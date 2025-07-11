#!/usr/bin/env bash
###############################################################################
#  run25qaInitial.sh  –  Run-25 Au+Au JET-DST quick QA  (trigger census edition)
#
#  • Builds run25auauCurrentDstRuns.txt from *.list inside dst_list/
#  • For each run queries DAQ replica for runtime, GL1-event count and all
#    GL1 triggers (with their ‘scaled’ value).
#  • Produces two tables:
#       1) per-run summary (identical to previous behaviour)
#       2) trigger census:
#            TriggerName | ON  | OFF | ABSENT
#         where  ON   = scaled > 0
#                OFF  = scaled = -1   (trigger deliberately off)
#                ABSENT = no row for that trigger in given run
#
#  2025-07-11
###############################################################################
set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob

########################  USER SETTINGS  ######################################
DST_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"
RUN_LIST_FILE="run25auauCurrentDstRuns.txt"
MIN_RUNTIME=300            # 5 min  (seconds)
###############################################################################

########################  COLOUR HELPERS  #####################################
ESC=$'\e['
RED=${ESC}0\;31m; YEL=${ESC}1\;33m; GRN=${ESC}0\;32m; BLU=${ESC}1\;34m; RST=${ESC}0m
say()   { printf "${BLU}[STEP]${RST} %s\n" "$*"; }
good()  { printf "${GRN}[ OK ]${RST} %s\n"  "$*"; }
warn()  { printf "${YEL}[WARN]${RST} %s\n" "$*"; }
fatal() { printf "${RED}✘ %s${RST}\n"     "$*" >&2; exit 2; }

trap 'fatal "Script aborted (line $LINENO) – check previous messages."' ERR

########################  PREREQUISITES  ######################################
command -v psql   >/dev/null || fatal "psql not found in \$PATH"
command -v column >/dev/null || fatal "column (util-linux) missing"
PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)

"${PSQL[@]}" -c "SELECT 1;" &>/dev/null \
  || fatal "Cannot reach DB – check network / Kerberos / .pgpass"
good "DB connectivity OK."

########################  1. BUILD RUN LIST  ##################################
say "Scanning $DST_LIST_DIR for JET-DST .list files"

runs=()
for f in "$DST_LIST_DIR"/{DST_JET,dst_jet}-*.list; do
  bn=${f##*/}; run=${bn#*-}; run=${run%.list}
  [[ $run =~ ^[0-9]{8}$ ]] && runs+=("$run") || warn "Ignored odd file: $bn"
done
((${#runs[@]})) || fatal "No valid JET-DST .list files found."

printf '%s\n' "${runs[@]}" | sort -u >"$RUN_LIST_FILE"
good "Run list saved to $RUN_LIST_FILE  ( $(wc -l <"$RUN_LIST_FILE") runs )"

########################  2. MAIN LOOP  #######################################
say "Starting per-run queries …"

header=$'runNumber\trunTime[s]\tGL1_evt\t<5min\tTriggers_ON(scaled)\n'
rows=()
shortCnt=0; totalTime=0; totalEvt=0

# --- trigger census structures ---------------------------------------------
declare -A trigON   # how many runs have trigger scaled>0
declare -A trigOFF  # how many runs have trigger scaled=-1
declare -A trigABS  # how many runs *lack* that trigger row

# helper: safe increment for associative arrays with “set -u”
inc() { local -n A=$1; local k=$2; A["$k"]=$(( ${A["$k"]:-0} + 1 )); }

idx=0
while read -r run; do
  ((++idx))
  say "[RUN $idx] Processing run $run"

  # ---------- runtime -------------------------------------------------------
  sql_rt="SELECT COALESCE(FLOOR(EXTRACT(EPOCH FROM (ertimestamp - brtimestamp)))::INT,0)
          FROM run WHERE runnumber=$run;"
  echo "[SQL] ${sql_rt//[$'\n\t']/ }"
  runtime=$("${PSQL[@]}" -c "$sql_rt" | tr -d '[:space:]')
  [[ -z $runtime ]] && runtime=0
  (( runtime < MIN_RUNTIME )) && { short="YES"; ((shortCnt++)); } || short=""

  # ---------- GL1 physics events -------------------------------------------
  sql_evt="SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
           FROM filelist
           WHERE runnumber=$run
             AND filename LIKE '%GL1_physics_gl1daq%.evt';"
  echo "[SQL] ${sql_evt//[$'\n\t']/ }"
  gl1evt=$("${PSQL[@]}" -c "$sql_evt" | tr -d '[:space:]')
  [[ -z $gl1evt ]] && gl1evt=0

  # ---------- triggers ------------------------------------------------------
  sql_trig="SELECT t.triggername, s.scaled
            FROM gl1_scalers s
            JOIN gl1_triggernames t
              ON s.index=t.index
             AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
            WHERE s.runnumber=$run;"
  echo "[SQL] ${sql_trig//[$'\n\t']/ }"
  mapfile -t trgRows < <("${PSQL[@]}" -c "$sql_trig")

  # track which triggers appeared at all for ABSENT counting
  declare -A seenInRun=()

  trigStr=""
  for row in "${trgRows[@]}"; do
    trig=${row%%$'\t'*}; scaled=${row##*$'\t'}
    seenInRun["$trig"]=1
    if [[ $scaled =~ ^-?[0-9]+$ ]]; then
      case $scaled in
        ''|*[!0-9-]*) ;;  # non-numeric, ignore
        *[!0-9]* )      ;; # safety
        *)
          if (( scaled > 0 )); then
            trigStr+="${trig}(${scaled}),"
            inc trigON  "$trig"
          elif (( scaled == -1 )); then
            inc trigOFF "$trig"
          fi
          ;;
      esac
    fi
  done
  # mark absences
  for trig in "${!trigON[@]}" "${!trigOFF[@]}"; do
    [[ -n ${seenInRun["$trig"]+set} ]] || inc trigABS "$trig"
  done
  trigStr=${trigStr%,}

  # ---------- accumulate totals --------------------------------------------
  totalTime=$(( totalTime + runtime ))
  totalEvt=$(( totalEvt  + gl1evt ))
  rows+=( "$run\t$runtime\t$gl1evt\t$short\t$trigStr" )
done <"$RUN_LIST_FILE"

########################  3. PER-RUN TABLE ####################################
echo -e "\n${GRN}========== Run-25 Jet-DST Per-Run Summary ==========${RST}"
printf '%s' "$header"
printf '%s\n' "${rows[@]}" | column -t -s $'\t'

########################  4. TRIGGER CENSUS ###################################
echo -e "\n${GRN}========== Trigger Census (per ${#rows[@]} runs) ==========${RST}"
printf 'TriggerName\tON\tOFF\tABSENT\n'
allTrig=($(printf '%s\n' "${!trigON[@]}" "${!trigOFF[@]}" "${!trigABS[@]}" | sort -u))
for t in "${allTrig[@]}"; do
  printf '%s\t%d\t%d\t%d\n' \
    "$t" \
    "${trigON[$t]:-0}" \
    "${trigOFF[$t]:-0}" \
    "${trigABS[$t]:-0}"
done | column -t -s $'\t'

########################  5. OVERALL STATS ####################################
echo "------------------------------------------------------------"
printf 'Runs analysed : %d\n'    "${#rows[@]}"
printf 'Short (<5 m)  : %d\n'    "$shortCnt"
printf 'Total runtime : %d s  (%.2f h)\n' "$totalTime" "$(bc -l <<<"$totalTime/3600")"
printf 'Total GL1 evt : %d\n'   "$totalEvt"
echo "============================================================"
good "Done."
###############################################################################
