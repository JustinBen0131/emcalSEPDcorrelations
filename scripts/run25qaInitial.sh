#!/usr/bin/env bash
###############################################################################
#  run25qaInitial.sh – Run‑25 Au+Au JET‑DST quick QA
#                     (runtime cut + trigger census & matrix + cut‑flow)
#
#  ① Builds run25auauCurrentDstRuns.txt from *.list in $DST_LIST_DIR
#  ② Per run queries DAQ replica for
#       • runtime      := ertimestamp – brtimestamp   [ s ]
#       • GL1 evt      := Σ( lastevent–firstevent+1 ) in *.evt
#       • every GL1 trigger’s “scaled” value
#  ③ Prints
#       • per‑run summary
#       • cut‑flow table  ( ≥5 min vs <5 min )
#       • list of <5 min runs
#       • trigger census  (ON|OFF|ABSENT)
#       • run × trigger matrix (✔/✖/–)
#
#  2025‑07‑11   –   sPHENIX offline QA utilities
###############################################################################
set -o errexit -o nounset -o pipefail
IFS=$'\n\t' ; shopt -s nullglob

########################  USER SETTINGS  ######################################
DST_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"
RUN_LIST_FILE="run25auauCurrentDstRuns.txt"
MIN_RUNTIME=300            # 5 min
###############################################################################

########################  COLOUR HELPERS  #####################################
ESC=$'\e['
RED=${ESC}0\;31m ; YEL=${ESC}1\;33m ; GRN=${ESC}0\;32m ; BLU=${ESC}1\;34m ; RST=${ESC}0m
say()   { printf "${BLU}[STEP]${RST} %s\n" "$*"; }
good()  { printf "${GRN}[ OK ]${RST} %s\n"  "$*"; }
warn()  { printf "${YEL}[WARN]${RST} %s\n" "$*"; }
fatal() { printf "${RED}✘ %s${RST}\n"     "$*" >&2; exit 2; }
trap 'fatal "Script aborted (line $LINENO) – check previous messages."' ERR

########################  BINARIES & HELPERS  #################################
command -v psql   >/dev/null || fatal "psql not found in \$PATH"
command -v column >/dev/null || fatal "`column` (util‑linux) missing"
command -v bc     >/dev/null || fatal "`bc` not found – install bc"

readonly PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)

# --- run SQL, always exit 0, return rows ------------------------------------
db() { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }

# --- numeric or zero (prevents (( … )) crashes) -----------------------------
num_or_zero() {
  [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf '0'
}

"${PSQL[@]}" -c 'SELECT 1;' &>/dev/null || \
  fatal "Cannot reach DAQ DB – check network / Kerberos / .pgpass"
good "DB connectivity OK."

########################  1. BUILD RUN LIST  ##################################
say "Scanning $DST_LIST_DIR for JET‑DST .list files"
runs=()
for f in "$DST_LIST_DIR"/{DST_JET,dst_jet}-*.list; do
  bn=${f##*/}; run=${bn#*-}; run=${run%.list}
  [[ $run =~ ^[0-9]{8}$ ]] && runs+=("$run") || warn "Ignored odd file: $bn"
done
((${#runs[@]})) || fatal "No valid JET‑DST .list files found."
printf '%s\n' "${runs[@]}" | sort -u >"$RUN_LIST_FILE"
good "Run list saved to $RUN_LIST_FILE  ( $(wc -l <"$RUN_LIST_FILE") runs )"

########################  2. MAIN LOOP  #######################################
say "Starting per‑run queries …"

header=$'runNumber\trunTime[s]\tGL1_evt\tON\tOFF\t<5min\tON_triggers\tOFF_triggers\n'
rows=()

declare -A onCnt offCnt absCnt status              # trigger meta
declare -a droppedRuns
totalTime=0 totalEvt=0 shortCnt=0 idx=0

while read -r run; do
  ((++idx))
  say "[RUN $idx] Processing run $run"

  # ---------- runtime (s) ---------------------------------------------------
  runtime_raw=$(db "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT
                    FROM run WHERE runnumber=$run;")
  runtime=$(num_or_zero "${runtime_raw//[[:space:]]/}")
  if (( runtime < MIN_RUNTIME )); then
    short='YES'; droppedRuns+=("$run"); ((shortCnt++))
  else
    short=''
  fi

  # ---------- GL1 physics events -------------------------------------------
  evt_raw=$(db "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
                FROM filelist
                WHERE runnumber=$run
                  AND filename LIKE '%GL1_physics_gl1daq%.evt';")
  gl1evt=$(num_or_zero "${evt_raw//[[:space:]]/}")

  # ---------- trigger scalers ----------------------------------------------
  mapfile -t trgRows < <(db "SELECT t.triggername, s.scaled
                             FROM gl1_scalers s
                             JOIN gl1_triggernames t
                               ON s.index=t.index
                              AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                             WHERE s.runnumber=$run;")

  declare -A seen=(); onList=(); offList=()
  for row in "${trgRows[@]}"; do
    IFS=$'\t' read -r trig scaled <<<"$row"
    seen["$trig"]=1
    if [[ $scaled =~ ^-?[0-9]+$ ]] && (( scaled > 0 )); then
      status["$run|$trig"]="✔"; (( onCnt["$trig"]++ ))
      onList+=("${trig}(${scaled})")
    else
      status["$run|$trig"]="✖"; (( offCnt["$trig"]++ ))
      offList+=("$trig")
    fi
  done
  # any trigger not seen this run but observed elsewhere → ABSENT
  for t in "${!onCnt[@]}" "${!offCnt[@]}"; do
    [[ ${seen[$t]+x} ]] || { status["$run|$t"]="–"; (( absCnt["$t"]++ )); }
  done

  totalTime=$(( totalTime + runtime ))
  totalEvt=$(( totalEvt  + gl1evt ))

  rows+=("$run\t$runtime\t$gl1evt\t${#onList[@]}\t${#offList[@]}\t$short\t\
$(IFS=,;echo "${onList[*]}")\t$(IFS=,;echo "${offList[*]}")")
done <"$RUN_LIST_FILE"

########################  3. PER‑RUN SUMMARY  #################################
echo -e "\n${GRN}==================== Per‑Run Summary ====================${RST}"
printf '%s' "$header"
printf '%s\n' "${rows[@]}" | column -t -s $'\t'

########################  4. CUT‑FLOW (runtime)  ##############################
passRuns=$(( ${#runs[@]} - shortCnt ))
passTime=0 passEvt=0
for row in "${rows[@]}"; do
  IFS=$'\t' read -r _run rt ev _ _ shortFlag _ _ <<<"$row"
  [[ $shortFlag == YES ]] && continue
  (( passTime += rt, passEvt += ev ))
done

echo -e "\n${GRN}==================== Cut‑flow (runtime) ===================${RST}"
printf 'Stage\tRuns\tEvents\tRuntime[h]\n'
printf 'Raw  \t%d\t%d\t%.2f\n'  "${#runs[@]}" "$totalEvt" "$(bc -l <<<"$totalTime/3600")"
printf '≥5 m \t%d\t%d\t%.2f\n'  "$passRuns"   "$passEvt"  "$(bc -l <<<"$passTime/3600")"

########################  5. DROPPED RUNS (<5 min) ###########################
echo -e "\n${YEL}Runs dropped for runtime < ${MIN_RUNTIME}s:${RST}"
((${#droppedRuns[@]})) && printf '%s\n' "${droppedRuns[@]}" | column || echo "(none)"

########################  6. TRIGGER CENSUS  ##################################
echo -e "\n${GRN}==================== Trigger Census ====================${RST}"
printf 'TriggerName\tON\tOFF\tABSENT\n'
for t in $(printf '%s\n' "${!onCnt[@]}" "${!offCnt[@]}" | sort -u); do
  printf '%s\t%d\t%d\t%d\n' "$t" "${onCnt[$t]:-0}" "${offCnt[$t]:-0}" "${absCnt[$t]:-0}"
done | column -t -s $'\t'

########################  7. RUN × TRIGGER MATRIX #############################
echo -e "\n${GRN}============== Run × Trigger Status Matrix ==============${RST}"
triggers=( $(printf '%s\n' "${!onCnt[@]}" "${!offCnt[@]}" | sort -u) )
printf 'runNumber'; for t in "${triggers[@]}"; do printf '\t%s' "$t"; done; echo
for run in "${runs[@]}"; do
  printf '%s' "$run"
  for t in "${triggers[@]}"; do
    printf '\t%s' "${status["$run|$t"]:-–}"
  done
  echo
done | column -t -s $'\t'

########################  8. OVERALL TOTALS  ##################################
echo -e "${GRN}------------------------------------------------------------${RST}"
printf 'Runs analysed : %d\n'   "${#runs[@]}"
printf 'Short (<5 m)  : %d\n'   "$shortCnt"
printf 'Total runtime : %d s  (%.2f h)\n' "$totalTime" "$(bc -l <<<"$totalTime/3600")"
printf 'Total GL1 evt : %d\n'   "$totalEvt"
echo -e "${GRN}============================================================${RST}"
good "Done."
