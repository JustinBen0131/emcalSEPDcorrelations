#!/usr/bin/env bash
###############################################################################
#  run25qaInitial.sh – Run‑25 Au+Au JET‑DST quick QA
#                     ▸ per‑run runtime / GL1 counts
#                     ▸ cut‑flow (§ ≥ 5 min)
#                     ▸ trigger activity summary
###############################################################################
set -o errexit -o nounset -o pipefail
IFS=$'\n\t'; shopt -s nullglob

########################  USER SETTINGS  ######################################
DST_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"
RUN_LIST_FILE="run25auauCurrentDstRuns.txt"
MIN_RUNTIME=300            # seconds (5 min)
TOP_N_TRIG=20              # show N most‑frequent triggers
###############################################################################

########################  COLOUR HELPERS  #####################################
ESC=$'\e['
RED=${ESC}0\;31m  ; YEL=${ESC}1\;33m  ; GRN=${ESC}0\;32m
BLU=${ESC}1\;34m  ; RST=${ESC}0m
say()   { printf "${BLU}[STEP]${RST} %s\n"  "$*"; }
good()  { printf "${GRN}[ OK ]${RST} %s\n"  "$*"; }
warn()  { printf "${YEL}[WARN]${RST} %s\n" "$*"; }
fatal() { printf "${RED}✘ %s${RST}\n"     "$*" >&2; exit 2; }
trap 'fatal "Script aborted (line $LINENO) – $BASH_COMMAND"' ERR

########################  BINARIES & HELPERS  #################################
for b in psql column bc; do command -v "$b" >/dev/null || fatal "$b not found"; done
readonly PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)

sql()          { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
num_or_zero()  { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }
safe_mapfile() { mapfile -t "$@" || true; }              # ignore empty result sets

sql 'SELECT 1;' >/dev/null || fatal "Cannot reach DAQ DB"
good "DB connectivity OK."

######################## 1. BUILD RUN LIST ####################################
say "Scanning $DST_LIST_DIR for .list files"
runs=()
for f in "$DST_LIST_DIR"/{DST_JET,dst_jet}-*.list; do
  rn=${f##*-}; rn=${rn%.list}
  [[ $rn =~ ^[0-9]{8}$ ]] && runs+=("$rn") || warn "Ignore: ${f##*/}"
done
((${#runs[@]})) || fatal "No JET‑DST lists found"
printf '%s\n' "${runs[@]}" | sort -u >"$RUN_LIST_FILE"
good "Run list ➔ $RUN_LIST_FILE  (${#runs[@]} runs)"

######################## 2. MAIN LOOP #########################################
say "Querying runs …"
header=$'Run\trunTime[s]\tGL1_evt\tON\tOFF\t<5m\tLegend'

rows=()
totalTime=0 totalEvt=0 shortCnt=0
dropped=()

for run in "${runs[@]}"; do
  say "• $run"

  # -- runtime ---------------------------------------------------------------
  rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT
            FROM run WHERE runnumber=$run;" | tr -d '[:space:]')
  rt=$(num_or_zero "$rt")
  flag=''; (( rt < MIN_RUNTIME )) && { flag='YES'; dropped+=("$run"); ((++shortCnt)); }

  # -- GL1 events ------------------------------------------------------------
  ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
            FROM filelist
            WHERE runnumber=$run
              AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d '[:space:]')
  ev=$(num_or_zero "$ev")

  # -- triggers --------------------------------------------------------------
  safe_mapfile trgRows < <(sql "SELECT scaled FROM gl1_scalers WHERE runnumber=$run;")
  on=0 off=0 legend='–'
  for s in "${trgRows[@]}"; do
    [[ $s =~ ^-?[0-9]+$ && $s -gt 0 ]] && (( ++on )) || (( ++off ))
  done
  (( on )) && legend='✔' || legend='✖'

  rows+=("$run"$'\t'"$rt"$'\t'"$ev"$'\t'"$on"$'\t'"$off"$'\t'"$flag"$'\t'"$legend")
  (( totalTime += rt, totalEvt += ev ))
done

######################## 3. PER‑RUN TABLE #####################################
echo -e "\n${GRN}================ Per‑run summary ================${RST}"
printf '%s\n' "$header" "${rows[@]}" | column -t -s $'\t'

######################## 4. CUT‑FLOW ##########################################
passRuns=$(( ${#runs[@]} - shortCnt ))
passTime=0 passEvt=0
for row in "${rows[@]}"; do
  IFS=$'\t' read -r _ rt ev _ _ flag _ <<<"$row"
  [[ $flag == YES ]] && continue
  (( passTime += rt, passEvt += ev ))
done

echo -e "\n${GRN}================ Cut‑flow (≥5 min) ===============${RST}"
printf "%-6s %5s %13s %12s\n" "Stage" "Runs" "Events" "Runtime[h]"
printf "%-6s %5d %13d %12.2f\n" "Raw" "${#runs[@]}" "$totalEvt" \
       "$(bc -l <<<"$totalTime/3600")"
printf "%-6s %5d %13d %12.2f\n" "Cut" "$passRuns"   "$passEvt"  \
       "$(bc -l <<<"$passTime/3600")"

######################## 5. DROPPED RUNS ######################################
echo -e "\n${YEL}Runs <5 min:${RST}"
((${#dropped[@]})) && printf '%s\n' "${dropped[@]}" | pr -3 -t || echo "(none)"

######################## 6. TRIGGER SUMMARY ###################################
echo -e "\n${GRN}================ Trigger totals =================${RST}"
printf '%-35s %8s %15s\n' "Trigger" "ONruns" "ΣScaled"
declare -A trigRun trigSum           # ← nounset‑safe ++/+= below

for run in "${runs[@]}"; do
  while IFS=$'\t' read -r name scaled; do
    [[ -z $name ]] && continue
    [[ $scaled =~ ^-?[0-9]+$ && $scaled -gt 0 ]] || continue
    trigRun["$name"]=$(( ${trigRun["$name"]:-0} + 1 ))
    trigSum["$name"]=$(( ${trigSum["$name"]:-0} + scaled ))
  done < <(sql "SELECT t.triggername, s.scaled
                FROM gl1_scalers s
                JOIN gl1_triggernames t
                  ON s.index = t.index
                 AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                WHERE s.runnumber = $run;")
done

for trg in "${!trigRun[@]}"; do
  printf '%-35s %8d %15d\n' "$trg" "${trigRun[$trg]}" "${trigSum[$trg]}"
done | sort -k2 -nr | head -"$TOP_N_TRIG"

######################## 7. OVERALL TOTALS ####################################
echo -e "${GRN}----------------------------------------------------${RST}"
printf 'Runs analysed : %d\n'   "${#runs[@]}"
printf 'Short (<5 m)  : %d\n'   "$shortCnt"
printf 'Total runtime : %d s (≈ %.2f h)\n' "$totalTime" "$(bc -l <<<"$totalTime/3600")"
printf 'Total GL1 evt : %d\n'   "$totalEvt"
echo -e "${GRN}====================================================${RST}"
good "Done."
