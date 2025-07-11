#!/usr/bin/env bash
##############################################################################
#  run_auau_run3_qa_submit.sh         – submit (or locally test) the sPHENIX
#                                       EMCAL×sEPD×MBD QA correlation job
#
#  Data‑set selector  (first argument, default = run24auau)
#  ───────────────────────────────────────────────────────────────────────────
#   run24auau   → Au+Au Run‑24  CALO‑DSTs  (file prefix DST_CALO_run2auau_…)
#   run25auau   → Au+Au Run‑25  JET‑DSTs   (file prefix DST_JET‑…)
#
#  Operational modes (second argument)
#  ───────────────────────────────────────────────────────────────────────────
#   local                     : process first DST of first .list (quick test)
#   condor                    : submit every run
#   condor firstTen           : stop when grand‑total jobs > MAX_JOBS
#   condor round <N>          : submit only runs listed in run_segments/
#   condorTest                : like “condor” but after first run only
#   splitRunList <file.txt>   : create segment files (no submission)
#
#  Project tree
#  ───────────────────────────────────────────────────────────────────────────
#   PROJECT_BASE/
#       dst_list/             • run‑wise .list files (input for this script)
#       tmp_condor_lists/     • auto‑clean split chunks
#       run_segments/         • segment*.txt created by splitRunList
#       log/  stdout/  error/ • Condor logs
#       run_auau_run3_qa.sh   • worker executable
##############################################################################
set -euo pipefail
shopt -s extglob      # for pattern matching convenience
IFS=$'\n\t'

########################  COLOUR & LOG HELPERS  ###############################
ESC=$'\e['
CLR_R=${ESC}0\;31m ; CLR_G=${ESC}0\;32m ; CLR_Y=${ESC}1\;33m
CLR_B=${ESC}1\;34m ; CLR_BLD=${ESC}1m   ; CLR_RST=${ESC}0m
say()   { printf "${CLR_B}➜${CLR_RST} %s\n" "$*"; }
good()  { printf "${CLR_G}%s${CLR_RST}\n"   "$*"; }
warn()  { printf "${CLR_Y}⚠ %s${CLR_RST}\n" "$*" >&2; }
fatal() { printf "${CLR_R}✘ %s${CLR_RST}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO)"' ERR

##############################################################################
# 0. DATA‑SET SELECTION  (first CLI argument, default = run24auau)
##############################################################################
DATASET=${1:-run24auau}
shift || true                      # leave $@ holding the operational mode

case "$DATASET" in
  run24auau)
    FILE_PREFIX="DST_CALO_run2auau_new_2024p007"
    LIST_PATTERN="${FILE_PREFIX}-000*.list"
    LIST_FMT="${FILE_PREFIX}-000%05d.list"    # printf pattern
    ;;

  run25auau)
    FILE_PREFIX="DST_JET"
    LIST_PATTERN="${FILE_PREFIX}-000*.list"
    LIST_FMT="${FILE_PREFIX}-%08d.list"
    ;;

  *)
    fatal "Unknown data‑set selector '$DATASET' – use run24auau or run25auau"
    ;;
esac

##############################################################################
# 1. USER‑TUNABLE CONSTANTS
##############################################################################
CHUNK_SIZE=5      # files per Condor job
MAX_JOBS=10000    # global cap for “condor firstTen”
##############################################################################

##############################################################################
# 2. PROJECT CONSTANTS  (edit paths here only if the project moves)
##############################################################################
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"

DST_LIST_DIR="${PROJECT_BASE}/dst_list"
TMP_LIST_DIR="${PROJECT_BASE}/tmp_condor_lists"
EXEC="${PROJECT_BASE}/run_auau_run3_qa.sh"       # worker
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"

LOGDIR="${PROJECT_BASE}/log"
OUTDIR="${PROJECT_BASE}/stdout"
ERRDIR="${PROJECT_BASE}/error"
mkdir -p "$TMP_LIST_DIR" "$LOGDIR" "$OUTDIR" "$ERRDIR"

# —— optional splitting of run lists into “rounds” ————————————————
RUN_SPLIT_DIR="${PROJECT_BASE}/run_segments"
SEGMENT_PREFIX="${RUN_SPLIT_DIR}/runSegment_${DATASET}_"
mkdir -p "$RUN_SPLIT_DIR"
##############################################################################

##############################################################################
# 3. OPERATIONAL MODE PARSING
##############################################################################
mode="${1:-}"          # local | condor | condorTest | splitRunList
limitSwitch="${2:-}"   # optional second token

case "$mode" in
  local|condor|condorTest|splitRunList) ;;
  '')  fatal "Missing operational mode  (local | condor | condorTest | splitRunList)" ;;
  *)   fatal "Unknown operational mode '$mode'" ;;
esac
##############################################################################

##############################################################################
# 4. HELPER: split_run_list
##############################################################################
split_run_list() {
  local master="$1"
  [[ -f "$master" ]] || { warn "Run‑list not found → $master"; return 1; }

  say  "Splitting $(basename "$master") → $RUN_SPLIT_DIR"
  local seg=1 jobs=0 current="${SEGMENT_PREFIX}${seg}.txt"
  : > "$current"

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue

    local listFile
    if [[ "$DATASET" == run24auau ]]; then
      rn=$(printf "%05d" "$rn")
      listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$rn")"
    else
      rn=$(printf "%08d" "$rn")
      listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$rn")"
    fi

    [[ -f "$listFile" ]] || { warn "  – list for run $rn missing – skipped"; continue; }
    local nFiles; nFiles=$(wc -l < "$listFile")
    (( nFiles )) || { warn "  – list for run $rn empty   – skipped"; continue; }

    local nJobs=$(( (nFiles + CHUNK_SIZE - 1) / CHUNK_SIZE ))

    if (( jobs + nJobs > MAX_JOBS )); then
      good "  [SEGMENT $seg] closed with $jobs jobs"
      (( ++seg ))
      current="${SEGMENT_PREFIX}${seg}.txt"
      : > "$current"
      jobs=0
    fi

    echo "$rn" >> "$current"
    (( jobs += nJobs ))
  done < "$master"

  good "  [SEGMENT $seg] closed with $jobs jobs"
  good "[OK] splitting finished – files are in $RUN_SPLIT_DIR"
}
##############################################################################

##############################################################################
# 5. EARLY‑EXIT MODE : splitRunList
##############################################################################
if [[ "$mode" == "splitRunList" ]]; then
  split_run_list "$limitSwitch"   # 2nd arg is the run file
  exit 0
fi
##############################################################################

##############################################################################
# 6. VERBOSITY + GLOBAL CAP
##############################################################################
VERBOSE=0
[[ "$mode" == "condorTest" || ( "$mode" == "condor" && "$limitSwitch" == "firstTen" ) ]] && VERBOSE=1
vecho() { (( VERBOSE )) && echo -e "${CLR_B}•${CLR_RST} $*"; }

jobCap=0
[[ "$mode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobCap=$MAX_JOBS
submitted=0
##############################################################################

##############################################################################
# 7. ROUND‑N SUPPORT  (optional external run list)
##############################################################################
runListFile=""
if [[ "$mode" == "condor" && "$limitSwitch" == "round" && "${3:-}" =~ ^[0-9]+$ ]]; then
  runListFile="${SEGMENT_PREFIX}${3}.txt"
  [[ -f "$runListFile" ]] || fatal "Segment file $runListFile not found"
  say  "Round ${3} selected → using run list $(basename "$runListFile")"
fi
##############################################################################

##############################################################################
# 8. BUILD MAIN ARRAYS  (runs[]  +  listFiles[])
##############################################################################
runs=()
listFiles=()

if [[ -n "$runListFile" ]]; then
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    if [[ "$DATASET" == run24auau ]]; then
      rn=$(printf "%05d" "$rn")
    else
      rn=$(printf "%08d" "$rn")
    fi
    lf="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$rn")"
    [[ -f "$lf" ]] || { warn "List for run $rn missing – skipped"; continue; }
    runs+=( "$rn" ); listFiles+=( "$lf" )
  done < "$runListFile"
else
  mapfile -t listFiles < <(ls "${DST_LIST_DIR}"/${LIST_PATTERN} 2>/dev/null | sort)
  (( ${#listFiles[@]} )) || fatal "No .list files found in ${DST_LIST_DIR}"
  for f in "${listFiles[@]}"; do
    bn=${f##*-}; runs+=( "${bn%.list}" )
  done
fi
##############################################################################

##############################################################################
# 9. LOCAL MODE – single quick test
##############################################################################
if [[ "$mode" == "local" ]]; then
  [[ -n "${runs[0]:-}" ]] || fatal "No runs available for local mode"

  # Use run# and files from first .list unless a run number is forced
  if [[ -n "$limitSwitch" ]]; then
    if [[ "$DATASET" == run24auau ]]; then
      runNumber=$(printf "%05d" "$limitSwitch")
    else
      runNumber=$(printf "%08d" "$limitSwitch")
    fi
  else
    runNumber="${runs[0]}"
  fi

  maxEvt="${3:-0}"
  listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumber")"
  [[ -s "$listFile" ]] || fatal "List‑file $listFile not found or empty"

  firstDST=$(head -n1 "$listFile")
  say  "Local test  –  run $runNumber"
  say  "First DST   : $firstDST"
  say  "maxEvt      : $maxEvt (0 → all)"

  tmpList=$(mktemp "${TMP_LIST_DIR}/local_${runNumber}_XXXX.list")
  echo "$firstDST" > "$tmpList"

  "${EXEC}" "$runNumber" "$tmpList" 0 "$CONDOR_OUT_BASE" "$maxEvt"

  rm -f "$tmpList"
  exit 0
fi
##############################################################################

##############################################################################
# 10. CONDOR / CONDORTEST LOOP
##############################################################################
for idx in "${!runs[@]}"; do
  run=${runs[$idx]}
  masterList=${listFiles[$idx]}

  vecho "Considering run $run  (list: $(basename "$masterList"))"

  # ~~ split into CHUNK_SIZE‑line files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rm -f "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${run}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null)
  (( ${#chunks[@]} )) || { warn "Empty run $run – skipped."; continue; }

  if (( jobCap && submitted + ${#chunks[@]} > jobCap )); then
    good "Adding run $run would exceed cap ($jobCap) – stopping."
    break
  fi

  chunkNo=0
  for listFile in "${chunks[@]}"; do
    (( ++chunkNo ))
    [[ -s "$listFile" ]] || { warn "Empty chunk $listFile"; continue; }

    firstDST=$(head -n1 "$listFile")
    tag=$(basename "${firstDST%.root}")

    subFile="${TMP_LIST_DIR}/${tag}.sub"
    cat > "$subFile" <<EOS
universe      = vanilla
executable    = $EXEC
arguments     = $run $listFile \$(Cluster) $CONDOR_OUT_BASE
log           = ${LOGDIR}/${tag}.log
output        = ${OUTDIR}/${tag}.out
error         = ${ERRDIR}/${tag}.err
request_memory= 1500MB
+JobFlavour   = "tomorrow"
queue
EOS

    if condor_submit "$subFile" >/dev/null; then
      (( ++submitted ))
      vecho "  submitted chunk $chunkNo/${#chunks[@]}"
    else
      warn "condor_submit failed for $subFile"
    fi
  done

  vecho "Completed run $run – jobs now at $submitted"
  [[ "$mode" == "condorTest" ]] && { vecho "condorTest done."; break; }
done

good "Grand‑total Condor jobs submitted: $submitted"
[[ $jobCap -gt 0 ]] && say  "Launch cap in effect          : $jobCap"
##############################################################################
