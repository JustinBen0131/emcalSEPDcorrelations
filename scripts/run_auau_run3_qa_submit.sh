#!/usr/bin/env bash
##############################################################################
#  run_auau_run3_qa_submit.sh         – submit (or locally test) the sPHENIX
#                                       EMCAL×sEPD×MBD QA correlation job
##############################################################################
set -euo pipefail
shopt -s extglob
IFS=$'\n\t'

########################  COLOUR & LOG HELPERS  ###############################
ESC=$'\e['
CLR_R=${ESC}0\;31m ; CLR_G=${ESC}0\;32m ; CLR_Y=${ESC}1\;33m
CLR_B=${ESC}1\;34m ; CLR_RST=${ESC}0m
say()   { printf "${CLR_B}➜${CLR_RST} %s\n" "$*"; }
good()  { printf "${CLR_G}%s${CLR_RST}\n"   "$*"; }
warn()  { printf "${CLR_Y}⚠ %s${CLR_RST}\n" "$*" >&2; }
fatal() { printf "${CLR_R}✘ %s${CLR_RST}\n" "$*" >&2; exit 1; }
trap 'fatal "Script aborted (line $LINENO)"' ERR

##############################################################################
# 0. DATA‑SET SELECTION
##############################################################################
DATASET=${1:-run24auau}        # run24auau | run25auau
shift || true                  # consume it

# ---------- NEW: optional DST‑type argument for run25auau -------------------
#   • dstjet      → DST_JET      (default if omitted)
#   • dstjetcalo  → DST_JETCALO
DSTTYPE=dstjet
if [[ "$DATASET" == run25auau && "${1:-}" =~ ^(dstjet|dstjetcalo)$ ]]; then
  DSTTYPE=${1,,}   # lower‑case for robustness
  shift            # consume the dst‑type token
fi
# ---------------------------------------------------------------------------

case "$DATASET" in
  run24auau)
    FILE_PREFIX="DST_CALO_run2auau_new_2024p007"
    LIST_PATTERN="${FILE_PREFIX}-000*.list"
    LIST_FMT="${FILE_PREFIX}-000%05d.list"
    PAD_FMT="%05d"
    ;;
  run25auau)
    # --------------------- NEW: derive prefix from $DSTTYPE -----------------
    case "$DSTTYPE" in
      dstjet)      FILE_PREFIX="DST_JET" ;;
      dstjetcalo)  FILE_PREFIX="DST_JETCALO" ;;
      *)           fatal "BUG: unhandled DSTTYPE ‘$DSTTYPE’" ;;
    esac
    # -----------------------------------------------------------------------
    LIST_PATTERN="${FILE_PREFIX}-000*.list"
    LIST_FMT="${FILE_PREFIX}-%08d.list"
    PAD_FMT="%08d"
    ;;
  *)
    fatal "Unknown data‑set selector '$DATASET' – use run24auau or run25auau"
    ;;
esac

##############################################################################
# 1. USER CONSTANTS
##############################################################################
CHUNK_SIZE=2
MAX_JOBS=10000
##############################################################################

##############################################################################
# 2. PATH CONSTANTS
##############################################################################
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
DST_LIST_DIR="${PROJECT_BASE}/dst_list"
TMP_LIST_DIR="${PROJECT_BASE}/tmp_condor_lists"
EXEC="${PROJECT_BASE}/run_auau_run3_qa.sh"
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
LOGDIR="${PROJECT_BASE}/log"
OUTDIR="${PROJECT_BASE}/stdout"
ERRDIR="${PROJECT_BASE}/error"
mkdir -p "$TMP_LIST_DIR" "$LOGDIR" "$OUTDIR" "$ERRDIR"

RUN_SPLIT_DIR="${PROJECT_BASE}/run_segments"
SEGMENT_PREFIX="${RUN_SPLIT_DIR}/runSegment_${DATASET}_"
mkdir -p "$RUN_SPLIT_DIR"
##############################################################################

##############################################################################
# 3. OPERATIONAL MODE PARSING
##############################################################################
mode="${1:-}"          # local | condor | condorTest | splitRunList
limitSwitch="${2:-}"   # optional 2nd token

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

  while IFS= read -r raw; do
    [[ -z "$raw" || "$raw" =~ ^# ]] && continue
    local runNumDec=$((10#$raw))
    local listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
    [[ -f "$listFile" ]] || { warn "  – list for run $raw missing – skipped"; continue; }

    local nFiles; nFiles=$(wc -l < "$listFile")
    (( nFiles )) || { warn "  – list for run $raw empty   – skipped"; continue; }

    local nJobs=$(( (nFiles + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    if (( jobs + nJobs > MAX_JOBS )); then
      good "  [SEGMENT $seg] closed with $jobs jobs"
      (( ++seg ))
      current="${SEGMENT_PREFIX}${seg}.txt"
      : > "$current"
      jobs=0
    fi

    printf "$PAD_FMT\n" "$runNumDec" >>"$current"
    (( jobs += nJobs ))
  done < "$master"

  good "  [SEGMENT $seg] closed with $jobs jobs"
  good "[OK] splitting finished – files are in $RUN_SPLIT_DIR"
}
##############################################################################

##############################################################################
# 5. EARLY‑EXIT: splitRunList
##############################################################################
if [[ "$mode" == "splitRunList" ]]; then
  split_run_list "$limitSwitch"
  exit 0
fi
##############################################################################

##############################################################################
# 6. VERBOSITY / CAP
##############################################################################
VERBOSE=0
[[ "$mode" == "condorTest" || ( "$mode" == "condor" && "$limitSwitch" == "firstTen" ) ]] && VERBOSE=1
vecho() { (( VERBOSE )) && echo -e "${CLR_B}•${CLR_RST} $*"; }

jobCap=0
[[ "$mode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobCap=$MAX_JOBS
submitted=0
##############################################################################

##############################################################################
# 7. ROUND‑N OR GOLDEN‑LIST SELECTION
##############################################################################
runListFile=""

if [[ "$mode" == "condor" && "$limitSwitch" == "round" && "${3:-}" =~ ^[0-9]+$ ]]; then
  runListFile="${SEGMENT_PREFIX}${3}.txt"
  [[ -f "$runListFile" ]] || fatal "Segment file $runListFile not found"
  say  "Round ${3} selected → using run list $(basename "$runListFile")"
fi

if [[ "$DATASET" == run25auau && -z "$runListFile" ]]; then
  for p in "${PROJECT_BASE}" .; do
    [[ -f "$p/run25GoldenRuns.txt" ]] && runListFile="$p/run25GoldenRuns.txt" && break
  done
  [[ -n "$runListFile" ]] && say  "run25auau selected – using golden run list $(basename "$runListFile")"
fi
##############################################################################

##############################################################################
# 8. BUILD MAIN ARRAYS
##############################################################################
runs=()
listFiles=()

if [[ -n "$runListFile" ]]; then
  while IFS= read -r raw; do
    [[ -z "$raw" || "$raw" =~ ^# ]] && continue
    runNumDec=$((10#$raw))
    runNumPad=$(printf "$PAD_FMT" "$runNumDec")
    listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
    [[ -f "$listFile" ]] || { warn "List for run $runNumPad missing – skipped"; continue; }
    runs+=( "$runNumPad" ); listFiles+=( "$listFile" )
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
# 9. LOCAL MODE
##############################################################################
if [[ "$mode" == "local" ]]; then
  [[ -n "${runs[0]:-}" ]] || fatal "No runs available for local mode"

  if [[ -n "$limitSwitch" ]]; then
    runNumDec=$((10#$limitSwitch))
    runNumber=$(printf "$PAD_FMT" "$runNumDec")
  else
    runNumber="${runs[0]}"
    runNumDec=$((10#$runNumber))
  fi

  maxEvt="${3:-0}"
  listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
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
  runPad=${runs[$idx]}
  runDec=$((10#$runPad))
  masterList=${listFiles[$idx]}

  vecho "Considering run $runPad  (list: $(basename "$masterList"))"

  rm -f "${TMP_LIST_DIR}/run${runPad}_chunk_"* 2>/dev/null || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${runPad}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${runPad}_chunk_"* 2>/dev/null)
  (( ${#chunks[@]} )) || { warn "Empty run $runPad – skipped."; continue; }

  if (( jobCap && submitted + ${#chunks[@]} > jobCap )); then
    good "Adding run $runPad would exceed cap ($jobCap) – stopping."
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
arguments     = $runPad $listFile \$(Cluster) $CONDOR_OUT_BASE
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

  vecho "Completed run $runPad – jobs now at $submitted"
  [[ "$mode" == "condorTest" ]] && { vecho "condorTest done."; break; }
done

good "Grand‑total Condor jobs submitted: $submitted"
[[ $jobCap -gt 0 ]] && say  "Launch cap in effect          : $jobCap"
##############################################################################
