#!/usr/bin/env bash
##############################################################################
#  run_auau_run3_qa_submit.sh
#
#  Modes
#  ───────────────────────────────────────────────────────────────────────────
#   local                     : process first DST of first .list (quick test)
#   condor                    : submit every run
#   condor firstTen           : stop when grand‑total jobs > MAX_JOBS
#   condor round  <N>         : submit only runs listed in run_segments/
#   condorTest                : like “condor” but after first run only
#   splitRunList <file.txt>   : create segment files (no submission)
#
#  ───────────────────────────────────────────────────────────────────────────
#   PROJECT_BASE
#     ├── dst_list/               (input *.list files – one per run)
#     ├── tmp_condor_lists/       (auto‑clean split chunks)
#     ├── run_segments/           (segment*.txt created by splitRunList)
#     ├── log/  stdout/  error/   (Condor logs)
#     └── run_auau_run3_qa.sh     (worker executable)
##############################################################################
set -euo pipefail

#################################### d##########################################
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
SEGMENT_PREFIX="${RUN_SPLIT_DIR}/runSegment_"
mkdir -p "$RUN_SPLIT_DIR"
##############################################################################

##############################################################################
# 3. ARGUMENT PARSING
##############################################################################
mode="${1:-}"
limitSwitch="${2:-}"

if [[ "$mode" != "local"        && \
      "$mode" != "condor"       && \
      "$mode" != "condorTest"   && \
      "$mode" != "splitRunList" ]]; then
  cat <<EOF
Usage:
  $0 local
  $0 condor [firstTen]
  $0 condor round <N>           (submit only segment N)
  $0 condorTest
  $0 splitRunList <runFile.txt> (create segment files)
EOF
  exit 1
fi
##############################################################################

##############################################################################
# 4. HELPER: split_run_list
##############################################################################
split_run_list() {
  local master="$1"
  [[ -f "$master" ]] || { echo "[ERROR] run‑list not found → $master"; return 1; }

  echo "[INFO] Splitting $(basename "$master") → $RUN_SPLIT_DIR"
  local seg=1 jobs=0 current="${SEGMENT_PREFIX}${seg}.txt"
  : > "$current"

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    rn=$(printf "%05d" "$rn")   # normalise to 5‑digits (e.g. 54128)

    local lf="${DST_LIST_DIR}/DST_CALO_run2auau_new_2024p007-000${rn}.list"
    [[ -f "$lf" ]] || { echo "[WARN] No list for run $rn – skipped"; continue; }

    local nFiles; nFiles=$(wc -l < "$lf")
    (( nFiles )) || { echo "[WARN] Empty list for run $rn – skipped"; continue; }

    local nJobs=$(( (nFiles + CHUNK_SIZE - 1) / CHUNK_SIZE ))

    # overflow guard – open new segment
    if (( jobs + nJobs > MAX_JOBS )); then
      echo "  [SEGMENT $seg] closed with $jobs jobs"
      (( ++seg ))
      current="${SEGMENT_PREFIX}${seg}.txt"
      : > "$current"
      jobs=0
    fi

    echo "$rn" >> "$current"
    (( jobs += nJobs ))
  done < "$master"

  echo "  [SEGMENT $seg] closed with $jobs jobs"
  echo "[OK] splitting finished – files are in $RUN_SPLIT_DIR"
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
vecho() { (( VERBOSE )) && echo "$*"; }

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
  [[ -f "$runListFile" ]] || { echo "[ERROR] Segment file $runListFile not found"; exit 2; }
  echo "[INFO] Round $3 selected → using run list $runListFile"
fi
##############################################################################

##############################################################################
# 8. BUILD MAIN ARRAYS  (runs[]  +  listFiles[])
##############################################################################
if [[ -n "$runListFile" ]]; then
  # — use an external list (segment file) —
  runs=()
  listFiles=()
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    rn=$(printf "%05d" "$rn")
    lf="${DST_LIST_DIR}/DST_CALO_run2auau_new_2024p007-000${rn}.list"
    [[ -f "$lf" ]] || { echo "[WARN] No list for run $rn – skipped"; continue; }
    runs+=( "$rn" )
    listFiles+=( "$lf" )
  done < "$runListFile"
else
  # — original directory scan —
  mapfile -t listFiles < <(ls "${DST_LIST_DIR}"/DST_CALO_run2auau_new_2024p007-000*.list 2>/dev/null | sort)
  if (( ${#listFiles[@]} == 0 )); then
    echo "[FATAL] No .list files found in ${DST_LIST_DIR}" >&2
    exit 2
  fi
  runs=()
  for f in "${listFiles[@]}"; do
    bn=${f##*-000}; runs+=( "${bn%.list}" )
  done
fi
##############################################################################

# ─────────────────────────────────────────────────────────────
# 9. LOCAL MODE – single quick test               (✱ modified)
# ─────────────────────────────────────────────────────────────
if [[ "$mode" == "local" ]]; then
  maxEvt="${2:-0}"                           # new (0 ⇒ all events)

  firstList="${listFiles[0]}"
  firstRun="${runs[0]}"
  firstDST="$(head -n1 "$firstList")" || { echo "[ERROR] Empty $firstList"; exit 3; }

  echo "[INFO] Local mode – Run=$firstRun"
  echo "[INFO] First DST  : $firstDST"

  tmpList=$(mktemp "${TMP_LIST_DIR}/local_${firstRun}_XXXX.list")
  echo "$firstDST" > "$tmpList"

  #               ↓   ↓           ↓            ↓                ↓ NEW
  "${EXEC}" "$firstRun" "$tmpList" 0 "$CONDOR_OUT_BASE" "$maxEvt"
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

  vecho "[VERBOSE] Considering run $run  (file: $(basename "$masterList"))"

  # ~~ split into CHUNK_SIZE‑line files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rm -f "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${run}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null)
  nChunks=${#chunks[@]}
  (( nChunks )) || { echo "[WARN] Empty run $run – skipping."; continue; }

  # —— global cap guard ————————————————————————————————————————————————
  if (( jobCap && submitted + nChunks > jobCap )); then
    echo "[INFO] Adding run $run would exceed cap ($jobCap) – stopping."
    break
  fi

  # ~~ submit each chunk ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  chunkNo=0
  for listFile in "${chunks[@]}"; do
    (( ++chunkNo ))
    [[ -s "$listFile" ]] || { echo "[ERROR] Empty chunk $listFile"; continue; }

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
queue
EOS

    if condor_submit "$subFile" >/dev/null; then
      (( ++submitted ))
      vecho "[VERBOSE]   submitted chunk $chunkNo/$nChunks"
    else
      echo "[ERROR] condor_submit failed for $subFile"
    fi
  done

  vecho "[VERBOSE] Completed run $run – jobs now at $submitted"
  [[ "$mode" == "condorTest" ]] && { vecho "[VERBOSE] condorTest done."; break; }
done

echo "[INFO] Grand‑total Condor jobs submitted: $submitted"
[[ $jobCap -gt 0 ]] && echo "[INFO] Cap for this launch: $jobCap"
##############################################################################
