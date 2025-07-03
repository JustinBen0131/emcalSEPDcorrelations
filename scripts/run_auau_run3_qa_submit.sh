#!/usr/bin/env bash
##############################################################################
# run_auau_run3_qa_submit.sh
#
#  Modes
#    local        : process first DST of first .list (quick sanity check)
#    condor       : submit every run (one Condor job per CHUNK_SIZE DSTs)
#    condorTest   : same as "condor" but stops after the first run
#    condor firstTen
#                 : keep adding whole runs until the grand‑total job count
#                   would exceed MAX_JOBS; then stop.
#
#  Project layout (fixed):
#      PROJECT_BASE = /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations
#      dst_list/    = *.list built by forceFileListSimCreation.sh
#      tmp_condor_lists/ – temporary split chunks (auto‑clean on each run)
#      log/ , stdout/ , error/  – Condor logs
##############################################################################
set -euo pipefail

########################  USER‑TUNABLE  #######################################
CHUNK_SIZE=5           # files per Condor job
MAX_JOBS=10000         # cap for "condor firstTen"
###############################################################################

#  PROJECT CONSTANTS  ---------------------------------------------------------
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
DST_LIST_DIR="${PROJECT_BASE}/dst_list"
TMP_LIST_DIR="${PROJECT_BASE}/tmp_condor_lists"
EXEC="${PROJECT_BASE}/run_auau_run3_qa.sh"          # worker
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"

LOGDIR="${PROJECT_BASE}/log"
OUTDIR="${PROJECT_BASE}/stdout"
ERRDIR="${PROJECT_BASE}/error"
mkdir -p "$TMP_LIST_DIR" "$LOGDIR" "$OUTDIR" "$ERRDIR"
# -----------------------------------------------------------------------------

########################  ARG PARSE  ##########################################
mode="${1:-}"
limitSwitch="${2:-}"

if [[ "$mode" != "local" && "$mode" != "condor" && "$mode" != "condorTest" ]]; then
  echo "Usage: $0 {local|condor|condorTest} [firstTen]" >&2
  exit 1
fi
# -----------------------------------------------------------------------------


########################  VERBOSITY HELPERS  ##################################
VERBOSE=0
[[ "$mode" == "condorTest" || ( "$mode" == "condor" && "$limitSwitch" == "firstTen" ) ]] && VERBOSE=1
vecho() { (( VERBOSE )) && echo "$*"; }
# -----------------------------------------------------------------------------


########################  GLOBAL JOB CAP  #####################################
jobCap=0
[[ "$mode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobCap=$MAX_JOBS
submitted=0
# -----------------------------------------------------------------------------

###############################################################################
#  BUILD LIST OF RUNS FROM THE FILENAMES IN dst_list/
###############################################################################
mapfile -t listFiles < <(ls "${DST_LIST_DIR}"/DST_CALO_run2auau_new_2024p007-000*.list 2>/dev/null | sort)
if (( ${#listFiles[@]} == 0 )); then
  echo "[FATAL] No .list files found in ${DST_LIST_DIR}" >&2
  exit 2
fi

# extract run numbers (00054128 …) → 54128
runs=()
for f in "${listFiles[@]}"; do
  bn=${f##*-000}; runs+=( "${bn%.list}" )
done
###############################################################################

########################  LOCAL MODE  #########################################
if [[ "$mode" == "local" ]]; then
  firstList="${listFiles[0]}"
  firstRun="${runs[0]}"
  firstDST="$(head -n1 "$firstList")"   || { echo "[ERROR] Empty $firstList" >&2; exit 3; }

  echo "[INFO] Local mode – Run=$firstRun"
  echo "[INFO] First DST  : $firstDST"

  tmpList=$(mktemp "${TMP_LIST_DIR}/local_${firstRun}_XXXX.list")
  echo "$firstDST" > "$tmpList"
  "${EXEC}" "$firstRun" "$tmpList" 0 "$CONDOR_OUT_BASE"
  rm -f "$tmpList"
  exit 0
fi
###############################################################################

########################  CONDOR / CONDORTEST  ################################
for idx in "${!runs[@]}"; do
  run=${runs[$idx]}
  masterList=${listFiles[$idx]}

  vecho "[VERBOSE] Considering run $run  (file: $(basename "$masterList"))"

  # split into N‑file chunks
  rm -f "${TMP_LIST_DIR}/run${run}_chunk_"* || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${run}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null)
  nChunks=${#chunks[@]}
  (( nChunks )) || { echo "[WARN] Empty run $run – skipping." >&2; continue; }

  # ---- global job‑cap logic -------------------------------------------------
  if (( jobCap > 0 )); then
    prospective=$((submitted + nChunks))
    vecho "[VERBOSE] Prospective total: $prospective / $jobCap"
    if (( prospective > jobCap )); then
      echo "[INFO] Job cap ($jobCap) would be exceeded – stopping at run $run."
      break
    fi
  fi
  # --------------------------------------------------------------------------

  # ---- submit every chunk ---------------------------------------------------
  chunkNo=0
  for listFile in "${chunks[@]}"; do
    ((++chunkNo))
    [[ -s "$listFile" ]] || { echo "[ERROR] Empty chunk $listFile" >&2; continue; }

    firstDST="$(head -n1 "$listFile")"
    tag="$(basename "${firstDST%.root}")"

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
      ((++submitted))
      vecho "[VERBOSE]   submitted chunk $chunkNo/$nChunks"
    else
      echo "[ERROR] condor_submit failed for $subFile" >&2
    fi
  done
  vecho "[VERBOSE] Completed run $run – jobs now at $submitted"

  [[ "$mode" == "condorTest" ]] && { vecho "[VERBOSE] condorTest done."; break; }
done

echo "[INFO] Grand‑total Condor jobs submitted: $submitted"
[[ $jobCap -gt 0 ]] && echo "[INFO] Cap for this run: $jobCap"
##############################################################################
