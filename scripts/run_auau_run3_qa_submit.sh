#!/usr/bin/env bash
##############################################################################
# run_auau_run3_qa_submit.sh
#
#  Modes
#  ─────
#    local        : run Fun4All_getJetTrigs.C once on the first DST file
#    condor       : submit 1 Condor job per <CHUNK_SIZE> files for **every** run
#    condorTest   : like "condor" but ONLY the first run is processed
#    condor firstTen
#                 : process successive whole runs until doing another run
#                   would raise the total number of jobs above <MAX_JOBS>.
#
#  Generalised submission policy (independent of CHUNK_SIZE / MAX_JOBS)
#  --------------------------------------------------------------------
#    1.  Chunking rule    – Lists are split into fixed‑size groups of N files
#                           (final group may contain < N). 1 group → 1 job.
#    2.  Atomic‑run rule  – A run is submitted in its entirety or not at all.
#    3.  Global‑limit     – Stop before any run that would breach MAX_JOBS.
#    4.  Test‑mode rule   – "condorTest" ignores the limit and submits every
#                           chunk for the first run only.
##############################################################################

set -euo pipefail

# =========  DEBUG / ERROR TRAPS  ============================================
trap 'echo "[FATAL] Script died on line $LINENO while executing: $BASH_COMMAND" >&2' ERR
[[ -n "${BASH_XTRACEFD:-}" ]] || { exec 9>&2; export BASH_XTRACEFD=9; }  # safe for set -x
DEBUG=${DEBUG:-0}        # export DEBUG=1 for shell trace
(( DEBUG )) && set -x
# ============================================================================

# ───────────── USER‑TUNABLE CONSTANTS ───────────────────────────────────────
CHUNK_SIZE=5            # N  ← files per Condor job
MAX_JOBS=10000          # Q  ← global cap for "condor firstTen"
# ────────────────────────────────────────────────────────────────────────────

mode="${1:-}"
limitSwitch="${2:-}"          # empty or 'firstTen'

if [[ "$mode" != "local" && "$mode" != "condor" && "$mode" != "condorTest" ]]; then
  echo "Usage: $0  {local | condor | condorTest}  [firstTen]" >&2
  exit 1
fi

# ───────────── VERBOSITY ────────────────────────────────────────────────────
VERBOSE=0
if [[ "$mode" == "condorTest" || ( "$mode" == "condor" && "$limitSwitch" == "firstTen" ) ]]; then
  VERBOSE=1
fi
vecho() { (( VERBOSE )) && echo "$@"; }

# ───────────── GLOBAL LIMIT SELECTION ───────────────────────────────────────
jobLimit=0                                # plain 'condor' → unlimited
[[ "$mode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobLimit=$MAX_JOBS
# ────────────────────────────────────────────────────────────────────────────

# ───────────── PATH SET‑UP (unchanged) ──────────────────────────────────────
USER="$(id -un)"
HOME="/sphenix/u/${USER}"
SCRATCH="${HOME}/scratch/TriggerAnalysis"
EXEC="${SCRATCH}/RunTriggerPlotter_Condor.sh"

CONDOR_BASE="/sphenix/tg/tg01/bulk/jbennett/TriggerAna"   # ROOT output
LOGDIR="${SCRATCH}/log"
OUTDIR="${SCRATCH}/stdout"
ERRDIR="${SCRATCH}/error"
mkdir -p "$LOGDIR" "$OUTDIR" "$ERRDIR"

TMP_LIST_DIR="${SCRATCH}/condor_lists"
MACRO_DIR="macro"
mkdir -p "$TMP_LIST_DIR"

vecho "[VERBOSE] Mode             : $mode ($limitSwitch)"
vecho "[VERBOSE] CHUNK_SIZE (N)   : $CHUNK_SIZE"
vecho "[VERBOSE] jobLimit (Q)     : $jobLimit   (0 → ∞ means no cap)"
vecho "[VERBOSE] Condor base path : $CONDOR_BASE"
vecho "[VERBOSE] Scratch path     : $SCRATCH"

runFile="${SCRATCH}/Final_RunNumbers_After_All_Cuts.txt"
[[ -s "$runFile" ]] || { echo "[ERROR] Run‑number file missing: $runFile" >&2; exit 1; }
readarray -t runs < "$runFile"
# ────────────────────────────────────────────────────────────────────────────

##############################################################################
#  LOCAL MODE (unchanged)
##############################################################################
if [[ "$mode" == "local" ]]; then
  firstRun="${runs[0]}"
  dstMaster="${SCRATCH}/dst_list/dst_jet_run2pp-000${firstRun}.list"
  [[ -s "$dstMaster" ]] || { echo "[ERROR] $dstMaster missing" >&2; exit 1; }

  firstFile="$(head -n1 "$dstMaster")"
  echo "[INFO] Local mode – first run: $firstRun"
  echo "[INFO] Using DST file: $firstFile"

  tmpList=$(mktemp "${TMP_LIST_DIR}/local_${firstRun}_XXXX.list")
  echo "$firstFile" > "$tmpList"

  outDir="${SCRATCH}/output_local/${firstRun}"
  mkdir -p "$outDir"
  rootOut="${outDir}/TrigPlot_local_${firstRun}.root"

  # environment
  set +u
  export PGHOST=localhost
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  export PGHOST=localhost
  set -u
  source /opt/sphenix/core/bin/setup_local.sh "$HOME/install"
  export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:"$HOME/install/include"

  root -b -l -q "macro/Fun4All_getJetTrigs.C(0, \"$tmpList\", \"$rootOut\")"
  echo "[INFO] Local run finished. Output → $rootOut"
  rm -f "$tmpList"
  exit 0
fi

##############################################################################
#  CONDOR / CONDORTEST MODE
##############################################################################
submitted=0
for run in "${runs[@]}"; do
  vecho "[VERBOSE] --------------"
  vecho "[VERBOSE] Considering run $run"

  masterList="${SCRATCH}/dst_list/dst_jet_run2pp-000${run}.list"
  if [[ ! -s "$masterList" ]]; then
    echo "[WARN] No list for run $run – skipping." >&2
    continue
  fi

  # split into chunks of CHUNK_SIZE lines
  rm -f "${TMP_LIST_DIR}/run${run}_chunk_"* || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${run}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${run}_chunk_"* 2>/dev/null)
  nChunks=${#chunks[@]}
  vecho "[VERBOSE] Chunks generated: $nChunks (size $CHUNK_SIZE)"
  if (( nChunks == 0 )); then
    echo "[ERROR] split produced zero chunks for run $run" >&2
    exit 2
  fi

  # ---- GLOBAL LIMIT CHECK --------------------------------------------------
  if (( jobLimit > 0 )); then
    prospective=$((submitted + nChunks))
    vecho "[VERBOSE] Prospective total jobs: $prospective / $jobLimit"
    if (( prospective > jobLimit )); then
      echo "[INFO] Job cap ($jobLimit) would be exceeded by adding run $run – stopping."
      break
    fi
  fi
  # --------------------------------------------------------------------------

  # ---- SUBMIT ALL CHUNKS FOR THIS RUN --------------------------------------
  chunkIdx=0
  for listFile in "${chunks[@]}"; do
    ((++chunkIdx))
    if [[ ! -s "$listFile" ]]; then
      echo "[ERROR] Empty chunk file: $listFile (index $chunkIdx)" >&2
      continue
    fi

    firstDST="$(head -n1 "$listFile")" || { echo "[ERROR] Cannot read $listFile" >&2; exit 3; }
    baseTag="$(basename "$firstDST" .root)"
    logFile="${LOGDIR}/${baseTag}.log"
    outFile="${OUTDIR}/${baseTag}.out"
    errFile="${ERRDIR}/${baseTag}.err"
    subFile="${SCRATCH}/TrigPlot_${baseTag}.sub"

    vecho "[VERBOSE]   [chunk $chunkIdx/$nChunks] Job tag : $baseTag"

    cat > "$subFile" <<EOL
universe      = vanilla
executable    = $EXEC
arguments     = $run $listFile \$(Cluster) $CONDOR_BASE
log           = $logFile
output        = $outFile
error         = $errFile
request_memory= 1000MB
queue
EOL

    if condor_submit "$subFile" >/dev/null; then
      ((++submitted))
    else
      echo "[ERROR] condor_submit failed for $subFile" >&2
    fi
  done
  vecho "[VERBOSE] Submitted $nChunks jobs for run $run"

  # ---- TEST‑MODE GUARD -----------------------------------------------------
  if [[ "$mode" == "condorTest" ]]; then
    vecho "[VERBOSE] condorTest complete – processed first run only."
    break
  fi
  # --------------------------------------------------------------------------
done

echo "[INFO] Grand total jobs submitted: $submitted"
if (( jobLimit > 0 )); then
  echo "[INFO] Job‑cap mode active – cap = $jobLimit."
fi

