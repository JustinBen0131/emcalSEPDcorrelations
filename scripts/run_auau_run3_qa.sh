#!/usr/bin/env bash
##############################################################################
# run_auau_run3_qa.sh
#  argv[1]  run number
#  argv[2]  list file (≤ N lines as produced by submitter)
#  argv[3]  Condor cluster‑ID (for unique filenames)
#  argv[4]  DEST_BASE — base directory for .root output
#            (if empty, defaults to $SCRATCH/output/<run>)
##############################################################################
set -euo pipefail

# ----- user paths -----------------------------------------------------------
USER="$(id -un)"
HOME="/sphenix/u/${USER}"
SCRATCH="/sphenix/u/${USER}/scratch/TriggerAnalysis"
MYINSTALL="${HOME}/install"
MACRO_DIR="${SCRATCH}/../macros"
# ----------------------------------------------------------------------------

runNumber="$1"
fileList="$2"
clusterID="${3:-0}"
destBase="${4:-}"

# ----- environment ----------------------------------------------------------
set +u
export PGHOST=localhost
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export PGHOST=localhost
set -u
source /opt/sphenix/core/bin/setup_local.sh "$HOME/install"
# ----------------------------------------------------------------------------

# ----- decide output directory ---------------------------------------------
if [[ -n "$destBase" ]]; then
  outDir="${destBase}/${runNumber}"
else
  outDir="${SCRATCH}/output/${runNumber}"
fi
mkdir -p "$outDir"
# ----------------------------------------------------------------------------

firstFile="$(head -n1 "$fileList")"
baseName="$(basename "$firstFile")"
rootOut="${outDir}/TrigPlot_run${runNumber}_c${clusterID}_${baseName%.*}.root"

echo "[INFO] $(date)  Run=$runNumber  Files=$(wc -l < "$fileList")"
echo "[INFO] Output → $rootOut"

root -b -l -q "macro/Fun4All_getJetTrigs.C(0, \"$fileList\", \"$rootOut\")"

echo "[INFO] Completed $(date)"

