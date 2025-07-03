#!/usr/bin/env bash
##############################################################################
# run_auau_run3_qa.sh
#   argv[1]  run number                    (e.g. 54128)
#   argv[2]  chunk .list                   (≤ CHUNK_SIZE lines)
#   argv[3]  Condor Cluster ID             (for unique file names)
#   argv[4]  DEST_BASE  (optional)         default → group scratch
#
#   Produces TrigPlot_run<run>_c<cluster>_<firstDST>.root in:
#       DEST_BASE/<run>/
##############################################################################
set -euo pipefail

########################  FIXED PATHS  ########################################
USER="$(id -un)"                       # still useful for local tests
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
MACRO_DIR="${PROJECT_BASE}/macro"
SCRATCH_BASE="/sphenix/u/${USER}/scratch/TriggerAnalysis"  # minimal local scratch

DEFAULT_DEST="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
################################################################################

runNumber="$1"; shift
fileList="$1"; shift
clusterID="${1:-0}"; shift
destBase="${1:-$DEFAULT_DEST}"

[[ -s "$fileList" ]] || { echo "[FATAL] Empty list file: $fileList" >&2; exit 2; }

#  -- Condor sets HOME wrong; force a sane env --------------------------------
set +u
export PGHOST=localhost
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export PGHOST=localhost
set -u
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/u/${USER}/install"
################################################################################

#  Output directory -----------------------------------------------------------
outDir="${destBase}/${runNumber}"
mkdir -p "$outDir"

firstFile="$(head -n1 "$fileList")"
baseTag="$(basename "${firstFile%.root}")"
rootOut="${outDir}/TrigPlot_run${runNumber}_c${clusterID}_${baseTag}.root"

echo "[INFO] $(date)  Run=$runNumber  Files=$(wc -l < "$fileList")"
echo "[INFO] Writing → $rootOut"

root -b -l -q "${MACRO_DIR}/Fun4All_emcalSEPDcorrelator.C(0, \"$fileList\", \"$rootOut\")"

echo "[INFO] Completed $(date)"
