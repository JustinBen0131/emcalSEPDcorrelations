#!/usr/bin/env bash
##############################################################################
# run_auau_run3_qa.sh
##############################################################################
set -euo pipefail

########################  FIXED PATHS  ########################################
USER="$(id -un)"                       # still useful for local tests
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
MACRO_DIR="${PROJECT_BASE}/macros"
SCRATCH_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"  # minimal local scratch

DEFAULT_DEST="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
################################################################################

runNumber="$1"; shift
fileList="$1"; shift
tag="$1"; shift
destBase="${1:-$DEFAULT_DEST}"; shift
evtMax="${1:-0}"

[[ -s "$fileList" ]] || { echo "[FATAL] Empty list file: $fileList" >&2; exit 2; }

#  -- Condor sets HOME wrong; force a sane env --------------------------------
set +u
export PGHOST=localhost
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export PGHOST=localhost
set -u
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/u/${USER}/install"

# ---- ROOT needs $HOME for $HOME/.root.mimes; define it if Condor wiped it ---
: "${HOME:=/sphenix/u/${USER}}"
export HOME
export ROOTENV_NO_HOME=1          # let ROOT skip ~/.root* if the file is absent
################################################################################


#  Output directory -----------------------------------------------------------
# If the path we received (`destBase`) already ends with the run number
# we use it as‑is; otherwise we append the run‑number folder once.
if [[ "${destBase##*/}" == "$runNumber" ]]; then
  outDir="${destBase}"
else
  outDir="${destBase}/${runNumber}"
fi
mkdir -p "$outDir"

firstFile="$(head -n1 "$fileList")"
baseTag="$(basename "${firstFile%.root}")"
rootOut="${outDir}/${tag}.root"

echo "[INFO] $(date)  Run=$runNumber  Files=$(wc -l < "$fileList")  (evtMax=$evtMax)"
echo "[INFO] Writing → $rootOut"

root -b -l -q \
    "${MACRO_DIR}/Fun4All_emcalSEPDcorrelator.C(${evtMax},\"${fileList}\",\"${rootOut}\",false)"

echo "[INFO] Completed $(date)"
