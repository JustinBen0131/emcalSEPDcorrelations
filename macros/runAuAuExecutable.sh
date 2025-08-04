#!/usr/bin/env bash
##############################################################################
#  runAuAuExecutable.sh – Condor wrapper, one job per input ROOT file
##############################################################################
set -euo pipefail

export RUN_LOCATION=sphenix          # force scratch tree

inputRoot="$1"; shift
outDir="$1"  ; shift
mkdir -p "${outDir}"

MACRO="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/macros/analyzeRun24or25auau.cpp"
buildDir="$(dirname "${MACRO}")/.aclic_build"

# ────────────────────────────────────────────────────────────────────────────
#  restore a usable environment for ROOT
# ────────────────────────────────────────────────────────────────────────────
USER="$(id -un)"

set +u
export PGHOST=localhost
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export PGHOST=localhost
set -u
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/u/${USER}/install"

: "${HOME:=/sphenix/u/${USER}}"
export HOME
export ROOTENV_NO_HOME=1
# ────────────────────────────────────────────────────────────────────────────

root -l -b -q -e "gSystem->SetBuildDir(\"${buildDir}\",kTRUE)" \
     -e ".L ${MACRO}+"                                          \
     -e "runOneQaPass(\"${inputRoot}\",\"${outDir}\")"
