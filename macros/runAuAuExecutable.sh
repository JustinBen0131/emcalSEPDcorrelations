#!/usr/bin/env bash
##############################################################################
#  runAuAuExecutable.sh – Condor wrapper, one job per input ROOT file
##############################################################################
set -euo pipefail

export RUN_LOCATION=sphenix          # force scratch tree

# ------------------------------------------------------------
#  Detect “final” mode  (no arguments OR first arg == --final)
# ------------------------------------------------------------
if [[ $# -eq 0 || ${1:-} == "--final" ]]; then
    finalMode="true"
else
    finalMode="false"
fi

if [[ "${finalMode}" == "false" ]]; then
    inputRoot="$1"; shift
    outDir="$1"  ; shift
    mkdir -p "${outDir}"
fi

MACRO="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/macros/analyzeRun24or25auau.cpp"

# one unique build folder per job (falls back to /tmp if the env‑var is absent)
buildDir="${_CONDOR_SCRATCH_DIR:-/tmp}/aclic_build_$$"
mkdir -p "${buildDir}"

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

if [[ "${finalMode}" == "true" ]]; then
    export COMBINED_ONLY=1          # ← tell the C++ macro to skip per‑run loops
    root -l -b -q \
         -e "gSystem->SetBuildDir(\"${buildDir}\",kTRUE)" \
         -e ".L ${MACRO}+" \
         -e "analyzeRun24or25auau(false,-1)"
else
    root -l -b -q \
         -e "gSystem->SetBuildDir(\"${buildDir}\",kTRUE)" \
         -e ".L ${MACRO}+" \
         -e "runOneQaPass(\"${inputRoot}\",\"${outDir}\")"
fi
