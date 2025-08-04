#!/usr/bin/env bash
##############################################################################
#  runAuAuExecutable.sh – Condor wrapper, one job per input ROOT file
##############################################################################
set -euo pipefail

export RUN_LOCATION=sphenix          # force scratch tree

inputRoot="$1"; shift
outDir="$1"  ; shift
mkdir -p "${outDir}"

MACRO="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/analyzeRun24or25auau.cpp"
buildDir="$(dirname "${MACRO}")/.aclic_build"

root -l -b -q -e "gSystem->SetBuildDir(\"${buildDir}\",kTRUE)" \
     -e ".L ${MACRO}+"                                          \
     -e "runOneQaPass(\"${inputRoot}\",\"${outDir}\")"

