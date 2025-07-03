#!/usr/bin/env bash
# forceFileListSimCreation.sh
#
# Purpose:
#   1) Build a G4Hits list file containing every .root file in:
#        /sphenix/lustre01/sphnxpro/mdc2/single_particle/g4hits/run0024/gamma_pt_200_40000
#      in ascending (sequential numeric) order.
#   2) Build a DST_CALO_CLUSTER list file containing every .root file in:
#        /sphenix/lustre01/sphnxpro/mdc2/single_particle/nopileup/calocluster/run0024/gamma_pt_200_40000
#      also in ascending order.
#   3) Write these two lists to:
#        /sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run24_type14_gamma_pt_200_40000
#      named G4Hits.list and DST_CALO_CLUSTER.list respectively.
#
# Usage:
#   ./forceFileListSimCreation.sh
#
# Note:
#   - We now use `find ... -name '*.root' | sort -V` instead of `ls -1v` to avoid
#     the "Argument list too long" error when many files are present.
#   - Make sure you have permission to read from the source directories 
#     and write to the output directory.

set -euo pipefail

##############################################################################
# PATHS
##############################################################################
G4HITS_SRC="/sphenix/lustre01/sphnxpro/mdc2/single_particle/g4hits/run0024/gamma_pt_200_40000"
CALO_SRC="/sphenix/lustre01/sphnxpro/mdc2/single_particle/nopileup/calocluster/run0024/gamma_pt_200_40000"

OUTDIR="/sphenix/u/patsfan753/scratch/PDCrun24pp/simListFiles/run24_type14_gamma_pt_200_40000"

G4HITS_LIST="${OUTDIR}/G4Hits.list"
CALO_LIST="${OUTDIR}/DST_CALO_CLUSTER.list"

##############################################################################
# CREATE OUTPUT DIRECTORY
##############################################################################
echo "[INFO] Creating output directory if not present: ${OUTDIR}"
mkdir -p "${OUTDIR}"

##############################################################################
# BUILD G4HITS LIST
##############################################################################
echo "[INFO] Building G4Hits list from: ${G4HITS_SRC}"
if [[ -d "${G4HITS_SRC}" ]]; then
  find "${G4HITS_SRC}" -type f -name '*.root' \
    | sort -V \
    > "${G4HITS_LIST}"

  echo "[INFO] Wrote $(wc -l < "${G4HITS_LIST}") lines to ${G4HITS_LIST}"
  echo "       (Sample entries from G4Hits.list)"
  head -n 5 "${G4HITS_LIST}"
else
  echo "[ERROR] Directory not found: ${G4HITS_SRC}"
  exit 1
fi
echo "------------------------------------------------------------------"

##############################################################################
# BUILD DST_CALO_CLUSTER LIST
##############################################################################
echo "[INFO] Building DST_CALO_CLUSTER list from: ${CALO_SRC}"
if [[ -d "${CALO_SRC}" ]]; then
  find "${CALO_SRC}" -type f -name '*.root' \
    | sort -V \
    > "${CALO_LIST}"

  echo "[INFO] Wrote $(wc -l < "${CALO_LIST}") lines to ${CALO_LIST}"
  echo "       (Sample entries from DST_CALO_CLUSTER.list)"
  head -n 5 "${CALO_LIST}"
else
  echo "[ERROR] Directory not found: ${CALO_SRC}"
  exit 1
fi
echo "------------------------------------------------------------------"

##############################################################################
# DONE
##############################################################################
echo "[DONE] File lists created:"
echo "  1) ${G4HITS_LIST}"
echo "  2) ${CALO_LIST}"
echo "------------------------------------------------------------------"


