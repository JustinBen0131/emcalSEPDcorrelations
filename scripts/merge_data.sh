#!/usr/bin/env bash
###############################################################################
# merge_data.sh
#
# (1)  ./merge_data.sh condor  [test|firstHalf]
#      Split the list of run‑number directories under CONDOR_OUT_BASE and
#      submit one Condor job per run.  Each job executes
#
#            hadd  output/<run>/output_<run>.root  \
#                 CONDOR_OUT_BASE/<run>/*.root
#
#      Sub‑modes:
#        test       – only the *first* run is merged (smoke test)
#        firstHalf  – submit ½ of the runs
#
# (2)  ./merge_data.sh addRuns [condor]
#      Merge every existing  output_<run>.root  into  output_total.root.
#      With the optional ‘condor’ it runs as one Condor job; otherwise
#      executes locally on the login node.
###############################################################################
set -euo pipefail

###############################################################################
# User configuration
###############################################################################
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/output"
RUN_MERGED_PREFIX="output"          # → output_<run>.root
TMP_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/tmp_run_lists"
CONDOR_STDOUT="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/stdout"
CONDOR_STDERR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/error"
CONDOR_LOGDIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/log"
HADD_WRAPPER="hadd_run_condor.sh"   # helper created on‑the‑fly

REQUEST_MEMORY="2000MB"             # raise if a single run > 2 GB in memory

###############################################################################
usage() {
  echo "Usage:"
  echo "  $0 condor  [test|firstHalf]"
  echo "  $0 addRuns [condor]"
  exit 1
}

[[ $# -eq 0 || $# -gt 2 ]] && usage
MODE="$1"; SUBMODE="${2:-}"
[[ "$MODE" != "condor" && "$MODE" != "addRuns" ]] && usage

mkdir -p "$OUTPUT_DIR" "$TMP_LIST_DIR" \
         "$CONDOR_STDOUT" "$CONDOR_STDERR" "$CONDOR_LOGDIR"

###############################################################################
# Helper that *only* loads sPHENIX and runs hadd        (nounset relaxed)
###############################################################################
cat > "$HADD_WRAPPER" <<'EOS'
#!/usr/bin/env bash
set -eo pipefail
set +u
export USER="$(id -un)"
export LOGNAME="$USER"
export HOME="/sphenix/u/$USER"
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/$USER/install"
set -u

LIST="$1"          # text file with input ROOTs (absolute paths)
OUT="$2"           # output ROOT
[[ -s "$LIST" ]] || { echo "[FATAL] empty list $LIST"; exit 2; }
echo "[hadd_run_condor] merging $(wc -l < "$LIST") files  →  $OUT"
hadd -v 3 -f "$OUT" @"$LIST"
EOS
chmod +x "$HADD_WRAPPER"

###############################################################################
# 1) CONDOR – one job per run‑directory
###############################################################################
if [[ "$MODE" == "condor" ]]; then
  # Enumerate every run directory that contains at least one ROOT fragment
  mapfile -t runs < <(find "$CONDOR_OUT_BASE" -mindepth 1 -maxdepth 1 -type d \
                      -printf "%f\n" | sort)

  [[ ${#runs[@]} -eq 0 ]] && { echo "[ERROR] no run directories found"; exit 1; }

  echo "[INFO] Found ${#runs[@]} run folders under $CONDOR_OUT_BASE"

  # Optional sub‑mode throttling
  if [[ "$SUBMODE" == "test" ]]; then
    runs=( "${runs[0]}" )
  elif [[ "$SUBMODE" == "firstHalf" ]]; then
    half=$(( ${#runs[@]} / 2 ))
    runs=( "${runs[@]:0:$half}" )
  elif [[ -n "$SUBMODE" ]]; then
    usage
  fi

  # Build one submit file --------------------------------------------------
  SUB="$TMP_LIST_DIR/merge_runs.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $HADD_WRAPPER
output     = $CONDOR_STDOUT/merge.\$(Cluster).\$(Process).out
error      = $CONDOR_STDERR/merge.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory = $REQUEST_MEMORY
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
stream_output = True
stream_error  = True
EOT

  # Generate one (list + queue line) per run -------------------------------
  i=0
  for run in "${runs[@]}"; do
    ((i++))
    inDir="$CONDOR_OUT_BASE/$run"
    outFile="$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root"

    # Build list file
    list="$TMP_LIST_DIR/in_${run}.txt"
    find "$inDir" -maxdepth 1 -type f -name "*.root" | sort > "$list"
    [[ ! -s "$list" ]] && { echo "[WARN] $run has no ROOT files – skipped"; continue; }

    # Pre‑submission clean‑up of old merged output
    rm -f "$outFile"

    echo "arguments = $list $outFile" >> "$SUB"
    echo "queue"                     >> "$SUB"
  done

  condor_submit "$SUB"
  echo "[INFO] Submitted $i Condor jobs for run‑by‑run merging"
  exit 0
fi

###############################################################################
# 2) ADDRUNS – merge output_<run>.root → output_total.root
###############################################################################
# Build list of per‑run outputs that already exist
mapfile -t runFiles < <(find "$OUTPUT_DIR" -maxdepth 1 -type f \
                        -name "${RUN_MERGED_PREFIX}_*.root" | sort)

[[ ${#runFiles[@]} -eq 0 ]] && { echo "[ERROR] no run‑level ROOTs found"; exit 1; }
[[ ${#runFiles[@]} -eq 1 ]] && { echo "[INFO] Only one run – nothing to add"; exit 0; }

LIST_ALL="$TMP_LIST_DIR/all_runs.txt"
printf "%s\n" "${runFiles[@]}" > "$LIST_ALL"
FINAL="$OUTPUT_DIR/${RUN_MERGED_PREFIX}_total.root"
rm -f "$FINAL"

if [[ "$SUBMODE" == "condor" ]]; then
  SUB="$TMP_LIST_DIR/final_merge.sub"; rm -f "$SUB"
  cat > "$SUB" <<EOT
universe   = vanilla
executable = $HADD_WRAPPER
output     = $CONDOR_STDOUT/final.\$(Cluster).\$(Process).out
error      = $CONDOR_STDERR/final.\$(Cluster).\$(Process).err
log        = $CONDOR_LOGDIR/final.\$(Cluster).\$(Process).log
request_memory = $REQUEST_MEMORY
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
arguments  = $LIST_ALL $FINAL
queue
EOT
  condor_submit "$SUB"
  echo "[INFO] Submitted one Condor job to create $(basename "$FINAL")"
  exit 0
fi

echo "[INFO] Local merge → $(basename "$FINAL")"
hadd -v 3 -f "$FINAL" @"$LIST_ALL"
echo "[DONE] $(ls -lh "$FINAL")"
