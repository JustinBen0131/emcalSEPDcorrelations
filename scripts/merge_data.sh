#!/usr/bin/env bash
###############################################################################
# merge_data.sh – highly verbose Condor hadd helper
#
# ▸ Purpose
#   1.  One “hadd” job per run directory      →   MODE = condor
#   2.  Merge all per‑run ROOT files          →   MODE = addRuns
#   3.  Merge a single run locally            →   MODE = local <runNumber>
#
# ▸ New in this version
#   • Busy‑run detection (unchanged from v2)
#   • **removeOtherJobs** switch:
#       – gathers every live HTCondor job (for $USER)
#       – deletes its chunk list, partial ROOT outputs and the job itself
#       – makes the merge phase immune to corrupted, half‑written files
###############################################################################

###############################################################################
# ---- 1. Strict mode + debug plumbing ----------------------------------------
###############################################################################
set -euo pipefail
IFS=$'\n\t'
(( ${DEBUG:-0} )) && set -x

trap 'err=$?; printf "\033[0;31m%s  ✘  line %d – cmd `%s` exited %d\033[0m\n" \
                "$(date "+%F %T")" "${BASH_LINENO[0]}" "$BASH_COMMAND" "$err" >&2' ERR

###############################################################################
# ---- 2. Colour helpers + logging wrappers -----------------------------------
###############################################################################
CLR_B="\033[1;34m"; CLR_G="\033[0;32m"; CLR_Y="\033[1;33m"; CLR_R="\033[0;31m"; CLR_RST="\033[0m"
ts()   { date "+%F %T"; }
say()  { printf "${CLR_B}%s  ➜  %s${CLR_RST}\n"  "$(ts)" "$*"; }
good() { printf "${CLR_G}%s  ✔  %s${CLR_RST}\n"  "$(ts)" "$*"; }
warn() { printf "${CLR_Y}%s  ⚠  %s${CLR_RST}\n"  "$(ts)" "$*" >&2; }
fatal(){ printf "${CLR_R}%s  ✘  %s${CLR_RST}\n"  "$(ts)" "$*" >&2; exit 1; }

###############################################################################
# ---- 3. User settings -------------------------------------------------------
###############################################################################
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/output"
RUN_MERGED_PREFIX="output"          # → output_<run>.root
TMP_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/tmp_run_lists"
CONDOR_STDOUT="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/stdout"
CONDOR_STDERR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/error"
CONDOR_LOGDIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/log"
HADD_WRAPPER="hadd_run_condor.sh"
REQUEST_MEMORY="2000MB"

mkdir -p "$OUTPUT_DIR" "$TMP_LIST_DIR" "$CONDOR_STDOUT" \
         "$CONDOR_STDERR" "$CONDOR_LOGDIR"

###############################################################################
# ---- 4. Usage & CLI parsing -------------------------------------------------
###############################################################################
usage() {
  cat <<EOF
Usage:
  $0 MODE [removeOtherJobs] [QUALIFIER]

  MODE
    condor   – merge every idle run via Condor (per‑run hadd)
               QUALIFIER = test | firstHalf
    addRuns  – grand‑total merge of per‑run outputs
               QUALIFIER = condor   (do the grand‑total on Condor)
    local    – merge a single run locally
               QUALIFIER = <runNumber>

  removeOtherJobs
      (optional, may follow any MODE)
      Purge every still‑running job in your Condor queue **before**
      the merge starts.  The purge:
        • deletes the affected *.list chunk files
        • deletes partial *.root files in \$CONDOR_OUT_BASE/<run>
        • runs 'condor_rm \$USER'

Environment:
  DEBUG=1   enable shell trace & extra logging
EOF
  exit 1
}

(( $# >= 1 )) || usage
MODE=$1; shift

#  -- optional purge switch ---------------------------------------------------
PURGE=0
if [[ ${1:-} == removeOtherJobs ]]; then
    PURGE=1
    shift
fi

SUBMODE=${1:-}

case "$MODE" in
  condor)     [[ -z "$SUBMODE" || "$SUBMODE" =~ ^(test|firstHalf)$ ]] || usage ;;
  addRuns)    [[ -z "$SUBMODE" || "$SUBMODE" == condor ]]              || usage ;;
  local)      [[ -n "$SUBMODE" && "$SUBMODE" =~ ^[0-9]{5,8}$ ]]        || usage ;;
  *)          usage ;;
esac

###############################################################################
# ---- 5. Build the tiny wrapper executed inside each Condor slot ------------
###############################################################################
cat > "$HADD_WRAPPER" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
LIST=$1; OUT=$2
[[ -s $LIST ]] || { echo "[FATAL] empty list $LIST"; exit 2; }

#  sPHENIX environment --------------------------------------------------------
set +u
export PGHOST=${PGHOST:-localhost}
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u

exec 1> >(stdbuf -oL cat) 2>&1          # live stdout/err streaming
echo "[wrapper] $(wc -l <"$LIST") inputs  →  $OUT"
hadd -v -v -v -f "$OUT" @"$LIST"
EOS
chmod +x "$HADD_WRAPPER"

###############################################################################
# ---- 6. Helpers -------------------------------------------------------------
###############################################################################
safe_find() {                                  # noisy, fault‑tolerant “find”
  local dir=$1 list=$2
  say "  • scanning $dir"
  (
    set +e +o pipefail
    find "$dir" -type f -name '*.root' -print 2> >(while read -r l; do warn "    find: $l"; done) |
      sort > "$list"
  )
  (( $? == 0 && -s $list )) || { warn "    ➜ no ROOT files"; return 1; }
  (( DEBUG )) && { say "    first 10 entries:"; head -n 10 "$list" | sed 's/^/      /'; }
}

# ---- 6a. busy‑run cache -----------------------------------------------------
declare -Ag busySet=()                 # busySet[run]=1   (global)

refresh_busy_runs() {
  busySet=()
  while read -r token; do busySet["$token"]=1; done < <(
      condor_q "$USER" -af Cmd Args 2>/dev/null |
      grep -Eo '[0-9]{5,8}' | sort -u
  )
}

# ---- 6b. heavy‑duty purge ---------------------------------------------------
purge_busy_jobs() {
  refresh_busy_runs
  (( ${#busySet[@]} )) || { good "No running Condor jobs – nothing to purge"; return; }

  say  "removeOtherJobs  –  purging $(printf '%d' "${#busySet[@]}") busy run(s)"

  # 1) remove chunk *.list files ------------------------------------------------
  mapfile -t busyLists < <(
      condor_q "$USER" -af Args 2>/dev/null |
      grep -Eo '/[^[:space:]]+tmp_condor_lists[^[:space:]]+\.list' | sort -u
  )
  if (( ${#busyLists[@]} )); then
      say "  • deleting ${#busyLists[@]} temporary list(s)"
      for f in "${busyLists[@]}"; do [[ -f $f ]] && rm -f "$f" && say "      – $(basename "$f")"; done
  else
      say "  • no chunk lists to delete"
  fi

  # 2) remove partial output ROOTs --------------------------------------------
  say "  • cleaning per‑run output directories"
  for run in "${!busySet[@]}"; do
      outDir="$CONDOR_OUT_BASE/$run"
      [[ -d $outDir ]] || continue
      n=$(find "$outDir" -type f -name '*.root' | wc -l)
      [[ $n -gt 0 ]] && { find "$outDir" -type f -name '*.root' -delete; say "      – run $run  ($n file(s) purged)"; }
  done

  # 3) remove Condor jobs ------------------------------------------------------
  say "  • condor_rm $USER"
  if condor_rm "$USER" >/dev/null 2>&1; then
      good "All Condor jobs removed"
  else
      warn "condor_rm failed – please check manually"
  fi
}

###############################################################################
# ---- 7. Optional pre‑merge purge -------------------------------------------
###############################################################################
if (( PURGE )); then purge_busy_jobs; fi
refresh_busy_runs                     # always refresh cache afterwards

###############################################################################
# ---- 8.  PER‑RUN MERGE (MODE = condor) --------------------------------------
###############################################################################
if [[ $MODE == condor ]]; then
  # housekeeping --------------------------------------------------------------
  say "Cleaning previous Condor text outputs"
  for d in "$CONDOR_STDOUT" "$CONDOR_STDERR" "$CONDOR_LOGDIR"; do
      [[ -d $d ]] && find "$d" -type f -delete
  done
  say "Removing stale per‑run ROOT files from $OUTPUT_DIR"
  find "$OUTPUT_DIR" -maxdepth 1 -type f -name "${RUN_MERGED_PREFIX}_????????.root" -delete

  # Step 1 – enumerate run directories ----------------------------------------
  say "Step 1 – enumerating run directories under $CONDOR_OUT_BASE"
  mapfile -t runs < <(find "$CONDOR_OUT_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

  # Step 2 – drop busy runs ----------------------------------------------------
  say "Step 2 – filtering out active Condor runs"
  runs=( $(for r in "${runs[@]}"; do [[ -z ${busySet[$r]+x} ]] && echo "$r"; done) )
  (( ${#busySet[@]} )) && say "    active: $(printf '%s ' "${!busySet[@]}")" || say "    none"
  (( ${#runs[@]} )) || { good "Nothing idle to merge"; exit 0; }

  say "Step 3 – ${#runs[@]} idle run(s) will be processed"
  case "$SUBMODE" in
      test)      runs=( "${runs[0]}" ); warn "TEST mode – only ${runs[0]}" ;;
      firstHalf) half=$(( ${#runs[@]} / 2 )); runs=( "${runs[@]:0:$half}" ); warn "FIRST‑HALF mode – $half run(s)" ;;
  esac

  # Condor submit description --------------------------------------------------
  SUB="$TMP_LIST_DIR/merge_runs.sub"; : > "$SUB"
  cat >>"$SUB" <<EOT
universe        = vanilla
executable      = $HADD_WRAPPER
output          = $CONDOR_STDOUT/merge.\$(Cluster).\$(Process).out
error           = $CONDOR_STDERR/merge.\$(Cluster).\$(Process).err
log             = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory  = $REQUEST_MEMORY
getenv          = True
stream_output   = True
stream_error    = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
EOT

  jobCnt=0
  for run in "${runs[@]}"; do
      say "── Run $run ───────────────────────────────────────────────────────"
      inDir=$CONDOR_OUT_BASE/$run
      list=$TMP_LIST_DIR/in_${run}.txt
      outFile=$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root
      safe_find "$inDir" "$list" || continue
      nFiles=$(wc -l <"$list")
      say "    will merge $nFiles file(s) → $outFile"
      printf 'arguments = %s %s\nqueue 1\n' "$list" "$outFile" >> "$SUB"
      (( ++jobCnt ))
  done

  (( jobCnt )) || fatal "Zero jobs created – aborting"
  (( DEBUG )) && { say "Submit description (first 10 lines):"; head -n 10 "$SUB" | sed 's/^/    /'; }
  say "Step 4 – submitting $jobCnt job(s) to Condor"
  condor_submit "$SUB"
  good "Condor submission finished"
  exit 0
fi

###############################################################################
# ---- 9.  SINGLE‑RUN LOCAL MERGE --------------------------------------------
###############################################################################
if [[ $MODE == local ]]; then
  run="$SUBMODE"
  [[ -n ${busySet[$run]+x} ]] && { warn "Run $run is still active in Condor – skipping local merge"; exit 0; }

  say "Local merge for run ${run}"
  inDir=$CONDOR_OUT_BASE/$run
  list=$TMP_LIST_DIR/in_${run}.txt
  outFile=$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root
  safe_find "$inDir" "$list" || fatal "No ROOT files found for run ${run}"
  [[ -f $outFile ]] && { say "  • removing old $outFile"; rm -f "$outFile"; }
  "$HADD_WRAPPER" "$list" "$outFile"
  good "Merge finished – $(ls -lh "$outFile")"
  exit 0
fi

###############################################################################
# ---- 10. GRAND‑TOTAL MERGE --------------------------------------------------
###############################################################################
say "Grand‑total stage – collecting per‑run outputs in $OUTPUT_DIR"
mapfile -t runFiles < <(
  find "$OUTPUT_DIR" -maxdepth 1 -type f -name "${RUN_MERGED_PREFIX}_*.root" |
  while read -r f; do
      [[ $f =~ _([0-9]{5,8})\.root$ ]] || continue
      r=${BASH_REMATCH[1]}
      [[ -z ${busySet[$r]+x} ]] && echo "$f"
  done | sort
)
(( ${#runFiles[@]} < 2 )) && { good "Nothing to add"; exit 0; }

LIST_ALL="$TMP_LIST_DIR/all_runs.txt"; printf '%s\n' "${runFiles[@]}" > "$LIST_ALL"
FINAL="$OUTPUT_DIR/${RUN_MERGED_PREFIX}_total.root"; rm -f "$FINAL"

if [[ $SUBMODE == condor ]]; then
  SUB="$TMP_LIST_DIR/final_merge.sub"; : >"$SUB"
  cat >>"$SUB" <<EOT
universe = vanilla
executable = $HADD_WRAPPER
output  = $CONDOR_STDOUT/final.\$(Cluster).\$(Process).out
error   = $CONDOR_STDERR/final.\$(Cluster).\$(Process).err
log     = $CONDOR_LOGDIR/final.\$(Cluster).\$(Process).log
request_memory  = $REQUEST_MEMORY
getenv          = True
stream_output   = True
stream_error    = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
arguments = $LIST_ALL $FINAL
queue 1
EOT
  (( DEBUG )) && { say "Submitting grand‑total job (DEBUG preview):"; sed 's/^/    /' "$SUB"; }
  condor_submit "$SUB"
  good "Grand‑total Condor job submitted"
else
  say "Running grand‑total merge locally → $(basename "$FINAL")"
  hadd -v 3 -f "$FINAL" @"$LIST_ALL"
  good "DONE – $(ls -lh "$FINAL")"
fi
