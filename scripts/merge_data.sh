#!/usr/bin/env bash
###############################################################################
# merge_data.sh – highly verbose Condor hadd helper
# Author: <you>
# Usage:
#   DEBUG=1 ./merge_data.sh condor [test|firstHalf]
#   DEBUG=1 ./merge_data.sh addRuns [condor]
###############################################################################

###############################################################################
# ---- 1. Strict mode + debug plumbing ----------------------------------------
###############################################################################
set -euo pipefail
IFS=$'\n\t'

(( ${DEBUG:-0} )) && set -x        # shell‑level trace if you want it

# global error hook – fires on *every* non‑zero exit status
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
# ---- 4. Usage guard ---------------------------------------------------------
###############################################################################
usage() {
  cat <<EOF
Usage:
  $0 condor [test|firstHalf]   # one Condor job per run
  $0 addRuns [condor]          # hadd run‑level outputs → total
Environment:
  DEBUG=1   enable shell trace & extra logging
EOF
  exit 1
}
[[ $# -lt 1 || $# -gt 2 ]] && usage
MODE=$1; SUBMODE=${2:-}
[[ $MODE != condor && $MODE != addRuns ]] && usage

###############################################################################
# ---- 5. Build the tiny wrapper executed inside each Condor slot ------------
###############################################################################
cat > "$HADD_WRAPPER" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
LIST=$1; OUT=$2
[[ -s $LIST ]] || { echo "[FATAL] empty list $LIST"; exit 2; }

# ── temporarily relax nounset so the sPHENIX env script can
#    reference variables that might be unset ──
set +u
export PGHOST=${PGHOST:-localhost}
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u

# flush every line so Condor can stream it back live
exec 1> >(stdbuf -oL cat) 2>&1

echo "[wrapper] $(wc -l <"$LIST") inputs  →  $OUT"
hadd -v -v -v -f "$OUT" @"$LIST"
EOS
chmod +x "$HADD_WRAPPER"

###############################################################################
# ---- 6. Helper: safe, noisy find --------------------------------------------
###############################################################################
safe_find() {
  local dir=$1 list_file=$2
  say "  • scanning $dir"
  # run find in a subshell with its own error handling
  (
    set +e +o pipefail
    find "$dir" -type f -name '*.root' -print 2> >(while read -r l; do warn "    find: $l"; done) |
      sort > "$list_file"
  )
  local rc=$?
  if (( rc != 0 )); then
      warn "    find exited $rc (ignored)"
  fi
  if [[ ! -s $list_file ]]; then
      warn "    ➜ NO root files found"
      return 1
  fi
  if (( DEBUG )); then
      say "    first 10 entries:"; head -n 10 "$list_file" | sed 's/^/      /'
  fi
  return 0
}

###############################################################################
# ---- 7.  PER‑RUN MERGE ------------------------------------------------------
###############################################################################
if [[ $MODE == condor ]]; then
  # ─────────────────────────────────────────────────────────────
  #  pre‑submission clean‑up
  #    • purge old stdout / stderr / log text files
  #    • remove any previous run‑merged ROOT files
  # ─────────────────────────────────────────────────────────────
  say "Cleaning previous Condor text outputs"
  for d in "$CONDOR_STDOUT" "$CONDOR_STDERR" "$CONDOR_LOGDIR"; do
      [[ -d $d ]] && find "$d" -type f -delete
  done

  say "Removing stale per‑run ROOT files from $OUTPUT_DIR"
  find "$OUTPUT_DIR" -maxdepth 1 -type f \
       -name "${RUN_MERGED_PREFIX}_????????.root" -delete

  say "Step 1 – enumerating run directories under $CONDOR_OUT_BASE"
  mapfile -t runs < <(find "$CONDOR_OUT_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

  say "Step 2 – filtering out active Condor runs"
  mapfile -t busy < <(condor_q "$USER" -af Cmd Args 2>/dev/null | grep -Eo '[0-9]{8}' | sort -u)
  if (( ${#busy[@]} )); then
      say "    active: ${busy[*]}"
  else
      say "    none"
  fi
  declare -A busySet; for r in "${busy[@]}";  do busySet[$r]=1; done
  runs=( $(for r in "${runs[@]}"; do [[ -z ${busySet[$r]+x} ]] && echo "$r"; done) )
  (( ${#runs[@]} )) || { good "Nothing idle to merge"; exit 0; }
  say "Step 3 – ${#runs[@]} idle run(s) will be processed"

  case $SUBMODE in
      test)      runs=( "${runs[0]}" ); warn "TEST mode – only ${runs[0]}" ;;
      firstHalf) half=$(( ${#runs[@]}/2 )); runs=( "${runs[@]:0:$half}" ); warn "FIRST‑HALF mode – $half run(s)" ;;
      "")        ;;
      *)         usage ;;
  esac

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

      if ! safe_find "$inDir" "$list"; then
          continue
      fi

      nFiles=$(wc -l <"$list")
      say "    will merge $nFiles file(s) → $outFile"

      printf 'arguments = %s %s\nqueue 1\n' "$list" "$outFile" >> "$SUB"
      (( ++jobCnt ))
  done

  (( jobCnt )) || fatal "Zero jobs created – aborting"

  if (( DEBUG )); then
      say "Submit description (first 10 lines):"
      head -n 10 "$SUB" | sed 's/^/    /'
  fi

  say "Step 4 – submitting $jobCnt job(s) to Condor"
  condor_submit "$SUB"
  good "Condor submission finished"
  exit 0
fi

###############################################################################
# ---- 8.  GRAND‑TOTAL MERGE --------------------------------------------------
###############################################################################
say "Grand‑total stage – collecting per‑run outputs in $OUTPUT_DIR"
mapfile -t runFiles < <(find "$OUTPUT_DIR" -maxdepth 1 -type f -name "${RUN_MERGED_PREFIX}_*.root" | sort)
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
    (( DEBUG )) && { say "Submitting grand‑total job (DEBUG preview):"; cat "$SUB" | sed 's/^/    /'; }
    condor_submit "$SUB"
    good "Grand‑total Condor job submitted"
else
    say "Running grand‑total merge locally → $(basename "$FINAL")"
    hadd -v 3 -f "$FINAL" @"$LIST_ALL"
    good "DONE – $(ls -lh "$FINAL")"
fi
