#!/usr/bin/env bash
##############################################################################
#  run_auau_run3_qa_submit.sh
#
#  Purpose
#  ───────
#  Submit – or locally dry-run – the sPHENIX **EMCAL × sEPD × MBD QA
#  correlation analysis** for Run-24/25 Au+Au data.  Jobs are either executed
#  directly on the login node (**local**) or dispatched to the Lab’s HTCondor
#  pool (**condor / condorTest**).
#
#  Quick synopsis
#  ──────────────
#    ./run_auau_run3_qa_submit.sh  <DATASET>  [QUALIFIER]  <MODE>  [OPTION] …
#
#    <DATASET>      run24auau      – uses Run-24 **CALO-DST** lists
#                   run25auau      – uses Run-25 JET-family lists *or*
#                                    Run-3 CALOFITTING lists (see QUALIFIER)
#
#    [QUALIFIER]    dstjet         – (default)  Run-25 JET DSTs
#                   dstjetcalo     – Run-25 JETCALO DSTs
#                   caloFitting    – **Run-3 CALOFITTING** sample
#
#    <MODE>         local          – run one chunk interactively for a sanity
#                                    check (prints the command it executes)
#                   condorTest     – submit *one* run to Condor for a smoke
#                                    test (verbose output)
#                   condor         – full Condor submission
#                   splitRunList   – utility: split a master run list into
#                                    manageable “round-n” segments
#
#    [OPTION]       firstTen       – with *condor*: cap the launch to 10 chunks
#                   round <N>      – with *condor*: use segment file # N
#                   <runID>        – with *local*: which run to test
#
#  Examples
#  ────────
#  • Run-24 quick local test (first DST only)
#      ./run_auau_run3_qa_submit.sh run24auau local
#
#  • Submit *all* Run-25 JET jobs to Condor (quiet mode)
#      ./run_auau_run3_qa_submit.sh run25auau condor
#
#  • Submit the first 10 CALOFITTING chunks for a smoke test
#      ./run_auau_run3_qa_submit.sh run25auau caloFitting condor firstTen
#
#  • Split a golden-run list into segment files of ≤10 000 Condor jobs each
#      ./run_auau_run3_qa_submit.sh run25auau splitRunList run25GoldenRuns.txt
#
#  How the qualifier logic works
#  ─────────────────────────────
#  ┌──────────┬──────────────────────────────────────────────────────────────┐
#  │ DATASET  │ 1st extra token → result                                    │
#  ├──────────┼──────────────────────────────────────────────────────────────┤
#  │ run24…   │ (none)                → CALO lists (Run-24 Y2Calib)         │
#  │ run25…   │ dstjet | dstjetcalo   → JET / JETCALO lists (Run-25)        │
#  │ run25…   │ caloFitting           → CALOFITTING lists (Run-3) *and*     │
#  │          │                        golden-run selector switches to      │
#  │          │                        **run3GoldenRuns.txt**               │
#  └──────────┴──────────────────────────────────────────────────────────────┘
#
#  Modes in detail
#  ───────────────
#  • local          : Runs a single chunk (first two DST files of the run)
#                     directly, useful for debugger/interactive gdb.
#  • condorTest     : Submits exactly one Condor job (first chunk of the
#                     first run) and prints progress messages.
#  • condor         : Walks every *.list file (or the supplied run list) and
#                     submits one Condor job per CHUNK_SIZE files.  Optional
#                     “firstTen” hard-caps the launch at 10 jobs.
#  • splitRunList   : Utility helper.  Takes a plain-text list of run numbers
#                     and produces runSegment_<DATASET>_<n>.txt files that each
#                     expand to ≤ MAX_JOBS Condor jobs.
#
#  Environment / paths
#  ───────────────────
#  PROJECT_BASE       = /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations
#  DST_LIST_DIR       = \$PROJECT_BASE/dst_list          (input *.list files)
#  TMP_LIST_DIR       = \$PROJECT_BASE/tmp_condor_lists  (auto-generated)
#  run_auau_run3_qa.sh – executable Condor wrapper (called per chunk)
#
##############################################################################

set -euo pipefail
shopt -s extglob
IFS=$'\n\t'

########################  COLOUR & LOG HELPERS  ###############################
ESC=$'\e['
CLR_R=${ESC}0\;31m ; CLR_G=${ESC}0\;32m ; CLR_Y=${ESC}1\;33m
CLR_B=${ESC}1\;34m ; CLR_RST=${ESC}0m
say()   { printf "${CLR_B}➜${CLR_RST} %s\n" "$*"; }
good()  { printf "${CLR_G}%s${CLR_RST}\n"   "$*"; }
warn()  { printf "${CLR_Y}⚠ %s${CLR_RST}\n" "$*" >&2; }
fatal() { printf "${CLR_R}✘ %s${CLR_RST}\n" "$*" >&2; exit 1; }
trap 'fatal "Script aborted (line $LINENO)"' ERR

##############################################################################
# 0. DATA‑SET & SPECIAL‑MODE SELECTION
##############################################################################
DATASET=${1:-run24auau}        # run24auau | run25auau
shift || true                  # always consume at least one token

# ---------------------------------------------------------------------------
# Recognise an *optional* special‑mode keyword *before* ordinary dst‑type:
#   run25auau caloFitting …
# This switches the submitter to CALOFITTING‑specific list handling.
# ---------------------------------------------------------------------------
CALOFIT=0
if [[ "$DATASET" == run25auau && "${1:-}" == caloFitting ]]; then
  CALOFIT=1
  shift                         # consume "caloFitting"
fi

# ---------------------------------------------------------------------------
# dst‑type selection for the standard JET* modes (only if CALOFIT == 0)
#   • dstjet      → DST_JET      (default)
#   • dstjetcalo  → DST_JETCALO
# ---------------------------------------------------------------------------
DSTTYPE=dstjet
if (( ! CALOFIT )) && [[ "$DATASET" == run25auau && "${1:-}" =~ ^(dstjet|dstjetcalo)$ ]]; then
  DSTTYPE=${1,,}                # lower‑case
  shift                         # consume the dst‑type token
fi
# ---------------------------------------------------------------------------

case "$DATASET" in
  run24auau)
    FILE_PREFIX="DST_CALO_run2auau_new_2024p007"
    LIST_PATTERN="${FILE_PREFIX}-000*.list"
    LIST_FMT="${FILE_PREFIX}-000%05d.list"
    PAD_FMT="%05d"
    ;;
  run25auau)
    if (( CALOFIT )); then
      FILE_PREFIX="DST_CALOFITTING_run3auau_new_newcdbtag_v006"
      LIST_PATTERN="${FILE_PREFIX}-*.list"          # no 000‑subdir for CALOFIT
      LIST_FMT="${FILE_PREFIX}-%08d.list"
      PAD_FMT="%08d"
    else
      case "$DSTTYPE" in
        dstjet)      FILE_PREFIX="DST_JET" ;;
        dstjetcalo)  FILE_PREFIX="DST_JETCALO" ;;
        *)           fatal "BUG: unhandled DSTTYPE ‘$DSTTYPE’" ;;
      esac
      LIST_PATTERN="${FILE_PREFIX}-000*.list"
      LIST_FMT="${FILE_PREFIX}-%08d.list"
      PAD_FMT="%08d"
    fi
    ;;
  *)
    fatal "Unknown data‑set selector '$DATASET' – use run24auau or run25auau"
    ;;
esac

##############################################################################
# 1. USER CONSTANTS
##############################################################################
CHUNK_SIZE=2
MAX_JOBS=10000
##############################################################################

##############################################################################
# 2. PATH CONSTANTS
##############################################################################
PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
DST_LIST_DIR="${PROJECT_BASE}/dst_list"
TMP_LIST_DIR="${PROJECT_BASE}/tmp_condor_lists"
EXEC="${PROJECT_BASE}/run_auau_run3_qa.sh"
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"

# --- new helper: expand to <base>/<run> and mkdir if needed -------------
outdir_for_run() {
  local run="$1"
  local dir="${CONDOR_OUT_BASE}/${run}"
  mkdir -p "$dir"          # make sure it exists (harmless if it already does)
  printf '%s' "$dir"
}
LOGDIR="${PROJECT_BASE}/log"
OUTDIR="${PROJECT_BASE}/stdout"
ERRDIR="${PROJECT_BASE}/error"
mkdir -p "$TMP_LIST_DIR" "$LOGDIR" "$OUTDIR" "$ERRDIR"

RUN_SPLIT_DIR="${PROJECT_BASE}/run_segments"
SEGMENT_PREFIX="${RUN_SPLIT_DIR}/runSegment_${DATASET}_"
mkdir -p "$RUN_SPLIT_DIR"

##############################################################################
# 3. OPERATIONAL MODE PARSING
##############################################################################
mode="${1:-}"          # local | condor | condorTest | splitRunList
limitSwitch="${2:-}"   # optional 2nd token

case "$mode" in
  local|condor|condorTest|splitRunList) ;;
  '')  fatal "Missing operational mode  (local | condor | condorTest | splitRunList)" ;;
  *)   fatal "Unknown operational mode '$mode'" ;;
esac

##############################################################################
# 4. HELPER: split_run_list
##############################################################################
split_run_list() {
  local master="$1"
  [[ -f "$master" ]] || { warn "Run-list not found → $master"; return 1; }

  say  "Splitting $(basename "$master") → $RUN_SPLIT_DIR"

  local seg=1 jobs=0              # <- define first
  local current="${SEGMENT_PREFIX}${seg}.txt"
  : > "$current"

  while IFS= read -r raw; do
    [[ -z "$raw" || "$raw" =~ ^# ]] && continue
    local runNumDec=$((10#$raw))
    local listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
    [[ -f "$listFile" ]] || { warn "  – list for run $raw missing – skipped"; continue; }

    local nFiles; nFiles=$(wc -l < "$listFile")
    (( nFiles )) || { warn "  – list for run $raw empty   – skipped"; continue; }

    local nJobs=$(( (nFiles + CHUNK_SIZE - 1) / CHUNK_SIZE ))
    if (( jobs + nJobs > MAX_JOBS )); then
      good "  [SEGMENT $seg] closed with $jobs jobs"
      (( ++seg ))
      current="${SEGMENT_PREFIX}${seg}.txt"
      : > "$current"
      jobs=0
    fi

    printf "$PAD_FMT\n" "$runNumDec" >>"$current"
    (( jobs += nJobs ))
  done < "$master"

  good "  [SEGMENT $seg] closed with $jobs jobs"
  good "[OK] splitting finished – files are in $RUN_SPLIT_DIR"
}

##############################################################################
# 5. EARLY‑EXIT: splitRunList
##############################################################################
if [[ "$mode" == "splitRunList" ]]; then
  split_run_list "$limitSwitch"
  exit 0
fi

##############################################################################
# 6. VERBOSITY / CAP
##############################################################################
VERBOSE=0
[[ "$mode" == "condorTest" || ( "$mode" == "condor" && "$limitSwitch" == "firstTen" ) ]] && VERBOSE=1
vecho() { (( VERBOSE )) && echo -e "${CLR_B}•${CLR_RST} $*"; }

# ──────────────────────────────────────────────────────────────────────────
#  Clean previous output before a new Condor campaign / test run
# --------------------------------------------------------------------------
if [[ "$mode" == "condorTest" || "$mode" == "condor" ]]; then
  say  "$mode – removing previous output"
  # – wipe every run-subfolder in the bulk tree
  rm -rf "${CONDOR_OUT_BASE:?}/"* || warn "Nothing to clean in ${CONDOR_OUT_BASE}"
  # – truncate per-job stdout / log / error directories
  rm -f  "${OUTDIR:?}/"* "${LOGDIR:?}/"* "${ERRDIR:?}/"* 2>/dev/null || true
fi

# (re‑)create the job IO directories unconditionally -----------------------
mkdir -p "$LOGDIR" "$OUTDIR" "$ERRDIR"
# ──────────────────────────────────────────────────────────────────────────

jobCap=0
[[ "$mode" == "condor" && "$limitSwitch" == "firstTen" ]] && jobCap=$MAX_JOBS
submitted=0

##############################################################################
# 7. ROUND‑N OR GOLDEN‑LIST SELECTION
##############################################################################
runListFile=""

# -- explicit “round N” selection -------------------------------------------
if [[ "$mode" == "condor" && "$limitSwitch" == "round" && "${3:-}" =~ ^[0-9]+$ ]]; then
  runListFile="${SEGMENT_PREFIX}${3}.txt"
  [[ -f "$runListFile" ]] || fatal "Segment file $runListFile not found"
  say  "Round ${3} selected → using run list $(basename "$runListFile")"
fi

# -- automatic golden list ---------------------------------------------------
if [[ "$DATASET" == run25auau && -z "$runListFile" ]]; then
  if (( CALOFIT )); then
    for p in "${PROJECT_BASE}" .; do
      [[ -f "$p/run3GoldenRuns.txt" ]] && runListFile="$p/run3GoldenRuns.txt" && break
    done
    [[ -n "$runListFile" ]] && \
      say "run25auau(caloFitting) – using golden run list $(basename "$runListFile")"
  else
    for p in "${PROJECT_BASE}" .; do
      [[ -f "$p/run25GoldenRuns.txt" ]] && runListFile="$p/run25GoldenRuns.txt" && break
    done
    [[ -n "$runListFile" ]] && \
      say "run25auau – using golden run list $(basename "$runListFile")"
  fi
fi

##############################################################################
# 8. BUILD MAIN ARRAYS (runs[]  &  listFiles[])
##############################################################################
runs=()
listFiles=()

if [[ -n "$runListFile" ]]; then
  while IFS= read -r raw; do
    [[ -z "$raw" || "$raw" =~ ^# ]] && continue
    runNumDec=$((10#$raw))
    runNumPad=$(printf "$PAD_FMT" "$runNumDec")
    listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
    [[ -f "$listFile" ]] || { warn "List for run $runNumPad missing – skipped"; continue; }
    runs+=( "$runNumPad" ); listFiles+=( "$listFile" )
  done < "$runListFile"
else
  mapfile -t listFiles < <(ls "${DST_LIST_DIR}"/${LIST_PATTERN} 2>/dev/null | sort)
  (( ${#listFiles[@]} )) || fatal "No .list files found in ${DST_LIST_DIR}"
  for f in "${listFiles[@]}"; do
    bn=${f##*-}; runs+=( "${bn%.list}" )
  done
fi

##############################################################################
# 9. LOCAL MODE  (final, working version)
##############################################################################
if [[ "$mode" == "local" ]]; then
  [[ -n "${runs[0]:-}" ]] || fatal "No runs available for local mode"

  runNumber="${runs[0]}"        # default = first run in list
  maxEvt=0                      # default = analyse everything

  if [[ -n "$limitSwitch" ]]; then
    if [[ "$limitSwitch" =~ ^[0-9]+$ && "$limitSwitch" -ge 100000 ]]; then
      # token looks like a real run number
      runNumber=$(printf "$PAD_FMT" "$limitSwitch")
      maxEvt="${3:-0}"          # 2‑nd numeric token becomes event cap
    else
      # token is the event cap itself
      maxEvt="$limitSwitch"
    fi
  fi

  runNumDec=$((10#$runNumber))                       # decimal copy
  listFile="${DST_LIST_DIR}/$(printf "$LIST_FMT" "$runNumDec")"
  [[ -s "$listFile" ]] || fatal "List‑file $listFile not found or empty"

  firstDST=$(head -n1 "$listFile")
  tag=$(basename "${firstDST%.root}")                # ← new: build tag

  say  "Local test  –  run $runNumber"
  say  "First DST   : $firstDST"
  say  "maxEvt      : $maxEvt (0 → all)"

  tmpList=$(mktemp "${TMP_LIST_DIR}/local_${runNumber}_XXXX.list")
  echo "$firstDST" > "$tmpList"

  "${EXEC}"  "$runNumber"  "$tmpList"  "${tag}"  "$(outdir_for_run "$runNumber")"  "$maxEvt"
  rm -f "$tmpList"
  exit 0
fi


##############################################################################
# 10. CONDOR / CONDORTEST LOOP
##############################################################################
for idx in "${!runs[@]}"; do
  runPad=${runs[$idx]}
  runDec=$((10#$runPad))
  masterList=${listFiles[$idx]}

  vecho "Considering run $runPad  (list: $(basename "$masterList"))"

  rm -f "${TMP_LIST_DIR}/run${runPad}_chunk_"* 2>/dev/null || true
  split -l "$CHUNK_SIZE" -d -a 3 "$masterList" "${TMP_LIST_DIR}/run${runPad}_chunk_"

  mapfile -t chunks < <(ls "${TMP_LIST_DIR}/run${runPad}_chunk_"* 2>/dev/null)
  (( ${#chunks[@]} )) || { warn "Empty run $runPad – skipped."; continue; }

  if (( jobCap && submitted + ${#chunks[@]} > jobCap )); then
    good "Adding run $runPad would exceed cap ($jobCap) – stopping."
    break
  fi

  chunkNo=0
  for listFile in "${chunks[@]}"; do
    (( ++chunkNo ))
    [[ -s "$listFile" ]] || { warn "Empty chunk $listFile"; continue; }

    firstDST=$(head -n1 "$listFile")
    tag=$(basename "${firstDST%.root}")

    subFile="${TMP_LIST_DIR}/${tag}.sub"
    cat > "$subFile" <<EOS
universe      = vanilla
executable    = $EXEC
arguments     = $runPad  $listFile  ${tag}  $(outdir_for_run $runPad)
log           = ${LOGDIR}/${tag}.log
output        = ${OUTDIR}/${tag}.out
error         = ${ERRDIR}/${tag}.err
request_memory= 3000MB
+JobFlavour   = "tomorrow"
queue
EOS
    if condor_submit "$subFile" >/dev/null; then
      (( ++submitted ))
      vecho "  submitted chunk $chunkNo/${#chunks[@]}"
    else
      warn "condor_submit failed for $subFile"
    fi
  done

  vecho "Completed run $runPad – jobs now at $submitted"
  [[ "$mode" == "condorTest" ]] && { vecho "condorTest done."; break; }
done

good "Grand‑total Condor jobs submitted: $submitted"
[[ $jobCap -gt 0 ]] && say  "Launch cap in effect          : $jobCap"
##############################################################################
