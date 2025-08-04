#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh — master driver for the Run‑24/25 Au+Au QA chain
#
#  ────────────────────────────  USAGE CHEAT‑SHEET  ─────────────────────────
#
#  LOCAL WORKSTATION (macOS / Linux)         ┃  SPHENIX ANALYSIS NODES
#  ───────────────────────────────────────── ┃  ───────────────────────────────
#  ./runAuAu.sh                  # all runs  ┃  … fromSPHENIXnode condorTest
#  ./runAuAu.sh testRun          # first     ┃  … fromSPHENIXnode condor
#  ./runAuAu.sh testCombined  N  # top‑N     ┃  … fromSPHENIXnode condorOnlyRuns
#  ./runAuAu.sh runFromCurrentCombined       ┃  … fromSPHENIXnode condorOnlyCombined
#                                            ┃
#     condorTest          → 1 run (largest) + hadd/QA
#     condor              → run‑by‑run jobs **plus** hadd/QA
#     condorOnlyRuns      → run‑by‑run jobs **only**
#     condorOnlyCombined  → hadd/QA **only**
#
#  Directory roots (auto‑selected via RUN_LOCATION):
#     local   → $HOME/Desktop/auauAnalysis/emcalSEPDcorrelations
#     sphenix → /sphenix/u/$USER/scratch/emcalSEPDcorrelations
#
##############################################################################
set -euo pipefail

# ───────────────────────────────  GLOBAL CONFIG  ───────────────────────────
macro="analyzeRun24or25auau.cpp"                 # ROOT macro to run
verbose="false"      # toggle with  -v / --verbose   or  VERBOSE=1
clean_output="true"  # local run: wipe old PNG/CSV before starting

# ‑‑ pretty terminal helpers
cBlu=$'\033[1;34m'; cYel=$'\033[1;33m'; cRed=$'\033[1;31m'; cEnd=$'\033[0m'
note() { printf "${cBlu}[INFO]${cEnd}  %s\n" "$*"; }
warn() { printf "${cYel}[WARN]${cEnd}  %s\n" "$*"; }
die () { printf "${cRed}[FATAL]${cEnd} %s\n" "$*" >&2; exit 2; }

# ────────────────────────────────  CLI PARSE  ──────────────────────────────
[[ "${VERBOSE:-0}" == 1 ]] && verbose="true"
if [[ ${1:-} =~ ^--?v(erbose)?$ ]]; then verbose="true"; shift; fi

mode="${1:-}"    # major mode token
test_arg="false" sample_arg="-1"

case "${mode}" in
  "") ;;  # desktop default
  fromLocalNode)      export RUN_LOCATION=local  ; shift ;;
  fromSPHENIXnode)    export RUN_LOCATION=sphenix
                      shift; mode="${1:-condor}" ; [[ $# -gt 0 ]] && shift ;;
  testRun)            test_arg="true"            ; shift ;;
  testCombined)       [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || die "Usage: testCombined <N>"
                      sample_arg="$2"            ; shift 2 ;;
  runFromCurrentCombined)
                      export COMBINED_ONLY=1     ; clean_output="false"; shift ;;
  condorTest|condor|condorOnlyRuns|condorOnlyCombined)
                      export RUN_LOCATION=sphenix; shift ;;
  *) ;;
esac
[[ $# -ge 1 ]] && export QA_ONLY="$1"            # optional QA filter

# ────────────────────────────────  PATH  SETUP  ────────────────────────────
if [[ "${RUN_LOCATION:-local}" == "local" ]]; then
    BASE="$HOME/Desktop/auauAnalysis/emcalSEPDcorrelations"
else
    BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
fi

INPUT_DIR="${BASE}/output"               # run‑by‑run ROOTs
OUTPUT_DIR="${BASE}/outputPlots"         # PNG/CSV per run
WRAPPER="${BASE}/macros/runAuAuExecutable.sh"

LOG_DIR="${BASE}/log";  OUT_DIR="${BASE}/stdout"; ERR_DIR="${BASE}/error"
SUBMIT_DIR="${BASE}/tmp_condor_submit"

# ─────────────────────────────  UTILITY SUB‑ROUTINES  ──────────────────────
_cleanup_condor_tree() {
    note "⏳  Cleaning submission tree"
    rm -rf "${SUBMIT_DIR:?}/"* "${OUTPUT_DIR:?}" 2>/dev/null || true
    find "${LOG_DIR}" "${OUT_DIR}" "${ERR_DIR}" -type f -delete 2>/dev/null || true
    mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${OUT_DIR}" "${ERR_DIR}" "${OUTPUT_DIR}"
}

_run_files() { ls "${INPUT_DIR}"/output_*.root 2>/dev/null | sort; }   # echo list

_pick_largest() {        # $1… list of files on stdin; echo largest
    awk '{print $0}' | xargs -I{} stat -c '%s %n' {} | sort -nr | head -1 | awk '{print $2}'
}

_submit_run_job() {      # arg1=file, arg2=runId
    local f="$1" run="$2" sub="${SUBMIT_DIR}/${run}.sub"
    cat > "$sub" <<EOS
universe        = vanilla
executable      = ${WRAPPER}
arguments       = ${f}  ${OUTPUT_DIR}/${run}
output          = ${OUT_DIR}/${run}.out
error           = ${ERR_DIR}/${run}.err
log             = ${LOG_DIR}/${run}.log
request_memory  = 4GB
+JobFlavour     = "tomorrow"
queue
EOS
    condor_submit "$sub" || die "condor_submit failed for run $run"
}

_wait_logs() { condor_wait "$@" || die "condor_wait reported failure"; }

_verify_outputs() {      # arg1=expectedCount
    local need=$1 tries=0 have
    note "🔍  Verifying PNG/CSV sub‑folders (${need} expected)"
    while :; do
        have=$(find "${OUTPUT_DIR}" -maxdepth 1 -type d -name '[0-9]*' | wc -l)
        note "    progress ${have}/${need}"
        (( have >= need )) && break
        (( tries++ > 60 )) && die "Timeout waiting for run outputs"
        sleep 30
    done
}

_submit_hadd_job() {
    local sub="${SUBMIT_DIR}/haddAll.sub"
    cat > "$sub" <<'EOS'
universe        = vanilla
executable      = /bin/bash
arguments       = -c '
set -euo pipefail
export RUN_LOCATION=sphenix
BASE=/sphenix/u/${USER}/scratch/emcalSEPDcorrelations
INPUT=${BASE}/output
OUT=${BASE}/output/output_ALL_COMBINED.root
echo "[HADD] merging → \$OUT"
hadd -f -k "\$OUT" "\$INPUT"/output_*.root
echo "[HADD] running combined QA pass"
root -l -b -q -e ".L ${BASE}/macros/analyzeRun24or25auau.cpp+" \
               -e "runOneQaPass(\"\$OUT\",\"${BASE}/outputPlots/Combined\")"'
output          = ${OUT_DIR}/hadd.out
error           = ${ERR_DIR}/hadd.err
log             = ${LOG_DIR}/hadd.log
request_memory  = 4GB
+JobFlavour     = "tomorrow"
queue
EOS
    condor_submit "$sub" || die "condor_submit failed for hadd job"
}

# ───────────────────────────────  CONDOR  BRANCH  ──────────────────────────
if [[ $mode =~ ^condor ]]; then
    [[ -x "${WRAPPER}" ]] || die "Wrapper ${WRAPPER} is missing or not executable"

    mapfile -t allRuns < <(_run_files)
    [[ ${#allRuns[@]} -gt 0 ]] || die "No ROOT files found in ${INPUT_DIR}"

    case "$mode" in
        condorTest)             # keep only largest file
            allRuns=( "$(_run_files | _pick_largest)" );;
        condorOnlyCombined)     # do not touch allRuns
            ;;
    esac

    _cleanup_condor_tree

    if [[ "$mode" != "condorOnlyCombined" ]]; then
        note "🚀  Submitting run‑by‑run jobs (${#allRuns[@]} total)"
        runLogs=()
        for f in "${allRuns[@]}"; do
            run=${f##*/}; run=${run%.root}; run=${run#output_}
            note "Submitting ${run}"
            _submit_run_job "$f" "$run"
            runLogs+=( "${LOG_DIR}/${run}.log" )
        done
        _wait_logs "${runLogs[@]}"

        [[ "$mode" == "condorOnlyRuns" ]] && { note "Run‑jobs done → exiting (condorOnlyRuns)"; exit 0; }

        _verify_outputs "${#allRuns[@]}"
    fi

    note "🧩  Submitting hadd + combined‑QA job"
    _submit_hadd_job
    note "✔  Workflow finished – monitor ${LOG_DIR}/hadd.log for progress"
    exit 0
fi

# ───────────────────────  DESKTOP / INTERACTIVE RUN  ───────────────────────
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"

build_dir="$(cd "$(dirname "${macro}")" && pwd)/.aclic_build"
mkdir -p "${build_dir}"

root_cmd=( root -l -b -q -e "gSystem->SetBuildDir(\"${build_dir}\",kTRUE)" \
           "${macro}+Ok(${test_arg},${sample_arg})" )
[[ "$verbose" == true ]] && note "+ ${root_cmd[*]}"

if [[ "$clean_output" == true && -d "${OUTPUT_DIR}" ]]; then
    note "🧹  Cleaning old local outputs under ${OUTPUT_DIR}"
    rm -rf "${OUTPUT_DIR:?}/"*
fi

find "$(dirname "$macro")" -maxdepth 1 -name "$(basename "${macro%.*}")_ACLiC_*" -delete
rm -f "$(dirname "$macro")/$(basename "${macro%.*}")"{.so,.d}

"${root_cmd[@]}" 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
