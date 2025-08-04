#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh  –  single‑entry driver for the Run‑24/25 Au+Au QA macro
#
#  ── Local laptop  ──────────────────────────────────────────────────────────
#     ./runAuAu.sh [qa_list]                  # analyse all runs
#     ./runAuAu.sh testRun [qa_list]          # first run only
#     ./runAuAu.sh testCombined <N> [qa_list] # top‑N runs, then hadd+QA
#     ./runAuAu.sh runFromCurrentCombined     # QA on existing hadd file
#
#  ── sPHENIX analysis nodes  ────────────────────────────────────────────────
#     ./runAuAu.sh fromSPHENIXnode condorTest # submit one job (smoke test)
#     ./runAuAu.sh fromSPHENIXnode condor     # full DAG – one job per run AND a final hadd/process job on all runs
#
#  ENVIRONMENT
#     RUN_LOCATION is set automatically:
#         local    → $HOME/Desktop/… tree
#         sphenix  → /sphenix/u/<user>/scratch/emcalSEPDcorrelations
#
#  LOG / STDOUT / STDERR (Condor)
#     /scratch/emcalSEPDcorrelations/{log|stdout|error}
##############################################################################
set -euo pipefail

# ────────────────────────────────────────────────────────────────────────────
# 0.  Global flags
# ────────────────────────────────────────────────────────────────────────────
macro="analyzeRun24or25auau.cpp"   # C++ macro to build/run
verbose="false"
clean_output="true"

if [[ "${VERBOSE:-0}" == 1 ]]; then verbose="true"; fi
if [[ ${1:-} == "--verbose" || ${1:-} == "-v" ]]; then
    verbose="true"; shift
fi

# ────────────────────────────────────────────────────────────────────────────
# 1.  MODE PARSER
# ────────────────────────────────────────────────────────────────────────────
mode="${1:-}"
test_arg="false"        # C++ parameter #1
sample_arg="-1"         # C++ parameter #2

case "${mode}" in
  "") ;;                                           # desktop – full suite
  fromLocalNode)      export RUN_LOCATION=local  ; shift ;;
  fromSPHENIXnode)
          export RUN_LOCATION=sphenix          # running on a sPHENIX node
          shift                                # drop the keyword itself
          mode="${1:-condor}"                  # take next token ⇢ condor / condorTest
          [[ $# -gt 0 ]] && shift              # consume it when present
          ;;
  testRun)            test_arg="true"            ; shift ;;
  testCombined)       [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] \
                         || { echo "Usage: … testCombined <N>"; exit 1; }
                      sample_arg="$2" ; shift 2 ;;
  runFromCurrentCombined)
                      export COMBINED_ONLY=1; clean_output="false"; shift ;;
  condorTest)         export RUN_LOCATION=sphenix; mode="condorTest"; shift ;;
  condor)             export RUN_LOCATION=sphenix; mode="condor"    ; shift ;;
  *) ;;
esac

# ────────────────────────────────────────────────────────────────────────────
# 2.  Helper logger functions
# ────────────────────────────────────────────────────────────────────────────
clr_blu=$'\033[1;34m'; clr_red=$'\033[1;31m'; clr_end=$'\033[0m'
note() { printf "${clr_blu}[INFO]${clr_end}  %s\n" "$*"; }
warn() { printf "${clr_red}[WARN]${clr_end}  %s\n" "$*"; }
die () { printf "${clr_red}[FATAL]${clr_end} %s\n" "$*" >&2; exit 2; }

# ────────────────────────────────────────────────────────────────────────────
# 3.  Condor submission pathway
# ────────────────────────────────────────────────────────────────────────────
if [[ "${mode}" == "condor" || "${mode}" == "condorTest" ]]; then
    note "Starting Condor submission (${mode})"

    PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
    EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh"

    SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"
    LOG_DIR="${PROJECT_BASE}/log"
    STDOUT_DIR="${PROJECT_BASE}/stdout"
    STDERR_DIR="${PROJECT_BASE}/error"

    INPUT_DIR="${PROJECT_BASE}/output"
    OUTPUT_DIR="${PROJECT_BASE}/outputPlots"

    [[ -x "${EXEC_WRAPPER}" ]] || die "Wrapper ${EXEC_WRAPPER} missing or not executable"

    mapfile -t roots < <(ls "${INPUT_DIR}"/output_*.root 2>/dev/null | sort)
    [[ ${#roots[@]} -gt 0 ]] || die "No ROOT files in ${INPUT_DIR}"

    # pick the ROOT file that has the *most* events (proxy = size)
    if [[ "${mode}" == "condorTest" ]]; then
        note "condorTest → selecting run with the highest statistics"

        largest=""
        maxSize=0
        for f in "${roots[@]}"; do
            # GNU/Linux uses “stat -c%s”, macOS/BSD uses “stat -f%z”
            sz=$( { stat -c%s "$f" 2>/dev/null || stat -f%z "$f"; } ) || sz=0
            (( sz > maxSize )) && { maxSize=$sz; largest="$f"; }
        done

        [[ -n "${largest}" ]] || die "Could not determine the largest ROOT file"
        roots=( "${largest}" )

        runPick=$(basename "${largest}")
        note "condorTest → will submit only ${runPick}  (size $((maxSize/1024/1024)) MB)"
    fi

    note "Cleaning previous submission artefacts"
    rm -rf "${SUBMIT_DIR:?}/"* \
           "${OUTPUT_DIR:?}" \
           "${LOG_DIR:?}/"* "${STDOUT_DIR:?}/"* "${STDERR_DIR:?}/"* 2>/dev/null || true              2>/dev/null || true
    mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}" "${OUTPUT_DIR}"

    # ──────────────────────────────────────────────────────────────────
    #  submit one job per run; wait; verify output; run hadd+QA
    # ──────────────────────────────────────────────────────────────────
    runLogs=()                          # HTCondor .log files for condor_wait
    expRuns=${#roots[@]}                # how many runs we expect to finish

    for rf in "${roots[@]}"; do
        run="$(basename "${rf}" .root)"; run="${run#output_}"
        sub="${SUBMIT_DIR}/${run}.sub"

        note "Submitting run ${run}"
        note "   input    : ${rf}"
        note "   out‑dir  : ${OUTPUT_DIR}/${run}"

        cat > "${sub}" <<EOS
        universe        = vanilla
        executable      = ${EXEC_WRAPPER}
        arguments       = ${rf}  ${OUTPUT_DIR}/${run}
        output          = ${STDOUT_DIR}/${run}.out
        error           = ${STDERR_DIR}/${run}.err
        log             = ${LOG_DIR}/${run}.log
        request_memory  = 4GB
        +JobFlavour     = "tomorrow"
        queue
EOS
        
        condor_submit "${sub}"              || die "condor_submit failed for ${run}"
        runLogs+=( "${LOG_DIR}/${run}.log" )
    done
    note "All ${expRuns} run-jobs submitted – exiting."
    exit 0
fi

# --------------------------------------------------------------------------
# 4.  Optional QA‑module filter on the command line
# --------------------------------------------------------------------------
if [[ $# -ge 1 ]]; then
    export QA_ONLY="$1"; shift
fi

# --------------------------------------------------------------------------
# 5.  Desktop / interactive execution (unchanged, just more noise if -v)
# --------------------------------------------------------------------------
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"
root_flags=(-l -b -q)

build_dir="$(mktemp -d "${TMPDIR:-/tmp}/aclic_build_XXXXXX")"
mkdir -p "${build_dir}"

root_cmd=(
  root "${root_flags[@]}"
  -e "gSystem->SetBuildDir(\"${build_dir}\",kTRUE)"
  "${macro}+Ok(${test_arg},${sample_arg})"
)
[[ "${verbose}" == "true" ]] && note "+ ${root_cmd[*]}"

output_root="${HOME}/Desktop/auauAnalysis/emcalSEPDcorrelations/output"
[[ "${RUN_LOCATION:-local}" == "local" ]] || clean_output="false"

if [[ "${clean_output}" == "true" && -d "${output_root}" ]]; then
    note "Cleaning old local output under ${output_root}"
    rm -rf "${output_root:?}/"*
fi

macro_dir="$(cd "$(dirname "${macro}")" && pwd)"
macro_base="$(basename "${macro%.*}")_cpp"

# ── wipe the previous ACLiC build area and recreate it ────────────────
rm -rf "${macro_dir}/.aclic_build"
mkdir -p "${macro_dir}/.aclic_build"

# ── remove any stray libraries from earlier compilations ──────────────
find "${macro_dir}" -maxdepth 1 -type f -name "${macro_base}_ACLiC_*" -delete
rm -f "${macro_dir}/${macro_base}".{so,d}

"${root_cmd[@]}" 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
