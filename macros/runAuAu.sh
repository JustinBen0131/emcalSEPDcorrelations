#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh  –  single‑entry driver for the Run‑24/25 Au+Au QA workflow
#
#  QUICK REFERENCE ──────────────────────────────────────────────────────────
#  Local machine
#    ./runAuAu.sh                           : analyse every run found locally
#    ./runAuAu.sh testRun                   : analyse the very first run only
#    ./runAuAu.sh testCombined <N>          : top‑N runs → hadd → full QA pass
#    ./runAuAu.sh runFromCurrentCombined    : QA pass on an existing COMBINED file
#
#  sPHENIX analysis node
#    ./runAuAu.sh fromSPHENIXnode condorTest               : smoke‑test (1 HTCondor job)
#    ./runAuAu.sh fromSPHENIXnode condor                   : full grid (1 job / run)
#    ./runAuAu.sh fromSPHENIXnode haddAndFinalize local    : **local** merge + combined QA
#      ├─ optional  haddAndFinalize        → submit only the merge + combined QA job (Condor)
#      ├─ optional  wipePrevPlots          → delete **all** PNGs under $OUTPUT_DIR before run
#      └─ examples
#             ./runAuAu.sh fromSPHENIXnode condor haddAndFinalize        # remote final pass
#             ./runAuAu.sh fromSPHENIXnode condorTest wipePrevPlots      # smoke‑test, clean PNGs
#
#  FLAGS / OPTIONS ──────────────────────────────────────────────────────────
#    haddAndFinalize          remote merge + combined QA (Condor)
#    local                    when placed **after haddAndFinalize**, do that step locally
#    wipePrevPlots            global wipe of $OUTPUT_DIR before submission
#    VERBOSE=1|2              extra shell diagnostics   (or add ‑v for level 1)
#
#  CLEAN‑UP MATRIX ──────────────────────────────────────────────────────────
#                        | tmp submit dir | log/out/err | outputPlots PNG tree
#    --------------------+---------------+-------------+----------------------
#    condor / condorTest | wiped always  | run‑IDs only| untouched
#    +wipePrevPlots      | ″             | ″           | **fully wiped**
#    haddAndFinalize     | ″             | generic     | untouched
#    haddAndFinalize local| N/A          | N/A         | untouched
#
#  PATHS  (auto‑detected) ───────────────────────────────────────────────────
#    RUN_LOCATION=local     → $HOME/Desktop/auauAnalysis/…
#    RUN_LOCATION=sphenix   → /sphenix/u/$USER/scratch/emcalSEPDcorrelations
#
#  HTCondor artefacts live in   $PROJECT_BASE/{log | stdout | error}
##############################################################################

set -euo pipefail

# ────────────────────────────────────────────────────────────────────────────
# 0.  Global flags
# ────────────────────────────────────────────────────────────────────────────
macro="analyzeRun24or25auau.cpp"   # C++ macro to build/run
verbose="false"
clean_output="true"
finalOnly="false"
localFinalize="false"          # ← run merge + QA locally
wipePlots="false"              # ← wipe outputPlots only if user adds ‘wipePrevPlots’

# ---------- numeric verbosity (environment or first positional) ----------
if [[ ${1:-} == VERBOSE=* ]]; then
    export VERBOSE="${1#VERBOSE=}"     # take the number after '='
    shift                              # drop this pseudo‑argument
fi
: "${VERBOSE:=0}"                      # default if nothing supplied

# ---------- script‑level noisy/quiet switch ------------------------------
if (( VERBOSE >= 1 )); then
    verbose="true"                     # enables extra shell messages
fi

# legacy --verbose / -v still sets VERBOSE=1
if [[ ${1:-} == "--verbose" || ${1:-} == "-v" ]]; then
    verbose="true"; export VERBOSE=1; shift
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
          export RUN_LOCATION=sphenix            # use scratch tree
          shift                                  # drop the keyword itself

          # ――― default submission mode ────────────────────────────────
          mode="condor"

          # ――― process all following keywords in any order ―───────────
          while [[ $# -gt 0 ]]; do
                case "$1" in
                    condor|condorTest)
                        mode="$1"
                        ;;
                    haddAndFinalize)
                        finalOnly="true"           # Condor final‑pass
                        ;;
                    local)
                        localFinalize="true"       # *local* final‑pass
                        ;;
                    wipePrevPlots)
                        wipePlots="true"
                        ;;
                    --)  shift; break ;;           # end‑of‑options marker
                    *)   break ;;                  # first non‑option → stop parsing
                esac
                shift
          done

          # if the user asked for the local variant, force a dedicated mode
          if [[ "${localFinalize}" == "true" ]]; then
                mode="localFinalize"
          fi
          ;;
  testRun)            test_arg="true"            ; shift ;;
  testCombined)       [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] \
                         || { echo "Usage: … testCombined <N>"; exit 1; }
                      sample_arg="$2" ; shift 2 ;;
  runFromCurrentCombined)
                      export COMBINED_ONLY=1; clean_output="false"; shift ;;
  condorTest)         export RUN_LOCATION=sphenix; mode="condorTest"; shift ;;
  condor)             export RUN_LOCATION=sphenix; mode="condor"    ; shift ;;
  haddAndFinalizeLocal)
                      export RUN_LOCATION=sphenix; localFinalize="true"; shift ;;
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

    # ------------------------------------------------------------------
    #  haddAndFinalize  →  exactly one job that merges all runs and
    #                      performs the combined QA pass
    # ------------------------------------------------------------------
   if [[ "${finalOnly}" == "true" ]]; then
        export HADD_AND_FINALIZE=1
        note "haddAndFinalize mode – submitting single finalize job"

        PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
        EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh"   # ← reuse one wrapper

        SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"
        LOG_DIR="${PROJECT_BASE}/log"
        STDOUT_DIR="${PROJECT_BASE}/stdout"
        STDERR_DIR="${PROJECT_BASE}/error"

        [[ -x "${EXEC_WRAPPER}" ]] || die "Wrapper ${EXEC_WRAPPER} missing or not executable"
        rm -rf "${SUBMIT_DIR:?}/"* 2>/dev/null || true
        for f in "${LOG_DIR}"/haddAndFinalize.* \
                 "${STDOUT_DIR}"/haddAndFinalize.* \
                 "${STDERR_DIR}"/haddAndFinalize.*; do
            [[ -e "$f" ]] && rm -f "$f"
        done
        mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}"

        # ⚠  deliberate: we do NOT touch \$OUTPUT_DIR here → existing plots stay intact

cat > "${SUBMIT_DIR}/haddAndFinalize.sub" <<EOS
        universe        = vanilla
        executable      = ${EXEC_WRAPPER}
        arguments       = --final
        output          = ${STDOUT_DIR}/haddAndFinalize.out
        error           = ${STDERR_DIR}/haddAndFinalize.err
        log             = ${LOG_DIR}/haddAndFinalize.log
        getenv          = True
        request_memory  = 8GB
        +JobFlavour     = "tomorrow"
        queue
EOS
        condor_submit "${SUBMIT_DIR}/haddAndFinalize.sub" \
            || die "condor_submit failed for finalize job"
        note "Finalize job submitted – exiting."
        exit 0
    fi

    note "Starting Condor submission (${mode})"

    PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
    EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh"

    SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"
    LOG_DIR="${PROJECT_BASE}/log"
    STDOUT_DIR="${PROJECT_BASE}/stdout"
    STDERR_DIR="${PROJECT_BASE}/error"

    INPUT_DIR="${PROJECT_BASE}/output"
    OUTPUT_DIR="/sphenix/tg/tg01/bulk/jbennett/GLOBAL_QA/outputPlots"

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

    # always clear temporary submit files
    rm -rf "${SUBMIT_DIR:?}/"* 2>/dev/null || true

    # wipe the entire PNG tree only when the user asked for it
    if [[ "${wipePlots}" == "true" && "${finalOnly}" != "true" ]]; then
        rm -rf "${OUTPUT_DIR:?}/"* 2>/dev/null || true
    fi

    # ── remove log/stdout/stderr files from any earlier attempt ────────────
    for rf in "${roots[@]}"; do
        run="$(basename "${rf}" .root)"; run="${run#output_}"
        rm -f  "${LOG_DIR}/${run}.log"    2>/dev/null || true
        rm -f  "${STDOUT_DIR}/${run}.out" 2>/dev/null || true
        rm -f  "${STDERR_DIR}/${run}.err" 2>/dev/null || true
    done

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
        getenv          = True             # ← pass entire env, incl. VERBOSE
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

# ──────────────────────────────────────────────────────────────────
# 3‑B.  haddAndFinalizeLocal  –  run merge + combined QA *locally*
# ──────────────────────────────────────────────────────────────────
if [[ "${localFinalize}" == "true" ]]; then
    note "haddAndFinalizeLocal – running merge + combined QA locally"
    PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
    EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh"
    [[ -x "${EXEC_WRAPPER}" ]] || die "Wrapper ${EXEC_WRAPPER} missing or not executable"
    export HADD_AND_FINALIZE=1
    "${EXEC_WRAPPER}" --final
    note "Local finalize completed."
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
