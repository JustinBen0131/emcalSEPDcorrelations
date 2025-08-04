#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh  –  single‑entry driver for the Run‑24/25 Au+Au QA macro
#
#  SUPPORTED INVOCATIONS
#  ───────────────────────────────────────────────────────────────────────────
#  DESKTOP  (macOS / personal laptop)
#  ───────────────────────────────────────────────────────────────────────────
#    ./runAuAu.sh [qa_list]                     # full pass over all runs
#    ./runAuAu.sh testRun [qa_list]             # first run only
#    ./runAuAu.sh testCombined <N> [qa_list]    # top‑N runs, then hadd+QA
#    ./runAuAu.sh runFromCurrentCombined        # *just* reruns QA on the
#                                               #   already‑existing hadd file
#
#  SPHENIX ANALYSIS NODES (scratch area ↔ Condor)
#  ───────────────────────────────────────────────────────────────────────────
#    ./runAuAu.sh fromSPHENIXnode condorTest    # 1 Condor job (smoke test)
#    ./runAuAu.sh fromSPHENIXnode condor        # one Condor job per run,
#                                               # then automatic hadd+QA
#
#  PATH POLICY
#  ───────────────────────────────────────────────────────────────────────────
#    • Environment variable RUN_LOCATION is set automatically:
#        local    → $HOME/Desktop/… tree
#        sphenix  → /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations
#    • Those paths are *immutable* – edit only in the C++ macro if required.
#
#  LOG / STDOUT / STDERR DESTINATION (Condor mode)
#    /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/log
#    /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/stdout
#    /sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/error
#
#  All other behaviour (QA filtering via QA_ONLY, verbosity via ‑v or
#  VERBOSE=1, cleaning rules, etc.) remains unchanged.
##############################################################################
set -euo pipefail

macro="analyzeRun24or25auau.cpp"         # C++ macro to build/run

# ────────────────────────────────────────────────────────────────────────────
# 0.  Global verbosity flag (CLI --verbose | -v  or  VERBOSE=1 env var)
# ────────────────────────────────────────────────────────────────────────────
verbose="false"
clean_output="true"
if [[ "${VERBOSE:-0}" == 1 ]]; then verbose="true"; fi
if [[ ${1:-} == "--verbose" || ${1:-} == "-v" ]]; then
    verbose="true"; shift         # drop the flag from $@
fi

# ────────────────────────────────────────────────────────────────────────────
# 1.  MODE PARSER  (empty | testRun | testCombined <N>)
# ────────────────────────────────────────────────────────────────────────────
mode="${1:-}"
test_arg="false"     # C++ parm #1
sample_arg="-1"      # C++ parm #2

case "${mode}" in
  "") ;;                                        # desktop – full suite
  fromLocalNode)                               # explicit alias; keeps default paths
        export RUN_LOCATION=local
        shift 1
        ;;
  fromSPHENIXnode)                             # ⇢ Condor launcher (all runs)
        export RUN_LOCATION=sphenix
        mode="condor"                          # remember for later logic
        shift 1
        ;;
  testRun)
        test_arg="true"
        shift 1
        ;;
  testCombined)
        [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Usage: … testCombined <N>"; exit 1; }
        sample_arg="$2"
        shift 2
        ;;
  runFromCurrentCombined)
        export COMBINED_ONLY=1
        clean_output="false"
        shift 1
        ;;
  condorTest)                                  # submit *one* Condor job
        export RUN_LOCATION=sphenix
        mode="condorTest"
        shift 1
        ;;
  condor)                                      # submit all Condor jobs
        export RUN_LOCATION=sphenix
        mode="condor"
        shift 1
        ;;
  *) ;;
esac


# ──────────────────────────────────────────────────────────────────────────
#  Condor (and condorTest) submission
# ──────────────────────────────────────────────────────────────────────────
if [[ "${mode}" == "condor" || "${mode}" == "condorTest" ]]; then
    PROJECT_BASE="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
    EXEC_WRAPPER="${PROJECT_BASE}/runAuAuExecutable.sh"

    SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"

    LOG_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/log"
    STDOUT_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/stdout"
    STDERR_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/error"

    INPUT_DIR="${PROJECT_BASE}/output"
    OUTPUT_DIR="${PROJECT_BASE}/outputPlots"

    # ── purge artefacts from any previous Condor campaign ──────────────────────
    rm -rf  "${SUBMIT_DIR:?}/"*            2>/dev/null || true          # stale .sub / .dag
    find    "${LOG_DIR}"    -type f -delete 2>/dev/null || true         # old .log files
    find    "${STDOUT_DIR}" -type f -delete 2>/dev/null || true         # old .out files
    find    "${STDERR_DIR}" -type f -delete 2>/dev/null || true         # old .err files
    rm -rf  "${OUTPUT_DIR:?}"              2>/dev/null || true          # previous PNG/CSV tree

    # ── recreate the cleaned directories ───────────────────────────────────────
    mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}" "${OUTPUT_DIR}"


    mapfile -t roots < <(ls "${INPUT_DIR}"/output_*.root 2>/dev/null | sort)
    [[ ${#roots[@]} -gt 0 ]] || { echo "[FATAL] no ROOT files in ${INPUT_DIR}"; exit 2; }

    [[ "${mode}" == "condorTest" ]] && roots=( "${roots[0]}" )   # 1st run only

    dag="${SUBMIT_DIR}/runAuAu.dag"; >"${dag}"
    jobIds=()

    for rf in "${roots[@]}"; do
        bn=$(basename "${rf}")
        run=${bn#output_}; run=${run%.root}
        sub="${SUBMIT_DIR}/${run}.sub"

        cat >"${sub}" <<EOS
universe      = vanilla
executable    = ${EXEC_WRAPPER}
arguments     = ${rf}  ${OUTPUT_DIR}/${run}
output        = ${STDOUT_DIR}/${run}.out
error         = ${STDERR_DIR}/${run}.err
log           = ${LOG_DIR}/${run}.log
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOS
        echo "JOB  J${run}  ${sub}" >>"${dag}"
        jobIds+=( "J${run}" )
    done

    # ----- post‑processing (hadd + combined QA) ---------------------------
    haddSub="${SUBMIT_DIR}/haddAll.sub"
    cat >"${haddSub}" <<'EOS'
universe      = vanilla
executable    = /bin/bash
arguments     = -c '
set -euo pipefail
export RUN_LOCATION=sphenix
BASE=/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations
INPUT=${BASE}/output
OUT=${BASE}/output/output_ALL_COMBINED.root
hadd -f -k "${OUT}" "${INPUT}"/output_*.root
root -l -b -q -e ".L ${BASE}/analyzeRun24or25auau.cpp+" \
                -e "runOneQaPass(\"${OUT}\",\"${BASE}/outputPlots/Combined\")"'
output        = ${STDOUT_DIR}/hadd.out
error         = ${STDERR_DIR}/hadd.err
log           = ${LOG_DIR}/hadd.log
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOS

    echo "JOB  HADD  ${haddSub}"  >>"${dag}"
    printf 'PARENT %s CHILD HADD\n' "${jobIds[@]}" >>"${dag}"

    condor_submit_dag -update_submit 0 "${dag}"
    echo "[OK] Condor DAG submitted."
    exit 0
fi


# ------------------------------------------------------------------
# Any remaining positional argument is treated as the QA module list
# ------------------------------------------------------------------
if [[ $# -ge 1 ]]; then
    export QA_ONLY="$1"          # e.g.  correlations,hcal
    shift
fi


# ────────────────────────────────────────────────────────────────────────────
# 2.  Silence duplicate‑rpath warnings from Apple ld
# ────────────────────────────────────────────────────────────────────────────
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"

# ────────────────────────────────────────────────────────────────────────────
# 3.  Build the ROOT command line
#     • use ‑q only in non‑verbose mode
# ────────────────────────────────────────────────────────────────────────────
root_flags=(-l -b -q)

# Define a dedicated build directory _once_ per run
build_dir="$(cd "$(dirname "${macro}")" && pwd)/.aclic_build"
mkdir -p "${build_dir}"

root_cmd=(
  root "${root_flags[@]}"
  -e "gSystem->SetBuildDir(\"${build_dir}\", kTRUE)"          # ← key line
  "${macro}+Ok(${test_arg},${sample_arg})"
)


# Optional: show the exact command we are about to run when verbose
[[ "${verbose}" == "true" ]] && echo "+ ${root_cmd[*]}" >&2

# ────────────────────────────────────────────────────────────────────────────
# 3.5  Clean previous output *and* stale ACLiC artefacts
# ────────────────────────────────────────────────────────────────────────────
output_root="${HOME}/Desktop/auauAnalysis/emcalSEPDcorrelations/output"
[[ "${RUN_LOCATION:-local}" == "local" ]] || clean_output="false"

# ❶ purge old QA PNG / CSV output  (unless combined‑only run)
if [[ "${clean_output}" == "true" && -d "${output_root}" ]]; then
    echo "Cleaning old output under ${output_root}" >&2
    rm -rf "${output_root:?}/"*
fi

# ❷ purge every ACLiC file that may linger from earlier builds
macro_dir="$(cd "$(dirname "${macro}")" && pwd)"
macro_base="$(basename "${macro%.*}")_cpp"

echo "Removing stale ACLiC artefacts in ${macro_dir}" >&2
find "${macro_dir}" -maxdepth 1 -type f -name "${macro_base}_ACLiC_*" -delete
rm -f "${macro_dir}/${macro_base}.so" "${macro_dir}/${macro_base}.d"

# ────────────────────────────────────────────────────────────────────────────
# 4.  Execute — filter only the duplicate-rpath warning
# ────────────────────────────────────────────────────────────────────────────
{
    "${root_cmd[@]}"
} 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
