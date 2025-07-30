#!/usr/bin/env bash
#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh — build & execute analyzeRun24or25auau.cpp
#
#  QUICK START
#  -----------
#      ./runAuAu.sh                       # all runs, all QA modules
#      ./runAuAu.sh correlations,jetqa    # all runs, CorrQA + JetQA only
#      ./runAuAu.sh testRun hcal,pi0      # first run, HCal + Pi0 QA
#      ./runAuAu.sh testCombined 10       # top‑10 runs merged, full QA suite
#      ./runAuAu.sh testCombined 10 correlations,hcal
#
#  FULL USAGE
#  ==========
#      ./runAuAu.sh [--verbose] [qa_list]
#      ./runAuAu.sh [--verbose] testRun [qa_list]
#      ./runAuAu.sh [--verbose] testCombined <N> [qa_list]
#      ./runAuAu.sh --help | -h
#
#  POSITIONAL ARGUMENTS
#      testRun                 analyse only the first discovered run
#      testCombined <N>        merge the N highest‑statistics runs, then analyse
#      qa_list                 (optional) comma‑separated list of QA‑module tags
#                              to execute, e.g.  correlations,hcal,jetqa
#                              — omit to run the complete QA suite.
#
#  OPTIONS
#      --verbose , -v          show every clang++ / linker command printed by ACLiC
#      --help    , -h          print this summary and exit
#
#  ENVIRONMENT
#      VERBOSE=1               same effect as --verbose
#      QA_ONLY                 alternative place to provide <qa_list>;
#                              the command‑line argument, if present, overrides it
#
#  NOTES
#      • All previous call patterns still work unchanged.
#      • Output PNG/CSV files under  $HOME/Desktop/auauAnalysis/…  are purged
#        automatically before each run.
#      • Duplicate‑rpath warnings from Apple ld are filtered out of the log.
##############################################################################
set -euo pipefail

macro="analyzeRun24or25auau.cpp"         # C++ macro to build/run

# ────────────────────────────────────────────────────────────────────────────
# 0.  Global verbosity flag (CLI --verbose | -v  or  VERBOSE=1 env var)
# ────────────────────────────────────────────────────────────────────────────
verbose="false"
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
  "") ;;
  testRun)
        test_arg="true"
        shift 1
        ;;
  testCombined)
        if [[ $# -lt 2 || ! "$2" =~ ^[0-9]+$ ]]; then
            echo "Usage: ./runAuAu.sh [--verbose] testCombined <N> [qa_list]" >&2
            exit 1
        fi
        sample_arg="$2"
        shift 2
        ;;
  *) ;;
esac

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

# ❶ purge old QA PNG / CSV output
if [[ -d "${output_root}" ]]; then
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
