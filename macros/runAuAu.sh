#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh — build & execute analyzeRun24or25auau.cpp
#
#  USAGE SUMMARY
#  =============
#      ./runAuAu.sh [--verbose]                        # analyse *all* runs
#      ./runAuAu.sh [--verbose] testRun               # only the first run
#      ./runAuAu.sh [--verbose] testCombined <N>      # N random runs + hadd
#      ./runAuAu.sh --help | -h                       # show this help
#
#  OPTIONS
#      --verbose  , -v  : show every clang++ / linker command ACLiC issues
#      --help     , -h  : print this usage summary and exit
#
#  ENVIRONMENT
#      VERBOSE=1        : same as --verbose flag
#
#  NOTES
#      • All previous call patterns still work.
#      • Duplicate‑rpath warnings from Apple ld are filtered out automatically.
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
  testRun)       test_arg="true" ;;
  testCombined)
        if [[ $# -lt 2 || ! "$2" =~ ^[0-9]+$ ]]; then
            echo "Usage: ./runAuAu.sh [--verbose] testCombined <N>" >&2
            exit 1
        fi
        sample_arg="$2"
        ;;
  *)  echo "Unknown option: ${mode}" >&2
      echo "Valid options:  (none) | testRun | testCombined <N>" >&2
      exit 1
      ;;
esac

# ────────────────────────────────────────────────────────────────────────────
# 2.  Silence duplicate‑rpath warnings from Apple ld
# ────────────────────────────────────────────────────────────────────────────
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"

# ────────────────────────────────────────────────────────────────────────────
# 3.  Build the ROOT command line
#     • use ‑q only in non‑verbose mode
# ────────────────────────────────────────────────────────────────────────────
root_flags=(-l -b)                # ‑l = no splash screen
if [[ "${verbose}" == "false" ]]; then
    root_flags+=(-q)              # hide extra ROOT chatter
fi

root_cmd=(root "${root_flags[@]}" \
          "${macro}+Ok(${test_arg},${sample_arg})")

# Optional: show the exact command we are about to run when verbose
[[ "${verbose}" == "true" ]] && echo "+ ${root_cmd[*]}" >&2

# ────────────────────────────────────────────────────────────────────────────
# 3.5  Clean previous output so every run starts with an empty directory
# ────────────────────────────────────────────────────────────────────────────
output_root="${HOME}/Desktop/auauAnalysis/emcalSEPDcorrelations/output"

if [[ -d "${output_root}" ]]; then
    echo "Cleaning old output under ${output_root}" >&2
    rm -rf "${output_root:?}/"*
fi

# ────────────────────────────────────────────────────────────────────────────
# 4.  Execute — filter only the duplicate-rpath warning
# ────────────────────────────────────────────────────────────────────────────
{
    "${root_cmd[@]}"
} 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
