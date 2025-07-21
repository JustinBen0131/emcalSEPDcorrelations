#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh — compile & execute analyzeRun24or25auau.cpp with clean output #
##############################################################################
set -euo pipefail                       # robust bash: fail on any error/undef

macro="analyzeRun24or25auau.cpp"        # C++ macro to build/run

# ────────────────────────────────────────────────────────────────────────────
# 1. ARGUMENT PARSER  (empty | testRun | testCombined <N>)
# ────────────────────────────────────────────────────────────────────────────
mode="${1:-}"          # first CLI token, if any
test_arg="false"       # C++ parm #1  (true → keep only first file)
sample_arg="-1"        # C++ parm #2  (‑1  → no random sampling)

case "${mode}" in
  "") ;;                                     # default: analyse all runs
  testRun)
      test_arg="true"                        # QA just the first input file
      ;;
  testCombined)
      if [[ $# -lt 2 || ! "$2" =~ ^[0-9]+$ ]]; then
          echo "Usage: ./runAuAu.sh testCombined <N>" >&2
          exit 1
      fi
      sample_arg="$2"                        # QA N random runs, then hadd
      ;;
  *)
      echo "Unknown option: ${mode}" >&2
      echo "Valid options:  (none) | testRun | testCombined <N>" >&2
      exit 1
      ;;
esac

# ────────────────────────────────────────────────────────────────────────────
# 2. Suppress duplicate‑library chatter (Apple ld ≥ Xcode 15)
# ────────────────────────────────────────────────────────────────────────────
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"

# ────────────────────────────────────────────────────────────────────────────
# 3. Build the ROOT command
#    (passes the two C++ bool/int parameters we just parsed)
# ────────────────────────────────────────────────────────────────────────────
root_cmd=(root -b -q "${macro}+Ok(${test_arg},${sample_arg})")

# ────────────────────────────────────────────────────────────────────────────
# 4. Run — strip *only* the duplicate‑rpath warning line; preserve exit code
# ────────────────────────────────────────────────────────────────────────────
{
  "${root_cmd[@]}"
} 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
