#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh — compile & execute analyzeRun24or25auau.cpp with clean output #
##############################################################################
set -euo pipefail                    # robust bash: fail on any error/undef

# ---------- 1. choose macro argument ----------------------------------------
macro="analyzeRun24or25auau.cpp"
test_arg="false"
if [[ ${1:-} == "testRun" ]]; then
  test_arg="true"    # run only the first available ROOT file
fi

# ---------- 2. suppress duplicate‑library chatter (Apple ld ≥ Xcode 15) -----
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"

# ---------- 3. build the ROOT command ---------------------------------------
root_cmd=(root -b -q "${macro}+Ok(${test_arg})")

# ---------- 4. run — strip *only* the duplicate‑rpath line ------------------
# keep exit status; show every other message unchanged
{
  "${root_cmd[@]}"
} 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
