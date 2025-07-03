#!/usr/bin/env bash
###############################################################################
# forceFileListSimCreation.sh – build .list files for “good” run numbers
#   • Reads  ../goodRunList_sEPD_run24auau.txt
#   • Writes ../dst_list/DST_CALO_run2auau_new_2024p007-<run>.list
#   • Searches   run_<low>_<high>/  (100-run bins) under BASE_DIR
#       ↳ if empty, retries one level higher
#   • If still no ROOT files, records the run and deletes any stale .list later
###############################################################################
set -euo pipefail
IFS=$'\n\t'

########################  COLOUR / LOG HELPERS  ###############################
ESC=$'\e['
CLR_RED=${ESC}0\;31m;  CLR_GRN=${ESC}0\;32m
CLR_YEL=${ESC}1\;33m;  CLR_BLU=${ESC}1\;34m
CLR_BLD=${ESC}1m;      CLR_RST=${ESC}0m

say()   { printf "${CLR_BLU}➜${CLR_RST} %s\n" "$*"; }
good()  { printf "${CLR_GRN}%s${CLR_RST}\n"   "$*"; }
warn()  { printf "${CLR_YEL}⚠ %s${CLR_RST}\n" "$*" >&2; }
fatal() { printf "${CLR_RED}✘ %s${CLR_RST}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO)"' ERR

############################  CONSTANTS  #####################################
SCRIPT_PWD="$(pwd)"
RUNLIST_FILE="${SCRIPT_PWD}/../goodRunList_sEPD_run24auau.txt"
LIST_DIR="${SCRIPT_PWD}/../dst_list"
BASE_DIR="/sphenix/lustre01/sphnxpro/physics/run2auau/caloy2calib/new_2024p007"
FILE_PREFIX="DST_CALO_run2auau_new_2024p007"

[[ -f "$RUNLIST_FILE" ]] || fatal "Run-list file not found: $RUNLIST_FILE"
mkdir -p "$LIST_DIR"

############################  PRE-FLIGHT  #####################################
say  "Working directory   : ${CLR_BLD}${SCRIPT_PWD}${CLR_RST}"
say  "Run-list file       : ${CLR_BLD}${RUNLIST_FILE}${CLR_RST}"
say  "Destination .list → : ${CLR_BLD}${LIST_DIR}${CLR_RST}"
say  "ROOT files base Dir : ${CLR_BLD}${BASE_DIR}${CLR_RST}"
echo

total_runs=$(grep -c '^[0-9]\+' "$RUNLIST_FILE")
good "Runs to process     : ${total_runs}"
echo

###########################  COUNTERS  ########################################
ok=0  miss_dir=0  miss_files=0  removed=0
inc() { : $(( $1 += 1 )); }          # safe under 'set -e'
declare -a cleanup_runs=()           # runs whose list files must be purged

###########################  MAIN LOOP  #######################################
shopt -s nullglob

while read -r run; do
    [[ -z "$run" || ! "$run" =~ ^[0-9]+$ ]] && continue

    run_pad=$(printf "%08d" "$run")
    low=$(( (run/100)*100 ))
    high=$(( low + 100 ))
    range_dir="${BASE_DIR}/run_$(printf '%08d_%08d' "$low" "$high")"

    say "▸ Run ${CLR_BLD}${run}${CLR_RST}  (range dir: ${range_dir##*/})"

    # -------- 1. directory exists? ------------------------------------------
    if [[ ! -d "$range_dir" ]]; then
        warn "  Range directory NOT found – skipping"
        inc miss_dir;  cleanup_runs+=("$run")
        continue
    fi

    # -------- 2. collect ROOT files (primary dir) ---------------------------
    files=( "${range_dir}/${FILE_PREFIX}-${run_pad}-"*.root )

    # -------- 3. fallback one level higher if none --------------------------
    if (( ${#files[@]} == 0 )); then
        fallback_dir=$(dirname "$range_dir")
        files=( "${fallback_dir}/${FILE_PREFIX}-${run_pad}-"*.root )
        if (( ${#files[@]} == 0 )); then
            warn "  No ROOT files in ${range_dir##*/} or its parent – skipping"
            inc miss_files;  cleanup_runs+=("$run")
            continue
        else
            warn "  Found files one directory up (${fallback_dir##*/})"
        fi
    fi

    # -------- 4. write list file --------------------------------------------
    list_file="${LIST_DIR}/${FILE_PREFIX}-${run_pad}.list"
    printf "%s\n" "${files[@]}" > "$list_file"
    good "  Wrote ${#files[@]} paths → ${list_file##*/}"
    inc ok
done < "$RUNLIST_FILE"

###########################  CLEAN-UP  ########################################
for bad_run in "${cleanup_runs[@]}"; do
    bad_pad=$(printf "%08d" "$bad_run")
    stale="${LIST_DIR}/${FILE_PREFIX}-${bad_pad}.list"
    if [[ -e "$stale" ]]; then
        rm -f "$stale"
        ((removed++))
        warn "  Removed empty/stale list → ${stale##*/}"
    fi
done

############################  SUMMARY  ########################################
echo
good "Finished:"
say  "  Successful lists : ${ok}"
[[ $miss_dir   -gt 0 ]] && warn "  Skipped (no dir)  : ${miss_dir}"
[[ $miss_files -gt 0 ]] && warn "  Skipped (no files): ${miss_files}"
[[ $removed    -gt 0 ]] && warn "  Stale lists purged: ${removed}"
echo "${CLR_BLD}All done.${CLR_RST}"
###############################################################################
