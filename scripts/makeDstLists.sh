#!/usr/bin/env bash
###############################################################################
#  makeDstLists.sh
#  ---------------
#  Build .list files for sPHENIX Run‑24 Au+Au CALO‑DSTs  *or*
#  Run‑25 Au+Au jet‑DSTs, depending on the first argument:
#
#      ./makeDstLists.sh run24auau   # ← current behaviour (good‑run file)
#      ./makeDstLists.sh run25auau   # ← auto‑scan jet DST repository
#
#  • List files are written to  $LIST_DIR  (wiped at start)
#  • One .list per run, containing the full path(s) to matching ROOT files
###############################################################################
set -euo pipefail
IFS=$'\n\t'

########################  COLOUR / LOG HELPERS  ###############################
ESC=$'\e['
CLR_RED=${ESC}0\;31m  ; CLR_GRN=${ESC}0\;32m
CLR_YEL=${ESC}1\;33m  ; CLR_BLU=${ESC}1\;34m
CLR_BLD=${ESC}1m      ; CLR_RST=${ESC}0m

say()   { printf "${CLR_BLU}➜${CLR_RST} %s\n" "$*"; }
good()  { printf "${CLR_GRN}%s${CLR_RST}\n"   "$*"; }
warn()  { printf "${CLR_YEL}⚠ %s${CLR_RST}\n" "$*" >&2; }
fatal() { printf "${CLR_RED}✘ %s${CLR_RST}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO)"' ERR

############################  CONSTANTS  #####################################
SCRIPT_PWD="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"

MODE=${1:-run24auau}               # default keeps old behaviour
case "$MODE" in
    run24auau)
        # ---- Run‑24 CALO‑DST settings -------------------------------------
        FILE_PREFIX="DST_CALO_run2auau_new_2024p007"
        BASE_DIR="/sphenix/lustre01/sphnxpro/physics/run2auau/caloy2calib/new_2024p007"
        RUNLIST_FILE="${SCRIPT_PWD}/../goodRunList_sEPD_run24auau.txt"
        [[ -f "$RUNLIST_FILE" ]] || fatal "Run‑list file not found: $RUNLIST_FILE"
        ;;
    run25auau)
        # ---- Run‑25 JET‑DST settings --------------------------------------
        FILE_PREFIX="DST_JET"
        BASE_DIR="/sphenix/tg/tg01/jets/vbailey/run25_jet_dsts/new_newcdbtag_v005/jet"
        RUNLIST_FILE=""                 # not used in this mode
        [[ -d "$BASE_DIR" ]] || fatal "Jet‑DST directory not found: $BASE_DIR"
        ;;
    *)
        fatal "Unknown mode '$MODE' – use 'run24auau' or 'run25auau'"
        ;;
esac

############################  PRE‑FLIGHT  #####################################
say  "Mode                  : ${CLR_BLD}${MODE}${CLR_RST}"
say  "Destination .list dir : ${CLR_BLD}${LIST_DIR}${CLR_RST}"
say  "ROOT file base dir    : ${CLR_BLD}${BASE_DIR}${CLR_RST}"

# Purge previous list files
rm -f "${LIST_DIR:?}/"* 2>/dev/null || true
mkdir -p "$LIST_DIR"
good "Cleared old list files in ${LIST_DIR}"
echo

###########################  HELPERS  #########################################
# increment a nameref counter safely under 'set -e'
inc() { local -n ref=$1; ref=$(( ref + 1 )); }

###########################  BUILD RUN SET  ###################################
declare -a run_numbers=()

if [[ "$MODE" == run24auau ]]; then
    mapfile -t run_numbers < <(grep -E '^[0-9]+' "$RUNLIST_FILE" | sort -u)
else
    # Extract 8‑digit run numbers from file names once
    while IFS= read -r -d '' f; do
        [[ $f =~ ([0-9]{8}) ]] && run_numbers+=("${BASH_REMATCH[1]}")
    done < <(find "$BASE_DIR" -maxdepth 1 -name "${FILE_PREFIX}-*.root" -print0)
    run_numbers=($(printf '%s\n' "${run_numbers[@]}" | sort -u))
fi

total_runs=${#run_numbers[@]}
[[ $total_runs -eq 0 ]] && fatal "No runs found to process"
good "Runs to process       : $total_runs"
echo

###########################  MAIN LOOP  #######################################
ok=0 miss_dir=0 miss_files=0 removed=0
declare -a cleanup_runs=()

shopt -s nullglob

for run in "${run_numbers[@]}"; do
    run_pad=$(printf "%08d" "$run")

    if [[ "$MODE" == run24auau ]]; then
        low=$(( (run/100)*100 ))
        high=$(( low + 100 ))
        range_dir="${BASE_DIR}/run_$(printf '%08d_%08d' "$low" "$high")"
        search_dirs=("$range_dir" "$(dirname "$range_dir")")
    else
        search_dirs=("$BASE_DIR")       # jet DSTs live flat in one dir
    fi

    say "▸ Run ${CLR_BLD}${run}${CLR_RST}"

    files=()
    for d in "${search_dirs[@]}"; do
        [[ -d "$d" ]] || continue
        files+=( "$d/${FILE_PREFIX}-${run_pad}-"*.root )
        (( ${#files[@]} )) && break     # stop at first dir with matches
    done

    if (( ${#files[@]} == 0 )); then
        warn "  No ROOT files found – skipping"
        inc miss_files; cleanup_runs+=("$run"); continue
    fi

    list_file="${LIST_DIR}/${FILE_PREFIX}-${run_pad}.list"
    printf "%s\n" "${files[@]}" > "$list_file"
    good "  Wrote ${#files[@]} paths → ${list_file##*/}"
    inc ok
done

###########################  CLEAN‑UP (stale)  ################################
for bad_run in "${cleanup_runs[@]}"; do
    stale="${LIST_DIR}/${FILE_PREFIX}-$(printf '%08d' "$bad_run").list"
    if [[ -e "$stale" ]]; then
        rm -f "$stale"; inc removed
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

