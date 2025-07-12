#!/usr/bin/env bash
###############################################################################
#  makeDstLists.sh  – build one “.list” per run with absolute paths to DST
#  ROOT files.
#
#  Modes:
#    run24auau → Run-24 Au+Au CALO-DSTs
#    run25auau → Run-25 Au+Au JET **and** JETCALO-DSTs (parallel)
#
#  No external run-number files are needed – runs are detected by scanning the
#  filenames themselves.
###############################################################################
set -euo pipefail
IFS=$'\n\t'

########################  COLOUR / LOG HELPERS  ###############################
esc=$'\e['
clr_red=${esc}0\;31m ; clr_grn=${esc}0\;32m
clr_yel=${esc}1\;33m ; clr_blu=${esc}1\;34m
clr_bld=${esc}1m     ; clr_rst=${esc}0m

say()   { printf "${clr_blu}➜${clr_rst} %s\n" "$*"; }
good()  { printf "${clr_grn}%s${clr_rst}\n"   "$*"; }
warn()  { printf "${clr_yel}⚠ %s${clr_rst}\n" "$*" >&2; }
fatal() { printf "${clr_red}✘ %s${clr_rst}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO)"' ERR

############################  CONSTANTS  #####################################
script_pwd="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
list_dir="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"

###########################  MODE SELECTION  ##################################
mode=${1:-run24auau}
declare -A prefix_dir     # key = prefix, value = base directory

case "$mode" in
  run24auau)
    # Run-24 CALO-DST repository (nested by run-range dirs)
    prefix_dir["DST_CALO_run2auau_new_2024p007"]="/sphenix/lustre01/sphnxpro/physics/run2auau/caloy2calib/new_2024p007"
    ;;

  run25auau)
    # Run-25 flat repositories – JET and JETCALO in parallel
    prefix_dir["DST_JET"]="/sphenix/tg/tg01/jets/vbailey/run25_jet_dsts/new_newcdbtag_v005/jet"
    prefix_dir["DST_JETCALO"]="/sphenix/tg/tg01/jets/vbailey/run25_jet_dsts/new_newcdbtag_v005/jetcalo"
    ;;

  *) fatal "Unknown mode “$mode” – use run24auau | run25auau" ;;
esac

# Check source directories exist
for d in "${prefix_dir[@]}"; do
  [[ -d "$d" ]] || fatal "Directory not found: $d"
done

############################  PRE-FLIGHT  #####################################
say "Mode                  : ${clr_bld}${mode}${clr_rst}"
say "Destination .list dir : ${clr_bld}${list_dir}${clr_rst}"
for k in "${!prefix_dir[@]}"; do
  say "Source for ${clr_bld}${k}${clr_rst}: ${prefix_dir[$k]}"
done

rm -rf "${list_dir:?}"/* 2>/dev/null || true
mkdir -p "$list_dir"
good "Output directory prepared"

###########################  COLLECT RUN IDs  #################################
declare -A run_set   # associative (run → 1) to deduplicate

for pfx in "${!prefix_dir[@]}"; do
  base=${prefix_dir[$pfx]}
  while IFS= read -r -d '' f; do
    [[ $f =~ ([0-9]{8}) ]] && run_set[${BASH_REMATCH[1]}]=1
  done < <(find "$base" -type f -name "${pfx}-????????-*.root" -print0)
done

runs=("${!run_set[@]}")
IFS=$'\n' runs=($(sort -n <<<"${runs[*]}")); IFS=$'\n\t'

(( ${#runs[@]} )) || fatal "No runs found to process"
good "Runs to process       : ${#runs[@]}"
echo

##############################  HELPERS  ######################################
inc() { local -n ref=$1; ref=$((ref+1)); }

############################  MAIN LOOP  ######################################
ok=0 miss=0 pruned=0
empty_lists=()

shopt -s nullglob

for run in "${runs[@]}"; do
  run_dec=$((10#$run))
  run8=$(printf "%08d" "$run_dec")
  say "▸ Run ${clr_bld}${run8}${clr_rst}"

  for pfx in "${!prefix_dir[@]}"; do
    base=${prefix_dir[$pfx]}
    mapfile -t files < <(find "$base" -type f -name "${pfx}-${run8}-*.root" -print)

    list_file="${list_dir}/${pfx}-${run8}.list"

    if (( ${#files[@]} )); then
      printf "%s\n" "${files[@]}" > "$list_file"
      good "  ${pfx}: wrote ${#files[@]} path(s)"
      inc ok
    else
      warn "  ${pfx}: no files found"
      empty_lists+=("$list_file"); inc miss
    fi
  done
done

############################  CLEAN-UP  #######################################
for lf in "${empty_lists[@]}"; do
  [[ -e "$lf" ]] || continue
  rm -f "$lf"; inc pruned
  warn "  Removed empty list $(basename "$lf")"
done

############################  SUMMARY  ########################################
echo
good "Finished:"
say "  Successful lists : $ok"
[[ $miss   -gt 0 ]] && warn "  Runs with no files: $miss"
[[ $pruned -gt 0 ]] && warn "  Empty lists pruned: $pruned"
echo "${clr_bld}All done.${clr_rst}"
###############################################################################
