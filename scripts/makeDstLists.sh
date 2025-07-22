#!/usr/bin/env bash
###############################################################################
#  makeDstLists.sh
#
#  Purpose
#  -------
#  Build one “*.list” file per run containing the absolute paths to DST ROOT
#  files.  These lists are subsequently consumed by downstream analysis or
#  GRID production tools.
#
#  High-level Workflow
#  -------------------
#    1.  Scan one or more DST repositories for files whose names encode the
#        run-number (e.g.   DST_JET-00068542-00017.root).
#    2.  For every run discovered, write a list file
#        $list_dir/<PREFIX>-<run>.list  containing all matching paths.
#    3.  (optional)  “caloFitting” workflow
#          • Retrieve Run-3 Au+Au run numbers that have CALOFITTING output.
#          • Perform a quick QA (runtime ≥ 5 min **and** GL1 events ≥ 1 × 10⁵).
#          • Keep the *golden* subset and create CALOFITTING list files.
#          • Optionally bypass CreateDstList.pl and build these lists by
#            traversing the raw production tree (see **forceFileList** mode).
#
#  Command-line Syntax
#  -------------------
#    ./makeDstLists.sh <mode> [caloFitting [forceFileList]]
#
#      <mode>            Required.  Selects the DST repositories to scan.
#                        ├─ run24auau   →  Run-24 Au+Au  CALO-DSTs
#                        └─ run25auau   →  Run-25 Au+Au  JET and JETCALO DSTs
#
#      caloFitting       Optional.  Enables the extended QA + CALOFITTING
#                        pipeline described above.
#
#      forceFileList     Optional **third** argument that is only honoured
#                        when *caloFitting* is also specified.  Instead of
#                        calling CreateDstList.pl, the script walks the
#                        directory
#                          /sphenix/lustre01/sphnxpro/production/run3auau/\
#                          physics/caloy2fitting/<tag>
#                        creates one list per run, and writes them to
#                        $list_dir with the same naming convention.
#
#  Typical Examples
#  ----------------
#    # Standard list building for Run-24 Au+Au CALO-DSTs
#    ./makeDstLists.sh run24auau
#
#    # Build Run-25 JET / JETCALO lists **plus** CALOFITTING QA & lists
#    ./makeDstLists.sh run25auau caloFitting
#
#    # Same as above, but force the CALOFITTING list files to be built from
#    # the production tree rather than via CreateDstList.pl
#    ./makeDstLists.sh run25auau caloFitting forceFileList
#
###############################################################################
set -euo pipefail
IFS=$'\n\t' ; shopt -s nullglob              # strict mode

########################  COLOUR / LOG HELPERS  ###############################
esc=$'\e['
clr_red=${esc}0\;31m ; clr_grn=${esc}0\;32m
clr_yel=${esc}1\;33m ; clr_blu=${esc}1\;34m
clr_bld=${esc}1m     ; clr_rst=${esc}0m

say()   { printf "${clr_blu}➜${clr_rst} %s\n" "$*"; }
good()  { printf "${clr_grn}%s${clr_rst}\n"   "$*"; }
warn()  { printf "${clr_yel}⚠ %s${clr_rst}\n" "$*" >&2; }
fatal() { printf "${clr_red}✘ %s${clr_rst}\n" "$*" >&2; exit 1; }

trap 'fatal "Script aborted (line $LINENO) – $BASH_COMMAND"' ERR

############################  CONSTANTS  #####################################
script_pwd="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
list_dir="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/dst_list"

###########################  ARGUMENT PARSING  ################################
mode=${1:-run24auau}
extra=${2:-}                       # e.g. "caloFitting" (case‑sensitive)

declare -A prefix_dir              # key = prefix, value = repository path

case "$mode" in
  run24auau)
    prefix_dir["DST_CALO_run2auau_new_2024p007"]="/sphenix/lustre01/sphnxpro/physics/run2auau/caloy2calib/new_2024p007"
    ;;
  run25auau)
    prefix_dir["DST_JET"]="/sphenix/tg/tg01/jets/vbailey/run25_jet_dsts/new_newcdbtag_v005/jet"
    prefix_dir["DST_JETCALO"]="/sphenix/tg/tg01/jets/vbailey/run25_jet_dsts/new_newcdbtag_v005/jetcalo"
    ;;
  *) fatal "Unknown mode “$mode” – use run24auau | run25auau" ;;
esac

#####################  PRE‑FLIGHT (repo existence, banner)  ###################
for p in "${prefix_dir[@]}"; do [[ -d "$p" ]] || fatal "Directory not found: $p"; done

say "Mode                  : ${clr_bld}${mode}${clr_rst}"
[[ -n $extra ]] && say "Extra workflow        : ${clr_bld}${extra}${clr_rst}"
say "Destination .list dir : ${clr_bld}${list_dir}${clr_rst}"
for k in "${!prefix_dir[@]}"; do say "Source for ${clr_bld}${k}${clr_rst}: ${prefix_dir[$k]}"; done

rm -rf "${list_dir:?}"/* 2>/dev/null || true
mkdir -p "$list_dir"
good "Output directory prepared"

###########################  COLLECT RUN IDs  #################################
declare -A run_set          # associative: run → 1 (dedup)

for pfx in "${!prefix_dir[@]}"; do
  base=${prefix_dir[$pfx]}
  while IFS= read -r -d '' f; do
    [[ $f =~ ([0-9]{8}) ]] && run_set[${BASH_REMATCH[1]}]=1
  done < <(find "$base" -type f -name "${pfx}-????????-*.root" -print0)
done

runs=("${!run_set[@]}") ; IFS=$'\n' runs=($(sort -n <<<"${runs[*]}")) ; IFS=$'\n\t'
(( ${#runs[@]} )) || fatal "No runs found to process"

good "Runs to process       : ${#runs[@]}"
echo

##############################  HELPERS  ######################################
inc() { local -n ref=$1; ref=$((ref+1)); }

############################  MAIN LOOP  ######################################
ok=0 miss=0 pruned=0
empty_lists=()

for run in "${runs[@]}"; do
  run8=$(printf "%08d" "$((10#$run))")
  say "▸ Run ${clr_bld}${run8}${clr_rst}"

  for pfx in "${!prefix_dir[@]}"; do
    base=${prefix_dir[$pfx]}
    mapfile -t files < <(find "$base" -type f -name "${pfx}-${run8}-*.root" -print)

    list_file="${list_dir}/${pfx}-${run8}.list"
    if (( ${#files[@]} )); then
      printf "%s\n" "${files[@]}" >"$list_file"
      good "  ${pfx}: wrote ${#files[@]} path(s)"
      inc ok
    else
      warn "  ${pfx}: no files found"
      empty_lists+=("$list_file") ; inc miss
    fi
  done
done

############################  CLEAN‑UP  #######################################
for lf in "${empty_lists[@]}"; do
  [[ -e "$lf" ]] || continue
  rm -f "$lf"; inc pruned
  warn "  Removed empty list $(basename "$lf")"
done

############################  SUMMARY  ########################################
echo
good "Finished (phase‑1 list building):"
say "  Successful lists : $ok"
[[ $miss   -gt 0 ]] && warn "  Runs with no files: $miss"
[[ $pruned -gt 0 ]] && warn "  Empty lists pruned: $pruned"

###############################################################################
#                        ░▒▓  EXTRA WORKFLOW ▓▒░                               #
###############################################################################
if [[ $extra == caloFitting ]]; then
  ###########################################################################
  # 0 . settings that you are likely to tweak occasionally                  #
  ###########################################################################
  tag="new_newcdbtag_v006"           # CreateDstList.pl  --tag
  dataset="run3auau"                 # CreateDstList.pl  --dataset
  calo_prefix="DST_CALOFITTING"
  min_runtime=300                    # s  (≥ 5 min)
  min_gl1_evt=100000                 # GL1 events cut
  top_trig=20                        # show N most‑frequent triggers
  run3_list="run3auau-${tag}.list"
  golden_txt="../run3GoldenRuns.txt"

  ###########################################################################
  # 1 . obtain Run‑3 Au+Au run numbers for CALOFITTING                      #
  ###########################################################################
  say "Collecting Run‑3 run numbers via CreateDstList.pl"
  CreateDstList.pl --tag "$tag" --dataset "$dataset" --printruns "$calo_prefix" \
                   >"$run3_list"
  good "Run‑list ➔ $run3_list  ($(wc -l <"$run3_list") runs)"

  ###########################################################################
  # 2 . quick QA (run‑time, GL1 counts, trigger scalers)                    #
  ###########################################################################
  PSQL=(psql -h sphnxdaqdbreplica -d daq -At -F $'\t' -q)
  sql()          { "${PSQL[@]}" -c "$1" 2>/dev/null || true; }
  num_or_zero()  { [[ $1 =~ ^-?[0-9]+$ ]] && printf '%s' "$1" || printf 0; }
  safe_mapfile() { mapfile -t "$@" || true; }

  mapfile -t run3Runs <"$run3_list"
  ((${#run3Runs[@]})) || fatal "No run3 runs retrieved from CreateDstList.pl"

  hdr=$'Run\trt[s]\tGL1ev\tTrigON\tTrigOFF\tFail\tNote'
  rows=(); golden=(); tot_rt=0; tot_ev=0; fail_cnt=0
  declare -A trig_run trig_sum trig_sd_sum trig_sd_cnt

  say "Querying DAQ DB …"
  for run in "${run3Runs[@]}"; do
    rt=$(sql "SELECT FLOOR(EXTRACT(EPOCH FROM (ertimestamp-brtimestamp)))::INT
              FROM run WHERE runnumber=$run;" | tr -d ' ') ; rt=$(num_or_zero "$rt")
    ev=$(sql "SELECT COALESCE(SUM(lastevent-firstevent+1),0)::BIGINT
              FROM filelist WHERE runnumber=$run
               AND filename LIKE '%GL1_physics_gl1daq%.evt';" | tr -d ' ')
    ev=$(num_or_zero "$ev")

    # -- trigger scalers + scaledown -----------------------------------------
    # We now read three columns:  triggername <TAB> index <TAB> scaled
    # and, for each index, fetch scaledownNN from gl1_scaledown.
    on=0  off=0
    while IFS=$'\t' read -r trg idx scaled; do
        [[ -z $trg ]] && continue                 # safety

        scaled=$(num_or_zero "$scaled")           # normalise
        idx=$(num_or_zero "$idx")                 # safety for scaledown query

        # obtain the scaledown factor for this bit in the current run
        sd=$(sql "SELECT scaledown${idx} FROM gl1_scaledown WHERE runnumber = $run;" | tr -d ' ')
        sd=$(num_or_zero "$sd")

        if (( scaled > 0 )); then
            (( ++on ))
            trig_run["$trg"]=$(( ${trig_run["$trg"]:-0} + 1 ))
            trig_sum["$trg"]=$(( ${trig_sum["$trg"]:-0} + scaled ))

            # only count *valid* scaledown factors ( > 0 )
            if (( sd > 0 )); then
                trig_sd_sum["$trg"]=$(( ${trig_sd_sum["$trg"]:-0} + sd ))
                trig_sd_cnt["$trg"]=$(( ${trig_sd_cnt["$trg"]:-0} + 1 ))
            fi
        else
            (( ++off ))
        fi
        
    done < <(sql "SELECT t.triggername, t.index, s.scaled
                   FROM gl1_scalers  s
                   JOIN gl1_triggernames t
                     ON s.index = t.index
                    AND s.runnumber BETWEEN t.runnumber AND t.runnumber_last
                  WHERE s.runnumber = $run;")



    note="GOOD"
    (( rt <  min_runtime )) && { note="SHORT"; }
    (( ev <  min_gl1_evt ))  && { note="${note/_/}+LOWEVT"; }
    [[ $note == GOOD ]] && golden+=("$run") || ((++fail_cnt))

    rows+=("$run"$'\t'"$rt"$'\t'"$ev"$'\t'"$on"$'\t'"$off"$'\t'"$note")
    (( tot_rt+=rt, tot_ev+=ev ))
  done

  ###########################################################################
  # 3 . print nice tables (per-run, cut-flow, scalers)                      #
  ###########################################################################
  echo -e "\n${clr_grn}================ Per-run summary ================${clr_rst}"
  printf '%s\n' "$hdr" "${rows[@]}" | column -t -s $'\t'

  ########################  CUT-FLOW BREAK-DOWN  ###############################
  # ── 1. numbers for the *full* sample ----------------------------------------
  all_runs=${#run3Runs[@]}
  all_rt=$tot_rt
  all_ev=$tot_ev

  # ── 2. after run-time ≥ ${min_runtime}s --------------------------------------
  timeRuns=0; timeRt=0; timeEv=0
  while IFS=$'\t' read -r _ rt ev _ _ note; do
      [[ $note == SHORT* ]] && continue        # drop the short ones
      ((timeRuns++, timeRt+=rt, timeEv+=ev))
  done <<<"$(printf '%s\n' "${rows[@]}")"

  # ── 3. after GL1 events ≥ ${min_gl1_evt} -------------------------------------
  eventRuns=0; eventRt=0; eventEv=0
  while IFS=$'\t' read -r _ rt ev _ _ note; do
      [[ $note != GOOD ]] && continue          # keep only the golden ones
      ((eventRuns++, eventRt+=rt, eventEv+=ev))
  done <<<"$(printf '%s\n' "${rows[@]}")"

  echo -e "\n${clr_grn}================ Cut-flow summary ===============${clr_rst}"
  printf '%-12s %6s %12s %12s\n' "Stage" "Runs" "GL1evt" "Runtime[h]"
  printf '%-12s %6d %12d %12.2f\n' "All"     "$all_runs"   "$all_ev"   "$(bc -l <<<"$all_rt/3600")"
  printf '%-12s %6d %12d %12.2f\n' "Time≥5m" "$timeRuns"   "$timeEv"   "$(bc -l <<<"$timeRt/3600")"
  printf '%-12s %6d %12d %12.2f\n' "Evt≥1e5" "$eventRuns"  "$eventEv"  "$(bc -l <<<"$eventRt/3600")"

  # keep_cnt is used later when we write the .list files ------------------------
  keep_cnt=$eventRuns

  ########################  TRIGGER TOTALS  ####################################
  echo -e "\n${clr_grn}================ Trigger totals =================${clr_rst}"

  # dynamic padding – find the longest trigger name and add 2 chars
  maxlen=7
  for trg in "${!trig_run[@]}"; do
      (( ${#trg} > maxlen )) && maxlen=${#trg}
  done
  (( maxlen += 2 ))

  hdr_fmt="%-${maxlen}s %8s %15s %10s\n"
  printf "$hdr_fmt" "Trigger" "ONruns" "ΣScaled" "AvgSD"

  for trg in "${!trig_run[@]}"; do
      cnt=${trig_sd_cnt[$trg]:-0}
      if (( cnt )); then
          avg=$(bc -l <<<"${trig_sd_sum[$trg]}/${cnt}")
          avg=$(printf "%.1f" "$avg")
      else
          avg="–"
      fi
      printf "$hdr_fmt" "$trg" "${trig_run[$trg]}" "${trig_sum[$trg]}" "$avg"
  done | sort -k2 -nr | head -"$top_trig"



  ###########################################################################
  # 4 . produce per‑run .list files for the golden sample                  #
  ###########################################################################

  # third positional argument decides how the .list files are built
  build_mode=${3:-create}        # create  →  use CreateDstList.pl   (default)
                                   # forceFileList → traverse caloy2fitting tree

  say "Generating ${calo_prefix} .list files for ${#golden[@]} golden runs  (mode=${build_mode})"

  if [[ $build_mode == forceFileList ]]; then
      # ----------------------------------------------------------------------
      # 4A. DIRECTORY TRAVERSAL BACK‑UP – build the list files ourselves
      #     Base      : /sphenix/lustre01/sphnxpro/production/run3auau/physics/caloy2fitting/<tag>
      #     Sub‑dirs  : run_<000NNN00>_<000NNN00> (100‑run buckets)
      #     Pattern   : DST_CALOFITTING_<dataset>_<tag>-<run8>-NNNNN.root
      # ----------------------------------------------------------------------
      calo_base="/sphenix/lustre01/sphnxpro/production/run3auau/physics/caloy2fitting/${tag}"

      for run in "${golden[@]}"; do
        run_dec=$((10#$run))                                  # strip any octal :contentReference[oaicite:0]{index=0}
        run8=$(printf "%08d" "$run_dec")                      # 8‑digit run number

        bucket_start=$(( (run_dec/100)*100 ))                 # integer division :contentReference[oaicite:1]{index=1}
        bucket_end=$(( bucket_start + 100 ))
        bucket_dir=$(printf "run_%08d_%08d" "$bucket_start" "$bucket_end")

        full_dir="${calo_base}/${bucket_dir}"
        [[ -d $full_dir ]] || { warn "  ${run}: directory ${bucket_dir} missing"; continue; }

        # collect all segments for this run, sort them numerically for reproducibility
        mapfile -t segs < <(find "$full_dir" -type f \
                     -name "${calo_prefix}_${dataset}_${tag}-${run8}-*.root" \
                     | sort -V)

        out_list="${list_dir}/${calo_prefix}_${dataset}_${tag}-${run8}.list"

        if (( ${#segs[@]} )); then
          printf '%s\n' "${segs[@]}" >"$out_list"
          good "  ${run}: wrote ${#segs[@]} path(s)"
        else
          warn "  ${run}: no CALOFITTING files found in ${bucket_dir}"
        fi
  done

  else
      # ----------------------------------------------------------------------
      # 4B. STANDARD PATH – use CreateDstList.pl exactly as before
      # ----------------------------------------------------------------------
      for run in "${golden[@]}"; do
          out_list="${list_dir}/${calo_prefix}-${run}.list"
          CreateDstList.pl --tag "$tag" --dataset "$dataset" \
                           --list "$run3_list" --run "$run" \
                           "$calo_prefix" >"$out_list"
          good "  ${run}: $(wc -l <"$out_list") path(s)"
      done
  fi

  printf '%s\n' "${golden[@]}" >"$golden_txt"
  good "Golden run list ➔ $golden_txt"

  say "CALOFITTING workflow complete."
fi

echo "${clr_bld}All done.${clr_rst}"
###############################################################################
