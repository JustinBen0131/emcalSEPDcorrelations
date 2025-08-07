#!/usr/bin/env bash
##############################################################################
#  runAuAuStages.sh  –  Au+Au Run‑24/25 **end‑to‑end QA driver**
#
#  Quick command map
#  ─────────────────
#  • stage0              : launch one Condor job **per run** (per‑run QA)
#  • stage1              : scan SEB/HCal, write *good‑run list* + bar chart
#  • stage2              : Condor “hadd” groups of ≤10 good runs each
#  • stage3 [mode] [...] : final merge **+ combined QA**
#        mode            : local  – merge/QA on login node
#                          condor – submit merge/QA to Condor (default)
#        skipStage2      : build combined file directly from good‑run list
#  • processOnly [mode] [qaList]
#        mode            : local  | condor   (default = local)
#        qaList          : comma‑separated subset of QA modules (see table)
#                          leave empty to run the **full** QA suite
#
#  Example cheat‑sheet
#  ───────────────────
#    ./runAuAuStages.sh stage0                    # full per‑run grid
#    ./runAuAuStages.sh stage1                    # build good‑run list
#    ./runAuAuStages.sh stage2                    # group hadd on Condor
#    ./runAuAuStages.sh stage3                    # final combined QA (Condor)
#    ./runAuAuStages.sh stage3 local skipStage2   # local merge, skip stage‑2
#    ./runAuAuStages.sh processOnly               # re‑run all QA on combined
#    ./runAuAuStages.sh processOnly correlations  # correlations QA only
#    ./runAuAuStages.sh processOnly condor pi0,jetqa
#
#  QA module keywords accepted by `processOnly`
#  ─────────────────────────────────────────────
#  | Keyword        | Module executed                    |
#  | -------------- | ---------------------------------- |
#  | correlations   | CALO × sEPD × MBD correlations     |
#  | hcal           | HCal (IHCal / OHCal / total)       |
#  | mbd            | MBD QA                             |
#  | sepd           | sEPD event‑plane / tile QA         |
#  | sepdother      | sEPD miscellaneous QA              |
#  | jetqa          | Jet QA (general & summary)         |
#  | eventqa        | Global event QA                    |
#  | triggerqa      | Trigger counters / rates           |
#  | pi0            | π⁰ invariant‑mass QA               |
#  | emcal          | EMCal occupancy / spectra          |
#  | vn             | vₙ analysis QA                      |
#
#  Pass multiple keywords comma‑separated, e.g.
#       QA_ONLY="correlations,pi0,jetqa"
#  (the wrapper sets this automatically for you when you use
#   `./runAuAuStages.sh processOnly [...] <list>`).
##############################################################################
set -euo pipefail

# ─────────  Verbosity control (same semantics as runAuAu.sh)  ─────────
# VERBOSE=0   → silent (default)
# VERBOSE=1   → extra informational messages
# VERBOSE=2   → shell trace (set -x) + all VERBOSE=1 output
# The first positional token ‘-v’ or ‘--verbose’ also forces VERBOSE=1.
# ----------------------------------------------------------------------
if [[ ${1:-} == VERBOSE=* ]]; then
    export VERBOSE="${1#VERBOSE=}"      # take numeric value after '='
    shift                               # drop this pseudo‑argument
fi
: "${VERBOSE:=0}"                       # default when nothing supplied

if [[ ${1:-} == "-v" || ${1:-} == "--verbose" ]]; then
    export VERBOSE=1
    shift
fi

# enable full shell trace for the most chatty level
if (( VERBOSE >= 2 )); then
    set -x
fi

# ───────────────  global paths (adapt only these when relocating) ───────────
PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
INPUT_DIR="${PROJECT_BASE}/output"                         # per‑run ROOTs
OUTPUT_DIR="/sphenix/tg/tg01/bulk/jbennett/GLOBAL_QA/outputPlots"
EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh" # already existing

TMP_BASE="${PROJECT_BASE}/tmp_online_stages"        # scratch for stages 1‑3
RUNLIST="${TMP_BASE}/runlist_stage2.txt"            # output of stage 1
GROUP_DIR="${TMP_BASE}/groups"                      # group hadds (stage 2)
mkdir -p "${TMP_BASE}" "${GROUP_DIR}"

# ───────────────  Condor bookkeeping (stage0) ───────────────
SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"
LOG_DIR="${PROJECT_BASE}/log"
STDOUT_DIR="${PROJECT_BASE}/stdout"
STDERR_DIR="${PROJECT_BASE}/error"

# ───────────────  helpers – pretty loggers ───────────────
clrB=$'\033[1;34m'; clrR=$'\033[1;31m'; clr0=$'\033[0m'
note(){ printf "${clrB}[INFO]${clr0}  %s\n" "$*"; }
warn(){ printf "${clrR}[WARN]${clr0}  %s\n" "$*"; }
die (){ printf "${clrR}[FATAL]${clr0} %s\n" "$*" >&2; exit 2; }
step(){ printf "\n${clrB}==========  %s  ==========${clr0}\n" "$*"; }

[[ -x "${EXEC_WRAPPER}" ]] || die "wrapper ${EXEC_WRAPPER} not found"

stage0(){   # optional arg: rescueBusy
  local action="${1:-}"

  # ───────────────────────────────────────────────────────────────────
  # rescueBusy: kill Stage‑0 Condor jobs (runAuAuExecutable.sh) and
  #             process those runs locally, sequentially
  # ───────────────────────────────────────────────────────────────────
  if [[ "$action" == rescueBusy ]]; then
    step "Stage‑0  –  rescueBusy (remove Condor per‑run jobs and run locally)"

    # 1) discover busy Stage‑0 runs from Condor
    declare -a busyRuns=()
    while read -r line || [[ -n "$line" ]]; do
      [[ -z "$line" ]] && continue
      # Args look like: "<INPUT_DIR>/output_XXXXXXXX.root  <OUTPUT_DIR>/XXXXXXXX"
      # Prefer extracting from the first token (input ROOT path)
      rf="${line%% *}"
      if [[ "$rf" =~ output_([0-9]{5,8})\.root ]]; then
        busyRuns+=( "${BASH_REMATCH[1]}" )
      elif [[ "$line" =~ /([0-9]{5,8})([[:space:]]|$) ]]; then
        busyRuns+=( "${BASH_REMATCH[1]}" )
      fi
    done < <(
      condor_q "$USER" \
        -constraint 'regexp("runAuAuExecutable.sh",Cmd) && (JobStatus == 1 || JobStatus == 2)' \
        -af Args 2>/dev/null
    ) || true

    if (( ${#busyRuns[@]} == 0 )); then
      note "No active Stage‑0 jobs to rescue."
      return
    fi

    # dedupe + show
    mapfile -t busyRuns < <(printf '%s\n' "${busyRuns[@]}" | sort -u)
    note "Active Stage‑0 run(s): $(printf '%s ' "${busyRuns[@]}")"

    # 2) remove those Condor jobs (do not touch any files)
    local busyExpr='regexp("runAuAuExecutable.sh",Cmd) && (JobStatus == 1 || JobStatus == 2)'
    local nKill
    nKill=$(condor_q "$USER" -constraint "$busyExpr" -af ClusterId ProcId 2>/dev/null | wc -l)
    if (( nKill > 0 )); then
      note "Removing $nKill Condor job(s) for Stage‑0"
      condor_rm "$USER" -constraint "$busyExpr" || warn "condor_rm returned non‑zero – continuing anyway"
    else
      note "No matching Condor jobs to remove"
    fi

    # 3) run the same work locally, one run after another
    mkdir -p "${OUTPUT_DIR}"
    for run in "${busyRuns[@]}"; do
      local rf="${INPUT_DIR}/output_${run}.root"
      local outdir="${OUTPUT_DIR}/${run}"
      if [[ ! -f "$rf" ]]; then
        warn "Input ROOT not found for run ${run} → ${rf}  (skipped)"
        continue
      fi
      mkdir -p "$outdir"
      note "Local Stage‑0 processing for run ${run}"
      "${EXEC_WRAPPER}" "$rf" "$outdir" || warn "Local Stage‑0 processing failed for run ${run}"
    done

    note "Stage‑0 rescueBusy finished."
    return
  fi

  # ───────────────────────────────────────────────────────────────────
  # Normal Stage‑0 submission path
  # ───────────────────────────────────────────────────────────────────
  step "Stage‑0  –  submit run‑by‑run Condor jobs"
  rm -rf "${SUBMIT_DIR:?}/"* 2>/dev/null || true
  mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}" "${OUTPUT_DIR}"

  mapfile -t roots < <(ls "${INPUT_DIR}"/output_*.root 2>/dev/null | sort)
  (( ${#roots[@]} )) || die "No per‑run ROOT files in ${INPUT_DIR}"

  for rf in "${roots[@]}"; do
      run=$(basename "${rf}" .root); run="${run#output_}"
      sub="${SUBMIT_DIR}/${run}.sub"
      cat >"${sub}" <<EOF
universe      = vanilla
executable    = ${EXEC_WRAPPER}
arguments     = ${rf}  ${OUTPUT_DIR}/${run}
log           = ${LOG_DIR}/${run}.log
output        = ${STDOUT_DIR}/${run}.out
error         = ${STDERR_DIR}/${run}.err
getenv        = True
request_memory= 1.5GB
+JobFlavour   = "tomorrow"
queue
EOF
      condor_submit "${sub}" || die "condor_submit failed for run ${run}"
  done
  note "All ${#roots[@]} run‑jobs submitted."
}


# ───────────────────────  stage-1  –  SEB / HCal diagnostics  ──────────────
stage1(){
  step "Stage-1  – SEB scan and good-run list"

  rm -f "${RUNLIST}"
  mkdir -p "${OUTPUT_DIR}/Combined"

  # --------  one-shot, *very* chatty C++ macro  ----------------------------
  cat >"${TMP_BASE}/stage1_seb.C"<<'CPP'
/*  Verbose SEB + HCal scanner
    ──────────────────────────
    • traverses every output_XXXXXXXX.root in <inDir>
    • shows per-run “MissingSEB.txt” status
    • builds a good-run list              → <runList>
    • draws MissingSEB_distribution.png   → <plotDir>/Combined
    • prints the full SEB summary table identical to the production macro
*/
#include <filesystem>
#include <fstream>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <regex>
#include <TH1.h>
#include <TCanvas.h>
namespace fs = std::filesystem;

void stage1_seb(const char* inDir,
                const char* plotDir,
                const char* runList)
{
   std::cout << "\n>>> Traversing directory   " << inDir << "\n";

   std::map<std::string,int>               sebCnt;     // SEB## → incidence
   std::map<std::string,std::vector<std::string>> badRunMap;
   std::vector<fs::path>                   runFiles;

   for (auto& e : fs::directory_iterator(inDir))
       if (e.is_regular_file()) runFiles.push_back(e.path());
   std::sort(runFiles.begin(), runFiles.end());

   std::ofstream good(runList);
   if(!good){ std::cerr<<"[FATAL] cannot write runList\n"; return; }

   /* ── Step-by-step over every ROOT file ─────────────────────────────── */
   for (const auto& f : runFiles) {
       std::string fname = f.filename().string();
       if (fname == "output_ALL_COMBINED.root") continue;

       std::smatch m;
       if (!std::regex_search(fname, m,
                              std::regex(R"(output_([0-9]{8})\.root)")))
           continue;
       std::string run = m[1].str();
       fs::path miss = fs::path(plotDir) / run / "MissingSEB.txt";

       std::cout << "  → inspecting " << fname << "  … ";

       if (fs::exists(miss) && fs::file_size(miss) > 0) {
           std::ifstream fMiss(miss);
           std::string tok;                     // header token (discard)
           fMiss >> tok;
           while (fMiss >> tok) {
               if (tok.rfind("SEB",0)!=0) continue;
               sebCnt[tok]++;  badRunMap[run].push_back(tok);
           }
           std::cout << "⚠️  missing-SEB file found (" << badRunMap[run].size()
                     << " SEB)\n";
       } else {
           std::cout << "✓ no missing SEB – good run\n";
           good << run << '\n';
       }
   }
   good.close();
   std::cout << "\n>>> Good-run list written to " << runList << "\n";

   /* ── Draw bar chart ────────────────────────────────────────────────── */
   TH1I h("h","Runs with missing SEB;SEB index;Number of runs",16,-0.5,15.5);
   for (int i=0;i<16;++i){
       char lab[8]; sprintf(lab,"SEB%02d",i);
       h.GetXaxis()->SetBinLabel(i+1,lab);
       h.SetBinContent(i+1, sebCnt[lab]);
   }
   h.SetFillColor(kAzure+1); h.SetBarWidth(0.8); h.SetBarOffset(0.1);
   TCanvas c("c","",800,500); c.SetGridy(); h.Draw("bar2");
   fs::create_directories(fs::path(plotDir)/"Combined");
   fs::path png = fs::path(plotDir)/"Combined"/"MissingSEB_distribution.png";
   c.SaveAs(png.c_str());
   std::cout << ">>> Bar chart saved to      " << png << "\n";

   /* ── Summary table (verbatim format of production macro) ───────────── */
   const std::size_t nBad = badRunMap.size();
   std::size_t nBad1=0, nBadMul=0;
   for (const auto& [_,v] : badRunMap) (v.size()==1 ? ++nBad1 : ++nBadMul);

   std::cout << "\n==========  Missing SEB summary  ==========\n";
   std::cout << "Runs with ≥1 missing SEB : " << nBad  << "\n"
             << "   ├─ exactly one SEB    : " << nBad1 << "\n"
             << "   └─ multiple SEBs      : " << nBadMul << "\n";
   std::cout << "SEB │ Runs\n----+-----\n";
   for (int i=0;i<16;++i){
       std::string key = "SEB"+std::to_string(i);
       int cnt = sebCnt.count(key)?sebCnt[key]:0;
       std::cout << std::left << std::setw(3) << key
                 << " │ " << cnt << "\n";
   }
   std::cout << "===========\n";
}
CPP

  # --- compile + execute; keep all stdout so the user sees the chatter ----
  root -l -b -q \
       -e "gSystem->SetBuildDir(\"${TMP_BASE}\",kTRUE)" \
       "${TMP_BASE}/stage1_seb.C+(\"${INPUT_DIR}\",\"${OUTPUT_DIR}\",\"${RUNLIST}\")" \
  || die "Stage-1 C++ macro failed"

  note "Stage-1 completed – detailed log printed above."
}


# ───────────────────────  stage-2  –  group hadd (Condor)  ────────────────
##
##  submit_group_hadd idx outfile infiles…
##  • writes a HTCondor submit file that directly calls `/usr/bin/hadd`
##  • submits that job and prints the cluster ID
##
submit_group_hadd() {
    local idx="$1"; shift
    local outFile="$1"; shift
    local -a inFiles=( "$@" )

    printf "      ↳ crafting group %02d  (%d file%s)\n" \
           "${idx}" "${#inFiles[@]}" "$([[ ${#inFiles[@]} -eq 1 ]] && echo "" || echo "s")"

    local submit="${TMP_BASE}/group_${idx}.sub"
    cat > "${submit}" <<EOF
universe        = vanilla
executable      = $(command -v hadd)
arguments       = -f ${outFile} ${inFiles[*]}
log             = ${TMP_BASE}/group_${idx}.log
output          = ${TMP_BASE}/group_${idx}.out
error           = ${TMP_BASE}/group_${idx}.err
getenv          = True
request_memory  = 1GB
+JobFlavour     = "tomorrow"
queue
EOF

    echo "      ↳ submitting group-${idx} …"
    local out
    if ! out=$(condor_submit "${submit}" 2>&1); then
        die "condor_submit failed for group ${idx}\n${out}"
    fi
    # extract "N job(s) submitted to cluster K"
    local jobs
    jobs="$(printf '%s\n' "${out}" | grep -Eo '[0-9]+ job(s)? submitted to cluster [0-9]+' || true)"
    echo "         ${jobs:-[no-confirmation-string]}"
}

stage2() {
    step "Stage-2  –  split good runs into ≤10-file batches and submit hadd jobs"

    [[ -f "${RUNLIST}" ]] || die "Run-list ${RUNLIST} missing – run stage1 first"
    mapfile -t runs < "${RUNLIST}"
    (( ${#runs[@]} )) || die "Run-list is empty"
    mkdir -p "${GROUP_DIR}"

    echo "    total good runs            : ${#runs[@]}"
    echo "    files per Condor hadd job   : 10"
    echo "    target directory for groups : ${GROUP_DIR}"
    echo

    local grpIdx=0
    local -a chunk=()

    for run in "${runs[@]}"; do
        local f="${INPUT_DIR}/output_${run}.root"
        if [[ ! -f "${f}" ]]; then
            warn "      ‼︎ ROOT file for run ${run} not found – skipped"
            continue
        fi

        chunk+=( "$f" )
        printf "      • queued %-14s (grp %02d, slot %d/10)\n" \
               "$(basename "$f")" "${grpIdx}" "${#chunk[@]}"

        if (( ${#chunk[@]} == 10 )); then
            echo "      → group complete – submit grp-${grpIdx}"
            submit_group_hadd "${grpIdx}" \
                               "${GROUP_DIR}/group_${grpIdx}.root" \
                               "${chunk[@]}"
            chunk=()
            (( ++grpIdx ))
        fi
    done

    # final partial
    if (( ${#chunk[@]} )); then
        echo "      → final partial group (${#chunk[@]} file(s)) – submit grp-${grpIdx}"
        submit_group_hadd "${grpIdx}" \
                           "${GROUP_DIR}/group_${grpIdx}.root" \
                           "${chunk[@]}"
    fi

    rm -f "${RUNLIST}"
    echo
    note "Stage-2 finished – all group-hadd jobs are now in the Condor queue."
}

# ───────────────────────  stage‑3  –  final hadd + QA  ─────────────────────
submit_final_hadd(){   # outfile  infile...
  local ofile="$1"; shift
  local sub="${TMP_BASE}/final_hadd.sub"
  cat >"${sub}"<<EOF
universe      = vanilla
executable    = /bin/bash
arguments     = -c "hadd -f ${ofile} $* && \
                    export COMBINED_ONLY=1 EXTERNAL_HADD=1 && \
                    ${EXEC_WRAPPER} --final"
log           = ${TMP_BASE}/final_hadd.log
output        = ${TMP_BASE}/final_hadd.out
error         = ${TMP_BASE}/final_hadd.err
getenv        = True
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOF
  condor_submit "${sub}" || die "condor_submit failed (final hadd)"
}

# ---------- helper: build list of input ROOT files --------------------------
# When skipStage2 is requested, we need to convert the good‑run list produced
# in Stage‑1 into full file paths under $INPUT_DIR.
collect_run_files() {          # arg 1 = run‑list file
    local list="$1"
    mapfile -t runs < "$list" || return 1
    (( ${#runs[@]} )) || return 1

    local f; local -a files=()
    for r in "${runs[@]}"; do
        f="${INPUT_DIR}/output_${r}.root"
        [[ -f "$f" ]] && files+=( "$f" ) || \
            warn "‼︎ run ${r} listed as good but file missing – skipped"
    done
    (( ${#files[@]} )) || return 1
    printf '%s\n' "${files[@]}"
}

# ---------- Stage‑3 – final merge + combined QA -----------------------------
# Usage:
#   ./runAuAuStages.sh stage3                  → normal path (needs Stage‑2)
#   ./runAuAuStages.sh stage3 skipStage2       → Condor merge, skip Stage‑2
#   ./runAuAuStages.sh stage3 local skipStage2 → local  merge, skip Stage‑2
stage3(){   # [local|condor]  [skipStage2]
  local mode="condor"
  local use_runlist="false"

  # first optional token
  case "$1" in
      local|condor)  mode="$1"; shift ;;
  esac
  # second optional token
  [[ ${1:-} == skipStage2 ]] && { use_runlist="true"; shift; }

  step "Stage‑3  –  final hadd & QA (${mode}${use_runlist:+, skipping Stage‑2})"

  # --------------------------------------------------------------------------
  # 1. Build input‑file array  →  inputs[@]
  # --------------------------------------------------------------------------
  declare -a inputs
  if [[ "$use_runlist" == "true" ]]; then
      [[ -f "$RUNLIST" ]] || die "Run‑list $RUNLIST missing – run stage1 first"
      mapfile -t inputs < <(collect_run_files "$RUNLIST") \
          || die "No valid ROOT files obtained from run list"
  else
      mapfile -t inputs < <(ls "$GROUP_DIR"/group_*.root 2>/dev/null | sort)
      (( ${#inputs[@]} )) || die "No group files – finish stage2 first or add skipStage2"
  fi

  COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"

  # --------------------------------------------------------------------------
  # 2. Local merge + QA
  # --------------------------------------------------------------------------
  if [[ "$mode" == "local" ]]; then
      note "Running final hadd locally"
      rm -f "$COMBINED"
      hadd -f "$COMBINED" "${inputs[@]}" || die "hadd failed"
      export COMBINED_ONLY=1 EXTERNAL_HADD=1
      "$EXEC_WRAPPER" --final
      [[ "$use_runlist" == "true" ]] || {
          rm -f "$GROUP_DIR"/group_*.root
          rmdir "$GROUP_DIR" 2>/dev/null || true
      }
      return
  fi

  # --------------------------------------------------------------------------
  # 3. Condor merge + QA
  # --------------------------------------------------------------------------
  local listfile="$TMP_BASE/hadd_input.lst"
  printf '%s\n' "${inputs[@]}" > "$listfile"

  local sub="$TMP_BASE/final_hadd.sub"
  cat > "$sub" <<EOF
universe      = vanilla
executable    = /bin/bash
transfer_input_files = $listfile
arguments     = -c "xargs -a $(basename "$listfile") hadd -f $COMBINED && \
                    export COMBINED_ONLY=1 EXTERNAL_HADD=1 && \
                    $EXEC_WRAPPER --final"
log           = $TMP_BASE/final_hadd.log
output        = $TMP_BASE/final_hadd.out
error         = $TMP_BASE/final_hadd.err
getenv        = True
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOF
  condor_submit "$sub" || die "condor_submit failed (final hadd)"
}

# ─────────────────────────────  main dispatch  ─────────────────────────────
case "${1:-}" in
  stage0)       shift; stage0 "$@" ;;
  stage1)       shift; stage1 "$@" ;;
  stage2)       shift; stage2 "$@" ;;
  stage3)       shift; stage3 "$@" ;;
  processOnly)
        # ────────────────────────────────────────────────────────────────
        # Usage examples
        #   ./runAuAuStages.sh processOnly                       → run all QA modules
        #   ./runAuAuStages.sh processOnly correlations          → correlations only
        #   ./runAuAuStages.sh processOnly correlations,pi0      → correlations + pi0
        #   ./runAuAuStages.sh processOnly local correlations    → same, executed locally
        #   ./runAuAuStages.sh processOnly condor correlations   → filtered QA, via Condor
        # ----------------------------------------------------------------
        shift                                # remove keyword ‘processOnly’
        execMode="local"                     # default = run here

        # optional first token  →  local | condor
        if [[ ${1:-} == local || ${1:-} == condor ]]; then
            execMode="$1"
            shift
        fi

        # optional second token  →  comma‑separated QA list
        # e.g.  correlations,pi0,jetqa
        if [[ $# -ge 1 && "$1" =~ ^[A-Za-z0-9_,]+$ ]]; then
            export QA_ONLY="${1,,}"          # lower‑case, pass to C++
            note "QA_ONLY filter applied → ${QA_ONLY}"
            shift
        fi

        COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"
        [[ -f "${COMBINED}" ]] || die "Combined file ${COMBINED} not found – run stage3 first"

        # Tell the C++ macro it should process only the combined file
        export COMBINED_ONLY=1
        export EXTERNAL_HADD=1

        if [[ "${execMode}" == "local" ]]; then
            "${EXEC_WRAPPER}" --final
        else
            sub="${TMP_BASE}/processOnly.sub"
            cat > "${sub}" <<EOF
universe      = vanilla
executable    = ${EXEC_WRAPPER}
arguments     = --final
output        = ${TMP_BASE}/processOnly.out
error         = ${TMP_BASE}/processOnly.err
log           = ${TMP_BASE}/processOnly.log
getenv        = True
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOF
            condor_submit "${sub}" || die "condor_submit failed (processOnly)"
        fi
        ;;
  *)
        echo "Usage: $0  stage0 | stage1 | stage2 | stage3 [local|condor] | processOnly"
        exit 1 ;;
esac
