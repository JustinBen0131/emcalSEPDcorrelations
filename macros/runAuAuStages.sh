#!/usr/bin/env bash
#################################################################################################
#  runAuAuStages.sh  –  Au+Au Run‑24/25 **end‑to‑end QA driver**
#
#  SYNOPSIS
#  --------
#    ./runAuAuStages.sh <command> [options]
#
#  COMMANDS (in‑house CLI)
#  -----------------------
#  • stage0 [rescueBusy]
#      Launch one Condor job **per run** (per‑run QA). With `rescueBusy`, kill
#      active Stage‑0 Condor jobs and process those runs locally, sequentially.
#
#  • stage1
#      Scan SEB/HCal; write the *good‑run list* and draw the
#      `Combined/MissingSEB_distribution.png` summary plot.
#
#  • stage2
#      Submit Condor “hadd” groups of ≤10 **good runs** each (uses the list from stage1).
#
#  • stage3 [local|condor] [skipStage2]
#      Final merge → `output_ALL_COMBINED.root` **and** run the combined QA.
#        - local  : do merge+QA on the login node.
#        - condor : submit the merge+QA to Condor (default).
#        - skipStage2 : build the combined file directly from the good‑run list
#                       produced by stage1 (i.e. no group files from stage2).
#
#  • processOnly [local|condor] [qaList] [triggerList]
#      Re‑run QA on the already‑combined file only.
#        - local|condor : where to execute (default = local).
#        - qaList       : comma‑separated subset of QA modules (see table below).
#                         If omitted → run the **full** QA suite.
#        - triggerList  : (optional) comma‑separated exact trigger names to limit
#                         processing (e.g. `MBD_NS_geq_2_vtx_lt_10,photon_8_plus_MBD_NS_geq_2_vtx_lt_10`).
#                         Internally forwarded as TRIGGER_ONLY. When active, TriggerQA
#                         is skipped (other modules still run for those triggers).
#
#  • processOnlyParrallel [ALL|qaList]
#      Submit one Condor job **per QA module** on the combined file.
#      Use `ALL` (default) or provide a comma‑separated `qaList`.
#      (Note: `TRIGGER_ONLY` is **not** applied in this parallel mode.)
#
#  QA MODULE KEYWORDS
#  ------------------
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
#  | vn             | vₙ analysis QA                     |
#
#  GLOBAL KNOBS (environment variables)
#  ------------------------------------
#  • VERBOSE=0|1|2      – 0: silent (default), 1: info logs, 2: +shell trace
#  • QA_ONLY=...        – comma‑separated QA modules (same tokens as the table)
#  • TRIGGER_ONLY=...   – comma‑separated trigger names to process (exact match)
#                         e.g. TRIGGER_ONLY="MBD_NS_geq_2_vtx_lt_10,photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
#TRIGGER_ONLY="MBD_NS_geq_2_vtx_lt_10,MBD_NS_geq_2_vtx_lt_30,MBD_NS_geq_2_vtx_lt_150" \
#VERBOSE=3 \
#./runAuAuStages.sh processOnlyParrallel correlations,hcal,mbd,sepd,sepdother,jetqa,eventqa,pi0,emcal,vn
#
#  QUICK EXAMPLES
#  --------------
#    ./runAuAuStages.sh stage0
#    ./runAuAuStages.sh stage0 rescueBusy
#    ./runAuAuStages.sh stage1
#    ./runAuAuStages.sh stage2
#    ./runAuAuStages.sh stage3                       # final QA via Condor
#    ./runAuAuStages.sh stage3 local skipStage2      # local merge from runlist + QA
#    ./runAuAuStages.sh processOnly                  # re‑run all QA on combined
#    ./runAuAuStages.sh processOnly correlations     # correlations only
#    ./runAuAuStages.sh processOnly condor pi0,jetqa # Condor; π0 + jet QA only
#    ./runAuAuStages.sh processOnly local emcal MBD_NS_geq_2_vtx_lt_10
#    ./runAuAuStages.sh processOnly condor correlations,pi0 \
#         MBD_NS_geq_2_vtx_lt_10,photon_8_plus_MBD_NS_geq_2_vtx_lt_10
####################################################################################################################
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
#include <TStyle.h>
#include <sstream>
#include <climits>
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

   /* ── Draw bar chart (fixed: key alignment, no stats box, dynamic title) ─── */
   const std::size_t nMissing = badRunMap.size();                      // runs with ≥1 missing SEB
   int runLow = INT_MAX, runHigh = -1;
   for (const auto &kv : badRunMap) {
        int r = std::stoi(kv.first);                                    // drop leading zeros
        if (r < runLow)  runLow  = r;
        if (r > runHigh) runHigh = r;
   }

   TH1I h("h","",16,-0.5,15.5);                                        // set title after filling
   for (int i=0;i<16;++i){
        char lab[8]; std::sprintf(lab,"SEB%02d",i);                     // pretty bin labels
        h.GetXaxis()->SetBinLabel(i+1, lab);

        // IMPORTANT: use non‑padded key to match how sebCnt[] is filled (e.g. "SEB0", ... "SEB15")
        std::string key = std::string("SEB") + std::to_string(i);
        int cnt = sebCnt.count(key) ? sebCnt.at(key) : 0;
        h.SetBinContent(i+1, cnt);
   }

   /* dynamic title: "Runs with Missing SEB (low → high, X runs with ≥ 1 missing)" */
   std::ostringstream ttl;
   if (nMissing > 0 && runLow != INT_MAX) {
        ttl << "Runs with Missing SEB (" << runLow
            << " #rightarrow " << runHigh
            << ", " << nMissing << " runs with #geq 1 missing)";
   } else {
        ttl << "Runs with Missing SEB (no runs with missing SEB)";
   }
   h.SetTitle( (ttl.str() + ";SEB index;Number of runs").c_str() );    // keep axis titles

   gStyle->SetOptStat(0);                                              // remove stats box
   h.SetFillColor(kAzure+1);
   h.SetBarWidth(0.8);
   h.SetBarOffset(0.1);

   TCanvas c("c","",800,500);
   c.SetGridy();
   h.Draw("bar2");

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
  stage0)
        note "[DISPATCH] stage0 | argv='${*}'"
        shift
        note "[DISPATCH] stage0 → calling: stage0 $*"
        stage0 "$@"
        ;;
  stage1)
        note "[DISPATCH] stage1 | argv='${*}'"
        shift
        note "[DISPATCH] stage1 → calling: stage1 $*"
        stage1 "$@"
        ;;
  stage2)
        note "[DISPATCH] stage2 | argv='${*}'"
        shift
        note "[DISPATCH] stage2 → calling: stage2 $*"
        stage2 "$@"
        ;;
  stage3)
        note "[DISPATCH] stage3 | argv='${*}'"
        shift
        note "[DISPATCH] stage3 → calling: stage3 $*"
        stage3 "$@"
        ;;
  quickEventCheck)
        # Usage:
        #   ./runAuAuStages.sh quickEventCheck <triggerName> [topN]
        # Example:
        #   ./runAuAuStages.sh quickEventCheck MBD_NS_geq_2_vtx_lt_10
        note "[DISPATCH] quickEventCheck | argv='${*}'"
        shift
        trig="${1:-}"
        topN="${2:-20}"
        note "[ARGS] trigger='${trig:-<empty>}'  topN='${topN}'"
        [[ -n "${trig}" ]] || die "Usage: $0 quickEventCheck <triggerName> [topN]"
        if ! [[ "${topN}" =~ ^[0-9]+$ ]]; then
          die "topN must be an integer (got '${topN}')"
        fi

        step "Quick Event Check – ${trig}"
        note "[SETUP] TMP_BASE='${TMP_BASE}'"
        note "[SETUP] INPUT_DIR='${INPUT_DIR}'"
        mkdir -p "${TMP_BASE}" || warn "[SETUP] mkdir -p '${TMP_BASE}' failed (continuing; may already exist)"

        # Emit a tiny ROOT macro that scans ${INPUT_DIR} and prints Top-N runs
        macro_path="${TMP_BASE}/quickEventCheck.C"
        note "[GEN] Writing ROOT macro to: ${macro_path}"
        cat > "${macro_path}" <<'CPP'
#include <filesystem>
#include <regex>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <cmath>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"

namespace fs = std::filesystem;

void quickEventCheck(const char* inDirC, const char* trigC, int topN = 20)
{
    const std::string inDir = inDirC ? inDirC : ".";
    const std::string trig  = trigC  ? trigC  : "";
    if (trig.empty()) {
        std::cerr << "[FATAL] trigger name is empty\n";
        return;
    }

    std::vector<std::pair<std::string, long long>> rows; // {run, scaledCount}
    const std::regex fileRe(R"(output_([0-9]{5,12})\.root)");

    for (const auto& e : fs::directory_iterator(inDir)) {
        if (!e.is_regular_file()) continue;
        const std::string fname = e.path().filename().string();

        std::smatch m;
        if (!std::regex_match(fname, m, fileRe)) continue;
        const std::string run = m[1].str();

        long long cnt = 0;
        std::unique_ptr<TFile> f(TFile::Open(e.path().c_str(), "READ"));
        if (f && !f->IsZombie()) {
            TDirectory* d = dynamic_cast<TDirectory*>(f->Get(trig.c_str()));
            if (d) {
                const std::string hname = "cnt_" + trig + "_scaled";
                TH1* h = dynamic_cast<TH1*>(d->Get(hname.c_str()));
                if (h) {
                    cnt = static_cast<long long>(std::llround(h->GetBinContent(1)));
                }
            }
        }
        rows.emplace_back(run, cnt);
    }

    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b){
                  if (a.second != b.second) return a.second > b.second; // desc by count
                  return a.first < b.first;                              // tie-break by run
              });

    const int n = std::min<int>(topN, rows.size());
    std::size_t runW = 3;
    for (const auto& r : rows) runW = std::max<std::size_t>(runW, r.first.size());

    std::cout << "\n==========  Quick Event Check  ==========\n";
    std::cout << "Directory : " << inDir << "\n";
    std::cout << "Trigger   : " << trig  << "\n";
    std::cout << "Histogram : cnt_" << trig << "_scaled\n\n";

    std::cout << std::left  << std::setw(6)      << "Rank"
              << std::setw(static_cast<int>(runW)+2) << "Run"
              << std::right << std::setw(16)     << "Scaled Counts" << "\n";

    std::cout << std::string(6 + (runW+2) + 16, '-') << "\n";
    for (int i = 0; i < n; ++i) {
        std::cout << std::left  << std::setw(6)                  << (i+1)
                  << std::setw(static_cast<int>(runW)+2)         << rows[i].first
                  << std::right << std::setw(16)                 << rows[i].second
                  << "\n";
    }
    std::cout << std::string(6 + (runW+2) + 16, '-') << "\n";
    std::cout << "Total runs scanned: " << rows.size() << "\n";
    std::cout << "=========================================\n";
}
CPP
        # Non-fatal confirmation of the generated file (do not change control flow)
        if [[ -s "${macro_path}" ]]; then
            note "[GEN] Macro created OK (size=$(wc -c < "${macro_path}") bytes)"
        else
            warn "[GEN] Macro file '${macro_path}' is missing or empty (ROOT will likely fail below)"
        fi

        note "[ROOT] Compiling and executing macro via ACLiC (C+)"
        note "[ROOT] Build dir: ${TMP_BASE}"
        note "[ROOT] Command: root -l -b -q -e \"gSystem->SetBuildDir(\\\"${TMP_BASE}\\\",kTRUE)\" \"${macro_path}+(\\\"${INPUT_DIR}\\\",\\\"${trig}\\\",${topN})\""

        # Compile + run the macro; print the table and exit
        root -l -b -q \
             -e "gSystem->SetBuildDir(\"${TMP_BASE}\",kTRUE)" \
             "${macro_path}+(\"${INPUT_DIR}\",\"${trig}\",${topN})" \
        || die "QuickEventCheck C++ macro failed"

        note "[DONE] quickEventCheck finished successfully"
        exit 0
        ;;
  processOnly)
        # ────────────────────────────────────────────────────────────────
        # Usage examples
        #   ./runAuAuStages.sh processOnly                       → run all QA modules
        #   ./runAuAuStages.sh processOnly correlations          → correlations only
        #   ./runAuAuStages.sh processOnly correlations,pi0      → correlations + pi0
        #   ./runAuAuStages.sh processOnly local correlations    → same, executed locally
        #   ./runAuAuStages.sh processOnly condor correlations   → filtered QA, via Condor
        # ----------------------------------------------------------------
        note "[DISPATCH] processOnly | argv='${*}'"
        shift                                # remove keyword ‘processOnly’
        execMode="local"                     # default = run here

        # optional first token  →  local | condor
        if [[ ${1:-} == local || ${1:-} == condor ]]; then
            execMode="$1"
            shift
        fi
        note "[ARGS] execMode='${execMode}'"

        # optional second token  →  comma‑separated QA list (lower‑cased)
        if [[ $# -ge 1 && "$1" =~ ^[A-Za-z0-9_,]+$ ]]; then
            export QA_ONLY="${1,,}"
            note "[FILTER] QA_ONLY='${QA_ONLY}'"
            shift
        else
            note "[FILTER] QA_ONLY (unset) – running full QA suite"
        fi

        # optional third token → trigger name(s), comma-separated
        if [[ $# -ge 1 && "$1" =~ ^[A-Za-z0-9_+,]+$ ]]; then
            export TRIGGER_ONLY="$1"
            note "[FILTER] TRIGGER_ONLY='${TRIGGER_ONLY}'"
            shift
        else
            note "[FILTER] TRIGGER_ONLY (unset) – all allowed triggers"
        fi

        COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"
        note "[CHECK] Looking for combined file: ${COMBINED}"
        [[ -f "${COMBINED}" ]] || die "Combined file ${COMBINED} not found – run stage3 first"

        # Tell the C++ macro it should process only the combined file
        export COMBINED_ONLY=1
        export EXTERNAL_HADD=1
        note "[ENV] COMBINED_ONLY=1  EXTERNAL_HADD=1  VERBOSE=${VERBOSE:-0}"

        if [[ "${execMode}" == "local" ]]; then
            note "[RUN] Local execution via: ${EXEC_WRAPPER} --final"
            "${EXEC_WRAPPER}" --final
        else
            note "[RUN] Submitting Condor job for processOnly"
            launcher="${TMP_BASE}/processOnly_launcher.sh"
            note "[FILE] Writing launcher: ${launcher}"
            cat > "${launcher}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export QA_ONLY="${QA_ONLY:-}"
export TRIGGER_ONLY="${TRIGGER_ONLY:-}"
export COMBINED_ONLY=1
export EXTERNAL_HADD=1
export VERBOSE="${VERBOSE:-0}"
exec "${EXEC_WRAPPER}" --final
EOF
            chmod +x "${launcher}"

            sub="${TMP_BASE}/processOnly.sub"
            note "[FILE] Writing submit file: ${sub}"
            cat > "${sub}" <<EOF
universe      = vanilla
executable    = ${launcher}
output        = ${TMP_BASE}/processOnly.out
error         = ${TMP_BASE}/processOnly.err
log           = ${TMP_BASE}/processOnly.log
stream_output = True
stream_error  = True
getenv        = True
request_memory= 4GB
+JobFlavour   = "tomorrow"
queue
EOF
            note "[CONDOR] condor_submit ${sub}"
            condor_submit "${sub}" || die "condor_submit failed (processOnly)"
        fi
        ;;
  processOnlyParrallel)
        # ────────────────────────────────────────────────────────────────
        # Parallel per‑module reprocessing on the already‑combined file.
        # Examples:
        #   ./runAuAuStages.sh processOnlyParrallel correlations,jetqa,pi0
        #   ./runAuAuStages.sh processOnlyParrallel ALL
        # ----------------------------------------------------------------
        note "[DISPATCH] processOnlyParrallel | argv='${*}'"
        shift  # remove keyword ‘processOnlyParrallel’

        COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"
        note "[CHECK] Looking for combined file: ${COMBINED}"
        [[ -f "${COMBINED}" ]] || die "Combined file ${COMBINED} not found – run stage3 first"

        # Full module list; 'ALL' expands to all of these
        all_modules=(correlations hcal mbd sepd sepdother jetqa eventqa triggerqa pi0 emcal vn)
        note "[INFO] Available modules: ${all_modules[*]}"

        # Parse requested modules
        req="${1:-ALL}"
        if [[ "${req^^}" == "ALL" ]]; then
            modules=( "${all_modules[@]}" )
            note "[ARGS] Requested modules: ALL → expanding to full set"
        else
            IFS=',' read -r -a modules <<< "${req,,}"
            note "[ARGS] Requested modules: ${modules[*]}"
        fi

        mkdir -p "${TMP_BASE}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}" || \
          warn "[SETUP] mkdir -p for tmp/log/stdout/stderr failed (continuing if already present)"

        # Submit one Condor job per module
        for m in "${modules[@]}"; do
            # validate module name
            if [[ ! " ${all_modules[*]} " =~ " ${m} " ]]; then
                warn "[FILTER] Unknown QA module '${m}' – skipped"
                continue
            fi

            note "[SUBMIT] Parallel processOnly job for module: ${m}"

            launcher="${TMP_BASE}/processOnly_${m}.sh"
            note "[FILE] Writing launcher: ${launcher}"
            cat > "${launcher}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export QA_ONLY="${m}"
export COMBINED_ONLY=1
export EXTERNAL_HADD=1
export VERBOSE="${VERBOSE:-0}"
exec "${EXEC_WRAPPER}" --final
EOF
            chmod +x "${launcher}"

            sub="${TMP_BASE}/processOnly_${m}.sub"
            note "[FILE] Writing submit file: ${sub}"
            cat > "${sub}" <<EOF
universe      = vanilla
executable    = ${launcher}
output        = ${STDOUT_DIR}/processOnly_${m}.out
error         = ${STDERR_DIR}/processOnly_${m}.err
log           = ${LOG_DIR}/processOnly_${m}.log
stream_output = True
stream_error  = True
getenv        = True
request_memory= 1.5GB
+JobFlavour   = "tomorrow"
queue
EOF
            note "[CONDOR] condor_submit ${sub}"
            condor_submit "${sub}" || warn "condor_submit failed for module ${m}"
        done
        ;;
  *)
        warn "[USAGE] Invalid or missing command: '${1:-<none>}'"
        echo "Usage: $0  stage0 [rescueBusy] | stage1 | stage2 | stage3 [local|condor] [skipStage2] | processOnly [local|condor] [qaList] [triggerList] | processOnlyParrallel [ALL|qaList] | quickEventCheck <triggerName> [topN]"
        exit 1
        ;;
esac
