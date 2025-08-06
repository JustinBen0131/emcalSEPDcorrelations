#!/usr/bin/env bash
##############################################################################
#  runAuAuStages.sh  –  4‑stage Au+Au Run‑24/25 QA pipeline
#
#    stage0  : submit run‑by‑run Condor grid  (per‑run QA)
#    stage1  : SEB diagnostics → good‑run list  + bar‑chart
#    stage2  : Condor group‑hadd jobs      (10 runs / job)
#    stage3  : final hadd  + combined QA   (local | condor)
##############################################################################
set -euo pipefail

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

# ──────────────────────────  stage‑0  –  per‑run Condor grid  ──────────────
stage0(){
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
request_memory= 4GB
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


# ────────────────────────────  stage‑2 helpers  ────────────────────────────
##
##  submit_group_hadd idx outfile infiles…
##  • builds a tiny wrapper shell script  (one per group)
##  • writes a HTCondor submit file that calls that wrapper
##  • submits the job and reports the cluster / proc id
##
submit_group_hadd() {
    local idx="$1"; shift                 # numeric group id 00,01,…
    local outFile="$1"; shift             # destination ROOT file
    local -a inFiles=( "$@" )             # all input ROOTs for this group

    printf "      ↳ crafting group %02d  (%d file%s)\n" \
           "${idx}" "${#inFiles[@]}" "$([[ ${#inFiles[@]} -eq 1 ]] && echo "" || echo "s")"

    # ------------------------------------------------------------------ #
    # 1.  Build the *wrapper* that Condor will actually execute
    # ------------------------------------------------------------------ #
    local wrapper="${TMP_BASE}/hadd_grp_${idx}.sh"
    cat > "${wrapper}" <<'EOS'
#!/bin/bash
set -euo pipefail
outFile="$1"; shift
echo "[hadd-grp$$] merging -> ${outFile}"
hadd -f "${outFile}" "$@"
EOS
    chmod +x "${wrapper}"

    # ------------------------------------------------------------------ #
    # 2.  Build the .sub file for Condor
    # ------------------------------------------------------------------ #
    local submit="${TMP_BASE}/group_${idx}.sub"
    cat > "${submit}" <<EOF
universe        = vanilla
executable      = ${wrapper}
arguments       = ${outFile} ${inFiles[*]}
log             = ${TMP_BASE}/group_${idx}.log
output          = ${TMP_BASE}/group_${idx}.out
error           = ${TMP_BASE}/group_${idx}.err
getenv          = True
request_memory  = 2GB
+JobFlavour     = "tomorrow"
queue
EOF

    # ------------------------------------------------------------------ #
    # 3.  Submit and report result
    # ------------------------------------------------------------------ #
    echo "      ↳ submitting group‑${idx} …"
    local submit_out
    if ! submit_out=$(condor_submit "${submit}" 2>&1); then
        die "condor_submit failed for group ${idx}\n${submit_out}"
    fi
    echo "         $(echo "${submit_out}" | grep -Eo '[0-9]+ job\(s\) submitted to cluster [0-9]+')"
}

# ───────────────────────────  stage‑2 main  ────────────────────────────────
stage2() {
    step "Stage‑2  –  split good runs into ≤10‑file batches and submit hadd jobs"

    # ── 0.  sanity checks ──────────────────────────────────────────────
    [[ -f "${RUNLIST}"     ]] || die "Run‑list ${RUNLIST} missing – run stage1 first"
    mapfile -t runs < "${RUNLIST}"
    (( ${#runs[@]} ))          || die "Run‑list is empty"
    mkdir -p "${GROUP_DIR}"

    echo "    total good runs            : ${#runs[@]}"
    echo "    files per Condor hadd job   : 10"
    echo "    target directory for groups : ${GROUP_DIR}"
    echo

    # ── 1.  walk the run list, chunk into groups of 10 ─────────────────
    local grpIdx=0
    local -a current=()
    for run in "${runs[@]}"; do
        local f="${INPUT_DIR}/output_${run}.root"
        if [[ ! -f "${f}" ]]; then
            warn "      ‼︎ ROOT file for run ${run} not found – skipped"
            continue
        fi

        current+=( "${f}" )
        printf "      • queued %-14s (grp %02d, slot %d/10)\n" \
               "$(basename "${f}")" "${grpIdx}" "${#current[@]}"

        if (( ${#current[@]} == 10 )); then
            echo "      → group complete – submit grp‑${grpIdx}"
            submit_group_hadd "${grpIdx}" \
                               "${GROUP_DIR}/group_${grpIdx}.root" \
                               "${current[@]}"
            current=()
            ((grpIdx++))
        fi
    done

    # ── 2.  final partial group (if any) ───────────────────────────────
    if (( ${#current[@]} )); then
        echo "      → final partial group (${#current[@]} file(s)) – submit grp‑${grpIdx}"
        submit_group_hadd "${grpIdx}" \
                           "${GROUP_DIR}/group_${grpIdx}.root" \
                           "${current[@]}"
    fi

    # ── 3.  cleanup & summary ──────────────────────────────────────────
    rm -f "${RUNLIST}"
    echo
    note "Stage‑2 finished – all group‑hadd jobs are now in the Condor queue."
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

stage3(){         # optional arg: local | condor (default condor)
  mode="${1:-condor}"
  step "Stage‑3  –  final hadd & QA (${mode})"
  mapfile -t groups < <(ls "${GROUP_DIR}"/group_*.root 2>/dev/null | sort)
  (( ${#groups[@]} )) || die "no group files – finish stage2 first"

  COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"
  if [[ "${mode}" == "local" ]]; then
      note "Running final hadd locally"
      rm -f "${COMBINED}"
      hadd -f "${COMBINED}" "${groups[@]}" || die "hadd failed"
      export COMBINED_ONLY=1 EXTERNAL_HADD=1
      "${EXEC_WRAPPER}" --final
      rm -f "${GROUP_DIR}"/group_*.root && rmdir "${GROUP_DIR}" 2>/dev/null || true
  else
      note "Submitting final hadd as Condor job"
      submit_final_hadd "${COMBINED}" "${groups[@]}"
  fi
}

# ─────────────────────────────  main dispatch  ─────────────────────────────
case "${1:-}" in
  stage0) shift; stage0 "$@" ;;
  stage1) shift; stage1 "$@" ;;
  stage2) shift; stage2 "$@" ;;
  stage3) shift; stage3 "$@" ;;
  *)  echo "Usage: $0  stage0 | stage1 | stage2 | stage3 [local|condor]"; exit 1 ;;
esac
