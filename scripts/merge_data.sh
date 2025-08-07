#!/usr/bin/env bash
###############################################################################
# merge_data.sh — sPHENIX Au+Au Run‑3 QA file‑merging helper
#
# OVERVIEW
# ▸ Each reconstruction “chunk” writes a ROOT file in
#       $CONDOR_OUT_BASE/<RUN>/<TAG>.root
# ▸ This script gathers those chunks and merges them with **hadd**, either
#   locally or on the HTCondor farm.  It also post‑processes the output by
#   re‑scaling trigger‑specific histograms with live/scaled counts from the
#   DAQ database.
#
# BASIC INVOCATION
#   ./merge_data.sh MODE [removeOtherJobs] [QUALIFIER]
#
#   MODE
#     condor             one Condor hadd job per run directory
#     condorStillRunning like condor, but chunks whose jobs are currently
#                        *IDLE* or *RUNNING* are **skipped** instead of killed
#     addRuns            merge the already‑produced per‑run ROOTs into a
#                        single “grand‑total” file
#     local <RUN>        merge a single run locally (quick debugging)
#
#   removeOtherJobs      (optional) purge *all* of your active HTCondor jobs
#                        before the merge; their chunk lists, partial ROOTs
#                        and the jobs themselves are deleted
#
#   QUALIFIER
#     condor             (only for addRuns) run the grand‑total merge on Condor
#     test               (only for condor / condorStillRunning) submit just the
#                        first idle run, useful for smoke tests
#     firstHalf          process only the first 50 % of idle runs
#     <runNumber>        required for MODE = local
#
# QUICK EXAMPLES
#   # Submit every idle run to Condor (one job per run)
#   ./merge_data.sh condor
#
#   # Resume a partially finished campaign; keep busy jobs alive
#   ./merge_data.sh condorStillRunning
#
#   # Same as above, but with shell tracing and verbose logging
#   DEBUG=1 ./merge_data.sh condorStillRunning
#
#   # Regenerate run 66484 locally, ignoring busy runs, with full debug output
#   DEBUG=1 ./merge_data.sh condorStillRunning local 00066484
#
#   # Produce the final grand‑total file on Condor
#   ./merge_data.sh addRuns condor
#
#   • Busy‑run detection
#   • **removeOtherJobs** switch:
#       – enumerates every live HTCondor job belonging to \$USER
#       – deletes its chunk list, partial ROOT outputs and the job itself
#       – guarantees the merge sees only complete, uncorrupted input.
#
# ENVIRONMENT
#   DEBUG=1          enables `set -x` and extra logging
#   PGHOST           overrides the host used for DAQ‑DB access
#
# NOTE
#   Individual chunk ROOTs are **never** deleted by default.  Only the merged
#   per‑run files and Condor log/output files are reconstructed on each launch.
###############################################################################

###############################################################################
# ---- 1. Strict mode + debug plumbing ----------------------------------------
###############################################################################
set -euo pipefail
IFS=$'\n\t'
DEBUG=${DEBUG:-0}        # ensure DEBUG always exists (0 if unset)
(( DEBUG )) && set -x

trap 'err=$?; printf "\033[0;31m%s  ✘  line %d – cmd `%s` exited %d\033[0m\n" \
                "$(date "+%F %T")" "${BASH_LINENO[0]}" "$BASH_COMMAND" "$err" >&2' ERR

###############################################################################
# ---- 2. Colour helpers + logging wrappers -----------------------------------
###############################################################################
CLR_B="\033[1;34m"; CLR_G="\033[0;32m"; CLR_Y="\033[1;33m"; CLR_R="\033[0;31m"; CLR_RST="\033[0m"
ts()   { date "+%F %T"; }
say()  { printf "${CLR_B}%s  ➜  %s${CLR_RST}\n"  "$(ts)" "$*"; }
good() { printf "${CLR_G}%s  ✔  %s${CLR_RST}\n"  "$(ts)" "$*"; }
warn() { printf "${CLR_Y}%s  ⚠  %s${CLR_RST}\n"  "$(ts)" "$*" >&2; }
fatal(){ printf "${CLR_R}%s  ✘  %s${CLR_RST}\n"  "$(ts)" "$*" >&2; exit 1; }

###############################################################################
# ---- 3. User settings -------------------------------------------------------
###############################################################################
CONDOR_OUT_BASE="/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations"
OUTPUT_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/output"
RUN_MERGED_PREFIX="output"          # → output_<run>.root
TMP_LIST_DIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/tmp_run_lists"
CONDOR_STDOUT="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/stdout"
CONDOR_STDERR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/error"
CONDOR_LOGDIR="/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/log"
QA_CMD_FILTER='regexp("run_auau_run3_qa.sh",Cmd)'
HADD_WRAPPER="hadd_run_condor.sh"
REQUEST_MEMORY="2000MB"
JOB_PRIO=100000


mkdir -p "$OUTPUT_DIR" "$TMP_LIST_DIR" "$CONDOR_STDOUT" \
         "$CONDOR_STDERR" "$CONDOR_LOGDIR"

###############################################################################
# ---- 4. Usage & CLI parsing -------------------------------------------------
###############################################################################
usage() {
  cat <<EOF
Usage:
  $0 MODE [removeOtherJobs] [QUALIFIER]

  MODE
    condor   – merge every idle run via Condor (per‑run hadd)
               QUALIFIER = test | firstHalf
    addRuns  – grand‑total merge of per‑run outputs
               QUALIFIER = condor   (do the grand‑total on Condor)
    local    – merge a single run locally
               QUALIFIER = <runNumber>

  removeOtherJobs
      (optional, may follow any MODE)
      Purge every still‑running job in your Condor queue **before**
      the merge starts.  The purge:
        • deletes the affected *.list chunk files
        • deletes partial *.root files in \$CONDOR_OUT_BASE/<run>
        • runs 'condor_rm \$USER'

Environment:
  DEBUG=1   enable shell trace & extra logging
EOF
  exit 1
}

(( $# >= 1 )) || usage
MODE=$1; shift

#  -- optional purge switch ---------------------------------------------------
PURGE=0
if [[ ${1:-} == removeOtherJobs ]]; then
    PURGE=1
    shift
fi

SUBMODE=${1:-}
RUNNUM=${2:-}          # only meaningful when SUBMODE == local

case "$MODE" in
  condor)
        [[ -z "$SUBMODE" || "$SUBMODE" =~ ^(test|firstHalf)$ ]] || usage
        ;;

  condorStillRunning)
        if [[ -z "$SUBMODE" || "$SUBMODE" =~ ^(test|firstHalf)$ ]]; then
            :
        elif [[ "$SUBMODE" == local && "$RUNNUM" =~ ^[0-9]{5,8}$ ]]; then
            :
        else
            usage
        fi
        ;;

  rescueBusy)          # NEW ─ kill active jobs and finish those runs locally
        [[ -z "$SUBMODE" ]] || usage
        ;;

  addRuns)
        [[ -z "$SUBMODE" || "$SUBMODE" == condor ]] || usage
        ;;

  local)
        [[ -n "$SUBMODE" && "$SUBMODE" =~ ^[0-9]{5,8}$ ]] || usage
        ;;

  *)
        usage
        ;;
esac


# Flag: should we keep busy runs and merely skip the still‑running chunk files?
SKIP_RUNNING=0
[[ $MODE == condorStillRunning ]] && SKIP_RUNNING=1

###############################################################################
# ---- 5. Build the tiny wrapper executed inside each Condor slot ------------
###############################################################################
cat > "$HADD_WRAPPER" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
[[ ${DEBUG_WRAPPER:-0} -ne 0 ]] && set -x         # export DEBUG_WRAPPER=1 for x‑trace
LIST=$1; OUT=$2
[[ -s $LIST ]] || { echo "[FATAL] empty list $LIST"; exit 2; }

#  sPHENIX environment --------------------------------------------------------
set +u
export PGHOST=${PGHOST:-localhost}
source /opt/sphenix/core/bin/sphenix_setup.sh -n
set -u

exec 1> >(stdbuf -oL cat) 2>&1          # live stdout/err streaming
echo -e "\e[1;34m[wrapper] $(wc -l <"$LIST") inputs  →  $OUT\e[0m"
echo "[wrapper] Input file list:"
nl -ba "$LIST"
hadd -v -v -v -f "$OUT" @"$LIST" || { echo -e "\e[0;31m[FATAL] hadd failed – abort\e[0m"; exit 3; }

# ──────────────────────────────────────────────────────────────────────
#  1.  determine run‑number from final file name (…_<RUN>.root)
#  2.  query DAQ DB for live / scaled counts (per trigger bit)
#  3.  rescale every histogram in the matching trigger directory
#      → prints OLD integral, FACTOR, NEW integral
# ──────────────────────────────────────────────────────────────────────
runNum=${OUT##*_}; runNum=${runNum%.root}
echo -e "\e[1;34m[wrapper] commencing trigger‑prescale scaling for run $runNum\e[0m"

export RUNNUM_ENV=$runNum        # pass to ROOT
export OUTFILE_ENV=$OUT

set -o pipefail                               # keep PIPESTATUS working
ROOT_LOG=$(mktemp /tmp/scale_${runNum}_XXXX.log)

root -l -b <<'EOF' 2>&1 | tee "$ROOT_LOG"
{
  //--------------------------------------------------------------------
  //                   verbose‑scaling.C  (inline)
  //--------------------------------------------------------------------
  #include <iostream>
  #include <iomanip>
  #include <map>
  #include <regex>
  #include <memory>
  using std::cout; using std::cerr; using std::endl;

  gErrorIgnoreLevel = kInfo;                 // print all ROOT messages
  if (std::getenv("DEBUG_ROOT")) gDebug = 1; // export DEBUG_ROOT=1 for ROOT‑level trace

  /* 0. ANSI helpers --------------------------------------------------- */
  const char* GRN="\033[0;32m", *YEL="\033[1;33m", *RED="\033[0;31m", *BLU="\033[1;34m", *RST="\033[0m";

  /* 1. trigger‑bit ↔ folder mapping ---------------------------------- */
  std::map<int,std::string> trigName = {
      {10, "MBD_NS_geq_2"},
      {11, "MBD_NS_geq_1"},
      {12, "MBD_NS_geq_2_vtx_lt_10"},
      {13, "MBD_NS_geq_2_vtx_lt_30"},
      {14, "MBD_NS_geq_2_vtx_lt_150"},
      {15, "MBD_NS_geq_1_vtx_lt_10"},
      {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
      {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
      {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
      {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
      {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
      {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
      {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
      {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
  };

  /* 2. DB query – collect live / scaled ------------------------------- */
  const char* envRun  = gSystem->Getenv("RUNNUM_ENV");                 // may be nullptr
  int         run     = -1;

  /*––– determine run number (strip leading 0s) ––––––––––––––––––––––––*/
  if (envRun && *envRun) {
      run = std::strtol(envRun, nullptr, 10);                          // 00066484 → 66484
  } else {
      std::string fn = gSystem->Getenv("OUTFILE_ENV") ? gSystem->Getenv("OUTFILE_ENV") : "";
      std::smatch  m;
      if (std::regex_search(fn, m, std::regex("([0-9]{5,8})\\.root$")))
          run = std::stoi(m[1]);
  }
  if (run <= 0) {
      cerr << RED << "[ERROR] cannot determine run number – abort" << RST << endl;
      return;
  }

  /*––– very early diagnostics –––––––––––––––––––––––––––––––––––––––––*/
  cout << BLU << "[debug] RUNNUM_ENV="   << (envRun ? envRun : "<unset>")
       << "  OUTFILE_ENV="              << (gSystem->Getenv("OUTFILE_ENV") ? gSystem->Getenv("OUTFILE_ENV") : "<unset>")
       << RST << endl;

  std::map<std::string,double> scaleFac;
  cout << BLU << "[scaling] querying DAQ DB for run " << run << RST << endl;

  /*––– open DB connection –––––––––––––––––––––––––––––––––––––––––––––*/
  std::unique_ptr<TSQLServer> db(
      TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro",""));

  if (!db || db->IsZombie()) {
      cerr << RED << "[ERROR] DB connect failed" << RST << endl;
      cerr << YEL << "[debug] PGHOST=" << (gSystem->Getenv("PGHOST") ? gSystem->Getenv("PGHOST") : "<unset>")
           << "  PGPORT=5432  DB=daq  USER=phnxro" << RST << endl;
      return;
  }
  cout << GRN << "[debug] DB connection OK – " << db->ServerInfo() << RST << endl;

  /*──– sanity: show first few tables (helps if schema changes) ––––––––*/
  {
      std::unique_ptr<TSQLResult> tbl(db->GetTables("daq", "%"));     // ← fixed
      cout << BLU << "[debug] available tables (first 10):" << RST << endl;
      for (int i = 0; tbl && i < 10; ++i) {
            std::unique_ptr<TSQLRow> row(tbl->Next());
            if (!row) break;
            cout << "   • " << row->GetField(0) << endl;
       }
  }

  /*––– main scaler query ––––––––––––––––––––––––––––––––––––––––––––––*/
  char q[512];
  std::sprintf(q,
      "SELECT s.index, s.live, s.scaled "
      "FROM gl1_scalers s WHERE s.runnumber=%d ORDER BY s.index;", run);
  cout << YEL << "[query] " << q << RST << endl;

  std::unique_ptr<TSQLResult> res(db->Query(q));
  if (!res || res->GetRowCount()==0) {
      cerr << RED << "[ERROR] query returned no rows" << RST << endl;

      /* extra diagnostics – show how many rows exist for that run */
      char qc[256];
      std::sprintf(qc,
        "SELECT COUNT(*) FROM gl1_scalers WHERE runnumber=%d;", run);

      std::unique_ptr<TSQLResult> rc(db->Query(qc));
      if (rc) {
          std::unique_ptr<TSQLRow> row(rc->Next());
          if (row)
              cout << YEL << "[debug] gl1_scalers rows for run " << run
                 << " = " << row->GetField(0) << RST << endl;
      }

      return;
  }

  /*––– consume scaler rows –‑ LOAD ALL SCALER DATA FIRST –––––––––––––‑*/
  std::map<int, std::pair<double,double>> scalerMap;   // idx → {live,scaled}
  while (auto row = res->Next()) {
      int idx = std::atoi(row->GetField(0));
      scalerMap[idx] = { std::atof(row->GetField(1)), std::atof(row->GetField(2)) };
      delete row;
  }

  /*––– now iterate over **every** trigger in trigName –––––––––––––––––*/
  for (const auto &tg : trigName)
  {
      int idx = tg.first;
      const std::string &dir = tg.second;

      auto itS = scalerMap.find(idx);
      if (itS == scalerMap.end()) {
          cout << YEL << "[skip] " << dir << " (no row in gl1_scalers)" << RST << endl;
          continue;
      }
      double live   = itS->second.first;
      double scaled = itS->second.second;

      /*── scaledown check – -1 means trigger OFF ––––––––––––––––––––––*/
      char qsd[256];
      std::sprintf(qsd,
          "SELECT scaledown%d FROM gl1_scaledown WHERE runnumber=%d;", idx, run);
      std::unique_ptr<TSQLResult> sdres(db->Query(qsd));
      double sdFactor = -1;
      if (sdres) {
          std::unique_ptr<TSQLRow> sdrow(sdres->Next());
          if (sdrow) sdFactor = std::atof(sdrow->GetField(0));
      }
      if (sdFactor < 0) {
          cout << YEL << "[skip] " << dir
               << " (scaledown = -1, no scaling applied)" << RST << endl;
          continue;
      }

      /*── active trigger – calculate scale factor –––––––––––––––––––––*/
      double factor = (scaled > 0) ? live / scaled : -1;
      scaleFac[dir] = factor;

      cout << std::setw(4) << idx << " → "
           << std::setw(40) << dir
           << " | live="   << std::setw(10) << live
           << " scaled="   << std::setw(10) << scaled
           << "  ⇒  factor=" << factor << endl;
  }

  if (scaleFac.empty()) {
      cerr << RED << "[ERROR] no scale factors found – abort" << RST << endl;
      return;
  }

  /* 3. open ROOT file ------------------------------------------------- */
  const char* fName = gSystem->Getenv("OUTFILE_ENV");
  TFile f(fName,"UPDATE");
  if (f.IsZombie()) { cerr << RED << "[ERROR] cannot open " << fName << RST << endl; return; }

  size_t nDirScaled=0, nHistScaled=0;

  for (auto& kv : scaleFac)
  {
      const std::string& dir = kv.first;
      double fac = kv.second;
      if (fac<=0) { cout << YEL << "[skip] " << dir << " (factor <=0)" << RST << endl; continue; }

      TDirectory* d = f.GetDirectory(dir.c_str());
      if (!d) { cout << YEL << "[skip] directory " << dir << " not present" << RST << endl; continue; }

      cout << GRN << "[scale] " << dir << "  factor=" << fac << RST << endl;
      ++nDirScaled;

      /* --- snapshot the list of keys --------------------------------- */
      std::vector<std::string> keyNames;
      {
          TIter itKey(d->GetListOfKeys());
          while (auto *k = static_cast<TKey*>(itKey()))
              keyNames.emplace_back(k->GetName());
      }

      /* --- now scale each histogram exactly once --------------------- */
      for (const auto &hname : keyNames)
      {
          TObject *obj = d->Get(hname.c_str());
          if (!obj || !obj->InheritsFrom("TH1")) { delete obj; continue; }
          TH1 *h = static_cast<TH1*>(obj);

          /*  Histograms to exclude from scaling
           *  – all counters   :  cnt_*              (unchanged)
           *  – master trigger :  h_MB_vs_Trigger    (unchanged)
           *  – any spectrum whose name **contains** “doNotScale_”
           *    e.g.  h_maxClusterEnergy_doNotScale_MBD_NS_geq_2_vtx_lt_150
           *           ↑───────────────────────────── added test
           */
          if (hname.rfind("cnt_",0)==0 ||
              hname=="h_MB_vs_Trigger" ||
              hname.find("doNotScale_")!=std::string::npos)
          {   delete h;  continue;   }

          double oldInt = h->Integral();
          if (oldInt==0) { delete h; continue; }

          h->Scale(fac);
          double newInt = h->Integral();

          cout << "   • " << std::left << std::setw(45) << hname
               << "  " << std::right << std::fixed << std::setprecision(1)
               << oldInt << "  →  " << newInt << endl;

          d->cd();
          h->Write(hname.c_str(), TObject::kOverwrite);
          delete h;
          ++nHistScaled;
      }
      
      f.cd();
  }

  f.Write(); f.Close();
  cout << BLU << "[summary] scaled " << nHistScaled << " histograms in "
       << nDirScaled << " trigger directories" << RST << endl;
}
.q
EOF
ROOT_RC=${PIPESTATUS[0]}                     # true exit code of ROOT
echo -e "\e[1;34m[wrapper] ROOT exit code = ${ROOT_RC}  •  full log ⇒ ${ROOT_LOG}\e[0m"

if (( ROOT_RC != 0 )); then
    echo -e "\e[0;31m[FATAL] ROOT scaling macro aborted – inspect ${ROOT_LOG}\e[0m"
    exit 4
fi

echo -e "\e[0;32m[wrapper] scaling finished – final file size: $(du -h "$OUT" | cut -f1)\e[0m"

EOS
chmod +x "$HADD_WRAPPER"

###############################################################################
# ---- 6. Helpers -------------------------------------------------------------
###############################################################################

# When SKIP_RUNNING = 1 build once‑per‑script a list of files still produced
SKIP_FILE="$TMP_LIST_DIR/skip_running_jobs.txt"
if (( SKIP_RUNNING )); then
    say "Building skip‑list for active Condor chunks → $(basename "$SKIP_FILE")"
    : > "$SKIP_FILE"
    condor_q "$USER" \
        -constraint "$QA_CMD_FILTER && (JobStatus == 1 || JobStatus == 2)" \
        -af Args 2>/dev/null |
      awk -v base="$CONDOR_OUT_BASE" '
          {
              run=$1; tag=$3;
              if (run ~ /^[0-9]+$/ && tag != "")
                  printf "%s/%s/%s.root\n", base, run, tag
          }' | sort -u > "$SKIP_FILE"
    lines=$(wc -l < "$SKIP_FILE")
    say "  • $lines file(s) will be ignored during hadd"
fi

safe_find() {                                  # robust, fault‑tolerant “find”
  local dir=$1 list=$2
  say "  • scanning $dir"

  ###########################################################################
  # 1. Skip‑list is NOT refreshed here.
  #    We keep the list that was built once at script start so that any
  #    chunk that was IDLE or RUNNING back then is excluded from every merge,
  #    even if the job finishes while we are scanning directories.
  ###########################################################################
  if (( SKIP_RUNNING )); then
      :
  fi

  ###########################################################################
  # 2. INITIAL CANDIDATE LIST  – all *.root presently on disk
  ###########################################################################
  (
    set +e +o pipefail
    find "$dir" -type f -name '*.root' -print 2> >(while read -r l; do warn "    find: $l"; done) |
      sort
  ) > "${list}.00_all"

  ###########################################################################
  # 3. SKIP LIST  – remove files whose chunks are still in Condor
  ###########################################################################
  if (( SKIP_RUNNING )) && [[ -s "$SKIP_FILE" ]]; then
      grep -F -v -f "$SKIP_FILE" "${list}.00_all" > "${list}.01_skip"
  else
      mv "${list}.00_all" "${list}.01_skip"
  fi

  ###########################################################################
  # 4. OPEN‑FILE FILTER  – drop ROOTs that are *currently* open for writing
  #    (covers the extremely slow‑write edge‑case).
  ###########################################################################
  if command -v lsof >/dev/null 2>&1; then
      lsof +D "$dir" 2>/dev/null | awk '/\.root$/ {print $9}' | sort -u > "${list}.open"
      if [[ -s "${list}.open" ]]; then
          grep -F -v -f "${list}.open" "${list}.01_skip" > "${list}"
      else
          mv "${list}.01_skip" "${list}"
      fi
      rm -f "${list}.open"
  else
      mv "${list}.01_skip" "${list}"
  fi
  rm -f "${list}.00_all" "${list}.01_skip"

  ###########################################################################
  # 5. FINAL SANITY  – abort run if nothing is left to merge
  ###########################################################################
  if [[ ! -s "$list" ]]; then
      warn "    ➜ no eligible ROOT files after safety filters"
      return 1
  fi

  (( DEBUG )) && { say "    first 10 entries:"; head -n 10 "$list" | sed 's/^/      /'; }
}



# ---- 6a. busy‑run cache -----------------------------------------------------
declare -Ag busySet=()                 # busySet[run]=1   (global)

refresh_busy_runs() {
  busySet=()

  # ---- A.  still‑running *chunk* jobs (run_auau_run3_qa.sh) ----------
  while read -r run; do
      [[ $run =~ ^[0-9]{5,8}$ ]] && busySet["$run"]=1
  done < <(
      condor_q "$USER" \
          -constraint 'regexp("run_auau_run3_qa.sh",Cmd) && (JobStatus == 1 || JobStatus == 2)' \
          -af Args 2>/dev/null |
      awk '{print $1}'
  )

  # ---- B.  still‑running *per‑run‑hadd* jobs (hadd_run_condor.sh) ----
  # Each Args line looks like:
  #   /path/to/in_00066749.txt  /path/to/output_00066749.root
  # We need just the 00066749 part; also prevent set ‑e from aborting
  # when the final read returns 1 (EOF).
  while read -r line || [[ -n $line ]]; do
      [[ $line =~ in_([0-9]{5,8})\.txt ]] && busySet["${BASH_REMATCH[1]}"]=1
  done < <(
      condor_q "$USER" \
          -constraint 'regexp("hadd_run_condor.sh",Cmd) && (JobStatus == 1 || JobStatus == 2)' \
          -af Args 2>/dev/null
  ) || true
}

# ---- 6b. heavy‑duty purge ---------------------------------------------------
purge_busy_jobs() {
  refresh_busy_runs
  (( ${#busySet[@]} )) || { good "No running Condor jobs – nothing to purge"; return; }

  say  "removeOtherJobs  –  purging $(printf '%d' "${#busySet[@]}") busy run(s)"

  # 1) remove chunk *.list files ------------------------------------------------
  mapfile -t busyLists < <(
      condor_q "$USER" -constraint "$QA_CMD_FILTER" -af Args 2>/dev/null |
      grep -Eo '/[^[:space:]]+tmp_condor_lists[^[:space:]]+\.list' | sort -u
  )
  if (( ${#busyLists[@]} )); then
      say "  • deleting ${#busyLists[@]} temporary list(s)"
      for f in "${busyLists[@]}"; do [[ -f $f ]] && rm -f "$f" && say "      – $(basename "$f")"; done
  else
      say "  • no chunk lists to delete"
  fi

  # 2) remove partial output ROOTs --------------------------------------------
  say "  • cleaning per‑run output directories"
  for run in "${!busySet[@]}"; do
      # build the list of ROOT files that the *still‑running* job(s) would create
      mapfile -t toDelete < <(
          condor_q "$USER" -constraint "$QA_CMD_FILTER" -af Args 2>/dev/null |
          awk -v r="$run" -v base="$CONDOR_OUT_BASE" '$1==r {printf "%s/%s/%s.root\n", base, $1, $3}'
      )
      (( ${#toDelete[@]} )) || { say "      – run $run  (no matching partial outputs)"; continue; }

      for f in "${toDelete[@]}"; do
          [[ -f $f ]] && rm -f "$f"
      done
      say "      – run $run  (${#toDelete[@]} file(s) purged)"
  done

  # 3) remove Condor jobs ------------------------------------------------------
  say "  • condor_rm $USER"
  if condor_rm "$USER" >/dev/null 2>&1; then
      good "All Condor jobs removed"
  else
      warn "condor_rm failed – please check manually"
  fi
}

###############################################################################
# ---- 7. Optional pre‑merge purge -------------------------------------------
###############################################################################
# ---- 7. Optional pre‑merge purge -------------------------------------------
if (( PURGE )); then          # -- purge requested
    purge_busy_jobs           #    ① delete tmp files + condor_rm
    busySet=()                #    ② pretend no runs are “busy” any more
else                          # -- normal path
    refresh_busy_runs         #    rebuild busySet from condor_q
fi

###############################################################################
# ---- 8.  PER‑RUN MERGE (MODE = condor) --------------------------------------
###############################################################################
if [[ ($MODE == condor || $MODE == condorStillRunning) && $SUBMODE != local ]]; then
  ###########################################################################
  # (Re)initialise work‑area for *every* submission, even when we are in
  # condorStillRunning mode.  We remove only the per‑run merged ROOT files
  # and the old Condor log/err/out to guarantee a clean slate; individual
  # segment ROOTs inside $CONDOR_OUT_BASE are NEVER touched.
  ###########################################################################
  say "Cleaning previous Condor text outputs"
  for d in "$CONDOR_STDOUT" "$CONDOR_STDERR" "$CONDOR_LOGDIR"; do
      [[ -d $d ]] && find "$d" -type f -delete
  done

  say "Removing stale per‑run ROOT files from $OUTPUT_DIR (will be regenerated)"
  find "$OUTPUT_DIR" -maxdepth 1 -type f -name "${RUN_MERGED_PREFIX}_????????.root" -delete


  # Step 1 – enumerate run directories ----------------------------------------
  say "Step 1 – enumerating run directories under $CONDOR_OUT_BASE"
  mapfile -t runs < <(find "$CONDOR_OUT_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

  # Step 2 – busy‑run handling -------------------------------------------------
  if (( SKIP_RUNNING )); then
      say "Step 2 – keeping runs that still have active chunks (individual files will be skipped)"
      # nothing is removed from the runs[] array
  else
      say "Step 2 – filtering out active Condor runs"
      runs=( $(for r in "${runs[@]}"; do [[ -z ${busySet[$r]+x} ]] && echo "$r"; done) )
      (( ${#busySet[@]} )) && say "    active: $(printf '%s ' "${!busySet[@]}")" || say "    none"
      (( ${#runs[@]} )) || { good "Nothing idle to merge"; exit 0; }
  fi

  say "Step 3 – ${#runs[@]} idle run(s) will be processed"
  case "$SUBMODE" in
      test)      runs=( "${runs[0]}" ); warn "TEST mode – only ${runs[0]}" ;;
      firstHalf) half=$(( ${#runs[@]} / 2 )); runs=( "${runs[@]:0:$half}" ); warn "FIRST‑HALF mode – $half run(s)" ;;
  esac

  # Condor submit description --------------------------------------------------
  SUB="$TMP_LIST_DIR/merge_runs.sub"; : > "$SUB"
  cat >>"$SUB" <<EOT
universe        = vanilla
executable      = $HADD_WRAPPER
output          = $CONDOR_STDOUT/merge.\$(Cluster).\$(Process).out
error           = $CONDOR_STDERR/merge.\$(Cluster).\$(Process).err
log             = $CONDOR_LOGDIR/merge.\$(Cluster).\$(Process).log
request_memory  = $REQUEST_MEMORY
priority        = $JOB_PRIO
getenv          = True
stream_output   = True
stream_error    = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
EOT

  jobCnt=0
  for run in "${runs[@]}"; do
      say "── Run $run ───────────────────────────────────────────────────────"
      inDir=$CONDOR_OUT_BASE/$run
      list=$TMP_LIST_DIR/in_${run}.txt
      outFile=$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root
      safe_find "$inDir" "$list" || continue
      nFiles=$(wc -l <"$list")
      say "    will merge $nFiles file(s) → $outFile"
      printf 'arguments = %s %s\nqueue 1\n' "$list" "$outFile" >> "$SUB"
      (( ++jobCnt ))
  done

  (( jobCnt )) || fatal "Zero jobs created – aborting"
  (( DEBUG )) && { say "Submit description (first 10 lines):"; head -n 10 "$SUB" | sed 's/^/    /'; }
  say "Step 4 – submitting $jobCnt job(s) to Condor"
  condor_submit "$SUB"
  good "Condor submission finished"
  exit 0
fi


###############################################################################
# ---- 9a.  RESCUE‑BUSY RUNS LOCALLY  (MODE = rescueBusy)  ---------------------
###############################################################################
if [[ $MODE == rescueBusy ]]; then
  refresh_busy_runs
  (( ${#busySet[@]} )) || { good "No active Condor jobs – nothing to rescue"; exit 0; }

  say "Rescuing ${#busySet[@]} run(s) still in Condor:  ${!busySet[@]}"
  # kill both chunk‑analysis and per‑run‑hadd jobs in one go
  busyExpr='regexp("run_auau_run3_qa.sh|hadd_run_condor.sh",Cmd) && (JobStatus == 1 || JobStatus == 2)'
  nKill=$(condor_q "$USER" -constraint "$busyExpr" -af ClusterId ProcId 2>/dev/null | wc -l)
  if (( nKill > 0 )); then
      say "  • removing $nKill active Condor job(s)"
      condor_rm "$USER" -constraint "$busyExpr" \
          || warn "condor_rm returned non‑zero – continuing anyway"
  else
      say "  • no active Condor jobs to remove"
  fi

  for run in "${!busySet[@]}"; do
      say "Local HADD for run ${run}"
      inDir=$CONDOR_OUT_BASE/$run
      list=$TMP_LIST_DIR/in_${run}.txt
      outFile=$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root

      safe_find "$inDir" "$list" || { warn "    no finished chunks yet – skipped"; continue; }
      [[ -f $outFile ]] && { say "    removing stale $outFile"; rm -f "$outFile"; }

      "$HADD_WRAPPER" "$list" "$outFile"
      good "Run ${run} merged locally  →  $(du -h "$outFile" | cut -f1)"
  done

  good "Rescue phase finished"
  exit 0
fi


###############################################################################
# ---- 9b.  SINGLE‑RUN LOCAL MERGE --------------------------------------------
###############################################################################
if [[ $MODE == local || ( $MODE == condorStillRunning && $SUBMODE == local ) ]]; then
  if [[ $MODE == local ]]; then
      run="$SUBMODE"         # ./merge_data.sh  local  <run>
  else
      run="$RUNNUM"          # ./merge_data.sh  condorStillRunning  local  <run>
  fi
  
  [[ -n ${busySet[$run]+x} ]] && { warn "Run $run is still active in Condor – skipping local merge"; exit 0; }

  say "Local merge for run ${run}"
  inDir=$CONDOR_OUT_BASE/$run
  list=$TMP_LIST_DIR/in_${run}.txt
  outFile=$OUTPUT_DIR/${RUN_MERGED_PREFIX}_${run}.root
  safe_find "$inDir" "$list" || fatal "No ROOT files found for run ${run}"
  [[ -f $outFile ]] && { say "  • removing old $outFile"; rm -f "$outFile"; }
  "$HADD_WRAPPER" "$list" "$outFile"
  good "Merge finished – $(ls -lh "$outFile")"
  exit 0
fi

###############################################################################
# ---- 10. GRAND‑TOTAL MERGE --------------------------------------------------
###############################################################################
say "Grand‑total stage – collecting per‑run outputs in $OUTPUT_DIR"
mapfile -t runFiles < <(
  find "$OUTPUT_DIR" -maxdepth 1 -type f -name "${RUN_MERGED_PREFIX}_*.root" |
  while read -r f; do
      [[ $f =~ _([0-9]{5,8})\.root$ ]] || continue
      r=${BASH_REMATCH[1]}
      [[ -z ${busySet[$r]+x} ]] && echo "$f"
  done | sort
)
(( ${#runFiles[@]} < 2 )) && { good "Nothing to add"; exit 0; }

LIST_ALL="$TMP_LIST_DIR/all_runs.txt"; printf '%s\n' "${runFiles[@]}" > "$LIST_ALL"
FINAL="$OUTPUT_DIR/${RUN_MERGED_PREFIX}_total.root"; rm -f "$FINAL"

if [[ $SUBMODE == condor ]]; then
  SUB="$TMP_LIST_DIR/final_merge.sub"; : >"$SUB"
  cat >>"$SUB" <<EOT
universe = vanilla
executable = $HADD_WRAPPER
output  = $CONDOR_STDOUT/final.\$(Cluster).\$(Process).out
error   = $CONDOR_STDERR/final.\$(Cluster).\$(Process).err
log     = $CONDOR_LOGDIR/final.\$(Cluster).\$(Process).log
request_memory  = $REQUEST_MEMORY
getenv          = True
stream_output   = True
stream_error    = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
arguments = $LIST_ALL $FINAL
queue 1
EOT
  (( DEBUG )) && { say "Submitting grand‑total job (DEBUG preview):"; sed 's/^/    /' "$SUB"; }
  condor_submit "$SUB"
  good "Grand‑total Condor job submitted"
else
  say "Running grand‑total merge locally → $(basename "$FINAL")"
  hadd -v 3 -f "$FINAL" @"$LIST_ALL"
  good "DONE – $(ls -lh "$FINAL")"
fi
