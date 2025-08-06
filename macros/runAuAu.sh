#!/usr/bin/env bash
##############################################################################
#  runAuAu.sh  –  single‑entry driver for the Run‑24/25 Au+Au QA workflow
#
#  QUICK REFERENCE ──────────────────────────────────────────────────────────
#  Local machine
#    ./runAuAu.sh                           : analyse every run found locally
#    ./runAuAu.sh testRun                   : analyse the very first run only
#    ./runAuAu.sh testCombined <N>          : top‑N runs → hadd → full QA pass
#    ./runAuAu.sh runFromCurrentCombined    : QA pass on an existing COMBINED file
#
#  sPHENIX analysis node
#    ./runAuAu.sh fromSPHENIXnode condorTest               : smoke‑test (1 HTCondor job)
#    ./runAuAu.sh fromSPHENIXnode condor                   : full grid (1 job / run)
#    ./runAuAu.sh fromSPHENIXnode haddAndFinalize local    : **local** merge + combined QA
#      ├─ optional  haddAndFinalize        → submit only the merge + combined QA job (Condor)
#      ├─ optional  wipePrevPlots          → delete **all** PNGs under $OUTPUT_DIR before run
#      └─ examples
#             ./runAuAu.sh fromSPHENIXnode condor haddAndFinalize        # remote final pass
#             ./runAuAu.sh fromSPHENIXnode condorTest wipePrevPlots      # smoke‑test, clean PNGs
#
#  FLAGS / OPTIONS ──────────────────────────────────────────────────────────
#    haddAndFinalize          remote merge + combined QA (Condor)
#    local                    when placed **after haddAndFinalize**, do that step locally
#    wipePrevPlots            global wipe of $OUTPUT_DIR before submission
#    VERBOSE=1|2              extra shell diagnostics   (or add ‑v for level 1)
#
#  CLEAN‑UP MATRIX ──────────────────────────────────────────────────────────
#                        | tmp submit dir | log/out/err | outputPlots PNG tree
#    --------------------+---------------+-------------+----------------------
#    condor / condorTest | wiped always  | run‑IDs only| untouched
#    +wipePrevPlots      | ″             | ″           | **fully wiped**
#    haddAndFinalize     | ″             | generic     | untouched
#    haddAndFinalize local| N/A          | N/A         | untouched
#
#  PATHS  (auto‑detected) ───────────────────────────────────────────────────
#    RUN_LOCATION=local     → $HOME/Desktop/auauAnalysis/…
#    RUN_LOCATION=sphenix   → /sphenix/u/$USER/scratch/emcalSEPDcorrelations
#
#  HTCondor artefacts live in   $PROJECT_BASE/{log | stdout | error}
##############################################################################
set -euo pipefail

# ────────────────────────────────────────────────────────────────────────────
# 0.  Global flags
# ────────────────────────────────────────────────────────────────────────────
macro="analyzeRun24or25auau.cpp"   # C++ macro to build/run
verbose="false"
clean_output="true"
finalOnly="false"
stageArg=""
localFinalize="false"          # ← run merge + QA locally
wipePlots="false"              # ← wipe outputPlots only if user adds ‘wipePrevPlots’

# ---------- numeric verbosity (environment or first positional) ----------
if [[ ${1:-} == VERBOSE=* ]]; then
    export VERBOSE="${1#VERBOSE=}"     # take the number after '='
    shift                              # drop this pseudo‑argument
fi
: "${VERBOSE:=0}"                      # default if nothing supplied

# ---------- script‑level noisy/quiet switch ------------------------------
if (( VERBOSE >= 1 )); then
    verbose="true"                     # enables extra shell messages
fi

# legacy --verbose / -v still sets VERBOSE=1
if [[ ${1:-} == "--verbose" || ${1:-} == "-v" ]]; then
    verbose="true"; export VERBOSE=1; shift
fi

# → enable bash trace when VERBOSE ≥2  (shows every shell command)
if (( VERBOSE >= 2 )); then
    set -x
fi

# ────────────────────────────────────────────────────────────────────────────
# 1.  MODE PARSER
# ────────────────────────────────────────────────────────────────────────────
mode="${1:-}"
test_arg="false"        # C++ parameter #1
sample_arg="-1"         # C++ parameter #2

# helper banner
clr_cyan=$'\033[1;36m'; clr_end=$'\033[0m'
step_banner() { printf "${clr_cyan}==========  %s  ==========${clr_end}\n" "$*"; }

step_banner "Argument parsing"
case "${mode}" in
  "") ;;                                           # desktop – full suite
  fromLocalNode)      export RUN_LOCATION=local  ; shift ;;
  fromSPHENIXnode)
          export RUN_LOCATION=sphenix            # use scratch tree
          shift                                  # drop the keyword itself

          # ――― default submission mode ────────────────────────────────
          mode="condor"

          # ――― process all following keywords in any order ―───────────
          while [[ $# -gt 0 ]]; do
                case "$1" in
                    condor|condorTest)
                        mode="$1"
                        ;;
                    haddAndFinalize)
                        finalOnly="true"           # activate multi‑stage mode
                        if [[ ${2:-} =~ ^stage[123]$ ]]; then
                            stageArg="$2"          # remember stage1 / stage2 / stage3
                            shift                  # consume the token
                        fi
                        ;;
                    local)
                        localFinalize="true"       # *local* final‑pass
                        ;;
                    wipePrevPlots)
                        wipePlots="true"
                        ;;
                    --)  shift; break ;;           # end‑of‑options marker
                    *)   break ;;                  # first non‑option → stop parsing
                esac
                shift
          done

          # if the user asked for the local variant, force a dedicated mode
          if [[ "${localFinalize}" == "true" ]]; then
                mode="localFinalize"
          fi
          ;;
  testRun)            test_arg="true"            ; shift ;;
  testCombined)       [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] \
                         || { echo "Usage: … testCombined <N>"; exit 1; }
                      sample_arg="$2" ; shift 2 ;;
  runFromCurrentCombined)
                      export COMBINED_ONLY=1; clean_output="false"; shift ;;
  condorTest)         export RUN_LOCATION=sphenix; mode="condorTest"; shift ;;
  condor)             export RUN_LOCATION=sphenix; mode="condor"    ; shift ;;
  haddAndFinalizeLocal)
                      export RUN_LOCATION=sphenix; localFinalize="true"; shift ;;
  *) ;;
esac

# ────────────────────────────────────────────────────────────────────────────
# 2.  Helper logger functions
# ────────────────────────────────────────────────────────────────────────────
clr_blu=$'\033[1;34m'; clr_red=$'\033[1;31m'; clr_end=$'\033[0m'
note() { printf "${clr_blu}[INFO]${clr_end}  %s\n" "$*"; }
warn() { printf "${clr_red}[WARN]${clr_end}  %s\n" "$*"; }
die () { printf "${clr_red}[FATAL]${clr_end} %s\n" "$*" >&2; exit 2; }



##############################################################################
#  stage1_prepare  – full Missing‑SEB diagnostics *and* good‑run list
#  Arguments :  $1 = INPUT_DIR   $2 = PLOTS_DIR   $3 = RUNLIST_OUT
#  Output    :
#      • ASCII summary table (stdout)
#      • PNG  bar chart   →  $PLOTS_DIR/Combined/MissingSEB_distribution.png
#      • Good‑run list   →  $RUNLIST_OUT      (one run‑ID per line)
##############################################################################
stage1_prepare() {
    local input_dir="$1"
    local plots_dir="$2"
    local runlist="$3"
    rm -f "${runlist}"; touch "${runlist}"

    declare -A sebCnt          # SEB## → counter
    declare -A badRuns         # runID → 1
    local goodCnt=0

    # ── scan per‑run MissingSEB.txt files ────────────────────────────────
    for f in "${input_dir}"/output_*.root; do
        [[ "$(basename "$f")" == "output_ALL_COMBINED.root" ]] && continue
        local run="${f##*/}"; run="${run#output_}"; run="${run%.root}"
        local miss="${plots_dir}/${run}/MissingSEB.txt"

        if [[ -s "${miss}" ]]; then
            badRuns["$run"]=1
            while read -r tok; do
                [[ "${tok}" == SEB* ]] || continue
                (( sebCnt["${tok}"]++ ))
            done < <(tail -n +2 "${miss}")      # skip header
        else
            echo "$run" >> "${runlist}"
            ((goodCnt++))
        fi
    done

    # ── textual summary (identical layout to C++ macro) ──────────────────
    step_banner "Missing SEB summary (Stage‑1)"
    local nBad=${#badRuns[@]}  nBad1=0  nBadMul=0
    for run in "${!badRuns[@]}"; do
        local nTok=$(grep -o 'SEB[0-9]\+' "${plots_dir}/${run}/MissingSEB.txt" | wc -l)
        (( nTok == 1 ? nBad1++ : nBadMul++ ))
    done
    printf "${clr_blu}Runs with ≥1 missing SEB : %s${clr_end}\n" "${nBad}"
    printf "   ├─ exactly one SEB    : %s\n" "${nBad1}"
    printf "   └─ multiple SEBs      : %s\n" "${nBadMul}"
    printf "SEB │ Runs\n"; printf -- '----+-----\n'
    for i in {0..15}; do
        key=$(printf "SEB%02d" "$i")
        val=${sebCnt[$key]:-0}
        printf '%-3s │ %d\n' "$key" "$val"
    done
    printf -- '===========\n'
    note "Good‑run list built – ${goodCnt} runs → ${runlist}"

    # ── bar chart via an inline ROOT macro ───────────────────────────────
    local png="${plots_dir}/Combined/MissingSEB_distribution.png"
    mkdir -p "$(dirname "${png}")"
    local setBins=""
    for i in {0..15}; do
        key=$(printf "SEB%02d" "$i")
        cnt=${sebCnt[$key]:-0}
        setBins+="h.SetBinContent($((i+1)),$cnt);"
    done

    root -l -b -q <<EOF >/dev/null 2>&1
{
   gROOT->SetBatch();
   TH1I h("h","Runs with missing SEB;SEB index;Number of runs",16,-0.5,15.5);
   ${setBins}
   for(int i=0;i<16;++i) h.GetXaxis()->SetBinLabel(i+1,Form("SEB%02d",i));
   h.SetFillColor(kAzure+1); h.SetBarWidth(0.8); h.SetBarOffset(0.1);
   TCanvas c("c","",800,500); c.SetGridy(); h.Draw("bar2");
   gSystem->mkdir("$(dirname "${png}")",true);
   c.SaveAs("${png}");
}
EOF
    note "Bar chart written → ${png}"
}


submit_group_hadd() {        # idx  outfile  infile...
    local idx="$1"; shift
    local ofile="$1"; shift
    local sub="${TMP_BASE}/group_${idx}.sub"
    cat > "$sub" <<EOF
universe        = vanilla
executable      = /bin/bash
arguments       = -c "hadd -f ${ofile} $*"
output          = ${TMP_BASE}/group_${idx}.out
error           = ${TMP_BASE}/group_${idx}.err
log             = ${TMP_BASE}/group_${idx}.log
getenv          = True
request_memory  = 2GB
+JobFlavour     = "tomorrow"
queue
EOF
    condor_submit "$sub" || die "condor_submit failed for group $idx"
}

submit_final_hadd() {        # outfile  infile...
    local ofile="$1"; shift
    local sub="${TMP_BASE}/final_hadd.sub"
    cat > "$sub" <<EOF
universe        = vanilla
executable      = /bin/bash
arguments       = -c "hadd -f ${ofile} $* && \
                      export COMBINED_ONLY=1 EXTERNAL_HADD=1 && \
                      root -l -b -q -e 'gSystem->SetBuildDir(\"${TMPDIR:-/tmp}\",kTRUE)' ${macro}+Ok(false,-1)"
output          = ${TMP_BASE}/final_hadd.out
error           = ${TMP_BASE}/final_hadd.err
log             = ${TMP_BASE}/final_hadd.log
getenv          = True
request_memory  = 4GB
+JobFlavour     = "tomorrow"
queue
EOF
    condor_submit "$sub" || die "condor_submit failed for final hadd"
}

##############################################################################
# perform_hadd  – build output_ALL_COMBINED.root with the classic ROOT hadd
# Skips runs that have MissingSEB.txt (i.e. “bad” runs).
##############################################################################
perform_hadd () {
    local input_dir="$1"        # e.g.  /sphenix/u/$USER/…/output
    local plots_dir="$2"        # e.g.  /sphenix/tg/tg01/…/outputPlots
    local combined="${input_dir}/output_ALL_COMBINED.root"

    note "hadd – collecting per‑run ROOT files"

    good_files=()                 # runs that will enter the merge
    declare -A sebCnt             # SEB##  → counter
    declare -A badRuns            # runID  → 1  (at least one missing SEB)

    # ── scan every per‑run ROOT file ──────────────────────────────────────
    for f in "${input_dir}"/output_*.root ; do
        [[ "$(basename "$f")" == "output_ALL_COMBINED.root" ]] && continue   # skip target file
        run="${f##*/}"; run="${run#output_}"; run="${run%.root}"

        missFile="${plots_dir}/${run}/MissingSEB.txt"
        if [[ -s "${missFile}" ]]; then
            badRuns["$run"]=1
            while read -r tok ; do
                [[ "${tok}" == SEB* ]] || continue
                (( sebCnt["${tok}"]++ ))
            done < <(tail -n +2 "${missFile}")      # skip header line
            warn "  ↳ run ${run} excluded (MissingSEB.txt not empty)"
            continue
        fi

        good_files+=( "${f}" )
    done

    # ── textual summary identical to the C++ macro ────────────────────────
    step_banner "Missing SEB summary"

    local nBad=${#badRuns[@]}
    local nBad1=0 nBadMul=0
    for run in "${!badRuns[@]}"; do
        nTok=$(grep -o 'SEB[0-9]\+' "${plots_dir}/${run}/MissingSEB.txt" | wc -l)
        (( nTok == 1 ? nBad1++ : nBadMul++ ))
    done

    printf "${clr_blu}Runs with ≥1 missing SEB : %s${clr_end}\n" "${nBad}"
    printf "   ├─ exactly one SEB    : %s\n" "${nBad1}"
    printf "   └─ multiple SEBs      : %s\n" "${nBadMul}"
    printf "SEB │ Runs\n"
    printf -- '----+-----\n'
    for i in {0..15}; do
        key=$(printf "SEB%02d" "$i")
        val=${sebCnt[$key]:-0}        # default to 0 when unset
        printf -- '%-3s │ %d\n' "$key" "$val"
    done
    printf -- '===========\n'

    # ── abort if we ended up with <2 good runs ────────────────────────────
    (( ${#good_files[@]} >= 2 )) || { warn "Need ≥2 good runs – aborting hadd"; return 1; }

    # ── size bookkeeping + merge ──────────────────────────────────────────
    total_mb=0
    for g in "${good_files[@]}"; do
        sz=$( { stat -c%s "$g" 2>/dev/null || stat -f%z "$g"; } ) || sz=0
        (( total_mb += sz/1024/1024 ))
    done
    note "Merging ${#good_files[@]} runs  (≈${total_mb} MB total)"
    rm -f "${combined}"
    hadd -f "${combined}" "${good_files[@]}" || die "hadd failed"
    note "Combined file created → ${combined}"
}


# ────────────────────────────────────────────────────────────────────────────
# 3.  Condor submission pathway
# ────────────────────────────────────────────────────────────────────────────
if [[ "${mode}" == "condor" || "${mode}" == "condorTest" ]]; then
    step_banner "Condor submission configuration"

    # ------------------------------------------------------------------
    #  haddAndFinalize  →  exactly one job that merges all runs and
    #                      performs the combined QA pass
    # ------------------------------------------------------------------
    if [[ "${finalOnly}" == "true" ]]; then
        # ------------------------------------------------------------------
        # Multi‑stage haddAndFinalize
        #    • stage1 : build run list          (SEB check only, no hadd)
        #    • stage2 : 10‑run group hadd jobs  (Condor)
        #    • stage3 : final hadd + QA pass    (local | condor)
        # ------------------------------------------------------------------
        stage="${stageArg:-stage1}"            # default when user omits token

        PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
        INPUT_DIR="${PROJECT_BASE}/output"
        OUTPUT_DIR="/sphenix/tg/tg01/bulk/jbennett/GLOBAL_QA/outputPlots"
        TMP_BASE="${PROJECT_BASE}/tmp_hadd_multistage"
        RUNLIST="${TMP_BASE}/runlist_stage2.txt"
        GROUP_DIR="${TMP_BASE}/groups"
        mkdir -p "${TMP_BASE}" "${GROUP_DIR}"

        case "${stage}" in
            stage1)
                step_banner "Stage‑1  –  SEB scan + good‑run list"
                stage1_prepare  "${INPUT_DIR}"  "${OUTPUT_DIR}"  "${RUNLIST}"
                note "Stage‑1 finished – launch Stage‑2 once the group jobs can run."
                ;;

            stage2)
                step_banner "Stage‑2  –  submit group hadd jobs"
                [[ -f "${RUNLIST}" ]] || die "Stage‑2: run list ${RUNLIST} missing"
                mapfile -t allRuns < "${RUNLIST}"
                (( ${#allRuns[@]} )) || die "Stage‑2: run list empty"

                grpIdx=0; grpRuns=()
                for r in "${allRuns[@]}"; do
                    grpRuns+=( "${INPUT_DIR}/output_${r}.root" )
                    if (( ${#grpRuns[@]} == 10 )); then
                        submit_group_hadd "${grpIdx}" "${GROUP_DIR}/group_${grpIdx}.root" "${grpRuns[@]}"
                        grpRuns=(); ((grpIdx++))
                    fi
                done
                if (( ${#grpRuns[@]} )); then
                    submit_group_hadd "${grpIdx}" "${GROUP_DIR}/group_${grpIdx}.root" "${grpRuns[@]}"
                fi
                rm -f "${RUNLIST}"
                note "Stage‑2 submitted – temporary run list removed."
                ;;

            stage3)
                step_banner "Stage‑3  –  final hadd and QA pass"
                shift || true                         # consume 'stage3'
                runLocal="no"
                [[ ${1:-} == local ]]  && { runLocal="yes"; shift; }
                [[ ${1:-} == condor ]] && { runLocal="no" ; shift; }

                mapfile -t grpFiles < <(ls "${GROUP_DIR}"/group_*.root 2>/dev/null | sort)
                (( ${#grpFiles[@]} )) || die "Stage‑3: no group files – finish Stage‑2 first"

                COMBINED="${INPUT_DIR}/output_ALL_COMBINED.root"

                if [[ "${runLocal}" == "yes" ]]; then
                    note "Local final hadd – merging ${#grpFiles[@]} group files"
                    rm -f "${COMBINED}"
                    hadd -f "${COMBINED}" "${grpFiles[@]}" || die "hadd failed"
                    export COMBINED_ONLY=1 EXTERNAL_HADD=1
                    root -l -b -q \
                         -e "gSystem->SetBuildDir(\"${TMPDIR:-/tmp}\",kTRUE)" \
                         "${macro}+Ok(false,-1)"
                else
                    note "Submitting final hadd as one Condor job"
                    submit_final_hadd "${COMBINED}" "${grpFiles[@]}"
                fi

                # clean‑up when combined file exists (local path) or when job finishes (condor)
                [[ -f "${COMBINED}" ]] && { rm -f "${GROUP_DIR}"/group_*.root; rmdir "${GROUP_DIR}" 2>/dev/null || true; }
                ;;

            *)
                die "Unknown haddAndFinalize stage '${stage}' (use stage1|stage2|stage3)"
                ;;
        esac
        exit 0
    fi

cat > "${SUBMIT_DIR}/haddAndFinalize.sub" <<EOS
        universe        = vanilla
        executable      = ${EXEC_WRAPPER}
        arguments       = --final
        output          = ${STDOUT_DIR}/haddAndFinalize.out
        error           = ${STDERR_DIR}/haddAndFinalize.err
        log             = ${LOG_DIR}/haddAndFinalize.log
        getenv          = True
        request_memory  = 8GB
        +JobFlavour     = "tomorrow"
        queue
EOS
        note "Submitting finalize HTCondor job"
        condor_submit "${SUBMIT_DIR}/haddAndFinalize.sub" \
            || die "condor_submit failed for finalize job"
        note "Finalize job submitted – exiting wrapper."
        exit 0
    fi

    note "Starting Condor submission (${mode})"

    PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
    EXEC_WRAPPER="${PROJECT_BASE}/macros/runAuAuExecutable.sh"

    SUBMIT_DIR="${PROJECT_BASE}/tmp_condor_submit"
    LOG_DIR="${PROJECT_BASE}/log"
    STDOUT_DIR="${PROJECT_BASE}/stdout"
    STDERR_DIR="${PROJECT_BASE}/error"

    INPUT_DIR="${PROJECT_BASE}/output"
    OUTPUT_DIR="/sphenix/tg/tg01/bulk/jbennett/GLOBAL_QA/outputPlots"

    [[ -x "${EXEC_WRAPPER}" ]] || die "Wrapper ${EXEC_WRAPPER} missing or not executable"

    mapfile -t roots < <(ls "${INPUT_DIR}"/output_*.root 2>/dev/null | sort)
    [[ ${#roots[@]} -gt 0 ]] || die "No ROOT files in ${INPUT_DIR}"
    note "Found ${#roots[@]} ROOT input files"

    # pick the ROOT file that has the *most* events (proxy = size)
    if [[ "${mode}" == "condorTest" ]]; then
        note "condorTest → selecting run with the highest statistics"

        largest=""
        maxSize=0
        for f in "${roots[@]}"; do
            sz=$( { stat -c%s "$f" 2>/dev/null || stat -f%z "$f"; } ) || sz=0
            (( sz > maxSize )) && { maxSize=$sz; largest="$f"; }
        done

        [[ -n "${largest}" ]] || die "Could not determine the largest ROOT file"
        roots=( "${largest}" )

        runPick=$(basename "${largest}")
        note "condorTest → will submit only ${runPick}  (size $((maxSize/1024/1024)) MB)"
    fi

    note "Cleaning previous submission artefacts"

    rm -rf "${SUBMIT_DIR:?}/"* 2>/dev/null || true

    if [[ "${wipePlots}" == "true" && "${finalOnly}" != "true" ]]; then
        warn "wipePrevPlots requested – deleting existing ${OUTPUT_DIR}"
        rm -rf "${OUTPUT_DIR:?}/"* 2>/dev/null || true
    fi

    for rf in "${roots[@]}"; do
        run="$(basename "${rf}" .root)"; run="${run#output_}"
        rm -f  "${LOG_DIR}/${run}.log"    2>/dev/null || true
        rm -f  "${STDOUT_DIR}/${run}.out" 2>/dev/null || true
        rm -f  "${STDERR_DIR}/${run}.err" 2>/dev/null || true
    done

    mkdir -p "${SUBMIT_DIR}" "${LOG_DIR}" "${STDOUT_DIR}" "${STDERR_DIR}" "${OUTPUT_DIR}"

    step_banner "Submitting per‑run jobs"
    runLogs=()
    expRuns=${#roots[@]}

    for rf in "${roots[@]}"; do
        run="$(basename "${rf}" .root)"; run="${run#output_}"
        sub="${SUBMIT_DIR}/${run}.sub"

        note "→ Submitting run ${run}"
        note "   input    : ${rf}"
        note "   out‑dir  : ${OUTPUT_DIR}/${run}"

cat > "${sub}" <<EOS
        universe        = vanilla
        executable      = ${EXEC_WRAPPER}
        arguments       = ${rf}  ${OUTPUT_DIR}/${run}
        output          = ${STDOUT_DIR}/${run}.out
        error           = ${STDERR_DIR}/${run}.err
        log             = ${LOG_DIR}/${run}.log
        getenv          = True
        request_memory  = 4GB
        +JobFlavour     = "tomorrow"
        queue
EOS
        
        condor_submit "${sub}" || die "condor_submit failed for ${run}"
        runLogs+=( "${LOG_DIR}/${run}.log" )
    done
    note "All ${expRuns} run‑jobs submitted – exiting submission block."
    exit 0
fi

# ──────────────────────────────────────────────────────────────────
# 3‑B.  haddAndFinalizeLocal  –  run merge + combined QA *locally*
# ──────────────────────────────────────────────────────────────────
if [[ "${localFinalize}" == "true" ]]; then
    step_banner "Local hadd + finalize"

    PROJECT_BASE="/sphenix/u/${USER}/scratch/emcalSEPDcorrelations"
    INPUT_DIR="${PROJECT_BASE}/output"
    OUTPUT_DIR="/sphenix/tg/tg01/bulk/jbennett/GLOBAL_QA/outputPlots"

    perform_hadd "${INPUT_DIR}" "${OUTPUT_DIR}" || exit 2

    # run the combined QA pass on the freshly‑hadded file
    export COMBINED_ONLY=1          # skip per‑run loops
    export EXTERNAL_HADD=1          # tells the C++ macro not to rebuild
    root -l -b -q \
         -e "gSystem->SetBuildDir(\"${TMPDIR:-/tmp}\",kTRUE)" \
         "${macro}+Ok(false,-1)"

    note "Local finalize completed."
    exit 0
fi

# --------------------------------------------------------------------------
# 4.  Optional QA‑module filter on the command line
# --------------------------------------------------------------------------
if [[ $# -ge 1 ]]; then
    export QA_ONLY="$1"; shift
    note "QA_ONLY filter applied → ${QA_ONLY}"
fi

# --------------------------------------------------------------------------
# 5.  Desktop / interactive execution
# --------------------------------------------------------------------------
step_banner "Interactive / desktop execution"
export LDFLAGS="${LDFLAGS:-} -Wl,-no_warn_duplicate_libraries"
root_flags=(-l -b -q)

build_dir="$(mktemp -d "${TMPDIR:-/tmp}/aclic_build_XXXXXX")"
mkdir -p "${build_dir}"

root_cmd=(
  root "${root_flags[@]}"
  -e "gSystem->SetBuildDir(\"${build_dir}\",kTRUE)"
  "${macro}+Ok(${test_arg},${sample_arg})"
)
[[ "${verbose}" == "true" ]] && note "+ ${root_cmd[*]}"

output_root="${HOME}/Desktop/auauAnalysis/emcalSEPDcorrelations/output"
[[ "${RUN_LOCATION:-local}" == "local" ]] || clean_output="false"

if [[ "${clean_output}" == "true" && -d "${output_root}" ]]; then
    warn "Cleaning old local output under ${output_root}"
    rm -rf "${output_root:?}/"*
fi

macro_dir="$(cd "$(dirname "${macro}")" && pwd)"
macro_base="$(basename "${macro%.*}")_cpp"

rm -rf "${macro_dir}/.aclic_build"
mkdir -p "${macro_dir}/.aclic_build"

find "${macro_dir}" -maxdepth 1 -type f -name "${macro_base}_ACLiC_*" -delete
rm -f "${macro_dir}/${macro_base}".{so,d}

note "Launching ROOT compilation / execution"
"${root_cmd[@]}" 2>&1 | grep -vE '^ld: warning: duplicate -rpath .+ ignored$'
note "ROOT macro completed"
