# CALO × sEPD × MBD – **Run‑25 Au+Au (CALOFITTING)**  
*(patsfan753 · July 2025 · BNL sPHENIX farm)*

End‑to‑end instructions to turn raw **CALOFITTING** DSTs into merged,
analysis‑ready ROOT files **and** publication‑quality flow‑vₙ plots.

---

## 0 · Folder Structure

```
emcalSEPDcorrelations/
├─ src/                 # C++ plugin
│   ├─ autogen.sh  configure.ac  Makefile.am
│   ├─ emcal_sepdCorrelator.{cc,h}
├─ scripts/             # pipeline helpers
│   ├─ makeDstLists.sh                  # step ①
│   ├─ run_auau_run3_qa_submit.sh       # step ② (+ splitRunList)
│   ├─ run_auau_run3_qa.sh              # executed on every Condor slot
│   └─ merge_data.sh                    # steps ③ & ④
├─ macros/
│   ├─ Fun4All_emcalSEPDcorrelator.C    # executes correlator + jets
│   └─ analyzeRun24or25auau.cpp         # offline plotting
└─ README.md            # you are here
```

---

## 1 · Build & Install the correlator plugin (one‑time)

```bash
cd src
prefix=${1:-$MYINSTALL}            # optional install path
./autogen.sh  --prefix="$prefix"
./configure   --prefix="$prefix"
make -j8 && make install
```

---

## 2 · Create per‑run DST lists *(~1 min)*

```bash
cd "$PROJECT_BASE"                       # e.g. /sphenix/u/$USER/scratch/emcalSEPDcorrelations
./makeDstLists.sh run25auau caloFitting forceFileList (having issues with createDSTList this works)
```

Results → `dst_list/DST_CALOFITTING_…-<run>.list`  
(one absolute path per line).

---

## 3 · Submit QA jobs to HTCondor

```bash
./run_auau_run3_qa_submit.sh run25auau caloFitting condor
```

> **Variants**  
> `local` · `condorTest` · `condor firstTen` · `condor round <N>`

Chunk size = 2 DSTs.  
Outputs →  
`/sphenix/tg/tg01/bulk/jbennett/emcalSEPDcorrelations/<run>/output_<tag>.root`

---

## 4 · Merge chunk files **per run**

```bash
DEBUG=1 ./merge_data.sh condor
```

Creates `output/output_<run>.root`

---

## 5 · Merge all runs into one file

```bash
# local
DEBUG=1 ./merge_data.sh addRuns
# or Condor
DEBUG=1 ./merge_data.sh addRuns condor
```

Result: `output/output_total.root`

---

## 6 · Offline analysis & plotting

```bash
# build once (ACLiC optimisation) then launch
chmod +x runAuAu.sh
./runAuAu.sh                   # default: all QA modules on all runs
```

### 6.1  Selecting QA modules at run‑time  

The last positional argument can be **a comma‑separated list of QA tags**.  
Only the listed modules run; everything else is skipped.  
If you omit the list, the macro behaves exactly as before (runs **all** modules).

| Tag            | QA module instantiated                            |
| -------------- | ------------------------------------------------- |
| `correlations` | detector‑correlation maps (CorrQA)                |
| `hcal`         | HCal QA (IHCal / OHCal / totalHCal)               |
| `mbd`          | MBD QA                                            |
| `sepd`         | sEPD QA (event‑plane, tile, …)                    |
| `sepdother`    | miscellaneous sEPD checks                         |
| `jetqa`        | jet‑trigger QA                                    |
| `eventqa`      | global event‑quality checks                       |
| `triggerqa`    | trigger‑counter summaries                         |
| `pi0`          | π⁰ / η invariant‑mass QA                          |
| `emcal`        | EMCal QA (non‑mass plots)                         |
| `vn`           | flow‑vₙ plots                                     |

> **Syntax**  
> `./runAuAu.sh [--verbose] [mode] [mode‑args] [qa_tag_1,qa_tag_2,…]`

### 6.2  Typical calls

```bash
# CorrQA and JetQA only, all runs
./runAuAu.sh correlations,jetqa

# first run only, HCal and π0 modules
./runAuAu.sh testRun hcal,pi0

# merge top‑5 runs, then CorrQA + HCal + vₙ plots
./runAuAu.sh testCombined 5 correlations,hcal,vn
```

Key outputs:

```
output/Combined/<trigger>/
   ├─ vn/FlowQA/…          (if ‘vn’ selected)
   │   ├─ v2_ALL_S_allCent_MB.png
   │   ├─ v2_cent20_40_MB.png
   │   └─ vbar2_ALL_S_vsCent_MB.png
   ├─ EMCal/…              (if ‘emcal’ selected)
   └─ InvariantMassSummary.csv
```


---

## 7 · Optional: split huge run lists

```bash
./run_auau_run3_qa_submit.sh run25auau splitRunList run25GoldenRuns.txt
./run_auau_run3_qa_submit.sh run25auau caloFitting condor round 2
```

`runSegment_run25auau_<n>.txt` limits submissions to ≤10 k jobs per round.

---

## 8 · Directory cheat‑sheet

```
$PROJECT_BASE/
├─ dst_list/                 # ≤1 file per run (step ①)
├─ tmp_condor_lists/         # per‑chunk lists
├─ log/  stdout/  error/     # Condor I/O
├─ output/
│   ├─ output_<run>.root     # after step ③
│   ├─ output_total.root     # after step ⑤
│   └─ Combined/…            # plots & CSV
└─ vn/FlowQA/                # inside Combined/<trigger>/
```

Happy analysing 🎉
