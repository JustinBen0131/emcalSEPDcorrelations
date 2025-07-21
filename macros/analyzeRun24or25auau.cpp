// analyzeRun24or25auau.cpp  – ROOT ≥ 6, C++17
// ===============================================================
//  • Pass‑0: catalogue every histogram in the file  (console + .txt)
//  • Pass‑1: modular QA (EMCal, HCal, sEPD, MBD, correlations,
//            #pi0 invariant‑mass spectra, …) with automatic centrality
//            slice replication and North/South map fusion.
// ===============================================================

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <random>
#include <thread>
#include <stdexcept>
#include <sstream>
#include <TStyle.h>
#include <TH2Poly.h>
#include <TFileMerger.h>
#include <filesystem>
#include <TH2.h>
#include <TH3.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TSequentialExecutor.hxx>    // defines ROOT::TSequentialExecutor
#include <Math/MinimizerOptions.h>         // defines ROOT::Math::MinimizerOptions
#include <TLegend.h>                       // full definition of TLegend
#include <iostream>
#include <TGraphErrors.h>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <cstdint>   // uintptr_t cast
#include <algorithm>
#include <cmath>

using std::string;
namespace fs = std::filesystem;

// ───────────────────────────────────────────────
//            ▼  USER SETTINGS  ▼
// ───────────────────────────────────────────────
namespace {
  /* directory that already contains all run‑merged ROOT files
   * produced by your Condor step:  output_00066XXX.root, …            */
  const fs::path kInputDir =
      "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/input/output";

  /* parent directory where PNGs/CSVs will appear:
   *   <kOutputDir>/<run‑number>/…    (per‑run)
   *   <kOutputDir>/Combined/…        (after hadd)                      */
  const fs::path kOutputDir =
      "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/output";
}
std::string kInputFile   = "";   // gets filled inside the loop
std::string kOutputBase  = "";   // ditto
std::set<string> kTriggersWanted{ "MBD_NandS_geq_2" };
/* ----------  PLOT CANVAS SIZES  -----------
 * Set a value ≤0 to keep the old automatic size
 * (nEta × px for width, nPhi × px for height).
 * -------------------------------------------------------- */
constexpr int kHCalCanvasW = 800;   // IHCal / OHCal width   [px]
constexpr int kHCalCanvasH = 1200;  // IHCal / OHCal height  [px]

constexpr int kEMCalCanvasW = 800;  // EMCal  hit‑map width  [px]
constexpr int kEMCalCanvasH = 1200; // EMCal  hit‑map height [px]
/*  <<<   #pi0‐fit master switch   >>>                                         *
 *  false  → spectra are drawn, but *no* TF1 fit is attempted and            *
 *           InvariantMassSummary.csv is left empty (except header).         *
 *  true   → run the Gaussian‑plus‑poly fit and fill CSV.                    */
constexpr bool kDoPi0Fit = true;
// ───────────────────────────────────────────────

// ---------- helper: collect *.root files, sorted -------------------
static std::vector<fs::path> listRunFiles(const fs::path& dir)
{
  std::vector<fs::path> v;
  std::regex re(R"(output_([0-9]{8})\.root)");
  for (auto& e : fs::directory_iterator(dir))
    if (e.is_regular_file() && e.path().extension() == ".root" &&
        std::regex_match(e.path().filename().string(), re))
      v.push_back(e.path());
  std::sort(v.begin(), v.end(),
            [](const fs::path& a, const fs::path& b){ return a.string() < b.string(); });
  return v;
}


// ╔══════════════════════════════════════════════╗
// ║                1.  LOGGING                   ║
// ╚══════════════════════════════════════════════╝
namespace term {
  constexpr const char* CLR_RST  = "\033[0m";
  constexpr const char* CLR_BOLD = "\033[1m";
  constexpr const char* CLR_CYAN = "\033[36m";
  constexpr const char* CLR_GRN  = "\033[32m";
  constexpr const char* CLR_YEL  = "\033[33m";
  constexpr const char* CLR_RED  = "\033[31m";
}

namespace log {
  inline void banner(const string& m)
  {
    std::cout << "\n" << term::CLR_BOLD << term::CLR_CYAN
              << "══════════════════════════════════════════════════════════════\n"
              << "  " << m << "\n"
              << "══════════════════════════════════════════════════════════════"
              << term::CLR_RST << "\n";
  }
  inline void info (const string& m){ std::cout<<term::CLR_CYAN<<m<<term::CLR_RST<<"\n"; }
  inline void ok   (const string& m){ std::cout<<term::CLR_GRN <<m<<term::CLR_RST<<"\n"; }
  inline void warn (const string& m){ std::cout<<term::CLR_YEL <<m<<term::CLR_RST<<"\n"; }
  inline void err  (const string& m){ std::cerr<<term::CLR_RED <<m<<term::CLR_RST<<"\n"; }

  // new ultra‑light trace channel
  inline void trace(const std::string& m)
  {
    std::cout << term::CLR_YEL << "[TRACE] " << m << term::CLR_RST << "\n";
  }
}

// helper --------------------------------------------------------------
inline void ensure_dir(const fs::path& p)
{
  std::error_code ec;
  fs::create_directories(p, ec);
}
inline string sf3(double x){
  std::ostringstream o; if(x==0){o<<"0";return o.str();}
  int e=int(floor(log10(fabs(x)))); o<<std::fixed<<std::setprecision(std::max(0,2-e))<<x;
  string s=o.str(); std::replace(s.begin(),s.end(),'.','p'); return s;
}



// ╔══════════════════════════════════════════════╗
// ║          2.  GEOMETRY  HELPERS               ║
// ╚══════════════════════════════════════════════╝
static inline int sector_from_idx(unsigned ieta,unsigned iphi)
{ if(iphi>=256) return -1; int base=iphi/8; return (ieta<48)?32+base:base; }

static inline int ib_from_idx(unsigned ieta,unsigned)
{
  const int table[12] = {5,4,3,2,1,0,0,1,2,3,4,5};
  return (ieta<96)?table[ieta/8]:-1;
}
inline bool isBadBoard(int sec,int ib){
  return ((sec==50&&ib==1)||(sec==4&&ib==1)||(sec==25&&ib==2));
}

static inline int hcal_sector_from_idx(unsigned, unsigned iphi){ return (iphi<64)?iphi/4:-1; }
static inline int hcal_plate_from_idx(unsigned ieta, unsigned){ return (ieta<24)?ieta/4:-1; }
inline bool isBadHcalPlate(int,int){ return false; }



// ╔══════════════════════════════════════════════╗
// ║ 3.  CENTRALITY‑TAG DISCOVERY                 ║
// ╚══════════════════════════════════════════════╝
using CentList = std::vector<string>;
CentList discoverSlices(TFile* f)
{
  std::set<string> tags;
  std::regex r(R"(_([0-9]+_[0-9]+)_.*$)");
  TIter itDir(f->GetListOfKeys());
  while(auto* kd=dynamic_cast<TKey*>(itDir())){
    if(strcmp(kd->GetClassName(),"TDirectoryFile")) continue;
    TDirectory* t=(TDirectory*)kd->ReadObj();
    TIter itH(t->GetListOfKeys());
    while(auto* kh=dynamic_cast<TKey*>(itH())){
      string n=kh->GetName(); std::smatch m;
      if(std::regex_search(n,m,r)) tags.insert(m[1]);
    }
  }
  if(tags.empty()) tags.insert("Inclusive");
  return {tags.begin(),tags.end()};
}

inline fs::path cPath(fs::path base, const string& slice, fs::path sub)
{
  // first enter the detector / analysis sub‑folder …
  base /= sub;

  // …then append the centrality tier *within* that tree
  base /= (slice == "Inclusive" ? "noCentralityDep"
                                : "Cent_" + slice);
  return base;
}



// ╔══════════════════════════════════════════════╗
// ║ 4.  #pi0 CUT KEY / PEAK FITTER                 ║
// ╚══════════════════════════════════════════════╝
struct CutKey{ float E,chi,asy,pLo,pHi; string trigger; };

bool decodeInvName(const string& n, CutKey& k)
{
  std::regex r(R"(mInv_pt([\-0-9\.]+)to([\-0-9\.]+)_E([0-9\.]+)_chi([0-9\.]+)_asy([0-9\.]+)_(.+))");
  std::smatch m; if(!std::regex_match(n,m,r)) return false;
  k.pLo=stof(m[1]); k.pHi=stof(m[2]); k.E=stof(m[3]); k.chi=stof(m[4]); k.asy=stof(m[5]); k.trigger=m[6];
  return true;
}

// ------------------------------------------------------------------
// Return the bin index whose centre is closest to x (1‑based)
// ------------------------------------------------------------------
static inline int binAt(const TH1* h, double x)
{
    const double bw = h->GetBinWidth(1);
    int idx = static_cast<int>(std::round((x - h->GetXaxis()->GetXmin())/bw)) + 1;
    return std::clamp(idx, 1, h->GetNbinsX());
}

// ╔══════════════════════════════════════════════╗
// ║ 5.  SAVE HELPERS  (with TRACE)               ║
// ╚══════════════════════════════════════════════╝
void save1D(TH1* h,const fs::path& p)
{
  log::trace("save1D → " + p.string());
  TCanvas c; h->SetStats(0); h->Draw();
  ensure_dir(p.parent_path()); c.SaveAs(p.string().c_str());
}
void save2D(TH2* h,const fs::path& p,const char* opt="COLZ")
{
  log::trace("save2D → " + p.string());
  TCanvas c; h->SetStats(0); h->Draw(opt);
  ensure_dir(p.parent_path()); c.SaveAs(p.string().c_str());
}



// ╔══════════════════════════════════════════════╗
// ║ 6.  NS MAP CACHE                             ║
// ╚══════════════════════════════════════════════╝
template<class MPair> class NSCache{
public:
  MPair& operator[](const string& k){ return _c[k]; }
  bool ready(const string& k){ return _c[k].n && _c[k].s; }
  MPair pop(const string& k){ MPair p=_c[k]; _c.erase(k); return p; }
  struct H {
        size_t operator()(const string& s) const noexcept {
            return std::hash<string>{}(s);
        }
  };
private:
  std::unordered_map<string,MPair,H> _c;
};
struct MapPair{ TH2* n=nullptr,*s=nullptr; };



// ╔══════════════════════════════════════════════╗
// ║ 7.  QA  BASE CLASS                           ║
// ╚══════════════════════════════════════════════╝
class QA{
public:
  QA(string trg, fs::path base, const CentList& c): trig(trg), root(base), slices(c){}
  virtual ~QA()=default;
  virtual bool process(TObject*)=0;
protected:
  string trig;
  fs::path root;
  const CentList& slices;
  static string sliceKey(const string& n){
    std::regex r(R"(_([0-9]+_[0-9]+)_)"); std::smatch m;
    return std::regex_search(n, m, r) ? m[1].str() : std::string("Inclusive");
  }
};



// ╔══════════════════════════════════════════════╗
// ║     #pi0   I N V A R I A N T ‑ M A S S   QA  ║
// ╚══════════════════════════════════════════════╝
class Pi0QA : public QA
{
 public:
  /* pT bins that appear in the file names – must match CutKey values */
  const std::vector<std::pair<float,float>> m_ptBins {
      {2,4},{4,6},{6,8},{8,10},{10,12},{12,15},{15,20},{20,30} };

  struct FitInfo {
    std::string slice;  double pLo, pHi;
    double mean, sigma, chi2;  int ndf;
  };
    
    
  struct FitPair {
            std::unique_ptr<TF1> total;
            std::unique_ptr<TF1> poly;

            /* default ctor/dtor */
            FitPair() = default;
            ~FitPair() = default;

            /* deep‑copy ctor */
            FitPair(const FitPair& other)
            {
                if (other.total) total.reset( static_cast<TF1*>( other.total->Clone() ) );
                if (other.poly ) poly .reset( static_cast<TF1*>( other.poly ->Clone() ) );
            }

            /* deep‑copy assignment */
            FitPair& operator=(const FitPair& other)
            {
                if (this != &other) {
                    total.reset( other.total ? static_cast<TF1*>( other.total->Clone() ) : nullptr );
                    poly .reset( other.poly  ? static_cast<TF1*>( other.poly ->Clone() )  : nullptr );
                }
                return *this;
            }

            /* move semantics remain the default */
            FitPair(FitPair&&) noexcept            = default;
            FitPair& operator=(FitPair&&) noexcept = default;
  };


  /* ------------------------------------------------------------------ *
   *  ctor – figure out the run label from the base directory            *
   *  …/<output>/<run>/<trigger>  →  run = last element of parent_path() *
   * ------------------------------------------------------------------ */
  Pi0QA(std::string           t,
        fs::path              b,
        const CentList&       s,
        std::ofstream&        csvFit_) :
        QA(std::move(t), std::move(b), s),
        csvFit(csvFit_)
  {
    /* run label */
    runID = root.parent_path().filename().string();

    fs::path p = root / "EMCal/invMassQA" / "Pi0SignalBackground.csv";
    ensure_dir(p.parent_path());
    csvSB.open(p);
    csvSB << "trigger,cent,pTlo,pThi,slice,windowSigma,sb,err\n";
  }

  ~Pi0QA() override
  {
      writeSummaryPanels();
      writeRunSummary();            // only fires once in the Combined pass
  }

  const auto& fitSummary() const { return _fitSummary; }

  // ----------------------------------------------------------------
  // MAIN ENTRY – called once per histogram
  // ----------------------------------------------------------------
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;
    std::string n = o->GetName();
      
    /* =============================================================== *
     * (A)  UNCUT 2‑D maps  –  handled here once and returned          *
     *      Histogram names:  Minv_vs_Asym_… , Minv_vs_chi2_… , …      *
     * =============================================================== */
   if (o->InheritsFrom(TH2::Class()) &&         // TH2F only
          n.rfind("Minv_vs_", 0) == 0)             // quick prefix test
   {
          /* optional safety: skip empty plots produced by empty events  */
          auto* h2 = static_cast<TH2*>(o);
          if (h2->GetEntries() <= 0 || h2->Integral() <= 0) return false;

          const std::string slice = sliceKey(n);            // "Inclusive" or "lo_hi"
          fs::path subDir = "EMCal/invMassQA/cutQA";        // new folder
          fs::path outPng = cPath(root, slice, subDir) / (n + ".png");
          ensure_dir(outPng.parent_path());

          TCanvas c;
          h2->SetStats(0);
          h2->Draw("COLZ");
          c.SaveAs(outPng.string().c_str());

          return true;                                      // ← histogram consumed
    }

      
    if (n.rfind("mInv_",0)!=0) return false;

    CutKey ck;              // parse name
    if (!decodeInvName(n, ck)) return false;

    const std::string slice = sliceKey(n);
    const bool        pTInt = (ck.pLo<0 || ck.pHi<0);

    //----------------------------------------------------------------
    // 0. Book‑keeping: where will the PNG go?
    //----------------------------------------------------------------
    const std::string combDir = "E"+sf3(ck.E)+"_Chi"+sf3(ck.chi)+"_Asym"+sf3(ck.asy);
    if (cutTag.empty()) {                          // first spectrum ever seen
          cutTag = combDir;
    } else if (combDir != cutTag) {                // ← a NEW cut combination
          writeSummaryPanels();                      // ➊ finish the previous cut
          _centralHists.clear();                     // ➋ reset caches
          _storedFit.clear();
          _storedEtaFit.clear();
          cutTag = combDir;                          // ➌ start the next cut
    }

    /* 1.  base directory = cut‑combination *inside* the proper centrality slice */
    fs::path subDir = "EMCal/invMassQA";
    subDir /= combDir;
    fs::path baseDir = cPath(root, slice, subDir);   // e.g. …/Cent_0_10/…/E2p00_Chi…/

    /* 2.  append the pT‑range only after the centrality folder has been inserted */
    if (!pTInt)
          baseDir /= ("pT_" + sf3(ck.pLo) + "_to_" + sf3(ck.pHi));

    ensure_dir(baseDir);                             // make sure it exists

    /* 3.  final PNG path */
    fs::path outPng = baseDir / (n + ".png");

    //----------------------------------------------------------------
    // 1.  Robust #pi0 peak search  → initial #mu, A
    //----------------------------------------------------------------
    TH1* h = static_cast<TH1*>(o);
    const double piFitLo = 0.05, piFitHi = 0.35;            // #pi0 window
    const int    iLoPi   = binAt(h,piFitLo),  iHiPi = binAt(h,piFitHi);

    int iMaxPi = iLoPi;
    double maxCntPi = 0;
    for(int i=iLoPi;i<=iHiPi;++i)
      if(h->GetBinContent(i)>maxCntPi){ maxCntPi=h->GetBinContent(i); iMaxPi=i; }

    const double piMu0     = h->GetBinCenter(iMaxPi);       // ~peak
    const double piAmp0    = maxCntPi;
    const double piSigma0  = 0.025;

    //----------------------------------------------------------------
    // 2.  Composite fit function  #pi0‑Gaus + poly‑2 background
    //----------------------------------------------------------------
    TF1 total("total","gaus(0)+pol2(3)",piFitLo,piFitHi);
    total.SetParNames("A","mu","sigma","c0","c1","c2");

    total.SetParameters(piAmp0, piMu0, 0.022,    // narrower #sigma start
                        1,       0,     0);      // flat background
    total.SetParLimits(0,   0,  1e9);            // A  ≥ 0
    total.SetParLimits(1,   0.10,  0.17);        // #mu  within window
    total.SetParLimits(2,   0.010, 0.060);       // #sigma  sensible range

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(3'000);

    /* --- bootstrap background ---------------------------------- */
    TF1 polyTmp("polyTmp","pol2",piFitLo,piFitHi);
    for(int ip = binAt(h,0.11); ip <= binAt(h,0.16); ++ip)
          h->SetBinError(ip, 1e9);            // exclude peak from χ²
    h->Fit(&polyTmp,"QRN0");
    for(int ip = binAt(h,0.11); ip <= binAt(h,0.16); ++ip)
          h->SetBinError(ip, std::sqrt(h->GetBinContent(ip)));

    total.SetParameters(piAmp0, piMu0, 0.022,
                        polyTmp.GetParameter(0),
                        polyTmp.GetParameter(1),
                        polyTmp.GetParameter(2));

    bool fitOK = (h->Fit(&total,"QRN0") == 0);

    /* -------------- #pi0 parameters -------------------------------- */
    double piMu=piMu0, piSig=piSigma0,
           piMuErr=0,   piSigErr=0;
    if(fitOK){
        piMu     = total.GetParameter(1);
        piSig    = total.GetParameter(2);
        piMuErr  = total.GetParError (1);
        piSigErr = total.GetParError (2);
    }

     /* ----------------------------------------------------------------
       * 3.  Robust η‑peak search & fit   (0.45 – 0.80 GeV)
       *     –  stage‑1:   find seed   (side‑band‑subtracted maximum)
       *     –  stage‑2:   background from side‑bands only (2nd‑order)
       *     –  stage‑3:   simultaneous   Gaus + Poly2 fit
      * ---------------------------------------------------------------- */
     double etaMu = 0, etaSig = 0, etaMuErr = 0, etaSigErr = 0;

     if (runID == "Combined")          // ← skip per‑run files
      {
          const double etaLo = 0.45, etaHi = 0.80;     // full η window
          const double side   = 0.03;                  // half‑width of blind zone
                                                       // around the candidate peak

          /* ---- (1) quick seed from raw spectrum ------------------------ */
          int iLo = binAt(h, etaLo),  iHi = binAt(h, etaHi);
          int iMax = iLo;  double maxCnt = 0.;
          for (int i = iLo; i <= iHi; ++i)
              if (h->GetBinContent(i) > maxCnt) { maxCnt = h->GetBinContent(i); iMax = i; }

          if (maxCnt > 0)                       // peak candidate found
          {
              etaMu  = h->GetBinCenter(iMax);
              etaSig = 0.040;                   // start slightly wider than π0

              /* ---- (2) background from *side‑bands only* ----------------- *
               *      exclude ±side around the seed to avoid bias            */
              TF1 bkg("bkg","pol2", etaLo, etaHi);

              // build an exclusion mask for χ² evaluation
              for (int i = iLo; i <= iHi; ++i) {
                  const double x = h->GetBinCenter(i);
                  const bool inPeak = (std::fabs(x - etaMu) < side);
                  h->SetBinError(i, inPeak ? 1e9 : std::sqrt(h->GetBinContent(i)));
              }
              h->Fit(&bkg, "QN0");              // quiet, no draw, store params

              // restore original errors
              for (int i = iLo; i <= iHi; ++i)
                  h->SetBinError(i, std::sqrt(h->GetBinContent(i)));

              /* ---- (3) combined fit – let background float, *
               *      but keep it close to the side‑band shape  */
              TF1 gEta("gEta", "gaus(0)+pol2(3)", etaLo, etaHi);
              gEta.SetParNames("A","#mu","#sigma","c0","c1","c2");
              gEta.SetParameters(maxCnt, etaMu, etaSig,
                                 bkg.GetParameter(0),
                                 bkg.GetParameter(1),
                                 bkg.GetParameter(2));

              // sensible limits
              gEta.SetParLimits(0,      0,      1e9);
              gEta.SetParLimits(1,   etaLo,   etaHi);
              gEta.SetParLimits(2,   0.020,   0.090);

              /* robust limits even if the nominal coefficient is ~0
               * give each parameter at least ±1 × 10⁻³ head‑room            */
              for (int ip = 3; ip <= 5; ++ip) {
                  const double p  = bkg.GetParameter(ip - 3);
                  const double dp = std::max(std::fabs(p)*0.20, 1e-3);
                  gEta.SetParLimits(ip, p - dp, p + dp);
              }
              
              const bool etaOK = (h->Fit(&gEta, "QRN0") == 0);   // Q:quiet R:range N:no‑store 0:draw suppressed
              if (etaOK) {
                  etaMu     = gEta.GetParameter(1);
                  etaMuErr  = gEta.GetParError (1);
                  etaSig    = gEta.GetParameter(2);
                  etaSigErr = gEta.GetParError (2);

                  /* cache for overview */
                  if (_storedEtaFit.count(slice) == 0)
                      _storedEtaFit[slice].reset(new TF1(gEta));
              } else {                          // fit failed → blanks
                  etaMu = etaSig = etaMuErr = etaSigErr = 0;
              }
          }
      }


    //----------------------------------------------------------------
    // 4.  Background TF1  (clone of the poly part)
    //----------------------------------------------------------------
    TF1 poly("bg","pol2",piFitLo,piFitHi);
    poly.SetParameters(total.GetParameter(3),
                       total.GetParameter(4),
                       total.GetParameter(5));
    poly.SetLineColor(kAzure+2); poly.SetLineWidth(2); poly.SetLineStyle(2);

    //----------------------------------------------------------------
    // 5.  Signal / Background CSV  (unchanged for #pi0)
    //----------------------------------------------------------------
    const std::vector<double> ws = {1.25,1.5,1.75,2.0,2.25};
    for(double w : ws){
      const int i1 = binAt(h, std::max(piMu-w*piSig, piFitLo));
      const int i2 = binAt(h, std::min(piMu+w*piSig, piFitHi));
      double S=0,B=0,sErr=0,bErr=0;
      for(int i=i1;i<=i2;++i){
        const double x  = h->GetBinCenter(i);
        const double bg = std::max(poly.Eval(x),0.);
        const double cnt= h->GetBinContent(i);
        B+=bg; S+=cnt-bg; sErr+=cnt; bErr+=bg;
      }
      bErr = std::sqrt(bErr); sErr = std::sqrt(sErr);
      double ratio=(B>0)?S/B:0, rErr=0;
      if(ratio>0) rErr=ratio*std::sqrt((sErr*sErr)/(S*S)+(bErr*bErr)/(B*B));
      csvSB << trig << ',' << ck.E << ',' << ck.chi << ',' << ck.asy << ','
            << slice << ',' << w << ',' << ratio << ',' << rErr << '\n';
    }

    //----------------------------------------------------------------
    // 6.  Pretty plot of this single spectrum -----------------------
    //----------------------------------------------------------------
      {
        TCanvas c; h->SetStats(0); h->Draw();
        poly.Draw("SAME"); total.Draw("SAME");
        if(_storedEtaFit.count(slice)) _storedEtaFit[slice]->Draw("SAME");

        /* ── one‑line cut‑combination label ─────────────────────────────── */
        {
            TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.038);  tl.SetTextAlign(13);
            tl.DrawLatex(0.14, 0.89,
                         Form("E > %.2f GeV   Asym #leq %.2f   #chi^{2} #leq %.2f",
                              ck.E, ck.asy, ck.chi));
        }

        TLegend leg(0.55,0.64,0.88,0.88);
        leg.SetBorderSize(0);
        leg.SetTextAlign(12);

        /* #pi0: write #mu‑line and #sigma‑line underneath one another */
        leg.AddEntry((TObject*)nullptr,
                     Form("#pi^{0}:  #mu = %.3f #pm %.3f GeV", piMu, piMuErr),
                     "");
        leg.AddEntry((TObject*)nullptr,
                     Form("          #sigma = %.3f #pm %.3f GeV", piSig, piSigErr),
                     "");

        /* η: do the same, but only if a peak was found */
        if (etaMu > 0) {
            leg.AddEntry((TObject*)nullptr,
                         Form("#eta:    #mu = %.3f #pm %.3f GeV", etaMu, etaMuErr),
                         "");
            leg.AddEntry((TObject*)nullptr,
                         Form("          #sigma = %.3f #pm %.3f GeV", etaSig, etaSigErr),
                         "");
        }

        leg.Draw();
        c.SaveAs(outPng.string().c_str());
      }

    //----------------------------------------------------------------
    // 7.  CSV & stored‑fit bookkeeping ------------------------------
    //----------------------------------------------------------------
    csvFit << trig << ',' << ck.E << ',' << ck.chi << ',' << ck.asy << ','
           << ck.pLo << ',' << ck.pHi << ','
           << piMu   << ',' << piMuErr  << ','
           << piSig  << ',' << piSigErr << ','
           << etaMu  << ',' << etaMuErr << ','
           << etaSig << ',' << etaSigErr << '\n';

    _fitSummary.emplace(n, FitInfo{slice,ck.pLo,ck.pHi,piMu,piSig,
                                   total.GetChisquare(),total.GetNDF()});

    /* --- keep a copy of the #pi0 fit (once per slice) --------------- */
    if(fitOK && _storedFit.count(slice)==0){
        _storedFit[slice].total.reset(new TF1(total));
        _storedFit[slice].poly .reset(new TF1(poly ));
    }

   if (fitOK) {                                 // ← pTInt condition removed
            auto &m = s_runPoints[slice];          // creates slice entry if missing
            if (m.count(runID) == 0)               // keep first occurrence only
                m[runID] = {piMu, piMuErr, piSig, piSigErr};
    }

    //----------------------------------------------------------------
    // 8.  Cache for overview canvases
    //----------------------------------------------------------------
    cacheForOverview(slice, n, h, pTInt);
    return true;
  }

 private:
  //---------------- helper: put clone into the right cache ----------
  void cacheForOverview(const std::string& slice,const std::string&,
                        TH1* src,bool pTIntegrated)
  {
    auto* cl = static_cast<TH1*>(src->Clone());
    cl->SetDirectory(nullptr);

    if (pTIntegrated)                                    // centrality‑summary
    {
        _centralHists.emplace(slice, cl);
    }
    else                                                  // pT‑binned summary
    {
        _ptHists[slice].push_back(cl);                      // ← new cache
    }
  }

  //---------------- write 2 × 3 overview + #mu,#sigma vs centrality ---------
  void writeSummaryPanels()
  {
    if(_centralHists.empty()) return;

    TCanvas cGrid("c_pi0Cent","#pi0 – all centralities",1800,1000);
    gStyle->SetOptTitle(0);
    cGrid.SetTopMargin(0.12);          // leave 12 % canvas height for the banner
    cGrid.Divide(3,2,0.01,0.01);

    std::vector<double> vC,vCerr,vMu,vMuErr,vSi,vSiErr;
    const double fitLo=0.05, fitHi=0.35;
    int pad=1;

    for(const auto& sl : slices){
      auto it=_centralHists.find(sl);
      if(it==_centralHists.end()) continue;
      TH1* h=it->second;
      cGrid.cd(pad++); h->SetStats(0);

      /* ---- pick stored fit or (rare) fall‑back quick fit ---------- */
      TF1 *fTot=nullptr,*fBg=nullptr;
      std::unique_ptr<TF1> tmpTot,tmpBg;
      if(_storedFit.count(sl)){
          fTot=_storedFit[sl].total.get();
          fBg =_storedFit[sl].poly .get();
      }else{
          tmpTot=std::make_unique<TF1>("fTmp","gaus(0)+pol2(3)",fitLo,fitHi);
          int iMax=h->GetMaximumBin();
          tmpTot->SetParameters(h->GetBinContent(iMax),
                                h->GetBinCenter (iMax),0.02,1,0,0);
          h->Fit(tmpTot.get(),"QRN0");
          tmpBg=std::make_unique<TF1>("fBgTmp","pol2",fitLo,fitHi);
          tmpBg->SetParameters(tmpTot->GetParameter(3),
                               tmpTot->GetParameter(4),
                               tmpTot->GetParameter(5));
          fTot=tmpTot.get(); fBg=tmpBg.get();
      }

      fTot->SetLineColor(kRed+1);   fTot->SetLineWidth(2);
      fBg ->SetLineColor(kBlue+2);  fBg ->SetLineWidth(2);
      fBg ->SetLineStyle(2);

      h->SetMaximum(1.15*h->GetMaximum());
      h->Draw(); fBg->Draw("SAME"); fTot->Draw("SAME");
      if(_storedEtaFit.count(sl)) _storedEtaFit[sl]->Draw("SAME");

      double mu=fTot->GetParameter(1), emu=fTot->GetParError(1);
      double si=fTot->GetParameter(2), esi=fTot->GetParError(2);

      /* ---- helper: ASCII‑only label (defined once per call) --------- */
      static const auto centLabel = [](const std::string& slice) {
            if (slice == "Inclusive") return std::string("Inclusive");
            const auto p = slice.find('_');
            const std::string lo = slice.substr(0, p);
            const std::string hi = slice.substr(p + 1);
            return "Centrality: " + lo + "-" + hi + " %";
      };

      const std::string lbl = centLabel(sl);

      /* ---------- very‑simple “top for first two, bottom for all others” logic ------------- */
      const double lm = gPad->GetLeftMargin();
      const double rm = gPad->GetRightMargin();
      const double tm = gPad->GetTopMargin();
      const double bm = gPad->GetBottomMargin();

      /*  Which slice is this?  Examples of sl: “0_10”, “10_20”, “20_30”, … */
      const bool isTopSlice = (sl == "30_40" || sl == "40_50" || sl == "50_60");

      /* 1) horizontal side – keep right‑hand corner but move ~7 % pad‑width left  */
      const double xText = 1.0 - rm - 0.42;    // shift left a touch

      /* 2) vertical anchor  (unchanged logic)                                     */
      const bool   putBottom = !isTopSlice;
      const double yAnchor   = putBottom ? bm + 0.07 : 1.0 - tm - 0.05;
      const double dy        = 0.063;

      /* 3) prettify trigger and run strings                                       */
      std::string runShort  = runID;                                // drop leading zeros
      if (std::all_of(runID.begin(), runID.end(), ::isdigit))
            runShort = std::to_string(std::stoi(runID));

      std::string trigLabel = trig;                                 // nicer trigger text
      if (trig == "MBD_NandS_geq_2") trigLabel = "MBD NS #geq 2";

      TLatex tx;  tx.SetNDC();  tx.SetTextSize(0.035);  tx.SetTextAlign(13);

      if (runID != "Combined")            // per‑run files: 5‑line block
      {
            const double yRun  = putBottom ? yAnchor + 4*dy : yAnchor;
            const double yTrig = putBottom ? yAnchor + 3*dy : yAnchor -   dy;
            const double yCent = putBottom ? yAnchor + 2*dy : yAnchor - 2*dy;
            const double yMu   = putBottom ? yAnchor +   dy : yAnchor - 3*dy;
            const double ySig  = putBottom ? yAnchor         : yAnchor - 4*dy;

            tx.DrawLatex(xText, yTrig, Form("Trigger: %s", trigLabel.c_str()));
            tx.DrawLatex(xText, yCent, lbl.c_str());                       // centrality
            tx.DrawLatex(xText, yMu,   Form("#mu = %.3f #pm %.3f GeV",  mu,  emu));
            tx.DrawLatex(xText, ySig,  Form("#sigma = %.3f #pm %.3f GeV", si, esi));
      }
      else                                 // combined file: keep original 3‑line block
      {
            const double y1 = putBottom ? yAnchor + 2*dy : yAnchor;
            const double y2 = putBottom ? yAnchor +   dy : yAnchor - dy;
            const double y3 = putBottom ? yAnchor         : yAnchor - 2*dy;

            tx.DrawLatex(xText, y1, lbl.c_str());
            tx.DrawLatex(xText, y2, Form("#mu = %.3f #pm %.3f GeV",  mu,  emu));
            tx.DrawLatex(xText, y3, Form("#sigma = %.3f #pm %.3f GeV", si, esi));
      }

      if(sl!="Inclusive"){
          int lo=std::stoi(sl.substr(0,sl.find('_')));
          int hi=std::stoi(sl.substr(sl.find('_')+1));
          vC.push_back(0.5*(lo+hi));   vCerr.push_back(0.5*(hi-lo));
          vMu.push_back(mu); vMuErr.push_back(emu);
          vSi.push_back(si); vSiErr.push_back(esi);
      }
    }
      fs::path pngGrid = root/"EMCal"/"invMassQA"/cutTag
                            /"Pi0Mass_AllCentrality.png";
      ensure_dir(pngGrid.parent_path());

    {   /* run‑number and cut combination header */
          std::string runShort = runID;
          if (std::all_of(runID.begin(), runID.end(), ::isdigit))
              runShort = std::to_string(std::stoi(runID));

          double eCut = 0, chiCut = 0, asyCut = 0;
          std::smatch m;
          if (std::regex_match(cutTag, m,
              std::regex(R"(E([0-9]+p[0-9]+)_Chi([0-9]+p[0-9]+)_Asym([0-9]+p[0-9]+))")))
          {
              auto p2d = [](const std::string& s){
                  return std::stod(std::regex_replace(s, std::regex("p"), "."));
              };
              eCut   = p2d(m[1]);   chiCut = p2d(m[2]);   asyCut = p2d(m[3]);
          }

          cGrid.cd();                                   // make *canvas* current
          TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035); tl.SetTextAlign(13);
          tl.DrawLatex(0.31, 0.99,
              Form("Run: %s   E > %.2f GeV   Asym #leq %.2f   #chi^{2} #leq %.2f",
                   runShort.c_str(), eCut, asyCut, chiCut));
    }

    cGrid.SaveAs(pngGrid.string().c_str());

    /* ---------- (B) #mu,#sigma versus centrality (#pi0 only, unchanged) --- */
    if(!vC.empty()){
      int n=vC.size();
      auto gMu = std::make_unique<TGraphErrors>(n,
                          vC.data(), vMu.data(), nullptr, vMuErr.data());   // no x‑errors
      auto gSi = std::make_unique<TGraphErrors>(n,
                          vC.data(), vSi.data(), nullptr, vSiErr.data());   // no x‑errors
      gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
      gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

      TCanvas cGS("c_mu_sigma_vs_cent",
                  "#pi^{0} peak position / width vs centrality",800,800);

      /* ---------- pad geometry: add a 2 % blank strip between μ‑ and σ‑panels ---- */
      const double padLeft  = 0.18;          // identical inner widths
      const double padRight = 0.04;          // symmetrical right margin

      const double gapFrac  = 0.02;          // 2 % of canvas → visual spacer
      const double fracBot  = 0.30;          // σ‑panel height (30 %)
      const double fracTop  = 1.0 - fracBot - gapFrac;

      /* ----------------------------- upper (μ) pad ------------------------------ */
      TPad *p1 = new TPad("p1", "upper",
                            0,                    gapFrac + fracBot,   // y‑low
                            1,                    1);                  // y‑high
      p1->SetBottomMargin(0.04);               // 4 % → small white gap
      p1->SetTopMargin   (0.04);
      p1->SetLeftMargin  (padLeft);
      p1->SetRightMargin (padRight);
      p1->Draw();
      p1->cd();

      gMu->SetTitle("; ;#mu_{#pi^{0}} (GeV/c^{2})");
      gMu->Draw("AP");                   // same x‑range for both pads
      gMu->GetXaxis()->SetLabelOffset(999);  // hide x‑labels & ticks in upper pad
      gMu->GetXaxis()->SetTitleOffset(999);
      gMu->GetXaxis()->SetTickLength(0);


      cGS.cd();

      /* ----------------------------- lower (σ) pad ------------------------------ */
      TPad *p2 = new TPad("p2", "lower",
                            0,                    0,                   // y‑low
                            1,                    fracBot);            // y‑high
      p2->SetTopMargin   (0.06);               // 4 % + 2 % gap visual balance
      p2->SetBottomMargin(0.38);
      p2->SetLeftMargin  (padLeft);
      p2->SetRightMargin (padRight);
      p2->Draw();
      p2->cd();

        gSi->SetTitle(";Centrality [%];#sigma_{#pi^{0}} (GeV/c^{2})");
        gSi->Draw("AP");

        /* ---------- axis fonts & ticks – scaled for the small pad ----------------- */
        gSi->GetXaxis()->SetNdivisions(506);
        gSi->GetXaxis()->SetTitleSize(0.09);
        gSi->GetXaxis()->SetLabelSize(0.07);

        gSi->GetYaxis()->SetTitleSize(0.09);
        gSi->GetYaxis()->SetLabelSize(0.07);
        gSi->GetYaxis()->SetTitleOffset(0.90);
        gSi->GetYaxis()->SetTickLength(0.035);

        /* ---------- run‑number + cut‑combination label ---------------------------- */
        {
            /* ---- (1) shorten run ID (strip leading zeros) ---- */
            std::string runShort = runID;
            if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                runShort = std::to_string(std::stoi(runID));

            /* ---- (2) parse the three cut values back from cutTag ---- */
            double eCut = 0, chiCut = 0, asyCut = 0;
            std::smatch m;
            if (std::regex_match(cutTag, m,
                    std::regex(R"(E([0-9]+p[0-9]+)_Chi([0-9]+p[0-9]+)_Asym([0-9]+p[0-9]+))"))) {
                auto p2d = [](const std::string& s){
                    std::string t = std::regex_replace(s, std::regex("p"), ".");
                    return std::stod(t);
                };
                eCut   = p2d(m[1].str());
                chiCut = p2d(m[2].str());
                asyCut = p2d(m[3].str());
            }

            /* ---- (3) draw label ---- */
            cGS.cd();                                         
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.025);

            tl.DrawLatex(0.4, 0.96, Form("Run:  %s", runShort.c_str()));
            tl.DrawLatex(0.4, 0.91,
                         Form("Cuts:  E > %.2f GeV   Asym #leq %.2f   #chi^{2} #leq %.2f",
                              eCut, asyCut, chiCut));
        }
                
        fs::path pngGraph = root/"EMCal"/"invMassQA"/cutTag
                             /"Pi0Mass_Sigma_vs_Centrality.png";
        ensure_dir(pngGraph.parent_path());          // ← create folder

        cGS.SaveAs(pngGraph.string().c_str());
        
        /* ---------- (C1) 2×3 grid of first six pT‑bin spectra (one canvas per centrality) --- */
        for (const auto& sl : slices)
        {
            if (sl == "Inclusive") continue;                       // skip inclusive slice
            auto itH = _ptHists.find(sl);
            if (itH == _ptHists.end() || itH->second.empty()) continue;

            const int nShow = std::min<int>(6, itH->second.size()); // show at most six bins
            TCanvas cGridPt(Form("c_pi0_ptGrid_%s", sl.c_str()),
                            Form("#pi^{0} invariant mass – Cent %s %% (first six p_{T} bins)", sl.c_str()),
                            1800, 1000);
            cGridPt.Divide(3, 2, 0.01, 0.01);

            for (int i = 0; i < nShow; ++i)
            {
                cGridPt.cd(i + 1);
                TH1* hPt = itH->second[i];
                hPt->SetStats(0);
                hPt->Draw();
            }

            fs::path pngGridPt = root / "EMCal" / "invMassQA" / cutTag / sl
                                 / "pTsummarizedMeanSigmaDistributions"
                                 / "Pi0Mass_First6pTbins.png";
            ensure_dir(pngGridPt.parent_path());
            cGridPt.SaveAs(pngGridPt.string().c_str());
        }

        /* ---------- (C2)  μ,σ  versus  pT   (one canvas per centrality) ---------- */
        for (const auto& sl : slices)
        {
            if (sl == "Inclusive") continue;
            auto it = _ptHists.find(sl);
            if (it == _ptHists.end()) continue;

            /* ---- collect μ,σ and errors ordered by pT centre ---- */
            std::map<double, FitInfo> byPt;
            for (const auto& [hName, fi] : _fitSummary)
                if (fi.slice == sl && fi.pLo >= 0 && fi.pHi >= 0)
                    byPt[0.5 * (fi.pLo + fi.pHi)] = fi;

            if (byPt.size() < 2) continue;

            const int n = byPt.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            int k = 0;
            for (const auto& [pt, fi] : byPt)
            {
                x[k]   = pt;
                yMu[k] = fi.mean;    eMu[k] = 0;
                ySi[k] = fi.sigma;   eSi[k] = 0;
                ++k;
            }

            auto gMu = std::make_unique<TGraphErrors>(n, x.data(), yMu.data(),
                                                      nullptr, eMu.data());
            auto gSi = std::make_unique<TGraphErrors>(n, x.data(), ySi.data(),
                                                      nullptr, eSi.data());
            gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
            gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

            TCanvas cPT(Form("c_mu_sigma_vs_pt_%s", sl.c_str()),
                        "#pi^{0} peak position / width vs p_{T}", 800, 800);

            const double padLeft = 0.18, padRight = 0.04, gap = 0.02, fracBot = 0.30;

            TPad* p1 = new TPad("p1", "", 0, gap + fracBot, 1, 1);
            p1->SetBottomMargin(0.04); p1->SetTopMargin(0.04);
            p1->SetLeftMargin(padLeft); p1->SetRightMargin(padRight);
            p1->Draw(); p1->cd();
            gMu->SetTitle("; ;m_{#pi^{0}}  (GeV)"); gMu->Draw("AP");
            gMu->GetXaxis()->SetLabelOffset(999); gMu->GetXaxis()->SetTitleOffset(999);

            cPT.cd();
            TPad* p2 = new TPad("p2", "", 0, 0, 1, fracBot);
            p2->SetTopMargin(0.06); p2->SetBottomMargin(0.38);
            p2->SetLeftMargin(padLeft); p2->SetRightMargin(padRight);
            p2->Draw(); p2->cd();
            gSi->SetTitle(";p_{T}  [GeV/#it{c}];#sigma_{#pi^{0}}  (GeV)");
            gSi->Draw("AP");
            gSi->GetXaxis()->SetNdivisions(506);
            gSi->GetXaxis()->SetTitleSize(0.09); gSi->GetXaxis()->SetLabelSize(0.07);
            gSi->GetYaxis()->SetTitleSize(0.09); gSi->GetYaxis()->SetLabelSize(0.07);
            gSi->GetYaxis()->SetTitleOffset(0.90); gSi->GetYaxis()->SetTickLength(0.035);

            fs::path pngPT = root / "EMCal" / "invMassQA" / cutTag / sl
                              / "pTsummarizedMeanSigmaDistributions"
                              / "Pi0Mass_Sigma_vs_pT.png";
            ensure_dir(pngPT.parent_path());
            cPT.SaveAs(pngPT.string().c_str());
        }
    }
  }

    //--------------------------------------------------------------------
    //  Write run‑by‑run π0‐mass summary   (called once, after “Combined”)
    //--------------------------------------------------------------------
  void writeRunSummary()
  {
        /* -------------------------------------------------------------
         * 1)  Guard clauses
         * ----------------------------------------------------------- */
        if (runID != "Combined") {
            std::cout << "[Pi0QA]  writeRunSummary() – skipped:  not in \"Combined\" pass\n";
            return;
        }
        if (s_runPoints.empty()) {
            std::cout << "[Pi0QA]  writeRunSummary() – WARNING:  s_runPoints is EMPTY\n";
            return;
        }

        // one PNG per (E,χ,asy) cut – avoid duplicates when several Pi0QA
        static std::unordered_set<std::string> s_done;
        if (s_done.count(cutTag)) {
            std::cout << "[Pi0QA]  writeRunSummary() – already written for cut "
                      << cutTag << '\n';
            return;
        }
        s_done.insert(cutTag);

        std::cout << "[Pi0QA]  ▶ building run‑summary for cut " << cutTag << '\n';
        std::cout << "          slices stored : " << s_runPoints.size() << '\n';

        /* -------------------------------------------------------------
         * 2)  Colour palette
         * ----------------------------------------------------------- */
        const int cols[] = {kBlue+1,kRed+1,kGreen+2,kMagenta+2,
                            kOrange+1,kCyan+2,kSpring+5,kPink+1};
        constexpr int nCols = sizeof(cols)/sizeof(int);

        /* -------------------------------------------------------------
         * 3)  Build one TGraphErrors per slice
         * ----------------------------------------------------------- */
        std::vector<TGraphErrors*> gMuList, gSiList;
        TLegend                    leg(0.12,0.73,0.42,0.88);  leg.SetBorderSize(0);

        int colourIdx = 0;
        for (const auto& [slice, mp] : s_runPoints)
        {
            /* --- collect run numbers (digits only) ------------------ */
            std::vector<int> runs;
            for (const auto& [runStr,_] : mp)
                if (std::all_of(runStr.begin(), runStr.end(), ::isdigit))
                    runs.push_back(std::stoi(runStr));

            if (runs.empty()) {
                std::cout << "          · slice \"" << slice << "\" skipped – no numeric runs\n";
                continue;
            }
            std::sort(runs.begin(), runs.end());

            /* --- fill arrays ---------------------------------------- */
            const int n = runs.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            for (int i = 0; i < n; ++i) {
                const auto& p = mp.at(std::to_string(runs[i]));
                x[i]  = runs[i];
                yMu[i]= p.mu;      eMu[i]= p.muErr;
                ySi[i]= p.sigma;   eSi[i]= p.sigmaErr;
            }

            auto gMu = new TGraphErrors(n, x.data(), yMu.data(), nullptr, eMu.data());
            auto gSi = new TGraphErrors(n, x.data(), ySi.data(), nullptr, eSi.data());

            const int col = cols[colourIdx++ % nCols];
            gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
            gSi->SetMarkerStyle(kFullCircle); gSi->SetLineWidth(2);
            gMu->SetMarkerColor(col); gMu->SetLineColor(col);
            gSi->SetMarkerColor(col); gSi->SetLineColor(col);

            gMuList.push_back(gMu);
            gSiList.push_back(gSi);

            const std::string lbl = (slice == "Inclusive") ?
                                     "Inclusive" : "Cent " + slice + " %";
            leg.AddEntry(gMu, lbl.c_str(), "pl");

            std::cout << "          · slice \"" << slice << "\" – "
                      << n << " runs, colour idx " << (colourIdx-1) << '\n';
        }

        if (gMuList.empty()) {
            std::cout << "[Pi0QA]  writeRunSummary() – ABORT:  all slices empty, no graph created\n";
            return;
        }

        /* -------------------------------------------------------------
         * 4)  Draw canvas (μ upper, σ lower)
         * ----------------------------------------------------------- */
        TCanvas cR("c_mu_sigma_vs_run_allCent",
                   "#pi^{0} peak position / width vs run", 900, 800);

        // ---- μ pad --------------------------------------------------
        TPad *p1 = new TPad("p1","", 0, 0.35, 1, 1);
        p1->SetBottomMargin(0.02);  p1->Draw();  p1->cd();
        gMuList.front()->SetTitle(";Run number;m_{#pi^{0}}  (GeV)");

        for (std::size_t i = 0; i < gMuList.size(); ++i)
            gMuList[i]->Draw(i == 0 ? "AP" : "P SAME");
        leg.Draw();

        // ---- σ pad --------------------------------------------------
        cR.cd();
        TPad *p2 = new TPad("p2","", 0, 0, 1, 0.32);
        p2->SetTopMargin(0.02);  p2->SetBottomMargin(0.30);
        p2->Draw();  p2->cd();

        gSiList.front()->SetTitle(";Run number;#sigma_{#pi^{0}}  (GeV)");
        for (std::size_t i = 0; i < gSiList.size(); ++i)
            gSiList[i]->Draw(i == 0 ? "AP" : "P SAME");

        /* -------------------------------------------------------------
         * 5)  Save PNG – catch I/O problems
         * ----------------------------------------------------------- */
        fs::path pngRun = root/"summary/pi0"/cutTag
                       /"CentralitySummaryHistograms"/"Pi0Mass_Sigma_vs_Run_AllCentrality.png";

        ensure_dir(pngRun.parent_path());

        try {
            cR.SaveAs(pngRun.string().c_str());
            std::cout << "[Pi0QA]  ✔  run‑summary written to "
                      << pngRun.string() << '\n';
        }
        catch (const std::exception& ex) {
            std::cerr << "[Pi0QA]  ERROR while saving \"" << pngRun.string()
                      << "\" – " << ex.what() << '\n';
        }
  }
  /* ---------------------------------------------------------------- */
  /*  data members                                                    */
  /* ---------------------------------------------------------------- */
  std::ofstream&                                   csvFit;
  std::ofstream                                    csvSB;

  std::unordered_map<std::string, TH1*>            _centralHists;

  std::unordered_map<std::string, FitPair>         _storedFit;
  std::unordered_map<std::string, std::unique_ptr<TF1>> _storedEtaFit;

  std::unordered_map<std::string,
                     std::vector<TH1*>>            _ptHists;
  std::unordered_map<std::string, FitInfo>         _fitSummary;

  /* run label of this instance */
  std::string runID;
  /* directory tag that identifies one (E , χ² , asym) cut‑combination */
  std::string cutTag;

  /* ---------- static: accumulate #pi0 points per slice & run ---------- */
  struct RunPoint { double mu, muErr, sigma, sigmaErr; };
  /*   slice → runID → RunPoint   */
  static inline std::unordered_map<
                std::string,
                std::unordered_map<std::string,RunPoint>>  s_runPoints;
  static inline bool s_summaryWritten = false;
};

// ─── Minimal fall‑back style helpers ────────────────────────────────────
inline void tidyAxes(TH1* h)
{
    if (!h) return;
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
}

inline void styleAxes(TH1* h, bool forSummary = false)
{
    if (!h) return;
    h->SetLineWidth(2);
    h->SetMarkerSize(1.0);
    if (forSummary) h->SetMarkerColor(kRed);
}

inline void setupPad(TVirtualPad* p)
{
    if (!p) return;
    p->SetLeftMargin (0.12);
    p->SetRightMargin(0.18);
    p->SetBottomMargin(0.12);
    p->SetTopMargin  (0.08);
}


/* ==================================================================== *
 * §‑1  Helpers added for run‑label and “auto‑tight” axis scaling       *
 * ==================================================================== */

/* Strip leading zeroes from a run‑directory name like “00044777” → “44777”.
   If the directory contains non‑digits (e.g. “Run_00044777”), the numeric
   suffix is preserved but zeroes are still removed.                     */
static std::string stripLeadingZeros(const std::string& runDir)
{
    std::smatch m;
    if (std::regex_search(runDir, m, std::regex(R"((\d+)$)")))
        return std::to_string(std::stoul(m[1].str()));   // 44777
    return runDir;                                       // fallback
}

static void tightenAxes(TH1* h)
{
    // Fine‑tune only what exists on a 1‑D histogram
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->GetXaxis()->SetTitleOffset(1.1F);
    h->GetYaxis()->SetTitleOffset(1.25F);
}

/* Tighten the displayed X/Y range so that the upper edge coincides with
   the last *non‑empty* bin, keeping the lower edge fixed at zero.       */
static void tightenAxes(TH2* h)
{
    auto* axX = h->GetXaxis();
    auto* axY = h->GetYaxis();

    const int lastX = h->GetNbinsX();
    const int lastY = h->GetNbinsY();

    int hiX = lastX;
    while (hiX > 1 && h->Integral(hiX, lastX, 1, lastY) == 0) --hiX;

    int hiY = lastY;
    while (hiY > 1 && h->Integral(1, lastX, hiY, lastY) == 0) --hiY;

    axX->SetRangeUser(0., axX->GetBinUpEdge(hiX));
    axY->SetRangeUser(0., axY->GetBinUpEdge(hiY));
}

/* Draw run‑number (without leading zeroes) in the upper‑left corner.    */
static void drawRunLabel(const std::string& runID)
{
    TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.035);
    tl.DrawLatex(0.04, 0.94, ("Run " + runID).c_str());
}

/* ──────────────────────────────────────────────────────────────────────────
 *  Correlation QA
 * ──────────────────────────────────────────────────────────────────────── */
struct NSPair { std::shared_ptr<TH2> n, s; };
static inline std::unordered_map<std::string, NSPair> g_nsCache;

class CorrQA : public QA
{
 public:
    /* ---------- static run-summary cache --------------------------- */
    using RunMap  = std::unordered_map<std::string, std::shared_ptr<TH2>>;
    using NameMap = std::unordered_map<std::string, RunMap>;
    
    CorrQA(const std::string& trig,
           const fs::path&    base,
           const CentList&    cent)
      : QA(trig, base, cent)
    {}

    ~CorrQA() override
    {
        writeCentralityOverviews();   // §4 below
        writeRunSummaries();          // implemented just below
    }

 private:
    static std::unordered_map<std::string, NameMap> s_cache;

    /* ------------------------------------------------------------------ *
     *  Produce per-run correlation summary panels.
     *  (Stub-version prevents compile errors – extend as required.)      *
     * ------------------------------------------------------------------ */
    void writeRunSummaries() {}

    // ─────────────────────────────── 1. per-histogram ──────────────────────
    bool process(TObject* o) override
    {
        /* ------------------------------------------------------------------ *
         * Accept both 2‑D (“…_vs_…”) maps  *and*  the new 1‑D Δη / Δφ spectra *
         * ------------------------------------------------------------------ */
        if (!o->InheritsFrom(TH1::Class())) return false;     // neither TH1 nor TH2 → skip

        const std::string hName = o->GetName();

        /* Reject Jet‑QA and event‑plane families – they are handled elsewhere */
        static const std::vector<std::string> skipPrefixes = {
            "h_leadEt_",          // jet leading/subleading E_T maps
            "h_maxJetEt_",        // jet spectra
            "h_jetEt_",           // 3‑D jet histos
            "h_Psi",              // event‑plane Ψₙ histos
            "p_R2_",              // resolution proxies
        };
        for (const auto& p : skipPrefixes)
            if (hName.rfind(p, 0) == 0)
                return false;                          // delegate to the right QA class

        const bool is2D = (hName.find("_vs_") != std::string::npos);
        const bool is1D = (hName.rfind("h_dEta_",0)==0) || (hName.rfind("h_dPhi_",0)==0);
        if (!(is2D || is1D)) return false;

        const std::string slice   = sliceKey(hName);          // Inclusive / 0_10 …
        const bool        hasCent = (slice != "Inclusive");

        /* ------------------------------------------------------------------ *
         * A.  Extract *clean* detector names                                  *
         * ------------------------------------------------------------------ */
        std::string tokA, tokB;
        if (is2D) {
            tokA = hName.substr( 2, hName.find("_vs_") - 2 );
            tokB = hName.substr( hName.find("_vs_") + 4 );
        } else {                           /* 1‑D  “h_dEta_<DET1>_<DET2>” …     */
            const std::string tail = hName.substr(7);          // drop “h_dEta_ / h_dPhi_”
            const std::size_t us   = tail.find('_');
            tokA = tail.substr(0, us);
            tokB = tail.substr(us + 1);
        }

        const std::string detA = canonicalDet(tokA);
        const std::string detB = canonicalDet(tokB);

        /* detect North/South once the canonical name is known --------------- */
        const bool isNorth = (tokA.find("_North")!=std::string::npos ||
                              tokB.find("_North")!=std::string::npos);
        const bool isSouth = (tokA.find("_South")!=std::string::npos ||
                              tokB.find("_South")!=std::string::npos);

        /* ------------------------------------------------------------------ *
         * B.  Directory group name construction                               *
         * ------------------------------------------------------------------ */
        const bool hasIH = (detA=="IHCal" || detB=="IHCal");
        const bool hasOH = (detA=="OHCal" || detB=="OHCal");

        std::string groupDir, hcalMode, otherDet;
        if (hasIH ^ hasOH) {                       /* exactly one HCal sub‑system  */
            hcalMode = hasIH ? "IHCal" : "OHCal";
            otherDet = (detA!=hcalMode)?detA:detB;

            /* put EMCal first to obtain  “EMCal_IHCal”  /  “EMCal_OHCal” */
            if (otherDet == "EMCal")
                groupDir = otherDet + "_" + hcalMode;
            else
                groupDir = hcalMode + "_" + otherDet;
        }
        else if (hasIH && hasOH) {                 /* both IHCal & OHCal present   */
            hcalMode = "totalHCal";
            otherDet = (detA!="IHCal" && detA!="OHCal") ? detA : detB;
            groupDir = hcalMode + "_" + otherDet;
        }
        else {                                     /* no HCal involved             */
            groupDir = (detA < detB) ? detA + "_" + detB
                                     : detB + "_" + detA;
        }

        /* ------------------------------------------------------------------ *
         * C.  Style & save the individual PNG (handles TH1 and TH2)           *
         * ------------------------------------------------------------------ */
        TH2* h2 = dynamic_cast<TH2*>(o);
        TH1* h1 = dynamic_cast<TH1*>(o);           /* works for TH1 *and* TH2      */

        tidyAxes(h1);               /* tidyAxes / styleAxes accept TH1 base class */
        styleAxes(h1,false);

        fs::path outDir = root / "correlations" / groupDir;
        if (hasCent) outDir /= ("Cent_" + slice);
        ensure_dir(outDir);

        fs::path pngFile = outDir / (o->GetName() + std::string(".png"));
        {
            TCanvas c("c_corr","",1100,800); setupPad(&c);

            if (h2) {                               /* 2‑D map */
                c.SetLogz();
                tightenAxes(h2);
                h2->Draw("COLZ");
            } else {                                /* 1‑D spectrum */
                tightenAxes(h1);
                h1->Draw();
            }
            drawRunLabel( stripLeadingZeros(root.parent_path().filename().string()) );
            c.SaveAs(pngFile.string().c_str());
        }

        /* ------------------------------------------------------------------
         * D.  Collect centrality-dependent clones for the overview canvas
         * ----------------------------------------------------------------- */
        if (hasCent)
        {
            const std::string baseKey = groupDir + "|" +
                                        stripCentSuffix(hName) ;  // "h_SEPD_vs_MBD"
            auto& v = m_centCache[baseKey];
            auto  cl = std::shared_ptr<TH2>(static_cast<TH2*>(h2->Clone()));
            cl->SetDirectory(nullptr); tidyAxes(cl.get()); styleAxes(cl.get(),false);
            v.emplace_back(slice, std::move(cl));
        }

        /* ─────────── existing E-F-G blocks stay exactly as before ───────── */
        handleNorthSouthPairing(hName, isNorth, isSouth, hasCent,
                                groupDir, slice, o);
        cacheForRunSummary(groupDir, hName, o);
        aggregateTotalHCal(hasIH, hasOH, otherDet, hName,
                           hasCent, slice, groupDir, o);
        return true;
    }

    /* ==================================================================== *
     * §0  Helper: canonical detector name                                  *
     * ==================================================================== */
    static std::string canonicalDet(const std::string& tok)
    {
        auto has = [&](const char* pat){ return tok.find(pat)!=std::string::npos; };

        if (has("IHCAL") || has("IHCal"))   return "IHCal";
        if (has("OHCAL") || has("OHCal"))   return "OHCal";
        if (has("CEMC"))                   return "EMCal";
        if (has("SEPD") || has("sEPD"))    return "sEPD";
        if (has("MBD"))                    return "MBD";

        /* fall-back: take the first token up to a digit or underscore */
        const std::size_t pos = tok.find_first_of("0123456789_");
        return tok.substr(0, pos);
    }

    /* strip “…_<lo>_<hi>_<trigger>” so all centrality clones map back
       to the *same* base name                                           */
    static std::string stripCentSuffix(const std::string& h)
    {
        const std::regex re("_(\\d{1,3})_(\\d{1,3})_.*$");
        return std::regex_replace(h, re, "");
    }

    /* ==================================================================== *
     * §1  North / South paired canvas (unchanged logic, moved to a helper) *
     * ==================================================================== */
    void handleNorthSouthPairing(const std::string& hName,bool isNorth,bool isSouth,
                                 bool hasCent, [[maybe_unused]] const std::string& groupDir,
                                 const std::string& slice,TObject* o)
    {
        const bool isNScandidate =
            (hName.find("SEPD")!=std::string::npos || hName.find("sEPD")!=std::string::npos) &&
            (isNorth ^ isSouth);

        if (!isNScandidate) return;

        std::string baseKey = std::regex_replace(hName,
                                                 std::regex("(_North|_South)"),
                                                 "");
        const std::string cacheKey = groupDir + "|" + baseKey + "|" + slice;
        auto& pair = g_nsCache[cacheKey];

        auto* cl = static_cast<TH2*>(o->Clone());
        cl->SetDirectory(nullptr);  tidyAxes(cl); styleAxes(cl,false);
        (isSouth ? pair.s : pair.n).reset(cl);

        if (pair.n && pair.s) {                          // both halves ready
            fs::path dir = root / "correlations" / groupDir;
            if (hasCent) dir /= ("Cent_" + slice);
            ensure_dir(dir);

            fs::path png = dir / (baseKey + std::string("_NS.png"));

            const double zMax = std::max(pair.n->GetMaximum(),
                                          pair.s->GetMaximum());
            pair.n->SetMaximum(zMax);  pair.s->SetMaximum(zMax);
            pair.n->SetMinimum(1);     pair.s->SetMinimum(1);

            TCanvas c("c_ns","",1200,600); c.Divide(2,1,0.01,0.01);
            c.cd(1); setupPad(gPad); gPad->SetLogz();
            tightenAxes(pair.s.get());
            pair.s->Draw("COLZ");
            drawRunLabel( stripLeadingZeros(root.parent_path().filename().string()) );

            c.cd(2); setupPad(gPad); gPad->SetLogz();
            tightenAxes(pair.n.get());
            pair.n->Draw("COLZ");
            drawRunLabel( stripLeadingZeros(root.parent_path().filename().string()) );
            c.SaveAs(png.string().c_str());

            g_nsCache.erase(cacheKey);
        }
    }

    /* ==================================================================== *
     * §2  Combined run-summary cache (unchanged, moved to helper)          *
     * ==================================================================== */
    void cacheForRunSummary(const std::string& groupDir,
                            const std::string& hName,
                            TObject* o)
    {
        const std::string runID = root.parent_path().filename().string();
        if (runID == "Combined") return;

        auto* cl = static_cast<TH2*>(o->Clone());
        cl->SetDirectory(nullptr); tidyAxes(cl); styleAxes(cl,true);

        s_cache[groupDir][hName][runID].reset(cl);
    }

    /* ==================================================================== *
     * §3  total-HCal aggregation (unchanged, moved to helper)              *
     * ==================================================================== */
    void aggregateTotalHCal(bool hasIH,bool hasOH,
                            const std::string& otherDet,
                            const std::string& hName,
                            bool hasCent,const std::string& slice,
                            [[maybe_unused]] const std::string& groupDir,
                            TObject* o)
    {
        if (!(hasIH ^ hasOH)) return;                 // both or none → skip

        const std::string totGroup = "totalHCal_" + otherDet;
        const std::string canonName =
            std::regex_replace(hName,
                               std::regex("(IHCAL|IHCal|OHCAL|OHCal)"),
                               "HCal");

        const std::string runID = root.parent_path().filename().string();
        const std::string aggKey = totGroup + "|" + canonName + "|" +
                                   slice + "|" + runID;

        struct Agg { std::shared_ptr<TH2> h; int parts = 0; };
        static std::unordered_map<std::string, Agg> agg;

        Agg& a = agg[aggKey];
        if (!a.h) {
            a.h.reset(static_cast<TH2*>(o->Clone()));
            a.h->SetDirectory(nullptr);
            a.h->SetName(canonName.c_str());
        } else {
            a.h->Add(static_cast<TH2*>(o));
        }
        if (++a.parts != 2) return;                  // wait for the other HCal

        tidyAxes(a.h.get()); styleAxes(a.h.get(), false);

        fs::path dir = root / "correlations" / totGroup;
        if (hasCent) dir /= ("Cent_" + slice);
        ensure_dir(dir);

        fs::path png = dir / (canonName + ".png");
        TCanvas cTot(("c_"+canonName).c_str(),"",1100,800); setupPad(&cTot);
        cTot.SetLogz();
        tightenAxes(a.h.get());
        a.h->Draw("COLZ");
        drawRunLabel( stripLeadingZeros(runID) );
        cTot.SaveAs(png.string().c_str());

        if (runID != "Combined") {
            auto* cl = static_cast<TH2*>(a.h->Clone());
            cl->SetDirectory(nullptr); styleAxes(cl,true);
            s_cache[totGroup][canonName][runID].reset(cl);
        }
    }

    /* ==================================================================== *
     * §4  Overview canvas with all centrality bins                         *
     * ==================================================================== */
    void writeCentralityOverviews()
    {
        for (auto& [key, vec] : m_centCache) {
            if (vec.empty()) continue;

            /* key format:  "groupDir|baseHistName"                          */
            const std::size_t pos = key.find('|');
            const std::string groupDir = key.substr(0, pos);
            const std::string baseHist = key.substr(pos+1);

            /* sort by user-provided centrality order                        */
            std::sort(vec.begin(), vec.end(),
                      [&](auto& a, auto& b)
                      {
                          const auto idx = [&](const std::string& s)
                          {
                              auto it = std::find(slices.begin(),
                                                  slices.end(), s);
                              return (it==slices.end())
                                      ? INT_MAX : std::distance(slices.begin(), it);
                          };
                          return idx(a.first) < idx(b.first);
                      });

            /* canvas geometry:  ≤6 slices → 2×3, else ceil(sqrt(N)) × same  */
            const int n = static_cast<int>(vec.size());
            int nCols = 3, nRows = 2;
            if (n > 6) {
                nCols = static_cast<int>(std::ceil(std::sqrt(n)));
                nRows = static_cast<int>(std::ceil(double(n)/nCols));
            }
            TCanvas c("c_overview","", nCols*550, nRows*500);
            c.Divide(nCols, nRows, 0.001, 0.001);

            for (int i = 0; i < n; ++i) {
                c.cd(i+1);  setupPad(gPad);
                gPad->SetLogz();
                tightenAxes(vec[i].second.get());
                vec[i].second->Draw("COLZ");
                drawRunLabel( stripLeadingZeros(root.parent_path().filename().string()) );
                TLatex tl; tl.SetNDC(); tl.SetTextSize(0.04);

                /* build human‑readable label  “low %  ≤ centrality < high %” */
                std::string label;
                if (vec[i].first == "Inclusive") {
                    label = "Inclusive";
                } else {
                    std::smatch m;
                    if (std::regex_match(vec[i].first, m, std::regex(R"((\d{1,3})_(\d{1,3}))"))) {
                        const int lo = std::stoi(m[1].str());
                        const int hi = std::stoi(m[2].str());
                        std::ostringstream oss;
                        oss << lo << "\\%\\;#leq\\;centrality\\;<\\;" << hi << "\\%";
                        label = oss.str();
                    } else {
                        label = vec[i].first;            // fallback – unexpected slice key
                    }
                }
                tl.DrawLatex(0.05, 0.85, label.c_str());
            }

            fs::path dir = root / "correlations" / groupDir;
            ensure_dir(dir);
            fs::path png = dir / (baseHist + std::string("_CentSummary.png"));
            c.SaveAs(png.string().c_str());
        }
        m_centCache.clear();
    }

    /* ==================================================================== *
     * data members                                                         *
     * ==================================================================== */
    using SliceClone = std::pair<std::string /*slice*/, std::shared_ptr<TH2>>;
    std::unordered_map<std::string, std::vector<SliceClone>> m_centCache;
};

/* static data */
std::unordered_map<std::string, CorrQA::NameMap> CorrQA::s_cache;


class EmcalQA : public QA
{
public:
  using QA::QA;               // forward all ctors from QA

  // ────────────────────────────────────────────────────────────────────
  //  summary panel (called once, after the very last histogram)
  // ────────────────────────────────────────────────────────────────────
  ~EmcalQA() override
  {
    if (_centralMaps.empty()) return;        // nothing to paint

    // 0) choose a perceptually‑uniform palette once for all pads
    gStyle->SetPalette(kBird);               // kViridis, kCool, …
    gStyle->SetNumberContours(100);

    // 1) impose the common Z‑range collected during processing
    if (_zMin < _zMax) {                     // at least one non‑empty map
      for (auto& kv : _centralMaps) {
        kv.second->SetMinimum(_zMin);
        kv.second->SetMaximum(_zMax);
      }
      Int_t   idx0 = gStyle->GetColorPalette(0);      // palette entry → index
      TColor* c0   = gROOT->GetColor(idx0);           // index → TColor*
      if (c0) c0->SetRGB(1.0, 1.0, 1.0);
    }

    // 2) compose the 2 × 3 overview canvas
    TCanvas c("c_emcalCent", "EMCal hit‑maps – all centralities",
              2100, 1200);
    c.Divide(3, 2, 0.01, 0.01);

    int pad = 1;
    for (const auto& kv : _centralMaps)        // iterate over the maps we really have
    {
        if (pad > 6) break;                      // 2 × 3 canvas → max 6 pads
        c.cd(pad++);
        kv.second->Draw("COLZ");

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
        const std::string lbl = (kv.first == "Inclusive")
                                ? "Inclusive"
                                : ("Cent " + kv.first);
        tl.DrawLatex(0.15, 0.85, lbl.c_str());
    }

    fs::path out = root / "EMCal" / "EMCalHitMap_AllCentrality.png";
    ensure_dir(out.parent_path());
    c.SaveAs(out.string().c_str());
  }

  // ────────────────────────────────────────────────────────────────────
  //  P R O C E S S   – called once for every histogram in “EMCal/*”
  // ────────────────────────────────────────────────────────────────────
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;
    const std::string n = o->GetName();
    if (n.rfind("h_EMC_", 0) != 0)          return false;   // not EMCal

    const std::string slice = sliceKey(n);
    const bool isMap        = (n.find("_EtaPhiMap_") != std::string::npos);

    // ── helper that turns a raw η‑φ TH2 into a ready‑to‑paint square PNG
    auto makePanel = [&](TH2* src) -> std::unique_ptr<TH2F>
    {
      /* 0)  constants -------------------------------------------------- */
      constexpr int nPhi = 256;                 // rows (Y)
      constexpr int nEta =  96;                 // cols (X)
      constexpr int px   =   6;                 // pixel‑size in PNG

      /* 1)  clone & bad‑board masking ---------------------------------- */
      std::unique_ptr<TH2> h(static_cast<TH2*>(src->Clone()));
      h->SetDirectory(nullptr);
      h->SetStats(0);
      h->SetContour(100);

      for (int ip = 0; ip < nPhi; ++ip)
        for (int ie = 0; ie < nEta; ++ie)
          if (isBadBoard(sector_from_idx(ie, ip),
                         ib_from_idx(ie, ip)))
            h->SetBinContent(h->FindBin(ip, ie), -9999.);   // masked → white

      /* 2)  rotate (η ↦ X, φ ↦ Y) -------------------------------------- */
      std::unique_ptr<TH2F> rot(new TH2F(("hRot_" + std::string(src->GetName())).c_str(),
                                         h->GetTitle(),
                                         nEta, 0, nEta,
                                         nPhi, 0, nPhi));
      for (int ip = 1; ip <= nPhi; ++ip)
        for (int ie = 1; ie <= nEta; ++ie)
          rot->SetBinContent(ie, ip, h->GetBinContent(ip, ie));

      rot->SetDirectory(nullptr);
      rot->SetMinimum(1.);                    // under‑flow colour = white
      rot->SetTitleOffset(0.9, "X"); rot->SetTitleOffset(1.4, "Y");
      rot->GetXaxis()->SetTitle("Tower #eta");
      rot->GetYaxis()->SetTitle("Tower #phi");
      rot->GetXaxis()->SetNdivisions(12, kFALSE);   // every 8 η
      rot->GetYaxis()->SetNdivisions(32, kFALSE);   // every 8 φ

      /* 2a) GLOBAL Z‑RANGE UPDATE  (ignore masked / empty) ------------- */
      {
        std::vector<double> vv; vv.reserve(nPhi * nEta);
        for (int y = 1; y <= nPhi; ++y)
          for (int x = 1; x <= nEta; ++x) {
            const double v = rot->GetBinContent(x, y);
            if (v > 0.) vv.push_back(v);
          }

        if (!vv.empty()) {
          std::sort(vv.begin(), vv.end());
          const double lo = vv.front();
          const double hi = vv[ static_cast<std::size_t>(vv.size() * kFracSaturation) ];
          _zMin = std::min(_zMin, lo);
          _zMax = std::max(_zMax, hi);
        }
      }
        
      std::string runLabel = root.parent_path().filename().string();   // “00066484” …
      if (std::all_of(runLabel.begin(), runLabel.end(), ::isdigit))    // numeric?
            runLabel = std::to_string(std::stoi(runLabel));              // strip leading 0s

      std::string baseTitle = src->GetTitle();                         // keep original
      if (baseTitle.empty()) baseTitle = src->GetName();               // fallback to name
      rot->SetTitle( (baseTitle + " (" + runLabel + ")").c_str() );

      /* 3)  save individual PNG (1 bin ⇔ 1 pixel) ---------------------- */
      const int cw = (kEMCalCanvasW > 0) ? kEMCalCanvasW : nEta * px;
      const int ch = (kEMCalCanvasH > 0) ? kEMCalCanvasH : nPhi * px;

      fs::path outPng = cPath(root, slice, "EMCal/otherQA")
                        / (src->GetName() + std::string(".png"));
      ensure_dir(outPng.parent_path());

      TCanvas c("c_emcal", "", cw, ch);
      c.SetRightMargin (0.17);
      c.SetFixedAspectRatio(false);           // allow resizing in a viewer
      c.SetLeftMargin  (0.14);
      c.SetBottomMargin(0.08);
      const double kTopPad = 0.055;          // 0.04 → 0.10 moves title upward
      c.SetTopMargin(kTopPad);
      c.SetFixedAspectRatio();                // lock 1 bin = 1 pixel

      rot->Draw("COLZ");

      /* 4)  ultra‑light grid (every 8 towers) -------------------------- */
      TLine l; l.SetLineColor(kBlack); l.SetLineWidth(1);
      for (int x = 0; x <= nEta; x += 8) l.DrawLine(x, 0, x, nPhi);
      for (int y = 0; y <= nPhi; y += 8) l.DrawLine(0, y, nEta, y);

      /* 5)  North / South divider (η = 48) ----------------------------- */
      TLine ns(48, 0, 48, nPhi);
      ns.SetLineColor(kBlack); ns.SetLineWidth(3);
      ns.Draw();

      /* 6)  sector / inner‑board labels (subtle grey) ------------------ */
      TLatex tx;
      tx.SetTextSize(0.020); tx.SetTextAlign(22);
      tx.SetTextColorAlpha(kGray + 2, 0.70);

      for (int s = 0; s < 64; ++s) {               // 32 per hemisphere
        const int    basePhi = (s % 32) * 8;
        const double yMid    = basePhi + 4;
        const double xSec    = (s < 32) ? 72 : 24;
        tx.DrawLatex(xSec, yMid, Form("S%d", s));
      }

      c.SaveAs(outPng.string().c_str());
      return rot;
    }; // makePanel

    // ── centralities: keep one rotated copy per slice for the overview
    if (isMap && o->InheritsFrom(TH2::Class())) {
      auto p = makePanel(static_cast<TH2*>(o));
      _centralMaps.emplace(slice, std::move(p));   // silently overwrites duplicates
    }
    // ── any other 2‑D EMCal histogram → PNG with default helper
    else if (o->InheritsFrom(TH2::Class())) {
      fs::path out = cPath(root, slice, "EMCal") / (n + ".png");
      save2D(static_cast<TH2*>(o), out);
    }
    // ── plain 1‑D EMCal histogram
    else {
      fs::path out = cPath(root, slice, "EMCal") / (n + ".png");
      save1D(static_cast<TH1*>(o), out);
    }
    return true;
  }

private:
  // ── global Z‑axis limits (updated by makePanel, applied in destructor)
  double _zMin { std::numeric_limits<double>::max()  };  // smallest >0 bin
  double _zMax { std::numeric_limits<double>::lowest() }; // 99.5 % quantile
  static constexpr double kFracSaturation = 0.995;        // keep 99.5 %

  // one rotated η‑φ map per centrality slice
  std::unordered_map<std::string, std::unique_ptr<TH2F>> _centralMaps;
};


// ─── HCal QA – proportional η–φ hit‑maps (IHCal / OHCal) ─────────────
class HcalQA : public QA
{
 public:
  using QA::QA;

  // ==================================================================
  //  helper: safe‑clone (detaches from any directory immediately)
  // ==================================================================
  static TH2* cloneDetach(const TH2* src, const char* newName = nullptr)
  {
      auto* c = static_cast<TH2*>(src->Clone(newName ? newName : src->GetName()));
      if (!c) throw std::runtime_error("TH2::Clone() returned nullptr!");
      c->SetDirectory(nullptr);          // *** critical: break ownership ***
      return c;
  }

  // ==================================================================
  //  per‑object processing entry point
  // ==================================================================
  bool process(TObject* o) override
  {
    /* -------------------------------------------------------------- *
     *  0.  Accept only IHCal / OHCal histograms                      *
     * -------------------------------------------------------------- */
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string n   = o->GetName();
    const bool isI        = n.rfind("h_IHCAL_", 0) == 0;
    const bool isO        = n.rfind("h_OHCAL_", 0) == 0;
    if (!isI && !isO) return false;

    const std::string slice = sliceKey(n);                 // Inclusive / Cent_x_y
    const bool isMap        = n.find("_EtaPhiMap_") != std::string::npos;

    log::trace("HcalQA  → processing \"" + n +
               "\"  slice=" + slice +
               (isMap ? "  (map)" : "  (scalar)"));

    /* -------------------------------------------------------------- *
     *  1.  Lambda that renders one η–φ map panel                     *
     * -------------------------------------------------------------- */
    auto makePanel = [&](TH2* src, const fs::path& outPng)
    {
      try {
          constexpr int nPhi = 64, nEta = 24, px = 18;

          // 1.1 working clone & bad‑plate masking --------------------
          std::unique_ptr<TH2> h( cloneDetach(src) );
          h->SetStats(0); h->SetContour(99);

          for (int ip = 0; ip < nPhi; ++ip)
            for (int ie = 0; ie < nEta; ++ie)
              if (isBadHcalPlate(hcal_sector_from_idx(ie, ip),
                                 hcal_plate_from_idx (ie, ip)))
                h->SetBinContent(h->FindBin(ip, ie), -9999.);

          // 1.2 rotate  (η → X, φ → Y) ------------------------------
          std::unique_ptr<TH2F> rot(
              static_cast<TH2F*>(cloneDetach(h.get(),
                                             ("hRot_" + std::string(src->GetName())).c_str())));
          rot->SetBins(nEta, 0, nEta, nPhi, 0, nPhi);

          for (int ip = 1; ip <= nPhi; ++ip)
            for (int ie = 1; ie <= nEta; ++ie)
              rot->SetBinContent(ie, ip, h->GetBinContent(ip, ie));

          // cosmetics …
          rot->SetMinimum(1.);
          rot->GetXaxis()->SetTitle("#eta index");
          rot->GetYaxis()->SetTitle("#phi index");

          // 1.3 canvas & save ---------------------------------------
          const int cw = (kHCalCanvasW > 0) ? kHCalCanvasW : nEta * px;
          const int ch = (kHCalCanvasH > 0) ? kHCalCanvasH : nPhi * px;

          TCanvas c("c_hcal", "", cw, ch);
          c.SetRightMargin(0.16);  c.SetLeftMargin(0.08);
          c.SetBottomMargin(0.08); c.SetTopMargin(0.055);
          c.SetFixedAspectRatio();

          rot->Draw("COLZ");
          gPad->Update();                 // make sure user coords are frozen

          // vertical grid (η index)
          for (int x = 0; x <= nEta; ++x) {
              auto *lx = new TLine(x, 0, x, nPhi);      // heap → canvas keeps it
              lx->SetLineColor(kBlack);
              lx->SetLineWidth((x % 8 == 0) ? 3 : 1);   // heavy every 8 towers
              lx->Draw("same");
          }

          // horizontal grid (φ index)
          for (int y = 0; y <= nPhi; ++y) {
              auto *ly = new TLine(0, y, nEta, y);
              ly->SetLineColor(kBlack);
              ly->SetLineWidth((y % 2 == 0) ? 3 : 1);   // heavy every 2 towers
              ly->Draw("same");
          }

          c.Modified(); c.Update();        // register all new primitives
          c.SaveAs(outPng.string().c_str());
          log::trace("HcalQA  → wrote " + outPng.string());
      }
      catch (const std::exception& ex) {
          log::warn(std::string("HcalQA  WARN  failed to save panel for \"")
                    + src->GetName() + "\": " + ex.what());
      }
    };
    /* ---------- end makePanel ------------------------------------ */

    /* -------------------------------------------------------------- *
     *  2.  Per‑arm outputs                                           *
     * -------------------------------------------------------------- */
    fs::path subDir = fs::path("HCal") / (isI ? "IHCal" : "OHCal");
    fs::path out    = cPath(root, slice, subDir) / (n + ".png");

    if (isMap && o->InheritsFrom(TH2::Class()))
        makePanel(static_cast<TH2*>(o), out);
    else if (o->InheritsFrom(TH2::Class()))
        save2D(static_cast<TH2*>(o), out);
    else
        save1D(static_cast<TH1*>(o), out);

    /* -------------------------------------------------------------- *
     *  3.  Build totalHCal maps once both IHCal & OHCal are present  *
     * -------------------------------------------------------------- */
    if (isMap && o->InheritsFrom(TH2::Class()))
    {
      struct Pair {
          std::unique_ptr<TH2> i, o;
          ~Pair() {               // safety: make sure no directory owns us
              if (i) i->SetDirectory(nullptr);
              if (o) o->SetDirectory(nullptr);
          }
      };
      static std::unordered_map<std::string, Pair> cache;   // key = slice|name

      std::string baseName = std::regex_replace(
                                   n,
                                   std::regex("^h_[IO]HCAL_"),   // unify prefix
                                   "h_HCAL_");

      const std::string key = slice + "|" + baseName;
      Pair& p = cache[key];

      try {
          if (isI) p.i.reset( cloneDetach(static_cast<TH2*>(o)) );
          if (isO) p.o.reset( cloneDetach(static_cast<TH2*>(o)) );
      }
      catch (const std::exception& ex) {
          log::warn("HcalQA  WARN  clone failed for \"" + n + "\": " + ex.what());
          return true;   // skip – don’t kill the job
      }

      if (p.i && p.o)                         // have both arms – combine now
      {
        auto tot = std::unique_ptr<TH2>(
                       cloneDetach(p.i.get(),
                                   (std::string(p.i->GetName())+"_tot").c_str()));
        tot->Add(p.o.get());
        tot->SetTitle("totalHcal");
        fs::path outTot = cPath(root, slice, fs::path("HCal") / "totalHcal")
                              / (baseName + "_total.png");
        makePanel(tot.get(), outTot);

        cache.erase(key);                     // free – Pair dtor detaches dirs
      }
    }

    return true;
  }
};



// ───────────────── Event‑plane observables (sEPD) ──────────────────────
class SepdPlaneQA : public QA
{
 public:
  SepdPlaneQA(std::string t, fs::path b, const CentList& s)
  : QA(std::move(t), std::move(b), s) {}

  bool process(TObject* o) override
  {
    if (!o) {
      log::err("[SepdPlaneQA] nullptr TObject received – skipping.");
      return false;
    }
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string n = o->GetName();
    log::trace("[SepdPlaneQA] Inspecting \"" + n + "\".");

    /* keep only the Ψₙ & resolution‑proxy histograms */
    static const std::vector<std::string> keys = {
      "h_Psi1_sEPD", "h_Psi2_sEPD", "h_Psi3_sEPD",
      "h_Psi1_res_vs_Qsum", "h_Psi2_res_vs_Qsum", "h_Psi3_res_vs_Qsum",
      "p_R2_vs_cent"
    };
    const bool match = std::any_of(keys.begin(), keys.end(),
                                   [&](const std::string& k){ return n.rfind(k,0)==0; });
    if (!match) return false;

    /* build output path & ensure directory exists */
    fs::path out = cPath(root, sliceKey(n), "sEPD/EventPlaneQA") / (n + ".png");
    try {
      ensure_dir(out.parent_path());
    } catch (const std::exception& e) {
      log::err("[SepdPlaneQA] Failed to create directory \"" +
               out.parent_path().string() + "\": " + e.what());
      return false;
    }

    /* save histogram – wrap in try/catch to surface ROOT errors */
    try {
      save1D(static_cast<TH1*>(o), out);   // works for TH1, TH2, TProfile
      log::ok("[SepdPlaneQA] Saved → " + out.string());
    } catch (const std::exception& e) {
      log::err("[SepdPlaneQA] Exception while saving \"" + n + "\": " + e.what());
      return false;
    }
    return true;
  }
};



// ─────────  South / North combiner for 2‑D hit‑maps  ────────────────────
template<class DERIVED>
class NSDetectorQA : public QA
{
 public:
  NSDetectorQA(std::string           t,
               std::filesystem::path b,
               const CentList&       s,
               NSCache<MapPair>&     c)
  : QA(std::move(t), std::move(b), s), _cache(c) {}

  bool process(TObject* o) override
  {
    if (!o) {
      log::err("[NSDetectorQA] nullptr TObject received – skipping.");
      return false;
    }
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string hName = o->GetName();
    if (!DERIVED::accept(hName)) return false;

    const std::string slice = sliceKey(hName);
    const bool south = hName.find("_South_") != std::string::npos;
    log::trace("[NSDetectorQA] " + std::string(south ? "South" : "North") +
               " arm histogram \"" + hName + "\" accepted.");

    /* ------------------------------------------------------------------ *
     *  scalar histograms (charge spectra etc.) → save immediately        *
     * ------------------------------------------------------------------ */
    if (!o->InheritsFrom(TH2::Class()))
    {
      fs::path out = cPath(root, slice, DERIVED::subdir) / (hName + ".png");
      try {
        ensure_dir(out.parent_path());
        save1D(static_cast<TH1*>(o), out);
        log::ok("[NSDetectorQA] Saved 1‑D histo → " + out.string());
      } catch (const std::exception& e) {
        log::err("[NSDetectorQA] Failed to save \"" + hName + "\": " + e.what());
        return false;
      }
      return true;
    }

    /* ------------------------------------------------------------------ *
     *  2‑D hit‑maps: cache until both arms are present                    *
     * ------------------------------------------------------------------ */
    MapPair& mp = _cache[trig + slice];

    // clone so the original can be deleted without affecting us
    TH2* hClone = static_cast<TH2*>(o->Clone());
    if (!hClone) {
      log::err("[NSDetectorQA] Clone failed for \"" + hName + "\" – skipping.");
      return false;
    }
    hClone->SetDirectory(nullptr);
    hClone->SetStats(0);

    south ? mp.s = hClone : mp.n = hClone;
    if (!_cache.ready(trig + slice)) {
      log::trace("[NSDetectorQA] Waiting for partner arm to arrive ("
                 + slice + ").");
      return true;                         // partner not yet seen
    }

    MapPair in = _cache.pop(trig + slice); // take ownership
    auto tidy  = [](TH2* h)
    {
      h->SetMinimum(0.);
      h->SetLineColor(kBlack);
      h->SetLineWidth(1);
    };
    tidy(in.s); tidy(in.n);

    /* ------------------------------------------------------------------ *
     *  finished S–N canvas                                               *
     * ------------------------------------------------------------------ */
    fs::path png = cPath(root, slice, DERIVED::subdir)
                 / (DERIVED::fileName(trig) + ".png");

    try {
      ensure_dir(png.parent_path());
    } catch (const std::exception& e) {
      log::err("[NSDetectorQA] Cannot create output dir \"" +
               png.parent_path().string() + "\": " + e.what());
      return false;
    }

    TCanvas c("c_hit", "", 1200, 600);
    c.Divide(2, 1, 0.01, 0.01);

    auto drawPad = [&](TH2* h, const char* ttl)
    {
      gPad->SetRightMargin(0.20);
      gPad->SetLeftMargin (0.10);
      gPad->SetBottomMargin(0.10);
      gPad->SetTopMargin  (0.08);

      h->SetTitle(ttl);
      h->GetZaxis()->SetTitle("Counts");
      h->GetZaxis()->SetTitleOffset(1.3);
      h->Draw("POLZ");                     // works for both hex & polar
    };

    try {
      c.cd(1); drawPad(in.s, DERIVED::titleSouth);
      c.cd(2); drawPad(in.n, DERIVED::titleNorth);
      c.SaveAs(png.string().c_str());
      log::ok("[NSDetectorQA] Combined S/N map saved → " + png.string());
    } catch (const std::exception& e) {
      log::err("[NSDetectorQA] Error while drawing/saving \"" +
               png.string() + "\": " + e.what());
      return false;
    }
    return true;
  }

 private:
  NSCache<MapPair>& _cache;
};

/* ───────────── MBD tag – unchanged (aside from minor comments) ───────── */
struct MBDTag
{
  static bool accept(const std::string& s)
  {
    return s.rfind("h_MBD_Hitmap_",0)==0     ||
           s.rfind("h_charge_MBD", 0)==0     ||
           s.rfind("h_Qsum_MBD"  , 0)==0;
  }
  static constexpr const char* subdir = "MBD/otherQA";
  static std::string fileName(const std::string& t)
  { return "MBD_Hitmap_NS_" + t; }
  static constexpr const char* titleSouth = "MBD South";
  static constexpr const char* titleNorth = "MBD North";
};

/* ───────────── sEPD tag – hit‑maps + ΣQ spectra ───────────── */
struct sEPDTag
{
  static bool accept(const std::string& s)
  {
    return s.rfind("h_sEPD_Hitmap_",0)==0   ||   // φ–r hit‑maps
           s.rfind("h_Qsum_sEPD"  ,0)==0    ||   // ΣQ spectra
           s.rfind("h_towerQ_SEPD",0)==0;        // ΣQ (tower) spectra
  }
  static constexpr const char* subdir = "sEPD/OtherQA";
  static std::string fileName(const std::string& t)
  { return "sEPD_Hitmap_NS_" + t; }
  static constexpr const char* titleSouth = "South";
  static constexpr const char* titleNorth = "North";
};

/*  concrete type aliases – unchanged public names  */
using MbdQA  = NSDetectorQA<MBDTag>;
using SepdQA = NSDetectorQA<sEPDTag>;




// ╔══════════════════════════════════════════════╗
// ║               E v e n t  Q A                 ║
// ╚══════════════════════════════════════════════╝
class EventQA : public QA
{
public:
    using QA::QA;

    // ────────────────────────────────────────────────────────────────
    // 1. Per‑histogram processing
    // ────────────────────────────────────────────────────────────────
    bool process(TObject* o) override
    {
        if (!o->InheritsFrom(TH1::Class())) return false;

        const std::string n = o->GetName();
        const bool isVz   = (n.rfind("h_vertexZ_"  ,0) == 0);
        const bool isCent = (n.rfind("h_centrality_",0) == 0);
        if (!isVz && !isCent) return false;

        /* “…/output/<RUN>/<trigger>/EventQA/…”  (always Inclusive) */
        const fs::path outPng = isVz
            ? root / "MBD" / "zVertex"   / "VertexZ.png"
            : root / "centrality" / "Centrality.png";
        ensure_dir(outPng.parent_path());

        std::unique_ptr<TH1> h(static_cast<TH1*>(o->Clone()));
        h->SetDirectory(nullptr);
        h->SetStats(0);

        const std::string runID = root.parent_path().filename().string();

        // ============================================================
        // (A)  primary‑vertex z  –  robust iterative Gaussian fit
        // ============================================================
        if (isVz)
        {
            /* --- STEP‑0 : robust seed from quantiles ---------------- */
            double probs[3] = {0.16, 0.50, 0.84};
            double q[3];
            h->GetQuantiles(3, q, probs);
            double mu    = q[1];
            double sigma = 0.5*(q[2]-q[0]);               // 68 % width
            if (sigma <= 0) sigma = h->GetRMS();
            if (sigma <= 0) sigma = 1;

            TF1 g("g","gaus", mu-3*sigma, mu+3*sigma);
            g.SetLineColor(kRed+1); g.SetLineWidth(2);
            g.SetParameters(h->GetMaximum(), mu, sigma);

            /* --- STEP‑1 : iterative 2.5 σ shrink until convergence -- */
            constexpr int    kMaxIter   = 5;
            constexpr double kNSigmaFit = 2.5;
            constexpr double kTol       = 1e-3;
            bool   fitOK = false;

            for (int it = 0; it < kMaxIter; ++it)
            {
                const double lo = mu - kNSigmaFit*sigma;
                const double hi = mu + kNSigmaFit*sigma;
                g.SetRange(lo, hi);
                g.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sigma);

                TFitResultPtr res = h->Fit(&g,"Q0RSLL");
                fitOK = (int)res == 0;
                if (!fitOK) break;

                const double muNew    = g.GetParameter(1);
                const double sigmaNew = std::fabs(g.GetParameter(2));

                const bool conv =
                       std::fabs(muNew   - mu)    < kTol*sigma &&
                       std::fabs(sigmaNew- sigma) < kTol*sigma;
                mu = muNew;  sigma = sigmaNew;
                if (conv) break;
            }

            const double muErr  = fitOK ? g.GetParError(1) : 0;
            const double sigErr = fitOK ? g.GetParError(2) : 0;

            /* --- per‑run PNG --------------------------------------- */
            {
                TCanvas c("c_vz","", 900, 600);
                h->Draw();
                if (fitOK) g.Draw("SAME");

                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.04);
                tx.DrawLatex(0.15,0.86,Form("#mu = %.2f #pm %.2f cm", mu,  muErr));
                tx.DrawLatex(0.15,0.80,Form("#sigma = %.2f #pm %.2f cm", sigma, sigErr));
                c.SaveAs(outPng.string().c_str());
            }

            /* --- cache for overlays & run‑summary ------------------ */
            if (runID != "Combined")
            {
                VzPoint& p = s_points[runID];
                p.mu = mu; p.muErr = muErr; p.sigma = sigma; p.sigmaErr = sigErr;
                p.hist.reset(static_cast<TH1*>(h->Clone()));
                p.hist->SetDirectory(nullptr);
            }

            /* ---------- NEW: keep event count for final table ------ */
            const long long nEvt = static_cast<long long>(h->GetEntries());
            auto& slot = s_evtCounts[runID];
            if (nEvt > slot) slot = nEvt;        // keep the largest if multiple triggers
        }

        // ============================================================
        // (B)  centrality spectrum – plain plot, normalised
        // ============================================================
        else
        {
            TCanvas c("c_cent","",900,600);
            h->Draw();

            /* print the run‑number in the upper‑right corner */
            TLatex tx;
            tx.SetNDC();           // Normalised device coordinates
            tx.SetTextAlign(31);   // right‑aligned, top‑aligned
            tx.SetTextSize(0.04);
            tx.DrawLatex(0.97,0.94, runID.c_str());

            c.SaveAs(outPng.string().c_str());

            if (runID != "Combined")
            {
                std::unique_ptr<TH1> cp(static_cast<TH1*>(h->Clone()));
                cp->SetDirectory(nullptr);
                if (cp->GetEntries() > 0) cp->Scale(1.0 / cp->GetEntries());
                s_centHists[runID] = std::move(cp);
            }
        }
        return true;
    }

    // ────────────────────────────────────────────────────────────────
    // 2. Final summary (executed once, after the “Combined” pass)
    // ────────────────────────────────────────────────────────────────
    ~EventQA() override
    {
        const std::string runID = root.parent_path().filename().string();
        if (runID != "Combined" || s_summaryWritten) return;
        s_summaryWritten = true;

        const fs::path outDir = root.parent_path();      // “…/Combined”

        /* ---------- helper: reproducible colour stream ------------ */
        auto nextColour = [](){
            static int idx = 0;
            static int palette[] = {kBlue+1,kRed+1,kGreen+2,kMagenta+2,
                                    kCyan+2,kOrange+1,kViolet,kAzure+2,
                                    kPink+1,kTeal+2};
            return palette[(idx++) % (sizeof(palette)/sizeof(int))];
        };

        /* (A) vertex‑Z overlay ------------------------------------- */
        if (!s_points.empty())
        {
            TCanvas c("c_vz_overlay","Primary‑vertex Z – all runs",900,600);
            TLegend leg(0.68,0.57,0.88,0.88); leg.SetBorderSize(0);

            bool first = true;
            for (auto& [run,p] : s_points)
            {
                Color_t col = nextColour();
                p.hist->SetLineColor(col); p.hist->SetLineWidth(2);
                p.hist->Draw(first ? "HIST" : "HIST SAME");
                leg.AddEntry(p.hist.get(), run.c_str(), "l");
                first = false;
            }
            leg.Draw();
            fs::path vzOverlay = root / "MBD" / "zVertex" / "VertexZ_AllRuns.png";
            ensure_dir(vzOverlay.parent_path());
            c.SaveAs(vzOverlay.string().c_str());
        }

        /* (B) centrality overlay ----------------------------------- */
        if (s_centHists.size() > 1)
        {
            TCanvas c("c_cent_overlay","Centrality – all runs",900,600);
            TLegend leg(0.5,0.45,0.75,0.65); leg.SetBorderSize(0);

            bool first = true;
            for (auto& [run,h] : s_centHists)
            {
                Color_t col = nextColour();
                h->SetLineColor(col); h->SetLineWidth(2);
                h->Draw(first ? "HIST" : "HIST SAME");
                leg.AddEntry(h.get(), run.c_str(), "l");
                first = false;
            }
            leg.Draw();
            fs::path centOverlay = root / "centrality" / "Centrality_AllRuns.png";
            ensure_dir(centOverlay.parent_path());
            c.SaveAs(centOverlay.string().c_str());
        }

        /* (C)  μ,σ  versus run number ------------------------------ */
        if (s_points.size() > 1)
        {
            std::vector<int> runs;
            for (auto& [r,_] : s_points)
                if (std::all_of(r.begin(),r.end(),::isdigit))
                    runs.push_back(std::stoi(r));
            std::sort(runs.begin(), runs.end());

            const int n = runs.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            for (int i = 0; i < n; ++i)
            {
                const auto& p = s_points[std::to_string(runs[i])];
                x[i]=runs[i]; yMu[i]=p.mu; eMu[i]=p.muErr;
                ySi[i]=p.sigma; eSi[i]=p.sigmaErr;
            }

            auto gMu = std::make_unique<TGraphErrors>(n,x.data(),yMu.data(),nullptr,eMu.data());
            auto gSi = std::make_unique<TGraphErrors>(n,x.data(),ySi.data(),nullptr,eSi.data());
            gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
            gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

            /* dynamic Y‑ranges */
            double minMu = yMu[0] - eMu[0], maxMu = yMu[0] + eMu[0];
            double maxSig = ySi[0] + eSi[0];
            for (int i = 1; i < n; ++i)
            {
                minMu  = std::min(minMu , yMu[i] - eMu[i]);
                maxMu  = std::max(maxMu , yMu[i] + eMu[i]);
                maxSig = std::max(maxSig, ySi[i] + eSi[i]);
            }
            double absMu = std::max(std::fabs(minMu), std::fabs(maxMu));
            if (absMu <= 0.) absMu = 0.01;
            gMu->SetMinimum(-1.10 * absMu);
            gMu->SetMaximum( 1.10 * absMu);
            gSi->SetMinimum(0.0);
            gSi->SetMaximum(1.10 * maxSig);

            TCanvas c("c_mu_sigma_vs_run","vertex‑Z  #mu,#sigma  vs run",900,800);

            TPad* p1 = new TPad("p1","",0,0.35,1,1);
            p1->SetBottomMargin(0.02); p1->Draw(); p1->cd();
            gMu->SetTitle(";Run number;#mu  [cm]");
            gMu->Draw("AP");

            c.cd();
            TPad* p2 = new TPad("p2","",0,0,1,0.32);
            p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30);
            p2->Draw(); p2->cd();
            gSi->SetTitle(";Run number;#sigma  [cm]");
            gSi->Draw("AP");

            fs::path pngRun = root / "MBD" / "zVertex" / "VertexZ_MeanSigma_vs_Run.png";
            ensure_dir(pngRun.parent_path());
            c.SaveAs(pngRun.string().c_str());
        }
    }

    // ────────────────────────────────────────────────────────────────
    // 3. Static containers (shared across ALL EventQA instances)
    // ────────────────────────────────────────────────────────────────
    struct VzPoint {
        double mu{}, muErr{}, sigma{}, sigmaErr{};
        std::unique_ptr<TH1> hist;
    };

    /* runID → fit results & histogram */
    static inline std::unordered_map<std::string,VzPoint>         s_points;
    /* runID → normalised centrality histogram */
    static inline std::unordered_map<std::string,std::unique_ptr<TH1>> s_centHists;

    /* ---------- runID → event count (entries in h_vertexZ_) -------- */
    static inline std::unordered_map<std::string,long long>       s_evtCounts;

    static inline bool  s_summaryWritten = false;

    /* accessor for final summary in main() */
    static const auto& eventCounts() { return s_evtCounts; }
};


// ────────────────────────────────────────────────────────────────────
//  Jet‑QA module
//      • “generalHistos” : every raw histogram & projection
//      • “summary”       : jet‑yield vs E_T   and   jet‑yield vs centrality
// ────────────────────────────────────────────────────────────────────
class JetQA : public QA
{
 public:
  using QA::QA;

  // =================================================================
  // 1. Per‑histogram processing (called many times per file)
  // =================================================================
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string n = o->GetName();
    const bool is1D =  n.rfind("h_maxJetEt_"         ,0) == 0;
    const bool is2D =  n.rfind("h_leadEt_vs_subEt_"  ,0) == 0;
    const bool is3D =  n.rfind("h_jetEt_area_nConst_",0) == 0;
    if (!is1D && !is2D && !is3D) return false;

    const std::string slice = sliceKey(n);                  // Inclusive / Cent_…
    const std::string rLab  = radiusTag(n);                 // r02 / r04 …

    /* -------- 1.1   *all* raw plots into …/generalHistos/ ---------- */
    fs::path baseGen = cPath(root, slice,
                             fs::path("jetQA/generalHistos") / rLab);
    ensure_dir(baseGen);

    /* the original plotting helpers are reused unchanged ------------- */
    if (is1D)          return save1Dplot(static_cast<TH1*>(o), baseGen, n);
    if (is2D)          return handle2D (static_cast<TH2*>(o), baseGen, n);
                       return handle3D (static_cast<TH3*>(o), baseGen, n);
  }

  // =================================================================
  // 2. Final summary – executed once at the very end of the run
  // =================================================================
  ~JetQA() override
  {
    /* ------------------------------------------------------------ *
     *  Collect every “h_maxJetEt_rXX_…Inclusive…” histogram that   *
     *  has statistics, integrate above a threshold and build:      *
     *     – yield(E_T)  per slice                                  *
     *     – yield vs centrality (cent = bin centre)                *
     * ------------------------------------------------------------ */
    constexpr double kEtMin = 20.0;                // GeV threshold
    std::map<std::string,std::unique_ptr<TH1>> maxEtPerSlice;   // slice→hist

    for (const TObject* obj : *gROOT->GetList())
    {
      const TH1* h = dynamic_cast<const TH1*>(obj);
      if (!h) continue;

      const std::string n = h->GetName();
      if (n.rfind("h_maxJetEt_",0)!=0) continue;         // keep only the 1‑D family
      if (h->Integral()<=0)      continue;

      const std::string slice = sliceKey(n);             // Inclusive / Cent_x_y
      /* clone under our ownership so original may disappear ---------- */
      maxEtPerSlice[slice].reset( static_cast<TH1*>(h->Clone()) );
      maxEtPerSlice[slice]->SetDirectory(nullptr);
    }

    if (maxEtPerSlice.empty()) return;          // nothing to summarise

    /* 2.1  yield(E_T) – one graph per slice ------------------------ */
    fs::path dirSum = root / "jetQA/summary";
    ensure_dir(dirSum);

    TCanvas cYield("c_yieldEt","Jet yield vs E_{T}",1100,850);
    TLegend leg(0.15,0.70,0.45,0.88); leg.SetBorderSize(0);

    int colList[]{kRed+1,kBlue+2,kGreen+2,kMagenta+2,kOrange+1};
    int iCol=0;

    std::vector<double> xCent, yYield;          // for yield‑vs‑cent graph

    for (auto& [slice, h] : maxEtPerSlice)
    {
        /* determine centrality bin centre ---------------------------- */
        double xC = 50.0;                                   // default = Inclusive
        if (slice.rfind("Cent_", 0) == 0)
        {
          std::smatch m;
          std::regex  re(R"(Cent_([0-9]+)_([0-9]+))");
          if (std::regex_match(slice, m, re))
            xC = 0.5 * (std::stod(m[1]) + std::stod(m[2]));
        }

        /* integrate Y(E_T>EtMin) ------------------------------------ */
        int    binMin = h->FindBin(kEtMin);
        double yield  = h->Integral(binMin, h->GetNbinsX());

        /* real event count for this slice --------------------------- */
        long long nEv = static_cast<long long>(h->GetEntries());
        if (nEv == 0)                     // skip empty slices (avoids div-by-zero)
          continue;

        yield /= nEv;                     // per-event yield

        xCent.push_back(xC);
        yYield.push_back(yield);

        /* build differential yield dN/dE_T --------------------------- */
        auto hDiff = std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone()));
        hDiff->Scale(1.0 / nEv, "width");       // safe: nEv > 0
        hDiff->SetLineColor(colList[iCol % 5]);
        hDiff->SetLineWidth(2);

        hDiff->SetTitle(Form("dN/dE_{T} – %s", slice.c_str()));
        hDiff->GetYaxis()->SetTitle("1/N_{ev}  dN/dE_{T}  [GeV^{-1}]");

        hDiff->Draw(iCol == 0 ? "HIST" : "HIST SAME");
        leg.AddEntry(hDiff.get(), slice.c_str(), "l");
        _owned1D.push_back(std::move(hDiff));
        ++iCol;
    }
    leg.Draw();
    cYield.SetLogy();
    cYield.SaveAs( (dirSum/"JetYield_vs_Et_AllSlices.png").string().c_str() );

    /* 2.2  yield vs centrality ------------------------------------- */
    if (xCent.size()>1)
    {
      auto gCent = std::make_unique<TGraphErrors>(xCent.size());
      for (std::size_t i=0;i<xCent.size();++i)
        gCent->SetPoint(i, xCent[i], yYield[i]);

      TCanvas cCent("c_yieldCent","Jet yield vs centrality",900,700);
      gCent->SetTitle(Form("Jet yield  E_{T}>%.0f GeV",kEtMin));
      gCent->GetXaxis()->SetTitle("centrality [%]");
      gCent->GetYaxis()->SetTitle("jets / event");
      gCent->SetMarkerStyle(kFullCircle); gCent->SetLineWidth(2);
      gCent->Draw("AP");

      cCent.SaveAs( (dirSum/"JetYield_vs_Centrality.png").string().c_str() );
      _ownedGraphs.push_back(std::move(gCent));
    }
  }

 private:
  // ======================= helper functions =========================
  static std::string radiusTag(const std::string& hname)
  {
    std::smatch m; std::regex re(R"(_(r[0-9]+|R[0-9]+)_)");
    return std::regex_search(hname,m,re) ? m[1].str() : "UnknownR";
  }

  bool save1Dplot(TH1* h, const fs::path& dir, const std::string& hname)
  {
    TCanvas c; c.SetLogy(); h->SetStats(0);
    h->SetTitle(makeTitle(hname).c_str());
    h->Draw();
    c.SaveAs( (dir/(hname+".png")).string().c_str() );
    return true;
  }
    
  // =============== 2‑D ==================================================
  bool handle2D(TH2* h, const fs::path& dir, const std::string& hname)
  {
    h->SetTitle(makeTitle(hname).c_str());
    TCanvas c; h->SetStats(0);
    h->Draw("COLZ");

    // y = x guideline
    const double xmax = h->GetXaxis()->GetXmax();
    TLine diag(0,0, xmax, xmax);
    diag.SetLineStyle(2); diag.SetLineWidth(2); diag.Draw();

    ensure_dir(dir);
    c.SaveAs((dir / (hname + ".png")).string().c_str());
    return true;
  }

  // =============== 3‑D ==================================================
  bool handle3D(TH3* h3, const fs::path& dir, const std::string& hname)
  {
    h3->SetTitle(makeTitle(hname).c_str());

    // main 3‑D view
    save3D(h3, dir / (hname + "_3D.png"));

    // orthogonal projections
    saveProjection(h3,"yx", dir / (hname + "_Et_vs_Area.png"));    // E_T vs A
    saveProjection(h3,"xz", dir / (hname + "_Et_vs_Nconst.png"));  // E_T vs N
    saveProjection(h3,"yz", dir / (hname + "_Area_vs_Nconst.png"));// A  vs N
    return true;
  }

  // ===== local helpers (only visible inside JetQA) ==================
  static void save3D(TH3* h, const fs::path& png)
  {
    log::trace("save3D → " + png.string());
    TCanvas c("c3D","",1200,1000);
    c.SetRightMargin(0.18);
    h->SetStats(0); h->SetContour(99);
    h->Draw("BOX2Z");                      // nice semi‑transparent boxes
    ensure_dir(png.parent_path());
    c.SaveAs(png.string().c_str());
  }

  static void saveProjection(TH3* h3,
                               const char* axes,
                               const fs::path& png)
    {
        TH1*   tmp = h3->Project3D(axes);                 // ROOT gives TH1*
        auto h2 = std::unique_ptr<TH2>(                  // take ownership
                      dynamic_cast<TH2*>(tmp));          // safe cast → TH2*

        if (!h2) {                                       // should never happen
            log::warn(std::string("Projection '")+axes+
                       "' of '"+h3->GetName()+"' is not TH2 – skipped");
            delete tmp;                                  // avoid leak
            return;
        }

        h2->SetDirectory(nullptr);                       // detach from gDirectory
        h2->SetStats(0);
        h2->SetTitle((std::string(h3->GetTitle())+
                     "  –  "+axes).c_str());

        save2D(h2.get(), png);
  }

  static std::string makeTitle(const std::string& hname)
  {
    // “h_jetEt_area_nConst_r02_MBD_NandS_geq_2”  →  “r02  (MBD_NandS_geq_2)”
    std::smatch m; std::regex re(R"(_(r[0-9]+|R[0-9]+).+?_(MBD.+))");
    return std::regex_search(hname,m,re) ? (m[1].str()+"  ("+m[2].str()+')')
                                         : hname;
  }
    
  // containers keeping produced objects alive
  std::vector<std::unique_ptr<TH1>>         _owned1D;
  std::vector<std::unique_ptr<TGraphErrors>> _ownedGraphs;
};


// ╔═══════════════════════════════════════════════════════════════════╗
// ║                V n   P l o t   P o s t ‑ P r o c e s s o r        ║
// ╚═══════════════════════════════════════════════════════════════════╝
class VnPlotQA : public QA
{
 public:
  VnPlotQA(std::string trig,
           std::filesystem::path base,
           const CentList& slices) :
      QA(std::move(trig), base, slices),
      _outDir(base / "vNana")
  {
    std::filesystem::create_directories(_outDir);
  }

  // ─────────────────────────────── 1. cache every TProfile ──────────
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TProfile::Class())) return false;

    /* accepted names
     *   – calorimeters        p_v<n>_<DET>_<lo>_<hi>_<trig>     (DET = CEMC_S …)
     *   – jets (new)          p_v<n>_JET_rXX_<lo>_<hi>_<trig>
     */
    static const std::regex re(
      R"(p_v([123])_((?:[A-Za-z0-9]+_[NS])|(?:JET_r[0-9]{2}))_([0-9]+)_([0-9]+)_(.+))");

    std::smatch m;
    const std::string h = o->GetName();
    if (!std::regex_match(h, m, re)) return false;

    const int         nHarm   = std::stoi(m[1]);        // 1 / 2 / 3
    const std::string detTag  = m[2];                   // CEMC_S …  OR  JET_r04
    const std::string centKey = m[3].str() + '_' + m[4].str();   // e.g. 10_20
    const std::string trigLab = m[5];

    _raw[trigLab][detTag][nHarm][centKey]
        .push_back(static_cast<TProfile*>(o));
    return true;
  }

  // ─────────────────────────────── 2. heavy lifting on exit ─────────
  ~VnPlotQA() override { writeCanvases(); }

 private:
  // =============================== helpers =============================
  static std::string detBase(const std::string& tag)
  {
      /* calorimeters -------------------------------------------------- */
      if (tag.find("CEMC")  != std::string::npos) return "EMCal";
      if (tag.find("IHCAL") != std::string::npos) return "IHCal";
      if (tag.find("OHCAL") != std::string::npos) return "OHCal";
      if (tag == "HCAL_S" || tag == "HCAL_N")      return "totalHCal";
      if (tag == "ALL_S"  || tag == "ALL_N")       return "totalCalo";

      /* jets ---------------------------------------------------------- */
      if (tag.rfind("JET_r",0)==0)                 return "jetvN";
      return tag;                                  // fallback
  }

  static std::string regionOf(const std::string& tag)
  {
      if (tag.rfind("JET_r",0)==0) return tag.substr(4,3); // r02 / r04 …
      return (tag.back()=='S') ? "South" : "North";
  }

  /* turn a TProfile → TGraphErrors, apply resolution if available ---- */
  std::unique_ptr<TGraphErrors>
  graphFromProf(const TProfile* p,
                int            nHarm,
                const std::string& centKey,
                const std::string& trigLab) const
  {
      const int nb = p->GetNbinsX();
      auto g = std::make_unique<TGraphErrors>(nb);

      /* obtain sub‑event resolution Rn(cent) once -------------------- */
      double Rn = 1.0;                       // default = no correction
      const std::string rName = Form("p_R%d_vs_cent_%s", nHarm, trigLab.c_str());
      if (auto* rProf = static_cast<TProfile*>(gROOT->FindObject(rName.c_str())))
      {
         const auto pos = centKey.find('_');
         const double cMid = 0.5*( std::stod(centKey.substr(0,pos)) +
                                   std::stod(centKey.substr(pos+1)) );
         const int bin = rProf->FindBin(cMid);
         const double val = rProf->GetBinContent(bin);
         if (val > 1e-6) Rn = val;
      }

      /* copy points --------------------------------------------------- */
      for (int i=1;i<=nb;++i)
      {
          const double xLo=p->GetXaxis()->GetBinLowEdge(i);
          const double xHi=p->GetXaxis()->GetBinUpEdge (i);
          const double x  =0.5*(xLo+xHi);
          const double ex =0.5*(xHi-xLo);

          double y  = p->GetBinContent(i) / Rn;      // corrected
          double ey = p->GetBinError  (i) / Rn;

          g->SetPoint     (i-1, x, y);
          g->SetPointError(i-1, ex, ey);
      }
      g->SetLineWidth(2); g->SetMarkerStyle(kFullCircle);
      return g;
  }

  /* averaged (pT‑integrated) v̅n with resolution ---------------------- */
  static std::pair<double,double> meanProf(const TProfile* p,double Rn)
  {
      double num=0,den=0,err2=0;
      for (int i=1;i<=p->GetNbinsX();++i)
      {
          const double y=p->GetBinContent(i)/Rn;
          const double e=p->GetBinError(i)/Rn;
          const double w=p->GetBinEntries(i);
          if (w<=0) continue;
          num+=w*y; den+=w; err2+=(w*e)*(w*e);
      }
      const double m=(den>0)?num/den:0;
      const double er=(den>0)?std::sqrt(err2)/den:0;
      return {m,er};
  }

  // ============================ canvas writer ==========================
  void writeCanvases()
  {
    /* --- regroup cache: split into detector base / region / harmonic -- */
    using ProfVec = std::vector<TProfile*>;
    using CentMap = std::map<std::string,ProfVec>;
    using HarmMap = std::map<int,CentMap>;
    using RegMap  = std::unordered_map<std::string,HarmMap>;  // region
    using DetMap  = std::unordered_map<std::string,RegMap>;   // detector

    std::unordered_map<std::string,DetMap> grp;               // by trigger

    for (auto& [trig, detMap] : _raw)
      for (auto& [detTag, hMap] : detMap)
      {
        const std::string det  = detBase(detTag);
        const std::string reg  = regionOf(detTag);
        RegMap& rmap = grp[trig][det];

        for (auto& [n,cMap] : hMap)
          for (auto& [cent,v] : cMap)
            rmap[reg][n][cent].insert(rmap[reg][n][cent].end(),
                                      v.begin(),v.end());
      }

    /* colour palette -------------------------------------------------- */
    const int colTbl[]{kRed+1,kBlue+2,kGreen+2,kMagenta+2,kCyan+2,kOrange+1};
    const int nCol = sizeof(colTbl)/sizeof(int);

    /* loop over regrouped structure ---------------------------------- */
    for (const auto& [trig, detMap] : grp)
      for (const auto& [det, regMap] : detMap)
        for (const auto& [reg, harmMap] : regMap)
          for (const auto& [n, centMap] : harmMap)
          {
            /* ===== (1) v_n(pT) : one plot per centrality ============= */
            for (const auto& [cent, vec] : centMap)
            {
              if (vec.empty()) continue;
              TCanvas c("c","",1100,850); c.SetGrid();

              auto g = graphFromProf(vec.front(), n, cent, trig);
              g->SetTitle(Form("v_{%d}(p_{T}) – %s %s (Cent %s %%)",
                               n,det.c_str(),reg.c_str(),cent.c_str()));
              g->Draw("AP");

              saveCanvas(c,{det,reg,Form("v%d",n),"Cent_"+cent},
                         Form("v%d_%s_%s_cent%s.png",
                              n,det.c_str(),reg.c_str(),cent.c_str()));
            }

            /* ===== (2) all‑cent plot ================================= */
            {
              TCanvas c("c_all","",1100,850); c.SetGrid();
              TLegend leg(0.15,0.70,0.45,0.88); leg.SetBorderSize(0);
              int colIdx=0; double yMax=0;

              for (const auto& [cent, vec] : centMap)
              {
                if (vec.empty()) continue;
                auto g = graphFromProf(vec.front(), n, cent, trig);
                int col = colTbl[colIdx++%nCol];
                g->SetLineColor(col); g->SetMarkerColor(col);
                g->SetTitle(Form("v_{%d}(p_{T}) – %s %s (all cent)",
                                 n,det.c_str(),reg.c_str()));
                g->Draw(colIdx==1?"APL":"PL SAME");
                leg.AddEntry(g.get(),Form("Cent %s %%",cent.c_str()),"pl");
                yMax = std::max(yMax,
                                *std::max_element(g->GetY(), g->GetY()+g->GetN()));
                _ownedGraphs.push_back(std::move(g));
              }
              if (yMax>0)
              {
                c.Update();
                if(auto* fr=static_cast<TH1*>(c.GetPrimitive("htemp")))
                  fr->SetMaximum(1.15*yMax);
              }
              leg.Draw();
              saveCanvas(c,{det,reg,Form("v%d",n),"summaryPlots"},
                         Form("v%d_%s_%s_allCent.png",
                              n,det.c_str(),reg.c_str()));
            }

            /* ===== (3) pT‑integrated v̅_n vs centrality ============== */
            {
              auto gCent = std::make_unique<TGraphErrors>();
              int ip=0;
              for (const auto& [cent, vec] : centMap)
              {
                if (vec.empty()) continue;

                const auto pos=cent.find('_');
                const double cMid = 0.5*( std::stod(cent.substr(0,pos)) +
                                          std::stod(cent.substr(pos+1)) );

                const std::string rName =
                    Form("p_R%d_vs_cent_%s",n,trig.c_str());
                double Rn=1.0;
                if (auto* r = static_cast<TProfile*>(gROOT->FindObject(rName.c_str())))
                Rn = r->GetBinContent(r->FindBin(cMid));

                const auto [mu,er] = meanProf(vec.front(),Rn);

                gCent->SetPoint     (ip, cMid, mu);
                gCent->SetPointError(ip, 0.5*(std::stod(cent.substr(pos+1))-
                                              std::stod(cent.substr(0,pos))), er);
                ++ip;
              }

              if (gCent->GetN()>0)
              {
                TCanvas c("c_cent","",1000,800); c.SetGrid();
                gCent->SetTitle(Form("v_{%d} vs centrality – %s %s",
                                     n,det.c_str(),reg.c_str()));
                gCent->SetMarkerStyle(kFullCircle); gCent->SetLineWidth(2);
                gCent->Draw("AP");
                _ownedGraphs.push_back(std::move(gCent));

                saveCanvas(c,{det,reg,Form("v%d",n),"summaryPlots"},
                           Form("vbar%d_%s_%s_vsCent.png",
                                n,det.c_str(),reg.c_str()));
              }
            }
          }
  }

  // ---------- save helper ---------------------------------------------
  void saveCanvas(TCanvas& c,
                  std::initializer_list<std::string> path,
                  const std::string& file) const
  {
      std::filesystem::path dir=_outDir;
      for (auto& p : path) dir/=p;
      std::filesystem::create_directories(dir);
      c.SaveAs((dir/file).c_str());
  }

  // ---------- data members --------------------------------------------
  std::filesystem::path _outDir;

  using ProfVec = std::vector<TProfile*>;
  std::unordered_map<
      std::string,                       // trigger
      std::unordered_map<
          std::string,                   // detTag
          std::map<
              int,                       // harmonic n
              std::map<std::string,ProfVec> > > > _raw;

  std::vector<std::unique_ptr<TGraphErrors>> _ownedGraphs;
  std::vector<std::unique_ptr<TProfile>>     _ownedProfiles;
};


/* ────────────────────────────────────────────────────────────────────
 *  9.  MAIN DRIVER  –  per‑run analysis in parallel + combined pass
 *      (process pool implementation – safe for all ROOT classes)
 * ────────────────────────────────────────────────────────────────── */

/* ------------------------------------------------------------------ */
/*  One complete QA pass for a single ROOT file                        */
/* ------------------------------------------------------------------ */
void runOneQaPass(const std::string& inFile,
                  const std::string& outBase)
{
  /* retain the existing globals – each worker process owns its copy */
  kInputFile  = inFile;
  kOutputBase = outBase;

  gStyle->SetOptStat(0);

  log::banner("sPHENIX Run‑24 Au+Au QA – Enhanced Macro");

  std::unique_ptr<TFile> in(TFile::Open(kInputFile.c_str(), "READ"));
  if (!in || in->IsZombie()) { log::err("Cannot open " + kInputFile); return; }
  log::ok("Input file opened");
    
  // ------------------------------------------------------------------
  // 0‑bis.  Quick statistics sanity‑check
  //        → skip this run if *every* histogram is empty
  // ------------------------------------------------------------------
  bool hasStatistics = false;                 // assume “all empty” for now
    {
      TIter itTop(in->GetListOfKeys());
      while (auto* kDir = dynamic_cast<TKey*>(itTop())) {
        if (strcmp(kDir->GetClassName(), "TDirectoryFile")) continue;

        const std::string trg = kDir->GetName();        // e.g. "MBD_NandS_geq_2"
        if (!kTriggersWanted.count(trg)) continue;      // ignore unwanted triggers

        TDirectory* dTrig = static_cast<TDirectory*>(kDir->ReadObj());
        TIter itH(dTrig->GetListOfKeys());

        while (auto* kHist = dynamic_cast<TKey*>(itH())) {
          std::unique_ptr<TObject> obj(kHist->ReadObj());   // RAII – autodelete
          if (!obj->InheritsFrom(TH1::Class())) continue;   // TH2/TH3 inherit TH1

          TH1* h = static_cast<TH1*>(obj.get());
          if (h->GetEntries() > 0 && h->Integral() > 0) {   // non‑empty hist found
            hasStatistics = true;
            break;
          }
        }
        if (hasStatistics) break;               // early exit if anything has stats
      }
    }

    if (!hasStatistics) {
      log::warn("RUN IS SKIPPED! -> every histogram in \"" + inFile +
                "\" is empty!! (zero entries / zero integral).");
      return;                                    // ← abort QA early → no output
  }
  // ------------------------------------------------------------------

  CentList slices = discoverSlices(in.get());
  slices.erase(
        std::remove_if(slices.begin(), slices.end(),
                       [](const std::string& s){ return s=="0_100"; }),
        slices.end());
  {
    std::ostringstream o; o << "Centrality slices: ";
    for (auto& s : slices) o << s << "  ";
    log::info(o.str());
  }

  /* -----------------------------  PASS‑0 : catalogue  ------------- */
  log::banner("Pass 0 – Catalogue");
  fs::path catTxt = fs::path(kOutputBase) / "AllHistogramNames.txt";
  ensure_dir(catTxt.parent_path());
  std::ofstream cat(catTxt);
  cat << "# Histogram catalogue for " << kInputFile << "\n";

  std::unordered_map<string,int> hCnt;
  TIter it0(in->GetListOfKeys());
  while (auto* kd = dynamic_cast<TKey*>(it0())) {
    if (strcmp(kd->GetClassName(), "TDirectoryFile")) continue;
    string trg = kd->GetName(); if (!kTriggersWanted.count(trg)) continue;
    TDirectory* d = static_cast<TDirectory*>(kd->ReadObj());
    TIter itH(d->GetListOfKeys());
    while (auto* kh = dynamic_cast<TKey*>(itH())) {
      cat << "[" << trg << "] " << kh->GetName() << "\n"; ++hCnt[trg];
    }
  }
  log::ok("Histogram list written → " + catTxt.string());

  /* -----------------------------  PASS‑1 : QA  -------------------- */
  log::banner("Pass 1 – QA Production");

  fs::path csvPath = fs::path(kOutputBase) / "InvariantMassSummary.csv";
  ensure_dir(csvPath.parent_path());
  std::ofstream csv(csvPath);
  csv << "trigger,E,Chi,Asym,pTlo,pThi,meanPi0,errPi0,sigmaPi0,errSigmaPi0,"
         "meanEta,errEta,sigmaEta,errSigmaEta\n";

  NSCache<MapPair> mbdCache, sepdCache;
  struct Cnt { int tot = 0, used = 0; };
  std::unordered_map<string,Cnt> stat;

  TIter itDir(in->GetListOfKeys());
  while (auto* kd = dynamic_cast<TKey*>(itDir())) {
    if (strcmp(kd->GetClassName(), "TDirectoryFile")) continue;
    string trg = kd->GetName(); if (!kTriggersWanted.count(trg)) continue;
    TDirectory* dTrig = static_cast<TDirectory*>(kd->ReadObj());

      for (auto& s : slices) {
        fs::path b = fs::path(kOutputBase) / trg;
        for (auto sub : { "correlations", "vNana", "centrality",
                            "EMCal/otherQA", "EMCal/invMassQA", "EMCal/invMassQA/cutQA",
                            "HCal/IHCal", "HCal/OHCal", "HCal/totalHCal",
                            "MBD/otherQA", "MBD/zVertex",
                            "sEPD/OtherQA", "sEPD/EventPlaneQA",
                            "jetQA/generalHistos", "jetQA/summary"})
            ensure_dir(b / sub);
      }

    /* QA modules ---------------------------------------------------- */
    std::vector<std::unique_ptr<QA>> qa;
    fs::path base = fs::path(kOutputBase) / trg;
    qa.emplace_back(std::make_unique<Pi0QA >(trg,base,slices,csv));
    qa.emplace_back(std::make_unique<CorrQA>(trg,base,slices));
    qa.emplace_back(std::make_unique<EmcalQA>(trg,base,slices));
    qa.emplace_back(std::make_unique<HcalQA >(trg,base,slices));
    qa.emplace_back(std::make_unique<MbdQA  >(trg,base,slices,mbdCache));
    qa.emplace_back(std::make_unique<SepdQA >(trg,base,slices,sepdCache));
    qa.emplace_back(std::make_unique<SepdPlaneQA>(trg,base,slices));
    qa.emplace_back(std::make_unique<JetQA >(trg,base,slices));
    qa.emplace_back(std::make_unique<VnPlotQA>(trg, base, slices));
    qa.emplace_back(std::make_unique<EventQA>(trg, base, slices));


    /* histogram loop ----------------------------------------------- */
    TIter itH(dTrig->GetListOfKeys());
    while (auto* kh = dynamic_cast<TKey*>(itH())) {
      log::trace("Handling [" + trg + "] \"" + string(kh->GetName()) + "\"");
      TObject* obj = kh->ReadObj();
      if (obj->InheritsFrom(TH1::Class()))
        static_cast<TH1*>(obj)->SetDirectory(nullptr);

      ++stat[trg].tot;
      for (auto& m : qa) if (m->process(obj)) { ++stat[trg].used; break; }
    }

    log::ok("Trigger " + trg + ": processed " +
            std::to_string(stat[trg].tot) + " objects");
  }

  csv.close();

  /* -----------------------------  Summary ------------------------- */
  log::banner("Summary");
  std::cout << term::CLR_BOLD
            << std::left  << std::setw(20) << "Trigger"
            << std::right << std::setw(12) << "Total"
            << std::setw(12)                << "Written"
            << term::CLR_RST << "\n";
  for (auto& [t,c] : stat)
    std::cout << std::left  << std::setw(20) << t
              << std::right << std::setw(12) << c.tot
              << std::setw(12)               << c.used << "\n";

  log::ok("All outputs under " + kOutputBase);
    
    /* ------------------------------------------------------------------ *
     *  RUN‑BY‑RUN EVENT SUMMARY  (vertex‑Z entries as proxy for events)   *
     * ------------------------------------------------------------------ */
    {
        using term::CLR_BOLD; using term::CLR_RST;
        const auto& ev = EventQA::eventCounts();
        if (!ev.empty())
        {
            std::vector<std::pair<int,long long>> rows;
            rows.reserve(ev.size());

            for (const auto& [runStr,n] : ev)
                if (std::all_of(runStr.begin(), runStr.end(), ::isdigit))
                    rows.emplace_back(std::stoi(runStr), n);

            std::sort(rows.begin(), rows.end(),
                      [](auto a, auto b){ return a.first < b.first; });

            log::banner("Run‑by‑run event statistics");
            std::cout << CLR_BOLD
                      << std::left  << std::setw(12) << "Run"
                      << std::right << std::setw(15) << "Events"
                      << CLR_RST << "\n";

            long long totalEv = 0;
            for (auto [run,n] : rows) {
                totalEv += n;
                std::cout << std::left  << std::setw(12) << run
                          << std::right << std::setw(15) << n   << "\n";
            }
            std::cout << CLR_BOLD
                      << std::left  << std::setw(12) << "TOTAL"
                      << std::right << std::setw(15) << totalEv
                      << CLR_RST << "\n\n";

            log::ok("Runs analysed : " + std::to_string(rows.size()) +
                    "   |   Total events : " + std::to_string(totalEv));
        }
    }
}


/* ------------------------------------------------------------------ */
/*  MAIN wrapper – fan-out over all runs with a preforked pool        */
/*    optional ‘runFilter’ → process just that one run         */
/* ------------------------------------------------------------------ */
void analyzeRun24or25auau(bool testRun = false, int nSample = -1)
{
    using fs::path;

    /* 0. discover input ROOT files --------------------------------- */
    std::vector<fs::path> runFiles = listRunFiles(kInputDir);
    if (runFiles.empty()) {
        log::err("No input files found in " + kInputDir.string());
        return;
    }

    if (nSample > 0 && static_cast<std::size_t>(nSample) < runFiles.size()) {
        std::mt19937 rng( static_cast<unsigned>(
            std::chrono::steady_clock::now().time_since_epoch().count()) );
        std::shuffle(runFiles.begin(), runFiles.end(), rng);
        runFiles.resize(nSample);
    }

    if (testRun && runFiles.size() > 1)
        runFiles.resize(1);            // old behaviour


//    /* 1. decide pool size (≤ physical cores, ≥ 1) ------------------ */
//    const std::size_t nWorkers =
//        std::min<std::size_t>(runFiles.size(),
//                              std::max<unsigned>(1, std::thread::hardware_concurrency()));
//
//    ROOT::TProcessExecutor exec(nWorkers);        // **preforked** process pool
//
//    log::banner("Launching " + std::to_string(runFiles.size())   +
//                " runs on "    + std::to_string(exec.GetPoolSize()) +
//                " parallel processes");

    /* --- run all QA passes in‑process: no fork, shared statics -------- */
    ROOT::TSequentialExecutor exec;                // single‑process executor

    log::banner("Running " + std::to_string(runFiles.size()) +
                " runs sequentially (shared memory)");
    
    
    /* 2. worker function – receives the *index* (not the path) ----- */
    auto worker = [&](unsigned int idx)->int
    {
        const path& f = runFiles[idx];

        /* extract 8‑digit run number from file name ---------------- */
        std::smatch    m;
        const std::regex reRun(R"(output_([0-9]{8})\.root)");
        const std::string fname = f.filename().string();
        std::regex_search(fname, m, reRun);
        const std::string run = m.empty() ? "Single" : m[1].str();

        const std::string outBase = (kOutputDir / run).string();

        /* pretty progress line – before --------------------------------*/
        log::info(std::string("▶  (") +
                  (idx + 1 < 10 ? " " : "") + std::to_string(idx + 1) + "/" +
                  std::to_string(runFiles.size()) + ")  Run " + run + "  –  start");

        const auto t0 = std::chrono::steady_clock::now();
        runOneQaPass(f.string(), outBase);
        const auto dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        /* pretty progress line – after ---------------------------------*/
        log::ok(std::string("✓  (") +
                (idx + 1 < 10 ? " " : "") + std::to_string(idx + 1) + "/" +
                std::to_string(runFiles.size()) + ")  Run " + run +
                "  –  done in " + std::to_string(dt).substr(0,5) + " s");
        return 0;                       // Map() insists on a return value
    };

    /* 3. fire the jobs – chunkSize = 1  → dynamic load balancing ---- */
    exec.Map(worker, ROOT::TSeqI(runFiles.size()));

    /* 4. optional: merge all runs and re‑run QA on the combined file */
    if (!testRun && runFiles.size() > 1) {
        const path combined = kInputDir / "output_ALL_COMBINED.root";
        log::banner("Hadd – building " + combined.string());

        // ------------------------------------------------------------
        // (a)  print a detailed file list with individual sizes
        // ------------------------------------------------------------
        log::info("Files to be merged (" + std::to_string(runFiles.size()) + " total):");
        std::uintmax_t totBytes = 0;
        for (const auto& f : runFiles) {
            const auto sz = fs::file_size(f);
            totBytes += sz;
            log::info("   + " + f.filename().string() +
                      "  (" + std::to_string(sz / 1'024'000) + " MB)");
        }
        log::info("   ------------------------------------------------");
        log::info("Accumulated input size : " +
                  std::to_string(totBytes / 1'024'000) + " MB");

        // ------------------------------------------------------------
        // (b)  run TFileMerger with progress timing
        // ------------------------------------------------------------
        const auto t0Hadd = std::chrono::steady_clock::now();

        TFileMerger merger(/*dryRun=*/false, /*verbose=*/true);
        merger.OutputFile(combined.c_str(), "RECREATE");

        for (const auto& f : runFiles) {
            log::trace("TFileMerger  ← adding  " + f.string());
            merger.AddFile(f.c_str());
        }

        if (!merger.Merge()) {
            log::err("TFileMerger failed – combined QA skipped");
            return;
        }

        const auto dHadd = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - t0Hadd).count();
        const auto outSize = fs::file_size(combined);

        log::ok("Combined ROOT file created in " +
                std::to_string(dHadd).substr(0,5) + " s,  size " +
                std::to_string(outSize / 1'024'000) + " MB");

        // ------------------------------------------------------------
        // (c)  re‑run the QA pass on the freshly merged file
        // ------------------------------------------------------------
        runOneQaPass(combined.string(),
                     (kOutputDir / "Combined").string());
    }
}
