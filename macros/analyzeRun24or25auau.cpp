// analyzeRun24or25auau.cpp  – ROOT ≥ 6, C++17
// ===============================================================
//  • Pass‑0: catalogue every histogram in the file  (console + .txt)
//  • Pass‑1: modular QA (EMCal, HCal, sEPD, MBD, correlations,
//            π0 invariant‑mass spectra, …) with automatic centrality
//            slice replication and North/South map fusion.
//  • Built‑in [TRACE] instrumentation to pinpoint run‑time crashes.
//  • Global toggle kDoPi0Fit to switch π0 mass fitting on/off
// ===============================================================

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <thread>
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
/*  <<<   π0‐fit master switch   >>>                                         *
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
  base /= (slice == "Inclusive" ? "noCentralityDep"
                                : "Cent_" + slice);
  return base / sub;
}



// ╔══════════════════════════════════════════════╗
// ║ 4.  π0 CUT KEY / PEAK FITTER                 ║
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
private:
  struct H{ size_t operator()(const string& s)const noexcept{ return std::hash<string>{}(s);} };
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
// ║     π0   I N V A R I A N T ‑ M A S S   QA    ║
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
    runID = root.parent_path().filename().string();  // "00066013", "Combined", …

    /* new CSV for S/B */
    fs::path p = root / "EMCal/pi0QA" / "Pi0SignalBackground.csv";
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
    if (n.rfind("mInv_",0)!=0) return false;          // not a π0 spectrum

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

    fs::path subDir = "EMCal/pi0QA";
    subDir /= combDir;                  // …/pi0QA/<cut>/[…]
    if(!pTInt) subDir/=("pT_"+sf3(ck.pLo)+"_to_"+sf3(ck.pHi));
    fs::path outPng = cPath(root,slice,subDir)/(n+".png");
    ensure_dir(outPng.parent_path());

    //----------------------------------------------------------------
    // 1.  Robust π0 peak search  → initial μ, A
    //----------------------------------------------------------------
    TH1* h = static_cast<TH1*>(o);
    const double piFitLo = 0.05, piFitHi = 0.35;            // π0 window
    const int    iLoPi   = binAt(h,piFitLo),  iHiPi = binAt(h,piFitHi);

    int iMaxPi = iLoPi;
    double maxCntPi = 0;
    for(int i=iLoPi;i<=iHiPi;++i)
      if(h->GetBinContent(i)>maxCntPi){ maxCntPi=h->GetBinContent(i); iMaxPi=i; }

    const double piMu0     = h->GetBinCenter(iMaxPi);       // ~peak
    const double piAmp0    = maxCntPi;
    const double piSigma0  = 0.025;

    //----------------------------------------------------------------
    // 2.  Composite fit function  π0‑Gaus + poly‑2 background
    //----------------------------------------------------------------
    TF1 total("total","gaus(0)+pol2(3)",piFitLo,piFitHi);
    total.SetParNames("A","mu","sigma","c0","c1","c2");

    total.SetParameters(piAmp0, piMu0, 0.022,    // narrower σ start
                        1,       0,     0);      // flat background
    total.SetParLimits(0,   0,  1e9);            // A  ≥ 0
    total.SetParLimits(1,   0.10,  0.17);        // μ  within window
    total.SetParLimits(2,   0.010, 0.060);       // σ  sensible range

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

    /* -------------- π0 parameters -------------------------------- */
    double piMu=piMu0, piSig=piSigma0,
           piMuErr=0,   piSigErr=0;
    if(fitOK){
        piMu     = total.GetParameter(1);
        piSig    = total.GetParameter(2);
        piMuErr  = total.GetParError (1);
        piSigErr = total.GetParError (2);
    }

    //----------------------------------------------------------------
    // 3.  Fixed‑background η peak search (0.45 – 0.80 GeV)
    //----------------------------------------------------------------
    double etaMu = 0, etaSig = 0, etaMuErr = 0, etaSigErr = 0;

    if (runID == "Combined")          // ← skip per‑run files
    {
          const double etaLo = 0.45, etaHi = 0.80;

          /* take the poly coefficients from the π0 fit (good first‑order BKG) */
          TF1 bkg("bkg","pol2",etaLo,etaHi);
          bkg.SetParameters(total.GetParameter(3),
                            total.GetParameter(4),
                            total.GetParameter(5));

          /* background‑subtracted max gives seed */
          int iLoEta = binAt(h, etaLo),  iHiEta = binAt(h, etaHi);
          int iMaxEta = iLoEta;  double maxCntEta = 0.;
          for (int i = iLoEta; i <= iHiEta; ++i)
          {
              double y = h->GetBinContent(i) -
                         std::max(bkg.Eval(h->GetBinCenter(i)), 0.);
              if (y > maxCntEta) { maxCntEta = y; iMaxEta = i; }
          }

          if (maxCntEta > 0)            // peak exists
          {
              etaMu  = h->GetBinCenter(iMaxEta);
              etaSig = 0.035;                          // rough width

              TF1 gEta("gEta","gaus(0)+pol2(3)",etaLo,etaHi);
              gEta.SetParNames("A","#mu","#sigma","c0","c1","c2");
              double ampEta = maxCntEta;
              gEta.SetParameters(ampEta, etaMu, etaSig,
                                 bkg.GetParameter(0),
                                 bkg.GetParameter(1),
                                 bkg.GetParameter(2));

              gEta.SetParLimits(0,      0,   1e9);
              gEta.SetParLimits(1,   etaLo,  etaHi);
              gEta.SetParLimits(2,   0.015,  0.080);

              /* fix background so the η‑Gaussian is stable */
              gEta.FixParameter(3, bkg.GetParameter(0));
              gEta.FixParameter(4, bkg.GetParameter(1));
              gEta.FixParameter(5, bkg.GetParameter(2));

              bool etaOK = (h->Fit(&gEta, "QRN0") == 0);
              if (etaOK)
              {
                  etaMu     = gEta.GetParameter(1);
                  etaMuErr  = gEta.GetParError (1);
                  etaSig    = gEta.GetParameter(2);
                  etaSigErr = gEta.GetParError (2);

                  /* store η fit once per slice for later overlays */
                  if (_storedEtaFit.count(slice) == 0)
                      _storedEtaFit[slice].reset(new TF1(gEta));
              }
              else
              {
                  etaMu = etaSig = etaMuErr = etaSigErr = 0;   // failed fit
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
    // 5.  Signal / Background CSV  (unchanged for π0)
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

      TLegend leg(0.55,0.64,0.88,0.88);
      leg.SetBorderSize(0);
      leg.SetTextAlign(12);                         // left‑align text

      /* π0: write μ‑line and σ‑line underneath one another */
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

    /* --- keep a copy of the π0 fit (once per slice) --------------- */
    if(fitOK && _storedFit.count(slice)==0){
        _storedFit[slice].total.reset(new TF1(total));
        _storedFit[slice].poly .reset(new TF1(poly ));
    }

    /* --- store run‑summary point (Inclusive, pT‑integrated) ------- */
    if(fitOK && pTInt && slice=="Inclusive" && s_runPoints.count(runID)==0){
        s_runPoints[runID] = {piMu, piMuErr, piSig, piSigErr};
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

    if(pTIntegrated)
      _centralHists.emplace(slice,cl);
  }

  //---------------- write 2 × 3 overview + μ,σ vs centrality ---------
  void writeSummaryPanels()
  {
    if(_centralHists.empty()) return;

    /* ---------- (A) 2 × 3 mass spectra panel -------------------- */
    TCanvas cGrid("c_pi0Cent","π0 – all centralities",1800,1000);
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

        /* ------------------------------------------------------------------ *
         *  Smart text placement:                                              *
         *     – If the peak sits on the left‑hand side (μ ≲ 0.22 GeV) the     *
         *        annotation is moved to the right (x ≈ 0.60) to avoid overlap *
         *     – Otherwise it stays on the default left margin (x ≈ 0.13).     *
         * ------------------------------------------------------------------ */
        const double xText = (mu < 0.22) ? 0.60 : 0.13;   // 0.22 GeV ~ middle of canvas

        TLatex tx;  tx.SetNDC();  tx.SetTextSize(0.04);
        tx.DrawLatex(xText, 0.86, lbl.c_str());
        tx.DrawLatex(xText, 0.78,
                     Form("#mu = %.3f #pm %.3f GeV",  mu,  emu));
        tx.DrawLatex(xText, 0.70,
                     Form("#sigma = %.3f #pm %.3f GeV", si, esi));


      if(sl!="Inclusive"){
          int lo=std::stoi(sl.substr(0,sl.find('_')));
          int hi=std::stoi(sl.substr(sl.find('_')+1));
          vC.push_back(0.5*(lo+hi));   vCerr.push_back(0.5*(hi-lo));
          vMu.push_back(mu); vMuErr.push_back(emu);
          vSi.push_back(si); vSiErr.push_back(esi);
      }
    }

    fs::path pngGrid = root/"EMCal/pi0QA"/cutTag/"Pi0Mass_AllCentrality.png";
    ensure_dir(pngGrid.parent_path()); cGrid.SaveAs(pngGrid.string().c_str());

    /* ---------- (B) μ,σ versus centrality (π0 only, unchanged) --- */
    if(!vC.empty()){
      int n=vC.size();
      auto gMu = std::make_unique<TGraphErrors>(n,
                          vC.data(), vMu.data(), nullptr, vMuErr.data());   // no x‑errors
      auto gSi = std::make_unique<TGraphErrors>(n,
                          vC.data(), vSi.data(), nullptr, vSiErr.data());   // no x‑errors
      gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
      gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

      TCanvas cGS("c_mu_sigma_vs_cent",
                  "π^{0} peak position / width vs centrality",800,800);

      TPad *p1=new TPad("p1","",0,0.35,1,1);
      p1->SetBottomMargin(0.02); p1->Draw(); p1->cd();
      gMu->SetTitle("π^{0} mass versus centrality;Centrality [%];m_{π^{0}} (GeV/c^{2})");
      gMu->Draw("AP");

      cGS.cd();
      TPad *p2=new TPad("p2","",0,0,1,0.32);
      p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30);
      p2->Draw(); p2->cd();
      gSi->SetTitle(";Centrality [%];σ_{π^{0}} (GeV/c^{2})");
      gSi->Draw("AP");

      fs::path pngGraph = root/"EMCal/pi0QA"/cutTag/"Pi0Mass_Sigma_vs_Centrality.png";

      cGS.SaveAs(pngGraph.string().c_str());
    }
  }

  //---------------- write run‑by‑run μ,σ summary (π0 only) ----------
  void writeRunSummary()
  {
      if (runID != "Combined" || s_runPoints.size() < 2) return;

      static std::unordered_set<std::string> done;   // one PNG per cut
      if (done.count(cutTag)) return;
      done.insert(cutTag);

      /* convert maps → sorted vectors (by run number) */
      std::vector<int> runs;
      for(const auto& [r,_]:s_runPoints)
          if(std::all_of(r.begin(),r.end(),::isdigit))
              runs.push_back(std::stoi(r));
      if(runs.empty()) return;
      std::sort(runs.begin(),runs.end());

      const int n=runs.size();
      std::vector<double> x(n),yMu(n),eMu(n),ySi(n),eSi(n);
      for(int i=0;i<n;++i){
          const auto& p=s_runPoints[std::to_string(runs[i])];
          x[i]=runs[i]; yMu[i]=p.mu; eMu[i]=p.muErr;
          ySi[i]=p.sigma; eSi[i]=p.sigmaErr;
      }
      auto gMu=new TGraphErrors(n,x.data(),yMu.data(),nullptr,eMu.data());
      auto gSi=new TGraphErrors(n,x.data(),ySi.data(),nullptr,eSi.data());
      gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
      gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

      TCanvas cR("c_mu_sigma_vs_run",
                 "π^{0} peak position / width vs run",900,800);

      TPad *p1=new TPad("p1","",0,0.35,1,1);
      p1->SetBottomMargin(0.02); p1->Draw(); p1->cd();
      gMu->SetTitle("π^{0} mass versus run;Run number;m_{π^{0}} (GeV/c^{2})");
      gMu->Draw("AP");

      cR.cd();
      TPad *p2=new TPad("p2","",0,0,1,0.32);
      p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30);
      p2->Draw(); p2->cd();
      gSi->SetTitle(";Run number;σ_{π^{0}} (GeV/c^{2})");
      gSi->Draw("AP");

      fs::path pngRun = root.parent_path() / (cutTag + "_Pi0Mass_Sigma_vs_Run.png");
      cR.SaveAs(pngRun.string().c_str());
  }

  /* ---------------------------------------------------------------- */
  /*  data members                                                    */
  /* ---------------------------------------------------------------- */
  std::ofstream&                                   csvFit;
  std::ofstream                                    csvSB;

  std::unordered_map<std::string, TH1*>            _centralHists;

  struct FitPair {
      std::unique_ptr<TF1> total;
      std::unique_ptr<TF1> poly;
  };
  std::unordered_map<std::string, FitPair>         _storedFit;
  std::unordered_map<std::string, std::unique_ptr<TF1>> _storedEtaFit;

  std::unordered_map<std::string,
                     std::vector<TH1*>>            _ptHists;
  std::unordered_map<std::string, FitInfo>         _fitSummary;

  /* run label of this instance */
  std::string runID;
  /* directory tag that identifies one (E , χ² , asym) cut‑combination */
  std::string cutTag;

  /* ---------- static: accumulate π0 points over all runs ---------- */
  struct RunPoint { double mu, muErr, sigma, sigmaErr; };
  static inline std::unordered_map<std::string, RunPoint> s_runPoints;
  static inline bool s_summaryWritten=false;
};



// ——— Detector–detector correlations ————————————
class CorrQA : public QA
{
public:
    using QA::QA;                               // inherit ctors

    // =====================   one histogram   ==========================
    bool process(TObject* o) override
    {
        if (!o->InheritsFrom(TH2::Class()))          return false;
        const std::string hName = o->GetName();
        if (hName.find("_vs_") == std::string::npos) return false;   // not a correlation

        const std::string slice   = sliceKey(hName);                 // Inclusive / Cent_…
        const bool        isSEPD  = (hName.find("SEPD") != std::string::npos ||
                                     hName.find("sEPD") != std::string::npos);
        const std::string subDir  = isSEPD ? "sEPD" : "Correlations";

        // ───────────────── tidy ranges & apply consistent style ─────────────────
        TH2* h2 = static_cast<TH2*>(o);
        tidyAxes (h2);
        styleAxes(h2, /*smallPad=*/false);

        // ───────────────── per‑run PNG  ───────────────────────────────
        {
            fs::path out = cPath(root, slice, subDir) / (h2->GetName() + std::string(".png"));
            ensure_dir(out.parent_path());

            TCanvas c("c_corr", "", 1100, 800);
            setupPad(&c);
            h2->Draw("COLZ");
            c.SaveAs(out.string().c_str());
        }

        // ───────────────── cache clone for summary sheet ──────────────
        const std::string runID = root.parent_path().filename().string();
        if (runID != "Combined")
        {
            auto* clone = static_cast<TH2*>(o->Clone());
            clone->SetDirectory(nullptr);
            tidyAxes (clone);
            styleAxes(clone, /*smallPad=*/true);   // default for summary grid
            s_cache[subDir][hName][runID].reset(clone);
        }
        return true;
    }

    // =====================   final summary   ==========================
    ~CorrQA() override
    {
        const std::string runID = root.parent_path().filename().string();
        if (runID != "Combined" || s_done) return;
        s_done = true;

        const fs::path baseOut = root.parent_path();   // “…/Combined/<trigger>”

        for (const auto& [subDir, byName] : s_cache)
            for (const auto& [hName, byRun] : byName)
            {
                if (byRun.empty()) continue;

                /* gather <run,hist> sorted by run number (as string) */
                std::vector<std::pair<std::string,TH2*>> runs;
                runs.reserve(byRun.size());
                for (const auto& [r,h] : byRun) runs.emplace_back(r, h.get());
                std::sort(runs.begin(), runs.end(),
                           [](auto& a, auto& b){ return a.first < b.first; });

                /* common log‑safe colour‑scale */
                double zMax = 0;
                for (const auto& [_,h] : runs) zMax = std::max(zMax, h->GetMaximum());
                for (const auto& [_,h] : runs) { h->SetMinimum(1); h->SetMaximum(zMax); }

                /* paginated 10×10 grids */
                const int perPage = 100;
                const int nPages  = (runs.size() + perPage - 1) / perPage;
                fs::path outDir   = baseOut / subDir / ("summary_" + hName);
                ensure_dir(outDir);

                for (int pg = 0; pg < nPages; ++pg)
                {
                    const int lo =  pg      * perPage;
                    const int hi = (pg + 1) * perPage;

                    TCanvas c(Form("c_%s_p%02d", hName.c_str(), pg+1), "", 3000, 3000);
                    c.Divide(10,10,0.000,0.000);

                    for (int i = lo; i < hi && i < (int)runs.size(); ++i)
                    {
                        c.cd(i - lo + 1);
                        setupPad(gPad);                       // pad‑specific margins/log‑Z
                        runs[i].second->Draw("COLZ");

                        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
                        tl.DrawLatex(0.02, 0.90, runs[i].first.c_str());
                    }
                    c.SaveAs((outDir / Form("page%02d.png", pg+1)).string().c_str());
                }
            }
    }

    // =====================   helpers   ================================
    static void tidyAxes(TH2* h)
    {
        int fx = 1, lx = h->GetNbinsX();
        while (fx <= lx && h->Integral(fx,fx,1,h->GetNbinsY()) == 0) ++fx;
        while (lx >= fx && h->Integral(lx,lx,1,h->GetNbinsY()) == 0) --lx;
        if (fx < lx) h->GetXaxis()->SetRange(fx,lx);

        int fy = 1, ly = h->GetNbinsY();
        while (fy <= ly && h->Integral(1,h->GetNbinsX(),fy,fy) == 0) ++fy;
        while (ly >= fy && h->Integral(1,h->GetNbinsX(),ly,ly) == 0) --ly;
        if (fy < ly) h->GetYaxis()->SetRange(fy,ly);
    }

    static void styleAxes(TH2* h, bool smallPad)
    {
        const double labSize = smallPad ? 0.028 : 0.030;
        const double titSize = smallPad ? 0.034 : 0.037;

        h->GetXaxis()->SetLabelSize(labSize);
        h->GetYaxis()->SetLabelSize(labSize);
        h->GetZaxis()->SetLabelSize(labSize);

        h->GetXaxis()->SetTitleSize(titSize);
        h->GetYaxis()->SetTitleSize(titSize);
        h->GetZaxis()->SetTitleSize(titSize);

        h->GetXaxis()->SetTitleOffset(1.15);
        h->GetYaxis()->SetTitleOffset(1.50);
        h->GetZaxis()->SetTitleOffset(1.20);

        // a subtle grid can be helpful but is optional:
        // h->SetContour(99);  // already default in many setups
    }

    /** common pad cosmetics: margins + log‑Z **/
    static void setupPad(TVirtualPad* pad)
    {
        pad->SetLogz();
        pad->SetRightMargin(0.16);
        pad->SetLeftMargin (0.12);
        pad->SetBottomMargin(0.12);
        pad->SetTopMargin  (0.06);
    }

    /* subDir → histName → runID → TH2 clone */
    using RunMap  = std::unordered_map<std::string,std::unique_ptr<TH2>>;
    using NameMap = std::unordered_map<std::string,RunMap>;
    static inline std::unordered_map<std::string,NameMap> s_cache;
    static inline bool s_done = false;
};



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
    for (const auto& sl : slices) {
      auto it = _centralMaps.find(sl);
      if (it == _centralMaps.end()) continue;

      c.cd(pad++);
      it->second->Draw("COLZ");

      TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
      tl.DrawLatex(0.15, 0.85,
                   (sl == "Inclusive" ? "Inclusive"
                                      : ("Cent " + sl)).c_str());
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

      /* 3)  save individual PNG (1 bin ⇔ 1 pixel) ---------------------- */
      const int cw = (kEMCalCanvasW > 0) ? kEMCalCanvasW : nEta * px;
      const int ch = (kEMCalCanvasH > 0) ? kEMCalCanvasH : nPhi * px;

      fs::path outPng = cPath(root, slice, "EMCal")
                      / (src->GetName() + std::string(".png"));
      ensure_dir(outPng.parent_path());

      TCanvas c("c_emcal", "", cw, ch);
      c.SetRightMargin (0.17);
      c.SetFixedAspectRatio(false);           // allow resizing in a viewer
      c.SetLeftMargin  (0.14);
      c.SetBottomMargin(0.08);
      c.SetTopMargin   (0.04);
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
 public: using QA::QA;

  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string n = o->GetName();
    const bool isI = n.rfind("h_IHCAL_",0)==0;
    const bool isO = n.rfind("h_OHCAL_",0)==0;
    if (!isI && !isO) return false;

    const std::string slice   = sliceKey(n);
    const bool        isMap   = (n.find("_EtaPhiMap_")!=std::string::npos);

    auto makePanel = [&](TH2* src, const fs::path& outPng)
    {
      /* 0) constants -------------------------------------------------- */
      constexpr int nPhi = 64;                      // rows     (Y)
      constexpr int nEta = 24;                      // columns  (X)
      constexpr int px   = 18;                      // pixel‑size (PNG)

      /* 1) local clone + bad‑plate masking ---------------------------- */
      std::unique_ptr<TH2> h(static_cast<TH2*>(src->Clone()));
      h->SetDirectory(nullptr); h->SetStats(0); h->SetContour(99);

      for (int ip = 0; ip < nPhi; ++ip)
        for (int ie = 0; ie < nEta; ++ie)
          if (isBadHcalPlate(hcal_sector_from_idx(ie,ip),
                             hcal_plate_from_idx(ie,ip)))
            h->SetBinContent(h->FindBin(ip,ie),-9999.);   // white holes

      /* 2) rotate  (η → X,  φ → Y) ----------------------------------- */
      std::unique_ptr<TH2F> rot(new TH2F(("hRot_"+std::string(src->GetName())).c_str(),
                                         h->GetTitle(),
                                         nEta, 0, nEta,         // X‑bins
                                         nPhi, 0, nPhi));       // Y‑bins
      for (int ip = 1; ip <= nPhi; ++ip)
        for (int ie = 1; ie <= nEta; ++ie)
          rot->SetBinContent(ie, ip, h->GetBinContent(ip, ie));

      rot->SetDirectory(nullptr);
      rot->SetMinimum(1.);
      rot->GetXaxis()->SetTitle("#eta index");
      rot->GetYaxis()->SetTitle("#phi index");

      /* --- make tick‑label fonts smaller & tidy titles ---------------- */
      const double kLabSize = 0.025;   // tick‑label font (default ≈ 0.04)
      const double kTitSize = 0.030;   // axis‑title font  (default ≈ 0.04)
      rot->GetXaxis()->SetLabelSize(kLabSize);
      rot->GetYaxis()->SetLabelSize(kLabSize);
      rot->GetZaxis()->SetLabelSize(kLabSize);   // palette scale
      rot->GetXaxis()->SetTitleSize(kTitSize);
      rot->GetYaxis()->SetTitleSize(kTitSize);
      rot->GetZaxis()->SetTitleSize(kTitSize);

      rot->GetXaxis()->SetNdivisions(nEta, kFALSE);
      rot->GetYaxis()->SetNdivisions(nPhi, kFALSE);

      /* --- canvas size: use user constant unless ≤0 ------------------ */
      const int cw = (kHCalCanvasW > 0) ? kHCalCanvasW : nEta * px;
      const int ch = (kHCalCanvasH > 0) ? kHCalCanvasH : nPhi * px;

      TCanvas c("c_hcal","",cw,ch);
      c.SetRightMargin(0.16);

      c.SetFixedAspectRatio(false);
      c.SetRightMargin(0.16);
      c.SetLeftMargin (0.08);
      c.SetBottomMargin(0.08);
      c.SetTopMargin  (0.04);
      c.SetFixedAspectRatio();                  // 1 bin ⇔ 1 pixel

      rot->Draw("COLZ");

      /* 4) grid ------------------------------------------------------- */
      TLine l; l.SetLineColor(kBlack);

      // thin: individual towers
      l.SetLineWidth(1);
      for (int x=1; x<nEta; ++x) l.DrawLine(x,0,x,nPhi);
      for (int y=1; y<nPhi; ++y) l.DrawLine(0,y,nEta,y);

      // thick: sector & plate structure
      l.SetLineWidth(4);
      for (int y=0; y<=nPhi; y+=2) l.DrawLine(0,y,nEta,y);     // 32 sectors
      l.DrawLine( 8,0, 8,nPhi);                                // inner‑boards
      l.DrawLine(16,0,16,nPhi);

      /* 5) save ------------------------------------------------------- */
      ensure_dir(outPng.parent_path());
      c.SaveAs(outPng.string().c_str());
    };
    /* ----------------------------------------------------------------- */

    // ---------- write outputs -----------------------------------------
    fs::path subDir = isI ? "IHCal" : "OHCal";
    fs::path out    = cPath(root, slice, subDir) / (n + ".png");

    if (isMap && o->InheritsFrom(TH2::Class()))
      makePanel(static_cast<TH2*>(o), out);
    else if (o->InheritsFrom(TH2::Class()))
      save2D(static_cast<TH2*>(o), out);
    else
      save1D(static_cast<TH1*>(o), out);

    return true;
  }
};


// ——— sEPD & MBD QA (NS combiner) ——————————————
template<class DERIVED> class NSDetectorQA : public QA{
public:
  NSDetectorQA(string t,fs::path b,const CentList& s,NSCache<MapPair>& c):
      QA(t,b,s),cache(c){}
  bool process(TObject* o) override
  {
    if(!o->InheritsFrom(TH1::Class())) return false;
    string n=o->GetName(); if(!DERIVED::accept(n)) return false;
    string sl=sliceKey(n);
    bool isSouth=n.find("_South_")!=string::npos;

    // A. arrival trace
    log::trace(string(DERIVED::subdir)+"  slice="+sl+
               "  hist=\""+n+"\"  south? "+(isSouth?"yes":"no"));

    /* scalar spectra */
    if(!o->InheritsFrom(TH2::Class())){
      auto save=[&](const string& slice){
        fs::path out=cPath(root,slice,DERIVED::subdir)/(n+".png");
        save1D(static_cast<TH1*>(o),out);
      };
      save(sl);
      return true;
    }

    // B. cache bookkeeping
    MapPair& mp=cache[trig+sl];
    TH2* clone=static_cast<TH2*>(o->Clone());
    clone->SetDirectory(nullptr); clone->SetStats(0);

    log::trace("   clone ptr = "+std::to_string((uintptr_t)clone));

    isSouth ? mp.s=clone : mp.n=clone;

    log::trace("   cache state  n="+std::to_string((uintptr_t)mp.n)+
               "  s="+std::to_string((uintptr_t)mp.s));

    if(!cache.ready(trig+sl)) return true;

    // C. both maps present → draw
    MapPair out=cache.pop(trig+sl);

    log::trace("   >> drawing combined canvas for slice "+sl+
               "  n="+std::to_string((uintptr_t)out.n)+
               "  s="+std::to_string((uintptr_t)out.s));

    // -----------------------------------------------------------------
    // ONE finished South–North canvas
    // -----------------------------------------------------------------
    auto save = [&](const std::string& slice)
      {
        fs::path outPng = cPath(root, slice, DERIVED::subdir)
                          / (DERIVED::fileName(trig) + ".png");

        /* -------- 1)  prepare uniform colour scale (log‑friendly) ------ */
        const double zMin = 0.;                                     // >0 for log
        const double zMax = std::max(out.s->GetMaximum(),
                                     out.n->GetMaximum());
        out.s->SetMinimum(zMin);  out.n->SetMinimum(zMin);
        out.s->SetMaximum(zMax);  out.n->SetMaximum(zMax);

        /* -------- 2)  canvas & two pads -------------------------------- */
        TCanvas c("c_mbd", "", 1200, 600);
        c.Divide(2, 1, 0.01, 0.01);

        auto drawPad = [&](TH2* h, const char* ttl)
        {
          gPad->SetRightMargin(0.20);          // leave room for palette + labels
          gPad->SetLeftMargin (0.10);
          gPad->SetBottomMargin(0.10);
          gPad->SetTopMargin  (0.08);

          h->SetTitle(ttl);
          h->GetZaxis()->SetTitle("Counts");
          h->GetZaxis()->SetTitleOffset(1.3);  // pull title away from palette
          h->Draw("POLZ");                     // hexagons + palette
        };

        c.cd(1); drawPad(out.s, DERIVED::titleSouth);
        c.cd(2); drawPad(out.n, DERIVED::titleNorth);

        ensure_dir(outPng.parent_path());
        c.SaveAs(outPng.string().c_str());
    };

    save(sl);
    return true;
  }
protected:
  NSCache<MapPair>& cache;
};

                                                                 
struct MBDTag{
  /* recognise every histogram that clearly belongs to the MBD detector
   * – hit‑maps (“h_MBD_…”) as well as charge & Q‑sum spectra           */
  static bool accept(const std::string& s)
  {
    return  s.rfind("h_MBD_"      , 0) == 0   ||   // hit‑maps
            s.rfind("h_charge_MBD", 0) == 0   ||   // per‑event ΣQ spectra
            s.rfind("h_Qsum_MBD"  , 0) == 0;        // centrality helper
  }
  static constexpr const char* subdir = "MBD";
  static std::string fileName(const std::string& t){ return "MBD_Hitmap_NS_" + t; }
  static constexpr const char* titleSouth = "MBD South";
  static constexpr const char* titleNorth = "MBD North";
};

struct sEPDTag{
  static bool accept(const std::string& s){
        /* accept both “sEPD” (hit‑maps, Q‑spectra, …) **and**
           “SEPD” (correlation histograms, etc.)               */
        return s.find("sEPD") != std::string::npos ||
               s.find("SEPD") != std::string::npos;
  }
  static constexpr const char* subdir="sEPD";
  static string fileName(const string& t){return "sEPD_Hitmap_NS_"+t;}
  static constexpr const char* titleSouth="sEPD South";
  static constexpr const char* titleNorth="sEPD North";
};

using MbdQA  = NSDetectorQA<MBDTag>;
using SepdQA = NSDetectorQA<sEPDTag>;

// ------------------------------------------------------------------
//  Jet‑QA module  –  handles 1‑D, 2‑D and 3‑D jet histograms
// ------------------------------------------------------------------
class JetQA : public QA
{
 public:
  using QA::QA;

  // == main entry ====================================================
  bool process(TObject* o) override
  {
    // accept only the three jet‑QA families --------------------------
    const std::string n = o->GetName();
    const bool is1D = (n.rfind("h_maxJetEt_"          ,0) == 0);
    const bool is2D = (n.rfind("h_leadEt_vs_subEt_"   ,0) == 0);
    const bool is3D = (n.rfind("h_jetEt_area_nConst_",0) == 0);
    if (!is1D && !is2D && !is3D) return false;

    const std::string slice = sliceKey(n);            // Inclusive / Cent_…
    const std::string rLab  = radiusTag(n);           // r02 / r04 / …

    fs::path outBase = cPath(root, slice, fs::path("jetQA") / rLab);

    if (is1D && o->InheritsFrom(TH1::Class()))
    {
      save1D(static_cast<TH1*>(o), outBase / (n + ".png"));
      return true;
    }
    if (is2D && o->InheritsFrom(TH2::Class()))
      return handle2D(static_cast<TH2*>(o), outBase, n);

    if (is3D && o->InheritsFrom(TH3::Class()))
      return handle3D(static_cast<TH3*>(o), outBase, n);

    return false;                         // should never reach here
  }

 private:
  // ---- helper: radius tag (“r02”, “r04”, …) ------------------------
  static std::string radiusTag(const std::string& hname)
  {
    std::smatch m; std::regex r(R"(_(r[0-9]+|R[0-9]+)_)");
    return std::regex_search(hname,m,r) ? m[1].str() : "UnknownR";
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
};

class VnPlotQA : public QA
{
 public:
    /* ------------------------------------------------------------------ *
     *  ctor – standard QA signature (+ auto vn sub‑folder)               *
     * ------------------------------------------------------------------ */
  VnPlotQA(std::string           trig,
             std::filesystem::path base,
             const CentList&       slices)
    : QA(std::move(trig), base, slices)     // initialise the QA base‑class
    , _outDir(base / "vn")                  // all vₙ plots go into …/<trigger>/vn
  {
      if (!std::filesystem::exists(_outDir))
        std::filesystem::create_directories(_outDir);
  }

  /* ------------------------------------------------------------------ *
   *  called once for *every* TObject in the input ROOT file            *
   * ------------------------------------------------------------------ */
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TProfile::Class())) return false;

    const std::string hname = o->GetName();
    /*  pattern:  p_v< n >_< DET >_< lo >_< hi >_< trigger >            */
    std::smatch m;
    static const std::regex re(
        R"(p_v([123])_([A-Za-z0-9_]+)_([0-9]+)_([0-9]+)_(.+))");

    if (!std::regex_match(hname, m, re)) return false;

    const int    nHarm = std::stoi(m[1]);               // 1,2,3
    const string det   = m[2];
    const string cent  = m[3].str() + '_' + m[4].str(); // "20_40"
    const string trig  = m[5];

    /* trigger → detector → harmonic → centrality → list<TProfile*> */
    _cache[trig][det][nHarm][cent]
        .push_back(static_cast<TProfile*>(o));
    return true;
  }

  /* ------------------------------------------------------------------ *
   *  destructor = single final plotting pass                           *
   * ------------------------------------------------------------------ */
  ~VnPlotQA() override { writeCanvases(); }

 private:
  using ProfileVec = std::vector<TProfile*>;

  /* ---- helper: convert one TProfile → TGraphErrors ----------------- */
  static std::unique_ptr<TGraphErrors> makeGraph(const TProfile* p)
  {
    const int nb = p->GetNbinsX();
    auto g = std::make_unique<TGraphErrors>(nb);

    for (int i = 1; i <= nb; ++i)
    {
      const double xLo = p->GetXaxis()->GetBinLowEdge(i);
      const double xHi = p->GetXaxis()->GetBinUpEdge (i);
      const double x   = 0.5 * (xLo + xHi);
      const double ex  = 0.5 * (xHi - xLo);

      const double y   = p->GetBinContent(i);
      const double ey  = p->GetBinError  (i);

      g->SetPoint     (i - 1, x,  y);
      g->SetPointError(i - 1, ex, ey);
    }
    g->SetLineWidth(2);
    g->SetMarkerStyle(kFullCircle);
    return g;
  }

  /* ---- heavy part: generate canvases and write PNGs ---------------- */
  void writeCanvases()
  {
    for (const auto& [trig, detMap] : _cache)
      for (const auto& [det, harmMap] : detMap)
        for (const auto& [n, centMap] : harmMap)          // n = 1,2,3
        {
          //----------------------------------------------------------------
          // (A) one canvas: *all* centralities, fixed {detector,n,trigger}
          //----------------------------------------------------------------
          {
            const std::string ttl = Form("v_{%d} vs p_{T} – %s (%s trigger)",
                                         n, det.c_str(), trig.c_str());

            TCanvas c(Form("c_v%d_%s_%s", n, det.c_str(), trig.c_str()),
                      "", 1200, 900);
            c.SetGrid();

            TLegend leg(0.15, 0.70, 0.45, 0.88); leg.SetBorderSize(0);

            int    col  = 1;
            double yMax = 0.;

            for (const auto& [cent, profList] : centMap)
            {
              if (profList.empty()) continue;
              auto g = makeGraph(profList.front());
              g->SetLineColor(col);
              g->SetMarkerColor(col);
              g->SetTitle(ttl.c_str());

              g->Draw(col == 1 ? "APL" : "PL SAME");
              leg.AddEntry(g.get(), Form("Cent %s %%", cent.c_str()), "pl");

              const int    npts = g->GetN();
              const double localMax =
                  TMath::MaxElement(npts, g->GetY());
              yMax = std::max(yMax, localMax);

              _ownedGraphs.push_back(std::move(g));
              ++col;
            }
              if (yMax > 0.) {
                  c.Update();                                            // ensure the frame exists
                  auto *frame = static_cast<TH1*>(c.GetPrimitive("htemp"));   // implicit axis frame
                  if (frame) frame->SetMaximum(1.15 * yMax);
              }

            leg.Draw();
            saveCanvas(c, det, Form("v%d_%s_allCent_%s.png",
                                    n, det.c_str(), trig.c_str()));
          }

          //----------------------------------------------------------------
          // (B) one canvas *per* centrality: compare detectors
          //----------------------------------------------------------------
          for (const auto& [centWanted, _junk] : centMap)
          {
            const std::string ttl =
                Form("v_{%d} vs p_{T} – Cent %s %% (%s trigger)",
                     n, centWanted.c_str(), trig.c_str());

            TCanvas c(Form("c_v%d_cent%s_%s", n,
                           centWanted.c_str(), trig.c_str()),
                      "", 1200, 900);
            c.SetGrid();

            TLegend leg(0.15, 0.70, 0.45, 0.88); leg.SetBorderSize(0);

            int col = 1;
            for (const auto& [det2, harmMap2] : detMap)
            {
              auto itH = harmMap2.find(n);
              if (itH == harmMap2.end()) continue;

              auto itC = itH->second.find(centWanted);
              if (itC == itH->second.end() || itC->second.empty()) continue;

              auto g = makeGraph(itC->second.front());
              g->SetLineColor(col);
              g->SetMarkerColor(col);
              g->SetTitle(ttl.c_str());

              g->Draw(col == 1 ? "APL" : "PL SAME");
              leg.AddEntry(g.get(), det2.c_str(), "pl");

              _ownedGraphs.push_back(std::move(g));
              ++col;
            }
            leg.Draw();
            saveCanvas(c, "centrality",
                       Form("v%d_cent%s_%s.png",
                            n, centWanted.c_str(), trig.c_str()));
          }
        } // … harmonic loop
  }

  /* ---- helper: save a canvas & create directory if needed ---------- */
  void saveCanvas(TCanvas& c,
                  const std::string& subDir,
                  const std::string& fileName) const
  {
    const std::filesystem::path dir = _outDir / "FlowQA" / subDir;
    if (!std::filesystem::exists(dir))
      std::filesystem::create_directories(dir);

    const std::string full = (dir / fileName).string();
    c.SaveAs(full.c_str());
  }

  /* ------------------------------------------------------------------ */
  /*  data members                                                      */
  /* ------------------------------------------------------------------ */
  std::filesystem::path _outDir;

  /*  trigger → detector → n(=1,2,3) → centrality → list<TProfile*>  */
  std::unordered_map<
      std::string,
      std::unordered_map<
          std::string,
          std::map<
              int,
              std::map<std::string, ProfileVec>>>> _cache;

  /* keep graphs alive until end‑of‑job */
  std::vector<std::unique_ptr<TGraphErrors>> _ownedGraphs;
};


class EventQA : public QA
{
public:
    using QA::QA;

    /* ==============================================================
     *  1.  Per‑histogram processing
     * ==============================================================*/
    bool process(TObject* o) override
    {
        if (!o->InheritsFrom(TH1::Class())) return false;

        const std::string n = o->GetName();
        const bool isVz   = (n.rfind("h_vertexZ_"  ,0) == 0);
        const bool isCent = (n.rfind("h_centrality_",0) == 0);
        if (!isVz && !isCent) return false;

        /*                                                                 *
         *   directory  “…/output/<RUN>/<trigger>/EventQA/…”               *
         *   (always Inclusive, no centrality slicing here)                *
         *-----------------------------------------------------------------*/
        const fs::path outPng =
            root / "EventQA" / (isVz ? "VertexZ.png" : "Centrality.png");
        ensure_dir(outPng.parent_path());

        std::unique_ptr<TH1> h(static_cast<TH1*>(o->Clone()));
        h->SetDirectory(nullptr);  h->SetStats(0);

        const std::string runID = root.parent_path().filename().string();

        /* =============================================================
         * (A)  primary‑vertex z – robust iterative Gaussian fit
         * =========================================================== */
        if (isVz)
        {
            //----------------------------------------------------------------
            //  STEP‑0   robust seed from quantiles (immune to tails)
            //----------------------------------------------------------------
            double probs[3] = {0.16, 0.50, 0.84};
            double q[3];
            h->GetQuantiles(3, q, probs);
            double mu     = q[1];
            double sigma  = 0.5*(q[2]-q[0]);          // ~68 % central width
            if (sigma <= 0) sigma = h->GetRMS();
            if (sigma <= 0) sigma = 1;                // safety

            TF1 g("g","gaus", mu-3*sigma, mu+3*sigma);
            g.SetLineColor(kRed+1); g.SetLineWidth(2);
            g.SetParameters(h->GetMaximum(), mu, sigma);

            //----------------------------------------------------------------
            //  STEP‑1   iterative 2.5 σ shrinking until convergence
            //----------------------------------------------------------------
            constexpr int    kMaxIter   = 5;
            constexpr double kNSigmaFit = 2.5;
            constexpr double kTol       = 1e-3;       // relative change
            bool     fitOK    = false;
            TFitResultPtr res;

            for (int it=0; it<kMaxIter; ++it)
            {
                const double lo = mu - kNSigmaFit*sigma;
                const double hi = mu + kNSigmaFit*sigma;
                g.SetRange(lo,hi);
                g.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sigma);

                res = h->Fit(&g,"Q0RSLL");            // Q:quiet 0:no draw R:range S:store LL:likelihood
                fitOK = (int)res == 0;
                if (!fitOK) break;

                const double muNew    = g.GetParameter(1);
                const double sigmaNew = std::fabs(g.GetParameter(2));

                const bool conv =
                       std::fabs(muNew   - mu)    < kTol*sigma &&
                       std::fabs(sigmaNew- sigma) < kTol*sigma;

                mu = muNew;  sigma = sigmaNew;
                if (conv) break;                    // converged
            }

            /* ---- fall‑back: if everything failed keep old values ---- */
            const double muErr  = fitOK ? g.GetParError(1) : 0;
            const double sigErr = fitOK ? g.GetParError(2) : 0;

            /* ---- per‑run PNG --------------------------------------- */
            {
                TCanvas c("c_vz","", 900, 600);
                h->Draw();
                if (fitOK) g.Draw("SAME");

                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.04);
                tx.DrawLatex(0.15,0.86,Form("#mu = %.2f #pm %.2f cm", mu,  muErr));
                tx.DrawLatex(0.15,0.80,Form("#sigma = %.2f #pm %.2f cm", sigma, sigErr));
                c.SaveAs(outPng.string().c_str());
            }

            /* ---- cache for global overlays / μ,σ vs run ------------ */
            if (runID != "Combined") {
                VzPoint& p = s_points[runID];
                p.mu = mu; p.muErr = muErr; p.sigma = sigma; p.sigmaErr = sigErr;
                p.hist.reset(static_cast<TH1*>(h->Clone()));
                p.hist->SetDirectory(nullptr);
            }
        }

        /* =============================================================
         * (B)  centrality spectrum – plain plot, normalised            *
         * =========================================================== */
        else
        {
            TCanvas c; h->Draw(); c.SaveAs(outPng.string().c_str());

            if (runID != "Combined") {
                std::unique_ptr<TH1> cp(static_cast<TH1*>(h->Clone()));
                cp->SetDirectory(nullptr);
                if (cp->GetEntries() > 0) cp->Scale(1.0 / cp->GetEntries());
                s_centHists[runID] = std::move(cp);
            }
        }
        return true;
    }


    /* ==============================================================
     *  2.  Final summary (executed once, after the Combined pass)
     * ==============================================================*/
    ~EventQA() override
    {
        const std::string runID = root.parent_path().filename().string();
        if (runID != "Combined" || s_summaryWritten) return;
        s_summaryWritten = true;

        const fs::path outDir = root.parent_path();           // “…/Combined”

        /* ---------- helper: reproducible colour stream ------------ */
        auto nextColour = [](){
            static int idx = 0;
            static int palette[] = {kBlue+1,kRed+1,kGreen+2,kMagenta+2,
                                    kCyan+2,kOrange+1,kViolet,kAzure+2,
                                    kPink+1,kTeal+2};
            return palette[(idx++)% (sizeof(palette)/sizeof(int))];
        };

        /*  (A) vertex‑Z overlay  ----------------------------------- */
        if (s_points.size() > 1)
        {
            TCanvas c("c_vz_overlay","Primary‑vertex Z – all runs",900,600);
            TLegend leg(0.68,0.57,0.88,0.88); leg.SetBorderSize(0);

            bool first = true;
            for (auto& [run,p] : s_points) {
                Color_t col = nextColour();
                p.hist->SetLineColor(col); p.hist->SetLineWidth(2);
                p.hist->Draw(first ? "HIST" : "HIST SAME");
                leg.AddEntry(p.hist.get(), run.c_str(), "l");
                first = false;
            }
            leg.Draw();
            c.SaveAs((outDir/"VertexZ_AllRuns.png").string().c_str());
        }

        /*  (B) centrality overlay   -------------------------------- */
        if (s_centHists.size() > 1)
        {
            TCanvas c("c_cent_overlay","Centrality – all runs",900,600);
            TLegend leg(0.68,0.57,0.88,0.88); leg.SetBorderSize(0);

            bool first = true;
            for (auto& [run,h] : s_centHists) {
                Color_t col = nextColour();
                h->SetLineColor(col); h->SetLineWidth(2);
                h->Draw(first ? "HIST" : "HIST SAME");
                leg.AddEntry(h.get(), run.c_str(), "l");
                first = false;
            }
            leg.Draw();
            c.SaveAs((outDir/"Centrality_AllRuns.png").string().c_str());
        }

        /*  (C) μ,σ versus run number  ------------------------------ */
        if (s_points.size() > 1)
        {
            std::vector<int> runs;
            for (auto& [r,_] : s_points)
                if (std::all_of(r.begin(),r.end(),::isdigit))
                    runs.push_back(std::stoi(r));
            std::sort(runs.begin(), runs.end());

            const int n = runs.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            for (int i=0;i<n;++i) {
                const auto& p = s_points[std::to_string(runs[i])];
                x[i]=runs[i]; yMu[i]=p.mu; eMu[i]=p.muErr;
                ySi[i]=p.sigma; eSi[i]=p.sigmaErr;
            }

            auto gMu = std::make_unique<TGraphErrors>(n,x.data(),yMu.data(),nullptr,eMu.data());
            auto gSi = std::make_unique<TGraphErrors>(n,x.data(),ySi.data(),nullptr,eSi.data());
            gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
            gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

            TCanvas c("c_mu_sigma_vs_run","vertex‑Z  μ,σ  vs run",900,800);

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

            c.SaveAs((outDir/"VertexZ_MeanSigma_vs_Run.png").string().c_str());
        }
    }

    /* ==============================================================
     *  3.  Static containers (shared across all EventQA instances)
     * ==============================================================*/
    struct VzPoint {
        double mu{}, muErr{}, sigma{}, sigmaErr{};
        std::unique_ptr<TH1> hist;
    };

    static inline std::unordered_map<std::string, VzPoint>        s_points;
    static inline std::unordered_map<std::string, std::unique_ptr<TH1>> s_centHists;
    static inline bool  s_summaryWritten = false;
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

  CentList slices = discoverSlices(in.get());
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

    /* directory skeleton for this trigger -------------------------- */
    for (auto& s : slices) {
      fs::path b = fs::path(kOutputBase) / trg;
      if (s != "Inclusive") b /= ("Cent_" + s);
      for (auto sub : { "Correlations", "EMCal/pi0QA", "IHCal", "OHCal",
                        "MBD", "sEPD", "jetQA", "EventQA" })
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
}


/* ------------------------------------------------------------------ */
/*  MAIN wrapper – fan‑out over all runs with a preforked pool        */
/* ------------------------------------------------------------------ */
void analyzeRun24or25auau()
{
    using fs::path;

    /* 0. discover input ROOT files --------------------------------- */
    std::vector<path> runFiles = listRunFiles(kInputDir);
    if (runFiles.empty()) {
        log::err("No input files found in " + kInputDir.string());
        return;
    }

    /* 1. decide pool size (≤ physical cores, ≥ 1) ------------------ */
    const std::size_t nWorkers =
        std::min<std::size_t>(runFiles.size(),
                              std::max<unsigned>(1, std::thread::hardware_concurrency()));

    ROOT::TProcessExecutor exec(nWorkers);        // **preforked** process pool

    log::banner("Launching " + std::to_string(runFiles.size())   +
                " runs on "    + std::to_string(exec.GetPoolSize()) +
                " parallel processes");

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
    if (runFiles.size() > 1) {
        const path combined = kInputDir / "output_ALL_COMBINED.root";
        log::banner("Hadd – building " + combined.string());

        TFileMerger merger(/*dryRun=*/false, /*verbose=*/true);
        merger.OutputFile(combined.c_str(), "RECREATE");
        for (const auto& f : runFiles) merger.AddFile(f.c_str());

        if (!merger.Merge()) {
            log::err("TFileMerger failed – combined QA skipped");
            return;
        }
        log::ok("Combined ROOT file written");

        runOneQaPass(combined.string(),
                     (kOutputDir / "Combined").string());
    }
}
