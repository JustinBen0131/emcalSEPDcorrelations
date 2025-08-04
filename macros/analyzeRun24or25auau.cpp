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
#include <TPaletteAxis.h>
#include <filesystem>
#include <TH2.h>
#include <TH3.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <TEllipse.h>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TSequentialExecutor.hxx>
#include <ROOT/TProcessExecutor.hxx>
#include <Math/MinimizerOptions.h>
#include <TLegend.h>
#include <iostream>
#include <TGraphErrors.h>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include "sepdgeomGrid.h"

using std::string;
namespace fs = std::filesystem;

// ───────────────────────────────────────────────
//            ▼  USER SETTINGS  ▼
// ───────────────────────────────────────────────
namespace {

    /* The environment variable RUN_LOCATION selects the platform:
     *     “local”   → mac / Desktop tree (default)
     *     “sphenix” → /sphenix/u/… scratch tree
     */
    const bool kOnSphenix = []{
        const char* env = std::getenv("RUN_LOCATION");
        return env && std::string(env) == "sphenix";
    }();

    /* Base of the whole analysis directory */
    const fs::path baseDir = kOnSphenix
        ? "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations"
        : "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations";

    const fs::path kInputDir  = kOnSphenix
            ? baseDir / "output"            // → output_*.root
            : baseDir / "input" / "output"; // → Desktop/…/input/output/*.root

    const fs::path kOutputDir = kOnSphenix
        ? baseDir / "outputPlots"       /* per‑run PNGs on sPHENIX node   */
        : baseDir / "outputPlots";           /* local Desktop output           */

}   // anonymous namespace

std::string kInputFile   = "";   // gets filled inside the loop
std::string kOutputBase  = "";   // ditto
std::set<std::string> kTriggersWanted{
    "MBD_NS_geq_2",
    "MBD_NS_geq_1",
    "MBD_NS_geq_2_vtx_lt_10",
    "MBD_NS_geq_2_vtx_lt_30",
    "MBD_NS_geq_2_vtx_lt_150",
    "MBD_NS_geq_1_vtx_lt_10",
    "photon_6_plus_MBD_NS_geq_2_vtx_lt_10",
    "photon_8_plus_MBD_NS_geq_2_vtx_lt_10",
    "photon_10_plus_MBD_NS_geq_2_vtx_lt_10",
    "photon_12_plus_MBD_NS_geq_2_vtx_lt_10",
    "photon_6_plus_MBD_NS_geq_2_vtx_lt_150",
    "photon_8_plus_MBD_NS_geq_2_vtx_lt_150",
    "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
    "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
};

inline std::string prettifyTrigger(std::string s)
{
    // 1.  dedicated multi‑character tokens  (order matters!)
    s = std::regex_replace(s, std::regex(R"(_plus_)")   , " + ");
    s = std::regex_replace(s, std::regex(R"(NS)")       , "N&S");
    s = std::regex_replace(s, std::regex(R"(geq)")      , "#geq");
    s = std::regex_replace(s, std::regex(R"(vtx_lt_([0-9]+))"), "vtx < $1");

    // 2.  remaining underscores  →  space
    std::replace(s.begin(), s.end(), '_', ' ');

    // 3.  collapse any accidental double spaces
    s = std::regex_replace(s, std::regex(R"(\s{2,})"), " ");

    return s;
}

/* -----------------------------------------------------------------
 *  One single map that the rest of the macro can query at will:
 *     kPrettyTrig.at("photon_10_plus_MBD_NS_geq_2_vtx_lt_150")
 *                              ↳  "photon 10 + MBD N&S #geq 2 vtx < 150"
 * ----------------------------------------------------------------- */
static const std::unordered_map<std::string,std::string> kPrettyTrig = []{
    std::unordered_map<std::string,std::string> m;
    for (const auto& raw : kTriggersWanted)
        m.emplace(raw, prettifyTrigger(raw));
    return m;                 // NB: executed once, at start‑up
}();

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

namespace ulog {
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
  ulog::trace("save1D → " + p.string());
  TCanvas c; h->SetStats(0); h->Draw();
  ensure_dir(p.parent_path()); c.SaveAs(p.string().c_str());
}
void save2D(TH2* h,const fs::path& p,const char* opt="COLZ")
{
  ulog::trace("save2D → " + p.string());
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
    /* ──────────────── ultra‑lightweight terminal logger ──────────────── */
    enum class Lvl { DBG, INFO, WARN, ERR };
    static void log(Lvl l, const std::string& m)
    {
        static const char* tag[]{"DBG","INF","WRN","ERR"};
        std::ostream& os = (l == Lvl::ERR) ? std::cerr : std::cout;
        os << "[Pi0QA] " << tag[static_cast<int>(l)] << "  " << m << '\n';
    }

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
        FitPair()  = default;
        ~FitPair() = default;

        /* deep‑copy ctor */
        FitPair(const FitPair& other)
        {
            if (other.total) total.reset( static_cast<TF1*>(other.total->Clone()) );
            if (other.poly ) poly .reset( static_cast<TF1*>(other.poly ->Clone()) );
        }

        /* deep‑copy assignment */
        FitPair& operator=(const FitPair& other)
        {
            if (this != &other) {
                total.reset(other.total ? static_cast<TF1*>(other.total->Clone()) : nullptr);
                poly .reset(other.poly  ? static_cast<TF1*>(other.poly ->Clone()) : nullptr);
            }
            return *this;
        }

        /* move semantics – default is fine */
        FitPair(FitPair&&) noexcept            = default;
        FitPair& operator=(FitPair&&) noexcept = default;
    };

    /* ------------------------------------------------------------------ *
     *  ctor – figure out the run label from the base directory            *
     *  …/<output>/<run>/<trigger>  →  run = last element of parent_path() *
     * ------------------------------------------------------------------ */
    Pi0QA(std::string t,
          fs::path    b,
          const CentList& s,
          std::ofstream&  csvFit_) :
        QA(std::move(t), std::move(b), s),
        csvFit(csvFit_)
    {
        try {
            /* run label */
            runID = root.parent_path().filename().string();
            log(Lvl::INFO,"initialising – runID = " + runID);

            /* open (create) S/B CSV */
            fs::path p = root / "EMCal/invMassQA" / "Pi0SignalBackground.csv";
            ensure_dir(p.parent_path());
            csvSB.open(p);
            if (!csvSB.is_open())
                throw std::runtime_error("unable to open " + p.string());

            csvSB << "trigger,cent,pTlo,pThi,slice,windowSigma,sb,err\n";
            log(Lvl::DBG,"opened S/B CSV at " + p.string());
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR,std::string("ctor error: ")+ex.what());
            throw;                       // unrecoverable – propagate
        }
    }

    ~Pi0QA() override
    {
        /* any exception here must be caught – dtor must not throw */
        try {
            log(Lvl::INFO,"writing summary panels …");
            writeSummaryPanels();
            log(Lvl::INFO,"writing run summary …");
            writeRunSummary();          // only fires once in the Combined pass
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR,std::string("dtor error: ")+ex.what());
        }
        log(Lvl::INFO,"Pi0QA destructed");
    }

    /* const accessor (unchanged) */
    const auto& fitSummary() const { return _fitSummary; }
    
    // -----------------------------------------------------------------
    // Helper performing the π0 peak fit.
    //  • returns true/false (= fitOK)
    //  • fills the caller‑provided TF1 `total` and all requested values
    // -----------------------------------------------------------------
    bool doPi0Fit(TH1* h,              // [in]  histogram to be fitted
                  TF1& total,          // [out] configured gaus+pol2 function
                  double& mu,          // [out] fitted μ
                  double& muErr,       // [out] μ error
                  double& sigma,       // [out] fitted σ
                  double& sigmaErr)    // [out] σ error
    {
        if (!h) return false;

        const double piFitLo = 0.05, piFitHi = 0.35;

        /* ---- seed determination -------------------------------------- */
        const int iLo = binAt(h, piFitLo);
        const int iHi = binAt(h, piFitHi);

        int iMax = iLo; double amp0 = 0.;
        for (int i = iLo; i <= iHi; ++i)
            if (h->GetBinContent(i) > amp0) { amp0 = h->GetBinContent(i); iMax = i; }

        const double mu0    = h->GetBinCenter(iMax);
        const double sigma0 = 0.025;

        log(Lvl::DBG, Form("π0 seed  amp=%g  mu=%g", amp0, mu0));

        /* ---- configure composite model -------------------------------- */
        total.SetRange(piFitLo, piFitHi);
        total.SetParNames("A","mu","sigma","c0","c1","c2");
        total.SetParameters(amp0, mu0, 0.022, 1, 0, 0);
        total.SetParLimits(0, 0, 1e9);
        total.SetParLimits(1, 0.10, 0.17);
        total.SetParLimits(2, 0.010, 0.060);

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(3'000);

        /* ---- pre‑fit background *only* -------------------------------- */
        TF1 polyTmp("polyTmp", "pol2", piFitLo, piFitHi);
        for (int ip = binAt(h, 0.11); ip <= binAt(h, 0.16); ++ip) h->SetBinError(ip, 1e9);
        h->Fit(&polyTmp, "QRN0");
        for (int ip = binAt(h, 0.11); ip <= binAt(h, 0.16); ++ip)
            h->SetBinError(ip, std::sqrt(h->GetBinContent(ip)));

        total.SetParameters(amp0, mu0, 0.022,
                            polyTmp.GetParameter(0),
                            polyTmp.GetParameter(1),
                            polyTmp.GetParameter(2));

        /* ---- final composite fit ------------------------------------- */
        const bool ok = (h->Fit(&total, "QRN0") == 0);
        log(Lvl::INFO, std::string("π0 fit ") + (ok ? "succeeded" : "FAILED"));

        if (ok)
        {
            mu       = total.GetParameter(1);
            sigma    = total.GetParameter(2);
            muErr    = total.GetParError (1);
            sigmaErr = total.GetParError (2);
        }
        else
        {
            mu = mu0; muErr = 0.;
            sigma = sigma0; sigmaErr = 0.;
        }

        log(Lvl::DBG, Form("π0  mu=%.5f±%.5f  sigma=%.5f±%.5f",
                           mu, muErr, sigma, sigmaErr));
        return ok;
    }
    

    bool doEtaFit(TH1* h,                    // [in]  invariant‑mass histogram
                  const std::string& slice,  // [in]  current cent‑slice key
                  double& mu, double& muErr, // [out] η‑peak μ and error
                  double& sigma, double& sigmaErr)          // [out] σ and error
    {
        mu = muErr = sigma = sigmaErr = 0.;          // default outputs

        if (runID != "Combined" || !h)               // only do it once per macro
            return false;

        const double etaLo = 0.45, etaHi = 0.80, side = 0.03;

        int iLo  = binAt(h, etaLo);
        int iHi  = binAt(h, etaHi);
        int iMax = iLo;  double maxCnt = 0.;

        for (int i = iLo; i <= iHi; ++i)
            if (h->GetBinContent(i) > maxCnt) { maxCnt = h->GetBinContent(i); iMax = i; }

        if (maxCnt <= 0) return false;               // nothing to fit

        const double mu0    = h->GetBinCenter(iMax);
        const double sigma0 = 0.040;

        /* --- local background pre‑fit ---------------------------------- */
        TF1 bkg("bkg","pol2",etaLo,etaHi);
        for (int i = iLo; i <= iHi; ++i)
        {
            const double x = h->GetBinCenter(i);
            const bool inPk = (std::fabs(x - mu0) < side);
            h->SetBinError(i, inPk ? 1e9 : std::sqrt(h->GetBinContent(i)));
        }
        h->Fit(&bkg,"QN0");
        for (int i = iLo; i <= iHi; ++i)
            h->SetBinError(i, std::sqrt(h->GetBinContent(i)));

        /* --- composite η‑model ----------------------------------------- */
        TF1 gEta("gEta","gaus(0)+pol2(3)",etaLo,etaHi);
        gEta.SetParNames("A","#mu","#sigma","c0","c1","c2");
        gEta.SetParameters(maxCnt, mu0, sigma0,
                           bkg.GetParameter(0),
                           bkg.GetParameter(1),
                           bkg.GetParameter(2));
        gEta.SetParLimits(0, 0, 1e9);
        gEta.SetParLimits(1, etaLo, etaHi);
        gEta.SetParLimits(2, 0.020, 0.090);
        for (int ip = 3; ip <= 5; ++ip)
        {
            const double p  = bkg.GetParameter(ip - 3);
            const double dp = std::max(std::fabs(p) * 0.20, 1e-3);
            gEta.SetParLimits(ip, p - dp, p + dp);
        }

        const bool ok = (h->Fit(&gEta, "QRN0") == 0);
        log(Lvl::INFO, std::string("η fit ") + (ok ? "succeeded" : "FAILED"));

        if (ok)
        {
            mu       = gEta.GetParameter(1);
            muErr    = gEta.GetParError (1);
            sigma    = gEta.GetParameter(2);
            sigmaErr = gEta.GetParError (2);

            // store first successful fit per slice
            if (_storedEtaFit.count(slice) == 0)
                _storedEtaFit[slice].reset(new TF1(gEta));
        }
        else
        {
            mu = muErr = sigma = sigmaErr = 0.;
        }

        return ok;
    }

    // -----------------------------------------------------------------------
    //  MAIN ENTRY – called once per histogram
    // -----------------------------------------------------------------------
    bool process(TObject* o) override
    {
        /* --------------------------------------------------------- *
         * 0. ultra‑early sanity checks (never throw from here)       *
         * --------------------------------------------------------- */
        if (!o)
        {
            log(Lvl::ERR,"process(): received nullptr – skipped");
            return false;
        }
        if (!o->InheritsFrom(TH1::Class()))
        {
            log(Lvl::DBG,Form("process(): \"%s\" is not TH1 – skipped",o->GetName()));
            return false;
        }

        const std::string n = o->GetName();
        log(Lvl::DBG,"process(): enter – \"" + n + '"');

        /* ========================================================= *
         *  Everything below is wrapped in a single try/catch guard  *
         * ========================================================= */
        try
        {
            // ================================================================
            // (A)  UNCUT 2‑D MAPS  – names start with “Minv_vs_”
            // ================================================================
            if (o->InheritsFrom(TH2::Class()) && n.rfind("Minv_vs_",0) == 0)
            {
                log(Lvl::DBG,"→ treating object as uncut 2‑D map");

                auto* h2 = static_cast<TH2*>(o);
                if (h2->GetEntries() <= 0 || h2->Integral() <= 0)
                {
                    log(Lvl::INFO,"empty 2‑D map – skipped");
                    return false;
                }

                // ----------------- 1. shrink axes ---------------------------
                auto shrinkAxis = [&](TH2* h, bool xAxis)
                {
                    const int nBins = xAxis ? h->GetNbinsX() : h->GetNbinsY();
                    int first =  nBins+1, last = -1;

                    for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
                    for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
                    {
                        if (h->GetBinContent(ix,iy) <= 0) continue;
                        int bin = xAxis ? ix : iy;
                        first   = std::min(first, bin);
                        last    = std::max(last , bin);
                    }

                    if (last >= first)
                    {
                        TAxis* ax = xAxis ? h->GetXaxis() : h->GetYaxis();
                        const double lo = ax->GetBinLowEdge(first);
                        const double hi = ax->GetBinUpEdge (last );
                        ax->SetRangeUser(lo,hi);
                        log(Lvl::DBG,Form("   axis %c shrunk  [%d,%d] → [%.3f,%.3f]",
                                          xAxis?'X':'Y',first,last,lo,hi));
                    }
                    else
                        log(Lvl::WARN,"   shrinkAxis(): all‑zero axis – kept default range");
                };

                shrinkAxis(h2,true);   // X
                shrinkAxis(h2,false);  // Y

                // ----------------- 2. draw & save ---------------------------
                const std::string slice = sliceKey(n);
                const fs::path subDir  = "EMCal/invMassQA/cutQA";
                const fs::path outPng  = cPath(root,slice,subDir)/(n + ".png");
                ensure_dir(outPng.parent_path());

                log(Lvl::INFO,"saving 2‑D map → " + outPng.string());
                try
                {
                    TCanvas c;
                    c.SetRightMargin(0.18);
                    c.SetLogz();
                    h2->SetStats(0);
                    h2->Draw("COLZ");
                    c.SaveAs(outPng.string().c_str());
                    log(Lvl::DBG,"2‑D map saved OK");
                }
                catch (const std::exception& ex)
                {
                    log(Lvl::ERR,std::string("SaveAs failed – ")+ex.what());
                }
                return true;    // ♦ histogram consumed
            }

            // ================================================================
            // (B)  INVARIANT‑MASS SPECTRA  – names start with “mInv_”
            // ================================================================
            if (n.rfind("mInv_",0) != 0)
            {
                log(Lvl::DBG,"not an mInv_* spectrum – skipped");
                return false;
            }

            CutKey ck;
            if (!decodeInvName(n, ck))
            {
                log(Lvl::WARN,"decodeInvName() failed – skipped");
                return false;
            }

            const std::string slice = sliceKey(n);
            const bool pTInt        = (ck.pLo < 0 || ck.pHi < 0);

            // -------- bookkeeping for (E,χ²,α) combination ------------------
            const std::string combDir =
                "E"+sf3(ck.E)+"_Chi"+sf3(ck.chi)+"_Asym"+sf3(ck.asy);

            if (cutTag.empty())
            {
                cutTag = combDir;
                log(Lvl::INFO,"cutTag initialised → " + cutTag);
            }
            else if (combDir != cutTag)
            {
                log(Lvl::INFO,"cutTag switch  " + cutTag + " → " + combDir);

                /* 1.  finish the *current* cut‑combination:
                       – centrality / pT grids
                       – run‑by‑run overlay (THIS was missing)                         */
                writeSummaryPanels();
                writeRunSummary();                         // <<< added line

                /* 2.  reset per‑cut caches _before_ we start the next combination    */
                _centralHists.clear(); _storedFit.clear(); _storedEtaFit.clear();

                /* 3.  activate the new cut‑combination                              */
                cutTag = combDir;
            }

            // -------- ensure output directory exists -----------------------
            const fs::path subDir = fs::path("EMCal/invMassQA") / combDir;
            const fs::path baseDir = cPath(root,slice,subDir);
            ensure_dir(baseDir);

            const fs::path outPng = baseDir / (n + ".png");
            log(Lvl::DBG,"spectrum PNG path → " + outPng.string());

            // --------------------------- fit --------------------------------
            TH1* h = static_cast<TH1*>(o);
            const double piFitLo = 0.05, piFitHi = 0.35;
            double piMu=0.,piMuErr=0.,piSig=0.,piSigErr=0.;

            TF1 total("total","gaus(0)+pol2(3)",piFitLo,piFitHi);
            const bool fitOK = doPi0Fit(h,total,piMu,piMuErr,piSig,piSigErr);

            double etaMu=0.,etaMuErr=0.,etaSig=0.,etaSigErr=0.;
            (void)doEtaFit(h,slice,etaMu,etaMuErr,etaSig,etaSigErr);

            TF1 poly("bg","pol2",piFitLo,piFitHi);
            poly.SetParameters(total.GetParameter(3),
                               total.GetParameter(4),
                               total.GetParameter(5));
            poly.SetLineColor(kAzure+2); poly.SetLineWidth(2); poly.SetLineStyle(2);

            // -------------- S/B ratios for CSV -----------------------------
            const std::vector<double> ws = {1.25,1.5,1.75,2.0,2.25};
            for (double w: ws)
            {
                const int i1 = binAt(h,std::max(piMu-w*piSig,piFitLo));
                const int i2 = binAt(h,std::min(piMu+w*piSig,piFitHi));

                double S=0,B=0,sErr=0,bErr=0;
                for (int i=i1;i<=i2;++i)
                {
                    const double x   = h->GetBinCenter(i);
                    const double bg  = std::max(poly.Eval(x),0.);
                    const double cnt = h->GetBinContent(i);
                    B+=bg; S+=cnt-bg; sErr+=cnt; bErr+=bg;
                }

                bErr = std::sqrt(bErr); sErr = std::sqrt(sErr);
                double ratio = (B>0) ? S/B : 0, rErr = 0;
                if (ratio>0) rErr = ratio*std::sqrt((sErr*sErr)/(S*S)+(bErr*bErr)/(B*B));

                csvSB << trig << ',' << ck.E << ',' << ck.chi << ',' << ck.asy << ','
                      << slice << ',' << w << ',' << ratio << ',' << rErr << '\n';
            }

            // -------------- render single spectrum -------------------------
            log(Lvl::INFO,"drawing spectrum → " + outPng.string());
            try
            {
                TCanvas c;                     // local canvas
                h->SetStats(0);
                h->Draw();                     // histogram first

                /* draw CLONES so the pad owns its own copies, not the stack objects */
                poly.DrawCopy("SAME");
                total.DrawCopy("SAME");
                if (_storedEtaFit.count(slice))
                    _storedEtaFit[slice]->DrawCopy("SAME");   // safe for stored fits

                // labels and legend (same code, only log around) -------------
                {
                    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.038); tl.SetTextAlign(13);

                    std::string centStr;
                    if (slice=="Inclusive") centStr="Inclusive";
                    else {
                        std::smatch m;
                        if (std::regex_match(slice,m,std::regex(R"((\d{1,3})_(\d{1,3}))")))
                            centStr="Cent "+m[1].str()+"-"+m[2].str()+" %";
                        else centStr=slice;
                    }

                    std::string pTStr = pTInt
                                      ? "p_{T}-integrated"
                                      : Form("%.2f < p_{T} < %.2f GeV/#it{c}", ck.pLo, ck.pHi);

                    /* line‑1: centrality slice */
                    tl.DrawLatex(0.14, 0.93, centStr.c_str());

                    /* line‑2: pT interval (or “pT‑integrated”) */
                    tl.DrawLatex(0.14, 0.88, pTStr.c_str());

                    /* line‑3: numeric cuts (E, asym, χ²) */
                    tl.DrawLatex(0.14, 0.83,
                                 Form("E #geq %.2f GeV   Asym < %.2f   #chi^{2} < %.2f",
                                      ck.E, ck.asy, ck.chi));

                    TLegend leg(0.55,0.64,0.88,0.88);
                    leg.SetBorderSize(0); leg.SetTextAlign(12);
                    leg.AddEntry((TObject*)nullptr,
                                 Form("#pi^{0}:  #mu = %.3f #pm %.3f GeV",piMu,piMuErr),"");
                    leg.AddEntry((TObject*)nullptr,
                                 Form("          #sigma = %.3f #pm %.3f GeV",piSig,piSigErr),"");
                    if (etaMu>0)
                    {
                        leg.AddEntry((TObject*)nullptr,
                                     Form("#eta:    #mu = %.3f #pm %.3f GeV",etaMu,etaMuErr),"");
                        leg.AddEntry((TObject*)nullptr,
                                     Form("          #sigma = %.3f #pm %.3f GeV",etaSig,etaSigErr),"");
                    }
                    leg.Draw();
                }

                c.SaveAs(outPng.string().c_str());
                log(Lvl::DBG,"spectrum PNG saved OK");
            }
            catch (const std::exception& ex)
            {
                log(Lvl::ERR,std::string("spectrum SaveAs failed – ")+ex.what());
            }

            /* ------------------------------------------------------------------ *
             *  Detach every TF1 from ALL ROOT lists that may own it before it    *
             *  goes out of scope – protects against later double deletion.       *
             *  >>> VERBOSE VERSION – prints every action so we can spot issues.  *
             * ------------------------------------------------------------------ */
            auto detachTF1 = [&](TF1* f, const char* tag)
            {
                if (!f) { log(Lvl::DBG, std::string("detachTF1(")+tag+"): f==nullptr"); return; }

                log(Lvl::DBG, Form("detachTF1(%s): pointer = %p", tag, (void*)f));

                auto tryRemove = [&](TCollection* lst, const char* lstName)
                {
                    if (!lst) return false;
                    if (lst->FindObject(f)) {
                        lst->Remove(f);
                        log(Lvl::DBG, std::string("   removed from ")+lstName);
                        return true;
                    }
                    return false;
                };

                bool removed = false;
                removed |= tryRemove(h->GetListOfFunctions(),            "hist‑function list");
                removed |= tryRemove(gPad ? gPad->GetListOfPrimitives()
                                          : nullptr,                     "pad primitives");
                removed |= tryRemove(gROOT ? gROOT->GetListOfFunctions()
                                           : nullptr,                    "gROOT function list");

                if (!removed)
                    log(Lvl::WARN, std::string("detachTF1(")+tag+
                                    "): pointer not found in any list – possible logic error");
            };

            /* call for both local functions with explicit tags */
            detachTF1(&poly , "poly");
            detachTF1(&total, "total");


            // ---------- CSV & caches (original logic, unchanged) -----------
            csvFit << trig << ',' << ck.E << ',' << ck.chi << ',' << ck.asy << ','
                   << ck.pLo << ',' << ck.pHi << ','
                   << slice  << ','                           /* centrality slice */
                   << piMu   << ',' << piMuErr  << ','
                   << piSig  << ',' << piSigErr << ','
                   << etaMu  << ',' << etaMuErr << ','
                   << etaSig << ',' << etaSigErr << '\n';

            _fitSummary.emplace(n, FitInfo{slice,ck.pLo,ck.pHi,piMu,piSig,
                                           total.GetChisquare(),total.GetNDF()});

            if (fitOK && _storedFit.count(slice)==0)
            {
                _storedFit[slice].total.reset(new TF1(total));
                _storedFit[slice].poly .reset(new TF1(poly ));
            }

            if (fitOK)
            {
                auto &m = s_runPoints[slice];
                if (m.count(runID)==0)
                    m[runID] = {piMu,piMuErr,piSig,piSigErr};
            }

            cacheForOverview(slice,n,h,pTInt);

            log(Lvl::DBG,"process(): exit OK for \"" + n + '"');
            return true;
        }
        catch (const std::exception& ex)
        {
            log(Lvl::ERR,std::string("process(): C++ exception – ")+ex.what());
            return false;
        }
        catch (...)
        {
            log(Lvl::ERR,"process(): UNKNOWN exception");
            return false;
        }
    }


 private:
    // -----------------------------------------------------------------------------
    // helper: clone histogram and cache it in the correct container
    // -----------------------------------------------------------------------------
    void cacheForOverview(const std::string& slice,
                          const std::string& /*hName*/,
                          TH1* src,
                          bool pTIntegrated)
    {
        try {
            if (!src) { log(Lvl::ERR,"cacheForOverview(): src==nullptr"); return; }
            auto* cl = static_cast<TH1*>(src->Clone());
            if (!cl)  { log(Lvl::ERR,"cacheForOverview(): Clone() returned nullptr"); return; }
            cl->SetDirectory(nullptr);

            if (pTIntegrated) {                          // centrality‑summary
                log(Lvl::DBG,"cacheForOverview(): central cache  slice="+slice);
                _centralHists.emplace(slice,cl);
            } else {                                     // pT‑binned summary
                log(Lvl::DBG,"cacheForOverview(): pT‑cache slice="+slice);
                _ptHists[slice].push_back(cl);
            }
        }
        catch(const std::exception& ex){
            log(Lvl::ERR,std::string("cacheForOverview(): exception – ")+ex.what());
        }
    }
    
    // -------------------------------------------------------------------------
    // helper: draw the 2×3 centrality grid and collect μ,σ vs centrality data
    // -------------------------------------------------------------------------
    void fillCentralityGrid(TCanvas& cGrid,
                            std::vector<double>& vC,  std::vector<double>& vCerr,
                            std::vector<double>& vMu, std::vector<double>& vMuErr,
                            std::vector<double>& vSi, std::vector<double>& vSiErr,
                            std::vector<std::string>& slicesDone)
    {
        const double fitLo = 0.05, fitHi = 0.35;
        int pad = 1;                                         // pad counter

        for (const auto& sl : slices)                        /* keep slice order */
        {
            auto it = _centralHists.find(sl);
            if (it == _centralHists.end()) {
                log(Lvl::WARN,"central slice \""+sl+"\" missing in cache");
                continue;
            }

            TH1* h = it->second;
            cGrid.cd(pad++);  h->SetStats(0);

            /* ---- pick stored fit or fall‑back quick fit ------------------- */
            TF1 *fTot = nullptr, *fBg = nullptr;
            std::unique_ptr<TF1> tmpTot, tmpBg;

            if (_storedFit.count(sl)) {
                fTot = _storedFit[sl].total.get();
                fBg  = _storedFit[sl].poly .get();
            } else {
                tmpTot = std::make_unique<TF1>("fTmp","gaus(0)+pol2(3)",fitLo,fitHi);
                int iMax = h->GetMaximumBin();
                tmpTot->SetParameters(h->GetBinContent(iMax),
                                      h->GetBinCenter (iMax),0.02,1,0,0);
                h->Fit(tmpTot.get(),"QRN0");
                tmpBg  = std::make_unique<TF1>("fBgTmp","pol2",fitLo,fitHi);
                tmpBg->SetParameters(tmpTot->GetParameter(3),
                                     tmpTot->GetParameter(4),
                                     tmpTot->GetParameter(5));
                fTot = tmpTot.get();  fBg = tmpBg.get();
            }

            fTot->SetLineColor(kRed+1);   fTot->SetLineWidth(2);
            fBg ->SetLineColor(kBlue+2);  fBg ->SetLineWidth(2);
            fBg ->SetLineStyle(2);

            h->SetMaximum(1.15 * h->GetMaximum());
            h->Draw();  fBg->Draw("SAME");  fTot->Draw("SAME");
            if (_storedEtaFit.count(sl)) _storedEtaFit[sl]->Draw("SAME");

            double mu  = fTot->GetParameter(1), emu = fTot->GetParError(1);
            double sig = fTot->GetParameter(2), esig = fTot->GetParError(2);

            /* ---- ASCII‑only centrality label ----------------------------- */
            static const auto centLabel = [](const std::string& slice){
                if (slice == "Inclusive") return std::string("Inclusive");
                const auto p = slice.find('_');
                return "Centrality: " +
                       slice.substr(0,p) + "-" + slice.substr(p+1) + " %";
            };
            const std::string lbl = centLabel(sl);

            /* ---- text block ---------------------------------------------- */
            const double lm = gPad->GetLeftMargin();
            const double rm = gPad->GetRightMargin();
            const double tm = gPad->GetTopMargin();
            const double bm = gPad->GetBottomMargin();
            const bool   isTopSlice = (sl == "30_40" || sl == "40_50" || sl == "50_60");
            const bool   putBottom  = !isTopSlice;
            const double yAnchor = putBottom ? bm + 0.07 : 1.0 - tm - 0.05;
            const double dy = 0.063;
            const double xText = 1.0 - rm - 0.42;

            std::string runShort = runID;
            if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                runShort = std::to_string(std::stoi(runID));

            std::string trigLabel;
            if (auto it = kPrettyTrig.find(trig); it != kPrettyTrig.end())
                trigLabel = it->second;
            else
                trigLabel = prettifyTrigger(trig);

            TLatex tx;  tx.SetNDC();  tx.SetTextSize(0.035);  tx.SetTextAlign(13);

            if (runID != "Combined") {
                const double yRun  = putBottom ? yAnchor + 4*dy : yAnchor;
                const double yTrig = putBottom ? yAnchor + 3*dy : yAnchor - dy;
                const double yCent = putBottom ? yAnchor + 2*dy : yAnchor - 2*dy;
                const double yMu   = putBottom ? yAnchor +   dy : yAnchor - 3*dy;
                const double ySig  = putBottom ? yAnchor       : yAnchor - 4*dy;

                tx.DrawLatex(xText, yTrig, Form("Trigger: %s", trigLabel.c_str()));
                tx.DrawLatex(xText, yCent, lbl.c_str());
                tx.DrawLatex(xText, yMu , Form("#mu = %.3f #pm %.3f GeV",   mu ,emu ));
                tx.DrawLatex(xText, ySig, Form("#sigma = %.3f #pm %.3f GeV",sig,esig));
            } else {
                const double y1 = putBottom ? yAnchor + 2*dy : yAnchor;
                const double y2 = putBottom ? yAnchor +   dy : yAnchor - dy;
                const double y3 = putBottom ? yAnchor       : yAnchor - 2*dy;
                tx.DrawLatex(xText, y1, lbl.c_str());
                tx.DrawLatex(xText, y2, Form("#mu = %.3f #pm %.3f GeV",   mu ,emu ));
                tx.DrawLatex(xText, y3, Form("#sigma = %.3f #pm %.3f GeV",sig,esig));
            }

            if (sl != "Inclusive") {
                int lo = std::stoi(sl.substr(0, sl.find('_')));
                int hi = std::stoi(sl.substr(sl.find('_') + 1));
                vC .push_back(0.5 * (lo + hi));   vCerr.push_back(0.5 * (hi - lo));
                vMu.push_back(mu);                vMuErr.push_back(emu);
                vSi.push_back(sig);               vSiErr.push_back(esig);
            }

            slicesDone.push_back(sl);                              // summary recap
        } /* end loop over slices */
    }

    // -------------------------------------------------------------------------
    // helper: draw the “μ,σ versus centrality” figure and save it to disk
    // -------------------------------------------------------------------------
    void plotMuSigmaVsCentrality(const std::vector<double>& vC,
                                 const std::vector<double>& vCerr,
                                 const std::vector<double>& vMu,
                                 const std::vector<double>& vMuErr,
                                 const std::vector<double>& vSi,
                                 const std::vector<double>& vSiErr)
    {
        (void)vCerr;
        
        if (vC.empty()) return;                          // nothing to do

        const int n = vC.size();
        auto gMu = std::make_unique<TGraphErrors>(n, vC.data(), vMu.data(),
                                                  nullptr, vMuErr.data());

        auto gSi = std::make_unique<TGraphErrors>(n, vC.data(), vSi.data(),
                                                  nullptr, vSiErr.data());
        gMu->SetMarkerStyle(kFullCircle);  gMu->SetLineWidth(2);
        gSi->SetMarkerStyle(kOpenCircle);  gSi->SetLineWidth(2);

        TCanvas cGS("c_mu_sigma_vs_cent",
                    "#pi^{0} peak position / width vs centrality", 800, 800);

        /* pad geometry (unchanged) */
        const double padLeft = 0.18, padRight = 0.04, gapFrac = 0.02, fracBot = 0.30;

        TPad* p1 = new TPad("p1","",0, gapFrac + fracBot, 1, 1);
        p1->SetBottomMargin(0.04);  p1->SetTopMargin(0.04);
        p1->SetLeftMargin(padLeft); p1->SetRightMargin(padRight);
        p1->Draw();  p1->cd();
        gMu->SetTitle("; ;#mu_{#pi^{0}} (GeV/c^{2})");
        gMu->Draw("AP");
        gMu->GetXaxis()->SetLabelOffset(999);
        gMu->GetXaxis()->SetTitleOffset(999);
        gMu->GetXaxis()->SetTickLength(0);
        cGS.cd();

        TPad* p2 = new TPad("p2","",0,0,1,fracBot);
        p2->SetTopMargin(0.06);  p2->SetBottomMargin(0.38);
        p2->SetLeftMargin(padLeft);  p2->SetRightMargin(padRight);
        p2->Draw();  p2->cd();
        gSi->SetTitle(";Centrality [%];#sigma_{#pi^{0}} (GeV/c^{2})");
        gSi->Draw("AP");
        gSi->GetXaxis()->SetNdivisions(506);
        gSi->GetXaxis()->SetTitleSize(0.09);  gSi->GetXaxis()->SetLabelSize(0.07);
        gSi->GetYaxis()->SetTitleSize(0.09);  gSi->GetYaxis()->SetLabelSize(0.07);
        gSi->GetYaxis()->SetTitleOffset(0.90); gSi->GetYaxis()->SetTickLength(0.035);

        /* run / cut / trigger annotation */
        {
            /* short run ID ------------------------------------------------- */
            std::string runShort = runID;
            if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                runShort = std::to_string(std::stoi(runID));

            /* numeric cuts from the tag ----------------------------------- */
            double eCut = 0, chiCut = 0, asyCut = 0;
            std::smatch m;
            if (std::regex_match(cutTag, m,
                    std::regex(R"(E([0-9]+p[0-9]+)_Chi([0-9]+p[0-9]+)_Asym([0-9]+p[0-9]+))")))
            {
                auto p2d = [](const std::string& s)
                           { return std::stod(std::regex_replace(s, std::regex("p"), ".")); };
                eCut   = p2d(m[1]);  chiCut = p2d(m[2]);  asyCut = p2d(m[3]);
            }

            /* prettified trigger label ------------------------------------ */
            std::string trigLabel;
            if (auto it = kPrettyTrig.find(trig); it != kPrettyTrig.end())
                trigLabel = it->second;
            else
                trigLabel = prettifyTrigger(trig);

            p1->cd();                    /* ensure we query the correct pad   */
            const double yMin = gMu->GetHistogram()->GetMinimum();
            const double yMax = gMu->GetHistogram()->GetMaximum();

            double yDataTop = *std::max_element(vMu.begin(), vMu.end());
            if (!vMuErr.empty())
                yDataTop += *std::max_element(vMuErr.begin(), vMuErr.end());

            const double headroom = (yMax - yDataTop) / (yMax - yMin);
            const bool   useBottomRight = (headroom < 0.12);   /* 12 % threshold */

            const double x0 = 0.5;                // safely inside right margin
            const double yStart = 0.45;
            const double dy = 0.042;                // line spacing

            p1->cd();
            TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.03);

            tl.DrawLatex(x0, yStart,
                         Form("Run:  %s", runShort.c_str()));
            tl.DrawLatex(x0, yStart - dy,
                         Form("Cuts:  E #geq %.2f GeV, #alpha < %.2f, #chi^{2} < %.2f",
                              eCut, asyCut, chiCut));
            tl.DrawLatex(x0, yStart - 2*dy,
                         Form("Trigger: %s", trigLabel.c_str()));
        }


        fs::path pngGraph = root / "EMCal" / "invMassQA" / cutTag
                           / "Pi0Mass_Sigma_vs_Centrality.png";
        ensure_dir(pngGraph.parent_path());
        cGS.SaveAs(pngGraph.string().c_str());
        log(Lvl::INFO,"μ,σ vs centrality PNG → " + pngGraph.string());
    }

    void buildPtPanels(std::vector<std::string>& ptPlotsDone)
    {
        for (const auto& sl : slices)
        {
            if (sl == "Inclusive") continue;

            auto itH = _ptHists.find(sl);
            if (itH == _ptHists.end() || itH->second.empty()) continue;

            /* (C1)  grid with the first six pT‑binned spectra ------------------ */
            const int nShow = std::min<int>(6, itH->second.size());
            TCanvas cGridPt(Form("c_pi0_ptGrid_%s", sl.c_str()),
                            Form("#pi^{0} invariant mass – Cent %s %% (first six p_{T} bins)",
                                 sl.c_str()),
                            1800, 1000);
            cGridPt.Divide(3, 2, 0.01, 0.01);

            for (int i = 0; i < nShow; ++i)
            {
                cGridPt.cd(i + 1);
                TH1* hPt = itH->second[i];
                hPt->SetStats(0);

                const double fitLo = 0.05, fitHi = 0.35;
                std::unique_ptr<TF1> fTot, fBg;

                int iMax = hPt->GetMaximumBin();
                fTot = std::make_unique<TF1>(Form("fTot_%d_%s", i, sl.c_str()),
                                             "gaus(0)+pol2(3)", fitLo, fitHi);
                fTot->SetParameters(hPt->GetBinContent(iMax),
                                    hPt->GetBinCenter (iMax), 0.022, 1, 0, 0);
                hPt->Fit(fTot.get(), "QRN0");

                fBg  = std::make_unique<TF1>(Form("fBg_%d_%s",  i, sl.c_str()),
                                             "pol2", fitLo, fitHi);
                fBg->SetParameters(fTot->GetParameter(3),
                                   fTot->GetParameter(4),
                                   fTot->GetParameter(5));

                hPt->Draw();
                fBg->SetLineColor(kBlue  + 2);  fBg->SetLineStyle(2);  fBg->SetLineWidth(2);
                fBg->DrawCopy("SAME");
                fTot->SetLineColor(kRed  + 1);  fTot->SetLineWidth(2);
                fTot->DrawCopy("SAME");
                if (_storedEtaFit.count(sl)) _storedEtaFit[sl]->Draw("SAME");

                CutKey tmpCK;
                if (decodeInvName(hPt->GetName(), tmpCK))
                {
                    std::string centLbl;
                    if (sl == "Inclusive")
                        centLbl = "Inclusive";
                    else
                        centLbl = Form("Cent %s %%", sl.c_str());

                    std::string ptLbl;
                    if (tmpCK.pLo >= 0 && tmpCK.pHi >= 0)
                        ptLbl = Form("%.2f < p_{T} < %.2f GeV/#it{c}",
                                     tmpCK.pLo, tmpCK.pHi);
                    else
                        ptLbl = "p_{T} IND";

                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045); tx.SetTextAlign(13);
                    tx.DrawLatex(0.14, 0.93, centLbl.c_str());
                    tx.DrawLatex(0.14, 0.88, ptLbl.c_str());
                }
            }

            const std::string slDir = (sl == "Inclusive" || sl == "noCentralityDep" ||
                                       sl.rfind("Cent_",0) == 0) ? sl : "Cent_" + sl;
            fs::path pngGridPt = root / "EMCal" / "invMassQA" / cutTag / slDir /
                                 "Pi0Mass_First6pTbins.png";
            ensure_dir(pngGridPt.parent_path());
            cGridPt.SaveAs(pngGridPt.string().c_str());
            log(Lvl::DBG, "pT‑grid PNG  → " + pngGridPt.string());
            ptPlotsDone.push_back(sl + " (grid)");

            /* (C2)  μ,σ versus pT ------------------------------------------- */
            std::map<double, FitInfo> byPt;
            for (const auto& [hName, fi] : _fitSummary)
                if (fi.slice == sl && fi.pLo >= 0 && fi.pHi >= 0)
                    byPt[0.5 * (fi.pLo + fi.pHi)] = fi;

            if (byPt.size() < 2) continue;

            int n = byPt.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            int k = 0;
            for (const auto& [pt, fi] : byPt)
            {
                x[k]  = pt;        yMu[k] = fi.mean;   eMu[k] = 0;
                ySi[k] = fi.sigma; eSi[k] = 0;         ++k;
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
            TPad* p1 = new TPad("p1","",0, gap + fracBot, 1, 1);
            p1->SetBottomMargin(0.04); p1->SetTopMargin(0.04);
            p1->SetLeftMargin(padLeft); p1->SetRightMargin(padRight);
            p1->Draw(); p1->cd();
            gMu->SetTitle("; ;m_{#pi^{0}}  (GeV)"); gMu->Draw("AP");
            gMu->GetXaxis()->SetLabelOffset(999);
            gMu->GetXaxis()->SetTitleOffset(999);

            cPT.cd();
            TPad* p2 = new TPad("p2","",0,0,1,fracBot);
            p2->SetTopMargin(0.06); p2->SetBottomMargin(0.38);
            p2->SetLeftMargin(padLeft); p2->SetRightMargin(padRight);
            p2->Draw(); p2->cd();
            gSi->SetTitle(";p_{T}  [GeV/#it{c}];#sigma_{#pi^{0}}  (GeV)");
            gSi->Draw("AP");
            gSi->GetXaxis()->SetNdivisions(506);
            gSi->GetXaxis()->SetTitleSize(0.09);  gSi->GetXaxis()->SetLabelSize(0.07);
            gSi->GetYaxis()->SetTitleSize(0.09);  gSi->GetYaxis()->SetLabelSize(0.07);
            gSi->GetYaxis()->SetTitleOffset(0.90); gSi->GetYaxis()->SetTickLength(0.035);

            fs::path pngPT = root / "EMCal" / "invMassQA" / cutTag / slDir /
                             "Pi0Mass_Sigma_vs_pT.png";
            ensure_dir(pngPT.parent_path());
            cPT.SaveAs(pngPT.string().c_str());
            log(Lvl::DBG, "μ,σ vs pT PNG → " + pngPT.string());
            ptPlotsDone.push_back(sl + " (μσ‑vs‑pT)");
        }
    }

    void writeSummaryPanels()
    {
        if (_centralHists.empty()) {
            log(Lvl::INFO,"writeSummaryPanels(): _centralHists is empty – nothing to do");
            return;
        }

        log(Lvl::INFO,"writeSummaryPanels(): building centrality overview");

        std::vector<std::string> slicesDone;          // ← for final terminal table
        std::vector<std::string> ptPlotsDone;

        try {
            TCanvas cGrid("c_pi0Cent","#pi0 – all centralities",1800,1000);
            gStyle->SetOptTitle(0);
            cGrid.SetTopMargin(0.12);
            cGrid.Divide(3,2,0.01,0.01);

            std::vector<double> vC,vCerr,vMu,vMuErr,vSi,vSiErr;
            const double fitLo=0.05, fitHi=0.35;
            int pad=1;

            fillCentralityGrid(cGrid,
                               vC, vCerr,
                               vMu, vMuErr,
                               vSi, vSiErr,
                               slicesDone);      // one call replaces the whole loop
            /* ---- add run + cut header, then save the 2×3 grid --------------- */
            {
                /* build stripped run label (no leading zeros) */
                std::string runShort = runID;
                if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                    runShort = std::to_string(std::stoi(runID));

                /* decode the E / χ² / α cuts from cutTag */
                double eCut = 0., chiCut = 0., asyCut = 0.;
                std::smatch m;
                if (std::regex_match(cutTag, m,
                        std::regex(R"(E([0-9]+p[0-9]+)_Chi([0-9]+p[0-9]+)_Asym([0-9]+p[0-9]+))")))
                {
                    auto p2d = [](const std::string& s)
                               { return std::stod(std::regex_replace(s, std::regex("p"), ".")); };
                    eCut   = p2d(m[1]);
                    chiCut = p2d(m[2]);
                    asyCut = p2d(m[3]);
                }

                /* draw the header once – use NDC so it sits above every pad */
                cGrid.cd();
                TLatex tl;  tl.SetNDC();  tl.SetTextSize(0.028);  tl.SetTextAlign(22);
                tl.DrawLatex(0.50, 0.96,
                    Form("Run %s   |   E #geq %.2f GeV, #alpha < %.2f, #chi^{2} < %.2f",
                         runShort.c_str(), eCut, asyCut, chiCut));
            }

            fs::path pngGrid = root/"EMCal"/"invMassQA"/cutTag
                               /"Pi0Mass_AllCentrality.png";
            ensure_dir(pngGrid.parent_path());
            cGrid.SaveAs(pngGrid.string().c_str());

            plotMuSigmaVsCentrality(vC, vCerr,
                                    vMu, vMuErr,
                                    vSi, vSiErr);   // one call replaces the whole block


            /* ----------- pT‑bin overview & μ,σ vs pT -------------------------- */
            buildPtPanels(ptPlotsDone);

            /* -------------- terminal summary table ------------------------ */
            std::ostringstream oss;
            oss << "writeSummaryPanels(): summary\n";
            oss << "  centrality slices rendered ("<<slicesDone.size()<<") : ";
            for (const auto& s: slicesDone) oss << s << ' ';
            oss << "\n  pT‑dependent canvases   ("<<ptPlotsDone.size()<<") : ";
            for (const auto& s: ptPlotsDone) oss << s << ' ';
            log(Lvl::INFO,oss.str());

        }
        catch(const std::exception& ex){
            log(Lvl::ERR,std::string("writeSummaryPanels(): exception – ")+ex.what());
        }
    }

    //--------------------------------------------------------------------
    //  Write run‑by‑run π0‑mass summary   (called once, after “Combined”)
    //--------------------------------------------------------------------
    void writeRunSummary()
    {
        /* ────────────────────────────────────────────────────────────────
         * 0.  ENTRY BANNER
         * ──────────────────────────────────────────────────────────────── */
        log(Lvl::INFO,
            "══════════════════════════════════════════════════════════════");
        log(Lvl::INFO,
            "writeRunSummary()  ➜  entering run‑summary builder");

        /* ----------------------------------------------------------------
         * 1.  Guard clauses
         * ---------------------------------------------------------------- */
        if (runID != "Combined") {
            log(Lvl::INFO,
                "writeRunSummary(): current pass = \"" + runID +
                "\"  →  run‑summary generation skipped (only generated in "
                "\"Combined\" pass)");
            return;
        }

        /* ----------------------------------------------------------------
         * 2.  Ensure we have data  (re‑hydrate from CSV files if needed)
         * ---------------------------------------------------------------- */
        if (s_runPoints.empty()) {
            log(Lvl::INFO,
                "writeRunSummary(): in‑memory cache empty – attempting to "
                "rebuild from per‑run CSV files");

            /* helper that parses one CSV and populates s_runPoints ---------- */
            auto loadCsv = [&](const fs::path& csvPath)
            {
                log(Lvl::DBG,"   · scanning " + csvPath.string());
                std::ifstream csv(csvPath);
                if (!csv) {
                    log(Lvl::DBG,"     ↳ file missing or unreadable – ignored");
                    return;
                }

                std::string line;
                std::getline(csv, line);               // discard header
                unsigned int nRows = 0;

                while (std::getline(csv, line)) {
                    if (line.empty()) continue;
                    std::stringstream ss(line);
                    std::string fld[15];               // first 15 columns are enough
                    for (int i = 0; i < 15 &&
                           std::getline(ss, fld[i], ','); ++i) {}

                    /* CSV columns (after earlier edit):
                       0 trig | 1 E | 2 Chi | 3 Asym | 4 pTlo | 5 pThi | 6 cent
                       7 μ | 8 μErr | 9 σ | 10 σErr | …                               */

                    if (fld[0] != trig)          continue;      // wrong trigger
                    if (std::stod(fld[4]) >= 0) continue;      // pT‑integrated only

                    const std::string& cent = fld[6];
                    const std::string runIDcsv =
                            csvPath.parent_path().filename().string();

                    s_runPoints[cent][runIDcsv] = {
                        std::stod(fld[7]), std::stod(fld[8]),
                        std::stod(fld[9]), std::stod(fld[10])
                    };
                    ++nRows;
                }
                log(Lvl::DBG, Form("     ↳ %u row(s) imported", nRows));
            };

            const fs::path runsRoot = root.parent_path().parent_path();
            log(Lvl::INFO,
                "   • searching run directories under " + runsRoot.string());
            unsigned int csvFound = 0;
            for (const auto& de : fs::directory_iterator(runsRoot)) {
                if (!de.is_directory())                           continue;
                if (de.path().filename() == "Combined")           continue;

                const fs::path csvPath = de.path() /
                                          "InvariantMassSummary.csv";
                if (fs::exists(csvPath)) ++csvFound;
                loadCsv(csvPath);
            }
            log(Lvl::INFO,
                Form("   • CSV scan finished  |  files found = %u", csvFound));

            if (s_runPoints.empty()) {
                log(Lvl::WARN,
                    "writeRunSummary(): still no data after CSV scan – "
                    "aborting run‑summary generation");
                return;
            }
        }
        /* ----------------------------------------------------------------
         * 3.  Duplicate‑prevention (per TRIGGER × CUT combination)
         * ---------------------------------------------------------------- */
        static std::unordered_set<std::string> s_done;
        const std::string tagKey = cutTag + "_" + trig;   // e.g.  E2p00_Chi1p00_…
        if (s_done.count(tagKey)) {
            log(Lvl::INFO,
                "writeRunSummary(): already executed for key \"" + tagKey +
                "\" – skipping duplicate");
            return;
        }
        s_done.insert(tagKey);

        log(Lvl::INFO,
            "▶ building run‑summary  |  cut = \"" + cutTag +
            "\"  |  slices = " + std::to_string(s_runPoints.size()));

        /* ----------------------------------------------------------------
         * 4.  Fixed colour palette (unchanged)
         * ---------------------------------------------------------------- */
        const int cols[] = {kBlue+1,kRed+1,kGreen+2,kMagenta+2,
                            kOrange+1,kCyan+2,kSpring+5,kPink+1};
        constexpr int nCols = sizeof(cols)/sizeof(int);

        try
        {
            std::vector<TGraphErrors*> gMuList, gSiList;
            std::vector<std::string>   slicesDone;
            /* legend – small font, transparent, top-right corner */
            TLegend leg(0.70, 0.7, 0.89, 0.9);        // x1,y1,x2,y2 in NDC
            leg.SetBorderSize(0);                       // no frame
            leg.SetFillStyle(0);                        // fully transparent
            leg.SetTextFont(42);
            leg.SetTextSize(0.03);                     // ≈ 70 % of previous size

            int colourIdx = 0;                          // running colour index used below

            std::vector<std::string> masterRuns;
            for (const auto& [sl, mp] : s_runPoints)
                for (const auto& [r, _] : mp) {
                    if (!std::all_of(r.begin(), r.end(), ::isdigit))     // skip non‑numeric IDs
                        continue;
                    if (std::find(masterRuns.begin(), masterRuns.end(), r) == masterRuns.end())
                        masterRuns.push_back(r);
                }

            std::sort(masterRuns.begin(), masterRuns.end(),
                      [](const std::string& a, const std::string& b)
                      { return std::stoi(a) < std::stoi(b); });


            std::unordered_map<std::string,int> runIdx;
            for (std::size_t i = 0; i < masterRuns.size(); ++i)
                runIdx[masterRuns[i]] = static_cast<int>(i);

            const auto& runLabels = masterRuns;        // used later for the custom axis

            log(Lvl::INFO,"   • regenerating per‑slice graphs");
            
            std::vector<std::string> sliceOrder;             // Inclusive first
            if (s_runPoints.count("Inclusive"))  sliceOrder.push_back("Inclusive");
            for (const auto& [sl,_] : s_runPoints)
                if (sl != "Inclusive") sliceOrder.push_back(sl);

            /* numeric sort for “Cent_x_y” etc. so 0‑10 precedes 10‑20 … */
            std::sort(sliceOrder.begin()+ (sliceOrder.front()=="Inclusive"),
                      sliceOrder.end(),
                      [](const std::string& a, const std::string& b)
                      {
                          auto num = [](const std::string& s)
                          {
                              std::size_t p = s.find_first_of("0123456789");
                              return (p==std::string::npos) ? 0 : std::stoi(s.substr(p));
                          };
                          return num(a) < num(b);
                      });

            const int   nSl    = static_cast<int>(sliceOrder.size());
            const double dx    = 0.8 / nSl;                 // bin width ~0.8 units
            std::unordered_map<std::string,double> sliceOffset;
            for (int i = 0; i < nSl; ++i)
                sliceOffset[sliceOrder[i]] = (i - (nSl-1)/2.0) * dx;


            for (const auto& [slice, mp] : s_runPoints)
            {
                /* ---------- prepare sorted run list ----------------------- */
                std::vector<std::pair<int,std::string>> runList;
                for (const auto& [runStr,_] : mp)
                    if (std::all_of(runStr.begin(), runStr.end(), ::isdigit))
                        runList.emplace_back(std::stoi(runStr), runStr);

                if (runList.empty()) {
                    log(Lvl::WARN,"     ↳ slice \"" + slice +
                                   "\" skipped – contains no numeric runs");
                    continue;
                }
                std::sort(runList.begin(), runList.end(),
                          [](auto& a, auto& b){ return a.first < b.first; });

                std::vector<double> xOverlay, xCentre, yMu, eMu, ySi, eSi;
                for (const auto& [runNum, runStr] : runList)
                {
                    const int idx         = runIdx[runStr];          // 0,1,2,…
                    const RunPoint& p     = mp.at(runStr);

                    const double off      = sliceOffset.at(slice);   // symmetric spread inside bin
                    xOverlay.push_back(idx + off);                   // overlay coordinate
                    xCentre .push_back(static_cast<double>(idx));    // exact bin centre

                    yMu.push_back(p.mu);      eMu.push_back(p.muErr);
                    ySi.push_back(p.sigma);   eSi.push_back(p.sigmaErr);
                }
                const int n = static_cast<int>(xOverlay.size());

                /* overlay graphs – smaller markers so they do not overlap */
                auto gMu = new TGraphErrors(n,xOverlay.data(),yMu.data(),nullptr,eMu.data());
                auto gSi = new TGraphErrors(n,xOverlay.data(),ySi.data(),nullptr,eSi.data());
                gMu->SetMarkerSize(0.70);
                gSi->SetMarkerSize(0.70);

                /* per-slice graphs – points sit centred in the run bin */
                auto gMuCent = new TGraphErrors(n,xCentre.data(),yMu.data(),nullptr,eMu.data());
                auto gSiCent = new TGraphErrors(n,xCentre.data(),ySi.data(),nullptr,eSi.data());

                const int col = cols[colourIdx++ % nCols];
                gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
                gSi->SetMarkerStyle(kFullCircle); gSi->SetLineWidth(2);
                gMu->SetMarkerColor(col);         gMu->SetLineColor(col);
                gSi->SetMarkerColor(col);         gSi->SetLineColor(col);

                gMuList.push_back(gMu);
                gSiList.push_back(gSi);

                const std::string lbl = (slice=="Inclusive")
                                            ? "Inclusive"
                                            : "Cent " + slice + " %";
                leg.AddEntry(gMu,lbl.c_str(),"pl");
                slicesDone.push_back(slice);

                /* ---------- verbose recap per slice ----------------------- */
                log(Lvl::DBG,
                    Form("     ↳ %-10s  runs = %3d  colourIdx = %2d",
                         slice.c_str(), n, colourIdx-1));

                /* ---------- per‑slice PNG -------------------------------- */
                {
                    /* ─────────────────────────────────────────────────────────────────────
                     * Per‑slice canvas – identical layout to overlay (1600×900, custom
                     * ticks, run‑number labels, etc.) but containing only this slice.
                     * ─────────────────────────────────────────────────────────────────── */
                    TCanvas cS(Form("c_mu_sigma_vs_run_%s",slice.c_str()),"",1600,900);

                    /* ---------- μ‑pad -------------------------------------------------- */
                    TPad *pTop = new TPad("pTop","",0,0.35,1,1);
                    pTop->SetBottomMargin(0.02); pTop->Draw(); pTop->cd();

                    gMuCent->SetTitle(";Run number;m_{#pi^{0}}  (GeV)");
                    gMuCent->SetMarkerStyle(gMu->GetMarkerStyle());
                    gMuCent->SetMarkerColor(gMu->GetMarkerColor());
                    gMuCent->SetLineColor  (gMu->GetLineColor());
                    gMuCent->Draw("AP");

                    /* custom x‑axis (baseline + ticks) – same code as overlay ------------ */
                    {
                        const int nPtsLoc = runLabels.size();
                        gMuCent->GetHistogram()->GetXaxis()->SetLimits(-0.5, nPtsLoc - 0.5);   // <<< align frame
                        TH1 *fr = gMuCent->GetHistogram();
                        fr->GetXaxis()->SetLabelOffset(999);
                        fr->GetXaxis()->SetTickLength(0);

                        gPad->Update();
                        const double y0    = gPad->GetUymin();
                        const double yTick = y0 + (gPad->GetUymax()-gPad->GetUymin())*0.03;

                        TLine *base = new TLine(-0.5, y0, nPtsLoc-0.5, y0);
                        base->SetLineWidth(1); base->Draw();

                        for (int i = 0; i <= nPtsLoc; ++i) {
                            const double xEdge = -0.5 + i;
                            TLine *tick = new TLine(xEdge, y0, xEdge, yTick);
                            tick->SetLineWidth(1); tick->Draw();
                        }
                        gPad->Modified(); gPad->Update();
                    }

                    /* ---------- σ‑pad -------------------------------------------------- */
                    cS.cd();
                    TPad *pBot = new TPad("pBot","",0,0,1,0.32);
                    pBot->SetTopMargin(0.02); pBot->SetBottomMargin(0.30);
                    pBot->Draw(); pBot->cd();

                    gSiCent->SetTitle("; ;#sigma_{#pi^{0}}  (GeV)");
                    gSiCent->SetMarkerStyle(gSi->GetMarkerStyle());
                    gSiCent->SetMarkerColor(gSi->GetMarkerColor());
                    gSiCent->SetLineColor  (gSi->GetLineColor());
                    gSiCent->Draw("AP");
                    
                    gSiCent->GetYaxis()->SetTitleSize(0.060);   // default ≈0.045
                    gSiCent->GetYaxis()->SetTitleOffset(0.90);  // move a bit closer

                    /* identical custom axis + run‑number labels ------------------------- */
                    {
                        const int nPtsLoc = runLabels.size();
                        gSiCent->GetHistogram()->GetXaxis()->SetLimits(-0.5, nPtsLoc - 0.5);   // <<< align frame
                        TH1 *fr = gSiCent->GetHistogram();
                        fr->GetXaxis()->SetLabelOffset(999);
                        fr->GetXaxis()->SetTickLength(0);

                        gPad->Update();
                        const double y0    = gPad->GetUymin();
                        const double yTick = y0 + (gPad->GetUymax()-gPad->GetUymin())*0.03;

                        TLine *base = new TLine(-0.5, y0, nPtsLoc-0.5, y0);
                        base->SetLineWidth(1); base->Draw();

                        for (int i = 0; i <= nPtsLoc; ++i) {
                            const double xEdge = -0.5 + i;
                            TLine *tick = new TLine(xEdge, y0, xEdge, yTick);
                            tick->SetLineWidth(1); tick->Draw();
                        }

                        /* centred run‑number labels (identical to overlay) */
                        const double x0   = gPad->GetLeftMargin();
                        const double xW   = 1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin();
                        const double yNDC = gPad->GetBottomMargin() * 0.75;

                        TLatex tl; tl.SetNDC(); tl.SetTextFont(42);
                        tl.SetTextAngle(55); tl.SetTextSize(0.042); tl.SetTextAlign(23);

                        for (int i = 0; i < nPtsLoc; ++i) {
                            const double xFrac = (i + 0.5) / static_cast<double>(nPtsLoc);
                            const double xNDC  = x0 + xFrac * xW;
                            const std::string lbl = std::to_string(std::stoi(runLabels[i]));
                            tl.DrawLatex(xNDC, yNDC, lbl.c_str());
                        }

                        gPad->Modified(); gPad->Update();
                    }

                    std::string subDir  = "";                                    // never create a sub‑folder
                    std::string fileTag = (slice=="Inclusive" || slice=="noCentralityDep")
                                            ? slice                              // “Inclusive”, “noCentralityDep”
                                            : "Cent_" + slice;                   // “Cent_0_10”, “Cent_10_20”, …


                    fs::path pngSlice = root / "EMCal" / "invMassQA" / cutTag;
                    if (!subDir.empty()) pngSlice /= subDir;        // add sub‑folder only if needed
                    ensure_dir(pngSlice);                           // make sure it exists

                    /* identical headline to overlay ------------------------------------ */
                    {
                        const std::string runLo = masterRuns.empty()
                                                    ? "?" : std::to_string(std::stoi(masterRuns.front()));
                        const std::string runHi = masterRuns.empty()
                                                    ? "?" : std::to_string(std::stoi(masterRuns.back()));
                        const int nRun = static_cast<int>(masterRuns.size());

                        const std::string sliceDesc =
                                (slice == "Inclusive" || slice == "noCentralityDep")
                                    ? slice
                                    : Form("Cent %s %%", slice.c_str());

                        const std::string canvTitle =
                            Form("Invariant Mass Summary %s #rightarrow %s (%d runs), %s",
                                 runLo.c_str(), runHi.c_str(), nRun, sliceDesc.c_str());

                        cS.cd();
                        TLatex head;
                        head.SetNDC();
                        head.SetTextFont(42);
                        head.SetTextAlign(22);
                        head.SetTextSize(0.050);
                        head.DrawLatex(0.50, 0.97, canvTitle.c_str());
                        cS.Modified();
                        cS.Update();
                    }

                    pngSlice /= "Pi0Mass_Sigma_vs_Run_" + fileTag + ".png";
                    cS.SaveAs(pngSlice.string().c_str());

                    log(Lvl::INFO,"       • per-slice PNG → " + pngSlice.string());
                }
            } // end slice loop

            if (gMuList.empty()) {
                log(Lvl::ERR,
                    "writeRunSummary(): no valid slices – cannot build overlay");
                return;
            }

            /* ----------------------------------------------------------------
             * 6.  Overlay canvas (μ top, σ bottom)   +   global PNG
             * ---------------------------------------------------------------- */
            log(Lvl::INFO,"   • drawing overlay canvas with all slices");

            /* build a descriptive canvas title – strip any leading zeros */
            const std::string runLo = masterRuns.empty()
                                        ? "?"
                                        : std::to_string(std::stoi(masterRuns.front()));
            const std::string runHi = masterRuns.empty()
                                        ? "?"
                                        : std::to_string(std::stoi(masterRuns.back()));
            const int nRun = static_cast<int>(masterRuns.size());

            const std::string canvTitle =
                Form("Invariant Mass Summary %s #rightarrow %s (%d runs), Overlay",
                     runLo.c_str(), runHi.c_str(), nRun);

            TCanvas cR("c_mu_sigma_vs_run_allCent", "", 1600, 900);
            
            /* -------------------------------------------------------------
             * draw the headline last so it is not hidden by the sub‑pads
             * ----------------------------------------------------------- */
            cR.cd();                 // back to the canvas (not to p1 / p2)

            TPad* p1 = new TPad("p1","",0,0.35,1,1);
            p1->SetBottomMargin(0.02); p1->Draw(); p1->cd();
            gMuList.front()->SetTitle(";Run number;m_{#pi^{0}}  (GeV)");

            const int nPts = runLabels.size();          // how many runs will be shown
            /* draw every slice, but ensure “Inclusive” is over‑plotted last so it remains
               visible (blue markers will sit on top of the others)                         */
            int idxInc = -1;
            for (std::size_t i = 0; i < slicesDone.size(); ++i)
                if (slicesDone[i] == "Inclusive") { idxInc = static_cast<int>(i); break; }

            if (idxInc >= 0) {
                gMuList[idxInc]->Draw("AP");                 // frame comes from Inclusive
                for (std::size_t i = 0; i < gMuList.size(); ++i)
                    if (static_cast<int>(i) != idxInc)
                        gMuList[i]->Draw("P SAME");          // other slices
                gMuList[idxInc]->Draw("P SAME");             // Inclusive on top
            } else {                                         // fallback
                for (std::size_t i = 0; i < gMuList.size(); ++i)
                    gMuList[i]->Draw(i == 0 ? "AP" : "P SAME");
            }

            /* expand the x–range so that the first marker sits at the first tick,
               not on the frame’s left border                                                */
            if (nPts > 0)
            {
                const double xMin = -0.5;               // centre of the first tick
                const double xMax =  nPts - 0.5;        // centre of the last tick

                gMuList.front()->GetHistogram()->GetXaxis()->SetLimits(xMin, xMax);
                gSiList.front()->GetHistogram()->GetXaxis()->SetLimits(xMin, xMax);
                gPad->Modified();                       // force the pad to redraw axes
            }

            leg.Draw();
            
            {
                TH1* fr = gMuList.front()->GetHistogram();
                if (fr) {
                    fr->GetXaxis()->SetLabelOffset(999);   // hide labels
                    fr->GetXaxis()->SetTickLength(0);      // suppress default ticks
                }

                gPad->Update();                            // now pad extents are known
                const double y0 = gPad->GetUymin();        // bottom edge in user coords

                /* -----------------------------------------------------------------
                 * Hand–built x‑axis: thin baseline plus primary ticks at every bin
                 * edge (−0.5, 0.5, 1.5 …).  Nothing is left to ROOT’s discretion.
                 * ----------------------------------------------------------------*/
                const double xMinUsr = -0.5;
                const double xMaxUsr =  nPts - 0.5;

                const double yBase   = y0;                                                // axis baseline (bottom of frame)
                const double yTickHi = y0 + (gPad->GetUymax() - gPad->GetUymin()) * 0.03; // tick length ≈ 3 % pad height

                /* baseline --------------------------------------------------------- */
                TLine *base = new TLine(xMinUsr, yBase, xMaxUsr, yBase);
                base->SetLineWidth(1);
                base->Draw();

                /* primary ticks ---------------------------------------------------- */
                for (int i = 0; i <= nPts; ++i) {      // i = 0 … nPts   →  (−0.5 … nPts−0.5)
                    const double xEdge = xMinUsr + i;
                    TLine *tick = new TLine(xEdge, yBase, xEdge, yTickHi);
                    tick->SetLineWidth(1);
                    tick->Draw();
                }

                gPad->Modified();
                gPad->Update();
            }
            
            cR.cd();
            TPad* p2 = new TPad("p2","",0,0,1,0.32);
            p2->SetTopMargin(0.02); p2->SetBottomMargin(0.30);
            p2->Draw(); p2->cd();
            gSiList.front()->SetTitle("; ;#sigma_{#pi^{0}}  (GeV)");
            int idxIncSi = -1;
            for (std::size_t i = 0; i < slicesDone.size(); ++i)
                if (slicesDone[i] == "Inclusive") { idxIncSi = static_cast<int>(i); break; }

            if (idxIncSi >= 0) {
                gSiList[idxIncSi]->Draw("AP");
                gSiList.front()->GetYaxis()->SetTitleSize(0.060);
                gSiList.front()->GetYaxis()->SetTitleOffset(0.55);
                
                for (std::size_t i = 0; i < gSiList.size(); ++i)
                    if (static_cast<int>(i) != idxIncSi)
                        gSiList[i]->Draw("P SAME");            // draw other slices
                gSiList[idxIncSi]->Draw("P SAME");             // Inclusive on top
            } else {                                           // fallback
                for (std::size_t i = 0; i < gSiList.size(); ++i)
                    gSiList[i]->Draw(i == 0 ? "AP" : "P SAME");
            }


            {
                const int nPts = runLabels.size();
                if (nPts == 0) {
                    log(Lvl::WARN,"writeRunSummary(): no run labels");
                } else {

                    /* -------------------------------------------------------------
                     * (1)  suppress the default numeric x–axis on the σ‑pad frame
                     * ----------------------------------------------------------- */
                    auto *fr = gSiList.front()->GetHistogram();
                    fr->GetXaxis()->SetLabelOffset(999);
                    fr->GetXaxis()->SetTickLength(0);

                    /* -------------------------------------------------------------
                     * (2)  tick marks only – no labels – with an invisible TGaxis
                     * ----------------------------------------------------------- */
                    gPad->Update();                             // make sure geometry is fixed
                    const double y0 = gPad->GetUymin();         // bottom edge of σ‑pad

                    /* -----------------------------------------------------------------
                     * Hand-built x-axis for σ-pad: baseline plus primary ticks at every
                     * bin edge (−0.5, 0.5, 1.5 …).  No automatic TGaxis involved.
                     * ----------------------------------------------------------------*/
                    const double xMinUsr = -0.5;
                    const double xMaxUsr =  nPts - 0.5;

                    const double yBase   = y0;                                                // axis baseline
                    const double yTickHi = y0 + (gPad->GetUymax() - gPad->GetUymin()) * 0.03; // tick length ≈3 %

                    /* baseline --------------------------------------------------------- */
                    TLine *baseσ = new TLine(xMinUsr, yBase, xMaxUsr, yBase);
                    baseσ->SetLineWidth(1);
                    baseσ->Draw();

                    /* primary ticks ---------------------------------------------------- */
                    for (int i = 0; i <= nPts; ++i) {
                        const double xEdge = xMinUsr + i;                                     // −0.5, 0.5, …
                        TLine *tickσ = new TLine(xEdge, yBase, xEdge, yTickHi);
                        tickσ->SetLineWidth(1);
                        tickσ->Draw();
                    }

                    gPad->Modified();
                    gPad->Update();

                    const double x0   = gPad->GetLeftMargin();
                    const double xW   = 1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin();
                    const double yNDC = gPad->GetBottomMargin() * 0.5;

                    TLatex tl;
                    tl.SetNDC();
                    tl.SetTextFont(42);
                    tl.SetTextAngle(55);
                    tl.SetTextSize(0.045);
                    tl.SetTextAlign(23);

                    for (int i = 0; i < nPts; ++i)
                    {
                        const double xFrac = (i + 0.5) / static_cast<double>(nPts);
                        const double xNDC  = x0 + xFrac * xW;
                        const std::string lbl = std::to_string(std::stoi(runLabels[i]));

                        tl.DrawLatex(xNDC, yNDC, lbl.c_str());
                    }
                    gPad->Modified();
                    gPad->Update();
                }
            }
            
            /* -------------------------------------------------------------
             * headline: draw it last so it sits on top of both pads
             * ----------------------------------------------------------- */
            cR.cd();
            TLatex head;
            head.SetNDC();
            head.SetTextFont(42);
            head.SetTextAlign(22);      // centred horizontally & vertically
            head.SetTextSize(0.050);
            head.DrawLatex(0.50, 0.97, canvTitle.c_str());   // y‑pos a bit below edge
            cR.Modified();
            cR.Update();


            fs::path pngRun = root / "EMCal" / "invMassQA" / cutTag /
                              "Pi0Mass_Sigma_vs_Run_AllCentrality.png";
            ensure_dir(pngRun.parent_path());
            cR.SaveAs(pngRun.string().c_str());
            
            log(Lvl::INFO,"   • global run‑summary PNG → " + pngRun.string());

            std::ostringstream oss;
            oss << "\n──────────  RUN‑SUMMARY CONTENTS  ──────────\n"
                << "Slices plotted  : " << slicesDone.size() << '\n'
                << "Slice list      : ";
            for (const auto& s : slicesDone) oss << s << ' ';
            oss << "\nOutput directory : " << pngRun.parent_path().string()
                << "\nPNG files        :"
                << "\n   • Overlay   : " << pngRun.filename().string()
                << "\n   • Per‑slice : Pi0Mass_Sigma_vs_Run_AllCentrality.png"
                   " (one for each slice)\n";
            log(Lvl::INFO,oss.str());
        }
        catch(const std::exception& ex)
        {
            log(Lvl::ERR,
                std::string("writeRunSummary(): fatal exception – ") + ex.what());
        }

        log(Lvl::INFO,
            "writeRunSummary()  ⇦  completed\n"
            "══════════════════════════════════════════════════════════════");
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

static void tightenAxes(TH2* h, double padFrac = 0.05)
{
    if (!h) return;

    TAxis* axX = h->GetXaxis();
    TAxis* axY = h->GetYaxis();

    const int nBX = axX->GetNbins();
    const int nBY = axY->GetNbins();

    int hiX = nBX;
    while (hiX > 1 && h->Integral(hiX, nBX, 1, nBY) == 0) --hiX;

    int hiY = nBY;
    while (hiY > 1 && h->Integral(1, nBX, hiY, nBY) == 0) --hiY;

    const double xLow = axX->GetBinLowEdge(1);
    const double yLow = axY->GetBinLowEdge(1);

    double xUp  = axX->GetBinUpEdge(hiX);
    double yUp  = axY->GetBinUpEdge(hiY);

    const double dx = (xUp - xLow) * padFrac;
    const double dy = (yUp - yLow) * padFrac;

    axX->SetRangeUser(xLow, xUp + dx);
    axY->SetRangeUser(yLow, yUp + dy);

    /* keep titles nicely centred */
    axX->CenterTitle(true);   axY->CenterTitle(true);
    axX->SetTitleOffset(1.1F); axY->SetTitleOffset(1.45F);
}
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
    
 enum class Lvl { DBG, INFO, WARN, ERR };
 static inline void log(Lvl lvl, const std::string& msg)
 {
    static const char* tag[]{"DBG","INF","WRN","ERR"};
    std::ostream& os = (lvl == Lvl::ERR) ? std::cerr : std::cout;
    os << "[CorrQA] " << tag[static_cast<int>(lvl)] << "  " << msg << '\n';
 }
 public:
    /* ---------- static run-summary cache --------------------------- */
    using RunMap  = std::unordered_map<std::string, std::shared_ptr<TH2>>;
    using NameMap = std::unordered_map<std::string, RunMap>;
    
    CorrQA(const std::string& trig,
           const fs::path&    base,
           const CentList&    cent)
      : QA(trig, base, cent),
        m_trig(trig)
    {}

    ~CorrQA() override
    {
        writeCentralityOverviews();
        writeRunSummaries();
        writeCalorimeterSummary();
    }

 private:
     static std::unordered_map<std::string, NameMap> s_cache;
     std::string m_trig;

     /* ─ calorimeter‑overview helper ─ */
     std::unordered_map<std::string, std::shared_ptr<TH2>> m_calSummary;

    void writeRunSummaries()
    {
        const std::string pass = root.parent_path().filename().string();

        /* execute only once – during the Combined pass */
        if (pass != "Combined") {
            log(Lvl::DBG,"writeRunSummaries(): pass = \"" + pass +
                          "\" – skipped (only runs in Combined)");
            return;
        }
        if (s_cache.empty()) {
            log(Lvl::DBG,"writeRunSummaries(): s_cache empty – nothing to summarise");
            return;
        }

        using HVec   = std::vector<std::shared_ptr<TH2>>;
        using RunMap = std::unordered_map<std::string, HVec>;     // runID → vec
        std::unordered_map<std::string, RunMap> groupRun;         // groupDir → …

        /* ---------------- collect all TH2 clones --------------------- */
        for (auto& [groupDir, nameMap] : s_cache)
            for (auto& [hName, runMap] : nameMap)
                for (auto& [runID, h] : runMap)
                    if (runID != "Combined" && h)
                        groupRun[groupDir][runID].push_back(h);

        if (groupRun.empty()) {
            log(Lvl::WARN,"writeRunSummaries(): no non‑empty histogram sets found");
            return;
        }

        /* ---------------- static geometry constants ------------------ */
        constexpr int nCols = 8, nRows = 8;
        constexpr int canW  = nCols * 350, canH = nRows * 350;
        constexpr int perPage = nCols * nRows;

        log(Lvl::INFO,"writeRunSummaries(): starting – "
                      + std::to_string(groupRun.size()) + " detector‑pair folders");

        /* ---------------- iterate over folders ----------------------- */
        for (auto& [groupDir, runMap] : groupRun)
        {
            fs::path baseDir = root / "correlations" / groupDir;
            ensure_dir(baseDir);
            log(Lvl::INFO,"   ↳ folder \"" + groupDir + "\"  (" +
                          std::to_string(runMap.size()) + " runs)");

            using RunH = std::pair<std::string /*runID*/, std::shared_ptr<TH2>>;
            std::vector<RunH> runs;
            for (auto& [runID, vec] : runMap)
                if (!vec.empty())
                    runs.emplace_back(runID, vec.front());      // take the first map for the run

            if (runs.empty()) {
                log(Lvl::WARN,"      • folder \"" + groupDir + "\" – no runs, skipped");
                continue;
            }

            std::size_t page = 0;
            for (std::size_t idx = 0; idx < runs.size(); idx += perPage)
            {
                ++page;
                const std::size_t nThis = std::min<std::size_t>(perPage,
                                                                runs.size() - idx);

                log(Lvl::INFO,Form("      • page %zu  (%zu runs)", page, nThis));

                try {
                    TCanvas c(Form("c_%s_page%zu", groupDir.c_str(), page),
                              "", canW, canH);
                    c.Divide(nCols, nRows, 0.001, 0.001);

                    for (std::size_t i = 0; i < nThis; ++i) {
                        c.cd(static_cast<int>(i) + 1);
                        gPad->SetLogz();
                        tightenAxes(runs[idx + i].second.get());
                        runs[idx + i].second->Draw("COLZ");

                        /* tiny run‑label centred below each sub‑pad ---------------- */
                        TLatex lab; lab.SetNDC();
                        lab.SetTextFont(42);
                        lab.SetTextAlign(23);      // centred horizontally
                        lab.SetTextSize(0.032);
                        lab.DrawLatex(0.50, 0.02,
                                      stripLeadingZeros(runs[idx + i].first).c_str());
                    }

                    TLatex header; header.SetNDC();
                    header.SetTextFont(42);
                    header.SetTextAlign(11);
                    header.SetTextSize(0.034);
                    header.DrawLatex(0.02, 0.97, groupDir.c_str());

                    fs::path png = baseDir /
                        (std::string("page") + std::to_string(page) + ".png");
                    c.SaveAs(png.string().c_str());

                    log(Lvl::INFO,"         ↳ saved " + png.string());
                }
                catch (const std::exception& ex) {
                    log(Lvl::ERR,"         ✖ ROOT exception on page "
                                 + std::to_string(page) + " – " + ex.what());
                }
                catch (...) {
                    log(Lvl::ERR,"         ✖ unknown exception on page "
                                 + std::to_string(page));
                }
            }
        }

        /* clear cache to free memory */
        log(Lvl::DBG,"writeRunSummaries(): clearing s_cache");
        s_cache.clear();
    }
    
    /* ------------------------------------------------------------------ *
     *  2×3 calorimeter summary (EMCal/IHCal/OHCal  ×  MBD/sEPD)          *
     * ------------------------------------------------------------------ */
    void writeCalorimeterSummary()
    {
        static const std::array<std::string,6> order = {
            "EMCal_MBD","IHCal_MBD","OHCal_MBD",
            "EMCal_sEPD","IHCal_sEPD","OHCal_sEPD"};

        /* stop early if none of the six inclusive maps are cached */
        bool any = false;
        for (const auto& k : order)
            if (m_calSummary.count(k)) { any = true; break; }
        if (!any) {
            log(Lvl::DBG,"writeCalorimeterSummary(): nothing cached – skipped");
            return;
        }

        fs::path outDir = root / "correlations" / "calorimeterSummary";
        ensure_dir(outDir);

        TCanvas c("c_calSum","", 3*550, 2*500);
        c.Divide(3,2,0.001,0.001);

        /* ── loop over the 6 pads ─────────────────────────────────────────── */
        for (std::size_t i = 0; i < order.size(); ++i) {
            c.cd(static_cast<int>(i)+1);
            setupPad(gPad);
            gPad->SetTopMargin(0.12);   /* room for per‑pad header */
            gPad->SetLogz();

            auto it = m_calSummary.find(order[i]);
            if (it != m_calSummary.end()) {
                tightenAxes(it->second.get());
                it->second->Draw("COLZ");
            } else {                                /* N/A pad */
                TLatex miss; miss.SetNDC(); miss.SetTextAlign(22);
                miss.SetTextFont(42); miss.SetTextSize(0.04);
                miss.DrawLatex(0.50,0.50,"N/A");
            }

            /* ── per‑pad title, e.g. “EMCal vs MBD” ───────────────────────── */
            const std::size_t us = order[i].find('_');
            const std::string tokA = order[i].substr(0, us);
            const std::string tokB = order[i].substr(us+1);
            const std::string padTitle =
                prettyDet(tokA) + "  vs  " + prettyDet(tokB);

            TLatex padTx; padTx.SetNDC();
            padTx.SetTextAlign(22); padTx.SetTextFont(42);
            padTx.SetTextSize(0.045);
            padTx.DrawLatex(0.50, 0.96, padTitle.c_str());
        }

        /* ── global labels: two lines, smaller font, tight in the corner ─── */
        c.cd();  /* main canvas pad */
        const std::string runID   = stripLeadingZeros(root.parent_path()
                                                      .filename().string());
        const std::string trigLbl = prettifyTrigger(m_trig);

        TLatex tx; tx.SetNDC();
        tx.SetTextAlign(11);         /* left‑justified, top‑aligned           */
        tx.SetTextFont(42);
        tx.SetTextSize(0.022);       /* smaller than per‑pad titles           */

        constexpr double x0 = 0.01;  /* far left edge of canvas               */
        constexpr double y0 = 0.96; /* very top                              */

        tx.DrawLatex(x0, y0, ("Run " + runID).c_str());          /* line 1   */
        tx.DrawLatex(x0, y0 - 0.020, trigLbl.c_str());           /* line 2   */


        fs::path png = outDir / "calorimeterSummary.png";
        c.SaveAs(png.string().c_str());
        log(Lvl::INFO,"writeCalorimeterSummary(): saved " + png.string());

        m_calSummary.clear();
    }

    // ─────────────────────────────── 1. per-histogram ──────────────────────
    bool process(TObject* o) override
    {
        /* ------------------------------------------------------------------ *
         * Accept both 2‑D (“…_vs_…”) maps  *and*  the new 1‑D Δη / Δφ spectra *
         * ------------------------------------------------------------------ */
        if (!o->InheritsFrom(TH1::Class())) return false;     // neither TH1 nor TH2 → skip

        const std::string hName = o->GetName();

        /* Reject Jet‑QA, resolution profiles, event‑plane histos, … */
        static const std::vector<std::string> skipPrefixes = {
            "h_leadEt_",          // jet leading/subleading E_T maps
            "h_maxJetEt_",        // jet spectra
            "h_jetEt_",           // 3‑D jet histos
            "h_Psi",              // event‑plane Ψₙ histos
            "p_R",                // ALL sub‑event‑resolution profiles p_R1_, p_R2_, p_R3_
            "Minv_vs_"
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
            /* special case:  sEPD South × sEPD North  → store under  “sEPD/NS”   */
            const bool isSEPD_NS =
                (detA == "sEPD" && detB == "sEPD") &&
                (tokA.find("_South") != std::string::npos ||
                 tokB.find("_South") != std::string::npos) &&
                (tokA.find("_North") != std::string::npos ||
                 tokB.find("_North") != std::string::npos);

            if (isSEPD_NS)
                groupDir = "sEPD/NS";              /* yields …/correlations/sEPD/NS */
            else
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

            /* Build a descriptive title
             *   – always show the two detector systems in “A vs B” order
             *   – append centrality information:
             *        Inclusive   →  “(Inclusive)”
             *        X_Y         →  “(Cent X–Y %)”
             */
            std::string plotTitle = std::regex_replace(groupDir,
                                                       std::regex("_"), " vs ");

            if (!hasCent) {                                   // centrality‑integrated (“Inclusive”)
                plotTitle += "  (Inclusive)";
            } else {                                          // explicit centrality slice “lo_hi”
                std::smatch m;
                if (std::regex_match(slice, m, std::regex(R"((\d{1,3})_(\d{1,3}))"))) {
                    plotTitle += "  (Cent " + m[1].str() + "-" + m[2].str() + " %)";
                }
            }


            if (h2) {                               // 2‑D map
                h2->SetTitle("");

                auto axisLabel = [](const std::string& tok,
                                    const std::string& det) -> std::string
                {
                    auto hasTok = [](const std::string& str, const std::string& tag)->bool
                    {
                        const std::string p = "_" + tag;
                        std::size_t pos = str.find(p);
                        while (pos != std::string::npos) {
                            std::size_t end = pos + p.size();
                            if (end == str.size() || str[end] == '_') return true;
                            pos = str.find(p, pos + 1);
                        }
                        return false;
                    };

                    const bool south = hasTok(tok, "South") || hasTok(tok, "S");
                    const bool north = hasTok(tok, "North") || hasTok(tok, "N");

                    if (det == "EMCal") return "#SigmaE_{CEMC}  [GeV]";
                    if (det == "IHCal") return "#SigmaE_{IHCal} [GeV]";
                    if (det == "OHCal") return "#SigmaE_{OHCal} [GeV]";
                    if (det == "HCal")  return "#SigmaE_{HCal}  [GeV]";
                    if (det == "MBD")   return "#SigmaQ_{MBD}   [Charge]";

                    if (det == "sEPD") {
                        if (south) return "#SigmaQ_{South} [Charge]";
                        if (north) return "#SigmaQ_{North} [Charge]";
                        return "#SigmaQ_{sEPD}  [Charge]";
                    }
                    return det;                                   // fallback
                };

                h2->GetXaxis()->SetTitle( axisLabel(tokA, detA).c_str() );
                h2->GetYaxis()->SetTitle( axisLabel(tokB, detB).c_str() );

                c.SetLogz();
                tightenAxes(h2);
                h2->Draw("COLZ");

                /* one centred header above the pad --------------------------- */
                TLatex tl;
                tl.SetNDC();
                tl.SetTextAlign(22);
                tl.SetTextFont(42);
                tl.SetTextSize(0.04);
                std::string header;
                if (detA == "sEPD" && detB == "sEPD")            // special N‑S correlation
                    header = "sEPD North - South Correlations";
                else
                    header = prettyDet(tokA) + "  vs  " + prettyDet(tokB);

                /* append pretty trigger in round brackets */
                header += "  (" + prettifyTrigger(m_trig) + ")";

                tl.DrawLatex(0.50, 0.96, header.c_str());
            }
            else {                                  // 1‑D spectrum (Δη / Δφ)
                h1->SetTitle("");

                /* figure‑out which variable we are drawing ------------------- */
                const char* xlab =
                    (hName.rfind("h_dEta_", 0) == 0) ? "#Delta#eta"
                                                     : "#Delta#phi";
                h1->GetXaxis()->SetTitle( xlab );
                h1->GetYaxis()->SetTitle( "Events" );

                tightenAxes(h1);
                h1->Draw();

                TLatex tl;
                tl.SetNDC();
                tl.SetTextAlign(22);
                tl.SetTextFont(42);
                tl.SetTextSize(0.04);
                std::string header = prettyDet(tokA) + "  vs  " + prettyDet(tokB) +
                                     "  (" + prettifyTrigger(m_trig) + ")";
                tl.DrawLatex(0.50, 0.96, header.c_str());
            }
            drawRunLabel( stripLeadingZeros(root.parent_path().filename().string()) );
            c.SaveAs(pngFile.string().c_str());

            /* cache six inclusive maps for the 2×3 calorimeter summary */
            if (!hasCent && h2) {
                static const std::unordered_set<std::string> want = {
                    "EMCal_MBD","IHCal_MBD","OHCal_MBD",
                    "EMCal_sEPD","IHCal_sEPD","OHCal_sEPD"};
                if (want.count(groupDir)) {
                    auto cl = std::shared_ptr<TH2>(static_cast<TH2*>(h2->Clone()));
                    cl->SetDirectory(nullptr); tidyAxes(cl.get()); styleAxes(cl.get(), false);
                    m_calSummary[groupDir] = std::move(cl);
                }
            }

            /* ─── mirror copy for unified “HCal” collection ───────────────────── */
            const bool involvesHCal =
                (groupDir.find("IHCal")     != std::string::npos) ||
                (groupDir.find("OHCal")     != std::string::npos) ||
                (groupDir.find("totalHCal") != std::string::npos);

            if (involvesHCal)
            {
                /* partner detector = first token in groupDir that does NOT contain “HCal” */
                std::string partner;
                {
                    std::stringstream ss(groupDir);
                    std::string tok;
                    while (std::getline(ss, tok, '_'))
                        if (tok.find("HCal") == std::string::npos) { partner = tok; break; }
                }

                if (!partner.empty())
                {
                    fs::path hcalDir = root / "correlations" / "HCal" / ("HCal_" + partner);
                    if (hasCent) hcalDir /= ("Cent_" + slice);        // keep centrality hierarchy
                    ensure_dir(hcalDir);

                    fs::path hcalPng = hcalDir / (o->GetName() + std::string(".png"));
                    c.SaveAs(hcalPng.string().c_str());               // second write
                }
            }
            /* ─────────────────────────────────────────────────────────────────── */

        }

        /* ------------------------------------------------------------------
         * D.  Collect centrality‑dependent clones for the overview canvas
         *     (only meaningful for 2‑D histograms)
         * ----------------------------------------------------------------- */
        if (hasCent && h2)                                     // <<< guard against 1‑D
        {
            const std::string baseKey = groupDir + "|" +
                                        stripCentSuffix(hName);  // "h_SEPD_vs_MBD"
            auto& v = m_centCache[baseKey];
            auto  cl = std::shared_ptr<TH2>(static_cast<TH2*>(h2->Clone()));
            cl->SetDirectory(nullptr);
            tidyAxes(cl.get());
            styleAxes(cl.get(), false);
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

    /* Pretty detector label that keeps possible North/South qualifiers */
    static std::string prettyDet(const std::string& tok)
    {
        std::string det = canonicalDet(tok);    // EMCal / sEPD / …

        /* match “…_North”, “…_South”, “…_N”, “…_S” only when they form
           a complete underscore‑separated token                               */
        auto hasToken = [](const std::string& str, const std::string& tag) -> bool
        {
            const std::string pat = "_" + tag;
            std::size_t pos = str.find(pat);
            while (pos != std::string::npos)
            {
                const std::size_t end = pos + pat.size();
                if (end == str.size() || str[end] == '_') return true; // exact token
                pos = str.find(pat, pos + 1);
            }
            return false;
        };

        const bool isN = hasToken(tok, "North") || hasToken(tok, "N");
        const bool isS = hasToken(tok, "South") || hasToken(tok, "S");

        if (isN)      det += " North";
        else if (isS) det += " South";
        return det;
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

            /* ------------------------------------------------------------ *
             * 1.  Produce a single “North + South” map by a straightforward
             *     bin‑by‑bin ADDITION, preserving statistics.
             * ------------------------------------------------------------ */
            std::unique_ptr<TH2> hTot(
                static_cast<TH2*>(pair.s->Clone(
                    (baseKey + std::string("_combined")).c_str())));
            hTot->Add(pair.n.get());                    // ⟵ real addition
            tidyAxes(hTot.get());  styleAxes(hTot.get(), false);

            fs::path pngTot = dir / (baseKey + std::string("_combined.png"));
            {
                TCanvas cTot(("c_" + baseKey + "_comb").c_str(), "", 1100, 800);
                setupPad(&cTot);
                cTot.SetLogz();
                tightenAxes(hTot.get());
                hTot->Draw("COLZ");
                drawRunLabel(stripLeadingZeros(root.parent_path().filename().string()));
                cTot.SaveAs(pngTot.string().c_str());
            }

            /* ------------------------------------------------------------ *
             * 2.  Keep the original side‑by‑side view for quick checks.
             * ------------------------------------------------------------ */
            fs::path pngNS = dir / (baseKey + std::string("_NS.png"));

            const double zMax = std::max(pair.n->GetMaximum(),
                                          pair.s->GetMaximum());
            pair.n->SetMaximum(zMax);  pair.s->SetMaximum(zMax);
            pair.n->SetMinimum(1);     pair.s->SetMinimum(1);

            TCanvas c("c_ns","",1200,600); c.Divide(2,1,0.01,0.01);
            c.cd(1); setupPad(gPad); gPad->SetLogz();
            tightenAxes(pair.s.get());
            pair.s->DrawCopy("COLZ");

            c.cd(2); setupPad(gPad); gPad->SetLogz();
            tightenAxes(pair.n.get());
            pair.n->DrawCopy("COLZ");

            drawRunLabel(stripLeadingZeros(root.parent_path().filename().string()));
            c.SaveAs(pngNS.string().c_str());

            g_nsCache.erase(cacheKey);                  // clean‑up
        }

    }

    /* ==================================================================== *
     * §2  Combined run‑summary cache – verbose, exception‑safe             *
     * ==================================================================== */
    void cacheForRunSummary(const std::string& groupDir,
                            const std::string& hName,
                            TObject* o)
    {
        const std::string runID = root.parent_path().filename().string();

        /* never cache the synthetic “Combined” pass itself */
        if (runID == "Combined") {
            log(Lvl::DBG,"cacheForRunSummary(): Combined pass – histogram \"" +
                          hName + "\" ignored");
            return;
        }

        /* guard against non‑TH2 objects – should never happen, but better safe */
        if (!o || !o->InheritsFrom(TH2::Class())) {
            const char* what = o ? o->IsA()->GetName() : "nullptr";
            log(Lvl::WARN,"cacheForRunSummary(): object \"" + hName +
                           "\" is " + what + ", expected TH2 – skipped");
            return;
        }

        try {
            auto* cl = static_cast<TH2*>(o->Clone());

            if (!cl) {
                log(Lvl::ERR,"cacheForRunSummary(): Clone() returned nullptr – \"" +
                              hName + "\" not cached");
                return;
            }

            cl->SetDirectory(nullptr);            /* detach from any TDirectory  */
            cl->SetBit(kCanDelete,false);         /* ROOT ownership protection   */

            tidyAxes(cl);
            styleAxes(cl,true);

            /* store in 3‑level cache:  detPair ▸ histName ▸ runID             */
            s_cache[groupDir][hName][runID].reset(cl);

            log(Lvl::DBG,"cacheForRunSummary(): cached \"" + hName +
                          "\" for run " + runID + " under \"" + groupDir + '"');
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR,"cacheForRunSummary(): std::exception while cloning \"" +
                          hName + "\" – " + ex.what());
        }
        catch (...) {
            log(Lvl::ERR,"cacheForRunSummary(): unknown exception while cloning \"" +
                          hName + '"');
        }
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

        /* ─── mirror copy for unified “HCal” collection ───────────────────── */
        {
            fs::path hDir = root / "correlations" / "HCal" / ("HCal_" + otherDet);
            if (hasCent) hDir /= ("Cent_" + slice);               // preserve cent‑subfolders
            ensure_dir(hDir);

            fs::path hPng = hDir / (canonName + std::string(".png"));
            cTot.SaveAs(hPng.string().c_str());                   // duplicate write
        }
        /* ─────────────────────────────────────────────────────────────────── */

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

            /* leave room for a one‑line header */
            c.SetTopMargin(0.12);

            /* split the canvas into regular pads */
            c.Divide(nCols, nRows, 0.001, 0.001);

            /* ── global header centred above the grid ───────────────────────── */
            {
                c.cd();                       /* main canvas pad (not a sub‑pad) */
                /* centred title -------------------------------------------------- */
                TLatex titleTx;  titleTx.SetNDC();
                titleTx.SetTextAlign(22);           /* centred */
                titleTx.SetTextFont(42);
                titleTx.SetTextSize(0.038);

                /* build the title from the *original* detector tokens so that any
                 * “North / South” qualifier survives the canonicalisation step      */
                std::string tokA, tokB;
                {
                    /* baseHist still carries the full detector tag, e.g.
                     *   h_SEPD_N_vs_CEMC_North    or   h_SEPD_S_vs_CEMC_South         */
                    const std::size_t vsPos = baseHist.find("_vs_");
                    if (vsPos != std::string::npos) {
                        tokA = baseHist.substr(2, vsPos - 2);   // drop leading “h_”
                        tokB = baseHist.substr(vsPos + 4);      // text after “_vs_”
                    } else {                                    // very rare fall‑back
                        const std::size_t us = groupDir.find('_');
                        tokA = groupDir.substr(0, us);
                        tokB = groupDir.substr(us + 1);
                    }
                }

                const std::string title =
                    prettyDet(tokA) + "  vs  " + prettyDet(tokB) +
                    "  -  Centrality overview";
                
                titleTx.DrawLatex(0.50, 0.97, title.c_str());

                /* single, larger run label (top‑left, once per canvas) ----------- */
                TLatex runTx; runTx.SetNDC();
                runTx.SetTextAlign(11);             /* left‑aligned */
                runTx.SetTextFont(42);
                runTx.SetTextSize(0.042);
                runTx.DrawLatex(0.02, 0.95,
                    ("Run " + stripLeadingZeros(root.parent_path().filename().string())).c_str());
            }

            for (int i = 0; i < n; ++i) {
                c.cd(i+1);  setupPad(gPad);
                gPad->SetTopMargin(0.15);
                gPad->SetLogz();
                tightenAxes(vec[i].second.get());
                vec[i].second->Draw("COLZ");
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
                        oss << "Centrality: " << lo << " #minus " << hi << " %";
                        label = oss.str();
                    } else {
                        label = vec[i].first;            // fallback – unexpected slice key
                    }
                }
                tl.DrawLatex(0.4, 0.18, label.c_str());
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

    fs::path out = root / "EMCal" / "otherQA" / "EMCalHitMap_AllCentrality.png";
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

    auto makePanel = [&](TH2* src,
                           const std::string& triggerName) -> std::unique_ptr<TH2F>
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
        
      rot->SetTitle("");    // canvas title left blank – header will be drawn later

      // ------------------------------------------------------------------
      // prepare header strings (used after the canvas exists)
      // ------------------------------------------------------------------
      std::string runID = root.parent_path().filename().string();
      if (std::all_of(runID.begin(), runID.end(), ::isdigit))
            runID = std::to_string(std::stoi(runID));              // strip leading zeros
      std::string runLabel  = "run " + runID;
      std::string trigLabel = prettifyTrigger(triggerName);


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

      /* draw header once a pad is active */
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize(0.028);          // smaller font, matches new HCal style
      tl.SetTextAlign(13);

      Long64_t nEvents = src->GetEntries();        // total entries in this map
      tl.DrawLatex(0.08, 0.99,
                     Form("EMCal Hitmap for %s (%s, nEvents = %lld)",
                          trigLabel.c_str(),
                          runLabel.c_str(),
                          static_cast<long long>(nEvents)));
      gPad->Update();

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

      /* ----------------------------------------------------------------
       *  missing‑SEB detector (run‑by‑run, Inclusive slice only)
       * ---------------------------------------------------------------- */
      {
            static std::unordered_set<std::string> s_doneRuns;   // run‑ID cache
            const std::string runKey = root.parent_path().filename().string();
            const bool isCombined    = (runKey == "Combined");

            if (!isCombined && slice == "Inclusive" && s_doneRuns.insert(runKey).second)
            {
                /* mapping  SEB → contiguous sector range (inclusive) */
                const std::pair<int,int> sebRange[16] = {
                    {16,19},  /* SEB0  */ {20,23},  /* SEB1  */ {28,31},  /* SEB2  */
                    {24,27},  /* SEB3  */ { 8,11},  /* SEB4  */ {12,15},  /* SEB5  */
                    { 4, 7},  /* SEB6  */ { 0, 3},  /* SEB7  */ {48,51},  /* SEB8  */
                    {52,55},  /* SEB9  */ {60,63},  /* SEB10 */ {56,59},  /* SEB11 */
                    {40,43},  /* SEB12 */ {44,47},  /* SEB13 */ {36,39},  /* SEB14 */
                    {32,35}   /* SEB15 */
                };

                /* helper – is *entire* sector empty? */
                auto sectorIsEmpty = [&](int sec)->bool
                {
                    double sum = 0.0;
                    for (int iphi = 0; iphi < 256; ++iphi)
                        for (int ieta = 0; ieta < 96; ++ieta)
                            if (sector_from_idx(ieta, iphi) == sec)
                                sum += src->GetBinContent(iphi + 1, ieta + 1);
                    return sum <= 0.0;
                };

                std::vector<int> missingSEB;
                for (int seb = 0; seb < 16; ++seb)
                {
                    const auto [lo, hi] = sebRange[seb];
                    bool empty = true;
                    for (int sec = lo; sec <= hi && empty; ++sec)
                        empty &= sectorIsEmpty(sec);
                    if (empty) missingSEB.push_back(seb);
                }

                if (!missingSEB.empty())
                {
                    const fs::path txt = baseDir / "output" / "MissingSEB.txt";
                    ensure_dir(txt.parent_path());

                    std::ofstream ofs(txt, std::ios::app);
                    ofs << runLabel;                       // first column = run number (no leading zeros)
                    for (int seb : missingSEB) ofs << '\t' << "SEB" << seb;
                    ofs << '\n';
              }
          }
      }

      c.SaveAs(outPng.string().c_str());
      return rot;

    }; // makePanel

    // ── centralities: keep one rotated copy per slice for the overview
    if (isMap && o->InheritsFrom(TH2::Class())) {
      auto p = makePanel(static_cast<TH2*>(o),
                           root.filename().string());   // raw trigger name
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
 /* ---------- global run‑level bookkeeping -------------------- */
 struct Range { int phiMin, phiMax, etaMin, etaMax; };
 static inline std::unordered_set<std::string>              s_runIDs;      // all runs seen
 static inline std::unordered_map<std::string,Range>        s_missing;     // run → first empty range
 static inline std::unordered_set<std::string>              s_scanned;     // guard: check once/run

 public:
  using QA::QA;

  static void writeMissingBinReport(const fs::path& txtPath);
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

    /* ---------- NEW: remember this run ID for the combined pass -------- */
    {
        std::string id = root.parent_path().filename().string();
        if (id != "Combined") HcalQA::s_runIDs.insert(id);
    }

      const std::string slice = sliceKey(n);                 // Inclusive / Cent_x_y
      const bool isMap        = n.find("_EtaPhiMap_") != std::string::npos;

      /* ---------- one‑time missing‑bin scan per run --------------- */
      if (isMap) {
          std::string runID = root.parent_path().filename().string();   // e.g. "0068408"
          if (runID != "Combined" && !s_scanned.count(runID)) {
              s_scanned.insert(runID);

              const TH2* h2 = static_cast<const TH2*>(o);
              const int nPhi = h2->GetNbinsX();
              const int nEta = h2->GetNbinsY();

              int phiMin =  999, phiMax = -1, etaMin =  999, etaMax = -1;
              for (int ip = 1; ip <= nPhi; ++ip)
                  for (int ie = 1; ie <= nEta; ++ie)
                      if (h2->GetBinContent(ip, ie) == 0) {
                          phiMin = std::min(phiMin, ip - 1);
                          phiMax = std::max(phiMax, ip - 1);
                          etaMin = std::min(etaMin, ie - 1);
                          etaMax = std::max(etaMax, ie - 1);
                      }

              if (phiMax >= 0 && etaMax >= 0)          // at least one empty cell detected
                  s_missing[runID] = {phiMin, phiMax, etaMin, etaMax};
          }
      }

      ulog::trace("HcalQA  → processing \"" + n +
               "\"  slice=" + slice +
               (isMap ? "  (map)" : "  (scalar)"));

    /* -------------------------------------------------------------- *
     *  1.  Lambda that renders one η–φ map panel                     *
     * -------------------------------------------------------------- */
    auto makePanel = [&](TH2* src,
                           const fs::path& outPng,
                           const std::string& triggerName)   // NEW: pass prettified trigger
    {
      try {
          ensure_dir(outPng.parent_path());
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
          rot->GetXaxis()->SetTitle("#eta");
          rot->GetYaxis()->SetTitle("#phi");

          // 1.3 canvas & save ---------------------------------------
          const int cw = (kHCalCanvasW > 0) ? kHCalCanvasW : nEta * px;
          const int ch = (kHCalCanvasH > 0) ? kHCalCanvasH : nPhi * px;

          TCanvas c("c_hcal", "", cw, ch);
          c.SetRightMargin(0.16);  c.SetLeftMargin(0.08);
          c.SetBottomMargin(0.08); c.SetTopMargin(0.055);
          c.SetFixedAspectRatio();

          rot->SetTitle("");
          rot->Draw("COLZ");
          
          /* ── add horizontal z-axis title (“counts”) ───────────────────────── */
          gPad->Update();                                   // palette now exists
          if (auto* pal = dynamic_cast<TPaletteAxis*>(
                  rot->GetListOfFunctions()->FindObject("palette")))
          {
              // 1. (optional) fine-tune palette position so we know exact NDC
              //    — keeps label from overlapping the plot frame
              pal->SetX1NDC(0.88);           // left  edge of palette
              pal->SetX2NDC(0.92);           // right edge of palette
              pal->SetY1NDC(0.12);           // bottom
              pal->SetY2NDC(0.92);           // top
              gPad->Modified(); gPad->Update();

              // 2. draw the horizontal label
              const double xMid = 0.5 * (pal->GetX1NDC() + pal->GetX2NDC());
              const double yTop = pal->GetY2NDC() + 0.02;   // a little above

              TLatex tz;
              tz.SetNDC();
              tz.SetTextSize(0.032);          // match your header size
              tz.SetTextAlign(22);            // centre, centred
              tz.DrawLatex(xMid, yTop, "counts");
          }

          /* ── descriptive title – plot type + run information ───────────── */
          std::string runID = root.parent_path().filename().string();          // "00067834" or "Combined"
          std::string runLabel;

          if (runID == "Combined")
          {
              std::vector<int> rNums;
              for (const auto& r : HcalQA::s_runIDs)
                  if (std::all_of(r.begin(), r.end(), ::isdigit))
                      rNums.push_back(std::stoi(r));

              if (!rNums.empty()) {
                  const auto [lo, hi] = std::minmax_element(rNums.begin(), rNums.end());

                  // TLatex escape sequence: "#rightarrow"
                  runLabel = Form("runs %d #rightarrow %d, %zu runs",
                                  *lo, *hi, rNums.size());
              } else {
                  runLabel = "combined";
              }
          }

          else  /* ordinary single run */
          {
              std::string cleanID = runID;
              if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                  cleanID = std::to_string(std::stoi(runID));   // drop leading 0s
              runLabel = "run " + cleanID;
          }

          /* what kind of map is this?  – inspect the histogram name */
          const std::string hName = src->GetName();
          const char* mapType =
                  (hName.find("_tot")  != std::string::npos) ? "TotalHCal Hitmap" :
                  (hName.find("IHCAL") != std::string::npos) ? "IHCal Hitmap"     :
                  (hName.find("OHCAL") != std::string::npos) ? "OHCal Hitmap"     :
                                                               "TotalHCal Hitmap";

          TLatex tl;
          tl.SetNDC();
          tl.SetTextSize(0.028);          // smaller font
          tl.SetTextAlign(13);

          std::string trigLabel  = prettifyTrigger(triggerName);
          Long64_t    nEvents    = src->GetEntries();          // histogram entry count

          tl.DrawLatex(0.08, 0.98,
                       Form("%s for %s (%s, nEvents = %lld)",
                            mapType,
                            trigLabel.c_str(),
                            runLabel.c_str(),
                            static_cast<long long>(nEvents)));

          gPad->Update();
          
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
          ulog::trace("HcalQA  → wrote " + outPng.string());
      }
      catch (const std::exception& ex) {
          ulog::warn(std::string("HcalQA  WARN  failed to save panel for \"")
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
        makePanel(static_cast<TH2*>(o), out, root.filename().string());
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
          ~Pair() {
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
          ulog::warn("HcalQA  WARN  clone failed for \"" + n + "\": " + ex.what());
          return true;
      }

        // build a *true* total‑HCal map once both IHCal & OHCal are ready
        if (p.i && p.o)
        {
            /* ------------------------------------------------------------ *
             * 1. create an empty clone that will hold the summed occupancies
             *    – we inherit binning & axes from IHCal, but clear contents
             * ------------------------------------------------------------ */
            auto tot = std::unique_ptr<TH2>(
                cloneDetach(p.i.get(), (baseName + "_tot").c_str()));
            tot->Reset("ICE");                           // zero all bins

            /* ------------------------------------------------------------ *
             * 2. loop over every bin and add IHCal+OHCal ONLY
             *    when neither side carries the “‑9999” masked‑plate marker
             * ------------------------------------------------------------ */
            const int nX = tot->GetNbinsX();
            const int nY = tot->GetNbinsY();

            for (int ix = 1; ix <= nX; ++ix)
                for (int iy = 1; iy <= nY; ++iy)
                {
                    const double cI = p.i->GetBinContent(ix, iy);
                    const double cO = p.o->GetBinContent(ix, iy);

                    if (cI < 0 || cO < 0) {               // at least one bad plate
                        tot->SetBinContent(ix, iy, -9999.);   // keep it masked
                    } else {
                        tot->SetBinContent(ix, iy, cI + cO);  // proper sum
                    }
                }

            tot->SetTitle("totalHCal");                   // picked up by makePanel()

            /* ------------------------------------------------------------ *
             * 3. write the PNG through the same helper as IHCal / OHCal
             * ------------------------------------------------------------ */
            fs::path outTot = cPath(root, slice, fs::path("HCal") / "totalHcal")
                                / (baseName + "_tot.png");
            makePanel(tot.get(), outTot, root.filename().string());

            cache.erase(key);                             // drop pair – not needed any more
        }
    }

    return true;
  }
};


void HcalQA::writeMissingBinReport(const fs::path& txtPath)
{
    ensure_dir(txtPath.parent_path());
    std::ofstream rep(txtPath);
    if (!rep) {
        ulog::warn("HcalQA  WARN  cannot open " + txtPath.string() + " for writing");
        return;
    }

    rep << "# Runs with empty IHCal / OHCal η–φ bins\n"
        << "# runID\tphi[min‑max]\teta[min‑max]\n";

    for (const auto& [run,r] : s_missing)
        rep << run << '\t'
            << r.phiMin << '-' << r.phiMax << '\t'
            << r.etaMin << '-' << r.etaMax << '\n';

    ulog::info("HcalQA  → wrote missing‑bin report to " + txtPath.string());
}



// ───────────────── Event‑plane observables (sEPD) ──────────────────────
class sEPDotherQA : public QA
{
 public:
  sEPDotherQA(std::string t, fs::path b, const CentList& s)
  : QA(std::move(t), std::move(b), s) {}

    bool process(TObject* o) override
    {
      if (!o) { ulog::err("[sEPDotherQA] nullptr TObject – skip"); return false; }
      if (!o->InheritsFrom(TH1::Class())) return false;

      const std::string n = o->GetName();
      ulog::trace("[sEPDotherQA] Inspecting \"" + n + "\"");

      static const std::vector<std::string> keys = {
          /* event‑plane QA – now explicit South & North histograms */
          "h_Psi1_sEPD_S", "h_Psi1_sEPD_N",
          "h_Psi2_sEPD_S", "h_Psi2_sEPD_N",
          "h_Psi3_sEPD_S", "h_Psi3_sEPD_N",
          "h_Psi1_res_vs_Qsum", "h_Psi2_res_vs_Qsum", "h_Psi3_res_vs_Qsum",
          "p_R1_vs_cent", "p_R2_vs_cent", "p_R3_vs_cent",
          /* ring‑level QA (unchanged) */
          "h_SEPD_RingOcc_South", "h_SEPD_RingOcc_North",
          "h_SEPD_RingQ_South",   "h_SEPD_RingQ_North"
      };
      const bool match = std::any_of(keys.begin(), keys.end(),
                                     [&](const std::string& k){ return n.rfind(k,0)==0; });
      if (!match) return false;

      /* 2. decide sub‑folder – tileQA vs EventPlaneQA -------------------- */
      const bool isTile =
           n.rfind("h_SEPD_RingOcc_",0)==0 || n.rfind("h_SEPD_RingQ_",0)==0;
      const std::string subdir = isTile ? "sEPD/tileQA" : "sEPD/EventPlaneQA";

      /* 3. build output path & create folders ---------------------------- */
      fs::path out = cPath(root, sliceKey(n), subdir) / (n + ".png");
      try { ensure_dir(out.parent_path()); }
      catch (const std::exception& e) {
        ulog::err("[sEPDotherQA] mkdir failed for \"" +
                 out.parent_path().string() + "\": " + e.what());
        return false;
      }

        /* ------------------------------------------------------------------
         *  Tile‑by‑tile QA
         *  ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑
         *  Collect the four ring histograms (Occ/Q  ×  South/North) that
         *  belong to the same trigger + centrality slice.  When all four
         *  are available create **one** summary canvas with **two pads only**:
         *
         *    pad 1 :  Occupancy  (South + North overlaid, legend)
         *    pad 2 :  ΣQ         (South + North overlaid, legend)
         *
         *  Non‑tile histograms (event‑plane QA) are still saved with save1D().
         * ------------------------------------------------------------------ */
        if (isTile)
        {
            const std::string slice   = sliceKey(n);                    // "Inclusive", "0_10", …
            const std::size_t usPos   = n.rfind('_');
            const std::string trigger = (usPos==std::string::npos) ? "UNKNOWN"
                                     : n.substr(usPos+1);               // last field
            const std::string key     = slice + "|" + trigger;

            enum { OCC_S, OCC_N, Q_S, Q_N };
            static std::unordered_map<std::string,
                                      std::array<std::shared_ptr<TH1>,4>> cache;

            auto& set = cache[key];
            if      (n.find("RingOcc_South")!=std::string::npos) set[OCC_S].reset( static_cast<TH1*>(o->Clone()) );
            else if (n.find("RingOcc_North")!=std::string::npos) set[OCC_N].reset( static_cast<TH1*>(o->Clone()) );
            else if (n.find("RingQ_South")  !=std::string::npos) set[Q_S  ].reset( static_cast<TH1*>(o->Clone()) );
            else if (n.find("RingQ_North")  !=std::string::npos) set[Q_N  ].reset( static_cast<TH1*>(o->Clone()) );

            if (!std::all_of(set.begin(), set.end(),
                             [](const auto& p){ return bool(p); }))
                return true;                     // still waiting for missing arms

            auto prep = [](TH1* h, int col)
            {
                h->SetDirectory(nullptr);
                h->SetLineColor(col);
                h->SetMarkerColor(col);
                h->SetLineWidth(2);
                h->SetStats(0);
            };

            /* ---------- build 1×2 summary canvas (bottom pads removed) ----- */
            TCanvas c("c_ringQA", "sEPD ring QA", 1200, 600);
            c.Divide(2,1,0.02,0.02);                       // two pads only

            const int clrS = kRed+1,  clrN = kBlue+2;

            /* pad 1 : Occupancy overlay ------------------------------------ */
            c.cd(1); gPad->SetGridy();
            prep(set[OCC_S].get(), clrS); prep(set[OCC_N].get(), clrN);

            double occIntS = set[OCC_S]->Integral();
            double occIntN = set[OCC_N]->Integral();
            if (occIntS > 0.) set[OCC_S]->Scale(1.0 / occIntS);
            if (occIntN > 0.) set[OCC_N]->Scale(1.0 / occIntN);

            set[OCC_S]->SetTitle("");
            set[OCC_S]->GetYaxis()->SetTitle("normalised hits");
            set[OCC_S]->Draw("hist");
            set[OCC_N]->Draw("hist same");
            TLegend leg1(0.55,0.2,0.88,0.4); leg1.SetBorderSize(0);
            leg1.AddEntry(set[OCC_S].get(),"South","l");
            leg1.AddEntry(set[OCC_N].get(),"North","l");
            leg1.Draw();

            /* pad 2 : ΣQ overlay ------------------------------------------- */
            c.cd(2); gPad->SetGridy();
            prep(set[Q_S].get(), clrS); prep(set[Q_N].get(), clrN);

            double qIntS = set[Q_S]->Integral();
            double qIntN = set[Q_N]->Integral();
            if (qIntS > 0.) set[Q_S]->Scale(1.0 / qIntS);
            if (qIntN > 0.) set[Q_N]->Scale(1.0 / qIntN);

            set[Q_S]->SetTitle("");
            set[Q_S]->GetYaxis()->SetTitle("normalised #SigmaQ");
            set[Q_S]->Draw("hist");
            set[Q_N]->Draw("hist same");
            TLegend leg2(0.55,0.2,0.88,0.4); leg2.SetBorderSize(0);
            leg2.AddEntry(set[Q_S].get(),"South","l");
            leg2.AddEntry(set[Q_N].get(),"North","l");
            leg2.Draw();

            /* global header ------------------------------------------------- */
            c.cd(0);
            /* header: detector + trigger + run + events */
            std::string runID = root.parent_path().filename().string();
            if (std::all_of(runID.begin(), runID.end(), ::isdigit))
                runID = std::to_string(std::stoi(runID));              // strip leading zeros

            unsigned long long nEvt =
                static_cast<unsigned long long>( set[OCC_S]->GetEntries()
                                               + set[OCC_N]->GetEntries() );

            std::string trigLabel = prettifyTrigger(trigger);

            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(13);
            tl.SetTextFont(42);
            tl.SetTextSize(0.028);          // smaller, consistent with HCal/EMCal
            tl.DrawLatex(0.08, 0.96,
                         Form("sEPD Ring QA for %s (run %s, nEvents = %llu)",
                              trigLabel.c_str(),
                              runID.c_str(),
                              static_cast<long long>(nEvt)));

            /* output -------------------------------------------------------- */
            fs::path dst = cPath(root, slice, "sEPD/tileQA")
                         / ("sEPD_RingQA_" + slice + "_" + trigger + ".png");
            ensure_dir(dst.parent_path());
            c.SaveAs(dst.string().c_str());
            ulog::ok("[sEPDotherQA] Saved combined tile QA → " + dst.string());

            cache.erase(key);                               // free memory
            return true;
        }


        /* ------------------------------------------------------------------
         *  Event‑plane (non‑tile) histograms – keep original behaviour
         * ----------------------------------------------------------------- */
        try {
          save1D(static_cast<TH1*>(o), out);
          ulog::ok("[sEPDotherQA] Saved → " + out.string());
        } catch (const std::exception& e) {
          ulog::err("[sEPDotherQA] save1D failed for \"" + n + "\": " + e.what());
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
      ulog::err("[NSDetectorQA] nullptr TObject received – skipping.");
      return false;
    }
    if (!o->InheritsFrom(TH1::Class())) return false;

    const std::string hName = o->GetName();
    if (!DERIVED::accept(hName)) return false;

    const std::string slice = sliceKey(hName);
    const bool south = hName.find("_South_") != std::string::npos;
    ulog::trace("[NSDetectorQA] " + std::string(south ? "South" : "North") +
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
        ulog::ok("[NSDetectorQA] Saved 1‑D histo → " + out.string());
      } catch (const std::exception& e) {
        ulog::err("[NSDetectorQA] Failed to save \"" + hName + "\": " + e.what());
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
      ulog::err("[NSDetectorQA] Clone failed for \"" + hName + "\" – skipping.");
      return false;
    }
    hClone->SetDirectory(nullptr);
    hClone->SetStats(0);

    south ? mp.s = hClone : mp.n = hClone;
    if (!_cache.ready(trig + slice)) {
      ulog::trace("[NSDetectorQA] Waiting for partner arm to arrive ("
                 + slice + ").");
      return true;                         // partner not yet seen
    }

    MapPair in = _cache.pop(trig + slice); // take ownership
    auto tidy = [](TH2* h)
    {
            /* smallest positive bin – later used for colour‑scale lower bound */
            double zMin = std::numeric_limits<double>::max();
            for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
                for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                    const double z = h->GetBinContent(ix,iy);
                    if (z > 0.0 && z < zMin) zMin = z;
                }

            /* — cosmetic: hide per‑bin outlines for sEPD (regular TH2),
                 but keep thin borders for MBD hex‑maps (TH2Poly) — */
            if (h->InheritsFrom(TH2Poly::Class())) {
                h->SetLineColor(kBlack);          // hexagon edges remain
                h->SetLineWidth(1);
            } else {
                h->SetLineColor(0);               // invisible → no “squiggly” rings
                h->SetLineWidth(0);
            }

            h->GetZaxis()->SetNdivisions(510);    // major + minor ticks on colour bar
      };
      tidy(in.s); tidy(in.n);

      /* ------------------------------------------------------------- *
       * Robust colour‑scale upper bound:
       *   – collect all positive bin contents
       *   – take the 99.5‑percentile (ignores the top 0.5 % outliers)
       *   – add 5 % head‑room
       * ------------------------------------------------------------- */
      auto quantile99 = [](TH2* h) -> double
      {
          std::vector<double> v;                          // positive entries only
          v.reserve(h->GetNbinsX() * h->GetNbinsY());

          for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
              for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                  double z = h->GetBinContent(ix, iy);
                  if (z > 0.0) v.push_back(z);
              }

          if (v.empty()) return 1.0;                      // completely empty map
          const std::size_t k = static_cast<std::size_t>(0.995 * v.size());
          std::nth_element(v.begin(), v.begin() + k, v.end());
          return v[k];
      };

      auto tinyPositive = [](TH2* h)->double
      {
          double m = std::numeric_limits<double>::max();
          for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
              for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                  double z = h->GetBinContent(ix,iy);
                  if (z > 0.0 && z < m) m = z;
              }
          return (m == std::numeric_limits<double>::max()) ? 1.0 : m;
      };

      /* --- pick common colour‑scale limits -------------------------------- */
      const bool isHex = in.s->InheritsFrom(TH2Poly::Class());   // both histos same type

      double zMin, zMax;

      if (isHex) {                                   /* MBD hexagonal maps */
          /* full range: start at zero, top = 105 % of real per‑arm maximum  */
          zMin = 0.0;
          zMax = 1.05 * std::max(in.s->GetMaximum(), in.n->GetMaximum());
          if (zMax <= zMin) zMax = zMin + 1.0;       // safety for empty maps
      } else {                                       /* regular TH2 (sEPD etc.) */
          zMin = std::min(tinyPositive(in.s), tinyPositive(in.n));
          zMax = 1.05 * std::max(quantile99(in.s), quantile99(in.n));
      }
      in.s->SetMinimum(zMin);
      in.n->SetMinimum(zMin);      /* per‑pad upper limit will be set in drawPad() */
      
      /* ------------------------------------------------------------------ *
       *  finished S–N canvas                                               *
       * ------------------------------------------------------------------ */
      fs::path png = cPath(root, slice, DERIVED::subdir)
                     / (DERIVED::fileName(trig) + ".png");

      /* ensure the directory exists */
      try
      {
          ensure_dir(png.parent_path());
      }
      catch (const std::exception& e)
      {
          ulog::err("[NSDetectorQA] Cannot create output dir \"" +
                   png.parent_path().string() + "\": " + e.what());
          return false;
      }

      /* ------------------------------------------------------------------ *
       *  extract the 8‑digit run number from kInputFile                    *
       * ------------------------------------------------------------------ */
      int runNumber = 0;                           // fall‑back when no match
      {
          std::smatch m;
          std::regex_search(kInputFile, m, std::regex(R"((\d{8}))"));
          if (!m.empty()) runNumber = std::stoi(m[1].str());
      }

      /* ------------------------------------------------------------------ *
       *  all drawing inside its own try/catch                              *
       * ------------------------------------------------------------------ */
      try
      {
          TCanvas c("c_hit", "", 1200, 600);
          c.Divide(2, 1, 0.00, 0.01);

          auto drawPad = [&](TH2* h, const char* ttl, bool withZ)
          {
              /* ------------------------------------------------------------
               * Decide view style at run‑time:
               *   • TH2Poly  → Cartesian hexagons (MBD)
               *   • TH2      → φ–r polar map      (sEPD)
               * ------------------------------------------------------------ */
              const bool isHex = h->InheritsFrom(TH2Poly::Class());

              if (isHex)                         /* ---------- MBD hex view ---------- */
              {
                  gPad->SetRightMargin(0.18);
                  gPad->SetLeftMargin (0.12);
                  gPad->SetBottomMargin(0.12);
                  gPad->SetTopMargin  (0.08);

                  h->SetTitle("");
                  h->GetZaxis()->SetTitle("Counts");
                  h->GetZaxis()->SetTitleOffset(1.30);
                  h->Draw("POLZ");
              }
              else                                /* ---------- sEPD polar view ---------- */
              /* ---------- sEPD polar view (enhanced) -------------------------------- */
              {
                  //--------------------------------------------------------------------
                  // ❶  Geometry constants & diagnostics
                  //--------------------------------------------------------------------
                  const double rInner = h->GetYaxis()->GetXmin();           // 0.15 cm
                  const double rOuter = h->GetYaxis()->GetXmax();           // 3.51 cm
                  const double edge   = rOuter;                             // detector radius

                  const int    nPhi   = h->GetNbinsX();
                  const int    nRing  = h->GetNbinsY();
                  const double dPhi   = 2.0 * TMath::Pi() / nPhi;
                  const double dR     = (rOuter - rInner) / nRing;

                  ulog::trace("      sEPD polar view:");
                  ulog::trace("         nPhi  = " + std::to_string(nPhi)  +
                             "  (Δφ = " + std::to_string(TMath::RadToDeg() * dPhi) + "°)");
                  ulog::trace("         nRing = " + std::to_string(nRing) +
                             "  (ΔR = " + std::to_string(dR) + " cm)");
                  ulog::trace("         rInner = " + std::to_string(rInner) + " cm"
                             "  |  rOuter = " + std::to_string(rOuter) + " cm");

                  //--------------------------------------------------------------------
                  // ❷  Pad style
                  //--------------------------------------------------------------------
                  /* South pad (withZ == false) expands right‑ward,
                   * North pad (withZ == true ) expands left‑ward but still keeps 18 % for the palette.
                   */
                  if (withZ) {
                      gPad->SetLeftMargin (0.06);   // slimmer left edge – gains ~4 % width
                      gPad->SetRightMargin(0.18);   // palette room unchanged
                  } else {
                      gPad->SetLeftMargin (0.10);   // original
                      gPad->SetRightMargin(0.18);   // slimmer right edge – gains ~12 % width
                  }
                  gPad->SetBottomMargin(0.10);
                  gPad->SetTopMargin  (0.08);
                  gPad->SetFixedAspectRatio();
                  gPad->SetLogz(0);                 // linear colour scale

                  //--------------------------------------------------------------------
                  // ❸  First draw (heat map) – ensures palette is created
                  //--------------------------------------------------------------------
                  const char* optFirst = withZ ? "COLZ POL AH"      : "COL POL AH";
                  const char* optSame  = withZ ? "same COLZ POL AH" : "same COL POL AH";
                  
                  h->SetTitle("");
                  h->GetZaxis()->SetTitle("Counts");
                  h->GetZaxis()->SetTitleOffset(1.30);

                  /* --- per‑pad upper bound : 99.5‑percentile + 5 % head‑room -------- */
                  {
                      std::vector<double> z;  z.reserve(h->GetNbinsX()*h->GetNbinsY());
                      for (int ix=1; ix<=h->GetNbinsX(); ++ix)
                          for (int iy=1; iy<=h->GetNbinsY(); ++iy) {
                              double v = h->GetBinContent(ix,iy);
                              if (v>0.0) z.push_back(v);
                          }
                      if (!z.empty()) {
                          std::size_t k = static_cast<std::size_t>(0.995*z.size());
                          std::nth_element(z.begin(), z.begin()+k, z.end());
                          h->SetMaximum(1.05 * z[k]);
                      }
                  }

                  h->SetMinimum(0.0);               // colour scale anchored at zero
                  h->Draw(optFirst);

                  //--------------------------------------------------------------------
                  // ❹  Draw empty frame (square) so ROOT keeps axes, then redraw map
                  //--------------------------------------------------------------------
                  gPad->DrawFrame(-edge, -edge, edge, edge);
                  h->Draw(optSame);

                  //--------------------------------------------------------------------
                  // ❺  Draw φ‑sector spokes and radial circles (grid)
                  //--------------------------------------------------------------------
                  static std::vector<TLine> spokes;
                  static std::vector<TEllipse> rings;
                  if (spokes.empty())
                  {
                      // φ spokes
                      for (int i = 0; i < nPhi; ++i)
                      {
                          const double a = i * dPhi;
                          spokes.emplace_back(0., 0., edge * std::cos(a), edge * std::sin(a));
                          spokes.back().SetLineColor(kBlack);
                          spokes.back().SetLineWidth(1);
                      }
                      // radial rings (dashed)
                      for (int j = 1; j <= nRing; ++j)
                      {
                          const double r = rInner + j * dR;
                          rings.emplace_back(0., 0., r, r);        // store ring
                          rings.back().SetFillStyle(0);            // no fill – outline only
                          rings.back().SetLineColor(kBlack);       // black outline
                          rings.back().SetLineStyle(1);            // solid
                          rings.back().SetLineWidth(1);            // skinny
                      }
                  }
                  for (auto& l : spokes) {
                      l.SetLineColor(kBlack);
                      l.SetLineStyle(1);
                      l.SetLineWidth(1);
                      l.Draw();
                  }
                  for (auto& e : rings) {
                      e.SetFillStyle(0);
                      e.SetLineColor(kBlack);
                      e.SetLineStyle(1);
                      e.SetLineWidth(1);
                      e.Draw();
                  }

                  //--------------------------------------------------------------------
                  // ❻  Geometry & count tables  (printed to stdout)
                  //--------------------------------------------------------------------
                  std::cout << "\n┌───────────────────── Radial Rings ─────────────────────┐\n"
                            << "│   j   │ r_in [cm] │ r_out [cm] │   Σ Counts  │\n"
                            << "├───────┼───────────┼────────────┼──────────────┤\n";

                  for (int j = 1; j <= nRing; ++j)
                  {
                      const double rLo = rInner + (j-1)*dR;
                      const double rHi = rInner +  j   *dR;
                      double sum = 0.0;
                      for (int i = 1; i <= nPhi; ++i) sum += h->GetBinContent(i, j);
                      std::cout << std::setw(5) << j
                                << "   │ " << std::setw(9) << std::fixed << std::setprecision(3) << rLo
                                << " │ "  << std::setw(10) << rHi
                                << " │ "  << std::setw(11) << std::fixed << std::setprecision(0) << sum
                                << " │\n";
                  }
                  std::cout << "└─────────────────────────────────────────────────────────┘\n";

                  std::cout << "\n┌────────────────── φ Sectors (Δφ = "
                            << std::setw(6) << std::fixed << std::setprecision(2)
                            << TMath::RadToDeg() * dPhi << "°) ──────────────────┐\n"
                            << "│   i   │  φ_low [°] │  φ_high [°] │   Σ Counts │\n"
                            << "├───────┼────────────┼─────────────┼────────────┤\n";

                  for (int i = 1; i <= nPhi; ++i)
                  {
                      const double φLo = (i-1)*dPhi * TMath::RadToDeg();
                      const double φHi =  i   *dPhi * TMath::RadToDeg();
                      double sum = 0.0;
                      for (int j = 1; j <= nRing; ++j) sum += h->GetBinContent(i, j);
                      std::cout << std::setw(5) << i
                                << "   │ " << std::setw(10) << std::fixed << std::setprecision(2) << φLo
                                << " │ "  << std::setw(11) << φHi
                                << " │ "  << std::setw(10) << std::fixed << std::setprecision(0) << sum
                                << " │\n";
                  }
                  std::cout << "└─────────────────────────────────────────────────────────┘\n";

                  //--------------------------------------------------------------------
                  // ❼  Canvas / pad metrics (useful when debugging aspect issues)
                  //--------------------------------------------------------------------
                  if (auto* c = gPad->GetCanvas())
                  {
                      ulog::trace("         canvas  WxH = "
                                 + std::to_string(c->GetWw()) + " × "
                                 + std::to_string(c->GetWh()) + " px");
                      ulog::trace("         pad      WxH = "
                                 + std::to_string(gPad->GetWw()) + " × "
                                 + std::to_string(gPad->GetWh()) + " px");
                  }
              }

              const unsigned long long nEvt =
                  static_cast<unsigned long long>( h->GetEntries() );

              std::string trigLabel = prettifyTrigger(trig);   // ‘trig’ is in outer scope

              char buf[160];
              snprintf(buf, sizeof(buf),
                       "%s for %s (run %d, nEvents = %llu)",
                       ttl,
                       trigLabel.c_str(),
                       runNumber,
                       nEvt);

              TLatex lbl;
              lbl.SetNDC();
              lbl.SetTextAlign(13);         // left‑aligned, same as HCal/EMCal
              lbl.SetTextFont(42);
              lbl.SetTextSize(0.028);       // smaller font
              lbl.DrawLatex(0.08, 0.98, buf);
          };

          /* draw both pads ------------------------------------------------ */
          c.cd(1); drawPad(in.s, DERIVED::titleSouth, false);
          c.cd(2); drawPad(in.n, DERIVED::titleNorth, true );

          c.SaveAs(png.string().c_str());
          ulog::ok("[NSDetectorQA] Combined S/N map saved → " + png.string());
      }
      catch (const std::exception& e)
      {
          ulog::err("[NSDetectorQA] Error while drawing/saving \"" +
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
    using QA::QA;   // inherit all constructors

    // ────────────────────────────────────────────────────────────────
    // 1. Per‑histogram processing
    // ────────────────────────────────────────────────────────────────
    bool process(TObject* o) override
    {
        /* -----------------------------------------------------------
         *  Guard clause – we only handle TH1 histograms
         * --------------------------------------------------------- */
        if (!o->InheritsFrom(TH1::Class()))
            return false;

        const std::string n = o->GetName();
        const bool isVz   = (n.rfind("h_vertexZ_"  ,0) == 0);
        const bool isCent = (n.rfind("h_centrality_",0) == 0);
        if (!isVz && !isCent)
            return false;                       // irrelevant histogram

        /* “…/output/<RUN>/<trigger>/EventQA/…”  (always Inclusive) */
        const fs::path outPng = isVz
          ? root / "MBD"        / "zVertex"    / "VertexZ.png"
          : root / "centrality" / "Centrality.png";
        ensure_dir(outPng.parent_path());

        /* clone → local ownership → detach from any directory */
        std::unique_ptr<TH1> h(static_cast<TH1*>(o->Clone()));
        h->SetDirectory(nullptr);
        h->SetStats(0);

        const std::string runID          = root.parent_path().filename().string();   // "44792384"  or "Combined"
        const bool        isCombinedPass = (runID == "Combined");

        /* ─── verbose header when we are analysing the COMBINED file ─── */
        if (isCombinedPass) {
            std::cout << "\n╔══════════════════════════════════════════════════════════════╗\n"
                         "║  EventQA  –  COMBINED pass : new histogram received          ║\n"
                         "╠══════════════════════════════════════════════════════════════╣\n"
                         "║  Trigger     : " << std::setw(20) << root.filename().string()          << " ║\n"
                         "║  Histogram   : " << std::setw(40) << n                                   << " ║\n"
                         "║  Type        : " << std::setw(12) << (isVz?  "vertex‑Z" : "centrality") << " ║\n"
                         "║  Entries     : " << std::setw(12) << h->GetEntries()                    << " ║\n"
                         "╚══════════════════════════════════════════════════════════════╝\n";
        }

        /* quick‑reject empty histograms (very low stats)  ---------- */
        if (h->GetEntries() < 20)
        {
            std::cout << "[EventQA‑DBG] skip   run=" << runID
                      << "  hist='" << n << "'  entries=" << h->GetEntries() << '\n';
            return false;
        }

        // ============================================================
        // (A)  primary‑vertex z  –  robust iterative Gaussian fit
        // ============================================================
        if (isVz)
        {
            const std::string trig = root.filename().string();   // trigger folder name
            std::cout << "\n[EventQA‑DBG] >>> vertex‑Z processing start <<<\n"
                      << "   runID          : " << runID           << '\n'
                      << "   trigger        : " << trig            << '\n'
                      << "   hist name      : " << n               << '\n'
                      << "   entries        : " << h->GetEntries() << '\n'
                      << "   output PNG     : " << outPng          << '\n';

            /* ------------------------------------------------------- *
             * STEP‑0 : robust seed from 16‑50‑84 % quantiles          *
             * ------------------------------------------------------- */
            double probs[3] = {0.16, 0.50, 0.84};
            double q[3];
            h->GetQuantiles(3, q, probs);
            double mu    = q[1];
            double sigma = 0.5*(q[2]-q[0]);          // 68 % inter‑quantile width
            if (sigma <= 0) sigma = h->GetRMS();
            if (sigma <= 0) sigma = 1.0;             // hard fallback
            std::cout << "[EventQA‑DBG] seed μ=" << mu << "  σ=" << sigma << '\n';

            TF1 g("g","gaus", mu-3*sigma, mu+3*sigma);
            g.SetLineColor(kRed+1);
            g.SetLineWidth(2);
            g.SetParameters(h->GetMaximum(), mu, sigma);

            /* ------------------------------------------------------- *
             * STEP‑1 : iterative 2.5 σ shrink until convergence       *
             * ------------------------------------------------------- */
            constexpr int    kMaxIter   = 5;
            constexpr double kNSigmaFit = 2.5;
            constexpr double kTol       = 1e-3;

            bool  fitOK    = false;
            int   fitStat  = -1;
            double chi2    = 0.0;

            for (int it = 0; it < kMaxIter; ++it)
            {
                const double lo = mu - kNSigmaFit*sigma;
                const double hi = mu + kNSigmaFit*sigma;
                g.SetRange(lo, hi);
                g.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sigma);

                TFitResultPtr res = h->Fit(&g,"Q0RSLL");
                fitStat           = static_cast<int>(res);
                fitOK             = (fitStat == 0);

                if (!fitOK)
                {
                    std::cout << "[EventQA‑DBG] iter " << it
                              << "  –  Fit **FAILED** (status " << fitStat << ")\n";
                    break;
                }

                const double muNew    = g.GetParameter(1);
                const double sigmaNew = std::fabs(g.GetParameter(2));

                const bool   conv =
                       std::fabs(muNew   - mu)    < kTol*sigma &&
                       std::fabs(sigmaNew- sigma) < kTol*sigma;

                mu    = muNew;
                sigma = sigmaNew;
                chi2  = g.GetChisquare();

                std::cout << "[EventQA‑DBG] iter " << it
                          << "  μ="  << mu
                          << "  σ="  << sigma
                          << "  χ²=" << chi2
                          << "  status=" << fitStat
                          << (conv ? "  (converged)" : "") << '\n';

                if (conv) break;
            }

            /* sanity filter – reject obviously wrong fits            */
            const bool saneParams =
                   std::fabs(mu) < 50.0   &&  sigma > 0.5 && sigma < 30.0;
            if (!saneParams)
            {
                std::cout << "[EventQA‑DBG] sanity check **FAILED**  |μ|<50 cm  && 0.5<σ<30 cm\n";
                fitOK = false;
            }

            const double muErr  = fitOK ? g.GetParError(1) : 0.0;
            const double sigErr = fitOK ? g.GetParError(2) : 0.0;

            /* ------------------------------------------------------- *
             * Per‑run diagnostic PNG                                  *
             * ------------------------------------------------------- */
            {
                TCanvas c("c_vz","", 900, 600);
                h->Draw();
                if (fitOK) g.Draw("SAME");

                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
                tx.DrawLatex(0.15,0.88,Form("#mu = %.2f #pm %.2f cm", mu,  muErr));
                tx.DrawLatex(0.15,0.84,Form("#sigma = %.2f #pm %.2f cm", sigma, sigErr));
                c.SaveAs(outPng.string().c_str());
            }

            /* ------------------------------------------------------- *
             * Decide whether to cache this fit for the summary        *
             * ------------------------------------------------------- */
            const long long nEvt  = static_cast<long long>(h->GetEntries());
            auto& maxEvtForRun    = s_evtCounts[runID];
            const bool firstForRun = (s_points.find(runID) == s_points.end());
            const bool takeThis    = fitOK && (firstForRun || nEvt > maxEvtForRun);

            std::cout << "[EventQA‑DBG] decision: "
                      << (takeThis ? "STORE" : "skip")
                      << "   (fitOK="  << fitOK
                      << ", first="   << firstForRun
                      << ", nEvt="    << nEvt
                      << ", maxEvt="  << maxEvtForRun
                      << ", combinedPass=" << isCombinedPass << ")\n";

            /* ---------- store or update counters ------------------ */
            if (!isCombinedPass && takeThis)
            {
                VzPoint& p = s_points[runID];
                p.mu        = mu;        p.muErr    = muErr;
                p.sigma     = sigma;     p.sigmaErr = sigErr;
                p.nEvt      = nEvt;

                p.hist.reset(static_cast<TH1*>(h->Clone()));
                p.hist->SetDirectory(nullptr);
                p.hist->Scale(1.0 / p.hist->GetEntries());   // normalise overlay

                maxEvtForRun = nEvt;      // update AFTER successful storage
                std::cout << "[EventQA‑DBG] stored as current BEST for run "
                          << runID << '\n';
            }
            else if (!isCombinedPass && nEvt > maxEvtForRun)
            {
                maxEvtForRun = nEvt;      // keep counter so later good fit can compare
                std::cout << "[EventQA‑DBG] updated event counter only → "
                          << maxEvtForRun << '\n';
            }

            std::cout << "[EventQA‑DBG] <<< vertex‑Z processing end <<<\n";
        }

        // ============================================================
        // (B)  centrality spectrum – plain plot, normalised
        // ============================================================
        else   /* isCent */
        {
            std::cout << "[EventQA‑DBG] centrality hist  run=" << runID
                      << "  entries=" << h->GetEntries() << '\n';

            TCanvas c("c_cent","",900,600);
            h->Draw();

            TLatex tx; tx.SetNDC(); tx.SetTextAlign(31); tx.SetTextSize(0.04);
            tx.DrawLatex(0.97,0.94, runID.c_str());
            c.SaveAs(outPng.string().c_str());

            if (!isCombinedPass)
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
        /* we are inside “…/output/Combined/<trigger>”                */
        const bool combinedHere = (root.parent_path().filename() == "Combined");
        if (!combinedHere || s_summaryWritten) return;
        s_summaryWritten = true;

        std::cout << "\n╔══════════════════════════════════════════════════════════════╗\n"
                     "║                EventQA  –  COMBINED SUMMARY                 ║\n"
                     "╠═══════════════════════╤═════════════════════════════════════╣\n"
                     "║  Vertex‑Z fits        │ " << std::setw(8) << s_points.size()    << "                                 ║\n"
                     "║  Centrality histograms│ " << std::setw(8) << s_centHists.size() << "                                 ║\n"
                     "╚═══════════════════════╧═════════════════════════════════════╝\n";
        std::cout << std::left;
        std::cout << "\n┌────────┬───────────┬──────────┬───────────┐\n"
                     "│  Run   │   μ (cm)  │ σ (cm)   │  Events   │\n"
                     "├────────┼───────────┼──────────┼───────────┤\n";
        for (const auto& [run,p] : s_points)
            std::cout << "│ " << std::setw(6)  << run
                      << " │ " << std::setw(9) << std::fixed << std::setprecision(2) << p.mu
                      << " │ " << std::setw(9) << p.sigma
                      << " │ " << std::setw(9) << p.nEvt << " │\n";
        std::cout << "└────────┴───────────┴──────────┴───────────┘\n\n";


        /* print table of stored run points for extra visibility      */
        for (auto& [run,p] : s_points)
            std::cout << "    run " << run << "  μ=" << p.mu
                      << "  σ=" << p.sigma << "  nEvt=" << p.nEvt << '\n';

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
            std::cout << "  – writing VertexZ_AllRuns overlay\n";
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
            std::cout << "    saved → " << vzOverlay << '\n';
        }

        /* (B) centrality overlay + paginated small‑multiples -------- */
        if (s_centHists.size() > 1)
        {
            /* ─── B‑1 : classic colour overlay (unchanged logic) ─── */
            std::cout << "  – writing Centrality_AllRuns overlay\n";
            TCanvas cOv("c_cent_overlay","Centrality – all runs",900,600);
            TLegend leg(0.50,0.45,0.75,0.65); leg.SetBorderSize(0);

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
            cOv.SaveAs(centOverlay.string().c_str());
            std::cout << "    saved → " << centOverlay << '\n';

            /* ─── B‑2 : 10×10 grid pages with one run per pad ───── */
            const int nPerPage = 100;              // 10 × 10
            std::vector<std::pair<std::string,TH1*>> ordered;
            ordered.reserve(s_centHists.size());

            /* reproducible order by numeric run number (if any) ---- */
            for (auto& [run,h] : s_centHists) ordered.emplace_back(run, h.get());
            std::sort(ordered.begin(), ordered.end(),
                      [](auto& a, auto& b)
                      {
                          const bool aNum = std::all_of(a.first.begin(),a.first.end(),::isdigit);
                          const bool bNum = std::all_of(b.first.begin(),b.first.end(),::isdigit);
                          if (aNum && bNum) return std::stoi(a.first) < std::stoi(b.first);
                          return a.first < b.first;
                      });

            const int nPages = static_cast<int>((ordered.size() + nPerPage - 1) / nPerPage);
            const fs::path outDir = root / "centrality";
            ensure_dir(outDir);

            for (int pg = 0; pg < nPages; ++pg)
            {
                const int start = pg * nPerPage;
                const int stop  = std::min<int>(start + nPerPage, ordered.size());
                const int nHere = stop - start;

                TCanvas cGrid( Form("c_cent_grid_%02d", pg+1),
                               "Centrality – per run", 2500, 2500 );
                cGrid.Divide(10,10, 0.001, 0.001);

                for (int idx = start; idx < stop; ++idx)
                {
                    cGrid.cd(idx - start + 1);
                    gPad->SetLeftMargin  (0.10);
                    gPad->SetRightMargin (0.02);
                    gPad->SetTopMargin   (0.05);
                    gPad->SetBottomMargin(0.10);

                    TH1* h = ordered[idx].second;
                    h->SetLineColor(kBlue+1); h->SetLineWidth(1);
                    h->Draw("HIST");

                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.05);
                    tx.SetTextAlign(33);    // top‑right
                    tx.DrawLatex(0.94, 0.90, ordered[idx].first.c_str());
                }

                std::string fname = Form("Centrality_AllRuns_grid_p%02d.png", pg+1);
                fs::path pagePng  = outDir / fname;
                cGrid.SaveAs(pagePng.string().c_str());
                std::cout << "    saved → " << pagePng << '\n';
            }
        }
        
        /* (C)  μ,σ  versus run number ------------------------------ */
        if (s_points.size() > 1)
        {
            std::cout << "  – writing μ,σ vs run summary plot\n";
            std::vector<int> runs;
            for (auto& [r,_] : s_points)
                if (std::all_of(r.begin(),r.end(),::isdigit))
                    runs.push_back(std::stoi(r));
            std::sort(runs.begin(), runs.end());

            const int n = runs.size();
            std::vector<double> x(n), yMu(n), eMu(n), ySi(n), eSi(n);
            for (int i = 0; i < n; ++i)
            {
                std::string runStr = Form("%08d", runs[i]);   // keep the 8‑digit run key
                const auto& p      = s_points.at(runStr);     // fetch existing entry
                x[i]   = runs[i];
                yMu[i] = p.mu;        eMu[i] = p.muErr;
                ySi[i] = p.sigma;     eSi[i] = p.sigmaErr;
            }

            auto gMu = std::make_unique<TGraphErrors>(n,x.data(),yMu.data(),
                                                      nullptr,eMu.data());
            auto gSi = std::make_unique<TGraphErrors>(n,x.data(),ySi.data(),
                                                      nullptr,eSi.data());
            gMu->SetMarkerStyle(kFullCircle); gMu->SetLineWidth(2);
            gSi->SetMarkerStyle(kOpenCircle); gSi->SetLineWidth(2);

            /* dynamic Y‑ranges ------------------------------------ */
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

            gMu->SetMinimum(-1.10 * absMu); gMu->SetMaximum(1.10 * absMu);
            gSi->SetMinimum(0.0);           gSi->SetMaximum(1.10 * maxSig);

            TCanvas c("c_mu_sigma_vs_run","vertex‑Z  #mu,#sigma  vs run",900,800);

            TPad* p1 = new TPad("p1","",0,0.35,1,1);
            p1->SetBottomMargin(0.02); p1->Draw(); p1->cd();
            /* ───── upper pad : mean μ ───────────────────────────────────────── */
            gMu->SetTitle("; ;#mu  [cm]");          // empty X‑title → only on lower pad
            gMu->GetXaxis()->SetLabelSize(0);       // hide run‑number labels here
            gMu->GetXaxis()->SetTickLength(0.02);   // keep short ticks
            gMu->GetYaxis()->SetTitleSize(0.05);
            gMu->GetYaxis()->SetLabelSize(0.04);
            gMu->Draw("AP");

            /* ───── lower pad : σ ───────────────────────────────────────────── */
            c.cd();
            TPad* p2 = new TPad("p2","",0,0,1,0.32);
            p2->SetTopMargin(0.02);
            p2->SetBottomMargin(0.35);              // more room for larger X‑axis text
            p2->Draw(); p2->cd();

            gSi->SetTitle(";Run number;#sigma  [cm]");
            gSi->GetXaxis()->SetTitleSize(0.06);    // larger labels / title
            gSi->GetXaxis()->SetLabelSize(0.05);
            gSi->GetYaxis()->SetTitleSize(0.05);    // match top panel size
            gSi->GetYaxis()->SetLabelSize(0.04);
            gSi->Draw("AP");

            fs::path pngRun = root / "MBD" / "zVertex" / "VertexZ_MeanSigma_vs_Run.png";
            ensure_dir(pngRun.parent_path());
            c.SaveAs(pngRun.string().c_str());
            std::cout << "    saved → " << pngRun << '\n';
        }

        std::cout << "[EventQA‑DBG] ===== summary complete =====\n";
    }

    // ────────────────────────────────────────────────────────────────
    // 3. Static containers
    // ────────────────────────────────────────────────────────────────
    struct VzPoint {
        double mu{}, muErr{}, sigma{}, sigmaErr{};
        long long nEvt{};
        std::unique_ptr<TH1> hist;
    };

    static inline std::unordered_map<std::string,VzPoint>              s_points;
    static inline std::unordered_map<std::string,std::unique_ptr<TH1>> s_centHists;
    static inline std::unordered_map<std::string,long long>            s_evtCounts;

    static inline bool  s_summaryWritten = false;

    /* accessor for external code if needed */
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

  enum class Lvl { DBG, INF, WRN, ERR };

  static void log(Lvl level, const std::string& msg)
  {
        switch (level) {
            case Lvl::DBG: ulog::trace(msg); break;
            case Lvl::INF: ulog::info (msg); break;
            case Lvl::WRN: ulog::warn (msg); break;
            case Lvl::ERR: ulog::err  (msg); break;
        }
  }

  // =================================================================
  // 1. Per‑histogram processing (called many times per file)
  // =================================================================
  bool process(TObject* o) override
  {
      //------------------------------------------------------------------
      // 0. Preliminary sanity checks
      //------------------------------------------------------------------
      if (!o) {
        log(Lvl::WRN,"process(): received nullptr – skipped");
        return false;
      }
      if (!o->InheritsFrom(TH1::Class())) {          // we process only TH1/2/3
        log(Lvl::DBG,Form("process(): \"%s\" is not TH1 – skipped",o->GetName()));
        return false;
      }

      const std::string n = o->GetName();            // histogram name (immutable)

      // Categorise once – fast string‑prefix tests (rfind(...,0) == 0 → starts with)
      const bool is1D = n.rfind("h_maxJetEt_"         ,0) == 0;
      const bool is2D = n.rfind("h_leadEt_vs_subEt_"  ,0) == 0;
      const bool is3D =   n.rfind("h_jetEt_area_nConst_",0) == 0
                       || n.rfind("h_jetEt_eta_phi_"     ,0) == 0     // legacy name
                       || n.rfind("h_jetEtEtaPhi_"       ,0) == 0;    // NEW name

      if (!is1D && !is2D && !is3D) {                 // not a Jet‑QA histogram
        log(Lvl::DBG,"process(): \""+n+
                      "\" not in jet‑QA family – skipped");
        return false;
      }

      log(Lvl::DBG,"process(): start \""+n+"\"  "
                   "(is1D="+std::to_string(is1D)+
                   ", is2D="+std::to_string(is2D)+
                   ", is3D="+std::to_string(is3D)+")");

      //------------------------------------------------------------------
      // 1. Everything inside a try/catch so that one bad histogram never
      //    aborts the entire QA pass.
      //------------------------------------------------------------------
      try
      {
        //----------------------------------------------------------------
        // 1.1  Build destination directory:  …/<slice>/jetQA/generalHistos/<rXX>
        //----------------------------------------------------------------
        const std::string slice = sliceKey(n);       // e.g. "Inclusive" or "Cent_0_10"
        const std::string rLab  = radiusTag(n);      // e.g. "r02", "r04", …

        fs::path baseGen =
            cPath(root, slice,
                  fs::path("jetQA/generalHistos") / rLab);

        ensure_dir(baseGen);                         // may throw on permission issues

        //----------------------------------------------------------------
        // 1.2  Delegate to the original helpers; propagate their return
        //      status so upstream code stays unchanged.
        //----------------------------------------------------------------
        bool ok = false;
        if (is1D)
          {
            TH1* h1 = static_cast<TH1*>(o);
            ok = save1Dplot(h1, baseGen, n);

            /* cache Inclusive h_maxJetEt_* for later summary */
            if (n.rfind("h_maxJetEt_",0) == 0 && slice == "Inclusive")
            {
                auto hClone = std::unique_ptr<TH1>(static_cast<TH1*>(h1->Clone()));
                hClone->SetDirectory(nullptr);
                _cachedMaxEt[slice] = std::move(hClone);
            }
        }
        else if (is2D)
          ok = handle2D  (static_cast<TH2*>(o), baseGen, n);
        else            // is3D
          ok = handle3D  (static_cast<TH3*>(o), baseGen, n);

          if (!ok)
            log(Lvl::WRN,"process(): helper reported failure for \""+n+"\"");

          /* ────────────────────────────────────────────────────────────────
           * NEW: keep one clone per centrality bin so we can build a 2×3
           *      overview grid at the very end of the run.
           * ──────────────────────────────────────────────────────────────── */
          if (slice.rfind("Cent_",0) == 0)                 /* skip “Inclusive” */
          {
              /* strip the “_lo_hi_” part so all six slices map to one key     */
              std::string base =
                  std::regex_replace(n, std::regex(R"(_\d{1,3}_\d{1,3}_)"), "_");

              auto cl = std::shared_ptr<TH1>(static_cast<TH1*>(o->Clone()));
              cl->SetDirectory(nullptr);
              _centCache[base].emplace_back(slice, std::move(cl));
          }

          log(Lvl::DBG,"process(): finished \""+n+"\"  (status=" +
                       std::string(ok ? "OK" : "FAIL") + ')');
          return ok;

      }
      catch (const std::exception& ex)
      {
        log(Lvl::ERR,std::string("process(): C++ exception for \"")+n+
                      "\" – "+ex.what());
        return false;                                 // safe fallback
      }
      catch (...)
      {
        log(Lvl::ERR,std::string("process(): UNKNOWN exception for \"")+n+"\"");
        return false;
      }
  }


  // =================================================================
  // 2. Final summary – executed once at the very end of the run
  // =================================================================
  ~JetQA() override
  {
      //----------------------------------------------------------------
      // 0. Early exit guards
      //----------------------------------------------------------------
      try {
        constexpr double kEtMin = 20.0;            // GeV threshold for the yield
        std::map<std::string,std::unique_ptr<TH1>> maxEtPerSlice;

        //----------------------------------------------------------------
        // 1. Harvest the cached Inclusive h_maxJetEt_* histograms
        //----------------------------------------------------------------
        if (_cachedMaxEt.empty()) {
            log(Lvl::WRN,"[JetQA::~JetQA] no cached h_maxJetEt_* histograms – summary skipped");
            return;
        }

        for (auto& [slice,hPtr] : _cachedMaxEt)
            if (hPtr && hPtr->Integral() > 0)
                maxEtPerSlice[slice].reset(static_cast<TH1*>(hPtr->Clone()));


        //----------------------------------------------------------------
        // 2. Prepare output directory
        //----------------------------------------------------------------
        fs::path dirSum = root / "jetQA/summary";
        try { ensure_dir(dirSum); }
        catch (const std::exception& ex) {
          log(Lvl::ERR,std::string("[JetQA::~JetQA] ensure_dir(")+dirSum.string()+
                        ") failed – "+ex.what());
          return;
        }

          //----------------------------------------------------------------
          // 3. Yield(E_T)  –  per‑slice PNGs **and** combined overlay
          //----------------------------------------------------------------
          TCanvas cYield("c_yieldEt","Jet yield vs E_{T}",1100,850);
          cYield.SetLogy();

          /* legend positioned high‑left, large markers for overlay */
          TLegend leg(0.18,0.68,0.48,0.88); leg.SetBorderSize(0);
          const int colList[]{kRed+1,kBlue+2,kGreen+2,kMagenta+2,kOrange+1};
          int iCol = 0;

          std::vector<double> xCent, yYield;                // for cent‑trend

          for (auto& [slice, h] : maxEtPerSlice)
          {
              /* --- 3.1 centrality coordinate -------------------------------- */
              double xC = 50.0;
              if (slice.rfind("Cent_",0) == 0) {
                  std::smatch m;
                  if (std::regex_match(slice,m,std::regex(R"(Cent_([0-9]+)_([0-9]+))")))
                      xC = 0.5*(std::stod(m[1])+std::stod(m[2]));
              }

              /* --- 3.2 integrate yield above kEtMin ------------------------- */
              const int    binMin = h->FindBin(kEtMin);
              double       yield  = h->Integral(binMin, h->GetNbinsX());
              const long long nEv = static_cast<long long>(h->GetEntries());
              if (nEv==0) continue;
              yield /= nEv;

              xCent.push_back(xC);  yYield.push_back(yield);

              /* --- 3.3 differential spectrum (clone, norm‑per‑event) ------- */
              auto hDiff = std::unique_ptr<TH1>(static_cast<TH1*>(h->Clone()));
              hDiff->Scale(1.0/nEv,"width");
              const int col = colList[iCol % (sizeof(colList)/sizeof(int))];
              hDiff->SetLineColor(col);  hDiff->SetLineWidth(2);
              hDiff->GetYaxis()->SetTitle("1/N_{ev}  dN/dE_{T}  [GeV^{-1}]");

              /* ---------- A) individual PNG for this slice ---------- */
              {
                  TCanvas cSlice(Form("c_yield_%s",slice.c_str()),
                                 Form("Jet yield vs E_{T} – %s",slice.c_str()),
                                 1100,850);
                  cSlice.SetLogy();
                  hDiff->Draw("HIST");
                  const fs::path pngSlice =
                         dirSum/(std::string("JetYield_vs_Et_")+slice+".png");
                  cSlice.SaveAs(pngSlice.string().c_str());
              }

              /* ---------- B) add to combined overlay ---------------- */
              hDiff->SetTitle("");                 // suppress title in overlay
              hDiff->Draw(iCol==0 ? "HIST" : "HIST SAME");
              leg.AddEntry(hDiff.get(), slice.c_str(),"l");

              _owned1D.push_back(std::move(hDiff));
              ++iCol;
        }

        /* draw legend & save overlay only if ≥1 slice was plotted */
        if (iCol>0) {
              leg.Draw();
              const fs::path pngYield = dirSum/"JetYield_vs_Et_AllSlices.png";
              cYield.SaveAs(pngYield.string().c_str());
        }


        //----------------------------------------------------------------
        // 4. Yield vs centrality (if ≥2 points)
        //----------------------------------------------------------------
        if (xCent.size() > 1)
        {
          auto gCent = std::make_unique<TGraphErrors>(xCent.size());
          for (std::size_t i=0;i<xCent.size();++i)
            gCent->SetPoint(i, xCent[i], yYield[i]);

          TCanvas cCent("c_yieldCent","Jet yield vs centrality",900,700);
          gCent->SetTitle(Form("Jet yield  E_{T}>%.0f GeV",kEtMin));
          gCent->GetXaxis()->SetTitle("centrality [%]");
          gCent->GetYaxis()->SetTitle("jets / event");
          gCent->SetMarkerStyle(kFullCircle);  gCent->SetLineWidth(2);
          gCent->Draw("AP");

          const fs::path pngCent = dirSum/"JetYield_vs_Centrality.png";
          try { cCent.SaveAs(pngCent.string().c_str()); }
          catch (const std::exception& ex) {
            log(Lvl::ERR,std::string("[JetQA::~JetQA] failed to save ")+pngCent.string()+
                          " – "+ex.what());
          }
          _ownedGraphs.push_back(std::move(gCent));
        }

        writeCentralityGrids();          /* NEW – 2×3 overview per histogram */

        log(Lvl::INF,"[JetQA::~JetQA] summary successfully written ("+
                                std::to_string(maxEtPerSlice.size())+
                                " slice(s))  and centrality grids");
      }
      catch (const std::exception& ex)
      {
        log(Lvl::ERR,std::string("[JetQA::~JetQA] C++ exception – ")+ex.what());
      }
      catch (...)
      {
        log(Lvl::ERR,"[JetQA::~JetQA] UNKNOWN exception");
      }
    }

 private:
    // ------------------------------------------------------------------
    //  radiusTag  – extract jet‑radius label (r02, r04 …) from the name
    // ------------------------------------------------------------------
    static std::string radiusTag(const std::string& hname)
    {
      ulog::trace("[JetQA] radiusTag()  →  hname=\"" + hname + '"');
      try {
        std::smatch m;  std::regex re(R"(_(r[0-9]+|R[0-9]+)_)");
        if (std::regex_search(hname, m, re)) {
          ulog::trace("[JetQA] radiusTag()  ←  \"" + m[1].str() + '"');
          return m[1].str();
        }
        ulog::warn("[JetQA] radiusTag()  –  pattern not found in \"" + hname + '"');
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] radiusTag() exception – ")+ex.what());
      }
      return "UnknownR";
    }

    // ------------------------------------------------------------------
    //  save1Dplot  – render 1‑D histogram with log‑Y
    // ------------------------------------------------------------------
    bool save1Dplot(TH1* h, const fs::path& dir, const std::string& hname)
    {
      ulog::trace("[JetQA] save1Dplot() enter – " + hname);
      if (!h) { ulog::err("[JetQA] save1Dplot() received nullptr"); return false; }

      try {
        ensure_dir(dir);
        TCanvas c; c.SetLogy();
        h->SetStats(0);
        h->SetTitle("");                 // avoid ROOT’s default title box
        h->Draw();

        TLatex ttl;
        ttl.SetNDC();                   // normalised device coordinates
        ttl.SetTextFont(42);            // plain, professional font
        ttl.SetTextAlign(23);           // centered horizontally, top‑aligned
        ttl.SetTextSize(0.045);
        ttl.DrawLatex(0.50, 0.94, makeTitle(hname).c_str());

        const fs::path png = dir/(hname + ".png");
        c.SaveAs(png.string().c_str());
        ulog::info("[JetQA] 1‑D plot saved → " + png.string());
        return true;
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] save1Dplot() exception – ")+ex.what());
        return false;
      }
    }

    // =============== 2‑D ==================================================
    bool handle2D(TH2* h, const fs::path& dir, const std::string& hname)
    {
      ulog::trace("[JetQA] handle2D() enter – " + hname);
      if (!h) { ulog::err("[JetQA] handle2D() nullptr"); return false; }

      try {
        ensure_dir(dir);
        TCanvas c;
        h->SetStats(0);
        h->SetTitle("");

        /* --- auto‑range: first & last non‑empty bins ------------------- */
        int fx = h->GetNbinsX()+1;                  // first   non‑empty X bin
        int fy = h->GetNbinsY()+1, ly = 0;          // full Y auto‑range
        double zMin = std::numeric_limits<double>::max();

        for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
              for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                  const double c = h->GetBinContent(ix,iy);
                  if (c <= 0) continue;
                  if (ix < fx) fx = ix;
                  if (iy < fy) fy = iy;  if (iy > ly) ly = iy;
                  if (c  < zMin) zMin = c;
              }

        /* X‑axis: keep full width on the right (lx = NbinsX) */
        if (fx <= h->GetNbinsX()) h->GetXaxis()->SetRange(fx, h->GetNbinsX());

        /* Y‑axis: still trim empty margins top & bottom */
        if (fy < ly) h->GetYaxis()->SetRange(fy, ly);

        /* Z‑axis minimum (colour scale) */
        if (zMin < std::numeric_limits<double>::max()) h->SetMinimum(zMin);

        /* ---------------------------------------------------------------
         *  Special treatment for   h_leadEt_vs_subEt_*   histograms
         *     • enable log‑Z colour scale for better dynamic range
         *     • suppress the diagonal guideline
         * -------------------------------------------------------------- */
        const bool isLeadSub = (hname.rfind("h_leadEt_vs_subEt_",0) == 0);
        if (isLeadSub) c.SetLogz();

        h->Draw("COLZ");

        /* show guideline only for non‑lead/sub histograms */
        if (!isLeadSub) {
              const double xmax = h->GetXaxis()->GetXmax();
              TLine diag(0,0, xmax, xmax);
              diag.SetLineStyle(2);
              diag.SetLineWidth(2);
              diag.Draw();
        }

        TLatex ttl; ttl.SetNDC();
        ttl.SetTextFont(42);
        ttl.SetTextAlign(23);
        ttl.SetTextSize(0.03);
        ttl.DrawLatex(0.50, 0.94, makeTitle(hname).c_str());


        const fs::path png = dir/(hname + ".png");
        c.SaveAs(png.string().c_str());
        ulog::info("[JetQA] 2‑D plot saved → " + png.string());
        return true;
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] handle2D() exception – ")+ex.what());
        return false;
      }
    }
    
    // ------------------------------------------------------------------
    //  saveJetQA2D – draw 2‑D histogram keeping the descriptive header
    // ------------------------------------------------------------------
    static void saveJetQA2D(TH2* h, const fs::path& png)
    {
        if (!h) { ulog::err("[JetQA] saveJetQA2D() nullptr"); return; }

        try {
            ensure_dir(png.parent_path());

            /* ── canvas identical to other Jet‑QA helpers ───────────────── */
            TCanvas c;
            h->SetStats(0);

            /* preserve the header that saveProjection() injected */
            const std::string header = h->GetTitle();
            h->SetTitle("");                     // suppress ROOT’s own title box

            /* automatic axis‑range trimming (same logic as handle2D) */
            int fx = h->GetNbinsX()+1, lx = 0;
            int fy = h->GetNbinsY()+1, ly = 0;
            double zMin = std::numeric_limits<double>::max();
            for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
                for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
                    const double cBin = h->GetBinContent(ix,iy);
                    if (cBin <= 0) continue;
                    if (ix < fx) fx = ix;  if (ix > lx) lx = ix;
                    if (iy < fy) fy = iy;  if (iy > ly) ly = iy;
                    if (cBin < zMin) zMin = cBin;
                }
            if (fx < lx) h->GetXaxis()->SetRange(fx,lx);
            if (fy < ly) h->GetYaxis()->SetRange(fy,ly);
            if (zMin < std::numeric_limits<double>::max()) h->SetMinimum(zMin);

            h->Draw("COLZ");

            /* draw the preserved header */
            TLatex ttl; ttl.SetNDC();
            ttl.SetTextFont(42);
            ttl.SetTextAlign(23);
            ttl.SetTextSize(0.03);
            ttl.DrawLatex(0.50, 0.94, header.c_str());

            c.SaveAs(png.string().c_str());
        }
        catch (const std::exception& ex) {
            ulog::err(std::string("[JetQA] saveJetQA2D() exception – ")+ex.what());
        }
    }

    // =============== 3‑D ==================================================
    bool handle3D(TH3* h3, const fs::path& dir, const std::string& hname)
    {
      ulog::trace("[JetQA] handle3D() enter – " + hname);
      if (!h3) { ulog::err("[JetQA] handle3D() nullptr"); return false; }

      try {
        ensure_dir(dir);
        h3->SetTitle(makeTitle(hname).c_str());

        // main 3‑D view
        save3D(h3, dir/(hname + "_3D.png"));

        const bool isEtaPhi =
                  hname.find("_eta_phi_") != std::string::npos   ||   // legacy
                  hname.find("EtaPhi_")    != std::string::npos;      // NEW

        if (!isEtaPhi) {
                saveProjection(h3,"yx", dir/(hname + "_Et_vs_Area.png"));    // E_T vs A
                saveProjection(h3,"xz", dir/(hname + "_Et_vs_Nconst.png"));  // E_T vs N
                saveProjection(h3,"yz", dir/(hname + "_Area_vs_Nconst.png"));/* A vs N */
        } else {
                /* Jet E_T × η × φ  →  two useful views */
                saveProjection(h3,"zx", dir/(hname + "_Et_vs_Eta.png"));     // E_T vs η
                saveProjection(h3,"zy", dir/(hname + "_Et_vs_Phi.png"));     // E_T vs φ
        }
        return true;
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] handle3D() exception – ")+ex.what());
        return false;
      }
    }

    // ------------------------------------------------------------------
    //  save3D – single 3‑D view (semi‑transparent boxes)
    // ------------------------------------------------------------------
    static void save3D(TH3* h, const fs::path& png)
    {
      if (!h) { ulog::err("[JetQA] save3D() nullptr"); return; }
      ulog::trace("[JetQA] save3D() → " + png.string());

      try {
        ensure_dir(png.parent_path());

        /* ── canvas with generous margins for colour bar and axis titles ── */
        TCanvas c("c3D","",1200,1000);
        c.SetLeftMargin  (0.14);
        c.SetBottomMargin(0.14);
        c.SetRightMargin (0.32);                 // extra room for palette

        /* ── histogram cosmetics ─────────────────────────────────────────── */
        h->SetStats(0);
        h->SetContour(99);

        /* shift axis titles away from tick labels so they do not overlap */
        h->GetXaxis()->SetTitleOffset(1.6);
        h->GetYaxis()->SetTitleOffset(2.0);
        h->GetZaxis()->SetTitleOffset(1.4);

        h->SetTitle("");          // suppress ROOT’s default title

        const bool isEtaPhi = (strstr(h->GetName(),"EtEtaPhi_")!=nullptr) ||
                                 (strstr(h->GetName(),"_eta_phi_")!=nullptr);

          /* For η–φ maps we first collapse the E_T axis so that the bar
             height represents the number of jets in each (η,φ) bin.          */
        std::unique_ptr<TH2> h2;          // only used for η–φ maps
        TH1* hDraw = h;                   // histogram actually drawn (TH2/TH3 inherit from TH1)

        if (isEtaPhi) {
              h2.reset( static_cast<TH2*>( h->Project3D("yx") ) );  // η vs φ
              h2->SetDirectory(nullptr);
              h2->SetStats(0);
              h2->SetContour(99);
              hDraw = h2.get();
        }

        hDraw->Draw("LEGO2");

        /* ── label the colour scale after the palette has been created ──── */
        gPad->Update();

        if (auto* pal =
                    dynamic_cast<TPaletteAxis*>(hDraw->GetListOfFunctions()
                                                  ->FindObject("palette")))
        {
                /* ── move palette farther right ───────────────────────── */
                const double x1 = 0.84;             // left edge
                const double x2 = 0.90;             // right edge
                pal->SetX1NDC(x1);
                pal->SetX2NDC(x2);

                /* ── tidy the vertical axis text (let us label outside) ─ */
                pal->GetAxis()->SetTitle("");       // no vertical title
                pal->GetAxis()->SetLabelSize(0.035);

                /* ── horizontal label centred above palette ───────────── */
                TLatex lab;
                lab.SetNDC();
                lab.SetTextFont(42);
                lab.SetTextAlign(23);               // centred horizontally
                lab.SetTextSize(0.032);
                const double xMid = 0.5*(x1+x2);
                const double yTop = pal->GetY2NDC() + 0.02;
                lab.DrawLatex(xMid, yTop, "Counts");
        }

        /* ── professional header via TLatex ─────────────────────────────── */
        {
              TLatex ttl; ttl.SetNDC();
              ttl.SetTextFont(42);
              ttl.SetTextAlign(23);
              ttl.SetTextSize(0.045);
              ttl.DrawLatex(0.50, 0.94, makeTitle(h->GetName()).c_str());
        }

        c.SaveAs(png.string().c_str());
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] save3D() exception – ")+ex.what());
      }
    }
    

    // ------------------------------------------------------------------
    //  saveProjection – 2‑D projection helper for TH3
    // ------------------------------------------------------------------
    static void saveProjection(TH3* h3,
                               const char* axes,
                               const fs::path& png)
    {
      if (!h3) { ulog::err("[JetQA] saveProjection() nullptr"); return; }
      ulog::trace(std::string("[JetQA] saveProjection(")+axes+") → "+png.string());

      try {
        TH1* tmp = h3->Project3D(axes);                    // ROOT returns TH1*
        auto h2  = std::unique_ptr<TH2>(dynamic_cast<TH2*>(tmp));

        if (!h2) {
          ulog::warn(std::string("[JetQA] projection '")+axes+
                    "' of '"+h3->GetName()+"' is not TH2 – skipped");
          delete tmp;
          return;
        }

        h2->SetDirectory(nullptr);                         // detach
        h2->SetStats(0);
        h2->SetTitle("");                                  // no “yx projection”
          
        /* --- auto‑range identical to handle2D() ------------------------ */
        int fx = h2->GetNbinsX()+1, lx = 0;
        int fy = h2->GetNbinsY()+1, ly = 0;
        double zMin = std::numeric_limits<double>::max();
        for (int ix = 1; ix <= h2->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= h2->GetNbinsY(); ++iy) {
                const double c = h2->GetBinContent(ix,iy);
                if (c <= 0) continue;
                if (ix < fx) fx = ix;  if (ix > lx) lx = ix;
                if (iy < fy) fy = iy;  if (iy > ly) ly = iy;
                if (c  < zMin) zMin = c;
            }
        if (fx < lx) h2->GetXaxis()->SetRange(fx,lx);
        if (fy < ly) h2->GetYaxis()->SetRange(fy,ly);
        if (zMin < std::numeric_limits<double>::max()) h2->SetMinimum(zMin);

        /* translate ROOT’s axis code (“yx”, “xz”, “yz”) → human‑readable text */
        const char* axesDesc =
                      strcmp(axes,"yx")==0 ? "E_{T} vs Area"          :
                      strcmp(axes,"xz")==0 ? "E_{T} vs N_{const}"     :
                      strcmp(axes,"yz")==0 ? "Area vs N_{const}"      :
                      strcmp(axes,"zx")==0 ? "E_{T} vs #eta"          :
                      strcmp(axes,"zy")==0 ? "E_{T} vs #phi"          :
                                             axes;                     /* fallback */

        /* build “axes‑only” title + radius/cent/trigger descriptor */
        std::string full = makeTitle(h3->GetName());         /* e.g.  "Jet E_T … (R = 0.2, …)" */
        std::string descriptor;
        size_t p = full.find('(');                           /* keep everything from '(' onward */
        if (p != std::string::npos) descriptor = full.substr(p);

        std::string finalTitle = std::string(axesDesc) + " " + descriptor;
        h2->SetTitle(finalTitle.c_str());

        saveJetQA2D(h2.get(), png);        // new helper that keeps the header text
      }
      catch (const std::exception& ex) {
        ulog::err(std::string("[JetQA] saveProjection() exception – ")+ex.what());
      }
    }

    // ------------------------------------------------------------------
    //  makeTitle – build a fully descriptive TLatex header
    //              e.g.  "Leading vs sub‑leading jet E_{T}
    //                     (R = 0.2, 20–30 %, MBD #N&S #geq 2 vtx < 10)"
    // ------------------------------------------------------------------
    static std::string makeTitle(const std::string& hname)
    {
        /* 1.  translate the histogram *kind* from its prefix ------------- */
        std::string kind = hname;                // fallback
        if      (hname.rfind("h_maxJetEt_"        ,0) == 0) kind = "Max jet E_{T}";
        else if (hname.rfind("h_leadEt_vs_subEt_" ,0) == 0) kind = "Leading vs subleading jet E_{T}";
        else if (hname.rfind("h_jetEt_area_nConst_",0) == 0)
                kind = "Jet E_{T} vs Area vs N_{const}";
        else if (hname.rfind("h_jetEtEtaPhi_"       ,0) == 0
              || hname.rfind("h_jetEt_eta_phi_"     ,0) == 0)
                kind = "Jet E_{T} vs #eta vs #phi";

        /* 2.  extract radius, optional centrality, trigger --------------- */
        // pattern:  _r02_(20)_(30)_<trigger>   OR   _r02_<trigger>
        std::smatch m;
        std::regex  re(R"(.*_(r[0-9]+|R[0-9]+)(?:_([0-9]+)_([0-9]+))?_(.+))");
        std::string radius = "R = ?";
        std::string cent   = "Centrality IND";
        std::string trig   = "?";

        if (std::regex_match(hname, m, re))
        {
            /* radius ---------------------------------------------------- */
            std::string rTag = m[1];
            int rInt = std::stoi(rTag.substr(1));          // "r02" -> 2
            std::ostringstream rs; rs.setf(std::ios::fixed); rs<<std::setprecision(1)<<rInt/10.0;
            radius = "R = " + rs.str();

            /* centrality (optional) ------------------------------------ */
            if (m[2].matched && m[3].matched)
                cent = m[2].str() + "#minus" + m[3].str() + " %";

            /* trigger --------------------------------------------------- */
            trig = prettifyTrigger(m[4]);
        }

        /* 3.  assemble final title --------------------------------------- */
        return kind + " (" + radius + ", " + cent + ", " + trig + ")";
    }
    
    /* =================================================================
     *  Centrality‑overview writer – one 2×3 canvas for each histogram
     *  kind × jet‑radius × trigger.  Saved under
     *        …/jetQA/generalHistos/<rXX>/centralitySummary/
     * ================================================================= */
    void writeCentralityGrids()
    {
        if (_centCache.empty()) return;

        /* desired ordering of the six bins */
        const std::array<std::string,6> order = {
            "Cent_0_10","Cent_10_20","Cent_20_40",
            "Cent_40_60","Cent_60_80","Cent_80_100"};

        for (auto& [base, vec] : _centCache)
        {
            if (vec.size() < 2) continue;                  /* nothing to plot */

            /* split out radius label and trigger for the title */
            std::smatch m;
            if (!std::regex_match(base, m,
                  std::regex(R"(.*_(r[0-9]+|R[0-9]+)_(.+))"))) continue;
            const std::string rLab  = m[1];
            const std::string trig  = prettifyTrigger(m[2]);

            fs::path dir = root / "jetQA/generalHistos" / rLab / "centralitySummary";
            ensure_dir(dir);

            TCanvas c(("c_"+base).c_str(),"",3*550,2*500);
            c.Divide(3,2,0.001,0.001);

            /* global header (plot type, trigger, radius) */
            c.cd();
            TLatex hd; hd.SetNDC(); hd.SetTextFont(42); hd.SetTextAlign(22);
            hd.SetTextSize(0.038);
            hd.DrawLatex(0.50,0.97,(makeTitle(base)+"  ["+trig+", "+rLab+']').c_str());

            /* draw six pads ------------------------------------------------- */
            int pad = 1;
            for (const auto& key : order)
            {
                c.cd(pad++);
                auto it = std::find_if(vec.begin(), vec.end(),
                                       [&](auto& p){return p.first==key;});
                if (it != vec.end()) {
                    TH1* h = it->second.get();
                    if (h->InheritsFrom(TH2::Class())||
                        h->InheritsFrom(TH3::Class()))
                        h->Draw("COLZ");
                    else
                        h->Draw("HIST");

                    TLatex lb; lb.SetNDC(); lb.SetTextFont(42);
                    lb.SetTextAlign(23); lb.SetTextSize(0.04);
                    lb.DrawLatex(0.50,0.92,
                                 (key.substr(5)+" %").c_str());   /* “lo_hi” → “lo %” */
                } else {
                    TLatex na; na.SetNDC(); na.SetTextFont(42);
                    na.SetTextAlign(22); na.SetTextSize(0.04);
                    na.DrawLatex(0.50,0.50,"N/A");
                }
            }

            c.SaveAs( (dir/(base+"_CentSummary.png")).string().c_str() );
        }
        _centCache.clear();
    }

    /* cache:  base‑histogram‑name  →  {slice, hist‑clone} */
    using SliceClone = std::pair<std::string,std::shared_ptr<TH1>>;
    std::unordered_map<std::string,std::vector<SliceClone>> _centCache;
    
    // containers keeping produced objects alive
    std::vector<std::unique_ptr<TH1>>          _owned1D;
    std::vector<std::unique_ptr<TGraphErrors>> _ownedGraphs;

    /* Inclusive h_maxJetEt_* histograms cached for the destructor */
    std::unordered_map<std::string,std::unique_ptr<TH1>> _cachedMaxEt;
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
    
  enum class Lvl { DBG, INF, WRN, ERR };
  static void log(Lvl l, const std::string& m)
  {
        static const char* tag[]{"DBG","INF","WRN","ERR"};
        std::ostream& os = (l == Lvl::ERR) ? std::cerr : std::cout;
        os << "[VnPlotQA] " << tag[static_cast<int>(l)] << "  " << m << '\n';
  }


  // ─────────────────────────────── 1. cache every TProfile ──────────
  bool process(TObject* o) override
  {
        /* ------------------------------------------------------------------
         *  0. Type gate – accept *only* TProfile objects
         * ------------------------------------------------------------------ */
       if (!o) {
          log(Lvl::ERR, "process(): received nullptr TObject* – aborting");
          return false;
       }

       /* ---------- 1.  Q‑vector scatter plots (TH2F) ------------------ */
       if (o->InheritsFrom(TH2::Class()) &&
          std::string(o->GetName()).rfind("h_Q", 0) == 0)
       {
          auto* h2 = static_cast<TH2*>(o);          // safe – we just checked
          h2->SetDirectory(nullptr);                // detach from ROOT dir

          TCanvas c("c_q", "Q‑vector", 800, 750);
          c.SetGrid();
          h2->SetStats(0);
          h2->Draw("colz");

          /* build …/vNana/Qvectors/<hist>.png  ----------------------- */
          fs::path dst = _outDir / "Qvectors"
                                   / (std::string(h2->GetName()) + ".png");
          try { std::filesystem::create_directories(dst.parent_path()); }
          catch (...) {/* directory already exists – ignore */}

          c.SaveAs(dst.string().c_str());
          log(Lvl::INF, "process(): stored Q‑vector scatter → " + dst.string());
          return true;                               // histogram handled
       }

      /* ---------- 2.  standard v_n TProfiles ------------------------ */
       if (!o->InheritsFrom(TProfile::Class())) return false;  // ignore others

        /* ------------------------------------------------------------------
         *  1. Histogram‑name validation
         *      accepted patterns:
         *          – calorimeters : p_v<n>_<DET>_<lo>_<hi>_<trig>   (DET = CEMC_S …)
         *          – jets        : p_v<n>_JET_rXX_<lo>_<hi>_<trig>
         * ------------------------------------------------------------------ */
        static const std::regex re(
            R"(p_v([123])_((?:[A-Za-z0-9]+_[NST])|(?:JET_r[0-9]{2}))_([0-9]+)_([0-9]+)_(.+))");


        const std::string hName = o->GetName();
        std::smatch       m;
        if (!std::regex_match(hName, m, re)) return false;          // silently skip

        /* ------------------------------------------------------------------
         *  2. Extract metadata from regex capture groups
         * ------------------------------------------------------------------ */
        int         nHarm   = 0;
        std::string detTag, centKey, trigLab;
        try {
            nHarm   = std::stoi(m[1]);                  // 1 / 2 / 3
            detTag  = m[2];                             // CEMC_S …  OR  JET_r04
            centKey = m[3].str() + '_' + m[4].str();    // e.g. 10_20
            trigLab = m[5];
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR, "process(): regex‑capture conversion failed – " +
                          std::string(ex.what()) + "  |  hist = \"" + hName + '"');
            return false;
        }

        /* ------------------------------------------------------------------
         *  3. Clone the profile so the pointer stays valid after walkDir()
         * ------------------------------------------------------------------ */
        std::unique_ptr<TProfile> owned;
        try {
            owned.reset( static_cast<TProfile*>(
                o->Clone( (hName + "_copy").c_str() )) );

            if (!owned) {
                log(Lvl::ERR, "process(): Clone() returned nullptr – hist \"" + hName + '"');
                return false;
            }

            owned->SetDirectory(nullptr);               // detach from any TDirectory
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR, "process(): Clone() threw – " + std::string(ex.what()) +
                          "  |  hist = \"" + hName + '"');
            return false;
        }

        /* ------------------------------------------------------------------
         *  4. Commit the stable pointer into the cache
         * ------------------------------------------------------------------ */
        try {
            TProfile* pStable = owned.get();             // will remain valid
            _ownedProfiles.push_back(std::move(owned));  // ownership kept here

            auto& vec = _raw[trigLab][detTag][nHarm][centKey];
            vec.push_back(pStable);

            log(Lvl::DBG, "process(): cached \"" + hName + "\"  →  trig=" + trigLab +
                          " det=" + detTag + " n=" + std::to_string(nHarm) +
                          " cent=" + centKey +
                          "   |   slice size now = " + std::to_string(vec.size()));
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR, "process(): exception while caching – " + std::string(ex.what()) +
                          "  |  hist = \"" + hName + '"');
            return false;
        }

        return true;   // histogram consumed successfully
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
      if (tag == "HCAL_S" || tag == "HCAL_N" || tag == "HCAL_T")
          return "totalHCal";
      if (tag == "ALL_S"  || tag == "ALL_N"  || tag == "ALL_T")
          return "totalCalo";

      /* jets ---------------------------------------------------------- */
      if (tag.rfind("JET_r",0)==0)                 return "jetvN";
      return tag;                                  // fallback
  }

    static std::string regionOf(const std::string& tag)
    {
        if (tag.rfind("JET_r",0)==0)              // jets keep the “rXX” code
            return tag.substr(4,3);               // e.g. r04
        const char suf = tag.back();              // last char = S / N / T
        if (suf == 'S') return "South";
        if (suf == 'N') return "North";
        return "Total";                           // *anything* else → combined arm
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
         const double val = rProf->GetBinContent(bin);     // ⟨cos n ΔΨ⟩ (can be ≤1)
         if (val > 0.0)                                    // guard against negatives / empty bin
              Rn = std::sqrt(val);                          // resolution  Rn = √⟨cos n ΔΨ⟩

      }

      /* copy points --------------------------------------------------- */
      for (int i=1;i<=nb;++i)
      {
          const double xLo = p->GetXaxis()->GetBinLowEdge(i);
          const double xHi = p->GetXaxis()->GetBinUpEdge (i);
          const double x   = 0.5*(xLo + xHi);
          const double ex  = 0.0;                          // suppress horizontal errors

          double y  = p->GetBinContent(i) / Rn;            // corrected
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
        log(Lvl::INF,"writeCanvases(): enter – regrouping cached profiles");

        /* ---------------- regroup cache:  trigger → detector → region ---------------- */
        using ProfVec = std::vector<TProfile*>;
        using CentMap = std::map<std::string,ProfVec>;          // key = "lo_hi"
        using HarmMap = std::map<int,CentMap>;                  // key = harmonic n
        using RegMap  = std::unordered_map<std::string,HarmMap>;/* region  (N / S / r04 …) */
        using DetMap  = std::unordered_map<std::string,RegMap>; /* detector (EMCal, jetvN) */

        std::unordered_map<std::string,DetMap> grp;             // regrouped by trigger

        /* ---------- 1. regroup profiling histograms --------------------- */
        try
        {
            for (auto& [trig, detMap] : _raw)
            for (auto& [detTag, hMap] : detMap)
            {
                const std::string det = detBase(detTag);
                const std::string reg = regionOf(detTag);
                RegMap& rmap         = grp[trig][det];

                for (auto& [n, cMap] : hMap)
                for (auto& [cent, v] : cMap)
                {
                    {
                        std::ostringstream oss;
                        oss << "  regroup: trig=" << trig
                            << " det=" << det << " reg=" << reg
                            << " n=" << n << " cent=" << cent
                            << " vec=" << v.size();
                        log(Lvl::DBG, oss.str());
                    }

                    auto& tgt = rmap[reg][n][cent];

                    /* ===== DEBUG pre‑state (allocation‑free) ===== */
                    {
                        char buf[256];
                        std::snprintf(buf, sizeof(buf),
                                      "    • tgt=%p  size=%zu  cap=%zu | v=%p  size=%zu  cap=%zu",
                                      static_cast<const void*>(&tgt), tgt.size(), tgt.capacity(),
                                      static_cast<const void*>(&v),   v.size(),   v.capacity());
                        log(Lvl::DBG, buf);
                    }


                    /* ===== alias‑check ===== */
                    if (&tgt == &v)
                    {
                        log(Lvl::ERR,"    ⛔ self‑insertion detected – skipping this vector");
                        continue;                                    // move on to next (cent,v) pair
                    }

                    /* ===== merge profile lists without risking self‑aliasing ========= */
                    try
                    {
                        if (v.empty())           // nothing to do
                            continue;

                        /* First encounter of this centrality --------------------------- *
                         * Just copy‑assign – no reallocation, no aliasing hazard.        */
                        if (tgt.empty())
                        {
                            tgt = v;             // deep copies the pointer list
                            log(Lvl::DBG, "    ✓ first copy – new size=1");
                            continue;
                        }

                        /* Subsequent fills – append safely ----------------------------- */
                        const std::size_t oldCap = tgt.capacity();
                        tgt.reserve(tgt.size() + v.size());   // single reallocation at most
                        tgt.insert(tgt.end(), v.begin(), v.end());

                        std::ostringstream ok;
                        ok << "    ✓ merged " << v.size()
                           << " profiles   new size=" << tgt.size()
                           << "   cap " << oldCap << " → " << tgt.capacity();
                        log(Lvl::DBG, ok.str());
                    }
                    catch (const std::exception& ex)
                    {
                        log(Lvl::ERR, std::string("    ❌ merging vectors failed – ") + ex.what());
                    }
                }
            }
        }
        catch(const std::exception& ex)
        {
            log(Lvl::ERR,std::string("regrouping failed – ")+ex.what());
            return;
        }

        /* ---------- 2. build missing “Total” regions (N ⊕ S) ------------- */
        for (auto& [trig, detMap] : grp)
        for (auto& [det , regMap] : detMap)
        {
            const bool hasN = regMap.count("North");
            const bool hasS = regMap.count("South");
            const bool hasT = regMap.count("Total");

            if (hasT || !(hasN && hasS)) continue;           // nothing to do

            auto& northMap = regMap["North"];
            auto& southMap = regMap["South"];
            auto& totMap   = regMap["Total"];                // will be created

            for (auto& [n, centN] : northMap)               // loop harmonics
            {
                auto& centS = southMap[n];
                for (auto& [centKey, vecN] : centN)          // loop centralities
                {
                    if (vecN.empty()) continue;
                    auto itS = centS.find(centKey);
                    if (itS == centS.end() || itS->second.empty()) continue;

                    /* 1. clone the first North profile ----------------------- */
                    auto* pTot = static_cast<TProfile*>(vecN.front()->Clone());
                    pTot->SetDirectory(nullptr);

                    /* 2. add the first South profile ------------------------- */
                    pTot->Add(itS->second.front());          // TH1::Add = weighted sum

                    /* 3. keep ownership so ROOT will not delete it ----------- */
                    _ownedProfiles.emplace_back(pTot);

                    /* 4. expose to the “Total” container --------------------- */
                    totMap[n][centKey].push_back(pTot);
                }
            }
        }

        /* ---------- 3a. colour palette ----------------------------------- */
        const int colTbl[] = {kRed+1,kBlue+2,kGreen+2,kMagenta+2,kCyan+2,kOrange+1};
        const int nCol     = sizeof(colTbl)/sizeof(int);

        /* ---------- 3b. iterate over regrouped structure ----------------- */
        for (const auto& [trig, detMap] : grp)
        {
            log(Lvl::INF,"trigger group: "+trig);
            for (const auto& [det, regMap] : detMap)
            for (const auto& [reg, harmMap] : regMap)
            for (const auto& [n, centMap]  : harmMap)
            {

                {
                    std::ostringstream oss;
                    oss << "  → det=" << det
                        << "  reg=" << reg
                        << "  n=" << n
                        << "  centSlices=" << centMap.size();
                    log(Lvl::INF, oss.str());
                }

                /* =========== (1)  v_n(pT)  – one canvas per centrality ============ */
                for (const auto& [cent, vec] : centMap)
                {
                    if (vec.empty()) { log(Lvl::WRN,"   slice "+cent+" is EMPTY"); continue; }

                    TCanvas c("c","",1100,850); c.SetGrid();

                    /* create the graph; skip the slice completely if it has no points */
                    auto gHold = graphFromProf(vec.front(), n, cent, trig);

                    if (!gHold || gHold->GetN() == 0)
                    {
                        log(Lvl::WRN, "      ↳ graph has zero points – canvas not saved for cent="
                                      + cent);
                        continue;                                   // nothing to draw
                    }

                    TGraphErrors* g = gHold.get();

                    /* ----------  human‑readable labels  -------------------------------- */
                    std::string centLab = "Centrality IND";
                    {
                        const size_t pos = cent.find('_');
                        if (pos != std::string::npos)
                            centLab = cent.substr(0,pos) + "–" + cent.substr(pos+1) + " %";
                    }
                    const std::string trigLab = prettifyTrigger(trig);

                    /* ----------  set new title  --------------------------------------- */
                    g->SetTitle(Form("%s v_{%d}, %s, %s; p_{T}^{tower}  [GeV]; v_{%d}",
                                     det.c_str(), n, centLab.c_str(), trigLab.c_str(), n));
                    g->Draw("AP");


                    /* ── dynamic Y‑axis: always show 0 and full error bars ────────────── */
                    double yMin = 0.0, yMax = 0.0;
                    for (int ip = 0; ip < g->GetN(); ++ip)
                    {
                        const double y  = g->GetY()[ip];
                        const double ey = g->GetEY()[ip];
                        yMin = std::min(yMin, y - ey);
                        yMax = std::max(yMax, y + ey);
                    }
                    /* ----------  attach axis titles & rescale frame ------------------- */
                    TH1* fr = g->GetHistogram();                       // robust – never nullptr
                    if (fr)
                    {
                        fr->SetMinimum(yMin < 0 ? 1.15 * yMin : 0.0);
                        fr->SetMaximum(1.15 * yMax);
                        fr->GetXaxis()->SetTitle("p_{T}^{tower}  [GeV]");
                        fr->GetYaxis()->SetTitle(Form("v_%d", n));
                        if (yMin < 0)
                        {
                            auto* z = new TLine(fr->GetXaxis()->GetXmin(), 0,
                                                fr->GetXaxis()->GetXmax(), 0);
                            z->SetLineStyle(kDashed); z->Draw();
                        }
                    }
                    c.Modified();  c.Update();

                    gHold.release();                               // ROOT canvas owns it now

                    saveCanvas(c, {det, reg, Form("v%d", n), "Cent_" + cent},
                               Form("v%d_%s_%s_cent%s.png",
                                    n, det.c_str(), reg.c_str(), cent.c_str()));
                    log(Lvl::DBG, "      saved per‑cent canvas for cent=" + cent);
                }

                /* =========== (2)  all‑centralities overlay ======================= */
                try
                {
                    log(Lvl::INF, "   ▶ building all‑centrality overlay for "
                                   + det + " " + reg + "  v" + std::to_string(n));

                    TCanvas cAll("c_all", "", 1100, 850);
                    cAll.SetGrid();

                    TLegend leg(0.15, 0.70, 0.45, 0.88);
                    leg.SetBorderSize(0);

                    int    colIdx = 0;
                    double yMax   = 0.0;
                    double yMin   = 0.0;

                    for (const auto& [cent, vec] : centMap)
                    {
                        std::ostringstream msg;
                        msg << "     • slice " << cent << "  vec.size=" << vec.size();
                        log(Lvl::DBG, msg.str());

                        if (vec.empty())
                        {
                            log(Lvl::WRN, "       ↳ empty slice – skipped");
                            continue;
                        }

                        std::unique_ptr<TGraphErrors> gPtr;
                        try
                        {
                            gPtr = graphFromProf(vec.front(), n, cent, trig);
                        }
                        catch (const std::exception& ex)
                        {
                            log(Lvl::ERR, "       ↳ graphFromProf failed – " + std::string(ex.what()));
                            continue;
                        }

                        const int col = colTbl[colIdx % nCol];
                        ++colIdx;

                        TGraphErrors* g = gPtr.get();
                        g->SetLineColor(col);
                        g->SetMarkerColor(col);
                        g->SetBit(kCanDelete, kFALSE);              // canvas must NOT delete it

                        g->SetTitle(Form("%s v_{%d}, Centrality IND, %s; p_{T}^{tower}  [GeV]; v_{%d}",
                                         det.c_str(), n, prettifyTrigger(trig).c_str(), n));

                        log(Lvl::DBG, "       ↳ drawing, colour=" + std::to_string(col));
                        g->Draw(colIdx == 1 ? "AP" : "P SAME");      // points only, no lines

                        leg.AddEntry(g, Form("Cent %s %%", cent.c_str()), "p");

                        double localMax = -1e9, localMin =  1e9;
                        for (int ip = 0; ip < g->GetN(); ++ip)
                        {
                            const double y  = g->GetY()[ip];
                            const double ey = g->GetEY()[ip];
                            localMax = std::max(localMax, y + ey);
                            localMin = std::min(localMin, y - ey);
                        }
                        yMax = std::max(yMax, localMax);
                        yMin = std::min(yMin, localMin);

                        _ownedGraphs.push_back(std::move(gPtr));
                    }

                    if (colIdx == 0)
                    {
                        log(Lvl::WRN, "   ⚠ no centrality slices drawn – overlay canvas skipped");
                        return;
                    }

                    /* ----------  pick the frame histogram from the FIRST graph -------- */
                    TH1* fr = nullptr;
                    for (auto obj : *cAll.GetListOfPrimitives())          // first TH1 is the frame
                        if ((fr = dynamic_cast<TH1*>(obj))) break;

                    if (fr)
                    {
                        fr->SetMaximum(1.15 * yMax);
                        fr->SetMinimum(yMin < 0 ? 1.15 * yMin : 0.0);
                        fr->GetXaxis()->SetTitle("p_{T}^{tower}  [GeV]");
                        fr->GetYaxis()->SetTitle(Form("v_%d", n));
                        if (yMin < 0)
                        {
                            auto* z = new TLine(fr->GetXaxis()->GetXmin(), 0,
                                                fr->GetXaxis()->GetXmax(), 0);
                            z->SetLineStyle(kDashed); z->Draw();
                        }
                    }
                    cAll.Modified();  cAll.Update();

                    leg.Draw();

                    saveCanvas(cAll,
                               {det, reg, Form("v%d", n), "summaryPlots"},
                               Form("v%d_%s_%s_allCent.png",
                                    n, det.c_str(), reg.c_str()));

                    log(Lvl::INF, "   ✔ overlay saved: v"
                                  + std::to_string(n) + "_" + det + "_" + reg + "_allCent.png");
                }
                catch (const std::exception& ex)
                {
                    log(Lvl::ERR, std::string("   ✖ overlay plot failed – ") + ex.what());
                }


                /* =========== (3)  pT‑integrated v̅_n versus centrality =========== */
                try
                {
                    auto gCent = std::make_unique<TGraphErrors>();
                    int ip = 0;

                    for (const auto& [cent, vec] : centMap)
                    {
                        if (vec.empty()) continue;

                        const size_t pos = cent.find('_');
                        const double cMid = 0.5 * ( std::stod(cent.substr(0,pos)) +
                                                    std::stod(cent.substr(pos+1)) );

                        const std::string rName =
                            Form("p_R%d_vs_cent_%s", n, trig.c_str());
                        double Rn = 1.0;                                   // default = no correction
                        if (auto* r = static_cast<TProfile*>(gROOT->FindObject(rName.c_str())))
                        {
                            const double val = r->GetBinContent(r->FindBin(cMid));
                            if (val > 0.0) Rn = std::sqrt(val);            // use √⟨cos n ΔΨ⟩
                        }

                        const auto [mu,er] = meanProf(vec.front(), Rn);

                        gCent->SetPoint     (ip, cMid, mu);
                        const double ex = 0.0;                       // no horizontal error bar
                        gCent->SetPointError(ip, ex, er);
                        
                        ++ip;
                    }

                    if (gCent->GetN() > 0)
                    {
                        TCanvas cCent("c_cent","",1000,800); cCent.SetGrid();
                        gCent->SetTitle(Form("%s #bar{v}_{%d} vs centrality, %s; centrality  [%%]; v_{%d}",
                                             det.c_str(), n, prettifyTrigger(trig).c_str(), n));
                        gCent->SetMarkerStyle(kFullCircle); gCent->SetLineWidth(2);
                        gCent->Draw("AP");

                        /* ── make sure lower edge shows 0 (or below if needed) -------------- */
                        double yMin = 0.0, yMax = 0.0;
                        for (int ip = 0; ip < gCent->GetN(); ++ip)
                        {
                            const double y  = gCent->GetY()[ip];
                            const double ey = gCent->GetEY()[ip];
                            yMin = std::min(yMin, y - ey);
                            yMax = std::max(yMax, y + ey);
                        }
                        TH1* fr = gCent->GetHistogram();
                        if (fr)
                        {
                            fr->SetMaximum(1.15 * yMax);
                            fr->SetMinimum(yMin < 0 ? 1.15 * yMin : 0.0);
                            fr->GetXaxis()->SetTitle("centrality  [%]");
                            fr->GetYaxis()->SetTitle(Form("v_%d", n));
                            if (yMin < 0)
                            {
                                auto* z = new TLine(fr->GetXaxis()->GetXmin(), 0,
                                                    fr->GetXaxis()->GetXmax(), 0);
                                z->SetLineStyle(kDashed); z->Draw();
                            }
                        }
                        cCent.Modified();  cCent.Update();

                        _ownedGraphs.push_back(std::move(gCent));

                        saveCanvas(cCent,{det,reg,Form("v%d",n),"summaryPlots"},
                                   Form("vbar%d_%s_%s_vsCent.png",
                                        n,det.c_str(),reg.c_str()));
                        
                        log(Lvl::DBG,"      saved v̅_n vs cent");
                    }
                    else
                        log(Lvl::WRN,"      no points for v̅_n vs cent – skipped");
                }
                catch(const std::exception& ex)
                {
                    log(Lvl::ERR,std::string("   v̅_n vs cent failed – ")+ex.what());
                }
            } /* end harmonic‑loop */
        }     /* end trigger‑loop */

      log(Lvl::INF,"writeCanvases(): exit OK");
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




// ╔══════════════════════════════════════════════╗
// ║               T r i g g e r  Q A             ║
// ╚══════════════════════════════════════════════╝
class TriggerQA : public QA
{
    enum class Lvl { DBG, INF, WRN, ERR };
    static void log(Lvl l, const std::string& m)
    {
        static const char* tag[]{"DBG","INF","WRN","ERR"};
        std::ostream& os = (l==Lvl::ERR) ? std::cerr : std::cout;
        os << "[TriggerQA] " << tag[static_cast<int>(l)] << "  " << m << '\n';
    }

public:
    using QA::QA;                       // inherit constructors

    // ────────────────────────────────────────────────────────────────
    // 1. per‑histogram processing – **verbose, exception‑safe version**
    // ────────────────────────────────────────────────────────────────
    bool process(TObject* o) override
    {
        /* ---------------------------------------------------------- *
         * 0.  Preliminary sanity checks                              *
         * ---------------------------------------------------------- */
        if (!o) {
            log(Lvl::ERR,"process(): received nullptr – skipped");
            return false;
        }
        if (!o->InheritsFrom(TH1::Class())) {          // we process only TH1/2
            log(Lvl::DBG,Form("process(): \"%s\" is not TH1 – skipped",o->GetName()));
            return false;
        }

        const std::string n = o->GetName();
        log(Lvl::DBG,"process(): start \""+n+'\"');

        try
        {
            /* ============================================================= *
             *  A.  Minimum‑bias × Trigger correlation matrix (TH2)          *
             * ============================================================= */
            if (n == "h_MB_vs_Trigger")
            {
                auto* h2 = dynamic_cast<TH2*>(o);
                if (!h2) {
                    log(Lvl::ERR,"process(): \""+n+"\" – dynamic_cast<TH2*> failed");
                    return false;
                }

                fs::path out = root / "triggerQA" / "MBTrigCorrelation.png";
                ensure_dir(out.parent_path());

                log(Lvl::DBG,"process(): rendering MB×Trigger correlation → "+out.string());

                TCanvas c("c_corr","MB vs Trigger decision",1200,700);
                c.SetRightMargin(0.18);
                h2->SetStats(0);
                h2->SetTitle("Minimum‑bias vs Trigger decision;Trigger key;Event category");
                c.SetLogz();
                h2->Draw("COLZ TEXT");

                /* ---- totals and accepted fraction ------------------------ */
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03); tx.SetTextAlign(13);
                const double padL=c.GetLeftMargin(), padT=c.GetTopMargin();
                const int nTrig = h2->GetNbinsX();
                long long totEvt=0, bothOK=0;

                for (int ix=1; ix<=nTrig; ++ix)
                {
                    bothOK += static_cast<long long>(h2->GetBinContent(ix,4));
                    for (int iy=1; iy<=4; ++iy)
                        totEvt += static_cast<long long>(h2->GetBinContent(ix,iy));
                }
                tx.DrawLatex(padL+0.02,1-padT-0.04,
                             Form("Total events : %lld",totEvt));
                tx.DrawLatex(padL+0.02,1-padT-0.08,
                             Form("Accepted (MB && Trig) : %lld  (%.2f%%)",
                                  bothOK, (totEvt ? 100.*bothOK/totEvt : 0.0)));

                c.SaveAs(out.string().c_str());
                log(Lvl::INF,"wrote "+out.string());
                return true;                                // finished branch A
            }

            /* ============================================================= *
             *  B.  Scalar counters & turn‑on spectra                         *
             * ============================================================= */

            const bool isCnt =
                  n.rfind("cnt_",0)                         == 0  ||   // raw/live/scaled bit counters
                  n.rfind("h_vtxRelToCut_",0)               == 0  ||   // vertex‑quality monitor
                  n.rfind("h_maxClusterEnergy_doNotScale_",0)== 0 ||   // turn‑on  (raw+live)
                  n.rfind("h_maxClusterE_EMC_",0)           == 0;      // turn‑on  (scaled)

            if (!isCnt) {
                log(Lvl::DBG,"process(): \""+n+"\" – not a counter/turn‑on histo");
                return false;
            }

            auto* h1 = dynamic_cast<TH1*>(o);
            if (!h1) {
                log(Lvl::ERR,"process(): \""+n+"\" – dynamic_cast<TH1*> failed");
                return false;
            }


            /* corrected code – four functional changes
             *   (1) strip *all* leading zeros after “…_lt_”
             *   (2) ignore empty histograms (Integral() == 0)
             *   (3) never overwrite an existing entry (keeps merged content)
             *   (4) log which keys are stored ─ easier debugging                */
            auto normalisedKey = [](std::string raw) -> std::string
            {
                /* “…_lt_00025” -> “…_lt_25”, “…_lt_010” -> “…_lt_10”, … */
                return std::regex_replace(raw, std::regex(R"(_lt_0+([1-9]\d*))"), "_lt_$1");
            };

            /* do ‑ N O T ‑ scale -------------------------------------------------- */
            if (n.rfind("h_maxClusterEnergy_doNotScale_",0)==0)
            {
                std::string trgName =
                    normalisedKey(n.substr(std::strlen("h_maxClusterEnergy_doNotScale_")));

                if (!_spectra.count(trgName))                       // ← no “Integral()” filter
                {
                    TH1* hClone = static_cast<TH1*>(h1->Clone());
                    hClone->SetDirectory(nullptr);
                    _spectra[trgName] = hClone;
                    log(Lvl::DBG,"process(): cached doNotScale \""+trgName+
                                  "\"  (Integral=" + std::to_string(h1->Integral()) + ')');
                }
            }

            /* s c a l e d  -------------------------------------------------------- */
            else if (n.rfind("h_maxClusterE_EMC_",0)==0)
            {
                std::string trgName =
                    normalisedKey(n.substr(std::strlen("h_maxClusterE_EMC_")));

                if (!_spectraScaled.count(trgName))
                {
                    TH1* hClone = static_cast<TH1*>(h1->Clone());
                    hClone->SetDirectory(nullptr);
                    _spectraScaled[trgName] = hClone;
                    log(Lvl::DBG,"process(): cached scaled \""+trgName+
                                  "\"  (Integral=" + std::to_string(h1->Integral()) + ')');
                }
            }

            /* ---- (B‑2) decide subfolder ---------------------------------- */
            std::string subdir;
            if (n.find("_raw")    != std::string::npos) subdir = "bitCounts/raw";
            else if (n.find("_live")   != std::string::npos) subdir = "bitCounts/live";
            else if (n.find("_scaled") != std::string::npos) subdir = "bitCounts/scaled";
            else if (n.rfind("h_vtxRelToCut_",0)==0)         subdir = "vtxCutCompliance";
            else                                            subdir = "misc";

            fs::path out = root / "triggerQA" / subdir / (n + ".png");
            ensure_dir(out.parent_path());

            /* ---- (B‑3) render the histogram ------------------------------ */
            log(Lvl::DBG,"process(): rendering counter → "+out.string());
            TCanvas c("c_cnt","",800,600);
            h1->SetStats(0);
            h1->SetLineWidth(2);
            const double yMax = 1.15 * h1->GetMaximum();
            h1->SetMaximum(yMax>0 ? yMax : 1.);
            h1->Draw("HIST TEXT00");

            c.SaveAs(out.string().c_str());
            log(Lvl::DBG,"saved "+out.string());
            return true;                                    // finished branch B
        }
        /* -------------------------------------------------------------- *
         * 4.  Robust exception handling                                   *
         * -------------------------------------------------------------- */
        catch (const std::exception& ex)
        {
            log(Lvl::ERR,std::string("process(): C++ exception for \"")+n+
                          "\" – "+ex.what());
            return false;
        }
        catch (...)
        {
            log(Lvl::ERR,std::string("process(): UNKNOWN exception for \"")+n+"\"");
            return false;
        }
    }

    /* ──────────────────────────────────────────────────────────────────
     * REPLACEMENT — paste this entire block verbatim
     * ────────────────────────────────────────────────────────────────── */
    ~TriggerQA() override
    {
        /* --------------------------------------------------------------
         * 0.  Entry banner & one‑time guard
         * -------------------------------------------------------------- */
        const std::string runID = root.parent_path().filename().string();   // e.g. "69420", "Combined"
        log(Lvl::DBG, "~TriggerQA()  –  enter  (runID = \"" + runID + "\")");

        if (runID != "Combined") {
            log(Lvl::DBG, "~TriggerQA()  –  non‑Combined pass → nothing to do");
            return;
        }
        if (_emitted) {
            log(Lvl::WRN, "~TriggerQA()  –  already executed once, skipping");
            return;
        }
        _emitted = true;

        /* --------------------------------------------------------------
         * 1.  Prepare unified output directory
         * -------------------------------------------------------------- */
        fs::path outDir = root.parent_path() / "turnOnCurveSummaries";
        try {
            ensure_dir(outDir);
            log(Lvl::INF, "Created / verified summary folder:  " + outDir.string());
        }
        catch (const std::exception& ex) {
            log(Lvl::ERR, std::string("ensure_dir(\"") + outDir.string() +
                          "\") failed – " + ex.what());
            return;
        }

        /* --------------------------------------------------------------
         * 2.  Define trigger groups (MB reference + photon hierarchy)
         * -------------------------------------------------------------- */
        static const std::vector<std::string> grp150 = {
            "MBD_NS_geq_2_vtx_lt_150",
            "photon_6_plus_MBD_NS_geq_2_vtx_lt_150",
            "photon_8_plus_MBD_NS_geq_2_vtx_lt_150",
            "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
            "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
        };
        static const std::vector<std::string> grp10  = {
            "MBD_NS_geq_2_vtx_lt_10",
            "photon_6_plus_MBD_NS_geq_2_vtx_lt_10",
            "photon_8_plus_MBD_NS_geq_2_vtx_lt_10",
            "photon_10_plus_MBD_NS_geq_2_vtx_lt_10",
            "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"
        };

        /* --------------------------------------------------------------
         * 3.  Generic overlay helper
         *     – works for *either* cache (do‑not‑scale / scaled)
         * -------------------------------------------------------------- */
        auto makeOverlay =
            [&](const std::vector<std::string>& trgList,
                const std::string& outName,
                const std::unordered_map<std::string,TH1*>& src)
        {
            log(Lvl::DBG, "makeOverlay(\"" + outName + "\") – starting");

            /* 3.1  Check reference MB trigger --------------------------- */
            const std::string& refKey = trgList.front();              // always MB entry
            if (!src.count(refKey)) {
                log(Lvl::WRN, "  MB reference \"" + refKey +
                               "\" missing – overlay skipped");
                return;
            }

            /* 3.2  Build list of available photon triggers ------------- */
            std::vector<std::string> avail;          // keeps original order
            avail.push_back(refKey);                 // MB reference always first

            for (std::size_t i = 1; i < trgList.size(); ++i) {
                if (src.count(trgList[i])) {
                    avail.push_back(trgList[i]);
                    log(Lvl::DBG, "  found histogram for \"" + trgList[i] + '"');
                } else {
                    log(Lvl::DBG, "  histogram for \"" + trgList[i] + "\" is missing");
                }
            }

            if (avail.size() == 1) {
                log(Lvl::WRN, "  no photon trigger histograms present → overlay skipped");
                return;
            }

            /* limiting palette length just in case */
            if (avail.size() > 5)
                log(Lvl::WRN, "  more than 5 curves – excess will re‑use last colour");

            TCanvas c(("c_"+outName).c_str(), "Trigger turn-on", 1100, 900);
            /* one column, two rows – small gap between pads             */
            c.Divide(1, 2, 0.005, 0.005);

            const int palette[5] = {kBlack, kBlue+2, kGreen+2, kRed+1, kMagenta+2};

            /* ----------------  LEFT pad : absolute spectra  ----------------- */
            c.cd(1); gPad->SetLogy();
            double globalMax = 0.0;
            for (auto& key : avail)
                globalMax = std::max(globalMax, src.at(key)->GetMaximum());

            for (std::size_t i = 0; i < avail.size(); ++i) {
                TH1* h = src.at(avail[i]);
                h->SetLineColor(palette[std::min<std::size_t>(i,4)]);
                h->SetLineWidth(2);
                h->SetMaximum(globalMax * 1.15);
                h->Draw(i == 0 ? "HIST" : "HIST SAME");
            }
            TLegend leg(0.52, 0.66, 0.88, 0.88); leg.SetBorderSize(0);
            for (std::size_t i = 0; i < avail.size(); ++i)
                leg.AddEntry(src.at(avail[i]), avail[i].c_str(), "l");
            leg.Draw();

            /* ----------------  RIGHT pad : ratios to MB --------------------- */
            c.cd(2); gPad->SetGridy();
            TH1* hRef = src.at(refKey);
            std::vector<std::unique_ptr<TH1>> ratios;

            for (std::size_t i = 1; i < avail.size(); ++i) {
                try {
                    ratios.emplace_back(
                        static_cast<TH1*>(src.at(avail[i])->Clone())
                    );
                }
                catch (...) {
                    log(Lvl::ERR, "  Clone() failed for \"" + avail[i] + '"');
                    continue;
                }

                ratios.back()->Divide(hRef);
                ratios.back()->SetLineColor(palette[std::min<std::size_t>(i,4)]);
                ratios.back()->SetLineWidth(2);
                ratios.back()->SetTitle(";E_{max}^{cluster}  [GeV];Ratio to MB trigger");

                if (ratios.back()->Integral() > 0.0)          // guard against an all‑zero ratio
                {
                    ratios.back()->SetMaximum(1.20);
                    ratios.back()->SetMinimum(0.00);
                    ratios.back()->Draw(i == 1 ? "HIST" : "HIST SAME");
                }
                else
                {
                    log(Lvl::WRN,"  ratio for \""+avail[i]+"\" is empty – skipped");
                }
            }

            TLegend leg2(0.52, 0.16, 0.88, 0.38); leg2.SetBorderSize(0);
            for (std::size_t i = 1; i < avail.size(); ++i)
                leg2.AddEntry(ratios[i-1].get(), avail[i].c_str(), "l");
            leg2.Draw();

            /* 3.4  Write PNG with robust error handling ---------------- */
            fs::path dst = outDir / outName;
            try {
                c.SaveAs(dst.string().c_str());
                log(Lvl::INF, "  wrote " + dst.string());
            }
            catch (const std::exception& ex) {
                log(Lvl::ERR, std::string("  SaveAs(\"") + dst.string() +
                              "\") failed – " + ex.what());
            }

            /* --------------------------------------------------------------
             *  individual photon‑pair overlays
             * -------------------------------------------------------------- */
            for (std::size_t i = 1; i < avail.size(); ++i)
            {
                const std::string& pho = avail[i];             // current photon trigger
                std::string tag = pho.substr(0, pho.find_first_of(" \t"));
                std::string oneName =
                    outName.substr(0,outName.find_last_of('.')) + "_" + tag + ".png";

                TCanvas c2(("c_"+oneName).c_str(), "Two‑curve turn‑on", 1100, 900);
                c2.Divide(1,2,0.005,0.005);

                /* top : spectra -------------------------------------------------- */
                c2.cd(1); gPad->SetLogy();
                double ymax = std::max(src.at(refKey)->GetMaximum(),
                                       src.at(pho)->GetMaximum());
                for (int pass=0; pass<2; ++pass)
                {
                    TH1* h = src.at(!pass ? refKey : pho);
                    h->SetLineColor(palette[pass]);
                    h->SetLineWidth(2);
                    h->SetMaximum(ymax*1.15);
                    h->Draw(pass==0 ? "HIST" : "HIST SAME");
                }
                TLegend l1(0.55,0.66,0.88,0.88); l1.SetBorderSize(0);
                l1.AddEntry(src.at(refKey), refKey.c_str(), "l");
                l1.AddEntry(src.at(pho),    pho.c_str(),    "l");
                l1.Draw();

                /* bottom : ratio ------------------------------------------------- */
                c2.cd(2); gPad->SetGridy();
                std::unique_ptr<TH1> r(
                    static_cast<TH1*>(src.at(pho)->Clone())
                );
                r->Divide(src.at(refKey));
                r->SetLineColor(palette[1]);
                r->SetLineWidth(2);
                r->SetTitle(";E_{cluster}^{max}  [GeV];Ratio to MB trigger");
                r->SetMaximum(1.20);  r->SetMinimum(0.00);
                r->Draw("HIST");

                fs::path dst2 = outDir / oneName;
                try {
                    c2.SaveAs(dst2.string().c_str());
                    log(Lvl::INF,"  wrote pair‑overlay " + dst2.string());
                }
                catch (const std::exception& ex) {
                    log(Lvl::ERR,"  SaveAs(\""+dst2.string()+"\") failed – "+ex.what());
                }
            }

            log(Lvl::DBG, "makeOverlay(\"" + outName + "\") – done");

        };

        /* --------------------------------------------------------------
         * 4.  Emit up to four summary PNGs
         * -------------------------------------------------------------- */
        makeOverlay(grp150, "TurnOn_doNotScale_vtx_lt_150.png", _spectra);
        makeOverlay(grp10 , "TurnOn_doNotScale_vtx_lt_10.png",  _spectra);
        makeOverlay(grp150, "TurnOn_scaled_vtx_lt_150.png",     _spectraScaled);
        makeOverlay(grp10 , "TurnOn_scaled_vtx_lt_10.png",      _spectraScaled);

        log(Lvl::INF, "~TriggerQA()  –  completed successfully");
    }

private:
    /* one‑per‑trigger cache of the “doNotScale” spectra */
    static inline std::unordered_map<std::string,TH1*> _spectra;        // doNotScale
    static inline std::unordered_map<std::string,TH1*> _spectraScaled;  // per‑event / scaled
    static inline bool _emitted = false;   // only once in Combined pass
};



/* ────────────────────────────────────────────────────────────────────
 *  9.  MAIN DRIVER  –  per‑run analysis in parallel + combined pass
 *      (process‑pool implementation – safe for all ROOT classes)
 * ────────────────────────────────────────────────────────────────── */

namespace  /* helpers stay local to this TU */ {

/* small convenience type from the original code -------------------- */
struct Cnt { int tot = 0, used = 0; };

/* ------------------------------------------------------------------ *
 *  QA‑module filter – controlled via environment variable QA_ONLY.
 *  Example:   export QA_ONLY="correlations,hcal,jetqa"
 * ------------------------------------------------------------------ */
static std::unordered_set<std::string> gQaFilter;

static bool wantQA(const std::string& tag)
{
    /* empty filter  →  accept every module */
    return gQaFilter.empty() || gQaFilter.count(tag);
}

/* one‑time initialiser – runs before main() ------------------------ */
struct _InitQaFilter_
{
    _InitQaFilter_()
    {
        const char* env = std::getenv("QA_ONLY");
        if (!env || !*env) return;                 // no list supplied
        std::stringstream ss(env);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            std::transform(tok.begin(), tok.end(), tok.begin(), ::tolower);
            gQaFilter.insert(tok);
        }
        if (!gQaFilter.empty()) {
            std::ostringstream o;
            o << "QA filter active → ";
            for (const auto& t : gQaFilter) o << t << ' ';
            ulog::info(o.str());
        }
    }
} _initQaFilter_;

/* ===================================================================
 * H‑0  :  open ROOT file + discover centrality slices
 * =================================================================== */
std::unique_ptr<TFile> openInputFile(const std::string& file)
{
    std::unique_ptr<TFile> in( TFile::Open(file.c_str(), "READ") );
    if (!in || in->IsZombie()) { ulog::err("Cannot open " + file); return nullptr; }
    ulog::ok("Input file opened");
    return in;
}

CentList discoverAndReportSlices(TFile* in)
{
    CentList slices = discoverSlices(in);
    slices.erase(std::remove_if(slices.begin(), slices.end(),
                                [](const std::string& s){ return s=="0_100"; }),
                 slices.end());
    {
        std::ostringstream o; o << "Centrality slices: ";
        for (auto& s : slices) o << s << "  ";
        ulog::info(o.str());
    }
    return slices;
}

/* ===================================================================
 * H‑1  :  PASS‑0 – histogram catalogue
 * =================================================================== */
void catalogueHistograms(TFile*             in,
                         const std::string& outBase,
                         const std::string& inFile)
{
    ulog::banner("Pass 0 – Catalogue");

    fs::path catTxt = fs::path(outBase) / "AllHistogramNames.txt";
    ensure_dir(catTxt.parent_path());
    std::ofstream cat(catTxt);
    cat << "# Histogram catalogue for " << inFile << "\n";

    std::unordered_map<string,int> hCnt;     /* preserved, still filled */
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
    ulog::ok("Histogram list written → " + catTxt.string());
}

/* ===================================================================
 * H‑2  :  PASS‑1 – full QA production
 *         Returns the two maps needed later for the summaries.
 *         ▶  This version is functionally identical to the original
 *            but adds exhaustive diagnostics and exception guards.
 * =================================================================== */
struct QaMaps {
    std::unordered_map<std::string,Cnt> stat;       // total/used histograms
    std::unordered_map<std::string,int> runsActive; // scaled‑trigger activity
};

QaMaps runQaProduction(TFile*              in,
                       const std::string&  outBase,
                       const CentList&     slices)
{
    ulog::banner("Pass 1 – QA Production  (outBase = \"" + outBase + "\")");

    //------------------------------------------------------------------
    // 0. CSV initialisation
    //------------------------------------------------------------------
    fs::path csvPath = fs::path(outBase) / "InvariantMassSummary.csv";
    ensure_dir(csvPath.parent_path());

    std::ofstream csv(csvPath);
    if (!csv) {
        ulog::err("Cannot open " + csvPath.string() + " for writing – aborting runQaProduction()");
        return {};
    }
    csv << "trigger,E,Chi,Asym,pTlo,pThi,cent,"
           "meanPi0,errPi0,sigmaPi0,errSigmaPi0,"
           "meanEta,errEta,sigmaEta,errSigmaEta\n";

    //------------------------------------------------------------------
    // 1. state managed across triggers
    //------------------------------------------------------------------
    NSCache<MapPair> mbdCache, sepdCache;
    QaMaps maps;                       // <‑‑ will be returned

    //------------------------------------------------------------------
    // 2. iterate over top‑level trigger directories (one “trg” loop)
    //    • emits start / end banners
    //    • traps every recoverable error
    //    • counts silent failures (errCnt)
    //------------------------------------------------------------------
    TIter itDir(in->GetListOfKeys());
    while (auto* kd = dynamic_cast<TKey*>(itDir()))
    {
        // ── guard: interested only in sub‑directories ──────────────────
        if (strcmp(kd->GetClassName(), "TDirectoryFile") != 0) continue;

        const std::string trg = kd->GetName();
        if (!kTriggersWanted.count(trg)) continue;

        const auto t0 = std::chrono::steady_clock::now();
        ulog::banner("Trigger \"" + trg + "\" — analysis START");

        try
        {
            //------------------------------------------------------------------
            // 2‑A. open trigger directory
            //------------------------------------------------------------------
            std::unique_ptr<TDirectory> trigDir(
                static_cast<TDirectory*>(kd->ReadObj()));

            if (!trigDir || trigDir->IsZombie()) {
                ulog::warn("  ↳ directory unreadable – skipping trigger \"" + trg + '"');
                continue;
            }

            //------------------------------------------------------------------
            // 2‑B. classify folder as “minimal” or “full”
            //------------------------------------------------------------------
            auto isMinimal = [](TDirectory* d)->bool
            {
                TIter it(d->GetListOfKeys());
                while (auto* k = dynamic_cast<TKey*>(it())) {
                    const std::string n = k->GetName();
                    const bool ok =  n == "h_MB_vs_Trigger"          ||
                                     n.rfind("cnt_",0)            == 0 ||
                                     n.rfind("h_vtxRelToCut_",0)  == 0;
                    if (!ok) return false;
                }
                return true;
            };

            //------------------------------------------------------------------
            // 2‑C. live‑trigger detection and directory routing
            //------------------------------------------------------------------
            auto* hCntScaled = dynamic_cast<TH1*>(
                    trigDir->Get(Form("cnt_%s_scaled", trg.c_str())));
            const bool hasLive = (hCntScaled && hCntScaled->GetBinContent(1) > 0);
            if (hasLive) ++maps.runsActive[trg];

            /* inactive triggers → …/scaledTriggersOFF_generalTrigQA/<trg>/ */
            fs::path trgBase = hasLive
                                 ? fs::path(outBase) / trg
                                 : fs::path(outBase) / "scaledTriggersOFF_generalTrigQA" / trg;

            /* decide whether only trigger‑QA is required */
            const bool minimal = isMinimal(trigDir.get()) || !hasLive;

            if (minimal)
            {
                /* MB‑only folder  or  scaled counts = 0 */
                ulog::trace("  ↳ scaled counts = 0 or minimal folder – writing to "
                           + trgBase.string());
                ensure_dir(trgBase / "triggerQA");
            }
            else
            {
                std::vector<std::string> sub = {
                    "correlations", "centrality",
                    "HCal/IHCal", "HCal/OHCal", "HCal/totalHCal",
                    "MBD/otherQA", "MBD/zVertex",
                    "sEPD/OtherQA", "sEPD/EventPlaneQA", "sEPD/tileQA",
                    "jetQA/generalHistos", "jetQA/summary", "triggerQA",
                    "vNana", "EMCal/otherQA",
                    "EMCal/invMassQA", "EMCal/invMassQA/cutQA"
                };

                for (const auto& sd : sub)
                    try {
                        ensure_dir(trgBase / sd);
                    }
                    catch (const std::exception& ex) {
                        ulog::warn("  ↳ cannot create \"" +
                                  (trgBase / sd).string() +
                                  "\" – " + ex.what());
                    }
            }


            //------------------------------------------------------------------
            // 2‑D. assemble QA module stack  (filter‑aware)
            //------------------------------------------------------------------
            std::vector<std::unique_ptr<QA>> qa;

            /* helper: add module only when wanted */
            auto push = [&](const std::string& tag, auto p)
                        { if (wantQA(tag)) qa.emplace_back(std::move(p)); };

            push("correlations", std::make_unique<CorrQA>(trg, trgBase, slices));
            push("hcal",         std::make_unique<HcalQA>(trg, trgBase, slices));
            push("mbd",          std::make_unique<MbdQA >(trg, trgBase, slices, mbdCache));
            push("sepd",         std::make_unique<SepdQA>(trg, trgBase, slices, sepdCache));
            push("sepdother",    std::make_unique<sEPDotherQA>(trg, trgBase, slices));
            push("jetqa",        std::make_unique<JetQA >(trg, trgBase, slices));
            push("eventqa",      std::make_unique<EventQA>(trg, trgBase, slices));
            push("triggerqa",    std::make_unique<TriggerQA>(trg, trgBase, slices));

            if (hasLive) {
                push("pi0",   std::make_unique<Pi0QA >(trg, trgBase, slices, csv));
                push("emcal", std::make_unique<EmcalQA>(trg, trgBase, slices));
                push("vn",    std::make_unique<VnPlotQA>(trg, trgBase, slices));
            }

            ulog::trace("  ↳ instantiated "
                       + std::to_string(qa.size()) + " QA modules");

            //------------------------------------------------------------------
            // 2‑E. recursive walk of all keys inside the trigger directory
            //------------------------------------------------------------------
            size_t errCnt = 0;

            auto walkDir = [&](TDirectory* dir, auto&& walk)->void {
                TIter it(dir->GetListOfKeys());
                while (auto* key = dynamic_cast<TKey*>(it()))
                {
                    // dive into sub‑directories first
                    if (strcmp(key->GetClassName(),"TDirectoryFile") == 0) {
                        walk(static_cast<TDirectory*>(key->ReadObj()), walk);
                        continue;
                    }

                    std::unique_ptr<TObject> obj(key->ReadObj());
                    if (!obj) { ++errCnt; continue; }

                    // detach TH1 from ROOT's dir ownership
                    if (obj->InheritsFrom(TH1::Class()))
                        static_cast<TH1*>(obj.get())->SetDirectory(nullptr);

                    ++maps.stat[trg].tot;
                    bool used = false;
                    for (auto& m : qa)
                        if (m->process(obj.get())) { used = true; break; }

                    if (used) ++maps.stat[trg].used;
                }
            };
            walkDir(trigDir.get(), walkDir);

            //------------------------------------------------------------------
            // 2‑F. end‑of‑trigger log line
            //------------------------------------------------------------------
            const auto dt = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - t0).count();

            ulog::ok("Trigger \"" + trg + "\"  – done in "
                    + std::to_string(dt).substr(0,5) + " s"
                    + "   |   used "
                    + std::to_string(maps.stat[trg].used) + "/"
                    + std::to_string(maps.stat[trg].tot) + " objects"
                    + (errCnt ? ("   |   silent‑failures = " + std::to_string(errCnt)) : ""));
        }
        catch (const std::exception& ex)
        {
            ulog::err("‼  Fatal exception while processing trigger \""
                     + trg + "\": " + ex.what());
        }
        ulog::banner("Trigger \"" + trg + "\" — analysis END");
    } // <- end trigger loop


    csv.close();
    ulog::info("Pass 1 – QA Production finished, CSV written → " + csvPath.string());
    return maps;  // ––––––––––––––––––––––––––––––––––––––––––––––––––
}


/* ===================================================================
 * H‑3  :  summary printouts   (three blocks, unchanged text)
 * =================================================================== */
void printScaledTriggerSummary(const std::unordered_map<std::string,int>& runsActive)
{
    ulog::banner("Scaled‑trigger activity summary");

    std::size_t trigCol = 0;
    for (const auto& [t,_] : runsActive) trigCol = std::max(trigCol, t.size());
    trigCol = std::max<std::size_t>(trigCol, 8);

    const std::string hRule(trigCol + 17, '=');

    std::cout << hRule << "\n"
              << term::CLR_BOLD
              << std::left  << std::setw(trigCol) << "Trigger"
              << " │ "
              << std::right << std::setw(12)      << "nRuns"
              << term::CLR_RST << "\n"
              << hRule << "\n";

    int totalRuns = 0;
    for (const auto& [t,n] : runsActive) {
        std::cout << std::left  << std::setw(trigCol) << t
                  << " │ "
                  << std::right << std::setw(12)     << n << "\n";
        totalRuns += n;
    }

    std::cout << hRule << "\n"
              << term::CLR_BOLD
              << std::left  << std::setw(trigCol) << "TOTAL"
              << " │ "
              << std::right << std::setw(12)     << totalRuns
              << term::CLR_RST << "\n"
              << hRule << "\n";
}

void printMbCorrelationSummary(TFile* in,
                               const std::unordered_map<std::string,Cnt>& stat)
{
    ulog::banner("Trigger ↔ MB correlation");

    std::size_t trigW = 0;
    for (const auto& [t,_] : stat) trigW = std::max(trigW, t.size());
    trigW = std::max<std::size_t>(trigW, 8);

    auto printRule = [&](char ch){ std::cout << std::string(trigW + 75, ch) << "\n"; };

    printRule('=');
    std::cout << term::CLR_BOLD
              << std::left  << std::setw(trigW) << "Trigger"
              << " │ " << std::right << std::setw(12) << "RAW"
              << " │ " << std::setw(12)               << "LIVE"
              << " │ " << std::setw(12)               << "SCALED"
              << " │ " << std::setw(12)               << "MB&&Trig"
              << " │ " << std::setw(12)               << "MB Only"
              << " │ " << std::setw(12)               << "Trig Only"
              << term::CLR_RST << "\n";
    printRule('-');

    long long totRaw = 0, totLive = 0, totScaled = 0,
              totBoth = 0, totMBOnly = 0, totTrigOnly = 0;

    TIter itTrig(in->GetListOfKeys());
    while (auto* kDir = dynamic_cast<TKey*>(itTrig())) {
        if (strcmp(kDir->GetClassName(),"TDirectoryFile")) continue;
        std::string trgName = kDir->GetName();
        if (!kTriggersWanted.count(trgName)) continue;

        TDirectory* dTrig = static_cast<TDirectory*>(kDir->ReadObj());

        auto getCnt = [&](const std::string& h)->long long {
            if (auto* h1 = dynamic_cast<TH1*>( dTrig->Get(h.c_str()) ))
                return static_cast<long long>(h1->GetBinContent(1));
            return 0LL;
        };
        const long long nRaw    = getCnt("cnt_" + trgName + "_raw");
        const long long nLive   = getCnt("cnt_" + trgName + "_live");
        const long long nScaled = getCnt("cnt_" + trgName + "_scaled");

        long long mb_and_trig = 0, mb_only = 0, trig_only = 0;
        if (auto* h2 = dynamic_cast<TH2*>( dTrig->Get("h_MB_vs_Trigger") )) {
            int xbin = -1;
            for (int ix = 1; ix <= h2->GetNbinsX(); ++ix)
                if (const char* lb = h2->GetXaxis()->GetBinLabel(ix);
                    lb && std::string(lb) == trgName) { xbin = ix; break; }
            if (xbin > 0) {
                trig_only   = static_cast<long long>(h2->GetBinContent(xbin,2));
                mb_only     = static_cast<long long>(h2->GetBinContent(xbin,3));
                mb_and_trig = static_cast<long long>(h2->GetBinContent(xbin,4));
            }
        }

        totRaw      += nRaw;
        totLive     += nLive;
        totScaled   += nScaled;
        totBoth     += mb_and_trig;
        totMBOnly   += mb_only;
        totTrigOnly += trig_only;

        const char* rowClr =
            (mb_and_trig > 0 && mb_only == 0 && trig_only == 0)
                ? term::CLR_GRN : term::CLR_CYAN;

        std::cout << rowClr
                  << std::left  << std::setw(trigW) << trgName
                  << " │ " << std::right << std::setw(12) << nRaw
                  << " │ " << std::setw(12)               << nLive
                  << " │ " << std::setw(12)               << nScaled
                  << " │ " << std::setw(12)               << mb_and_trig
                  << " │ " << std::setw(12)               << mb_only
                  << " │ " << std::setw(12)               << trig_only
                  << term::CLR_RST << "\n";
    }

    printRule('=');
    std::cout << term::CLR_BOLD
              << std::left  << std::setw(trigW) << "TOTAL"
              << " │ " << std::right << std::setw(12) << totRaw
              << " │ " << std::setw(12)               << totLive
              << " │ " << std::setw(12)               << totScaled
              << " │ " << std::setw(12)               << totBoth
              << " │ " << std::setw(12)               << totMBOnly
              << " │ " << std::setw(12)               << totTrigOnly
              << term::CLR_RST << "\n";
    printRule('=');

    ulog::ok("All outputs under " + kOutputBase);
}

void printRunEventSummary()
{
    using term::CLR_BOLD; using term::CLR_RST;
    const auto& ev = EventQA::eventCounts();
    if (ev.empty()) return;

    std::vector<std::pair<int,long long>> rows;
    rows.reserve(ev.size());

    for (const auto& [runStr,n] : ev)
        if (std::all_of(runStr.begin(), runStr.end(), ::isdigit))
            rows.emplace_back(std::stoi(runStr), n);

    std::sort(rows.begin(), rows.end(),
              [](auto a, auto b){ return a.first < b.first; });

    ulog::banner("Run‑by‑run event statistics");
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

    ulog::ok("Runs analysed : " + std::to_string(rows.size()) +
            "   |   Total events : " + std::to_string(totalEv));
}

/* ===================================================================
 * H‑4  :  helper that prints the MB‑/trigger‑summary + events
 * =================================================================== */
void printAllSummaries(TFile* in, const QaMaps& maps)
{
    printScaledTriggerSummary(maps.runsActive);
    printMbCorrelationSummary(in, maps.stat);
    printRunEventSummary();
}

} // anonymous namespace



/* ------------------------------------------------------------------ */
/*  One complete QA pass for a single ROOT file                        */
/* ------------------------------------------------------------------ */
void runOneQaPass(const std::string& inFile,
                  const std::string& outBase)
{
    /* original globals remain exactly as before ------------------- */
    kInputFile  = inFile;
    kOutputBase = outBase;
    gStyle->SetOptStat(0);

    ulog::banner("sPHENIX Run‑24 Au+Au QA – Enhanced Macro");

    std::unique_ptr<TFile> in = openInputFile(kInputFile);
    if (!in) return;

    CentList slices = discoverAndReportSlices(in.get());

    catalogueHistograms(in.get(), kOutputBase, kInputFile);

    QaMaps maps = runQaProduction(in.get(), kOutputBase, slices);

    printAllSummaries(in.get(), maps);
}

static std::vector<fs::path> discoverInputRuns()
{
    std::vector<fs::path> runFiles = listRunFiles(kInputDir);
    if (runFiles.empty())
        ulog::err("No input files found in " + kInputDir.string());
    return runFiles;                        // calling code checks empty()
}


static void selectTopNRuns(int                     nSample,
                           std::vector<fs::path>&  runFiles)
{
    if (nSample <= 0) return;               // nothing to do

    struct RunStat { fs::path file; long long nEvt; std::string run; };
    std::vector<RunStat> ranked;   ranked.reserve(runFiles.size());

    ulog::banner("Run pre‑scan – counting entries in h_vertexZ_<trigger>");

    for (const auto& f : runFiles)
    {
        std::smatch m;
        std::string fname = f.filename().string();
        std::regex_search(fname, m, std::regex(R"(output_([0-9]{8})\.root)"));
        std::string runStr = m.empty() ? "UNKNOWN" : m[1].str();

        std::unique_ptr<TFile> tf( TFile::Open(f.c_str(), "READ") );
        if (!tf || tf->IsZombie()) {
            ulog::warn("   ↳ " + runStr + "  – file unreadable, skipped");
            continue;
        }

        long long bestInFile = 0;
        TIter itTop( tf->GetListOfKeys() );
        while (auto* kDir = dynamic_cast<TKey*>(itTop())) {

            if (strcmp(kDir->GetClassName(), "TDirectoryFile")) continue;
            const std::string trg = kDir->GetName();
            if (!kTriggersWanted.count(trg)) continue;

            TDirectory* dTrig = static_cast<TDirectory*>(kDir->ReadObj());
            const std::string hName = "h_vertexZ_" + trg;

            if (auto* h = dynamic_cast<TH1*>( dTrig->Get(hName.c_str()) ))
                bestInFile = std::max(bestInFile,
                              static_cast<long long>(h->GetEntries()));
        }

        ulog::trace("   ↳ run " + runStr +
                   "   entries = " + std::to_string(bestInFile));

        if (bestInFile > 0) ranked.push_back( {f, bestInFile, runStr} );
    }

    if (ranked.empty()) {
        ulog::err("No non‑empty runs found – aborting sampling.");
        runFiles.clear();
        return;
    }

    std::sort(ranked.begin(), ranked.end(),
              [](const RunStat& a, const RunStat& b){ return a.nEvt > b.nEvt; });

    ulog::banner("Top‑run ranking (h_vertexZ entries)");
    std::cout << term::CLR_BOLD
              << std::left  << std::setw(6)  << "Rank"
              << std::setw(12) << "Run"
              << std::right << std::setw(15) << "Entries"
              << term::CLR_RST << "\n";

    int rankIdx = 1;
    for (const auto& r : ranked) {
        std::cout << std::left  << std::setw(6)  << rankIdx++
                  << std::setw(12) << r.run
                  << std::right << std::setw(15) << r.nEvt << "\n";
        if (rankIdx > nSample && nSample > 0) break;
    }

    if (static_cast<int>(ranked.size()) > nSample)
        ranked.resize(nSample);

    ulog::info("Selected " + std::to_string(ranked.size()) +
              " run(s) with the highest statistics for combined pass.");

    runFiles.clear();
    for (auto& r : ranked) runFiles.push_back( std::move(r.file) );
}


static void processRunsSequentially(const std::vector<fs::path>& runFiles,
                                    bool                         testRun)
{
    using fs::path;

    std::vector<fs::path> runs = runFiles;   // local copy (unaltered)

    if (testRun && runs.size() > 1) runs.resize(1);

    ROOT::TSequentialExecutor exec;
    ulog::banner("Running " + std::to_string(runs.size()) +
                " runs sequentially (shared memory)");

    auto worker = [&](unsigned int idx)->int
    {
        const path& f = runs[idx];

        std::smatch    m;
        const std::regex reRun(R"(output_([0-9]{8})\.root)");
        const std::string fname = f.filename().string();
        std::regex_search(fname, m, reRun);
        const std::string run = m.empty() ? "Single" : m[1].str();

        const std::string outBase = (kOutputDir / run).string();

        ulog::info(std::string("▶  (") +
                  (idx + 1 < 10 ? " " : "") + std::to_string(idx + 1) + "/" +
                  std::to_string(runs.size()) + ")  Run " + run + "  –  start");

        const auto t0 = std::chrono::steady_clock::now();
        runOneQaPass(f.string(), outBase);
        const auto dt = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t0).count();

        ulog::ok(std::string("✓  (") +
                (idx + 1 < 10 ? " " : "") + std::to_string(idx + 1) + "/" +
                std::to_string(runs.size()) + ")  Run " + run +
                "  –  done in " + std::to_string(dt).substr(0,5) + " s");
        return 0;
    };

    exec.Map(worker, ROOT::TSeqI(runs.size()));
}


static void mergeRunsAndReprocess(const std::vector<fs::path>& runFiles,
                                  bool                         testRun)
{
    using fs::path;
    if (testRun || runFiles.size() <= 1) return;

    /* ---------- dump HCal missing‑bin report -------------------- */
    HcalQA::writeMissingBinReport( kOutputDir / "MissingHCalBins.txt" );

    std::unordered_map<std::string, std::vector<std::string>> badRunMap;
    std::unordered_map<std::string, int>                      sebCount;

    {   /* 4 a. read MissingSEB.txt  +  bar‑chart output */
        const fs::path sebTxt  = baseDir / "output" / "MissingSEB.txt";
        const fs::path sebPng  = baseDir / "output" / "MissingSEB_distribution.png";

        std::ifstream miss( sebTxt.string() );
        if (!miss) {
            ulog::warn("MissingSEB.txt not found – merging all runs");
        } else {
            std::string line;
            while (std::getline(miss, line)) {
                if (line.empty()) continue;
                std::istringstream iss(line);
                std::string run;  iss >> run;
                std::string seb;
                while (iss >> seb) {
                    badRunMap[run].push_back(seb);
                    ++sebCount[seb];                 // accumulate counts
                }
            }
        }

        /* ---------- create bar‑chart (SEB00 … SEB15) ------------------ */
        TH1I hSEB("hMissingSEB",
                  "Runs with missing SEB;SEB index;Number of runs",
                  16, -0.5, 15.5);                   // 16 integer bins 0…15

        for (const auto& [seb,cnt] : sebCount) {
            int idx = -1;
            try { idx = std::stoi( seb.substr(3) ); }  // "SEB7" → 7
            catch (...) { continue; }
            if (idx >= 0 && idx < 16) hSEB.SetBinContent(idx + 1, cnt);
        }
        hSEB.SetFillColor(kAzure + 1);
        hSEB.SetBarWidth(0.8);
        hSEB.SetBarOffset(0.1);
        hSEB.GetXaxis()->SetTickLength(0);
        hSEB.LabelsOption("h");                       // horizontal x‑labels

        /* label the x‑axis bins “SEB00” … “SEB15” */
        for (int i = 1; i <= 16; ++i)
            hSEB.GetXaxis()->SetBinLabel(i,
                   Form("SEB%02d", i - 1));

        TCanvas cSEB("cMissingSEB","",800,500);
        gPad->SetGridy();
        hSEB.Draw("bar2");

        ensure_dir( sebPng.parent_path() );
        cSEB.SaveAs( sebPng.string().c_str() );

        ulog::ok("Missing SEB distribution plot saved → " + sebPng.string());
    }


    /* 4 b. filter list */
    std::vector<fs::path> mergeFiles;
    for (const auto& f : runFiles) {
        std::smatch m;
        const std::string fname = f.filename().string();
        if (std::regex_search(fname, m, std::regex(R"(output_([0-9]{8})\.root)"))) {
            const std::string run = m[1].str();
            if (badRunMap.count(run) == 0) mergeFiles.push_back(f);
        }
    }

    /* 4 c. terminal summary */
    ulog::banner("Missing SEB summary");
    std::size_t nBad = badRunMap.size(), nBad1 = 0, nBadMul = 0;
    for (const auto& [_,v] : badRunMap) (v.size()==1 ? ++nBad1 : ++nBadMul);

    std::cout << term::CLR_BOLD
              << "Runs with ≥1 missing SEB : " << nBad << "\n"
              << "   ├─ exactly one SEB    : " << nBad1 << "\n"
              << "   └─ multiple SEBs      : " << nBadMul << "\n"
              << term::CLR_RST << std::endl;

    if (!sebCount.empty()) {
        std::size_t w = 0;
        for (const auto& [seb,_] : sebCount) w = std::max(w, seb.size());
        w = std::max<std::size_t>(w, 5);

        std::cout << std::left << std::setw(w) << "SEB"
                  << " │ " << "Runs\n"
                  << std::string(w + 7, '-') << "\n";
        for (const auto& [seb,c] : sebCount)
            std::cout << std::left << std::setw(w) << seb
                      << " │ " << c << "\n";
        std::cout << std::string(w + 7, '=') << std::endl;
    }

    /* 4 d. hadd */
    if (mergeFiles.size() < 2) {
        ulog::warn("Skipping hadd – need ≥2 good runs, have "
                  + std::to_string(mergeFiles.size()));
        return;
    }

    const path combined = kInputDir / "output_ALL_COMBINED.root";
    ulog::banner("Hadd – building " + combined.string());

    ulog::info("Files to be merged (" + std::to_string(mergeFiles.size()) + " total):");
    std::uintmax_t totBytes = 0;
    for (const auto& f : mergeFiles) {
        const auto sz = fs::file_size(f);
        totBytes += sz;
        ulog::info("   + " + f.filename().string() +
                  "  (" + std::to_string(sz / 1'024'000) + " MB)");
    }
    ulog::info("   ------------------------------------------------");
    ulog::info("Accumulated input size : " +
              std::to_string(totBytes / 1'024'000) + " MB");

    const auto t0Hadd = std::chrono::steady_clock::now();

    TFileMerger merger(/*dryRun=*/false, /*verbose=*/true);
    merger.OutputFile(combined.c_str(), "RECREATE");
    for (const auto& f : mergeFiles) merger.AddFile(f.c_str());

    if (!merger.Merge()) {
        ulog::err("TFileMerger failed – combined QA skipped");
        return;
    }

    const auto dHadd   = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - t0Hadd).count();
    const auto outSize = fs::file_size(combined);

    ulog::ok("Combined ROOT file created in " +
            std::to_string(dHadd).substr(0,5) + " s,  size " +
            std::to_string(outSize / 1'024'000) + " MB");

    runOneQaPass(combined.string(),
                 (kOutputDir / "Combined").string());
}


void analyzeRun24or25auau(bool testRun = false, int nSample = -1)
{
    /* ── special mode: analyse existing combined file only ───────── */
    const bool combinedOnly = (std::getenv("COMBINED_ONLY") != nullptr);
    if (combinedOnly)
    {
        fs::path combined = kInputDir / "output_ALL_COMBINED.root";
        if (!fs::exists(combined)) {
            ulog::err("COMBINED_ONLY set but " + combined.string() + " not found");
            return;
        }
        ulog::banner("Combined‑only mode → " + combined.string());
        runOneQaPass(combined.string(), (kOutputDir / "Combined").string());
        return;                                // skip per‑run processing
    }

    /* 0. discover input ROOT files --------------------------------- */
    std::vector<fs::path> runFiles = discoverInputRuns();
    if (runFiles.empty()) return;

    /* 1. optional Top‑N sampling ----------------------------------- */
    selectTopNRuns(nSample, runFiles);
    if (runFiles.empty()) return;          // sampling aborted due to no stats

    /* 2. process every run sequentially */
    processRunsSequentially(runFiles, testRun);

    /* 3. optional merge + re‑run on combined file ------------------ */
    mergeRunsAndReprocess(runFiles, testRun);
}
