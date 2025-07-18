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
#include <TStyle.h>
#include <TH2Poly.h>
#include <TFileMerger.h>
#include <filesystem>
#include <TH3.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <cstdint>   // uintptr_t cast

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

std::array<double,8> fitPi0(TH1* h,bool& ok)
{
  double lo=h->GetXaxis()->GetXmin(), hi=h->GetXaxis()->GetXmax();
  TF1* f=new TF1("f","[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]*exp(-0.5*((x-[7])/[8])**2)",lo,hi);
  f->SetParameters(1,0,0,h->GetMaximum(),0.135,0.01,h->GetMaximum()/5.,0.55,0.02);
  f->SetParLimits(4,0.11,0.16); f->SetParLimits(5,0.005,0.05);
  f->SetParLimits(7,0.45,0.70); f->SetParLimits(8,0.01,0.06);
  ok=(h->Fit(f,"QNRS")==0);
  return {f->GetParameter(4),f->GetParError(4),
          f->GetParameter(5),f->GetParError(5),
          f->GetParameter(7),f->GetParError(7),
          f->GetParameter(8),f->GetParError(8)};
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
// ║ 8.  QA MODULES                              ║
// ╚══════════════════════════════════════════════╝
class Pi0QA : public QA
{
public:
  /* pT bins that appear in the file names – must match CutKey values   */
  const std::vector<std::pair<float,float>> m_ptBins {
      {2,4},{4,6},{6,8},{8,10},{10,12},{12,15},{15,20},{20,30} };

  struct FitInfo {
    std::string slice;  double pLo, pHi;         // bin identifiers
    double mean, sigma;                          // Gaussian parameters
    double chi2; int ndf;                        // fit quality
  };

  Pi0QA(std::string t, fs::path b, const CentList& s, std::ofstream& csv)
      : QA(std::move(t), std::move(b), s), csv(csv) {}

  // --------------------------------------------------------------------
  // D E S T R U C T O R   – creates the multi‑panel summary canvases
  // --------------------------------------------------------------------
  ~Pi0QA() override
  {
    // -------- 2×3 canvas: one pad per centrality (pT‑integrated) ------
    if (!_centralHists.empty()) {
      TCanvas cCent("c_pi0Cent", "Pi0 mass – all centralities", 1800, 1000);
      cCent.Divide(3, 2, 0.01, 0.01);

      int pad = 1;
      for (const auto& sl : slices) {
        auto it = _centralHists.find(sl);
        if (it == _centralHists.end()) continue;

        cCent.cd(pad++);
        it->second->SetStats(0);
        it->second->Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
        std::string centTxt = (sl == "Inclusive")
                                  ? "Inclusive"
                                  : ("Centrality: " +
                                     sl.substr(0, sl.find('_')) + " to " +
                                     sl.substr(sl.find('_') + 1) + " %");
        tl.DrawLatex(0.12, 0.85, centTxt.c_str());
      }
      fs::path out = root / "EMCal/pi0QA" / "Pi0Mass_AllCentrality.png";
      ensure_dir(out.parent_path());
      cCent.SaveAs(out.string().c_str());
    }

    // -------- 2×3 canvas per centrality: one pad per pT bin ----------
    for (const auto& [slice, vec] : _ptHists) {
      const int padsPerPage = 6;
      int page = 0, padInPage = 0;
      std::unique_ptr<TCanvas> cPT;

      auto newCanvas = [&](int pg){
        std::string cname = "c_pi0PT_" + slice + "_p" + std::to_string(pg);
        cPT.reset(new TCanvas(cname.c_str(), cname.c_str(), 1800, 1000));
        cPT->Divide(3, 2, 0.01, 0.01);
        padInPage = 0;
      };

      newCanvas(page);

      for (size_t i = 0; i < vec.size(); ++i) {
        if (!vec[i]) continue;                    // pT bin absent

        if (padInPage == padsPerPage) {          // new page
          fs::path out = root / "EMCal/pi0QA" /
                         ("Pi0Mass_pT_" + slice +
                          "_page" + std::to_string(page) + ".png");
          ensure_dir(out.parent_path());
          cPT->SaveAs(out.string().c_str());
          ++page;
          newCanvas(page);
        }

        cPT->cd(++padInPage);
        vec[i]->SetStats(0);
        vec[i]->Draw();

        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
        const auto& bin = m_ptBins[i];
        tl.DrawLatex(0.10, 0.85,
                       Form("%.0f #leq p_{T}^{M#gamma#gamma} < %.0f GeV",
                            bin.first, bin.second));
      }

      /* write last (or single) canvas */
      if (padInPage) {
        fs::path out = root / "EMCal/pi0QA" /
                       ("Pi0Mass_pT_" + slice +
                        "_page" + std::to_string(page) + ".png");
        ensure_dir(out.parent_path());
        cPT->SaveAs(out.string().c_str());
      }
    }
  }

  const auto& fitSummary() const { return _fitSummary; }

  // --------------------------------------------------------------------
  // P R O C E S S   – called for every TObject in the ROOT file
  // --------------------------------------------------------------------
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;
    std::string n = o->GetName();
    if (n.rfind("mInv_", 0) != 0) return false;          // not π0 spectrum

    CutKey ck;
    if (!decodeInvName(n, ck)) return false;

    const std::string slice = sliceKey(n);
    const bool isPtInt      = (ck.pLo < 0 || ck.pHi < 0);

    // ----- build output sub‑directory path identical to before --------
    fs::path sub = "EMCal/pi0QA";
    sub /= ("E"    + sf3(ck.E)   +
            "_Chi" + sf3(ck.chi) +
            "_Asym" + sf3(ck.asy));
    if (!isPtInt)
      sub /= ("pT_" + sf3(ck.pLo) + "_to_" + sf3(ck.pHi));

    fs::path outPng = cPath(root, slice, sub) / (n + ".png");
    ensure_dir(outPng.parent_path());

    // -------------------- perform the single‑Gaussian fit -------------
    TH1* h = static_cast<TH1*>(o);
    const double lo = 0.08, hi = 0.20;                // fixed fit region

    TF1 fitFun("fPi0", "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])**2)",
               lo, hi);
    fitFun.SetParameters(1, 0, 0, h->GetMaximum(), 0.135, 0.01);
    fitFun.SetParLimits(4, 0.12, 0.14);               // μ window
    fitFun.SetParLimits(5, 0.005, 0.03);              // σ upper bound

    bool fitOK = (h->Fit(&fitFun, "QNRS") == 0);

    TF1 polyBg("polyBg", "[0]+[1]*x+[2]*x*x", lo, hi);
    polyBg.SetParameters(fitFun.GetParameter(0),
                         fitFun.GetParameter(1),
                         fitFun.GetParameter(2));
    polyBg.SetLineColor(kBlue+1); polyBg.SetLineWidth(2);

    TF1 gausSig("gausSig", "[0]*exp(-0.5*((x-[1])/[2])**2)", lo, hi);
    gausSig.SetParameters(fitFun.GetParameter(3),
                            fitFun.GetParameter(4),
                            fitFun.GetParameter(5));
    gausSig.SetLineColor(kRed);   gausSig.SetLineWidth(2);

    // ----------------------------- draw & save ------------------------
    TCanvas c1;
    h->SetStats(0);
    h->Draw();
    polyBg.Draw("SAME"); gausSig.Draw("SAME");

    TLegend leg(0.55, 0.70, 0.88, 0.88); leg.SetBorderSize(0);
    leg.AddEntry(&gausSig, "Gaussian signal",      "l");   // red
    leg.AddEntry(&polyBg,  "Polynomial background","l");   // blue
    leg.Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.04);
    tl.DrawLatex(0.55, 0.63,
                   Form("Gaussian: %.4f #pm %.4f GeV",
                        gausSig.GetParameter(1),
                        gausSig.GetParameter(2)));
    c1.SaveAs(outPng.string().c_str());

    // ----------------------------- CSV + map --------------------------
    if (fitOK) {
      csv << trig  << ',' << ck.E   << ',' << ck.chi << ',' << ck.asy << ','
          << ck.pLo << ',' << ck.pHi << ','
          << gausSig.GetParameter(1) << ','         // mean
          << gausSig.GetParameter(2) << ','         // sigma
          << fitFun.GetChisquare()   << ','         // χ²
          << fitFun.GetNDF()         << '\n';       // ndf

      _fitSummary.emplace(n, FitInfo{slice, ck.pLo, ck.pHi,
                                     gausSig.GetParameter(1),
                                     gausSig.GetParameter(2),
                                     fitFun.GetChisquare(),
                                     fitFun.GetNDF()});
    }

    // ------ cache histograms for overview canvases --------------------
    if (isPtInt) {
      auto* cl = static_cast<TH1*>(h->Clone(("__cl_" + n).c_str()));
      cl->SetDirectory(nullptr);
      _centralHists.emplace(slice, cl);
    } else {
      // find the matching pT bin index
      auto it = std::find_if(m_ptBins.begin(), m_ptBins.end(),
                             [&](auto& p){ return fabs(p.first-ck.pLo)<1e-3
                                                && fabs(p.second-ck.pHi)<1e-3; });
      if (it != m_ptBins.end()) {
        size_t idx = std::distance(m_ptBins.begin(), it);
        auto& vec  = _ptHists[slice];
        if (vec.size() < m_ptBins.size()) vec.resize(m_ptBins.size(), nullptr);
        auto* cl = static_cast<TH1*>(h->Clone(("__cl_" + n).c_str()));
        cl->SetDirectory(nullptr);
        vec[idx] = cl;
      }
    }
    return true;
  }

private:
  std::ofstream&                                   csv;
  std::unordered_map<std::string, TH1*>            _centralHists;  // pT‑int
  std::unordered_map<std::string,
                     std::vector<TH1*>>            _ptHists;       // by pT
  std::unordered_map<std::string, FitInfo>         _fitSummary;    // all fits
};


// ——— Detector–detector correlations ————————————
class CorrQA : public QA{
public: using QA::QA;
  bool process(TObject* o) override
  {
    if(!o->InheritsFrom(TH2::Class())) return false;
    string n=o->GetName(); if(n.find("_vs_")==string::npos) return false;
    string sl=sliceKey(n);

    auto save=[&](const string& slice){
        fs::path out = cPath(root, slice, "Correlations") / (n + ".png");

        // make a local, non‑owned pointer for clarity
        TH2* h2 = static_cast<TH2*>(o);
        h2->SetStats(0);

        // ── tighten   X‑range ───────────────────────────────────
        int firstX = 1, lastX = h2->GetNbinsX();
        while (firstX <= lastX &&
               h2->Integral(firstX, firstX, 1, h2->GetNbinsY()) == 0) ++firstX;
        while (lastX  >= firstX &&
               h2->Integral(lastX,  lastX,  1, h2->GetNbinsY()) == 0) --lastX;
        if (firstX < lastX) h2->GetXaxis()->SetRange(firstX, lastX);

        // ── tighten   Y‑range ───────────────────────────────────
        int firstY = 1, lastY = h2->GetNbinsY();
        while (firstY <= lastY &&
               h2->Integral(1, h2->GetNbinsX(), firstY, firstY) == 0) ++firstY;
        while (lastY  >= firstY &&
               h2->Integral(1, h2->GetNbinsX(), lastY,  lastY ) == 0) --lastY;
        if (firstY < lastY) h2->GetYaxis()->SetRange(firstY, lastY);

        // ── draw with log‑Z scale ───────────────────────────────
        TCanvas c("c", "", 1100, 800);
        c.SetLogz();
        h2->Draw("COLZ");

        ensure_dir(out.parent_path());
        c.SaveAs(out.string().c_str());
    };

    save(sl);
    return true;
  }
};

// ─── EMCal QA – η‑φ maps with bad‑board masking + 2×3 overview ───────
class EmcalQA : public QA
{
public:
  using QA::QA;

  /* summary panel shown at the end of the job */
  ~EmcalQA() override
  {
    if (_centralMaps.empty()) return;

    TCanvas c("c_emcalCent", "EMCal hit‑maps – all centralities", 2100, 1200);
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

  // ------------------------------------------------------------------
  //  P R O C E S S   – called for every histogram in the EMCal folder
  // ------------------------------------------------------------------
  bool process(TObject* o) override
  {
    if (!o->InheritsFrom(TH1::Class())) return false;
    const std::string n = o->GetName();
    if (n.rfind("h_EMC_", 0) != 0)      return false;

    const std::string slice   = sliceKey(n);
    const bool        isMap   = (n.find("_EtaPhiMap_") != std::string::npos);

    auto makePanel = [&](TH2* src) -> std::unique_ptr<TH2F>
    {
        /* 0)  constants -------------------------------------------------- */
        constexpr int nPhi = 256;                // rows  (Y)
        constexpr int nEta =  96;                // cols  (X)
        constexpr int px   =   6;                // pixel‑size in the PNG

        /* 1)  local clone + bad‑board masking ---------------------------- */
        std::unique_ptr<TH2> h(static_cast<TH2*>(src->Clone()));
        h->SetDirectory(nullptr);   h->SetStats(0);   h->SetContour(100);

        for (int ip = 0; ip < nPhi; ++ip)
          for (int ie = 0; ie < nEta; ++ie)
            if (isBadBoard(sector_from_idx(ie, ip),
                           ib_from_idx(ie, ip)))
              h->SetBinContent(h->FindBin(ip, ie), -9999.);   // hole = white

        /* 2)  rotate  (η → X,  φ → Y) ----------------------------------- */
        std::unique_ptr<TH2F> rot(new TH2F(("hRot_"+std::string(src->GetName())).c_str(),
                                           h->GetTitle(),
                                           nEta, 0, nEta,       // X
                                           nPhi, 0, nPhi));     // Y
        for (int ip = 1; ip <= nPhi; ++ip)
          for (int ie = 1; ie <= nEta; ++ie)
            rot->SetBinContent(ie, ip, h->GetBinContent(ip, ie));

        rot->SetDirectory(nullptr);
        rot->SetMinimum(1.);                         // under‑flow = white
        rot->SetTitleOffset(0.9,"X"); rot->SetTitleOffset(0.9,"Y");
        rot->GetXaxis()->SetTitle("Tower #eta");
        rot->GetYaxis()->SetTitle("Tower #phi");
        rot->GetXaxis()->SetNdivisions(12,kFALSE);   // every 8 η
        rot->GetYaxis()->SetNdivisions(32,kFALSE);   // every 8 φ

        /* 3)  square‑pixel canvas --------------------------------------- */
        const int cw = nEta * px;                    // 96⋅6 = 576 px
        const int ch = nPhi * px;                    // 256⋅6 = 1536 px

        fs::path outPng = cPath(root, slice, "EMCal") / (src->GetName() + std::string(".png"));
        ensure_dir(outPng.parent_path());

        TCanvas c("c_emcal","",cw,ch);
        c.SetRightMargin(0.17);
        c.SetLeftMargin (0.08);
        c.SetBottomMargin(0.08);
        c.SetTopMargin  (0.04);
        c.SetFixedAspectRatio();                     // 1 bin ⇔ 1 pixel

        rot->Draw("COLZ");

        /* 4)  ultra‑light grid (every 8 towers) -------------------------- */
        TLine l;
        l.SetLineColor(kBlack); l.SetLineWidth(1);
        for (int x = 0; x <= nEta; x += 8) { l.DrawLine(x, 0, x, nPhi); }
        for (int y = 0; y <= nPhi; y += 8) { l.DrawLine(0, y, nEta, y); }

        /* 5)  North / South divider (η = 48) ----------------------------- */
        TLine ns(48, 0, 48, nPhi); ns.SetLineColor(kBlack); ns.SetLineWidth(3);
        ns.Draw();

        /* 6)  sector / inner‑board labels (subtle grey) ------------------ */
        TLatex tx; tx.SetTextSize(0.020); tx.SetTextAlign(22); tx.SetTextColorAlpha(kGray+2,0.7);

        for (int s = 0; s < 64; ++s) {                      // 32 per hemisphere
          const int basePhi = (s % 32) * 8;
          const double yMid = basePhi + 4;
          const double xSec = (s < 32) ? 72 : 24;
          tx.DrawLatex(xSec, yMid, Form("S%d", s));
        }

        c.SaveAs(outPng.string().c_str());
        return rot;
    };

      
    if (isMap && o->InheritsFrom(TH2::Class())) {
      auto p = makePanel(static_cast<TH2*>(o));
      if (_centralMaps.find(slice) == _centralMaps.end())
        _centralMaps.emplace(slice, std::move(p));

    } else if (o->InheritsFrom(TH2::Class())) {
      fs::path out = cPath(root, slice, "EMCal") / (n + ".png");
      save2D(static_cast<TH2*>(o), out);

    } else {                                            // plain 1‑D
      fs::path out = cPath(root, slice, "EMCal") / (n + ".png");
      save1D(static_cast<TH1*>(o), out);
    }
    return true;
  }

private:
  std::unordered_map<std::string,
                     std::unique_ptr<TH2F>> _centralMaps;
};


// ——— HCal QA ————————————————————————————————————
class HcalQA : public QA{
public: using QA::QA;
  bool process(TObject* o) override
  {
    if(!o->InheritsFrom(TH1::Class())) return false;
    string n=o->GetName();
    bool isI=n.rfind("h_IHCAL_",0)==0, isO=n.rfind("h_OHCAL_",0)==0;
    if(!isI && !isO) return false;
    string sl=sliceKey(n); bool etaPhi=n.find("_EtaPhiMap_")!=string::npos;

    auto save=[&](const string& slice){
      fs::path sub=isI?"IHCal":"OHCal";
      fs::path out=cPath(root,slice,sub)/(n+".png");
      if(etaPhi && o->InheritsFrom(TH2::Class()))
      {
        auto* h=static_cast<TH2*>(o);
        std::unique_ptr<TH2> h2(static_cast<TH2*>(h->Clone()));
        h2->SetDirectory(nullptr); h2->SetStats(0); h2->SetContour(99); h2->SetMinimum(1.);
        for(int iphi=0;iphi<64;++iphi)
          for(int ieta=0;ieta<24;++ieta)
            if(isBadHcalPlate(hcal_sector_from_idx(ieta,iphi),
                              hcal_plate_from_idx(ieta,iphi)))
              h2->SetBinContent(h2->FindBin(iphi,ieta),-9999.);
          // replacement – identical hit‑map but with the sector/board/tower grid
          TCanvas c("c","IHCal/OHCal hit‑map",1300,800);
          c.SetRightMargin(0.16);
          c.SetBottomMargin(0.07);
          c.SetLeftMargin(0.08);

          /* axis labelling identical to the online monitor */
          h2->GetXaxis()->SetTitle("#eta index");
          h2->GetYaxis()->SetTitle("#phi index");
          h2->GetXaxis()->SetNdivisions(24, kFALSE);   // 24 η bins
          h2->GetYaxis()->SetNdivisions(64, kFALSE);   // 64 φ bins
          h2->Draw("COLZ");

          /* ── thick horizontal lines: 32 sectors (every 2 φ towers) ── */
          for (int s = 0; s < 32; ++s) {
            TLine *ls = new TLine(0, (s+1)*2, 24, (s+1)*2);
            ls->SetLineColor(kBlack); ls->SetLineWidth(4); ls->Draw();
          }

          /* ── thick vertical lines at η = 8 and 16 (inner‑board groups) ── */
          TLine lB1(8,  0, 8, 64);  lB1.SetLineColor(kBlack); lB1.SetLineWidth(4); lB1.Draw();
          TLine lB2(16, 0, 16, 64); lB2.SetLineColor(kBlack); lB2.SetLineWidth(4); lB2.Draw();

          /* ── fine grid: individual towers (1 × 1) ── */
          for (int iphi = 0; iphi < 64; ++iphi) {
            TLine *ly = new TLine(0, iphi+1, 24, iphi+1);
            ly->SetLineColor(kBlack); ly->SetLineWidth(1); ly->Draw();
          }
          for (int ieta = 0; ieta < 24; ++ieta) {
            TLine *lx = new TLine(ieta+1, 0, ieta+1, 64);
            lx->SetLineColor(kBlack); lx->SetLineWidth(1); lx->Draw();
          }

          ensure_dir(out.parent_path());
          c.SaveAs(out.string().c_str());
      }
      else if(o->InheritsFrom(TH2::Class())) save2D(static_cast<TH2*>(o),out);
      else                                   save1D(static_cast<TH1*>(o),out);
    };
    save(sl);
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

    auto save=[&](const string& slice){
      fs::path outPng=cPath(root,slice,DERIVED::subdir)/(DERIVED::fileName(trig)+".png");
      TCanvas c("c","",1100,600); c.Divide(2,1,0.01,0.01);
      c.cd(1); out.s->SetTitle(DERIVED::titleSouth); out.s->Draw("COL POLZ");
      c.cd(2); out.n->SetTitle(DERIVED::titleNorth); out.n->Draw("COL POLZ");
      ensure_dir(outPng.parent_path()); c.SaveAs(outPng.string().c_str());
    };
    save(sl);
    return true;
  }
protected:
  NSCache<MapPair>& cache;
};

                                                                 
struct MBDTag{
  // accept only true‑MBD hit‑maps – histogram name must START with “h_MBD_”
  static bool accept(const std::string& s)
  {
    return s.rfind("h_MBD_", 0) == 0;   // good MBD names: h_MBD_…
  }
  static constexpr const char* subdir = "MBD";
  static std::string fileName(const std::string& t){ return "MBD_Hitmap_NS_" + t; }
  static constexpr const char* titleSouth = "MBD South";
  static constexpr const char* titleNorth = "MBD North";
};

struct sEPDTag{
  static bool accept(const string& s){return s.find("sEPD")!=string::npos;}
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

  static void saveProjection(TH3* h3, const char* axes,
                             const fs::path& png)
  {
    std::unique_ptr<TH2> h2(h3->Project3D(axes));  // ROOT owns → copy
    h2->SetDirectory(nullptr); h2->SetStats(0);
    h2->SetTitle((std::string(h3->GetTitle())+"  –  "+axes).c_str());
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


// ╔══════════════════════════════════════════════╗
// ║ 9.  MAIN DRIVER  – run‑by‑run + combined     ║
// ╚══════════════════════════════════════════════╝
void analyzeRun24or25auau()
{
  gStyle->SetOptStat(0);          // ROOT style once is enough
  gStyle->SetPalette(kViridis);
  // ------------------------------------------------------------
  // 0.  Collect all “output_########.root” files in kInputDir
  // ------------------------------------------------------------
  std::vector<fs::path> runFiles = listRunFiles(kInputDir);

  /*  If the directory is empty we fall back to the original
   *  single‑file behaviour: whatever kInputFile already points to
   *  (this lets you keep testing on one file if you want).       */
  if (runFiles.empty()) {
    if (kInputFile.empty()) {
      log::err("No input files found and kInputFile is empty – nothing to do.");
      return;
    }
    runFiles.push_back(kInputFile);          // “single‑file” mode
  }

  // ------------------------------------------------------------
  // 1.  Lambda = ONE COMPLETE QA PASS  (the old body unchanged)
  // ------------------------------------------------------------
  auto runOneQaPass = [&]()
  {
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

    // ─── Pass‑0 catalogue ──────────────────────────────────────
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

    // ─── Pass‑1 QA ─────────────────────────────────────────────
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

      // directory skeleton for this trigger
      for (auto& s : slices) {
        fs::path b = fs::path(kOutputBase) / trg;
        if (s != "Inclusive") b /= ("Cent_" + s);
        for (auto sub : { "Correlations", "EMCal/pi0QA", "IHCal", "OHCal",
                          "MBD", "sEPD", "jetQA" })
          ensure_dir(b / sub);
      }

      // QA modules
      std::vector<std::unique_ptr<QA>> qa;
      fs::path base = fs::path(kOutputBase) / trg;
      qa.emplace_back(std::make_unique<Pi0QA >(trg,base,slices,csv));
      qa.emplace_back(std::make_unique<CorrQA>(trg,base,slices));
      qa.emplace_back(std::make_unique<EmcalQA>(trg,base,slices));
      qa.emplace_back(std::make_unique<HcalQA >(trg,base,slices));
      qa.emplace_back(std::make_unique<MbdQA  >(trg,base,slices,mbdCache));
      qa.emplace_back(std::make_unique<SepdQA >(trg,base,slices,sepdCache));
      qa.emplace_back(std::make_unique<JetQA >(trg,base,slices));

      // histogram loop
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

    // ─── Summary ───────────────────────────────────────────────
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
  };   // <‑‑ end lambda runOneQaPass
  // ------------------------------------------------------------
  // 2.  Run the QA pass for every individual run file
  // ------------------------------------------------------------
  std::regex reRun(R"(output_([0-9]{8})\.root)");
  for (const auto& f : runFiles) {
    std::smatch m;
    std::string fname = f.filename().string();
    std::regex_search(fname.cbegin(), fname.cend(), m, reRun);
    const std::string run = m.empty() ? "Single" : m[1].str();

    kInputFile  = f.string();
    kOutputBase = (kOutputDir / run).string();

    log::info("▶  Run " + run + "  →  " + kOutputBase);
    runOneQaPass();
  }

  // ------------------------------------------------------------
  // 3.  Merge (“hadd”) all run files and do one final pass
  // ------------------------------------------------------------
  if (runFiles.size() > 1) {
    fs::path combined = kInputDir / "output_ALL_COMBINED.root";
    log::banner("Hadd – building " + combined.string());

    TFileMerger merger(false, true);           // print summary
    merger.OutputFile(combined.c_str(), "RECREATE");
    for (const auto& f : runFiles) merger.AddFile(f.string().c_str());

    if (!merger.Merge()) {
      log::err("TFileMerger failed – combined QA skipped");
      return;
    }
    log::ok("Combined ROOT file written");

    kInputFile  = combined.string();
    kOutputBase = (kOutputDir / "Combined").string();

    log::banner("Final QA pass on combined file");
    runOneQaPass();
  }
}
