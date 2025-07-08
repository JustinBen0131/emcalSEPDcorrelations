// analyzeRun24auau.cpp   (ROOT ≥ 6, c++17)
// ======================================================================
//  * catalogue every histogram in the file (prints + text‑file)
//  * then run QA / invariant‑mass analysis in modular detector classes
//  * one output tree per trigger as described above
// ======================================================================

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH2Poly.h>
#include <unordered_map>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <map>
#include <vector>
#include <tuple>
#include <array>
#include <memory>
#include <iomanip>
#include <exception>

using std::string;
namespace fs = std::filesystem;

// ----------------------------------------------------------------------
//               ----------  USER CONFIGURATION  ----------
// ----------------------------------------------------------------------
const string kInputFile =
    "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/"
    "emcal_sepd_analysis_run54128_c0_DST_CALO_run2auau_new_2024p007-00054128-00000.root";

const string kOutputBase =
    "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/output";

/* <<<  LIST THE TRIGGERS YOU WANT TO PROCESS HERE  >>>               */
std::set<string> kTriggersWanted { "MBD_NandS_geq_2" };   // <-- edit here if needed
// ----------------------------------------------------------------------

std::map<string,int> kTriggerColours {
  {"MBD_NandS_geq_2", kBlack},
  {"MBD_NandS_geq_1", kBlue+1},
  {"Photon_10_GeV"  , kRed+1},
  {"Photon_15_GeV"  , kGreen+2}
};
std::map<string,string> kPrettyTriggerName {
  {"MBD_NandS_geq_2", "MBD ≥ 2"},
  {"MBD_NandS_geq_1", "MBD ≥ 1"},
  {"Photon_10_GeV"  , "γ‑10"},
  {"Photon_15_GeV"  , "γ‑15"}
};

// ----------------------------------------------------------------------
//                         LOGGING HELPERS
// ----------------------------------------------------------------------
namespace log {
  constexpr const char* kClrReset  = "\033[0m";
  constexpr const char* kClrRed    = "\033[31m";
  constexpr const char* kClrGreen  = "\033[32m";
  constexpr const char* kClrYellow = "\033[33m";
  constexpr const char* kClrCyan   = "\033[36m";

  inline void info(const string& m)
  { std::cout<<kClrCyan  <<"[INFO] "<<m<<kClrReset<<"\n"; }

  inline void ok  (const string& m)
  { std::cout<<kClrGreen <<"[OK]   "<<m<<kClrReset<<"\n"; }

  inline void warn(const string& m)
  { std::cout<<kClrYellow<<"[WARN] "<<m<<kClrReset<<"\n"; }

  inline void err (const string& m)
  { std::cerr<<kClrRed   <<"[ERR]  "<<m<<kClrReset<<"\n"; }
}

// ----------------------------------------------------------------------
//                     GENERIC UTILITY FUNCTIONS
// ----------------------------------------------------------------------
inline void ensure_dir(const fs::path& p)
{
  std::error_code ec;
  fs::create_directories(p, ec);
}

inline string sf3(double x)
{
  std::ostringstream o;
  if(x==0.) { o<<"0"; return o.str(); }
  int e = int(std::floor(std::log10(std::fabs(x))));
  o<<std::fixed<<std::setprecision(std::max(0,2-e))<<x;
  string s=o.str(); std::replace(s.begin(),s.end(),'.','p'); return s;
}

// ------------ EMCal geometry helpers (unchanged) ----------------------
static inline int sector_from_idx(unsigned ieta,unsigned iphi)
{ if(iphi>=256) return -1; int base=iphi/8; return (ieta<48)?32+base:base; }

static inline int ib_from_idx(unsigned ieta,unsigned)
{
  if(ieta< 8) return 5; else if(ieta<16) return 4; else if(ieta<24) return 3;
  else if(ieta<32) return 2; else if(ieta<40) return 1; else if(ieta<48) return 0;
  else if(ieta<56) return 0; else if(ieta<64) return 1; else if(ieta<72) return 2;
  else if(ieta<80) return 3; else if(ieta<88) return 4; else if(ieta<96) return 5;
  return -1;
}

// ----------------------------------------------------------------------
//                    BAD‑INTERFACE‑BOARD  (CEMC)
// ----------------------------------------------------------------------
inline bool isBadBoard(int sector, int ib)
{
  return ((sector==50 && ib==1) ||
          (sector==4  && ib==1) ||
          (sector==25 && ib==2));
}

// ---------- invariant‑mass helpers (same as before) -------------------
struct CutKey {
  float E,chi,asy,pLo,pHi; string trigger;
  string to_string() const
  {
    std::ostringstream o; o<<"E"<<sf3(E)<<"_Chi"<<sf3(chi)<<"_Asym"<<sf3(asy);
    if(pLo>=0&&pHi>=0) o<<"/pT_"<<sf3(pLo)<<"_to_"<<sf3(pHi); return o.str();
  }
};

bool decodeInvName(const string& n, CutKey& o)
{
  std::regex r(R"(mInv_pt([\-0-9\.]+)to([\-0-9\.]+)_E([0-9\.]+)_chi([0-9\.]+)_asy([0-9\.]+)_(.+))");
  std::smatch m; if(!std::regex_match(n,m,r)) return false;
  o.pLo=std::stof(m[1]); o.pHi=std::stof(m[2]); o.E =std::stof(m[3]);
  o.chi=std::stof(m[4]); o.asy=std::stof(m[5]); o.trigger=m[6]; return true;
}

std::tuple<bool,std::array<double,8>> fitInvariantMass(TH1* h)
{
  double lo = h->GetXaxis()->GetXmin();
  double hi = h->GetXaxis()->GetXmax();

  TF1* f = new TF1("fTotal",
        "[0]+[1]*x+[2]*x*x+[3]*exp(-0.5*((x-[4])/[5])**2)+[6]*exp(-0.5*((x-[7])/[8])**2)",
        lo, hi);

  f->SetParameters(1,0,0,
                   h->GetMaximum(),0.135,0.01,
                   h->GetMaximum()/5.,0.55,0.02);

  f->SetParLimits(5,0.005,0.05);
  f->SetParLimits(8,0.01,0.06);
  f->SetParLimits(4,0.11,0.16);
  f->SetParLimits(7,0.45,0.70);

  int status = h->Fit(f,"QNR");
  if(status!=0) log::warn("Fit did not fully converge (status "+std::to_string(status)+")");

  std::array<double,8> p{{f->GetParameter(4),f->GetParError(4),
                          f->GetParameter(5),f->GetParError(5),
                          f->GetParameter(7),f->GetParError(7),
                          f->GetParameter(8),f->GetParError(8)}};
  return {true,p};
}

// -------- generic drawing helpers  ------------------------------------
void saveEtaPhiMapWithGrid(TH2* h,const fs::path& out)
{
  TCanvas c("c","",900,500); h->SetStats(0); h->Draw("COLZ");
  for(int iphi=0; iphi<=256; iphi+=8){ TLine l(iphi,0,iphi,96); l.SetLineStyle(3); l.SetLineColor(kWhite); l.Draw(); }
  for(int ieta=0; ieta<=96; ieta+=8){  TLine l(0,ieta,256,ieta); l.SetLineStyle(3); l.SetLineColor(kWhite); l.Draw(); }
  TLine eq(0,48,256,48); eq.SetLineColor(kBlack); eq.SetLineWidth(2); eq.Draw();
  ensure_dir(out.parent_path()); c.SaveAs(out.string().c_str());
}
void saveHist1D(TH1* h,const fs::path& out){ TCanvas c; h->SetStats(0); h->Draw(); ensure_dir(out.parent_path()); c.SaveAs(out.string().c_str()); }
void saveHist2D(TH2* h,const fs::path& out){ TCanvas c; h->SetStats(0); h->Draw("COLZ"); ensure_dir(out.parent_path()); c.SaveAs(out.string().c_str()); }

// ======================================================================
//                         QA CLASS HIERARCHY
// ======================================================================
class TriggerQA {
public:
  TriggerQA(const string& trg, const fs::path& out) :
    trigger(trg), outDir(out) {}
  virtual ~TriggerQA() = default;
  virtual bool process(TObject* obj) = 0;

protected:
  string   trigger;
  fs::path outDir;
};

// --------------------------- π0 / correlations -------------------------
class Pi0QA : public TriggerQA {
public:
  Pi0QA(const string& trg, const fs::path& out, std::ofstream& csvStream) :
      TriggerQA(trg,out), csv(csvStream) {}

  bool process(TObject* obj) override
  {
    string hname = obj->GetName();
    if(hname.rfind("mInv_",0)!=0 || !obj->InheritsFrom(TH1::Class())) return false;

    CutKey ck; if(!decodeInvName(hname,ck)) return false;

    fs::path dirPath = outDir /
                       ("E"+sf3(ck.E)+"_Chi"+sf3(ck.chi)+"_Asym"+sf3(ck.asy));
    if(ck.pLo>=0&&ck.pHi>=0) dirPath /= ("pT_"+sf3(ck.pLo)+"_to_"+sf3(ck.pHi));
    ensure_dir(dirPath);

    TH1* h = static_cast<TH1*>(obj);
    auto [ok,pars] = fitInvariantMass(h);

    TCanvas c; h->SetStats(0); h->Draw();
    if(ok){
      TF1* f = h->GetFunction("fTotal"); f->SetLineColor(kRed); f->Draw("SAME");
      csv << ck.trigger << "," << ck.E << "," << ck.chi << "," << ck.asy << ","
          << ck.pLo << "," << ck.pHi << ","
          << pars[0] << "," << pars[1] << "," << pars[2] << "," << pars[3] << ","
          << pars[4] << "," << pars[5] << "," << pars[6] << "," << pars[7] << "\n";
    }
    fs::path pngPath = dirPath / (hname + ".png");
    c.SaveAs(pngPath.string().c_str());
    log::info("Correlations saved " + pngPath.string());
    return true;
  }
private:
  std::ofstream& csv;
};

// ---------------------------  EMCal QA  --------------------------------
class EmcalQA : public TriggerQA {
public:
  using TriggerQA::TriggerQA;

private:
  void drawEtaPhiWithMask(TH2* h, const fs::path& outPNG)
  {
    std::unique_ptr<TH2> h2(dynamic_cast<TH2*>(h->Clone(Form("%s_clone",h->GetName()))));
    h2->SetDirectory(nullptr);

    // --- zero out bad boards
    for(int phi=0; phi<256; ++phi){
      for(int eta=0; eta<96; ++eta){
        int sector=sector_from_idx(eta,phi);
        int ib    =ib_from_idx(eta,phi);
        if(isBadBoard(sector,ib))
          h2->SetBinContent(h2->GetXaxis()->FindBin(phi),
                            h2->GetYaxis()->FindBin(eta), -9999.);
      }
    }
    h2->SetContour(99); h2->SetMinimum(1.0);

    gStyle->SetOptStat(0);
    TCanvas c("cCEMC","CEMC η‑φ",1600,1200);
    c.SetRightMargin(0.15);
    h2->SetTitle(Form("%s – bad IBs masked", h->GetTitle()));
    h2->GetXaxis()->SetTitle("Tower ϕ index");
    h2->GetYaxis()->SetTitle("Tower η index");
    h2->Draw("COLZ");

    for(int i=0;i<=256;i+=8){ TLine l(i,0,i,96); l.Draw(); }
    for(int j=0;j<=96; j+=8){ TLine l(0,j,256,j); l.Draw(); }

    TLatex latSec; latSec.SetTextSize(0.018); latSec.SetTextAlign(22); latSec.SetTextColor(kRed);
    TLatex latIB ; latIB .SetTextSize(0.025); latIB .SetTextAlign(22); latIB .SetTextColor(kRed);
    for(int s=0;s<64;++s){
      double xC=(s%32)*8+4;
      double yC=(s<32)? 72. : 24.;
      latSec.DrawLatex(xC,yC,Form("S%d",s));
      for(int ib=0; ib<6; ++ib){
        double yIB=(s<32)? (48+ib*8+4) : ((5-ib)*8+4);
        latIB.DrawLatex(xC,yIB,Form("%d",ib));
      }
    }

    ensure_dir(outPNG.parent_path());
    c.SaveAs(outPNG.string().c_str());
  }

public:
  bool process(TObject* obj) override
  {
    if(!obj->InheritsFrom(TH1::Class())) return false;

    string hname=obj->GetName();
    bool etaPhi = (hname.find("h_EMC_EtaPhiMap_")==0);
    if(!etaPhi && hname.find("h_EMC_")!=0) return false;

    fs::path outPng = outDir/(hname+".png");
    try{
      if(etaPhi && obj->InheritsFrom(TH2::Class()))
        drawEtaPhiWithMask(static_cast<TH2*>(obj),outPng);
      else if(obj->InheritsFrom(TH2::Class()))
        saveHist2D(static_cast<TH2*>(obj),outPng);
      else
        saveHist1D(static_cast<TH1*>(obj),outPng);

      log::info("EMCal saved "+outPng.string());
    }catch(const std::exception& e){
      log::err("EMCal failed for "+hname+": "+e.what());
    }
    return true;
  }
};

// --------------------------  HCal QA  ----------------------------------
class HcalQA : public TriggerQA {
public:
  using TriggerQA::TriggerQA;

  bool process(TObject* obj) override
  {
    if(!obj->InheritsFrom(TH1::Class())) return false;

    string hname = obj->GetName();
    bool isI = (hname.find("h_IHCAL_")==0);
    bool isO = (hname.find("h_OHCAL_")==0);
    if(!isI && !isO) return false;

    fs::path subDir = outDir/(isI? "IHCal":"OHCal");
    ensure_dir(subDir);
    fs::path outPng = subDir/(hname+".png");

    if(obj->InheritsFrom(TH2::Class()))
      saveHist2D(static_cast<TH2*>(obj),outPng);
    else
      saveHist1D(static_cast<TH1*>(obj),outPng);

    log::info(string(isI?"IHCal":"OHCal")+" saved "+outPng.string());
    return true;
  }
};

// --------------------------  sEPD QA  ----------------------------------
class SepdQA : public TriggerQA {
public:
  using TriggerQA::TriggerQA;
  bool process(TObject* obj) override
  {
    string hname=obj->GetName();
    if(hname.find("sEPD")==string::npos) return false;

    fs::path outPng = outDir/(hname+".png");
    if(obj->InheritsFrom(TH2::Class()))
      saveHist2D(static_cast<TH2*>(obj),outPng);
    else
      saveHist1D(static_cast<TH1*>(obj),outPng);

    log::info("sEPD saved "+outPng.string());
    return true;
  }
};

// ---------------------------  MBD QA  ----------------------------------
class MbdQA : public TriggerQA {
public:
  using TriggerQA::TriggerQA;

  bool process(TObject* obj) override
  {
    //--------------------------------------------------------------------
    // 0.  Reject anything that is not MBD‑related
    //--------------------------------------------------------------------
    const std::string hname = obj->GetName();
    if (hname.find("MBD") == std::string::npos) return false;

    //--------------------------------------------------------------------
    // 1.  One‑dimensional MBD histograms → original handling
    //--------------------------------------------------------------------
    if (!obj->InheritsFrom(TH2::Class())) {
      fs::path outPng = outDir / (hname + ".png");
      saveHist1D(static_cast<TH1*>(obj), outPng);
      log::info("MBD saved " + outPng.string());
      return true;
    }

    //--------------------------------------------------------------------
    // 2.  Two‑dimensional hit‑maps: cache until both North & South exist
    //--------------------------------------------------------------------
    TH2*  h2  = static_cast<TH2*>(obj);
    auto* hCl = static_cast<TH2*>(h2->Clone(Form("%s_clone", h2->GetName())));
    hCl->SetDirectory(nullptr);                        // keep after file closes
    hCl->SetStats(0);                                  // no stats box

    const bool isSouth = (hname.find("_South_") != std::string::npos);
    const std::string cacheKey = trigger;              // one cache per trigger

    struct Pair { TH2* south = nullptr; TH2* north = nullptr; };
    static std::unordered_map<std::string, Pair> cache;

    Pair& p = cache[cacheKey];
    if (isSouth) p.south = hCl; else p.north = hCl;

    // --- If the partner map is not here yet, just return --------------
    if (!p.south || !p.north) {
      log::info("MBD cached " + hname + " (waiting for partner map)");
      return true;
    }

    //--------------------------------------------------------------------
    // 3.  Both maps present → make the side‑by‑side canvas
    //--------------------------------------------------------------------
    fs::path outPng = outDir / "MBD_Hitmap_SouthNorth_" + trigger + ".png";
    gStyle->SetOptStat(0);

    TCanvas c("cMBD", "MBD South & North Hit‑Maps", 1100, 600);
    c.Divide(2, 1, 0.01, 0.01);                     // minimal gaps

    c.cd(1);
    p.south->Draw("POLY COLZ");                     // hex bin‑colours
    p.south->SetTitle("MBD South Hit‑map");

    c.cd(2);
    p.north->Draw("POLY COLZ");
    p.north->SetTitle("MBD North Hit‑map");

    ensure_dir(outPng.parent_path());
    c.SaveAs(outPng.string().c_str());
    log::ok("MBD combined map saved " + outPng.string());

    //--------------------------------------------------------------------
    // 4.  Also save *individual* PNGs exactly as before
    //--------------------------------------------------------------------
    saveHist2D(p.south, outDir / (p.south->GetName() + std::string(".png")));
    saveHist2D(p.north, outDir / (p.north->GetName() + std::string(".png")));

    //--------------------------------------------------------------------
    // 5.  Clean‑up the cache entry
    //--------------------------------------------------------------------
    delete p.south;
    delete p.north;
    cache.erase(cacheKey);
    return true;
  }
};

// --- definition of the static cache map
std::unordered_map<std::string, MbdQA::Pair> MbdQA::cache;


// ======================================================================
//                         MAIN DRIVER
// ======================================================================
void analyzeRun24auau()
{
  try {
    // -------------------------------------------------------- open file
    std::unique_ptr<TFile> in(TFile::Open(kInputFile.c_str(),"READ"));
    if(!in || in->IsZombie())
      throw std::runtime_error("Cannot open "+kInputFile);
    log::ok("Opened "+kInputFile);

    // --------------------------------------------------- Pass 0 catalogue
    fs::path listFile = fs::path(kOutputBase)/"AllHistogramNames.txt";
    ensure_dir(listFile.parent_path());
    std::ofstream list(listFile);
    list << "# Histogram catalogue for " << kInputFile << "\n";

    TIter itDir0(in->GetListOfKeys());
    while(auto* kdir=dynamic_cast<TKey*>(itDir0())){
      if(strcmp(kdir->GetClassName(),"TDirectoryFile")) continue;
      string trig=kdir->GetName();
      TDirectory* d=(TDirectory*)kdir->ReadObj();
      TIter itH(d->GetListOfKeys());
      while(auto* kh=dynamic_cast<TKey*>(itH())){
        string hname=kh->GetName();
        std::cout<<"["<<trig<<"] "<<hname<<"\n";
        list<<"["<<trig<<"] "<<hname<<"\n";
      }
    }
    log::ok("Histogram list written to "+listFile.string());

    // ------------------------------------------ CSV for invariant mass
    fs::path csvPath = fs::path(kOutputBase)/"InvariantMassSummary.csv";
    ensure_dir(csvPath.parent_path());
    std::ofstream csv(csvPath);
    csv<<"trigger,E,Chi,Asym,pTlo,pThi,meanPi0,errPi0,sigmaPi0,errSigmaPi0,"
          "meanEta,errEta,sigmaEta,errSigmaEta\n";

    // ------------------------------------------------------ Pass 1 QA
    TIter itDir(in->GetListOfKeys());
    while(auto* kdir=dynamic_cast<TKey*>(itDir())){
      if(strcmp(kdir->GetClassName(),"TDirectoryFile")) continue;
      string trig=kdir->GetName();
      if(!kTriggersWanted.count(trig)){
        log::info("Skipping trigger "+trig+" (not requested)"); continue;
      }

      log::info("Processing trigger "+trig);
      TDirectory* dTrig=(TDirectory*)kdir->ReadObj();

      // -------- create trigger output skeleton
      fs::path trigBase = fs::path(kOutputBase)/trig;
      ensure_dir(trigBase/"Correlations");
      ensure_dir(trigBase/"EMCal");
      ensure_dir(trigBase/"IHCal");
      ensure_dir(trigBase/"OHCal");
      ensure_dir(trigBase/"MBD");
      ensure_dir(trigBase/"sEPD");

      // -------- instantiate QA handlers
      std::vector<std::unique_ptr<TriggerQA>> qa;
      qa.emplace_back(std::make_unique<Pi0QA>(trig,trigBase/"Correlations",csv));
      qa.emplace_back(std::make_unique<EmcalQA>(trig,trigBase/"EMCal"));
      qa.emplace_back(std::make_unique<HcalQA>(trig,trigBase));
      qa.emplace_back(std::make_unique<SepdQA>(trig,trigBase/"sEPD"));
      qa.emplace_back(std::make_unique<MbdQA >(trig,trigBase/"MBD"));

      // -------- histogram loop
      TIter itH(dTrig->GetListOfKeys());
      while(auto* kh=dynamic_cast<TKey*>(itH())){
        std::unique_ptr<TObject> obj(kh->ReadObj());
        bool handled=false;
        for(auto& h : qa)
          if(h->process(obj.get())) { handled=true; break; }

        if(!handled){
          fs::path misc = trigBase/"Misc"; ensure_dir(misc);
          if(obj->InheritsFrom(TH2::Class()))
            saveHist2D(static_cast<TH2*>(obj.get()), misc/(string(obj->GetName())+".png"));
          else if(obj->InheritsFrom(TH1::Class()))
            saveHist1D(static_cast<TH1*>(obj.get()), misc/(string(obj->GetName())+".png"));
          log::warn("Saved unmatched histogram → Misc/"+string(obj->GetName())+".png");
        }
      } // hist loop
    } // trigger loop

    csv.close();
    log::ok("Finished – outputs in "+kOutputBase);
  }
  catch(const std::exception& e){
    log::err("Fatal: "+string(e.what()));
  }
}

// ----------------------------------------------------------------------
//  Allow `root -l -b -q analyzeRun24auau.cpp`
// ----------------------------------------------------------------------
void run_analyzeRun24auau(){ analyzeRun24auau(); }
run_analyzeRun24auau();
