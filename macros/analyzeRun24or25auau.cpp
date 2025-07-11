// analyzeRun24auau.cpp  – ROOT ≥ 6, C++17
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

#include <filesystem>
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
const string kInputFile =
    "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/"
    "emcal_sepd_analysis_run54280_c0_DST_CALO_run2auau_new_2024p007-00054280-00000.root";

const string kOutputBase =
    "/Users/patsfan753/Desktop/auauAnalysis/emcalSEPDcorrelations/output";

std::set<string> kTriggersWanted{ "MBD_NandS_geq_2" };

/*  <<<   π0‐fit master switch   >>>                                         *
 *  false  → spectra are drawn, but *no* TF1 fit is attempted and            *
 *           InvariantMassSummary.csv is left empty (except header).         *
 *  true   → run the Gaussian‑plus‑poly fit and fill CSV.                    */
constexpr bool kDoPi0Fit = false;
// ───────────────────────────────────────────────


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
// ║ 6.  NS MAP CACHE                            ║
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
// ║ 7.  QA  BASE CLASS                          ║
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

// ——— π0 QA —————————————————————————————————————————
class Pi0QA : public QA{
public:
  Pi0QA(string t,fs::path b,const CentList& s,std::ofstream& csv):
        QA(t,b,s),csv(csv){}
  bool process(TObject* o) override
  {
    if(!o->InheritsFrom(TH1::Class())) return false;
    string n=o->GetName(); if(n.rfind("mInv_",0)!=0) return false;

    CutKey ck; if(!decodeInvName(n,ck)) return false;
    string sl=sliceKey(n);

    fs::path sub = "EMCal/pi0QA";
    sub/=("E"+sf3(ck.E)+"_Chi"+sf3(ck.chi)+"_Asym"+sf3(ck.asy));
    if(ck.pLo>=0&&ck.pHi>=0) sub/=("pT_"+sf3(ck.pLo)+"_to_"+sf3(ck.pHi));

    auto save=[&](const string& slice)
    {
      fs::path out=cPath(root,slice,sub)/(n+".png");
      ensure_dir(out.parent_path());
      TH1* h=static_cast<TH1*>(o);

      bool ok=false;
      std::array<double,8> p{};
      if(kDoPi0Fit){ p = fitPi0(h, ok); }

      TCanvas c; h->SetStats(0); h->Draw();
      if(kDoPi0Fit && ok){ h->GetFunction("f")->SetLineColor(kRed); }
      c.SaveAs(out.string().c_str());

      if(kDoPi0Fit && ok){
        csv<<trig<<","<<ck.E<<","<<ck.chi<<","<<ck.asy<<","
           <<ck.pLo<<","<<ck.pHi<<","<<p[0]<<","<<p[1]<<","
           <<p[2]<<","<<p[3]<<","<<p[4]<<","<<p[5]<<","
           <<p[6]<<","<<p[7]<<"\n";
      }
    };

    save(sl);
    return true;
  }
private:
  std::ofstream& csv;
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
      fs::path out=cPath(root,slice,"Correlations")/(n+".png");
      save2D(static_cast<TH2*>(o),out);
    };
    save(sl);
    return true;
  }
};

// ——— EMCal QA ————————————————————————————————
class EmcalQA : public QA{
public: using QA::QA;
  bool process(TObject* o) override
  {
    if(!o->InheritsFrom(TH1::Class())) return false;
    string n=o->GetName(); if(n.rfind("h_EMC_",0)!=0) return false;
    string sl=sliceKey(n);
    bool etaPhi=n.find("_EtaPhiMap_")!=string::npos;

    auto save=[&](const string& slice){
      fs::path out=cPath(root,slice,"EMCal")/(n+".png");
      if(etaPhi && o->InheritsFrom(TH2::Class()))
      {
        auto* h=static_cast<TH2*>(o);
        std::unique_ptr<TH2> h2(static_cast<TH2*>(h->Clone()));
        h2->SetDirectory(nullptr); h2->SetStats(0); h2->SetContour(99);
        for(int iphi=0;iphi<256;++iphi)
          for(int ieta=0;ieta<96;++ieta)
            if(isBadBoard(sector_from_idx(ieta,iphi),ib_from_idx(ieta,iphi)))
              h2->SetBinContent(h2->FindBin(iphi,ieta),-9999.);
        TCanvas c("c","",1600,1200); c.SetRightMargin(0.15);
        h2->Draw("COLZ");
        TLine l(0,48,256,48); l.SetLineColor(kWhite); l.SetLineWidth(2); l.Draw();
        ensure_dir(out.parent_path()); c.SaveAs(out.string().c_str());
      }
      else if(o->InheritsFrom(TH2::Class())) save2D(static_cast<TH2*>(o),out);
      else                                   save1D(static_cast<TH1*>(o),out);
    };
    save(sl);
    return true;
  }
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
        TCanvas c("c","",1300,800); c.SetRightMargin(0.16); h2->Draw("COLZ");
        ensure_dir(out.parent_path()); c.SaveAs(out.string().c_str());
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
      c.cd(1); out.s->SetTitle(DERIVED::titleSouth); out.s->Draw(("COL POLZ");
      c.cd(2); out.n->SetTitle(DERIVED::titleNorth); out.n->Draw(("COL POLZ");
      ensure_dir(outPng.parent_path()); c.SaveAs(outPng.string().c_str());
    };
    save(sl);
    return true;
  }
protected:
  NSCache<MapPair>& cache;
};

struct MBDTag{
  // accept only true‑MBD histograms: must contain “MBD” and NOT “sEPD”
  static bool accept(const string& s)
  {
    return s.find("MBD")  != string::npos &&   // MBD present
           s.find("sEPD") == string::npos;     // but no sEPD substring
  }
  static constexpr const char* subdir = "MBD";
  static string fileName(const string& t){ return "MBD_Hitmap_NS_" + t; }
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



// ╔══════════════════════════════════════════════╗
// ║ 9.  MAIN DRIVER                              ║
// ╚══════════════════════════════════════════════╝
void analyzeRun24or25auau()
{
  gStyle->SetOptStat(0);

  log::banner("sPHENIX Run‑24 Au+Au QA – Enhanced Macro");

  std::unique_ptr<TFile> in(TFile::Open(kInputFile.c_str(),"READ"));
  if(!in||in->IsZombie()){ log::err("Cannot open "+kInputFile); return; }
  log::ok("Input file opened");

  CentList slices = discoverSlices(in.get());
  {
    std::ostringstream o; o<<"Centrality slices: ";
    for(auto& s:slices) o<<s<<"  ";
    log::info(o.str());
  }

  // ─── Pass 0 catalogue ─────────────────────────
  log::banner("Pass 0 – Catalogue");
  fs::path catTxt=fs::path(kOutputBase)/"AllHistogramNames.txt";
  ensure_dir(catTxt.parent_path());
  std::ofstream cat(catTxt);
  cat<<"# Histogram catalogue for "<<kInputFile<<"\n";

  std::unordered_map<string,int> hCount;
  TIter it0(in->GetListOfKeys());
  while(auto* kd=dynamic_cast<TKey*>(it0())){
    if(strcmp(kd->GetClassName(),"TDirectoryFile")) continue;
    string trg=kd->GetName(); if(!kTriggersWanted.count(trg)) continue;
    TDirectory* d=(TDirectory*)kd->ReadObj();
    TIter itH(d->GetListOfKeys());
    while(auto* kh=dynamic_cast<TKey*>(itH())){
      cat<<"["<<trg<<"] "<<kh->GetName()<<"\n"; ++hCount[trg];
    }
  }
  log::ok("Histogram list written → "+catTxt.string());

  // ─── Pass 1 QA ────────────────────────────────
  log::banner("Pass 1 – QA Production");

  fs::path csvPath=fs::path(kOutputBase)/"InvariantMassSummary.csv";
  ensure_dir(csvPath.parent_path());
  std::ofstream csv(csvPath);
  csv<<"trigger,E,Chi,Asym,pTlo,pThi,meanPi0,errPi0,sigmaPi0,errSigmaPi0,"
        "meanEta,errEta,sigmaEta,errSigmaEta\n";

  NSCache<MapPair> mbdCache, sepdCache;

  struct Counter{ int hTotal=0,hUsed=0; };
  std::unordered_map<string,Counter> trigStat;

  TIter itDir(in->GetListOfKeys());
  while(auto* kd=dynamic_cast<TKey*>(itDir())){
    if(strcmp(kd->GetClassName(),"TDirectoryFile")) continue;
    string trg=kd->GetName(); if(!kTriggersWanted.count(trg)) continue;
    TDirectory* dTrig=(TDirectory*)kd->ReadObj();

    // set up output directories
    for(auto& s:slices){
      fs::path b=fs::path(kOutputBase)/trg;
      if(s!="Inclusive") b/=("Cent_"+s);
      for(auto sub:{"Correlations","EMCal/pi0QA","IHCal","OHCal","MBD","sEPD"})
        ensure_dir(b/sub);
    }

    // instantiate QA modules
    std::vector<std::unique_ptr<QA>> qa;
    fs::path base=fs::path(kOutputBase)/trg;
    qa.emplace_back(std::make_unique<Pi0QA >(trg,base,slices,csv));
    qa.emplace_back(std::make_unique<CorrQA>(trg,base,slices));
    qa.emplace_back(std::make_unique<EmcalQA>(trg,base,slices));
    qa.emplace_back(std::make_unique<HcalQA >(trg,base,slices));
    qa.emplace_back(std::make_unique<MbdQA  >(trg,base,slices,mbdCache));
    qa.emplace_back(std::make_unique<SepdQA >(trg,base,slices,sepdCache));

    // histogram loop
    TIter itH(dTrig->GetListOfKeys());
    while(auto* kh=dynamic_cast<TKey*>(itH())){
      log::trace("Handling ["+trg+"] histogram \""+string(kh->GetName())+"\"");
      TObject* obj=kh->ReadObj();
      if(obj->InheritsFrom(TH1::Class()))
        static_cast<TH1*>(obj)->SetDirectory(nullptr);

      ++trigStat[trg].hTotal;
      for(auto& m:qa) if(m->process(obj)){ ++trigStat[trg].hUsed; break; }
    }

    log::ok("Trigger "+trg+": processed "+
            std::to_string(trigStat[trg].hTotal)+" objects");
  }

  csv.close();

  // ─── Summary ──────────────────────────────────
  log::banner("Summary");
  std::cout<<term::CLR_BOLD
           <<std::left<<std::setw(20)<<"Trigger"
           <<std::right<<std::setw(12)<<"Total"
           <<std::setw(12)<<"Written"<<term::CLR_RST<<"\n";
  for(auto& [t,c]:trigStat)
    std::cout<<std::left<<std::setw(20)<<t
             <<std::right<<std::setw(12)<<c.hTotal
             <<std::setw(12)<<c.hUsed<<"\n";

  log::ok("All outputs under "+kOutputBase);
}

