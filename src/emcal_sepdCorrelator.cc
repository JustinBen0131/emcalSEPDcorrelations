//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator
//  Implementation file  – no duplicated definitions
//==========================================================================
#include "emcal_sepdCorrelator.h"
//––– Fun4All / PHOOL -------------------------------------------------------
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <jetbase/JetContainer.h>
#include <array>
//––– ROOT & CLHEP ----------------------------------------------------------
#include <TProfile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH2Poly.h>
#include <CLHEP/Vector/ThreeVector.h>
//––– CDB access ------------------------------------------------------------
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
//––– sPHENIX objects -------------------------------------------------------
#include <globalvertex/GlobalVertex.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawClusterUtility.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <epd/EpdGeom.h>
#include <epd/EpdReco.h>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasClassifier.h>   // optional but handy

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/Eventplaneinfov1.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

// Standard C++ -------------------------------------------------------------
#include <atomic>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>


#ifdef _OPENMP
  #include <omp.h>
#endif

//–––––––– helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#define CLR_BLUE   "\033[1;34m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_RESET  "\033[0m"

#undef  LOG
#define LOG(lvl, colour, msg)                                           \
  do {                                                                  \
    if (static_cast<int>(Verbosity()) >= static_cast<int>(lvl))         \
      std::cout << colour << msg << CLR_RESET << std::endl;             \
  } while (false)

/** Always print, independent of Verbosity() */
#define PROGRESS(MSG)                                                     \
    do {                                                                  \
        if (static_cast<int>(Verbosity()) >= 1)                           \
            std::cout << CLR_CYAN << MSG << CLR_RESET << std::endl;       \
    } while (false)
//==========================================================================
//  ctor
//==========================================================================
emcal_sepdCorrelator::emcal_sepdCorrelator(const std::string& outFile)
: SubsysReco("emcal_sepdCorrelator"),
  Outfile(outFile)
{
  if (Outfile.empty())
  {
    std::cerr << "[FATAL] output filename is empty.\n";
    std::exit(EXIT_FAILURE);
  }
}

emcal_sepdCorrelator::~emcal_sepdCorrelator()
{
  
}

int emcal_sepdCorrelator::Init(PHCompositeNode* topNode)
{
  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – starting");

  /* 0.  book-keeping & QA histograms --------------------------------- */
  out = new TFile(Outfile.c_str(), "RECREATE");
  LOG(1, CLR_GREEN, "[Init] opened output file: " << Outfile);

  trigAna = new TriggerAnalyzer();
  LOG(1, CLR_GREEN, "[Init] booking scalar QA histograms …");
  createHistos_Data();

  /* 1.  optional DST node-tree dump ---------------------------------- */
  if (Verbosity() >= 2)           // ← adjust threshold as desired
  {
    std::cout << CLR_CYAN
              << "\n[Init] ── DST node-tree dump ────────────────────────────"
              << CLR_RESET << std::endl;

    /* depth-first walk implemented with a std::function so that the
       lambda can recurse without shadowing problems                       */
    std::function<void(PHCompositeNode*, int)> dumpTree =
      [&](PHCompositeNode* node, int depth)
    {
      if (!node) return;

      /* print this node ------------------------------------------------ */
      const std::string indent(depth * 3, ' ');
      std::cout << indent << node->getName()
                << " (" << node->getType() << ")";

      /* for IO-data nodes also print the contained class name ---------- */
      if (node->getType() == "PHIODataNode")
        std::cout << " <" << node->getClass() << '>';

      std::cout << '\n';

      /* iterate over children with the *real* PHOOL API ---------------- */
      PHNodeIterator it(node);
      auto& kids = it.ls();                       // PHPointerList<PHNode>
      for (size_t i = 0; i < kids.length(); ++i)
      {
        PHNode* child = kids[i];
        if (!child) continue;

        if (auto* comp = dynamic_cast<PHCompositeNode*>(child))
        {
          dumpTree(comp, depth + 1);              // recurse
        }
        else
        {
          const std::string ind2((depth + 1) * 3, ' ');
          std::cout << ind2 << child->getName()
                    << " (" << child->getType() << ")";
          if (child->getType() == "PHIODataNode")
            std::cout << " <" << child->getClass() << '>';
          std::cout << '\n';
        }
      }
    };

    /* pick the correct root: use argument if non-null, else global ---- */
    PHCompositeNode* root = topNode
                              ? topNode
                              : Fun4AllServer::instance()->topNode();

    dumpTree(root, 0);

    std::cout << CLR_CYAN
              << "[Init] ──────────────────────────────────────────────────\n"
              << CLR_RESET << std::endl;
  }

  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – done");
  return Fun4AllReturnCodes::EVENT_OK;
}

int emcal_sepdCorrelator::InitRun(PHCompositeNode* topNode)
{
  /* 0. banner -------------------------------------------------------- */
  const uint64_t run = recoConsts::instance()->get_uint64Flag("TIMESTAMP", 0);
  LOG(1, CLR_BLUE, "[InitRun] ------------------------------------------------------------");
  LOG(1, CLR_BLUE, "[InitRun] Starting InitRun  –  TIMESTAMP = " << run);

  /* 1. geometry‑dependent hit‑maps – only once per job --------------- */
  if (!m_mapsBooked)
  {
    LOG(1, CLR_GREEN, "[InitRun] booking hit‑maps …");
    bookShapeHitMaps(topNode);
    m_mapsBooked = true;
  }

  /* 2. (lazy) SEPD mapping will be done the first time we see data --- */
  m_sepdMapReady = false;         // force rebuild after run change

  /* 3. sanity‑check user centrality edges --------------------------- */
  if (m_centEdges.empty())
    LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges vector is EMPTY");
  else
  {
    bool mono = std::is_sorted(m_centEdges.begin(), m_centEdges.end());
    if (!mono)
      LOG(0, CLR_YELLOW, "[InitRun] WARNING: centrality edges not monotonic");
  }
  /* ------------------------------------------------------------------ */
  /* Flow‑vn: allocate one FlowAcc vector per detector (8 pT bins)      */
  /* ------------------------------------------------------------------ */
  if (m_flowAcc.empty())                               // only first InitRun
  {
      const std::vector<std::string> dets =
          {"CEMC","IHCAL","OHCAL","HCAL","ALL"};
      for (const auto& d : dets)
        m_flowAcc[d].assign(m_ptBins.size(), {});         // 8 × zero‑initialised
  }
    
  LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
  return Fun4AllReturnCodes::EVENT_OK;
}


void emcal_sepdCorrelator::bookShapeHitMaps(PHCompositeNode* topNode)
{
  auto* mbdg = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
  auto* epdg = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");

  if (!mbdg || !epdg)
    std::cerr << "[WARN] MbdGeom or EpdGeom missing – hit‑maps will be empty\n";

  for (const auto& kv : triggerNameMap)
  {
    const std::string& trig = kv.second;
    if (Verbosity() > 1)
      std::cout << CLR_BLUE << "  ├─ trigger \"" << trig
                << "\" – booking hit‑maps" << CLR_RESET << std::endl;

    HistMap& H = qaHistogramsByTrigger[trig];

    /* ensure (or create) ROOT sub‑dir */
    TDirectory* d = out->GetDirectory(trig.c_str());
    if (!d) d = out->mkdir(trig.c_str());
    d->cd();

    /* MBD */
    H["h_MBD_Hitmap_South_" + trig] = makeMbdHitmap("h_MBD_Hitmap_South_" + trig, mbdg, 0);
    H["h_MBD_Hitmap_North_" + trig] = makeMbdHitmap("h_MBD_Hitmap_North_" + trig, mbdg, 1);

    /* sEPD */
    H["h_sEPD_Hitmap_South_" + trig] = makeEpdHitmap("h_sEPD_Hitmap_South_" + trig, epdg, 0);
    H["h_sEPD_Hitmap_North_" + trig] = makeEpdHitmap("h_sEPD_Hitmap_North_" + trig, epdg, 1);

    /* EMCal / HCal (η,φ) maps */
    H["h_EMC_EtaPhiMap_" + trig]  = new TH2F(("h_EMC_EtaPhiMap_"  + trig).c_str(),
                                             "CEMC tower map;#phi (0...255);#eta (0...95)",
                                             256, 0, 256, 96, 0, 96);

    H["h_IHCAL_EtaPhiMap_" + trig] = new TH2F(("h_IHCAL_EtaPhiMap_" + trig).c_str(),
                                              "IHCAL tower map;#phi (0...63);#eta (0...23)",
                                              64, 0, 64, 24, 0, 24);

    H["h_OHCAL_EtaPhiMap_" + trig] = new TH2F(("h_OHCAL_EtaPhiMap_" + trig).c_str(),
                                              "OHCAL tower map;#phi (0...63);#eta (0...23)",
                                              64, 0, 64, 24, 0, 24);
      
    /* –– centrality‑binned clones –– */
    auto cloneHitMap = [&](const std::string& base, TH2F* src)
    {
        for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
        {
          const int lo = m_centEdges[i], hi = m_centEdges[i + 1];
          std::ostringstream n;  n << base << '_' << lo << '_' << hi << '_' << trig;
          H[n.str()] = dynamic_cast<TH2F*>(src->Clone(n.str().c_str()));
          H[n.str()]->Reset();                       // start empty
          H[n.str()]->SetDirectory(nullptr);         // detach from file
        }
    };

    /* create centrality clones for every calorimeter map */
    cloneHitMap("h_EMC_EtaPhiMap",  static_cast<TH2F*>(H["h_EMC_EtaPhiMap_"  + trig]));
    cloneHitMap("h_IHCAL_EtaPhiMap",static_cast<TH2F*>(H["h_IHCAL_EtaPhiMap_" + trig]));
    cloneHitMap("h_OHCAL_EtaPhiMap",static_cast<TH2F*>(H["h_OHCAL_EtaPhiMap_" + trig]));

  }

  out->cd();
}

void emcal_sepdCorrelator::bookTowerAndClusterQA(const std::string& trig, HistMap& H)
{
  const int nbE = 200; const double eMax = 50.;
  for (auto& ci : m_caloInfo)
  {
    const std::string label = std::get<2>(ci);
    H["h_towerE_" + label] =
        new TH1F(("h_towerE_" + label + "_" + trig).c_str(),
                 (label + " tower E;E [GeV]").c_str(),
                 nbE, 0, eMax);

    if (label == "CEMC")
      H["h_clusterE_EMC"] =
          new TH1F(("h_clusterE_EMC_" + trig).c_str(),
                   "EMC cluster E;E [GeV]", nbE, 0, eMax);
  }
}

void emcal_sepdCorrelator::bookChargeQA(const std::string& trig, HistMap& H)
{
  const int nbQ = 400; const double qMax = 600.;
  H["h_towerQ_SEPD"] = new TH1F(("h_towerQ_SEPD_" + trig).c_str(),
                                "sEPD tower charge;Q [ADC]", nbQ, 0, qMax);
  H["h_charge_MBD"]  = new TH1F(("h_charge_MBD_" + trig).c_str(),
                                "MBD PMT charge sum;Q [ADC]", nbQ, 0, qMax);
}

void emcal_sepdCorrelator::bookEnergyChargeCorrel(const std::string& trig,
                                                  HistMap&           H)
{
  /* helper that produces an auto‑extending TH2F ------------------------- */
  auto book2 = [&](const char* n, const char* t,
                   int nx, double x0, double x1,
                   int ny, double y0, double y1)
  {
    TH2F* h = new TH2F(n, t, nx, x0, x1, ny, y0, y1);
    h->SetCanExtend(TH1::kAllAxes);
    return h;
  };

  /* axis presets – identical to the original implementation ------------- */
  const int    nC   = 240;               // charge axis   (5 ADC  per bin)
  const double cMax = 1200.;             // covers ΣQ for central Au+Au
  const int    nE   = 240;               // energy axis   (5 GeV per bin)
  const double eMax = 1200.;             // covers ΣEt for all subsystems

  /* --------------------------------------------------------------------
   * 1)  keep the *original* un‑binned histograms (exact names preserved)
   * ------------------------------------------------------------------ */
  H["h_SEPD_vs_CEMC"] = book2(("h_SEPD_vs_CEMC_" + trig).c_str(),
                              "sEPD Q vs CEMC #SigmaE",
                              nC, 0, cMax, nE, 0, eMax);
  H["h_SEPD_vs_IHCAL"] = book2(("h_SEPD_vs_IHCAL_" + trig).c_str(),
                               "sEPD Q vs IHCAL #SigmaE",
                               nC, 0, cMax, nE, 0, eMax);
  H["h_SEPD_vs_OHCAL"] = book2(("h_SEPD_vs_OHCAL_" + trig).c_str(),
                               "sEPD Q vs OHCAL #SigmaE",
                               nC, 0, cMax, nE, 0, eMax);
  H["h_SEPD_vs_MBD"]   = book2(("h_SEPD_vs_MBD_"  + trig).c_str(),
                               "sEPD Q vs MBD #SigmaQ",
                               nC, 0, cMax, nC, 0, cMax);

  H["h_MBD_vs_CEMC"]   = book2(("h_MBD_vs_CEMC_" + trig).c_str(),
                               "MBD #SigmaQ vs CEMC #SigmaE",
                               nC, 0, cMax, nE, 0, eMax);
  H["h_MBD_vs_IHCAL"]  = book2(("h_MBD_vs_IHCAL_" + trig).c_str(),
                               "MBD #SigmaQ vs IHCAL #SigmaE",
                               nC, 0, cMax, nE, 0, eMax);
  H["h_MBD_vs_OHCAL"]  = book2(("h_MBD_vs_OHCAL_" + trig).c_str(),
                               "MBD #SigmaQ vs OHCAL #SigmaE",
                               nC, 0, cMax, nE, 0, eMax);

  H["h_SEPD_S_vs_CEMC_South"] = book2(("h_SEPD_S_vs_CEMC_South_" + trig).c_str(),
                                      "#SigmaQ_{sEPD South}  vs  #SigmaEt_{CEMC #eta<0}",
                                      nC, 0, cMax, nE, 0, eMax);
  H["h_SEPD_N_vs_CEMC_North"] = book2(("h_SEPD_N_vs_CEMC_North_" + trig).c_str(),
                                      "#SigmaQ_{sEPD North}  vs  #SigmaEt_{CEMC #eta>0}",
                                      nC, 0, cMax, nE, 0, eMax);

  /* --------------------------------------------------------------------
   * 2)  centrality‑binned clones – one per {lo,hi} range the user gave
   * ------------------------------------------------------------------ */
  auto addClone = [&](const std::string& base,
                      int lo, int hi,
                      const char* title,
                      int nx,double x0,double x1,
                      int ny,double y0,double y1)
  {
    std::ostringstream key;
    key << base << '_' << lo << '_' << hi << '_' << trig;
    H[key.str()] = book2(key.str().c_str(), title, nx,x0,x1, ny,y0,y1);

    /* remember the index so we can look it up quickly at fill time */
    std::ostringstream tag; tag << '_' << lo << '_' << hi;
    m_centIdxCache[tag.str()] = 1;        // value unused – we only need the key
  };

  /* loop over consecutive edges: [e0,e1), [e1,e2), … ------------------- */
  for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
  {
    const int lo = m_centEdges[i];
    const int hi = m_centEdges[i + 1];

    addClone("h_SEPD_vs_CEMC", lo,hi,
             "sEPD Q vs CEMC #SigmaE", nC,0,cMax, nE,0,eMax);
    addClone("h_SEPD_vs_IHCAL",lo,hi,
             "sEPD Q vs IHCAL #SigmaE",nC,0,cMax, nE,0,eMax);
    addClone("h_SEPD_vs_OHCAL",lo,hi,
             "sEPD Q vs OHCAL #SigmaE",nC,0,cMax, nE,0,eMax);
    addClone("h_SEPD_vs_MBD",  lo,hi,
             "sEPD Q vs MBD #SigmaQ",  nC,0,cMax, nC,0,cMax);

    addClone("h_MBD_vs_CEMC", lo,hi,
             "MBD #SigmaQ vs CEMC #SigmaE", nC,0,cMax, nE,0,eMax);
    addClone("h_MBD_vs_IHCAL",lo,hi,
             "MBD #SigmaQ vs IHCAL #SigmaE",nC,0,cMax, nE,0,eMax);
    addClone("h_MBD_vs_OHCAL",lo,hi,
             "MBD #SigmaQ vs OHCAL #SigmaE",nC,0,cMax, nE,0,eMax);

    addClone("h_SEPD_S_vs_CEMC_South",lo,hi,
             "#SigmaQ_{sEPD South}  vs  #SigmaEt_{CEMC #eta<0}",
             nC,0,cMax, nE,0,eMax);
    addClone("h_SEPD_N_vs_CEMC_North",lo,hi,
             "#SigmaQ_{sEPD North}  vs  #SigmaEt_{CEMC #eta>0}",
             nC,0,cMax, nE,0,eMax);
  }
}


/* ----------------------------------------------------------------------
 * bookPi0MassSpectra – π0 invariant–mass spectra
 *                      (global  +  centrality‑tagged clones)
 * -------------------------------------------------------------------- */
void emcal_sepdCorrelator::bookPi0MassSpectra(const std::string& trig,
                                              HistMap&           H)
{
  const int    nM   = 150;
  const double mMax = 1.5;    // [GeV/c²]

  auto addHist = [&](const std::string& baseKey)
  {
    /* 1) always keep the original (centrality‑independent) spectrum */
    const std::string hNameGlobal = baseKey + "_" + trig;
    H[hNameGlobal] = new TH1F(hNameGlobal.c_str(),
                              "m_{#gamma#gamma};GeV/c^{2}",
                              nM, 0., mMax);

    /* 2) add one clone for every user‑defined {lo,hi} percentile bin */
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    {
      const int lo = m_centEdges[i];
      const int hi = m_centEdges[i + 1];

      std::ostringstream name;
      name << baseKey << '_' << lo << '_' << hi << '_' << trig;

      H[name.str()] = new TH1F(name.str().c_str(),
                               "m_{#gamma#gamma};GeV/c^{2}",
                               nM, 0., mMax);
    }
  };

  /* ---- pT‑binned spectra ----------------------------------------------- */
  for (auto pt : m_ptBins)
    for (float Emin : m_minClusE)
      for (float chi : m_chi2Cuts)
        for (float a : m_asymCuts)
          addHist(invKey(pt.first, pt.second, Emin, chi, a));

  /* ---- inclusive (all‑pT) spectra -------------------------------------- */
  for (float Emin : m_minClusE)
    for (float chi : m_chi2Cuts)
      for (float a : m_asymCuts)
        addHist(invKey(-1, -1, Emin, chi, a));      // pT = −1 sentinel
}


void emcal_sepdCorrelator::bookFlowQA(const std::string& trig, HistMap& H)
{
  /* one TProfile per {detector, harmonic, cent bin} filled versus pT‑bin index */
  const int nPt = m_ptBins.size();

  auto make = [&](const std::string& det, int n, int lo, int hi)
  {
    std::ostringstream name;
    name << "p_v" << n << '_' << det << '_' << lo << '_' << hi << '_' << trig;
    auto* p = new TProfile(name.str().c_str(),
                           Form("v_{%d} (%s);p_{T} bin index;v_{%d}",n,det.c_str(),n),
                           nPt, 0, nPt, "s");
    p->SetStats(0);
    H[name.str()] = p;
  };

  for (int n : {2,3})
    for (const auto& det : {"CEMC","IHCAL","OHCAL","HCAL","ALL"})
      for (std::size_t i=0;i+1<m_centEdges.size();++i)
        make(det,n,m_centEdges[i],m_centEdges[i+1]);
}


//==========================================================================
//  bookEventPlaneCentralityQA – charge spectra, ψn distributions
//                              + sub‑event resolution proxies  (n = 1,2,3)
//==========================================================================
void
emcal_sepdCorrelator::bookEventPlaneCentralityQA(const std::string& trig,
                                                 HistMap&           H)
{
  /* 1. ΣQ spectra ---------------------------------------------------- */
  H["h_Qsum_MBD"]  = new TH1F(("h_Qsum_MBD_"  + trig).c_str(),
                              "MBD #SigmaQ;#SigmaQ_{MBD} [ADC]",
                              600, 0, 1200);

  H["h_Qsum_sEPD"] = new TH1F(("h_Qsum_sEPD_" + trig).c_str(),
                              "sEPD #SigmaQ;#SigmaQ_{sEPD} [ADC]",
                              600, 0, 1200);

  /* 2. detector‑to‑detector ΣQ map ---------------------------------- */
  H["h_Qsum_MBD_vs_sEPD"] =
      new TH2F(("h_Qsum_MBD_vs_sEPD_" + trig).c_str(),
               "#Sigma Q_{MBD} vs #Sigma Q_{sEPD};#Sigma Q_{MBD};#Sigma Q_{sEPD}",
               300, 0, 1200, 300, 0, 1200);

  /* 3. ψ n distributions (South, n = 1,2,3) ------------------------- */
  H["h_Psi1_sEPD"] =
      new TH1F(("h_Psi1_sEPD_" + trig).c_str(),
               "sEPD #Psi_{1};#Psi_{1} [rad]", 120, -TMath::Pi(), TMath::Pi());

  H["h_Psi2_sEPD"] =
      new TH1F(("h_Psi2_sEPD_" + trig).c_str(),
               "sEPD #Psi_{2};#Psi_{2} [rad]", 120, -TMath::Pi(), TMath::Pi());

  H["h_Psi3_sEPD"] =
      new TH1F(("h_Psi3_sEPD_" + trig).c_str(),
               "sEPD #Psi_{3};#Psi_{3} [rad]", 120, -TMath::Pi(), TMath::Pi());

  /* 4. sub‑event resolution proxies  ⟨cos n(Ψ^N–Ψ^S)⟩ vs ΣQ --------- */
  const int nBinsQ = 12;              // same binning for all three
  H["h_Psi1_res_vs_Qsum"] =
      new TProfile(("h_Psi1_res_vs_Qsum_" + trig).c_str(),
                   "cos (#Psi_{1}^{N}-#Psi_{1}^{S}) vs #Sigma Q_{sEPD};"
                   "#Sigma Q_{sEPD};#LT cosΔ#Psi #GT",
                   nBinsQ, 0, 1200, "s");

  H["h_Psi2_res_vs_Qsum"] =
      new TProfile(("h_Psi2_res_vs_Qsum_" + trig).c_str(),
                   "cos 2(#Psi_{2}^{N}-#Psi_{2}^{S}) vs #Sigma Q_{sEPD};"
                   "#Sigma Q_{sEPD};#LT cos2Δ#Psi #GT",
                   nBinsQ, 0, 1200, "s");

  H["h_Psi3_res_vs_Qsum"] =
      new TProfile(("h_Psi3_res_vs_Qsum_" + trig).c_str(),
                   "cos 3(#Psi_{3}^{N}-#Psi_{3}^{S}) vs #Sigma Q_{sEPD};"
                   "#Sigma Q_{sEPD};#LT cos3Δ#Psi #GT",
                   nBinsQ, 0, 1200, "s");
}


// ----------------------------------------------------------------------
// Book jet QA spectra:
//   (A) 1‑D   max‑E_T
//   (B) 3‑D   E_T • area • Nconst           (+ optional 2‑D lead‑vs‑sub)
//   (C) jet‑flow profiles  v_n^{jet}(pT)    (n = 1,2,3)
// ----------------------------------------------------------------------
void
emcal_sepdCorrelator::bookJetQA(const std::string& trig, HistMap& H)
{
  const int nbEt = 200;  const double etMax = 200.;   // 1 GeV/bin
  const int nbA  =  80;  const double aMax  = 0.80;   // ΔA = 0.01
  const int nbN  = 200;  const double nMax  = 200.;   // 1 constituent/bin

  for (const auto& r : kJetRadii)            // r.first = "r02", …
  {
    /* -------- (A) max‑jet‑E_T (global) ----------------------------- */
    const std::string base1 = "h_maxJetEt_" + std::string(r.first);
    H[base1 + "_" + trig] =
        new TH1F((base1 + "_" + trig).c_str(),
                 ("max jet E_{T} (" + std::string(r.first) +
                  ");E_{T} [GeV]").c_str(),
                 nbEt, 0, etMax);

    /* -------- (B) E_T–area–Nconst (global) ------------------------- */
    const std::string base3 = "h_jetEt_area_nConst_" + std::string(r.first);
    H[base3 + "_" + trig] =
        new TH3F((base3 + "_" + trig).c_str(),
                 ("Jet E_{T} vs area vs N_{const} (" + std::string(r.first) +
                  ");E_{T} [GeV];Area;N_{const}").c_str(),
                 nbEt, 0, etMax, nbA, 0, aMax, nbN, 0, nMax);

    /* -------- OPTIONAL 2‑D lead‑vs‑sub plot (global) --------------- */
    const std::string base2 = "h_leadEt_vs_subEt_" + std::string(r.first);
    H[base2 + "_" + trig] =
        new TH2F((base2 + "_" + trig).c_str(),
                 ("Leading vs sub‑leading jet E_{T} (" + std::string(r.first) +
                  ");E_{T}^{lead} [GeV];E_{T}^{sub} [GeV]").c_str(),
                 nbEt, 0, etMax, nbEt, 0, etMax);

    /* -------- centrality clones for the three objects -------------- */
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    {
      const int lo = m_centEdges[i], hi = m_centEdges[i + 1];
      std::ostringstream n1, n2, n3;
      n1 << base1 << '_' << lo << '_' << hi << '_' << trig;
      n2 << base2 << '_' << lo << '_' << hi << '_' << trig;
      n3 << base3 << '_' << lo << '_' << hi << '_' << trig;

      H[n1.str()] = new TH1F(n1.str().c_str(),
                             ("max jet E_{T} (" + std::string(r.first) +
                              ");E_{T} [GeV]").c_str(),
                             nbEt, 0, etMax);

      H[n2.str()] = new TH2F(n2.str().c_str(),
                             ("Leading vs sub‑leading jet E_{T} (" +
                              std::string(r.first) +
                              ");E_{T}^{lead} [GeV];E_{T}^{sub} [GeV]").c_str(),
                             nbEt, 0, etMax, nbEt, 0, etMax);

      H[n3.str()] = new TH3F(n3.str().c_str(),
                             ("Jet E_{T} vs area vs N_{const} (" +
                              std::string(r.first) +
                              ");E_{T} [GeV];Area;N_{const}").c_str(),
                             nbEt, 0, etMax, nbA, 0, aMax, nbN, 0, nMax);
    }

    /* =================================================================
     * (C)  v_n^{jet}(p_T)  profiles   (n = 1,2,3)
     *     – one per {harmonic, centrality slice, trigger}
     * ================================================================= */
    const int nPtBins = 40;                     // 5 GeV bins up to 200 GeV
    for (int n : {1, 2, 3})
      for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      {
        const int lo = m_centEdges[i], hi = m_centEdges[i + 1];
        std::ostringstream nm;
        nm << "p_v" << n << "_JET_" << r.first << '_'   // e.g. p_v2_JET_r02_20_40_MB
           << lo << '_' << hi << '_' << trig;

        H[nm.str()] = new TProfile(
            nm.str().c_str(),
            Form("v_{%d}^{jet} (%s);p_{T}^{jet} [GeV];v_{%d}",
                 n, r.first, n),
            nPtBins, 0, 200, "s");
      }
  } // radii loop
}


void emcal_sepdCorrelator::createHistos_Data()
{
  for (const auto& kv : triggerNameMap)
  {
    const std::string trig = kv.second;
    if (Verbosity() > 1)
      std::cout << CLR_BLUE << "  ├─ trigger \"" << trig
                << "\" – booking scalar QA" << CLR_RESET << std::endl;

    out->mkdir(trig.c_str())->cd();
    HistMap& H = qaHistogramsByTrigger[trig];

    bookTowerAndClusterQA(trig, H);
    bookChargeQA(trig, H);
    bookFlowQA(trig, H);
    bookEnergyChargeCorrel(trig, H);
    bookPi0MassSpectra(trig, H);
    bookEventPlaneCentralityQA(trig, H);
    bookJetQA(trig, H);
    H["h_vertexZ"] = new TH1F(("h_vertexZ_" + trig).c_str(),
                                "Primary vertex z;z_{vtx} [cm]",
                                240, -60., 60.);
      
    H["h_centrality"] =
          new TH1F(("h_centrality_" + trig).c_str(),
                   "Centrality percentile (MBD);centrality [%];Events",
                   100, 0., 100.);
    out->cd();
  }
}

void emcal_sepdCorrelator::buildSepdChannelMap()
{
  const std::size_t nChan = m_sepd->size();
  m_epdKey.assign(nChan, std::numeric_limits<unsigned>::max());

  const std::string mapName = "SEPD_CHANNELMAP";
  const std::string field   = "epd_channel_map";
  CDBTTree tree{ CDBInterface::instance()->getUrl(mapName) };

  std::size_t nMapped = 0;
  for (std::size_t ch = 0; ch < nChan; ++ch)
  {
    const int tile = tree.GetIntValue(ch, field);   // 0…511 or 999
    if (tile == 999) continue;                      // empty slot
    ++nMapped;
    const unsigned arm = (ch >= 384) ? 1u : 0u;     // 0 = S, 1 = N
    const unsigned id  = arm * 256u + tile;
    m_epdKey[ch] = TowerInfoDefs::encode_epd(id);
  }

  const double frac = 100. * nMapped / nChan;
  LOG(1, CLR_GREEN, "[SEPD‑map] mapped " << nMapped << " / " << nChan
                                         << " channels (" << std::fixed
                                         << std::setprecision(1) << frac << "%)");
  m_sepdMapReady = true;
}

int emcal_sepdCorrelator::process_event(PHCompositeNode* topNode)
{
  /* ------------------------------------------------------------------ */
  /* 0. Banner & running counter                                        */
  /* ------------------------------------------------------------------ */
  ++event_count;
  PROGRESS("=================================   event " << std::setw(3) << event_count << "    ===================================== "
           "(Verb=" << Verbosity() << ")");

  /* ------------------------------------------------------------------ */
  /* 1. Mandatory node check                                            */
  /* ------------------------------------------------------------------ */
  LOG(4, CLR_BLUE, "  [process_event] – node sanity");

  if (!fetchNodes(topNode))
  {
    LOG(4, CLR_YELLOW, "    mandatory node missing → ABORTEVENT");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 2. On-the-fly SEPD map build                                       */
  /* ------------------------------------------------------------------ */
  if (!m_sepdMapReady)
  {
    LOG(5, CLR_BLUE, "    building SEPD channel map …");
    m_sepd = findNode::getClass<TowerInfoContainer>(topNode,
                                                    "TOWERINFO_CALIB_SEPD");
    if (!m_sepd)
    {
      LOG(4, CLR_YELLOW,
          "    SEPD container not yet present – event skipped");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    buildSepdChannelMap();
  }

  /* ------------------------------------------------------------------ */
  /* 3. Trigger decoding                                                */
  /* ------------------------------------------------------------------ */
  trigAna->decodeTriggers(topNode);

  std::vector<std::string> activeTrig;
  if (Verbosity() >= 4) LOG(4, CLR_BLUE, "    Trigger matrix:");

  for (auto& [bitname, key] : triggerNameMap)
  {
    ++m_trigStat[key].tested;
    const bool fired = trigAna->didTriggerFire(bitname);
    if (fired) { activeTrig.push_back(key); ++m_trigStat[key].fired; }

    if (Verbosity() >= 4)
      std::cout << "      • " << std::left << std::setw(25) << bitname
                << " → "
                << (fired ? CLR_GREEN "FIRED" : CLR_YELLOW "–")
                << CLR_RESET << '\n';
  }

  if (activeTrig.empty())
  {
    ++m_evtNoTrig;
    LOG(4, CLR_YELLOW, "    no configured trigger fired – skip");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 4. Vertex-z QA & cut                                               */
  /* ------------------------------------------------------------------ */
  for (const auto& t : activeTrig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_vertexZ"])->Fill(m_vz);

  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
  {
    LOG(4, CLR_YELLOW, "    |vz| = " << std::fabs(m_vz)
                       << " cm ≥ cut (" << m_vzCut << ") – skip");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 5. Detector-level QA                                               */
  /* ------------------------------------------------------------------ */
  LOG(5, CLR_BLUE, "    running detector-level QA");
  doCaloQA(activeTrig);
  doSepdQA(activeTrig);

  if (!std::isfinite(m_vz) || std::abs(m_vz) > 60.0)
  {
      LOG(4, CLR_YELLOW, "    Vertex z = " << m_vz << " out of calibration bounds – skip");
      return Fun4AllReturnCodes::ABORTEVENT;
  }
  /* ------------------------------------------------------------------ */
  /* 6. Centrality lookup & diagnostics                                 */
  /* ------------------------------------------------------------------ */

  /* grab the two nodes we need --------------------------------------- */
  CentralityInfo   *_central        = findNode::getClass<CentralityInfo>(topNode,"CentralityInfo");
  MinimumBiasInfo  *_minimumbiasinfo= findNode::getClass<MinimumBiasInfo>(topNode,"MinimumBiasInfo");

  /* bail out gracefully if either of them is missing ----------------- */
  if (!_central)
  {
      LOG(4, CLR_YELLOW, "    CentralityInfo node **missing** – event skipped");
      return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* make sure the mbd_NS bin is there -------------------------------- */
  if (!_central->has_centrality_bin(CentralityInfo::PROP::mbd_NS))
  {
      LOG(4, CLR_YELLOW, "    CentralityInfo mbd_NS not filled yet – skip");
      return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!_minimumbiasinfo)
  {
      LOG(4, CLR_YELLOW, "    MinimumBiasInfo node **missing** – event skipped");
      return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* now it’s safe to fetch it ---------------------------------------- */
  float centile = _central->get_centrality_bin(CentralityInfo::PROP::mbd_NS);

  /* translate the percentile to an integer “centrality bin” ---------- */
  if (!std::isfinite(centile) || centile < 0.f)
  {
    LOG(4, CLR_YELLOW,
          "    mbd_NS centile invalid – defaulting to full 0-100 %");
    m_centBin = -1;                                   // minimum-bias
  }
  else
  {
    m_centBin = static_cast<int>(centile);
    LOG(5, CLR_GREEN, "    centrality bin = " << m_centBin << '%');
  }
    
  if (!_minimumbiasinfo)
  {
    LOG(4, CLR_YELLOW, "    MinimumBiasInfo node **missing** – event skipped");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* translate the percentile to an integer “centrality bin” ---------- */
  if (!std::isfinite(centile) || centile < 0.f)
  {
    LOG(4, CLR_YELLOW,
          "    mbd_NS centile invalid – defaulting to full 0-100 %");
    m_centBin = -1;                                  // → treat as minimum-bias
  }
  else
  {
    m_centBin = static_cast<int>(centile);
    LOG(5, CLR_GREEN, "    centrality bin = " << m_centBin << '%');
  }
  /* fill centrality histogram */
  if (std::isfinite(centile) && centile >= 0.f && centile <= 100.f)
      for (const auto& t : activeTrig)
        static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_centrality"])
            ->Fill(centile);
    
  doMbdQA (activeTrig);
  doPi0QA (activeTrig);
  fillCorrelations(activeTrig);
  fillFlowHists(activeTrig);

  if (doJetQA(topNode, activeTrig) == Fun4AllReturnCodes::ABORTRUN)
  {
    LOG(2, CLR_YELLOW, "    doJetQA requested ABORTRUN");
    return Fun4AllReturnCodes::ABORTRUN;
  }

  LOG(4, CLR_GREEN, "  [process_event] finished OK");
  return Fun4AllReturnCodes::EVENT_OK;
}


//==========================================================================
//  fetchNodes – check presence of all required nodes & cache pointers
//==========================================================================
bool emcal_sepdCorrelator::fetchNodes(PHCompositeNode* top)
{
    
  /* ------------------------------------------------------------------ */
  /* 0.  Reject non‑minimum‑bias events up‑front                        */
  /* ------------------------------------------------------------------ */
  m_isMinBias = false;   // reset per event

  MinimumBiasInfo* mbInfo =
         findNode::getClass<MinimumBiasInfo>(top,"MinimumBiasInfo");

  if (!mbInfo)
  {
       LOG(1, CLR_YELLOW,
           "  – MinimumBiasInfo node **missing** → skip event");
       return false;                       // ABORTEVENT at caller
  }

  m_isMinBias = mbInfo->isAuAuMinimumBias();

  if (!m_isMinBias)
  {
       LOG(3, CLR_CYAN,
           "  – event is NOT Au+Au minimum‑bias → skip event");
       return false;                       // nothing else to do
  }
    
  /* ––– primary vertex –––––––––––––––––––––––––––––––––––––––––––––––– */
  GlobalVertexMap* vmap = findNode::getClass<GlobalVertexMap>(top,"GlobalVertexMap");
  m_vtx = nullptr; m_vx = m_vy = m_vz = 0.;

  if (!vmap)          { LOG(1, CLR_YELLOW, "  – GlobalVertexMap node **missing** → skip event"); return false; }
  if (vmap->empty())  { LOG(2, CLR_YELLOW, "  – GlobalVertexMap is **empty** → skip event");     return false; }

  m_vtx = vmap->begin()->second;
  if (!m_vtx)         { LOG(1, CLR_YELLOW, "  – vertex pointer null → skip event");              return false; }

  m_vx = m_vtx->get_x();
  m_vy = m_vtx->get_y();
  m_vz = m_vtx->get_z();

  /* ––– calorimeter towers & geometry –––––––––––––––––––––––––––––––– */
  m_calo.clear();
  for (const auto& ci : m_caloInfo)
  {
    const std::string node = std::get<0>(ci),
                      geo  = std::get<1>(ci),
                      lbl  = std::get<2>(ci);

    auto* tw = findNode::getClass<TowerInfoContainer>(top, node);
    auto* ge = findNode::getClass<RawTowerGeomContainer>(top, geo);
    if (!tw || !ge)
    { LOG(2, CLR_YELLOW, "  – missing " << lbl << " nodes → skip"); return false; }

    m_calo[lbl] = { tw, ge, 0. };
  }

  /* ––– remaining detectors –––––––––––––––––––––––––––––––––––––––––– */
  m_sepd     = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_SEPD");
  m_mbdpmts  = findNode::getClass<MbdPmtContainer  >(top, "MbdPmtContainer");
  m_mbdgeom  = findNode::getClass<MbdGeom         >(top, "MbdGeom");
  m_epdgeom  = findNode::getClass<EpdGeom         >(top, "TOWERGEOM_EPD");
  m_epmap    = findNode::getClass<EventplaneinfoMap>(top, "EventplaneinfoMap");
  m_clus     = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");

  const bool ok_sepd = (m_sepd && m_epdgeom);
  const bool ok_mbd  = (m_mbdpmts && m_mbdgeom);

  if (!ok_sepd || !ok_mbd)
  {
    LOG(2, CLR_YELLOW, "  – missing mandatory SEPD and/or MBD nodes → skip event");
    return false;
  }
  return true;
}

TH2Poly* emcal_sepdCorrelator::makeMbdHitmap(const std::string& name,
                                             const MbdGeom*    geom,
                                             int               arm)   // 0 = S, 1 = N
{
  LOG(4, CLR_BLUE, "[makeMbdHitmap] =================================================");
  LOG(4, CLR_BLUE, "[makeMbdHitmap] name=\"" << name << "\"  arm=" << arm);

  /* ------------------------------------------------------------------ */
  /* 0. Create empty container                                          */
  /* ------------------------------------------------------------------ */
  auto *h = new TH2Poly(name.c_str(), ";x (cm);y (cm)", 0, 0, 0, 0);
  if (!geom)
  {
    LOG(4, CLR_YELLOW, "[makeMbdHitmap] geom==nullptr → return empty TH2Poly");
    return h;
  }

  /* ------------------------------------------------------------------ */
  /* 1. Collect PMT centres for the requested arm                       */
  /* ------------------------------------------------------------------ */
  std::vector<std::pair<double,double>> pos;  pos.reserve(64);

  const unsigned totPmt = 128;             // 64 PMTs per arm
  for (unsigned ip = 0; ip < totPmt; ++ip)
  {
    if (geom->get_arm(ip) != arm)
    {
      LOG(6, CLR_YELLOW, "  ip=" << ip << " → other arm – skip");
      continue;
    }

    const double cx = geom->get_x(ip);
    const double cy = geom->get_y(ip);
    if (std::isnan(cx) || std::isnan(cy))
    {
      LOG(6, CLR_YELLOW, "  ip=" << ip << " → NaN centre – skip");
      continue;
    }

    pos.emplace_back(cx, cy);
    LOG(6, CLR_GREEN, "  ip=" << ip << " centre=(" << cx << "," << cy << ")");
  }

  LOG(5, CLR_CYAN, "[makeMbdHitmap] collected " << pos.size()
                   << " centres for arm " << arm);
  if (pos.empty())
  {
    LOG(4, CLR_YELLOW, "[makeMbdHitmap] no valid PMTs – return empty TH2Poly");
    return h;
  }

  /* ------------------------------------------------------------------ */
  /* 2. Derive safe hexagon radius (flat-to-flat = √3·r)                */
  /* ------------------------------------------------------------------ */
  double dMin = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < pos.size(); ++i)
    for (std::size_t j = i + 1; j < pos.size(); ++j)
      dMin = std::min(
        dMin,
        std::hypot(pos[i].first  - pos[j].first,
                   pos[i].second - pos[j].second));

  if (!std::isfinite(dMin) || dMin <= 0.)
  {
    LOG(4, CLR_YELLOW,
        "[makeMbdHitmap] failed to compute dMin (dMin=" << dMin
        << ") – return empty TH2Poly");
    return h;
  }

  const double r = 0.97 * dMin / std::sqrt(3.0);   // 3 % safety margin
  LOG(5, CLR_CYAN, "[makeMbdHitmap] dMin="<<dMin<<"  → hex-r="<<r);

  /* ------------------------------------------------------------------ */
  /* 3. Create one regular flat-top hexagon per PMT centre              */
  /* ------------------------------------------------------------------ */
  h->SetFloat();                                    // allow bin touch summing
  double x[6], y[6];

  std::size_t hexCnt = 0;
  for (const auto& [cx, cy] : pos)
  {
    for (int k = 0; k < 6; ++k)
    {
      const double ang = k * M_PI / 3.0;            // 0°,60°,120°…
      x[k] = cx + r * std::cos(ang);
      y[k] = cy + r * std::sin(ang);
    }
    h->AddBin(6, x, y);
    ++hexCnt;

    LOG(6, CLR_GREEN,
        "  hex#" << hexCnt << " centre=(" << cx << "," << cy << ") added");
  }

  LOG(4, CLR_GREEN, "[makeMbdHitmap] completed – " << hexCnt
                    << " hexagons booked for arm " << arm);
  LOG(4, CLR_BLUE,  "[makeMbdHitmap] =================================================");
  return h;
}


TH2F*
emcal_sepdCorrelator::makeEpdHitmap(const std::string& name,
                                    EpdGeom* /*geom*/, int /*arm*/)
{
  // 24 × φ bins (15° each), 16 radial rings (tile 0 occupies ring 0)
  auto* h = new TH2F(name.c_str(),
                     ";#varphi  [rad];r  [cm]",
                     24, 0, 2*TMath::Pi(),
                     16, 0.15, 3.55);
  h->SetCanExtend(TH1::kAllAxes);   // keeps old “auto‑extend” behaviour
  return h;
}

//==========================================================================
//  doCaloQA – tower and cluster spectra   +   vn (n = 2,3) accumulators
//==========================================================================
void emcal_sepdCorrelator::doCaloQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doCaloQA] processing calorimeter towers");

  /* ------------------------------------------------------------------ */
  /* 0.  clear per‑event FlowAcc containers                             */
  /* ------------------------------------------------------------------ */
  for (auto& kv : m_flowAcc)
    for (auto& acc : kv.second) acc.reset();

  /* ------------------------------------------------------------------ */
  /* 1.  loop over calorimeter subsystems                               */
  /* ------------------------------------------------------------------ */
  for (auto& ck : m_calo)
  {
    const std::string&  lbl   = ck.first;          // "CEMC", "IHCAL", …
    TowerInfoContainer* twC   = ck.second.tw;
    double&             sumE  = ck.second.sumE;
    std::size_t         nHit  = 0;

    for (unsigned ch = 0; ch < twC->size(); ++ch)
    {
      auto* tw = twC->get_tower_at_channel(ch); if (!tw) continue;
      const double e = tw->get_energy();        if (e <= 0) continue;

      /* ---------- bookkeeping & 1‑D spectra ------------------------- */
      sumE += e; ++nHit;
      for (const auto& t : trig)
        static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerE_" + lbl])->Fill(e);

      /* ---------- geometry look‑up --------------------------------- */
      const unsigned key  = twC->encode_key(ch);
      const unsigned iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
      const unsigned ieta = TowerInfoDefs::getCaloTowerEtaBin(key);

      RawTowerGeom* tg  = ck.second.g->get_tower_geometry(key);
      const double  eta = tg ? tg->get_eta() : 0.;
      const double  et  = e / std::cosh(eta);

      /* ---------- arm‑separated ΣE​T -------------------------------- */
      if (lbl == "CEMC")       (eta < 0 ? m_cemcEt_arm[0] : m_cemcEt_arm[1]) += et;
      else if (lbl == "IHCAL") (eta < 0 ? m_ihcalEt_arm[0] : m_ihcalEt_arm[1]) += et;
      else if (lbl == "OHCAL") (eta < 0 ? m_ohcalEt_arm[0] : m_ohcalEt_arm[1]) += et;

      /* ---------- η–φ hit‑maps (global + centrality) ---------------- */
      std::string Hmap;
      if      (lbl == "CEMC")  Hmap = "h_EMC_EtaPhiMap_";
      else if (lbl == "IHCAL") Hmap = "h_IHCAL_EtaPhiMap_";
      else if (lbl == "OHCAL") Hmap = "h_OHCAL_EtaPhiMap_";
      else                     Hmap.clear();                // safety

      if (!Hmap.empty())
      {
        for (const auto& t : trig)
        {
          /* global map */
          static_cast<TH2F*>(qaHistogramsByTrigger[t][Hmap + t])
              ->Fill(iphi, ieta, e);

          /* centrality‑tagged clone */
          if (m_centBin >= 0)
          {
            int lo = 0, hi = 100;
            for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
              if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
              { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }

            std::ostringstream k;
            k << Hmap.substr(0, Hmap.size() - 1)     // drop trailing ‘_’
              << '_' << lo << '_' << hi << '_' << t;

            auto it = qaHistogramsByTrigger[t].find(k.str());
            if (it != qaHistogramsByTrigger[t].end())
              static_cast<TH2F*>(it->second)->Fill(iphi, ieta, e);
          }
        }
      }

      /* ----------------------------------------------------------------
       *  v​n accumulators  (harmonics n = 1, 2, 3) – one entry / pT bin
       * -------------------------------------------------------------- */
      int ptBin = -1;
      for (std::size_t b = 0; b < m_ptBins.size(); ++b)
        if (et >= m_ptBins[b].first && et < m_ptBins[b].second)
        { ptBin = static_cast<int>(b); break; }

      if (ptBin >= 0)
      {
        const double phi = tg ? tg->get_phi() : 0.;

        const double c1 = std::cos(phi),  s1 = std::sin(phi);
        const double c2 = c1 * c1 - s1 * s1;          // cos 2φ
        const double s2 = 2.0 * c1 * s1;              // sin 2φ
        const double c3 = c1 * c2 - s1 * s2;          // cos 3φ
        const double s3 = s1 * c2 + c1 * s2;          // sin 3φ

        auto accumulate = [&](const std::string& det)
        {
          auto& a = m_flowAcc[det][ptBin];
          a.sumW  += et;
          a.qx[1] += et * c1;  a.qy[1] += et * s1;
          a.qx[2] += et * c2;  a.qy[2] += et * s2;
          a.qx[3] += et * c3;  a.qy[3] += et * s3;
        };

        accumulate(lbl);                         // native subsystem
        if (lbl == "IHCAL" || lbl == "OHCAL") accumulate("HCAL");
        accumulate("ALL");
      }
    } /* tower loop */

    LOG(3, CLR_GREEN, "    " << lbl << " : "
                             << nHit << " fired, ΣE = " << sumE);
  }   /* detector loop */

  /* ------------------------------------------------------------------ */
  /* 2.  EMC clusters – simple E spectrum                               */
  /* ------------------------------------------------------------------ */
  if (!m_clus) return;

  RawClusterContainer::ConstRange cr = m_clus->getClusters();
  LOG(3, CLR_GREEN, "    EMC clusters : " << std::distance(cr.first, cr.second));

  for (auto it = cr.first; it != cr.second; ++it)
    for (const auto& t : trig)
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_clusterE_EMC"])
          ->Fill(it->second->get_energy());
}



//==========================================================================
//  doSepdQA – sEPD charge, hit‑map & event‑plane QA
//            (Ψ1, Ψ2, Ψ3 for both arms;   South planes used as reference)
//==========================================================================
void
emcal_sepdCorrelator::doSepdQA(const std::vector<std::string>& trig)
{
  LOG(2, CLR_BLUE, "[doSepdQA] ===============================================");
  if (!m_sepd || !m_epdgeom)
  { LOG(1, CLR_YELLOW, "[doSepdQA] SEPD nodes missing – skip"); return; }

  /* 1. initialise ---------------------------------------------------- */
  m_sepdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  /*  Q‑vector accumulators  ----------------------------------------- */
  float qx1S=0., qy1S=0., qx1N=0., qy1N=0.;   // harmonic 1
  float qx2S=0., qy2S=0., qx2N=0., qy2N=0.;   // harmonic 2
  float qx3S=0., qy3S=0., qx3N=0., qy3N=0.;   // harmonic 3

  /* 2. tower loop ---------------------------------------------------- */
  const auto nChan = m_sepd->size();
  for (unsigned ch = 0; ch < nChan; ++ch)
  {
    if (ch >= 744) continue;                           // electronics limit
    auto* ti = m_sepd->get_tower_at_channel(ch); if (!ti) continue;
    const double w = ti->get_energy();               if (w <= 0) continue;

    const unsigned key = m_epdKey[ch];
    if (key == std::numeric_limits<unsigned>::max()) continue;

    const int    arm = TowerInfoDefs::get_epd_arm(key);  // 0 = S, 1 = N
    const double phi = m_epdgeom->get_phi(key);
    const double r   = m_epdgeom->get_r  (key);

    /* hit‑maps ------------------------------------------------------ */
    const std::string hp = (arm == 0) ? "h_sEPD_Hitmap_South_" : "h_sEPD_Hitmap_North_";
    const double phPlot = (phi < 0) ? phi + 2 * M_PI : phi;  // [0,2π) for plotting
    for (const auto& t : trig)
      static_cast<TH2F*>(qaHistogramsByTrigger[t][hp + t])->Fill(phPlot, r, w);

    /* harmonic‑1 / 2 / 3 Q‑vectors --------------------------------- */
    const double c1 = std::cos(      phi), s1 = std::sin(      phi);
    const double c2 = std::cos(2.0 * phi), s2 = std::sin(2.0 * phi);
    const double c3 = std::cos(3.0 * phi), s3 = std::sin(3.0 * phi);

    if (arm == 0)          /* South (reference) */
    {
      qx1S += w * c1;  qy1S += w * s1;
      qx2S += w * c2;  qy2S += w * s2;
      qx3S += w * c3;  qy3S += w * s3;
      ++nFiredS;
    }
    else                   /* North */
    {
      qx1N += w * c1;  qy1N += w * s1;
      qx2N += w * c2;  qy2N += w * s2;
      qx3N += w * c3;  qy3N += w * s3;
      ++nFiredN;
    }

    /* bookkeeping --------------------------------------------------- */
    m_sepdQ          += w;
    m_sepdQ_arm[arm] += w;
  }

  /* 3. ΣQ spectrum --------------------------------------------------- */
  for (const auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerQ_SEPD"])->Fill(m_sepdQ);

  /* 4. ψn (n=1,2,3) -------------------------------------------------- */
  // --- try to fetch EventPlaneReco map for ψ2 (legacy support) -------
  bool usedMap = false;
  if (m_epmap && !m_epmap->empty())
  {
    auto* epdS = m_epmap->get(EventplaneinfoMap::sEPDS);
    auto* epdN = m_epmap->get(EventplaneinfoMap::sEPDN);
    if (epdS && epdN)
    {
      const auto q2S = epdS->get_qvector(2);
      const auto q2N = epdN->get_qvector(2);
      if ((q2S.first || q2S.second) && (q2N.first || q2N.second))
      {
        Eventplaneinfov1 h;
        m_psi2_S = h.GetPsi(q2S.first, q2S.second, 2);
        m_psi2_N = h.GetPsi(q2N.first, q2N.second, 2);
        usedMap  = true;
      }
    }
  }

  /* tower‑based fallback (all harmonics) ----------------------------- */
  if (!usedMap)
  {
    m_psi2_S = 0.5 * std::atan2(qy2S, qx2S);
    m_psi2_N = 0.5 * std::atan2(qy2N, qx2N);
  }

  m_psi1_S = (std::abs(qx1S) < 1e-9 && std::abs(qy1S) < 1e-9)
             ? 0. : std::atan2(qy1S, qx1S);
  m_psi1_N = (std::abs(qx1N) < 1e-9 && std::abs(qy1N) < 1e-9)
             ? 0. : std::atan2(qy1N, qx1N);

  m_psi3_S = (std::abs(qx3S) < 1e-9 && std::abs(qy3S) < 1e-9)
             ? 0. : (1. / 3.) * std::atan2(qy3S, qx3S);
  m_psi3_N = (std::abs(qx3N) < 1e-9 && std::abs(qy3N) < 1e-9)
             ? 0. : (1. / 3.) * std::atan2(qy3N, qx3N);

  /* 5. fill EP QA histograms ---------------------------------------- */
  const double cos1dPsi = std::cos(    m_psi1_N - m_psi1_S);
  const double cos2dPsi = std::cos(2 * (m_psi2_N - m_psi2_S));
  const double cos3dPsi = std::cos(3 * (m_psi3_N - m_psi3_S));

  for (const auto& t : trig)
  {
    auto& H = qaHistogramsByTrigger[t];

    static_cast<TH1F*>(H["h_Psi1_sEPD"])->Fill(m_psi1_S);
    static_cast<TH1F*>(H["h_Psi2_sEPD"])->Fill(m_psi2_S);
    static_cast<TH1F*>(H["h_Psi3_sEPD"])->Fill(m_psi3_S);

    static_cast<TProfile*>(H["h_Psi1_res_vs_Qsum"])->Fill(m_sepdQ, cos1dPsi);
    static_cast<TProfile*>(H["h_Psi2_res_vs_Qsum"])->Fill(m_sepdQ, cos2dPsi);
    static_cast<TProfile*>(H["h_Psi3_res_vs_Qsum"])->Fill(m_sepdQ, cos3dPsi);
  }

  /* 6. verbose summary ---------------------------------------------- */
  LOG(3, CLR_GREEN,
      "    SEPD ΣQ=" << m_sepdQ
      << "  Ψ1(S)=" << m_psi1_S
      << "  Ψ2(S)=" << m_psi2_S
      << "  Ψ3(S)=" << m_psi3_S
      << "  hits(S,N)=" << nFiredS << ',' << nFiredN);
  LOG(4, CLR_BLUE, "[doSepdQA] ===============================================");
}



//==========================================================================
//  doMbdQA – integrated charge + hex hit‑map
//==========================================================================
void emcal_sepdCorrelator::doMbdQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doMbdQA]");

  m_mbdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  const unsigned nPmt = static_cast<unsigned>(m_mbdpmts->get_npmt());
  for (unsigned ip = 0; ip < nPmt; ++ip)
  {
    MbdPmtHit* p = m_mbdpmts->get_pmt(ip);
    const double q = p->get_q(); if (q <= 0) continue;

    double cx = m_mbdgeom->get_x(ip);
    double cy = m_mbdgeom->get_y(ip);
    const std::string key = (m_mbdgeom->get_arm(ip) == 0 ? "h_MBD_Hitmap_South_"
                                                         : "h_MBD_Hitmap_North_");

    for (auto& t : trig)
      static_cast<TH2Poly*>(qaHistogramsByTrigger[t][key + t])->Fill(cx, cy, q);

    m_mbdQ += q;
    m_mbdQ_arm[m_mbdgeom->get_arm(ip)] += q;
    ++(m_mbdgeom->get_arm(ip) ? nFiredN : nFiredS);
  }

  for (auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_charge_MBD"])->Fill(m_mbdQ);

  LOG(3, CLR_GREEN, "    MBD #SigmaQ = " << m_mbdQ
                     << "  (South PMTs: " << nFiredS
                     << ", North PMTs: " << nFiredN << ")");

  fillCentralityQA(trig);
}

//==========================================================================
//  fillCentralityQA – #SigmaQ histos & correlation map (MBD ↔ sEPD)
//==========================================================================
void emcal_sepdCorrelator::fillCentralityQA(const std::vector<std::string>& trig)
{
  for (auto& t : trig)
  {
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Qsum_MBD"        ])->Fill(m_mbdQ);
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Qsum_sEPD"       ])->Fill(m_sepdQ);
    static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_Qsum_MBD_vs_sEPD"])->Fill(m_mbdQ, m_sepdQ);
  }
}

/* ----------------------------------------------------------------------
 * doPi0QA – γγ invariant‑mass spectra (global + centrality‑tagged)
 *            (improved verbosity)
 * -------------------------------------------------------------------- */
void emcal_sepdCorrelator::doPi0QA(const std::vector<std::string>& trig)
{
  /* ------------------------------------------------------------------ */
  /* 0)  Entrance message                                               */
  /* ------------------------------------------------------------------ */
  LOG(3, CLR_BLUE,
      "[doPi0QA] entered  – "
      << (m_clus ? m_clus->size() : 0)
      << " EMC clusters available in this event");

  /* ---------- early exits ------------------------------------------- */
  if (!m_clus || m_clus->size() < 2)
  {
    LOG(4, CLR_YELLOW, "    [doPi0QA] < 2 clusters – nothing to do");
    return;
  }

  /* ---------- build a filtered cluster cache ------------------------ */
  const float EminMin = *std::min_element(m_minClusE.begin(), m_minClusE.end());
  const float chi2Max = *std::max_element(m_chi2Cuts.begin(), m_chi2Cuts.end());
  const float asymMax = *std::max_element(m_asymCuts.begin(), m_asymCuts.end());

  struct Clu { TLorentzVector v; float E, pt, chi; };
  std::vector<Clu> cl; cl.reserve(m_clus->size());

  std::size_t nRejectedE   = 0;
  std::size_t nRejectedChi = 0;

  for (auto [it, end] = m_clus->getClusters(); it != end; ++it)
  {
    const RawCluster* c = it->second;
    if (c->get_energy() < EminMin) { ++nRejectedE;   continue; }
    if (c->get_chi2()  > chi2Max)  { ++nRejectedChi; continue; }

    const auto eVec = RawClusterUtility::GetEVec(*c, {m_vx, m_vy, m_vz});
    cl.push_back({{}, static_cast<float>(c->get_energy()),
                       static_cast<float>(eVec.perp()),
                       static_cast<float>(c->get_chi2())});
    cl.back().v.SetPtEtaPhiE(cl.back().pt,
                             eVec.pseudoRapidity(),
                             eVec.phi(),
                             cl.back().E);
  }

  if (cl.size() < 2)
  {
    LOG(4, CLR_YELLOW,
        "    [doPi0QA] Cluster filtering left "
        << cl.size() << " usable clusters (rejected "
        << nRejectedE << " for E, " << nRejectedChi
        << " for χ²) – aborting");
    return;
  }

  LOG(3, CLR_GREEN,
      "    [doPi0QA] kept " << cl.size() << " clusters  ("
      << nRejectedE << " rejected by E, "
      << nRejectedChi << " rejected by χ²)");

  /* ---------- determine centrality slice ---------------------------- */
  int lo = 0, hi = 100;
  if (m_centBin >= 0)
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      if (m_centBin >= m_centEdges[i] &&
          m_centBin <  m_centEdges[i + 1])
      { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }

  std::ostringstream tagSS; tagSS << '_' << lo << '_' << hi;   // "_20_40"
  const std::string centTag = tagSS.str();

  LOG(4, CLR_CYAN,
      "    [doPi0QA] centrality = " << m_centBin
      << "%  → slice " << lo << "–" << hi << '%');

  /* ---------- loop over cluster pairs ------------------------------- */
  const std::size_t totalPairs = cl.size() * (cl.size() - 1) / 2;
  std::atomic<std::size_t> pairCnt{0};
  constexpr std::size_t reportEvery = 50'000;

#ifdef _OPENMP
  #pragma omp parallel default(shared)
#endif
  {
    std::map<std::string, CutStat> evtStatLocal;

#ifdef _OPENMP
    #pragma omp for schedule(dynamic,256)
#endif
    for (std::size_t idx = 0; idx < totalPairs; ++idx)
    {
      /* triangular index → (i,j) */
      const std::size_t i = static_cast<std::size_t>((std::sqrt(8.0*idx+1)-1)/2);
      const std::size_t j = idx - i*(i+1)/2 + i + 1;

      const Clu &c1 = cl[i], &c2 = cl[j];
      const float asym = std::fabs(c1.E - c2.E) / (c1.E + c2.E);
      if (asym > asymMax) { ++pairCnt; continue; }            // global veto

      /* ---- helper that fills both histogram flavours ---------------- */
      auto fillBoth = [&](const std::string& baseKey,
                          const std::vector<std::string>& trigList,
                          float mInv)
      {
        for (const auto& t : trigList)
        {
          auto& H = qaHistogramsByTrigger[t];

          /* global spectrum */
          const std::string gKey = baseKey + "_" + t;
          auto itG = H.find(gKey);
          if (itG != H.end())
            static_cast<TH1F*>(itG->second)->Fill(mInv);

          /* centrality‑tagged clone */
          const std::string cKey = baseKey + centTag + "_" + t;
          auto itC = H.find(cKey);
          if (itC != H.end())
            static_cast<TH1F*>(itC->second)->Fill(mInv);
        }
      };

      /* ---- inclusive (all‑pT) spectra ------------------------------ */
      for (float Emin : m_minClusE)
      for (float chiMx: m_chi2Cuts)
      for (float aMx  : m_asymCuts)
      {
        const std::string key = statKey(-1,-1,Emin,chiMx,aMx);
        auto& st = evtStatLocal[key]; ++st.tested;

        bool pass = (c1.E >= Emin && c2.E >= Emin) &&
                    (c1.chi <= chiMx && c2.chi <= chiMx) &&
                    (asym   <= aMx);

        if (!pass)
        {
          if      (c1.E < Emin || c2.E < Emin) ++st.failE;
          else if (c1.chi > chiMx || c2.chi > chiMx) ++st.failChi;
          else                                     ++st.failAsy;
          continue;
        }

        const float mInv = (c1.v + c2.v).M();
        fillBoth(invKey(-1,-1,Emin,chiMx,aMx), trig, mInv);
        ++st.passed;
      }

      /* ---- pT‑binned spectra -------------------------------------- */
      for (const auto& pb : m_ptBins)
      {
        const float ptLo = pb.first, ptHi = pb.second;
        if (c1.pt < ptLo || c1.pt >= ptHi ||
            c2.pt < ptLo || c2.pt >= ptHi) continue;

        for (float Emin : m_minClusE)
        for (float chiMx: m_chi2Cuts)
        for (float aMx  : m_asymCuts)
        {
          const std::string key = statKey(ptLo,ptHi,Emin,chiMx,aMx);
          auto& st = evtStatLocal[key]; ++st.tested;

          bool pass = (c1.E >= Emin && c2.E >= Emin) &&
                      (c1.chi <= chiMx && c2.chi <= chiMx) &&
                      (asym   <= aMx);

          if (!pass)
          {
            if      (c1.E < Emin || c2.E < Emin) ++st.failE;
            else if (c1.chi > chiMx || c2.chi > chiMx) ++st.failChi;
            else                                     ++st.failAsy;
            continue;
          }

          const float mInv = (c1.v + c2.v).M();
          fillBoth(invKey(ptLo,ptHi,Emin,chiMx,aMx), trig, mInv);
          ++st.passed;
        }
      }

      /* ---- periodic progress ticker -------------------------------- */
      if ((++pairCnt % reportEvery) == 0)
      {
#ifdef _OPENMP
        if (omp_get_thread_num() == 0)
#endif
          std::cout << CLR_CYAN << "    [doPi0QA] processed "
                    << pairCnt << " / " << totalPairs << " pairs ("
                    << std::fixed << std::setprecision(1)
                    << 100.0 * pairCnt / totalPairs << "%)\r"
                    << CLR_RESET << std::flush;
      }
    } // end pair loop

    /* ---- merge local stats ----------------------------------------- */
#ifdef _OPENMP
    #pragma omp critical
#endif
    {
      for (const auto& kv : evtStatLocal)
      {
        m_evtStat[kv.first].tested  += kv.second.tested;
        m_evtStat[kv.first].failE   += kv.second.failE;
        m_evtStat[kv.first].failChi += kv.second.failChi;
        m_evtStat[kv.first].failAsy += kv.second.failAsy;
        m_evtStat[kv.first].passed  += kv.second.passed;
      }
    }
  } // end parallel region

  /* ---------- summary & exit ---------------------------------------- */
  const std::size_t processedPairs = pairCnt.load();
  LOG(3, CLR_GREEN,
      "    [doPi0QA] completed – "
      << processedPairs << " / " << totalPairs
      << " pairs processed (" << std::fixed << std::setprecision(1)
      << (100.0 * processedPairs / totalPairs) << "%)"
      << " – centrality slice " << lo << "–" << hi << '%');
}


/* ----------------------------------------------------------------------
 * fillCorrelations – ΣE / ΣQ detector‑level correlations
 *                    (writes to global‑AND‑centrality histograms)
 * -------------------------------------------------------------------- */
void emcal_sepdCorrelator::fillCorrelations(const std::vector<std::string>& trig)
{
  /* —— 0. Basic book‑keeping ———————————————————————————————— */
  static std::uint64_t callCtr = 0;   // stays local to this TU
  ++callCtr;
  LOG(3, CLR_BLUE, "[fillCorrelations] call #" << callCtr
                          << "  (triggers=" << trig.size() << ")  BEGIN");

  const double cemc  = m_calo["CEMC" ].sumE;
  const double ihcal = m_calo["IHCAL"].sumE;
  const double ohcal = m_calo["OHCAL"].sumE;

  /* —— 1. Determine centrality slice this event belongs to ———————— */
  int lo = 0, hi = 100;                      // default = ‘all events’
  if (m_centBin >= 0)
  {
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
      { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }
  }
  const std::string tag = '_' + std::to_string(lo) + '_' + std::to_string(hi);

  LOG(4, CLR_CYAN, "    centrality bin = " << m_centBin
                       << "  → slice " << lo << "–" << hi << " %");

  /* —— 2. Auxiliary lambda that *safely* writes to TH2F ——————————— */
  auto safeFill = [&](auto* obj, double x, double y,
                      const std::string& hName, const std::string& trg)
  {
    if (!obj)
    {
      LOG(1, CLR_YELLOW, "      [WARN] histogram '" << hName
                               << "' missing for trigger '" << trg << '\'');
      return false;
    }
    static_cast<TH2F*>(obj)->Fill(x, y);
    return true;
  };

  /* counters for diagnostic summary */
  std::unordered_map<std::string, std::size_t> binsFilled;

  /* —— 3. Loop over every *active* trigger ———————————————— */
  for (const auto& t : trig)
  {
    auto& H = qaHistogramsByTrigger[t];

    /* — 3.1 arm‑matched maps —————————————————————— */
    binsFilled[t] += safeFill(H["h_SEPD_S_vs_CEMC_South"], m_sepdQ_arm[0],
                              m_cemcEt_arm[0],
                              "h_SEPD_S_vs_CEMC_South", t);
    binsFilled[t] += safeFill(H["h_SEPD_N_vs_CEMC_North"], m_sepdQ_arm[1],
                              m_cemcEt_arm[1],
                              "h_SEPD_N_vs_CEMC_North", t);

    /* — 3.2 global (all events) maps ————————————— */
    binsFilled[t] += safeFill(H["h_SEPD_vs_CEMC" ], m_sepdQ, cemc ,
                              "h_SEPD_vs_CEMC", t);
    binsFilled[t] += safeFill(H["h_SEPD_vs_IHCAL"], m_sepdQ, ihcal,
                              "h_SEPD_vs_IHCAL", t);
    binsFilled[t] += safeFill(H["h_SEPD_vs_OHCAL"], m_sepdQ, ohcal,
                              "h_SEPD_vs_OHCAL", t);
    binsFilled[t] += safeFill(H["h_SEPD_vs_MBD"  ], m_sepdQ, m_mbdQ,
                              "h_SEPD_vs_MBD", t);

    binsFilled[t] += safeFill(H["h_MBD_vs_CEMC" ], m_mbdQ, cemc ,
                              "h_MBD_vs_CEMC", t);
    binsFilled[t] += safeFill(H["h_MBD_vs_IHCAL"], m_mbdQ, ihcal,
                              "h_MBD_vs_IHCAL", t);
    binsFilled[t] += safeFill(H["h_MBD_vs_OHCAL"], m_mbdQ, ohcal,
                              "h_MBD_vs_OHCAL", t);

    /* — 3.3 centrality‑tagged clones ———————————— */
    auto tryCent = [&](const std::string& base,
                       double x, double y)
    {
      const std::string hkey = base + tag + '_' + t;
      auto it = H.find(hkey);
      if (it == H.end()) return false;
      static_cast<TH2F*>(it->second)->Fill(x, y);
      return true;
    };

    binsFilled[t] += tryCent("h_SEPD_vs_CEMC",  m_sepdQ, cemc);
    binsFilled[t] += tryCent("h_SEPD_vs_IHCAL", m_sepdQ, ihcal);
    binsFilled[t] += tryCent("h_SEPD_vs_OHCAL", m_sepdQ, ohcal);
    binsFilled[t] += tryCent("h_SEPD_vs_MBD",   m_sepdQ, m_mbdQ);

    binsFilled[t] += tryCent("h_MBD_vs_CEMC",   m_mbdQ, cemc);
    binsFilled[t] += tryCent("h_MBD_vs_IHCAL",  m_mbdQ, ihcal);
    binsFilled[t] += tryCent("h_MBD_vs_OHCAL",  m_mbdQ, ohcal);

    binsFilled[t] += tryCent("h_SEPD_S_vs_CEMC_South", m_sepdQ_arm[0],
                             m_cemcEt_arm[0]);
    binsFilled[t] += tryCent("h_SEPD_N_vs_CEMC_North", m_sepdQ_arm[1],
                             m_cemcEt_arm[1]);
  } // trigger loop

  /* —— 4. Human‑readable one‑line summary ——————————————————— */
  LOG(3, CLR_GREEN,
      "    ΣE(CEMC)=" << cemc
      << "  ΣE(IHCAL)=" << ihcal
      << "  ΣE(OHCAL)=" << ohcal
      << "  ΣQ(MBD)="   << m_mbdQ
      << "  ΣQ(sEPD)="  << m_sepdQ
      << "  slice="     << lo << "–" << hi << " %");

  if (Verbosity() >= 4)
  {
    for (const auto& [trg, n] : binsFilled)
      LOG(4, CLR_CYAN, "      trigger '" << trg
                         << "': filled " << n << " correlation‑bins");
  }

  LOG(3, CLR_BLUE, "[fillCorrelations] call #" << callCtr << "  END");
}



//==========================================================================
//  trivial helpers
//==========================================================================
std::string emcal_sepdCorrelator::invKey(float ptLo, float ptHi,
                                         float e, float chi, float a)
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1)
      << "mInv_pt" << ptLo << "to" << ptHi
      << "_E" << e << "_chi" << chi << "_asy" << a;
  return oss.str();
}

std::string emcal_sepdCorrelator::statKey(float ptLo, float ptHi,
                                          float Emin, float chi, float a)
{
  std::ostringstream o; o << std::fixed << std::setprecision(1);
  if (ptLo < 0) o << "allPt";
  else          o << "pt" << ptLo << "to" << ptHi;
  o << "_E" << Emin << "_chi" << chi << "_asy" << a;
  return o.str();
}


// ----------------------------------------------------------------------
// Return the highest transverse energy of the jets in a container
// ----------------------------------------------------------------------
float
emcal_sepdCorrelator::getMaxJetEt(JetContainer* jets) const
{
  if (!jets) return 0.f;

  float maxEt = 0.f;
  for (const Jet* jet : *jets)               // ← iterate directly
  {
    if (!jet) continue;
    const float et = static_cast<float>(jet->get_et());
    if (et > maxEt) maxEt = et;
  }
  return maxEt;
}


// ----------------------------------------------------------------------
//  doJetQA – fill jet QA histograms
//            (max‑E_T, shape, lead–sub  +  v_n^{jet})
// ----------------------------------------------------------------------
int
emcal_sepdCorrelator::doJetQA(PHCompositeNode*                topNode,
                              const std::vector<std::string>& trig)
{
  LOG(4, CLR_BLUE, "  [doJetQA] – entering");

  /* ------------------------------------------------------------------ */
  /* 0. Centrality slice tag & reference planes                         */
  /* ------------------------------------------------------------------ */
  int lo = 0, hi = 100;
  if (m_centBin >= 0)
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
      { lo = m_centEdges[i];  hi = m_centEdges[i + 1]; break; }

  const std::string tag = '_' + std::to_string(lo) + '_' + std::to_string(hi);

  /* reference angles from sEPD South (index = harmonic) */
  const double psi[4] = {0., m_psi1_S, m_psi2_S, m_psi3_S};

  /* ------------------------------------------------------------------ */
  /* 1.  Scan every configured jet radius                               */
  /* ------------------------------------------------------------------ */
  std::unordered_map<std::string, float>                     maxJetEt;
  std::unordered_map<std::string, std::array<const Jet*, 2>> topTwoJets;

  for (const auto& [radKey, nodeName] : kJetRadii)
  {
    LOG(5, CLR_BLUE, "    checking jet node '" << nodeName << "' (" << radKey << ')');

    JetContainer* jets = findNode::getClass<JetContainer>(topNode, nodeName);
    if (!jets)
    {
      LOG(4, CLR_YELLOW, "    missing JetContainer '" << nodeName << "' – ABORTRUN");
      return Fun4AllReturnCodes::ABORTRUN;
    }

    const Jet* lead = nullptr;
    const Jet* sub  = nullptr;

    /* ---------- loop over *all* jets:                                *
     *            – select lead/sub for QA shapes                       *
     *            – apply |η|>3 gap and fill v_n profiles               */
    for (const Jet* j : *jets)
    {
      if (!j) continue;

      /* keep track of two highest‑E_T jets (no η gap here) */
      if (!lead || j->get_et() > lead->get_et()) { sub = lead;  lead = j; }
      else if (!sub  || j->get_et() > sub ->get_et()) { sub = j; }

      /* === jet flow ================================================ */
      if (std::abs(j->get_eta()) < 3.0) continue;   // η‑gap vs. SEPD

      const double phi = j->get_phi();
      const double pt  = j->get_pt();
      double vn[4] = {0.};
      for (int n : {1, 2, 3})
        vn[n] = std::cos(n * (phi - psi[n]));

      /* fill one profile per active trigger ------------------------- */
      for (const auto& t : trig)
      {
        auto& H = qaHistogramsByTrigger[t];

        std::ostringstream k1, k2, k3;
        k1 << "p_v1_JET_" << radKey << tag << '_' << t;
        k2 << "p_v2_JET_" << radKey << tag << '_' << t;
        k3 << "p_v3_JET_" << radKey << tag << '_' << t;

        if (auto it = H.find(k1.str()); it != H.end())
          static_cast<TProfile*>(it->second)->Fill(pt, vn[1]);
        if (auto it = H.find(k2.str()); it != H.end())
          static_cast<TProfile*>(it->second)->Fill(pt, vn[2]);
        if (auto it = H.find(k3.str()); it != H.end())
          static_cast<TProfile*>(it->second)->Fill(pt, vn[3]);
      }
    } // jet loop

    const float etMax = lead ? static_cast<float>(lead->get_et()) : 0.f;
    maxJetEt  [radKey] = etMax;
    topTwoJets[radKey] = {lead, sub};

    LOG(5, CLR_GREEN, "      max E_T = " << etMax << " GeV");
  } // radii loop (jet scan)

  /* ------------------------------------------------------------------ */
  /* 2.  Fill shape / max‑E_T histograms (uses topTwoJets)              */
  /* ------------------------------------------------------------------ */
  for (const auto& [radKey, etMax] : maxJetEt)
  {
    const std::string base1 = "h_maxJetEt_"          + radKey;
    const std::string base2 = "h_leadEt_vs_subEt_"   + radKey;
    const std::string base3 = "h_jetEt_area_nConst_" + radKey;

    for (const auto& t : trig)
    {
      auto& H = qaHistogramsByTrigger[t];

      /* ---- 2.1 1‑D max‑E_T ---------------------------------------- */
      if (auto* hG = dynamic_cast<TH1F*>(H[base1 + "_" + t])) hG->Fill(etMax);
      if (auto it = H.find(base1 + tag + "_" + t); it != H.end())
        static_cast<TH1F*>(it->second)->Fill(etMax);

      /* ---- 2.2 jet‑shape histograms ------------------------------- */
      const auto& jets = topTwoJets[radKey];

      for (const Jet* j : jets)                       // 3‑D shape
      {
        if (!j) continue;
        static_cast<TH3F*>(H[base3 + "_" + t])
            ->Fill(j->get_et(),
                   j->get_property(Jet::PROPERTY::prop_area),
                   static_cast<int>(j->size_comp()));

        if (auto it3 = H.find(base3 + tag + "_" + t); it3 != H.end())
          static_cast<TH3F*>(it3->second)
              ->Fill(j->get_et(),
                     j->get_property(Jet::PROPERTY::prop_area),
                     static_cast<int>(j->size_comp()));
      }

      /* 2‑D lead‑vs‑sub: need both jets ----------------------------- */
      if (jets[0] && jets[1])
      {
        static_cast<TH2F*>(H[base2 + "_" + t])
            ->Fill(jets[0]->get_et(), jets[1]->get_et());

        if (auto it2 = H.find(base2 + tag + "_" + t); it2 != H.end())
          static_cast<TH2F*>(it2->second)
              ->Fill(jets[0]->get_et(), jets[1]->get_et());
      }
    } // trigger loop
  }   // radii loop (hist fill)

  LOG(4, CLR_GREEN, "  [doJetQA] – finished OK");
  return Fun4AllReturnCodes::EVENT_OK;
}


void emcal_sepdCorrelator::fillFlowHists(const std::vector<std::string>& trig)
{
  /* centrality slice -------------------------------------------------- */
  int lo=0, hi=100;
  if (m_centBin>=0)
    for (std::size_t i=0;i+1<m_centEdges.size();++i)
      if (m_centBin>=m_centEdges[i] && m_centBin<m_centEdges[i+1])
      { lo=m_centEdges[i]; hi=m_centEdges[i+1]; break; }

  const double cPsi2 = std::cos(2*m_psi2_S), sPsi2 = std::sin(2*m_psi2_S);
  const double cPsi3 = std::cos(3*m_psi3_S), sPsi3 = std::sin(3*m_psi3_S);

  for (const auto& [det, vec] : m_flowAcc)
  {
    for (std::size_t ib=0; ib<vec.size(); ++ib)
    {
      const auto& a = vec[ib];
      if (a.sumW<=0) continue;

      const double v2 = (a.qx[2]*cPsi2 + a.qy[2]*sPsi2)/a.sumW;
      const double v3 = (a.qx[3]*cPsi3 + a.qy[3]*sPsi3)/a.sumW;

      for (const auto& t : trig)
      {
        auto& H = qaHistogramsByTrigger[t];
        std::ostringstream k2,k3;
        k2<<"p_v2_"<<det<<'_'<<lo<<'_'<<hi<<'_'<<t;
        k3<<"p_v3_"<<det<<'_'<<lo<<'_'<<hi<<'_'<<t;
        static_cast<TProfile*>(H[k2.str()])->Fill(int(ib)+0.5, v2);
        static_cast<TProfile*>(H[k3.str()])->Fill(int(ib)+0.5, v3);
      }
    }
  }
}


//==========================================================================
//  ResetEvent / Reset / End – unchanged logic, single definition
//==========================================================================
int emcal_sepdCorrelator::ResetEvent(PHCompositeNode*)
{
  /* wipe every per‑event cache so NOTHING bleeds into the next event */
  std::fill(std::begin(m_sepdQ_arm),  std::end(m_sepdQ_arm),  0.);
  std::fill(std::begin(m_mbdQ_arm),   std::end(m_mbdQ_arm),   0.);
  std::fill(std::begin(m_cemcEt_arm), std::end(m_cemcEt_arm), 0.);
  std::fill(std::begin(m_ihcalEt_arm),std::end(m_ihcalEt_arm),0.);
  std::fill(std::begin(m_ohcalEt_arm),std::end(m_ohcalEt_arm),0.);

  m_sepdQ  = 0.;
  m_mbdQ   = 0.;
  m_psi1_N = 0.;
  m_psi1_S = 0.;
  m_psi2_N = 0.;
  m_psi2_S = 0.;
  m_psi3_N = 0.;
  m_psi3_S = 0.;
  m_centBin = -1;
  for (auto& kv : m_calo) kv.second.sumE = 0.;

  m_evtStat.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

int emcal_sepdCorrelator::Reset(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

// --------------------------------------------------------------------------
//  End – enhanced diagnostics, robust against dangling pointers
// --------------------------------------------------------------------------
int emcal_sepdCorrelator::End(PHCompositeNode*)
{
  auto warn = [&](const std::string& m)
  { std::cerr << CLR_YELLOW << "[End] " << m << CLR_RESET << '\n'; };

  auto info = [&](int lvl, const std::string& m)
  { LOG(lvl, CLR_GREEN, "[End] " << m); };

  //--------------------------------------------------------------------
  // 1. Basic checks on the output file pointer
  //--------------------------------------------------------------------
  if (!out)           { warn("TFile* 'out' is nullptr – nothing to write");  return Fun4AllReturnCodes::ABORTEVENT; }
  if (!out->IsOpen()) { warn("Output file is *not* open ("+std::string(out->GetName())+")"); return Fun4AllReturnCodes::ABORTEVENT; }

  //--------------------------------------------------------------------
  // 2. Write histograms trigger‑by‑trigger
  //--------------------------------------------------------------------
  std::size_t nHistExpected = 0, nHistWritten = 0;

  for (auto& [trig, hMap] : qaHistogramsByTrigger)
  {
    TDirectory* dir = out->GetDirectory(trig.c_str());
    if (!dir) dir = out->mkdir(trig.c_str());
    dir->cd();

    std::size_t nExpThis = 0, nWrtThis = 0;

    for (auto& [key, obj] : hMap)
    {
      ++nHistExpected; ++nExpThis;

      TH1* h = dynamic_cast<TH1*>(obj);
      if (!h)              { warn("Object '"+key+"' (trigger "+trig+") is not TH1 – skipped"); continue; }
      if (h->GetEntries()==0)
      { if (Verbosity()>1) warn("Histogram '"+key+"' (trigger "+trig+") has 0 entries – skipped"); continue; }

      try
      {
        if (h->Write("", TObject::kOverwrite) > 0) { ++nHistWritten; ++nWrtThis; }
        else warn("Write() returned 0 for '"+key+"' (trigger "+trig+")");
      }
      catch (const std::exception& e)
      { warn("Exception while writing '"+key+"' – "+std::string(e.what())); }
    }
    info(1, "trigger '"+trig+"': "+std::to_string(nWrtThis)+" / "
               +std::to_string(nExpThis)+" histograms written");
    out->cd();
  }

  //--------------------------------------------------------------------
  // 3.  Human‑readable summary  (must run *before* the file is deleted)
  //--------------------------------------------------------------------
  if (Verbosity() > 0)
  {
    std::cout << "\n\033[1mHistogram summary\033[0m\n"
              << "\033[1mTrigger                        │ Histogram                           │  Entries\033[0m\n"
              << "-------------------------------------------------------------------------------\n";

    for (const auto& [trig, hMap] : qaHistogramsByTrigger)
      for (const auto& [key, obj]  : hMap)
        if (const TH1* h = dynamic_cast<const TH1*>(obj))
        {
            if (h->GetEntries() == 0) continue;          // <‑‑ only list filled ones
            std::cout << std::left << std::setw(30) << trig << " │ "
            << std::setw(32) << key  << " │ "
            << std::right<< std::setw(10)
            << static_cast<Long64_t>(h->GetEntries()) << '\n';
        }
    if (!m_totStat.empty())
    {
      std::cout << "-------------------------------------------------------------------------------\n"
                << "\n\033[1mπ0‑QA cut summary (all events)\033[0m\n"
                << "\033[1mcut‑key                                    │ tested        passed    eff[%]\033[0m\n"
                << "--------------------------------------------------------------------------\n";
      for (const auto& [key, s] : m_totStat)
      {
        const double eff = s.tested ? 100.*s.passed / s.tested : 0.;
        std::cout << std::left  << std::setw(44) << key << " │ "
                  << std::right << std::setw(12) << s.tested
                  << std::setw(12) << s.passed
                  << std::setw(9)  << std::fixed << std::setprecision(2) << eff << '\n';
      }
      std::cout << "--------------------------------------------------------------------------\n";
    }
  }

  //--------------------------------------------------------------------
  // 4.  Write footer & close the file
  //--------------------------------------------------------------------
  if (Verbosity() >= 1)
      std::cout << "\nOutput ROOT file →  " << out->GetName() << "\n\n";
    
  info(1, "writing TFile footer and closing ("+std::to_string(nHistWritten)
           +" / "+std::to_string(nHistExpected)+" objects written)");

  try      { out->Write(); out->Close(); }
  catch (const std::exception& e)
  { warn("Exception during TFile::Write/Close – "+std::string(e.what())); }

  delete out; out = nullptr;          // safe: we no longer dereference histos
  info(0, "Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}


void emcal_sepdCorrelator::Print(const std::string&) const
{ /* nothing to print beyond ROOT histograms */ }
