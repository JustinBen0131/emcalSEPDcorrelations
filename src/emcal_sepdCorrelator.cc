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
#include <TKey.h>
#include <CLHEP/Vector/ThreeVector.h>
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
#include <ffarawobjects/Gl1Packet.h>
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
#include <algorithm>   // std::clamp
#include <cmath>       // std::cosh, std::hypot, std::fmod
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <tuple>

#ifdef _OPENMP
  #include <omp.h>
#endif

//–––––––– helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#define CLR_BLUE   "\033[1;34m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_MAGENTA "\033[1;35m"
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

// ----------------------------------------------------------------------
//  getCentralitySlice  – returns true if the current event falls into
//                        one of the user–defined centrality intervals
//                        (lo ≤ cent < hi).  On success it also builds
//                        the “_lo_hi” tag that all QA routines use
//                        when addressing centrality‑specific histograms.
//
//  • lo,hi   – returned lower / upper edge (undefined if the function
//              returns false)
//  • tag     – empty if the function returns false
// ----------------------------------------------------------------------
bool emcal_sepdCorrelator::getCentralitySlice(int& lo,
                                              int& hi,
                                              std::string& tag) const
{
  lo = hi = -1;
  for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
    { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }

  const bool hasSlice = (lo >= 0);
  tag.clear();
  if (hasSlice)
  {
    std::ostringstream ss;
    ss << '_' << lo << '_' << hi;
    tag = ss.str();
  }
  return hasSlice;
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
    

  //---------------------------------------------------------------------------
  //  Correlation map between Minimum‑Bias decision and trigger status
  //  ─────────────────────────────────────────────────────────────────────────
  //
  //  Y‑axis categories (Event category):
  //    1  →  ¬MB  &  ¬Trig   — event fails both cuts
  //    2  →  ¬MB  &   Trig   — event rejected **only** by MB cut
  //    3  →   MB  &  ¬Trig   — event rejected **only** by trigger logic
  //    4  →   MB  &   Trig   — event accepted by *both* cuts
  //---------------------------------------------------------------------------
  {
        out->mkdir("triggerQA")->cd();                      // separate folder
        const int nTrig = triggerNameMap.size();

        h_MBTrigCorr = new TH2I("h_MB_vs_Trigger",
                                "Minimum‑bias vs Trigger decision;"
                                "Trigger key;"
                                "Event category",
                                nTrig, 0.5, nTrig + 0.5,          // X: trigger
                                4,      0.5, 4.5);                // Y: category

        // Label Y‑bins with the categories explained above
        h_MBTrigCorr->GetYaxis()->SetBinLabel(1, "¬MB  ¬Trig");
        h_MBTrigCorr->GetYaxis()->SetBinLabel(2, "¬MB   Trig");
        h_MBTrigCorr->GetYaxis()->SetBinLabel(3, " MB  ¬Trig");
        h_MBTrigCorr->GetYaxis()->SetBinLabel(4, " MB   Trig");

        // Label X‑bins with trigger keys and cache their indices
        int ib = 1;
        for (auto& [bitName, key] : triggerNameMap)
        {
            h_MBTrigCorr->GetXaxis()->SetBinLabel(ib, key.c_str());
            m_trigBin[key] = ib++;                      // cache index
        }
        out->cd();                                     // back to root
  }


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
        const std::vector<std::string> dets = {
          "CEMC_S","CEMC_N","CEMC_T",
          "IHCAL_S","IHCAL_N","IHCAL_T",
          "OHCAL_S","OHCAL_N","OHCAL_T",
          "HCAL_S","HCAL_N","HCAL_T",
          "ALL_S","ALL_N","ALL_T"
        };
        for (const auto& d : dets)
            m_flowAcc[d].assign(m_ptBins.size(), {});       // one FlowAcc per pT‑bin
  }
  LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
  return Fun4AllReturnCodes::EVENT_OK;
}

// Return the vertex‑cut contained in “…_vtx_lt_<N>”
// ‑1  ➜  no explicit cut in the name.
static int
extractVtxCut(const std::string& trigName)
{
  static const std::regex re(R"(vtx_lt_(\d+))");
  std::smatch m;
  return std::regex_search(trigName, m, re) ? std::stoi(m[1]) : -1;
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
      
    H["h_SEPD_RingOcc_South_" + trig] =
          new TH1F(("h_SEPD_RingOcc_South_" + trig).c_str(),
                   "sEPD South – ring occupancy;ring index;hits / event",
                   16, -0.5, 15.5);

    H["h_SEPD_RingOcc_North_" + trig] =
          new TH1F(("h_SEPD_RingOcc_North_" + trig).c_str(),
                   "sEPD North – ring occupancy;ring index;hits / event",
                   16, -0.5, 15.5);

    H["h_SEPD_RingQ_South_" + trig] =
          new TH1F(("h_SEPD_RingQ_South_" + trig).c_str(),
                   "sEPD South – ΣQ per ring;ring index;ΣQ  [ADC]",
                   16, -0.5, 15.5);

    H["h_SEPD_RingQ_North_" + trig] =
          new TH1F(("h_SEPD_RingQ_North_" + trig).c_str(),
                   "sEPD North – ΣQ per ring;ring index;ΣQ  [ADC]",
                   16, -0.5, 15.5);

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
      
    /* generic centrality‑clone helper – works for TH2F, TH2Poly, … */
    auto cloneHitMap = [&](const std::string& base, TObject* src)
      {
        for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
        {
          const int lo = m_centEdges[i], hi = m_centEdges[i + 1];

          std::ostringstream n;
          n << base << '_' << lo << '_' << hi << '_' << trig;   // final key

          TObject* c = src->Clone(n.str().c_str());             // deep copy
          if (auto* h = dynamic_cast<TH1*>(c))
          {
            h->Reset();                 // start empty
            h->SetDirectory(nullptr);   // detach from current directory
          }
          H[n.str()] = c;               // register in map
        }
    };

    /* ── calorimeters  ─────────────────────────── */
    cloneHitMap("h_EMC_EtaPhiMap",   H["h_EMC_EtaPhiMap_"   + trig]);
    cloneHitMap("h_IHCAL_EtaPhiMap", H["h_IHCAL_EtaPhiMap_" + trig]);
    cloneHitMap("h_OHCAL_EtaPhiMap", H["h_OHCAL_EtaPhiMap_" + trig]);
    cloneHitMap("h_MBD_Hitmap_South", H["h_MBD_Hitmap_South_" + trig]);
    cloneHitMap("h_MBD_Hitmap_North", H["h_MBD_Hitmap_North_" + trig]);
    cloneHitMap("h_sEPD_Hitmap_South", H["h_sEPD_Hitmap_South_" + trig]);
    cloneHitMap("h_sEPD_Hitmap_North", H["h_sEPD_Hitmap_North_" + trig]);
    cloneHitMap("h_SEPD_RingOcc_South", H["h_SEPD_RingOcc_South_" + trig]);
    cloneHitMap("h_SEPD_RingOcc_North", H["h_SEPD_RingOcc_North_" + trig]);
    cloneHitMap("h_SEPD_RingQ_South",   H["h_SEPD_RingQ_South_"   + trig]);
    cloneHitMap("h_SEPD_RingQ_North",   H["h_SEPD_RingQ_North_"   + trig]);
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
    {
        /* per‑cluster spectrum (already there) */
        H["h_clusterE_EMC"] =
            new TH1F(("h_clusterE_EMC_" + trig).c_str(),
                     "EMC cluster E;E [GeV]", nbE, 0, eMax);

        /* >>> NEW – per‑event maximum‑cluster‑energy <<< */
        H["h_maxClusterE_EMC"] =
            new TH1F(("h_maxClusterE_EMC_" + trig).c_str(),
                     "max EMC cluster E per event;E_{max} [GeV];Events",
                     nbE, 0, eMax);

        H["h_maxClusterEnergy_doNotScale_" + trig] =
            new TH1F(("h_maxClusterEnergy_doNotScale_" + trig).c_str(),
                     "max EMC cluster E per event (doNotScale);E_{max} [GeV];Events",
                     nbE, 0, eMax);

        /* centrality‑tagged clones ------------------------------------ */
        for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
        {
          const int lo = m_centEdges[i]   ;
          const int hi = m_centEdges[i+1] ;

          std::ostringstream n;
          n << "h_maxClusterE_EMC_" << lo << '_' << hi << '_' << trig;

          H[n.str()] = new TH1F(n.str().c_str(),
                                "max EMC cluster E per event;E_{max} [GeV];Events",
                                nbE, 0, eMax);
        }
     }
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
    
  H["h_IHCAL_vs_CEMC"] = book2(("h_IHCAL_vs_CEMC_" + trig).c_str(),
                                 "IHCAL #SigmaE vs CEMC #SigmaE",
                                 nE, 0, eMax,   nE, 0, eMax);

  H["h_OHCAL_vs_CEMC"] = book2(("h_OHCAL_vs_CEMC_" + trig).c_str(),
                                 "OHCAL #SigmaE vs CEMC #SigmaE",
                                 nE, 0, eMax,   nE, 0, eMax);
    
    
  H["h_dEta_CEMC_IHCAL"] = new TH1F(("h_dEta_CEMC_IHCAL_" + trig).c_str(),
                                      "#Delta#eta (IHCAL – CEMC);#Delta#eta;Events",
                                      120, -6.0, 6.0);     // 0.1‑wide bins
  H["h_dPhi_CEMC_IHCAL"] = new TH1F(("h_dPhi_CEMC_IHCAL_" + trig).c_str(),
                                      "#Delta#phi (IHCAL – CEMC);#Delta#phi;Events",
                                      128, -TMath::Pi(), TMath::Pi());   // 5° bins

  /* --- ΣQ(sEPD South) × ΣQ(sEPD North) ----------------------- */
  H["h_SEPD_S_vs_SEPD_N"] =
        book2(("h_SEPD_S_vs_SEPD_N_" + trig).c_str(),
              "#SigmaQ_{sEPD South}  vs  #SigmaQ_{sEPD North};"
              "#SigmaQ_{South} [ADC];#SigmaQ_{North} [ADC]",
              nC, 0, cMax,   // X‑axis = South arm
              nC, 0, cMax);  // Y‑axis = North arm
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
    
  /* helper that clones a TH1F for a given centrality slice --------------- */
  auto clone1D = [&](const std::string& base,
                       int lo, int hi,
                       const char*   title,
                       int nb, double loX, double hiX)
    {
        std::ostringstream key;           // final histogram key
        key << base << '_' << lo << '_' << hi << '_' << trig;
        H[key.str()] = new TH1F(key.str().c_str(), title, nb, loX, hiX);

        /* allow quick lookup at fill time (value is unused) */
        std::ostringstream tag; tag << '_' << lo << '_' << hi;
        m_centIdxCache[tag.str()] = 1;
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
    addClone("h_SEPD_S_vs_SEPD_N", lo,hi,
               "#SigmaQ_{sEPD South}  vs  #SigmaQ_{sEPD North}",
               nC,0,cMax, nC,0,cMax);
      
    addClone("h_IHCAL_vs_CEMC", lo,hi,
               "IHCAL #SigmaE vs CEMC #SigmaE",
               nE,0,eMax, nE,0,eMax);
    addClone("h_OHCAL_vs_CEMC", lo,hi,
               "OHCAL #SigmaE vs CEMC #SigmaE",
               nE,0,eMax, nE,0,eMax);

    clone1D("h_dEta_CEMC_IHCAL", lo,hi,
              "#Delta#eta (IHCAL – CEMC);#Delta#eta;Events",
              120,-6.0,6.0);

    clone1D("h_dPhi_CEMC_IHCAL", lo,hi,
              "#Delta#phi (IHCAL – CEMC);#Delta#phi;Events",
              128,-TMath::Pi(),TMath::Pi());
  }
}


/* ----------------------------------------------------------------------
 * bookPi0MassSpectra – π0 invariant‑mass spectra
 *                      (global  +  centrality‑tagged clones)
 * -------------------------------------------------------------------- */
void emcal_sepdCorrelator::bookPi0MassSpectra(const std::string& trig,
                                              HistMap&           H)
{
  const int    nM   = 150;
  const double mMax = 1.5;    // [GeV/c²]

  /* ---- helper for the existing 1‑D spectra ------------------------------ */
  auto addHist1D = [&](const std::string& baseKey)
  {
    /* 1) centrality‑independent */
    const std::string hNameGlobal = baseKey + "_" + trig;
    H[hNameGlobal] = new TH1F(hNameGlobal.c_str(),
                              "m_{#gamma#gamma};GeV/c^{2}",
                              nM, 0., mMax);

    /* 2) centrality‑tagged clones */
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

  /* ---------- NEW ► helper for the three 2‑D correlation plots ---------- */
  auto addHist2D = [&](const std::string& baseKey,
                       int    nx, double xmin, double xmax,
                       int    ny, double ymin, double ymax,
                       const  char* xTitle,
                       const  char* yTitle)
  {
    /* global */
    const std::string gName = baseKey + "_" + trig;
    H[gName] = new TH2F(gName.c_str(),
                        (std::string(xTitle) + ";" + yTitle).c_str(),
                        nx, xmin, xmax,
                        ny, ymin, ymax);

    /* centrality‑tagged */
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    {
      const int lo = m_centEdges[i];
      const int hi = m_centEdges[i + 1];

      std::ostringstream name;
      name << baseKey << '_' << lo << '_' << hi << '_' << trig;

      H[name.str()] = new TH2F(name.str().c_str(),
                               (std::string(xTitle) + ";" + yTitle).c_str(),
                               nx, xmin, xmax,
                               ny, ymin, ymax);
    }
  };
  /* --------------------------------------------------------------------- */

  /* ---- pT‑binned 1‑D spectra ------------------------------------------ */
  for (auto pt : m_ptBins)
    for (float Emin : m_minClusE)
      for (float chi : m_chi2Cuts)
        for (float a : m_asymCuts)
          addHist1D(invKey(pt.first, pt.second, Emin, chi, a));

  /* ---- inclusive 1‑D spectra ------------------------------------------ */
  for (float Emin : m_minClusE)
    for (float chi : m_chi2Cuts)
      for (float a : m_asymCuts)
        addHist1D(invKey(-1, -1, Emin, chi, a));      // pT = −1 sentinel

  /* ---------- NEW ► un‑cut 2‑D correlation plots (once per trigger) ---- */
  addHist2D("Minv_vs_Asym",
            50,  0.0, 1.0,     // |E1−E2|/(E1+E2)
            120, 0.0, 0.6,     // mInv
            "|E_{1}-E_{2}|/(E_{1}+E_{2})",
            "m_{#gamma#gamma} (GeV/c^{2})");

  addHist2D("Minv_vs_chi2",
            60,  0.0, 6.0,     // (χ²₁+χ²₂)/2
            120, 0.0, 0.6,
            "(#chi^{2}_{1}+ #chi^{2}_{2})/2",
            "m_{#gamma#gamma} (GeV/c^{2})");

  addHist2D("Minv_vs_Eavg",
            120, 0.0, 20.0,    // (E1+E2)/2
            120, 0.0, 0.6,
            "(E_{1}+E_{2})/2  (GeV)",
            "m_{#gamma#gamma} (GeV/c^{2})");
}


void emcal_sepdCorrelator::bookFlowQA(const std::string& trig, HistMap& H)
{
    /* one TProfile per {detector, harmonic, centrality} filled versus *real* pT */
    std::vector<double> ptEdge;                         // lower‑edge array
    ptEdge.reserve(m_ptBins.size()+1);
    ptEdge.push_back(m_ptBins.front().first);
    for (auto& b : m_ptBins) ptEdge.push_back(b.second);

    auto make = [&](const std::string& det, int n, int lo, int hi)
    {
      std::ostringstream name;
      name << "p_v" << n << '_' << det << '_' << lo << '_' << hi << '_' << trig;
      auto* p = new TProfile(name.str().c_str(),
                             Form("v_{%d} (%s);p_{T}^{tower} [GeV];v_{%d}",
                                  n,det.c_str(),n),
                             ptEdge.size()-1, &ptEdge[0], "s");
      p->SetStats(0);
      H[name.str()] = p;
    };

    /* detectors:  South, North, Total (South+North) */
    static const std::array<const char*,15> detKeys = {
        "CEMC_S","CEMC_N","CEMC_T",
        "IHCAL_S","IHCAL_N","IHCAL_T",
        "OHCAL_S","OHCAL_N","OHCAL_T",
        "HCAL_S","HCAL_N","HCAL_T",
        "ALL_S","ALL_N","ALL_T"
    };

    /* centrality‑binned profiles */
    for (int n : {1,2,3})
      for (const auto& det : detKeys)
        for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
          make(det, n, m_centEdges[i], m_centEdges[i + 1]);

    /* minimum‑bias (0–100 %) profiles */
    for (int n : {1,2,3})
      for (const auto& det : detKeys)
        make(det, n, 0, 100);

    /* ----------  store ⟨cos n ΔΨ⟩ (n = 1,2,3) vs centrality -------------- */
    {
      std::vector<double> centEdgeD(m_centEdges.begin(), m_centEdges.end());

      for (int n : {1, 2, 3})
      {
        std::ostringstream key;   key  << "p_R" << n << "_vs_cent_" << trig;
        std::ostringstream titl;  titl << "cos ";
        if (n > 1) titl << n;
        titl << "(#Psi_{" << n << "}^{N}-#Psi_{" << n << "}^{S}) vs centrality;"
             << "centrality bin [%];#LT cos";
        if (n > 1) titl << ' ' << n;
        titl << "Δ#Psi #GT";

        auto* p = new TProfile(key.str().c_str(), titl.str().c_str(),
                               centEdgeD.size() - 1,
                               centEdgeD.data(), "s");
        p->SetStats(0);
        H[key.str()] = p;
      }
   }
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

  /* 3. ψ n distributions  (South & North, n = 1,2,3) ---------------- */
  H["h_Psi1_sEPD_S"] = new TH1F(("h_Psi1_sEPD_S_" + trig).c_str(),
                                  "sEPD South  #Psi_{1};#Psi_{1} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());
  H["h_Psi1_sEPD_N"] = new TH1F(("h_Psi1_sEPD_N_" + trig).c_str(),
                                  "sEPD North  #Psi_{1};#Psi_{1} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());

  H["h_Psi2_sEPD_S"] = new TH1F(("h_Psi2_sEPD_S_" + trig).c_str(),
                                  "sEPD South  #Psi_{2};#Psi_{2} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());
  H["h_Psi2_sEPD_N"] = new TH1F(("h_Psi2_sEPD_N_" + trig).c_str(),
                                  "sEPD North  #Psi_{2};#Psi_{2} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());

  H["h_Psi3_sEPD_S"] = new TH1F(("h_Psi3_sEPD_S_" + trig).c_str(),
                                  "sEPD South  #Psi_{3};#Psi_{3} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());
  H["h_Psi3_sEPD_N"] = new TH1F(("h_Psi3_sEPD_N_" + trig).c_str(),
                                  "sEPD North  #Psi_{3};#Psi_{3} [rad]",
                                  120, -TMath::Pi(), TMath::Pi());

  /* --- Q‑vector scatter plots for detailed debugging --------------- */
  for (auto& [det, vec] : m_flowAcc)                // every detector key
      for (int n : {1,2,3})
      {
          std::ostringstream nm;
          nm << "h_Q" << n << '_' << det << '_' << trig;
          H[nm.str()] = new TH2F(nm.str().c_str(),
                                 Form("Q_{%d} (%s);Q_{x};Q_{y}", n, det.c_str()),
                                 200, -500., 500.,
                                 200, -500., 500.);
  }


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
// Book jet QA spectra
//   (A) 1‑D   max‑E_T                            – h_maxJetEt_*
//   (B) 3‑D   E_T × area × Nconst                – h_jetEt_area_nConst_*
//   (C) 2‑D   leading‑ vs sub‑leading‑jet E_T    – h_leadEt_vs_subEt_*
//   (D) jet‑flow profiles  v_n^{jet}(pT) (n = 1,2,3) – p_v{n}_JET_*
// ----------------------------------------------------------------------
void emcal_sepdCorrelator::bookJetQA(const std::string& trig, HistMap& H)
{
  /* ---------- axis configuration ------------------------------------ */
  const int    nbEt  = 200;  const double etMax  = 200.;   // 1 GeV/bin
  const int    nbA   =  80;  const double aMax   = 0.80;   // ΔA = 0.01
  const int    nbN   = 200;  const double nMax   = 200.;   // 1 constituent/bin
  const int    nPtPr =  40;                                 // 5 GeV bins

  for (const auto& r : kJetRadii)              // r.first = "r02", …
  {
      const std::string base1 = "h_maxJetEt_"            + std::string(r.first);
      const std::string base2 = "h_leadEt_vs_subEt_"     + std::string(r.first);
      const std::string base3 = "h_jetEt_area_nConst_"   + std::string(r.first);
      const std::string base4 = "h_jetEtEtaPhi_"         + std::string(r.first);   // NEW

      /* ---------- global histograms ---------- */
      H[base1 + "_" + trig] = new TH1F(
          (base1 + "_" + trig).c_str(),
          ("max jet E_{T} (" + std::string(r.first) + ");E_{T} [GeV]").c_str(),
          nbEt, 0, etMax);

      H[base2 + "_" + trig] = new TH2F(
          (base2 + "_" + trig).c_str(),
          ("Leading vs sub‑leading jet E_{T} (" + std::string(r.first) +
           ");E_{T}^{lead} [GeV];E_{T}^{sub} [GeV]").c_str(),
          nbEt, 0, etMax, nbEt, 0, etMax);

      H[base3 + "_" + trig] = new TH3F(
          (base3 + "_" + trig).c_str(),
          ("Jet E_{T} vs area vs N_{const} (" + std::string(r.first) +
           ");E_{T} [GeV];Area;N_{const}").c_str(),
          nbEt, 0, etMax, nbA, 0, aMax, nbN, 0, nMax);

      /* NEW ►  E_{T} × η × φ */
      const int nbEta = 120;  const double etaMin = -6.0, etaMax = 6.0;
      const int nbPhi = 128;  const double phiMin = -TMath::Pi(), phiMax =  TMath::Pi();
      H[base4 + "_" + trig] = new TH3F(
          (base4 + "_" + trig).c_str(),
          ("Jet E_{T} vs #eta vs #phi (" + std::string(r.first) +
           ");#eta;#phi [rad];E_{T} [GeV]").c_str(),
          nbEta, etaMin, etaMax,
          nbPhi, phiMin, phiMax,
          nbEt , 0      , etMax);

      /* ----------------------------------------------------------------
       * (A) – (D)   ONE clone per user‑defined centrality slice
       * ----------------------------------------------------------------*/
      for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      {
        const int lo = m_centEdges[i];
        const int hi = m_centEdges[i + 1];

        std::ostringstream n1,n2,n3,n4;
        n1 << base1 << '_' << lo << '_' << hi << '_' << trig;
        n2 << base2 << '_' << lo << '_' << hi << '_' << trig;
        n3 << base3 << '_' << lo << '_' << hi << '_' << trig;
        n4 << base4 << '_' << lo << '_' << hi << '_' << trig;   // NEW

        H[n1.str()] = new TH1F(n1.str().c_str(),
            ("max jet E_{T} (" + std::string(r.first) + ");E_{T} [GeV]").c_str(),
            nbEt, 0, etMax);

        H[n2.str()] = new TH2F(n2.str().c_str(),
            ("Leading vs sub‑leading jet E_{T} (" + std::string(r.first) +
             ");E_{T}^{lead} [GeV];E_{T}^{sub} [GeV]").c_str(),
            nbEt, 0, etMax, nbEt, 0, etMax);

        H[n3.str()] = new TH3F(n3.str().c_str(),
            ("Jet E_{T} vs area vs N_{const} (" + std::string(r.first) +
             ");E_{T} [GeV];Area;N_{const}").c_str(),
            nbEt, 0, etMax, nbA, 0, aMax, nbN, 0, nMax);

        /* clone of the NEW histogram */
        H[n4.str()] = new TH3F(n4.str().c_str(),
            ("Jet E_{T} vs #eta vs #phi (" + std::string(r.first) +
             ");#eta;#phi [rad];E_{T} [GeV]").c_str(),
            nbEta, etaMin, etaMax,
            nbPhi, phiMin, phiMax,
            nbEt , 0      , etMax);
      }


    /* ----------------------------------------------------------------
     * (D)  v_n^{jet}(p_T)  – one profile per {n, cent‑slice, trigger}
     * ----------------------------------------------------------------*/
    for (int n : {1,2,3})
      for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      {
        const int lo = m_centEdges[i]   ;
        const int hi = m_centEdges[i+1] ;

        std::ostringstream nm;
        nm << "p_v" << n << "_JET_" << r.first << '_'
           << lo << '_' << hi << '_' << trig;

        H[nm.str()] = new TProfile(
            nm.str().c_str(),
            Form("v_{%d}^{jet} (%s);p_{T}^{jet} [GeV];v_{%d}",
                 n, r.first, n),
            nPtPr, 0, 200, "s");
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
      
    H["cnt_"+trig+"_raw"] =
          new TH1I(("cnt_"+trig+"_raw").c_str(),
                   (trig+" – raw bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);

    H["cnt_"+trig+"_live"] =
          new TH1I(("cnt_"+trig+"_live").c_str(),
                   (trig+" – live bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);

    H["cnt_"+trig+"_scaled"] =
          new TH1I(("cnt_"+trig+"_scaled").c_str(),
                   (trig+" – scaled bit fired;flag;Events").c_str(),
                   1, 0.5, 1.5);
      
      /* ── vertex‑cut compliance histogram ─────────────────────────────── */
    if (int vCut = extractVtxCut(trig); vCut > 0)
    {
        auto* h = new TH1I(Form("h_vtxRelToCut_%s", trig.c_str()),
                           Form("|z_{vtx}| vs trigger‑cut = %d cm;relation;Events", vCut),
                           3, 0.5, 3.5);
        h->GetXaxis()->SetBinLabel(1, "< cut");
        h->GetXaxis()->SetBinLabel(2, "= cut");
        h->GetXaxis()->SetBinLabel(3, "> cut");
        H[h->GetName()] = h;
    }
      
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

// ------------------------------------------------------------------
//  Turn a 64‑bit word into “{ 10,11,14 }” for human‑friendly printout
// ------------------------------------------------------------------
std::string
emcal_sepdCorrelator::bitsetToList(uint64_t word) const
{
  std::ostringstream os;
  os << '{';
  bool first = true;
  for (unsigned b = 0; b < 64; ++b)
    if ((word >> b) & 1U)
    {
      if (!first) os << ',';
      os << ' ' << b;
      first = false;
    }
  if (!first) os << ' ';
  os << '}';
  return os.str();
}

// ------------------------------------------------------------------
//  One‑line + (optionally) detailed trigger summary per event
// ------------------------------------------------------------------
void
emcal_sepdCorrelator::printTriggerSummary(const std::vector<std::string>& active,
                                          uint64_t wRaw, uint64_t wLive, uint64_t wScaled) const
{
  if (Verbosity() < 3) return;                  // show only when asked for it

  std::cout << CLR_CYAN
            << "Evt " << event_count
            << " | vtx = (" << std::fixed << std::setprecision(2)
            << m_vx << ',' << m_vy << ',' << m_vz << ") cm"
            << " | MB=" << (m_isMinBias ? "✓" : "✗")
            << " | scaled bits = " << bitsetToList(wScaled) << '\n';

  /* list accepted vs. rejected trigger keys -------------------------- */
  std::ostringstream acc, rej;
  for (auto& [bit, key] : triggerNameMap)
    ((std::find(active.begin(), active.end(), key) != active.end())
        ? acc : rej) << ' ' << key;

  std::cout << "        | accepted triggers:"
            << (acc.str().empty() ? " –" : acc.str()) << '\n'
            << "        | rejected triggers:"
            << (rej.str().empty() ? " –" : rej.str()) << CLR_RESET << std::endl;

  /* full decision table for heavy‑debugging -------------------------- */
  if (Verbosity() >= 6)
  {
    std::cout << "bit │ key                                   │ raw live scaled  decision\n"
                 "────┼───────────────────────────────────────┼────────────────────────\n";
    for (auto& [bit, key] : triggerNameMap)
    {
      const bool r = (wRaw    >> bit) & 1U;
      const bool l = (wLive   >> bit) & 1U;
      const bool s = (wScaled >> bit) & 1U;
      std::cout << std::setw(3) << bit << " │ "
                << std::left  << std::setw(31) << key << " │ "
                << (r ? "✓" : "✗") << "    "
                << (l ? "✓" : "✗") << "    "
                << (s ? "✓" : "✗") << "      "
                << (s ? "ACCEPT" : "reject") << '\n';
    }
    std::cout << std::endl;
  }
}

// ------------------------------------------------------------------
//  firstEventCuts – returns true iff
//        (Minimum‑bias  &&  ≥1 scaled bit fired  &&  |vz| < m_vzCut)
//  Also fills bookkeeping histograms and the per‑trigger vertex QA.
// ------------------------------------------------------------------
bool
emcal_sepdCorrelator::firstEventCuts(PHCompositeNode* topNode,
                                     std::vector<std::string>& activeTrig)
{
  /* 0. clear output vector ----------------------------------------- */
  activeTrig.clear();

  /* 1. minimum‑bias flag (cached earlier in fetchNodes) ------------ */
  const bool isMB = m_isMinBias;
  LOG(2, CLR_BLUE, "[firstEventCuts] event " << event_count
         << "  –  isMB = " << std::boolalpha << isMB);

  /* 2. obtain GL1 trigger packet ----------------------------------- */
  uint64_t wScaled = 0, wLive = 0, wRaw = 0;

  auto* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1)                     // fallback: DST name “14001”
      gl1 = findNode::getClass<Gl1Packet>(topNode, "14001");

  if (gl1)
  {
    wScaled = gl1->lValue(0, "ScaledVector");
    wLive   = gl1->lValue(0, "LiveVector");
    wRaw    = gl1->lValue(0, "TriggerVector");

    LOG(3, CLR_BLUE, "  GL1  raw=0x" << std::hex << wRaw
                     << "  live=0x"   << wLive
                     << "  scaled=0x" << wScaled << std::dec);
  }
  else
    LOG(1, CLR_YELLOW,
        "  GL1Packet node missing – assuming all trigger bits = 0");

  /* 3. decode **once** per word – reuse inside the loop ------------ */
  const auto bitsScaled = extractTriggerBits(wScaled, event_count);
  const auto bitsLive   = extractTriggerBits(wLive  , event_count);
  const auto bitsRaw    = extractTriggerBits(wRaw   , event_count);

  /* 4. per‑trigger loop -------------------------------------------- */
  for (const auto& [bitIdx, key] : triggerNameMap)
  {
    ++m_trigStat[key].tested;

    const bool firedScaled = checkTriggerCondition(bitsScaled, bitIdx);
    const bool firedLive   = checkTriggerCondition(bitsLive  , bitIdx);
    const bool firedRaw    = checkTriggerCondition(bitsRaw   , bitIdx);

    /* 4a. scalar counters ------------------------------------------ */
    auto safeFill = [&](const std::string& hname)
    {
      auto& H = qaHistogramsByTrigger[key];
      if (auto it = H.find(hname); it != H.end())
          static_cast<TH1I*>(it->second)->Fill(1);
    };
    if (firedRaw   ) safeFill("cnt_" + key + "_raw");
    if (firedLive  ) safeFill("cnt_" + key + "_live");
    if (firedScaled) safeFill("cnt_" + key + "_scaled");

    /* 4b. MB × Trigger correlation map ----------------------------- */
    if (h_MBTrigCorr)
    {
      const int cat = isMB ? (firedScaled ? 4 : 3)
                           : (firedScaled ? 2 : 1);
      h_MBTrigCorr->Fill(m_trigBin[key], cat);
    }

    /* 4c. accept trigger, do vertex‑cut QA ------------------------- */
    if (!firedScaled) continue;

    activeTrig.push_back(key);
    ++m_trigStat[key].fired;

    if (int vCut = extractVtxCut(key); vCut > 0)
    {
      const double vzAbs = std::fabs(m_vz);
      const int bin = (vzAbs < vCut) ? 1 : ((vzAbs == vCut) ? 2 : 3);

      auto& H = qaHistogramsByTrigger[key];
      if (auto it = H.find("h_vtxRelToCut_" + key); it != H.end())
        static_cast<TH1I*>(it->second)->Fill(bin);
    }
  } /* trigger loop */

  /* 4d. human‑readable summary ------------------------------------ */
  printTriggerSummary(activeTrig, wRaw, wLive, wScaled);

  /* 5. global vertex‑z veto --------------------------------------- */
  bool pass = isMB && !activeTrig.empty();
  if (pass && m_useVzCut && std::fabs(m_vz) >= m_vzCut) pass = false;

  LOG(2, CLR_BLUE, "  → firstEventCuts(): " << (pass ? "PASS" : "FAIL"));
  return pass;
}



// ======================================================================
//  process_event – one‑event driver, with streamlined logging
// ======================================================================
int emcal_sepdCorrelator::process_event(PHCompositeNode* topNode)
{
  /* ------------------------------------------------------------------ */
  /* 0.  Banner & running counter                                       */
  /* ------------------------------------------------------------------ */
  ++event_count;
  PROGRESS("=================================   event "
           << std::setw(4) << event_count
           << "   =====================================  (Verb="
           << Verbosity() << ')');

  /* ------------------------------------------------------------------ */
  /* 1.  Mandatory node check                                           */
  /* ------------------------------------------------------------------ */
  LOG(4, CLR_BLUE, "  [process_event] – node sanity");
  if (!fetchNodes(topNode))
  {
    LOG(4, CLR_YELLOW,
        "    mandatory node missing OR not Au+Au minimum‑bias →  ABORTEVENT");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /* ------------------------------------------------------------------ */
  /* 2.  Build the SEPD channel map on the fly (first event only)       */
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
    buildSepdChannelMap();     // sets m_sepdMapReady = true
  }

  /* ------------------------------------------------------------------ */
  /* 3.  Combined MB‑and‑Trigger gate (fills h_MBTrigCorr)              */
  /* ------------------------------------------------------------------ */
  std::vector<std::string> activeTrig;
  if (!firstEventCuts(topNode, activeTrig))
    {
        ++m_evtNoTrig;                   // keep previous statistics
        LOG(4, CLR_YELLOW,
            "    event rejected by MB/Trigger gate – skip");
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    /* 4.  Vertex‑z QA  -------------------------------------- */
    for (const auto& t : activeTrig) {
        static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_vertexZ"])
            ->Fill(m_vz);
    }


    /* ------------------------------------------------------------------ */
    /* 5.  Centrality lookup & diagnostics                                */
    /*      (must precede detector‑level QA so centrality clones fill)    */
    /* ------------------------------------------------------------------ */
    CentralityInfo* central =
          findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

    if (!central)
    {
        LOG(4, CLR_YELLOW,
            "    CentralityInfo node missing – skip");
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    const float centile =
        central->get_centrality_bin(CentralityInfo::PROP::mbd_NS);

    if (!std::isfinite(centile) || centile < 0.f)
    {
        LOG(4, CLR_YELLOW,
          "    mbd_NS centile invalid – treating as minimum‑bias (0 – 100 %)");
        m_centBin = -1;                          // minimum‑bias
    }
    else
    {
        m_centBin = static_cast<int>(centile);
        LOG(5, CLR_GREEN, "    centrality bin = " << m_centBin << '%');
    }

    /* centrality histogram (filled once the value is validated) -------- */
    if (centile >= 0.f && centile <= 100.f) {
        for (const auto& t : activeTrig) {
            static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_centrality"])
                ->Fill(centile);
        }
    }

    /* guard: vz must be reasonable for centrality calibration ---------- */
    if (!std::isfinite(m_vz) || std::abs(m_vz) > 60.0)
    {
        LOG(4, CLR_YELLOW,
            "    Vertex‑z (" << m_vz
            << " cm) outside calibration bounds – skip event");
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    /* ------------------------------------------------------------------ */
    /* 6.  Detector‑level QA                                              */
    /* ------------------------------------------------------------------ */
    LOG(5, CLR_BLUE, "    running detector‑level QA");
    doCaloQA(topNode, activeTrig);
    doSepdQA(activeTrig);      // SEPD charge maps & ψ_n
    doMbdQA(activeTrig);
    doPi0QA(activeTrig);
    fillCorrelations(activeTrig);
    fillFlowHists(activeTrig);

    if (doJetQA(topNode, activeTrig) == Fun4AllReturnCodes::ABORTRUN)
    {
      LOG(2, CLR_YELLOW, "    doJetQA requested ABORTRUN");
      return Fun4AllReturnCodes::ABORTRUN;
    }

    LOG(4, CLR_GREEN, "  [process_event] – completed OK");
    return Fun4AllReturnCodes::EVENT_OK;
}


//==========================================================================
//  fetchNodes – check presence of all required nodes & cache pointers
//==========================================================================
bool emcal_sepdCorrelator::fetchNodes(PHCompositeNode* top)
{
    
  /* ------------------------------------------------------------------ */
  /* 0.  Minimum‑bias information         */
  /* ------------------------------------------------------------------ */
  m_isMinBias = false;                                // reset per event

  if (auto* mbInfo = findNode::getClass<MinimumBiasInfo>(top,
                                                           "MinimumBiasInfo"))
        m_isMinBias = mbInfo->isAuAuMinimumBias();
  else
        LOG(1, CLR_YELLOW,
            "  – MinimumBiasInfo node missing (treating as !MB)");

  /* NOTE: Actual rejection is done in firstEventCuts(), so we only
   *       cache m_isMinBias here and keep processing. */

    
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

TH2F*
emcal_sepdCorrelator::makeEpdHitmap(const std::string& name,
                                    EpdGeom* /*geom*/, int /*arm*/)
{
    static const Double_t rEdge[17] =   // cm – inner radius of each ring
      { 0.15, 0.35, 0.55, 0.75, 0.95, 1.15, 1.35, 1.55,
        1.75, 1.95, 2.15, 2.35, 2.55, 2.75, 2.95, 3.15, 3.55 };

    auto* h = new TH2F(name.c_str(),
                       ";#varphi  [deg];r  [cm]",
                       24, 0., 360.,             // 24 × 15°
                       16, rEdge);               // *** non‑linear *** radii
    // Allow the uniform φ‑axis to grow if ever needed, but keep the
    // variable‑bin r‑axis fixed – this prevents the ROOT ExtendAxis
    // warning and preserves overflow statistics.
    h->SetCanExtend(TH1::kXaxis);
    return h;
}

// ==========================================================================
//  doSepdQA – charge / hit‑maps (global & centrality) + Ψ1,Ψ2,Ψ3 QA
//             now with correct polar filling (tile‑0 = 30°)
// ==========================================================================
void emcal_sepdCorrelator::doSepdQA(const std::vector<std::string>& trig)
{
  LOG(2, CLR_BLUE, "[doSepdQA] ──────────────────────────────────────────────");

  /* ------------------------------------------------------------------ */
  /* 0. Mandatory nodes present?                                        */
  /* ------------------------------------------------------------------ */
  if (!m_sepd || !m_epdgeom)
  {
    LOG(1, CLR_YELLOW,
        "[doSepdQA] mandatory SEPD nodes missing "
        "(m_sepd=" << m_sepd << ", m_epdgeom=" << m_epdgeom << ") – SKIP");
    return;
  }

  /* ------------------------------------------------------------------ */
  /* Helper lambdas                                                     */
  /* ------------------------------------------------------------------ */
  auto warnOnce = [this](const std::string& key)
  {
    static std::unordered_set<std::string> issued;
    if (issued.insert(key).second)
      LOG(1, CLR_YELLOW, "      [WARN] histogram key \"" << key
                                  << "\" is missing – first occurrence");
  };

  auto safeFillH1 = [](TObject* h, double x)
    { return h ? (static_cast<TH1*>(h)->Fill(x), true) : false; };

  /*  Works for TH2* *or* TProfile*.  Falls back silently if unknown.   */
  auto safeFillH2 = [](TObject* h, double x, double y, double w = 1.)
    {
      if (!h) return false;
      if (auto* p = dynamic_cast<TProfile*>(h)) { p->Fill(x, y, w); return true; }
      if (auto* h2 = dynamic_cast<TH2*>(h))     { h2->Fill(x, y, w); return true; }
      return false;
    };
    

  auto fillPolar = [](TH2* h, double phi, double r, double w,
                        unsigned epdKey)
    {
        if (!h) return;

        const double twoPi   = 2.*TMath::Pi();
        const double rad2deg = 180.0 / TMath::Pi();      // ROOT “POL” expects degrees
        auto wrapDeg = [twoPi, rad2deg](double a)
        {
            a = std::fmod(a, twoPi);                     // 0 ≤ φ < 2π  (rad)
            if (a < 0) a += twoPi;
            return a * rad2deg;                          // convert to degrees
        };

        const int ring = TowerInfoDefs::get_epd_rbin(epdKey);

        if (ring == 0)                                   // tile‑0 spans 30°
        {
            const double dphi = TMath::Pi() / 12.;       // 15° in rad
            h->Fill(wrapDeg(phi)       , r, w);
            h->Fill(wrapDeg(phi + dphi), r, w);
        }
        else
        {
            h->Fill(wrapDeg(phi), r, w);
        }
    };

  /* ------------------------------------------------------------------ */

  /* ------------------------------------------------------------------ */
  /* 1. Per‑event initialisation                                        */
  /* ------------------------------------------------------------------ */
  m_sepdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  float qxS[4] = {0}, qyS[4] = {0};
  float qxN[4] = {0}, qyN[4] = {0};

  /* ------------------------------------------------------------------ */
  /* 2. Determine centrality slice                                      */
  /* ------------------------------------------------------------------ */
  int         lo , hi ;
  std::string sliceTag;
  const bool  hasSlice = getCentralitySlice(lo,hi,sliceTag);

  /* ------------------------------------------------------------------ */
  /* 3. Channel loop  +  radial QA counters                              */
  /* ------------------------------------------------------------------ */
  int radialCnt[2][16] = {};          // [arm][ring]

  for (unsigned ch = 0, nChan = m_sepd->size(); ch < nChan; ++ch)
  {
      if (ch >= 744) continue;                         // hardware guard

      auto* ti = m_sepd->get_tower_at_channel(ch); if (!ti) continue;
      const double w = ti->get_energy();
      if (w <= 0.) continue;

      const unsigned key = m_epdKey[ch];
      if (key == std::numeric_limits<unsigned>::max()) continue;

      const int    arm  = TowerInfoDefs::get_epd_arm(key);   // 0=S,1=N
      const double phi  = m_epdgeom->get_phi(key);
      const int    ring = TowerInfoDefs::get_epd_rbin(key);  // 0 … 15
      const double r    = 0.25 + 0.20 * ring;               // centre of radial bin
      const double phiPlot = (phi < 0) ? phi + 2*M_PI : phi; // [0,2π)

      ++radialCnt[arm][ring];                               // QA counter
      
      /* ─── per‑ring QA fills ─────────────────────────────────────── */
      for (const std::string& t : trig)
      {
        auto& H = qaHistogramsByTrigger[t];

        /* global histograms ------------------------------------------------ */
        const std::string hOcc = (arm == 0)
            ? "h_SEPD_RingOcc_South_" + t
            : "h_SEPD_RingOcc_North_" + t;

        const std::string hQ   = (arm == 0)
            ? "h_SEPD_RingQ_South_"   + t
            : "h_SEPD_RingQ_North_"   + t;

        if (auto* ho = dynamic_cast<TH1*>(H[hOcc])) ho->Fill(ring);
        if (auto* hq = dynamic_cast<TH1*>(H[hQ ]))  hq->Fill(ring, w);

          if (hasSlice) {
            const std::string tag   = sliceTag + '_' + t;

            const std::string hOccC = (arm == 0)
                ? "h_SEPD_RingOcc_South_" + tag
                : "h_SEPD_RingOcc_North_" + tag;

            const std::string hQC   = (arm == 0)
                ? "h_SEPD_RingQ_South_"   + tag
                : "h_SEPD_RingQ_North_"   + tag;

            if (auto it = H.find(hOccC); it != H.end())
                static_cast<TH1*>(it->second)->Fill(ring);
            if (auto it = H.find(hQC); it != H.end())
                static_cast<TH1*>(it->second)->Fill(ring, w);
          }
      }


      /* ---- 3a. Hit‑maps (global + slice) ---------------------------- */
      const std::string baseKey = (arm == 0)
                                ? "h_sEPD_Hitmap_South_"
                                : "h_sEPD_Hitmap_North_";

      for (const std::string& t : trig)
      {
        auto& H = qaHistogramsByTrigger[t];

        // ---------- global hit‑map ------------------------------------
        if (auto* h = dynamic_cast<TH2*>(H[baseKey + t]))
          fillPolar(h, phiPlot, r, 1.0, key);
        else
          warnOnce(baseKey + t);

        // ---------- centrality slice ----------------------------------
        if (hasSlice) {
            const std::string keyC = baseKey.substr(0, baseKey.size() - 1)
                                   + sliceTag + '_' + t;
            auto itC = H.find(keyC);
            if (itC != H.end())
              if (auto* hC = dynamic_cast<TH2*>(itC->second))
                fillPolar(hC, phiPlot, r, 1.0, key);
         }
      }

      /* ---- 3b. Q‑vector sums (n = 1,2,3) ---------------------------- */
      const double c1 = std::cos(phi),       s1 = std::sin(phi);
      const double c2 = std::cos(2*phi),     s2 = std::sin(2*phi);
      const double c3 = std::cos(3*phi),     s3 = std::sin(3*phi);

      float* qx = (arm==0) ? qxS : qxN;
      float* qy = (arm==0) ? qyS : qyN;

      qx[1]+=w*c1; qy[1]+=w*s1;
      qx[2]+=w*c2; qy[2]+=w*s2;
      qx[3]+=w*c3; qy[3]+=w*s3;
      (arm==0 ? ++nFiredS : ++nFiredN);

      /* ---- 3c. ΣQ bookkeeping --------------------------------------- */
      m_sepdQ          += w;
      m_sepdQ_arm[arm] += w;
    } // end channel loop

    /* ------------------------------------------------------------------ */
    /* 3d.  RADIAL QA print‑out  (requires Verbosity() ≥ 1)               */
    /* ------------------------------------------------------------------ */
    if (Verbosity() > 0)
    {
      auto printArm = [&](int arm, const char* lbl)
      {
        std::ostringstream os;
        os << lbl << " : [";
        int tot = 0;
        for (int r = 0; r < 16; ++r)
        {
          os << std::setw(4) << radialCnt[arm][r] << (r==15?"]":" ");
          tot += radialCnt[arm][r];
        }
        os << "  (sum = " << tot << ")";
        LOG(1, CLR_BLUE, os.str());
      };
      printArm(0, "SEPD South");
      printArm(1, "SEPD North");
  }


  /* ------------------------------------------------------------------ */
  /* 4. ΣQ spectrum (per‑event)                                         */
  /* ------------------------------------------------------------------ */
  for (const std::string& t : trig)
    safeFillH1(qaHistogramsByTrigger[t]["h_towerQ_SEPD"], m_sepdQ);

  /* ------------------------------------------------------------------ */
  /* 5. Event‑plane reconstruction (Ψ1,Ψ2,Ψ3)                           */
  /*      … (unchanged) …                                               */
  /* ------------------------------------------------------------------ */
  bool usedMap=false;
  if (m_epmap && !m_epmap->empty())
  {
    auto* epdS = m_epmap->get(EventplaneinfoMap::sEPDS);
    auto* epdN = m_epmap->get(EventplaneinfoMap::sEPDN);
    if (epdS && epdN)
    {
      const auto q2S = epdS->get_qvector(2);
      const auto q2N = epdN->get_qvector(2);
      if ((q2S.first||q2S.second) && (q2N.first||q2N.second))
      {
        Eventplaneinfov1 h;
        m_psi2_S = h.GetPsi(q2S.first,q2S.second,2);
        m_psi2_N = h.GetPsi(q2N.first,q2N.second,2);
        usedMap  = true;
      }
    }
  }

  if (!usedMap)    // fallback: tower‑based
  {
    m_psi2_S = 0.5 * std::atan2(qyS[2],qxS[2]);
    m_psi2_N = 0.5 * std::atan2(qyN[2],qxN[2]);
  }

  m_psi1_S = (std::hypot(qxS[1],qyS[1])<1e-9) ? 0. : std::atan2(qyS[1],qxS[1]);
  m_psi1_N = (std::hypot(qxN[1],qyN[1])<1e-9) ? 0. : std::atan2(qyN[1],qxN[1]);

  m_psi3_S = (std::hypot(qxS[3],qyS[3])<1e-9) ? 0.
                                              : (1./3.)*std::atan2(qyS[3],qxS[3]);
  m_psi3_N = (std::hypot(qxN[3],qyN[3])<1e-9) ? 0.
                                              : (1./3.)*std::atan2(qyN[3],qxN[3]);

  /* ------------------------------------------------------------------ */
  /* 6. QA histograms (Ψn + sub‑event resolution proxies)               */
  /* ------------------------------------------------------------------ */
  const double cos1 = std::cos( m_psi1_N - m_psi1_S);
  const double cos2 = std::cos(2*(m_psi2_N - m_psi2_S));
  const double cos3 = std::cos(3*(m_psi3_N - m_psi3_S));

  for (const std::string& t : trig)
  {
      auto& H = qaHistogramsByTrigger[t];

      /* event‑plane shapes, South & North */
      safeFillH1(H["h_Psi1_sEPD_S"], m_psi1_S);
      safeFillH1(H["h_Psi1_sEPD_N"], m_psi1_N);
      safeFillH1(H["h_Psi2_sEPD_S"], m_psi2_S);
      safeFillH1(H["h_Psi2_sEPD_N"], m_psi2_N);
      safeFillH1(H["h_Psi3_sEPD_S"], m_psi3_S);
      safeFillH1(H["h_Psi3_sEPD_N"], m_psi3_N);

      /* sub‑event correlation versus total charge */
      safeFillH2(H["h_Psi1_res_vs_Qsum"], m_sepdQ, cos1);
      safeFillH2(H["h_Psi2_res_vs_Qsum"], m_sepdQ, cos2);
      safeFillH2(H["h_Psi3_res_vs_Qsum"], m_sepdQ, cos3);
  }


  /* ------------------------------------------------------------------ */
  /* 7. Summary                                                         */
  /* ------------------------------------------------------------------ */
  LOG(3, CLR_GREEN,
        "    ΣQ=" << m_sepdQ
        << "  Ψ1S=" << m_psi1_S << "  Ψ2S=" << m_psi2_S
        << "  Ψ3S=" << m_psi3_S << "  hits(S,N)=" << nFiredS << ',' << nFiredN
        << (usedMap ? "  (Ψ2 from EventplaneinfoMap)" : "  (Ψ2 from towers)"));
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



//==========================================================================
//  doMbdQA – integrated charge  +  hex‑hitmaps (global *and* centrality)
//==========================================================================
void emcal_sepdCorrelator::doMbdQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doMbdQA]");

  /* ------------------------------------------------------------------ */
  /* 0. Event‑level bookkeeping                                         */
  /* ------------------------------------------------------------------ */
  m_mbdQ      = 0.;
  std::size_t nFiredS = 0,            // PMT multiplicities
              nFiredN = 0;

  /* ------------------------------------------------------------------ */
  /* 1. Determine this event’s centrality slice once                    */
  /* ------------------------------------------------------------------ */
  int         lo , hi ;
  std::string sliceTag;
  const bool  hasSlice = getCentralitySlice(lo,hi,sliceTag);

  /* ------------------------------------------------------------------ */
  /* 2. Loop over all PMTs                                              */
  /* ------------------------------------------------------------------ */
  const unsigned nPmt = static_cast<unsigned>(m_mbdpmts->get_npmt());
  for (unsigned ip = 0; ip < nPmt; ++ip)
  {
    MbdPmtHit* p = m_mbdpmts->get_pmt(ip);
    if (!p) continue;

    const double q = p->get_q();        // calibrated charge
    if (q <= 0.) continue;              // skip empty PMTs

    const double cx = m_mbdgeom->get_x(ip);
    const double cy = m_mbdgeom->get_y(ip);

    /* decide which hit‑map this PMT belongs to --------------------- */
    const bool isSouth = (m_mbdgeom->get_arm(ip) == 0);
    const std::string baseKey = isSouth ?
                                "h_MBD_Hitmap_South_" :
                                "h_MBD_Hitmap_North_";

    for (const auto& t : trig)
    {
      /* ---- (a) global map --------------------------------------- */
      static_cast<TH2Poly*>(qaHistogramsByTrigger[t][baseKey + t])
          ->Fill(cx, cy, q);

      /* ---- (b) centrality‑tagged clone (only if cent info valid) - */
        if (hasSlice)
        {
          const std::string keyC =
              baseKey.substr(0, baseKey.size() - 1) + sliceTag + '_' + t;

          auto it = qaHistogramsByTrigger[t].find(keyC);
          if (it != qaHistogramsByTrigger[t].end())
            static_cast<TH2Poly*>(it->second)->Fill(cx, cy, q);
        }
    }

    /* ---- book‑keeping -------------------------------------------- */
    m_mbdQ += q;
    m_mbdQ_arm[isSouth ? 0 : 1] += q;
    isSouth ? ++nFiredS : ++nFiredN;
  } // PMT loop

  /* ------------------------------------------------------------------ */
  /* 3. Per‑event ΣQ spectrum                                          */
  /* ------------------------------------------------------------------ */
  for (const auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_charge_MBD"])
        ->Fill(m_mbdQ);

  /* ------------------------------------------------------------------ */
  /* 4. Verbose summary + downstream correlations                       */
  /* ------------------------------------------------------------------ */
  LOG(3, CLR_GREEN, "    MBD ΣQ = " << m_mbdQ
                     << "  (South PMTs: " << nFiredS
                     << ", North PMTs: " << nFiredN << ")");

  fillCentralityQA(trig);              // keeps existing behaviour
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

//----------------------------------------------------------------------
//  accumulateFlowContribution – single‑tower q‑vector bookkeeping
//----------------------------------------------------------------------
void emcal_sepdCorrelator::accumulateFlowContribution(const std::string& calorimeter,
                                                      unsigned           ieta,
                                                      double             et,
                                                      double             phi,
                                                      int                ptBin)
{
  if (et <= 0) return;                        // safety ‑ should never fire

  /*……… harmonic basis ……………………………………………………………………………………*/
  const double c1 = std::cos(phi),  s1 = std::sin(phi);
  const double c2 =  c1*c1 - s1*s1;          // cos 2φ
  const double s2 =  2.0*c1*s1;              // sin 2φ
  const double c3 =  c1*c2 - s1*s2;          // cos 3φ
  const double s3 =  s1*c2 + c1*s2;          // sin 3φ

  /*……… helper: safely add to a detector key …………………………………*/
  auto push = [&](const std::string& detKey)
  {
    if (!m_flowAcc.count(detKey) ||
        static_cast<std::size_t>(ptBin) >= m_flowAcc[detKey].size())
    {
      static std::unordered_set<std::string> warned;
      if (warned.insert(detKey).second)
        LOG(1, CLR_YELLOW, "      [WARN] detector \"" << detKey
                               << "\" missing or pT‑bin out of range");
      return;
    }

    auto& acc = m_flowAcc[detKey][ptBin];
    acc.sumW  += et;
    acc.qx[1] += et*c1;  acc.qy[1] += et*s1;
    acc.qx[2] += et*c2;  acc.qy[2] += et*s2;
    acc.qx[3] += et*c3;  acc.qy[3] += et*s3;
  };

  /*……… route to the correct detector buckets ……………………………*/
  if (calorimeter == "CEMC")
  {
        bool south = isSouthCEMC(ieta);
        push(south ? "CEMC_S" : "CEMC_N");
        push("CEMC_T");                       // combined S+N
        push(south ? "ALL_S"  : "ALL_N");
        push("ALL_T");                        // combined S+N
  }
  else if (calorimeter == "IHCAL")
  {
      bool south = isSouthHCal(ieta);
      push(south ? "IHCAL_S" : "IHCAL_N");
      push("IHCAL_T");
      push(south ? "HCAL_S"  : "HCAL_N");
      push("HCAL_T");
      push(south ? "ALL_S"   : "ALL_N");
      push("ALL_T");
  }
  else if (calorimeter == "OHCAL")
  {
      bool south = isSouthHCal(ieta);
      push(south ? "OHCAL_S" : "OHCAL_N");
      push("OHCAL_T");
      push(south ? "HCAL_S"  : "HCAL_N");
      push("HCAL_T");
      push(south ? "ALL_S"   : "ALL_N");
      push("ALL_T");
  }
}


//==========================================================================
//  doCaloQA – tower / cluster QA  +  per‑event vₙ (n = 1,2,3) accumulators
//==========================================================================
void emcal_sepdCorrelator::doCaloQA(
        PHCompositeNode*                   topNode,
        const std::vector<std::string>&    trig)
{
  LOG(3, CLR_BLUE, "[doCaloQA]  ─ processing calorimeter towers");
  /* ------------------------------------------------------------------ */
  /* Centrality slice helper                                            */
  /* ------------------------------------------------------------------ */
  int         centLo , centHi ;
  std::string centTag;
  const bool  hasCentSlice = getCentralitySlice(centLo,centHi,centTag);

  //------------------------------------------------------------------------
  // (0)  Reset per‑event flow accumulators
  //------------------------------------------------------------------------
  for (auto& [det, vec] : m_flowAcc)
    for (auto& acc : vec) acc.reset();

  //------------------------------------------------------------------------
  // (1)  Tower loop – iterate over every configured calorimeter subsystem
  //------------------------------------------------------------------------
  std::vector<std::pair<double,double>> cemcPos;          // (η,φ) cache

  for (auto& ck : m_calo)
  {
    /* handy aliases --------------------------------------------------- */
    const std::string&  label = ck.first;                 // "CEMC", "IHCAL", …
    auto*               twCon = ck.second.tw;             // TowerInfoContainer
    double&             sumE  = ck.second.sumE;           // ΣE in this det
    std::size_t         nHit  = 0;                        // fired towers

    /* tower loop ------------------------------------------------------ */
    for (unsigned ch = 0; ch < twCon->size(); ++ch)
    {
        auto* tw = twCon->get_tower_at_channel(ch);
        if (!tw)                        continue;
        
        const double e = tw->get_energy();
        if (e < m_towMinE)              continue;           // noise suppression
        
        /* 1‑A) simple book‑keeping + spectrum fill ---------------------- */
        sumE += e; ++nHit;
        for (const auto& t : trig)
            static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerE_" + label])
            ->Fill(e);
        
        /* 1‑B) geometry look‑up ---------------------------------------- */
        const unsigned key  = twCon->encode_key(ch);
        const unsigned iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
        const unsigned ieta = TowerInfoDefs::getCaloTowerEtaBin(key);
        
        RawTowerGeom* tg   = ck.second.g->get_tower_geometry(key);
        const double  eta  = tg ? tg->get_eta() : 0.0;
        const double  et   = e  / std::cosh(eta);           // transverse energy
        
        /* 1‑C) ΣE bookkeeping per arm ---------------------------------- */
        if (label == "CEMC")
        {
            isSouthCEMC(ieta) ? m_cemcEt_arm[0] += et
            : m_cemcEt_arm[1] += et;
        }
        else if (label == "IHCAL")
        {
            isSouthHCal(ieta) ? m_ihcalEt_arm[0] += et
            : m_ihcalEt_arm[1] += et;
        }
        else if (label == "OHCAL")
        {
            isSouthHCal(ieta) ? m_ohcalEt_arm[0] += et
            : m_ohcalEt_arm[1] += et;
        }
        
        /* 1‑D) η–φ hit‑map fills (global + centrality) ------------------ */
        std::string hMapPrefix;
        if      (label == "CEMC")  hMapPrefix = "h_EMC_EtaPhiMap_";
        else if (label == "IHCAL") hMapPrefix = "h_IHCAL_EtaPhiMap_";
        else if (label == "OHCAL") hMapPrefix = "h_OHCAL_EtaPhiMap_";
        
        if (!hMapPrefix.empty())
        {
            for (const auto& t : trig)
            {
                /* global map */
                static_cast<TH2F*>(qaHistogramsByTrigger[t][hMapPrefix + t])
                ->Fill(iphi, ieta, e);
                
                if (hasCentSlice)
                {
                    std::ostringstream k;
                    k << hMapPrefix.substr(0, hMapPrefix.size()-1)   // drop trailing '_'
                      << centTag << '_' << t;

                    auto it = qaHistogramsByTrigger[t].find(k.str());
                    if (it != qaHistogramsByTrigger[t].end())
                        static_cast<TH2F*>(it->second)->Fill(iphi, ieta, e);
                }
            }
        }
        
        /* ----------------------------------------------------------------
         * (A)  vₙ accumulators   (harmonics n = 1,2,3)
         *      – one entry per pT‑bin  ×  detector
         * ----------------------------------------------------------------*/
        int ptBin = -1;
        for (std::size_t b = 0; b < m_ptBins.size(); ++b)
            if (et >= m_ptBins[b].first && et < m_ptBins[b].second)
            { ptBin = static_cast<int>(b); break; }
        
        if (et < 0)
        {
            LOG(20, CLR_YELLOW, "      [WARN] negative et=" << et
                << " – tower skipped");
            continue;
        }
        
        if (ptBin < 0)          // clamp under/overflow to nearest bin edge
        {
            ptBin = (et < m_ptBins.front().first) ? 0
            : static_cast<int>(m_ptBins.size() - 1);
            
            static bool warned = false;
            if (!warned)
            {
                LOG(20, CLR_YELLOW, "      [INFO] et=" << et
                    << " outside configured pT range – clamped "
                    << "to bin " << ptBin);
                warned = true;
            }
        }
        
        double phi = tg ? tg->get_phi() : 0.0;        // (-π,π]
        if (phi < 0) phi += 2.*M_PI;                  //  [0,2π)
        
        /* cache CEMC tower positions for later Δη/Δφ comparisons -------- */
        if (label == "CEMC")
            cemcPos.emplace_back(eta, phi);
        
        /* Δη/Δφ resolution plots: nearest CEMC tower to an IHCAL tower -- */
        if (label == "IHCAL" && !cemcPos.empty())
        {
            double bestDEta = 999.0, bestDPhi = 999.0;
            for (const auto& ep : cemcPos)
            {
                const double dEta = eta - ep.first;
                const double dPhi = TVector2::Phi_mpi_pi(phi - ep.second);
                if (std::hypot(dEta, dPhi) < std::hypot(bestDEta, bestDPhi))
                {
                    bestDEta = dEta;  bestDPhi = dPhi;
                }
            }
            
            for (const auto& t : trig)
            {
                /* global */
                static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_dEta_CEMC_IHCAL"])
                ->Fill(bestDEta);
                static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_dPhi_CEMC_IHCAL"])
                ->Fill(bestDPhi);
                
                /* centrality‑tagged */
                if (hasCentSlice)
                {
                    const std::string kEta = "h_dEta_CEMC_IHCAL" + centTag + '_' + t;
                    const std::string kPhi = "h_dPhi_CEMC_IHCAL" + centTag + '_' + t;
                    
                    if (auto it = qaHistogramsByTrigger[t].find(kEta);
                        it != qaHistogramsByTrigger[t].end())
                        static_cast<TH1F*>(it->second)->Fill(bestDEta);
                    
                    if (auto it = qaHistogramsByTrigger[t].find(kPhi);
                        it != qaHistogramsByTrigger[t].end())
                        static_cast<TH1F*>(it->second)->Fill(bestDPhi);
                }
            }
        }
        
        accumulateFlowContribution(label, ieta, et, phi, ptBin);
    }

    LOG(3, CLR_GREEN, "    " << label << " : "
                             << nHit << " towers  |  ΣE = " << sumE);
  }   // end detector loop

  //------------------------------------------------------------------------
  // (2)  EMC cluster QA – simple E spectrum + per‑event max cluster‑E
  //------------------------------------------------------------------------
  if (!m_clus) return;

  float maxClusE = 0.f;
  const auto clRange = m_clus->getClusters();

  LOG(3, CLR_GREEN, "    EMC clusters : "
                    << std::distance(clRange.first, clRange.second));

  for (auto it = clRange.first; it != clRange.second; ++it)
  {
    const float eClus = static_cast<float>(it->second->get_energy());
    for (const auto& t : trig)
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_clusterE_EMC"])
          ->Fill(eClus);

    if (eClus > maxClusE) maxClusE = eClus;
  }

  /* ------------------------------------------------------------------ */
  /* (2‑A)  fill *max cluster‑E* histograms (+ optional verbose)        */
  /* ------------------------------------------------------------------ */
  if (maxClusE > 0.f)
  {
    /* determine centrality slice (once) ------------------------------ */
    int         lo , hi ;
    std::string tag;
    const bool  hasSlice = getCentralitySlice(lo,hi,tag);

    /* (A‑1) global & centrality‑tagged histograms -------------------- */
    for (const auto& t : trig)
    {
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_maxClusterE_EMC"])
          ->Fill(maxClusE);

        if (hasSlice)
        {
          const std::string key = "h_maxClusterE_EMC" + tag + '_' + t;
          auto&             H   = qaHistogramsByTrigger[t];
          if (auto it = H.find(key); it != H.end())
            static_cast<TH1F*>(it->second)->Fill(maxClusE);
        }

      /* extra diagnostics if Verbosity()>5 --------------------------- */
      if (Verbosity() > 5)
        LOG(6, CLR_MAGENTA, "      [maxClusE] trg=" << t
                            << "  centSlice=" << lo << "–" << hi
                            << "  Emax=" << maxClusE);
    }

    /* (A‑2)  do‑not‑scale reference histograms for turn‑on curves ---- */
    bool rawMbd150 = false;
    bool photonScaled[4] = {false,false,false,false};

    auto* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1) gl1 = findNode::getClass<Gl1Packet>(topNode, "14001");
    if (gl1)
    {
      const uint64_t wRaw    = gl1->lValue(0, "TriggerVector");
      const uint64_t wScaled = gl1->lValue(0, "ScaledVector");

      rawMbd150       = (wRaw    >> 14) & 0x1ULL;  // MBD_NS_geq_2_vtx_lt_150
      photonScaled[0] = (wScaled >> 20) & 0x1ULL;  // photon_6
      photonScaled[1] = (wScaled >> 21) & 0x1ULL;  // photon_8
      photonScaled[2] = (wScaled >> 22) & 0x1ULL;  // photon_10
      photonScaled[3] = (wScaled >> 23) & 0x1ULL;  // photon_12
    }

    if (rawMbd150)
    {
      /* (i) baseline MBD reference ----------------------------------- */
      {
        const std::string trg = "MBD_NS_geq_2_vtx_lt_150";
        auto& H = qaHistogramsByTrigger[trg];
        auto  it = H.find("h_maxClusterEnergy_doNotScale_" + trg);
        if (it != H.end())
          static_cast<TH1F*>(it->second)->Fill(maxClusE);
      }

      /* (ii) rare photon triggers ------------------------------------ */
      static const char* keyList[4] =
      {
        "photon_6_plus_MBD_NS_geq_2_vtx_lt_150",
        "photon_8_plus_MBD_NS_geq_2_vtx_lt_150",
        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
        "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
      };

      for (int i = 0; i < 4; ++i)
      {
        if (!photonScaled[i]) continue;              // live‑bit required
        const std::string trg = keyList[i];
        auto& H = qaHistogramsByTrigger[trg];
        auto  it = H.find("h_maxClusterEnergy_doNotScale_" + trg);
        if (it != H.end())
          static_cast<TH1F*>(it->second)->Fill(maxClusE);

        if (Verbosity() > 5)
          LOG(6, CLR_MAGENTA, "      [maxClusE] filled doNotScale for " << trg
                              << "  Emax=" << maxClusE);
      }
    }
  } // end if(maxClusE>0)
}



/* ----------------------------------------------------------------------
 * doPi0QA – γγ invariant‑mass spectra (global + centrality‑tagged)
 * -------------------------------------------------------------------- */
void emcal_sepdCorrelator::doPi0QA(const std::vector<std::string>& trig)
{
  /* 0)  Entrance message ------------------------------------------------ */
  LOG(3, CLR_BLUE,
      "[doPi0QA] entered  – "
      << (m_clus ? m_clus->size() : 0)
      << " EMC clusters available in this event");

  /* early exits --------------------------------------------------------- */
  if (!m_clus || m_clus->size() < 2)
  { LOG(4, CLR_YELLOW, "    [doPi0QA] < 2 clusters – nothing to do"); return; }

  /* build filtered cluster cache --------------------------------------- */
  const float EminMin = *std::min_element(m_minClusE.begin(), m_minClusE.end());
  const float chi2Max = *std::max_element(m_chi2Cuts.begin(), m_chi2Cuts.end());
  const float asymMax = *std::max_element(m_asymCuts.begin(), m_asymCuts.end());

  struct Clu { TLorentzVector v; float E, pt, chi; };
  std::vector<Clu> cl;  cl.reserve(m_clus->size());

  std::size_t nRejectedE = 0, nRejectedChi = 0;

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

  /* determine centrality slice ----------------------------------------- */
  int         lo , hi ;
  std::string centTag;
  const bool  hasSlice = getCentralitySlice(lo,hi,centTag);
  if (!hasSlice) { LOG(4, CLR_YELLOW, "    [doPi0QA] event outside slices – skip"); return; }


  LOG(4, CLR_CYAN,
      "    [doPi0QA] centrality = " << m_centBin
      << "%  → slice " << lo << "–" << hi << '%');

  /* loop over cluster pairs -------------------------------------------- */
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
      if (asym > asymMax) { ++pairCnt; continue; }      // global veto

      /* ---------- unconditional 2‑D correlation fills (leading‑cluster vars) ---- */
      {
          const float mInvPair = (c1.v + c2.v).M();

          /* --- pick the leading‑E cluster and its χ² --------------------------- */
          const bool  leadIsC1 = (c1.E >= c2.E);
          const float eLead    = leadIsC1 ? c1.E   : c2.E;   // leading‑cluster energy
          const float chiLead  = leadIsC1 ? c1.chi : c2.chi; // χ² of that same cluster

          /* helper lambda – unchanged ------------------------------------------ */
          auto fill2D = [&](const std::string& prefix,
                            float x, float y,
                            const std::string& trg)
          {
            auto& H = qaHistogramsByTrigger[trg];

            const std::string gKey = prefix + "_" + trg;
            if (auto it = H.find(gKey); it != H.end())
              static_cast<TH2F*>(it->second)->Fill(x, y);

            const std::string cKey = prefix + centTag + "_" + trg;
            if (auto it = H.find(cKey); it != H.end())
              static_cast<TH2F*>(it->second)->Fill(x, y);
          };

          for (const auto& t : trig)
          {
            fill2D("Minv_vs_Asym", asym,    mInvPair, t); // unchanged
            fill2D("Minv_vs_chi2", chiLead, mInvPair, t); // now leading‑cluster χ²
            fill2D("Minv_vs_Eavg", eLead,   mInvPair, t); // now leading‑cluster E
          }
      }

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

  /* summary & exit ---------------------------------------------------- */
  const std::size_t processedPairs = pairCnt.load();
  LOG(3, CLR_GREEN,
      "    [doPi0QA] completed – "
      << processedPairs << " / " << totalPairs
      << " pairs processed (" << std::fixed << std::setprecision(1)
      << (100.0 * processedPairs / totalPairs) << "%)"
      << " – centrality slice " << lo << "–" << hi << '%');
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
  int         lo , hi ;
  std::string tag;
  const bool  hasSlice = getCentralitySlice(lo,hi,tag);

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
      
      
    binsFilled[t] += safeFill(H["h_SEPD_S_vs_SEPD_N"], m_sepdQ_arm[0],
                                m_sepdQ_arm[1],
                                "h_SEPD_S_vs_SEPD_N", t);
    binsFilled[t] += safeFill(H["h_MBD_vs_CEMC" ], m_mbdQ, cemc ,
                              "h_MBD_vs_CEMC", t);
    binsFilled[t] += safeFill(H["h_MBD_vs_IHCAL"], m_mbdQ, ihcal,
                                "h_MBD_vs_IHCAL", t);
    binsFilled[t] += safeFill(H["h_MBD_vs_OHCAL"], m_mbdQ, ohcal,
                                "h_MBD_vs_OHCAL", t);

    /* --- NEW global EMCal ↔ HCal maps ---------------------------------- */
    binsFilled[t] += safeFill(H["h_IHCAL_vs_CEMC"], ihcal, cemc,
                                "h_IHCAL_vs_CEMC", t);
    binsFilled[t] += safeFill(H["h_OHCAL_vs_CEMC"], ohcal, cemc,
                                "h_OHCAL_vs_CEMC", t);

    auto tryCent = [&](const std::string& base,
                         double x, double y)
    {
        if (!hasSlice) return false;
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
        
    binsFilled[t] += tryCent("h_SEPD_S_vs_SEPD_N",
                                 m_sepdQ_arm[0], m_sepdQ_arm[1]);
    binsFilled[t] += tryCent("h_IHCAL_vs_CEMC", ihcal, cemc);
    binsFilled[t] += tryCent("h_OHCAL_vs_CEMC", ohcal, cemc);
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
int emcal_sepdCorrelator::doJetQA(PHCompositeNode*                topNode,
                                  const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doJetQA] – entering (centBin=" << m_centBin << ")");

  /* ------------------------------------------------------------------ */
  /* 0.  Work out the centrality slice that this event belongs to       */
  /* ------------------------------------------------------------------ */
  int         lo , hi ;
  std::string sliceTag;
  const bool  hasSlice = getCentralitySlice(lo,hi,sliceTag);

  const double psi[4] = {0., m_psi1_S, m_psi2_S, m_psi3_S};

  /* ------------------------------------------------------------------ */
  /* helpers (warn‑once + safe TH1/2/3 fill)                            */
  /* ------------------------------------------------------------------ */
  auto warnOnce = [this](const std::string& key)
  {
        static std::unordered_set<std::string> issued;
        if (issued.insert(key).second)
            LOG(1, CLR_YELLOW, "      [WARN] histogram key \"" << key
                                   << "\" is missing – first occurrence");
  };

  auto safeFillH1 = [&](TObject* o, double x)
  { if (auto* h = dynamic_cast<TH1*>(o)) h->Fill(x); else warnOnce(o ? o->GetName() : "null"); };

  /* ------------------------------------------------------------------ */
  /* 1.  Scan every configured jet radius                               */
  /* ------------------------------------------------------------------ */
  struct TwoJets { const Jet* lead = nullptr; const Jet* sub = nullptr; };

  std::unordered_map<std::string,float>   maxEt;
  std::unordered_map<std::string,TwoJets> bestPair;   // keyed by radKey

  for (const auto& [radKey, nodeName] : kJetRadii)
  {
    JetContainer* jets = findNode::getClass<JetContainer>(topNode, nodeName);
    if (!jets)
    {
      LOG(0, CLR_YELLOW, "      [FATAL] JetContainer \"" << nodeName
                         << "\" missing – ABORTRUN");
      return Fun4AllReturnCodes::ABORTRUN;
    }

    TwoJets pair;
    for (const Jet* j : *jets)
    {
      if (!j) continue;

      /* ---- leading & sub‑leading ---------------------------------- */
      if (!pair.lead || j->get_et() > pair.lead->get_et())
      { pair.sub = pair.lead; pair.lead = j; }
      else if (!pair.sub || j->get_et() > pair.sub->get_et())
      { pair.sub = j; }

      /* ---- v_n^{jet} (only if the event is in a slice) ------------ */
      if (hasSlice && std::abs(j->get_eta()) >= 3.0)
      {
        const double phi = j->get_phi();
        const double pt  = j->get_pt();

        for (int n : {1,2,3})
        {
          const double vn = std::cos(n * (phi - psi[n]));
          for (const auto& t : trig)
          {
            std::ostringstream key;
            key << "p_v" << n << "_JET_" << radKey << sliceTag << '_' << t;

            auto& H = qaHistogramsByTrigger[t];
            if (auto it = H.find(key.str()); it != H.end())
              static_cast<TProfile*>(it->second)->Fill(pt, vn);
          }
        }
      }
    }

    maxEt   [radKey] = pair.lead ? pair.lead->get_et() : 0.f;
    bestPair[radKey] = pair;

    LOG(4, CLR_GREEN, "      radius " << radKey
                        << "  jets=" << jets->size()
                        << "  maxE_T=" << maxEt[radKey]);
  }

    /* ------------------------------------------------------------------ */
    /* 2.  Histogram fills  +  verbose geometry diagnostics               */
    /* ------------------------------------------------------------------ */
    constexpr double etaMinDbg = -6.0, etaMaxDbg = 6.0;
    constexpr double phiMinDbg = -M_PI, phiMaxDbg =  M_PI;

    for (const auto& [radKey, etMax] : maxEt)
    {
      LOG(4, CLR_CYAN,
          "    [JetQA] radius " << radKey
          << "  maxE_T=" << etMax
          << "  (η range " << etaMinDbg << " … " << etaMaxDbg
          << ", φ range " << phiMinDbg << " … " << phiMaxDbg << " rad)");

      const std::string b1 = "h_maxJetEt_"          + radKey;
      const std::string b2 = "h_leadEt_vs_subEt_"   + radKey;
      const std::string b3 = "h_jetEt_area_nConst_" + radKey;
      const std::string b4 = "h_jetEtEtaPhi_"       + radKey;   // NEW

      const TwoJets& J = bestPair.at(radKey);

      /* fallback jet area if FastJet did not store one ----------------- */
      const double Rguess       = (radKey.size()>1 && radKey[0]=='r')
                                ? 0.1 * std::stod(radKey.substr(1)) : 0.4;
      const double areaFallback = M_PI * Rguess * Rguess;   // π R²

      for (const std::string& trg : trig)
      {
        auto& H = qaHistogramsByTrigger[trg];

        /* ---- 2.1  max‑E_T histograms -------------------------------- */
        safeFillH1(H[b1 + "_" + trg], etMax);                 // global
        if (hasSlice) safeFillH1(H[b1 + sliceTag + '_' + trg], etMax);

        /* ---- helpers for TH2 / TH3 ---------------------------------- */
        auto fill3 = [&](const std::string& key,double x,double y,double z)
        {
          TObject* o = H.count(key) ? H[key] : nullptr;
          if (auto* h = dynamic_cast<TH3F*>(o)) h->Fill(x,y,z); else warnOnce(key);
        };
        auto fill2 = [&](const std::string& key,double x,double y)
        {
          TObject* o = H.count(key) ? H[key] : nullptr;
          if (auto* h = dynamic_cast<TH2F*>(o)) h->Fill(x,y); else warnOnce(key);
        };

        /* ---- 2.2  jet‑shape histograms ------------------------------ */
        for (const Jet* j : {J.lead, J.sub})
        {
          if (!j) continue;

          const bool isLead = (j == J.lead);
          const double et   = j->get_et();
          const int    nC   = static_cast<int>(j->size_comp());
          const double area = areaFallback;                    // FastJet area not stored

          const double aOK  = std::clamp(area , 0.0, 0.80);
          const double nCOK = std::clamp<double>(nC, 0.0, 200);

          /* area–Nconst shape (existing) */
          fill3(b3 + "_" + trg, et, aOK, nCOK);               // global
          if (hasSlice) fill3(b3 + sliceTag + '_' + trg, et, aOK, nCOK);

          const double eta = j->get_eta();
          double phi = j->get_phi();                          // [-π,π] inclusive ⇒ make upper edge exclusive
          if (phi >= M_PI) phi -= 2.0 * M_PI;                 // shift +π into the first bin (−π)


          fill3(b4 + "_" + trg, eta, phi, et);                // global
          if (hasSlice) fill3(b4 + sliceTag + '_' + trg, eta, phi, et);

          /* -------- additional diagnostics (controlled by Verbosity) -------- */
          LOG(5, CLR_MAGENTA,
              "        " << (isLead ? "[lead] " : "[sub ] ")
              << "E_T=" << et << " GeV"
              << "  η="  << eta
              << "  φ="  << phi
              << "  area=" << aOK
              << "  nConst=" << nC);
          if (eta < etaMinDbg || eta > etaMaxDbg || phi < phiMinDbg || phi > phiMaxDbg)
            LOG(5, CLR_YELLOW,
                "        ↳ WARNING: (η,φ) outside booking range – histogram under/overflow");
        } // jet loop

        /* ---- 2.3  leading‑vs‑sub‑leading correlation ---------------- */
        if (J.lead && J.sub)
        {
          const double lEt = J.lead->get_et();
          const double sEt = J.sub ->get_et();

          fill2(b2 + "_" + trg, lEt, sEt);               // global
          if (hasSlice) fill2(b2 + sliceTag + '_' + trg, lEt, sEt);

          LOG(6, CLR_GREEN,
              "        Lead–Sub pair filled  (E_T^lead=" << lEt
              << ",  E_T^sub=" << sEt << ")");
        }
      } // trigger loop
    }   // radius loop

  LOG(3, CLR_GREEN, "  [doJetQA] – completed OK");
  return Fun4AllReturnCodes::EVENT_OK;
}


//––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
//  fillFlowHists – stores raw v₁, v₂, v₃  +  ⟨cos 2 ΔΨ⟩  vs centrality
//                  with lightweight, self‑throttling sanity checks
//––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
void
emcal_sepdCorrelator::fillFlowHists(const std::vector<std::string>& trig)
{
  /* (0) map this event to a <lo,hi> centrality slice ----------------- */
  int         lo , hi ;
  std::string sliceTag;
  const bool  hasSlice = getCentralitySlice(lo,hi,sliceTag);
  if (!hasSlice) return;     // skip events outside user slices

  /* -------- DEBUG: store the exact Q‑vectors that feed v_n --------- */
  for (const auto& [det, vec] : m_flowAcc)
      for (std::size_t ib = 0; ib < vec.size(); ++ib)
      {
          const FlowAcc& a = vec[ib];
          if (a.sumW <= 0.) continue;

          for (int n : {1,2,3})
          {
              const double qx = a.qx[n];
              const double qy = a.qy[n];

              for (const std::string& t : trig)
              {
                  std::ostringstream key;
                  key << "h_Q" << n << '_' << det << '_' << t;
                  auto& H = qaHistogramsByTrigger[t];
                  if (auto it = H.find(key.str()); it != H.end())
                      static_cast<TH2F*>(it->second)->Fill(qx, qy);
              }
          }
      }


  if (Verbosity() >= 3)
    LOG(3, CLR_CYAN, "[fillFlowHists] cent=" << m_centBin
                       << "  → slice " << lo << "–" << hi << '%');

  /* (1) harmonic basis for both SEPD sub‑event planes ---------------- */
  const double cPsi1S = std::cos(m_psi1_S),  sPsi1S = std::sin(m_psi1_S);
  const double cPsi2S = std::cos(2.*m_psi2_S), sPsi2S = std::sin(2.*m_psi2_S);
  const double cPsi3S = std::cos(3.*m_psi3_S), sPsi3S = std::sin(3.*m_psi3_S);

  const double cPsi1N = std::cos(m_psi1_N),  sPsi1N = std::sin(m_psi1_N);
  const double cPsi2N = std::cos(2.*m_psi2_N), sPsi2N = std::sin(2.*m_psi2_N);
  const double cPsi3N = std::cos(3.*m_psi3_N), sPsi3N = std::sin(3.*m_psi3_N);

  /* ------------------------------------------------------------------ *
   *  helper that prints **once** per anomaly type                      *
   * ------------------------------------------------------------------ */
  auto warnOnce = [&](const std::string& tag, const std::string& msg)
  {
    static std::unordered_set<std::string> told;
    if (told.insert(tag).second && Verbosity() >= 6)
      LOG(6, CLR_YELLOW, "[fillFlowHists] WARNING: " << msg);
  };

  /* (2) detector → pT‑bin loop -------------------------------------- */
  for (const auto& [det, vec] : m_flowAcc)
  {
    for (std::size_t ib = 0; ib < vec.size(); ++ib)
    {
      const FlowAcc& a = vec[ib];
      if (a.sumW <= 0.) continue;

      const double ptCtr = 0.5 * (m_ptBins[ib].first + m_ptBins[ib].second);

      /* -------------  AUTOCORRELATION‑SAFE PROJECTION  ----------------
       *  - South‑arm towers use the NORTH event‑plane and vice‑versa.
       *  - Combined “_T” buckets are skipped now and filled later from
       *    the weighted average of their *_S / *_N counterparts.
       * --------------------------------------------------------------- */
      if (det.size() >= 2 && det.rfind("_T") == det.size() - 2)
          continue;                                     // handle *_T after this loop

      const bool southArm = (det.rfind("_S") == det.size() - 2);

      const double c1 = southArm ? cPsi1N : cPsi1S;
      const double s1 = southArm ? sPsi1N : sPsi1S;
      const double c2 = southArm ? cPsi2N : cPsi2S;
      const double s2 = southArm ? sPsi2N : sPsi2S;
      const double c3 = southArm ? cPsi3N : cPsi3S;
      const double s3 = southArm ? sPsi3N : sPsi3S;

      const double v1 = (a.qx[1]*c1 + a.qy[1]*s1) / a.sumW;
      const double v2 = (a.qx[2]*c2 + a.qy[2]*s2) / a.sumW;
      const double v3 = (a.qx[3]*c3 + a.qy[3]*s3) / a.sumW;

      /* (2a)   sanity checks (printed only once per flavour) --------- */
      auto chk = [&](double v, const char* lab)
      {
          if (!std::isfinite(v))
          {
            std::ostringstream msg;
            msg << "non‑finite " << lab << " for det=" << det << ", bin=" << ib;
            warnOnce(det + lab + "_nan", msg.str());
          }
          else if (std::fabs(v) > 1.5)
          {
            std::ostringstream msg;
            msg << '|' << lab << "| = " << v
                << " > 1.5  (det=" << det << ", bin=" << ib << ")";
            warnOnce(det + lab + "_big", msg.str());
          }
      };
      chk(v1, "v1");  chk(v2, "v2");  chk(v3, "v3");

      /* (3) fill per‑trigger profiles (unit weight) ------------------ */
      for (const std::string& t : trig)
      {
        auto& H = qaHistogramsByTrigger[t];

        const std::string k1 = Form("p_v1_%s_%d_%d_%s", det.c_str(), lo, hi, t.c_str());
        const std::string k2 = Form("p_v2_%s_%d_%d_%s", det.c_str(), lo, hi, t.c_str());
        const std::string k3 = Form("p_v3_%s_%d_%d_%s", det.c_str(), lo, hi, t.c_str());

        if (auto* p = dynamic_cast<TProfile*>(H[k1])) p->Fill(ptCtr, v1, 1.0);
        if (auto* p = dynamic_cast<TProfile*>(H[k2])) p->Fill(ptCtr, v2, 1.0);
        if (auto* p = dynamic_cast<TProfile*>(H[k3])) p->Fill(ptCtr, v3, 1.0);
      }

        if (Verbosity() >= 7)
        {
            const double ptLo = m_ptBins[ib].first;   // lower edge of this pT bin
            const double ptHi = m_ptBins[ib].second;  // upper edge of this pT bin
            LOG(7, CLR_MAGENTA, "  det=" << det << "  ib=" << ib
                            << "  pT=(" << ptLo << "–" << ptHi << ")"
                            << "  v=( " << v1 << ", " << v2 << ", " << v3 << " )"
                            << "  Σw=" << a.sumW);
        }
     }
  }

  /* (4) store ⟨cos n ΔΨ⟩ (n = 1,2,3) for later Rₙ extraction ----------- */
  const double cos1 = std::cos( m_psi1_N - m_psi1_S );
  const double cos2 = std::cos( 2.*(m_psi2_N - m_psi2_S) );
  const double cos3 = std::cos( 3.*(m_psi3_N - m_psi3_S) );

  auto chkCos = [&](double c, int n)
  {
      if (!std::isfinite(c) || std::fabs(c) > 1.0)
        warnOnce(Form("cos%d_out_of_range", n),
                 Form("cos %dΔΨ = %.3f outside [-1,1]", n, c));
  };
  chkCos(cos1,1);  chkCos(cos2,2);  chkCos(cos3,3);

  for (const std::string& t : trig)
  {
      if (auto* p = dynamic_cast<TProfile*>(qaHistogramsByTrigger[t]["p_R1_vs_cent_" + t]))
        p->Fill(static_cast<double>(lo), cos1, 1.0);
      if (auto* p = dynamic_cast<TProfile*>(qaHistogramsByTrigger[t]["p_R2_vs_cent_" + t]))
        p->Fill(static_cast<double>(lo), cos2, 1.0);
      if (auto* p = dynamic_cast<TProfile*>(qaHistogramsByTrigger[t]["p_R3_vs_cent_" + t]))
        p->Fill(static_cast<double>(lo), cos3, 1.0);
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
  out->cd("triggerQA");
  if (h_MBTrigCorr && h_MBTrigCorr->GetEntries() > 0) h_MBTrigCorr->Write();

  //--------------------------------------------------------------------
  // 4.  Write footer & close the file
  //--------------------------------------------------------------------
  if (Verbosity() >= 1)
      std::cout << "\nOutput ROOT file →  " << out->GetName() << "\n\n";
    
  info(1, "writing TFile footer and closing ("+std::to_string(nHistWritten)
           +" / "+std::to_string(nHistExpected)+" objects written)");

  try      { out->Close(); }
  catch (const std::exception& e)
  { warn("Exception during TFile::Write/Close – "+std::string(e.what())); }

  delete out; out = nullptr;          // safe: we no longer dereference histos
  info(0, "Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}


void emcal_sepdCorrelator::Print(const std::string&) const
{ /* nothing to print beyond ROOT histograms */ }
