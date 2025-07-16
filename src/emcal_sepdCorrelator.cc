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

//––– ROOT & CLHEP ----------------------------------------------------------
#include <TProfile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH2Poly.h>
#include <CLHEP/Vector/ThreeVector.h>
//––– CDB access ------------------------------------------------------------
#include <cdbobjects/CDBTTree.h>        // <-- defines CDBTTree
#include <ffamodules/CDBInterface.h>     // <-- gives you CDBInterface::instance()


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
  do { std::cout << CLR_CYAN << MSG << CLR_RESET << std::endl; } while (false)

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

/* ======================================================================
 *  Init – one-time module setup
 *    • books QA histograms
 *    • (optional) prints the full DST node tree (Verbosity ≥ 2)
 * ====================================================================*/
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


//======================================================================
//  InitRun – geometry‑dependent booking (only once per run)
//======================================================================
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

  LOG(1, CLR_BLUE, "[InitRun] InitRun completed successfully");
  return Fun4AllReturnCodes::EVENT_OK;
}



//==========================================================================
//  bookShapeHitMaps – hex (MBD) & polar (sEPD) hit‑maps, one per trigger
//==========================================================================
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
  }
  out->cd();
}

//–––––––––––––––––––– helper book‑ers (tower/cluster, charge, correlations,
//                                  π0 spectra, EP/centrality) ––––––––––––
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

/* ----------------------------------------------------------------------
 * bookEnergyChargeCorrel  – detector–detector ΣE / ΣQ correlation maps
 *                           (global   +   centrality‑tagged clones)
 * -------------------------------------------------------------------- */
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


void emcal_sepdCorrelator::bookEventPlaneCentralityQA(const std::string& trig, HistMap& H)
{
  /* 1. #SigmaQ spectra */
  H["h_Qsum_MBD"]  = new TH1F(("h_Qsum_MBD_"  + trig).c_str(),
                              "MBD #SigmaQ;#SigmaQ_{MBD} [ADC]",   600, 0, 1200);
  H["h_Qsum_sEPD"] = new TH1F(("h_Qsum_sEPD_" + trig).c_str(),
                              "sEPD #SigmaQ;#SigmaQ_{sEPD} [ADC]", 600, 0, 1200);

  /* 2. detector‑to‑detector #SigmaQ map */
  H["h_Qsum_MBD_vs_sEPD"] = new TH2F(("h_Qsum_MBD_vs_sEPD_" + trig).c_str(),
                                     "#Sigma Q_{MBD} vs #Sigma Q_{sEPD};#Sigma Q_{MBD};#Sigma Q_{sEPD}",
                                     300, 0, 1200, 300, 0, 1200);

  /* 3. #Psi₂ distributions & resolution proxy */
  H["h_Psi2_sEPD"] =
      new TH1F(("h_Psi2_sEPD_" + trig).c_str(), "sEPD #Psi_{2};#Psi_{2} [rad]",
               120, -TMath::Pi(), TMath::Pi());

  H["h_Psi2_res_vs_Qsum"] =
      new TProfile(("h_Psi2_res_vs_Qsum_" + trig).c_str(),
                   "cos 2(#Psi_{N}-#Psi_{S}) vs #Sigma Q_{sEPD};#Sigma Q_{sEPD};#langle cos2Δ#Psi #rangle",
                   12, 0, 1200, "s");
}

// ----------------------------------------------------------------------
// Book max‑jet‑E_T spectra (global + centrality‑tagged clones)
// ----------------------------------------------------------------------
void
emcal_sepdCorrelator::bookJetQA(const std::string& trig, HistMap& H)       // <<< NEW
{
  const int nbEt = 200;            // 1 GeV per bin
  const double etMax = 200.;

  for (const auto& r : kJetRadii)
  {
    const std::string base = std::string("h_maxJetEt_") + r.first;

    /* global histogram */
    H[base + "_" + trig] =
        new TH1F((base + "_" + trig).c_str(),
                 ("max jet E_{T} ("+std::string(r.first)+");E_{T} [GeV]").c_str(),
                 nbEt, 0, etMax);

    /* centrality‑tagged clones */
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
    {
      int lo = m_centEdges[i];
      int hi = m_centEdges[i + 1];
      std::ostringstream n;
      n << base << '_' << lo << '_' << hi << '_' << trig;
      H[n.str()] =
          new TH1F(n.str().c_str(),
                   ("max jet E_{T} ("+std::string(r.first)+");E_{T} [GeV]").c_str(),
                   nbEt, 0, etMax);
    }
  }
}


//==========================================================================
//  createHistos_Data – scalar QA & QA maps
//==========================================================================
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
    bookEnergyChargeCorrel(trig, H);
    bookPi0MassSpectra(trig, H);
    bookEventPlaneCentralityQA(trig, H);
    bookJetQA(trig, H);
    H["h_vertexZ"] = new TH1F(("h_vertexZ_" + trig).c_str(),
                                "Primary vertex z;z_{vtx} [cm]",
                                240, -60., 60.);
      
    out->cd();
  }
}

// ----------------------------------------------------------------------
//  buildSepdChannelMap – executed once, after TOWERINFO_CALIB_SEPD exists
// ----------------------------------------------------------------------
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

//==========================================================================
//  process_event – orchestration with exhaustive verbosity
//==========================================================================
int emcal_sepdCorrelator::process_event(PHCompositeNode* topNode)
{
  /* ------------------------------------------------------------------ */
  /* 0. Banner & running counter                                        */
  /* ------------------------------------------------------------------ */
  ++event_count;
  PROGRESS("=================================   event " << std::setw(3) << event_count << "===================================== "
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
    
  /* ------------------------------------------------------------------ */
  /* 7. Remaining QA                                                    */
  /* ------------------------------------------------------------------ */
  doMbdQA (activeTrig);
  doPi0QA (activeTrig);
  fillCorrelations(activeTrig);

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
//  doCaloQA – tower and cluster spectra
//==========================================================================
void emcal_sepdCorrelator::doCaloQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doCaloQA] processing calorimeter towers");

  for (auto& ck : m_calo)
  {
    const std::string& lbl = ck.first;
    TowerInfoContainer* twC = ck.second.tw;
    double& sumE = ck.second.sumE;
    std::size_t nHit = 0;

    for (unsigned ch = 0; ch < twC->size(); ++ch)
    {
      auto* tw = twC->get_tower_at_channel(ch); if (!tw) continue;
      const double e = tw->get_energy();        if (e <= 0) continue;

      sumE += e; ++nHit;

      for (auto& t : trig)
        static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerE_" + lbl])->Fill(e);

      /* η‑φ maps & transverse‑energy arm sums */
      unsigned int key  = twC->encode_key(ch);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(key);

      RawTowerGeom* tg   = ck.second.g->get_tower_geometry(key);
      const double eta_t = tg ? tg->get_eta() : 0.;
      const double et    = e / std::cosh(eta_t);

      if (lbl == "CEMC")
      {
        if (eta_t < 0) m_cemcEt_arm[0] += et;
        else           m_cemcEt_arm[1] += et;
      }
      else if (lbl == "IHCAL")
      {
        if (eta_t < 0) m_ihcalEt_arm[0] += et;
        else           m_ihcalEt_arm[1] += et;
      }
      else if (lbl == "OHCAL")
      {
        if (eta_t < 0) m_ohcalEt_arm[0] += et;
        else           m_ohcalEt_arm[1] += et;
      }

      for (auto& t : trig)
      {
        if      (lbl == "CEMC")
          static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_EMC_EtaPhiMap_"  + t])->Fill(iphi, ieta, e);
        else if (lbl == "IHCAL")
          static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_IHCAL_EtaPhiMap_" + t])->Fill(iphi, ieta, e);
        else if (lbl == "OHCAL")
          static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_OHCAL_EtaPhiMap_" + t])->Fill(iphi, ieta, e);
      }
    } // tower loop

    LOG(3, CLR_GREEN, "    " << lbl << " : " << nHit << " fired, #SigmaE = " << sumE);
  }

  /* EMC clusters */
  if (!m_clus) return;
  RawClusterContainer::ConstRange cr = m_clus->getClusters();
  LOG(3, CLR_GREEN, "    EMC clusters : " << std::distance(cr.first, cr.second));

  for (auto it = cr.first; it != cr.second; ++it)
    for (auto& t : trig)
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_clusterE_EMC"])
          ->Fill(it->second->get_energy());
}

//==========================================================================
//  doSepdQA – sEPD charge, hit-map & event-plane QA   (✱ EXTENDED VERBOSITY ✱)
//==========================================================================
void emcal_sepdCorrelator::doSepdQA(const std::vector<std::string>& trig)
{
  /* ------------------------------------------------------------------ */
  /* 0. Banner + sanity on pointers                                     */
  /* ------------------------------------------------------------------ */
  LOG(2, CLR_BLUE, "[doSepdQA] ===============================================");
  LOG(2, CLR_BLUE, "[doSepdQA] m_sepd="   << (m_sepd   ? "OK":"nullptr")
                       << "  m_epdgeom="  << (m_epdgeom? "OK":"nullptr")
                       << "  m_epmap="    << (m_epmap  ? "OK":"nullptr"));

  if (!m_sepd || !m_epdgeom)
  {
    LOG(1, CLR_YELLOW, "[doSepdQA] SEPD container or geometry missing – skip");
    return;
  }

  /* ------------------------------------------------------------------ */
  /* 1. Initialise accumulators                                         */
  /* ------------------------------------------------------------------ */
  m_sepdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;
  double qxS = 0., qyS = 0.;
  double qxN = 0., qyN = 0.;

  /* ------------------------------------------------------------------ */
  /* 2. Tower loop                                                      */
  /* ------------------------------------------------------------------ */
  const auto nChan = m_sepd->size();
  LOG(2, CLR_CYAN, "[doSepdQA] iterating over " << nChan << " SEPD channels");

  for (unsigned ch = 0; ch < nChan; ++ch)
  {
    if (ch >= 744) continue;
    auto *ti = m_sepd->get_tower_at_channel(ch);
    if (!ti)
    {
      LOG(6, CLR_YELLOW, "  ch="<<ch<<" → nullptr – continue");
      continue;
    }

    const double w = ti->get_energy();
    if (w <= 0)
    {
      LOG(11, CLR_YELLOW, "  ch="<<ch<<" → E="<<w<<" ≤ 0 – continue");
      continue;
    }

    const unsigned key = m_epdKey[ch];
    if (key == std::numeric_limits<unsigned>::max())
    {
      LOG(11, CLR_YELLOW, "  ch="<<ch<<" → unmapped (999) – continue");
      continue;                                   // empty electronics slot
    }

    const int    arm = TowerInfoDefs::get_epd_arm(key); // 0=S,1=N
    const double r   = m_epdgeom->get_r  (key);
    const double phi = m_epdgeom->get_phi(key);

    LOG(11, CLR_GREEN, "  ch="<<ch
                 <<" arm="<<arm<<" r="<<r<<" phi="<<phi<<" w="<<w);

    /* ---- hit-map fill --------------------------------------------- */
    const std::string hpfx = (arm == 0)
                             ? "h_sEPD_Hitmap_South_" : "h_sEPD_Hitmap_North_";
    const double ph = (phi < 0) ? phi + 2.*M_PI : phi;   // 0…2π

    for (const auto& t : trig)
      static_cast<TH2F*>(qaHistogramsByTrigger[t][hpfx + t])
          ->Fill(ph, r, w);

    /* ---- Q-vector -------------------------------------------------- */
    const double c2 = std::cos(2 * phi);
    const double s2 = std::sin(2 * phi);
    if (arm == 0)
    {
      qxS += w * c2;  qyS += w * s2;  ++nFiredS;
      LOG(11, CLR_CYAN, "    -> South: qxS="<<qxS<<"  qyS="<<qyS
                           <<"  nS="<<nFiredS);
    }
    else
    {
      qxN += w * c2;  qyN += w * s2;  ++nFiredN;
      LOG(11, CLR_CYAN, "    -> North: qxN="<<qxN<<"  qyN="<<qyN
                           <<"  nN="<<nFiredN);
    }

    m_sepdQ          += w;
    m_sepdQ_arm[arm] += w;
    LOG(11, CLR_CYAN, "    -> running ΣQ = "<<m_sepdQ);
  } // channel loop

  LOG(3, CLR_GREEN, "[doSepdQA] tower loop finished:"
                    <<" ΣQ="<<m_sepdQ
                    <<"  SouthHits="<<nFiredS
                    <<"  NorthHits="<<nFiredN);

    /* ------------------------------------------------------------------ */
    /* 3. Scalar ΣQ histogram fill                                        */
    /* ------------------------------------------------------------------ */
    for (const auto& t : trig)
    {
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerQ_SEPD"])
          ->Fill(m_sepdQ);
    }

    /* ------------------------------------------------------------------ */
    /* 4. Event-plane angles from EventplaneinfoMap (⇒ Ψ2)                */
    /* ------------------------------------------------------------------ */

    bool usedMap = false;

    /* ---- 4.1  Map sanity --------------------------------------------- */
    if (!m_epmap)
    {
      LOG(5, CLR_YELLOW,
          "[doSepdQA] EventplaneinfoMap pointer is NULL – "
          "EventPlaneReco probably never ran (check module order)");
    }
    else if (m_epmap->empty())
    {
      LOG(5, CLR_YELLOW,
          "[doSepdQA] EventplaneinfoMap is EMPTY – "
          "all arms failed quality cuts or EventPlaneReco executed too early");
    }
    else
    {
      LOG(11, CLR_BLUE,
          "[doSepdQA] EventplaneinfoMap contains South & North objects?");

      /* ---- 4.2  South / North arm lookup ----------------------------- */
      auto* epdS = m_epmap->get(EventplaneinfoMap::sEPDS);
      auto* epdN = m_epmap->get(EventplaneinfoMap::sEPDN);

      LOG(11, CLR_BLUE,
          "              South = " << (epdS ? "present" : "MISSING")
          << "  |  North = "      << (epdN ? "present" : "MISSING"));

      if (!epdS || !epdN)
      {
        if (!epdS)
          LOG(5, CLR_YELLOW,
              "              South arm object absent – "
              "did it fail <min_tiles>/<qsum> cuts in EventPlaneReco?");
        if (!epdN)
          LOG(5, CLR_YELLOW,
              "              North arm object absent – "
              "did it fail <min_tiles>/<qsum> cuts in EventPlaneReco?");
      }
      else
      {
        /* ---- 4.3  Q-vector content check ----------------------------- */
        const auto q2S = epdS->get_qvector(2);
        const auto q2N = epdN->get_qvector(2);

        LOG(11, CLR_BLUE,
            "              Q2-vectors  S=(" << q2S.first << ", " << q2S.second
                                           << ")  N=(" << q2N.first << ", " << q2N.second << ')');

        const bool badS = (q2S.first == 0. && q2S.second == 0.);
        const bool badN = (q2N.first == 0. && q2N.second == 0.);

        if (badS || badN)
        {
          if (badS)
            LOG(5, CLR_YELLOW,
                "              South Q2 = (0,0) – quality cut inside EventPlaneReco");
          if (badN)
            LOG(5, CLR_YELLOW,
                "              North Q2 = (0,0) – quality cut inside EventPlaneReco");
        }
        else
        {
          /* ---- 4.4  All good → compute Ψ2 ----------------------------- */
          Eventplaneinfov1 helper;          // supplies GetPsi()
          m_psi2_S = helper.GetPsi(q2S.first, q2S.second, 2);
          m_psi2_N = helper.GetPsi(q2N.first, q2N.second, 2);
          usedMap  = true;
        }
      }
    }

    /* ---- 4.5  Fallback ------------------------------------------------ */
    if (!usedMap)
    {
      LOG(5, CLR_YELLOW,
          "[doSepdQA] EP-map unusable – falling back to tower-based Q-vectors");
      m_psi2_S = 0.5 * std::atan2(qyS, qxS);
      m_psi2_N = 0.5 * std::atan2(qyN, qxN);
    }

    LOG(4, CLR_GREEN,
        "[doSepdQA] Ψ₂(S) = " << m_psi2_S << "   Ψ₂(N) = " << m_psi2_N);

  const double cos2dPsi = std::cos(2 * (m_psi2_N - m_psi2_S));

  /* ------------------------------------------------------------------ */
  /* 5. EP histograms                                                   */
  /* ------------------------------------------------------------------ */
  for (const auto& t : trig)
  {
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Psi2_sEPD"])
        ->Fill(m_psi2_S);
    static_cast<TProfile*>(qaHistogramsByTrigger[t]["h_Psi2_res_vs_Qsum"])
        ->Fill(m_sepdQ, cos2dPsi);
  }

  /* ------------------------------------------------------------------ */
  /* 6. Footer message                                                  */
  /* ------------------------------------------------------------------ */
  LOG(3, CLR_GREEN, "    SEPD ΣQ="<<m_sepdQ
                     <<"  (South hits="<<nFiredS
                     <<", North hits="<<nFiredN<<")");
  LOG(4, CLR_BLUE,  "[doSepdQA] ===============================================");
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
//  doJetQA – fill jet QA histograms with exhaustive diagnostics
// ----------------------------------------------------------------------
int emcal_sepdCorrelator::doJetQA(PHCompositeNode*              topNode,
                                  const std::vector<std::string>& trig)
{
  LOG(4, CLR_BLUE, "  [doJetQA] – entering");

  /* ------------------------------------------------------------------ */
  /* 1.  Gather max-E_T for every configured jet radius                 */
  /* ------------------------------------------------------------------ */
  std::unordered_map<std::string,float> maxJetEt;   // key = "R02", "R04", …

  for (const auto& [radKey, nodeName] : kJetRadii)
  {
    LOG(5, CLR_BLUE, "    checking jet node '" << nodeName << "' (" << radKey << ')');

    JetContainer* jets = findNode::getClass<JetContainer>(topNode, nodeName);
    if (!jets)
    {
      LOG(4, CLR_YELLOW,
          "    missing JetContainer '" << nodeName << "' – ABORTRUN");
      return Fun4AllReturnCodes::ABORTRUN;
    }

    const float etMax = getMaxJetEt(jets);
    maxJetEt[radKey]  = etMax;

    LOG(5, CLR_GREEN, "      max E_T = " << etMax << " GeV");
  }

  /* ------------------------------------------------------------------ */
  /* 2.  Determine centrality slice tag                                 */
  /* ------------------------------------------------------------------ */
  int lo = 0, hi = 100;
  if (m_centBin >= 0)
    for (std::size_t i = 0; i + 1 < m_centEdges.size(); ++i)
      if (m_centBin >= m_centEdges[i] && m_centBin < m_centEdges[i + 1])
      { lo = m_centEdges[i]; hi = m_centEdges[i + 1]; break; }

  const std::string tag = '_' + std::to_string(lo) + '_' + std::to_string(hi);

  LOG(5, CLR_CYAN, "    centrality slice = " << lo << "–" << hi << " %  (tag '" << tag << "')");

  /* ------------------------------------------------------------------ */
  /* 3.  Fill histograms                                                */
  /* ------------------------------------------------------------------ */
  for (const auto& [radKey, etMax] : maxJetEt)
  {
    const std::string base = "h_maxJetEt_" + radKey;

    for (const auto& t : trig)
    {
      auto& H = qaHistogramsByTrigger[t];

      /* global histogram – must exist */
      auto* hGlobal = dynamic_cast<TH1F*>(H[base + "_" + t]);
      if (hGlobal) hGlobal->Fill(etMax);
      else LOG(4, CLR_YELLOW, "      [WARN] missing hist '" << base << '_' << t << '\'');

      /* centrality-tagged clone – may or may not exist                   */
      const std::string cKey = base + tag + "_" + t;
      auto it = H.find(cKey);
      if (it != H.end())
      {
        dynamic_cast<TH1F*>(it->second)->Fill(etMax);
        LOG(6, CLR_GREEN, "        • wrote " << radKey << " to '" << cKey << '\'');
      }
      else if (Verbosity() >= 6)
        LOG(6, CLR_BLUE,  "        • no cent-clone '" << cKey << '\'');
    }
  }

  LOG(4, CLR_GREEN, "  [doJetQA] – finished OK");
  return Fun4AllReturnCodes::EVENT_OK;
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
  m_psi2_S = 0.;
  m_psi2_N = 0.;
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
          std::cout << std::left << std::setw(30) << trig << " │ "
                    << std::setw(32) << key  << " │ "
                    << std::right<< std::setw(10)
                    << static_cast<Long64_t>(h->GetEntries()) << '\n';

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
  // 4.  Write footer & close the file (histograms already written)
  //--------------------------------------------------------------------
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
