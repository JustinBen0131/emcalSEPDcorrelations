//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator – implementation
//==========================================================================

#include "emcal_sepdCorrelator.h"

//––– Fun4All / PHOOL -------------------------------------------------------
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <TProfile.h>
#include <CLHEP/Vector/ThreeVector.h>

//––– Other ROOT headers ----------------------------------------------------
#include <TDirectory.h>
#include <TSystem.h>
#include <TMath.h>
#include <globalvertex/GlobalVertex.h>   // full definition of GlobalVertex
#include <calobase/TowerInfo.h>   // full definition of TowerInfo
#include <iostream>
#include <mbd/MbdPmtHit.h>
#include <calobase/RawCluster.h>    // full definition of RawClust
#include <calobase/RawClusterUtility.h>
#include <calobase/TowerInfoDefs.h>
#include <epd/EpdGeom.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

//–––––––– helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#define COUT_BLUE(MSG) \
    if (Verbosity() > 0) std::cout << "\033[1;34m" << MSG << "\033[0m\n"


#define CLR_BLUE   "\033[1;34m"
#define CLR_CYAN   "\033[1;36m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_RESET  "\033[0m"

/** Print `msg` if Verbosity() ≥ lvl, decorated with the chosen color */
#define LOG(lvl, colour, msg)                                  \
  do { if (Verbosity() >= (lvl))                               \
         std::cout << colour << msg << CLR_RESET << std::endl; \
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

//==========================================================================
//  Init – one‑time module setup
//==========================================================================
int emcal_sepdCorrelator::Init(PHCompositeNode* topNode)
{
  (void) topNode;                              // silence unused‑parameter warning
  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – starting");

  out     = new TFile(Outfile.c_str(),"RECREATE");
  LOG(1, CLR_GREEN, "[Init] opened output file: " << Outfile);

  trigAna = new TriggerAnalyzer();

  /* 1) scalar QA -------------------------------------------------------- */
  LOG(1, CLR_GREEN, "[Init] booking scalar QA histograms …");
  createHistos_Data();

  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – done");
  return Fun4AllReturnCodes::EVENT_OK;
}

int emcal_sepdCorrelator::InitRun(PHCompositeNode* topNode)
{
  if (m_mapsBooked) return Fun4AllReturnCodes::EVENT_OK;   // already done

  LOG(1, CLR_GREEN, "[InitRun] geometry is present – booking hit‑maps …");
  bookShapeHitMaps(topNode);
  m_mapsBooked = true;

  return Fun4AllReturnCodes::EVENT_OK;
}


//==========================================================================
//  bookShapeHitMaps – hex (MBD) & polar (sEPD) hit‑maps, one per trigger
//==========================================================================
void emcal_sepdCorrelator::bookShapeHitMaps(PHCompositeNode* topNode)
{
  auto* mbdg = findNode::getClass<MbdGeom>(topNode,"MbdGeom");
  auto* epdg = findNode::getClass<EpdGeom>(topNode,"TOWERGEOM_EPD");
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
    H["h_MBD_Hitmap_South_"+trig] = makeMbdHitmap("h_MBD_Hitmap_South_"+trig,
                                                    mbdg, 0);
    H["h_MBD_Hitmap_North_"+trig] = makeMbdHitmap("h_MBD_Hitmap_North_"+trig,
                                                    mbdg, 1);

    /* sEPD */
    H["h_sEPD_Hitmap_South_"+trig] = makeEpdHitmap("h_sEPD_Hitmap_South_"+trig,epdg,0);
    H["h_sEPD_Hitmap_North_"+trig] = makeEpdHitmap("h_sEPD_Hitmap_North_"+trig,epdg,1);

    /* EMCal map (256 × 96) */
    H["h_EMC_EtaPhiMap_"+trig] =
          new TH2F(("h_EMC_EtaPhiMap_"+trig).c_str(),
                   "CEMC tower map;#phi (0–255);#eta (0–95)",
                   256,0,256, 96,0,96);

    /* IHCal map ( 64 × 24) */
    H["h_IHCAL_EtaPhiMap_"+trig] =
          new TH2F(("h_IHCAL_EtaPhiMap_"+trig).c_str(),
                   "IHCAL tower map;#phi (0–63);#eta (0–23)",
                   64,0,64, 24,0,24);

    /* OHCal map ( 64 × 24) */
    H["h_OHCAL_EtaPhiMap_"+trig] =
          new TH2F(("h_OHCAL_EtaPhiMap_"+trig).c_str(),
                   "OHCAL tower map;#phi (0–63);#eta (0–23)",
                   64,0,64, 24,0,24);
  }
  out->cd();
}


//–––– helper: tower + cluster spectra ––––––––––––––––––––––––––––––––––––
void emcal_sepdCorrelator::bookTowerAndClusterQA(const std::string& trig,
                                                 HistMap& H)
{
  const int nbE = 200; const double eMax = 50.;
  for (auto& ci : m_caloInfo)
  {
    const std::string label = std::get<2>(ci);
    H["h_towerE_"+label] = new TH1F(("h_towerE_"+label+"_"+trig).c_str(),
                                    (label+" tower E;E [GeV]").c_str(),
                                    nbE,0,eMax);
    if (label=="CEMC")
      H["h_clusterE_EMC"] = new TH1F(("h_clusterE_EMC_"+trig).c_str(),
                                     "EMC cluster E;E [GeV]",nbE,0,eMax);
  }
}

//–––– helper: single‑variable charge spectra –––––––––––––––––––––––––––––
void emcal_sepdCorrelator::bookChargeQA(const std::string& trig, HistMap& H)
{
  const int nbQ = 400; const double qMax = 600.;
  H["h_towerQ_SEPD"] = new TH1F(("h_towerQ_SEPD_"+trig).c_str(),
                                "sEPD tower charge;Q [ADC]",nbQ,0,qMax);
  H["h_charge_MBD"]  = new TH1F(("h_charge_MBD_"+trig).c_str(),
                                "MBD PMT charge sum;Q [ADC]",nbQ,0,qMax);
}

//–––– helper: ΣE / ΣQ correlations –––––––––––––––––––––––––––––––––––––––
void emcal_sepdCorrelator::bookEnergyChargeCorrel(const std::string& trig,
                                                  HistMap& H)
{
  auto book2=[&](const char* n,const char* t,
                 int nx,double x0,double x1,
                 int ny,double y0,double y1)
  { return new TH2F(n,t,nx,x0,x1,ny,y0,y1); };

  const int nC=120; const double cMax=60.;
  const int nE=120; const double eMax=60.;

  H["h_SEPD_vs_CEMC"]=book2(("h_SEPD_vs_CEMC_"+trig).c_str(),"sEPD Q vs CEMC ΣE",
                            nC,0,cMax,nE,0,eMax);
  H["h_SEPD_vs_IHCAL"]=book2(("h_SEPD_vs_IHCAL_"+trig).c_str(),"sEPD Q vs IHCAL ΣE",
                             nC,0,cMax,nE,0,eMax);
  H["h_SEPD_vs_OHCAL"]=book2(("h_SEPD_vs_OHCAL_"+trig).c_str(),"sEPD Q vs OHCAL ΣE",
                             nC,0,cMax,nE,0,eMax);
  H["h_SEPD_vs_MBD"] =book2(("h_SEPD_vs_MBD_"+trig).c_str(),"sEPD Q vs MBD ΣQ",
                             nC,0,cMax,nC,0,cMax);

  H["h_MBD_vs_CEMC"] =book2(("h_MBD_vs_CEMC_"+trig).c_str(),"MBD ΣQ vs CEMC ΣE",
                             nC,0,cMax,nE,0,eMax);
  H["h_MBD_vs_IHCAL"]=book2(("h_MBD_vs_IHCAL_"+trig).c_str(),"MBD ΣQ vs IHCAL ΣE",
                             nC,0,cMax,nE,0,eMax);
  H["h_MBD_vs_OHCAL"]=book2(("h_MBD_vs_OHCAL_"+trig).c_str(),"MBD ΣQ vs OHCAL ΣE",
                             nC,0,cMax,nE,0,eMax);
    
  // --- NEW per-arm maps ---------------------------------------------------
  H["h_SEPD_S_vs_CEMC_South"] = book2(("h_SEPD_S_vs_CEMC_South_"+trig).c_str(),
                                        "ΣQ_{sEPD South}  vs  ΣEₜ_{CEMC η<0}",
                                        nC,0,cMax, nE,0,eMax);
  H["h_SEPD_N_vs_CEMC_North"] = book2(("h_SEPD_N_vs_CEMC_North_"+trig).c_str(),
                                        "ΣQ_{sEPD North}  vs  ΣEₜ_{CEMC η>0}",
                                        nC,0,cMax, nE,0,eMax);
}

//–––– helper: π0 invariant‑mass spectra ––––––––––––––––––––––––––––––––––
void emcal_sepdCorrelator::bookPi0MassSpectra(const std::string& trig,
                                              HistMap& H)
{
  const int nM=150; const double mMax=1.5;
  for (auto pt : m_ptBins)
    for (float Emin: m_minClusE)
      for (float chi: m_chi2Cuts)
        for (float a: m_asymCuts)
        {
          const std::string k = invKey(pt.first,pt.second,Emin,chi,a)+"_"+trig;
          H[k] = new TH1F(k.c_str(),"m_{#gamma#gamma};GeV/c^{2}",nM,0,mMax);
        }
    
    
    // --- NEW: inclusive (pT‑independent) spectra ---------------------------
  for (float Emin: m_minClusE)
    for (float chi: m_chi2Cuts)
      for (float a  : m_asymCuts)
      {
          const std::string k = invKey(-1,-1,Emin,chi,a)+"_"+trig; // pt = ‑1 : sentinel
          H[k] = new TH1F(k.c_str(),"m_{#gamma#gamma};GeV/c^{2}",nM,0,mMax);
      }
}

// -------------------------------------------------------------------------
// bookEventPlaneCentralityQA – all histograms used below             EP-CENT
// -------------------------------------------------------------------------
void emcal_sepdCorrelator::bookEventPlaneCentralityQA(const std::string& trig,
                                                      HistMap& H)
{
  /* --- 1. plain ΣQ --------------------------------------------------- */
  H["h_Qsum_MBD"]  = new TH1F(("h_Qsum_MBD_"+trig).c_str(),
                              "MBD ΣQ;ΣQ_{MBD} [ADC]",  600,0,1200);
  H["h_Qsum_sEPD"] = new TH1F(("h_Qsum_sEPD_"+trig).c_str(),
                              "sEPD ΣQ;ΣQ_{sEPD} [ADC]",600,0,1200);

  /* --- 2. cross‑detector map ----------------------------------------- */
  H["h_Qsum_MBD_vs_sEPD"] = new TH2F(("h_Qsum_MBD_vs_sEPD_"+trig).c_str(),
                                     "ΣQ_{MBD} vs ΣQ_{sEPD};ΣQ_{MBD};ΣQ_{sEPD}",
                                     300,0,1200, 300,0,1200);

  /* --- 3. event‑plane distributions ---------------------------------- */
  H["h_Psi2_sEPD"]   = new TH1F(("h_Psi2_sEPD_"+trig).c_str(),
                                "sEPD Ψ_{2};Ψ_{2} [rad]", 120,-TMath::Pi(),TMath::Pi());
  H["h_Psi2_res_vs_Qsum"] = new TProfile(("h_Psi2_res_vs_Qsum_"+trig).c_str(),
                                  "cos 2(Ψ_{N}-Ψ_{S}) vs ΣQ_{sEPD};ΣQ_{sEPD};⟨cos2ΔΨ⟩",
                                  12,0,1200,"s");
}


//==========================================================================
//  createHistos_Data  – delegate to four tiny helpers
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

    bookTowerAndClusterQA  (trig, H);
    bookChargeQA           (trig, H);
    bookEnergyChargeCorrel (trig, H);
    bookPi0MassSpectra     (trig, H);
    bookEventPlaneCentralityQA(trig, H);   // <<< EP‑CENT QA
      
    out->cd();
  }
}



//==========================================================================
//  process_event  – orchestration only
//==========================================================================
int emcal_sepdCorrelator::process_event(PHCompositeNode* topNode)
{
  ++event_count;
  PROGRESS("[event " << std::setw(9) << event_count << "]");

  if (!fetchNodes(topNode)) return Fun4AllReturnCodes::ABORTEVENT;

  trigAna->decodeTriggers(topNode);

  // ─── interrogate every configured trigger ────────────────
  std::vector<std::string> act;                 // active trigger keys
  if (Verbosity() >= 3)
    std::cout << CLR_BLUE << "    Trigger status:\n";

  for (auto& kv : triggerNameMap)
  {
    const std::string &bitname = kv.first;      // e.g. "MBD N&S >= 2"
    const std::string &key     = kv.second;     // e.g. "MBD_NandS_geq_2"
    ++m_trigStat[key].tested;

    const bool fired = trigAna->didTriggerFire(bitname);
    if (fired) { act.push_back(key); ++m_trigStat[key].fired; }

    if (Verbosity() >= 3)
      std::cout << "      • " << std::left << std::setw(25) << bitname
                << " → " << (fired ? CLR_GREEN "FIRED" : CLR_YELLOW "–")
                << CLR_RESET << '\n';
  }

  if (act.empty()) {
    ++m_evtNoTrig;
    if (Verbosity() >= 3)
      std::cout << CLR_YELLOW
                << "      → event skipped: no configured trigger fired\n"
                << CLR_RESET;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  doCaloQA(act);
  doSepdQA(act);
  doMbdQA (act);
  doPi0QA (act);
  fillCorrelations(act);
    
  return Fun4AllReturnCodes::EVENT_OK;
}

bool emcal_sepdCorrelator::fetchNodes(PHCompositeNode* top)
{
  /* ── grab the primary vertex ───────────────────────────────────────── */
  GlobalVertexMap* vmap = findNode::getClass<GlobalVertexMap>(top,"GlobalVertexMap");

  m_vtx = nullptr; m_vx = m_vy = m_vz = 0.;

  if (!vmap) {
    LOG(1, CLR_YELLOW, "  – GlobalVertexMap node **missing** → skip event");
    return false;
  }
  if (vmap->empty()) {
    LOG(2, CLR_YELLOW, "  – GlobalVertexMap is **empty** → skip event");
    return false;
  }

  m_vtx = vmap->begin()->second;          // cache the pointer
  if (!m_vtx) {                           // extra safety
    LOG(1, CLR_YELLOW, "  – vertex pointer null → skip event");
    return false;
  }

  m_vx = m_vtx->get_x();
  m_vy = m_vtx->get_y();
  m_vz = m_vtx->get_z();

  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut) {
    LOG(2, CLR_YELLOW, "  – |vz| = " << std::fabs(m_vz)
                                     << " cm exceeds cut (" << m_vzCut << ") → skip");
    return false;
  }
    
  /* calorimeter nodes -------------------------------------------------- */
  m_calo.clear();
  for (const auto& ci : m_caloInfo)
  {
    const std::string node = std::get<0>(ci),
                      geo  = std::get<1>(ci),
                      lbl  = std::get<2>(ci);

    auto* tw = findNode::getClass<TowerInfoContainer>   (top, node);
    auto* ge = findNode::getClass<RawTowerGeomContainer>(top, geo);
    if (!tw || !ge)
    { LOG(2, CLR_YELLOW, "  – missing " << lbl << " nodes → skip"); return false; }

    m_calo[lbl] = { tw, ge, 0. };
  }

  /* remaining detectors ------------------------------------------------ */

  // --- sEPD calibrated tower container ---------------------------------
  m_sepd = findNode::getClass<TowerInfoContainer>(top, "TOWERINFO_CALIB_SEPD");
  if (!m_sepd)
  {
     LOG(2, CLR_YELLOW, "  – node \"TOWERINFO_CALIB_SEPD\" **missing**");
  }

  // --- MBD PMT hits -----------------------------------------------------
  m_mbdpmts = findNode::getClass<MbdPmtContainer>(top, "MbdPmtContainer");
  if (!m_mbdpmts)
  {
      LOG(2, CLR_YELLOW, "  – node \"MbdPmtContainer\" **missing**");
  }

  // --- MBD offline geometry (RUN node) ---------------------------------
  m_mbdgeom = findNode::getClass<MbdGeom>(top, "MbdGeom");
  if (!m_mbdgeom)
  {
      LOG(2, CLR_YELLOW, "  – node \"MbdGeom\" **missing** (MBD geometry)");
  }

  // --- sEPD offline geometry (RUN node) --------------------------------
  m_epdgeom = findNode::getClass<EpdGeom>(top, "TOWERGEOM_EPD");
  if (!m_epdgeom)
  {
      LOG(2, CLR_YELLOW, "  – node \"TOWERGEOM_EPD\" **missing** (sEPD geometry)");
  }

  // --- EMC cluster container -------------------------------------------
  m_clus = findNode::getClass<RawClusterContainer>(top, "CLUSTERINFO_CEMC");

  // ---------------------------------------------------------------------
  // Require BOTH calibrated tower data and corresponding detector geometries
  // ---------------------------------------------------------------------
  const bool ok_sepd = (m_sepd && m_epdgeom);
  const bool ok_mbd  = (m_mbdpmts && m_mbdgeom);

  if (!ok_sepd || !ok_mbd)
  {
      LOG(2, CLR_YELLOW, "  – missing mandatory SEPD and/or MBD nodes → skip event");
  }
  return (ok_sepd && ok_mbd);

}

// ════════════════════════════════════════════════════════════════════════
//  Geometry helpers – TH2Poly booking
// ════════════════════════════════════════════════════════════════════════
TH2Poly* emcal_sepdCorrelator::makeMbdHitmap(const std::string& name,
                                             MbdGeom* geom,
                                             int       arm)   // 0 = S, 1 = N
{
  auto* h = new TH2Poly();                 // ROOT‑safe constructor
  h->SetNameTitle(name.c_str(), ";x (cm);y (cm)");

  LOG(2, CLR_BLUE,
        "    ↳ makeMbdHitmap(\"" << name << "\")  arm=" << (arm ? "North" : "South"));
    
  if (!geom)
  {
    LOG(1, CLR_YELLOW, "      [WARN] MbdGeom is nullptr – map left empty");
    return h;
  }

  std::size_t nAdded = 0;
  const double r = 3.0;                     // hex-radius [cm] ≈ actual tile size

  for (unsigned ip = 0; ip < 128; ++ip)     // 64 PMTs per arm → 128 total
  {
        if (geom->get_arm(ip) != arm) continue;

        const double cx = geom->get_x(ip);
        const double cy = geom->get_y(ip);
        if (std::isnan(cx) || std::isnan(cy))  continue;

        double x[7]{}, y[7]{};
        for (int k = 0; k < 6; ++k)
        {
            const double ang = TMath::Pi()/6. + k * TMath::Pi()/3.;
            x[k] = cx + r * std::cos(ang);
            y[k] = cy + r * std::sin(ang);
        }
        x[6] = x[0];  y[6] = y[0];

        h->AddBin(6, x, y);
        ++nAdded;
  }

  LOG(3, CLR_GREEN, "      → " << nAdded << " hex‑bins booked");
  return h;
}

TH2Poly* emcal_sepdCorrelator::makeEpdHitmap(const std::string& name,
                                             EpdGeom* geom,
                                             int       arm)   // 0 = S, 1 = N
{
  auto* h = new TH2Poly();
  h->SetNameTitle(name.c_str(), ";x (cm);y (cm)");
  LOG(2, CLR_BLUE, "    ↳ makeEpdHitmap(\"" << name << "\")  arm=" << (arm?"North":"South"));

  if (!geom)
  {
    LOG(1, CLR_YELLOW, "      [WARN] EpdGeom is nullptr – map left empty");
    return h;
  }

  /* one polygon per physical tile (256 tiles per arm) ------------------ */
  std::size_t nAdded = 0;

  for (unsigned tile = 0; tile < 256; ++tile)
  {
    const unsigned id = arm*256 + tile;     // packed ID expected by EpdGeom

    /* the offline geometry delivers the 4 tile corners in local order   *
     * (ix = 0…3).  Corner 4 closes the polygon.                         */
    double x[5]{}, y[5]{};
    /* tile centre returned by EpdGeom ----------------------------------- */
    const double r_cen  = geom->get_r (id);
    const double phi_cen= geom->get_phi(id);
    const double x_cen  = r_cen * std::cos(phi_cen);
    const double y_cen  = r_cen * std::sin(phi_cen);

    /* build a tiny square (≈4 cm side) around the centre ---------------- */
    const double d = 2.0;               // half‑side length  [cm]
    x[0] = x_cen - d;  y[0] = y_cen - d;
    x[1] = x_cen + d;  y[1] = y_cen - d;
    x[2] = x_cen + d;  y[2] = y_cen + d;
    x[3] = x_cen - d;  y[3] = y_cen + d;
    x[4] = x[0];       y[4] = y[0];

    /* guard against uninitialised tiles (beam-pipe hole, etc.) */
    if (std::isnan(x[0]) || std::isnan(y[0])) continue;

    h->AddBin(4, x, y);
    ++nAdded;
  }

  LOG(3, CLR_GREEN, "      → " << nAdded << " tile-bins booked");
  return h;
}

// ════════════════════════════════════════════════════════════════════════
//  doCaloQA – tower and cluster spectra
// ════════════════════════════════════════════════════════════════════════
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
            const double e = tw->get_energy();         if (e <= 0) continue;
            sumE += e; ++nHit;
            
            for (auto& t : trig)
                static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerE_" + lbl])->Fill(e);
            
            unsigned int key  = twC->encode_key(ch);                // universal helper
            unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
            unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(key);
            
            
            // ------------------------------------------------------------------
            //  NEW  :  transverse energy and arm-specific sums
            // ------------------------------------------------------------------
            RawTowerGeom* tg = ck.second.g->get_tower_geometry(key);
            const double eta_t = tg ? tg->get_eta() : 0.;    // falls back to 0 if null
            const double et    = e / std::cosh(eta_t);       // E_T = E cosh⁻¹ η

            if (lbl == "CEMC") {
              if (eta_t < 0) m_cemcEt_arm[0] += et;   // South
              else           m_cemcEt_arm[1] += et;   // North
            }
            else if (lbl == "IHCAL") {
              if (eta_t < 0) m_ihcalEt_arm[0] += et;
              else           m_ihcalEt_arm[1] += et;
            }
            else if (lbl == "OHCAL") {
              if (eta_t < 0) m_ohcalEt_arm[0] += et;
              else           m_ohcalEt_arm[1] += et;
            }

            for (auto& t : trig)
            {
              if      (lbl == "CEMC")
                static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_EMC_EtaPhiMap_"+t])
                    ->Fill(iphi, ieta, e);
              else if (lbl == "IHCAL")
                static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_IHCAL_EtaPhiMap_"+t])
                    ->Fill(iphi, ieta, e);
              else if (lbl == "OHCAL")
                static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_OHCAL_EtaPhiMap_"+t])
                    ->Fill(iphi, ieta, e);
            }
        }
        
        LOG(3, CLR_GREEN, "    " << lbl << " : " << nHit << " fired, ΣE = " << sumE);
    }
    
    /* EMC clusters ------------------------------------------------------- */
    if (!m_clus) return;
    RawClusterContainer::ConstRange cr = m_clus->getClusters();
    std::size_t nClus = std::distance(cr.first, cr.second);
    LOG(3, CLR_GREEN, "    EMC clusters : " << nClus);
    
    for (auto it = cr.first; it != cr.second; ++it)
    {
        const RawCluster* cl = it->second;
        for (auto& t : trig)
            static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_clusterE_EMC"])
            ->Fill(cl->get_energy());
    }
}

// ════════════════════════════════════════════════════════════════════════
//  doSepdQA – sEPD charge, hit-map  *and* event-plane QA          EP-CENT
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::doSepdQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doSepdQA]");

  /* ------------------------------------------------------------------ *
   * 1) running sums                                                    *
   * ------------------------------------------------------------------ */
  m_sepdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  /* Q-vectors for 2-nd harmonic (South / North) ----------------------- */
  double qxS = 0., qyS = 0.;
  double qxN = 0., qyN = 0.;

  for (unsigned ch = 0; ch < m_sepd->size(); ++ch)
  {
    auto* ti = m_sepd->get_tower_at_channel(ch); if (!ti) continue;
    const double w = ti->get_energy();           if (w <= 0) continue;

    /* build the official packed key once ------------------------------ */
    const unsigned key = TowerInfoDefs::encode_epd(ch);

    const int  arm  = TowerInfoDefs::get_epd_arm(key);   // 0 = South, 1 = North
    const double r   = m_epdgeom->get_r  (key);
    const double phi = m_epdgeom->get_phi(key);
    const double x   = r * std::cos(phi);
    const double y   = r * std::sin(phi);

    /* --- hit-map fill ------------------------------------------------ */
    const std::string hpfx = (arm==0 ? "h_sEPD_Hitmap_South_" : "h_sEPD_Hitmap_North_");
    for (auto& t : trig)
      static_cast<TH2Poly*>(qaHistogramsByTrigger[t][hpfx + t])->Fill(x, y, w);

    /* --- Q-vector accumulation -------------------------------------- */
    const double c2 = std::cos(2*phi), s2 = std::sin(2*phi);
    if (arm==0) { qxS += w*c2;  qyS += w*s2;  ++nFiredS; }
    else        { qxN += w*c2;  qyN += w*s2;  ++nFiredN; }

    m_sepdQ      += w;
    m_sepdQ_arm[arm] += w;      // NEW – arm‑specific sum
  }

  /* ------------------------------------------------------------------ *
   * 2) scalar ΣQ histogram                                             *
   * ------------------------------------------------------------------ */
  for (auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerQ_SEPD"])->Fill(m_sepdQ);

  /* ------------------------------------------------------------------ *
   * 3) event-plane angles & resolution proxy                           *
   * ------------------------------------------------------------------ */
  m_psi2_S = 0.5 * std::atan2(qyS, qxS);      // South   Ψ₂
  m_psi2_N = 0.5 * std::atan2(qyN, qxN);      // North   Ψ₂
  const double cos2dPsi = std::cos( 2*(m_psi2_N - m_psi2_S) );

  for (auto& t : trig)
  {
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Psi2_sEPD"])
        ->Fill(m_psi2_S);                     // store South Ψ₂
    auto  prof = dynamic_cast<TProfile*>( qaHistogramsByTrigger[t]["h_Psi2_res_vs_Qsum"] );
    if (prof)  prof->Fill(m_sepdQ, cos2dPsi);
  }

  /* ------------------------------------------------------------------ *
   * 4) verbose printout                                                *
   * ------------------------------------------------------------------ */
  LOG(3, CLR_GREEN, "    SEPD ΣQ = " << m_sepdQ
                     << "  (South hits: " << nFiredS
                     << ", North hits: " << nFiredN << ")");
}


// -------------------------------------------------------------------------
// fillCentralityQA – fills ΣQ histos & the correlation             EP-CENT
// -------------------------------------------------------------------------
void emcal_sepdCorrelator::fillCentralityQA(const std::vector<std::string>& trig)
{
  for (auto& t : trig)
  {
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Qsum_MBD"] )->Fill(m_mbdQ);
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_Qsum_sEPD"])->Fill(m_sepdQ);
    static_cast<TH2F*>(qaHistogramsByTrigger[t]["h_Qsum_MBD_vs_sEPD"])
        ->Fill(m_mbdQ, m_sepdQ);
  }
}

// ════════════════════════════════════════════════════════════════════════
//  doMbdQA – integrated charge + hex hit‑map
// ════════════════════════════════════════════════════════════════════════
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
          static_cast<TH2Poly*>(qaHistogramsByTrigger[t][key + t])
              ->Fill(cx, cy, q);

    m_mbdQ          += q;
    m_mbdQ_arm[ m_mbdgeom->get_arm(ip) ] += q;   // NEW
    ++(m_mbdgeom->get_arm(ip) ? nFiredN : nFiredS);
  }

  for (auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_charge_MBD"])->Fill(m_mbdQ);

  LOG(3, CLR_GREEN, "    MBD ΣQ = " << m_mbdQ
                     << "  (South PMTs: " << nFiredS
                     << ", North PMTs: " << nFiredN << ")");
    
  fillCentralityQA(trig);
}

// ════════════════════════════════════════════════════════════════════════
//  doPi0QA – γγ invariant‑mass spectra with live progress feedback
//           and minor speed‑ups that preserve all physics logic
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::doPi0QA(const std::vector<std::string>& trig)
{
  if (!m_clus || m_clus->size() < 2) return;

  /* ------------------------------------------------------------------ *
   * 1) build a pre‑filtered local cache of clusters                    *
   * ------------------------------------------------------------------ */
  const float  EminMin  = *std::min_element(m_minClusE.begin(),  m_minClusE.end());
  const float  chi2Max  = *std::max_element(m_chi2Cuts.begin(),  m_chi2Cuts.end());
  const float  asymMax  = *std::max_element(m_asymCuts.begin(),  m_asymCuts.end());

  struct Clu { TLorentzVector v; float E, pt, chi; };
  std::vector<Clu> cl; cl.reserve(m_clus->size());

  for (auto [it, end] = m_clus->getClusters(); it != end; ++it)
  {
    const RawCluster* c = it->second;
    if (c->get_energy() < EminMin || c->get_chi2() > chi2Max) continue;  // early veto

    const auto eVec = RawClusterUtility::GetEVec(*c, {m_vx, m_vy, m_vz});
    cl.push_back({{}, static_cast<float>(c->get_energy()),
                       static_cast<float>(eVec.perp()),
                       static_cast<float>(c->get_chi2())});

    cl.back().v.SetPtEtaPhiE(cl.back().pt,
                             eVec.pseudoRapidity(),
                             eVec.phi(),
                             cl.back().E);
  }
  if (cl.size() < 2) return;               // nothing left

  /* ------------------------------------------------------------------ *
   * 2) announce workload                                               *
   * ------------------------------------------------------------------ */
  const std::size_t totalPairs = (cl.size() * (cl.size() - 1)) / 2;
  std::cout << CLR_CYAN << "    [doPi0QA] will analyse "
            << totalPairs << " cluster pairs" << CLR_RESET << std::endl;

  /* ------------------------------------------------------------------ *
   * 3) parallel pair scan                                              *
   * ------------------------------------------------------------------ */
  constexpr std::size_t reportEvery = 50'000;
  std::atomic<std::size_t> pairCnt{0};

  /* local (thread‑private) statistics & histogram pointers ------------ */
  #ifdef _OPENMP
    #pragma omp parallel default(shared)
  #endif
  {
    std::map<std::string,CutStat> evtStatLocal;

    #ifdef _OPENMP
        #pragma omp for schedule(dynamic,256)
    #endif
    for (std::size_t idx = 0; idx < cl.size() * (cl.size() - 1) / 2; ++idx)
    {
      std::size_t i = static_cast<std::size_t>(
          (std::sqrt(8.0 * idx + 1) - 1) / 2);       // invert triangular index
      std::size_t j = idx - i * (i + 1) / 2 + i + 1;

      const Clu &c1 = cl[i], &c2 = cl[j];
      const float asym = std::fabs(c1.E - c2.E) / (c1.E + c2.E);
      if (asym > asymMax) { ++pairCnt; continue; }   // global asym veto

      /* --- inclusive block ---------------------------------------- */
      for (float Emin : m_minClusE)
      for (float chiMx : m_chi2Cuts)
      for (float aMx : m_asymCuts)
      {
        const std::string keyInc = statKey(-1,-1,Emin,chiMx,aMx);
        auto &st = evtStatLocal[keyInc]; ++st.tested;

        if (c1.E < Emin || c2.E < Emin)        continue;
        if (c1.chi > chiMx || c2.chi > chiMx)  continue;
        if (asym   > aMx)                      continue;

        const float mInv = (c1.v + c2.v).M();
        for (auto &t : trig)
          static_cast<TH1F*>( qaHistogramsByTrigger[t]
                    [invKey(-1,-1,Emin,chiMx,aMx)+"_"+t] )->Fill(mInv);
        ++st.passed;
      }

      /* --- pT‑binned block ---------------------------------------- */
      for (const auto &pb : m_ptBins)
      {
        const float ptLo = pb.first, ptHi = pb.second;
        const std::string kAll = statKey(ptLo,ptHi,0,0,0);
        auto &stAll = evtStatLocal[kAll]; ++stAll.tested;

        if (c1.pt < ptLo || c1.pt >= ptHi ||
            c2.pt < ptLo || c2.pt >= ptHi) continue;

        for (float Emin : m_minClusE)
        for (float chiMx : m_chi2Cuts)
        for (float aMx  : m_asymCuts)
        {
          const std::string key = statKey(ptLo,ptHi,Emin,chiMx,aMx);
          auto &st = evtStatLocal[key]; ++st.tested;

          if (c1.E < Emin || c2.E < Emin)      continue;
          if (c1.chi > chiMx || c2.chi > chiMx)continue;
          if (asym   > aMx)                    continue;

          const float mInv = (c1.v + c2.v).M();
          for (auto &t : trig)
            static_cast<TH1F*>( qaHistogramsByTrigger[t]
                    [invKey(ptLo,ptHi,Emin,chiMx,aMx)+"_"+t] )->Fill(mInv);
          ++st.passed;
        }
      }

      /* ---- live progress (single thread) -------------------------- */
    #ifdef _OPENMP
      if ((++pairCnt % reportEvery) == 0 && omp_get_thread_num() == 0)
    #else
      if ((++pairCnt % reportEvery) == 0)
    #endif
      {
        const double pct = 100.0 * static_cast<double>(pairCnt) /
                           static_cast<double>(totalPairs);
        std::cout << CLR_CYAN << "    [doPi0QA] processed "
                  << pairCnt << " / " << totalPairs
                  << " (" << std::fixed << std::setprecision(1) << pct << "%)\r"
                  << CLR_RESET << std::flush;
      }
    }   // omp for

    /* --- merge thread‑local statistics ----------------------------- */
  #ifdef _OPENMP
    #pragma omp critical
  #endif
    {
      for (auto &kv : evtStatLocal)
      {
        m_evtStat[kv.first].tested += kv.second.tested;
        m_evtStat[kv.first].passed += kv.second.passed;
      }
    }
  }     // omp parallel

  /* clear the “\r” progress line and print final summary -------------- */
  std::cout << CLR_CYAN << "    [doPi0QA] finished "
            << totalPairs << " / " << totalPairs << " (100.0%)            "
            << CLR_RESET << std::endl;

  /* ------------------------------------------------------------------ *
   * 4) optional per‑event summary (unchanged)                          *
   * ------------------------------------------------------------------ */
  if (Verbosity() >= 2)
  {
    std::cout << CLR_GREEN
              << "    π0‑QA summary (this event)\n"
              << "    cut‑key                                    "
                 "tested   passed   eff[%]\n"
              << "    -----------------------------------------------------------\n";
    for (const auto& kv : m_evtStat)
    {
      const auto& s = kv.second;
      const double eff = s.tested ? 100. * s.passed / s.tested : 0.;
      std::cout << "    " << std::left << std::setw(40) << kv.first
                << std::right << std::setw(8) << s.tested
                << std::setw(9) << s.passed
                << std::setw(9) << std::fixed << std::setprecision(1) << eff
                << '\n';
      m_totStat[kv.first].tested += s.tested;
      m_totStat[kv.first].passed += s.passed;
    }
    std::cout << CLR_RESET;
  }
}



// ════════════════════════════════════════════════════════════════════════
//  fillCorrelations – ΣE / ΣQ detector‑level correlations
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::fillCorrelations(const std::vector<std::string>& trig)
{
  const double cemc = m_calo["CEMC" ].sumE;
  const double ihc  = m_calo["IHCAL"].sumE;
  const double ohc  = m_calo["OHCAL"].sumE;

  // ------------------------------------------------------------------
  //  side-matched correlations  (transverse energy in calorimeters)
  // ------------------------------------------------------------------
    for (auto& t : trig)
    {
      auto& H = qaHistogramsByTrigger[t];

      // --- South (η<0) ----------------------------------------------------
      static_cast<TH2F*>(H["h_SEPD_S_vs_CEMC_South"])
          ->Fill(m_sepdQ_arm[0], m_cemcEt_arm[0]);

      // --- North (η>0) ----------------------------------------------------
      static_cast<TH2F*>(H["h_SEPD_N_vs_CEMC_North"])
          ->Fill(m_sepdQ_arm[1], m_cemcEt_arm[1]);

      // keep the old global plots for cross-checks
      static_cast<TH2F*>(H["h_SEPD_vs_CEMC"])
          ->Fill(m_sepdQ, m_calo["CEMC"].sumE);
      static_cast<TH2F*>(H["h_SEPD_vs_IHCAL"])
          ->Fill(m_sepdQ, m_calo["IHCAL"].sumE);
      static_cast<TH2F*>(H["h_SEPD_vs_OHCAL"])
          ->Fill(m_sepdQ, m_calo["OHCAL"].sumE);
      static_cast<TH2F*>(H["h_SEPD_vs_MBD"])
          ->Fill(m_sepdQ, m_mbdQ);

      static_cast<TH2F*>(H["h_MBD_vs_CEMC"])
          ->Fill(m_mbdQ, m_calo["CEMC"].sumE);
      static_cast<TH2F*>(H["h_MBD_vs_IHCAL"])
          ->Fill(m_mbdQ, m_calo["IHCAL"].sumE);
      static_cast<TH2F*>(H["h_MBD_vs_OHCAL"])
          ->Fill(m_mbdQ, m_calo["OHCAL"].sumE);
    }


  LOG(3, CLR_BLUE, "  [fillCorrelations]  ΣE(CEMC)=" << cemc
                         << "  ΣE(IHCAL)=" << ihc
                         << "  ΣE(OHCAL)=" << ohc);
}

// ════════════════════════════════════════════════════════════════════════
//  trivial helpers
// ════════════════════════════════════════════════════════════════════════
std::string emcal_sepdCorrelator::invKey(float ptLo,float ptHi,
                                         float e,float chi,float a)
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1)
      << "mInv_pt" << ptLo << "to" << ptHi
      << "_E"   << e
      << "_chi" << chi
      << "_asy" << a;
  return oss.str();
}


std::string emcal_sepdCorrelator::statKey(float ptLo,float ptHi,
                                          float Emin,float chi,float a)
{
  std::ostringstream o;
  o<<std::fixed<<std::setprecision(1);
  if(ptLo<0) o<<"allPt";
  else       o<<"pt"<<ptLo<<"to"<<ptHi;
  o<<"_E"<<Emin<<"_chi"<<chi<<"_asy"<<a;
  return o.str();
}
//==========================================================================
//  ResetEvent / End / Reset  (unchanged – keep your original)
//==========================================================================
// ‑‑‑ in emcal_sepdCorrelator.cc  (completely replace the old method) ‑‑‑
int emcal_sepdCorrelator::ResetEvent(PHCompositeNode*)
{
  //--------------------------------------------------------------------
  // wipe every per‑event cache so NOTHING bleeds into the next event
  //--------------------------------------------------------------------
  std::fill(std::begin(m_sepdQ_arm),  std::end(m_sepdQ_arm),  0.);
  std::fill(std::begin(m_mbdQ_arm),   std::end(m_mbdQ_arm),   0.);
  std::fill(std::begin(m_cemcEt_arm), std::end(m_cemcEt_arm), 0.);
  std::fill(std::begin(m_ihcalEt_arm),std::end(m_ihcalEt_arm),0.);
  std::fill(std::begin(m_ohcalEt_arm),std::end(m_ohcalEt_arm),0.);

  m_sepdQ   = 0.;
  m_mbdQ    = 0.;
  m_psi2_S  = 0.;
  m_psi2_N  = 0.;
  for (auto& kv : m_calo) kv.second.sumE = 0.;     // CEMC/IHCAL/OHCAL

  m_evtStat.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

int emcal_sepdCorrelator::Reset     (PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }

int emcal_sepdCorrelator::End(PHCompositeNode*)
{
  COUT_BLUE("Writing output…");
  if(!out||!out->IsOpen()) return Fun4AllReturnCodes::ABORTEVENT;
  for(auto& tk:qaHistogramsByTrigger){
    TDirectory* d=out->mkdir(tk.first.c_str()); d->cd();
    for(auto& hk:tk.second){
      TH1* h=dynamic_cast<TH1*>(hk.second);
      if(h && h->GetEntries()) h->Write();
    }
    out->cd();
  }
  out->Write();
  out->Close();
  delete out;
  out = nullptr;

  /* ------------ console table (only if Verbosity() > 0) -------------- */
  if (Verbosity() > 0)
  {
      std::cout << "\n\033[1mHistogram summary\033[0m\n"
                << "\033[1mTrigger                        │ Histogram                           │  Entries\033[0m\n"
                << "-------------------------------------------------------------------------------\n";

      for (const auto& tk : qaHistogramsByTrigger)
        for (const auto& hk : tk.second)
        {
          TH1* h = dynamic_cast<TH1*>(hk.second);
          if (!h) continue;
          std::cout << std::left  << std::setw(30) << tk.first << " │ "
                    << std::setw(32) << hk.first   << " │ "
                    << std::right << std::setw(10) << static_cast<Long64_t>(h->GetEntries())
                    << '\n';
        }

      std::cout << "-------------------------------------------------------------------------------\n";
  }

  // ---------- π0‑QA run summary ----------------------------------
  if(Verbosity()>0 && !m_totStat.empty())
  {
      std::cout << "\n\033[1mπ0‑QA cut summary (all events)\033[0m\n"
                << "\033[1mcut‑key                                    │ "
                   "tested        passed    eff[%]\033[0m\n"
                << "--------------------------------------------------------------------------\n";
      for(const auto& kv : m_totStat)
      {
        const auto& s = kv.second;
        const double eff = s.tested ? 100.*s.passed/s.tested : 0.;
        std::cout << std::left << std::setw(44) << kv.first << " │ "
                  << std::right<< std::setw(12)<< s.tested
                  << std::setw(12)<< s.passed
                  << std::setw(9) << std::fixed << std::setprecision(2) << eff
                  << '\n';
      }
      std::cout << "--------------------------------------------------------------------------\n";
  }
  COUT_BLUE("Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}

void emcal_sepdCorrelator::Print(const std::string&) const { /* nothing */ }
