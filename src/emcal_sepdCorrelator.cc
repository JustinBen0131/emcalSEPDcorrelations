//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator – implementation
//==========================================================================

#include "emcal_sepdCorrelator.h"

//––– Fun4All / PHOOL -------------------------------------------------------
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>

//––– Other ROOT headers ----------------------------------------------------
#include <TDirectory.h>
#include <TSystem.h>
#include <TMath.h>
#include <globalvertex/GlobalVertex.h>   // full definition of GlobalVertex
#include <calobase/TowerInfo.h>   // full definition of TowerInfo
#include <iostream>

//–––––––– helpers ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#define COUT_BLUE(MSG) \
    if (Verbosity() > 0) std::cout << "\033[1;34m" << MSG << "\033[0m\n"


#define CLR_BLUE   "\033[1;34m"
#define CLR_GREEN  "\033[1;32m"
#define CLR_YELLOW "\033[1;33m"
#define CLR_RESET  "\033[0m"

/** Print `msg` if Verbosity() ≥ lvl, decorated with the chosen color */
#define LOG(lvl, colour, msg)                                  \
  do { if (Verbosity() >= (lvl))                               \
         std::cout << colour << msg << CLR_RESET << std::endl; \
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

//==========================================================================
//  Init – one‑time module setup
//==========================================================================
int emcal_sepdCorrelator::Init(PHCompositeNode* topNode)
{
  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – starting");

  out     = new TFile(Outfile.c_str(),"RECREATE");
  LOG(1, CLR_GREEN, "[Init] opened output file: " << Outfile);

  trigAna = new TriggerAnalyzer();

  /* 1) scalar QA -------------------------------------------------------- */
  LOG(1, CLR_GREEN, "[Init] booking scalar QA histograms …");
  createHistos_Data();

  /* 2) detector‑shape hit‑maps ----------------------------------------- */
  LOG(1, CLR_GREEN, "[Init] booking detector hit‑maps …");
  bookShapeHitMaps(topNode);

  LOG(1, CLR_BLUE, "[Init] emcal_sepdCorrelator – done");
  return Fun4AllReturnCodes::EVENT_OK;
}

//==========================================================================
//  bookShapeHitMaps – hex (MBD) & polar (sEPD) hit‑maps, one per trigger
//==========================================================================
void emcal_sepdCorrelator::bookShapeHitMaps(PHCompositeNode* topNode)
{
  auto* mbdg = findNode::getClass<MbdGeom>(topNode,"MbdGeom");
  auto* epdg = findNode::getClass<EpdGeom>(topNode,"EpdGeom");
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
    H["h_sEPD_Hitmap_South_"+trig] = makeEpdHitmap("h_sEPD_Hitmap_South_"+trig,0);
    H["h_sEPD_Hitmap_North_"+trig] = makeEpdHitmap("h_sEPD_Hitmap_North_"+trig,1);
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

    out->cd();
  }
}



//==========================================================================
//  process_event  – orchestration only
//==========================================================================
int emcal_sepdCorrelator::process_event(PHCompositeNode* topNode)
{
  ++event_count;

  if (!fetchNodes(topNode)) return Fun4AllReturnCodes::ABORTEVENT;

  trigAna->decodeTriggers(topNode);
  std::vector<std::string> act;
  for (auto& kv : triggerNameMap)
    if (trigAna->didTriggerFire(kv.first)) act.push_back(kv.second);
  if (act.empty()) return Fun4AllReturnCodes::ABORTEVENT;

  doCaloQA(act);
  doSepdQA(act);
  doMbdQA (act);
  doPi0QA (act);
  fillCorrelations(act);
    
  return Fun4AllReturnCodes::EVENT_OK;
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

//––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
TH2Poly* emcal_sepdCorrelator::makeEpdHitmap(const std::string& name,
                                             int /*arm*/)
{
  auto* h = new TH2Poly();
  h->SetNameTitle(name.c_str(), ";x (cm);y (cm)");
    
  LOG(2, CLR_BLUE, "    ↳ makeEpdHitmap(\"" << name << "\")");

  constexpr int    NR   = 16;        // radial rings
  constexpr int    NPhi = 12;        // azimuthal sectors
  constexpr double Rmax = 15.5;      // cm
  const double     dR   = Rmax / NR;

  std::size_t nAdded = 0;

  for (int ir = 0; ir < NR; ++ir)
  {
    const double r0 = ir * dR;
    const double r1 = (ir + 1) * dR;

    for (int ip = 0; ip < NPhi; ++ip)
    {
      const double p0 = ip       * TMath::TwoPi() / NPhi;
      const double p1 = (ip + 1) * TMath::TwoPi() / NPhi;

      double x[5]{}, y[5]{};
      x[0] = r0 * std::cos(p0);  y[0] = r0 * std::sin(p0);
      x[1] = r1 * std::cos(p0);  y[1] = r1 * std::sin(p0);
      x[2] = r1 * std::cos(p1);  y[2] = r1 * std::sin(p1);
      x[3] = r0 * std::cos(p1);  y[3] = r0 * std::sin(p1);
      x[4] = x[0];               y[4] = y[0];

      h->AddBin(4, x, y);
      ++nAdded;
    }
  }

  LOG(3, CLR_GREEN, "      → " << nAdded << " polar‑bins booked");
  return h;
}

// ════════════════════════════════════════════════════════════════════════
//  fetchNodes – fills per‑event caches, applies vertex cut
// ════════════════════════════════════════════════════════════════════════
bool emcal_sepdCorrelator::fetchNodes(PHCompositeNode* top)
{
  LOG(3, CLR_BLUE, "[fetchNodes] grabbing event nodes");

  auto* vmap = findNode::getClass<GlobalVertexMap>(top,"GlobalVertexMap");
  if (!vmap || vmap->empty())
  { LOG(2, CLR_YELLOW, "  – no vertex → skip event"); return false; }

  m_vz = vmap->begin()->second->get_z();
  if (m_useVzCut && std::fabs(m_vz) >= m_vzCut)
  { LOG(2, CLR_YELLOW, "  – |vz| = " << m_vz << " cm > cut → skip"); return false; }

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
  m_sepd     = findNode::getClass<TowerInfoContainer>(top,"TOWERS_SEPD");
  m_mbdpmts  = findNode::getClass<MbdPmtContainer>   (top,"MbdPmtContainer");
  m_mbdgeom  = findNode::getClass<MbdGeom>           (top,"MbdGeom");
  m_epdgeom  = findNode::getClass<EpdGeom>           (top,"EpdGeom");
  m_clus     = findNode::getClass<RawClusterContainer>(top,"CLUSTERINFO_CEMC");

  const bool ok = (m_sepd && m_mbdpmts && m_mbdgeom && m_epdgeom);
  if (!ok) LOG(2, CLR_YELLOW, "  – missing SEPD/MBD geometry → skip");
  return ok;
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
    }

    LOG(3, CLR_GREEN, "    " << lbl << " : " << nHit << " fired, ΣE = " << sumE);
  }

  /* EMC clusters ------------------------------------------------------- */
  if (!m_clus) return;
  std::size_t nClus = m_clus->size();
  LOG(3, CLR_GREEN, "    EMC clusters : " << nClus);

  for (const auto* cl : *m_clus)
    for (auto& t : trig)
      static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_clusterE_EMC"])
          ->Fill(cl->get_energy());
}

// ════════════════════════════════════════════════════════════════════════
//  doSepdQA – integrated charge + polar hit‑map
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::doSepdQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doSepdQA]");

  m_sepdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  for (unsigned ch = 0; ch < m_sepd->size(); ++ch)
  {
    auto* ti = m_sepd->get_tower_at_channel(ch); if (!ti) continue;
    const double q = ti->get_energy();           if (q <= 0) continue;

    const int arm  = (ch < 256 ? 0 : 1);         // 0 S / 1 N
    const int tile = ch % 256;
    double x, y, z;  m_epdgeom->GetGlobalPosition(arm, tile, x, y, z);

    const std::string key = (arm == 0 ? "h_sEPD_Hitmap_South_"
                                      : "h_sEPD_Hitmap_North_");

    for (auto& t : trig)
      static_cast<TH2Poly*>(qaHistogramsByTrigger[t][key + t])->Fill(x, y, q);

    m_sepdQ += q;
    ++(arm ? nFiredN : nFiredS);
  }

  for (auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_towerQ_SEPD"])->Fill(m_sepdQ);

  LOG(3, CLR_GREEN, "    SEPD ΣQ = " << m_sepdQ
                     << "  (South hits: " << nFiredS
                     << ", North hits: " << nFiredN << ")");
}

// ════════════════════════════════════════════════════════════════════════
//  doMbdQA – integrated charge + hex hit‑map
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::doMbdQA(const std::vector<std::string>& trig)
{
  LOG(3, CLR_BLUE, "  [doMbdQA]");

  m_mbdQ = 0.;
  std::size_t nFiredS = 0, nFiredN = 0;

  for (unsigned ip = 0; ip < m_mbdpmts->get_npmt(); ++ip)
  {
    auto* p = m_mbdpmts->get_pmt(ip);
    const double q = p->get_q(); if (q <= 0) continue;
    
    double cx = m_mbdgeom->get_x(ip);
    double cy = m_mbdgeom->get_y(ip);
    const std::string key = (p->get_arm() == 0 ? "h_MBD_Hitmap_South_"
                                                : "h_MBD_Hitmap_North_");

    for (auto& t : trig)
          static_cast<TH2Poly*>(qaHistogramsByTrigger[t][key + t])
              ->Fill(cx, cy, q);

    m_mbdQ += q;
    ++(p->get_arm() ? nFiredN : nFiredS);
  }

  for (auto& t : trig)
    static_cast<TH1F*>(qaHistogramsByTrigger[t]["h_charge_MBD"])->Fill(m_mbdQ);

  LOG(3, CLR_GREEN, "    MBD ΣQ = " << m_mbdQ
                     << "  (South PMTs: " << nFiredS
                     << ", North PMTs: " << nFiredN << ")");
}

// ════════════════════════════════════════════════════════════════════════
//  doPi0QA – γγ invariant‑mass spectra  (summary prints only)
// ════════════════════════════════════════════════════════════════════════
void emcal_sepdCorrelator::doPi0QA(const std::vector<std::string>& trig)
{
  if (!m_clus || m_clus->size() < 2) return;

  const std::size_t nClus = m_clus->size();
  LOG(3, CLR_BLUE, "  [doPi0QA] " << nClus << " clusters in event");

  struct Clu { TLorentzVector v; float E, pt, chi; };
  std::vector<Clu> cl; cl.reserve(nClus);

  for (const auto* c : *m_clus)
  {
    cl.push_back({{}, c->get_energy(), c->get_pt(), c->get_chi2()});
    cl.back().v.SetPtEtaPhiE(cl.back().pt, c->get_eta(), c->get_phi(), cl.back().E);
  }

  /* count accepted pairs per pT‑bin for a concise summary */
  std::map<std::pair<float,float>, std::size_t> pairCounter;

  for (std::size_t i = 0; i < cl.size(); ++i)
    for (std::size_t j = i + 1; j < cl.size(); ++j)
    {
      const float e1 = cl[i].E,  e2 = cl[j].E,
                  pt1 = cl[i].pt, pt2 = cl[j].pt,
                  chi1 = cl[i].chi, chi2 = cl[j].chi,
                  asym = std::fabs(e1 - e2) / (e1 + e2),
                  mInv = (cl[i].v + cl[j].v).M();

      for (auto pb : m_ptBins)
      {
        if (pt1 < pb.first || pt1 >= pb.second) continue;
        if (pt2 < pb.first || pt2 >= pb.second) continue;

        bool accepted = false;

        for (float Emin  : m_minClusE)
          if (e1 >= Emin && e2 >= Emin)
            for (float chiMax : m_chi2Cuts)
              if (chi1 < chiMax && chi2 < chiMax)
                for (float aMax : m_asymCuts)
                  if (asym < aMax)
                  {
                    const std::string base = invKey(pb.first,pb.second,
                                                    Emin,chiMax,aMax);

                    for (auto& t : trig)
                      static_cast<TH1F*>(qaHistogramsByTrigger[t][base + "_" + t])
                          ->Fill(mInv);
                    accepted = true;
                  }

        if (accepted) ++pairCounter[pb];
      }
    }

  /* concise table ------------------------------------------------------ */
  if (Verbosity() >= 2 && !pairCounter.empty())
  {
    std::cout << CLR_GREEN << "    π0 pairs accepted (per pT bin):";
    for (const auto& kv : pairCounter)
      std::cout << "  [" << kv.first.first << "–" << kv.first.second
                << " GeV] " << kv.second;
    std::cout << CLR_RESET << std::endl;
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

  for (auto& t : trig)
  {
    auto& H = qaHistogramsByTrigger[t];

    static_cast<TH2F*>(H["h_SEPD_vs_CEMC"])->Fill(m_sepdQ, cemc);
    static_cast<TH2F*>(H["h_SEPD_vs_IHCAL"])->Fill(m_sepdQ, ihc);
    static_cast<TH2F*>(H["h_SEPD_vs_OHCAL"])->Fill(m_sepdQ, ohc);
    static_cast<TH2F*>(H["h_SEPD_vs_MBD" ])->Fill(m_sepdQ, m_mbdQ);

    static_cast<TH2F*>(H["h_MBD_vs_CEMC"])->Fill(m_mbdQ, cemc);
    static_cast<TH2F*>(H["h_MBD_vs_IHCAL"])->Fill(m_mbdQ, ihc);
    static_cast<TH2F*>(H["h_MBD_vs_OHCAL"])->Fill(m_mbdQ, ohc);
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


//==========================================================================
//  ResetEvent / End / Reset  (unchanged – keep your original)
//==========================================================================
int emcal_sepdCorrelator::ResetEvent(PHCompositeNode*) { return Fun4AllReturnCodes::EVENT_OK; }
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
  out->Write(); out->Close(); delete out; out=nullptr;
  COUT_BLUE("Done.");
  return Fun4AllReturnCodes::EVENT_OK;
}

void emcal_sepdCorrelator::Print(const std::string&) const { /* nothing */ }
