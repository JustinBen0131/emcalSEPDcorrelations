// Tell Emacs this is C++   -*- C++ -*-
#ifndef EMCALSEPDCORRELATOR_H
#define EMCALSEPDCORRELATOR_H
//==========================================================================
//  EMCal × sEPD × MBD correlator – headers
//==========================================================================

//––– Framework ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <phool/PHCompositeNode.h>

//––– ROOT base ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TLorentzVector.h>

//––– sPHENIX objects ––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawClusterContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/src_epdReco/EpdGeom.h"
#include <centrality/CentralityInfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

//––– STL ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <jetbase/JetContainer.h>

// --------------------------------------------------------------------------
//  Helper alias: maps histogram name → ROOT object*
// --------------------------------------------------------------------------
using HistMap = std::map<std::string, TObject*>;
class CentralityInfo;
class MinimumBiasInfo;
class Fun4AllHistoManager;
class PHCompositeNode;

// ==========================================================================
//  CLASS DCaloTowerCalibECLARATION
// ==========================================================================
class emcal_sepdCorrelator : public SubsysReco
{
 public:
  static constexpr std::array<std::pair<const char*, const char*>, 1> kJetRadii {{
        {"r02", "AntiKt_TowerInfo_HIRecoSeedsRaw_r02"}
  }};
      
  // ---------- construction / destruction ----------------------------------
  explicit emcal_sepdCorrelator(const std::string& out = "caloTreeData.root");
  ~emcal_sepdCorrelator() override = default;

  // ---------- Fun4All hooks ------------------------------------------------
  int Init          (PHCompositeNode*) override;
  int InitRun       (PHCompositeNode*) override;
  int process_event (PHCompositeNode*) override;
  int ResetEvent    (PHCompositeNode*) override;
  int End           (PHCompositeNode*) override;
  int Reset         (PHCompositeNode*) override;
  void Print(const std::string& what = "ALL") const override;

  // ---------- user configuration ------------------------------------------
  void setVerbose  (int level)   { Verbosity(level);        }
  void setRunNumber(int r)       { m_runNumber = r;         }
  void setVzCut    (double c)    { m_vzCut = std::fabs(c);  }
  void enableVzCut (bool f = true) { m_useVzCut = f;        }
  void setCentralityEdges(const std::vector<int>& e) { m_centEdges = e; }

 private:
  // ======================================================================
  // 1) One‑time booking helpers
  // ======================================================================
  void      createHistos_Data();                          // main booker
  TH2Poly*  makeMbdHitmap(const std::string&, const MbdGeom*, int arm);
  TH2F*     makeEpdHitmap(const std::string& name, EpdGeom* geom, int arm);
  bool      m_sepdMapReady {false};
  void      buildSepdChannelMap();  
  // ======================================================================
  // 2) Per‑event helpers (called in process_event)
  // ======================================================================
  bool fetchNodes (PHCompositeNode*);                     // guards + cache
  void doCaloQA   (const std::vector<std::string>&);
  void doSepdQA   (const std::vector<std::string>&);
  void doMbdQA    (const std::vector<std::string>&);
  void doPi0QA    (const std::vector<std::string>&);
  void fillCorrelations (const std::vector<std::string>&);

  // ======================================================================
  // 3) Configuration & run‑time state
  // ======================================================================
  // --- run‑wide -----------------------------------------------------------
  int         m_runNumber   = -1;
  bool        verbose       = true;
  double      m_vzCut       = 10.;        // [cm]
  bool        m_useVzCut    = true;
  const GlobalVertex* m_vtx {nullptr};
  double m_vx {0.}, m_vy {0.}, m_vz {0.};
  std::vector<unsigned> m_epdKey;
  std::vector<int> m_centEdges {0,10,20,30,40,50,60};
  int                         m_centBin   = -1;   // 0…99
  std::map<std::string,int>   m_centIdxCache;         // "0_10" → 0, etc.

  std::string       Outfile;              // ROOT output file name
  TFile*            out      = nullptr;
  TriggerAnalyzer*  trigAna  = nullptr;
  std::size_t       event_count = 0;

  // --- trigger bookkeeping -----------------------------------------------
    std::map<std::string, std::string> triggerNameMap {
        {"MBD N&S >= 2", "MBD_NandS_geq_2"}};
  std::map<std::string, HistMap>     qaHistogramsByTrigger;

  // --- analysis cuts ------------------------------------------------------
  const std::vector<float>               m_asymCuts {0.7f};
  const std::vector<float>               m_chi2Cuts {4.f};
  const std::vector<float>               m_minClusE {2.f};
  const std::vector<std::pair<float,float>> m_ptBins {
        {2,4},{4,6},{6,8},{8,10},{10,12},{12,15},{15,20},{20,30} };

  // --- detector lists -----------------------------------------------------
  const std::vector<std::tuple<std::string,std::string,std::string>> m_caloInfo {
        {"TOWERINFO_CALIB_CEMC_RETOWER",   "TOWERGEOM_CEMC",   "CEMC"},
        {"TOWERINFO_CALIB_HCALIN", "TOWERGEOM_HCALIN", "IHCAL"},
        {"TOWERINFO_CALIB_HCALOUT","TOWERGEOM_HCALOUT","OHCAL"} };

  // --- run‑time caches ----------------------------------------------------
  struct CaloCache {
    TowerInfoContainer*    tw = nullptr;
    RawTowerGeomContainer* g  = nullptr;
    double                 sumE = 0.;
  };
  std::map<std::string, CaloCache> m_calo;         // keyed by "CEMC"/…

  TowerInfoContainer*  m_sepd     = nullptr;
  MbdPmtContainer*     m_mbdpmts  = nullptr;
  MbdGeom*             m_mbdgeom  = nullptr;
  EpdGeom*             m_epdgeom  = nullptr;
  EventplaneinfoMap*   m_epmap    = nullptr;   // pointer to EventplaneinfoMap
  RawClusterContainer* m_clus     = nullptr;

  // --- helpers that book sets of histograms -------------------------------
  void bookShapeHitMaps           (PHCompositeNode* topNode);
  void bookTowerAndClusterQA      (const std::string& trig, HistMap& H);
  void bookChargeQA               (const std::string& trig, HistMap& H);
  void bookEnergyChargeCorrel     (const std::string& trig, HistMap& H);
  void bookPi0MassSpectra         (const std::string& trig, HistMap& H);
  void bookEventPlaneCentralityQA (const std::string& trig, HistMap& H);
  void fillCentralityQA           (const std::vector<std::string>& trig);
  void fillEventPlaneQA           (const std::vector<std::string>& trig);

  // --- per‑event scalars ---------------------------------------------------
  double m_sepdQ  = 0.,  m_mbdQ  = 0.;     // integrated charges
  double m_psi2_N = 0.,  m_psi2_S = 0.;    // event‑plane angles (rad)

  // --- per‑arm caches (0 = South / η<0, 1 = North / η>0) ------------------
  double m_sepdQ_arm  [2] {0., 0.};
  double m_mbdQ_arm   [2] {0., 0.};
  double m_cemcEt_arm [2] {0., 0.};
  double m_ihcalEt_arm[2] {0., 0.};
  double m_ohcalEt_arm[2] {0., 0.};

  // --- utility ------------------------------------------------------------
  static std::string invKey (float ptLo,float ptHi,
                             float minE,float maxChi,float maxAsy);
  static std::string statKey(float ptLo,float ptHi,
                             float Emin,float chiMax,float aMax);

  bool m_mapsBooked = false;          // hit‑maps booked in InitRun()

  struct CutStat {
    std::size_t tested  = 0;
    std::size_t failE   = 0;
    std::size_t failChi = 0;
    std::size_t failAsy = 0;
    std::size_t passed  = 0;
  };
  std::map<std::string,CutStat> m_evtStat;  // per‑event
  std::map<std::string,CutStat> m_totStat;  // run‑wide

  struct TrigStat {
    std::size_t tested = 0;   // decodeTriggers() succeeded
    std::size_t fired  = 0;   // trigger accepted
  };
  std::map<std::string,TrigStat> m_trigStat;
  std::size_t m_evtNoTrig = 0;

  void  bookJetQA (const std::string& trig, HistMap& H);
  int   doJetQA   (PHCompositeNode* topNode, const std::vector<std::string>& trig);
  float getMaxJetEt(JetContainer* jets) const;
    
  // ======================================================================
  // 4) Static mapping helpers
  // ======================================================================
  /** Map EMCal tower indices (ieta,iphi) → sector 0–63 */
  static inline int sector_from_idx(unsigned int ieta, unsigned int iphi)
  {
    if (iphi >= 256) return -1;
    const int base = iphi / 8;                 // 8 φ bins / sector slice
    return (ieta < 48) ? 32 + base : base;     // bottom vs. top half
  }

  /** Map tower indices (ieta,iphi) → inner‑barrel number 0–5 */
  static inline int ib_from_idx(unsigned int ieta, unsigned int /*iphi*/)
  {
    if      (ieta <  8) return 5;
    else if (ieta < 16) return 4;
    else if (ieta < 24) return 3;
    else if (ieta < 32) return 2;
    else if (ieta < 40) return 1;
    else if (ieta < 48) return 0;
    else if (ieta < 56) return 0;
    else if (ieta < 64) return 1;
    else if (ieta < 72) return 2;
    else if (ieta < 80) return 3;
    else if (ieta < 88) return 4;
    else if (ieta < 96) return 5;
    return -1;
  }
};

// --------------------------------------------------------------------------
#endif  // EMCALSEPDCORRELATOR_H
