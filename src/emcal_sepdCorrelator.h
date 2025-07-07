// Tell Emacs this is C++   -*- C++ -*-
#ifndef EMCALSEPDCORRELATOR_H
#define EMCALSEPDCORRELATOR_H
//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator – headers
//  ------------------------------------------------------------------
//  PUBLIC  :  unchanged Fun4All module interface.
//  PRIVATE :  book‑once helpers  |  per‑event helpers  |  caches.
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
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>        // ← new: brings in get_eta()
#include <calobase/RawClusterContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include <epd/EpdGeom.h>

//––– STL ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using HistMap = std::map<std::string,TObject*>;

//==========================================================================
//                               CLASS
//==========================================================================
class emcal_sepdCorrelator : public SubsysReco
{
 public:
  explicit emcal_sepdCorrelator(const std::string& out = "caloTreeData.root");
  ~emcal_sepdCorrelator() override = default;

  // Fun4All hooks ---------------------------------------------------------
  int  Init            (PHCompositeNode*) override;
  int  InitRun         (PHCompositeNode*) override;
  int  process_event   (PHCompositeNode*) override;
  int  ResetEvent      (PHCompositeNode*) override;
  int  End             (PHCompositeNode*) override;
  int  Reset           (PHCompositeNode*) override;
  void Print           (const std::string& what = "ALL") const override;

  // user bits -------------------------------------------------------------
  void setVerbose(int level) { Verbosity(level); }
  void setRunNumber(int r)         { m_runNumber = r;   }
  void setVzCut(double c)          { m_vzCut = std::fabs(c); }
  void enableVzCut(bool f=true)    { m_useVzCut = f;    }

 private:
  //=======================================================================
  // 1)  One‑time booking helpers
  //=======================================================================
  void createHistos_Data();                                  // main booker
  TH2Poly* makeMbdHitmap(const std::string&, MbdGeom*, int arm);
  TH2Poly* makeEpdHitmap(const std::string& name, EpdGeom* geom, int arm);

  //=======================================================================
  // 2)  Per‑event helpers (called in process_event)
  //=======================================================================
  bool fetchNodes (PHCompositeNode*);                        // guards + cache
  void doCaloQA   (const std::vector<std::string>&);
  void doSepdQA   (const std::vector<std::string>&);
  void doMbdQA    (const std::vector<std::string>&);
  void doPi0QA    (const std::vector<std::string>&);
  void fillCorrelations (const std::vector<std::string>&);

  //=======================================================================
  // 3)  Configuration & state
  //=======================================================================
  // ––– run‑wide ----------------------------------------------------------
  int         m_runNumber  = -1;
  bool        verbose      = true;
  double      m_vzCut      = 10.;      // [cm]
  bool        m_useVzCut   = true;
  const GlobalVertex* m_vtx {nullptr};
  double m_vx {0.}, m_vy {0.}, m_vz {0.};
    
  std::string Outfile;                 // ROOT output

  TFile*            out   = nullptr;
  TriggerAnalyzer*  trigAna = nullptr;
  std::size_t       event_count = 0;

  // ––– trigger map (unchanged) ------------------------------------------
  std::map<std::string,std::string> triggerNameMap {
     {"MBD N&S >= 2", "MBD_NandS_geq_2"}
  };
  std::map<std::string,
           std::map<std::string,TObject*>> qaHistogramsByTrigger;

  // ––– cut tables (unchanged) -------------------------------------------
  const std::vector<float>               m_asymCuts   {0.5f,0.7f};
  const std::vector<float>               m_chi2Cuts   {4.f};
  const std::vector<float>               m_minClusE   {1.f,2.f};
  const std::vector<std::pair<float,float>> m_ptBins {
        {2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10},
        {10,12},{12,15},{15,20},{20,30} };

  // ––– calorimeter convenience list -------------------------------------
  const std::vector<std::tuple<std::string,std::string,std::string>> m_caloInfo {
        {"TOWERINFO_CALIB_CEMC",   "TOWERGEOM_CEMC",   "CEMC"},
        {"TOWERINFO_CALIB_HCALIN", "TOWERGEOM_HCALIN", "IHCAL"},
        {"TOWERINFO_CALIB_HCALOUT","TOWERGEOM_HCALOUT","OHCAL"} };

  // ––– run‑time caches ---------------------------------------------------
  struct CaloCache { TowerInfoContainer* tw=nullptr;
                     RawTowerGeomContainer* g=nullptr;
                     double sumE=0.; };
  std::map<std::string,CaloCache>  m_calo;  // "CEMC" …
  TowerInfoContainer*  m_sepd    = nullptr;
  MbdPmtContainer*     m_mbdpmts = nullptr;
  MbdGeom*             m_mbdgeom = nullptr;
  EpdGeom*             m_epdgeom = nullptr;
  RawClusterContainer* m_clus    = nullptr;

  void bookShapeHitMaps        (PHCompositeNode* topNode);          ///< NEW
  void bookTowerAndClusterQA   (const std::string& trig, HistMap& H);
  void bookChargeQA            (const std::string& trig, HistMap& H);
  void bookEnergyChargeCorrel  (const std::string& trig, HistMap& H);
  void bookPi0MassSpectra      (const std::string& trig, HistMap& H);
  void bookEventPlaneCentralityQA (const std::string& trig, HistMap& H);
  void   fillCentralityQA         (const std::vector<std::string>& trig);
  void   fillEventPlaneQA         (const std::vector<std::string>& trig);
    
  double m_sepdQ = 0., m_mbdQ = 0.;
  double m_psi2_N = 0., m_psi2_S = 0.;        ///< Ψ₂ from North/South sEPD

  //----------------------------------------------------------------
  // --- per‑arm ( 0 = South / η < 0 , 1 = North / η > 0 ) caches
  //----------------------------------------------------------------
  double m_sepdQ_arm [2] {0., 0.};   ///< ΣQ  sEPD  (ADC counts)
  double m_mbdQ_arm  [2] {0., 0.};   ///< ΣQ  MBD   (ADC counts)
  double m_cemcEt_arm[2] {0., 0.};   ///< ΣEₜ CEMC  (GeV)
  double m_ihcalEt_arm[2]{0., 0.};   ///< ΣEₜ IHCAL (GeV)
  double m_ohcalEt_arm[2]{0., 0.};   ///< ΣEₜ OHCAL (GeV)

  // ––– utility -----------------------------------------------------------
  static std::string invKey(float ptLo,float ptHi,
                              float minE,float maxChi,float maxAsy);
  static std::string statKey(float ptLo,float ptHi,
                               float Emin,float chiMax,float aMax);
  };
  bool m_mapsBooked = false;

  struct CutStat {
  std::size_t tested  = 0;   ///< # candidate pairs tested
  std::size_t passed  = 0;   ///< # pairs that survived all cuts
  };
  /// per‑event scratch pad (cleared in ResetEvent)
  std::map<std::string,CutStat> m_evtStat;
  /// run‑wide accumulation (written in End())
  std::map<std::string,CutStat> m_totStat;


  // ───── inside the class (private or public, as you prefer) ──────────────
  struct TrigStat {               // run‑wide bookkeeping
    std::size_t tested = 0;       // events where decodeTriggers() succeeded
    std::size_t fired  = 0;       // events accepted by this trigger
  };
  std::map<std::string,TrigStat>  m_trigStat;
  std::size_t m_evtNoTrig = 0;    // events rejected because no trigger fired

  //==========================================================================
  //──────────────── mapping helpers ─────────────────────────────────────────
  /** Map EMCal tower (ieta,iphi) to  sector 0–63 */
  static inline int sector_from_idx(unsigned int ieta, unsigned int iphi)
  {
    if (iphi >= 256) return -1;
    const int base = iphi / 8;                 // 8 φ bins per sector slice
    return (ieta < 48) ? 32 + base             // bottom half
                     :           base;       // top half
  }

  /** Map tower (ieta,iphi) to  IB number 0–5 */
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


#endif  // EMCALSEPDCORRELATOR_H
