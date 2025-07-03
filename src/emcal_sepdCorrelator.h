// Tell Emacs this is C++   -*- C++ -*-
#ifndef EMCALSEPDCORRELATOR_H
#define EMCALSEPDCORRELATOR_H
//==========================================================================
//  sPHENIX EMCal × sEPD × MBD correlator – headers
//  Author:  <your name>          (world‑class clean‑room version)
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
#include <calobase/RawTowerGeomContainer.h>
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
  static TH2Poly* makeMbdHitmap(const std::string&, MbdGeom*, int arm);
  static TH2Poly* makeEpdHitmap(const std::string&,            int arm);

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
  double      m_vzCut      = 30.;      // [cm]
  bool        m_useVzCut   = true;
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
    
  double m_sepdQ = 0., m_mbdQ = 0., m_vz = 0.;

  // ––– utility -----------------------------------------------------------
  static std::string invKey(float ptLo,float ptHi,
                            float minE,float maxChi,float maxAsy);
};
//==========================================================================

#endif  // EMCALSEPDCORRELATOR_H
