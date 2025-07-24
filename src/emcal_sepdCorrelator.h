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
#include <TH3F.h> 
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
#include <epd/EpdGeom.h>
#include <epd/EpdReco.h>
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
#include <bitset>          // ← for extractTriggerBits()
#include <TVector2.h>      // ← TVector2::Phi_mpi_pi in doCaloQA()
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
  ~emcal_sepdCorrelator() override;

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
  bool m_isMinBias {false};
  // ======================================================================
  // 2) Per‑event helpers (called in process_event)
  // ======================================================================
  bool fetchNodes (PHCompositeNode*);                     // guards + cache
  void accumulateFlowContribution(const std::string& calorimeter,
                                    unsigned           ieta,
                                    double             et,
                                    double             phi,
                                    int                ptBin);
    
  void doCaloQA   (PHCompositeNode* topNode,
                   const std::vector<std::string>&);
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
  double      m_towMinE   {0.050};
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

  // --- trigger bookkeeping if using TriggerAnalyzer package -----------------------------------------------
//  std::map<std::string, std::string> triggerNameMap {
//        {"MBD N&S >= 2", "MBD_NandS_geq_2"}};
    
  // --- trigger bookkeeping :  bit‑index  →  human‑readable key ----------
  std::map<int, std::string> triggerNameMap = {
      {10, "MBD_NS_geq_2"},
      {11, "MBD_NS_geq_1"},
      {12, "MBD_NS_geq_2_vtx_lt_10"},
      {13, "MBD_NS_geq_2_vtx_lt_30"},
      {14, "MBD_NS_geq_2_vtx_lt_150"},
      {15, "MBD_NS_geq_1_vtx_lt_10"},
      {16, "photon_6_plus_MBD_NS_geq_2_vtx_lt_10"},
      {17, "photon_8_plus_MBD_NS_geq_2_vtx_lt_10"},
      {18, "photon_10_plus_MBD_NS_geq_2_vtx_lt_10"},
      {19, "photon_12_plus_MBD_NS_geq_2_vtx_lt_10"},
      {20, "photon_6_plus_MBD_NS_geq_2_vtx_lt_150"},
      {21, "photon_8_plus_MBD_NS_geq_2_vtx_lt_150"},
      {22, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"},
      {23, "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"}
  };
    
  std::map<std::string, HistMap>     qaHistogramsByTrigger;
  // ── trigger QA helpers ────────────────────────────────────────────────
  TH2I* h_MBTrigCorr   = nullptr;                 // 2‑D map: MinBias × Trigger
  std::unordered_map<std::string,int> m_trigBin;  // trigger‑key → x‑bin index

  // first‑event gate: MB + trigger selection (declared here, defined in .cc)
  bool firstEventCuts(PHCompositeNode*   topNode,
                        std::vector<std::string>& activeTrig);
  // --- analysis cuts ------------------------------------------------------
  const std::vector<float>               m_asymCuts {0.5f, 0.7f};
  const std::vector<float>               m_chi2Cuts {1.f, 4.f};
  const std::vector<float>               m_minClusE {2.f};
  const std::vector<std::pair<float,float>> m_ptBins {
        {2,4},{4,6},{6,8},{8,10},{10,12},{12,15},{15,20},{20,30} };

  // --- detector lists -----------------------------------------------------
  const std::vector<std::tuple<std::string,std::string,std::string>> m_caloInfo {
        {"TOWERINFO_CALIB_CEMC",   "TOWERGEOM_CEMC",   "CEMC"},
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

  double m_psi1_N = 0.,  m_psi1_S = 0.;    // Ψ1 North / South
  double m_psi3_N = 0.,  m_psi3_S = 0.;;   // Ψ3 North / South
    
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
    
  // --- vn flow helpers ------------------------------------------------------
  static constexpr int kMaxHarm = 3;                       // we need n=2,3
  struct FlowAcc
    {
        double sumW{0.};
        double qx[4]{}, qy[4]{};

        void reset()
        {
          sumW = 0.;
          std::fill(std::begin(qx), std::end(qx), 0.);
          std::fill(std::begin(qy), std::end(qy), 0.);
        }
  };
  std::map<std::string, std::vector<FlowAcc>>  m_flowAcc;  // "CEMC" → 8 pT bins

  void   bookFlowQA(const std::string& trig, HistMap& H);
  void   fillFlowHists(const std::vector<std::string>& trig);
  // --------------------------------------------------------------
  //  Tower-index ⇒ hemisphere lookup (inline, header‑only)
  // --------------------------------------------------------------
  inline bool isSouthCEMC (unsigned ieta) { return ieta < 48; }
  inline bool isNorthCEMC (unsigned ieta) { return ieta >= 48; }

  inline bool isSouthHCal(unsigned ieta) { return ieta < 12; }   // IHCAL & OHCAL
  inline bool isNorthHCal(unsigned ieta) { return ieta >= 12; }
    
  // emcal_sepdCorrelator.h  (private section)
  void  printTriggerSummary(const std::vector<std::string>& active,
                               uint64_t wRaw,uint64_t wLive,uint64_t wScaled) const;
  std::string bitsetToList(uint64_t word) const;

    
    /*
     following two functions are for seperate raw trigger bit QA not using triggerAnalyzer
     */
  inline std::vector<int> extractTriggerBits(uint64_t b_gl1_scaledvec, [[maybe_unused]]int entry) {
        std::vector<int> trig_bits;
        std::bitset<64> bits(b_gl1_scaledvec);
//        if (verbose) {
//            std::cout << "Processing entry " << entry << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;
//        }
//        
        for (unsigned int bit = 0; bit < 64; bit++) {
            if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
                trig_bits.push_back(bit);
            }
        }
        return trig_bits;
  }

  // Inline function to check trigger condition
  inline bool checkTriggerCondition(const std::vector<int> &trig_bits, int inputBit) {
        for (const int &bit : trig_bits) {
            if (bit == inputBit) {
//                if (verbose) {
//                    std::cout << "  Trigger condition met with bit: " << bit << std::endl;
//                }
                
                return true;
            }
        }
//        if (verbose) {
//            std::cout << "  No relevant trigger conditions met." << std::endl;
//        }
        
        return false;
  }

  // --------------------------------------------------------------------------

};
#endif  // EMCALSEPDCORRELATOR_H
