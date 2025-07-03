// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EMCALSEPDCORRELATOR_H
#define EMCALSEPDCORRELATOR_H

#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include <iomanip>
#include <TTree.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <unordered_map>

//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

class PHCompositeNode;
class Fun4AllHistoManager;
class TFile;
class RawCluster;
class TowerInfoContainer;
class TH1F;
class TH2F;

class emcal_sepdCorrelator : public SubsysReco{
public:
    
    // A constructor that takes two std::string arguments
    emcal_sepdCorrelator(const std::string &dataOutFile = "caloTreeData.root");
    
    ~emcal_sepdCorrelator() override;
    
    int Init(PHCompositeNode *topNode) override;
    
    int process_event(PHCompositeNode *topNode) override;
    
    int ResetEvent(PHCompositeNode *topNode) override;
    
    int End(PHCompositeNode *topNode) override;
    
    int Reset(PHCompositeNode * /*topNode*/) override;
    
    void Print(const std::string &what = "ALL") const override;
    
    void setGenEvent(int eventGet)     {getEvent = eventGet;}

    /// Turn verbose mode on or off
    void setVerbose(bool v) { verbose = v; }
    
    void setRunNumber(int runnumber) { m_runNumber = runnumber; }
    
    /// change the numerical cut (negative values are abs‑ed)
    void setVzCut(double cut)          { m_vzCut   = std::fabs(cut); }
    /// turn the cut on/off from the macro
    void enableVzCut(bool enable=true) { m_useVzCut = enable;        }

private:
    
    int m_runNumber = -1;
    std::map<std::string, double> m_scaleFactors;
    
    TFile *out     = nullptr;  // data

    // 3) Filenames:
    std::string Outfile;     // data output file
    
    int getEvent;
    TriggerAnalyzer* trigAna{nullptr};
    std::size_t event_count = 0;
    
    std::map<std::string, std::string> triggerNameMap = {
        {"MBD N&S >= 1",          "MBD_NandS_geq_1"}
    };
    
    // Pointer to the active trigger name map for the current run
    std::map<int, std::string>* activeTriggerNameMap = nullptr;
    
    static constexpr std::array<std::pair<const char*, const char*>, 2> kJetRadii {{
        {"r03", "AntiKt_unsubtracted_r03"},
        {"r05", "AntiKt_unsubtracted_r05"}
    }};
    
    bool   verbose = true;

    float m_vertex;
    double m_vx, m_vy, m_vz;
    
    double m_vzCut   = 30.0;   // [cm]  default threshold
    bool   m_useVzCut = true;  // enable/disable flag

    std::map<std::string, std::map<std::string, TObject*>> qaHistogramsByTrigger;
    void createHistos_Data();
    
    void checkMbdAndFillNewHists(
            PHCompositeNode*                      topNode,
            const std::unordered_map<std::string,float>& jetEtByRadius);
    
    inline std::string formatFloatForFilename(float value) {
        std::ostringstream ss;
        // Increase the precision to handle more decimal places accurately
        ss << std::fixed << std::setprecision(3) << value;
        std::string str = ss.str();
        size_t dotPos = str.find('.');
        if (dotPos != std::string::npos) {
            // Replace '.' with "point"
            str = str.substr(0, dotPos) + "point" + str.substr(dotPos + 1);
        }
        // Remove trailing zeros and 'point' for whole numbers
        if (value == static_cast<int>(value)) {
            size_t pointPos = str.find("point");
            if (pointPos != std::string::npos) {
                str.erase(pointPos);
            }
        } else {
            // Remove trailing zeros for decimal values
            str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        }
        return str;
    }
};

#endif  

