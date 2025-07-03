#include "emcal_sepdCorrelator.h"
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//Fun4All
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <unordered_map>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>

//Tower stuff
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <chrono>  // Include for timing

//GL1 Information
#include <ffarawobjects/Gl1Packet.h>

//for cluster vertex correction
#include <CLHEP/Geometry/Point3D.h>

//for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>
// caloEvalStack for cluster to truth matching
#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <jetbase/Jetv1.h>
#include <jetbase/Jetv2.h>
#include <jetbase/JetContainer.h>

#include <filesystem>
#include <fstream>
#include <locale>

#define ANSI_COLOR_RED_BOLD "\033[1;31m"
#define ANSI_COLOR_BLUE_BOLD "\033[1;34m"
#define ANSI_COLOR_GREEN_BOLD "\033[1;32m"
#define ANSI_COLOR_RESET "\033[0m"


//____________________________________________________________________________
emcal_sepdCorrelator::emcal_sepdCorrelator(const std::string& dataOutFile)
  : SubsysReco("emcal_sepdCorrelator")
  , Outfile(dataOutFile)
{
  // ------------------------------------------------------------------
  // Basic sanity‑check: the caller must supply a non‑empty filename.
  // Abort early if they forgot to pass one.
  // ------------------------------------------------------------------
  if (Outfile.empty())
  {
    std::cerr << "\n[ERROR] emcal_sepdCorrelator constructed with an empty "
                 "output‑file name!\n"
                 "        You must pass the desired .root path, e.g.\n"
                 "        new emcal_sepdCorrelator(\"/path/to/file.root\");\n"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "[DEBUG] emcal_sepdCorrelator::emcal_sepdCorrelator() ctor called\n"
            << "        Data output will go to: " << Outfile << std::endl;
}

//____________________________________________________________________________..
emcal_sepdCorrelator::~emcal_sepdCorrelator() {
    std::cout << "[DEBUG] emcal_sepdCorrelator::~emcal_sepdCorrelator() destructor called." << std::endl;
}

//____________________________________________________________________________..
int emcal_sepdCorrelator::Init(PHCompositeNode* /*topNode*/) {
    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Initializing emcal_sepdCorrelator -- RUNNING Init" << ANSI_COLOR_RESET << std::endl;
    }
    
    if (verbose) {
        std::cout << "[INFO] Running in DATA mode." << std::endl;
    }
    out = new TFile(Outfile.c_str(),"RECREATE");

    trigAna = new TriggerAnalyzer();

    createHistos_Data();

    //so that the histos actually get written out
    Fun4AllServer *se = Fun4AllServer::instance();
    if (verbose) {
        se -> Print("NODETREE");
    }
    std::cout << "TriggerAnalyzer::Init(PHCompositeNode *topNode) Initializing" << std::endl;


    // If we get here, we have at least one of them set to true
    return Fun4AllReturnCodes::EVENT_OK;
}




//____________________________________________________________________________..
void emcal_sepdCorrelator::createHistos_Data() {
    std::cout << "[DEBUG] Entering emcal_sepdCorrelator::createHistos()..." << std::endl;

    for (const auto& kv : triggerNameMap)
    {
        // kv.first  -> DB name  (unused here except for logging)
        // kv.second -> shortName
        const std::string& triggerName = kv.second;
        
        if (verbose)
        {
            std::cout << ANSI_COLOR_BLUE_BOLD
            << "Creating histograms for trigger: " << triggerName
            << ANSI_COLOR_RESET << std::endl;
        }
        
        // Create a directory for the current trigger
        TDirectory* triggerDir = out->mkdir(triggerName.c_str());
        if (!triggerDir) {
            std::cerr << "[ERROR] Failed to create directory for trigger: " << triggerName << std::endl;
            exit(EXIT_FAILURE);
        }
        triggerDir->cd(); // Set the current directory to the trigger directory
        
        
        std::map<std::string, TObject*>& qaHistograms = qaHistogramsByTrigger[triggerName];
        // Helper functions for creating histograms with logging
        auto createHistogram = [&](const std::string& name, const std::string& title, int bins, double xMin, double xMax) {
            // Check if a histogram with the same name already exists
            if (out->Get(name.c_str()) != nullptr) {
                std::cerr << "\n[ERROR] Duplicate histogram detected: " << name << std::endl;
                std::cerr << "A histogram with this name already exists in the output file." << std::endl;
                std::cerr << "Aborting Fun4All macro to avoid further conflicts." << std::endl;
                exit(EXIT_FAILURE);
            }
            
            if (verbose) {
                std::cout << ANSI_COLOR_RED_BOLD << "Creating histogram: " << name << " - " << title << ANSI_COLOR_RESET << std::endl;
            }
            
            TH1F* hist = new TH1F(name.c_str(), title.c_str(), bins, xMin, xMax);
            hist->SetDirectory(out); // Ensure it is linked to the output file
            return hist;
        };
        

        for (const auto& r : kJetRadii)
        {
            const std::string tag = std::string("_") + r.first;           // "_r03", "_r06"

            qaHistograms["h_leadingJetET" + tag + "_" + triggerName] =
                createHistogram("h_leadingJetET" + tag + "_" + triggerName,
                                "Leading Jet E_{T} (" + std::string(r.first) + "); Jet E_{T} [GeV]",
                                50, 0, 50);

            qaHistograms["h_leadingJetET" + tag + "_NewTriggerFilling_doNotScale_" + triggerName] =
                createHistogram("h_leadingJetET" + tag + "_NewTriggerFilling_doNotScale_" + triggerName,
                                "Leading Jet E_{T} (" + std::string(r.first) + "); Jet E_{T} [GeV]",
                                50, 0, 50);
        }
    }
}


namespace {
    inline float getMaxJetEt(JetContainer* jets)
    {
        if (!jets) return 0.f;
        float m = 0.f;
        for (const auto* j : *jets) m = std::max(m, j->get_et());
        return m;
    }
}


//____________________________________________________________________________..
int emcal_sepdCorrelator::process_event(PHCompositeNode *topNode)
{
    event_count++;

    std::cout << "\n========== Processing emcal_sepdCorrelator -- Event " << event_count << " ==========\n";

    // 1) Decode triggers
    if (trigAna)
    {
        if (verbose)
        {
            std::cout << "[DEBUG] About to call trigAna->decodeTriggers()"
                      << " on topNode=" << topNode << std::endl;
        }
        trigAna->decodeTriggers(topNode);
    }
    else
    {
        std::cerr << "[ERROR] No emcal_sepdCorrelator pointer!\n";
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    // 2) Build a vector of fired trigger *short* names using the DB-names in triggerNameMap.
    std::vector<std::string> activeTriggerNames;
    activeTriggerNames.reserve(triggerNameMap.size());

    for (const auto &kv : triggerNameMap)
    {
        const std::string &dbTriggerName   = kv.first;   // how it's known in the DB
        const std::string &histFriendlyStr = kv.second;  // short name for histograms

        // 3) Check if this DB trigger fired
        if (trigAna->didTriggerFire(dbTriggerName))
        {
            // 4) Save the *short* name for future reference
            activeTriggerNames.push_back(histFriendlyStr);

            if (verbose)
            {
                std::cout << "Trigger fired: \"" << dbTriggerName
                          << "\" => short name \"" << histFriendlyStr << "\"" << std::endl;
            }
        }
    }
    
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    m_vx = m_vy = m_vz = 0; // Initialize vertex coordinates to zero

    // Check if the GlobalVertexMap node is missing
    if (!vertexmap) {
        std::cout << "Error: GlobalVertexMap node is missing." << std::endl;
    } else {
        if (verbose) {
            std::cout << "GlobalVertexMap node found." << std::endl;
        }
        // Check if the vertex map is empty
        if (vertexmap->empty()) {
            if (verbose) {
                std::cout << "Warning: GlobalVertexMap is empty." << std::endl;
            }
            return Fun4AllReturnCodes::ABORTEVENT;
        }

        // Access the first vertex in the map
        GlobalVertex* vtx = vertexmap->begin()->second;

        if (vtx) {
            // Retrieve vertex coordinates
            m_vx = vtx->get_x();
            m_vy = vtx->get_y();
            m_vz = vtx->get_z();
            if (verbose) {
                std::cout << "Vertex coordinates retrieved: "
                          << "x = " << m_vx << ", "
                          << "y = " << m_vy << ", "
                          << "z = " << m_vz << std::endl;
            }
            if (m_useVzCut && std::abs(m_vz) >= m_vzCut)
            {
                if (verbose)
                {
                    std::cout << "Skipping event: |m_vz| = "
                              << std::abs(m_vz) << " ≥ " << m_vzCut << " cm (configured)\n";
                }
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            else if (verbose)
            {
                std::cout << "Vertex within cut.\n";
            }

        } else if (verbose) {
            std::cout << "Warning: Vertex object is null." << std::endl;
        }
    }

    // -----------------------------------------------------------------------
    // Build a map: radius‑tag  ->  max jet Et in this event
    std::unordered_map<std::string, float> maxJetEt;
    for (const auto& r : kJetRadii)
    {
        auto* jets = findNode::getClass<JetContainer>(topNode, r.second);
        if (!jets)
        {
            if (verbose)
                std::cout << "Aborting run: missing jet container " << r.second << " …\n";
            return Fun4AllReturnCodes::ABORTRUN;
        }
        maxJetEt[r.first] = getMaxJetEt(jets);
    }
    // ---------------------------------------------------------------------------
    // Loop over active trigger names  (UNCHANGED)
    for (const std::string& firedShortName : activeTriggerNames)
    {
        auto& qaHistograms = qaHistogramsByTrigger[firedShortName];

        // Fill the scaled (“normal”) histograms – one per radius
        for (const auto& kv : maxJetEt)               // kv.first = "r03"/"r06"
        {
            const std::string histName =
                "h_leadingJetET_" + kv.first + "_" + firedShortName;

            TH1F* h = static_cast<TH1F*>(qaHistograms[histName]);
            if (!h)
            {
                std::cerr << "Error: Histogram " << histName << " not found.\n";
                continue;
            }
            h->Fill(kv.second);

            if (verbose)
                std::cout << "Filled " << histName << " with " << kv.second << '\n';
        }
    }
    // ----------  NEW: fill the “doNotScale” histograms -------------------------
    checkMbdAndFillNewHists(topNode, maxJetEt);
    // ---------------------------------------------------------------------------


    return Fun4AllReturnCodes::EVENT_OK;
}



//____________________________________________________________________________
void emcal_sepdCorrelator::checkMbdAndFillNewHists(
    PHCompositeNode*                            topNode,
    const std::unordered_map<std::string,float>& jetEtByRadius)
{
    if (!trigAna)
    {
        std::cerr << "[ERROR] No emcal_sepdCorrelator pointer!\n";
        return;
    }

    // Make sure trigger bits for *this* event are decoded
    trigAna->decodeTriggers(topNode);

    // --------------------------------------------------------------------
    // We always demand the *RAW* MBD bit; if it is absent we do nothing.
    // --------------------------------------------------------------------
    const std::string mbdDbName    = "MBD N&S >= 1";
    const std::string mbdShortName = "MBD_NandS_geq_1";

    if (!trigAna->checkRawTrigger(mbdDbName))
    {
        if (verbose)
            std::cout << "[INFO] Raw MBD bit did not fire – skipping do‑not‑scale hists.\n";
        return;
    }

    // --------------------------------------------------------------------
    // 1) Fill MBD’s own “do not scale” histograms (one per radius)
    // --------------------------------------------------------------------
    for (const auto& r : kJetRadii)   // r.first = "r03"/"r06"
    {
        const std::string tag   = "_" + std::string(r.first);   // "_r03"
        const std::string hName = "h_leadingJetET" + tag +
                                  "_NewTriggerFilling_doNotScale_" +
                                  mbdShortName;

        auto& histMap = qaHistogramsByTrigger[mbdShortName];
        auto  it      = histMap.find(hName);
        if (it == histMap.end() || !(it->second))
        {
            if (verbose)
                std::cerr << "[WARNING] Histogram " << hName << " not booked.\n";
            continue;
        }

        if (auto* h = dynamic_cast<TH1F*>(it->second))
        {
            h->Fill(jetEtByRadius.at(r.first));
            if (verbose)
                std::cout << "[INFO] Filled " << hName << " with "
                        << jetEtByRadius.at(r.first) << '\n';
        }
    }

    // --------------------------------------------------------------------
    // 2) For every *other* trigger we now require:
    //       (a) raw‑MBD bit  **and**
    //       (b) *scaled* version of the rare trigger (didTriggerFire)
    // --------------------------------------------------------------------
    for (const auto& kv : triggerNameMap)
    {
        const std::string& dbTriggerName   = kv.first;   // DB name
        const std::string& histFriendlyStr = kv.second;  // short name

        if (dbTriggerName == mbdDbName) continue;        // already handled

        // need the scaled (live) bit for the rare trigger
        if (!trigAna->didTriggerFire(dbTriggerName)) continue;

        for (const auto& r : kJetRadii)
        {
            const std::string tag   = "_" + std::string(r.first);
            const std::string hName = "h_leadingJetET" + tag +
                                      "_NewTriggerFilling_doNotScale_" +
                                      histFriendlyStr;

            auto& histMap = qaHistogramsByTrigger[histFriendlyStr];
            auto  it      = histMap.find(hName);
            if (it == histMap.end() || !(it->second))
            {
                if (verbose)
                    std::cerr << "[WARNING] Histogram " << hName << " not booked.\n";
                continue;
            }

            if (auto* h = dynamic_cast<TH1F*>(it->second))
            {
                h->Fill(jetEtByRadius.at(r.first));
            }
        }
    }
}

    


int emcal_sepdCorrelator::ResetEvent(PHCompositeNode* /*topNode*/)
{

    if (verbose) {
        std::cout << ANSI_COLOR_BLUE_BOLD << "Resetting event..." << ANSI_COLOR_RESET << std::endl;
    }
    
    return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________
int emcal_sepdCorrelator::End(PHCompositeNode* /*topNode*/)
{
    if (verbose)
    {
        std::cout << ANSI_COLOR_BLUE_BOLD
                  << "emcal_sepdCorrelator::End() – finishing up and writing histograms"
                  << ANSI_COLOR_RESET << std::endl;
    }

    // ------------------------------------------------------------------ //
    // 1) Sanity‑check the output file
    // ------------------------------------------------------------------ //
    if (!out || !out->IsOpen())
    {
        std::cerr << ANSI_COLOR_RED_BOLD
                  << "[ERROR] Output file is not open – nothing to write!"
                  << ANSI_COLOR_RESET << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }
    out->cd();

    // ------------------------------------------------------------------ //
    // 2) Write the QA histograms you actually filled
    // ------------------------------------------------------------------ //
    std::size_t totalWritten = 0;

    for (auto& trigKV : qaHistogramsByTrigger)
    {
        const std::string& trigName   = trigKV.first;   // e.g. "Jet_8_GeV_plus_MBD_NS_geq_1"
        auto&              histMap    = trigKV.second;  // std::map<std::string,TObject*>

        // Make sure the directory exists
        TDirectory* dir = out->GetDirectory(trigName.c_str());
        if (!dir) dir = out->mkdir(trigName.c_str());
        dir->cd();

        for (auto& kv : histMap)
        {
            const std::string& hName = kv.first;
            TObject*           obj   = kv.second;

            if (!obj) continue;                            // null pointer guard
            TH1* h = dynamic_cast<TH1*>(obj);
            if (!h)        continue;                       // not a TH1 – ignore
            if (h->GetEntries() == 0) continue;            // skip empty hists

            if (h->Write() == 0)
            {
                std::cerr << ANSI_COLOR_RED_BOLD
                          << "[ERROR] Failed to write histogram " << hName
                          << ANSI_COLOR_RESET << std::endl;
            }
            else if (verbose)
            {
                std::cout << "[INFO] Wrote " << hName
                          << " (" << h->GetEntries() << " entries)"
                          << " to directory " << trigName << '\n';
            }
            ++totalWritten;
        }
    }

    if (verbose)
    {
        std::cout << "[INFO] Wrote a total of " << totalWritten
                  << " non‑empty histograms.\n";
    }

    // ------------------------------------------------------------------ //
    // 3) Close and clean up
    // ------------------------------------------------------------------ //
    out->cd();
    out->Write();     // flush TFile metadata (directory structure, etc.)
    out->Close();
    delete out;
    out = nullptr;

    if (verbose)
    {
        std::cout << ANSI_COLOR_GREEN_BOLD
                  << "emcal_sepdCorrelator::End() completed successfully."
                  << ANSI_COLOR_RESET << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int emcal_sepdCorrelator::Reset(PHCompositeNode* /*topNode*/) {
 std::cout << "emcal_sepdCorrelator::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void emcal_sepdCorrelator::Print(const std::string &what) const {
  std::cout << "emcal_sepdCorrelator::Print(const std::string &what) const Printing info for " << what << std::endl;
}

