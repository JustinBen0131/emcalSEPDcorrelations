#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <caloreco/CaloTowerStatus.h>

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/DetermineTowerRho.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <jetbackground/TowerRho.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include "/sphenix/user/patsfan753/tutorials/tutorials/CaloDataAnaRun24pp/clusterIsoCopy_src/ClusterIso.h"

#include <calotreegen/caloTreeGen.h>
#include "/sphenix/u/patsfan753/scratch/TriggerAnalysis/src/JetTriggerPlotter.h"
#include <Calo_Calib.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libEMCalSEPD.so)


static const bool WANT_VERBOSE = false;
 
void Fun4All_emcalSEPDcorrelator(const int nEvents = 0,
                         const char *listFile = "input_files.list",
                         const char *inName = "commissioning.root")
{

    Fun4AllServer *se = Fun4AllServer::instance();
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Fun4AllServer instance acquired: " << se << std::endl;
    }
    gSystem->Load("libg4dst");
    
    // Basic run config
    recoConsts *rc = recoConsts::instance();
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Setting CDB_GLOBALTAG to 'ProdA_2024'..." << std::endl;
    }
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    
    // Read the first filename to extract the run number
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Attempting to open listFile: " << listFile << std::endl;
    }
    std::ifstream infile(listFile);
    std::string firstFilename;
    if (!infile.is_open())
    {
        std::cerr << "[ERROR] Could not open input file list: " << listFile << std::endl;
        return;
    }
    if (!std::getline(infile, firstFilename))
    {
        std::cerr << "[ERROR] Input file list is empty: " << listFile << std::endl;
        return;
    }
    std::cout << "[DEBUG] First filename read: " << firstFilename << std::endl;
    
    int runnumber = -1;  // Declare at an outer scope
    // Extract run number from the first filename
    std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstFilename);
    runnumber = runseg.first;
    int segnumber = runseg.second;
    if (WANT_VERBOSE) {
        std::cout << "[DEBUG] Extracted run: " << runnumber
        << " segment: " << segnumber << std::endl;
    }

    if (runnumber <= 0)
    {
        std::cerr << "[ERROR] Invalid run number extracted from first file: "
        << runnumber << ". Exiting..." << std::endl;
        return;
    }
    rc->set_uint64Flag("TIMESTAMP", runnumber);
    
    
    
    auto* correlator = new emcal_sepdCorrelator(inName);    // writes to rootOut
    correlator->setVzCut(30.0);
    correlator->enableVzCut();    // (re)enable – default true
    correlator->setVerbose(false);
    se->registerSubsystem(correlator);
    
    TriggerRunInfoReco *triggerruninforeco = new TriggerRunInfoReco();
    triggerruninforeco->Verbosity(0);
    se->registerSubsystem(triggerruninforeco);

    Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTcalo");


    infile.clear();
    infile.seekg(0, std::ios::beg);

    // Read all filenames and add them to the input manager
    std::string filename;
    while (std::getline(infile, filename))
    {
        if (filename.empty()) continue;
        in->AddFile(filename.c_str());
        if (WANT_VERBOSE) {
            std::cout << "[INFO] Added input file: " << filename << std::endl;
        }

    }
    infile.close();
    se->registerInputManager(in);

    se->run(nEvents);
    se->End();

    // Now done => exit
    if (WANT_VERBOSE)
    {
      std::cout << "[DEBUG] Done => exiting.\n";
    }
    gSystem->Exit(0);
}


#endif

