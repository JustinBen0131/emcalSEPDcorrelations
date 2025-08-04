//======================================================================
//  Fun4All_emcalSEPDcorrelator.C  – ana.495‑compatible driver
//  --------------------------------------------------------------------
//  * No dependency on CaloGeomInit or MbdGeomReco (not shipped with ana).
//  * Lots of run‑time sanity checks to pinpoint problems quickly.
//  * Fails hard (std::runtime_error) on any unrecoverable condition.
//======================================================================
#pragma once
#if defined(__CINT__) || defined(__CLING__)
  R__ADD_INCLUDE_PATH($OFFLINE_MAIN/include)
#endif
#if defined(__CLING__)
  #pragma cling add_include_path("$ENV{OFFLINE_MAIN}/include")
#endif

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All / sPHENIX ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <calobase/TowerInfoContainer.h>
#include <mbd/MbdPmtContainer.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>         // for PHIODataNode

#include <ffamodules/FlagHandler.h>
#include <ffamodules/CDBInterface.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <calobase/RawTowerGeom.h>
#include <caloreco/RawTowerCalibration.h>
#include <jetbase/Jet.h>
#include <jetbase/FastJetOptions.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <epd/EpdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerCalib.h>
#include <eventplaneinfo/EventPlaneReco.h>
#include <eventplaneinfo/Eventplaneinfo.h>
#include <centrality/CentralityReco.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>
#include <calotrigger/MinimumBiasClassifier.h>   // optional but handy
#include <zdcinfo/ZdcReco.h>
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <jetbase/FastJetOptions.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/CopyAndSubtractJets.h>

// analysis module
#include "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/src/emcal_sepdCorrelator.h"

// C / C++
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <Calo_Calib.C>

//–––– ROOT libraries –––––––––––––––––––––––––––––––––––––––––––––––––––
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libeventplaneinfo.so)
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY( libzdcinfo.so )
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(/sphenix/user/patsfan753/install/lib/libEMCalSEPD.so)

//======================================================================
//  Convenience helpers
//======================================================================
namespace detail
{
  /// Throw a nicely formatted exception on unrecoverable error
  [[noreturn]] void bail(const std::string& msg)
  {
    std::ostringstream oss;
    oss << "\n[FATAL] Fun4All_emcalSEPDcorrelator :: " << msg << '\n';
    throw std::runtime_error(oss.str());
  }

  /// Trim whitespace from both ends (for robust list‑file parsing)
  inline std::string trim(std::string s)
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  }
}

namespace detail
{
  inline FastJetAlgoSub* fjAlgo(const float R)
  {
    FastJetOptions o({Jet::ANTIKT, Jet::SRC::VOID, R, /*ptmin*/0.0, /*verbosity*/0});
    return new FastJetAlgoSub(o);
  }
}

//======================================================================
//  The actual steering macro
//======================================================================
void Fun4All_emcalSEPDcorrelator(const int   nEvents   =  0,
                                 const char* listFile  = "input_files.list",
                                 const char* outRoot   = "TrigPlot.root",
                                 const bool  verbose   = false)
{
  //--------------------------------------------------------------------
  // 0.  Banner & basic environment sanity
  //--------------------------------------------------------------------
  std::cout << "\n>>> Fun4All_emcalSEPDcorrelator – ana.495 driver <<<\n"
            << "    Input list : " << listFile  << '\n'
            << "    Output file: " << outRoot   << '\n'
            << "    nEvents    : " << nEvents   << (nEvents==0? " (all)\n":"\n");

  Fun4AllServer* se = Fun4AllServer::instance();
  if (!se) detail::bail("unable to obtain Fun4AllServer instance!");

  //--------------------------------------------------------------------
  // 1.  Parse the file list & determine run / segment
  //--------------------------------------------------------------------
  std::ifstream list(listFile);
  if (!list.is_open())
      detail::bail("cannot open input list \"" + std::string(listFile) + "\"");

  std::vector<std::string> files;
  for (std::string line; std::getline(list, line); )
  {
      line = detail::trim(line);
      if (!line.empty()) files.emplace_back(line);
  }
  if (files.empty())
      detail::bail("input list \"" + std::string(listFile) + "\" is empty");

  const std::string& firstFile = files.front();
  const auto [run, seg]        = Fun4AllUtils::GetRunSegment(firstFile);
  if (run <= 0)
      detail::bail("failed to extract run number from first file: " + firstFile);

  if (verbose)
      std::cout << "[INFO] Run=" << run << "  Seg=" << seg
                << "  (" << files.size() << " files)\n";

  //--------------------------------------------------------------------
  // 2.  Global run flags
  //--------------------------------------------------------------------
  recoConsts* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG","newcdbtag");
  rc->set_uint64Flag("TIMESTAMP",     run);
  CDBInterface::instance() -> Verbosity(1);
    
  std::unique_ptr<FlagHandler> flag = std::make_unique<FlagHandler>();
  se->registerSubsystem(flag.get());
    
  auto* inDST = new Fun4AllDstInputManager("DSTcalofitting");
  for (const auto& f : files) inDST->AddFile(f);
  se->registerInputManager(inDST);

    
  for (const std::string& det : {"CEMC","HCALIN","HCALOUT"})
    {
      auto *geom = new CaloGeomMapping(("Geom_"+det).c_str());
      geom->set_detector_name(det);          // << detector tag
      geom->set_UseDetailedGeometry(false);   // (optional, but nice)
      se->registerSubsystem(geom);
  }
//  CaloGeomMapping* geomMap = new CaloGeomMapping("CEMC_GeomFiller");
//  geomMap->set_detector_name("CEMC");
//  geomMap->set_UseDetailedGeometry(true);   // we want the 8-vertex blocks
//  geomMap->Verbosity(0);
//  se->registerSubsystem(geomMap);           // register *before* anything that uses it
//
    
  //////////////////////////////
  // set statuses on raw towers
  std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  se->registerSubsystem(statusEMC);
    
  CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
  statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
  statusHCalIn->set_time_cut(2);
  se->registerSubsystem(statusHCalIn);

  CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
  statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
  statusHCALOUT->set_time_cut(2);
  se->registerSubsystem(statusHCALOUT);


  ////////////////////
  // Calibrate towers
  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  se->registerSubsystem(calibEMC);

  std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  se->registerSubsystem(calibOHCal);

  std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  se->registerSubsystem(calibIHCal);
    
    
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  ClusterBuilder->set_UseAltZVertex(1); // Use MBD Vertex for vertex-based corrections
  se->registerSubsystem(ClusterBuilder);

  std::cout << "Calibrating sEPD" << std::endl;
  std::unique_ptr<EpdReco> epdreco = std::make_unique<EpdReco>();
  epdreco->Verbosity(0);
  se->registerSubsystem(epdreco.get());

  // // MBD/BBC Reconstruction
  std::cout << "Calibrating MBD" << std::endl;
  std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
  se->registerSubsystem(mbdreco.release());

  std::cout << "Calibrating ZDC" << std::endl;
  std::unique_ptr<ZdcReco> zdcreco = std::make_unique<ZdcReco>();
  zdcreco->set_zdc1_cut(0.0);
  zdcreco->set_zdc2_cut(0.0);
  se->registerSubsystem(zdcreco.get());
    
  std::cout << "Retrieving Vtx Info" << std::endl;
  std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
  se->registerSubsystem(gvertex.get());
    
  std::cout << "building minbias classifier" << std::endl;
  auto* mb = new MinimumBiasClassifier();
  mb->Verbosity(0);
  mb->setOverwriteScale(
        "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
  se->registerSubsystem(mb);

    
  std::cout << "building centrality classifier" << std::endl;
  auto* cent = new CentralityReco();
  cent->Verbosity(0);
  cent->setOverwriteScale(
        "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
  se->registerSubsystem(cent);
    
    
  std::cout << "building EP info" << std::endl;
  std::unique_ptr<EventPlaneReco> epreco = std::make_unique<EventPlaneReco>();
  epreco->set_sepd_epreco(true);
  epreco->Verbosity(0);
  se->registerSubsystem(epreco.get());


    
    //--------------------------------------------------------------------
    // 3d)  HI‑style tower‑jet background subtraction (+ jet reco)
    //--------------------------------------------------------------------
    {
      auto banner = [](const std::string &m)
      { std::cout << "\n[BG‑SUB] " << m << std::endl; };

      const int vLvl = 1;            // common verbosity

      /* (i) – one retower is enough for all radii ----------------------- */
      banner("(i)  Retowering EMCal (0.025×0.025) …");
      auto *rcemc = new RetowerCEMC("RetowerCEMC");
      rcemc->set_towerinfo(true);
      rcemc->set_frac_cut(0.50);
      rcemc->set_towerNodePrefix("TOWERINFO_CALIB");
      rcemc->Verbosity(vLvl);
      se->registerSubsystem(rcemc);

      /* (ii‑½) – tower‑by‑tower subtraction, shared by every radius ----- */
      auto *st = new SubtractTowers("SubtractTowers_Sub1");
      st->set_towerinfo(true);
      st->set_towerNodePrefix("TOWERINFO_CALIB");
      st->Verbosity(vLvl);
      se->registerSubsystem(st);

      /* loop over all radii defined in the header ----------------------- */
      for (const auto &[radKey, rawNodeName] :
           emcal_sepdCorrelator::kJetRadii)            // e.g. {"r02","AntiKt_TowerInfo_HIRecoSeedsRaw_r02"}
      {
        const int   D = std::stoi(std::string(radKey).substr(1)); // "r02"→2, "r04"→4 …
        const float R = 0.1f * D;                                 // 0.2, 0.4, …

        /* helper names that carry the radius tag ----------------------- */
        const std::string rawRecoName = "HIRecoSeedsRaw_" + std::string(radKey);
        const std::string subRecoName = "HIRecoSeedsSub_" + std::string(radKey);
        const std::string algoRaw     = "AntiKt_TowerInfo_" + rawRecoName;
        const std::string algoSub = rawNodeName;

        /* (ii) – raw‑seed jets ----------------------------------------- */
        banner("(ii)  Reconstructing raw Anti‑kT jets R=" + std::to_string(R));
        auto *seedReco = new JetReco(rawRecoName);
        seedReco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB"));
        seedReco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB"));
        seedReco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB"));
        seedReco->add_algo(detail::fjAlgo(R), algoRaw);          // **critical name**
        seedReco->set_algo_node("AntiKt_TowerInfo");
        seedReco->set_input_node("TOWERINFO_CALIB");
        seedReco->Verbosity(vLvl);
        se->registerSubsystem(seedReco);

        /* (iii‑a) – UE from raw jets  → TowerInfoBackground_Sub1 -------- */
        auto *dtb1 = new DetermineTowerBackground("DetTowerBkg_Sub1_"+radKey);
        dtb1->SetBackgroundOutputName("TowerInfoBackground_Sub1");
        dtb1->SetSeedType(0);
        dtb1->SetSeedJetD(D);
        dtb1->set_towerinfo(true);
        dtb1->set_towerNodePrefix("TOWERINFO_CALIB");
        dtb1->Verbosity(vLvl);
        se->registerSubsystem(dtb1);

        /* (iii‑b) – UE for tower subtraction  → TowerInfoBackground_Sub2 */
        auto *dtb2 = new DetermineTowerBackground("DetTowerBkg_Sub2_"+radKey);
        dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
        dtb2->SetSeedType(0);
        dtb2->SetSeedJetD(D);
        dtb2->set_towerinfo(true);
        dtb2->set_towerNodePrefix("TOWERINFO_CALIB");
        dtb2->Verbosity(vLvl);
        se->registerSubsystem(dtb2);

        /* (v) – jets on UE‑subtracted towers ---------------------------- */
        banner("(v)   Reconstructing UE‑subtracted jets R=" + std::to_string(R));
        auto *subReco = new JetReco(subRecoName);
        subReco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1,"TOWERINFO_CALIB"));
        subReco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1,"TOWERINFO_CALIB"));
        subReco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1,"TOWERINFO_CALIB"));
        subReco->add_algo(detail::fjAlgo(R), algoSub);
        subReco->set_algo_node("AntiKt_TowerInfo");
        subReco->set_input_node("TOWER");
        subReco->Verbosity(vLvl);
        se->registerSubsystem(subReco);

        /* (vi‑a) – UE from subtracted jets  → TowerInfoBackground_Sub3 -- */
        auto *dtb3 = new DetermineTowerBackground("DetTowerBkg_Sub3_"+radKey);
        dtb3->SetBackgroundOutputName("TowerInfoBackground_Sub3");
        dtb3->SetSeedType(1);
        dtb3->SetSeedJetD(D);
        dtb3->set_towerinfo(true);
        dtb3->set_towerNodePrefix("TOWERINFO_CALIB");
        dtb3->Verbosity(vLvl);
        se->registerSubsystem(dtb3);

        /* (vi‑b) – copy & residual UE subtraction ----------------------- */
        auto *casj = new CopyAndSubtractJets("CopyAndSubtractJets_"+radKey);
        casj->set_towerinfo(true);
        casj->set_towerNodePrefix("TOWERINFO_CALIB");
        casj->Verbosity(3);
        se->registerSubsystem(casj);
      } // radius loop
    }

    
  // 3e) Run‑information helper (optional but handy)
  auto* trigInfo = new TriggerRunInfoReco();
  trigInfo->Verbosity(verbose ? 1 : 0);
  se->registerSubsystem(trigInfo);

  // 3f) User analysis module – must come *last*
  auto* correl = new emcal_sepdCorrelator(outRoot);
  correl->setVzCut(10.);
  correl->enableVzCut(true);
  correl->setVerbose(10);
  se->registerSubsystem(correl);
    

  //--------------------------------------------------------------------
  // 5.  Run
  //--------------------------------------------------------------------
  try
  {
    if (verbose) std::cout << "[INFO] Starting event loop …\n";
    se->run(nEvents);
    se->End();
    if (verbose) std::cout << "[INFO] Finished successfully.\n";
  }
  catch (const std::exception& e)
  {
    detail::bail(std::string("exception in Fun4All: ") + e.what());
  }

//  //--------------------------------------------------------------------
//  // 6.  Clean exit
//  //--------------------------------------------------------------------
//  gSystem->Exit(0);
}

#endif   // ROOT_VERSION guard
