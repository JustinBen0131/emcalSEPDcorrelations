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
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>         // for PHIODataNode

#include <ffamodules/CDBInterface.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <epd/EpdReco.h>
#include <mbd/MbdReco.h>
#include <globalvertex/GlobalVertexReco.h>

// new – high‑level reconstruction
#include <eventplaneinfo/EventPlaneReco.h>
#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>   // optional but handy

#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>

// analysis module
#include "/sphenix/u/patsfan753/scratch/emcalSEPDcorrelations/src/emcal_sepdCorrelator.h"

// C / C++
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

//–––– ROOT libraries –––––––––––––––––––––––––––––––––––––––––––––––––––
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libeventplaneinfo.so)
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
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
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP",     run);
  PHRandomSeed();

  //--------------------------------------------------------------------
  // 3.  Register reconstruction / analysis subsystems  (⟨strict order⟩)
  //--------------------------------------------------------------------

//    // ------------------------------------------------------------------
//    // erase sEPD container if it exists AND its size != 744  -------------
//    // (older DST_CALO files contain only 721 channels)
//    // ------------------------------------------------------------------
//    auto* top = se->topNode();              // Fun4AllServer already exists
//    TowerInfoContainer* sepd = findNode::getClass<TowerInfoContainer>(
//                                 top, "TOWERINFO_CALIB_SEPD");
//
//    if (sepd && sepd->size() != 744)
//    {
//      // locate the parent composite node ("SEPD") that owns the data node
//      PHNodeIterator it(top);
//      auto* sepdNode = dynamic_cast<PHCompositeNode*>(
//                         it.findFirst("PHCompositeNode", "SEPD"));
//
//      if (sepdNode)
//      {
//        PHNodeIterator it2(sepdNode);
//        if (auto* dataNode = it2.findFirst("PHIODataNode", "TOWERINFO_CALIB_SEPD"))
//        {
//          sepdNode->removeNode(dataNode);
//          std::cout << "[INFO] removed stale sEPD node (size "
//                    << sepd->size() << ") – will rebuild with EpdReco\n";
//        }
//      }
//    }
//    //  -- 3a) (Re)build a correct sEPD container -------------------------
//    auto* epdReco = new EpdReco();      // Verbosity, etc. optional
//    se->registerSubsystem(epdReco);
//
//
//  if (!findNode::getClass<MbdPmtContainer>(top,"MbdPmtContainer"))
//  {
//      se->registerSubsystem( new MbdReco() );
//  }

  // 3b) Primary‑vertex finder (needs MBD & EPD information)
  auto gvr = new GlobalVertexReco();
  se->registerSubsystem( gvr );

  // 3c) Centrality determination (needs vertex + detector charges)
  auto cent = new CentralityReco();
  cent->setOverwriteScale ("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_54280.root");
  cent->setOverwriteVtx   ("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_54280.root");
  cent->setOverwriteDivs  ("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/divs/cdb_centrality_54280.root");
  se->registerSubsystem(cent);

  auto mbclass = new MinimumBiasClassifier();
  mbclass->Verbosity(0);
  se->registerSubsystem(mbclass);

  // 3d) Event‑plane reconstruction (needs calibrated sEPD/MBD)
  auto epreco = new EventPlaneReco();
  epreco->set_sepd_epreco(true);         // build sEPD Q‑vector
  se->registerSubsystem(epreco);

  // 3e) Run‑information helper (optional but handy)
  auto* trigInfo = new TriggerRunInfoReco();
  trigInfo->Verbosity(verbose ? 1 : 0);
  se->registerSubsystem(trigInfo);

  // 3f) User analysis module – must come *last*
  auto* correl = new emcal_sepdCorrelator(outRoot);
  correl->setVzCut(30.);
  correl->enableVzCut(true);
  correl->setVerbose(10);
  se->registerSubsystem(correl);


  //--------------------------------------------------------------------
  // 4.  Input manager
  //--------------------------------------------------------------------
  auto* inDST = new Fun4AllDstInputManager("DSTcalo");
  for (const auto& f : files) inDST->AddFile(f);
  se->registerInputManager(inDST);

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

  //--------------------------------------------------------------------
  // 6.  Clean exit
  //--------------------------------------------------------------------
  gSystem->Exit(0);
}

#endif   // ROOT_VERSION guard
