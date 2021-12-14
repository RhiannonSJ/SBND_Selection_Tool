#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"

using namespace cppsecrets;
using namespace selection;

int MainTest(const char *config){
  
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;

  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the class ConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get variables from config
  std::string input_location  = "";
  std::string input_filename  = "";
  std::string exceptions_file = "";
  std::string stats_location  = "";
  unsigned int total_files = 0;
  unsigned int detector = 0; // 0 = sbnd, 1 = uboone, 2 = icarus
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputFileLocation",input_location);
  p->getValue("InputFileName",    input_filename);
  p->getValue("ExceptionsFile",   exceptions_file);
  p->getValue("StatFileLocation", stats_location); 
  p->getValue("Detector",         detector);
  p->getValue("TotalFiles",       total_files);
  p->getValue("MinXFid",          minx_fid);
  p->getValue("MinYFid",          miny_fid);
  p->getValue("MinZFid",          minz_fid);
  p->getValue("MaxXFid",          maxx_fid);
  p->getValue("MaxYFid",          maxy_fid);
  p->getValue("MaxZFid",          maxz_fid);
  p->getValue("MinXAV",           minx_av);
  p->getValue("MinYAV",           miny_av);
  p->getValue("MinZAV",           minz_av);
  p->getValue("MaxXAV",           maxx_av);
  p->getValue("MaxYAV",           maxy_av);
  p->getValue("MaxZAV",           maxz_av);

  //------------------------------------------------------------------------------------------
  //                                    Initialise
  //------------------------------------------------------------------------------------------

  // Get the active and fiducial geometry objects
  Geometry fiducial(minx_fid,miny_fid,minz_fid,maxx_fid,maxy_fid,maxz_fid,true);
  Geometry active(minx_av,miny_av,minz_av,maxx_av,maxy_av,maxz_av,false);
  PlaneList planes = active.GetExternalPlaneList();

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;

  int start = static_cast<int>(time(NULL));
  double pot = 0.; 

  std::vector<int> exceptions;
  FillExceptions(exceptions_file.c_str(),exceptions);

  // COUNTERS
  unsigned int nue_true  = 0; 

  unsigned int ccinc_cc0pi = 0;
  unsigned int ccinc_cc1pi = 0;
  unsigned int ccinc_ccoth = 0;
  unsigned int ccinc_nc0pi = 0;
  unsigned int ccinc_nc1pi = 0;
  unsigned int ccinc_ncoth = 0;
  unsigned int ccinc_nue   = 0;
  unsigned int ccinc_true  = 0; 
  unsigned int ccinc_sig   = 0; 
  unsigned int ccinc_sel   = 0; 

  unsigned int cc0pi_cc0pi = 0;
  unsigned int cc0pi_cc1pi = 0;
  unsigned int cc0pi_ccoth = 0;
  unsigned int cc0pi_nc0pi = 0;
  unsigned int cc0pi_nc1pi = 0;
  unsigned int cc0pi_ncoth = 0;
  unsigned int cc0pi_nue   = 0;
  unsigned int cc0pi_true  = 0; 
  unsigned int cc0pi_sig   = 0; 
  unsigned int cc0pi_sel   = 0; 

  unsigned int cc1pi_cc0pi = 0;
  unsigned int cc1pi_cc1pi = 0;
  unsigned int cc1pi_ccoth = 0;
  unsigned int cc1pi_nc0pi = 0;
  unsigned int cc1pi_nc1pi = 0;
  unsigned int cc1pi_ncoth = 0;
  unsigned int cc1pi_nue   = 0;
  unsigned int cc1pi_true  = 0; 
  unsigned int cc1pi_sig   = 0; 
  unsigned int cc1pi_sel   = 0; 

  unsigned int ncinc_cc0pi = 0;
  unsigned int ncinc_cc1pi = 0;
  unsigned int ncinc_ccoth = 0;
  unsigned int ncinc_nc0pi = 0;
  unsigned int ncinc_nc1pi = 0;
  unsigned int ncinc_ncoth = 0;
  unsigned int ncinc_nue   = 0;
  unsigned int ncinc_true  = 0; 
  unsigned int ncinc_sig   = 0; 
  unsigned int ncinc_sel   = 0; 

  unsigned int nc0pi_cc0pi = 0;
  unsigned int nc0pi_cc1pi = 0;
  unsigned int nc0pi_ccoth = 0;
  unsigned int nc0pi_nc0pi = 0;
  unsigned int nc0pi_nc1pi = 0;
  unsigned int nc0pi_ncoth = 0;
  unsigned int nc0pi_nue   = 0;
  unsigned int nc0pi_true  = 0; 
  unsigned int nc0pi_sig   = 0; 
  unsigned int nc0pi_sel   = 0; 

  unsigned int nc1pi_cc0pi = 0;
  unsigned int nc1pi_cc1pi = 0;
  unsigned int nc1pi_ccoth = 0;
  unsigned int nc1pi_nc0pi = 0;
  unsigned int nc1pi_nc1pi = 0;
  unsigned int nc1pi_ncoth = 0;
  unsigned int nc1pi_nue   = 0;
  unsigned int nc1pi_true  = 0; 
  unsigned int nc1pi_sig   = 0; 
  unsigned int nc1pi_sel   = 0;

  unsigned int all_tracks_contained   = 0;
  unsigned int precuts_passed         = 0;
  unsigned int ccinc_passed           = 0;
  unsigned int ccinc_topology_passed  = 0;

  TopologyMap nc_map      = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc_map      = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_map   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map   = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc0pi2p_map = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap nc0pi_map   = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_map   = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nue_map     = GeneralAnalysisHelper::GetNuETopologyMap();

  std::vector< TopologyMap > maps({cc_map, cc0pi_map, cc1pi_map, cc0pi2p_map, nc_map, nc0pi_map, nc1pi_map, nue_map});
  std::map< std::string, std::map<std::string, unsigned int> > topology_rates;

  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    for(const Event &e : events){
      // Boolean to check that we have selected a CC Inclusive using the method laid out in PIDMain
      bool cc_inclusive_passed = GeneralAnalysisHelper::PassedCCInclusive(e,detector);

      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial() || 
         !GeneralAnalysisHelper::MaxOneLongEscapingTrack(e) || 
         !GeneralAnalysisHelper::MinOneRecoTrack(e)) continue;
      precuts_passed++;

      // Now we are looking at selected CC Inclusive events, 
      // start trying to identify particle types
      if(cc_inclusive_passed){
        ccinc_passed++;
        if(e.CheckRecoTopology(maps[0])){
          ccinc_topology_passed++;
          if(e.CheckMCTopology(maps[1]))      ccinc_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) ccinc_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) ccinc_ccoth++;
          else if(e.CheckMCTopology(maps[5])) ccinc_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) ccinc_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) ccinc_ncoth++;
          else if(e.CheckMCTopology(maps[7])) ccinc_nue++;
        }
        if(e.CheckRecoTopology(maps[1])){
          if(e.CheckMCTopology(maps[1]))      cc0pi_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) cc0pi_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) cc0pi_ccoth++;
          else if(e.CheckMCTopology(maps[5])) cc0pi_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) cc0pi_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) cc0pi_ncoth++;
          else if(e.CheckMCTopology(maps[7])) cc0pi_nue++;
        }
        if(e.CheckRecoTopology(maps[2])){
          if(e.CheckMCTopology(maps[1])) cc1pi_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) cc1pi_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) cc1pi_ccoth++;
          else if(e.CheckMCTopology(maps[5])) cc1pi_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) cc1pi_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) cc1pi_ncoth++;
          else if(e.CheckMCTopology(maps[7])) cc1pi_nue++;
        }
        if(e.CheckRecoTopology(maps[4])){
          if(e.CheckMCTopology(maps[1]))      ncinc_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) ncinc_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) ncinc_ccoth++;
          else if(e.CheckMCTopology(maps[5])) ncinc_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) ncinc_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) ncinc_ncoth++;
          else if(e.CheckMCTopology(maps[7])) ncinc_nue++;
        }
        if(e.CheckRecoTopology(maps[5])){
          if(e.CheckMCTopology(maps[1]))      nc0pi_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) nc0pi_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) nc0pi_ccoth++;
          else if(e.CheckMCTopology(maps[5])) nc0pi_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) nc0pi_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) nc0pi_ncoth++;
          else if(e.CheckMCTopology(maps[7])) nc0pi_nue++;
        }
        if(e.CheckRecoTopology(maps[6])){
          if(e.CheckMCTopology(maps[1]))      nc1pi_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) nc1pi_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) nc1pi_ccoth++;
          else if(e.CheckMCTopology(maps[5])) nc1pi_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) nc1pi_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) nc1pi_ncoth++;
          else if(e.CheckMCTopology(maps[7])) nc1pi_nue++;
        }
      }
      // Overall efficiencies 
      if(e.CheckMCTopology(maps[0])){
        ccinc_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[0])) ccinc_sig++;
      }
      if(e.CheckMCTopology(maps[1])){
        cc0pi_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[1])) cc0pi_sig++;
      }
      if(e.CheckMCTopology(maps[2])) {
        cc1pi_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[2])) cc1pi_sig++;
      }
      if(e.CheckMCTopology(maps[4])){
        ncinc_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[4])) ncinc_sig++;
      }
      if(e.CheckMCTopology(maps[5])){
        nc0pi_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[5])) nc0pi_sig++;
      }
      if(e.CheckMCTopology(maps[6])){
        nc1pi_true++;
        if(cc_inclusive_passed && e.CheckRecoTopology(maps[6])) nc1pi_sig++;
      }
      if(e.CheckMCTopology(maps[7])) nue_true++;

      // Overall purities
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[0])) ccinc_sel++;
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[1])) cc0pi_sel++;
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[2])) cc1pi_sel++;
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[4])) ncinc_sel++;
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[5])) nc0pi_sel++;
      if(cc_inclusive_passed && e.CheckRecoTopology(maps[6])) nc1pi_sel++;
    }
  }
  std::vector<std::string> reco_topology_names{"CC~Inc.", "CC~0\\pi", "CC~1\\pi", "NC~Inc.", "NC~0\\pi", "NC~1\\pi"};

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"full_topological_breakdown.txt");

  file << "==============================================================================================================" << std::endl;
  //file << " Total number of events with all tracks contained : " << all_tracks_contained << std::endl;
  file << " Total number of events passing the reconstructable cuts  : " << precuts_passed << std::endl;
  file << " Total number of events passing the CC Inclusive cuts     : " << ccinc_passed << std::endl;
  file << " Total number of events selected as CC Inclusive topology : " << ccinc_topology_passed << std::endl;
  file << std::setw(12) << "True \\ Reco" << "||" <<  std::setw(10) << " CC Inc. " << std::setw(10) << " CC 0Pi " << std::setw(10) << " CC 1Pi " << std::setw(10) << " NC Inc. " << std::setw(10) << " NC 0Pi " << std::setw(10) << " NC 1Pi " << std::endl;
  file << std::setw(12) << " CC 0Pi "     << "||" << std::setw(10) << ccinc_cc0pi << std::setw(10) << cc0pi_cc0pi << std::setw(10) << cc1pi_cc0pi << std::setw(10) << ncinc_cc0pi << std::setw(10) << nc0pi_cc0pi << std::setw(10) << nc1pi_cc0pi << std::endl; 
  file << std::setw(12) << " CC 1Pi "     << "||" << std::setw(10) << ccinc_cc1pi << std::setw(10) << cc0pi_cc1pi << std::setw(10) << cc1pi_cc1pi << std::setw(10) << ncinc_cc1pi << std::setw(10) << nc0pi_cc1pi << std::setw(10) << nc1pi_cc1pi <<  std::endl;  
  file << std::setw(12) << " CC Other "   << "||" << std::setw(10) << ccinc_ccoth << std::setw(10) << cc0pi_ccoth << std::setw(10) << cc1pi_ccoth << std::setw(10) << ncinc_ccoth << std::setw(10) << nc0pi_ccoth << std::setw(10) << nc1pi_ccoth << std::endl; 
  file << std::setw(12) << " NC 0Pi "     << "||" << std::setw(10) << ccinc_nc0pi << std::setw(10) << cc0pi_nc0pi << std::setw(10) << cc1pi_nc0pi << std::setw(10) << ncinc_nc0pi << std::setw(10) << nc0pi_nc0pi << std::setw(10) << nc1pi_nc0pi << std::endl;  
  file << std::setw(12) << " NC 1Pi "     << "||" << std::setw(10) << ccinc_nc1pi << std::setw(10) << cc0pi_nc1pi << std::setw(10) << cc1pi_nc1pi << std::setw(10) << ncinc_nc1pi << std::setw(10) << nc0pi_nc1pi << std::setw(10) << nc1pi_nc1pi  << std::endl;  
  file << std::setw(12) << " NC Other "   << "||" << std::setw(10) << ccinc_ncoth << std::setw(10) << cc0pi_ncoth << std::setw(10) << cc1pi_ncoth << std::setw(10) << ncinc_ncoth << std::setw(10) << nc0pi_ncoth << std::setw(10) << nc1pi_ncoth  << std::endl;  
  file << "==============================================================================================================" << std::endl;
  file << " CC Inc.    true       : " << ccinc_true << std::endl; 
  file << " CC 0Pi     true       : " << cc0pi_true << std::endl; 
  file << " CC 1Pi     true       : " << cc1pi_true << std::endl; 
  file << " NC Inc.    true       : " << ncinc_true << std::endl; 
  file << " NC 0Pi     true       : " << nc0pi_true << std::endl; 
  file << " NC 1Pi     true       : " << nc1pi_true << std::endl; 
  file << " Nu E       true       : " << nue_true   << std::endl; 
  file << "===================================================" << std::endl;
  file << " CC Inc.    efficiency : " << ccinc_sig/double(ccinc_true) << std::endl; 
  file << " CC 0Pi     efficiency : " << cc0pi_sig/double(cc0pi_true) << std::endl; 
  file << " CC 1Pi     efficiency : " << cc1pi_sig/double(cc1pi_true) << std::endl; 
  file << " NC Inc.    efficiency : " << ncinc_sig/double(ncinc_true) << std::endl; 
  file << " NC 0Pi     efficiency : " << nc0pi_sig/double(nc0pi_true) << std::endl; 
  file << " NC 1Pi     efficiency : " << nc1pi_sig/double(nc1pi_true) << std::endl; 
  file << "===================================================" << std::endl;
  file << " CC Inc.   purity      : " << ccinc_sig/double(ccinc_sel)  << std::endl; 
  file << " CC 0Pi    purity      : " << cc0pi_sig/double(cc0pi_sel)  << std::endl; 
  file << " CC 1Pi    purity      : " << cc1pi_sig/double(cc1pi_sel)  << std::endl; 
  file << " NC Inc.   purity      : " << ncinc_sig/double(ncinc_sel)  << std::endl; 
  file << " NC 0Pi    purity      : " << nc0pi_sig/double(nc0pi_sel)  << std::endl; 
  file << " NC 1Pi    purity      : " << nc1pi_sig/double(nc1pi_sel)  << std::endl; 
  file << "===================================================" << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

