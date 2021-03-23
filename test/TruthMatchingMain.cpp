#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"

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

  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  EventSelectionTool::EventList all_events;
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    selection::LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    all_events.insert(all_events.end(), events.begin(), events.end());
  }
  std::cout << " Total number of events in sample: " << all_events.size() << std::endl;

  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Files to hold particle statistics
  ofstream all_file;
  all_file.open(stats_location+"particle_stats.txt");

  ofstream mis_id_file;
  mis_id_file.open(stats_location+"mis_identification_stats.txt");

  GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(all_events, all_file, false);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(all_events, cc_signal_map, "CC Inclusive",  all_file, false);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(all_events, cc0pi_signal_map, "CC 0 Pi",    all_file, false);

  GeneralAnalysisHelper::FillGeneralParticleMisIdStatisticsFile(all_events, mis_id_file, false);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(all_events, cc_signal_map, "CC Inclusive",  mis_id_file, false);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(all_events, cc0pi_signal_map, "CC 0 Pi",    mis_id_file, false);

  ofstream all_file_tex;
  all_file_tex.open(stats_location+"particle_stats.tex");

  ofstream mis_id_file_tex;
  mis_id_file_tex.open(stats_location+"mis_identification_stats.tex");

  GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(all_events, all_file_tex, true);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(all_events, cc_signal_map, "CC Inclusive",  all_file_tex, true);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(all_events, cc0pi_signal_map, "CC 0 Pi",    all_file_tex, true);

  GeneralAnalysisHelper::FillGeneralParticleMisIdStatisticsFile(all_events, mis_id_file_tex, true);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(all_events, cc_signal_map, "CC Inclusive",  mis_id_file_tex, true);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(all_events, cc0pi_signal_map, "CC 0 Pi",    mis_id_file_tex, true);

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

