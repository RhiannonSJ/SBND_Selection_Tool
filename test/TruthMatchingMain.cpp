#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
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

using namespace selection;

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/";
  std::string plots_location  = "../Output_Selection_Tool/plots/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  // Load the events into the event list
  for( unsigned int i = 0; i < 200; ++i ){

    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    ss << "/hepstore/rjones/Samples/FNAL/210518_analysis_sample_200/7700210_" << i << "/output_file.root";
    //ss << "/hepstore/rjones/Samples/FNAL/210518_analysis_sample_200/7726975_" << i << "/output_file.root";
    //ss << "/hepstore/rjones/Samples/FNAL/150518_analysis_sample/7864704_" << i << "/output_file.root";
    //ss << "/hepstore/rjones/Samples/FNAL/sbn_workshop_0318_new/4883618_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
    
    std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;

  }

  // Load the events into the event list
  for( unsigned int i = 0; i < 300; ++i ){

    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    //ss << "/hepstore/rjones/Samples/FNAL/210518_analysis_sample_200/7700210_" << i << "/output_file.root";
    ss << "/hepstore/rjones/Samples/FNAL/210518_analysis_sample_200/7726975_" << i << "/output_file.root";
    //ss << "/hepstore/rjones/Samples/FNAL/150518_analysis_sample/7864704_" << i << "/output_file.root";
    //ss << "/hepstore/rjones/Samples/FNAL/sbn_workshop_0318_new/4883618_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
    
    std::cout << "Loaded file " << std::setw(4) << 200 + i << '\r' << flush;

  }

  std::cout << std::endl;
  
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

  GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(events, all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, nc_signal_map, "NC Inclusive",  all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc_signal_map, "CC Inclusive",  all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc0pi_signal_map, "CC 0 Pi",    all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc1pi_signal_map, "CC 1 Pi+/-", all_file);

  GeneralAnalysisHelper::FillGeneralParticleMisIdStatisticsFile(events, mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, nc_signal_map, "NC Inclusive",  mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc_signal_map, "CC Inclusive",  mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc0pi_signal_map, "CC 0 Pi",    mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc1pi_signal_map, "CC 1 Pi+/-", mis_id_file);

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
