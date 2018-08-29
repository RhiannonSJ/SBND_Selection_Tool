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
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

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
  std::string stats_location = "../Output_Selection_Tool/statistics/breakdown/track_pid_testing/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){

    //if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/hepstore/rjones/Samples/FNAL/290818_analysis_sample/11206561_"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/old_220518_ana_files/8110339_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    
    //std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
  }
  std::cout << std::endl;

  // COUNTERS
  unsigned int ccinc_cc0pi = 0;
  unsigned int ccinc_cc1pi = 0;
  unsigned int ccinc_ccoth = 0;
  unsigned int ccinc_ncoth = 0;
  unsigned int ccinc_true  = 0; 
  unsigned int ccinc_sig   = 0; 
  unsigned int ccinc_sel   = 0; 
  
  unsigned int cc0pi_cc0pi = 0;
  unsigned int cc0pi_cc1pi = 0;
  unsigned int cc0pi_ccoth = 0;
  unsigned int cc0pi_ncoth = 0;
  unsigned int cc0pi_true  = 0; 
  unsigned int cc0pi_sig   = 0; 
  unsigned int cc0pi_sel   = 0; 
  
  unsigned int cc1pi_cc0pi = 0;
  unsigned int cc1pi_cc1pi = 0;
  unsigned int cc1pi_ccoth = 0;
  unsigned int cc1pi_ncoth = 0;
  unsigned int cc1pi_true  = 0; 
  unsigned int cc1pi_sig   = 0; 
  unsigned int cc1pi_sel   = 0; 

  unsigned int all_tracks_contained   = 0;
  unsigned int max_one_escaping_track = 0;

  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();

  std::vector< TopologyMap > maps({cc_signal_map, cc0pi_signal_map, cc1pi_signal_map});  

  // First, ensure all tracks are contained
  for(const Event &e : events){

    //if(!e.AllContained()) continue;
    //all_tracks_contained++;

    // Check the true vertex is in the fiducial volume
    if(e.IsSBNDTrueFiducial()){
      if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;
      max_one_escaping_track++;
      
      if(e.CheckRecoTopology(maps[0])){
        if(e.CheckMCTopology(maps[1]))     ccinc_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) ccinc_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) ccinc_ccoth++;
        else ccinc_ncoth++;
      }
      if(e.CheckRecoTopology(maps[1])){
        if(e.CheckMCTopology(maps[1]))     cc0pi_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) cc0pi_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) cc0pi_ccoth++;
        else cc0pi_ncoth++;
      }
      if(e.CheckRecoTopology(maps[2])){
        if(e.CheckMCTopology(maps[1]))     cc1pi_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) cc1pi_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) cc1pi_ccoth++;
        else cc1pi_ncoth++;
      }

      // Overall efficiencies 
      if(e.CheckMCTopology(maps[0])){
        ccinc_true++;
        if(e.CheckRecoTopology(maps[0])) ccinc_sig++;
      }
      if(e.CheckMCTopology(maps[1])){
        cc0pi_true++;
        if(e.CheckRecoTopology(maps[1])) cc0pi_sig++;
      }
      if(e.CheckMCTopology(maps[2])){
        cc1pi_true++;
        if(e.CheckRecoTopology(maps[2])) cc1pi_sig++;
      }
      // Overall purities
      if(e.CheckRecoTopology(maps[0])) ccinc_sel++;
      if(e.CheckRecoTopology(maps[1])) cc0pi_sel++;
      if(e.CheckRecoTopology(maps[2])) cc1pi_sel++;
    }
  }

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"raquel_topology_breakdown.txt");

  file << "================================================================" << std::endl;
  //file << " Total number of events with all tracks contained : " << all_tracks_contained << std::endl;
  file << " Total number of events with maximum one escaping track : " << max_one_escaping_track << std::endl;
  file << std::setw(12) << "True \\ Reco" << "||" <<  std::setw(10) << " CC Inc. " << std::setw(10) << " CC 0Pi " << std::setw(10) << " CC 1Pi " << std::endl;
  file << std::setw(12) << " CC 0Pi "     << "||" << std::setw(10) << ccinc_cc0pi << std::setw(10) << cc0pi_cc0pi << std::setw(10) << cc1pi_cc0pi << std::endl;  
  file << std::setw(12) << " CC 1Pi "     << "||" << std::setw(10) << ccinc_cc1pi << std::setw(10) << cc0pi_cc1pi << std::setw(10) << cc1pi_cc1pi << std::endl;  
  file << std::setw(12) << " CC Other "   << "||" << std::setw(10) << ccinc_ccoth << std::setw(10) << cc0pi_ccoth << std::setw(10) << cc1pi_ccoth << std::endl;  
  file << std::setw(12) << " NC "         << "||" << std::setw(10) << ccinc_ncoth << std::setw(10) << cc0pi_ncoth << std::setw(10) << cc1pi_ncoth << std::endl;  
  file << "================================================================" << std::endl;
  file << " CC Inc. efficiency : " << ccinc_sig/double(ccinc_true) << std::endl; 
  file << " CC 0Pi  efficiency : " << cc0pi_sig/double(cc0pi_true) << std::endl; 
  file << " CC 1Pi  efficiency : " << cc1pi_sig/double(cc1pi_true) << std::endl; 
  file << "================================================================" << std::endl;
  file << " CC Inc. purity     : " << ccinc_sig/double(ccinc_sel)  << std::endl; 
  file << " CC 0Pi  purity     : " << cc0pi_sig/double(cc0pi_sel)  << std::endl; 
  file << " CC 1Pi  purity     : " << cc1pi_sig/double(cc1pi_sel)  << std::endl; 
  file << "================================================================" << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
