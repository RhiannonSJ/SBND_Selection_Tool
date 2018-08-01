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
  std::string stats_location = "../Output_Selection_Tool/statistics/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap numu_signal_map  = GeneralAnalysisHelper::GetNuMuTopologyMap();
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){

    if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    ss << "/hepstore/rjones/Samples/FNAL/old_220518_ana_files/8110339_" << i << "/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
    
    //std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
 
  /*
   *
   *      TESTING
   *
   */
  // Counter for total number of escaping tracks
  unsigned int too_many_escape         = 0;
  unsigned int too_many_true_contained = 0;
  unsigned int too_many_true_escaping  = 0;

  for(const Event &e : events){
    //Counter for event-based track counting
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
      too_many_escape++;
      if(e.IsSBNDTrueFiducial()) too_many_true_contained++;
      else too_many_true_escaping++;
    }
  
  }

  /*
   *
   *      OUTPUTS
   *
   */

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"containment_testing.txt");


  file << "=====================================================================" << std::endl;
  file << " Total number of events                                       : " << events.size()   << std::endl;
  file << " Number of events with more than one escaping track           : " << too_many_escape << std::endl;
  file << " Events with more than one escaping and true vertex contained : " << too_many_true_contained << std::endl;
  file << " Events with more than one escaping and true vertex escaping  : " << too_many_true_escaping  << std::endl;
  file << "=====================================================================" << std::endl;

  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
