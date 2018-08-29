#include "../include/EventSelectionTool.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/Event.h"
#include <iostream>
#include <iomanip>
#include <fstream>
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
  std::cout << " Start: Local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Running all files " << std::endl;
  
  std::string evd_location = "../Output_Selection_Tool/evd/";
  
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
  
  TopologyMap cc0pi   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc0pi2p = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
 
  // Initialise the file to hold file and event ids for different topologies 
  ofstream file;
  file.open(evd_location+"event_display_ids.txt");
  file << std::endl;
  file << "-------------------------------------------------------------------------------------------------" << std::endl;
  file << std::endl;
  file << std::setw(16) << " Type " << std::setw(16) << "MC, Reco, Both" << std::setw(8) << " File " << std::setw(8) << " Event " << std::endl;
  file << std::endl;
  file << "-------------------------------------------------------------------------------------------------" << std::endl;
  file << std::endl;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.CheckMCTopology(cc0pi2p) && e.IsSBNDTrueFiducial()) {
      file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " MC " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;

      if(e.CheckRecoTopology(cc0pi2p))
        file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " Both " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
    }
    if(e.CheckRecoTopology(cc0pi2p) && e.IsSBNDTrueFiducial())
      file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " Reco " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
  }
 
  file.close();

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
