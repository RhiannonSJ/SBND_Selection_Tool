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
  std::string tables_location = "../Output_Selection_Tool/tables/";
  std::string plots_location  = "../Output_Selection_Tool/plots/";

  //------------------------------------------------------------------------------------------
  //                                        Counters
  //------------------------------------------------------------------------------------------

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
  for( unsigned int i = 0; i < 500; ++i ){
  
    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    ss << "/hepstore/rjones/Samples/FNAL/sbn_workshop_0318_new/4883618_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
  }
  
  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc_signal_map) && e.CheckRecoTopology(cc_signal_map)) correctly_reconstructed_cc++;
    if(e.CheckRecoTopology(cc_signal_map)) reco_topology_cc++;
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc_signal_map))   true_topology_cc++;

  }
  
  std::cout << "===========================================================" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "===========================================================" << std::endl;

  /*
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.58, 0.68, 0.88, 0.88 );
  l->SetBorderSize(0);

  l->AddEntry( h_true_energy,      " True Monte Carlo",    "l" );
  l->AddEntry( h_reco_energy,      " Selected ",           "l" );
  l->AddEntry( h_good_energy,      " Correctly selected ", "l" );
  */

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
