#include "../include/EventSelectionTool.h"
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
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc0pi_signal_map;
  TopologyMap cc1pi_signal_map;
  TopologyMap ccpi0_signal_map;
 
  std::vector< int > mu;
  std::vector< int > pi;
  std::vector< int > pi0;
  
  mu.push_back( 13 );
  pi.push_back( 211 );
  pi.push_back(-211 );
  pi0.push_back( 111 );
  
  cc0pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc0pi_signal_map.insert( std::make_pair( pi,  0 ) );
  cc0pi_signal_map.insert( std::make_pair( pi0, 0 ) );

  cc1pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc1pi_signal_map.insert( std::make_pair( pi,  1 ) );
  
  ccpi0_signal_map.insert( std::make_pair( mu,  1 ) );
  ccpi0_signal_map.insert( std::make_pair( pi0, 1 ) );
  
  for( unsigned int i = 0; i < 500; ++i ){
  
    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    ss << "/hepstore/rjones/Samples/FNAL/sbn_workshop_0318/5103689_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events, i);
  }

  // Initialise the file to hold file and event ids for different topologies 
  ofstream file;
  file.open("event_display_ids.txt");
  file << std::endl;
  file << "-------------------------------------------------------------------------" << std::endl;
  file << std::endl;
  file << std::setw(8) << " Type " << std::setw(8) << " File " << std::setw(8) << " Event " << std::endl;
  file << std::endl;
  file << "-------------------------------------------------------------------------" << std::endl;
  file << std::endl;

  // Counter for CC 1pi protons
  unsigned int protons_reco = 0;
  unsigned int protons_true = 0;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.CheckMCTopology(cc0pi_signal_map) && e.CheckRecoTopology(cc0pi_signal_map)) {
      file << std::setw(8) << " 0 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    else if(e.CheckRecoTopology(cc0pi_signal_map)) {
      file << std::setw(8) << " 1 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    if(e.CheckMCTopology(cc0pi_signal_map) && !e.CheckRecoTopology(cc0pi_signal_map)) {
      file << std::setw(8) << " 2 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    
    if(e.CheckMCTopology(cc1pi_signal_map) && e.CheckRecoTopology(cc1pi_signal_map)) {
      
      protons_true = e.CountMCParticlesWithPdg(2212);
      protons_reco = e.CountRecoParticlesWithPdg(2212);
     
      if(e.GetPhysicalProcess() == 3) { 
        file << std::setw(8) << " 9 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
      }
      // If we are missing many protons
      if( protons_true >= 4 && protons_reco == 0){
        file << std::setw(8) << " 3p " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
      }
      file << std::setw(8) << " 3 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;

    }
    else if(e.CheckRecoTopology(cc1pi_signal_map)) {
      file << std::setw(8) << " 4 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    if(e.CheckMCTopology(cc1pi_signal_map) && !e.CheckRecoTopology(cc1pi_signal_map)) {
      file << std::setw(8) << " 5 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    
    if(e.CheckMCTopology(ccpi0_signal_map) && e.CheckRecoTopology(ccpi0_signal_map)) {
      file << std::setw(8) << " 6 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
    else if(e.CheckRecoTopology(ccpi0_signal_map)) {
      file << std::setw(8) << " 7 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    } 
    if(e.CheckMCTopology(ccpi0_signal_map) && !e.CheckRecoTopology(ccpi0_signal_map)) {
      file << std::setw(8) << " 8 " << std::setw(8) << e.GetFileNumber() << std::setw(8) << e.GetID() << std::endl;
    }
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
