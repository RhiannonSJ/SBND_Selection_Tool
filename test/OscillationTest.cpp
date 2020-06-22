#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
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
#include "TROOT.h"
#include "TAxis.h"

using namespace selection;

void LoadAllEvents(EventSelectionTool::EventList &events, 
                   const unsigned int &start_file, 
                   const unsigned int &end_file, 
                   const int &start_time, 
                   double &pot, 
                   std::vector<unsigned int> &exceptions);

int MainTest(const unsigned int &start_file = 0, const unsigned int &end_file = 1){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string file_location  = "/sbnd/data/users/rsjones/Output_Selection_Tool/files/valor/cut75cm/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  // Maps
  TopologyMap ccinc = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc2pi = GeneralAnalysisHelper::GetCC2PiTopologyMap();
  TopologyMap ccpi0 = GeneralAnalysisHelper::GetCCPi0TopologyMap();

  // Variables to get
  bool iscc, isnc;  // Charged or neutral current event current
  int nu_pdg, mode; // Neutrino pdg code and scattering mode
  double baseline; // Baseline of the neutrino from the flux
  double enu_true, enu_reco, enu_reco_mc, qsqr; // Neutrino kinematics, truth and reco
  double mu_momentum, mu_cos_z; // Muon kinematics, reco
  int nkaons, npip, npim, npi0; // Particle counting, truth
  // Enumeration (all CC):
  //   -1 = undefined
  //    0 = 0pi
  //    1 = 1chpi
  //    2 = 2chpi
  //    3 = 1pi0
  //    4 = other
  int true_topology = -1;
  int reco_topology = -1;
  
  double pot; // Subrun information

  TTree *t_run    = new TTree("valor_tree"," Tree to hold variables needed for the VALOR analysis");
  TTree *t_subrun = new TTree("subrun_tree"," Tree to hold subrun variables, such as POT");
  t_run->SetDirectory(0);
  t_subrun->SetDirectory(0);
 
  TH1D *h_ccinc_e_true = new TH1D("h_ccinc_e_true","Neutrino energy distribution for a true CC Inclusive final state",100, 0, 3);
  TH1D *h_ccinc_e_sig  = new TH1D("h_ccinc_e_sig", "Neutrino energy distribution for a signal CC Inclusive final state",100, 0, 3);

  t_run->Branch("iscc",          &iscc,          "iscc/O");
  t_run->Branch("isnc",          &isnc,          "isnc/O");
  t_run->Branch("baseline",      &baseline,      "baseline/D");
  t_run->Branch("nu_pdg",        &nu_pdg,        "nu_pdg/I");
  t_run->Branch("enu_true",      &enu_true,      "enu_true/D");
  t_run->Branch("enu_reco_mc",   &enu_reco_mc,   "enu_reco_mc/D");
  t_run->Branch("enu_reco",      &enu_reco,      "enu_reco/D");
  t_run->Branch("mu_momentum",   &mu_momentum,   "mu_momentum/D");
  t_run->Branch("mu_cos_z",      &mu_cos_z,      "mu_cos_z/D");
  t_run->Branch("qsqr",          &qsqr,          "qsqr/D");
  t_run->Branch("mode",          &mode,          "mode/I");
  t_run->Branch("nkaons",        &nkaons,        "nkaons/I");
  t_run->Branch("npip",          &npip,          "npip/I");
  t_run->Branch("npim",          &npim,          "npim/I");
  t_run->Branch("npi0",          &npi0,          "npi0/I");
  t_run->Branch("true_topology", &true_topology, "true_topology/I");
  t_run->Branch("reco_topology", &reco_topology, "reco_topology/I");

  t_subrun->Branch("pot",        &pot,           "pot/D");
  
  int start = static_cast<int>(time(NULL));
  std::vector<unsigned int> exceptions;
  exceptions.clear();

  // Read in txt file of list of empty input directories
  std::fstream exception_file("exceptions.txt");
  if(!exception_file)
    std::cout << " No exceptions file given" << std::endl;
  else{
    std::string s_exc;
    while (std::getline(exception_file, s_exc)) {
      unsigned int i_exc;
      std::istringstream ss_exc(s_exc);
      ss_exc >> i_exc;
      exceptions.push_back(i_exc);
      ss_exc.str(std::string());
      s_exc.clear();
    }
  }

  /*
  std::cout << " Skipping files in directory : " << std::endl;
  for(const int & ex : exceptions)
    std::cout << " - " << ex << " - ";
  std::cout << std::endl;
  */

  LoadAllEvents(events, start_file, end_file, start, pot, exceptions);

  // Counter to quantify how many events have strange particle energies in them
  int bad_events = 0;

  // Loop over events and perform vertexing study
  for(const Event &e : events){
    bool bad_event = false;

    // TPC track criteria
    if(!e.IsSBNDRecoFiducial()) continue;
    if(!GeneralAnalysisHelper::MaxOneLongEscapingTrack(e)) continue;
    iscc     = e.GetIsCC();
    isnc     = !e.GetIsCC();
    nu_pdg   = e.GetNeutrinoPdgCode();
    enu_true = e.GetTrueNuEnergy();
    qsqr     = e.GetTrueNuQ2();
    mode     = e.GetPhysicalProcess();
    baseline = e.GetBaseline();

    // Start the counters
    nkaons      = 0;
    npip        = 0;
    npim        = 0;
    npi0        = 0;
    enu_reco    = 0.;
    enu_reco_mc = 0.;
    
    // Particles 
    ParticleList mc   = e.GetMCParticleList();
    ParticleList reco = e.GetRecoParticleList();

    // Neutrino vertex is within the ficudial border
    if(e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(ccinc)){
        if(!bad_event){
          for(const Particle &p : reco){
            if(p.GetKineticEnergy() <= 0.) continue;
            if(p.GetKineticEnergy() > 10) { // Higher than 10 GeV
              bad_event = true;
              bad_events++;
              std::cerr << " Badly defined energy of the particle (>10GeV), skipping event " << std::endl;
              break;
            }
            if(p.GetPdgCode() == 13){
              mu_momentum = p.GetModulusMomentum();
              mu_cos_z    = p.GetCosTheta();
            } 
            double mass = 0.;
            if(p.GetFromRecoTrack() && GeneralAnalysisHelper::ParticleHasAMatch(e,p) >= 0) {
              if(GeneralAnalysisHelper::GetBestMCParticle(e,p).GetPdgCode() == 211  || 
                  GeneralAnalysisHelper::GetBestMCParticle(e,p).GetPdgCode() == -211 ||
                  GeneralAnalysisHelper::GetBestMCParticle(e,p).GetPdgCode() == 13){
                mass = GeneralAnalysisHelper::GetBestMCParticle(e,p).GetMass();
              }
            }
            // Calculate the reconstructed energy from the reco kinetic (visible) + mass of outgoing particles 
            // Do not include hadron mass since the interaction didn't *produce* a hadron, it was *with* a hadron
            enu_reco += (p.GetKineticEnergy() + mass);
          }
          for(const Particle &p : mc){
            if(p.GetPdgCode() ==  311 || p.GetPdgCode() == -321 || p.GetPdgCode() == 321) nkaons++;
            if(p.GetPdgCode() ==  211) npip++;
            if(p.GetPdgCode() == -211) npim++;
            if(p.GetPdgCode() ==  111) npi0++; 
            
            // MC reconstructed neutrino energy 
            double mass = 0.;
            if(p.GetPdgCode() == 211  || 
               p.GetPdgCode() == -211 ||
               p.GetPdgCode() == 13){
              mass = p.GetMass();
            } // Pdgcodes
            enu_reco_mc += (p.GetKineticEnergy() + mass);
          } // Particles

          if(e.CheckRecoTopology(cc0pi))
            reco_topology = 0;
          else if(e.CheckRecoTopology(cc1pi))
            reco_topology = 1;
          else if(e.CheckRecoTopology(cc2pi))
            reco_topology = 2;
          else if(e.CheckRecoTopology(ccpi0))
            reco_topology = 3;
          else
            reco_topology = 4;
          if(e.CheckMCTopology(cc0pi))
            true_topology = 0;
          else if(e.CheckMCTopology(cc1pi))
            true_topology = 1;
          else if(e.CheckMCTopology(cc2pi))
            true_topology = 2;
          else if(e.CheckMCTopology(ccpi0))
            true_topology = 3;
          else
            true_topology = 4;
          t_run->Fill();
        } // Bad event
      } // Topology
      // Check all true CC Inclusive events and plot the reconstructed energy of them
      if(e.CheckMCTopology(ccinc)){
        if(bad_event) continue;
        double enu_mc_reco = 0.;
        for(const Particle &p : e.GetMCParticleList()){
          // Calculate the reconstructed energy from the reco kinetic (visible) + the true particle mass (cheating)
          double mass = 0.;
          if(p.GetPdgCode() == 211  || 
             p.GetPdgCode() == -211 ||
             p.GetPdgCode() == 13)
            mass = p.GetMass();
          enu_mc_reco += (p.GetKineticEnergy() + mass);
        } // Particles
        h_ccinc_e_true->Fill(enu_mc_reco);
        if(e.CheckRecoTopology(ccinc)){
          h_ccinc_e_sig->Fill(enu_mc_reco);
        } // Signal CCInclusive
      } // Topology
    } // Fiducial
  } // Events
  std::cout << " Total events with particles with bad energy : " << bad_events << std::endl;
  // Print the total pot from all the samples
  std::cout << " Total POT in the samples is: " << pot << std::endl;
  t_subrun->Fill();

  // Write histograms to file and external files
  TH1D *h_ccinc_e_eff = (TH1D*)h_ccinc_e_sig->Clone("h_ccinc_e_eff");
  h_ccinc_e_eff->Divide(h_ccinc_e_true);

  // Output TFile
  TFile f((file_location+"ccinc_selection_180620_"+std::to_string(start_file)+"-"+std::to_string(end_file)+".root").c_str(), "RECREATE");

  h_ccinc_e_true->Write();
  h_ccinc_e_sig->Write();
  h_ccinc_e_eff->Write();
  t_run->Write();
  t_subrun->Write();

  f.Write();
  f.Close();

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

} // MainTest

void LoadAllEvents(EventSelectionTool::EventList &events, 
                   const unsigned int &start_file, 
                   const unsigned int &end_file, 
                   const int &start_time, 
                   double &pot, 
                   std::vector<unsigned int> &exceptions){
  double total_pot = 0;
  std::vector<unsigned int>::iterator it;
  // Load the events into the event list
  for( unsigned int i = start_file; i < end_file; ++i ){
    it = std::find(exceptions.begin(), exceptions.end(),i);
    if(it != exceptions.end()) continue;
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    //name = "/sbnd/app/users/rsjones/LArSoft_v08_54_00/test/output_file.root";
    name = "/pnfs/sbnd/persistent/users/rsjones/sbnd_selection_170620/selection/"+std::to_string(i)+"/output_file.root";
    //name = "/pnfs/sbnd/persistent/users/rsjones/mcp0.9_neutrino_with_subrun/selection/"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    double temp_pot = 0.;
    unsigned int total_files = end_file - start_file;
    EventSelectionTool::LoadEventList(file_name, events, i, temp_pot);
    EventSelectionTool::GetTimeLeft(start_time,total_files,i-start_file);
    total_pot += temp_pot;
  }
  std::cout << std::endl;
  pot = total_pot;
} // LoadAllEvents
