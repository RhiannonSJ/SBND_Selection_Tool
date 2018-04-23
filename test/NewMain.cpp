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
  std::cout << " Start: Local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string file_location = "../Output_Selection_Tool/plots/";

  // Counters
  unsigned int correctly_reconstructed_cc  = 0;
  unsigned int true_topology_cc            = 0;
  unsigned int reco_topology_cc            = 0;

  unsigned int correctly_reconstructed_nc  = 0;
  unsigned int true_topology_nc            = 0;
  unsigned int reco_topology_nc            = 0;

  unsigned int correctly_reconstructed_0pi = 0;
  unsigned int true_topology_0pi           = 0;
  unsigned int reco_topology_0pi           = 0;
  unsigned int true_ccqe                   = 0;
  unsigned int true_ccmec                  = 0;
  unsigned int true_ccdis                  = 0;
  unsigned int true_ccres                  = 0;
  unsigned int true_cccoh                  = 0;
  unsigned int reco_ccqe                   = 0;
  unsigned int reco_ccmec                  = 0;
  unsigned int reco_ccdis                  = 0;
  unsigned int reco_ccres                  = 0;
  unsigned int reco_cccoh                  = 0;
  unsigned int good_ccqe                   = 0;
  unsigned int good_ccmec                  = 0;
  unsigned int good_ccdis                  = 0;
  unsigned int good_ccres                  = 0;
  unsigned int good_cccoh                  = 0;
  unsigned int proton_multiplicity_reco    = 0;
  unsigned int proton_multiplicity_good    = 0;
  float proton_energy_sum_reco             = 0.;
  float proton_energy_sum_good             = 0.;

  unsigned int correctly_reconstructed_1pi = 0;
  unsigned int true_topology_1pi           = 0;
  unsigned int reco_topology_1pi           = 0;

  unsigned int correctly_reconstructed_pi0 = 0;
  unsigned int true_topology_pi0           = 0;
  unsigned int reco_topology_pi0           = 0;

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
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
  std::cout << " After loading events: Local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Neutrino energy histograms
  TH1F *h_reco_energy           = new TH1F("h_reco_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 3);
  TH1F *h_good_energy           = new TH1F("h_good_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 3);
  TH1F *h_true_energy           = new TH1F("h_true_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 3);
  TH2F *h_reco_energy_cos       = new TH2F("h_reco_energy_cos",      "#nu_{#mu} CC 0#pi neutino energy and opening angle", 50, 0, 3, 50, -1, 1);
  TH2F *h_good_energy_cos       = new TH2F("h_good_energy_cos",      "#nu_{#mu} CC 0#pi neutino energy and opening angle", 50, 0, 3, 50, -1, 1);
  TH2F *h_reco_energy_proton    = new TH2F("h_reco_energy_protons",  "#nu_{#mu} CC 0#pi neutino energy and proton multiplicity", 50, 0, 3, 5, 0, 5);
  TH2F *h_good_energy_proton    = new TH2F("h_good_energy_protons",  "#nu_{#mu} CC 0#pi neutino energy and proton multiplicity", 50, 0, 3, 5, 0, 5);
  TH2F *h_reco_energy_proton_e  = new TH2F("h_reco_energy_proton_e", "#nu_{#mu} CC 0#pi neutino energy and proton energy", 50, 0, 3, 50, 0, 3);
  TH2F *h_good_energy_proton_e  = new TH2F("h_good_energy_proton_e", "#nu_{#mu} CC 0#pi neutino energy and proton energy", 50, 0, 3, 50, 0, 3);

  // Get CC 0pi true and reconstructed neutrino energies
  std::vector<float> true_neutrino_energy, reco_neutrino_energy, good_neutrino_energy, good_cos_theta, reco_cos_theta, good_proton_multiplicity, reco_proton_multiplicity, good_proton_energy_sum, reco_proton_energy_sum;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc_signal_map) && e.CheckRecoTopology(cc_signal_map)) correctly_reconstructed_cc++;
    if(e.CheckRecoTopology(cc_signal_map)) reco_topology_cc++;
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc_signal_map))   true_topology_cc++;

    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(nc_signal_map) && e.CheckRecoTopology(nc_signal_map)) correctly_reconstructed_nc++;
    if(e.CheckRecoTopology(nc_signal_map)) reco_topology_nc++;
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(nc_signal_map))   true_topology_nc++;
    
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc1pi_signal_map) && e.CheckRecoTopology(cc1pi_signal_map)) correctly_reconstructed_1pi++;
    if(e.CheckRecoTopology(cc1pi_signal_map)) reco_topology_1pi++;
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc1pi_signal_map))   true_topology_1pi++;
    
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(ccpi0_signal_map) && e.CheckRecoTopology(ccpi0_signal_map)) correctly_reconstructed_pi0++;
    if(e.CheckRecoTopology(ccpi0_signal_map)) reco_topology_pi0++;
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(ccpi0_signal_map))   true_topology_pi0++;
    
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc0pi_signal_map) && e.CheckRecoTopology(cc0pi_signal_map)) {
      
      // Counter
      correctly_reconstructed_0pi++;
      
      // Get the well selected reconstructed neutrino energy
      ParticleList parts       = e.GetRecoParticleList();
      unsigned int n_particles = e.GetRecoParticleList().size();
      
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 2212) {
        proton_multiplicity_good++;
        proton_energy_sum_good += parts[i].GetEnergy();
      }
      good_proton_multiplicity.push_back(proton_multiplicity_good);
      good_proton_energy_sum.push_back(proton_energy_sum_good);
      
      // Get reconstructed energy
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 13) {
        good_neutrino_energy.push_back(CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e));

        TVector3 z;
        z[0] = 0;
        z[1] = 0;
        z[2] = 1;

        float p    = parts[i].GetMomentum().Mag();
        float cth  = (1/p) * (parts[i].GetMomentum()).Dot(z);
        good_cos_theta.push_back(cth);

      }
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) good_ccqe++;
        if(e.GetNuanceCode() == 10)   good_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) good_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) good_ccres++;
        if(e.GetNuanceCode() == 1097) good_cccoh++;
          
      }
    }
    if(e.CheckRecoTopology(cc0pi_signal_map)) {
      
      // Counter
      reco_topology_0pi++;

      // Get the well selected reconstructed neutrino energy
      ParticleList parts       = e.GetRecoParticleList();
      unsigned int n_particles = e.GetRecoParticleList().size();
      
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 2212){
        proton_multiplicity_reco++;
        proton_energy_sum_reco += parts[i].GetEnergy();
      }
      reco_proton_multiplicity.push_back(proton_multiplicity_reco);
      reco_proton_energy_sum.push_back(proton_energy_sum_reco);
      
      // Get reconstructed energy
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 13) {
        reco_neutrino_energy.push_back(CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e));
        
        TVector3 z;
        z[0] = 0;
        z[1] = 0;
        z[2] = 1;

        float p    = parts[i].GetMomentum().Mag();
        float cth  = (1/p) * (parts[i].GetMomentum()).Dot(z);
        reco_cos_theta.push_back(cth);
      }
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) reco_ccqe++;
        if(e.GetNuanceCode() == 10)   reco_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) reco_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) reco_ccres++;
        if(e.GetNuanceCode() == 1097) reco_cccoh++;
          
      }
    }
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(cc0pi_signal_map)) {
     
      // Counter
      true_topology_0pi++;
    
      // Get true neutrino energy
      true_neutrino_energy.push_back(e.GetTrueNuEnergy());
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) true_ccqe++;
        if(e.GetNuanceCode() == 10)   true_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) true_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) true_ccres++;
        if(e.GetNuanceCode() == 1097) true_cccoh++;
          
      }
    }
  }
  for( unsigned int i = 0; i < true_neutrino_energy.size(); ++i) h_true_energy->Fill(true_neutrino_energy[i]);
  for( unsigned int i = 0; i < reco_neutrino_energy.size(); ++i) {
    h_reco_energy->Fill(reco_neutrino_energy[i]);
    h_reco_energy_proton->Fill(reco_neutrino_energy[i],reco_proton_multiplicity[i]);
    h_reco_energy_proton_e->Fill(reco_neutrino_energy[i],reco_proton_energy_sum[i]);
    h_reco_energy_cos->Fill(reco_neutrino_energy[i],reco_cos_theta[i]);
  }
  for( unsigned int i = 0; i < good_neutrino_energy.size(); ++i) {
    h_good_energy->Fill(good_neutrino_energy[i]);
    h_good_energy_proton->Fill(good_neutrino_energy[i],good_proton_multiplicity[i]);
    h_good_energy_proton_e->Fill(good_neutrino_energy[i],good_proton_energy_sum[i]);
    h_good_energy_cos->Fill(good_neutrino_energy[i],good_cos_theta[i]);
  }

  std::cout << "===========================================================" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed NC events     : " << correctly_reconstructed_nc/double(true_topology_nc) * 100                         << std::endl;
  std::cout << " Purity of reconstructed NC events                   : " << (correctly_reconstructed_nc)/double(reco_topology_nc) * 100    << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC events     : " << correctly_reconstructed_cc/double(true_topology_cc) * 100                         << std::endl;
  std::cout << " Purity of reconstructed CC events                   : " << (correctly_reconstructed_cc)/double(reco_topology_cc) * 100    << std::endl;
  std::cout << " CC signal                                           : " << correctly_reconstructed_cc << std::endl;
  std::cout << " CC selected                                         : " << reco_topology_cc           << std::endl;
  std::cout << " CC true                                             : " << true_topology_cc           << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC 0Pi events : " << correctly_reconstructed_0pi/double(true_topology_0pi) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC 0Pi events               : " << (correctly_reconstructed_0pi)/double(reco_topology_0pi) * 100 << std::endl;
  std::cout << " CC 0pi signal                                       : " << correctly_reconstructed_0pi << std::endl;
  std::cout << " CC 0pi selected                                     : " << reco_topology_0pi           << std::endl;
  std::cout << " CC 0pi true                                         : " << true_topology_0pi           << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC 1Pi events : " << correctly_reconstructed_1pi/double(true_topology_1pi) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC 1Pi events               : " << (correctly_reconstructed_1pi)/double(reco_topology_1pi) * 100 << std::endl; 
  std::cout << " CC 1pi signal                                       : " << correctly_reconstructed_1pi << std::endl;
  std::cout << " CC 1pi selected                                     : " << reco_topology_1pi           << std::endl;
  std::cout << " CC 1pi true                                         : " << true_topology_1pi           << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC Pi0 events : " << correctly_reconstructed_pi0/double(true_topology_pi0) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC Pi0 events               : " << (correctly_reconstructed_pi0)/double(reco_topology_pi0) * 100 << std::endl; 
  std::cout << " CC pi0 signal                                       : " << correctly_reconstructed_pi0 << std::endl;
  std::cout << " CC pi0 selected                                     : " << reco_topology_pi0           << std::endl;
  std::cout << " CC pi0 true                                         : " << true_topology_pi0           << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "===========================================================" << std::endl;

  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.58, 0.68, 0.88, 0.88 );
  l->SetBorderSize(0);

  l->AddEntry( h_true_energy,      " True Monte Carlo",    "l" );
  l->AddEntry( h_reco_energy,      " Selected ",           "l" );
  l->AddEntry( h_good_energy,      " Correctly selected ", "l" );
    
  h_true_energy->SetLineColor(2);
  h_true_energy->SetStats(kFALSE);
  h_reco_energy->SetLineColor(4);
  h_reco_energy->SetStats(kFALSE);
  h_reco_energy->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  h_good_energy->SetLineColor(6);
  h_good_energy->SetStats(kFALSE);

  h_reco_energy->Draw();
  h_good_energy->Draw("same");
  h_true_energy->Draw("same");
  l->Draw();

  c->SaveAs((file_location+"cc0pi_nu_energy.root").c_str());
  c->Clear();

  h_reco_energy_cos->SetStats(kFALSE);
  h_reco_energy_cos->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_cos->GetYaxis()->SetTitle("cos#theta_{#mu}");
  h_reco_energy_cos->Draw("colz");

  c->SaveAs((file_location+"cc0pi_reco_energy_cos.root").c_str());
  c->Clear();
  
  h_good_energy_cos->SetStats(kFALSE);
  h_good_energy_cos->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_cos->GetYaxis()->SetTitle("cos#theta_{#mu}");
  h_good_energy_cos->Draw("colz");
 
  c->SaveAs((file_location+"cc0pi_good_energy_cos.root").c_str());
  c->Clear();
  
  h_reco_energy_proton->SetStats(kFALSE);
  h_reco_energy_proton->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_proton->GetYaxis()->SetTitle("Proton multiplicity");
  h_reco_energy_proton->Draw("colz");
 
  c->SaveAs((file_location+"cc0pi_reco_energy_proton.root").c_str());
  c->Clear();
  
  h_good_energy_proton->SetStats(kFALSE);
  h_good_energy_proton->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_proton->GetYaxis()->SetTitle("Proton multiplicity");
  h_good_energy_proton->Draw("colz");
 
  c->SaveAs((file_location+"cc0pi_good_energy_proton.root").c_str());
  c->Clear();
  
  h_reco_energy_proton_e->SetStats(kFALSE);
  h_reco_energy_proton_e->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_proton_e->GetYaxis()->SetTitle("Proton energy sum [GeV]");
  h_reco_energy_proton_e->Draw("colz");
 
  c->SaveAs((file_location+"cc0pi_reco_proton_energies.root").c_str());
  c->Clear();
  
  h_good_energy_proton_e->SetStats(kFALSE);
  h_good_energy_proton_e->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_proton_e->GetYaxis()->SetTitle("Proton energy sum [GeV]");
  h_good_energy_proton_e->Draw("colz");
 
  c->SaveAs((file_location+"cc0pi_good_proton_energies.root").c_str());
  c->Clear();
 
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
