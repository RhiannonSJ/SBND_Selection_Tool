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
  std::cout << " Running all files " << std::endl;
  
  // Counters
  unsigned int correctly_reconstructed_cc     = 0;
  unsigned int true_topology_cc               = 0;
  unsigned int reco_topology_cc               = 0;

  unsigned int correctly_reconstructed_other  = 0;
  unsigned int true_topology_other            = 0;
  unsigned int reco_topology_other            = 0;

  unsigned int correctly_reconstructed_nc     = 0;
  unsigned int true_topology_nc               = 0;
  unsigned int reco_topology_nc               = 0;
                                              
  unsigned int correctly_reconstructed_0pi    = 0;
  unsigned int true_topology_0pi              = 0;
  unsigned int reco_topology_0pi              = 0;
  unsigned int true_ccqe                      = 0;
  unsigned int true_ccmec                     = 0;
  unsigned int true_ccdis                     = 0;
  unsigned int true_ccres                     = 0;
  unsigned int true_cccoh                     = 0;
  unsigned int reco_ccqe                      = 0;
  unsigned int reco_ccmec                     = 0;
  unsigned int reco_ccdis                     = 0;
  unsigned int reco_ccres                     = 0;
  unsigned int reco_cccoh                     = 0;
  unsigned int good_ccqe                      = 0;
  unsigned int good_ccmec                     = 0;
  unsigned int good_ccdis                     = 0;
  unsigned int good_ccres                     = 0;
  unsigned int good_cccoh                     = 0;
  unsigned int proton_multiplicity_reco       = 0;
  unsigned int proton_multiplicity_good       = 0;
  float proton_energy_sum_reco                = 0.;
  float proton_energy_sum_good                = 0.;
                                              
  unsigned int correctly_reconstructed_1pi    = 0;
  unsigned int true_topology_1pi              = 0;
  unsigned int reco_topology_1pi              = 0;
                                              
  unsigned int correctly_reconstructed_pi0    = 0;
  unsigned int true_topology_pi0              = 0;
  unsigned int reco_topology_pi0              = 0;

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc_signal_map;
  TopologyMap nc_signal_map;
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
  
  cc_signal_map.insert( std::make_pair( mu,  1 ) );
  nc_signal_map.insert( std::make_pair( mu,  0 ) );
  
  cc0pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc0pi_signal_map.insert( std::make_pair( pi,  0 ) );
  cc0pi_signal_map.insert( std::make_pair( pi0, 0 ) );

  cc1pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc1pi_signal_map.insert( std::make_pair( pi,  1 ) );
  cc1pi_signal_map.insert( std::make_pair( pi0, 0 ) );
  
  ccpi0_signal_map.insert( std::make_pair( mu,  1 ) );
  ccpi0_signal_map.insert( std::make_pair( pi,  0 ) );
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

  // Neutrino energy histograms
  TH1F *h_energy_eff            = new TH1F("h_energy_eff",           "#nu_{#mu} CC 0#pi neutrino energy efficiency",50,0, 2.5);
  TH1F *h_reco_energy           = new TH1F("h_reco_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 2.5);
  TH1F *h_good_energy           = new TH1F("h_good_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 2.5);
  TH1F *h_true_energy           = new TH1F("h_true_energy",          "#nu_{#mu} CC 0#pi neutrino energy",100,0, 2.5);
  TH1F *h_reco_energy_mc        = new TH1F("h_reco_energy_mc",       "#nu_{#mu} CC 0#pi Reco - MC neutrino energy",100,-2, 2.5);
  TH1F *h_good_energy_mc        = new TH1F("h_good_energy_mc",       "#nu_{#mu} CC 0#pi Reco - MC neutrino energy",100,-2, 2.5);
  TH1F *h_true_energy_mc        = new TH1F("h_true_energy_mc",       "#nu_{#mu} CC 0#pi Reco - MC neutrino energy",100,-2, 2.5);
  TH1F *h_reco_energy_mu        = new TH1F("h_reco_energy_mu",       "#nu_{#mu} CC 0#pi muon energy",    100,0, 2.5);
  TH1F *h_good_energy_mu        = new TH1F("h_good_energy_mu",       "#nu_{#mu} CC 0#pi muon energy",    100,0, 2.5);
  TH1F *h_true_energy_mu        = new TH1F("h_true_energy_mu",       "#nu_{#mu} CC 0#pi muon energy",    100,0, 2.5);
  TH2F *h_reco_energy_cos       = new TH2F("h_reco_energy_cos",      "#nu_{#mu} CC 0#pi neutino energy and opening angle", 50, 0, 3, 50, -1, 1);
  TH2F *h_good_energy_cos       = new TH2F("h_good_energy_cos",      "#nu_{#mu} CC 0#pi neutino energy and opening angle", 50, 0, 3, 50, -1, 1);
  TH2F *h_reco_energy_proton    = new TH2F("h_reco_energy_protons",  "#nu_{#mu} CC 0#pi neutino energy and proton multiplicity", 50, 0, 3, 5, 0, 5);
  TH2F *h_good_energy_proton    = new TH2F("h_good_energy_protons",  "#nu_{#mu} CC 0#pi neutino energy and proton multiplicity", 50, 0, 3, 5, 0, 5);
  TH2F *h_reco_energy_proton_e  = new TH2F("h_reco_energy_proton_e", "#nu_{#mu} CC 0#pi neutino energy and proton energy", 50, 0, 3, 50, 0, 3);
  TH2F *h_good_energy_proton_e  = new TH2F("h_good_energy_proton_e", "#nu_{#mu} CC 0#pi neutino energy and proton energy", 50, 0, 3, 50, 0, 3);

  // Get CC 0pi true and reconstructed neutrino energies
  std::vector<float> mc_true_neutrino_energy, true_neutrino_energy, reco_neutrino_energy, good_neutrino_energy, good_cos_theta, reco_cos_theta, good_proton_multiplicity, reco_proton_multiplicity, good_proton_energy_sum, reco_proton_energy_sum, true_muon_energy, reco_muon_energy, good_muon_energy, efficiency_energy;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.CheckMCTopology(cc_signal_map) && e.CheckRecoTopology(cc_signal_map)) correctly_reconstructed_cc++;
    if(e.CheckRecoTopology(cc_signal_map)) reco_topology_cc++;
    if(e.CheckMCTopology(cc_signal_map))   true_topology_cc++;

    if(e.CheckMCTopology(nc_signal_map) && e.CheckRecoTopology(nc_signal_map)) correctly_reconstructed_nc++;
    if(e.CheckRecoTopology(nc_signal_map)) reco_topology_nc++;
    if(e.CheckMCTopology(nc_signal_map))   true_topology_nc++;
    
    if(e.CheckMCTopology(cc1pi_signal_map) && e.CheckRecoTopology(cc1pi_signal_map)) correctly_reconstructed_1pi++;
    if(e.CheckRecoTopology(cc1pi_signal_map)) reco_topology_1pi++;
    if(e.CheckMCTopology(cc1pi_signal_map))   true_topology_1pi++;
    
    if(e.CheckMCTopology(ccpi0_signal_map) && e.CheckRecoTopology(ccpi0_signal_map)) correctly_reconstructed_pi0++;
    if(e.CheckRecoTopology(ccpi0_signal_map)) reco_topology_pi0++;
    if(e.CheckMCTopology(ccpi0_signal_map))   true_topology_pi0++;
    
    if(e.CheckMCTopology(cc0pi_signal_map) && e.CheckRecoTopology(cc0pi_signal_map)) {
      
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
        good_neutrino_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]));
        good_muon_energy.push_back(parts[i].GetEnergy());

        efficiency_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]) / e.GetTrueNuEnergy());

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
        reco_muon_energy.push_back(parts[i].GetEnergy());
        reco_neutrino_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]));
        
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
    if(e.CheckMCTopology(cc0pi_signal_map)) {
     
      // Counter
      true_topology_0pi++;
    
      ParticleList parts       = e.GetMCParticleList();
      unsigned int n_particles = e.GetMCParticleList().size();

      // Get the MC true neutrino energy
      mc_true_neutrino_energy.push_back(e.GetTrueNuEnergy());

      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 13) {
        true_neutrino_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]));
        true_muon_energy.push_back(parts[i].GetEnergy());
      }

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
  for( unsigned int i = 0; i < true_neutrino_energy.size(); ++i) {
    h_true_energy->Fill(true_neutrino_energy[i]);
    h_true_energy_mc->Fill(true_neutrino_energy[i] - mc_true_neutrino_energy[i]);
    h_true_energy_mu->Fill(true_muon_energy[i]);
  }
  for( unsigned int i = 0; i < reco_neutrino_energy.size(); ++i) {
    h_reco_energy->Fill(reco_neutrino_energy[i]);
    h_reco_energy_mc->Fill(reco_neutrino_energy[i] - mc_true_neutrino_energy[i]);
    h_reco_energy_mu->Fill(reco_muon_energy[i]);
    h_reco_energy_proton->Fill(reco_neutrino_energy[i],reco_proton_multiplicity[i]);
    h_reco_energy_proton_e->Fill(reco_neutrino_energy[i],reco_proton_energy_sum[i]);
    h_reco_energy_cos->Fill(reco_neutrino_energy[i],reco_cos_theta[i]);
  }
  for( unsigned int i = 0; i < good_neutrino_energy.size(); ++i) {
    h_energy_eff->Fill(good_neutrino_energy[i],efficiency_energy[i]);
    h_good_energy->Fill(good_neutrino_energy[i]);
    h_good_energy_mu->Fill(good_muon_energy[i]);
    h_good_energy_mc->Fill(good_neutrino_energy[i] - mc_true_neutrino_energy[i]);
    h_good_energy_proton->Fill(good_neutrino_energy[i],good_proton_multiplicity[i]);
    h_good_energy_proton_e->Fill(good_neutrino_energy[i],good_proton_energy_sum[i]);
    h_good_energy_cos->Fill(good_neutrino_energy[i],good_cos_theta[i]);
  }

  std::cout << "===========================================================" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed NC events     : " << correctly_reconstructed_nc/double(true_topology_nc) * 100                         << std::endl;
  std::cout << " Purity of reconstructed NC events                   : " << (correctly_reconstructed_nc)/double(reco_topology_nc) * 100    << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed CC events     : " << correctly_reconstructed_cc/double(true_topology_cc) * 100                         << std::endl;
  std::cout << " Purity of reconstructed CC events                   : " << (correctly_reconstructed_cc)/double(reco_topology_cc) * 100    << std::endl; 
  std::cout << " Number of CC events reconstructed                   : " << reco_topology_cc << std::endl; 
  std::cout << " Number of CC events correctly reconstructed         : " << correctly_reconstructed_cc << std::endl; 
  std::cout << " Number of CC events true                            : " << true_topology_cc << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed CC 0Pi events : " << correctly_reconstructed_0pi/double(true_topology_0pi) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC 0Pi events               : " << (correctly_reconstructed_0pi)/double(reco_topology_0pi) * 100 << std::endl;
  std::cout << " Number of CC 0pi events reconstructed               : " << reco_topology_0pi << std::endl; 
  std::cout << " Number of CC 0pi events correctly reconstructed     : " << correctly_reconstructed_0pi << std::endl; 
  std::cout << " Number of CC 0pi events true                        : " << true_topology_0pi << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed CC 1Pi events : " << correctly_reconstructed_1pi/double(true_topology_1pi) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC 1Pi events               : " << (correctly_reconstructed_1pi)/double(reco_topology_1pi) * 100 << std::endl; 
  std::cout << " Number of CC 1pi events reconstructed               : " << reco_topology_1pi << std::endl; 
  std::cout << " Number of CC 1pi events correctly reconstructed     : " << correctly_reconstructed_1pi << std::endl; 
  std::cout << " Number of CC 1pi events true                        : " << true_topology_1pi << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed CC Pi0 events : " << correctly_reconstructed_pi0/double(true_topology_pi0) * 100                       << std::endl;
  std::cout << " Purity of reconstructed CC Pi0 events               : " << (correctly_reconstructed_pi0)/double(reco_topology_pi0) * 100 << std::endl; 
  std::cout << " Number of CC pi0 events reconstructed               : " << reco_topology_pi0 << std::endl; 
  std::cout << " Number of CC pi0 events correctly reconstructed     : " << correctly_reconstructed_pi0 << std::endl; 
  std::cout << " Number of CC pi0 events true                        : " << true_topology_pi0 << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Efficiency of correctly reconstructed CC Other      : " << correctly_reconstructed_other/double(true_topology_other) * 100                         << std::endl;
  std::cout << " Purity of reconstructed CC Other events             : " << (correctly_reconstructed_other)/double(reco_topology_other) * 100    << std::endl; 
  std::cout << " Number of CC Other events reconstructed             : " << reco_topology_other << std::endl; 
  std::cout << " Number of CC Other events correctly reconstructed   : " << correctly_reconstructed_other << std::endl; 
  std::cout << " Number of CC Other events true                      : " << true_topology_other << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "===========================================================" << std::endl;

  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.58, 0.68, 0.88, 0.88 );
  l->SetBorderSize(0);

  l->AddEntry( h_true_energy,      " True Monte Carlo",    "l" );
  l->AddEntry( h_reco_energy,      " Selected ",           "l" );
  l->AddEntry( h_good_energy,      " Signal ",             "l" );
    
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

  c->SaveAs("plots/cc0pi_nu_energy.root");
  c->Clear();

  h_true_energy_mc->SetLineColor(2);
  h_true_energy_mc->SetStats(kFALSE);
  h_reco_energy_mc->SetLineColor(4);
  h_reco_energy_mc->SetStats(kFALSE);
  h_reco_energy_mc->GetXaxis()->SetTitle("Reconstructed - MC Neutrino energy [GeV]");
  h_good_energy_mc->SetLineColor(6);
  h_good_energy_mc->SetStats(kFALSE);

  h_true_energy_mc->Draw();
  h_good_energy_mc->Draw("same");
  h_reco_energy_mc->Draw("same");
  l->Draw();

  c->SaveAs("plots/cc0pi_nu_energy_mc.root");
  c->Clear();

  h_true_energy_mu->SetLineColor(2);
  h_true_energy_mu->SetStats(kFALSE);
  h_true_energy_mu->GetXaxis()->SetTitle("Muon Energy [GeV]");
  h_reco_energy_mu->SetLineColor(4);
  h_reco_energy_mu->SetStats(kFALSE);
  h_good_energy_mu->SetLineColor(6);
  h_good_energy_mu->SetStats(kFALSE);

  h_reco_energy_mu->Draw();
  h_good_energy_mu->Draw("same");
  h_true_energy_mu->Draw("same");
  l->Draw();

  c->SaveAs("plots/cc0pi_mu_energy.root");
  c->Clear();

  h_energy_eff->SetStats(kFALSE);
  h_energy_eff->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  h_energy_eff->GetYaxis()->SetTitle("Efficiency");

  h_energy_eff->Draw();

  c->SaveAs("plots/cc0pi_nu_energy_eff.root");
  c->Clear();

  h_reco_energy_cos->SetStats(kFALSE);
  h_reco_energy_cos->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_cos->GetYaxis()->SetTitle("cos#theta_{#mu}");
  h_reco_energy_cos->Draw("colz");

  c->SaveAs("plots/cc0pi_reco_energy_cos.root");
  c->Clear();
  
  h_good_energy_cos->SetStats(kFALSE);
  h_good_energy_cos->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_cos->GetYaxis()->SetTitle("cos#theta_{#mu}");
  h_good_energy_cos->Draw("colz");
 
  c->SaveAs("plots/cc0pi_good_energy_cos.root");
  c->Clear();
  
  h_reco_energy_proton->SetStats(kFALSE);
  h_reco_energy_proton->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_proton->GetYaxis()->SetTitle("Proton multiplicity");
  h_reco_energy_proton->Draw("colz");
 
  c->SaveAs("plots/cc0pi_reco_energy_proton.root");
  c->Clear();
  
  h_good_energy_proton->SetStats(kFALSE);
  h_good_energy_proton->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_proton->GetYaxis()->SetTitle("Proton multiplicity");
  h_good_energy_proton->Draw("colz");
 
  c->SaveAs("plots/cc0pi_good_energy_proton.root");
  c->Clear();
  
  h_reco_energy_proton_e->SetStats(kFALSE);
  h_reco_energy_proton_e->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_reco_energy_proton_e->GetYaxis()->SetTitle("Proton energy sum [GeV]");
  h_reco_energy_proton_e->Draw("colz");
 
  c->SaveAs("plots/cc0pi_reco_proton_energies.root");
  c->Clear();
  
  h_good_energy_proton_e->SetStats(kFALSE);
  h_good_energy_proton_e->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_good_energy_proton_e->GetYaxis()->SetTitle("Proton energy sum [GeV]");
  h_good_energy_proton_e->Draw("colz");
 
  c->SaveAs("plots/cc0pi_good_proton_energies.root");
  c->Clear();
 
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
