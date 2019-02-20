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
#include "TROOT.h"
#include "TAxis.h"

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
  std::string plots_location  = "../Output_Selection_Tool/plots/vertices/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  // Maps
  TopologyMap cc0pi_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
 
  // Initiate histograms
  TH1D *h_truedirection_mu = new TH1D("h_truedirection_mu","", 50, -1, 1);
  TH1D *h_truedirection_pi = new TH1D("h_truedirection_pi","", 50, -1, 1);

  TH1D *h_truemomentum_mu = new TH1D("h_truemomentum_mu","", 50, 0, 2);
  TH1D *h_truemomentum_pi = new TH1D("h_truemomentum_pi","", 50, 0, 2);

  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
  
  // Loop over events and perform vertexing study
  for(const Event &e : events){

    // Particles 
    ParticleList mc = e.GetMCParticleList();

    // TPC track criteria
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;

    // Neutrino vertex is within the ficudial border
    if(e.IsSBNDTrueFiducial()){
      // True muon-pion comparison of z-angle
      for(const Particle &p : mc){
        if(p.GetPdgCode() == 13){
          h_truedirection_mu->Fill(p.GetCosTheta());
          h_truemomentum_mu->Fill(p.GetModulusMomentum());
        }
        if(p.GetPdgCode() == 211 || p.GetPdgCode() == -211){
          h_truedirection_pi->Fill(p.GetCosTheta());
          h_truemomentum_pi->Fill(p.GetModulusMomentum());
        }
      }
    }
  }

  h_truedirection_mu->Scale(1./double(h_truedirection_mu->Integral()));
  h_truedirection_pi->Scale(1./double(h_truedirection_pi->Integral()));
  h_truedirection_mu->SetFillColor(kSpring-3);
  h_truedirection_pi->SetFillColor(kOrange+7);
  h_truedirection_mu->SetFillStyle(3001);
  h_truedirection_pi->SetFillStyle(3001);
  h_truedirection_mu->SetLineColor(kSpring-3);
  h_truedirection_pi->SetLineColor(kOrange+7);

  h_truemomentum_mu->Scale(1./double(h_truemomentum_mu->Integral()));
  h_truemomentum_pi->Scale(1./double(h_truemomentum_pi->Integral()));
  h_truemomentum_mu->SetFillColor(kSpring-3);
  h_truemomentum_pi->SetFillColor(kOrange+7);
  h_truemomentum_mu->SetFillStyle(3001);
  h_truemomentum_pi->SetFillStyle(3001);
  h_truemomentum_mu->SetLineColor(kSpring-3);
  h_truemomentum_pi->SetLineColor(kOrange+7);

  TCanvas *c      = new TCanvas("c", "", 800, 600);
  TLegend *l      = new TLegend(0.22,0.68,0.52,0.88);
  
  // Legend
  l->AddEntry(h_truedirection_mu,  "#mu direction",  "f");
  l->AddEntry(h_truedirection_pi,  "#pi direction",  "f");
  l->SetLineWidth(0);
  l->SetTextAlign(12);
  l->SetTextFont(133);

  // Max y
  double max_direction = 1.1 * std::max(h_truedirection_mu->GetMaximum(), h_truedirection_pi->GetMaximum());

  // Draw plots
  h_truedirection_mu->SetStats(0);
  h_truedirection_mu->GetXaxis()->SetTitle("Track angle, cos(#theta)");
  h_truedirection_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_truedirection_mu->GetXaxis()->SetTitleFont(132);
  h_truedirection_mu->GetXaxis()->SetLabelFont(132);
  h_truedirection_mu->GetYaxis()->SetTitleFont(132);
  h_truedirection_mu->GetYaxis()->SetLabelFont(132);
  h_truedirection_mu->GetYaxis()->SetTitleOffset(1.2);
  h_truedirection_mu->GetYaxis()->SetMaxDigits(3);
  h_truedirection_mu->GetYaxis()->SetRangeUser(0, max_direction);
  h_truedirection_mu->Draw("hist");
  h_truedirection_pi->Draw("same hist");
  l->Draw("same");

  c->SaveAs((plots_location+"track_direction.root").c_str());
  c->Clear();
  l->Clear();

  // Legend
  TLegend *l_mom = new TLegend(0.58,0.68,0.88,0.88);
  l_mom->AddEntry(h_truemomentum_mu,  "#mu momentum",  "f");
  l_mom->AddEntry(h_truemomentum_pi,  "#pi momentum",  "f");
  l_mom->SetLineWidth(0);
  l_mom->SetTextAlign(12);
  l_mom->SetTextFont(133);

  // Max y
  double max_momentum = 1.1 * std::max(h_truemomentum_mu->GetMaximum(), h_truemomentum_pi->GetMaximum());

  // Draw plots
  h_truemomentum_mu->SetStats(0);
  h_truemomentum_mu->GetXaxis()->SetTitle("Track momentum, [GeV]");
  h_truemomentum_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_truemomentum_mu->GetXaxis()->SetTitleFont(132);
  h_truemomentum_mu->GetXaxis()->SetLabelFont(132);
  h_truemomentum_mu->GetYaxis()->SetTitleFont(132);
  h_truemomentum_mu->GetYaxis()->SetLabelFont(132);
  h_truemomentum_mu->GetYaxis()->SetTitleOffset(1.2);
  h_truemomentum_mu->GetYaxis()->SetMaxDigits(3);
  h_truemomentum_mu->GetYaxis()->SetRangeUser(0, max_momentum);
  h_truemomentum_mu->Draw("hist");
  h_truemomentum_pi->Draw("same hist");
  l_mom->Draw("same");

  c->SaveAs((plots_location+"track_momentum.root").c_str());
  c->Clear();
  l_mom->Clear();

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
