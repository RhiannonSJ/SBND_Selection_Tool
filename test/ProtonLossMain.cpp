#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <time.h>
#include <algorithm>
#include <vector>
#include <string>
#include "TVector3.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"

using namespace selection;

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<int> &exceptions);

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string plots = "/sbnd/data/users/rsjones/Output_Selection_Tool/plots/mcp2/cut75cm/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 824;
  double pot = 0.; 

  std::vector<int> exceptions;
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

  LoadAllEvents(events, total_files, start, pot, exceptions);
  
  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  /*
   * Proton loss and mis-identification energy and length study
   *
   * Loop over events
   *   Loop over list of reconstructed particles
   *      If true == 2212 && reco == 2212: Matched
   *        Get kinetic energy of true
   *        Get length of true
   *      If reco != 2212 && true == 2212: Mis-identified true proton
   *        Get kinetic energy of true
   *        Get length of true
   *   Loop over list of MC particles
   *      If true == 2212
   *        Loop over reconstructed particles and see if any have corresponding MC ID
   *          If not: missed proton
   *            Get kinetic energy of true
   *            Get length of true
   *  
   */
  TopologyMap cc    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0 = GeneralAnalysisHelper::GetCCPi0TopologyMap();

  std::vector<float> signal_energy, signal_length, background_energy, background_length, missed_energy, missed_length;
  std::vector<int>   signal_hits, background_hits, missed_hits;
  unsigned int total           = 0;
  unsigned int missed          = 0;
  unsigned int missed_below_25 = 0;
  unsigned int missed_below_10 = 0;

  // Proton multiplicity counters
  unsigned int ccqe  = 0;
  unsigned int ccmec = 0;
  unsigned int ccres = 0;
  unsigned int ccdis = 0;
  unsigned int ccoth = 0;
  unsigned int ncany = 0;
  
  unsigned int ccqe_reco = 0;
  unsigned int ccmec_reco = 0;
  unsigned int ccres_reco = 0;
  unsigned int ccdis_reco = 0;
  unsigned int ccoth_reco = 0;
  unsigned int ncany_reco = 0;
  
  TH1D *h_ccqe  = new TH1D("h_ccqe",  "Proton multiplicity", 10, 0, 9);
  TH1D *h_ccmec = new TH1D("h_ccmec", "Proton multiplicity", 10, 0, 9);
  TH1D *h_ccres = new TH1D("h_ccres", "Proton multiplicity", 10, 0, 9);
  TH1D *h_ccdis = new TH1D("h_ccdis", "Proton multiplicity", 10, 0, 9);
  TH1D *h_ccoth = new TH1D("h_ccoth", "Proton multiplicity", 10, 0, 9);
  TH1D *h_nc    = new TH1D("h_nc",    "Proton multiplicity", 10, 0, 9);

  TH1D *hr_ccqe  = new TH1D("hr_ccqe",  "Proton multiplicity", 10, 0, 9);
  TH1D *hr_ccmec = new TH1D("hr_ccmec", "Proton multiplicity", 10, 0, 9);
  TH1D *hr_ccres = new TH1D("hr_ccres", "Proton multiplicity", 10, 0, 9);
  TH1D *hr_ccdis = new TH1D("hr_ccdis", "Proton multiplicity", 10, 0, 9);
  TH1D *hr_ccoth = new TH1D("hr_ccoth", "Proton multiplicity", 10, 0, 9);
  TH1D *hr_nc    = new TH1D("hr_nc",    "Proton multiplicity", 10, 0, 9);

  for(const Event &e : events){

    if(!e.IsSBNDTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) continue;
    
    ParticleList reco_particles = e.GetRecoParticleList(); 
    ParticleList true_particles = e.GetMCParticleList();

    for(Particle &p_reco : reco_particles){
      // If the particle is a reconstructed track (not pi0)
      if(p_reco.GetFromRecoTrack() && GeneralAnalysisHelper::ParticleHasAMatch(e, p_reco) >= 0){
        Particle p_match = GeneralAnalysisHelper::GetBestMCParticle(e,p_reco);
        if(p_match.GetPdgCode() == 2212 && p_match.GetNumberOfHits() >= 5){
          if(p_reco.GetPdgCode() == 2212 && p_reco.GetNumberOfHits() >= 5){
            // True and reconstructed proton: signal
            signal_energy.push_back(p_match.GetKineticEnergy());
            signal_length.push_back(p_match.GetLength());
            signal_hits.push_back(p_match.GetNumberOfHits());
          }
          else if(p_reco.GetPdgCode() != 2212 && p_reco.GetNumberOfHits() >= 5){
            // True and mis-indentified proton
            background_energy.push_back(p_match.GetKineticEnergy());
            background_length.push_back(p_match.GetLength());
            background_hits.push_back(p_match.GetNumberOfHits());
          }
        }
      }
    }
    for(Particle &p_true : true_particles){
      if(p_true.GetPdgCode() == 2212 && !GeneralAnalysisHelper::HasBeenReconstructed(e, p_true) && p_true.GetNumberOfHits() >= 5){
        // True proton, not reconstructed
        missed_energy.push_back(p_true.GetKineticEnergy());
        missed_length.push_back(p_true.GetLength());
        missed_hits.push_back(p_true.GetNumberOfHits());
        missed++;
        if(p_true.GetNumberOfHits() < 25) missed_below_25++;
        if(p_true.GetNumberOfHits() < 10) missed_below_10++;
      }
      else if(p_true.GetPdgCode() == 2212 && GeneralAnalysisHelper::HasBeenReconstructed(e, p_true) && p_true.GetNumberOfHits() >= 5)
        total++;
    }
    if(e.CheckMCTopology(cc0pi)){
      for(const Particle &mc : true_particles){
        if(mc.GetPdgCode() == 2212 && mc.GetNumberOfHits() >= 5){
          if(e.GetPhysicalProcess() == 1 && e.GetIsCC())  ccqe++;
          else if(e.GetPhysicalProcess() == 3 && e.GetIsCC())  ccdis++;
          else if(e.GetPhysicalProcess() == 4 && e.GetIsCC())  ccres++;
          else if(e.GetPhysicalProcess() == 10 && e.GetIsCC()) ccmec++;
          else if(e.GetIsCC()) ccoth++;
          if(!e.GetIsCC()) ncany++;
        }
      }
      h_ccqe->Fill(ccqe);
      h_ccmec->Fill(ccmec);
      h_ccres->Fill(ccres);
      h_ccdis->Fill(ccdis);
      h_ccoth->Fill(ccoth);
      h_nc->Fill(ncany);
    }
    if(e.CheckRecoTopology(cc0pi)){
      for(const Particle &reco : reco_particles){
        if(reco.GetPdgCode() == 2212 && reco.GetNumberOfHits() >= 5){
          if(e.GetPhysicalProcess() == 1 && e.GetIsCC())  ccqe_reco++;
          else if(e.GetPhysicalProcess() == 3 && e.GetIsCC())  ccdis_reco++;
          else if(e.GetPhysicalProcess() == 4 && e.GetIsCC())  ccres_reco++;
          else if(e.GetPhysicalProcess() == 10 && e.GetIsCC()) ccmec_reco++;
          else if(e.GetIsCC()) ccoth_reco++;
          if(!e.GetIsCC()) ncany_reco++;
        }
      }
      hr_ccqe->Fill(ccqe_reco);
      hr_ccmec->Fill(ccmec_reco);
      hr_ccres->Fill(ccres_reco);
      hr_ccdis->Fill(ccdis_reco);
      hr_ccoth->Fill(ccoth_reco);
      hr_nc->Fill(ncany_reco);
    }
  }
  std::vector<TH1D*> signal;
  signal.push_back(h_ccqe);
  signal.push_back(h_ccmec);
  signal.push_back(h_ccres);
  signal.push_back(h_ccdis);
  signal.push_back(h_ccoth);
  signal.push_back(h_nc);

  std::vector<TH1D*> signal_reco;
  signal_reco.push_back(hr_ccqe);
  signal_reco.push_back(hr_ccmec);
  signal_reco.push_back(hr_ccres);
  signal_reco.push_back(hr_ccdis);
  signal_reco.push_back(hr_ccoth);
  signal_reco.push_back(hr_nc);

  /*
   *
   * Fill
   *
   */
  TH1D *h_signal_energy = new TH1D("h_signal_energy","Correctly reconstructed proton kinetic energies",40,0,0.6);
  TH1D *h_signal_length = new TH1D("h_signal_length","Correctly reconstructed proton lengths",40,0,15);
  TH1D *h_signal_hits   = new TH1D("h_signal_hits",  "Correctly reconstructed proton hits",40,0,40);

  TH1D *h_missed_energy = new TH1D("h_missed_energy","Missed proton kinetic energies",40,0,0.6);
  TH1D *h_missed_length = new TH1D("h_missed_length","Missed proton lengths",40,0,15);
  TH1D *h_missed_hits   = new TH1D("h_missed_hits",  "Missed proton hits",40,0,40);

  TH1D *h_background_energy = new TH1D("h_background_energy","Mis-identified proton kinetic energies",40,0,0.6);
  TH1D *h_background_length = new TH1D("h_background_length","Mis-identified proton lengths",40,0,15);
  TH1D *h_background_hits   = new TH1D("h_background_hits",  "Mis-identified proton hits",40,0,40);


  for(unsigned int i = 0; i < signal_energy.size(); ++i){
    h_signal_energy->Fill(signal_energy[i]);
    h_signal_length->Fill(signal_length[i]);
    h_signal_hits->Fill(signal_hits[i]);
  }
  for(unsigned int i = 0; i < background_energy.size(); ++i){
    h_background_energy->Fill(background_energy[i]);
    h_background_length->Fill(background_length[i]);
    h_background_hits->Fill(background_hits[i]);
  }
  for(unsigned int i = 0; i < missed_energy.size(); ++i){
    h_missed_energy->Fill(missed_energy[i]);
    h_missed_length->Fill(missed_length[i]);
    h_missed_hits->Fill(missed_hits[i]);
  }

  /*
   *
   * Stack
   *
   */
  TCanvas *c = new TCanvas("c","",700,900);
  TLegend *l = new TLegend( 0.68, 0.53, 0.88, 0.88 );
  THStack *hsT = new THStack("hsT","True CC 0#pi proton multiplicity pre-FSI");
  THStack *hsR = new THStack("hsR","Reco CC 0#pi proton multiplicity pre-FSI");

  int pal[6];
  pal[0]  = kPink+5;
  pal[1]  = kOrange+1;
  pal[2]  = kAzure+7;
  pal[3]  = kTeal-5;
  pal[4]  = kYellow-3;
  pal[5]  = kViolet+6;

  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTitleFont(132, "X");
  gStyle->SetTitleFont(132, "Y");
  gStyle->SetLabelFont(132, "X");
  gStyle->SetLabelFont(132, "Y");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(6, pal);
  
  double norm      = 0;
  double norm_reco = 0;
  
  assert(signal.size() == signal_reco.size());
  for (unsigned int i = 0; i < signal.size(); ++i){
      norm      += signal[i]->Integral();
      norm_reco += signal_reco[i]->Integral();
  }

  for ( unsigned int i = 0; i < signal.size(); ++i ){
    if(norm > 0){
      signal[i]->SetFillColor(pal[i]);
      signal[i]->SetFillStyle(3001);
      signal[i]->SetLineColor(pal[i]);
      signal[i]->Scale(1/norm);
    }
    if(norm_reco > 0){
      signal_reco[i]->SetFillColor(pal[i]);
      signal_reco[i]->SetFillStyle(3001);
      signal_reco[i]->SetLineColor(pal[i]);
      signal_reco[i]->Scale(1/norm_reco);
    }
  }
  
  for ( unsigned int i = 0; i < signal.size(); ++i ){
    hsT->Add(signal[i]);
    hsR->Add(signal_reco[i]);
  }


  /*
   *
   * Draw
   *
   */
  l->AddEntry(h_ccqe,  "CC QE",    "f");
  l->AddEntry(h_ccres, "CC RES",   "f");
  l->AddEntry(h_ccmec, "CC MEC",   "f");
  l->AddEntry(h_ccdis, "CC DIS",   "f");
  l->AddEntry(h_ccoth, "CC Other", "f");
  l->AddEntry(h_nc,    "NC",       "f");
  l->SetLineWidth(0);
  l->SetTextAlign(22);
  l->SetTextFont(132);
  
  gROOT->ForceStyle();
  c->UseCurrentStyle();
  hsT->Draw("hist");
  l->Draw("same");
  const char *x_label = "Proton multiplicity";
  const char *y_label = "Normalised event rate";
  hsT->GetXaxis()->SetTitle(x_label);   
  hsT->GetYaxis()->SetTitle(y_label);   
  
  c->SaveAs((plots+"true_multiplicity_pre_fsi.png").c_str());
  c->SaveAs((plots+"true_multiplicity_pre_fsi.root").c_str());

  l->Clear();
  c->Clear();
  
  l->AddEntry(hr_ccqe,  "CC QE",    "f");
  l->AddEntry(hr_ccres, "CC RES",   "f");
  l->AddEntry(hr_ccmec, "CC MEC",   "f");
  l->AddEntry(hr_ccdis, "CC DIS",   "f");
  l->AddEntry(hr_ccoth, "CC Other", "f");
  l->AddEntry(hr_nc,    "NC",       "f");
  l->SetLineWidth(0);
  l->SetTextAlign(22);
  l->SetTextFont(132);
  
  gROOT->ForceStyle();
  c->UseCurrentStyle();
  hsR->Draw("hist");
  l->Draw("same");
  hsR->GetXaxis()->SetTitle(x_label);   
  hsR->GetYaxis()->SetTitle(y_label);   
  
  c->SaveAs((plots+"reco_multiplicity_pre_fsi.png").c_str());
  c->SaveAs((plots+"reco_multiplicity_pre_fsi.root").c_str());

  l->Clear();
  c->Clear();
  
  
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  l->AddEntry( h_signal_energy,     " True, reconstructed proton", "l" );
  l->AddEntry( h_background_energy, " True, misidentified proton", "l" );
  l->AddEntry( h_missed_energy,     " True, missed proton", "l" );
  
  h_signal_energy->SetLineColor(2);
  h_signal_energy->SetStats(kFALSE);
  h_signal_energy->Scale(1/double(signal_energy.size()));
  h_background_energy->SetLineColor(4);
  h_background_energy->SetStats(kFALSE);
  h_background_energy->GetXaxis()->SetTitle("Proton Kinetic Energy [GeV]");
  h_background_energy->Scale(1/double(background_energy.size()));
  h_missed_energy->SetLineColor(8);
  h_missed_energy->SetStats(kFALSE);
  h_missed_energy->Scale(1/double(missed_energy.size()));

  float max_y_e = std::max(h_background_energy->GetMaximum(), std::max(h_signal_energy->GetMaximum(),h_missed_energy->GetMaximum()));
  
  h_background_energy->GetYaxis()->SetRangeUser(0,max_y_e*1.1);
  h_background_energy->Draw("hist");
  h_missed_energy->Draw("hist same");
  h_signal_energy->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"proton_kinetic_energy.root").c_str());
  c->SaveAs((plots+"proton_kinetic_energy.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_length,     " True, reconstructed proton", "l" );
  l->AddEntry( h_background_length, " True, misidentified proton", "l" );
  l->AddEntry( h_missed_length,     " True, missed proton", "l" );
    
  h_signal_length->SetLineColor(2);
  h_signal_length->SetStats(kFALSE);
  h_signal_length->Scale(1/double(signal_length.size()));
  h_background_length->SetLineColor(4);
  h_background_length->SetStats(kFALSE);
  h_background_length->GetXaxis()->SetTitle("Proton Length [cm]");
  h_background_length->Scale(1/double(background_length.size()));
  h_missed_length->SetLineColor(8);
  h_missed_length->SetStats(kFALSE);
  h_missed_length->Scale(1/double(missed_length.size()));

  float max_y_l = std::max(h_background_length->GetMaximum(), std::max(h_signal_length->GetMaximum(),h_missed_length->GetMaximum()));
  h_background_length->GetYaxis()->SetRangeUser(0,max_y_l*1.1);
  h_background_length->Draw("hist");
  h_missed_length->Draw("hist same");
  h_signal_length->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"proton_length.root").c_str());
  c->SaveAs((plots+"proton_length.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_hits,     " True, reconstructed proton", "l" );
  l->AddEntry( h_background_hits, " True, misidentified proton", "l" );
  l->AddEntry( h_missed_hits,     " True, missed proton", "l" );
  
  h_signal_hits->SetLineColor(2);
  h_signal_hits->SetStats(kFALSE);
  h_signal_hits->Scale(1/double(signal_hits.size()));
  h_background_hits->SetLineColor(4);
  h_background_hits->SetStats(kFALSE);
  h_background_hits->GetXaxis()->SetTitle("Proton Hits [GeV]");
  h_background_hits->Scale(1/double(background_hits.size()));
  h_missed_hits->SetLineColor(8);
  h_missed_hits->SetStats(kFALSE);
  h_missed_hits->Scale(1/double(missed_hits.size()));

  float max_y_h = std::max(h_background_hits->GetMaximum(), std::max(h_signal_hits->GetMaximum(),h_missed_hits->GetMaximum()));
  
  h_background_hits->GetYaxis()->SetRangeUser(0,max_y_h*1.1);
  h_background_hits->Draw("hist");
  h_missed_hits->Draw("hist same");
  h_signal_hits->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"proton_hits.root").c_str());
  c->SaveAs((plots+"proton_hits.png").c_str());
  c->Clear();

  std::cout << " Total number of reconstructed protons :                " << total << std::endl;
  std::cout << " Missed protons :                                       " << missed << std::endl;
  std::cout << " Missed protons with less than 25 hits :                " << missed_below_25 << std::endl;
  std::cout << " Missed protons with less than 10 hits  :               " << missed_below_10 << std::endl;
  std::cout << " Percentage of missed protons with less than 25 hits :  " << missed_below_25 / double(missed) << std::endl;
  std::cout << " Percentage of missed protons with less than 10 hits  : " << missed_below_10 / double(missed) << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end)       << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<int> &exceptions) {
  double total_pot = 0;
  std::vector<int>::iterator it;
  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    it = std::find(exceptions.begin(), exceptions.end(),i);
    if(it != exceptions.end()) continue;
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/pnfs/sbnd/persistent/users/rsjones/sbnd_selection_150620/selection/"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    double temp_pot = 0.;
    EventSelectionTool::LoadEventList(file_name, events, i, temp_pot);
    EventSelectionTool::GetTimeLeft(start_time,total_files,i);
    total_pot += temp_pot;
  }
  std::cout << std::endl;
  pot = total_pot;
} // LoadAllEvents
