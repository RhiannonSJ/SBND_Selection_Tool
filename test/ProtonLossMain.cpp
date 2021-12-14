#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"

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

using namespace cppsecrets;
using namespace selection;

int MainTest(const char *config){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the class ConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get variables from config
  std::string input_location  = "";
  std::string input_filename  = "";
  std::string exceptions_file = "";
  std::string stats_location  = "";
  std::string plots_location  = "";
  unsigned int total_files = 0;
  unsigned int detector = 0; // 0 = sbnd, 1 = uboone, 2 = icarus
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputFileLocation",input_location);
  p->getValue("InputFileName",    input_filename);
  p->getValue("ExceptionsFile",   exceptions_file);
  p->getValue("StatFileLocation", stats_location); 
  p->getValue("PlotFileLocation", plots_location); 
  p->getValue("TotalFiles",       total_files);
  p->getValue("Detector",         detector);
  p->getValue("MinXFid",          minx_fid);
  p->getValue("MinYFid",          miny_fid);
  p->getValue("MinZFid",          minz_fid);
  p->getValue("MaxXFid",          maxx_fid);
  p->getValue("MaxYFid",          maxy_fid);
  p->getValue("MaxZFid",          maxz_fid);
  p->getValue("MinXAV",           minx_av);
  p->getValue("MinYAV",           miny_av);
  p->getValue("MinZAV",           minz_av);
  p->getValue("MaxXAV",           maxx_av);
  p->getValue("MaxYAV",           maxy_av);
  p->getValue("MaxZAV",           maxz_av);

  //------------------------------------------------------------------------------------------
  //                                    Initialise
  //------------------------------------------------------------------------------------------

  // Get the active and fiducial geometry objects
  Geometry fiducial(minx_fid,miny_fid,minz_fid,maxx_fid,maxy_fid,maxz_fid,true);
  Geometry active(minx_av,miny_av,minz_av,maxx_av,maxy_av,maxz_av,false);
  PlaneList planes = active.GetExternalPlaneList();

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;

  int start = static_cast<int>(time(NULL));
  double pot = 0.; 

  std::vector<int> exceptions;
  FillExceptions(exceptions_file.c_str(),exceptions);

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

  // First, ensure all tracks are contained
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
    // Now loop over the events
    for(const Event &e : events){
      bool cc_inclusive_passed = GeneralAnalysisHelper::PassedCCInclusive(e,detector);

      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial() || 
         !GeneralAnalysisHelper::MaxOneLongEscapingTrack(e) || 
         !GeneralAnalysisHelper::MinOneRecoTrack(e) ||
         !cc_inclusive_passed) continue;

      ParticleList reco_particles = e.GetRecoParticleList(); 
      ParticleList true_particles = e.GetMCParticleList();

      for(Particle &p_reco : reco_particles){
        // If the particle is a reconstructed track (not pi0)
        if(p_reco.GetFromRecoTrack() && GeneralAnalysisHelper::ParticleHasAMatch(e, p_reco) >= 0){
          Particle p_match = GeneralAnalysisHelper::GetBestMCParticle(e,p_reco);
          if(p_match.GetPdgCode() == 2212 && p_match.GetNumberOfHits() >= 5 && p_match.GetKineticEnergy() > 0.021){
            if(p_reco.GetPdgCode() == 2212 && p_reco.GetNumberOfHits() >= 5 && p_reco.GetKineticEnergy() > 0.021){
              // True and reconstructed proton: signal
              signal_energy.push_back(p_match.GetKineticEnergy());
              signal_length.push_back(p_match.GetLength());
              signal_hits.push_back(p_match.GetNumberOfHits());
            }
            else if(p_reco.GetPdgCode() != 2212 && p_reco.GetNumberOfHits() >= 5 && p_reco.GetKineticEnergy() > 0.021){
              // True and mis-indentified proton
              background_energy.push_back(p_match.GetKineticEnergy());
              background_length.push_back(p_match.GetLength());
              background_hits.push_back(p_match.GetNumberOfHits());
            }
          }
        }
      }
      for(Particle &p_true : true_particles){
        if(p_true.GetPdgCode() == 2212 && !GeneralAnalysisHelper::HasBeenReconstructed(e, p_true) && p_true.GetNumberOfHits() >= 5 && p_true.GetKineticEnergy() > 0.021){
          // True proton, not reconstructed
          missed_energy.push_back(p_true.GetKineticEnergy());
          missed_length.push_back(p_true.GetLength());
          missed_hits.push_back(p_true.GetNumberOfHits());
          missed++;
          if(p_true.GetNumberOfHits() < 25) missed_below_25++;
          if(p_true.GetNumberOfHits() < 10) missed_below_10++;
        }
        else if(p_true.GetPdgCode() == 2212 && GeneralAnalysisHelper::HasBeenReconstructed(e, p_true) && p_true.GetNumberOfHits() >= 5 && p_true.GetKineticEnergy() > 0.021)
          total++;
      }
    }
  }
  TH1D *h_signal_energy = new TH1D("h_signal_energy","Correctly reconstructed proton kinetic energies",80,0,0.6);
  TH1D *h_signal_length = new TH1D("h_signal_length","Correctly reconstructed proton lengths",80,0,15);
  TH1D *h_signal_hits   = new TH1D("h_signal_hits",  "Correctly reconstructed proton hits",80,0,160);

  TH1D *h_missed_energy = new TH1D("h_missed_energy","Missed proton kinetic energies",80,0,0.6);
  TH1D *h_missed_length = new TH1D("h_missed_length","Missed proton lengths",80,0,15);
  TH1D *h_missed_hits   = new TH1D("h_missed_hits",  "Missed proton hits",80,0,160);

  TH1D *h_background_energy = new TH1D("h_background_energy","Mis-identified proton kinetic energies",80,0,0.6);
  TH1D *h_background_length = new TH1D("h_background_length","Mis-identified proton lengths",80,0,15);
  TH1D *h_background_hits   = new TH1D("h_background_hits",  "Mis-identified proton hits",80,0,160);


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

  TCanvas *c = new TCanvas("c","",900,900);
  c->SetLeftMargin  (0.138796 );
  c->SetRightMargin (0.0334448);
  c->SetBottomMargin(0.132404 );
  c->SetTopMargin   (0.0734);

  TLegend *l = new TLegend(0.15, 0.94, 0.97, 0.99);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.045);
  l->SetNColumns(3);
  l->SetMargin(0.2);
  l->SetTextFont(132);
  
  // Style
  std::vector<TH1D*> h_signal     = {h_signal_energy,h_signal_length,h_signal_hits};
  std::vector<TH1D*> h_background = {h_background_energy,h_background_length,h_background_hits};
  std::vector<TH1D*> h_missed     = {h_missed_energy,h_missed_length,h_missed_hits};

  SetHistogramStyle(h_signal,     3001, 1, 905, 905, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_background, 3001, 1, 867, 867, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_missed,     3001, 1, 835, 835, 2, 132, 1, 1.2, false);

  l->AddEntry( h_signal_energy,     " Reconstructed", "f" );
  l->AddEntry( h_background_energy, " Misidentified", "f" );
  l->AddEntry( h_missed_energy,     " Missed", "f" );
  
  h_signal_energy->Scale(1/double(h_signal_energy->Integral()));
  h_background_energy->Scale(1/double(h_background_energy->Integral()));
  h_missed_energy->Scale(1/double(h_missed_energy->Integral()));
  
  h_background_energy->GetYaxis()->SetTitle("Normalised event rate");
  h_background_energy->GetXaxis()->SetTitle("Proton Kinetic Energy [GeV]");
  
  float max_y_e = std::max(h_background_energy->GetMaximum(), std::max(h_signal_energy->GetMaximum(),h_missed_energy->GetMaximum()));
  
  h_background_energy->GetYaxis()->SetRangeUser(0,max_y_e*1.1);
  h_background_energy->SetTitle("");
  h_background_energy->Draw("hist");
  h_missed_energy->Draw("hist same");
  h_signal_energy->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"ccpassed/proton_kinetic_energy.root").c_str());
  c->SaveAs((plots_location+"ccpassed/proton_kinetic_energy.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_length,     " Reconstructed", "f" );
  l->AddEntry( h_background_length, " Misidentified", "f" );
  l->AddEntry( h_missed_length,     " Missed", "f" );
    
  h_signal_length->Scale(1/double(h_signal_length->Integral()));
  h_background_length->Scale(1/double(h_background_length->Integral()));
  h_missed_length->Scale(1/double(h_missed_length->Integral()));
  
  h_background_length->GetYaxis()->SetTitle("Normalised event rate");
  h_background_length->GetXaxis()->SetTitle("Proton Length [cm]");

  float max_y_l = std::max(h_background_length->GetMaximum(), std::max(h_signal_length->GetMaximum(),h_missed_length->GetMaximum()));
  h_background_length->GetYaxis()->SetRangeUser(0,max_y_l*1.1);
  h_background_length->SetTitle("");
  h_background_length->Draw("hist");
  h_missed_length->Draw("hist same");
  h_signal_length->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"ccpassed/proton_length.root").c_str());
  c->SaveAs((plots_location+"ccpassed/proton_length.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_hits,     " Reconstructed", "f" );
  l->AddEntry( h_background_hits, " Misidentified", "f" );
  l->AddEntry( h_missed_hits,     " Missed", "f" );
  
  h_signal_hits->Scale(1/double(h_signal_hits->Integral()));
  h_background_hits->Scale(1/double(h_background_hits->Integral()));
  h_missed_hits->Scale(1/double(h_missed_hits->Integral()));
  
  h_background_hits->GetYaxis()->SetTitle("Normalised event rate");
  h_background_hits->GetXaxis()->SetTitle("Proton Hits [GeV]");

  float max_y_h = std::max(h_background_hits->GetMaximum(), std::max(h_signal_hits->GetMaximum(),h_missed_hits->GetMaximum()));
  
  h_background_hits->GetYaxis()->SetRangeUser(0,max_y_h*1.1);
  h_background_hits->SetTitle("");
  h_background_hits->Draw("hist");
  h_missed_hits->Draw("hist same");
  h_signal_hits->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"ccpassed/proton_hits.root").c_str());
  c->SaveAs((plots_location+"ccpassed/proton_hits.png").c_str());
  c->Clear();

  // Stats file
  std::ofstream stats;
  stats.open(stats_location+"proton_loss.txt");

  stats << " Total number of reconstructed protons :                " << total << std::endl;
  stats << " Missed protons :                                       " << missed << std::endl;
  stats << " Missed protons with less than 25 hits :                " << missed_below_25 << std::endl;
  stats << " Missed protons with less than 10 hits  :               " << missed_below_10 << std::endl;
  stats << " Percentage of missed protons with less than 25 hits :  " << missed_below_25 / double(missed) << std::endl;
  stats << " Percentage of missed protons with less than 10 hits  : " << missed_below_10 / double(missed) << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end)       << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
