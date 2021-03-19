#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"

using namespace cppsecrets;
using namespace selection;

int MainTest(const char *config){
  
  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
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

  // Setup the histograms
  TH2D *h_dedx_resrange_mu = new TH2D("h_dedx_resrange_mu", "",250,0,20,250,0,20);
  TH2D *h_dedx_resrange_pr = new TH2D("h_dedx_resrange_pr", "",250,0,20,250,0,20);
  TH2D *h_dedx_resrange_pi = new TH2D("h_dedx_resrange_pi", "",250,0,20,250,0,20);

  TH2D *h_dedx_resrange_mu_mc = new TH2D("h_dedx_resrange_mu_mc", "",250,0,20,250,0,20);
  TH2D *h_dedx_resrange_pr_mc = new TH2D("h_dedx_resrange_pr_mc", "",250,0,20,250,0,20);
  TH2D *h_dedx_resrange_pi_mc = new TH2D("h_dedx_resrange_pi_mc", "",250,0,20,250,0,20);

  // Counters
  unsigned int fNMuons   = 0;
  unsigned int fNPions   = 0;
  unsigned int fNProtons = 0;

  unsigned int fNMuonsMC   = 0;
  unsigned int fNPionsMC   = 0;
  unsigned int fNProtonsMC = 0;

  // First, ensure all tracks are contained
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    selection::LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active, true);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
    // Now loop over the events
    for(const Event &e : events){
      // First, fill the truth-level dEdx
      ParticleList mcList   = e.GetMCParticleList();


      ParticleList recoList = e.GetRecoParticleList();
      // Now loop over the reconstructed particles
      for(const Particle &p: recoList){
        // Make sure particle is a reconstructed track and is fully contained
        if(!p.GetFromRecoTrack()) continue;
        if(!p.GetTrackContained()) continue;
        if(p.GetMomentum().Mag() < 0.2) continue; // 200 MeV threshold

        // Now loop over the dEdx and residual ranges and fill
        if(p.GetdEdx().size() == 0) continue;
        assert(p.GetdEdx().size() == p.GetResidualRange().size());
        for(unsigned int e = 0; e < p.GetdEdx().size(); ++e){
        
          double dEdx     = p.GetdEdx().at(e); // MeV
          double resRange = p.GetResidualRange().at(e); // cm
          if(abs(p.GetPdgCode()) == 13){
            h_dedx_resrange_mu->Fill(resRange, dEdx);
            fNMuons++;
          }
          if(abs(p.GetPdgCode()) == 211){
            h_dedx_resrange_pi->Fill(resRange, dEdx);
            fNPions++;
          }
          if(abs(p.GetPdgCode()) == 2212){
            h_dedx_resrange_pr->Fill(resRange, dEdx);
            fNProtons++;
          }
        }
      } // Reco particles
    } // Events
  } // Files


  // Print counters
  std::cout << "Total number of muons: " << fNMuons << ", protons: " << fNProtons << ", pions: " << fNPions << std::endl;

  // Format the histograms for the overlay
  h_dedx_resrange_mu  ->SetMarkerColor(kViolet+5);
  h_dedx_resrange_pi  ->SetMarkerColor(kOrange+5);
  h_dedx_resrange_pr  ->SetMarkerColor(kTeal-6);
  h_dedx_resrange_mu  ->SetMarkerStyle(7);
  h_dedx_resrange_pi  ->SetMarkerStyle(7);
  h_dedx_resrange_pr  ->SetMarkerStyle(7);

  h_dedx_resrange_mu->GetXaxis()->SetTitle("Residual Range [cm]");
  h_dedx_resrange_mu->GetYaxis()->SetTitle("dE/dx [MeV / cm]");
  h_dedx_resrange_mu->GetXaxis()->SetTitleSize(0.06);
  h_dedx_resrange_mu->GetYaxis()->SetTitleSize(0.06);
  h_dedx_resrange_mu->GetXaxis()->SetLabelSize(0.05);
  h_dedx_resrange_mu->GetYaxis()->SetLabelSize(0.05);
  h_dedx_resrange_mu->GetXaxis()->SetTitleFont(132);
  h_dedx_resrange_mu->GetYaxis()->SetTitleFont(132);
  h_dedx_resrange_mu->GetXaxis()->SetLabelFont(132);
  h_dedx_resrange_mu->GetYaxis()->SetLabelFont(132);
  h_dedx_resrange_mu->SetStats(0);

  // Make the canvas
  // Open ROOT file for writing
  TFile fFile((plots_location+"testdEdx.root").c_str(),"RECREATE");

  TCanvas *fC = new TCanvas("c","",900,900);
  fC->SetLeftMargin(0.13);
  fC->SetBottomMargin(0.13);
  fC->SetTopMargin(0.08);
  fC->SetRightMargin(0.025);

  TLegend *fL = new TLegend(0.14,0.93,0.98,0.99);
  fL->SetTextFont(132);
  fL->SetTextSize(0.05);
  fL->SetBorderSize(0);
  fL->SetNColumns(3);

  fL->AddEntry(h_dedx_resrange_mu, "#color[885]{Muons}",   "p");
  fL->AddEntry(h_dedx_resrange_pi, "#color[805]{Pions}",   "p");
  fL->AddEntry(h_dedx_resrange_pr, "#color[834]{Protons}", "p");

  h_dedx_resrange_mu->Draw();
  h_dedx_resrange_pi->Draw("same");
  h_dedx_resrange_pr->Draw("same");
  fL->Draw("same");

  h_dedx_resrange_mu  ->Write();
  h_dedx_resrange_pi  ->Write();
  h_dedx_resrange_pr  ->Write();

  fC->Write();
  fFile.Close();

  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()


