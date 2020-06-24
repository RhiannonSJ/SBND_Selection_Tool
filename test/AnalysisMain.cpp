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

  TopologyMap nc_map      = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc_map      = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_map   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map   = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc0pi2p_map = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap nc0pi_map   = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_map   = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nue_map     = GeneralAnalysisHelper::GetNuETopologyMap();

  std::vector< TopologyMap > maps({cc_map, cc0pi_map, cc1pi_map, cc0pi2p_map, nc_map, nc0pi_map, nc1pi_map, nue_map});  

  // COUNTERS
  unsigned int total_true_ccinc = 0;
  unsigned int total_true_ncinc = 0;
  
  unsigned int reco_fid_true_ccinc = 0;
  unsigned int reco_fid_true_ncinc = 0;
  
  unsigned int reco_true_fid_true_ccinc = 0;
  unsigned int reco_true_fid_true_ncinc = 0;
  
  unsigned int reco_true_fid_max_1_escapes_true_ccinc = 0;
  unsigned int reco_true_fid_max_1_escapes_true_ncinc = 0;
  
  // HISTOGRAMS
  TH1D *h_longest_diff_mu = new TH1D("h_longest_diff_mu", "Percentage length difference between 2 longest tracks",40,0,1);
  TH1D *h_longest_diff_pr = new TH1D("h_longest_diff_pr", "Percentage length difference between 2 longest tracks",40,0,1);
  TH1D *h_longest_diff_pi = new TH1D("h_longest_diff_pi", "Percentage length difference between 2 longest tracks",40,0,1);

  TH1D *h_longest_length_mu = new TH1D("h_longest_length_mu", "Longest track length",40,0,700);
  TH1D *h_longest_length_pr = new TH1D("h_longest_length_pr", "Longest track length",40,0,700);
  TH1D *h_longest_length_pi = new TH1D("h_longest_length_pi", "Longest track length",40,0,700);

  TH1D *h_length_mu = new TH1D("h_length_mu", "Track length",40,0,500);
  TH1D *h_length_pr = new TH1D("h_length_pr", "Track length",40,0,500);
  TH1D *h_length_pi = new TH1D("h_length_pi", "Track length",40,0,500);

  TH1D *h_chi2_pr_mu = new TH1D("h_chi2_pr_mu", "#Chi^{2}_{proton}",40,0,120);
  TH1D *h_chi2_pr_pr = new TH1D("h_chi2_pr_pr", "#Chi^{2}_{proton}",40,0,120);
  TH1D *h_chi2_pr_pi = new TH1D("h_chi2_pr_pi", "#Chi^{2}_{proton}",40,0,120);

  TH1D *h_chi2_mu_mu = new TH1D("h_chi2_mu_mu", "#Chi^{2}_{muon}",40,0,70);
  TH1D *h_chi2_mu_pr = new TH1D("h_chi2_mu_pr", "#Chi^{2}_{muon}",40,0,70);
  TH1D *h_chi2_mu_pi = new TH1D("h_chi2_mu_pi", "#Chi^{2}_{muon}",40,0,70);
  
  TH1D *h_chi2_ratio_mu = new TH1D("h_chi2_ratio_mu", "#Chi^{2}_{muon}/#Chi^{2}_{proton}",40,0,0.5);
  TH1D *h_chi2_ratio_pr = new TH1D("h_chi2_ratio_pr", "#Chi^{2}_{muon}/#Chi^{2}_{proton}",40,0,0.5);
  TH1D *h_chi2_ratio_pi = new TH1D("h_chi2_ratio_pi", "#Chi^{2}_{muon}/#Chi^{2}_{proton}",40,0,0.5);

  TH2D *h_chi2_mu_chi2_pr_mu = new TH2D("h_chi2_mu_chi2_pr_mu", "#Chi^{2}_{muon} vs #Chi^{2}_{proton}",40,0,100,40,0,100);
  TH2D *h_chi2_mu_chi2_pr_pr = new TH2D("h_chi2_mu_chi2_pr_pr", "#Chi^{2}_{muon} vs #Chi^{2}_{proton}",40,0,100,40,0,100);
  TH2D *h_chi2_mu_chi2_pr_pi = new TH2D("h_chi2_mu_chi2_pr_pi", "#Chi^{2}_{muon} vs #Chi^{2}_{proton}",40,0,100,40,0,100);

  std::vector<TH1D*> h_muon   = {h_longest_diff_mu, h_longest_length_mu, h_length_mu, h_chi2_pr_mu, h_chi2_mu_mu, h_chi2_ratio_mu};
  std::vector<TH1D*> h_pion   = {h_longest_diff_pi, h_longest_length_pi, h_length_pi, h_chi2_pr_pi, h_chi2_mu_pi, h_chi2_ratio_pi};
  std::vector<TH1D*> h_proton = {h_longest_diff_pr, h_longest_length_pr, h_length_pr, h_chi2_pr_pr, h_chi2_mu_pr, h_chi2_ratio_pr};
  std::vector<TH2D*> h_2d     = {h_chi2_mu_chi2_pr_mu, h_chi2_mu_chi2_pr_pr, h_chi2_mu_chi2_pr_pi};
  
  // First, ensure all tracks are contained
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
    for(const Event &e : events){
      // Fill the counters
      if(e.CheckMCTopology(cc_map)) {
        total_true_ccinc++;
        if(e.IsRecoFiducial()){
          reco_fid_true_ccinc++;
          if(e.IsTrueFiducial()){
            reco_true_fid_true_ccinc++;
            if(GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
              reco_true_fid_max_1_escapes_true_ccinc++;
            }
          }
        }
      }
      else if(e.CheckMCTopology(nc_map)) {
        total_true_ncinc++;
        if(e.IsRecoFiducial()){
          reco_fid_true_ncinc++;
          if(e.IsTrueFiducial()){
            reco_true_fid_true_ncinc++;
            if(GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
              reco_true_fid_max_1_escapes_true_ncinc++;
            }
          }
        }
      }
      // Define and fill the histograms
      // Longest and second longest track lengths
      double longest = -std::numeric_limits<double>::max();
      double second  = -std::numeric_limits<double>::max();
      int longest_id = -1;
      int second_id  = -1;
      bool min_2_tracks = false;
      for(const Particle &p : e.GetMCParticleList()){
        if(p.GetLength() > longest){
          longest_id = p.GetMCId();
          longest = p.GetLength();
        }
      }
      for(const Particle &p : e.GetMCParticleList()){
        if(p.GetMCId() != longest_id && p.GetLength() > second){
          second_id = p.GetMCId();
          second = p.GetLength();
        }
      }

      // Get the fractional difference between the longest and second longest track lengths
      double diff = -999.;
      if(longest_id != -1 && second_id != -1)
        diff = (longest - second)/longest;

      // Now fill the appropriate histograms
      for(const Particle &p : e.GetMCParticleList()){
        // Muon
        if(abs(p.GetPdgCode()) == 13){
          h_length_mu->Fill(p.GetLength());
          if(p.GetMCId() == longest_id){
            h_longest_diff_mu->Fill(diff);
            h_longest_length_mu->Fill(longest);
          }
        }
        // Pion
        else if(abs(p.GetPdgCode()) == 211){
          h_length_pi->Fill(p.GetLength());
          if(p.GetMCId() == longest_id){
            h_longest_diff_pi->Fill(diff);
            h_longest_length_pi->Fill(longest);
          }
        }
        // Proton
        else if(abs(p.GetPdgCode()) == abs(2212)){
          h_length_pr->Fill(p.GetLength());
          if(p.GetMCId() == longest_id){
            h_longest_diff_pr->Fill(diff);
            h_longest_length_pr->Fill(longest);
          }
        }
      }

      // Need reconstructed particles for chi2 variables
      // Though take the PDGCode from truth
      for(const Particle &p : e.GetRecoParticleList()){
        if(!p.GetFromRecoTrack()) continue;
        if(!p.GetHasCalorimetry()) continue;
        if(GeneralAnalysisHelper::ParticleHasAMatch(e,p) == -1) continue;
        Particle mcp = GeneralAnalysisHelper::GetBestMCParticle(e,p);
        double chi2_ratio = p.GetChi2Mu()/p.GetChi2P();
        // Muon
        if(abs(mcp.GetPdgCode()) == 13){
          h_chi2_pr_mu->Fill(p.GetChi2P());
          h_chi2_mu_mu->Fill(p.GetChi2Mu());
          h_chi2_ratio_mu->Fill(chi2_ratio);
          h_chi2_mu_chi2_pr_mu->Fill(p.GetChi2Mu(),p.GetChi2P());
        }
        // Pion
        if(abs(mcp.GetPdgCode()) == 211){
          h_chi2_pr_pi->Fill(p.GetChi2P());
          h_chi2_mu_pi->Fill(p.GetChi2Mu());
          h_chi2_ratio_pi->Fill(chi2_ratio);
          h_chi2_mu_chi2_pr_pi->Fill(p.GetChi2Mu(),p.GetChi2P());
        }
        // Proton
        if(abs(mcp.GetPdgCode()) == 2212){
          h_chi2_pr_pr->Fill(p.GetChi2P());
          h_chi2_mu_pr->Fill(p.GetChi2Mu());
          h_chi2_ratio_pr->Fill(chi2_ratio);
          h_chi2_mu_chi2_pr_pr->Fill(p.GetChi2Mu(),p.GetChi2P());
        }
      }
    }
  }

  // Write histograms
  // Set the style of the muon, pion and proton histograms
  SetHistogramStyle(h_muon,   3001, 1, 905, 905, 2, 132, 1.2, 1.5, true);
  SetHistogramStyle(h_pion,   3001, 1, 801, 801, 2, 132, 1.2, 1.5, true);
  SetHistogramStyle(h_proton, 3001, 1, 867, 867, 2, 132, 1.2, 1.5, true);

  // Test
  TCanvas *c = new TCanvas("c","",900,900);
  c->SetTopMargin(0.0320233);
  c->SetBottomMargin(0.0946143);
  c->SetLeftMargin(0.113181);
  c->SetRightMargin(0.0320233);

  TLegend * l = new TLegend(0.355,0.681,0.977,0.929);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextFont(132);

  // Longest diffs
  l->AddEntry(h_longest_diff_mu,"Muons","f");
  l->AddEntry(h_longest_diff_pi,"Pions","f");
  l->AddEntry(h_longest_diff_pr,"Protons","f");

  double max_diff = 1.1 * std::max(h_longest_diff_mu->GetMaximum(), std::max(h_longest_diff_pi->GetMaximum(),h_longest_diff_pr->GetMaximum()));
  h_longest_diff_mu->GetYaxis()->SetRangeUser(0.,max_diff);
  h_longest_diff_mu->SetTitle("");
  h_longest_diff_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_longest_diff_mu->GetXaxis()->SetTitle("Fractional length difference between two longest tracks [cm]");
  h_longest_diff_mu->Draw("hist");
  h_longest_diff_pi->Draw("hist same");
  h_longest_diff_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"longest_length_differences.root").c_str());
  c->SaveAs((plots_location+"longest_length_differences.png").c_str());
  l->Clear();
  c->Clear();

  // Longest lengths
  l->AddEntry(h_longest_length_mu,"Muons","f");
  l->AddEntry(h_longest_length_pi,"Pions","f");
  l->AddEntry(h_longest_length_pr,"Protons","f");

  double max_longest = 1.1 * std::max(h_longest_length_mu->GetMaximum(), std::max(h_longest_length_pi->GetMaximum(),h_longest_length_pr->GetMaximum()));
  h_longest_length_mu->GetYaxis()->SetRangeUser(0.,max_longest);
  h_longest_length_mu->SetTitle("");
  h_longest_length_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_longest_length_mu->GetXaxis()->SetTitle("Longest track length [cm]");
  h_longest_length_mu->Draw("hist");
  h_longest_length_pi->Draw("hist same");
  h_longest_length_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"longest_lengths.root").c_str());
  c->SaveAs((plots_location+"longest_lengths.png").c_str());
  l->Clear();
  c->Clear();

  // Track lengths
  l->AddEntry(h_length_mu,"Muons","f");
  l->AddEntry(h_length_pi,"Pions","f");
  l->AddEntry(h_length_pr,"Protons","f");

  double max_length = 1.1 * std::max(h_length_mu->GetMaximum(), std::max(h_length_pi->GetMaximum(),h_length_pr->GetMaximum()));
  h_length_mu->GetYaxis()->SetRangeUser(0.,max_length);
  h_length_mu->SetTitle("");
  h_length_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_length_mu->GetXaxis()->SetTitle("Track length [cm]");
  h_length_mu->Draw("hist");
  h_length_pi->Draw("hist same");
  h_length_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_lengths.root").c_str());
  c->SaveAs((plots_location+"track_lengths.png").c_str());
  l->Clear();
  c->Clear();

  // Chi2 muon
  l->AddEntry(h_chi2_mu_mu,"Muons","f");
  l->AddEntry(h_chi2_mu_pi,"Pions","f");
  l->AddEntry(h_chi2_mu_pr,"Protons","f");

  double max_chi2_mu = 1.1 * std::max(h_chi2_mu_mu->GetMaximum(), std::max(h_chi2_mu_pi->GetMaximum(),h_chi2_mu_pr->GetMaximum()));
  h_chi2_mu_mu->GetYaxis()->SetRangeUser(0.,max_chi2_mu);
  h_chi2_mu_mu->SetTitle("");
  h_chi2_mu_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_chi2_mu_mu->GetXaxis()->SetTitle("#Chi^{2}_{#mu}");
  h_chi2_mu_mu->Draw("hist");
  h_chi2_mu_pi->Draw("hist same");
  h_chi2_mu_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_mus.root").c_str());
  c->SaveAs((plots_location+"track_chi2_mus.png").c_str());
  l->Clear();
  c->Clear();

  // Chi2 proton
  l->AddEntry(h_chi2_pr_mu,"Muons","f");
  l->AddEntry(h_chi2_pr_pi,"Pions","f");
  l->AddEntry(h_chi2_pr_pr,"Protons","f");

  double max_chi2_pr = 1.1 * std::max(h_chi2_pr_mu->GetMaximum(), std::max(h_chi2_pr_pi->GetMaximum(),h_chi2_pr_pr->GetMaximum()));
  h_chi2_pr_mu->GetYaxis()->SetRangeUser(0.,max_chi2_pr);
  h_chi2_pr_mu->SetTitle("");
  h_chi2_pr_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_chi2_pr_mu->GetXaxis()->SetTitle("#Chi^{2}_{proton}");
  h_chi2_pr_mu->Draw("hist");
  h_chi2_pr_pi->Draw("hist same");
  h_chi2_pr_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_prs.root").c_str());
  c->SaveAs((plots_location+"track_chi2_prs.png").c_str());
  l->Clear();
  c->Clear();

  // Chi2 muon/proton
  l->AddEntry(h_chi2_ratio_mu,"Muons","f");
  l->AddEntry(h_chi2_ratio_pi,"Pions","f");
  l->AddEntry(h_chi2_ratio_pr,"Protons","f");

  double max_chi2_ratio = 1.1 * std::max(h_chi2_ratio_mu->GetMaximum(), std::max(h_chi2_ratio_pi->GetMaximum(),h_chi2_ratio_pr->GetMaximum()));
  h_chi2_ratio_mu->GetYaxis()->SetRangeUser(0.,max_chi2_ratio);
  h_chi2_ratio_mu->SetTitle("");
  h_chi2_ratio_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_chi2_ratio_mu->GetXaxis()->SetTitle("#Chi^{2}_{#mu} / #Chi^{2}_{proton}");
  h_chi2_ratio_mu->Draw("hist");
  h_chi2_ratio_pi->Draw("hist same");
  h_chi2_ratio_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_ratios.root").c_str());
  c->SaveAs((plots_location+"track_chi2_ratios.png").c_str());
  l->Clear();
  c->Clear();

  // Write statistics
  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"analysis_cuts.txt");

  file << "=====================================================================" << std::endl;
  file << " Total POT used to generate this sample: " << pot << std::endl;
  file << std::setw(20) << "Cut/True top."    << "||";
  file <<  std::setw(15) << " CC Inc. ";
  file << std::setw(15) << " NC Inc. ";
  file << std::setw(15) << " CCInc. purity" << std::endl;

  file << std::setw(20) << " Total "          << "||";
  file << std::setw(15) << total_true_ccinc;
  file << std::setw(15) << total_true_ncinc;
  file << std::setw(15) << total_true_ccinc / static_cast<double>(total_true_ccinc+total_true_ncinc) << std::endl;

  file << std::setw(20) << " Reco fiducial "  << "||";
  file << std::setw(15) << reco_fid_true_ccinc;
  file << std::setw(15) << reco_fid_true_ncinc;
  file << std::setw(15) << reco_fid_true_ccinc/static_cast<double>(reco_fid_true_ccinc+reco_fid_true_ncinc) << std::endl; 

  file << std::setw(20) << " True fiducial "  << "||";
  file << std::setw(15) << reco_true_fid_true_ccinc;
  file << std::setw(15) << reco_true_fid_true_ncinc;
  file << std::setw(15) << reco_true_fid_true_ccinc/static_cast<double>(reco_true_fid_true_ccinc+reco_true_fid_true_ncinc) << std::endl; 

  file << std::setw(20) << " Max. 1 escapes " << "||";
  file << std::setw(15) << reco_true_fid_max_1_escapes_true_ccinc;
  file << std::setw(15) << reco_true_fid_max_1_escapes_true_ncinc;
  file << std::setw(15) << reco_true_fid_max_1_escapes_true_ccinc/static_cast<double>(reco_true_fid_max_1_escapes_true_ccinc+reco_true_fid_max_1_escapes_true_ncinc) << std::endl; 

  file << "=====================================================================" << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

