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
  double totPOT = 0.;
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
  p->getValue("TotalPOT",         totPOT);
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

  TopologyMap nc_map      = GeneralAnalysisHelper::GetNuMuNCTopologyMap();
  TopologyMap cc_map      = GeneralAnalysisHelper::GetCCIncTopologyMap();

  // COUNTERS
  unsigned int total_true_ccinc = 0;
  unsigned int total_true_ncinc = 0;
  unsigned int total_truenuoth = 0;
  
  unsigned int true_fid_true_ccinc = 0;
  unsigned int true_fid_true_ncinc = 0;
  unsigned int true_fid_truenuoth = 0;
  
  unsigned int reco_true_fid_true_ccinc = 0;
  unsigned int reco_true_fid_true_ncinc = 0;
  unsigned int reco_true_fid_truenuoth = 0;
  
  unsigned int reco_true_fid_max_1_escapes_true_ccinc = 0;
  unsigned int reco_true_fid_max_1_escapes_true_ncinc = 0;
  unsigned int reco_true_fid_max_1_escapes_truenuoth = 0;
 
  unsigned int reco_true_fid_max_1_escapes_min_1_track_true_ccinc = 0;
  unsigned int reco_true_fid_max_1_escapes_min_1_track_true_ncinc = 0;
  unsigned int reco_true_fid_max_1_escapes_min_1_track_truenuoth = 0;

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

  TH1D *h_chi2_pr_mu = new TH1D("h_chi2_pr_mu", "#chi^{2}_{proton}",40,0,120);
  TH1D *h_chi2_pr_pr = new TH1D("h_chi2_pr_pr", "#chi^{2}_{proton}",40,0,120);
  TH1D *h_chi2_pr_pi = new TH1D("h_chi2_pr_pi", "#chi^{2}_{proton}",40,0,120);

  TH1D *h_chi2_mu_mu = new TH1D("h_chi2_mu_mu", "#chi^{2}_{muon}",40,0,70);
  TH1D *h_chi2_mu_pr = new TH1D("h_chi2_mu_pr", "#chi^{2}_{muon}",40,0,70);
  TH1D *h_chi2_mu_pi = new TH1D("h_chi2_mu_pi", "#chi^{2}_{muon}",40,0,70);
  
  TH1D *h_chi2_ratio_mu = new TH1D("h_chi2_ratio_mu", "#chi^{2}_{muon}/#chi^{2}_{proton}",40,0,0.5);
  TH1D *h_chi2_ratio_pr = new TH1D("h_chi2_ratio_pr", "#chi^{2}_{muon}/#chi^{2}_{proton}",40,0,0.5);
  TH1D *h_chi2_ratio_pi = new TH1D("h_chi2_ratio_pi", "#chi^{2}_{muon}/#chi^{2}_{proton}",40,0,0.5);

  TH2D *h_chi2_mu_chi2_pr_mu = new TH2D("h_chi2_mu_chi2_pr_mu", "#chi^{2}_{muon} vs #chi^{2}_{proton}",40,0,100,40,0,100);
  TH2D *h_chi2_mu_chi2_pr_pr = new TH2D("h_chi2_mu_chi2_pr_pr", "#chi^{2}_{muon} vs #chi^{2}_{proton}",40,0,100,40,0,100);
  TH2D *h_chi2_mu_chi2_pr_pi = new TH2D("h_chi2_mu_chi2_pr_pi", "#chi^{2}_{muon} vs #chi^{2}_{proton}",40,0,100,40,0,100);

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
      if(e.CheckMCNeutrino(14) && e.CheckMCTopology(cc_map)) {
        total_true_ccinc++;
        if(e.IsTrueFiducial()){
          true_fid_true_ccinc++;
          if(e.IsRecoFiducial()){
            reco_true_fid_true_ccinc++;
            if(GeneralAnalysisHelper::MaxOneLongEscapingTrack(e)){
              reco_true_fid_max_1_escapes_true_ccinc++;
              if(GeneralAnalysisHelper::MinOneRecoTrack(e)){
                reco_true_fid_max_1_escapes_min_1_track_true_ccinc++;
              }
            }
          }
        }
      }
      else if(e.CheckMCNeutrino(14)) {
        total_true_ncinc++;
        if(e.IsTrueFiducial()){
          true_fid_true_ncinc++;
          if(e.IsRecoFiducial()){
            reco_true_fid_true_ncinc++;
            if(GeneralAnalysisHelper::MaxOneLongEscapingTrack(e)){
              reco_true_fid_max_1_escapes_true_ncinc++;
              if(GeneralAnalysisHelper::MinOneRecoTrack(e)){
                reco_true_fid_max_1_escapes_min_1_track_true_ncinc++;
              }
            }
          }
        }
      }
      else if(e.CheckMCNeutrino(12) || e.CheckMCNeutrino(-12) || e.CheckMCNeutrino(-14)){
        total_truenuoth++;
        if(e.IsTrueFiducial()){
          true_fid_truenuoth++;
          if(e.IsRecoFiducial()){
            reco_true_fid_truenuoth++;
            if(GeneralAnalysisHelper::MaxOneLongEscapingTrack(e)){
              reco_true_fid_max_1_escapes_truenuoth++;
              if(GeneralAnalysisHelper::MinOneRecoTrack(e)){
                reco_true_fid_max_1_escapes_min_1_track_truenuoth++;
              }
            }
          }
        }
      }
      else{
        std::cout << " Somehow haven't covered everything..." << std::endl;
      }
      // Define and fill the histograms
      // Longest and second longest track lengths
      double longest = -std::numeric_limits<double>::max();
      double second  = -std::numeric_limits<double>::max();
      int longest_id = -1;
      int second_id  = -1;
      GeneralAnalysisHelper::LongestMCTrackLength(e,longest);
      GeneralAnalysisHelper::LongestMCTrackID(e,longest_id);
      for(const Particle &p : e.GetMCParticleList()){
        if(p.ID() != longest_id && p.GetLength() > second){
          second_id = p.ID();
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
          if(p.ID() == longest_id){
            h_longest_diff_mu->Fill(diff);
            h_longest_length_mu->Fill(longest);
          }
        }
        // Pion
        else if(abs(p.GetPdgCode()) == 211){
          h_length_pi->Fill(p.GetLength());
          if(p.ID() == longest_id){
            h_longest_diff_pi->Fill(diff);
            h_longest_length_pi->Fill(longest);
          }
        }
        // Proton
        else if(abs(p.GetPdgCode()) == abs(2212)){
          h_length_pr->Fill(p.GetLength());
          if(p.ID() == longest_id){
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
  SetHistogramStyle(h_muon,   3001, 1, 905, 905, 2, 132, 1, 1.2, true);
  SetHistogramStyle(h_pion,   3001, 1, 801, 801, 2, 132, 1, 1.2, true);
  SetHistogramStyle(h_proton, 3001, 1, 867, 867, 2, 132, 1, 1.2, true);

  // Test
  TCanvas *c = new TCanvas("c","",900,900);
  c->SetLeftMargin  (0.138796 );
  c->SetRightMargin (0.0334448);
  c->SetBottomMargin(0.132404 );
  c->SetTopMargin   (0.0365854);
  
  TLegend * l = new TLegend(0.27,0.88,0.99,0.95);
  l->SetNColumns(3);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.05);
  l->SetTextFont(132);

  // Longest diffs
  l->AddEntry(h_longest_diff_mu,"Longest: #mu","f");
  l->AddEntry(h_longest_diff_pi,"Longest: #pi","f");
  l->AddEntry(h_longest_diff_pr,"Longest: pr","f");

  double max_diff = 1.1 * std::max(h_longest_diff_mu->GetMaximum(), std::max(h_longest_diff_pi->GetMaximum(),h_longest_diff_pr->GetMaximum()));
  h_longest_diff_mu->GetYaxis()->SetRangeUser(0.,max_diff);
  h_longest_diff_mu->SetTitle("");
  h_longest_diff_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_longest_diff_mu->GetXaxis()->SetTitle("#Delta L / L_{longest}, between two longest tracks");
  h_longest_diff_mu->Draw("hist");
  h_longest_diff_pi->Draw("hist same");
  h_longest_diff_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"longest_length_differences.root").c_str());
  c->SaveAs((plots_location+"longest_length_differences.png").c_str());
  c->SaveAs((plots_location+"longest_length_differences.pdf").c_str());
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
  c->SaveAs((plots_location+"longest_lengths.pdf").c_str());
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
  c->SaveAs((plots_location+"track_lengths.pdf").c_str());
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
  h_chi2_mu_mu->GetXaxis()->SetTitle("#chi^{2}_{#mu}");
  h_chi2_mu_mu->Draw("hist");
  h_chi2_mu_pi->Draw("hist same");
  h_chi2_mu_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_mus.root").c_str());
  c->SaveAs((plots_location+"track_chi2_mus.png").c_str());
  c->SaveAs((plots_location+"track_chi2_mus.pdf").c_str());
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
  h_chi2_pr_mu->GetXaxis()->SetTitle("#chi^{2}_{proton}");
  h_chi2_pr_mu->Draw("hist");
  h_chi2_pr_pi->Draw("hist same");
  h_chi2_pr_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_prs.root").c_str());
  c->SaveAs((plots_location+"track_chi2_prs.png").c_str());
  c->SaveAs((plots_location+"track_chi2_prs.pdf").c_str());
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
  h_chi2_ratio_mu->GetXaxis()->SetTitle("#chi^{2}_{#mu} / #chi^{2}_{proton}");
  h_chi2_ratio_mu->Draw("hist");
  h_chi2_ratio_pi->Draw("hist same");
  h_chi2_ratio_pr->Draw("hist same");
  l->Draw();

  c->SaveAs((plots_location+"track_chi2_ratios.root").c_str());
  c->SaveAs((plots_location+"track_chi2_ratios.png").c_str());
  c->SaveAs((plots_location+"track_chi2_ratios.pdf").c_str());
  l->Clear();
  c->Clear();

  // Write statistics
  // POT scaling factor
  double potScale = totPOT / static_cast<double>(pot);
  std::cout << " Total POT: " << totPOT << ", constructed with: " << pot << " POT --> scale = " << potScale << std::endl;
  
  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"analysis_cuts.txt");

  file << "====================================================================================" << std::endl;
  file << " Total POT used to generate this sample: " << pot << std::endl;
  file << std::setw(20) << "Cut/True top."    << "||";
  file <<  std::setw(15) << " CC Inc. ";
  file << std::setw(15) << " NC Inc. ";
  file << std::setw(15) << " Nue Inc. ";
  file << std::setw(15) << " CCInc. purity" << std::endl;

  file << std::setw(20) << " Total "          << "||";
  file << std::setw(15) << static_cast<int>(total_true_ccinc*potScale);
  file << std::setw(15) << static_cast<int>(total_true_ncinc*potScale);
  file << std::setw(15) << static_cast<int>(total_truenuoth*potScale);
  file << std::setw(15) << total_true_ccinc / static_cast<double>(total_true_ccinc+total_true_ncinc+total_truenuoth) << std::endl;

  file << std::setw(20) << " True fiducial "  << "||";
  file << std::setw(15) << static_cast<int>(true_fid_true_ccinc*potScale);
  file << std::setw(15) << static_cast<int>(true_fid_true_ncinc*potScale);
  file << std::setw(15) << static_cast<int>(true_fid_truenuoth*potScale);
  file << std::setw(15) << true_fid_true_ccinc/static_cast<double>(true_fid_true_ccinc+true_fid_true_ncinc+true_fid_truenuoth) << std::endl; 

  file << std::setw(20) << " Reco fiducial "  << "||";
  file << std::setw(15) << static_cast<int>(reco_true_fid_true_ccinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_true_ncinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_truenuoth*potScale);
  file << std::setw(15) << reco_true_fid_true_ccinc/static_cast<double>(reco_true_fid_true_ccinc+reco_true_fid_true_ncinc+reco_true_fid_truenuoth) << std::endl; 

  file << std::setw(20) << " Max. 1 escapes " << "||";
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_true_ccinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_true_ncinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_truenuoth*potScale);
  file << std::setw(15) << reco_true_fid_max_1_escapes_true_ccinc/static_cast<double>(reco_true_fid_max_1_escapes_true_ccinc+reco_true_fid_max_1_escapes_true_ncinc+reco_true_fid_max_1_escapes_truenuoth) << std::endl; 

  file << std::setw(20) << " Min. 1 track " << "||";
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_true_ccinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_true_ncinc*potScale);
  file << std::setw(15) << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_truenuoth*potScale);
  file << std::setw(15) << reco_true_fid_max_1_escapes_min_1_track_true_ccinc/static_cast<double>(reco_true_fid_max_1_escapes_min_1_track_true_ccinc+reco_true_fid_max_1_escapes_min_1_track_true_ncinc+reco_true_fid_max_1_escapes_min_1_track_truenuoth) << std::endl; 

  file << "====================================================================================" << std::endl;

  ofstream fileTeX;
  fileTeX.open(stats_location+"analysis_cuts.tex");
  // TeX file header
  fileTeX << "\\begin{table}[h!] " << std::endl;
  fileTeX << "  \\centering " << std::endl;
  fileTeX << "  \\renewcommand{\\arraystretch}{1.4}" << std::endl;
  fileTeX << "  \\begin{tabular}{ m{3.4cm} * {4}{ >{\\centering\\arraybackslash}m{2.5cm} } }" << endl;
  fileTeX << "    \\hline" << endl;
  fileTeX << "    True topology~$\\rightarrow$ & \\multirow{2}{*}{$\\nu_{\\mu}$~CC} & \\multirow{2}{*}{$\\nu_{\\mu}$~NC} & \\multirow{2}{*}{$\\nu,\\bar{\\nu}$~Other} & \\multirow{2}{*}{$\\nu_{\\mu}$~CC~Purity} \\\\" << std::endl; 
  fileTeX << "    $\\downarrow$~Cut applied & & & \\\\" << std::endl;
  fileTeX << "    \\hline" << endl;

  fileTeX << "    Total & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_true_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_true_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_truenuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << total_true_ccinc / static_cast<double>(total_true_ccinc+total_true_ncinc+total_truenuoth)*100 << "~\\% \\\\ " << std::endl;

  fileTeX << "    True fiducial & "; 
  fileTeX << "    \\num{ " << static_cast<int>(true_fid_true_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(true_fid_true_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(true_fid_truenuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << true_fid_true_ccinc/static_cast<double>(true_fid_true_ccinc+true_fid_true_ncinc+true_fid_truenuoth)*100 << "~\\% \\\\ " << std::endl;

  fileTeX << "    Reco fiducial & "; 
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_true_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_true_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_truenuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << reco_true_fid_true_ccinc/static_cast<double>(reco_true_fid_true_ccinc+reco_true_fid_true_ncinc+reco_true_fid_truenuoth)*100 << "~\\% \\\\ "  << std::endl;

  fileTeX << "    Max. 1 escapes & "; 
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_true_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_true_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_truenuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << reco_true_fid_max_1_escapes_true_ccinc/static_cast<double>(reco_true_fid_max_1_escapes_true_ccinc+reco_true_fid_max_1_escapes_true_ncinc+reco_true_fid_max_1_escapes_truenuoth)*100 << "~\\% \\\\ "  << std::endl;

  fileTeX << "    Min. 1 track & "; 
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_true_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_true_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(reco_true_fid_max_1_escapes_min_1_track_truenuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << reco_true_fid_max_1_escapes_min_1_track_true_ccinc/static_cast<double>(reco_true_fid_max_1_escapes_min_1_track_true_ccinc+reco_true_fid_max_1_escapes_min_1_track_true_ncinc+reco_true_fid_max_1_escapes_min_1_track_truenuoth)*100 << "~\\% \\\\" << std::endl; 
  fileTeX << "    \\hline" << endl;

  // Tex file end
  fileTeX << "  \\end{tabular}"<< endl;
  fileTeX << "\\end{table}" << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

