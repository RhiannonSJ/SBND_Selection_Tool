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

void SetHistErrors(TH1D *h, TH1D *den);

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
  unsigned int detector = 0; // 0 = sbnd, 1 = uboone, 2 = icarus
  double totPOT = 0.;
  double diff_cut = 0;
  double length_cut = 0;
  double longest_cut = 0;
  double chi2p_cut = 0;
  double chi2mu_cut = 0;
  double chi2ratio_cut = 0;
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
  p->getValue("PreDiffCut",       diff_cut);
  p->getValue("PreLengthCut",     length_cut);
  p->getValue("PreLongestCut",    longest_cut);
  p->getValue("PreChi2PCut",      chi2p_cut);
  p->getValue("PreChi2MuCut",     chi2mu_cut);
  p->getValue("PreChi2RatioCut",  chi2ratio_cut);

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
  TopologyMap cc0pi_map   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map   = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc0pi2p_map = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap nc0pi_map   = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_map   = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nue_map     = GeneralAnalysisHelper::GetNuECCTopologyMap();

  std::vector< TopologyMap > maps({cc_map, cc0pi_map, cc1pi_map, cc0pi2p_map, nc_map, nc0pi_map, nc1pi_map, nue_map});  

  // HISTOGRAMS

  // Longest/One escaping
  TH1D *h_one_escapes_enu = new TH1D("h_one_escapes_enu", "One track escapes",40,0,3);
  TH1D *h_one_escapes_pmu = new TH1D("h_one_escapes_pmu", "One track escapes",40,0,3);
  TH1D *h_one_escapes_pth = new TH1D("h_one_escapes_pth", "One track escapes",40,-1,1);

  TH1D *h_long_escapes_enu = new TH1D("h_long_escapes_enu", "Longest track escapes",40,0,3);
  TH1D *h_long_escapes_pmu = new TH1D("h_long_escapes_pmu", "Longest track escapes",40,0,3);
  TH1D *h_long_escapes_pth = new TH1D("h_long_escapes_pth", "Longest track escapes",40,-1,1);

  // Length cut/none escaping
  TH1D *h_none_escape_enu = new TH1D("h_none_escape_enu", "No tracks escape",40,0,3);
  TH1D *h_none_escape_pmu = new TH1D("h_none_escape_pmu", "No tracks escape",40,0,3);
  TH1D *h_none_escape_pth = new TH1D("h_none_escape_pth", "No tracks escape",40,-1,1);

  TH1D *h_min_length_enu = new TH1D("h_min_length_enu", "Minimum track length / No tracks escape",40,0,3);
  TH1D *h_min_length_pmu = new TH1D("h_min_length_pmu", "Minimum track length / No tracks escape",40,0,3);
  TH1D *h_min_length_pth = new TH1D("h_min_length_pth", "Minimum track length / No tracks escape",40,-1,1);

  // Diff/min2
  TH1D *h_min2_enu = new TH1D("h_min2_enu", "Minimum 2 tracks",40,0,3);
  TH1D *h_min2_pmu = new TH1D("h_min2_pmu", "Minimum 2 tracks",40,0,3);
  TH1D *h_min2_pth = new TH1D("h_min2_pth", "Minimum 2 tracks",40,-1,1);

  TH1D *h_diff_enu = new TH1D("h_diff_enu", "Length difference cuts / Minimum 2 tracks",40,0,3);
  TH1D *h_diff_pmu = new TH1D("h_diff_pmu", "Length difference cuts / Minimum 2 tracks",40,0,3);
  TH1D *h_diff_pth = new TH1D("h_diff_pth", "Length difference cuts / Minimum 2 tracks",40,-1,1);

  // Chi2 & length/ 1 track or no diff
  TH1D *h_1track_nodiff_enu = new TH1D("h_1track_nodiff_enu", "1 track no diff",40,0,3);
  TH1D *h_1track_nodiff_pmu = new TH1D("h_1track_nodiff_pmu", "1 track no diff",40,0,3);
  TH1D *h_1track_nodiff_pth = new TH1D("h_1track_nodiff_pth", "1 track no diff",40,-1,1);

  TH1D *h_chi2_length_enu = new TH1D("h_chi2_length_enu", "Longest length & #chi^{2}_{muon} & #chi^{2}_{proton} cuts / 1 track or didn't pass diff",40,0,3);
  TH1D *h_chi2_length_pmu = new TH1D("h_chi2_length_pmu", "Longest length & #chi^{2}_{muon} & #chi^{2}_{proton} cuts / 1 track or didn't pass diff ",40,0,3);
  TH1D *h_chi2_length_pth = new TH1D("h_chi2_length_pth", "Longest length & #chi^{2}_{muon} & #chi^{2}_{proton} cuts / 1 track or didn't pass diff ",40,-1,1);

  // COUNTERS
  unsigned int total_ccinc  = 0;
  unsigned int total_ncinc  = 0;
  unsigned int total_nuoth = 0;
  
  unsigned int total_fiducial_cuts_ccinc  = 0;
  unsigned int total_fiducial_cuts_ncinc  = 0;
  unsigned int total_fiducial_cuts_nuoth = 0;

  unsigned int one_escaping_ccinc = 0;
  unsigned int one_escaping_ncinc = 0;
  unsigned int one_escaping_nuoth = 0;

  unsigned int longest_escaping_ccinc = 0;
  unsigned int longest_escaping_ncinc = 0;
  unsigned int longest_escaping_nuoth = 0;

  // If all tracks are contained
  // Considering muons
  unsigned int track_length_cut_ccinc = 0;
  unsigned int track_length_cut_ncinc = 0;
  unsigned int track_length_cut_nuoth = 0;

  unsigned int track_length_min2_ccinc = 0;
  unsigned int track_length_min2_ncinc = 0;
  unsigned int track_length_min2_nuoth = 0;

  unsigned int track_length_diff_cut_ccinc = 0;
  unsigned int track_length_diff_cut_ncinc = 0;
  unsigned int track_length_diff_cut_nuoth = 0;

  // Considering protons
  unsigned int chi2_cuts_ccinc = 0;
  unsigned int chi2_cuts_ncinc = 0;
  unsigned int chi2_cuts_nuoth = 0;

  // Considering protons
  unsigned int all_cuts_ccinc = 0;
  unsigned int all_cuts_ncinc = 0;
  unsigned int all_cuts_nuoth = 0;

  // First, ensure all tracks are contained
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    for(const Event &e : events){
      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial()) continue;
      if(e.CheckMCTopology(cc_map))
        total_ccinc++;
      else if(e.CheckMCNeutrino(14) && e.CheckMCTopology(nc_map))
        total_ncinc++;
      else if(e.CheckMCNeutrino(12) || e.CheckMCNeutrino(-12) || e.CheckMCNeutrino(-14))
        total_nuoth++;
      else{
        std::cout << " Missed something" << std::endl;
      }

      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial() || 
         !GeneralAnalysisHelper::MaxOneLongEscapingTrack(e) || 
         !GeneralAnalysisHelper::MinOneRecoTrack(e)) continue;

      double enu = e.GetTrueNuEnergy();

      // Get the longest and second longest tracks
      double longest = -std::numeric_limits<double>::max();
      double second  = -std::numeric_limits<double>::max();
      int longest_id = -1;
      int second_id  = -1;
      GeneralAnalysisHelper::LongestRecoTrackLength(e,longest);
      GeneralAnalysisHelper::LongestRecoTrackID(e,longest_id);
      for(const Particle &p : e.GetRecoParticleList()){
        if(!p.GetFromRecoTrack()) continue;
        if(p.ID() != longest_id && p.GetLength() > second){
          second_id = p.ID();
          second = p.GetLength();
        }
      }

      // Get the fractional difference between the longest and second longest track lengths
      double diff = -999.;
      bool min_2_tracks = (longest_id != -1 && second_id != -1);
      if(min_2_tracks)
        diff = (longest - second)/longest;

      if(e.CheckMCTopology(cc_map)){
        total_fiducial_cuts_ccinc++;
        double pmu = 0.;
        double pth = 0.;
        std::vector<float> pmus, pths;
        GeneralAnalysisHelper::GetMCModulusMomentumWithPdg(e,13,pmus);
        GeneralAnalysisHelper::GetMCCosThetaWithPdg(e,13,pths);
        if(pmus.size() > 1 || pths.size() > 1){
          std::cerr << " Error: Found more than 2 muons in a no-pileup sample " << std::endl;
        }
        pmu = pmus.at(0);
        pth = pths.at(0);

        // If 1 track escapes
        if(GeneralAnalysisHelper::NumberEscapingTracks(e) == 1){
          one_escaping_ccinc++;
          all_cuts_ccinc++;
          h_one_escapes_enu->Fill(enu);
          h_one_escapes_pmu->Fill(pmu);
          h_one_escapes_pth->Fill(pth);
          for(const Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            if(p.ID() == longest_id && p.GetOneEndTrackContained() && longest > longest_cut){
              h_long_escapes_enu->Fill(enu);
              h_long_escapes_pmu->Fill(pmu);
              h_long_escapes_pth->Fill(pth);
              longest_escaping_ccinc++;
              break;
            }
          }
        }
        // If no tracks escape
        // This is where the PID cuts actually matter
        // Start trying to find the muon
        else{
          h_none_escape_enu->Fill(enu);
          h_none_escape_pmu->Fill(pmu);
          h_none_escape_pth->Fill(pth);
          // Not considering protons first
          bool track_cut_passed = false;
          for(Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            if(p.GetLength() > length_cut){
              track_length_cut_ccinc++; 
              track_cut_passed = true;
              break;
            }
          }
          if(track_cut_passed){
            h_min_length_enu->Fill(enu);
            h_min_length_pmu->Fill(pmu);
            h_min_length_pth->Fill(pth);
            if(min_2_tracks){
              h_min2_enu->Fill(enu);
              h_min2_pmu->Fill(pmu);
              h_min2_pth->Fill(pth);
              track_length_min2_ccinc++;
              if(diff > diff_cut){
                h_diff_enu->Fill(enu);
                h_diff_pmu->Fill(pmu);
                h_diff_pth->Fill(pth);
                track_length_diff_cut_ccinc++;
                all_cuts_ccinc++;
              }
            }
            if((min_2_tracks && diff <= diff_cut) || !min_2_tracks){
              h_1track_nodiff_enu->Fill(enu);
              h_1track_nodiff_pmu->Fill(pmu);
              h_1track_nodiff_pth->Fill(pth);
              bool found_p_mu = false;
              bool found_ratio = false;
              // Now consider chi2 variables
              for(const Particle &p : e.GetRecoParticleList()){
                if(!p.GetFromRecoTrack()) continue;
                //Check for clear protons
                if(detector != 2){
                  if(found_p_mu == false && ((p.GetChi2Mu()/p.GetChi2P()) < chi2ratio_cut || (p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    h_chi2_length_enu->Fill(enu);
                    h_chi2_length_pmu->Fill(pmu);
                    h_chi2_length_pth->Fill(pth);
                    found_p_mu = true;
                    chi2_cuts_ccinc++;
                    all_cuts_ccinc++;
                  }
                }
                else{
                  if(found_p_mu == false && ((p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    h_chi2_length_enu->Fill(enu);
                    h_chi2_length_pmu->Fill(pmu);
                    h_chi2_length_pth->Fill(pth);
                    found_p_mu = true;
                    chi2_cuts_ccinc++;
                    all_cuts_ccinc++;
                  }
                }
              }
            }
          }
        }
      }
      else if(e.CheckMCNeutrino(14) && e.CheckMCTopology(nc_map)) {
        total_fiducial_cuts_ncinc++;
        if(GeneralAnalysisHelper::NumberEscapingTracks(e) == 1){
          one_escaping_ncinc++;
          for(const Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            all_cuts_ncinc++;
            if(p.ID() == longest_id && p.GetOneEndTrackContained() && longest > longest_cut){
              longest_escaping_ncinc++;
              break;
            }
          }
        }
        else{
          bool track_cut_passed = false;
          for(Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            if(p.GetLength() > length_cut){
              track_length_cut_ncinc++; 
              track_cut_passed = true;
              break;
            }
          }
          if(track_cut_passed){
            if(min_2_tracks){
              track_length_min2_ncinc++;
              if(diff > diff_cut){
                track_length_diff_cut_ncinc++;
                all_cuts_ncinc++;
              }
            }
            if((min_2_tracks && diff <= diff_cut) || !min_2_tracks){
              // Now consider chi2 variables
              bool found_p_mu = false;
              for(const Particle &p : e.GetRecoParticleList()){
                if(!p.GetFromRecoTrack()) continue;
                //Check for clear protons
                if(detector != 2){
                  if(found_p_mu == false && ((p.GetChi2Mu()/p.GetChi2P()) < chi2ratio_cut || (p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    found_p_mu = true;
                    chi2_cuts_ncinc++;
                    all_cuts_ncinc++;
                  }
                }
                else{
                  if(found_p_mu == false && ((p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    found_p_mu = true;
                    chi2_cuts_ncinc++;
                    all_cuts_ncinc++;
                  }
                }
              }
            }
          }
        }
      }
      else if(e.CheckMCNeutrino(12) || e.CheckMCNeutrino(-12) || e.CheckMCNeutrino(-14)){
        total_fiducial_cuts_nuoth++;
        if(GeneralAnalysisHelper::NumberEscapingTracks(e) == 1){
          one_escaping_nuoth++;
          for(const Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            all_cuts_nuoth++;
            if(p.ID() == longest_id && p.GetOneEndTrackContained() && longest > longest_cut){
              longest_escaping_nuoth++;
              break;
            }
          }
        }
        else{
          bool track_cut_passed = false;
          for(Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            if(p.GetLength() > length_cut){
              track_length_cut_nuoth++; 
              track_cut_passed = true;
              break;
            }
          }
          if(track_cut_passed){
            if(min_2_tracks){
              track_length_min2_nuoth++;
              if(diff > diff_cut){
                track_length_diff_cut_nuoth++;
                all_cuts_nuoth++;
              }
            }
            if((min_2_tracks && diff <= diff_cut) || !min_2_tracks){
              // Now consider chi2 variables
              bool found_p_mu = false;
              for(const Particle &p : e.GetRecoParticleList()){
                if(!p.GetFromRecoTrack()) continue;
                //Check for clear protons
                if(detector != 2){
                  if(found_p_mu == false && ((p.GetChi2Mu()/p.GetChi2P()) < chi2ratio_cut || (p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    found_p_mu = true;
                    chi2_cuts_nuoth++;
                    all_cuts_nuoth++;
                  }
                }
                else{
                  if(found_p_mu == false && ((p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)){
                    found_p_mu = true;
                    chi2_cuts_nuoth++;
                    all_cuts_nuoth++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Construct and write efficiency histograms
  TH1D *h_eff_min_length_enu = (TH1D*) h_min_length_enu->Clone("h_eff_min_length_enu");
  TH1D *h_eff_min_length_pmu = (TH1D*) h_min_length_pmu->Clone("h_eff_min_length_pmu");
  TH1D *h_eff_min_length_pth = (TH1D*) h_min_length_pth->Clone("h_eff_min_length_pth");

  TH1D *h_eff_diff_enu = (TH1D*) h_diff_enu->Clone("h_eff_diff_enu");
  TH1D *h_eff_diff_pmu = (TH1D*) h_diff_pmu->Clone("h_eff_diff_pmu");
  TH1D *h_eff_diff_pth = (TH1D*) h_diff_pth->Clone("h_eff_diff_pth");

  TH1D *h_eff_chi2_length_enu = (TH1D*) h_chi2_length_enu->Clone("h_eff_chi2_length_enu");
  TH1D *h_eff_chi2_length_pmu = (TH1D*) h_chi2_length_pmu->Clone("h_eff_chi2_length_pmu");
  TH1D *h_eff_chi2_length_pth = (TH1D*) h_chi2_length_pth->Clone("h_eff_chi2_length_pth");

  // Now divide through by the appropraite histogram and set the style
  h_eff_min_length_enu->Divide(h_none_escape_enu);
  h_eff_min_length_pmu->Divide(h_none_escape_pmu);
  h_eff_min_length_pth->Divide(h_none_escape_pth);

  h_eff_diff_enu->Divide(h_min2_enu);
  h_eff_diff_pmu->Divide(h_min2_pmu);
  h_eff_diff_pth->Divide(h_min2_pth);

  h_eff_chi2_length_enu->Divide(h_1track_nodiff_enu);
  h_eff_chi2_length_pmu->Divide(h_1track_nodiff_pmu);
  h_eff_chi2_length_pth->Divide(h_1track_nodiff_pth);

  std::vector<TH1D*> h_min_length  = {h_eff_min_length_enu,  h_eff_min_length_pmu,  h_eff_min_length_pth};
  std::vector<TH1D*> h_diff        = {h_eff_diff_enu,        h_eff_diff_pmu,        h_eff_diff_pth};
  std::vector<TH1D*> h_chi2_length = {h_eff_chi2_length_enu, h_eff_chi2_length_pmu, h_eff_chi2_length_pth};

  SetHistogramStyle(h_min_length,  3001, 1, 801, 801, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_diff,        3001, 1, 867, 867, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_chi2_length, 3001, 1, 835, 835, 2, 132, 1, 1.2, false);

  TCanvas *c = new TCanvas("c","",900,900);
  c->SetLeftMargin  (0.138796 );
  c->SetRightMargin (0.0334448);
  c->SetBottomMargin(0.132404 );
  c->SetTopMargin   (0.063 );
  
  gStyle->SetTitleFont(132,"t");

  h_eff_min_length_enu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_eff_min_length_enu->Draw("hist");

  c->SaveAs((plots_location+"enu_efficiency_min_length.root").c_str());
  c->SaveAs((plots_location+"enu_efficiency_min_length.png").c_str());
  c->SaveAs((plots_location+"enu_efficiency_min_length.pdf").c_str());
  c->Clear();

  h_eff_min_length_pmu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_eff_min_length_pmu->Draw("hist");

  c->SaveAs((plots_location+"pmu_efficiency_min_length.root").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_min_length.png").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_min_length.pdf").c_str());
  c->Clear();

  h_eff_min_length_pth->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_pth->GetXaxis()->SetTitle("True muon scattering angle");
  h_eff_min_length_pth->Draw("hist");

  c->SaveAs((plots_location+"pth_efficiency_min_length.root").c_str());
  c->SaveAs((plots_location+"pth_efficiency_min_length.png").c_str());
  c->SaveAs((plots_location+"pth_efficiency_min_length.pdf").c_str());
  c->Clear();

  h_eff_diff_enu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_diff_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_eff_diff_enu->Draw("hist");

  c->SaveAs((plots_location+"enu_efficiency_diff.root").c_str());
  c->SaveAs((plots_location+"enu_efficiency_diff.png").c_str());
  c->SaveAs((plots_location+"enu_efficiency_diff.pdf").c_str());
  c->Clear();

  h_eff_diff_pmu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_diff_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_eff_diff_pmu->Draw("hist");

  c->SaveAs((plots_location+"pmu_efficiency_diff.root").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_diff.png").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_diff.pdf").c_str());
  c->Clear();

  h_eff_diff_pth->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_diff_pth->GetXaxis()->SetTitle("True muon scattering angle");
  h_eff_diff_pth->Draw("hist");

  c->SaveAs((plots_location+"pth_efficiency_diff.root").c_str());
  c->SaveAs((plots_location+"pth_efficiency_diff.png").c_str());
  c->SaveAs((plots_location+"pth_efficiency_diff.pdf").c_str());
  c->Clear();

  h_eff_chi2_length_enu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_chi2_length_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_eff_chi2_length_enu->Draw("hist");

  c->SaveAs((plots_location+"enu_efficiency_chi2_length.root").c_str());
  c->SaveAs((plots_location+"enu_efficiency_chi2_length.png").c_str());
  c->SaveAs((plots_location+"enu_efficiency_chi2_length.pdf").c_str());
  c->Clear();

  h_eff_chi2_length_pmu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_chi2_length_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_eff_chi2_length_pmu->Draw("hist");

  c->SaveAs((plots_location+"pmu_efficiency_chi2_length.root").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_chi2_length.png").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_chi2_length.pdf").c_str());
  c->Clear();

  h_eff_chi2_length_pth->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_chi2_length_pth->GetXaxis()->SetTitle("True muon scattering angle");
  h_eff_chi2_length_pth->Draw("hist");

  c->SaveAs((plots_location+"pth_efficiency_chi2_length.root").c_str());
  c->SaveAs((plots_location+"pth_efficiency_chi2_length.png").c_str());
  c->SaveAs((plots_location+"pth_efficiency_chi2_length.pdf").c_str());
  c->Clear();

  // Now overlay the histograms
  SetHistogramStyle(h_min_length,  0, 1, 801, 801, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_diff,        0, 2, 867, 867, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_chi2_length, 0, 1, 835, 835, 2, 132, 1, 1.2, false);

  // Not going to print title so change canvas margin
  c->SetTopMargin(0.0734);

  TLegend * l = new TLegend(0.13,0.94,0.97,0.99);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.034);
  l->SetNColumns(3);
  l->SetMargin(0.2);
  l->SetTextFont(132);

  l->AddEntry(h_eff_min_length_enu,  "One track min length","l");
  l->AddEntry(h_eff_diff_enu,        "#Delta L / L_{longest}","l");
  l->AddEntry(h_eff_chi2_length_enu, "#chi^{2} & length cuts","l");

  h_eff_min_length_enu->SetTitle("");
  h_eff_min_length_enu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_eff_min_length_enu->GetYaxis()->SetRangeUser(0.2,1.1);
  h_eff_min_length_enu->Draw("hist");
  h_eff_diff_enu->Draw("hist same");
  h_eff_chi2_length_enu->Draw("hist same");
  l->Draw("same");

  c->SaveAs((plots_location+"enu_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"enu_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"enu_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  l->AddEntry(h_eff_min_length_pmu,  "One track min length","l");
  l->AddEntry(h_eff_diff_pmu,        "#Delta L / L_{longest}","l");
  l->AddEntry(h_eff_chi2_length_pmu, "#chi^{2} & length cuts","l");

  h_eff_min_length_pmu->SetTitle("");
  h_eff_min_length_pmu->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_eff_min_length_pmu->GetYaxis()->SetRangeUser(0.2,1.1);
  h_eff_min_length_pmu->Draw("hist");
  h_eff_diff_pmu->Draw("hist same");
  h_eff_chi2_length_pmu->Draw("hist same");
  l->Draw("same");

  c->SaveAs((plots_location+"pmu_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"pmu_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  l->AddEntry(h_eff_min_length_pth,  "One track min length","l");
  l->AddEntry(h_eff_diff_pth,        "#Delta L / L_{longest}","l");
  l->AddEntry(h_eff_chi2_length_pth, "#chi^{2} & length cuts","l");

  h_eff_min_length_pth->SetTitle("");
  h_eff_min_length_pth->GetYaxis()->SetTitle("Efficiency w.r.t previous cut");
  h_eff_min_length_pth->GetXaxis()->SetTitle("True muon scattering angle");
  h_eff_min_length_pth->GetYaxis()->SetRangeUser(0.2,1.1);
  h_eff_min_length_pth->Draw("hist");
  h_eff_diff_pth->Draw("hist same");
  h_eff_chi2_length_pth->Draw("hist same");
  l->Draw("same");

  c->SaveAs((plots_location+"pth_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"pth_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"pth_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  // Cumulative efficiencies
  TH1D *h_cum_total_enu = (TH1D*) h_one_escapes_enu->Clone("h_cum_total_enu");
  TH1D *h_cum_total_pmu = (TH1D*) h_one_escapes_pmu->Clone("h_cum_total_pmu");
  TH1D *h_cum_total_pth = (TH1D*) h_one_escapes_pth->Clone("h_cum_total_pth");

  TH1D *h_cum_one_escapes_enu = (TH1D*) h_one_escapes_enu->Clone("h_cum_one_escapes_enu");
  TH1D *h_cum_one_escapes_pmu = (TH1D*) h_one_escapes_pmu->Clone("h_cum_one_escapes_pmu");
  TH1D *h_cum_one_escapes_pth = (TH1D*) h_one_escapes_pth->Clone("h_cum_one_escapes_pth");

  TH1D *h_cum_min_length_enu = (TH1D*) h_min_length_enu->Clone("h_cum_min_length_enu");
  TH1D *h_cum_min_length_pmu = (TH1D*) h_min_length_pmu->Clone("h_cum_min_length_pmu");
  TH1D *h_cum_min_length_pth = (TH1D*) h_min_length_pth->Clone("h_cum_min_length_pth");

  TH1D *h_cum_diff_enu = (TH1D*) h_diff_enu->Clone("h_cum_diff_enu");
  TH1D *h_cum_diff_pmu = (TH1D*) h_diff_pmu->Clone("h_cum_diff_pmu");
  TH1D *h_cum_diff_pth = (TH1D*) h_diff_pth->Clone("h_cum_diff_pth");

  TH1D *h_cum_chi2_length_enu = (TH1D*) h_chi2_length_enu->Clone("h_cum_chi2_length_enu");
  TH1D *h_cum_chi2_length_pmu = (TH1D*) h_chi2_length_pmu->Clone("h_cum_chi2_length_pmu");
  TH1D *h_cum_chi2_length_pth = (TH1D*) h_chi2_length_pth->Clone("h_cum_chi2_length_pth");

  // Add one and none escaping together to get the denominator
  h_cum_total_enu->Add(h_none_escape_enu);
  h_cum_total_pmu->Add(h_none_escape_pmu);
  h_cum_total_pth->Add(h_none_escape_pth);

  h_cum_diff_enu->Add(h_cum_one_escapes_enu);
  h_cum_diff_pmu->Add(h_cum_one_escapes_pmu);
  h_cum_diff_pth->Add(h_cum_one_escapes_pth);

  h_cum_chi2_length_enu->Add(h_cum_diff_enu);
  h_cum_chi2_length_pmu->Add(h_cum_diff_pmu);
  h_cum_chi2_length_pth->Add(h_cum_diff_pth);

  // Now divide by the totals
  h_cum_one_escapes_enu->Divide(h_cum_total_enu);
  h_cum_one_escapes_pmu->Divide(h_cum_total_pmu);
  h_cum_one_escapes_pth->Divide(h_cum_total_pth);

  SetHistErrors(h_cum_one_escapes_enu,h_cum_total_enu);
  SetHistErrors(h_cum_one_escapes_pmu,h_cum_total_pmu);
  SetHistErrors(h_cum_one_escapes_pth,h_cum_total_pth);
 
  h_cum_diff_enu->Divide(h_cum_total_enu);
  h_cum_diff_pmu->Divide(h_cum_total_pmu);
  h_cum_diff_pth->Divide(h_cum_total_pth);

  SetHistErrors(h_cum_diff_enu,h_cum_total_enu);
  SetHistErrors(h_cum_diff_pmu,h_cum_total_pmu);
  SetHistErrors(h_cum_diff_pth,h_cum_total_pth);
  
  h_cum_chi2_length_enu->Divide(h_cum_total_enu);
  h_cum_chi2_length_pmu->Divide(h_cum_total_pmu);
  h_cum_chi2_length_pth->Divide(h_cum_total_pth);

  SetHistErrors(h_cum_chi2_length_enu,h_cum_total_enu);
  SetHistErrors(h_cum_chi2_length_pmu,h_cum_total_pmu);
  SetHistErrors(h_cum_chi2_length_pth,h_cum_total_pth);
  
  std::vector<TH1D*> h_cum_one_escapes = {h_cum_one_escapes_enu, h_cum_one_escapes_pmu, h_cum_one_escapes_pth};
  std::vector<TH1D*> h_cum_diff        = {h_cum_diff_enu,        h_cum_diff_pmu,        h_cum_diff_pth};
  std::vector<TH1D*> h_cum_chi2_length = {h_cum_chi2_length_enu, h_cum_chi2_length_pmu, h_cum_chi2_length_pth};

  SetHistogramStyle(h_cum_one_escapes, 0, 1, 905, 905, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_cum_diff,        0, 1, 867, 867, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_cum_chi2_length, 0, 1, 835, 835, 2, 132, 1, 1.2, false);
  
  l->SetX1NDC(0.13);
  l->SetY1NDC(0.94);
  l->SetX2NDC(0.97);
  l->SetY2NDC(0.99);
  
  // Now overlay the histograms
  l->AddEntry(h_cum_one_escapes_enu, "One track escapes","ep");
  l->AddEntry(h_cum_diff_enu,        "#Delta L / L_{longest}","ep");
  l->AddEntry(h_cum_chi2_length_enu, "#chi^{2} & length cuts","ep");

  h_cum_one_escapes_enu->SetTitle("");
  h_cum_one_escapes_enu->GetYaxis()->SetTitle("Cumulative efficiency");
  h_cum_one_escapes_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_cum_one_escapes_enu->GetYaxis()->SetRangeUser(0.0,1.1);
  h_cum_one_escapes_enu->Draw("E1 X0");
  h_cum_diff_enu->Draw("E1 X0 same");
  h_cum_chi2_length_enu->Draw("E1 X0 same");
  l->Draw("same");

  c->SaveAs((plots_location+"enu_cumulative_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"enu_cumulative_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"enu_cumulative_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();
  
  l->AddEntry(h_cum_one_escapes_pmu, "One track escapes","ep");
  l->AddEntry(h_cum_diff_pmu,        "#Delta L / L_{longest}","ep");
  l->AddEntry(h_cum_chi2_length_pmu, "#chi^{2} & length cuts","ep");

  h_cum_one_escapes_pmu->SetTitle("");
  h_cum_one_escapes_pmu->GetYaxis()->SetTitle("Cumulative efficiency");
  h_cum_one_escapes_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_cum_one_escapes_pmu->GetYaxis()->SetRangeUser(0.0,1.1);
  h_cum_one_escapes_pmu->Draw("E1 X0");
  h_cum_diff_pmu->Draw("E1 X0 same");
  h_cum_chi2_length_pmu->Draw("E1 X0 same");
  l->Draw("same");

  c->SaveAs((plots_location+"pmu_cumulative_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"pmu_cumulative_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"pmu_cumulative_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  l->AddEntry(h_cum_one_escapes_pth, "One track escapes","ep");
  l->AddEntry(h_cum_diff_pth,        "#Delta L / L_{longest}","ep");
  l->AddEntry(h_cum_chi2_length_pth, "#chi^{2} & length cuts","ep");

  h_cum_one_escapes_pth->SetTitle("");
  h_cum_one_escapes_pth->GetYaxis()->SetTitle("Cumulative efficiency");
  h_cum_one_escapes_pth->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h_cum_one_escapes_pth->GetYaxis()->SetRangeUser(0.0,1.1);
  h_cum_one_escapes_pth->Draw("E1 X0");
  h_cum_diff_pth->Draw("E1 X0 same");
  h_cum_chi2_length_pth->Draw("E1 X0 same");
  l->Draw("same");

  c->SaveAs((plots_location+"pth_cumulative_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"pth_cumulative_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"pth_cumulative_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  // Write statistics
  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"preselection_cut_rates.txt");

  file << "===================================================================================================" << std::endl;
  file << " Total POT used to generate this sample: " << pot << std::endl;
  file << std::setw(20) << "Cut/True top."    << "||";
  file << std::setw(15) << " CC Inc. ";
  file << std::setw(15) << " NC Inc. ";
  file << std::setw(15) << " Nue Inc. ";
  file << std::setw(15) << " CCInc. purity" << std::endl;

  file << std::setw(20) << " Total " << "||";
  file << std::setw(15) << total_ccinc;
  file << std::setw(15) << total_ncinc;
  file << std::setw(15) << total_nuoth;
  file << std::setw(15) << total_ccinc / static_cast<double>(total_ccinc+total_ncinc+total_nuoth) << std::endl;

  file << std::setw(20) << " Total fiducial " << "||";
  file << std::setw(15) << total_fiducial_cuts_ccinc;
  file << std::setw(15) << total_fiducial_cuts_ncinc;
  file << std::setw(15) << total_fiducial_cuts_nuoth;
  file << std::setw(15) << total_fiducial_cuts_ccinc / static_cast<double>(total_fiducial_cuts_ccinc+total_fiducial_cuts_ncinc+total_fiducial_cuts_nuoth) << std::endl;

  file << std::setw(20) << " All cuts "   << "||";
  file << std::setw(15) << all_cuts_ccinc;
  file << std::setw(15) << all_cuts_ncinc;
  file << std::setw(15) << all_cuts_nuoth;
  file << std::setw(15) << all_cuts_ccinc / static_cast<double>(all_cuts_ccinc+all_cuts_ncinc+all_cuts_nuoth) << std::endl;


  file << "===================================================================================================" << std::endl;

  file << std::setw(20) << " One escaping "    << "||";
  file << std::setw(15) << one_escaping_ccinc;
  file << std::setw(15) << one_escaping_ncinc;
  file << std::setw(15) << one_escaping_nuoth;
  file << std::setw(15) << one_escaping_ccinc / static_cast<double>(one_escaping_ccinc+one_escaping_ncinc+one_escaping_nuoth) << std::endl;

  file << "---------------------------------------------------------------------------------------------------" << std::endl;

  file << std::setw(20) << " None escaping "    << "||";
  file << std::setw(15) << total_fiducial_cuts_ccinc-one_escaping_ccinc;
  file << std::setw(15) << total_fiducial_cuts_ncinc-one_escaping_ncinc;
  file << std::setw(15) << total_fiducial_cuts_nuoth-one_escaping_nuoth;
  file << std::setw(15) << (total_fiducial_cuts_ccinc-one_escaping_ccinc) / static_cast<double>((total_fiducial_cuts_ccinc-one_escaping_ccinc)+(total_fiducial_cuts_ncinc-one_escaping_ncinc)+(total_fiducial_cuts_nuoth-one_escaping_nuoth)) << std::endl;

  file << "---------------------------------------------------------------------------------------------------" << std::endl;
  file << " If all tracks are contained: "                                        << std::endl;
  file << "---------------------------------------------------------------------------------------------------" << std::endl;

  file << std::setw(20) << " 1 track > 10cm "   << "||";
  file << std::setw(15) << track_length_cut_ccinc;
  file << std::setw(15) << track_length_cut_ncinc;
  file << std::setw(15) << track_length_cut_nuoth;
  file << std::setw(15) << track_length_cut_ccinc / static_cast<double>(track_length_cut_ccinc+track_length_cut_ncinc+track_length_cut_nuoth) << std::endl;

  file << std::setw(20) << " Min 2 tracks "   << "||";
  file << std::setw(15) << track_length_min2_ccinc;
  file << std::setw(15) << track_length_min2_ncinc;
  file << std::setw(15) << track_length_min2_nuoth;
  file << std::setw(15) << track_length_min2_ccinc / static_cast<double>(track_length_min2_ccinc+track_length_min2_ncinc+track_length_min2_nuoth) << std::endl;

  file << "---------------------------------------------------------------------------------------------------" << std::endl;
  file << " If all tracks are contained & track cut passed: "                     << std::endl;
  file << "---------------------------------------------------------------------------------------------------" << std::endl;

  file << std::setw(20) << " Diff "   << "||";
  file << std::setw(15) << track_length_diff_cut_ccinc;
  file << std::setw(15) << track_length_diff_cut_ncinc;
  file << std::setw(15) << track_length_diff_cut_nuoth;
  file << std::setw(15) << track_length_diff_cut_ccinc / static_cast<double>(track_length_diff_cut_ccinc+track_length_diff_cut_ncinc+track_length_diff_cut_nuoth) << std::endl;

  file << std::setw(20) << " Chi2 & length "   << "||";
  file << std::setw(15) << chi2_cuts_ccinc;
  file << std::setw(15) << chi2_cuts_ncinc;
  file << std::setw(15) << chi2_cuts_nuoth;
  file << std::setw(15) << chi2_cuts_ccinc / static_cast<double>(chi2_cuts_ccinc+chi2_cuts_ncinc+chi2_cuts_nuoth) << std::endl;

  file << "===================================================================================================" << std::endl;

  double potScale = totPOT / static_cast<double>(pot);
  std::cout << " Total POT: " << totPOT << ", constructed with: " << pot << " POT --> scale = " << potScale << std::endl;
  
  double total_escaping = one_escaping_ccinc+one_escaping_ncinc+one_escaping_nuoth;
  double total_diff_cut = track_length_diff_cut_ccinc+track_length_diff_cut_ncinc+track_length_diff_cut_nuoth;
  double total_chi2     = chi2_cuts_ccinc+chi2_cuts_ncinc+chi2_cuts_nuoth;

  ofstream fileTeX;
  fileTeX.open(stats_location+"preselection_cut_rates.tex");
  // TeX file header
  fileTeX << "\\begin{table}[h!] " << std::endl;
  fileTeX << "  \\centering " << std::endl;
  fileTeX << "  \\renewcommand{\\arraystretch}{1.4}" << std::endl;
  fileTeX << "  \\begin{tabular}{ m{3.1cm} * {2}{ >{\\centering\\arraybackslash}m{1.6cm} } *{1}{ >{\\centering\\arraybackslash}m{1.3cm} } * {2}{ >{\\centering\\arraybackslash}m{2.7cm} } }" << endl;
  fileTeX << "    \\hline" << endl;
  fileTeX << "    True topology~$\\rightarrow$ & \\multirow{2}{*}{$\\nu_{\\mu}$~CC} & \\multirow{2}{*}{$\\nu_{\\mu}$~NC} & \\multirow{2}{*}{$\\nu,\\bar{\\nu}$~Other} & \\multirow{2}{*}{$\\nu_{\\mu}$~CC~Pur.} & \\multirow{2}{*}{$\\nu_{\\mu}$~CC~Eff.} \\\\" << std::endl; 
  fileTeX << "    $\\downarrow$~Cut applied & & & \\\\" << std::endl;
  fileTeX << "    \\hline" << endl;

  fileTeX << "    Total & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_fiducial_cuts_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_fiducial_cuts_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(total_fiducial_cuts_nuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << total_fiducial_cuts_ccinc / static_cast<double>(total_fiducial_cuts_ccinc+total_fiducial_cuts_ncinc+total_fiducial_cuts_nuoth)*100 <<  "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(total_fiducial_cuts_ccinc,total_fiducial_cuts_ccinc+total_fiducial_cuts_ncinc+total_fiducial_cuts_nuoth)*100 << "~\\% & ";
  fileTeX << std::fixed << std::setprecision(2) << total_fiducial_cuts_ccinc / static_cast<double>(total_fiducial_cuts_ccinc)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(total_fiducial_cuts_ccinc,total_fiducial_cuts_ccinc)*100 <<"~\\% \\\\ " << std::endl;
  
  fileTeX << "    \\hdashline" << endl;
  
  fileTeX << "    One escapes & ";
  fileTeX << "    \\num{ " << static_cast<int>(one_escaping_ccinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(one_escaping_ncinc*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>(one_escaping_nuoth*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << one_escaping_ccinc / static_cast<double>(total_escaping)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(one_escaping_ccinc,total_escaping)*100 <<"~\\% & ";
  fileTeX << std::fixed << std::setprecision(2) << one_escaping_ccinc / static_cast<double>(total_fiducial_cuts_ccinc)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(one_escaping_ccinc,total_fiducial_cuts_ccinc)*100 <<"~\\% \\\\ ";
  fileTeX << std::endl;
 
  fileTeX << "    + $\\Delta L / L_{\\text{longest}}$ & ";
  fileTeX << "    \\num{ " << static_cast<int>((track_length_diff_cut_ccinc+one_escaping_ccinc)*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>((track_length_diff_cut_ncinc+one_escaping_ncinc)*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>((track_length_diff_cut_nuoth+one_escaping_nuoth)*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << (track_length_diff_cut_ccinc+one_escaping_ccinc) / static_cast<double>(total_escaping+total_diff_cut)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(track_length_diff_cut_ccinc+one_escaping_ccinc,total_escaping+total_diff_cut)*100 <<"~\\% & ";
  fileTeX << std::fixed << std::setprecision(2) << (track_length_diff_cut_ccinc+one_escaping_ccinc) / static_cast<double>(total_fiducial_cuts_ccinc)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(track_length_diff_cut_ccinc+one_escaping_ccinc,total_fiducial_cuts_ccinc)*100 <<"~\\% \\\\ ";
  fileTeX << std::endl;
  
  fileTeX << "    + $\\chi^{2}$ \\& length & ";
  fileTeX << "    \\num{ " << static_cast<int>((chi2_cuts_ccinc+track_length_diff_cut_ccinc+one_escaping_ccinc)*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>((chi2_cuts_ncinc+track_length_diff_cut_ncinc+one_escaping_ncinc)*potScale) << "} & ";
  fileTeX << "    \\num{ " << static_cast<int>((chi2_cuts_nuoth+track_length_diff_cut_nuoth+one_escaping_nuoth)*potScale) << "} & ";
  fileTeX << std::fixed << std::setprecision(2) << (chi2_cuts_ccinc+track_length_diff_cut_ccinc+one_escaping_ccinc) / static_cast<double>(total_escaping+total_diff_cut+total_chi2)*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(chi2_cuts_ccinc+track_length_diff_cut_ccinc+one_escaping_ccinc,total_escaping+total_diff_cut+total_chi2)*100 << "~\\% & ";
  fileTeX << std::fixed << std::setprecision(2) << (chi2_cuts_ccinc+track_length_diff_cut_ccinc+one_escaping_ccinc) / static_cast<double>(total_fiducial_cuts_ccinc)*100  << "~$\\pm$~" << std::fixed << std::setprecision(2) << GeneralAnalysisHelper::RatioUncertainty(chi2_cuts_ccinc+track_length_diff_cut_ccinc+one_escaping_ccinc,total_fiducial_cuts_ccinc)*100 << "~\\% \\\\ ";
  fileTeX << std::endl;
  
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

// Set the errors on the efficiency plots given the numerator and denomenator
void SetHistErrors(TH1D *h, TH1D *den){

  // Make sure the binning is the same
  if(h->GetNbinsX() != den->GetNbinsX()){
    std::cerr << " Error: Number of bins in the effiency histogram does not match that of the denomenator: " << std::endl;
    std::cerr << "   Histogram: " << h->GetNbinsX() << ", denomenator: " << den->GetNbinsX() << std::endl;
    exit(1);  
  }

  // Classical approach doesn't work at the limits ( e = 0,1 )
  //
  // Instead use Bayesean Approach from: https://cds.cern.ch/record/1010669/files/0701199.pdf?version=1
  //
  for(int b = 1; b <= h->GetNbinsX(); ++b){
    double num = h->GetBinContent(b) * den->GetBinContent(b);
    double denom = den->GetBinContent(b);
    double VarA = ((num+1)*(num+2))/static_cast<double>((denom+2)*(denom+3));
    double VarB = TMath::Power(num+1,2)/static_cast<double>(TMath::Power(denom+2,2));
    double sig  = TMath::Sqrt(VarA - VarB);
    h->SetBinError(b,sig);
//    std::cout << " Selected: " << num << ", total: " << denom << std::endl;
//    std::cout << " Eff: " << h->GetBinContent(b) << " +/- " << sig << std::endl;
  }
}
