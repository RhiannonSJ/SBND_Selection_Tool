#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include <stdexcept>
#include "TROOT.h"
#include "TMath.h"
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

using namespace selection;

// The plane object
class Plane{

  public  :
    // Constructor 
    Plane(const TVector3 &V, const TVector3 &A, const TVector3 &B);

    // Getters
    TVector3 GetV() const;
    TVector3 GetA() const;
    TVector3 GetB() const;
    TVector3 GetUnitA() const;
    TVector3 GetUnitB() const;
    TVector3 GetUnitN() const;
    float    GetAlpha() const;
    float    GetBeta() const;

  private :
    // Member variables
    // Caps    = vector
    // No Caps = unit vector
    TVector3 m_a, m_b, m_n, m_A, m_B, m_V;
    float m_alpha, m_beta;
};

// Distance between particle start point and the point at which it crosses an (infinite) plane
float GetDistanceFromParticleToPlane(const Plane &plane, const Particle &particle);
bool  CheckIfParticleIntersectsPlane(const Plane &plane, const Particle &particle);
bool  IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane);

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/";
  std::string plots_location = "../Output_Selection_Tool/plots/escaping_track/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap numu_signal_map  = GeneralAnalysisHelper::GetNuMuTopologyMap();
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  // Define the fiducial planes of the detector
  //  The order: +/-x, +/-y, +/-z
  //
  // To take into account the fiducial border define 
  //    Half with fiducial   = Half width  - width fiducial border 
  //    Half height fiducial = Half height - height fiducial border
  //    Half length fiducial = Half length - length fiducial border
  float w = 200; // half full detector width (x)
  float h = 200; // half full detector height (y)
  float l = 250; // half full detector length (z)
  float w_border = 10; // width fiducial border
  float h_border = 20; // height fiducial border
  float l_border = 10; // length fiducial border
  float fid_w = w - w_border; // fiducial half width
  float fid_h = h - h_border; // fiducial half height
  float fid_l = l - l_border; // fiducial half length

  // Define the ficudial planes of the detector
  std::vector<Plane> planes;
  planes.emplace_back(TVector3( fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0,  fid_l));
  planes.emplace_back(TVector3(-fid_w, 0, l),         TVector3(0, fid_h, 0), TVector3(0, 0, -fid_l));
  planes.emplace_back(TVector3(0,  fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0, -fid_l));
  planes.emplace_back(TVector3(0, -fid_h, l),         TVector3(fid_w, 0, 0), TVector3(0, 0,  fid_l));
  planes.emplace_back(TVector3(0, 0, l_border),       TVector3(fid_w, 0, 0), TVector3(0,  fid_h, 0));
  planes.emplace_back(TVector3(0, 0, (2*l)-l_border), TVector3(fid_w, 0, 0), TVector3(0, -fid_h, 0));

  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){

   // if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/hepstore/rjones/Samples/FNAL/290818_analysis_sample/11206561_"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/old_220518_ana_files/8110339_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    
    //std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
 
  /*
   *
   *      TESTING
   *
   */
  // Counter for total number of escaping tracks
  unsigned int max_one_escapes         = 0;
  unsigned int too_many_escape         = 0;
  unsigned int too_many_true_contained = 0;
  unsigned int too_many_true_escaping  = 0;

  for(const Event &e : events){
    //Counter for event-based track counting
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
      too_many_escape++;
      if(e.IsSBNDTrueFiducial()) too_many_true_contained++;
      else too_many_true_escaping++;
    }
    else{
      max_one_escapes++;
    }
  }

  // Plot histograms for events with only 1 escaping tracks
  //    Plot for both scenarios when the escaping track is the longest and isn't the longest
  //    - Length of the escaping particle
  //    - Distance from closest face
  //      - x, y, z
  //    Also get the percentage of time that the escaping track is the longest
  TH1D *h_length_muon            = new TH1D("h_length_muon","Length of the escaping track", 50, 0, 400);
  TH1D *h_length_not_muon        = new TH1D("h_length_not_muon","Length of the escaping track", 50, 0, 400);
  TH1D *h_length_longest         = new TH1D("h_length_longest","Length of the escaping track", 50, 0, 400);
  TH1D *h_length_not             = new TH1D("h_length_not","Length of the escaping track", 50, 0, 400);
  TH1D *h_dist_x_longest         = new TH1D("h_dist_x_longest","Distance of the neutrino vertex from the nearest X face", 50, 0, 200);
  TH1D *h_dist_x_not             = new TH1D("h_dist_x_not","Distance of the neutrino vertex from the nearest X face", 50, 0, 200);
  TH1D *h_dist_y_longest         = new TH1D("h_dist_y_longest","Distance of the neutrino vertex from the nearest Y face", 50, 0, 200);
  TH1D *h_dist_y_not             = new TH1D("h_dist_y_not","Distance of the neutrino vertex from the nearest Y face", 50, 0, 200);
  TH1D *h_dist_z_longest         = new TH1D("h_dist_z_longest","Distance of the neutrino vertex from the nearest Z face", 50, 0, 200);
  TH1D *h_dist_z_not             = new TH1D("h_dist_z_not","Distance of the neutrino vertex from the nearest Z face", 50, 0, 200);
  TH1D *h_closest_longest        = new TH1D("h_closest_longest","Distance of the neutrino vertex from the nearest fiducial border", 50, 0, 200);
  TH1D *h_closest_not            = new TH1D("h_closest_not","Distance of the neutrino vertex from the nearest fiducial border", 50, 0, 200);
  TH1D *h_closest_muon           = new TH1D("h_closest_muon","Distance of the neutrino vertex from the nearest fiducial border", 50, 0, 200);
  TH1D *h_closest_not_muon       = new TH1D("h_closest_not_muon","Distance of the neutrino vertex from the nearest fiducial border", 50, 0, 200);
  TH2D *h_length_dist_muon       = new TH2D("h_length_dist_muon", "Length of the escaping track vs. distance from closest border", 50,0,200,20,0,400);
  TH2D *h_length_dist_not        = new TH2D("h_length_dist_not", "Length of the escaping track vs. distance from closest border", 50,0,200,20,0,400);
  TH1D *h_escaping_dist_muon     = new TH1D("h_escaping_dist_muon","Distance of the neutrino vertex from the escaping fiducial border", 100, 0, 500);
  TH1D *h_escaping_dist_not_muon = new TH1D("h_escaping_dist_not_muon","Distance of the neutrino vertex from the escaping fiducial border", 100, 0, 500);

  // Tree for TMVA calculation
  bool signal;
  float length, chi2p, chi2mu, dist;
  TFile f("../Output_Selection_Tool/tmva/trees/TMVA_input.root","RECREATE");
  TTree *escaping_muon_pid = new TTree("escaping_muon_pid", "TTree to hold variables which might distinguish an escaping track as a muon");
  TTree *all_particle_pid  = new TTree("all_particle_pid",  "TTree to hold variables which might distinguish an escaping track as a muon");
  escaping_muon_pid->Branch("signal", &signal, "signal/O");
  escaping_muon_pid->Branch("length", &length, "length/F");
  escaping_muon_pid->Branch("chi2p",  &chi2p,  "chi2p/F");
  escaping_muon_pid->Branch("chi2mu", &chi2mu, "chi2mu/F");
  escaping_muon_pid->Branch("dist",   &dist,   "dist/F");

  /*
   *
   * LENGTH
   *
   */
  unsigned int longest_escapes = 0;
  unsigned int true_muon_escapes = 0;
  unsigned int longest_over_100_escapes = 0;
  unsigned int events_with_1_escaping_track = 0;
  unsigned int longest_escaping_true_muon = 0;

  double escaping_track_length = -1.;

  for(const Event &e : events){
    bool longest_track_escapes = false;
    bool muon_escapes = false;
    // Only look at events with 1 escaping track
    if(GeneralAnalysisHelper::NumberEscapingTracks(e) != 1) continue;
    events_with_1_escaping_track++;
    
    // Now plot some things
    // Find the ID of the longest track
    double max = -std::numeric_limits<double>::max();
    int max_id = -1;
    for(const Particle &p : e.GetRecoParticleList()){
      if(p.GetFromRecoTrack() && p.GetLength() > max){
        max_id = p.GetMCParticleIdHits();
        max = p.GetLength();
      }
    }
    if(max_id == -1) continue;

    // Plot
    for(const Particle &p : e.GetRecoParticleList()){
      if(p.GetFromRecoTrack() && p.GetMCParticleIdHits() == max_id && p.GetOneEndTrackContained()){
        h_length_longest->Fill(p.GetLength());
        longest_escapes++;
        longest_track_escapes = true;
        for(const Particle &mcp : e.GetMCParticleList()){
          if(mcp.GetMCId() == max_id && mcp.GetPdgCode() == 13)
            longest_escaping_true_muon++;
        }
        if(p.GetLength() > 100) longest_over_100_escapes++;
      }
      else if(p.GetFromRecoTrack() && p.GetOneEndTrackContained()) 
        h_length_not->Fill(p.GetLength());
          
      if(p.GetFromRecoTrack() && p.GetOneEndTrackContained()){
        for(const Particle &mcp : e.GetMCParticleList()){
          if(mcp.GetPdgCode() == 13 && mcp.GetOneEndTrackContained()){
            h_length_muon->Fill(p.GetLength());
            escaping_track_length = p.GetLength();
            true_muon_escapes++;
            muon_escapes = true;
          }
          else if(mcp.GetOneEndTrackContained()){
            h_length_not_muon->Fill(p.GetLength());
            escaping_track_length = p.GetLength();
          }
        }
      }
    }
    
    /*
     *
     *  VERTEX DISTANCE
     *
     */
    // Fiducial volume parameters 
    float min_fid_x = e.GetMinimumFiducialDimensions()[0];
    float min_fid_y = e.GetMinimumFiducialDimensions()[1];
    float min_fid_z = e.GetMinimumFiducialDimensions()[2];
    float max_fid_x = e.GetMaximumFiducialDimensions()[0];
    float max_fid_y = e.GetMaximumFiducialDimensions()[1];
    float max_fid_z = e.GetMaximumFiducialDimensions()[2];
    
    // Neutrino vertex
    float nu_vtx_x = e.GetRecoNuVertex()[0];
    float nu_vtx_y = e.GetRecoNuVertex()[1];
    float nu_vtx_z = e.GetRecoNuVertex()[2];

    // Find out which the neutrino vertex is closest to, minimum or maximum
    //    Call the minimum distance delta(x,y,z)
    float delta_x = std::min(std::abs(max_fid_x - nu_vtx_x), std::abs(nu_vtx_x - min_fid_x));
    float delta_y = std::min(std::abs(max_fid_y - nu_vtx_y), std::abs(nu_vtx_y - min_fid_y));
    float delta_z = std::min(std::abs(max_fid_z - nu_vtx_z), std::abs(nu_vtx_z - min_fid_z));

    float delta_min = std::min(delta_x, std::min(delta_y, delta_z));

    if(longest_track_escapes){
      h_dist_x_longest->Fill(delta_x);
      h_dist_y_longest->Fill(delta_y);
      h_dist_z_longest->Fill(delta_z);
      h_closest_longest->Fill(delta_min);
    }
    else{
      h_dist_x_not->Fill(delta_x);
      h_dist_y_not->Fill(delta_y);
      h_dist_z_not->Fill(delta_z);
      h_closest_not->Fill(delta_min);
    }
    if(muon_escapes){
      h_length_dist_muon->Fill(delta_min, escaping_track_length);
      h_closest_muon->Fill(delta_min);
    }
    else{
      h_length_dist_not->Fill(delta_min, escaping_track_length);
      h_closest_not_muon->Fill(delta_min);
    }

    /*
     *
     *      VERTEX DISTANCE FROM ESCAPING TRACK LOCATION
     *
     */
    
    // Find the location at which the escaping track leaves the TPC
    for(const Particle &p : e.GetRecoParticleList()){
      // Check that we are looking at the escaping track
      if(p.GetFromRecoTrack() && p.GetOneEndTrackContained()){
        float distance_to_intersection_point = -std::numeric_limits<float>::max();
        // Loop over the fiducial planes and find out which the escaping particle passed through
        unsigned int i = 0;
        for(const Plane &plane : planes){
          i++;
          if(!CheckIfParticleIntersectsPlane(plane, p)) continue;
          distance_to_intersection_point = GetDistanceFromParticleToPlane(plane,p);
          break;
        }
        if(distance_to_intersection_point < 0 ) {
          std::cout << " Start : " << p.GetVertex()[0] << ", " << p.GetVertex()[1] << ", " << p.GetVertex()[2] << std::endl;
          std::cout << " End   : " << p.GetEnd()[0]    << ", " << p.GetEnd()[1]    << ", " << p.GetEnd()[2]    << std::endl;
        }
        if(muon_escapes) {
          h_escaping_dist_muon->Fill(distance_to_intersection_point);

          // Fill the TMVA tree 
          signal = 1;
          length = p.GetLength();
          chi2p  = p.GetChi2P();
          chi2mu = p.GetChi2Mu();
          dist   = distance_to_intersection_point;
        }
        else {
          h_escaping_dist_not_muon->Fill(distance_to_intersection_point);

          // Fill the TMVA tree 
          signal = 0;
          length = p.GetLength();
          chi2p  = p.GetChi2P();
          chi2mu = p.GetChi2Mu();
          dist   = distance_to_intersection_point;
        }
        escaping_muon_pid->Fill();
      }
    }
  }

  escaping_muon_pid->Write();
  f.Close();

  /*
   *
   *      OUTPUTS
   *
   */
  // Length vs dist 2D plot
  TCanvas *c = new TCanvas("c","",800,600);

  h_length_dist_muon->GetXaxis()->SetTitle("Distance to closest fiducial border [cm]");
  h_length_dist_muon->GetYaxis()->SetTitle("Length of the escaping track [cm]");
  h_length_dist_muon->Draw("colz");
  h_length_dist_muon->SetStats(0);

  c->SaveAs((plots_location+"length_dist_muon.png").c_str());
  c->SaveAs((plots_location+"length_dist_muon.root").c_str());

  c->Clear();
  
  h_length_dist_not->GetXaxis()->SetTitle("Distance to closest fiducial border [cm]");
  h_length_dist_not->GetYaxis()->SetTitle("Length of the escaping track [cm]");
  h_length_dist_not->Draw("colz");
  h_length_dist_not->SetStats(0);

  c->SaveAs((plots_location+"length_dist_not.png").c_str());
  c->SaveAs((plots_location+"length_dist_not.root").c_str());

  c->Clear();

  TLegend *leg = new TLegend(0.43,0.68,0.88,0.88);
  leg->AddEntry(h_escaping_dist_muon, "The true muon escapes", "l");
  leg->AddEntry(h_escaping_dist_not_muon, "The true muon doesn't escape", "l");
  leg->SetLineWidth(0);

  h_escaping_dist_muon->SetLineColor(2);
  h_escaping_dist_not_muon->SetLineColor(4);
  
  double int_edist_mu = h_escaping_dist_muon->Integral();
  double int_edist_not_mu = h_escaping_dist_not_muon->Integral();
  h_escaping_dist_muon->Scale(1./int_edist_mu);
  h_escaping_dist_not_muon->Scale(1./int_edist_not_mu);
 
  double maxedistmu = 1.1*std::max(h_escaping_dist_muon->GetMaximum(), h_escaping_dist_not_muon->GetMaximum());
  h_escaping_dist_muon->GetYaxis()->SetRangeUser(0,maxedistmu);
  h_escaping_dist_not_muon->GetYaxis()->SetRangeUser(0,maxedistmu);
  h_escaping_dist_muon->GetXaxis()->SetTitle("Distance to escaping fiducial border [cm]");
  h_escaping_dist_muon->Draw("hist");
  h_escaping_dist_not_muon->Draw("hist same");
  leg->Draw("same");
  h_escaping_dist_muon->SetStats(0);
  
  c->SaveAs((plots_location+"escaping_dist_muon.png").c_str());
  c->SaveAs((plots_location+"escaping_dist_muon.root").c_str());

  leg->Clear();
  c->Clear();
  
  leg->AddEntry(h_closest_muon, "The true muon escapes", "l");
  leg->AddEntry(h_closest_not_muon, "The true muon doesn't escape", "l");
  
  h_closest_muon->SetLineColor(2);
  h_closest_not_muon->SetLineColor(4);

  double int_dist_mu = h_closest_muon->Integral();
  double int_dist_not_mu = h_closest_not_muon->Integral();
  h_closest_muon->Scale(1./int_dist_mu);
  h_closest_not_muon->Scale(1./int_dist_not_mu);
 
  double maxdistmu = 1.1*std::max(h_closest_muon->GetMaximum(), h_closest_not_muon->GetMaximum());
  h_closest_muon->GetYaxis()->SetRangeUser(0,maxdistmu);
  h_closest_not_muon->GetYaxis()->SetRangeUser(0,maxdistmu);
  h_closest_muon->GetXaxis()->SetTitle("Distance to closest fiducial border [cm]");
  h_closest_muon->Draw("hist");
  h_closest_not_muon->Draw("hist same");
  leg->Draw("same");
  h_closest_muon->SetStats(0);

  c->SaveAs((plots_location+"closest_fid_distance_of_escaping_muon.png").c_str());
  c->SaveAs((plots_location+"closest_fid_distance_of_escaping_muon.root").c_str());

  leg->Clear();
  c->Clear();
  
  leg->AddEntry(h_length_muon, "The true muon escapes", "l");
  leg->AddEntry(h_length_not_muon, "The true muon doesn't escape", "l");
  
  h_length_muon->SetLineColor(2);
  h_length_not_muon->SetLineColor(4);

  double int_mu = h_length_muon->Integral();
  double int_not_mu = h_length_not_muon->Integral();
  h_length_muon->Scale(1./int_mu);
  h_length_not_muon->Scale(1./int_not_mu);
  
  double maxmu = 1.1*std::max(h_length_muon->GetMaximum(), h_length_not_muon->GetMaximum());
  h_length_muon->GetYaxis()->SetRangeUser(0.1,maxmu);
  h_length_not_muon->GetYaxis()->SetRangeUser(0.1,maxmu);
  h_length_muon->GetXaxis()->SetTitle("Length [cm]");
  h_length_muon->Draw("hist");
  h_length_not_muon->Draw("hist same");
  leg->Draw("same");
  h_length_muon->SetStats(0);

  c->SaveAs((plots_location+"length_of_escaping_muon.png").c_str());
  c->SaveAs((plots_location+"length_of_escaping_muon.root").c_str());

  leg->Clear();
  c->Clear();

  h_length_longest->SetLineColor(2);
  h_length_not->SetLineColor(4);

  leg->AddEntry(h_length_longest, "The longest track escapes", "l");
  leg->AddEntry(h_length_not, "The longest track doesn't escape", "l");
  leg->SetLineWidth(0);

  double int_long     = h_length_longest->Integral();
  double int_not_long = h_length_not->Integral();
  h_length_longest->Scale(1./int_long);
  h_length_not->Scale(1./int_not_long);
  
  double max = 1.1*std::max(h_length_longest->GetMaximum(), h_length_not->GetMaximum());
  h_length_longest->GetYaxis()->SetRangeUser(0,max);
  h_length_not->GetYaxis()->SetRangeUser(0,max);
  h_length_longest->GetXaxis()->SetTitle("Length [cm]");
  h_length_longest->Draw("hist");
  h_length_not->Draw("hist same");
  leg->Draw("same");
  h_length_longest->SetStats(0);

  c->SaveAs((plots_location+"length_of_escaping_track.png").c_str());
  c->SaveAs((plots_location+"length_of_escaping_track.root").c_str());

  c->Clear();

  h_dist_x_longest->SetLineColor(2);
  h_dist_x_not->SetLineColor(4);

  double int_dist_x     = h_dist_x_longest->Integral();
  double int_dist_x_not = h_dist_x_not->Integral();
  h_dist_x_longest->Scale(1./int_dist_x);
  h_dist_x_not->Scale(1./int_dist_x_not);
  
  double maxx = 1.1*std::max(h_dist_x_longest->GetMaximum(), h_dist_x_not->GetMaximum());
  h_dist_x_longest->GetYaxis()->SetRangeUser(0,maxx);
  h_dist_x_not->GetYaxis()->SetRangeUser(0,maxx);
  h_dist_x_longest->GetXaxis()->SetTitle("#Delta X [cm]");
  h_dist_x_longest->Draw();
  h_dist_x_not->Draw("same");
  leg->Draw("same");
  h_dist_x_longest->SetStats(0);

  c->SaveAs((plots_location+"delta_x.png").c_str());
  c->SaveAs((plots_location+"delta_x.root").c_str());

  c->Clear();

  h_dist_y_longest->SetLineColor(2);
  h_dist_y_not->SetLineColor(4);

  double int_dist_y     = h_dist_y_longest->Integral();
  double int_dist_y_not = h_dist_y_not->Integral();
  h_dist_y_longest->Scale(1./int_dist_y);
  h_dist_y_not->Scale(1./int_dist_y_not);
  
  double maxy = 1.1*std::max(h_dist_y_longest->GetMaximum(), h_dist_y_not->GetMaximum());
  h_dist_y_longest->GetYaxis()->SetRangeUser(0,maxy);
  h_dist_y_not->GetYaxis()->SetRangeUser(0,maxy);
  h_dist_y_longest->GetXaxis()->SetTitle("#Delta Y [cm]");
  h_dist_y_longest->Draw();
  h_dist_y_not->Draw("same");
  leg->Draw("same");
  h_dist_y_longest->SetStats(0);

  c->SaveAs((plots_location+"delta_y.png").c_str());
  c->SaveAs((plots_location+"delta_y.root").c_str());

  c->Clear();

  h_dist_z_longest->SetLineColor(2);
  h_dist_z_not->SetLineColor(4);

  double int_dist_z     = h_dist_z_longest->Integral();
  double int_dist_z_not = h_dist_z_not->Integral();
  h_dist_z_longest->Scale(1./int_dist_z);
  h_dist_z_not->Scale(1./int_dist_z_not);
  
  double maxz = 1.1*std::max(h_dist_z_longest->GetMaximum(), h_dist_z_not->GetMaximum());
  h_dist_z_longest->GetYaxis()->SetRangeUser(0,maxz);
  h_dist_z_not->GetYaxis()->SetRangeUser(0,maxz);
  h_dist_z_longest->GetXaxis()->SetTitle("#Delta Z [cm]");
  h_dist_z_longest->Draw();
  h_dist_z_not->Draw("same");
  leg->Draw("same");
  h_dist_z_longest->SetStats(0);

  c->SaveAs((plots_location+"delta_z.png").c_str());
  c->SaveAs((plots_location+"delta_z.root").c_str());

  c->Clear();

  h_closest_longest->SetLineColor(2);
  h_closest_not->SetLineColor(4);

  double maxclosest = 1.1*std::max(h_closest_longest->GetMaximum(), h_closest_not->GetMaximum());
  h_closest_longest->GetYaxis()->SetRangeUser(0,maxclosest);
  h_closest_not->GetYaxis()->SetRangeUser(0,maxclosest);
  h_closest_longest->GetXaxis()->SetTitle("Closest fiducial border distance [cm]");
  h_closest_longest->Draw();
  h_closest_not->Draw("same");
  leg->Draw("same");
  h_closest_longest->SetStats(0);

  c->SaveAs((plots_location+"delta_min.png").c_str());
  c->SaveAs((plots_location+"delta_min.root").c_str());

  c->Clear();

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"containment_testing.txt");

  file << "=====================================================================" << std::endl;
  file << " Total number of events                                             : " << events.size()   << std::endl;
  file << " Number of events with maximum one escaping track                   : " << max_one_escapes << std::endl;
  file << " Number of events with exactly one escaping track                   : " << events_with_1_escaping_track << std::endl;
  file << " Number of events with more than one escaping track                 : " << too_many_escape << std::endl;
  file << " Events with more than one escaping and true vertex contained       : " << too_many_true_contained << std::endl;
  file << " Events with more than one escaping and true vertex escaping        : " << too_many_true_escaping  << std::endl;
  file << " ------------------------------------------------------------------- " << std::endl;
  file << " For events with only 1 escaping track : " << std::endl;
  file << " ------------------------------------------------------------------- " << std::endl;
  file << " Percentage of events where the true muon escapes                   : " << 100 * true_muon_escapes/double(events_with_1_escaping_track) << std::endl;
  file << " Percentage of events where the longest track escapes               : " << 100 * longest_escapes/double(events_with_1_escaping_track) << std::endl;
  file << " Percentage of events with longest escaping and > 100 cm            : " << 100 * longest_over_100_escapes/double(events_with_1_escaping_track) << std::endl;
  file << " Percentage of events with longest escaping and longest true muon   : " << 100 * longest_escaping_true_muon/double(events_with_1_escaping_track) << std::endl;
  file << " Percentage of longest escaping track events with longest true muon : " << 100 * longest_escaping_true_muon/double(longest_escapes) << std::endl;
  file << "=====================================================================" << std::endl;

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

// Plane object constructor
Plane::Plane(const TVector3 &V, const TVector3 &A, const TVector3 &B) :
  m_A(A),
  m_B(B),
  m_V(V) {
    m_a = ((1/m_A.Mag()) * m_A);
    m_b = ((1/m_B.Mag()) * m_B);
    m_n = m_a.Cross(m_b);  
    m_alpha = m_A.Mag();
    m_beta  = m_B.Mag();

}

// Plane object getters
TVector3 Plane::GetV() const {return m_V;}
TVector3 Plane::GetA() const {return m_A;}
TVector3 Plane::GetB() const {return m_B;}
TVector3 Plane::GetUnitA() const {return m_a;}
TVector3 Plane::GetUnitB() const {return m_b;} 
TVector3 Plane::GetUnitN() const {return m_n;}
float Plane::GetAlpha() const {return m_alpha;}
float Plane::GetBeta()  const {return m_beta;}


// Get distance to plane 
float GetDistanceFromParticleToPlane(const Plane &plane, const Particle &particle){

  // Get the value of the unit vector of the particle dotted with the normal to the plane
  // If this is zero, the particle is parallel so throw an exception and catch it in the main
  // then continue looping through the list of planes
  TVector3 track_direction   = (particle.GetEnd() - particle.GetVertex()).Unit();
  float direction_from_plane = track_direction.Dot(plane.GetUnitN());

  if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()){
    throw std::domain_error("The particle is parallel to the plane");
  }

  return (1/direction_from_plane)*((plane.GetV() - particle.GetVertex()).Dot(plane.GetUnitN()));
}


bool CheckIfParticleIntersectsPlane(const Plane &plane, const Particle &particle){

  float d = -std::numeric_limits<float>::max();
  try{
    // Will throw exception if the particle is parallel to the plane
    d = GetDistanceFromParticleToPlane(plane, particle);
  }
  // If caught, the particle is parallel to the plane and does not intersect
  catch(const std::domain_error&) {return false;}

  if(d < 0 || d > particle.GetLength()) return false;
  TVector3 track_direction    = (particle.GetEnd() - particle.GetVertex()).Unit();
  TVector3 intersection_point = particle.GetVertex() + d * track_direction; 

  return IsProjectedPointInPlaneBounds(intersection_point, plane);
}

bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane){
  // Check if the point lies within the bound plane
  return (std::abs((point-plane.GetV()).Dot(plane.GetUnitA())) <= plane.GetAlpha() && std::abs((point-plane.GetV()).Dot(plane.GetUnitB())) <= plane.GetBeta());
}
