#include "../include/EventSelectionTool.h"
#include <iostream>
#include <numeric>
#include "TLeaf.h"
#include "TBranch.h"
#include "TVector3.h"
#include <algorithm>
#include <iterator>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <stdexcept>

namespace selection{

  double EventSelectionTool::GetPOT(TTree *subrun, const unsigned int det){
    // Get the pot for the individual file
    unsigned int n = subrun->GetEntries();
    double pot = 0.;
    for(unsigned int j = 0; j < n; ++j){
      if(det == 0){
        subrun->GetEntry(j);
        TBranch *b_pot = subrun->GetBranch("subrun_pot");
        pot += b_pot->GetLeaf("subrun_pot")->GetValue();
      }
      // There is a bug in the uboone and icarus POT so need to divide by 10
      else{
        subrun->GetEntry(j);
        TBranch *b_pot = subrun->GetBranch("subrun_pot");
        pot += b_pot->GetLeaf("subrun_pot")->GetValue()*0.1;
      }
    }
    return pot;
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionTool::GetTimeLeft(const int start_time, const int total, const unsigned int i){
   
    int now = static_cast<int>(time(NULL));
    int diff = now - start_time;
    int n_left = total - (i+1);
    int time_left = std::round(diff/double(i+1) * n_left);
    int seconds_per_minute = 60; 
    int seconds_per_hour = 3600; 
    int seconds_per_day = 86400;
    int days_left = time_left / (seconds_per_day);
    int hours_left = (time_left - (days_left * seconds_per_day)) / (seconds_per_hour);
    int minutes_left = (time_left - (days_left * seconds_per_day) - (hours_left * seconds_per_hour)) / (seconds_per_minute);
    int seconds_left = time_left - (days_left * seconds_per_day) - (hours_left * seconds_per_hour) - (minutes_left * seconds_per_minute);
    if(i%2){
      std::cout << " Estimated time left: " << std::setw(5) << days_left << " days, ";
      std::cout                             << std::setw(5) << hours_left << " hours, ";
      std::cout                             << std::setw(5) << minutes_left << " minutes, ";
      std::cout                             << std::setw(5) << seconds_left << " seconds." << '\r' << flush;
    }
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionTool::LoadEventList(const std::string &file_name, EventList &event_list, const int &file, double &pot, const Geometry &fid, const Geometry &av, const bool &runEverything ){
 
    TFile f(file_name.c_str());
    if(!f.IsOpen()){
      std::cerr << " Error opening file " << std::endl;
      exit(1);
    }

    TTree *t_event    = (TTree*) f.Get("selectionTree/event_tree");
    TTree *t_subrun   = (TTree*) f.Get("selectionTree/subrun_tree");
    TTree *t_particle = (TTree*) f.Get("selectionTree/particle_tree");
    TTree *t_track    = (TTree*) f.Get("selectionTree/recotrack_tree");
    TTree *t_shower   = (TTree*) f.Get("selectionTree/recoshower_tree");
   
    UIdToTrackListMap all_tracks;
    LoadTracks(t_track,all_tracks, av, runEverything);  

    UIdToShowerListMap all_showers;
    LoadShowers(t_shower,all_showers,t_track->GetEntries());  

    UIdToParticleListMap all_mcparticles;
    LoadMCParticles(t_particle,all_mcparticles,av);

    TBranch *b_event_id        = t_event->GetBranch("event_id");
    TBranch *b_time_now        = t_event->GetBranch("time_now");
    TBranch *b_r_vertex        = t_event->GetBranch("r_vertex");
    TBranch *b_t_vertex        = t_event->GetBranch("t_vertex");
    TBranch *b_t_interaction   = t_event->GetBranch("t_interaction");
    TBranch *b_t_scatter       = t_event->GetBranch("t_scatter");
    TBranch *b_t_baseline      = t_event->GetBranch("t_baseline");
    TBranch *b_t_iscc          = t_event->GetBranch("t_iscc");
    TBranch *b_t_nu_pdgcode    = t_event->GetBranch("t_nu_pdgcode");
    TBranch *b_t_charged_pions = t_event->GetBranch("t_charged_pions");
    TBranch *b_t_neutral_pions = t_event->GetBranch("t_neutral_pions");
    TBranch *b_t_vertex_energy = t_event->GetBranch("t_vertex_energy");
    TBranch *b_t_neutrino_qsqr = t_event->GetBranch("t_qsqr");
    TBranch *b_detector        = t_subrun->GetBranch("detector_enum");
    
    t_subrun->GetEntry(0);
    unsigned int det = b_detector->GetLeaf("detector_enum")->GetValue();
   
    pot = EventSelectionTool::GetPOT(t_subrun,det);
    
    unsigned int n_events = t_event->GetEntries();

    for(unsigned int j = 0; j < n_events; ++j){

      ParticleList mcparticles;
      ParticleList recoparticles;
      TrackList    tracks;
      ShowerList   showers;

      TVector3 r_vertex, t_vertex;
      unsigned int interaction, pions_ch, pions_neu, scatter;
      int neutrino_pdg;
      bool iscc(false);
      float neu_energy;
      float neu_qsqr;
      float baseline;

      t_event->GetEntry(j);

      int event_id = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now = b_time_now->GetLeaf("time_now")->GetValue();
      r_vertex[0]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(0);
      r_vertex[1]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(1);
      r_vertex[2]  = b_r_vertex->GetLeaf("r_vertex")->GetValue(2);
      t_vertex[0]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(0);
      t_vertex[1]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(1);
      t_vertex[2]  = b_t_vertex->GetLeaf("t_vertex")->GetValue(2);
      interaction  = b_t_interaction->GetLeaf("t_interaction")->GetValue();
      scatter      = b_t_scatter->GetLeaf("t_scatter")->GetValue();
      baseline     = b_t_baseline->GetLeaf("t_baseline")->GetValue();
      iscc         = b_t_iscc->GetLeaf("t_iscc")->GetValue();
      neutrino_pdg = b_t_nu_pdgcode->GetLeaf("t_nu_pdgcode")->GetValue();
      pions_ch     = b_t_charged_pions->GetLeaf("t_charged_pions")->GetValue();
      pions_neu    = b_t_neutral_pions->GetLeaf("t_neutral_pions")->GetValue();
      neu_energy   = b_t_vertex_energy->GetLeaf("t_vertex_energy")->GetValue();
      neu_qsqr     = b_t_neutrino_qsqr->GetLeaf("t_qsqr")->GetValue();
    
      const UniqueId event_identification(event_id,time_now);

      const auto track_itr = all_tracks.find(event_identification);
      if(track_itr != all_tracks.end())
        tracks = track_itr->second;

      const auto shower_itr = all_showers.find(event_identification);
      if(shower_itr != all_showers.end())
        showers = shower_itr->second;

      const auto mcparticle_itr = all_mcparticles.find(event_identification);
      if(mcparticle_itr != all_mcparticles.end())
        mcparticles = mcparticle_itr->second;
      
      if(tracks.size() != 0) {
//        EventSelectionTool::GetRecoParticleFromTrackSBND(tracks, recoparticles, av);
        if(det == 0)
          EventSelectionTool::GetRecoParticleFromTrackSBND(tracks, recoparticles, av);
        if(det == 1)
          EventSelectionTool::GetRecoParticleFromTrackMicroBooNE(tracks, recoparticles, av);
        if(det == 2)
          EventSelectionTool::GetRecoParticleFromTrackICARUS(tracks, recoparticles, av);
      }
      if(showers.size() != 0) EventSelectionTool::GetRecoParticleFromShower(showers, r_vertex, recoparticles, av);
     
      // Check if any particles should be flipped
      EventSelectionTool::CheckAndFlip(r_vertex, recoparticles);

      event_list.emplace_back(mcparticles, recoparticles, interaction, scatter, neutrino_pdg, pions_ch, pions_neu, iscc, t_vertex, r_vertex, neu_energy, neu_qsqr, file, event_id, baseline, fid);
    }
  }
  //------------------------------------------------------------------------------------------ 
  void EventSelectionTool::LoadMCParticles(TTree *mcparticle_tree, UIdToParticleListMap &mcparticles, const Geometry &g){
  
    TBranch *b_event_id = mcparticle_tree->GetBranch("event_id");
    TBranch *b_time_now = mcparticle_tree->GetBranch("time_now");
    TBranch *b_id       = mcparticle_tree->GetBranch("p_id");
    TBranch *b_n_hits   = mcparticle_tree->GetBranch("p_n_hits");
    TBranch *b_pdgcode  = mcparticle_tree->GetBranch("p_pdgcode");
    TBranch *b_status   = mcparticle_tree->GetBranch("p_status");
    TBranch *b_mass     = mcparticle_tree->GetBranch("p_mass");
    TBranch *b_energy   = mcparticle_tree->GetBranch("p_energy");
    TBranch *b_vertex   = mcparticle_tree->GetBranch("p_vertex");
    TBranch *b_end      = mcparticle_tree->GetBranch("p_end");
    TBranch *b_momentum = mcparticle_tree->GetBranch("p_momentum");
    
    unsigned int n_entries = mcparticle_tree->GetEntries();

      for(unsigned int i = 0; i < n_entries; ++i){
    
      mcparticle_tree->GetEntry(i);
      
      int event_id          = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now          = b_time_now->GetLeaf("time_now")->GetValue();
      const UniqueId id(event_id,time_now);
      
      double temp_vertex[3];
      double temp_end[3];
      double temp_momentum[3];
      
      int p_id              = b_id->GetLeaf("p_id")->GetValue();
      int pdgcode           = b_pdgcode->GetLeaf("p_pdgcode")->GetValue();
      int statuscode        = b_status->GetLeaf("p_status")->GetValue();
      int n_hits            = b_n_hits->GetLeaf("p_n_hits")->GetValue();
      float mass            = b_mass->GetLeaf("p_mass")->GetValue();
      float energy          = b_energy->GetLeaf("p_energy")->GetValue();
      temp_vertex[0]        = b_vertex->GetLeaf("p_vertex")->GetValue(0);
      temp_vertex[1]        = b_vertex->GetLeaf("p_vertex")->GetValue(1);
      temp_vertex[2]        = b_vertex->GetLeaf("p_vertex")->GetValue(2);
      temp_end[0]           = b_end->GetLeaf("p_end")->GetValue(0);
      temp_end[1]           = b_end->GetLeaf("p_end")->GetValue(1);
      temp_end[2]           = b_end->GetLeaf("p_end")->GetValue(2);
      temp_momentum[0]      = b_momentum->GetLeaf("p_momentum")->GetValue(0);
      temp_momentum[1]      = b_momentum->GetLeaf("p_momentum")->GetValue(1);
      temp_momentum[2]      = b_momentum->GetLeaf("p_momentum")->GetValue(2);
 
      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      TVector3 momentum(temp_momentum);

      mcparticles[id].emplace_back(p_id, pdgcode, statuscode, n_hits, mass, energy, vertex, end, momentum, g);
    }
  }
  //------------------------------------------------------------------------------------------ 
  void EventSelectionTool::LoadShowers(TTree *shower_tree, UIdToShowerListMap &showers, const int &n){
  
    TBranch *b_event_id   = shower_tree->GetBranch("event_id");
    TBranch *b_time_now   = shower_tree->GetBranch("time_now");
    TBranch *b_n_hits     = shower_tree->GetBranch("sh_n_hits");
    TBranch *b_vertex     = shower_tree->GetBranch("sh_start");
    TBranch *b_direction  = shower_tree->GetBranch("sh_direction");
    TBranch *b_open_angle = shower_tree->GetBranch("sh_open_angle");
    TBranch *b_length     = shower_tree->GetBranch("sh_length");
    TBranch *b_energy     = shower_tree->GetBranch("sh_energy");
    
    unsigned int n_entries = shower_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      shower_tree->GetEntry(i);

      int event_id      = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now      = b_time_now->GetLeaf("time_now")->GetValue();
      const UniqueId id(event_id,time_now);
      const int pid = n+i; // n = number of tracks
      
      double temp_vertex[3];
      double temp_direction[3];
      
      temp_vertex[0]    = b_vertex->GetLeaf("sh_start")->GetValue(0);
      temp_vertex[1]    = b_vertex->GetLeaf("sh_start")->GetValue(1);
      temp_vertex[2]    = b_vertex->GetLeaf("sh_start")->GetValue(2);
      temp_direction[0] = b_direction->GetLeaf("sh_direction")->GetValue(0);
      temp_direction[1] = b_direction->GetLeaf("sh_direction")->GetValue(1);
      temp_direction[2] = b_direction->GetLeaf("sh_direction")->GetValue(2);
      float open_angle  = b_open_angle->GetLeaf("sh_open_angle")->GetValue();
      float length      = b_length->GetLeaf("sh_length")->GetValue();
      float energy      = b_energy->GetLeaf("sh_energy")->GetValue();
      int n_hits        = b_n_hits->GetLeaf("sh_n_hits")->GetValue();
 
      TVector3 vertex(temp_vertex);
      TVector3 direction(temp_direction);

      showers[id].emplace_back(pid, n_hits, vertex, direction, open_angle, length, energy);
    } 
  }
  //------------------------------------------------------------------------------------------ 
  void EventSelectionTool::LoadTracks(TTree *track_tree, UIdToTrackListMap &tracks, const Geometry &g, const bool &runEverything ){

    // Get the outermost faces of the detector
    float min_x = *std::min_element(g.GetMinX().begin(),g.GetMinX().end());
    float min_y = *std::min_element(g.GetMinY().begin(),g.GetMinY().end());
    float min_z = *std::min_element(g.GetMinZ().begin(),g.GetMinZ().end());
    float max_x = *std::max_element(g.GetMaxX().begin(),g.GetMaxX().end());
    float max_y = *std::max_element(g.GetMaxY().begin(),g.GetMaxY().end());
    float max_z = *std::max_element(g.GetMaxZ().begin(),g.GetMaxZ().end());

    TBranch *b_event_id         = track_tree->GetBranch("event_id");
    TBranch *b_time_now         = track_tree->GetBranch("time_now");
    TBranch *b_id_charge        = track_tree->GetBranch("tr_id_charge");
    TBranch *b_id_energy        = track_tree->GetBranch("tr_id_energy");
    TBranch *b_id_hits          = track_tree->GetBranch("tr_id_hits");
    TBranch *b_n_hits           = track_tree->GetBranch("tr_n_hits");
    TBranch *b_vertex           = track_tree->GetBranch("tr_vertex");
    TBranch *b_end              = track_tree->GetBranch("tr_end");
    TBranch *b_pida             = track_tree->GetBranch("tr_pida");
    TBranch *b_chi2_mu          = track_tree->GetBranch("tr_chi2_mu");
    TBranch *b_chi2_pi          = track_tree->GetBranch("tr_chi2_pi");
    TBranch *b_chi2_pr          = track_tree->GetBranch("tr_chi2_pr");
    TBranch *b_chi2_ka          = track_tree->GetBranch("tr_chi2_ka");
    TBranch *b_length           = track_tree->GetBranch("tr_length");
    TBranch *b_kinetic_energy   = track_tree->GetBranch("tr_kinetic_energy");
    TBranch *b_mcs_mom_muon     = track_tree->GetBranch("tr_mcs_mom_muon");
    TBranch *b_range_mom_muon   = track_tree->GetBranch("tr_range_mom_muon");
    TBranch *b_range_mom_proton = track_tree->GetBranch("tr_range_mom_proton");
    TBranch *b_size             = track_tree->GetBranch("tr_dedx_size");      
    TBranch *b_residual_range;
    TBranch *b_dedx; 

    if(!runEverything){
      // Do not process the following branches
      track_tree->SetBranchStatus("tr_dedx",0);      
      track_tree->SetBranchStatus("tr_residual_range",0); 
    }
    else{
      b_residual_range = track_tree->GetBranch("tr_residual_range"); 
      b_dedx           = track_tree->GetBranch("tr_dedx");
    }
    
    unsigned int n_entries = track_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
      track_tree->GetEntry(i);

      int event_id         = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now         = b_time_now->GetLeaf("time_now")->GetValue();
      const UniqueId id(event_id,time_now);
      
      double temp_vertex[3];
      double temp_end[3];

      int n_hits             = b_n_hits->GetLeaf("tr_n_hits")->GetValue();
      int id_charge          = b_id_charge->GetLeaf("tr_id_charge")->GetValue();
      int id_energy          = b_id_energy->GetLeaf("tr_id_energy")->GetValue();
      int id_hits            = b_id_hits->GetLeaf("tr_id_hits")->GetValue();
      temp_vertex[0]         = b_vertex->GetLeaf("tr_vertex")->GetValue(0);
      temp_vertex[1]         = b_vertex->GetLeaf("tr_vertex")->GetValue(1);
      temp_vertex[2]         = b_vertex->GetLeaf("tr_vertex")->GetValue(2);
      temp_end[0]            = b_end->GetLeaf("tr_end")->GetValue(0);
      temp_end[1]            = b_end->GetLeaf("tr_end")->GetValue(1);
      temp_end[2]            = b_end->GetLeaf("tr_end")->GetValue(2);
      float pida             = b_pida->GetLeaf("tr_pida")->GetValue();
      float chi2_mu          = b_chi2_mu->GetLeaf("tr_chi2_mu")->GetValue();
      float chi2_pi          = b_chi2_pi->GetLeaf("tr_chi2_pi")->GetValue();
      float chi2_pr          = b_chi2_pr->GetLeaf("tr_chi2_pr")->GetValue();
      float chi2_ka          = b_chi2_ka->GetLeaf("tr_chi2_ka")->GetValue();
      float length           = b_length->GetLeaf("tr_length")->GetValue();
      float kinetic_energy   = b_kinetic_energy->GetLeaf("tr_kinetic_energy")->GetValue();
      float mcs_mom_muon     = b_mcs_mom_muon->GetLeaf("tr_mcs_mom_muon")->GetValue();
      float range_mom_muon   = b_range_mom_muon->GetLeaf("tr_range_mom_muon")->GetValue();
      float range_mom_proton = b_range_mom_proton->GetLeaf("tr_range_mom_proton")->GetValue();

      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      
      float vertex_x = temp_vertex[0];                        
      float vertex_y = temp_vertex[1];                        
      float vertex_z = temp_vertex[2];                        
      float end_x    = temp_end[0];                        
      float end_y    = temp_end[1];                        
      float end_z    = temp_end[2];                        
                                                                                   
      bool does_vtx_escape = 
        (vertex_x < min_x || vertex_x > max_x ||
         vertex_y < min_y || vertex_y > max_y ||
         vertex_z < min_z || vertex_z > max_z);
      bool does_end_escape = 
        (end_x    < min_x || end_x    > max_x ||
         end_y    < min_y || end_y    > max_y ||
         end_z    < min_z || end_z    > max_z);
      bool not_contained = does_vtx_escape || does_end_escape;

      bool one_end_escapes = true;
      if(does_vtx_escape && does_end_escape) one_end_escapes   = false;
      if(!does_vtx_escape && !does_end_escape) one_end_escapes = false;

      std::vector<float> dedx;
      std::vector<float> residual_range;

      if(runEverything){
        unsigned int n_dedx = b_size->GetLeaf("tr_dedx_size")->GetValue(); // Get the number of entries for the dedx & residual range branches
        for(unsigned int j = 0; j < n_dedx; ++j){
          dedx.push_back(b_dedx->GetLeaf("tr_dedx")->GetValue(j));
          residual_range.push_back(b_residual_range->GetLeaf("tr_residual_range")->GetValue(j));
        }
      }

      if(n_hits < 5) continue;
      tracks[id].emplace_back(i, id_charge, id_energy, id_hits, n_hits, pida, chi2_mu, chi2_pi, chi2_pr, chi2_ka, length, kinetic_energy, mcs_mom_muon, range_mom_muon, range_mom_proton,vertex, end, !not_contained, one_end_escapes, dedx, residual_range);
    
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  float EventSelectionTool::GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end){
    // Get the value of the unit vector of the particle dotted with the normal to the plane
    // If this is zero, the particle is parallel so throw an exception and catch it in the main
    // then continue looping through the list of planes
    TVector3 track_direction   = (end - vtx).Unit();
    float direction_from_plane = track_direction.Dot(plane.GetUnitN());

    if(std::abs(direction_from_plane) <= std::numeric_limits<float>::epsilon()){
      throw std::domain_error("The track is parallel to the plane");
    }

    return (1/direction_from_plane)*((plane.GetV() - vtx).Dot(plane.GetUnitN()));
  }

  //------------------------------------------------------------------------------------------ 

  float EventSelectionTool::GetDistanceFromParticleToPlane(const Plane &plane, const Particle &particle){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = particle.GetVertex();
    TVector3 end = particle.GetEnd();
    return EventSelectionTool::GetDistanceToPlane(plane, vtx, end);
  }

  //------------------------------------------------------------------------------------------ 

  float EventSelectionTool::GetDistanceFromTrackToPlane(const Plane &plane, const Track &track){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = track.m_vertex;
    TVector3 end = track.m_end;
    return EventSelectionTool::GetDistanceToPlane(plane, vtx, end);
  }

  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionTool::CheckIfParticleIntersectsPlane(const Plane &plane, const Particle &particle){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = particle.GetVertex();
    TVector3 end = particle.GetEnd();
    float length = particle.GetLength();
    return EventSelectionTool::CheckIfIntersectsPlane(plane, vtx, end, length);
  }
  
  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionTool::CheckIfTrackIntersectsPlane(const Plane &plane, const Track &track){
    // Get the vertex and end of the particle and pass to the distance function
    TVector3 vtx = track.m_vertex;
    TVector3 end = track.m_end;
    float length = track.m_length;
    return EventSelectionTool::CheckIfIntersectsPlane(plane, vtx, end, length);
  }

  //------------------------------------------------------------------------------------------ 
  
  bool EventSelectionTool::CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length){
    float d = -std::numeric_limits<float>::max();
    try{
      // Will throw exception if the particle is parallel to the plane
      d = GetDistanceToPlane(plane, vtx, end);
    }
    // If caught, the particle is parallel to the plane and does not intersect
    catch(const std::domain_error&) {return false;}

    if(d < 0 || d > length) return false;
    TVector3 track_direction    = (end - vtx).Unit();
    TVector3 intersection_point = vtx + d * track_direction; 

    return IsProjectedPointInPlaneBounds(intersection_point, plane);
  }

  //------------------------------------------------------------------------------------------ 

  bool EventSelectionTool::IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane){
    // Check if the point lies within the bound plane
    return (std::abs((point-plane.GetV()).Dot(plane.GetUnitA())) <= plane.GetAlpha() && std::abs((point-plane.GetV()).Dot(plane.GetUnitB())) <= plane.GetBeta());
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionTool::CheckAndFlip(const TVector3 &vtx, ParticleList &particles){
  
    // Loop over reconstructed particles
    // Check if the end point of the particle is closer to the neutrino vertex than the start
    // Flip if true
    for(Particle &p : particles){
      // Make sure the particle we are looking at is a reconstructed track
      if(!p.GetFromRecoTrack()) continue;
    
      // If the end is closer to the neutrino vertex than the start, flip it
      float nu_vtx_dist = (vtx - p.GetVertex()).Mag(); 
      float nu_end_dist = (vtx - p.GetEnd()).Mag(); 
      if(nu_vtx_dist > nu_end_dist) p.FlipTrack();
    }
  }

  //------------------------------------------------------------------------------------------ 
  /*
  void EventSelectionTool::GetRecoParticleFromTrack(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(trk.m_one_end_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

    // Check if the track is contained and passes the distance cut to the escaping border
    bool contained_and_passes_distance_cut = false;
    if(exactly_one_escapes){
      for(unsigned int i = 0; i < track_list.size(); ++i){
        const Track &trk(track_list[i]);
        // Only one end contained and the escaping track's length is greater than 100 cm
        if(trk.m_one_end_contained && trk.m_length >= 100){
          // Find out if the neutrino vertex is far enough from the escaping face
          float distance_to_intersection_point = -std::numeric_limits<float>::max();
          // Loop over the fiducial planes and find out which the escaping particle passed through
          PlaneList planes = g.GetExternalPlaneList();
          for(const Plane &plane : planes){
            if(!EventSelectionTool::CheckIfTrackIntersectsPlane(plane, trk)) continue;
            distance_to_intersection_point = EventSelectionTool::GetDistanceFromTrackToPlane(plane,trk);
            // Make sure the distance to the escaping border is big enough 
            // Run Contianment Main to determine this value from mu-pi comparison plot
            if(distance_to_intersection_point > 75){
              contained_and_passes_distance_cut = true;
              break;
            }
          }
        }
      }
    }
    
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &candidate(track_list[i]);
      // Get the lengths of the tracks and find the longest track and compare to the rest of
      // the lengths
      if(candidate.m_length > longest_track_length) {
        longest_track_length = candidate.m_length;
        longest_track_id     = i;
      }
    }

    bool always_longest(true);
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      // Find out if the longest track is always 1.5x longer than all the others in the event
      const Track &track(track_list[id]);
      if(track.m_length*1.5 >= longest_track_length && id != longest_track_id) always_longest = false;
    }
  
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;
    std::vector<int> used_ids;
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(contained_and_passes_distance_cut){
        // If one end is contained and the neutrino vertex is more than 75 cm from the escaping border
        if(track.m_one_end_contained) 
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        else if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        else
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        else if(EventSelectionTool::GetMuonByChi2Muon(track) == 13 || (id == longest_track_id && always_longest))
          mu_candidates.push_back(id);
        else
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
      }
    }
    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.emplace_back(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, muon.m_id, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range,g);
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_track_id == id && always_longest) {
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionTool::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else {
        std::cerr << "Haven't found the muon from the candidates" << std::endl;
        throw 10;
      }
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.emplace_back(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, muon.m_id, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range,g);
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range, g); 
    } 
  }*/
  
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromTrackSBND(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g){
    // Define any pre-calculated cuts
    double diff_cut      = 0.65;
    double longest_cut   = 100.;
    double length_cut    = 10.;
    double chi2p_cut     = 87.;
    double chi2mu_cut    = 14.;
    double chi2ratio_cut = 0.075;
    GetRecoParticleFromTrack(track_list, recoparticle_list, g, diff_cut, length_cut, longest_cut, chi2p_cut, chi2mu_cut, chi2ratio_cut, 0);
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromTrackMicroBooNE(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g){
    // Define any pre-calculated cuts
    double diff_cut      = 0.65;
    double longest_cut   = 100.;
    double length_cut    = 10.;
    double chi2p_cut     = 85.;
    double chi2mu_cut    = 19.;
    double chi2ratio_cut = 0.075;
    GetRecoParticleFromTrack(track_list, recoparticle_list, g, diff_cut, length_cut, longest_cut, chi2p_cut, chi2mu_cut, chi2ratio_cut, 1);
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromTrackICARUS(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g){
    // Define any pre-calculated cuts
    double diff_cut      = 0.65;
    double longest_cut   = 100.;
    double length_cut    = 10.;
    double chi2p_cut     = 20.;
    double chi2mu_cut    = 26.;
    double chi2ratio_cut = 1.;
    GetRecoParticleFromTrack(track_list, recoparticle_list, g, diff_cut, length_cut, longest_cut, chi2p_cut, chi2mu_cut, chi2ratio_cut, 2);
  }
  
  //------------------------------------------------------------------------------------------ 

  void EventSelectionTool::GetRecoParticleFromTrack(const TrackList &track_list, 
                                                    ParticleList &recoparticle_list, 
                                                    const Geometry &g,
                                                    const double &diff_cut,
                                                    const double &length_cut,
                                                    const double &longest_cut,
                                                    const double &chi2p_cut,
                                                    const double &chi2mu_cut,
                                                    const double &chi2ratio_cut,
                                                    const unsigned int &det){
    // Assign ridiculously short length to initiate the longest track length
    float longest  = -std::numeric_limits<float>::max();
    float second   = -std::numeric_limits<float>::max();
    int longest_id = -1;
    int second_id  = -1;
  
    // Check if exactly 1 track escapes
    bool one_long_escapes = false;
    unsigned int n_escaping = 0;
    int ntracks = track_list.size();
    for(int i = 0; i < ntracks; ++i){
      const Track &trk(track_list[i]);
      if(trk.m_one_end_contained && trk.m_length >= 100) n_escaping++;
    }
    if(n_escaping == 1) one_long_escapes = true;

    // Find the longest and second longest tracks
    for(int i = 0; i < ntracks; ++i){
      const Track &candidate(track_list[i]);
      // Find the longest track
      if(candidate.m_length > longest) {
        longest    = candidate.m_length;
        longest_id = i;
      }
    }
    for(int i = 0; i < ntracks; ++i){
      const Track &candidate(track_list[i]);
      // Find the second longest track
      if(candidate.m_length > second  && i != longest_id) {
        second    = candidate.m_length;
        second_id = i;
      }
    }
    // Get the fractional difference between the longest and second longest track lengths
    double diff = -999.;
    bool min_2_tracks = (longest_id != -1 && second_id != -1);
    if(min_2_tracks)
      diff = (longest - second)/longest;

    // Muon candidates 
    std::vector<unsigned int> mu_candidates;
    std::vector<int> used_ids;
    // Loop over track list
    for(int id = 0; id < ntracks; ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(one_long_escapes){
        // If one end is contained and the neutrino vertex is more than 75 cm from the escaping border
        if(track.m_one_end_contained) 
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 13, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        else if(EventSelectionTool::GetProtonByChi2Proton(track,det) == 2212 || track.m_length < length_cut)
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        else
          recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g); 
      }
      else{
        if(det == 0 || det == 1){
          // If the Chi2 Proton hypothesis gives proton, call the track a proton
          // Otherwise, call it a muon candidate
          if(EventSelectionTool::GetProtonByChi2Proton(track,det) == 2212 || track.m_length < length_cut)
            recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
          else if(EventSelectionTool::GetMuonByChi2Proton(track,det) == 13 || (id == longest_id && diff > diff_cut) || (id == longest_id && longest > longest_cut))
            mu_candidates.push_back(id);
          else
            recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        }
        else{
          if((id == longest_id && diff > diff_cut) || (id == longest_id && longest > longest_cut))
            mu_candidates.push_back(id);
          else if(EventSelectionTool::GetProtonByChi2Proton(track,det) == 2212 || track.m_length < length_cut)
            recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
          else
            recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range,g);
        
        }
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.emplace_back(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, muon.m_id, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range,g);
      return;
    }
    
    // If more than one muon candidate exists
    bool foundTheMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      int id = mu_candidates[i];
      const Track &candidate(track_list[id]);
      if(longest_id == id){
        muonID = id;
        foundTheMuon = true;
        break;
      }
    }
    if(!foundTheMuon) {
      // Find the smallest chi^2 under the muon hypothesis
      muonID = EventSelectionTool::GetMuonByChi2(track_list, mu_candidates);
      if(muonID != std::numeric_limits<unsigned int>::max()) foundTheMuon = true;
      else {
        std::cerr << "Haven't found the muon from the candidates" << std::endl;
        throw 10;
      }
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.emplace_back(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, muon.m_id, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_mcs_mom_muon, muon.m_range_mom_muon, muon.m_range_mom_proton,  muon.m_length, muon.m_vertex, muon.m_end, muon.m_chi2_pr, muon.m_chi2_mu, muon.m_chi2_pi, muon.m_dedx, muon.m_residual_range,g);
    
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.emplace_back(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, track.m_id, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_mcs_mom_muon, track.m_range_mom_muon, track.m_range_mom_proton,  track.m_length, track.m_vertex, track.m_end, track.m_chi2_pr, track.m_chi2_mu, track.m_chi2_pi, track.m_dedx, track.m_residual_range, g); 
    }
  }
  
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list, const Geometry &g){

    // New method
    // Calculate distance of closest approach between 2 photon showers 
    // There must be at least 2 showers to try and find a pi0

    std::vector<unsigned int> used_photon;
    ShowerList unused_showers;
    std::vector<unsigned int>::iterator it;
    unused_showers.clear();
    used_photon.clear();

    if(shower_list.size() == 1) unused_showers.push_back(shower_list[0]);
    else{
      // Loop over showers
      for(unsigned int i = 0; i < shower_list.size(); ++i){

        // Vector to hold 'j' values of candidate photon to pair with current photon
        // Vector to hold to position at which the photons meet to call it the point at which
        // the pi0 decayed
        // Vector to hold the distance of closest approach, in order to minimise this
        int id = shower_list[i].m_id;
        std::vector<unsigned int> candidate_id_for_pair;
        std::vector<TVector3> decay_point;
        std::vector<float> c_distance;
        std::vector<float> pi0_energy;
        std::vector<float> pi0_inv_mass_diff;
        std::vector<int> total_hits;
        candidate_id_for_pair.clear();
        decay_point.clear();
        c_distance.clear();
        total_hits.clear();

        for( unsigned int j = i+1; j < shower_list.size(); ++j){
          // If we are only looking at a single photon, continue
          if(i==j) continue;

          // If the photon has already been assigned to a pi0
          for(unsigned int k = 0; k < used_photon.size(); ++k) if(used_photon[k] == j) continue;

          // Get the distance of closest approach of the current two showers we are looking at
          TVector3 dir_1, dir_2, start_1, start_2, link;

          dir_1   = shower_list[i].m_direction;
          dir_2   = shower_list[j].m_direction;
          start_1 = shower_list[i].m_vertex;
          start_2 = shower_list[j].m_vertex;
          link    = start_1 - start_2;

          float a = dir_1.Dot(dir_1); 
          float b = dir_1.Dot(dir_2); 
          float c = dir_2.Dot(dir_2); 
          float d = dir_1.Dot(link); 
          float e = dir_2.Dot(link);
          float denomenator = a*c - b*b;

          // Get the invariant mass of the current 2 showers
          float energy_1   = shower_list[i].m_energy;
          float energy_2   = shower_list[j].m_energy;

          float cos_theta  = (b / (std::sqrt(a)*std::sqrt(c)));  
          float inv_mass   = std::sqrt(2*energy_1*energy_2*(1-cos_theta));
          float pi0_mass   = 134.97; // MeV

          float mass_diff  = std::abs(pi0_mass/1000. - inv_mass);

          // If the lines are parallel
          if(denomenator == 0) continue;

          float m_closest = (b*e - c*d)/denomenator;
          float n_closest = (a*e - b*d)/denomenator;

          TVector3 d_closest  = link + ((b*e - c*d)*dir_1 - (a*e - b*d)*dir_2)*(1/denomenator);
          float mag_d_closest =  d_closest.Mag();

          TVector3 d_middle   = (start_1 + m_closest*dir_1) + 0.5*d_closest;

          // If the distance of closest approach is smaller than 15 cm 
          // or the invariant mass of the photons is within 20% of the pion mass
          // call it a candidate
          if(mag_d_closest < 15 && mass_diff * (1./pi0_mass) < 0.2) { 
            candidate_id_for_pair.push_back(j);
            decay_point.push_back(d_closest);
            c_distance.push_back(mag_d_closest);
            total_hits.push_back(shower_list[i].m_n_hits + shower_list[j].m_n_hits);
            pi0_inv_mass_diff.push_back(mass_diff);
            pi0_energy.push_back(energy_1 + energy_2);
          }
        } // Inner shower list

        if(candidate_id_for_pair.size() == 0) continue;
        if(candidate_id_for_pair.size() == 1){
          // If the location with respect to the neutrino vertex at which the photons 
          // were produced by the candidate pi0 is more than 15 cm continue
          if((reco_vertex - decay_point[0]).Mag() > 15) continue;

          // push back a pi0 corresponding to the two photons
          recoparticle_list.emplace_back(111, id, total_hits[0], reco_vertex, decay_point[0], pi0_energy[0],g);
          used_photon.push_back(i);
          used_photon.push_back(candidate_id_for_pair[0]);

        } 
        else{
          // Find the minimum distance
          std::vector<float>::iterator min      = std::min_element(c_distance.begin(), c_distance.end());
          TVector3 best_decay_point             = decay_point[std::distance(c_distance.begin(), min)];
          int best_total_hits                   = total_hits[std::distance(c_distance.begin(), min)];
          float best_pi0_energy                 = pi0_energy[std::distance(c_distance.begin(), min)];
          recoparticle_list.emplace_back(111, id, best_total_hits, reco_vertex, best_decay_point, best_pi0_energy,g);
          used_photon.push_back(candidate_id_for_pair[std::distance(c_distance.begin(), min)]);
          used_photon.push_back(i);
        } // candidates
      } // Shower loop
      // Find any unused showers
      if(shower_list.size() > used_photon.size() && used_photon.size() != 0){
        for(unsigned int i = 0; i < shower_list.size(); ++i) {
          // If the current shower is not in the used photon list
          // Push it to unused showers
          it = std::find(used_photon.begin(), used_photon.end(), i);
          if(it == used_photon.end())
            unused_showers.push_back(shower_list[i]);
        }
      } // Used photons
    } // If more than 1 shower

    // For every unused shower, check the distance of the shower vertex from the neutrino vertex
    // If it is <= 15 cm, call it an electron
    // If it is > 15 cm, call it a photon 
    for(unsigned int i = 0; i < unused_showers.size(); ++i) {
      TVector3 shower_start = unused_showers[i].m_vertex;
      double conversion_length = (shower_start - reco_vertex).Mag();
      if(conversion_length <= 15) recoparticle_list.emplace_back(11, unused_showers[i].m_id, unused_showers[i].m_n_hits, unused_showers[i].m_vertex, (unused_showers[i].m_length * unused_showers[i].m_direction ), unused_showers[i].m_energy,g);
      else recoparticle_list.emplace_back(22, unused_showers[i].m_id, unused_showers[i].m_n_hits, unused_showers[i].m_vertex, (unused_showers[i].m_length * unused_showers[i].m_direction ), unused_showers[i].m_energy,g);
    } // Unused showers
  }
//------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByChi2(const Track &track){
    // Push the chi2 values onto a vector to find the minimum
    // NOT MUON
    std::map<int, float> chi2_map;

    chi2_map.insert(std::map<int, float>::value_type(211,  track.m_chi2_pi));
    chi2_map.insert(std::map<int, float>::value_type(321,  track.m_chi2_ka));
    chi2_map.insert(std::map<int, float>::value_type(2212, track.m_chi2_pr));

    float min_chi2 = std::numeric_limits<float>::max();
    int best_pdg   = std::numeric_limits<int>::max();

    for(std::map<int, float>::const_iterator it = chi2_map.begin(); it != chi2_map.end(); ++it){
    
      if(it->second < min_chi2){
        min_chi2 = it->second; 
        best_pdg = it->first;
      }
    } 
    return best_pdg;
  }

  //------------------------------------------------------------------------------------------ 
 
  int EventSelectionTool::GetMuonByChi2(const TrackList &tracks, const std::vector<unsigned int> &mu_candidates){

    // Loop over muon candidates and find smallest corresponding chi2_mu
    float min_chi2_mu = std::numeric_limits<float>::max();
    unsigned int min_chi2_id = std::numeric_limits<unsigned int>::max();

    // Find the smallest chi^2 under the muon hypothesis
    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
      unsigned int id = mu_candidates[i];
      const Track &candidate(tracks[id]);
      if(candidate.m_chi2_mu < min_chi2_mu) {
        min_chi2_mu = candidate.m_chi2_mu; 
        min_chi2_id = id;
      }
    }
    return min_chi2_id;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetProtonByChi2Proton(const Track &track, const unsigned int &det){
    double cut = 0.;
    if(det == 0 || det == 1)
      cut = 87.;
    else
      cut = 20.;

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr < cut) return 2212;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetMuonByChi2MuonProtonRatio(const Track &track, const unsigned int &det){
    double cut = 0.;
    if(det == 0 || det == 1)
      cut = 0.075;
    else
      cut = 1.;

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_mu/track.m_chi2_pr <= cut) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  int EventSelectionTool::GetMuonByChi2Proton(const Track &track, const unsigned int &det){
    double cut = 0.;
    if(det == 0)
      cut = 87.;
    else if(det == 1)
      cut = 85.;
    else
      cut = 20.;

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr >= cut) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetMuonByChi2Muon(const Track &track, const unsigned int &det){
    double cut = 0.;
    if(det == 0)
      cut = 14.;
    else if(det == 1)
      cut = 19.;
    else
      cut = 26.;

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_mu < cut) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Track::Track(const int id, const int mc_id_charge, const int mc_id_energy, const int mc_id_hits, const int n_hits, const float pida, const float chi2_mu, const float chi2_pi, const float chi2_pr, const float chi2_ka, const float length, const float kinetic_energy, const float mcs_momentum_muon, const float range_momentum_muon, const float range_momentum_proton, const TVector3 &vertex, const TVector3 &end, const bool &contained, const bool &one_end_contained, const std::vector<float> &dedx, const std::vector<float> &residual_range) :
    m_id(id),
    m_mc_id_charge(mc_id_charge),
    m_mc_id_energy(mc_id_energy),
    m_mc_id_hits(mc_id_hits),
    m_n_hits(n_hits),
    m_pida(pida),
    m_chi2_mu(chi2_mu),
    m_chi2_pi(chi2_pi),
    m_chi2_pr(chi2_pr),
    m_chi2_ka(chi2_ka),
    m_length(length),
    m_kinetic_energy(kinetic_energy),
    m_mcs_mom_muon(mcs_momentum_muon),
    m_range_mom_muon(range_momentum_muon),
    m_range_mom_proton(range_momentum_proton),
    m_vertex(vertex),
    m_end(end),
    m_contained(contained),
    m_one_end_contained(one_end_contained),
    m_dedx(dedx),
    m_residual_range(residual_range) {}

  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Shower::Shower(const int id, const int n_hits, const TVector3 &vertex, const TVector3 &direction, const float open_angle, const float length, const float energy) :
    m_id(id),
    m_n_hits(n_hits),
    m_vertex(vertex),
    m_direction(direction),
    m_open_angle(open_angle),
    m_length(length),
    m_energy(energy/1000.) {}

} // namespace: selection
//------------------------------------------------------------------------------------------ 
