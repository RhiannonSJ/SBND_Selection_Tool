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

namespace selection{
 
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
 
  void EventSelectionTool::LoadEventList(const std::string &file_name, EventList &event_list, const int &file){
 
    TFile f(file_name.c_str());
    TTree *t_event    = (TTree*) f.Get("event_tree");
    TTree *t_particle = (TTree*) f.Get("particle_tree");
    TTree *t_track    = (TTree*) f.Get("recotrack_tree");
    TTree *t_shower   = (TTree*) f.Get("recoshower_tree");

    TBranch *b_event_id        = t_event->GetBranch("event_id");
    TBranch *b_time_now        = t_event->GetBranch("time_now");
    TBranch *b_r_vertex        = t_event->GetBranch("r_vertex");
    TBranch *b_t_vertex        = t_event->GetBranch("t_vertex");
    TBranch *b_t_interaction   = t_event->GetBranch("t_interaction");
    TBranch *b_t_scatter       = t_event->GetBranch("t_scatter");
    TBranch *b_t_iscc          = t_event->GetBranch("t_iscc");
    TBranch *b_t_nu_pdgcode    = t_event->GetBranch("t_nu_pdgcode");
    TBranch *b_t_charged_pions = t_event->GetBranch("t_charged_pions");
    TBranch *b_t_neutral_pions = t_event->GetBranch("t_neutral_pions");
    TBranch *b_t_vertex_energy = t_event->GetBranch("t_vertex_energy");
    
    unsigned int n_events = t_event->GetEntries();

    unsigned int start_tracks      = 0;
    unsigned int start_showers     = 0;
    unsigned int start_mcparticles = 0;

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
      iscc         = b_t_iscc->GetLeaf("t_iscc")->GetValue();
      neutrino_pdg = b_t_nu_pdgcode->GetLeaf("t_nu_pdgcode")->GetValue();
      pions_ch     = b_t_charged_pions->GetLeaf("t_charged_pions")->GetValue();
      pions_neu    = b_t_neutral_pions->GetLeaf("t_neutral_pions")->GetValue();
      neu_energy   = b_t_vertex_energy->GetLeaf("t_vertex_energy")->GetValue();
   
      std::pair<int,int> event_identification(event_id,time_now);

      EventSelectionTool::GetTrackList(start_tracks, t_track, event_identification, tracks);
      EventSelectionTool::GetShowerList(start_showers, t_shower, event_identification, showers);
      EventSelectionTool::GetMCParticleList(start_mcparticles, t_particle, event_identification, mcparticles);
      EventSelectionTool::GetRecoParticleFromTrackRaquel(tracks, recoparticles);
      EventSelectionTool::GetRecoParticleFromShower(showers, r_vertex, recoparticles);
      
      event_list.push_back(Event(mcparticles, recoparticles, interaction, scatter, neutrino_pdg, pions_ch, pions_neu, iscc, t_vertex, r_vertex, neu_energy, file, event_id));

      start_tracks      += tracks.size();
      start_showers     += showers.size();
      start_mcparticles += mcparticles.size();
    }
  }

  //------------------------------------------------------------------------------------------ 
 
  void EventSelectionTool::GetUniqueEventList(TTree *event_tree, UniqueEventIdList &unique_event_list){
  
    TBranch *b_event_id = event_tree->GetBranch("event_id");
    TBranch *b_time_now = event_tree->GetBranch("time_now");
    
    unsigned int n_events = event_tree->GetEntries();

    for(unsigned int i = 0; i < n_events; ++i){
   
      event_tree->GetEntry(i);

      int event_id_a = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now_a = b_time_now->GetLeaf("time_now")->GetValue();
  
      bool shouldAdd(true);
    
      for(unsigned int j = 0; j < n_events; ++j){
        event_tree->GetEntry(j);

        int event_id_b = b_event_id->GetLeaf("event_id")->GetValue();
        int time_now_b = b_time_now->GetLeaf("time_now")->GetValue();
        
        if (event_id_a == event_id_b && time_now_a == time_now_b && i != j){
          shouldAdd = false;
          break;
        }
      }
       
      if (shouldAdd)
        unique_event_list.push_back(std::pair<int, int>(event_id_a, time_now_a));
    }        
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetTrackList(unsigned int start, TTree *track_tree, const std::pair<int, int> &unique_event, TrackList &track_list){
   
    TBranch *b_event_id       = track_tree->GetBranch("event_id");
    TBranch *b_time_now       = track_tree->GetBranch("time_now");
    TBranch *b_id_charge      = track_tree->GetBranch("tr_id_charge");
    TBranch *b_id_energy      = track_tree->GetBranch("tr_id_energy");
    TBranch *b_id_hits        = track_tree->GetBranch("tr_id_hits");
    TBranch *b_n_hits         = track_tree->GetBranch("tr_n_hits");
    TBranch *b_vertex         = track_tree->GetBranch("tr_vertex");
    TBranch *b_end            = track_tree->GetBranch("tr_end");
    TBranch *b_pida           = track_tree->GetBranch("tr_pida");
    TBranch *b_chi2_mu        = track_tree->GetBranch("tr_chi2_mu");
    TBranch *b_chi2_pi        = track_tree->GetBranch("tr_chi2_pi");
    TBranch *b_chi2_pr        = track_tree->GetBranch("tr_chi2_pr");
    TBranch *b_chi2_ka        = track_tree->GetBranch("tr_chi2_ka");
    TBranch *b_length         = track_tree->GetBranch("tr_length");
    TBranch *b_kinetic_energy = track_tree->GetBranch("tr_kinetic_energy");
    
    unsigned int n_entries = track_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      track_tree->GetEntry(i);

      int event_id         = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now         = b_time_now->GetLeaf("time_now")->GetValue();
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;
      
      double temp_vertex[3];
      double temp_end[3];

      int id_charge        = b_id_charge->GetLeaf("tr_id_charge")->GetValue();
      int id_energy        = b_id_energy->GetLeaf("tr_id_energy")->GetValue();
      int id_hits          = b_id_hits->GetLeaf("tr_id_hits")->GetValue();
      int n_hits           = b_n_hits->GetLeaf("tr_n_hits")->GetValue();
      temp_vertex[0]       = b_vertex->GetLeaf("tr_vertex")->GetValue(0);
      temp_vertex[1]       = b_vertex->GetLeaf("tr_vertex")->GetValue(1);
      temp_vertex[2]       = b_vertex->GetLeaf("tr_vertex")->GetValue(2);
      temp_end[0]          = b_end->GetLeaf("tr_end")->GetValue(0);
      temp_end[1]          = b_end->GetLeaf("tr_end")->GetValue(1);
      temp_end[2]          = b_end->GetLeaf("tr_end")->GetValue(2);
      float pida           = b_pida->GetLeaf("tr_pida")->GetValue();
      float chi2_mu        = b_chi2_mu->GetLeaf("tr_chi2_mu")->GetValue();
      float chi2_pi        = b_chi2_pi->GetLeaf("tr_chi2_pi")->GetValue();
      float chi2_pr        = b_chi2_pr->GetLeaf("tr_chi2_pr")->GetValue();
      float chi2_ka        = b_chi2_ka->GetLeaf("tr_chi2_ka")->GetValue();
      float length         = b_length->GetLeaf("tr_length")->GetValue();
      float kinetic_energy = b_kinetic_energy->GetLeaf("tr_kinetic_energy")->GetValue();

      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      
      float vertex_x = temp_vertex[0];                        
      float vertex_y = temp_vertex[1];                        
      float vertex_z = temp_vertex[2];                        
      float end_x    = temp_end[0];                        
      float end_y    = temp_end[1];                        
      float end_z    = temp_end[2];                        
                                                                                   
      // Co-ordinate offset in cm
      int sbnd_half_length_x = 400;
      int sbnd_half_length_y = 400;
      int sbnd_half_length_z = 500;
      
      int sbnd_offset_x = 200;
      int sbnd_offset_y = 200;
      int sbnd_offset_z = 0;

      int sbnd_border_x = 10;
      int sbnd_border_y = 20;
      int sbnd_border_z = 10;

      bool contained = true;

      if (    (vertex_x > (sbnd_half_length_x - sbnd_offset_x - sbnd_border_x)) 
           || (vertex_x < (-sbnd_offset_x + sbnd_border_x))          
           || (vertex_y > (sbnd_half_length_y - sbnd_offset_y - sbnd_border_y)) 
           || (vertex_y < (-sbnd_offset_y + sbnd_border_y))          
           || (vertex_z > (sbnd_half_length_z - sbnd_offset_z - sbnd_border_z)) 
           || (vertex_z < (-sbnd_offset_z + sbnd_border_z))
           || (end_x    > (sbnd_half_length_x - sbnd_offset_x - sbnd_border_x)) 
           || (end_x    < (-sbnd_offset_x + sbnd_border_x))          
           || (end_y    > (sbnd_half_length_y - sbnd_offset_y - sbnd_border_y)) 
           || (end_y    < (-sbnd_offset_y + sbnd_border_y))          
           || (end_z    > (sbnd_half_length_z - sbnd_offset_z - sbnd_border_z)) 
           || (end_z    < (-sbnd_offset_z + sbnd_border_z))) contained = false; 

      track_list.push_back(Track(id_charge, id_energy, id_hits, n_hits, pida, chi2_mu, chi2_pi, chi2_pr, chi2_ka, length, kinetic_energy, vertex, end, contained));
    
    } 
  }
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetShowerList(unsigned int start, TTree *shower_tree, const std::pair<int, int> &unique_event, ShowerList &shower_list){
  
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
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;
      
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

      shower_list.push_back(Shower(n_hits, vertex, direction, open_angle, length, energy));
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetMCParticleList(unsigned int start, TTree *mcparticle_tree, const std::pair<int, int> &unique_event, ParticleList &mcparticle_list){
    
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
      
      if(event_id != unique_event.first || time_now != unique_event.second) continue;

      double temp_vertex[3];
      double temp_end[3];
      double temp_momentum[3];
      
      int id                = b_id->GetLeaf("p_id")->GetValue();
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

      mcparticle_list.push_back(Particle(id, pdgcode, statuscode, n_hits, mass, energy, vertex, end, momentum));
      
    }
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionTool::GetRecoParticleFromTrackRaquel(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
  
    // Get the longest track id
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

    // Check if exactly 1 track escapes and that it is the longest track
    bool only_longest_escapes = false;
    bool none_escape = true;
    bool longest_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(i == longest_track_id) longest_escapes = true;
      if(!trk.m_contained) n_escaping++;
    }
    if(n_escaping == 1 && longest_escapes) only_longest_escapes = true;
    if(n_escaping != 0) none_escape = false;

    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(only_longest_escapes){
        if(!track.m_contained)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
      }
      else if(!none_escape) continue;
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetMuonByChi2Proton(track) == 13 || (id == longest_track_id && always_longest) || !track.m_contained)
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
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
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
    
    for(unsigned int id = 0; id < mu_candidates.size(); ++id){
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
    } 
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionTool::GetRecoParticleFromTrackChi2P(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(!trk.m_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

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
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(exactly_one_escapes){
        if(!track.m_contained) 
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetMuonByChi2Proton(track) == 13 || (id == longest_track_id && always_longest) || !track.m_contained)
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
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
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
    
    for(unsigned int id = 0; id < mu_candidates.size(); ++id){
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
    } 
  }

  //------------------------------------------------------------------------------------------ 

  void EventSelectionTool::GetRecoParticleFromTrack1Escaping(const TrackList &track_list, ParticleList &recoparticle_list){

    // Assign ridiculously short length to initiate the longest track length
    float longest_track_length      = -std::numeric_limits<float>::max();
    unsigned int longest_track_id   =  std::numeric_limits<unsigned int>::max();
   
    // Check if exactly 1 track escapes
    bool exactly_one_escapes = false;
    unsigned int n_escaping = 0;
    for(unsigned int i = 0; i < track_list.size(); ++i){
      const Track &trk(track_list[i]);
      if(!trk.m_contained) n_escaping++;
    }
    if(n_escaping == 1) exactly_one_escapes = true;

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
    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
      const Track &track(track_list[id]);
      // If exactly one particle escapes, call it the muon
      // Then identify protons
      // Then everything else
      if(exactly_one_escapes){
        if(!track.m_contained) 
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 13, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
      }
      else{
        // If the Chi2 Proton hypothesis gives proton, call the track a proton
        // Otherwise, call it a muon candidate
        if(EventSelectionTool::GetProtonByChi2Proton(track) == 2212)
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, 2212, track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
        else if(EventSelectionTool::GetPdgByChi2MuonCandidate(track) == 13 || (id == longest_track_id && always_longest) || !track.m_contained)
          mu_candidates.push_back(id);
        else
          recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end));
      }
    }

    // If the muon was found by length, this will return
    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
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
      else throw 10;
    }

    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(muon.m_mc_id_charge, muon.m_mc_id_energy, muon.m_mc_id_hits, 13, muon.m_n_hits, muon.m_kinetic_energy, muon.m_length, muon.m_vertex, muon.m_end));
    
    for(unsigned int id = 0; id < mu_candidates.size(); ++id){
      if(id == muonID) continue;
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(track.m_mc_id_charge, track.m_mc_id_energy, track.m_mc_id_hits, EventSelectionTool::GetPdgByChi2(track), track.m_n_hits, track.m_kinetic_energy, track.m_length, track.m_vertex, track.m_end)); 
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list){
 
    // To start with, 
    //  If an event has 2 showers
    //  Push back a pi0 at that point
    //if(shower_list.size() == 2) recoparticle_list.push_back(Particle(111,reco_vertex,reco_vertex));
    
    // New method
    // Calculate distance of closest approach between 2 photon showers 

    // There must be at least 2 showers
    if(shower_list.size() <= 2) return;

    std::vector<unsigned int> used_photon;
    used_photon.clear();

    // Loop over showers
    for(unsigned int i = 0; i < shower_list.size(); ++i){

      // Vector to hold 'j' values of candidate photon to pair with current photon
      // Vector to hold to position at which the photons meet to call it the point at which
      // the pi0 decayed
      // Vector to hold the distance of closest approach, in order to minimise this
      std::vector<unsigned int> candidate_id_for_pair;
      std::vector<TVector3> decay_point;
      std::vector<float> c_distance;
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
        float energy_1  = shower_list[i].m_energy;
        float energy_2  = shower_list[j].m_energy;

        float cos_theta = (b / (std::sqrt(a)*std::sqrt(c)));  
        float inv_mass  = std::sqrt(2*energy_1*energy_2*(1-cos_theta));

        float pi0_mass       = 134.97; // MeV
        
        // If the lines are parallel
        if(denomenator == 0) continue;

        float m_closest = (b*e - c*d)/denomenator;
        float n_closest = (a*e - b*d)/denomenator;

        TVector3 d_closest  = link + ((b*e - c*d)*dir_1 - (a*e - b*d)*dir_2)*(1/denomenator);
        float mag_d_closest =  d_closest.Mag();

        TVector3 d_middle   = (start_1 + m_closest*dir_1) + 0.5*d_closest;

        // If the distance of closest approach is smaller than 15 cm call it a candidate
        if(mag_d_closest < 15) { 
          candidate_id_for_pair.push_back(j);
          decay_point.push_back(d_closest);
          c_distance.push_back(mag_d_closest);
          total_hits.push_back(shower_list[i].m_n_hits + shower_list[j].m_n_hits);
        }
      }
  
      if(candidate_id_for_pair.size() == 0) continue;
      if(candidate_id_for_pair.size() == 1){
      
        // If the location with respect to the neutrino vertex at which the photons 
        // were produced by the candidate pi0 is more than 15 cm continue
        if((reco_vertex - decay_point[0]).Mag() > 15) continue;

        // push back a pi0 corresponding to the two photons
        recoparticle_list.push_back(Particle(111, total_hits[0], reco_vertex,decay_point[0]));
      } 
      else{
      
        // Find the minimum distance
        std::vector<float>::iterator min = std::min_element(c_distance.begin(), c_distance.end());
        TVector3 best_decay_point = decay_point[std::distance(c_distance.begin(), min)];
        int best_total_hits       = total_hits[std::distance(c_distance.begin(), min)];
        recoparticle_list.push_back(Particle(111,best_total_hits, reco_vertex, best_decay_point));
        used_photon.push_back(candidate_id_for_pair[std::distance(c_distance.begin(), min)]);
      
      }
    }
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
  
  int EventSelectionTool::GetProtonByChi2Proton(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr < 80) return 2212;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetMuonByChi2Proton(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_pr > 80) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByChi2MuonCandidate(const Track &track){

    // Limit based on particle gun plots of muon and proton vs proton chi^2
    if(track.m_chi2_mu < 16) return 13;
    return std::numeric_limits<int>::max();
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByPIDA(const Track &track){

    // Muon
    if(track.m_pida >= 5 && track.m_pida < 9) return 13;

    //Pion
    if(track.m_pida >= 9 && track.m_pida < 13) return 211;
  
    //Kaon
    if(track.m_pida >= 13 && track.m_pida < 13.5) return 321;
  
    //Proton
    if(track.m_pida >= 13) return 2212;

    return std::numeric_limits<int>::max();
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByPIDAStrict(const Track &track){

    // Muon
    if(track.m_pida >= 7.5  && track.m_pida < 8) return 13;

    //Pion
    if(track.m_pida >= 8    && track.m_pida < 9) return 211;
  
    //Kaon
    if(track.m_pida >= 12.9 && track.m_pida < 13.5) return 321;
  
    //Proton
    if(track.m_pida >= 16.9 && track.m_pida < 17.4) return 2212;

    return std::numeric_limits<int>::max();
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Track::Track(const int mc_id_charge, const int mc_id_energy, const int mc_id_hits, const int n_hits, const float pida, const float chi2_mu, const float chi2_pi, const float chi2_pr, const float chi2_ka, const float length, const float kinetic_energy, const TVector3 &vertex, const TVector3 &end, const bool &contained) :
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
    m_vertex(vertex),
    m_end(end),
    m_contained(contained) {}

  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Shower::Shower(const int n_hits, const TVector3 &vertex, const TVector3 &direction, const float open_angle, const float length, const float energy) :
    m_n_hits(n_hits),
    m_vertex(vertex),
    m_direction(direction),
    m_open_angle(open_angle),
    m_length(length),
    m_energy(energy/3.) {}

} // namespace: selection
