#ifndef EVENT_SELECTION_TOOL_H
#define EVENT_SELECTION_TOOL_H

#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "Event.h"
#include "Particle.h"
#include "Plane.h"

namespace selection{
 
  /**
   * @brief  EventSelectionTool helper class
   */
  class EventSelectionTool {

    private : 
      class Track;
      class Shower;

    public : 

      typedef std::vector<Plane>                 PlaneList;
      typedef std::vector<Particle>              ParticleList;
      typedef std::vector<Event>                 EventList;
      typedef std::vector<Track>                 TrackList;
      typedef std::vector<Shower>                ShowerList;
      typedef std::pair<int,int>                 UniqueId;
      typedef std::map<UniqueId, TrackList>      UIdToTrackListMap;
      typedef std::map<UniqueId, ShowerList>     UIdToShowerListMap;
      typedef std::map<UniqueId, ParticleList>   UIdToParticleListMap;
      
      /**
       * @brief get the pot corresponding to each individual file
       *
       * @param  subrun tree with the some of the subrun information from larsoft
       * @param  det detector enumeration
       *
       * @return pot from the current file
       *
       */
      static double GetPOT(TTree *subrun, const unsigned int det);

      /**
       * @brief  load the list of events to analyse from the root file
       *
       * @param  file_name name of the root file to access
       * @param  event_list vector of events to fill
       * @param  file number of the file of the current event
       * @param  pot from the current file
       * @param  fid fiducial volume of this detector
       * @param  av active volume of this detector
       * @param  runEverything whether or not to get the really hefty branches
       *
       */
      static void LoadEventList(const std::string &file_name, EventList &event_list, const int &file, double &pot, const Geometry &fid, const Geometry &av, const bool &runEverything );

      /**
       * @brief  Output the length of time left in the running
       *
       * @param  start_time time at the start of the run
       * @param  total number of events
       * @param  iteration
       *
       */
      static void GetTimeLeft(const int start_time, const int total, const unsigned int i);
      
      /**
       * @brief  Load mcparticles
       *
       * @param  mcparticle_tree  
       * @param  map between unique event identifier and list of mcparticle objects
       * @param  g active volume of the detector
       *
       */
      static void LoadMCParticles(TTree *mcparticle_tree, UIdToParticleListMap &mcparticles, const Geometry &g);

      /**
       * @brief  Load showers
       *
       * @param  shower_tree  
       * @param  map between unique event identifier and list of shower objects
       * @param  n number of tracks to help construct unique particle id
       *
       */
      static void LoadShowers(TTree *shower_tree, UIdToShowerListMap &showers, const int &n);

      /**
       * @brief  Load tracks
       *
       * @param  track_tree  
       * @param  map between unique event identifier and list of track objects
       * @param  g active volume of the detector
       *
       */
      static void LoadTracks(TTree *track_tree, UIdToTrackListMap &tracks, const Geometry &g, const bool &runEverything );
      
      /**
       * @brief  The the distance between a particle and a defined plane
       *
       * @param  plane 
       * @param  particle
       * 
       * @return the distance between the particle and the plane
       */
      static float GetDistanceFromParticleToPlane(const Plane &plane, const Particle &particle);

      /**
       * @brief  The the distance between a particle and a defined plane
       *
       * @param  plane 
       * @param  track
       * 
       * @return the distance between the particle and the plane
       */
      static float GetDistanceFromTrackToPlane(const Plane &plane, const Track &track);

      /**
       * @brief  The the distance to a defined plane
       *
       * @param  plane 
       * @param  vtx
       * @param  end
       * 
       * @return the distance from the neutrino vertex to the plane
       */
      static float GetDistanceToPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end);
      
      /**
       * @brief  Check if the particle intersects a given plane
       *
       * @param  plane
       * @param  particle
       *
       * @return true or false
       */
      static bool CheckIfParticleIntersectsPlane(const Plane &plane, const Particle &particle);

      /**
       * @brief  Check if the track intersects a given plane
       *
       * @param  plane
       * @param  track
       *
       * @return true or false
       */
      static bool CheckIfTrackIntersectsPlane(const Plane &plane, const Track &track);

      /**
       * @brief  Check if it intersects a given plane
       *
       * @param  plane
       * @param  vtx
       * @param  end
       * @param  length
       *
       * @return true or false
       */
      static bool CheckIfIntersectsPlane(const Plane &plane, const TVector3 &vtx, const TVector3 &end, const float &length);

      /**
       * @brief  Check if the projected point is within the bounds of the given plane
       *
       * @param  point
       * @param  plane
       *
       * @return true or false
       */
      static bool IsProjectedPointInPlaneBounds(const TVector3 &point, const Plane &plane);
      
    private :

      /**
       * @brief  check if a particle should be flipped and flip it
       *
       * @param  vtx the neutrino interaction vertex
       * @param  particles the particles in the event
       *
       */
      static void CheckAndFlip(const TVector3 &vtx, ParticleList &particles);
      
      /**
       * @brief  get a list of reconstructed particles from track objects in old version
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       * @param  g active volume geometry of the outer detector
       * @param  diff_cut length difference cutbetween longest two tracks
       * @param  longest_cut longest track length cut
       * @param  chi2p_cut cut on chi2 proton
       * @param  chi2mu_cut cut on chi2 muon
       * @param  chi2ratio_cut cut on the ratio chi2mu/chi2p
       * @param  icarus whether we are calling the ICARUS PID or not
       *
       */
//      static void GetRecoParticleFromTrack(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g);
      static void GetRecoParticleFromTrack(const TrackList &track_list, 
                                           ParticleList &recoparticle_list, 
                                           const Geometry &g,
                                           const double &diff_cut,
                                           const double &length_cut,
                                           const double &longest_cut,
                                           const double &chi2p_cut,
                                           const double &chi2mu_cut,
                                           const double &chi2ratio_cut,
                                           const unsigned int &det);
      
      /**
       * @brief  get a list of reconstructed particles from track objects in SBND
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       * @param  g active volume geometry of the outer detector
       *
       */
      static void GetRecoParticleFromTrackSBND(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g);
      
      /**
       * @brief  get a list of reconstructed particles from track objects in MicroBooNE
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       * @param  g active volume geometry of the outer detector
       *
       */
      static void GetRecoParticleFromTrackMicroBooNE(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g);
      
      /**
       * @brief  get a list of reconstructed particles from track objects in ICARUS
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       * @param  g active volume geometry of the outer detector
       *
       */
      static void GetRecoParticleFromTrackICARUS(const TrackList &track_list, ParticleList &recoparticle_list, const Geometry &g);
      
      /**
       * @brief  get a list of reconstructed particles from track objects using original method
       *
       * @param  track_list list of tracks in the event
       * @param  recoparticle_list particle list to fill
       * @param  g active volume geometry of the outer detector
       *
       */
      static void GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list, const Geometry &g);

      /**
       * @brief  get the particle id based on its chi2 value
       *
       * @param  track the track to find the pdg of
       *
       * @return pdg
       *
       */
      static int GetPdgByChi2(const Track &track);
  
      /**
       * @brief  get the best muon candidate based on the smalleset chi2_mu value
       *
       * @param  tracks the track list to loop over
       * @param  mu_candidates the ids of the muon candidates from the track list
       *
       * @return best_id
       *
       */
      static int GetMuonByChi2(const TrackList &tracks, const std::vector<unsigned int> &mu_candidates);

      /**
       * @brief  get the whether the particle is a muon under the chi^2 proton hypothesis
       *
       * @param  track the track to find the pdg of
       *
       * @return pdg
       *
       */
      static int GetMuonByChi2Proton(const Track &track, const unsigned int &det);
      
      /**
       * @brief  get the whether the particle is a proton under the chi^2 proton hypothesis
       *
       * @param  track the track to find the pdg of
       * @param  detector
       *
       * @return pdg
       *
       */
      static int GetProtonByChi2Proton(const Track &track, const unsigned int &det);
      
      /**
       * @brief  get the whether the particle is a muon candidate from the chi2muon/chi2proton ratio
       *
       * @param  track the track to find the pdg of
       * @param  detector
       *
       * @return pdg
       *
       */
      static int GetMuonByChi2MuonProtonRatio(const Track &track, const unsigned int &det);

      /**
       * @brief  get the whether the particle is a muon candidate under the chi^2 muon hypothesis
       *
       * @param  detector
       * @param  track the track to find the pdg of
       *
       * @return pdg
       *
       */
      static int GetMuonByChi2Muon(const Track &track, const unsigned int &det);
      
      /**
       * @brief  Track class 
       */
      class Track{
      
        public : 
          
          /**
           * @brief  Constructor
           *
           * @param  id unique ID
           * @param  mc_id_charge mc TrackID corresponding to MCParticle using charge 
           * @param  mc_id_energy mc TrackID corresponding to MCParticle using energy
           * @param  mc_id_hits mc TrackID corresponding to MCParticle using hits
           * @param  pida pida value
           * @param  chi2_mu chi squared value for the muon fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_pi chi squared value for the pion fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_pr chi squared value for the proton fit of the reconstructed dEdx to the expected distribution
           * @param  chi2_ka chi squared value for the kaon fit of the reconstructed dEdx to the expected distribution
           * @param  length track length
           * @param  kinetic_energy track kinetic energy
           * @param  vertex vertex of the track
           * @param  end end point of the track
           * @param  contained whether or not the reconstructed track is contained within the  fiducial volume
           * @param  one_end_contained whether or not the reconstructed track has one end contained within the  fiducial volume
           *
           */
          Track(const int id, const int mc_id_charge, const int mc_id_energy, const int mc_id_hits, const int n_hits, const float pida, const float chi2_mu, const float chi2_pi, const float chi2_pr, const float chi2_ka, const float length, const float kinetic_energy, const float mcs_momentum_muon, const float range_momentum_muon, const float range_momentum_proton, const TVector3 &vertex, const TVector3 &end, const bool &contained, const bool &one_end_contained, const std::vector<float> &dedx, const std::vector<float> &residual_range);

          // Member variables
          int      m_id;                       ///< unique ID
          int      m_mc_id_charge;             ///< mc TrackID corresponding to MCParticle using charge
          int      m_mc_id_energy;             ///< mc TrackID corresponding to MCParticle using energy
          int      m_mc_id_hits;               ///< mc TrackID corresponding to MCParticle using hits
          int      m_n_hits;                   ///< number of hits in the track
          float    m_pida;                     ///< pida value
          float    m_chi2_mu;                  ///< chi squared fit to the muon expected dEdx
          float    m_chi2_pi;                  ///< chi squared fit to the pion expected dEdx
          float    m_chi2_pr;                  ///< chi squared fit to the proton expected dEdx 
          float    m_chi2_ka;                  ///< chi squared fit to the kaon expected dEdx
          float    m_length;                   ///< length of the track
          float    m_kinetic_energy;           ///< kinetic energy of the track
          float    m_mcs_mom_muon;             ///< multiple coulomb scattering momentum is the particle is an escaping muon
          float    m_range_mom_muon;           ///< range momentum if the particle is a contained muon 
          float    m_range_mom_proton;         ///< range momentum if the particle is a contained proton
          TVector3 m_vertex;                   ///< vertex of the track         
          TVector3 m_end;                      ///< end of the track
          bool     m_contained;                ///< whether or not the reconstructed track is contained
          bool     m_one_end_contained;        ///< whether or not the reconstructed track has one contained end
          std::vector<float> m_dedx;           ///< vector of the dedx distribution of the reconstructed track
          std::vector<float> m_residual_range; ///< vector of the residual range distribution of the track
      
      }; // Track
      
      /**
       * @brief  Shower class 
       */
      class Shower{
      
        public : 
          
          /**
           * @brief  Constructor
           *
           * @param  id unique ID
           * @param  vertex vertex of the shower
           * @param  direction direction of the shower
           * @param  open_angle opening angle at the vertex of the shower
           * @param  length length of the shower
           *
           */
          Shower(const int id, const int n_hits, const TVector3 &vertex, const TVector3 &direction, const float open_angle, const float length, const float energy);

          // Member variables
          int      m_id;         ///< unique ID
          int      m_n_hits;     ///< number of hits in the shower
          TVector3 m_vertex;     ///< vertex of the shower 
          TVector3 m_direction;  ///< direction of the shower
          float    m_open_angle; ///< opening angle at the vertex of the shower
          float    m_length;     ///< length of the shower
          float    m_energy;     ///< energy of the shower

      }; // Shower
  }; // EventSelectionTool
} // namespace: selection
#endif
