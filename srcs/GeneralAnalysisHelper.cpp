#include "../include/GeneralAnalysisHelper.h"

namespace selection{

  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNuMuTopologyMap() {
    TopologyMap signal_map_numu;
    signal_map_numu.insert(TopologyMap::value_type({14},1));
    return signal_map_numu;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNCTopologyMap() {
    TopologyMap signal_map_nc;
    signal_map_nc.insert(TopologyMap::value_type({13},0));
    return signal_map_nc;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNC0PiTopologyMap() {
    TopologyMap signal_map_nc0pi;
    signal_map_nc0pi.insert(TopologyMap::value_type({13},0));
    signal_map_nc0pi.insert(TopologyMap::value_type({211, -211, 111},0));
    return signal_map_nc0pi;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNC1PiTopologyMap() {
    TopologyMap signal_map_nc1pi;
    signal_map_nc1pi.insert(TopologyMap::value_type({13},0));
    signal_map_nc1pi.insert(TopologyMap::value_type({211, -211},1));
    signal_map_nc1pi.insert(TopologyMap::value_type({111},0));
    return signal_map_nc1pi;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNC2PiTopologyMap() {
    TopologyMap signal_map_nc2pi;
    signal_map_nc2pi.insert(TopologyMap::value_type({13},0));
    signal_map_nc2pi.insert(TopologyMap::value_type({211, -211},2));
    signal_map_nc2pi.insert(TopologyMap::value_type({111},0));
    return signal_map_nc2pi;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCCIncTopologyMap() {
    TopologyMap signal_map_cc_inc;
    signal_map_cc_inc.insert(TopologyMap::value_type({13},1));
    return signal_map_cc_inc;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0PiTopologyMap() {
    TopologyMap signal_map_cc_0pi;
    signal_map_cc_0pi.insert(TopologyMap::value_type({13},1));
    signal_map_cc_0pi.insert(TopologyMap::value_type({211, -211, 111},0));
    return signal_map_cc_0pi;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0Pi1PTopologyMap() {
    TopologyMap signal_map_cc_0pi_1p;
    signal_map_cc_0pi_1p.insert(TopologyMap::value_type({13},1));
    signal_map_cc_0pi_1p.insert(TopologyMap::value_type({211, -211, 111},0));
    signal_map_cc_0pi_1p.insert(TopologyMap::value_type({2212},1));
    return signal_map_cc_0pi_1p;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0Pi2PTopologyMap() {
    TopologyMap signal_map_cc_0pi_2p;
    signal_map_cc_0pi_2p.insert(TopologyMap::value_type({13},1));
    signal_map_cc_0pi_2p.insert(TopologyMap::value_type({211, -211, 111},0));
    signal_map_cc_0pi_2p.insert(TopologyMap::value_type({2212},2));
    return signal_map_cc_0pi_2p;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0Pi3PTopologyMap() {
    TopologyMap signal_map_cc_0pi_3p;
    signal_map_cc_0pi_3p.insert(TopologyMap::value_type({13},1));
    signal_map_cc_0pi_3p.insert(TopologyMap::value_type({211, -211, 111},0));
    signal_map_cc_0pi_3p.insert(TopologyMap::value_type({2212},3));
    return signal_map_cc_0pi_3p;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0Pi5PTopologyMap() {
    TopologyMap signal_map_cc_0pi_5p;
    signal_map_cc_0pi_5p.insert(TopologyMap::value_type({13},1));
    signal_map_cc_0pi_5p.insert(TopologyMap::value_type({211, -211, 111},0));
    signal_map_cc_0pi_5p.insert(TopologyMap::value_type({2212},5));
    return signal_map_cc_0pi_5p;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC1PiTopologyMap() { 
    TopologyMap signal_map_cc_1pi;
    signal_map_cc_1pi.insert(TopologyMap::value_type({13},1));
    signal_map_cc_1pi.insert(TopologyMap::value_type({211, -211},1));
    signal_map_cc_1pi.insert(TopologyMap::value_type({111},0));
    return signal_map_cc_1pi;
  }
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC2PiTopologyMap() { 
    TopologyMap signal_map_cc_2pi;
    signal_map_cc_2pi.insert(TopologyMap::value_type({13},1));
    signal_map_cc_2pi.insert(TopologyMap::value_type({211, -211},2));
    signal_map_cc_2pi.insert(TopologyMap::value_type({111},0));
    return signal_map_cc_2pi;
  }
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCCPi0TopologyMap() {
    TopologyMap signal_map_cc_pi0;
    signal_map_cc_pi0.insert(TopologyMap::value_type({13},1));
    signal_map_cc_pi0.insert(TopologyMap::value_type({111},1));
    signal_map_cc_pi0.insert(TopologyMap::value_type({211, -211},0));
    return signal_map_cc_pi0;
  } 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCCPi0PiCTopologyMap() {
    TopologyMap signal_map_cc_pi0_pic;
    signal_map_cc_pi0_pic.insert(TopologyMap::value_type({13},1));
    signal_map_cc_pi0_pic.insert(TopologyMap::value_type({111},1));
    signal_map_cc_pi0_pic.insert(TopologyMap::value_type({211, -211},1));
    return signal_map_cc_pi0_pic;
  } 


  //----------------------------------------------------------------------------------------
  //      DO NOT USE ON RECO 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNuETopologyMap() {
    TopologyMap signal_map_nue;
    signal_map_nue.insert(TopologyMap::value_type({12},1));
    return signal_map_nue;
  } 

  //----------------------------------------------------------------------------------------

  bool GeneralAnalysisHelper::MinOneRecoTrack(const Event &e){
    for(const Particle &p : e.GetRecoParticleList()){
      if(!p.GetFromRecoTrack()) continue;
      return true;
    }
    return false;
  }

  //----------------------------------------------------------------------------------------

  bool GeneralAnalysisHelper::PassedCCInclusive(const Event &e, const unsigned int &det){
    // Define cut variables per detector
    double diff_cut      = 0.7;
    double length_cut    = 10.;
    double longest_cut   = 50.;
    double escaping_cut  = 90.;
    double chi2p_cut     = 0.;
    double chi2mu_cut    = 0.;
    double chi2ratio_cut = 0.;
    if(det == 0){
      chi2p_cut     = 65.;
      chi2mu_cut    = 19.;
      chi2ratio_cut = 0.08;
    }
    if(det == 1){
      chi2p_cut     = 65.;
      chi2mu_cut    = 21.;
      chi2ratio_cut = 0.08;
      diff_cut      = 0.78;
    }
    if(det == 2){
      chi2p_cut     = 15.;
      chi2mu_cut    = 26.;
      chi2ratio_cut = 1.;
    }

    if(!e.IsRecoFiducial() || !MaxOneLongEscapingTrack(e) || !MinOneRecoTrack(e)) return false;

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

    // If 1 track escapes return true
    if(GeneralAnalysisHelper::NumberEscapingTracks(e) == 1)
      return true;

    // If no tracks escape
    else{
      // Not considering protons first
      bool track_cut_passed = false;
      for(Particle &p : e.GetRecoParticleList()){
        if(!p.GetFromRecoTrack()) continue;
        if(p.GetLength() > length_cut){
          track_cut_passed = true;
          break;
        }
      }
      if(track_cut_passed){
        // If we have two tracks and the longest two are very different in length, pass
        if(min_2_tracks && diff > diff_cut)
          return true;

        // Otherwise look at other variables
        else if((min_2_tracks && diff <= diff_cut) || !min_2_tracks){
          // Now consider chi2 variables
          for(const Particle &p : e.GetRecoParticleList()){
            if(!p.GetFromRecoTrack()) continue;
            //Check for clear protons
            if(det != 2){
              if((p.GetChi2Mu()/p.GetChi2P()) < chi2ratio_cut || (p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)
                return true;
            }
            else if((p.GetChi2P() > chi2p_cut && p.GetChi2Mu() < chi2mu_cut) || longest > longest_cut)
              return true;

          }
        }
      }
    }
    // Otherwise, we haven't passed any cuts - return false
    return false;
  }

  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::NumberEscapingTracks(const Event &e){
    unsigned int escaping_tracks = 0;
    int i = -1;
    for(const Particle &p : e.GetRecoParticleList()){
      i++;
      // Make sure the particle is a reconstructed track and check if it escapes
      if(p.GetFromRecoTrack() && p.GetOneEndTrackContained()) {
        escaping_tracks++;
      }
    }
    return escaping_tracks;
  }

  //----------------------------------------------------------------------------------------

  bool GeneralAnalysisHelper::MaxOneEscapingTrack(const Event &e){
    if(GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) return false;
    return true;
  }


  //----------------------------------------------------------------------------------------

  bool GeneralAnalysisHelper::MaxOneLongEscapingTrack(const Event &e){
    if(GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) return false;
    if(GeneralAnalysisHelper::NumberEscapingTracks(e) == 0) return true;
    double escaping_particle_length = 0;
    for(const Particle &p : e.GetRecoParticleList()){
      if(p.GetFromRecoTrack() && p.GetOneEndTrackContained()){
        escaping_particle_length = p.GetLength();
        break;
      }
    }
    if(escaping_particle_length >= 100) return true;
    return false;
  }

  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::TopologyStatistics(const Event &e, const TopologyMap signal_map_topology, double & count_true, double & count_signal, double & count_selected){
    if(e.CheckMCTopology(signal_map_topology)) count_true++;
    if(e.CheckRecoTopology(signal_map_topology)) count_selected++;
    if(e.CheckMCTopology(signal_map_topology) && e.CheckRecoTopology(signal_map_topology)) count_signal++;
  }

  //----------------------------------------------------------------------------------------
  
  ParticleMatrix GeneralAnalysisHelper::TopologyMatrix(const Event &e, ParticleMatrix &count_true_topology, ParticleMatrix &count_signal_topology, ParticleMatrix &count_selected_topology){
    std::vector<TopologyMap> topology_vector(5);
    topology_vector[0] = GeneralAnalysisHelper::GetNCTopologyMap();
    topology_vector[1] = GeneralAnalysisHelper::GetCCIncTopologyMap();
    topology_vector[2] = GeneralAnalysisHelper::GetCC0PiTopologyMap();
    topology_vector[3] = GeneralAnalysisHelper::GetCC1PiTopologyMap();
    topology_vector[4] = GeneralAnalysisHelper::GetCCPi0TopologyMap();

    for(unsigned int i=0; i < topology_vector.size(); ++i ){
      for(unsigned int j=0; j < topology_vector.size(); ++j ){
        if (e.CheckMCTopology(topology_vector[i])) count_true_topology[i][j]++;
        if (e.CheckRecoTopology(topology_vector[i])) count_selected_topology[i][j]++;
        if (e.CheckMCTopology(topology_vector[i]) && e.CheckRecoTopology(topology_vector[i])) count_signal_topology[i][j]++;
      }
    }
    return count_signal_topology;
  }
      
  //------------------------------------------------------------------------------------------ 
  
  bool GeneralAnalysisHelper::HasBeenReconstructed(const Event &e, const Particle &p){
    
    // Check if a true particle has a corresponding reconstructed particle
    // Only needs to happen once, may happen more than once but this is so that the user
    // can quickly check if the mc particle has been reconstructed
    //
    // Check that we are looking at an MC particle
    int true_id = p.GetMCId();
    ParticleList reco_particles = e.GetRecoParticleList();
    for(const Particle &p_reco : reco_particles){
      if(p_reco.GetFromRecoTrack() && GeneralAnalysisHelper::ParticleHasAMatch(e, p_reco) >= 0){
        if(GeneralAnalysisHelper::GetBestMCParticle(e,p_reco).GetMCId() == true_id) return true;
      }
    }
    return false;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int GeneralAnalysisHelper::ParticleHasAMatch(const Event &e, const Particle &p){
    // Starting from hits (since this is the chosen best method) find out if there is a match
    ParticleList particles = e.GetMCParticleList();
    for(const Particle &part : particles){
      if(part.GetMCId() == p.GetMCParticleIdHits())        return 0;
      else if(part.GetMCId() == p.GetMCParticleIdCharge()) return 1;
      else if(part.GetMCId() == p.GetMCParticleIdEnergy()) return 2;
    }
    return -1;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle GeneralAnalysisHelper::GetMCParticleCharge(const Event &e, const Particle &particle) {
    int charge_id = particle.GetMCParticleIdCharge();
    return GetMCParticle(charge_id, e.GetMCParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle GeneralAnalysisHelper::GetMCParticleEnergy(const Event &e, const Particle &particle) {
    int energy_id = particle.GetMCParticleIdEnergy();
    return GetMCParticle(energy_id, e.GetMCParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle GeneralAnalysisHelper::GetMCParticleHits(const Event &e, const Particle &particle) {
    int hits_id = particle.GetMCParticleIdHits();
    return GetMCParticle(hits_id, e.GetMCParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle GeneralAnalysisHelper::GetMCParticle(const int id, const ParticleList &particle_list) {
    for(const Particle &p : particle_list) {
      if(p.GetMCId() == id) return p;
    }
    std::cerr << "GetMCParticle" << std::endl;
    throw 8;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle GeneralAnalysisHelper::GetBestMCParticle(const Event &e, const Particle &particle) {
    /*
     * If the reconstructed particle has a match by
     *    0 == hits
     *    1 == charge
     *    2 == energy
     * get the MCParticle using the relevant method in order of preference (0 -> 2)
     */
    if(GeneralAnalysisHelper::ParticleHasAMatch(e, particle) == 0){
      return GeneralAnalysisHelper::GetMCParticleHits(e, particle);
    }
    else if(GeneralAnalysisHelper::ParticleHasAMatch(e, particle) == 1){
      return GeneralAnalysisHelper::GetMCParticleCharge(e, particle);
    }
    else if(GeneralAnalysisHelper::ParticleHasAMatch(e, particle) == 2){
      return GeneralAnalysisHelper::GetMCParticleEnergy(e, particle);
    }
    std::cerr << "GetBestMCParticle" << std::endl;
    throw 9;
  }
  
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg){
    // Check if the event is a selected event
    if(e.IsTrueFiducial() && e.CheckRecoTopology(topology)){
      return GeneralAnalysisHelper::CountMatchedParticles(e, e.GetRecoParticleList(), pdg);
    }
    else return 0;  
  }
      
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg){
    // Check if the event is a signal event
    if(e.IsTrueFiducial() && e.CheckRecoTopology(topology) && e.CheckMCTopology(topology)){
      return GeneralAnalysisHelper::CountMatchedParticles(e, e.GetRecoParticleList(), pdg);
    }
    else return 0;
  }
      
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesAll(const Event &e, const int pdg){
    return GeneralAnalysisHelper::CountMatchedParticles(e, e.GetRecoParticleList(), pdg);
  }
      
  //----------------------------------------------------------------------------------------
  
  unsigned int GeneralAnalysisHelper::CountMatchedParticles(const Event &e, const ParticleList &particle_list, const int pdg){
    unsigned int matched_particles = 0;
    // Loop over given particle list, if MatchedParticle, add to counter
    for(const Particle &p : particle_list){
      if(p.GetPdgCode() == pdg){
        if(GeneralAnalysisHelper::MatchedParticle(e,p)) matched_particles++;
      }
    }
    return matched_particles;
  }
  
  //----------------------------------------------------------------------------------------
  
  bool GeneralAnalysisHelper::MatchedParticle(const Event &e, const Particle &p){
    /**
     *
     * If the reconstructed particle has more than 5 hits and 
     *  If the reconstructed particle has been matched by
     *    0 == hits
     *    1 == charge
     *    2 == energy
     *    and the reconstructed pdgcode is the same at the truth pdgcode
     *    MATCHED PARTICLE == TRUE
     */
    if(p.GetNumberOfHits() >= 5){
      if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 0      && abs(GeneralAnalysisHelper::GetMCParticleHits(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
      else if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 1 && abs(GeneralAnalysisHelper::GetMCParticleCharge(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
      else if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 2 && abs(GeneralAnalysisHelper::GetMCParticleEnergy(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
      else return false;
    }
    else return false;
  }
  
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMCParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg){
    if(e.IsTrueFiducial() && e.CheckRecoTopology(topology)){
      return e.CountMCParticlesWithPdg(pdg);
    }
    return 0;
  }

  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMCParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg){
    if(e.CheckRecoTopology(topology) && e.CheckMCTopology(topology) && e.IsTrueFiducial()){
      return e.CountMCParticlesWithPdg(pdg);
    }
    return 0;
  }

  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountRecoParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg){
    if(e.IsTrueFiducial() && e.CheckRecoTopology(topology)){
      return e.CountRecoParticlesWithPdg(pdg);
    }
    return 0;
  }

  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountRecoParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg){
    if(e.CheckRecoTopology(topology) && e.CheckMCTopology(topology) && e.IsTrueFiducial()){
      return e.CountRecoParticlesWithPdg(pdg);
    }
    return 0;
  }
  
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int true_pdg, const int reco_pdg){
    if(e.IsTrueFiducial() && e.CheckRecoTopology(topology)){
      return GeneralAnalysisHelper::CountMisMatchedParticles(e, true_pdg, reco_pdg);
    }
    return 0;
  }
  
  //----------------------------------------------------------------------------------------
  
  unsigned int GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int true_pdg, const int reco_pdg){
    if(e.CheckRecoTopology(topology) && e.CheckMCTopology(topology) && e.IsTrueFiducial()){
      return GeneralAnalysisHelper::CountMisMatchedParticles(e, true_pdg, reco_pdg);
    }
    return 0;
  }
  
  //----------------------------------------------------------------------------------------
  
  unsigned int GeneralAnalysisHelper::CountMisMatchedParticles(const Event &e, const int true_pdg, const int reco_pdg){
    /*
     * Loop over given event's reconstructed particle list, if the current particle 
     * misidentification corresponds to 
     *    truth         = true_pdg
     *    reconstructed = reco_pdg
    * add to mismatched_particles counter
    *
    */
    unsigned int mismatched_particles = 0;
    ParticleList particles = e.GetRecoParticleList();
    for(const Particle &p : particles){
      if(p.GetPdgCode() == reco_pdg && GeneralAnalysisHelper::ParticleHasAMatch(e, p) >= 0){
        if(GeneralAnalysisHelper::GetBestMCParticle(e, p).GetPdgCode() == true_pdg && GeneralAnalysisHelper::GetBestMCParticle(e, p).GetNumberOfHits() >= 5 && p.GetNumberOfHits() >= 5) mismatched_particles++;
      }
    }
    return mismatched_particles;
  }

  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(const EventList &ev_list, const TopologyMap &topology, const std::string &topology_name, std::ofstream &os){

    // Counters for each particle type
    unsigned int mc_selected_muons     = 0;
    unsigned int mc_selected_protons   = 0;
    unsigned int mc_selected_pions     = 0;
    unsigned int mc_signal_muons       = 0;
    unsigned int mc_signal_protons     = 0;
    unsigned int mc_signal_pions       = 0;
    unsigned int reco_selected_muons   = 0;
    unsigned int reco_selected_protons = 0;
    unsigned int reco_selected_pions   = 0;
    unsigned int reco_signal_muons     = 0;
    unsigned int reco_signal_protons   = 0;
    unsigned int reco_signal_pions     = 0;
    unsigned int selected_muons        = 0;
    unsigned int selected_protons      = 0;
    unsigned int selected_pions        = 0;
    unsigned int signal_muons          = 0;
    unsigned int signal_protons        = 0;
    unsigned int signal_pions          = 0;

    for(const Event &e : ev_list){
      
      if(!e.IsTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) continue;

      mc_selected_muons     += GeneralAnalysisHelper::CountMCParticlesByTopologySelected(e, topology, 13);
      mc_selected_pions     += GeneralAnalysisHelper::CountMCParticlesByTopologySelected(e, topology, 211);
      mc_selected_pions     += GeneralAnalysisHelper::CountMCParticlesByTopologySelected(e, topology, -211);
      mc_selected_protons   += GeneralAnalysisHelper::CountMCParticlesByTopologySelected(e, topology, 2212);
      
      mc_signal_muons       += GeneralAnalysisHelper::CountMCParticlesByTopologySignal(e, topology, 13);
      mc_signal_pions       += GeneralAnalysisHelper::CountMCParticlesByTopologySignal(e, topology, 211);
      mc_signal_pions       += GeneralAnalysisHelper::CountMCParticlesByTopologySignal(e, topology, -211);
      mc_signal_protons     += GeneralAnalysisHelper::CountMCParticlesByTopologySignal(e, topology, 2212);

      reco_selected_muons   += GeneralAnalysisHelper::CountRecoParticlesByTopologySelected(e, topology, 13);
      reco_selected_pions   += GeneralAnalysisHelper::CountRecoParticlesByTopologySelected(e, topology, 211);
      reco_selected_pions   += GeneralAnalysisHelper::CountRecoParticlesByTopologySelected(e, topology, -211);
      reco_selected_protons += GeneralAnalysisHelper::CountRecoParticlesByTopologySelected(e, topology, 2212);
      
      reco_signal_muons     += GeneralAnalysisHelper::CountRecoParticlesByTopologySignal(e, topology, 13);
      reco_signal_pions     += GeneralAnalysisHelper::CountRecoParticlesByTopologySignal(e, topology, 211);
      reco_signal_pions     += GeneralAnalysisHelper::CountRecoParticlesByTopologySignal(e, topology, -211);
      reco_signal_protons   += GeneralAnalysisHelper::CountRecoParticlesByTopologySignal(e, topology, 2212);

      selected_muons        += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 13);
      selected_pions        += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 211);
      selected_pions        += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, -211);
      selected_protons      += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 2212);
                            
      signal_muons          += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 13);
      signal_pions          += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 211);
      signal_pions          += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, -211);
      signal_protons        += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 2212);
    }

    os << "    " << topology_name                                                                                               << std::endl;  
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << std::setw(16) << "Proton" << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Selected event, MC particles";
    os << std::setw(16) << mc_selected_muons;
    os << std::setw(16) << mc_selected_pions;
    os << std::setw(16) << mc_selected_protons;
    os << std::endl;
    os << std::setw(35) << "Selected event, Reco particles";
    os << std::setw(16) << reco_selected_muons;
    os << std::setw(16) << reco_selected_pions;
    os << std::setw(16) << reco_selected_protons;
    os << std::endl;
    os << std::setw(35) << "Selected event, Matched particles";
    os << std::setw(16) << selected_muons;
    os << std::setw(16) << selected_pions;
    os << std::setw(16) << selected_protons;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Signal event, MC particles";
    os << std::setw(16) << mc_signal_muons;
    os << std::setw(16) << mc_signal_pions;
    os << std::setw(16) << mc_signal_protons;
    os << std::endl;
    os << std::setw(35) << "Signal event, Reco particles";
    os << std::setw(16) << reco_signal_muons;
    os << std::setw(16) << reco_signal_pions;
    os << std::setw(16) << reco_signal_protons;
    os << std::endl;
    os << std::setw(35) << "Signal event, Matched particles";
    os << std::setw(16) << signal_muons;
    os << std::setw(16) << signal_pions;
    os << std::setw(16) << signal_protons;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Selected, Efficiency";
    os << std::setw(16) << std::setprecision(5) << 100 * selected_muons/double(mc_selected_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * selected_pions/double(mc_selected_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * selected_protons/double(mc_selected_protons);
    os << std::endl;
    os << std::setw(35) << "Selected, Purity";
    os << std::setw(16) << std::setprecision(5) << 100 * selected_muons/double(reco_selected_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * selected_pions/double(reco_selected_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * selected_protons/double(reco_selected_protons);
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Signal, Efficiency";
    os << std::setw(16) << std::setprecision(5) << 100 * signal_muons/double(mc_signal_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * signal_pions/double(mc_signal_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * signal_protons/double(mc_signal_protons);
    os << std::endl;
    os << std::setw(35) << "Signal, Purity";
    os << std::setw(16) << std::setprecision(5) << 100 * signal_muons/double(reco_signal_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * signal_pions/double(reco_signal_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * signal_protons/double(reco_signal_protons);
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
  }
  
  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(const EventList &ev_list, const TopologyMap &topology, const std::string &topology_name, std::ofstream &os){
  
    /* 
     * Counters for particle type and mis-identified counterpart
     *    left particle:  reconstructed
     *    right particle: true
     *
     * e.g. muon_charged_pion = charged pion has been misidentified as a muon
     *
     */
    unsigned int signal_muon        = 0;
    unsigned int signal_muon_pion   = 0;
    unsigned int signal_muon_proton = 0;
    unsigned int signal_pion        = 0;
    unsigned int signal_pion_muon   = 0;
    unsigned int signal_pion_proton = 0;
    unsigned int signal_proton      = 0;
    unsigned int signal_proton_muon = 0;
    unsigned int signal_proton_pion = 0;
    
    unsigned int selected_muon        = 0;
    unsigned int selected_muon_pion   = 0;
    unsigned int selected_muon_proton = 0;
    unsigned int selected_pion        = 0;
    unsigned int selected_pion_muon   = 0;
    unsigned int selected_pion_proton = 0;
    unsigned int selected_proton      = 0;
    unsigned int selected_proton_muon = 0;
    unsigned int selected_proton_pion = 0;
    
    for(const Event &e : ev_list){
      if(!e.IsTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) continue;
      
      selected_muon        += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 13, 13);
      selected_muon_pion   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 211, 13);
      selected_muon_pion   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, -211, 13);
      selected_muon_proton += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 2212, 13);
      
      signal_muon          += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 13, 13);
      signal_muon_pion     += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 211, 13);
      signal_muon_pion     += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, -211, 13);
      signal_muon_proton   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 2212, 13);
      
      selected_pion        += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 211, 211);
      selected_pion        += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, -211, 211);
      selected_pion_muon   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 13, 211);
      selected_pion_muon   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 13, -211);
      selected_pion_proton += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 2212, 211);
      selected_pion_proton += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 2212, -211);
      
      signal_pion          += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 211, 211);
      signal_pion          += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, -211, 211);
      signal_pion_muon     += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 13, 211);
      signal_pion_muon     += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 13, -211);
      signal_pion_proton   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 2212, 211);
      signal_pion_proton   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 2212, -211);
      
      selected_proton      += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 2212, 2212);
      selected_proton_muon += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 13, 2212);
      selected_proton_pion += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, 211, 2212);
      selected_proton_pion += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySelected(e, topology, -211, 2212);
      
      signal_proton        += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 2212, 2212);
      signal_proton_muon   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 13, 2212);
      signal_proton_pion   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, 211, 2212);
      signal_proton_pion   += GeneralAnalysisHelper::CountMisMatchedParticlesByTopologySignal(e, topology, -211, 2212);
    }

    os << "    " << topology_name                                                                                               << std::endl;  
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Selected event";
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "True/Reco" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << std::setw(16) << "Proton" << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Muon";
    os << std::setw(16) << selected_muon;
    os << std::setw(16) << selected_pion_muon;
    os << std::setw(16) << selected_proton_muon;
    os << std::endl;
    os << std::setw(35) << "Charged pion";
    os << std::setw(16) << selected_muon_pion;
    os << std::setw(16) << selected_pion;
    os << std::setw(16) << selected_proton_pion;
    os << std::endl;
    os << std::setw(35) << "Proton";
    os << std::setw(16) << selected_muon_proton;
    os << std::setw(16) << selected_pion_proton;
    os << std::setw(16) << selected_proton;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Signal event";
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "True/Reco" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << std::setw(16) << "Proton" << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Muon";
    os << std::setw(16) << signal_muon;
    os << std::setw(16) << signal_pion_muon;
    os << std::setw(16) << signal_proton_muon;
    os << std::endl;
    os << std::setw(35) << "Charged pion";
    os << std::setw(16) << signal_muon_pion;
    os << std::setw(16) << signal_pion;
    os << std::setw(16) << signal_proton_pion;
    os << std::endl;
    os << std::setw(35) << "Proton";
    os << std::setw(16) << signal_muon_proton;
    os << std::setw(16) << signal_pion_proton;
    os << std::setw(16) << signal_proton;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
  }

  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(const EventList &ev_list, std::ofstream &os){
   
    unsigned int mc_muons           = 0;
    unsigned int mc_protons         = 0;
    unsigned int mc_charged_pions   = 0;
    
    unsigned int reco_muons         = 0;
    unsigned int reco_protons       = 0;
    unsigned int reco_charged_pions = 0;
    
    unsigned int muons              = 0;
    unsigned int protons            = 0;
    unsigned int charged_pions      = 0;
    
    for(const Event &e : ev_list){
      if(!e.IsTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) continue;
      
      mc_muons           += e.CountMCParticlesWithPdg(13);
      mc_charged_pions   += e.CountMCParticlesWithPdg(211);
      mc_charged_pions   += e.CountMCParticlesWithPdg(-211);
      mc_protons         += e.CountMCParticlesWithPdg(2212);
      
      reco_muons         += e.CountRecoParticlesWithPdg(13);
      reco_charged_pions += e.CountRecoParticlesWithPdg(211);
      reco_charged_pions += e.CountRecoParticlesWithPdg(-211);
      reco_protons       += e.CountRecoParticlesWithPdg(2212);

      muons              += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 13);
      charged_pions      += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 211);
      charged_pions      += GeneralAnalysisHelper::CountMatchedParticlesAll(e, -211);
      protons            += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 2212);
    }

    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << "    All events"                                                                                               << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << std::setw(35) << "" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << "Proton" << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << std::setw(35) << "MC Particles";
    os << std::setw(16) << mc_muons;
    os << std::setw(16) << mc_charged_pions;
    os << std::setw(16) << mc_protons;
    os << std::endl;
    os << std::setw(35) << "Reco Particles";
    os << std::setw(16) << reco_muons;
    os << std::setw(16) << reco_charged_pions;
    os << std::setw(16) << reco_protons;
    os << std::endl;
    os << std::setw(35) << "Matched Particles";
    os << std::setw(16) << muons;
    os << std::setw(16) << charged_pions;
    os << std::setw(16) << protons;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"           << std::endl;
    os << std::setw(35) << "Efficiency";
    os << std::setw(16) << std::setprecision(5) << 100 * muons/double(mc_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * charged_pions/double(mc_charged_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * protons/double(mc_protons);
    os << std::endl;
    os << std::setw(35) << "Purity";
    os << std::setw(16) << std::setprecision(5) << 100 * muons/double(reco_muons);
    os << std::setw(16) << std::setprecision(5) << 100 * charged_pions/double(reco_charged_pions);
    os << std::setw(16) << std::setprecision(5) << 100 * protons/double(reco_protons);
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
  
  }
  
  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::FillGeneralParticleMisIdStatisticsFile(const EventList &ev_list, std::ofstream &os){
    /* 
     * Counters for particle type and mis-identified counterpart
     *    left particle:  reconstructed
     *    right particle: true
     *
     * e.g. muon_charged_pion = charged pion has been misidentified as a muon
     *
     */
    unsigned int muon        = 0;
    unsigned int muon_pion   = 0;
    unsigned int muon_proton = 0;
    unsigned int pion        = 0;
    unsigned int pion_muon   = 0;
    unsigned int pion_proton = 0;
    unsigned int proton      = 0;
    unsigned int proton_muon = 0;
    unsigned int proton_pion = 0;
    
    for(const Event &e : ev_list){
      if(!e.IsTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) > 1) continue;
      
      muon          += GeneralAnalysisHelper::CountMisMatchedParticles(e, 13, 13);
      muon_pion     += GeneralAnalysisHelper::CountMisMatchedParticles(e, 211, 13);
      muon_pion     += GeneralAnalysisHelper::CountMisMatchedParticles(e, -211, 13);
      muon_proton   += GeneralAnalysisHelper::CountMisMatchedParticles(e, 2212, 13);
      
      pion          += GeneralAnalysisHelper::CountMisMatchedParticles(e, 211, 211);
      pion          += GeneralAnalysisHelper::CountMisMatchedParticles(e, -211, 211);
      pion_muon     += GeneralAnalysisHelper::CountMisMatchedParticles(e, 13, 211);
      pion_muon     += GeneralAnalysisHelper::CountMisMatchedParticles(e, 13, -211);
      pion_proton   += GeneralAnalysisHelper::CountMisMatchedParticles(e, 2212, 211);
      pion_proton   += GeneralAnalysisHelper::CountMisMatchedParticles(e, 2212, -211);
      
      proton        += GeneralAnalysisHelper::CountMisMatchedParticles(e, 2212, 2212);
      proton_muon   += GeneralAnalysisHelper::CountMisMatchedParticles(e, 13, 2212);
      proton_pion   += GeneralAnalysisHelper::CountMisMatchedParticles(e, 211, 2212);
      proton_pion   += GeneralAnalysisHelper::CountMisMatchedParticles(e, -211, 2212);
    }
    
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << "    All events"                                                                           << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << std::setw(35) << "True/Reco" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << std::setw(16) << "Proton" << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    os << std::setw(35) << "Muon";
    os << std::setw(16) << muon;
    os << std::setw(16) << pion_muon;
    os << std::setw(16) << proton_muon;
    os << std::endl;
    os << std::setw(35) << "Charged pion";
    os << std::setw(16) << muon_pion;
    os << std::setw(16) << pion;
    os << std::setw(16) << proton_pion;
    os << std::endl;
    os << std::setw(35) << "Proton";
    os << std::setw(16) << muon_proton;
    os << std::setw(16) << pion_proton;
    os << std::setw(16) << proton;
    os << std::endl;
    os << "-------------------------------------------------------------------------------------"    << std::endl;
    
  }
  
  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetRecoLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths) {
    LengthWithPdg(pdg, e.GetRecoParticleList(), lengths);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetMCLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths) {
    LengthWithPdg(pdg, e.GetMCParticleList(), lengths);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::LengthWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &lengths) {
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg ) lengths.push_back(particle_list[i].GetLength());
    }
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::LongestMCTrackID(const Event &e, int &id) {
    double longest = -std::numeric_limits<double>::max();
    int longest_id = -1;
    for(const Particle & p : e.GetMCParticleList()){
      if(p.GetLength() > longest){
        longest_id = p.ID();
        longest = p.GetLength();
      }
    }
    id = longest_id;
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::LongestRecoTrackID(const Event &e, int &id) {
    double longest = -std::numeric_limits<double>::max();
    int longest_id = -1;
    for(const Particle &p : e.GetRecoParticleList()){
      if(!p.GetFromRecoTrack()) continue;
      if(p.GetLength() > longest){
        longest_id = p.ID();
        longest = p.GetLength();
      }
    }
    id = longest_id;
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::LongestMCTrackLength(const Event &e, double &length) {
    double longest = -std::numeric_limits<double>::max();
    int longest_id = -1;
    LongestMCTrackID(e,longest_id);
    if(longest_id != -1){
      for(const Particle & p : e.GetMCParticleList()){
        if(p.ID() == longest_id){
          length = p.GetLength();
          break;
        }
      }
    }
    else
      length = longest;
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::LongestRecoTrackLength(const Event &e, double &length) {
    double longest = -std::numeric_limits<double>::max();
    int longest_id = -1;
    LongestRecoTrackID(e,longest_id);
    if(longest_id != -1){
      for(const Particle &p : e.GetRecoParticleList()){
        if(!p.GetFromRecoTrack()) continue;
        if(p.ID() == longest_id){
          length = p.GetLength();
          break;
        }
      }
    }
    else
      length = longest;
  }
  
  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::GetRecoCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas) {
    CosThetaWithPdg(pdg, e.GetRecoParticleList(), cos_thetas);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetMCCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas) {
    CosThetaWithPdg(pdg, e.GetMCParticleList(), cos_thetas);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::CosThetaWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &cos_thetas) {
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) cos_thetas.push_back(particle_list[i].GetCosTheta());
    }
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetMCEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies) {
    EnergyWithPdg(pdg, e.GetMCParticleList(), energies);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetRecoEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies) {
    EnergyWithPdg(pdg, e.GetRecoParticleList(), energies);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::EnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &energies) {
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) energies.push_back(particle_list[i].GetEnergy());
    }
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetMCKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies) {
    KineticEnergyWithPdg(pdg, e.GetMCParticleList(), kinetic_energies);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies) {
    KineticEnergyWithPdg(pdg, e.GetRecoParticleList(), kinetic_energies);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &kinetic_energies) {
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) kinetic_energies.push_back(particle_list[i].GetKineticEnergy());
    }
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetMCModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod) {
    ModulusMomentumWithPdg(pdg, e.GetMCParticleList(), momentum_mod);
  }

  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::GetRecoModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod) {
    ModulusMomentumWithPdg(pdg, e.GetRecoParticleList(), momentum_mod);
  }
  
  //----------------------------------------------------------------------------------------

  void GeneralAnalysisHelper::ModulusMomentumWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &momentum_mod) {
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) momentum_mod.push_back(particle_list[i].GetModulusMomentum());
    }
  }
} // selection
