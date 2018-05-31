#include "../include/Event.h"
#include "../include/EventSelectionTool.h"
namespace selection{
  
  Event::Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const int neutrino_pdg, const unsigned int charged_pi, const unsigned int neutral_pi, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex, const float neutrino_energy) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_nuance(nuance),
    m_nu_pdg(neutrino_pdg),
    m_charged_pi(charged_pi),
    m_neutral_pi(neutral_pi),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex), 
    m_neutrino_energy(neutrino_energy) {
    
      // Co-ordinate offset in cm
      m_sbnd_half_length_x = 400;
      m_sbnd_half_length_y = 400;
      m_sbnd_half_length_z = 500;
      
      m_sbnd_offset_x = 200;
      m_sbnd_offset_y = 200;
      m_sbnd_offset_z = 0;

      m_sbnd_border_x = 10;
      m_sbnd_border_y = 20;
      m_sbnd_border_z = 10;

    }


  //------------------------------------------------------------------------------------------ 
    
  unsigned int Event::CountMCParticlesWithPdg(const int pdg) const{
 
    return this->CountParticlesWithPdg(pdg, m_mc_particles);

  }

  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountRecoParticlesWithPdg(const int pdg) const{
    
    return this->CountParticlesWithPdg(pdg, m_reco_particles);

  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckMCTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckRecoTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_reco_particles);

  }
  
  //------------------------------------------------------------------------------------------ 

  Particle Event::GetMostEnergeticRecoParticle() const{
 
    return this->GetMostEnergeticParticle(m_reco_particles);

  }
  
  //------------------------------------------------------------------------------------------ 

  Particle Event::GetMostEnergeticTrueParticle() const{
 
    return this->GetMostEnergeticParticle(m_mc_particles);

  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetMCParticleList() const{
  
    return m_mc_particles;

  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetRecoParticleList() const{
  
    return m_reco_particles;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNuanceCode() const{
  
    return m_nuance;
  
  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNeutrinoPdgCode() const{
  
    return m_nu_pdg;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNChargedPions() const{
  
    return m_charged_pi;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNNeutralPions() const{
  
    return m_neutral_pi;

  }

  //------------------------------------------------------------------------------------------ 
  
  bool Event::IsSBNDTrueFiducial() const{
       
    // Check the neutrino interaction vertex is within the fiducial volume      
     float nu_vertex_x = m_mc_vertex[0];                        
     float nu_vertex_y = m_mc_vertex[1];                        
     float nu_vertex_z = m_mc_vertex[2];                        
                                                                                 
     if (    (nu_vertex_x > (m_sbnd_half_length_x - m_sbnd_offset_x - m_sbnd_border_x)) 
          || (nu_vertex_x < (-m_sbnd_offset_x + m_sbnd_border_x))          
          || (nu_vertex_y > (m_sbnd_half_length_y - m_sbnd_offset_y - m_sbnd_border_y)) 
          || (nu_vertex_y < (-m_sbnd_offset_y + m_sbnd_border_y))          
          || (nu_vertex_z > (m_sbnd_half_length_z - m_sbnd_offset_z - m_sbnd_border_z)) 
          || (nu_vertex_z < (-m_sbnd_offset_z + m_sbnd_border_z))) return false; 

     return true;
  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetPhysicalProcess() const{

    // QEL
    if(m_nuance == 0 
    || m_nuance == 1001 
    || m_nuance == 1002) return 0;
    // MEC
    else if(m_nuance == 10) return 1;
    // RES
    else if(m_nuance == 1 
         || m_nuance == 1003 
         || m_nuance == 1004
         || m_nuance == 1005
         || m_nuance == 1006
         || m_nuance == 1007
         || m_nuance == 1008
         || m_nuance == 1009
         || m_nuance == 1010) return 2;
    // DIS
    else if(m_nuance == 2
         || m_nuance == 1091) return 3;
    // COH
    else if(m_nuance == 1097) return 4;
    // Non RES 1pi
    else if(m_charged_pi + m_neutral_pi == 1) return 5;
    // Other
    else return 6;
  
  }
      
  //------------------------------------------------------------------------------------------ 

  bool Event::GetIsCC() const{
  
    return m_is_cc;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetMCNuVertex() const{
  
    return m_mc_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetRecoNuVertex() const{
  
    return m_reco_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 
  
  float Event::GetTrueNuEnergy() const{
  
    return m_neutrino_energy;

  }
  
  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const{
    unsigned int particle_counter = 0;
    for(unsigned int i = 0; i < particle_list.size(); ++i) if(particle_list[i].GetPdgCode() == pdg && particle_list[i].GetNumberOfHits() > 5) particle_counter++;
    return particle_counter;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Event::CheckTopology(const TopologyMap &topology, const ParticleList &particle_list) const{
    // Loop over the map
    for( TopologyMap::const_iterator it = topology.begin(); it != topology.end(); ++it ){
      // Define temporary variables for the current map element
      std::vector< int > pdg_codes = it->first; 
      int n_total                  = it->second;
      // Count the number of particles in the current event with the same PDG codes 
      // as given by the chosen topology
      int counter = 0;
      // Loop over particles in current event
      for(unsigned int i = 0; i < pdg_codes.size(); ++i){
        counter += this->CountParticlesWithPdg(pdg_codes[i], particle_list);
      }
      if(counter != n_total) return false;
    }
    return true;
  }
  
  //------------------------------------------------------------------------------------------
 
  Particle Event::GetMostEnergeticParticle(const ParticleList &particle_list) const{
    float highest_energy   = -std::numeric_limits<float>::max();
    unsigned int energy_id = std::numeric_limits<unsigned int >::max();
    for(unsigned int i = 0; i < particle_list.size(); ++i){
      if(!particle_list[i].GetHasCalorimetry()) continue;
      if(particle_list[i].GetEnergy() > highest_energy) energy_id = i;
    }
    return particle_list[energy_id];
  }
  
} // selection
