#include "../include/Event.h"
#include "../include/EventSelectionTool.h"
namespace selection{
  
  Event::Event(const ParticleList &mc_particles, 
               const ParticleList &reco_particles, 
               const unsigned int interaction, 
               const unsigned int scatter, 
               const int neutrino_pdg, 
               const unsigned int charged_pi, 
               const unsigned int neutral_pi, 
               const bool is_cc, 
               const TVector3 &mc_vertex, 
               const TVector3 &reco_vertex, 
               const float neutrino_energy, 
               const float neutrino_qsqr, 
               const int &file, 
               const int &id, 
               const float &baseline,
               const Geometry &g) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_interaction(interaction),
    m_scatter(scatter),
    m_nu_pdg(neutrino_pdg),
    m_charged_pi(charged_pi),
    m_neutral_pi(neutral_pi),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex), 
    m_neutrino_energy(neutrino_energy),
    m_neutrino_qsqr(neutrino_qsqr),
    m_file(file),
    m_id(id),
    m_baseline(baseline),
    m_geometry(g){}

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
  
  int Event::GetFileId() const{
    return m_file;
  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetId() const{
    return m_id;
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
  
  int Event::GetInteractionType() const{
  
    return m_interaction;
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int Event::GetPhysicalProcess() const{
  
    return m_scatter;
  
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
  
  TVector3 Event::GetMinimumFiducialDimensions() const{
    return TVector3(*std::min_element(m_geometry.GetMinX().begin(),m_geometry.GetMinX().end()),
                    *std::min_element(m_geometry.GetMinY().begin(),m_geometry.GetMinY().end()),
                    *std::min_element(m_geometry.GetMinZ().begin(),m_geometry.GetMinZ().end()));
  }

  //------------------------------------------------------------------------------------------ 
  
  /*
  std::pair<float, float> Event::GetCentralFiducialDimensions() const{
    std::pair<float, float> central(-m_sbnd_border_x_max1, m_sbnd_border_x_min2);
    return central;
  }*/

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Event::GetMaximumFiducialDimensions() const{
    return TVector3(*std::max_element(m_geometry.GetMaxX().begin(),m_geometry.GetMaxX().end()),
                    *std::max_element(m_geometry.GetMaxY().begin(),m_geometry.GetMaxY().end()),
                    *std::max_element(m_geometry.GetMaxZ().begin(),m_geometry.GetMaxZ().end()));
  }
  
  //------------------------------------------------------------------------------------------ 
  
  bool Event::IsTrueFiducial() const{
       
    // Check the neutrino interaction vertex is within the fiducial volume      
    float vertex_x = m_mc_vertex[0];                        
    float vertex_y = m_mc_vertex[1];                        
    float vertex_z = m_mc_vertex[2];

    bool in_fiducial = false;
    for(unsigned int n = 0; n < m_geometry.GetNTPCs(); ++n){
      float min_x = m_geometry.GetMinX().at(n);
      float min_y = m_geometry.GetMinY().at(n);
      float min_z = m_geometry.GetMinZ().at(n);
      float max_x = m_geometry.GetMaxX().at(n);
      float max_y = m_geometry.GetMaxY().at(n);
      float max_z = m_geometry.GetMaxZ().at(n);
      
      if(vertex_x >= min_x && vertex_x <= max_x &&
         vertex_y >= min_y && vertex_y <= max_y &&
         vertex_z >= min_z && vertex_z <= max_z) in_fiducial = true;
    }
    return in_fiducial;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Event::IsRecoFiducial() const{
       
    // Check the neutrino interaction vertex is within the fiducial volume      
    float vertex_x = m_reco_vertex[0];                        
    float vertex_y = m_reco_vertex[1];                        
    float vertex_z = m_reco_vertex[2];                
   
    bool in_fiducial = false;
    for(unsigned int n = 0; n < m_geometry.GetNTPCs(); ++n){
      float min_x = m_geometry.GetMinX().at(n);
      float min_y = m_geometry.GetMinY().at(n);
      float min_z = m_geometry.GetMinZ().at(n);
      float max_x = m_geometry.GetMaxX().at(n);
      float max_y = m_geometry.GetMaxY().at(n);
      float max_z = m_geometry.GetMaxZ().at(n);

      if(vertex_x >= min_x && vertex_x <= max_x &&
         vertex_y >= min_y && vertex_y <= max_y &&
         vertex_z >= min_z && vertex_z <= max_z) in_fiducial = true;
    }
    return in_fiducial;
  }

  //------------------------------------------------------------------------------------------ 
 
  bool Event::AllRecoContained() const{
    for(const Particle &p : m_reco_particles){
      if(p.GetFromRecoTrack() && !p.GetTrackContained()) return false;
    }
    return true;
  }
  
  //------------------------------------------------------------------------------------------ 

  int Event::NumberOfEscapingRecoParticles() const{
    return this->Event::NumberOfEscapingParticles(this->Event::GetRecoParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
      
  int Event::NumberOfEscapingMCParticles() const{
    return this->Event::NumberOfEscapingParticles(this->Event::GetMCParticleList());
  }
  
  //------------------------------------------------------------------------------------------ 
      
  int Event::NumberOfEscapingParticles(const ParticleList &particles) const{
    int escaping = 0;
    for(const Particle &p : particles){
      if(!p.GetFromRecoTrack() || !p.GetHasCalorimetry()) continue;
      if(p.GetOneEndTrackContained()) escaping++;
    }
    return escaping;
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
  
  float Event::GetTrueNuQ2() const{
  
    return m_neutrino_qsqr;

  }
  
  //------------------------------------------------------------------------------------------ 
  
  float Event::GetBaseline() const{
    
    return m_baseline;

  }

  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const{
    unsigned int particle_counter = 0;
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg && particle_list[i].GetNumberOfHits() >= 5) {
        if(pdg == 2212) {// Make sure the proton energy is realistic for reconstruction
          if(particle_list[i].GetKineticEnergy() > 0.021) particle_counter++;
        }
        else particle_counter++;
      }
    }
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
