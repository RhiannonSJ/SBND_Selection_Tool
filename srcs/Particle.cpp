#include <string>
#include <vector>
#include <cmath>
#include "TVector3.h"
#include "../include/Particle.h"

namespace selection{

  Particle::Particle(const int mc_id,
                     const int pdg, 
                     const int status, 
                     const int n_hits,  
                     const float mass, 
                     const float energy,
                     const TVector3 &vertex,
                     const TVector3 &end,
                     const TVector3 &momentum,
                     const Geometry &g) : 
    m_mc_id(mc_id),
    m_pdg(pdg),
    m_status(status),
    m_n_hits(n_hits),
    m_mass(mass),
    m_energy(energy),
    m_has_calorimetry(true),
    m_from_reco_track(false),
    m_vertex(vertex),
    m_end(end),
    m_momentum(momentum),
    m_geometry(g){
      m_length = sqrt(pow(end[0] - vertex[0], 2) + pow(end[1] - vertex[1], 2) + pow(end[0] - vertex[0], 2));
      m_id = m_mc_id;
    }

  //------------------------------------------------------------------------------------------ 
  
  Particle::Particle(const int mc_id_charge, 
                     const int mc_id_energy, 
                     const int mc_id_hits,
                     const int id,
                     const int pdg, 
                     const int n_hits, 
                     const float kinetic_energy, 
                     const float mcs_momentum_muon, 
                     const float range_momentum_muon, 
                     const float range_momentum_proton, 
                     const float length, 
                     const TVector3 &vertex, 
                     const TVector3 &end, 
                     const float &chi2p, 
                     const float &chi2mu, 
                     const float &chi2pi, 
                     const std::vector<float> &dedx, 
                     const std::vector<float> &residual_range, 
                     const Geometry &g) :
    m_mc_id_charge(mc_id_charge),
    m_mc_id_energy(mc_id_energy),
    m_mc_id_hits(mc_id_hits),
    m_id(id),
    m_pdg(pdg),
    m_n_hits(n_hits),
    m_kinetic_energy(kinetic_energy),
    m_mcs_mom_muon(mcs_momentum_muon),
    m_range_mom_muon(range_momentum_muon),
    m_range_mom_proton(range_momentum_proton),
    m_length(length),
    m_has_calorimetry(true),
    m_from_reco_track(true),
    m_vertex(vertex),
    m_end(end),
    m_chi2p(chi2p),
    m_chi2mu(chi2mu),
    m_chi2pi(chi2pi),
    m_dedx(dedx),
    m_residual_range(residual_range),
    m_geometry(g){
      // Set member variables
      m_mass   = this->GetMassFromPdg(pdg);
      m_energy = m_mass + (kinetic_energy/1000.);
      
      // Get the magnitude of the momentum
      double momentum_magnitude = sqrt(pow(m_energy,2) - pow(m_mass,2));
      m_momentum = momentum_magnitude * (m_end - m_vertex)*(1/double((m_end-m_vertex).Mag()));
    }

  //------------------------------------------------------------------------------------------ 

  Particle::Particle(const int pdg,
                     const int id,
                     const int n_hits, 
                     const TVector3 &vertex, 
                     const TVector3 &end, 
                     const float &energy,
                     const Geometry &g) :
    m_pdg(pdg),
    m_id(id),
    m_n_hits(n_hits),
    m_mass(this->GetMassFromPdg(pdg)),
    m_has_calorimetry(true),
    m_from_reco_track(false),
    m_vertex(vertex),
    m_end(end),
    m_energy(energy),
    m_geometry(g){
      m_length = sqrt(pow(end[0] - vertex[0], 2) + pow(end[1] - vertex[1], 2) + pow(end[0] - vertex[0], 2));
    }

  //------------------------------------------------------------------------------------------ 
  
  int Particle::ID() const{
    return m_id;
  }
  float Particle::GetMassFromPdg(const int pdg) const{
    switch(abs(pdg)){ 
      case 211: 
        return 0.1395701; 
      case 111:
        return 0.1349766;
      case 13:
        return 0.1056583;
      case 11:
        return 0.0005489;
      case 22:
        return 0.0;
      case 2212:
        return 0.9382720;
      case 321:
        return 0.4936770;
      default:
        throw 2;
    }
  }

  //------------------------------------------------------------------------------------------ 
  
  int Particle::GetPdgCode() const{return m_pdg;}

  //------------------------------------------------------------------------------------------ 
  
  int Particle::GetStatusCode() const{
    // Make sure we are looking at an mcparticle
    if(m_has_calorimetry && !m_from_reco_track){
      return m_status;
    }
    std::cerr << "GetStatusCode" << std::endl;
    throw 5;
  }
      
  //------------------------------------------------------------------------------------------ 
  
  int Particle::GetNumberOfHits() const{return m_n_hits;}

  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetMass() const{return m_mass;}

  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetEnergy() const{
    if(!m_has_calorimetry) {
      std::cerr << " No calorimetry " << std::endl;
      throw 1;
    }
    return m_energy;
  }

  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetKineticEnergy() const{
    if(!m_has_calorimetry) {
      std::cerr << " No calorimetry " << std::endl;
      throw 1;
    }
    return m_energy-m_mass;
  }

  //------------------------------------------------------------------------------------------ 
 
  float Particle::GetChi2P() const{
    if(!m_from_reco_track) {
      std::cerr << " Not a track " << std::endl;
      throw 1;
    }
    return m_chi2p;
  }

  //------------------------------------------------------------------------------------------ 
 
  float Particle::GetChi2Mu() const{
    if(!m_from_reco_track) {
      std::cerr << " Not a track " << std::endl;
      throw 1;
    }
    return m_chi2mu;
  }

  //------------------------------------------------------------------------------------------ 
 
  float Particle::GetChi2Pi() const{
    if(!m_from_reco_track) {
      std::cerr << " Not a track " << std::endl;
      throw 1;
    }
    return m_chi2pi;
  }

  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetLength() const{return m_length;}

  //------------------------------------------------------------------------------------------ 
 
  std::vector<float> Particle::GetdEdx() const{
    if(!m_from_reco_track) {
      std::cerr << " Not a track " << std::endl;
      throw 1;
    }
    return m_dedx;
  }

  //------------------------------------------------------------------------------------------ 
 
  std::vector<float> Particle::GetResidualRange() const{
    if(!m_from_reco_track) {
      std::cerr << " Not a track " << std::endl;
      throw 1;
    }
    return m_residual_range;
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetVertex() const{return m_vertex;}

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetEnd() const{return m_end;}

  //------------------------------------------------------------------------------------------ 
  
  void Particle::FlipTrack(){
    TVector3 temp;
    temp     = m_vertex;
    m_vertex = m_end;
    m_end    = temp;
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetMomentum() const{
    if(!m_has_calorimetry) {
      std::cerr << " No calorimetry " << std::endl;
      throw 1;
    }
    if(!m_from_reco_track) return m_momentum;
    else if(this->Particle::GetPdgCode() == 13 && this->Particle::GetOneEndTrackContained() && this->Particle::GetLength() >= 100)
      return this->Particle::GetMCSMomentumMuon();
    else if(this->Particle::GetPdgCode() == 13 && !this->Particle::GetTrackContained() && this->Particle::GetLength() >= 50)
      return this->Particle::GetRangeMomentumMuon();
    else if(this->Particle::GetPdgCode() == 2212 && !this->Particle::GetTrackContained() && this->Particle::GetLength() >= 50)
      return this->Particle::GetRangeMomentumProton();
    else
      return m_momentum;
  }

  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetModulusMomentum() const{
    if(!m_has_calorimetry) {
      std::cerr << " No calorimetry " << std::endl;
      throw 1;
    }
    if(!m_from_reco_track) return m_momentum.Mag();
    else if(this->Particle::GetPdgCode() == 13 && this->Particle::GetOneEndTrackContained() && this->Particle::GetLength() >= 100)
      return m_mcs_mom_muon;
    else if(this->Particle::GetPdgCode() == 13 && this->Particle::GetTrackContained() && this->Particle::GetLength() >= 50)
      return m_range_mom_muon;
    else if(this->Particle::GetPdgCode() == 2212 && this->Particle::GetTrackContained() && this->Particle::GetLength() >= 50)
      return m_range_mom_proton;
    else
      return m_momentum.Mag();
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetMCSMomentumMuon() const{
    return m_mcs_mom_muon * (1. / (this->Particle::GetEnd() - this->Particle::GetVertex()).Mag()) * (this->Particle::GetEnd() - this->Particle::GetVertex());
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetRangeMomentumMuon() const{
    return m_range_mom_muon * (1. / (this->Particle::GetEnd() - this->Particle::GetVertex()).Mag()) * (this->Particle::GetEnd() - this->Particle::GetVertex());
  }

  //------------------------------------------------------------------------------------------ 
  
  TVector3 Particle::GetRangeMomentumProton() const{
    return m_range_mom_proton * (1. / (this->Particle::GetEnd() - this->Particle::GetVertex()).Mag()) * (this->Particle::GetEnd() - this->Particle::GetVertex());
  }
  
  //------------------------------------------------------------------------------------------ 

  int Particle::GetMCId() const{ 
    if(m_has_calorimetry && !m_from_reco_track){
      return m_mc_id;
    } 
    std::cout << "GetMCId" << std::endl;
    throw 5;
  }

  //------------------------------------------------------------------------------------------ 

  int Particle::GetMCParticleIdCharge() const{
    // Check that the particle has been made using reconstructed tracks
    if(m_has_calorimetry && m_from_reco_track){
      return m_mc_id_charge;
    }
    
    std::cout << "GetMCIdCharge" << std::endl;

    throw 6;
  }

  //------------------------------------------------------------------------------------------ 

  int Particle::GetMCParticleIdEnergy() const{
    // Check that the particle has been made using reconstructed tracks
    if(m_has_calorimetry && m_from_reco_track){
      return m_mc_id_energy;
    }
    
    std::cout << "GetMCIdEnergy" << std::endl;

    throw 6;
  }

  //------------------------------------------------------------------------------------------ 

  int Particle::GetMCParticleIdHits() const{
    // Check that the particle has been made using reconstructed tracks
    if(m_has_calorimetry && m_from_reco_track){
      return m_mc_id_hits;
    }
    
    std::cout << "GetMCIdHits" << std::endl;

    throw 6;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Particle::GetHasCalorimetry() const{return m_has_calorimetry;}

  //------------------------------------------------------------------------------------------ 
  
  bool Particle::GetFromRecoTrack() const{return m_from_reco_track;}
  
  //------------------------------------------------------------------------------------------ 
  
  float Particle::GetCosTheta() const{
    TVector3 u_z(0,0,1);
    return (m_end - m_vertex).Dot( u_z )/((m_end - m_vertex).Mag());
  }

  //------------------------------------------------------------------------------------------ 
 
  bool Particle::GetTrackContained() const{

    if(!m_from_reco_track && !m_has_calorimetry) 
      std::cerr << "Containment can only be checked for reconstructed tracks or MCParticles" << std::endl;

    // Check the neutrino interaction vertex is within the fiducial volume
    float vertex_x = m_vertex[0];                        
    float vertex_y = m_vertex[1];                        
    float vertex_z = m_vertex[2];                        
    float end_x    = m_end[0];                        
    float end_y    = m_end[1];                        
    float end_z    = m_end[2];                        
   
    // Get the outermost faces of the detector
    float min_x = *std::min_element(m_geometry.GetMinX().begin(),m_geometry.GetMinX().end());
    float min_y = *std::min_element(m_geometry.GetMinY().begin(),m_geometry.GetMinY().end());
    float min_z = *std::min_element(m_geometry.GetMinZ().begin(),m_geometry.GetMinZ().end());
    float max_x = *std::max_element(m_geometry.GetMaxX().begin(),m_geometry.GetMaxX().end());
    float max_y = *std::max_element(m_geometry.GetMaxY().begin(),m_geometry.GetMaxY().end());
    float max_z = *std::max_element(m_geometry.GetMaxZ().begin(),m_geometry.GetMaxZ().end());

    if(vertex_x < min_x || vertex_x > max_x ||
       vertex_y < min_y || vertex_y > max_y ||
       vertex_z < min_z || vertex_z > max_z ||
       end_x    < min_x || end_x    > max_x ||
       end_y    < min_y || end_y    > max_y ||
       end_z    < min_z || end_z    > max_z) return false;

    return true;
  }
  
  //------------------------------------------------------------------------------------------ 
 
  bool Particle::GetOneEndTrackContained() const{

    if(!m_from_reco_track && !m_has_calorimetry) 
      std::cerr << "Containment can only be checked for reconstructed tracks or MCParticles" << std::endl;

    // Check the neutrino interaction vertex is within the fiducial volume
    float vertex_x = m_vertex[0];                        
    float vertex_y = m_vertex[1];                        
    float vertex_z = m_vertex[2];                        
    float end_x    = m_end[0];                        
    float end_y    = m_end[1];                        
    float end_z    = m_end[2];                        
   
    // Get the outermost faces of the detector
    float min_x = *std::min_element(m_geometry.GetMinX().begin(),m_geometry.GetMinX().end());
    float min_y = *std::min_element(m_geometry.GetMinY().begin(),m_geometry.GetMinY().end());
    float min_z = *std::min_element(m_geometry.GetMinZ().begin(),m_geometry.GetMinZ().end());
    float max_x = *std::max_element(m_geometry.GetMaxX().begin(),m_geometry.GetMaxX().end());
    float max_y = *std::max_element(m_geometry.GetMaxY().begin(),m_geometry.GetMaxY().end());
    float max_z = *std::max_element(m_geometry.GetMaxZ().begin(),m_geometry.GetMaxZ().end());

    bool does_vtx_escape = 
      (vertex_x < min_x || vertex_x > max_x ||
       vertex_y < min_y || vertex_y > max_y ||
       vertex_z < min_z || vertex_z > max_z);
    bool does_end_escape = 
      (end_x    < min_x || end_x    > max_x ||
       end_y    < min_y || end_y    > max_y ||
       end_z    < min_z || end_z    > max_z);

    if(does_vtx_escape && does_end_escape) return false;
    if(!does_vtx_escape && !does_end_escape) return false;

    return true;
  }
} // Selection
