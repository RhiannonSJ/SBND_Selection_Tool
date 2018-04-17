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

    for(unsigned int i = 0; i < particle_list.size(); ++i) if(particle_list[i].GetPdgCode() == pdg) particle_counter++;

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
  
  Particle Event::GetMCParticleCharge(const Particle &particle) const{
    int charge_id = particle.GetMCParticleIdCharge();
    return GetMCParticle(charge_id, m_mc_particles);
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle Event::GetMCParticleEnergy(const Particle &particle) const{
    int energy_id = particle.GetMCParticleIdEnergy();
    return GetMCParticle(energy_id, m_mc_particles);
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle Event::GetMCParticleHits(const Particle &particle) const{
    int hits_id = particle.GetMCParticleIdHits();
    return GetMCParticle(hits_id, m_mc_particles);
  }
  
  //------------------------------------------------------------------------------------------ 
  
  Particle Event::GetMCParticle(const int id, const ParticleList &particle_list ) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetMCId() == id) return particle_list[i];
    }
    throw 8;
  }
  
 //------------------------------------------------------------------------------------------------


  float Event::GetRecoLengthWithPdg(const int pdg) const{
    return this -> LengthWithPdg( pdg, m_reco_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCLengthWithPdg(const int pdg) const{
    return this -> LengthWithPdg( pdg, m_mc_particles );
  }

  //------------------------------------------------------------------------------------------------

  float Event::LengthWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg ){ return particle_list[i].GetLength();}
    }
    return 0.;
  }

  //------------------------------------------------------------------------------------------------

  float Event::GetRecoCosThetaWithPdg(const int pdg) const{
    return this -> CosThetaWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCCosThetaWithPdg(const int pdg) const{
    return this -> CosThetaWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------

  float Event::CosThetaWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) {
	return particle_list[i].GetCosTheta();
      }
    }
    return 0.;
  }
 //------------------------------------------------------------------------------------------------


  float Event::GetMCEnergyWithPdg(const int pdg) const{
    return this -> EnergyWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoEnergyWithPdg(const int pdg) const{
    return this -> EnergyWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::EnergyWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetEnergy();
    }
    return 0.;
 }

  //------------------------------------------------------------------------------------------------


  float Event::GetMCKineticEnergyWithPdg(const int pdg) const{
    return this -> KineticEnergyWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoKineticEnergyWithPdg(const int pdg) const{
    return this -> KineticEnergyWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetKineticEnergy();
    }
    return 0.;
 }
  //------------------------------------------------------------------------------------------------


  float Event::GetMCModulusMomentumWithPdg(const int pdg) const{
    return this ->ModulusMomentumWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoModulusMomentumWithPdg(const int pdg) const{
    return this -> ModulusMomentumWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::ModulusMomentumWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetModulusMomentum();
    }
    return 0.;
 }

  //----------------------------------------------
  //---------------------------------------------------------------------------------------------

  float Event::GetDeltaEnergy( const ParticleList & particle_list ) const{
    TLorentzVector pion, p2;
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
	if( GetNuanceCode() == 1003 || GetNuanceCode() == 1009 ){
	  if( particle_list[i].GetPdgCode() ==  211 ){
	    pion.SetE( particle_list[i].GetEnergy() );
	    pion.SetVect( particle_list[i].GetMomentum() );
	  }else if( particle_list[i].GetPdgCode() == 2212 ) {
	    p2.SetE( particle_list[i].GetEnergy() );
	    p2.SetVect( particle_list[i].GetMomentum() );
	  }
	  
       	} else if( GetNuanceCode() == 1005 ){
	  if( particle_list[i].GetPdgCode() == 211 ){
	    pion.SetE( particle_list[i].GetEnergy() );
	    pion.SetVect( particle_list[i].GetMomentum() );
	  }else if( particle_list[i].GetPdgCode() == 2112 ) {
	    p2.SetE( particle_list[i].GetEnergy() );
	    p2.SetVect( particle_list[i].GetMomentum() );
	  }
	  
      }
      }
      return (pion + p2).M();
  }
    
  //-----------------------------------------------------------------------------------------------
  float Event::GetMCDeltaEnergy( ) const{
     return this -> GetDeltaEnergy( m_mc_particles );
  }
  //---------------------------------------------------------------------------------------------                                                     

  float Event::GetDeltaEnergy_p( const ParticleList & particle_list ) const{
    TLorentzVector pion, p2;                                                                                          
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      	if( particle_list[i].GetPdgCode() ==  211 ){
	  pion.SetE( particle_list[i].GetEnergy() );
	  pion.SetVect( particle_list[i].GetMomentum() );
	}else if( particle_list[i].GetPdgCode() == 2212 ) {
	  p2.SetE( particle_list[i].GetEnergy() );
	  p2.SetVect( particle_list[i].GetMomentum() );
	}
    }
    return (pion + p2).M();
  }

  
  //-------------------------------------------------------------------------------------------------

  float Event::GetCC0piNeutrinoEnergy( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
	// Get the values needed
	e    = particle_list[i].GetEnergy();
	p    = particle_list[i].GetMomentum().Mag();
	cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
	
	reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_p*m_p - m_n*m_n)*0.5);
	return reco;
      }
    }
    return reco;
  
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetCC1piNeutrinoEnergy( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_D   =   1.232;   // Delta mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    // Assuming D+ production ( xsec(n) bigger in Ar ). Modified m_n to consider the D++ production.
    m_n=(22*m_n + 18*m_p)/40.;
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
	// Get the values needed
	e    = particle_list[i].GetEnergy();
	p    = particle_list[i].GetMomentum().Mag();
	cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
	
	reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_D*m_D - m_n*m_n)*0.5);
	return reco;
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetCC1piNeutrinoEnergyMethod2( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float m_pi  = 0.13957;   // Pion mass, GeV
    float reco=0,p_mu, p_pi, e_mu, e_pi;   // track variables

    // Vector of z direction
    TVector3  z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;
    double angle_mu, angle_pi; //cos angle_mu_pi;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      for(unsigned int j = 0; j < particle_list.size(); ++j) {
	if( particle_list[i].GetPdgCode() == 13 ){
	  e_mu    = particle_list[i].GetEnergy();
	  p_mu    = particle_list[i].GetMomentum().Mag();
	  angle_mu = particle_list[i].GetCosTheta();
	  if( particle_list[j].GetPdgCode() == 211 ){
	    e_pi    = particle_list[j].GetEnergy();
	    p_pi    = particle_list[j].GetMomentum().Mag();
	    angle_pi =  particle_list[j].GetCosTheta();
	    reco = (m_mu*m_mu+m_pi*m_pi-2*m_p*(e_mu+e_pi)+2*p_mu*p_pi)/(2*(e_mu+e_pi-p_mu*angle_mu-p_pi*angle_pi-m_p));
	    return reco;
	    }
	}
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetMCCC0piNeutrinoEnergy(  ) const{
    return this -> GetCC0piNeutrinoEnergy( m_mc_particles ); 
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetMCCC1piNeutrinoEnergy(  ) const{
    return this -> GetCC1piNeutrinoEnergy( m_mc_particles );   
  }
   //-------------------------------------------------------------------------------------------------
 
  float Event::GetRecoCC0piNeutrinoEnergy(  ) const{
    return this -> GetCC0piNeutrinoEnergy( m_reco_particles ); 
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetRecoCC1piNeutrinoEnergy(  ) const{
    return this -> GetCC1piNeutrinoEnergy( m_reco_particles );   
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetRecoCC1piNeutrinoEnergyMethod2(  ) const{
    return this -> GetCC1piNeutrinoEnergyMethod2( m_reco_particles );   
  }
  //-------------------------------------------------------------------------------------------------
 
  Particle Event::GetMostEnergeticParticle(const ParticleList &particle_list) const{

    float highest_energy   = -std::numeric_limits<float>::max();
    unsigned int energy_id = std::numeric_limits<unsigned int >::max();

    for(unsigned int i = 0; i < particle_list.size(); ++i){
    
      if(!particle_list[i].GetHasCalorimetry()) continue;

      if(particle_list[i].GetEnergy() > highest_energy) energy_id = i;
    }

    return particle_list[energy_id];

  }
  
  double Event::Efficiency( const std::vector< double > & CountMC, const std::vector< double > & CountTReco, const std::vector< double > & CountReco, const TopologyMap &topology  ) const {
   ofstream rfile ;
   rfile.open( "results.txt" ) ;

   for( int i = 0; i<5; ++i ){
    rfile << "__________________________________________________________"                                                     << "\n";
    rfile                                                                                                                     << "\n";
    rfile << "                 TOPOLOGY NUMBER " << i                                                                         << "\n";
    rfile << "__________________________________________________________"                                                     << "\n";
    rfile << "CountMC = "                       << CountMC[i]                                                                 << "\n";
    rfile << "CountReco = "                     << CountReco[i]                                                               << "\n";
    rfile << "SameCount = "                     << CountTReco[i]                                                              << "\n";
    rfile << "Background = "                    << CountReco[i] - CountTReco[i]                                               << "\n";
    rfile << "Correct Reconstructed Events[%]=" << ( CountTReco[i] / CountMC[i] ) * 100                                       << "\n";
    rfile << "Purity[%]="                       << (( CountTReco[i] ) / CountReco[i] ) * 100                                  << "\n";
    rfile << "Background_Rejection[%]="         << (1-( CountReco[i] - CountTReco[i] ) / ( CountMC[0]+CountMC[1]-CountMC[i] ) ) * 100 << "\n";
    rfile << "__________________________________________________________"                                                     << "\n";
 }

   if( topology == GeneralAnalysisHelper::GetNCTopologyMap() ) return ( CountTReco[0] / CountMC[0] ) * 100 ;
   if( topology == GeneralAnalysisHelper::GetCCIncTopologyMap() ) return ( CountTReco[1] / CountMC[1] ) * 100 ;
   if( topology == GeneralAnalysisHelper::GetCC0PiTopologyMap() ) return ( CountTReco[2] / CountMC[2] ) * 100 ;
   if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) return ( CountTReco[3] / CountMC[3] ) * 100 ;
   if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() ) return ( CountTReco[4] / CountMC[4] ) * 100 ;
   return 0;
  }

  void Event::SaveTopologyMatrix( const ParticleMatrix & Count_MC_Topology, const ParticleMatrix & Count_TReco_Topology, const ParticleMatrix & Count_Reco_Topology ) const {
   ofstream TMfile ;
   TMfile.open( "TopologyMatrix.txt" ) ;
    TMfile                                                                                                                       << "\n";
    TMfile << "____________________________________________________________"                                                     << "\n";
    TMfile                                                                                                                       << "\n";
    TMfile << "  TOPOLOGY MATRIX - TRUE RECO  (#_TReco / #_Total_MC) : "                                                         << "\n";
    TMfile << "____________________________________________________________"                                                     << "\n";
    for( unsigned int i = 0 ; i < 5; ++i ){
      TMfile << "(";
      for( unsigned int k = 0 ; k < 5 ; ++k ) {
	if( Count_TReco_Topology[i][k]!=0 ){
	  TMfile << ( Count_TReco_Topology[i][k] / Count_MC_Topology[i][k] ) * 100 << "   ,   ";
	} else TMfile << "   --   ";
    }
      TMfile <<")"                                                                                                               << "\n";
    }
    
  }

  void Event::EventInformationParticles( std::string event_file, const int event_number ) const {
    ofstream efile ;
    efile.open( event_file , std::ofstream::app);

    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << " EVENT NUMBER =                                            " << event_number                 << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "TRUE EVENTS      : "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "   muons         : " << CountMCParticlesWithPdg(13)                                          << "\n";
    efile << "   pi+/-         : " << CountMCParticlesWithPdg(211) + CountMCParticlesWithPdg(-211)         << "\n";
    efile << "   pi0           : " << CountMCParticlesWithPdg(111)                                         << "\n";
    efile << "   protons       : " << CountMCParticlesWithPdg(2212)                                        << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "TRUE TOPOLOGY    : "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    if(CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()))     efile << "   NC            : " << "TRUE"                  << "\n";
    if(CheckMCTopology(GeneralAnalysisHelper::GetCCIncTopologyMap())) efile << "   ccincl.       : " << "TRUE"                  << "\n";
    if(CheckMCTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap())) efile << "   cc0pi         : " << "TRUE"                  << "\n";
    if(CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap())) efile << "   cc1pi+/-      : " << "TRUE"                  << "\n";
    if(CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap())) efile << "   cc1pi0        : " << "TRUE"                  << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << " SELECTED EVENTS :                                         "                                 << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "   muons         : " << CountRecoParticlesWithPdg(13)                                        << "\n";
    efile << "   pi+/-         : " << CountRecoParticlesWithPdg(211) + CountRecoParticlesWithPdg(-211)     << "\n";
    efile << "   pi0           : " << CountRecoParticlesWithPdg(111)                                       << "\n";
    efile << "   protons       : " << CountRecoParticlesWithPdg(2212)                                      << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "SELECTED TOPOLOGY: "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    if(CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()))     efile << "   NC            : " << "TRUE"                << "\n";
    if(CheckRecoTopology(GeneralAnalysisHelper::GetCCIncTopologyMap())) efile << "   ccincl.       : " << "TRUE"                << "\n";
    if(CheckRecoTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap())) efile << "   cc0pi         : " << "TRUE"                << "\n";
    if(CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap())) efile << "   cc1pi+/-      : " << "TRUE"                << "\n";
    if(CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap())) efile << "   cc1pi0        : " << "TRUE"                << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";


  }

  void Event::EventProperties( const TopologyMap &topology, std::string event_file, const int event_number ) const {
    ofstream lfile, afile, Kfile ;
    lfile.open( event_file+="_length.txt" , std::ofstream::app);
    afile.open( event_file+="_angle.txt"  , std::ofstream::app);
    Kfile.open( event_file+="_KineticEnergy.txt" , std::ofstream::app);

    if( CheckMCTopology( topology ) ) { /** Change the topology here **/
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "EVENT NUMBER      : " << event_number                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "LENGTH INFORMATION: "                                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "TRUE EVENTS       : "                                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "Muon length         : " << GetMCLengthWithPdg(  13  )         << "\n";
       lfile << "pi+/- length        : " << GetMCLengthWithPdg( 211  )  << "\n";
       lfile << "pi0   length        : " << GetMCLengthWithPdg( 111  )  << "\n";
       lfile << "p length            : " << GetMCLengthWithPdg( 2212 ) << "\n";
       
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "EVENT NUMBER       : " << event_number                                                       << "\n";
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "ANGULAR INFORMATION: "                                                                       << "\n";
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "TRUE EVENTS        : "                                                                       << "\n"; 
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "Muon angle       : " << GetMCCosThetaWithPdg(  13  )  << "\n";
       afile << "Pion angle       : " << GetMCCosThetaWithPdg( 211  )  << "\n";
       afile << "pi0   angle      : " << GetMCCosThetaWithPdg( 111  )  << "\n";
       afile << "Proton angle     : " << GetMCCosThetaWithPdg( 2212 )  << "\n";

       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "EVENT NUMBER       : " << event_number                                                       << "\n";
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "Muon Kenergy       : " << GetMCKineticEnergyWithPdg(  13  )  << "\n";
       Kfile << "Pion Kenergy       : " << GetMCKineticEnergyWithPdg( 211  )  << "\n";
       
       if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy     : " << GetMCKineticEnergyWithPdg( 2212 )  << "\n";

       if( CheckRecoTopology( topology ) ) {

	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 lfile << "SIGNAL EVENTS     : "                                                                        << "\n";
	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()) )lfile << "Muon length          : " << GetRecoLengthWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )lfile << "pi+/- length      : "  << GetRecoLengthWithPdg( 211  ) << "\n";
	 if( CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )lfile << "pi0   length      : "  << GetRecoLengthWithPdg( 111  ) << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 )lfile << "p length            : " << GetRecoLengthWithPdg( 2212 )  << "\n";
	 

	 afile << "-----------------------------------------------------------"                                 << "\n";
	 afile << "SIGNAL EVENTS      : "                                                                       << "\n"; 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap())) afile << "Muon angle              : " << GetRecoCosThetaWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )afile << "Pion angle         : " << GetRecoCosThetaWithPdg( 211  )  << "\n";
	 if( CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )afile << "pi0   angle        : " << GetRecoCosThetaWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) afile << "Proton angle        : " << GetRecoCosThetaWithPdg( 2212 )  << "\n";
 	 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) Kfile << "Muon Kenergy         : " << GetRecoKineticEnergyWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )Kfile << "Pion Kenergy       : " << GetRecoKineticEnergyWithPdg( 211  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy      : " << GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
       }

     }
    if( CheckRecoTopology( topology ) ) { // Change topology here

	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 lfile << "SELECTED EVENTS   : "                                                                        << "\n";
	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) lfile << "Muon length       : " << GetRecoLengthWithPdg(  13  )         << "\n";
	 if( CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )lfile << "pi+/- length      : " << GetRecoLengthWithPdg( 211  )  << "\n";
	 if( CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )lfile << "pi0   length      : " << GetRecoLengthWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) lfile << "p length          : " << GetRecoLengthWithPdg( 2212 ) << "\n";
	 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 afile << "SELECTED EVENTS      : "                                                                     << "\n"; 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) afile << "Muon angle         : " << GetRecoCosThetaWithPdg(  13  )  << "\n";
	 if( CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )afile << "Pion angle         : " << GetRecoCosThetaWithPdg( 211  )  << "\n";
	 if( CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )afile << "pi0   angle        : " << GetRecoCosThetaWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) afile << "Proton angle       : " << GetRecoCosThetaWithPdg( 2212 )  << "\n";

	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) Kfile << "Muon Kenergy       : " << GetRecoKineticEnergyWithPdg(  13  )  << "\n";
	 if( CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )Kfile << "Pion Kenergy       : " << GetRecoKineticEnergyWithPdg( 211  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy     : " << GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
     }

  }

} // selection
