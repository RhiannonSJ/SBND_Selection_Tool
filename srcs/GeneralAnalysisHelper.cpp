#include "../include/GeneralAnalysisHelper.h"

namespace selection{

  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNCTopologyMap() {
    TopologyMap signal_map_nc;
    signal_map_nc.insert(TopologyMap::value_type({13},0));
    return signal_map_nc;
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
  TopologyMap GeneralAnalysisHelper::GetCC1PiTopologyMap() { 
    TopologyMap signal_map_cc_1pi;
    signal_map_cc_1pi.insert(TopologyMap::value_type({13},1));
    signal_map_cc_1pi.insert(TopologyMap::value_type({211, -211},1));
    return signal_map_cc_1pi;
  }
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCCPi0TopologyMap() {
    TopologyMap signal_map_cc_pi0;
    signal_map_cc_pi0.insert(TopologyMap::value_type({13},1));
    signal_map_cc_pi0.insert(TopologyMap::value_type({111},1));
    return signal_map_cc_pi0;
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
  
  int GeneralAnalysisHelper::ParticleHasAMatch(const Event &e, const Particle &p){
    // Starting from hits (since this is the chosen best method) find out if there is a match
    ParticleList particles = e.GetMCParticleList();
    for(const Particle &part : particles){
      if(part.GetMCId() == p.GetMCParticleIdHits()) return 0;
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
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetMCId() == id) return particle_list[i];
    }
    throw 8;
  }
  
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesByTopologyTrue(const Event &e, const TopologyMap &topology, const int pdg){
    // Check if the event is a true event
    if(e.IsSBNDTrueFiducial() && e.CheckMCTopology(topology)){
      return GeneralAnalysisHelper::CountMatchedParticles(e, e.GetRecoParticleList(), pdg);
    }
    else return 0;  
  }
      
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg){
    // Check if the event is a selected event
    if(e.CheckRecoTopology(topology)){
      return GeneralAnalysisHelper::CountMatchedParticles(e, e.GetRecoParticleList(), pdg);
    }
    else return 0;  
  }
      
  //----------------------------------------------------------------------------------------

  unsigned int GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg){
    // Check if the event is a signal event
    if(e.IsSBNDTrueFiducial() && e.CheckRecoTopology(topology) && e.CheckMCTopology(topology)){
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
    for(unsigned int i = 0; i < particle_list.size(); ++i){
      if(particle_list[i].GetPdgCode() == abs(pdg)){
        if(GeneralAnalysisHelper::MatchedParticle(e,particle_list[i])) matched_particles++;
      }
    }
    return matched_particles;
  }
  
  //----------------------------------------------------------------------------------------
  
  bool GeneralAnalysisHelper::MatchedParticle(const Event &e, const Particle &p){
    if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 0      && abs(GeneralAnalysisHelper::GetMCParticleHits(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
    else if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 1 && abs(GeneralAnalysisHelper::GetMCParticleCharge(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
    else if(GeneralAnalysisHelper::ParticleHasAMatch(e, p) == 2 && abs(GeneralAnalysisHelper::GetMCParticleEnergy(e, p).GetPdgCode()) == p.GetPdgCode()) return true;
    else return false;
  }
  
  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(const EventList &ev_list, const TopologyMap &topology, const std::string &topology_name, std::ofstream &os){

    // Counters for each particle type
    unsigned int total_true_muons         = 0;
    unsigned int total_true_protons       = 0;
    unsigned int total_true_charged_pions = 0;
    unsigned int total_reco_muons         = 0;
    unsigned int total_reco_protons       = 0;
    unsigned int total_reco_charged_pions = 0;
    unsigned int true_muons               = 0;
    unsigned int true_protons             = 0;
    unsigned int true_charged_pions       = 0;
    unsigned int selected_muons           = 0;
    unsigned int selected_protons         = 0;
    unsigned int selected_charged_pions   = 0;
    unsigned int signal_muons             = 0;
    unsigned int signal_protons           = 0;
    unsigned int signal_charged_pions     = 0;

    for(unsigned int i = 0; i < ev_list.size(); ++i){                                                                                                           
      Event e(ev_list[i]);

      total_true_muons         += e.CountMCParticlesWithPdg(13);
      total_true_charged_pions += e.CountMCParticlesWithPdg(211);
      total_true_charged_pions += e.CountMCParticlesWithPdg(-211);
      total_true_protons       += e.CountMCParticlesWithPdg(2212);

      total_reco_muons         += e.CountRecoParticlesWithPdg(13);
      total_reco_charged_pions += e.CountRecoParticlesWithPdg(211);
      total_reco_charged_pions += e.CountRecoParticlesWithPdg(-211);
      total_reco_protons       += e.CountRecoParticlesWithPdg(2212);

      true_muons               += GeneralAnalysisHelper::CountMatchedParticlesByTopologyTrue(e, topology, 13);
      true_charged_pions       += GeneralAnalysisHelper::CountMatchedParticlesByTopologyTrue(e, topology, 211);
      true_charged_pions       += GeneralAnalysisHelper::CountMatchedParticlesByTopologyTrue(e, topology, -211);
      true_protons             += GeneralAnalysisHelper::CountMatchedParticlesByTopologyTrue(e, topology, 2212);
      
      selected_muons           += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 13);
      selected_charged_pions   += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 211);
      selected_charged_pions   += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, -211);
      selected_protons         += GeneralAnalysisHelper::CountMatchedParticlesByTopologySelected(e, topology, 2212);
      
      signal_muons             += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 13);
      signal_charged_pions     += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 211);
      signal_charged_pions     += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, -211);
      signal_protons           += GeneralAnalysisHelper::CountMatchedParticlesByTopologySignal(e, topology, 2212);
    }

    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << "    " << topology_name                                                                                                                           << std::endl;  
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << std::setw(16) << "Proton" << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "Truth event";
    os << std::setw(16) << true_muons;
    os << std::setw(16) << true_charged_pions;
    os << std::setw(16) << true_protons;
    os << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "Selected event";
    os << std::setw(16) << selected_muons;
    os << std::setw(16) << selected_charged_pions;
    os << std::setw(16) << selected_protons;
    os << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "Signal event";
    os << std::setw(16) << signal_muons;
    os << std::setw(16) << signal_charged_pions;
    os << std::setw(16) << signal_protons;
    os << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
  }

  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(const EventList &ev_list, std::ofstream &os){
    
    unsigned int muons         = 0;
    unsigned int protons       = 0;
    unsigned int charged_pions = 0;
    unsigned int neutral_pions = 0;
    
    for(unsigned int i = 0; i < ev_list.size(); ++i){                                                                                                           
      Event e(ev_list[i]);
      muons         += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 13);
      charged_pions += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 211) + GeneralAnalysisHelper::CountMatchedParticlesAll(e, -211);
      protons       += GeneralAnalysisHelper::CountMatchedParticlesAll(e, 2212);
    }

    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << "    All events"                                                                                                                                  << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "" << std::setw(16) << "Muon" << std::setw(16) << "Charged pion" << std::setw(16) << "Proton" << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
    os << std::setw(16) << "";
    os << std::setw(16) << muons;
    os << std::setw(16) << charged_pions;
    os << std::setw(16) << protons;
    os << std::endl;
    os << "---------------------------------------------------------------------------------"                                                               << std::endl;
  
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
  
  //----------------------------------------------------------------------------------------
 
  double GeneralAnalysisHelper::Efficiency(const std::vector< double > & count_mc, const std::vector< double > & count_signal, const std::vector< double > & count_selected, const TopologyMap &topology) {
    ofstream rfile;
    rfile.open( "../Output_Selection_Tool/statistics/results.txt" );

    for( int i = 0; i<5; ++i ){
      rfile << "__________________________________________________________"                                                             << "\n";
      rfile                                                                                                                             << "\n";
      rfile << "                 TOPOLOGY NUMBER " << i                                                                                 << "\n";
      rfile << "__________________________________________________________"                                                             << "\n";
      rfile << "Count MC = "                      << count_mc[i]                                                                        << "\n";
      rfile << "Count Selected = "                << count_selected[i]                                                                  << "\n";
      rfile << "Count Signal = "                  << count_signal[i]                                                                    << "\n";
      rfile << "Background = "                    << count_selected[i] - count_signal[i]                                                << "\n";
      rfile << "Correct Reconstructed Events[%]=" << ( count_signal[i] / count_mc[i] ) * 100                                            << "\n";
      rfile << "Purity[%]="                       << (( count_signal[i] ) / count_selected[i] ) * 100                                   << "\n";
      rfile << "Background_Rejection[%]="         << (1-( count_selected[i] - count_signal[i] ) / ( count_mc[0]+count_mc[1]-count_mc[i] ) ) * 100 << "\n";
      rfile << "__________________________________________________________"                                                             << "\n";
    }

    if( topology == GeneralAnalysisHelper::GetNCTopologyMap() )    return ( count_signal[0] / count_mc[0] ) * 100;
    if( topology == GeneralAnalysisHelper::GetCCIncTopologyMap() ) return ( count_signal[1] / count_mc[1] ) * 100;
    if( topology == GeneralAnalysisHelper::GetCC0PiTopologyMap() ) return ( count_signal[2] / count_mc[2] ) * 100;
    if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) return ( count_signal[3] / count_mc[3] ) * 100;
    if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() ) return ( count_signal[4] / count_mc[4] ) * 100;
    return 0;
  }
  
  //------------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::SaveTopologyMatrix(const ParticleMatrix &count_mc_topology, const ParticleMatrix &count_signal_topology, const ParticleMatrix &count_selected_topology) {
   ofstream TMfile ;
   TMfile.open( "../Output_Selection_Tool/statistics/TopologyMatrix.txt" ) ;
    TMfile                                                                   << "\n";
    TMfile << "____________________________________________________________" << "\n";
    TMfile                                                                   << "\n";
    TMfile << "  TOPOLOGY MATRIX - TRUE RECO  (#_TReco / #_Total_MC) : "     << "\n";
    TMfile << "____________________________________________________________" << "\n";
    for( unsigned int i = 0 ; i < 5; ++i ){
      TMfile << "(";
      for( unsigned int k = 0 ; k < 5 ; ++k ) {
        if( count_signal_topology[i][k]!=0 ){
          TMfile << ( count_signal_topology[i][k] / count_mc_topology[i][k] ) * 100 << "   ,   ";
	      } 
        else TMfile << "   --   ";
      }
      TMfile <<")"                                                           << "\n";
    }
  }
  
  //------------------------------------------------------------------------------------------
 
  void GeneralAnalysisHelper::EventInformationParticles(const Event &e, std::string event_file, const int event_number) {
    ofstream efile ;
    efile.open( event_file , std::ofstream::app);

    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << " EVENT NUMBER =                                            " << event_number                       << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << "TRUE EVENTS      : "                                                                               << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << "   muons         : " << e.CountMCParticlesWithPdg(13)                                              << "\n";
    efile << "   pi+/-         : " << e.CountMCParticlesWithPdg(211) + e.CountMCParticlesWithPdg(-211)           << "\n";
    efile << "   pi0           : " << e.CountMCParticlesWithPdg(111)                                             << "\n";
    efile << "   protons       : " << e.CountMCParticlesWithPdg(2212)                                            << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << "TRUE TOPOLOGY    : "                                                                               << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    if(e.CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()))    efile << "   NC           : " << "TRUE"   << "\n";
    if(e.CheckMCTopology(GeneralAnalysisHelper::GetCCIncTopologyMap())) efile << "   ccincl.       : " << "TRUE"   << "\n";
    if(e.CheckMCTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap())) efile << "   cc0pi         : " << "TRUE"   << "\n";
    if(e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap())) efile << "   cc1pi+/-      : " << "TRUE"   << "\n";
    if(e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap())) efile << "   cc1pi0        : " << "TRUE"   << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << " SELECTED EVENTS :                                         "                                       << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << "   muons         : " << e.CountRecoParticlesWithPdg(13)                                            << "\n";
    efile << "   pi+/-         : " << e.CountRecoParticlesWithPdg(211) + e.CountRecoParticlesWithPdg(-211)       << "\n";
    efile << "   pi0           : " << e.CountRecoParticlesWithPdg(111)                                           << "\n";
    efile << "   protons       : " << e.CountRecoParticlesWithPdg(2212)                                          << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    efile << "SELECTED TOPOLOGY: "                                                                               << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";
    if(e.CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()))    efile << "   NC           : " << "TRUE" << "\n";
    if(e.CheckRecoTopology(GeneralAnalysisHelper::GetCCIncTopologyMap())) efile << "   ccincl.       : " << "TRUE" << "\n";
    if(e.CheckRecoTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap())) efile << "   cc0pi         : " << "TRUE" << "\n";
    if(e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap())) efile << "   cc1pi+/-      : " << "TRUE" << "\n";
    if(e.CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap())) efile << "   cc1pi0        : " << "TRUE" << "\n";
    efile << "-----------------------------------------------------------"                                       << "\n";

  }
  
  //------------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::EventProperties(const Event &e, const TopologyMap &topology, std::string event_file, const int event_number) {
    ofstream lfile, afile, Kfile ;
    lfile.open( event_file+="_length.txt" , std::ofstream::app);
    afile.open( event_file+="_angle.txt"  , std::ofstream::app);
    Kfile.open( event_file+="_KineticEnergy.txt" , std::ofstream::app);

    /*
    if( e.CheckMCTopology( topology ) ) { // Change the topology here
      lfile << "-----------------------------------------------------------" << "\n";
      lfile << "EVENT NUMBER      : " << event_number                        << "\n";
      lfile << "-----------------------------------------------------------" << "\n";
      lfile << "LENGTH INFORMATION: "                                        << "\n";
      lfile << "-----------------------------------------------------------" << "\n";
      lfile << "TRUE EVENTS       : "                                        << "\n";
      lfile << "-----------------------------------------------------------" << "\n";
      lfile << "Muon length         : " << GeneralAnalysisHelper::GetMCLengthWithPdg(  13  )  << "\n";
      lfile << "pi+/- length        : " << GeneralAnalysisHelper::GetMCLengthWithPdg( 211  )  << "\n";
      lfile << "pi0   length        : " << GeneralAnalysisHelper::GetMCLengthWithPdg( 111  )  << "\n";
      lfile << "p length            : " << GeneralAnalysisHelper::GetMCLengthWithPdg( 2212 )  << "\n";
      
      afile << "-----------------------------------------------------------" << "\n";
      afile << "EVENT NUMBER       : " << event_number                       << "\n";
      afile << "-----------------------------------------------------------" << "\n";
      afile << "ANGULAR INFORMATION: "                                       << "\n";
      afile << "-----------------------------------------------------------" << "\n";
      afile << "TRUE EVENTS        : "                                       << "\n"; 
      afile << "-----------------------------------------------------------" << "\n";
      afile << "Muon angle       : " << GeneralAnalysisHelper::GetMCCosThetaWithPdg(  13  )  << "\n";
      afile << "Pion angle       : " << GeneralAnalysisHelper::GetMCCosThetaWithPdg( 211  )  << "\n";
      afile << "pi0   angle      : " << GeneralAnalysisHelper::GetMCCosThetaWithPdg( 111  )  << "\n";
      afile << "Proton angle     : " << GeneralAnalysisHelper::GetMCCosThetaWithPdg( 2212 )  << "\n";

      Kfile << "-----------------------------------------------------------" << "\n";
      Kfile << "EVENT NUMBER       : " << event_number                       << "\n";
      Kfile << "-----------------------------------------------------------" << "\n";
      Kfile << "KINETIC ENERGY [GeV] : "                                     << "\n";
      Kfile << "-----------------------------------------------------------" << "\n";
      Kfile << "TRUE EVENTS        : "                                       << "\n"; 
      Kfile << "-----------------------------------------------------------" << "\n";
      Kfile << "Muon Kenergy       : " << GeneralAnalysisHelper::GetMCKineticEnergyWithPdg(  13  )  << "\n";
      Kfile << "Pion Kenergy       : " << GeneralAnalysisHelper::GetMCKineticEnergyWithPdg( 211  )  << "\n";
       
      if( e.CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy     : " << GeneralAnalysisHelper::GetMCKineticEnergyWithPdg( 2212 )  << "\n";
      if( e.CheckRecoTopology( topology ) ) {
         lfile << "-----------------------------------------------------------"     << "\n";
         lfile << "SIGNAL EVENTS     : "                                            << "\n";
         lfile << "-----------------------------------------------------------"     << "\n";
         if( !e.CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()) )lfile     << "Muon length       : " << GeneralAnalysisHelper::GetRecoLengthWithPdg(  13  ) << "\n";
         if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )lfile   << "pi+/- length      : " << GeneralAnalysisHelper::GetRecoLengthWithPdg( 211  ) << "\n";
         if( e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )lfile   << "pi0   length      : " << GeneralAnalysisHelper::GetRecoLengthWithPdg( 111  ) << "\n";
         if( e.CountMCParticlesWithPdg(2212)!=0 )lfile                                << "p length          : " << GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) << "\n";

         afile << "-----------------------------------------------------------"     << "\n";
         afile << "SIGNAL EVENTS      : "                                           << "\n"; 
         afile << "-----------------------------------------------------------"     << "\n";
         if( !e.CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap())) afile     << "Muon angle        : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg(  13  )  << "\n";
         if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) ) afile  << "Pion angle        : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 211  )  << "\n";
         if( e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) ) afile  << "Pi0 angle         : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 111  )  << "\n";
         if( e.CountMCParticlesWithPdg(2212)!=0 ) afile                               << "Proton angle      : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 2212 )  << "\n";
         
         Kfile << "-----------------------------------------------------------"     << "\n";
         Kfile << "KINETIC ENERGY [GeV] : "                                         << "\n";
         Kfile << "-----------------------------------------------------------"     << "\n";
         Kfile << "TRUE EVENTS        : "                                           << "\n"; 
         Kfile << "-----------------------------------------------------------"     << "\n";
         if( !e.CheckMCTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) Kfile    << "Muon Kenergy      : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg(  13  )  << "\n";
         if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) ) Kfile  << "Pion Kenergy      : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg( 211  )  << "\n";
         if( e.CountMCParticlesWithPdg(2212)!=0 ) Kfile                               << "Proton Kenergy    : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
      }
    }
    if( e.CheckRecoTopology( topology ) ) { // Change topology here

      lfile << "-----------------------------------------------------------"                                 << "\n";
      lfile << "SELECTED EVENTS   : "                                                                        << "\n";
      lfile << "-----------------------------------------------------------"                                 << "\n";
      if( !e.CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) lfile  << "Muon length       : "   << GeneralAnalysisHelper::GetRecoLengthWithPdg(  13  ) << "\n";
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) )lfile << "pi+/- length      : "   << GeneralAnalysisHelper::GetRecoLengthWithPdg( 211  ) << "\n";
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) )lfile << "pi0   length      : "   << GeneralAnalysisHelper::GetRecoLengthWithPdg( 111  ) << "\n";
      if( e.CountMCParticlesWithPdg(2212)!=0 ) lfile                               << "p length          : "   << GeneralAnalysisHelper::GetRecoLengthWithPdg( 2212 ) << "\n";
      
      afile << "-----------------------------------------------------------"                                 << "\n";
      afile << "SELECTED EVENTS      : "                                                                     << "\n"; 
      afile << "-----------------------------------------------------------"                                 << "\n";
      if( !e.CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) afile   << "Muon angle         : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg(  13  )  << "\n";
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) ) afile << "Pion angle         : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 211  )  << "\n";
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) ) afile << "pi0   angle        : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 111  )  << "\n";
      if( e.CountMCParticlesWithPdg(2212)!=0 ) afile                                << "Proton angle       : " << GeneralAnalysisHelper::GetRecoCosThetaWithPdg( 2212 )  << "\n";

      Kfile << "-----------------------------------------------------------"                                 << "\n";
      Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
      Kfile << "-----------------------------------------------------------"                                 << "\n";
      Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
      Kfile << "-----------------------------------------------------------"                                 << "\n";
      if( !e.CheckRecoTopology(GeneralAnalysisHelper::GetNCTopologyMap()) ) Kfile   << "Muon Kenergy       : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg(  13  )  << "\n";
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) ) Kfile << "Pion Kenergy       : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg( 211  )  << "\n";
      if( e.CountMCParticlesWithPdg(2212)!=0 ) Kfile                                << "Proton Kenergy     : " << GeneralAnalysisHelper::GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
    }
  */
  }
} // selection
