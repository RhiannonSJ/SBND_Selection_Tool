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
    if(e.CheckMCTopology(signal_map_topology) == 1) count_true++;
    if(e.CheckRecoTopology(signal_map_topology) == 1) count_selected++;
    if(e.CheckMCTopology(signal_map_topology) == 1) count_signal++;
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
        if (e.CheckMCTopology(topology_vector[i]) == 1) count_true_topology[i][j]++;
        if (e.CheckRecoTopology(topology_vector[i]) == 1) count_selected_topology[i][j]++;
        if (e.CheckMCTopology(topology_vector[i]) == 1) count_signal_topology[i][j]++;
      }
    }
    return count_signal_topology;
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
    rfile.open( "results.txt" );

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
   TMfile.open( "TopologyMatrix.txt" ) ;
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
