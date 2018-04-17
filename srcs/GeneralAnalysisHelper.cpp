#include "../include/GeneralAnalysisHelper.h"

namespace selection{

  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetNCTopologyMap() {return signal_map_nc;} 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCCIncTopologyMap() {return signal_map_cc_inc;} 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC0PiTopologyMap() {return signal_map_CC_0pi;} 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetCC1PiTopologyMap() {return signal_map_cc_1pi;} 
  //----------------------------------------------------------------------------------------
  TopologyMap GeneralAnalysisHelper::GetPi0TopologyMap() {return signal_map_cc_pi0;} 
  //----------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::SetTopologies() {

    // Particle pgd code vectors 
    std::vector< int > mu, pi, pi0, p;
    mu.push_back(13);
    pi.push_back(211);
    pi.push_back(-211);
    pi0.push_back(111);
    p.push_back(2212);
    
    /*                                                               
     * NC topology  ( Topology number 0 )                            
     * TopologyMap signal_map_nc;                                    
     */
    signal_map_nc.insert(std::make_pair(mu, 0));

    /*                                                                 
     * CC inclusive topology ( Topology number 1 )                     
     * TopologyMap signal_map_cc_inc;                            
     */
    signal_map_cc_inc.insert(std::make_pair(mu, 1));

    /*                                                                 
     * CC0pi topology ( Topology number 2 )                            
     * TopologyMap signal_map_cc_0pi;                                  
     */
    signal_map_cc_0pi.insert(std::make_pair(mu,  1));
    signal_map_cc_0pi.insert(std::make_pair(pi,  0));
    signal_map_cc_0pi.insert(std::make_pair(pi0, 0));

    /*                                                                 
     * CC1pi+/- topology  ( Topology number 3 )                        
     * TopologyMap signal_map_cc_1pi;                                  
     */
    signal_map_cc_1pi.insert(std::make_pair(mu, 1));
    signal_map_cc_1pi.insert(std::make_pair(pi, 1));
    
    /*                                                                 
     * CC1pi0 topology ( Topology number 4 )                           
     * TopologyMap signal_map_cc_pi0;                                  
     */
    signal_map_cc_pi0.insert(std::make_pair(mu,  1));
    signal_map_cc_pi0.insert(std::make_pair(pi0, 1));
  }
  
  //-------------------------------------------------------------------------------------------------
  
  void GeneralAnalysisHelper::TopologyStatistics(const Event &e, const TopologyMap signal_map_topology, double & count_true, double & count_signal, double & count_selected){
    if(e.CheckMCTopology(signal_map_topology) == 1) count_true++;
    if(e.CheckRecoTopology(signal_map_topology) == 1 && e.GetRecoCC0piNeutrinoEnergy() > 0 && e.GetRecoCC1piNeutrinoEnergy() > 0) count_selected++;
    if(e.CheckMCTopology(signal_map_topology) == 1 && e.CheckRecoTopology(signal_map_topology) == 1 && e.GetRecoCC0piNeutrinoEnergy() > 0 && e.GetRecoCC1piNeutrinoEnergy() > 0) count_signal++;
  }

  //----------------------------------------------------------------------------------------
  
  ParticleMatrix GeneralAnalysisHelper::TopologyMatrix(const Event &e, ParticleMatrix &count_true_topology, ParticleMatrix &count_signal_topology, ParticleMatrix &count_selected_topology){
    std::vector<TopologyMap> topology_vector(5);
    topology_vector[0] = signal_map_nc;
    topology_vector[1] = signal_map_cc_inc;
    topology_vector[2] = signal_map_cc_0pi;
    topology_vector[3] = signal_map_cc_1pi;
    topology_vector[4] = signal_map_cc_pi0;

    for(unsigned int i=0; i < topology_vector.size(); ++i ){
      for(unsigned int j=0; j < topology_vector.size(); ++j ){
        if (e.CheckMCTopology(topology_vector[i]) == 1) count_true_topology[i][j]++;
        if (e.CheckRecoTopology(topology_vector[i]) == 1 && e.GetRecoCC0piNeutrinoEnergy() > 0 && e.GetRecoCC1piNeutrinoEnergy() > 0) count_selected_topology[i][j]++;
        if (e.CheckMCTopology(topology_vector[i]) == 1 && e.CheckRecoTopology(topology_vector[j]) == 1 && e.GetRecoCC0piNeutrinoEnergy() > 0 && e.GetRecoCC1piNeutrinoEnergy() > 0) count_signal_topology[i][j]++;
      }
    }
    return count_signal_topology;
  }
 
  //-----------------------------------------------------------------------------------------------
  

} // selection
