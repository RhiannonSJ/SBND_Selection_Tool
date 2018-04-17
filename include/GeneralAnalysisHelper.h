#ifndef GENERAL_ANALYSIS_HELPER_H
#define GENERAL_ANALYSIS_HELPER_H

#include <string>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "EventSelectionTool.h"
#include "Event.h"
#include "Particle.h"

namespace selection{
  
  /**
   * @brief  GeneralAnalysisHelper helper class
   */
  class GeneralAnalysisHelper {

    private : 
      
      TopologyMap signal_map_nc;
      TopologyMap signal_map_cc_inc;
      TopologyMap signal_map_cc_0pi;
      TopologyMap signal_map_cc_1pi;
      TopologyMap signal_map_cc_pi0;
      
    public : 

      typedef std::vector<Particle> ParticleList;
      typedef std::vector<Event>    EventList;

      /**
       * @brief  Get NC topology map
       */
      static TopologyMap GetNCTopologyMap();

      /**
       * @brief  Get CC inclusive topology map
       */
      static TopologyMap GetCCIncTopologyMap();

      /**
       * @brief  Get CC 0Pi topology map
       */
      static TopologyMap GetCC0PiTopologyMap();

      /**
       * @brief  Get CC 1Pi topology map
       */
      static TopologyMap GetCC1PiTopologyMap();

      /**
       * @brief  Get CC 1Pi0 topology map
       */
      static TopologyMap GetCCPi0TopologyMap();

      /**                                                              
       * @brief  Define topologies
       */
      static void SetTopologies();

      /**                                                              
       * @brief  Gives the number of MC, Reco and Coincidences for a given topology                                                           
       *
       * @param  e current event
       * @param  signal_map_topology chosen topology
       * @param  count_true number of true events with chosen topology
       * @param  count_signal number of signal events with chosen topology
       * @param  count selected number of selected events with chosen topology
       */
      static void TopologyStatistics(const Event &e, const TopologyMap signal_map_topology, double &count_true, double &count_signal, double &count_selected);
     
      /**                                                              
       * @brief  Obtains the Topology matrix for a specific set of events                                                                    
       *
       * @param  e current event
       * @param  count_true_topology number of true events with chosen topology
       * @param  count_signal_topology number of signal events with chosen topology
       * @param  count selected_topology number of selected events with chosen topology
       *
       * @return Matrix of topology-based statistics
       */
      static ParticleMatrix TopologyMatrix(const Event &e, ParticleMatrix &count_true_topology, ParticleMatrix &count_signal_topology, ParticleMatrix &count_selected_topology);
  


  }; // GeneralAnalysisHelper
} // namespace: selection
#endif
