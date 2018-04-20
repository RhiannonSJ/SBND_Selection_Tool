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
      
      /**                                                              
       * @brief  Get the track lengths for a given pdg
       *
       * @param  pdg
       * @param  particle_list
       * @param  lengths
       */
      static void LengthWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the track opening angles for a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  cos_thetas
       */
      static void CosThetaWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &cos_thetas);
      
      /**                                                              
       * @brief  Returns the energies of particles with a given pdg 
       *
       * @param  pdg
       * @param  particle_list
       * @param  energies
       */
      static void EnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &energies);

      /**                                                              
       * @brief  Returns the Kinetic energies of particles with a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  kinetic_energies
       */
      static void KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &kinetic_energies);

      /**                                                              
       * @brief  Returns the magnitude of the momentum of particles with a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  momentum_mod
       */
      static void ModulusMomentumWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &momentum_mod);

      
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
  
      /**                                                              
       * @brief  Get the longest MC tracks with a given pdg     
       *
       * @param  event
       * @param  pdg                                                 
       * @param  lengths 
       */
      static void GetMCLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the longest reconstructed tracks with a given pdg
       *
       * @param  event
       * @param  pdg                                                    
       * @param  lengths 
       */
      static void GetRecoLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the cos thetas for given MC pdg
       *
       * @param  event
       * @param  pdg
       * @param  cos_thetas
       */
      static void GetMCCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas);

      /**                                                              
       * @brief  Get the cos theta with reco pdg
       *
       * @param  event
       * @param  pdg
       * @param  cos_thetas
       */
      static void GetRecoCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas);
      
      /**                                                              
       * @brief  Get the MC energies with a given pdg
       *
       * @param  event
       * @param  pdg
       * @param  energies
       */
      static void GetMCEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies);

      /**                                                              
       * @brief  Get the reco energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  energies
       */
      static void GetRecoEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies);

      /**                                                              
       * @brief  Get the MC kinetic energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  kinetic_energies
       */
      static void GetMCKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies);

      /**                                                              
       * @brief  Get the reco kinetic energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  kinetic_energies
       */
      static void GetRecoKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies);
      
      /**                                                              
       * @brief  Get the MC modulus of the momentum with pdg
       *
       * @param  event
       * @param  pdg
       * @param  momentum_mod
       */
      static void GetMCModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod);
      
      /**                                                              
       * @brief  Get the reco modulus of the momentum with pdg
       *
       * @param  event
       * @param  pdg
       * @param  momentum_mod
       */
      static void GetRecoModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod);
      
      /**
       * @brief  Calculates the Efficiency, Purity, Background Rejection
       *         Parameters for Efficiency calculation ( MC, signal and selected ) for a given topology : 
       *         0-> No muon, 
       *         1 -> CCinclusive,
       *         2-> CC0pi, 
       *         3-> CC1pi+/-,
       *         4-> CC1pi0
       *
       * @param  count_mc
       * @param  count_signal 
       * @param  count_selected 
       **/
      static double Efficiency(const std::vector< double > &count_mc, const std::vector< double > &count_signal, const std::vector< double > &count_selected, const TopologyMap &topology);

      /**
       * @brief  Save Topology Matrix into a file
       *         BackGround Study : topology mis identification table 
       *         0-> No muon, 
       *         1 -> CCinclusive,
       *         2-> CC0pi, 
       *         3-> CC1pi+/-,
       *         4-> CC1pi0
       *
       * @param  count_mc_topology 
       * @param  count_signal_topology 
       * @param  count_selected_topology 
       **/
      static void SaveTopologyMatrix(const ParticleMatrix &count_mc_topology, const ParticleMatrix &count_signal_topology, const ParticleMatrix &count_selected_topology);

      /**
       * @brief Saves the event information in a file ( types of particles in the event 
       * and topology for True and Selected
       **/
      static void EventInformationParticles(const Event &e, const std::string name, const int event_number);

      /**
       * @brief Saves the event characteristics in a file ( length , angle and kinetic energy ) 
       * for the selected topology
       **/
      static void EventProperties(const Event &e, const TopologyMap &topology, std::string event_file, const int event_number);

  }; // GeneralAnalysisHelper
} // namespace: selection
#endif
