#ifndef CC1PI_ANALYSIS_HELPER_H
#define CC1PI_ANALYSIS_HELPER_H

#include <string>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "GeneralAnalysisHelper.h"
#include "EventSelectionTool.h"
#include "Event.h"
#include "Particle.h"

namespace selection{
  
  /**
   * @brief  CC1piAnalysisHelper helper class
   */
  class CC1piAnalysisHelper {

    public : 

      typedef std::vector<Particle> ParticleList;
      typedef std::vector<Event>    EventList;
      
      /**                                                              
       * @brief  Gives the number of times the Muon is the longest track and the pion or proton is the second longest track.                                     
       *
       * @param  e current event
       * @param  signal_map_topology the chosen topology 
       * @param  count_longest count the number of times the muon is the longest
       * @param  count_second_longest count the number of times the proton or pion is the second longest
       *
       * @return Particle matrix with statistics based on length and pdg
       */
      static ParticleMatrix LengthBasedStatistics(const Event &e, const TopologyMap signal_map_topology, ParticleMatrix &count_longest, ParticleMatrix &count_second_longest);
      
      /**                                                              
       * @brief Returns the MC momentum transfer for cc1pi     
       * @param pdg                        
       */
      float GetMCQ2WithPdg(const Event &e, const int pdg) const;
      /**                                                              
       * @brief Returns the Reco momentum transfer for cc1pi     
       * @param pdg                        
       */
      float GetRecoQ2WithPdg(const Event &e, const int pdg) const;

      /**                                                              
       * @brief Returns MC Energy of the particle with longest track for cc1p                         
       */
      float GetMCEnergyLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      float GetRecoEnergyLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      float GetMCEnergySecondLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      float GetRecoEnergySecondLongest(const Event &e) const;

      float GetMCKineticEnergyLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      float GetRecoKineticEnergyLongest(const Event &e) const;
      /**                                                              
       * @brief Returns MC kinetic Energy  of the particle with the second longest track for cc1p    
       */
      float GetMCKineticEnergySecondLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco kinetic Energy  of the particle with the second longest track for cc1p    
       */
      float GetRecoKineticEnergySecondLongest(const Event &e) const;
      /**                                                              
       * @brief Returns MC momentum module  of the particle with longest track for cc1p    
       */
      float GetMCModulusMomentumLongest(const Event &e) const;
      /**                                                              
       * @brief Returns Reco momentum module  of the particle with longest track for cc1p    
       */
      float GetRecoModulusMomentumLongest(const Event &e) const;
      /**                                                              
       * @brief Returns momentum module  of the particle with the second longest track for cc1p    
       */
      float GetMCModulusMomentumSecondLongest(const Event &e) const;
      /**                                                              
       * @brief Returns momentum module of the particle with the second longest track for cc1p    
       */
      float GetRecoModulusMomentumSecondLongest(const Event &e) const;

    private : 

      /**                                                              
       * @brief Returns the energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetEnergyLongest(const Event &e, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetEnergySecondLongest(const Event &e, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the kinetic energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetKineticEnergyLongest(const Event &e, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the kinetic energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetKineticEnergySecondLongest(const Event &e, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the momentum module of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      static float GetModulusMomentumLongest(const Event &e, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the momentum module of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      static float GetModulusMomentumSecondLongest(const Event &e, const ParticleList &particle_list) const;

  }; // CC1piAnalysisHelper
} // namespace: selection
#endif
