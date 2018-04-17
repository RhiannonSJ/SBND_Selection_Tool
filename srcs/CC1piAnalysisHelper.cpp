#include "../include/CC1piAnalysisHelper.h"

namespace selection{
  
  //-----------------------------------------------------------------------------------------------
  
  ParticleMatrix CC1piAnalysisHelper::LengthBasedStatistics(const Event &e, const TopologyMap signal_map_topology, ParticleMatrix &count_longest, ParticleMatrix &count_second_longest){

    if(e.CheckMCTopology(signal_map_topology) == 1){
      //LONGEST
      if(e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 2212 ) && e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 211 ) ) {count_longest[0][0]++;}
      else if(e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg( 2212 ) && e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg(13) ) {count_longest[0][1]++;}
      else if(e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg( 211 ) && e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg(13) ) {count_longest[0][2]++;}

      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1 ){
        if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) ){count_longest[1][0]++ ;}
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_longest[1][1]++;} 
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_longest[1][2]++;}
      }

      if( e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ) {
        if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) ){count_longest[2][0]++ ;}
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_longest[2][1]++;} 
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_longest[2][2]++;}
      }

      //SECOND LONGEST
      if( e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 211 ) && e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg( 2212 ) ) {count_second_longest[0][1]++;} 
      else if( e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg( 211 ) && e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg( 13 ) ) {count_second_longest[0][1]++;} 
      else if( e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg( 13 ) && e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 2212 ) ) {count_second_longest[0][0]++;} 
      else if( e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg( 13 ) && e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 211 ) ) {count_second_longest[0][0]++;} 
      else if( e.GetMCLengthWithPdg( 13 ) > e.GetMCLengthWithPdg( 2212 ) && e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg( 211 ) ) {count_second_longest[0][2]++; std::cout<<"Hello"<<std::endl;} 
      else if( e.GetMCLengthWithPdg( 211 ) > e.GetMCLengthWithPdg( 2212 ) && e.GetMCLengthWithPdg( 2212 ) > e.GetMCLengthWithPdg( 13 ) ) {count_second_longest[0][2]++; std::cout<<"Hello"<<std::endl;}
	
      if( e.CheckRecoTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1 ){
        if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[1][1]++;}
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_second_longest[1][1]++;}
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[1][0]++;}
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) ) {count_second_longest[1][0]++;}
        else if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) ) {count_second_longest[1][2]++;} 
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_second_longest[1][2]++;}
          
      }
      if( e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ) {
        if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[2][1]++;}
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) && e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_second_longest[2][1]++;}
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 13 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) ) {count_second_longest[2][0]++;}
        else if( e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) && e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 211 ) ) {count_second_longest[2][0]++;}
        else if( e.GetRecoLengthWithPdg( 13 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 211 ) ) {count_second_longest[2][2]++;} 
        else if( e.GetRecoLengthWithPdg( 211 ) > e.GetRecoLengthWithPdg( 2212 ) && e.GetRecoLengthWithPdg( 2212 ) > e.GetRecoLengthWithPdg( 13 ) ) {count_second_longest[2][2]++;} 
      }
    }
    return count_longest;
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCQ2WithPdg(const Event &e, const int pdg) const{
    return -( e.GetMCModulusMomentumWithPdg( 13 )-e.GetMCModulusMomentumWithPdg( 211 ))*(e.GetMCModulusMomentumWithPdg( 13 )-e.GetMCModulusMomentumWithPdg( 211 ));
  }

  //------------------------------------------------------------------------------------------------


  float CC1piAnalysisHelper::GetRecoQ2WithPdg(const Event &e, const int pdg) const{
    return -( e.GetRecoModulusMomentumWithPdg( 13 )-e.GetRecoModulusMomentumWithPdg( 211 ))*(e.GetRecoModulusMomentumWithPdg( 13 )-e.GetRecoModulusMomentumWithPdg( 211 ));
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetEnergyLongest(const Event &e, const ParticleList &particle_list) const {
 
    float MaxEnergy = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxEnergy = particle_list[i].GetEnergy();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCEnergyLongest(const Event &e) const {
    return GetEnergyLongest(e, e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoEnergyLongest(const Event &e) const {
    return GetEnergyLongest(e.GetRecoParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetEnergySecondLongest(const Event &e, const ParticleList &particle_list) const {
    float SecondMaxEnergy = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if (MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111){
          SecondMaxEnergy = particle_list[i].GetEnergy();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxEnergy;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCEnergySecondLongest(const Event &e) const {
    return GetEnergySecondLongest(e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoEnergySecondLongest(const Event &e) const {
    return GetEnergySecondLongest(e.GetRecoParticleList());
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetKineticEnergyLongest(const Event &e, const ParticleList & particle_list) const {
    float MaxKineticEnergy = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxKineticEnergy = particle_list[i].GetKineticEnergy();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxKineticEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCKineticEnergyLongest(const Event &e) const {
    return GetKineticEnergyLongest(e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoKineticEnergyLongest(const Event &e) const {
    return GetKineticEnergyLongest(e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetKineticEnergySecondLongest(const Event &e, const ParticleList & particle_list) const {
    float SecondMaxKineticEnergy = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
          SecondMaxKineticEnergy = particle_list[i].GetKineticEnergy();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxKineticEnergy;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCKineticEnergySecondLongest(const Event &e) const {
    return GetKineticEnergySecondLongest(e.GetMCParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoKineticEnergySecondLongest(const Event &e) const {
    return GetKineticEnergySecondLongest(e.GetRecoParticleList());
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetModulusMomentumLongest(const Event &e, const ParticleList & particle_list) const {
    float MaxMomentum = 0.;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxMomentum = particle_list[i].GetModulusMomentum();
          MaxLength = particle_list[i].GetLength();
        }
      }
    }
    return MaxMomentum;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCModulusMomentumLongest(const Event &e) const {
    return GetModulusMomentumLongest(e.GetMCParticleList());
  }


  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoModulusMomentumLongest(const Event &e) const {
    return GetModulusMomentumLongest(e.GetRecoParticleList());
  }
  
  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetModulusMomentumSecondLongest(const Event &e, const ParticleList & particle_list) const {
    float SecondMaxModulusMomentum = 0.;
    float SecondMaxLength = 0;
    float MaxLength = 0.;
    if( e.CheckMCTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
        if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
          MaxLength = particle_list[i].GetLength();
        }
        else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
          SecondMaxModulusMomentum = particle_list[i].GetModulusMomentum();
          SecondMaxLength = particle_list[i].GetLength();
        }
      }
      return SecondMaxModulusMomentum;
    }
    return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetMCModulusMomentumSecondLongest(const Event &e) const {
    return GetModulusMomentumSecondLongest(e.GetMCParticleList());
  }


  //-----------------------------------------------------------------------------------------------

  float CC1piAnalysisHelper::GetRecoModulusMomentumSecondLongest(const Event &e) const {
    return GetModulusMomentumSecondLongest(e.GetRecoParticleList());
  }


  
} // selection
