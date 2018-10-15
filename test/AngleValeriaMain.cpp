 #include "../include/CC0piAnalysisHelper.h"                                                    
 #include "../include/GeneralAnalysisHelper.h"                                                  
 #include "../include/EventSelectionTool.h"                                                     
 #include "../include/Event.h"                                                                  
 #include <iostream>                                                                            
 #include <sstream>                                                                             
 #include <numeric>                                                                             
 #include <time.h>                                                                              
 #include <algorithm>                                                                           
 #include "TVector3.h"                                                                          
 #include "TH1.h"                                                                               
 #include "TH2.h"                                                                               
 #include "TCanvas.h"                                                                           
 #include "TLegend.h"                                                                           
 #include "TLatex.h"                                                                            
 #include "TStyle.h"                                                                            
 #include "TColor.h"                                                                            
 #include "TObjArray.h"                                                                         
 #include "THStack.h"
 #include "TProfile.h"

 using namespace selection;                                                                     
                                                                                                
 int MainTest(){ 
                                                                                                
   time_t rawtime;                                                                              
   struct tm * timeinfo;                                                                        
   time (&rawtime);                                                                             
   timeinfo = localtime (&rawtime);                                                             
   std::cout << "-----------------------------------------------------------" << std::endl;     
   std::cout << " Start local time and date:  " << asctime(timeinfo)          << std::endl;     
   std::cout << "-----------------------------------------------------------" << std::endl;     
                                                                                                
   // Output file location                                                                      
   std::string momenta = "../Output_Selection_Tool/plots/angle/momenta_distributions/";   
   std::string stats_location= "../Output_Selection_Tool/statistics/";
   std::string lengthpath = "../Output_Selection_Tool/plots/angle/length_distributions/";
   std::string cosinetheta="../Output_Selection_Tool/plots/angle/costheta_distributions/";
   
   //Load events                                         
   
   // Initialise event list and the topology maps                                               
   EventSelectionTool::EventList events;
   int start = static_cast<int>(time(NULL));                                                    
   unsigned int total = 500;                                                                   
                                                                                               
   // Load the events into the event list                                                       
   for( unsigned int i = 0; i < total; ++i ){                                                   
                                                                                                
     std::string name;                                                                          
     name.clear();       
     char file_name[1024];                                                                      
     name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
     //name = "/hepstore/rjones/Samples/FNAL/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
     //name = "/pnfs/sbnd/persistent/users/rsjones/analysis_sample/120918_ana_files/11509725_"+std::to_string(i)+"/output_file.root";
     strcpy( file_name, name.c_str() );                                                         
    
     EventSelectionTool::LoadEventList(file_name, events, i);                                                                                                                              
     EventSelectionTool::GetTimeLeft(start,total,i);                                            
   }                                                                                            
   std::cout << std::endl;  

   TopologyMap cc_signal_map         = GeneralAnalysisHelper::GetCCIncTopologyMap(); 
   TopologyMap cc0pi2p_signal_map    = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap(); 

std::vector <float> signal_costheta_mu, signal_costheta_p1, signal_costheta_p2, signal_costheta_p1p2;
std::vector <float> reco_costheta_mu, reco_costheta_p1, reco_costheta_p2, reco_costheta_p1p2;
std::vector <float> true_costheta_mu, true_costheta_p1, true_costheta_p2, true_costheta_p1p2; 
std::vector <float> signal_momentum_mu, signal_momentum_p1, signal_momentum_p2;                
std::vector <float> reco_momentum_mu, reco_momentum_p1, reco_momentum_p2;                      
std::vector <float> true_momentum_mu, true_momentum_p1, true_momentum_p2;
std::vector <float> signal_length_mu, signal_length_p1, signal_length_p2;                
std::vector <float> reco_length_mu, reco_length_p1, reco_length_p2;                      
std::vector <float> true_length_mu, true_length_p1, true_length_p2;


int tottrue      = 0; // counter of total CC true events
int totreco      = 0; // counter of total CC reconstructed events
int cc0pi2pmc    = 0; // counter of cc0pi2p true events
int cc0pi2pre    = 0; // counter of cc0pi2p reconstructed events  
int cc0pi2ps1    = 0; // counter of cc0pi2p signal events by muons
int cc0pi2ps2    = 0; // counter of cc0pi2p signal events by protons
int mutrue       = 0; // counter of true muons
int mureco       = 0; // counter of reconstructed particles with escaping tracks 
int musignal     = 0; // counter of signal muons
int nproton      = 0; // counter for comparing true proton momenta
int nprotonreco  = 0; // counter for comparing reco proton momenta
int nprotonsig   = 0; // counter for comparing signal proton momenta
int mucount      = 0; // counter for true muons
int mucountr     = 0; // counter for reco muons
int mucounts     = 0; // counter for signal muons
int ptotreco     = 0; // counter of reconstructed particles with not escaping & w/ at least 5 hits tracks 
int ptottrue     = 0; // counter of true protons 
int ptotsig      = 0; // counter of signal protons
int truev        = 0; // counter of true events considering both protons and muons  
int ntotevents   = 0; // counter of total number of events 
int recopevent   = 0; // counter of total number of reconstructed events
float momento; 
float momento2;
float momento3; 
float costheta;
float costheta2;
float costheta3;
float length;
float length2;
float length3; 
std:: vector <float> mo;
std:: vector <float> res;
std:: vector <float> zezi;
std:: vector <float> mores, mores2, mores3;
ofstream file;
file.open(stats_location+"costheta.txt");

// loop over the events 

for(const Event &e : events){
  if(e.IsSBNDTrueFiducial()){
    if(bool max = GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
    ++ntotevents;
    ParticleList true_particles = e.GetMCParticleList();  
    ParticleList reco_particles = e.GetRecoParticleList();                                    
    if(e.CheckRecoTopology(cc_signal_map)){
      ++totreco;
    }
    if(e.CheckMCTopology(cc_signal_map)){
    ++tottrue;
    }

// cc0pi2P true events

    if(e.CheckMCTopology(cc0pi2p_signal_map)){
      ++cc0pi2pmc;  
      
      //loop over the true particles
  
      for(Particle &p_true : true_particles){
        if(p_true.GetPdgCode() == 13){
          ++mutrue; 
          ++mucount;
          momento3=p_true.GetModulusMomentum();
          length3=p_true.GetLength();
          costheta3=p_true.GetCosTheta();
        }
        else if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021){
          if(p_true.GetModulusMomentum()>0.00001 && p_true.GetModulusMomentum() < 100) {
            ++nproton; 
            ++ptottrue;
            mo.push_back(p_true.GetModulusMomentum());
            res.push_back(p_true.GetCosTheta());   
            mores.push_back(p_true.GetLength());
          }
        }
      }

      true_momentum_mu.push_back(momento3);
      true_costheta_mu.push_back(costheta3);
      true_length_mu.push_back(length3);
      if(mo[0]<mo[1] && mo[1]!=0 && mo[0]!=0){
       true_costheta_p1.push_back(res[0]);
       true_costheta_p2.push_back(res[1]);
       true_momentum_p1.push_back(mo[0]);
       true_momentum_p2.push_back(mo[1]);
       true_length_p1.push_back(mores[0] );
       true_length_p2.push_back(mores[1]);
      }
       else if(mo[1]<mo[0] && mo[1]!=0 && mo[0]!= 0){
       true_costheta_p1.push_back(res[1]);
       true_costheta_p2.push_back(res[0]);
       true_momentum_p1.push_back(mo[1]);
       true_momentum_p2.push_back(mo[0]);
       true_length_p1.push_back(mores[1] );
       true_length_p2.push_back(mores[0]);
       }
       mo.clear();
       res.clear();
       mores.clear();
       nproton=0;
       mucount=0;
    } 
   
    // CC0Pi2P reconstructed events
 
    if(e.CheckRecoTopology(cc0pi2p_signal_map)){
      ++cc0pi2pre;
      if(unsigned int k=e.CountRecoParticlesWithPdg(2212)==2){
  
        // loop over the reconstructed particles  
     
        for(Particle &p_reco : reco_particles){        
          if(p_reco.GetFromRecoTrack()){
            if(p_reco.GetPdgCode()==13){         
              ++mureco;
              ++mucountr;
              momento2=p_reco.GetModulusMomentum();
              length2=p_reco.GetLength();
             costheta2= p_reco.GetCosTheta();
           
              //get the MC particle by hits 
           
              int phits = p_reco.GetMCParticleIdHits();           
              for(Particle &p_true : true_particles){
                int mcid = p_true.GetMCId();  
                if(p_true.GetPdgCode()==13){ 
                  if(mcid==phits){             
                    ++musignal;
                    ++mucounts;
                    costheta=p_reco.GetCosTheta();
                    momento=p_reco.GetModulusMomentum();
                    length=p_reco.GetLength();
                  }
                }
              }
            }
            else if(p_reco.GetPdgCode()==2212 && p_reco.GetKineticEnergy()>=0.021){
              ++ptotreco;
              ++nprotonreco; 
              //           if(p_reco.GetModulusMomentum()>=0.250){ 
              mo.push_back(p_reco.GetModulusMomentum());
              res.push_back(p_reco.GetCosTheta());
              mores3.push_back(p_reco.GetLength());
              //get the MC particle by hits 
              int  phits = p_reco.GetMCParticleIdHits();           
              for(Particle &p_true : true_particles){
                int mcid = p_true.GetMCId();  
                if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021){ 
                  if(mcid==phits){
                    ++nprotonsig;
                    ++ptotsig;
                    zezi.push_back(p_reco.GetModulusMomentum());         
                    mores.push_back(p_reco.GetCosTheta());
                    mores2.push_back(p_reco.GetLength());
                      
                  } 
                }              
              }
            }
            }
    } //closes the loop on reco particles
    reco_length_mu.push_back(length2);
     reco_costheta_mu.push_back(costheta2);
     reco_momentum_mu.push_back(momento2); 
       if(mo[0]<mo[1]){
       reco_costheta_p1.push_back(res[0]);
       reco_costheta_p2.push_back(res[1]);
       reco_momentum_p1.push_back(mo[0]);
       reco_momentum_p2.push_back(mo[1]);
       reco_length_p1.push_back(mores3[0]);
       reco_length_p2.push_back(mores3[1]);
       }
       else if(mo[1]<=mo[0]){
       reco_costheta_p1.push_back(res[1]);
       reco_costheta_p2.push_back(res[0]);
       reco_momentum_p1.push_back(mo[1]);
       reco_momentum_p2.push_back(mo[0]);
       reco_length_p1.push_back(mores3[1]);
       reco_length_p2.push_back(mores3[0]);
       }
      
     if(nprotonsig==2 && mucounts==1){
        ++cc0pi2ps2;
        signal_costheta_mu.push_back(costheta);
        signal_momentum_mu.push_back(momento);
        signal_length_mu.push_back(length);
      if(zezi[0]<=zezi[1]){
       signal_costheta_p1.push_back(mores[0]);
       signal_costheta_p2.push_back(mores[1]);
       signal_momentum_p1.push_back(zezi[0]);                                                
       signal_momentum_p2.push_back(zezi[1]); 
       signal_length_p1.push_back(mores2[0]);
       signal_length_p2.push_back(mores2[1]);
      }
       else if(zezi[1]<zezi[0]){
       signal_momentum_p1.push_back(zezi[1]);
       signal_momentum_p2.push_back(zezi[0]);
       signal_costheta_p1.push_back(mores[1]);
       signal_costheta_p2.push_back(mores[0]); 
       signal_length_p1.push_back(mores2[1]);
       signal_length_p2.push_back(mores2[0]);
       }
      }
     res.clear(); 
     nprotonreco=0;
     mucountr=0;
     zezi.clear();
     mores.clear();
     nprotonsig=0;
     mucounts=0;
     mores2.clear();
     mores3.clear();
  }
  }  
 } // end of the condition on reconstructed topology map
}
} // end of the loop on the event         

// writing on the file all the useful information
file << "il numero di eventi generati è :      " << ntotevents << endl;
file << endl;
file << "==================== EVENTI SIMULATI =================" << endl;
file << endl;
file << "Topologia CC " << tottrue << endl;
file << "Topologia CC0pi2p  " << cc0pi2pmc << endl;
file << endl;
file << " ======= Particelle prodotte negli eventi simulati ====== " << endl;
file << endl;
file << "Muoni " << mutrue << endl; 
file << "Protoni " << ptottrue << endl;
file << endl;
file << "======================================================" << endl;
file << endl;
file << "================ EVENTI RICOSTRUITI ==================" << endl;
file << endl;
file << "Topologia CC " << totreco << endl;
file << "Topologia CC0pi2p " << cc0pi2pre << endl; 
file << "Eventi contati " << recopevent << endl;
file << endl; 
file << "===== Particelle prodotte negli eventi ricostruiti ===== " << endl; 
file << endl;
file << "Muoni" << mureco << endl;
file << "Protoni " << ptotreco << endl;
file << "======================================================" << endl;
file << endl;
file << "================== EVENTI DI SEGNALE ==================" << endl;
file << endl;
file << "Eventi contati" << cc0pi2ps2 << endl;
file << endl; 
file << " ===== Particelle contate negli eventi di segnale ===== " << endl;
file << endl;
file << "Muoni " << musignal << endl; 
file << "Protoni " << ptotsig << endl; 
file << "=======================================================" << endl;
file << endl; 
file << "_____ true events ____" << endl;
file << "dimensioni impulsi theta mu, p1 e p2 " << true_momentum_mu.size() << "   " << true_momentum_p1.size() << "     " << true_momentum_p2.size() << endl;
file << "dimensioni cos theta mu, p1 e p2 " << true_costheta_mu.size();
file <<"            " <<  true_costheta_p1.size(); 
file << "            " << true_costheta_p2.size() << endl;
for(unsigned int i=0; i<true_costheta_p1.size(); i++){
  file <<  "impulso 1" << true_momentum_p1[i]  << " cos theta 1 = " << true_costheta_p1[i] << endl; 
  file << " impulso 2" << true_momentum_p2[i] << " cos theta 2 = " << true_costheta_p2[i] << endl; 
}

file << "_____ reco events ____"<< endl;
file << "dimensioni impulsi theta mu, p1 e p2 " << reco_momentum_mu.size() << "   " << reco_momentum_p1.size() << "     " << reco_momentum_p2.size() << endl;
file << "dimensioni cos theta mu, p1 e p2 " << reco_costheta_mu.size();
file <<"            " <<  reco_costheta_p1.size(); 
file << "            " << reco_costheta_p2.size() << endl;
for(unsigned int i=0; i<reco_costheta_p1.size(); i++){
  file << "impulso " << reco_momentum_mu[i] << " cos theta mu = " << reco_costheta_mu[i] << endl;
  file <<  "impulso " << reco_momentum_p1[i] << " cos theta 1 = " << reco_costheta_p1[i] << endl; 
  file << "impulso " << reco_momentum_p2[i] <<" cos theta 2 = " << reco_costheta_p2[i] << endl; 
}

file << "_____ signal events ____"<< endl;
file << "dimensioni impulsi theta mu, p1 e p2 " << signal_momentum_mu.size() << "   " << signal_momentum_p1.size() << "     " << signal_momentum_p2.size() << endl;
file << "dimensioni cos theta mu, p1 e p2 " << signal_costheta_mu.size();
file <<"            " <<  signal_costheta_p1.size(); 
file << "            " << signal_costheta_p2.size() << endl;
for(unsigned int i=0; i<signal_costheta_p1.size(); i++){
  file << "impulso " << signal_momentum_mu[i] << " cos theta mu = " << signal_costheta_mu[i] << endl;
  file <<  "impulso " << signal_momentum_p1[i] << " cos theta 1 = " << signal_costheta_p1[i] << endl; 
 file << "impulso " << signal_momentum_p2[i] <<" cos theta 2 = " << signal_costheta_p2[i] << endl; 
}


/*
//efficiency and purity
//
double eff;
double pur; 

eff = ((double) cc0pi2ps2)/(truev); 
pur = ((double) cc0pi2ps2)/(recopevent);

file << "Efficienza = " << eff << endl;
file << endl;
file << "Purity = " << pur << endl;
file << endl; 
*/


// fill histograms

TH1D *h_reco_mu = new TH1D("h_reco_mu", "Muon momentum distribution", 40, 0 ,5);                                        
TH1D *h_true_mu = new TH1D("h_true_mu", "Muon momentum distribution" , 40, 0, 5);                                    
TH1D *h_true_p1 = new TH1D("h_true_p1", "p1 momentum distribution", 40, 0, 2);                                       
TH1D *h_reco_p1 = new TH1D("h_reco_p1", "p1 momentum distribution", 40, 0, 2);                                         
TH1D *h_true_p2 = new TH1D("h_true_p2", "p2 momentum distribution", 40, 0, 2);                                       
TH1D *h_reco_p2 = new TH1D("h_reco_p2", "p2 momentum distribution", 40, 0, 2);                                       
TH1D *h_sign_mu = new TH1D("h_sign_mu","Muon momentum distribution", 40,0,5);                                            
TH1D *h_sign_p1 = new TH1D("h_sign_p1", "p1 momentum distribution", 40, 0, 2);                                       
TH1D *h_sign_p2 = new TH1D("h_sign_p2", "p2 momentum distribution", 40, 0, 2);                                       
for(unsigned int i=0; i<reco_momentum_mu.size(); i++){                                         
  h_reco_mu->Fill(reco_momentum_mu[i]);                                                        
}                                                                                              
for(unsigned int i=0; i<true_momentum_mu.size(); i++){                                        
  h_true_mu->Fill(true_momentum_mu[i]);                                                       
}                                                                                             
for(unsigned int i=0; i<true_momentum_p1.size(); i++){                                         
  h_true_p1->Fill(true_momentum_p1[i]);                                                         
}                                                                                             
for(unsigned int i=0; i<reco_momentum_p1.size(); i++){                                         
  h_reco_p1->Fill(reco_momentum_p1[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<true_momentum_p2.size(); i++){                                         
  h_true_p2->Fill(true_momentum_p2[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<reco_momentum_p2.size(); i++){                                         
  h_reco_p2->Fill(reco_momentum_p2[i]);                                                         
}                                                                                              
for(unsigned int i=0; i<signal_momentum_mu.size(); i++){                                      
  h_sign_mu->Fill(signal_momentum_mu[i]);                                                      
}                                                                                              
for(unsigned int i=0; i<signal_momentum_p1.size(); i++){                                      
  h_sign_p1->Fill(signal_momentum_p1[i]);                                                      
}                                                                                 
for(unsigned int i=0; i<signal_momentum_p2.size(); i++){                                       
  h_sign_p2->Fill(signal_momentum_p2[i]);                                                      
} 

TH1D *a_reco_mu = new TH1D("a_reco_mu", "Muon cosine distribution", 40, -1 , 1);
TH1D *a_true_mu = new TH1D("a_true_mu", "Muon cosine distribution"  , 40, -1, 1);
TH1D *a_true_p1 = new TH1D("a_true_p1", " p1 cosine distribution", 40, -1, 1);
TH1D *a_reco_p1 = new TH1D("a_reco_p1", "p1 cosine distribution", 40, -1, 1);
TH1D *a_true_p2 = new TH1D("a_true_p2", " p2 cosine distribution", 40, -1, 1);
TH1D *a_reco_p2 = new TH1D("a_reco_p2", "p2 cosine distribution", 40, -1, 1);
TH1D *a_sign_mu = new TH1D("a_sign_mu", "Muon cosine distribution", 40, -1, 1);
TH1D *a_sign_p1 = new TH1D("a_sign_p1", "p1 cosine distribution", 40, -1, 1);
TH1D *a_sign_p2 = new TH1D("a_sign_p2", "p2 cosine distribution", 40, -1, 1);

for(unsigned int i=0; i<reco_costheta_mu.size(); i++){
  a_reco_mu->Fill(reco_costheta_mu[i]);
}
for(unsigned int i=0; i<true_costheta_mu.size(); i++){
  a_true_mu->Fill(true_costheta_mu[i]);
}
for(unsigned int i=0; i<true_costheta_p1.size(); i++){
 a_true_p1->Fill(true_costheta_p1[i]);
}
for(unsigned int i=0; i<reco_costheta_p1.size(); i++){
 a_reco_p1->Fill(reco_costheta_p1[i]);
}
for(unsigned int i=0; i<true_costheta_p2.size(); i++){
 a_true_p2->Fill(true_costheta_p2[i]);
}
for(unsigned int i=0; i<reco_costheta_p2.size(); i++){
 a_reco_p2->Fill(reco_costheta_p2[i]);
}
for(unsigned int i=0; i<signal_costheta_mu.size(); i++){
 a_sign_mu->Fill(signal_costheta_mu[i]);
}
for(unsigned int i=0; i<signal_costheta_p1.size(); i++){
 a_sign_p1->Fill(signal_costheta_p1[i]);
}

for(unsigned int i=0; i<signal_costheta_p2.size(); i++){
 a_sign_p2->Fill(signal_costheta_p2[i]);
}



TH1D *l_reco_mu = new TH1D("l_reco_mu", "Muon length distribution", 40, 0 , 40);
TH1D *l_true_mu = new TH1D("l_true_mu", "Muon length distribution"  , 40, 0, 40);
TH1D *l_true_p1 = new TH1D("l_true_p1", " p1 length distribution", 40, 0, 40);
TH1D *l_reco_p1 = new TH1D("l_reco_p1", "p1 length distribution", 40, 0, 40);
TH1D *l_true_p2 = new TH1D("l_true_p2", " p2 length distribution", 40, 0, 40);
TH1D *l_reco_p2 = new TH1D("l_reco_p2", "p2 length distribution", 40, 0, 40);
TH1D *l_sign_mu = new TH1D("l_sign_mu", "Muon length distribution", 40, 0, 40);
TH1D *l_sign_p1 = new TH1D("l_sign_p1", "p1 length distribution", 40, 0, 40);
TH1D *l_sign_p2 = new TH1D("l_sign_p2", "p2 length distribution", 40, 0, 40);

for(unsigned int i=0; i<reco_length_mu.size(); i++){
  l_reco_mu->Fill(reco_length_mu[i]);
}
for(unsigned int i=0; i<true_length_mu.size(); i++){
  l_true_mu->Fill(true_length_mu[i]);
}
for(unsigned int i=0; i<true_length_p1.size(); i++){
 l_true_p1->Fill(true_length_p1[i]);
}
for(unsigned int i=0; i<reco_length_p1.size(); i++){
 l_reco_p1->Fill(reco_length_p1[i]);
}
for(unsigned int i=0; i<true_length_p2.size(); i++){
 l_true_p2->Fill(true_length_p2[i]);
}
for(unsigned int i=0; i<reco_length_p2.size(); i++){
 l_reco_p2->Fill(reco_length_p2[i]);
}
for(unsigned int i=0; i<signal_length_mu.size(); i++){
 l_sign_mu->Fill(signal_length_mu[i]);
}
for(unsigned int i=0; i<signal_length_p1.size(); i++){
 l_sign_p1->Fill(signal_length_p1[i]);
}

for(unsigned int i=0; i<signal_length_p2.size(); i++){
 l_sign_p2->Fill(signal_length_p2[i]);
}


file << "=====================================================================" << endl;
file << "enties for the histogram of true mu " << h_true_mu->GetEntries() << endl;
file << "entries for the histogram of reco mu " << h_reco_mu->GetEntries() << endl;
file << "entries for the histogram of signal mu " << h_sign_mu->GetEntries() << endl;
file << "=====================================================================" << endl;
file << "=====================================================================" << endl;
file << "enties for the histogram of true mu costheta" << a_true_mu->GetEntries() << endl;
file << "entries for the histogram of reco mu costheta " << a_reco_mu->GetEntries() << endl;
file << "entries for the histogram of signal mu costheta" << a_sign_mu->GetEntries() << endl;
file << "=====================================================================" << endl;
file << "=====================================================================" << endl;
file << "entries for the histogram of true p1 " << h_true_p1->GetEntries() << endl;
file << " entries for the histogram of reco p1 " << h_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p1 " << h_sign_p1->GetEntries() << endl; 
file << "=====================================================================" << endl;
file << "enties for the histogram of true p1 costheta" << a_true_p1->GetEntries() << endl;
file << "entries for the histogram of reco p1 costheta " << a_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p1 costheta" << a_sign_p1->GetEntries() << endl;
file << "=====================================================================" << endl;
file << "=====================================================================" << endl;
file << "entries for the histogram of true p2 " << h_true_p1->GetEntries() << endl;
file << " entries for the histogram of reco p2 " << h_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p2 " << h_sign_p1->GetEntries() << endl; 
file << "=====================================================================" << endl;
file << "enties for the histogram of true p2 costheta" << a_true_p1->GetEntries() << endl;
file << "entries for the histogram of reco p2 costheta " << a_reco_p1->GetEntries() << endl;
file << "entries for the histogram of signal p2 costheta" << a_sign_p1->GetEntries() << endl;
file << "=====================================================================" << endl;




// Drawing

TCanvas *c = new TCanvas();                                                                    
gStyle->SetOptStat(0);
h_reco_mu->SetLineColorAlpha(kBlue, 0.30);
h_reco_mu->SetLineWidth(2);
h_reco_mu->SetFillColorAlpha(kBlue, 0.20);
h_reco_mu->SetFillStyle(3001);
h_true_mu->SetLineColorAlpha(kRed, 0.30);
h_true_mu->SetLineWidth(2);
h_true_mu->SetFillColorAlpha(kRed, 0.20);
h_true_mu->SetFillStyle(3001);
h_sign_mu->SetLineColorAlpha(kGreen, 0.30);
h_sign_mu->SetLineWidth(2);
h_sign_mu->SetFillColorAlpha(kGreen, 0.20);
h_sign_mu->SetFillStyle(3001);
h_true_mu->Scale(1/(h_true_mu->GetEntries()));                                                 
h_reco_mu->Scale(1/(h_reco_mu->GetEntries()));                                                 
h_sign_mu->Scale(1/(h_sign_mu->GetEntries()));                                                 
h_sign_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_reco_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_true_mu->GetYaxis()->SetRangeUser(0, 0.35);                                                   
h_sign_mu->GetYaxis()->SetTitle("Events");                                                     
h_sign_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_reco_mu->GetYaxis()->SetTitle("Events");                                                     
h_reco_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_true_mu->GetYaxis()->SetTitle("Events");                                                     
h_true_mu->GetXaxis()->SetTitle("p_{\\mu} (GeV/c)");                                             
h_reco_mu->Draw("HIST");                                                                       
h_true_mu->Draw("HIST SAME");                                                                  
h_sign_mu->Draw("HIST SAME"); 
TLegend *legend = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend->AddEntry(h_true_mu,"true","f");                                                            
legend->AddEntry(h_reco_mu, "reco", "f");
legend->AddEntry(h_sign_mu, "signal", "f");
legend->Draw();  
c->Update();
c->SaveAs((momenta+"Muon_Histogram_momentum.root").c_str());                                       
c->SaveAs((momenta+"Muon_Histogram_momentum.png").c_str());                                        

TCanvas *a1 = new TCanvas();                                                                    
gStyle->SetOptStat(0);    
a_reco_mu->SetLineColorAlpha(kBlue, 0.30);
a_reco_mu->SetLineWidth(2);
a_reco_mu->SetFillColorAlpha(kBlue, 0.20);
a_reco_mu->SetFillStyle(3001);
a_true_mu->SetLineColorAlpha(kRed, 0.30);
a_true_mu->SetLineWidth(2);
a_true_mu->SetFillColorAlpha(kRed, 0.20);
a_true_mu->SetFillStyle(3001);
a_sign_mu->SetLineColorAlpha(kGreen, 0.30);
a_sign_mu->SetLineWidth(2);
a_sign_mu->SetFillColorAlpha(kGreen, 0.20);
a_sign_mu->SetFillStyle(3001);
a_true_mu->Scale(1/(a_true_mu->GetEntries()));                                                 
a_reco_mu->Scale(1/(a_reco_mu->GetEntries()));                                                 
a_sign_mu->Scale(1/(a_sign_mu->GetEntries()));                                                 
a_sign_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_reco_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_true_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_sign_mu->GetYaxis()->SetTitle("Events");                                                     
a_sign_mu->GetXaxis()->SetTitle("cos \\theta_{\\mu} ");                                             
a_reco_mu->GetYaxis()->SetTitle("Events");                                                     
a_reco_mu->GetXaxis()->SetTitle("cos \\theta_{\\mu} ");                                             
a_true_mu->GetYaxis()->SetTitle("Events");                                                     
a_true_mu->GetXaxis()->SetTitle("cos \\theta_{\\mu} ");                                             
a_reco_mu->Draw("HIST");                                                                       
a_true_mu->Draw("HIST SAME");                                                                  
a_sign_mu->Draw("HIST SAME"); 
TLegend *legenda = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legenda->AddEntry(a_true_mu,"true","f");                                                            
legenda->AddEntry(a_reco_mu, "reco", "f");
legenda->AddEntry(a_sign_mu, "signal", "f");
legenda->Draw();  
a1->Update();
a1->SaveAs((cosinetheta+"Muon_Histogram_costheta.root").c_str());                                       
a1->SaveAs((cosinetheta+"Muon_Histogram_costheta.png").c_str());                                        


TCanvas *b1 = new TCanvas();                                                                    
gStyle->SetOptStat(0);    
l_reco_mu->SetLineColorAlpha(kBlue, 0.30);
l_reco_mu->SetLineWidth(2);
l_reco_mu->SetFillColorAlpha(kBlue, 0.20);
l_reco_mu->SetFillStyle(3001);
l_true_mu->SetLineColorAlpha(kRed, 0.30);
l_true_mu->SetLineWidth(2);
l_true_mu->SetFillColorAlpha(kRed, 0.20);
l_true_mu->SetFillStyle(3001);
l_sign_mu->SetLineColorAlpha(kGreen, 0.30);
l_sign_mu->SetLineWidth(2);
l_sign_mu->SetFillColorAlpha(kGreen, 0.20);
l_sign_mu->SetFillStyle(3001);
l_true_mu->Scale(1/(l_true_mu->GetEntries()));                                                 
l_reco_mu->Scale(1/(l_reco_mu->GetEntries()));                                                 
l_sign_mu->Scale(1/(l_sign_mu->GetEntries()));                                                 
l_sign_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_reco_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_true_mu->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_sign_mu->GetYaxis()->SetTitle("Events");                                                     
l_sign_mu->GetXaxis()->SetTitle("length (cm) ");                                             
l_reco_mu->GetYaxis()->SetTitle("Events");                                                     
l_reco_mu->GetXaxis()->SetTitle("length (cm)");                                             
l_true_mu->GetYaxis()->SetTitle("Events");                                                     
l_true_mu->GetXaxis()->SetTitle("length(cm) ");                                             
l_reco_mu->Draw("HIST");                                                                       
l_true_mu->Draw("HIST SAME");                                                                  
l_sign_mu->Draw("HIST SAME"); 
TLegend *legendab1 = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legendab1->AddEntry(l_true_mu,"true","f");                                                            
legendab1->AddEntry(l_reco_mu, "reco", "f");
legendab1->AddEntry(l_sign_mu, "signal", "f");
legendab1->Draw();  
b1->Update();
b1->SaveAs((lengthpath+"Muon_Histogram_length.root").c_str());                                       
b1->SaveAs((lengthpath+"Muon_Histogram_length.png").c_str());                                        





TCanvas *c2 = new TCanvas();                                                                    
gStyle->SetOptStat(0);                                                                     
h_reco_p1->SetLineColorAlpha(kBlue, 0.30);
h_reco_p1->SetLineWidth(2);
h_reco_p1->SetFillColorAlpha(kBlue, 0.20);
h_reco_p1->SetFillStyle(3001);
h_true_p1->SetLineColorAlpha(kRed, 0.30);
h_true_p1->SetLineWidth(2);
h_true_p1->SetFillColorAlpha(kRed, 0.20);
h_true_p1->SetFillStyle(3001);
h_sign_p1->SetLineColorAlpha(kGreen, 0.30);
h_sign_p1->SetLineWidth(2);
h_sign_p1->SetFillColorAlpha(kGreen, 0.20);
h_sign_p1->SetFillStyle(3001);


h_true_p1->Scale(1/(h_true_p1->GetEntries()));                                                 
h_reco_p1->Scale(1/(h_reco_p1->GetEntries()));                                                
h_sign_p1->Scale(1/(h_sign_p1->GetEntries()));                                               
h_sign_p1->GetYaxis()->SetRangeUser(0, 0.25);                                                   
h_reco_p1->GetYaxis()->SetRangeUser(0, 0.25);                                                   
h_true_p1->GetYaxis()->SetRangeUser(0, 0.25);                                                   
h_sign_p1->GetYaxis()->SetTitle("Events");                                                     
h_sign_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_reco_p1->GetYaxis()->SetTitle("Events");                                                     
h_reco_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_true_p1->GetYaxis()->SetTitle("Events");                                                     
h_true_p1->GetXaxis()->SetTitle("p_{p_{1}} (GeV/c)");                                             
h_reco_p1->Draw("HIST");                                                                       
h_true_p1->Draw("HIST SAME");                                                                  
h_sign_p1->Draw("HIST SAME"); 
TLegend *legend2 = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend2->AddEntry(h_true_p1,"true","f");                                                            
legend2->AddEntry(h_reco_p1, "reco", "f");
legend2->AddEntry(h_sign_p1, "signal", "f");
legend2->Draw();  
c2->Update();
c2->SaveAs((momenta+"p1_Histogram_momentum.root").c_str());                                       
c2->SaveAs((momenta+"p1_Histogram_momentum.png").c_str());                                        

TCanvas *a2 = new TCanvas();                                                                    
gStyle->SetOptStat(0);

a_reco_p1->SetLineColorAlpha(kBlue, 0.30);
a_reco_p1->SetLineWidth(2);
a_reco_p1->SetFillColorAlpha(kBlue, 0.20);
a_reco_p1->SetFillStyle(3001);
a_true_p1->SetLineColorAlpha(kRed, 0.30);
a_true_p1->SetLineWidth(2);
a_true_p1->SetFillColorAlpha(kRed, 0.20);
a_true_p1->SetFillStyle(3001);
a_sign_p1->SetLineColorAlpha(kGreen, 0.30);
a_sign_p1->SetLineWidth(2);
a_sign_p1->SetFillColorAlpha(kGreen, 0.20);
a_sign_p1->SetFillStyle(3001);

a_true_p1->Scale(1/(a_true_p1->GetEntries()));                                                 
a_reco_p1->Scale(1/(a_reco_p1->GetEntries()));                                                 
a_sign_p1->Scale(1/(a_sign_p1->GetEntries()));                                                 
a_sign_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_reco_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_true_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
a_sign_p1->GetYaxis()->SetTitle("Events");                                                     
a_sign_p1->GetXaxis()->SetTitle("cos \\theta_{p_{1}} ");                                             
a_reco_p1->GetYaxis()->SetTitle("Events");                                                     
a_reco_p1->GetXaxis()->SetTitle("cos \\theta_{p_{1}} ");                                             
a_true_p1->GetYaxis()->SetTitle("Events");                                                     
a_true_p1->GetXaxis()->SetTitle("cos \\theta_{p_{1}} ");                                             
a_reco_p1->Draw("HIST");                                                                       
a_true_p1->Draw("HIST SAME");                                                                  
a_sign_p1->Draw("HIST SAME"); 
TLegend *legenda1 = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legenda1->AddEntry(a_true_p1,"true","f");                                                            
legenda1->AddEntry(a_reco_p1, "reco", "f");
legenda1->AddEntry(a_sign_p1, "signal", "f");
legenda1->Draw();  
a2->Update();
a2->SaveAs((cosinetheta+"p1_Histogram_costheta.root").c_str());                                       
a2->SaveAs((cosinetheta+"p1_Histogram_costheta.png").c_str());                                        


TCanvas *b2 = new TCanvas();                                                                    
gStyle->SetOptStat(0);    
l_reco_p1->SetLineColorAlpha(kBlue, 0.30);
l_reco_p1->SetLineWidth(2);
l_reco_p1->SetFillColorAlpha(kBlue, 0.20);
l_reco_p1->SetFillStyle(3001);
l_true_p1->SetLineColorAlpha(kRed, 0.30);
l_true_p1->SetLineWidth(2);
l_true_p1->SetFillColorAlpha(kRed, 0.20);
l_true_p1->SetFillStyle(3001);
l_sign_p1->SetLineColorAlpha(kGreen, 0.30);
l_sign_p1->SetLineWidth(2);
l_sign_p1->SetFillColorAlpha(kGreen, 0.20);
l_sign_p1->SetFillStyle(3001);
l_true_p1->Scale(1/(l_true_p1->GetEntries()));                                                 
l_reco_p1->Scale(1/(l_reco_p1->GetEntries()));                                                 
l_sign_p1->Scale(1/(l_sign_p1->GetEntries()));                                                 
l_sign_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_reco_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_true_p1->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_sign_p1->GetYaxis()->SetTitle("Events");                                                     
l_sign_p1->GetXaxis()->SetTitle("length (cm) ");                                             
l_reco_p1->GetYaxis()->SetTitle("Events");                                                     
l_reco_p1->GetXaxis()->SetTitle("length (cm)");                                             
l_true_p1->GetYaxis()->SetTitle("Events");                                                     
l_true_p1->GetXaxis()->SetTitle("length(cm) ");                                             
l_reco_p1->Draw("HIST");                                                                       
l_true_p1->Draw("HIST SAME");                                                                  
l_sign_p1->Draw("HIST SAME"); 
TLegend *legendab2 = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legendab2->AddEntry(l_true_p1,"true","f");                                                            
legendab2->AddEntry(l_reco_p1, "reco", "f");
legendab2->AddEntry(l_sign_p1, "signal", "f");
legendab2->Draw();  
b2->Update();
b2->SaveAs((lengthpath+"p1_Histogram_length.root").c_str());                                       
b2->SaveAs((lengthpath+"p1_Histogram_length.png").c_str());                                        






TCanvas *c3 = new TCanvas();                                                                    
gStyle->SetOptStat(0);                                                                      

h_reco_p2->SetLineColorAlpha(kBlue, 0.30);
h_reco_p2->SetLineWidth(2);
h_reco_p2->SetFillColorAlpha(kBlue, 0.20);
h_reco_p2->SetFillStyle(3001);
h_true_p2->SetLineColorAlpha(kRed, 0.30);
h_true_p2->SetLineWidth(2);
h_true_p2->SetFillColorAlpha(kRed, 0.20);
h_true_p2->SetFillStyle(3001);
h_sign_p2->SetLineColorAlpha(kGreen, 0.30);
h_sign_p2->SetLineWidth(2);
h_sign_p2->SetFillColorAlpha(kGreen, 0.20);
h_sign_p2->SetFillStyle(3001);
h_true_p2->Scale(1/(h_true_p2->GetEntries()));                                                 
h_reco_p2->Scale(1/(h_reco_p2->GetEntries()));                                                 
h_sign_p2->Scale(1/(h_sign_p2->GetEntries()));                                                 
h_sign_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
h_reco_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
h_true_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
h_sign_p2->GetYaxis()->SetTitle("Events");                                                     
h_sign_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_reco_p2->GetYaxis()->SetTitle("Events");                                                     
h_reco_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_true_p2->GetYaxis()->SetTitle("Events");                                                     
h_true_p2->GetXaxis()->SetTitle("p_{p_{2}} (GeV/c)");                                             
h_reco_p2->Draw("HIST");                                                                       
h_true_p2->Draw("HIST SAME");                                                                  
h_sign_p2->Draw("HIST SAME"); 
TLegend *legend3 = new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend3->AddEntry(h_true_p2,"true","f");                                                            
legend3->AddEntry(h_reco_p2, "reco", "f");
legend3->AddEntry(h_sign_p2, "signal", "f");
legend3->Draw();  
c3->Update();
c3->SaveAs((momenta+"p2_Histogram_momentum.root").c_str());                                       
c3->SaveAs((momenta+"p2_Histogram_momentum.png").c_str());                                        

TCanvas *a3 = new TCanvas();                                                                    
gStyle->SetOptStat(0);                                                           

a_reco_p2->SetLineColorAlpha(kBlue, 0.30);
a_reco_p2->SetLineWidth(2);
a_reco_p2->SetFillColorAlpha(kBlue, 0.20);
a_reco_p2->SetFillStyle(3001);
a_true_p2->SetLineColorAlpha(kRed, 0.30);
a_true_p2->SetLineWidth(2);
a_true_p2->SetFillColorAlpha(kRed, 0.20);
a_true_p2->SetFillStyle(3001);
a_sign_p2->SetLineColorAlpha(kGreen, 0.30);
a_sign_p2->SetLineWidth(2);
a_sign_p2->SetFillColorAlpha(kGreen, 0.20);
a_sign_p2->SetFillStyle(3001);

a_true_p2->Scale(1/(a_true_p2->GetEntries()));                                                 
a_reco_p2->Scale(1/(a_reco_p2->GetEntries()));                                                 
a_sign_p2->Scale(1/(a_sign_p2->GetEntries()));                                                 
a_sign_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
a_reco_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
a_true_p2->GetYaxis()->SetRangeUser(0, 0.2);                                                   
a_sign_p2->GetYaxis()->SetTitle("Events");                                                     
a_sign_p2->GetXaxis()->SetTitle("cos \\theta_{p_{2}} ");                                             
a_reco_p2->GetYaxis()->SetTitle("Events");                                                     
a_reco_p2->GetXaxis()->SetTitle("cos \\theta_{p_{2}} ");                                             
a_true_p2->GetYaxis()->SetTitle("Events");                                                     
a_true_p2->GetXaxis()->SetTitle("cos \\theta_{p_{2}} ");                                             
a_reco_p2->Draw("HIST");                                                                       
a_true_p2->Draw("HIST SAME");                                                                  
a_sign_p2->Draw("HIST SAME"); 
TLegend *legenda2 = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legenda2->AddEntry(a_true_p2,"true","f");                                                            
legenda2->AddEntry(a_reco_p2, "reco", "f");
legenda2->AddEntry(a_sign_p2, "signal", "f");
legenda1->Draw();  
a3->Update();
a3->SaveAs((cosinetheta+"p2_Histogram_costheta.root").c_str());                                       
a3->SaveAs((cosinetheta+"p2_Histogram_costheta.png").c_str());                                        

TCanvas *b3 = new TCanvas();                                                                    
gStyle->SetOptStat(0);    
l_reco_p2->SetLineColorAlpha(kBlue, 0.30);
l_reco_p2->SetLineWidth(2);
l_reco_p2->SetFillColorAlpha(kBlue, 0.20);
l_reco_p2->SetFillStyle(3001);
l_true_p2->SetLineColorAlpha(kRed, 0.30);
l_true_p2->SetLineWidth(2);
l_true_p2->SetFillColorAlpha(kRed, 0.20);
l_true_p2->SetFillStyle(3001);
l_sign_p2->SetLineColorAlpha(kGreen, 0.30);
l_sign_p2->SetLineWidth(2);
l_sign_p2->SetFillColorAlpha(kGreen, 0.20);
l_sign_p2->SetFillStyle(3001);
l_true_p2->Scale(1/(l_true_p2->GetEntries()));                                                 
l_reco_p2->Scale(1/(l_reco_p2->GetEntries()));                                                 
l_sign_p2->Scale(1/(l_sign_p2->GetEntries()));                                                 
l_sign_p2->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_reco_p2->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_true_p2->GetYaxis()->SetRangeUser(0, 0.15);                                                   
l_sign_p2->GetYaxis()->SetTitle("Events");                                                     
l_sign_p2->GetXaxis()->SetTitle("length (cm) ");                                             
l_reco_p2->GetYaxis()->SetTitle("Events");                                                     
l_reco_p2->GetXaxis()->SetTitle("length (cm)");                                             
l_true_p2->GetYaxis()->SetTitle("Events");                                                     
l_true_p2->GetXaxis()->SetTitle("length(cm) ");                                             
l_reco_p2->Draw("HIST");                                                                       
l_true_p2->Draw("HIST SAME");                                                                  
l_sign_p2->Draw("HIST SAME"); 
TLegend *legendab3 = new TLegend(0.1, 0.9, 0.3, 0.7);                                             
legendab3->AddEntry(l_true_p2,"true","f");                                                            
legendab3->AddEntry(l_reco_p2, "reco", "f");
legendab3->AddEntry(l_sign_p2, "signal", "f");
legendab3->Draw();  
b3->Update();
b3->SaveAs((lengthpath+"p2_Histogram_length.root").c_str());                                       
b3->SaveAs((lengthpath+"p2_Histogram_length.png").c_str());                                        



file.close();

return 0;
}
