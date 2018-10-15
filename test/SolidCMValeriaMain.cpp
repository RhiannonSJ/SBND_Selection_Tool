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
#include "TGraph.h"
#include "TCanvas.h"                                                                           
 #include "TLegend.h"                                                                           
 #include "TLatex.h"                                                                            
 #include "TStyle.h"                                                                            
 #include "TColor.h"                                                                            
 #include "TObjArray.h"                                                                         
 #include "THStack.h"
 #include "TProfile.h"
 #include "TGraph.h"
 #include <math.h>
#include <vector>
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
   std::string angle = "../Output_Selection_Tool/plots/solid/angle_distributions/CM_angle/";
   std::string stats_location= "../Output_Selection_Tool/statistics/";
                                                                                                
   //Load events                                         
   
   // Initialise event list and the topology maps                                               
   EventSelectionTool::EventList events;
   int start = static_cast<int>(time(NULL));                                                    
  
   unsigned int total = 500 ;                                                                    
                                                                                               
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
int cc0pi2pmc    = 0; // counter of cc0pi2p true events
int cc0pi2pre    = 0; // counter of cc0pi2p reconstructed events  
int cc0pi2ps2    = 0; // counter of cc0pi2p signal events by protons
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
Double_t momento; 
Double_t momento2;
Double_t momento3; 
std:: vector <float> mo;
std:: vector <Double_t> res, res2;
std:: vector <Double_t> zezi;
std:: vector <Double_t> mores;
TVector3 neutrino, momentum, n, q, qs, nu_sign, momentum_sign, momentum_p, n_sign;
TVector3 p1, p2, p1r, p2r, p1s, p2s;
std:: vector <Double_t> angle_true;
std:: vector <Double_t> angle_reco;
std:: vector <Double_t> angle_sign;
std:: vector <int> interactions; 
interactions.push_back(1); // QE 
interactions.push_back(3); // DIS 
interactions.push_back(4); // RESONANCE
interactions.push_back(5); // COHERENT 
interactions.push_back(10);// MEC 
std:: vector < vector <Double_t> >  angle_int_sign;
std:: vector <Double_t> angle_QE_sign, angle_DIS_sign, angle_RES_sign, angle_COH_sign, angle_MEC_sign;
angle_int_sign.push_back(angle_QE_sign);
angle_int_sign.push_back(angle_DIS_sign);
angle_int_sign.push_back(angle_RES_sign);
angle_int_sign.push_back(angle_COH_sign);
angle_int_sign.push_back(angle_MEC_sign);
std:: vector < vector <Double_t> >  angle_int_reco;
std:: vector <Double_t> angle_QE_reco, angle_DIS_reco, angle_RES_reco, angle_COH_reco, angle_MEC_reco;
angle_int_reco.push_back(angle_QE_reco);
angle_int_reco.push_back(angle_DIS_reco);
angle_int_reco.push_back(angle_RES_reco);
angle_int_reco.push_back(angle_COH_reco);
angle_int_reco.push_back(angle_MEC_reco); 
std:: vector < vector <Double_t> >  angle_int_true;
std:: vector <Double_t> angle_QE_true, angle_DIS_true, angle_RES_true, angle_COH_true, angle_MEC_true;
angle_int_true.push_back(angle_QE_true);
angle_int_true.push_back(angle_DIS_true);
angle_int_true.push_back(angle_RES_true);
angle_int_true.push_back(angle_COH_true);
angle_int_true.push_back(angle_MEC_true);

Double_t angletrue, anglereco, anglesign;
ofstream file;
file.open(stats_location+"solidCM.txt");

// loop over the events 
for(const Event &e : events){
  if(e.IsSBNDTrueFiducial()){
    if(bool max = GeneralAnalysisHelper::MaxOneEscapingTrack(e)){
    ++ntotevents;
    ParticleList true_particles = e.GetMCParticleList();  
    ParticleList reco_particles = e.GetRecoParticleList();                                    

    // cc0pi2P true events
    
    if(e.CheckMCTopology(cc0pi2p_signal_map)){
      ++cc0pi2pmc;  
      //loop over the true particles

      if(unsigned int j=e.CountMCParticlesWithPdg(2212)==2){  
        for(Particle &p_true : true_particles){
          float nu_energy = CC0piAnalysisHelper::GetMCCC0piNeutrinoEnergy(e);
          neutrino(0)=0;
          neutrino(1)=0;
          neutrino(2)=nu_energy; 
          if(p_true.GetPdgCode()==13){
            ++mucount;
            momentum=p_true.GetMomentum();
          }
          else if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021 && p_true.GetModulusMomentum()>=0.250){
            ++nproton; 
            ++ptottrue; 
            momentum=p_true.GetMomentum();
            mo.push_back(p_true.GetModulusMomentum());       
            for(int i=0; i<3; i++){
              res.push_back(momentum(i));
            }
          }
        }
      }
      if(mo[0]<mo[1] && mo[1]!=0 && mo[0]!=0){
        for(unsigned int i=0; i<3; i++){
          p1(i)=res[i];
        }
        for(unsigned int i=0; i<3; i++){
          p2(i)=res[i+3];
        }
      }
      else if(mo[1]<mo[0] && mo[1]!=0 && mo[0]!= 0){
        for(unsigned int i=0; i<3; i++){
          p2(i)=res[i];
        } 
        for(unsigned int i=0; i<3; i++){
          p1(i)=res[i+3];
        }
      }
      for(unsigned int i=0; i<3; i++){
        q(i)=neutrino(i)-momentum(i);
      }
      for(unsigned int i=0; i<3; i++){
        n(i)=p2(i)-q(i);
      }
      angle_true.push_back(n.Angle(p1));
      angletrue=n.Angle(p1);
        for(unsigned int i=0; i<interactions.size(); i++){
          if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
             angle_int_true[i].push_back(angletrue); 
          } 
        }
     
      
      
      mo.clear(); 
      res.clear();
      nproton=0;
      mucount=0;
    }
    // CC0Pi2P reconstructed events

    if(e.CheckRecoTopology(cc0pi2p_signal_map)){
      ++cc0pi2pre;
      if(unsigned int k=e.CountRecoParticlesWithPdg(2212)==2){
        // loop over the reconstructed particles 
        
        for(Particle &p_reco : reco_particles){       
          float nu_energy = CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e);
          neutrino(0)=0;
          neutrino(1)=0;
          neutrino(2)=nu_energy;
          for(Particle &p_true:true_particles){
            int mcid = p_true.GetMCId();
            int phits = p_reco.GetMCParticleIdHits();                 
            if(mcid==phits){            
              float nu_energy_s = CC0piAnalysisHelper::GetMCCC0piNeutrinoEnergy(e);
              nu_sign(0)=0;
              nu_sign(1)=0;
              nu_sign(2)=nu_energy_s;    
            }       
          }          
          if(p_reco.GetFromRecoTrack()){
            if(p_reco.GetPdgCode()==13){         
              ++mucountr;           
              momentum=p_reco.GetMomentum();

           //get the MC particle by hits 
          
              int phits = p_reco.GetMCParticleIdHits();           
              for(Particle &p_true : true_particles){
                int mcid = p_true.GetMCId();
                if(mcid==phits){            
                  if(p_true.GetPdgCode()==13){ 
                    ++mucounts;
                    momentum_sign=p_reco.GetMomentum();
                    
                  }
                }
              }
            }
            else if(p_reco.GetPdgCode()==2212 && p_reco.GetKineticEnergy()>=0.021 && p_reco.GetModulusMomentum()>=0.250){
              ++ptotreco;
              ++nprotonreco; 
              mo.push_back(p_reco.GetModulusMomentum());
              momentum=p_reco.GetMomentum();
              for(int i=0; i<3; i++){
                res.push_back(momentum(i));
              }
              //get the MC particle by hits
              
              int  phits = p_reco.GetMCParticleIdHits();           
              for(Particle &p_true : true_particles){
                int mcid = p_true.GetMCId(); 
                    if(p_true.GetPdgCode()==2212 && p_true.GetKineticEnergy()>=0.021 && p_true.GetModulusMomentum()>=0.250){ 
                      if(mcid==phits){
                        ++nprotonsig;
                        ++ptotsig;
                        zezi.push_back(p_reco.GetModulusMomentum());         
                        momentum_p=p_reco.GetMomentum();
                        for(int i=0; i<3; i++){
                          res2.push_back(momentum_p(i));
                        }
                      } 
                    }              
              }
            }
          }
        } //closes the loop on reco particles
        if(mo[0]<mo[1]){
          for(unsigned int i=0; i<3; i++){
            p1r(i)=res[i];
          }
          for(unsigned int i=0; i<3; i++){
            p2r(i)=res[i+3];
          }
        }
        else if(mo[1]<=mo[0]){
          for(unsigned int i=0; i<3; i++){
            p2r(i)=res[i];     
          }
          for(unsigned int i=0; i<3; i++){
            p1r(i)=res[i+3];
          }  
        }
        for(unsigned int i=0; i<3; i++){
          q(i)=neutrino(i)-momentum(i);
        }
        for(unsigned int i=0; i<3; i++){
          n(i)=p2r(i)-q(i);
        }      
        angle_reco.push_back(n.Angle(p1r));
        anglereco=n.Angle(p1r);
        for(unsigned int i=0; i<interactions.size(); i++){
          if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
             angle_int_reco[i].push_back(anglereco); 
          } 
        }
        
        if(nprotonsig==2 && mucounts==1){
          ++cc0pi2ps2;
          if(zezi[0]<=zezi[1]){
            for(unsigned int i=0; i<3; i++){
              p1s(i)=res2[i];
            }
            for(unsigned int i=0; i<3; i++){
              p2s(i)=res2[i+3];
            }
            
          }
          else if(zezi[1]<zezi[0]){
            for(unsigned int i=0; i<3; i++){
              p2s(i)=res2[i];
            }
            for(unsigned int i=0; i<3; i++){
              p1s(i)=res2[i+3];
            }
          }
          for(unsigned int i=0; i<3; i++){
            qs(i)=nu_sign(i)-momentum_sign(i);
          }
          for(unsigned int i=0; i<3; i++){
            n_sign(i)=p2s(i)-qs(i);
          }
          angle_sign.push_back(n_sign.Angle(p1s));
          anglesign=n_sign.Angle(p1s);
       for(unsigned int i=0; i<interactions.size(); i++){
          if (unsigned int j=e.GetPhysicalProcess()==interactions[i]){
             angle_int_sign[i].push_back(anglesign); 
          } 
        }
        }
       
        res.clear(); 
        nprotonreco=0;
        mucountr=0;
        zezi.clear();
        nprotonsig=0;
        mucounts=0;
        res2.clear();

      }
    }  
    } // end of the condition on reconstructed topology map

  }

} // end of the loop on the event         
/*for(unsigned int k=0; k< angle_sign.size(); k++){
         cout << "angolo = " << angle_sign[k] << endl;
}
for(unsigned int i=0; i<interactions.size(); i++){

//for(unsigned int j=0; j<angle_int_true[i].size(); j++){
  cout <<  "dimensioni angolo in funzione interazione = " << i << " " << angle_int_true[i].size() << endl;
 //}         
}
  cout << "dimensioni angolo tot = " << angle_true.size() << endl; 
 

for(unsigned int i=0; i<angle_true.size(); i++){
 cout << "angle true = " << angle_true[i] << " "; 
} 
cout << endl;


for(unsigned int i=0; i<angle_reco.size(); i++){
 cout << "angle reco = " << angle_reco[i] << " "; 
} 
cout << endl;
 
for(unsigned int i=0; i<angle_sign.size(); i++){
 cout << "angle sign = " << angle_sign[i] << " "; 
} 
cout << endl;
 */
// ============================================
//               Fill histograms
// ============================================


// ____________ cosine of solid angle distributions _______________


TH1D *a_reco = new TH1D("a_reco", " angle distribution", 40, -1,1);                                        
TH1D *a_true = new TH1D("a_true", " angle distribution" , 40, -1, 1);                                    
TH1D *a_sign = new TH1D("a_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_true.size(); i++){                                       
  a_true->Fill(cos(angle_true[i]));                                                      
}

for(unsigned int i=0; i<angle_reco.size(); i++){                                       
  a_reco->Fill(cos(angle_reco[i]));                                                      
}

for(unsigned int i=0; i<angle_sign.size(); i++){                                       
  a_sign->Fill(cos(angle_sign[i]));                                                      
}

//QE
TH1D *qe_reco = new TH1D("qe_reco", " angle distribution", 40, -1,1);                                        
TH1D *qe_true = new TH1D("qe_true", " angle distribution" , 40, -1, 1);                                    
TH1D *qe_sign = new TH1D("qe_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_int_true[0].size(); i++){                                       
  qe_true->Fill(cos(angle_int_true[0][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[0].size(); i++){                                       
  qe_reco->Fill(cos(angle_int_reco[0][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[0].size(); i++){                                       
  qe_sign->Fill(cos(angle_int_sign[0][i]));                                                      
}




// DIS
TH1D *dis_reco = new TH1D("dis_reco", " angle distribution", 40, -1,1);                                        
TH1D *dis_true = new TH1D("dis_true", " angle distribution" , 40, -1, 1);                                    
TH1D *dis_sign = new TH1D("dis_sign", " angle distribution", 40, -1, 1);                                       

for(unsigned int i=0; i<angle_int_true[1].size(); i++){                                       
  dis_true->Fill(cos(angle_int_true[1][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[1].size(); i++){                                       
  dis_reco->Fill(cos(angle_int_reco[1][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[1].size(); i++){                                       
  dis_sign->Fill(cos(angle_int_sign[1][i]));                                                      
}

// RES
TH1D *res_reco = new TH1D("res_reco", " angle distribution", 40, -1,1);                                        
TH1D *res_true = new TH1D("res_true", " angle distribution" , 40, -1, 1);                                    
TH1D *res_sign = new TH1D("res_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[2].size(); i++){                                       
  res_true->Fill(cos(angle_int_true[2][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[2].size(); i++){                                       
  res_reco->Fill(cos(angle_int_reco[2][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[2].size(); i++){                                       
 res_sign->Fill(cos(angle_int_sign[2][i]));                                                      
}


// COH
TH1D *coh_reco = new TH1D("coh_reco", " angle distribution", 40, -1,1);                                        
TH1D *coh_true = new TH1D("coh_true", " angle distribution" , 40, -1, 1);                                    
TH1D *coh_sign = new TH1D("coh_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[3].size(); i++){                                       
  coh_true->Fill(cos(angle_int_true[3][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[3].size(); i++){                                       
  coh_reco->Fill(cos(angle_int_reco[3][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[3].size(); i++){                                       
  coh_sign->Fill(cos(angle_int_sign[3][i]));                                                      
}


// MEC

TH1D *mec_reco = new TH1D("mec_reco", " angle distribution", 40, -1,1);                                        
TH1D *mec_true = new TH1D("mec_true", " angle distribution" , 40, -1, 1);                                    
TH1D *mec_sign = new TH1D("mec_sign", " angle distribution", 40, -1, 1);                                       


for(unsigned int i=0; i<angle_int_true[4].size(); i++){                                       
  mec_true->Fill(cos(angle_int_true[4][i]));                                                      
}

for(unsigned int i=0; i<angle_int_reco[4].size(); i++){                                       
  mec_reco->Fill(cos(angle_int_reco[4][i]));                                                      
}

for(unsigned int i=0; i<angle_int_sign[4].size(); i++){                                       
  mec_sign->Fill(cos(angle_int_sign[4][i]));                                                      
}




file << "=====================================================================" << endl;
file << "enties for the histogram of true mu " << a_true->GetEntries() << endl;
file << "entries for the histogram of reco mu " << a_reco->GetEntries() << endl;
file << "entries for the histogram of signal mu " << a_sign->GetEntries() << endl;
file << "=====================================================================" << endl;


cout << 0 << endl;

//=========================================================
//                        Drawing
//=========================================================

// cosine of solid angle in CM frame distribution 

TCanvas *a = new TCanvas();                                                                    
gStyle->SetOptStat(0);
a_true->SetLineColorAlpha(kBlue, 0.30);
a_true->SetLineWidth(2);
a_true->SetFillColorAlpha(kBlue, 0.20);
a_true->SetFillStyle(3001);
a_reco->SetLineColorAlpha(kRed, 0.30);
a_reco->SetLineWidth(2);
a_reco->SetFillColorAlpha(kRed, 0.20);
a_reco->SetFillStyle(3001);
a_sign->SetLineColorAlpha(kGreen, 0.30);
a_sign->SetLineWidth(2);
a_sign->SetFillColorAlpha(kGreen, 0.20);
a_sign->SetFillStyle(3001);
a_true->Scale(1/(a_true->GetEntries()));                                                 
a_reco->Scale(1/(a_reco->GetEntries()));                                                 
a_sign->Scale(1/(a_sign->GetEntries()));                                                 
a_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
a_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
a_true->GetYaxis()->SetRangeUser(0, 0.18);                                                  
a_sign->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
a_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
a_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
a_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
a_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
a_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
a_reco->Draw("HIST");                                                                       
a_true->Draw("HIST SAME");                                                                  
a_sign->Draw("HIST SAME"); 
TLegend *legendn= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legendn->AddEntry(a_true,"true","f");                                                            
legendn->AddEntry(a_reco, "reco", "f");
legendn->AddEntry(a_sign, "signal", "f");
legendn->Draw();  
a->Update();
a->SaveAs((angle+"Angle_cm_Histogram.root").c_str());                                       
a->SaveAs((angle+"Angle_cm_Histogram.png").c_str());                                        

// true events 
//
TCanvas *b = new TCanvas();                                                                    
gStyle->SetOptStat(0); 
qe_true->SetFillColorAlpha(4, 0.5);
qe_true->Scale(1/(a_true->GetEntries()));                                                 
qe_true->GetYaxis()->SetRangeUser(0, 0.18);                                                   
qe_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
qe_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

dis_true->SetFillColorAlpha(2, 0.5);
dis_true->Scale(1/(a_true->GetEntries()));                                                 
dis_true->GetYaxis()->SetRangeUser(0, 0.18);                                                   
dis_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
dis_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

res_true->SetFillColorAlpha(5, 0.5);
res_true->Scale(1/(a_true->GetEntries()));                                                 
res_true->GetYaxis()->SetRangeUser(0, 0.18);                                                   
res_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
res_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

coh_true->SetFillColorAlpha(6, 0.5);
coh_true->Scale(1/(a_true->GetEntries()));                                                 
coh_true->GetYaxis()->SetRangeUser(0, 0.18);                                                   
coh_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
coh_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

mec_true->SetFillColorAlpha(3, 0.5);
mec_true->Scale(1/(a_true->GetEntries()));                                                 
mec_true->GetYaxis()->SetRangeUser(0, 0.18);                                                   
mec_true->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
mec_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             


THStack *hs_true=new THStack("hs_true", "cos(\\theta^{,}_{p1n}) true" ); 
hs_true->Add(qe_true);
hs_true->Add(res_true);
hs_true->Add(coh_true);
hs_true->Add(dis_true);
hs_true->Add(mec_true);
hs_true->Write();
hs_true->Draw("HIST");
hs_true->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");
hs_true->GetYaxis()->SetTitle(" Arbitrary  units");

TLegend *legend_true= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_true->AddEntry(qe_true,"QE","f");                                                            
legend_true->AddEntry(res_true, "RES", "f");
legend_true->AddEntry(coh_true, "COH", "f");
legend_true->AddEntry(dis_true, "DIS" , "f");
legend_true->AddEntry(mec_true, "MEC", "f");
legend_true->Draw("SAME");  
b->Update();
b->SaveAs((angle+"Angle_cm_interactions_true.root").c_str());                                       
b->SaveAs((angle+"Angle_cm_interactions_true.png").c_str());                                        

// reco events 

TCanvas *c = new TCanvas();
gStyle->SetOptStat(0);
qe_reco->SetLineWidth(2);
qe_reco->SetFillColorAlpha(4, 0.5);
qe_reco->Scale(1/(a_reco->GetEntries()));                                                 
qe_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
qe_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
qe_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

dis_reco->SetFillColorAlpha(2, 0.5);
dis_reco->Scale(1/(a_reco->GetEntries()));                                                 
dis_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
dis_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
dis_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

res_reco->SetFillColorAlpha(5, 0.5);
res_reco->Scale(1/(a_reco->GetEntries()));                                                 
res_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
res_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
res_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

coh_reco->SetFillColorAlpha(6, 0.5);
coh_reco->Scale(1/(a_reco->GetEntries()));                                                 
coh_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
coh_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
coh_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             

mec_reco->SetFillColorAlpha(3, 0.5);
mec_reco->Scale(1/(a_reco->GetEntries()));                                                 
mec_reco->GetYaxis()->SetRangeUser(0, 0.18);                                                   
mec_reco->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
mec_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             




THStack *hs_reco=new THStack("hs_reco","cos(\\theta^{,}_{p1n}) reco" );
hs_reco->Add(qe_reco);
hs_reco->Add(res_reco);
hs_reco->Add(coh_reco);
hs_reco->Add(dis_reco);
hs_reco->Add(mec_reco);
hs_reco->Draw("HIST"); 
hs_reco->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");
hs_reco->GetYaxis()->SetTitle(" Arbitrary  units");


TLegend *legend_reco= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_reco->AddEntry(qe_reco,"QE","f");                                                            
legend_reco->AddEntry(res_reco, "RES", "f");
legend_reco->AddEntry(coh_reco, "COH", "f");
legend_reco->AddEntry(dis_reco, "DIS" , "f");
legend_reco->AddEntry(mec_reco, "MEC", "f");
legend_reco->Draw("SAME");  
c->Update();
c->SaveAs((angle+"Angle_cm_interactions_reco.root").c_str());                                       
c->SaveAs((angle+"Angle_cm_interactions_reco.png").c_str());                                        


// reco events 

TCanvas *d = new TCanvas();
gStyle->SetOptStat(0);
qe_sign->SetLineWidth(2);
qe_sign->SetFillColorAlpha(4, 0.5);
qe_sign->Scale(1/(a_sign->GetEntries()));                                                 
qe_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
qe_sign->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
qe_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
qe_sign->Draw("HIST");

dis_sign->SetFillColorAlpha(2, 0.5);
dis_sign->Scale(1/(a_sign->GetEntries()));                                                 
dis_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
dis_sign->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
dis_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
dis_sign->Draw("HIST SAME");

res_sign->SetFillColorAlpha(5, 0.5);
res_sign->Scale(1/(a_sign->GetEntries()));                                                 
res_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
res_sign->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
res_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
res_sign->Draw("HIST SAME");

coh_sign->SetFillColorAlpha(6, 0.5);
coh_sign->Scale(1/(a_sign->GetEntries()));                                                 
coh_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
coh_sign->GetYaxis()->SetTitle(" Arbitrary  units");                                                     
coh_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             



mec_sign->SetFillColorAlpha(3, 0.5);
mec_sign->Scale(1/(a_sign->GetEntries()));                                                 
mec_sign->GetYaxis()->SetRangeUser(0, 0.18);                                                   
mec_sign->GetYaxis()->SetTitle("Arbitrary  units");                                                     
mec_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");                                             
mec_sign->Draw("HIST SAME");



THStack *hs_sign=new THStack("hs_sign","cos(\\theta^{,}_{p1n}) signal");
hs_sign->Add(qe_sign);
hs_sign->Add(res_sign);
hs_sign->Add(coh_sign);
hs_sign->Add(dis_sign);
hs_sign->Add(mec_sign);
hs_sign->Draw("HIST");
hs_sign->GetXaxis()->SetTitle("cos(\\theta^{,}_{p1n})");
hs_sign->GetYaxis()->SetTitle(" Arbitrary  units");

TLegend *legend_sign= new TLegend(0.9, 0.7, 0.7, 0.9);                                             
legend_sign->AddEntry(qe_sign,"QE","f");                                                            
legend_sign->AddEntry(res_sign, "RES", "f");
legend_sign->AddEntry(coh_sign, "COH", "f");
legend_sign->AddEntry(dis_sign, "DIS" , "f");
legend_sign->AddEntry(mec_sign, "MEC", "f");
legend_sign->Draw("SAME");  
d->Update();
d->SaveAs((angle+"Angle_cm_interactions_sign.root").c_str());                                       
d->SaveAs((angle+"Angle_cm_interactions_sign.png").c_str());                                        



file.close();

return 0;
}
