#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <numeric>
#include <time.h>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TAxis.h"

using namespace cppsecrets;
using namespace selection;

int MainTest(const char *config){

  time_t rawtime;
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTime(rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;


  //------------------------------------------------------------------------------------------
  //                                    Configure
  //------------------------------------------------------------------------------------------
  // Create object of the class ConfigReader
  // Parse the configuration file
  // Dump map on the console after parsing it
  ConfigReader* p = ConfigReader::getInstance();
  p->parseFile(config);
  std::cout << " Variables from configuration file: " << std::endl;
  p->dumpFileValues();
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Get variables from config
  std::string input_location  = "";
  std::string input_filename  = "";
  std::string exceptions_file = "";
  std::string file_location   = "";
  unsigned int total_files = 0;
  unsigned int detector = 0; // 0 = sbnd, 1 = uboone, 2 = icarus
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputFileLocation",input_location);
  p->getValue("InputFileName",    input_filename);
  p->getValue("ExceptionsFile",   exceptions_file);
  p->getValue("TreeFileLocation", file_location); 
  p->getValue("Detector",         detector);
  p->getValue("TotalFiles",       total_files);
  p->getValue("MinXFid",          minx_fid);
  p->getValue("MinYFid",          miny_fid);
  p->getValue("MinZFid",          minz_fid);
  p->getValue("MaxXFid",          maxx_fid);
  p->getValue("MaxYFid",          maxy_fid);
  p->getValue("MaxZFid",          maxz_fid);
  p->getValue("MinXAV",           minx_av);
  p->getValue("MinYAV",           miny_av);
  p->getValue("MinZAV",           minz_av);
  p->getValue("MaxXAV",           maxx_av);
  p->getValue("MaxYAV",           maxy_av);
  p->getValue("MaxZAV",           maxz_av);

  //------------------------------------------------------------------------------------------
  //                                    Initialise
  //------------------------------------------------------------------------------------------

  // Get the active and fiducial geometry objects
  Geometry fiducial(minx_fid,miny_fid,minz_fid,maxx_fid,maxy_fid,maxz_fid,true);
  Geometry active(minx_av,miny_av,minz_av,maxx_av,maxy_av,maxz_av,false);
  PlaneList planes = active.GetExternalPlaneList();

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;

  int start = static_cast<int>(time(NULL));
  double pot = 0.; 

  std::vector<int> exceptions;
  FillExceptions(exceptions_file.c_str(),exceptions);

  // Maps
  TopologyMap ccinc = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc2pi = GeneralAnalysisHelper::GetCC2PiTopologyMap();
  TopologyMap ccpi0 = GeneralAnalysisHelper::GetCCPi0TopologyMap();

  // Variables to get
  bool iscc, isnc;  // Charged or neutral current event current
  int nu_pdg, mode; // Neutrino pdg code and scattering mode
  double baseline; // Baseline of the neutrino from the flux
  double enu_true, enu_reco, enu_reco_mc, qsqr; // Neutrino kinematics, truth and reco
  double mu_momentum, mu_cos_z; // Muon kinematics, reco
  int nkaons, npip, npim, npi0; // Particle counting, truth
  // Enumeration (all CC):
  //   -1 = undefined
  //    0 = 0pi
  //    1 = 1chpi
  //    2 = 2chpi
  //    3 = 1pi0
  //    4 = other
  int true_topology = -1;
  int reco_topology = -1;

  TTree *t_run    = new TTree("valor_tree"," Tree to hold variables needed for the VALOR analysis");
  TTree *t_subrun = new TTree("subrun_tree"," Tree to hold subrun variables, such as POT");
  t_run->SetDirectory(0);
  t_subrun->SetDirectory(0);

  t_run->Branch("iscc",          &iscc,          "iscc/O");
  t_run->Branch("isnc",          &isnc,          "isnc/O");
  t_run->Branch("baseline",      &baseline,      "baseline/D");
  t_run->Branch("nu_pdg",        &nu_pdg,        "nu_pdg/I");
  t_run->Branch("enu_true",      &enu_true,      "enu_true/D");
  t_run->Branch("enu_reco_mc",   &enu_reco_mc,   "enu_reco_mc/D");
  t_run->Branch("enu_reco",      &enu_reco,      "enu_reco/D");
  t_run->Branch("mu_momentum",   &mu_momentum,   "mu_momentum/D");
  t_run->Branch("mu_cos_z",      &mu_cos_z,      "mu_cos_z/D");
  t_run->Branch("qsqr",          &qsqr,          "qsqr/D");
  t_run->Branch("mode",          &mode,          "mode/I");
  t_run->Branch("nkaons",        &nkaons,        "nkaons/I");
  t_run->Branch("npip",          &npip,          "npip/I");
  t_run->Branch("npim",          &npim,          "npim/I");
  t_run->Branch("npi0",          &npi0,          "npi0/I");
  t_run->Branch("true_topology", &true_topology, "true_topology/I");
  t_run->Branch("reco_topology", &reco_topology, "reco_topology/I");

  t_subrun->Branch("pot",        &pot,           "pot/D");

  // Counter to quantify how many events have strange particle energies in them
  int bad_events = 0;

  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    for(const Event &e : events){
      bool cc_inclusive_passed = GeneralAnalysisHelper::PassedCCInclusive(e,detector);

      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial() || 
         !GeneralAnalysisHelper::MaxOneLongEscapingTrack(e) || 
         !GeneralAnalysisHelper::MinOneRecoTrack(e)) continue;

      // Now we are looking at selected CC Inclusive events, 
      // start trying to identify particle types
      if(cc_inclusive_passed){
        iscc     = e.GetIsCC();
        isnc     = !e.GetIsCC();
        nu_pdg   = e.GetNeutrinoPdgCode();
        enu_true = e.GetTrueNuEnergy();
        qsqr     = e.GetTrueNuQ2();
        mode     = e.GetPhysicalProcess();
        baseline = e.GetBaseline();

        // Start the counters
        nkaons      = 0;
        npip        = 0;
        npim        = 0;
        npi0        = 0;
        enu_reco    = 0.;
        enu_reco_mc = 0.;

        // Particles 
        ParticleList mc   = e.GetMCParticleList();
        ParticleList reco = e.GetRecoParticleList();

        if(e.CheckRecoTopology(ccinc)){
          for(const Particle &p : reco){
            double mass = 0.;
            double KE   = 0.;
            KE = sqrt(pow(p.GetModulusMomentum(),2) + pow(p.GetMass(),2)) - p.GetMass();
            if(KE <= 0. || isnan(KE)) continue;
            if(KE > 10) { // Higher than 10 GeV
              bad_events++;
              break;
            }
            if(p.GetPdgCode() == 13){
              mu_momentum = p.GetModulusMomentum();
              mu_cos_z    = p.GetCosTheta();
            } 
            if(p.GetFromRecoTrack()) {
              // Calculated KE from momentum and mass for every particle

              if(p.GetPdgCode() == 211  || 
                  p.GetPdgCode() == -211 ||
                  p.GetPdgCode() == 13){
                mass = p.GetMass();
              }
            }
            // Calculate the reconstructed energy from the reco kinetic (visible) + mass of outgoing particles 
            // Do not include hadron mass since the interaction didn't *produce* a hadron, it was *with* a hadron
            enu_reco += KE + mass;
            if(isnan(enu_reco)){
              std::cout << " Found NAN: " << enu_reco << std::endl;
              enu_reco = 0;
            }
          }
          for(const Particle &p : mc){
            if(p.GetPdgCode() ==  311 || p.GetPdgCode() == -321 || p.GetPdgCode() == 321) nkaons++;
            if(p.GetPdgCode() ==  211) npip++;
            if(p.GetPdgCode() == -211) npim++;
            if(p.GetPdgCode() ==  111) npi0++; 

            // MC reconstructed neutrino energy 
            double mass = 0.;
            double KE   = 0.;
            // Calculated KE from momentum and mass for every particle
            KE = sqrt(pow(p.GetModulusMomentum(),2) + pow(p.GetMass(),2)) - p.GetMass();

            if(p.GetPdgCode() == 211  || 
                p.GetPdgCode() == -211 ||
                p.GetPdgCode() == 13){
              mass = p.GetMass();
            } // Pdgcodes
            enu_reco_mc += KE + mass;
            //enu_reco_mc += (p.GetKineticEnergy() + mass);
          } // Particles

          if(e.CheckRecoTopology(cc0pi))
            reco_topology = 0;
          else if(e.CheckRecoTopology(cc1pi))
            reco_topology = 1;
          else if(e.CheckRecoTopology(cc2pi))
            reco_topology = 2;
          else if(e.CheckRecoTopology(ccpi0))
            reco_topology = 3;
          else
            reco_topology = 4;
          if(e.CheckMCTopology(cc0pi))
            true_topology = 0;
          else if(e.CheckMCTopology(cc1pi))
            true_topology = 1;
          else if(e.CheckMCTopology(cc2pi))
            true_topology = 2;
          else if(e.CheckMCTopology(ccpi0))
            true_topology = 3;
          else
            true_topology = 4;

          t_run->Fill();
        } // Chosen Topology
      } // CCInc Precuts
    } // Events
  } // Files
  std::cout << " Total events with particles with energy > 10GeV (skipped) : " << bad_events << std::endl;
  // Print the total pot from all the samples
  std::cout << " Total POT in the samples is: " << pot << std::endl;
  t_subrun->Fill();

  // Output TFile
  TFile f((file_location+"ccinc_selection_final.root").c_str(), "RECREATE");

  t_run->Write();
  t_subrun->Write();

  f.Write();
  f.Close();

  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest

