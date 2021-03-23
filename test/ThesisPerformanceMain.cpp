#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Geometry.h"
#include "../include/Plane.h"
#include "../include/Event.h"
#include "../include/Particle.h"
#include "../include/Setup.h"
#include "../include/ConfigReader.h"

using namespace cppsecrets;
using namespace selection;

void SetHistErrors(TH1D *h, TH1D *den);

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
  std::string stats_location  = "";
  std::string plots_location  = "";
  unsigned int total_files = 0;
  unsigned int detector = 0; // 0 = sbnd, 1 = uboone, 2 = icarus
  double totPOT = 0.;
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputFileLocation",input_location);
  p->getValue("InputFileName",    input_filename);
  p->getValue("ExceptionsFile",   exceptions_file);
  p->getValue("StatFileLocation", stats_location); 
  p->getValue("PlotFileLocation", plots_location); 
  p->getValue("Detector",         detector);
  p->getValue("TotalFiles",       total_files);
  p->getValue("TotalPOT",         totPOT);
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

  // HISTOGRAMS

  // True CC Inclusive 
  TH1D *h_true_ccinc_enu = new TH1D("h_true_ccinc_enu", "",40,0,3);
  TH1D *h_true_ccinc_pmu = new TH1D("h_true_ccinc_pmu", "",40,0,3);
  TH1D *h_true_ccinc_pth = new TH1D("h_true_ccinc_pth", "",40,-1,1);

  // True CC Inclusive 
  TH1D *h_true_cc0pi_enu = new TH1D("h_true_cc0pi_enu", "",40,0,3);
  TH1D *h_true_cc0pi_pmu = new TH1D("h_true_cc0pi_pmu", "",40,0,3);
  TH1D *h_true_cc0pi_pth = new TH1D("h_true_cc0pi_pth", "",40,-1,1);

  // Signal CC Inclusive 
  TH1D *h_signal_ccinc_enu = new TH1D("h_signal_ccinc_enu", "",40,0,3);
  TH1D *h_signal_ccinc_pmu = new TH1D("h_signal_ccinc_pmu", "",40,0,3);
  TH1D *h_signal_ccinc_pth = new TH1D("h_signal_ccinc_pth", "",40,-1,1);

  // Signal CC Inclusive 
  TH1D *h_signal_cc0pi_enu = new TH1D("h_signal_cc0pi_enu", "",40,0,3);
  TH1D *h_signal_cc0pi_pmu = new TH1D("h_signal_cc0pi_pmu", "",40,0,3);
  TH1D *h_signal_cc0pi_pth = new TH1D("h_signal_cc0pi_pth", "",40,-1,1);

  // COUNTERS
  unsigned int nue_true  = 0; 

  unsigned int ccinc_cc0pi = 0;
  unsigned int ccinc_cc1pi = 0;
  unsigned int ccinc_ccoth = 0;
  unsigned int ccinc_nc0pi = 0;
  unsigned int ccinc_nc1pi = 0;
  unsigned int ccinc_ncoth = 0;
  unsigned int ccinc_nue   = 0;
  unsigned int ccinc_true  = 0; 
  unsigned int ccinc_sig   = 0; 
  unsigned int ccinc_sel   = 0; 

  unsigned int cc0pi_cc0pi = 0;
  unsigned int cc0pi_cc1pi = 0;
  unsigned int cc0pi_ccoth = 0;
  unsigned int cc0pi_nc0pi = 0;
  unsigned int cc0pi_nc1pi = 0;
  unsigned int cc0pi_ncoth = 0;
  unsigned int cc0pi_nue   = 0;
  unsigned int cc0pi_true  = 0; 
  unsigned int cc0pi_sig   = 0; 
  unsigned int cc0pi_sel   = 0; 

  unsigned int all_tracks_contained   = 0;
  unsigned int precuts_passed         = 0;
  unsigned int ccinc_passed           = 0;
  unsigned int ccinc_topology_passed  = 0;

  TopologyMap nc_map        = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc_map        = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_map     = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map     = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc1pinpi0_map = GeneralAnalysisHelper::GetCC1PiNPi0TopologyMap();
  TopologyMap cc0pi2p_map   = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap nc0pi_map     = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_map     = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nc1pinpi0_map = GeneralAnalysisHelper::GetNC1PiNPi0TopologyMap();
  TopologyMap ncnpi_map     = GeneralAnalysisHelper::GetNCNPiTopologyMap();
  TopologyMap nue_map       = GeneralAnalysisHelper::GetNuECCTopologyMap();
  TopologyMap nuenc_map     = GeneralAnalysisHelper::GetNuECCTopologyMap();
  std::vector<TopologyMap> cc_for_other{cc0pi_map, cc1pinpi0_map};
  std::vector<TopologyMap> nc_for_other{nc0pi_map, nc1pinpi0_map};
  TopologyMap ccoth_map   = GeneralAnalysisHelper::GetOtherTopologyMap(1,cc_for_other);
  TopologyMap ncoth_map   = GeneralAnalysisHelper::GetOtherTopologyMap(2,nc_for_other);
  
  std::vector< TopologyMap > reco_maps({cc_map, cc0pi_map});
  std::vector< TopologyMap > true_maps({cc0pi_map, cc1pi_map, ccoth_map, nc0pi_map, nc1pi_map, ncoth_map, nue_map});
  std::vector<std::string> reco_topology_names{"CC~Inc.", "CC~$0\\pi$"};
  std::vector<std::string> count_features{"CC~$0\\pi$", "CC~$1\\pi^{\\pm}$", "CC~Oth.", "NC~$0\\pi$", "NC~$1\\pi^{\\pm}$", "NC~Oth.", "$\\nu_{e}$", "True", "Selected", "Signal"};
  std::vector<std::string> reco_topology_names_notex{"CC Inc.", "CC 0Pi"};
  std::vector<std::string> count_features_notex{"CC 0Pi", "CC 1Pi", "CC Oth.", "NC 0pi", "NC 1Pi", "NC Oth.", "Nue", "True", "Selected", "Signal"};
  std::map< std::string, std::map<std::string, unsigned int> > topology_count_rates, topology_count_rates_notex;

  for(const std::string &f : count_features){
    std::map< std::string, unsigned int > count_rates;
    for(const std::string &t : reco_topology_names){
      count_rates.emplace(t,0); 
    }
    topology_count_rates.emplace(f,count_rates);
  }
  for(const std::string &f : count_features_notex){
    std::map< std::string, unsigned int > count_rates;
    for(const std::string &t : reco_topology_names_notex){
      count_rates.emplace(t,0); 
    }
    topology_count_rates_notex.emplace(f,count_rates);
  }

  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    selection::LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    for(const Event &e : events){
    
      // Boolean to check that we have selected a CC Inclusive using the method laid out in PIDMain
      bool cc_inclusive_passed = GeneralAnalysisHelper::PassedCCInclusive(e,detector);

      if(!e.IsRecoFiducial() || 
         !e.IsTrueFiducial() || 
         !GeneralAnalysisHelper::MaxOneLongEscapingTrack(e) || 
         !GeneralAnalysisHelper::MinOneRecoTrack(e)) continue;
      precuts_passed++;

      double enu = e.GetTrueNuEnergy();
      double pmu = 0.;
      double pth = 0.;
      // Now we are looking at selected CC Inclusive events, 
      // start trying to identify particle types
      if(cc_inclusive_passed){
        ccinc_passed++;
        for(unsigned int j = 0; j < reco_maps.size(); ++j){
          if(e.CheckRecoTopology(reco_maps.at(j))){          
            if(j == 0) //CCInc
              ccinc_topology_passed++;
            for(unsigned int k = 0; k < true_maps.size(); ++k){
              TString currMC(count_features_notex.at(k).data());
              if(currMC.Contains("NC")){
                if(e.CheckMCTopology(true_maps.at(k)) && e.CheckMCNeutrino(14)){
                  topology_count_rates.at(count_features.at(k)).at(reco_topology_names.at(j))++;
                  topology_count_rates_notex.at(count_features_notex.at(k)).at(reco_topology_names_notex.at(j))++;
                }
              }
              else if(currMC.Contains("Nue")){
                if(e.CheckMCNeutrino(12) || e.CheckMCNeutrino(-12)){
                  topology_count_rates.at(count_features.at(k)).at(reco_topology_names.at(j))++;
                  topology_count_rates_notex.at(count_features_notex.at(k)).at(reco_topology_names_notex.at(j))++;
                }
              }
              else{
                if(e.CheckMCTopology(true_maps.at(k))){
                  topology_count_rates.at(count_features.at(k)).at(reco_topology_names.at(j))++;
                  topology_count_rates_notex.at(count_features_notex.at(k)).at(reco_topology_names_notex.at(j))++;
                }
              }
            }
          }
        }
      }

      // Overall efficiencies 
      for(unsigned int j = 0; j < reco_maps.size(); ++j){
        if(e.CheckMCTopology(reco_maps.at(j))){
          topology_count_rates.at("True").at(reco_topology_names.at(j))++;
          topology_count_rates_notex.at("True").at(reco_topology_names_notex.at(j))++;

          if(j == 0 || j == 1){
            std::vector<float> pmus, pths;
            GeneralAnalysisHelper::GetMCModulusMomentumWithPdg(e,13,pmus);
            GeneralAnalysisHelper::GetMCCosThetaWithPdg(e,13,pths);
            if(pmus.size() > 1 || pths.size() > 1){
              std::cerr << " Error: Found more than 2 muons in a no-pileup sample " << std::endl;
            }
            pmu = pmus.at(0);
            pth = pths.at(0);
          }

          // CCInc true hist
          if(j == 0){
            h_true_ccinc_enu->Fill(enu);
            h_true_ccinc_pmu->Fill(pmu);
            h_true_ccinc_pth->Fill(pth);
          }
          
            // CC0pi true hist
          if(j == 1){
            h_true_cc0pi_enu->Fill(enu);
            h_true_cc0pi_pmu->Fill(pmu);
            h_true_cc0pi_pth->Fill(pth);
          }

          if(cc_inclusive_passed && e.CheckRecoTopology(reco_maps.at(j))){ 
            topology_count_rates.at("Signal").at(reco_topology_names.at(j))++;
            topology_count_rates_notex.at("Signal").at(reco_topology_names_notex.at(j))++;
            // CCInc signal hist
            if(j == 0){
              h_signal_ccinc_enu->Fill(enu);
              h_signal_ccinc_pmu->Fill(pmu);
              h_signal_ccinc_pth->Fill(pth);
            }

            // CC0pi signal hist
            if(j == 1){
              h_signal_cc0pi_enu->Fill(enu);
              h_signal_cc0pi_pmu->Fill(pmu);
              h_signal_cc0pi_pth->Fill(pth);
            }
          }
        }
      }

      // Overall purities
      for(unsigned int j = 0; j < reco_maps.size(); ++j){
        if(cc_inclusive_passed && e.CheckRecoTopology(reco_maps.at(j))){
          topology_count_rates.at("Selected").at(reco_topology_names.at(j))++;
          topology_count_rates_notex.at("Selected").at(reco_topology_names_notex.at(j))++;
        }
      }
    }
  }

  // POT scaling factor
  double potScale = totPOT / static_cast<double>(pot);
  
  // Histograms
  // Efficiency
  TH1D *h_eff_ccinc_enu = static_cast<TH1D*>(h_signal_ccinc_enu->Clone("h_eff_ccinc_enu"));
  TH1D *h_eff_cc0pi_enu = static_cast<TH1D*>(h_signal_cc0pi_enu->Clone("h_eff_cc0pi_enu"));

  TH1D *h_eff_ccinc_pmu = static_cast<TH1D*>(h_signal_ccinc_pmu->Clone("h_eff_ccinc_pmu"));
  TH1D *h_eff_cc0pi_pmu = static_cast<TH1D*>(h_signal_cc0pi_pmu->Clone("h_eff_cc0pi_pmu"));

  TH1D *h_eff_ccinc_pth = static_cast<TH1D*>(h_signal_ccinc_pth->Clone("h_eff_ccinc_pth"));
  TH1D *h_eff_cc0pi_pth = static_cast<TH1D*>(h_signal_cc0pi_pth->Clone("h_eff_cc0pi_pth"));

  h_eff_ccinc_enu->Divide(h_true_ccinc_enu);
  h_eff_cc0pi_enu->Divide(h_true_cc0pi_enu);
                                         
  h_eff_ccinc_pmu->Divide(h_true_ccinc_pmu);
  h_eff_cc0pi_pmu->Divide(h_true_cc0pi_pmu);
                                         
  h_eff_ccinc_pth->Divide(h_true_ccinc_pth);
  h_eff_cc0pi_pth->Divide(h_true_cc0pi_pth);

  SetHistErrors(h_eff_ccinc_enu, h_true_ccinc_enu);
  SetHistErrors(h_eff_cc0pi_enu, h_true_cc0pi_enu);

  SetHistErrors(h_eff_ccinc_pmu, h_true_ccinc_pmu);
  SetHistErrors(h_eff_cc0pi_pmu, h_true_cc0pi_pmu);

  SetHistErrors(h_eff_ccinc_pth, h_true_ccinc_pth);
  SetHistErrors(h_eff_cc0pi_pth, h_true_cc0pi_pth);

  // Style
  std::vector<TH1D*> h_inc  = {h_eff_ccinc_enu,  h_eff_ccinc_pmu,  h_eff_ccinc_pth};
  std::vector<TH1D*> h_0pi  = {h_eff_cc0pi_enu,  h_eff_cc0pi_pmu,  h_eff_cc0pi_pth};

  SetHistogramStyle(h_inc,  0, 1, 905, 905, 2, 132, 1, 1.2, false);
  SetHistogramStyle(h_0pi,  0, 1, 867, 867, 2, 132, 1, 1.2, false);

  TLegend * l = new TLegend(0.24,0.94,0.85,0.99);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.045);
  l->SetNColumns(2);
  l->SetMargin(0.2);
  l->SetTextFont(132);

  TCanvas *c = new TCanvas("c","",900,900);
  c->SetLeftMargin  (0.138796 );
  c->SetRightMargin (0.0334448);
  c->SetBottomMargin(0.132404 );
  c->SetTopMargin   (0.0734 );

  // Efficiencies
  l->AddEntry(h_eff_ccinc_enu,"#nu_{#mu} CC Inc.","ep");
  l->AddEntry(h_eff_cc0pi_enu,"#nu_{#mu} CC 0#pi","ep");

  h_eff_ccinc_enu->SetTitle("");
  h_eff_ccinc_enu->GetYaxis()->SetTitle("Efficiency");
  h_eff_ccinc_enu->GetXaxis()->SetTitle("True neutrino energy [GeV]");
  h_eff_ccinc_enu->GetYaxis()->SetRangeUser(0.0,1.1);
  h_eff_ccinc_enu->Draw("E1 X0");
  h_eff_cc0pi_enu->Draw("E1 X0 same");
  l->Draw("same");
  
  c->SaveAs((plots_location+"ccpassed/enu_pid_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"ccpassed/enu_pid_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"ccpassed/enu_pid_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();
  
  l->AddEntry(h_eff_ccinc_pmu,"#nu_{#mu} CC Inc.","ep");
  l->AddEntry(h_eff_cc0pi_pmu,"#nu_{#mu} CC 0#pi","ep");

  h_eff_ccinc_pmu->SetTitle("");
  h_eff_ccinc_pmu->GetYaxis()->SetTitle("Efficiency");
  h_eff_ccinc_pmu->GetXaxis()->SetTitle("True muon momentum [GeV]");
  h_eff_ccinc_pmu->GetYaxis()->SetRangeUser(0.0,1.1);
  h_eff_ccinc_pmu->Draw("E1 X0");
  h_eff_cc0pi_pmu->Draw("E1 X0 same");
  l->Draw("same");
  
  c->SaveAs((plots_location+"ccpassed/pmu_pid_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"ccpassed/pmu_pid_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"ccpassed/pmu_pid_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();
  
  l->AddEntry(h_eff_ccinc_pth,"#nu_{#mu} CC Inc.","ep");
  l->AddEntry(h_eff_cc0pi_pth,"#nu_{#mu} CC 0#pi","ep");

  h_eff_ccinc_pth->SetTitle("");
  h_eff_ccinc_pth->GetYaxis()->SetTitle("Efficiency");
  h_eff_ccinc_pth->GetXaxis()->SetTitle("True cos#theta_{#mu}");
  h_eff_ccinc_pth->GetYaxis()->SetRangeUser(0.0,1.1);
  h_eff_ccinc_pth->Draw("E1 X0");
  h_eff_cc0pi_pth->Draw("E1 X0 same");
  l->Draw("same");
  
  c->SaveAs((plots_location+"ccpassed/pth_pid_efficiency_overlay.root").c_str());
  c->SaveAs((plots_location+"ccpassed/pth_pid_efficiency_overlay.png").c_str());
  c->SaveAs((plots_location+"ccpassed/pth_pid_efficiency_overlay.pdf").c_str());
  c->Clear();
  l->Clear();

  // Files to hold particle statistics
  ofstream file,fileTeX;
  file.open(stats_location+"thesis_topological_breakdown.txt");
  fileTeX.open(stats_location+"thesis_topological_breakdown.tex");

  // Text file header
  file << "=============================================================================" << std::endl;
  file << " POT scaling factor for the current detector              : " << potScale << std::endl;
  file << "-----------------------------------------------------------------------------" << std::endl;
  file << " Total number of events passing the reconstructable cuts  : " << precuts_passed << std::endl;
  file << " Total number of events passing the CC Inclusive cuts     : " << ccinc_passed << std::endl;
  file << " Total number of events selected as CC Inclusive topology : " << ccinc_topology_passed << std::endl;
  file << "-----------------------------------------------------------------------------" << std::endl;
  
  // TeX file header
  fileTeX << "\\begin{table}[h!] " << std::endl;
  fileTeX << "  \\small " << std::endl;
  fileTeX << "  \\centering " << std::endl;
  fileTeX << "  \\renewcommand{\\arraystretch}{1.4}" << std::endl;
  fileTeX << "  \\begin{tabular}{ m{2cm} * {" << reco_topology_names.size() << "}{ >{\\centering\\arraybackslash}m{3cm} } }" << endl;
  fileTeX << "    \\hline" << endl;
  
  //Text file first row
  file << std::setw(12) << "True \\ Reco" << " || ";
  for(const std::string &t : reco_topology_names_notex){
    file << std::setw(10) << t;
  }
  file << std::endl;
  file << "-----------------------------------------------------------------------------" << std::endl;

  //Tex file first row
  fileTeX << "    Reco~$\\rightarrow$ &"; 
  for(unsigned int i = 0; i < reco_topology_names.size() - 1; ++i){
    fileTeX << "\\multirow{2}{*}{" << reco_topology_names.at(i) << "} &";
  }
  fileTeX << "\\multirow{2}{*}{" << reco_topology_names.at(reco_topology_names.size()-1) << "} \\\\" << std::endl;
  fileTeX << "    $\\downarrow$~True &";
  for(unsigned int i = 0; i < reco_topology_names.size() - 1; ++i){
    fileTeX << " & ";
  }
  fileTeX << " \\\\" << std::endl;
  fileTeX << "    \\hline" << endl;

  // Text file contents
  for(const std::string &f : count_features_notex){
    if(f == "True"){
      file << "-----------------------------------------------------------------------------" << std::endl;
    }
    file << std::setw(12) << f << " || ";
    for(const std::string &t : reco_topology_names_notex){
      file << std::setw(10) << std::fixed << std::setprecision(2) << static_cast<int>(topology_count_rates_notex.at(f).at(t)*potScale);
    }
    file << std::endl;
  }
  file << "-----------------------------------------------------------------------------" << std::endl;
  
  // Tex file contents
  for(const std::string &f : count_features){
    if(f == "True"){
      fileTeX << "    \\hline" << std::endl;
    }
    fileTeX << "    " << f << " & ";
    for(unsigned int i = 0; i < reco_topology_names.size()-1; ++i){
      const std::string t = reco_topology_names.at(i);
      fileTeX << "\\num{ " << std::fixed << std::setprecision(2) << static_cast<int>(topology_count_rates.at(f).at(t)*potScale) << " } &";
    }
    fileTeX << "\\num{ " << std::fixed << std::setprecision(2) << static_cast<int>(topology_count_rates.at(f).at(reco_topology_names.at(reco_topology_names.size()-1))*potScale) << " } \\\\ " << std::endl;
  }
  fileTeX << "    \\hline" << std::endl;

  // Text file efficiency
  file << std::setw(12) << "Efficiency " << " || ";
  for(const std::string &t : reco_topology_names_notex){
    double eff = topology_count_rates_notex.at("Signal").at(t) / static_cast<double>(topology_count_rates_notex.at("True").at(t));
    double err = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates_notex.at("Signal").at(t),topology_count_rates_notex.at("True").at(t));
    file << std::setw(10) << std::fixed << std::setprecision(2) << eff << " +/-" << std::fixed << std::setprecision(2) << err;
  }
  file << std::endl;
  file << std::setw(12) << "Purity" << " || ";
  for(const std::string &t : reco_topology_names_notex){
    double pur = topology_count_rates_notex.at("Signal").at(t) / static_cast<double>(topology_count_rates_notex.at("Selected").at(t));
    double err = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates_notex.at("Signal").at(t),topology_count_rates_notex.at("Selected").at(t));
    file << std::setw(10) << std::fixed << std::setprecision(2) << pur << " +/-" << std::fixed << std::setprecision(2) << err;
  }
  file << std::endl;

  // Tex file efficiency
  fileTeX << "    Efficiency & ";
  for(unsigned int i = 0; i < reco_topology_names.size()-1; ++i){
    const std::string t = reco_topology_names.at(i);
    double eff = topology_count_rates.at("Signal").at(t) / static_cast<double>(topology_count_rates.at("True").at(t));
    double err = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates.at("Signal").at(t),topology_count_rates.at("True").at(t));
    fileTeX << std::fixed << std::setprecision(2) << eff*100 << "~$\\pm$~" <<  std::fixed << std::setprecision(2) <<err*100 << "~\\% &";
  }
  const std::string te = reco_topology_names.at(reco_topology_names.size()-1);
  double effEnd = topology_count_rates.at("Signal").at(te) / static_cast<double>(topology_count_rates.at("True").at(te));
  double errEnd = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates.at("Signal").at(te),topology_count_rates.at("True").at(te));
  fileTeX << std::fixed << std::setprecision(2) << effEnd*100 << "~$\\pm$~" <<  std::fixed << std::setprecision(2) <<errEnd*100 << "~\\% \\\\" << std::endl;

  fileTeX << "    Purity & ";
  for(unsigned int i = 0; i < reco_topology_names.size()-1; ++i){
    const std::string t = reco_topology_names.at(i);
    double pur = topology_count_rates.at("Signal").at(t) / static_cast<double>(topology_count_rates.at("Selected").at(t));
    double err = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates.at("Signal").at(t),topology_count_rates.at("Selected").at(t));
    fileTeX << std::fixed << std::setprecision(2) << pur*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << err*100 << "~\\% &";
  }
  double purEnd = topology_count_rates.at("Signal").at(te) / static_cast<double>(topology_count_rates.at("Selected").at(te));
  errEnd = GeneralAnalysisHelper::RatioUncertainty(topology_count_rates.at("Signal").at(te),topology_count_rates.at("Selected").at(te));
  fileTeX << std::fixed << std::setprecision(2) << purEnd*100 << "~$\\pm$~" << std::fixed << std::setprecision(2) << errEnd*100 << "~\\% \\\\" << std::endl;
  fileTeX << "    \\hline" << std::endl;

  // Text file end
  file << "=============================================================================" << std::endl;
  
  // Tex file end
  fileTeX << "  \\end{tabular}"<< endl;
  fileTeX << "\\end{table}" << std::endl;
  
  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

// Set the errors on the efficiency plots given the numerator and denomenator
void SetHistErrors(TH1D *h, TH1D *den){

  // Make sure the binning is the same
  if(h->GetNbinsX() != den->GetNbinsX()){
    std::cerr << " Error: Number of bins in the effiency histogram does not match that of the denomenator: " << std::endl;
    std::cerr << "   Histogram: " << h->GetNbinsX() << ", denomenator: " << den->GetNbinsX() << std::endl;
    exit(1);  
  }

  // Classical approach doesn't work at the limits ( e = 0,1 )
  //
  // Instead use Bayesean Approach from: https://cds.cern.ch/record/1010669/files/0701199.pdf?version=1
  //
  for(int b = 1; b <= h->GetNbinsX(); ++b){
    double num = h->GetBinContent(b) * den->GetBinContent(b);
    double denom = den->GetBinContent(b);
    h->SetBinError(b,GeneralAnalysisHelper::RatioUncertainty(num, denom));
    
//    std::cout << " Selected: " << num << ", total: " << denom << std::endl;
//    std::cout << " Eff: " << h->GetBinContent(b) << " +/- " << sig << std::endl;
  }
}
