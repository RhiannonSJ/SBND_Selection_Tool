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

  unsigned int cc1pi_cc0pi = 0;
  unsigned int cc1pi_cc1pi = 0;
  unsigned int cc1pi_ccoth = 0;
  unsigned int cc1pi_nc0pi = 0;
  unsigned int cc1pi_nc1pi = 0;
  unsigned int cc1pi_ncoth = 0;
  unsigned int cc1pi_nue   = 0;
  unsigned int cc1pi_true  = 0; 
  unsigned int cc1pi_sig   = 0; 
  unsigned int cc1pi_sel   = 0; 

  unsigned int ncinc_cc0pi = 0;
  unsigned int ncinc_cc1pi = 0;
  unsigned int ncinc_ccoth = 0;
  unsigned int ncinc_nc0pi = 0;
  unsigned int ncinc_nc1pi = 0;
  unsigned int ncinc_ncoth = 0;
  unsigned int ncinc_nue   = 0;
  unsigned int ncinc_true  = 0; 
  unsigned int ncinc_sig   = 0; 
  unsigned int ncinc_sel   = 0; 

  unsigned int nc0pi_cc0pi = 0;
  unsigned int nc0pi_cc1pi = 0;
  unsigned int nc0pi_ccoth = 0;
  unsigned int nc0pi_nc0pi = 0;
  unsigned int nc0pi_nc1pi = 0;
  unsigned int nc0pi_ncoth = 0;
  unsigned int nc0pi_nue   = 0;
  unsigned int nc0pi_true  = 0; 
  unsigned int nc0pi_sig   = 0; 
  unsigned int nc0pi_sel   = 0; 

  unsigned int nc1pi_cc0pi = 0;
  unsigned int nc1pi_cc1pi = 0;
  unsigned int nc1pi_ccoth = 0;
  unsigned int nc1pi_nc0pi = 0;
  unsigned int nc1pi_nc1pi = 0;
  unsigned int nc1pi_ncoth = 0;
  unsigned int nc1pi_nue   = 0;
  unsigned int nc1pi_true  = 0; 
  unsigned int nc1pi_sig   = 0; 
  unsigned int nc1pi_sel   = 0;

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
  std::vector<TopologyMap> cc_for_other{cc0pi_map, cc1pinpi0_map};
  std::vector<TopologyMap> nc_for_other{nc0pi_map, nc1pinpi0_map};
  TopologyMap ccoth_map   = GeneralAnalysisHelper::GetOtherTopologyMap(1,cc_for_other);
  TopologyMap ncoth_map   = GeneralAnalysisHelper::GetOtherTopologyMap(2,nc_for_other);
  
  std::vector< TopologyMap > reco_maps({cc_map, cc0pi_map, cc1pi_map, nc_map, nc0pi_map, ncnpi_map});
  std::vector< TopologyMap > true_maps({cc0pi_map, cc1pi_map, ccoth_map, nc0pi_map, nc1pi_map, ncoth_map, nue_map});
  std::vector<std::string> reco_topology_names{"CC~Inc.", "CC~$0\\pi$", "CC~$1\\pi^{\\pm}$", "NC~Inc.", "NC~$0\\pi$", "NC~n$\\pi^{\\pm,0}$"};
  std::vector<std::string> count_features{"CC~$0\\pi$", "CC~$1\\pi^{\\pm}$", "CC~Oth.", "NC~$0\\pi$", "NC~$1\\pi^{\\pm}$", "NC~Oth.", "$\\nu_{e}$", "True", "Signal", "Selected"};
  std::vector<std::string> reco_topology_names_notex{"CC Inc.", "CC 0Pi", "CC 1Pi", "NC Inc.", "NC 0Pi", "NC NPi"};
  std::vector<std::string> count_features_notex{"CC 0Pi", "CC 1Pi", "CC Oth.", "NC 0pi", "NC 1Pi", "NC Oth.", "Nue", "True", "Signal", "Selected"};
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

      // Now we are looking at selected CC Inclusive events, 
      // start trying to identify particle types
      if(cc_inclusive_passed){
        ccinc_passed++;
        for(unsigned int j = 0; j < reco_maps.size(); ++j){
          if(e.CheckRecoTopology(reco_maps.at(j))){          
            if(j == 0) //CCInc
              ccinc_topology_passed++;
            for(unsigned int k = 0; k < true_maps.size(); ++k){
              if(e.CheckMCTopology(true_maps.at(k))){
                topology_count_rates.at(count_features.at(k)).at(reco_topology_names.at(j))++;
                topology_count_rates_notex.at(count_features_notex.at(k)).at(reco_topology_names_notex.at(j))++;
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
          if(cc_inclusive_passed && e.CheckRecoTopology(reco_maps.at(j))){ 
            topology_count_rates.at("Signal").at(reco_topology_names.at(j))++;
            topology_count_rates_notex.at("Signal").at(reco_topology_names_notex.at(j))++;
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
  
  // Files to hold particle statistics
  ofstream file,fileTeX;
  file.open(stats_location+"full_topological_breakdown.txt");
  fileTeX.open(stats_location+"full_topological_breakdown.tex");

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
  fileTeX << "  \\begin{tabular}{ m{1.5cm} * {" << reco_topology_names.size() << "}{ >{\\centering\\arraybackslash}m{2cm} } }" << endl;
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
      file << std::setw(10) << std::setprecision(5) << static_cast<int>(topology_count_rates_notex.at(f).at(t)*potScale);
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
      fileTeX << "\\num{ " << std::setprecision(5) << static_cast<int>(topology_count_rates.at(f).at(t)*potScale) << " } &";
    }
    fileTeX << "\\num{ " << std::setprecision(5) << static_cast<int>(topology_count_rates.at(f).at(reco_topology_names.at(reco_topology_names.size()-1))*potScale) << " } \\\\ " << std::endl;
  }
  fileTeX << "    \\hline" << std::endl;

  // Text file efficiency
  file << std::setw(12) << "Efficiency " << " || ";
  for(const std::string &t : reco_topology_names_notex){
    double eff = topology_count_rates_notex.at("Signal").at(t) / static_cast<double>(topology_count_rates_notex.at("True").at(t));
    file << std::setw(10) << std::setprecision(4) << eff;
  }
  file << std::endl;
  file << std::setw(12) << "Purity" << " || ";
  for(const std::string &t : reco_topology_names_notex){
    double pur = topology_count_rates_notex.at("Signal").at(t) / static_cast<double>(topology_count_rates_notex.at("Selected").at(t));
    file << std::setw(10) << std::setprecision(4) << pur;
  }
  file << std::endl;

  // Tex file efficiency
  fileTeX << "    Efficiency & ";
  for(unsigned int i = 0; i < reco_topology_names.size()-1; ++i){
    const std::string t = reco_topology_names.at(i);
    double eff = topology_count_rates.at("Signal").at(t) / static_cast<double>(topology_count_rates.at("True").at(t));
    fileTeX << std::setprecision(4) << eff*100 << "~\\% &";
  }
  const std::string te = reco_topology_names.at(reco_topology_names.size()-1);
  double effEnd = topology_count_rates.at("Signal").at(te) / static_cast<double>(topology_count_rates.at("True").at(te));
  fileTeX << std::setprecision(4) << effEnd*100 << "~\\% \\\\" << std::endl;

  fileTeX << "    Purity & ";
  for(unsigned int i = 0; i < reco_topology_names.size()-1; ++i){
    const std::string t = reco_topology_names.at(i);
    double pur = topology_count_rates.at("Signal").at(t) / static_cast<double>(topology_count_rates.at("Selected").at(t));
    fileTeX << std::setprecision(4) << pur*100 << "~\\% &";
  }
  double purEnd = topology_count_rates.at("Signal").at(te) / static_cast<double>(topology_count_rates.at("Selected").at(te));
  fileTeX << std::setprecision(4) << purEnd*100 << "~\\% \\\\" << std::endl;
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

