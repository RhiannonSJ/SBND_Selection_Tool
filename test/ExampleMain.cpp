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
  std::string plots_location  = "";
  unsigned int total_files = 0;
  std::vector<double> minx_fid, miny_fid, minz_fid;
  std::vector<double> maxx_fid, maxy_fid, maxz_fid;
  std::vector<double> minx_av, miny_av, minz_av;
  std::vector<double> maxx_av, maxy_av, maxz_av;

  p->getValue("InputFileLocation",input_location);
  p->getValue("InputFileName",    input_filename);
  p->getValue("ExceptionsFile",   exceptions_file);
  p->getValue("StatFileLocation", stats_location); 
  p->getValue("PlotFileLocation", plots_location); 
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

  // First, ensure all tracks are contained
  for( unsigned int i = 0; i < total_files; ++i ){
    EventSelectionTool::EventList events;
    LoadAllEventsInFile(input_location, input_filename, events, i, pot, exceptions, fiducial, active);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
    // Now loop over the events
    for(const Event &e : events){
    }
  }
  std::cout << "-----------------------------------------------------------" << std::endl;
  time_t rawtime_end;
  GetTime(rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  GetTotalTime(rawtime, rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()


