#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include "../include/Plane.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <time.h>
#include <ctime>
#include <stdexcept>
#include "TROOT.h"
#include "TMath.h"
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

namespace selection{

  /**
   * @brief  Load a single file and fill the event list 
   *
   * @param  input_loc location of the input file subdirectories
   * @param  input_name name of the root file to open
   * @param  events empty event list to fill
   * @param  file_index index of the current file
   * @param  pot empty pot object to add to
   * @param  exceptions list of files to skip (if they are empty or broken)
   * @param  fid fiducial volume geometry object
   * @param  av active volume geometry object
   *
   */
  void LoadAllEventsInFile(const std::string &input_loc, 
                           const std::string &input_name, 
                           EventSelectionTool::EventList &events, 
                           const unsigned int &file_index, 
                           double &pot, 
                           std::vector<int> &exceptions,
                           const Geometry &fid,
                           const Geometry &av);

  /**
   * @brief  Fill the list of files to skip
   *
   * @param  file input exception list
   * @param  exceptions empty vector to fill
   *
   */
  void FillExceptions(const char *file_name, std::vector<int> &exceptions);

  /**
   * @brief  Get the time to monitor the efficiency of the programme
   *
   * @param  rawtime time object to read
   *
   */
  void GetTime(time_t &rawtime);
  
  /**
   * @brief  Get the total time taken and print nicely
   *
   * @param  starttime start time object to read
   * @param  endtime end time object to read
   *
   */
  void GetTotalTime(time_t &starttime, time_t &endtime);
} // namespace

