#include "../include/Setup.h"

namespace selection{
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void LoadAllEventsInFile(const std::string &input_loc, 
      const std::string &input_name, 
      EventSelectionTool::EventList &events, 
      const unsigned int &file_index, 
      double &pot, 
      std::vector<int> &exceptions,
      const Geometry &fid,
      const Geometry &av,
      const bool &runEverything){
    std::vector<int>::const_iterator it = std::find(exceptions.begin(), exceptions.end(),file_index);
    if(it != exceptions.end()) return;

    // Get the filenames
    const std::string name = input_loc+"/"+std::to_string(file_index)+"/"+input_name;
    //const std::string name = "/pnfs/sbnd/persistent/users/rsjones/sbnd_selection_doms_files_190620/selection/"+std::to_string(file_index)+"/selection_tree.root";

    double temp_pot = 0.;
    EventSelectionTool::LoadEventList(name, events, file_index, temp_pot, fid, av, runEverything);
    pot += temp_pot;
  } // LoadAllEventsInFile
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void LoadAllEventsInFile(const std::string &input_loc, 
      const std::string &input_name, 
      EventSelectionTool::EventList &events, 
      const unsigned int &file_index, 
      double &pot, 
      std::vector<int> &exceptions,
      const Geometry &fid,
      const Geometry &av){
  
    LoadAllEventsInFile(input_loc, input_name, events, file_index, pot, exceptions, fid, av, false);

  } // LoadAllEventsInFile
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void FillExceptions(const char *file_name, std::vector<int> &exceptions){
    // Read in txt file of list of empty input directories
    std::fstream file(file_name);
    if(!file){
      std::cout << " No exceptions file given" << std::endl;
      return;
    }
    else{
      std::string s_exc;
      while (std::getline(file, s_exc)) {
        unsigned int i_exc;
        std::istringstream ss_exc(s_exc);
        ss_exc >> i_exc;
        exceptions.push_back(i_exc);
        ss_exc.str(std::string());
        s_exc.clear();
      }
    }
  } // FillExceptions
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void GetTime(time_t &rawtime){
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    std::cout << " Local time and date:  " << asctime(timeinfo)                << std::endl;
  } // GetTime
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void GetTotalTime(time_t &starttime, time_t &endtime){
    double timediff = difftime(endtime,starttime);

    // Calculate the time in a nice format
    int seconds_per_minute = 60; 
    int seconds_per_hour   = 3600; 
    int hours              = timediff / seconds_per_hour;
    int minutes            = (timediff - (hours * seconds_per_hour)) / (seconds_per_minute);
    int seconds            = timediff - (hours * seconds_per_hour) - (minutes * seconds_per_minute);

    // Now print
    std::cout << " Total time taken: "    << std::setw(4) << hours << " hours, ";
    std::cout                             << std::setw(4) << minutes << " minutes, ";
    std::cout                             << std::setw(4) << seconds << " seconds." << std::endl;
  } // GetTotalTime
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  void SetHistogramStyle(std::vector<TH1D*> hists, 
      const int &fillstyle, 
      const int &linestyle, 
      const int &fillcolour,
      const int &linecolour,
      const int &linewidth,
      const int &font,
      const double &xoffset,
      const double &yoffset,
      const bool &scale){
    // Loop over the histograms and get styling!
    for(TH1D *h : hists){
      h->SetStats(0);
      h->SetFillStyle(fillstyle);
      h->SetLineStyle(linestyle);
      h->SetFillColor(fillcolour);
      h->SetLineColor(linecolour);
      h->SetLineWidth(linewidth);
      h->SetTitleFont(font);
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleFont(font);
      h->GetXaxis()->SetLabelFont(font);
      h->GetYaxis()->SetTitleFont(font);
      h->GetYaxis()->SetLabelFont(font);
      h->GetXaxis()->SetTitleOffset(xoffset);
      h->GetYaxis()->SetTitleOffset(yoffset);
      h->GetYaxis()->SetMaxDigits(3);
      if(!scale) continue;
      if(h->Integral() > 0)
        h->Scale(1/h->Integral());
      else
        std::cout << " Can't scale this histogram: The integral is 0" << std::endl;
    }
  }
  // --------------------------------------------------------------------------------------------------------------------------------------------------
}
