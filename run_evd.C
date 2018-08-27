{
/*
 * Example list of source and main files from RooUnfold_Kinematics
 *
 * NB: The main is compiled last since it depends on everything else
 *
 *
gROOT->ProcessLine(".L source_files/all_classes.cpp+");
gROOT->ProcessLine(".L source_files/cov_func.cpp+");
gROOT->ProcessLine(".L source_files/unfolding_functions.cpp+");
gROOT->ProcessLine(".L source_files/cross_section.cpp+");
gROOT->ProcessLine(".L source_files/slicing.cpp+");
gROOT->ProcessLine(".L source_files/run_all_per_signal.cpp+");
gROOT->ProcessLine(".L source_files/main_cc0pi.cpp+");
*
*
*/

gROOT->ProcessLine(".L srcs/Particle.cpp+");
gROOT->ProcessLine(".L srcs/Event.cpp+");
gROOT->ProcessLine(".L srcs/EventSelectionTool.cpp+");
gROOT->ProcessLine(".L srcs/GeneralAnalysisHelper.cpp+");
gROOT->ProcessLine(".L test/EvdMain.cpp+");
}
