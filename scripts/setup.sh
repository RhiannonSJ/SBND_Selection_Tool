source /grid/fermiapp/products/sbnd/setup_sbnd.sh
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup && retval="$?"
setup root v6_18_04b -q e17:prof
cd $SELECTION_TOOL
