source /grid/fermiapp/products/sbnd/setup_sbnd.sh
source /cvmfs/larsoft.opensciencegrid.org/products/setup && retval="$?"
setup root v6_18_04b -q e17:prof
cd $SELECTION_TOOL
