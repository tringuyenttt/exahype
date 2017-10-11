# set binary and configuration file
EXABINARY=$(ls ./ExaHyPE-*) # don't forget trailing slash ie. ./ExaHype_bla
EXACONF=../MHD_AlfenWave.exahype

# setup simulation parameters, apply workaround
export exahype_parameter_workaround="yes"
export exahype_initialdata=$(grep constants $EXACONF | sed -n "s/^.*initialdata:\(\S*\)}.*/\1/p")

# make a log file
exec > >(awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' |  tee run-$(hostname).log) 2>&1

echo "This is $0 on $(hostname) at $(date)"
echo "Running with exahype_parameter_workaround and initialdata=$exahype_initialdata"
echo "This will be a great fun"
echo

# make sure stuff is there which is needed by the exahype binary
mkdir -p output

# delete old output
rm *vtk *log-file

# run the stuff
time $EXABINARY $EXACONF

# pack all the output files
tar cvfz results.tar.gzip *.vtk
