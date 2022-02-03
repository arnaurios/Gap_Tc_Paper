#!/bin/bash

# CHOOSE DENSITIES AND TEMPERATURES TO BE RUN AND SUBMIT A BUNCH OF JOBS
dens=(0.02 0.032 0.03 0.04 0.06 0.08 0.10 0.12 0.128 0.14 0.16 0.18 0.20 0.24 0.28 0.32 0.36 0.40 0.44 0.48);
original_folder=$PWD;

# length of the arrays
nden=${#dens[*]};

# loop over densities and temperatures
it=0; 
while test $it -lt $nden
do
# Get into folder and run
	runfolder="/Users/arnau/Volumes/Surrey/project/work/NSCL_project/SCGF/neumat/N3LO_noSRG_TBFcorrelated_intreg/g21a_rho${dens[it]}"
	tofolder="$PWD/g2a_rho${dens[it]}/"
	echo $runfolder $tofolder
	if [ ! -e "$tofolder" ]; then
    # Control will enter here if $DIRECTORY doesn't exist
            mkdir $tofolder
        fi

	cd $runfolder
#	cp wim*.dat $tofolder
#	cp *.eps $tofolder
	cp momdis* $tofolder

	
	cd $tofolder
	
	cd $original_folder
    it=`expr $it + 1`
done
