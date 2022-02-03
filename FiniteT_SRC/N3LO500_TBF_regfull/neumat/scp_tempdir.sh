#!/bin/bash

# CHOOSE DENSITIES AND TEMPERATURES TO BE RUN AND SUBMIT A BUNCH OF JOBS
dens=(0.03);
tfinal=(0.3);


original_folder=$PWD;
# length of the arrays
nden=${#dens[*]};

user=arnaurios
server=goeppert.fqa.ub.edu
dir_base=/home/arnaurios/project/NSCL_project/SCGF/neumat/N3LO500_noSRG_TBFcorrelated_intreg/auto_extrap_gII_disp/

# loop over densities and temperatures
it=0;
while test $it -lt $nden
do
# Get into folder and run
    for temp in ${tfinal[@]}
    do
      if [[ $temp == 0 ]]; then
        runfolder="g2a_rho${dens[it]: -4}"
        basefolder="g2a_rho${dens[it]}"
        dir_root=$dir_base"data"
      else
        runfolder="g2a_rho${dens[it]: -4}_t$temp"
        basefolder="g2a_rho${dens[it]}_t$temp"
        dir_root=$dir_base"data_tf"
      fi

      echo $basefolder
      echo $runfolder

      diraux=$dir_root/$basefolder
      echo $diraux
      if ssh $user@$server [ -d "$diraux" ]; then

        echo $diraux
        if [ ! -d "$runfolder" ]; then
           mkdir $runfolder
        else
           rm -r $runfolder
           mkdir $runfolder
        fi

        # rsync -aP $user@$server:$diraux/*.eps $runfolder
        rsync -aP $user@$server:$diraux/*.pdf $runfolder
        rsync -aP $user@$server:$diraux/wim_t0.dat $runfolder
      fi

    done
    it=`expr $it + 1`
done
