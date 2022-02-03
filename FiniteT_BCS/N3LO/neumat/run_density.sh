#!/bin/bash
# CHOOSE DENSITY
dens=(0.30) # 0.40 0.50 0.60 0.70 0.80 0.90 1.00);
dens=(0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00);
dens=(1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40);
#temp_array="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1"
#itemp=1

original_folder=$PWD;


for rho in ${dens[@]}
do
  echo $rho
  echo $temp_array
  echo $itemp

    input_file=inp.$rho
    echo $input_file
    if [ -f "$input_file" ]; then
       rm $input_file
    fi
    #echo 1,$rho >> $input_file
    echo $rho,0.,1 >> $input_file
    #echo $itemp $temp_array >> $input_file
    echo 0.0001,0.1,16 >> $input_file
    echo 1S0 >> $input_file
    #echo T >> $input_file

    cat $input_file
    ./bcs.x < $input_file
#    python3 fig_denominator_temp.py $rho

    savefolder="gaps/kf_"$rho
    echo $savefolder
    if [ ! -d "$savefolder" ]; then
       mkdir $savefolder
    else
       rm -r $savefolder
       mkdir $savefolder
    fi

    cp *.dat $savefolder"/"

done
