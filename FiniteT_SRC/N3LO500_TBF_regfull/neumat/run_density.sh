#!/bin/bash
# CHOOSE DENSITY
dens=(.010 .015 0.02 0.03 .032 0.04 0.05 0.06 0.08);

original_folder=$PWD;
#temp_array=("0")
for rho in ${dens[@]}
do
  echo $rho
  echo g2a_rho${rho}_t*/
  temp_array="0"
  itemp=1
  for f in g2a_rho${rho}_t*/
    do [ -d "$f" ]
      echo $f is indeed a folder ;
      temp=${f:13:3}
      echo $temp
      #temp_array+=("${temp_array[@]}" "{$temp}")
      temp_array=$temp_array" $temp"
      itemp=$((itemp + 1))
    done
    echo $temp_array
    echo $itemp

    input_file=inp.$rho
    echo $input_file
    if [ -f "$input_file" ]; then
       rm $input_file
    fi
    echo 1,$rho >> $input_file
    echo 1S0 >> $input_file
    echo T >> $input_file
    echo $itemp,$temp_array >> $input_file

    cat $input_file
    ./bcs.x < $input_file
    python3 fig_denominator_temp.py $rho


    savefolder="gaps/"$rho
    echo $savefolder
    if [ ! -d "$savefolder" ]; then
       mkdir $savefolder
    else
       rm -r $savefolder
       mkdir $savefolder
    fi

    cp *.dat $savefolder"/"
    cp wim_den$rho.pdf $savefolder"/"

done
