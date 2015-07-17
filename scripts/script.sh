#!/bin/bash
echo '#!/usr/bin/env gnuplot' > plot.gnu
echo "set size ratio 1" >> plot.gnu
echo "set cbrange [0.8:1.6]" >> plot.gnu
echo "set pm3d map" >> plot.gnu
#echo "set palette gray" >> plot.gnu
#for i in {0..500..10}     # replace 1 and 10 with any number to like to. For more check the syntax  
for (( i=0; i <= $1; i+=$2 ))
do
#	echo "splot '/Users/niktre/Documents/MPIP/LB_codes/Jasna/output_drop_2015_07_16/density${i}' matrix ">>plot.gnu
	echo "splot '/Users/niktre/Documents/MPIP/LB_codes/Jasna/output/density${i}' matrix ">>plot.gnu
	echo "pause 0.1" >>plot.gnu      # one second interval between plots
#        echo "pause -1" >>plot.gnu    # press enter for each transition 
done

