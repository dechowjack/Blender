#!/bin/bash
start=$1
counter=0
while [ $counter -le 199 ]
do
	run=$(($start+$counter))
        printdir="Pix$run"
        cpdir="Inputs/Pix$run"
        echo "Running $printdir"
        cp "$PBS_O_WORKDIR/$cpdir"/*.txt ./
        /apps/julia/1.1.1/bin/julia Blender_Algorithm_v5.0.jl . >log.txt
        cp *.txt "$PBS_O_WORKDIR/$cpdir"/ 
	rm *.txt
        ((counter++))
done
