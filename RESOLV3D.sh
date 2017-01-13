#!/bin/bash

if [ "$#" -ne 7 ]
then
	echo "Invalid parameters."
	echo "Usage: sh runCONCEPT3D <input file> <output prefix> <iterations>"
else
	echo "Running CONCEPT3D with following parameters:"
	echo "===> Working directory: ./"
	echo "===> Input file: " ${1}
	echo "===> Output prefix: " ${2}
	echo "===> Iterations: " ${3}
	echo "===> Cutoff: " ${4}
	echo "===> Genome: " ${5}
	echo "===> Key type: " ${6}
	echo "===> Background genes: " ${7}
	echo `date`
	R --no-save --no-restore --max-ppsize=500000 --slave --args ${1} ${2} ${3} ${4} ${5} ${6} ${7} < runRESOLV3D.R > runRESOLV3D.log 2>runRESOLV3D.err 
	echo "===> Finished"
	echo `date`
fi
