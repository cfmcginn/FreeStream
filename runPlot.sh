#!/bin/bash

dir=output

#.5076 GeV-1 to .1fm

spaceFactorNum=1000
spaceFactorDenom=5076

spaceGeVInv=0.5948442576


spaceGeVInv=$(echo "$spaceGeVInv * $spaceFactorNum" | bc -l)
echo $spaceGeVInv
spaceGeVInv=$(echo "$spaceGeVInv / $spaceFactorDenom" | bc -l)
echo $spaceGeVInv

DATE=`date +%Y%m%d`
./bin/processED.exe $dir $spaceGeVInv 0.0 0.0 1

rm -f pdfDir/$DATE/*.gif

bash /home/cfmcginn/CUBHIG/slides/allDirToTex.sh pdfDir/$DATE "" 1 "Chris McGinn" 1
