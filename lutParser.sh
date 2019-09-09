#!/bin/bash

outFileName=input/tempToEnergyDensityLUT.txt

DATE=`date +%Y%m%d`

temps=()
energyDensities=()

for i in logdir/LUT/logFS*Point*.log
do
#    echo $i
    temp=$(grep Step0 $i)
    temp=${temp#* }
    while [[ $temp == " "* ]]
    do
	temp=${temp# *}
    done
    temp=${temp#* }
    while [[ $temp == " "* ]]
    do
	temp=${temp# *}
    done

    temp=$(echo $temp)
    
    while [[ $temp == *" "* ]]
    do
	temp=${temp% *}
    done

    energyDensity=${i#*_PointMag}
    energyDensity=${energyDensity%.log}
    energyDensity=$(echo $energyDensity | sed -e "s@p@.@g")
    
#    echo $temp $energyDensity

    tempCount=$(echo $temp | wc -c)
    energyDensityCount=$(echo $energyDensity | wc -c)

    if [[ $tempCount -gt 0 ]]
    then
	if [[ $energyDensityCount -gt 0 ]]
	then
	    temps+=($temp)
	    energyDensities+=($energyDensity)
	fi
    fi
done

if [[ -f $outFileName ]]
then
    rm $outFileName
fi

pos=0
for i in "${temps[@]}"
do
    echo $i,"${energyDensities[$pos]}" >> $outFileName
    pos=$((pos + 1))
done
