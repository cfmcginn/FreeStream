#!/bin/bash

if [[ $# -ge 2 ]]
then
    echo "USAGE: bash run.sh <OVERRIDE QUERY-default 0>. exit 1"
    exit 1
fi

doOverride=0
if [[ $# -eq 1 ]]
then
    if [[ $1 -eq 0 ]]
    then
	doOverride=$1
    elif [[ $1 -eq 1 ]]
    then
	doOverride=$1
    else
	echo "Given argument $1 must be 0 or 1. exit 1"
	exit 1
    fi	
fi

#check we can correctly grab number of cpu
count=$(lscpu | grep "CPU(s):" | wc -l)

if [[ $count -le 0 ]]
then
    echo "lscpu appears to not return nCPU. Recommend running ./bin/FS.exe manually. exit 1"
    exit 1
fi

FSCount=$(ps | grep FS.exe | wc -l)
if [[ $FSCount -ge 1 ]]
then
    echo "YOU HAVE RUNNING FS INSTANCES. PLEASE TERMINATE, THEN RUN. exit 1"
    exit 1
fi

make

DATE=`date +%Y%m%d`

if [[ $doOverride -eq 0 ]]
then
    goodInput=0
    while [[ $goodInput -eq 0 ]]
    do
	read -n1 -r -p "Delete output dir? [y/n]..." inVal
	if [[ $inVal == "y" ]]
	then
	    echo " Deleting output"	
	    rm -f output/*.dat
	    goodInput=1
	elif [[ $inVal == "n" ]]
	then
	    echo " Not deleting output"
	    goodInput=1
	else
	    echo " Please answer y/n"
	fi	
    done
fi

#determine number of cpus
nCPU=$(lscpu | grep -n "CPU(s):" | head -n 1)
nCPU=${nCPU#*"CPU(s):"}

while [[ $nCPU == " "* ]]
do
    nCPU=${nCPU#" "}
done
while [[ $nCPU == *" "* ]]
do
    nCPU=${nCPU%" "*}
done

coresToSave=2
totalNCPU=$nCPU
nCPU=$((nCPU - $coresToSave))
if [[ $nCPU -lt 1 ]]
then
    nCPU=1
fi

echo "TotalCPU: $totalNCPU, capping at $nCPU"
#Set this to keep computer from melting
export OMP_NUM_THREADS=$nCPU

nLattice=(40)
spacing=0.1487110644
#spacing=0.5076

pointMag=.0012094
pointMagStr=$pointMag

#Lets initialize input energy density
doPointSrc=0
doIPGlasma=1
doDefault=0
boolArr=($doPointSrc $doIPGlasma $doDefault)
boolTot=0
for i in "${boolArr[@]}"
do
    if [[ $i -eq 1 ]]
    then
	dummy=0
    elif [[ $i -eq 0 ]]
    then
	dummy=0
    else
	echo "Please select 0 or 1 for boolean inputs. exit 1"
	exit 1
    fi
	
    boolTot=$((boolTot + $i))
done

if [[ $boolTot -le 0 ]]
then
    echo "No boolean initialized to 1. exit 1"
    exit 1
elif [[ $boolTot -ge 2 ]]
then
    echo "Multiple boolean initialized to 1. exit 1"
    exit 1
fi
     
if [[ $pointMagStr == *"."* ]]
then
    pointMagStr=$(echo $pointMagStr | sed -e "s@\.@p@g")

    if [[ $pointMagStr == "p"* ]]
    then
	pointMagStr=0$pointMagStr
    fi
fi

mkdir -p logdir/$DATE

for i in "${nLattice[@]}"
do
    cp data/paramsTEMPLATE.txt data/params.txt
    sed -i -e "s@NLATTICE@$i@g" data/params.txt
    sed -i -e "s@SPACINGGEVINV@$spacing@g" data/params.txt

    logName=N"$i"_
    if [[ $doPointSrc -eq 1 ]]
    then
	logName="$logName"PointMag$pointMagStr
    elif [[ $doIPGlasma -eq 1 ]]
    then
	logName="$logName"IPGlasma
    elif [[ $doDefault -eq 1 ]] 
    then
	logName="$logName"Default
    fi
	
    start=$SECONDS
    #DO THIS FIRST OR BUGS
    ./bin/initE.exe >& logdir/$DATE/logInitE_$logName.log

    #OVERRIDE ENERGY DENSITY WITH X
    #    ./bin/InitED.exe $i 1.0 0.7
    #OVERRIDE ENERGY DENSITY WITH POINT SOURCE 
    if [[ $doPointSrc -eq 1 ]]
    then
	./bin/initPointSource.exe $i $pointMag
    elif [[ $doIPGlasma -eq 1 ]]
    then
	./bin/initIPGlasma.exe 0 9 10 5 6 0.1
	./bin/rescaleInitED.exe data/params.txt input/inited.dat input/tempToEnergyDensityLUT.txt
	#	exit 1
    fi

    ./bin/checkAndResetInitLattice.exe $i input/inited.dat 

    val=$?
    if [[ $val -eq 1 ]]
    then
	exit 1
    fi
    
    #copy over the initial conditions for later comparison
    cp input/inited.dat output/inited--0.100_$logName.dat
    #Process & background
    ./bin/FS.exe data/params.txt $logName >& logdir/$DATE/logFS_$logName.log &

    #Following is for inf handling and early termination
    count=$(grep "Done" logdir/$DATE/logFS_$logName.log | wc -l)
    while [[ $count -le 0 ]]
    do
	sleep 5

	checkInf=$(grep "problem here"  logdir/$DATE/logFS_$logName.log | wc -l)
	if [[ $checkInf -gt 0 ]]
	then
	    echo "INF ERROR, BREAKING"

	    psStr=$(ps | grep FS)
	    psStr=${psStr%" pts"*}
	    while [[ $psStr == " " ]]
	    do
		psStr=${psStr# }
		psStr=${psStr% }
	    done
	    kill $psStr
	    duration=$(( SECONDS - start ))    
	    echo "Lattice N=$i took $duration seconds on $nCPU nCPU before inf failure..."
	    exit 1
	fi
	    	
        count=$(grep "Done" logdir/$DATE/logFS_$logName.log | wc -l)
	count2=$(grep "TERMINATING" logdir/$DATE/logFS_$logName.log | wc -l)
	count=$((count + $count2))
    done

    duration=$(( SECONDS - start ))    
    echo "Lattice N=$i took $duration seconds on $nCPU nCPU..."
done

wait
echo "Running FS.exe complete!"
    
