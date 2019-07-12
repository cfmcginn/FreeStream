#!/bin.bash

#check we can correctly grab number of cpu
count=$(lscpu | grep "CPU(s):" | wc -l)

if [[ $count -le 0 ]]
then
    echo "lscpu appears to not return nCPU. Recommend running ./bin/FS.exe manually. exit 1"
    exit 1
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

totalNCPU=$nCPU
nCPU=$((nCPU - 4))
if [[ $nCPU -lt 1 ]]
then
    nCPU=1
fi

echo "TotalCPU: $totalNCPU, capping at $nCPU"
#Set this to keep computer from melting
export OMP_NUM_THREADS=$nCPU

nLattice=(80)

mkdir -p logdir

for i in "${nLattice[@]}"
do
    cp data/paramsTEMPLATE.txt data/params.txt
    sed -i -e "s@NLATTICE@$i@g" data/params.txt

    start=$SECONDS
    #DO THIS FIRST OR BUGS
    ./bin/initE.exe >& logdir/logInitE_N$i.log
    #OVERRIDE ENERGY DENSITY   
#    ./bin/InitED.exe $i 1.0 0.7
    ./bin/FS.exe data/params.txt >& logdir/logFS_N$i.log &

    count=$(grep "Done" logdir/logFS_N$i.log | wc -l)

    while [[ $count -le 0 ]]
    do
	sleep 5

	checkInf=$(grep "problem here"  logdir/logFS_N$i.log | wc -l)
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
	    	
        count=$(grep "Done" logdir/logFS_N$i.log | wc -l)
    done

    duration=$(( SECONDS - start ))    
    echo "Lattice N=$i took $duration seconds on $nCPU nCPU..."
done

wait
echo "Running FS.exe complete!"
    
