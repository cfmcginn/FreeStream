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
nCPU=$((nCPU - 2))
if [[ $nCPU -lt 1 ]]
then
    nCPU=1
fi

echo "TotalCPU: $totalNCPU, capping at $nCPU"
#Set this to keep computer from melting
export OMP_NUM_THREADS=$nCPU

nLattice=(10 20)

mkdir -p logdir

for i in "${nLattice[@]}"
do
    cp data/paramsTEMPLATE.txt data/params.txt
    sed -i -e "s@NLATTICE@$i@g" data/params.txt

    start=$SECONDS
    #DO THIS FIRST OR BUGS
    ./bin/initE.exe >& logdir/logInitE_N$i.log
    ./bin/FS.exe data/params.txt >& logdir/logFS_N$i.log
    duration=$(( SECONDS - start ))

    echo "Lattice N=$i took $duration seconds..."
done

echo "Running FS.exe complete!"
    
