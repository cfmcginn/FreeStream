#!/bin/bash

startVal=0.0002094

for i in `seq 0 100`
do

    cp runTEMPLATE.sh run.sh
    sed -i -e "s@POINTMAGIN@$startVal@g" run.sh

    bash run.sh 1
    startVal=$(echo $startVal + .00001 | bc -l)
done

echo ""
echo "RUN2 COMPLETE!"
echo ""
