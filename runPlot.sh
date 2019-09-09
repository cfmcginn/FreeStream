#!/bin/bash

dir=output

#.5076 GeV-1 to .1fm

spaceFactorNum=1000
spaceFactorDenom=5076

#spaceGeVInv=0.5948442576
spaceGeVInv=0.1487110644


spaceGeVInv=$(echo "$spaceGeVInv * $spaceFactorNum" | bc -l)
echo $spaceGeVInv
spaceGeVInv=$(echo "$spaceGeVInv / $spaceFactorDenom" | bc -l)
echo $spaceGeVInv

DATE=`date +%Y%m%d`
#./bin/processED.exe $dir $spaceGeVInv 9.5 5.5 0
./bin/processED.exe $dir $spaceGeVInv 0 0 0

rm -f pdfDir/$DATE/*.gif

comboArr=()
comboSize=0

for i in pdfDir/$DATE/inited*.pdf
do
    if [[ $i == *"--"* ]]
    then
#	echo $i
	tempStr=${i#*_}
	tempStr=${tempStr%_*}
	comboArr+=($tempStr)
	comboSize=$((comboSize + 1))
    fi
done

if [[ $comboSize -ge 2 ]]
then
    for i in pdfDir/$DATE/inited*${comboArr[0]}*.pdf
    do
	tempStr=$(echo $i | sed -e "s@${comboArr[0]}@COMBO@g")
	tempStr2=$(echo $i | sed -e "s@${comboArr[0]}@${comboArr[1]}@g")
	pdfjam $i $tempStr2 --nup 2x1 --landscape --outfile $tempStr
    done
fi


#mkdir -p pdfTEMP
#for j in "${comboArr[@]}"
#do
#    for i in pdfDir/$DATE/inited*$j*.pdf
#    do
#	mv $i pdfTEMP
#    done
#done
   
bash /home/cfmcginn/CUBHIG/slides/allDirToTex.sh pdfDir/$DATE "" 1 "Chris McGinn" 1

#mv pdfTEMP/*.pdf pdfDir/$DATE
