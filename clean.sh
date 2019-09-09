#!/bin/bash

DATE=`date +%Y%m%d`

for i in pdfDir/$DATE/*.pdf
do
    if [[ -f $i ]]
    then
	rm $i
    fi
done

for i in output/*.dat
do
    if [[ -f $i ]]
    then
	rm $i
    fi
done

for i in /home/cfmcginn/CUBHIG/slides/$DATE*/figures/*.pdf
do
    if [[ -f $i ]]
    then
	rm $i
    fi
done
