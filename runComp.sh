#!/bin/bash

DATE=`date +%Y%m%d`

bigFile=/home/cfmcginn/Samples/IPGlasmaFiles20190814/IPGlasma_flat_useNucleus0_grid1024_g2mu0.10_m0.15_run00000.root

freeFile=/home/cfmcginn/CUBHIG/tempPlots2019/Aug16plots/processedEDHist_20190816.root

./bin/compareIPFS.exe $bigFile $freeFile


rm -f pdfDir/$DATE/h2_evt00000_t00001_$DATE.png
rm -f pdfDir/$DATE/h2_evt00000_t00036_$DATE.png

bash /home/cfmcginn/CUBHIG/slides/allDirToTex.sh pdfDir/$DATE/ "" 1 "Chris McGinn" 1
