//Author Chris McGinn
//cpp
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/returnFileList.h"
#include "include/stringUtil.h"
#include "include/plotUtilities.h"

int processED(std::string inDir, double latticeSpacing=0.1/*in fm*/, double xOffset = 0.0/*in fm*/, double yOffset = 0.0/*in fm*/, bool doCentralLightCone = false, std::string sortStr = "inited-")
{
  if(!checkDir(inDir)){
    std::cout << "Input directory \'" << inDir << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::vector<std::string> fileList = returnFileList(inDir, ".dat");
  std::vector<double> timeStamps;
  if(fileList.size() == 0){
    std::cout << "Input directory \'" << inDir << "\' contains no valid .dat files. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;  

  //Lets filter out everything that isnt energy density
  unsigned int pos = 0;
  while(pos < fileList.size()){
    if(fileList[pos].find(sortStr) == std::string::npos) fileList.erase(fileList.begin()+pos);
    else ++pos;
  }

  std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Now extract time from each file name and sort
  for(auto const & file : fileList){
    std::string time = file;
    time.replace(0, time.find("-")+1, "");
    time.replace(time.find(".dat"), 4, "");
    timeStamps.push_back(std::stod(time));
  }

  unsigned int fI = 0;
  while(fI < fileList.size()-1){
    bool isGood = true;
    for(unsigned int fI2 = fI+1; fI2 < fileList.size(); ++fI2){
      if(timeStamps[fI] > timeStamps[fI2]){
	double tempTimeStamp = timeStamps[fI];
	std::string tempFileName = fileList[fI];

	timeStamps[fI] = timeStamps[fI2];
	fileList[fI] = fileList[fI2];

	timeStamps[fI2] = tempTimeStamp;
	fileList[fI2] = tempFileName;
	
	isGood = false;
      }      
    }

    if(isGood) ++fI;    
  }
  
  //Lets check our inputs
  std::cout << "Processing " << fileList.size() << " files..." << std::endl;
  for(auto const & file : fileList){
    std::cout << " " << file << std::endl;
  }
  std::cout << "Corresponding to times: " << std::endl;
  for(auto const & time : timeStamps){
    std::cout << " " << prettyString(time, 3, false) << " fm/c" << std::endl;
  }

  double baseTime = timeStamps[0];
  
  //Process the first file to get a handle on the binning
  Int_t total = 0;

  std::ifstream inFile(fileList[0].c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}    
    std::vector<std::string> tempVect = strToVect(tempStr);
    if(tempVect.size() == 0) continue;

    total = tempVect.size();
  }  

  Int_t nX = TMath::Sqrt(total);
  if(nX*nX != total){
    std::cout << "Total " << total << " is not equal to nX*nX " << nX << "*" << nX << ". return 1" << std::endl;
    return 1;
  }

  std::cout << "nX x nY: " << nX << "x" << nX << std::endl;
  inFile.close();

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  std::vector<std::vector<double> > filePreProcessVect;
  Double_t minVal = 1000000000000000;
  Double_t maxVal = -minVal;
  for(auto const & file : fileList){
    inFile.open(file.c_str());

    filePreProcessVect.push_back({});
    while(std::getline(inFile, tempStr)){
      while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}
      std::vector<std::string> tempVect = strToVect(tempStr);
      if(tempVect.size() == 0) continue;            
      
      for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
	Double_t val = std::stod(tempVect[vI]);
	filePreProcessVect[filePreProcessVect.size()-1].push_back(val);
	if(val > maxVal) maxVal = val;
	if(val < minVal) minVal = val;
      }
    }
    
    inFile.close();
  }

  Double_t lowHigh = (nX*latticeSpacing)/2;

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(18);
  
  pos = 0;
  for(auto const & file : filePreProcessVect){
    TCanvas* canv_p = new TCanvas("canv_p", "", 900, 900);
    canv_p->SetTopMargin(0.05);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetRightMargin(0.18);
    canv_p->SetBottomMargin(0.12);

    TH2D* hist_p = new TH2D("hist_h", ";x (fm);y (fm)", nX, -lowHigh+xOffset, lowHigh+xOffset, nX, -lowHigh+yOffset, lowHigh+yOffset);
    centerTitles(hist_p);

    TH2D* lightCone_p = new TH2D("lightCone_p", "", 10*nX, -lowHigh+xOffset, lowHigh+xOffset, 10*nX, -lowHigh+yOffset, lowHigh+yOffset);
    
    std::vector<double> contents;
    for(unsigned int vI = 0; vI < file.size(); ++vI){
      hist_p->SetBinContent((vI/nX)+1, (vI%nX)+1, file[vI]);

      contents.push_back(file[vI]);
    }

    std::sort(std::begin(contents), std::end(contents));
    
    Double_t minVal2 = contents[0]/2.;
    Int_t maxPos = TMath::Min(997*contents.size()/1000, contents.size()-1);
    Double_t maxVal2 = contents[maxPos];

    if(TMath::Abs(contents[0] - contents[maxPos]) < 0.000001){
      maxVal2 = contents[contents.size()-1] * 1.1;
    }
    
    hist_p->SetMaximum(maxVal);
    hist_p->SetMinimum(minVal);
    hist_p->SetMaximum(maxVal2);
    hist_p->SetMinimum(minVal2);
    hist_p->SetMaximum(0.0003);
    hist_p->SetMinimum(0.0);

    //IF DOING SURF DO ADDITIONAL OFFSETS
    std::string drawOpt = "COLZ";

    if(isStrSame(drawOpt, "SURF2")){
      hist_p->GetXaxis()->SetTitleOffset(2.0);
      hist_p->GetYaxis()->SetTitleOffset(2.0);
    }
    //    hist_p->DrawCopy("SURF2");
    //    gStyle->SetPalette(kRainBow);
    hist_p->DrawCopy(drawOpt.c_str());

    label_p->DrawLatex(0.25, .97, ("#bf{t = " + prettyString(timeStamps[pos], 5, false) + " fm/c}").c_str());

    if(doCentralLightCone){
      double radius = timeStamps[pos] - baseTime;

      if(radius > 0.001){ 
	for(Int_t bIX = 0; bIX < lightCone_p->GetXaxis()->GetNbins(); ++bIX){
	  Double_t xLow = lightCone_p->GetXaxis()->GetBinLowEdge(bIX+1);
	  Double_t xHigh = lightCone_p->GetXaxis()->GetBinLowEdge(bIX+2);
	  
	  for(Int_t bIY = 0; bIY < lightCone_p->GetYaxis()->GetNbins(); ++bIY){
	    Double_t yLow = lightCone_p->GetYaxis()->GetBinLowEdge(bIY+1);
	    Double_t yHigh = lightCone_p->GetYaxis()->GetBinLowEdge(bIY+2);

	    bool goodBin = false;

	    if(TMath::Abs(xOffset - xLow) < radius){
	      Double_t y1 = yOffset + TMath::Sqrt(radius*radius - (xLow - xOffset)*(xLow - xOffset));
	      Double_t y3 = yOffset - TMath::Sqrt(radius*radius - (xLow - xOffset)*(xLow - xOffset));

	      if(y1 > yLow && y1 < yHigh) goodBin = true;
	      if(y3 > yLow && y3 < yHigh) goodBin = true;
	    }
	    
	    if(TMath::Abs(xOffset - xHigh) < radius){
	      Double_t y2 = yOffset + TMath::Sqrt(radius*radius - (xHigh - xOffset)*(xHigh - xOffset));
	      Double_t y4 = yOffset - TMath::Sqrt(radius*radius - (xHigh - xOffset)*(xHigh - xOffset));

	      if(y2 > yLow && y2 < yHigh) goodBin = true;
	      if(y4 > yLow && y4 < yHigh) goodBin = true;
	    }

	    if(goodBin) lightCone_p->SetBinContent(bIX+1, bIY+1, hist_p->GetMaximum());
	    else lightCone_p->SetBinContent(bIX+1, bIY+1, 0);
	    lightCone_p->SetBinError(bIX+1, bIY+1, 0);
	  }
	}

	//	gStyle->SetPalette(kRust);
	lightCone_p->DrawCopy("SAME");

	std::cout << "INT: " << lightCone_p->Integral() << std::endl;
      }
    }
    
    gStyle->SetOptStat(0);
    
    std::string saveName = fileList[pos];
    while(saveName.find("/") != std::string::npos){saveName.replace(0, saveName.find("/")+1, "");}
    saveName = "pdfDir/" + dateStr + "/" + saveName;
    saveName.replace(saveName.find(".dat"), 4, "");
    while(saveName.find(".") != std::string::npos){saveName.replace(saveName.find("."), 1, "p");}
    saveName = saveName + "_" + dateStr + ".pdf";    
    quietSaveAs(canv_p, saveName);

    saveName.replace(saveName.find(".pdf"), 4, ".gif");    
    quietSaveAs(canv_p, saveName);

    delete hist_p;
    delete lightCone_p;
    delete canv_p;
    
    ++pos;
  }

  delete label_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 2 || argc > 6){
    std::cout << "Usage: ./bin/processED.exe <inDir> <latticeSpacing=0.1 default> <xOffset=0.0 default> <yOffset=0.0 default> <doCentralLightCone=false default>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += processED(argv[1]);
  else if(argc == 3) retVal += processED(argv[1], std::stod(argv[2]));
  else if(argc == 4) retVal += processED(argv[1], std::stod(argv[2]), std::stod(argv[3]));
  else if(argc == 5) retVal += processED(argv[1], std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]));
  else if(argc == 6) retVal += processED(argv[1], std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stoi(argv[5]));
  return retVal;
}
