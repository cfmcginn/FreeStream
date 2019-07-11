//cpp
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/returnFileList.h"
#include "include/stringUtil.h"
#include "include/plotUtilities.h"

int processED(std::string inDir)
{
  if(!checkDir(inDir)){
    std::cout << "Input directory \'" << inDir << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList = returnFileList(inDir, ".dat");
  std::vector<double> timeStamps;
  if(fileList.size() == 0){
    std::cout << "Input directory \'" << inDir << "\' contains no valid .dat files. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  //Lets filter out everything that isnt energy density
  unsigned int pos = 0;
  while(pos < fileList.size()){
    if(fileList[pos].find("inited") == std::string::npos) fileList.erase(fileList.begin()+pos);
    else ++pos;
  }
  
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

  pos = 0;
  for(auto const & file : filePreProcessVect){
    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetRightMargin(0.12);
    canv_p->SetBottomMargin(0.12);

    TH2D* hist_p = new TH2D("hist_h", ";x;y", nX, -1, 1, nX, -1, 1);
    
    for(unsigned int vI = 0; vI < file.size(); ++vI){
      hist_p->SetBinContent(vI%nX, vI/nX, file[vI]);
    }    

    hist_p->SetMaximum(maxVal);
    hist_p->SetMinimum(minVal);
    hist_p->DrawCopy("COLZ");

    gStyle->SetOptStat(0);
    
    std::string saveName = fileList[pos];
    while(saveName.find("/") != std::string::npos){saveName.replace(0, saveName.find("/")+1, "");}
    saveName = "pdfDir/" + dateStr + "/" + saveName;
    saveName.replace(saveName.find(".dat"), 4, "");
    while(saveName.find(".") != std::string::npos){saveName.replace(saveName.find("."), 1, "p");}
    saveName = saveName + "_" + dateStr + ".pdf";    
    quietSaveAs(canv_p, saveName);

    delete hist_p;
    delete canv_p;
    ++pos;
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/processED.exe <inDir>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += processED(argv[1]);
  return retVal;
}