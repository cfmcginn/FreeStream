//Author Chris McGinn
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
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

int compRootToDat(std::string inRootFile, std::string inDatFile, bool doPrint = false)
{
  if(!checkFile(inRootFile) || inRootFile.find(".root") == std::string::npos){
    std::cout << "Given inRootFile \'" << inRootFile << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  if(!checkFile(inDatFile) || inDatFile.find(".dat") == std::string::npos){
    std::cout << "Given inDatFile \'" << inDatFile << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inRootFile.c_str(), "READ");
  TH2D* rootHist_p = (TH2D*)inFile_p->Get("subsetHistWithBuffer_h");

  const Int_t nBinsRoot = rootHist_p->GetXaxis()->GetNbins();

  std::ifstream inFile(inDatFile.c_str());
  std::string tempStr;
  std::vector<std::string> datVect;
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;

    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}
    
    datVect = strToVect(tempStr);
  }  
  inFile.close();

  const Int_t nBinsDat = TMath::Sqrt(datVect.size());
  if(nBinsDat*nBinsDat != datVect.size()){
    std::cout << "datVect " << datVect.size() << " not sqrt. return 1" << std::endl;
    return 1;
  }
  if(nBinsDat != nBinsRoot){
    std::cout << "nBinsDat " << nBinsDat << " not equal to nBinsRoot " << nBinsRoot << ". return 1" << std::endl;
    return 1;
  }

  const Int_t nMaxBins = 2000;
  Int_t nBinsX = -1;
  Double_t binsX[nMaxBins];
  Int_t nBinsY = -1;
  Double_t binsY[nMaxBins];

  for(Int_t bIX = 0; bIX < rootHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    ++nBinsX;
    binsX[nBinsX] = rootHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
  }

  for(Int_t bIY = 0; bIY < rootHist_p->GetYaxis()->GetNbins()+1; ++bIY){
    ++nBinsY;
    binsY[nBinsY] = rootHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
  }

  TH2D* datHist_p = new TH2D("datHist_h", ";x (fm);y (fm)", nBinsX, binsX, nBinsY, binsY);

  Double_t pos = 0;
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    for(Int_t bIY = 0; bIY < nBinsY; ++bIY){
      Double_t val = std::stod(datVect[pos]);
      
      datHist_p->SetBinContent(bIX+1, bIY+1, val);
      datHist_p->SetBinError(bIX+1, bIY+1, 0.0);

      ++pos;
    }
  }

  setSumW2({rootHist_p, datHist_p});
  centerTitles({rootHist_p, datHist_p});

  if(doPrint){
    const double min = 0.001;
    
    for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
      for(Int_t bIY = 0; bIY < nBinsY; ++bIY){
	double rootContent = rootHist_p->GetBinContent(bIX+1, bIY);
	double datContent = datHist_p->GetBinContent(bIX+1, bIY);

	if(rootContent < min && datContent < min) continue;
	if(TMath::Abs(rootContent - datContent) < min) continue;

	std::cout << "Discrepancy: " << rootContent << ", " << datContent << std::endl;
      }
    }
  }
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 900, 900);
  canv_p->SetTopMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetRightMargin(0.14);
  canv_p->SetBottomMargin(0.12);
  
  datHist_p->Divide(rootHist_p);
  datHist_p->DrawCopy("COLZ");

  gStyle->SetOptStat(0);

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  std::string saveName = "pdfDir/" + dateStr + "/compRoot_" + dateStr + ".pdf";  
  quietSaveAs(canv_p, saveName);
  delete canv_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;  
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 4){
    std::cout << "Usage: ./bin/compRootToDat.exe <inRootFile> <inDatFile> <doPrint=0 default>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 3) retVal += compRootToDat(argv[1], argv[2]);
  else if(argc == 4) retVal += compRootToDat(argv[1], argv[2], std::stoi(argv[3]));
  return retVal;
}
