//Author Chris McGinn

//cpp
#include <fstream>
#include <iostream>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/macroHistToSubsetHist.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"

int initIPGlasma(const int timeInit = 0, const double xLow = 5, const double xHigh = 10, const double yLow = 5, const double yHigh = 10, const double zeroBuffer = 0.5)
{   
  if(xHigh <= xLow){
    std::cout << "xHigh " << xHigh << " is <= xLow " << xLow << ". return 1" << std::endl;
    return 1;
  }
  if(yHigh <= yLow){
    std::cout << "yHigh " << yHigh << " is <= yLow " << yLow << ". return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  const std::string ipGlasmaFile = "input/IPGlasma_flat_useNucleus0_grid1024_g2mu0.10_m0.15_run00000.root";
  
  if(!checkFile(ipGlasmaFile)){
    std::cout << "Input ipGlasmaFile \'" << ipGlasmaFile << "\' is not found. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(ipGlasmaFile.c_str(), "READ");
  std::vector<std::string> histList = returnRootFileContentsList(inFile_p, "TH2D");

  const std::string timeInitStr = std::to_string(timeInit);

  const std::string histStrBase = "h2_evt00000_t00000";
  std::string histStrInit = histStrBase.substr(0, histStrBase.size()-timeInitStr.size()) + timeInitStr;

  bool histFound = false;
  
  for(auto const & hist : histList){
    if(isStrSame(hist, histStrInit)){
      histFound = true;
      break;
    }
  }

  if(!histFound){
    std::cout << "Requested initial time hist \'" << histStrInit << "\' is not found in file \'" << ipGlasmaFile << "\'. return 1" << std::endl;
    return 1;
  }

  std::cout << "Found requested initial time \'" << histStrInit << "\'. Processing for x " << xLow << "-" << xHigh << "fm, y " << yLow << "-" << yHigh << " fm." << std::endl;

  //CreateBinning for our subhist
  const Int_t nBinsMax = 2000;
  Int_t nBinsX = -1;
  Double_t binsX[nBinsMax];
  Int_t nBinsY = -1;
  Double_t binsY[nBinsMax];

  Int_t xPosLow = -1;
  Int_t xPosHigh = -1;
  Int_t yPosLow = -1;
  Int_t yPosHigh = -1;

  TH2D* initHist_p = (TH2D*)inFile_p->Get(histStrInit.c_str());

  bool badBound = false;
  if(xLow < initHist_p->GetXaxis()->GetBinLowEdge(1)){
    std::cout << "XLow " << xLow << " is less than inhist low edge " << initHist_p->GetXaxis()->GetBinLowEdge(1) << ". return 1" << std::endl;
    badBound = true;
  }
  if(xHigh > initHist_p->GetXaxis()->GetBinLowEdge(initHist_p->GetXaxis()->GetNbins()+1)){
    std::cout << "XHigh " << xHigh << " is less than inhist high edge " << initHist_p->GetXaxis()->GetBinLowEdge(initHist_p->GetXaxis()->GetNbins()+1) << ". return 1" << std::endl;
    badBound = true;
  }  
  if(yLow < initHist_p->GetYaxis()->GetBinLowEdge(1)){
    std::cout << "YLow " << yLow << " is less than inhist low edge " << initHist_p->GetYaxis()->GetBinLowEdge(1) << ". return 1" << std::endl;
    badBound = true;
  }
  if(yHigh > initHist_p->GetYaxis()->GetBinLowEdge(initHist_p->GetYaxis()->GetNbins()+1)){
    std::cout << "YHigh " << yHigh << " is less than inhist high edge " << initHist_p->GetYaxis()->GetBinLowEdge(initHist_p->GetYaxis()->GetNbins()+1) << ". return 1" << std::endl;
    badBound = true;
  }

  if(badBound){
    inFile_p->Close();
    delete inFile_p;
    return 1;
  }

  Double_t minDelta = initHist_p->GetXaxis()->GetBinLowEdge(initHist_p->GetNbinsX()+1) - initHist_p->GetXaxis()->GetBinLowEdge(1);
  for(Int_t bIX = 0; bIX < initHist_p->GetXaxis()->GetNbins(); ++bIX){
    if(initHist_p->GetXaxis()->GetBinLowEdge(bIX+2) - initHist_p->GetXaxis()->GetBinLowEdge(bIX+1) < minDelta){
      minDelta = initHist_p->GetXaxis()->GetBinLowEdge(bIX+2) - initHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
    }
  }
  for(Int_t bIY = 0; bIY < initHist_p->GetYaxis()->GetNbins(); ++bIY){
    if(initHist_p->GetYaxis()->GetBinLowEdge(bIY+2) - initHist_p->GetYaxis()->GetBinLowEdge(bIY+1) < minDelta){
      minDelta = initHist_p->GetYaxis()->GetBinLowEdge(bIY+2) - initHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
    }
  }
  minDelta /= 2.;

  for(Int_t bIX = 0; bIX < initHist_p->GetXaxis()->GetNbins(); ++bIX){
    Double_t lowE = initHist_p->GetXaxis()->GetBinLowEdge(bIX+1);

    if(TMath::Abs(lowE - xLow) < minDelta) xPosLow = bIX;
    if(TMath::Abs(lowE - xHigh) < minDelta) xPosHigh = bIX;

    if(xPosLow >= 0 && xPosHigh >= 0) break;
  }

  for(Int_t bIY = 0; bIY < initHist_p->GetYaxis()->GetNbins(); ++bIY){
    Double_t lowE = initHist_p->GetYaxis()->GetBinLowEdge(bIY+1);

    if(TMath::Abs(lowE - yLow) < minDelta) yPosLow = bIY;
    if(TMath::Abs(lowE - yHigh) < minDelta) yPosHigh = bIY;

    if(yPosLow >= 0 && yPosHigh >= 0) break;
  }  
  
  std::cout << "xPosLow: " << xLow << ", " << xPosLow << std::endl;
  std::cout << "xPosHigh: " << xHigh << ", " << xPosHigh << std::endl;
  std::cout << "yPosLow: " << yLow << ", " << yPosLow << std::endl;
  std::cout << "yPosHigh: " << yHigh << ", " << yPosHigh << std::endl;

  for(Int_t bIX = xPosLow; bIX < xPosHigh+1; ++bIX){
    ++nBinsX;
    binsX[nBinsX] = initHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
  }
  for(Int_t bIY = yPosLow; bIY < yPosHigh+1; ++bIY){
    ++nBinsY;
    binsY[nBinsY] = initHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
  }

  Double_t totalED = 0.0;
  Double_t totalN = 0.0;
  for(Int_t bIX = 0; bIX < initHist_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < initHist_p->GetYaxis()->GetNbins(); ++bIY){
      totalED += initHist_p->GetBinContent(bIX+1, bIY+1);
      ++totalN;
    }    
  }

  Double_t aveED = totalED/totalN;
  std::cout << "AVERAGE ED: " << aveED << std::endl;
  const Double_t factorToED = 1;//targetAveED/aveED;

  for(Int_t bIX = 0; bIX < initHist_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < initHist_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t newVal = initHist_p->GetBinContent(bIX+1, bIY+1)*factorToED;
      Double_t newErr = initHist_p->GetBinError(bIX+1, bIY+1)*factorToED;

      initHist_p->SetBinContent(bIX+1, bIY+1, newVal);
      initHist_p->SetBinError(bIX+1, bIY+1, newErr);
    }
  }
      
  
  TH2D* subsetHist_p = new TH2D("subsetHist_h", ";x (fm);y (fm)", nBinsX, binsX, nBinsY, binsY);
  setSumW2(subsetHist_p);
  centerTitles(subsetHist_p);
  macroHistToSubsetHist(initHist_p, subsetHist_p, true);

  std::cout << "Precheck 1: " << subsetHist_p->GetBinContent(1, 1) << std::endl;
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.12);
  canv_p->SetLeftMargin(0.12);
  canv_p->SetBottomMargin(0.12);

  subsetHist_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  std::string saveName = "pdfDir/" + dateStr + "/initIPGlasma_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  std::cout << "CREATING nLattice: " << nBinsX << "x" << nBinsY << "..." << std::endl;
  Double_t edgeDelta = binsX[1] - binsX[0];
  Int_t nLatticeBuffer = zeroBuffer/edgeDelta;
  std::cout << " Additional buffer #pm " << nLatticeBuffer << "x" << edgeDelta << "=" << edgeDelta*nLatticeBuffer << " (Requested " << zeroBuffer << ")." << std::endl;

  Int_t nBinsXFull = -1;
  Double_t binsXFull[nBinsMax];
  Int_t nBinsYFull = -1;
  Double_t binsYFull[nBinsMax];

  for(Int_t bIX = 0; bIX < nLatticeBuffer*2 + nBinsX + 1; ++bIX){
    ++nBinsXFull;

    if(bIX < nLatticeBuffer) binsXFull[nBinsXFull] = binsX[0] - edgeDelta*(nLatticeBuffer - bIX);
    else if(bIX >= nBinsX+1+nLatticeBuffer) binsXFull[nBinsXFull] = binsXFull[nBinsXFull-1] + edgeDelta;
    else binsXFull[nBinsXFull] = binsX[bIX - nLatticeBuffer];   
  }

  std::cout << "BinsX: " << nBinsX << std::endl;
  for(Int_t bIX = 0; bIX < nBinsX; ++bIX){
    std::cout << binsX[bIX] << ", ";
  }
  std::cout << binsX[nBinsX] << "." << std::endl;
  std::cout << std::endl;

  std::cout << "BinsXFull: " << nBinsXFull << std::endl;
  for(Int_t bIX = 0; bIX < nBinsXFull; ++bIX){
    std::cout << binsXFull[bIX] << ", ";
  }
  std::cout << binsXFull[nBinsXFull] << "." << std::endl;
  std::cout << std::endl;

  for(Int_t bIY = 0; bIY < nLatticeBuffer*2 + nBinsY + 1; ++bIY){
    ++nBinsYFull;

    if(bIY < nLatticeBuffer) binsYFull[nBinsYFull] = binsY[0] - edgeDelta*(nLatticeBuffer - bIY);
    else if(bIY >= nBinsY+1+nLatticeBuffer) binsYFull[nBinsYFull] = binsYFull[nBinsYFull-1] + edgeDelta;
    else binsYFull[nBinsYFull] = binsY[bIY - nLatticeBuffer];   
  }

  TH2D* subsetHistWithBuffer_p = new TH2D("subsetHistWithBuffer_h", ";x (fm);y (fm)", nBinsXFull, binsXFull, nBinsYFull, binsYFull);
  setSumW2(subsetHistWithBuffer_p);
  centerTitles(subsetHistWithBuffer_p);

  for(Int_t bIX = 0; bIX < nBinsXFull; ++bIX){
    for(Int_t bIY = 0; bIY < nBinsYFull; ++bIY){
      if(bIX < nLatticeBuffer) subsetHistWithBuffer_p->SetBinContent(bIX+1, bIY+1, 0.0);
      else if(bIY < nLatticeBuffer) subsetHistWithBuffer_p->SetBinContent(bIX+1, bIY+1, 0.0);
      else if(bIX > nLatticeBuffer+nBinsX) subsetHistWithBuffer_p->SetBinContent(bIX+1, bIY+1, 0.0);
      else if(bIY > nLatticeBuffer+nBinsX) subsetHistWithBuffer_p->SetBinContent(bIX+1, bIY+1, 0.0);
      else subsetHistWithBuffer_p->SetBinContent(bIX+1, bIY+1, subsetHist_p->GetBinContent(bIX+1-nLatticeBuffer, bIY+1-nLatticeBuffer));

      subsetHistWithBuffer_p->SetBinError(bIX+1, bIY+1, 0.0);
    }    
  }

  std::cout << "Precheck 2: " << subsetHistWithBuffer_p->GetBinContent(1, 1) << std::endl;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.12);
  canv_p->SetLeftMargin(0.12);
  canv_p->SetBottomMargin(0.12);

  subsetHistWithBuffer_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);
  saveName = "pdfDir/" + dateStr + "/initIPGlasmaWithBuffer_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  std::ofstream outFile("input/inited.dat");

  for(Int_t bIX = 0; bIX < subsetHistWithBuffer_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < subsetHistWithBuffer_p->GetYaxis()->GetNbins(); ++bIY){
      outFile << subsetHistWithBuffer_p->GetBinContent(bIX+1, bIY+1) << "\t";
    }
  }
  
  outFile.close();

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  saveName = "output/" + dateStr + "/initIPGlasma_" + dateStr + ".root";
  TFile* outFile_p = new TFile(saveName.c_str(), "RECREATE");

  subsetHist_p->Write("", TObject::kOverwrite);
  subsetHistWithBuffer_p->Write("", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;
  
  delete subsetHistWithBuffer_p;
  delete subsetHist_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 1 || argc > 7){
    std::cout << "Usage: ./bin/initIPGlasma.exe <timeInit=0 default> <xLow=5 default> <xHigh=10 default> <yLow=5 default> <yHigh=10 default> <zeroBuffer=0.5 default>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 1) retVal += initIPGlasma();
  else if(argc == 2) retVal += initIPGlasma(std::stoi(argv[1]));
  else if(argc == 3) retVal += initIPGlasma(std::stoi(argv[1]), std::stod(argv[2]));
  else if(argc == 4) retVal += initIPGlasma(std::stoi(argv[1]), std::stod(argv[2]), std::stod(argv[3]));
  else if(argc == 5) retVal += initIPGlasma(std::stoi(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]));
  else if(argc == 6) retVal += initIPGlasma(std::stoi(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]));
  else if(argc == 7) retVal += initIPGlasma(std::stoi(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]));
  return retVal;
}
