//Author: Chris McGinn (2019.08.15)
//cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/macroHistToSubsetHist.h"
#include "include/plotUtilities.h"
#include "include/returnRootFileContentsList.h"
#include "include/stringUtil.h"

int compareIPFS(std::string inCompFileName)
{
  if(!checkFile(inCompFileName) || inCompFileName.find(".root") == std::string::npos){
    std::cout << "Given inCompFileName \'" << inCompFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }
  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* inFile_p = new TFile(inCompFileName.c_str(), "READ");
  TH1D* timeHist_p = (TH1D*)inFile_p->Get("h1_time");
  std::vector<std::string> th2List = returnRootFileContentsList(inFile_p, "TH2D");
  std::vector<int> th2Times;
  
  for(auto const & th2 : th2List){
    std::string timeStr = th2;
    while(timeStr.find("_t") != std::string::npos){timeStr.replace(0, timeStr.find("_t")+2, "");}
    th2Times.push_back(std::stoi(timeStr));
  }

  for(unsigned int tI = 0; tI < th2List.size()-1; ++tI){
    for(unsigned int tI2 = tI+1; tI2 < th2List.size(); ++tI2){
      if(th2Times[tI] > th2Times[tI2]){
	int tempTime = th2Times[tI];
	std::string tempName = th2List[tI];

	th2Times[tI] = th2Times[tI2];
	th2List[tI] = th2List[tI2];

	th2Times[tI] = tempTime;
	th2List[tI] = tempName;
      }
    }
  }

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  const Int_t nMaxBins = 2000;
  Int_t nBinsX = 0;
  Double_t binsX[nMaxBins];
  Int_t nBinsY = 0;
  Double_t binsY[nMaxBins];

  TLine* line_p = new TLine();
  line_p->SetLineColor(kRed);
  line_p->SetLineWidth(4);
  //  line_p->SetLineStyle(2);
  
  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(42);

  
  const Double_t minSubVal = 5.0;
  const Double_t maxSubVal = 10.0;

  for(unsigned int tI = 0; tI < th2List.size(); ++tI){
    Double_t time = timeHist_p->GetBinContent(tI+1);

    TCanvas* canv_p = new TCanvas("canv_p", "", 1500*2, 1500);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.01);
    canv_p->SetBottomMargin(0.01);

    canv_p->Divide(2, 1);
    
    canv_p->cd();
    canv_p->cd(1);

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    
    TH2D* hist_p = (TH2D*)inFile_p->Get(th2List[tI].c_str());
    hist_p->GetXaxis()->SetTitle("x (fm)");
    hist_p->GetYaxis()->SetTitle("y (fm)");
    centerTitles(hist_p);

    hist_p->DrawCopy("COLZ");    
    gStyle->SetOptStat(0);

    label_p->DrawLatex(0.2, 0.96, ("t = " + prettyString(time, 4, false) + " fm/c").c_str());
    
    line_p->DrawLine(minSubVal, minSubVal, minSubVal, maxSubVal);
    line_p->DrawLine(maxSubVal, minSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, maxSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, minSubVal, maxSubVal, minSubVal);

    nBinsX = -1;
    for(Int_t bIX = 0; bIX < hist_p->GetXaxis()->GetNbins()+1; ++bIX){
      if(hist_p->GetXaxis()->GetBinLowEdge(bIX+1) >= minSubVal && hist_p->GetXaxis()->GetBinLowEdge(bIX+1) < maxSubVal){
	++nBinsX;
	binsX[nBinsX] = hist_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }      
    }
    
    nBinsY = -1;
    for(Int_t bIY = 0; bIY < hist_p->GetYaxis()->GetNbins()+1; ++bIY){
      if(hist_p->GetYaxis()->GetBinLowEdge(bIY+1) >= minSubVal && hist_p->GetYaxis()->GetBinLowEdge(bIY+1) < maxSubVal){
	++nBinsY;
	binsY[nBinsY] = hist_p->GetYaxis()->GetBinLowEdge(bIY+1);
      }      
    }

    if(tI == 0){
      std::cout << "SUB BIN CHECK: " << std::endl;
      std::cout << "nBinsX=" << nBinsX << std::endl;
      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	std::cout << binsX[bIX] << ", ";      
      }
      std::cout << std::endl;
      std::cout << "nBinsY=" << nBinsY << std::endl;
      for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
	std::cout << binsY[bIY] << ", ";      
      }
      std::cout << std::endl;
    }

    TH2D* hist2_p = new TH2D("hist2_p", ";x (fm);y (fm)", nBinsX, binsX, nBinsY, binsY);
    centerTitles(hist2_p);
    macroHistToSubsetHist(hist_p, hist2_p);        
    
    if(tI == 0) std::cout << "nBinsX x nBinsY: " << hist_p->GetXaxis()->GetNbins() << " x " << hist_p->GetYaxis()->GetNbins() << std::endl;

    canv_p->cd();
    canv_p->cd(2);

    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    //    hist2_p->SetMaximum(hist_p->GetMaximum());
    hist2_p->DrawCopy("COLZ");

    line_p->DrawLine(minSubVal, minSubVal, minSubVal, maxSubVal);
    line_p->DrawLine(maxSubVal, minSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, maxSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, minSubVal, maxSubVal, minSubVal);

    gStyle->SetOptStat(0);
    
    std::string saveName = "pdfDir/" + dateStr + "/" + th2List[tI] + "_" + dateStr + ".png";
    quietSaveAs(canv_p, saveName);
    delete canv_p;

    delete hist2_p;
  }

  delete line_p;
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])  
{
  if(argc != 2){
    std::cout << "Usage: ./bin/compareIPFS.exe <inCompFileName>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += compareIPFS(argv[1]);
  return retVal;
}
