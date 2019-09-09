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

int compareIPFS(std::string inCompFileName, std::string inFreeStreamName = "")
{
  if(!checkFile(inCompFileName) || inCompFileName.find(".root") == std::string::npos){
    std::cout << "Given inCompFileName \'" << inCompFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(inFreeStreamName.size() != 0){
    if(!checkFile(inFreeStreamName) || inFreeStreamName.find(".root") == std::string::npos){
      std::cout << "Given inFreeStreamName \'" << inFreeStreamName << "\' is invalid. return 1" << std::endl;
      return 1;
    }
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

  const unsigned int nTH2 = 2;
  //const unsigned int nTH2 = th2List.size();

  Double_t marginTop = 0.05;
  Double_t marginOther = 0.14;
  
  for(unsigned int tI = 0; tI < nTH2; ++tI){
    Double_t time = timeHist_p->GetBinContent(tI+1);

    TCanvas* canv_p = new TCanvas("canv_p", "", 1000*5./3., 1000);
    canv_p->SetTopMargin(0.00);
    canv_p->SetRightMargin(0.00);
    canv_p->SetLeftMargin(0.00);
    canv_p->SetBottomMargin(0.00);    

    Double_t canvHeight = canv_p->GetWh();
    Double_t canvWidth = canv_p->GetWw();

    Double_t topMarginPix = canvHeight*marginTop;
    Double_t bottomMarginPix = canvHeight*marginOther;

    Double_t leftMarginPix = canvWidth*marginOther*0.6;
    Double_t rightMarginPix = canvWidth*marginOther*0.6;

    std::cout << "HEIGHT and WIDTH: " << canvHeight << ", " << canvWidth << std::endl;   

    const Int_t nPads = 3;
    TPad* pads[nPads];
    for(Int_t pI = 0; pI < nPads; ++pI){
      pads[pI] = NULL;
    }
    
    pads[0] = new TPad("pad0", "", 0.0, 0.0, 0.6, 1.0);
    pads[1] = new TPad("pad1", "", 0.6, 0.6667, 0.8, 1.0);
    pads[2] = new TPad("pad2", "", 0.6, 0.3333, 0.8, 0.6667);

    canv_p->cd();
    pads[0]->Draw("SAME");
    pads[0]->cd();
    pads[0]->SetTopMargin(marginTop);
    pads[0]->SetRightMargin(marginOther);
    pads[0]->SetLeftMargin(marginOther);
    pads[0]->SetBottomMargin(marginOther);
    
    TH2D* hist_p = (TH2D*)inFile_p->Get(th2List[tI].c_str());
    hist_p->GetXaxis()->SetTitle("x (fm/c)");
    hist_p->GetYaxis()->SetTitle("y (fm/c)");
    centerTitles(hist_p);
    hist_p->DrawCopy("COLZ");

    line_p->DrawLine(minSubVal, minSubVal, minSubVal, maxSubVal);
    line_p->DrawLine(maxSubVal, minSubVal, maxSubVal, maxSubVal);

    line_p->DrawLine(minSubVal, maxSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, minSubVal, maxSubVal, minSubVal);


    label_p->DrawLatex(0.2, 0.96, ("t = " + prettyString(time, 4, false) + " fm/c, IPGlasma Evolution").c_str());
    
    Double_t minVal = hist_p->GetXaxis()->GetBinLowEdge(1);
    Double_t maxVal = hist_p->GetXaxis()->GetBinLowEdge(hist_p->GetXaxis()->GetNbins()+2);
    
    Double_t pixPerX = (canvWidth*0.6 - rightMarginPix - leftMarginPix)/(maxVal - minVal);
    Double_t pixPerY = (canvHeight - bottomMarginPix - topMarginPix)/(maxVal - minVal);

    Double_t xLow = leftMarginPix/(canvWidth) + pixPerX*(minSubVal - minVal)/(canvWidth);
    Double_t xHigh = leftMarginPix/(canvWidth) + pixPerX*(maxSubVal - minVal)/(canvWidth);

    Double_t yLow = bottomMarginPix/(canvHeight) + pixPerY*(minSubVal - minVal)/canvHeight;
    Double_t yHigh = bottomMarginPix/(canvHeight) + pixPerY*(maxSubVal - minVal)/canvHeight;

    
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
    
    if(tI == 0 && false){
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
    pads[1]->Draw("SAME");
    pads[1]->cd();

    pads[1]->SetTopMargin(marginTop/0.3333);
    pads[1]->SetRightMargin(marginOther);
    pads[1]->SetLeftMargin(marginOther);
    pads[1]->SetBottomMargin(marginOther);

    //    hist2_p->SetMaximum(hist_p->GetMaximum());
    hist2_p->DrawCopy("COLZ");
    label_p->DrawLatex(0.2, 0.88, "Zoom");
    
    line_p->DrawLine(minSubVal, minSubVal, minSubVal, maxSubVal);
    line_p->DrawLine(maxSubVal, minSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, maxSubVal, maxSubVal, maxSubVal);
    line_p->DrawLine(minSubVal, minSubVal, maxSubVal, minSubVal);

    /*
    line_p->DrawLineNDC(0, 0, 0, 1);
    line_p->DrawLineNDC(1, 0, 1, 1);

    line_p->DrawLineNDC(0, 1, 1, 1);
    line_p->DrawLineNDC(0, 0, 1, 0);
    */
    
    pixPerX = (canvWidth/2. - rightMarginPix - leftMarginPix)/(maxSubVal - minSubVal);
    pixPerY = (canvHeight - bottomMarginPix - topMarginPix)/(maxSubVal - minSubVal);

    Double_t xLow2 = 0.6 + marginOther*0.2;
    Double_t xHigh2 = 0.8 - marginOther*0.2;

    Double_t yLow2 = 0.6667 + marginOther*0.3333;
    Double_t yHigh2 = 1.0 - topMarginPix/(canvHeight);

    std::cout << "Xlow, ylow, xhigh, yhigh: " << xLow << ", " << yLow << ", " << xHigh << ", " << yHigh << std::endl;

    std::cout << "Xlow2, ylow2, xhigh2, yhigh2: " << xLow2 << ", " << yLow2 << ", " << xHigh2 << ", " << yHigh2 << std::endl;

    Double_t m = (yLow2 - yLow)/(xHigh2 - xHigh);
    Double_t b = yLow2 - m*xHigh2;
    Double_t yLow3 = m*xLow2 + b;

    m = (yHigh2 - yHigh)/(xHigh2 - xHigh);
    b = yHigh2 - m*xHigh2;
    Double_t yHigh3 = m*xLow2 + b;
    
    
    canv_p->cd();
    line_p->SetLineWidth(2);
    line_p->DrawLineNDC(xLow, yLow, xLow2, yLow2);
    line_p->DrawLineNDC(xLow, yHigh, xLow2, yHigh2);

    line_p->DrawLineNDC(xHigh, yLow, xLow2, yLow3);
    line_p->DrawLineNDC(xHigh, yHigh, xLow2, yHigh3);


    line_p->DrawLineNDC(xLow2, yLow3, xHigh2, yLow2);
    line_p->SetLineStyle(10);
    line_p->DrawLineNDC(xLow2, yHigh3, xHigh2, yHigh2);

    line_p->SetLineStyle(1);
    line_p->SetLineWidth(4);

    gStyle->SetOptStat(0);
    
    if(inFreeStreamName.size() != 0){
      TFile* inFile2_p = new TFile(inFreeStreamName.c_str(), "READ");
      TH2D* hist3_p = (TH2D*)inFile2_p->Get("inited_0p250_h");
      TH2D* hist4_p = new TH2D("hist4_h", ";x (fm/c);y (fm/c)", nBinsX, binsX, nBinsY, binsY);
      std::cout << hist4_p->GetNbinsX() << std::endl;
      
      centerTitles(hist4_p);
      canv_p->cd();
      pads[2]->Draw("SAME");
      pads[2]->cd();
      
      pads[2]->SetTopMargin(marginTop/0.3333);
      pads[2]->SetRightMargin(marginOther);
      pads[2]->SetLeftMargin(marginOther);
      pads[2]->SetBottomMargin(marginOther);
      

      macroHistToSubsetHist(hist3_p, hist4_p);

      hist4_p->DrawCopy("COLZ");
      delete hist4_p;

      label_p->DrawLatex(0.2, 0.88, "Freestream");

      
      inFile2_p->Close();
      delete inFile2_p;
    }

    std::string saveName = "pdfDir/" + dateStr + "/" + th2List[tI] + "_" + dateStr + ".png";
    quietSaveAs(canv_p, saveName);
    for(Int_t pI = 0; pI < nPads; ++pI){      
      if(pads[pI] != NULL) delete pads[pI];
    }

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
  if(argc < 2 || argc > 3){
    std::cout << "Usage: ./bin/compareIPFS.exe <inCompFileName> <inFreeStreamName=\"\" default>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 2) retVal += compareIPFS(argv[1]);
  else if(argc == 3) retVal += compareIPFS(argv[1], argv[2]);
  return retVal;
}
