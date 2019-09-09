//cpp
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStyle.h"

int getSquareFactor(int inVal)
{
  for(int i = 1; i < inVal; ++i){
    int square = i*i;
    if(square == inVal) return i;
  }

  std::cout << "Square factor not found, return -1" << std::endl;
  return -1;
}


std::vector<std::string> getVectFromStr(std::string inStr)
{
  inStr = inStr + ",";
  while(inStr.find(",,") != std::string::npos){
    inStr.replace(inStr.find(",,"), 2, ",");
  }  
  
  std::vector<std::string> retVect;
  while(inStr.find(",") != std::string::npos){
    retVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  return retVect;
}

int compDat(const std::string inFileName1, const std::string inFileName2)
{
  if(!checkFile(inFileName1) || inFileName1.find(".dat") == std::string::npos){
    std::cout << "inFileName1 \'" << inFileName1 << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  if(!checkFile(inFileName2) || inFileName2.find(".dat") == std::string::npos){
    std::cout << "inFileName2 \'" << inFileName2 << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  std::ifstream inFile(inFileName1.c_str());
  std::string tempStr;
  std::vector<std::string> fullVect1, fullVect2;
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;
    while(tempStr.find("  ") != std::string::npos){tempStr.replace(tempStr.find("  "), 2, " ");} 
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, " ");} 
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");} 
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");} 

    if(tempStr.size() == 0) continue;    

    if(fullVect1.size() != 0) std::cout << "WARNING: VECTOR ALREADY HAS CONTENTS - BE CAREFUL" << std::endl;

    std::vector<std::string> tempVect = getVectFromStr(tempStr);
    fullVect1.insert(fullVect1.end(), tempVect.begin(), tempVect.end());
  }
  inFile.close();
  inFile.open(inFileName2);
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;
    while(tempStr.find("  ") != std::string::npos){tempStr.replace(tempStr.find("  "), 2, " ");} 
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, " ");} 
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");} 
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");} 

    if(tempStr.size() == 0) continue;    

    if(fullVect2.size() != 0) std::cout << "WARNING: VECTOR ALREADY HAS CONTENTS - BE CAREFUL" << std::endl;

    std::vector<std::string> tempVect = getVectFromStr(tempStr);
    fullVect2.insert(fullVect2.end(), tempVect.begin(), tempVect.end());
  }
  inFile.close();

  //PostProcessing of the fullVectors
  for(unsigned int i = 0; i < fullVect1.size(); ++i){
    tempStr = fullVect1[i];

    if(tempStr.find("e") != std::string::npos){
      std::string tempStr2 = tempStr.substr(tempStr.find("e")+1, tempStr.size());
      if(std::stoi(tempStr2) <= -9) tempStr = "0";
    }    

    if(tempStr.find(".") != std::string::npos){
      std::string tempStr2 = "";
      if(tempStr.find("e") != std::string::npos){
	tempStr2 = tempStr.substr(tempStr.find("e")+1, tempStr.size());
	tempStr = tempStr.substr(0, tempStr.find("e"));
      }
      
      unsigned int posDot = tempStr.find(".");
      std::string lastChar = "";
      while(tempStr.size() - posDot > 3){
	lastChar = tempStr.substr(tempStr.size()-1, 1);
	tempStr = tempStr.substr(0, tempStr.size()-1);
      }
      if(lastChar.size() != 0){
	if(std::string("56789").find(lastChar) != std::string::npos){
	  lastChar = "";
	}
      }
      
      if(tempStr2.size() != 0) tempStr = tempStr + "e" + tempStr2;
    }

    fullVect1[i] = tempStr;
  }

  for(unsigned int i = 0; i < fullVect2.size(); ++i){
    tempStr = fullVect2[i];

    if(tempStr.find("e") != std::string::npos){
      std::string tempStr2 = tempStr.substr(tempStr.find("e")+1, tempStr.size());
      if(std::stoi(tempStr2) <= -9) tempStr = "0";
    }    

    if(tempStr.find(".") != std::string::npos){
      std::string tempStr2 = "";
      if(tempStr.find("e") != std::string::npos){
	tempStr2 = tempStr.substr(tempStr.find("e")+1, tempStr.size());
	tempStr = tempStr.substr(0, tempStr.find("e"));
      }
      
      unsigned int posDot = tempStr.find(".");
      std::string lastChar = "";
      while(tempStr.size() - posDot > 3){
	lastChar = tempStr.substr(tempStr.size()-1, 1);
	tempStr = tempStr.substr(0, tempStr.size()-1);
      }
      if(lastChar.size() != 0){
	if(std::string("56789").find(lastChar) != std::string::npos){
	  lastChar = "";
	}
      }
      
      if(tempStr2.size() != 0) tempStr = tempStr + "e" + tempStr2;
    }

    fullVect2[i] = tempStr;
  }


  unsigned int squareFactor1 = getSquareFactor(fullVect1.size());
  unsigned int squareFactor2 = getSquareFactor(fullVect2.size());
  unsigned int maxSize = 0;
  for(unsigned int i = 0; i < fullVect1.size(); ++i){
    if(fullVect1[i].size() > maxSize) maxSize = fullVect1[i].size();
  }
  for(unsigned int i = 0; i < fullVect2.size(); ++i){
    if(fullVect2[i].size() > maxSize) maxSize = fullVect2[i].size();
  }
  ++maxSize;
    
  std::cout << "Map of \'" << squareFactor1 << "x" << squareFactor1 << "\' in file \'" << inFileName1 << "\'." << std::endl;
  for(unsigned int i = 0; i < fullVect1.size(); ++i){
    tempStr = fullVect1[i];
    while(tempStr.size() < maxSize){tempStr = tempStr + " ";}
    std::cout << tempStr;
    if(i%squareFactor1 == squareFactor1-1){
      std::cout << std::endl;
      std::cout << std::endl;
    }
  }

  std::cout << "Map of \'" << squareFactor2 << "x" << squareFactor2 << "\' in file \'" << inFileName2 << "\'." << std::endl;
  for(unsigned int i = 0; i < fullVect2.size(); ++i){
    tempStr = fullVect2[i];
    while(tempStr.size() < maxSize){tempStr = tempStr + " ";}
    std::cout << tempStr;
    if(i%squareFactor2 == squareFactor2-1){
      std::cout << std::endl;
      std::cout << std::endl;
    }
  }

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);;

  TH2D* hist1_p = nullptr;
  TH2D* hist2_p = nullptr;
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.14);
  hist1_p = new TH2D("hist1_p", ";X;Y", squareFactor1, 0-0.5, squareFactor1-0.5, squareFactor1, 0-0.5, squareFactor1-0.5);
  centerTitles(hist1_p);
  
  for(unsigned int i = 0; i < fullVect1.size(); ++i){
    hist1_p->SetBinContent(i/squareFactor1 + 1, i%squareFactor1 + 1, std::stod(fullVect1[i]));
    hist1_p->SetBinError(i/squareFactor1 + 1, i%squareFactor1 + 1, 0.0);
  }

  hist1_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);

  std::string saveName = "pdfDir/" + dateStr + "/canv1_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.14);

  hist2_p = new TH2D("hist2_p", ";X;Y", squareFactor2, 0-0.5, squareFactor2-0.5, squareFactor2, 0-0.5, squareFactor2-0.5);
  centerTitles(hist2_p);

  for(unsigned int i = 0; i < fullVect2.size(); ++i){
    hist2_p->SetBinContent(i/squareFactor2 + 1, i%squareFactor2 + 1, std::stod(fullVect2[i]));
    hist2_p->SetBinError(i/squareFactor2 + 1, i%squareFactor2 + 1, 0.0);
  }

  hist2_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);

  saveName = "pdfDir/" + dateStr + "/canv2_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  std::vector<double> maxFactor; // exclude 0;

  for(Int_t bIX = 0; bIX < hist1_p->GetNbinsX(); ++bIX){
    for(Int_t bIY = 0; bIY < hist1_p->GetNbinsY(); ++bIY){

      Double_t val1 = hist1_p->GetBinContent(bIX+1, bIY+1);
      Double_t val2 = hist2_p->GetBinContent(bIX+1, bIY+1);

      if(val1 >= TMath::Power(10, -20)){
	if(val2 >= TMath::Power(10, -20)){

	  maxFactor.push_back(val1/val2);
	  
	}
      }
    }
  }

  std::sort(std::begin(maxFactor), std::end(maxFactor));

  std::cout << "Max ratio diff: " << maxFactor[0] << "-" << maxFactor[maxFactor.size()-1] << std::endl;

  canv_p = new TCanvas("canv_p", "", 450, 450);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.14);

  hist1_p->Divide(hist2_p);

  std::vector<double> vals;
  for(Int_t bIX = 0; bIX < hist1_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < hist1_p->GetYaxis()->GetNbins(); ++bIY){
      vals.push_back(hist1_p->GetBinContent(bIX+1, bIY+1));
    }
  }

  std::sort(std::begin(vals), std::end(vals));
  double normVal = vals[vals.size()/2];
  for(Int_t bIX = 0; bIX < hist1_p->GetXaxis()->GetNbins(); ++bIX){
    for(Int_t bIY = 0; bIY < hist1_p->GetYaxis()->GetNbins(); ++bIY){
      double val = hist1_p->GetBinContent(bIX+1, bIY+1)/normVal;
      hist1_p->SetBinContent(bIX+1, bIY+1, val);
    }
  }

  hist1_p->DrawCopy("COLZ");
  gStyle->SetOptStat(0);

  saveName = "pdfDir/" + dateStr + "/canvDiv_" + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  delete canv_p;

  
  delete hist1_p;
  delete hist2_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/compDat.exe <inFileName1> <inFileName2>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += compDat(argv[1], argv[2]);
  return retVal;
}
