//Author Chris McGinn

//cpp
#include <iostream>
#include <vector>

//ROOT
#include "TFile.h"
#include "TH2D.h"

//Local
#include "include/checkMakeDir.h"
#include "include/returnRootFileContentsList.h"

int initIPGlasma(int time)
{
  const std::string ipGlasmaFile = "input/IPGlasma_flat_useNucleus0_grid1024_g2mu0.10_m0.15_run00000.root";
  
  if(!checkFile(ipGlasmaFile)){
    std::cout << "Input ipGlasmaFile \'" << ipGlasmaFile << "\' is not found. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(ipGlasmaFile.c_str(), "READ");
  std::vector<std::string> histList = returnRootFileContentsList(inFile_p, "TH2D");

  for(auto const & hist : histList){
    std::cout << hist << std::endl;  
  }
  
  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/initIPGlasma.exe <time>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += initIPGlasma(std::stoi(argv[1]));
  return retVal;
}
