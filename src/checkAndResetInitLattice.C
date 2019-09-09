//Author: Chris McGinn

//cpp
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"

std::vector<std::string> vectFromInFile(std::string inFileName)
{
  std::vector<std::string> retVect = {};
  if(!checkFile(inFileName)) return retVect;

  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");}
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}

    if(tempStr.size() == 0) continue;
    
    std::vector<std::string> tempVect = strToVect(tempStr);
    retVect.insert(retVect.begin(), tempVect.begin(), tempVect.end());
  }
  inFile.close();
  
  return retVect;
}

int checkAndResetInitLattice(const int nLat, const std::string inInitEDName)
{
  const std::string mainStr = "inited";
  if(!fileCheckWithMsg(inInitEDName, mainStr + ".dat", "1")) return 1;
  std::vector<std::string> mainFullVect = vectFromInFile(inInitEDName);
  
  int fileNLat = std::sqrt(mainFullVect.size());
  if(fileNLat*fileNLat != mainFullVect.size()){
    std::cout << "Given file has size " << mainFullVect.size() << " with no integer square!!!! return 1" << std::endl;
    std::cout << " (int)sqrt(" << mainFullVect.size() << ")=" << fileNLat << "; " << fileNLat << "*" << fileNLat << "=" << fileNLat*fileNLat << "!=" << mainFullVect.size() << std::endl;
    return 1;
  }
  else if(fileNLat != nLat){
    std::cout << "Given nLattice=" << nLat << " is not equal to that foundin file \'" << inInitEDName << "\', " << fileNLat << ". return 1" << std::endl;
    return 1;
  }
  
  std::vector<std::string> otherInit = {"pi", "pixx", "pixy", "piyy", "ux", "uy"};
  std::vector<std::string> otherFileList;
  for(auto const & other : otherInit){
    std::cout << inInitEDName << ", " << mainStr << std::endl;
    otherFileList.push_back(inInitEDName);
    otherFileList[otherFileList.size()-1].replace(otherFileList[otherFileList.size()-1].find(mainStr), mainStr.size(), "init" + other);
    std::cout << otherFileList[otherFileList.size()-1] << std::endl;
  }


  std::cout << "Other files to check: " << std::endl;
  for(unsigned int fI = 0; fI < otherFileList.size(); ++fI){
    std::cout << fI << "/" << otherFileList.size() << ": " << otherFileList[fI] << std::endl;
  }
  
  for(unsigned int fI = 0; fI < otherFileList.size(); ++fI){
    if(!fileCheckWithMsg(otherFileList[fI], ".dat", "1")){
      std::cout << "WARNING: File \'" << otherFileList[fI] << "\' is missing. we will initialize it w/ " << nLat << "x" << nLat << " lattice of zeroes" << std::endl;

      std::ofstream outFile(otherFileList[fI].c_str());
      for(int bIX = 0; bIX < nLat; ++bIX){
	for(int bIY = 0; bIY < nLat; ++bIY){
	  outFile << "0\t";
	}
      }
      outFile.close();
    }
    else{
      std::vector<std::string> otherFullVect = vectFromInFile(otherFileList[fI]);
      if(otherFullVect.size() != mainFullVect.size()){
	std::cout << "Other initializing file \'" << otherFileList[fI] << "\' is size=" << otherFullVect.size() << ", not matching main size of " << mainFullVect.size() << ". Will be reset to matched size, all zero." << std::endl;

	std::ofstream outFile(otherFileList[fI].c_str());
	for(int bIX = 0; bIX < nLat; ++bIX){
	  for(int bIY = 0; bIY < nLat; ++bIY){
	    outFile << "0\t";
	  }
	}
	outFile.close();
      
      }      
    }
  }
  
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/checkAndResetInitLattice.exe <nLat> <inInitEDName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkAndResetInitLattice(std::stoi(argv[1]), argv[2]);
  return retVal;
}
