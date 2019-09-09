//cpp
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"



int checkED(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".dat") == std::string::npos){
    std::cout << "inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fullVect;
    
  std::ifstream inFile(inFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");}
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}

    if(tempStr.size() == 0) continue;

    std::vector<std::string> tempVect = strToVect(tempStr);
    fullVect.insert(fullVect.end(), tempVect.begin(), tempVect.end());
  }
  inFile.close();

  int rootVal = std::sqrt(fullVect.size());
  std::cout << "sqrt(" << fullVect.size() << ")=" << rootVal << ", " << rootVal << "*" << rootVal << "=" << rootVal*rootVal << std::endl;

  for(unsigned int sI = 0; sI < fullVect.size(); ++sI){
    if(isStrSame(fullVect[sI], "0")){
      std::cout << "FOUND 0: " << sI << std::endl;
    }
  }
  
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkED.exe <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkED(argv[1]);
  return retVal;
}
