//cpp
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"
#include "include/returnFileList.h"

std::vector<std::string> getVectFromStr(std::string inStr)
{
  inStr = inStr + ",";
  while(inStr.find(",,") != std::string::npos){inStr.replace(inStr.find(",,"), 2, ",");}
  if(inStr.size() == 1) return {};

  std::vector<std::string> tempVect;
  while(inStr.find(",") != std::string::npos){
    tempVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  
  return tempVect;
}

int checkInitFiles(std::string dirPath)
{
  if(!checkDir(dirPath)){
    std::cout << "Given dirPath \'" << dirPath << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<std::string> fileList = returnFileList(dirPath, "input/init");
  std::vector<int> sizes;
  unsigned int pos = 0;
  while(pos < fileList.size()){
    if(fileList[pos].find("ORIG") != std::string::npos) fileList.erase(fileList.begin()+pos);
    else ++pos;
  }
  
  std::cout << "Processing..." << std::endl;
  for(auto const & file : fileList){
    std::ifstream inFile(file.c_str());
    std::string tempStr;
    std::vector<std::string> fullVect;
    
    while(std::getline(inFile, tempStr)){
      if(tempStr.size() == 0) continue;

      while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");}
      while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}
      while(tempStr.find(",,") != std::string::npos){tempStr.replace(tempStr.find(",,"), 2, ",");}
      if(tempStr.size() == 1) continue;
      
      std::vector<std::string> tempVect = getVectFromStr(tempStr);
      fullVect.insert(fullVect.end(), tempVect.begin(), tempVect.end());
    }    
    
    sizes.push_back(fullVect.size());    
    inFile.close();
  }

  std::cout << "Size results: " << std::endl;
  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " " << fileList[fI] << ": " << sizes[fI] << std::endl;
  }
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkInitFiles.exe <dirPath>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += checkInitFiles(argv[1]);
  return retVal;
}
