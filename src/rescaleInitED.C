//Author: Chris McGinn

//cpp
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//Local
#include "include/checkMakeDir.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

//ROOT
#include "TDatime.h"
#include "TMath.h"

int rescaleInitED(std::string inParams, std::string inInitED, std::string inLUT)
{
  //Check files input are good
  bool goodFile = fileCheckWithMsg(inParams, ".txt", "1");
  goodFile = goodFile && fileCheckWithMsg(inInitED, ".dat", "1");
  goodFile = goodFile && fileCheckWithMsg(inLUT, ".txt", "1");
  if(!goodFile) return 1;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;
  
  //Lets get the temperature from the inParams file
  std::ifstream inFile(inParams.c_str());
  std::string tempStr = "";
  double temp = -1.0;
  
  while(std::getline(inFile, tempStr)){
    if(tempStr.find("TSTART") == std::string::npos) continue;

    tempStr.replace(0, tempStr.find("TSTART")+6, "");
    if(tempStr.size() == 0) continue;
    while(tempStr.substr(0, 1).find(" ") != std::string::npos){
      tempStr.replace(0, 1, "");
      if(tempStr.size() == 0) break;      
    }
    while(tempStr.substr(tempStr.size()-1, 1).find(" ") != std::string::npos){
      tempStr.replace(tempStr.size()-1, 1, "");
      if(tempStr.size() == 0) break;      
    }
    if(tempStr.size() == 0) continue;

    temp = std::stod(tempStr);
  }
  inFile.close();

  //Check temperature and either return or output for debug
  if(temp < 0){
    std::cout << "Couldn't find valid temperature, \'" << temp << "\' less than 0. return 1" << std::endl;
    return 1;
  }  
  std::cout << "Processing for temperature: " << prettyString(temp, 5, false) << std::endl;

  //Now lets build a pair of vectors from the LUT (Note: Lut is created by looking at central temperature reports from logs in processing point sources on an odd lattice)
  std::vector<double> temps;
  std::vector<double> energyDensities;
  inFile.open(inLUT.c_str());
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;
    if(tempStr.find(",") == std::string::npos) continue;

    temps.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
    tempStr.replace(0, tempStr.find(",")+1, "");
    energyDensities.push_back(std::stod(tempStr));
  }  
  inFile.close();

  //Do a LUT sort by temp
  for(int tI = 0; tI < (int)(temps.size()-1); ++tI){
    for(int tI2 = tI+1; tI2 < (int)temps.size(); ++tI2){
      if(temps[tI] > temps[tI2]){
	double tempTemp = temps[tI];
	double tempED = energyDensities[tI];

	temps[tI] = temps[tI2];
	energyDensities[tI] = energyDensities[tI2];

      	temps[tI2] = tempTemp;
	energyDensities[tI2] = tempED;
      }
    }    
  }

  //Test that the temperature exists within the lut range. If not, return with message. DO NOT EXTRAPOLATE DUMBDUMB YOU WILL MESS IT UP
  if(temp < temps[0] || temp > temps[temps.size()-1]){
    std::cout << "Temperature from \'" << inParams << "\', " << prettyString(temp, 5, false) << ", is outside range of validity from LUT, " << prettyString(temps[0], 5, false) << "-" << prettyString(temps[temps.size()-1], 5, false) << ". return 1" << std::endl;
    return 1;
  }

  //Lets find the closest temperature in the lut
  double delta = temps[temps.size()-1] - temps[0];
  int closestTempPos = -1;
  for(unsigned int tI = 0; tI < temps.size(); ++tI){
    if(TMath::Abs(temps[tI] - temp) < delta){
      delta = TMath::Abs(temps[tI] - temp);
      closestTempPos = tI;
    }
    //Think a break here would be safe but this code is fast so let it ride
  }
  
  //Ok, now we need to process the energy density file to rescale
  
  std::vector<double> initEDLatticeFlat;
  inFile.open(inInitED.c_str());
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue;

    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, ",");}
    while(tempStr.find("\t") != std::string::npos){tempStr.replace(tempStr.find("\t"), 1, ",");}

    std::vector<std::string> tempVect = strToVect(tempStr);
    std::vector<double> tempVect2;
    for(auto const & val : tempVect){
      tempVect2.push_back(std::stod(val));
    }
    initEDLatticeFlat.insert(initEDLatticeFlat.end(), tempVect2.begin(), tempVect2.end());
  }  
  inFile.close();

  //Lets just check that the lattice is square
  int size = TMath::Sqrt(initEDLatticeFlat.size());
  if(size*size != (int)initEDLatticeFlat.size()){
    std::cout << "Lattice size is not a square. return 1." << std::endl;
    return 1;
  }

  std::vector<double> initEDLatticeFlat2 = initEDLatticeFlat;
  std::sort(initEDLatticeFlat2.begin(), initEDLatticeFlat2.end());
  double median = initEDLatticeFlat2[initEDLatticeFlat2.size()/2 - 1];
  if(initEDLatticeFlat2.size()%2 == 0) median = (median + initEDLatticeFlat2[initEDLatticeFlat2.size()/2])/2.;
  double rescaleFactor = energyDensities[closestTempPos]/median;

  std::cout << "RESCALING BY " << energyDensities[closestTempPos] << "/" << median << "=" << rescaleFactor << std::endl;
  
  std::string outFileOrig = inInitED;
  outFileOrig.replace(outFileOrig.find(".dat"), 4, ("_PreRescale_" + dateStr + ".dat").c_str());
  std::ofstream outFile(outFileOrig.c_str());
  for(auto const & val : initEDLatticeFlat){
    outFile << val << "\t";
  }
  outFile.close();

  std::string outFileNew = inInitED;
  outFileNew.replace(outFileNew.find(".dat"), 4, ("_PostRescale_" + dateStr + ".dat").c_str());
  outFile.open(outFileNew.c_str());
  for(auto const & val : initEDLatticeFlat){
    outFile << val*rescaleFactor << "\t";
  }
  outFile.close();

  outFile.open(inInitED.c_str());
  for(auto const & val : initEDLatticeFlat){
    outFile << val*rescaleFactor << "\t";
  }
  outFile.close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/rescaleInitED.exe <inParams> <inInitED> <inLUT>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += rescaleInitED(argv[1], argv[2], argv[3]);
  return retVal;
}
