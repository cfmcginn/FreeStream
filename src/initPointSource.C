//Based on InitED, Original Author: Megan Byres, University Colorado Boulder
//initPointSource for a point source calc. by Chris McGinn

//cpp 
#include <fstream>
#include <string>

//ROOT
#include "TF2.h"
#include "TH2D.h"
#include "TMath.h"

//Local
#include "include/getLinBins.h"
#include "include/plotUtilities.h"

int InitED(const int nLattice)
{
  // some constants
  const double e0 = 1.4; //based on peak in gaus case

  std::ofstream txtFile; // file to put text values of function in     
  txtFile.open("input/inited.dat");

  const int maxNLattice = 1000; // Define this just to set array size limits
  if(nLattice > maxNLattice){
    std::cout << "Given nLattice \'" << nLattice << "\' is greater than maxNLattice \'" << maxNLattice << "\'. return 1" << std::endl;
    return 1;
  }

  std::cout << "Creating lattice " << nLattice << "x" << nLattice << "." << std::endl;
  std::cout << "Point source of " << prettyString(e0, 2, false) << " magnitude." << std::endl;
  
  for(int i = 0; i < nLattice; ++i){
    for(int j = 0; j < nLattice; ++j){
      if(i == nLattice/2 && j == nLattice/2) txtFile << prettyString(e0, 10, false) << "\t";
      else txtFile << "0.0" << "\t";
    }
    //    txtFile << std::endl;
  }
  txtFile << std::endl;
  txtFile.close();

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/initPointSource.exe <nLattice> . return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += InitED(std::stoi(argv[1]));
  return retVal;  
}
