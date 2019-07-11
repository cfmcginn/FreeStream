//Original Author: Megan Byres, University Colorado Boulder
//Modification by Chris McGinn

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

int InitED(const int nLattice, const double widthFromZeroFM, const double sigma)
{
  // some constants
  const double e0 = 0.00150022 * TMath::Power (140.*widthFromZeroFM/(5.*nLattice), 4);
  const double SCAL = 1. / 0.3989423 / 0.00630285;
  // create histogram
  TH2D* h2 = NULL;
  h2 = new TH2D("h2","energy", nLattice, -widthFromZeroFM, widthFromZeroFM, nLattice, -widthFromZeroFM, widthFromZeroFM);
  h2->Sumw2();
  // create function
  TF2* gauss = NULL;
  gauss = new TF2("gauss", "xygaus(0)", -10, 10, -10, 10);
  gauss->SetParameters(1, 0, sigma, 0, sigma);
  std::ofstream txtFile; // file to put text values of function in     
  txtFile.open("input/inited.dat");

  const int maxNLattice = 1000; // Define this just to set array size limits
  if(nLattice > maxNLattice){
    std::cout << "Given nLattice \'" << nLattice << "\' is greater than maxNLattice \'" << maxNLattice << "\'. return 1" << std::endl;
    return 1;
  }

  Double_t bins[maxNLattice+1];

  getLinBins(-widthFromZeroFM, widthFromZeroFM, nLattice, bins);


  std::cout << "Creating lattice " << nLattice << "x" << nLattice << ", from -" << (Int_t)widthFromZeroFM << " to " << (Int_t)widthFromZeroFM << " fm. Points at: " << std::endl;
  for(Int_t bI = 0; bI < nLattice; ++bI){
    std::cout << (bins[bI] + bins[bI+1])/2. << ", ";
  }
  std::cout << "fm in x and y. Gaussian width of " << prettyString(sigma, 2, false) << " fm." << std::endl;
  
  
  for(int i = 0; i < nLattice; ++i){
    double valX = (bins[i] + bins[i+1])/2.;
      
    for(int j = 0; j < nLattice; ++j){
      double valY = (bins[j] + bins[j+1])/2.;
      h2->SetBinContent(i, j, gauss->Eval(valX, valY)*e0*SCAL); //set to function at that point * both constants
    }
  }

  h2->Scale(1./h2->Integral()); // normalize it

  for(int i = 0; i < nLattice; ++i){
    for(int j = 0; j < nLattice; ++j){
      txtFile << prettyString(h2->GetBinContent(i, j), 10, false) << "\t";
    }
    //    txtFile << std::endl;
  }
  txtFile << std::endl;
  txtFile.close();

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "Usage: ./bin/InitED.exe <nLattice> <widthFromZeroFM> <sigmaFM>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += InitED(std::stoi(argv[1]), std::stod(argv[2]), std::stod(argv[3]));
  return retVal;  
}
