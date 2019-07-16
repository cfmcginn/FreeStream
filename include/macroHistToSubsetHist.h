#ifndef MACROHISTTOSUBSETHIST_H
#define MACROHISTTOSUBSETHIST_H

#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"

bool macroHistToSubsetHist(TH1D* macroHist_p, TH1D* subsetHist_p, bool doSumW2 = false)
{
  std::vector<double> macroBins, subsetBins;

  for(Int_t bIX = 0; bIX < macroHist_p->GetNbinsX()+1; ++bIX){
    macroBins.push_back(macroHist_p->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIX = 0; bIX < subsetHist_p->GetNbinsX()+1; ++bIX){
    subsetBins.push_back(subsetHist_p->GetBinLowEdge(bIX+1));
  }

  //Check that subset bins Delta > 1
  bool allSubsetLargeDelta = true;
  for(unsigned int sI = 0; sI < subsetBins.size()-1; ++sI){
    if(subsetBins[sI+1] - subsetBins[sI] < .001){
      allSubsetLargeDelta = false;
      break;
    } 
  }

  if(!allSubsetLargeDelta){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have delta > .001" << std::endl;
    std::cout << " Subset hist: ";
    for(auto const & subsetVal : subsetBins){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  //Check that macrobins are well contained by subset bins
  bool allBinsGood = true;
  for(auto const & subVal : subsetBins){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBins){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }

  if(!allBinsGood){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have macro hist match" << std::endl;
    std::cout << " Macro name: " << macroHist_p->GetName() << std::endl;
    std::cout << " subset name: " << subsetHist_p->GetName() << std::endl;
    std::cout << " Macro hist: ";
    for(auto const & macroVal : macroBins){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset hist: ";
    for(auto const & subsetVal : subsetBins){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }
  else{
    for(unsigned int sI = 0; sI < subsetBins.size(); ++sI){
      Double_t val = 0;
      Double_t err = 0;
      
      for(Int_t bIX = 0; bIX < macroHist_p->GetNbinsX(); ++bIX){
	if(subsetBins[sI] <= macroHist_p->GetBinCenter(bIX+1) && macroHist_p->GetBinCenter(bIX+1) < subsetBins[sI+1]){
	  val += macroHist_p->GetBinContent(bIX+1);
	  err = TMath::Sqrt(err*err + macroHist_p->GetBinError(bIX+1)*macroHist_p->GetBinError(bIX+1));
	}
      }

      subsetHist_p->SetBinContent(sI+1, val);
      if(doSumW2) subsetHist_p->SetBinError(sI+1, err);
      else subsetHist_p->SetBinError(sI+1, TMath::Sqrt(val));
    }
  }
 
  return allBinsGood;
}


bool macroHistToSubsetHist(TH2D* macroHist_p, TH2D* subsetHist_p, bool doSumW2 = false, bool doSumLowY = false, bool doSumHighY = false)
{
  std::vector<double> macroBinsX, subsetBinsX, macroBinsY, subsetBinsY;

  for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    macroBinsX.push_back(macroHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIX = 0; bIX < subsetHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    subsetBinsX.push_back(subsetHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins()+1; ++bIY){
    macroBinsY.push_back(macroHist_p->GetYaxis()->GetBinLowEdge(bIY+1));
  }
  for(Int_t bIY = 0; bIY < subsetHist_p->GetYaxis()->GetNbins()+1; ++bIY){
    subsetBinsY.push_back(subsetHist_p->GetYaxis()->GetBinLowEdge(bIY+1));
  }

  //Check that subset bins Delta > 1
  bool allSubsetLargeDelta = true;
  for(unsigned int sI = 0; sI < subsetBinsX.size()-1; ++sI){
    if(subsetBinsX[sI+1] - subsetBinsX[sI] < .001){
      allSubsetLargeDelta = false;
      break;
    } 
  }
  for(unsigned int sI = 0; sI < subsetBinsY.size()-1; ++sI){
    if(subsetBinsY[sI+1] - subsetBinsY[sI] < .001){
      allSubsetLargeDelta = false;
      break;
    } 
  }

  if(!allSubsetLargeDelta){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have delta > .001" << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset histY: ";
    for(auto const & subsetVal : subsetBinsY){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  //Check that macrobins are well contained by subset bins
  bool allBinsGood = true;
  for(auto const & subVal : subsetBinsX){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBinsX){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }
  for(auto const & subVal : subsetBinsY){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBinsY){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }


  if(!allBinsGood){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have macro hist match" << std::endl;
    std::cout << " Macro name: " << macroHist_p->GetName() << std::endl;
    std::cout << " subset name: " << subsetHist_p->GetName() << std::endl;
    std::cout << " Macro histX: ";
    for(auto const & macroVal : macroBinsX){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Macro histY: ";
    for(auto const & macroVal : macroBinsY){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset histY: ";
    for(auto const & subsetVal : subsetBinsY){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }
  else{
    for(unsigned int sIX = 0; sIX < subsetBinsX.size(); ++sIX){
      for(unsigned int sIY = 0; sIY < subsetBinsY.size(); ++sIY){
	Double_t val = 0;
	Double_t err = 0;

	std::vector<Int_t> binsX, binsY;
      
	for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins(); ++bIX){
	  if(subsetBinsX[sIX] <= macroHist_p->GetXaxis()->GetBinCenter(bIX+1) && macroHist_p->GetXaxis()->GetBinCenter(bIX+1) < subsetBinsX[sIX+1]) binsX.push_back(bIX);
	}

	for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins(); ++bIY){
	  if(doSumLowY && sIY == 0){
	    if(macroHist_p->GetYaxis()->GetBinCenter(bIY+1) < subsetBinsY[sIY+1]) binsY.push_back(bIY);
	  }
	  else if(doSumHighY && sIY == subsetBinsY.size()-1){
	    if(subsetBinsY[sIY] <= macroHist_p->GetYaxis()->GetBinCenter(bIY+1)) binsY.push_back(bIY);
	  }
	  else{
	    if(subsetBinsY[sIY] <= macroHist_p->GetYaxis()->GetBinCenter(bIY+1) && macroHist_p->GetYaxis()->GetBinCenter(bIY+1) < subsetBinsY[sIY+1]) binsY.push_back(bIY);
	  }
	}


	for(unsigned int bIX = 0; bIX < binsX.size(); ++bIX){
	  for(unsigned int bIY = 0; bIY < binsY.size(); ++bIY){
	    val += macroHist_p->GetBinContent(binsX[bIX]+1, binsY[bIY]+1);
	    err = TMath::Sqrt(err*err + macroHist_p->GetBinError(binsX[bIX]+1, binsY[bIY]+1)*macroHist_p->GetBinError(binsX[bIX]+1, binsY[bIY]+1));
	  }
	}

	subsetHist_p->SetBinContent(sIX+1, sIY+1, val);
	if(doSumW2) subsetHist_p->SetBinError(sIX+1, sIY+1, err);
	else subsetHist_p->SetBinError(sIX+1, sIY+1, TMath::Sqrt(val));
      }
    }
  }
 
  return allBinsGood;
}


bool macroHistToSubsetHistX(TH2D* macroHist_p, TH1D* subsetHist_p, bool doSumW2 = false)
{
  std::vector<double> macroBinsX, subsetBinsX, macroBinsY;

  for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    macroBinsX.push_back(macroHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIX = 0; bIX < subsetHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    subsetBinsX.push_back(subsetHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins()+1; ++bIY){
    macroBinsY.push_back(macroHist_p->GetYaxis()->GetBinLowEdge(bIY+1));
  }

  //Check that subset bins Delta > 1
  bool allSubsetLargeDelta = true;
  for(unsigned int sI = 0; sI < subsetBinsX.size()-1; ++sI){
    if(subsetBinsX[sI+1] - subsetBinsX[sI] < .001){
      allSubsetLargeDelta = false;
      break;
    } 
  }

  if(!allSubsetLargeDelta){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have delta > .001" << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  //Check that macrobins are well contained by subset bins
  bool allBinsGood = true;
  for(auto const & subVal : subsetBinsX){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBinsX){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }

  if(!allBinsGood){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have macro hist match" << std::endl;
    std::cout << " Macro name: " << macroHist_p->GetName() << std::endl;
    std::cout << " subset name: " << subsetHist_p->GetName() << std::endl;
    std::cout << " Macro histX: ";
    for(auto const & macroVal : macroBinsX){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }
  else{
    for(unsigned int sIX = 0; sIX < subsetBinsX.size(); ++sIX){
      Double_t val = 0;
      Double_t err = 0;

      std::vector<Int_t> binsX;
      
      for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins(); ++bIX){
	if(subsetBinsX[sIX] <= macroHist_p->GetXaxis()->GetBinCenter(bIX+1) && macroHist_p->GetXaxis()->GetBinCenter(bIX+1) < subsetBinsX[sIX+1]) binsX.push_back(bIX);
      }
      
      for(unsigned int bIX = 0; bIX < binsX.size(); ++bIX){
	for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins(); ++bIY){
	  val += macroHist_p->GetBinContent(binsX[bIX]+1, bIY+1);
	  err = TMath::Sqrt(err*err + macroHist_p->GetBinError(binsX[bIX]+1, bIY+1)*macroHist_p->GetBinError(binsX[bIX]+1, bIY+1));
	}
      }

      subsetHist_p->SetBinContent(sIX+1, val);
      if(doSumW2) subsetHist_p->SetBinError(sIX+1, err);
      else subsetHist_p->SetBinError(sIX+1, TMath::Sqrt(val));      
    }    
  }
  
  return allBinsGood;
}



bool macroHistToSubsetHistY(TH2D* macroHist_p, TH1D* subsetHist_p, bool doSumW2 = false)
{
  std::vector<double> macroBinsX, subsetBinsX, macroBinsY;

  for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    macroBinsX.push_back(macroHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIX = 0; bIX < subsetHist_p->GetXaxis()->GetNbins()+1; ++bIX){
    subsetBinsX.push_back(subsetHist_p->GetXaxis()->GetBinLowEdge(bIX+1));
  }
  for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins()+1; ++bIY){
    macroBinsY.push_back(macroHist_p->GetYaxis()->GetBinLowEdge(bIY+1));
  }

  //Check that subset bins Delta > 1
  bool allSubsetLargeDelta = true;
  for(unsigned int sI = 0; sI < subsetBinsX.size()-1; ++sI){
    if(subsetBinsX[sI+1] - subsetBinsX[sI] < .001){
      allSubsetLargeDelta = false;
      break;
    } 
  }

  if(!allSubsetLargeDelta){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have delta > .001" << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  //Check that macrobins are well contained by subset bins
  bool allBinsGood = true;
  for(auto const & subVal : subsetBinsX){
    bool binHasMatch = false;
    for(auto const & macroVal : macroBinsY){
      if(TMath::Abs(macroVal - subVal) < 0.5){
	binHasMatch = true;
	break;
      }
    }

    if(!binHasMatch){
      allBinsGood = false;
      break;
    }
  }

  if(!allBinsGood){
    std::cout << "Error in macroHistToSubsetHist: Not all bin boundaries in subset hist have macro hist match" << std::endl;
    std::cout << " Macro name: " << macroHist_p->GetName() << std::endl;
    std::cout << " subset name: " << subsetHist_p->GetName() << std::endl;
    std::cout << " Macro histY: ";
    for(auto const & macroVal : macroBinsY){
      std::cout << macroVal << ",";
    }
    std::cout << std::endl;
    std::cout << " Subset histX: ";
    for(auto const & subsetVal : subsetBinsX){
      std::cout << subsetVal << ",";
    }
    std::cout << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }
  else{
    for(unsigned int sIX = 0; sIX < subsetBinsX.size(); ++sIX){
      Double_t val = 0;
      Double_t err = 0;

      std::vector<Int_t> binsY;
      
      for(Int_t bIY = 0; bIY < macroHist_p->GetYaxis()->GetNbins(); ++bIY){
	if(subsetBinsX[sIX] <= macroHist_p->GetYaxis()->GetBinCenter(bIY+1) && macroHist_p->GetYaxis()->GetBinCenter(bIY+1) < subsetBinsX[sIX+1]) binsY.push_back(bIY);
      }


      for(Int_t bIX = 0; bIX < macroHist_p->GetXaxis()->GetNbins(); ++bIX){
	for(unsigned int bIY = 0; bIY < binsY.size(); ++bIY){
	  val += macroHist_p->GetBinContent(bIX+1, binsY[bIY]+1);
	  err = TMath::Sqrt(err*err + macroHist_p->GetBinError(bIX+1, binsY[bIY]+1)*macroHist_p->GetBinError(bIX+1, binsY[bIY]+1));
	}
      }

      subsetHist_p->SetBinContent(sIX+1, val);
      if(doSumW2) subsetHist_p->SetBinError(sIX+1, err);
      else subsetHist_p->SetBinError(sIX+1, TMath::Sqrt(val));      
    }
  }
 
  return allBinsGood;
}


void test()
{
  TRandom3* randGen_p = new TRandom3(0);

  TH1D* test1_h = new TH1D("test1_h", "", 10, 0, 10);
  TH1D* test2_h = new TH1D("test2_h", "", 5, 0, 10);

  for(unsigned int i = 0; i < 500; ++i){
    Double_t val = randGen_p->Gaus(5, 1.0);

    test1_h->Fill(val);
  }

  macroHistToSubsetHist(test1_h, test2_h);


  test1_h->Print("ALL");
  test2_h->Print("ALL");

  delete test1_h;
  delete test2_h;

  delete randGen_p;
      
  return;
}

#endif
