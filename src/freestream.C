//original author Jamie Nagle
//Modified for use w/ makefile Chris McGinn

//cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/plotUtilities.h"

int freestream(int itime = 2, int ievent = 0)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  if(itime == 1) std::cout << "Selected time evo = 0.1 fm/c" << std::endl;
  if(itime == 2) std::cout << "Selected time evo = 0.2 fm/c" << std::endl;  
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // author:     j.nagle
  // date:       created 07/02/2019
  // pseudocode:
  // read in t=0+ file and map each point to a circle of radius cT to calculate free streaming evolution
  // compare with t = T file

  char fname[100];
  const std::string ipGlasmaName = "input/IPGlasma_flat_useNucleus0_grid1024_g2mu0.10_m0.15_run00000.root";
  if(!checkFile(ipGlasmaName)){
    std::cout << "ipGlasmaFile \'" << ipGlasmaName << "\' is not found. return 1" << std::endl;
    return 1;
  }

  sprintf(fname, ipGlasmaName.c_str());
  TFile *fin = new TFile(fname, "READ");
  //  if (!fin) std::cout << "Could not open filename = " << fname << std::endl;

  char hname[100];
  sprintf(hname,"h2_evt0000%d_t00000",ievent);							     
  TH2D *h0 = static_cast <TH2D *> (fin->Get(hname)); // t = 1.0e-5
  TH2D *h2;
  sprintf(hname,"h2_evt0000%d_t00038",ievent);							       
  if (itime == 1) h2 = static_cast <TH2D *> (fin->Get(hname)); // t = 0.1 fm/c  
  sprintf(hname,"h2_evt0000%d_t00071",ievent);
  if (itime == 2) h2 = static_cast <TH2D *> (fin->Get(hname)); // t = 0.2 fm/c

  int nbins   = h0->GetNbinsX();
  double low  = h0->GetXaxis()->GetBinCenter(1) + h0->GetXaxis()->GetBinWidth(1)/2.0;
  double high = h0->GetXaxis()->GetBinCenter(nbins) + h0->GetXaxis()->GetBinWidth(nbins)/2.0;

  TH2D *hf = (TH2D *) h0->Clone();
  hf->Reset();
  // new TH2D("hf","hf",nbins,low,high,nbins,low,high);

  
  const static int NR = 1000; // number of points around the circle
  double dx[NR];
  double dy[NR];
  double ctau = 0.2;
  if (itime == 1) ctau = 0.1;
  if (itime == 2) ctau = 0.2;
  
  double dphi = 2.0*TMath::Pi() / ((double) NR);
  for (int ir=0;ir<NR;ir++) {
    dx[ir] = ctau * TMath::Sin(dphi * (double) ir);
    dy[ir] = ctau * TMath::Cos(dphi * (double) ir);
  }

  double zooml = 5.0;
  double zoomh = 10.0;
  // based on the zoom, figure out the range!
  int lowbin = h0->GetXaxis()->FindBin(zooml-2.0*ctau);
  int highbin = h0->GetXaxis()->FindBin(zoomh+2.0*ctau);

  if (lowbin  < 1) lowbin = 1;
  if (highbin > h0->GetNbinsX()) highbin = h0->GetNbinsX();
  
  // loop over all bins
  for (int ix=lowbin;ix<=highbin;ix++) {
    for (int iy=lowbin;iy<=highbin;iy++) {

      double weight = h0->GetBinContent(ix,iy);
      // check around a circle
      for (int ir=0;ir<NR;ir++) {
	int newx = hf->GetXaxis()->FindBin( h0->GetXaxis()->GetBinCenter(ix) + dx[ir] );
	int newy = hf->GetYaxis()->FindBin( h0->GetYaxis()->GetBinCenter(iy) + dy[ir] );
	hf->SetBinContent(newx,newy, hf->GetBinContent(newx, newy) + weight );
      }
    }
  }

  // integrate only over zoom range - without buffer zone of 2 * ctau
  lowbin = h0->GetXaxis()->FindBin(zooml);
  highbin = h0->GetXaxis()->FindBin(zoomh);

  hf->Scale(h2->Integral(lowbin,highbin,lowbin,highbin) / hf->Integral(lowbin,highbin,lowbin,highbin) );
  //  hf->Scale((1.0 / ctau) * (1.0 / (double) NR));
    
  h0->GetXaxis()->SetRangeUser(zooml,zoomh);
  h0->GetYaxis()->SetRangeUser(zooml,zoomh);  
  h2->GetXaxis()->SetRangeUser(zooml,zoomh);
  h2->GetYaxis()->SetRangeUser(zooml,zoomh);  
  hf->GetXaxis()->SetRangeUser(zooml,zoomh);
  hf->GetYaxis()->SetRangeUser(zooml,zoomh);  

  std::cout << h0->Integral() << std::endl;
  std::cout << h2->Integral() << std::endl;
  
  TCanvas *c = new TCanvas("c","c",10,10,2700,770);
  c->Divide(4,1);
  c->cd(1);
  gPad->SetRightMargin(0.15);
  h0->GetZaxis()->SetRangeUser(0.0,h0->GetMaximum());
  h0->DrawCopy("colz");
  c->cd(2);
  gPad->SetRightMargin(0.15);  
  h2->GetZaxis()->SetRangeUser(0.0,h2->GetMaximum());
  h2->DrawCopy("colz");
  c->cd(3);
  gPad->SetRightMargin(0.15);  
  hf->GetZaxis()->SetRangeUser(0.0,hf->GetMaximum());
  hf->DrawCopy("colz");
  c->cd(4);
  gPad->SetRightMargin(0.15);  
  TH2D *hr = (TH2D *) hf->Clone();
  hr->Reset();
  hr->Divide(hf,h2,1.0,1.0);
  //  hr->Add(hf,h2,1.0,-1.0);
  hr->DrawCopy("colz");

  TH2D *hmap = new TH2D("hmap","hmap",120,0.0,120.0,50,0.0,3.0);
  TProfile *hmapprof = new TProfile("hmapprof","hmapprof",120,0.0,120.0);  
  // ratio difference versus value 
  for (int ix=lowbin; ix<=highbin; ix++) {
    for (int iy=lowbin; iy<=highbin; iy++) {
      hmap->Fill(h2->GetBinContent(ix,iy) , hf->GetBinContent(ix,iy)/h2->GetBinContent(ix,iy));
      hmapprof->Fill(h2->GetBinContent(ix,iy) , hf->GetBinContent(ix,iy)/h2->GetBinContent(ix,iy));      
    }
  }
  hmap->DrawCopy("colz");
  hmapprof->SetMarkerStyle(20);
  hmapprof->DrawCopy("p,same");

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  std::string saveName = "pdfDir/" + dateStr + "/freestream_" + dateStr + ".pdf";
  
  quietSaveAs(c, saveName);
  
  delete hmap;
  delete c;
  
  fin->Close();
  delete fin;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 1 || argc > 3){
    std::cout << "Usage: ./bin/freestream.exe <itime=2 default> <ievent=0 default>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  if(argc == 1) retVal += freestream();
  else if(argc == 2) retVal += freestream(std::stoi(argv[1]));
  else if(argc == 3) retVal += freestream(std::stoi(argv[1]), std::stoi(argv[2]));
  return retVal;
}
