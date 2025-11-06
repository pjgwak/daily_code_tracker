//
//  createMCWeightingFuncs.C
//  
//
//  Created by Javier Martin Blanco on 28/06/16.
//
//

#include <stdio.h>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include "FitFunctions.h"


using namespace std;
void createMCWeightingFuncs(const char* fileData, const char* fileMC, const char* partName="JPsi", const char* collName="PP", const char* rapName="006", const char* partType = "prompt", Int_t polOrder = 3, Int_t nRand=100);

void customiseLegend(TLegend& legend)
{
  legend.SetBorderSize(0);
  legend.SetTextSize(0.04);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(0);
  legend.SetFillStyle(1001);
}

void createAll(bool isPrompt=true){
  if (isPrompt)
  {
    createMCWeightingFuncs("Distributions_data_prompt/cPP_pt_rap0006_cent0100.root", "Distributions_mc_prompt/cPP_pt_rap0006_cent0100.root","JPsi","PP","006","prompt");
    createMCWeightingFuncs("Distributions_data_prompt/cPbPb_pt_rap0006_cent0100.root", "Distributions_mc_prompt/cPbPb_pt_rap0006_cent0100.root","JPsi","PbPb","006","prompt");
    
    createMCWeightingFuncs("Distributions_data_prompt/cPP_pt_rap0612_cent0100.root", "Distributions_mc_prompt/cPP_pt_rap0612_cent0100.root","JPsi","PP","0612","prompt");
    createMCWeightingFuncs("Distributions_data_prompt/cPbPb_pt_rap0612_cent0100.root", "Distributions_mc_prompt/cPbPb_pt_rap0612_cent0100.root","JPsi","PbPb","0612","prompt");
    
    createMCWeightingFuncs("Distributions_data_prompt/cPP_pt_rap1218_cent0100.root", "Distributions_mc_prompt/cPP_pt_rap1218_cent0100.root","JPsi","PP","1218","prompt");
    createMCWeightingFuncs("Distributions_data_prompt/cPbPb_pt_rap1218_cent0100.root", "Distributions_mc_prompt/cPbPb_pt_rap1218_cent0100.root","JPsi","PbPb","1218","prompt");
    
    createMCWeightingFuncs("Distributions_data_prompt/cPP_pt_rap1824_cent0100.root", "Distributions_mc_prompt/cPP_pt_rap1824_cent0100.root","JPsi","PP","1824","prompt");
    createMCWeightingFuncs("Distributions_data_prompt/cPbPb_pt_rap1824_cent0100.root", "Distributions_mc_prompt/cPbPb_pt_rap1824_cent0100.root","JPsi","PbPb","1824","prompt");
    
//    createMCWeightingFuncs("Distributions_data_prompt/cPP_pt_rap0024_cent0100.root", "Distributions_mc_prompt/cPP_pt_rap0024_cent0100.root","JPsi","PP","024","prompt");
//    createMCWeightingFuncs("Distributions_data_prompt/cPbPb_pt_rap0024_cent0100.root", "Distributions_mc_prompt/cPbPb_pt_rap0024_cent0100.root","JPsi","PbPb","024","prompt");
  }
  else
  {
    createMCWeightingFuncs("Distributions_data_nonprompt/cPP_pt_rap0006_cent0100.root", "Distributions_mc_nonprompt/cPP_pt_rap0006_cent0100.root","JPsi","PP","006","nonprompt");
    createMCWeightingFuncs("Distributions_data_nonprompt/cPbPb_pt_rap0006_cent0100.root", "Distributions_mc_nonprompt/cPbPb_pt_rap0006_cent0100.root","JPsi","PbPb","006","nonprompt");
    
    createMCWeightingFuncs("Distributions_data_nonprompt/cPP_pt_rap0612_cent0100.root", "Distributions_mc_nonprompt/cPP_pt_rap0612_cent0100.root","JPsi","PP","0612","nonprompt");
    createMCWeightingFuncs("Distributions_data_nonprompt/cPbPb_pt_rap0612_cent0100.root", "Distributions_mc_nonprompt/cPbPb_pt_rap0612_cent0100.root","JPsi","PbPb","0612","nonprompt");
    
    createMCWeightingFuncs("Distributions_data_nonprompt/cPP_pt_rap1218_cent0100.root", "Distributions_mc_nonprompt/cPP_pt_rap1218_cent0100.root","JPsi","PP","1218","nonprompt");
    createMCWeightingFuncs("Distributions_data_nonprompt/cPbPb_pt_rap1218_cent0100.root", "Distributions_mc_nonprompt/cPbPb_pt_rap1218_cent0100.root","JPsi","PbPb","1218","nonprompt");
    
    createMCWeightingFuncs("Distributions_data_nonprompt/cPP_pt_rap1824_cent0100.root", "Distributions_mc_nonprompt/cPP_pt_rap1824_cent0100.root","JPsi","PP","1824","nonprompt");
    createMCWeightingFuncs("Distributions_data_nonprompt/cPbPb_pt_rap1824_cent0100.root", "Distributions_mc_nonprompt/cPbPb_pt_rap1824_cent0100.root","JPsi","PbPb","1824","nonprompt");
    
//    createMCWeightingFuncs("Distributions_data_nonprompt/cPP_pt_rap0024_cent0100.root", "Distributions_mc_nonprompt/cPP_pt_rap0024_cent0100.root","JPsi","PP","024","nonprompt");
//    createMCWeightingFuncs("Distributions_data_nonprompt/cPbPb_pt_rap0024_cent0100.root", "Distributions_mc_nonprompt/cPbPb_pt_rap0024_cent0100.root","JPsi","PbPb","024","nonprompt");
  }
}


void createMCWeightingFuncs(const char* fileData, const char* fileMC, const char* partName, const char* collName, const char* rapName, const char* partType, Int_t polOrder, Int_t nRand)
// partName = "JPsi" or "Psi2S", collName = "PP" or "PbPb"
// rapName = "Mid" or "Fwd"
// nRand is the number of random data distributions that we want to generate from the original one
// polOrder is the order of the polinomial to be used to fit the DATA/MC distribution (supported 1 or 2)
{
  // Input parameters checks

  // Read the input files and extract the corresponding distributions
  TFile* fData = new TFile(fileData,"READ");
  TFile* fMC = new TFile(fileMC,"READ");
  if (!fData || !fMC)
  {
    cout << "[ERROR]: One of the two input files was not found" << endl;
    return;
  }
  
  TString hName("");
  if (!strcmp(collName,"PP")) hName = "PPDistribution";
  else hName = "PbPbDistribution";
//  TString hName(Form("hdistr_%s_%s_%s",partName,collName,rapName));
  
  TGraphErrors* histoData = static_cast<TGraphErrors*>(fData->FindObjectAny(hName.Data()));
  if (!histoData)
  {
    cout << "Data histo " << hName.Data() << " not found in file " << fileData << endl;
    return;
  }

  TGraphErrors* histoMC = static_cast<TGraphErrors*>(fMC->FindObjectAny(hName.Data()));
  if (!histoMC)
  {
    cout << "MC histo " << hName.Data() << " not found in file " << fileMC << endl;
    return;
  }
  
  // Retrieve histos for 0-2.4 for normalisation
  TString sDataFile(fileData);
  sDataFile.Remove(sDataFile.First("/")+1,sDataFile.Sizeof());
  TFile* fNormData = new TFile(Form("%s%s",sDataFile.Data(),hName.Contains("PP") ? "cPP_pt_rap0024_cent0100.root" : "cPbPb_pt_rap0024_cent0100.root"),"READ");
  TFile* fNormDataLowPt = new TFile(Form("%s%s",sDataFile.Data(),hName.Contains("PP") ? "cPP_pt_rap1824_cent0100.root" : "cPbPb_pt_rap1824_cent0100.root"),"READ");
  sDataFile.ReplaceAll("data","mc");
  TFile* fNormMC = new TFile(Form("%s%s",sDataFile.Data(),hName.Contains("PP") ? "cPP_pt_rap0024_cent0100.root" : "cPbPb_pt_rap0024_cent0100.root"),"READ");
  TFile* fNormMCLowPt = new TFile(Form("%s%s",sDataFile.Data(),hName.Contains("PP") ? "cPP_pt_rap1824_cent0100.root" : "cPbPb_pt_rap1824_cent0100.root"),"READ");
  if (!fNormData || !fNormMC)
  {
    cout << "[ERROR]: One of the two input normalisation files was not found" << endl;
    return;
  }
  if (!fNormDataLowPt || !fNormMCLowPt)
  {
    cout << "[ERROR]: One of the two input normalisation files for low pt was not found" << endl;
    return;
  }
  
  TGraphErrors* histoNormData = static_cast<TGraphErrors*>(fNormData->FindObjectAny(hName.Data()));
  if (!histoData)
  {
    cout << "Normalisation data histo " << hName.Data() << " not found in file " << fNormData->GetName() << endl;
    return;
  }
  TGraphErrors* histoNormDataLowPt = static_cast<TGraphErrors*>(fNormDataLowPt->FindObjectAny(hName.Data()));
  if (!histoData)
  {
    cout << "Normalisation data histo " << hName.Data() << " not found in file " << fNormDataLowPt->GetName() << endl;
    return;
  }
  
  TGraphErrors* histoNormMC = static_cast<TGraphErrors*>(fNormMC->FindObjectAny(hName.Data()));
  if (!histoMC)
  {
    cout << "Normalisation MC histo " << hName.Data() << " not found in file " << fNormMC->GetName() << endl;
    return;
  }
  TGraphErrors* histoNormMCLowPt = static_cast<TGraphErrors*>(fNormMCLowPt->FindObjectAny(hName.Data()));
  if (!histoMC)
  {
    cout << "Normalisation MC histo " << hName.Data() << " not found in file " << fNormMCLowPt->GetName() << endl;
    return;
  }
  
  // Calculate normalisation factors for data and MC
  double normMC(0.), normData(0.);
  for (int i = 0 ; i < histoNormData->GetN() ; i++)
  {
    double x, y;
    double xMC, yMC;
    histoNormData->GetPoint(i,x,y);
    normData += y;
    histoNormMC->GetPoint(i,xMC,yMC);
    normMC += yMC;
  }
  
  for (int i = 0 ; i < histoNormDataLowPt->GetN() ; i++) // Add also Njpsi at low pt to normalisations
  {
    double x, y;
    double xMC, yMC;
    histoNormDataLowPt->GetPoint(i,x,y);
    if (x<6.5)
    {
      normData += y;
      histoNormMCLowPt->GetPoint(i,xMC,yMC);
      normMC += yMC;
    }
    else break;
  }

//  cout << "Total NJpsi in DATA = " << normData << endl;
//  cout << "Total NJpsi in MC = " << normMC << endl;
  
  normData = 1./normData;
  normMC = 1./normMC;
  
  // Create output objects
  TObjArray* aFunctions = new TObjArray();
  aFunctions->SetOwner(kTRUE);
  TObjArray* aRatios = new TObjArray();
  aRatios->SetOwner(kTRUE);
  TObjArray* aDataHistos = new TObjArray();
  aDataHistos->SetOwner(kTRUE);
  
  Double_t* exLow = histoData->GetEXlow();
  Double_t* exHigh = histoData->GetEXhigh();
  Double_t* exyLow = histoData->GetEYlow();
  Double_t* eyHigh = histoData->GetEYhigh();
  Double_t* x = histoData->GetX();
  Double_t* y = histoData->GetY();

  Int_t nPoints = histoData->GetN();
  Double_t xAxis[nPoints+1];
  for (int i = 0; i<nPoints ; i++)
  {
    xAxis[i] = x[i]-exLow[i];
//    cout << i << " ; " << xAxis[i] << endl;
  }
  xAxis[nPoints] = x[nPoints-1]+exHigh[nPoints-1];
//  cout << nPoints << " ; " << xAxis[nPoints] << endl;
  
  TH1D* histoDataClone = new TH1D(Form("hdistrData_%s_%s_nominal",partName,collName),"",nPoints,xAxis);
//  histoDataClone->Sumw2(kTRUE);
  for (int i = 0; i<nPoints ; i++)
  {
    histoDataClone->SetBinContent(i+1,y[i]);
    histoDataClone->SetBinError(i+1,exyLow[i]);
  }

  // Normalise the initial DATA and MC distributions
//  Double_t norm = 1./histoDataClone->Integral();
//  histoDataClone->Scale(norm); // Normalise histo to 1
  for (int i = 0; i<nPoints ; i++) // Normalise histo to 1
  {
    histoDataClone->SetBinContent(i+1,histoDataClone->GetBinContent(i+1)*normData);
    histoDataClone->SetBinError(i+1,histoDataClone->GetBinError(i+1)*normData);
  }
  histoDataClone->SetLineColor(1);
  histoDataClone->SetMarkerColor(1);
  histoDataClone->SetStats(0);
  histoDataClone->GetXaxis()->SetTitle("p_{T}");
  histoDataClone->GetYaxis()->SetTitle("dN/dp_{T}");
  aDataHistos->Add(histoDataClone); // First object in the arrays will correspond to the nominal results
  
  
  Double_t* exLowMC = histoMC->GetEXlow();
  Double_t* exHighMC = histoMC->GetEXhigh();
  Double_t* exyLowMC = histoMC->GetEYlow();
  Double_t* eyHighMC = histoMC->GetEYhigh();
  Double_t* xMC = histoMC->GetX();
  Double_t* yMC = histoMC->GetY();
  
  Int_t nPointsMC = histoMC->GetN();
  Double_t xAxisMC[nPointsMC+1];
  for (int i = 0; i<nPointsMC ; i++)
  {
    xAxisMC[i] = xMC[i]-exLowMC[i];
//    cout << xAxisMC[i] << endl;
  }
  xAxisMC[nPointsMC] = xMC[nPointsMC-1]+exHighMC[nPointsMC-1];
//  cout << xAxisMC[nPointsMC] << endl;
  
  TH1D* histoMCClone = new TH1D(Form("hdistrMC_%s_%s_nominal",partName,collName),"",nPointsMC,xAxisMC);
//  histoMCClone->Sumw2(kTRUE);
  for (int i = 0; i<nPointsMC ; i++)
  {
    histoMCClone->SetBinContent(i+1,yMC[i]);
    histoMCClone->SetBinError(i+1,exyLowMC[i]);
  }

//  norm = 1./histoMCClone->Integral();
//  histoMCClone->Scale(norm); // Normalise histo to 1
  for (int i = 0; i<nPointsMC ; i++) // Normalise histo to 1
  {
    histoMCClone->SetBinContent(i+1,histoMCClone->GetBinContent(i+1)*normMC);
    histoMCClone->SetBinError(i+1,histoMCClone->GetBinError(i+1)*normMC);
  }
  histoMCClone->SetLineColor(4);
  histoMCClone->SetMarkerColor(4);
  histoMCClone->SetStats(0);
  histoMCClone->GetXaxis()->SetTitle("p_{T}");
  histoMCClone->GetYaxis()->SetTitle("dN/dp_{T}");
  
  // Loop on nRand to get the randomised data distributions
  TH1* histoDataClone2(0x0);
  for (Int_t n=0 ; n<nRand ; n++)
  {
    histoDataClone2 = static_cast<TH1D*>(histoDataClone->Clone(Form("hdistr_%s_%s_n%d",partName,collName,n)));
    histoDataClone2->Reset();
    
    Int_t nBins = histoDataClone2->GetNbinsX();
    for (Int_t i=1; i<=nBins ; i++)
    {
      Double_t newYield = histoDataClone->GetBinContent(i);
      Double_t newError = histoDataClone->GetBinError(i);
      
//      if (i != nBins) // This is done to not randomise the last pt bin, which has already a big statistical uncertainty on the data/mc ratio
//      {
        TF1* fRand = new TF1("fRand","TMath::Gaus(x,[0],[1],1)",newYield-6.*newError,newYield+6.*newError); // Function to randomise yields in each bin
        fRand->SetNpx(10000);
        fRand->SetParameter(0,newYield);
        fRand->SetParameter(1,newError);
        
        newYield = fRand->GetRandom();
        while (newYield<=0.) {newYield = fRand->GetRandom();} // This introduces a bias, but avoid to have 0, so infinite data/mc ratios
//      }
      histoDataClone2->SetBinContent(i,newYield);
      histoDataClone2->SetBinError(i,newError);
      
      histoDataClone2->SetLineColor(2);
      histoDataClone2->SetMarkerColor(2);
      
    }
    
//    norm = 1./histoDataClone->Integral();
//    for (Int_t i=1; i<=nBins ; i++) // Normalise histo to 1
//    {
//      histoMCClone->SetBinContent(i,histoMCClone->GetBinContent(i)*normMC);
//      histoMCClone->SetBinError(i,histoMCClone->GetBinError(i)*normMC);
//    }
//    histoDataClone->Scale(norm); // Normalise randomised histo to 1
    aDataHistos->Add(histoDataClone2);
  }
  
  // Loop on randomised histos, ratio DATA/MC and fitting function
  TIter nextHisto(aDataHistos);
  histoDataClone = 0x0;
  histoDataClone2 = 0x0;
  Int_t nh(0);
  while ( (histoDataClone = static_cast<TH1D*>(nextHisto())) )
  {
    histoDataClone2 = static_cast<TH1D*>(histoDataClone->Clone(nh == 0 ? Form("hratio_%s_%s_nominal",partName,collName) : Form("hratio_%s_%s_n%d",partName,collName,nh-1)));
    Int_t nBins = histoDataClone2->GetNbinsX();
    for (Int_t i=1; i<=nBins ; i++) // Divide Data/MC
    {
      double dataY = histoDataClone2->GetBinContent(i);
      histoDataClone2->SetBinContent(i,dataY/histoMCClone->GetBinContent(i));
      double errorY = histoDataClone2->GetBinContent(i)*TMath::Sqrt(TMath::Power(histoDataClone2->GetBinError(i)/dataY,2.));
      histoDataClone2->SetBinError(i,errorY);
      if (dataY-errorY<=0) {
        cout << "dataY:errorY " << dataY << " " << errorY << endl;
        errorY = dataY;
      }
      histoDataClone2->SetBinError(i,errorY);
    }
//    histoDataClone2->Divide(histoMCClone);
    histoDataClone2->SetLineColor(nh==0 ? 1 : 2);
    histoDataClone2->SetMarkerColor(nh==0 ? 1 : 2);
    histoDataClone2->GetXaxis()->SetTitle("p_{T}");
    histoDataClone2->GetYaxis()->SetTitle("dN/dp_{T}(DATA)/dN/dp_{T}(MC)");
    histoDataClone2->GetYaxis()->SetRangeUser(0.0, 2.5);
    
    TString fName(nh==0 ? "wFunction_nominal" : Form("wFunction_n%d",nh-1));
    TF1* fWeight(0x0);
    if (polOrder==1) fWeight= new TF1(fName.Data(),Pol1,6.5,50.,2);
    else if (polOrder==2) fWeight= new TF1(fName.Data(),Pol2,6.5,50.,3);
    else if (polOrder==3) {
      double xmin = (xAxis[0]<=6.5 ? 6.5 : 3);
      fWeight= new TF1(fName.Data(),"TMath::Erf(x/[0]+[1])+[2]",xmin,50.);

      if (!strcmp(partType,"prompt")) fWeight->SetParameters(-3.5, 4.0, 2.0);
      else fWeight->SetParameters(22, 0.3, 0.2);

      if (!strcmp(collName,"PbPb") && !strcmp(rapName,"1824")) {
        if (!strcmp(partType,"prompt")) histoDataClone2->GetYaxis()->SetRangeUser(0.0, 2.5);
        else histoDataClone2->GetYaxis()->SetRangeUser(0.0, 2.5);
        fWeight->SetParameters(-15, 0.1, 2.0);
      }
    }
    fWeight->SetLineColor(nh==0 ? kBlack : kBlue );
    int status = 0, fitcount = 0;
    do {
      status = histoDataClone2->Fit(fWeight, "N");
      Double_t *pars = static_cast<Double_t*>(fWeight->GetParameters());
      pars[0] = pars[0]+1.1*(fitcount%2==0 ? 1 : -1);
      pars[1] = pars[1]+1.3*(fitcount%2==0 ? 1 : -1);
      pars[2] = pars[2]+1.1*(fitcount%2==0 ? 1 : -1);
      fWeight->SetParameters(pars);
      fitcount++;
    } while (status!=0 && fitcount<5);

    if (status!=0) {
      cout << "FIT FAILED: " << fName << endl;
    }
    
    aFunctions->Add(fWeight);
    aRatios->Add(histoDataClone2);
    
    nh++;
  }
  
  // Create canvases with the randomised data and weight functions
  TCanvas* cData = new TCanvas("canv","canv",600,600);
  TLegend* lData = new TLegend(0.6002865,0.6919831,0.7679083,0.8565401);
  customiseLegend(*lData);
  nextHisto.Reset();
  nh = 0;
  while ( (histoDataClone = static_cast<TH1D*>(nextHisto())) )
  {
    if (nh == 0)
    {
      lData->AddEntry(histoDataClone,"Nominal data","l");
      lData->AddEntry(histoMCClone,"Nominal MC","l");
    }
    histoDataClone->DrawClone("same");
    nh++;
  }
  lData->Draw("same");
  nextHisto.Reset(); //Following lines are just to plot the nominal histo on top of the others
  histoDataClone = static_cast<TH1D*>(nextHisto());
  histoDataClone->DrawClone("same");
  histoMCClone->DrawClone("same"); // Draw also the MC distribution
  
  TCanvas* cFuncs = new TCanvas("canvf","canvf",600,600);
  TLegend* lFuncs = new TLegend(0.6002865,0.6919831,0.7679083,0.8565401);
  customiseLegend(*lFuncs);
  TIter nextFunc(aFunctions);
  TF1* fDummy(0x0);
  nh = 0;
  while ( (fDummy = static_cast<TF1*>(nextFunc())) )
  {
    if (nh == 0)
    {
      lFuncs->AddEntry(fDummy,"Nominal","l");
      TH1D* hDummy = static_cast<TH1D*>((aRatios->At(0))->Clone()); // Get ratio histo for the axes label
      Double_t ymax = hDummy->GetMaximum();
      Double_t ymin = hDummy->GetMinimum();
      hDummy->Reset("ICES");
      hDummy->GetYaxis()->SetRangeUser(ymin-0.3*ymin,ymax+0.3*ymax);
      hDummy->Draw();
    }
    else fDummy->DrawClone("same");
    nh++;
  }
  lFuncs->Draw("same");
  nextFunc.Reset(); //Following lines are just to plot the nominal function on top of the others
  fDummy = static_cast<TF1*>(nextFunc());
  fDummy->DrawClone("same");
  
  TCanvas* cRatios = new TCanvas("canvr","canvr",600,600);
  TLegend* lRatios = new TLegend(0.6002865,0.6919831,0.7679083,0.8565401);
  customiseLegend(*lRatios);
  TIter nextRatio(aRatios);
  nh = 0;
  while ( (histoDataClone = static_cast<TH1D*>(nextRatio())) )
  {
    if (nh == 0) lRatios->AddEntry(histoDataClone,"Nominal","l");
    histoDataClone->DrawClone("same");
    nh++;
  }
  lRatios->Draw("same");
  nextRatio.Reset(); //Following lines are just to plot the nominal histo on top of the others
  histoDataClone = static_cast<TH1D*>(nextRatio());
  histoDataClone->DrawClone("same");
  
  
  // Save arrays with results to file
  TFile* fSave = new TFile(Form("weights_%s_%s_%s_%s.root",partName,collName,rapName,partType),"RECREATE");
  cData->Write("dataDistributions",TObject::kOverwrite | TObject::kSingleKey);
  cFuncs->Write("weightingFuncs",TObject::kOverwrite | TObject::kSingleKey);
  cRatios->Write("ratiosDataMC",TObject::kOverwrite | TObject::kSingleKey);
  histoMCClone->Write(histoMCClone->GetName(), TObject::kOverwrite | TObject::kSingleKey);
  aDataHistos->Write("hRandData", TObject::kOverwrite | TObject::kSingleKey);
  aRatios->Write("DataOverMC", TObject::kOverwrite | TObject::kSingleKey);
  aFunctions->Write("wFunctions", TObject::kOverwrite | TObject::kSingleKey);
  
  fSave->Close();
  
  // Close files and delete pointers
  fData->Close(); delete fData;
  fMC->Close(); delete fMC;
}
