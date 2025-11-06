#ifndef __CINT__
#endif
#include "initOniaTree.C"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TProfile.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TChain.h"
#include "TCut.h"
#include "TLegendEntry.h"
#include "TPaveStats.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooBinning.h"
#include "RooHist.h"

#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include "FitFunctions.h"

#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TString.h"
#include "signal.h"

void Plotting()
{
  TH1::AddDirectory(kFALSE);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetOptTitle(0); // at least most of the time
  gStyle->SetOptStat("nemr"); // most of the time, sometimes "nemriou" might be useful to display name, 
  //number of entries, mean, rms, integral, overflow and underflow
  gStyle->SetOptFit(0); // set to 1 only if you want to display fit results
  
  //=================================== Histo for all cut variables ==========================================
  using namespace RooFit;
  using namespace RooStats;
 
  void setWSRange(RooWorkspace *ws);
  void defineMassBkg(RooWorkspace *ws);
  void defineMassSig(RooWorkspace *ws);    

  double RError(double,double,double,double);
  double WeightRatioError(double,double,double,double);
  
  bool is3D =1;
  bool significance = 0;
  bool isHI=0; 
  bool Prompt = 0;
  bool nonPrompt = 1;

  bool savePlot=0;
  bool Binning =0;
  bool ctauFit = 1;
  

  float ptmin=3.0; float ptmax=30.0; float ymin=1.6; float ymax=2.4; bool absRapidity=true; 
  double threshold=0.8; double efficiency=0.9; int centmin=0; int centmax=200;

  ymin=1.6; ymax=2.4; ptmin=3.0; ptmax=30.0; 
  //ymin=1.6; ymax=2.4; ptmin=3.0; ptmax=6.5; 
  //ymin=1.6; ymax=2.4; ptmin=6.5; ptmax=12.0; 
  //ymin=1.6; ymax=2.4; ptmin=12.0; ptmax=30.0; 
  
  //ymin=0.0; ymax=1.6; ptmin=6.5; ptmax=30.0;  
  //ymin=0.0; ymax=1.6; ptmin=6.5; ptmax=9.0;  
  //ymin=0.0; ymax=1.6; ptmin=9.0; ptmax=12.0;  
  //ymin=0.0; ymax=1.6; ptmin=12.0; ptmax=15.0;  
  //ymin=0.0; ymax=1.6; ptmin=15.0; ptmax=20.0;  
  //ymin=0.0; ymax=1.6; ptmin=20.0; ptmax=30.0;  
  
  int range1 = 0; int fbin = 0;
  int range2 = 0; int lbin = 0;
  float xVal = 0.0;

  if(significance){
    range1 = -10;
    range2 = 50;
    fbin = -10; lbin = 50;
    xVal = 11.0;
  }else{
    range1 = -1;
    range2 = 3;
    fbin = -5; lbin = 10;
    xVal = 0.5;
  }

  float lCut=-10.0;
  TLegend *leg = new TLegend(0.38,0.5,0.87,0.75);
  if (isHI)
    leg->SetHeader("PYTHIA: PbPb #sqrt{s_{NN}} = 5.02 TeV");
  else
    leg->SetHeader("PYTHIA: pp #sqrt{s} = 5.02 TeV");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetMargin(0.15);
  leg->SetTextSize(0.032);
  leg->SetTextColor(2);

  TCut defaultCut = NULL;
  //  TCut defaultCut = "Reco_QQ_sign==0&&Gen_QQ_4mom.M()>3.09&&Gen_QQ_4mom.M()<3.10&&Reco_QQ_4mom.M()>2.9&&Reco_QQ_4mom.M()<3.2&&Reco_QQ_ctauTrue>-10";
  defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>2.9&&Reco_QQ_4mom.M()<3.2&&Reco_QQ_ctau>-500&&Reco_QQ_ctauTrue>-10";
  if(significance)
    defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>2.85&&Reco_QQ_4mom.M()<3.25&&Reco_QQ_ctau>-500&&Reco_QQ_ctauErr>-500&&Reco_QQ_ctauTrue>-10";
  if(is3D){
    defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>2.85&&Reco_QQ_4mom.M()<3.25&&Reco_QQ_ctau3D>-500&&Reco_QQ_ctauTrue3D>-10";
    if(significance)
      defaultCut = "Reco_QQ_sign==0&&Reco_QQ_4mom.M()>2.85&&Reco_QQ_4mom.M()<3.25&&Reco_QQ_ctau3D>-500&&Reco_QQ_ctauErr3D>-500&&Reco_QQ_ctauTrue3D>-10";
  }
  
  TCut ptCut = Form("Reco_QQ_4mom.Pt()>%3.1f&&Reco_QQ_4mom.Pt()<%3.1f",ptmin,ptmax);
  TCut rapCut;
  if (absRapidity)
    rapCut = Form("abs(Reco_QQ_4mom.Rapidity())>%3.1f&&abs(Reco_QQ_4mom.Rapidity())<%3.1f",ymin,ymax);
  else
    rapCut = Form("Reco_QQ_4mom.Rapidity()>%3.1f&&Reco_QQ_4mom.Rapidity()<%3.1f",ymin,ymax);
  TCut centCut;
  if (isHI)
    centCut = Form("Centrality>=%d&&Centrality<%d",centmin,centmax);
  
  // new Acceptance
  TCut IsPlAccept = "(abs(Reco_QQ_mupl_4mom->Eta()) < 2.4 && ( (abs(Reco_QQ_mupl_4mom->Eta()) < 1.2 && Reco_QQ_mupl_4mom->Pt() >= 3.5 ) || (1.2 <= abs(Reco_QQ_mupl_4mom->Eta()) && abs(Reco_QQ_mupl_4mom->Eta()) < 2.1 && Reco_QQ_mupl_4mom->Pt() >= 5.77-1.89*abs(Reco_QQ_mupl_4mom->Eta() )) || (2.1 <= abs(Reco_QQ_mupl_4mom->Eta()) && Reco_QQ_mupl_4mom->Pt() >= 1.8 )))";
  TCut IsMiAccept = "(abs(Reco_QQ_mumi_4mom->Eta()) < 2.4 && ( (abs(Reco_QQ_mumi_4mom->Eta()) < 1.2 && Reco_QQ_mumi_4mom->Pt() >= 3.5 ) || (1.2 <= abs(Reco_QQ_mumi_4mom->Eta()) && abs(Reco_QQ_mumi_4mom->Eta()) < 2.1 && Reco_QQ_mumi_4mom->Pt() >= 5.77-1.89*abs(Reco_QQ_mumi_4mom->Eta() )) || (2.1 <= abs(Reco_QQ_mumi_4mom->Eta()) && Reco_QQ_mumi_4mom->Pt() >= 1.8 )))";

  TCut MuPlID = "((Reco_QQ_mupl_isGoodMuon==1) && Reco_QQ_mupl_nTrkWMea> 5 &&  Reco_QQ_mupl_nPixWMea > 0 && TMath::Abs(Reco_QQ_mupl_dxy) < 0.3 && TMath::Abs(Reco_QQ_mupl_dz) < 20 && Reco_QQ_VtxProb > 0.01 )";
  TCut MuMiID = "((Reco_QQ_mumi_isGoodMuon==1) && Reco_QQ_mumi_nTrkWMea> 5 &&  Reco_QQ_mumi_nPixWMea > 0 && TMath::Abs(Reco_QQ_mumi_dxy) < 0.3 && TMath::Abs(Reco_QQ_mumi_dz) < 20 && Reco_QQ_VtxProb > 0.01 )";
  
  
  unsigned int trigBit;
  if (isHI)
    trigBit=1; // DoubleMu0_HighQ
  else
    trigBit=1; // DoubleMu0_HighQ
  TCut trigCut = Form("(HLTriggers&%u)==%u&&(Reco_QQ_trig&%u)==%u",trigBit,trigBit,trigBit,trigBit);

  TString fname;
  TString fnamepng;
  TString outfname;
  TString plotDir;
  
  if(!is3D){plotDir = "Plots/2015";}
  if(!is3D && significance){plotDir = "Plots/2015/significance";}
  if(is3D){plotDir = "Plots/2015/3D";}
  if(is3D && Binning){plotDir = "Plots/2015/3D/Binning";}
  if(is3D && ctauFit){plotDir = "Plots/2015/3D/ctauFit";}
  if(is3D && significance){plotDir = "Plots/2015/3D/significance";}
  
  if (isHI) {
    if(Prompt){
      fname = Form("/Jpsi_PbPb_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.pdf",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
      fnamepng = Form("/Jpsi_PbPb_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.png",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
      outfname = Form("/Jpsi_PbPb_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.root",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
    }
    if(nonPrompt){
      fname = Form("/Jpsi_PbPb_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.pdf",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
      fnamepng = Form("/Jpsi_PbPb_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.png",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
      outfname = Form("/Jpsi_PbPb_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f_Cent_%d-%d.root",efficiency,ymin,ymax,ptmin,ptmax,int(centmin/2.0),int(centmax/2.0));
    }
      
  } 
  else {
    if(Prompt){
      fname = Form("/Jpsi_pp_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",efficiency,ymin,ymax,ptmin,ptmax);
      fnamepng = Form("/Jpsi_pp_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.png",efficiency,ymin,ymax,ptmin,ptmax);
      outfname = Form("/Jpsi_pp_XcheckEff_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.root",efficiency,ymin,ymax,ptmin,ptmax);
    }
    if(nonPrompt){
      fname = Form("/Jpsi_pp_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.pdf",efficiency,ymin,ymax,ptmin,ptmax);
      fnamepng = Form("/Jpsi_pp_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.png",efficiency,ymin,ymax,ptmin,ptmax);
      outfname = Form("/Jpsi_pp_XcheckEff_NP_%3.1f_Rap_%3.1f-%3.1f_Pt_%3.1f-%3.1f.root",efficiency,ymin,ymax,ptmin,ptmax);
    }
  }
  std::cout << fname << std::endl;

  std::cout << "default: " << defaultCut.GetTitle() << std::endl;
  std::cout << "pt cut: " << ptCut.GetTitle() << std::endl;
  std::cout << "rapidity cut: " << rapCut.GetTitle() << std::endl;
  std::cout << "centrality cut: " << centCut.GetTitle() << std::endl;
  std::cout << "trigger bit: " << trigCut.GetTitle() << std::endl;

  TLatex *lpt;
  if (ptmin==0.0)
    lpt = new TLatex(0.38,0.845,Form("p_{T} < %3.1f GeV/c",ptmax));
  else{
    lpt = new TLatex(0.38,0.845,Form("%3.1f < p_{T} < %3.1f GeV/c",ptmin,ptmax));
    lpt->SetTextSize(0.033);
  }
  TLatex *lrap;
  if (absRapidity){
    if (ymin==0.0){
      lrap = new TLatex(0.38,0.77,Form("|y| < %3.1f",ymax));
      lrap->SetTextSize(0.033);
    }
    else{
      lrap = new TLatex(0.38,0.77,Form("%3.1f < |y| < %3.1f",ymin,ymax));
      lrap->SetTextSize(0.033);
    }
  }
  else {
    if (ymin==0.0)
      lrap = new TLatex(0.38,0.77,Form("y < %3.1f",ymax));
    else
      lrap = new TLatex(0.38,0.77,Form("%3.1f < y < %3.1f",ymin,ymax));
  }

  lpt->SetNDC(kTRUE);
  lrap->SetNDC(kTRUE);
  
  TLatex *lcent = NULL;
  if (isHI) {
    lcent = new TLatex(0.38,0.69,Form("Cent. %d-%d%%",int(centmin/2.0),int(centmax/2.0)));
    lcent->SetTextSize(0.033);
    lcent->SetNDC(kTRUE);
  }
  
  TChain *iTree = NULL;
  iTree = new TChain("hionia/myTree");
  if (isHI) {
    if(Prompt)
      {
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/Jpsi/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
      }
    if(nonPrompt)
      {
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	//iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/PbPb/BJpsi/OniaTree_Pythia8_BJpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
      }
  }
  else{
    if(Prompt)
      iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/pp/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
    if(nonPrompt)
      iTree->Add("/home/abdulla/CMS/quarkoniaData/Charm2015/rootfiles/MC/pp/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
  }
  
  Float_t         muPlPt;
  Float_t         muMiPt;
  Float_t         muPlEta;
  Float_t         muMiEta;
  Float_t         QQctau3D;
  Float_t         QQctau3DErr;
  
  double invariantMass;
  double JpsiPt;
  double JpsiRap;
  double vProb;

  float muPlDxy;
  float muPlDz;
  int muPlNhits;
  int muPlNPxlLayers;
  int muPlNTrkLayers;
  float muPlInnerChi;
  float muPlGlobalChi;
  int muPlNMuonhits;
  bool muPlGoodMu;
  bool muPlTrkArb;
  bool muPlTMOneStaTight;
  bool muPlHighPurity;

  float muMiDxy;
  float muMiDz;
  int muMiNhits;
  int muMiNPxlLayers;
  int muMiNTrkLayers;
  float muMiInnerChi;
  float muMiGlobalChi;
  int muMiNMuonhits;
  bool muMiGoodMu;
  bool muMiTrkArb;
  bool muMiTMOneStaTight;
  bool muMiHighPurity;
  
  
  int nx = 14, ny = 15000;
  //double xBins[22] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,30.0}; //19
  double xBins[15] = {3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,15.0,18.0,22.0,30.0}; //19
  
  //double xBins[4] = {3.0,6.5,12.0,30.0}; //19
  //double xBins[6] = {6.5,9.0,12.0,15.0,20.0,30.0}; //19
  
  //double xBins[13] = {6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,15.0,18.0,21.0,24.0,30.0}; //13
  //double xBins[8] = {6.0,9.0,12.0,15.0,18.0,21.0,24.0,30.0}; //19
  //TH2D *hPR2D = new TH2D("hPR2D","hPR2D",nx,xBins,ny,-5.0,10.0);
  
  //TH1F *hPR = new TH1F("hPR","hPR;Inv_mass",100,-2.5,3.2);
  TH1F *hPR = new TH1F("hPR","hPR;p_{T}",nx,xBins);
  //TH1F *hPR = new TH1F("hPR","hPR;l_{J/#psi} (mm);counts",1500,-5,10);
  hPR->Sumw2();
  hPR->SetLineColor(1);
  //TH1F *hPRfunc = new TH1F("hPRfunc","hPRfunc;Inv_mass",100,-2.8,3.2);
  TH1F *hPRfunc = new TH1F("hPRfunc","hPRfunc;p_{T}",nx,xBins);
  //TH1F *hPRfunc = new TH1F("hPRfunc","hPRfunc;l_{J/#psi} (mm);counts",1500,-5,10);
  hPRfunc->Sumw2();
  hPRfunc->SetLineColor(2);
  //TH1F* hPRCut = (TH1F*) hPR->Clone("hPRCut");
  
  TFile *Fitfile = new TFile("Plots/2015/3D/ctauFit/Jpsi_pp_eff_0.9_Rap_1.6-2.4_Pt_3.0-30.0.root");
  TGraph *FitGr = (TGraph*)Fitfile->Get("Graph");
  TF1 *func1 = (TF1*)FitGr->GetFunction("FitFn");
  
  new TCanvas;
  TAxis *XaxctauFW = FitGr->GetXaxis();
  XaxctauFW->SetLimits(0.0,30.0);
  TAxis *YaxctauFW = FitGr->GetYaxis();
  YaxctauFW->SetLimits(0.0,0.1);
  FitGr->GetXaxis()->SetRangeUser(0.0,30.0);
  FitGr->GetYaxis()->SetRangeUser(0.0,0.1);
  //FitGr->Draw("AP");
  //FitGr->Draw();
  func1->Draw();

  initOniaTree(iTree);
  Long64_t nentries = iTree->GetEntries();
  nentries = 500000;          
  cout<<"Total Entries in data  Tree  "<<"  "<<nentries<< "====="<<endl;
  Long64_t nbytes = 0, nb = 0;
 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = iTree->GetEntry(jentry);   
    
    for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
      TLorentzVector *qq4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
      TLorentzVector *qq4mupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
      TLorentzVector *qq4mumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
      
      invariantMass = qq4mom->M();
      JpsiPt = qq4mom->Pt();
      JpsiRap = qq4mom->Rapidity();

      muPlPt = qq4mupl->Pt();
      muMiPt = qq4mumi->Pt();
      muPlEta = qq4mupl->Eta();
      muMiEta = qq4mumi->Eta();

      QQctau3D = Reco_QQ_ctau3D[iQQ];
      //QQctauTrue3D = Reco_QQ_ctauTrue3D[iQQ];
      
      vProb = Reco_QQ_VtxProb[iQQ];
      muPlDxy=Reco_QQ_mupl_dxy[iQQ];
      muPlDz=Reco_QQ_mupl_dz[iQQ];
      muPlNPxlLayers=Reco_QQ_mupl_nPixWMea[iQQ];
      muPlNTrkLayers=Reco_QQ_mupl_nTrkWMea[iQQ];
      muPlTMOneStaTight=Reco_QQ_mupl_TMOneStaTight[iQQ];
      muPlGoodMu = Reco_QQ_mupl_isGoodMuon[iQQ];
      muPlHighPurity = Reco_QQ_mupl_highPurity[iQQ];

      muMiDxy=Reco_QQ_mumi_dxy[iQQ];
      muMiDz=Reco_QQ_mumi_dz[iQQ];
      muMiNPxlLayers=Reco_QQ_mumi_nPixWMea[iQQ];
      muMiNTrkLayers=Reco_QQ_mumi_nTrkWMea[iQQ];
      muMiTMOneStaTight=Reco_QQ_mumi_TMOneStaTight[iQQ];
      muMiGoodMu = Reco_QQ_mumi_isGoodMuon[iQQ];
      muMiHighPurity = Reco_QQ_mumi_highPurity[iQQ];
      
      bool defaultpass =0, ptpass =0, rappass =0, PosPass=0, NegPass=0, AllCut=0 ,PosIn=0, NegIn=0,TriggerPass=0;
      
      //if(Reco_QQ_sign[iQQ]==0 && invariantMass>2.85 && invariantMass<3.25 && Reco_QQ_ctau3D[iQQ]>-500 && Reco_QQ_ctauTrue3D[iQQ]>-10) {defaultpass =1;}
      if(Reco_QQ_sign[iQQ]==0 && invariantMass>2.85 && invariantMass<3.25 ) {defaultpass =1;}
      if(!defaultpass) continue;
      
      if(JpsiPt>ptmin && JpsiPt<ptmax)  {ptpass =1;}
      if( fabs(JpsiRap)>ymin && fabs(JpsiRap) <ymax )   {rappass =1;}

      if(!ptpass) continue;
      if(!rappass) continue; 
      //centCut = Form("Centrality>=%d&&Centrality<%d",centmin,centmax);
      
      if( (HLTriggers&(ULong64_t)pow(2,0))==pow(2,0) && (Reco_QQ_trig[iQQ]&(ULong64_t)pow(2,0))==pow(2,0) )
	{TriggerPass =1;}
      
      if (fabs(muPlEta) < 2.4 && ( (fabs(muPlEta) < 1.2 && muPlPt >= 3.5 ) || (1.2 <= fabs(muPlEta) && fabs(muPlEta) < 2.1 && muPlPt >= 5.77-1.89*fabs(muPlEta )) || (2.1 <= fabs(muPlEta) && muPlPt >= 1.8 )) )  {PosIn =1;}
	
      if(fabs(muMiEta) < 2.4 && ( (fabs(muMiEta) < 1.2 && muMiPt >= 3.5 ) || (1.2 <= fabs(muMiEta) && fabs(muMiEta) < 2.1 && muMiPt >= 5.77-1.89*fabs(muMiEta )) || (2.1 <= fabs(muMiEta) && muMiPt >= 1.8 )) )  {NegIn =1;}

      if( (muPlGoodMu==1) && muPlNTrkLayers> 5 &&  muPlNPxlLayers > 0 && TMath::Abs(muPlDxy) < 0.3 && TMath::Abs(muPlDz) < 20 && vProb > 0.01){PosPass=1;}

      if( (muMiGoodMu==1) && muMiNTrkLayers> 5 &&  muMiNPxlLayers > 0 && TMath::Abs(muMiDxy) < 0.3 && TMath::Abs(muMiDz) < 20 ){NegPass=1;}

      //if((Reco_QQ_mupl_isGoodMuon==1) && Reco_QQ_mupl_nTrkWMea> 5 &&  Reco_QQ_mupl_nPixWMea > 0 && TMath::fabs(Reco_QQ_mupl_dxy) < 0.3 && TMath::fabs(Reco_QQ_mupl_dz) < 20 && Reco_QQ_VtxProb > 0.01 )  {PosPass =1;}
      //if((Reco_QQ_mumi_isGoodMuon==1) && Reco_QQ_mumi_nTrkWMea> 5 &&  Reco_QQ_mumi_nPixWMea > 0 && TMath::fabs(Reco_QQ_mumi_dxy) < 0.3 && TMath::fabs(Reco_QQ_mumi_dz) < 20 && Reco_QQ_VtxProb > 0.01 )  {NegPass =1;}
      
      if(TriggerPass &&PosIn && NegIn && PosPass && NegPass) {AllCut = 1;}
      if(!AllCut) continue;
      
      //hPR->Fill(QQctau3D);
      hPR->Fill(JpsiPt);
      
      //for pp fwdrap:         cut.dMuon.ctauCut = "( ctau < (0.0128 + (0.288/pt)) )";
      //for pbpb fwdrap:     cut.dMuon.ctauCut = "( ctau < (0.0147 + (0.279/pt)) )";
      //for pp midrap:         cut.dMuon.ctauCut = "( ctau < (0.0097 + (0.249/pt)) )";
      //for pbpb midrap:     cut.dMuon.ctauCut = "( ctau < (0.0134 + (0.216/pt)) )";
      if( QQctau3D < (0.0128 + (0.288/JpsiPt) ) ) {hPRfunc->Fill(JpsiPt);}
      
    }
  }
  
  //TH1F* ratio = (TH1F*)hPR->Clone();      
  //ratio->Divide(hPRfunc);

  TH1F* ratio = (TH1F*)hPRfunc->Clone();      
  ratio->Divide(hPR);
  
  ratio->SetStats(0);
  //gStyle->SetFitFormat("1.3f");
  
  TCanvas *rt = new TCanvas();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  ratio->SetMarkerColor(kGreen+3);
  ratio->SetLineColor(0);
  ratio->SetMarkerStyle(29);
  ratio->SetMarkerSize(1.3);
  
  XaxctauFW = ratio->GetXaxis();
  XaxctauFW->SetLimits(0.0,30.0);
  YaxctauFW = ratio->GetYaxis();
  YaxctauFW->SetLimits(0.5,1.5);
  ratio->GetXaxis()->SetRangeUser(0.0,30.0);
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->Draw();
 
  //*/
  
  
  lpt->Draw();
  lrap->Draw();
  leg->Draw();
  
  if (isHI)
    lcent->Draw();

  if (savePlot) {
    rt->SaveAs(plotDir+fname);
    rt->SaveAs(plotDir+fnamepng);
  }

  return;

  /*
  //TH1D *ctauMid_pT = new TH1D("ctauMid_pT","ctauMid_pT",90,0,30);
  //TH1D *ctauFW_pT = new TH1D("ctauFW_pT","ctauFW_pT",90,0,30);
  
  //double ctauFrw[5] = {0.04,0.03,0.02,0.02,0.01};
  double ctauFrw[4] = {0.07,0.04,0.03};
  double ctauMd[6] = {0.04,0.03,0.02,0.02,0.01};
  
  double Fbins[4] = {3,6.5,12.0,30.0};
  double Mbins[6] = {6.5,9.0,12.0,15.0,20.0,30.0};
  double x1bins[4] = {0.0};
  double x2bins[6] = {0.0};
  
  for(int i=0; i<3; i++){
    x1bins[i] = (Fbins[i] + Fbins[i+1])/2.0;
    cout<<x1bins[i]<<endl;
    //ctauFW_pT->SetBinContent(x1bins[i], ctauFrw[i-1]);
  } 

  for(int i=0; i<5; i++){
    x2bins[i] = (Mbins[i] + Mbins[i+1])/2.0;
    //ctauMid_pT->SetBinContent(x2bins[i], ctauMd[i]);
  } 
  
  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();

  TGraph *ctauFW_pT = new TGraph(3,x1bins,ctauFrw);
  
  TF1 *FitFn = new TF1("FitFn",Pol2,0.0,30.0,3);
  FitFn->SetParNames("C","a","b");
  FitFn->SetLineWidth(2.0);
  FitFn->SetLineColor(2);
  ctauFW_pT->Fit("FitFn","", "",3.0,30.0);

  ctauFW_pT->SetLineColor(2);    
  ctauFW_pT->SetLineWidth(3); 
  ctauFW_pT->SetMarkerSize(1.5);
  ctauFW_pT->GetYaxis()->SetTitle("c#tau");
  ctauFW_pT->GetXaxis()->SetTitle("p_{T} bin");

  TAxis *XaxctauFW = ctauFW_pT->GetXaxis();
  XaxctauFW->SetLimits(0.0,30.0);
  TAxis *YaxctauFW = ctauFW_pT->GetYaxis();
  YaxctauFW->SetLimits(0.0,0.1);
  ctauFW_pT->GetXaxis()->SetRangeUser(0.0,30.0);
  ctauFW_pT->GetYaxis()->SetRangeUser(0.0,0.1);
  ctauFW_pT->Draw("AP*");

  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  
  TF1 *FitFn2 = new TF1("FitFn2",Exp,0.0,30.0,2);
  TGraph *ctauMid_pT = new TGraph(5,x2bins,ctauMd);
  ctauMid_pT->SetLineColor(2);    
  ctauMid_pT->SetLineWidth(3); 
  ctauMid_pT->SetMarkerSize(1.5);
  ctauMid_pT->GetYaxis()->SetTitle("c#tau");
  ctauMid_pT->GetXaxis()->SetTitle("p_{T} bin");
  
  ctauMid_pT->Fit("FitFn2","", "",3.0,30.0);
  FitFn2->SetParNames("C1","n");

  TAxis *XaxctauMid = ctauMid_pT->GetXaxis();
  XaxctauMid->SetLimits(0.0,30.0);
  TAxis *YaxctauMid = ctauMid_pT->GetYaxis();
  YaxctauMid->SetLimits(0.0,0.1);
  ctauMid_pT->GetXaxis()->SetRangeUser(0.0,30.0);
  ctauMid_pT->GetYaxis()->SetRangeUser(0.0,0.1);
  ctauMid_pT->Draw("AP*");
  
  TLatex *t = new TLatex();
  t->SetNDC(); t->SetTextAlign(12);t->SetTextColor(1);
  t->SetTextSize(0.03);
  
  //t->DrawLatex(0.135,0.21," 10  pp               6.5-30   0.0-1.6");
  //lh->Draw("same");
  
  TF1 *myfunc1=ctauFW_pT->GetFunction("FitFn");
  TF1 *myfunc2=ctauMid_pT->GetFunction("FitFn2");
  
  TH1F *funHist = new TH1F("funHist","funHist",30,0.0,30.0);
  TH1F *funHist2 = new TH1F("funHist2","funHist2",30,0.0,30.0);
  
  funHist->Eval(myfunc1);
  funHist2->Eval(myfunc2);
  //funHist = myfunc->GetHistogram();
  //root > TH1F h("h","",100,0,10);
  //root > h.Eval(f1);

  cout<<"funHist:Entries "<<funHist->GetEntries()<<endl;
  for(int i=1;i<=funHist->GetNbinsX();++i)
    {
      cout<<"bin:pT "<<i<<"ctau: "<<funHist->GetBinContent(i)<<endl<<endl;
    }

  for(int i=1;i<=funHist2->GetNbinsX();++i)
    {
      cout<<"bin:pT "<<i<<"ctau: "<<funHist2->GetBinContent(i)<<endl<<endl;
    }

  new TCanvas;
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridy();
  TAxis *XaxHistMid = funHist->GetXaxis();
  XaxHistMid->SetLimits(0.0,30.0);
  TAxis *YaxHistMid = funHist->GetYaxis();
  YaxHistMid->SetLimits(0.0,0.1);
  funHist->GetXaxis()->SetRangeUser(0.0,30.0);
  funHist->GetYaxis()->SetRangeUser(0.0,0.1);
  funHist->Draw();
  funHist2->SetLineColor(2);
  funHist2->Draw("same");
  //gPad->SaveAs("Plots/ctauMid.png");   
  //gPad->SaveAs("plots/ctauMid.pdf");   
  */

}
