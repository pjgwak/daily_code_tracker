#ifndef sigExtrChecks_C
#define sigExtrChecks_C

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "TLine.h"
#include "TObjArray.h"
#include "RooMCStudy.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooDataSet.h"

#include <vector>

#include "Macros/Utilities/resultUtils.h"

using namespace std;
using namespace RooFit;

struct poi {
  Char_t name[64];
  float val;
  float err;
};

void sigExtrChecks(
                    const char* workDirName,
                    const char* DSTag="DATA",
                    double nSamples = 500,
                    int check2Run=-1) {
  // workDirName: usual tag where to look for files in Output
  // module2Run: variable to choose which check to run (-1=all, 0=integrals check only, 1=toy mc check only )
  
  TString saveName("");
  bool runIntCheck = false;
  bool runToyMCCheck = false;
  if ( check2Run == -1 )
  {
    saveName = "sigExtrAllChecksResults";
    runIntCheck = true;
    runToyMCCheck = true;
  }
  else if ( check2Run == 0 )
  {
    saveName = "sigExtrPDFsChecksResults";
    runIntCheck = true;
  }
  else if ( check2Run == 1 )
  {
    saveName = "sigExtrToyMCChecksResults";
    runToyMCCheck = true;
  }
  else
  {
    cout << "Error: invalid option for check2Run parameter. Valid options are -1=all, 0=integrals check only, 1=toy mc check only" << endl;
    return;
  }
    
  //   bin edges
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  // model names
  Char_t jpsiName[128], psipName[128], bkgName[128];
  // collision system
  Char_t collSystem[8];
  
  // list of files
  vector<TString> theFiles = fileList(workDirName,"",DSTag);
  
  //______ Create arrays and histos to store results
  int nPart = 2;
  const char* aParticle[2] = {"Jpsi","Psi2S"};
  
  TObjArray* aDiffs = new TObjArray();
  aDiffs->SetOwner(kTRUE);
  for ( int i = 0 ; i < nPart ; i++)
  {
    TH1* hDiff = new TH1D(Form("hDiff_%s",aParticle[i]),"",theFiles.size(),0,theFiles.size());
    hDiff->GetYaxis()->SetTitle("(#int{DATA} - #int{PDF_{BKG}} - #int{PDF_{SIG}}) / #int{PDF_{SIG}}");
    hDiff->GetYaxis()->SetTitleOffset(1.15);
    hDiff->SetStats(0);
    hDiff->SetDirectory(0);
    hDiff->SetMarkerColor(1);
    hDiff->SetMarkerStyle(20);
    hDiff->SetMarkerSize(1);
    hDiff->SetLineColor(1);
    
    aDiffs->Add(hDiff);
  }
  
  TH1* hNJpsiToy = new TH1D("hNJpsiToy","",theFiles.size(),0,theFiles.size());
  hNJpsiToy->GetYaxis()->SetTitle("N_{J#psi}");
  hNJpsiToy->GetYaxis()->SetTitleOffset(1.15);
  hNJpsiToy->SetStats(0);
  hNJpsiToy->SetDirectory(0);
  hNJpsiToy->SetMarkerColor(2);
  hNJpsiToy->SetMarkerStyle(20);
  hNJpsiToy->SetMarkerSize(1);
  hNJpsiToy->SetLineColor(2);
  
  TH1* hNJpsiNom = new TH1D("hNJpsiNom","",theFiles.size(),0,theFiles.size());
  hNJpsiNom->GetYaxis()->SetTitle("N_{J#psi}");
  hNJpsiNom->GetYaxis()->SetTitleOffset(1.15);
  hNJpsiNom->SetStats(0);
  hNJpsiNom->SetDirectory(0);
  hNJpsiNom->SetMarkerColor(1);
  hNJpsiNom->SetMarkerStyle(20);
  hNJpsiNom->SetMarkerSize(1);
  hNJpsiNom->SetLineColor(1);
  
  TH1* hNJpsiDiff = new TH1D("hNJpsiDiff","",theFiles.size(),0,theFiles.size());
  hNJpsiDiff->GetYaxis()->SetTitle("N_{J#psi}^{nom} - N_{J#psi}^{toy}");
  hNJpsiDiff->GetYaxis()->SetTitleOffset(1.15);
  hNJpsiDiff->SetStats(0);
  hNJpsiDiff->SetDirectory(0);
  hNJpsiDiff->SetMarkerColor(1);
  hNJpsiDiff->SetMarkerStyle(20);
  hNJpsiDiff->SetMarkerSize(1);
  hNJpsiDiff->SetLineColor(1);
  
  TH1* hRPsiPJpsiToy = new TH1D("hRPsiPJpsiToy","",theFiles.size(),0,theFiles.size());
  hRPsiPJpsiToy->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J#psi}");
  hRPsiPJpsiToy->GetYaxis()->SetTitleOffset(1.15);
  hRPsiPJpsiToy->SetStats(0);
  hRPsiPJpsiToy->SetDirectory(0);
  hRPsiPJpsiToy->SetMarkerColor(2);
  hRPsiPJpsiToy->SetMarkerStyle(20);
  hRPsiPJpsiToy->SetMarkerSize(1);
  hRPsiPJpsiToy->SetLineColor(2);
  
  TH1* hRPsiPJpsiNom = new TH1D("hRPsiPJpsiNom","",theFiles.size(),0,theFiles.size());
  hRPsiPJpsiNom->GetYaxis()->SetTitle("N_{#psi(2S)}/N_{J#psi}");
  hRPsiPJpsiNom->GetYaxis()->SetTitleOffset(1.15);
  hRPsiPJpsiNom->SetStats(0);
  hRPsiPJpsiNom->SetDirectory(0);
  hRPsiPJpsiNom->SetMarkerColor(1);
  hRPsiPJpsiNom->SetMarkerStyle(20);
  hRPsiPJpsiNom->SetMarkerSize(1);
  hRPsiPJpsiNom->SetLineColor(1);
  
  TH1* hRPsiPJpsiDiff = new TH1D("hRPsiPJpsiDiff","",theFiles.size(),0,theFiles.size());
  hRPsiPJpsiDiff->GetYaxis()->SetTitle("(N_{#psi(2S)}/N_{J#psi})^{nom} - (N_{#psi(2S)}/N_{J#psi})^{toy}");
  hRPsiPJpsiDiff->GetYaxis()->SetTitleOffset(1.15);
  hRPsiPJpsiDiff->SetStats(0);
  hRPsiPJpsiDiff->SetDirectory(0);
  hRPsiPJpsiDiff->SetMarkerColor(1);
  hRPsiPJpsiDiff->SetMarkerStyle(20);
  hRPsiPJpsiDiff->SetMarkerSize(1);
  hRPsiPJpsiDiff->SetLineColor(1);

  
  TObjArray* aToys = new TObjArray();
  aToys->SetOwner(kTRUE);
  //_____
  
  
  int cnt=0;
  for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) //Start loop on files
  {
    cout << "Parsing file " << cnt << " / " << theFiles.size() << ": " << *it << endl;
    
    // parse the file name to get info
    anabin thebin = binFromFile(*it);
    ptmin = thebin.ptbin().low();
    ptmax = thebin.ptbin().high();
    ymin = thebin.rapbin().low();
    ymax = thebin.rapbin().high();
    centmin = thebin.centbin().low();
    centmax = thebin.centbin().high();
    strcpy(collSystem, (it->Index("PbPb")>0) ? "PbPb" : "PP");
    
    TString binName(Form("%s-Pt[%.1f,%.1f]-Y[%.1f,%.1f]-C[%.1f,%.1f]",collSystem,ptmin,ptmax,ymin,ymax,centmin,centmax));
    
    TFile *f = TFile::Open(*it); RooWorkspace *ws = NULL;
    if (!f)
    {
      cout << "Error, file " << *it << " does not exist." << endl;
      continue;
    }
    else
    {
      ws = (RooWorkspace*) f->Get("workspace");
      if (!ws)
      {
        cout << "Error, workspace not found in " << *it << "." << endl;
        f->Close();
        delete f;
        continue;
      }
    }
    
    // get the models and variables from the ws
    TString pdfName(Form("pdfMASS_Tot_%s",collSystem));
    RooFitResult *fitresult = static_cast<RooFitResult*>(ws->genobj(Form("fitResult_%s",pdfName.Data())));
    
    RooAddPdf *model = static_cast<RooAddPdf*>(pdfFromWS(ws, Form("_%s",collSystem), "pdfMASS_Tot")); //Get the total pdf
    if(!model)
    {
      cout << "Error: Missing total model" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooAbsPdf *model_sig = pdfFromWS(ws, Form("_%s",collSystem), "pdfMASSTot_Jpsi");
    if(!model_sig)
    {
      cout << "Error: Missing signal model" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooAbsPdf *model_sigPsip = pdfFromWS(ws, Form("_%s",collSystem), "pdfMASSTot_Psi2S");
    if(!model_sigPsip)
    {
      cout << "Error: Missing signal model" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooAbsPdf *model_bkg = pdfFromWS(ws, Form("_%s",collSystem), "pdfMASSTot_Bkg");
    if(!model_bkg)
    {
      cout << "Error: Missing background model" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooRealVar *invMass = ws->var("invMass");
    if(!invMass)
    {
      cout << "Error: Missing invariant mass variable" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooRealVar *normPsi = ws->var(Form("N_Jpsi_%s",collSystem));
    if(!normPsi)
    {
      cout << "Error: Missing nof psi variable" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooRealVar *normPsiP = static_cast<RooRealVar*>(ws->function(Form("N_Psi2S_%s",collSystem)));
    if(!normPsiP)
    {
      cout << "Error: Missing nof psi variable" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooRealVar *normBkg = ws->var(Form("N_Bkg_%s",collSystem));
    if(!normBkg)
    {
      cout << "Error: Missing nof bkg variable" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    RooRealVar *Psi2SOverJpsi = ws->var(Form("RFrac2Svs1S_%s",collSystem));
    if(!Psi2SOverJpsi)
    {
      cout << "Error: Missing yields ratio variable" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    const char* token = Form("_%s",collSystem);
    RooAbsData *dat = dataFromWS(ws, token, Form("dOS_%s", DSTag));
    if(!dat)
    {
      cout << "Error: Missing data" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    // create a clone of the workspace to work in
    RooWorkspace* wsC = static_cast<RooWorkspace*>(ws->Clone("workspaceClone"));
    if(!wsC)
    {
      cout << "Error: Not able to clone workspace" << endl;
      
      delete ws;
      delete f;
      continue;
    }
    
    //_____________________________________________________________________
    //____________Computation of Int{DATA} - Int{PDF_Bkg} =? Int{PDF_Sig}__
    //_____________________________________________________________________
    if (runIntCheck)
    {
      for ( int i = 0 ; i < nPart ; i++ )
      {
        double mMin(0.), mMax(0.);
        const char* range = "";
        bool isJpsi;
        if ( strcmp(aParticle[i],"Jpsi") )
        {
          mMin = 3.4;
          mMax = 4.5;
          range = "invMass>3.4";
          isJpsi = false;
        }
        else
        {
          mMin = 2.2;
          mMax = 3.4;
          range = "invMass<3.4";
          isJpsi = true;
        }
        invMass->setRange("r1",mMin,mMax);
        
        double intData = dat->reduce(range)->sumEntries();
        
        RooAbsReal* intSig(0x0);
        double nSig(0);
        if (isJpsi)
        {
          intSig = model_sig->createIntegral(*invMass,NormSet(*invMass),Range("r1"));
          nSig = normPsi->getVal()*intSig->getVal();
        }
        else
        {
          intSig = model_sigPsip->createIntegral(*invMass,NormSet(*invMass),Range("r1"));
          nSig = normPsiP->getVal()*intSig->getVal();
        }
        double nSigErr = intSig->getPropagatedError(*fitresult);
        if ( nSigErr == 0 ) nSigErr = (isJpsi ? normPsi->getError() : normPsiP->getError());
        else nSigErr *= (isJpsi ? normPsi->getVal() : normPsiP->getVal());
        
        RooAbsReal* intBkg = model_bkg->createIntegral(*invMass,NormSet(*invMass),Range("r1"));
        double nBkg = normBkg->getVal()*intBkg->getVal();
        double nBkgErr = intBkg->getPropagatedError(*fitresult);
        if ( nBkgErr == 0 ) nBkgErr = normBkg->getError();
        else nBkgErr *= normBkg->getVal();
        
        double relDiff = (intData - nBkg - nSig)/nSig;
        double relDiffErr = TMath::Sqrt( TMath::Power(TMath::Sqrt(intData)/nSig,2) + TMath::Power(nBkgErr/nSig,2) /*+ TMath::Power(2.*nSigErr/nSig,2)*/ );
        
        TH1* hDiff = static_cast<TH1*>(aDiffs->FindObject(Form("hDiff_%s",aParticle[i])));
        hDiff->SetBinContent(cnt+1,relDiff);
        hDiff->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
        hDiff->SetBinError(cnt+1,relDiffErr);
      }
    }
    //___________________________________________________________________
    //___________________________________________________________________
    
    
    //______________________________________________
    //_____________________ Toy MC study ___________
    //______________________________________________
    if (runToyMCCheck)
    {
      double nJpsi = normPsi->getVal(); // Nof Jpsi from the nominal fit in DATA
      double nJpsiErr = normPsi->getError();
      
      double nPsi2S = normPsiP->getVal(); // Nof Psi2S from the nominal fit in DATA
      double nPsi2SErr = normPsiP->getError();
      
      double nPsi2SOverJpsi = Psi2SOverJpsi->getVal(); // Ratio Psi2S over Jpsi from the nominal fit in DATA
      double nPsi2SOverJpsiErr = Psi2SOverJpsi->getError();
      
      
      double nData = dat->sumEntries();
      
      hNJpsiNom->SetBinContent(cnt+1,nJpsi);
      hNJpsiNom->SetBinError(cnt+1,nJpsiErr);
      hNJpsiNom->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      hRPsiPJpsiNom->SetBinContent(cnt+1,nPsi2SOverJpsi);
      hRPsiPJpsiNom->SetBinError(cnt+1,nPsi2SOverJpsiErr);
      hRPsiPJpsiNom->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      wsC->cd();
      if(!wsC->loadSnapshot(Form("%s_parIni",pdfName.Data()))) //Get the snapshot with the initial parameters of the original fit
      {
        cout << "Error loading initial parameters workspace" << endl;
        delete wsC;
        delete ws;
        delete f;
        continue;
      }
      
      RooAbsPdf* model_init = pdfFromWS(wsC, Form("_%s",collSystem), "pdfMASS_Tot"); //Get the total pdf with the initial parameters for fitting
      RooRealVar *invMass_init = wsC->var("invMass");
      RooRealVar *normPsi_init = wsC->var(Form("N_Jpsi_%s",collSystem));
      RooRealVar *normPsi2S_init = static_cast<RooRealVar*>(wsC->function(Form("N_Psi2S_%s",collSystem)));
      RooRealVar *Psi2SOverJpsi_init = wsC->var(Form("RFrac2Svs1S_%s",collSystem));
      
      RooAddPdf *sumPDF = new RooAddPdf("sumPDF","",RooArgList(*model_sig,*model_sigPsip,*model_bkg),RooArgList(*normPsi,*normPsiP,*normBkg)); // Create total pdf from the sum of the pdfs (otherwise it does not work, Idk why)
      
      // Instantiate MC study manager
      RooMCStudy* mcstudy = new RooMCStudy(*sumPDF,*invMass_init,FitModel(*model_init),Extended(),Silence(),FitOptions(NumCPU(8), Extended(kTRUE)));//
      
      // Generate and fit nSamples samples of Poisson(nData) events
      mcstudy->generateAndFit(nSamples,nData,kFALSE);
      
      //_____ Plot and fit NJpsi, NPsi2S and NPsi2S/NJpsi distributions
      TCanvas* cdistrJpsi = new TCanvas(Form("cNJpsiDistr_%s",binName.Data()));
      const RooDataSet* jpsiDataSet = &(mcstudy->fitParDataSet());
      
      normPsi_init->setRange("FitNJpsiWindow",nJpsi-4*nJpsiErr,nJpsi+4*nJpsiErr);
      
      RooRealVar meanJPsi("meanJPsi","average NJPsi",nJpsi,nJpsi-(0.05*nJpsi),nJpsi+(0.05*nJpsi));
      RooRealVar sigmaJPsi("sigmaJPsi","average sigmaNJpsi",nJpsiErr,0.9*nJpsiErr,1.1*nJpsiErr);
      RooGaussian gausJpsi("gausJpsi","gaussian fit",*normPsi_init,meanJPsi,sigmaJPsi);
      RooFitResult* fitResultJpsi = gausJpsi.fitTo(*(RooAbsData*)jpsiDataSet,Range("FitNJpsiWindow"));
      
      RooRealVar meanJPsiNom("meanJPsiNom","average NJPsi",nJpsi);
      RooRealVar sigmaJPsiNom("sigmaJPsiNom","average sigmaNJpsi",nJpsiErr);
      RooGaussian gausJpsiNom("gausJpsiNom","gaussian fit",*normPsi_init,meanJPsiNom,sigmaJPsiNom);
      
      RooPlot* frameJpsi = normPsi_init->frame(Bins(100),Range(nJpsi-5*nJpsiErr,nJpsi+6*nJpsiErr));
      frameJpsi->SetTitle(binName.Data());
      jpsiDataSet->plotOn(frameJpsi);
      gausJpsi.plotOn(frameJpsi,LineColor(kRed), LineStyle(1), Precision(1e-4),Range(frameJpsi->GetXaxis()->GetXmin(),frameJpsi->GetXaxis()->GetXmax()));
      gausJpsiNom.plotOn(frameJpsi,LineColor(kBlue), LineStyle(1), Precision(1e-4),Range(frameJpsi->GetXaxis()->GetXmin(),frameJpsi->GetXaxis()->GetXmax()));
      frameJpsi->Draw();
      
      TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.036); float dy = 0.04;
      t->DrawLatex(0.45, 0.85-dy, Form("#color[4]{N_{J/#psi}^{nom} = %0.2f #pm %0.2f}",nJpsi,nJpsiErr)); dy+=0.050;
      t->DrawLatex(0.45, 0.85-dy, Form("#color[2]{N_{J/#psi}^{toy} = %0.2f(%0.2f) #pm %0.2f(%0.2f)}",meanJPsi.getVal(),meanJPsi.getError(),sigmaJPsi.getVal(),sigmaJPsi.getError())); dy+=0.050;
      t->Draw("same");
      
      cdistrJpsi->Update();
      cdistrJpsi->SaveAs(Form("cNJpsiDistr_%s.pdf",binName.Data()));
      aToys->Add(cdistrJpsi);
      
      hNJpsiToy->SetBinContent(cnt+1,meanJPsi.getVal());
      hNJpsiToy->SetBinError(cnt+1,sigmaJPsi.getVal());
      hNJpsiToy->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      hNJpsiDiff->SetBinContent(cnt+1,nJpsi - meanJPsi.getVal());
      hNJpsiDiff->SetBinError(cnt+1,meanJPsi.getError());
      hNJpsiDiff->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      
      TCanvas* cdistrPsi2SOverJpsi = new TCanvas(Form("cPsi2SOverJpsiDistr_%s",binName.Data()));
      const RooDataSet* psi2sOverJpsiDataSet = &(mcstudy->fitParDataSet());
      
      Psi2SOverJpsi_init->setRange("FitRatioWindow",nPsi2SOverJpsi-4*nPsi2SOverJpsiErr,nPsi2SOverJpsi+4*nPsi2SOverJpsiErr);
      
      RooRealVar meanNPsi2SOverJPsi("meanNPsi2SOverJPsi","average NJPsi",nPsi2SOverJpsi,nPsi2SOverJpsi-(0.05*nPsi2SOverJpsi),nPsi2SOverJpsi+(0.05*nPsi2SOverJpsi));
      RooRealVar sigmaNPsi2SOverJPsi("sigmaNPsi2SOverJPsi","average sigmaNJpsi",nPsi2SOverJpsiErr,0.9*nPsi2SOverJpsiErr,1.1*nPsi2SOverJpsiErr);
      RooGaussian gausNPsi2SOverJPsi("gausNPsi2SOverJPsi","gaussian fit",*Psi2SOverJpsi_init,meanNPsi2SOverJPsi,sigmaNPsi2SOverJPsi);
      RooFitResult* fitResultNPsi2SOverJPsi = gausNPsi2SOverJPsi.fitTo(*(RooAbsData*)psi2sOverJpsiDataSet,Range("FitRatioWindow"));
      
      RooRealVar meanNPsi2SOverJPsiNom("meanNPsi2SOverJPsiNom","average NJPsi",nPsi2SOverJpsi);
      RooRealVar sigmaNPsi2SOverJPsiNom("sigmaNPsi2SOverJPsiNom","average sigmaNJpsi",nPsi2SOverJpsiErr);
      RooGaussian gausNPsi2SOverJPsiNom("gausNPsi2SOverJPsiNom","gaussian fit",*Psi2SOverJpsi_init,meanNPsi2SOverJPsiNom,sigmaNPsi2SOverJPsiNom);
      
      RooPlot* frameNPsi2SOverJPsi = Psi2SOverJpsi_init->frame(Bins(100),Range(nPsi2SOverJpsi-5*nPsi2SOverJpsiErr,nPsi2SOverJpsi+6*nPsi2SOverJpsiErr));
      frameNPsi2SOverJPsi->SetTitle(binName.Data());
      psi2sOverJpsiDataSet->plotOn(frameNPsi2SOverJPsi);
      gausNPsi2SOverJPsi.plotOn(frameNPsi2SOverJPsi,LineColor(kRed), LineStyle(1), Precision(1e-4),Range(frameNPsi2SOverJPsi->GetXaxis()->GetXmin(),frameNPsi2SOverJPsi->GetXaxis()->GetXmax()));
      gausNPsi2SOverJPsiNom.plotOn(frameNPsi2SOverJPsi,LineColor(kBlue), LineStyle(1), Precision(1e-4),Range(frameNPsi2SOverJPsi->GetXaxis()->GetXmin(),frameNPsi2SOverJPsi->GetXaxis()->GetXmax()));
      frameNPsi2SOverJPsi->Draw();
      
      TLatex *t2 = new TLatex(); t2->SetNDC(); t2->SetTextSize(0.036); dy = 0.04;
      t2->DrawLatex(0.35, 0.85-dy, Form("#color[4]{(N_{#psi(2S)}/N_{J#psi})^{nom} = %0.4f #pm %0.4f}",nPsi2SOverJpsi,nPsi2SOverJpsiErr)); dy+=0.050;
      t2->DrawLatex(0.35, 0.85-dy, Form("#color[2]{(N_{#psi(2S)}/N_{J#psi})^{toy} = %0.4f(%0.4f) #pm %0.4f(%0.4f)}",meanNPsi2SOverJPsi.getVal(),meanNPsi2SOverJPsi.getError(),sigmaNPsi2SOverJPsi.getVal(),sigmaNPsi2SOverJPsi.getError())); dy+=0.050;
      t2->Draw("same");
      
      cdistrPsi2SOverJpsi->Update();
      cdistrPsi2SOverJpsi->SaveAs(Form("cPsi2OverJpsiSDistr_%s.pdf",binName.Data()));
      aToys->Add(cdistrPsi2SOverJpsi);
      
      hRPsiPJpsiToy->SetBinContent(cnt+1,meanNPsi2SOverJPsi.getVal());
      hRPsiPJpsiToy->SetBinError(cnt+1,sigmaNPsi2SOverJPsi.getVal());
      hRPsiPJpsiToy->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      hRPsiPJpsiDiff->SetBinContent(cnt+1,nPsi2SOverJpsi - meanNPsi2SOverJPsi.getVal());
      hRPsiPJpsiDiff->SetBinError(cnt+1,meanNPsi2SOverJPsi.getError());
      hRPsiPJpsiDiff->GetXaxis()->SetBinLabel(cnt+1,binName.Data());
      
      //___ Save toy MC to file
      mcstudy->SetName(Form("MCStudy_%s",binName.Data()));
      
      TString fileName(it->Data());
      fileName.Remove(0,fileName.Last('/')+1);
      TFile* fSaveMCStudy = new TFile(Form("MCTOY_%s",fileName.Data()),"RECREATE");
      mcstudy->Write();
      fSaveMCStudy->Close(); delete fSaveMCStudy;
      delete mcstudy;
      //_____
    }
    //___________________________________________________________________
    //___________________________________________________________________
    
    
    delete ws;
    delete wsC;
    delete f;
    
    
    cnt++;
  } //End loop on the files
  
  
  //___ Draw and Save
  if (runIntCheck)
  {
    for ( int i = 0 ; i < nPart ; i++ )
    {
      TCanvas* c = new TCanvas(Form("relDiff_%s",aParticle[i]),"",90,116,1265,535);
      c->Range(-3.690909,-0.01066472,33.30606,0.01252061);
      c->SetFillColor(0);
      c->SetBorderMode(0);
      c->SetBorderSize(2);
      c->SetRightMargin(0.1163896);
      c->SetTopMargin(0.03732809);
      c->SetBottomMargin(0.1630648);
      c->SetFrameBorderMode(0);
      c->SetFrameBorderMode(0);
      
      TH1* hDiff = static_cast<TH1*>(aDiffs->FindObject(Form("hDiff_%s",aParticle[i])));
      hDiff->Draw("p");
      
      TLine* l = new TLine(0.,0.,hDiff->GetXaxis()->GetXmax(),0.);
      l->SetLineWidth(3);
      l->Draw("same");
      
      c->SaveAs(Form("%sSigExtrXCheck.pdf",aParticle[i]));
      aDiffs->Add(c);
    }
  }
  
  TFile* fSave = new TFile(Form("%s.root",saveName.Data()),"RECREATE");
  if (runIntCheck) aDiffs->Write("relDiffs", TObject::kOverwrite | TObject::kSingleKey);
  if (runToyMCCheck)
  {
    TH1* hNJpsiNomClone = static_cast<TH1*>(hNJpsiNom->Clone("hNJpsiNomClone"));
    TH1* hNJpsiToyClone = static_cast<TH1*>(hNJpsiToy->Clone("hNJpsiToyClone"));
    aToys->Add(hNJpsiNom);
    aToys->Add(hNJpsiToy);
    
    TCanvas* c1 = new TCanvas("TestJpsiYield","",90,116,1265,535);
    c1->Range(-3.690909,-0.01066472,33.30606,0.01252061);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetRightMargin(0.1163896);
    c1->SetTopMargin(0.03732809);
    c1->SetBottomMargin(0.1630648);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
    
    TLegend* l1 = new TLegend(0.1284695,0.7362205,0.5463918,0.9271654);
    l1->SetBorderSize(0);
    l1->SetLineColor(1);
    l1->SetLineStyle(1);
    l1->SetLineWidth(1);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->AddEntry(hNJpsiNomClone,"Nominal","lp");
    l1->AddEntry(hNJpsiToyClone,"Toy MC","lp");
    
    hNJpsiNomClone->Draw("p");
    hNJpsiToyClone->Draw("samep");
    l1->Draw("same");
    
    c1->SaveAs("NJpsiToyMCCheck.pdf");
    aToys->Add(c1);
    
    
    TH1* hRPsiPJpsiNomClone = static_cast<TH1*>(hRPsiPJpsiNom->Clone("hRPsiPJpsiNomClone"));
    TH1* hRPsiPJpsiToyClone = static_cast<TH1*>(hRPsiPJpsiToy->Clone("hRPsiPJpsiToyClone"));
    aToys->Add(hRPsiPJpsiNom);
    aToys->Add(hRPsiPJpsiToy);
    
    TCanvas* c2 = new TCanvas("TestYieldsRatio","",90,116,1265,535);
    c2->Range(-3.690909,-0.01066472,33.30606,0.01252061);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    c2->SetRightMargin(0.1163896);
    c2->SetTopMargin(0.03732809);
    c2->SetBottomMargin(0.1630648);
    c2->SetFrameBorderMode(0);
    c2->SetFrameBorderMode(0);
    
    TLegend* l2 = new TLegend(0.1284695,0.7362205,0.5463918,0.9271654);
    l2->SetBorderSize(0);
    l2->SetLineColor(1);
    l2->SetLineStyle(1);
    l2->SetLineWidth(1);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->AddEntry(hRPsiPJpsiNomClone,"Nominal","lp");
    l2->AddEntry(hRPsiPJpsiToyClone,"Toy MC","lp");
    
    hRPsiPJpsiNomClone->Draw("p");
    hRPsiPJpsiToyClone->Draw("samep");
    l2->Draw("same");
    
    c2->SaveAs("PsiPOverJpsiToyMCCheck.pdf");
    aToys->Add(c2);
    
    
    TH1* hNJpsiDiffClone = static_cast<TH1*>(hNJpsiDiff->Clone("hNJpsiDiffClone"));
    aToys->Add(hNJpsiDiff);

    TCanvas* c3 = new TCanvas("TestJpsiYieldDiff","",90,116,1265,535);
    c3->Range(-3.690909,-0.01066472,33.30606,0.01252061);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetBorderSize(2);
    c3->SetRightMargin(0.1163896);
    c3->SetTopMargin(0.03732809);
    c3->SetBottomMargin(0.1630648);
    c3->SetFrameBorderMode(0);
    c3->SetFrameBorderMode(0);
    
    hNJpsiDiffClone->Draw("p");
    
    TLine* l0 = new TLine(0.,0.,hNJpsiDiffClone->GetXaxis()->GetXmax(),0.);
    l0->SetLineWidth(3);
    l0->Draw("same");
    
    c3->SaveAs("NJpsiToyMCCheckDiffs.pdf");
    aToys->Add(c3);
    
    
    TH1* hRPsiPJpsiDiffClone = static_cast<TH1*>(hRPsiPJpsiDiff->Clone("hRPsiPJpsiDiffClone"));
    aToys->Add(hRPsiPJpsiDiff);
    
    TCanvas* c4 = new TCanvas("TestYieldsRatioDiff","",90,116,1265,535);
    c4->Range(-3.690909,-0.01066472,33.30606,0.01252061);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetBorderSize(2);
    c4->SetRightMargin(0.1163896);
    c4->SetTopMargin(0.03732809);
    c4->SetBottomMargin(0.1630648);
    c4->SetFrameBorderMode(0);
    c4->SetFrameBorderMode(0);
    
    hRPsiPJpsiDiffClone->Draw("p");
    
    TLine* l01 = new TLine(0.,0.,hRPsiPJpsiDiffClone->GetXaxis()->GetXmax(),0.);
    l01->SetLineWidth(3);
    l01->Draw("same");
    
    c4->SaveAs("PsiPOverJpsiToyMCCheckDiffs.pdf");
    aToys->Add(c4);
    
    
    aToys->Write("aToys", TObject::kOverwrite | TObject::kSingleKey);
  }
  fSave->Close();
}

#endif // #ifndef sigExtrChecks_C
