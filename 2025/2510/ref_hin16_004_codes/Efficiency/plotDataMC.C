#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"

#include "FitFunctions.h"
#include "../Fitter/Macros/CMS/CMS_lumi.C"
#include "../Fitter/Macros/CMS/tdrstyle.C"
#include "../Fitter/Macros/Utilities/texUtils.h"

#include <iostream>
#include <fstream>

using namespace std;

void setErr(TH1F *hist);
void setErrWt(TH1F *hist1, TH1F *hist2);
void fixCentPp(TH1F *hist);
TH1F* integrateHist(TH1F *hist);

void plotDataMC() {

  //plot Data_MC comparison, data points from the fits

  TString colltag, deptag, raptag, cuttag, objtag;
  for (int icoll=0; icoll<1; icoll++) {
    colltag = (icoll==0) ? "pp" : "pbpb";
    
    for (int iobj=0; iobj<1; iobj++) {
      objtag = (iobj==0) ? "jpsi" : "psi2s";
      
      TString name = "files/histos_" + objtag + "_" + colltag + ".root";
      TFile *fjpsi = new TFile(name);
            
      for (int idep=0; idep<1; idep++) { // idep=2 -> integrated
	deptag = (idep==0) ? "pt" : "cent";
	
	for (int irap=0; irap<2; irap++) {
	  raptag = (irap==0) ? "mid" : "fwd";
	  
	  for (int icut=4; icut<5; icut++) { // icut=3 -> pt-integrated cut efficiency, icut=4 -> pt-dep cut efficiency
	    if(icut==0) cuttag = "_";
	    else if (icut==1 || icut==3) cuttag = "cut_";
	    else cuttag = "ptdepcut_";
	    
	    setTDRStyle();
	    gStyle->SetEndErrorSize(3);
	    
	    TString hname;
	    hname = "hnum" + cuttag + deptag + raptag;
	    TH1F *hjpsinum = (TH1F*) fjpsi->Get(hname);
	    hname = (icut<3) ? "hden_" + deptag + raptag : "hnum_" + deptag + raptag;
	    TH1F *hjpsiden = (TH1F*) fjpsi->Get(hname);
	    
	    if (idep==2) {
	      hjpsinum = integrateHist(hjpsinum);
	      hjpsiden = integrateHist(hjpsiden);
	    }
	    
	    if (icoll==0 && idep==1) { // centrality for pp... fill all the bins
	      fixCentPp(hjpsinum);
	      fixCentPp(hjpsiden);
	    }
	    
	    const int nbins_ptmid = 5;
	    const float bins_ptmid[nbins_ptmid+1] = {6.5, 9, 12, 15, 20, 30};
	    const int nbins_ptfwd = 3;
	    const float bins_ptfwd[nbins_ptfwd+1] = {3, 6.5, 12, 30};
	    TH1F *hjpsiMid;
	    TH1F *hjpsiFwd;
	    
	    int nbins = hjpsinum->GetNbinsX();
	    cout<<" nbibs: "<<nbins<<endl;
	    
	    if(irap==0)
	      hjpsiMid = new TH1F("hjpsiMid","hjpsiMid",nbins_ptmid,bins_ptmid);
	    if(irap==1)
	      hjpsiFwd = new TH1F("hjpsiFwd","hjpsiFwd",nbins_ptfwd,bins_ptfwd);
	    
	    nbins = (irap==0) ? 5 : 3;
	    std::vector<double> ptMid = {58653.734,49509.714,18024.597,8912.3154,3088.4191}; // pp data Jpsi
	    std::vector<double> ptMidErr = {250.80012,230.37362,137.7950,97.365753,57.877162};
	    
	    std::vector<double> ptFwd = {80647.242,64563.777,8649.4736,};
	    std::vector<double> ptFwdErr = {305.99594, 265.74569,97.767669};
	    
	    if(icoll == 0 && iobj==1){
	      ptMid = {2277.2612,2057.1032,776.04791,399.03088,147.03543}; // pp data Psi2s
	      ptMidErr = {63.425487,57.534194,34.686260,25.240642,15.657732};
	      
	      ptFwd = {2893.1272,2243.7263,372.88653};
	      ptFwdErr = {89.004768,66.136161,26.768508};
	    }
	    
	    if(icoll == 1 && iobj==0){
	      ptMid = {7890.9941,7001.1967,2821.2739,1414.6502,642.99639}; // PbPb data jpsi
	      ptMidErr = {105.10430,90.561050,55.901050,39.513504,26.796726};
	      
	      ptFwd = {9658.9218,7334.6621,1232.4993};
	      ptFwdErr = {172.92077,115.03877,38.281230};
	    }
	    
	    if(icoll == 1 && iobj==1){
	      ptMid = {104.72235,119.38916,61.021923,28.467804,10.201046}; // PbPb data jpsi
	      ptMidErr = {47.020084,26.660514,13.654704,9.0622406,5.9582667};
	      
	      ptFwd = {63.292587,119.99708,36.432670,};
	      ptFwdErr = {118.25901,53.776134,12.202820};
	    }
	    
	    for (int i=0; i<nbins; i++) {
	      if(irap==0){
		hjpsiMid->SetBinContent(i+1,ptMid[i]);
		hjpsiMid->SetBinError(i+1, ptMidErr[i]);

		cout<<" ptMid[i] "<<ptMid[i]<<endl;
	      }
	      
	      if(irap==1){
		hjpsiFwd->SetBinContent(i+1,ptFwd[i]);
		hjpsiFwd->SetBinError(i+1, ptFwdErr[i]);
	      }
	    }
	    
	    gStyle->SetOptTitle(0); // at least most of the time                                                            
	    gStyle->SetOptStat(0); // most of the time, sometimes "nemriou" might be useful to display name,  
	    //number of entries, mean, rms, integral, overflow and underflow                                     
	    gStyle->SetOptFit(1);
	    
	    hjpsinum->Sumw2();
 	    TLatex tl; TString cname;
	    TString Compame = "DataMC";
	    TCanvas *c2 = new TCanvas();
	    
	    if(irap==0){
	      hjpsiMid->GetYaxis()->SetTitle("dN/dp_{T}");
	      hjpsiMid->GetXaxis()->SetTitle((idep==1) ? "Centrality bin" : "p_{T}");
	      hjpsiMid->Scale(1/hjpsiMid->Integral());
	      hjpsinum->Scale(1/hjpsinum->Integral());
	      hjpsinum->SetMarkerColor(kRed);
	      hjpsinum->SetLineColor(kRed);
	      hjpsiMid->Draw();
	      hjpsinum->Draw("E1 same");
	    }
	    if(irap==1){
	      hjpsiFwd->GetYaxis()->SetTitle("dN/dp_{T}");
	      hjpsiFwd->GetXaxis()->SetTitle((idep==1) ? "Centrality bin" : "p_{T}");
	      hjpsiFwd->Scale(1/hjpsiFwd->Integral());
	      hjpsinum->Scale(1/hjpsinum->Integral());
	      hjpsiFwd->Draw();
	      hjpsinum->SetMarkerColor(kRed);
	      hjpsinum->SetLineColor(kRed);
	      hjpsinum->Draw("E1 same");
	    }
	    
	    cname = "files/" + Compame + colltag + "_"  + objtag + "_" + deptag + "_" + raptag + "_" + cuttag;
	    double yshift =0; if (icoll==1 && irap==1) yshift=0.3;
	    if (icut>=3) yshift=0;
	    TLegend *tleg = new TLegend(0.61,0.667,0.92,0.78);
	    tleg->SetBorderSize(0);
	    if(irap==0)
	      tleg->AddEntry(hjpsiMid,Form("Data: %s", iobj==0 ? "J/#psi" : "#psi(2S)"),"lp");
	    if(irap==1)
	      tleg->AddEntry(hjpsiFwd,Form("Data: %s", iobj==0 ? "J/#psi" : "#psi(2S)"),"lp");
	    tleg->AddEntry(hjpsinum,Form("MC: %s", iobj==0 ? "J/#psi" : "#psi(2S)"),"lp");
	    tleg->Draw();
	    
	    tl.DrawLatex((irap==0) ? 9.4 : 5.2,1.3, colltag + TString(", ") + ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") + ((icut==0) ? "no #font[12]{l}_{J/#psi}^{3D} cut" : ((icut==1 || icut==3) ? "cut #font[12]{l}_{J/#psi}^{3D} cut" : "pt-dep #font[12]{l}_{J/#psi}^{3D} cut")));	    
	    c2->SaveAs(cname + ".root");
	    c2->SaveAs(cname + ".png");
	    c2->SaveAs(cname + ".pdf");
	    
	    TH1F *hjpsiDataPbPb;
	    if(irap==0) hjpsiDataPbPb = (TH1F*)hjpsiMid->Clone();
	    if(irap==1) hjpsiDataPbPb = (TH1F*)hjpsiFwd->Clone();
	    hjpsiDataPbPb->SetDirectory(0);
	    hjpsiDataPbPb->SetName("hjpsiDataPbPb");
	    hjpsiDataPbPb->Divide(hjpsinum);
	    hjpsiDataPbPb->GetYaxis()->SetTitle("DATA/MC");
	    
	    TAxis *YaxRatio = hjpsiDataPbPb->GetYaxis();
	    YaxRatio->SetLimits(0.4,1.4);
	    hjpsiDataPbPb->GetYaxis()->SetRangeUser(0.4,1.4);
	    
	    TCanvas *c3 = new TCanvas();
	    hjpsiDataPbPb->Draw();
	    
	    TString CompameR = "DataMC_Ratio";
	    cname = "files/" + CompameR + colltag + "_"  + objtag  + "_" + deptag + "_" + raptag + "_" + cuttag;
	    yshift =0; if (icoll==1 && irap==1) yshift=0.3;
	    if (icut>=3) yshift=0;
	    TLegend *tlegR = new TLegend(0.6,0.667,0.92,0.78);
	    tlegR->SetBorderSize(0);
	    if(irap==0)
	      tlegR->AddEntry(hjpsiMid,Form("Data/MC Ratio: %s", iobj==0 ? "J/#psi" : "#psi(2S)"),"p");
	    if(irap==1)
	      tlegR->AddEntry(hjpsiFwd,Form("Data/MC Ratio: %s", iobj==0 ? "J/#psi" : "#psi(2S)"),"p");
	    tlegR->Draw();
	    
	    tl.DrawLatex((irap==0) ? 9.4 : 5.2,1.3, colltag + TString(", ") + ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") + ((icut==0) ? "no #font[12]{l}_{J/#psi}^{3D} cut" : ((icut==1 || icut==3) ? "cut #font[12]{l}_{J/#psi}^{3D} cut" : "pt-dep #font[12]{l}_{J/#psi}^{3D} cut")));
	    TF1 *FitFn = new TF1("FitFn",(irap==0) ? Pol2 : Pol1,0.0,30.0,(irap==0) ? 3 : 2); 
	    FitFn->SetLineWidth(2.0);
	    FitFn->SetLineColor(2);
	    hjpsiDataPbPb->Fit("FitFn","", "",(irap==0) ? 6.5 : 3.0,30.0);
	    FitFn->Draw("same");
	    
	    c3->SaveAs(cname + ".root");
	    c3->SaveAs(cname + ".png");
	    c3->SaveAs(cname + ".pdf");
	    
	    CompameR = "Ratio";
	    cname = "files/" + CompameR + colltag + "_"  + objtag + "_" + deptag + "_" + raptag + "_" + cuttag;
	    
	    TFile *outf = new TFile(cname + ".root","RECREATE");
	    hjpsiDataPbPb->Write();
	    FitFn->Write();
	    outf->Close();
	    
	    //delete c1;  
	    //delete c2;  
	    //delete c3;  
	    //delete hjpsiDataPbPb;
	    //delete hjpsinum;    
	    //delete FitFn;
	    //delete hjpsiMid;
	    //delete hjpsiFwd;
	    
	  } // icut loop (without / with ctau cut)
	} // irap loop (mid / fwd)
      } // idep loop (pt / centrality)
    }
  } // icoll loop (pp / pbpb)
  
}

void setErr(TH1F *hist) {
  int nbins = hist->GetNbinsX();
  for (int i=1; i<nbins+1; i++) {
    hist->SetBinError(i,sqrt(hist->GetBinContent(i)));
  }
}

void setErrWt(TH1F *h1, TH1F *h2) {
  // for details see https://root.cern.ch/root/html532/src/TH1.cxx.html#2643
  Double_t b1,b2,w,ErrWt;
  int nbins = h1->GetNbinsX();
  for (int i=1; i<nbins+1; i++) {
    
    b1  = h1->GetBinContent(i);
    b2  = h2->GetBinContent(i);
    w = b1/b2;
    
    Double_t e1 = h1->GetBinError(i);
    Double_t e2 = h2->GetBinError(i);
    ErrWt = TMath::Abs( ( (1-2*w)*e1*e1 + w*w*e2*e2 )/(b2*b2) );
    h1->SetBinError(i,ErrWt);
  }
}

void fixCentPp(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   float y = hist->GetBinContent(1);
   float dy = hist->GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist->SetBinContent(i,y);
      hist->SetBinError(i,dy);
   }
}

TH1F* integrateHist(TH1F *hist) {
   TString name = hist->GetName(); name = name + "_int";
   TString title = hist->GetTitle(); title = title + " integrated";
   double integral, integralerror;
   int nbins = hist->GetNbinsX();
   integral = hist->IntegralAndError(1,nbins,integralerror);
   TH1F *ans = new TH1F(name, title, 1, hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinUpEdge(nbins));
   ans->SetBinContent(1,integral);
   ans->SetBinError(1,integralerror);
   return ans;
}
