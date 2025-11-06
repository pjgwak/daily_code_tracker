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

void plotCompare() {
  TFile *fjpsi_pp = new TFile("files/histos_jpsi_pp.root");
  TFile *fpsi2s_pp = new TFile("files/histos_psi2s_pp.root");
  TFile *fnpjpsi_pp = new TFile("files/histos_npjpsi_pp.root");
  TFile *fjpsi_pbpb = new TFile("files/histos_jpsi_pbpb.root");
  TFile *fpsi2s_pbpb = new TFile("files/histos_psi2s_pbpb.root");
  TFile *fnpjpsi_pbpb = new TFile("files/histos_npjpsi_pbpb.root");

  TFile *fjpsi_pp_Wt = new TFile("filesDataMC/histos_jpsi_pp.root");
  TFile *fpsi2s_pp_Wt  = new TFile("filesDataMC/histos_psi2s_pp.root");
  TFile *fnpjpsi_pp_Wt = new TFile("filesDataMC/histos_npjpsi_pp.root");
  TFile *fjpsi_pbpb_Wt = new TFile("filesDataMC/histos_jpsi_pbpb.root");
  TFile *fpsi2s_pbpb_Wt = new TFile("filesDataMC/histos_psi2s_pbpb.root");
  TFile *fnpjpsi_pbpb_Wt = new TFile("filesDataMC/histos_npjpsi_pbpb.root");
  
  // compare two efficiencies, one nominal and other weighted MC efficiency
  
  TString colltag, deptag, raptag, cuttag;
  for (int icoll=0; icoll<1; icoll++) {
    colltag = (icoll==0) ? "pp" : "pbpb";
    
    TString name = "files/histos_jpsi_" + colltag + ".root";
    TFile *fjpsi = new TFile(name);
    name = "files/histos_psi2s_" + colltag + ".root";
    TFile *fpsi2s = new TFile(name);
    name = "files/histos_npjpsi_" + colltag + ".root";
    TFile *fnpjpsi = new TFile(name);

    TString name_Wt = "filesDataMC/histos_jpsi_" + colltag + ".root";
    TFile *fjpsi_Wt  = new TFile(name);
    name = "filesDataMC/histos_psi2s_" + colltag + ".root";
    TFile *fpsi2s_Wt = new TFile(name);
    name = "filesDataMC/histos_npjpsi_" + colltag + ".root";
    TFile *fnpjpsi_Wt = new TFile(name);
    
    for (int idep=1; idep<2; idep++) { // idep=2 -> integrated                                                                              
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
	  TString hname_Wt;
	  
	  TCanvas *c1 = new TCanvas();
	  TLatex tl; TString cname;
	  TString effname = (icut<3) ? "singleff_" : "ctaucuteff_";
	  cname = "files/Comp/" + effname + colltag + "_" + deptag + "_" + raptag + "_" + cuttag;
	  
	  gStyle->SetOptTitle(0); // at least most of the time                                                                                
          gStyle->SetOptStat(0); // most of the time, sometimes "nemriou" might be useful to display name,                                  
          //number of entries, mean, rms, integral, overflow and underflow                                                                   
          gStyle->SetOptFit(0);

	  // now, let's draw simple ratios of efficiencies: psi(2S)/J/psi
	  // but do it only once
	  //delete c1;
	  //delete tg_jpsi; delete tg_psi2s; delete tg_npjpsi;
	  //continue;
	  
	  hname = "hnum" + cuttag + deptag + raptag;
	  TH1F *hjpsipp = (TH1F*) fjpsi_pp->Get(hname); hjpsipp->SetName(TString(hjpsipp->GetName()) + "_jpsipp");
	  TH1F *hpsi2spp = (TH1F*) fpsi2s_pp->Get(hname); hpsi2spp->SetName(TString(hpsi2spp->GetName()) + "_psi2spp");
	  TH1F *hjpsipbpb = (TH1F*) fjpsi_pbpb->Get(hname); hjpsipbpb->SetName(TString(hjpsipbpb->GetName()) + "_jpsipbpb");
	  TH1F *hpsi2spbpb = (TH1F*) fpsi2s_pbpb->Get(hname); hpsi2spbpb->SetName(TString(hpsi2spbpb->GetName()) + "_psi2spbpb");
	  hname = "hden_" + deptag + raptag;
	  TH1F *hjpsidenpp = (TH1F*) fjpsi_pp->Get(hname); hjpsidenpp->SetName(TString(hjpsidenpp->GetName()) + "_jpsidenpp");
	  TH1F *hpsi2sdenpp = (TH1F*) fpsi2s_pp->Get(hname); hpsi2sdenpp->SetName(TString(hpsi2sdenpp->GetName()) + "_psi2sdenpp");
	  TH1F *hjpsidenpbpb = (TH1F*) fjpsi_pbpb->Get(hname); hjpsipbpb->SetName(TString(hjpsidenpbpb->GetName()) + "_jpsidenpbpb");
	  TH1F *hpsi2sdenpbpb = (TH1F*) fpsi2s_pbpb->Get(hname); hpsi2sdenpbpb->SetName(TString(hpsi2sdenpbpb->GetName()) + "_psi2sdenpbpb");

	  hname_Wt = "hnum" + cuttag + deptag + raptag;
	  TH1F *hjpsipp_Wt = (TH1F*) fjpsi_pp_Wt->Get(hname_Wt); hjpsipp_Wt->SetName(TString(hjpsipp_Wt->GetName()) + "_jpsipp_Wt");
	  TH1F *hpsi2spp_Wt = (TH1F*) fpsi2s_pp_Wt->Get(hname_Wt); hpsi2spp_Wt->SetName(TString(hpsi2spp_Wt->GetName()) + "_psi2spp_Wt");
	  TH1F *hjpsipbpb_Wt = (TH1F*) fjpsi_pbpb_Wt->Get(hname_Wt); hjpsipbpb_Wt->SetName(TString(hjpsipbpb_Wt->GetName()) + "_jpsipbpb_Wt");
	  TH1F *hpsi2spbpb_Wt = (TH1F*) fpsi2s_pbpb_Wt->Get(hname_Wt); hpsi2spbpb_Wt->SetName(TString(hpsi2spbpb_Wt->GetName()) + "_psi2spbpb_Wt");
	  hname_Wt = "hden_" + deptag + raptag;
	  TH1F *hjpsidenpp_Wt = (TH1F*) fjpsi_pp_Wt->Get(hname_Wt); hjpsidenpp_Wt->SetName(TString(hjpsidenpp_Wt->GetName()) + "_jpsidenpp_Wt");
	  TH1F *hpsi2sdenpp_Wt = (TH1F*) fpsi2s_pp_Wt->Get(hname_Wt); hpsi2sdenpp_Wt->SetName(TString(hpsi2sdenpp_Wt->GetName()) + "_psi2sdenpp_Wt");
	  TH1F *hjpsidenpbpb_Wt = (TH1F*) fjpsi_pbpb_Wt->Get(hname_Wt); hjpsipbpb_Wt->SetName(TString(hjpsidenpbpb_Wt->GetName()) + "_jpsidenpbpb_Wt");
	  TH1F *hpsi2sdenpbpb_Wt = (TH1F*) fpsi2s_pbpb_Wt->Get(hname_Wt); hpsi2sdenpbpb_Wt->SetName(TString(hpsi2sdenpbpb_Wt->GetName()) + "_psi2sdenpbpb_Wt");
	  if (idep==2) {
	    hjpsipp = integrateHist(hjpsipp);
	    hpsi2spp = integrateHist(hpsi2spp);
	    hjpsipbpb = integrateHist(hjpsipbpb);
	    hpsi2spbpb = integrateHist(hpsi2spbpb);
	    hjpsidenpp = integrateHist(hjpsidenpp);
	    hpsi2sdenpp = integrateHist(hpsi2sdenpp);
	    hjpsidenpbpb = integrateHist(hjpsidenpbpb);
	    hpsi2sdenpbpb = integrateHist(hpsi2sdenpbpb);
	  }
	  
	  if (idep==1) { // centrality for pp... fill all the bins
	    fixCentPp(hjpsipp);
	    fixCentPp(hpsi2spp);
	    fixCentPp(hjpsidenpp);
	    fixCentPp(hpsi2sdenpp);
	  }
	  
	  if (idep==2) {
	    hjpsipp_Wt = integrateHist(hjpsipp_Wt);
	    hpsi2spp_Wt = integrateHist(hpsi2spp_Wt);
	    hjpsipbpb_Wt = integrateHist(hjpsipbpb_Wt);
	    hpsi2spbpb_Wt = integrateHist(hpsi2spbpb_Wt);
	    hjpsidenpp_Wt = integrateHist(hjpsidenpp_Wt);
	    hpsi2sdenpp_Wt = integrateHist(hpsi2sdenpp_Wt);
	    hjpsidenpbpb_Wt = integrateHist(hjpsidenpbpb_Wt);
	    hpsi2sdenpbpb_Wt = integrateHist(hpsi2sdenpbpb_Wt);
	  }
	  
	  if (idep==1) { // centrality for pp_Wt... fill all the bins
	    fixCentPp(hjpsipp_Wt);
	    fixCentPp(hpsi2spp_Wt);
	    fixCentPp(hjpsidenpp_Wt);
	    fixCentPp(hpsi2sdenpp_Wt);
	  }
	  
	  hjpsipp->Sumw2(true); hjpsidenpp->Sumw2(true);
	  hpsi2spp->Sumw2(true); hpsi2sdenpp->Sumw2(true);
	  hjpsipp->Divide(hjpsipp,hjpsidenpp,1,1,"B");
	  hpsi2spp->Divide(hpsi2spp,hpsi2sdenpp,1,1,"B");
	  hjpsipbpb->Divide(hjpsipbpb,hjpsidenpbpb,1,1,"B");
	  hpsi2spbpb->Divide(hpsi2spbpb,hpsi2sdenpbpb,1,1,"B");
	  
	  hjpsipp_Wt->Sumw2(true); hjpsidenpp_Wt->Sumw2(true);
	  hpsi2spp_Wt->Sumw2(true); hpsi2sdenpp_Wt->Sumw2(true);
	  
	  //TH1F *hjpsippDp_Wt = (TH1F*)hjpsipp_Wt->Clone();
          //hjpsippDp_Wt->SetDirectory(0);
          //hjpsippDp_Wt->SetName("hjpsippDp_Wt");
	  
	  hjpsipp_Wt->Divide(hjpsipp_Wt,hjpsidenpp_Wt,1,1,"B");
	  hpsi2spp_Wt->Divide(hpsi2spp_Wt,hpsi2sdenpp_Wt,1,1,"B");
	  hjpsipbpb_Wt->Divide(hjpsipbpb_Wt,hjpsidenpbpb_Wt,1,1,"B");
	  hpsi2spbpb_Wt->Divide(hpsi2spbpb_Wt,hpsi2sdenpbpb_Wt,1,1,"B");
	  
	  hpsi2spp->Divide(hjpsipp);
	  hpsi2spbpb->Divide(hjpsipbpb);
	  hpsi2spp_Wt->Divide(hjpsipp_Wt);
	  hpsi2spbpb_Wt->Divide(hjpsipbpb_Wt);
	  
	  hpsi2spp->SetMarkerColor(kBlack);
	  hpsi2spp->SetLineColor(kBlack);
	  hpsi2spp_Wt->SetMarkerColor(kRed);
	  hpsi2spp_Wt->SetLineColor(kRed);
	  
	  hpsi2spbpb->SetMarkerColor(kBlue);
	  hpsi2spbpb->SetLineColor(kBlue);
	  hpsi2spbpb_Wt->SetMarkerColor(kGreen);
	  hpsi2spbpb_Wt->SetLineColor(kGreen);
	  
	  //hpsi2spp->Divide(hpsi2spp_Wt);
	  //hpsi2spbpb->Divide(hpsi2spbpb_Wt);
	  
	  cname = "files/Comp/simpleratio_" + deptag + "_" + raptag + "_" + cuttag;
	  
	  TCanvas *c2 = new TCanvas();
	  if (idep<2) {
	    // plot
	    TH1F *haxes = new TH1F("haxes","haxes",1,0,(idep==1) ? 100 : 30);
	    haxes->GetYaxis()->SetTitle("Eff(#psi(2S)) / Eff(J/#psi)");
            haxes->GetYaxis()->SetRangeUser(0.5,1.5);
            haxes->SetBinContent(1,1);
            haxes->Draw();
	    haxes->GetXaxis()->SetTitle((idep==1) ? "Centrality bin" : "p_{T}");
	    hpsi2spp->Draw("E1 same");
	    hpsi2spp_Wt->Draw("E1 same");
	    hpsi2spbpb->Draw("E1 same");
	    hpsi2spbpb_Wt->Draw("E1 same");
	    
	    TLegend *tleg2 = new TLegend(0.7,0.17,0.89,0.31);
	    tleg2->SetBorderSize(0);
	    tleg2->AddEntry(hpsi2spp,"pp ratio","lp");
	    tleg2->AddEntry(hpsi2spp_Wt,"pp ratio weighted","lp");
	    tleg2->AddEntry(hpsi2spbpb,"pbpb ratio","lp");
	    tleg2->AddEntry(hpsi2spbpb_Wt,"pbpb ratio weighted","lp");
	    tleg2->Draw();
	    
	    tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
			 + ((icut==0) ? "no #font[12]{l}_{J/#psi}^{3D} cut" : ((icut==1) ? "cut #font[12]{l}_{J/#psi}^{3D} cut" : "pt-dep #font[12]{l}_{J/#psi}^{3D} cut")));
	    
	    c2->SaveAs(cname + ".root");
	    c2->SaveAs(cname + ".png");
	    c2->SaveAs(cname + ".pdf");

	    delete haxes;
	  }
	  
	  // at last, the double ratio
	  hpsi2spbpb->Divide(hpsi2spp);
	  hpsi2spbpb_Wt->Divide(hpsi2spp_Wt);
	  
	  hpsi2spbpb->SetMarkerColor(kBlack);
	  hpsi2spbpb->SetLineColor(kBlack);
	  hpsi2spbpb_Wt->SetMarkerColor(kRed);
	  hpsi2spbpb_Wt->SetLineColor(kRed);
          
	  cname = "files/Comp/doubleratio_" + deptag + "_" + raptag + "_" + cuttag;
	  TCanvas *c3 = new TCanvas();
	  	  
	  if (idep<2) {
	    // plot
	    TH1F *haxes = new TH1F("haxes","haxes",1,0,(idep==1) ? 100 : 30);
            haxes->GetYaxis()->SetTitle("[Eff(#psi(2S)) / Eff(J/#psi)]_{PbPb} / [Eff(#psi(2S)) / Eff(J/#psi)]_{pp}");
            haxes->GetYaxis()->SetRangeUser(0.5,1.5);
            haxes->SetBinContent(1,1);
	    haxes->GetXaxis()->SetTitle((idep==1) ? "Centrality bin" : "p_{T}");
	    haxes->GetYaxis()->SetTitleSize(0.03);
	    haxes->GetYaxis()->SetTitleOffset(2);
	    haxes->Draw();
	    hpsi2spbpb->Draw("E1 same");
	    hpsi2spbpb_Wt->Draw("E1 same");

	    TLegend *tleg4 = new TLegend(0.7,0.17,0.89,0.31);
	    tleg4->SetBorderSize(0);
	    tleg4->AddEntry(hpsi2spbpb,"DR","lp");
	    tleg4->AddEntry(hpsi2spbpb_Wt,"DR weighted","lp");
	    tleg4->Draw();
	    
	    tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
			 + ((icut==0) ? "no #font[12]{l}_{J/#psi}^{3D} cut" : ((icut==1) ? "cut #font[12]{l}_{J/#psi}^{3D} cut" : "pt-dep #font[12]{l}_{J/#psi}^{3D} cut")));
	    c3->SaveAs(cname + ".root");
	    c3->SaveAs(cname + ".png");
	    c3->SaveAs(cname + ".pdf");
	    
	    delete haxes;
	    
	  }
	  
	  // clean behind ourselves
	  //delete c1;
	  //delete tg_jpsi; delete tg_psi2s; delete tg_npjpsi;
	  //delete hjpsinum;	       
	  //delete hpsi2snum; 	    
	  //delete hnpjpsinum;
	  //delete hjpsiden;
	  //delete hpsi2sden; 
	  //delete hnpjpsiden;	       
	  //delete hname;
	  
	} // icut loop (without / with ctau cut)
      } // irap loop (mid / fwd)
    } // idep loop (pt / centrality)
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
