#include "../Macros/CMS/CMS_lumi.C"
#include "../Macros/CMS/tdrstyle.C"
#include "datasheet.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>



using namespace std;

const bool isPreliminary = false;

void SetLegendStyle(TLegend &leg) ;
int SetColor(const int idx, const bool bkg=false) ;
int SetMarker(const int i) ;
void SetGraphStyle(TGraph *g, const int idx, const bool bkg, const int i) ;
void SetGraphStyle(TGraph *g, const int idx, const bool bkg, const double alphaV, const int i);

void comparePsi2SJpsi(
    const int compare=0, //0: vs. cent, 1: vs. pT
    const int drawopt=0 //0: 2 graphs each in a same region, 1: all 4 graphs
    ) ;

void compareYnSJpsi(
    const bool doprompt,
    const int compare, //0: vs. cent, 1: vs. pT
    const int drawopt //0: Y(1,2,3S), J/psi, 1: Y(1,2S), J/psi, psi(2S), 2: Y(1), J/psi, psi(2S)
    ) ;

void extendLowPtFor14005(
    const bool doprompt
    ) ;



//////
void compareRAA(
    const bool doprompt=false,
    const int compare=0, //0: all, 1: psi2S+D meson, 2: all except B
    const bool v2=false, //true: extended x-axis range
    const bool logx=false
    ) {

  setTDRStyle();

  const int jpsicolor = 0;
  const int jpsimarker = 0;
  const int psi2scolor = 6;
  const int psi2smarker = 1;
  const int Dcolor=2;
  const int Dmarker=1;
  const int Bcolor=8;
  const int Bmarker=1;
  const int Chgmarker=1;

  // Other results
  TGraphAsymmErrors *DmesonRAA_cent_0_100_1_syst = new TGraphAsymmErrors(5,DmesonRAA_fx3001,DmesonRAA_fy3001,DmesonRAA_felx3001,DmesonRAA_fehx3001,DmesonRAA_fely3001,DmesonRAA_fehy3001);
  TGraphAsymmErrors *DmesonRAA_cent_0_100_1 = new TGraphAsymmErrors(5,DmesonRAA_fx3001,DmesonRAA_fy3001,DmesonRAA_felx3001,DmesonRAA_fehx3001,DmesonRAA_fey3001,DmesonRAA_fey3001);
  DmesonRAA_cent_0_100_1->SetName("DmesonRAA_cent_0_100_1");
  DmesonRAA_cent_0_100_1_syst->SetName("DmesonRAA_cent_0_100_1_syst");
  
  TGraphAsymmErrors *DmesonRAA_cent_0_100_2_syst = new TGraphAsymmErrors(9,DmesonRAA_fx3002,DmesonRAA_fy3002,DmesonRAA_felx3002,DmesonRAA_fehx3002,DmesonRAA_fely3002,DmesonRAA_fehy3002);
  TGraphAsymmErrors *DmesonRAA_cent_0_100_2 = new TGraphAsymmErrors(9,DmesonRAA_fx3002,DmesonRAA_fy3002,DmesonRAA_felx3002,DmesonRAA_fehx3002,DmesonRAA_fey3002,DmesonRAA_fey3002);
  DmesonRAA_cent_0_100_2->SetName("DmesonRAA_cent_0_100_2");
  DmesonRAA_cent_0_100_2_syst->SetName("DmesonRAA_cent_0_100_2_syst");
  
  SetGraphStyle(DmesonRAA_cent_0_100_1, Dcolor, 0, Dmarker);
  SetGraphStyle(DmesonRAA_cent_0_100_1_syst, Dcolor, 1, Dmarker);
  SetGraphStyle(DmesonRAA_cent_0_100_2, Dcolor, 0, Dmarker);
  SetGraphStyle(DmesonRAA_cent_0_100_2_syst, Dcolor, 1, Dmarker);

  TGraphAsymmErrors *BmesonRAA_cent_0_100_syst = new TGraphAsymmErrors(5,BmesonRAA_fx3002,BmesonRAA_fy3002,BmesonRAA_felx3002,BmesonRAA_fehx3002,BmesonRAA_fely3002,BmesonRAA_fehy3002);
  TGraphAsymmErrors *BmesonRAA_cent_0_100 = new TGraphAsymmErrors(5,BmesonRAA_fx3002,BmesonRAA_fy3002,BmesonRAA_felx3002,BmesonRAA_fehx3002,BmesonRAA_fey3002,BmesonRAA_fey3002);
  BmesonRAA_cent_0_100->SetName("BmesonRAA_cent_0_100");
  BmesonRAA_cent_0_100_syst->SetName("BmesonRAA_cent_0_100_syst");

  SetGraphStyle(BmesonRAA_cent_0_100, Bcolor, 0, Bmarker);
  SetGraphStyle(BmesonRAA_cent_0_100_syst, Bcolor, 1, Bmarker);

  TH1D *hChgRAA_cent_0_100 = new TH1D("h_ChgRAA_cent_0_100",";p_{T} [GeV/c];R_{AA}",37,chgRAA_xAxis);
  TH1D *hChgRAA_cent_0_100_syst = new TH1D("h_ChgRAA_cent_0_100_syst",";p_{T} [GeV/c];R_{AA}",37,chgRAA_xAxis);
  for (int i=0; i<37; i++) {
    hChgRAA_cent_0_100->SetBinContent(i+1,chgRAA_yAxis[i]);
    hChgRAA_cent_0_100->SetBinError(i+1,chgRAA_yError[i]);
    hChgRAA_cent_0_100_syst->SetBinContent(i+1,chgRAA_yAxis[i]);
    hChgRAA_cent_0_100_syst->SetBinError(i+1,chgRAA_yErrorSyst[i]);
  }
  TGraphErrors *ChgRAA_cent_0_100 = new TGraphErrors(hChgRAA_cent_0_100);
  TGraphErrors *ChgRAA_cent_0_100_syst = new TGraphErrors(hChgRAA_cent_0_100_syst);
  ChgRAA_cent_0_100->SetName("ChgRAA_cent_0_100");
  ChgRAA_cent_0_100_syst->SetName("ChgRAA_cent_0_100_syst");

  ChgRAA_cent_0_100->SetMarkerStyle(SetMarker(Chgmarker));
  ChgRAA_cent_0_100->SetMarkerColor(kBlack);
  ChgRAA_cent_0_100->SetLineColor(kBlack);
  ChgRAA_cent_0_100->SetMarkerSize(1.1);
  ChgRAA_cent_0_100_syst->SetMarkerStyle(SetMarker(Chgmarker));
  ChgRAA_cent_0_100_syst->SetMarkerColor(kBlack);
  ChgRAA_cent_0_100_syst->SetLineColor(kBlack);
  ChgRAA_cent_0_100_syst->SetMarkerSize(1.1);
  ChgRAA_cent_0_100_syst->SetFillColorAlpha(kBlack, 0.15);

  // 16025 results
  string pr_fname[3] = {
    "./root/result_JPsi_RAA_pt_prompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_prompt_1_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_prompt_2_linearX_accEffCorr.root"
  };
  string np_fname[3] = {
    "./root/result_JPsi_RAA_pt_nonprompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_nonprompt_1_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_nonprompt_2_linearX_accEffCorr.root"
  };
  
  TFile inputf( doprompt ? pr_fname[2].c_str() : np_fname[2].c_str() );
  TFile inputf0( doprompt ? pr_fname[0].c_str() : np_fname[0].c_str() ); // to grab the lowest 3 pt bins when drawing integrated RAA vs. pt

  // Load histograms
  TCanvas *c1 = (TCanvas*)inputf.Get("c1");
  TGraphAsymmErrors *bin_0 = (TGraphAsymmErrors *)c1->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst = (TGraphAsymmErrors *)c1->FindObject("bin_0_syst");
  bin_0->SetName("bin_0");
  bin_0_syst->SetName("bin_0_syst");
  TH1F *haxes = (TH1F*)c1->FindObject("haxes");

  TCanvas *c0 = (TCanvas*)inputf0.Get("c1");
  TGraphAsymmErrors *bin_lowpt_ = (TGraphAsymmErrors *)c0->FindObject("bin_3");
  TGraphAsymmErrors *bin_lowpt_syst_ = (TGraphAsymmErrors *)c0->FindObject("bin_3_syst");

  int npoint_lowpt = bin_lowpt_->GetN();
  Double_t *x_lowpt = bin_lowpt_->GetX();
  Double_t *y_lowpt = bin_lowpt_->GetY();
  Double_t *exlow_lowpt = bin_lowpt_->GetEXlow();
  Double_t *eylow_lowpt = bin_lowpt_->GetEYlow();
  Double_t *exhigh_lowpt = bin_lowpt_->GetEXhigh();
  Double_t *eyhigh_lowpt = bin_lowpt_->GetEYhigh();
  Double_t *eylow_syst_lowpt = bin_lowpt_syst_->GetEYlow();
  Double_t *eyhigh_syst_lowpt = bin_lowpt_syst_->GetEYhigh();

  Double_t g_x_lowpt[3], g_y_lowpt[3], g_exl_lowpt[3], g_exh_lowpt[3], g_eyl_lowpt[3], g_eyh_lowpt[3], g_eyl_syst_lowpt[3], g_eyh_syst_lowpt[3];
  for (int i=0; i<3; i++) {
    g_x_lowpt[i] = x_lowpt[i];
    g_y_lowpt[i] = y_lowpt[i];
    g_exl_lowpt[i] = exlow_lowpt[i];
    g_exh_lowpt[i] = exhigh_lowpt[i];
    g_eyl_lowpt[i] = eylow_lowpt[i];
    g_eyh_lowpt[i] = eyhigh_lowpt[i];
    g_eyl_syst_lowpt[i] = eylow_syst_lowpt[i];
    g_eyh_syst_lowpt[i] = eyhigh_syst_lowpt[i];
  }
  
  TGraphAsymmErrors *bin_lowpt = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_lowpt,g_eyh_lowpt);
  TGraphAsymmErrors *bin_lowpt_syst = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_syst_lowpt,g_eyh_syst_lowpt);
  bin_lowpt->SetName("bin_lowpt");
  bin_lowpt_syst->SetName("bin_lowpt_syst");
  SetGraphStyle(bin_lowpt, jpsicolor, 0, 4);
  SetGraphStyle(bin_lowpt_syst, jpsicolor, 1, 4);
  
  // 16004 psi2S
  string psi2s_fname[2] = {
    "./root16004/result_Psi2S_RAA_cent_prompt_16004_3_accEffCorr.root",
    "./root16004/result_Psi2S_RAA_pt_prompt_16004_3_linearX_accEffCorr.root"
  };
  
  TFile inputf_psi2s( psi2s_fname[1].c_str() );

  TCanvas *c1_psi2s = (TCanvas*)inputf_psi2s.Get("c1");
  TGraphAsymmErrors *bin_0_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0_syst");
  bin_0_psi2s->SetName("bin_0_psi2s");
  bin_0_syst_psi2s->SetName("bin_0_syst_psi2s");
  SetGraphStyle(bin_0_psi2s, psi2scolor, 0, psi2smarker);
  SetGraphStyle(bin_0_syst_psi2s, psi2scolor, 1, psi2smarker);
  
  TGraphAsymmErrors *bin_1_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1");
  TGraphAsymmErrors *bin_1_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1_syst");
  bin_1_psi2s->SetName("bin_1_psi2s");
  bin_1_syst_psi2s->SetName("bin_1_syst_psi2s");
  SetGraphStyle(bin_1_psi2s, psi2scolor+1, 0, psi2smarker);
  SetGraphStyle(bin_1_syst_psi2s, psi2scolor+1, 1, psi2smarker);


  double xpos1 = haxes->GetXaxis()->GetBinLowEdge(1);
  int lastbin = haxes->GetNbinsX();
  double xpos2 = v2 ? 400 : 60;
//  double xpos2 = haxes->GetXaxis()->GetBinLowEdge(lastbin) + haxes->GetXaxis()->GetBinWidth(lastbin);
  double xaxis[2] = {0};
  xaxis[0] = logx ? 1 : xpos1;
  xaxis[1] = xpos2;

  TH1F *gaxes = new TH1F("gaxes","",1,xaxis);
  gaxes->GetXaxis()->CenterTitle();
  gaxes->GetXaxis()->SetTitle(haxes->GetXaxis()->GetTitle());
  if (logx) {
    gaxes->GetXaxis()->SetTitleSize(0.05);
    gaxes->GetXaxis()->SetTitleOffset(1.15);
    gaxes->GetXaxis()->SetLabelSize(0.04);
    gaxes->GetXaxis()->SetLabelOffset(0);
  }
  gaxes->GetYaxis()->SetTitle(haxes->GetYaxis()->GetTitle());
  gaxes->GetYaxis()->SetRangeUser(0, 1.7);

  TGraphAsymmErrors *bin_1, *bin_1_syst, *bin_2, *bin_2_syst, *bin_3, *bin_3_syst;

  // Draw graphs into 1 canvas
  gStyle->SetEndErrorSize(5);
  TCanvas canv("canv","canv",600,600);
  canv.SetLogx(logx);
  gaxes->Draw("");
  canv.Update();
  if (compare==1) { 
    DmesonRAA_cent_0_100_1_syst->Draw("5 same");
    DmesonRAA_cent_0_100_1->Draw("p same");
    DmesonRAA_cent_0_100_2_syst->Draw("5 same");
    DmesonRAA_cent_0_100_2->Draw("p same");
    bin_0_syst_psi2s->Draw("5 same");
    bin_0_psi2s->Draw("p same");
    bin_1_syst_psi2s->Draw("5 same");
    bin_1_psi2s->Draw("p same");
  }
  else {
    ChgRAA_cent_0_100_syst->Draw("5 same");
    ChgRAA_cent_0_100->Draw("p same");
    DmesonRAA_cent_0_100_1_syst->Draw("5 same");
    DmesonRAA_cent_0_100_1->Draw("p same");
    DmesonRAA_cent_0_100_2_syst->Draw("5 same");
    DmesonRAA_cent_0_100_2->Draw("p same");
    if (compare==0) {
      BmesonRAA_cent_0_100_syst->Draw("5 same");
      BmesonRAA_cent_0_100->Draw("p same");
    }
  }

  bin_0_syst->Draw("5 same");
  bin_0->Draw("p same");
  bin_lowpt_syst->Draw("5 same");
  bin_lowpt->Draw("p same");
  
  canv.Update();

  TLegend leg(0.22,0.71,0.55,0.89,NULL,"brNDC");
  SetLegendStyle(leg);
  leg.SetTextSize(0.033);
  leg.SetMargin(0.15);
  leg.SetNColumns(2);

  TLegend leg1(0.215,0.64,0.55,0.805,NULL,"brNDC");
  SetLegendStyle(leg1);
  leg1.SetTextSize(0.033);

  TLegend leg2(0.215,0.64,0.55,0.805,NULL,"brNDC");
  SetLegendStyle(leg2);
  leg2.SetTextSize(0.033);

  if (doprompt) {
    if (compare==1) { // psi2s+D
      leg.SetY1(0.71);
      leg.SetY2(0.89);
      leg1.SetY1(0.60);
      leg1.SetY2(0.70);
    } 
  } else {
    if (compare==0) { // all
      leg.SetY1(0.74);
      leg.SetY2(0.89);
      leg1.SetY1(0.62);
      leg1.SetY2(0.72);

      leg2.SetX1(0.45);
      leg2.SetX2(0.75);
      leg2.SetY1(0.62);
      leg2.SetY2(0.72);
    }
    else if (compare==2) { // B is not drawn
      leg.SetY1(0.77);
      leg.SetY2(0.89);
      leg1.SetY1(0.64);
      leg1.SetY2(0.74);

      leg2.SetX1(0.44);
      leg2.SetX2(0.75);
      leg2.SetY1(0.64);
      leg2.SetY2(0.74);
    }
  };

  TLegendEntry *entry;

  // vs. pT
  if (doprompt && compare==1) { // this is only for pr
    leg.SetHeader("Hidden charm");
    entry = leg.AddEntry((TObject*)0,"#scale[1.0]{Prompt J/#psi}","");
    entry = leg.AddEntry((TObject*)0,"","");
    entry = leg.AddEntry("bin_lowpt","#scale[0.83]{1.8 < |y| < 2.4}","p");
    entry = leg.AddEntry("bin_0","#scale[0.83]{|y| < 2.4}","p");
    entry = leg.AddEntry((TObject*)0,"#scale[1.0]{Prompt #psi(2S)}","");
    entry = leg.AddEntry((TObject*)0,"","");
    entry = leg.AddEntry("bin_1_psi2s","#scale[0.83]{1.6 < |y| < 2.4}","p");
    entry = leg.AddEntry("bin_0_psi2s","#scale[0.83]{|y| < 1.6}","p");
    leg1.SetHeader("Open charm");
    entry = leg1.AddEntry("DmesonRAA_cent_0_100_1","#scale[0.83]{D^{0}} #scale[0.67]{HIN-16-001}","p");
    entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 1}","");
  }
  else if (!doprompt && compare!=1) { // this is only for np
    leg.SetHeader("Open beauty");
    entry = leg.AddEntry((TObject*)0,"#scale[1.0]{J/#psi from b hadrons}","");
    entry = leg.AddEntry((TObject*)0,"","");
    entry = leg.AddEntry("bin_lowpt","#scale[0.83]{1.8 < |y| < 2.4}","p");
    entry = leg.AddEntry("bin_0","#scale[0.83]{|y| < 2.4}","p");
    if (compare==0) {
      entry = leg.AddEntry("BmesonRAA_cent_0_100","#scale[0.83]{B^{+}} #scale[0.75]{|y| < 2.4}, #scale[0.67]{arXiv:1705.04727}","p");
    }
    leg1.SetHeader("Open charm");
    entry = leg1.AddEntry("DmesonRAA_cent_0_100_1","#scale[0.83]{D^{0}} #scale[0.67]{HIN-16-001}","p");
    entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 1}","");
    leg2.SetHeader("Light flavor");
    entry = leg2.AddEntry("ChgRAA_cent_0_100","#scale[0.83]{Charged hadrons} #scale[0.67]{JHEP 04 (2017)}","p");
    leg2.AddEntry((TObject*)0,"#scale[0.75]{|#eta| < 1}","");
  }

  // TEXTS
  TLatex tl;
  tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.044);
  tl.SetTextSize(0.037);
  if (v2) tl.DrawLatex(0.70,0.18,"Cent. 0-100%");
  else tl.DrawLatex(0.64,0.53,"Cent. 0-100%");

  // x-axis is always vs. pT -> use 106 option
  CMS_lumi( (TPad*) gPad, 106, 33, "", isPreliminary); 

  // draw a line at RAA=1
  TLine line; line.SetLineColor(kBlack); line.SetNDC();
  line.DrawLine(xpos1,1,xpos2,1);
  canv.Update();
  leg.Draw("same");
  leg1.Draw("same");
  if (!doprompt) {
    leg2.Draw("same");
  }

  // Draw upper limits 
  TArrow arrow; arrow.SetNDC(); arrow.SetArrowSize(0.03); arrow.SetLineWidth(2);
  arrow.SetOption("<-|");
  if (doprompt && compare==1) {
    arrow.SetLineColor(SetColor(psi2scolor+1));
    arrow.DrawArrow(psi2S_pt_fwd_fx, 0.01, psi2S_pt_fwd_fx, psi2S_pt_fwd_fy);
  }

  // Global uncertainty boxes
  TBox *box0, *box_psi2s_mid, *box_psi2s_fwd;
  double dx= (logx ? 5 : 1.25);
  if (v2) dx = 40;
  box0 = new TBox(xpos2-dx, Jpsi_pt_rap0024_global[0], xpos2, Jpsi_pt_rap0024_global[1]);
  box0->SetLineColor(kRed+2);
  box0->SetFillColorAlpha(kRed,0.35);
  box0->Draw("lf");
  if (compare==1) {
    box_psi2s_mid = new TBox(xpos2-dx*2, psi2S_pt_mid_global[0], xpos2-dx, psi2S_pt_mid_global[1]);
    box_psi2s_mid->SetFillColorAlpha(SetColor(psi2scolor),0.7);
    box_psi2s_fwd = new TBox(xpos2-dx*3, psi2S_pt_fwd_global[0], xpos2-dx*2, psi2S_pt_fwd_global[1]);
    box_psi2s_fwd->SetFillColorAlpha(SetColor(psi2scolor+1),0.7);
    box_psi2s_mid->Draw("lf");
    box_psi2s_fwd->Draw("lf");
  }

  canv.RedrawAxis();

  // Save canvas
  string outputname = "RAA_pt_";
  outputname += Form("%s",doprompt?"pr_":"np_");
  outputname += "integrated_";
  if (compare==0) outputname += "all";
  else if (compare==1) outputname += "D";
  else if (compare==2) outputname += "excludeB";
  if (logx) outputname += Form("_log_%s.pdf", v2?"v2":"v1");
  else outputname += Form("_lin_%s.pdf", v2?"v2":"v1");

  cout << outputname << endl;
  canv.SaveAs(outputname.c_str());

  return ;
}

void comparePsi2SJpsi(
    const int compare, //0: vs. cent, 1: vs. pT
    const int drawopt //0: 2 graphs each in a same region, 1: all 4 graphs
    ) {

  setTDRStyle();

  // Set style indexes
  const int jpsicolor = 0;
  const int jpsimarker = 0;
  const int psi2scolor = 6;
  const int psi2smarker = 1;

  // 16004 results
  string jpsi_fname[2] = {
    "./root16004/result_JPsi_RAA_cent_prompt_16004_3_accEffCorr.root",
    "./root16004/result_JPsi_RAA_pt_prompt_16004_3_linearX_accEffCorr.root"
  };
  string psi2s_fname[2] = {
    "./root16004/result_Psi2S_RAA_cent_prompt_16004_3_accEffCorr.root",
    "./root16004/result_Psi2S_RAA_pt_prompt_16004_3_linearX_accEffCorr.root"
  };
  
  TFile inputf_jpsi( jpsi_fname[compare].c_str() );
  TFile inputf_psi2s( psi2s_fname[compare].c_str() );

  // Load histograms
  TCanvas *c1_jpsi = (TCanvas*)inputf_jpsi.Get("c1");
  TGraphAsymmErrors *bin_0_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0_syst");
  bin_0_jpsi->SetName("bin_0_jpsi");
  bin_0_syst_jpsi->SetName("bin_0_syst_jpsi");
  SetGraphStyle(bin_0_jpsi, jpsicolor, 0, jpsimarker);
  SetGraphStyle(bin_0_syst_jpsi, jpsicolor, 1, jpsimarker);

  TGraphAsymmErrors *bin_1_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_1");
  TGraphAsymmErrors *bin_1_syst_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_1_syst");
  bin_1_jpsi->SetName("bin_1_jpsi");
  bin_1_syst_jpsi->SetName("bin_1_syst_jpsi");
  SetGraphStyle(bin_1_jpsi, jpsicolor+drawopt, 0, jpsimarker);
  SetGraphStyle(bin_1_syst_jpsi, jpsicolor+drawopt, 1, jpsimarker);

  TCanvas *c1_psi2s = (TCanvas*)inputf_psi2s.Get("c1");
  TGraphAsymmErrors *bin_0_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0_syst");
  bin_0_psi2s->SetName("bin_0_psi2s");
  bin_0_syst_psi2s->SetName("bin_0_syst_psi2s");
  SetGraphStyle(bin_0_psi2s, psi2scolor, 0, psi2smarker);
  SetGraphStyle(bin_0_syst_psi2s, psi2scolor, 1, psi2smarker);
  
  TGraphAsymmErrors *bin_1_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1");
  TGraphAsymmErrors *bin_1_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1_syst");
  bin_1_psi2s->SetName("bin_1_psi2s");
  bin_1_syst_psi2s->SetName("bin_1_syst_psi2s");
  SetGraphStyle(bin_1_psi2s, psi2scolor+drawopt, 0, psi2smarker);
  SetGraphStyle(bin_1_syst_psi2s, psi2scolor+drawopt, 1, psi2smarker);

  // Dummy objects for drawing
  TH1F *haxes = (TH1F*)c1_psi2s->FindObject("haxes");
  TH1F *bin_0_jpsi_clone = (TH1F*)haxes->Clone("bin_0_jpsi_clone");
  bin_0_jpsi_clone->SetMarkerColor(kBlack);
  bin_0_jpsi_clone->SetMarkerStyle(SetMarker(jpsimarker));
  TH1F *bin_0_psi2s_clone = (TH1F*)haxes->Clone("bin_0_psi2s_clone");
  bin_0_psi2s_clone->SetMarkerColor(kBlack);
  bin_0_psi2s_clone->SetMarkerStyle(SetMarker(psi2smarker));
  
  if (SetMarker(jpsimarker) == kFullStar) bin_0_jpsi_clone->SetMarkerSize(2.3);
  else if (SetMarker(jpsimarker) == kFullDiamond) bin_0_jpsi_clone->SetMarkerSize(2.2);
  else if (SetMarker(jpsimarker) == kFullCross) bin_0_jpsi_clone->SetMarkerSize(2.0);
  else bin_0_jpsi_clone->SetMarkerSize(1.5);
  if (SetMarker(psi2smarker) == kFullStar) bin_0_psi2s_clone->SetMarkerSize(2.3);
  else if (SetMarker(psi2smarker) == kFullDiamond) bin_0_psi2s_clone->SetMarkerSize(2.2);
  else if (SetMarker(psi2smarker) == kFullCross) bin_0_psi2s_clone->SetMarkerSize(2.0);
  else bin_0_psi2s_clone->SetMarkerSize(1.5);

  // Draw graphs into 1 canvas
  gStyle->SetEndErrorSize(5);
  TCanvas canv("canv","canv",600,600);
  haxes->GetYaxis()->SetRangeUser(0,1.5);
  haxes->Draw("");
  TCanvas canv1("canv1","canv1",600,600); // for drawopt==0 case
  haxes->Draw("");
  
  if (drawopt==0) {
    canv.cd();
    bin_0_syst_jpsi->Draw("5 same");
    bin_0_jpsi->Draw("p same");
    bin_0_syst_psi2s->Draw("5 same");
    bin_0_psi2s->Draw("p same");
    canv.Update();
    
    canv1.cd();
    bin_1_syst_jpsi->Draw("5 same");
    bin_1_jpsi->Draw("p same");
    bin_1_syst_psi2s->Draw("5 same");
    bin_1_psi2s->Draw("p same");
    canv1.Update();
  } else if (drawopt==1) {
    canv.cd();
    bin_0_syst_jpsi->Draw("5 same");
    bin_0_jpsi->Draw("p same");
    bin_0_syst_psi2s->Draw("5 same");
    bin_0_psi2s->Draw("p same");
    bin_1_syst_jpsi->Draw("5 same");
    bin_1_jpsi->Draw("p same");
    bin_1_syst_psi2s->Draw("5 same");
    bin_1_psi2s->Draw("p same");
    canv.Update();
  }

  // Set TLegend
  //drawopt==0
  TLegend leg(0.57,0.50,0.80,0.64,NULL,"brNDC");
  SetLegendStyle(leg);
  leg.SetMargin(0.25);
  leg.SetTextSize(0.05);
  //drawopt==1
  TLegend leg1(0.20,0.67,0.40,0.89,NULL,"brNDC");
  SetLegendStyle(leg1);
  leg1.SetNColumns(2);
  leg1.SetMargin(0.5);
  leg1.SetTextSize(0.032);
  
  TLegendEntry *entry;

  // TEXTS
  TLatex tl;
  tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.044);

  // Draw upper limits 
  TArrow arrow; arrow.SetNDC(); arrow.SetArrowSize(0.03); arrow.SetLineWidth(2);
  arrow.SetOption("<-|");
  
  if (drawopt==0) { // 2 graphs in each canvas
    canv.cd();
    entry = leg.AddEntry("bin_0_jpsi", "Prompt J/#psi","p");
    entry = leg.AddEntry("bin_0_psi2s","Prompt #psi(2S)","p");
    if (compare==0) { // vs. cent
      canv.cd();
      tl.SetTextColor(kBlack);
      tl.DrawLatex(0.21,0.82,"#splitline{|y| < 1.6}{6.5 < p_{T} < 30 GeV/c}");
      arrow.SetLineColor(SetColor(psi2scolor+drawopt));
      arrow.DrawArrow(psi2S_cent_mid_fx, 0.01, psi2S_cent_mid_fx, psi2S_cent_mid_fy); 
      canv1.cd();
      tl.SetTextColor(kBlack);
      tl.DrawLatex(0.21,0.82,"#splitline{1.6 < |y| < 2.4}{3 < p_{T} < 30 GeV/c}");
      arrow.DrawArrow(psi2S_cent_fwd_fx[0], 0.01, psi2S_cent_fwd_fx[0], psi2S_cent_fwd_fy[0]);
      arrow.DrawArrow(psi2S_cent_fwd_fx[1], 0.01, psi2S_cent_fwd_fx[1], psi2S_cent_fwd_fy[1]);
    } else if (compare==1) { // vs. pT
      canv.cd();
      tl.SetTextColor(kBlack);
      tl.DrawLatex(0.21,0.82,"#splitline{Cent. 0-100%}{|y| < 1.6}");
      canv1.cd();
      tl.SetTextColor(kBlack);
      tl.DrawLatex(0.21,0.82,"#splitline{Cent. 0-100%}{1.6 < |y| < 2.4}");
      arrow.SetLineColor(SetColor(psi2scolor+drawopt));
      arrow.DrawArrow(psi2S_pt_fwd_fx, 0.01, psi2S_pt_fwd_fx, psi2S_pt_fwd_fy);
    }
  }
  else if (drawopt==1) { // all 4 graphs in 1 canvas
    canv.cd();
    if (compare==0) { // vs. cent
      entry = leg1.AddEntry((TObject*)0, "|y| < 1.6, 6.5 < p_{T} < 30 GeV/c","");
      entry = leg1.AddEntry((TObject*)0,"","");
      entry = leg1.AddEntry("bin_0_jpsi", "Prompt J/#psi","p");
      entry = leg1.AddEntry("bin_0_psi2s","Prompt #psi(2S)","p");
      entry = leg1.AddEntry((TObject*)0,"1.6 < |y| < 2.4, 3 < p_{T} < 30 GeV/c","");
      entry = leg1.AddEntry((TObject*)0,"","");
      entry = leg1.AddEntry("bin_1_jpsi", "Prompt J/#psi","p");
      entry = leg1.AddEntry("bin_1_psi2s","Prompt #psi(2S)","p");
      arrow.SetLineColor(SetColor(psi2scolor));
      arrow.DrawArrow(psi2S_cent_mid_fx, 0.01, psi2S_cent_mid_fx, psi2S_cent_mid_fy); 
      arrow.SetLineColor(SetColor(psi2scolor+drawopt));
      arrow.DrawArrow(psi2S_cent_fwd_fx[0], 0.01, psi2S_cent_fwd_fx[0], psi2S_cent_fwd_fy[0]);
      arrow.DrawArrow(psi2S_cent_fwd_fx[1], 0.01, psi2S_cent_fwd_fx[1], psi2S_cent_fwd_fy[1]);
    } else if (compare==1) { // vs. pT
      entry = leg1.AddEntry((TObject*)0, "|y| < 1.6, Cent. 0-100%","");
      entry = leg1.AddEntry((TObject*)0,"","");
      entry = leg1.AddEntry("bin_0_jpsi", "Prompt J/#psi","p");
      entry = leg1.AddEntry("bin_0_psi2s","Prompt #psi(2S)","p");
      entry = leg1.AddEntry((TObject*)0,"1.6 < |y| < 2.4, Cent. 0-100%","");
      entry = leg1.AddEntry((TObject*)0,"","");
      entry = leg1.AddEntry("bin_1_jpsi", "Prompt J/#psi","p");
      entry = leg1.AddEntry("bin_1_psi2s","Prompt #psi(2S)","p");
      arrow.SetLineColor(SetColor(psi2scolor+drawopt));
      arrow.DrawArrow(psi2S_pt_fwd_fx, 0.01, psi2S_pt_fwd_fx, psi2S_pt_fwd_fy);
    }
  }
 
  canv.cd();
  if (compare == 0) // vs. cent
    CMS_lumi( (TPad*) gPad, 1061, 33, "", isPreliminary ); 
  else
    CMS_lumi( (TPad*) gPad, 106, 33, "", isPreliminary ); 
  canv1.cd();
  if (compare == 0) // vs. cent
    CMS_lumi( (TPad*) gPad, 1061, 33, "", isPreliminary ); 
  else
    CMS_lumi( (TPad*) gPad, 106, 33, "", isPreliminary ); 

  // draw a line at RAA=1
  double xpos1 = haxes->GetXaxis()->GetBinLowEdge(1);
  int lastbin = haxes->GetNbinsX();
  double xpos2 = haxes->GetXaxis()->GetBinLowEdge(lastbin) + haxes->GetXaxis()->GetBinWidth(lastbin);

  // Global uncertainty boxes
  TBox box;
  if (compare==0) { // vs. cent
    double x1=xpos2, x2=xpos2-20, x3=xpos2, x4=xpos2-20;
    canv.cd();
    box.SetFillColorAlpha(kRed-7,0.7);
    box.DrawBox(x2-20, Jpsi_cent_mid_global[0], x2, Jpsi_cent_mid_global[1]);
    box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
    box.DrawBox(x1-20, psi2S_cent_mid_global[0], x1, psi2S_cent_mid_global[1]);
    
    if (drawopt==0) { // 2 graphs in each canvas
      canv1.cd();
      box.SetFillColorAlpha(kRed-7,0.7);
    } else if (drawopt==1) { // all 4 graphs in 1 canvas
      x3 = x2-20;
      x4 = x3-20;
      box.SetFillColorAlpha(SetColor(jpsicolor+drawopt),0.7);
    }
    box.DrawBox(x4-20, Jpsi_cent_fwd_global[0], x4, Jpsi_cent_fwd_global[1]);
    box.SetFillColorAlpha(SetColor(psi2scolor+drawopt),0.7);
    box.DrawBox(x3-20, psi2S_cent_fwd_global[0], x3, psi2S_cent_fwd_global[1]);
  } else if (compare==1) { // vs. pT
    double x1=xpos2, x2=x1-1.25, x3=xpos2, x4=x3-1.25;
    canv.cd();
    box.SetFillColorAlpha(kRed-7,0.7);
    box.DrawBox(x2-1.25, Jpsi_pt_mid_global[0], x2, Jpsi_pt_mid_global[1]);
    box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
    box.DrawBox(x1-1.25, psi2S_pt_mid_global[0], x1, psi2S_pt_mid_global[1]);
    
    if (drawopt==0) { // 2 graphs in each canvas
      canv1.cd();
      box.SetFillColorAlpha(kRed-7,0.7);
    } else if (drawopt==1) { // all 4 graphs in 1 canvas
      x3 = x2-1.25;
      x4 = x3-1.25;
      box.SetFillColorAlpha(SetColor(jpsicolor+drawopt),0.7);
    }
    box.DrawBox(x4-1.25, Jpsi_pt_fwd_global[0], x4, Jpsi_pt_fwd_global[1]);
    box.SetFillColorAlpha(SetColor(psi2scolor+drawopt),0.7);
    box.DrawBox(x3-1.25, psi2S_pt_fwd_global[0], x3, psi2S_pt_fwd_global[1]);
  }

  TLine line; line.SetLineColor(kBlack); line.SetNDC();
  canv.cd();
  line.DrawLine(xpos1,1,xpos2,1);
  canv.Update();
  canv1.cd();
  line.DrawLine(xpos1,1,xpos2,1);
  canv1.Update();

  if (drawopt==0) {
    canv.cd();
    leg.Draw();
    canv.RedrawAxis();
    canv.SaveAs(Form("Psi2S_RAA_midrap_%s.pdf",compare==0?"cent":"pt"));

    canv1.cd();
    leg.Draw();
    canv1.RedrawAxis();
    canv1.SaveAs(Form("Psi2S_RAA_forward_%s.pdf",compare==0?"cent":"pt"));
  } else if (drawopt==1) {
    canv.cd();
    leg1.Draw();
    canv.RedrawAxis();
    canv.SaveAs(Form("Psi2S_RAA_allRegions_%s.pdf",compare==0?"cent":"pt"));
  }

}

void compareYnSJpsi(
    const bool doprompt,
    const int compare, //0: vs. cent, 1: vs. pT
    const int drawopt //0: Y(1,2,3S), J/psi, 1: Y(1,2,3S), J/psi, psi(2S) mid, 2: Y(1), J/psi, psi(2S)
    ) {

  setTDRStyle();

  // Set style indexes (avoid marker==4 cause it's for low pT J/psi RAA)
  const int jpsicolor = 0;
  const int jpsimarker = 0;
  const int psi2scolor = 6;
  const int psi2smarker = 1;
  const int y1scolor = 3;
  const int y1smarker = 2;

  // 16025 results
  string jpsi_cent_fname[2] = {
    "./root/result_JPsi_RAA_cent_prompt_0_accEffCorr.root",
    "./root/result_JPsi_RAA_cent_prompt_1_accEffCorr.root"
  };
  string jpsi_pt_fname[2] = {
    "./root/result_JPsi_RAA_pt_prompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_prompt_2_linearX_accEffCorr.root"
  };
  string npjpsi_cent_fname[2] = {
    "./root/result_JPsi_RAA_cent_nonprompt_0_accEffCorr.root",
    "./root/result_JPsi_RAA_cent_nonprompt_1_accEffCorr.root"
  };
  string npjpsi_pt_fname[2] = {
    "./root/result_JPsi_RAA_pt_nonprompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_nonprompt_2_linearX_accEffCorr.root"
  };
  // 16004 results
  string psi2s_fname[2] = {
    "./root16004/result_Psi2S_RAA_cent_prompt_16004_3_accEffCorr.root",
    "./root16004/result_Psi2S_RAA_pt_prompt_16004_3_linearX_accEffCorr.root"
  };

  TFile inputf_jpsi( compare==0 ? jpsi_cent_fname[1].c_str() : jpsi_pt_fname[1].c_str() );
  TFile inputf_jpsi_lowpt( compare==0 ? jpsi_cent_fname[0].c_str() : jpsi_pt_fname[0].c_str() );
  TFile inputf_npjpsi( compare==0 ? npjpsi_cent_fname[1].c_str() : npjpsi_pt_fname[1].c_str() );
  TFile inputf_npjpsi_lowpt( compare==0 ? npjpsi_cent_fname[0].c_str() : npjpsi_pt_fname[0].c_str() );
  TFile inputf_psi2s( psi2s_fname[compare].c_str() );

  // Load histograms
  // J/psi RAA: 0-2.4, 6.5-50, 0-100%
  TCanvas *c1_jpsi;
  if (doprompt) c1_jpsi = (TCanvas*)inputf_jpsi.Get("c1");
  else c1_jpsi = (TCanvas*)inputf_npjpsi.Get("c1");
  TGraphAsymmErrors *bin_0_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0_syst");
  bin_0_jpsi->SetName("bin_0_jpsi");
  bin_0_syst_jpsi->SetName("bin_0_syst_jpsi");
  SetGraphStyle(bin_0_jpsi, jpsicolor, 0, jpsimarker);
  SetGraphStyle(bin_0_syst_jpsi, jpsicolor, 1, jpsimarker);

  // J/psi RAA: 1.8-2.4, 3-6.5, 0-100%
  TCanvas *c1_jpsi_lowpt;
  if (doprompt) c1_jpsi_lowpt = (TCanvas*)inputf_jpsi_lowpt.Get("c1");
  else c1_jpsi_lowpt = (TCanvas*)inputf_npjpsi_lowpt.Get("c1");
  TGraphAsymmErrors *bin_lowpt_ = (TGraphAsymmErrors *)c1_jpsi_lowpt->FindObject("bin_3");
  TGraphAsymmErrors *bin_lowpt_syst_ = (TGraphAsymmErrors *)c1_jpsi_lowpt->FindObject("bin_3_syst");

  int npoint_lowpt = bin_lowpt_->GetN();
  Double_t *x_lowpt = bin_lowpt_->GetX();
  Double_t *y_lowpt = bin_lowpt_->GetY();
  Double_t *exlow_lowpt = bin_lowpt_->GetEXlow();
  Double_t *eylow_lowpt = bin_lowpt_->GetEYlow();
  Double_t *exhigh_lowpt = bin_lowpt_->GetEXhigh();
  Double_t *eyhigh_lowpt = bin_lowpt_->GetEYhigh();
  Double_t *eylow_syst_lowpt = bin_lowpt_syst_->GetEYlow();
  Double_t *eyhigh_syst_lowpt = bin_lowpt_syst_->GetEYhigh();

  Double_t g_x_lowpt[3], g_y_lowpt[3], g_exl_lowpt[3], g_exh_lowpt[3], g_eyl_lowpt[3], g_eyh_lowpt[3], g_eyl_syst_lowpt[3], g_eyh_syst_lowpt[3];
  for (int i=0; i<3; i++) {
    g_x_lowpt[i] = x_lowpt[i];
    g_y_lowpt[i] = y_lowpt[i];
    g_exl_lowpt[i] = exlow_lowpt[i];
    g_exh_lowpt[i] = exhigh_lowpt[i];
    g_eyl_lowpt[i] = eylow_lowpt[i];
    g_eyh_lowpt[i] = eyhigh_lowpt[i];
    g_eyl_syst_lowpt[i] = eylow_syst_lowpt[i];
    g_eyh_syst_lowpt[i] = eyhigh_syst_lowpt[i];
  }
  
  TGraphAsymmErrors *bin_lowpt = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_lowpt,g_eyh_lowpt);
  TGraphAsymmErrors *bin_lowpt_syst = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_syst_lowpt,g_eyh_syst_lowpt);
  bin_lowpt->SetName("bin_lowpt");
  bin_lowpt_syst->SetName("bin_lowpt_syst");
  SetGraphStyle(bin_lowpt, jpsicolor, 0, 4);
  SetGraphStyle(bin_lowpt_syst, jpsicolor, 1, 4);
 
  // psi(2S) RAA: 0-1.6, 6.5-30, 0-100%
  TCanvas *c1_psi2s = (TCanvas*)inputf_psi2s.Get("c1");
  TGraphAsymmErrors *bin_0_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_0_syst");
  bin_0_psi2s->SetName("bin_0_psi2s");
  bin_0_syst_psi2s->SetName("bin_0_syst_psi2s");
  SetGraphStyle(bin_0_psi2s, psi2scolor, 0, psi2smarker);
  SetGraphStyle(bin_0_syst_psi2s, psi2scolor, 1, psi2smarker);
  
  // psi(2S) RAA: 1.6-2.4, 3-50, 0-100%
  TGraphAsymmErrors *bin_1_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1");
  TGraphAsymmErrors *bin_1_syst_psi2s = (TGraphAsymmErrors *)c1_psi2s->FindObject("bin_1_syst");
  bin_1_psi2s->SetName("bin_1_psi2s");
  bin_1_syst_psi2s->SetName("bin_1_syst_psi2s");
  SetGraphStyle(bin_1_psi2s, psi2scolor+1, 0, psi2smarker);
  SetGraphStyle(bin_1_syst_psi2s, psi2scolor+1, 1, psi2smarker);

  // Y(1,2,3S) graphs
  TGraphAsymmErrors *Y1SRAA, *Y1SRAA_syst;
  TGraphAsymmErrors *Y2SRAA, *Y2SRAA_syst;
  TGraphAsymmErrors *Y3SRAA;
  if (compare==0) {
    Y1SRAA = new TGraphAsymmErrors(9,Y1SRAA_cent_fx,Y1SRAA_cent_fy,Y1SRAA_cent_fex,Y1SRAA_cent_fex,Y1SRAA_cent_fey_stat,Y1SRAA_cent_fey_stat);
    Y1SRAA_syst = new TGraphAsymmErrors(9,Y1SRAA_cent_fx,Y1SRAA_cent_fy,Y1SRAA_cent_fex,Y1SRAA_cent_fex,Y1SRAA_cent_fely_syst,Y1SRAA_cent_fehy_syst);
    Y2SRAA = new TGraphAsymmErrors(9,Y2SRAA_cent_fx,Y2SRAA_cent_fy,Y2SRAA_cent_fex,Y2SRAA_cent_fex,Y2SRAA_cent_fey_stat,Y2SRAA_cent_fey_stat);
    Y2SRAA_syst = new TGraphAsymmErrors(9,Y2SRAA_cent_fx,Y2SRAA_cent_fy,Y2SRAA_cent_fex,Y2SRAA_cent_fex,Y2SRAA_cent_fely_syst,Y2SRAA_cent_fehy_syst);
  } else if (compare==1) {
    Y1SRAA = new TGraphAsymmErrors(6,Y1SRAA_pt_fx,Y1SRAA_pt_fy,Y1SRAA_pt_fex,Y1SRAA_pt_fex,Y1SRAA_pt_fey_stat,Y1SRAA_pt_fey_stat);
    Y1SRAA_syst = new TGraphAsymmErrors(6,Y1SRAA_pt_fx,Y1SRAA_pt_fy,Y1SRAA_pt_fex,Y1SRAA_pt_fex,Y1SRAA_pt_fey_syst,Y1SRAA_pt_fey_syst);
    Y2SRAA = new TGraphAsymmErrors(3,Y2SRAA_pt_fx,Y2SRAA_pt_fy,Y2SRAA_pt_fex,Y2SRAA_pt_fex,Y2SRAA_pt_fey_stat,Y2SRAA_pt_fey_stat);
    Y2SRAA_syst = new TGraphAsymmErrors(3,Y2SRAA_pt_fx,Y2SRAA_pt_fy,Y2SRAA_pt_fex,Y2SRAA_pt_fex,Y2SRAA_pt_fey_syst,Y2SRAA_pt_fey_syst);
  }
  Y1SRAA->SetName("Y1SRAA");
  Y1SRAA_syst->SetName("Y1SRAA_syst");
  SetGraphStyle(Y1SRAA, y1scolor, 0, y1smarker);
  SetGraphStyle(Y1SRAA_syst, y1scolor, 1, y1smarker);
  Y2SRAA->SetName("Y2SRAA");
  Y2SRAA_syst->SetName("Y2SRAA_syst");
  SetGraphStyle(Y2SRAA, y1scolor+1, 0, y1smarker);
  SetGraphStyle(Y2SRAA_syst, y1scolor+1, 1, y1smarker);
  Y3SRAA = new TGraphAsymmErrors(); //for drawing, set some color, marker
  Y3SRAA->SetName("Y3SRAA");
  SetGraphStyle(Y3SRAA, y1scolor+2, 0, y1smarker);
  Y3SRAA->SetLineWidth(2);

  // Dummy objects for drawing
  TH1F *haxes = (TH1F*)c1_jpsi->FindObject("haxes");
  
  // Draw graphs
  gStyle->SetEndErrorSize(5);
  TCanvas canv("canv","canv",600,600);
  haxes->GetYaxis()->SetRangeUser(0,1.7);
  haxes->Draw("");
  
  //0: Y(1,2,3S), J/psi, 1: Y(1,2S), J/psi, psi(2S), 2: Y(1), J/psi, psi(2S)
  if (drawopt==0) {
    Y1SRAA_syst->Draw("5 same");
    Y1SRAA->Draw("p same");
    Y2SRAA_syst->Draw("5 same");
    Y2SRAA->Draw("p same");
    Y3SRAA->Draw("p same");
    bin_0_syst_jpsi->Draw("5 same");
    bin_0_jpsi->Draw("p same");
    if (compare==1) { // vs. pt
      bin_lowpt_syst->Draw("5 same");
      bin_lowpt->Draw("p same");
    }
  } else if (drawopt==1) {
    Y1SRAA_syst->Draw("5 same");
    Y1SRAA->Draw("p same");
    Y2SRAA_syst->Draw("5 same");
    Y2SRAA->Draw("p same");
    Y3SRAA->Draw("p same");
    bin_0_syst_jpsi->Draw("5 same");
    bin_0_jpsi->Draw("p same");
    if (compare==1) { // vs. pt
      bin_lowpt_syst->Draw("5 same");
      bin_lowpt->Draw("p same");
    }
    if (doprompt) {
      bin_0_syst_psi2s->Draw("5 same");
      bin_0_psi2s->Draw("p same");
    }
  } else if (drawopt==2) {
    Y1SRAA_syst->Draw("5 same");
    Y1SRAA->Draw("p same");
    bin_0_syst_jpsi->Draw("5 same");
    bin_0_jpsi->Draw("p same");
    if (compare==1) { // vs. pt
      bin_lowpt_syst->Draw("5 same");
      bin_lowpt->Draw("p same");
    }
    if (doprompt) {
      bin_0_syst_psi2s->Draw("5 same");
      bin_0_psi2s->Draw("p same");
      bin_1_syst_psi2s->Draw("5 same");
      bin_1_psi2s->Draw("p same");
    }
  }
  canv.Update();

  // Set TLegend
  TLegend leg1(0.23,0.715,0.63,0.85,NULL,"brNDC");
  SetLegendStyle(leg1);
  leg1.SetMargin(0.07);
  leg1.SetTextSize(0.033);
  if (compare==1) {
    leg1.SetX2(0.55);
    leg1.SetMargin(0.2);
    leg1.SetNColumns(2); // vs. pT
  }

  TLegend leg(0.23,0.615,0.63,0.705,NULL,"brNDC");
  SetLegendStyle(leg);
  leg.SetMargin(0.2);
  leg.SetTextSize(0.033);
  if (drawopt==0 || drawopt==1) leg.SetNColumns(3);
  
  TLegend leg2(0.23,0.715,0.63,0.85,NULL,"brNDC");
  SetLegendStyle(leg2);
  leg2.SetMargin(0.07);
  leg2.SetTextSize(0.033);

  if (compare==0) { // vs. cent
    if (doprompt) {
      leg1.SetY1(0.72);
      leg1.SetY2(0.89);
      leg.SetY1(0.63);
      leg.SetY2(0.70);
      leg2.SetY1(0.60);
      leg2.SetY2(0.625);
    } else {
      leg1.SetY1(0.77);
      leg1.SetY2(0.89);
      leg.SetY1(0.66);
      leg.SetY2(0.74);
      leg2.SetY1(0.63);
      leg2.SetY2(0.655);
    }
  } else if (compare==1) { // vs. pT
    if (doprompt) {
      leg1.SetY1(0.72);
      leg1.SetY2(0.89);
      leg.SetY1(0.63);
      leg.SetY2(0.69);
      leg2.SetY1(0.60);
      leg2.SetY2(0.625);
    } else {
      leg1.SetY1(0.77);
      leg1.SetY2(0.88);
      leg.SetY1(0.65);
      leg.SetY2(0.72);
      leg2.SetY1(0.62);
      leg2.SetY2(0.645);
    }
  }
 
  TLegendEntry *entry;

  // TEXTS
  TLatex tl;
  tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.037);
  tl.SetTextColor(kBlack);

  // Draw upper limits 
  TArrow arrow; arrow.SetNDC(); arrow.SetArrowSize(0.03); arrow.SetLineWidth(2);
  arrow.SetOption("<-|");
  
  // draw a line at RAA=1
  double xpos1 = haxes->GetXaxis()->GetBinLowEdge(1);
  int lastbin = haxes->GetNbinsX();
  double xpos2 = haxes->GetXaxis()->GetBinLowEdge(lastbin) + haxes->GetXaxis()->GetBinWidth(lastbin);
  if (compare==0) xpos2=420;

  // Global uncertainty boxes
  TBox box;

  if (compare == 0) // vs. cent
    CMS_lumi( (TPad*) gPad, 1061, 33, "", isPreliminary ); 
  else
    CMS_lumi( (TPad*) gPad, 106, 33, "", isPreliminary ); 
  
  if (compare==0) { // vs. cent
    if (drawopt==0) {
      if (doprompt) leg1.SetHeader("#scale[1.15]{Charmonia}");
      else leg1.SetHeader("#scale[1.15]{Open beauty}");
      if (doprompt) {
        //entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi} #scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","p");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","");
      } else {
        //entry = leg1.AddEntry("bin_0_jpsi", "#scale[0.83]{Nonprompt J/#psi} #scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","p");
        entry = leg1.AddEntry("bin_0_jpsi", "#scale[0.83]{J/#psi from b hadrons}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","");
      }
      //if (doprompt) leg.SetHeader("#scale[1.15]{Bottomonia} #scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}");
      //else leg.SetHeader("#scale[1.15]{Hidden beauty} #scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}");
      if (doprompt) leg.SetHeader("#scale[1.15]{Bottomonia}"); 
      else leg.SetHeader("#scale[1.15]{Hidden beauty}");
      entry = leg.AddEntry("Y1SRAA", "#scale[0.83]{#varUpsilon(1S)}","p");
      entry = leg.AddEntry("Y2SRAA", "#scale[0.83]{#varUpsilon(2S)}","p");
      entry = leg.AddEntry("Y3SRAA", "#scale[0.83]{#varUpsilon(3S)}","l");
      entry = leg2.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}","");
    } else if (drawopt==1) {
      if (doprompt) leg1.SetHeader("#scale[1.15]{Charmonia}");
      else leg1.SetHeader("#scale[1.15]{Open beauty}");
      if (doprompt) {
        //entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi} #scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","p");
        //entry = leg1.AddEntry("bin_0_psi2s","#scale[0.83]{Prompt #psi(2S)} #scale[0.75]{|y| < 1.6, 6.5 < p_{T} < 30 GeV/c}","p");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","");
        entry = leg1.AddEntry("bin_0_psi2s","#scale[0.83]{Prompt #psi(2S)}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 1.6, 6.5 < p_{T} < 30 GeV/c}","");
      } else {
        //entry = leg1.AddEntry("bin_0_jpsi", "#scale[0.83]{Nonprompt J/#psi} #scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","p");
        entry = leg1.AddEntry("bin_0_jpsi", "#scale[0.83]{J/#psi from b hadrons}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","");
      }
      if (doprompt) leg.SetHeader("#scale[1.15]{Bottomonia}"); 
      else leg.SetHeader("#scale[1.15]{Hidden beauty}");
      entry = leg.AddEntry("Y1SRAA", "#scale[0.83]{#varUpsilon(1S)}","p");
      entry = leg.AddEntry("Y2SRAA", "#scale[0.83]{#varUpsilon(2S)}","p");
      entry = leg.AddEntry("Y3SRAA", "#scale[0.83]{#varUpsilon(3S)}","l");
      entry = leg2.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}","");
    } else if (drawopt==2) {
      tl.SetTextSize(0.050);
      if (doprompt) tl.DrawLatex(0.56,0.57,"Prompt J/#psi, #psi(2S)");
      else tl.DrawLatex(0.62,0.57,"J/#psi from b hadrons");
      entry = leg.AddEntry("bin_0_jpsi", "J/#psi #scale[0.8]{|y| < 2.4, 6.5 < p_{T} < 50 GeV/c}","p");
      if (doprompt) {
        entry = leg.AddEntry("bin_0_psi2s","#psi(2S) #scale[0.8]{|y| < 1.6, 6.5 < p_{T} < 30 GeV/c}","p");
        entry = leg.AddEntry("bin_1_psi2s","#psi(2S) #scale[0.8]{(1.6 < |y| < 2.4, 3 < p_{T} < 50 GeV/c)}","p");
      }
      entry = leg.AddEntry("Y1SRAA", "#varUpsilon(1S) #scale[0.8]{|y| < 2.4, p_{T} < 30 GeV/c}","p");
    }
  } else if (compare==1) { // vs. pT
    tl.SetTextSize(0.037);
    tl.DrawLatex(0.64,0.515,"Cent. 0-100%");
    if (drawopt==0) {
      if (doprompt) leg1.SetHeader("#scale[1.15]{Charmonia}");
      else leg1.SetHeader("#scale[1.15]{Open beauty}");
      if (doprompt) {
        //entry = leg1.AddEntry("bin_lowpt", "#scale[0.83]{Prompt J/#psi} #scale[0.75]{1.8 < |y| < 2.4}","p");
        //entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi} #scale[0.75]{|y| < 2.4}","p");
        entry = leg1.AddEntry("bin_lowpt", "#scale[0.83]{Prompt J/#psi}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{1.8 < |y| < 2.4}","");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{Prompt J/#psi}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4}","");
      } else {
        //entry = leg1.AddEntry("bin_lowpt", "#scale[0.83]{Nonprompt J/#psi} #scale[0.75]{1.8 < |y| < 2.4}","p");
        //entry = leg1.AddEntry("bin_0_jpsi", "#scale[0.83]{Nonprompt J/#psi} #scale[0.75]{|y| < 2.4}","p");
        entry = leg1.AddEntry("bin_lowpt", "#scale[0.83]{J/#psi from b hadrons}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{1.8 < |y| < 2.4}","");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.83]{J/#psi from b hadrons}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4}","");
      }
      if (doprompt) leg.SetHeader("#scale[1.15]{Bottomonia}"); 
      else leg.SetHeader("#scale[1.15]{Hidden beauty}");
      entry = leg.AddEntry("Y1SRAA", "#scale[0.83]{#varUpsilon(1S)}","p");
      entry = leg.AddEntry("Y2SRAA", "#scale[0.83]{#varUpsilon(2S)}","p");
      entry = leg.AddEntry("Y3SRAA", "#scale[0.83]{#varUpsilon(3S)}","l");
      entry = leg2.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}","");
    } else if (drawopt==1) {
      if (doprompt) leg1.SetHeader("#scale[1.15]{Charmonia}");
      else leg1.SetHeader("#scale[1.15]{Open beauty}");
      if (doprompt) {
        entry = leg1.AddEntry((TObject*)0, "#scale[0.83]{Prompt J/#psi}","");
        entry = leg1.AddEntry((TObject*)0,"","");
        entry = leg1.AddEntry("bin_lowpt","#scale[0.75]{1.8 < |y| < 2.4}","p");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.75]{|y| < 2.4}","p");
        entry = leg1.AddEntry((TObject*)0,"#scale[0.83]{Prompt #psi(2S)}","");
        entry = leg1.AddEntry((TObject*)0,"","");
        entry = leg1.AddEntry("bin_0_psi2s","#scale[0.75]{|y| < 1.6}","p");
      } else {
        entry = leg1.AddEntry((TObject*)0, "#scale[0.83]{J/#psi from b hadrons}","");
        entry = leg1.AddEntry((TObject*)0,"","");
        entry = leg1.AddEntry("bin_lowpt","#scale[0.75]{1.8 < |y| < 2.4}","p");
        entry = leg1.AddEntry("bin_0_jpsi","#scale[0.75]{|y| < 2.4}","p");
      }
      if (doprompt) leg.SetHeader("#scale[1.15]{Bottomonia}"); 
      else leg.SetHeader("#scale[1.15]{Hidden beauty}");
      entry = leg.AddEntry("Y1SRAA", "#scale[0.83]{#varUpsilon(1S)}","p");
      entry = leg.AddEntry("Y2SRAA", "#scale[0.83]{#varUpsilon(2S)}","p");
      entry = leg.AddEntry("Y3SRAA", "#scale[0.83]{#varUpsilon(3S)}","l");
      entry = leg2.AddEntry((TObject*)0,"#scale[0.75]{|y| < 2.4, p_{T} < 30 GeV/c}, #scale[0.67]{HIN-16-023}","");
    } else if (drawopt==2) {
      tl.SetTextSize(0.050);
      if (doprompt) tl.DrawLatex(0.56,0.57,"Prompt J/#psi, #psi(2S)");
      else tl.DrawLatex(0.62,0.57,"J/#psi from b hadrons");
      entry = leg.AddEntry("bin_lowpt", "J/#psi #scale[0.8]{1.8 < |y| < 2.4}","p");
      entry = leg.AddEntry("bin_0_jpsi", "J/#psi #scale[0.8]{|y| < 2.4}","p");
      if (doprompt) {
        entry = leg.AddEntry("bin_0_psi2s","#psi(2S) #scale[0.8]{|y| < 1.6}","p");
        entry = leg.AddEntry("bin_1_psi2s","#psi(2S) #scale[0.8]{1.6 < |y| < 2.4}","p");
      }
      entry = leg.AddEntry("Y1SRAA", "#varUpsilon(1S) #scale[0.8]{|y| < 2.4}","p");
    }
  }

  // Draw C.L.
  if (compare==0) { // vs. cent
    double dx=20;
    double x1=xpos2, x2=x1-dx;
    if (drawopt==0) {
      arrow.SetLineColor(SetColor(y1scolor+2));
      arrow.DrawArrow(Y3SRAA_cent_fx[0], 0.0, Y3SRAA_cent_fx[0], Y3SRAA_cent_fehy[0]); 
      arrow.DrawArrow(Y3SRAA_cent_fx[1], 0.0, Y3SRAA_cent_fx[1], Y3SRAA_cent_fehy[1]); 

      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      if (doprompt) box.DrawBox(x2, Jpsi_cent_rap0024_global[0], x1, Jpsi_cent_rap0024_global[1]);
      else box.DrawBox(x2, NPJpsi_cent_rap0024_global[0], x1, NPJpsi_cent_rap0024_global[1]);
    } else if (drawopt==1) {
      arrow.SetLineColor(SetColor(y1scolor+2));
      arrow.DrawArrow(Y3SRAA_cent_fx[0], 0.0, Y3SRAA_cent_fx[0], Y3SRAA_cent_fehy[0]); 
      arrow.DrawArrow(Y3SRAA_cent_fx[1], 0.0, Y3SRAA_cent_fx[1], Y3SRAA_cent_fehy[1]); 

      if (doprompt) {
        arrow.SetLineColor(SetColor(psi2scolor));
        arrow.DrawArrow(psi2S_cent_mid_fx, 0.01, psi2S_cent_mid_fx, psi2S_cent_mid_fy);
      }
      
      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      if (doprompt) box.DrawBox(x2, Jpsi_cent_rap0024_global[0], x1, Jpsi_cent_rap0024_global[1]);
      else box.DrawBox(x2, NPJpsi_cent_rap0024_global[0], x1, NPJpsi_cent_rap0024_global[1]);
      if (doprompt) {
        box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
        box.DrawBox(x2-dx, psi2S_cent_mid_global[0], x2, psi2S_cent_mid_global[1]);
      }
    } else if (drawopt==2) {
      if (doprompt) {
        arrow.SetLineColor(SetColor(psi2scolor));
        arrow.DrawArrow(psi2S_cent_mid_fx, 0.01, psi2S_cent_mid_fx, psi2S_cent_mid_fy); 
        arrow.SetLineColor(SetColor(psi2scolor+1));
        arrow.DrawArrow(psi2S_cent_fwd_fx[0], 0.01, psi2S_cent_fwd_fx[0], psi2S_cent_fwd_fy[0]);
        arrow.DrawArrow(psi2S_cent_fwd_fx[1], 0.01, psi2S_cent_fwd_fx[1], psi2S_cent_fwd_fy[1]);
      }
      
      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      if (doprompt) box.DrawBox(x2, Jpsi_cent_rap0024_global[0], x1, Jpsi_cent_rap0024_global[1]);
      else box.DrawBox(x2, NPJpsi_cent_rap0024_global[0], x1, NPJpsi_cent_rap0024_global[1]);
      if (doprompt) {
        box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
        box.DrawBox(x2-dx, psi2S_cent_mid_global[0], x2, psi2S_cent_mid_global[1]);
        box.SetFillColorAlpha(SetColor(psi2scolor+1),0.7);
        box.DrawBox(x2-dx*2, psi2S_cent_fwd_global[0], x2-dx, psi2S_cent_fwd_global[1]);
      }
    }
  } else if (compare==1) { // vs. pT
    double dx=2.5;
    double x1=xpos2, x2=x1-dx;
    if (drawopt==0) {
      arrow.SetLineColor(SetColor(y1scolor+2));
      arrow.DrawArrow(Y3SRAA_pt_fx[0], 0.0, Y3SRAA_pt_fx[0], Y3SRAA_pt_fehy[0]); 
      arrow.DrawArrow(Y3SRAA_pt_fx[1], 0.0, Y3SRAA_pt_fx[1], Y3SRAA_pt_fehy[1]); 

      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      box.DrawBox(x2, Jpsi_pt_rap0024_global[0], x1, Jpsi_pt_rap0024_global[1]);
    } else if (drawopt==1) {
      arrow.SetLineColor(SetColor(y1scolor+2));
      arrow.DrawArrow(Y3SRAA_pt_fx[0], 0.0, Y3SRAA_pt_fx[0], Y3SRAA_pt_fehy[0]); 
      arrow.DrawArrow(Y3SRAA_pt_fx[1], 0.0, Y3SRAA_pt_fx[1], Y3SRAA_pt_fehy[1]); 

      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      box.DrawBox(x2, Jpsi_pt_rap0024_global[0], x1, Jpsi_pt_rap0024_global[1]);
      if (doprompt) {
        box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
        box.DrawBox(x2-dx, psi2S_pt_mid_global[0], x2, psi2S_pt_mid_global[1]);
      }
    } else if (drawopt==2) {
      if (doprompt) {
        arrow.SetLineColor(SetColor(psi2scolor+1));
        arrow.DrawArrow(psi2S_pt_fwd_fx, 0.01, psi2S_pt_fwd_fx, psi2S_pt_fwd_fy);
      }
      
      box.SetFillColorAlpha(SetColor(jpsicolor),0.7);
      box.DrawBox(x2, Jpsi_pt_rap0024_global[0], x1, Jpsi_pt_rap0024_global[1]);
      if (doprompt) {
        box.SetFillColorAlpha(SetColor(psi2scolor),0.7);
        box.DrawBox(x2-dx*2, psi2S_pt_mid_global[0], x2-dx, psi2S_pt_mid_global[1]);
        box.SetFillColorAlpha(SetColor(psi2scolor+1),0.7);
        box.DrawBox(x2-dx*3, psi2S_pt_fwd_global[0], x2-dx*2, psi2S_pt_fwd_global[1]);
      }
    }
  }

  TLine line; line.SetLineColor(kBlack); line.SetNDC();
  line.DrawLine(xpos1,1,xpos2,1);
  canv.Update();

  leg.Draw();
  leg1.Draw();
  leg2.Draw();
  canv.RedrawAxis();
  if (doprompt) canv.SaveAs(Form("Charmonia_Bottomonia_RAA_opt%d_%s.pdf",drawopt,compare==0?"cent":"pt"));
  else canv.SaveAs(Form("NPJpsi_Bottomonia_RAA_opt%d_%s.pdf",drawopt,compare==0?"cent":"pt"));

}


void extendLowPtFor14005(
    const bool doprompt
    ) {

  setTDRStyle();

  // Set style indexes (avoid marker==4 cause it's for low pT J/psi RAA)
  const int hin14005color = -2;
  const int hin14005marker = 0;
  const int jpsicolor = 0;
  const int jpsimarker = 0;

  // 16025 results
  string jpsi_pt_fname[2] = {
    "./root/result_JPsi_RAA_pt_prompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_prompt_14005_2_linearX_accEffCorr.root"
  };
  string npjpsi_pt_fname[2] = {
    "./root/result_JPsi_RAA_pt_nonprompt_0_linearX_accEffCorr.root",
    "./root/result_JPsi_RAA_pt_nonprompt_14005_2_linearX_accEffCorr.root"
  };
  
  TFile inputf_jpsi(jpsi_pt_fname[1].c_str() );
  TFile inputf_jpsi_lowpt(jpsi_pt_fname[0].c_str() );
  TFile inputf_npjpsi(npjpsi_pt_fname[1].c_str() );
  TFile inputf_npjpsi_lowpt(npjpsi_pt_fname[0].c_str() );

  // Load histograms
  // J/psi RAA: 0-2.4, 6.5-50, 0-100%
  TCanvas *c1_jpsi;
  if (doprompt) c1_jpsi = (TCanvas*)inputf_jpsi.Get("c1");
  else c1_jpsi = (TCanvas*)inputf_npjpsi.Get("c1");
  TGraphAsymmErrors *bin_0_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0");
  TGraphAsymmErrors *bin_0_syst_jpsi = (TGraphAsymmErrors *)c1_jpsi->FindObject("bin_0_syst");
  bin_0_jpsi->SetName("bin_0_jpsi");
  bin_0_syst_jpsi->SetName("bin_0_syst_jpsi");
  SetGraphStyle(bin_0_jpsi, jpsicolor, 0, jpsimarker);
  SetGraphStyle(bin_0_syst_jpsi, jpsicolor, 1, jpsimarker);

  // J/psi RAA: 1.8-2.4, 3-6.5, 0-100%
  TCanvas *c1_jpsi_lowpt;
  if (doprompt) c1_jpsi_lowpt = (TCanvas*)inputf_jpsi_lowpt.Get("c1");
  else c1_jpsi_lowpt = (TCanvas*)inputf_npjpsi_lowpt.Get("c1");
  TGraphAsymmErrors *bin_lowpt_ = (TGraphAsymmErrors *)c1_jpsi_lowpt->FindObject("bin_3");
  TGraphAsymmErrors *bin_lowpt_syst_ = (TGraphAsymmErrors *)c1_jpsi_lowpt->FindObject("bin_3_syst");

  int npoint_lowpt = bin_lowpt_->GetN();
  Double_t *x_lowpt = bin_lowpt_->GetX();
  Double_t *y_lowpt = bin_lowpt_->GetY();
  Double_t *exlow_lowpt = bin_lowpt_->GetEXlow();
  Double_t *eylow_lowpt = bin_lowpt_->GetEYlow();
  Double_t *exhigh_lowpt = bin_lowpt_->GetEXhigh();
  Double_t *eyhigh_lowpt = bin_lowpt_->GetEYhigh();
  Double_t *eylow_syst_lowpt = bin_lowpt_syst_->GetEYlow();
  Double_t *eyhigh_syst_lowpt = bin_lowpt_syst_->GetEYhigh();

  Double_t g_x_lowpt[3], g_y_lowpt[3], g_exl_lowpt[3], g_exh_lowpt[3], g_eyl_lowpt[3], g_eyh_lowpt[3], g_eyl_syst_lowpt[3], g_eyh_syst_lowpt[3];
  for (int i=0; i<3; i++) {
    g_x_lowpt[i] = x_lowpt[i];
    g_y_lowpt[i] = y_lowpt[i];
    g_exl_lowpt[i] = exlow_lowpt[i];
    g_exh_lowpt[i] = exhigh_lowpt[i];
    g_eyl_lowpt[i] = eylow_lowpt[i];
    g_eyh_lowpt[i] = eyhigh_lowpt[i];
    g_eyl_syst_lowpt[i] = eylow_syst_lowpt[i];
    g_eyh_syst_lowpt[i] = eyhigh_syst_lowpt[i];
  }
  
  TGraphAsymmErrors *bin_lowpt = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_lowpt,g_eyh_lowpt);
  TGraphAsymmErrors *bin_lowpt_syst = new TGraphAsymmErrors(3,g_x_lowpt,g_y_lowpt,g_exl_lowpt,g_exh_lowpt,g_eyl_syst_lowpt,g_eyh_syst_lowpt);
  bin_lowpt->SetName("bin_lowpt");
  bin_lowpt_syst->SetName("bin_lowpt_syst");
  SetGraphStyle(bin_lowpt, jpsicolor, 0, 4);
  SetGraphStyle(bin_lowpt_syst, jpsicolor, 1, 4);

  // 14005 J/psi RAA: 0-2.4, 6.5-30, 0-100%
  TGraphErrors *jpsi_14005;
  TGraphErrors *jpsi_14005_syst;
  if (doprompt) {
    jpsi_14005 = new TGraphErrors(6,prJpsi_14005_pt_rap0024_x,prJpsi_14005_pt_rap0024_y,prJpsi_14005_pt_rap0024_ex_stat,prJpsi_14005_pt_rap0024_ey_stat);
    jpsi_14005_syst = new TGraphErrors(6,prJpsi_14005_pt_rap0024_x,prJpsi_14005_pt_rap0024_y,prJpsi_14005_pt_rap0024_ex_syst,prJpsi_14005_pt_rap0024_ey_syst);
  } else {
    jpsi_14005 = new TGraphErrors(6,npJpsi_14005_pt_rap0024_x,npJpsi_14005_pt_rap0024_y,npJpsi_14005_pt_rap0024_ex_stat,npJpsi_14005_pt_rap0024_ey_stat);
    jpsi_14005_syst = new TGraphErrors(6,npJpsi_14005_pt_rap0024_x,npJpsi_14005_pt_rap0024_y,npJpsi_14005_pt_rap0024_ex_syst,npJpsi_14005_pt_rap0024_ey_syst);
  }
  jpsi_14005->SetName("jpsi_14005");
  jpsi_14005_syst->SetName("jpsi_14005_syst");
  SetGraphStyle(jpsi_14005, hin14005color, 0, hin14005marker);
  SetGraphStyle(jpsi_14005_syst, hin14005color, 1, hin14005marker);

  // J/psi RAA: 1.6-2.4, 3-6.5, 0-100%
  TGraphErrors *jpsi_14005_lowpt;
  TGraphErrors *jpsi_14005_syst_lowpt;
  if (doprompt) {
    jpsi_14005_lowpt = new TGraphErrors(3,prJpsi_14005_pt_rap1624_x,prJpsi_14005_pt_rap1624_y,prJpsi_14005_pt_rap1624_ex_stat,prJpsi_14005_pt_rap1624_ey_stat);
    jpsi_14005_syst_lowpt = new TGraphErrors(3,prJpsi_14005_pt_rap1624_x,prJpsi_14005_pt_rap1624_y,prJpsi_14005_pt_rap1624_ex_syst,prJpsi_14005_pt_rap1624_ey_syst);
  } else {
    jpsi_14005_lowpt = new TGraphErrors(3,npJpsi_14005_pt_rap1624_x,npJpsi_14005_pt_rap1624_y,npJpsi_14005_pt_rap1624_ex_stat,npJpsi_14005_pt_rap1624_ey_stat);
    jpsi_14005_syst_lowpt = new TGraphErrors(3,npJpsi_14005_pt_rap1624_x,npJpsi_14005_pt_rap1624_y,npJpsi_14005_pt_rap1624_ex_syst,npJpsi_14005_pt_rap1624_ey_syst);
  }
  jpsi_14005_lowpt->SetName("jpsi_14005_lowpt");
  jpsi_14005_syst_lowpt->SetName("jpsi_14005_syst_lowpt");
  SetGraphStyle(jpsi_14005_lowpt, hin14005color, 0, 4);
  SetGraphStyle(jpsi_14005_syst_lowpt, hin14005color, 1, 4);
 

  // Dummy objects for drawing
  TH1F *haxes = (TH1F*)c1_jpsi->FindObject("haxes");
  
  // Draw graphs
  gStyle->SetEndErrorSize(5);
  TCanvas canv("canv","canv",600,600);
  haxes->GetYaxis()->SetRangeUser(0,1.5);
  haxes->Draw("");
  
  jpsi_14005_syst->Draw("5 same");
  jpsi_14005->Draw("p same"); 
  jpsi_14005_syst_lowpt->Draw("5 same");
  jpsi_14005_lowpt->Draw("p same"); 
  bin_0_syst_jpsi->Draw("5 same");
  bin_0_jpsi->Draw("p same");
  bin_lowpt_syst->Draw("5 same");
  bin_lowpt->Draw("p same");
   
  canv.Update();

  // Set TLegend
  //drawopt==0
  TLegend leg(0.27,0.51,0.67,0.71,NULL,"brNDC");
//  TLegend leg(0.25,0.56,0.67,0.74,NULL,"brNDC");
  SetLegendStyle(leg);
  leg.SetTextSize(0.036);
  leg.SetNColumns(2);
  leg.SetMargin(0.2);
  
  TLegendEntry *entry;

  // TEXTS
  TLatex tl;
  tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.037);
  tl.SetTextColor(kBlack);

  // draw a line at RAA=1
  double xpos1 = haxes->GetXaxis()->GetBinLowEdge(1);
  int lastbin = haxes->GetNbinsX();
  double xpos2 = haxes->GetXaxis()->GetBinLowEdge(lastbin) + haxes->GetXaxis()->GetBinWidth(lastbin);

  // x-axis is always vs. pT -> use 106 option
  CMS_lumi( (TPad*) gPad, 106, 33, "", isPreliminary ); 
  
  tl.SetTextSize(0.037);
  tl.DrawLatex(0.20,0.795,"Cent. 0-100%");
  
  tl.SetTextSize(0.057); 
  if (doprompt) tl.DrawLatex(0.20,0.85,"Prompt J/#psi");
  else tl.DrawLatex(0.20,0.85,"J/#psi from b hadrons");

  entry = leg.AddEntry((TObject*)0, "#splitline{#sqrt[]{s_{NN}} = 2.76 TeV}{#scale[0.6]{(Eur.Phys.J. C77 (2017) no.4, 252)}}","");
  entry = leg.AddEntry((TObject*)0,"","");
  entry = leg.AddEntry("jpsi_14005_lowpt", "1.6 < |y| < 2.4","p");
  entry = leg.AddEntry("jpsi_14005", "|y| < 2.4","p");
  entry = leg.AddEntry((TObject*)0, "#sqrt[]{s_{NN}} = 5.02 TeV","");
  entry = leg.AddEntry((TObject*)0,"","");
  entry = leg.AddEntry("bin_lowpt", "1.8 < |y| < 2.4","p");
  entry = leg.AddEntry("bin_0_jpsi", "|y| < 2.4","p");

//  entry = leg.AddEntry((TObject*)0, "#splitline{#sqrt[]{s_{NN}} = 2.76 TeV}{#scale[0.6]{(Eur.Phys.J. C77 (2017) no.4, 252)}}","");
//  entry = leg.AddEntry("jpsi_14005_lowpt", "1.6 < |y| < 2.4","p");
//  entry = leg.AddEntry("jpsi_14005", "|y| < 2.4","p");
//  entry = leg.AddEntry((TObject*)0, "#sqrt[]{s_{NN}} = 5.02 TeV","");
//  entry = leg.AddEntry("bin_lowpt", "1.8 < |y| < 2.4","p");
//  entry = leg.AddEntry("bin_0_jpsi", "|y| < 2.4","p");

  // Global uncertainty boxes
  TBox box;
  TBox *box1, *box2;

  double dx=1.25;
  double x1=xpos2, x2=x1-dx;

  box1 = new TBox(x2, Jpsi_14005_pt_rap0024_global[0], x1, Jpsi_14005_pt_rap0024_global[1]);
  box1->SetLineColor(kBlack);
  box1->Draw("l");
  box.SetLineColor(kBlack);
  box.SetFillColorAlpha(13, 1);
  box.DrawBox(x2, Jpsi_14005_pt_rap0024_global[0], x1, Jpsi_14005_pt_rap0024_global[1]);

  box2 = new TBox(x2-dx, Jpsi_pt_rap0024_global[0], x2, Jpsi_pt_rap0024_global[1]);
  box2->SetLineColor(kRed+2);
  box2->Draw("l");
  box.SetLineColor(kRed+2);
  box.SetFillColorAlpha(kRed-7,0.5);
  box.DrawBox(x2-dx, Jpsi_pt_rap0024_global[0], x2, Jpsi_pt_rap0024_global[1]);
  
  TLine line; line.SetLineColor(kBlack); line.SetNDC();
  line.DrawLine(xpos1,1,xpos2,1);
  canv.Update();

  leg.Draw();
  canv.RedrawAxis();
  canv.SaveAs(Form("result_JPsi_RAA_pt_%s_2RapRanges_1_linearX_accEffCorr.pdf", doprompt? "prompt" : "nonprompt"));

}


void SetLegendStyle(TLegend &leg) {
  leg.SetMargin(0.1);
  leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.03);
  leg.SetLineColor(1);
  leg.SetLineStyle(1);
  leg.SetLineWidth(1);
  leg.SetFillColor(10);
  leg.SetFillStyle(0);
}

int SetColor(const int idx, const bool bkg) {
  const int palette[] =    {kRed+2, kOrange+9, kOrange+1, kSpring+2, kTeal+2, kGreen+3, kAzure+7, kBlue+1, kViolet+1, kViolet+2};
  const int palettebkg[] = {kRed-9, kOrange-9, kOrange-4, kSpring-9, kTeal-8, kGreen-9, kAzure-9, kBlue-9, kMagenta-10, kViolet-9};
  //const int palette[] =    {kRed+2, kOrange+9, kOrange+1, kOrange+5, kGreen+3, kTeal+2, kAzure+7, kBlue+1, kViolet+1, kViolet+2};
  //const int palettebkg[] = {kRed-9, kOrange-9, kOrange-4, kYellow-9, kGreen-9, kTeal-8, kAzure-9, kBlue-9, kMagenta-10, kViolet-9};
  const int size = sizeof(palette)/sizeof(int);

  if (idx>=0) {
    if (bkg) return palettebkg[static_cast<int>(idx % size)];
    else return palette[static_cast<int>(idx % size)];
  } else {
    return 13-idx;
  }
}

int SetMarker(const int i) {
  if (i==0) return kFullSquare;
  else if (i==1) return kFullCircle;
  else if (i==2) return kFullStar;
  else if (i==3) return kFullCross;
  else if (i==4) return kFullDiamond;
  else if (i==5) return kOpenSquare;
  else if (i==6) return kOpenCircle;
  else if (i==7) return kOpenStar;
  else if (i==8) return kFullTriangleDown;
  else return kOpenCross;
}

void SetGraphStyle(TGraph *g, const int idx, const bool bkg, const int i) {
  if (bkg) g->SetFillColorAlpha(SetColor(idx,bkg), 0.3);

  g->SetLineColor(SetColor(idx));
  
  g->SetMarkerColor(SetColor(idx));
  g->SetMarkerStyle(SetMarker(i));
  if (SetMarker(i) == kFullStar) g->SetMarkerSize(2.3);
  else if (SetMarker(i) == kFullDiamond) g->SetMarkerSize(2.2);
  else if (SetMarker(i) == kFullCross) g->SetMarkerSize(2.0);
  else g->SetMarkerSize(1.5);

}
void SetGraphStyle(TGraph *g, const int idx, const bool bkg, const double alphaV, const int i) {
  if (bkg) g->SetFillColorAlpha(SetColor(idx,bkg), alphaV);

  g->SetLineColor(SetColor(idx));
  
  g->SetMarkerColor(SetColor(idx));
  g->SetMarkerStyle(SetMarker(i));
  if (SetMarker(i) == kFullStar) g->SetMarkerSize(2.3);
  else if (SetMarker(i) == kFullDiamond) g->SetMarkerSize(2.2);
  else if (SetMarker(i) == kFullCross) g->SetMarkerSize(2.0);
  else g->SetMarkerSize(1.5);

}
