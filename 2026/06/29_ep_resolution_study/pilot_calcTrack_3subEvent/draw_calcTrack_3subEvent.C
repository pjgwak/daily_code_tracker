#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::cout;
using std::endl;

void ApplyStyle()
{
  gROOT->SetBatch(kTRUE);
  const char *rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  if (!gSystem->AccessPathName(rootlogon))
    gROOT->ProcessLine(Form(".X %s", rootlogon));

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
}

bool HasBranch(TTree *t, const char *name)
{
  return t && t->GetBranch(name);
}

void SetupPad()
{
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);
  gPad->SetTicks(1, 1);
}

void SetupHist(TH1 *h, const char *xTitle, int color = kBlack)
{
  if (!h)
    return;
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitle("Entries");
  h->GetYaxis()->CenterTitle(true);
  const double maxY = h->GetMaximum();
  h->SetMaximum(maxY > 0 ? maxY * 1.6 : 1.0);
}

void DrawPanelLabel(const char *text, double x = 0.62, double y = 0.85)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.03);
  latex.DrawLatex(x, y, text);
}

void SaveCanvas(TCanvas *c, const TString &outDir, const TString &sample, const TString &name)
{
  c->SaveAs(Form("%s/plots_%s_%s.pdf", outDir.Data(), sample.Data(), name.Data()));
}

TH1 *DrawTreeHist(TTree *t, const char *expr, const char *histSpec, const char *cut,
                  const char *xTitle, const char *label, int color = kBlack,
                  double labelX = 0.62, double labelY = 0.85)
{
  if (!t)
    return nullptr;
  t->Draw(Form("%s>>%s", expr, histSpec), cut, "goff");
  TString histName(histSpec);
  histName.Resize(histName.First('('));
  TH1 *h = (TH1 *)gDirectory->Get(histName);
  if (!h)
  {
    cout << "[WARN] missing histogram: " << histName << endl;
    return nullptr;
  }
  SetupHist(h, xTitle, color);
  h->Draw("hist");
  DrawPanelLabel(label, labelX, labelY);
  return h;
}

void DrawTreeHist2D(TTree *t, const char *expr, const char *histSpec, const char *cut,
                    const char *xTitle, const char *yTitle, const char *label)
{
  if (!t)
    return;
  gPad->SetRightMargin(0.16);
  t->Draw(Form("%s>>%s", expr, histSpec), cut, "goff");
  TString histName(histSpec);
  histName.Resize(histName.First('('));
  TH2 *h = (TH2 *)gDirectory->Get(histName);
  if (!h)
  {
    cout << "[WARN] missing histogram: " << histName << endl;
    return;
  }
  h->SetStats(0);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetYaxis()->CenterTitle(true);
  h->Draw("colz");
  DrawPanelLabel(label, 0.18, 0.86);
}

void SetHistStyle(TH1 *h, int color, int marker, int width = 2)
{
  if (!h)
    return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(width);
  h->SetStats(0);
}


TH1 *MakeTreeHist(TTree *t, const char *expr, const char *histSpec, const char *cut,
                  const char *xTitle, int color)
{
  if (!t)
    return nullptr;
  t->Draw(Form("%s>>%s", expr, histSpec), cut, "goff");
  TString histName(histSpec);
  histName.Resize(histName.First('('));
  TH1 *h = (TH1 *)gDirectory->Get(histName);
  if (!h)
    return nullptr;
  SetupHist(h, xTitle, color);
  h->SetMaximum(h->GetMaximum() > 0 ? h->GetMaximum() * 1.35 : 1.0);
  return h;
}

void draw_calcTrack_3subEvent(
    const char *inputFile = "./roots/pr/skim_PbPb2023_isMC1_PR_EP_3Sub_job.root",
    const char *sampleArg = "pr")
{
  ApplyStyle();

  TString sample(sampleArg);
  TString outDir = Form("figs/%s", sample.Data());
  if (gSystem->AccessPathName(outDir))
    gSystem->mkdir(outDir, true);

  TFile *f = TFile::Open(inputFile, "READ");
  if (!f || f->IsZombie())
  {
    cout << "[ERROR] cannot open input file: " << inputFile << endl;
    return;
  }

  TTree *t = (TTree *)f->Get("hionia/myTree");
  if (!t)
  {
    cout << "[ERROR] hionia/myTree not found in " << inputFile << endl;
    return;
  }

  TH1 *hCos2TrkHFp = (TH1 *)f->Get("hCos2TrkHFp");
  TH1 *hCos2TrkHFm = (TH1 *)f->Get("hCos2TrkHFm");
  TH1 *hCos2HFpHFm = (TH1 *)f->Get("hCos2HFpHFm");
  TH1 *hResolution3Sub = (TH1 *)f->Get("hResolution3Sub");
  TH2 *hQ2xy = (TH2 *)f->Get("hQ2xy");

  const char *validTrk = "has_trk_ep && ep>-90 && ep_rec>-90 && ep_flat>-90";
  const char *validHF = "ep_hfp2>-90 && ep_hfm2>-90";
  const char *validAll = "has_trk_ep && ep_flat>-90 && ep_hfp2>-90 && ep_hfm2>-90";
  const char *deltaTrkHFp = "0.5*atan2(sin(2*(ep_flat-ep_hfp2)),cos(2*(ep_flat-ep_hfp2)))";
  const char *deltaTrkHFm = "0.5*atan2(sin(2*(ep_flat-ep_hfm2)),cos(2*(ep_flat-ep_hfm2)))";
  const char *deltaHFpHFm = "0.5*atan2(sin(2*(ep_hfp2-ep_hfm2)),cos(2*(ep_hfp2-ep_hfm2)))";
  const char *deltaCalcInputTrk = "0.5*atan2(sin(2*(ep_flat-rpAng[3])),cos(2*(ep_flat-rpAng[3])))";

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.03);

  // Same items as draw_ep_angles.py: raw, recentered, and flattened tracker EP.
  TCanvas *cEpGroup = new TCanvas("cEpGroup", "tracker EP angle group", 1500, 500);
  cEpGroup->Divide(3, 1);
  cEpGroup->cd(1); SetupPad(); DrawTreeHist(t, "ep", "h_ep(100,-1.5,1.5)", validTrk, "#Psi_{2,trk} (rad)", "Raw EP Angle", kBlack, 0.65, 0.85);
  cEpGroup->cd(2); SetupPad(); DrawTreeHist(t, "ep_rec", "h_ep_rec(100,-1.5,1.5)", validTrk, "#Psi_{2,trk} (rad)", "Recentering", kBlack, 0.68, 0.85);
  cEpGroup->cd(3); SetupPad(); DrawTreeHist(t, "ep_flat", "h_ep_flat(100,-1.5,1.5)", validTrk, "#Psi_{2,trk} (rad)", "Recentering + Flattening", kBlack, 0.50, 0.85);
  SaveCanvas(cEpGroup, outDir, sample, "ep_group");

  // Same items as draw_q_vectors.py plus the 2D Q-vector map from the macro output.
  TCanvas *cQGroup = new TCanvas("cQGroup", "Q-vector distributions", 1500, 500);
  cQGroup->Divide(3, 1);
  cQGroup->cd(1); SetupPad(); DrawTreeHist(t, "Q2x_tree", "h_qx(100,-0.5,0.5)", "has_trk_ep", "Q_{2,x}", "Q_{2,x}", kBlack, 0.72, 0.85);
  cQGroup->cd(2); SetupPad(); DrawTreeHist(t, "Q2y_tree", "h_qy(100,-0.5,0.5)", "has_trk_ep", "Q_{2,y}", "Q_{2,y}", kBlack, 0.72, 0.85);
  cQGroup->cd(3); SetupPad();
  if (hQ2xy)
  {
    gPad->SetRightMargin(0.16);
    hQ2xy->SetStats(0);
    hQ2xy->GetXaxis()->SetTitle("Q_{2,x}");
    hQ2xy->GetXaxis()->CenterTitle(true);
    hQ2xy->GetYaxis()->SetTitle("Q_{2,y}");
    hQ2xy->GetYaxis()->CenterTitle(true);
    hQ2xy->Draw("colz");
    DrawPanelLabel("Q_{2,x} vs Q_{2,y}", 0.18, 0.86);
  }
  SaveCanvas(cQGroup, outDir, sample, "q_vector");

  // Additional angle diagnostics not covered by the Python scripts.
  TCanvas *cHF = new TCanvas("cHF", "HF EP angle distributions", 1500, 500);
  cHF->Divide(3, 1);
  cHF->cd(1); SetupPad(); DrawTreeHist(t, "ep_hfp2", "h_ep_hfp2(100,-1.5,1.5)", validHF, "#Psi_{2,HF+} (rad)", "HF+ input", kRed + 1, 0.68, 0.85);
  cHF->cd(2); SetupPad(); DrawTreeHist(t, "ep_hfm2", "h_ep_hfm2(100,-1.5,1.5)", validHF, "#Psi_{2,HF-} (rad)", "HF- input", kBlue + 1, 0.68, 0.85);
  cHF->cd(3); SetupPad(); DrawTreeHist(t, "rpAng[2]", "h_ep_hf2(100,-1.5,1.5)", "nEP>2", "#Psi_{2,HF} (rad)", "HF combined input", kGreen + 2, 0.55, 0.85);
  SaveCanvas(cHF, outDir, sample, "hf_angles");

  TCanvas *cEpOverlay = new TCanvas("cEpOverlay", "EP angle overlay", 900, 700);
  SetupPad();
  TH1 *hEpAng = MakeTreeHist(t, "rpAng[2]", "h_overlay_epang(100,-1.5,1.5)", "nEP>2", "#Psi_{2} (rad)", kGreen + 2);
  TH1 *hHFp = MakeTreeHist(t, "ep_hfp2", "h_overlay_hfp(100,-1.5,1.5)", validHF, "#Psi_{2} (rad)", kRed + 1);
  TH1 *hHFm = MakeTreeHist(t, "ep_hfm2", "h_overlay_hfm(100,-1.5,1.5)", validHF, "#Psi_{2} (rad)", kBlue + 1);
  TH1 *hCalTrack = MakeTreeHist(t, "ep_flat", "h_overlay_caltrk(100,-1.5,1.5)", validTrk, "#Psi_{2} (rad)", kBlack);
  double overlayMax = 0.0;
  for (TH1 *h : {hEpAng, hHFp, hHFm, hCalTrack})
    if (h)
      overlayMax = std::max(overlayMax, h->GetMaximum());
  if (hEpAng)
  {
    hEpAng->SetTitle(";#Psi_{2} (rad);Entries");
    hEpAng->SetMaximum(overlayMax > 0 ? overlayMax * 1.2 : 1.0);
    hEpAng->Draw("hist");
  }
  if (hHFp)
    hHFp->Draw(hEpAng ? "hist same" : "hist");
  if (hHFm)
    hHFm->Draw((hEpAng || hHFp) ? "hist same" : "hist");
  if (hCalTrack)
    hCalTrack->Draw((hEpAng || hHFp || hHFm) ? "hist same" : "hist");
  TLegend *legOverlay = new TLegend(0.16, 0.70, 0.48, 0.88);
  legOverlay->SetBorderSize(0);
  legOverlay->SetFillStyle(0);
  if (hEpAng)
    legOverlay->AddEntry(hEpAng, "EPAng", "l");
  if (hHFp)
    legOverlay->AddEntry(hHFp, "HF+", "l");
  if (hHFm)
    legOverlay->AddEntry(hHFm, "HF-", "l");
  if (hCalTrack)
    legOverlay->AddEntry(hCalTrack, "CalTrack", "l");
  legOverlay->Draw();
  DrawPanelLabel("EP angle overlay", 0.60, 0.85);
  SaveCanvas(cEpOverlay, outDir, sample, "ep_overlay");

  TCanvas *cDelta = new TCanvas("cDelta", "folded EP angle differences", 1500, 500);
  cDelta->Divide(3, 1);
  cDelta->cd(1); SetupPad(); DrawTreeHist(t, deltaTrkHFp, "h_d_trk_hfp(100,-1.6,1.6)", validAll, "#Delta#Psi_{2}(trk,HF+) (rad)", "trk - HF+", kRed + 1, 0.68, 0.85);
  cDelta->cd(2); SetupPad(); DrawTreeHist(t, deltaTrkHFm, "h_d_trk_hfm(100,-1.6,1.6)", validAll, "#Delta#Psi_{2}(trk,HF-) (rad)", "trk - HF-", kBlue + 1, 0.68, 0.85);
  cDelta->cd(3); SetupPad(); DrawTreeHist(t, deltaHFpHFm, "h_d_hfp_hfm(100,-1.6,1.6)", validHF, "#Delta#Psi_{2}(HF+,HF-) (rad)", "HF+ - HF-", kGreen + 2, 0.66, 0.85);
  SaveCanvas(cDelta, outDir, sample, "angle_differences");

  TCanvas *cTrkCompare = new TCanvas("cTrkCompare", "calculated tracker vs input tracker EP", 1500, 500);
  cTrkCompare->Divide(3, 1);
  cTrkCompare->cd(1); SetupPad(); DrawTreeHist(t, "ep_flat", "h_ep_flat_cmp(100,-1.5,1.5)", validTrk, "#Psi_{2,trk}^{calc} (rad)", "Calculated tracker", kBlack, 0.58, 0.85);
  cTrkCompare->cd(2); SetupPad(); DrawTreeHist(t, "rpAng[3]", "h_ep_trackmid_input(100,-1.5,1.5)", "nEP>3", "#Psi_{2,trk}^{input} (rad)", "Input trackmid2", kMagenta + 1, 0.60, 0.85);
  cTrkCompare->cd(3); SetupPad(); DrawTreeHist(t, deltaCalcInputTrk, "h_d_calc_input_trk(100,-1.6,1.6)", "has_trk_ep && ep_flat>-90 && nEP>3", "#Delta#Psi_{2}(calc,input) (rad)", "calc - input track", kOrange + 7, 0.54, 0.85);
  SaveCanvas(cTrkCompare, outDir, sample, "trk_input_compare");

  TCanvas *cAngle2D = new TCanvas("cAngle2D", "2D EP angle comparisons", 1500, 500);
  cAngle2D->Divide(3, 1);
  cAngle2D->cd(1); SetupPad(); DrawTreeHist2D(t, "ep_hfp2:ep_flat", "h2_trk_hfp(80,-1.5,1.5,80,-1.5,1.5)", validAll, "#Psi_{2,trk}^{calc}", "#Psi_{2,HF+}", "trk vs HF+");
  cAngle2D->cd(2); SetupPad(); DrawTreeHist2D(t, "ep_hfm2:ep_flat", "h2_trk_hfm(80,-1.5,1.5,80,-1.5,1.5)", validAll, "#Psi_{2,trk}^{calc}", "#Psi_{2,HF-}", "trk vs HF-");
  cAngle2D->cd(3); SetupPad(); DrawTreeHist2D(t, "rpAng[3]:ep_flat", "h2_trk_input(80,-1.5,1.5,80,-1.5,1.5)", "has_trk_ep && ep_flat>-90 && nEP>3", "#Psi_{2,trk}^{calc}", "#Psi_{2,trk}^{input}", "calc vs input trk");
  SaveCanvas(cAngle2D, outDir, sample, "angle_2d");

  if (hCos2TrkHFp && hCos2TrkHFm && hCos2HFpHFm && hResolution3Sub)
  {
    const double meanTrkHFp = hResolution3Sub->GetBinContent(1);
    const double meanTrkHFm = hResolution3Sub->GetBinContent(2);
    const double meanHFpHFm = hResolution3Sub->GetBinContent(3);
    const double rTrk = hResolution3Sub->GetBinContent(4);
    const double entries = hResolution3Sub->GetBinContent(5);

    SetHistStyle(hCos2TrkHFp, kRed + 1, 20);
    SetHistStyle(hCos2TrkHFm, kBlue + 1, 21);
    SetHistStyle(hCos2HFpHFm, kGreen + 2, 22);

    TCanvas *cCorr = new TCanvas("cCorr", "3-subevent correlations", 900, 700);
    SetupPad();
    double maxY = std::max({hCos2TrkHFp->GetMaximum(), hCos2TrkHFm->GetMaximum(), hCos2HFpHFm->GetMaximum()});
    hCos2TrkHFp->SetTitle(";cos 2(#Psi_{A}-#Psi_{B});Entries");
    hCos2TrkHFp->SetMaximum(maxY > 0 ? maxY * 1.35 : 1.0);
    hCos2TrkHFp->Draw("hist");
    hCos2TrkHFm->Draw("hist same");
    hCos2HFpHFm->Draw("hist same");

    TLegend *legCorr = new TLegend(0.16, 0.70, 0.58, 0.88);
    legCorr->SetBorderSize(0);
    legCorr->SetFillStyle(0);
    legCorr->AddEntry(hCos2TrkHFp, Form("Trk-HF+  mean = %.4f", meanTrkHFp), "l");
    legCorr->AddEntry(hCos2TrkHFm, Form("Trk-HF-  mean = %.4f", meanTrkHFm), "l");
    legCorr->AddEntry(hCos2HFpHFm, Form("HF+-HF-  mean = %.4f", meanHFpHFm), "l");
    legCorr->Draw();
    latex.SetTextSize(0.040);
    latex.DrawLatex(0.63, 0.83, Form("R_{trk} = %.4f", rTrk));
    latex.DrawLatex(0.63, 0.77, Form("Entries = %.0f", entries));
    SaveCanvas(cCorr, outDir, sample, "cos2_correlations");

    TCanvas *cRes = new TCanvas("cRes", "3-subevent summary", 900, 650);
    SetupPad();
    cRes->SetBottomMargin(0.20);
    TH1 *hResDraw = (TH1 *)hResolution3Sub->Clone("hResDraw");
    hResDraw->SetStats(0);
    hResDraw->SetLineColor(kBlack);
    hResDraw->SetFillColor(kAzure - 9);
    hResDraw->SetBarWidth(0.65);
    hResDraw->SetBarOffset(0.15);
    hResDraw->SetTitle(";Quantity;Value");
    hResDraw->SetMinimum(0.0);
    hResDraw->SetMaximum(std::max(1.05, rTrk * 1.25));
    hResDraw->LabelsOption("v", "X");
    hResDraw->Draw("bar2 text0");
    SaveCanvas(cRes, outDir, sample, "resolution_summary");

    cout << std::setprecision(12)
         << "R_Trk = sqrt(" << meanTrkHFp << " * " << meanTrkHFm << " / " << meanHFpHFm << ") = " << rTrk << endl;
  }
  else
  {
    cout << "[WARN] 3-subevent histograms are missing; skipped correlation and summary plots" << endl;
  }

  cout << "Input: " << inputFile << endl;
  cout << "Sample: " << sample << endl;
  cout << "saved in: " << outDir << endl;
  f->Close();
}
