#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "jpsi_raa_values.h"



namespace
{
const char *kRootlogonPath = "/data/users/pjgwak/input_files/rootlogon.C";

TGraphErrors *makeCentGraph(const double *x, const double *xerr, const double *y, const double *yerr, int nBins, bool zeroXErr)
{
  double xerrLocal[8] = {0.0};
  double yLocal[8] = {0.0};
  double yerrLocal[8] = {0.0};
  for (int i = 0; i < nBins; ++i)
  {
    xerrLocal[i] = zeroXErr ? 0.0 : xerr[i];
    yLocal[i] = y[nBins - 1 - i];
    yerrLocal[i] = yerr[nBins - 1 - i];
  }
  return new TGraphErrors(nBins, x, yLocal, xerrLocal, yerrLocal);
}

void styleGraph(TGraphErrors *g, Color_t color, Style_t markerStyle)
{
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(1.6);
  g->SetLineWidth(3);
}

void styleSysGraph(TGraphErrors *g, Color_t color)
{
  g->SetLineColor(color);
  g->SetFillColorAlpha(color, 0.14);
  g->SetMarkerSize(0);
}

void drawTextBlock(TCanvas *c)
{
  (void)c;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.96, 0.935, "PbPb 5.02 TeV");
  latex.SetTextAlign(11);
  latex.SetTextFont(72);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.19, 0.935, "CMS");
  latex.SetTextFont(42);
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.285, 0.935, "Internal");
}

void drawComparison(TGraphErrors *gCurrent,
                    TGraphErrors *gOld, TGraphErrors *gOldSys,
                    const char *legendCurrent, const char *legendOld,
                    const char *outName)
{
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gCurrent->SetTitle("");
  gCurrent->GetXaxis()->SetTitle("<N_{Part}>");
  gCurrent->GetXaxis()->CenterTitle();
  gCurrent->GetYaxis()->SetTitle("R_{AA}");
  gCurrent->GetYaxis()->CenterTitle();
  gCurrent->GetXaxis()->SetTitleSize(0.05);
  gCurrent->GetYaxis()->SetTitleSize(0.05);
  gCurrent->GetXaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetTitleOffset(1.20);
  gCurrent->GetXaxis()->SetLimits(0.0, 400.0);
  gCurrent->SetMinimum(0.0);
  gCurrent->SetMaximum(1.44);
  gCurrent->Draw("AP");
  gCurrent->Draw("P SAME");
  gOld->Draw("P SAME");
  if (gOldSys)
    gOldSys->Draw("5");
  gCurrent->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.73, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gCurrent, legendCurrent, "pe");
  leg->AddEntry(gOld, legendOld, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

void drawComparisonWithTwoOlds(TGraphErrors *gCurrent,
                               TGraphErrors *gOldA, TGraphErrors *gOldASys,
                               TGraphErrors *gOldB, TGraphErrors *gOldBSys,
                               const char *legendCurrent,
                               const char *legendOldA, const char *legendOldB,
                               const char *outName)
{
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gCurrent->SetTitle("");
  gCurrent->GetXaxis()->SetTitle("<N_{Part}>");
  gCurrent->GetXaxis()->CenterTitle();
  gCurrent->GetYaxis()->SetTitle("R_{AA}");
  gCurrent->GetYaxis()->CenterTitle();
  gCurrent->GetXaxis()->SetTitleSize(0.05);
  gCurrent->GetYaxis()->SetTitleSize(0.05);
  gCurrent->GetXaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetTitleOffset(1.20);
  gCurrent->GetXaxis()->SetLimits(0.0, 400.0);
  gCurrent->SetMinimum(0.0);
  gCurrent->SetMaximum(1.44);
  gCurrent->Draw("AP");
  gCurrent->Draw("P SAME");
  gOldA->Draw("P SAME");
  if (gOldASys)
    gOldASys->Draw("5");
  gOldB->Draw("P SAME");
  if (gOldBSys)
    gOldBSys->Draw("5");
  gCurrent->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.73, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gCurrent, legendCurrent, "pe");
  leg->AddEntry(gOldA, legendOldA, "pe");
  leg->AddEntry(gOldB, legendOldB, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

void drawComparisonWithThreeOlds(TGraphErrors *gCurrent,
                                 TGraphErrors *gOldA, TGraphErrors *gOldASys,
                                 TGraphErrors *gOldB, TGraphErrors *gOldBSys,
                                 TGraphErrors *gOldC, TGraphErrors *gOldCSys,
                                 const char *legendCurrent,
                                 const char *legendOldA, const char *legendOldB, const char *legendOldC,
                                 const char *outName)
{
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gCurrent->SetTitle("");
  gCurrent->GetXaxis()->SetTitle("<N_{Part}>");
  gCurrent->GetXaxis()->CenterTitle();
  gCurrent->GetYaxis()->SetTitle("R_{AA}");
  gCurrent->GetYaxis()->CenterTitle();
  gCurrent->GetXaxis()->SetTitleSize(0.05);
  gCurrent->GetYaxis()->SetTitleSize(0.05);
  gCurrent->GetXaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetTitleOffset(1.20);
  gCurrent->GetXaxis()->SetLimits(0.0, 400.0);
  gCurrent->SetMinimum(0.0);
  gCurrent->SetMaximum(1.44);
  gCurrent->Draw("AP");
  gCurrent->Draw("P SAME");
  gOldA->Draw("P SAME");
  if (gOldASys)
    gOldASys->Draw("5");
  gOldB->Draw("P SAME");
  if (gOldBSys)
    gOldBSys->Draw("5");
  gOldC->Draw("P SAME");
  if (gOldCSys)
    gOldCSys->Draw("5");
  gCurrent->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.61, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gCurrent, legendCurrent, "pe");
  leg->AddEntry(gOldA, legendOldA, "pe");
  leg->AddEntry(gOldB, legendOldB, "pe");
  leg->AddEntry(gOldC, legendOldC, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

} // namespace

void compare_Npart_Jpsi(bool isSys = true)
{
  gROOT->Macro(kRootlogonPath);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("./figs", true);

  // Current Run-2 points: TnPL2L3 only.
  TGraphErrors *gRun2MidPr = makeCentGraph(jpsi_raa::kMidNpart, jpsi_raa::kMidNpartXErr, jpsi_raa::kTnPL2L3CentMidPr, jpsi_raa::kTnPL2L3CentMidPrStat, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gRun2FwdPr = makeCentGraph(jpsi_raa::kFwdNpart, jpsi_raa::kFwdNpartXErr, jpsi_raa::kTnPL2L3CentFwdPr, jpsi_raa::kTnPL2L3CentFwdPrStat, jpsi_raa::kNCentFwd, !isSys);
  TGraphErrors *gRun2MidNp = makeCentGraph(jpsi_raa::kMidNpart, jpsi_raa::kMidNpartXErr, jpsi_raa::kTnPL2L3CentMidNp, jpsi_raa::kTnPL2L3CentMidNpStat, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gRun2FwdNp = makeCentGraph(jpsi_raa::kFwdNpart, jpsi_raa::kFwdNpartXErr, jpsi_raa::kTnPL2L3CentFwdNp, jpsi_raa::kTnPL2L3CentFwdNpStat, jpsi_raa::kNCentFwd, !isSys);

  styleGraph(gRun2MidPr, kAzure + 2, 20);
  styleGraph(gRun2FwdPr, kAzure + 2, 20);
  styleGraph(gRun2MidNp, kRed + 1, 21);
  styleGraph(gRun2FwdNp, kRed + 1, 21);

  double xHinMidPrAbsY0p6[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xErrHinMidPrAbsY0p6[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yHinMidPrAbsY0p6[6] = {0.736, 0.599, 0.456, 0.412, 0.332, 0.241};
  double yErrHinMidPrAbsY0p6[6] = {0.043, 0.034, 0.021, 0.016, 0.012, 0.008};
  double ySysHinMidPrAbsY0p6[6] = {0.092, 0.052, 0.033, 0.029, 0.021, 0.015};
  double yHinMidPrAbsY1p6[6] = {0.704, 0.615, 0.461, 0.397, 0.330, 0.251};
  double yErrHinMidPrAbsY1p6[6] = {0.025, 0.021, 0.013, 0.010, 0.007, 0.005};
  double ySysHinMidPrAbsY1p6[6] = {0.087, 0.053, 0.035, 0.029, 0.024, 0.020};

  double xHinFwdPrAbsY1p8to2p4[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xErrHinFwdPrAbsY1p8to2p4[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yHinFwdPrAbsY1p8to2p4[6] = {0.697, 0.566, 0.514, 0.564, 0.443, 0.328};
  double yErrHinFwdPrAbsY1p8to2p4[6] = {0.064, 0.053, 0.042, 0.051, 0.037, 0.034};
  double ySysHinFwdPrAbsY1p8to2p4[6] = {0.111, 0.080, 0.070, 0.082, 0.071, 0.060};
  double yHinFwdPrLowPtAbsY1p8to2p4[6] = {0.672, 0.514, 0.484, 0.349, 0.270, 0.212};
  double yErrHinFwdPrLowPtAbsY1p8to2p4[6] = {0.049, 0.039, 0.029, 0.020, 0.015, 0.012};
  double ySysHinFwdPrLowPtAbsY1p8to2p4[6] = {0.084, 0.049, 0.038, 0.026, 0.020, 0.016};
  double xHinFwdPrAbsY1p6to2p4[3] = {32.7, 160.3, 311.5};
  double xErrHinFwdPrAbsY1p6to2p4[3] = {4.0, 4.0, 4.0};
  double yHinFwdPrAbsY1p6to2p4[3] = {0.719, 0.558, 0.372};
  double yErrHinFwdPrAbsY1p6to2p4[3] = {0.024, 0.016, 0.010};
  double ySysHinFwdPrAbsY1p6to2p4[3] = {0.091, 0.067, 0.058};

  double xHinMidNpAbsY0p6[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xErrHinMidNpAbsY0p6[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yHinMidNpAbsY0p6[6] = {0.705, 0.621, 0.538, 0.531, 0.456, 0.356};
  double yErrHinMidNpAbsY0p6[6] = {0.042, 0.036, 0.025, 0.021, 0.016, 0.011};
  double ySysHinMidNpAbsY0p6[6] = {0.088, 0.052, 0.038, 0.036, 0.025, 0.020};

  double xHinFwdNpAbsY1p8to2p4[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xErrHinFwdNpAbsY1p8to2p4[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yHinFwdNpAbsY1p8to2p4[6] = {0.737, 0.637, 0.600, 0.490, 0.424, 0.358};
  double yErrHinFwdNpAbsY1p8to2p4[6] = {0.055, 0.048, 0.035, 0.028, 0.023, 0.019};
  double ySysHinFwdNpAbsY1p8to2p4[6] = {0.094, 0.073, 0.049, 0.039, 0.034, 0.029};
  double yHinFwdNpLowPtAbsY1p8to2p4[6] = {0.894, 0.899, 0.809, 0.626, 0.638, 0.560};
  double yErrHinFwdNpLowPtAbsY1p8to2p4[6] = {0.082, 0.084, 0.067, 0.057, 0.052, 0.058};
  double ySysHinFwdNpLowPtAbsY1p8to2p4[6] = {0.140, 0.137, 0.133, 0.126, 0.129, 0.170};

  TGraphErrors *gHinMidPrAbsY0p6 = new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY0p6, isSys ? xErrHinMidPrAbsY0p6 : nullptr, yErrHinMidPrAbsY0p6);
  TGraphErrors *gHinMidPrAbsY0p6Sys = isSys ? new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY0p6, xErrHinMidPrAbsY0p6, ySysHinMidPrAbsY0p6) : nullptr;
  TGraphErrors *gHinMidPrAbsY1p6 = new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY1p6, isSys ? xErrHinMidPrAbsY0p6 : nullptr, yErrHinMidPrAbsY1p6);
  TGraphErrors *gHinMidPrAbsY1p6Sys = isSys ? new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY1p6, xErrHinMidPrAbsY0p6, ySysHinMidPrAbsY1p6) : nullptr;
  TGraphErrors *gHinFwdPrAbsY1p8to2p4 = new TGraphErrors(6, xHinFwdPrAbsY1p8to2p4, yHinFwdPrAbsY1p8to2p4, isSys ? xErrHinFwdPrAbsY1p8to2p4 : nullptr, yErrHinFwdPrAbsY1p8to2p4);
  TGraphErrors *gHinFwdPrAbsY1p8to2p4Sys = isSys ? new TGraphErrors(6, xHinFwdPrAbsY1p8to2p4, yHinFwdPrAbsY1p8to2p4, xErrHinFwdPrAbsY1p8to2p4, ySysHinFwdPrAbsY1p8to2p4) : nullptr;
  TGraphErrors *gHinFwdPrLowPtAbsY1p8to2p4 = new TGraphErrors(6, xHinFwdPrAbsY1p8to2p4, yHinFwdPrLowPtAbsY1p8to2p4, isSys ? xErrHinFwdPrAbsY1p8to2p4 : nullptr, yErrHinFwdPrLowPtAbsY1p8to2p4);
  TGraphErrors *gHinFwdPrLowPtAbsY1p8to2p4Sys = isSys ? new TGraphErrors(6, xHinFwdPrAbsY1p8to2p4, yHinFwdPrLowPtAbsY1p8to2p4, xErrHinFwdPrAbsY1p8to2p4, ySysHinFwdPrLowPtAbsY1p8to2p4) : nullptr;
  TGraphErrors *gHinFwdPrAbsY1p6to2p4 = new TGraphErrors(3, xHinFwdPrAbsY1p6to2p4, yHinFwdPrAbsY1p6to2p4, isSys ? xErrHinFwdPrAbsY1p6to2p4 : nullptr, yErrHinFwdPrAbsY1p6to2p4);
  TGraphErrors *gHinFwdPrAbsY1p6to2p4Sys = isSys ? new TGraphErrors(3, xHinFwdPrAbsY1p6to2p4, yHinFwdPrAbsY1p6to2p4, xErrHinFwdPrAbsY1p6to2p4, ySysHinFwdPrAbsY1p6to2p4) : nullptr;
  TGraphErrors *gHinMidNpAbsY0p6 = new TGraphErrors(6, xHinMidNpAbsY0p6, yHinMidNpAbsY0p6, isSys ? xErrHinMidNpAbsY0p6 : nullptr, yErrHinMidNpAbsY0p6);
  TGraphErrors *gHinMidNpAbsY0p6Sys = isSys ? new TGraphErrors(6, xHinMidNpAbsY0p6, yHinMidNpAbsY0p6, xErrHinMidNpAbsY0p6, ySysHinMidNpAbsY0p6) : nullptr;
  TGraphErrors *gHinFwdNpAbsY1p8to2p4 = new TGraphErrors(6, xHinFwdNpAbsY1p8to2p4, yHinFwdNpAbsY1p8to2p4, isSys ? xErrHinFwdNpAbsY1p8to2p4 : nullptr, yErrHinFwdNpAbsY1p8to2p4);
  TGraphErrors *gHinFwdNpAbsY1p8to2p4Sys = isSys ? new TGraphErrors(6, xHinFwdNpAbsY1p8to2p4, yHinFwdNpAbsY1p8to2p4, xErrHinFwdNpAbsY1p8to2p4, ySysHinFwdNpAbsY1p8to2p4) : nullptr;
  TGraphErrors *gHinFwdNpLowPtAbsY1p8to2p4 = new TGraphErrors(6, xHinFwdNpAbsY1p8to2p4, yHinFwdNpLowPtAbsY1p8to2p4, isSys ? xErrHinFwdNpAbsY1p8to2p4 : nullptr, yErrHinFwdNpLowPtAbsY1p8to2p4);
  TGraphErrors *gHinFwdNpLowPtAbsY1p8to2p4Sys = isSys ? new TGraphErrors(6, xHinFwdNpAbsY1p8to2p4, yHinFwdNpLowPtAbsY1p8to2p4, xErrHinFwdNpAbsY1p8to2p4, ySysHinFwdNpLowPtAbsY1p8to2p4) : nullptr;

  styleGraph(gHinMidPrAbsY0p6, kBlack, 24);
  styleGraph(gHinMidPrAbsY1p6, kGreen + 3, 25);
  styleGraph(gHinFwdPrAbsY1p8to2p4, kBlack, 24);
  styleGraph(gHinFwdPrLowPtAbsY1p8to2p4, kGreen + 3, 25);
  styleGraph(gHinFwdPrAbsY1p6to2p4, kMagenta + 1, 26);
  styleGraph(gHinMidNpAbsY0p6, kBlack, 25);
  styleGraph(gHinFwdNpAbsY1p8to2p4, kBlack, 25);
  styleGraph(gHinFwdNpLowPtAbsY1p8to2p4, kGreen + 3, 24);
  if (gHinMidPrAbsY0p6Sys) styleSysGraph(gHinMidPrAbsY0p6Sys, kGray + 1);
  if (gHinMidPrAbsY1p6Sys) styleSysGraph(gHinMidPrAbsY1p6Sys, kGreen - 7);
  if (gHinFwdPrAbsY1p8to2p4Sys) styleSysGraph(gHinFwdPrAbsY1p8to2p4Sys, kGray + 1);
  if (gHinFwdPrLowPtAbsY1p8to2p4Sys) styleSysGraph(gHinFwdPrLowPtAbsY1p8to2p4Sys, kGreen - 7);
  if (gHinFwdPrAbsY1p6to2p4Sys) styleSysGraph(gHinFwdPrAbsY1p6to2p4Sys, kMagenta - 9);
  if (gHinMidNpAbsY0p6Sys) styleSysGraph(gHinMidNpAbsY0p6Sys, kGray + 1);
  if (gHinFwdNpAbsY1p8to2p4Sys) styleSysGraph(gHinFwdNpAbsY1p8to2p4Sys, kGray + 1);
  if (gHinFwdNpLowPtAbsY1p8to2p4Sys) styleSysGraph(gHinFwdNpLowPtAbsY1p8to2p4Sys, kGreen - 7);

  // HIN-16-025 reference points are hard-coded below to keep each legacy binning explicit.
  drawComparisonWithTwoOlds(gRun2MidPr, gHinMidPrAbsY0p6, gHinMidPrAbsY0p6Sys, gHinMidPrAbsY1p6, gHinMidPrAbsY1p6Sys,
                 "PR (6.5 < p_{T} < 40, |y| < 1.6)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 50, |y| < 0.6)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 30, |y| < 1.6)",
                 "compare_mid_Npart_Jpsi_PR");
  drawComparisonWithThreeOlds(gRun2FwdPr, gHinFwdPrAbsY1p8to2p4, gHinFwdPrAbsY1p8to2p4Sys, gHinFwdPrLowPtAbsY1p8to2p4, gHinFwdPrLowPtAbsY1p8to2p4Sys, gHinFwdPrAbsY1p6to2p4, gHinFwdPrAbsY1p6to2p4Sys,
                 "PR (3.5 < p_{T} < 40, 1.6 < |y| < 2.4)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 50, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, 3.0 < p_{T} < 6.5, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, 3.0 < p_{T} < 30, 1.6 < |y| < 2.4)",
                 "compare_fwd_Npart_Jpsi_PR");
  drawComparison(gRun2MidNp, gHinMidNpAbsY0p6, gHinMidNpAbsY0p6Sys,
                 "NP (6.5 < p_{T} < 40, |y| < 1.6)", "NP (HIN-16-025, 6.5 < p_{T} < 50, |y| < 0.6)",
                 "compare_mid_Npart_Jpsi_NP");
  drawComparisonWithTwoOlds(gRun2FwdNp, gHinFwdNpAbsY1p8to2p4, gHinFwdNpAbsY1p8to2p4Sys, gHinFwdNpLowPtAbsY1p8to2p4, gHinFwdNpLowPtAbsY1p8to2p4Sys,
                 "NP (3.5 < p_{T} < 40, 1.6 < |y| < 2.4)",
                 "NP (HIN-16-025, 6.5 < p_{T} < 50, 1.8 < |y| < 2.4)",
                 "NP (HIN-16-025, 3.0 < p_{T} < 6.5, 1.8 < |y| < 2.4)",
                 "compare_fwd_Npart_Jpsi_NP");
}
