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
const char *kOldRootDir = "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/Macros/final_Results/roots";
const char *kSystRootPath = "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/Macros/syst_summary_Jpsi/syst_roots/total_syst.root";
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

TGraphErrors *makeSysCentGraph(const double *x, const double *xerr, const double *y, TH1D *hFrac, int nBins)
{
  double yLocal[8] = {0.0};
  double yerr[8] = {0.0};
  for (int i = 0; i < nBins; ++i)
  {
    yLocal[i] = y[nBins - 1 - i];
    yerr[i] = hFrac ? yLocal[i] * hFrac->GetBinContent(i + 1) : 0.0;
  }
  return new TGraphErrors(nBins, x, yLocal, xerr, yerr);
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

void drawTextBlock(const char *line1, const char *line2, TCanvas *c)
{
  (void)c;
  (void)line1;
  (void)line2;
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

void drawComparison(TGraphErrors *gNew, TGraphErrors *gNewSys,
                    TGraphErrors *gAlt,
                    TGraphErrors *gOld, TGraphErrors *gOldSys,
                    const char *legendNew, const char *legendAlt, const char *legendOld,
                    const char *line1, const char *line2,
                    const char *outName)
{
  (void)gNew;
  (void)gNewSys;
  (void)legendNew;
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gAlt->SetTitle("");
  gAlt->GetXaxis()->SetTitle("<N_{Part}>");
  gAlt->GetXaxis()->CenterTitle();
  gAlt->GetYaxis()->SetTitle("R_{AA}");
  gAlt->GetYaxis()->CenterTitle();
  gAlt->GetXaxis()->SetTitleSize(0.05);
  gAlt->GetYaxis()->SetTitleSize(0.05);
  gAlt->GetXaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetTitleOffset(1.20);
  gAlt->GetXaxis()->SetLimits(0.0, 400.0);
  gAlt->SetMinimum(0.0);
  gAlt->SetMaximum(1.44);
  gAlt->Draw("AP");
  gAlt->Draw("P SAME");
  gOld->Draw("P SAME");
  if (gOldSys)
    gOldSys->Draw("5");
  gAlt->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.73, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gAlt, legendAlt, "pe");
  leg->AddEntry(gOld, legendOld, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(line1, line2, c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

void drawComparisonWithTwoOlds(TGraphErrors *gNew, TGraphErrors *gNewSys,
                               TGraphErrors *gAlt,
                               TGraphErrors *gOldA, TGraphErrors *gOldASys,
                               TGraphErrors *gOldB, TGraphErrors *gOldBSys,
                               const char *legendNew, const char *legendAlt,
                               const char *legendOldA, const char *legendOldB,
                               const char *line1, const char *line2,
                               const char *outName)
{
  (void)gNew;
  (void)gNewSys;
  (void)legendNew;
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gAlt->SetTitle("");
  gAlt->GetXaxis()->SetTitle("<N_{Part}>");
  gAlt->GetXaxis()->CenterTitle();
  gAlt->GetYaxis()->SetTitle("R_{AA}");
  gAlt->GetYaxis()->CenterTitle();
  gAlt->GetXaxis()->SetTitleSize(0.05);
  gAlt->GetYaxis()->SetTitleSize(0.05);
  gAlt->GetXaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetTitleOffset(1.20);
  gAlt->GetXaxis()->SetLimits(0.0, 400.0);
  gAlt->SetMinimum(0.0);
  gAlt->SetMaximum(1.44);
  gAlt->Draw("AP");
  gAlt->Draw("P SAME");
  gOldA->Draw("P SAME");
  if (gOldASys)
    gOldASys->Draw("5");
  gOldB->Draw("P SAME");
  if (gOldBSys)
    gOldBSys->Draw("5");
  gAlt->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.73, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gAlt, legendAlt, "pe");
  leg->AddEntry(gOldA, legendOldA, "pe");
  leg->AddEntry(gOldB, legendOldB, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(line1, line2, c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

void drawComparisonWithThreeOlds(TGraphErrors *gNew, TGraphErrors *gNewSys,
                                 TGraphErrors *gAlt,
                                 TGraphErrors *gOldA, TGraphErrors *gOldASys,
                                 TGraphErrors *gOldB, TGraphErrors *gOldBSys,
                                 TGraphErrors *gOldC, TGraphErrors *gOldCSys,
                                 const char *legendNew, const char *legendAlt,
                                 const char *legendOldA, const char *legendOldB, const char *legendOldC,
                                 const char *line1, const char *line2,
                                 const char *outName)
{
  (void)gNew;
  (void)gNewSys;
  (void)legendNew;
  TCanvas *c = new TCanvas(outName, "", 800, 800);
  c->cd();
  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->Draw();
  pad->cd();

  gAlt->SetTitle("");
  gAlt->GetXaxis()->SetTitle("<N_{Part}>");
  gAlt->GetXaxis()->CenterTitle();
  gAlt->GetYaxis()->SetTitle("R_{AA}");
  gAlt->GetYaxis()->CenterTitle();
  gAlt->GetXaxis()->SetTitleSize(0.05);
  gAlt->GetYaxis()->SetTitleSize(0.05);
  gAlt->GetXaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetTitleOffset(1.20);
  gAlt->GetXaxis()->SetLimits(0.0, 400.0);
  gAlt->SetMinimum(0.0);
  gAlt->SetMaximum(1.44);
  gAlt->Draw("AP");
  gAlt->Draw("P SAME");
  gOldA->Draw("P SAME");
  if (gOldASys)
    gOldASys->Draw("5");
  gOldB->Draw("P SAME");
  if (gOldBSys)
    gOldBSys->Draw("5");
  gOldC->Draw("P SAME");
  if (gOldCSys)
    gOldCSys->Draw("5");
  gAlt->Draw("P SAME");

  TLegend *leg = new TLegend(0.42, 0.61, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0175);
  leg->SetEntrySeparation(0.010);
  leg->AddEntry(gAlt, legendAlt, "pe");
  leg->AddEntry(gOldA, legendOldA, "pe");
  leg->AddEntry(gOldB, legendOldB, "pe");
  leg->AddEntry(gOldC, legendOldC, "pe");
  leg->Draw();

  TLine *line = new TLine(0.0, 1.0, 400.0, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(line1, line2, c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

} // namespace

void compare_Npart_Jpsi(bool isSys = true)
{
  gROOT->Macro(kRootlogonPath);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("./figs", true);

  TFile *fSys = TFile::Open(kSystRootPath, "READ");
  TFile *fPrMidOld = TFile::Open(Form("%s/RAA_PR_JPsi_HIN_16_025_mid_Npart.root", kOldRootDir), "READ");
  TFile *fPrFwdOld = TFile::Open(Form("%s/RAA_PR_JPsi_HIN_16_025_fwd_Npart.root", kOldRootDir), "READ");
  TFile *fNpMidOld = TFile::Open(Form("%s/RAA_NP_JPsi_HIN_16_025_Npart.root", kOldRootDir), "READ");
  TFile *fNpFwdOld = TFile::Open(Form("%s/RAA_NP_JPsi_HIN_16_025_fwd_Npart.root", kOldRootDir), "READ");

  TH1D *hSysMidPr = fSys ? static_cast<TH1D *>(fSys->Get("mid_cent_PR")) : nullptr;
  TH1D *hSysFwdPr = fSys ? static_cast<TH1D *>(fSys->Get("fwd_cent_PR")) : nullptr;
  TH1D *hSysMidNp = fSys ? static_cast<TH1D *>(fSys->Get("mid_cent_NP")) : nullptr;
  TH1D *hSysFwdNp = fSys ? static_cast<TH1D *>(fSys->Get("fwd_cent_NP")) : nullptr;

  TGraphErrors *gMidPr = makeCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidPr, jpsi_raa::kCentMidPrErr, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gFwdPr = makeCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdPr, jpsi_raa::kCentFwdPrErr, jpsi_raa::kNCentFwd, !isSys);
  TGraphErrors *gMidNp = makeCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidNp, jpsi_raa::kCentMidNpErr, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gFwdNp = makeCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdNp, jpsi_raa::kCentFwdNpErr, jpsi_raa::kNCentFwd, !isSys);
  TGraphErrors *gMidPrL2L3 = makeCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidPrL2L3, jpsi_raa::kCentMidPrErrL2L3, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gFwdPrL2L3 = makeCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdPrL2L3, jpsi_raa::kCentFwdPrErrL2L3, jpsi_raa::kNCentFwd, !isSys);
  TGraphErrors *gMidNpL2L3 = makeCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidNpL2L3, jpsi_raa::kCentMidNpErrL2L3, jpsi_raa::kNCentMid, !isSys);
  TGraphErrors *gFwdNpL2L3 = makeCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdNpL2L3, jpsi_raa::kCentFwdNpErrL2L3, jpsi_raa::kNCentFwd, !isSys);

  TGraphErrors *gMidPrSys = isSys ? makeSysCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidPr, hSysMidPr, jpsi_raa::kNCentMid) : nullptr;
  TGraphErrors *gFwdPrSys = isSys ? makeSysCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdPr, hSysFwdPr, jpsi_raa::kNCentFwd) : nullptr;
  TGraphErrors *gMidNpSys = isSys ? makeSysCentGraph(jpsi_raa::kNpartMid, jpsi_raa::kNpartXErrMid, jpsi_raa::kCentMidNp, hSysMidNp, jpsi_raa::kNCentMid) : nullptr;
  TGraphErrors *gFwdNpSys = isSys ? makeSysCentGraph(jpsi_raa::kNpartFwd, jpsi_raa::kNpartXErrFwd, jpsi_raa::kCentFwdNp, hSysFwdNp, jpsi_raa::kNCentFwd) : nullptr;

  styleGraph(gMidPr, kBlue + 2, 20);
  styleGraph(gFwdPr, kBlue + 2, 20);
  styleGraph(gMidNp, kRed + 2, 21);
  styleGraph(gFwdNp, kRed + 2, 21);
  styleGraph(gMidPrL2L3, kAzure + 2, 20);
  styleGraph(gFwdPrL2L3, kAzure + 2, 20);
  styleGraph(gMidNpL2L3, kRed + 1, 21);
  styleGraph(gFwdNpL2L3, kRed + 1, 21);
  if (gMidPrSys) styleSysGraph(gMidPrSys, kBlue - 9);
  if (gFwdPrSys) styleSysGraph(gFwdPrSys, kBlue - 9);
  if (gMidNpSys) styleSysGraph(gMidNpSys, kRed - 9);
  if (gFwdNpSys) styleSysGraph(gFwdNpSys, kRed - 9);

  double xMidOld[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xMidOldErr[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yMidPrOld[6] = {0.736, 0.599, 0.456, 0.412, 0.332, 0.241};
  double yMidPrOldErr[6] = {0.043, 0.034, 0.021, 0.016, 0.012, 0.008};
  double yMidPrOldSys[6] = {0.092, 0.052, 0.033, 0.029, 0.021, 0.015};
  double yMidPrOldWide[6] = {0.704, 0.615, 0.461, 0.397, 0.330, 0.251};
  double yMidPrOldWideErr[6] = {0.025, 0.021, 0.013, 0.010, 0.007, 0.005};
  double yMidPrOldWideSys[6] = {0.087, 0.053, 0.035, 0.029, 0.024, 0.020};

  double xFwdPrOld[6] = {21.9, 86.9, 131.4, 189.2, 264.2, 358.8};
  double xFwdPrOldErr[6] = {4.3, 4.3, 4.3, 4.3, 4.3, 4.3};
  double yFwdPrOld[6] = {0.697, 0.566, 0.514, 0.564, 0.443, 0.328};
  double yFwdPrOldErr[6] = {0.064, 0.053, 0.042, 0.051, 0.037, 0.034};
  double yFwdPrOldSys[6] = {0.111, 0.080, 0.070, 0.082, 0.071, 0.060};
  double yFwdPrOldLow[6] = {0.672, 0.514, 0.484, 0.349, 0.270, 0.212};
  double yFwdPrOldLowErr[6] = {0.049, 0.039, 0.029, 0.020, 0.015, 0.012};
  double yFwdPrOldLowSys[6] = {0.084, 0.049, 0.038, 0.026, 0.020, 0.016};
  double xFwdPrOldWide[3] = {32.7, 160.3, 311.5};
  double xFwdPrOldWideErr[3] = {4.0, 4.0, 4.0};
  double yFwdPrOldWide[3] = {0.719, 0.558, 0.372};
  double yFwdPrOldWideErr[3] = {0.024, 0.016, 0.010};
  double yFwdPrOldWideSys[3] = {0.091, 0.067, 0.058};

  double yMidNpOld[6] = {0.705, 0.621, 0.538, 0.531, 0.456, 0.356};
  double yMidNpOldErr[6] = {0.042, 0.036, 0.025, 0.021, 0.016, 0.011};
  double yMidNpOldSys[6] = {0.088, 0.052, 0.038, 0.036, 0.025, 0.020};

  double yFwdNpOld[6] = {0.737, 0.637, 0.600, 0.490, 0.424, 0.358};
  double yFwdNpOldErr[6] = {0.055, 0.048, 0.035, 0.028, 0.023, 0.019};
  double yFwdNpOldSys[6] = {0.094, 0.073, 0.049, 0.039, 0.034, 0.029};
  double yFwdNpOldLow[6] = {0.894, 0.899, 0.809, 0.626, 0.638, 0.560};
  double yFwdNpOldLowErr[6] = {0.082, 0.084, 0.067, 0.057, 0.052, 0.058};
  double yFwdNpOldLowSys[6] = {0.140, 0.137, 0.133, 0.126, 0.129, 0.170};

  TGraphErrors *gMidPrOld = new TGraphErrors(6, xMidOld, yMidPrOld, isSys ? xMidOldErr : nullptr, yMidPrOldErr);
  TGraphErrors *gMidPrOldSys = isSys ? new TGraphErrors(6, xMidOld, yMidPrOld, xMidOldErr, yMidPrOldSys) : nullptr;
  TGraphErrors *gMidPrOldWide = new TGraphErrors(6, xMidOld, yMidPrOldWide, isSys ? xMidOldErr : nullptr, yMidPrOldWideErr);
  TGraphErrors *gMidPrOldWideSys = isSys ? new TGraphErrors(6, xMidOld, yMidPrOldWide, xMidOldErr, yMidPrOldWideSys) : nullptr;
  TGraphErrors *gFwdPrOld = new TGraphErrors(6, xFwdPrOld, yFwdPrOld, isSys ? xFwdPrOldErr : nullptr, yFwdPrOldErr);
  TGraphErrors *gFwdPrOldSys = isSys ? new TGraphErrors(6, xFwdPrOld, yFwdPrOld, xFwdPrOldErr, yFwdPrOldSys) : nullptr;
  TGraphErrors *gFwdPrOldLow = new TGraphErrors(6, xFwdPrOld, yFwdPrOldLow, isSys ? xFwdPrOldErr : nullptr, yFwdPrOldLowErr);
  TGraphErrors *gFwdPrOldLowSys = isSys ? new TGraphErrors(6, xFwdPrOld, yFwdPrOldLow, xFwdPrOldErr, yFwdPrOldLowSys) : nullptr;
  TGraphErrors *gFwdPrOldWide = new TGraphErrors(3, xFwdPrOldWide, yFwdPrOldWide, isSys ? xFwdPrOldWideErr : nullptr, yFwdPrOldWideErr);
  TGraphErrors *gFwdPrOldWideSys = isSys ? new TGraphErrors(3, xFwdPrOldWide, yFwdPrOldWide, xFwdPrOldWideErr, yFwdPrOldWideSys) : nullptr;
  TGraphErrors *gMidNpOld = new TGraphErrors(6, xMidOld, yMidNpOld, isSys ? xMidOldErr : nullptr, yMidNpOldErr);
  TGraphErrors *gMidNpOldSys = isSys ? new TGraphErrors(6, xMidOld, yMidNpOld, xMidOldErr, yMidNpOldSys) : nullptr;
  TGraphErrors *gFwdNpOld = new TGraphErrors(6, xMidOld, yFwdNpOld, isSys ? xMidOldErr : nullptr, yFwdNpOldErr);
  TGraphErrors *gFwdNpOldSys = isSys ? new TGraphErrors(6, xMidOld, yFwdNpOld, xMidOldErr, yFwdNpOldSys) : nullptr;
  TGraphErrors *gFwdNpOldLow = new TGraphErrors(6, xMidOld, yFwdNpOldLow, isSys ? xMidOldErr : nullptr, yFwdNpOldLowErr);
  TGraphErrors *gFwdNpOldLowSys = isSys ? new TGraphErrors(6, xMidOld, yFwdNpOldLow, xMidOldErr, yFwdNpOldLowSys) : nullptr;

  styleGraph(gMidPrOld, kBlack, 24);
  styleGraph(gMidPrOldWide, kGreen + 3, 25);
  styleGraph(gFwdPrOld, kBlack, 24);
  styleGraph(gFwdPrOldLow, kGreen + 3, 25);
  styleGraph(gFwdPrOldWide, kMagenta + 1, 26);
  styleGraph(gMidNpOld, kBlack, 25);
  styleGraph(gFwdNpOld, kBlack, 25);
  styleGraph(gFwdNpOldLow, kGreen + 3, 24);
  if (gMidPrOldSys) styleSysGraph(gMidPrOldSys, kGray + 1);
  if (gMidPrOldWideSys) styleSysGraph(gMidPrOldWideSys, kGreen - 7);
  if (gFwdPrOldSys) styleSysGraph(gFwdPrOldSys, kGray + 1);
  if (gFwdPrOldLowSys) styleSysGraph(gFwdPrOldLowSys, kGreen - 7);
  if (gFwdPrOldWideSys) styleSysGraph(gFwdPrOldWideSys, kMagenta - 9);
  if (gMidNpOldSys) styleSysGraph(gMidNpOldSys, kGray + 1);
  if (gFwdNpOldSys) styleSysGraph(gFwdNpOldSys, kGray + 1);
  if (gFwdNpOldLowSys) styleSysGraph(gFwdNpOldLowSys, kGreen - 7);

  drawComparisonWithTwoOlds(gMidPr, gMidPrSys, gMidPrL2L3, gMidPrOld, gMidPrOldSys, gMidPrOldWide, gMidPrOldWideSys,
                 "", "PR (6.5 < p_{T} < 40, |y| < 1.6)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 50, |y| < 0.6)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 30, |y| < 1.6)",
                 "6.5 < p_{T} < 40 GeV/c", "|y| < 1.6", "compare_mid_Npart_Jpsi_PR");
  drawComparisonWithThreeOlds(gFwdPr, gFwdPrSys, gFwdPrL2L3, gFwdPrOld, gFwdPrOldSys, gFwdPrOldLow, gFwdPrOldLowSys, gFwdPrOldWide, gFwdPrOldWideSys,
                 "", "PR (3.5 < p_{T} < 40, 1.6 < |y| < 2.4)",
                 "PR (HIN-16-025, 6.5 < p_{T} < 50, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, 3.0 < p_{T} < 6.5, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, 3.0 < p_{T} < 30, 1.6 < |y| < 2.4)",
                 "3.5 < p_{T} < 40 GeV/c", "1.6 < |y| < 2.4", "compare_fwd_Npart_Jpsi_PR");
  drawComparison(gMidNp, gMidNpSys, gMidNpL2L3, gMidNpOld, gMidNpOldSys,
                 "", "NP (6.5 < p_{T} < 40, |y| < 1.6)", "NP (HIN-16-025, 6.5 < p_{T} < 50, |y| < 0.6)",
                 "6.5 < p_{T} < 40 GeV/c", "|y| < 1.6", "compare_mid_Npart_Jpsi_NP");
  drawComparisonWithTwoOlds(gFwdNp, gFwdNpSys, gFwdNpL2L3, gFwdNpOld, gFwdNpOldSys, gFwdNpOldLow, gFwdNpOldLowSys,
                 "", "NP (3.5 < p_{T} < 40, 1.6 < |y| < 2.4)",
                 "NP (HIN-16-025, 6.5 < p_{T} < 50, 1.8 < |y| < 2.4)",
                 "NP (HIN-16-025, 3.0 < p_{T} < 6.5, 1.8 < |y| < 2.4)",
                 "3.5 < p_{T} < 40 GeV/c", "1.6 < |y| < 2.4", "compare_fwd_Npart_Jpsi_NP");
}
