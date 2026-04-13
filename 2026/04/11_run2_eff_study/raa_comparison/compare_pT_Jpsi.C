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

TGraphErrors *makePtGraph(const double *edges, const double *y, const double *yerr, int nBins, bool zeroXErr)
{
  double x[8] = {0.0};
  double xerr[8] = {0.0};
  for (int i = 0; i < nBins; ++i)
  {
    x[i] = jpsi_raa::pt_center(edges, i);
    xerr[i] = zeroXErr ? 0.0 : jpsi_raa::pt_half_width(edges, i);
  }
  return new TGraphErrors(nBins, x, y, xerr, yerr);
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
                    double xMax, const char *outName)
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
  gCurrent->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gCurrent->GetXaxis()->CenterTitle();
  gCurrent->GetYaxis()->SetTitle("R_{AA}");
  gCurrent->GetYaxis()->CenterTitle();
  gCurrent->GetXaxis()->SetTitleSize(0.05);
  gCurrent->GetYaxis()->SetTitleSize(0.05);
  gCurrent->GetXaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetTitleOffset(1.20);
  gCurrent->GetXaxis()->SetLimits(0.0, xMax);
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

  TLine *line = new TLine(0.0, 1.0, xMax, 1.0);
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
                               double xMax, const char *outName)
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
  gCurrent->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gCurrent->GetXaxis()->CenterTitle();
  gCurrent->GetYaxis()->SetTitle("R_{AA}");
  gCurrent->GetYaxis()->CenterTitle();
  gCurrent->GetXaxis()->SetTitleSize(0.05);
  gCurrent->GetYaxis()->SetTitleSize(0.05);
  gCurrent->GetXaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetLabelSize(0.042);
  gCurrent->GetYaxis()->SetTitleOffset(1.20);
  gCurrent->GetXaxis()->SetLimits(0.0, xMax);
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

  TLine *line = new TLine(0.0, 1.0, xMax, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

} // namespace

void compare_pT_Jpsi(bool isSys = true)
{
  gROOT->Macro(kRootlogonPath);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("./figs", true);

  // Current Run-2 points: TnPL2L3 only.
  TGraphErrors *gRun2MidPr = makePtGraph(jpsi_raa::kPtMidBinEdges, jpsi_raa::kTnPL2L3PtMidPr, jpsi_raa::kTnPL2L3PtMidPrStat, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gRun2FwdPr = makePtGraph(jpsi_raa::kPtFwdBinEdges, jpsi_raa::kTnPL2L3PtFwdPr, jpsi_raa::kTnPL2L3PtFwdPrStat, jpsi_raa::kNPtFwd, !isSys);
  TGraphErrors *gRun2MidNp = makePtGraph(jpsi_raa::kPtMidBinEdges, jpsi_raa::kTnPL2L3PtMidNp, jpsi_raa::kTnPL2L3PtMidNpStat, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gRun2FwdNp = makePtGraph(jpsi_raa::kPtFwdBinEdges, jpsi_raa::kTnPL2L3PtFwdNp, jpsi_raa::kTnPL2L3PtFwdNpStat, jpsi_raa::kNPtFwd, !isSys);

  styleGraph(gRun2MidPr, kAzure + 2, 20);
  styleGraph(gRun2FwdPr, kAzure + 2, 20);
  styleGraph(gRun2MidNp, kRed + 1, 21);
  styleGraph(gRun2FwdNp, kRed + 1, 21);

  double xHinMidPrAbsY0p6[6] = {7.5, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xErrHinMidPrAbsY0p6[6] = {1.0, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yHinMidPrAbsY0p6[6] = {0.335, 0.347, 0.330, 0.350, 0.359, 0.427};
  double yErrHinMidPrAbsY0p6[6] = {0.019, 0.018, 0.015, 0.012, 0.017, 0.029};
  double ySysHinMidPrAbsY0p6[6] = {0.038, 0.025, 0.019, 0.018, 0.017, 0.022};
  double xHinMidPrAbsY1p6[5] = {7.75, 10.5, 13.5, 17.5, 25.0};
  double xErrHinMidPrAbsY1p6[5] = {1.25, 1.5, 1.5, 2.5, 5.0};
  double yHinMidPrAbsY1p6[5] = {0.353, 0.337, 0.366, 0.356, 0.473};
  double yErrHinMidPrAbsY1p6[5] = {0.008, 0.007, 0.010, 0.011, 0.021};
  double ySysHinMidPrAbsY1p6[5] = {0.035, 0.021, 0.021, 0.023, 0.033};
  double xHinMidNpAbsY0p6[6] = {7.5, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xErrHinMidNpAbsY0p6[6] = {1.0, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yHinMidNpAbsY0p6[6] = {0.476, 0.534, 0.432, 0.466, 0.456, 0.494};
  double yErrHinMidNpAbsY0p6[6] = {0.027, 0.028, 0.020, 0.016, 0.022, 0.034};
  double ySysHinMidNpAbsY0p6[6] = {0.057, 0.043, 0.028, 0.024, 0.022, 0.025};
  double xHinFwdPrAbsY1p8to2p4[10] = {3.75, 5.0, 6.0, 7.0, 8.0, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xErrHinFwdPrAbsY1p8to2p4[10] = {0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yHinFwdPrAbsY1p8to2p4[10] = {0.553, 0.365, 0.328, 0.334, 0.301, 0.355, 0.291, 0.317, 0.358, 0.298};
  double yErrHinFwdPrAbsY1p8to2p4[10] = {0.047, 0.024, 0.022, 0.020, 0.021, 0.026, 0.020, 0.018, 0.033, 0.042};
  double ySysHinFwdPrAbsY1p8to2p4[10] = {0.113, 0.050, 0.036, 0.028, 0.024, 0.027, 0.020, 0.019, 0.030, 0.034};
  double xHinFwdPrAbsY1p6to2p4[3] = {4.75, 9.25, 21.0};
  double xErrHinFwdPrAbsY1p6to2p4[3] = {1.75, 2.75, 9.0};
  double yHinFwdPrAbsY1p6to2p4[3] = {0.470, 0.329, 0.348};
  double yErrHinFwdPrAbsY1p6to2p4[3] = {0.017, 0.008, 0.014};
  double ySysHinFwdPrAbsY1p6to2p4[3] = {0.079, 0.029, 0.024};
  double xHinFwdNpAbsY1p8to2p4[10] = {3.75, 5.0, 6.0, 7.0, 8.0, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xErrHinFwdNpAbsY1p8to2p4[10] = {0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yHinFwdNpAbsY1p8to2p4[10] = {0.739, 0.635, 0.517, 0.576, 0.483, 0.393, 0.420, 0.425, 0.357, 0.391};
  double yErrHinFwdNpAbsY1p8to2p4[10] = {0.063, 0.042, 0.034, 0.034, 0.034, 0.029, 0.030, 0.025, 0.033, 0.056};
  double ySysHinFwdNpAbsY1p8to2p4[10] = {0.261, 0.111, 0.060, 0.060, 0.055, 0.037, 0.033, 0.028, 0.033, 0.036};

  TGraphErrors *gHinMidPrAbsY0p6 = new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY0p6, isSys ? xErrHinMidPrAbsY0p6 : nullptr, yErrHinMidPrAbsY0p6);
  TGraphErrors *gHinMidPrAbsY0p6Sys = isSys ? new TGraphErrors(6, xHinMidPrAbsY0p6, yHinMidPrAbsY0p6, xErrHinMidPrAbsY0p6, ySysHinMidPrAbsY0p6) : nullptr;
  TGraphErrors *gHinMidPrAbsY1p6 = new TGraphErrors(5, xHinMidPrAbsY1p6, yHinMidPrAbsY1p6, isSys ? xErrHinMidPrAbsY1p6 : nullptr, yErrHinMidPrAbsY1p6);
  TGraphErrors *gHinMidPrAbsY1p6Sys = isSys ? new TGraphErrors(5, xHinMidPrAbsY1p6, yHinMidPrAbsY1p6, xErrHinMidPrAbsY1p6, ySysHinMidPrAbsY1p6) : nullptr;
  TGraphErrors *gHinMidNpAbsY0p6 = new TGraphErrors(6, xHinMidNpAbsY0p6, yHinMidNpAbsY0p6, isSys ? xErrHinMidNpAbsY0p6 : nullptr, yErrHinMidNpAbsY0p6);
  TGraphErrors *gHinMidNpAbsY0p6Sys = isSys ? new TGraphErrors(6, xHinMidNpAbsY0p6, yHinMidNpAbsY0p6, xErrHinMidNpAbsY0p6, ySysHinMidNpAbsY0p6) : nullptr;
  TGraphErrors *gHinFwdPrAbsY1p8to2p4 = new TGraphErrors(10, xHinFwdPrAbsY1p8to2p4, yHinFwdPrAbsY1p8to2p4, isSys ? xErrHinFwdPrAbsY1p8to2p4 : nullptr, yErrHinFwdPrAbsY1p8to2p4);
  TGraphErrors *gHinFwdPrAbsY1p8to2p4Sys = isSys ? new TGraphErrors(10, xHinFwdPrAbsY1p8to2p4, yHinFwdPrAbsY1p8to2p4, xErrHinFwdPrAbsY1p8to2p4, ySysHinFwdPrAbsY1p8to2p4) : nullptr;
  TGraphErrors *gHinFwdPrAbsY1p6to2p4 = new TGraphErrors(3, xHinFwdPrAbsY1p6to2p4, yHinFwdPrAbsY1p6to2p4, isSys ? xErrHinFwdPrAbsY1p6to2p4 : nullptr, yErrHinFwdPrAbsY1p6to2p4);
  TGraphErrors *gHinFwdPrAbsY1p6to2p4Sys = isSys ? new TGraphErrors(3, xHinFwdPrAbsY1p6to2p4, yHinFwdPrAbsY1p6to2p4, xErrHinFwdPrAbsY1p6to2p4, ySysHinFwdPrAbsY1p6to2p4) : nullptr;
  TGraphErrors *gHinFwdNpAbsY1p8to2p4 = new TGraphErrors(10, xHinFwdNpAbsY1p8to2p4, yHinFwdNpAbsY1p8to2p4, isSys ? xErrHinFwdNpAbsY1p8to2p4 : nullptr, yErrHinFwdNpAbsY1p8to2p4);
  TGraphErrors *gHinFwdNpAbsY1p8to2p4Sys = isSys ? new TGraphErrors(10, xHinFwdNpAbsY1p8to2p4, yHinFwdNpAbsY1p8to2p4, xErrHinFwdNpAbsY1p8to2p4, ySysHinFwdNpAbsY1p8to2p4) : nullptr;

  styleGraph(gHinMidPrAbsY0p6, kBlack, 24);
  styleGraph(gHinMidPrAbsY1p6, kGreen + 3, 25);
  styleGraph(gHinMidNpAbsY0p6, kBlack, 25);
  styleGraph(gHinFwdPrAbsY1p8to2p4, kBlack, 24);
  styleGraph(gHinFwdPrAbsY1p6to2p4, kMagenta + 1, 26);
  styleGraph(gHinFwdNpAbsY1p8to2p4, kBlack, 25);
  if (gHinMidPrAbsY0p6Sys) styleSysGraph(gHinMidPrAbsY0p6Sys, kGray + 1);
  if (gHinMidPrAbsY1p6Sys) styleSysGraph(gHinMidPrAbsY1p6Sys, kGreen - 7);
  if (gHinMidNpAbsY0p6Sys) styleSysGraph(gHinMidNpAbsY0p6Sys, kGray + 1);
  if (gHinFwdPrAbsY1p8to2p4Sys) styleSysGraph(gHinFwdPrAbsY1p8to2p4Sys, kGray + 1);
  if (gHinFwdPrAbsY1p6to2p4Sys) styleSysGraph(gHinFwdPrAbsY1p6to2p4Sys, kMagenta - 9);
  if (gHinFwdNpAbsY1p8to2p4Sys) styleSysGraph(gHinFwdNpAbsY1p8to2p4Sys, kGray + 1);

  // HIN-16-025 reference points are hard-coded below to keep labels and binning explicit.
  drawComparisonWithTwoOlds(gRun2MidPr, gHinMidPrAbsY0p6, gHinMidPrAbsY0p6Sys, gHinMidPrAbsY1p6, gHinMidPrAbsY1p6Sys,
                 "PR (Cent. 0-90%, |y| < 1.6)",
                 "PR (HIN-16-025, Cent. 0-100%, |y| < 0.6)",
                 "PR (HIN-16-025, Cent. 0-100%, |y| < 1.6)",
                 40.0, "compare_mid_pT_Jpsi_PR");
  drawComparisonWithTwoOlds(gRun2FwdPr, gHinFwdPrAbsY1p8to2p4, gHinFwdPrAbsY1p8to2p4Sys, gHinFwdPrAbsY1p6to2p4, gHinFwdPrAbsY1p6to2p4Sys,
                 "PR (Cent. 0-90%, 1.6 < |y| < 2.4)",
                 "PR (HIN-16-025, Cent. 0-100%, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, Cent. 0-100%, 1.6 < |y| < 2.4)",
                 40.0, "compare_fwd_pT_Jpsi_PR");
  drawComparison(gRun2MidNp, gHinMidNpAbsY0p6, gHinMidNpAbsY0p6Sys,
                 "NP (Cent. 0-90%, |y| < 1.6)", "NP (HIN-16-025, Cent. 0-100%, |y| < 0.6)",
                 40.0, "compare_mid_pT_Jpsi_NP");
  drawComparison(gRun2FwdNp, gHinFwdNpAbsY1p8to2p4, gHinFwdNpAbsY1p8to2p4Sys,
                 "NP (Cent. 0-90%, 1.6 < |y| < 2.4)", "NP (HIN-16-025, Cent. 0-100%, 1.8 < |y| < 2.4)",
                 40.0, "compare_fwd_pT_Jpsi_NP");
}
