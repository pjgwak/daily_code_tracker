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

TGraphErrors *makeSysPtGraph(const double *edges, const double *y, TH1D *hFrac, int nBins)
{
  double x[8] = {0.0};
  double xerr[8] = {0.0};
  double yerr[8] = {0.0};
  for (int i = 0; i < nBins; ++i)
  {
    x[i] = jpsi_raa::pt_center(edges, i);
    xerr[i] = jpsi_raa::pt_half_width(edges, i);
    yerr[i] = hFrac ? y[i] * hFrac->GetBinContent(i + 1) : 0.0;
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
                    double xMax, const char *outName)
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
  gAlt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gAlt->GetXaxis()->CenterTitle();
  gAlt->GetYaxis()->SetTitle("R_{AA}");
  gAlt->GetYaxis()->CenterTitle();
  gAlt->GetXaxis()->SetTitleSize(0.05);
  gAlt->GetYaxis()->SetTitleSize(0.05);
  gAlt->GetXaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetTitleOffset(1.20);
  gAlt->GetXaxis()->SetLimits(0.0, xMax);
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

  TLine *line = new TLine(0.0, 1.0, xMax, 1.0);
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
                               double xMax, const char *outName)
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
  gAlt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gAlt->GetXaxis()->CenterTitle();
  gAlt->GetYaxis()->SetTitle("R_{AA}");
  gAlt->GetYaxis()->CenterTitle();
  gAlt->GetXaxis()->SetTitleSize(0.05);
  gAlt->GetYaxis()->SetTitleSize(0.05);
  gAlt->GetXaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetLabelSize(0.042);
  gAlt->GetYaxis()->SetTitleOffset(1.20);
  gAlt->GetXaxis()->SetLimits(0.0, xMax);
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

  TLine *line = new TLine(0.0, 1.0, xMax, 1.0);
  line->SetLineStyle(7);
  line->Draw();
  drawTextBlock(line1, line2, c);
  c->SaveAs(Form("./figs/%s.pdf", outName));
}

} // namespace

void compare_pT_Jpsi(bool isSys = true)
{
  gROOT->Macro(kRootlogonPath);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gSystem->mkdir("./figs", true);

  TFile *fSys = TFile::Open(kSystRootPath, "READ");
  TFile *fPrMidOld = TFile::Open(Form("%s/RAA_PR_Jpsi_HIN_16_025_mid_pT.root", kOldRootDir), "READ");
  TFile *fPrFwdOld = TFile::Open(Form("%s/RAA_PR_Jpsi_HIN_16_025_fwd_pT.root", kOldRootDir), "READ");
  TFile *fNpMidOld = TFile::Open(Form("%s/RAA_NP_Jpsi_HIN_16_025_mid_pT.root", kOldRootDir), "READ");
  TFile *fNpFwdOld = TFile::Open(Form("%s/RAA_NP_Jpsi_HIN_16_025_fwd_pT.root", kOldRootDir), "READ");

  TH1D *hSysMidPr = fSys ? static_cast<TH1D *>(fSys->Get("mid_pt_PR")) : nullptr;
  TH1D *hSysFwdPr = fSys ? static_cast<TH1D *>(fSys->Get("fwd_pt_PR")) : nullptr;
  TH1D *hSysMidNp = fSys ? static_cast<TH1D *>(fSys->Get("mid_pt_NP")) : nullptr;
  TH1D *hSysFwdNp = fSys ? static_cast<TH1D *>(fSys->Get("fwd_pt_NP")) : nullptr;

  TGraphErrors *gMidPr = makePtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidPr, jpsi_raa::kPtMidPrErr, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gFwdPr = makePtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdPr, jpsi_raa::kPtFwdPrErr, jpsi_raa::kNPtFwd, !isSys);
  TGraphErrors *gMidNp = makePtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidNp, jpsi_raa::kPtMidNpErr, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gFwdNp = makePtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdNp, jpsi_raa::kPtFwdNpErr, jpsi_raa::kNPtFwd, !isSys);
  TGraphErrors *gMidPrL2L3 = makePtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidPrL2L3, jpsi_raa::kPtMidPrErrL2L3, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gFwdPrL2L3 = makePtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdPrL2L3, jpsi_raa::kPtFwdPrErrL2L3, jpsi_raa::kNPtFwd, !isSys);
  TGraphErrors *gMidNpL2L3 = makePtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidNpL2L3, jpsi_raa::kPtMidNpErrL2L3, jpsi_raa::kNPtMid, !isSys);
  TGraphErrors *gFwdNpL2L3 = makePtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdNpL2L3, jpsi_raa::kPtFwdNpErrL2L3, jpsi_raa::kNPtFwd, !isSys);

  TGraphErrors *gMidPrSys = isSys ? makeSysPtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidPr, hSysMidPr, jpsi_raa::kNPtMid) : nullptr;
  TGraphErrors *gFwdPrSys = isSys ? makeSysPtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdPr, hSysFwdPr, jpsi_raa::kNPtFwd) : nullptr;
  TGraphErrors *gMidNpSys = isSys ? makeSysPtGraph(jpsi_raa::kPtMidEdges, jpsi_raa::kPtMidNp, hSysMidNp, jpsi_raa::kNPtMid) : nullptr;
  TGraphErrors *gFwdNpSys = isSys ? makeSysPtGraph(jpsi_raa::kPtFwdEdges, jpsi_raa::kPtFwdNp, hSysFwdNp, jpsi_raa::kNPtFwd) : nullptr;

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

  double xMidOld[6] = {7.5, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xMidOldErr[6] = {1.0, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yMidPrOld[6] = {0.335, 0.347, 0.330, 0.350, 0.359, 0.427};
  double yMidPrOldErr[6] = {0.019, 0.018, 0.015, 0.012, 0.017, 0.029};
  double yMidPrOldSys[6] = {0.038, 0.025, 0.019, 0.018, 0.017, 0.022};
  double xMidPrOldWide[5] = {7.75, 10.5, 13.5, 17.5, 25.0};
  double xMidPrOldWideErr[5] = {1.25, 1.5, 1.5, 2.5, 5.0};
  double yMidPrOldWide[5] = {0.353, 0.337, 0.366, 0.356, 0.473};
  double yMidPrOldWideErr[5] = {0.008, 0.007, 0.010, 0.011, 0.021};
  double yMidPrOldWideSys[5] = {0.035, 0.021, 0.021, 0.023, 0.033};
  double yMidNpOld[6] = {0.476, 0.534, 0.432, 0.466, 0.456, 0.494};
  double yMidNpOldErr[6] = {0.027, 0.028, 0.020, 0.016, 0.022, 0.034};
  double yMidNpOldSys[6] = {0.057, 0.043, 0.028, 0.024, 0.022, 0.025};
  double xFwdPrOld[10] = {3.75, 5.0, 6.0, 7.0, 8.0, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xFwdPrOldErr[10] = {0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yFwdPrOld[10] = {0.553, 0.365, 0.328, 0.334, 0.301, 0.355, 0.291, 0.317, 0.358, 0.298};
  double yFwdPrOldErr[10] = {0.047, 0.024, 0.022, 0.020, 0.021, 0.026, 0.020, 0.018, 0.033, 0.042};
  double yFwdPrOldSys[10] = {0.113, 0.050, 0.036, 0.028, 0.024, 0.027, 0.020, 0.019, 0.030, 0.034};
  double xFwdPrOldWide[3] = {4.75, 9.25, 21.0};
  double xFwdPrOldWideErr[3] = {1.75, 2.75, 9.0};
  double yFwdPrOldWide[3] = {0.470, 0.329, 0.348};
  double yFwdPrOldWideErr[3] = {0.017, 0.008, 0.014};
  double yFwdPrOldWideSys[3] = {0.079, 0.029, 0.024};
  double xMidNpOld[6] = {7.5, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xMidNpOldErr[6] = {1.0, 0.5, 0.75, 2.0, 2.5, 5.0};
  double xFwdNpOld[10] = {3.75, 5.0, 6.0, 7.0, 8.0, 9.0, 10.25, 13.0, 17.5, 25.0};
  double xFwdNpOldErr[10] = {0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 2.0, 2.5, 5.0};
  double yFwdNpOld[10] = {0.739, 0.635, 0.517, 0.576, 0.483, 0.393, 0.420, 0.425, 0.357, 0.391};
  double yFwdNpOldErr[10] = {0.063, 0.042, 0.034, 0.034, 0.034, 0.029, 0.030, 0.025, 0.033, 0.056};
  double yFwdNpOldSys[10] = {0.261, 0.111, 0.060, 0.060, 0.055, 0.037, 0.033, 0.028, 0.033, 0.036};

  TGraphErrors *gMidPrOld = new TGraphErrors(6, xMidOld, yMidPrOld, isSys ? xMidOldErr : nullptr, yMidPrOldErr);
  TGraphErrors *gMidPrOldSys = isSys ? new TGraphErrors(6, xMidOld, yMidPrOld, xMidOldErr, yMidPrOldSys) : nullptr;
  TGraphErrors *gMidPrOldWide = new TGraphErrors(5, xMidPrOldWide, yMidPrOldWide, isSys ? xMidPrOldWideErr : nullptr, yMidPrOldWideErr);
  TGraphErrors *gMidPrOldWideSys = isSys ? new TGraphErrors(5, xMidPrOldWide, yMidPrOldWide, xMidPrOldWideErr, yMidPrOldWideSys) : nullptr;
  TGraphErrors *gFwdPrOld = new TGraphErrors(10, xFwdPrOld, yFwdPrOld, isSys ? xFwdPrOldErr : nullptr, yFwdPrOldErr);
  TGraphErrors *gFwdPrOldSys = isSys ? new TGraphErrors(10, xFwdPrOld, yFwdPrOld, xFwdPrOldErr, yFwdPrOldSys) : nullptr;
  TGraphErrors *gFwdPrOldWide = new TGraphErrors(3, xFwdPrOldWide, yFwdPrOldWide, isSys ? xFwdPrOldWideErr : nullptr, yFwdPrOldWideErr);
  TGraphErrors *gFwdPrOldWideSys = isSys ? new TGraphErrors(3, xFwdPrOldWide, yFwdPrOldWide, xFwdPrOldWideErr, yFwdPrOldWideSys) : nullptr;
  TGraphErrors *gMidNpOld = new TGraphErrors(6, xMidNpOld, yMidNpOld, isSys ? xMidNpOldErr : nullptr, yMidNpOldErr);
  TGraphErrors *gMidNpOldSys = isSys ? new TGraphErrors(6, xMidNpOld, yMidNpOld, xMidNpOldErr, yMidNpOldSys) : nullptr;
  TGraphErrors *gFwdNpOld = new TGraphErrors(10, xFwdNpOld, yFwdNpOld, isSys ? xFwdNpOldErr : nullptr, yFwdNpOldErr);
  TGraphErrors *gFwdNpOldSys = isSys ? new TGraphErrors(10, xFwdNpOld, yFwdNpOld, xFwdNpOldErr, yFwdNpOldSys) : nullptr;

  styleGraph(gMidPrOld, kBlack, 24);
  styleGraph(gMidPrOldWide, kGreen + 3, 25);
  styleGraph(gFwdPrOld, kBlack, 24);
  styleGraph(gFwdPrOldWide, kMagenta + 1, 26);
  styleGraph(gMidNpOld, kBlack, 25);
  styleGraph(gFwdNpOld, kBlack, 25);
  if (gMidPrOldSys) styleSysGraph(gMidPrOldSys, kGray + 1);
  if (gMidPrOldWideSys) styleSysGraph(gMidPrOldWideSys, kGreen - 7);
  if (gFwdPrOldSys) styleSysGraph(gFwdPrOldSys, kGray + 1);
  if (gFwdPrOldWideSys) styleSysGraph(gFwdPrOldWideSys, kMagenta - 9);
  if (gMidNpOldSys) styleSysGraph(gMidNpOldSys, kGray + 1);
  if (gFwdNpOldSys) styleSysGraph(gFwdNpOldSys, kGray + 1);

  drawComparisonWithTwoOlds(gMidPr, gMidPrSys, gMidPrL2L3, gMidPrOld, gMidPrOldSys, gMidPrOldWide, gMidPrOldWideSys,
                 "", "PR (Cent. 0-90%, |y| < 1.6)",
                 "PR (HIN-16-025, Cent. 0-100%, |y| < 0.6)",
                 "PR (HIN-16-025, Cent. 0-100%, |y| < 1.6)",
                 "Cent. 0-90%", "|y| < 1.6", 40.0, "compare_mid_pT_Jpsi_PR");
  drawComparisonWithTwoOlds(gFwdPr, gFwdPrSys, gFwdPrL2L3, gFwdPrOld, gFwdPrOldSys, gFwdPrOldWide, gFwdPrOldWideSys,
                 "", "PR (Cent. 0-90%, 1.6 < |y| < 2.4)",
                 "PR (HIN-16-025, Cent. 0-100%, 1.8 < |y| < 2.4)",
                 "PR (HIN-16-025, Cent. 0-100%, 1.6 < |y| < 2.4)",
                 "Cent. 0-90%", "1.6 < |y| < 2.4", 40.0, "compare_fwd_pT_Jpsi_PR");
  drawComparison(gMidNp, gMidNpSys, gMidNpL2L3, gMidNpOld, gMidNpOldSys,
                 "", "NP (Cent. 0-90%, |y| < 1.6)", "NP (HIN-16-025, Cent. 0-100%, |y| < 0.6)",
                 "Cent. 0-90%", "|y| < 1.6", 40.0, "compare_mid_pT_Jpsi_NP");
  drawComparison(gFwdNp, gFwdNpSys, gFwdNpL2L3, gFwdNpOld, gFwdNpOldSys,
                 "", "NP (Cent. 0-90%, 1.6 < |y| < 2.4)", "NP (HIN-16-025, Cent. 0-100%, 1.8 < |y| < 2.4)",
                 "Cent. 0-90%", "1.6 < |y| < 2.4", 40.0, "compare_fwd_pT_Jpsi_NP");
}
