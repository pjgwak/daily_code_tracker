#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TBox.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>

namespace
{
TString macroDirectory()
{
  TString macroPath = __FILE__;
  if (!gSystem->IsAbsoluteFileName(macroPath))
    macroPath = TString(gSystem->WorkingDirectory()) + "/" + macroPath;
  return gSystem->DirName(macroPath);
}

TString outputPath(const char *outputFileName)
{
  return macroDirectory() + "/" + outputFileName;
}
} // namespace

void raw_yield_mid()
{
  TStopwatch executionTimer;
  executionTimer.Start();
  struct ExecutionTimerGuard {
    TStopwatch *stopwatch;
    explicit ExecutionTimerGuard(TStopwatch *timer) : stopwatch(timer) {}
    ~ExecutionTimerGuard() {
      stopwatch->Stop();
      stopwatch->Print();
    }
  } executionTimerGuard(&executionTimer);

  // =====  CMS plot style =====
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // ===== OO raw yields =====
  // OO, ctau subrange method (updated bins)
  const int nominalBinCount = 8;

  double nominalPtLow[nominalBinCount] = {7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  double nominalPtHigh[nominalBinCount] = {8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 30.0};

  // nominal raw yields
  double nominalYields[nominalBinCount] = {
      6253.0, 5021.8, 3571.2, 4250.2,
      2208.2, 1185.7, 1081.6, 602.8};

  // statistical errors
  double nominalStatErrors[nominalBinCount] = {
      95.1, 81.1, 67.1, 72.3,
      51.4, 37.2, 35.6, 27.3};

  // systematics not provided
  double nominalSysErrors[nominalBinCount];
  for (int binIndex = 0; binIndex < nominalBinCount; ++binIndex)
    nominalSysErrors[binIndex] = 0.0; // sys=0

  double nominalPtCenters[nominalBinCount], nominalPtHalfWidths[nominalBinCount];
  for (int binIndex = 0; binIndex < nominalBinCount; ++binIndex)
  {
    const double binWidth = nominalPtHigh[binIndex] - nominalPtLow[binIndex];
    nominalPtCenters[binIndex] = 0.5 * (nominalPtLow[binIndex] + nominalPtHigh[binIndex]);
    nominalPtHalfWidths[binIndex] = 0.5 * binWidth;
    nominalYields[binIndex] /= binWidth;
    nominalStatErrors[binIndex] /= binWidth;
    nominalSysErrors[binIndex] /= binWidth;
  }

  // ===== OO |y|<1.6 (2d fit) =====
  // y = 0.00-1.60
  const int fitBinCount = 8;
  double fitPtLow[fitBinCount] = {7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  double fitPtHigh[fitBinCount] = {8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 30.0};

  double fitYields[fitBinCount] = {6231.41197, 4948.96553, 3533.20990, 4211.54399, 2222.78147, 1156.57651, 1058.15308, 623.83447};
  double fitStatErrors[fitBinCount] = {65.98198, 55.39028, 46.58958, 50.12010, 36.09343, 25.75811, 24.75737, 18.96578};

  double fitSysErrors[fitBinCount];
  for (int binIndex = 0; binIndex < fitBinCount; ++binIndex)
    fitSysErrors[binIndex] = 0.0; // sys=0

  double fitPtCenters[fitBinCount], fitPtHalfWidths[fitBinCount];
  for (int binIndex = 0; binIndex < fitBinCount; ++binIndex)
  {
    const double binWidth = fitPtHigh[binIndex] - fitPtLow[binIndex];
    fitPtCenters[binIndex] = 0.5 * (fitPtLow[binIndex] + fitPtHigh[binIndex]);
    fitPtHalfWidths[binIndex] = 0.5 * (fitPtHigh[binIndex] - fitPtLow[binIndex]);
    fitYields[binIndex] /= binWidth;
    fitStatErrors[binIndex] /= binWidth;
    fitSysErrors[binIndex] /= binWidth;
  }

  // ===== Determine the displayed yield range =====
  double minimumYield = 1.0e30;
  double maximumYield = 0.0;
  for (int binIndex = 0; binIndex < nominalBinCount; ++binIndex)
  {
    const double lowerYield = nominalYields[binIndex] - nominalStatErrors[binIndex];
    const double upperYield = nominalYields[binIndex] + nominalStatErrors[binIndex];
    if (lowerYield > 0.0 && lowerYield < minimumYield)
      minimumYield = lowerYield;
    if (upperYield > maximumYield)
      maximumYield = upperYield;
  }
  for (int binIndex = 0; binIndex < fitBinCount; ++binIndex)
  {
    const double lowerYield = fitYields[binIndex] - fitStatErrors[binIndex];
    const double upperYield = fitYields[binIndex] + fitStatErrors[binIndex];
    if (lowerYield > 0.0 && lowerYield < minimumYield)
      minimumYield = lowerYield;
    if (upperYield > maximumYield)
      maximumYield = upperYield;
  }

  // ===== Graphs =====
  // nominal: systematic band
  TGraphAsymmErrors *nominalSystematicGraph = new TGraphAsymmErrors(nominalBinCount, nominalPtCenters, nominalYields, nominalPtHalfWidths, nominalPtHalfWidths, nominalSysErrors, nominalSysErrors);
  nominalSystematicGraph->SetFillStyle(0);
  nominalSystematicGraph->SetLineColor(kRed + 1);
  nominalSystematicGraph->SetLineWidth(2);

  // nominal: statistical errors
  TGraphErrors *nominalStatisticalGraph = new TGraphErrors(nominalBinCount, nominalPtCenters, nominalYields, nominalPtHalfWidths, nominalStatErrors);
  nominalStatisticalGraph->SetMarkerStyle(24);
  nominalStatisticalGraph->SetMarkerSize(1.7);
  nominalStatisticalGraph->SetMarkerColor(kRed + 1);
  nominalStatisticalGraph->SetLineColor(kRed + 1);
  nominalStatisticalGraph->SetLineWidth(3);

  // OO |y|<1.6: sys (currently 0)
  TGraphAsymmErrors *fitSystematicGraph = new TGraphAsymmErrors(fitBinCount, fitPtCenters, fitYields, fitPtHalfWidths, fitPtHalfWidths, fitSysErrors, fitSysErrors);
  fitSystematicGraph->SetFillColorAlpha(kAzure + 2, 0.20);
  fitSystematicGraph->SetFillStyle(1001);
  fitSystematicGraph->SetLineColor(kAzure + 2);
  fitSystematicGraph->SetLineWidth(3);

  // OO |y|<1.6: stat
  TGraphErrors *fitStatisticalGraph = new TGraphErrors(fitBinCount, fitPtCenters, fitYields, fitPtHalfWidths, fitStatErrors);
  fitStatisticalGraph->SetMarkerStyle(20);
  fitStatisticalGraph->SetMarkerSize(1.7);
  fitStatisticalGraph->SetMarkerColor(kAzure + 2);
  fitStatisticalGraph->SetLineColor(kAzure + 2);
  fitStatisticalGraph->SetLineWidth(3);

  // ===== Draw =====
  TCanvas *canvas = new TCanvas("canvas", "OO raw yields", 800, 800);
  canvas->SetTitle("");
  TPad plotPad("plotPad", "plotPad", 0.0, 0.0, 1.0, 1.0);
  plotPad.SetTopMargin(0.08);
  plotPad.SetBottomMargin(0.13);
  plotPad.SetLogy();
  plotPad.Draw();
  plotPad.cd();
  gStyle->SetOptStat(0);
  nominalSystematicGraph->SetTitle("");
  nominalSystematicGraph->Draw("A2");
  nominalSystematicGraph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  nominalSystematicGraph->GetYaxis()->SetTitle("dN/dp_{T}");
  nominalSystematicGraph->GetXaxis()->CenterTitle();
  nominalSystematicGraph->GetYaxis()->CenterTitle();
  nominalSystematicGraph->GetXaxis()->SetTitleSize(0.05);
  nominalSystematicGraph->GetYaxis()->SetTitleSize(0.05);
  nominalSystematicGraph->GetXaxis()->SetLabelSize(0.042);
  nominalSystematicGraph->GetYaxis()->SetLabelSize(0.042);
  nominalSystematicGraph->GetYaxis()->SetTitleOffset(1.25);
  nominalSystematicGraph->GetXaxis()->SetLimits(nominalPtLow[0], fitPtHigh[fitBinCount - 1]);
  nominalSystematicGraph->SetMinimum(0.8 * minimumYield);
  nominalSystematicGraph->SetMaximum(4 * maximumYield);

  nominalSystematicGraph->Draw("A2");
  fitSystematicGraph->Draw("2 SAME");
  nominalStatisticalGraph->Draw("P SAME");
  fitStatisticalGraph->Draw("P SAME");

  // ------------------------------------------------------------------
  // labels (match ctau_pr.C style/positions)
  // ------------------------------------------------------------------
  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.032);
    label.SetTextFont(42);
    label.SetTextAlign(31);
    label.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.04);
    label.SetTextFont(72);
    label.DrawLatex(0.19, 0.935, "CMS Internal");
  }
  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    // label.SetTextColor(kRed + 1);
    label.DrawLatex(0.19, 0.865, "|y| < 1.6");
  }

  TLegend *legend = new TLegend(0.7, 0.8, 0.90, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.024);
  legend->SetEntrySeparation(0.01);
  legend->AddEntry(fitStatisticalGraph, "OO (2D fit)", "lp");
  legend->AddEntry(nominalStatisticalGraph, "OO (Nominal)", "lp");
  legend->Draw();

  canvas->SaveAs(outputPath("raw_yield_OO_mid.pdf"));
  // canvas->SaveAs(outputPath("plot_OO_mid.png"));

  // ===== Draw the raw yields and their 2D-fit / nominal ratio =====
  double yieldRatios[fitBinCount], yieldRatioStatErrors[fitBinCount];
  for (int binIndex = 0; binIndex < fitBinCount; ++binIndex)
  {
    yieldRatios[binIndex] = fitYields[binIndex] / nominalYields[binIndex];
    yieldRatioStatErrors[binIndex] = yieldRatios[binIndex] *
        std::sqrt(std::pow(fitStatErrors[binIndex] / fitYields[binIndex], 2) +
                  std::pow(nominalStatErrors[binIndex] / nominalYields[binIndex], 2));
  }

  TGraphAsymmErrors *ratioCanvasNominalSystematicGraph =
      static_cast<TGraphAsymmErrors *>(nominalSystematicGraph->Clone("ratioCanvasNominalSystematicGraphMid"));
  TGraphAsymmErrors *ratioCanvasFitSystematicGraph =
      static_cast<TGraphAsymmErrors *>(fitSystematicGraph->Clone("ratioCanvasFitSystematicGraphMid"));
  TGraphErrors *ratioCanvasNominalStatisticalGraph =
      static_cast<TGraphErrors *>(nominalStatisticalGraph->Clone("ratioCanvasNominalStatisticalGraphMid"));
  TGraphErrors *ratioCanvasFitStatisticalGraph =
      static_cast<TGraphErrors *>(fitStatisticalGraph->Clone("ratioCanvasFitStatisticalGraphMid"));

  TGraphErrors *yieldRatioGraph = new TGraphErrors(
      fitBinCount, fitPtCenters, yieldRatios, fitPtHalfWidths, yieldRatioStatErrors);
  yieldRatioGraph->SetName("yieldRatioGraphMid");
  yieldRatioGraph->SetTitle("");
  yieldRatioGraph->SetMarkerStyle(20);
  yieldRatioGraph->SetMarkerSize(1.7);
  yieldRatioGraph->SetMarkerColor(kAzure + 2);
  yieldRatioGraph->SetLineColor(kAzure + 2);
  yieldRatioGraph->SetLineWidth(3);

  TCanvas *ratioCanvas = new TCanvas("ratioCanvasMid", "OO raw yields and ratio", 800, 800);
  ratioCanvas->SetTitle("");

  TPad *ratioTopPad = new TPad("ratioTopPadMid", "ratioTopPadMid", 0.0, 0.25, 1.0, 1.0);
  ratioTopPad->SetBottomMargin(0.00001);
  ratioTopPad->SetTopMargin(0.08);
  ratioTopPad->SetLogy();
  ratioTopPad->Draw();
  ratioTopPad->cd();

  ratioCanvasNominalSystematicGraph->Draw("A2");
  ratioCanvasNominalSystematicGraph->GetXaxis()->SetTitle("");
  ratioCanvasNominalSystematicGraph->GetXaxis()->SetLabelSize(0.0);
  ratioCanvasNominalSystematicGraph->GetYaxis()->SetTitle("dN/dp_{T}");
  ratioCanvasNominalSystematicGraph->GetYaxis()->CenterTitle();
  ratioCanvasNominalSystematicGraph->GetYaxis()->SetTitleSize(0.05);
  ratioCanvasNominalSystematicGraph->GetYaxis()->SetLabelSize(0.042);
  ratioCanvasNominalSystematicGraph->GetYaxis()->SetTitleOffset(1.25);
  ratioCanvasNominalSystematicGraph->GetXaxis()->SetLimits(nominalPtLow[0], fitPtHigh[fitBinCount - 1]);
  ratioCanvasNominalSystematicGraph->SetMinimum(0.8 * minimumYield);
  ratioCanvasNominalSystematicGraph->SetMaximum(4 * maximumYield);
  ratioCanvasFitSystematicGraph->Draw("2 SAME");
  ratioCanvasNominalStatisticalGraph->Draw("P SAME");
  ratioCanvasFitStatisticalGraph->Draw("P SAME");

  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.032);
    label.SetTextFont(42);
    label.SetTextAlign(31);
    label.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.04);
    label.SetTextFont(72);
    label.DrawLatex(0.19, 0.935, "CMS Internal");
  }
  {
    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.03);
    label.SetTextFont(42);
    label.DrawLatex(0.19, 0.865, "|y| < 1.6");
  }

  TLegend *ratioCanvasLegend = new TLegend(0.7, 0.8, 0.90, 0.90);
  ratioCanvasLegend->SetBorderSize(0);
  ratioCanvasLegend->SetFillStyle(0);
  ratioCanvasLegend->SetTextSize(0.024);
  ratioCanvasLegend->SetEntrySeparation(0.01);
  ratioCanvasLegend->AddEntry(ratioCanvasFitStatisticalGraph, "OO (2D fit)", "lp");
  ratioCanvasLegend->AddEntry(ratioCanvasNominalStatisticalGraph, "OO (Nominal)", "lp");
  ratioCanvasLegend->Draw();

  ratioCanvas->cd();
  TPad *ratioBottomPad = new TPad("ratioBottomPadMid", "ratioBottomPadMid", 0.0, 0.0, 1.0, 0.25);
  ratioBottomPad->SetTopMargin(0.00001);
  ratioBottomPad->SetBottomMargin(0.4);
  ratioBottomPad->Draw();
  ratioBottomPad->cd();

  yieldRatioGraph->Draw("AP");
  yieldRatioGraph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  yieldRatioGraph->GetXaxis()->CenterTitle();
  yieldRatioGraph->GetXaxis()->SetTitleSize(0.15);
  yieldRatioGraph->GetXaxis()->SetLabelSize(0.10);
  yieldRatioGraph->GetXaxis()->SetLimits(nominalPtLow[0], fitPtHigh[fitBinCount - 1]);
  yieldRatioGraph->GetYaxis()->SetTitle("2D/Nom.");
  yieldRatioGraph->GetYaxis()->CenterTitle();
  yieldRatioGraph->GetYaxis()->SetNdivisions(505);
  yieldRatioGraph->GetYaxis()->SetTitleSize(0.12);
  yieldRatioGraph->GetYaxis()->SetLabelSize(0.10);
  yieldRatioGraph->GetYaxis()->SetTitleOffset(0.55);
  yieldRatioGraph->SetMinimum(0.85);
  yieldRatioGraph->SetMaximum(1.15);

  TLine ratioReferenceLine(nominalPtLow[0], 1.0, fitPtHigh[fitBinCount - 1], 1.0);
  ratioReferenceLine.SetLineStyle(2);
  ratioReferenceLine.SetLineWidth(2);
  ratioReferenceLine.Draw("SAME");
  yieldRatioGraph->Draw("P SAME");

  ratioCanvas->SaveAs(outputPath("raw_yield_OO_mid_ratio.pdf"));
}