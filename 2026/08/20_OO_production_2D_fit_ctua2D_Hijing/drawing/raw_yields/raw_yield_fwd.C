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

void raw_yield_fwd()
{
  TStopwatch executionTimer;
  executionTimer.Start();
  struct ExecutionTimerGuard
  {
    TStopwatch *stopwatch;
    explicit ExecutionTimerGuard(TStopwatch *timer) : stopwatch(timer) {}
    ~ExecutionTimerGuard()
    {
      stopwatch->Stop();
      stopwatch->Print();
    }
  } executionTimerGuard(&executionTimer);

  // =====  CMS plot style =====
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // ===== OO raw yields (12 bins) =====
  // ctau subrange method
  const int nominalBinCount = 12;

  double nominalPtLow[nominalBinCount] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                              7.0, 8.0, 9.0, 10.0, 12.0, 14.0};

  double nominalPtHigh[nominalBinCount] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                               8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  // nominal raw yields
  double nominalYields[nominalBinCount] = {
      25932.8, 24533.1, 19598.2, 13779.8,
      10321.4, 6640.5, 4467.7, 2796.3,
      1771.4, 2060.3, 943.6, 952.2};

  // statistical errors
  double nominalStatErrors[nominalBinCount] = {
      660.2, 474.7, 363.4, 260.2,
      191.1, 134.8, 99.4, 72.3,
      56.6, 59.0, 39.3, 38.0};

  // systematics not provided
  double nominalSysErrors[nominalBinCount] = {
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0};

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

  // ===== OO 1.6 < |y| < 2.4 (2d fit) =====
  // y = 1.60-2.40
  const int fitBinCount = 12;
  double fitPtLow[fitBinCount] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0};
  double fitPtHigh[fitBinCount] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  double fitYields[fitBinCount] = {24599.43971, 26167.35345, 19606.33220, 14111.38178, 10079.75989, 6589.37522, 4368.08358, 2780.10135, 1728.71869, 2024.21631, 896.68853, 923.16679};
  double fitStatErrors[fitBinCount] = {330.66030, 270.30127, 207.43002, 162.53041, 115.63698, 85.06409, 64.41203, 49.29955, 38.67589, 42.76550, 25.90626, 25.99896};

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
  nominalStatisticalGraph->SetMarkerSize(1.5);
  nominalStatisticalGraph->SetMarkerColor(kRed + 1);
  nominalStatisticalGraph->SetLineColor(kRed + 1);
  nominalStatisticalGraph->SetLineWidth(2);

  // OO 1.6 < |y| < 2.4: sys (currently 0)
  TGraphAsymmErrors *fitSystematicGraph = new TGraphAsymmErrors(fitBinCount, fitPtCenters, fitYields, fitPtHalfWidths, fitPtHalfWidths, fitSysErrors, fitSysErrors);
  fitSystematicGraph->SetFillStyle(0); // no fill
  fitSystematicGraph->SetLineColor(kAzure + 2);
  fitSystematicGraph->SetLineWidth(2);

  // OO 1.6 < |y| < 2.4: stat
  TGraphErrors *fitStatisticalGraph = new TGraphErrors(fitBinCount, fitPtCenters, fitYields, fitPtHalfWidths, fitStatErrors);
  fitStatisticalGraph->SetMarkerStyle(20);
  fitStatisticalGraph->SetMarkerSize(1.4);
  fitStatisticalGraph->SetMarkerColor(kAzure + 2);
  fitStatisticalGraph->SetLineColor(kAzure + 2);
  fitStatisticalGraph->SetLineWidth(2);

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
  nominalSystematicGraph->GetXaxis()->SetLimits(nominalPtLow[0], nominalPtHigh[nominalBinCount - 1]);
  nominalSystematicGraph->SetMinimum(0.8 * minimumYield);
  nominalSystematicGraph->SetMaximum(4 * maximumYield);

  nominalSystematicGraph->Draw("A2");
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
    label.DrawLatex(0.19, 0.865, "1.6 < |y| < 2.4");
  }

  TLegend *legend = new TLegend(0.7, 0.8, 0.90, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.024);
  legend->SetEntrySeparation(0.01);
  legend->AddEntry(fitStatisticalGraph, "OO (2D fit)", "lp");
  legend->AddEntry(nominalStatisticalGraph, "OO (Nominal)", "lp");
  legend->Draw();

  canvas->SaveAs(outputPath("raw_yield_OO_fwd.pdf"));
  // canvas->SaveAs(outputPath("plot_OO_fwd.png"));
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
      static_cast<TGraphAsymmErrors *>(nominalSystematicGraph->Clone("ratioCanvasNominalSystematicGraphFwd"));
  TGraphErrors *ratioCanvasNominalStatisticalGraph =
      static_cast<TGraphErrors *>(nominalStatisticalGraph->Clone("ratioCanvasNominalStatisticalGraphFwd"));
  TGraphErrors *ratioCanvasFitStatisticalGraph =
      static_cast<TGraphErrors *>(fitStatisticalGraph->Clone("ratioCanvasFitStatisticalGraphFwd"));

  TGraphErrors *yieldRatioGraph = new TGraphErrors(
      fitBinCount, fitPtCenters, yieldRatios, fitPtHalfWidths, yieldRatioStatErrors);
  yieldRatioGraph->SetName("yieldRatioGraphFwd");
  yieldRatioGraph->SetTitle("");
  yieldRatioGraph->SetMarkerStyle(20);
  yieldRatioGraph->SetMarkerSize(1.4);
  yieldRatioGraph->SetMarkerColor(kAzure + 2);
  yieldRatioGraph->SetLineColor(kAzure + 2);
  yieldRatioGraph->SetLineWidth(2);

  TCanvas *ratioCanvas = new TCanvas("ratioCanvasFwd", "OO raw yields and ratio", 800, 800);
  ratioCanvas->SetTitle("");

  TPad *ratioTopPad = new TPad("ratioTopPadFwd", "ratioTopPadFwd", 0.0, 0.25, 1.0, 1.0);
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
  ratioCanvasNominalSystematicGraph->GetXaxis()->SetLimits(nominalPtLow[0], nominalPtHigh[nominalBinCount - 1]);
  ratioCanvasNominalSystematicGraph->SetMinimum(0.8 * minimumYield);
  ratioCanvasNominalSystematicGraph->SetMaximum(4 * maximumYield);
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
    label.DrawLatex(0.19, 0.865, "1.6 < |y| < 2.4");
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
  TPad *ratioBottomPad = new TPad("ratioBottomPadFwd", "ratioBottomPadFwd", 0.0, 0.0, 1.0, 0.25);
  ratioBottomPad->SetTopMargin(0.00001);
  ratioBottomPad->SetBottomMargin(0.4);
  ratioBottomPad->Draw();
  ratioBottomPad->cd();

  yieldRatioGraph->Draw("AP");
  yieldRatioGraph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  yieldRatioGraph->GetXaxis()->CenterTitle();
  yieldRatioGraph->GetXaxis()->SetTitleSize(0.15);
  yieldRatioGraph->GetXaxis()->SetLabelSize(0.10);
  yieldRatioGraph->GetXaxis()->SetLimits(nominalPtLow[0], nominalPtHigh[nominalBinCount - 1]);
  yieldRatioGraph->GetYaxis()->SetTitle("2D/Nom.");
  yieldRatioGraph->GetYaxis()->CenterTitle();
  yieldRatioGraph->GetYaxis()->SetNdivisions(505);
  yieldRatioGraph->GetYaxis()->SetTitleSize(0.12);
  yieldRatioGraph->GetYaxis()->SetLabelSize(0.10);
  yieldRatioGraph->GetYaxis()->SetTitleOffset(0.55);
  yieldRatioGraph->SetMinimum(0.85);
  yieldRatioGraph->SetMaximum(1.15);

  TLine ratioReferenceLine(nominalPtLow[0], 1.0, nominalPtHigh[nominalBinCount - 1], 1.0);
  ratioReferenceLine.SetLineStyle(2);
  ratioReferenceLine.SetLineWidth(2);
  ratioReferenceLine.Draw("SAME");
  yieldRatioGraph->Draw("P SAME");

  ratioCanvas->SaveAs(outputPath("raw_yield_OO_fwd_ratio.pdf"));
}
