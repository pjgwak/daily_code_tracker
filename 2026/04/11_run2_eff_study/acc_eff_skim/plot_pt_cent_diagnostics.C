#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace
{
TH1D *ProjectPt(TH2D *hist, double centLow, double centHigh, const char *name)
{
  if (!hist)
    return nullptr;

  TAxis *yAxis = hist->GetYaxis();
  const int firstBin = yAxis->FindFixBin(centLow + 1e-6);
  const int lastBin = yAxis->FindFixBin(centHigh - 1e-6);
  TH1D *proj = hist->ProjectionX(name, firstBin, lastBin, "e");
  proj->SetDirectory(nullptr);
  return proj;
}

TH1D *BuildEffProjection(TH2D *num, TH2D *den, double centLow, double centHigh,
                         const char *name, const char *title)
{
  if (!num || !den)
    return nullptr;

  TH1D *projNum = ProjectPt(num, centLow, centHigh, Form("%s_num", name));
  TH1D *projDen = ProjectPt(den, centLow, centHigh, Form("%s_den", name));
  if (!projNum || !projDen)
    return nullptr;

  TH1D *eff = static_cast<TH1D *>(projNum->Clone(name));
  eff->SetDirectory(nullptr);
  eff->SetTitle(title);
  eff->Divide(projNum, projDen, 1.0, 1.0);

  delete projNum;
  delete projDen;
  return eff;
}

void StyleHist(TH1D *hist, int color, int markerStyle)
{
  if (!hist)
    return;

  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
}

void StyleHist2D(TH2D *hist)
{
  if (!hist)
    return;

  hist->SetTitleSize(0.05, "XY");
  hist->SetLabelSize(0.04, "XY");
  hist->SetTitleOffset(1.0, "X");
  hist->SetTitleOffset(1.15, "Y");
  hist->SetMarkerSize(1.35);
}

void ComputeVisibleRange(TH2D *hist, double &minValue, double &maxValue)
{
  minValue = std::numeric_limits<double>::max();
  maxValue = -std::numeric_limits<double>::max();

  if (!hist)
    return;

  for (int ix = 1; ix <= hist->GetNbinsX(); ++ix)
  {
    for (int iy = 1; iy <= hist->GetNbinsY(); ++iy)
    {
      const double value = hist->GetBinContent(ix, iy);
      if (!std::isfinite(value) || value == 0.0)
        continue;
      minValue = std::min(minValue, value);
      maxValue = std::max(maxValue, value);
    }
  }

  if (minValue == std::numeric_limits<double>::max())
  {
    minValue = 0.0;
    maxValue = 1.0;
  }
}

void SetRatioZRange(TH2D *ratio)
{
  double minValue = 0.0;
  double maxValue = 1.0;
  ComputeVisibleRange(ratio, minValue, maxValue);

  if (maxValue <= minValue)
  {
    minValue -= 0.05;
    maxValue += 0.05;
  }

  const double span = maxValue - minValue;
  const double padding = std::max(0.02, 0.12 * span);
  ratio->SetMinimum(minValue - padding);
  ratio->SetMaximum(maxValue + padding);
}

void WriteCountTable(TH2D *histPtW, TH2D *histNoPtW, const std::string &txtPath)
{
  if (!histPtW || !histNoPtW)
    return;

  std::ofstream out(txtPath);
  if (!out)
  {
    std::cout << "[ERROR] failed to write table: " << txtPath << "\n";
    return;
  }

  out << "# pt_low pt_high cent_low cent_high ptw_count noptw_count ratio_ptw_over_noptw\n";
  out << std::setprecision(10);

  for (int ix = 1; ix <= histPtW->GetNbinsX(); ++ix)
  {
    const double ptLow = histPtW->GetXaxis()->GetBinLowEdge(ix);
    const double ptHigh = histPtW->GetXaxis()->GetBinUpEdge(ix);

    for (int iy = 1; iy <= histPtW->GetNbinsY(); ++iy)
    {
      const double centLow = histPtW->GetYaxis()->GetBinLowEdge(iy);
      const double centHigh = histPtW->GetYaxis()->GetBinUpEdge(iy);
      const double valuePtW = histPtW->GetBinContent(ix, iy);
      const double valueNoPtW = histNoPtW->GetBinContent(ix, iy);
      const double ratio = (valueNoPtW != 0.0) ? (valuePtW / valueNoPtW) : 0.0;

      out << ptLow << ' ' << ptHigh << ' '
          << centLow << ' ' << centHigh << ' '
          << valuePtW << ' ' << valueNoPtW << ' ' << ratio << '\n';
    }
  }
}

void SaveCanvas(TCanvas &canvas, const std::string &basePath)
{
  canvas.SaveAs((basePath + ".png").c_str());
  canvas.SaveAs((basePath + ".pdf").c_str());
}

void Draw2DComparison(TFile *inputFile,
                      const char *histNamePtW,
                      const char *histNameNoPtW,
                      const std::string &basePath,
                      bool isEfficiency)
{
  TH2D *histPtW = dynamic_cast<TH2D *>(inputFile->Get(histNamePtW));
  TH2D *histNoPtW = dynamic_cast<TH2D *>(inputFile->Get(histNameNoPtW));
  if (!histPtW || !histNoPtW)
  {
    std::cout << "[ERROR] missing 2D histogram pair: "
              << histNamePtW << ", " << histNameNoPtW << "\n";
    return;
  }

  TH2D *ratio = static_cast<TH2D *>(histPtW->Clone(Form("%s_ratio", histNamePtW)));
  ratio->SetDirectory(nullptr);
  ratio->SetTitle("p_{T} reweight / no p_{T} reweight");
  ratio->Divide(histNoPtW);
  StyleHist2D(histPtW);
  StyleHist2D(histNoPtW);
  StyleHist2D(ratio);

  if (isEfficiency)
  {
    histPtW->SetMinimum(0.0);
    histNoPtW->SetMinimum(0.0);
  }
  SetRatioZRange(ratio);

  std::string safeTag = histNamePtW;
  std::replace(safeTag.begin(), safeTag.end(), '/', '_');
  TCanvas canvas(Form("c_%s", safeTag.c_str()), "", 2400, 850);
  canvas.Divide(3, 1);

  canvas.cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.13);
  gStyle->SetPaintTextFormat(isEfficiency ? ".3f" : ".3e");
  histPtW->Draw(isEfficiency ? "COLZ TEXT" : "COLZ TEXT");

  canvas.cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.13);
  gStyle->SetPaintTextFormat(isEfficiency ? ".3f" : ".3e");
  histNoPtW->Draw(isEfficiency ? "COLZ TEXT" : "COLZ TEXT");

  canvas.cd(3);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.13);
  gStyle->SetPaintTextFormat(".3f");
  ratio->Draw("COLZ TEXT");

  SaveCanvas(canvas, basePath);
  delete ratio;
}

void DrawProjectionComparison(TFile *inputFile,
                              const char *denNamePtW,
                              const char *denNameNoPtW,
                              const char *numNamePtW,
                              const char *numNameNoPtW,
                              const std::string &quantity,
                              const std::string &basePath)
{
  TH2D *denPtW = dynamic_cast<TH2D *>(inputFile->Get(denNamePtW));
  TH2D *denNoPtW = dynamic_cast<TH2D *>(inputFile->Get(denNameNoPtW));
  TH2D *numPtW = dynamic_cast<TH2D *>(inputFile->Get(numNamePtW));
  TH2D *numNoPtW = dynamic_cast<TH2D *>(inputFile->Get(numNameNoPtW));
  if (!denPtW || !denNoPtW || !numPtW || !numNoPtW)
  {
    std::cout << "[ERROR] missing projection inputs for " << basePath << "\n";
    return;
  }

  const double centRanges[4][2] = {{0.0, 180.0}, {0.0, 20.0}, {20.0, 60.0}, {60.0, 180.0}};
  const char *centLabels[4] = {"0-90% inclusive", "0-10%", "10-30%", "30-90%"};
  TH1D *histPtWKeep[4] = {nullptr, nullptr, nullptr, nullptr};
  TH1D *histNoPtWKeep[4] = {nullptr, nullptr, nullptr, nullptr};

  std::string safeTag = basePath;
  std::replace(safeTag.begin(), safeTag.end(), '/', '_');
  TCanvas canvas(Form("c_proj_%s", safeTag.c_str()), "", 1600, 1200);
  canvas.Divide(2, 2);

  for (int i = 0; i < 4; ++i)
  {
    TH1D *histPtW = nullptr;
    TH1D *histNoPtW = nullptr;
    if (quantity == "den")
    {
      histPtW = ProjectPt(denPtW, centRanges[i][0], centRanges[i][1],
                          Form("proj_den_ptw_%d_%s", i, safeTag.c_str()));
      histNoPtW = ProjectPt(denNoPtW, centRanges[i][0], centRanges[i][1],
                            Form("proj_den_noptw_%d_%s", i, safeTag.c_str()));
    }
    else if (quantity == "num")
    {
      histPtW = ProjectPt(numPtW, centRanges[i][0], centRanges[i][1],
                          Form("proj_num_ptw_%d_%s", i, safeTag.c_str()));
      histNoPtW = ProjectPt(numNoPtW, centRanges[i][0], centRanges[i][1],
                            Form("proj_num_noptw_%d_%s", i, safeTag.c_str()));
    }
    else
    {
      histPtW = BuildEffProjection(numPtW, denPtW, centRanges[i][0], centRanges[i][1],
                                   Form("proj_eff_ptw_%d_%s", i, safeTag.c_str()),
                                   ";p_{T} (GeV/c);Efficiency");
      histNoPtW = BuildEffProjection(numNoPtW, denNoPtW, centRanges[i][0], centRanges[i][1],
                                     Form("proj_eff_noptw_%d_%s", i, safeTag.c_str()),
                                     ";p_{T} (GeV/c);Efficiency");
    }

    if (!histPtW || !histNoPtW)
      continue;

    histPtWKeep[i] = histPtW;
    histNoPtWKeep[i] = histNoPtW;

    StyleHist(histPtW, kRed + 1, 20);
    StyleHist(histNoPtW, kBlue + 1, 24);

    const double maxY = std::max(histPtW->GetMaximum(), histNoPtW->GetMaximum());
    histPtW->SetMaximum((quantity == "eff") ? std::min(1.05, std::max(0.2, maxY * 1.25)) : maxY * 1.25);
    histPtW->SetMinimum(0.0);
    histNoPtW->SetMinimum(0.0);
    histPtW->SetTitle(Form(";p_{T} (GeV/c);%s", quantity == "eff" ? "Efficiency" : "Weighted counts"));

    canvas.cd(i + 1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    histPtW->Draw("E1");
    histNoPtW->Draw("E1 SAME");

    TLegend legend(0.52, 0.72, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(histPtW, "p_{T} reweight on", "lep");
    legend.AddEntry(histNoPtW, "p_{T} reweight off", "lep");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.045);
    label.DrawLatex(0.15, 0.92, centLabels[i]);
  }

  SaveCanvas(canvas, basePath);

  for (int i = 0; i < 4; ++i)
  {
    delete histPtWKeep[i];
    delete histNoPtWKeep[i];
  }
}
} // namespace

void plot_pt_cent_diagnostics(
    const char *inputPath = "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/skim_roots/eff_PbPb2018_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
    const char *outputDir = "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/plots/pt_cent_diagnostics")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".3g");

  if (gSystem->AccessPathName(outputDir) && gSystem->mkdir(outputDir, true) != 0)
  {
    std::cout << "[ERROR] failed to create output directory: " << outputDir << "\n";
    return;
  }

  TFile *inputFile = TFile::Open(inputPath, "READ");
  if (!inputFile || inputFile->IsZombie())
  {
    std::cout << "[ERROR] failed to open input ROOT file: " << inputPath << "\n";
    return;
  }

  Draw2DComparison(inputFile, "hist_eff2d_den_mid_ptw", "hist_eff2d_den_mid_noptw",
                   std::string(outputDir) + "/mid_den_2d", false);
  Draw2DComparison(inputFile, "hist_eff2d_num_mid_ptw", "hist_eff2d_num_mid_noptw",
                   std::string(outputDir) + "/mid_num_2d", false);
  Draw2DComparison(inputFile, "hist_eff2d_mid_ptw", "hist_eff2d_mid_noptw",
                   std::string(outputDir) + "/mid_eff_2d", true);
  Draw2DComparison(inputFile, "hist_eff2d_den_fwd_ptw", "hist_eff2d_den_fwd_noptw",
                   std::string(outputDir) + "/fwd_den_2d", false);
  Draw2DComparison(inputFile, "hist_eff2d_num_fwd_ptw", "hist_eff2d_num_fwd_noptw",
                   std::string(outputDir) + "/fwd_num_2d", false);
  Draw2DComparison(inputFile, "hist_eff2d_fwd_ptw", "hist_eff2d_fwd_noptw",
                   std::string(outputDir) + "/fwd_eff_2d", true);

  WriteCountTable(dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_den_mid_ptw")),
                  dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_den_mid_noptw")),
                  std::string(outputDir) + "/mid_den_2d.txt");
  WriteCountTable(dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_num_mid_ptw")),
                  dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_num_mid_noptw")),
                  std::string(outputDir) + "/mid_num_2d.txt");
  WriteCountTable(dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_den_fwd_ptw")),
                  dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_den_fwd_noptw")),
                  std::string(outputDir) + "/fwd_den_2d.txt");
  WriteCountTable(dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_num_fwd_ptw")),
                  dynamic_cast<TH2D *>(inputFile->Get("hist_eff2d_num_fwd_noptw")),
                  std::string(outputDir) + "/fwd_num_2d.txt");

  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_mid_ptw", "hist_eff2d_den_mid_noptw",
                           "hist_eff2d_num_mid_ptw", "hist_eff2d_num_mid_noptw",
                           "den", std::string(outputDir) + "/mid_den_proj");
  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_mid_ptw", "hist_eff2d_den_mid_noptw",
                           "hist_eff2d_num_mid_ptw", "hist_eff2d_num_mid_noptw",
                           "num", std::string(outputDir) + "/mid_num_proj");
  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_mid_ptw", "hist_eff2d_den_mid_noptw",
                           "hist_eff2d_num_mid_ptw", "hist_eff2d_num_mid_noptw",
                           "eff", std::string(outputDir) + "/mid_eff_proj");
  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_fwd_ptw", "hist_eff2d_den_fwd_noptw",
                           "hist_eff2d_num_fwd_ptw", "hist_eff2d_num_fwd_noptw",
                           "den", std::string(outputDir) + "/fwd_den_proj");
  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_fwd_ptw", "hist_eff2d_den_fwd_noptw",
                           "hist_eff2d_num_fwd_ptw", "hist_eff2d_num_fwd_noptw",
                           "num", std::string(outputDir) + "/fwd_num_proj");
  DrawProjectionComparison(inputFile,
                           "hist_eff2d_den_fwd_ptw", "hist_eff2d_den_fwd_noptw",
                           "hist_eff2d_num_fwd_ptw", "hist_eff2d_num_fwd_noptw",
                           "eff", std::string(outputDir) + "/fwd_eff_proj");

  inputFile->Close();
  delete inputFile;

  std::cout << "[INFO] plots saved under " << outputDir << "\n";
}
