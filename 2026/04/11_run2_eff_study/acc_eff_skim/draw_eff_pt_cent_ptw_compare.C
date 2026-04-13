#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

using std::cout;
using std::endl;

namespace
{
TString resolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}

void styleHist(TH1D *hist, Color_t color, Style_t markerStyle, Style_t lineStyle)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(1.5);
  hist->SetLineStyle(lineStyle);
  hist->SetLineWidth(3);
}

void configureAxes(TH1D *hist, const char *xTitle)
{
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle("Efficiency");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.25);
}

void drawCmsLabel(const TString &sampleLabel, const TString &regionLabel)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);

  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.935, "PbPb 2018 5.02 TeV");

  tx.SetTextSize(0.040);
  tx.SetTextFont(72);
  tx.SetTextAlign(11);
  tx.DrawLatex(0.18, 0.935, "CMS Internal");

  tx.SetTextFont(42);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.18, 0.865, sampleLabel);
  tx.DrawLatex(0.18, 0.815, regionLabel);
}

TH1D *loadHist(TFile *file, const char *histName, const char *cloneName, bool required = true)
{
  TH1D *hist = dynamic_cast<TH1D *>(file->Get(histName));
  if (!hist)
  {
    if (required)
      cout << "[ERROR] missing histogram: " << histName << endl;
    else
      cout << "[WARN] optional histogram not found: " << histName << endl;
    return nullptr;
  }

  TH1D *cloned = static_cast<TH1D *>(hist->Clone(cloneName));
  cloned->SetDirectory(nullptr);
  return cloned;
}

TH1D *loadHistWithFallback(TFile *file, const char *preferredName,
                           const char *fallbackName, const char *cloneName, bool required = true)
{
  TH1D *hist = dynamic_cast<TH1D *>(file->Get(preferredName));
  if (!hist && fallbackName)
    hist = dynamic_cast<TH1D *>(file->Get(fallbackName));
  if (!hist)
  {
    if (required)
      cout << "[ERROR] missing histogram: " << preferredName << endl;
    else
      cout << "[WARN] missing histogram: " << preferredName << " (fallback: " << fallbackName << ")" << endl;
    return nullptr;
  }

  TH1D *cloned = static_cast<TH1D *>(hist->Clone(cloneName));
  cloned->SetDirectory(nullptr);
  return cloned;
}

TH1D *makeFractionHist(TH1D *numerator, TH1D *denominator, const char *cloneName)
{
  TH1D *hist = static_cast<TH1D *>(numerator->Clone(cloneName));
  hist->SetDirectory(nullptr);
  hist->Reset("ICES");
  for (int ibin = 1; ibin <= numerator->GetNbinsX(); ++ibin)
  {
    const double den = denominator->GetBinContent(ibin);
    if (den <= 0.0)
      continue;
    hist->SetBinContent(ibin, numerator->GetBinContent(ibin) / den);
  }
  return hist;
}

void savePtWComparePlot(TH1D *histOff, TH1D *histOn, const TString &outPath,
                        const TString &sampleLabel, const TString &regionLabel,
                        const char *xTitle)
{
  TCanvas canv("c_ptw_compare", "", 800, 800);
  canv.cd();

  TPad pad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad.SetTopMargin(0.08);
  pad.SetBottomMargin(0.13);
  pad.SetLeftMargin(0.13);
  pad.SetRightMargin(0.04);
  pad.Draw();
  pad.cd();

  gPad->SetTicks(1, 1);
  gStyle->SetOptStat(0);

  std::unique_ptr<TH1D> frame(static_cast<TH1D *>(histOff->Clone(Form("%s_frame", histOff->GetName()))));
  frame->Reset("ICES");
  configureAxes(frame.get(), xTitle);
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, std::max(histOff->GetMaximum(), histOn->GetMaximum()) * 1.35));
  frame->Draw("AXIS");

  histOff->Draw("E1 SAME");
  histOn->Draw("E1 SAME");

  drawCmsLabel(sampleLabel, regionLabel);

  TLegend leg(0.56, 0.78, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.AddEntry(histOff, "p_{T} weight OFF", "lp");
  leg.AddEntry(histOn, "p_{T} weight ON", "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveFwdCentDebugPlot(TH1D *effFwdOff, TH1D *effFwdOn,
                          TH1D *effHighOff, TH1D *effHighOn,
                          TH1D *effLowOff, TH1D *effLowOn,
                          TH1D *lowFracOff, TH1D *lowFracOn,
                          const TString &outPath, const TString &sampleLabel)
{
  TCanvas canv("c_cent_fwd_debug", "", 850, 900);
  canv.cd();

  TPad top("top", "top", 0.0, 0.36, 1.0, 1.0);
  top.SetTopMargin(0.08);
  top.SetBottomMargin(0.02);
  top.SetLeftMargin(0.13);
  top.SetRightMargin(0.04);
  top.Draw();

  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.36);
  bot.SetTopMargin(0.02);
  bot.SetBottomMargin(0.28);
  bot.SetLeftMargin(0.13);
  bot.SetRightMargin(0.04);
  bot.Draw();

  top.cd();
  gPad->SetTicks(1, 1);
  std::unique_ptr<TH1D> frameTop(static_cast<TH1D *>(effFwdOff->Clone("frameTop")));
  frameTop->Reset("ICES");
  configureAxes(frameTop.get(), "Centrality");
  frameTop->GetXaxis()->SetLabelSize(0.0);
  frameTop->GetXaxis()->SetTitleSize(0.0);
  frameTop->SetMinimum(0.0);
  frameTop->SetMaximum(std::max(1.05, std::max(effFwdOn->GetMaximum(), effHighOn->GetMaximum()) * 1.35));
  frameTop->Draw("AXIS");

  styleHist(effFwdOff, kBlack, 24, 2);
  styleHist(effFwdOn, kBlack, 20, 1);
  styleHist(effHighOff, kAzure + 2, 24, 2);
  styleHist(effHighOn, kAzure + 2, 20, 1);
  styleHist(effLowOff, kOrange + 7, 25, 2);
  styleHist(effLowOn, kOrange + 7, 21, 1);

  effFwdOff->Draw("E1 SAME");
  effFwdOn->Draw("E1 SAME");
  effHighOff->Draw("E1 SAME");
  effHighOn->Draw("E1 SAME");
  effLowOff->Draw("E1 SAME");
  effLowOn->Draw("E1 SAME");

  drawCmsLabel(sampleLabel, "1.6 < |y| < 2.4 centrality breakdown");

  TLegend legTop(0.47, 0.63, 0.92, 0.90);
  legTop.SetBorderSize(0);
  legTop.SetFillStyle(0);
  legTop.SetTextSize(0.026);
  legTop.AddEntry(effFwdOff, "combined 3.5 < p_{T} < 40, ptW OFF", "lp");
  legTop.AddEntry(effFwdOn, "combined 3.5 < p_{T} < 40, ptW ON", "lp");
  legTop.AddEntry(effHighOff, "high p_{T} 6.5 < p_{T} < 40, ptW OFF", "lp");
  legTop.AddEntry(effHighOn, "high p_{T} 6.5 < p_{T} < 40, ptW ON", "lp");
  legTop.AddEntry(effLowOff, "low p_{T} 3.5 < p_{T} < 6.5, ptW OFF", "lp");
  legTop.AddEntry(effLowOn, "low p_{T} 3.5 < p_{T} < 6.5, ptW ON", "lp");
  legTop.Draw();

  bot.cd();
  gPad->SetTicks(1, 1);
  std::unique_ptr<TH1D> frameBot(static_cast<TH1D *>(lowFracOff->Clone("frameBot")));
  frameBot->Reset("ICES");
  configureAxes(frameBot.get(), "Centrality");
  frameBot->GetYaxis()->SetTitle("Low-p_{T} fraction");
  frameBot->GetYaxis()->SetTitleOffset(1.25);
  frameBot->SetMinimum(0.0);
  frameBot->SetMaximum(1.0);
  frameBot->Draw("AXIS");

  styleHist(lowFracOff, kGreen + 2, 24, 2);
  styleHist(lowFracOn, kMagenta + 1, 20, 1);
  lowFracOff->Draw("E1 SAME");
  lowFracOn->Draw("E1 SAME");

  TLine line;
  line.SetLineColor(kGray + 1);
  line.SetLineStyle(2);
  line.DrawLine(frameBot->GetXaxis()->GetXmin(), 0.5, frameBot->GetXaxis()->GetXmax(), 0.5);

  TLegend legBot(0.56, 0.72, 0.92, 0.90);
  legBot.SetBorderSize(0);
  legBot.SetFillStyle(0);
  legBot.SetTextSize(0.028);
  legBot.AddEntry(lowFracOff, "low-p_{T} denominator fraction, ptW OFF", "lp");
  legBot.AddEntry(lowFracOn, "low-p_{T} denominator fraction, ptW ON", "lp");
  legBot.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

void draw_eff_pt_cent_ptw_compare(bool isNCollW = true, bool isGenW = true, bool isTnPW = true)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString outDir = baseDir + "/outputs/figs/eff_pt_cent_ptw_compare";
  const TString ptDir = outDir + "/pt";
  const TString centDir = outDir + "/cent";

  gSystem->mkdir(outDir, true);
  gSystem->mkdir(ptDir, true);
  gSystem->mkdir(centDir, true);

  const TString commonTagOff = Form("_ncollW%d_genW%d_ptW0_tnpW%d",
                                    isNCollW ? 1 : 0, isGenW ? 1 : 0, isTnPW ? 1 : 0);
  const TString commonTagOn = Form("_ncollW%d_genW%d_ptW1_tnpW%d",
                                   isNCollW ? 1 : 0, isGenW ? 1 : 0, isTnPW ? 1 : 0);

  struct SampleInput {
    TString shortLabel;
    TString displayLabel;
  };

  const std::vector<SampleInput> samples = {
      {"PR", "Prompt MC"},
      {"NP", "Nonprompt MC"},
  };

  struct HistCase {
    const char *histName;
    const char *suffix;
    const char *regionLabel;
    const char *xTitle;
    bool isCent;
  };

  const std::vector<HistCase> histCases = {
      {"hist_eff_mid", "mid", "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c", "p_{T} (GeV/c)", false},
      {"hist_eff_fwd", "fwd", "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent0010", "pt_mid_cent0010", "0-10%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent1020", "pt_mid_cent1020", "10-20%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent2030", "pt_mid_cent2030", "20-30%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent3040", "pt_mid_cent3040", "30-40%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent4050", "pt_mid_cent4050", "40-50%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_mid_cent5090", "pt_mid_cent5090", "50-90%, |y| < 1.6", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_fwd_cent0010", "pt_fwd_cent0010", "0-10%, 1.6 < |y| < 2.4", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_fwd_cent1030", "pt_fwd_cent1030", "10-30%, 1.6 < |y| < 2.4", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_fwd_cent3050", "pt_fwd_cent3050", "30-50%, 1.6 < |y| < 2.4", "p_{T} (GeV/c)", false},
      {"hist_eff_pt_fwd_cent5090", "pt_fwd_cent5090", "50-90%, 1.6 < |y| < 2.4", "p_{T} (GeV/c)", false},
      {"hist_eff_cent_mid", "cent_mid", "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c", "Centrality", true},
      {"hist_eff_cent_fwd", "cent_fwd", "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c", "Centrality", true},
  };

  for (const auto &sample : samples)
  {
    const TString pathOff = inputDir + "/eff_PbPb2018_isMC1_" + sample.shortLabel + commonTagOff + "_Dimuon_MiniAOD.root";
    const TString pathOn = inputDir + "/eff_PbPb2018_isMC1_" + sample.shortLabel + commonTagOn + "_Dimuon_MiniAOD.root";

    std::unique_ptr<TFile> fOff(TFile::Open(pathOff, "READ"));
    std::unique_ptr<TFile> fOn(TFile::Open(pathOn, "READ"));
    if (!fOff || fOff->IsZombie())
    {
      cout << "[ERROR] failed to open pT-off file: " << pathOff << endl;
      continue;
    }
    if (!fOn || fOn->IsZombie())
    {
      cout << "[ERROR] failed to open pT-on file: " << pathOn << endl;
      continue;
    }

    for (const auto &histCase : histCases)
    {
      std::unique_ptr<TH1D> hOff(loadHistWithFallback(fOff.get(), histCase.histName, nullptr,
                                                      TString::Format("%s_%s_off", sample.shortLabel.Data(), histCase.suffix),
                                                      false));
      std::unique_ptr<TH1D> hOn(loadHistWithFallback(fOn.get(), histCase.histName, nullptr,
                                                     TString::Format("%s_%s_on", sample.shortLabel.Data(), histCase.suffix),
                                                     false));
      if (!hOff || !hOn)
        continue;

      styleHist(hOff.get(), kGray + 2, 24, 2);
      styleHist(hOn.get(), kRed + 1, 20, 1);

      TString sampleTag = sample.shortLabel;
      sampleTag.ToLower();
      const TString subDir = histCase.isCent ? centDir : ptDir;
      const TString outPath = subDir + "/eff_" + sampleTag + "_ptw_compare_" + histCase.suffix +
                              Form("_ncollW%d_genW%d_tnpW%d",
                                   isNCollW ? 1 : 0, isGenW ? 1 : 0, isTnPW ? 1 : 0);
      savePtWComparePlot(hOff.get(), hOn.get(), outPath, sample.displayLabel, histCase.regionLabel, histCase.xTitle);
    }
  }
}
