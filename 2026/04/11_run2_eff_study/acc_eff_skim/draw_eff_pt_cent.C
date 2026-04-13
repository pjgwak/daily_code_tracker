#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
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

void styleHist(TH1D *hist, Color_t color, Style_t markerStyle, double markerSize)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(markerSize);
  hist->SetLineWidth(3);
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

void configureAxes(TH1D *hist, const char *xTitle, const char *yTitle)
{
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.25);
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
  // Detach from the TFile so canvases remain valid after the file is closed.
  return cloned;
}

TH1D *loadHistWithFallback(TFile *file, const char *preferredName,
                           const char *fallbackName, const char *cloneName)
{
  TH1D *hist = dynamic_cast<TH1D *>(file->Get(preferredName));
  if (!hist && fallbackName)
    hist = dynamic_cast<TH1D *>(file->Get(fallbackName));
  if (!hist)
  {
    cout << "[WARN] missing histogram: " << preferredName;
    if (fallbackName)
      cout << " (fallback: " << fallbackName << ")";
    cout << endl;
    return nullptr;
  }

  TH1D *cloned = static_cast<TH1D *>(hist->Clone(cloneName));
  cloned->SetDirectory(nullptr);
  return cloned;
}

void saveSinglePlot(TH1D *hist, const TString &outPath, const TString &sampleLabel,
                    const TString &regionLabel, const char *xTitle, const char *yTitle)
{
  TCanvas canv("c_eff_single", "", 800, 800);
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

  configureAxes(hist, xTitle, yTitle);
  hist->SetMinimum(0.0);
  hist->SetMaximum(std::max(1.05, hist->GetMaximum() * 1.35));
  // E1 visualizes the statistical errors already stored in the input histogram.
  hist->Draw("E1");

  drawCmsLabel(sampleLabel, regionLabel);

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveOverlayPlot(TH1D *histPr, TH1D *histNp, const TString &outPath, const TString &regionLabel,
                     const char *xTitle)
{
  TCanvas canv("c_eff_overlay", "", 800, 800);
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

  std::unique_ptr<TH1D> frame(static_cast<TH1D *>(histPr->Clone(Form("%s_frame", histPr->GetName()))));
  frame->Reset("ICES");
  configureAxes(frame.get(), xTitle, "Efficiency");
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, std::max(histPr->GetMaximum(), histNp->GetMaximum()) * 1.35));
  frame->Draw("AXIS");

  // The drawing macro does not recalculate uncertainties; it only displays them.
  histPr->Draw("E1 SAME");
  histNp->Draw("E1 SAME");

  drawCmsLabel("Prompt / Nonprompt MC", regionLabel);

  TLegend leg(0.58, 0.78, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.AddEntry(histPr, "Prompt", "lp");
  leg.AddEntry(histNp, "Nonprompt", "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveSummaryPlot(const std::vector<std::pair<TH1D *, TH1D *>> &histPairs,
                     const std::vector<TString> &regionLabels,
                     const TString &outPath, const char *xTitle,
                     const TString &headerLabel,
                     bool compactLegend = false,
                     double fixedYMax = -1.0)
{
  if (histPairs.empty() || histPairs.size() != regionLabels.size())
    return;

  TCanvas canv("c_eff_summary", "", 800, 800);
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

  double maxY = 0.0;
  for (const auto &pair : histPairs)
  {
    if (pair.first)
      maxY = std::max(maxY, pair.first->GetMaximum());
    if (pair.second)
      maxY = std::max(maxY, pair.second->GetMaximum());
  }

  std::unique_ptr<TH1D> frame(static_cast<TH1D *>(histPairs.front().first->Clone("summary_frame")));
  frame->Reset("ICES");
  configureAxes(frame.get(), xTitle, "Efficiency");
  frame->SetMinimum(0.0);
  frame->SetMaximum((fixedYMax > 0.0) ? fixedYMax : std::max(1.05, maxY * 1.35));
  frame->Draw("AXIS");

  const Color_t caseColors[] = {kBlack, kRed + 1, kAzure + 1, kGreen + 2, kMagenta + 1,
                                kOrange + 7, kCyan + 2, kPink + 7, kGray + 2, kYellow + 2};
  const Style_t prMarkers[] = {20, 21, 22, 23, 29, 33, 34, 47, 43, 41};
  const Style_t npMarkers[] = {24, 25, 26, 32, 30, 27, 28, 46, 42, 40};

  const size_t nLabels = regionLabels.size();
  const int nColumns = (nLabels > 6 || compactLegend) ? 2 : 1;
  const double textSize = compactLegend ? 0.0165 : ((nLabels > 8) ? 0.020 : ((nLabels > 5) ? 0.022 : 0.025));
  const double x1 = compactLegend ? 0.60 : ((nColumns == 2) ? 0.34 : 0.56);
  const double x2 = compactLegend ? 0.965 : 0.94;
  const double y1 = compactLegend ? 0.56 : ((nLabels > 8) ? 0.48 : ((nLabels > 5) ? 0.54 : 0.66));
  const double y2 = 0.90;
  TLegend leg(x1, y1, x2, y2);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(textSize);
  leg.SetNColumns(nColumns);
  if (compactLegend) {
    leg.SetColumnSeparation(-0.03);
    leg.SetMargin(0.12);
    leg.SetEntrySeparation(0.04);
  }

  for (size_t i = 0; i < histPairs.size(); ++i)
  {
    TH1D *histPr = histPairs[i].first;
    TH1D *histNp = histPairs[i].second;
    if (!histPr || !histNp)
      continue;

    styleHist(histPr, caseColors[i % 10], prMarkers[i % 10], 1.45);
    styleHist(histNp, caseColors[i % 10], npMarkers[i % 10], 1.45);
    histPr->SetLineStyle(1);
    histNp->SetLineStyle(7);
    histPr->SetLineWidth(3);
    histNp->SetLineWidth(3);

    histPr->Draw("E1 SAME");
    histNp->Draw("E1 SAME");

    leg.AddEntry(histPr, "PR " + regionLabels[i], "lp");
    leg.AddEntry(histNp, "NP " + regionLabels[i], "lp");
  }

  drawCmsLabel("Prompt / Nonprompt MC", headerLabel);

  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

void draw_eff_pt_cent(bool isNCollW = true, bool isGenW = true, bool isPtW = true, bool isTnPW = true)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString outDir = baseDir + "/outputs/figs/eff_pt_cent";
  const TString ptDir = outDir + "/pt";
  const TString ptPrDir = ptDir + "/pr";
  const TString ptNpDir = ptDir + "/np";
  const TString ptCompareDir = ptDir + "/compare";
  const TString centDir = outDir + "/cent";
  const TString centPrDir = centDir + "/pr";
  const TString centNpDir = centDir + "/np";
  const TString centCompareDir = centDir + "/compare";
  const TString summaryDir = outDir + "/summary";
  const TString summaryPtDir = summaryDir + "/pt";
  const TString summaryCentDir = summaryDir + "/cent";
  const TString weightTag = Form("_ncollW%d_genW%d_ptW%d_tnpW%d",
                                 isNCollW ? 1 : 0, isGenW ? 1 : 0,
                                 isPtW ? 1 : 0, isTnPW ? 1 : 0);

  gSystem->mkdir(outDir, true);
  gSystem->mkdir(ptDir, true);
  gSystem->mkdir(ptPrDir, true);
  gSystem->mkdir(ptNpDir, true);
  gSystem->mkdir(ptCompareDir, true);
  gSystem->mkdir(centDir, true);
  gSystem->mkdir(centPrDir, true);
  gSystem->mkdir(centNpDir, true);
  gSystem->mkdir(centCompareDir, true);
  gSystem->mkdir(summaryDir, true);
  gSystem->mkdir(summaryPtDir, true);
  gSystem->mkdir(summaryCentDir, true);

  const TString prPath =
      inputDir + "/eff_PbPb2018_isMC1_PR" + weightTag + "_Dimuon_MiniAOD.root";
  const TString npPath =
      inputDir + "/eff_PbPb2018_isMC1_NP" + weightTag + "_Dimuon_MiniAOD.root";

  std::unique_ptr<TFile> fPr(TFile::Open(prPath, "READ"));
  std::unique_ptr<TFile> fNp(TFile::Open(npPath, "READ"));
  if (!fPr || fPr->IsZombie())
  {
    cout << "[ERROR] failed to open prompt file: " << prPath << endl;
    return;
  }
  if (!fNp || fNp->IsZombie())
  {
    cout << "[ERROR] failed to open nonprompt file: " << npPath << endl;
    return;
  }

  std::unique_ptr<TH1D> hPrPtMid(loadHist(fPr.get(), "hist_eff_mid", "hPrPtMid"));
  std::unique_ptr<TH1D> hPrPtFwd(loadHist(fPr.get(), "hist_eff_fwd", "hPrPtFwd"));
  std::unique_ptr<TH1D> hNpPtMid(loadHist(fNp.get(), "hist_eff_mid", "hNpPtMid"));
  std::unique_ptr<TH1D> hNpPtFwd(loadHist(fNp.get(), "hist_eff_fwd", "hNpPtFwd"));
  if (!hPrPtMid || !hPrPtFwd || !hNpPtMid || !hNpPtFwd)
    return;

  std::unique_ptr<TH1D> hPrCentMid(loadHist(fPr.get(), "hist_eff_cent_mid", "hPrCentMid", false));
  std::unique_ptr<TH1D> hPrCentFwd(loadHist(fPr.get(), "hist_eff_cent_fwd", "hPrCentFwd", false));
  std::unique_ptr<TH1D> hNpCentMid(loadHist(fNp.get(), "hist_eff_cent_mid", "hNpCentMid", false));
  std::unique_ptr<TH1D> hNpCentFwd(loadHist(fNp.get(), "hist_eff_cent_fwd", "hNpCentFwd", false));

  styleHist(hPrPtMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrPtFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpPtMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpPtFwd.get(), kOrange + 7, 25, 1.7);

  const TString midLabel = "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c";
  const TString fwdLabel = "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c";

  saveSinglePlot(hPrPtMid.get(), ptPrDir + "/eff_pt_pr_mid" + weightTag, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hPrPtFwd.get(), ptPrDir + "/eff_pt_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hNpPtMid.get(), ptNpDir + "/eff_pt_np_mid" + weightTag, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hNpPtFwd.get(), ptNpDir + "/eff_pt_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency");

  saveOverlayPlot(hPrPtMid.get(), hNpPtMid.get(), ptCompareDir + "/eff_pt_pr_np_mid" + weightTag, midLabel, "p_{T} (GeV/c)");
  saveOverlayPlot(hPrPtFwd.get(), hNpPtFwd.get(), ptCompareDir + "/eff_pt_pr_np_fwd" + weightTag, fwdLabel, "p_{T} (GeV/c)");

  std::vector<std::unique_ptr<TH1D>> ptSummaryStore;
  std::vector<std::pair<TH1D *, TH1D *>> ptSummaryPairs;
  std::vector<TString> ptSummaryLabels;

  ptSummaryStore.emplace_back(static_cast<TH1D *>(hPrPtMid->Clone("hPrPtMidSummary")));
  ptSummaryStore.back()->SetDirectory(nullptr);
  ptSummaryStore.emplace_back(static_cast<TH1D *>(hNpPtMid->Clone("hNpPtMidSummary")));
  ptSummaryStore.back()->SetDirectory(nullptr);
  ptSummaryPairs.emplace_back(ptSummaryStore[ptSummaryStore.size() - 2].get(),
                              ptSummaryStore[ptSummaryStore.size() - 1].get());
  ptSummaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");

  ptSummaryStore.emplace_back(static_cast<TH1D *>(hPrPtFwd->Clone("hPrPtFwdSummary")));
  ptSummaryStore.back()->SetDirectory(nullptr);
  ptSummaryStore.emplace_back(static_cast<TH1D *>(hNpPtFwd->Clone("hNpPtFwdSummary")));
  ptSummaryStore.back()->SetDirectory(nullptr);
  ptSummaryPairs.emplace_back(ptSummaryStore[ptSummaryStore.size() - 2].get(),
                              ptSummaryStore[ptSummaryStore.size() - 1].get());
  ptSummaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");

  saveSummaryPlot(ptSummaryPairs, ptSummaryLabels,
                  summaryPtDir + "/eff_pt_pr_np_summary" + weightTag,
                  "p_{T} (GeV/c)", "p_{T} (GeV/c)");

  struct PtCentCase {
    const char *suffix;
    const char *label;
    const char *summaryLabel;
  };
  const std::vector<PtCentCase> ptCentCases = {
      {"pt_mid_cent0010", "0-10%, |y| < 1.6", "mid 0-10%"},
      {"pt_mid_cent1020", "10-20%, |y| < 1.6", "mid 10-20%"},
      {"pt_mid_cent2030", "20-30%, |y| < 1.6", "mid 20-30%"},
      {"pt_mid_cent3040", "30-40%, |y| < 1.6", "mid 30-40%"},
      {"pt_mid_cent4050", "40-50%, |y| < 1.6", "mid 40-50%"},
      {"pt_mid_cent5090", "50-90%, |y| < 1.6", "mid 50-90%"},
      {"pt_fwd_cent0010", "0-10%, 1.6 < |y| < 2.4", "fwd 0-10%"},
      {"pt_fwd_cent1030", "10-30%, 1.6 < |y| < 2.4", "fwd 10-30%"},
      {"pt_fwd_cent3050", "30-50%, 1.6 < |y| < 2.4", "fwd 30-50%"},
      {"pt_fwd_cent5090", "50-90%, 1.6 < |y| < 2.4", "fwd 50-90%"},
  };
  std::vector<std::pair<TH1D *, TH1D *>> ptCentSummaryPairs;
  std::vector<TString> ptCentSummaryLabels;
  for (const auto &centCase : ptCentCases)
  {
    const TString histName = TString::Format("hist_eff_%s", centCase.suffix);
    std::unique_ptr<TH1D> hPr(loadHist(fPr.get(), histName, TString::Format("%s_pr", histName.Data()), false));
    std::unique_ptr<TH1D> hNp(loadHist(fNp.get(), histName, TString::Format("%s_np", histName.Data()), false));
    if (!hPr || !hNp)
      continue;

    styleHist(hPr.get(), kAzure + 2, 24, 1.7);
    styleHist(hNp.get(), kOrange + 7, 25, 1.7);

    const TString suffix = centCase.suffix;
    const TString label = centCase.label;
    saveSinglePlot(hPr.get(), ptPrDir + "/eff_pt_pr_" + suffix + weightTag, "Prompt MC", label, "p_{T} (GeV/c)", "Efficiency");
    saveSinglePlot(hNp.get(), ptNpDir + "/eff_pt_np_" + suffix + weightTag, "Nonprompt MC", label, "p_{T} (GeV/c)", "Efficiency");
    saveOverlayPlot(hPr.get(), hNp.get(), ptCompareDir + "/eff_pt_pr_np_" + suffix + weightTag, label, "p_{T} (GeV/c)");
    ptCentSummaryPairs.emplace_back(hPr.release(), hNp.release());
    ptCentSummaryLabels.push_back(centCase.summaryLabel);
  }
  if (!ptCentSummaryPairs.empty())
    saveSummaryPlot(ptCentSummaryPairs, ptCentSummaryLabels,
                    summaryPtDir + "/eff_pt_pr_np_cent_summary" + weightTag,
                    "p_{T} (GeV/c)", "p_{T} (GeV/c)");

  auto saveRegionSummary = [&](const std::vector<const char *> &suffixes, const TString &outName, const TString &header) {
    std::vector<std::pair<TH1D *, TH1D *>> pairs;
    std::vector<TString> labels;
    std::vector<std::unique_ptr<TH1D>> store;
    for (const char *suffix : suffixes)
    {
      const TString histName = TString::Format("hist_eff_%s", suffix);
      std::unique_ptr<TH1D> hPr(loadHist(fPr.get(), histName, TString::Format("%s_pr", histName.Data()), false));
      std::unique_ptr<TH1D> hNp(loadHist(fNp.get(), histName, TString::Format("%s_np", histName.Data()), false));
      if (!hPr || !hNp)
        continue;
      styleHist(hPr.get(), kAzure + 2, 24, 1.7);
      styleHist(hNp.get(), kOrange + 7, 25, 1.7);
      labels.push_back(suffix);
      store.push_back(std::move(hPr));
      store.push_back(std::move(hNp));
      pairs.emplace_back(store[store.size() - 2].get(), store[store.size() - 1].get());
    }
    if (!pairs.empty())
      saveSummaryPlot(pairs, labels, summaryPtDir + outName + weightTag, "p_{T} (GeV/c)", header);
  };
  saveRegionSummary({"pt_mid_cent0010", "pt_mid_cent1020", "pt_mid_cent2030", "pt_mid_cent3040", "pt_mid_cent4050", "pt_mid_cent5090"},
                    "/eff_pt_pr_np_cent_summary_mid", "|y| < 1.6");
  saveRegionSummary({"pt_fwd_cent0010", "pt_fwd_cent1030", "pt_fwd_cent3050", "pt_fwd_cent5090"},
                    "/eff_pt_pr_np_cent_summary_fwd", "1.6 < |y| < 2.4");

  std::vector<std::pair<TH1D *, TH1D *>> centSummaryPairs;
  std::vector<TString> centSummaryLabels;
  if (hPrCentMid && hPrCentFwd && hNpCentMid && hNpCentFwd)
  {
    styleHist(hPrCentMid.get(), kAzure + 2, 24, 1.7);
    styleHist(hPrCentFwd.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpCentMid.get(), kOrange + 7, 25, 1.7);
    styleHist(hNpCentFwd.get(), kOrange + 7, 25, 1.7);

    saveSinglePlot(hPrCentMid.get(), centPrDir + "/eff_cent_pr_mid" + weightTag, "Prompt MC", midLabel, "Centrality", "Efficiency");
    saveSinglePlot(hPrCentFwd.get(), centPrDir + "/eff_cent_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "Centrality", "Efficiency");
    saveSinglePlot(hNpCentMid.get(), centNpDir + "/eff_cent_np_mid" + weightTag, "Nonprompt MC", midLabel, "Centrality", "Efficiency");
    saveSinglePlot(hNpCentFwd.get(), centNpDir + "/eff_cent_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "Centrality", "Efficiency");

    saveOverlayPlot(hPrCentMid.get(), hNpCentMid.get(), centCompareDir + "/eff_cent_pr_np_mid" + weightTag, midLabel, "Centrality");
    saveOverlayPlot(hPrCentFwd.get(), hNpCentFwd.get(), centCompareDir + "/eff_cent_pr_np_fwd" + weightTag, fwdLabel, "Centrality");
    centSummaryPairs.push_back({hPrCentMid.get(), hNpCentMid.get()});
    centSummaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");
    centSummaryPairs.push_back({hPrCentFwd.get(), hNpCentFwd.get()});
    centSummaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");
  }
  else
  {
    cout << "[WARN] centrality efficiency histograms are not present in the ROOT files." << endl;
    cout << "[WARN] pT efficiency plots were produced, but centrality plots were skipped." << endl;
  }

  if (!centSummaryPairs.empty())
    saveSummaryPlot(centSummaryPairs, centSummaryLabels,
                    summaryCentDir + "/eff_cent_pr_np_summary" + weightTag,
                    "Centrality", "Centrality", true, 1.8);
}
