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

void drawCmsLabel(const TString &sampleLabel, const TString &regionLabel)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);

  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.935, "pp 5.02 TeV");

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

std::unique_ptr<TH1D> cloneDetached(TH1D *hist, const TString &name)
{
  std::unique_ptr<TH1D> cloned(static_cast<TH1D *>(hist->Clone(name)));
  cloned->SetDirectory(nullptr);
  return cloned;
}

void saveSinglePlot(TH1D *hist, const TString &outPath, const TString &sampleLabel,
                    const TString &regionLabel, const char *xTitle, const char *yTitle)
{
  TCanvas canv("c_eff_single_pp", "", 800, 800);
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
  hist->Draw("E1");

  drawCmsLabel(sampleLabel, regionLabel);

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveOverlayPlot(TH1D *histPr, TH1D *histNp, const TString &outPath,
                     const TString &regionLabel, const char *xTitle)
{
  TCanvas canv("c_eff_overlay_pp", "", 800, 800);
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
                     const TString &headerLabel)
{
  if (histPairs.empty() || histPairs.size() != regionLabels.size())
    return;

  TCanvas canv("c_eff_summary_pp", "", 800, 800);
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

  std::unique_ptr<TH1D> frame(static_cast<TH1D *>(histPairs.front().first->Clone("summary_frame_pp")));
  frame->Reset("ICES");
  configureAxes(frame.get(), xTitle, "Efficiency");
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, maxY * 1.35));
  frame->Draw("AXIS");

  const Color_t caseColors[] = {kBlack, kRed + 1, kAzure + 1, kGreen + 2, kMagenta + 1,
                                kOrange + 7, kCyan + 2, kPink + 7};
  const Style_t prMarkers[] = {20, 21, 22, 23, 29, 33, 34, 47};
  const Style_t npMarkers[] = {24, 25, 26, 32, 30, 27, 28, 46};

  const int nColumns = (regionLabels.size() > 4) ? 2 : 1;
  TLegend leg((nColumns == 2) ? 0.36 : 0.56, 0.64, 0.94, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize((regionLabels.size() > 4) ? 0.022 : 0.025);
  leg.SetNColumns(nColumns);

  for (size_t i = 0; i < histPairs.size(); ++i)
  {
    TH1D *histPr = histPairs[i].first;
    TH1D *histNp = histPairs[i].second;
    if (!histPr || !histNp)
      continue;

    styleHist(histPr, caseColors[i % 8], prMarkers[i % 8], 1.45);
    styleHist(histNp, caseColors[i % 8], npMarkers[i % 8], 1.45);
    histPr->SetLineStyle(1);
    histNp->SetLineStyle(7);
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

void draw_eff_pp(bool isNCollW = true, bool isGenW = true, bool isPtW = true, bool isTnPW = true)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString outDir = baseDir + "/outputs/figs/eff_pp";
  const TString ptDir = outDir + "/pt";
  const TString ptSingleDir = ptDir + "/single";
  const TString ptOverlayDir = ptDir + "/overlay";
  const TString yDir = outDir + "/y";
  const TString ySingleDir = yDir + "/single";
  const TString yOverlayDir = yDir + "/overlay";
  const TString summaryDir = outDir + "/summary";
  const TString weightTag = Form("_ncollW%d_genW%d_ptW%d_tnpW%d",
                                 isNCollW ? 1 : 0, isGenW ? 1 : 0,
                                 isPtW ? 1 : 0, isTnPW ? 1 : 0);

  gSystem->mkdir(outDir, true);
  gSystem->mkdir(ptDir, true);
  gSystem->mkdir(ptSingleDir, true);
  gSystem->mkdir(ptOverlayDir, true);
  gSystem->mkdir(yDir, true);
  gSystem->mkdir(ySingleDir, true);
  gSystem->mkdir(yOverlayDir, true);
  gSystem->mkdir(summaryDir, true);

  const TString prPath =
      inputDir + "/eff_pp5p02TeV_isMC1_PR" + weightTag + "_Dimuon_MiniAOD.root";
  const TString npPath =
      inputDir + "/eff_pp5p02TeV_isMC1_NP" + weightTag + "_Dimuon_MiniAOD.root";

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

  styleHist(hPrPtMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrPtFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpPtMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpPtFwd.get(), kOrange + 7, 25, 1.7);

  const TString midLabel = "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c";
  const TString fwdLabel = "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c";

  saveSinglePlot(hPrPtMid.get(), ptSingleDir + "/eff_pr_mid" + weightTag, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hPrPtFwd.get(), ptSingleDir + "/eff_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hNpPtMid.get(), ptSingleDir + "/eff_np_mid" + weightTag, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency");
  saveSinglePlot(hNpPtFwd.get(), ptSingleDir + "/eff_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency");

  saveOverlayPlot(hPrPtMid.get(), hNpPtMid.get(), ptOverlayDir + "/eff_pr_np_mid" + weightTag, midLabel, "p_{T} (GeV/c)");
  saveOverlayPlot(hPrPtFwd.get(), hNpPtFwd.get(), ptOverlayDir + "/eff_pr_np_fwd" + weightTag, fwdLabel, "p_{T} (GeV/c)");

  std::vector<std::unique_ptr<TH1D>> ptSummaryStore;
  std::vector<std::pair<TH1D *, TH1D *>> ptSummaryPairs;
  std::vector<TString> ptSummaryLabels;

  ptSummaryStore.push_back(cloneDetached(hPrPtMid.get(), "summary_pr_mid"));
  ptSummaryStore.push_back(cloneDetached(hNpPtMid.get(), "summary_np_mid"));
  ptSummaryPairs.emplace_back(ptSummaryStore[ptSummaryStore.size() - 2].get(),
                              ptSummaryStore[ptSummaryStore.size() - 1].get());
  ptSummaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");

  ptSummaryStore.push_back(cloneDetached(hPrPtFwd.get(), "summary_pr_fwd"));
  ptSummaryStore.push_back(cloneDetached(hNpPtFwd.get(), "summary_np_fwd"));
  ptSummaryPairs.emplace_back(ptSummaryStore[ptSummaryStore.size() - 2].get(),
                              ptSummaryStore[ptSummaryStore.size() - 1].get());
  ptSummaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");

  saveSummaryPlot(ptSummaryPairs, ptSummaryLabels,
                  summaryDir + "/eff_pt_pr_np_summary" + weightTag,
                  "p_{T} (GeV/c)", "p_{T} summary");
}
