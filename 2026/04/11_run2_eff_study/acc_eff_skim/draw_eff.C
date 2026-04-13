// Read efficiency files in:
// /data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/acc_eff_skim/skim_roots/
//
// Draw:
// 1. Individual efficiency plots for prompt/nonprompt MC
// 2. Prompt vs nonprompt overlays on the same canvas

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

void drawCmsLabel(const TString &sampleLabel, const TString &rapidityLabel)
{
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.040);
    tx.SetTextFont(72);
    tx.DrawLatex(0.19, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.030);
    tx.SetTextFont(42);
    tx.DrawLatex(0.19, 0.865, sampleLabel);
    tx.DrawLatex(0.19, 0.815, rapidityLabel);
  }
}

void configureAxes(TH1D *hist, const char *yTitle)
{
  hist->SetTitle("");
  if (!hist->GetXaxis()->GetTitle() || TString(hist->GetXaxis()->GetTitle()).Length() == 0)
    hist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hist->GetYaxis()->SetTitle(yTitle);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.25);
}

void saveSinglePlot(TH1D *hist, const TString &outPath, const TString &sampleLabel,
                    const TString &selectionLabel, const char *xTitle,
                    const char *yTitle, const char *legendLabel)
{
  TCanvas canv("c_eff_single", "", 800, 800);
  canv.cd();

  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.SetLeftMargin(0.13);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  pad1.cd();

  gPad->SetTicks(1, 1);
  gStyle->SetOptStat(0);

  configureAxes(hist, yTitle);
  hist->GetXaxis()->SetTitle(xTitle);
  hist->SetMinimum(0.0);
  hist->SetMaximum(std::max(1.05, hist->GetMaximum() * 1.35));
  // E1 visualizes the statistical errors already stored in the input histogram.
  hist->Draw("E1");

  drawCmsLabel(sampleLabel, selectionLabel);

  TLegend leg(0.58, 0.80, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.AddEntry(hist, legendLabel, "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveOverlayPlot(TH1D *histPr, TH1D *histNp, const TString &outPath,
                     const TString &selectionLabel, const char *xTitle)
{
  TCanvas canv("c_eff_overlay", "", 800, 800);
  canv.cd();

  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.SetLeftMargin(0.13);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  pad1.cd();

  gPad->SetTicks(1, 1);
  gStyle->SetOptStat(0);

  std::unique_ptr<TH1D> frame(static_cast<TH1D *>(histPr->Clone(Form("%s_frame", histPr->GetName()))));
  frame->Reset("ICES");
  configureAxes(frame.get(), "Efficiency");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, std::max(histPr->GetMaximum(), histNp->GetMaximum()) * 1.35));
  frame->Draw("AXIS");

  // The drawing macro does not recalculate uncertainties; it only displays them.
  histPr->Draw("E1 SAME");
  histNp->Draw("E1 SAME");

  drawCmsLabel("Prompt / Nonprompt MC", selectionLabel);

  TLegend leg(0.56, 0.76, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.AddEntry(histPr, "Prompt", "lp");
  leg.AddEntry(histNp, "Nonprompt", "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

TH1D *loadHist(TFile *file, const char *histName, const char *cloneName)
{
  TH1D *hist = dynamic_cast<TH1D *>(file->Get(histName));
  if (!hist)
    return nullptr;

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
    return nullptr;

  TH1D *cloned = static_cast<TH1D *>(hist->Clone(cloneName));
  cloned->SetDirectory(nullptr);
  return cloned;
}

void saveSummaryPlot(const std::vector<std::pair<TH1D *, TH1D *>> &histPairs,
                     const std::vector<TString> &selectionLabels,
                     const TString &outPath, const char *xTitle)
{
  if (histPairs.empty() || histPairs.size() != selectionLabels.size())
    return;

  TCanvas canv("c_eff_summary", "", 800, 800);
  canv.cd();

  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.SetLeftMargin(0.13);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  pad1.cd();

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
  configureAxes(frame.get(), "Efficiency");
  frame->GetXaxis()->SetTitle(xTitle);
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, maxY * 1.35));
  frame->Draw("AXIS");

  const Color_t caseColors[] = {kBlack, kRed + 1, kAzure + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                                kCyan + 2, kPink + 7, kGray + 2, kYellow + 2};
  const Style_t prMarkers[] = {20, 21, 22, 23, 29, 33, 34, 47, 43, 41};
  const Style_t npMarkers[] = {24, 25, 26, 32, 30, 27, 28, 46, 42, 40};

  const size_t nLabels = selectionLabels.size();
  const int nColumns = (nLabels > 6) ? 2 : 1;
  const double textSize = (nLabels > 8) ? 0.020 : ((nLabels > 5) ? 0.022 : 0.025);
  const double x1 = (nColumns == 2) ? 0.34 : 0.56;
  const double x2 = 0.94;
  const double y1 = (nLabels > 8) ? 0.48 : ((nLabels > 5) ? 0.54 : 0.66);
  const double y2 = 0.90;
  TLegend leg(x1, y1, x2, y2);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(textSize);
  leg.SetNColumns(nColumns);

  for (size_t i = 0; i < histPairs.size(); ++i)
  {
    TH1D *histPr = histPairs[i].first;
    TH1D *histNp = histPairs[i].second;
    if (!histPr || !histNp)
      continue;

    styleHist(histPr, caseColors[i % 10], prMarkers[i % 10], 1.35);
    styleHist(histNp, caseColors[i % 10], npMarkers[i % 10], 1.35);
    histPr->SetLineStyle(1);
    histNp->SetLineStyle(7);

    histPr->Draw("E1 SAME");
    histNp->Draw("E1 SAME");

    leg.AddEntry(histPr, "PR " + selectionLabels[i], "lp");
    leg.AddEntry(histNp, "NP " + selectionLabels[i], "lp");
  }

  drawCmsLabel("Prompt / Nonprompt MC", "Summary");

  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

void draw_eff(bool isNCollW = true, bool isGenW = true, bool isPtW = false, bool isTnPW = true,
              const char *extraTag = "", const char *figSubdir = "")
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString figsDir = baseDir + "/outputs/figs";
  const TString tag = extraTag ? TString(extraTag) : TString("");
  TString effFigsDir = figsDir + "/eff";
  if (figSubdir && TString(figSubdir).Length() > 0)
    effFigsDir = figsDir + "/" + TString(figSubdir);
  const TString ptDir = effFigsDir + "/pt";
  const TString summaryDir = effFigsDir + "/summary";
  const TString countsDir = effFigsDir + "/counts";
  const TString countsPtDir = countsDir + "/pt";
  const TString effPtDir = ptDir + "/single";
  const TString overlayPtDir = ptDir + "/overlay";
  const TString weightTag = Form("_ncollW%d_genW%d_ptW%d_tnpW%d%s",
                                 isNCollW ? 1 : 0, isGenW ? 1 : 0,
                                 isPtW ? 1 : 0, isTnPW ? 1 : 0, tag.Data());

  gSystem->mkdir(figsDir, true);
  gSystem->mkdir(effFigsDir, true);
  gSystem->mkdir(ptDir, true);
  gSystem->mkdir(summaryDir, true);
  gSystem->mkdir(countsDir, true);
  gSystem->mkdir(countsPtDir, true);
  gSystem->mkdir(effPtDir, true);
  gSystem->mkdir(overlayPtDir, true);

  const TString prPath =
      inputDir + "/eff_PbPb2018_isMC1_PR" + weightTag + "_Dimuon_MiniAOD.root";
  const TString npPath =
      inputDir + "/eff_PbPb2018_isMC1_NP" + weightTag + "_Dimuon_MiniAOD.root";

  std::unique_ptr<TFile> fPr(TFile::Open(prPath, "READ"));
  std::unique_ptr<TFile> fNp(TFile::Open(npPath, "READ"));
  if (!fPr || fPr->IsZombie())
  {
    cout << "[ERROR] failed to open prompt efficiency file: " << prPath << endl;
    return;
  }
  if (!fNp || fNp->IsZombie())
  {
    cout << "[ERROR] failed to open nonprompt efficiency file: " << npPath << endl;
    return;
  }

  std::unique_ptr<TH1D> hPrMid(loadHist(fPr.get(), "hist_eff_mid", "hist_eff_mid_pr"));
  std::unique_ptr<TH1D> hPrFwd(loadHist(fPr.get(), "hist_eff_fwd", "hist_eff_fwd_pr"));
  std::unique_ptr<TH1D> hNpMid(loadHist(fNp.get(), "hist_eff_mid", "hist_eff_mid_np"));
  std::unique_ptr<TH1D> hNpFwd(loadHist(fNp.get(), "hist_eff_fwd", "hist_eff_fwd_np"));
  std::unique_ptr<TH1D> hPrCentMid(loadHist(fPr.get(), "hist_eff_cent_mid", "hist_eff_cent_mid_pr"));
  std::unique_ptr<TH1D> hPrCentFwd(loadHist(fPr.get(), "hist_eff_cent_fwd", "hist_eff_cent_fwd_pr"));
  std::unique_ptr<TH1D> hNpCentMid(loadHist(fNp.get(), "hist_eff_cent_mid", "hist_eff_cent_mid_np"));
  std::unique_ptr<TH1D> hNpCentFwd(loadHist(fNp.get(), "hist_eff_cent_fwd", "hist_eff_cent_fwd_np"));
  std::unique_ptr<TH1D> hPrDenMid(loadHist(fPr.get(), "hist_eff_den_mid", "hist_eff_den_mid_pr"));
  std::unique_ptr<TH1D> hPrDenFwd(loadHist(fPr.get(), "hist_eff_den_fwd", "hist_eff_den_fwd_pr"));
  std::unique_ptr<TH1D> hPrNumMid(loadHist(fPr.get(), "hist_eff_num_mid", "hist_eff_num_mid_pr"));
  std::unique_ptr<TH1D> hPrNumFwd(loadHist(fPr.get(), "hist_eff_num_fwd", "hist_eff_num_fwd_pr"));
  std::unique_ptr<TH1D> hNpDenMid(loadHist(fNp.get(), "hist_eff_den_mid", "hist_eff_den_mid_np"));
  std::unique_ptr<TH1D> hNpDenFwd(loadHist(fNp.get(), "hist_eff_den_fwd", "hist_eff_den_fwd_np"));
  std::unique_ptr<TH1D> hNpNumMid(loadHist(fNp.get(), "hist_eff_num_mid", "hist_eff_num_mid_np"));
  std::unique_ptr<TH1D> hNpNumFwd(loadHist(fNp.get(), "hist_eff_num_fwd", "hist_eff_num_fwd_np"));
  if (!hPrMid || !hPrFwd || !hNpMid || !hNpFwd ||
      !hPrCentMid || !hPrCentFwd || !hNpCentMid || !hNpCentFwd ||
      !hPrDenMid || !hPrDenFwd || !hPrNumMid || !hPrNumFwd ||
      !hNpDenMid || !hNpDenFwd || !hNpNumMid || !hNpNumFwd)
  {
    cout << "[ERROR] missing one or more efficiency histograms in input ROOT files." << endl;
    return;
  }

  styleHist(hPrMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpFwd.get(), kOrange + 7, 25, 1.7);
  styleHist(hPrCentMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrCentFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpCentMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpCentFwd.get(), kOrange + 7, 25, 1.7);
  styleHist(hPrDenMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrDenFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrNumMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrNumFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpDenMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpDenFwd.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpNumMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpNumFwd.get(), kOrange + 7, 25, 1.7);

  const TString midLabel = "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c";
  const TString fwdLabel = "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c";

  saveSinglePlot(hPrDenMid.get(), countsPtDir + "/eff_den_pr_mid" + weightTag, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Counts", "Efficiency denominator");
  saveSinglePlot(hPrDenFwd.get(), countsPtDir + "/eff_den_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Counts", "Efficiency denominator");
  saveSinglePlot(hPrNumMid.get(), countsPtDir + "/eff_num_pr_mid" + weightTag, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Counts", "Efficiency numerator");
  saveSinglePlot(hPrNumFwd.get(), countsPtDir + "/eff_num_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Counts", "Efficiency numerator");
  saveSinglePlot(hNpDenMid.get(), countsPtDir + "/eff_den_np_mid" + weightTag, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Counts", "Efficiency denominator");
  saveSinglePlot(hNpDenFwd.get(), countsPtDir + "/eff_den_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Counts", "Efficiency denominator");
  saveSinglePlot(hNpNumMid.get(), countsPtDir + "/eff_num_np_mid" + weightTag, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Counts", "Efficiency numerator");
  saveSinglePlot(hNpNumFwd.get(), countsPtDir + "/eff_num_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Counts", "Efficiency numerator");

  saveSinglePlot(hPrMid.get(), effPtDir + "/eff_pr_mid" + weightTag, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency", "Efficiency");
  saveSinglePlot(hPrFwd.get(), effPtDir + "/eff_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency", "Efficiency");
  saveSinglePlot(hNpMid.get(), effPtDir + "/eff_np_mid" + weightTag, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Efficiency", "Efficiency");
  saveSinglePlot(hNpFwd.get(), effPtDir + "/eff_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Efficiency", "Efficiency");

  saveOverlayPlot(hPrMid.get(), hNpMid.get(), overlayPtDir + "/eff_pr_np_mid" + weightTag, midLabel, "p_{T} (GeV/c)");
  saveOverlayPlot(hPrFwd.get(), hNpFwd.get(), overlayPtDir + "/eff_pr_np_fwd" + weightTag, fwdLabel, "p_{T} (GeV/c)");

  std::vector<std::unique_ptr<TH1D>> ptSummaryStore;
  std::vector<std::pair<TH1D *, TH1D *>> ptSummaryPairs;
  std::vector<TString> ptSummaryLabels;
  ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrMid->Clone("summary_pr_mid"))));
  ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpMid->Clone("summary_np_mid"))));
  ptSummaryPairs.push_back({ptSummaryStore[ptSummaryStore.size() - 2].get(), ptSummaryStore[ptSummaryStore.size() - 1].get()});
  ptSummaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");
  ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrFwd->Clone("summary_pr_fwd"))));
  ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpFwd->Clone("summary_np_fwd"))));
  ptSummaryPairs.push_back({ptSummaryStore[ptSummaryStore.size() - 2].get(), ptSummaryStore[ptSummaryStore.size() - 1].get()});
  ptSummaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");

  struct PlotCase {
    const char *suffix;
    const char *label;
    const char *summaryLabel;
    const char *xTitle;
  };
  const std::vector<PlotCase> ptCases = {
      {"pt_mid_cent0010", "0-10%, |y| < 1.6", "mid 0-10%", "p_{T} (GeV/c)"},
      {"pt_mid_cent1020", "10-20%, |y| < 1.6", "mid 10-20%", "p_{T} (GeV/c)"},
      {"pt_mid_cent2030", "20-30%, |y| < 1.6", "mid 20-30%", "p_{T} (GeV/c)"},
      {"pt_mid_cent3040", "30-40%, |y| < 1.6", "mid 30-40%", "p_{T} (GeV/c)"},
      {"pt_mid_cent4050", "40-50%, |y| < 1.6", "mid 40-50%", "p_{T} (GeV/c)"},
      {"pt_mid_cent5090", "50-90%, |y| < 1.6", "mid 50-90%", "p_{T} (GeV/c)"},
      {"pt_fwd_cent0010", "0-10%, 1.6 < |y| < 2.4", "fwd 0-10%", "p_{T} (GeV/c)"},
      {"pt_fwd_cent1030", "10-30%, 1.6 < |y| < 2.4", "fwd 10-30%", "p_{T} (GeV/c)"},
      {"pt_fwd_cent3050", "30-50%, 1.6 < |y| < 2.4", "fwd 30-50%", "p_{T} (GeV/c)"},
      {"pt_fwd_cent5090", "50-90%, 1.6 < |y| < 2.4", "fwd 50-90%", "p_{T} (GeV/c)"},
  };

  for (const auto &plotCase : ptCases)
  {
    const TString histEffName = TString::Format("hist_eff_%s", plotCase.suffix);
    const TString histDenName = TString::Format("hist_eff_den_%s", plotCase.suffix);
    const TString histNumName = TString::Format("hist_eff_num_%s", plotCase.suffix);

    std::unique_ptr<TH1D> hPrEff(loadHist(fPr.get(), histEffName, TString::Format("%s_pr", histEffName.Data())));
    std::unique_ptr<TH1D> hNpEff(loadHist(fNp.get(), histEffName, TString::Format("%s_np", histEffName.Data())));
    std::unique_ptr<TH1D> hPrDen(loadHist(fPr.get(), histDenName, TString::Format("%s_pr", histDenName.Data())));
    std::unique_ptr<TH1D> hNpDen(loadHist(fNp.get(), histDenName, TString::Format("%s_np", histDenName.Data())));
    std::unique_ptr<TH1D> hPrNum(loadHist(fPr.get(), histNumName, TString::Format("%s_pr", histNumName.Data())));
    std::unique_ptr<TH1D> hNpNum(loadHist(fNp.get(), histNumName, TString::Format("%s_np", histNumName.Data())));
    if (!hPrEff || !hNpEff || !hPrDen || !hNpDen || !hPrNum || !hNpNum)
    {
      cout << "[WARNING] missing one or more histograms for " << plotCase.suffix << endl;
      continue;
    }

    styleHist(hPrEff.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpEff.get(), kOrange + 7, 25, 1.7);
    styleHist(hPrDen.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpDen.get(), kOrange + 7, 25, 1.7);
    styleHist(hPrNum.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpNum.get(), kOrange + 7, 25, 1.7);

    const TString suffix = plotCase.suffix;
    const TString label = plotCase.label;
    saveSinglePlot(hPrDen.get(), countsPtDir + "/eff_den_pr_" + suffix + weightTag, "Prompt MC", label, plotCase.xTitle, "Counts", "Efficiency denominator");
    saveSinglePlot(hNpDen.get(), countsPtDir + "/eff_den_np_" + suffix + weightTag, "Nonprompt MC", label, plotCase.xTitle, "Counts", "Efficiency denominator");
    saveSinglePlot(hPrNum.get(), countsPtDir + "/eff_num_pr_" + suffix + weightTag, "Prompt MC", label, plotCase.xTitle, "Counts", "Efficiency numerator");
    saveSinglePlot(hNpNum.get(), countsPtDir + "/eff_num_np_" + suffix + weightTag, "Nonprompt MC", label, plotCase.xTitle, "Counts", "Efficiency numerator");
    saveSinglePlot(hPrEff.get(), effPtDir + "/eff_pr_" + suffix + weightTag, "Prompt MC", label, plotCase.xTitle, "Efficiency", "Efficiency");
    saveSinglePlot(hNpEff.get(), effPtDir + "/eff_np_" + suffix + weightTag, "Nonprompt MC", label, plotCase.xTitle, "Efficiency", "Efficiency");
    saveOverlayPlot(hPrEff.get(), hNpEff.get(), overlayPtDir + "/eff_pr_np_" + suffix + weightTag, label, plotCase.xTitle);

    ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrEff->Clone(TString::Format("summary_pr_%s", plotCase.suffix)))));
    ptSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpEff->Clone(TString::Format("summary_np_%s", plotCase.suffix)))));
    ptSummaryPairs.push_back({ptSummaryStore[ptSummaryStore.size() - 2].get(), ptSummaryStore[ptSummaryStore.size() - 1].get()});
    ptSummaryLabels.push_back(plotCase.summaryLabel);
  }

  saveSummaryPlot(ptSummaryPairs, ptSummaryLabels, summaryDir + "/eff_pr_np_pt_summary" + weightTag, "p_{T} (GeV/c)");

  const std::vector<PlotCase> centCases = {};

  std::vector<std::unique_ptr<TH1D>> centSummaryStore;
  std::vector<std::pair<TH1D *, TH1D *>> centSummaryPairs;
  std::vector<TString> centSummaryLabels;
  centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrCentMid->Clone("summary_cent_pr_mid"))));
  centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpCentMid->Clone("summary_cent_np_mid"))));
  centSummaryPairs.push_back({centSummaryStore[centSummaryStore.size() - 2].get(), centSummaryStore[centSummaryStore.size() - 1].get()});
  centSummaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");
  centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrCentFwd->Clone("summary_cent_pr_fwd"))));
  centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpCentFwd->Clone("summary_cent_np_fwd"))));
  centSummaryPairs.push_back({centSummaryStore[centSummaryStore.size() - 2].get(), centSummaryStore[centSummaryStore.size() - 1].get()});
  centSummaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");

  for (const auto &plotCase : centCases)
  {
    const TString histEffName = TString::Format("hist_eff_%s", plotCase.suffix);
    const TString histDenName = TString::Format("hist_eff_den_%s", plotCase.suffix);
    const TString histNumName = TString::Format("hist_eff_num_%s", plotCase.suffix);

    std::unique_ptr<TH1D> hPrEff(loadHist(fPr.get(), histEffName, TString::Format("%s_pr", histEffName.Data())));
    std::unique_ptr<TH1D> hNpEff(loadHist(fNp.get(), histEffName, TString::Format("%s_np", histEffName.Data())));
    std::unique_ptr<TH1D> hPrDen(loadHist(fPr.get(), histDenName, TString::Format("%s_pr", histDenName.Data())));
    std::unique_ptr<TH1D> hNpDen(loadHist(fNp.get(), histDenName, TString::Format("%s_np", histDenName.Data())));
    std::unique_ptr<TH1D> hPrNum(loadHist(fPr.get(), histNumName, TString::Format("%s_pr", histNumName.Data())));
    std::unique_ptr<TH1D> hNpNum(loadHist(fNp.get(), histNumName, TString::Format("%s_np", histNumName.Data())));
    if (!hPrEff || !hNpEff || !hPrDen || !hNpDen || !hPrNum || !hNpNum)
    {
      cout << "[WARNING] missing one or more histograms for " << plotCase.suffix << endl;
      continue;
    }

    styleHist(hPrEff.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpEff.get(), kOrange + 7, 25, 1.7);
    styleHist(hPrDen.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpDen.get(), kOrange + 7, 25, 1.7);
    styleHist(hPrNum.get(), kAzure + 2, 24, 1.7);
    styleHist(hNpNum.get(), kOrange + 7, 25, 1.7);

    const TString suffix = plotCase.suffix;
    const TString label = plotCase.label;
    centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hPrEff->Clone(TString::Format("summary_cent_pr_%s", plotCase.suffix)))));
    centSummaryStore.push_back(std::unique_ptr<TH1D>(static_cast<TH1D *>(hNpEff->Clone(TString::Format("summary_cent_np_%s", plotCase.suffix)))));
    centSummaryPairs.push_back({centSummaryStore[centSummaryStore.size() - 2].get(), centSummaryStore[centSummaryStore.size() - 1].get()});
    centSummaryLabels.push_back(plotCase.summaryLabel);
  }

  saveSummaryPlot(centSummaryPairs, centSummaryLabels, summaryDir + "/eff_pr_np_cent_summary" + weightTag, "Centrality");

}
