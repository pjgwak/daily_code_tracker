#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
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

void drawCmsLabel(const TString &collisionLabel, const TString &sampleLabel, const TString &regionLabel)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);

  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.935, collisionLabel);

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

double integratedAcceptance(TH1D *num, TH1D *den)
{
  if (!num || !den)
    return 0.0;

  const double denIntegral = den->Integral();
  if (denIntegral <= 0.0)
    return 0.0;

  return num->Integral() / denIntegral;
}

void saveSinglePlot(TH1D *hist, const TString &outPath,
                    const TString &collisionLabel, const TString &sampleLabel,
                    const TString &regionLabel, const char *xTitle, const char *yTitle,
                    double integratedAcc)
{
  TCanvas canv("c_acc_single", "", 800, 800);
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

  TLine guide(hist->GetXaxis()->GetXmin(), integratedAcc,
              hist->GetXaxis()->GetXmax(), integratedAcc);
  guide.SetLineColor(kGray + 2);
  guide.SetLineStyle(7);
  guide.SetLineWidth(3);
  guide.Draw("SAME");

  drawCmsLabel(collisionLabel, sampleLabel, regionLabel);

  TLegend leg(0.56, 0.70, 0.94, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.024);
  leg.AddEntry(hist, "Acceptance by p_{T} bin", "lp");
  leg.AddEntry(&guide, "Integrated acc", "l");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveOverlayPlot(TH1D *histPr, TH1D *histNp, const TString &outPath,
                     const TString &collisionLabel, const TString &regionLabel,
                     const char *xTitle, const char *yTitle,
                     double integratedPr, double integratedNp)
{
  TCanvas canv("c_acc_overlay", "", 800, 800);
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
  configureAxes(frame.get(), xTitle, yTitle);
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, std::max(histPr->GetMaximum(), histNp->GetMaximum()) * 1.35));
  frame->Draw("AXIS");

  histPr->Draw("E1 SAME");
  histNp->Draw("E1 SAME");

  TLine linePr(histPr->GetXaxis()->GetXmin(), integratedPr,
               histPr->GetXaxis()->GetXmax(), integratedPr);
  linePr.SetLineColor(histPr->GetLineColor());
  linePr.SetLineStyle(7);
  linePr.SetLineWidth(3);
  linePr.Draw("SAME");

  TLine lineNp(histNp->GetXaxis()->GetXmin(), integratedNp,
               histNp->GetXaxis()->GetXmax(), integratedNp);
  lineNp.SetLineColor(histNp->GetLineColor());
  lineNp.SetLineStyle(3);
  lineNp.SetLineWidth(3);
  lineNp.Draw("SAME");

  drawCmsLabel(collisionLabel, "Prompt / Nonprompt MC", regionLabel);

  TLegend leg(0.56, 0.66, 0.94, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.024);
  leg.AddEntry(histPr, "Prompt", "lp");
  leg.AddEntry(histNp, "Nonprompt", "lp");
  leg.AddEntry(&linePr, "PR integrated", "l");
  leg.AddEntry(&lineNp, "NP integrated", "l");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveSummaryPlot(const std::vector<std::pair<TH1D *, TH1D *>> &histPairs,
                     const std::vector<TString> &regionLabels,
                     const TString &outPath, const TString &collisionLabel,
                     const char *xTitle, const char *yTitle, const TString &headerLabel)
{
  if (histPairs.empty() || histPairs.size() != regionLabels.size())
    return;

  TCanvas canv("c_acc_summary", "", 800, 800);
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
  double minX = 1e9;
  double maxX = -1e9;
  for (const auto &pair : histPairs)
  {
    if (pair.first)
    {
      maxY = std::max(maxY, pair.first->GetMaximum());
      minX = std::min(minX, pair.first->GetXaxis()->GetBinLowEdge(1));
      maxX = std::max(maxX, pair.first->GetXaxis()->GetBinUpEdge(pair.first->GetNbinsX()));
    }
    if (pair.second)
    {
      maxY = std::max(maxY, pair.second->GetMaximum());
      minX = std::min(minX, pair.second->GetXaxis()->GetBinLowEdge(1));
      maxX = std::max(maxX, pair.second->GetXaxis()->GetBinUpEdge(pair.second->GetNbinsX()));
    }
  }

  if (minX >= maxX)
    return;

  std::unique_ptr<TH1D> frame(new TH1D("summary_frame_acc", "", 100, minX, maxX));
  configureAxes(frame.get(), xTitle, yTitle);
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, maxY * 1.35));
  frame->Draw("AXIS");

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

    const bool isMidRegion = regionLabels[i].Contains("|y|<1.6");
    const Style_t regionMarker = isMidRegion ? 27 : 24;

    styleHist(histPr, kAzure + 2, regionMarker, 1.45);
    styleHist(histNp, kRed + 1, regionMarker, 1.45);
    histPr->SetLineStyle(1);
    histNp->SetLineStyle(7);
    histPr->Draw("E1 SAME");
    histNp->Draw("E1 SAME");

    leg.AddEntry(histPr, "PR " + regionLabels[i], "lp");
    leg.AddEntry(histNp, "NP " + regionLabels[i], "lp");
  }

  drawCmsLabel(collisionLabel, "Prompt / Nonprompt MC", headerLabel);
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

void draw_acc(bool isPbPb = false, bool isGenW = true, bool isPtW = true)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString collLabel = isPbPb ? "PbPb2018" : "pp2018";
  const TString cmsLabel = isPbPb ? "PbPb 5.02 TeV" : "pp 5.02 TeV";
  const TString outDir = baseDir + "/outputs/figs/" + (isPbPb ? TString("acc_pbpb") : TString("acc_pp"));
  const TString ptDir = outDir + "/pt";
  const TString ptSingleDir = ptDir + "/single";
  const TString ptOverlayDir = ptDir + "/overlay";
  const TString summaryDir = outDir + "/summary";
  const TString weightTag = Form("_ncollW0_genW%d_ptW%d", isGenW ? 1 : 0, isPtW ? 1 : 0);

  gSystem->mkdir(outDir, true);
  gSystem->mkdir(ptDir, true);
  gSystem->mkdir(ptSingleDir, true);
  gSystem->mkdir(ptOverlayDir, true);
  gSystem->mkdir(summaryDir, true);

  const TString prPath =
      inputDir + "/acc_" + collLabel + "_ppInput_isMC1_PR" + weightTag + ".root";
  const TString npPath =
      inputDir + "/acc_" + collLabel + "_ppInput_isMC1_NP" + weightTag + ".root";

  std::unique_ptr<TFile> fPr(TFile::Open(prPath, "READ"));
  std::unique_ptr<TFile> fNp(TFile::Open(npPath, "READ"));
  if (!fPr || fPr->IsZombie())
  {
    cout << "[ERROR] failed to open prompt acceptance file: " << prPath << endl;
    return;
  }
  if (!fNp || fNp->IsZombie())
  {
    cout << "[ERROR] failed to open nonprompt acceptance file: " << npPath << endl;
    return;
  }

  std::unique_ptr<TH1D> hPrMid(loadHist(fPr.get(), "hist_acc_mid", "hPrAccMid"));
  std::unique_ptr<TH1D> hPrFwd(loadHist(fPr.get(), "hist_acc_fwd", "hPrAccFwd"));
  std::unique_ptr<TH1D> hNpMid(loadHist(fNp.get(), "hist_acc_mid", "hNpAccMid"));
  std::unique_ptr<TH1D> hNpFwd(loadHist(fNp.get(), "hist_acc_fwd", "hNpAccFwd"));
  std::unique_ptr<TH1D> hPrNumMid(loadHist(fPr.get(), "hist_acc_num_mid", "hPrAccNumMid"));
  std::unique_ptr<TH1D> hPrDenMid(loadHist(fPr.get(), "hist_acc_den_mid", "hPrAccDenMid"));
  std::unique_ptr<TH1D> hPrNumFwd(loadHist(fPr.get(), "hist_acc_num_fwd", "hPrAccNumFwd"));
  std::unique_ptr<TH1D> hPrDenFwd(loadHist(fPr.get(), "hist_acc_den_fwd", "hPrAccDenFwd"));
  std::unique_ptr<TH1D> hNpNumMid(loadHist(fNp.get(), "hist_acc_num_mid", "hNpAccNumMid"));
  std::unique_ptr<TH1D> hNpDenMid(loadHist(fNp.get(), "hist_acc_den_mid", "hNpAccDenMid"));
  std::unique_ptr<TH1D> hNpNumFwd(loadHist(fNp.get(), "hist_acc_num_fwd", "hNpAccNumFwd"));
  std::unique_ptr<TH1D> hNpDenFwd(loadHist(fNp.get(), "hist_acc_den_fwd", "hNpAccDenFwd"));
  if (!hPrMid || !hPrFwd || !hNpMid || !hNpFwd)
    return;
  if (!hPrNumMid || !hPrDenMid || !hPrNumFwd || !hPrDenFwd || !hNpNumMid || !hNpDenMid || !hNpNumFwd || !hNpDenFwd)
    return;

  const double prMidIntegrated = integratedAcceptance(hPrNumMid.get(), hPrDenMid.get());
  const double prFwdIntegrated = integratedAcceptance(hPrNumFwd.get(), hPrDenFwd.get());
  const double npMidIntegrated = integratedAcceptance(hNpNumMid.get(), hNpDenMid.get());
  const double npFwdIntegrated = integratedAcceptance(hNpNumFwd.get(), hNpDenFwd.get());

  styleHist(hPrMid.get(), kAzure + 2, 27, 1.7);
  styleHist(hPrFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpMid.get(), kRed + 1, 27, 1.7);
  styleHist(hNpFwd.get(), kRed + 1, 24, 1.7);

  const TString midLabel = "|y| < 1.6, 6.5 < p_{T} < 40 GeV/c";
  const TString fwdLabel = "1.6 < |y| < 2.4, 3.5 < p_{T} < 40 GeV/c";

  saveSinglePlot(hPrMid.get(), ptSingleDir + "/acc_pr_mid" + weightTag,
                 cmsLabel, "Prompt MC", midLabel, "p_{T} (GeV/c)", "Acceptance", prMidIntegrated);
  saveSinglePlot(hPrFwd.get(), ptSingleDir + "/acc_pr_fwd" + weightTag,
                 cmsLabel, "Prompt MC", fwdLabel, "p_{T} (GeV/c)", "Acceptance", prFwdIntegrated);
  saveSinglePlot(hNpMid.get(), ptSingleDir + "/acc_np_mid" + weightTag,
                 cmsLabel, "Nonprompt MC", midLabel, "p_{T} (GeV/c)", "Acceptance", npMidIntegrated);
  saveSinglePlot(hNpFwd.get(), ptSingleDir + "/acc_np_fwd" + weightTag,
                 cmsLabel, "Nonprompt MC", fwdLabel, "p_{T} (GeV/c)", "Acceptance", npFwdIntegrated);

  saveOverlayPlot(hPrMid.get(), hNpMid.get(), ptOverlayDir + "/acc_pr_np_mid" + weightTag,
                  cmsLabel, midLabel, "p_{T} (GeV/c)", "Acceptance", prMidIntegrated, npMidIntegrated);
  saveOverlayPlot(hPrFwd.get(), hNpFwd.get(), ptOverlayDir + "/acc_pr_np_fwd" + weightTag,
                  cmsLabel, fwdLabel, "p_{T} (GeV/c)", "Acceptance", prFwdIntegrated, npFwdIntegrated);

  std::vector<std::unique_ptr<TH1D>> summaryStore;
  std::vector<std::pair<TH1D *, TH1D *>> summaryPairs;
  std::vector<TString> summaryLabels;

  summaryStore.push_back(cloneDetached(hPrMid.get(), "summary_pr_mid"));
  summaryStore.push_back(cloneDetached(hNpMid.get(), "summary_np_mid"));
  summaryPairs.emplace_back(summaryStore[summaryStore.size() - 2].get(),
                            summaryStore[summaryStore.size() - 1].get());
  summaryLabels.push_back("|y|<1.6, 6.5<p_{T}<40");

  summaryStore.push_back(cloneDetached(hPrFwd.get(), "summary_pr_fwd"));
  summaryStore.push_back(cloneDetached(hNpFwd.get(), "summary_np_fwd"));
  summaryPairs.emplace_back(summaryStore[summaryStore.size() - 2].get(),
                            summaryStore[summaryStore.size() - 1].get());
  summaryLabels.push_back("1.6<|y|<2.4, 3.5<p_{T}<40");

  saveSummaryPlot(summaryPairs, summaryLabels,
                  summaryDir + "/acc_summary" + weightTag,
                  cmsLabel, "Observable", "Acceptance", "Acceptance summary");
}
