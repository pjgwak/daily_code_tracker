// Read acceptance files in:
// /data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/acc_eff_skim/skim_roots/
//
// Draw:
// 1. Individual acceptance plots for prompt/nonprompt MC in mid/forward rapidity
// 2. Prompt vs nonprompt overlays on the same canvas for each rapidity region

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <iostream>
#include <memory>

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
                    const TString &rapidityLabel, const char *yTitle, const char *legendLabel)
{
  TCanvas canv("c_acc_single", "", 800, 800);
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
  hist->SetMinimum(0.0);
  hist->SetMaximum(std::max(1.05, hist->GetMaximum() * 1.35));
  hist->Draw("E1");

  drawCmsLabel(sampleLabel, rapidityLabel);

  TLegend leg(0.58, 0.80, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.AddEntry(hist, legendLabel, "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}

void saveOverlayPlot(TH1D *histPr, TH1D *histNp, const TString &outPath, const TString &rapidityLabel)
{
  TCanvas canv("c_acc_overlay", "", 800, 800);
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
  configureAxes(frame.get(), "Acceptance");
  frame->SetMinimum(0.0);
  frame->SetMaximum(std::max(1.05, std::max(histPr->GetMaximum(), histNp->GetMaximum()) * 1.35));
  frame->Draw("AXIS");

  histPr->Draw("E1 SAME");
  histNp->Draw("E1 SAME");

  drawCmsLabel("Prompt / Nonprompt MC", rapidityLabel);

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
  return cloned;
}
} // namespace

void draw_acc(bool isNCollW = false, bool isGenW = true, bool isPtW = true)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir = resolveBaseDir();
  const TString inputDir = baseDir + "/skim_roots";
  const TString figsDir = baseDir + "/outputs/figs";
  const TString accFigsDir = figsDir + "/acc";
  const TString effPrDir = accFigsDir + "/acceptance/pr";
  const TString effNpDir = accFigsDir + "/acceptance/np";
  const TString effCompareDir = accFigsDir + "/acceptance/compare";
  const TString denPrDir = accFigsDir + "/denominator/pr";
  const TString denNpDir = accFigsDir + "/denominator/np";
  const TString numPrDir = accFigsDir + "/numerator/pr";
  const TString numNpDir = accFigsDir + "/numerator/np";
  const TString weightTag = Form("_ncollW%d_genW%d_ptW%d",
                                 isNCollW ? 1 : 0, isGenW ? 1 : 0, isPtW ? 1 : 0);

  gSystem->mkdir(figsDir, true);
  gSystem->mkdir(accFigsDir, true);
  gSystem->mkdir(effPrDir, true);
  gSystem->mkdir(effNpDir, true);
  gSystem->mkdir(effCompareDir, true);
  gSystem->mkdir(denPrDir, true);
  gSystem->mkdir(denNpDir, true);
  gSystem->mkdir(numPrDir, true);
  gSystem->mkdir(numNpDir, true);

  const TString prPath =
      inputDir + "/acc_OO2025_isMC1_PR" + weightTag + "_Dimuon_MiniAOD_Private_MC.root";
  const TString npPath =
      inputDir + "/acc_OO2025_isMC1_NP" + weightTag + "_Dimuon_MiniAOD_Private_MC.root";

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

  std::unique_ptr<TH1D> hPrMid(loadHist(fPr.get(), "hist_acc_mid", "hist_acc_mid_pr"));
  std::unique_ptr<TH1D> hPrFwd(loadHist(fPr.get(), "hist_acc_fwd", "hist_acc_fwd_pr"));
  std::unique_ptr<TH1D> hNpMid(loadHist(fNp.get(), "hist_acc_mid", "hist_acc_mid_np"));
  std::unique_ptr<TH1D> hNpFwd(loadHist(fNp.get(), "hist_acc_fwd", "hist_acc_fwd_np"));
  std::unique_ptr<TH1D> hPrDenMid(loadHist(fPr.get(), "hist_den_mid", "hist_den_mid_pr"));
  std::unique_ptr<TH1D> hPrDenFwd(loadHist(fPr.get(), "hist_den_fwd", "hist_den_fwd_pr"));
  std::unique_ptr<TH1D> hPrNumMid(loadHist(fPr.get(), "hist_num_mid", "hist_num_mid_pr"));
  std::unique_ptr<TH1D> hPrNumFwd(loadHist(fPr.get(), "hist_num_fwd", "hist_num_fwd_pr"));
  std::unique_ptr<TH1D> hNpDenMid(loadHist(fNp.get(), "hist_den_mid", "hist_den_mid_np"));
  std::unique_ptr<TH1D> hNpDenFwd(loadHist(fNp.get(), "hist_den_fwd", "hist_den_fwd_np"));
  std::unique_ptr<TH1D> hNpNumMid(loadHist(fNp.get(), "hist_num_mid", "hist_num_mid_np"));
  std::unique_ptr<TH1D> hNpNumFwd(loadHist(fNp.get(), "hist_num_fwd", "hist_num_fwd_np"));
  if (!hPrMid || !hPrFwd || !hNpMid || !hNpFwd ||
      !hPrDenMid || !hPrDenFwd || !hPrNumMid || !hPrNumFwd ||
      !hNpDenMid || !hNpDenFwd || !hNpNumMid || !hNpNumFwd)
  {
    cout << "[ERROR] missing one or more acceptance histograms in input ROOT files." << endl;
    return;
  }

  styleHist(hPrMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpFwd.get(), kOrange + 7, 25, 1.7);
  styleHist(hPrDenMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrDenFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrNumMid.get(), kAzure + 2, 24, 1.7);
  styleHist(hPrNumFwd.get(), kAzure + 2, 24, 1.7);
  styleHist(hNpDenMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpDenFwd.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpNumMid.get(), kOrange + 7, 25, 1.7);
  styleHist(hNpNumFwd.get(), kOrange + 7, 25, 1.7);

  const TString midLabel = "|y| < 1.6";
  const TString fwdLabel = "1.6 < |y| < 2.4";

  saveSinglePlot(hPrDenMid.get(), denPrDir + "/acc_den_pr_mid" + weightTag, "Prompt MC", midLabel, "Counts", "Acceptance denominator");
  saveSinglePlot(hPrDenFwd.get(), denPrDir + "/acc_den_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "Counts", "Acceptance denominator");
  saveSinglePlot(hPrNumMid.get(), numPrDir + "/acc_num_pr_mid" + weightTag, "Prompt MC", midLabel, "Counts", "Acceptance numerator");
  saveSinglePlot(hPrNumFwd.get(), numPrDir + "/acc_num_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "Counts", "Acceptance numerator");
  saveSinglePlot(hNpDenMid.get(), denNpDir + "/acc_den_np_mid" + weightTag, "Nonprompt MC", midLabel, "Counts", "Acceptance denominator");
  saveSinglePlot(hNpDenFwd.get(), denNpDir + "/acc_den_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "Counts", "Acceptance denominator");
  saveSinglePlot(hNpNumMid.get(), numNpDir + "/acc_num_np_mid" + weightTag, "Nonprompt MC", midLabel, "Counts", "Acceptance numerator");
  saveSinglePlot(hNpNumFwd.get(), numNpDir + "/acc_num_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "Counts", "Acceptance numerator");

  saveSinglePlot(hPrMid.get(), effPrDir + "/acc_pr_mid" + weightTag, "Prompt MC", midLabel, "Acceptance", "Acceptance");
  saveSinglePlot(hPrFwd.get(), effPrDir + "/acc_pr_fwd" + weightTag, "Prompt MC", fwdLabel, "Acceptance", "Acceptance");
  saveSinglePlot(hNpMid.get(), effNpDir + "/acc_np_mid" + weightTag, "Nonprompt MC", midLabel, "Acceptance", "Acceptance");
  saveSinglePlot(hNpFwd.get(), effNpDir + "/acc_np_fwd" + weightTag, "Nonprompt MC", fwdLabel, "Acceptance", "Acceptance");

  saveOverlayPlot(hPrMid.get(), hNpMid.get(), effCompareDir + "/acc_pr_np_mid" + weightTag, midLabel);
  saveOverlayPlot(hPrFwd.get(), hNpFwd.get(), effCompareDir + "/acc_pr_np_fwd" + weightTag, fwdLabel);
}
