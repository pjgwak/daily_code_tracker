
#include <algorithm>
#include <cmath>
#include <memory>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

void mass_ttree_draw_only(float ptLow = 3.5, float ptHigh = 50,
                          float yLow = 0, float yHigh = 2.4,
                          int bins = 90)
{
  const char *pdName = "PromptReco";
  const char *runLabel = "404337_404474";
  const char *triggerName = "L1SingleMu0_Open";
  const char *dataRoot = "/data/users/pjgwak/work/daily_code_tracker/2026/05/27_PbPb_run_muon_performance/mass_fitting/PromptReco/404337_404468/roodataset_roots/PromptReco/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_PromptReco_404337_404474_SingleMu0_Open.root";

  gStyle->SetOptStat(0);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  auto inputFile = std::unique_ptr<TFile>(TFile::Open(dataRoot, "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    printf("Cannot open input file:\n  %s\n", dataRoot);
    return;
  }

  auto tree = dynamic_cast<TTree *>(inputFile->Get("RooTreeDataStore_dataset_"));
  if (!tree) {
    printf("Cannot find TTree RooTreeDataStore_dataset_\n");
    return;
  }

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("mass", 1);
  tree->SetBranchStatus("pt", 1);
  tree->SetBranchStatus("y", 1);

  double mass = 0.0;
  double pt = 0.0;
  double y = 0.0;
  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("y", &y);

  auto formatTag = [](double value) { return TString::Format("%.2f", value); };
  const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
  const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
  const TString figTag = yTag + "_" + ptTag;

  TH1D hMass("hMass", "", bins, 2.6, 3.5);
  hMass.Sumw2();

  const Long64_t nEntries = tree->GetEntries();
  Long64_t selected = 0;
  const Long64_t reportEvery = 5000000;
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (!(mass > 2.6 && mass < 3.5 &&
          pt > ptLow && pt < ptHigh &&
          std::abs(y) > yLow && std::abs(y) < yHigh)) {
      continue;
    }
    hMass.Fill(mass);
    ++selected;
    if ((selected % reportEvery) == 0) {
      printf("Selected %lld entries\n", selected);
    }
  }

  printf("Draw-only selected entries: %lld / %lld\n", selected, nEntries);
  printf("Cut: mass > 2.6 && mass < 3.5 && pt > %g && pt < %g && abs(y) > %g && abs(y) < %g\n",
         ptLow, ptHigh, yLow, yHigh);

  const TString figDir = TString::Format("figs/%s/mass_draw", yTag.Data());
  const TString rootDir = TString::Format("roots/%s/mass_draw", yTag.Data());
  gSystem->mkdir(figDir, true);
  gSystem->mkdir(rootDir, true);

  TCanvas canvas("canvas", "mass draw only", 800, 700);
  canvas.SetLeftMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.08);
  canvas.SetBottomMargin(0.12);

  hMass.SetMarkerStyle(20);
  hMass.SetMarkerSize(0.7);
  hMass.SetLineColor(kBlack);
  hMass.GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
  hMass.GetYaxis()->SetTitle("Events");
  hMass.GetXaxis()->CenterTitle();
  hMass.GetYaxis()->CenterTitle();
  hMass.GetYaxis()->SetTitleOffset(1.15);
  hMass.SetMaximum(hMass.GetMaximum() * 1.35);
  hMass.Draw("E");

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.935, Form("PbPb %s run %s", pdName, runLabel));

  tx.SetTextAlign(11);
  tx.SetTextFont(72);
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.21, 0.935, "CMS Internal");

  tx.SetTextFont(42);
  tx.SetTextSize(0.030);
  double xtext = 0.18;
  double y0 = 0.86;
  double dy = -0.045;
  int k = 0;
  tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
  if (yLow == 0)
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
  else
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
  tx.DrawLatex(xtext, y0 + dy * k++, triggerName);
  // tx.DrawLatex(xtext, y0 + dy * k++, Form("Entries = %lld", selected));

  const TString pdfName = Form("%s/mass_draw_%s.pdf", figDir.Data(), figTag.Data());
  const TString pngName = Form("%s/mass_draw_%s.png", figDir.Data(), figTag.Data());
  const TString rootName = Form("%s/mass_draw_%s.root", rootDir.Data(), figTag.Data());
  canvas.SaveAs(pdfName);
  canvas.SaveAs(pngName);

  TCanvas logCanvas("logCanvas", "mass draw only log", 800, 700);
  logCanvas.SetLeftMargin(0.13);
  logCanvas.SetRightMargin(0.04);
  logCanvas.SetTopMargin(0.08);
  logCanvas.SetBottomMargin(0.12);
  logCanvas.SetLogy(true);
  hMass.SetMinimum(0.5);
  hMass.SetMaximum(hMass.GetMaximum() * 20.0);
  hMass.Draw("E");
  tx.DrawLatex(0.96, 0.935, Form("PbPb %s run %s", pdName, runLabel));
  tx.SetTextAlign(11);
  tx.SetTextFont(72);
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.18, 0.935, "CMS Internal");
  tx.SetTextFont(42);
  tx.SetTextSize(0.030);
  k = 0;
  tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
  if (yLow == 0)
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
  else
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
  tx.DrawLatex(xtext, y0 + dy * k++, triggerName);
  tx.DrawLatex(xtext, y0 + dy * k++, Form("Entries = %lld", selected));

  const TString logPdfName = Form("%s/mass_draw_%s_log.pdf", figDir.Data(), figTag.Data());
  const TString logPngName = Form("%s/mass_draw_%s_log.png", figDir.Data(), figTag.Data());
  logCanvas.SaveAs(logPdfName);
  logCanvas.SaveAs(logPngName);

  TFile out(rootName, "RECREATE");
  hMass.Write("hMass");
  out.Close();

  printf("Saved:\n");
  printf("  %s\n", pdfName.Data());
  printf("  %s\n", pngName.Data());
  printf("  %s\n", logPdfName.Data());
  printf("  %s\n", logPngName.Data());
  printf("  %s\n", rootName.Data());
}
