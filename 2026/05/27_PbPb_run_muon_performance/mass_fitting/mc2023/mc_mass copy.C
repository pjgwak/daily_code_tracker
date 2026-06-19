#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TParameter.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using namespace RooFit;

void mc_mass(float ptLow = 6.5, float ptHigh = 50,
          float yLow = 0.0, float yHigh = 2.4,
          int bins = 90,
          int nSignalCBComponentsArg = 2,
          int nSignalGaussComponentsArg = 2,
          const char *triggerName = "HLT_HIL1SingleMu0_Open_v")
{
  gStyle->SetOptStat(0);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  const char *recoLabel = "McPrompt2023";
  const TString mcRoot = TString::Format(
      "/data/users/pjgwak/work/daily_code_tracker/2026/05/27_PbPb_run_muon_performance/"
      "skimming/roodataset_roots/%s/"
      "RooDataSet_miniAOD_isMC1_PR_globalOff_Jpsi_cent0_200_"
      "Effw0_Accw0_PtW0_TnP0_%s_%s.root",
      recoLabel, recoLabel, triggerName);

  auto inputFile = std::unique_ptr<TFile>(TFile::Open(mcRoot));
  if (!inputFile || inputFile->IsZombie()) {
    printf("Cannot open input file:\n  %s\n", mcRoot.Data());
    return;
  }
  auto inputData = dynamic_cast<RooDataSet *>(inputFile->Get("dataset"));
  if (!inputData) {
    printf("Cannot find RooDataSet named 'dataset' in:\n  %s\n", mcRoot.Data());
    return;
  }

  TString cut = Form("mass > 2.6 && mass < 3.5 && "
                     "pt > %g && pt < %g && abs(y) > %g && abs(y) < %g",
                     ptLow, ptHigh, yLow, yHigh);
  auto data = std::unique_ptr<RooDataSet>(
      static_cast<RooDataSet *>(inputData->reduce(cut)));
  if (!data || data->numEntries() <= 0) {
    printf("No entries after cut:\n  %s\n", cut.Data());
    return;
  }
  printf("Prompt MC mass entries: %lld\n", (Long64_t)data->numEntries());
  printf("Input: %s\n", mcRoot.Data());
  printf("Cut: %s\n", cut.Data());

  RooRealVar mass("mass", "m_{#mu#mu}", 2.60, 3.50, "GeV");
  const int nSignalCBComponents = std::clamp(nSignalCBComponentsArg, 0, 2);
  const int nSignalGaussComponents = std::clamp(nSignalGaussComponentsArg, 0, 2);
  if (nSignalCBComponents + nSignalGaussComponents <= 0) {
    printf("Invalid prompt MC mass model: no signal components enabled\n");
    return;
  }

  const double meanSeed = 3.095;
  const double sigmaG1Seed = 0.030;
  const double sigmaG2Seed = 0.095;
  const double sigmaCB1Seed = 0.045;
  const double alphaCB1Seed = 1.6;
  const double nCB1Seed = 3.0;
  const double sigmaCB2Seed = 0.060;
  const double alphaCB2Seed = -2.0;
  const double nCB2Seed = 4.0;
  const double fracG1Seed = 0.25;
  const double fracCB1Seed = 0.45;
  const double fracCB2Seed = 0.20;

  RooRealVar mean("mean", "Jpsi mean", meanSeed, 3.05, 3.15);
  RooRealVar sigmaG1("sigmaG1", "G1 sigma", sigmaG1Seed, 0.001, 0.300);
  RooRealVar sigmaG2("sigmaG2", "G2 sigma", sigmaG2Seed, 0.001, 0.500);
  RooGaussian sigG1("sigG1", "signal Gaussian 1", mass, mean, sigmaG1);
  RooGaussian sigG2("sigG2", "signal Gaussian 2", mass, mean, sigmaG2);

  RooRealVar sigmaCB1("sigmaCB1", "CB left sigma", sigmaCB1Seed, 0.001, 0.400);
  RooRealVar alphaCB1("alphaCB1", "CB left alpha", alphaCB1Seed, 0.20, 8.0);
  RooRealVar nCB1("nCB1", "CB left n", nCB1Seed, 1.01, 80.0);
  RooCBShape sigCB1("sigCB1", "signal CB left", mass, mean, sigmaCB1, alphaCB1, nCB1);

  RooRealVar sigmaCB2("sigmaCB2", "CB right sigma", sigmaCB2Seed, 0.001, 0.20);
  RooRealVar alphaCB2("alphaCB2", "CB right alpha", alphaCB2Seed, -8.0, -0.20);
  RooRealVar nCB2("nCB2", "CB right n", nCB2Seed, 1.01, 80.0);
  RooCBShape sigCB2("sigCB2", "signal CB right", mass, mean, sigmaCB2, alphaCB2, nCB2);

  RooRealVar fracG1("fracG1", "G1 fraction", fracG1Seed, 1e-4, 1.0 - 1e-4);
  RooRealVar fracCB1("fracCB1", "CB1 fraction", fracCB1Seed, 1e-4, 1.0 - 1e-4);
  RooRealVar fracCB2("fracCB2", "CB2 fraction", fracCB2Seed, 1e-4, 1.0 - 1e-4);

  RooArgList sigList;
  RooArgList sigFracList;
  if (nSignalGaussComponents >= 1) sigList.add(sigG1);
  if (nSignalCBComponents >= 1) sigList.add(sigCB1);
  if (nSignalCBComponents >= 2) sigList.add(sigCB2);
  if (nSignalGaussComponents >= 2) sigList.add(sigG2);
  if (sigList.getSize() >= 2) sigFracList.add(fracG1);
  if (sigList.getSize() >= 3) sigFracList.add(fracCB1);
  if (sigList.getSize() >= 4) sigFracList.add(fracCB2);

  std::unique_ptr<RooAddPdf> sigMassMulti;
  RooAbsPdf *sigMassPtr = dynamic_cast<RooAbsPdf *>(sigList.at(0));
  if (sigList.getSize() >= 2) {
    sigMassMulti = std::make_unique<RooAddPdf>(
        "sigMass", "configurable prompt MC signal mass", sigList, sigFracList, true);
    sigMassPtr = sigMassMulti.get();
  }
  RooAbsPdf &sigMass = *sigMassPtr;

  auto result = std::unique_ptr<RooFitResult>(
      sigMass.fitTo(*data, Extended(false), Save(true), PrintLevel(-1),
                    Strategy(1), Offset(true), NumCPU(12), RecoverFromUndefinedRegions(1.0)));
  if (!result || result->status() != 0 || result->covQual() < 2) {
    result.reset(sigMass.fitTo(*data, Extended(false), Save(true), PrintLevel(-1),
                               Strategy(2), Offset(true), NumCPU(12), RecoverFromUndefinedRegions(1.0)));
  }
  if (!result) {
    printf("Fit failed: RooFitResult is null\n");
    return;
  }
  printf("\nPrompt MC mass result:\n");
  result->Print("v");

  auto pullChi2 = [](const RooHist *pull) {
    double chi2 = 0.0;
    int n = 0;
    if (!pull) return std::pair<double, int>(chi2, n);
    for (int i = 0; i < pull->GetN(); ++i) {
      double x = 0.0;
      double y = 0.0;
      pull->GetPoint(i, x, y);
      if (!std::isfinite(y)) continue;
      chi2 += y * y;
      ++n;
    }
    return std::pair<double, int>(chi2, n);
  };

  auto frame = std::unique_ptr<RooPlot>(mass.frame(Bins(bins), Title("")));
  data->plotOn(frame.get(), Name("dataMass"));
  sigMass.plotOn(frame.get(), LineColor(kBlack), LineWidth(2), Name("fitMass"));
  if (nSignalGaussComponents >= 1) {
    sigMass.plotOn(frame.get(), Components(sigG1), LineStyle(kDashed), LineColor(kBlue + 1), Name("sigG1"));
  }
  if (nSignalCBComponents >= 1) {
    sigMass.plotOn(frame.get(), Components(sigCB1), LineStyle(kDashed), LineColor(kGreen + 2), Name("sigCB1"));
  }
  if (nSignalCBComponents >= 2) {
    sigMass.plotOn(frame.get(), Components(sigCB2), LineStyle(kDashed), LineColor(kMagenta + 1), Name("sigCB2"));
  }
  if (nSignalGaussComponents >= 2) {
    sigMass.plotOn(frame.get(), Components(sigG2), LineStyle(kDashed), LineColor(kRed + 1), Name("sigG2"));
  }
  double ymax = -1e300;
  if (auto *hdata = dynamic_cast<RooHist *>(frame->getHist("dataMass"))) {
    for (int i = 0; i < hdata->GetN(); ++i) {
      double xval = 0.0;
      double yval = 0.0;
      hdata->GetPoint(i, xval, yval);
      ymax = std::max(ymax, yval);
    }
  }
  if (ymax > 0.0 && ymax < 1e300) frame->SetMaximum(ymax * 1.8);
  frame->SetMinimum(0.0);

  auto pull = frame->pullHist("dataMass", "fitMass");
  const auto chi = pullChi2(pull);
  const int nFloat = result->floatParsFinal().getSize();
  const int ndf = std::max(1, chi.second - nFloat);
  const double chi2OverNdf = chi.first / ndf;
  auto pullFrame = std::unique_ptr<RooPlot>(mass.frame(Bins(bins), Title("")));
  pullFrame->addPlotable(pull, "P");
  pullFrame->SetMinimum(-5.0);
  pullFrame->SetMaximum(5.0);

  auto formatTag = [](double value) { return TString::Format("%.2f", value); };
  const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
  const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
  const TString trigTag(triggerName);
  const TString figTag = yTag + "_" + ptTag + "_" + trigTag;
  const TString figDir = TString::Format("figs/%s/mc_mass", yTag.Data());
  const TString rootDir = TString::Format("roots/%s/mc_mass", yTag.Data());
  gSystem->mkdir(figDir, true);
  gSystem->mkdir(rootDir, true);

  TCanvas canvas("canvas", "prompt MC mass fit", 800, 800);
  TPad top("top", "top", 0.0, 0.25, 1.0, 1.0);
  TPad bottom("bottom", "bottom", 0.0, 0.0, 1.0, 0.25);
  top.SetBottomMargin(0.00001);
  top.SetTopMargin(0.08);
  top.SetLeftMargin(0.13);
  bottom.SetTopMargin(0.00001);
  bottom.SetBottomMargin(0.40);
  bottom.SetLeftMargin(0.13);
  top.Draw();
  bottom.Draw();

  top.cd();
  frame->GetYaxis()->SetTitle("Events");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.042);
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw("e");
  TLegend legend(0.52, 0.66, 0.78, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.03);
  legend.AddEntry(frame->findObject("dataMass"), "Prompt MC", "lep");
  legend.AddEntry(frame->findObject("fitMass"), "Signal fit", "l");
  if (nSignalGaussComponents >= 1) legend.AddEntry(frame->findObject("sigG1"), "G1", "l");
  if (nSignalCBComponents >= 1) legend.AddEntry(frame->findObject("sigCB1"), "CB left", "l");
  if (nSignalCBComponents >= 2) legend.AddEntry(frame->findObject("sigCB2"), "CB right", "l");
  if (nSignalGaussComponents >= 2) legend.AddEntry(frame->findObject("sigG2"), "G2", "l");
  legend.Draw("same");

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.18, 0.86, "Prompt MC J/#psi #rightarrow #mu^{+}#mu^{-}");
  if (yLow == 0)
    tx.DrawLatex(0.18, 0.81, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
  else
    tx.DrawLatex(0.18, 0.81, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
  tx.DrawLatex(0.18, 0.76, triggerName);
  tx.DrawLatex(0.18, 0.71, Form("status/cov = %d/%d", result->status(), result->covQual()));
  tx.DrawLatex(0.18, 0.66, Form("#chi^{2}/ndf = %.1f/%d = %.2f",
                                 chi.first, ndf, chi2OverNdf));

  TLatex tp;
  tp.SetNDC();
  tp.SetTextSize(0.024);
  tp.SetTextFont(42);
  double xtext = 0.78;
  double y0 = 0.62;
  double dy = -0.042;
  int k = 0;
  auto printVar = [&](const char *title, const RooRealVar &var) {
    tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), var.getError()));
  };
  printVar("#mu", mean);
  if (nSignalGaussComponents >= 1) printVar("#sigma_{G1}", sigmaG1);
  if (nSignalCBComponents >= 1) printVar("#sigma_{CB1}", sigmaCB1);
  if (nSignalCBComponents >= 2) {
    printVar("#sigma_{CB2}", sigmaCB2);
    printVar("#alpha_{CB2}", alphaCB2);
    printVar("n_{CB2}", nCB2);
  }

  bottom.cd();
  pullFrame->GetYaxis()->SetTitle("Pull");
  pullFrame->GetYaxis()->CenterTitle();
  pullFrame->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
  pullFrame->GetXaxis()->CenterTitle();
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetYaxis()->SetTitleSize(0.12);
  pullFrame->GetYaxis()->SetTitleOffset(0.40);
  pullFrame->GetYaxis()->SetLabelSize(0.10);
  pullFrame->GetXaxis()->SetTitleSize(0.14);
  pullFrame->GetXaxis()->SetTitleOffset(1.00);
  pullFrame->GetXaxis()->SetLabelSize(0.10);
  pullFrame->GetXaxis()->SetLabelOffset(0.02);
  pullFrame->Draw();
  TLine zero(mass.getMin(), 0.0, mass.getMax(), 0.0);
  TLine plus3(mass.getMin(), 3.0, mass.getMax(), 3.0);
  TLine minus3(mass.getMin(), -3.0, mass.getMax(), -3.0);
  plus3.SetLineStyle(kDashed);
  minus3.SetLineStyle(kDashed);
  plus3.SetLineColor(kRed + 1);
  minus3.SetLineColor(kRed + 1);
  zero.Draw("same");
  plus3.Draw("same");
  minus3.Draw("same");

  const TString pdfName = Form("%s/mc_mass_fit_%s.pdf", figDir.Data(), figTag.Data());
  const TString pngName = Form("%s/mc_mass_fit_%s.png", figDir.Data(), figTag.Data());
  const TString rootName = Form("%s/mc_mass_model_%s.root", rootDir.Data(), figTag.Data());
  canvas.SaveAs(pdfName);
  canvas.SaveAs(pngName);

  TFile out(rootName, "RECREATE");
  result->Write("mcMassFitResult");
  sigMass.Write("mcMassModel");
  TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
  TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
  TParameter<int>("massBins", bins).Write();
  TParameter<double>("massChi2", chi.first).Write();
  TParameter<int>("massNdf", ndf).Write();
  TParameter<double>("massChi2OverNdf", chi2OverNdf).Write();
  TParameter<int>("massFitStatus", result->status()).Write();
  TParameter<int>("massFitCovQual", result->covQual()).Write();
  TParameter<int>("mcMassBkgIncluded", 0).Write();
  auto writePar = [&](RooRealVar &v) {
    TParameter<double>(v.GetName(), v.getVal()).Write();
    TParameter<double>(Form("%sErr", v.GetName()), v.getError()).Write();
  };
  writePar(mean);
  writePar(sigmaG1);
  writePar(sigmaG2);
  writePar(sigmaCB1);
  writePar(alphaCB1);
  writePar(nCB1);
  writePar(sigmaCB2);
  writePar(alphaCB2);
  writePar(nCB2);
  writePar(fracG1);
  writePar(fracCB1);
  writePar(fracCB2);
  out.Close();

  printf("\nSaved:\n");
  printf("  %s\n", pdfName.Data());
  printf("  %s\n", pngName.Data());
  printf("  %s\n", rootName.Data());
}
