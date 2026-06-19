#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooCrystalBall.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using namespace RooFit;

void mc_mass_binned(float ptLow = 3.5, float ptHigh = 50,
                    float yLow = 0.0, float yHigh = 2.4,
                    int bins = 1000,
                    int nSignalCBComponentsArg = 2,
                    int nSignalGaussComponentsArg = 2,
                    const char *triggerName = "HLT_HIL1SingleMu0_Open_v",
                    bool useDoubleSidedCB = false)
{
  gStyle->SetOptStat(0);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  const char *recoLabel = "McPrompt2023";
  const char *recoDisplayLabel = "MC Prompt 2023";
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
  mass.setBins(bins);
  RooDataHist binnedData("binnedData", "binned prompt MC mass", RooArgSet(mass), *data);
  const int nSignalCBComponents = std::clamp(nSignalCBComponentsArg, 0, 2);
  const int nSignalGaussComponents = std::clamp(nSignalGaussComponentsArg, 0, 3);
  if (nSignalCBComponents + nSignalGaussComponents <= 0) {
    printf("Invalid prompt MC mass model: no signal components enabled\n");
    return;
  }

  const double meanSeed = useDoubleSidedCB ? 3.097146 : 3.095;

  RooRealVar mean("mean", "Jpsi mean", meanSeed, 3.05, 3.15);
  const double sigmaG1Seed = useDoubleSidedCB ? 0.01929 : 0.030;
  const double sigmaG2Seed = useDoubleSidedCB ? 0.03679 : 0.095;
  const double sigmaG3Seed = useDoubleSidedCB ? 0.06203 : 0.150;
  RooRealVar sigmaG1("sigmaG1", "G1 sigma", sigmaG1Seed, 0.001, 0.080);
  RooRealVar sigmaG2Scale("sigmaG2Scale", "G2/G1 sigma scale",
                          sigmaG2Seed / sigmaG1Seed, 1.001, 20.0);
  RooFormulaVar sigmaG2("sigmaG2", "@0*@1", RooArgList(sigmaG1, sigmaG2Scale));
  RooRealVar sigmaG3Scale("sigmaG3Scale", "G3/G2 sigma scale",
                          sigmaG3Seed / sigmaG2Seed, 1.001, 20.0);
  RooFormulaVar sigmaG3("sigmaG3", "@0*@1", RooArgList(sigmaG2, sigmaG3Scale));

  const double sigmaDSCB1LSeed = useDoubleSidedCB ? 0.06336 : 0.045;
  const double sigmaDSCB1RSeed = useDoubleSidedCB ? 0.03755 : 0.060;
  const double sigmaDSCB2LSeed = 0.0275;
  const double sigmaDSCB2RSeed = 0.0780;
  RooRealVar sigmaCB1("sigmaCB1", useDoubleSidedCB ? "DSCB1 left sigma" : "CB1 sigma",
                      sigmaDSCB1LSeed, 0.001, 0.300);
  RooRealVar sigmaCB2Scale("sigmaCB2Scale", useDoubleSidedCB ? "DSCB1 right/left sigma scale" : "CB2/CB1 sigma scale",
                           sigmaDSCB1RSeed / sigmaDSCB1LSeed, 0.10, 10.0);
  RooFormulaVar sigmaCB2("sigmaCB2", "@0*@1", RooArgList(sigmaCB1, sigmaCB2Scale));
  RooRealVar sigmaDSCB2LScale("sigmaDSCB2LScale", "DSCB2 left/DSCB1 left sigma scale",
                              sigmaDSCB2LSeed / sigmaDSCB1LSeed, 0.10, 10.0);
  RooFormulaVar sigmaDSCB2L("sigmaDSCB2L", "@0*@1", RooArgList(sigmaCB1, sigmaDSCB2LScale));
  RooRealVar sigmaDSCB2RScale("sigmaDSCB2RScale", "DSCB2 right/DSCB1 right sigma scale",
                              sigmaDSCB2RSeed / sigmaDSCB1RSeed, 0.10, 10.0);
  RooFormulaVar sigmaDSCB2R("sigmaDSCB2R", "@0*@1", RooArgList(sigmaCB2, sigmaDSCB2RScale));

  RooGaussian sigG1("sigG1", "signal Gaussian 1", mass, mean, sigmaG1);
  RooGaussian sigG2("sigG2", "signal Gaussian 2", mass, mean, sigmaG2);
  RooGaussian sigG3("sigG3", "signal Gaussian 3", mass, mean, sigmaG3);

  RooRealVar alphaDSCBL("alphaDSCBL", useDoubleSidedCB ? "shared DSCB left alpha" : "CB left alpha",
                        useDoubleSidedCB ? 1.2 : 1.6, 0.20, 6.0);
  RooRealVar nDSCBL("nDSCBL", useDoubleSidedCB ? "shared DSCB left n" : "CB left n",
                    useDoubleSidedCB ? 2.0 : 3.0, 1.01, 80.0);
  RooRealVar alphaDSCBR("alphaDSCBR", useDoubleSidedCB ? "shared DSCB right alpha" : "CB right |alpha|",
                        useDoubleSidedCB ? 1.7 : 2.0, 0.20, 6.0);
  RooFormulaVar alphaCBRight("alphaCBRight", useDoubleSidedCB ? "@0" : "-@0",
                             RooArgList(alphaDSCBR));
  RooRealVar nDSCBR("nDSCBR", useDoubleSidedCB ? "shared DSCB right n" : "CB right n",
                    useDoubleSidedCB ? 2.0 : 4.0, 1.01, 80.0);
  RooRealVar alphaDSCB2LScale("alphaDSCB2LScale", "DSCB2 left alpha scale", 1.0, 0.25, 4.0);
  RooFormulaVar alphaDSCB2L("alphaDSCB2L", "@0*@1", RooArgList(alphaDSCBL, alphaDSCB2LScale));
  // RooRealVar nDSCB2LScale("nDSCB2LScale", "DSCB2 left n offset scale", 1.0, 0.05, 20.0);
  // RooFormulaVar nDSCB2L("nDSCB2L", "1.0+(@0-1.0)*@1", RooArgList(nDSCBL, nDSCB2LScale));
  RooRealVar nDSCB2L("nDSCB2L", "DSCB2 left n offset scale", 1.0, 0.05, 20.0);
  RooRealVar alphaDSCB2RScale("alphaDSCB2RScale", "DSCB2 right alpha scale", 1.0, 0.25, 4.0);
  RooFormulaVar alphaDSCB2R("alphaDSCB2R", "@0*@1", RooArgList(alphaDSCBR, alphaDSCB2RScale));
  // RooRealVar nDSCB2RScale("nDSCB2RScale", "DSCB2 right n offset scale", 1.0, 0.05, 20.0);
  // RooFormulaVar nDSCB2R("nDSCB2R", "1.0+(@0-1.0)*@1", RooArgList(nDSCBR, nDSCB2RScale));
  RooRealVar nDSCB2R("nDSCB2R", "DSCB2 right n", 1.0, 0.05, 20.0);
  RooCBShape sigCB1("sigCB1", "signal CB left", mass, mean, sigmaCB1, alphaDSCBL, nDSCBL);
  RooCBShape sigCB2("sigCB2", "signal CB right", mass, mean, sigmaCB2, alphaCBRight, nDSCBR); // shared n
  RooCrystalBall sigDSCB("sigDSCB", "signal double-sided Crystal Ball", mass, mean,
                         sigmaCB1, sigmaCB2, alphaDSCBL, nDSCBL, alphaDSCBR, nDSCBL);
  RooCrystalBall sigDSCB2("sigDSCB2", "signal double-sided Crystal Ball 2", mass, mean, sigmaDSCB2L, sigmaDSCB2R, alphaDSCB2L, nDSCB2L, alphaDSCB2R, nDSCB2L);

  RooRealVar frac1("frac1", "recursive fraction 1", useDoubleSidedCB ? 0.38 : 0.33,
                   1e-5, 1.0 - 1e-5);
  RooRealVar frac2("frac2", "recursive fraction 2", useDoubleSidedCB ? 0.24 : 0.45,
                   1e-5, 1.0 - 1e-5);
  RooRealVar frac3("frac3", "recursive fraction 3", useDoubleSidedCB ? 0.47 : 0.30,
                   1e-5, 1.0 - 1e-5);
  RooRealVar frac4("frac4", "recursive fraction 4", 0.50, 1e-5, 1.0 - 1e-5);

  RooArgList sigList;
  RooArgList sigFracList;
  if (useDoubleSidedCB) {
    if (nSignalCBComponents >= 1) sigList.add(sigDSCB);
    if (nSignalCBComponents >= 2) sigList.add(sigDSCB2);
    if (nSignalGaussComponents >= 1) sigList.add(sigG1);
    if (nSignalGaussComponents >= 2) sigList.add(sigG2);
    if (nSignalGaussComponents >= 3) sigList.add(sigG3);
  } else {
    if (nSignalGaussComponents >= 1) sigList.add(sigG1);
    if (nSignalCBComponents >= 1) sigList.add(sigCB1);
    if (nSignalCBComponents >= 2) sigList.add(sigCB2);
    if (nSignalGaussComponents >= 2) sigList.add(sigG2);
    if (nSignalGaussComponents >= 3) sigList.add(sigG3);
  }
  if (sigList.getSize() >= 2) sigFracList.add(frac1);
  if (sigList.getSize() >= 3) sigFracList.add(frac2);
  if (sigList.getSize() >= 4) sigFracList.add(frac3);
  if (sigList.getSize() >= 5) sigFracList.add(frac4);
  printf("Model: nDSCB/CB=%d, nGaussian=%d, nComponents=%d, useDoubleSidedCB=%d\n",
         nSignalCBComponents, nSignalGaussComponents, sigList.getSize(),
         useDoubleSidedCB ? 1 : 0);
  printf("Component order:");
  for (int i = 0; i < sigList.getSize(); ++i) {
    printf(" %s", sigList.at(i)->GetName());
  }
  printf("\n");

  std::unique_ptr<RooAddPdf> sigMassMulti;
  RooAbsPdf *sigMassPtr = dynamic_cast<RooAbsPdf *>(sigList.at(0));
  if (sigList.getSize() >= 2) {
    sigMassMulti = std::make_unique<RooAddPdf>(
        "sigMass", "configurable prompt MC signal mass", sigList, sigFracList, true);
    sigMassPtr = sigMassMulti.get();
  }
  RooAbsPdf &sigMass = *sigMassPtr;
  const double integratedBinningPrecision = 1e-6;
  const double plotPrecision = 1e-6;

  auto isGoodFit = [](const RooFitResult *fit) {
    return fit && fit->status() == 0 && fit->covQual() == 3;
  };
  auto printFitAttempt = [](const char *label, const RooFitResult *fit) {
    if (!fit) {
      printf("%s: null RooFitResult\n", label);
      return;
    }
    printf("%s: status=%d, covQual=%d, edm=%.6g, minNll=%.6g\n",
           label, fit->status(), fit->covQual(), fit->edm(), fit->minNll());
  };

  auto isBetterFit = [&](const RooFitResult *candidate, const RooFitResult *current) {
    if (!candidate) return false;
    if (!current) return true;
    if (isGoodFit(candidate) != isGoodFit(current)) return isGoodFit(candidate);
    if (candidate->covQual() != current->covQual()) return candidate->covQual() > current->covQual();
    const bool candidateStatusOk = candidate->status() == 0;
    const bool currentStatusOk = current->status() == 0;
    if (candidateStatusOk != currentStatusOk) return candidateStatusOk;
    return candidate->minNll() < current->minNll();
  };

  std::vector<std::unique_ptr<RooFitResult>> fitResults;
  std::vector<int> fitAttemptNumbers;
  fitResults.emplace_back(sigMass.fitTo(binnedData, Extended(false), Save(true), PrintLevel(-1), Strategy(1), Offset(true), NumCPU(12), Minimizer("Minuit2", "migrad"), IntegrateBins(integratedBinningPrecision), RecoverFromUndefinedRegions(1.0)));
  fitAttemptNumbers.push_back(1);
  printFitAttempt("Fit attempt 1 (Strategy 1)", fitResults.back().get());
  if (!isGoodFit(fitResults.back().get())) {
    printf("Refitting because attempt 1 did not satisfy status=0 and covQual=3\n");
    fitResults.emplace_back(sigMass.fitTo(binnedData, Extended(false), Save(true), PrintLevel(-1),
                                          Strategy(2), Offset(true), NumCPU(12), Minimizer("Minuit2", "migrad"),
                                          IntegrateBins(integratedBinningPrecision),
                                          RecoverFromUndefinedRegions(1.0)));
    fitAttemptNumbers.push_back(2);
    printFitAttempt("Fit attempt 2 (Strategy 2)", fitResults.back().get());
  }
  if (!isGoodFit(fitResults.back().get())) {
    printf("Refitting once more from attempt 2 parameters because final quality is still not status=0 and covQual=3\n");
    fitResults.emplace_back(sigMass.fitTo(binnedData, Extended(false), Save(true), PrintLevel(-1),
                                          Strategy(2), Offset(true), NumCPU(12), Minimizer("Minuit2", "migrad"),
                                          IntegrateBins(integratedBinningPrecision),
                                          RecoverFromUndefinedRegions(1.0)));
    fitAttemptNumbers.push_back(3);
    printFitAttempt("Fit attempt 3 (Strategy 2)", fitResults.back().get());
  }
  int bestFitIndex = -1;
  for (int i = 0; i < static_cast<int>(fitResults.size()); ++i) {
    if (isBetterFit(fitResults[i].get(),
                    bestFitIndex >= 0 ? fitResults[bestFitIndex].get() : nullptr)) {
      bestFitIndex = i;
    }
  }
  if (bestFitIndex < 0) {
    printf("Fit failed: no RooFitResult was produced\n");
    return;
  }
  if (bestFitIndex != static_cast<int>(fitResults.size()) - 1) {
    printf("Using fit attempt %d as final result because it has the best status/covQual/minNll\n",
           fitAttemptNumbers[bestFitIndex]);
    auto params = std::unique_ptr<RooArgSet>(sigMass.getParameters(binnedData));
    params->assignValueOnly(fitResults[bestFitIndex]->floatParsFinal());
  }
  const int fitAttempt = fitAttemptNumbers[bestFitIndex];
  auto result = std::move(fitResults[bestFitIndex]);
  if (!result) {
    printf("Fit failed: RooFitResult is null\n");
    return;
  }
  if (!isGoodFit(result.get())) {
    printf("WARNING: final fit still failed quality requirement status=0 and covQual=3\n");
  }
  printf("\nPrompt MC binned mass result:\n");
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
  binnedData.plotOn(frame.get(), Name("dataMass"));
  sigMass.plotOn(frame.get(), LineColor(kBlack), LineWidth(2),
                 Precision(plotPrecision), Name("fitMass"));
  if (nSignalGaussComponents >= 1) {
    sigMass.plotOn(frame.get(), Components(sigG1),
                   LineStyle(kDashed), LineColor(kBlue + 1), LineWidth(2),
                   Precision(plotPrecision), Name("sigG1"));
  }
  if (useDoubleSidedCB && nSignalCBComponents >= 1) {
    sigMass.plotOn(frame.get(), Components(sigDSCB),
                   LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(2),
                   Precision(plotPrecision), Name("sigDSCB"));
  } else if (nSignalCBComponents >= 1) {
    sigMass.plotOn(frame.get(), Components(sigCB1),
                   LineStyle(kDashed), LineColor(kGreen + 2), LineWidth(2),
                   Precision(plotPrecision), Name("sigCB1"));
  }
  if (useDoubleSidedCB && nSignalCBComponents >= 2) {
    sigMass.plotOn(frame.get(), Components(sigDSCB2),
                   LineStyle(kDashed), LineColor(kMagenta + 1), LineWidth(2),
                   Precision(plotPrecision), Name("sigDSCB2"));
  }
  if (!useDoubleSidedCB && nSignalCBComponents >= 2) {
    sigMass.plotOn(frame.get(), Components(sigCB2),
                   LineStyle(kDashed), LineColor(kMagenta + 1), LineWidth(2),
                   Precision(plotPrecision), Name("sigCB2"));
  }
  if (nSignalGaussComponents >= 2) {
    sigMass.plotOn(frame.get(), Components(sigG2),
                   LineStyle(kDashed), LineColor(kRed + 1), LineWidth(2),
                   Precision(plotPrecision), Name("sigG2"));
  }
  if (nSignalGaussComponents >= 3) {
    sigMass.plotOn(frame.get(), Components(sigG3),
                   LineStyle(kDashed), LineColor(kOrange + 7), LineWidth(2),
                   Precision(plotPrecision), Name("sigG3"));
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
  double maxAbsPull = 0.0;
  int nPullGt3 = 0;
  if (pull) {
    for (int i = 0; i < pull->GetN(); ++i) {
      double x = 0.0;
      double y = 0.0;
      pull->GetPoint(i, x, y);
      if (!std::isfinite(y)) continue;
      maxAbsPull = std::max(maxAbsPull, std::abs(y));
      if (std::abs(y) > 3.0) ++nPullGt3;
    }
  }
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
  const TString figDir = TString::Format("figs/%s/mc_mass_binned", yTag.Data());
  const TString rootDir = TString::Format("roots/%s/mc_mass_binned", yTag.Data());
  gSystem->mkdir(figDir, true);
  gSystem->mkdir(rootDir, true);

  TCanvas canvas("canvas", "prompt MC binned mass fit", 800, 800);
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
  TLegend legend(0.50, 0.66, 0.74, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.03);
  legend.AddEntry(frame->findObject("dataMass"), "Prompt MC", "lep");
  legend.AddEntry(frame->findObject("fitMass"), "Signal fit", "l");
  if (nSignalGaussComponents >= 1) legend.AddEntry(frame->findObject("sigG1"), "G1", "l");
  if (useDoubleSidedCB && nSignalCBComponents >= 1) legend.AddEntry(frame->findObject("sigDSCB"), "DSCB1", "l");
  if (useDoubleSidedCB && nSignalCBComponents >= 2) legend.AddEntry(frame->findObject("sigDSCB2"), "DSCB2", "l");
  if (!useDoubleSidedCB && nSignalCBComponents >= 1) legend.AddEntry(frame->findObject("sigCB1"), "CB left", "l");
  if (!useDoubleSidedCB && nSignalCBComponents >= 2) legend.AddEntry(frame->findObject("sigCB2"), "CB right", "l");
  if (nSignalGaussComponents >= 2) legend.AddEntry(frame->findObject("sigG2"), "G2", "l");
  if (nSignalGaussComponents >= 3) legend.AddEntry(frame->findObject("sigG3"), "G3", "l");
  legend.Draw("same");

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, Form("PbPb %s", recoDisplayLabel));
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.04);
    tx.SetTextFont(72);
    tx.DrawLatex(0.25, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    double xtext = 0.175;
    double y0 = 0.865;
    double dy = -0.05;
    int k = 0;
    tx.DrawLatex(xtext, y0 + dy * k++, "Prompt MC J/#psi #rightarrow #mu^{+}#mu^{-} (binned)");
    if (yLow == 0)
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f",
                                             ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f",
                                             ptLow, ptHigh, yLow, yHigh));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%s", triggerName));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("Attempt %d : status=%d, covQual=%d",
                                           fitAttempt, result->status(), result->covQual()));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("#chi^{2}/ndf = %.1f/%d = %.2f",
                                           chi.first, ndf, chi2OverNdf));
  }
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.024);
    tp.SetTextFont(42);
    double xtext = 0.74;
    double y0 = 0.87;
    double dy = -0.045;
    int k = 0;
    auto printVar = [&](const char *title, const RooAbsReal &var) {
      if (auto *rrv = dynamic_cast<const RooRealVar *>(&var)) {
        if (rrv->isConstant())
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
        else
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g",
                                                 title, rrv->getVal(), rrv->getError()));
      } else {
        const double err = result ? var.getPropagatedError(*result) : 0.0;
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g",
                                               title, var.getVal(), err));
      }
    };
    printVar("#mu", mean);
    if (nSignalGaussComponents >= 1) printVar("#sigma_{G1}", sigmaG1);
    if (nSignalCBComponents >= 1) {
      printVar(useDoubleSidedCB ? "#sigma_{1L}" : "#sigma_{CB1}", sigmaCB1);
      if (useDoubleSidedCB) {
        printVar("#sigma_{1R}", sigmaCB2);
        printVar("#alpha_{L}", alphaDSCBL);
        printVar("#alpha_{R}", alphaDSCBR);
        printVar("n_{L}", nDSCBL);
        printVar("n_{R}", nDSCBR);
      }
    }
    if (useDoubleSidedCB && nSignalCBComponents >= 2) {
      printVar("#sigma_{2L}", sigmaDSCB2L);
      printVar("#sigma_{2R}", sigmaDSCB2R);
      printVar("#alpha_{2L}", alphaDSCB2L);
      printVar("#alpha_{2R}", alphaDSCB2R);
    }
    if (!useDoubleSidedCB && nSignalCBComponents >= 2) {
      printVar("#sigma_{CB2}", sigmaCB2);
      printVar("#alpha_{CB2}", alphaCBRight);
      printVar("n_{CB2}", nDSCBR);
    }
    if (nSignalGaussComponents >= 2) printVar("#sigma_{G2}", sigmaG2);
    if (nSignalGaussComponents >= 3) printVar("#sigma_{G3}", sigmaG3);
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

  const TString pdfName = Form("%s/mc_mass_binned_fit_%s.pdf", figDir.Data(), figTag.Data());
  const TString pngName = Form("%s/mc_mass_binned_fit_%s.png", figDir.Data(), figTag.Data());
  const TString rootName = Form("%s/mc_mass_binned_model_%s.root", rootDir.Data(), figTag.Data());
  canvas.SaveAs(pdfName);
  canvas.SaveAs(pngName);

  TFile out(rootName, "RECREATE");
  result->Write("mcMassBinnedFitResult");
  sigMass.Write("mcMassBinnedModel");
  TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
  TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
  TParameter<int>("useDoubleSidedCB", useDoubleSidedCB ? 1 : 0).Write();
  TParameter<int>("massBins", bins).Write();
  TParameter<int>("mcMassBinnedFit", 1).Write();
  TParameter<double>("integratedBinningPrecision", integratedBinningPrecision).Write();
  TParameter<double>("plotPrecision", plotPrecision).Write();
  TParameter<double>("massChi2", chi.first).Write();
  TParameter<int>("massNdf", ndf).Write();
  TParameter<double>("massChi2OverNdf", chi2OverNdf).Write();
  TParameter<double>("massMaxAbsPull", maxAbsPull).Write();
  TParameter<int>("massNPullGt3", nPullGt3).Write();
  TParameter<int>("massFitStatus", result->status()).Write();
  TParameter<int>("massFitCovQual", result->covQual()).Write();
  TParameter<int>("massFitAttempt", fitAttempt).Write();
  TParameter<int>("massFitGood", isGoodFit(result.get()) ? 1 : 0).Write();
  TParameter<int>("mcMassBkgIncluded", 0).Write();
  auto writePar = [&](RooAbsReal &v) {
    double err = 0.0;
    if (auto *rrv = dynamic_cast<RooRealVar *>(&v)) {
      err = rrv->getError();
    } else if (result) {
      err = v.getPropagatedError(*result);
    }
    TParameter<double>(v.GetName(), v.getVal()).Write();
    TParameter<double>(Form("%sErr", v.GetName()), err).Write();
  };
  writePar(mean);
  if (nSignalGaussComponents >= 1) writePar(sigmaG1);
  if (useDoubleSidedCB) {
    if (nSignalCBComponents >= 1) {
      writePar(sigmaCB1);
      writePar(alphaDSCBL);
      writePar(nDSCBL);
      writePar(sigmaCB2);
      writePar(alphaDSCBR);
      writePar(nDSCBR);
    }
    if (nSignalCBComponents >= 2) {
      writePar(sigmaDSCB2L);
      writePar(alphaDSCB2L);
      writePar(nDSCB2L);
      writePar(sigmaDSCB2R);
      writePar(alphaDSCB2R);
      writePar(nDSCB2R);
      writePar(alphaDSCB2LScale);
      // writePar(nDSCB2LScale);
      writePar(alphaDSCB2RScale);
      // writePar(nDSCB2RScale);
    }
  } else {
    if (nSignalCBComponents >= 1) {
      writePar(sigmaCB1);
      writePar(alphaDSCBL);
      writePar(nDSCBL);
    }
    if (nSignalCBComponents >= 2) {
      writePar(sigmaCB2);
      writePar(alphaCBRight);
      writePar(nDSCBR);
    }
  }
  if (nSignalGaussComponents >= 2) writePar(sigmaG2);
  if (nSignalGaussComponents >= 3) writePar(sigmaG3);
  if (sigList.getSize() >= 2) writePar(frac1);
  if (sigList.getSize() >= 3) writePar(frac2);
  if (sigList.getSize() >= 4) writePar(frac3);
  if (sigList.getSize() >= 5) writePar(frac4);
  double remainingCoeff = 1.0;
  for (int i = 0; i < sigList.getSize(); ++i) {
    double coeff = remainingCoeff;
    if (i < sigFracList.getSize()) {
      auto *frac = dynamic_cast<RooAbsReal *>(sigFracList.at(i));
      const double fracVal = frac ? frac->getVal() : 0.0;
      coeff *= fracVal;
      remainingCoeff *= (1.0 - fracVal);
    }
    TParameter<double>(Form("componentCoeff%d", i + 1), coeff).Write();
  }
  if (nSignalGaussComponents >= 2) writePar(sigmaG2Scale);
  if (nSignalGaussComponents >= 3) writePar(sigmaG3Scale);
  if (useDoubleSidedCB || nSignalCBComponents >= 2) writePar(sigmaCB2Scale);
  if (useDoubleSidedCB && nSignalCBComponents >= 2) {
    writePar(sigmaDSCB2LScale);
    writePar(sigmaDSCB2RScale);
  }
  out.Close();

  printf("\nSaved:\n");
  printf("  %s\n", pdfName.Data());
  printf("  %s\n", pngName.Data());
  printf("  %s\n", rootName.Data());
}
