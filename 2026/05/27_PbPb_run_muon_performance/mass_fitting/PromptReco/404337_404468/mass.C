#include <memory>
#include <algorithm>
#include <cmath>
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
#include "RooBernstein.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using namespace RooFit;

void mass(float ptLow = 6.5, float ptHigh = 50,
                             float yLow = 0, float yHigh = 2.4,
                             int bkgOrder = 1,
                             int bins = 90,
                             int nSignalCBComponentsArg = 0,
                             int nSignalGaussComponentsArg = 2,
                             int nBkgExpComponentsArg = 0,
                             int nBkgChebyOrderArg = 3,
                             bool fixSignalShapeFromMc = false)
{
  // ------------------------------------------------------------------
  // run and trigger configuration
  // ------------------------------------------------------------------
  const char *pdName = "PromptReco";
  const char *runLabel = "404337_404474";
  const char *triggerName = "L1SingleMu0_Open";
  // 
  // Bit21: HLT_HIMinimumBiasHF1AND_v
  // Bit22: HLT_HIMinimumBiasHF1ANDZDC1nOR_v


  gStyle->SetOptStat(0);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // ------------------------------------------------------------------
  // load input dataset
  // ------------------------------------------------------------------
  // const char *dataRoot = "/data/users/pjgwak/work/daily_code_tracker/2026/05/27_PbPb_run_muon_performance/skimming/roodataset_roots/PromptReco/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_PromptReco_404337_404468_SingleMu0_Open.root";
  const char *dataRoot = "/data/users/pjgwak/work/daily_code_tracker/2026/05/27_PbPb_run_muon_performance/mass_fitting/PromptReco/404337_404468/roodataset_roots/PromptReco/TestRooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_PromptReco_404337_404474_SingleMu0_Open.root";
      
  auto inputFile = std::unique_ptr<TFile>(TFile::Open(dataRoot));
  if (!inputFile || inputFile->IsZombie()) {
    printf("Cannot open input file:\n  %s\n", dataRoot);
    return;
  }
  auto inputData = dynamic_cast<RooDataSet *>(inputFile->Get("dataset"));
  if (!inputData) {
    printf("Cannot find RooDataSet named 'dataset'\n");
    return;
  }

  // ------------------------------------------------------------------
  // select kinematic bin
  // ------------------------------------------------------------------
  TString cut = Form("mass > 2.6 && mass < 3.5 && "
                    //  "ctau3D > -6 && ctau3D < 8 && "
                     "pt > %g && pt < %g && abs(y) > %g && abs(y) < %g",
                     ptLow, ptHigh, yLow, yHigh);
  auto data = std::unique_ptr<RooDataSet>(
      static_cast<RooDataSet *>(inputData->reduce(cut)));
  if (!data || data->numEntries() <= 0) {
    printf("No entries after cut:\n  %s\n", cut.Data());
    return;
  }
  printf("Mass tune entries: %lld\n", (Long64_t)data->numEntries());
  printf("Cut: %s\n", cut.Data());

  RooRealVar mass("mass", "m_{#mu#mu}", 2.60, 3.50, "GeV");

  // ------------------------------------------------------------------
  // build output labels and read MC mass seeds
  // ------------------------------------------------------------------
  auto formatTag = [](double value) { return TString::Format("%.2f", value); };
  const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
  const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
  const TString figTag = yTag + "_" + ptTag;

  const TString mcMassName = Form("roots/%s/mc_mass/mc_mass_model_%s.root",
                                  yTag.Data(), figTag.Data());
  auto mcMassFile = std::unique_ptr<TFile>(TFile::Open(mcMassName));
  auto readMcMass = [&](const char *name, double fallback) {
    if (!mcMassFile || mcMassFile->IsZombie()) return fallback;
    auto *param = dynamic_cast<TParameter<double> *>(mcMassFile->Get(name));
    return param ? param->GetVal() : fallback;
  };
  auto readMcMassInt = [&](const char *name, int fallback) {
    if (!mcMassFile || mcMassFile->IsZombie()) return fallback;
    auto *param = dynamic_cast<TParameter<int> *>(mcMassFile->Get(name));
    return param ? param->GetVal() : fallback;
  };
  if (mcMassFile && !mcMassFile->IsZombie()) {
    printf("MC mass template: %s\n", mcMassName.Data());
  } else {
    printf("MC mass template not found, using defaults:\n  %s\n", mcMassName.Data());
  }
  auto clampValue = [](double value, double low, double high) {
    if (!std::isfinite(value)) return 0.5 * (low + high);
    return std::min(std::max(value, low), high);
  };
  auto rangeAround = [&](double value, double hardLow, double hardHigh,
                         double relLow, double relHigh,
                         double minLowGap, double minHighGap) {
    value = clampValue(value, hardLow + minLowGap, hardHigh - minHighGap);
    double low = value - std::max(minLowGap, std::abs(value) * relLow);
    double high = value + std::max(minHighGap, std::abs(value) * relHigh);
    low = std::max(hardLow, low);
    high = std::min(hardHigh, high);
    if (low >= value) low = std::max(hardLow, value - minLowGap);
    if (high <= value) high = std::min(hardHigh, value + minHighGap);
    return std::pair<double, double>(low, high);
  };
  auto initInRange = [&](double value, const std::pair<double, double> &range) {
    const double eps = std::max(1e-9, 1e-4 * (range.second - range.first));
    return clampValue(value, range.first + eps, range.second - eps);
  };

  // ------------------------------------------------------------------
  // configure signal and background components
  // ------------------------------------------------------------------
  const int nSignalCBComponents = std::clamp(
      nSignalCBComponentsArg >= 0
          ? nSignalCBComponentsArg
          : readMcMassInt("nSignalCBComponents", 2),
      0, 2);
  const int nSignalGaussComponents = std::clamp(
      nSignalGaussComponentsArg >= 0
          ? nSignalGaussComponentsArg
          : readMcMassInt("nSignalGaussComponents", 2),
      0, 2);
  if (nSignalCBComponents + nSignalGaussComponents <= 0) {
    printf("Invalid signal mass model: no signal components enabled\n");
    return;
  }
  const int nBkgExpComponents = std::clamp(nBkgExpComponentsArg, 0, 1);
  const int nBkgChebyOrder = std::clamp(
      nBkgChebyOrderArg >= 0 ? nBkgChebyOrderArg : bkgOrder, 0, 6);
  if (nBkgExpComponents > 0 && nBkgChebyOrder > 0) {
    printf("Invalid background mass model: choose exp or Cheby, not both\n");
    return;
  }

  // ------------------------------------------------------------------
  // build signal mass model
  // ------------------------------------------------------------------
  const double meanSeed = readMcMass("mean", 3.095);
  const auto meanRange = rangeAround(meanSeed, 3.05, 3.15, 0.003, 0.003, 0.010, 0.010);
  RooRealVar mean("mean", "J/#psi mean", initInRange(meanSeed, meanRange),
                  meanRange.first, meanRange.second);
  const double sigmaG1Seed = readMcMass("sigmaG1", 0.030);
  const auto sigmaG1Range = rangeAround(sigmaG1Seed, 0.001, 0.1, 0.75, 3.0, 0.008, 0.030);
  RooRealVar sigmaG1("sigmaG1", "G1 sigma", initInRange(sigmaG1Seed, sigmaG1Range),
                     sigmaG1Range.first, sigmaG1Range.second);
  const double sigmaG2Seed = readMcMass("sigmaG2", 0.095);
  const auto sigmaG2Range = rangeAround(sigmaG2Seed, 0.001, 0.15, 0.75, 3.0, 0.015, 0.060);
  RooRealVar sigmaG2("sigmaG2", "G2 sigma", initInRange(sigmaG2Seed, sigmaG2Range),
                     sigmaG2Range.first, sigmaG2Range.second);
  RooGaussian sigG1("sigG1", "signal Gaussian 1", mass, mean, sigmaG1);
  RooGaussian sigG2("sigG2", "signal Gaussian 2", mass, mean, sigmaG2);

  const double sigmaCB1Seed = readMcMass("sigmaCB1", 0.045);
  const auto sigmaCB1Range = rangeAround(sigmaCB1Seed, 0.001, 0.400, 0.75, 3.0, 0.010, 0.050);
  RooRealVar sigmaCB1("sigmaCB1", "CB left sigma",
                      initInRange(sigmaCB1Seed, sigmaCB1Range),
                      sigmaCB1Range.first, sigmaCB1Range.second);
  const double alphaCB1Seed = readMcMass("alphaCB1", 1.6);
  const auto alphaCB1Range = rangeAround(alphaCB1Seed, 0.01, 30.0, 0.75, 3.0, 0.20, 2.0);
  RooRealVar alphaCB1("alphaCB1", "CB left alpha",
                      initInRange(alphaCB1Seed, alphaCB1Range),
                      alphaCB1Range.first, alphaCB1Range.second);
  const double nCB1Seed = readMcMass("nCB1", 2.6);
  const auto nCB1Range = rangeAround(nCB1Seed, 0.01, 5.0, 0.75, 5.0, 0.50, 10.0);
  RooRealVar nCB1("nCB1", "CB left n", initInRange(nCB1Seed, nCB1Range),
                  nCB1Range.first, nCB1Range.second);
  RooCBShape sigCB1("sigCB1", "signal CB left", mass, mean,
                    sigmaCB1, alphaCB1, nCB1);

  const double sigmaCB2Seed = readMcMass("sigmaCB2", 0.060);
  const auto sigmaCB2Range = rangeAround(sigmaCB2Seed, 0.001, 0.500, 0.75, 3.0, 0.010, 0.060);
  RooRealVar sigmaCB2("sigmaCB2", "CB right sigma",
                      initInRange(sigmaCB2Seed, sigmaCB2Range),
                      sigmaCB2Range.first, sigmaCB2Range.second);
  const double alphaCB2Seed = readMcMass("alphaCB2", -2.0);
  const double alphaCB2Init = clampValue(alphaCB2Seed, -30.0 + 1e-3, -0.01 - 1e-3);
  const auto alphaCB2Range = std::pair<double, double>(
      std::max(-30.0, alphaCB2Init - std::max(0.20, std::abs(alphaCB2Init) * 3.0)),
      std::min(-0.01, alphaCB2Init + std::max(2.0, std::abs(alphaCB2Init) * 0.75)));
  RooRealVar alphaCB2("alphaCB2", "CB right alpha",
                      initInRange(alphaCB2Init, alphaCB2Range),
                      alphaCB2Range.first, alphaCB2Range.second);
  const double nCB2Seed = readMcMass("nCB2", 4.0);
  const auto nCB2Range = rangeAround(nCB2Seed, 0.01, 500.0, 0.75, 5.0, 0.50, 10.0);
  RooRealVar nCB2("nCB2", "CB right n", initInRange(nCB2Seed, nCB2Range),
                  nCB2Range.first, nCB2Range.second);
  RooCBShape sigCB2("sigCB2", "signal CB right", mass, mean,
                    sigmaCB2, alphaCB2, nCB2);

  RooRealVar fracG1("fracG1", "G1 fraction",
                    clampValue(readMcMass("fracG1", 0.25), 1e-4, 1.0 - 1e-4),
                    1e-4, 1.0 - 1e-4);
  RooRealVar fracCB1("fracCB1", "CB1 fraction",
                     clampValue(readMcMass("fracCB1", 0.45), 1e-4, 1.0 - 1e-4),
                     1e-4, 1.0 - 1e-4);
  RooRealVar fracCB2("fracCB2", "CB2 fraction",
                     clampValue(readMcMass("fracCB2", 0.20), 1e-4, 1.0 - 1e-4),
                     1e-4, 1.0 - 1e-4);
  if (fixSignalShapeFromMc && mcMassFile && !mcMassFile->IsZombie()) {
    for (auto *v : {&sigmaG2, &sigmaCB1, &alphaCB1, &nCB1,
                    &sigmaCB2, &alphaCB2, &nCB2, &fracG1, &fracCB1, &fracCB2}) {
      v->setConstant(true);
    }
  }
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
        "sigMass", "configurable signal mass", sigList, sigFracList, true);
    sigMassPtr = sigMassMulti.get();
  }
  RooAbsPdf &sigMass = *sigMassPtr;

  // ------------------------------------------------------------------
  // build background mass model and total PDF
  // ------------------------------------------------------------------
  RooRealVar bkgP0("bkgP0", "bkg p0", 0.0);
  bkgP0.setConstant(true);
  std::vector<std::unique_ptr<RooRealVar>> bkgPars;
  RooArgList bkgList;
  for (int i = 1; i <= nBkgChebyOrder; ++i) {
    bkgPars.push_back(std::make_unique<RooRealVar>(
        Form("bkgP%d", i), Form("bkg p%d", i), 0.0, -1.0, 1.0));
    bkgList.add(*bkgPars.back());
  }
  RooRealVar bkgLambda("bkgLambda", "background exponential slope",
                       -1.0, -20.0, -1e-4);
  std::unique_ptr<RooAbsPdf> bkgMassOwned;
  if (nBkgExpComponents >= 1) {
    bkgMassOwned = std::make_unique<RooExponential>(
        "bkgMass", "exponential background", mass, bkgLambda);
  } else {
    bkgMassOwned = std::make_unique<RooChebychev>(
        "bkgMass", "Chebychev background", mass, bkgList);
  }
  RooAbsPdf &bkgMass = *bkgMassOwned;

  RooRealVar nSig("nSig", "signal yield", 0.08 * data->numEntries(),
                  0.0, 0.8 * data->numEntries());
  RooRealVar nBkg("nBkg", "background yield", 0.92 * data->numEntries(),
                  0.0, 1.5 * data->numEntries());
  RooAddPdf model("massModel", "mass model", RooArgList(sigMass, bkgMass),
                  RooArgList(nSig, nBkg));

  // ------------------------------------------------------------------
  // fit mass model
  // ------------------------------------------------------------------
  RooArgSet externalConstraints;
  RooRealVar sigmaG1ConstraintMean("sigmaG1ConstraintMean",
                                   "sigmaG1 MC constraint mean",
                                   sigmaG1Seed);
  const double sigmaG1McErr = readMcMass("sigmaG1Err", 0.0);
  const double sigmaG1ConstraintWidth =
      std::max({sigmaG1McErr, 0.10 * std::abs(sigmaG1Seed), 0.001});
  RooRealVar sigmaG1ConstraintSigma("sigmaG1ConstraintSigma",
                                    "sigmaG1 MC constraint sigma",
                                    sigmaG1ConstraintWidth, 1e-9, 1.0);
  sigmaG1ConstraintMean.setConstant(true);
  sigmaG1ConstraintSigma.setConstant(true);
  RooGaussian sigmaG1Constraint("sigmaG1Constraint", "sigmaG1 MC constraint",
                                sigmaG1, sigmaG1ConstraintMean,
                                sigmaG1ConstraintSigma);
  if (fixSignalShapeFromMc && mcMassFile && !mcMassFile->IsZombie()) {
    externalConstraints.add(sigmaG1Constraint);
  }

  auto result = std::unique_ptr<RooFitResult>(
      model.fitTo(*data, Extended(true), Save(true), PrintLevel(-1),
                  Strategy(2), Offset(true), SumW2Error(false),
                  ExternalConstraints(externalConstraints)));
  printf("\nMass tune result:\n");
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
  const int nFloat = result->floatParsFinal().getSize();

  // ------------------------------------------------------------------
  // draw mass fit
  // ------------------------------------------------------------------
  auto frame = std::unique_ptr<RooPlot>(
      mass.frame(Bins(bins), Title("")));
  data->plotOn(frame.get(), Name("dataMass"));
  model.plotOn(frame.get(), LineColor(kBlack), LineWidth(2), Name("fitMass"));
  model.plotOn(frame.get(), Components(sigMass), LineStyle(kDashed),
               LineColor(kBlue + 1), LineWidth(2), Name("sigMass"));
  model.plotOn(frame.get(), Components(bkgMass), LineStyle(kDashed),
               LineColor(kRed + 1), LineWidth(2), Name("bkgMass"));
  double ymax = -1e300;
  if (auto *hdata = dynamic_cast<RooHist *>(frame->getHist("dataMass"))) {
    for (int i = 0; i < hdata->GetN(); ++i) {
      double xval = 0.0;
      double yval = 0.0;
      hdata->GetPoint(i, xval, yval);
      if (yval > ymax) ymax = yval;
    }
  }
  if (ymax > 0.0 && ymax < 1e300) frame->SetMaximum(ymax * 1.8);
  frame->SetMinimum(0.0);
  auto pull = frame->pullHist("dataMass", "fitMass");
  const auto chi = pullChi2(pull);
  const int ndf = std::max(1, chi.second - nFloat);

  auto pullFrame = std::unique_ptr<RooPlot>(mass.frame(Bins(bins), Title("")));
  pullFrame->addPlotable(pull, "P");
  pullFrame->SetMinimum(-5.0);
  pullFrame->SetMaximum(5.0);

  const TString figDir = TString::Format("figs/%s/mass", yTag.Data());
  const TString rootDir = TString::Format("roots/%s/mass", yTag.Data());
  gSystem->mkdir(figDir, true);
  gSystem->mkdir(rootDir, true);
  TCanvas canvas("canvas", "mass tuned fit", 800, 800);
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
  TLegend legend(0.50, 0.70, 0.74, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.03);
  legend.AddEntry(frame->findObject("dataMass"), "Data", "lep");
  legend.AddEntry(frame->findObject("fitMass"), "Fit", "l");
  legend.AddEntry(frame->findObject("sigMass"), "Signal", "l");
  legend.AddEntry(frame->findObject("bkgMass"), "Background", "l");
  legend.Draw("same");
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, Form("PbPb %s run %s", pdName, runLabel));
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.04);
    tx.SetTextFont(72);
    tx.DrawLatex(0.19, 0.935, "CMS Internal");
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
    tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
    if (yLow == 0)
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f",
                                             ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f",
                                             ptLow, ptHigh, yLow, yHigh));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%s", triggerName));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("Status : MINIMIZE=%d, covQual=%d",
                                           result->status(), result->covQual()));
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
    auto printVar = [&](const char *title, const RooRealVar &var) {
      if (var.isConstant())
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
      else
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g",
                                               title, var.getVal(), var.getError()));
    };
    printVar("N_{sig}", nSig);
    printVar("N_{bkg}", nBkg);
    printVar("#mu", mean);
    if (nSignalGaussComponents >= 1) printVar("#sigma_{G1}", sigmaG1);
    if (nSignalCBComponents >= 1) {
      printVar("#sigma_{CB1}", sigmaCB1);
      // printVar("#alpha_{CB1}", alphaCB1);
      // printVar("n_{CB1}", nCB1);
    }
    if (nSignalCBComponents >= 2) {
      printVar("#sigma_{CB2}", sigmaCB2);
      printVar("#alpha_{CB2}", alphaCB2);
      printVar("n_{CB2}", nCB2);
    }
    if (nSignalGaussComponents >= 2) printVar("#sigma_{G2}", sigmaG2);
  }

  // ------------------------------------------------------------------
  // draw pull and print summary
  // ------------------------------------------------------------------
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

  // ------------------------------------------------------------------
  // save figures and fit result
  // ------------------------------------------------------------------
  const TString pdfName = Form("%s/mass_fit_%s.pdf", figDir.Data(), figTag.Data());
  const TString pngName = Form("%s/mass_fit_%s.png", figDir.Data(), figTag.Data());
  const TString rootName = Form("%s/mass_model_%s.root", rootDir.Data(), figTag.Data());
  canvas.SaveAs(pdfName);
  canvas.SaveAs(pngName);

  TFile out(rootName, "RECREATE");
  result->Write("massFitResult");
  model.Write("massModel");
  TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
  TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
  TParameter<int>("nBkgExpComponents", nBkgExpComponents).Write();
  TParameter<int>("nBkgChebyOrder", nBkgChebyOrder).Write();
  TParameter<int>("bkgOrder", nBkgChebyOrder).Write();
  TParameter<int>("massBins", bins).Write();
  TParameter<double>("massChi2", chi.first).Write();
  TParameter<int>("massNdf", ndf).Write();
  TParameter<int>("massFitStatus", result->status()).Write();
  TParameter<int>("massFitCovQual", result->covQual()).Write();
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
  writePar(nSig);
  writePar(nBkg);
  TParameter<double>("bkgP0", bkgP0.getVal()).Write();
  if (nBkgExpComponents >= 1) writePar(bkgLambda);
  for (auto &p : bkgPars) writePar(*p);
  out.Close();

  printf("\nSaved:\n");
  printf("  %s\n", pdfName.Data());
  printf("  %s\n", pngName.Data());
  printf("  %s\n", rootName.Data());
}
