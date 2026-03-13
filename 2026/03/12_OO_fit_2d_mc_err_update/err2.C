#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooLandau.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "TH1.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"
#include "RooHist.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace RooFit;

static std::pair<double, int> err2_chi2_from_pull(RooHist *hpull)
{
  double chi2 = 0.0;
  int n = 0;
  if (!hpull)
    return {0.0, 0};
  double x = 0.0, y = 0.0;
  for (int i = 0; i < hpull->GetN(); ++i)
  {
    hpull->GetPoint(i, x, y);
    if (!std::isfinite(y))
      continue;
    chi2 += y * y;
    ++n;
  }
  return {chi2, n};
}

static void apply_logy_auto_range(RooPlot *plot, const char *histName, double topScale = 20.0, double bottomScale = 1e3)
{
  if (!plot)
    return;
  auto *hist = dynamic_cast<RooHist *>(plot->getHist(histName));
  if (!hist)
    return;

  double peak = 0.0;
  double minPositive = std::numeric_limits<double>::infinity();
  for (int i = 0; i < hist->GetN(); ++i)
  {
    double x = 0.0, y = 0.0;
    hist->GetPoint(i, x, y);
    if (!std::isfinite(y) || y <= 0.0)
      continue;
    peak = std::max(peak, y);
    minPositive = std::min(minPositive, y);
  }
  if (peak <= 0.0)
    return;

  double ymin = peak / bottomScale;
  if (std::isfinite(minPositive))
    ymin = std::min(ymin, 0.5 * minPositive);
  ymin = std::max(ymin, 1e-3);
  const double ymax = std::max(peak * topScale, ymin * 10.0);
  plot->SetMinimum(ymin);
  plot->SetMaximum(ymax);
}

static void err2_report_params_near_limit(const RooFitResult *fr, const char *tag, double relTol = 0.03)
{
  if (!fr || relTol <= 0.0)
    return;
  const RooArgList &pars = fr->floatParsFinal();
  for (int i = 0; i < pars.getSize(); ++i)
  {
    auto *rrv = dynamic_cast<RooRealVar *>(pars.at(i));
    if (!rrv)
      continue;
    const double lo = rrv->getMin();
    const double hi = rrv->getMax();
    if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo))
      continue;
    const double span = hi - lo;
    const double v = rrv->getVal();
    const double dLow = (v - lo) / span;
    const double dHigh = (hi - v) / span;
    if (dLow < relTol || dHigh < relTol)
    {
      printf("[RangeCheck][%s] %s = %.6g near bound [%.6g, %.6g]\n",
             tag ? tag : "fit", rrv->GetName(), v, lo, hi);
    }
  }
}

void err2(float ptLow = 14, float ptHigh = 20, float yLow = 1.6, float yHigh = 2.4)
{
  bool isWeight = false;

  // ------------------------------------------------------------------
  // model control
  // ------------------------------------------------------------------
  enum ErrPdfChoice
  {
    kErrPdfAnalytic = 0,
    kErrPdfHist = 1,
  };
  // Choose how many ctau-error components to use in each (pt, y) bin.
  int bkgErrPdfOpt = kErrPdfHist;
  int sigErrPdfOpt = kErrPdfAnalytic; // kErrPdfHist, kErrPdfAnalytic
  const int histPdfInterpolationOrder = 1;
  int nBkgTimeErrGaussComponents = 2;
  int nBkgTimeErrLandauComponents = 1;
  int nBkgTimeErrLognormalComponents = 1;
  int nSigTimeErrGaussComponents = 2;
  int nSigTimeErrLandauComponents = 1;
  int nSigTimeErrLognormalComponents = 1;
  bool useStagedTimeErrFitBkg = true;
  bool useStagedTimeErrFitSig = true;

  if (yLow == 1.6f)
  {
    if (ptLow == 1.0f && ptHigh == 2.0f)
    {
      nBkgTimeErrGaussComponents = 1;     // 0~2
      nBkgTimeErrLandauComponents = 2;    // 0~2
      nBkgTimeErrLognormalComponents = 1; // 0~1

      nSigTimeErrGaussComponents = 1;     // 0~2
      nSigTimeErrLandauComponents = 0;    // 0~2
      nSigTimeErrLognormalComponents = 1; // 0~1
    }
  } else
  {
    if (ptLow == 7.0f && ptHigh == 8.0f)
    {
      nBkgTimeErrGaussComponents = 1;     // 0~2
      nBkgTimeErrLandauComponents = 1;    // 0~2
      nBkgTimeErrLognormalComponents = 1; // 0~1

      nSigTimeErrGaussComponents = 1;     // 0~2
      nSigTimeErrLandauComponents = 1;    // 0~2
      nSigTimeErrLognormalComponents = 1; // 0~1
    }
  }

  nBkgTimeErrGaussComponents = std::clamp(nBkgTimeErrGaussComponents, 0, 2);
  nBkgTimeErrLandauComponents = std::clamp(nBkgTimeErrLandauComponents, 0, 2);
  nBkgTimeErrLognormalComponents = std::clamp(nBkgTimeErrLognormalComponents, 0, 1);
  nSigTimeErrGaussComponents = std::clamp(nSigTimeErrGaussComponents, 0, 2);
  nSigTimeErrLandauComponents = std::clamp(nSigTimeErrLandauComponents, 0, 2);
  nSigTimeErrLognormalComponents = std::clamp(nSigTimeErrLognormalComponents, 0, 1);

  if (nBkgTimeErrGaussComponents + nBkgTimeErrLandauComponents + nBkgTimeErrLognormalComponents <= 0)
  {
    std::cerr << "ERROR: at least one background timeErr component is required." << std::endl;
    return;
  }
  if (nSigTimeErrGaussComponents + nSigTimeErrLandauComponents + nSigTimeErrLognormalComponents <= 0)
  {
    std::cerr << "ERROR: at least one signal timeErr component is required." << std::endl;
    return;
  }
  const bool bkgUseGaus1 = (nBkgTimeErrGaussComponents >= 1);
  const bool bkgUseGaus2 = (nBkgTimeErrGaussComponents >= 2);
  const bool bkgUseLandau1 = (nBkgTimeErrLandauComponents >= 1);
  const bool bkgUseLandau2 = (nBkgTimeErrLandauComponents >= 2);
  const bool bkgUseLognormal = (nBkgTimeErrLognormalComponents >= 1);
  const int nBkgTimeErrComponents =
      static_cast<int>(bkgUseGaus1) + static_cast<int>(bkgUseGaus2) + static_cast<int>(bkgUseLandau1) +
      static_cast<int>(bkgUseLandau2) + static_cast<int>(bkgUseLognormal);
  const bool sigUseGaus1 = (nSigTimeErrGaussComponents >= 1);
  const bool sigUseGaus2 = (nSigTimeErrGaussComponents >= 2);
  const bool sigUseLandau1 = (nSigTimeErrLandauComponents >= 1);
  const bool sigUseLandau2 = (nSigTimeErrLandauComponents >= 2);
  const bool sigUseLognormal = (nSigTimeErrLognormalComponents >= 1);
  const int nSigTimeErrComponents =
      static_cast<int>(sigUseGaus1) + static_cast<int>(sigUseGaus2) + static_cast<int>(sigUseLandau1) +
      static_cast<int>(sigUseLandau2) + static_cast<int>(sigUseLognormal);

  // ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
  const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
  const char *DSET_NAME = "dataset";

  TFile *inputFile = TFile::Open(DATA_ROOT);
  if (!inputFile || inputFile->IsZombie())
  {
    std::cerr << "ERROR: cannot open input data file: " << DATA_ROOT << std::endl;
    return;
  }

  RooDataSet *inputData = dynamic_cast<RooDataSet *>(inputFile->Get(DSET_NAME));
  if (!inputData)
  {
    std::cerr << "ERROR: RooDataSet '" << DSET_NAME << "' not found in " << DATA_ROOT << std::endl;
    return;
  }

  TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)",
                          ptLow, ptHigh, yLow, yHigh);
  auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
  if (!dataBasic || dataBasic->numEntries() <= 0)
  {
    std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
    return;
  }

  auto *massTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
  auto *timeTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
  auto *timeErrTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DErr"));
  if (!massTmp || !timeTmp || !timeErrTmp)
  {
    std::cerr << "ERROR: required variables mass/ctau3D/ctau3DErr are missing in dataset." << std::endl;
    return;
  }

  auto quantileRange = [&](RooDataSet &ds, RooRealVar &var, double qLo, double qHi, bool positiveOnly)
  {
    std::vector<double> vals;
    vals.reserve(ds.numEntries());
    for (int i = 0; i < ds.numEntries(); ++i)
    {
      ds.get(i);
      const double v = var.getVal();
      if (!std::isfinite(v))
        continue;
      if (positiveOnly && v <= 0.0)
        continue;
      vals.push_back(v);
    }
    if (vals.empty())
      return std::make_pair(var.getMin(), var.getMax());
    std::sort(vals.begin(), vals.end());
    auto qAt = [&](double q)
    {
      long long idx = llround(q * (vals.size() - 1));
      if (idx < 0)
        idx = 0;
      if (idx >= static_cast<long long>(vals.size()))
        idx = static_cast<long long>(vals.size()) - 1;
      return vals[static_cast<size_t>(idx)];
    };
    return std::make_pair(qAt(qLo), qAt(qHi));
  };

  const auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.001, 0.999, false);
  auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.995, true);
  if (errRange.first < 1e-6)
    errRange.first = 1e-6;

  TString cutAll = Form(
      "(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
      "(ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)",
      ptLow, ptHigh, yLow, yHigh, ctRange.first, ctRange.second, errRange.first, errRange.second);
  auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutAll)));
  if (!dataSel || dataSel->numEntries() <= 0)
  {
    std::cerr << "ERROR: no entries after final selection: " << cutAll << std::endl;
    return;
  }

  const double sidebandLeftMax = 2.9;
  const double sidebandRightMin = 3.2;
  TString cutSideband = Form("(mass >= 2.6 && mass < %g) || (mass > %g && mass <= 3.5)",
                             sidebandLeftMax, sidebandRightMin);
  auto dataSB = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSel->reduce(cutSideband)));
  if (!dataSB || dataSB->numEntries() <= 0)
  {
    std::cerr << "ERROR: no entries after sideband selection: " << cutSideband << std::endl;
    return;
  }

  auto formatTag = [](double value) { return TString::Format("%.2f", value); };
  const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
  const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
  const TString figTag = yTag + "_" + ptTag;
  const TString figDir = TString::Format("figs/%s/err2", yTag.Data());
  const TString resultDir = TString::Format("roots/%s/err2", yTag.Data());
  const TString massFileName = TString::Format("roots/%s/mass/mass_model_%s.root", yTag.Data(), figTag.Data());
  auto figName = [&](const char *name)
  {
    return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
  };
  const TString resultFileName = TString::Format("%s/err2_model_%s.root", resultDir.Data(), figTag.Data());

  gSystem->mkdir(figDir, true);
  gSystem->mkdir(resultDir, true);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
  RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
  RooRealVar &obs_timeErr = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
  obs_mass.SetTitle("mass");
  obs_mass.setUnit("GeV/c^{2}");
  obs_time.SetTitle("#font[12]{l}_{J/#psi}");
  obs_time.setUnit("mm");
  obs_timeErr.SetTitle("event-by-event ctau error");
  obs_timeErr.setUnit("mm");
  obs_time.setRange(ctRange.first, ctRange.second);
  obs_timeErr.setRange(errRange.first, errRange.second);
  obs_timeErr.setMin(errRange.first);
  const int timeErrPlotBins = std::max(2, obs_timeErr.getBins());

  RooDataSet *data = dataSB.get();
  auto timeErrData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_timeErr))));
  auto dataSR = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSel->reduce(Form("(mass >= %g && mass <= %g)", sidebandLeftMax, sidebandRightMin))));

  const double errSpan = errRange.second - errRange.first;
  const double gaus1MeanInit = errRange.first + 0.14 * errSpan;
  const double gaus2MeanInit = errRange.first + 0.24 * errSpan;
  const double tailMpvInit = errRange.first + 0.40 * errSpan;
  const double gaus1SigmaInit = std::max(0.025 * errSpan, 5e-4);
  const double gaus2SigmaInit = std::max(0.055 * errSpan, 1e-3);
  const double tailWidthInit = std::max(0.10 * errSpan, 1e-3);
  const double tail2MpvInit = errRange.first + 0.62 * errSpan;
  const double tail2WidthInit = std::max(0.18 * errSpan, 2e-3);
  const double lognM0Init = errRange.first + 0.72 * errSpan;
  const double lognKInit = 0.35;
  const double errMeanMin = std::max(1e-6, errRange.first - 0.35 * errSpan);
  const double errMeanMax = errRange.second + 0.80 * errSpan;
  const double gaus1SigmaMax = std::max(0.40 * errSpan, 6e-3);
  const double gaus2SigmaMax = std::max(0.70 * errSpan, 1e-2);
  const double tailWidthMax = std::max(1.20 * errSpan, 2e-2);
  const double tail2WidthMax = std::max(1.80 * errSpan, 3e-2);
  const auto fracToRatio = [](double frac)
  {
    const double f = std::clamp(frac, 1e-6, 1.0 - 1e-6);
    return f / (1.0 - f);
  };

  RooRealVar timeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", gaus1MeanInit, errMeanMin, errMeanMax);
  RooRealVar timeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", gaus1SigmaInit, 5e-4, gaus1SigmaMax);
  RooRealVar timeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", gaus2MeanInit, errMeanMin, errMeanMax);
  RooRealVar timeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", gaus2SigmaInit, 1e-3, gaus2SigmaMax);
  RooRealVar timeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", tailMpvInit, errMeanMin, errMeanMax);
  RooRealVar timeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", tailWidthInit, 1e-3, tailWidthMax);
  RooRealVar timeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", tail2MpvInit, errMeanMin, errMeanMax);
  RooRealVar timeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", tail2WidthInit, 1e-3, tail2WidthMax);
  RooRealVar timeErrLognM0("timeErrLognM0", "timeErrLognM0", lognM0Init, errMeanMin, errMeanMax);
  RooRealVar timeErrLognK("timeErrLognK", "timeErrLognK", lognKInit, 0.01, 3.0);
  RooRealVar timeErrCore1FracRatio("timeErrCore1FracRatio", "timeErrCore1FracRatio", fracToRatio(0.60), 1e-3, 1e3);
  RooFormulaVar timeErrCore1Frac("timeErrCore1Frac", "@0/(1.0+@0)", RooArgList(timeErrCore1FracRatio));
  RooRealVar timeErrTailFracRatio("timeErrTailFracRatio", "timeErrTailFracRatio", fracToRatio(0.08), 1e-4, 1e4);
  RooFormulaVar timeErrTailFrac("timeErrTailFrac", "@0/(1.0+@0)", RooArgList(timeErrTailFracRatio));
  RooRealVar timeErrTail2FracRatio("timeErrTail2FracRatio", "timeErrTail2FracRatio", fracToRatio(0.05), 1e-4, 1e4);
  RooFormulaVar timeErrTail2Frac("timeErrTail2Frac", "@0/(1.0+@0)", RooArgList(timeErrTail2FracRatio));
  RooRealVar timeErrLognFracRatio("timeErrLognFracRatio", "timeErrLognFracRatio", fracToRatio(0.03), 1e-4, 1e4);
  RooFormulaVar timeErrLognFrac("timeErrLognFrac", "@0/(1.0+@0)", RooArgList(timeErrLognFracRatio));

  RooGaussian timeErrGaus1("timeErrGaus1", "timeErrGaus1", obs_timeErr, timeErrGaus1Mean, timeErrGaus1Sigma);
  RooGaussian timeErrGaus2("timeErrGaus2", "timeErrGaus2", obs_timeErr, timeErrGaus2Mean, timeErrGaus2Sigma);
  RooLandau timeErrTail("timeErrTail", "timeErrTail", obs_timeErr, timeErrTailMpv, timeErrTailWidth);
  RooLandau timeErrTail2("timeErrTail2", "timeErrTail2", obs_timeErr, timeErrTail2Mpv, timeErrTail2Width);
  RooLognormal timeErrLogn("timeErrLogn", "timeErrLogn", obs_timeErr, timeErrLognM0, timeErrLognK);

  RooArgList timeErrPdfList;
  RooArgList timeErrFracList;
  if (bkgUseLandau1)
  {
    timeErrPdfList.add(timeErrTail);
    if (nBkgTimeErrComponents > 1)
      timeErrFracList.add(timeErrTailFrac);
  }
  if (bkgUseLandau2)
  {
    timeErrPdfList.add(timeErrTail2);
    if (static_cast<int>(timeErrPdfList.getSize()) < nBkgTimeErrComponents)
      timeErrFracList.add(timeErrTail2Frac);
  }
  if (bkgUseLognormal)
  {
    timeErrPdfList.add(timeErrLogn);
    if (static_cast<int>(timeErrPdfList.getSize()) < nBkgTimeErrComponents)
      timeErrFracList.add(timeErrLognFrac);
  }
  if (bkgUseGaus1)
  {
    timeErrPdfList.add(timeErrGaus1);
    if (static_cast<int>(timeErrPdfList.getSize()) < nBkgTimeErrComponents)
      timeErrFracList.add(timeErrCore1Frac);
  }
  if (bkgUseGaus2)
    timeErrPdfList.add(timeErrGaus2);

  RooAddPdf timeErrPdfAnalytic("timeErrPdf", "timeErrPdf", timeErrPdfList, timeErrFracList, true);
  RooAddPdf timeErrPdfCoreOnly(
      "timeErrPdfCoreOnly", "timeErrPdfCoreOnly",
      RooArgList(timeErrGaus1, timeErrGaus2),
      RooArgList(timeErrCore1Frac),
      true);

  std::unique_ptr<RooDataHist> bkgErrHistData;
  std::unique_ptr<RooHistPdf> bkgErrHistPdf;
  RooAbsPdf *timeErrPdf = &timeErrPdfAnalytic;
  std::unique_ptr<RooFitResult> timeErrResult;
  if (bkgErrPdfOpt == kErrPdfHist)
  {
    bkgErrHistData = std::make_unique<RooDataHist>("bkgErrHistData", "", RooArgSet(obs_timeErr), *timeErrData);
    bkgErrHistPdf = std::make_unique<RooHistPdf>(
        "bkgErrHistPdf", "bkgErrHistPdf", RooArgSet(obs_timeErr), *bkgErrHistData, histPdfInterpolationOrder);
    timeErrPdf = bkgErrHistPdf.get();
  }
  else
  {
    if (useStagedTimeErrFitBkg && bkgUseGaus2)
      timeErrPdfCoreOnly.fitTo(*timeErrData, PrintLevel(-1), SumW2Error(isWeight));
    timeErrResult = std::unique_ptr<RooFitResult>(
        timeErrPdfAnalytic.fitTo(*timeErrData, Save(), PrintLevel(-1), SumW2Error(isWeight)));
    err2_report_params_near_limit(timeErrResult.get(), "bkg_timeErr");
  }

  auto readIntParam = [](TFile &f, const char *name, int fallback = -999) {
    auto *param = dynamic_cast<TParameter<int> *>(f.Get(name));
    return param ? param->GetVal() : fallback;
  };
  auto readDoubleParam = [](TFile &f, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN()) {
    auto *param = dynamic_cast<TParameter<double> *>(f.Get(name));
    return param ? param->GetVal() : fallback;
  };

  std::unique_ptr<TFile> massFile(TFile::Open(massFileName));
  if (!massFile || massFile->IsZombie())
  {
    std::cerr << "ERROR: cannot open mass model file: " << massFileName << std::endl;
    return;
  }

  const int nBkgExpComponents = std::clamp(readIntParam(*massFile, "nBkgExpComponents", 0), 0, 1);
  const int nBkgChebyOrder = std::clamp(readIntParam(*massFile, "nBkgChebyOrder", 0), 0, 6);
  const double nBkgMassVal = readDoubleParam(*massFile, "Nbkg", std::max(1.0, 0.5 * dataSel->numEntries()));
  RooConstVar bkg_mass_lambda("bkg_mass_lambda", "bkg_mass_lambda", readDoubleParam(*massFile, "bkg_mass_lambda", -1.0));
  RooConstVar bkg_mass_p1("bkg_mass_p1", "bkg_mass_p1", readDoubleParam(*massFile, "bkg_mass_p1", 0.0));
  RooConstVar bkg_mass_p2("bkg_mass_p2", "bkg_mass_p2", readDoubleParam(*massFile, "bkg_mass_p2", 0.0));
  RooConstVar bkg_mass_p3("bkg_mass_p3", "bkg_mass_p3", readDoubleParam(*massFile, "bkg_mass_p3", 0.0));
  RooConstVar bkg_mass_p4("bkg_mass_p4", "bkg_mass_p4", readDoubleParam(*massFile, "bkg_mass_p4", 0.0));
  RooConstVar bkg_mass_p5("bkg_mass_p5", "bkg_mass_p5", readDoubleParam(*massFile, "bkg_mass_p5", 0.0));
  RooConstVar bkg_mass_p6("bkg_mass_p6", "bkg_mass_p6", readDoubleParam(*massFile, "bkg_mass_p6", 0.0));
  RooArgList chebyCoeffList;
  if (nBkgChebyOrder >= 1) chebyCoeffList.add(bkg_mass_p1);
  if (nBkgChebyOrder >= 2) chebyCoeffList.add(bkg_mass_p2);
  if (nBkgChebyOrder >= 3) chebyCoeffList.add(bkg_mass_p3);
  if (nBkgChebyOrder >= 4) chebyCoeffList.add(bkg_mass_p4);
  if (nBkgChebyOrder >= 5) chebyCoeffList.add(bkg_mass_p5);
  if (nBkgChebyOrder >= 6) chebyCoeffList.add(bkg_mass_p6);
  std::unique_ptr<RooAbsPdf> bkgMassPdf;
  if (nBkgExpComponents == 1)
    bkgMassPdf = std::make_unique<RooExponential>("bkg_mass", "bkg_mass", obs_mass, bkg_mass_lambda);
  else
    bkgMassPdf = std::make_unique<RooChebychev>("bkg_mass", "bkg_mass", obs_mass, chebyCoeffList);

  const double massSb1Low = 2.6, massSb1High = 2.8;
  const double massSb2Low = 3.3, massSb2High = 3.5;
  const double massSigLow = 2.8, massSigHigh = 3.3;
  obs_mass.setRange("SR", massSigLow, massSigHigh);
  obs_mass.setRange("SB1", massSb1Low, massSb1High);
  obs_mass.setRange("SB2", massSb2Low, massSb2High);
  auto iBkgSR = std::unique_ptr<RooAbsReal>(bkgMassPdf->createIntegral(obs_mass, Range("SR")));
  auto iBkgSB1 = std::unique_ptr<RooAbsReal>(bkgMassPdf->createIntegral(obs_mass, Range("SB1")));
  auto iBkgSB2 = std::unique_ptr<RooAbsReal>(bkgMassPdf->createIntegral(obs_mass, Range("SB2")));
  const double fracBkgSR = iBkgSR ? iBkgSR->getVal() : 0.0;
  const double fracBkgSB = (iBkgSB1 ? iBkgSB1->getVal() : 0.0) + (iBkgSB2 ? iBkgSB2->getVal() : 0.0);
  const double NBkg_SR = nBkgMassVal * fracBkgSR;
  const double NBkg_SB = nBkgMassVal * fracBkgSB;
  const double scaleBkg = NBkg_SB > 0.0 ? NBkg_SR / NBkg_SB : 0.0;
  printf("[MassBkgScale] nBkgMass=%.1f frac(SR)=%.6f frac(SB)=%.6f => NBkg_SR=%.1f NBkg_SB=%.1f scale=%.6f\n",
         nBkgMassVal, fracBkgSR, fracBkgSB, NBkg_SR, NBkg_SB, scaleBkg);

  auto hErrSigSR = std::unique_ptr<TH1>(dataSR ? dataSR->createHistogram("hErrSigSR", obs_timeErr, Binning(timeErrPlotBins, errRange.first, errRange.second)) : nullptr);
  auto hErrBkgSB = std::unique_ptr<TH1>(dataSB->createHistogram("hErrBkgSB", obs_timeErr, Binning(timeErrPlotBins, errRange.first, errRange.second)));
  std::unique_ptr<RooDataHist> sigErrData;
  if (hErrSigSR && hErrBkgSB)
  {
    hErrSigSR->Add(hErrBkgSB.get(), -scaleBkg);
    for (int i = 1; i <= hErrSigSR->GetNbinsX(); ++i)
    {
      if (hErrSigSR->GetBinContent(i) < 0.0)
      {
        hErrSigSR->SetBinContent(i, 0.0);
        hErrSigSR->SetBinError(i, 0.0);
      }
    }
    sigErrData = std::make_unique<RooDataHist>("sigErrData", "", RooArgSet(obs_timeErr), hErrSigSR.get());
  }
  if (!sigErrData || sigErrData->sumEntries() <= 0.0)
  {
    std::cerr << "ERROR: empty SR-SB*scale err distribution" << std::endl;
    return;
  }

  RooRealVar sigErrGaus1Mean("sigErrGaus1Mean", "sigErrGaus1Mean", gaus1MeanInit, errMeanMin, errMeanMax);
  RooRealVar sigErrGaus1Sigma("sigErrGaus1Sigma", "sigErrGaus1Sigma", gaus1SigmaInit, 5e-4, gaus1SigmaMax);
  RooRealVar sigErrGaus2Mean("sigErrGaus2Mean", "sigErrGaus2Mean", gaus2MeanInit, errMeanMin, errMeanMax);
  RooRealVar sigErrGaus2Sigma("sigErrGaus2Sigma", "sigErrGaus2Sigma", gaus2SigmaInit, 1e-3, gaus2SigmaMax);
  RooRealVar sigErrTailMpv("sigErrTailMpv", "sigErrTailMpv", tailMpvInit, errMeanMin, errMeanMax);
  RooRealVar sigErrTailWidth("sigErrTailWidth", "sigErrTailWidth", tailWidthInit, 1e-3, tailWidthMax);
  RooRealVar sigErrTail2Mpv("sigErrTail2Mpv", "sigErrTail2Mpv", tail2MpvInit, errMeanMin, errMeanMax);
  RooRealVar sigErrTail2Width("sigErrTail2Width", "sigErrTail2Width", tail2WidthInit, 1e-3, tail2WidthMax);
  RooRealVar sigErrLognM0("sigErrLognM0", "sigErrLognM0", lognM0Init, errMeanMin, errMeanMax);
  RooRealVar sigErrLognK("sigErrLognK", "sigErrLognK", lognKInit, 0.01, 3.0);
  RooRealVar sigErrCore1FracRatio("sigErrCore1FracRatio", "sigErrCore1FracRatio", fracToRatio(0.60), 1e-3, 1e3);
  RooFormulaVar sigErrCore1Frac("sigErrCore1Frac", "@0/(1.0+@0)", RooArgList(sigErrCore1FracRatio));
  RooRealVar sigErrTailFracRatio("sigErrTailFracRatio", "sigErrTailFracRatio", fracToRatio(0.08), 1e-4, 1e4);
  RooFormulaVar sigErrTailFrac("sigErrTailFrac", "@0/(1.0+@0)", RooArgList(sigErrTailFracRatio));
  RooRealVar sigErrTail2FracRatio("sigErrTail2FracRatio", "sigErrTail2FracRatio", fracToRatio(0.05), 1e-4, 1e4);
  RooFormulaVar sigErrTail2Frac("sigErrTail2Frac", "@0/(1.0+@0)", RooArgList(sigErrTail2FracRatio));
  RooRealVar sigErrLognFracRatio("sigErrLognFracRatio", "sigErrLognFracRatio", fracToRatio(0.03), 1e-4, 1e4);
  RooFormulaVar sigErrLognFrac("sigErrLognFrac", "@0/(1.0+@0)", RooArgList(sigErrLognFracRatio));

  RooGaussian sigErrGaus1("sigErrGaus1", "sigErrGaus1", obs_timeErr, sigErrGaus1Mean, sigErrGaus1Sigma);
  RooGaussian sigErrGaus2("sigErrGaus2", "sigErrGaus2", obs_timeErr, sigErrGaus2Mean, sigErrGaus2Sigma);
  RooLandau sigErrTail("sigErrTail", "sigErrTail", obs_timeErr, sigErrTailMpv, sigErrTailWidth);
  RooLandau sigErrTail2("sigErrTail2", "sigErrTail2", obs_timeErr, sigErrTail2Mpv, sigErrTail2Width);
  RooLognormal sigErrLogn("sigErrLogn", "sigErrLogn", obs_timeErr, sigErrLognM0, sigErrLognK);
  RooArgList sigTimeErrPdfList;
  RooArgList sigTimeErrFracList;
  if (sigUseLandau1)
  {
    sigTimeErrPdfList.add(sigErrTail);
    if (nSigTimeErrComponents > 1)
      sigTimeErrFracList.add(sigErrTailFrac);
  }
  if (sigUseLandau2)
  {
    sigTimeErrPdfList.add(sigErrTail2);
    if (static_cast<int>(sigTimeErrPdfList.getSize()) < nSigTimeErrComponents)
      sigTimeErrFracList.add(sigErrTail2Frac);
  }
  if (sigUseLognormal)
  {
    sigTimeErrPdfList.add(sigErrLogn);
    if (static_cast<int>(sigTimeErrPdfList.getSize()) < nSigTimeErrComponents)
      sigTimeErrFracList.add(sigErrLognFrac);
  }
  if (sigUseGaus1)
  {
    sigTimeErrPdfList.add(sigErrGaus1);
    if (static_cast<int>(sigTimeErrPdfList.getSize()) < nSigTimeErrComponents)
      sigTimeErrFracList.add(sigErrCore1Frac);
  }
  if (sigUseGaus2)
    sigTimeErrPdfList.add(sigErrGaus2);

  RooAddPdf sigTimeErrPdfAnalytic("sigTimeErrPdf", "sigTimeErrPdf", sigTimeErrPdfList, sigTimeErrFracList, true);
  RooAddPdf sigTimeErrPdfCoreOnly(
      "sigTimeErrPdfCoreOnly", "sigTimeErrPdfCoreOnly",
      RooArgList(sigErrGaus1, sigErrGaus2),
      RooArgList(sigErrCore1Frac),
      true);
  std::unique_ptr<RooHistPdf> sigErrHistPdf;
  RooAbsPdf *sigTimeErrPdf = &sigTimeErrPdfAnalytic;
  std::unique_ptr<RooFitResult> sigTimeErrResult;
  if (sigErrPdfOpt == kErrPdfHist)
  {
    sigErrHistPdf = std::make_unique<RooHistPdf>(
        "sigErrHistPdf", "sigErrHistPdf", RooArgSet(obs_timeErr), *sigErrData, histPdfInterpolationOrder);
    sigTimeErrPdf = sigErrHistPdf.get();
  }
  else
  {
    if (useStagedTimeErrFitSig && sigUseGaus2)
      sigTimeErrPdfCoreOnly.fitTo(*sigErrData, PrintLevel(-1), SumW2Error(isWeight));
    sigTimeErrResult = std::unique_ptr<RooFitResult>(
        sigTimeErrPdfAnalytic.fitTo(*sigErrData, Save(), PrintLevel(-1), SumW2Error(isWeight)));
    err2_report_params_near_limit(sigTimeErrResult.get(), "sig_timeErr");
  }

  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  {
    return fr ? fr->findObject(n) : nullptr;
  };

  TCanvas *cTimeErr = new TCanvas("timeErrModel", "timeErrModel", 800, 800);
  TPad *timeErrPad1 = new TPad("timeErrPad1", "timeErrPad1", 0.0, 0.25, 1.0, 1.0);
  timeErrPad1->SetBottomMargin(0.00001);
  timeErrPad1->SetTopMargin(0.08);
  timeErrPad1->SetLogy();
  timeErrPad1->Draw();
  timeErrPad1->cd();

  RooPlot *timeErrPlot = obs_timeErr.frame(Range(errRange.first, errRange.second), Title(""));
  if (isWeight)
    timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), DataError(RooAbsData::SumW2), Name("data"));
  else
    timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), Name("data"));
  timeErrPdf->plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
  if (bkgErrPdfOpt == kErrPdfAnalytic)
  {
    if (bkgUseLandau1)
      timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
    if (bkgUseLandau2)
      timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
    if (bkgUseLognormal)
      timeErrPdf->plotOn(timeErrPlot, Components(timeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
    if (bkgUseGaus1)
      timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
    if (bkgUseGaus2)
      timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
  }
  apply_logy_auto_range(timeErrPlot, "data");
  timeErrPlot->GetYaxis()->SetTitle("Events");
  timeErrPlot->GetYaxis()->SetTitleOffset(1.6);
  timeErrPlot->GetXaxis()->SetTitle("");
  timeErrPlot->Draw("e");

  TLegend timeErrLeg(0.50, 0.66, 0.74, 0.89);
  timeErrLeg.SetBorderSize(0);
  timeErrLeg.SetFillStyle(0);
  timeErrLeg.SetTextSize(0.03);
  if (auto *o = findObj(timeErrPlot, "data"))
    timeErrLeg.AddEntry(o, "Data", "lep");
  if (auto *o = findObj(timeErrPlot, "model"))
    timeErrLeg.AddEntry(o, bkgErrPdfOpt == kErrPdfHist ? "RooHistPdf" : "Fit", "l");
  if (bkgErrPdfOpt == kErrPdfAnalytic && bkgUseLandau1)
    if (auto *o = findObj(timeErrPlot, "tail"))
    timeErrLeg.AddEntry(o, "Landau tail", "l");
  if (bkgErrPdfOpt == kErrPdfAnalytic && bkgUseLandau2)
    if (auto *o = findObj(timeErrPlot, "tail2"))
    timeErrLeg.AddEntry(o, "Landau tail 2", "l");
  if (bkgErrPdfOpt == kErrPdfAnalytic && bkgUseLognormal)
    if (auto *o = findObj(timeErrPlot, "logn"))
    timeErrLeg.AddEntry(o, "Log-normal tail", "l");
  if (bkgErrPdfOpt == kErrPdfAnalytic && bkgUseGaus1)
    if (auto *o = findObj(timeErrPlot, "gaus1"))
    timeErrLeg.AddEntry(o, "Gauss 1", "l");
  if (bkgErrPdfOpt == kErrPdfAnalytic && bkgUseGaus2)
    if (auto *o = findObj(timeErrPlot, "gaus2"))
    timeErrLeg.AddEntry(o, "Gauss 2", "l");
  timeErrLeg.Draw("same");

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, "OO #sqrt{s_{NN}} = 5.36 TeV (9 nb^{-1})");
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
    double xtext = 0.19, y0 = 0.865, dy = -0.05;
    int k = 0;
    tx.DrawLatex(xtext, y0 + dy * k++, "Mass sideband data");
    if (yLow == 0)
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    int minimize = -1;
    int hesse = -1;
    if (timeErrResult)
    {
      for (UInt_t i = 0, n = timeErrResult->numStatusHistory(); i < n; ++i)
      {
        const char *lab = timeErrResult->statusLabelHistory(i);
        if (lab)
        {
          if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
            minimize = timeErrResult->statusCodeHistory(i);
          else if (!strcmp(lab, "HESSE"))
            hesse = timeErrResult->statusCodeHistory(i);
        }
      }
    }
    if (bkgErrPdfOpt == kErrPdfHist)
      tx.DrawLatex(0.19, 0.765, Form("Status : RooHistPdf template (%d)", histPdfInterpolationOrder));
    else
      tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
  }
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.024);
    tp.SetTextFont(42);
    double xtext = 0.72, y0 = 0.87, dy = -0.04;
    int k = 0;
    auto printVar = [&](const char *title, RooAbsReal &var)
    {
      if (auto *rrv = dynamic_cast<RooRealVar *>(&var))
      {
        if (rrv->isConstant())
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
        else
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
        return;
      }
      const double err = timeErrResult ? var.getPropagatedError(*timeErrResult) : 0.0;
      if (err > 0.0 && std::isfinite(err))
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
      else
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
    };
    if (bkgErrPdfOpt == kErrPdfHist)
    {
      tp.DrawLatex(xtext, y0 + dy * k++, Form("bkgErrPdfOpt = %d", bkgErrPdfOpt));
      tp.DrawLatex(xtext, y0 + dy * k++, Form("interp = %d", histPdfInterpolationOrder));
      tp.DrawLatex(xtext, y0 + dy * k++, Form("bins = %d", timeErrPlotBins));
    }
    else
    {
      if (bkgUseGaus1)
      {
        printVar("#mu_{1}", timeErrGaus1Mean);
        printVar("#sigma_{1}", timeErrGaus1Sigma);
      }
      if (bkgUseGaus2)
      {
        printVar("#mu_{2}", timeErrGaus2Mean);
        printVar("#sigma_{2}", timeErrGaus2Sigma);
      }
      if (bkgUseLandau1)
      {
        printVar("mpv_{L}", timeErrTailMpv);
        printVar("#sigma_{L}", timeErrTailWidth);
        if (nBkgTimeErrComponents > 1)
          printVar("f_{tail}", timeErrTailFrac);
      }
      if (bkgUseLandau2)
      {
        printVar("mpv_{L2}", timeErrTail2Mpv);
        printVar("#sigma_{L2}", timeErrTail2Width);
        if (nBkgTimeErrComponents > 1)
          printVar("f_{tail2}", timeErrTail2Frac);
      }
      if (bkgUseLognormal)
      {
        printVar("m_{LN}", timeErrLognM0);
        printVar("k_{LN}", timeErrLognK);
        if (nBkgTimeErrComponents > 1)
          printVar("f_{LN}", timeErrLognFrac);
      }
      if (bkgUseGaus1 && bkgUseGaus2)
        printVar("f_{G1}", timeErrCore1Frac);
    }
  }

  cTimeErr->cd();
  TPad *timeErrPad2 = new TPad("timeErrPad2", "timeErrPad2", 0.0, 0.0, 1.0, 0.25);
  timeErrPad2->SetTopMargin(0.00001);
  timeErrPad2->SetBottomMargin(0.4);
  timeErrPad2->Draw();
  timeErrPad2->cd();

  RooPlot *timeErrPullPlot = obs_timeErr.frame(Range(errRange.first, errRange.second), Title(""));
  RooHist *timeErrPull = timeErrPlot->pullHist("data", "model");
  if (timeErrPull)
    timeErrPullPlot->addPlotable(timeErrPull, "P");
  timeErrPullPlot->GetYaxis()->SetTitle("Pull");
  timeErrPullPlot->GetXaxis()->SetTitle("ctau3D error [mm]");
  timeErrPullPlot->GetXaxis()->CenterTitle();
  timeErrPullPlot->SetMinimum(-8);
  timeErrPullPlot->SetMaximum(8);
  timeErrPullPlot->GetYaxis()->SetNdivisions(505);
  timeErrPullPlot->GetYaxis()->SetTitleSize(0.12);
  timeErrPullPlot->GetYaxis()->SetLabelSize(0.10);
  timeErrPullPlot->GetXaxis()->SetTitleSize(0.15);
  timeErrPullPlot->GetXaxis()->SetLabelSize(0.10);
  timeErrPullPlot->Draw();

  auto timeErrChi = err2_chi2_from_pull(timeErrPull);
  const int timeErrNPar = (bkgErrPdfOpt == kErrPdfHist || !timeErrResult) ? 0 : timeErrResult->floatParsFinal().getSize();
  const int timeErrNdf = std::max(1, timeErrChi.second - timeErrNPar);
  const double timeErrPvalue = TMath::Prob(timeErrChi.first, timeErrNdf);
  {
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.085);
    tc.SetTextFont(42);
    tc.SetTextAlign(33);
    tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", timeErrChi.first, timeErrNdf, timeErrPvalue));
  }

  TLine timeErrLine(errRange.first, 0.0, errRange.second, 0.0);
  timeErrLine.SetLineStyle(2);
  timeErrLine.Draw("same");

  cTimeErr->Print(figName("errBkg"));
  delete cTimeErr;

  TCanvas *cSigTimeErr = new TCanvas("sigTimeErrModel", "sigTimeErrModel", 800, 800);
  TPad *sigTimeErrPad1 = new TPad("sigTimeErrPad1", "sigTimeErrPad1", 0.0, 0.25, 1.0, 1.0);
  sigTimeErrPad1->SetBottomMargin(0.00001);
  sigTimeErrPad1->SetTopMargin(0.08);
  sigTimeErrPad1->SetLogy();
  sigTimeErrPad1->Draw();
  sigTimeErrPad1->cd();

  RooPlot *sigTimeErrPlot = obs_timeErr.frame(Range(errRange.first, errRange.second), Title(""));
  if (isWeight)
    sigErrData->plotOn(sigTimeErrPlot, DataError(RooAbsData::SumW2), Name("data"));
  else
    sigErrData->plotOn(sigTimeErrPlot, Name("data"));
  sigTimeErrPdf->plotOn(sigTimeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
  if (sigErrPdfOpt == kErrPdfAnalytic)
  {
    if (sigUseLandau1)
      sigTimeErrPdf->plotOn(sigTimeErrPlot, Components(sigErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
    if (sigUseLandau2)
      sigTimeErrPdf->plotOn(sigTimeErrPlot, Components(sigErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
    if (sigUseLognormal)
      sigTimeErrPdf->plotOn(sigTimeErrPlot, Components(sigErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
    if (sigUseGaus1)
      sigTimeErrPdf->plotOn(sigTimeErrPlot, Components(sigErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
    if (sigUseGaus2)
      sigTimeErrPdf->plotOn(sigTimeErrPlot, Components(sigErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
  }
  apply_logy_auto_range(sigTimeErrPlot, "data");
  sigTimeErrPlot->GetYaxis()->SetTitle("Events");
  sigTimeErrPlot->GetYaxis()->SetTitleOffset(1.6);
  sigTimeErrPlot->GetXaxis()->SetTitle("");
  sigTimeErrPlot->Draw("e");

  TLegend sigTimeErrLeg(0.50, 0.66, 0.74, 0.89);
  sigTimeErrLeg.SetBorderSize(0);
  sigTimeErrLeg.SetFillStyle(0);
  sigTimeErrLeg.SetTextSize(0.03);
  if (auto *o = findObj(sigTimeErrPlot, "data"))
    sigTimeErrLeg.AddEntry(o, "Data", "lep");
  if (auto *o = findObj(sigTimeErrPlot, "model"))
    sigTimeErrLeg.AddEntry(o, sigErrPdfOpt == kErrPdfHist ? "RooHistPdf" : "Fit", "l");
  if (sigErrPdfOpt == kErrPdfAnalytic && sigUseLandau1)
    if (auto *o = findObj(sigTimeErrPlot, "tail"))
    sigTimeErrLeg.AddEntry(o, "Landau tail", "l");
  if (sigErrPdfOpt == kErrPdfAnalytic && sigUseLandau2)
    if (auto *o = findObj(sigTimeErrPlot, "tail2"))
    sigTimeErrLeg.AddEntry(o, "Landau tail 2", "l");
  if (sigErrPdfOpt == kErrPdfAnalytic && sigUseLognormal)
    if (auto *o = findObj(sigTimeErrPlot, "logn"))
    sigTimeErrLeg.AddEntry(o, "Log-normal tail", "l");
  if (sigErrPdfOpt == kErrPdfAnalytic && sigUseGaus1)
    if (auto *o = findObj(sigTimeErrPlot, "gaus1"))
    sigTimeErrLeg.AddEntry(o, "Gauss 1", "l");
  if (sigErrPdfOpt == kErrPdfAnalytic && sigUseGaus2)
    if (auto *o = findObj(sigTimeErrPlot, "gaus2"))
    sigTimeErrLeg.AddEntry(o, "Gauss 2", "l");
  sigTimeErrLeg.Draw("same");

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, "OO #sqrt{s_{NN}} = 5.36 TeV (9 nb^{-1})");
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
    double xtext = 0.19, y0 = 0.865, dy = -0.05;
    int k = 0;
    tx.DrawLatex(xtext, y0 + dy * k++, "Signal region - scaled sideband");
    if (yLow == 0)
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    int minimize = -1;
    int hesse = -1;
    if (sigTimeErrResult)
    {
      for (UInt_t i = 0, n = sigTimeErrResult->numStatusHistory(); i < n; ++i)
      {
        const char *lab = sigTimeErrResult->statusLabelHistory(i);
        if (lab)
        {
          if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
            minimize = sigTimeErrResult->statusCodeHistory(i);
          else if (!strcmp(lab, "HESSE"))
            hesse = sigTimeErrResult->statusCodeHistory(i);
        }
      }
    }
    if (sigErrPdfOpt == kErrPdfHist)
      tx.DrawLatex(0.19, 0.765, Form("Status : RooHistPdf template (%d)", histPdfInterpolationOrder));
    else
      tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
  }
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.024);
    tp.SetTextFont(42);
    double xtext = 0.72, y0 = 0.87, dy = -0.04;
    int k = 0;
    auto printVar = [&](const char *title, RooAbsReal &var)
    {
      if (auto *rrv = dynamic_cast<RooRealVar *>(&var))
      {
        if (rrv->isConstant())
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
        else
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
        return;
      }
      const double err = sigTimeErrResult ? var.getPropagatedError(*sigTimeErrResult) : 0.0;
      if (err > 0.0 && std::isfinite(err))
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
      else
        tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
    };
    if (sigErrPdfOpt == kErrPdfHist)
    {
      tp.DrawLatex(xtext, y0 + dy * k++, Form("sigErrPdfOpt = %d", sigErrPdfOpt));
      tp.DrawLatex(xtext, y0 + dy * k++, Form("interp = %d", histPdfInterpolationOrder));
      tp.DrawLatex(xtext, y0 + dy * k++, Form("bins = %d", timeErrPlotBins));
    }
    else
    {
      if (sigUseGaus1)
      {
        printVar("#mu_{1}", sigErrGaus1Mean);
        printVar("#sigma_{1}", sigErrGaus1Sigma);
      }
      if (sigUseGaus2)
      {
        printVar("#mu_{2}", sigErrGaus2Mean);
        printVar("#sigma_{2}", sigErrGaus2Sigma);
      }
      if (sigUseLandau1)
      {
        printVar("mpv_{L}", sigErrTailMpv);
        printVar("#sigma_{L}", sigErrTailWidth);
        if (nSigTimeErrComponents > 1)
          printVar("f_{tail}", sigErrTailFrac);
      }
      if (sigUseLandau2)
      {
        printVar("mpv_{L2}", sigErrTail2Mpv);
        printVar("#sigma_{L2}", sigErrTail2Width);
        if (nSigTimeErrComponents > 1)
          printVar("f_{tail2}", sigErrTail2Frac);
      }
      if (sigUseLognormal)
      {
        printVar("m_{LN}", sigErrLognM0);
        printVar("k_{LN}", sigErrLognK);
        if (nSigTimeErrComponents > 1)
          printVar("f_{LN}", sigErrLognFrac);
      }
      if (sigUseGaus1 && sigUseGaus2)
        printVar("f_{G1}", sigErrCore1Frac);
    }
  }

  cSigTimeErr->cd();
  TPad *sigTimeErrPad2 = new TPad("sigTimeErrPad2", "sigTimeErrPad2", 0.0, 0.0, 1.0, 0.25);
  sigTimeErrPad2->SetTopMargin(0.00001);
  sigTimeErrPad2->SetBottomMargin(0.4);
  sigTimeErrPad2->Draw();
  sigTimeErrPad2->cd();

  RooPlot *sigTimeErrPullPlot = obs_timeErr.frame(Range(errRange.first, errRange.second), Title(""));
  RooHist *sigTimeErrPull = sigTimeErrPlot->pullHist("data", "model");
  if (sigTimeErrPull)
    sigTimeErrPullPlot->addPlotable(sigTimeErrPull, "P");
  sigTimeErrPullPlot->GetYaxis()->SetTitle("Pull");
  sigTimeErrPullPlot->GetXaxis()->SetTitle("ctau3D error [mm]");
  sigTimeErrPullPlot->GetXaxis()->CenterTitle();
  sigTimeErrPullPlot->SetMinimum(-8);
  sigTimeErrPullPlot->SetMaximum(8);
  sigTimeErrPullPlot->GetYaxis()->SetNdivisions(505);
  sigTimeErrPullPlot->GetYaxis()->SetTitleSize(0.12);
  sigTimeErrPullPlot->GetYaxis()->SetLabelSize(0.10);
  sigTimeErrPullPlot->GetXaxis()->SetTitleSize(0.15);
  sigTimeErrPullPlot->GetXaxis()->SetLabelSize(0.10);
  sigTimeErrPullPlot->Draw();

  auto sigTimeErrChi = err2_chi2_from_pull(sigTimeErrPull);
  const int sigTimeErrNPar = (sigErrPdfOpt == kErrPdfHist || !sigTimeErrResult) ? 0 : sigTimeErrResult->floatParsFinal().getSize();
  const int sigTimeErrNdf = std::max(1, sigTimeErrChi.second - sigTimeErrNPar);
  const double sigTimeErrPvalue = TMath::Prob(sigTimeErrChi.first, sigTimeErrNdf);
  {
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.085);
    tc.SetTextFont(42);
    tc.SetTextAlign(33);
    tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", sigTimeErrChi.first, sigTimeErrNdf, sigTimeErrPvalue));
  }

  TLine sigTimeErrLine(errRange.first, 0.0, errRange.second, 0.0);
  sigTimeErrLine.SetLineStyle(2);
  sigTimeErrLine.Draw("same");

  cSigTimeErr->Print(figName("errSig"));
  delete cSigTimeErr;

  TFile outFile(resultFileName, "RECREATE");
  TParameter<double>("errLow", errRange.first).Write();
  TParameter<double>("errHigh", errRange.second).Write();
  TParameter<int>("timeErrPlotBins", timeErrPlotBins).Write();
  TParameter<int>("bkgErrPdfOpt", bkgErrPdfOpt).Write();
  TParameter<int>("sigErrPdfOpt", sigErrPdfOpt).Write();
  TParameter<int>("nBkgTimeErrGaussComponents", nBkgTimeErrGaussComponents).Write();
  TParameter<int>("nBkgTimeErrLandauComponents", nBkgTimeErrLandauComponents).Write();
  TParameter<int>("nBkgTimeErrLognormalComponents", nBkgTimeErrLognormalComponents).Write();
  TParameter<int>("nBkgTimeErrComponents", nBkgTimeErrComponents).Write();
  TParameter<int>("nSigTimeErrGaussComponents", nSigTimeErrGaussComponents).Write();
  TParameter<int>("nSigTimeErrLandauComponents", nSigTimeErrLandauComponents).Write();
  TParameter<int>("nSigTimeErrLognormalComponents", nSigTimeErrLognormalComponents).Write();
  TParameter<int>("nSigTimeErrComponents", nSigTimeErrComponents).Write();
  TParameter<double>("scaleBkg", scaleBkg).Write();
  if (timeErrResult)
    timeErrResult->Write("timeErrResult");
  if (sigTimeErrResult)
    sigTimeErrResult->Write("sigTimeErrResult");
  if (bkgErrPdfOpt == kErrPdfAnalytic)
  {
    timeErrGaus1Mean.Write();
    timeErrGaus1Sigma.Write();
    timeErrGaus2Mean.Write();
    timeErrGaus2Sigma.Write();
    timeErrTailMpv.Write();
    timeErrTailWidth.Write();
    timeErrTail2Mpv.Write();
    timeErrTail2Width.Write();
    timeErrLognM0.Write();
    timeErrLognK.Write();
    timeErrCore1FracRatio.Write();
    timeErrTailFracRatio.Write();
    timeErrTail2FracRatio.Write();
    timeErrLognFracRatio.Write();
    timeErrCore1Frac.Write();
    timeErrTailFrac.Write();
    timeErrTail2Frac.Write();
    timeErrLognFrac.Write();
  }
  if (sigErrPdfOpt == kErrPdfAnalytic)
  {
    sigErrGaus1Mean.Write();
    sigErrGaus1Sigma.Write();
    sigErrGaus2Mean.Write();
    sigErrGaus2Sigma.Write();
    sigErrTailMpv.Write();
    sigErrTailWidth.Write();
    sigErrTail2Mpv.Write();
    sigErrTail2Width.Write();
    sigErrLognM0.Write();
    sigErrLognK.Write();
    sigErrCore1FracRatio.Write();
    sigErrTailFracRatio.Write();
    sigErrTail2FracRatio.Write();
    sigErrLognFracRatio.Write();
    sigErrCore1Frac.Write();
    sigErrTailFrac.Write();
    sigErrTail2Frac.Write();
    sigErrLognFrac.Write();
  }
  outFile.Close();

  const TString figErrBkg = figName("errBkg");
  const TString figErrSig = figName("errSig");
  std::cout << "------------------ FIT RESULT FOR ERR BKG --------------" << std::endl;
  if (timeErrResult)
    timeErrResult->Print("v");
  std::cout << "------------------ FIT RESULT FOR ERR SIG --------------" << std::endl;
  if (sigTimeErrResult)
    sigTimeErrResult->Print("v");
  std::cout << "[FIG] err bkg fit : " << figErrBkg << std::endl;
  std::cout << "[FIG] err sig fit : " << figErrSig << std::endl;
}
