/* RooFit tutorial
https://github.com/cofitzpa/roofit_tutorial_solutions/blob/master/roofit_tutorial_solution.C
 *
 * Highlights some of the basic features of RooFit by making
 * fits to the data provided in dataset.root, increasing in complexity
 *
 * This tutorial written by Conor Fitzpatrick
 * conor.fitzpatrick@cern.ch
 * with thanks to Wouter Verkerke and David Kirkby
 * If you spot any bugs or have any problems working through this, let me know...
 *
 */

#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TParameter.h"
#include "TMath.h"
#include "TString.h"
#include "RooHist.h"
#include "saved_fit_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace RooFit;

struct ScopedMacroTimer
{
	TStopwatch sw;
	TString label;

	ScopedMacroTimer(const char *macroName, float ptLow, float ptHigh, float yLow, float yHigh)
		: label(TString::Format("%s(pt=%.2f-%.2f, |y|=%.2f-%.2f)", macroName, ptLow, ptHigh, yLow, yHigh))
	{
		sw.Start();
		std::cout << "[TIMER] start " << label << std::endl;
	}

	~ScopedMacroTimer()
	{
		sw.Stop();
		std::cout << Form("[TIMER] done %s | real %.2f s | cpu %.2f s", label.Data(), sw.RealTime(), sw.CpuTime()) << std::endl;
	}
};

static std::pair<double, int> chi2_from_pull(RooHist *hpull)
{
	double chi2 = 0.;
	int n = 0;
	if (!hpull)
		return {0., 0};
	double x = 0., y = 0.;
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

void ctau_pr(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4, bool drawFromSavedFit = false, bool publish = false)
{
	ScopedMacroTimer timer("ctau_pr", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;

	/************* SOME BASICS FIRST ***************
	 * First things first, make sure you have access to ROOT.
	 * To do this on lxplus, type "SetupProject Gaudi ROOT"
	 *
	 * Next, make sure you've got a folder containing this
	 * macro and the dataset.root ntuple.
	 *
	 * Open dataset.root in a TBrowser and have a look at
	 * the ntuple. In it you'll see two columns labelled mass, time.
	 * This is the same as any other ntuple, there's nothing RooFit specific about it.
	 *
	 * The first few bits of code in this tutorial will load the dataset from
	 * ntuple and plot them in the manner of RooFit. This is to get you started.
	 * Once you know what is going on in the first few lines,
	 * run this macro in the same directory as the one containing dataset.root:
	 * "root -x roofit_tutorial.C"
	 * If you see the plots produced, you're ready to read on...
	 */

	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	enum ErrPdfChoice
	{
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
	};
	int nResolutionComponents = 3; // 1~4 components
	int errPdfOpt = kErrPdfHist;
	const int histPdfInterpolationOrder = 1;
	double errPrefitDataFraction = 0.5;
	double promptRetryPrefitDataFraction = 0.7;
	double promptFinalPrefitDataFraction = 0.9;
	double promptTime1ScaleSeed = 1.0;
	double promptTime2DeltaSeed = 0.8;
	double promptTime3DeltaSeed = 1.7;
	double promptTime4DeltaSeed = 2.5;
	double promptYieldFrac1Seed = 0.40;
	double promptYieldFrac2Seed = 0.30;
	double promptYieldFrac3Seed = 0.10;
	double promptYieldFrac4Seed = 0.05;
	double promptFinalTime1ScaleSeed = promptTime1ScaleSeed;
	double promptFinalTime2DeltaSeed = promptTime2DeltaSeed;
	double promptFinalTime3DeltaSeed = promptTime3DeltaSeed;
	double promptFinalTime4DeltaSeed = promptTime4DeltaSeed;
	double promptFinalYieldFrac1Seed = promptYieldFrac1Seed;
	double promptFinalYieldFrac2Seed = promptYieldFrac2Seed;
	double promptFinalYieldFrac3Seed = promptYieldFrac3Seed;
	double promptFinalYieldFrac4Seed = promptYieldFrac4Seed;
	int nTimeErrGaussComponents = 2;
	int nTimeErrLandauComponents = 1;
	int nTimeErrLognormalComponents = 1;
	if (yLow == 1.6f)
	{
		if (ptLow == 140.0f && ptHigh == 200.0f)
		{
			nResolutionComponents = 2;
		}
		else if (ptLow == 12.0f && ptHigh == 14.0f)
		{
			promptTime1ScaleSeed = 0.85;
			promptTime2DeltaSeed = 1.00;
			promptTime3DeltaSeed = 4.50;
			promptYieldFrac1Seed = 0.82;
			promptYieldFrac2Seed = 0.17;
			promptYieldFrac3Seed = 0.008;
			promptYieldFrac4Seed = 0.0;
			promptFinalTime1ScaleSeed = 0.85;
			promptFinalTime2DeltaSeed = 0.95;
			promptFinalTime3DeltaSeed = 6.00;
			promptFinalYieldFrac1Seed = 0.75;
			promptFinalYieldFrac2Seed = 0.20;
			promptFinalYieldFrac3Seed = 0.05;
			promptFinalYieldFrac4Seed = 0.0;
			promptRetryPrefitDataFraction = 0.85;
			promptFinalPrefitDataFraction = 1.0;
		}
	} else if(yLow==0.0f) {
		if (ptLow == 12.0f && ptHigh == 14.0f)
		{
			nResolutionComponents = 2;
			// nTimeErrGaussComponents = 1;
			// nTimeErrLandauComponents = 1;
			// nTimeErrLognormalComponents = 1;
		}
		else if (ptLow == 16.0f && ptHigh == 20.0f)
		{
			promptTime1ScaleSeed = 0.87;
			promptTime2DeltaSeed = 0.73;
			promptTime3DeltaSeed = 8.00;
			promptYieldFrac1Seed = 0.87;
			promptYieldFrac2Seed = 0.13;
			promptYieldFrac3Seed = 0.0025;
			promptYieldFrac4Seed = 0.0;
			promptFinalTime1ScaleSeed = 0.82;
			promptFinalTime2DeltaSeed = 0.70;
			promptFinalTime3DeltaSeed = 10.0;
			promptFinalYieldFrac1Seed = 0.75;
			promptFinalYieldFrac2Seed = 0.20;
			promptFinalYieldFrac3Seed = 0.05;
			promptFinalYieldFrac4Seed = 0.0;
			promptRetryPrefitDataFraction = 0.85;
			promptFinalPrefitDataFraction = 1.0;
		}
	}
	nResolutionComponents = std::clamp(nResolutionComponents, 1, 4);
	nTimeErrGaussComponents = std::clamp(nTimeErrGaussComponents, 0, 2);
	nTimeErrLandauComponents = std::clamp(nTimeErrLandauComponents, 0, 2);
	nTimeErrLognormalComponents = std::clamp(nTimeErrLognormalComponents, 0, 1);
	if (nTimeErrGaussComponents + nTimeErrLandauComponents + nTimeErrLognormalComponents <= 0)
	{
		std::cerr << "ERROR: at least one prompt timeErr component is required." << std::endl;
		return;
	}
	const bool useGaus1 = (nTimeErrGaussComponents >= 1);
	const bool useGaus2 = (nTimeErrGaussComponents >= 2);
	const bool useLandau1 = (nTimeErrLandauComponents >= 1);
	const bool useLandau2 = (nTimeErrLandauComponents >= 2);
	const bool useLognormal = (nTimeErrLognormalComponents >= 1);
	const int nTimeErrComponents =
			static_cast<int>(useGaus1) + static_cast<int>(useGaus2) +
			static_cast<int>(useLandau1) + static_cast<int>(useLandau2) +
			static_cast<int>(useLognormal);

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	// First we open the actual RooDataSet used in the prompt ctau analysis.
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_PR_globalOn_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	// Trim plotting/fit ranges with robust quantiles to suppress extreme outliers.
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

	// ------------------------------------------------------------------
	// determine fit ranges
	// ------------------------------------------------------------------
	const auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.0001, 0.9999, false);
	auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.999, true);
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

	// ------------------------------------------------------------------
	// output and observable setup
	// ------------------------------------------------------------------
	auto formatTag = [](double value)
	{
		return TString::Format("%.2f", value);
	};
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figDir = TString::Format("figs/%s/ctau_pr", yTag.Data());
	const TString resultDir = TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", resultDir.Data(), figTag.Data());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit && !load_saved_fit_file(savedFitFile, resolutionFileName, "prompt ctau"))
		return;
	if (drawFromSavedFit)
		nResolutionComponents = std::clamp(read_saved_int_param(savedFitFile.get(), "nResolutionComponents", nResolutionComponents), 1, 4);

	auto *massVar = static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	auto *timeVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	auto *timeErrVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	if (!massVar || !timeVar || !timeErrVar)
	{
		std::cerr << "ERROR: required observables mass/ctau3D/ctau3DErr are missing in selected dataset." << std::endl;
		return;
	}
	timeErrVar->setRange(errRange.first, errRange.second);
	timeErrVar->setMin(errRange.first);
	timeErrVar->setMax(errRange.second);
	int timeErrPlotBins = std::max(2, timeErrVar->getBins() * 2);
	int timeErrTrimBins = timeErrPlotBins;
	auto timeErrScanData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSel->reduce(RooArgSet(*timeErrVar))));
	auto timeErrScanHist = std::unique_ptr<TH1>(timeErrScanData ? timeErrScanData->createHistogram(
			"timeErrScanHist", *timeErrVar, Binning(timeErrTrimBins)) : nullptr);
	if (timeErrScanHist)
	{
		const double errBinWidth = timeErrScanHist->GetBinWidth(1);
		const int peakBin = timeErrScanHist->GetMaximumBin();
		for (int i = peakBin + 1; i <= timeErrScanHist->GetNbinsX(); ++i)
		{
			if (timeErrScanHist->GetBinContent(i) <= 0.0)
			{
				errRange.second = timeErrScanHist->GetXaxis()->GetBinLowEdge(i);
				break;
			}
		}
		const int rebinned = static_cast<int>(std::lround((errRange.second - errRange.first) / errBinWidth));
		timeErrTrimBins = std::max(1, rebinned);
	}
	if (errRange.second <= errRange.first)
	{
		std::cerr << "ERROR: invalid ctau3DErr range after peak-side trimming: ["
		          << errRange.first << ", " << errRange.second << "]" << std::endl;
		return;
	}
	TString cutErrTrim = Form("ctau3DErr >= %g && ctau3DErr < %g", errRange.first, errRange.second);
	dataSel.reset(static_cast<RooDataSet *>(dataSel->reduce(cutErrTrim)));
	if (!dataSel || dataSel->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after ctau3DErr trimming: " << cutErrTrim << std::endl;
		return;
	}
	timeErrPlotBins = timeErrTrimBins;

	massVar = static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	timeVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	timeErrVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	if (!massVar || !timeVar || !timeErrVar)
	{
		std::cerr << "ERROR: required observables mass/ctau3D/ctau3DErr are missing after trimming." << std::endl;
		return;
	}
	RooRealVar obs_mass("mass", "mass", 2.6, 3.5);
	RooRealVar obs_time("ctau3D", "ctau3D", ctRange.first, ctRange.second);
	RooRealVar obs_timeErr("ctau3DErr", "ctau3DErr", errRange.first, errRange.second);
	obs_mass.setBins(massVar->getBins());
	obs_time.setBins(timeVar->getBins());
	obs_timeErr.setBins(timeErrVar->getBins());
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_mass.setRange(2.6, 3.5);
	obs_mass.setMin(2.6);
	obs_mass.setMax(3.5);
	obs_time.SetTitle("#font[12]{l}_{J/#psi}");
	obs_time.setUnit("mm");
	obs_timeErr.SetTitle("event-by-event ctau error");
	obs_timeErr.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	obs_timeErr.setRange(errRange.first, errRange.second);
	obs_timeErr.setMin(errRange.first);
	obs_timeErr.setMax(errRange.second);
	const int massPlotBins = std::max(2, obs_mass.getBins() * 2);
	const int timePlotBins = std::max(2, obs_time.getBins() * 2);

	RooDataSet *data = dataSel.get();

	// ------------------------------------------------------------------
	// build ctau-error model
	// ------------------------------------------------------------------
	// Use an analytic core+tail model for the per-event error distribution.
	auto timeErrData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_timeErr))));
	RooArgSet condObs(obs_timeErr);
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
	auto fracToRatio = [](double frac)
	{
		const double eps = 1e-6;
		const double f = std::clamp(frac, eps, 1.0 - eps);
		return f / (1.0 - f);
	};
	RooRealVar timeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", gaus1MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", gaus1SigmaInit, 5e-4, std::max(0.15 * errSpan, 2e-3));
	RooRealVar timeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", gaus2MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", gaus2SigmaInit, 1e-3, std::max(0.30 * errSpan, 4e-3));
	RooRealVar timeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", tailMpvInit, errRange.first, errRange.second);
	RooRealVar timeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", tailWidthInit, 5e-4, std::max(0.60 * errSpan, 6e-3));
	RooRealVar timeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", tail2MpvInit, errRange.first, errRange.second);
	RooRealVar timeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", tail2WidthInit, 5e-4, std::max(1.00 * errSpan, 1.0e-2));
	RooRealVar timeErrLognM0("timeErrLognM0", "timeErrLognM0", lognM0Init, std::max(errRange.first, 1e-6), errRange.second);
	RooRealVar timeErrLognK("timeErrLognK", "timeErrLognK", lognKInit, 0.05, 0.95);
	RooRealVar timeErrCore1FracRatio("timeErrCore1FracRatio", "timeErrCore1FracRatio", fracToRatio(0.60), 1e-4, 1e4);
	RooFormulaVar timeErrCore1Frac("timeErrCore1Frac", "@0/(1.0+@0)", RooArgList(timeErrCore1FracRatio));
	RooRealVar timeErrTailFracRatio("timeErrTailFracRatio", "timeErrTailFracRatio", fracToRatio(0.08), 1e-5, 1e5);
	RooFormulaVar timeErrTailFrac("timeErrTailFrac", "@0/(1.0+@0)", RooArgList(timeErrTailFracRatio));
	RooRealVar timeErrTail2FracRatio("timeErrTail2FracRatio", "timeErrTail2FracRatio", fracToRatio(0.05), 1e-5, 1e5);
	RooFormulaVar timeErrTail2Frac("timeErrTail2Frac", "@0/(1.0+@0)", RooArgList(timeErrTail2FracRatio));
	RooRealVar timeErrLognFracRatio("timeErrLognFracRatio", "timeErrLognFracRatio", fracToRatio(0.03), 1e-5, 1e5);
	RooFormulaVar timeErrLognFrac("timeErrLognFrac", "@0/(1.0+@0)", RooArgList(timeErrLognFracRatio));
	RooGaussian timeErrGaus1("timeErrGaus1", "timeErrGaus1", obs_timeErr, timeErrGaus1Mean, timeErrGaus1Sigma);
	RooGaussian timeErrGaus2("timeErrGaus2", "timeErrGaus2", obs_timeErr, timeErrGaus2Mean, timeErrGaus2Sigma);
	RooLandau timeErrTail("timeErrTail", "timeErrTail", obs_timeErr, timeErrTailMpv, timeErrTailWidth);
	RooLandau timeErrTail2("timeErrTail2", "timeErrTail2", obs_timeErr, timeErrTail2Mpv, timeErrTail2Width);
	RooLognormal timeErrLogn("timeErrLogn", "timeErrLogn", obs_timeErr, timeErrLognM0, timeErrLognK);
	RooArgList timeErrPdfList;
	RooArgList timeErrFracList;
	if (useLandau1)
	{
		timeErrPdfList.add(timeErrTail);
		if (nTimeErrComponents > 1)
			timeErrFracList.add(timeErrTailFrac);
	}
	if (useLandau2)
	{
		timeErrPdfList.add(timeErrTail2);
		if (static_cast<int>(timeErrPdfList.getSize()) < nTimeErrComponents)
			timeErrFracList.add(timeErrTail2Frac);
	}
	if (useLognormal)
	{
		timeErrPdfList.add(timeErrLogn);
		if (static_cast<int>(timeErrPdfList.getSize()) < nTimeErrComponents)
			timeErrFracList.add(timeErrLognFrac);
	}
	if (useGaus1)
	{
		timeErrPdfList.add(timeErrGaus1);
		if (static_cast<int>(timeErrPdfList.getSize()) < nTimeErrComponents)
			timeErrFracList.add(timeErrCore1Frac);
	}
	if (useGaus2)
		timeErrPdfList.add(timeErrGaus2);
	std::unique_ptr<RooAddPdf> timeErrPdfAnalytic = std::make_unique<RooAddPdf>(
			"timeErrPdf", "timeErrPdf", timeErrPdfList, timeErrFracList, true);
	std::unique_ptr<RooDataHist> timeErrHistData;
	std::unique_ptr<TH1> timeErrHistTemplate;
	std::unique_ptr<RooHistPdf> timeErrPdfHist;
	RooAbsPdf *timeErrPdf = timeErrPdfAnalytic.get();
	RooFitResult *timeErrResult = nullptr;
	std::unique_ptr<RooFitResult> savedTimeErrResult;
	if (errPdfOpt == kErrPdfHist)
	{
		timeErrHistTemplate = std::unique_ptr<TH1>(timeErrData->createHistogram(
				"ctau_pr_hTimeErr", obs_timeErr, Binning(timeErrPlotBins, errRange.first, errRange.second)));
		if (!timeErrHistTemplate)
		{
			std::cerr << "ERROR: failed to build ctau3DErr histogram template." << std::endl;
			return;
		}
		timeErrHistData = std::make_unique<RooDataHist>(
				"timeErrHistData", "timeErrHistData",
				RooArgSet(obs_timeErr), timeErrHistTemplate.get());
		timeErrPdfHist = std::make_unique<RooHistPdf>(
				"timeErrHistPdf", "timeErrHistPdf",
				RooArgSet(obs_timeErr), *timeErrHistData, histPdfInterpolationOrder);
		timeErrPdf = timeErrPdfHist.get();
	}
	else
	{
		if (drawFromSavedFit)
		{
			savedTimeErrResult = clone_saved_fit_result(savedFitFile.get(), "timeErrResult");
			timeErrResult = savedTimeErrResult.get();
			if (!timeErrResult)
			{
				std::cerr << "ERROR: timeErrResult not found in saved prompt ctau file: " << resolutionFileName << std::endl;
				return;
			}
			apply_saved_fit_result(timeErrResult, *timeErrPdfAnalytic, RooArgSet(obs_timeErr));
		}
		else
			timeErrResult = timeErrPdfAnalytic->fitTo(*timeErrData, Save(true), PrintLevel(-1), SumW2Error(isWeight), PrefitDataFraction(errPrefitDataFraction));
		if (!drawFromSavedFit && timeErrResult && timeErrResult->status() != 0)
		{
			std::cout << "[WARN] timeErr fit did not converge (status=" << timeErrResult->status()
								<< "), retrying once." << std::endl;
			delete timeErrResult;
			timeErrResult = timeErrPdfAnalytic->fitTo(*timeErrData, Save(true), PrintLevel(-1), SumW2Error(isWeight), PrefitDataFraction(errPrefitDataFraction));
		}
		timeErrGaus1Mean.setConstant(true);
		timeErrGaus1Sigma.setConstant(true);
		timeErrGaus2Mean.setConstant(true);
		timeErrGaus2Sigma.setConstant(true);
		timeErrTailMpv.setConstant(true);
		timeErrTailWidth.setConstant(true);
		timeErrTail2Mpv.setConstant(true);
		timeErrTail2Width.setConstant(true);
		timeErrLognM0.setConstant(true);
		timeErrLognK.setConstant(true);
		timeErrCore1FracRatio.setConstant(true);
		timeErrTailFracRatio.setConstant(true);
		timeErrTail2FracRatio.setConstant(true);
		timeErrLognFracRatio.setConstant(true);
	}

	// ------------------------------------------------------------------
	// draw ctau-error model
	// ------------------------------------------------------------------
	// Save a control plot of the per-event error distribution and its fitted model.
	// RooPlot stores drawn objects by name; use that to build legends safely.
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
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errRange.first, errRange.second), DataError(RooAbsData::SumW2), Name("data"));
	else
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errRange.first, errRange.second), Name("data"));
	timeErrPdf->plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (errPdfOpt == kErrPdfAnalytic && useLandau1)
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
	if (errPdfOpt == kErrPdfAnalytic && useLandau2)
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
	if (errPdfOpt == kErrPdfAnalytic && useLognormal)
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
	if (errPdfOpt == kErrPdfAnalytic && useGaus1)
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
	if (errPdfOpt == kErrPdfAnalytic && useGaus2)
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
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
		timeErrLeg.AddEntry(o, errPdfOpt == kErrPdfHist ? "RooHistPdf" : "Fit", "l");
	if (errPdfOpt == kErrPdfAnalytic && useLandau1)
		if (auto *o = findObj(timeErrPlot, "tail"))
		timeErrLeg.AddEntry(o, "Landau tail", "l");
	if (errPdfOpt == kErrPdfAnalytic && useLandau2)
		if (auto *o = findObj(timeErrPlot, "tail2"))
		timeErrLeg.AddEntry(o, "Landau tail 2", "l");
	if (errPdfOpt == kErrPdfAnalytic && useLognormal)
		if (auto *o = findObj(timeErrPlot, "logn"))
		timeErrLeg.AddEntry(o, "Log-normal tail", "l");
	if (errPdfOpt == kErrPdfAnalytic && useGaus1)
		if (auto *o = findObj(timeErrPlot, "gaus1"))
		timeErrLeg.AddEntry(o, "Gauss 1", "l");
	if (errPdfOpt == kErrPdfAnalytic && useGaus2)
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
			tx.DrawLatex(xtext, y0 + dy * k++, "Prompt J/#psi MC");
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
			const int status = timeErrResult ? timeErrResult->status() : -1;
			int hesse = -1;
			if (timeErrResult)
			{
				for (UInt_t i = 0, n = timeErrResult->numStatusHistory(); i < n; ++i)
				{
					const char *lab = timeErrResult->statusLabelHistory(i);
					if (lab && TString(lab) == "HESSE")
					{
						hesse = timeErrResult->statusCodeHistory(i);
						break;
					}
				}
			}
			if (!publish)
			{
				if (errPdfOpt == kErrPdfHist)
					tx.DrawLatex(0.19, 0.765, Form("Status : RooHistPdf template (%d)", histPdfInterpolationOrder));
				else
					tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
			}
		}
		if (!publish)
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
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
		};
			if (errPdfOpt == kErrPdfHist)
			{
				tp.DrawLatex(xtext, y0 + dy * k++, Form("errPdfOpt = %d", errPdfOpt));
				tp.DrawLatex(xtext, y0 + dy * k++, Form("interp = %d", histPdfInterpolationOrder));
				tp.DrawLatex(xtext, y0 + dy * k++, Form("bins = %d", obs_timeErr.getBins()));
			}
			else
			{
				if (useGaus1)
				{
					printVar("#mu_{1}", timeErrGaus1Mean);
					printVar("#sigma_{1}", timeErrGaus1Sigma);
				}
				if (useGaus2)
				{
					printVar("#mu_{2}", timeErrGaus2Mean);
					printVar("#sigma_{2}", timeErrGaus2Sigma);
				}
				if (useLandau1)
				{
					printVar("mpv_{L}", timeErrTailMpv);
					printVar("#sigma_{L}", timeErrTailWidth);
					if (nTimeErrComponents > 1)
						printVar("f_{tail}", timeErrTailFrac);
				}
				if (useLandau2)
				{
					printVar("mpv_{L2}", timeErrTail2Mpv);
					printVar("#sigma_{L2}", timeErrTail2Width);
					if (nTimeErrComponents > 1)
						printVar("f_{tail2}", timeErrTail2Frac);
				}
				if (useLognormal)
				{
					printVar("m_{LN}", timeErrLognM0);
					printVar("k_{LN}", timeErrLognK);
					if (nTimeErrComponents > 1)
						printVar("f_{LN}", timeErrLognFrac);
				}
				if (useGaus1 && nTimeErrComponents > 1)
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

	auto timeErrChi = chi2_from_pull(timeErrPull);
	const int timeErrNPar = (errPdfOpt == kErrPdfAnalytic && timeErrResult) ? timeErrResult->floatParsFinal().getSize() : 0;
	const int timeErrNdf = std::max(1, timeErrChi.second - timeErrNPar);
	const double timeErrPvalue = TMath::Prob(timeErrChi.first, timeErrNdf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", timeErrChi.first, timeErrNdf));
	}

	TLine timeErrLine(errRange.first, 0.0, errRange.second, 0.0);
	timeErrLine.SetLineStyle(2);
	timeErrLine.Draw("same");

	cTimeErr->Print(figName("timeerr_model"));
	delete cTimeErr;

	// ------------------------------------------------------------------
	// quick data preview
	// ------------------------------------------------------------------
	// Next, let's plot the RooDataSet in each dimension to see what
	// the distributions look like:

	// Create a TCanvas
	TCanvas *c = new TCanvas("data", "data", 1024, 512);
	// Chop it in two to show both dimensions:
	c->Divide(2);

	// First we'll plot the mass on the left hand side
	c->cd(1);
	// Adjust the margins as the root default is terrible
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	// Now the important bit: A RooPlot is made in the mass dimension:
	RooPlot *massplot = obs_mass.frame();

	// We plot the data on this RooPlot:
	if (isWeight)
		data->plotOn(massplot, Binning(massPlotBins), DataError(RooAbsData::SumW2));
	else
		data->plotOn(massplot, Binning(massPlotBins));
	massplot->GetYaxis()->SetTitleOffset(1.6);
	massplot->Draw();

	// Now we plot the time distribution on the RHS:
	c->cd(2);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	// Logarithmic axis as this looks better:
	gPad->SetLogy();
	RooPlot *timeplot = obs_time.frame();
	if (isWeight)
		data->plotOn(timeplot, Binning(timePlotBins), Name("data_preview"), DataError(RooAbsData::SumW2));
	else
		data->plotOn(timeplot, Binning(timePlotBins), Name("data_preview"));
	apply_logy_auto_range(timeplot, "data_preview");
	timeplot->GetYaxis()->SetTitleOffset(1.6);
	timeplot->Draw();
	c->Draw();

	/***********************TUTORIAL STARTS HERE************************
	 * Your task is to try and fit to this dataset to extract the signal
	 * yield and lifetime.
	 * We'll start with the signal yield, fitting to the mass distribution.
	 * Then we'll try to extract the signal lifetime from the time distribution
	 * Finally, we'll simultaneously fit to both time and mass to see if it
	 * improves the errors.
	 *
	 * The tutorial is written in such a way that you should uncomment
	 * lines one at a time and fill them in as you go along. This way
	 * the code can be re-run at each step to see what happens.
	 */

	// ------------------------------------------------------------------
	// build ctau resolution model
	// ------------------------------------------------------------------
	// Use only lifetime fit (mass fit skipped on purpose).
	// THE CTAU TIME PDF
	// Modeled as a sum of up to four Gaussian resolution components.
	// Use per-event ctau errors as sigma scale factors while keeping the mean scale fixed to 1.
	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", 1.0);
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", 0.0, -0.5, 0.5);
	ctauTime1Mean.setConstant(true);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale", promptTime1ScaleSeed, 0.5, 2.0);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", 0.0, -0.5, 0.5);
	ctauTime2Mean.setConstant(true);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", promptTime2DeltaSeed, 0.001, 5.0);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", 0.0, -0.5, 0.5);
	ctauTime3Mean.setConstant(true);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", promptTime3DeltaSeed, 0.001, 10.0);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", 0.0, -0.5, 0.5);
	ctauTime4Mean.setConstant(true);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", promptTime4DeltaSeed, 0.05, 15.0);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale, ctauMeanScale, obs_timeErr);

	// ------------------------------------------------------------------
	// build ctau fit model
	// ------------------------------------------------------------------
		// Yield parameters for the ctau Gaussian components.
		RooRealVar Nctau1("Nctau1", "Nctau1", 0.0, data->numEntries());
		RooRealVar Nctau2("Nctau2", "Nctau2", 0.0, data->numEntries());
		RooRealVar Nctau3("Nctau3", "Nctau3", 0.0, data->numEntries());
		RooRealVar Nctau4("Nctau4", "Nctau4", 0.0, data->numEntries());
		const auto seedPromptYield = [&](double frac) {
			return std::clamp(data->numEntries() * frac, 0.0, static_cast<double>(data->numEntries()));
		};
		const auto resetPromptTimeSeeds = [&]() {
			ctauTime1Scale.setVal(promptTime1ScaleSeed);
			ctauTime2Delta.setVal(promptTime2DeltaSeed);
			ctauTime3Delta.setVal(promptTime3DeltaSeed);
			ctauTime4Delta.setVal(promptTime4DeltaSeed);
			Nctau1.setVal(seedPromptYield(promptYieldFrac1Seed));
			Nctau2.setVal(seedPromptYield(promptYieldFrac2Seed));
			Nctau3.setVal(seedPromptYield(promptYieldFrac3Seed));
			Nctau4.setVal(seedPromptYield(promptYieldFrac4Seed));
		};
		const auto resetPromptTimeSeedsFinal = [&]() {
			ctauTime1Scale.setVal(promptFinalTime1ScaleSeed);
			ctauTime2Delta.setVal(promptFinalTime2DeltaSeed);
			ctauTime3Delta.setVal(promptFinalTime3DeltaSeed);
			ctauTime4Delta.setVal(promptFinalTime4DeltaSeed);
			Nctau1.setVal(seedPromptYield(promptFinalYieldFrac1Seed));
			Nctau2.setVal(seedPromptYield(promptFinalYieldFrac2Seed));
			Nctau3.setVal(seedPromptYield(promptFinalYieldFrac3Seed));
			Nctau4.setVal(seedPromptYield(promptFinalYieldFrac4Seed));
		};
		resetPromptTimeSeeds();

		// THE LIFETIME FIT
		// The total ctau model is the prompt resolution mixture conditioned on the fitted error PDF.
		RooArgList promptTimePdfList;
		RooArgList promptTimeYieldList;
		promptTimePdfList.add(ctauTime1);
		promptTimeYieldList.add(Nctau1);
		if (nResolutionComponents >= 2)
		{
			promptTimePdfList.add(ctauTime2);
			promptTimeYieldList.add(Nctau2);
		}
		if (nResolutionComponents >= 3)
		{
			promptTimePdfList.add(ctauTime3);
			promptTimeYieldList.add(Nctau3);
		}
		else
		{
			Nctau3.setVal(0.0);
			Nctau3.setConstant(true);
			ctauTime3Mean.setConstant(true);
			ctauTime3Delta.setConstant(true);
		}
		if (nResolutionComponents >= 4)
		{
			promptTimePdfList.add(ctauTime4);
			promptTimeYieldList.add(Nctau4);
		}
		else
		{
			Nctau4.setVal(0.0);
			Nctau4.setConstant(true);
			ctauTime4Mean.setConstant(true);
			ctauTime4Delta.setConstant(true);
		}
		if (nResolutionComponents < 2)
		{
			Nctau2.setVal(0.0);
			Nctau2.setConstant(true);
			ctauTime2Mean.setConstant(true);
			ctauTime2Delta.setConstant(true);
		}
	auto time_pdf_core = std::make_unique<RooAddPdf>(
			"time_pdf_core", "time_pdf_core",
			promptTimePdfList, promptTimeYieldList);
	RooProdPdf time_pdf("time_pdf", "time_pdf", RooArgSet(*timeErrPdf), Conditional(RooArgSet(*time_pdf_core), RooArgSet(obs_time)));
	RooFitResult *time_result = nullptr;
	std::unique_ptr<RooFitResult> savedTimeResult;
	if (drawFromSavedFit)
	{
		savedTimeResult = clone_saved_fit_result(savedFitFile.get(), "timeResult");
		time_result = savedTimeResult.get();
		if (!time_result)
		{
			std::cerr << "ERROR: timeResult not found in saved prompt ctau file: " << resolutionFileName << std::endl;
			return;
		}
		apply_saved_fit_result(time_result, time_pdf, RooArgSet(obs_time, obs_timeErr));
	}
	else
		time_result = time_pdf.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
									 PrintLevel(-1), PrefitDataFraction(errPrefitDataFraction),
									 RecoverFromUndefinedRegions(1.0));
	if (!drawFromSavedFit && time_result && time_result->status() != 0)
	{
		std::cout << "[WARN] prompt time fit did not converge (status=" << time_result->status()
							<< "), retrying with reset seeds and stronger fit options." << std::endl;
		delete time_result;
		resetPromptTimeSeeds();
		time_result = time_pdf.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
									 PrintLevel(-1),
									 PrefitDataFraction(promptRetryPrefitDataFraction),
									 Strategy(1), Offset(true),
									 RecoverFromUndefinedRegions(1.0));
	}
	if (!drawFromSavedFit && time_result && time_result->status() != 0)
	{
		std::cout << "[WARN] prompt time fit still did not converge (status=" << time_result->status()
							<< "), retrying with Strategy(2)." << std::endl;
		delete time_result;
		resetPromptTimeSeedsFinal();
		time_result = time_pdf.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
									 PrintLevel(-1),
									 PrefitDataFraction(promptFinalPrefitDataFraction),
									 Strategy(2), Offset(true),
									 RecoverFromUndefinedRegions(1.0));
	}
	if (time_result)
		time_result->Print();

		const double nCtauTotal = std::max(1e-12,
			Nctau1.getVal() +
			(nResolutionComponents >= 2 ? Nctau2.getVal() : 0.0) +
			(nResolutionComponents >= 3 ? Nctau3.getVal() : 0.0) +
			(nResolutionComponents >= 4 ? Nctau4.getVal() : 0.0));
		const double ctauFrac1 = Nctau1.getVal() / nCtauTotal;
		const double ctauFrac2 = Nctau2.getVal() / nCtauTotal;
		const double ctauFrac3 = Nctau3.getVal() / nCtauTotal;
	if (!drawFromSavedFit)
	{
		TFile resolutionFile(resolutionFileName, "RECREATE");
		if (!resolutionFile.IsZombie())
		{
			if (time_result)
				time_result->Write("timeResult");
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			TParameter<double>("ctauMeanScale", ctauMeanScale.getVal()).Write();
			TParameter<double>("ctauTime1Mean", ctauTime1Mean.getVal()).Write();
			TParameter<double>("ctauTime1Scale", ctauTime1Scale.getVal()).Write();
			if (nResolutionComponents >= 2)
			{
				TParameter<double>("ctauTime2Mean", ctauTime2Mean.getVal()).Write();
				TParameter<double>("ctauTime2Scale", ctauTime2Scale.getVal()).Write();
				TParameter<double>("ctauFrac1", ctauFrac1).Write();
			}
			if (nResolutionComponents >= 3)
			{
				TParameter<double>("ctauTime3Mean", ctauTime3Mean.getVal()).Write();
				TParameter<double>("ctauTime3Scale", ctauTime3Scale.getVal()).Write();
				TParameter<double>("ctauFrac2", ctauFrac2).Write();
			}
			if (nResolutionComponents >= 4)
			{
				TParameter<double>("ctauTime4Mean", ctauTime4Mean.getVal()).Write();
				TParameter<double>("ctauTime4Scale", ctauTime4Scale.getVal()).Write();
				TParameter<double>("ctauFrac3", ctauFrac3).Write();
			}
			resolutionFile.Write();
		}
		else
		{
			std::cerr << "ERROR: cannot create resolution file: " << resolutionFileName << std::endl;
		}
	}
	if (drawFromSavedFit)
		std::cout << "[PlotOnly] Loaded saved prompt ctau fit and left ROOT file unchanged: " << resolutionFileName << std::endl;
	else
		std::cout << "Saved ctau prompt resolution parameters to " << resolutionFileName << std::endl;

	RooProdPdf time_pdf_component1(
			"time_pdf_component1",
			"time_pdf_component1",
			RooArgSet(*timeErrPdf),
			Conditional(RooArgSet(ctauTime1), RooArgSet(obs_time)));
	RooProdPdf time_pdf_component2(
			"time_pdf_component2",
			"time_pdf_component2",
			RooArgSet(*timeErrPdf),
			Conditional(RooArgSet(ctauTime2), RooArgSet(obs_time)));
	RooProdPdf time_pdf_component3(
			"time_pdf_component3",
			"time_pdf_component3",
			RooArgSet(*timeErrPdf),
			Conditional(RooArgSet(ctauTime3), RooArgSet(obs_time)));
	RooProdPdf time_pdf_component4(
			"time_pdf_component4",
			"time_pdf_component4",
			RooArgSet(*timeErrPdf),
			Conditional(RooArgSet(ctauTime4), RooArgSet(obs_time)));

	// ------------------------------------------------------------------
	// draw ctau fit
	// ------------------------------------------------------------------
	TCanvas *cLifetime = new TCanvas("c_lifetime", "c_lifetime", 800, 800);
	TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
	pad1->SetBottomMargin(0.00001);
	pad1->SetTopMargin(0.08);
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();

	RooPlot *timefitplot = obs_time.frame(Range(ctRange.first, ctRange.second), Title(""));
	if (isWeight)
		data->plotOn(timefitplot, Binning(timePlotBins), DataError(RooAbsData::SumW2), Name("data"));
	else
		data->plotOn(timefitplot, Binning(timePlotBins), Name("data"));
	time_pdf.plotOn(timefitplot, LineColor(kBlack), LineWidth(2), Name("model"));
	time_pdf_component1.plotOn(
			timefitplot,
			LineColor(kBlue + 1),
			LineStyle(kDashed),
			LineWidth(2),
			Normalization(Nctau1.getVal(), RooAbsReal::NumEvent),
			Name("ctau_component_1"));
	if (nResolutionComponents >= 2)
	{
		time_pdf_component2.plotOn(
				timefitplot,
				LineColor(kRed + 1),
				LineStyle(kDashed),
				LineWidth(2),
				Normalization(Nctau2.getVal(), RooAbsReal::NumEvent),
				Name("ctau_component_2"));
	}
		if (nResolutionComponents >= 3)
		{
			time_pdf_component3.plotOn(
					timefitplot,
					LineColor(kMagenta + 1),
					LineStyle(kDashed),
					LineWidth(2),
					Normalization(Nctau3.getVal(), RooAbsReal::NumEvent),
					Name("ctau_component_3"));
		}
		if (nResolutionComponents >= 4)
		{
			time_pdf_component4.plotOn(
					timefitplot,
					LineColor(kOrange + 7),
					LineStyle(kDashed),
					LineWidth(2),
					Normalization(Nctau4.getVal(), RooAbsReal::NumEvent),
					Name("ctau_component_4"));
		}

	apply_logy_auto_range(timefitplot, "data");
	timefitplot->GetYaxis()->SetTitle("Events");
	timefitplot->GetYaxis()->SetTitleOffset(1.6);
	timefitplot->GetXaxis()->SetTitle("");
	timefitplot->Draw("e");

	TLegend leg(0.50, 0.70, 0.72, 0.89);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.03);
	if (auto *o = findObj(timefitplot, "data"))
		leg.AddEntry(o, "Data", "lep");
	if (auto *o = findObj(timefitplot, "model"))
		leg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(timefitplot, "ctau_component_1"))
		leg.AddEntry(o, "Core 1", "l");
	if (nResolutionComponents >= 2)
		if (auto *o = findObj(timefitplot, "ctau_component_2"))
		leg.AddEntry(o, "Core 2", "l");
		if (nResolutionComponents >= 3)
			if (auto *o = findObj(timefitplot, "ctau_component_3"))
				leg.AddEntry(o, "Core 3", "l");
		if (nResolutionComponents >= 4)
			if (auto *o = findObj(timefitplot, "ctau_component_4"))
				leg.AddEntry(o, "Core 4", "l");
	leg.Draw("same");

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
			tx.DrawLatex(xtext, y0 + dy * k++, "Prompt J/#psi MC");
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
		const int status = time_result ? time_result->status() : -1;
		int hesse = -1;
		if (time_result)
		{
			for (UInt_t i = 0, n = time_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = time_result->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = time_result->statusCodeHistory(i);
					break;
				}
			}
		}
		if (!publish)
			if (!publish)
				tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
	}
	if (!publish)
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.73, y0 = 0.87, dy = -0.045;
		int k = 0;
		auto printVar = [&](const char *title, RooAbsReal &var)
		{
			auto *rrv = dynamic_cast<RooRealVar *>(&var);
			if (rrv && rrv->isConstant())
			{
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
				return;
			}
			const double err = time_result ? var.getPropagatedError(*time_result) : 0.0;
			tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
			};
			printVar("N_{1}", Nctau1);
			if (nResolutionComponents >= 2)
				printVar("N_{2}", Nctau2);
			printVar("#mu_{1}", ctauTime1Mean);
			printVar("s_{1}", ctauTime1Scale);
			if (nResolutionComponents >= 2)
			{
				printVar("s_{2}", ctauTime2Scale);
				printVar("#Delta s_{21}", ctauTime2Delta);
			}
			if (nResolutionComponents >= 3)
			{
				printVar("N_{3}", Nctau3);
				printVar("s_{3}", ctauTime3Scale);
				printVar("#Delta s_{32}", ctauTime3Delta);
			}
			if (nResolutionComponents >= 4)
			{
				printVar("N_{4}", Nctau4);
				printVar("s_{4}", ctauTime4Scale);
				printVar("#Delta s_{43}", ctauTime4Delta);
			}
		}

	cLifetime->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
	pad2->SetTopMargin(0.00001);
	pad2->SetBottomMargin(0.4);
	pad2->Draw();
	pad2->cd();

	RooPlot *fpull = obs_time.frame(Range(ctRange.first, ctRange.second), Title(""));
	RooHist *hpull = timefitplot->pullHist("data", "model");
	if (hpull)
		fpull->addPlotable(hpull, "P");
	fpull->GetYaxis()->SetTitle("Pull");
	fpull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} [mm]");
	fpull->GetXaxis()->CenterTitle();
	fpull->SetMinimum(-8);
	fpull->SetMaximum(8);
	fpull->GetYaxis()->SetNdivisions(505);
	fpull->GetYaxis()->SetTitleSize(0.12);
	fpull->GetYaxis()->SetLabelSize(0.10);
	fpull->GetXaxis()->SetTitleSize(0.15);
	fpull->GetXaxis()->SetLabelSize(0.10);
	fpull->Draw();

	// Count only the floating parameters of the ctau fit itself in chi2/ndf.
	auto chiM = chi2_from_pull(hpull);
	int npar = time_result ? time_result->floatParsFinal().getSize() : 0;
	int ndf = std::max(1, chiM.second - npar);
	double pvalue = TMath::Prob(chiM.first, ndf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", chiM.first, ndf));
	}

	TLine line(ctRange.first, 0.0, ctRange.second, 0.0);
	line.SetLineStyle(2);
	line.Draw("same");

	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

	// ------------------------------------------------------------------
	// print final result
	// ------------------------------------------------------------------
	std::cout << "------------------ FIT RESULT SUMMARY --------------------" << std::endl;
	std::cout << "Prompt resolution components used in PR fit: " << nResolutionComponents << std::endl;
	std::cout << "errPdfOpt in PR fit: " << errPdfOpt
			 << (errPdfOpt == kErrPdfHist ? " (RooHistPdf)" : " (analytic)") << std::endl;
	std::cout << "histPdfInterpolationOrder in PR fit: " << histPdfInterpolationOrder << std::endl;
	std::cout << "TimeErr components used in PR fit: G=" << nTimeErrGaussComponents
			 << " L=" << nTimeErrLandauComponents
			 << " LN=" << nTimeErrLognormalComponents << std::endl;
	std::cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << std::endl;
	if (timeErrResult)
	{
		std::cout << "------------------ FIT RESULT FOR TIME ERR ---------------" << std::endl;
		timeErrResult->Print("v");
	}
	else if (errPdfOpt == kErrPdfHist)
	{
		std::cout << "RooHistPdf template mode: no analytic err fit result." << std::endl;
	}
	time_result->Print();
	const TString figTimeErr = figName("timeerr_model");
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_pr err fit : " << figTimeErr << std::endl;
	std::cout << "[FIG] ctau_pr lifetime fit : " << figLifetime << std::endl;
	// time_pdf.Print("V");
}
