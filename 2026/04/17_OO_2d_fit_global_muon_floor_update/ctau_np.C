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
#include "RooMsgService.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
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
#include <array>
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

void ctau_np(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4, bool drawFromSavedFit = false, bool publish = false)
{
	ScopedMacroTimer timer("ctau_np", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	bool isSilence = true; // turn off tracing messages
	enum ErrPdfChoice
	{
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
	};
	int errPdfOpt = kErrPdfHist;
	const int histPdfInterpolationOrder = 1;
	double errPrefitDataFraction = 0.5;
	int nTimeErrGaussComponents = 2;
	int nTimeErrLandauComponents = 1;
	int nTimeErrLognormalComponents = 0;

	// Choose how many SingleSided signal components to use in each (pt, y) bin.
	int nSignalSSComponents = 2;
	if (yLow == 0.0f)
	{
		if (ptLow == 200.0f && ptHigh == 350.0f)
		{
			nSignalSSComponents = 2; // 1~3
		}
	}
	nTimeErrGaussComponents = std::clamp(nTimeErrGaussComponents, 0, 2);
	nTimeErrLandauComponents = std::clamp(nTimeErrLandauComponents, 0, 2);
	nTimeErrLognormalComponents = std::clamp(nTimeErrLognormalComponents, 0, 1);
	if (nTimeErrGaussComponents + nTimeErrLandauComponents + nTimeErrLognormalComponents <= 0)
	{
		std::cerr << "ERROR: at least one NP timeErr component is required." << std::endl;
		return;
	}
	const bool useGaus1 = (nTimeErrGaussComponents >= 1);
	const bool useGaus2 = (nTimeErrGaussComponents >= 2);
	const bool useLandau1 = (nTimeErrLandauComponents >= 1);
	const int nTimeErrComponents =
			static_cast<int>(useGaus1) + static_cast<int>(useGaus2) + static_cast<int>(useLandau1);
	nSignalSSComponents = std::clamp(nSignalSSComponents, 1, 3);
	if (isSilence)
	{
		RooMsgService::instance().getStream(0).removeTopic(RooFit::Tracing);
		RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
		RooMsgService::instance().getStream(1).removeTopic(RooFit::Tracing);
		RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
	}

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	// First we open the actual RooDataSet used in the prompt ctau analysis.
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_NP_globalOn_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	// ------------------------------------------------------------------
	// determine fit ranges
	// ------------------------------------------------------------------
	const auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.0001, 0.999, false);
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
	const TString figDir = TString::Format("figs/%s/ctau_np", yTag.Data());
	const TString prResultDir = TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString resultDir = TString::Format("roots/%s/ctau_np", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString npModelFileName = TString::Format("%s/ctau_np_model_%s.root", resultDir.Data(), figTag.Data());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit && !load_saved_fit_file(savedFitFile, npModelFileName, "nonprompt ctau"))
		return;
	if (drawFromSavedFit)
		nSignalSSComponents = std::clamp(read_saved_int_param(savedFitFile.get(), "nSignalSSComponents", nSignalSSComponents), 1, 3);

	auto *obs_time_scan = static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	auto *obs_timeErr_scan = static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	if (!obs_time_scan || !obs_timeErr_scan)
	{
		std::cerr << "ERROR: required variables ctau3D/ctau3DErr are missing before trimming." << std::endl;
		return;
	}
	obs_time_scan->setRange(ctRange.first, ctRange.second);
	obs_timeErr_scan->setRange(errRange.first, errRange.second);
	obs_timeErr_scan->setMin(errRange.first);
	obs_timeErr_scan->setMax(errRange.second);
	const int timePlotBins = std::max(2, obs_time_scan->getBins() * 2);
	int timeErrPlotBins = std::max(2, obs_timeErr_scan->getBins() * 2);
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", prResultDir.Data(), figTag.Data());

	int timeErrTrimBins = timeErrPlotBins;
	auto timeErrScanData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSel->reduce(RooArgSet(*obs_timeErr_scan))));
	auto timeErrScanHist = std::unique_ptr<TH1>(timeErrScanData ? timeErrScanData->createHistogram(
			"timeErrScanHist", *obs_timeErr_scan, Binning(timeErrTrimBins)) : nullptr);
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

	auto *massVar = static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	auto *timeVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	auto *timeErrVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	if (!massVar || !timeVar || !timeErrVar)
	{
		std::cerr << "ERROR: required observables mass/ctau3D/ctau3DErr are missing in selected dataset." << std::endl;
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

	// Keep the decay constants away from the numerically unstable tau -> 0 limit.
	const double absoluteLifetimeFloor = 1e-2;
	double lifetimeFloorCoeff = 0.25;
	if (yLow == 1.6f && ptHigh <= 4.0f)
		lifetimeFloorCoeff = 0.05;
	double signalLifetimeFloor = std::max(1.0 * errRange.first, absoluteLifetimeFloor);
	double bkgLifetimeFloor = std::max(1.0 * errRange.first, absoluteLifetimeFloor);
	double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	RooDataSet *data = dataSel.get();

	// ------------------------------------------------------------------
	// build ctau-error model
	// ------------------------------------------------------------------
	// Use an analytic core+tail model for the per-event error distribution.
	auto timeErrData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_timeErr))));
	const double errSpan = errRange.second - errRange.first;
	const double gaus1MeanInit = errRange.first + 0.14 * errSpan;
	const double gaus2MeanInit = errRange.first + 0.24 * errSpan;
	const double tailMpvInit = errRange.first + 0.40 * errSpan;
	const double gaus1SigmaInit = std::max(0.025 * errSpan, 5e-4);
	const double gaus2SigmaInit = std::max(0.055 * errSpan, 1e-3);
	const double tailWidthInit = std::max(0.10 * errSpan, 1e-3);
	RooRealVar timeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", gaus1MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", gaus1SigmaInit, 5e-4, std::max(0.15 * errSpan, 2e-3));
	RooRealVar timeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", gaus2MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", gaus2SigmaInit, 1e-3, std::max(0.30 * errSpan, 4e-3));
	RooRealVar timeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", tailMpvInit, errRange.first, errRange.second);
	RooRealVar timeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", tailWidthInit, 1e-3, std::max(0.50 * errSpan, 5e-3));
	RooRealVar timeErrCore1Frac("timeErrCore1Frac", "timeErrCore1Frac", 0.60, 0.05, 0.95);
	RooRealVar timeErrTailFrac("timeErrTailFrac", "timeErrTailFrac", 0.08, 0.001, 0.30);
	RooGaussian timeErrGaus1("timeErrGaus1", "timeErrGaus1", obs_timeErr, timeErrGaus1Mean, timeErrGaus1Sigma);
	RooGaussian timeErrGaus2("timeErrGaus2", "timeErrGaus2", obs_timeErr, timeErrGaus2Mean, timeErrGaus2Sigma);
	RooLandau timeErrTail("timeErrTail", "timeErrTail", obs_timeErr, timeErrTailMpv, timeErrTailWidth);
	RooArgList timeErrPdfList;
	RooArgList timeErrFracList;
	if (useLandau1)
	{
		timeErrPdfList.add(timeErrTail);
		if (nTimeErrComponents > 1)
			timeErrFracList.add(timeErrTailFrac);
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
				"ctau_np_hTimeErr", obs_timeErr, Binning(timeErrPlotBins, errRange.first, errRange.second)));
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
				std::cerr << "ERROR: timeErrResult not found in saved nonprompt ctau file: " << npModelFileName << std::endl;
				return;
			}
			apply_saved_fit_result(timeErrResult, *timeErrPdfAnalytic, RooArgSet(obs_timeErr));
		}
		else
			timeErrResult = timeErrPdfAnalytic->fitTo(
					*timeErrData,
					Save(true),
					PrintLevel(-1),
					SumW2Error(isWeight),
					PrefitDataFraction(errPrefitDataFraction));
		timeErrGaus1Mean.setConstant(true);
		timeErrGaus1Sigma.setConstant(true);
		timeErrGaus2Mean.setConstant(true);
		timeErrGaus2Sigma.setConstant(true);
		timeErrTailMpv.setConstant(true);
		timeErrTailWidth.setConstant(true);
		timeErrCore1Frac.setConstant(true);
		timeErrTailFrac.setConstant(true);
	}

	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	// ------------------------------------------------------------------
	// draw ctau-error model
	// ------------------------------------------------------------------
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
		tx.DrawLatex(xtext, y0 + dy * k++, "Nonprompt J/#psi MC");
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
		if (errPdfOpt == kErrPdfHist)
		{
			if (!publish)
				tx.DrawLatex(0.19, 0.765, Form("Status : RooHistPdf template (%d)", histPdfInterpolationOrder));
		}
		else
		{
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
		auto printVar = [&](const char *title, RooRealVar &var)
		{
			if (var.isConstant())
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), var.getError()));
		};
		if (errPdfOpt == kErrPdfHist)
		{
			tp.DrawLatex(xtext, y0 + dy * k++, Form("errPdfOpt = %d", errPdfOpt));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("interp = %d", histPdfInterpolationOrder));
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
	// load prompt resolution model from ctau_pr
	// ------------------------------------------------------------------
	TFile *resolutionFile = TFile::Open(resolutionFileName);
	if (!resolutionFile || resolutionFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open prompt resolution file: " << resolutionFileName << std::endl;
		return;
	}
	auto readDoubleParam = [&](const char *name) -> double
	{
		auto *param = dynamic_cast<TParameter<double> *>(resolutionFile->Get(name));
		if (!param)
		{
			std::cerr << "ERROR: missing double parameter '" << name << "' in " << resolutionFileName << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		return param->GetVal();
	};
	auto *nCompParam = dynamic_cast<TParameter<int> *>(resolutionFile->Get("nResolutionComponents"));
	if (!nCompParam)
	{
		std::cerr << "ERROR: missing parameter 'nResolutionComponents' in " << resolutionFileName << std::endl;
		return;
	}
	const int nResolutionComponents = nCompParam->GetVal();
	if (nResolutionComponents < 1 || nResolutionComponents > 4)
	{
		std::cerr << "ERROR: unsupported number of resolution components: " << nResolutionComponents << std::endl;
		return;
	}
	const double ctauMeanScaleVal = readDoubleParam("ctauMeanScale");
	const double ctauTime1MeanVal = readDoubleParam("ctauTime1Mean");
	const double ctauTime1ScaleVal = readDoubleParam("ctauTime1Scale");
	double ctauTime2MeanVal = 0.0;
	double ctauTime2ScaleVal = 0.0;
	double ctauFrac1Val = 0.0;
	double ctauTime3MeanVal = 0.0;
	double ctauTime3ScaleVal = 0.0;
	double ctauFrac2Val = 0.0;
	double ctauTime4MeanVal = 0.0;
	double ctauTime4ScaleVal = 0.0;
	double ctauFrac3Val = 0.0;
	if (nResolutionComponents >= 2)
	{
		ctauTime2MeanVal = readDoubleParam("ctauTime2Mean");
		ctauTime2ScaleVal = readDoubleParam("ctauTime2Scale");
		ctauFrac1Val = readDoubleParam("ctauFrac1");
	}
	if (nResolutionComponents >= 3)
	{
		ctauTime3MeanVal = readDoubleParam("ctauTime3Mean");
		ctauTime3ScaleVal = readDoubleParam("ctauTime3Scale");
		ctauFrac2Val = readDoubleParam("ctauFrac2");
	}
	if (nResolutionComponents >= 4)
	{
		ctauTime4MeanVal = readDoubleParam("ctauTime4Mean");
		ctauTime4ScaleVal = readDoubleParam("ctauTime4Scale");
		ctauFrac3Val = readDoubleParam("ctauFrac3");
	}
	if (!std::isfinite(ctauMeanScaleVal) || !std::isfinite(ctauTime1MeanVal) || !std::isfinite(ctauTime1ScaleVal) ||
			(nResolutionComponents >= 2 && (!std::isfinite(ctauTime2MeanVal) || !std::isfinite(ctauTime2ScaleVal) || !std::isfinite(ctauFrac1Val))) ||
			(nResolutionComponents >= 3 && (!std::isfinite(ctauTime3MeanVal) || !std::isfinite(ctauTime3ScaleVal) || !std::isfinite(ctauFrac2Val))) ||
			(nResolutionComponents >= 4 && (!std::isfinite(ctauTime4MeanVal) || !std::isfinite(ctauTime4ScaleVal) || !std::isfinite(ctauFrac3Val))))
	{
		std::cerr << "ERROR: invalid prompt resolution parameters in " << resolutionFileName << std::endl;
		return;
	}
	RooFitResult *resolutionResult = dynamic_cast<RooFitResult *>(resolutionFile->Get("timeResult"));
	auto readResolutionFitError = [&](const char *name) -> double
	{
		if (!resolutionResult)
			return std::numeric_limits<double>::quiet_NaN();
		auto *var = dynamic_cast<RooRealVar *>(resolutionResult->floatParsFinal().find(name));
		if (!var)
			var = dynamic_cast<RooRealVar *>(resolutionResult->constPars().find(name));
		return var ? var->getError() : std::numeric_limits<double>::quiet_NaN();
	};
	auto readResolutionYield = [&](const char *name) -> std::pair<double, double>
	{
		if (!resolutionResult)
			return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
		auto *var = dynamic_cast<RooRealVar *>(resolutionResult->floatParsFinal().find(name));
		if (!var)
			var = dynamic_cast<RooRealVar *>(resolutionResult->constPars().find(name));
		return var ? std::make_pair(var->getVal(), var->getError())
							 : std::make_pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	};
	auto resolutionFractionError = [&](int idx) -> double
	{
		std::array<std::pair<double, double>, 4> y = {
				readResolutionYield("Nctau1"),
				readResolutionYield("Nctau2"),
				readResolutionYield("Nctau3"),
				readResolutionYield("Nctau4")};
		double total = 0.0;
		for (int i = 0; i < nResolutionComponents; ++i)
		{
			if (!std::isfinite(y[i].first))
				return std::numeric_limits<double>::quiet_NaN();
			total += y[i].first;
		}
		if (!(total > 0.0))
			return std::numeric_limits<double>::quiet_NaN();
		double var = 0.0;
		for (int i = 0; i < nResolutionComponents; ++i)
		{
			if (!(std::isfinite(y[i].second) && y[i].second > 0.0))
				return std::numeric_limits<double>::quiet_NaN();
			const double deriv = (i == idx ? total - y[idx].first : -y[idx].first) / (total * total);
			var += deriv * deriv * y[i].second * y[i].second;
		}
		return std::sqrt(var);
	};
	resolutionFile->Close();
	std::cout << "Loaded ctau prompt resolution from " << resolutionFileName << std::endl;

	// ------------------------------------------------------------------
	// build signal ctau model
	// ------------------------------------------------------------------
	const double signalLifetimeInit = std::max(0.08, 1.5 * signalLifetimeFloor);
	const double signalLifetimeCeil = std::max(signalLifetimeInit * 5.0, maxStableLifetime);
	RooRealVar signal_lifetime(
			"signal_lifetime", "signal_lifetime",
			signalLifetimeInit,
			signalLifetimeFloor + 1e-2,
			signalLifetimeCeil);

	double maxPromptScaleSaved = std::max(ctauTime1ScaleVal, 0.0);
	if (nResolutionComponents >= 2 && std::isfinite(ctauTime2ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime2ScaleVal);
	if (nResolutionComponents >= 3 && std::isfinite(ctauTime3ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime3ScaleVal);
	if (nResolutionComponents >= 4 && std::isfinite(ctauTime4ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime4ScaleVal);

	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", ctauMeanScaleVal);
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", ctauTime1MeanVal);
	RooRealVar ctauTime1Scale(
			"ctauTime1Scale", "ctauTime1Scale",
			ctauTime1ScaleVal,
			std::max(0.3, 0.5 * ctauTime1ScaleVal),
			std::max(3.0, 2.0 * ctauTime1ScaleVal));
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", ctauTime2MeanVal);
	const double ctauTime2DeltaVal = std::max(0.05, ctauTime2ScaleVal - ctauTime1ScaleVal);
	RooRealVar ctauTime2Delta(
			"ctauTime2Delta", "ctauTime2Delta",
			ctauTime2DeltaVal,
			0.05,
			std::max(3.0, 2.0 * ctauTime2DeltaVal));
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac1("ctauFrac1", "ctauFrac1", ctauFrac1Val);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", ctauTime3MeanVal);
	const double ctauTime3DeltaVal = std::max(0.05, ctauTime3ScaleVal - ctauTime2ScaleVal);
	RooRealVar ctauTime3Delta(
			"ctauTime3Delta", "ctauTime3Delta",
			ctauTime3DeltaVal,
			0.05,
			std::max(3.0, 2.0 * ctauTime3DeltaVal));
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac2("ctauFrac2", "ctauFrac2", ctauFrac2Val);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", ctauTime4MeanVal);
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	RooRealVar ctauTime4Delta(
			"ctauTime4Delta", "ctauTime4Delta",
			ctauTime4DeltaVal,
			0.05,
			std::max(3.0, 2.0 * ctauTime4DeltaVal));
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac3("ctauFrac3", "ctauFrac3", ctauFrac3Val);
	ctauTime1Mean.setConstant(true);
	ctauTime2Mean.setConstant(true);
	ctauTime3Mean.setConstant(true);
	ctauTime4Mean.setConstant(true);
	RooArgSet timeConstraints;
	std::vector<std::unique_ptr<RooRealVar>> timeConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> timeConstraintPdfs;
	auto constraintSigmaFloor = [&](double central, double relFloor, double absFloor)
	{
		return std::max(absFloor, relFloor * std::abs(central));
	};
	auto addConstraint = [&](const char *baseName, RooRealVar &var, double central, double sigma)
	{
		if (!(std::isfinite(central) && std::isfinite(sigma) && sigma > 0.0))
			return;
		var.setVal(central);
		const TString meanName = TString::Format("%s_mean", baseName);
		const TString sigmaName = TString::Format("%s_sigma", baseName);
		const TString pdfName = TString::Format("%s_constraint", baseName);
		timeConstraintConsts.push_back(std::make_unique<RooRealVar>(meanName, meanName, central));
		timeConstraintConsts.back()->setConstant(true);
		timeConstraintConsts.push_back(std::make_unique<RooRealVar>(
				sigmaName, sigmaName, sigma, 1e-9, std::max(10.0 * sigma, 1e-8)));
		timeConstraintConsts.back()->setConstant(true);
		timeConstraintPdfs.push_back(std::make_unique<RooGaussian>(
				pdfName, pdfName, var,
				*timeConstraintConsts[timeConstraintConsts.size() - 2],
				*timeConstraintConsts.back()));
		timeConstraints.add(*timeConstraintPdfs.back());
	};
	auto constraintSigmaFromFit = [&](const char *fitName, double central, double relFloor, double absFloor)
	{
		double sigma = constraintSigmaFloor(central, relFloor, absFloor);
		const double fitErr = readResolutionFitError(fitName);
		if (std::isfinite(fitErr) && fitErr > 0.0)
			sigma = std::max(sigma, fitErr);
		return sigma;
	};
	auto fractionConstraintSigma = [&](int idx, double central)
	{
		double sigma = constraintSigmaFloor(central, 0.10, 0.02);
		const double fitErr = resolutionFractionError(idx);
		if (std::isfinite(fitErr) && fitErr > 0.0)
			sigma = std::max(sigma, fitErr);
		return sigma;
	};
	addConstraint("ctauTime1Scale", ctauTime1Scale, ctauTime1ScaleVal,
			constraintSigmaFromFit("ctauTime1Scale", ctauTime1ScaleVal, 0.10, 0.05));
	if (nResolutionComponents >= 2)
	{
		addConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal,
				constraintSigmaFromFit("ctauTime2Delta", ctauTime2DeltaVal, 0.10, 0.05));
		addConstraint("ctauFrac1", ctauFrac1, ctauFrac1Val,
				fractionConstraintSigma(0, ctauFrac1Val));
	}
	if (nResolutionComponents >= 3)
	{
		addConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal,
				constraintSigmaFromFit("ctauTime3Delta", ctauTime3DeltaVal, 0.10, 0.05));
		addConstraint("ctauFrac2", ctauFrac2, ctauFrac2Val,
				fractionConstraintSigma(1, ctauFrac2Val));
	}
	if (nResolutionComponents >= 4)
	{
		addConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal,
				constraintSigmaFromFit("ctauTime4Delta", ctauTime4DeltaVal, 0.10, 0.05));
		addConstraint("ctauFrac3", ctauFrac3, ctauFrac3Val,
				fractionConstraintSigma(2, ctauFrac3Val));
	}

	std::unique_ptr<RooAddModel> promptResolutionModel;
	RooResolutionModel *time_resolution_ptr = &ctauTime1;
	if (nResolutionComponents >= 2)
	{
		RooArgList resolutionPdfList;
		RooArgList resolutionFracList;
		resolutionPdfList.add(ctauTime1);
		resolutionPdfList.add(ctauTime2);
		resolutionFracList.add(ctauFrac1);
		if (nResolutionComponents >= 3)
		{
			resolutionPdfList.add(ctauTime3);
			resolutionFracList.add(ctauFrac2);
		}
		if (nResolutionComponents >= 4)
		{
			resolutionPdfList.add(ctauTime4);
			resolutionFracList.add(ctauFrac3);
		}
		promptResolutionModel = std::make_unique<RooAddModel>(
				"promptResolutionModel", "promptResolutionModel",
				resolutionPdfList, resolutionFracList);
		time_resolution_ptr = promptResolutionModel.get();
	}
	RooResolutionModel &time_resolution = *time_resolution_ptr;
	const double resolutionDrivenLifetimeFloor =
		std::max(absoluteLifetimeFloor, lifetimeFloorCoeff * maxPromptScaleSaved * std::max(errRange.second, absoluteLifetimeFloor));
	signalLifetimeFloor = std::max(signalLifetimeFloor, resolutionDrivenLifetimeFloor);
	bkgLifetimeFloor = std::max(bkgLifetimeFloor, resolutionDrivenLifetimeFloor);
	maxStableLifetime = std::max(maxStableLifetime, 20.0 * bkgLifetimeFloor);

	// Build up to three positive lifetime components that share the prompt resolution kernel.
	RooDecay signal_ss1_time("signal_ss1_time", "signal_ss1_time", obs_time, signal_lifetime, time_resolution, RooDecay::SingleSided);
	const double signalLifetime2Floor = std::max(0.75 * signalLifetimeFloor, 5e-3);
	const double signalLifetime2Init = std::max(0.20, 2.0 * signalLifetime2Floor);
	const double signalLifetime2Ceil = std::max(signalLifetime2Init * 5.0, maxStableLifetime);
	RooRealVar signal2_lifetime(
			"signal2_lifetime", "signal2_lifetime",
			signalLifetime2Init,
			signalLifetime2Floor + 2e-2,
			signalLifetime2Ceil);
	RooDecay signal_ss2_time("signal_ss2_time", "signal_ss2_time", obs_time, signal2_lifetime, time_resolution, RooDecay::SingleSided);
	const double signalLifetime3Floor = std::max(1.00 * signalLifetimeFloor, 5e-3);
	const double signalLifetime3Init = std::max(0.40, 2.0 * signalLifetime3Floor);
	const double signalLifetime3Ceil = std::max(signalLifetime3Init * 5.0, maxStableLifetime);
	RooRealVar signal3_lifetime(
			"signal3_lifetime", "signal3_lifetime",
			signalLifetime3Init,
			signalLifetime3Floor + 1e-2,
			signalLifetime3Ceil);
	RooDecay signal_ss3_time("signal_ss3_time", "signal_ss3_time", obs_time, signal3_lifetime, time_resolution, RooDecay::SingleSided);
	RooRealVar NsignalSS1("NsignalSS1", "NsignalSS1", 0.5 * data->numEntries(), 0.0 * data->numEntries(), std::max(1, data->numEntries()));
	RooRealVar NsignalSS2("NsignalSS2", "NsignalSS2", 0.3 * data->numEntries(), 0.0, std::max(1, data->numEntries()));
	RooRealVar NsignalSS3("NsignalSS3", "NsignalSS3", 0.2 * data->numEntries(), 0.0, std::max(1, data->numEntries()));
	auto setSeedIfFinite = [](RooRealVar &var, double value)
	{
		if (!std::isfinite(value))
			return;
		var.setVal(std::clamp(value, var.getMin(), var.getMax()));
	};
	auto applyLifetimeSeedDefaults = [&]()
	{
		double tau1Seed = 0.30;
		double tau2Seed = 0.46;
		double tau3Seed = 0.80;
		double frac1Seed = 0.60;
		double frac2Seed = 0.25;
		if (yLow == 0.0f)
		{
			tau1Seed = 0.34;
			tau2Seed = 0.48;
			tau3Seed = 0.90;
			frac1Seed = 0.62;
			frac2Seed = 0.24;
		}
		else if (ptHigh <= 2.0f)
		{
			tau1Seed = 0.05;
			tau2Seed = 0.45;
			tau3Seed = 0.90;
			frac1Seed = 0.25;
			frac2Seed = 0.20;
		}
		else if (ptHigh <= 6.0f)
		{
			tau1Seed = 0.18;
			tau2Seed = 0.42;
			tau3Seed = 0.85;
			frac1Seed = 0.45;
			frac2Seed = 0.25;
		}
		setSeedIfFinite(signal_lifetime, tau1Seed);
		if (nSignalSSComponents >= 2)
			setSeedIfFinite(signal2_lifetime, tau2Seed);
		if (nSignalSSComponents >= 3)
			setSeedIfFinite(signal3_lifetime, tau3Seed);
		const double totalEntries = std::max(1.0, static_cast<double>(data->numEntries()));
		NsignalSS1.setVal(std::clamp(frac1Seed * totalEntries, NsignalSS1.getMin(), NsignalSS1.getMax()));
		if (nSignalSSComponents >= 2)
		{
			const double frac2Norm = (nSignalSSComponents >= 3) ? frac2Seed : (1.0 - frac1Seed);
			NsignalSS2.setVal(std::clamp(frac2Norm * totalEntries, NsignalSS2.getMin(), NsignalSS2.getMax()));
		}
		if (nSignalSSComponents >= 3)
		{
			const double frac3Seed = std::max(0.05, 1.0 - frac1Seed - frac2Seed);
			NsignalSS3.setVal(std::clamp(frac3Seed * totalEntries, NsignalSS3.getMin(), NsignalSS3.getMax()));
		}
	};
	auto enforceOrderedLifetimeSeeds = [&]()
	{
		const double tauGap12 = 0.05;
		const double tauGap23 = 0.08;
		setSeedIfFinite(signal_lifetime, signal_lifetime.getVal());
		if (nSignalSSComponents >= 2)
		{
			const double tau2Min = std::min(signal2_lifetime.getMax(), std::max(signal2_lifetime.getMin(), signal_lifetime.getVal() + tauGap12));
			setSeedIfFinite(signal2_lifetime, std::max(signal2_lifetime.getVal(), tau2Min));
		}
		if (nSignalSSComponents >= 3)
		{
			const double tau3Min = std::min(signal3_lifetime.getMax(), std::max(signal3_lifetime.getMin(), signal2_lifetime.getVal() + tauGap23));
			setSeedIfFinite(signal3_lifetime, std::max(signal3_lifetime.getVal(), tau3Min));
		}
	};
	applyLifetimeSeedDefaults();
	enforceOrderedLifetimeSeeds();
	std::unique_ptr<RooAddPdf> signal_time;
		if (nSignalSSComponents == 1)
		{
			NsignalSS2.setVal(0.0);
			NsignalSS2.setConstant(true);
			NsignalSS3.setVal(0.0);
			NsignalSS3.setConstant(true);
			signal_time = std::make_unique<RooAddPdf>(
					"signal_time", "signal_time",
					RooArgList(signal_ss1_time),
					RooArgList(NsignalSS1));
		}
	else if (nSignalSSComponents == 2)
	{
		NsignalSS3.setVal(0.0);
		NsignalSS3.setConstant(true);
		signal_time = std::make_unique<RooAddPdf>(
				"signal_time", "signal_time",
				RooArgList(signal_ss1_time, signal_ss2_time),
				RooArgList(NsignalSS1, NsignalSS2));
	}
	else
	{
		signal_time = std::make_unique<RooAddPdf>(
				"signal_time", "signal_time",
				RooArgList(signal_ss1_time, signal_ss2_time, signal_ss3_time),
				RooArgList(NsignalSS1, NsignalSS2, NsignalSS3));
	}

	// ------------------------------------------------------------------
	// fit the signal ctau model
	// ------------------------------------------------------------------
	RooProdPdf time_pdf("time_pdf", "time_pdf", RooArgSet(*timeErrPdf), Conditional(RooArgSet(*signal_time), RooArgSet(obs_time)));
	const double ctauTime1MeanSeed = ctauTime1Mean.getVal();
	const double ctauTime1ScaleSeed = ctauTime1Scale.getVal();
	const double ctauTime2DeltaSeed = ctauTime2Delta.getVal();
	const double ctauTime3DeltaSeed = ctauTime3Delta.getVal();
	const double ctauTime4DeltaSeed = ctauTime4Delta.getVal();
	const double ctauFrac1Seed = ctauFrac1.getVal();
	const double ctauFrac2Seed = ctauFrac2.getVal();
	const double ctauFrac3Seed = ctauFrac3.getVal();
	const double signalLifetimeSeed = signal_lifetime.getVal();
	const double signalLifetime2Seed = signal2_lifetime.getVal();
	const double signalLifetime3Seed = signal3_lifetime.getVal();
	const double signalYield1Seed = NsignalSS1.getVal();
	const double signalYield2Seed = NsignalSS2.getVal();
	const double signalYield3Seed = NsignalSS3.getVal();
	auto resetTimeFitSeeds = [&]()
	{
		setSeedIfFinite(ctauTime1Mean, ctauTime1MeanSeed);
		setSeedIfFinite(ctauTime1Scale, ctauTime1ScaleSeed);
		if (nResolutionComponents >= 2)
		{
			setSeedIfFinite(ctauTime2Delta, ctauTime2DeltaSeed);
			setSeedIfFinite(ctauFrac1, ctauFrac1Seed);
		}
		if (nResolutionComponents >= 3)
		{
			setSeedIfFinite(ctauTime3Delta, ctauTime3DeltaSeed);
			setSeedIfFinite(ctauFrac2, ctauFrac2Seed);
		}
		if (nResolutionComponents >= 4)
		{
			setSeedIfFinite(ctauTime4Delta, ctauTime4DeltaSeed);
			setSeedIfFinite(ctauFrac3, ctauFrac3Seed);
		}
		setSeedIfFinite(signal_lifetime, signalLifetimeSeed);
		setSeedIfFinite(NsignalSS1, signalYield1Seed);
		if (nSignalSSComponents >= 2)
		{
			setSeedIfFinite(signal2_lifetime, signalLifetime2Seed);
			if (!NsignalSS2.isConstant())
				setSeedIfFinite(NsignalSS2, signalYield2Seed);
		}
		if (nSignalSSComponents >= 3)
		{
			setSeedIfFinite(signal3_lifetime, signalLifetime3Seed);
			if (!NsignalSS3.isConstant())
				setSeedIfFinite(NsignalSS3, signalYield3Seed);
		}
		enforceOrderedLifetimeSeeds();
	};
	RooFitResult *time_result = nullptr;
	std::unique_ptr<RooFitResult> savedTimeResult;
	if (drawFromSavedFit)
	{
		savedTimeResult = clone_saved_fit_result(savedFitFile.get(), "timeResult");
		time_result = savedTimeResult.get();
		if (time_result)
		{
			apply_saved_fit_result(time_result, time_pdf, RooArgSet(obs_time, obs_timeErr));
		}
		else
		{
			auto applySavedValue = [&](const char *name, RooRealVar &var)
			{
				const double value = read_saved_double_param(savedFitFile.get(), name, std::numeric_limits<double>::quiet_NaN());
				if (std::isfinite(value))
					var.setVal(value);
			};
			applySavedValue("signal_lifetime", signal_lifetime);
			if (nSignalSSComponents >= 2)
				applySavedValue("signal2_lifetime", signal2_lifetime);
			if (nSignalSSComponents >= 3)
				applySavedValue("signal3_lifetime", signal3_lifetime);
			applySavedValue("NsignalSS1", NsignalSS1);
			if (nSignalSSComponents >= 2)
				applySavedValue("NsignalSS2", NsignalSS2);
			if (nSignalSSComponents >= 3)
				applySavedValue("NsignalSS3", NsignalSS3);
			std::cout << "[PlotOnly] timeResult not found; loaded saved nonprompt ctau TParameters from "
								<< npModelFileName << std::endl;
		}
	}
	else
	{
		struct TimeFitAttempt
		{
			int strategy;
			double prefitFraction;
			const char *label;
		};
		const std::array<TimeFitAttempt, 3> attempts = {{
			{1, errPrefitDataFraction, "nominal"},
			{1, std::max(errPrefitDataFraction, 0.70), "prefit70"},
			{2, std::max(errPrefitDataFraction, 0.85), "strategy2"}
		}};
		for (size_t i = 0; i < attempts.size(); ++i)
		{
			const auto &attempt = attempts[i];
			if (i > 0)
				resetTimeFitSeeds();
			if (time_result)
			{
				delete time_result;
				time_result = nullptr;
			}
			time_result = time_pdf.fitTo(
					*data,
					Strategy(attempt.strategy),
					Extended(),
					Save(),
					SumW2Error(isWeight),
					Offset(true),
					ExternalConstraints(timeConstraints),
					PrefitDataFraction(attempt.prefitFraction),
					RecoverFromUndefinedRegions(1.0));
			const bool fitOk = time_result && time_result->status() == 0 && time_result->covQual() >= 2;
			if (fitOk)
				break;
			if (time_result && i + 1 < attempts.size())
			{
				std::cout << "[WARN] np time fit attempt '" << attempt.label
						  << "' ended with status=" << time_result->status()
						  << " covQual=" << time_result->covQual()
						  << ", retrying." << std::endl;
			}
		}
	}
	if (time_result)
		time_result->Print();

	if (!drawFromSavedFit)
	{
		TFile modelFile(npModelFileName, "RECREATE");
		if (!modelFile.IsZombie())
		{
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
			TParameter<double>("ctauMeanScale", ctauMeanScale.getVal()).Write();
			TParameter<double>("ctauTime1Mean", ctauTime1Mean.getVal()).Write();
			TParameter<double>("ctauTime1Scale", ctauTime1Scale.getVal()).Write();
			if (nResolutionComponents >= 2)
			{
				TParameter<double>("ctauTime2Mean", ctauTime1Mean.getVal()).Write();
				TParameter<double>("ctauTime2Scale", ctauTime2Scale.getVal()).Write();
				TParameter<double>("ctauFrac1", ctauFrac1.getVal()).Write();
			}
			if (nResolutionComponents >= 3)
			{
				TParameter<double>("ctauTime3Mean", ctauTime1Mean.getVal()).Write();
				TParameter<double>("ctauTime3Scale", ctauTime3Scale.getVal()).Write();
				TParameter<double>("ctauFrac2", ctauFrac2.getVal()).Write();
			}
			if (nResolutionComponents >= 4)
			{
				TParameter<double>("ctauTime4Mean", ctauTime1Mean.getVal()).Write();
				TParameter<double>("ctauTime4Scale", ctauTime4Scale.getVal()).Write();
				TParameter<double>("ctauFrac3", ctauFrac3.getVal()).Write();
			}
			TParameter<double>("signal_lifetime", signal_lifetime.getVal()).Write();
			TParameter<double>("absoluteLifetimeFloor", absoluteLifetimeFloor).Write();
			TParameter<double>("lifetimeFloorCoeff", lifetimeFloorCoeff).Write();
			TParameter<double>("maxPromptScaleSaved", maxPromptScaleSaved).Write();
			TParameter<double>("resolutionDrivenLifetimeFloor", resolutionDrivenLifetimeFloor).Write();
			TParameter<double>("signalLifetimeFloor", signalLifetimeFloor).Write();
			TParameter<double>("signalLifetimeFloorApplied", signalLifetimeFloor).Write();
			TParameter<double>("NsignalSS1", NsignalSS1.getVal()).Write();
			if (nSignalSSComponents >= 2)
			{
				TParameter<double>("signal2_lifetime", signal2_lifetime.getVal()).Write();
				TParameter<double>("signalLifetime2Floor", signalLifetime2Floor).Write();
				TParameter<double>("NsignalSS2", NsignalSS2.getVal()).Write();
			}
			if (nSignalSSComponents >= 3)
			{
				TParameter<double>("signal3_lifetime", signal3_lifetime.getVal()).Write();
				TParameter<double>("signalLifetime3Floor", signalLifetime3Floor).Write();
				TParameter<double>("NsignalSS3", NsignalSS3.getVal()).Write();
			}
			modelFile.Write();
		}
		else
		{
			std::cerr << "ERROR: cannot create NP model file: " << npModelFileName << std::endl;
		}
	}
	if (drawFromSavedFit)
		std::cout << "[PlotOnly] Loaded saved nonprompt ctau fit and left ROOT file unchanged: " << npModelFileName << std::endl;
	else
		std::cout << "Saved ctau NP model parameters to " << npModelFileName << std::endl;

	// Build component-only PDFs for plotting with the same conditional structure.
	std::unique_ptr<RooProdPdf> time_pdf_ss1;
	std::unique_ptr<RooProdPdf> time_pdf_ss2;
	std::unique_ptr<RooProdPdf> time_pdf_ss3;
	time_pdf_ss1 = std::make_unique<RooProdPdf>(
			"time_pdf_ss1",
			"time_pdf_ss1",
			RooArgSet(*timeErrPdf),
			Conditional(RooArgSet(signal_ss1_time), RooArgSet(obs_time)));
	if (nSignalSSComponents >= 2)
	{
		time_pdf_ss2 = std::make_unique<RooProdPdf>(
				"time_pdf_ss2",
				"time_pdf_ss2",
				RooArgSet(*timeErrPdf),
				Conditional(RooArgSet(signal_ss2_time), RooArgSet(obs_time)));
	}
	if (nSignalSSComponents == 3)
	{
		time_pdf_ss3 = std::make_unique<RooProdPdf>(
				"time_pdf_ss3",
				"time_pdf_ss3",
				RooArgSet(*timeErrPdf),
				Conditional(RooArgSet(signal_ss3_time), RooArgSet(obs_time)));
	}

	// ------------------------------------------------------------------
	// draw ctau fit
	// ------------------------------------------------------------------
	TCanvas *cLifetime = new TCanvas("cLifetime", "cLifetime", 800, 800);
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
	time_pdf_ss1->plotOn(
			timefitplot,
			LineColor(kBlue + 1),
			LineStyle(kDashed),
			LineWidth(2),
			Normalization(NsignalSS1.getVal(), RooAbsReal::NumEvent),
			Name("ss1_component"));
	if (time_pdf_ss2)
	{
		time_pdf_ss2->plotOn(
				timefitplot,
				LineColor(kRed + 1),
				LineStyle(kDashed),
				LineWidth(2),
				Normalization(NsignalSS2.getVal(), RooAbsReal::NumEvent),
				Name("ss2_component"));
	}
	if (time_pdf_ss3)
	{
		time_pdf_ss3->plotOn(
				timefitplot,
				LineColor(kMagenta + 1),
				LineStyle(kDashed),
				LineWidth(2),
				Normalization(NsignalSS3.getVal(), RooAbsReal::NumEvent),
				Name("ss3_component"));
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
	if (auto *o = findObj(timefitplot, "ss1_component"))
		leg.AddEntry(o, "SS1", "l");
	if (time_pdf_ss2)
		if (auto *o = findObj(timefitplot, "ss2_component"))
			leg.AddEntry(o, "SS2", "l");
	if (time_pdf_ss3)
		if (auto *o = findObj(timefitplot, "ss3_component"))
			leg.AddEntry(o, "SS3", "l");
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
		tx.DrawLatex(xtext, y0 + dy * k++, "Nonprompt J/#psi MC");
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
		auto printVar = [&](const char *title, const RooAbsReal &var)
		{
			const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
			if (rrv && rrv->isConstant())
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
			else if (rrv)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
		};
		auto printFloatingVar = [&](const char *title, const RooAbsReal &var)
		{
			const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
			if (rrv && rrv->isConstant())
				return;
			printVar(title, var);
		};
			printVar("N_{SS1}", NsignalSS1);
			if (nSignalSSComponents >= 2)
				printVar("N_{SS2}", NsignalSS2);
			if (nSignalSSComponents == 3)
				printVar("N_{SS3}", NsignalSS3);
			printVar("#tau_{SS1}", signal_lifetime);
		if (nSignalSSComponents >= 2)
			printVar("#tau_{SS2}", signal2_lifetime);
		if (nSignalSSComponents == 3)
			printVar("#tau_{SS3}", signal3_lifetime);
			printFloatingVar("#sigma_{1}^{res}", ctauTime1Scale);
			if (nResolutionComponents >= 2)
				printFloatingVar("#Delta s_{21}^{res}", ctauTime2Delta);
			if (nResolutionComponents >= 3)
				printFloatingVar("#Delta s_{32}^{res}", ctauTime3Delta);
			if (nResolutionComponents >= 4)
				printFloatingVar("#Delta s_{43}^{res}", ctauTime4Delta);
			if (nResolutionComponents >= 2)
				printFloatingVar("f_{res1}", ctauFrac1);
			if (nResolutionComponents >= 3)
				printFloatingVar("f_{res2}", ctauFrac2);
			if (nResolutionComponents >= 4)
				printFloatingVar("f_{res3}", ctauFrac3);
	}

	// ------------------------------------------------------------------
	// draw pull and print summary
	// ------------------------------------------------------------------
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

	cout << "------------------ FIT RESULT SUMMARY --------------------" << endl;
	cout << "errPdfOpt in NP fit: " << errPdfOpt
			 << (errPdfOpt == kErrPdfHist ? " (RooHistPdf)" : " (analytic)") << endl;
	cout << "histPdfInterpolationOrder in NP fit: " << histPdfInterpolationOrder << endl;
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	if (timeErrResult)
	{
		cout << "------------------ FIT RESULT FOR TIME ERR ---------------" << endl;
		timeErrResult->Print("v");
	}
	if (time_result)
		time_result->Print("v");
	{
		auto printFloorVsValue2Sigma = [&](const char *name, RooAbsReal &var, double floor)
		{
			double err = 0.0;
			bool hasErr = false;
			if (auto *rrv = dynamic_cast<RooRealVar *>(&var))
			{
				err = rrv->getError();
				hasErr = std::isfinite(err) && err > 0.0;
			}
			else if (time_result)
			{
				err = var.getPropagatedError(*time_result);
				hasErr = std::isfinite(err) && err > 0.0;
			}
			const double value2Sigma = var.getVal() + (hasErr ? 2.0 * err : 0.0);
			std::cout << Form("[FLOOR] %s : floor=%.6g , value+2sigma=%.6g", name, floor, value2Sigma);
			if (std::isfinite(floor))
				std::cout << Form(" , margin=%.6g", value2Sigma - floor);
			std::cout << std::endl;
		};
		std::cout << "---------------- FLOOR VS VALUE+2SIGMA -------------" << std::endl;
		printFloorVsValue2Sigma("signal_lifetime", signal_lifetime, signalLifetimeFloor);
		if (nSignalSSComponents >= 2)
			printFloorVsValue2Sigma("signal2_lifetime", signal2_lifetime, signalLifetime2Floor);
		if (nSignalSSComponents >= 3)
			printFloorVsValue2Sigma("signal3_lifetime", signal3_lifetime, signalLifetime3Floor);
	}
	cout << "Prompt resolution components used in NP fit: " << nResolutionComponents << endl;
	cout << "SingleSided components used in NP fit: " << nSignalSSComponents << endl;
	cout << "TimeErr components used in NP fit: G=" << nTimeErrGaussComponents
			 << " L=" << nTimeErrLandauComponents
			 << " LN=" << nTimeErrLognormalComponents << endl;
	const TString figTimeErr = figName("timeerr_model");
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_np err fit : " << figTimeErr << std::endl;
	std::cout << "[FIG] ctau_np lifetime fit : " << figLifetime << std::endl;
}
