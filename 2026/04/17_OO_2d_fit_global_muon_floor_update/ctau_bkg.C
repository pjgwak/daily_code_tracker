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
#include "RooLognormal.h"
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

void ctau_bkg(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4, bool isSilence = true, bool drawFromSavedFit = false, bool publish = false)
{
	ScopedMacroTimer timer("ctau_bkg", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;
	double errPrefitDataFraction = 0.5;
	double bkgRetryPrefitDataFraction = 0.8;
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	enum ErrPdfChoice
	{
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
	};
	const int overrideErrPdfOpt = -1;
	const int histPdfInterpolationOrder = 1;
	int errPdfOpt = kErrPdfAnalytic;
	int nTimeErrGaussComponents = 2;
	int nTimeErrLandauComponents = 2;
	int nTimeErrLognormalComponents = 1;

	// Choose how many background components to use in each (pt, y) bin.
	int nSignalSSComponents = 1;
	int nSignalFlipComponents = 1;
	int nSignalDSComponents = 1;
	if (yLow == 0.0f)
	{
		if (ptLow >= 10.0f)
		{
			nSignalSSComponents = 1; // 0~3
			nSignalFlipComponents = 0; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
	}
	else if (yLow == 1.6f)
	{
		if (ptLow == 1.0f && ptHigh == 2.0f)
		{
			nSignalSSComponents = 2; // 0~3
			nSignalFlipComponents = 1; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
		if (ptLow == 2.0f && ptHigh == 3.0f)
		{
			nSignalSSComponents = 2; // 0~3
			nSignalFlipComponents = 1; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
		if (ptLow == 3.0f && ptHigh == 4.0f)
		{
			nSignalSSComponents = 2; // 0~3
			nSignalFlipComponents = 1; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
		// if (ptLow == 8.0f && ptHigh == 9.0f)
		// {
		// 	nSignalSSComponents = 1; // 0~3
		// 	nSignalFlipComponents = 1; // 0~3
		// 	nSignalDSComponents = 1; // 0~3
		// }
		// if (ptLow == 14.0f && ptHigh == 20.0f)
		// {
		// 	nSignalSSComponents = 1; // 0~3
		// 	nSignalFlipComponents = 0; // 0~3
		// 	nSignalDSComponents = 0; // 0~3
		// 	nTimeErrGaussComponents = 2;
		// 	nTimeErrLandauComponents = 2;
		// 	nTimeErrLognormalComponents = 1;
		// }
	}
	nSignalSSComponents = std::clamp(nSignalSSComponents, 0, 3);
	nSignalFlipComponents = std::clamp(nSignalFlipComponents, 0, 3);
	nSignalDSComponents = std::clamp(nSignalDSComponents, 0, 3);
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
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOn_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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
	auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.001, 0.999, false);
	auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.995, true);
	if (errRange.first < 1e-6)
		errRange.first = 1e-6;

	if (yLow == 0.0f)
	{
		if (ptLow >= 10.0f)
		{
			// dummy
		}
	}
	else if (yLow == 1.6f)
	{
		if (ptLow == 1.0f && ptHigh == 2.0f)
		{
			ctRange.first = -4;
			ctRange.second = 6;
		}
		if (ptLow == 2.0f && ptHigh == 3.0f)
		{
			ctRange.first = -3;
			ctRange.second = 3;
		}
	}

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

	// ------------------------------------------------------------------
	// output and observable setup
	// ------------------------------------------------------------------
	auto formatTag = [](double value)
	{
		return TString::Format("%.2f", value);
	};
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString figDir = TString::Format("figs/%s/ctau_bkg", yTag.Data());
	const TString prResultDir = TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString bkgResultDir = TString::Format("roots/%s/ctau_bkg", yTag.Data());
	const TString errFileName = TString::Format("roots/%s/err2/err2_model_%s.root", yTag.Data(), figTag.Data());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(bkgResultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

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
	const int timePlotBins = std::max(2, obs_time.getBins() * 2);
	const int timeErrPlotBins = std::max(2, obs_timeErr.getBins() * 2);
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", prResultDir.Data(), figTag.Data());
	const TString fitResultFileName = TString::Format("%s/ctau_bkg_fitresult_%s.root", bkgResultDir.Data(), figTag.Data());
	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit && !load_saved_fit_file(savedFitFile, fitResultFileName, "background ctau"))
		return;

	const double absoluteLifetimeFloor = 1e-2;
	double lifetimeFloorCoeff = 0.25;
	if (yLow == 1.6f && ptHigh <= 4.0f)
		lifetimeFloorCoeff = 0.05;
	double bkgLifetimeFloor = std::max(1.0 * errRange.first, absoluteLifetimeFloor);
	double bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 2e-2);
	const double bkgSymLifetimeFloor = std::max(1.5 * errRange.first, 1e-2);
	double bkgFlipLifetimeFloor = std::max(0.75 * bkgLifetimeFloor, 5e-3);
	double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	RooDataSet *data = dataSB.get();

	// ------------------------------------------------------------------
	// build ctau-error model from err2.C outputs
	// ------------------------------------------------------------------
	auto timeErrData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_timeErr))));
	auto openFile = [](const TString &path) -> TFile * {
		TFile *f = TFile::Open(path);
		if (!f || f->IsZombie())
		{
			std::cerr << "ERROR: cannot open " << path << std::endl;
			return nullptr;
		}
		return f;
	};
	auto readSavedIntParam = [](TFile &f, const char *name, int fallback = -999) {
		auto *param = dynamic_cast<TParameter<int> *>(f.Get(name));
		return param ? param->GetVal() : fallback;
	};
	auto readSavedDoubleParam = [](TFile &f, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN()) {
		auto *param = dynamic_cast<TParameter<double> *>(f.Get(name));
		return param ? param->GetVal() : fallback;
	};
	auto readSavedErrValue = [](TFile &f, RooFitResult *fr, const char *name, double fallback) {
		if (auto *param = dynamic_cast<TParameter<double> *>(f.Get(name)))
			return param->GetVal();
		if (auto *var = dynamic_cast<RooRealVar *>(f.Get(name)))
			return var->getVal();
		if (auto *abs = dynamic_cast<RooAbsReal *>(f.Get(name)))
			return abs->getVal();
		if (fr)
		{
			if (auto *var = dynamic_cast<RooRealVar *>(fr->floatParsFinal().find(name)))
				return var->getVal();
			if (auto *var = dynamic_cast<RooRealVar *>(fr->constPars().find(name)))
				return var->getVal();
		}
		return fallback;
	};
	std::unique_ptr<TFile> errFile(openFile(errFileName));
	if (!errFile)
		return;
	const double errFitLow = readSavedDoubleParam(*errFile, "errLow", errRange.first);
	const double errFitHigh = readSavedDoubleParam(*errFile, "errHigh", errRange.second);
	obs_timeErr.setRange(errFitLow, errFitHigh);
	obs_timeErr.setMin(errFitLow);
	errPdfOpt = overrideErrPdfOpt >= 0 ? overrideErrPdfOpt : readSavedIntParam(*errFile, "bkgErrPdfOpt", kErrPdfAnalytic);
	nTimeErrGaussComponents = std::clamp(readSavedIntParam(*errFile, "nBkgTimeErrGaussComponents", 1), 0, 2);
	nTimeErrLandauComponents = std::clamp(readSavedIntParam(*errFile, "nBkgTimeErrLandauComponents", 1), 0, 2);
	nTimeErrLognormalComponents = std::clamp(readSavedIntParam(*errFile, "nBkgTimeErrLognormalComponents", 0), 0, 1);
	if (nTimeErrGaussComponents + nTimeErrLandauComponents + nTimeErrLognormalComponents <= 0)
	{
		std::cerr << "ERROR: at least one background timeErr component is required." << std::endl;
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
	RooFitResult *timeErrResult = dynamic_cast<RooFitResult *>(errFile->Get("timeErrResult"));
	if (errPdfOpt == kErrPdfAnalytic && !timeErrResult)
	{
		std::cerr << "ERROR: timeErrResult not found in " << errFileName << std::endl;
		return;
	}
	const double errUpper = std::max(errFitHigh, 1e-3);
	RooRealVar timeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", readSavedErrValue(*errFile, timeErrResult, "timeErrGaus1Mean", 0.02), errFitLow, errUpper);
	RooRealVar timeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", readSavedErrValue(*errFile, timeErrResult, "timeErrGaus1Sigma", 0.005), 1e-6, errUpper);
	RooRealVar timeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", readSavedErrValue(*errFile, timeErrResult, "timeErrGaus2Mean", timeErrGaus1Mean.getVal()), errFitLow, errUpper);
	RooRealVar timeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", readSavedErrValue(*errFile, timeErrResult, "timeErrGaus2Sigma", timeErrGaus1Sigma.getVal()), 1e-6, errUpper);
	RooRealVar timeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", readSavedErrValue(*errFile, timeErrResult, "timeErrTailMpv", 0.02), errFitLow, errUpper);
	RooRealVar timeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", readSavedErrValue(*errFile, timeErrResult, "timeErrTailWidth", 0.002), 1e-6, errUpper);
	RooRealVar timeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", readSavedErrValue(*errFile, timeErrResult, "timeErrTail2Mpv", timeErrTailMpv.getVal()), errFitLow, errUpper);
	RooRealVar timeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", readSavedErrValue(*errFile, timeErrResult, "timeErrTail2Width", timeErrTailWidth.getVal()), 1e-6, errUpper);
	RooRealVar timeErrLognM0("timeErrLognM0", "timeErrLognM0", readSavedErrValue(*errFile, timeErrResult, "timeErrLognM0", 0.02), errFitLow, errUpper);
	RooRealVar timeErrLognK("timeErrLognK", "timeErrLognK", readSavedErrValue(*errFile, timeErrResult, "timeErrLognK", 0.35), 0.01, 3.0);
	RooRealVar timeErrCore1FracRatio("timeErrCore1FracRatio", "timeErrCore1FracRatio", readSavedErrValue(*errFile, timeErrResult, "timeErrCore1FracRatio", 1.5), 1e-3, 1e3);
	RooFormulaVar timeErrCore1Frac("timeErrCore1Frac", "@0/(1.0+@0)", RooArgList(timeErrCore1FracRatio));
	RooRealVar timeErrTailFracRatio("timeErrTailFracRatio", "timeErrTailFracRatio", readSavedErrValue(*errFile, timeErrResult, "timeErrTailFracRatio", 0.08 / 0.92), 1e-4, 1e4);
	RooFormulaVar timeErrTailFrac("timeErrTailFrac", "@0/(1.0+@0)", RooArgList(timeErrTailFracRatio));
	RooRealVar timeErrTail2FracRatio("timeErrTail2FracRatio", "timeErrTail2FracRatio", readSavedErrValue(*errFile, timeErrResult, "timeErrTail2FracRatio", 0.05 / 0.95), 1e-4, 1e4);
	RooFormulaVar timeErrTail2Frac("timeErrTail2Frac", "@0/(1.0+@0)", RooArgList(timeErrTail2FracRatio));
	RooRealVar timeErrLognFracRatio("timeErrLognFracRatio", "timeErrLognFracRatio", readSavedErrValue(*errFile, timeErrResult, "timeErrLognFracRatio", 0.03 / 0.97), 1e-4, 1e4);
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
	RooAddPdf timeErrPdf("timeErrPdf", "timeErrPdf", timeErrPdfList, timeErrFracList, true);
	auto fixVar = [](RooRealVar &var) { var.setConstant(true); };
	fixVar(timeErrGaus1Mean);
	fixVar(timeErrGaus1Sigma);
	fixVar(timeErrGaus2Mean);
	fixVar(timeErrGaus2Sigma);
	fixVar(timeErrTailMpv);
	fixVar(timeErrTailWidth);
	fixVar(timeErrTail2Mpv);
	fixVar(timeErrTail2Width);
	fixVar(timeErrLognM0);
	fixVar(timeErrLognK);
	fixVar(timeErrCore1FracRatio);
	fixVar(timeErrTailFracRatio);
	fixVar(timeErrTail2FracRatio);
	fixVar(timeErrLognFracRatio);

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
	RooPlot *timeErrPlot = obs_timeErr.frame(Range(errFitLow, errFitHigh), Title(""));
	if (isWeight)
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errFitLow, errFitHigh), DataError(RooAbsData::SumW2), Name("data"));
	else
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errFitLow, errFitHigh), Name("data"));
	timeErrPdf.plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (useLandau1)
		timeErrPdf.plotOn(timeErrPlot, Components(timeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
	if (useLandau2)
		timeErrPdf.plotOn(timeErrPlot, Components(timeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
	if (useLognormal)
		timeErrPdf.plotOn(timeErrPlot, Components(timeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
	if (useGaus1)
		timeErrPdf.plotOn(timeErrPlot, Components(timeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
	if (useGaus2)
		timeErrPdf.plotOn(timeErrPlot, Components(timeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
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
		timeErrLeg.AddEntry(o, "Fit", "l");
	if (useLandau1)
		if (auto *o = findObj(timeErrPlot, "tail"))
		timeErrLeg.AddEntry(o, "Landau tail", "l");
	if (useLandau2)
		if (auto *o = findObj(timeErrPlot, "tail2"))
		timeErrLeg.AddEntry(o, "Landau tail 2", "l");
	if (useLognormal)
		if (auto *o = findObj(timeErrPlot, "logn"))
		timeErrLeg.AddEntry(o, "Log-normal tail", "l");
	if (useGaus1)
		if (auto *o = findObj(timeErrPlot, "gaus1"))
		timeErrLeg.AddEntry(o, "Gauss 1", "l");
	if (useGaus2)
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

	cTimeErr->cd();
	TPad *timeErrPad2 = new TPad("timeErrPad2", "timeErrPad2", 0.0, 0.0, 1.0, 0.25);
	timeErrPad2->SetTopMargin(0.00001);
	timeErrPad2->SetBottomMargin(0.4);
	timeErrPad2->Draw();
	timeErrPad2->cd();

	RooPlot *timeErrPullPlot = obs_timeErr.frame(Range(errFitLow, errFitHigh), Title(""));
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
	const int timeErrNPar = timeErrResult ? timeErrResult->floatParsFinal().getSize() : 0;
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

	TLine timeErrLine(errFitLow, 0.0, errFitHigh, 0.0);
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
	resolutionFile->Close();
	std::cout << "Loaded ctau prompt resolution from " << resolutionFileName << std::endl;

	// ------------------------------------------------------------------
	// build background ctau model
	// ------------------------------------------------------------------
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
	ctauFrac1.setConstant(true);
	ctauTime3Mean.setConstant(true);
	ctauFrac2.setConstant(true);
	ctauTime4Mean.setConstant(true);
	ctauFrac3.setConstant(true);

	std::unique_ptr<RooAddModel> promptResolutionModel;
	RooResolutionModel *time_resolution_ptr = &ctauTime1;
	RooAbsPdf *bkg_time_plb_ptr = &ctauTime1;
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
			resolutionPdfList, resolutionFracList
		);
		time_resolution_ptr = promptResolutionModel.get();
		bkg_time_plb_ptr = promptResolutionModel.get();
	}
	RooResolutionModel &time_resolution = *time_resolution_ptr;
	RooAbsPdf &bkg_time_plb = *bkg_time_plb_ptr;
	double maxPromptScaleSaved = std::max(ctauTime1ScaleVal, 0.0);
	if (nResolutionComponents >= 2 && std::isfinite(ctauTime2ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime2ScaleVal);
	if (nResolutionComponents >= 3 && std::isfinite(ctauTime3ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime3ScaleVal);
	if (nResolutionComponents >= 4 && std::isfinite(ctauTime4ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime4ScaleVal);
	const double resolutionDrivenLifetimeFloor =
		std::max(absoluteLifetimeFloor, lifetimeFloorCoeff * maxPromptScaleSaved * std::max(errFitHigh, absoluteLifetimeFloor));
	bkgLifetimeFloor = std::max(bkgLifetimeFloor, resolutionDrivenLifetimeFloor);
	bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 2e-2);
	bkgFlipLifetimeFloor = std::max(std::max(0.75 * bkgLifetimeFloor, 5e-3), resolutionDrivenLifetimeFloor);
	maxStableLifetime = std::max(maxStableLifetime, 20.0 * bkgLifetimeFloor);

	const double bkgLifetimeInit = std::max(0.05, 1.5 * bkgLifetimeFloor);
	const double bkgLifetimeCeil = std::max(0.50, maxStableLifetime);
	RooRealVar bkg_ss_lifetime(
		"bkg_ss_lifetime", "bkg_ss_lifetime",
		bkgLifetimeInit,
		bkgLifetimeFloor + 1e-2,
		bkgLifetimeCeil);

	const double bkgLifetime2Init = std::max(0.04, 1.5 * bkgLifetime2Floor);
	const double bkgLifetime2Ceil = std::max(0.50, maxStableLifetime);
	RooRealVar bkg_ss_lifetime2(
		"bkg_ss_lifetime2", "bkg_ss_lifetime2",
		bkgLifetime2Init,
		bkgLifetime2Floor + 2e-2,
		bkgLifetime2Ceil);
	const double bkgLifetime3Floor = std::max(2.0 * bkgLifetimeFloor, 3e-2);
	const double bkgLifetime3Init = std::max(0.10, 1.5 * bkgLifetime3Floor);
	const double bkgLifetime3Ceil = std::max(0.80, maxStableLifetime);
	RooRealVar bkg_ss_lifetime3(
		"bkg_ss_lifetime3", "bkg_ss_lifetime3",
		bkgLifetime3Init,
		bkgLifetime3Floor + 1e-2,
		bkgLifetime3Ceil);

	const double bkgSymLifetimeInit = std::max(0.04, 1.4 * bkgSymLifetimeFloor);
	const double bkgSymLifetimeCeil = std::max(0.30, 0.5 * maxStableLifetime);
	const double bkgSymLifetimeMinGap = std::max(0.30 * bkgSymLifetimeFloor, 2e-3);
	const double bkgSymLifetimeMaxGap = std::max(bkgSymLifetimeCeil - bkgSymLifetimeFloor, 2.0 * bkgSymLifetimeMinGap);
	RooRealVar bkg_ds_lifetime(
		"bkg_ds_lifetime", "bkg_ds_lifetime",
		bkgSymLifetimeInit,
		bkgSymLifetimeFloor + bkgSymLifetimeMinGap,
		bkgSymLifetimeCeil);
	const double bkgSymLifetime2Floor = std::max(2.0 * bkgSymLifetimeFloor, 2e-2);
	const double bkgSymLifetime2Init = std::max(0.07, 1.4 * bkgSymLifetime2Floor);
	const double bkgSymLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgSymLifetime2MinGap = std::max(0.25 * bkgSymLifetime2Floor, 3e-3);
	const double bkgSymLifetime2MaxGap = std::max(bkgSymLifetime2Ceil - bkgSymLifetime2Floor, 2.0 * bkgSymLifetime2MinGap);
	RooRealVar bkg_ds_lifetime2(
		"bkg_ds_lifetime2", "bkg_ds_lifetime2",
		bkgSymLifetime2Init,
		bkgSymLifetime2Floor + bkgSymLifetime2MinGap,
		bkgSymLifetime2Ceil);
	const double bkgSymLifetime3Floor = std::max(3.0 * bkgSymLifetimeFloor, 4e-2);
	const double bkgSymLifetime3Init = std::max(0.12, 1.3 * bkgSymLifetime3Floor);
	const double bkgSymLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetime3MinGap = std::max(0.20 * bkgSymLifetime3Floor, 4e-3);
	const double bkgSymLifetime3MaxGap = std::max(bkgSymLifetime3Ceil - bkgSymLifetime3Floor, 2.0 * bkgSymLifetime3MinGap);
	RooRealVar bkg_ds_lifetime3(
		"bkg_ds_lifetime3", "bkg_ds_lifetime3",
		bkgSymLifetime3Init,
		bkgSymLifetime3Floor + bkgSymLifetime3MinGap,
		bkgSymLifetime3Ceil);

	const double bkgFlipLifetimeInit = std::max(0.05, 1.2 * bkgFlipLifetimeFloor);
	const double bkgFlipLifetimeCeil = std::max(0.50, maxStableLifetime);
	RooRealVar bkg_flip_lifetime(
		"bkg_flip_lifetime", "bkg_flip_lifetime",
		bkgFlipLifetimeInit,
		bkgFlipLifetimeFloor + 1e-2,
		bkgFlipLifetimeCeil);
	const double bkgFlipLifetime2Floor = std::max(1.5 * bkgFlipLifetimeFloor, 1e-2);
	const double bkgFlipLifetime2Init = std::max(0.08, 1.3 * bkgFlipLifetime2Floor);
	const double bkgFlipLifetime2Ceil = std::max(0.60, maxStableLifetime);
	RooRealVar bkg_flip_lifetime2(
		"bkg_flip_lifetime2", "bkg_flip_lifetime2",
		bkgFlipLifetime2Init,
		bkgFlipLifetime2Floor + 2e-2,
		bkgFlipLifetime2Ceil);
	const double bkgFlipLifetime3Floor = std::max(2.5 * bkgFlipLifetimeFloor, 2e-2);
	const double bkgFlipLifetime3Init = std::max(0.12, 1.2 * bkgFlipLifetime3Floor);
	const double bkgFlipLifetime3Ceil = std::max(0.80, maxStableLifetime);
	RooRealVar bkg_flip_lifetime3(
		"bkg_flip_lifetime3", "bkg_flip_lifetime3",
		bkgFlipLifetime3Init,
		bkgFlipLifetime3Floor + 1e-2,
		bkgFlipLifetime3Ceil);

	RooDecay bkg_time_ss("bkg_time_ss", "bkg_time_ss", obs_time, bkg_ss_lifetime, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss2("bkg_time_ss2", "bkg_time_ss2", obs_time, bkg_ss_lifetime2, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss3("bkg_time_ss3", "bkg_time_ss3", obs_time, bkg_ss_lifetime3, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ds("bkg_time_ds", "bkg_time_ds", obs_time, bkg_ds_lifetime, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds2("bkg_time_ds2", "bkg_time_ds2", obs_time, bkg_ds_lifetime2, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds3("bkg_time_ds3", "bkg_time_ds3", obs_time, bkg_ds_lifetime3, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_flip("bkg_time_flip", "bkg_time_flip", obs_time, bkg_flip_lifetime, time_resolution, RooDecay::Flipped);
	RooDecay bkg_time_flip2("bkg_time_flip2", "bkg_time_flip2", obs_time, bkg_flip_lifetime2, time_resolution, RooDecay::Flipped);
	RooDecay bkg_time_flip3("bkg_time_flip3", "bkg_time_flip3", obs_time, bkg_flip_lifetime3, time_resolution, RooDecay::Flipped);

	RooRealVar Nplb("Nplb", "Nplb", 0.15 * data->numEntries(), 0.0, 7.0 * std::max(1, data->numEntries()));
	RooRealVar Nss1("Nss1", "Nss1", 0.30 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nss2("Nss2", "Nss2", 0.20 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nss3("Nss3", "Nss3", 0.10 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nds1("Nds1", "Nds1", 0.20 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nds2("Nds2", "Nds2", 0.10 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nds3("Nds3", "Nds3", 0.05 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nflip1("Nflip1", "Nflip1", 0.15 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nflip2("Nflip2", "Nflip2", 0.08 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	RooRealVar Nflip3("Nflip3", "Nflip3", 0.04 * data->numEntries(), 0.0, 1.0 * std::max(1, data->numEntries()));
	auto resetBkgSeeds = [&]() {
		if (yLow == 1.6f && ptLow == 12.0f && ptHigh == 14.0f)
		{
			ctauTime1Scale.setVal(0.86);
			if (nResolutionComponents >= 2)
				ctauTime2Delta.setVal(0.94);
			if (nResolutionComponents >= 3)
				ctauTime3Delta.setVal(3.57);
			Nplb.setVal(40.0);
			Nss1.setVal(100.0);
			Nds1.setVal(30.0);
			Nflip1.setVal(8.0);
			bkg_ss_lifetime.setVal(0.49);
			bkg_ds_lifetime.setVal(0.20);
			bkg_flip_lifetime.setVal(0.11);
		}
	};
	resetBkgSeeds();

	RooArgList bkgTimePdfList;
	RooArgList bkgTimeYieldList;
	bkgTimePdfList.add(bkg_time_plb);
	bkgTimeYieldList.add(Nplb);
	if (nSignalSSComponents >= 1)
	{
		bkgTimePdfList.add(bkg_time_ss);
		bkgTimeYieldList.add(Nss1);
	}
	if (nSignalSSComponents >= 2)
	{
		bkgTimePdfList.add(bkg_time_ss2);
		bkgTimeYieldList.add(Nss2);
	}
	if (nSignalSSComponents >= 3)
	{
		bkgTimePdfList.add(bkg_time_ss3);
		bkgTimeYieldList.add(Nss3);
	}
	if (nSignalDSComponents >= 1)
	{
		bkgTimePdfList.add(bkg_time_ds);
		bkgTimeYieldList.add(Nds1);
	}
	if (nSignalDSComponents >= 2)
	{
		bkgTimePdfList.add(bkg_time_ds2);
		bkgTimeYieldList.add(Nds2);
	}
	if (nSignalDSComponents >= 3)
	{
		bkgTimePdfList.add(bkg_time_ds3);
		bkgTimeYieldList.add(Nds3);
	}
	if (nSignalFlipComponents >= 1)
	{
		bkgTimePdfList.add(bkg_time_flip);
		bkgTimeYieldList.add(Nflip1);
	}
	if (nSignalFlipComponents >= 2)
	{
		bkgTimePdfList.add(bkg_time_flip2);
		bkgTimeYieldList.add(Nflip2);
	}
	if (nSignalFlipComponents >= 3)
	{
		bkgTimePdfList.add(bkg_time_flip3);
		bkgTimeYieldList.add(Nflip3);
	}
	RooAddPdf bkg_time("bkg_time", "bkg_time", bkgTimePdfList, bkgTimeYieldList);
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
	if (nResolutionComponents >= 2)
		addConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal,
			constraintSigmaFloor(ctauTime2DeltaVal, 0.10, 0.05));
	if (nResolutionComponents >= 3)
		addConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal,
			constraintSigmaFloor(ctauTime3DeltaVal, 0.10, 0.05));
	if (nResolutionComponents >= 4)
		addConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal,
			constraintSigmaFloor(ctauTime4DeltaVal, 0.10, 0.05));

	// ------------------------------------------------------------------
	// fit background ctau model
	// ------------------------------------------------------------------
	RooProdPdf time_pdf("time_pdf", "time_pdf", RooArgSet(timeErrPdf), Conditional(RooArgSet(bkg_time), RooArgSet(obs_time)));
	RooFitResult *time_result = nullptr;
	std::unique_ptr<RooFitResult> savedTimeResult;
	if (drawFromSavedFit)
	{
		savedTimeResult = clone_saved_fit_result(savedFitFile.get(), "timeResult");
		time_result = savedTimeResult.get();
		if (!time_result)
		{
			std::cerr << "ERROR: timeResult not found in saved background ctau file: " << fitResultFileName << std::endl;
			return;
		}
		apply_saved_fit_result(time_result, time_pdf, RooArgSet(obs_time, obs_timeErr));
	}
	else
	{
		if (timeConstraints.getSize() > 0)
			time_result = time_pdf.fitTo(
				*data, Extended(), Save(),
				// Strategy(2),
				SumW2Error(isWeight), PrefitDataFraction(errPrefitDataFraction),
				ExternalConstraints(timeConstraints), RecoverFromUndefinedRegions(1.0));
		else
			time_result = time_pdf.fitTo(
				*data, Extended(), Save(),
				// Strategy(2),
				SumW2Error(isWeight), PrefitDataFraction(errPrefitDataFraction), RecoverFromUndefinedRegions(1.0));
	}
	if (!drawFromSavedFit && time_result && (time_result->status() != 0 || time_result->covQual() < 2))
	{
		delete time_result;
		resetBkgSeeds();
		if (timeConstraints.getSize() > 0)
			time_result = time_pdf.fitTo(
				*data, Extended(), Save(), SumW2Error(isWeight), PrefitDataFraction(bkgRetryPrefitDataFraction),
				ExternalConstraints(timeConstraints),
				Strategy(2), Offset(true),
				RecoverFromUndefinedRegions(1.0));
		else
			time_result = time_pdf.fitTo(
				*data, Extended(), Save(), SumW2Error(isWeight), PrefitDataFraction(bkgRetryPrefitDataFraction),
				Strategy(2), Offset(true),
				RecoverFromUndefinedRegions(1.0));
		std::cout << "[OK] ctau_bkg plain refit done" << std::endl;
	}
	if (time_result)
		time_result->Print();

	RooProdPdf time_pdf_ss1(
		"time_pdf_ss1", "time_pdf_ss1",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ss), RooArgSet(obs_time)));
	RooProdPdf time_pdf_plb(
		"time_pdf_plb", "time_pdf_plb",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_plb), RooArgSet(obs_time)));
	RooProdPdf time_pdf_ss2(
		"time_pdf_ss2", "time_pdf_ss2",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ss2), RooArgSet(obs_time)));
	RooProdPdf time_pdf_ss3(
		"time_pdf_ss3", "time_pdf_ss3",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ss3), RooArgSet(obs_time)));
	RooProdPdf time_pdf_ds(
		"time_pdf_ds", "time_pdf_ds",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ds), RooArgSet(obs_time)));
	RooProdPdf time_pdf_ds2(
		"time_pdf_ds2", "time_pdf_ds2",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ds2), RooArgSet(obs_time)));
	RooProdPdf time_pdf_ds3(
		"time_pdf_ds3", "time_pdf_ds3",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_ds3), RooArgSet(obs_time)));
	RooProdPdf time_pdf_flip(
		"time_pdf_flip", "time_pdf_flip",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_flip), RooArgSet(obs_time)));
	RooProdPdf time_pdf_flip2(
		"time_pdf_flip2", "time_pdf_flip2",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_flip2), RooArgSet(obs_time)));
	RooProdPdf time_pdf_flip3(
		"time_pdf_flip3", "time_pdf_flip3",
		RooArgSet(timeErrPdf),
		Conditional(RooArgSet(bkg_time_flip3), RooArgSet(obs_time)));

	// ------------------------------------------------------------------
	// draw background ctau fit
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
	time_pdf_plb.plotOn(timefitplot, LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2),
		Normalization(Nplb.getVal(), RooAbsReal::NumEvent), Name("plb_component"));
	if (nSignalSSComponents >= 1)
		time_pdf_ss1.plotOn(timefitplot, LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2),
			Normalization(Nss1.getVal(), RooAbsReal::NumEvent), Name("ss1_component"));
	if (nSignalSSComponents >= 2)
		time_pdf_ss2.plotOn(timefitplot, LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2),
			Normalization(Nss2.getVal(), RooAbsReal::NumEvent), Name("ss2_component"));
	if (nSignalSSComponents >= 3)
		time_pdf_ss3.plotOn(timefitplot, LineColor(kAzure + 2), LineStyle(kDashed), LineWidth(2),
			Normalization(Nss3.getVal(), RooAbsReal::NumEvent), Name("ss3_component"));
	if (nSignalDSComponents >= 1)
		time_pdf_ds.plotOn(timefitplot, LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2),
			Normalization(Nds1.getVal(), RooAbsReal::NumEvent), Name("ds1_component"));
	if (nSignalDSComponents >= 2)
		time_pdf_ds2.plotOn(timefitplot, LineColor(kViolet + 1), LineStyle(kDashed), LineWidth(2),
			Normalization(Nds2.getVal(), RooAbsReal::NumEvent), Name("ds2_component"));
	if (nSignalDSComponents >= 3)
		time_pdf_ds3.plotOn(timefitplot, LineColor(kPink + 7), LineStyle(kDashed), LineWidth(2),
			Normalization(Nds3.getVal(), RooAbsReal::NumEvent), Name("ds3_component"));
	if (nSignalFlipComponents >= 1)
		time_pdf_flip.plotOn(timefitplot, LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2),
			Normalization(Nflip1.getVal(), RooAbsReal::NumEvent), Name("flip1_component"));
	if (nSignalFlipComponents >= 2)
		time_pdf_flip2.plotOn(timefitplot, LineColor(kOrange - 3), LineStyle(kDashed), LineWidth(2),
			Normalization(Nflip2.getVal(), RooAbsReal::NumEvent), Name("flip2_component"));
	if (nSignalFlipComponents >= 3)
		time_pdf_flip3.plotOn(timefitplot, LineColor(kYellow + 2), LineStyle(kDashed), LineWidth(2),
			Normalization(Nflip3.getVal(), RooAbsReal::NumEvent), Name("flip3_component"));
	apply_logy_auto_range(timefitplot, "data");
	timefitplot->GetYaxis()->SetTitle("Events");
	timefitplot->GetYaxis()->SetTitleOffset(1.6);
	timefitplot->GetXaxis()->SetTitle("");
	timefitplot->Draw("e");

	TLegend leg(0.50, 0.66, 0.74, 0.89);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.03);
	if (auto *o = findObj(timefitplot, "data"))
		leg.AddEntry(o, "Data", "lep");
	if (auto *o = findObj(timefitplot, "model"))
		leg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(timefitplot, "plb_component"))
		leg.AddEntry(o, "PLB", "l");
	if (nSignalSSComponents >= 1)
		if (auto *o = findObj(timefitplot, "ss1_component"))
			leg.AddEntry(o, "SS1", "l");
	if (nSignalSSComponents >= 2)
		if (auto *o = findObj(timefitplot, "ss2_component"))
			leg.AddEntry(o, "SS2", "l");
	if (nSignalSSComponents >= 3)
		if (auto *o = findObj(timefitplot, "ss3_component"))
			leg.AddEntry(o, "SS3", "l");
	if (nSignalDSComponents >= 1)
		if (auto *o = findObj(timefitplot, "ds1_component"))
			leg.AddEntry(o, "DS1", "l");
	if (nSignalDSComponents >= 2)
		if (auto *o = findObj(timefitplot, "ds2_component"))
			leg.AddEntry(o, "DS2", "l");
	if (nSignalDSComponents >= 3)
		if (auto *o = findObj(timefitplot, "ds3_component"))
			leg.AddEntry(o, "DS3", "l");
	if (nSignalFlipComponents >= 1)
		if (auto *o = findObj(timefitplot, "flip1_component"))
			leg.AddEntry(o, "Flip1", "l");
	if (nSignalFlipComponents >= 2)
		if (auto *o = findObj(timefitplot, "flip2_component"))
			leg.AddEntry(o, "Flip2", "l");
	if (nSignalFlipComponents >= 3)
		if (auto *o = findObj(timefitplot, "flip3_component"))
			leg.AddEntry(o, "Flip3", "l");
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
		double xtext = 0.73, y0 = 0.87, dy = -0.04;
		int k = 0;
		auto printVar = [&](const char *title, const RooAbsReal &var)
		{
			const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
			if (rrv && rrv->isConstant())
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
			else if (rrv)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
			else if (time_result)
			{
				const double err = var.getPropagatedError(*time_result);
				if (err > 0.0 && std::isfinite(err))
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
				else
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
			}
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
		};
		printVar("N_{PLB}", Nplb);
		if (nSignalSSComponents >= 1)
			printVar("N_{SS1}", Nss1);
		if (nSignalSSComponents >= 2)
			printVar("N_{SS2}", Nss2);
		if (nSignalSSComponents >= 3)
			printVar("N_{SS3}", Nss3);
		if (nSignalDSComponents >= 1)
			printVar("N_{DS1}", Nds1);
		if (nSignalDSComponents >= 2)
			printVar("N_{DS2}", Nds2);
		if (nSignalDSComponents >= 3)
			printVar("N_{DS3}", Nds3);
		if (nSignalFlipComponents >= 1)
			printVar("N_{Flip1}", Nflip1);
		if (nSignalFlipComponents >= 2)
			printVar("N_{Flip2}", Nflip2);
		if (nSignalFlipComponents >= 3)
			printVar("N_{Flip3}", Nflip3);
		printVar("#sigma_{1}^{res}", ctauTime1Scale);
		if (nResolutionComponents >= 2 && !ctauTime2Delta.isConstant())
			printVar("#Delta s_{21}^{res}", ctauTime2Delta);
		if (nResolutionComponents >= 3 && !ctauTime3Delta.isConstant())
			printVar("#Delta s_{32}^{res}", ctauTime3Delta);
		if (nResolutionComponents >= 4 && !ctauTime4Delta.isConstant())
			printVar("#Delta s_{43}^{res}", ctauTime4Delta);
		if (nSignalSSComponents >= 1)
			printVar("#tau_{SS1}", bkg_ss_lifetime);
		if (nSignalSSComponents >= 2)
			printVar("#tau_{SS2}", bkg_ss_lifetime2);
		if (nSignalSSComponents >= 3)
			printVar("#tau_{SS3}", bkg_ss_lifetime3);
		if (nSignalDSComponents >= 1)
			printVar("#tau_{DS1}", bkg_ds_lifetime);
		if (nSignalDSComponents >= 2)
			printVar("#tau_{DS2}", bkg_ds_lifetime2);
		if (nSignalDSComponents >= 3)
			printVar("#tau_{DS3}", bkg_ds_lifetime3);
		if (nSignalFlipComponents >= 1)
			printVar("#tau_{Flip1}", bkg_flip_lifetime);
		if (nSignalFlipComponents >= 2)
			printVar("#tau_{Flip2}", bkg_flip_lifetime2);
		if (nSignalFlipComponents >= 3)
			printVar("#tau_{Flip3}", bkg_flip_lifetime3);
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

	if (!drawFromSavedFit)
	{
		TFile *fitResultFile = TFile::Open(fitResultFileName, "RECREATE");
		if (!fitResultFile || fitResultFile->IsZombie())
		{
			std::cerr << "ERROR: cannot create output fit result file: " << fitResultFileName << std::endl;
			return;
		}
		if (timeErrResult)
			timeErrResult->Write("timeErrResult");
		if (time_result)
			time_result->Write("timeResult");
		TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
		TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
		TParameter<int>("nSignalFlipComponents", nSignalFlipComponents).Write();
		TParameter<int>("nSignalDSComponents", nSignalDSComponents).Write();
		TParameter<double>("ctRangeMin", ctRange.first).Write();
		TParameter<double>("ctRangeMax", ctRange.second).Write();
		TParameter<double>("errRangeMin", errRange.first).Write();
		TParameter<double>("errRangeMax", errRange.second).Write();
		TParameter<double>("absoluteLifetimeFloor", absoluteLifetimeFloor).Write();
		TParameter<double>("lifetimeFloorCoeff", lifetimeFloorCoeff).Write();
		TParameter<double>("maxPromptScaleSaved", maxPromptScaleSaved).Write();
		TParameter<double>("resolutionDrivenLifetimeFloor", resolutionDrivenLifetimeFloor).Write();
		TParameter<double>("bkgLifetimeFloorApplied", bkgLifetimeFloor).Write();
		TParameter<double>("bkgLifetime2FloorApplied", bkgLifetime2Floor).Write();
		TParameter<double>("bkgSymLifetimeFloorApplied", bkgSymLifetimeFloor).Write();
		TParameter<double>("bkgFlipLifetimeFloorApplied", bkgFlipLifetimeFloor).Write();
		fitResultFile->Write();
		fitResultFile->Close();
		std::cout << "Saved ctau background fit results to " << fitResultFileName << std::endl;
	}
	else
		std::cout << "[PlotOnly] Loaded saved background ctau fit and left ROOT file unchanged: " << fitResultFileName << std::endl;

	cout << "------------------ FIT RESULT SUMMARY --------------------" << endl;
	cout << "errPdfOpt in BKG fit: " << errPdfOpt
			 << (errPdfOpt == kErrPdfHist ? " (RooHistPdf)" : " (analytic)") << endl;
	cout << "histPdfInterpolationOrder in BKG fit: " << histPdfInterpolationOrder << endl;
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	if (timeErrResult)
	{
		cout << "------------------ FIT RESULT FOR TIME ERR ---------------" << endl;
		timeErrResult->Print("v");
	}
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
		if (nSignalSSComponents >= 1)
			printFloorVsValue2Sigma("bkg_ss_lifetime", bkg_ss_lifetime, bkgLifetimeFloor);
		if (nSignalSSComponents >= 2)
			printFloorVsValue2Sigma("bkg_ss_lifetime2", bkg_ss_lifetime2, bkgLifetime2Floor);
		if (nSignalSSComponents >= 3)
			printFloorVsValue2Sigma("bkg_ss_lifetime3", bkg_ss_lifetime3, bkgLifetime3Floor);
		if (nSignalDSComponents >= 1)
			printFloorVsValue2Sigma("bkg_ds_lifetime", bkg_ds_lifetime, bkgSymLifetimeFloor);
		if (nSignalDSComponents >= 2)
			printFloorVsValue2Sigma("bkg_ds_lifetime2", bkg_ds_lifetime2, bkgSymLifetime2Floor);
		if (nSignalDSComponents >= 3)
			printFloorVsValue2Sigma("bkg_ds_lifetime3", bkg_ds_lifetime3, bkgSymLifetime3Floor);
		if (nSignalFlipComponents >= 1)
			printFloorVsValue2Sigma("bkg_flip_lifetime", bkg_flip_lifetime, bkgFlipLifetimeFloor);
		if (nSignalFlipComponents >= 2)
			printFloorVsValue2Sigma("bkg_flip_lifetime2", bkg_flip_lifetime2, bkgFlipLifetime2Floor);
		if (nSignalFlipComponents >= 3)
			printFloorVsValue2Sigma("bkg_flip_lifetime3", bkg_flip_lifetime3, bkgFlipLifetime3Floor);
	}
	cout << "Prompt resolution components used in BKG fit: " << nResolutionComponents << endl;
	cout << "Background components used in BKG fit: SS=" << nSignalSSComponents
			 << " Flip=" << nSignalFlipComponents
			 << " DS=" << nSignalDSComponents << endl;
	cout << "TimeErr components used in BKG fit: G=" << nTimeErrGaussComponents
			 << " L=" << nTimeErrLandauComponents
			 << " LN=" << nTimeErrLognormalComponents << endl;
	const TString figTimeErr = figName("timeerr_model");
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_bkg err fit : " << figTimeErr << std::endl;
	std::cout << "[FIG] ctau_bkg lifetime fit : " << figLifetime << std::endl;
}
