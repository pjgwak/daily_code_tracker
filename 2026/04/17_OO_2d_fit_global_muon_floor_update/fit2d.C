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
#include "RooCrystalBall.h"
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
#include "TString.h"
#include "TParameter.h"
#include "TMath.h"
#include "RooHist.h"
#include "RooLognormal.h"
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

void fit2d(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4,
	bool fixResolutionParams = false, bool drawFromSavedFit = false, bool publish = false){
	ScopedMacroTimer timer("fit2d", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;

	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	bool isSilence = true;
	enum ErrPdfChoice {
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
		kErrPdfHistPrMc = 2,
	};
	const int overrideBkgErrPdfOpt = -1; // -1: follow err2 setting
	const int overrideSigErrPdfOpt = -1; // ditto
	const int histPdfInterpolationOrder = 1;
	auto isHistErrPdfChoice = [](int choice) {
		return choice == kErrPdfHist || choice == kErrPdfHistPrMc;
	};
	auto errPdfChoiceLabel = [&](int choice) -> const char * {
		switch (choice)
		{
		case kErrPdfHist:
			return "RooHistPdf";
		case kErrPdfHistPrMc:
			return "RooHistPdf(PR_MC_ROOT)";
		default:
			return "analytic";
		}
	};
	if (isSilence)
	{
		RooMsgService::instance().getStream(0).removeTopic(RooFit::Tracing);
		RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
		RooMsgService::instance().getStream(1).removeTopic(RooFit::Tracing);
		RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
	}

	// First we open the actual RooDataSet used in the prompt ctau analysis.
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOn_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
	const char *PR_MC_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_PR_globalOn_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
	const char *DSET_NAME = "dataset";
	TFile *inputFile = TFile::Open(DATA_ROOT);
	if (!inputFile || inputFile->IsZombie()) {
		std::cerr << "ERROR: cannot open input data file: " << DATA_ROOT << std::endl;
		return;
	}

	RooDataSet *inputData = dynamic_cast<RooDataSet *>(inputFile->Get(DSET_NAME));
	if (!inputData) {
		std::cerr << "ERROR: RooDataSet '" << DSET_NAME << "' not found in " << DATA_ROOT << std::endl;
		return;
	}

	TFile *prMcFile = TFile::Open(PR_MC_ROOT);
	if (!prMcFile || prMcFile->IsZombie()) {
		std::cerr << "ERROR: cannot open prompt-MC file: " << PR_MC_ROOT << std::endl;
		return;
	}
	RooDataSet *prMcInput = dynamic_cast<RooDataSet *>(prMcFile->Get(DSET_NAME));
	if (!prMcInput) {
		std::cerr << "ERROR: RooDataSet '" << DSET_NAME << "' not found in " << PR_MC_ROOT << std::endl;
		return;
	}

	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)",
		ptLow, ptHigh, yLow, yHigh);
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	auto *massTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
	auto *timeTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
	auto *timeErrTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DErr"));
	if (!massTmp || !timeTmp || !timeErrTmp) {
		std::cerr << "ERROR: required variables mass/ctau3D/ctau3DErr are missing in dataset." << std::endl;
		return;
	}

	auto quantileRange = [&](RooDataSet &ds, RooRealVar &var, double qLo, double qHi, bool positiveOnly) {
		std::vector<double> vals;
		vals.reserve(ds.numEntries());
		for (int i = 0; i < ds.numEntries(); ++i) {
			ds.get(i);
			const double v = var.getVal();
			if (!std::isfinite(v)) continue;
			if (positiveOnly && v <= 0.0) continue;
			vals.push_back(v);
		}
		if (vals.empty()) return std::make_pair(var.getMin(), var.getMax());
		std::sort(vals.begin(), vals.end());
		auto qAt = [&](double q) {
			long long idx = llround(q * (vals.size() - 1));
			if (idx < 0) idx = 0;
			if (idx >= static_cast<long long>(vals.size())) idx = static_cast<long long>(vals.size()) - 1;
			return vals[static_cast<size_t>(idx)];
		};
		return std::make_pair(qAt(qLo), qAt(qHi));
	};

	auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.001, 0.999, false);
	auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.995, true);
	if (errRange.first < 1e-6) errRange.first = 1e-6;

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
			ctRange.first = -6;
			ctRange.second = 8;
		}
		if (ptLow == 2.0f && ptHigh == 3.0f)
		{
			ctRange.first = -3;
			ctRange.second = 3;
		}
	}
	const double errMin = 0.0;
	auto formatTag = [](double value) { return TString::Format("%.2f", value); };
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString baseFigTag = yTag + "_" + ptTag;
	const TString figTag = baseFigTag;
	const TString figDir = TString::Format("figs/%s/fit2d", yTag.Data());
	const TString fitRootDir = TString::Format("roots/%s/fit2d", yTag.Data());
	const TString fitRootName = TString::Format("%s/fit2d_result_%s.root", fitRootDir.Data(), figTag.Data());
	const TString errFileName = TString::Format("roots/%s/err2/err2_model_%s.root", yTag.Data(), baseFigTag.Data());
	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit && !load_saved_fit_file(savedFitFile, fitRootName, "2D"))
		return;

	auto readErrFileDouble = [](const TString &path, const char *name, double fallback) {
		std::unique_ptr<TFile> f(TFile::Open(path));
		if (!f || f->IsZombie()) return fallback;
		auto *param = dynamic_cast<TParameter<double> *>(f->Get(name));
		return param ? param->GetVal() : fallback;
	};
	const double errFitLow = readErrFileDouble(errFileName, "errLow", errMin);
	const double errFitHigh = readErrFileDouble(errFileName, "errHigh", errRange.second);

	TString cutAll = Form(
		"(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
		"(ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)",
		ptLow, ptHigh, yLow, yHigh, ctRange.first, ctRange.second, errFitLow, errFitHigh
	);
	auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutAll)));
	if (!dataSel || dataSel->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after final selection: " << cutAll << std::endl;
		return;
	}

	auto *massVar = static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	auto *timeVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	auto *timeErrVar = static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	if (!massVar || !timeVar || !timeErrVar) {
		std::cerr << "ERROR: required observables mass/ctau3D/ctau3DErr are missing in selected dataset." << std::endl;
		return;
	}
	RooRealVar obs_mass("mass", "mass", 2.6, 3.5);
	RooRealVar obs_time("ctau3D", "ctau3D", ctRange.first, ctRange.second);
	RooRealVar obs_timeErr("ctau3DErr", "ctau3DErr", errFitLow, errFitHigh);
	obs_mass.setBins(massVar->getBins());
	obs_time.setBins(timeVar->getBins());
	obs_timeErr.setBins(timeErrVar->getBins());
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_mass.setRange(2.6, 3.5);
	obs_mass.setMin(2.6);
	obs_mass.setMax(3.5);
	obs_time.SetTitle("ctau3D");
	obs_time.setUnit("mm");
	obs_timeErr.SetTitle("event-by-event ctau error");
	obs_timeErr.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	obs_timeErr.setRange(errFitLow, errFitHigh);
	obs_timeErr.setMin(errFitLow);
	const int timePlotBins = std::max(2, obs_time.getBins() * 2);
	const bool isLowPtForwardBin = (yLow == 1.6f && ptHigh <= 3.0f);

	// Keep the decay constants away from the numerically unstable tau -> 0 limit.
	const double absoluteLifetimeFloor = 1e-2;
	double lifetimeFloorCoeff = 0.25;
	if (yLow == 1.6f && ptHigh <= 4.0f)
		lifetimeFloorCoeff = 0.05;
	double signalLifetimeFloor = std::max(1.0 * errRange.first, absoluteLifetimeFloor);
	double bkgLifetimeFloor = std::max(1.0 * errRange.first, absoluteLifetimeFloor);
	double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	RooDataSet *data = dataSel.get();

	const double jpsiMassLow = 2.9;
	const double jpsiMassHigh = 3.3;
	const double sidebandLeftMax = 2.8;
	const double sidebandRightMin = 3.3;
	TString cutSignalRegion = Form("(mass >= %g && mass <= %g)", jpsiMassLow, jpsiMassHigh);
	TString cutSideband = Form("(mass >= 2.6 && mass < %g) || (mass > %g && mass <= 3.5)", sidebandLeftMax, sidebandRightMin);
	auto dataSB = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(cutSideband)));
	if (!dataSB || dataSB->numEntries() <= 0) dataSB.reset(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_mass, obs_time, obs_timeErr))));

	const int timeErrPlotBinsDefault = std::max(2, obs_timeErr.getBins() * 2);
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(fitRootDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	const TString massFileName = TString::Format("roots/%s/mass/mass_model_%s.root", yTag.Data(), baseFigTag.Data());
	const TString prFileName = TString::Format("roots/%s/ctau_pr/ctau_resolution_%s.root", yTag.Data(), baseFigTag.Data());
	const TString npFileName = TString::Format("roots/%s/ctau_np/ctau_np_model_%s.root", yTag.Data(), baseFigTag.Data());
	const TString bkgFileName = TString::Format("roots/%s/ctau_bkg/ctau_bkg_fitresult_%s.root", yTag.Data(), baseFigTag.Data());

	auto openFile = [](const TString &path) -> TFile * {
		TFile *f = TFile::Open(path);
		if (!f || f->IsZombie()) {
			std::cerr << "ERROR: cannot open " << path << std::endl;
			return nullptr;
		}
		return f;
	};
	auto readIntParam = [](TFile &f, const char *name, int fallback = -999) {
		auto *param = dynamic_cast<TParameter<int> *>(f.Get(name));
		return param ? param->GetVal() : fallback;
	};
	auto readDoubleParam = [](TFile &f, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN()) {
		auto *param = dynamic_cast<TParameter<double> *>(f.Get(name));
		return param ? param->GetVal() : fallback;
	};
	auto readResultValue = [](RooFitResult &fr, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN()) {
		auto *var = dynamic_cast<RooRealVar *>(fr.floatParsFinal().find(name));
		if (!var) var = dynamic_cast<RooRealVar *>(fr.constPars().find(name));
		return var ? var->getVal() : fallback;
	};
	auto readResultError = [](RooFitResult &fr, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN()) {
		auto *var = dynamic_cast<RooRealVar *>(fr.floatParsFinal().find(name));
		if (!var) var = dynamic_cast<RooRealVar *>(fr.constPars().find(name));
		if (!var)
			return fallback;
		const double err = var->getError();
		return (std::isfinite(err) && err > 0.0) ? err : fallback;
	};
	auto readBkgResultValue = [&](RooFitResult &fr, const char *primaryName, const char *legacyName,
		double fallback = std::numeric_limits<double>::quiet_NaN()) {
		const double primary = readResultValue(fr, primaryName, std::numeric_limits<double>::quiet_NaN());
		if (std::isfinite(primary))
			return primary;
		if (legacyName && legacyName[0] != '\0')
		{
			const double legacy = readResultValue(fr, legacyName, std::numeric_limits<double>::quiet_NaN());
			if (std::isfinite(legacy))
				return legacy;
		}
		return fallback;
	};

	std::unique_ptr<TFile> massFile(openFile(massFileName));
	std::unique_ptr<TFile> prFile(openFile(prFileName));
	std::unique_ptr<TFile> npFile(openFile(npFileName));
	std::unique_ptr<TFile> bkgFile(openFile(bkgFileName));
	std::unique_ptr<TFile> errFile(openFile(errFileName));
	if (!massFile || !prFile || !npFile || !bkgFile || !errFile) return;
	const int timeErrPlotBins = std::max(2, readIntParam(*errFile, "timeErrPlotBins", timeErrPlotBinsDefault));
	const double savedScaleBkg = readDoubleParam(*errFile, "scaleBkg", std::numeric_limits<double>::quiet_NaN());
	const int savedBkgErrPdfOpt = readIntParam(*errFile, "bkgErrPdfOpt", kErrPdfAnalytic);
	const int savedSigErrPdfOpt = readIntParam(*errFile, "sigErrPdfOpt", kErrPdfAnalytic);
	const int bkgErrPdfOpt = overrideBkgErrPdfOpt >= 0 ? overrideBkgErrPdfOpt : savedBkgErrPdfOpt;
	const int sigErrPdfOpt = overrideSigErrPdfOpt >= 0 ? overrideSigErrPdfOpt : savedSigErrPdfOpt;

	RooFitResult *bkgErrResult = dynamic_cast<RooFitResult *>(errFile->Get("timeErrResult"));
	RooFitResult *signalErrResult = dynamic_cast<RooFitResult *>(errFile->Get("sigTimeErrResult"));
	if (!bkgErrResult || (sigErrPdfOpt == kErrPdfAnalytic && !signalErrResult)) {
		std::cerr << "ERROR: timeErrResult/sigTimeErrResult not found in " << errFileName << std::endl;
		return;
	}
	if (!dataSB || dataSB->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after sideband selection: " << cutSideband << std::endl;
		return;
	}
	auto readErrVar = [&](RooFitResult &fr, const char *name, double fallback = 0.0) {
		return readResultValue(fr, name, fallback);
	};
	obs_timeErr.setRange(errFitLow, errFitHigh);
	obs_timeErr.setMin(errFitLow);
	auto readErrStoredValue = [&](const char *name, RooFitResult *fr, double fallback = 0.0) {
		if (auto *param = dynamic_cast<TParameter<double> *>(errFile->Get(name)))
			return param->GetVal();
		if (auto *var = dynamic_cast<RooRealVar *>(errFile->Get(name)))
			return var->getVal();
		if (auto *abs = dynamic_cast<RooAbsReal *>(errFile->Get(name)))
			return abs->getVal();
		return fr ? readErrVar(*fr, name, fallback) : fallback;
	};
	const int nErrBkgGauss = std::clamp(readIntParam(*errFile, "nBkgTimeErrGaussComponents", 1), 0, 2);
	const int nErrBkgLandau = std::clamp(readIntParam(*errFile, "nBkgTimeErrLandauComponents", 1), 0, 2);
	const int nErrBkgLogn = std::clamp(readIntParam(*errFile, "nBkgTimeErrLognormalComponents", 0), 0, 1);
	const int nErrSigGauss = std::clamp(readIntParam(*errFile, "nSigTimeErrGaussComponents", 1), 0, 2);
	const int nErrSigLandau = std::clamp(readIntParam(*errFile, "nSigTimeErrLandauComponents", 1), 0, 2);
	const int nErrSigLogn = std::clamp(readIntParam(*errFile, "nSigTimeErrLognormalComponents", 0), 0, 1);
	const bool useErrBkgGaus1 = (nErrBkgGauss >= 1);
	const bool useErrBkgGaus2 = (nErrBkgGauss >= 2);
	const bool useErrBkgLandau1 = (nErrBkgLandau >= 1);
	const bool useErrBkgLandau2 = (nErrBkgLandau >= 2);
	const bool useErrBkgLogn = (nErrBkgLogn >= 1);
	const int nErrBkgComponents = static_cast<int>(useErrBkgGaus1) + static_cast<int>(useErrBkgGaus2) +
		static_cast<int>(useErrBkgLandau1) + static_cast<int>(useErrBkgLandau2) + static_cast<int>(useErrBkgLogn);
	const bool useErrSigGaus1 = (nErrSigGauss >= 1);
	const bool useErrSigGaus2 = (nErrSigGauss >= 2);
	const bool useErrSigLandau1 = (nErrSigLandau >= 1);
	const bool useErrSigLandau2 = (nErrSigLandau >= 2);
	const bool useErrSigLogn = (nErrSigLogn >= 1);
	const int nErrSigComponents = static_cast<int>(useErrSigGaus1) + static_cast<int>(useErrSigGaus2) +
		static_cast<int>(useErrSigLandau1) + static_cast<int>(useErrSigLandau2) + static_cast<int>(useErrSigLogn);

	auto bkgErrData = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSB->reduce(RooArgSet(obs_timeErr))));
	auto dataSR = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(cutSignalRegion)));
	auto prMcSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(prMcInput->reduce(cutAll)));
	auto prMcErrData = std::unique_ptr<RooDataSet>();
	if (!dataSR || dataSR->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after signal-region selection: " << cutSignalRegion << std::endl;
		return;
	}
	if (!prMcSel || prMcSel->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after prompt-MC selection: " << cutAll << std::endl;
		return;
	}
	auto *prMcTimeErr = dynamic_cast<RooRealVar *>(prMcSel->get()->find("ctau3DErr"));
	if (!prMcTimeErr) {
		std::cerr << "ERROR: ctau3DErr is missing in selected prompt-MC dataset" << std::endl;
		return;
	}
	prMcTimeErr->setRange(errFitLow, errFitHigh);
	prMcTimeErr->setMin(errFitLow);
	prMcErrData.reset(static_cast<RooDataSet *>(prMcSel->reduce(RooArgSet(*prMcTimeErr))));
	if (!prMcErrData || prMcErrData->numEntries() <= 0) {
		std::cerr << "ERROR: no prompt-MC ctau3DErr entries available after selection" << std::endl;
		return;
	}

	const double errUpper = std::max(errFitHigh, 1e-3);
	RooRealVar bkgTimeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", readErrStoredValue("timeErrGaus1Mean", bkgErrResult, 0.02), errFitLow, errUpper);
	RooRealVar bkgTimeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", readErrStoredValue("timeErrGaus1Sigma", bkgErrResult, 0.005), 1e-6, errUpper);
	RooRealVar bkgTimeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", readErrStoredValue("timeErrGaus2Mean", bkgErrResult, bkgTimeErrGaus1Mean.getVal()), errFitLow, errUpper);
	RooRealVar bkgTimeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", readErrStoredValue("timeErrGaus2Sigma", bkgErrResult, bkgTimeErrGaus1Sigma.getVal()), 1e-6, errUpper);
	RooRealVar bkgTimeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", readErrStoredValue("timeErrTailMpv", bkgErrResult, 0.02), errFitLow, errUpper);
	RooRealVar bkgTimeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", readErrStoredValue("timeErrTailWidth", bkgErrResult, 0.002), 1e-6, errUpper);
	RooRealVar bkgTimeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", readErrStoredValue("timeErrTail2Mpv", bkgErrResult, bkgTimeErrTailMpv.getVal()), errFitLow, errUpper);
	RooRealVar bkgTimeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", readErrStoredValue("timeErrTail2Width", bkgErrResult, bkgTimeErrTailWidth.getVal()), 1e-6, errUpper);
	RooRealVar bkgTimeErrLognM0("timeErrLognM0", "timeErrLognM0", readErrStoredValue("timeErrLognM0", bkgErrResult, 0.02), errFitLow, errUpper);
	RooRealVar bkgTimeErrLognK("timeErrLognK", "timeErrLognK", readErrStoredValue("timeErrLognK", bkgErrResult, 0.35), 0.01, 3.0);
	RooRealVar bkgTimeErrCore1FracRatio("timeErrCore1FracRatio", "timeErrCore1FracRatio", readErrStoredValue("timeErrCore1FracRatio", bkgErrResult, 1.5), 1e-3, 1e3);
	RooFormulaVar bkgTimeErrCore1Frac("timeErrCore1Frac", "@0/(1.0+@0)", RooArgList(bkgTimeErrCore1FracRatio));
	RooRealVar bkgTimeErrTailFracRatio("timeErrTailFracRatio", "timeErrTailFracRatio", readErrStoredValue("timeErrTailFracRatio", bkgErrResult, 0.08 / 0.92), 1e-4, 1e4);
	RooFormulaVar bkgTimeErrTailFrac("timeErrTailFrac", "@0/(1.0+@0)", RooArgList(bkgTimeErrTailFracRatio));
	RooRealVar bkgTimeErrTail2FracRatio("timeErrTail2FracRatio", "timeErrTail2FracRatio", readErrStoredValue("timeErrTail2FracRatio", bkgErrResult, 0.05 / 0.95), 1e-4, 1e4);
	RooFormulaVar bkgTimeErrTail2Frac("timeErrTail2Frac", "@0/(1.0+@0)", RooArgList(bkgTimeErrTail2FracRatio));
	RooRealVar bkgTimeErrLognFracRatio("timeErrLognFracRatio", "timeErrLognFracRatio", readErrStoredValue("timeErrLognFracRatio", bkgErrResult, 0.03 / 0.97), 1e-4, 1e4);
	RooFormulaVar bkgTimeErrLognFrac("timeErrLognFrac", "@0/(1.0+@0)", RooArgList(bkgTimeErrLognFracRatio));
	RooGaussian bkgTimeErrGaus1("gausBkg", "gausBkg", obs_timeErr, bkgTimeErrGaus1Mean, bkgTimeErrGaus1Sigma);
	RooGaussian bkgTimeErrGaus2("gausBkg2", "gausBkg2", obs_timeErr, bkgTimeErrGaus2Mean, bkgTimeErrGaus2Sigma);
	RooLandau bkgTimeErrTail("landauBkg", "landauBkg", obs_timeErr, bkgTimeErrTailMpv, bkgTimeErrTailWidth);
	RooLandau bkgTimeErrTail2("landauBkg2", "landauBkg2", obs_timeErr, bkgTimeErrTail2Mpv, bkgTimeErrTail2Width);
	RooLognormal bkgTimeErrLogn("lognBkg", "lognBkg", obs_timeErr, bkgTimeErrLognM0, bkgTimeErrLognK);
	RooArgList bkgErrPdfList;
	RooArgList bkgErrFracList;
	if (useErrBkgLandau1) {
		bkgErrPdfList.add(bkgTimeErrTail);
		if (nErrBkgComponents > 1) bkgErrFracList.add(bkgTimeErrTailFrac);
	}
	if (useErrBkgLandau2) {
		bkgErrPdfList.add(bkgTimeErrTail2);
		if (static_cast<int>(bkgErrPdfList.getSize()) < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrTail2Frac);
	}
	if (useErrBkgLogn) {
		bkgErrPdfList.add(bkgTimeErrLogn);
		if (static_cast<int>(bkgErrPdfList.getSize()) < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrLognFrac);
	}
	if (useErrBkgGaus1) {
		bkgErrPdfList.add(bkgTimeErrGaus1);
		if (static_cast<int>(bkgErrPdfList.getSize()) < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrCore1Frac);
	}
	if (useErrBkgGaus2) bkgErrPdfList.add(bkgTimeErrGaus2);
	RooAddPdf bkgErrPdf("pdfErr", "pdfErr", bkgErrPdfList, bkgErrFracList, true);

	RooRealVar sigTimeErrGaus1Mean("sigErrGaus1Mean", "sigErrGaus1Mean", readErrStoredValue("sigErrGaus1Mean", signalErrResult, 0.02), errFitLow, errUpper);
	RooRealVar sigTimeErrGaus1Sigma("sigErrGaus1Sigma", "sigErrGaus1Sigma", readErrStoredValue("sigErrGaus1Sigma", signalErrResult, 0.003), 1e-6, errUpper);
	RooRealVar sigTimeErrGaus2Mean("sigErrGaus2Mean", "sigErrGaus2Mean", readErrStoredValue("sigErrGaus2Mean", signalErrResult, sigTimeErrGaus1Mean.getVal()), errFitLow, errUpper);
	RooRealVar sigTimeErrGaus2Sigma("sigErrGaus2Sigma", "sigErrGaus2Sigma", readErrStoredValue("sigErrGaus2Sigma", signalErrResult, sigTimeErrGaus1Sigma.getVal()), 1e-6, errUpper);
	RooRealVar sigTimeErrTailMpv("sigErrTailMpv", "sigErrTailMpv", readErrStoredValue("sigErrTailMpv", signalErrResult, 0.02), errFitLow, errUpper);
	RooRealVar sigTimeErrTailWidth("sigErrTailWidth", "sigErrTailWidth", readErrStoredValue("sigErrTailWidth", signalErrResult, 0.002), 1e-6, errUpper);
	RooRealVar sigTimeErrTail2Mpv("sigErrTail2Mpv", "sigErrTail2Mpv", readErrStoredValue("sigErrTail2Mpv", signalErrResult, sigTimeErrTailMpv.getVal()), errFitLow, errUpper);
	RooRealVar sigTimeErrTail2Width("sigErrTail2Width", "sigErrTail2Width", readErrStoredValue("sigErrTail2Width", signalErrResult, sigTimeErrTailWidth.getVal()), 1e-6, errUpper);
	RooRealVar sigTimeErrLognM0("sigErrLognM0", "sigErrLognM0", readErrStoredValue("sigErrLognM0", signalErrResult, 0.02), errFitLow, errUpper);
	RooRealVar sigTimeErrLognK("sigErrLognK", "sigErrLognK", readErrStoredValue("sigErrLognK", signalErrResult, 0.35), 0.01, 3.0);
	RooRealVar sigTimeErrCore1FracRatio("sigErrCore1FracRatio", "sigErrCore1FracRatio", readErrStoredValue("sigErrCore1FracRatio", signalErrResult, 1.5), 1e-3, 1e3);
	RooFormulaVar sigTimeErrCore1Frac("sigErrCore1Frac", "@0/(1.0+@0)", RooArgList(sigTimeErrCore1FracRatio));
	RooRealVar sigTimeErrTailFracRatio("sigErrTailFracRatio", "sigErrTailFracRatio", readErrStoredValue("sigErrTailFracRatio", signalErrResult, 0.08 / 0.92), 1e-4, 1e4);
	RooFormulaVar sigTimeErrTailFrac("sigErrTailFrac", "@0/(1.0+@0)", RooArgList(sigTimeErrTailFracRatio));
	RooRealVar sigTimeErrTail2FracRatio("sigErrTail2FracRatio", "sigErrTail2FracRatio", readErrStoredValue("sigErrTail2FracRatio", signalErrResult, 0.05 / 0.95), 1e-4, 1e4);
	RooFormulaVar sigTimeErrTail2Frac("sigErrTail2Frac", "@0/(1.0+@0)", RooArgList(sigTimeErrTail2FracRatio));
	RooRealVar sigTimeErrLognFracRatio("sigErrLognFracRatio", "sigErrLognFracRatio", readErrStoredValue("sigErrLognFracRatio", signalErrResult, 0.03 / 0.97), 1e-4, 1e4);
	RooFormulaVar sigTimeErrLognFrac("sigErrLognFrac", "@0/(1.0+@0)", RooArgList(sigTimeErrLognFracRatio));
	RooGaussian sigTimeErrGaus1("gausSig", "gausSig", obs_timeErr, sigTimeErrGaus1Mean, sigTimeErrGaus1Sigma);
	RooGaussian sigTimeErrGaus2("gausSig2", "gausSig2", obs_timeErr, sigTimeErrGaus2Mean, sigTimeErrGaus2Sigma);
	RooLandau sigTimeErrTail("landauSig", "landauSig", obs_timeErr, sigTimeErrTailMpv, sigTimeErrTailWidth);
	RooLandau sigTimeErrTail2("landauSig2", "landauSig2", obs_timeErr, sigTimeErrTail2Mpv, sigTimeErrTail2Width);
	RooLognormal sigTimeErrLogn("lognSig", "lognSig", obs_timeErr, sigTimeErrLognM0, sigTimeErrLognK);
	RooArgList sigErrPdfList;
	RooArgList sigErrFracList;
	if (useErrSigLandau1) {
		sigErrPdfList.add(sigTimeErrTail);
		if (nErrSigComponents > 1) sigErrFracList.add(sigTimeErrTailFrac);
	}
	if (useErrSigLandau2) {
		sigErrPdfList.add(sigTimeErrTail2);
		if (static_cast<int>(sigErrPdfList.getSize()) < nErrSigComponents) sigErrFracList.add(sigTimeErrTail2Frac);
	}
	if (useErrSigLogn) {
		sigErrPdfList.add(sigTimeErrLogn);
		if (static_cast<int>(sigErrPdfList.getSize()) < nErrSigComponents) sigErrFracList.add(sigTimeErrLognFrac);
	}
	if (useErrSigGaus1) {
		sigErrPdfList.add(sigTimeErrGaus1);
		if (static_cast<int>(sigErrPdfList.getSize()) < nErrSigComponents) sigErrFracList.add(sigTimeErrCore1Frac);
	}
	if (useErrSigGaus2) sigErrPdfList.add(sigTimeErrGaus2);
	std::unique_ptr<RooAddPdf> signalErrPdf = std::make_unique<RooAddPdf>("pdfErrSig", "pdfErrSig", sigErrPdfList, sigErrFracList, true);

	auto fixVar = [](RooRealVar &var) { var.setConstant(true); };
	fixVar(bkgTimeErrGaus1Mean); fixVar(bkgTimeErrGaus1Sigma); fixVar(bkgTimeErrGaus2Mean); fixVar(bkgTimeErrGaus2Sigma);
	fixVar(bkgTimeErrTailMpv); fixVar(bkgTimeErrTailWidth); fixVar(bkgTimeErrTail2Mpv); fixVar(bkgTimeErrTail2Width);
	fixVar(bkgTimeErrLognM0); fixVar(bkgTimeErrLognK); fixVar(bkgTimeErrCore1FracRatio); fixVar(bkgTimeErrTailFracRatio);
	fixVar(bkgTimeErrTail2FracRatio); fixVar(bkgTimeErrLognFracRatio);
	fixVar(sigTimeErrGaus1Mean); fixVar(sigTimeErrGaus1Sigma); fixVar(sigTimeErrGaus2Mean); fixVar(sigTimeErrGaus2Sigma);
	fixVar(sigTimeErrTailMpv); fixVar(sigTimeErrTailWidth); fixVar(sigTimeErrTail2Mpv); fixVar(sigTimeErrTail2Width);
	fixVar(sigTimeErrLognM0); fixVar(sigTimeErrLognK); fixVar(sigTimeErrCore1FracRatio); fixVar(sigTimeErrTailFracRatio);
	fixVar(sigTimeErrTail2FracRatio); fixVar(sigTimeErrLognFracRatio);

	auto *bkgTimeResult = dynamic_cast<RooFitResult *>(bkgFile->Get("timeResult"));
	if (!bkgTimeResult) {
		std::cerr << "ERROR: timeResult not found in " << bkgFileName << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// mass model from mass.C
	// ------------------------------------------------------------------
	const int nSignalGaussComponents = std::clamp(readIntParam(*massFile, "nSignalGaussComponents", 0), 0, 2);
	const int nSignalCBComponents = std::clamp(readIntParam(*massFile, "nSignalCBComponents", 0), 0, 2);
	const int nBkgExpComponents = std::clamp(readIntParam(*massFile, "nBkgExpComponents", 0), 0, 1);
	const int nBkgChebyOrder = std::clamp(readIntParam(*massFile, "nBkgChebyOrder", 0), 0, 6);
	if (nSignalGaussComponents + nSignalCBComponents <= 0) {
		std::cerr << "ERROR: invalid saved mass signal configuration." << std::endl;
		return;
	}

	RooConstVar signal_mass_mean("signal_mass_mean","signal_mass_mean",readDoubleParam(*massFile,"signal_mass_mean",3.096));
	RooConstVar signal_mass_sigma("signal_mass_sigma","signal_mass_sigma",readDoubleParam(*massFile,"signal_mass_sigma",0.03));
	RooConstVar signal_mass_sigma2("signal_mass_sigma2","signal_mass_sigma2",readDoubleParam(*massFile,"signal_mass_sigma2",0.06));
	RooConstVar signal_mass_cb_sigma("signal_mass_cb_sigma","signal_mass_cb_sigma",readDoubleParam(*massFile,"signal_mass_cb_sigma",0.035));
	RooConstVar signal_mass_cb_sigma2("signal_mass_cb_sigma2","signal_mass_cb_sigma2",readDoubleParam(*massFile,"signal_mass_cb_sigma2",0.055));
	RooConstVar signal_mass_cb_alpha("signal_mass_cb_alpha","signal_mass_cb_alpha",readDoubleParam(*massFile,"signal_mass_cb_alpha",1.5));
	RooConstVar signal_mass_cb_alpha2("signal_mass_cb_alpha2","signal_mass_cb_alpha2",readDoubleParam(*massFile,"signal_mass_cb_alpha2",2.0));
	RooConstVar signal_mass_cb_n("signal_mass_cb_n","signal_mass_cb_n",readDoubleParam(*massFile,"signal_mass_cb_n",3.0));
	RooConstVar signal_mass_cb_n2("signal_mass_cb_n2","signal_mass_cb_n2",readDoubleParam(*massFile,"signal_mass_cb_n2",4.0));
	RooConstVar signal_mass_frac1("signal_mass_frac1","signal_mass_frac1",readDoubleParam(*massFile,"signal_mass_frac1",0.65));
	RooConstVar signal_mass_frac2("signal_mass_frac2","signal_mass_frac2",readDoubleParam(*massFile,"signal_mass_frac2",0.20));
	RooConstVar signal_mass_frac3("signal_mass_frac3","signal_mass_frac3",readDoubleParam(*massFile,"signal_mass_frac3",0.10));

	RooGaussian signal_mass_gaus("signal_mass_gaus","signal_mass_gaus",obs_mass,signal_mass_mean,signal_mass_sigma);
	RooGaussian signal_mass_gaus2("signal_mass_gaus2","signal_mass_gaus2",obs_mass,signal_mass_mean,signal_mass_sigma2);
	std::unique_ptr<RooAbsPdf> signal_mass_cb_owned;
	if (nSignalCBComponents >= 2)
	{
		signal_mass_cb_owned = std::make_unique<RooCrystalBall>(
			"signal_mass_cb","signal_mass_cb",
			obs_mass,signal_mass_mean,
			signal_mass_cb_sigma,signal_mass_cb_sigma2,
			signal_mass_cb_alpha,signal_mass_cb_n,
			signal_mass_cb_alpha2,signal_mass_cb_n2);
	}
	else if (nSignalCBComponents >= 1)
	{
		signal_mass_cb_owned = std::make_unique<RooCBShape>(
			"signal_mass_cb","signal_mass_cb",
			obs_mass,signal_mass_mean,signal_mass_cb_sigma,signal_mass_cb_alpha,signal_mass_cb_n);
	}
	RooAbsPdf *signal_mass_cb = signal_mass_cb_owned.get();

	RooArgList signalMassPdfList;
	if (nSignalGaussComponents >= 1) signalMassPdfList.add(signal_mass_gaus);
	if (nSignalCBComponents >= 1) signalMassPdfList.add(*signal_mass_cb);
	if (nSignalGaussComponents >= 2) signalMassPdfList.add(signal_mass_gaus2);
	RooArgList signalMassFracList;
	if (signalMassPdfList.getSize() >= 2) signalMassFracList.add(signal_mass_frac1);
	if (signalMassPdfList.getSize() >= 3) signalMassFracList.add(signal_mass_frac2);
	if (signalMassPdfList.getSize() >= 4) signalMassFracList.add(signal_mass_frac3);
	std::unique_ptr<RooAbsPdf> signal_mass_owned;
	RooAbsPdf *signal_mass_pdf = nullptr;
	if (signalMassPdfList.getSize() == 1 && nSignalGaussComponents >= 1) signal_mass_pdf = static_cast<RooAbsPdf *>(&signal_mass_gaus);
	if (signalMassPdfList.getSize() == 1 && nSignalGaussComponents == 0) signal_mass_pdf = signal_mass_cb;
	if (!signal_mass_pdf) {
		signal_mass_owned = std::make_unique<RooAddPdf>("signal_mass","signal_mass",signalMassPdfList,signalMassFracList,true);
		signal_mass_pdf = signal_mass_owned.get();
	}

	RooConstVar bkg_mass_lambda("bkg_mass_lambda","bkg_mass_lambda",readDoubleParam(*massFile,"bkg_mass_lambda",-1.0));
	RooConstVar bkg_mass_p1("bkg_mass_p1","bkg_mass_p1",readDoubleParam(*massFile,"bkg_mass_p1",0.0));
	RooConstVar bkg_mass_p2("bkg_mass_p2","bkg_mass_p2",readDoubleParam(*massFile,"bkg_mass_p2",0.0));
	RooConstVar bkg_mass_p3("bkg_mass_p3","bkg_mass_p3",readDoubleParam(*massFile,"bkg_mass_p3",0.0));
	RooConstVar bkg_mass_p4("bkg_mass_p4","bkg_mass_p4",readDoubleParam(*massFile,"bkg_mass_p4",0.0));
	RooConstVar bkg_mass_p5("bkg_mass_p5","bkg_mass_p5",readDoubleParam(*massFile,"bkg_mass_p5",0.0));
	RooConstVar bkg_mass_p6("bkg_mass_p6","bkg_mass_p6",readDoubleParam(*massFile,"bkg_mass_p6",0.0));
	RooArgList chebyCoeffList;
	if (nBkgChebyOrder >= 1) chebyCoeffList.add(bkg_mass_p1);
	if (nBkgChebyOrder >= 2) chebyCoeffList.add(bkg_mass_p2);
	if (nBkgChebyOrder >= 3) chebyCoeffList.add(bkg_mass_p3);
	if (nBkgChebyOrder >= 4) chebyCoeffList.add(bkg_mass_p4);
	if (nBkgChebyOrder >= 5) chebyCoeffList.add(bkg_mass_p5);
	if (nBkgChebyOrder >= 6) chebyCoeffList.add(bkg_mass_p6);
	std::unique_ptr<RooAbsPdf> bkg_mass_owned;
	if (nBkgExpComponents == 1)
		bkg_mass_owned = std::make_unique<RooExponential>("bkg_mass","bkg_mass",obs_mass,bkg_mass_lambda);
	else
		bkg_mass_owned = std::make_unique<RooChebychev>("bkg_mass","bkg_mass",obs_mass,chebyCoeffList);
	RooAbsPdf *bkg_mass_pdf = bkg_mass_owned.get();

	// ------------------------------------------------------------------
	// ctau prompt resolution from ctau_pr.C
	// ------------------------------------------------------------------
	const int nResolutionComponents = std::clamp(readIntParam(*prFile, "nResolutionComponents", 1), 1, 4);
	RooConstVar ctauMeanScale("ctauMeanScale","ctauMeanScale",readDoubleParam(*prFile,"ctauMeanScale",1.0));
	RooConstVar ctauTime1Mean("ctauTime1Mean","ctauTime1Mean",readDoubleParam(*prFile,"ctauTime1Mean",0.0));
	const double ctauTime1ScaleVal = readDoubleParam(*prFile,"ctauTime1Scale",1.0);
	RooRealVar ctauTime1Scale("ctauTime1Scale","ctauTime1Scale",
		ctauTime1ScaleVal,
		std::max(0.3, 0.5 * ctauTime1ScaleVal),
		std::max(3.0, 2.0 * ctauTime1ScaleVal));
	RooConstVar ctauTime2Mean("ctauTime2Mean","ctauTime2Mean",readDoubleParam(*prFile,"ctauTime2Mean",0.0));
	const double ctauTime2ScaleVal = readDoubleParam(*prFile,"ctauTime2Scale",1.0);
	const double ctauTime2DeltaVal = std::max(0.001, ctauTime2ScaleVal - ctauTime1ScaleVal);
	RooRealVar ctauTime2Delta("ctauTime2Delta","ctauTime2Delta",
		ctauTime2DeltaVal,
		0.001,
		std::max(3.0, 2.0 * ctauTime2DeltaVal));
	RooFormulaVar ctauTime2Scale("ctauTime2Scale","@0+@1",RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooConstVar ctauTime3Mean("ctauTime3Mean","ctauTime3Mean",readDoubleParam(*prFile,"ctauTime3Mean",0.0));
	const double ctauTime3ScaleVal = readDoubleParam(*prFile,"ctauTime3Scale",1.0);
	const double ctauTime3DeltaVal = std::max(0.001, ctauTime3ScaleVal - ctauTime2ScaleVal);
	RooRealVar ctauTime3Delta("ctauTime3Delta","ctauTime3Delta",
		ctauTime3DeltaVal,
		0.001,
		std::max(3.0, 2.0 * ctauTime3DeltaVal));
	RooFormulaVar ctauTime3Scale("ctauTime3Scale","@0+@1",RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooConstVar ctauTime4Mean("ctauTime4Mean","ctauTime4Mean",readDoubleParam(*prFile,"ctauTime4Mean",0.0));
	const double ctauTime4ScaleVal = readDoubleParam(*prFile,"ctauTime4Scale",1.0);
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	RooRealVar ctauTime4Delta("ctauTime4Delta","ctauTime4Delta",
		ctauTime4DeltaVal,
		0.05,
		std::max(3.0, 2.0 * ctauTime4DeltaVal));
	RooFormulaVar ctauTime4Scale("ctauTime4Scale","@0+@1",RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooConstVar ctauFrac1("ctauFrac1","ctauFrac1",readDoubleParam(*prFile,"ctauFrac1",0.5));
	RooConstVar ctauFrac2("ctauFrac2","ctauFrac2",readDoubleParam(*prFile,"ctauFrac2",0.2));
	RooConstVar ctauFrac3("ctauFrac3","ctauFrac3",readDoubleParam(*prFile,"ctauFrac3",0.1));

	RooGaussModel ctauTime1("ctauTime1","ctauTime1",obs_time,ctauTime1Mean,ctauTime1Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime2("ctauTime2","ctauTime2",obs_time,ctauTime2Mean,ctauTime2Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime3("ctauTime3","ctauTime3",obs_time,ctauTime3Mean,ctauTime3Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime4("ctauTime4","ctauTime4",obs_time,ctauTime4Mean,ctauTime4Scale,ctauMeanScale,obs_timeErr);
	double maxPromptScaleSaved = std::max(ctauTime1ScaleVal, 0.0);
	if (nResolutionComponents >= 2 && std::isfinite(ctauTime2ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime2ScaleVal);
	if (nResolutionComponents >= 3 && std::isfinite(ctauTime3ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime3ScaleVal);
	if (nResolutionComponents >= 4 && std::isfinite(ctauTime4ScaleVal))
		maxPromptScaleSaved = std::max(maxPromptScaleSaved, ctauTime4ScaleVal);
	const double resolutionDrivenLifetimeFloor =
		std::max(absoluteLifetimeFloor, lifetimeFloorCoeff * maxPromptScaleSaved * std::max(errFitHigh, absoluteLifetimeFloor));
	signalLifetimeFloor = std::max(signalLifetimeFloor, resolutionDrivenLifetimeFloor);
	bkgLifetimeFloor = std::max(bkgLifetimeFloor, resolutionDrivenLifetimeFloor);
	maxStableLifetime = std::max(maxStableLifetime, 20.0 * bkgLifetimeFloor);
	if (fixResolutionParams)
	{
		ctauTime1Scale.setVal(ctauTime1ScaleVal);
		ctauTime1Scale.setConstant(true);
		ctauTime2Delta.setVal(ctauTime2DeltaVal);
		ctauTime2Delta.setConstant(true);
		ctauTime3Delta.setVal(ctauTime3DeltaVal);
		ctauTime3Delta.setConstant(true);
		ctauTime4Delta.setVal(ctauTime4DeltaVal);
		ctauTime4Delta.setConstant(true);
	}

	std::unique_ptr<RooAddModel> promptResolutionModel;
	RooResolutionModel *promptResolutionPtr = &ctauTime1;
	if (nResolutionComponents >= 2) {
		RooArgList resolutionPdfList;
		RooArgList resolutionFracList;
		resolutionPdfList.add(ctauTime1);
		resolutionPdfList.add(ctauTime2);
		resolutionFracList.add(ctauFrac1);
		if (nResolutionComponents >= 3) {
			resolutionPdfList.add(ctauTime3);
			resolutionFracList.add(ctauFrac2);
		}
		if (nResolutionComponents >= 4) {
			resolutionPdfList.add(ctauTime4);
			resolutionFracList.add(ctauFrac3);
		}
		promptResolutionModel = std::make_unique<RooAddModel>("promptResolutionModel","promptResolutionModel",resolutionPdfList,resolutionFracList);
		promptResolutionPtr = promptResolutionModel.get();
	}
	RooResolutionModel &promptResolution = *promptResolutionPtr;
	RooAbsPdf *promptTimePdf = promptResolutionPtr;

	// ------------------------------------------------------------------
	// nonprompt signal model from ctau_np.C
	// ------------------------------------------------------------------
	const int nSignalSSComponentsSaved = std::clamp(readIntParam(*npFile, "nSignalSSComponents", 1), 1, 3);
	const double npYield1 = readDoubleParam(*npFile, "NsignalSS1", 1.0);
	const double npYield2 = readDoubleParam(*npFile, "NsignalSS2", 0.0);
	const double npYield3 = readDoubleParam(*npFile, "NsignalSS3", 0.0);
	const double signalLifetime1Val = readDoubleParam(*npFile,"signal_lifetime",0.1);
	const double signalLifetime2Val = readDoubleParam(*npFile,"signal2_lifetime",0.2);
	const double signalLifetime3Val = readDoubleParam(*npFile,"signal3_lifetime",0.4);
	const double signalLifetimeFloorSaved = std::max(
		readDoubleParam(*npFile, "signalLifetimeFloor", signalLifetimeFloor),
		signalLifetimeFloor);
	const double signalLifetime2FloorSaved = std::max(
		readDoubleParam(*npFile, "signalLifetime2Floor", std::max(0.75 * signalLifetimeFloorSaved, 1e-2)),
		std::max(0.75 * signalLifetimeFloorSaved, 1e-2));
	const double signalLifetime3FloorSaved = std::max(
		readDoubleParam(*npFile, "signalLifetime3Floor", std::max(1.00 * signalLifetimeFloorSaved, 1e-2)),
		std::max(1.00 * signalLifetimeFloorSaved, 1e-2));
	const bool useSignalSS2 = (nSignalSSComponentsSaved >= 2) && std::isfinite(signalLifetime2Val) && signalLifetime2Val > signalLifetime2FloorSaved;
	const bool useSignalSS3 = useSignalSS2 && (nSignalSSComponentsSaved >= 3) && std::isfinite(signalLifetime3Val) && signalLifetime3Val > signalLifetime3FloorSaved;
	const int nSignalSSComponents = 1 + (useSignalSS2 ? 1 : 0) + (useSignalSS3 ? 1 : 0);
	const double signalLifetimeInit = std::max(std::max(signalLifetime1Val, signalLifetimeFloorSaved + 1e-3), std::max(0.08, 1.5 * signalLifetimeFloorSaved));
	const double signalLifetimeCeil = std::max(signalLifetimeInit * 5.0, maxStableLifetime);
	RooRealVar signal_lifetime(
		"signal_lifetime", "signal_lifetime",
		signalLifetimeInit,
		signalLifetimeFloorSaved + 1e-2,
		signalLifetimeCeil);
	const double signalLifetime2Floor = signalLifetime2FloorSaved;
	const double signalLifetime2Init = std::max(std::max(signalLifetime2Val, signalLifetime2Floor + 1e-3), std::max(0.20, 2.0 * signalLifetime2Floor));
	const double signalLifetime2Ceil = std::max(signalLifetime2Init * 5.0, maxStableLifetime);
	RooRealVar signal2_lifetime(
		"signal2_lifetime", "signal2_lifetime",
		signalLifetime2Init,
		signalLifetime2Floor + 2e-2,
		signalLifetime2Ceil);
	const double signalLifetime3Floor = signalLifetime3FloorSaved;
	const double signalLifetime3Init = std::max(std::max(signalLifetime3Val, signalLifetime3Floor + 1e-3), std::max(0.40, 2.0 * signalLifetime3Floor));
	const double signalLifetime3Ceil = std::max(signalLifetime3Init * 5.0, maxStableLifetime);
	RooRealVar signal3_lifetime(
		"signal3_lifetime", "signal3_lifetime",
		signalLifetime3Init,
		signalLifetime3Floor + 1e-2,
		signalLifetime3Ceil);
	const double npYieldSum = std::max(1e-6, npYield1 + npYield2 + npYield3);
	RooRealVar NsignalSS1("NsignalSS1","NsignalSS1",
		std::max(npYield1, 1e-6), 0.0, std::max(5.0 * npYieldSum, 10.0));
	RooRealVar NsignalSS2("NsignalSS2","NsignalSS2",
		std::max(npYield2, 1e-6), 0.0, std::max(5.0 * npYieldSum, 10.0));
	RooRealVar NsignalSS3("NsignalSS3","NsignalSS3",
		std::max(npYield3, 1e-6), 0.0, std::max(5.0 * npYieldSum, 10.0));
	RooDecay signal_ss1_time("signal_ss1_time","signal_ss1_time",obs_time,signal_lifetime,promptResolution,RooDecay::SingleSided);
	RooDecay signal_ss2_time("signal_ss2_time","signal_ss2_time",obs_time,signal2_lifetime,promptResolution,RooDecay::SingleSided);
	RooDecay signal_ss3_time("signal_ss3_time","signal_ss3_time",obs_time,signal3_lifetime,promptResolution,RooDecay::SingleSided);
	std::unique_ptr<RooAddPdf> signal_np_owned;
	RooAbsPdf *signal_np_time = nullptr;
	if (nSignalSSComponents == 1) {
		signal_np_time = &signal_ss1_time;
	} else if (nSignalSSComponents == 2) {
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time","signal_np_time",RooArgList(signal_ss1_time,signal_ss2_time),RooArgList(NsignalSS1,NsignalSS2));
		signal_np_time = signal_np_owned.get();
	} else {
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time","signal_np_time",RooArgList(signal_ss1_time,signal_ss2_time,signal_ss3_time),RooArgList(NsignalSS1,NsignalSS2,NsignalSS3));
		signal_np_time = signal_np_owned.get();
	}

	RooRealVar bFraction("bFraction","bFraction",0.15,0.0,1.0);
	RooAddPdf signal_time("signal_time","signal_time",RooArgList(*signal_np_time,*promptTimePdf),RooArgList(bFraction),true);

	// ------------------------------------------------------------------
	// background ctau model from ctau_bkg.C
	// ------------------------------------------------------------------
	const int nBkgSSComponents = std::clamp(readIntParam(*bkgFile, "nSignalSSComponents", 0), 0, 3);
	const int nBkgFlipComponents = std::clamp(readIntParam(*bkgFile, "nSignalFlipComponents", 0), 0, 3);
	const int nBkgDSComponents = std::clamp(readIntParam(*bkgFile, "nSignalDSComponents", 0), 0, 3);

	RooRealVar bkg_lifetime("bkg_lifetime","bkg_lifetime",
		readBkgResultValue(*bkgTimeResult,"bkg_ss_lifetime","bkg_lifetime",std::max(0.05, 1.5 * bkgLifetimeFloor)));
	RooRealVar bkg_lifetime2("bkg_lifetime2","bkg_lifetime2",
		readBkgResultValue(*bkgTimeResult,"bkg_ss_lifetime2","bkg_lifetime2",std::max(0.04, 1.5 * std::max(1.5 * bkgLifetimeFloor, 2e-2))));
	RooRealVar bkg_lifetime3("bkg_lifetime3","bkg_lifetime3",
		readBkgResultValue(*bkgTimeResult,"bkg_ss_lifetime3","bkg_lifetime3",std::max(0.10, 1.5 * std::max(2.0 * bkgLifetimeFloor, 3e-2))));
	RooRealVar bkg_sym_lifetime("bkg_sym_lifetime","bkg_sym_lifetime",
		readBkgResultValue(*bkgTimeResult,"bkg_ds_lifetime","bkg_sym_lifetime",std::max(0.04, 1.4 * std::max(1.5 * errRange.first, 1e-2))));
	RooRealVar bkg_sym_lifetime2("bkg_sym_lifetime2","bkg_sym_lifetime2",
		readBkgResultValue(*bkgTimeResult,"bkg_ds_lifetime2","bkg_sym_lifetime2",std::max(0.07, 1.4 * std::max(2.0 * std::max(1.5 * errRange.first, 1e-2), 2e-2))));
	RooRealVar bkg_sym_lifetime3("bkg_sym_lifetime3","bkg_sym_lifetime3",
		readBkgResultValue(*bkgTimeResult,"bkg_ds_lifetime3","bkg_sym_lifetime3",std::max(0.12, 1.3 * std::max(3.0 * std::max(1.5 * errRange.first, 1e-2), 4e-2))));
	RooRealVar bkg_flip_lifetime("bkg_flip_lifetime","bkg_flip_lifetime",
		readBkgResultValue(*bkgTimeResult,"bkg_flip_lifetime","",std::max(0.05, 1.2 * std::max(0.75 * bkgLifetimeFloor, 5e-3))));
	RooRealVar bkg_flip_lifetime2("bkg_flip_lifetime2","bkg_flip_lifetime2",
		readBkgResultValue(*bkgTimeResult,"bkg_flip_lifetime2","",std::max(0.08, 1.3 * std::max(1.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 1e-2))));
	RooRealVar bkg_flip_lifetime3("bkg_flip_lifetime3","bkg_flip_lifetime3",
		readBkgResultValue(*bkgTimeResult,"bkg_flip_lifetime3","",std::max(0.12, 1.2 * std::max(2.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 2e-2))));
	bkg_lifetime.setConstant(false);
	bkg_lifetime2.setConstant(false);
	bkg_lifetime3.setConstant(false);
	bkg_sym_lifetime.setConstant(false);
	bkg_sym_lifetime2.setConstant(false);
	bkg_sym_lifetime3.setConstant(false);
	bkg_flip_lifetime.setConstant(false);
	bkg_flip_lifetime2.setConstant(false);
	bkg_flip_lifetime3.setConstant(false);

	const double bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 2e-2);
	const double bkgLifetime3Floor = std::max(2.0 * bkgLifetimeFloor, 3e-2);
	const double bkgSymLifetimeFloor = std::max(1.5 * errRange.first, 1e-2);
	const double bkgSymLifetime2Floor = std::max(2.0 * bkgSymLifetimeFloor, 2e-2);
	const double bkgSymLifetime3Floor = std::max(3.0 * bkgSymLifetimeFloor, 4e-2);
	const double bkgFlipLifetimeFloor = std::max(0.75 * bkgLifetimeFloor, 5e-3);
	const double bkgFlipLifetime2Floor = std::max(1.5 * bkgFlipLifetimeFloor, 1e-2);
	const double bkgFlipLifetime3Floor = std::max(2.5 * bkgFlipLifetimeFloor, 2e-2);
	const double bkgSymLifetimeMinGap = std::max(0.30 * bkgSymLifetimeFloor, 2e-3);
	const double bkgSymLifetime2MinGap = std::max(0.25 * bkgSymLifetime2Floor, 3e-3);
	const double bkgSymLifetime3MinGap = std::max(0.20 * bkgSymLifetime3Floor, 4e-3);
	bkg_lifetime.setRange(bkgLifetimeFloor + 1e-2, std::max(0.50, maxStableLifetime));
	bkg_lifetime2.setRange(bkgLifetime2Floor + 2e-2, std::max(0.50, maxStableLifetime));
	bkg_lifetime3.setRange(bkgLifetime3Floor + 1e-2, std::max(0.80, maxStableLifetime));
	bkg_sym_lifetime.setRange(bkgSymLifetimeFloor + bkgSymLifetimeMinGap, std::max(0.30, 0.5 * maxStableLifetime));
	bkg_sym_lifetime2.setRange(bkgSymLifetime2Floor + bkgSymLifetime2MinGap, std::max(0.50, maxStableLifetime));
	bkg_sym_lifetime3.setRange(bkgSymLifetime3Floor + bkgSymLifetime3MinGap, std::max(0.80, maxStableLifetime));
	bkg_flip_lifetime.setRange(bkgFlipLifetimeFloor + 1e-2, std::max(0.50, maxStableLifetime));
	bkg_flip_lifetime2.setRange(bkgFlipLifetime2Floor + 2e-2, std::max(0.60, maxStableLifetime));
	bkg_flip_lifetime3.setRange(bkgFlipLifetime3Floor + 1e-2, std::max(0.80, maxStableLifetime));

	RooDecay bkg_time_ss("bkg_time_ss","bkg_time_ss",obs_time,bkg_lifetime,promptResolution,RooDecay::SingleSided);
	RooDecay bkg_time_ss2("bkg_time_ss2","bkg_time_ss2",obs_time,bkg_lifetime2,promptResolution,RooDecay::SingleSided);
	RooDecay bkg_time_ss3("bkg_time_ss3","bkg_time_ss3",obs_time,bkg_lifetime3,promptResolution,RooDecay::SingleSided);
	RooDecay bkg_time_ds("bkg_time_ds","bkg_time_ds",obs_time,bkg_sym_lifetime,promptResolution,RooDecay::DoubleSided);
	RooDecay bkg_time_ds2("bkg_time_ds2","bkg_time_ds2",obs_time,bkg_sym_lifetime2,promptResolution,RooDecay::DoubleSided);
	RooDecay bkg_time_ds3("bkg_time_ds3","bkg_time_ds3",obs_time,bkg_sym_lifetime3,promptResolution,RooDecay::DoubleSided);
	RooDecay bkg_time_flip("bkg_time_flip","bkg_time_flip",obs_time,bkg_flip_lifetime,promptResolution,RooDecay::Flipped);
	RooDecay bkg_time_flip2("bkg_time_flip2","bkg_time_flip2",obs_time,bkg_flip_lifetime2,promptResolution,RooDecay::Flipped);
	RooDecay bkg_time_flip3("bkg_time_flip3","bkg_time_flip3",obs_time,bkg_flip_lifetime3,promptResolution,RooDecay::Flipped);

	const double NplbVal = readResultValue(*bkgTimeResult, "Nplb", 1.0);
	const double Nss1Val = readResultValue(*bkgTimeResult, "Nss1", 0.0);
	const double Nss2Val = readResultValue(*bkgTimeResult, "Nss2", 0.0);
	const double Nss3Val = readResultValue(*bkgTimeResult, "Nss3", 0.0);
	const double Nds1Val = readResultValue(*bkgTimeResult, "Nds1", 0.0);
	const double Nds2Val = readResultValue(*bkgTimeResult, "Nds2", 0.0);
	const double Nds3Val = readResultValue(*bkgTimeResult, "Nds3", 0.0);
	const double Nflip1Val = readResultValue(*bkgTimeResult, "Nflip1", 0.0);
	const double Nflip2Val = readResultValue(*bkgTimeResult, "Nflip2", 0.0);
	const double Nflip3Val = readResultValue(*bkgTimeResult, "Nflip3", 0.0);
	const double bkgYieldSum = std::max(1e-12,
		NplbVal +
		(nBkgSSComponents >= 1 ? Nss1Val : 0.0) +
		(nBkgSSComponents >= 2 ? Nss2Val : 0.0) +
		(nBkgSSComponents >= 3 ? Nss3Val : 0.0) +
		(nBkgDSComponents >= 1 ? Nds1Val : 0.0) +
		(nBkgDSComponents >= 2 ? Nds2Val : 0.0) +
		(nBkgDSComponents >= 3 ? Nds3Val : 0.0) +
		(nBkgFlipComponents >= 1 ? Nflip1Val : 0.0) +
		(nBkgFlipComponents >= 2 ? Nflip2Val : 0.0) +
		(nBkgFlipComponents >= 3 ? Nflip3Val : 0.0));

	std::vector<std::unique_ptr<RooRealVar>> bkgCoeffVars;
	RooArgList bkgTimePdfList;
	RooArgList bkgTimeCoeffList;
	auto addBkgComponent = [&](RooAbsPdf &pdf, const char *name, double yield) {
		bkgTimePdfList.add(pdf);
		bkgCoeffVars.push_back(std::make_unique<RooRealVar>(
			name, name, std::max(yield, 1e-6), 0.0, std::max(5.0 * bkgYieldSum, 10.0)));
		bkgTimeCoeffList.add(*bkgCoeffVars.back());
	};
	addBkgComponent(*promptTimePdf, "Nplb", NplbVal);
	if (nBkgSSComponents >= 1) addBkgComponent(bkg_time_ss, "Nss1", Nss1Val);
	if (nBkgSSComponents >= 2) addBkgComponent(bkg_time_ss2, "Nss2", Nss2Val);
	if (nBkgSSComponents >= 3) addBkgComponent(bkg_time_ss3, "Nss3", Nss3Val);
	if (nBkgDSComponents >= 1) addBkgComponent(bkg_time_ds, "Nds1", Nds1Val);
	if (nBkgDSComponents >= 2) addBkgComponent(bkg_time_ds2, "Nds2", Nds2Val);
	if (nBkgDSComponents >= 3) addBkgComponent(bkg_time_ds3, "Nds3", Nds3Val);
	if (nBkgFlipComponents >= 1) addBkgComponent(bkg_time_flip, "Nflip1", Nflip1Val);
	if (nBkgFlipComponents >= 2) addBkgComponent(bkg_time_flip2, "Nflip2", Nflip2Val);
	if (nBkgFlipComponents >= 3) addBkgComponent(bkg_time_flip3, "Nflip3", Nflip3Val);

	std::unique_ptr<RooAddPdf> bkg_time_owned;
	RooAbsPdf *bkg_time_pdf = nullptr;
	if (bkgTimePdfList.getSize() == 1) {
		bkg_time_pdf = promptTimePdf;
	} else {
		bkg_time_owned = std::make_unique<RooAddPdf>("bkg_time","bkg_time",bkgTimePdfList,bkgTimeCoeffList);
		bkg_time_pdf = bkg_time_owned.get();
	}

	// ------------------------------------------------------------------
	// full 2D fit: all shapes fixed, only Nsig/Nbkg and bFraction free
	// ------------------------------------------------------------------
	RooRealVar Nsig("Nsig","Nsig",readDoubleParam(*massFile,"Nsig",std::max(1.0, 0.5 * data->numEntries())),0.0,2.0 * data->numEntries());
	RooRealVar Nbkg("Nbkg","Nbkg",readDoubleParam(*massFile,"Nbkg",std::max(1.0, 0.5 * data->numEntries())),0.0,2.0 * data->numEntries());
	Nsig.setConstant(true);
	Nbkg.setConstant(true);

	auto floor_empty_template_bins = [](TH1 *hist)
	{
		if (!hist)
			return;
		double minPositive = std::numeric_limits<double>::infinity();
		double total = 0.0;
		for (int i = 1; i <= hist->GetNbinsX(); ++i)
		{
			const double content = hist->GetBinContent(i);
			if (content > 0.0)
			{
				minPositive = std::min(minPositive, content);
				total += content;
			}
		}
		const double floor = std::max({0.2, std::isfinite(minPositive) ? 0.1 * minPositive : 0.0, 1e-7 * total});
		for (int i = 1; i <= hist->GetNbinsX(); ++i)
		{
			if (hist->GetBinContent(i) <= 0.0)
			{
				hist->SetBinContent(i, floor);
				hist->SetBinError(i, 0.0);
			}
		}
	};
	auto make_smoothed_template = [&](TH1 *source, const char *name) -> std::unique_ptr<TH1>
	{
		if (!source)
			return nullptr;
		const int rebinFactor = source->GetNbinsX() >= 80 ? 4 : 2;
		std::unique_ptr<TH1> templ;
		if (source->GetNbinsX() % rebinFactor == 0)
			templ = std::unique_ptr<TH1>(source->Rebin(rebinFactor, name));
		else
			templ = std::unique_ptr<TH1>(static_cast<TH1 *>(source->Clone(name)));
		templ->SetDirectory(nullptr);
		templ->Smooth(1);
		floor_empty_template_bins(templ.get());
		return templ;
	};
	std::unique_ptr<RooAbsData> signalErrData;
	std::unique_ptr<TH1> signalErrPdfTemplateHist;
	{
		if (sigErrPdfOpt == kErrPdfHistPrMc) {
			auto hErrSigPrMc = std::unique_ptr<TH1>(prMcErrData->createHistogram(
				"fit2d_hErrSigPrMc", obs_timeErr, Binning(timeErrPlotBins, errFitLow, errFitHigh)));
			if (!hErrSigPrMc) {
				std::cerr << "ERROR: failed to build prompt-MC err histogram" << std::endl;
				return;
			}
			signalErrData = std::make_unique<RooDataHist>("signalErrDataPrMc", "", RooArgSet(obs_timeErr), hErrSigPrMc.get());
		} else {
			const double scaleBkg = std::isfinite(savedScaleBkg) ? savedScaleBkg : 0.0;
			if (scaleBkg <= 0.0) {
				std::cerr << "ERROR: invalid saved scaleBkg in " << errFileName << std::endl;
				return;
			}
			auto hErrSigSR = std::unique_ptr<TH1>(dataSR->createHistogram("fit2d_hErrSigSR", obs_timeErr, Binning(timeErrPlotBins, errFitLow, errFitHigh)));
			auto hErrBkgSB = std::unique_ptr<TH1>(dataSB->createHistogram("fit2d_hErrBkgSB", obs_timeErr, Binning(timeErrPlotBins, errFitLow, errFitHigh)));
			if (!hErrSigSR || !hErrBkgSB) {
				std::cerr << "ERROR: failed to build err histograms for signal-region subtraction" << std::endl;
				return;
			}
			hErrSigSR->Add(hErrBkgSB.get(), -scaleBkg);
			for (int i = 1; i <= hErrSigSR->GetNbinsX(); ++i)
			{
				if (hErrSigSR->GetBinContent(i) < 0.0)
				{
					hErrSigSR->SetBinContent(i, 0.0);
					hErrSigSR->SetBinError(i, 0.0);
				}
			}
			signalErrData = std::make_unique<RooDataHist>("signalErrDataSub", "", RooArgSet(obs_timeErr), hErrSigSR.get());
			if (sigErrPdfOpt == kErrPdfHist) {
				signalErrPdfTemplateHist = make_smoothed_template(hErrSigSR.get(), "fit2d_hErrSigPdfTemplate");
			}
		}
	}
	if (!signalErrData || signalErrData->sumEntries() <= 0.0) {
		std::cerr << "ERROR: empty signal ctau3DErr template/input distribution" << std::endl;
		return;
	}
	std::unique_ptr<RooDataHist> bkgErrHistData;
	std::unique_ptr<TH1> bkgErrHistTemplate;
	std::unique_ptr<RooHistPdf> bkgErrPdfHist;
	RooAbsPdf *bkgErrPdfPtr = &bkgErrPdf;
	if (isHistErrPdfChoice(bkgErrPdfOpt)) {
		bkgErrHistTemplate = std::unique_ptr<TH1>(bkgErrData->createHistogram(
			"fit2d_hErrBkgTemplate", obs_timeErr, Binning(timeErrPlotBins, errFitLow, errFitHigh)));
		if (!bkgErrHistTemplate) {
			std::cerr << "ERROR: failed to build background ctau3DErr histogram template" << std::endl;
			return;
		}
		bkgErrHistData = std::make_unique<RooDataHist>(
			"bkgErrHistData", "bkgErrHistData",
			RooArgSet(obs_timeErr), bkgErrHistTemplate.get());
		bkgErrPdfHist = std::make_unique<RooHistPdf>(
			"bkgErrHistPdf", "bkgErrHistPdf",
			RooArgSet(obs_timeErr), *bkgErrHistData, histPdfInterpolationOrder);
		bkgErrPdfPtr = bkgErrPdfHist.get();
	}
	std::unique_ptr<RooDataHist> signalErrHistTemplate;
	std::unique_ptr<RooHistPdf> signalErrPdfHist;
	RooAbsPdf *signalErrPdfPtr = signalErrPdf.get();
	if (isHistErrPdfChoice(sigErrPdfOpt)) {
		auto *signalErrHist = dynamic_cast<RooDataHist *>(signalErrData.get());
		if (!signalErrHist) {
			std::cerr << "ERROR: signalErrData is not a RooDataHist for RooHistPdf mode" << std::endl;
			return;
		}
		if (signalErrPdfTemplateHist)
			signalErrHistTemplate = std::make_unique<RooDataHist>(
				"signalErrHistTemplate", "signalErrHistTemplate",
				RooArgSet(obs_timeErr), signalErrPdfTemplateHist.get());
		else
			signalErrHistTemplate = std::make_unique<RooDataHist>(*signalErrHist);
		signalErrPdfHist = std::make_unique<RooHistPdf>(
			"signalErrHistPdf", "signalErrHistPdf",
			RooArgSet(obs_timeErr), *signalErrHistTemplate, histPdfInterpolationOrder);
		signalErrPdfPtr = signalErrPdfHist.get();
	}

	std::vector<std::unique_ptr<RooRealVar>> constraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> constraintPdfs;
	RooArgSet constraints;
	auto constraintSigmaFloor = [&](double central, double fitErr, double relFloor, double absFloor)
	{
		double sigma = absFloor;
		if (std::isfinite(central))
			sigma = std::max(sigma, relFloor * std::abs(central));
		if (std::isfinite(fitErr) && fitErr > 0.0)
			sigma = std::max(sigma, fitErr);
		return sigma;
	};
	auto constraintSigmaBounded = [&](double central, double fitErr,
		double relFloor, double absFloor, double relCeil, double absCeil)
	{
		double sigma = constraintSigmaFloor(central, fitErr, relFloor, absFloor);
		double upper = absCeil;
		if (std::isfinite(central))
			upper = std::max(upper, relCeil * std::abs(central));
		if (std::isfinite(upper) && upper > 0.0)
			sigma = std::min(sigma, upper);
		return sigma;
	};
	auto addConstraint = [&](const char *baseName, RooRealVar &var, double central, double sigma)
	{
		if (!(std::isfinite(central) && std::isfinite(sigma) && sigma > 0.0))
			return;
		var.setVal(central);
		const TString meanName = TString::Format("%s_mean", baseName);
		const TString sigmaName = TString::Format("%s_sigma", baseName);
		const TString pdfName = TString::Format("%s_constraint", baseName);
		constraintConsts.push_back(std::make_unique<RooRealVar>(meanName, meanName, central));
		constraintConsts.back()->setConstant(true);
		constraintConsts.push_back(std::make_unique<RooRealVar>(
			sigmaName, sigmaName, sigma, 1e-9, std::max(10.0 * sigma, 1e-8)));
		constraintConsts.back()->setConstant(true);
		constraintPdfs.push_back(std::make_unique<RooGaussian>(pdfName, pdfName, var,
			*constraintConsts[constraintConsts.size() - 2], *constraintConsts.back()));
		constraints.add(*constraintPdfs.back());
	};
	if (!fixResolutionParams)
	{
		addConstraint("ctauTime1Scale", ctauTime1Scale, ctauTime1ScaleVal,
			constraintSigmaFloor(ctauTime1ScaleVal, std::numeric_limits<double>::quiet_NaN(), 0.10, 0.05));
		if (nResolutionComponents >= 2)
			addConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal,
				constraintSigmaFloor(ctauTime2DeltaVal, std::numeric_limits<double>::quiet_NaN(), 0.10, 0.05));
		if (nResolutionComponents >= 3)
			addConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal,
				constraintSigmaFloor(ctauTime3DeltaVal, std::numeric_limits<double>::quiet_NaN(), 0.10, 0.05));
		if (nResolutionComponents >= 4)
			addConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal,
				constraintSigmaFloor(ctauTime4DeltaVal, std::numeric_limits<double>::quiet_NaN(), 0.10, 0.05));
	}
	addConstraint("NsignalSS1", NsignalSS1, npYield1,
		constraintSigmaFloor(npYield1, std::numeric_limits<double>::quiet_NaN(), 0.25, 1.0));
	if (useSignalSS2)
		addConstraint("NsignalSS2", NsignalSS2, npYield2,
			constraintSigmaFloor(npYield2, std::numeric_limits<double>::quiet_NaN(), 0.30, 1.0));
	if (useSignalSS3)
		addConstraint("NsignalSS3", NsignalSS3, npYield3,
			constraintSigmaFloor(npYield3, std::numeric_limits<double>::quiet_NaN(), 0.35, 1.0));
	for (auto &coeff : bkgCoeffVars)
	{
		const TString coeffName = coeff->GetName();
		const bool isFlipCoeff = coeffName.BeginsWith("Nflip");
		const bool isDsCoeff = coeffName.BeginsWith("Nds");
		const bool isPromptCoeff = coeffName == "Nplb";
		const double coeffSigma = constraintSigmaBounded(
			readResultValue(*bkgTimeResult, coeffName.Data(), coeff->getVal()),
			readResultError(*bkgTimeResult, coeffName.Data(), std::numeric_limits<double>::quiet_NaN()),
			isFlipCoeff ? 0.10 : (isDsCoeff ? 0.15 : 0.25),
			isFlipCoeff ? 5.0 : (isDsCoeff ? 3.0 : 10.0),
			isPromptCoeff ? 0.35 : (isFlipCoeff ? 0.50 : (isDsCoeff && isLowPtForwardBin ? 0.25 : 0.60)),
			isPromptCoeff ? 1500.0 : (isFlipCoeff ? 80.0 : (isDsCoeff && isLowPtForwardBin ? 60.0 : 500.0)));
		addConstraint(coeffName.Data(), *coeff,
			readResultValue(*bkgTimeResult, coeffName.Data(), coeff->getVal()),
			coeffSigma);
	}
	addConstraint("signal_lifetime", signal_lifetime, signalLifetime1Val,
		constraintSigmaFloor(signalLifetime1Val, std::numeric_limits<double>::quiet_NaN(), 0.15, 0.01));
	if (useSignalSS2)
		addConstraint("signal2_lifetime", signal2_lifetime, signalLifetime2Val,
			constraintSigmaFloor(signalLifetime2Val, std::numeric_limits<double>::quiet_NaN(), 0.15, 0.005));
	if (useSignalSS3)
		addConstraint("signal3_lifetime", signal3_lifetime, signalLifetime3Val,
			constraintSigmaFloor(signalLifetime3Val, std::numeric_limits<double>::quiet_NaN(), 0.15, 0.005));
	addConstraint("bkg_lifetime", bkg_lifetime,
		readBkgResultValue(*bkgTimeResult, "bkg_ss_lifetime", "bkg_lifetime", bkg_lifetime.getVal()),
		constraintSigmaFloor(bkg_lifetime.getVal(),
			dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime"))->getError() :
			dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime"))->getError() : 0.02,
			0.10, 0.01));
	if (nBkgSSComponents >= 2)
		addConstraint("bkg_lifetime2", bkg_lifetime2,
			readBkgResultValue(*bkgTimeResult, "bkg_ss_lifetime2", "bkg_lifetime2", bkg_lifetime2.getVal()),
			constraintSigmaFloor(bkg_lifetime2.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime2")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime2"))->getError() :
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime2")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime2"))->getError() : 0.02,
				0.10, 0.01));
	if (nBkgSSComponents >= 3)
		addConstraint("bkg_lifetime3", bkg_lifetime3,
			readBkgResultValue(*bkgTimeResult, "bkg_ss_lifetime3", "bkg_lifetime3", bkg_lifetime3.getVal()),
			constraintSigmaFloor(bkg_lifetime3.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime3")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ss_lifetime3"))->getError() :
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime3")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_lifetime3"))->getError() : 0.02,
				0.10, 0.015));
	if (nBkgDSComponents >= 1)
	{
		const double bkgSymLifetimeErr =
			dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime"))->getError() :
			dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime"))->getError() : 0.02;
		const double bkgSymLifetimeSigma = isLowPtForwardBin
			? constraintSigmaBounded(bkg_sym_lifetime.getVal(), bkgSymLifetimeErr, 0.10, 0.01, 0.50, 0.20)
			: constraintSigmaFloor(bkg_sym_lifetime.getVal(), bkgSymLifetimeErr, 0.10, 0.01);
		addConstraint("bkg_sym_lifetime", bkg_sym_lifetime,
			readBkgResultValue(*bkgTimeResult, "bkg_ds_lifetime", "bkg_sym_lifetime", bkg_sym_lifetime.getVal()),
			bkgSymLifetimeSigma);
	}
	if (nBkgDSComponents >= 2)
		addConstraint("bkg_sym_lifetime2", bkg_sym_lifetime2,
			readBkgResultValue(*bkgTimeResult, "bkg_ds_lifetime2", "bkg_sym_lifetime2", bkg_sym_lifetime2.getVal()),
			constraintSigmaFloor(bkg_sym_lifetime2.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime2")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime2"))->getError() :
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime2")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime2"))->getError() : 0.02,
				0.10, 0.015));
	if (nBkgDSComponents >= 3)
		addConstraint("bkg_sym_lifetime3", bkg_sym_lifetime3,
			readBkgResultValue(*bkgTimeResult, "bkg_ds_lifetime3", "bkg_sym_lifetime3", bkg_sym_lifetime3.getVal()),
			constraintSigmaFloor(bkg_sym_lifetime3.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime3")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_ds_lifetime3"))->getError() :
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime3")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_sym_lifetime3"))->getError() : 0.02,
				0.10, 0.02));
	if (nBkgFlipComponents >= 1)
		addConstraint("bkg_flip_lifetime", bkg_flip_lifetime,
			readResultValue(*bkgTimeResult, "bkg_flip_lifetime", bkg_flip_lifetime.getVal()),
			constraintSigmaBounded(bkg_flip_lifetime.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime"))->getError() : 0.02,
				0.15, 0.01, 0.75, 0.12));
	if (nBkgFlipComponents >= 2)
		addConstraint("bkg_flip_lifetime2", bkg_flip_lifetime2,
			readResultValue(*bkgTimeResult, "bkg_flip_lifetime2", bkg_flip_lifetime2.getVal()),
			constraintSigmaBounded(bkg_flip_lifetime2.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime2")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime2"))->getError() : 0.02,
				0.15, 0.015, 0.75, 0.20));
	if (nBkgFlipComponents >= 3)
		addConstraint("bkg_flip_lifetime3", bkg_flip_lifetime3,
			readResultValue(*bkgTimeResult, "bkg_flip_lifetime3", bkg_flip_lifetime3.getVal()),
			constraintSigmaBounded(bkg_flip_lifetime3.getVal(),
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime3")) ?
					dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_lifetime3"))->getError() : 0.02,
				0.15, 0.02, 0.75, 0.30));

	// ===== Plotting Only Model =====
	// only used to plot prompt and nonprompt components separately
	// These are not used for fitting.
	RooProdPdf prompt_core("prompt_core","prompt_core",RooArgSet(*signal_mass_pdf, *signalErrPdfPtr),Conditional(RooArgSet(*promptTimePdf),RooArgSet(obs_time)));
	RooProdPdf np_core("np_core","np_core",RooArgSet(*signal_mass_pdf, *signalErrPdfPtr),Conditional(RooArgSet(*signal_np_time),RooArgSet(obs_time)));
	// ===============================

	RooProdPdf signal_core("signal_core","signal_core",RooArgSet(*signal_mass_pdf, *signalErrPdfPtr),Conditional(RooArgSet(signal_time),RooArgSet(obs_time)));
	RooProdPdf bkg_core("bkg_core","bkg_core",RooArgSet(*bkg_mass_pdf, *bkgErrPdfPtr),Conditional(RooArgSet(*bkg_time_pdf),RooArgSet(obs_time)));
	RooAddPdf model("model","model",RooArgList(signal_core,bkg_core),RooArgList(Nsig,Nbkg));
	RooFitResult *model_result = nullptr;
	std::unique_ptr<RooFitResult> savedModelResult;
	if (drawFromSavedFit)
	{
		savedModelResult = clone_saved_fit_result(savedFitFile.get(), "modelResult");
		model_result = savedModelResult.get();
		if (!model_result)
		{
			std::cerr << "ERROR: modelResult not found in saved 2D fit file: " << fitRootName << std::endl;
			return;
		}
		apply_saved_fit_result(model_result, model, RooArgSet(obs_mass, obs_time, obs_timeErr));
	}
	else if (constraints.getSize() > 0)
		model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
			ExternalConstraints(constraints), RecoverFromUndefinedRegions(1.0));
	else
		model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
			RecoverFromUndefinedRegions(1.0));
	if (!drawFromSavedFit && model_result && (model_result->status() != 0 || model_result->covQual() < 2))
	{
		delete model_result;
		if (constraints.getSize() > 0)
			model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
				ExternalConstraints(constraints), RecoverFromUndefinedRegions(1.0));
		else
			model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight),
				RecoverFromUndefinedRegions(1.0));
		std::cout << "[OK] fit2d plain refit done" << std::endl;
	}

	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	TCanvas *cSigTimeErr = new TCanvas("sigTimeErrModel", "sigTimeErrModel", 800, 800);
	TPad *sigTimeErrPad1 = new TPad("sigTimeErrPad1", "sigTimeErrPad1", 0.0, 0.25, 1.0, 1.0);
	sigTimeErrPad1->SetBottomMargin(0.00001);
	sigTimeErrPad1->SetTopMargin(0.08);
	sigTimeErrPad1->SetLogy();
	sigTimeErrPad1->Draw();
	sigTimeErrPad1->cd();

	RooPlot *sigTimeErrPlot = obs_timeErr.frame(Range(errFitLow, errFitHigh), Title(""));
	if (isWeight)
		signalErrData->plotOn(sigTimeErrPlot, DataError(RooAbsData::SumW2), Name("data"));
	else
		signalErrData->plotOn(sigTimeErrPlot, Name("data"));
	signalErrPdfPtr->plotOn(sigTimeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (sigErrPdfOpt == kErrPdfAnalytic) {
		if (useErrSigLandau1) signalErrPdfPtr->plotOn(sigTimeErrPlot, Components(sigTimeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
		if (useErrSigLandau2) signalErrPdfPtr->plotOn(sigTimeErrPlot, Components(sigTimeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
		if (useErrSigLogn) signalErrPdfPtr->plotOn(sigTimeErrPlot, Components(sigTimeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
		if (useErrSigGaus1) signalErrPdfPtr->plotOn(sigTimeErrPlot, Components(sigTimeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
		if (useErrSigGaus2) signalErrPdfPtr->plotOn(sigTimeErrPlot, Components(sigTimeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
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
		sigTimeErrLeg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(sigTimeErrPlot, "tail"))
		sigTimeErrLeg.AddEntry(o, "Landau tail", "l");
	if (useErrSigLandau2)
		if (auto *o = findObj(sigTimeErrPlot, "tail2"))
		sigTimeErrLeg.AddEntry(o, "Landau 2", "l");
	if (auto *o = findObj(sigTimeErrPlot, "gaus1"))
		sigTimeErrLeg.AddEntry(o, "Gauss 1", "l");
	if (useErrSigGaus2)
		if (auto *o = findObj(sigTimeErrPlot, "gaus2"))
		sigTimeErrLeg.AddEntry(o, "Gauss 2", "l");
	if (useErrSigLogn)
		if (auto *o = findObj(sigTimeErrPlot, "logn"))
		sigTimeErrLeg.AddEntry(o, "Log-normal tail", "l");
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
		tx.DrawLatex(xtext, y0 + dy * k++, sigErrPdfOpt == kErrPdfHistPrMc ? "Prompt MC" : "SR - scaled SB");
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
		const int status = signalErrResult ? signalErrResult->status() : -1;
		int hesse = -1;
		if (signalErrResult)
		{
			for (UInt_t i = 0, n = signalErrResult->numStatusHistory(); i < n; ++i)
			{
				const char *lab = signalErrResult->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = signalErrResult->statusCodeHistory(i);
					break;
				}
			}
		}
		if (!publish)
			if (!publish)
				tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
	}

	cSigTimeErr->cd();
	TPad *sigTimeErrPad2 = new TPad("sigTimeErrPad2", "sigTimeErrPad2", 0.0, 0.0, 1.0, 0.25);
	sigTimeErrPad2->SetTopMargin(0.00001);
	sigTimeErrPad2->SetBottomMargin(0.4);
	sigTimeErrPad2->Draw();
	sigTimeErrPad2->cd();

	RooPlot *sigTimeErrPullPlot = obs_timeErr.frame(Range(errFitLow, errFitHigh), Title(""));
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

	auto sigTimeErrChi = chi2_from_pull(sigTimeErrPull);
	const int sigTimeErrNPar = signalErrResult ? signalErrResult->floatParsFinal().getSize() : 0;
	const int sigTimeErrNdf = std::max(1, sigTimeErrChi.second - sigTimeErrNPar);
	const double sigTimeErrPvalue = TMath::Prob(sigTimeErrChi.first, sigTimeErrNdf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", sigTimeErrChi.first, sigTimeErrNdf));
	}

	TLine sigTimeErrLine(errFitLow, 0.0, errFitHigh, 0.0);
	sigTimeErrLine.SetLineStyle(2);
	sigTimeErrLine.Draw("same");
	cSigTimeErr->Print(figName("errSig"));
	delete cSigTimeErr;

	TCanvas *cTimeErr = new TCanvas("timeErrModel", "timeErrModel", 800, 800);
	TPad *timeErrPad1 = new TPad("timeErrPad1", "timeErrPad1", 0.0, 0.25, 1.0, 1.0);
	timeErrPad1->SetBottomMargin(0.00001);
	timeErrPad1->SetTopMargin(0.08);
	timeErrPad1->SetLogy();
	timeErrPad1->Draw();
	timeErrPad1->cd();

	RooPlot *timeErrPlot = obs_timeErr.frame(Range(errFitLow, errFitHigh), Title(""));
	if (isWeight)
		bkgErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errFitLow, errFitHigh), DataError(RooAbsData::SumW2), Name("data"));
	else
		bkgErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins, errFitLow, errFitHigh), Name("data"));
	bkgErrPdfPtr->plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (bkgErrPdfOpt == kErrPdfAnalytic) {
		if (useErrBkgGaus1) bkgErrPdfPtr->plotOn(timeErrPlot, Components(bkgTimeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
		if (useErrBkgGaus2) bkgErrPdfPtr->plotOn(timeErrPlot, Components(bkgTimeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
		if (useErrBkgLandau1) bkgErrPdfPtr->plotOn(timeErrPlot, Components(bkgTimeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
		if (useErrBkgLandau2) bkgErrPdfPtr->plotOn(timeErrPlot, Components(bkgTimeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
		if (useErrBkgLogn) bkgErrPdfPtr->plotOn(timeErrPlot, Components(bkgTimeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
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
		timeErrLeg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(timeErrPlot, "tail"))
		timeErrLeg.AddEntry(o, "Landau tail", "l");
	if (useErrBkgLandau2)
		if (auto *o = findObj(timeErrPlot, "tail2"))
		timeErrLeg.AddEntry(o, "Landau 2", "l");
	if (auto *o = findObj(timeErrPlot, "gaus1"))
		timeErrLeg.AddEntry(o, "Gauss 1", "l");
	if (useErrBkgGaus2)
		if (auto *o = findObj(timeErrPlot, "gaus2"))
		timeErrLeg.AddEntry(o, "Gauss 2", "l");
	if (useErrBkgLogn)
		if (auto *o = findObj(timeErrPlot, "logn"))
		timeErrLeg.AddEntry(o, "Log-normal tail", "l");
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
		const int status = bkgErrResult ? bkgErrResult->status() : -1;
		int hesse = -1;
		if (bkgErrResult)
		{
			for (UInt_t i = 0, n = bkgErrResult->numStatusHistory(); i < n; ++i)
			{
				const char *lab = bkgErrResult->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = bkgErrResult->statusCodeHistory(i);
					break;
				}
			}
		}
		if (!publish)
			if (!publish)
				tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
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
	const int timeErrNPar = bkgErrResult ? bkgErrResult->floatParsFinal().getSize() : 0;
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
	cTimeErr->Print(figName("errBkg"));
	delete cTimeErr;

	TCanvas *cMass = new TCanvas("cMass", "cMass", 800, 800);
	TPad *massPad1 = new TPad("massPad1", "massPad1", 0.0, 0.25, 1.0, 1.0);
	massPad1->SetBottomMargin(0.00001);
	massPad1->SetTopMargin(0.08);
	massPad1->Draw();
	massPad1->cd();

	RooPlot *massfitplot = obs_mass.frame(Title(""));
	if (isWeight)
		data->plotOn(massfitplot, DataError(RooAbsData::SumW2), Name("data"));
	else
		data->plotOn(massfitplot, Name("data"));
	model.plotOn(massfitplot, LineColor(kBlack), LineWidth(2), Name("model"));
	model.plotOn(massfitplot, Components(signal_core), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Name("signal_component"));
	model.plotOn(massfitplot, Components(bkg_core), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("bkg_component"));
	double massYmax = -1e300;
	if (auto *hdata = dynamic_cast<RooHist *>(massfitplot->getHist("data")))
	{
		for (int i = 0; i < hdata->GetN(); ++i)
		{
			double xval = 0.0, yval = 0.0;
			hdata->GetPoint(i, xval, yval);
			if (yval > massYmax)
				massYmax = yval;
		}
	}
	if (massYmax > 0.0 && massYmax < 1e300)
		massfitplot->SetMaximum(massYmax * 1.8);
	massfitplot->GetYaxis()->SetTitle("Events");
	massfitplot->GetYaxis()->SetTitleOffset(1.6);
	massfitplot->GetXaxis()->SetTitle("");
	massfitplot->Draw("e");

	TLegend massLeg(0.50, 0.72, 0.74, 0.89);
	massLeg.SetBorderSize(0);
	massLeg.SetFillStyle(0);
	massLeg.SetTextSize(0.03);
	if (auto *o = findObj(massfitplot, "data"))
		massLeg.AddEntry(o, "Data", "lep");
	if (auto *o = findObj(massfitplot, "model"))
		massLeg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(massfitplot, "signal_component"))
		massLeg.AddEntry(o, "Signal", "l");
	if (auto *o = findObj(massfitplot, "bkg_component"))
		massLeg.AddEntry(o, "Background", "l");
	massLeg.Draw("same");

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
		tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
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
		const int status = model_result ? model_result->status() : -1;
		int hesse = -1;
		if (model_result)
		{
			for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = model_result->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = model_result->statusCodeHistory(i);
					break;
				}
			}
		}
		if (!publish)
			tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.74, y0 = 0.87, dy = -0.045;
		int k = 0;
		auto printReal = [&](const char *title, RooAbsReal &var)
		{
			if (auto *rrv = dynamic_cast<RooRealVar *>(&var))
			{
				if (rrv->isConstant())
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
				else
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
				return;
			}
			if (model_result)
			{
				const double err = var.getPropagatedError(*model_result);
				if (err > 0.0 && std::isfinite(err))
				{
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
					return;
				}
			}
			tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
		};
		auto printValue = [&](const char *title, double val, double err = -1.0, bool isFixed = false)
		{
			if (isFixed)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, val));
			else if (err >= 0.0)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, val, err));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, val));
		};
		if (publish)
		{
			printReal("f_{B}", bFraction);
			printReal("N_{sig}", Nsig);
			printReal("N_{bkg}", Nbkg);
		}
		else
		{
			printReal("N_{sig}", Nsig);
			printReal("N_{bkg}", Nbkg);
			printReal("f_{B}", bFraction);
			printReal("m_{0}", signal_mass_mean);
			if (nSignalGaussComponents >= 1)
				printReal("#sigma_{G1}", signal_mass_sigma);
			if (nSignalGaussComponents >= 2)
				printReal("#sigma_{G2}", signal_mass_sigma2);
			if (nSignalCBComponents >= 1)
			{
				printValue("#sigma_{CB1}", signal_mass_cb_sigma.getVal(), -1.0, true);
				printValue("#alpha_{CB1}", signal_mass_cb_alpha.getVal(), -1.0, true);
				printValue("n_{CB1}", signal_mass_cb_n.getVal(), -1.0, true);
			}
			if (nSignalCBComponents >= 2)
			{
				printValue("#alpha_{CB2}", signal_mass_cb_alpha2.getVal(), -1.0, true);
				printValue("n_{CB2}", signal_mass_cb_n2.getVal(), -1.0, true);
				printValue("#sigma_{CB2}", signal_mass_cb_sigma2.getVal(), -1.0, true);
			}
			if (signalMassPdfList.getSize() >= 2)
				printValue("f_{sig1}", signal_mass_frac1.getVal(), -1.0, true);
			if (signalMassPdfList.getSize() >= 3)
				printValue("f_{sig2}", signal_mass_frac2.getVal(), -1.0, true);
			if (signalMassPdfList.getSize() >= 4)
				printValue("f_{sig3}", signal_mass_frac3.getVal(), -1.0, true);
			if (nBkgExpComponents == 1)
				printValue("#lambda_{bkg}", bkg_mass_lambda.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 1)
				printValue("p_{1}", bkg_mass_p1.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 2)
				printValue("p_{2}", bkg_mass_p2.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 3)
				printValue("p_{3}", bkg_mass_p3.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 4)
				printValue("p_{4}", bkg_mass_p4.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 5)
				printValue("p_{5}", bkg_mass_p5.getVal(), -1.0, true);
			if (nBkgChebyOrder >= 6)
				printValue("p_{6}", bkg_mass_p6.getVal(), -1.0, true);
		}
	}

	cMass->cd();
	TPad *massPad2 = new TPad("massPad2", "massPad2", 0.0, 0.0, 1.0, 0.25);
	massPad2->SetTopMargin(0.00001);
	massPad2->SetBottomMargin(0.4);
	massPad2->Draw();
	massPad2->cd();

	RooPlot *massPullPlot = obs_mass.frame(Title(""));
	RooHist *massPull = massfitplot->pullHist("data", "model");
	if (massPull)
		massPullPlot->addPlotable(massPull, "P");
	massPullPlot->GetYaxis()->SetTitle("Pull");
	massPullPlot->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
	massPullPlot->GetXaxis()->CenterTitle();
	massPullPlot->SetMinimum(-8);
	massPullPlot->SetMaximum(8);
	massPullPlot->GetYaxis()->SetNdivisions(505);
	massPullPlot->GetYaxis()->SetTitleSize(0.12);
	massPullPlot->GetYaxis()->SetLabelSize(0.10);
	massPullPlot->GetXaxis()->SetTitleSize(0.15);
	massPullPlot->GetXaxis()->SetLabelSize(0.10);
	massPullPlot->Draw();

	auto massChi = chi2_from_pull(massPull);
	int massNpar = model_result ? model_result->floatParsFinal().getSize() : 0;
	int massNdf = std::max(1, massChi.second - massNpar);
	double massPvalue = TMath::Prob(massChi.first, massNdf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", massChi.first, massNdf));
	}
	TLine massLine(obs_mass.getMin(), 0.0, obs_mass.getMax(), 0.0);
	massLine.SetLineStyle(2);
	massLine.Draw("same");
	cMass->Print(figName("mass_fit"));
	delete cMass;

	TCanvas *cLifetime = new TCanvas("cLifetime", "cLifetime", 800, 800);
	TPad *timePad1 = new TPad("timePad1", "timePad1", 0.0, 0.25, 1.0, 1.0);
	timePad1->SetBottomMargin(0.00001);
	timePad1->SetTopMargin(0.08);
	timePad1->SetLogy();
	timePad1->Draw();
	timePad1->cd();

	RooPlot *timefitplot = obs_time.frame(Range(ctRange.first, ctRange.second), Title(""));
	if (isWeight)
		data->plotOn(timefitplot, Binning(timePlotBins), DataError(RooAbsData::SumW2), Name("data"));
	else
		data->plotOn(timefitplot, Binning(timePlotBins), Name("data"));
	model.plotOn(timefitplot, LineColor(kBlack), LineWidth(2), Name("model"));
	const double npYield = Nsig.getVal() * bFraction.getVal();
	const double prYield = std::max(0.0, Nsig.getVal() - npYield);
	prompt_core.plotOn(timefitplot, LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2),
		Normalization(prYield, RooAbsReal::NumEvent), Name("pr_component"));
	np_core.plotOn(timefitplot, LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2),
		Normalization(npYield, RooAbsReal::NumEvent), Name("np_component"));
	bkg_core.plotOn(timefitplot, LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2),
		Normalization(Nbkg.getVal(), RooAbsReal::NumEvent), Name("bkg_component"));
	apply_logy_auto_range(timefitplot, "data");
	timefitplot->GetYaxis()->SetTitle("Events");
	timefitplot->GetYaxis()->SetTitleOffset(1.6);
	timefitplot->GetXaxis()->SetTitle("");
	timefitplot->Draw("e");

	TLegend timeLeg(0.50, 0.72, 0.74, 0.89);
	timeLeg.SetBorderSize(0);
	timeLeg.SetFillStyle(0);
	timeLeg.SetTextSize(0.03);
	if (auto *o = findObj(timefitplot, "data"))
		timeLeg.AddEntry(o, "Data", "lep");
	if (auto *o = findObj(timefitplot, "model"))
		timeLeg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(timefitplot, "pr_component"))
		timeLeg.AddEntry(o, "PR", "l");
	if (auto *o = findObj(timefitplot, "np_component"))
		timeLeg.AddEntry(o, "NP", "l");
	if (auto *o = findObj(timefitplot, "bkg_component"))
		timeLeg.AddEntry(o, "Bkg", "l");
	timeLeg.Draw("same");

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
		tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
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
		const int status = model_result ? model_result->status() : -1;
		int hesse = -1;
		if (model_result)
		{
			for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = model_result->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = model_result->statusCodeHistory(i);
					break;
				}
			}
		}
		if (!publish)
			tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.019);
		tp.SetTextFont(42);
		double xtext = 0.71, y0 = 0.88, dy = -0.034;
		int k = 0;
		auto printReal = [&](const char *title, RooAbsReal &var)
		{
			if (auto *rrv = dynamic_cast<RooRealVar *>(&var))
			{
				if (rrv->isConstant())
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
				else
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
				return;
			}
			if (model_result)
			{
				const double err = var.getPropagatedError(*model_result);
				if (err > 0.0 && std::isfinite(err))
				{
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
					return;
				}
			}
			tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
		};
		auto printValue = [&](const char *title, double val, double err = -1.0, bool isFixed = false)
		{
			if (isFixed)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, val));
			else if (err >= 0.0)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, val, err));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, val));
		};
		if (publish)
		{
			printReal("f_{B}", bFraction);
			printReal("N_{sig}", Nsig);
			printReal("N_{bkg}", Nbkg);
		}
		else
		{
			printReal("f_{B}", bFraction);
			printReal("N_{sig}", Nsig);
			printReal("N_{bkg}", Nbkg);
			printReal("ctau s_{1}", ctauTime1Scale);
			if (nResolutionComponents >= 2)
				printReal("#Deltactau s_{2}", ctauTime2Delta);
			if (nResolutionComponents >= 3)
				printReal("#Deltactau s_{3}", ctauTime3Delta);
			if (nResolutionComponents >= 4)
				printReal("#Deltactau s_{4}", ctauTime4Delta);
			printReal("#tau_{NP1}", signal_lifetime);
			if (useSignalSS2)
				printReal("#tau_{NP2}", signal2_lifetime);
			if (useSignalSS3)
				printReal("#tau_{NP3}", signal3_lifetime);
			printReal("#tau_{Bkg}", bkg_lifetime);
			if (nBkgSSComponents >= 2)
				printReal("#tau_{Bkg2}", bkg_lifetime2);
			if (nBkgSSComponents >= 3)
				printReal("#tau_{Bkg3}", bkg_lifetime3);
			if (nBkgFlipComponents >= 1)
				printReal("#tau_{flip}", bkg_flip_lifetime);
			if (nBkgFlipComponents >= 2)
				printReal("#tau_{flip2}", bkg_flip_lifetime2);
			if (nBkgFlipComponents >= 3)
				printReal("#tau_{flip3}", bkg_flip_lifetime3);
			printValue("#tau_{sym}", bkg_sym_lifetime.getVal(), -1.0, bkg_sym_lifetime.isConstant());
		}
	}

	cLifetime->cd();
	TPad *timePad2 = new TPad("timePad2", "timePad2", 0.0, 0.0, 1.0, 0.25);
	timePad2->SetTopMargin(0.00001);
	timePad2->SetBottomMargin(0.4);
	timePad2->Draw();
	timePad2->cd();

	RooPlot *timePullPlot = obs_time.frame(Range(ctRange.first, ctRange.second), Title(""));
	RooHist *timePull = timefitplot->pullHist("data", "model");
	if (timePull)
		timePullPlot->addPlotable(timePull, "P");
	timePullPlot->GetYaxis()->SetTitle("Pull");
	timePullPlot->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} [mm]");
	timePullPlot->GetXaxis()->CenterTitle();
	timePullPlot->SetMinimum(-8);
	timePullPlot->SetMaximum(8);
	timePullPlot->GetYaxis()->SetNdivisions(505);
	timePullPlot->GetYaxis()->SetTitleSize(0.12);
	timePullPlot->GetYaxis()->SetLabelSize(0.10);
	timePullPlot->GetXaxis()->SetTitleSize(0.15);
	timePullPlot->GetXaxis()->SetLabelSize(0.10);
	timePullPlot->Draw();

	auto timeChi = chi2_from_pull(timePull);
	int timeNpar = model_result ? model_result->floatParsFinal().getSize() : 0;
	int timeNdf = std::max(1, timeChi.second - timeNpar);
	double timePvalue = TMath::Prob(timeChi.first, timeNdf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", timeChi.first, timeNdf));
	}
	TLine timeLine(ctRange.first, 0.0, ctRange.second, 0.0);
	timeLine.SetLineStyle(2);
	timeLine.Draw("same");
	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

	if (!drawFromSavedFit)
	{
		std::unique_ptr<TFile> fitOut(TFile::Open(fitRootName, "RECREATE"));
		if (!fitOut || fitOut->IsZombie())
		{
			std::cerr << "ERROR: cannot create fit2d output file: " << fitRootName << std::endl;
		}
		else
		{
			int hesse = -1;
			if (model_result)
			{
				for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
				{
					const char *lab = model_result->statusLabelHistory(i);
					if (lab && TString(lab) == "HESSE")
					{
						hesse = model_result->statusCodeHistory(i);
						break;
					}
				}
				model_result->Write("modelResult");
			}
			TParameter<double>("ptLow", ptLow).Write();
			TParameter<double>("ptHigh", ptHigh).Write();
			TParameter<double>("yLow", yLow).Write();
			TParameter<double>("yHigh", yHigh).Write();
			TParameter<double>("bFraction", bFraction.getVal()).Write();
			TParameter<double>("bFractionErr", bFraction.getError()).Write();
			TParameter<double>("Nsig", Nsig.getVal()).Write();
			TParameter<double>("Nbkg", Nbkg.getVal()).Write();
			TParameter<int>("fixResolutionParams", fixResolutionParams ? 1 : 0).Write();
			TParameter<double>("maxPromptScaleSaved", maxPromptScaleSaved).Write();
			TParameter<double>("errFitHigh", errFitHigh).Write();
			TParameter<double>("lifetimeFloorCoeff", lifetimeFloorCoeff).Write();
			TParameter<double>("absoluteLifetimeFloor", absoluteLifetimeFloor).Write();
			TParameter<double>("resolutionDrivenLifetimeFloor", resolutionDrivenLifetimeFloor).Write();
			TParameter<double>("signalLifetimeFloorApplied", signalLifetimeFloor).Write();
			TParameter<double>("bkgLifetimeFloorApplied", bkgLifetimeFloor).Write();
			TParameter<double>("sigTimeErrChi2", sigTimeErrChi.first).Write();
			TParameter<int>("sigTimeErrNdf", sigTimeErrNdf).Write();
			TParameter<double>("sigTimeErrPvalue", sigTimeErrPvalue).Write();
			TParameter<double>("timeErrChi2", timeErrChi.first).Write();
			TParameter<int>("timeErrNdf", timeErrNdf).Write();
			TParameter<double>("timeErrPvalue", timeErrPvalue).Write();
			TParameter<double>("massChi2", massChi.first).Write();
			TParameter<int>("massNdf", massNdf).Write();
			TParameter<double>("massPvalue", massPvalue).Write();
			TParameter<double>("timeChi2", timeChi.first).Write();
			TParameter<int>("timeNdf", timeNdf).Write();
			TParameter<double>("timePvalue", timePvalue).Write();
			TParameter<int>("fitStatus", model_result ? model_result->status() : -1).Write();
			TParameter<int>("hesseStatus", hesse).Write();
			fitOut->Write();
		}
	}
	else
		std::cout << "[PlotOnly] Loaded saved 2D fit and left ROOT file unchanged: " << fitRootName << std::endl;

	const TString figErrSig = figName("errSig");
	const TString figErrBkg = figName("errBkg");
	const TString figMass = figName("mass_fit");
	const TString figLifetime = figName("lifetime_fit");
	cout << "----------------- FIT RESULT FOR THE 2D MODEL ------------" << endl;
	cout << "bkgErrPdfOpt in fit2d: " << bkgErrPdfOpt
			 << " (" << errPdfChoiceLabel(bkgErrPdfOpt) << ")" << endl;
	cout << "sigErrPdfOpt in fit2d: " << sigErrPdfOpt
			 << " (" << errPdfChoiceLabel(sigErrPdfOpt) << ")" << endl;
	cout << "histPdfInterpolationOrder in fit2d: " << histPdfInterpolationOrder << endl;
	if (bkgErrResult)
	{
		cout << "---------------- FIT RESULT FOR BKG ERR PDF --------------" << endl;
		bkgErrResult->Print("v");
	}
	if (signalErrResult)
	{
		cout << "---------------- FIT RESULT FOR SIG ERR PDF --------------" << endl;
		signalErrResult->Print("v");
	}
	model_result->Print("v");
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
			else if (model_result)
			{
				err = var.getPropagatedError(*model_result);
				hasErr = std::isfinite(err) && err > 0.0;
			}
			const double value2Sigma = var.getVal() + (hasErr ? 2.0 * err : 0.0);
			std::cout << Form("[FLOOR] %s : floor=%.6g , value+2sigma=%.6g", name, floor, value2Sigma);
			if (std::isfinite(floor))
				std::cout << Form(" , margin=%.6g", value2Sigma - floor);
			std::cout << std::endl;
		};
		std::cout << "---------------- FLOOR VS VALUE+2SIGMA -------------" << std::endl;
		printFloorVsValue2Sigma("signal_lifetime", signal_lifetime, signalLifetimeFloorSaved);
		if (useSignalSS2)
			printFloorVsValue2Sigma("signal2_lifetime", signal2_lifetime, signalLifetime2Floor);
		if (useSignalSS3)
			printFloorVsValue2Sigma("signal3_lifetime", signal3_lifetime, signalLifetime3Floor);
		printFloorVsValue2Sigma("bkg_lifetime", bkg_lifetime, bkgLifetimeFloor);
		if (nBkgSSComponents >= 2)
			printFloorVsValue2Sigma("bkg_lifetime2", bkg_lifetime2, bkgLifetime2Floor);
		if (nBkgSSComponents >= 3)
			printFloorVsValue2Sigma("bkg_lifetime3", bkg_lifetime3, bkgLifetime3Floor);
		if (nBkgDSComponents >= 1)
			printFloorVsValue2Sigma("bkg_sym_lifetime", bkg_sym_lifetime, bkgSymLifetimeFloor);
		if (nBkgDSComponents >= 2)
			printFloorVsValue2Sigma("bkg_sym_lifetime2", bkg_sym_lifetime2, bkgSymLifetime2Floor);
		if (nBkgDSComponents >= 3)
			printFloorVsValue2Sigma("bkg_sym_lifetime3", bkg_sym_lifetime3, bkgSymLifetime3Floor);
		if (nBkgFlipComponents >= 1)
			printFloorVsValue2Sigma("bkg_flip_lifetime", bkg_flip_lifetime, bkgFlipLifetimeFloor);
		if (nBkgFlipComponents >= 2)
			printFloorVsValue2Sigma("bkg_flip_lifetime2", bkg_flip_lifetime2, bkgFlipLifetime2Floor);
		if (nBkgFlipComponents >= 3)
			printFloorVsValue2Sigma("bkg_flip_lifetime3", bkg_flip_lifetime3, bkgFlipLifetime3Floor);
	}
	std::cout << "[FIG] fit2d signal err fit : " << figErrSig << std::endl;
	std::cout << "[FIG] fit2d background err fit : " << figErrBkg << std::endl;
	std::cout << "[FIG] fit2d mass fit : " << figMass << std::endl;
	std::cout << "[FIG] fit2d lifetime fit : " << figLifetime << std::endl;
}
