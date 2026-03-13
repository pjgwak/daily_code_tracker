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
#include "TString.h"
#include "TParameter.h"
#include "TMath.h"
#include "RooHist.h"
#include "RooLognormal.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace RooFit;

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

void fit2d(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4){
	bool isWeight = false;
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	enum ErrPdfChoice
	{
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
		kErrPdfHistPrMc = 2,
	};
	const int overrideBkgErrPdfOpt = -1; // -1: automatic choice, default
	const int overrideSigErrPdfOpt = 2; // set to kErrPdfHistPrMc to build sig err RooHistPdf directly from PR_MC_ROOT
	const int histPdfInterpolationOrder = 1;
	auto isHistErrPdfChoice = [](int choice) {
		return choice == kErrPdfHist || choice == kErrPdfHistPrMc;
	};
	auto errPdfChoiceLabel = [&](int choice) -> const char * {
		switch (choice) {
		case kErrPdfHist:
			return "RooHistPdf";
		case kErrPdfHistPrMc:
			return "RooHistPdf(PR_MC_ROOT)";
		default:
			return "analytic";
		}
	};

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
	const char *PR_MC_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	const auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.001, 0.999, false);
	auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.995, true);
	if (errRange.first < 1e-6) errRange.first = 1e-6;
	const double errMin = 0.0;
	auto formatTag = [](double value) { return TString::Format("%.2f", value); };
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString figDir = TString::Format("figs/%s/fit2d", yTag.Data());
	const TString fitRootDir = TString::Format("roots/%s/fit2d", yTag.Data());
	const TString fitRootName = TString::Format("%s/fit2d_result_%s.root", fitRootDir.Data(), figTag.Data());
	const TString errFileName = TString::Format("roots/%s/err2/err2_model_%s.root", yTag.Data(), figTag.Data());

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

	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	RooRealVar &obs_timeErr = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_time.SetTitle("ctau3D");
	obs_time.setUnit("mm");
	obs_timeErr.SetTitle("event-by-event ctau error");
	obs_timeErr.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	obs_timeErr.setRange(errFitLow, errFitHigh);
	obs_timeErr.setMin(errFitLow);
	const int timePlotBins = std::max(2, obs_time.getBins());

	// Keep the decay constants away from the numerically unstable tau -> 0 limit.
	const double signalLifetimeFloor = std::max(1.0 * errRange.first, 1e-2);
	const double bkgLifetimeFloor = std::max(1.0 * errRange.first, 1e-2);
	const double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	RooDataSet *data = dataSel.get();

	const double sidebandLeftMax = 2.9;
	const double sidebandRightMin = 3.2;
	TString cutSignalRegion = Form("(mass >= %g && mass <= %g)", sidebandLeftMax, sidebandRightMin);
	TString cutSideband = Form("(mass >= 2.6 && mass < %g) || (mass > %g && mass <= 3.5)", sidebandLeftMax, sidebandRightMin);
	auto dataSB = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(data->reduce(cutSideband)));
	if (!dataSB || dataSB->numEntries() <= 0) dataSB.reset(static_cast<RooDataSet *>(data->reduce(RooArgSet(obs_mass, obs_time, obs_timeErr))));

	const int timeErrPlotBinsDefault = std::max(2, obs_timeErr.getBins());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(fitRootDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	const TString massFileName = TString::Format("roots/%s/mass/mass_model_%s.root", yTag.Data(), figTag.Data());
	const TString prFileName = TString::Format("roots/%s/ctau_pr/ctau_resolution_%s.root", yTag.Data(), figTag.Data());
	const TString npFileName = TString::Format("roots/%s/ctau_np/ctau_np_model_%s.root", yTag.Data(), figTag.Data());
	const TString bkgFileName = TString::Format("roots/%s/ctau_bkg/ctau_bkg_fitresult_%s.root", yTag.Data(), figTag.Data());

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

	std::unique_ptr<TFile> massFile(openFile(massFileName));
	std::unique_ptr<TFile> prFile(openFile(prFileName));
		std::unique_ptr<TFile> npFile(openFile(npFileName));
			std::unique_ptr<TFile> bkgFile(openFile(bkgFileName));
			std::unique_ptr<TFile> errFile(openFile(errFileName));
			if (!massFile || !prFile || !npFile || !bkgFile || !errFile) return;
			const int timeErrPlotBins = std::max(2, readIntParam(*errFile, "timeErrPlotBins", timeErrPlotBinsDefault));
			const int savedBkgErrPdfOpt = readIntParam(*errFile, "bkgErrPdfOpt", kErrPdfAnalytic);
			const int savedSigErrPdfOpt = readIntParam(*errFile, "sigErrPdfOpt", kErrPdfAnalytic);
			const int bkgErrPdfOpt = overrideBkgErrPdfOpt >= 0 ? overrideBkgErrPdfOpt : savedBkgErrPdfOpt;
			const int sigErrPdfOpt = overrideSigErrPdfOpt >= 0 ? overrideSigErrPdfOpt : savedSigErrPdfOpt;
			const int nErrBkgGauss = std::clamp(readIntParam(*errFile, "nBkgTimeErrGaussComponents", 2), 0, 2);
			const int nErrBkgLandau = std::clamp(readIntParam(*errFile, "nBkgTimeErrLandauComponents", 2), 0, 2);
			const int nErrBkgLogn = std::clamp(readIntParam(*errFile, "nBkgTimeErrLognormalComponents", 1), 0, 1);
			const int nErrSigGauss = std::clamp(readIntParam(*errFile, "nSigTimeErrGaussComponents", 2), 0, 2);
			const int nErrSigLandau = std::clamp(readIntParam(*errFile, "nSigTimeErrLandauComponents", 2), 0, 2);
			const int nErrSigLogn = std::clamp(readIntParam(*errFile, "nSigTimeErrLognormalComponents", 1), 0, 1);
			const double savedScaleBkg = readDoubleParam(*errFile, "scaleBkg", std::numeric_limits<double>::quiet_NaN());

		RooFitResult *bkgErrResult = dynamic_cast<RooFitResult *>(errFile->Get("timeErrResult"));
		RooFitResult *signalErrResult = dynamic_cast<RooFitResult *>(errFile->Get("sigTimeErrResult"));
		if ((bkgErrPdfOpt == kErrPdfAnalytic && !bkgErrResult) ||
		    (sigErrPdfOpt == kErrPdfAnalytic && !signalErrResult)) {
			std::cerr << "ERROR: required analytic err fit results not found in " << errFileName << std::endl;
			return;
		}
	if (!dataSB || dataSB->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after sideband selection: " << cutSideband << std::endl;
		return;
	}
		auto readErrVar = [&](RooFitResult *fr, const char *name, double fallback = 0.0) {
			return fr ? readResultValue(*fr, name, fallback) : fallback;
		};
		obs_timeErr.setRange(errFitLow, errFitHigh);
		obs_timeErr.setMin(errFitLow);
		const bool useErrBkgGaus1 = (nErrBkgGauss >= 1);
		const bool useErrBkgGaus2 = (nErrBkgGauss >= 2);
		const bool useErrBkgLandau1 = (nErrBkgLandau >= 1);
		const bool useErrBkgLandau2 = (nErrBkgLandau >= 2);
		const bool useErrBkgLogn = (nErrBkgLogn >= 1);
		const int nErrBkgComponents =
			static_cast<int>(useErrBkgGaus1) + static_cast<int>(useErrBkgGaus2) +
			static_cast<int>(useErrBkgLandau1) + static_cast<int>(useErrBkgLandau2) +
			static_cast<int>(useErrBkgLogn);
		const bool useErrSigGaus1 = (nErrSigGauss >= 1);
		const bool useErrSigGaus2 = (nErrSigGauss >= 2);
		const bool useErrSigLandau1 = (nErrSigLandau >= 1);
		const bool useErrSigLandau2 = (nErrSigLandau >= 2);
		const bool useErrSigLogn = (nErrSigLogn >= 1);
		const int nErrSigComponents =
			static_cast<int>(useErrSigGaus1) + static_cast<int>(useErrSigGaus2) +
			static_cast<int>(useErrSigLandau1) + static_cast<int>(useErrSigLandau2) +
			static_cast<int>(useErrSigLogn);
		auto fracToRatio = [](double frac) {
			const double f = std::clamp(frac, 1e-6, 1.0 - 1e-6);
			return f / (1.0 - f);
		};
		if ((bkgErrPdfOpt == kErrPdfAnalytic && nErrBkgComponents <= 0) ||
		    (sigErrPdfOpt == kErrPdfAnalytic && nErrSigComponents <= 0)) {
			std::cerr << "ERROR: invalid saved err2 component configuration in " << errFileName << std::endl;
			return;
		}

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
		RooRealVar bkgTimeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", readErrVar(bkgErrResult, "timeErrGaus1Mean", 0.02), errFitLow, errUpper);
		RooRealVar bkgTimeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", readErrVar(bkgErrResult, "timeErrGaus1Sigma", 0.005), 1e-6, errUpper);
		RooRealVar bkgTimeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", readErrVar(bkgErrResult, "timeErrGaus2Mean", bkgTimeErrGaus1Mean.getVal()), errFitLow, errUpper);
		RooRealVar bkgTimeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", readErrVar(bkgErrResult, "timeErrGaus2Sigma", bkgTimeErrGaus1Sigma.getVal()), 1e-6, errUpper);
			RooRealVar bkgTimeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", readErrVar(bkgErrResult, "timeErrTailMpv", 0.02), errFitLow, errUpper);
			RooRealVar bkgTimeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", readErrVar(bkgErrResult, "timeErrTailWidth", 0.002), 1e-6, errUpper);
			RooRealVar bkgTimeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", readErrVar(bkgErrResult, "timeErrTail2Mpv", bkgTimeErrTailMpv.getVal()), errFitLow, errUpper);
			RooRealVar bkgTimeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", readErrVar(bkgErrResult, "timeErrTail2Width", bkgTimeErrTailWidth.getVal()), 1e-6, errUpper);
			RooRealVar bkgTimeErrLognM0("timeErrLognM0", "timeErrLognM0", readErrVar(bkgErrResult, "timeErrLognM0", 0.02), errFitLow, errUpper);
			RooRealVar bkgTimeErrLognK("timeErrLognK", "timeErrLognK", readErrVar(bkgErrResult, "timeErrLognK", 0.35), 0.05, 0.95);
			RooRealVar bkgTimeErrCore1FracRatio("timeErrCore1FracRatio", "timeErrCore1FracRatio", readErrVar(bkgErrResult, "timeErrCore1FracRatio", fracToRatio(0.6)), 1e-6, 1e6);
			RooFormulaVar bkgTimeErrCore1Frac("timeErrCore1Frac", "@0/(1.0+@0)", RooArgList(bkgTimeErrCore1FracRatio));
			RooRealVar bkgTimeErrTailFracRatio("timeErrTailFracRatio", "timeErrTailFracRatio", readErrVar(bkgErrResult, "timeErrTailFracRatio", fracToRatio(0.08)), 1e-6, 1e6);
			RooFormulaVar bkgTimeErrTailFrac("timeErrTailFrac", "@0/(1.0+@0)", RooArgList(bkgTimeErrTailFracRatio));
			RooRealVar bkgTimeErrTail2FracRatio("timeErrTail2FracRatio", "timeErrTail2FracRatio", readErrVar(bkgErrResult, "timeErrTail2FracRatio", fracToRatio(0.05)), 1e-6, 1e6);
			RooFormulaVar bkgTimeErrTail2Frac("timeErrTail2Frac", "@0/(1.0+@0)", RooArgList(bkgTimeErrTail2FracRatio));
			RooRealVar bkgTimeErrLognFracRatio("timeErrLognFracRatio", "timeErrLognFracRatio", readErrVar(bkgErrResult, "timeErrLognFracRatio", fracToRatio(0.03)), 1e-6, 1e6);
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
			if (bkgErrPdfList.getSize() < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrTail2Frac);
		}
		if (useErrBkgLogn) {
			bkgErrPdfList.add(bkgTimeErrLogn);
			if (bkgErrPdfList.getSize() < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrLognFrac);
		}
		if (useErrBkgGaus1) {
			bkgErrPdfList.add(bkgTimeErrGaus1);
			if (bkgErrPdfList.getSize() < nErrBkgComponents) bkgErrFracList.add(bkgTimeErrCore1Frac);
		}
		if (useErrBkgGaus2)
			bkgErrPdfList.add(bkgTimeErrGaus2);
		std::unique_ptr<RooAddPdf> bkgErrPdfAnalytic = std::make_unique<RooAddPdf>(
			"pdfErr", "pdfErr", bkgErrPdfList, bkgErrFracList, true);
	std::unique_ptr<RooDataHist> bkgErrHistData;
	std::unique_ptr<RooHistPdf> bkgErrPdfHist;
	RooAbsPdf *bkgErrPdf = bkgErrPdfAnalytic.get();

		RooRealVar sigTimeErrGaus1Mean("sigErrGaus1Mean", "sigErrGaus1Mean", readErrVar(signalErrResult, "sigErrGaus1Mean", 0.02), errFitLow, errUpper);
		RooRealVar sigTimeErrGaus1Sigma("sigErrGaus1Sigma", "sigErrGaus1Sigma", readErrVar(signalErrResult, "sigErrGaus1Sigma", 0.003), 1e-6, errUpper);
		RooRealVar sigTimeErrGaus2Mean("sigErrGaus2Mean", "sigErrGaus2Mean", readErrVar(signalErrResult, "sigErrGaus2Mean", sigTimeErrGaus1Mean.getVal()), errFitLow, errUpper);
		RooRealVar sigTimeErrGaus2Sigma("sigErrGaus2Sigma", "sigErrGaus2Sigma", readErrVar(signalErrResult, "sigErrGaus2Sigma", sigTimeErrGaus1Sigma.getVal()), 1e-6, errUpper);
			RooRealVar sigTimeErrTailMpv("sigErrTailMpv", "sigErrTailMpv", readErrVar(signalErrResult, "sigErrTailMpv", 0.02), errFitLow, errUpper);
			RooRealVar sigTimeErrTailWidth("sigErrTailWidth", "sigErrTailWidth", readErrVar(signalErrResult, "sigErrTailWidth", 0.002), 1e-6, errUpper);
			RooRealVar sigTimeErrTail2Mpv("sigErrTail2Mpv", "sigErrTail2Mpv", readErrVar(signalErrResult, "sigErrTail2Mpv", sigTimeErrTailMpv.getVal()), errFitLow, errUpper);
			RooRealVar sigTimeErrTail2Width("sigErrTail2Width", "sigErrTail2Width", readErrVar(signalErrResult, "sigErrTail2Width", sigTimeErrTailWidth.getVal()), 1e-6, errUpper);
			RooRealVar sigTimeErrLognM0("sigErrLognM0", "sigErrLognM0", readErrVar(signalErrResult, "sigErrLognM0", 0.02), errFitLow, errUpper);
			RooRealVar sigTimeErrLognK("sigErrLognK", "sigErrLognK", readErrVar(signalErrResult, "sigErrLognK", 0.35), 0.05, 0.95);
			RooRealVar sigTimeErrCore1FracRatio("sigErrCore1FracRatio", "sigErrCore1FracRatio", readErrVar(signalErrResult, "sigErrCore1FracRatio", fracToRatio(0.6)), 1e-6, 1e6);
			RooFormulaVar sigTimeErrCore1Frac("sigErrCore1Frac", "@0/(1.0+@0)", RooArgList(sigTimeErrCore1FracRatio));
			RooRealVar sigTimeErrTailFracRatio("sigErrTailFracRatio", "sigErrTailFracRatio", readErrVar(signalErrResult, "sigErrTailFracRatio", fracToRatio(0.08)), 1e-6, 1e6);
			RooFormulaVar sigTimeErrTailFrac("sigErrTailFrac", "@0/(1.0+@0)", RooArgList(sigTimeErrTailFracRatio));
			RooRealVar sigTimeErrTail2FracRatio("sigErrTail2FracRatio", "sigErrTail2FracRatio", readErrVar(signalErrResult, "sigErrTail2FracRatio", fracToRatio(0.05)), 1e-6, 1e6);
			RooFormulaVar sigTimeErrTail2Frac("sigErrTail2Frac", "@0/(1.0+@0)", RooArgList(sigTimeErrTail2FracRatio));
			RooRealVar sigTimeErrLognFracRatio("sigErrLognFracRatio", "sigErrLognFracRatio", readErrVar(signalErrResult, "sigErrLognFracRatio", fracToRatio(0.03)), 1e-6, 1e6);
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
			if (sigErrPdfList.getSize() < nErrSigComponents) sigErrFracList.add(sigTimeErrTail2Frac);
		}
		if (useErrSigLogn) {
			sigErrPdfList.add(sigTimeErrLogn);
			if (sigErrPdfList.getSize() < nErrSigComponents) sigErrFracList.add(sigTimeErrLognFrac);
		}
		if (useErrSigGaus1) {
			sigErrPdfList.add(sigTimeErrGaus1);
			if (sigErrPdfList.getSize() < nErrSigComponents) sigErrFracList.add(sigTimeErrCore1Frac);
		}
		if (useErrSigGaus2)
			sigErrPdfList.add(sigTimeErrGaus2);
		std::unique_ptr<RooAddPdf> signalErrPdfAnalytic = std::make_unique<RooAddPdf>("pdfErrSig", "pdfErrSig",
			sigErrPdfList, sigErrFracList, true);
	std::unique_ptr<RooDataHist> signalErrHistTemplate;
	std::unique_ptr<RooHistPdf> signalErrPdfHist;
	RooAbsPdf *signalErrPdf = signalErrPdfAnalytic.get();

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
	if (nBkgExpComponents + nBkgChebyOrder <= 0) {
		std::cerr << "ERROR: no background component configured in mass.C output." << std::endl;
		return;
	}
	if (nBkgExpComponents > 0 && nBkgChebyOrder > 0) {
		std::cerr << "ERROR: inconsistent background model configuration in mass.C output." << std::endl;
		return;
	}

	RooRealVar signal_mass_mean("signal_mass_mean","signal_mass_mean",readDoubleParam(*massFile,"signal_mass_mean",3.096),3.05,3.15);
	RooRealVar signal_mass_sigma("signal_mass_sigma","signal_mass_sigma",readDoubleParam(*massFile,"signal_mass_sigma",0.03),0.001,0.080);
	const double sigmaFloor = 1e-4;
	const double signal_mass_sigma2_init = std::max(signal_mass_sigma.getVal() + 0.03, sigmaFloor + 1e-4);
	RooRealVar signal_mass_sigma2_log_offset(
			"signal_mass_sigma2_log_offset", "log(signal_mass_sigma2 - floor)",
			std::log(std::max(readDoubleParam(*massFile,"signal_mass_sigma2",signal_mass_sigma2_init) - sigmaFloor, 1e-4)),
			std::log(1e-4),
			std::log(1.0));
	RooFormulaVar signal_mass_sigma2("signal_mass_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_sigma2_log_offset));
	RooFormulaVar signal_mass_sigma_delta2("signal_mass_sigma_delta2","@0-@1",RooArgList(signal_mass_sigma2,signal_mass_sigma));
	RooRealVar signal_mass_cb_sigma_base("signal_mass_cb_sigma_base","signal_mass_cb_sigma_base",readDoubleParam(*massFile,"signal_mass_cb_sigma",0.035),0.008,0.080);
	std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_from_gaus;
	std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_delta_from_gaus;
	std::unique_ptr<RooRealVar> signal_mass_cb_sigma_log_offset;
	RooAbsReal *signal_mass_cb_sigma = &signal_mass_cb_sigma_base;
	RooAbsReal *signal_mass_cb_sigma_delta = nullptr;
	if (nSignalGaussComponents >= 1) {
		const double signal_mass_cb_sigma_init = std::max(signal_mass_sigma.getVal() + 0.005, sigmaFloor + 1e-4);
		signal_mass_cb_sigma_log_offset = std::make_unique<RooRealVar>(
				"signal_mass_cb_sigma_log_offset", "log(signal_mass_cb_sigma - floor)",
				std::log(std::max(signal_mass_cb_sigma_init - sigmaFloor, 1e-4)),
				std::log(1e-4),
				std::log(1.0));
		signal_mass_cb_sigma_from_gaus = std::make_unique<RooFormulaVar>(
				"signal_mass_cb_sigma", Form("%g + exp(@0)", sigmaFloor),
				RooArgList(*signal_mass_cb_sigma_log_offset));
		signal_mass_cb_sigma_delta_from_gaus = std::make_unique<RooFormulaVar>("signal_mass_cb_sigma_delta","@0-@1",RooArgList(*signal_mass_cb_sigma_from_gaus,signal_mass_sigma));
		signal_mass_cb_sigma = signal_mass_cb_sigma_from_gaus.get();
		signal_mass_cb_sigma_delta = signal_mass_cb_sigma_delta_from_gaus.get();
		signal_mass_cb_sigma_log_offset->setVal(std::log(std::max(readDoubleParam(*massFile,"signal_mass_cb_sigma", signal_mass_cb_sigma->getVal()) - sigmaFloor, 1e-4)));
	} else {
		signal_mass_cb_sigma_delta = &signal_mass_cb_sigma_base;
	}
	const double signal_mass_cb_sigma2_init = std::max(signal_mass_cb_sigma->getVal() + 0.02, sigmaFloor + 1e-4);
	RooRealVar signal_mass_cb_sigma2_log_offset(
			"signal_mass_cb_sigma2_log_offset", "log(signal_mass_cb_sigma2 - floor)",
			std::log(std::max(signal_mass_cb_sigma2_init - sigmaFloor, 1e-4)),
			std::log(1e-4),
			std::log(1.0));
	RooFormulaVar signal_mass_cb_sigma2("signal_mass_cb_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_cb_sigma2_log_offset));
	RooFormulaVar signal_mass_cb_sigma_delta2("signal_mass_cb_sigma_delta2","@0-@1",RooArgList(signal_mass_cb_sigma2,*signal_mass_cb_sigma));
	RooRealVar signal_mass_cb_alpha("signal_mass_cb_alpha","signal_mass_cb_alpha",readDoubleParam(*massFile,"signal_mass_cb_alpha",1.5),0.1,10.0);
	RooRealVar signal_mass_cb_alpha2("signal_mass_cb_alpha2","signal_mass_cb_alpha2",readDoubleParam(*massFile,"signal_mass_cb_alpha2",2.0),0.1,10.0);
	RooRealVar signal_mass_cb_n("signal_mass_cb_n","signal_mass_cb_n",readDoubleParam(*massFile,"signal_mass_cb_n",3.0),0.1,20.0);
	RooRealVar signal_mass_cb_n2("signal_mass_cb_n2","signal_mass_cb_n2",readDoubleParam(*massFile,"signal_mass_cb_n2",4.0),0.2,10.0);
	RooRealVar signal_mass_frac_ratio1("signal_mass_frac_ratio1","signal_mass_frac_ratio1",readDoubleParam(*massFile,"signal_mass_frac_ratio1",1.86),1e-3,1e3);
	RooFormulaVar signal_mass_frac1("signal_mass_frac1","@0/(1.0+@0)",RooArgList(signal_mass_frac_ratio1));
	RooRealVar signal_mass_frac_ratio2("signal_mass_frac_ratio2","signal_mass_frac_ratio2",readDoubleParam(*massFile,"signal_mass_frac_ratio2",0.25),1e-3,1e3);
	RooFormulaVar signal_mass_frac2("signal_mass_frac2","@0/(1.0+@0)",RooArgList(signal_mass_frac_ratio2));
	RooRealVar signal_mass_frac_ratio3("signal_mass_frac_ratio3","signal_mass_frac_ratio3",readDoubleParam(*massFile,"signal_mass_frac_ratio3",0.111111),1e-3,1e3);
	RooFormulaVar signal_mass_frac3("signal_mass_frac3","@0/(1.0+@0)",RooArgList(signal_mass_frac_ratio3));

	if (nSignalGaussComponents >= 1)
		signal_mass_sigma.setVal(readDoubleParam(*massFile,"signal_mass_sigma",signal_mass_sigma.getVal()));
	signal_mass_sigma2_log_offset.setVal(std::log(std::max(readDoubleParam(*massFile,"signal_mass_sigma2",signal_mass_sigma2.getVal()) - sigmaFloor, 1e-4)));
	signal_mass_cb_sigma_base.setVal(readDoubleParam(*massFile,"signal_mass_cb_sigma",signal_mass_cb_sigma_base.getVal()));
	signal_mass_cb_sigma2_log_offset.setVal(std::log(std::max(readDoubleParam(*massFile,"signal_mass_cb_sigma2",signal_mass_cb_sigma2.getVal()) - sigmaFloor, 1e-4)));
	signal_mass_cb_alpha.setVal(readDoubleParam(*massFile,"signal_mass_cb_alpha",signal_mass_cb_alpha.getVal()));
	signal_mass_cb_alpha2.setVal(readDoubleParam(*massFile,"signal_mass_cb_alpha2",signal_mass_cb_alpha2.getVal()));
	signal_mass_cb_n.setVal(readDoubleParam(*massFile,"signal_mass_cb_n",signal_mass_cb_n.getVal()));
	signal_mass_cb_n2.setVal(readDoubleParam(*massFile,"signal_mass_cb_n2",signal_mass_cb_n2.getVal()));
	signal_mass_frac_ratio1.setVal(readDoubleParam(*massFile,"signal_mass_frac_ratio1",signal_mass_frac_ratio1.getVal()));
	signal_mass_frac_ratio2.setVal(readDoubleParam(*massFile,"signal_mass_frac_ratio2",signal_mass_frac_ratio2.getVal()));
	signal_mass_frac_ratio3.setVal(readDoubleParam(*massFile,"signal_mass_frac_ratio3",signal_mass_frac_ratio3.getVal()));
	signal_mass_mean.setVal(readDoubleParam(*massFile,"signal_mass_mean",signal_mass_mean.getVal()));

	const auto fixMassVar = [](RooRealVar &var) {
		var.setConstant(true);
	};
	fixMassVar(signal_mass_mean);
	fixMassVar(signal_mass_sigma);
	fixMassVar(signal_mass_cb_sigma_base);
	fixMassVar(signal_mass_cb_alpha);
	fixMassVar(signal_mass_cb_alpha2);
	fixMassVar(signal_mass_cb_n);
	fixMassVar(signal_mass_cb_n2);
	fixMassVar(signal_mass_cb_sigma2_log_offset);
	fixMassVar(signal_mass_sigma2_log_offset);
	fixMassVar(signal_mass_frac_ratio1);
	fixMassVar(signal_mass_frac_ratio2);
	fixMassVar(signal_mass_frac_ratio3);
	if (signal_mass_cb_sigma_log_offset)
		signal_mass_cb_sigma_log_offset->setConstant(true);
	

	RooGaussian signal_mass_gaus("signal_mass_gaus","signal_mass_gaus",obs_mass,signal_mass_mean,signal_mass_sigma);
	RooGaussian signal_mass_gaus2("signal_mass_gaus2","signal_mass_gaus2",obs_mass,signal_mass_mean,signal_mass_sigma2);
	std::unique_ptr<RooAbsPdf> signal_mass_cb_owned;
	if (nSignalCBComponents >= 2) {
		signal_mass_cb_owned = std::make_unique<RooCrystalBall>(
			"signal_mass_cb","signal_mass_cb",
			obs_mass,signal_mass_mean,
			*signal_mass_cb_sigma,signal_mass_cb_sigma2,
			signal_mass_cb_alpha,signal_mass_cb_n,
			signal_mass_cb_alpha2,signal_mass_cb_n2);
	} else if (nSignalCBComponents >= 1) {
		signal_mass_cb_owned = std::make_unique<RooCBShape>(
			"signal_mass_cb","signal_mass_cb",
			obs_mass,signal_mass_mean,*signal_mass_cb_sigma,signal_mass_cb_alpha,signal_mass_cb_n);
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

	const double chebyRawLimit = 4.0;
	RooRealVar bkg_mass_p1_raw("bkg_mass_p1_raw","bkg_mass_p1_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooRealVar bkg_mass_p2_raw("bkg_mass_p2_raw","bkg_mass_p2_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooRealVar bkg_mass_p3_raw("bkg_mass_p3_raw","bkg_mass_p3_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooRealVar bkg_mass_p4_raw("bkg_mass_p4_raw","bkg_mass_p4_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooRealVar bkg_mass_p5_raw("bkg_mass_p5_raw","bkg_mass_p5_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooRealVar bkg_mass_p6_raw("bkg_mass_p6_raw","bkg_mass_p6_raw",0.0,-chebyRawLimit,chebyRawLimit);
	RooFormulaVar bkg_mass_p1("bkg_mass_p1","tanh(@0)",RooArgList(bkg_mass_p1_raw));
	RooFormulaVar bkg_mass_p2("bkg_mass_p2","tanh(@0)",RooArgList(bkg_mass_p2_raw));
	RooFormulaVar bkg_mass_p3("bkg_mass_p3","tanh(@0)",RooArgList(bkg_mass_p3_raw));
	RooFormulaVar bkg_mass_p4("bkg_mass_p4","tanh(@0)",RooArgList(bkg_mass_p4_raw));
	RooFormulaVar bkg_mass_p5("bkg_mass_p5","tanh(@0)",RooArgList(bkg_mass_p5_raw));
	RooFormulaVar bkg_mass_p6("bkg_mass_p6","tanh(@0)",RooArgList(bkg_mass_p6_raw));
	RooRealVar bkg_mass_lambda("bkg_mass_lambda","bkg_mass_lambda",readDoubleParam(*massFile,"bkg_mass_lambda",-1.0),-10.0,-1e-4);
	if (nBkgChebyOrder >= 1) bkg_mass_p1_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p1_raw",bkg_mass_p1_raw.getVal()));
	if (nBkgChebyOrder >= 2) bkg_mass_p2_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p2_raw",bkg_mass_p2_raw.getVal()));
	if (nBkgChebyOrder >= 3) bkg_mass_p3_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p3_raw",bkg_mass_p3_raw.getVal()));
	if (nBkgChebyOrder >= 4) bkg_mass_p4_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p4_raw",bkg_mass_p4_raw.getVal()));
	if (nBkgChebyOrder >= 5) bkg_mass_p5_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p5_raw",bkg_mass_p5_raw.getVal()));
	if (nBkgChebyOrder >= 6) bkg_mass_p6_raw.setVal(readDoubleParam(*massFile,"bkg_mass_p6_raw",bkg_mass_p6_raw.getVal()));
	if (nBkgChebyOrder >= 1) fixMassVar(bkg_mass_p1_raw);
	if (nBkgChebyOrder >= 2) fixMassVar(bkg_mass_p2_raw);
	if (nBkgChebyOrder >= 3) fixMassVar(bkg_mass_p3_raw);
	if (nBkgChebyOrder >= 4) fixMassVar(bkg_mass_p4_raw);
	if (nBkgChebyOrder >= 5) fixMassVar(bkg_mass_p5_raw);
	if (nBkgChebyOrder >= 6) fixMassVar(bkg_mass_p6_raw);
	fixMassVar(bkg_mass_lambda);
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
		0.5,
		3.0);
	RooConstVar ctauTime2Mean("ctauTime2Mean","ctauTime2Mean",readDoubleParam(*prFile,"ctauTime2Mean",0.0));
	const double ctauTime2ScaleVal = readDoubleParam(*prFile,"ctauTime2Scale",1.0);
	const double ctauTime2DeltaVal = std::max(0.05, ctauTime2ScaleVal - ctauTime1ScaleVal);
	RooRealVar ctauTime2Delta("ctauTime2Delta","ctauTime2Delta",
		ctauTime2DeltaVal,
		0.05,
		5.0);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale","@0+@1",RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooConstVar ctauTime3Mean("ctauTime3Mean","ctauTime3Mean",readDoubleParam(*prFile,"ctauTime3Mean",0.0));
	const double ctauTime3ScaleVal = readDoubleParam(*prFile,"ctauTime3Scale",1.0);
	const double ctauTime3DeltaVal = std::max(0.05, ctauTime3ScaleVal - ctauTime2ScaleVal);
	RooRealVar ctauTime3Delta("ctauTime3Delta","ctauTime3Delta",
		ctauTime3DeltaVal,
		0.05,
		10.0);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale","@0+@1",RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooConstVar ctauTime4Mean("ctauTime4Mean","ctauTime4Mean",readDoubleParam(*prFile,"ctauTime4Mean",0.0));
	const double ctauTime4ScaleVal = readDoubleParam(*prFile,"ctauTime4Scale",1.0);
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	RooRealVar ctauTime4Delta("ctauTime4Delta","ctauTime4Delta",
		ctauTime4DeltaVal,
		0.05,
		15.0);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale","@0+@1",RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooConstVar ctauFrac1("ctauFrac1","ctauFrac1",readDoubleParam(*prFile,"ctauFrac1",0.5));
	RooConstVar ctauFrac2("ctauFrac2","ctauFrac2",readDoubleParam(*prFile,"ctauFrac2",0.2));
	RooConstVar ctauFrac3("ctauFrac3","ctauFrac3",readDoubleParam(*prFile,"ctauFrac3",0.1));

	RooGaussModel ctauTime1("ctauTime1","ctauTime1",obs_time,ctauTime1Mean,ctauTime1Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime2("ctauTime2","ctauTime2",obs_time,ctauTime2Mean,ctauTime2Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime3("ctauTime3","ctauTime3",obs_time,ctauTime3Mean,ctauTime3Scale,ctauMeanScale,obs_timeErr);
	RooGaussModel ctauTime4("ctauTime4","ctauTime4",obs_time,ctauTime4Mean,ctauTime4Scale,ctauMeanScale,obs_timeErr);

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
	const double signalLifetimeFloorSaved = readDoubleParam(*npFile, "signalLifetimeFloor", signalLifetimeFloor);
	const double signalLifetime2FloorSaved = readDoubleParam(*npFile, "signalLifetime2Floor", std::max(0.75 * signalLifetimeFloorSaved, 1e-2));
	const double signalLifetime3FloorSaved = readDoubleParam(*npFile, "signalLifetime3Floor", std::max(1.00 * signalLifetimeFloorSaved, 1e-2));
	const double signalLogLifetimeOffsetVal = readDoubleParam(
		*npFile,
		"signal_log_lifetime_offset",
		std::log(std::max(signalLifetime1Val - signalLifetimeFloorSaved, 1e-3)));
	const double signal2LogLifetimeOffsetVal = readDoubleParam(
		*npFile,
		"signal2_log_lifetime_offset",
		std::log(std::max(signalLifetime2Val - signalLifetime2FloorSaved, 1e-3)));
	const double signal3LogLifetimeOffsetVal = readDoubleParam(
		*npFile,
		"signal3_log_lifetime_offset",
		std::log(std::max(signalLifetime3Val - signalLifetime3FloorSaved, 1e-3)));
	const bool useSignalSS2 = (nSignalSSComponentsSaved >= 2) && std::isfinite(signalLifetime2Val) && signalLifetime2Val > signalLifetime2FloorSaved;
	const bool useSignalSS3 = useSignalSS2 && (nSignalSSComponentsSaved >= 3) && std::isfinite(signalLifetime3Val) && signalLifetime3Val > signalLifetime3FloorSaved;
	const int nSignalSSComponents = 1 + (useSignalSS2 ? 1 : 0) + (useSignalSS3 ? 1 : 0);
	const double signalLifetimeInit = std::max(std::max(signalLifetime1Val, signalLifetimeFloorSaved + 1e-3), std::max(0.08, 1.5 * signalLifetimeFloorSaved));
	const double signalLifetimeCeil = std::max(signalLifetimeInit * 5.0, maxStableLifetime);
	RooRealVar signal_log_lifetime_offset(
		"signal_log_lifetime_offset", "log(signal_lifetime - floor)",
		signalLogLifetimeOffsetVal,
		std::log(1e-2),
		std::log(std::max(signalLifetimeCeil - signalLifetimeFloorSaved, 2e-4)));
	RooFormulaVar signal_lifetime(
		"signal_lifetime",
		Form("%g + exp(@0)", signalLifetimeFloorSaved),
		RooArgList(signal_log_lifetime_offset));
	const double signalLifetime2Floor = signalLifetime2FloorSaved;
	const double signalLifetime2Init = std::max(std::max(signalLifetime2Val, signalLifetime2Floor + 1e-3), std::max(0.20, 2.0 * signalLifetime2Floor));
	const double signalLifetime2Ceil = std::max(signalLifetime2Init * 5.0, maxStableLifetime);
	RooRealVar signal2_log_lifetime_offset(
		"signal2_log_lifetime_offset", "log(signal2_lifetime - floor)",
		signal2LogLifetimeOffsetVal,
		std::log(1e-2),
		std::log(std::max(signalLifetime2Ceil - signalLifetime2Floor, 2e-4)));
	RooFormulaVar signal2_lifetime(
		"signal2_lifetime",
		Form("%g + exp(@0)", signalLifetime2Floor),
		RooArgList(signal2_log_lifetime_offset));
	const double signalLifetime3Floor = signalLifetime3FloorSaved;
	const double signalLifetime3Init = std::max(std::max(signalLifetime3Val, signalLifetime3Floor + 1e-3), std::max(0.40, 2.0 * signalLifetime3Floor));
	const double signalLifetime3Ceil = std::max(signalLifetime3Init * 5.0, maxStableLifetime);
	RooRealVar signal3_log_lifetime_offset(
		"signal3_log_lifetime_offset", "log(signal3_lifetime - floor)",
		signal3LogLifetimeOffsetVal,
		std::log(1e-2),
		std::log(std::max(signalLifetime3Ceil - signalLifetime3Floor, 2e-4)));
	RooFormulaVar signal3_lifetime(
		"signal3_lifetime",
		Form("%g + exp(@0)", signalLifetime3Floor),
		RooArgList(signal3_log_lifetime_offset));
	RooConstVar NsignalSS1("NsignalSS1","NsignalSS1",npYield1);
	RooConstVar NsignalSS2("NsignalSS2","NsignalSS2",npYield2);
	RooConstVar NsignalSS3("NsignalSS3","NsignalSS3",npYield3);

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

	RooRealVar bkg_log_lifetime_offset("bkg_log_lifetime_offset","bkg_log_lifetime_offset",readResultValue(*bkgTimeResult,"bkg_log_lifetime_offset",std::log(std::max(0.05, 1.5 * bkgLifetimeFloor) - bkgLifetimeFloor)));
	RooRealVar bkg_log_lifetime2_offset("bkg_log_lifetime2_offset","bkg_log_lifetime2_offset",readResultValue(*bkgTimeResult,"bkg_log_lifetime2_offset",std::log(std::max(0.04, 1.5 * std::max(1.5 * bkgLifetimeFloor, 2e-2)) - std::max(1.5 * bkgLifetimeFloor, 2e-2))));
	RooRealVar bkg_log_lifetime3_offset("bkg_log_lifetime3_offset","bkg_log_lifetime3_offset",readResultValue(*bkgTimeResult,"bkg_log_lifetime3_offset",std::log(std::max(0.10, 1.5 * std::max(2.0 * bkgLifetimeFloor, 3e-2)) - std::max(2.0 * bkgLifetimeFloor, 3e-2))));
	RooRealVar bkg_sym_log_lifetime_offset("bkg_sym_log_lifetime_offset","bkg_sym_log_lifetime_offset",readResultValue(*bkgTimeResult,"bkg_sym_log_lifetime_offset",std::log(std::max(0.04, 1.4 * std::max(1.5 * errRange.first, 1e-2)) - std::max(1.5 * errRange.first, 1e-2))));
	RooRealVar bkg_sym_log_lifetime2_offset("bkg_sym_log_lifetime2_offset","bkg_sym_log_lifetime2_offset",readResultValue(*bkgTimeResult,"bkg_sym_log_lifetime2_offset",std::log(std::max(0.07, 1.4 * std::max(2.0 * std::max(1.5 * errRange.first, 1e-2), 2e-2)) - std::max(2.0 * std::max(1.5 * errRange.first, 1e-2), 2e-2))));
	RooRealVar bkg_sym_log_lifetime3_offset("bkg_sym_log_lifetime3_offset","bkg_sym_log_lifetime3_offset",readResultValue(*bkgTimeResult,"bkg_sym_log_lifetime3_offset",std::log(std::max(0.12, 1.3 * std::max(3.0 * std::max(1.5 * errRange.first, 1e-2), 4e-2)) - std::max(3.0 * std::max(1.5 * errRange.first, 1e-2), 4e-2))));
	RooRealVar bkg_flip_log_lifetime_offset("bkg_flip_log_lifetime_offset","bkg_flip_log_lifetime_offset",readResultValue(*bkgTimeResult,"bkg_flip_log_lifetime_offset",std::log(std::max(0.05, 1.2 * std::max(0.75 * bkgLifetimeFloor, 5e-3)) - std::max(0.75 * bkgLifetimeFloor, 5e-3))));
	RooRealVar bkg_flip_log_lifetime2_offset("bkg_flip_log_lifetime2_offset","bkg_flip_log_lifetime2_offset",readResultValue(*bkgTimeResult,"bkg_flip_log_lifetime2_offset",std::log(std::max(0.08, 1.3 * std::max(1.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 1e-2)) - std::max(1.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 1e-2))));
	RooRealVar bkg_flip_log_lifetime3_offset("bkg_flip_log_lifetime3_offset","bkg_flip_log_lifetime3_offset",readResultValue(*bkgTimeResult,"bkg_flip_log_lifetime3_offset",std::log(std::max(0.12, 1.2 * std::max(2.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 2e-2)) - std::max(2.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 2e-2))));
	bkg_log_lifetime_offset.setConstant(false);
	bkg_log_lifetime2_offset.setConstant(false);
	bkg_log_lifetime3_offset.setConstant(false);
	bkg_sym_log_lifetime_offset.setConstant(true);
	bkg_sym_log_lifetime2_offset.setConstant(true);
	bkg_sym_log_lifetime3_offset.setConstant(true);
	bkg_flip_log_lifetime_offset.setConstant(false);
	bkg_flip_log_lifetime2_offset.setConstant(false);
	bkg_flip_log_lifetime3_offset.setConstant(false);

	const double bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 2e-2);
	const double bkgLifetime3Floor = std::max(2.0 * bkgLifetimeFloor, 3e-2);
	const double bkgSymLifetimeFloor = std::max(1.5 * errRange.first, 1e-2);
	const double bkgSymLifetime2Floor = std::max(2.0 * bkgSymLifetimeFloor, 2e-2);
	const double bkgSymLifetime3Floor = std::max(3.0 * bkgSymLifetimeFloor, 4e-2);
	const double bkgFlipLifetimeFloor = std::max(0.75 * bkgLifetimeFloor, 5e-3);
	const double bkgFlipLifetime2Floor = std::max(1.5 * bkgFlipLifetimeFloor, 1e-2);
	const double bkgFlipLifetime3Floor = std::max(2.5 * bkgFlipLifetimeFloor, 2e-2);
	const double bkgLifetimeCeil = std::max(0.50, maxStableLifetime);
	const double bkgLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetimeCeil = std::max(0.30, 0.5 * maxStableLifetime);
	const double bkgSymLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgSymLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetimeMinGap = std::max(0.30 * bkgSymLifetimeFloor, 2e-3);
	const double bkgSymLifetime2MinGap = std::max(0.25 * bkgSymLifetime2Floor, 3e-3);
	const double bkgSymLifetime3MinGap = std::max(0.20 * bkgSymLifetime3Floor, 4e-3);
	bkg_log_lifetime_offset.setRange(
		std::log(1e-4),
		std::log(std::max(bkgLifetimeCeil - bkgLifetimeFloor, 2e-4)));
	bkg_log_lifetime2_offset.setRange(
		std::log(1e-4),
		std::log(std::max(bkgLifetime2Ceil - bkgLifetime2Floor, 2e-4)));
	bkg_log_lifetime3_offset.setRange(
		std::log(1e-4),
		std::log(std::max(bkgLifetime3Ceil - bkgLifetime3Floor, 2e-4)));
	bkg_sym_log_lifetime_offset.setRange(
		std::log(bkgSymLifetimeMinGap),
		std::log(std::max(bkgSymLifetimeCeil - bkgSymLifetimeFloor, 2.0 * bkgSymLifetimeMinGap)));
	bkg_sym_log_lifetime2_offset.setRange(
		std::log(bkgSymLifetime2MinGap),
		std::log(std::max(bkgSymLifetime2Ceil - bkgSymLifetime2Floor, 2.0 * bkgSymLifetime2MinGap)));
	bkg_sym_log_lifetime3_offset.setRange(
		std::log(bkgSymLifetime3MinGap),
		std::log(std::max(bkgSymLifetime3Ceil - bkgSymLifetime3Floor, 2.0 * bkgSymLifetime3MinGap)));
	bkg_flip_log_lifetime_offset.setRange(
		std::log(1e-4),
		std::log(std::max(std::max(0.50, maxStableLifetime) - bkgFlipLifetimeFloor, 2e-4)));
	bkg_flip_log_lifetime2_offset.setRange(
		std::log(1e-4),
		std::log(std::max(std::max(0.60, maxStableLifetime) - bkgFlipLifetime2Floor, 2e-4)));
	bkg_flip_log_lifetime3_offset.setRange(
		std::log(1e-4),
		std::log(std::max(std::max(0.80, maxStableLifetime) - bkgFlipLifetime3Floor, 2e-4)));
	RooFormulaVar bkg_lifetime("bkg_lifetime",Form("%g + exp(@0)", bkgLifetimeFloor),RooArgList(bkg_log_lifetime_offset));
	RooFormulaVar bkg_lifetime2("bkg_lifetime2",Form("%g + exp(@0)", bkgLifetime2Floor),RooArgList(bkg_log_lifetime2_offset));
	RooFormulaVar bkg_lifetime3("bkg_lifetime3",Form("%g + exp(@0)", bkgLifetime3Floor),RooArgList(bkg_log_lifetime3_offset));
	RooFormulaVar bkg_sym_lifetime("bkg_sym_lifetime",Form("%g + exp(@0)", bkgSymLifetimeFloor),RooArgList(bkg_sym_log_lifetime_offset));
	RooFormulaVar bkg_sym_lifetime2("bkg_sym_lifetime2",Form("%g + exp(@0)", bkgSymLifetime2Floor),RooArgList(bkg_sym_log_lifetime2_offset));
	RooFormulaVar bkg_sym_lifetime3("bkg_sym_lifetime3",Form("%g + exp(@0)", bkgSymLifetime3Floor),RooArgList(bkg_sym_log_lifetime3_offset));
	RooFormulaVar bkg_flip_lifetime("bkg_flip_lifetime",Form("%g + exp(@0)", bkgFlipLifetimeFloor),RooArgList(bkg_flip_log_lifetime_offset));
	RooFormulaVar bkg_flip_lifetime2("bkg_flip_lifetime2",Form("%g + exp(@0)", bkgFlipLifetime2Floor),RooArgList(bkg_flip_log_lifetime2_offset));
	RooFormulaVar bkg_flip_lifetime3("bkg_flip_lifetime3",Form("%g + exp(@0)", bkgFlipLifetime3Floor),RooArgList(bkg_flip_log_lifetime3_offset));

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

	std::vector<std::unique_ptr<RooConstVar>> bkgCoeffVars;
	RooArgList bkgTimePdfList;
	RooArgList bkgTimeCoeffList;
	auto addBkgComponent = [&](RooAbsPdf &pdf, const char *name, double yield) {
		bkgTimePdfList.add(pdf);
		bkgCoeffVars.push_back(std::make_unique<RooConstVar>(name, name, std::max(0.0, yield)));
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

	std::unique_ptr<RooAbsData> signalErrData;
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
		}
	}
	if (!signalErrData || signalErrData->sumEntries() <= 0.0) {
		std::cerr << "ERROR: empty signal ctau3DErr template/input distribution" << std::endl;
		return;
	}
		if (isHistErrPdfChoice(bkgErrPdfOpt)) {
			bkgErrHistData = std::make_unique<RooDataHist>(
				"bkgErrHistData", "bkgErrHistData",
				RooArgSet(obs_timeErr), *bkgErrData);
			bkgErrPdfHist = std::make_unique<RooHistPdf>(
				"bkgErrHistPdf", "bkgErrHistPdf",
				RooArgSet(obs_timeErr), *bkgErrHistData, histPdfInterpolationOrder);
			bkgErrPdf = bkgErrPdfHist.get();
		}
		if (isHistErrPdfChoice(sigErrPdfOpt)) {
			auto *signalErrHist = dynamic_cast<RooDataHist *>(signalErrData.get());
			if (!signalErrHist) {
				std::cerr << "ERROR: signalErrData is not a RooDataHist for RooHistPdf mode" << std::endl;
				return;
		}
		signalErrHistTemplate = std::make_unique<RooDataHist>(*signalErrHist);
		signalErrPdfHist = std::make_unique<RooHistPdf>(
			"signalErrHistPdf", "signalErrHistPdf",
			RooArgSet(obs_timeErr), *signalErrHistTemplate, histPdfInterpolationOrder);
		signalErrPdf = signalErrPdfHist.get();
	}

	std::vector<std::unique_ptr<RooRealVar>> constraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> constraintPdfs;
	RooArgSet constraints;
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
	auto addLogLifetimeConstraint = [&](const char *baseName, RooRealVar &offsetVar, double floor, double central, double sigma)
	{
		if (!(std::isfinite(floor) && std::isfinite(central) && std::isfinite(sigma) && sigma > 0.0))
			return;
		if (central <= floor)
			central = floor + std::max(1e-3, 0.25 * floor);
		const double low = std::max(floor + 1e-6, central - sigma);
		const double high = central + sigma;
		const double centralOffset = std::log(std::max(central - floor, 1e-6));
		const double lowOffset = std::log(std::max(low - floor, 1e-6));
		const double highOffset = std::log(std::max(high - floor, 1e-6));
		const double sigmaOffset = std::max(std::abs(highOffset - lowOffset) * 0.5, 1e-3);
		addConstraint(baseName, offsetVar, centralOffset, sigmaOffset);
	};

	addConstraint("ctauTime1Scale", ctauTime1Scale, ctauTime1ScaleVal, std::max(0.05, 0.10 * ctauTime1ScaleVal));
	if (nResolutionComponents >= 2)
		addConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal, std::max(0.05, 0.10 * ctauTime2DeltaVal));
	if (nResolutionComponents >= 3)
		addConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal, std::max(0.05, 0.10 * ctauTime3DeltaVal));
	if (nResolutionComponents >= 4)
		addConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal, std::max(0.05, 0.10 * ctauTime4DeltaVal));
	addLogLifetimeConstraint("signal_lifetime", signal_log_lifetime_offset, signalLifetimeFloor, signalLifetime1Val, std::max(0.01, 0.15 * signalLifetime1Val));
	if (useSignalSS2)
		addLogLifetimeConstraint("signal2_lifetime", signal2_log_lifetime_offset, signalLifetime2Floor, signalLifetime2Val, std::max(0.005, 0.15 * signalLifetime2Val));
	if (useSignalSS3)
		addLogLifetimeConstraint("signal3_lifetime", signal3_log_lifetime_offset, signalLifetime3Floor, signalLifetime3Val, std::max(0.005, 0.15 * signalLifetime3Val));
	addConstraint("bkg_log_lifetime_offset", bkg_log_lifetime_offset,
		readResultValue(*bkgTimeResult, "bkg_log_lifetime_offset", bkg_log_lifetime_offset.getVal()),
		std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime_offset")) ?
			dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime_offset"))->getError() : 0.1));
	if (nBkgSSComponents >= 2)
		addConstraint("bkg_log_lifetime2_offset", bkg_log_lifetime2_offset,
			readResultValue(*bkgTimeResult, "bkg_log_lifetime2_offset", bkg_log_lifetime2_offset.getVal()),
			std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime2_offset")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime2_offset"))->getError() : 0.1));
	if (nBkgSSComponents >= 3)
		addConstraint("bkg_log_lifetime3_offset", bkg_log_lifetime3_offset,
			readResultValue(*bkgTimeResult, "bkg_log_lifetime3_offset", bkg_log_lifetime3_offset.getVal()),
			std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime3_offset")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_log_lifetime3_offset"))->getError() : 0.1));
	if (nBkgFlipComponents >= 1)
		addConstraint("bkg_flip_log_lifetime_offset", bkg_flip_log_lifetime_offset,
			readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime_offset", bkg_flip_log_lifetime_offset.getVal()),
			std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime_offset")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime_offset"))->getError() : 0.1));
	if (nBkgFlipComponents >= 2)
		addConstraint("bkg_flip_log_lifetime2_offset", bkg_flip_log_lifetime2_offset,
			readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime2_offset", bkg_flip_log_lifetime2_offset.getVal()),
			std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime2_offset")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime2_offset"))->getError() : 0.1));
	if (nBkgFlipComponents >= 3)
		addConstraint("bkg_flip_log_lifetime3_offset", bkg_flip_log_lifetime3_offset,
			readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime3_offset", bkg_flip_log_lifetime3_offset.getVal()),
			std::max(1e-3, dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime3_offset")) ?
				dynamic_cast<RooRealVar *>(bkgTimeResult->floatParsFinal().find("bkg_flip_log_lifetime3_offset"))->getError() : 0.1));

	// ===== Plotting Only Model =====
	// only used to plot prompt and nonprompt components separately
	// These are not used for fitting.
	RooProdPdf prompt_core("prompt_core","prompt_core",RooArgSet(*signal_mass_pdf, *signalErrPdf),Conditional(RooArgSet(*promptTimePdf),RooArgSet(obs_time)));
	RooProdPdf np_core("np_core","np_core",RooArgSet(*signal_mass_pdf, *signalErrPdf),Conditional(RooArgSet(*signal_np_time),RooArgSet(obs_time)));
	// ===============================

	RooProdPdf signal_core("signal_core","signal_core",RooArgSet(*signal_mass_pdf, *signalErrPdf),Conditional(RooArgSet(signal_time),RooArgSet(obs_time)));
	RooProdPdf bkg_core("bkg_core","bkg_core",RooArgSet(*bkg_mass_pdf, *bkgErrPdf),Conditional(RooArgSet(*bkg_time_pdf),RooArgSet(obs_time)));
	RooAddPdf model("model","model",RooArgList(signal_core,bkg_core),RooArgList(Nsig,Nbkg));
	RooFitResult *model_result = nullptr;
	if (constraints.getSize() > 0)
		model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight), ExternalConstraints(constraints));
	else
		model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight));
	if (model_result && (model_result->status() != 0 || model_result->covQual() < 2))
	{
		const double refitSigmaScale = 1.0;
		RooArgSet retryConstraints;
		if (constraints.getSize() > 0)
			retryConstraints.add(constraints);
		RooArgList floatPars(model_result->floatParsFinal());
		auto modelVars = std::unique_ptr<RooArgSet>(model.getVariables());
		std::vector<std::unique_ptr<RooConstVar>> refitConsts;
		std::vector<std::unique_ptr<RooGaussian>> refitPdfs;
		for (int i = 0; i < floatPars.getSize(); ++i)
		{
			auto *src = dynamic_cast<RooRealVar *>(floatPars.at(i));
			if (!src || !modelVars)
				continue;
			auto *dst = dynamic_cast<RooRealVar *>(modelVars->find(src->GetName()));
			if (!dst)
				continue;
			const double val = src->getVal();
			const double err = refitSigmaScale * src->getError();
			if (!(err > 0.0))
				continue;
			refitConsts.push_back(std::make_unique<RooConstVar>(Form("refit_%s_mean", dst->GetName()), "", val));
			refitConsts.push_back(std::make_unique<RooConstVar>(Form("refit_%s_sigma", dst->GetName()), "", err));
			refitPdfs.push_back(std::make_unique<RooGaussian>(Form("refit_constr_%s", dst->GetName()), "",
				*dst, *refitConsts[refitConsts.size() - 2], *refitConsts.back()));
			retryConstraints.add(*refitPdfs.back());
		}
		if (retryConstraints.getSize() > 0)
		{
			delete model_result;
			model_result = model.fitTo(*data, Extended(), Save(), SumW2Error(isWeight), ExternalConstraints(retryConstraints));
			std::cout << "[OK] fit2d constrained refit done" << std::endl;
		}
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
	signalErrPdf->plotOn(sigTimeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (sigErrPdfOpt == kErrPdfAnalytic) {
		if (nErrSigLandau >= 1) signalErrPdf->plotOn(sigTimeErrPlot, Components(sigTimeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
		if (nErrSigLandau >= 2) signalErrPdf->plotOn(sigTimeErrPlot, Components(sigTimeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
		if (useErrSigLogn) signalErrPdf->plotOn(sigTimeErrPlot, Components(sigTimeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
		if (nErrSigGauss >= 1) signalErrPdf->plotOn(sigTimeErrPlot, Components(sigTimeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
		if (nErrSigGauss >= 2) signalErrPdf->plotOn(sigTimeErrPlot, Components(sigTimeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
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
		sigTimeErrLeg.AddEntry(o, isHistErrPdfChoice(sigErrPdfOpt) ? errPdfChoiceLabel(sigErrPdfOpt) : "Fit", "l");
	if (sigErrPdfOpt == kErrPdfAnalytic) {
		if (auto *o = findObj(sigTimeErrPlot, "tail"))
			sigTimeErrLeg.AddEntry(o, "Landau tail", "l");
		if (nErrSigLandau >= 2)
			if (auto *o = findObj(sigTimeErrPlot, "tail2"))
				sigTimeErrLeg.AddEntry(o, "Landau 2", "l");
		if (auto *o = findObj(sigTimeErrPlot, "gaus1"))
			sigTimeErrLeg.AddEntry(o, "Gauss 1", "l");
		if (nErrSigGauss >= 2)
			if (auto *o = findObj(sigTimeErrPlot, "gaus2"))
				sigTimeErrLeg.AddEntry(o, "Gauss 2", "l");
		if (useErrSigLogn)
			if (auto *o = findObj(sigTimeErrPlot, "logn"))
				sigTimeErrLeg.AddEntry(o, "Log-normal tail", "l");
	}
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
		int minimize = -1;
		int hesse = -1;
		if (signalErrResult)
		{
			for (UInt_t i = 0, n = signalErrResult->numStatusHistory(); i < n; ++i)
			{
				const char *lab = signalErrResult->statusLabelHistory(i);
				if (lab)
				{
					if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
						minimize = signalErrResult->statusCodeHistory(i);
					else if (!strcmp(lab, "HESSE"))
						hesse = signalErrResult->statusCodeHistory(i);
				}
			}
		}
		if (isHistErrPdfChoice(sigErrPdfOpt))
			tx.DrawLatex(0.19, 0.765, Form("Status : %s template (%d)", errPdfChoiceLabel(sigErrPdfOpt), histPdfInterpolationOrder));
		else
			tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
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
	const int sigTimeErrNPar = (isHistErrPdfChoice(sigErrPdfOpt) || !signalErrResult) ? 0 : signalErrResult->floatParsFinal().getSize();
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
		bkgErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), DataError(RooAbsData::SumW2), Name("data"));
	else
		bkgErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), Name("data"));
	bkgErrPdf->plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (bkgErrPdfOpt == kErrPdfAnalytic) {
		if (nErrBkgGauss >= 1) bkgErrPdf->plotOn(timeErrPlot, Components(bkgTimeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
		if (nErrBkgGauss >= 2) bkgErrPdf->plotOn(timeErrPlot, Components(bkgTimeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
		if (nErrBkgLandau >= 1) bkgErrPdf->plotOn(timeErrPlot, Components(bkgTimeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
		if (nErrBkgLandau >= 2) bkgErrPdf->plotOn(timeErrPlot, Components(bkgTimeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
		if (useErrBkgLogn) bkgErrPdf->plotOn(timeErrPlot, Components(bkgTimeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
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
		timeErrLeg.AddEntry(o, isHistErrPdfChoice(bkgErrPdfOpt) ? errPdfChoiceLabel(bkgErrPdfOpt) : "Fit", "l");
	if (bkgErrPdfOpt == kErrPdfAnalytic) {
		if (auto *o = findObj(timeErrPlot, "tail"))
			timeErrLeg.AddEntry(o, "Landau tail", "l");
		if (nErrBkgLandau >= 2)
			if (auto *o = findObj(timeErrPlot, "tail2"))
				timeErrLeg.AddEntry(o, "Landau 2", "l");
		if (auto *o = findObj(timeErrPlot, "gaus1"))
			timeErrLeg.AddEntry(o, "Gauss 1", "l");
		if (nErrBkgGauss >= 2)
			if (auto *o = findObj(timeErrPlot, "gaus2"))
				timeErrLeg.AddEntry(o, "Gauss 2", "l");
		if (useErrBkgLogn)
			if (auto *o = findObj(timeErrPlot, "logn"))
				timeErrLeg.AddEntry(o, "Log-normal tail", "l");
	}
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
		if (bkgErrResult)
		{
			for (UInt_t i = 0, n = bkgErrResult->numStatusHistory(); i < n; ++i)
			{
				const char *lab = bkgErrResult->statusLabelHistory(i);
				if (lab)
				{
					if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
						minimize = bkgErrResult->statusCodeHistory(i);
					else if (!strcmp(lab, "HESSE"))
						hesse = bkgErrResult->statusCodeHistory(i);
				}
			}
		}
		if (isHistErrPdfChoice(bkgErrPdfOpt))
			tx.DrawLatex(0.19, 0.765, Form("Status : %s template (%d)", errPdfChoiceLabel(bkgErrPdfOpt), histPdfInterpolationOrder));
		else
			tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
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
	const int timeErrNPar = isHistErrPdfChoice(bkgErrPdfOpt) ? 0 : (bkgErrResult ? bkgErrResult->floatParsFinal().getSize() : 0);
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
		int minimize = -1;
		int hesse = -1;
		if (model_result)
		{
			for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = model_result->statusLabelHistory(i);
				if (lab)
				{
					if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
						minimize = model_result->statusCodeHistory(i);
					else if (!strcmp(lab, "HESSE"))
						hesse = model_result->statusCodeHistory(i);
				}
			}
		}
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.73, y0 = 0.87, dy = -0.045;
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
			if (model_result)
			{
				const double err = var.getPropagatedError(*model_result);
				if (err > 0.0 && std::isfinite(err))
				{
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
					return;
				}
			}
			tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
		};
		printVar("N_{sig}", Nsig);
		printVar("N_{bkg}", Nbkg);
		printVar("f_{B}", bFraction);
		printVar("s_{1}", ctauTime1Scale);
		if (nResolutionComponents >= 2)
		{
			printVar("s_{2}", ctauTime2Scale);
			printVar("#Delta s_{21}", ctauTime2Delta);
		}
		if (nResolutionComponents >= 3)
		{
			printVar("s_{3}", ctauTime3Scale);
			printVar("#Delta s_{32}", ctauTime3Delta);
		}
		if (nResolutionComponents >= 4)
		{
			printVar("s_{4}", ctauTime4Scale);
			printVar("#Delta s_{43}", ctauTime4Delta);
		}
		printVar("#tau_{sig1}", signal_lifetime);
		if (useSignalSS2)
			printVar("#tau_{sig2}", signal2_lifetime);
		if (useSignalSS3)
			printVar("#tau_{sig3}", signal3_lifetime);
		printVar("#tau_{bkg1}", bkg_lifetime);
		if (nBkgSSComponents >= 2)
			printVar("#tau_{bkg2}", bkg_lifetime2);
		if (nBkgSSComponents >= 3)
			printVar("#tau_{bkg3}", bkg_lifetime3);
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", massChi.first, massNdf, massPvalue));
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
		int minimize = -1;
		int hesse = -1;
		if (model_result)
		{
			for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = model_result->statusLabelHistory(i);
				if (lab)
				{
					if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
						minimize = model_result->statusCodeHistory(i);
					else if (!strcmp(lab, "HESSE"))
						hesse = model_result->statusCodeHistory(i);
				}
			}
		}
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.73, y0 = 0.87, dy = -0.045;
		int k = 0;
		auto printVal = [&](const char *title, double val, double err = -1.0)
		{
			if (err >= 0.0)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, val, err));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, val));
		};
		printVal("f_{B}", bFraction.getVal(), bFraction.getError());
		printVal("N_{sig}", Nsig.getVal());
		printVal("N_{bkg}", Nbkg.getVal());
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", timeChi.first, timeNdf, timePvalue));
	}
	TLine timeLine(ctRange.first, 0.0, ctRange.second, 0.0);
	timeLine.SetLineStyle(2);
	timeLine.Draw("same");
	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

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
			TParameter<int>("bkgErrPdfOpt", bkgErrPdfOpt).Write();
			TParameter<int>("sigErrPdfOpt", sigErrPdfOpt).Write();
			TParameter<int>("histPdfInterpolationOrder", histPdfInterpolationOrder).Write();
			TParameter<int>("fitStatus", model_result ? model_result->status() : -1).Write();
			TParameter<int>("hesseStatus", hesse).Write();
			fitOut->Write();
		}
	}

	cout << "----------------- FIT RESULT FOR THE 2D MODEL ------------" << endl;
	cout << "bkgErrPdfOpt in fit2d: " << bkgErrPdfOpt
			 << " (" << errPdfChoiceLabel(bkgErrPdfOpt) << ")" << endl;
	cout << "sigErrPdfOpt in fit2d: " << sigErrPdfOpt
			 << " (" << errPdfChoiceLabel(sigErrPdfOpt) << ")" << endl;
	model_result->Print("v");
	const TString figErrSig = figName("errSig");
	const TString figErrBkg = figName("errBkg");
	const TString figMass = figName("mass_fit");
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] fit2d err sig fit : " << figErrSig << std::endl;
	std::cout << "[FIG] fit2d err bkg fit : " << figErrBkg << std::endl;
	std::cout << "[FIG] fit2d mass fit : " << figMass << std::endl;
	std::cout << "[FIG] fit2d lifetime fit : " << figLifetime << std::endl;
}
