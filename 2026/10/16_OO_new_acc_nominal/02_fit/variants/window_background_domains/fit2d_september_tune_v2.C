#include "../../../include/rule_override.h"
#include "../../../include/prompt_window_transport.h"
#include "../../../include/projection_diagnostics.h"
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
#include "RooMultiVarGaussian.h"
#include "TDecompChol.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooCrystalBall.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "RooAddModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooMinimizer.h"
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
#include "TStopwatch.h"
#include "TSystem.h"
#include "TString.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TMath.h"
#include "../../../include/saved_fit_helpers.h"
#include "../../../include/fit_config.h"
#include "RooHist.h"
#include "RooLognormal.h"
#include "RooMsgService.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace RooFit;

// The objective includes auxiliary-only coordinates absent from the data PDF.
inline std::unique_ptr<RooArgSet> objective_parameters(RooAbsPdf &model, const RooArgSet &constraints, const RooArgSet &observables)
{
	auto out = std::unique_ptr<RooArgSet>(model.getParameters(observables));
	for (auto *arg : constraints)
	{
		auto *pdf = dynamic_cast<RooAbsPdf *>(arg);
		if (!pdf)
			throw std::runtime_error("Constraint is not PDF");
		std::unique_ptr<RooArgSet> params(pdf->getParameters(observables));
		out->add(*params, true);
	}
	return out;
}
inline void restore_objective_result(const RooFitResult &result, RooArgSet &params, bool restoreErrors = true)
{
	for (auto *arg : result.floatParsFinal())
	{
		auto *source = dynamic_cast<const RooRealVar *>(arg);
		auto *target = dynamic_cast<RooRealVar *>(params.find(arg->GetName()));
		if (!source || !target)
			throw std::runtime_error(std::string("Missing objective float coordinate: ") + arg->GetName());
		if (!std::isfinite(source->getVal()) || source->getVal() < target->getMin() || source->getVal() > target->getMax())
			throw std::runtime_error("Saved coordinate outside unchanged domain");
		target->setVal(source->getVal());
		target->setConstant(false);
		if (restoreErrors && std::isfinite(source->getError()) && source->getError() >= 0)
			target->setError(source->getError());
	}
	for (auto *arg : result.constPars())
		if (auto *target = dynamic_cast<RooRealVar *>(params.find(arg->GetName())))
			if (auto *source = dynamic_cast<const RooRealVar *>(arg))
			{
				target->setVal(source->getVal());
				target->setConstant(true);
			}
}
inline void write_objective_state(const char *path, const RooArgSet &params)
{
	std::ofstream out(path);
	out << std::setprecision(17) << "name,value,error,min,max,constant\n";
	for (auto *arg : params)
		if (auto *v = dynamic_cast<RooRealVar *>(arg))
			out << v->GetName() << ',' << v->getVal() << ',' << v->getError() << ',' << v->getMin() << ',' << v->getMax() << ',' << v->isConstant() << '\n';
}

struct Fit2dScopedMacroTimer
{
	TStopwatch sw;
	TString label;

	Fit2dScopedMacroTimer(const char *macroName, float ptLow, float ptHigh, float yLow, float yHigh)
			: label(TString::Format("%s(pt=%.2f-%.2f, |y|=%.2f-%.2f)", macroName, ptLow, ptHigh, yLow, yHigh))
	{
		sw.Start();
		std::cout << "[TIMER] start " << label << std::endl;
	}

	~Fit2dScopedMacroTimer()
	{
		sw.Stop();
		std::cout << Form("[TIMER] done %s | real %.2f s | cpu %.2f s", label.Data(), sw.RealTime(), sw.CpuTime()) << std::endl;
	}
};
// 2026-07-31 22:06 KST: count only projection-relevant floating parameters in each reduced chi2 NDF.
struct ProjectionParameterSummary
{
	std::vector<std::string> massOnly;
	std::vector<std::string> lifetimeOnly;
	std::vector<std::string> common;
	std::vector<std::string> unknown;
};

static bool starts_with(const std::string &value, const char *prefix)
{
	return value.rfind(prefix, 0) == 0;
}

static ProjectionParameterSummary classify_projection_parameters(const RooFitResult *result)
{
	ProjectionParameterSummary summary;
	if (!result)
		return summary;

	const RooArgList &floating = result->floatParsFinal();
	for (int i = 0; i < floating.getSize(); ++i)
	{
		const RooAbsArg *arg = floating.at(i);
		if (!arg)
			continue;
		const std::string name = arg->GetName();
		if (name == "Nsig" || name == "Nbkg")
			summary.common.push_back(name);
		else if (starts_with(name, "signal_mass_") || starts_with(name, "bkg_mass_"))
			summary.massOnly.push_back(name);
		else if (name == "bFraction" || starts_with(name, "npCompositionFrac") || name == "bkgPositiveMassSlope" || starts_with(name, "bkgCompositionFrac") || starts_with(name, "ctau") || starts_with(name, "bkgCtau") ||
						 name.find("lifetime") != std::string::npos || starts_with(name, "NsignalSS") ||
						 name == "Nplb" || starts_with(name, "Nss") || starts_with(name, "Nds") || starts_with(name, "Nflip"))
			summary.lifetimeOnly.push_back(name);
		else
			summary.unknown.push_back(name);
	}
	return summary;
}

static std::string join_parameter_names(const std::vector<std::string> &names)
{
	std::string joined;
	for (const auto &name : names)
	{
		if (!joined.empty())
			joined += ",";
		joined += name;
	}
	return joined;
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

void fit2d_september_tune_v2(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4,
					 bool fixResolutionParams = false, bool drawFromSavedFit = false, bool publish = true,
					 const char *outputDir = "", const char *zoomOutputDir = "",
                     double zoomLow = -2.0, double zoomHigh = 3.0)
{
    // A nonempty zoom directory is strictly a saved-fit, lifetime-PDF-only replay.
    const bool zoomOnly = zoomOutputDir[0];
    if (zoomOnly) {
        drawFromSavedFit = true;
        publish = true;
        if (outputDir[0] || !std::isfinite(zoomLow) || !std::isfinite(zoomHigh) || zoomLow >= zoomHigh)
            throw std::runtime_error("Invalid saved-fit zoom arguments");
    }
	Fit2dScopedMacroTimer timer("fit2d", ptLow, ptHigh, yLow, yHigh);
	// KST - Both regions use main prefit inputs; no tuning result is loaded as a fit seed.
	const double unset = std::numeric_limits<double>::quiet_NaN();
	// KST - Support both configured regions; retain actual bin arguments for cuts and output tags.
	const bool isMid = yLow == 0.f && yHigh == 1.6f;
	const bool isForward = yLow == 1.6f && yHigh == 2.4f;
	const int midEdges[] = {7,8,9,10,12,14,16,20,30};
	const int forwardEdges[] = {1,2,3,4,5,6,7,8,9,10,12,14,20};
	const int *edges = isMid ? midEdges : forwardEdges;
	const int nBins = isMid ? 8 : 12;
	bool supported = false;
	for (int i=0;i<nBins;++i) supported |= ptLow == edges[i] && ptHigh == edges[i+1];
	if (!supported || (!isMid && !isForward))
		throw std::runtime_error("Use one of the eight mid or twelve forward bins");
	// Mid bins >=14 reuse the legacy 14-20 numerical window, not its parameter policy.
	// KST - settingPt* selects policy only; ptLow/ptHigh still select the real data and calibration files.
	const float settingPtLow = isMid && ptLow >= 14 ? 14.f : ptLow;
	const float settingPtHigh = isMid && ptLow >= 14 ? 20.f : ptHigh;
	if (fixResolutionParams)
		throw std::runtime_error("The configured models require fixResolutionParams=false");
    OOFinalFitPolicy policy = oo_final_policy(isMid, ptLow);
    const OORuleOverride windowRule = oo_load_rule_override(Form("y%.2f_%.2f",yLow,yHigh),ptLow,ptHigh);
    // Only the declared September-strength sensitivities are accepted. Means,
    // native covariance correlations, resolution domains and no-scale model stay intact.
    const std::vector<std::string> allowedTuneKeys = {"ct_low", "ct_high", "width_sigma_multiplier",
        "prompt_fraction_cov_scale", "prompt_full_cov_scale", "np_cov_scale", "background_penalty_scale",
        "profile_mass", "profile_bkg", "np_longest_free"};
    for (const auto &kv : windowRule.values) {
        if (std::find(allowedTuneKeys.begin(), allowedTuneKeys.end(), kv.first) == allowedTuneKeys.end())
            throw std::runtime_error("Unsupported September tuning key: " + kv.first);
        if (!std::isfinite(kv.second)) throw std::runtime_error("Nonfinite September tuning value");
        if ((kv.first == "profile_mass" || kv.first == "profile_bkg" || kv.first == "np_longest_free")
            && kv.second != 0. && kv.second != 1.) throw std::runtime_error("Boolean tuning key must be 0 or 1");
    }
    const double widthSigmaMultiplier = windowRule.get("width_sigma_multiplier", 1.);
    const double promptFullCovarianceScale = windowRule.get("prompt_full_cov_scale",100.);
    policy.fractionCovarianceScale = windowRule.get("prompt_fraction_cov_scale", policy.fractionCovarianceScale);
    policy.backgroundPenaltyScale = windowRule.get("background_penalty_scale", policy.backgroundPenaltyScale);
    const double npCovarianceScale = windowRule.get("np_cov_scale", isMid ? 1. : (ptLow<2?175.:ptLow<3?700.:ptLow<10?1.:2.));
    const bool npLongestLifetimeFree = windowRule.flag("np_longest_free", isMid || (isForward && ptLow >= 5));
    for (double v : {widthSigmaMultiplier, promptFullCovarianceScale, policy.fractionCovarianceScale, policy.backgroundPenaltyScale, npCovarianceScale})
        if (!(v > 0.)) throw std::runtime_error("Prior sigma/covariance multipliers must be positive");
    const bool highPt = policy.highPt;
    const bool directHighNP = highPt;
    double npLifetimeLow = highPt ? 0.1 : unset;
    double npLifetimeHigh = highPt ? 1.0 : unset;
    double np2LifetimeLow = unset, np2LifetimeHigh = unset;
    const double npLifetimeInit = highPt ? 0.35 : unset;
    const bool constrainNP = isForward || !highPt;
    const bool constrainBackgroundLifetimes = true;
    const bool constrainPrompt = true;
    const bool constrainPromptFractions = !isMid;
    const bool floatNativePromptMean = policy.floatPromptMean;
    const bool fullPromptCovariance = policy.fullPromptCovariance;
    const bool profileBackground = windowRule.flag("profile_bkg", isMid || highPt);
    const bool profileNPShape = highPt;
    const bool conditionalBackground = false;
    const bool profileMass = windowRule.flag("profile_mass", policy.profileMass);
    policy.profileMass = profileMass;
    const bool transportNP = true;
    const bool recoveryFit = false;
    const bool tuneMidWidths = isMid;
    const double promptGapMin = 1e-4;
    const int widthPolicy = policy.widthPolicy;
    std::cout << "[FIT_POLICY] version=2 group=" << policy.group << std::endl;

	if (publish)
		drawFromSavedFit = true;
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	bool isWeight = false;

	// Keep the decay constants away from the numerically unstable tau -> 0 limit.
	const double signalLifetimeFloor = 2e-3;
	const double defaultBkgLifetimeFloor = 5e-3;

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	const char *DATA_ROOT = "inputs/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	// basic cut
	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)",
													ptLow, ptHigh, yLow, yHigh);
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	// temporal variable before applying cut
	auto *massTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
	auto *timeTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau"));
	if (!massTmp || !timeTmp)
	{
		std::cerr << "ERROR: required variables mass/ctau are missing in dataset." << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// ctau range cut
	// remove 0.1% of from the peak.
	// ------------------------------------------------------------------
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

	// KST - The mapped pT policy reproduces the selected mid windows and the original forward windows.
	std::pair<double, double> ctRange{-6.0, 6.0};
	if (settingPtLow <= 2)
		ctRange = {-20.0, 20.0};
	if (settingPtLow == 5 || settingPtLow == 6 || settingPtLow == 14)
		ctRange = {-0.5, 2.0};
	// Windows selected in the recorded pilot fits.
	if (settingPtLow == 6)
		ctRange = {-0.4, 2.0};
	if (settingPtLow == 7)
		ctRange = {-1.8, 3.0};
	if (settingPtLow == 8 || settingPtLow == 9)
		ctRange = {-0.5, 2.0};
	if (settingPtLow == 10)
		ctRange = {-2.0, 4.0};
	if (settingPtLow == 12)
		ctRange = {-0.5, 4.0};
    if(windowRule.active) {
        if(!windowRule.has("ct_low") || !windowRule.has("ct_high")) throw std::runtime_error("Window rule requires both bounds");
        ctRange={windowRule.get("ct_low",0),windowRule.get("ct_high",0)};
        if(!std::isfinite(ctRange.first)||!std::isfinite(ctRange.second)||ctRange.first>=ctRange.second) throw std::runtime_error("Invalid ctau window");
        std::cout<<"[WINDOW_RULE] "<<windowRule.rule<<" ctau=["<<ctRange.first<<","<<ctRange.second<<"]"<<std::endl;
    }

	// ------------------------------------------------------------------
	// preapre name tags
	// ------------------------------------------------------------------
	auto formatTag = [](double value)
	{ return TString::Format("%.2f", value); };
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figTag = yTag + "_" + ptTag;
	TString figBaseDir = "figs";
	if (publish)
		figBaseDir = "figs_publish";
	if (outputDir[0])
		figBaseDir = TString(outputDir) + "/" + figBaseDir;
	const TString figDir = zoomOnly ? TString(zoomOutputDir)
        : TString::Format("%s/%s/fit2d", figBaseDir.Data(), yTag.Data());
	const TString fitRootDir = (outputDir[0] ? TString(outputDir) + "/" : TString("")) + TString::Format("roots/%s/fit2d", yTag.Data());
	const TString fitRootName = TString::Format("%s/fit2d_result_%s.root", fitRootDir.Data(), figTag.Data());
	if (drawFromSavedFit)
	{
		TFile saved(fitRootName, "READ");
		auto *low = saved.Get<TParameter<double>>("appliedCtauLow");
		auto *high = saved.Get<TParameter<double>>("appliedCtauHigh");
		if (saved.IsZombie() || !low || !high || low->GetVal() != ctRange.first || high->GetVal() != ctRange.second)
			throw std::runtime_error("Saved fit window differs from the current main configuration");
		auto checkFlag = [&](const char *name, int expected)
		{
			if (auto *flag = saved.Get<TParameter<int>>(name))
				if (flag->GetVal() != expected)
					throw std::runtime_error(std::string("Saved fit policy differs from main: ") + name);
		};
		auto *backend = saved.Get<TNamed>("configuredMinimizer");
		if (!backend || std::string(backend->GetTitle()) != "Minuit2")
			throw std::runtime_error("Saved fit minimizer differs from the common backend");
		checkFlag("stableGaussianConvolution", 0);
		checkFlag("bFractionFixedForProfile", 0);
		checkFlag("widthPolicy", widthPolicy);
		checkFlag("promptConstraintsApplied", constrainPrompt);
		checkFlag("promptFractionPriorApplied", constrainPromptFractions && settingPtLow < 5 && !fullPromptCovariance);
		checkFlag("bkgLifetimeConstraintEnabled", constrainBackgroundLifetimes);
		checkFlag("npConstraintPolicy", constrainNP ? 1 : 0);
		checkFlag("floatNativePromptMean", floatNativePromptMean);
		checkFlag("floatNPLifetimeScale", 0);
		checkFlag("conditionalBackgroundPositiveSlope", conditionalBackground);
		checkFlag("backgroundShapeProfiled", profileBackground);
		checkFlag("massShapeProfiled", profileMass);
		checkFlag("massTailParametersProfiled", 0);
        if (auto *tail = saved.Get<TParameter<int>>("npLongestLifetimeFree"))
            if (tail->GetVal()!=npLongestLifetimeFree) throw std::runtime_error("Saved longest NP lifetime policy mismatch");
		if (auto *lo = saved.Get<TParameter<double>>("npFirstLifetimeLower"))
		{
			auto *hi = saved.Get<TParameter<double>>("npFirstLifetimeUpper");
			if (!hi)
				throw std::runtime_error("Incomplete saved NP lifetime bounds");
			npLifetimeLow = lo->GetVal();
			npLifetimeHigh = hi->GetVal();
		}
		if (auto *lo = saved.Get<TParameter<double>>("npSecondLifetimeLower"))
		{
			auto *hi = saved.Get<TParameter<double>>("npSecondLifetimeUpper");
			if (!hi)
				throw std::runtime_error("Incomplete saved second NP lifetime bounds");
			np2LifetimeLow = lo->GetVal();
			np2LifetimeHigh = hi->GetVal();
		}
	}

	// All 12 forward bins use native widths directly, with no additional
	// common width multiplier. This flag only records the disabled policy.
	const bool commonPromptScale = false;
	const bool independentNativeSigmas = widthPolicy == 5;
	// Default Minuit function precision; only the accepted stopping tolerance is set.
	const double fitTolerance = isForward ? 1e-3 : (policy.functionPrecision > 0 ? 1e-1 : 1e-3);

	// ------------------------------------------------------------------
	// apply cut
	// ------------------------------------------------------------------
	TString cutAll = Form(
			"(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
			"(ctau >= %g && ctau <= %g)",
			ptLow, ptHigh, yLow, yHigh, ctRange.first, ctRange.second);

	auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutAll)));
	if (!dataSel || dataSel->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after final selection: " << cutAll << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// connect observable
	// ------------------------------------------------------------------
	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau"));
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_mass.setRange(2.6, 3.5);
	obs_mass.setMin(2.6);
	obs_mass.setMax(3.5);
	obs_mass.setRange("massFit", 2.6, 3.5);
	obs_time.SetTitle("ctau");
	obs_time.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	obs_time.setBins(10000, "cache");
	obs_time.setRange("ctauFull", ctRange.first, ctRange.second);
	int timePlotBins = std::max(400, int(std::ceil((ctRange.second - ctRange.first) / .02)));
	// Preserve the existing forward display; mid retains at least 400 bins.
	if (!isMid && settingPtLow == 14.0f && settingPtHigh == 20.0f)
		timePlotBins = 200;
	if (yLow >= 1.6f && settingPtHigh <= 5.0f)
		timePlotBins = std::max(timePlotBins, 250);

	RooDataSet *data = dataSel.get();

	// ------------------------------------------------------------------
	// Load saved 1D prefit results for the 2D model
	//   massFile : signal and background mass PDFs, component counts, and mass fit result
	//   prFile   : prompt ctau resolution means, widths, and component fractions
	//   npFile   : nonprompt signal ctau lifetimes and fixed component yields
	//   bkgFile  : background ctau PDFs, yields, and their saved fit result
	// ------------------------------------------------------------------

	// Prepare output directories and apply the common plotting style before
	// any canvases are created.
	gSystem->mkdir(figDir, true);
	if (!zoomOnly) gSystem->mkdir(fitRootDir, true);
	const TString diagnosticDir = fitRootDir + "/diagnostics_" + figTag;
	if (!zoomOnly) gSystem->mkdir(diagnosticDir, true);
	auto diagnosticPath = [&](const char *name)
	{ return diagnosticDir + "/" + name; };

	gROOT->Macro(gSystem->Getenv("OO_ROOTLOGON") ? gSystem->Getenv("OO_ROOTLOGON") : "include/rootlogon.C");

	// Saved prefit file paths.  All four files must refer to this same bin tag.
    auto windowPrefitBase=[](const char*env){return gSystem->Getenv(env)?TString(gSystem->Getenv(env))+"/":TString("");};
	const TString massFileName = windowPrefitBase("OO_MASS_PREFIT_DIR") + TString::Format("roots/%s/mass/mass_model_%s.root", yTag.Data(), figTag.Data());
	const TString prFileName = windowPrefitBase("OO_PR_PREFIT_DIR") + TString::Format("roots/%s/ctau_pr/ctau_resolution_%s.root", yTag.Data(), figTag.Data());
	const TString npFileName = windowPrefitBase("OO_NP_PREFIT_DIR") + TString::Format("roots/%s/ctau_np/ctau_np_model_%s.root", yTag.Data(), figTag.Data());
	const TString bkgFileName = windowPrefitBase("OO_BKG_PREFIT_DIR") + TString::Format("roots/%s/ctau_bkg/ctau_bkg_fitresult_%s.root", yTag.Data(), figTag.Data());

	// File-opening helper
	auto openFile = [](const TString &path) -> TFile *
	{
		TFile *f = TFile::Open(path);
		if (!f || f->IsZombie())
		{
			std::cerr << "ERROR: cannot open " << path << std::endl;
			return nullptr;
		}
		return f;
	};

	// Saved scalar metadata readers
	auto readIntParam = [](TFile &f, const char *name, int fallback = -999)
	{
		auto *param = dynamic_cast<TParameter<int> *>(f.Get(name));
		if (param)
			return param->GetVal();
		return fallback;
	};
	auto readDoubleParam = [](TFile &f, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN())
	{
		auto *param = dynamic_cast<TParameter<double> *>(f.Get(name));
		if (param)
			return param->GetVal();
		return fallback;
	};

	// RooFitResult readers.  Look in both the floating and constant parameter
	// lists because a saved prefit can intentionally fix a shape parameter.
	auto readResultValue = [](RooFitResult &fr, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN())
	{
		auto *var = dynamic_cast<RooRealVar *>(fr.floatParsFinal().find(name));
		if (!var)
			var = dynamic_cast<RooRealVar *>(fr.constPars().find(name));
		if (var)
			return var->getVal();
		return fallback;
	};
	auto readResultError = [](RooFitResult &fr, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN())
	{
		auto *var = dynamic_cast<RooRealVar *>(fr.floatParsFinal().find(name));
		if (!var)
			var = dynamic_cast<RooRealVar *>(fr.constPars().find(name));
		if (var)
			return var->getError();
		return fallback;
	};

	// Open every required input
	std::unique_ptr<TFile> massFile(openFile(massFileName));
	std::unique_ptr<TFile> prFile(openFile(prFileName));
	std::unique_ptr<TFile> npFile(openFile(npFileName));
	std::unique_ptr<TFile> bkgFile(openFile(bkgFileName));
	if (!massFile || !prFile || !npFile || !bkgFile)
		return;

    if (!drawFromSavedFit) {
        require_prefit_quality(*massFile, "fit_result");
        require_prefit_quality(*prFile, "timeResult");
        require_prefit_quality(*npFile, "timeResult");
        require_prefit_quality(*bkgFile, "timeResult");
        TFile mcCalibration(Form("roots/%s/mc_mass/mc_mass_model_%s.root", yTag.Data(), figTag.Data()), "READ");
        require_prefit_quality(mcCalibration, "fit_result");
    }

	// Rebuild saved topology before constructing PDFs; calibration counts can differ.
	if (drawFromSavedFit)
	{
		TFile savedCounts(fitRootName, "READ");
		if (savedCounts.IsZombie())
			throw std::runtime_error("Missing final result for replay");
		if (readIntParam(savedCounts, "fixedMassPolicyVersion", -1) != (profileMass ? 0 : 1) || readIntParam(savedCounts, "coreFixedMixturePolicyVersion", -1) != (widthPolicy == 2 ? 1 : 0) || readIntParam(savedCounts, (profileNPShape && !isForward) ? "npIdentifiableShapeMCConstraint" : "jointNativeNPMcConstraint", 0) != (constrainNP ? 1 : 0))
			throw std::runtime_error("Replay requires fixed-mass constant-mean policy metadata");
		if (readDoubleParam(savedCounts, "promptFractionMCErrorScale", -1.0) != oo_prompt_fraction_error_scale(policy, isMid, ptLow) ||
				readDoubleParam(savedCounts, "symmetricBkgDomainPriorSigmas", -1.0) != 3.0 ||
				readIntParam(savedCounts, "fixedPromptFractions", -1) != (settingPtLow >= 5.0f ? 1 : 0))
			throw std::runtime_error("Saved result has a different common constraint policy");
		if (readIntParam(savedCounts, "widthPolicy", -1) != widthPolicy ||
				readIntParam(savedCounts, "shareSignalBackgroundResolution", -1) != 0)
			throw std::runtime_error("Saved resolution policy differs from the common policy");
		if (readIntParam(savedCounts, "nResolutionComponents", -1) != readIntParam(*prFile, "nResolutionComponents", -2) ||
				readIntParam(savedCounts, "nSignalGaussComponents", -1) != readIntParam(*massFile, "nSignalGaussComponents", -2))
			throw std::runtime_error("Saved final topology differs from the native calibration");
		for (auto entry : {std::make_pair("nSignalCBComponents", massFile.get()), std::make_pair("nSignalSSComponents", npFile.get())})
		{
			const int saved = readIntParam(savedCounts, entry.first, -1);
			if (saved >= 0 && saved != readIntParam(*entry.second, entry.first, -1))
				throw std::runtime_error("Replay calibration component count differs from saved final result");
		}
	}

	// Saved fit results used below for parameter values and uncertainties.
	auto *massFitResult = dynamic_cast<RooFitResult *>(massFile->Get("fit_result"));
	auto *npTimeResult = dynamic_cast<RooFitResult *>(npFile->Get("timeResult"));
	auto *bkgTimeResult = dynamic_cast<RooFitResult *>(bkgFile->Get("timeResult"));
	if (!bkgTimeResult)
	{
		std::cerr << "ERROR: timeResult not found in " << bkgFileName << std::endl;
		return;
	}

	// Shared lifetime bounds used while rebuilding the saved background model.
	const double bkgLifetimeFloor = readDoubleParam(*bkgFile, "bkgLifetimeFloor", defaultBkgLifetimeFloor);
	const double lifetimeRangeFloor = 1e-2;
	const double maxStableLifetime = std::max(profileNPShape ? 6.0 : ctRange.second - ctRange.first, 20.0 * lifetimeRangeFloor);

	// ------------------------------------------------------------------
	// build mass model from mass.C
	// ------------------------------------------------------------------
	const int nSignalGaussComponents = std::clamp(readIntParam(*massFile, "nSignalGaussComponents", 0), 0, 2);
	const int nSignalCBComponents = std::clamp(readIntParam(*massFile, "nSignalCBComponents", 0), 0, 2);
	const int nBkgExpComponents = std::clamp(readIntParam(*massFile, "nBkgExpComponents", 0), 0, 1);
	const int nBkgChebyOrder = std::clamp(readIntParam(*massFile, "nBkgChebyOrder", 0), 0, 6);

	if (nSignalGaussComponents + nSignalCBComponents <= 0)
	{
		std::cerr << "ERROR: invalid saved mass signal configuration." << std::endl;
		return;
	}
	if (nBkgExpComponents + nBkgChebyOrder <= 0)
	{
		std::cerr << "ERROR: no background component configured in mass.C output." << std::endl;
		return;
	}
	if (nBkgExpComponents > 0 && nBkgChebyOrder > 0)
	{
		std::cerr << "ERROR: inconsistent background model configuration in mass.C output." << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// Mass signal shape parameters. Widths are expressed through logarithmic
	// offsets so every fitted width remains strictly positive.
	// ------------------------------------------------------------------
	RooRealVar signal_mass_mean("signal_mass_mean", "signal_mass_mean", readDoubleParam(*massFile, "signal_mass_mean", 3.096), 3.05, 3.15);
	RooRealVar signal_mass_sigma("signal_mass_sigma", "signal_mass_sigma", readDoubleParam(*massFile, "signal_mass_sigma", 0.03), 0.001, 0.080);

	const double sigmaFloor = 1e-4;
	const double signal_mass_sigma2_init = std::max(signal_mass_sigma.getVal() + 0.03, sigmaFloor + 1e-4);
	RooRealVar signal_mass_sigma2_log_offset(
			"signal_mass_sigma2_log_offset", "log(signal_mass_sigma2 - floor)",
			std::log(std::max(readDoubleParam(*massFile, "signal_mass_sigma2", signal_mass_sigma2_init) - sigmaFloor, 1e-4)),
			std::log(1e-4),
			std::log(5.0));
	RooFormulaVar signal_mass_sigma2("signal_mass_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_sigma2_log_offset));
	RooFormulaVar signal_mass_sigma_delta2("signal_mass_sigma_delta2", "@0-@1", RooArgList(signal_mass_sigma2, signal_mass_sigma));
	RooRealVar signal_mass_cb_sigma_base("signal_mass_cb_sigma_base", "signal_mass_cb_sigma_base", readDoubleParam(*massFile, "signal_mass_cb_sigma", 0.035), 0.008, 0.080);
	std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_from_gaus;
	std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_delta_from_gaus;
	std::unique_ptr<RooRealVar> signal_mass_cb_sigma_log_offset;
	RooAbsReal *signal_mass_cb_sigma = &signal_mass_cb_sigma_base;
	RooAbsReal *signal_mass_cb_sigma_delta = nullptr;

	// When a Gaussian core is present, parameterize the CB width as another
	// positive log-offset.  Without that core, use the independent base width.
	if (nSignalGaussComponents >= 1)
	{
		const double signal_mass_cb_sigma_init = std::max(signal_mass_sigma.getVal() + 0.005, sigmaFloor + 1e-4);
		signal_mass_cb_sigma_log_offset = std::make_unique<RooRealVar>(
				"signal_mass_cb_sigma_log_offset", "log(signal_mass_cb_sigma - floor)",
				std::log(std::max(signal_mass_cb_sigma_init - sigmaFloor, 1e-4)),
				std::log(1e-4),
				std::log(5.0));
		signal_mass_cb_sigma_from_gaus = std::make_unique<RooFormulaVar>(
				"signal_mass_cb_sigma", Form("%g + exp(@0)", sigmaFloor),
				RooArgList(*signal_mass_cb_sigma_log_offset));
		signal_mass_cb_sigma_delta_from_gaus = std::make_unique<RooFormulaVar>("signal_mass_cb_sigma_delta", "@0-@1", RooArgList(*signal_mass_cb_sigma_from_gaus, signal_mass_sigma));
		signal_mass_cb_sigma = signal_mass_cb_sigma_from_gaus.get();
		signal_mass_cb_sigma_delta = signal_mass_cb_sigma_delta_from_gaus.get();
		signal_mass_cb_sigma_log_offset->setVal(std::log(std::max(readDoubleParam(*massFile, "signal_mass_cb_sigma", signal_mass_cb_sigma->getVal()) - sigmaFloor, 1e-4)));
	}
	else
	{
		signal_mass_cb_sigma_delta = &signal_mass_cb_sigma_base;
	}

	const double signal_mass_cb_sigma2_init = std::max(signal_mass_cb_sigma->getVal() + 0.02, sigmaFloor + 1e-4);
	RooRealVar signal_mass_cb_sigma2_log_offset(
			"signal_mass_cb_sigma2_log_offset", "log(signal_mass_cb_sigma2 - floor)",
			std::log(std::max(signal_mass_cb_sigma2_init - sigmaFloor, 1e-4)),
			std::log(1e-4),
			std::log(5.0));
	RooFormulaVar signal_mass_cb_sigma2("signal_mass_cb_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_cb_sigma2_log_offset));
	RooFormulaVar signal_mass_cb_sigma_delta2("signal_mass_cb_sigma_delta2", "@0-@1", RooArgList(signal_mass_cb_sigma2, *signal_mass_cb_sigma));

	// Crystal-Ball tail parameters and recursive signal-mixture ratios are
	// seeded from mass.C. The ratio transform below enforces physical fractions.
	RooRealVar signal_mass_cb_alpha("signal_mass_cb_alpha", "signal_mass_cb_alpha", readDoubleParam(*massFile, "signal_mass_cb_alpha", 1.5), 0.01, 50.0);
	RooRealVar signal_mass_cb_alpha2("signal_mass_cb_alpha2", "signal_mass_cb_alpha2", readDoubleParam(*massFile, "signal_mass_cb_alpha2", 2.0), 0.01, 50.0);
	RooRealVar signal_mass_cb_n("signal_mass_cb_n", "signal_mass_cb_n", readDoubleParam(*massFile, "signal_mass_cb_n", 3.0), 0.01, 100.0);
	RooRealVar signal_mass_cb_n2("signal_mass_cb_n2", "signal_mass_cb_n2", readDoubleParam(*massFile, "signal_mass_cb_n2", 4.0), 0.01, 100.0);

	RooRealVar signal_mass_frac_ratio1("signal_mass_frac_ratio1", "signal_mass_frac_ratio1", readDoubleParam(*massFile, "signal_mass_frac_ratio1", 1.86), 1e-3, 1e3);
	RooFormulaVar signal_mass_frac1("signal_mass_frac1", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio1));
	RooRealVar signal_mass_frac_ratio2("signal_mass_frac_ratio2", "signal_mass_frac_ratio2", readDoubleParam(*massFile, "signal_mass_frac_ratio2", 0.25), 1e-3, 1e3);
	RooFormulaVar signal_mass_frac2("signal_mass_frac2", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio2));
	RooRealVar signal_mass_frac_ratio3("signal_mass_frac_ratio3", "signal_mass_frac_ratio3", readDoubleParam(*massFile, "signal_mass_frac_ratio3", 0.111111), 1e-3, 1e3);
	RooFormulaVar signal_mass_frac3("signal_mass_frac3", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio3));

	// Reapply saved prefit central values after constructing the variables.
	// This keeps the declared bounds separate from the starting values.
	if (nSignalGaussComponents >= 1)
		signal_mass_sigma.setVal(readDoubleParam(*massFile, "signal_mass_sigma", signal_mass_sigma.getVal()));
	signal_mass_sigma2_log_offset.setVal(std::log(std::max(readDoubleParam(*massFile, "signal_mass_sigma2", signal_mass_sigma2.getVal()) - sigmaFloor, 1e-4)));
	signal_mass_cb_sigma_base.setVal(readDoubleParam(*massFile, "signal_mass_cb_sigma", signal_mass_cb_sigma_base.getVal()));
	signal_mass_cb_sigma2_log_offset.setVal(std::log(std::max(readDoubleParam(*massFile, "signal_mass_cb_sigma2", signal_mass_cb_sigma2.getVal()) - sigmaFloor, 1e-4)));
	signal_mass_cb_alpha.setVal(readDoubleParam(*massFile, "signal_mass_cb_alpha", signal_mass_cb_alpha.getVal()));
	signal_mass_cb_alpha2.setVal(readDoubleParam(*massFile, "signal_mass_cb_alpha2", signal_mass_cb_alpha2.getVal()));
	signal_mass_cb_n.setVal(readDoubleParam(*massFile, "signal_mass_cb_n", signal_mass_cb_n.getVal()));
	signal_mass_cb_n2.setVal(readDoubleParam(*massFile, "signal_mass_cb_n2", signal_mass_cb_n2.getVal()));
	signal_mass_frac_ratio1.setVal(readDoubleParam(*massFile, "signal_mass_frac_ratio1", signal_mass_frac_ratio1.getVal()));
	signal_mass_frac_ratio2.setVal(readDoubleParam(*massFile, "signal_mass_frac_ratio2", signal_mass_frac_ratio2.getVal()));
	signal_mass_frac_ratio3.setVal(readDoubleParam(*massFile, "signal_mass_frac_ratio3", signal_mass_frac_ratio3.getVal()));
	signal_mass_mean.setVal(readDoubleParam(*massFile, "signal_mass_mean", signal_mass_mean.getVal()));

	// declare mass signal PDFs
	RooGaussian signal_mass_gaus("signal_mass_gaus", "signal_mass_gaus", obs_mass, signal_mass_mean, signal_mass_sigma);
	RooGaussian signal_mass_gaus2("signal_mass_gaus2", "signal_mass_gaus2", obs_mass, signal_mass_mean, signal_mass_sigma2);
	std::unique_ptr<RooAbsPdf> signal_mass_cb_owned;
	if (nSignalCBComponents >= 2)
	{
		signal_mass_cb_owned = std::make_unique<RooCrystalBall>(
				"signal_mass_cb", "signal_mass_cb",
				obs_mass, signal_mass_mean,
				*signal_mass_cb_sigma, signal_mass_cb_sigma2,
				signal_mass_cb_alpha, signal_mass_cb_n,
				signal_mass_cb_alpha2, signal_mass_cb_n2);
	}
	else if (nSignalCBComponents >= 1)
	{
		signal_mass_cb_owned = std::make_unique<RooCBShape>(
				"signal_mass_cb", "signal_mass_cb",
				obs_mass, signal_mass_mean, *signal_mass_cb_sigma, signal_mass_cb_alpha, signal_mass_cb_n);
	}
	RooAbsPdf *signal_mass_cb = signal_mass_cb_owned.get();

	// list PDFs and mass fractions
	RooArgList signalMassPdfList;
	if (nSignalGaussComponents >= 1)
		signalMassPdfList.add(signal_mass_gaus);
	if (nSignalCBComponents >= 1)
		signalMassPdfList.add(*signal_mass_cb);
	if (nSignalGaussComponents >= 2)
		signalMassPdfList.add(signal_mass_gaus2);

	RooArgList signalMassFracList;
	if (signalMassPdfList.getSize() >= 2)
		signalMassFracList.add(signal_mass_frac1);
	if (signalMassPdfList.getSize() >= 3)
		signalMassFracList.add(signal_mass_frac2);
	if (signalMassPdfList.getSize() >= 4)
		signalMassFracList.add(signal_mass_frac3);
	std::unique_ptr<RooAbsPdf> signal_mass_owned;

	// combined mass signal model
	RooAbsPdf *signal_mass_pdf = nullptr;
	if (signalMassPdfList.getSize() == 1 && nSignalGaussComponents >= 1)
		signal_mass_pdf = static_cast<RooAbsPdf *>(&signal_mass_gaus);
	if (signalMassPdfList.getSize() == 1 && nSignalGaussComponents == 0)
		signal_mass_pdf = signal_mass_cb;
	if (!signal_mass_pdf)
	{
		signal_mass_owned = std::make_unique<RooAddPdf>("signal_mass", "signal_mass", signalMassPdfList, signalMassFracList, true);
		signal_mass_pdf = signal_mass_owned.get();
	}

	// ------------------------------------------------------------------
	// Mass background shape: either one exponential or a saved Chebychev
	// ------------------------------------------------------------------
	RooRealVar bkg_mass_p1("bkg_mass_p1", "bkg_mass_p1", readDoubleParam(*massFile, "bkg_mass_p1", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p2("bkg_mass_p2", "bkg_mass_p2", readDoubleParam(*massFile, "bkg_mass_p2", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p3("bkg_mass_p3", "bkg_mass_p3", readDoubleParam(*massFile, "bkg_mass_p3", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p4("bkg_mass_p4", "bkg_mass_p4", readDoubleParam(*massFile, "bkg_mass_p4", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p5("bkg_mass_p5", "bkg_mass_p5", readDoubleParam(*massFile, "bkg_mass_p5", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p6("bkg_mass_p6", "bkg_mass_p6", readDoubleParam(*massFile, "bkg_mass_p6", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_lambda("bkg_mass_lambda", "bkg_mass_lambda", readDoubleParam(*massFile, "bkg_mass_lambda", -1.0), -10.0, -1e-4);

	RooArgList chebyCoeffList;
	if (nBkgChebyOrder >= 1)
		chebyCoeffList.add(bkg_mass_p1);
	if (nBkgChebyOrder >= 2)
		chebyCoeffList.add(bkg_mass_p2);
	if (nBkgChebyOrder >= 3)
		chebyCoeffList.add(bkg_mass_p3);
	if (nBkgChebyOrder >= 4)
		chebyCoeffList.add(bkg_mass_p4);
	if (nBkgChebyOrder >= 5)
		chebyCoeffList.add(bkg_mass_p5);
	if (nBkgChebyOrder >= 6)
		chebyCoeffList.add(bkg_mass_p6);
	// declare mass bkg model
	std::unique_ptr<RooAbsPdf> bkg_mass_owned;
	if (nBkgExpComponents == 1)
		bkg_mass_owned = std::make_unique<RooExponential>("bkg_mass", "bkg_mass", obs_mass, bkg_mass_lambda);
	else
		bkg_mass_owned = std::make_unique<RooChebychev>("bkg_mass", "bkg_mass", obs_mass, chebyCoeffList);
	RooAbsPdf *bkg_mass_pdf = bkg_mass_owned.get();

	// ------------------------------------------------------------------
	// ctau prompt resolution from ctau_pr.C
	// ------------------------------------------------------------------
	const int nResolutionCalibrationComponents = std::clamp(readIntParam(*prFile, "nResolutionComponents", 4), 1, 4);
	const int nResolutionComponents = nResolutionCalibrationComponents;
	auto fitRangeAround = [](double center, double hardLow, double hardHigh, double relWidth, double absWidth)
	{
		if (!std::isfinite(center))
			center = 0.5 * (hardLow + hardHigh);
		center = std::clamp(center, hardLow, hardHigh);
		const double width = std::max(absWidth, relWidth * std::abs(center));
		double low = std::max(hardLow, center - width);
		double high = std::min(hardHigh, center + width);
		if (!(low < high))
		{
			low = hardLow;
			high = hardHigh;
		}
		return std::make_pair(low, high);
	};

	// Native mean is the sole shared fit coordinate when released. For pT
	// 2-3 and 12-14 this is equivalent to the former shift parameterization;
	// changing coordinates does not remove the free mean degree of freedom.
	const double nativePromptMean = readDoubleParam(*prFile, "ctauTime1Mean", 0.);
	if (floatNativePromptMean)
	{
		for (int i = 2; i <= nResolutionComponents; ++i)
			if (std::abs(readDoubleParam(*prFile, Form("ctauTime%dMean", i), 0.) - nativePromptMean) > 1e-12)
				throw std::runtime_error("Shared native mean requires common MC component means");
	}
	std::unique_ptr<RooAbsReal> ctauTime1MeanOwned;
	if (fullPromptCovariance)
		ctauTime1MeanOwned = std::make_unique<RooRealVar>("ctauTime1Mean", "ctauTime1Mean", readDoubleParam(*prFile, "ctauTime1Mean", 0.0), -0.01, 0.01);
	else if (floatNativePromptMean)
		ctauTime1MeanOwned = std::make_unique<RooRealVar>("ctauTime1Mean", "shared native prompt mean (mm)", nativePromptMean, nativePromptMean - .05, nativePromptMean + .05);
	else
		ctauTime1MeanOwned = std::make_unique<RooConstVar>("ctauTime1Mean", "ctauTime1Mean", readDoubleParam(*prFile, "ctauTime1Mean", 0.0));
	RooAbsReal &ctauTime1Mean = *ctauTime1MeanOwned;
	const double ctauTime1ScaleVal = readDoubleParam(*prFile, "ctauTime1Scale", 1.0);
	const auto ctauTime1ScaleRange = fitRangeAround(ctauTime1ScaleVal, 0.005, 0.20, 1.0, 0.015);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale",
														std::clamp(ctauTime1ScaleVal, ctauTime1ScaleRange.first, ctauTime1ScaleRange.second),
														ctauTime1ScaleRange.first,
														ctauTime1ScaleRange.second);
	std::unique_ptr<RooAbsReal> ctauTime2MeanOwned;
	if (fullPromptCovariance)
		ctauTime2MeanOwned = std::make_unique<RooFormulaVar>("ctauTime2Mean", "@0", RooArgList(ctauTime1Mean));
	else if (floatNativePromptMean)
		ctauTime2MeanOwned = std::make_unique<RooFormulaVar>("ctauTime2Mean", "@0", RooArgList(ctauTime1Mean));
	else
		ctauTime2MeanOwned = std::make_unique<RooConstVar>("ctauTime2Mean", "ctauTime2Mean", readDoubleParam(*prFile, "ctauTime2Mean", 0.0));
	RooAbsReal &ctauTime2Mean = *ctauTime2MeanOwned;
	const double ctauTime2ScaleVal = readDoubleParam(*prFile, "ctauTime2Scale", 1.0);
	const double ctauTime2DeltaVal = std::max(0.001, ctauTime2ScaleVal - ctauTime1ScaleVal);
	const auto ctauTime2DeltaRange = fitRangeAround(ctauTime2DeltaVal, 0.001, 0.30, 2.0, 0.040);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta",
														std::clamp(ctauTime2DeltaVal, ctauTime2DeltaRange.first, ctauTime2DeltaRange.second),
														ctauTime2DeltaRange.first,
														ctauTime2DeltaRange.second);
	std::unique_ptr<RooAbsReal> ctauTime2ScaleOwned;
	if (independentNativeSigmas)
		ctauTime2ScaleOwned = std::make_unique<RooRealVar>("ctauTime2Scale", "Gaussian sigma 2 (mm)", ctauTime2ScaleVal, .001, 1.);
	else
		ctauTime2ScaleOwned = std::make_unique<RooFormulaVar>("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooAbsReal &ctauTime2Scale = *ctauTime2ScaleOwned;
	std::unique_ptr<RooAbsReal> ctauTime3MeanOwned;
	if (fullPromptCovariance)
		ctauTime3MeanOwned = std::make_unique<RooFormulaVar>("ctauTime3Mean", "@0", RooArgList(ctauTime1Mean));
	else if (floatNativePromptMean)
		ctauTime3MeanOwned = std::make_unique<RooFormulaVar>("ctauTime3Mean", "@0", RooArgList(ctauTime1Mean));
	else
		ctauTime3MeanOwned = std::make_unique<RooConstVar>("ctauTime3Mean", "ctauTime3Mean", readDoubleParam(*prFile, "ctauTime3Mean", 0.0));
	RooAbsReal &ctauTime3Mean = *ctauTime3MeanOwned;
	const double ctauTime3ScaleVal = readDoubleParam(*prFile, "ctauTime3Scale", 1.0);
	const double ctauTime3DeltaVal = std::max(0.001, ctauTime3ScaleVal - ctauTime2ScaleVal);
	const auto ctauTime3DeltaRange = fitRangeAround(ctauTime3DeltaVal, 0.001, 0.60, 2.0, 0.080);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta",
														std::clamp(ctauTime3DeltaVal, ctauTime3DeltaRange.first, ctauTime3DeltaRange.second),
														ctauTime3DeltaRange.first,
														ctauTime3DeltaRange.second);
	std::unique_ptr<RooAbsReal> ctauTime3ScaleOwned;
	if (independentNativeSigmas)
		ctauTime3ScaleOwned = std::make_unique<RooRealVar>("ctauTime3Scale", "Gaussian sigma 3 (mm)", ctauTime3ScaleVal, .001, 1.);
	else
		ctauTime3ScaleOwned = std::make_unique<RooFormulaVar>("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooAbsReal &ctauTime3Scale = *ctauTime3ScaleOwned;
	std::unique_ptr<RooAbsReal> ctauTime4MeanOwned;
	if (fullPromptCovariance)
		ctauTime4MeanOwned = std::make_unique<RooFormulaVar>("ctauTime4Mean", "@0", RooArgList(ctauTime1Mean));
	else if (floatNativePromptMean)
		ctauTime4MeanOwned = std::make_unique<RooFormulaVar>("ctauTime4Mean", "@0", RooArgList(ctauTime1Mean));
	else
		ctauTime4MeanOwned = std::make_unique<RooConstVar>("ctauTime4Mean", "ctauTime4Mean", readDoubleParam(*prFile, "ctauTime4Mean", 0.0));
	RooAbsReal &ctauTime4Mean = *ctauTime4MeanOwned;
	const double ctauTime4ScaleVal = readDoubleParam(*prFile, "ctauTime4Scale", 1.0);
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	const auto ctauTime4DeltaRange = fitRangeAround(ctauTime4DeltaVal, 0.05, 1.20, 1.5, 0.200);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta",
														std::clamp(ctauTime4DeltaVal, ctauTime4DeltaRange.first, ctauTime4DeltaRange.second),
														ctauTime4DeltaRange.first,
														ctauTime4DeltaRange.second);
	std::unique_ptr<RooAbsReal> ctauTime4ScaleOwned;
	if (independentNativeSigmas)
		ctauTime4ScaleOwned = std::make_unique<RooRealVar>("ctauTime4Scale", "Gaussian sigma 4 (mm)", ctauTime4ScaleVal, .001, 1.);
	else
		ctauTime4ScaleOwned = std::make_unique<RooFormulaVar>("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooAbsReal &ctauTime4Scale = *ctauTime4ScaleOwned;

	// init and limits for prompt resolution sigmas
	if (independentNativeSigmas)
	{
		RooRealVar *widths[] = {&ctauTime1Scale, dynamic_cast<RooRealVar *>(&ctauTime2Scale), dynamic_cast<RooRealVar *>(&ctauTime3Scale), dynamic_cast<RooRealVar *>(&ctauTime4Scale)};
		const double settingsPt8[][4] = {
				{0.021701926429079071, 0.0006860952528097039, 0.008844259464591251, 0.035377037858365004},
				{0.040291707987239675, 0.0012738016446594807, 0.016420216005948982, 0.065680864023795929},
				{0.088061697746611503, 0.0027840253249295085, 0.035888081470954153, 0.14355232588381661},
				{0.27280791998720766, 0.0086246822116821387, 0.11117833415605684, 0.44471333662422735}};
		// note - pT2: accepted default-function-precision study, sigma upper
		// bounds +5%; direct starts and limits in mm. No added resolution scale.
		const double settingsPt2[][4] = {
				{0.042646956108431382, 0.00402233622283827, 0.023083212174472189, 0.0969494911327832},
				{0.093792713655341947, 0.0042225305681415082, 0.041772963783148861, 0.17544644788922523},
				{0.26044771950337353, 0.046773623380517271, 0.083888147047522676, 0.35233021759959526},
				{0.50101053042516475, 0.089852919322086638, 0.23396686422312157, 0.98266082973711066}};
		const double (*settings)[4] = settingPtLow == 2 ? settingsPt2 : settingsPt8;
		for (int i = 0; i < 4; ++i)
		{
			widths[i]->setRange(settings[i][2], settings[i][3]);
			widths[i]->setVal(settings[i][0]);
			widths[i]->setError(settings[i][1]);
		}
	}

	// Recursive fractions map the saved mixture coefficients to component weights.
	const double ctauFrac1Val = readDoubleParam(*prFile, "ctauFrac1", 0.5);
	const double ctauFrac2Val = readDoubleParam(*prFile, "ctauFrac2", 0.2);
	const double ctauFrac3Val = readDoubleParam(*prFile, "ctauFrac3", 0.1);
	double compFrac1Val = ctauFrac1Val;
	if (nResolutionComponents == 1)
		compFrac1Val = 1.0;
	double compFrac2Val = 0.0;
	if (nResolutionComponents >= 2)
	{
		compFrac2Val = 1.0 - ctauFrac1Val;
		if (nResolutionComponents != 2)
			compFrac2Val *= ctauFrac2Val;
	}
	double compFrac3Val = 0.0;
	if (nResolutionComponents >= 3)
	{
		compFrac3Val = (1.0 - ctauFrac1Val) * (1.0 - ctauFrac2Val);
		if (nResolutionComponents != 3)
			compFrac3Val *= ctauFrac3Val;
	}
	RooRealVar ctauRecursiveFrac1("ctauRecursiveFrac1", "", ctauFrac1Val, 0., 1.);
	RooRealVar ctauRecursiveFrac2("ctauRecursiveFrac2", "", ctauFrac2Val, 0., 1.);
	RooRealVar ctauRecursiveFrac3("ctauRecursiveFrac3", "", ctauFrac3Val, 0., 1.);
	ctauRecursiveFrac1.setConstant(true);
	ctauRecursiveFrac2.setConstant(true);
	ctauRecursiveFrac3.setConstant(true);
	RooFormulaVar ctauFrac1("ctauFrac1", "@0", RooArgList(ctauRecursiveFrac1));
	RooFormulaVar ctauFrac2("ctauFrac2", "(1-@0)*@1", RooArgList(ctauRecursiveFrac1, ctauRecursiveFrac2));
	RooFormulaVar ctauFrac3("ctauFrac3", "(1-@0)*(1-@1)*@2", RooArgList(ctauRecursiveFrac1, ctauRecursiveFrac2, ctauRecursiveFrac3));

	auto makeResolution = [&](const char *name, RooAbsReal &mean, RooAbsReal &width) -> std::unique_ptr<RooGaussModel>
	{
		return std::make_unique<RooGaussModel>(name, name, obs_time, mean, width);
	};
	// Build the Gaussian components, then use either one component or their sum.
	auto ctauTime1Owned = makeResolution("ctauTime1", ctauTime1Mean, ctauTime1Scale);
	RooGaussModel &ctauTime1 = *ctauTime1Owned;
	auto ctauTime2Owned = makeResolution("ctauTime2", ctauTime2Mean, ctauTime2Scale);
	RooGaussModel &ctauTime2 = *ctauTime2Owned;
	auto ctauTime3Owned = makeResolution("ctauTime3", ctauTime3Mean, ctauTime3Scale);
	RooGaussModel &ctauTime3 = *ctauTime3Owned;
	auto ctauTime4Owned = makeResolution("ctauTime4", ctauTime4Mean, ctauTime4Scale);
	RooGaussModel &ctauTime4 = *ctauTime4Owned;

	// Preserve native PR calibration coordinates, then condition the same physical
	// Gaussian mixture on this target window. NP uses its full-space kernel.
	const auto prTransportCalibration = oo_prompt_transport::calibrationWindow(figTag.Data(), prFileName.Data());
	RooAbsReal *prTransportMeans[] = {&ctauTime1Mean, &ctauTime2Mean, &ctauTime3Mean, &ctauTime4Mean};
	RooAbsReal *prTransportWidths[] = {&ctauTime1Scale, &ctauTime2Scale, &ctauTime3Scale, &ctauTime4Scale};
	RooAbsPdf *prTransportGaussians[] = {&ctauTime1, &ctauTime2, &ctauTime3, &ctauTime4};
	RooAbsReal *prTransportNativeFirst[] = {&ctauFrac1, &ctauFrac2, &ctauFrac3};
	RooArgList prTransportGaussianList, prTransportMeanList, prTransportWidthList, prTransportNativeList;
	for (int i = 0; i < nResolutionComponents; ++i) {
		prTransportGaussianList.add(*prTransportGaussians[i]);
		prTransportMeanList.add(*prTransportMeans[i]);
		prTransportWidthList.add(*prTransportWidths[i]);
		if (i < nResolutionComponents - 1) prTransportNativeList.add(*prTransportNativeFirst[i]);
	}
	oo_prompt_transport::Model promptTransport(obs_time, prTransportGaussianList, prTransportMeanList,
		prTransportWidthList, prTransportNativeList, prTransportCalibration, ctRange.first, ctRange.second);
	RooResolutionModel &promptResolution = *promptTransport.resolution;
	RooAbsPdf *promptTimePdf = promptTransport.prompt.get();

	// ------------------------------------------------------------------
	// ctau prompt resolution - for ctau bkg model (prompt like background)
	// ------------------------------------------------------------------
	const int nBkgResolutionComponents = std::clamp(readIntParam(*bkgFile, "nResolutionComponents", nResolutionComponents), 1, 4);
	double bkgSavedMean = readResultValue(*bkgTimeResult, "ctauTime1Mean", ctauTime1Mean.getVal());
	// Preserve independently fitted Gaussian origins as fixed offsets from mean1.
	// Equal saved means give zero offsets and exactly the original shared model.
	double bkgSavedMeanOffsets[4] = {0., 0., 0., 0.};
	for (int i = 0; i < nBkgResolutionComponents; ++i) {
		const double savedMean = readResultValue(*bkgTimeResult, Form("ctauTime%dMean", i + 1), bkgSavedMean);
		if (!std::isfinite(savedMean)) throw std::runtime_error("Nonfinite saved background Gaussian mean");
		bkgSavedMeanOffsets[i] = savedMean - bkgSavedMean;
	}
	RooRealVar bkgCtauSharedMean("bkgCtauSharedMean", "shared background origin / mm", bkgSavedMean, -.5, .5);
	std::unique_ptr<RooAbsReal> bkgCtauTime1MeanOwned;
	if (profileBackground)
		bkgCtauTime1MeanOwned = std::make_unique<RooFormulaVar>("bkgCtauTime1Mean", "@0", RooArgList(bkgCtauSharedMean));
	else
		bkgCtauTime1MeanOwned = std::make_unique<RooConstVar>("bkgCtauTime1Mean", "bkgCtauTime1Mean", readResultValue(*bkgTimeResult, "ctauTime1Mean", ctauTime1Mean.getVal()));
	RooAbsReal &bkgCtauTime1Mean = *bkgCtauTime1MeanOwned;
	double bkgCtauTime1ScaleVal = readResultValue(*bkgTimeResult, "ctauTime1Scale", ctauTime1ScaleVal);
	const auto bkgCtauTime1ScaleRange = fitRangeAround(bkgCtauTime1ScaleVal, 0.005, 0.80, 0.75, 0.030);
	RooRealVar bkgCtauTime1Scale("bkgCtauTime1Scale", "bkgCtauTime1Scale",
															 std::clamp(bkgCtauTime1ScaleVal, bkgCtauTime1ScaleRange.first, bkgCtauTime1ScaleRange.second),
															 bkgCtauTime1ScaleRange.first,
															 bkgCtauTime1ScaleRange.second);
	std::unique_ptr<RooAbsReal> bkgCtauTime2MeanOwned;
	if (profileBackground)
		bkgCtauTime2MeanOwned = std::make_unique<RooFormulaVar>("bkgCtauTime2Mean", Form("@0+(%.17g)", bkgSavedMeanOffsets[1]), RooArgList(bkgCtauSharedMean));
	else
		bkgCtauTime2MeanOwned = std::make_unique<RooConstVar>("bkgCtauTime2Mean", "bkgCtauTime2Mean", readResultValue(*bkgTimeResult, "ctauTime2Mean", ctauTime2Mean.getVal()));
	RooAbsReal &bkgCtauTime2Mean = *bkgCtauTime2MeanOwned;
	const double bkgCtauTime2ScaleVal = readResultValue(*bkgTimeResult, "ctauTime2Scale", ctauTime2ScaleVal);
	double bkgCtauTime2DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime2Delta", bkgCtauTime2ScaleVal - bkgCtauTime1ScaleVal));
	const auto bkgCtauTime2DeltaRange = fitRangeAround(bkgCtauTime2DeltaVal, 0.001, 1.00, 1.00, 0.060);
	RooRealVar bkgCtauTime2Delta("bkgCtauTime2Delta", "bkgCtauTime2Delta",
															 std::clamp(bkgCtauTime2DeltaVal, bkgCtauTime2DeltaRange.first, bkgCtauTime2DeltaRange.second),
															 bkgCtauTime2DeltaRange.first,
															 bkgCtauTime2DeltaRange.second);
	RooFormulaVar bkgCtauTime2Scale("bkgCtauTime2Scale", "@0+@1", RooArgList(bkgCtauTime1Scale, bkgCtauTime2Delta));
	std::unique_ptr<RooAbsReal> bkgCtauTime3MeanOwned;
	if (profileBackground)
		bkgCtauTime3MeanOwned = std::make_unique<RooFormulaVar>("bkgCtauTime3Mean", Form("@0+(%.17g)", bkgSavedMeanOffsets[2]), RooArgList(bkgCtauSharedMean));
	else
		bkgCtauTime3MeanOwned = std::make_unique<RooConstVar>("bkgCtauTime3Mean", "bkgCtauTime3Mean", readResultValue(*bkgTimeResult, "ctauTime3Mean", ctauTime3Mean.getVal()));
	RooAbsReal &bkgCtauTime3Mean = *bkgCtauTime3MeanOwned;
	const double bkgCtauTime3ScaleVal = readResultValue(*bkgTimeResult, "ctauTime3Scale", ctauTime3ScaleVal);
	double bkgCtauTime3DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime3Delta", bkgCtauTime3ScaleVal - bkgCtauTime2ScaleVal));
	const auto bkgCtauTime3DeltaRange = fitRangeAround(bkgCtauTime3DeltaVal, 0.001, 2.00, 1.00, 0.120);
	RooRealVar bkgCtauTime3Delta("bkgCtauTime3Delta", "bkgCtauTime3Delta",
															 std::clamp(bkgCtauTime3DeltaVal, bkgCtauTime3DeltaRange.first, bkgCtauTime3DeltaRange.second),
															 bkgCtauTime3DeltaRange.first,
															 bkgCtauTime3DeltaRange.second);
	RooFormulaVar bkgCtauTime3Scale("bkgCtauTime3Scale", "@0+@1", RooArgList(bkgCtauTime2Scale, bkgCtauTime3Delta));
	std::unique_ptr<RooAbsReal> bkgCtauTime4MeanOwned;
	if (profileBackground)
		bkgCtauTime4MeanOwned = std::make_unique<RooFormulaVar>("bkgCtauTime4Mean", Form("@0+(%.17g)", bkgSavedMeanOffsets[3]), RooArgList(bkgCtauSharedMean));
	else
		bkgCtauTime4MeanOwned = std::make_unique<RooConstVar>("bkgCtauTime4Mean", "bkgCtauTime4Mean", readResultValue(*bkgTimeResult, "ctauTime4Mean", ctauTime4Mean.getVal()));
	RooAbsReal &bkgCtauTime4Mean = *bkgCtauTime4MeanOwned;
	const double bkgCtauTime4ScaleVal = readResultValue(*bkgTimeResult, "ctauTime4Scale", ctauTime4ScaleVal);
	double bkgCtauTime4DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime4Delta", bkgCtauTime4ScaleVal - bkgCtauTime3ScaleVal));
	const auto bkgCtauTime4DeltaRange = fitRangeAround(bkgCtauTime4DeltaVal, 0.001, 3.00, 1.00, 0.250);
	RooRealVar bkgCtauTime4Delta("bkgCtauTime4Delta", "bkgCtauTime4Delta",
															 std::clamp(bkgCtauTime4DeltaVal, bkgCtauTime4DeltaRange.first, bkgCtauTime4DeltaRange.second),
															 bkgCtauTime4DeltaRange.first,
															 bkgCtauTime4DeltaRange.second);
	RooFormulaVar bkgCtauTime4Scale("bkgCtauTime4Scale", "@0+@1", RooArgList(bkgCtauTime3Scale, bkgCtauTime4Delta));

	// Saved recursive fractions define the Gaussian-mixture weights.
	double bkgSavedPhysical[3] = {readResultValue(*bkgTimeResult, "ctauFrac1", compFrac1Val), readResultValue(*bkgTimeResult, "ctauFrac2", compFrac2Val), readResultValue(*bkgTimeResult, "ctauFrac3", compFrac3Val)};
	double bkgSavedRecursive[3] = {.5, .5, .5}, bkgRemain = 1.;
	if (profileBackground)
		for (int i = 0; i < nBkgResolutionComponents - 1; i++)
		{
			if (bkgSavedPhysical[i] < 0 || bkgSavedPhysical[i] > bkgRemain + 1e-10)
				throw std::runtime_error("Invalid background weights");
			bkgSavedRecursive[i] = std::clamp(bkgSavedPhysical[i] / bkgRemain, 0., 1.);
			bkgRemain -= bkgSavedPhysical[i];
		}
	RooRealVar bkgCtauRecursiveFrac1("bkgCtauRecursiveFrac1", "", bkgSavedRecursive[0], 0., 1.), bkgCtauRecursiveFrac2("bkgCtauRecursiveFrac2", "", bkgSavedRecursive[1], 0., 1.), bkgCtauRecursiveFrac3("bkgCtauRecursiveFrac3", "", bkgSavedRecursive[2], 0., 1.);
	std::unique_ptr<RooAbsReal> bkgCtauFrac1Owned;
	if (profileBackground)
		bkgCtauFrac1Owned = std::make_unique<RooFormulaVar>("bkgCtauFrac1", "@0", RooArgList(bkgCtauRecursiveFrac1));
	else
		bkgCtauFrac1Owned = std::make_unique<RooConstVar>("bkgCtauFrac1", "bkgCtauFrac1", bkgSavedPhysical[0]);
	RooAbsReal &bkgCtauFrac1 = *bkgCtauFrac1Owned;
	std::unique_ptr<RooAbsReal> bkgCtauFrac2Owned;
	if (profileBackground)
		bkgCtauFrac2Owned = std::make_unique<RooFormulaVar>("bkgCtauFrac2", "(1-@0)*@1", RooArgList(bkgCtauRecursiveFrac1, bkgCtauRecursiveFrac2));
	else
		bkgCtauFrac2Owned = std::make_unique<RooConstVar>("bkgCtauFrac2", "bkgCtauFrac2", bkgSavedPhysical[1]);
	RooAbsReal &bkgCtauFrac2 = *bkgCtauFrac2Owned;
	std::unique_ptr<RooAbsReal> bkgCtauFrac3Owned;
	if (profileBackground)
		bkgCtauFrac3Owned = std::make_unique<RooFormulaVar>("bkgCtauFrac3", "(1-@0)*(1-@1)*@2", RooArgList(bkgCtauRecursiveFrac1, bkgCtauRecursiveFrac2, bkgCtauRecursiveFrac3));
	else
		bkgCtauFrac3Owned = std::make_unique<RooConstVar>("bkgCtauFrac3", "bkgCtauFrac3", bkgSavedPhysical[2]);
	RooAbsReal &bkgCtauFrac3 = *bkgCtauFrac3Owned;

	// Gaussian resolution terms for the prompt-like background component.
	auto bkgCtauTime1Owned = makeResolution("bkgCtauTime1", bkgCtauTime1Mean, bkgCtauTime1Scale);
	RooGaussModel &bkgCtauTime1 = *bkgCtauTime1Owned;
	auto bkgCtauTime2Owned = makeResolution("bkgCtauTime2", bkgCtauTime2Mean, bkgCtauTime2Scale);
	RooGaussModel &bkgCtauTime2 = *bkgCtauTime2Owned;
	auto bkgCtauTime3Owned = makeResolution("bkgCtauTime3", bkgCtauTime3Mean, bkgCtauTime3Scale);
	RooGaussModel &bkgCtauTime3 = *bkgCtauTime3Owned;
	auto bkgCtauTime4Owned = makeResolution("bkgCtauTime4", bkgCtauTime4Mean, bkgCtauTime4Scale);
	RooGaussModel &bkgCtauTime4 = *bkgCtauTime4Owned;

	// Assemble the requested number of Gaussian resolution components.
	std::unique_ptr<RooAddModel> bkgPromptResolutionModel;
	RooResolutionModel *bkgPromptResolutionPtr = &bkgCtauTime1;
	if (nBkgResolutionComponents >= 2)
	{
		RooArgList bkgResolutionPdfList;
		RooArgList bkgResolutionFracList;
		bkgResolutionPdfList.add(bkgCtauTime1);
		bkgResolutionPdfList.add(bkgCtauTime2);
		bkgResolutionFracList.add(bkgCtauFrac1);
		if (nBkgResolutionComponents >= 3)
		{
			bkgResolutionPdfList.add(bkgCtauTime3);
			bkgResolutionFracList.add(bkgCtauFrac2);
		}
		if (nBkgResolutionComponents >= 4)
		{
			bkgResolutionPdfList.add(bkgCtauTime4);
			bkgResolutionFracList.add(bkgCtauFrac3);
		}
		bkgPromptResolutionModel = std::make_unique<RooAddModel>("bkgPromptResolutionModel", "bkgPromptResolutionModel", bkgResolutionPdfList, bkgResolutionFracList);
		bkgPromptResolutionPtr = bkgPromptResolutionModel.get();
	}
	// PDFs used by the prompt-like background and its convolved tail components.
	RooResolutionModel &bkgPromptResolution = *bkgPromptResolutionPtr;
	RooAbsPdf *bkgPromptTimePdf = bkgPromptResolutionPtr;

	// Apply the selected native, common-scale, or ordered-width policy.
	RooRealVar *promptWidths[] = {&ctauTime1Scale, independentNativeSigmas ? dynamic_cast<RooRealVar *>(&ctauTime2Scale) : &ctauTime2Delta, independentNativeSigmas ? dynamic_cast<RooRealVar *>(&ctauTime3Scale) : &ctauTime3Delta, independentNativeSigmas ? dynamic_cast<RooRealVar *>(&ctauTime4Scale) : &ctauTime4Delta};
	for (int i = 0; i < 4; ++i)
		promptWidths[i]->setConstant(fullPromptCovariance || (widthPolicy == 1 && i < nResolutionComponents - 2) || (widthPolicy == 2 && i != nResolutionComponents - 1));
	if (profileBackground)
	{
		bkgCtauTime1Scale.setRange(.001, 2.);
		bkgCtauTime1Scale.setConstant(false);
		bkgCtauTime2Delta.setRange(0., 2.);
		bkgCtauTime2Delta.setConstant(false);
		bkgCtauTime3Delta.setRange(0., 2.);
		bkgCtauTime3Delta.setConstant(false);
		bkgCtauTime4Delta.setRange(0., 2.);
		bkgCtauTime4Delta.setConstant(false);
	}
	else
	{
		bkgCtauTime1Scale.setConstant(true);
		bkgCtauTime2Delta.setConstant(true);
		bkgCtauTime3Delta.setConstant(true);
		bkgCtauTime4Delta.setConstant(true);
	}

	// ------------------------------------------------------------------
	// nonprompt signal model from ctau_np.C
	// ------------------------------------------------------------------
	// import NP component count, lifetimes, floors, and yields from the pure-exponential ctau_np prefit.
	// np fit uses simple exp functions. Here we used RooDecay convoluted with Gauss
	const int nSignalSSComponentsSaved = std::clamp(readIntParam(*npFile, "nSignalSSComponents", 1), 1, 3);

	const double signalLifetime1Val = readDoubleParam(*npFile, "signal_lifetime", 0.1);
	const double signalLifetime2Val = readDoubleParam(*npFile, "signal2_lifetime", 0.2);
	const double signalLifetime3Val = readDoubleParam(*npFile, "signal3_lifetime", 0.4);
	const double signalLifetimeFloorSaved = readDoubleParam(*npFile, "signalLifetimeFloor", signalLifetimeFloor);
	const double signalLifetime2FloorSaved = readDoubleParam(*npFile, "signalLifetime2Floor", std::max(0.75 * signalLifetimeFloorSaved, 5e-3));
	const double signalLifetime3FloorSaved = readDoubleParam(*npFile, "signalLifetime3Floor", std::max(1.00 * signalLifetimeFloorSaved, 5e-3));

	const double signalLogLifetimeOffsetVal = readDoubleParam(
			*npFile, "signal_log_lifetime_offset",
			std::log(std::max(signalLifetime1Val - signalLifetimeFloorSaved, 1e-3)));
	const double signal2LogLifetimeOffsetVal = readDoubleParam(
			*npFile, "signal2_log_lifetime_offset",
			std::log(std::max(signalLifetime2Val - signalLifetime2FloorSaved, 1e-3)));
	const double signal3LogLifetimeOffsetVal = readDoubleParam(
			*npFile, "signal3_log_lifetime_offset",
			std::log(std::max(signalLifetime3Val - signalLifetime3FloorSaved, 1e-3)));

	const bool useSignalSS2 = (nSignalSSComponentsSaved >= 2);
	const bool useSignalSS3 = (nSignalSSComponentsSaved >= 3);
	const int nSignalSSComponents = nSignalSSComponentsSaved;

	// Component 1: nominal single-sided lifetime, convolved with the prompt
	// resolution transferred from ctau_pr.C.
	const double signalLifetimeInit = std::max(0.45, 1.5 * signalLifetimeFloorSaved);
	const double signalLifetimeCeil = std::max(signalLifetimeInit * 5.0, maxStableLifetime);
	const double signalLogLifetimeFitErr = npTimeResult
																						 ? std::max(0.05, readResultError(*npTimeResult, "signal_log_lifetime_offset", 0.15))
																						 : 0.15;
	double signalLogLifetimeMin = -4.0;
	double signalLogLifetimeMax = std::log(std::max(signalLifetimeCeil - signalLifetimeFloorSaved, 1e-4));
	RooRealVar signal_log_lifetime_offset(
			"signal_log_lifetime_offset", "log(signal_lifetime - floor)",
			std::clamp(signalLogLifetimeOffsetVal, signalLogLifetimeMin, signalLogLifetimeMax),
			signalLogLifetimeMin, signalLogLifetimeMax);
	RooFormulaVar signal_lifetime(
			"signal_lifetime",
			Form("%g + exp(@0)", signalLifetimeFloorSaved),
			RooArgList(signal_log_lifetime_offset));

	// Component 2: independent floor and log offset, as in ctau_np.C.
	const double signalLifetime2Init = std::max(0.32, 2.0 * signalLifetime2FloorSaved);
	const double signalLifetime2Ceil = std::max(signalLifetime2Init * 5.0, maxStableLifetime);
	RooRealVar signal2_log_lifetime_offset(
			"signal2_log_lifetime_offset", "log(signal2_lifetime - floor)",
			signal2LogLifetimeOffsetVal, -6.0,
			std::log(std::max(signalLifetime2Ceil - signalLifetime2FloorSaved, 1e-4)));
	RooFormulaVar signal2_lifetime(
			"signal2_lifetime",
			Form("%g + exp(@0)", signalLifetime2FloorSaved),
			RooArgList(signal2_log_lifetime_offset));

	// Component 3: independent floor and log offset, as in ctau_np.C.
	const double signalLifetime3Init = std::max(0.40, 2.0 * signalLifetime3FloorSaved);
	const double signalLifetime3Ceil = std::max(signalLifetime3Init * 5.0, maxStableLifetime);
	RooRealVar signal3_log_lifetime_offset(
			"signal3_log_lifetime_offset", "log(signal3_lifetime - floor)",
			signal3LogLifetimeOffsetVal, -6.0,
			std::log(std::max(signalLifetime3Ceil - signalLifetime3FloorSaved, 2e-4)));
	RooFormulaVar signal3_lifetime(
			"signal3_lifetime",
			Form("%g + exp(@0)", signalLifetime3FloorSaved),
			RooArgList(signal3_log_lifetime_offset));

	// Initial flags: the native joint MC constraint later floats all active NP coordinates.
	if (useSignalSS2)
		signal2_log_lifetime_offset.setConstant(true);
	if (useSignalSS3)
		signal3_log_lifetime_offset.setConstant(true);

	// Extended yields imported from the saved ctau_np.C fit.
	const double npYield1 = readDoubleParam(*npFile, "NsignalSS1", 1.0);
	const double npYield2 = readDoubleParam(*npFile, "NsignalSS2", 0.0);
	const double npYield3 = readDoubleParam(*npFile, "NsignalSS3", 0.0);
	RooRealVar NsignalSS1("NsignalSS1", "NsignalSS1", npYield1, 0., 2. * (npYield1 + npYield2 + npYield3));
	RooRealVar NsignalSS2("NsignalSS2", "NsignalSS2", npYield2, 0., 2. * (npYield1 + npYield2 + npYield3));
	RooRealVar NsignalSS3("NsignalSS3", "NsignalSS3", npYield3, 0., 2. * (npYield1 + npYield2 + npYield3));

	// Component PDFs corresponding to the three lifetime parameterizations.
	RooDecay signal_ss1_time("signal_ss1_time", "signal_ss1_time", obs_time, signal_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss2_time("signal_ss2_time", "signal_ss2_time", obs_time, signal2_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss3_time("signal_ss3_time", "signal_ss3_time", obs_time, signal3_lifetime, promptResolution, RooDecay::SingleSided);

	// Transport native truth-window yields to reconstructed-window mixture weights.
	const double npTransportLow = transportNP ? readDoubleParam(*npFile, "npFitCtauLow") : 0.;
	const double npTransportHigh = transportNP ? readDoubleParam(*npFile, "npFitCtauHigh") : 1.;
	if (!(std::isfinite(npTransportLow) && std::isfinite(npTransportHigh) && npTransportLow >= 0 && npTransportHigh > npTransportLow))
		throw std::runtime_error("Missing truth calibration window metadata");
	if (transportNP)
		obs_time.setRange("npTransportTarget", ctRange.first, ctRange.second);

	// Avoid the artificial +/-1e6 range, which overflows RooGaussModel integrals.
	RooRealVar *npNativeYieldVars[] = {&NsignalSS1, &NsignalSS2, &NsignalSS3};
	RooAbsReal *npLifetimeVars[] = {&signal_lifetime, &signal2_lifetime, &signal3_lifetime};
	RooAbsPdf *npDecayVars[] = {&signal_ss1_time, &signal_ss2_time, &signal_ss3_time};
	std::vector<std::unique_ptr<RooAbsReal>> npTransportIntegrals;
	std::vector<std::unique_ptr<RooFormulaVar>> npTransportWeights;
	RooArgList npTransportWeightList;
	if (transportNP)
		for (int i = 0; i < nSignalSSComponents; ++i)
		{
			npTransportIntegrals.emplace_back(npDecayVars[i]->createIntegral(RooArgSet(obs_time), Range("npTransportTarget")));
			RooAbsReal *target = npTransportIntegrals.back().get();
			npTransportIntegrals.push_back(std::make_unique<RooFormulaVar>(
					Form("npTransportFullIntegral%d", i + 1), "2*@0", RooArgList(*npLifetimeVars[i])));
			RooAbsReal *universe = npTransportIntegrals.back().get();
			npTransportWeights.push_back(std::make_unique<RooFormulaVar>(Form("signal_np_transport_q%d", i + 1),
																																	 Form("@0*@1/@2/(exp(-%.17g/@3)-exp(-%.17g/@3))", npTransportLow, npTransportHigh),
																																	 RooArgList(*npNativeYieldVars[i], *target, *universe, *npLifetimeVars[i])));
			npTransportWeightList.add(*npTransportWeights.back());
		}

	// Without an MC prior, only mixture fractions are identifiable: the common
	// normalization of the native MC yields cancels from the NP data PDF.
	RooRealVar npCompositionFrac1("npCompositionFrac1", "NP component 1 fraction", npYield1 / std::max(1e-12, npYield1 + npYield2 + npYield3), 0., 1.);
	RooRealVar npCompositionFrac2("npCompositionFrac2", "NP recursive component 2 fraction", npYield2 / std::max(1e-12, npYield2 + npYield3), 0., 1.);
	npCompositionFrac1.setConstant(constrainNP || nSignalSSComponents < 2);
	npCompositionFrac2.setConstant(constrainNP || nSignalSSComponents < 3);
	// Assemble the requested number of saved nonprompt lifetime components.
	// The joint MC constraint below floats the active native composition.
	std::unique_ptr<RooAddPdf> signal_np_owned;
	RooAbsPdf *signal_np_time = nullptr;
	if (nSignalSSComponents == 1)
	{
		signal_np_time = &signal_ss1_time;
	}
	else if (nSignalSSComponents == 2)
	{
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time", RooArgList(signal_ss1_time, signal_ss2_time), !constrainNP ? RooArgList(npCompositionFrac1) : (transportNP ? npTransportWeightList : RooArgList(NsignalSS1, NsignalSS2)), !constrainNP);
		signal_np_time = signal_np_owned.get();
	}
	else
	{
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time", RooArgList(signal_ss1_time, signal_ss2_time, signal_ss3_time), !constrainNP ? RooArgList(npCompositionFrac1, npCompositionFrac2) : (transportNP ? npTransportWeightList : RooArgList(NsignalSS1, NsignalSS2, NsignalSS3)), !constrainNP);
		signal_np_time = signal_np_owned.get();
	}

	// User-authorized bin-dependent optimizer starts and bounds.
	// These do not change the two-region prompt-width or MC-fraction policy.
	double configuredBSeed = 0.1, bFractionMin = 0.06, bFractionMax = 0.60;
	if (settingPtLow == 2.0f && settingPtHigh == 3.0f)
	{
		configuredBSeed = 0.10;
		bFractionMin = 0.08;
		bFractionMax = 0.16;
	}
	else if (settingPtLow == 8.0f && settingPtHigh == 9.0f)
	{
		configuredBSeed = 0.32;
		bFractionMin = 0.20;
		bFractionMax = 0.35;
	}
	else if (settingPtLow == 9.0f && settingPtHigh == 10.0f)
	{
		configuredBSeed = 0.34;
		bFractionMin = 0.16;
		bFractionMax = 0.40;
	}
	if (!constrainNP || (settingPtLow == 3 || settingPtLow == 12 || settingPtLow == 4 || settingPtLow == 6 || settingPtLow == 10 || settingPtLow == 14) || settingPtLow == 1 || transportNP)
	{
		bFractionMin = 0.;
		bFractionMax = 1.;
	}
	if (fullPromptCovariance)
		configuredBSeed = 0.2;
	if (isMid) { configuredBSeed = 0.2; bFractionMin = 0.; bFractionMax = 1.; }
    if (isForward) { configuredBSeed=.2; bFractionMin=0.; bFractionMax=1.; }
	const double bFractionSeed = configuredBSeed;
	RooRealVar bFraction(
			"bFraction", "bFraction",
			std::clamp(bFractionSeed, bFractionMin, bFractionMax),
			bFractionMin, bFractionMax);

	// Signal ctau PDF used by signal_core in the final 2D model.
	RooAddPdf signal_time("signal_time", "signal_time", RooArgList(*signal_np_time, *promptTimePdf), RooArgList(bFraction), true);

	// ------------------------------------------------------------------
	// background ctau model from ctau_bkg.C
	// combines prompt-like resolution, single-sided, double-sided,
	// and flipped lifetime components.
	// ------------------------------------------------------------------
	// Component configuration
	const int nBkgSSComponents = std::clamp(readIntParam(*bkgFile, "nSignalSSComponents", 0), 0, 3);
	const int nBkgFlipComponents = std::clamp(readIntParam(*bkgFile, "nSignalFlipComponents", 0), 0, 3);
	const int nBkgDSComponents = std::clamp(readIntParam(*bkgFile, "nSignalDSComponents", 0), 0, 3);

	// Saved log-lifetime seeds, with bin-specific fallback initializations.
	double bkgLogLifetimeOffsetVal = readResultValue(*bkgTimeResult, "bkg_log_lifetime_offset", std::log(std::max(0.05, 1.5 * bkgLifetimeFloor) - bkgLifetimeFloor));
	double bkgLogLifetime2OffsetVal = readResultValue(*bkgTimeResult, "bkg_log_lifetime2_offset", std::log(std::max(0.04, 1.5 * std::max(1.5 * bkgLifetimeFloor, 2e-2)) - std::max(1.5 * bkgLifetimeFloor, 2e-2)));
	double bkgLogLifetime3OffsetVal = readResultValue(*bkgTimeResult, "bkg_log_lifetime3_offset", std::log(std::max(0.10, 1.5 * std::max(2.0 * bkgLifetimeFloor, 3e-2)) - std::max(2.0 * bkgLifetimeFloor, 3e-2)));
	double bkgSymLogLifetimeOffsetVal = readResultValue(*bkgTimeResult, "bkg_sym_log_lifetime_offset", std::log(std::max(0.04, 1.4 * 1e-2) - 1e-2));
	double bkgFlipLogLifetimeOffsetVal = readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime_offset", std::log(std::max(0.05, 1.2 * std::max(0.75 * bkgLifetimeFloor, 5e-3)) - std::max(0.75 * bkgLifetimeFloor, 5e-3)));
	double bkgFlipLogLifetime2OffsetVal = readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime2_offset", std::log(std::max(0.08, 1.3 * std::max(1.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 1e-2)) - std::max(1.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 1e-2)));

	// Floating log offsets for positive background lifetimes.
	RooRealVar bkg_log_lifetime_offset("bkg_log_lifetime_offset", "bkg_log_lifetime_offset", bkgLogLifetimeOffsetVal);
	RooRealVar bkg_log_lifetime2_offset("bkg_log_lifetime2_offset", "bkg_log_lifetime2_offset", bkgLogLifetime2OffsetVal);
	RooRealVar bkg_log_lifetime3_offset("bkg_log_lifetime3_offset", "bkg_log_lifetime3_offset", bkgLogLifetime3OffsetVal);
	RooRealVar bkg_sym_log_lifetime_offset("bkg_sym_log_lifetime_offset", "bkg_sym_log_lifetime_offset", bkgSymLogLifetimeOffsetVal);
	RooRealVar bkg_sym_log_lifetime2_offset("bkg_sym_log_lifetime2_offset", "bkg_sym_log_lifetime2_offset", readResultValue(*bkgTimeResult, "bkg_sym_log_lifetime2_offset", std::log(std::max(0.07, 1.4 * 2e-2) - 2e-2)));
	RooRealVar bkg_sym_log_lifetime3_offset("bkg_sym_log_lifetime3_offset", "bkg_sym_log_lifetime3_offset", readResultValue(*bkgTimeResult, "bkg_sym_log_lifetime3_offset", std::log(std::max(0.12, 1.3 * 4e-2) - 4e-2)));
	const double bkgSymLogLifetimeCenter = bkg_sym_log_lifetime_offset.getVal();
	const double bkgSymLogLifetime2Center = bkg_sym_log_lifetime2_offset.getVal();
	const double bkgSymLogLifetime3Center = bkg_sym_log_lifetime3_offset.getVal();
	RooRealVar bkg_flip_log_lifetime_offset("bkg_flip_log_lifetime_offset", "bkg_flip_log_lifetime_offset", bkgFlipLogLifetimeOffsetVal);
	RooRealVar bkg_flip_log_lifetime2_offset("bkg_flip_log_lifetime2_offset", "bkg_flip_log_lifetime2_offset", bkgFlipLogLifetime2OffsetVal);
	RooRealVar bkg_flip_log_lifetime3_offset("bkg_flip_log_lifetime3_offset", "bkg_flip_log_lifetime3_offset", readResultValue(*bkgTimeResult, "bkg_flip_log_lifetime3_offset", std::log(std::max(0.12, 1.2 * std::max(2.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 2e-2)) - std::max(2.5 * std::max(0.75 * bkgLifetimeFloor, 5e-3), 2e-2))));

	// fix or float
	bkg_log_lifetime_offset.setConstant(false);
	bkg_log_lifetime2_offset.setConstant(false);
	bkg_log_lifetime3_offset.setConstant(false);
	bkg_sym_log_lifetime_offset.setConstant(false);
	bkg_sym_log_lifetime2_offset.setConstant(false);
	bkg_sym_log_lifetime3_offset.setConstant(false);
	bkg_flip_log_lifetime_offset.setConstant(false);
	bkg_flip_log_lifetime2_offset.setConstant(false);
	bkg_flip_log_lifetime3_offset.setConstant(false);

	// Physical lifetime floors and the allowed ranges for each PDFs.
	const double bkgLifetime2Floor = readDoubleParam(*bkgFile, "bkgLifetime2Floor", std::max(1.5 * bkgLifetimeFloor, 1e-2));
	const double bkgLifetime3Floor = readDoubleParam(*bkgFile, "bkgLifetime3Floor", std::max(2.0 * bkgLifetimeFloor, 1.5e-2));
	const double bkgSymLifetimeFloor = readDoubleParam(*bkgFile, "bkgSymLifetimeFloor", 5e-3);
	const double bkgSymLifetime2Floor = readDoubleParam(*bkgFile, "bkgSymLifetime2Floor", std::max(2.0 * bkgSymLifetimeFloor, 1e-2));
	const double bkgSymLifetime3Floor = readDoubleParam(*bkgFile, "bkgSymLifetime3Floor", std::max(3.0 * bkgSymLifetimeFloor, 2e-2));
	const double bkgFlipLifetimeFloor = readDoubleParam(*bkgFile, "bkgFlipLifetimeFloor", 2e-3);
	const double bkgFlipLifetime2Floor = readDoubleParam(*bkgFile, "bkgFlipLifetime2Floor", std::max(1.5 * bkgFlipLifetimeFloor, 5e-3));
	const double bkgFlipLifetime3Floor = readDoubleParam(*bkgFile, "bkgFlipLifetime3Floor", std::max(2.5 * bkgFlipLifetimeFloor, 1e-2));
	const double bkgLifetimeCeil = std::max(0.50, maxStableLifetime);
	const double bkgLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetimeCeil = std::max(0.30, 0.5 * maxStableLifetime);
	const double bkgSymLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgSymLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetimeMinGap = std::max(1.00 * bkgSymLifetimeFloor, 5e-3);
	const double bkgSymLifetime2MinGap = std::max(0.50 * bkgSymLifetime2Floor, 6e-3);
	const double bkgSymLifetime3MinGap = std::max(0.35 * bkgSymLifetime3Floor, 8e-3);
	const double bkgSymLogLifetimeWidth = 3.0 * std::max(readResultError(*bkgTimeResult, "bkg_sym_log_lifetime_offset", 0.15), 0.15);
	const double bkgSymLogLifetime2Width = 3.0 * std::max(readResultError(*bkgTimeResult, "bkg_sym_log_lifetime2_offset", 0.20), 0.20);
	const double bkgSymLogLifetime3Width = 3.0 * std::max(readResultError(*bkgTimeResult, "bkg_sym_log_lifetime3_offset", 0.25), 0.25);

	// Helpers to keep log offsets within their physical fit ranges.
	auto clampLogOffsetToRange = [](RooRealVar &offset)
	{
		if (offset.getVal() < offset.getMin())
			offset.setVal(offset.getMin());
		else if (offset.getVal() > offset.getMax())
			offset.setVal(offset.getMax());
	};
	auto setLogLifetimeRange = [&](RooRealVar &offset, double floor, double hardLow, double hardHigh)
	{
		const double lowGap = std::max(hardLow - floor, 1e-4);
		const double highGap = std::max(hardHigh - floor, 2.0 * lowGap);
		offset.setRange(std::log(lowGap), std::log(highGap));
		clampLogOffsetToRange(offset);
	};
	auto setLogLifetimeRangeAround = [&](RooRealVar &offset, double floor, double hardLow, double hardHigh, double center, double width)
	{
		const double lowGap = std::max(hardLow - floor, 1e-4);
		const double highGap = std::max(hardHigh - floor, 2.0 * lowGap);
		const double hardLowOffset = std::log(lowGap);
		const double hardHighOffset = std::log(highGap);
		double low = hardLowOffset;
		double high = hardHighOffset;
		if (std::isfinite(center) && std::isfinite(width) && width > 0.0)
		{
			low = std::max(hardLowOffset, center - width);
			high = std::min(hardHighOffset, center + width);
		}
		if (!(std::isfinite(low) && std::isfinite(high) && low < high))
		{
			low = hardLowOffset;
			high = hardHighOffset;
		}
		offset.setRange(low, high);
		clampLogOffsetToRange(offset);
	};

	if (profileBackground)
	{
		RooRealVar *freeBkgOffsets[] = {&bkg_log_lifetime_offset, &bkg_log_lifetime2_offset, &bkg_log_lifetime3_offset, &bkg_sym_log_lifetime_offset, &bkg_sym_log_lifetime2_offset, &bkg_sym_log_lifetime3_offset, &bkg_flip_log_lifetime_offset, &bkg_flip_log_lifetime2_offset, &bkg_flip_log_lifetime3_offset};
		double freeBkgFloors[] = {bkgLifetimeFloor, bkgLifetime2Floor, bkgLifetime3Floor, bkgSymLifetimeFloor, bkgSymLifetime2Floor, bkgSymLifetime3Floor, bkgFlipLifetimeFloor, bkgFlipLifetime2Floor, bkgFlipLifetime3Floor};
		for (int i = 0; i < 9; i++)
			setLogLifetimeRange(*freeBkgOffsets[i], freeBkgFloors[i], freeBkgFloors[i] + 1e-4, profileNPShape ? 12.0 : std::max(10., 2. * (ctRange.second - ctRange.first)));
	}
	else
	{
		// Apply the ranges before constructing the lifetime formula variables.
		setLogLifetimeRange(bkg_log_lifetime_offset, bkgLifetimeFloor, 0.02, bkgLifetimeCeil);
		setLogLifetimeRange(bkg_log_lifetime2_offset, bkgLifetime2Floor, 0.05, bkgLifetime2Ceil);
		setLogLifetimeRange(bkg_log_lifetime3_offset, bkgLifetime3Floor, 0.10, bkgLifetime3Ceil);
		setLogLifetimeRangeAround(bkg_sym_log_lifetime_offset, bkgSymLifetimeFloor, 0.02,
															bkgSymLifetimeCeil, bkgSymLogLifetimeCenter, bkgSymLogLifetimeWidth);
		setLogLifetimeRangeAround(bkg_sym_log_lifetime2_offset, bkgSymLifetime2Floor, 0.05,
															bkgSymLifetime2Ceil, bkgSymLogLifetime2Center, bkgSymLogLifetime2Width);
		setLogLifetimeRangeAround(bkg_sym_log_lifetime3_offset, bkgSymLifetime3Floor, 0.10,
															bkgSymLifetime3Ceil, bkgSymLogLifetime3Center, bkgSymLogLifetime3Width);
		setLogLifetimeRange(bkg_flip_log_lifetime_offset, bkgFlipLifetimeFloor, 0.02,
												std::max(0.50, maxStableLifetime));
		setLogLifetimeRange(bkg_flip_log_lifetime2_offset, bkgFlipLifetime2Floor, 0.05,
												std::max(0.60, maxStableLifetime));
		setLogLifetimeRange(bkg_flip_log_lifetime3_offset, bkgFlipLifetime3Floor, 0.10,
												std::max(0.80, maxStableLifetime));
	}
	// Matched-window calibrations are empirical fitted PDFs. Preserve their saved
	// active lifetime coordinates and domains; window width must not clip a lifetime.
	RooRealVar *bkgNativeOffsets[] = {&bkg_log_lifetime_offset, &bkg_log_lifetime2_offset, &bkg_log_lifetime3_offset, &bkg_sym_log_lifetime_offset, &bkg_sym_log_lifetime2_offset, &bkg_sym_log_lifetime3_offset, &bkg_flip_log_lifetime_offset, &bkg_flip_log_lifetime2_offset, &bkg_flip_log_lifetime3_offset};
	const bool bkgNativeActive[] = {nBkgSSComponents>=1,nBkgSSComponents>=2,nBkgSSComponents>=3,nBkgDSComponents>=1,nBkgDSComponents>=2,nBkgDSComponents>=3,nBkgFlipComponents>=1,nBkgFlipComponents>=2,nBkgFlipComponents>=3};
	double bkgPreRestoreOffsets[9] = {};
	int bkgNativeRestoreChanges = 0;
	for(int i=0;i<9;++i) {
		if(!bkgNativeActive[i]) continue;
		auto *v=bkgNativeOffsets[i];
		auto *saved=dynamic_cast<const RooRealVar*>(bkgTimeResult->floatParsFinal().find(v->GetName()));
		if(!saved) saved=dynamic_cast<const RooRealVar*>(bkgTimeResult->constPars().find(v->GetName()));
		if(!saved||!std::isfinite(saved->getVal())||!std::isfinite(saved->getMin())||!std::isfinite(saved->getMax())||!(saved->getMin()<saved->getMax())||saved->getVal()<saved->getMin()||saved->getVal()>saved->getMax())
			throw std::runtime_error(std::string("Invalid saved background lifetime domain: ")+v->GetName());
		bkgPreRestoreOffsets[i]=v->getVal();
		if(std::abs(v->getVal()-saved->getVal())>1e-12) ++bkgNativeRestoreChanges;
		v->setRange(saved->getMin(),saved->getMax());
		v->setVal(saved->getVal());
	}
	// Positive lifetime variables used by the convolved background PDFs.
	RooFormulaVar bkg_lifetime("bkg_lifetime", Form("%g + exp(@0)", bkgLifetimeFloor), RooArgList(bkg_log_lifetime_offset));
	RooFormulaVar bkg_lifetime2("bkg_lifetime2", Form("%g + exp(@0)", bkgLifetime2Floor), RooArgList(bkg_log_lifetime2_offset));
	RooFormulaVar bkg_lifetime3("bkg_lifetime3", Form("%g + exp(@0)", bkgLifetime3Floor), RooArgList(bkg_log_lifetime3_offset));
	RooFormulaVar bkg_sym_lifetime("bkg_sym_lifetime", Form("%g + exp(@0)", bkgSymLifetimeFloor), RooArgList(bkg_sym_log_lifetime_offset));
	RooFormulaVar bkg_sym_lifetime2("bkg_sym_lifetime2", Form("%g + exp(@0)", bkgSymLifetime2Floor), RooArgList(bkg_sym_log_lifetime2_offset));
	RooFormulaVar bkg_sym_lifetime3("bkg_sym_lifetime3", Form("%g + exp(@0)", bkgSymLifetime3Floor), RooArgList(bkg_sym_log_lifetime3_offset));
	RooFormulaVar bkg_flip_lifetime("bkg_flip_lifetime", Form("%g + exp(@0)", bkgFlipLifetimeFloor), RooArgList(bkg_flip_log_lifetime_offset));
	RooFormulaVar bkg_flip_lifetime2("bkg_flip_lifetime2", Form("%g + exp(@0)", bkgFlipLifetime2Floor), RooArgList(bkg_flip_log_lifetime2_offset));
	RooFormulaVar bkg_flip_lifetime3("bkg_flip_lifetime3", Form("%g + exp(@0)", bkgFlipLifetime3Floor), RooArgList(bkg_flip_log_lifetime3_offset));

	// Convolve each lifetime PDFs with the background prompt resolution.
	RooDecay bkg_time_ss("bkg_time_ss", "bkg_time_ss", obs_time, bkg_lifetime, bkgPromptResolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss2("bkg_time_ss2", "bkg_time_ss2", obs_time, bkg_lifetime2, bkgPromptResolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss3("bkg_time_ss3", "bkg_time_ss3", obs_time, bkg_lifetime3, bkgPromptResolution, RooDecay::SingleSided);
	RooDecay bkg_time_ds("bkg_time_ds", "bkg_time_ds", obs_time, bkg_sym_lifetime, bkgPromptResolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds2("bkg_time_ds2", "bkg_time_ds2", obs_time, bkg_sym_lifetime2, bkgPromptResolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds3("bkg_time_ds3", "bkg_time_ds3", obs_time, bkg_sym_lifetime3, bkgPromptResolution, RooDecay::DoubleSided);
	RooDecay bkg_time_flip("bkg_time_flip", "bkg_time_flip", obs_time, bkg_flip_lifetime, bkgPromptResolution, RooDecay::Flipped);
	RooDecay bkg_time_flip2("bkg_time_flip2", "bkg_time_flip2", obs_time, bkg_flip_lifetime2, bkgPromptResolution, RooDecay::Flipped);
	RooDecay bkg_time_flip3("bkg_time_flip3", "bkg_time_flip3", obs_time, bkg_flip_lifetime3, bkgPromptResolution, RooDecay::Flipped);

	// keep weak yields away from the zero-seed covariance singularity in this bin.
	double NplbVal = readResultValue(*bkgTimeResult, "Nplb", 1.0);
	double Nss1Val = readResultValue(*bkgTimeResult, "Nss1", 0.0);
	double Nss2Val = readResultValue(*bkgTimeResult, "Nss2", 0.0);
	const double Nss3Val = readResultValue(*bkgTimeResult, "Nss3", 0.0);
	double Nds1Val = readResultValue(*bkgTimeResult, "Nds1", 0.0);
	const double Nds2Val = readResultValue(*bkgTimeResult, "Nds2", 0.0);
	const double Nds3Val = readResultValue(*bkgTimeResult, "Nds3", 0.0);
	double Nflip1Val = readResultValue(*bkgTimeResult, "Nflip1", 0.0);
	double Nflip2Val = readResultValue(*bkgTimeResult, "Nflip2", 0.0);
	const double Nflip3Val = readResultValue(*bkgTimeResult, "Nflip3", 0.0);
	double bkgYieldSum = NplbVal;
	if (nBkgSSComponents >= 1)
		bkgYieldSum += Nss1Val;
	if (nBkgSSComponents >= 2)
		bkgYieldSum += Nss2Val;
	if (nBkgSSComponents >= 3)
		bkgYieldSum += Nss3Val;
	if (nBkgDSComponents >= 1)
		bkgYieldSum += Nds1Val;
	if (nBkgDSComponents >= 2)
		bkgYieldSum += Nds2Val;
	if (nBkgDSComponents >= 3)
		bkgYieldSum += Nds3Val;
	if (nBkgFlipComponents >= 1)
		bkgYieldSum += Nflip1Val;
	if (nBkgFlipComponents >= 2)
		bkgYieldSum += Nflip2Val;
	if (nBkgFlipComponents >= 3)
		bkgYieldSum += Nflip3Val;
	bkgYieldSum = std::max(1e-12, bkgYieldSum);

	// Rebuild the extended background mixture and Gaussian yield constraints.
	std::vector<std::unique_ptr<RooAbsReal>> bkgCoeffVars;
	std::vector<std::string> bkgComponentNames;
	std::vector<double> bkgComponentSeeds;
	RooArgList bkgTimePdfList;
	RooArgList bkgTimeCoeffList;
	RooArgSet bkgCoeffConstraints;
	std::vector<std::unique_ptr<RooConstVar>> bkgCoeffConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> bkgCoeffConstraintPdfs;
	auto addBkgComponent = [&](RooAbsPdf &pdf, const char *name, double yield)
	{
		const double init = std::max(1e-6, yield);
		bkgTimePdfList.add(pdf);
		if (profileBackground)
		{
			bkgComponentNames.emplace_back(name);
			bkgComponentSeeds.push_back(init);
			return;
		}
		const double sigma = std::max({0.25 * init, 0.05 * bkgYieldSum, 1.0});
		const double lower = 1e-6;
		const double upper = std::max(10.0 * bkgYieldSum, 10.0 * init);
		bkgCoeffVars.push_back(std::make_unique<RooRealVar>(name, name, init, lower, upper));
		bkgTimeCoeffList.add(*bkgCoeffVars.back());
		bkgCoeffConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_mean", name), "", init));
		bkgCoeffConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_sigma", name), "", sigma));
		bkgCoeffConstraintPdfs.push_back(std::make_unique<RooGaussian>(
				Form("%s_constraint", name), "", *bkgCoeffVars.back(),
				*bkgCoeffConstraintConsts[bkgCoeffConstraintConsts.size() - 2],
				*bkgCoeffConstraintConsts.back()));
		bkgCoeffConstraints.add(*bkgCoeffConstraintPdfs.back());
	};
	addBkgComponent(*bkgPromptTimePdf, "Nplb", NplbVal);
	if (nBkgSSComponents >= 1)
		addBkgComponent(bkg_time_ss, "Nss1", Nss1Val);
	if (nBkgSSComponents >= 2)
		addBkgComponent(bkg_time_ss2, "Nss2", Nss2Val);
	if (nBkgSSComponents >= 3)
		addBkgComponent(bkg_time_ss3, "Nss3", Nss3Val);
	if (nBkgDSComponents >= 1)
		addBkgComponent(bkg_time_ds, "Nds1", Nds1Val);
	if (nBkgDSComponents >= 2)
		addBkgComponent(bkg_time_ds2, "Nds2", Nds2Val);
	if (nBkgDSComponents >= 3)
		addBkgComponent(bkg_time_ds3, "Nds3", Nds3Val);
	if (nBkgFlipComponents >= 1)
		addBkgComponent(bkg_time_flip, "Nflip1", Nflip1Val);
	if (nBkgFlipComponents >= 2)
		addBkgComponent(bkg_time_flip2, "Nflip2", Nflip2Val);
	if (nBkgFlipComponents >= 3)
		addBkgComponent(bkg_time_flip3, "Nflip3", Nflip3Val);

	std::vector<std::unique_ptr<RooRealVar>> bkgCompositionFractions;
	RooArgList bkgCompositionCoords;
	if (profileBackground)
	{
		double seedRemain = 0.;
		for (double seed : bkgComponentSeeds)
			seedRemain += seed;
		for (int i = 0; i < bkgComponentSeeds.size() - 1; i++)
		{
			double value = bkgComponentSeeds[i] / seedRemain;
			seedRemain -= bkgComponentSeeds[i];
			bkgCompositionFractions.push_back(std::make_unique<RooRealVar>(Form("bkgCompositionFrac%d", i + 1), "", std::clamp(value, 0., 1.), 0., 1.));
			bkgCompositionFractions.back()->setError(.03);
			bkgCompositionCoords.add(*bkgCompositionFractions.back());
		}
		for (int i = 0; i < bkgComponentSeeds.size(); i++)
		{
			std::string expr = "1";
			RooArgList args;
			int n = std::min(i + 1, int(bkgComponentSeeds.size()) - 1);
			for (int j = 0; j < n; j++)
			{
				args.add(*bkgCompositionFractions[j]);
				expr += (j == i ? "*@" : "*(1-@") + std::to_string(j) + (j == i ? "" : ")");
			}
			bkgCoeffVars.push_back(std::make_unique<RooFormulaVar>(bkgComponentNames[i].c_str(), expr.c_str(), args));
			bkgTimeCoeffList.add(*bkgCoeffVars.back());
		}
	}
	// Effective mass dependence of the positive-tail mixture odds, not of tau.
	// Distinct background sources can have different mass spectra; their mixture
	// then need not factorize in mass and ctau (ATLAS arXiv:2309.17177, Sec. 3.3).
	// log[wSS/(1-wSS)] = log[A/C] + k*(m-3.05): a first-order log-odds model.
	// All weights are normalized at each mass; the B_mass marginal is preserved.
	// Historical extension retained as an inactive implementation only;
	// conditionalBackground=false in every production bin. Physical
	// correlations may still be present. Historical references and plots:
	// tuning/background_mass_logic_20260911/README.md. Intercept priors retained.
	RooRealVar bkgPositiveMassSlope("bkgPositiveMassSlope", "positive-tail log-odds mass slope / GeV", settingPtLow == 6 ? 1. : 0., -5., 5.);
	RooAddPdf bkg_time_reference("bkg_time_reference", "zero-slope reference", bkgTimePdfList, bkgTimeCoeffList);
	bkg_time_reference.fixCoefNormalization(RooArgSet(obs_time));
	std::vector<std::unique_ptr<RooFormulaVar>> conditionalBkgWeights, conditionalBkgFractions;
	RooArgList conditionalBkgWeightList, conditionalBkgFractionList;
	if (conditionalBackground)
		for (int i = 0; i < bkgTimeCoeffList.getSize(); ++i)
		{
			bool positive = starts_with(bkgTimeCoeffList.at(i)->GetName(), "Nss");
			conditionalBkgWeights.push_back(std::make_unique<RooFormulaVar>(Form("conditionalBkgWeight%d", i), positive ? "@0*exp(@1*(@2-3.05))" : "@0", (positive ? RooArgList(*bkgTimeCoeffList.at(i), bkgPositiveMassSlope, obs_mass) : RooArgList(*bkgTimeCoeffList.at(i)))));
			conditionalBkgWeightList.add(*conditionalBkgWeights.back());
		}
	std::string denominator = "(";
	for (int i = 0; i < conditionalBkgWeightList.getSize(); ++i)
	{
		if (i)
			denominator += "+";
		denominator += "@" + std::to_string(i);
	}
	denominator += ")";
	for (int i = 0; i < conditionalBkgWeightList.getSize() - 1; ++i)
	{
		conditionalBkgFractions.push_back(std::make_unique<RooFormulaVar>(Form("conditionalBkgFraction%d", i), ("@" + std::to_string(i) + "/" + denominator).c_str(), conditionalBkgWeightList));
		conditionalBkgFractionList.add(*conditionalBkgFractions.back());
	}
	// Background ctau PDF used by bkg_core in the final 2D model.
	std::unique_ptr<RooAddPdf> bkg_time_owned;
	RooAbsPdf *bkg_time_pdf = nullptr;
	if (bkgTimePdfList.getSize() == 1)
	{
		bkg_time_pdf = bkgPromptTimePdf;
	}
	else
	{
		bkg_time_owned = std::make_unique<RooAddPdf>("bkg_time", "bkg_time", bkgTimePdfList, conditionalBackground ? conditionalBkgFractionList : bkgTimeCoeffList);
		if (conditionalBackground)
			bkg_time_owned->fixCoefNormalization(RooArgSet(obs_time));
		bkg_time_pdf = bkg_time_owned.get();
	}

	// ------------------------------------------------------------------
	// Final extended 2D model
	//
	//   model(m, ctau) = Nsig * [S_mass(m) * S_ctau(ctau)]
	//                   + Nbkg * [B_mass(m) * B_ctau(ctau | m)]
	// The conditional term reduces to B_ctau(ctau) when no slope is active.
	// Nsig and Nbkg are the only overall event yields.
	// ------------------------------------------------------------------
	const double nEntries = data->numEntries();
	const double defaultYield = std::max(1.0, 0.5 * nEntries);
	const double maxYield = 2.0 * nEntries;
	const double nSigInit = readDoubleParam(*massFile, "Nsig", defaultYield);
	const double nBkgInit = readDoubleParam(*massFile, "Nbkg", defaultYield);

	// Extended signal and background yields.
	RooRealVar Nsig("Nsig", "Nsig", nSigInit, 0.0, maxYield);
	RooRealVar Nbkg("Nbkg", "Nbkg", nBkgInit, 0.0, maxYield);
	Nsig.setConstant(false);
	Nbkg.setConstant(false);

	// prompt_core and np_core kept for plotting (prompt/nonprompt component)
	RooProdPdf prompt_core(
			"prompt_core", "prompt_core",
			RooArgSet(*signal_mass_pdf, *promptTimePdf));
	RooProdPdf np_core(
			"np_core", "np_core",
			RooArgSet(*signal_mass_pdf, *signal_np_time));

	// Two fitted 2D components: signal (mass x signal ctau) and background.
	RooProdPdf signal_core(
			"signal_core", "signal_core",
			RooArgSet(*signal_mass_pdf, signal_time));
	std::unique_ptr<RooProdPdf> bkgCoreOwned;
	if (conditionalBackground)
		bkgCoreOwned = std::make_unique<RooProdPdf>("bkg_core", "mass PDF times ctau PDF conditional on mass", *bkg_mass_pdf, Conditional(*bkg_time_pdf, RooArgSet(obs_time)));
	else
		bkgCoreOwned = std::make_unique<RooProdPdf>("bkg_core", "bkg_core", RooArgSet(*bkg_mass_pdf, *bkg_time_pdf));
	RooProdPdf &bkg_core = *bkgCoreOwned;

	// Extended sum used by fitTo(): yields multiply the corresponding PDFs.
	RooAddPdf model(
			"model", "model",
			RooArgList(signal_core, bkg_core),
			RooArgList(Nsig, Nbkg));

	auto auditConditionalBackground = [&](const char *phase)
	{
		if (!conditionalBackground)
			return;
		double savedM = obs_mass.getVal(), savedT = obs_time.getVal(), savedK = bkgPositiveMassSlope.getVal();
		RooArgSet timeNorm(obs_time), jointNorm(obs_mass, obs_time), massNorm(obs_mass);
		auto marginal = std::unique_ptr<RooAbsReal>(bkg_core.createIntegral(timeNorm, NormSet(jointNorm)));
		double maxMarginal = 0., maxZero = 0., minValue = 1e300;
		for (double k : {-5., 0., 5.})
		{
			bkgPositiveMassSlope.setVal(k);
			for (int j = 0; j <= 18; j++)
			{
				obs_mass.setVal(2.6 + .9 * j / 18.);
				double expected = bkg_mass_pdf->getVal(&massNorm), actual = marginal->getVal();
				maxMarginal = std::max(maxMarginal, std::abs(actual / expected - 1.));
				for (int q = 0; q <= 200; q++)
				{
					obs_time.setVal(ctRange.first + (ctRange.second - ctRange.first) * q / 200.);
					double v = bkg_core.getVal(&jointNorm);
					if (!std::isfinite(v) || v < 0.)
						throw std::runtime_error("Invalid conditional background density");
					minValue = std::min(minValue, v);
					if (k == 0.)
					{
						double old = bkg_time_reference.getVal(&timeNorm), now = bkg_time_pdf->getVal(&timeNorm);
						maxZero = std::max(maxZero, std::abs(now - old) / std::max(1e-14, std::abs(old)));
					}
				}
			}
		}
		obs_mass.setVal(savedM);
		obs_time.setVal(savedT);
		bkgPositiveMassSlope.setVal(savedK);
		std::ofstream o(diagnosticPath(Form("conditional_audit_%s_pt%.0f_%.0f.json", phase, double(ptLow), double(ptHigh))));
		o << std::setprecision(17) << "{\"max_mass_marginal_relative_error\":" << maxMarginal << ",\"max_zero_slope_relative_error\":" << maxZero << ",\"min_density\":" << minValue << "}\n";
		if (maxMarginal > 1e-7 || maxZero > 1e-7)
			throw std::runtime_error("Conditional background normalization audit failed");
	};
	auditConditionalBackground("prefit");

	// ------------------------------------------------------------------
	// implement Gaussian constraint
	// ------------------------------------------------------------------
	std::vector<std::unique_ptr<RooConstVar>> massConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> massConstraintPdfs;
	RooArgSet massConstraints;
	auto addMassConstraint = [&](const char *baseName, RooRealVar &var, double central, double sigma)
	{
		if (!(std::isfinite(central) && std::isfinite(sigma) && sigma > 0.0))
			return;
		if (std::isfinite(var.getMin()))
			central = std::max(central, var.getMin());
		if (std::isfinite(var.getMax()))
			central = std::min(central, var.getMax());
		var.setVal(central);
		massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_mean", baseName), "", central));
		massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_sigma", baseName), "", sigma));
		massConstraintPdfs.push_back(std::make_unique<RooGaussian>(
				Form("%s_constraint", baseName), "",
				var, *massConstraintConsts[massConstraintConsts.size() - 2], *massConstraintConsts.back()));
		massConstraints.add(*massConstraintPdfs.back());
	};
	if(isForward && ptLow>=10){
        auto *cal=prFile->Get<RooFitResult>("timeResult");
        auto *mean=cal?dynamic_cast<const RooRealVar*>(cal->floatParsFinal().find("ctauTime1Mean")):nullptr;
        if(!mean || mean->getError()<=0)throw std::runtime_error("Missing MC mean error");
        addMassConstraint("promptMeanMC",dynamic_cast<RooRealVar&>(ctauTime1Mean),mean->getVal(),8.*mean->getError());
    }
	auto addPositiveMassConstraint = [&](const char *baseName, RooRealVar &offsetVar, double floor, double central, double sigma)
	{
		if (!(std::isfinite(floor) && std::isfinite(central) && std::isfinite(sigma) && sigma > 0.0))
			return;
		if (central <= floor)
			central = floor + std::max(1e-4, 0.25 * std::max(floor, 1e-4));
		const double low = std::max(floor + 1e-6, central - sigma);
		const double high = central + sigma;
		const double centralOffset = std::log(std::max(central - floor, 1e-6));
		const double lowOffset = std::log(std::max(low - floor, 1e-6));
		const double highOffset = std::log(std::max(high - floor, 1e-6));
		const double sigmaOffset = std::max(1e-3, 0.5 * std::abs(highOffset - lowOffset));
		addMassConstraint(baseName, offsetVar, centralOffset, sigmaOffset);
	};
	auto massErr = [&](const char *name, double fallback)
	{
		double err = readDoubleParam(*massFile, Form("%s_err", name), fallback);
		if (massFitResult)
			err = readResultError(*massFitResult, name, err);
		return err;
	};
	// Nsig and Nbkg are fitted once in the final likelihood; no same-data yield priors.
	// Study alternative: native joint NP MC covariance for every active component.
	std::unique_ptr<RooMultiVarGaussian> npJointConstraint;
	if ((isForward || !profileNPShape) && constrainNP)
	{
		// Propagate the native joint MC covariance of active lifetimes and mixture yields.
		if (!npTimeResult || npTimeResult->status() != 0 || npTimeResult->covQual() != 3)
			throw std::runtime_error("Joint NP constraint requires a valid MC calibration");
		RooRealVar *npLifetimes[] = {&signal_log_lifetime_offset, &signal2_log_lifetime_offset, &signal3_log_lifetime_offset};
		RooRealVar *npYields[] = {&NsignalSS1, &NsignalSS2, &NsignalSS3};
        // Mid: release the MC-longest lifetime, retaining the covariance of other coordinates.
        const double nativeLifetimes[] = {signalLifetime1Val, signalLifetime2Val, signalLifetime3Val};
        const int longestNP = std::max_element(nativeLifetimes, nativeLifetimes+nSignalSSComponents)-nativeLifetimes;
		RooArgList npCoordinates, npCenters;
		std::vector<int> npIndices;
		for (int i = 0; i < nSignalSSComponents; ++i)
		{
			for (auto *v : {npLifetimes[i], npYields[i]})
			{
                if (npLongestLifetimeFree && v == npLifetimes[longestNP]) { v->setConstant(false); continue; }
				const int index = npTimeResult->floatParsFinal().index(v->GetName());
				if (index < 0)
					throw std::runtime_error("Missing active NP covariance coordinate");
				auto *saved = dynamic_cast<const RooRealVar *>(npTimeResult->floatParsFinal().at(index));
				v->setVal(saved->getVal());
				if (policy.calibrationSteps) v->setError(saved->getError());
				v->setConstant(false);
				npCoordinates.add(*v);
				npIndices.push_back(index);
				massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("npJoint_%s_mean", v->GetName()), "", saved->getVal()));
				npCenters.add(*massConstraintConsts.back());
			}
		}
		TMatrixDSym npJointCovariance(npIndices.size());
		for (int i = 0; i < int(npIndices.size()); ++i)
			for (int j = 0; j < int(npIndices.size()); ++j)
				npJointCovariance(i, j) = npCovarianceScale * npTimeResult->covarianceMatrix()(npIndices[i], npIndices[j]);
		TDecompChol npCovarianceCheck(npJointCovariance);
		if (!npCovarianceCheck.Decompose())
			throw std::runtime_error("NP MC covariance is not positive definite");
		npJointConstraint = std::make_unique<RooMultiVarGaussian>("npJointConstraint", "joint NP MC constraint", npCoordinates, npCenters, npJointCovariance);
		massConstraints.add(*npJointConstraint);
	}

	if (!constrainNP)
	{
		NsignalSS1.setConstant(true);
		NsignalSS2.setConstant(true);
		NsignalSS3.setConstant(true);
		signal_log_lifetime_offset.setConstant(false);
		signal2_log_lifetime_offset.setConstant(nSignalSSComponents < 2);
		signal3_log_lifetime_offset.setConstant(nSignalSSComponents < 3);
		if (npJointConstraint || massConstraints.find("npJointConstraint"))
			throw std::runtime_error("NP constraint remained enabled in unconstrained profile");
	}
	// Alternative common low-pT policy: propagate native MC fraction covariance.
	std::unique_ptr<RooMultiVarGaussian> promptFractionConstraint;
	if (!fullPromptCovariance && settingPtLow < 5.0f && nResolutionComponents > 1)
	{
		auto *cal = prFile->Get<RooFitResult>("timeResult");
		if (!cal || cal->status() != 0 || cal->covQual() != 3)
			throw std::runtime_error("Invalid prompt fraction MC covariance");
		RooRealVar *vars[] = {&ctauRecursiveFrac1, &ctauRecursiveFrac2, &ctauRecursiveFrac3};
		RooArgList coordinates, centers;
		std::vector<int> indices;
		for (int i = 0; i < nResolutionComponents - 1; ++i)
		{
			int index = cal->floatParsFinal().index(Form("ctauFrac%dVar", i + 1));
			if (index < 0)
				throw std::runtime_error("Missing native MC fraction coordinate");
			auto *v = dynamic_cast<const RooRealVar *>(cal->floatParsFinal().at(index));
			vars[i]->setVal(v->getVal());
			vars[i]->setConstant(false);
			coordinates.add(*vars[i]);
			indices.push_back(index);
			massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("promptFractionMC%d_mean", i + 1), "", v->getVal()));
			centers.add(*massConstraintConsts.back());
		}
		TMatrixDSym covariance(indices.size());
		for (int i = 0; i < int(indices.size()); ++i)
			for (int j = 0; j < int(indices.size()); ++j)
				covariance(i, j) = policy.fractionCovarianceScale * cal->covarianceMatrix()(indices[i], indices[j]);
		TDecompChol check(covariance);
		if (!check.Decompose())
			throw std::runtime_error("Non-positive prompt fraction covariance");
		promptFractionConstraint = std::make_unique<RooMultiVarGaussian>("promptFractionConstraint", "native MC fraction covariance", coordinates, centers, covariance);
		if (constrainPromptFractions)
			massConstraints.add(*promptFractionConstraint);
	}

	auto bkgTimeErr = [&](const char *name, double fallback)
	{
		double err = fallback;
		if (bkgTimeResult)
			err = readResultError(*bkgTimeResult, name, err);
		return err;
	};

	// Standard mass policy: keep every saved signal/background shape parameter fixed.
	signal_mass_mean.setConstant(true);
	signal_mass_sigma.setConstant(true);
	signal_mass_sigma2_log_offset.setConstant(true);
	signal_mass_cb_sigma2_log_offset.setConstant(true);
	signal_mass_cb_alpha.setConstant(true);
	signal_mass_cb_n.setConstant(true);
	signal_mass_frac_ratio1.setConstant(true);
	signal_mass_frac_ratio2.setConstant(true);
	signal_mass_frac_ratio3.setConstant(true);
	if (signal_mass_cb_sigma_log_offset)
		signal_mass_cb_sigma_log_offset->setConstant(true);
	else
		signal_mass_cb_sigma_base.setConstant(true);
	if (nSignalCBComponents >= 2)
	{
		signal_mass_cb_alpha2.setConstant(true);
		signal_mass_cb_n2.setConstant(true);
	}

	// Same Gaussian error prescription; apply only to floating width coordinates.
	const double absoluteWidthPenalty[] = {0.010, 0.020, 0.040, 0.080};
	const double absoluteWidthCenters[] = {ctauTime1ScaleVal, ctauTime2DeltaVal, ctauTime3DeltaVal, ctauTime4DeltaVal};
	for (int i = 0; i < nResolutionComponents; ++i)
	{
		if (constrainPrompt && !fullPromptCovariance && !promptWidths[i]->isConstant())
			addMassConstraint(promptWidths[i]->GetName(), *promptWidths[i],
												absoluteWidthCenters[i], widthSigmaMultiplier * (isForward ? oo_integrated_width_sigma(ptLow,i,nResolutionComponents) : absoluteWidthPenalty[i] * oo_prompt_penalty_scale(isMid, ptLow, i)));
	}
	// Independent prompt MC: one native full covariance for mean, widths and recursive fractions.
	std::unique_ptr<RooMultiVarGaussian> promptFullConstraint;
	if (fullPromptCovariance)
	{
		auto *cal = prFile->Get<RooFitResult>("timeResult");
		if (!cal || cal->status() != 0 || cal->covQual() != 3)
			throw std::runtime_error("Invalid prompt full MC calibration");
		RooRealVar *fractionVars[] = {&ctauRecursiveFrac1, &ctauRecursiveFrac2, &ctauRecursiveFrac3};
		RooArgList coordinates, centers;
		std::vector<int> indices;
		auto addCoordinate = [&](RooRealVar &variable, const char *nativeName)
		{
			int index = cal->floatParsFinal().index(nativeName);
			if (index < 0)
				throw std::runtime_error(std::string("Missing prompt covariance coordinate: ") + nativeName);
			auto *v = dynamic_cast<const RooRealVar *>(cal->floatParsFinal().at(index));
			if (!v || !std::isfinite(v->getError()) || v->getError() <= 0)
				throw std::runtime_error("Invalid native prompt parameter");
			variable.setRange(v->getMin(), v->getMax());
			variable.setVal(v->getVal());
			variable.setConstant(false);
			coordinates.add(variable);
			indices.push_back(index);
			massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("promptFull_%s_center", nativeName), "", v->getVal()));
			centers.add(*massConstraintConsts.back());
			std::cout << "[PromptCovMapping] " << nativeName << " -> " << variable.GetName() << " index=" << index << std::endl;
		};
		addCoordinate(dynamic_cast<RooRealVar &>(ctauTime1Mean), "ctauTime1Mean");
		for (int i = 0; i < nResolutionComponents; ++i)
			addCoordinate(*promptWidths[i], promptWidths[i]->GetName());
		for (int i = 0; i < nResolutionComponents - 1; ++i)
			addCoordinate(*fractionVars[i], Form("ctauFrac%dVar", i + 1));
		if (int(indices.size()) != cal->floatParsFinal().getSize())
			throw std::runtime_error("Prompt covariance contains unmapped floating coordinates");
		TMatrixDSym covariance(indices.size());
		for (int i = 0; i < int(indices.size()); ++i)
			for (int j = 0; j < int(indices.size()); ++j)
				covariance(i, j) = promptFullCovarianceScale * cal->covarianceMatrix()(indices[i], indices[j]);
		TDecompChol check(covariance);
		if (!check.Decompose())
			throw std::runtime_error("Prompt full covariance is not positive definite");
		promptFullConstraint = std::make_unique<RooMultiVarGaussian>("promptFullConstraint", "native full prompt MC covariance", coordinates, centers, covariance);
		massConstraints.add(*promptFullConstraint);
	}
	if (!profileBackground)
	{
		if (constrainBackgroundLifetimes)
		{
			auto addBkgLogOffsetConstraint = [&](const char *name, RooRealVar &var, double fallbackSigma)
			{
				if (!bkgTimeResult->floatParsFinal().find(name))
					throw std::runtime_error(std::string("Missing active background calibration: ") + name);
				double constraintCenter = readResultValue(*bkgTimeResult, name, var.getVal());
				addMassConstraint(name, var, constraintCenter, std::max(fallbackSigma, bkgTimeErr(name, fallbackSigma)) * policy.backgroundPenaltyScale);
			};
			if (nBkgSSComponents >= 1) addBkgLogOffsetConstraint("bkg_log_lifetime_offset", bkg_log_lifetime_offset, 0.12);
			if (nBkgSSComponents >= 2) addBkgLogOffsetConstraint("bkg_log_lifetime2_offset", bkg_log_lifetime2_offset, 0.16);
			if (nBkgSSComponents >= 3) addBkgLogOffsetConstraint("bkg_log_lifetime3_offset", bkg_log_lifetime3_offset, 0.20);
			if (nBkgDSComponents >= 1) addBkgLogOffsetConstraint("bkg_sym_log_lifetime_offset", bkg_sym_log_lifetime_offset, 0.15);
			if (nBkgDSComponents >= 2) addBkgLogOffsetConstraint("bkg_sym_log_lifetime2_offset", bkg_sym_log_lifetime2_offset, 0.20);
			if (nBkgDSComponents >= 3) addBkgLogOffsetConstraint("bkg_sym_log_lifetime3_offset", bkg_sym_log_lifetime3_offset, 0.25);
			if (nBkgFlipComponents >= 1) addBkgLogOffsetConstraint("bkg_flip_log_lifetime_offset", bkg_flip_log_lifetime_offset, 0.12);
			if (nBkgFlipComponents >= 2) addBkgLogOffsetConstraint("bkg_flip_log_lifetime2_offset", bkg_flip_log_lifetime2_offset, 0.16);
			if (nBkgFlipComponents >= 3) addBkgLogOffsetConstraint("bkg_flip_log_lifetime3_offset", bkg_flip_log_lifetime3_offset, 0.20);
		}
		for (auto *arg : bkgCoeffConstraints)
			massConstraints.add(*arg);
	}
	if (nBkgExpComponents == 1)
	{
		bkg_mass_lambda.setConstant(true);
	}
	else
	{
		bkg_mass_p1.setConstant(true);
		bkg_mass_p2.setConstant(true);
		bkg_mass_p3.setConstant(true);
		bkg_mass_p4.setConstant(true);
		bkg_mass_p5.setConstant(true);
		bkg_mass_p6.setConstant(true);
	}

	if (profileMass)
	{

		// Conventional joint mass/lifetime profile: float mass mean, widths,
		// mixture and background shape; CB tail parameters remain fixed.
		// Remove same-data mass-yield priors to avoid counting mass information twice.
		std::unique_ptr<RooArgSet> massShapeParameters(model.getParameters(RooArgSet(obs_mass, obs_time)));
		for (auto *arg : *massShapeParameters)
		{
			auto *v = dynamic_cast<RooRealVar *>(arg);
			if (!v)
				continue;
			std::string n = v->GetName();
			if ((starts_with(n, "signal_mass_") || starts_with(n, "bkg_mass_")) &&
					(n.find("alpha") == std::string::npos && n != "signal_mass_cb_n" && n != "signal_mass_cb_n2"))
				v->setConstant(false);
		}

		for (const char *name : {"Nsig_mass_constraint", "Nbkg_mass_constraint"})
		{
			auto *c = massConstraints.find(name);
			if (c)
				massConstraints.remove(*c);
		}
		bFraction.setRange(0., 1.);
		bFraction.setVal(0.2);
	}
    // Mass-fixing study: preserve the profiled control's B domain and seed.
    // These numerical settings are independent of the mass-shape treatment.
    if (!profileMass && ((isForward && ptLow >= 8 && ptLow < 10) || (isMid && ptLow < 12))) {
        bFraction.setRange(0., 1.);
        bFraction.setVal(0.2);
    }
	if (promptGapMin >= 0)
		for (int i = 1; i < nResolutionComponents; ++i)
			if (!promptWidths[i]->isConstant())
			{
				if (promptGapMin >= promptWidths[i]->getMax())
					throw std::runtime_error("Prompt gap lower bound exceeds its domain");
				promptWidths[i]->setMin(promptGapMin);
			}
	if (profileBackground)
	{
		auto parameters = std::unique_ptr<RooArgSet>(model.getParameters(RooArgSet(obs_mass, obs_time)));
		for (auto *arg : *parameters)
			if (auto *variable = dynamic_cast<RooRealVar *>(arg))
				if (starts_with(variable->GetName(), "bkgCtau"))
					variable->setConstant(true);
	}

    {
        std::ofstream constants(diagnosticPath("SEPTEMBER_TUNING_PRIOR_CONSTANTS.csv"));
        constants << std::setprecision(17) << "name,value\n";
        for (auto &v : massConstraintConsts) constants << v->GetName() << ',' << v->getVal() << '\n';
        for (auto &v : bkgCoeffConstraintConsts) constants << v->GetName() << ',' << v->getVal() << '\n';
    }
	// ------------------------------------------------------------------
	// fit
	// ------------------------------------------------------------------
	// Optimizer controls only: apply after MC constraints have seeded their
	// coordinates; neither their centers nor covariance matrices are changed.
	if (!std::isnan(npLifetimeLow))
	{
		if (!(npLifetimeLow > signalLifetimeFloorSaved && npLifetimeHigh > npLifetimeLow))
			throw std::runtime_error("NP lifetime bounds must be increasing and above the native floor");
		signal_log_lifetime_offset.setRange(std::log(npLifetimeLow - signalLifetimeFloorSaved), std::log(npLifetimeHigh - signalLifetimeFloorSaved));
		signal_log_lifetime_offset.setVal(std::clamp(signal_log_lifetime_offset.getVal(), signal_log_lifetime_offset.getMin() + 1e-8, signal_log_lifetime_offset.getMax() - 1e-8));
	}
	if (!std::isnan(npLifetimeInit))
	{
		if (!(npLifetimeInit > signalLifetimeFloorSaved))
			throw std::runtime_error("NP lifetime start must exceed its floor");
		double start = std::log(npLifetimeInit - signalLifetimeFloorSaved);
		if (!(start > signal_log_lifetime_offset.getMin() && start < signal_log_lifetime_offset.getMax()))
			throw std::runtime_error("NP lifetime start must lie inside its bounds");
		signal_log_lifetime_offset.setVal(start);
	}
	if (!std::isnan(np2LifetimeLow))
	{
		if (!(np2LifetimeLow > signalLifetime2FloorSaved && np2LifetimeHigh > np2LifetimeLow))
			throw std::runtime_error("Invalid second NP lifetime bounds");
		signal2_log_lifetime_offset.setRange(std::log(np2LifetimeLow - signalLifetime2FloorSaved), std::log(np2LifetimeHigh - signalLifetime2FloorSaved));
		signal2_log_lifetime_offset.setVal(std::clamp(signal2_log_lifetime_offset.getVal(), signal2_log_lifetime_offset.getMin() + 1e-8, signal2_log_lifetime_offset.getMax() - 1e-8));
	}

	// initial values for pT2/pT8 native-sigma tests from their supplied
	// physical fit coordinates. Starts only: every listed coordinate remains free.
	if (false && (independentNativeSigmas || (isForward && ptLow == 2)) && !drawFromSavedFit)
	{
		struct Start
		{
			const char *name;
			double value, error;
		};
		const Start startsPt8[] = {
				{"Nbkg", 3469.2833327580274, 127.60448224582501},
				{"Nds1", 238.81332919685883, 56.647000047892476},
				{"Nplb", 1808.2760151834948, 361.43597696497932},
				{"Nsig", 2665.5070500666393, 124.41748736866361},
				{"NsignalSS1", 89438.103408266237, 1075.5023395481694},
				{"NsignalSS2", 3295.2151469520177, 1034.5408093478322},
				{"Nss1", 119.65222334047809, 55.581442921048072},
				{"Nss2", 240.53835960417862, 62.160259725216477},
				{"bFraction", 0.25409841819534029, 0.012197602464351823},
				{"bkg_log_lifetime2_offset", -1.0610115769530539, 0.13744603624998486},
				{"bkg_log_lifetime_offset", -3.3513073079619811, 1.2443827168250259},
				{"bkg_mass_p1", -0.37041125623390642, 0.031279260348821031},
				{"bkg_mass_p2", 0.011021407885570978, 0.047259851097424595},
				{"bkg_sym_log_lifetime_offset", -1.6300310584836539, 0.25062321197227599},
				{"ctauTime1Scale", 0.021701926429079071, 0.0006860952528097039},
				{"ctauTime2Scale", 0.040291707987239675, 0.0012738016446594807},
				{"ctauTime3Scale", 0.088061697746611503, 0.0027840253249295085},
				{"ctauTime4Scale", 0.27280791998720766, 0.0086246822116821387},
				{"signal2_log_lifetime_offset", -2.0406727644135789, 0.19059564798342499},
				{"signal_log_lifetime_offset", -0.9666765593134804, 0.0066261320660649381},
				{"signal_mass_cb_sigma2_log_offset", -2.8809885413516101, 0.32711549089114889},
				{"signal_mass_cb_sigma_log_offset", -2.7799211604072123, 0.19332988114714689},
				{"signal_mass_frac_ratio1", 3.8972344536759094, 6.3203812906595891},
				{"signal_mass_mean", 3.0922649536781308, 0.0018911976609792269},
				{"signal_mass_sigma", 0.042599332028115833, 0.0036316304612638191},
				{"bkg_log_lifetime3_offset", -2.5532945886480176, 0.19994119838789759},
				{"bkg_sym_log_lifetime2_offset", -2.9957322735555527, 0.19631682761423064},
				{"bkg_sym_log_lifetime3_offset", -2.5257286443102083, 0.24539603451778791},
				{"bkg_flip_log_lifetime_offset", -4.4504146879311941, 0.11998868022388054},
				{"bkg_flip_log_lifetime2_offset", -3.6394844453378266, 0.15997289707827589},
				{"bkg_flip_log_lifetime3_offset", -3.534123910913308, 0.19994650920097734}};
		const Start startsPt2[] = {
				{"Nbkg", 291804.02631470823, 470.48568375426112},
				{"Nds1", 20993.525005207131, 4615.8936089517974},
				{"Nplb", 181283.27181072527, 39773.759251721742},
				{"Nsig", 24082.348938888696, 266.6998976591567},
				{"NsignalSS1", 135758.15583856832, 4126.167919952044},
				{"NsignalSS2", 75975.818130692875, 4118.7943253232297},
				{"Nss1", 1630.4229307477613, 477.82887305350698},
				{"Nss2", 943.35771198181203, 234.85469298311529},
				{"bFraction", 0.1225746456447807, 0.0085457011174172037},
				{"bkg_log_lifetime2_offset", 0.60722923467094658, 0.064884961405173303},
				{"bkg_log_lifetime_offset", -1.6358618047203206, 0.2863714663212209},
				{"bkg_sym_log_lifetime_offset", -0.72270085461756262, 0.0175006609615867},
				{"ctauTime1Mean", nativePromptMean + 0.00091044056546419195, 0.0012574069482733066},
				{"ctauRecursiveFrac1", 0.2303952299229097, 0.049959860045862808},
				{"ctauRecursiveFrac2", 0.84690155582736049, 0.026671575720543406},
				{"ctauRecursiveFrac3", 0.65123807329683059, 0.18977499209610138},
				{"ctauTime1Scale", 0.042646956108431382, 0.00402233622283827},
				{"ctauTime2Scale", 0.093792713655341947, 0.0042225305681415082},
				{"ctauTime3Scale", 0.26044771950337353, 0.046773623380517271},
				{"ctauTime4Scale", 0.50101053042516475, 0.089852919322086638},
				{"signal2_log_lifetime_offset", -1.4306098952324557, 0.026663249925991339},
				{"signal_log_lifetime_offset", -0.65723425256440515, 0.010345629046771943}};
		const Start *starts = settingPtLow == 2 ? startsPt2 : startsPt8;
		const size_t startCount = settingPtLow == 2 ? sizeof(startsPt2) / sizeof(Start) : sizeof(startsPt8) / sizeof(Start);
		auto parameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
		for (size_t index = 0; index < startCount; ++index)
		{
			const auto &start = starts[index];
			auto *target = dynamic_cast<RooRealVar *>(parameters->find(start.name));
			if (!target || target->isConstant())
				continue;
			// KST - Only mid 8-9 uses native NP-MC starts to avoid underflow of its new joint constraint.
			if (isMid && settingPtLow == 8 && (starts_with(start.name, "NsignalSS") ||
			    (!starts_with(start.name, "bkg") && std::string(start.name).find("log_lifetime_offset") != std::string::npos)))
			{
				// Forward NP-MC yields/lifetimes are many mid-MC sigma from the new
				// constraint center and make its PDF underflow. Start at that center.
				std::cout << "[MID_NP_START] " << start.name << " forward=" << start.value
				          << " mid=" << target->getVal() << std::endl;
				if (auto *saved = dynamic_cast<const RooRealVar *>(npTimeResult->floatParsFinal().find(start.name)))
					target->setError(saved->getError());
				continue;
			}
			if (start.value < target->getMin() || start.value > target->getMax())
			{
				// KST - Preserve the forward guard; mid may retain a valid native seed without changing its bounds.
				if (!isMid)
					throw std::runtime_error(std::string("Native sigma start out of bounds: ") + start.name);
				// A forward numerical seed may be outside the mid-prefit-derived domain.
				// Retain the valid native mid start; keep bounds and priors unchanged.
				std::cout << "[MID_START_FALLBACK] " << start.name << " forward=" << start.value
				          << " mid=" << target->getVal() << " range=[" << target->getMin()
				          << "," << target->getMax() << "]" << std::endl;
				continue;
			}
			target->setVal(start.value);
			if (start.error > 0)
				target->setError(start.error);
		}
		if (settingPtLow == 2 && !zoomOnly)
			write_objective_state(diagnosticPath("NATIVE_MAIN_INITIAL_PARAMETERS.csv"), *parameters);
	}
	const double npSecondLifetimeInitial = signal2_lifetime.getVal();
	const double npFractionInitial = npCompositionFrac1.getVal();
	if (!zoomOnly)
	{
		std::ofstream audit(diagnosticPath("EXTERNAL_CONSTRAINTS.txt"));
		for (auto *arg : massConstraints)
			audit << arg->GetName() << '\n';
	}
	const double npFirstLifetimeInitial = signal_lifetime.getVal();
	// Reproduce the selected starts while keeping B free in the likelihood.
	// Native-sigma bins already received their complete starting coordinates above.
	double bInitial = unset;
	if (!isMid && std::isnan(bInitial))
	{
		if (settingPtLow == 6)
			bInitial = 0.1;
		if (settingPtLow == 8 || settingPtLow == 9)
			bInitial = independentNativeSigmas ? bFraction.getVal() : 0.5;
		if (settingPtLow == 10)
			bInitial = 0.3389741915833851;
		if (settingPtLow == 12)
			bInitial = 0.406727298050799;
	}
	if (!std::isnan(bInitial))
	{
		if (bInitial <= bFraction.getMin() || bInitial >= bFraction.getMax())
			throw std::runtime_error("B starting value must lie strictly inside its fit domain");
		bFraction.setVal(bInitial);
	}

    if (policy.calibrationSteps && !drawFromSavedFit) {
        auto parameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
        for (auto *arg : *parameters) {
            auto *v = dynamic_cast<RooRealVar *>(arg);
            if (!v || v->isConstant()) continue;
            for (auto *file : {prFile.get(), npFile.get(), bkgFile.get()}) {
                auto *r = file->Get<RooFitResult>("timeResult");
                auto *saved = r ? dynamic_cast<const RooRealVar *>(r->floatParsFinal().find(v->GetName())) : nullptr;
                if (saved && saved->getError()>0) {v->setError(saved->getError());break;}
            }
        }
    }
    if (!isMid && policy.calibrationSteps && !drawFromSavedFit) {
        auto parameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
        unsigned count=0;const auto *starts=oo_f1_initial_values(count);
        for (unsigned i=0;i<count;++i) {
            auto *v=dynamic_cast<RooRealVar *>(parameters->find(starts[i].name));
            if (!v || v->isConstant() || starts[i].value<v->getMin() || starts[i].value>v->getMax())
                throw std::runtime_error(std::string("Invalid F1 numerical start: ")+starts[i].name);
            v->setVal(starts[i].value);v->setError(starts[i].error);
        }
    }
	double configuredBInitial = bFraction.getVal();
	int fitMinimizeStatus = -1;
	int fitHesseStatus = -1;
	std::unique_ptr<TFile> saved2dFitFile;

	// Fit-result ownership: modelResultHolder owns newly created or restored
	// results, while model_result is the non-owning pointer used below.
	std::unique_ptr<RooFitResult> modelResultHolder;
	RooFitResult *model_result = nullptr;

	// publish mode: not fit, just re-draw plots
	if (drawFromSavedFit)
	{
		if (!load_saved_fit_file(saved2dFitFile, fitRootName, "2D fit"))
			return;
		modelResultHolder = clone_saved_fit_result(saved2dFitFile.get(), "modelResult");
		if (!modelResultHolder)
		{
			std::cerr << "ERROR: modelResult not found in saved 2D fit file: " << fitRootName << std::endl;
			return;
		}
		auto parameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
		restore_objective_result(*modelResultHolder, *parameters);
		model_result = modelResultHolder.get();
		std::cout << "[PlotOnly] Loaded saved 2D fit and skipped fit: " << fitRootName << std::endl;
	}
	else if (isForward && std::isfinite(fitTolerance))
    {
        auto parameters=objective_parameters(model,massConstraints,RooArgSet(obs_mass,obs_time));
        struct Coordinate {RooRealVar *var;double value,error;bool constant;};
        std::vector<Coordinate> initial;
        for(auto *arg:*parameters)if(auto *v=dynamic_cast<RooRealVar *>(arg))
            initial.push_back({v,v->getVal(),v->getError(),v->isConstant()});
        std::unique_ptr<RooAbsReal> rawObjective(model.createNLL(*data,EvalBackend("cpu"),Extended(true),Offset(false),ExternalConstraints(massConstraints)));
        std::ofstream attempts(diagnosticPath("TIGHT_MINIMIZER_ATTEMPTS.csv"));
        attempts << std::setprecision(17) << "start_id,initial_B,nuisance_shift,attempt,tolerance,minimize,hesse,status,covQual,edm,raw_nll,B\n";
        std::unique_ptr<RooFitResult> best;
        double bestNll=std::numeric_limits<double>::infinity();
        int bestMinimize=-1,bestHesse=-1;
        const double bStarts[]={.15,.2,.65,.2,.2};
        const double shifts[]={0.,0.,0.,-.5,.5};
        for(int start=0;start<5;++start) {
            for(auto &x:initial) {
                x.var->setConstant(x.constant);
                x.var->setError(x.error);
                double value=x.value;
                const std::string name=x.var->GetName();
                if(!x.constant && name!="bFraction" && name!="Nsig" && name!="Nbkg"
                   && std::isfinite(x.error) && x.error>0) {
                    value+=shifts[start]*x.error;
                    const double margin=1e-8*std::max(1.,std::abs(x.value));
                    value=std::clamp(value,x.var->getMin()+margin,x.var->getMax()-margin);
                }
                x.var->setVal(value);
            }
            bFraction.setVal(bStarts[start]);
            write_objective_state(diagnosticPath(Form("INITIAL_START_%d.csv",start)),*parameters);
            std::unique_ptr<RooAbsReal> nll(model.createNLL(*data,EvalBackend("cpu"),Extended(true),Offset(true),ExternalConstraints(massConstraints)));
            RooMinimizer minimizer(*nll);
            minimizer.setMinimizerType(policy.minimizer);
            minimizer.setStrategy(2);
            minimizer.fitter()->Config().MinimizerOptions().SetPrecision(policy.functionPrecision);
            minimizer.setEps(fitTolerance);minimizer.setPrintLevel(1);
            minimizer.optimizeConst(0);minimizer.setRecoverFromNaNStrength(1.0);
            for(int attempt=0;attempt<2;++attempt) {
                const int ms=minimizer.minimize(policy.minimizer,"migrad");
                const int hs=minimizer.hesse();
                std::unique_ptr<RooFitResult> result(minimizer.save());
                const double raw=rawObjective->getVal();
                const bool valid=ms==0&&hs==0&&result->status()==0&&result->covQual()==3
                    &&std::isfinite(result->edm())&&result->edm()>=0&&result->edm()<1e-3&&std::isfinite(raw);
                attempts << start << ',' << bStarts[start] << ',' << shifts[start] << ',' << attempt << ','
                         << fitTolerance << ',' << ms << ',' << hs << ',' << result->status() << ','
                         << result->covQual() << ',' << result->edm() << ',' << raw << ',' << bFraction.getVal() << '\n';
                attempts.flush();
                if(valid && raw<bestNll) {
                    best.reset(static_cast<RooFitResult *>(result->Clone()));bestNll=raw;
                    bestMinimize=ms;bestHesse=hs;configuredBInitial=bStarts[start];
                }
                if(valid&&result->edm()<1e-6)break;
            }
        }
        if(!best)throw std::runtime_error("No valid minimum among five configured start slots");
        restore_objective_result(*best,*parameters);
        modelResultHolder=std::move(best);model_result=modelResultHolder.get();
        fitMinimizeStatus=bestMinimize;fitHesseStatus=bestHesse;
        std::cout << "[COMMON_MULTISTART] chosen B seed=" << configuredBInitial
                  << " raw NLL=" << std::setprecision(17) << bestNll << std::endl;
    }
	else if (isMid && std::isfinite(fitTolerance))
    {
        // Identical parameter treatment in every start. Nominal B is never an input.
        auto parameters=objective_parameters(model,massConstraints,RooArgSet(obs_mass,obs_time));
        struct InitialCoordinate {RooRealVar *var;double value,error;bool constant;};
        std::vector<InitialCoordinate> initial;
        for(auto *arg:*parameters) if(auto *v=dynamic_cast<RooRealVar *>(arg))
            initial.push_back({v,v->getVal(),v->getError(),v->isConstant()});
        std::unique_ptr<RooAbsReal> rawObjective(model.createNLL(*data,EvalBackend("cpu"),Extended(true),Offset(false),ExternalConstraints(massConstraints)));
        std::ofstream attempts(diagnosticPath("TIGHT_MINIMIZER_ATTEMPTS.csv"));
        attempts << std::setprecision(17) << "initial_B,attempt,tolerance,minimize,hesse,status,covQual,edm,raw_nll,B\n";
        std::unique_ptr<RooFitResult> best;
        double bestNll=std::numeric_limits<double>::infinity();
        bool bestValid=false;
        int bestMinimize=-1,bestHesse=-1;
        for(double initialB:{0.15,0.2,0.65}) {
            for(auto &x:initial) {x.var->setVal(x.value);x.var->setError(x.error);x.var->setConstant(x.constant);}
            bFraction.setVal(initialB);
            std::unique_ptr<RooAbsReal> nll(model.createNLL(*data,EvalBackend("cpu"),Extended(true),Offset(true),ExternalConstraints(massConstraints)));
            RooMinimizer minimizer(*nll);
            minimizer.setMinimizerType(policy.minimizer);
            minimizer.setStrategy(2);
            if(policy.functionPrecision>0)minimizer.fitter()->Config().MinimizerOptions().SetPrecision(policy.functionPrecision);
            minimizer.setEps(fitTolerance);minimizer.setPrintLevel(1);
            minimizer.optimizeConst(0);minimizer.setRecoverFromNaNStrength(1.0);
            for(int attempt=0;attempt<2;++attempt) {
                const int ms=minimizer.minimize(policy.minimizer,"migrad");
                const int hs=minimizer.hesse();
                std::unique_ptr<RooFitResult> result(minimizer.save());
                const double raw=rawObjective->getVal();
                const bool valid=ms==0&&hs==0&&result->status()==0&&result->covQual()==3&&std::isfinite(result->edm())&&result->edm()>=0&&result->edm()<1e-3&&std::isfinite(raw);
                attempts << initialB << ',' << attempt << ',' << fitTolerance << ',' << ms << ',' << hs << ',' << result->status() << ',' << result->covQual() << ',' << result->edm() << ',' << raw << ',' << bFraction.getVal() << '\n';
                attempts.flush();
                if(std::isfinite(raw)&&(!best||(valid&&!bestValid)||(valid==bestValid&&raw<bestNll))) {
                    best.reset(static_cast<RooFitResult *>(result->Clone()));bestNll=raw;bestValid=valid;
                    bestMinimize=ms;bestHesse=hs;configuredBInitial=initialB;
                }
                if(valid&&result->edm()<std::max(1e-6,1e-3*fitTolerance))break;
            }
        }
        if(!best)throw std::runtime_error("No finite result from the three initial B values");
        restore_objective_result(*best,*parameters);
        modelResultHolder=std::move(best);model_result=modelResultHolder.get();
        fitMinimizeStatus=bestMinimize;fitHesseStatus=bestHesse;
        std::cout << "[MULTISTART] chosen B seed=" << configuredBInitial << " raw NLL=" << std::setprecision(17) << bestNll << " valid=" << bestValid << std::endl;
    }
	else if (std::isfinite(fitTolerance))
	{
		// Explicit RooMinimizer control: RooFit::fitTo resets the global tolerance.
		std::unique_ptr<RooAbsReal> nll(model.createNLL(*data, EvalBackend("cpu"), Extended(true), Offset(true), ExternalConstraints(massConstraints)));
		RooMinimizer minimizer(*nll);
		minimizer.setMinimizerType(policy.minimizer);
		minimizer.setStrategy(2);
        if (policy.functionPrecision > 0) minimizer.fitter()->Config().MinimizerOptions().SetPrecision(policy.functionPrecision);
		minimizer.setEps(fitTolerance);
		minimizer.setPrintLevel(1);
		minimizer.optimizeConst(0);
		minimizer.setRecoverFromNaNStrength(1.0);
		std::ofstream attempts(diagnosticPath("TIGHT_MINIMIZER_ATTEMPTS.csv"));
		attempts << std::setprecision(17) << "attempt,tolerance,minimize,hesse,status,covQual,edm,nll\n";
		std::unique_ptr<RooFitResult> bestValidAttempt;
		for (int attempt = 0; attempt < 2; ++attempt)
		{
			fitMinimizeStatus = minimizer.minimize(policy.minimizer, "migrad");
			fitHesseStatus = minimizer.hesse();
			modelResultHolder.reset(minimizer.save());
			model_result = modelResultHolder.get();
			attempts << attempt << ',' << fitTolerance << ',' << fitMinimizeStatus << ',' << fitHesseStatus << ',' << model_result->status() << ',' << model_result->covQual() << ',' << model_result->edm() << ',' << model_result->minNll() << '\n';
            attempts.flush();
            if (fitMinimizeStatus == 0 && fitHesseStatus == 0 && model_result->status() == 0
                && model_result->covQual() == 3 && std::isfinite(model_result->edm())
                && model_result->edm() >= 0 && std::isfinite(model_result->minNll())
                && (!bestValidAttempt || model_result->minNll() < bestValidAttempt->minNll()))
                bestValidAttempt.reset(static_cast<RooFitResult *>(model_result->Clone()));
			if (fitMinimizeStatus == 0 && fitHesseStatus == 0 && model_result->status() == 0 && model_result->covQual() == 3 && model_result->edm() < std::max(1e-6, 1e-3 * fitTolerance))
				break;
		}
        if (bestValidAttempt) {
            auto parameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
            restore_objective_result(*bestValidAttempt, *parameters);
            modelResultHolder = std::move(bestValidAttempt);
            model_result = modelResultHolder.get();
            fitMinimizeStatus = 0;
            fitHesseStatus = 0;
        }

	}
	else if (recoveryFit)
	{
		// 2026-09-10 03:22 KST - Regenerate recovery initialization from fresh prefits, never a saved fit input.
		std::unique_ptr<RooArgSet> coldObjective(model.getParameters(RooArgSet(obs_mass, obs_time)));
		for (auto *arg : massConstraints)
		{
			auto *constraint = dynamic_cast<RooAbsPdf *>(arg);
			if (!constraint)
				throw std::runtime_error("Recovery constraint is not a PDF");
			std::unique_ptr<RooArgSet> auxiliary(constraint->getParameters(RooArgSet(obs_mass, obs_time)));
			coldObjective->add(*auxiliary, true);
		}
		std::unique_ptr<RooArgSet> coldSnapshot(static_cast<RooArgSet *>(coldObjective->snapshot()));
		if (!coldSnapshot)
			throw std::runtime_error("Cannot snapshot recovery initial objective");
		auto fitPrecursor = [&]()
		{
			return std::unique_ptr<RooFitResult>(model.fitTo(*data, EvalBackend("cpu"), Save(true), Extended(true), Offset(true),
																											 RecoverFromUndefinedRegions(1.0), Strategy(2), Minimizer("Minuit2"), Optimize(false), PrintLevel(1), ExternalConstraints(massConstraints)));
		};
		std::ofstream precursorAudit(diagnosticPath("RECOVERY_PRECURSOR_ATTEMPTS.csv"));
		precursorAudit << std::setprecision(17) << "attempt,status,covQual,edm,raw_nll\n";
		std::unique_ptr<RooFitResult> precursor;
		for (int attempt = 0; attempt < 2; ++attempt)
		{
			precursor = fitPrecursor();
			if (!precursor)
				throw std::runtime_error("Fresh recovery precursor returned no result");
			std::unique_ptr<RooAbsReal> precursorRaw(model.createNLL(*data, EvalBackend("cpu"), Extended(true), Offset(false), ExternalConstraints(massConstraints)));
			double value = precursorRaw->getVal();
			precursorAudit << attempt << ',' << precursor->status() << ',' << precursor->covQual() << ',' << precursor->edm() << ',' << value << '\n';
			precursorAudit.flush();
			{
				TFile output(diagnosticPath(Form("recovery_precursor_attempt%d.root", attempt)), "RECREATE");
				precursor->Write("modelResult");
				TParameter<double>("rawFinalNll", value).Write();
			}
			if (attempt == 1 || !(precursor->status() != 0 || precursor->covQual() < 3 || precursor->edm() > 1e-3))
				break;
		}
		// 2026-09-10 03:22 KST - Restore cold auxiliary starts before the original model-only warm transfer.
		for (auto *arg : *coldSnapshot)
			if (auto *saved = dynamic_cast<RooRealVar *>(arg))
			{
				auto *target = dynamic_cast<RooRealVar *>(coldObjective->find(saved->GetName()));
				if (!target)
					throw std::runtime_error("Missing cold recovery objective coordinate");
				target->setVal(saved->getVal());
				target->setError(saved->getError());
				target->setConstant(saved->isConstant());
			}
		auto *warm = precursor.get();
		std::unique_ptr<RooArgSet> pars(model.getParameters(RooArgSet(obs_mass, obs_time)));
		for (auto *arg : *pars)
			if (auto *v = dynamic_cast<RooRealVar *>(arg))
			{
				if (v->isConstant())
					continue;
				auto *w = dynamic_cast<const RooRealVar *>(warm->floatParsFinal().find(v->GetName()));
				if (!w)
					throw std::runtime_error("Missing floating warm coordinate");
				v->setVal(w->getVal());
				v->setError(w->getError());
			}
		std::unique_ptr<RooAbsReal> objective(model.createNLL(*data, EvalBackend("cpu"), Extended(true), Offset(true), ExternalConstraints(massConstraints)));
		std::unique_ptr<RooAbsReal> raw(model.createNLL(*data, EvalBackend("cpu"), Extended(true), Offset(false), ExternalConstraints(massConstraints)));
		std::cout << "[RECOVERY_WARM_NLL] " << std::setprecision(17) << raw->getVal() << std::endl;
		if (!std::isfinite(raw->getVal()))
			throw std::runtime_error("Nonfinite warm NLL");
		RooMinimizer minimizer(*objective);
		minimizer.setMinimizerType("Minuit2");
		minimizer.setStrategy(2);
		minimizer.setPrintLevel(1);
		minimizer.setEps(.001);
		minimizer.setMaxFunctionCalls(200000);
		minimizer.setMaxIterations(10000);
		minimizer.setOffsetting(true);
		int mig = minimizer.migrad();
		int hes = minimizer.hesse();
		modelResultHolder.reset(minimizer.save());
		model_result = modelResultHolder.get();
		std::cout << "[RECOVERY_CODES] migrad=" << mig << " hesse=" << hes << " rawFinal=" << raw->getVal() << std::endl;
	}

	else if (profileNPShape)
	{
		std::unique_ptr<RooAbsReal> nll(model.createNLL(*data, EvalBackend("legacy"), Extended(true), Offset(false), ExternalConstraints(massConstraints)));
		RooMinimizer minimizer(*nll);
		minimizer.setMinimizerType("Minuit2");
		minimizer.setStrategy(2);
		minimizer.setEps(.001);
		minimizer.setPrintLevel(1);
		minimizer.optimizeConst(0);
		std::ofstream attempts(diagnosticPath("FIT_ATTEMPTS.csv"));
		attempts << std::setprecision(17) << "attempt,status,covQual,edm,minimize,hesse,raw_nll,valid\n";
		double bestNll = std::numeric_limits<double>::infinity();
		bool bestValid = false;
		for (int attempt = 0; attempt < 2; ++attempt)
		{
			int ms = minimizer.minimize("Minuit2", "migrad"), hs = minimizer.hesse();
			auto result = std::unique_ptr<RooFitResult>(minimizer.save());
			double raw;
			{
				std::unique_ptr<RooAbsReal> r(model.createNLL(*data, EvalBackend("legacy"), Extended(true), Offset(false), ExternalConstraints(massConstraints)));
				raw = r->getVal();
			}
			bool valid = result->status() == 0 && result->covQual() == 3 && ms == 0 && hs == 0 && std::isfinite(raw);
			attempts << attempt << ',' << result->status() << ',' << result->covQual() << ',' << result->edm() << ',' << ms << ',' << hs << ',' << raw << ',' << valid << '\n';
			attempts.flush();
			{
				TFile archive(diagnosticPath(Form("fit_attempt%d.root", attempt)), "RECREATE");
				result->Write("modelResult");
				TParameter<double>("rawFinalNll", raw).Write();
				TParameter<int>("minimizeStatus", ms).Write();
				TParameter<int>("hesseStatus", hs).Write();
			}
			if (std::isfinite(raw) && ((valid && !bestValid) || (valid == bestValid && raw < bestNll)))
			{
				bestNll = raw;
				bestValid = valid;
				modelResultHolder.reset(static_cast<RooFitResult *>(result->Clone()));
				fitMinimizeStatus = ms;
				fitHesseStatus = hs;
			}
			if (valid && result->edm() < 1e-3)
				break;
		}
		if (!modelResultHolder)
			throw std::runtime_error("No finite fit attempt");
		{
			auto objective = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
			restore_objective_result(*modelResultHolder, *objective);
			write_objective_state(diagnosticPath("FINAL_OBJECTIVE_PARAMETERS.csv"), *objective);
		}
		model_result = modelResultHolder.get();
	}
	else
	{
		// Fit mode: first perform the nominal likelihood minimization, then keep
		// its RooFitResult alive in modelResultHolder for diagnostics and retries.

		// Nominal extended 2D likelihood fit
		if (massConstraints.getSize() > 0)
		{
			model_result = model.fitTo(*data,
																 Save(true),
																 Extended(true),
																 Offset(true),
																 RecoverFromUndefinedRegions(1.0),
																 Strategy(2),
																 Minimizer("Minuit2"),
																 Optimize(false),
																 PrintLevel(1),
																 ExternalConstraints(massConstraints));
		}
		else
		{
			model_result = model.fitTo(*data,
																 Save(true),
																 Extended(true),
																 Offset(true),
																 RecoverFromUndefinedRegions(1.0),
																 Strategy(2),
																 Minimizer("Minuit2"),
																 Optimize(false),
																 PrintLevel(1));
		}
		modelResultHolder.reset(model_result);
		model_result = modelResultHolder.get();

		// Retry once when the first pass did not converge cleanly.  The existing
		// forward pT 1-2 exception controls the retry offset convention.
		if (model_result && (model_result->status() != 0 ||
												 model_result->covQual() < 3 || model_result->edm() > 1e-3))
		{
			if (massConstraints.getSize() > 0)
			{
				modelResultHolder.reset(model.fitTo(*data,
																						Save(true),
																						Extended(true),
																						Offset(true),
																						RecoverFromUndefinedRegions(1.0),
																						Strategy(2),
																						Minimizer("Minuit2"),
																						Optimize(false),
																						PrintLevel(1),
																						ExternalConstraints(massConstraints)));
			}
			else
			{
				modelResultHolder.reset(model.fitTo(*data,
																						Save(true),
																						Extended(true),
																						Offset(true),
																						RecoverFromUndefinedRegions(1.0),
																						Strategy(2),
																						Minimizer("Minuit2"),
																						Optimize(false),
																						PrintLevel(1)));
			}
			model_result = modelResultHolder.get();
		}
	}

	// Extract the minimizer and HESSE status from RooFit's history
	if (model_result && (!profileNPShape || drawFromSavedFit) && (!std::isfinite(fitTolerance) || drawFromSavedFit))
	{
		for (UInt_t i = 0, n = model_result->numStatusHistory(); i < n; ++i)
		{
			const char *lab = model_result->statusLabelHistory(i);
			if (!lab)
				continue;
			if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
				fitMinimizeStatus = model_result->statusCodeHistory(i);
			else if (!strcmp(lab, "HESSE"))
				fitHesseStatus = model_result->statusCodeHistory(i);
		}
	}

	if (drawFromSavedFit && saved2dFitFile)
	{
		fitMinimizeStatus = readIntParam(*saved2dFitFile, "minimizeStatus", fitMinimizeStatus);
		fitHesseStatus = readIntParam(*saved2dFitFile, "hesseStatus", fitHesseStatus);
	}

	if (!zoomOnly)
	{
		// KST - Save model-only and full-objective coordinates to distinguish physical from auxiliary bounds.
		auto modelParameters = std::unique_ptr<RooArgSet>(model.getParameters(RooArgSet(obs_mass, obs_time)));
		write_objective_state(diagnosticPath("FINAL_MODEL_PARAMETERS.csv"), *modelParameters);
		auto objectiveParameters = objective_parameters(model, massConstraints, RooArgSet(obs_mass, obs_time));
		write_objective_state(diagnosticPath("FINAL_OBJECTIVE_PARAMETERS.csv"), *objectiveParameters);
	}

	// Propagate the final 2D covariance matrix to the reported nonprompt
	double bFractionErr = -1.0;
	if (model_result)
		bFractionErr = bFraction.getPropagatedError(*model_result);

	// Report B and its covariance error in the actual fitted ctau window only.
	// Range changes require sufficient signal retention; never extrapolate B.
	if (outputDir[0])
	{
		std::ofstream audit(diagnosticPath("FRACTION_WINDOW_AUDIT.csv"));
		audit << std::setprecision(17) << "fit_low,fit_high,B_in_fit_window,B_fit_error\n"
					<< ctRange.first << ',' << ctRange.second << ','
					<< bFraction.getVal() << ',' << bFractionErr << '\n';
	}

	// ------------------------------------------------------------------
	// draw - 1D mass
	// ------------------------------------------------------------------
	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		if (fr)
			return fr->findObject(n);
		return nullptr;
	};

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
	model.plotOn(massfitplot, Precision(1e-6), LineColor(kBlack), LineWidth(2), Name("model"));
	model.plotOn(massfitplot, Precision(1e-6), Components(signal_core), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Name("signal_component"));
	model.plotOn(massfitplot, Precision(1e-6), Components(bkg_core), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("bkg_component"));
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
	if (!publish)
	{
		TLatex tx;
		tx.SetNDC();
		tx.SetTextSize(0.03);
		tx.SetTextFont(42);
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", fitMinimizeStatus, fitHesseStatus));
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
		if (publish)
		{
			printVar("N_{sig}", Nsig);
			printVar("N_{bkg}", Nbkg);
			printVar("f_{B}", bFraction);
		}
		else
		{
			printVar("N_{sig}", Nsig);
			printVar("N_{bkg}", Nbkg);
			printVar("f_{B}", bFraction);
			printVar("s_{1}", ctauTime1Scale);
			if (nResolutionComponents >= 2)
			{
				printVar("s_{2}", ctauTime2Scale);
				if (!independentNativeSigmas)
					printVar("#Delta s_{21}", ctauTime2Delta);
			}
			if (nResolutionComponents >= 3)
			{
				printVar("s_{3}", ctauTime3Scale);
				if (!independentNativeSigmas)
					printVar("#Delta s_{32}", ctauTime3Delta);
			}
			if (nResolutionComponents >= 4)
			{
				printVar("s_{4}", ctauTime4Scale);
				if (!independentNativeSigmas)
					printVar("#Delta s_{43}", ctauTime4Delta);
			}
			printVar("#tau_{sig1}", signal_lifetime);
			if (useSignalSS2)
				printVar("#tau_{sig2}", signal2_lifetime);
			printVar("#tau_{bkg1}", bkg_lifetime);
			if (nBkgSSComponents >= 2)
				printVar("#tau_{bkg2}", bkg_lifetime2);
			if (nBkgSSComponents >= 3)
				printVar("#tau_{bkg3}", bkg_lifetime3);
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

	// mass uses mass-only plus common parameters; lifetime-only parameters are not subtracted.
	auditConditionalBackground("final");
	const ProjectionParameterSummary projectionParameters = classify_projection_parameters(model_result);
	if (!projectionParameters.unknown.empty())
	{
		std::cerr << "ERROR: unclassified floating fit parameters: "
							<< join_parameter_names(projectionParameters.unknown) << std::endl;
		delete cMass;
		return;
	}
	int totalFloatNpar = 0;
	if (model_result)
		totalFloatNpar = model_result->floatParsFinal().getSize();
	const int massNpar = projectionParameters.massOnly.size() + projectionParameters.common.size();
	const int timeNpar = projectionParameters.lifetimeOnly.size() + projectionParameters.common.size();
	int massDataBins = prefit_chi2_bins(*massfitplot);
	const int massNdf = massDataBins - massNpar;
	if (massNdf <= 0)
	{
		std::cerr << "ERROR: non-positive mass NDF: bins=" << massDataBins
							<< " npar=" << massNpar << std::endl;
		delete cMass;
		return;
	}
	const double massChi2Ndf = massfitplot->chiSquare("model", "data", massNpar);
	const double massChi2 = massChi2Ndf * massNdf;
	const double massPvalue = TMath::Prob(massChi2, massNdf);
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", massChi2, massNdf));
	}
	TLine massLine(obs_mass.getMin(), 0.0, obs_mass.getMax(), 0.0);
	massLine.SetLineStyle(2);
	massLine.Draw("same");
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};

	if (!zoomOnly) cMass->Print(figName("mass_fit"));
	delete cMass;

	// ------------------------------------------------------------------
	// draw - 1D ctau
	// ------------------------------------------------------------------
	// note - Display windows only: retain the full fit selection, PDF normalization,
	// data binning and projection chi2. Zoom the axes after constructing the curves.
	std::pair<double, double> ctPlotRange = ctRange;
	if (settingPtLow == 2 && settingPtHigh == 3)
		ctPlotRange = {-5.0, 5.0};
	if (settingPtLow == 3 && settingPtHigh == 4)
		ctPlotRange = {-4.0, 4.0};
	if (settingPtLow == 4 && settingPtHigh == 5)
		ctPlotRange = {-2.1, 3.5};
	if (settingPtLow == 10 && settingPtHigh == 12)
		ctPlotRange = {-1.0, 2.0};
	if (settingPtLow == 12 && settingPtHigh == 14)
		ctPlotRange = {-0.5, 2.0};
	if (zoomOnly) ctPlotRange = {zoomLow, zoomHigh};
	ctPlotRange.first = std::max(ctPlotRange.first, ctRange.first);
	ctPlotRange.second = std::min(ctPlotRange.second, ctRange.second);
	if (ctPlotRange.first >= ctPlotRange.second)
		ctPlotRange = ctRange;
    if (zoomOnly) {
        obs_time.setRange("signalZoom", ctPlotRange.first, ctPlotRange.second);
        const RooArgSet timeNorm(obs_time);
        for (RooAbsPdf *pdf : {promptTimePdf, signal_np_time}) {
            std::unique_ptr<RooAbsReal> fraction(pdf->createIntegral(timeNorm, NormSet(timeNorm), Range("signalZoom")));
            const double retained = fraction->getVal();
            std::cout << "[ZOOM] " << pdf->GetName() << " retained=" << 100*retained << "%" << std::endl;
            if (!std::isfinite(retained) || retained < .99)
                throw std::runtime_error("Widen the display window to retain at least 99% of each PR/NP signal");
        }
    }
	const bool zoomCtauPlot = ctPlotRange != ctRange;
    if (zoomOnly) {
        const int visibleBins=oo_diagnostics::display_bins(*data,obs_time,ctPlotRange.first,ctPlotRange.second);
        timePlotBins=std::max(timePlotBins, int(std::ceil(visibleBins*(ctRange.second-ctRange.first)/(ctPlotRange.second-ctPlotRange.first))));
        std::cout << "[ZOOM_BINNING] visible bins >= " << visibleBins << ", full-range bins=" << timePlotBins << std::endl;
    }

	std::cout << "[PLOT] ctau display=[" << ctPlotRange.first << "," << ctPlotRange.second
						<< "], fit=[" << ctRange.first << "," << ctRange.second
						<< "]; chi2 uses full fit range" << std::endl;
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
	// One interpolation precision and full fitted-yield normalization for all bins.
    const double ctauCurvePrecision = 1e-6;
    model.plotOn(timefitplot, Precision(ctauCurvePrecision),
                 Normalization(Nsig.getVal() + Nbkg.getVal(), RooAbsReal::NumEvent),
                 LineColor(kBlack), LineWidth(2), Name("model"));
	const double npYield = Nsig.getVal() * bFraction.getVal();
	const double prYield = std::max(0.0, Nsig.getVal() - npYield);
	prompt_core.plotOn(timefitplot, Precision(ctauCurvePrecision), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2),
										 Normalization(prYield, RooAbsReal::NumEvent), Name("pr_component"));
	np_core.plotOn(timefitplot, Precision(ctauCurvePrecision), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2),
								 Normalization(npYield, RooAbsReal::NumEvent), Name("np_component"));
	bkg_core.plotOn(timefitplot, Precision(ctauCurvePrecision), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2),
									Normalization(Nbkg.getVal(), RooAbsReal::NumEvent), Name("bkg_component"));
	apply_logy_auto_range(timefitplot, "data");
	timefitplot->GetYaxis()->SetTitle("Events");
	timefitplot->GetYaxis()->SetTitleOffset(1.6);
	timefitplot->GetXaxis()->SetTitle("");
	if (zoomCtauPlot)
		timefitplot->GetXaxis()->SetLimits(ctPlotRange.first, ctPlotRange.second);
	timefitplot->Draw("e");

	TLegend timeLeg(0.19, 0.52, 0.43, 0.725);
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
	if (!publish)
	{
		TLatex tx;
		tx.SetNDC();
		tx.SetTextSize(0.03);
		tx.SetTextFont(42);
		const int minimize = fitMinimizeStatus;
		const int hesse = fitHesseStatus;
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
		if (publish)
		{
			printVal("N_{sig}", Nsig.getVal());
			printVal("N_{bkg}", Nbkg.getVal());
			printVal("f_{B}", bFraction.getVal(), bFractionErr);
		}
		else
		{
			printVal("f_{B}", bFraction.getVal(), bFractionErr);
			printVal("N_{sig}", Nsig.getVal());
			printVal("N_{bkg}", Nbkg.getVal());
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
	if (zoomCtauPlot)
		timePullPlot->GetXaxis()->SetLimits(ctPlotRange.first, ctPlotRange.second);
	timePullPlot->Draw();

	// lifetime uses lifetime-only plus common parameters; mass-only parameters are not subtracted.
	int timeDataBins = prefit_chi2_bins(*timefitplot);
	const int timeNdf = timeDataBins - timeNpar;
	if (timeNdf <= 0)
	{
		std::cerr << "ERROR: non-positive lifetime NDF: bins=" << timeDataBins
							<< " npar=" << timeNpar << std::endl;
		delete cLifetime;
		return;
	}
	const double timeChi2Ndf = timefitplot->chiSquare("model", "data", timeNpar);
	const double timeChi2 = timeChi2Ndf * timeNdf;
	const double timePvalue = TMath::Prob(timeChi2, timeNdf);
	std::cout << "[NDF] total floating parameters=" << totalFloatNpar << std::endl;
	std::cout << "[NDF] mass-only=" << projectionParameters.massOnly.size()
						<< " [" << join_parameter_names(projectionParameters.massOnly) << "]" << std::endl;
	std::cout << "[NDF] lifetime-only=" << projectionParameters.lifetimeOnly.size()
						<< " [" << join_parameter_names(projectionParameters.lifetimeOnly) << "]" << std::endl;
	std::cout << "[NDF] common=" << projectionParameters.common.size()
						<< " [" << join_parameter_names(projectionParameters.common) << "]" << std::endl;
	std::cout << "[NDF] mass bins=" << massDataBins << " npar=" << massNpar << " ndf=" << massNdf << std::endl;
	std::cout << "[NDF] lifetime bins=" << timeDataBins << " npar=" << timeNpar << " ndf=" << timeNdf << std::endl;
	std::cout << Form("[FIT] mass chi2/ndf = %.1f/%d = %.3f, p = %.3g",
										massChi2, massNdf, massChi2Ndf, massPvalue)
						<< std::endl;
	std::cout << Form("[FIT] ctau chi2/ndf = %.1f/%d = %.3f, p = %.3g",
										timeChi2, timeNdf, timeChi2Ndf, timePvalue)
						<< std::endl;
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form(zoomCtauPlot ? "#chi^{2}/ndf (full range) = %.1f/%d" : "#chi^{2}/ndf = %.1f/%d", timeChi2, timeNdf));
	}
	TLine timeLine(ctPlotRange.first, 0.0, ctPlotRange.second, 0.0);
	timeLine.SetLineStyle(2);
	timeLine.Draw("same");
	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

    if (gSystem->Getenv("OO_DIAGNOSTICS_DIR")) {
        std::unique_ptr<RooAbsPdf> timeProjection(model.createProjection(RooArgSet(obs_mass)));
        std::unique_ptr<RooAbsPdf> massProjection(model.createProjection(RooArgSet(obs_time)));
        oo_diagnostics::run(*massProjection, *data, obs_mass, massNpar, obs_mass.getBins(), "fit2d_mass_"+figTag, false);
        oo_diagnostics::run(*timeProjection, *data, obs_time, timeNpar, timePlotBins, "fit2d_ctau_"+figTag, true,
            {{promptTimePdf, prYield, "PR", kRed+1}, {signal_np_time, npYield, "NP", kBlue+1}});
    }

    if (!drawFromSavedFit && (!model_result || fitMinimizeStatus != 0 || fitHesseStatus != 0
        || model_result->status() != 0 || model_result->covQual() != 3
        || !std::isfinite(model_result->edm()) || model_result->edm() < 0 || model_result->edm() >= 1e-3))
        throw std::runtime_error("September tuning final failed strict 0/3/EDM<1e-3 gate; diagnostics preserved");
	// ------------------------------------------------------------------
	// save fit result and projection diagnostics
	// ------------------------------------------------------------------
	if (!drawFromSavedFit)
	{
		std::unique_ptr<TFile> fitOut(TFile::Open(fitRootName, "RECREATE"));
		if (!fitOut || fitOut->IsZombie())
		{
			std::cerr << "ERROR: cannot create fit2d output file: " << fitRootName << std::endl;
		}
		else
		{
			if (model_result)
				model_result->Write("modelResult");
			TParameter<int>("fixedMassPolicyVersion", 1).Write();
			TParameter<int>("jointNativeNPMcConstraint", constrainNP && (isForward || !profileNPShape) ? 1 : 0).Write();
			TParameter<int>("npIdentifiableShapeMCConstraint", constrainNP && profileNPShape && !isForward ? 1 : 0).Write();
			TParameter<double>("promptFractionMCErrorScale", constrainPromptFractions ? (fullPromptCovariance ? 1.0 : (settingPtLow < 5.0f ? 2.0 : 0.0)) : 0.0).Write();
			TParameter<double>("symmetricBkgDomainPriorSigmas", 3.0).Write();
			TParameter<int>("coreFixedMixturePolicyVersion", widthPolicy == 2 ? 1 : 0).Write();
			TParameter<int>("twoRegionPolicyVersion", (isMid || fullPromptCovariance) ? 0 : 1).Write();
			TParameter<double>("appliedCtauLow", ctRange.first).Write();
			TParameter<double>("appliedCtauHigh", ctRange.second).Write();
			TParameter<int>("appliedTimePlotBins", timePlotBins).Write();
			if (std::isfinite(fitTolerance))
				TParameter<double>("explicitMinimizerTolerance", fitTolerance).Write();
			TParameter<int>("stableGaussianConvolution", 0).Write();
			TParameter<int>("bFractionFixedForProfile", bFraction.isConstant() ? 1 : 0).Write();
			TParameter<double>("bFractionInitial", configuredBInitial).Write();
            if (isMid) TParameter<int>("minimizerStartCount",3).Write();
			TParameter<double>("bFractionLower", bFraction.getMin()).Write();
			TParameter<double>("bFractionUpper", bFraction.getMax()).Write();
			if (nSignalSSComponents >= 2)
			{
				TParameter<double>("npSecondLifetimeInitial", npSecondLifetimeInitial).Write();
				TParameter<double>("npSecondLifetimeLower", signalLifetime2FloorSaved + std::exp(signal2_log_lifetime_offset.getMin())).Write();
				TParameter<double>("npSecondLifetimeUpper", signalLifetime2FloorSaved + std::exp(signal2_log_lifetime_offset.getMax())).Write();
				TParameter<double>("npSecondLifetimeFinal", signal2_lifetime.getVal()).Write();
				TParameter<double>("npFractionInitial", npFractionInitial).Write();
			}
			TParameter<double>("npFirstLifetimeInitial", npFirstLifetimeInitial).Write();
			TParameter<double>("npFirstLifetimeLower", signalLifetimeFloorSaved + std::exp(signal_log_lifetime_offset.getMin())).Write();
			TParameter<double>("npFirstLifetimeUpper", signalLifetimeFloorSaved + std::exp(signal_log_lifetime_offset.getMax())).Write();
			TParameter<double>("npFirstLifetimeFinal", signal_lifetime.getVal()).Write();
			TParameter<double>("npFirstLogLifetimeMCError", readResultError(*npTimeResult, "signal_log_lifetime_offset", 0.0)).Write();
			TParameter<double>("twoRegionPtBoundary", 5.0).Write();
			TParameter<int>("fixedPromptFractions", settingPtLow >= 5.0f ? 1 : 0).Write();
			TParameter<double>("promptBroadestAbsolutePenalty", (fullPromptCovariance || !constrainPrompt || commonPromptScale) ? 0.0 : absoluteWidthPenalty[nResolutionComponents - 1]).Write();
			TParameter<int>("fitPolicyVersion", 2).Write();
			TNamed("fitPolicyGroup", policy.group).Write();
			TParameter<double>("configuredFunctionPrecision", policy.functionPrecision).Write();
			TParameter<double>("configuredMinimizerTolerance", fitTolerance).Write();
			TParameter<int>("widthPolicy", widthPolicy).Write();
            if (isForward) {
                TParameter<double>("promptMeanMCErrorMultiplier",ptLow>=10?8.:0.).Write();
                TParameter<double>("coreWidthSigmaMultiplier",oo_integrated_core_strength(ptLow)).Write();
                TParameter<double>("promptFractionAdditionalCovarianceScale",0.).Write();
                TParameter<int>("backgroundMassPolicy",0).Write();
                TParameter<int>("fwdFourGroupPolicyVersion",1).Write();
            }
			TParameter<int>("shareSignalBackgroundResolution", 0).Write();
			TParameter<int>("nResolutionCalibrationComponents", nResolutionCalibrationComponents).Write();
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
			TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
			TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
			TParameter<double>("ptLow", ptLow).Write();
			TParameter<double>("ptHigh", ptHigh).Write();
			TParameter<double>("yLow", yLow).Write();
			TParameter<double>("yHigh", yHigh).Write();
			TParameter<int>("directFitVersion", 1).Write();
			TParameter<int>("npDataComponentsOverride", -1).Write();
			TParameter<int>("nSignalSSCalibrationComponents", nSignalSSComponentsSaved).Write();
			TParameter<int>("promptConstraintsApplied", constrainPrompt ? 1 : 0).Write();
			TParameter<int>("promptFractionPriorApplied", constrainPromptFractions && settingPtLow < 5 && !fullPromptCovariance ? 1 : 0).Write();
			TParameter<int>("promptWidthPriorsApplied", constrainPrompt && !fullPromptCovariance && !commonPromptScale ? 1 : 0).Write();
			int promptWidthPriorCount = 0;
			for (int i = 0; i < nResolutionComponents; ++i)
				promptWidthPriorCount += constrainPrompt && !fullPromptCovariance && !promptWidths[i]->isConstant();
			TParameter<int>("promptWidthPriorCount", promptWidthPriorCount).Write();
			TParameter<int>("bkgLifetimeConstraintEnabled", constrainBackgroundLifetimes ? 1 : 0).Write();
			TParameter<int>("bkgLifetimeConstraintsApplied", (!profileBackground && constrainBackgroundLifetimes) ? 1 : 0).Write();
			TParameter<int>("npMCConstraintApplied", constrainNP ? 1 : 0).Write();
			TParameter<int>("npConstraintPolicy", constrainNP ? 1 : 0).Write();
			TParameter<int>("npLifetimeMCConstraintApplied", constrainNP ? 1 : 0).Write();
			TParameter<int>("npCompositionOnlyMCConstraintApplied", 0).Write();

			TParameter<int>("floatNativePromptMean", floatNativePromptMean ? 1 : 0).Write();
			// Legacy output names for derived diagnostics only: data mean
			// minus fixed MC mean and its fitted error (zero in fixed bins).
			// These TParameter objects are not RooFit likelihood coordinates.
			TParameter<double>("ctauDataMeanShift", floatNativePromptMean ? ctauTime1Mean.getVal() - nativePromptMean : 0.).Write();
			TParameter<double>("ctauDataMeanShiftErr", floatNativePromptMean ? dynamic_cast<RooRealVar &>(ctauTime1Mean).getError() : 0.).Write();
			TParameter<int>("floatNPLifetimeScale", 0).Write();
			TParameter<double>("npDataLifetimeScale", 1.0).Write();
			TParameter<double>("npDataLifetimeScaleErr", 0.0).Write();
			TParameter<int>("massShapeProfiled", profileMass ? 1 : 0).Write();
			TParameter<int>("massTailParametersProfiled", 0).Write();
			TParameter<double>("massTailPowerLower", 0.01).Write();
			TParameter<int>("massTailPowersFixedToMC", 0).Write();
			TParameter<double>("promptWidthGapMin", promptGapMin).Write();
			TParameter<int>("conditionalBackgroundPositiveSlope", conditionalBackground ? 1 : 0).Write();
			TParameter<int>("backgroundAndYieldPriorsReleased", profileBackground ? 1 : 0).Write();
			TParameter<double>("backgroundYieldUpperNSigma", -1).Write();
			TParameter<int>("nBkgChebyOrder", nBkgChebyOrder).Write();
			TParameter<int>("nBkgExpComponents", nBkgExpComponents).Write();
            write_oo_policy_metadata(*fitOut, policy, isMid, ptLow, nResolutionComponents, nBkgSSComponents, nBkgDSComponents, nBkgFlipComponents, nBkgResolutionComponents);
            auto actualD = [](const char *n, double v) { TParameter<double>(n,v).Write("",TObject::kOverwrite); };
            auto actualI = [](const char *n, int v) { TParameter<int>(n,v).Write("",TObject::kOverwrite); };
            const bool actualNP = massConstraints.find("npJointConstraint");
            const bool actualFraction = massConstraints.find("promptFractionConstraint");
            const bool actualFullPrompt = massConstraints.find("promptFullConstraint");
            actualI("septemberTuningVersion",2);
            TNamed("septemberTuningDefinition", "Native MC centers unchanged; scalar Gaussian sigma/native covariance multipliers only; optional declared profiling; no overall resolution scale").Write();
            actualD("septTune_width_sigma_multiplier",widthSigmaMultiplier);
            actualD("septTune_prompt_full_cov_scale",promptFullCovarianceScale);
            actualD("septTune_prompt_fraction_cov_scale",policy.fractionCovarianceScale);
            actualD("septTune_np_cov_scale",npCovarianceScale);
            actualD("septTune_background_penalty_scale",policy.backgroundPenaltyScale);
            actualI("septTune_profile_mass",profileMass);
            actualI("septTune_profile_bkg",profileBackground);
            actualI("septTune_np_longest_free",npLongestLifetimeFree);
            actualD("promptJointMCCovarianceScale",actualFullPrompt?promptFullCovarianceScale:0.);
            actualD("promptFractionMCCovarianceScale",actualFraction?policy.fractionCovarianceScale:0.);
            actualD("promptFractionMCErrorScale",actualFraction?std::sqrt(policy.fractionCovarianceScale):(actualFullPrompt?std::sqrt(promptFullCovarianceScale):0.));
            actualD("npJointMCCovarianceScale",actualNP?npCovarianceScale:0.);
            int actualWidthCount=0, actualBkgLifetimeCount=0, actualBkgYieldCount=0;
            double actualBroadestSigma=0.;
            for (int i=0;i<4;++i) {
                const bool active=i<nResolutionComponents && massConstraints.find(Form("%s_constraint",promptWidths[i]->GetName()));
                const double sigma=active?widthSigmaMultiplier*(isForward?oo_integrated_width_sigma(ptLow,i,nResolutionComponents):absoluteWidthPenalty[i]*oo_prompt_penalty_scale(isMid,ptLow,i)):0.;
                actualD(Form("promptWidthPriorSigma%d",i+1),sigma);
                actualWidthCount+=active;
                if(i==nResolutionComponents-1) actualBroadestSigma=sigma;
            }
            for (auto *arg : massConstraints) {
                const std::string name=arg->GetName();
                if (starts_with(name,"bkg_") && name.find("_constraint")!=std::string::npos) ++actualBkgLifetimeCount;
                if (name=="Nplb_constraint" || starts_with(name,"Nss") || starts_with(name,"Nds") || starts_with(name,"Nflip")) ++actualBkgYieldCount;
            }
            actualD("promptBroadestAbsolutePenalty",actualBroadestSigma);
            actualI("promptWidthPriorsApplied",actualWidthCount>0);
            actualI("promptWidthPriorCount",actualWidthCount);
            actualI("promptFractionPriorApplied",actualFraction);
            actualI("promptFractionAnyPriorApplied",actualFraction || actualFullPrompt);
            actualI("promptFullCovariancePriorApplied",actualFullPrompt);
            actualI("fixedPromptFractions",ctauRecursiveFrac1.isConstant() && ctauRecursiveFrac2.isConstant() && ctauRecursiveFrac3.isConstant());
            actualI("backgroundShapeProfiled",profileBackground);
            int resolutionFloat=0, tailFloat=0, compositionFloat=0;
            for(auto *a:model_result->floatParsFinal()) {
                const std::string n=a->GetName();
                if(starts_with(n,"bkgCtau")) ++resolutionFloat;
                if(starts_with(n,"bkg_") && n.find("lifetime")!=std::string::npos) ++tailFloat;
                if(starts_with(n,"bkgComposition") || n=="Nplb" || starts_with(n,"Nss") || starts_with(n,"Nds") || starts_with(n,"Nflip")) ++compositionFloat;
            }
            actualI("backgroundResolutionFloatCount",resolutionFloat);
            actualI("backgroundLifetimeFloatCount",tailFloat);
            actualI("backgroundCompositionFloatCount",compositionFloat);
            actualI("massShapeProfiled",profileMass);
            actualI("fixedMassPolicyVersion",profileMass?0:1);
            actualI("bkgLifetimeConstraintsApplied",actualBkgLifetimeCount>0);
            actualI("activeBackgroundLifetimeConstraintCount",actualBkgLifetimeCount);
            actualI("activeBackgroundYieldConstraintCount",actualBkgYieldCount);
            actualD("backgroundLifetimePenaltyMultiplier",actualBkgLifetimeCount?policy.backgroundPenaltyScale:0.);
            actualI("npLongestLifetimeFree",npLongestLifetimeFree);
            actualI("npConstrainedCoordinateCount",actualNP?2*nSignalSSComponents-(npLongestLifetimeFree?1:0):0);
            actualI("minimizerStartCount",isForward?5:3);
            actualD("configuredMinimizerTolerance",fitTolerance);
			TParameter<int>("fullNativePromptCovariance", fullPromptCovariance ? 1 : 0).Write();

			TParameter<int>("commonPromptResolutionScalePolicy", commonPromptScale ? 1 : 0).Write();
			TParameter<int>("independentNativeSigmaPolicy", independentNativeSigmas ? 1 : 0).Write();
			if (conditionalBackground)
			{
				TParameter<double>("bkgPositiveMassSlope", bkgPositiveMassSlope.getVal()).Write();
				TParameter<double>("bkgPositiveMassSlopeErr", bkgPositiveMassSlope.getError()).Write();
			}
			std::unique_ptr<RooAbsReal> rawNll(model.createNLL(*data, Extended(true), Offset(false), ExternalConstraints(massConstraints)));
			if (!std::isfinite(rawNll->getVal()))
				throw std::runtime_error("Nonfinite final likelihood");
			TParameter<double>("rawFinalNll", rawNll->getVal()).Write();

			TParameter<int>("backgroundNativeLifetimeDomainApplied",1).Write();
			TParameter<int>("backgroundLifetimeSeedClampCount",bkgNativeRestoreChanges).Write();
			for(int i=0;i<9;++i) if(bkgNativeActive[i]) {
				TParameter<double>(Form("native_%s_min",bkgNativeOffsets[i]->GetName()),bkgNativeOffsets[i]->getMin()).Write();
				TParameter<double>(Form("native_%s_max",bkgNativeOffsets[i]->GetName()),bkgNativeOffsets[i]->getMax()).Write();
				TParameter<double>(Form("legacy_seed_%s",bkgNativeOffsets[i]->GetName()),bkgPreRestoreOffsets[i]).Write();
			}
			TParameter<int>("bkgMeanOffsetsApplied", profileBackground ? 1 : 0).Write();
			TParameter<double>("bkgMeanBaseSaved", bkgSavedMean).Write();
			for (int i = 0; i < nBkgResolutionComponents; ++i)
				TParameter<double>(Form("bkgMeanOffset%d", i + 1), bkgSavedMeanOffsets[i]).Write();
			TParameter<int>("promptWindowTransportApplied", 1).Write();
			TParameter<double>("prCalibrationCtauLow", prTransportCalibration.low).Write();
			TParameter<double>("prCalibrationCtauHigh", prTransportCalibration.high).Write();
			TNamed("prCalibrationPrefitMD5", prTransportCalibration.md5.c_str()).Write();
			for (int i = 0; i < nResolutionComponents; ++i) {
				TParameter<double>(Form("prPhysicalFraction%d", i + 1), promptTransport.fullWeights[i]->getVal()/promptTransport.fullSum->getVal()).Write();
				TParameter<double>(Form("prWindowFraction%d", i + 1), promptTransport.windowWeights[i]->getVal()/promptTransport.windowSum->getVal()).Write();
			}
			if (transportNP)
			{
				TParameter<int>("npTransportApplied", 1).Write();
				TParameter<double>("npTransportCalibrationLow", npTransportLow).Write();
				TParameter<double>("npTransportCalibrationHigh", npTransportHigh).Write();
				double qsum = 0.;
				for (auto &q : npTransportWeights)
					qsum += q->getVal();
				for (int i = 0; i < nSignalSSComponents; ++i)
					TParameter<double>(Form("npTransportFraction%d", i + 1), npTransportWeights[i]->getVal() / qsum).Write();
			}
			windowRule.write(*fitOut);
            TNamed("windowMassPrefit", massFileName.Data()).Write();
            TNamed("windowBackgroundPrefit", bkgFileName.Data()).Write();
            TParameter<double>("bFraction", bFraction.getVal()).Write();
			TParameter<double>("bFractionErr", bFractionErr).Write();
			TParameter<double>("Nsig", Nsig.getVal()).Write();
			TParameter<double>("Nbkg", Nbkg.getVal()).Write();
			TParameter<int>("projectionNdfVersion", 3).Write();
			TParameter<int>("totalFloatNpar", totalFloatNpar).Write();
			TParameter<int>("massPullBins", massDataBins).Write();
			TParameter<int>("timePullBins", timeDataBins).Write();
			TParameter<int>("massNpar", massNpar).Write();
			TParameter<int>("timeNpar", timeNpar).Write();
			TNamed("massOnlyNparNames", join_parameter_names(projectionParameters.massOnly).c_str()).Write();
			TNamed("lifetimeOnlyNparNames", join_parameter_names(projectionParameters.lifetimeOnly).c_str()).Write();
			TNamed("commonNparNames", join_parameter_names(projectionParameters.common).c_str()).Write();
			TParameter<double>("massChi2", massChi2).Write();
			TParameter<int>("massNdf", massNdf).Write();
			TParameter<double>("massPvalue", massPvalue).Write();
			TParameter<double>("timeChi2", timeChi2).Write();
			TParameter<int>("timeNdf", timeNdf).Write();
			TParameter<double>("timePvalue", timePvalue).Write();
			int fitStatus = -1;
			if (model_result)
				fitStatus = model_result->status();
			TParameter<int>("fitStatus", fitStatus).Write();
			TParameter<int>("hesseStatus", fitHesseStatus).Write();
			TParameter<int>("minimizeStatus", fitMinimizeStatus).Write();
			fitOut->Write();
			// Replace the old dispatcher's sidecar when this macro produces a new result.
			std::ofstream provenance(fitRootName + ".provenance.json");
			provenance << std::setprecision(17)
								 << "{\"producer\":\"fit2d.C\",\"direct_fit_version\":1,"
								 << "\"pt_low\":" << settingPtLow << ",\"pt_high\":" << settingPtHigh
								 << ",\"y_low\":" << yLow << ",\"y_high\":" << yHigh
								 << ",\"ctau_low\":" << ctRange.first << ",\"ctau_high\":" << ctRange.second
								 << ",\"b_fraction\":" << bFraction.getVal()
								 << ",\"b_fraction_error\":" << bFractionErr << "}\n";
		}
	}
	else
		std::cout << "[PlotOnly] Left saved 2D ROOT file unchanged: " << fitRootName << std::endl;

	// ------------------------------------------------------------------
	// print final result
	// ------------------------------------------------------------------
	cout << "----------------- FIT RESULT FOR THE 2D MODEL ------------" << endl;
	if (model_result)
		model_result->Print("v");
	const TString figMass = figName("mass_fit");
	const TString figLifetime = figName("lifetime_fit");
	if (!zoomOnly) std::cout << "[FIG] fit2d mass fit : " << figMass << std::endl;
	std::cout << "[FIG] fit2d lifetime fit : " << figLifetime << std::endl;
}