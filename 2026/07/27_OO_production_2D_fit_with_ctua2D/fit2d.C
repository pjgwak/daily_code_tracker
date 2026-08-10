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
#include "saved_fit_helpers.h"
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
		if (name == "Nsig" || name == "Nbkg" || name == "bFraction")
			summary.common.push_back(name);
		else if (starts_with(name, "signal_mass_") || starts_with(name, "bkg_mass_"))
			summary.massOnly.push_back(name);
		else if (starts_with(name, "ctau") || starts_with(name, "bkgCtau") ||
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

void fit2d(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4,
					 bool fixResolutionParams = false, bool drawFromSavedFit = false, bool publish = false)
{
	Fit2dScopedMacroTimer timer("fit2d", ptLow, ptHigh, yLow, yHigh);
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
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	std::pair<double, double> ctRange{-6.0, 6.0};
	if (yLow == 0.0f)
		ctRange = {-1.0, 3.0};
	if (yLow >= 1.6f && ptLow < 2.0f && ptHigh <= 2.0f)
		ctRange = {-20.0, 20.0};
	else if (yLow >= 1.6f && std::abs(ptLow - 2.0f) < 1e-3f && std::abs(ptHigh - 3.0f) < 1e-3f)
		ctRange = {-8.0, 10.0};
	else if (yLow >= 1.6f && ptHigh >= 12.0f)
		ctRange = {-2.0, 4.0};

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
	const TString figDir = TString::Format("%s/%s/fit2d", figBaseDir.Data(), yTag.Data());
	const TString fitRootDir = TString::Format("roots/%s/fit2d", yTag.Data());
	const TString fitRootName = TString::Format("%s/fit2d_result_%s.root", fitRootDir.Data(), figTag.Data());

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
	int timePlotBins = std::max(100, obs_time.getBins());
	if (yLow >= 1.6f && ptHigh <= 5.0f)
		timePlotBins = 250;

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
	gSystem->mkdir(fitRootDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	// Saved prefit file paths.  All four files must refer to this same bin tag.
	const TString massFileName = TString::Format("roots/%s/mass/mass_model_%s.root", yTag.Data(), figTag.Data());
	const TString prFileName = TString::Format("roots/%s/ctau_pr/ctau_resolution_%s.root", yTag.Data(), figTag.Data());
	const TString npFileName = TString::Format("roots/%s/ctau_np/ctau_np_model_%s.root", yTag.Data(), figTag.Data());
	const TString bkgFileName = TString::Format("roots/%s/ctau_bkg/ctau_bkg_fitresult_%s.root", yTag.Data(), figTag.Data());

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
	const double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * lifetimeRangeFloor);

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
	const int nResolutionComponents = std::clamp(readIntParam(*prFile, "nResolutionComponents", 1), 1, 4);
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

	RooConstVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", readDoubleParam(*prFile, "ctauTime1Mean", 0.0));
	const double ctauTime1ScaleVal = readDoubleParam(*prFile, "ctauTime1Scale", 1.0);
	const auto ctauTime1ScaleRange = fitRangeAround(ctauTime1ScaleVal, 0.005, 0.20, 1.0, 0.015);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale",
														std::clamp(ctauTime1ScaleVal, ctauTime1ScaleRange.first, ctauTime1ScaleRange.second),
														ctauTime1ScaleRange.first,
														ctauTime1ScaleRange.second);
	RooConstVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", readDoubleParam(*prFile, "ctauTime2Mean", 0.0));
	const double ctauTime2ScaleVal = readDoubleParam(*prFile, "ctauTime2Scale", 1.0);
	const double ctauTime2DeltaVal = std::max(0.001, ctauTime2ScaleVal - ctauTime1ScaleVal);
	const auto ctauTime2DeltaRange = fitRangeAround(ctauTime2DeltaVal, 0.001, 0.30, 2.0, 0.040);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta",
														std::clamp(ctauTime2DeltaVal, ctauTime2DeltaRange.first, ctauTime2DeltaRange.second),
														ctauTime2DeltaRange.first,
														ctauTime2DeltaRange.second);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooConstVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", readDoubleParam(*prFile, "ctauTime3Mean", 0.0));
	const double ctauTime3ScaleVal = readDoubleParam(*prFile, "ctauTime3Scale", 1.0);
	const double ctauTime3DeltaVal = std::max(0.001, ctauTime3ScaleVal - ctauTime2ScaleVal);
	const auto ctauTime3DeltaRange = fitRangeAround(ctauTime3DeltaVal, 0.001, 0.60, 2.0, 0.080);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta",
														std::clamp(ctauTime3DeltaVal, ctauTime3DeltaRange.first, ctauTime3DeltaRange.second),
														ctauTime3DeltaRange.first,
														ctauTime3DeltaRange.second);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooConstVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", readDoubleParam(*prFile, "ctauTime4Mean", 0.0));
	const double ctauTime4ScaleVal = readDoubleParam(*prFile, "ctauTime4Scale", 1.0);
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	const auto ctauTime4DeltaRange = fitRangeAround(ctauTime4DeltaVal, 0.05, 1.20, 1.5, 0.200);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta",
														std::clamp(ctauTime4DeltaVal, ctauTime4DeltaRange.first, ctauTime4DeltaRange.second),
														ctauTime4DeltaRange.first,
														ctauTime4DeltaRange.second);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));

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
	RooConstVar ctauFrac1("ctauFrac1", "ctauFrac1", compFrac1Val);
	RooConstVar ctauFrac2("ctauFrac2", "ctauFrac2", compFrac2Val);
	RooConstVar ctauFrac3("ctauFrac3", "ctauFrac3", compFrac3Val);

	// Build the Gaussian components, then use either one component or their sum.
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale);
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale);
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale);
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale);

	// build ctau resolution model
	std::unique_ptr<RooAddModel> promptResolutionModel;
	RooResolutionModel *promptResolutionPtr = &ctauTime1;
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
		promptResolutionModel = std::make_unique<RooAddModel>("promptResolutionModel", "promptResolutionModel", resolutionPdfList, resolutionFracList);
		promptResolutionPtr = promptResolutionModel.get();
	}
	RooResolutionModel &promptResolution = *promptResolutionPtr;
	RooAbsPdf *promptTimePdf = promptResolutionPtr;

	// ------------------------------------------------------------------
	// ctau prompt resolution - for ctau bkg model (prompt like background)
	// ------------------------------------------------------------------
	const int nBkgResolutionComponents = std::clamp(readIntParam(*bkgFile, "nResolutionComponents", nResolutionComponents), 1, 4);
	RooConstVar bkgCtauTime1Mean("bkgCtauTime1Mean", "bkgCtauTime1Mean", readResultValue(*bkgTimeResult, "ctauTime1Mean", ctauTime1Mean.getVal()));
	double bkgCtauTime1ScaleVal = readResultValue(*bkgTimeResult, "ctauTime1Scale", ctauTime1ScaleVal);
	const auto bkgCtauTime1ScaleRange = fitRangeAround(bkgCtauTime1ScaleVal, 0.005, 0.80, 0.75, 0.030);
	RooRealVar bkgCtauTime1Scale("bkgCtauTime1Scale", "bkgCtauTime1Scale",
														 std::clamp(bkgCtauTime1ScaleVal, bkgCtauTime1ScaleRange.first, bkgCtauTime1ScaleRange.second),
														 bkgCtauTime1ScaleRange.first,
														 bkgCtauTime1ScaleRange.second);
	RooConstVar bkgCtauTime2Mean("bkgCtauTime2Mean", "bkgCtauTime2Mean", readResultValue(*bkgTimeResult, "ctauTime2Mean", ctauTime2Mean.getVal()));
	const double bkgCtauTime2ScaleVal = readResultValue(*bkgTimeResult, "ctauTime2Scale", ctauTime2ScaleVal);
	double bkgCtauTime2DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime2Delta", bkgCtauTime2ScaleVal - bkgCtauTime1ScaleVal));
	const auto bkgCtauTime2DeltaRange = fitRangeAround(bkgCtauTime2DeltaVal, 0.001, 1.00, 1.00, 0.060);
	RooRealVar bkgCtauTime2Delta("bkgCtauTime2Delta", "bkgCtauTime2Delta",
														 std::clamp(bkgCtauTime2DeltaVal, bkgCtauTime2DeltaRange.first, bkgCtauTime2DeltaRange.second),
														 bkgCtauTime2DeltaRange.first,
														 bkgCtauTime2DeltaRange.second);
	RooFormulaVar bkgCtauTime2Scale("bkgCtauTime2Scale", "@0+@1", RooArgList(bkgCtauTime1Scale, bkgCtauTime2Delta));
	RooConstVar bkgCtauTime3Mean("bkgCtauTime3Mean", "bkgCtauTime3Mean", readResultValue(*bkgTimeResult, "ctauTime3Mean", ctauTime3Mean.getVal()));
	const double bkgCtauTime3ScaleVal = readResultValue(*bkgTimeResult, "ctauTime3Scale", ctauTime3ScaleVal);
	double bkgCtauTime3DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime3Delta", bkgCtauTime3ScaleVal - bkgCtauTime2ScaleVal));
	const auto bkgCtauTime3DeltaRange = fitRangeAround(bkgCtauTime3DeltaVal, 0.001, 2.00, 1.00, 0.120);
	RooRealVar bkgCtauTime3Delta("bkgCtauTime3Delta", "bkgCtauTime3Delta",
														 std::clamp(bkgCtauTime3DeltaVal, bkgCtauTime3DeltaRange.first, bkgCtauTime3DeltaRange.second),
														 bkgCtauTime3DeltaRange.first,
														 bkgCtauTime3DeltaRange.second);
	RooFormulaVar bkgCtauTime3Scale("bkgCtauTime3Scale", "@0+@1", RooArgList(bkgCtauTime2Scale, bkgCtauTime3Delta));
	RooConstVar bkgCtauTime4Mean("bkgCtauTime4Mean", "bkgCtauTime4Mean", readResultValue(*bkgTimeResult, "ctauTime4Mean", ctauTime4Mean.getVal()));
	const double bkgCtauTime4ScaleVal = readResultValue(*bkgTimeResult, "ctauTime4Scale", ctauTime4ScaleVal);
	double bkgCtauTime4DeltaVal = std::max(0.001, readResultValue(*bkgTimeResult, "ctauTime4Delta", bkgCtauTime4ScaleVal - bkgCtauTime3ScaleVal));
	const auto bkgCtauTime4DeltaRange = fitRangeAround(bkgCtauTime4DeltaVal, 0.001, 3.00, 1.00, 0.250);
	RooRealVar bkgCtauTime4Delta("bkgCtauTime4Delta", "bkgCtauTime4Delta",
															 std::clamp(bkgCtauTime4DeltaVal, bkgCtauTime4DeltaRange.first, bkgCtauTime4DeltaRange.second),
															 bkgCtauTime4DeltaRange.first,
															 bkgCtauTime4DeltaRange.second);
	RooFormulaVar bkgCtauTime4Scale("bkgCtauTime4Scale", "@0+@1", RooArgList(bkgCtauTime3Scale, bkgCtauTime4Delta));

	// Saved recursive fractions define the Gaussian-mixture weights.
	RooConstVar bkgCtauFrac1("bkgCtauFrac1", "bkgCtauFrac1", readResultValue(*bkgTimeResult, "ctauFrac1", compFrac1Val));
	RooConstVar bkgCtauFrac2("bkgCtauFrac2", "bkgCtauFrac2", readResultValue(*bkgTimeResult, "ctauFrac2", compFrac2Val));
	RooConstVar bkgCtauFrac3("bkgCtauFrac3", "bkgCtauFrac3", readResultValue(*bkgTimeResult, "ctauFrac3", compFrac3Val));

	// Gaussian resolution terms for the prompt-like background component.
	RooGaussModel bkgCtauTime1("bkgCtauTime1", "bkgCtauTime1", obs_time, bkgCtauTime1Mean, bkgCtauTime1Scale);
	RooGaussModel bkgCtauTime2("bkgCtauTime2", "bkgCtauTime2", obs_time, bkgCtauTime2Mean, bkgCtauTime2Scale);
	RooGaussModel bkgCtauTime3("bkgCtauTime3", "bkgCtauTime3", obs_time, bkgCtauTime3Mean, bkgCtauTime3Scale);
	RooGaussModel bkgCtauTime4("bkgCtauTime4", "bkgCtauTime4", obs_time, bkgCtauTime4Mean, bkgCtauTime4Scale);

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


	// fix resolution parameters
	if (fixResolutionParams)
	{
		ctauTime1Scale.setConstant(true);
		ctauTime2Delta.setConstant(true);
		ctauTime3Delta.setConstant(true);
		ctauTime4Delta.setConstant(true);
	}
	if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
	{
		bkgCtauTime1Scale.setConstant(false);
		bkgCtauTime2Delta.setConstant(false);
		bkgCtauTime3Delta.setConstant(false);
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
	if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
	{
		const double signalLogLifetimeHalfWidth = std::max(0.15, 3.0 * signalLogLifetimeFitErr);
		signalLogLifetimeMin = std::max(signalLogLifetimeMin, signalLogLifetimeOffsetVal - signalLogLifetimeHalfWidth);
		signalLogLifetimeMax = std::min(signalLogLifetimeMax, signalLogLifetimeOffsetVal + signalLogLifetimeHalfWidth);
	}
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

	// Let NP tau1 float; retain the saved SS2/SS3 lifetime shapes.
	if (useSignalSS2)
		signal2_log_lifetime_offset.setConstant(true);
	if (useSignalSS3)
		signal3_log_lifetime_offset.setConstant(true);

	// Extended yields imported from the saved ctau_np.C fit.
	const double npYield1 = readDoubleParam(*npFile, "NsignalSS1", 1.0);
	const double npYield2 = readDoubleParam(*npFile, "NsignalSS2", 0.0);
	const double npYield3 = readDoubleParam(*npFile, "NsignalSS3", 0.0);
	RooConstVar NsignalSS1("NsignalSS1", "NsignalSS1", npYield1);
	RooConstVar NsignalSS2("NsignalSS2", "NsignalSS2", npYield2);
	RooConstVar NsignalSS3("NsignalSS3", "NsignalSS3", npYield3);

	// Component PDFs corresponding to the three lifetime parameterizations.
	RooDecay signal_ss1_time("signal_ss1_time", "signal_ss1_time", obs_time, signal_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss2_time("signal_ss2_time", "signal_ss2_time", obs_time, signal2_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss3_time("signal_ss3_time", "signal_ss3_time", obs_time, signal3_lifetime, promptResolution, RooDecay::SingleSided);

	// Assemble the requested number of saved nonprompt lifetime components.
	// Their shapes and extended yields remain fixed; only bFraction floats.
	std::unique_ptr<RooAddPdf> signal_np_owned;
	RooAbsPdf *signal_np_time = nullptr;
	if (nSignalSSComponents == 1)
	{
		signal_np_time = &signal_ss1_time;
	}
	else if (nSignalSSComponents == 2)
	{
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time", RooArgList(signal_ss1_time, signal_ss2_time), RooArgList(NsignalSS1, NsignalSS2));
		signal_np_time = signal_np_owned.get();
	}
	else
	{
		signal_np_owned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time", RooArgList(signal_ss1_time, signal_ss2_time, signal_ss3_time), RooArgList(NsignalSS1, NsignalSS2, NsignalSS3));
		signal_np_time = signal_np_owned.get();
	}

	// Mix the fixed nonprompt template with the prompt template through bFraction.
	double bFractionSeed = 0.10;
	if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
		bFractionSeed = 0.08;
	else if (yLow == 1.6f && ptLow == 2.0f && ptHigh == 3.0f)
		bFractionSeed = 0.11;
	double bFractionMin = 0.06;
	double bFractionMax = 0.6;
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
	if (yLow == 0.0f && ptLow == 9.0f && ptHigh == 10.0f)
	{
		bkgLogLifetimeOffsetVal = -0.77138;
		bkgLogLifetime2OffsetVal = -9.2072;
		bkgSymLogLifetimeOffsetVal = -1.1749;
		bkgFlipLogLifetimeOffsetVal = -4.4021;
		bkgFlipLogLifetime2OffsetVal = -9.1986;
	}
	else if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
	{
		bkgLogLifetimeOffsetVal = -9.2103;
		bkgLogLifetime2OffsetVal = -0.23172;
		bkgLogLifetime3OffsetVal = 2.8166;
		bkgSymLogLifetimeOffsetVal = 1.4301;
		bkgFlipLogLifetimeOffsetVal = -3.0782;
		bkgFlipLogLifetime2OffsetVal = -0.43428;
	}

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
	const double bkgSymLogLifetimeWidth = std::max(2.0 * readResultError(*bkgTimeResult, "bkg_sym_log_lifetime_offset", 0.15), 0.12);
	const double bkgSymLogLifetime2Width = std::max(2.0 * readResultError(*bkgTimeResult, "bkg_sym_log_lifetime2_offset", 0.20), 0.16);
	const double bkgSymLogLifetime3Width = std::max(2.0 * readResultError(*bkgTimeResult, "bkg_sym_log_lifetime3_offset", 0.25), 0.20);

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
	if (yLow == 0.0f && ptLow == 9.0f && ptHigh == 10.0f)
	{
		NplbVal = 306.41;
		Nss1Val = 277.69;
		Nss2Val = 270.05;
		Nds1Val = 106.14;
		Nflip1Val = 132.72;
		Nflip2Val = 50.38;
	}
	else if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
	{
		NplbVal = 168700.0;
		Nss1Val = 35766.0;
		Nss2Val = 13467.0;
		Nds1Val = 6783.3;
		Nflip1Val = 36092.0;
		Nflip2Val = 14918.0;
	}
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
	std::vector<std::unique_ptr<RooRealVar>> bkgCoeffVars;
	RooArgList bkgTimePdfList;
	RooArgList bkgTimeCoeffList;
	RooArgSet bkgCoeffConstraints;
	std::vector<std::unique_ptr<RooConstVar>> bkgCoeffConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> bkgCoeffConstraintPdfs;
	auto addBkgComponent = [&](RooAbsPdf &pdf, const char *name, double yield)
	{
		const double init = std::max(1e-6, yield);
		bkgTimePdfList.add(pdf);
		bkgCoeffVars.push_back(std::make_unique<RooRealVar>(name, name, init, 1e-6, std::max(10.0 * bkgYieldSum, 10.0 * init)));
		bkgTimeCoeffList.add(*bkgCoeffVars.back());
		const double sigma = std::max({0.25 * init, 0.05 * bkgYieldSum, 1.0});
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

	// Background ctau PDF used by bkg_core in the final 2D model.
	std::unique_ptr<RooAddPdf> bkg_time_owned;
	RooAbsPdf *bkg_time_pdf = nullptr;
	if (bkgTimePdfList.getSize() == 1)
	{
		bkg_time_pdf = bkgPromptTimePdf;
	}
	else
	{
		bkg_time_owned = std::make_unique<RooAddPdf>("bkg_time", "bkg_time", bkgTimePdfList, bkgTimeCoeffList);
		bkg_time_pdf = bkg_time_owned.get();
	}

	// ------------------------------------------------------------------
	// Final extended 2D model
	//
	//   model(m, ctau) = Nsig * [S_mass(m) * S_ctau(ctau)]
	//                   + Nbkg * [B_mass(m) * B_ctau(ctau)]
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
	RooProdPdf bkg_core(
			"bkg_core", "bkg_core",
			RooArgSet(*bkg_mass_pdf, *bkg_time_pdf));

	// Extended sum used by fitTo(): yields multiply the corresponding PDFs.
	RooAddPdf model(
			"model", "model",
			RooArgList(signal_core, bkg_core),
			RooArgList(Nsig, Nbkg));

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
	addMassConstraint("Nsig_mass", Nsig, nSigInit,
										std::max(1.0, massErr("Nsig", readDoubleParam(*massFile, "Nsig_err", std::max(1.0, 0.05 * data->numEntries())))));
	addMassConstraint("Nbkg_mass", Nbkg, readDoubleParam(*massFile, "Nbkg", std::max(1.0, 0.5 * data->numEntries())),
										std::max(1.0, massErr("Nbkg", readDoubleParam(*massFile, "Nbkg_err", std::max(1.0, 0.05 * data->numEntries())))));
	// NP tau1 is free in the 2D fit, with the saved 1D uncertainty preventing
	// excursions into the numerically unstable short-lifetime region.
	addMassConstraint("signal_log_lifetime_offset", signal_log_lifetime_offset, signalLogLifetimeOffsetVal, signalLogLifetimeFitErr);
	auto bkgTimeErr = [&](const char *name, double fallback)
	{
		double err = fallback;
		if (bkgTimeResult)
			err = readResultError(*bkgTimeResult, name, err);
		return err;
	};

	// fix mass parameters
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

	// add liftime constraint
	addMassConstraint("ctauTime1Scale", ctauTime1Scale, ctauTime1ScaleVal, 0.010);
	if (nResolutionComponents >= 2)
		addMassConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal, 0.020);
	if (nResolutionComponents >= 3)
		addMassConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal, 0.040);
	if (nResolutionComponents >= 4)
		addMassConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal, 0.080);
	if (yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f)
	{
		addMassConstraint("bkgCtauTime1Scale", bkgCtauTime1Scale, bkgCtauTime1ScaleVal, std::max(0.005, bkgTimeErr("ctauTime1Scale", 0.005)));
		addMassConstraint("bkgCtauTime2Delta", bkgCtauTime2Delta, bkgCtauTime2DeltaVal, std::max(0.010, bkgTimeErr("ctauTime2Delta", 0.010)));
		addMassConstraint("bkgCtauTime3Delta", bkgCtauTime3Delta, bkgCtauTime3DeltaVal, std::max(0.015, bkgTimeErr("ctauTime3Delta", 0.015)));
		addMassConstraint("bkgCtauTime4Delta", bkgCtauTime4Delta, bkgCtauTime4DeltaVal, std::max(0.080, bkgTimeErr("ctauTime4Delta", 0.080)));
	}
	auto addBkgLogOffsetConstraint = [&](const char *name, RooRealVar &var, double fallbackSigma)
	{
		double constraintCenter = readResultValue(*bkgTimeResult, name, var.getVal());
		if ((yLow == 0.0f && ptLow == 9.0f && ptHigh == 10.0f) ||
				(yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f))
			constraintCenter = var.getVal();
		addMassConstraint(name, var, constraintCenter, std::max(fallbackSigma, bkgTimeErr(name, fallbackSigma)));
	};
	addBkgLogOffsetConstraint("bkg_log_lifetime_offset", bkg_log_lifetime_offset, 0.12);
	addBkgLogOffsetConstraint("bkg_log_lifetime2_offset", bkg_log_lifetime2_offset, 0.16);
	addBkgLogOffsetConstraint("bkg_log_lifetime3_offset", bkg_log_lifetime3_offset, 0.20);
	addBkgLogOffsetConstraint("bkg_sym_log_lifetime_offset", bkg_sym_log_lifetime_offset, 0.15);
	addBkgLogOffsetConstraint("bkg_sym_log_lifetime2_offset", bkg_sym_log_lifetime2_offset, 0.20);
	addBkgLogOffsetConstraint("bkg_sym_log_lifetime3_offset", bkg_sym_log_lifetime3_offset, 0.25);
	addBkgLogOffsetConstraint("bkg_flip_log_lifetime_offset", bkg_flip_log_lifetime_offset, 0.12);
	addBkgLogOffsetConstraint("bkg_flip_log_lifetime2_offset", bkg_flip_log_lifetime2_offset, 0.16);
	addBkgLogOffsetConstraint("bkg_flip_log_lifetime3_offset", bkg_flip_log_lifetime3_offset, 0.20);
	for (auto *arg : bkgCoeffConstraints)
		massConstraints.add(*arg);
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

	// ------------------------------------------------------------------
	// fit
	// ------------------------------------------------------------------
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
		apply_saved_fit_result(modelResultHolder.get(), model, RooArgSet(obs_mass, obs_time));
		model_result = modelResultHolder.get();
		std::cout << "[PlotOnly] Loaded saved 2D fit and skipped fit: " << fitRootName << std::endl;
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
												Offset(yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f),
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
												Offset(yLow == 1.6f && ptLow == 1.0f && ptHigh == 2.0f),
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
	if (model_result)
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

	// Propagate the final 2D covariance matrix to the reported nonprompt
	double bFractionErr = -1.0;
	if (model_result)
		bFractionErr = bFraction.getPropagatedError(*model_result);

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
	RooHist *massDataForChi2 = dynamic_cast<RooHist *>(massfitplot->getHist("data"));
	int massDataBins = obs_mass.getBins();
	if (massDataForChi2)
		massDataBins = massDataForChi2->GetN();
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

	cMass->Print(figName("mass_fit"));
	delete cMass;

	// ------------------------------------------------------------------
	// draw - 1D ctau
	// ------------------------------------------------------------------
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
	timePullPlot->Draw();

	// lifetime uses lifetime-only plus common parameters; mass-only parameters are not subtracted.
	RooHist *timeDataForChi2 = dynamic_cast<RooHist *>(timefitplot->getHist("data"));
	int timeDataBins = timePlotBins;
	if (timeDataForChi2)
		timeDataBins = timeDataForChi2->GetN();
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", timeChi2, timeNdf));
	}
	TLine timeLine(ctRange.first, 0.0, ctRange.second, 0.0);
	timeLine.SetLineStyle(2);
	timeLine.Draw("same");
	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

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
			TParameter<double>("bFractionErr", bFractionErr).Write();
			TParameter<double>("Nsig", Nsig.getVal()).Write();
			TParameter<double>("Nbkg", Nbkg.getVal()).Write();
			TParameter<int>("projectionNdfVersion", 2).Write();
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
			TParameter<int>("hesseStatus", hesse).Write();
			fitOut->Write();
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
	std::cout << "[FIG] fit2d mass fit : " << figMass << std::endl;
	std::cout << "[FIG] fit2d lifetime fit : " << figLifetime << std::endl;
}
