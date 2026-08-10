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
#include "RooParametricStepFunction.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TArrayD.h"
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
#include "saved_fit_helpers.h"
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

void ctau_bkg(float ptLow = 2, float ptHigh = 3, float yLow = 1.6, float yHigh = 2.4,
	bool isSilence = true, bool drawFromSavedFit = false, bool publish = false)
{
	ScopedMacroTimer timer("ctau_bkg", ptLow, ptHigh, yLow, yHigh);
	(void)isSilence;
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;

	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	int nSignalSSComponents = 2;
	int nSignalFlipComponents = 2;
	int nSignalDSComponents = 1;
	if (yLow == 0.0f)
	{
		if (ptLow == 20.0f && ptHigh == 30.0f)
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
			nSignalSSComponents = 3; // 0~3
			nSignalFlipComponents = 2; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
		if (ptLow == 2.0f && ptHigh == 3.0f)
		{
			nSignalSSComponents = 3; // 0~3
			nSignalFlipComponents = 2; // 0~3
			nSignalDSComponents = 1; // 0~3
		}
		// if (ptLow == 14.0f && ptHigh == 20.0f)
		// {
		// 	nSignalSSComponents = 1; // 0~3
		// 	nSignalFlipComponents = 0; // 0~3
		// 	nSignalDSComponents = 0; // 0~3
		// }
	}
	nSignalSSComponents = std::clamp(nSignalSSComponents, 0, 3);
	nSignalFlipComponents = std::clamp(nSignalFlipComponents, 0, 3);
	nSignalDSComponents = std::clamp(nSignalDSComponents, 0, 3);

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

	// basic cuts
	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)",
		ptLow, ptHigh, yLow, yHigh);
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	// temporal observable - to decide ctau fit range
	auto *massTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
	auto *timeTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau"));
	if (!massTmp || !timeTmp)
	{
		std::cerr << "ERROR: required variables mass/ctau are missing in dataset." << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// set ctau fit range
	// ------------------------------------------------------------------
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

	const auto ctRangeQuantile = quantileRange(*dataBasic, *timeTmp, 0.0005, 0.9995, false);
	// The lowest forward-pT bin needs the validated wide ctau fit range.
	auto ctRange = ctRangeQuantile;
	if (ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f)
		ctRange = std::make_pair(-20.0, 20.0);

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
	// make sideband dataset
	// ------------------------------------------------------------------
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
	// preapre name tags
	// ------------------------------------------------------------------
	auto formatTag = [](double value)
	{
		return TString::Format("%.2f", value);
	};
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figBaseDir = publish ? "figs_publish" : "figs";
	const TString figDir = TString::Format("%s/%s/ctau_bkg", figBaseDir.Data(), yTag.Data());
	const TString prResultDir = TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString bkgResultDir = TString::Format("roots/%s/ctau_bkg", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(bkgResultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	// ------------------------------------------------------------------
	// connect observable
	// ------------------------------------------------------------------
	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau"));
		obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_time.SetTitle("#font[12]{l}_{J/#psi}");
	obs_time.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	const int timePlotBins = std::max(2, obs_time.getBins());
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", prResultDir.Data(), figTag.Data());
	const TString fitResultFileName = TString::Format("%s/ctau_bkg_fitresult_%s.root", bkgResultDir.Data(), figTag.Data());


	RooDataSet *data = dataSB.get();

		auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	// ------------------------------------------------------------------
	// Load the prompt-resolution calibration saved by ctau_pr.
	// ------------------------------------------------------------------
	// The file contains one required Gaussian core and up to three optional
	// wider cores, selected by nResolutionComponents.
	TFile *resolutionFile = TFile::Open(resolutionFileName);
	if (!resolutionFile || resolutionFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open prompt resolution file: " << resolutionFileName << std::endl;
		return;
	}

	// Read a named double parameter and mark missing values as invalid.
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

	// Check how many Gaussian cores should be reconstructed.
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

	// Core 1 is always present.
	const double ctauMeanScaleVal = readDoubleParam("ctauMeanScale");
	const double ctauTime1MeanVal = readDoubleParam("ctauTime1Mean");
	const double ctauTime1ScaleVal = readDoubleParam("ctauTime1Scale");

	// Later cores and their stored recursive mixture fractions are optional.
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

	// Reject incomplete or non-finite calibrations before building the model.
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

	// RooAddModel needs absolute weights; its final core receives the remainder.
	const double compFrac1Val = (nResolutionComponents == 1) ? 1.0 : ctauFrac1Val;
	const double compFrac2Val = (nResolutionComponents >= 2)
		? ((nResolutionComponents == 2) ? (1.0 - ctauFrac1Val) : (1.0 - ctauFrac1Val) * ctauFrac2Val)
		: 0.0;
	const double compFrac3Val = (nResolutionComponents >= 3)
		? ((nResolutionComponents == 3) ? (1.0 - ctauFrac1Val) * (1.0 - ctauFrac2Val)
			: (1.0 - ctauFrac1Val) * (1.0 - ctauFrac2Val) * ctauFrac3Val)
		: 0.0;

	// set initial value
	const double ctauScale1Seed = std::max(ctauTime1ScaleVal, 1e-3);
	const double ctauScale2Seed = (nResolutionComponents >= 2) ? std::max(ctauTime2ScaleVal, ctauScale1Seed) : ctauScale1Seed;
	const double ctauScale3Seed = (nResolutionComponents >= 3) ? std::max(ctauTime3ScaleVal, ctauScale2Seed) : ctauScale2Seed;
	const double ctauScale4Seed = (nResolutionComponents >= 4) ? std::max(ctauTime4ScaleVal, ctauScale3Seed) : ctauScale3Seed;

	// scale to set conservative lifetime floors.
	double bkgLifetimeFloor = 1e-2;
	double bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 2e-2);
	double bkgSymLifetimeFloor = 1e-2;
	double bkgFlipLifetimeFloor = 5e-3;
	if ((ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f) || nResolutionComponents >= 4)
	{
		const double floorScaleSeed = ctauScale4Seed;
		bkgLifetimeFloor = std::max(1e-2, 0.35 * floorScaleSeed);
		bkgLifetime2Floor = std::max(1.5 * bkgLifetimeFloor, 8e-2);
		bkgSymLifetimeFloor = std::max(1.5 * bkgLifetimeFloor, 5e-2);
		bkgFlipLifetimeFloor = std::max(0.75 * bkgLifetimeFloor, 3e-2);
	}
	const double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	// ------------------------------------------------------------------
	// Build the prompt ctau-resolution parameters imported from ctau_pr.
	//
	// The model contains one to four Gaussian cores. Each later core is
	// expressed as a positive width increment relative to the preceding core,
	// preserving the physical ordering sigma_1 <= sigma_2 <= ... . Imported
	// means and mixture weights are fixed; only widths vary with constraints.
	// ------------------------------------------------------------------
	// Kept for compatibility with the ctau_pr output parameter schema.
	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", ctauMeanScaleVal);

	// Core 1: the narrowest Gaussian resolution component.
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", ctauTime1MeanVal);
	RooRealVar ctauTime1Scale(
		"ctauTime1Scale", "ctauTime1Scale",
		ctauTime1ScaleVal,
		std::max(0.35 * ctauTime1ScaleVal, 1e-3),
		std::max(2.0 * ctauTime1ScaleVal, 5e-3));
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale);

	// The low-pT forward bin uses the validated wider delta limits.
	double ctauTime2DeltaMax = std::max(0.10, 4.0 * std::max(1e-6, ctauTime2ScaleVal - ctauTime1ScaleVal));
	double ctauTime3DeltaMax = std::max(0.20, 4.0 * std::max(1e-6, ctauTime3ScaleVal - ctauTime2ScaleVal));
	double ctauTime4DeltaMax = std::max(0.40, 4.0 * std::max(1e-6, ctauTime4ScaleVal - ctauTime3ScaleVal));
	if (ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f)
	{
		ctauTime2DeltaMax = 0.70;
		ctauTime3DeltaMax = 1.50;
		ctauTime4DeltaMax = 2.00;
	}

	// Core 2: a non-negative increment makes its width no smaller than Core 1.
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", ctauTime2MeanVal);
	RooRealVar ctauTime2Delta(
		"ctauTime2Delta", "ctauTime2Delta",
		std::max(1e-6, ctauTime2ScaleVal - ctauTime1ScaleVal),
		1e-6,
		ctauTime2DeltaMax);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale);
	// Absolute weight of Core 1 in RooAddModel; the final core gets the remainder.
	RooRealVar ctauFrac1("ctauFrac1", "ctauFrac1", compFrac1Val);

	// Core 3
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", ctauTime3MeanVal);
	RooRealVar ctauTime3Delta(
		"ctauTime3Delta", "ctauTime3Delta",
		std::max(1e-6, ctauTime3ScaleVal - ctauTime2ScaleVal),
		1e-6,
		ctauTime3DeltaMax);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale);
	RooRealVar ctauFrac2("ctauFrac2", "ctauFrac2", compFrac2Val);

	// Core 4
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", ctauTime4MeanVal);
	RooRealVar ctauTime4Delta(
		"ctauTime4Delta", "ctauTime4Delta",
		std::max(1e-6, ctauTime4ScaleVal - ctauTime3ScaleVal),
		1e-6,
		ctauTime4DeltaMax);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale);
	RooRealVar ctauFrac3("ctauFrac3", "ctauFrac3", compFrac3Val);

	// Freeze the prompt calibration except for the width parameters constrained below.
	ctauTime1Mean.setConstant(true);
	ctauTime1Scale.setConstant(false);
	ctauTime2Mean.setConstant(true);
	ctauTime2Delta.setConstant(false);
	ctauFrac1.setConstant(true);
	ctauTime3Mean.setConstant(true);
	ctauTime3Delta.setConstant(false);
	ctauFrac2.setConstant(true);
	ctauTime4Mean.setConstant(true);
	ctauTime4Delta.setConstant(false);
	ctauFrac3.setConstant(true);

	// ------------------------------------------------------------------
	// Constrain floating resolution widths near the ctau_pr calibration.
	// ------------------------------------------------------------------
	RooArgSet externalConstraints;

	// Own the fixed Gaussian means, widths, and PDFs passed to fitTo().
	std::vector<std::unique_ptr<RooRealVar>> constraintMeans;
	std::vector<std::unique_ptr<RooRealVar>> constraintSigmas;
	std::vector<std::unique_ptr<RooGaussian>> constraintPdfs;

	// Clamp each seed to its bounds, then add a fixed-centre Gaussian constraint.
	auto addResConstraint = [&](RooRealVar &var,
		const char *name,
		double center,
		double sigma)
	{
		if (!std::isfinite(center) || !std::isfinite(sigma) || sigma <= 0.0)
			return;

		// Keep the initial value strictly inside the allowed fit range.
		const double low = var.getMin();
		const double high = var.getMax();
		const double eps = std::max(1e-9, 1e-4 * (high - low));
		var.setVal(std::min(std::max(center, low + eps), high - eps));

		constraintMeans.push_back(std::make_unique<RooRealVar>(
			Form("%s_constraint_mean", name), Form("%s constraint mean", name), center));
		constraintSigmas.push_back(std::make_unique<RooRealVar>(
			Form("%s_constraint_sigma", name), Form("%s constraint sigma", name),
			sigma, 1e-9, std::max(10.0 * sigma, 1e-8)));
		constraintMeans.back()->setConstant(true);
		constraintSigmas.back()->setConstant(true);

		constraintPdfs.push_back(std::make_unique<RooGaussian>(
			Form("%s_constraint", name), Form("%s constraint", name),
			var, *constraintMeans.back(), *constraintSigmas.back()));
		externalConstraints.add(*constraintPdfs.back());
	};

	// Later cores are constrained through their positive width increments.
	const double ctauTime2DeltaVal = std::max(1e-6, ctauTime2ScaleVal - ctauTime1ScaleVal);
	const double ctauTime3DeltaVal = std::max(1e-6, ctauTime3ScaleVal - ctauTime2ScaleVal);
	const double ctauTime4DeltaVal = std::max(1e-6, ctauTime4ScaleVal - ctauTime3ScaleVal);

	// Use tighter resolution constraints only in the lowest forward-pT bin.
	double ctauTime1ConstraintSigma = std::max(0.012, 0.35 * ctauTime1ScaleVal);
	double ctauTime2ConstraintSigma = std::max(0.010, 0.35 * ctauTime2DeltaVal);
	double ctauTime3ConstraintSigma = std::max(0.020, 0.35 * ctauTime3DeltaVal);
	double ctauTime4ConstraintSigma = std::max(0.040, 0.35 * ctauTime4DeltaVal);
	if (ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f)
	{
		ctauTime1ConstraintSigma = std::max(0.004, 0.12 * ctauTime1ScaleVal);
		ctauTime2ConstraintSigma = std::max(0.010, 0.30 * ctauTime2DeltaVal);
		ctauTime3ConstraintSigma = std::max(0.020, 0.30 * ctauTime3DeltaVal);
		ctauTime4ConstraintSigma = std::max(0.040, 0.30 * ctauTime4DeltaVal);
	}

	// Apply only the constraints for enabled resolution cores.
	addResConstraint(ctauTime1Scale, "ctauTime1Scale", ctauTime1ScaleVal, ctauTime1ConstraintSigma);
	if (nResolutionComponents >= 2)
		addResConstraint(ctauTime2Delta, "ctauTime2Delta", ctauTime2DeltaVal, ctauTime2ConstraintSigma);
	if (nResolutionComponents >= 3)
		addResConstraint(ctauTime3Delta, "ctauTime3Delta", ctauTime3DeltaVal, ctauTime3ConstraintSigma);
	if (nResolutionComponents >= 4)
		addResConstraint(ctauTime4Delta, "ctauTime4Delta", ctauTime4DeltaVal, ctauTime4ConstraintSigma);


	// ------------------------------------------------------------------
	// Assemble the prompt ctau-resolution PDF
	//
	// A single Gaussian core is already a RooResolutionModel. For two to four
	// cores, RooAddModel combines the enabled Gaussian cores with their fixed
	// imported mixture weights. RooAddModel assigns the remaining probability
	// to its final core, so the fraction list has one fewer entry than the PDF
	// list.
	// ------------------------------------------------------------------
	std::unique_ptr<RooAddModel> promptResolutionModel;

	// Use the single core directly unless a multi-core mixture is required.
	RooResolutionModel *time_resolution_ptr = &ctauTime1;
	RooAbsPdf *bkg_time_plb_ptr = &ctauTime1;

	if (nResolutionComponents >= 2)
	{
		// Add the mandatory first two cores and their first absolute weight.
		RooArgList resolutionPdfList;
		RooArgList resolutionFractionList;
		resolutionPdfList.add(ctauTime1);
		resolutionPdfList.add(ctauTime2);
		resolutionFractionList.add(ctauFrac1);

		// Append optional wider cores in increasing-width order.
		if (nResolutionComponents >= 3)
		{
			resolutionPdfList.add(ctauTime3);
			resolutionFractionList.add(ctauFrac2);
		}
		if (nResolutionComponents >= 4)
		{
			resolutionPdfList.add(ctauTime4);
			resolutionFractionList.add(ctauFrac3);
		}

		promptResolutionModel = std::make_unique<RooAddModel>(
			"promptResolutionModel",
			"promptResolutionModel",
			resolutionPdfList,
			resolutionFractionList);

		// The mixture is used both for decay convolution and the prompt PDF term.
		time_resolution_ptr = promptResolutionModel.get();
		bkg_time_plb_ptr = promptResolutionModel.get();
	}

	RooResolutionModel &time_resolution = *time_resolution_ptr;
	RooAbsPdf &bkg_time_plb = *bkg_time_plb_ptr;

	// ------------------------------------------------------------------
	// build ctau bkg model
	// ------------------------------------------------------------------
	// Background-lifetime ranges and special-bin seeds
	// Start from the common configuration, then override only validated bins.
	double bkgLifetimeInit = std::max(0.02, 1.2 * bkgLifetimeFloor);
	double bkgLifetimeCeil = std::max(0.50, maxStableLifetime);
	double bkgLifetime2Init = std::max(0.04, 1.5 * bkgLifetime2Floor);
	double bkgLifetime2Ceil = std::max(0.50, maxStableLifetime);

	double bkgLifetime3Floor = std::max(2.0 * bkgLifetimeFloor, 3e-2);
	if ((ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f) || nResolutionComponents >= 4)
		bkgLifetime3Floor = std::max(2.0 * bkgLifetimeFloor, 0.15);
	double bkgLifetime3Init = std::max(0.10, 1.5 * bkgLifetime3Floor);
	double bkgLifetime3Ceil = std::max(0.80, maxStableLifetime);

	double bkgSymLifetimeInit = std::max(0.04, 1.4 * bkgSymLifetimeFloor);
	double bkgSymLifetimeCeil = std::max(0.30, 0.5 * maxStableLifetime);

	double bkgFlipLifetimeInit = std::max(0.05, 1.2 * bkgFlipLifetimeFloor);
	double bkgFlipLifetimeCeil = std::max(0.50, maxStableLifetime);

	// The lowest forward-pT bin has a separately validated broad lifetime range.
	if (ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f)
	{
		bkgLifetimeInit = 0.18;
		bkgLifetimeCeil = 5.0;
		bkgLifetime2Init = 0.50;
		bkgLifetime2Ceil = 10.0;
		bkgLifetime3Init = 3.0;
		bkgLifetime3Ceil = 80.0;
		bkgSymLifetimeInit = 1.40;
		bkgSymLifetimeCeil = 20.0;
		bkgFlipLifetimeInit = 0.28;
		bkgFlipLifetimeCeil = 8.0;
	}

	// Allow the weak prompt and flipped lifetimes to explore the wider range in the forward bins
	double bkgLifetimeOffsetMin = 3e-3;
	double bkgLifetimeOffsetMax = 3e-2;
	double bkgFlipLifetimeOffsetMin = 1e-2;
	double bkgFlipLifetimeOffsetMax = 1.5e-1;
	if (ptLow < 2.0f && ptHigh <= 2.0f && yLow >= 1.6f)
	{
		bkgLifetimeOffsetMin = 1e-3;
		bkgLifetimeOffsetMax = bkgLifetimeCeil - bkgLifetimeFloor;
		bkgFlipLifetimeOffsetMin = 1e-4;
		bkgFlipLifetimeOffsetMax = bkgFlipLifetimeCeil - bkgFlipLifetimeFloor;
	}
	else if (yLow == 1.6f && ptLow == 4.0f && ptHigh == 5.0f)
	{
		bkgLifetimeOffsetMin = 1e-6;
		bkgLifetimeOffsetMax = bkgLifetimeCeil - bkgLifetimeFloor;
		bkgFlipLifetimeOffsetMin = 1e-4;
		bkgFlipLifetimeOffsetMax = bkgFlipLifetimeCeil - bkgFlipLifetimeFloor;
	}
	else if (yLow == 1.6f &&
		((ptLow == 3.0f && ptHigh == 4.0f) || ptLow == 9.0f || ptLow >= 12.0f))
	{
		bkgLifetimeOffsetMin = 1e-4;
		bkgLifetimeOffsetMax = bkgLifetimeCeil - bkgLifetimeFloor;
		bkgFlipLifetimeOffsetMin = 1e-4;
		bkgFlipLifetimeOffsetMax = bkgFlipLifetimeCeil - bkgFlipLifetimeFloor;
	}

	// Constrain bkg_sym_lifetime = floor + exp(log_offset) away from its floor.
	// The gap bounds prevent an unstable near-zero lifetime and fit-boundary sticking.
	// this is the value to avoid nan error of RooDecay. Got from experience for OO collisions
	const double bkgSymLifetimeMinGap = std::max(0.30 * bkgSymLifetimeFloor, 2e-3);
	const double bkgSymLifetimeMaxGap = std::max(bkgSymLifetimeCeil - bkgSymLifetimeFloor, 2.0 * bkgSymLifetimeMinGap);
	const double bkgSymLifetime2Floor = std::max(2.0 * bkgSymLifetimeFloor, 2e-2);
	const double bkgSymLifetime2Init = std::max(0.07, 1.4 * bkgSymLifetime2Floor);
	const double bkgSymLifetime2Ceil = std::max(0.50, maxStableLifetime);
	const double bkgSymLifetime2MinGap = std::max(0.25 * bkgSymLifetime2Floor, 3e-3);
	const double bkgSymLifetime2MaxGap = std::max(bkgSymLifetime2Ceil - bkgSymLifetime2Floor, 2.0 * bkgSymLifetime2MinGap);
	const double bkgSymLifetime3Floor = std::max(3.0 * bkgSymLifetimeFloor, 4e-2);
	const double bkgSymLifetime3Init = std::max(0.12, 1.3 * bkgSymLifetime3Floor);
	const double bkgSymLifetime3Ceil = std::max(0.80, maxStableLifetime);
	const double bkgSymLifetime3MinGap = std::max(0.20 * bkgSymLifetime3Floor, 4e-3);
	const double bkgSymLifetime3MaxGap = std::max(bkgSymLifetime3Ceil - bkgSymLifetime3Floor, 2.0 * bkgSymLifetime3MinGap);
	const double bkgFlipLifetime2Floor = std::max(1.5 * bkgFlipLifetimeFloor, 1e-2);
	const double bkgFlipLifetime2Init = std::max(0.08, 1.3 * bkgFlipLifetime2Floor);
	const double bkgFlipLifetime2Ceil = std::max(0.60, maxStableLifetime);
	const double bkgFlipLifetime3Floor = std::max(2.5 * bkgFlipLifetimeFloor, 2e-2);
	const double bkgFlipLifetime3Init = std::max(0.12, 1.2 * bkgFlipLifetime3Floor);
	const double bkgFlipLifetime3Ceil = std::max(0.80, maxStableLifetime);

	// Generic lifetime and yield seeds.
	double bkgLogLifetimeSeed = std::log(std::max(bkgLifetimeInit - bkgLifetimeFloor, 1e-4));
	double bkgLogLifetime2Seed = std::log(std::max(bkgLifetime2Init - bkgLifetime2Floor, 1e-4));
	double bkgSymLogLifetimeSeed = std::log(std::max(bkgSymLifetimeInit - bkgSymLifetimeFloor, bkgSymLifetimeMinGap));
	double bkgFlipLogLifetimeSeed = std::log(std::max(bkgFlipLifetimeInit - bkgFlipLifetimeFloor, 1e-4));
	double bkgFlipLogLifetime2Seed = std::log(std::max(bkgFlipLifetime2Init - bkgFlipLifetime2Floor, 1e-4));
	double nplbFraction = 0.60;
	double nss1Fraction = 0.10;
	double nss2Fraction = 0.05;
	double nds1Seed = 0.05 * data->numEntries();
	double nflip1Seed = 0.08 * data->numEntries();
	double nflip2Fraction = 0.05;

	// Use special seeds only for the bins
	if (yLow == 1.6f && ptLow == 3.0f && ptHigh == 4.0f)
	{
		bkgLogLifetimeSeed = -2.7490;
		bkgLogLifetime2Seed = -2.3435;
		bkgSymLogLifetimeSeed = -2.5666;
		bkgFlipLogLifetimeSeed = 1.25;
		bkgFlipLogLifetime2Seed = -3.1420;
		nplbFraction = 0.15;
		nss1Fraction = 0.30;
		nss2Fraction = 0.20;
		nds1Seed = 0.20 * data->numEntries();
		nflip1Seed = 1.0;
		nflip2Fraction = 0.08;
	}
	else if (yLow == 1.6f && ptLow == 5.0f && ptHigh == 6.0f)
	{
		bkgLogLifetimeSeed = -5.80914;
		bkgLogLifetime2Seed = -1.21487;
		bkgSymLogLifetimeSeed = -0.912962;
		bkgFlipLogLifetimeSeed = -2.99405;
		bkgFlipLogLifetime2Seed = -3.95034;
		nplbFraction = 0.7320;
		nss1Fraction = 0.0855;
		nss2Fraction = 0.0454;
		nds1Seed = 0.0579 * data->numEntries();
		nflip1Seed = 0.0477 * data->numEntries();
		nflip2Fraction = 0.0317;
	}
	else if (yLow == 1.6f && ptLow == 9.0f && ptHigh == 10.0f)
	{
		bkgLogLifetimeSeed = -9.20;
		bkgLogLifetime2Seed = -1.103;
		bkgSymLogLifetimeSeed = -1.273;
		bkgFlipLogLifetimeSeed = -9.20;
		bkgFlipLogLifetime2Seed = -1.696;
		nplbFraction = 0.52;
		nss1Fraction = 0.156;
		nss2Fraction = 0.162;
		nds1Seed = 1.0;
		nflip1Seed = 0.107 * data->numEntries();
		nflip2Fraction = 0.057;
	}
	else if (yLow == 1.6f && ptLow == 14.0f && ptHigh == 20.0f)
	{
		bkgLogLifetimeSeed = -9.182;
		bkgLogLifetime2Seed = -1.265;
		bkgSymLogLifetimeSeed = -1.719;
		bkgFlipLogLifetimeSeed = -3.417;
		bkgFlipLogLifetime2Seed = -4.168;
		nplbFraction = 0.291;
		nss1Fraction = 0.166;
		nss2Fraction = 0.396;
		nds1Seed = 1.0;
		nflip1Seed = 0.0695 * data->numEntries();
		nflip2Fraction = 0.0783;
	}

	// parameter - single sided bkg model
	RooRealVar bkg_log_lifetime_offset(
		"bkg_log_lifetime_offset", "log(bkg_lifetime - floor)",
		bkgLogLifetimeSeed,
		std::log(bkgLifetimeOffsetMin),
		std::log(std::max(bkgLifetimeOffsetMax, 2e-3)));
	RooFormulaVar bkg_lifetime(
		"bkg_lifetime",
		Form("%g + exp(@0)", bkgLifetimeFloor),
		RooArgList(bkg_log_lifetime_offset));
	RooRealVar bkg_log_lifetime2_offset(
		"bkg_log_lifetime2_offset", "log(bkg_lifetime2 - floor)",
		bkgLogLifetime2Seed,
		std::log(1e-4),
		std::log(std::max(bkgLifetime2Ceil - bkgLifetime2Floor, 2e-4)));
	RooFormulaVar bkg_lifetime2(
		"bkg_lifetime2",
		Form("%g + exp(@0)", bkgLifetime2Floor),
		RooArgList(bkg_log_lifetime2_offset));
	RooRealVar bkg_log_lifetime3_offset(
		"bkg_log_lifetime3", "log(bkg_lifetime3 - floor)",
		std::log(std::max(bkgLifetime3Init - bkgLifetime3Floor, 1e-4)),
		std::log(1e-4),
		std::log(std::max(bkgLifetime3Ceil - bkgLifetime3Floor, 2e-4)));
	RooFormulaVar bkg_lifetime3(
		"bkg_lifetime3",
		Form("%g + exp(@0)", bkgLifetime3Floor),
		RooArgList(bkg_log_lifetime3_offset));

	// parameter - doubled sided bkg model
	RooRealVar bkg_sym_log_lifetime_offset(
		"bkg_sym_log_lifetime_offset", "log(bkg_sym_lifetime - floor)",
		bkgSymLogLifetimeSeed,
		std::log(bkgSymLifetimeMinGap),
		std::log(bkgSymLifetimeMaxGap));
	RooFormulaVar bkg_sym_lifetime(
		"bkg_sym_lifetime",
		Form("%g + exp(@0)", bkgSymLifetimeFloor),
		RooArgList(bkg_sym_log_lifetime_offset));
	RooRealVar bkg_sym_log_lifetime2_offset(
		"bkg_sym_log_lifetime2_offset", "log(bkg_sym_lifetime2 - floor)",
		std::log(std::max(bkgSymLifetime2Init - bkgSymLifetime2Floor, bkgSymLifetime2MinGap)),
		std::log(bkgSymLifetime2MinGap),
		std::log(bkgSymLifetime2MaxGap));
	RooFormulaVar bkg_sym_lifetime2(
		"bkg_sym_lifetime2",
		Form("%g + exp(@0)", bkgSymLifetime2Floor),
		RooArgList(bkg_sym_log_lifetime2_offset));
	RooRealVar bkg_sym_log_lifetime3_offset(
		"bkg_sym_log_lifetime3_offset", "log(bkg_sym_lifetime3 - floor)",
		std::log(std::max(bkgSymLifetime3Init - bkgSymLifetime3Floor, bkgSymLifetime3MinGap)),
		std::log(bkgSymLifetime3MinGap),
		std::log(bkgSymLifetime3MaxGap));
	RooFormulaVar bkg_sym_lifetime3(
		"bkg_sym_lifetime3",
		Form("%g + exp(@0)", bkgSymLifetime3Floor),
		RooArgList(bkg_sym_log_lifetime3_offset));

	// parameter - flipped bkg model
	RooRealVar bkg_flip_log_lifetime_offset(
		"bkg_flip_log_lifetime_offset", "log(bkg_flip_lifetime - floor)",
		bkgFlipLogLifetimeSeed,
		std::log(bkgFlipLifetimeOffsetMin),
		std::log(std::max(bkgFlipLifetimeOffsetMax, 2e-4)));
	RooFormulaVar bkg_flip_lifetime(
		"bkg_flip_lifetime",
		Form("%g + exp(@0)", bkgFlipLifetimeFloor),
		RooArgList(bkg_flip_log_lifetime_offset));
	RooRealVar bkg_flip_log_lifetime2_offset(
		"bkg_flip_log_lifetime2_offset", "log(bkg_flip_lifetime2 - floor)",
		bkgFlipLogLifetime2Seed,
		std::log(5e-3),
		std::log(1.0));
	RooFormulaVar bkg_flip_lifetime2(
		"bkg_flip_lifetime2",
		Form("%g + exp(@0)", bkgFlipLifetime2Floor),
		RooArgList(bkg_flip_log_lifetime2_offset));
	RooRealVar bkg_flip_log_lifetime3_offset(
		"bkg_flip_log_lifetime3_offset", "log(bkg_flip_lifetime3 - floor)",
		std::log(std::max(bkgFlipLifetime3Init - bkgFlipLifetime3Floor, 1e-4)),
		std::log(1e-4),
		std::log(std::max(bkgFlipLifetime3Ceil - bkgFlipLifetime3Floor, 2e-4)));
	RooFormulaVar bkg_flip_lifetime3(
		"bkg_flip_lifetime3",
		Form("%g + exp(@0)", bkgFlipLifetime3Floor),
		RooArgList(bkg_flip_log_lifetime3_offset));

	// declaer bkg model
	RooDecay bkg_time_ss("bkg_time_ss", "bkg_time_ss", obs_time, bkg_lifetime, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss2("bkg_time_ss2", "bkg_time_ss2", obs_time, bkg_lifetime2, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ss3("bkg_time_ss3", "bkg_time_ss3", obs_time, bkg_lifetime3, time_resolution, RooDecay::SingleSided);
	RooDecay bkg_time_ds("bkg_time_ds", "bkg_time_ds", obs_time, bkg_sym_lifetime, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds2("bkg_time_ds2", "bkg_time_ds2", obs_time, bkg_sym_lifetime2, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_ds3("bkg_time_ds3", "bkg_time_ds3", obs_time, bkg_sym_lifetime3, time_resolution, RooDecay::DoubleSided);
	RooDecay bkg_time_flip("bkg_time_flip", "bkg_time_flip", obs_time, bkg_flip_lifetime, time_resolution, RooDecay::Flipped);
	RooDecay bkg_time_flip2("bkg_time_flip2", "bkg_time_flip2", obs_time, bkg_flip_lifetime2, time_resolution, RooDecay::Flipped);
	RooDecay bkg_time_flip3("bkg_time_flip3", "bkg_time_flip3", obs_time, bkg_flip_lifetime3, time_resolution, RooDecay::Flipped);

	// fraction of bkg models - used Nevt to always make them positive
	const double bkgYieldMax = 1.20 * std::max(1, data->numEntries());
	RooRealVar Nplb("Nplb", "Nplb", nplbFraction * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nss1("Nss1", "Nss1", nss1Fraction * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nss2("Nss2", "Nss2", nss2Fraction * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nss3("Nss3", "Nss3", 0.05 * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nds1("Nds1", "Nds1", nds1Seed, 0.0, bkgYieldMax);
	RooRealVar Nds2("Nds2", "Nds2", 0.02 * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nds3("Nds3", "Nds3", 0.02 * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nflip1("Nflip1", "Nflip1", nflip1Seed, 0.0, bkgYieldMax);
	RooRealVar Nflip2("Nflip2", "Nflip2", nflip2Fraction * data->numEntries(), 0.0, bkgYieldMax);
	RooRealVar Nflip3("Nflip3", "Nflip3", 0.02 * data->numEntries(), 0.0, bkgYieldMax);

	// list PDFs and fractions to consider possible combinations
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

	// declare bkg model
	RooAddPdf bkg_time("bkg_time", "bkg_time", bkgTimePdfList, bkgTimeYieldList);
	RooAbsPdf &time_pdf = bkg_time;

	// ------------------------------------------------------------------
	// Fit or restore the background ctau model.
	// ------------------------------------------------------------------
	std::unique_ptr<TFile> savedBkgFitFile;
	std::unique_ptr<RooFitResult> timeResultHolder;

	// helper to keep the fit options identical for the baseline fit and its retry
	// only change the strategy
	auto fitTimeModel = [&](int strategy)
	{
		return std::unique_ptr<RooFitResult>(
			time_pdf.fitTo(
				*data,
				Extended(),
				Save(),
				SumW2Error(isWeight),
				Strategy(strategy),
				Offset(true),
				RecoverFromUndefinedRegions(1.0),
				ExternalConstraints(externalConstraints)));
	};

	if (drawFromSavedFit)
	{
		// publish mode: restore the saved parameters without fitting again.
		if (!load_saved_fit_file(savedBkgFitFile, fitResultFileName, "background ctau"))
			return;
		timeResultHolder = clone_saved_fit_result(savedBkgFitFile.get(), "timeResult");
		if (!timeResultHolder)
		{
			std::cerr << "ERROR: timeResult not found in saved background ctau file: " << fitResultFileName << std::endl;
			return;
		}
		apply_saved_fit_result(timeResultHolder.get(), time_pdf, RooArgSet(obs_time));
		std::cout << "[PlotOnly] Loaded saved background ctau fit and skipped fit: " << fitResultFileName << std::endl;
	}
	else
	{
		// Start with the nominal strategy;
		timeResultHolder = fitTimeModel(1);

		// Retry once if convergence is insufficient.
		if (timeResultHolder &&
			(timeResultHolder->status() != 0 || timeResultHolder->covQual() < 3 || timeResultHolder->edm() > 1e-3))
		{
			timeResultHolder = fitTimeModel(2);
		}
	}

	RooFitResult *time_result = timeResultHolder.get();
	if (time_result)
		time_result->Print();

	// Save the completed fit session before creating plots.
	if (!drawFromSavedFit)
	{
		TFile fitResultFile(fitResultFileName, "RECREATE");
		if (!fitResultFile.IsZombie())
		{
			if (time_result)
				time_result->Write("timeResult");

			// Store the model configuration needed to interpret the saved fit.
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
			TParameter<int>("nSignalFlipComponents", nSignalFlipComponents).Write();
			TParameter<int>("nSignalDSComponents", nSignalDSComponents).Write();
			TParameter<double>("bkgLifetimeFloor", bkgLifetimeFloor).Write();
			TParameter<double>("bkgLifetime2Floor", bkgLifetime2Floor).Write();
			TParameter<double>("bkgLifetime3Floor", bkgLifetime3Floor).Write();
			TParameter<double>("bkgSymLifetimeFloor", bkgSymLifetimeFloor).Write();
			TParameter<double>("bkgSymLifetime2Floor", bkgSymLifetime2Floor).Write();
			TParameter<double>("bkgSymLifetime3Floor", bkgSymLifetime3Floor).Write();
			TParameter<double>("bkgFlipLifetimeFloor", bkgFlipLifetimeFloor).Write();
			TParameter<double>("bkgFlipLifetime2Floor", bkgFlipLifetime2Floor).Write();
			TParameter<double>("bkgFlipLifetime3Floor", bkgFlipLifetime3Floor).Write();
			TParameter<double>("ctRangeMin", ctRange.first).Write();
			TParameter<double>("ctRangeMax", ctRange.second).Write();
			fitResultFile.Write();

			std::cout << "Saved ctau background fit results to " << fitResultFileName << std::endl;
		}
		else
		{
			std::cerr << "ERROR: cannot create output fit result file: " << fitResultFileName << std::endl;
		}
	}
	else
	{
		std::cout << "[PlotOnly] Left saved background ctau ROOT file unchanged: " << fitResultFileName << std::endl;
	}

	// PDF objects for plotting
	RooAbsPdf &time_pdf_ss1 = bkg_time_ss;
	RooAbsPdf &time_pdf_plb = bkg_time_plb;
	RooAbsPdf &time_pdf_ss2 = bkg_time_ss2;
	RooAbsPdf &time_pdf_ss3 = bkg_time_ss3;
	RooAbsPdf &time_pdf_ds = bkg_time_ds;
	RooAbsPdf &time_pdf_ds2 = bkg_time_ds2;
	RooAbsPdf &time_pdf_ds3 = bkg_time_ds3;
	RooAbsPdf &time_pdf_flip = bkg_time_flip;
	RooAbsPdf &time_pdf_flip2 = bkg_time_flip2;
	RooAbsPdf &time_pdf_flip3 = bkg_time_flip3;

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
	if (!publish)
	{
		TLatex tx;
		tx.SetNDC();
		tx.SetTextSize(0.03);
		tx.SetTextFont(42);
		int minimize = -1;
		int hesse = -1;
		if (time_result)
		{
			for (UInt_t i = 0, n = time_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = time_result->statusLabelHistory(i);
				if (lab)
				{
					if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
						minimize = time_result->statusCodeHistory(i);
					else if (!strcmp(lab, "HESSE"))
						hesse = time_result->statusCodeHistory(i);
				}
			}
		}
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
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
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
			}
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
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
			printVar("resScale_{1}", ctauTime1Scale);
		if (nSignalSSComponents >= 1)
			printVar("#tau_{SS1}", bkg_lifetime);
		if (nSignalSSComponents >= 2)
			printVar("#tau_{SS2}", bkg_lifetime2);
		if (nSignalSSComponents >= 3)
			printVar("#tau_{SS3}", bkg_lifetime3);
		if (nSignalDSComponents >= 1)
			printVar("#tau_{DS1}", bkg_sym_lifetime);
		if (nSignalDSComponents >= 2)
			printVar("#tau_{DS2}", bkg_sym_lifetime2);
		if (nSignalDSComponents >= 3)
			printVar("#tau_{DS3}", bkg_sym_lifetime3);
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

	const int npar = time_result ? time_result->floatParsFinal().getSize() : 0;
	RooHist *hdataForChi2 = dynamic_cast<RooHist *>(timefitplot->getHist("data"));
	const int ndataBins = hdataForChi2 ? hdataForChi2->GetN() : timePlotBins;
	const int ndf = std::max(1, ndataBins - npar);
	const double chi2ndf = timefitplot->chiSquare("model", "data", npar);
	const double chi2 = chi2ndf * ndf;
	const double pvalue = TMath::Prob(chi2, ndf);
	std::cout << Form("Background ctau chi2/ndf: %.1f/%d = %.3f, p=%.3g",
		chi2, ndf, chi2ndf, pvalue) << std::endl;
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", chi2, ndf));
	}

	TLine line(ctRange.first, 0.0, ctRange.second, 0.0);
	line.SetLineStyle(2);
	line.Draw("same");

	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;



	// ------------------------------------------------------------------
	// print final result
	// ------------------------------------------------------------------
	cout << "------------------ FIT RESULT SUMMARY --------------------" << endl;
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	time_result->Print("v");
	cout << "Prompt resolution components used in BKG fit: " << nResolutionComponents << endl;
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_bkg lifetime fit : " << figLifetime << std::endl;
}
