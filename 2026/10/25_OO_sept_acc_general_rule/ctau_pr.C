#include "utils/prefit_jitter.h"
#include "utils/robust_chi2.h"
#include "utils/projection_diagnostics.h"
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
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
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
#include "TMath.h"
#include "utils/saved_fit_helpers.h"
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

void ctau_pr(float ptLow = 2, float ptHigh = 3, float yLow = 1.6, float yHigh = 2.4,
	bool drawFromSavedFit = false, bool publish = false)
{
	ScopedMacroTimer timer("ctau_pr", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;

	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	int nResolutionComponents = 3; // 1~4 components

	if (ptHigh <= 14.0f)
	{
		nResolutionComponents = 4;
	}
	// Study hook: OO_PR_NRES overrides the prompt Gaussian count (prefit model-selection scan).
	if (gSystem->Getenv("OO_PR_NRES")) nResolutionComponents = std::clamp(std::atoi(gSystem->Getenv("OO_PR_NRES")), 1, 4);
	
	if (yLow == 1.6f)
	{
		if (ptLow == 1.0f && ptHigh == 2.0f) // dummy
		{
			// nResolutionComponents = 4;
		}
	}
	nResolutionComponents = std::clamp(nResolutionComponents, 1, 4);


	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	const char *DATA_ROOT = "inputs/roodataset_roots/RooDataSet_miniAOD_isMC1_PR_Hijing_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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

	// 24 study hook: OO_PR_QUANT="qlo,qhi" changes the MC quantile window of the prompt prefit.
	double prQuantLo = 0.0005, prQuantHi = 0.9995;
	if (const char *q = gSystem->Getenv("OO_PR_QUANT")) sscanf(q, "%lf,%lf", &prQuantLo, &prQuantHi);
	const auto ctRange = quantileRange(*dataBasic, *timeTmp, prQuantLo, prQuantHi, false);

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
	// preapre name tags
	// ------------------------------------------------------------------
	auto formatTag = [](double value)
	{
		return TString::Format("%.2f", value);
	};
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	// Study hook: OO_PREFIT_OUT=<dir> writes this prefit (ROOT + figures) under <dir>/ instead of the default roots/, figs/.
	const TString prefitOut = gSystem->Getenv("OO_PREFIT_OUT") ? TString(gSystem->Getenv("OO_PREFIT_OUT")) + "/" : TString("");
	const TString figBaseDir = prefitOut + (publish ? "figs_publish" : "figs");
	const TString figDir = TString::Format("%s/%s/ctau_pr", figBaseDir.Data(), yTag.Data());
	const TString resultDir = prefitOut + TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", resultDir.Data(), figTag.Data());
	
	std::unique_ptr<TFile> savedResolutionFile;
	if (drawFromSavedFit)
	{
		if (!load_saved_fit_file(savedResolutionFile, resolutionFileName, "prompt ctau"))
			return;
		nResolutionComponents = std::clamp(read_saved_int_param(savedResolutionFile.get(), "nResolutionComponents", nResolutionComponents), 1, 4);
	}
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("include/rootlogon.C");

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
	const int massPlotBins = std::max(2, obs_mass.getBins() * 2);
	const int timePlotBins = std::max(2, obs_time.getBins());

	RooDataSet *data = dataSel.get();

	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	// ------------------------------------------------------------------
	// build ctau resolution model
	// ------------------------------------------------------------------
	// Prompt J/#psi ctau resolution is described by up to four Gaussian cores.
	// All cores share one fitted mean. Their widths are constructed as a
	// positive cumulative sum, ensuring sigma_1 < sigma_2 < sigma_3 < sigma_4.
	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", 1.0);

	// initial values
	double time1MeanSeed = 0.0;
	double time1ScaleSeed = 0.015;
	double time2DeltaSeed = 0.010;
	double time3DeltaSeed = 0.010;
	double time4DeltaSeed = 0.10;
	double frac1Seed = 0.65;
	double frac2Seed = 0.70;
	double frac3Seed = 0.70;

	// Bin-specific seeds from stable mid-rapidity fits.
	if (yLow == 0.0f && yHigh == 1.6f)
	{
		if (ptLow == 16.0f && ptHigh == 20.0f)
		{
			time1MeanSeed = 1.2739e-05;
			time1ScaleSeed = 1.2069e-02;
			time2DeltaSeed = 0.0137335;
			time3DeltaSeed = 0.0395093;
			time4DeltaSeed = 0.142132;
			frac1Seed = 6.1895e-02;
			frac2Seed = 6.9552e-01;
			frac3Seed = 9.5632e-01;
		}
	}

	// Bin-specific seeds from stable forward-rapidity fits.
	// fwd
	if (yLow == 1.6f && yHigh == 2.4f)
	{
		if (ptLow == 1.0f && ptHigh == 2.0f)
		{
			time1MeanSeed = 6.226e-4;
			time1ScaleSeed = 0.06851;
			time2DeltaSeed = 0.06103;
			time3DeltaSeed = 0.1281;
			time4DeltaSeed = 0.4865;
			frac1Seed = 0.4036;
			frac2Seed = 0.7600;
			frac3Seed = 0.8707;
		}
		else if (ptLow == 2.0f && ptHigh == 3.0f)
		{
			time1MeanSeed = 7.073e-5;
			time1ScaleSeed = 0.04538;
			time2DeltaSeed = 0.03539;
			time3DeltaSeed = 0.07643;
			time4DeltaSeed = 0.3227;
			frac1Seed = 0.4341;
			frac2Seed = 0.7762;
			frac3Seed = 0.8701;
		}
		else if (ptLow == 3.0f && ptHigh == 4.0f)
		{
			time1MeanSeed = 2.08e-4;
			time1ScaleSeed = 0.0338;
			time2DeltaSeed = 0.0277;
			time3DeltaSeed = 0.0611;
			time4DeltaSeed = 0.243;
			frac1Seed = 0.425;
			frac2Seed = 0.802;
			frac3Seed = 0.858;
		}
	}
	

	// Core 1: common mean and the narrowest resolution width.
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", time1MeanSeed, -0.01, 0.01);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale", time1ScaleSeed, 0.005, 0.20);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale);

	// Core 2: same mean, with a positive increment over the Core-1 width.
	RooFormulaVar ctauTime2Mean("ctauTime2Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", time2DeltaSeed, 1e-4, 0.30);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale);

	// Core 3: same mean, with a further positive width increment.
	RooFormulaVar ctauTime3Mean("ctauTime3Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", time3DeltaSeed, 0.003, 0.20);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale);

	// Core 4: broadest optional tail core, also sharing the common mean.
	RooFormulaVar ctauTime4Mean("ctauTime4Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", time4DeltaSeed, 0.05, 1.50);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale);

	// fractions
	RooRealVar ctauFrac1Var("ctauFrac1Var", "ctauFrac1Var", frac1Seed, 0.05, 0.98);
	RooRealVar ctauFrac2Var("ctauFrac2Var", "ctauFrac2Var", frac2Seed, 0.001, 0.98);
	RooRealVar ctauFrac3Var("ctauFrac3Var", "ctauFrac3Var", frac3Seed, 0.01, 0.99);

	// fill the list of PDF and fractions
	RooArgList promptTimePdfList;
	RooArgList promptTimeFracList;
	promptTimePdfList.add(ctauTime1);

	if (nResolutionComponents >= 2)
	{
		promptTimePdfList.add(ctauTime2);
		promptTimeFracList.add(ctauFrac1Var);
	}
	else
	{
		ctauFrac1Var.setVal(0.0);
	}

	if (nResolutionComponents >= 3)
	{
		promptTimePdfList.add(ctauTime3);
		promptTimeFracList.add(ctauFrac2Var);
	}
	else
	{
		ctauFrac2Var.setVal(0.0);
	}

	if (nResolutionComponents >= 4)
	{
		promptTimePdfList.add(ctauTime4);
		promptTimeFracList.add(ctauFrac3Var);
	}
	else
	{
		ctauFrac3Var.setVal(0.0);
	}

	// ctau prompt resolution model
	RooAddPdf time_pdf("time_pdf", "time_pdf", promptTimePdfList, promptTimeFracList, true);
	
	
	// ------------------------------------------------------------------
	// fit prompt ctau resolution
	// ------------------------------------------------------------------
	// In plot-only mode, restore the saved resolution parameters and skip the
	// minimization. Otherwise, fit the selected Gaussian-mixture model to data.
	std::unique_ptr<RooFitResult> timeResultHolder;
	RooFitResult *time_result = nullptr;

	if (drawFromSavedFit)
	{
		// The saved file stores absolute widths, while this model fits positive
		// width increments. Convert each saved width back to its increment.
		ctauTime1Mean.setVal(read_saved_double_param(savedResolutionFile.get(), "ctauTime1Mean", ctauTime1Mean.getVal()));
		ctauTime1Scale.setVal(read_saved_double_param(savedResolutionFile.get(), "ctauTime1Scale", ctauTime1Scale.getVal()));

		if (nResolutionComponents >= 2)
		{
			const double savedScale = read_saved_double_param(savedResolutionFile.get(), "ctauTime2Scale", ctauTime2Scale.getVal());
			ctauTime2Delta.setVal(std::max(1e-6, savedScale - ctauTime1Scale.getVal()));
			ctauFrac1Var.setVal(read_saved_double_param(savedResolutionFile.get(), "ctauFrac1", ctauFrac1Var.getVal()));
		}

		if (nResolutionComponents >= 3)
		{
			const double savedScale = read_saved_double_param(savedResolutionFile.get(), "ctauTime3Scale", ctauTime3Scale.getVal());
			ctauTime3Delta.setVal(std::max(1e-6, savedScale - ctauTime2Scale.getVal()));
			ctauFrac2Var.setVal(read_saved_double_param(savedResolutionFile.get(), "ctauFrac2", ctauFrac2Var.getVal()));
		}

		if (nResolutionComponents >= 4)
		{
			const double savedScale = read_saved_double_param(savedResolutionFile.get(), "ctauTime4Scale", ctauTime4Scale.getVal());
			ctauTime4Delta.setVal(std::max(1e-6, savedScale - ctauTime3Scale.getVal()));
			ctauFrac3Var.setVal(read_saved_double_param(savedResolutionFile.get(), "ctauFrac3", ctauFrac3Var.getVal()));
		}

		std::cout << "[PlotOnly] Loaded saved prompt ctau parameters and skipped fit: "
		          << resolutionFileName << std::endl;
	}
	else
	{
		// Default optimizer setup: absolute NLL (Offset=false) with Strategy 2.
		// The two validated exceptional bins override only the required setting.
		bool firstFitOffset = false;
		int firstFitStrategy = 2;
		bool retryFitOffset = true;

		if (yLow == 0.0f && ptLow == 12.0f && ptHigh == 14.0f)
		{
			firstFitOffset = true;
			retryFitOffset = false;
		}
		else if (yLow == 1.6f && ptLow == 5.0f && ptHigh == 6.0f)
		{
			firstFitStrategy = 1;
			retryFitOffset = false;
		}

		oo_prefit_jitter(time_pdf, RooArgSet(obs_time));
		// First fit uses the bin-specific setup chosen above.
		timeResultHolder.reset(
			time_pdf.fitTo(*data, Save(), Minimizer("Minuit2", "migrad"), SumW2Error(isWeight), Offset(firstFitOffset), Strategy(firstFitStrategy)));

		// Retry once only when the minimizer or covariance quality is insufficient.
		const bool fitNeedsRetry = timeResultHolder &&
			(timeResultHolder->status() != 0 || timeResultHolder->covQual() < 3 || timeResultHolder->edm() > 1e-3);
		if (fitNeedsRetry)
		{
			timeResultHolder.reset(
				time_pdf.fitTo(*data, Save(), Minimizer("Minuit2", "migrad"), SumW2Error(isWeight), Offset(retryFitOffset), Strategy(2)));
		}

		// 2026-10-15 prefit robustness: if both bin-specific attempts fail, try the remaining
		// (Strategy, Offset) combinations in a fixed order and keep the first valid result.
		auto validPrefit = [](const RooFitResult *r) { return r && r->status() == 0 && r->covQual() == 3 && r->edm() <= 1e-3; };
		for (int strategy : {2, 1})
			for (bool offset : {true, false})
			{
				if (validPrefit(timeResultHolder.get())) break;
				timeResultHolder.reset(time_pdf.fitTo(*data, Save(), Minimizer("Minuit2", "migrad"), SumW2Error(isWeight), Offset(offset), Strategy(strategy)));
				std::cout << "[PREFIT_ROBUST] strategy=" << strategy << " offset=" << offset << " status=" << timeResultHolder->status()
				          << " covQual=" << timeResultHolder->covQual() << std::endl;
			}

		time_result = timeResultHolder.get();
		if (time_result)
			time_result->Print();
	}

	// for log: ctauFrac is recursive fraction from fit
	const double ctauFrac1 = (nResolutionComponents >= 2) ? ctauFrac1Var.getVal() : 1.0;
	const double ctauFrac2 = (nResolutionComponents >= 3) ? ctauFrac2Var.getVal() : 0.0;
	const double ctauFrac3 = (nResolutionComponents >= 4) ? ctauFrac3Var.getVal() : 0.0;

	// for log: transform recursive fraction to simple fractions of components
	// nEntries * compFrac = number of events in each PDF
	const double compFrac1 = (nResolutionComponents == 1) ? 1.0 : ctauFrac1;
	const double compFrac2 = (nResolutionComponents >= 2) ? ((nResolutionComponents == 2) ? (1.0 - ctauFrac1) : (1.0 - ctauFrac1) * ctauFrac2) : 0.0;
	const double compFrac3 = (nResolutionComponents >= 3) ? ((nResolutionComponents == 3) ? (1.0 - ctauFrac1) * (1.0 - ctauFrac2) : (1.0 - ctauFrac1) * (1.0 - ctauFrac2) * ctauFrac3) : 0.0;
	const double compFrac4 = (nResolutionComponents >= 4) ? ((1.0 - ctauFrac1) * (1.0 - ctauFrac2) * (1.0 - ctauFrac3)) : 0.0;
	const double nEntries = data->numEntries();

	// save outputs
	if (!drawFromSavedFit)
	{
		TFile resolutionFile(resolutionFileName, "RECREATE");
		if (!resolutionFile.IsZombie())
		{
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			// 24 study: offset-free NLL at the saved minimum, comparable across starts.
			if (!drawFromSavedFit)
			{
				std::unique_ptr<RooAbsReal> rawNll(time_pdf.createNLL(*data, Offset(false)));
				TParameter<double>("rawPrefitNll", rawNll->getVal()).Write();
			}
			TParameter<double>("prCalibrationLow", ctRange.first).Write();
			TParameter<double>("prCalibrationHigh", ctRange.second).Write();
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
			if (time_result)
				time_result->Write("timeResult");
			resolutionFile.Write();
		}
		else
		{
			std::cerr << "ERROR: cannot create resolution file: " << resolutionFileName << std::endl;
		}
		std::cout << "Saved ctau prompt resolution parameters to " << resolutionFileName << std::endl;
	}

	RooAbsPdf &time_pdf_component1 = ctauTime1;
	RooAbsPdf &time_pdf_component2 = ctauTime2;
	RooAbsPdf &time_pdf_component3 = ctauTime3;
	RooAbsPdf &time_pdf_component4 = ctauTime4;

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
	time_pdf.plotOn(timefitplot, Precision(1e-6), LineColor(kBlack), LineWidth(2), Name("model"));
	time_pdf_component1.plotOn(
			timefitplot, Precision(1e-6),
			LineColor(kBlue + 1),
			LineStyle(kDashed),
			LineWidth(2),
			Normalization(compFrac1 * nEntries, RooAbsReal::NumEvent),
			Name("ctau_component_1"));
	if (nResolutionComponents >= 2)
	{
		time_pdf_component2.plotOn(
				timefitplot, Precision(1e-6),
				LineColor(kRed + 1),
				LineStyle(kDashed),
				LineWidth(2),
				Normalization(compFrac2 * nEntries, RooAbsReal::NumEvent),
				Name("ctau_component_2"));
	}
		if (nResolutionComponents >= 3)
		{
			time_pdf_component3.plotOn(
					timefitplot, Precision(1e-6),
					LineColor(kMagenta + 1),
					LineStyle(kDashed),
					LineWidth(2),
					Normalization(compFrac3 * nEntries, RooAbsReal::NumEvent),
					Name("ctau_component_3"));
		}
		if (nResolutionComponents >= 4)
		{
			time_pdf_component4.plotOn(
					timefitplot, Precision(1e-6),
					LineColor(kOrange + 7),
					LineStyle(kDashed),
					LineWidth(2),
					Normalization(compFrac4 * nEntries, RooAbsReal::NumEvent),
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
			if (nResolutionComponents >= 2)
				printVar("f_{1}", ctauFrac1Var);
			if (nResolutionComponents >= 3)
				printVar("f_{2}", ctauFrac2Var);
				printVar("mean", ctauTime1Mean);
				printVar("#sigma_{1}", ctauTime1Scale);
				if (nResolutionComponents >= 2)
				{
					printVar("#Delta #sigma_{21}", ctauTime2Delta);
				}
				if (nResolutionComponents >= 4)
				{
					printVar("f_{3}", ctauFrac3Var);
				}
				if (nResolutionComponents >= 3)
				{
					printVar("#Delta #sigma_{32}", ctauTime3Delta);
				}
				if (nResolutionComponents >= 4)
				{
					printVar("#Delta #sigma_{43}", ctauTime4Delta);
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
	const int npar = time_result ? time_result->floatParsFinal().getSize() : saved_prefit_npar(resolutionFileName, "timeResult");
	RooHist *hdataForChi2 = dynamic_cast<RooHist *>(timefitplot->getHist("data"));
	const int ndataBins = prefit_chi2_bins(*timefitplot);
	const int ndf = ndataBins - npar;
	const double chi2ndf = timefitplot->chiSquare("model", "data", npar);
	record_prefit_quality(resolutionFileName, chi2ndf, ndataBins, npar, drawFromSavedFit);
	// Binning-robust check: equal-population bins inside the populated range (utils/robust_chi2.h).
	{
		const auto robust = oo_robust_chi2(time_pdf, *data, obs_time, npar);
		oo_robust_chi2_toy(time_pdf, *data, obs_time, npar);
		if (!drawFromSavedFit) oo_write_robust_chi2(resolutionFileName, robust);
	}
	const double chi2 = chi2ndf * ndf;
	const double pvalue = TMath::Prob(chi2, ndf);
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
	oo_diagnostics::run(time_pdf, *data, obs_time, npar, timePlotBins, "ctau_pr_" + figTag, true);


	// ------------------------------------------------------------------
	// print final result
	// ------------------------------------------------------------------
	cout << "------------------ FIT RESULT SUMMARY --------------------" << endl;
	cout << "Prompt resolution components used in PR fit: " << nResolutionComponents << endl;
	cout << Form("Prompt ctau chi2/ndf: %.1f/%d = %.3f, p=%.3g", chi2, ndf, chi2ndf, pvalue) << endl;
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	if (time_result)
		time_result->Print();
	else
		cout << "[PlotOnly] fit result print skipped; using saved prompt ctau parameters." << endl;
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_pr lifetime fit : " << figLifetime << std::endl;
	// time_pdf.Print("V");
}
