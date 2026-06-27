#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooResolutionModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TParameter.h"
#include "TString.h"
#include "TAxis.h"
#include "RooHist.h"
#include "saved_fit_helpers.h"
#include <algorithm>
#include <cmath>
#include <cstring>
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

static std::pair<double, double> quantileRange(RooDataSet &ds, RooRealVar &var, double qLo, double qHi, bool positiveOnly)
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
}

void ctau_np(
	float ptLow = 6.5, float ptHigh = 9,
	float yLow = 0, float yHigh = 1.6,
	float muPtCut = 0.0,
	int PR = 1, // kept for compatibility: 0=PR, 1=NP, 2=Inc.
	float ctauCut = 0.08,
	int nSignalSSComponents = 2, // 1~3
	bool drawFromSavedFit = false,
	bool publish = false,
	int ssOverride = -1,
	bool retryFailedFit = true)
{
	ScopedMacroTimer timer("ctau_np", ptLow, ptHigh, yLow, yHigh);
	if (publish)
		drawFromSavedFit = true;
	bool isWeight = false;
	(void)muPtCut;
	(void)PR;
	(void)ctauCut;

	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	if (yLow == 1.6f && ptLow >= 3.0f && ptHigh <= 6.5f)
		nSignalSSComponents = 3;
	nSignalSSComponents = std::clamp(nSignalSSComponents, 1, 3);
	if (ssOverride >= 1)
		nSignalSSComponents = std::clamp(ssOverride, 1, 3);

	// ------------------------------------------------------------------
	// load nonprompt reco-MC dataset
	// ------------------------------------------------------------------
	const char *DATA_ROOT = "/data/users/pjgwak/work/raa_pb18/run2_raa_pbpb2018/skimmedFiles/OniaRooDataSet_isMC1_BtoJpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_260303_root634.root";
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

	TString accCut = "(((abs(eta1) <= 1.2) && (pt1 >= 3.5)) || ((abs(eta2) <= 1.2) && (pt2 >= 3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)))";
	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && (recoQQsign==0) && %s",
		ptLow, ptHigh, yLow, yHigh, accCut.Data());
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	auto *massTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
	auto *timeTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
	if (!massTmp || !timeTmp)
	{
		std::cerr << "ERROR: required variables mass/ctau3D are missing in dataset." << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// determine fit ranges
	// ------------------------------------------------------------------
	auto ctRangeRaw = quantileRange(*dataBasic, *timeTmp, 0.01, 0.999, false);
	ctRangeRaw.first = 0.05; // hard cut, NP decay have > 0
	const auto ctRange = std::make_pair(std::max(0.0, ctRangeRaw.first), ctRangeRaw.second);
	TString cutAll = Form(
		"%s && !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes) && "
		"(ctau3D >= %g && ctau3D <= %g)",
		cutBasic.Data(), ctRange.first, ctRange.second);
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
	const TString figBaseDir = publish ? "figs_publish" : "figs";
	const TString figDir = TString::Format("%s/%s/ctau_np", figBaseDir.Data(), yTag.Data());
	const TString prResultDir = TString::Format("roots/%s/ctau_pr", yTag.Data());
	const TString resultDir = TString::Format("roots/%s/ctau_np", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", prResultDir.Data(), figTag.Data());
	const TString modelFileName = TString::Format("%s/ctau_np_model_%s.root", resultDir.Data(), figTag.Data());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit)
	{
		if (!load_saved_fit_file(savedFitFile, modelFileName, "nonprompt ctau"))
			return;
		nSignalSSComponents = std::clamp(read_saved_int_param(savedFitFile.get(), "nSignalSSComponents", nSignalSSComponents), 1, 3);
	}

	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_time.SetTitle("#font[12]{l}_{J/#psi}");
	obs_time.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	const int timePlotBins = 500;
	obs_time.setBins(timePlotBins);
	RooDataSet *data = dataSel.get();
	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	// ------------------------------------------------------------------
	// load prompt resolution model from ctau_pr
	// ------------------------------------------------------------------
	TFile *resolutionFile = TFile::Open(resolutionFileName);
	if (!resolutionFile || resolutionFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open prompt resolution file: " << resolutionFileName << std::endl;
		std::cerr << "       Run ctau_pr for this (pt, y) bin before ctau_np." << std::endl;
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

	const double compFrac1Val = (nResolutionComponents == 1) ? 1.0 : ctauFrac1Val;
	const double compFrac2Val = (nResolutionComponents >= 2)
		? ((nResolutionComponents == 2) ? (1.0 - ctauFrac1Val) : (1.0 - ctauFrac1Val) * ctauFrac2Val)
		: 0.0;
	const double compFrac3Val = (nResolutionComponents >= 3)
		? ((nResolutionComponents == 3) ? (1.0 - ctauFrac1Val) * (1.0 - ctauFrac2Val)
										: (1.0 - ctauFrac1Val) * (1.0 - ctauFrac2Val) * ctauFrac3Val)
		: 0.0;

	// ------------------------------------------------------------------
	// build signal ctau model: reference-style resolution convolution
	// ------------------------------------------------------------------
	const double signalLifetimeFloor = 5e-3;
	const double bkgLifetimeFloor = 1e-2;
	const double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

	const double signalLifetimeInit = std::max(0.08, 1.5 * signalLifetimeFloor);
	const double signalLifetimeCeil = std::max(signalLifetimeInit * 5.0, maxStableLifetime);
	RooRealVar signal_log_lifetime_offset(
		"signal_log_lifetime_offset", "log(signal_lifetime - floor)",
		std::log(std::max(signalLifetimeInit - signalLifetimeFloor, 1e-4)),
		std::log(1e-2),
		std::log(std::max(signalLifetimeCeil - signalLifetimeFloor, 2e-4)));
	RooFormulaVar signal_lifetime(
		"signal_lifetime",
		Form("%g + exp(@0)", signalLifetimeFloor),
		RooArgList(signal_log_lifetime_offset));

	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", ctauMeanScaleVal);
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", ctauTime1MeanVal);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale", ctauTime1ScaleVal);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", ctauTime2MeanVal);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", std::max(1e-6, ctauTime2ScaleVal - ctauTime1ScaleVal));
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale);
	RooRealVar ctauFrac1("ctauFrac1", "ctauFrac1", compFrac1Val);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", ctauTime3MeanVal);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", std::max(1e-6, ctauTime3ScaleVal - ctauTime2ScaleVal));
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale);
	RooRealVar ctauFrac2("ctauFrac2", "ctauFrac2", compFrac2Val);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", ctauTime4MeanVal);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", std::max(1e-6, ctauTime4ScaleVal - ctauTime3ScaleVal));
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale);
	RooRealVar ctauFrac3("ctauFrac3", "ctauFrac3", compFrac3Val);
	ctauTime1Mean.setConstant(true);
	ctauTime1Scale.setConstant(true);
	ctauTime2Mean.setConstant(true);
	ctauTime2Delta.setConstant(true);
	ctauFrac1.setConstant(true);
	ctauTime3Mean.setConstant(true);
	ctauTime3Delta.setConstant(true);
	ctauFrac2.setConstant(true);
	ctauTime4Mean.setConstant(true);
	ctauTime4Delta.setConstant(true);
	ctauFrac3.setConstant(true);

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

	RooDecay signal_ss1_time("signal_ss1_time", "signal_ss1_time", obs_time, signal_lifetime, time_resolution, RooDecay::SingleSided);
	const double signalLifetime2Floor = std::max(0.75 * signalLifetimeFloor, 5e-3);
	const double signalLifetime2Init = std::max(0.20, 2.0 * signalLifetime2Floor);
	const double signalLifetime2Ceil = std::max(signalLifetime2Init * 5.0, maxStableLifetime);
	RooRealVar signal2_log_lifetime_offset(
		"signal2_log_lifetime_offset", "log(signal2_lifetime - floor)",
		std::log(std::max(signalLifetime2Init - signalLifetime2Floor, 1e-4)),
		std::log(1e-2),
		std::log(std::max(signalLifetime2Ceil - signalLifetime2Floor, 2e-4)));
	RooFormulaVar signal2_lifetime(
		"signal2_lifetime",
		Form("%g + exp(@0)", signalLifetime2Floor),
		RooArgList(signal2_log_lifetime_offset));
	RooDecay signal_ss2_time("signal_ss2_time", "signal_ss2_time", obs_time, signal2_lifetime, time_resolution, RooDecay::SingleSided);
	const double signalLifetime3Floor = std::max(1.00 * signalLifetimeFloor, 5e-3);
	const double signalLifetime3Init = std::max(0.40, 2.0 * signalLifetime3Floor);
	const double signalLifetime3Ceil = std::max(signalLifetime3Init * 5.0, maxStableLifetime);
	RooRealVar signal3_log_lifetime_offset(
		"signal3_log_lifetime_offset", "log(signal3_lifetime - floor)",
		std::log(std::max(signalLifetime3Init - signalLifetime3Floor, 1e-4)),
		std::log(1e-2),
		std::log(std::max(signalLifetime3Ceil - signalLifetime3Floor, 2e-4)));
	RooFormulaVar signal3_lifetime(
		"signal3_lifetime",
		Form("%g + exp(@0)", signalLifetime3Floor),
		RooArgList(signal3_log_lifetime_offset));
	RooDecay signal_ss3_time("signal_ss3_time", "signal_ss3_time", obs_time, signal3_lifetime, time_resolution, RooDecay::SingleSided);
	RooRealVar NsignalSS1("NsignalSS1", "NsignalSS1", 0.50 * data->numEntries(), 0.0, 5.0 * std::max(1, data->numEntries()));
	RooRealVar NsignalSS2("NsignalSS2", "NsignalSS2", 0.30 * data->numEntries(), 0.0, 5.0 * std::max(1, data->numEntries()));
	RooRealVar NsignalSS3("NsignalSS3", "NsignalSS3", 0.20 * data->numEntries(), 0.0, 5.0 * std::max(1, data->numEntries()));
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
	// fit or restore the signal ctau model
	// ------------------------------------------------------------------
	if (drawFromSavedFit)
	{
		signal_log_lifetime_offset.setVal(read_saved_double_param(savedFitFile.get(), "signal_log_lifetime_offset", signal_log_lifetime_offset.getVal()));
		NsignalSS1.setVal(read_saved_double_param(savedFitFile.get(), "NsignalSS1", NsignalSS1.getVal()));
		if (nSignalSSComponents >= 2)
		{
			signal2_log_lifetime_offset.setVal(read_saved_double_param(savedFitFile.get(), "signal2_log_lifetime_offset", signal2_log_lifetime_offset.getVal()));
			NsignalSS2.setVal(read_saved_double_param(savedFitFile.get(), "NsignalSS2", NsignalSS2.getVal()));
		}
		if (nSignalSSComponents >= 3)
		{
			signal3_log_lifetime_offset.setVal(read_saved_double_param(savedFitFile.get(), "signal3_log_lifetime_offset", signal3_log_lifetime_offset.getVal()));
			NsignalSS3.setVal(read_saved_double_param(savedFitFile.get(), "NsignalSS3", NsignalSS3.getVal()));
		}
	}
	std::unique_ptr<RooFitResult> timeResultHolder;
	RooFitResult *time_result = nullptr;
	if (!drawFromSavedFit)
	{
		timeResultHolder.reset(signal_time->fitTo(*data, Extended(), Save(), SumW2Error(isWeight), Strategy(2), PrintLevel(0), Minimizer("Minuit2", "migrad")));
		time_result = timeResultHolder.get();
		if (retryFailedFit && (!time_result || time_result->status() != 0 || time_result->covQual() < 2))
		{
			const int status = time_result ? time_result->status() : -999;
			const int covQual = time_result ? time_result->covQual() : -999;
			std::cout << "[Retry] ctau_np fit status=" << status << " covQual=" << covQual << "; retry default minimizer" << std::endl;
			timeResultHolder.reset(signal_time->fitTo(*data, Extended(), Save(), SumW2Error(isWeight), Strategy(2), PrintLevel(0)));
			time_result = timeResultHolder.get();
		}
		if (time_result)
			time_result->Print();
	}
	else
		std::cout << "[PlotOnly] Loaded saved nonprompt ctau parameters and skipped fit: " << modelFileName << std::endl;

	if (!drawFromSavedFit)
	{
		TFile modelFile(modelFileName, "RECREATE");
		if (!modelFile.IsZombie())
		{
			TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
			TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
			TParameter<double>("ctauMeanScale", ctauMeanScale.getVal()).Write();
			TParameter<double>("ctauTime1Mean", ctauTime1Mean.getVal()).Write();
			TParameter<double>("ctauTime1Scale", ctauTime1Scale.getVal()).Write();
			if (nResolutionComponents >= 2)
			{
				TParameter<double>("ctauTime2Mean", ctauTime2Mean.getVal()).Write();
				TParameter<double>("ctauTime2Scale", ctauTime2Scale.getVal()).Write();
				TParameter<double>("ctauFrac1", ctauFrac1.getVal()).Write();
			}
			if (nResolutionComponents >= 3)
			{
				TParameter<double>("ctauTime3Mean", ctauTime3Mean.getVal()).Write();
				TParameter<double>("ctauTime3Scale", ctauTime3Scale.getVal()).Write();
				TParameter<double>("ctauFrac2", ctauFrac2.getVal()).Write();
			}
			if (nResolutionComponents >= 4)
			{
				TParameter<double>("ctauTime4Mean", ctauTime4Mean.getVal()).Write();
				TParameter<double>("ctauTime4Scale", ctauTime4Scale.getVal()).Write();
				TParameter<double>("ctauFrac3", ctauFrac3.getVal()).Write();
			}
			TParameter<double>("signal_lifetime", signal_lifetime.getVal()).Write();
			TParameter<double>("signalLifetimeFloor", signalLifetimeFloor).Write();
			TParameter<double>("signal_log_lifetime_offset", signal_log_lifetime_offset.getVal()).Write();
			TParameter<double>("NsignalSS1", NsignalSS1.getVal()).Write();
			if (nSignalSSComponents >= 2)
			{
				TParameter<double>("signal2_lifetime", signal2_lifetime.getVal()).Write();
				TParameter<double>("signalLifetime2Floor", signalLifetime2Floor).Write();
				TParameter<double>("signal2_log_lifetime_offset", signal2_log_lifetime_offset.getVal()).Write();
				TParameter<double>("NsignalSS2", NsignalSS2.getVal()).Write();
			}
			if (nSignalSSComponents >= 3)
			{
				TParameter<double>("signal3_lifetime", signal3_lifetime.getVal()).Write();
				TParameter<double>("signalLifetime3Floor", signalLifetime3Floor).Write();
				TParameter<double>("signal3_log_lifetime_offset", signal3_log_lifetime_offset.getVal()).Write();
				TParameter<double>("NsignalSS3", NsignalSS3.getVal()).Write();
			}
			if (time_result)
				time_result->Write("fitResult");
			modelFile.Write();
		}
		else
		{
			std::cerr << "ERROR: cannot create NP model file: " << modelFileName << std::endl;
		}
		std::cout << "Saved ctau NP model parameters to " << modelFileName << std::endl;
	}
	else
		std::cout << "[PlotOnly] Left saved ctau NP ROOT file unchanged: " << modelFileName << std::endl;

	RooAbsPdf *time_pdf_ss1 = &signal_ss1_time;
	RooAbsPdf *time_pdf_ss2 = nSignalSSComponents >= 2 ? static_cast<RooAbsPdf *>(&signal_ss2_time) : nullptr;
	RooAbsPdf *time_pdf_ss3 = nSignalSSComponents == 3 ? static_cast<RooAbsPdf *>(&signal_ss3_time) : nullptr;

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
	signal_time->plotOn(timefitplot, LineColor(kBlack), LineWidth(2), Name("model"));
	time_pdf_ss1->plotOn(timefitplot, LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Normalization(NsignalSS1.getVal(), RooAbsReal::NumEvent), Name("ss1_component"));
	if (time_pdf_ss2)
		time_pdf_ss2->plotOn(timefitplot, LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Normalization(NsignalSS2.getVal(), RooAbsReal::NumEvent), Name("ss2_component"));
	if (time_pdf_ss3)
		time_pdf_ss3->plotOn(timefitplot, LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Normalization(NsignalSS3.getVal(), RooAbsReal::NumEvent), Name("ss3_component"));
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
		tx.DrawLatex(0.96, 0.935, "pp #sqrt{s} = 5.02 TeV (28.0 pb^{-1})");
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
		auto printFitVar = [&](const char *title, const RooAbsReal &var)
		{
			const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
			if (rrv)
			{
				if (rrv->isConstant())
					return;
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
				return;
			}
			if (time_result)
			{
				const double err = var.getPropagatedError(*time_result);
				if (err > 0.0 && std::isfinite(err))
				{
					tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
					return;
				}
			}
			tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
		};
		printFitVar("N_{SS1}", NsignalSS1);
		if (nSignalSSComponents >= 2)
			printFitVar("N_{SS2}", NsignalSS2);
		if (nSignalSSComponents == 3)
			printFitVar("N_{SS3}", NsignalSS3);
		printFitVar("#tau^{raw}_{SS1}", signal_log_lifetime_offset);
		if (nSignalSSComponents >= 2)
			printFitVar("#tau^{raw}_{SS2}", signal2_log_lifetime_offset);
		if (nSignalSSComponents == 3)
			printFitVar("#tau^{raw}_{SS3}", signal3_log_lifetime_offset);
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
	std::cout << Form("Nonprompt ctau chi2/ndf: %.1f/%d = %.3f, p=%.3g",
		chiM.first, ndf, ndf > 0 ? chiM.first / ndf : 0.0, pvalue) << std::endl;
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
	std::cout << "------------------ FIT RESULT SUMMARY --------------------" << std::endl;
	std::cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << std::endl;
	if (time_result)
		time_result->Print("v");
	else
		std::cout << "[PlotOnly] fit result print skipped; using saved nonprompt ctau parameters." << std::endl;
	std::cout << "Prompt resolution components used in NP fit: " << nResolutionComponents << std::endl;
	std::cout << "SingleSided components used in NP fit: " << nSignalSSComponents << std::endl;
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_np lifetime fit : " << figLifetime << std::endl;
}
