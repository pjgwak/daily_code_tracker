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
#include "TParameter.h"
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

static RooParametricStepFunction *makeStepPdfFromHist(
	const char *pdfName,
	RooRealVar &x,
	TH1 *h,
	RooArgList &ownedCoeffs,
	double floorFrac = 1e-6,
	double minCoeff = 1e-8,
	double maxCoeff = 1.0)
{
	if (!h)
		return nullptr;

	const int nBins = h->GetNbinsX();
	if (nBins < 2)
		return nullptr;

	std::vector<double> boundaries(nBins + 1);
	for (int i = 1; i <= nBins; ++i)
		boundaries[i - 1] = h->GetXaxis()->GetBinLowEdge(i);
	boundaries[nBins] = h->GetXaxis()->GetBinUpEdge(nBins);

	TArrayD boundaryArray(static_cast<int>(boundaries.size()), boundaries.data());

	double totalYield = 0.0;
	double minWidth = std::numeric_limits<double>::infinity();
	for (int i = 1; i <= nBins; ++i)
	{
		double c = h->GetBinContent(i);
		if (!std::isfinite(c) || c < 0.0)
			c = 0.0;
		totalYield += c;
		minWidth = std::min(minWidth, h->GetBinWidth(i));
	}
	if (!(totalYield > 0.0) || !(minWidth > 0.0))
		return nullptr;

	ownedCoeffs.removeAll();
	double maxHeight = minCoeff;
	double usedArea = 0.0;
	for (int i = 1; i <= nBins - 1; ++i)
	{
		double c = h->GetBinContent(i);
		if (!std::isfinite(c) || c <= 0.0)
			c = floorFrac * totalYield;

		const double width = h->GetBinWidth(i);
		double frac = c / totalYield;
		double hgt = frac / width;
		hgt = std::max(minCoeff, hgt);
		maxHeight = std::max(maxHeight, hgt);
		usedArea += hgt * width;

		auto *coef = new RooRealVar(
			Form("%s_coef_%02d", pdfName, i - 1),
			Form("%s_coef_%02d", pdfName, i - 1),
			hgt,
			minCoeff,
			std::max(maxCoeff, 5.0 * maxHeight));
		coef->setConstant(true);
		ownedCoeffs.add(*coef);
	}

	const double remainingArea = std::max(0.0, 1.0 - usedArea);
	const double impliedLastHeight = remainingArea / h->GetBinWidth(nBins);
	std::cout << "[DEBUG] " << pdfName
	          << " totalYield = " << totalYield
	          << ", usedArea = " << usedArea
	          << ", remainingArea = " << remainingArea
	          << ", impliedLastHeight = " << impliedLastHeight
	          << std::endl;

	return new RooParametricStepFunction(
		pdfName, pdfName,
		x,
		ownedCoeffs,
		boundaryArray,
		nBins);
}

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

void ctau_pr(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4)
{
	bool isWeight = false;
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	enum TimeErrModelPlaceholder
	{
		kTimeErrDefault = 0,
	};
	enum ErrPdfChoice
	{
		kErrPdfAnalytic = 0,
		kErrPdfHist = 1,
	};
	int nResolutionComponents = 3; // 1~4 components
	int errPdfOpt = kErrPdfHist; // 0: analytic, 1: RooHistPdf
	const int histPdfInterpolationOrder = 1;
	const TimeErrModelPlaceholder timeErrModel = kTimeErrDefault; // placeholder only
	bool useStagedTimeErrFit = false;
	
	if (yLow == 1.6f)
	{
		if (ptLow == 14.0f && ptHigh == 20.0f)
		{
			nResolutionComponents = 2;
			useStagedTimeErrFit = false;
		}
	}
	nResolutionComponents = std::clamp(nResolutionComponents, 1, 4);

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	// First we open the actual RooDataSet used in the prompt ctau analysis.
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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
	const auto ctRange = quantileRange(*dataBasic, *timeTmp, 0.001, 0.999, false);
	auto errRange = quantileRange(*dataBasic, *timeErrTmp, 0.001, 0.995, true);
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
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	RooRealVar &obs_time = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3D"));
	RooRealVar &obs_timeErr = *static_cast<RooRealVar *>(dataSel->get()->find("ctau3DErr"));
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	obs_time.SetTitle("#font[12]{l}_{J/#psi}");
	obs_time.setUnit("mm");
	obs_timeErr.SetTitle("event-by-event ctau error");
	obs_timeErr.setUnit("mm");
	obs_time.setRange(ctRange.first, ctRange.second);
	obs_timeErr.setRange(errRange.first, errRange.second);
	obs_timeErr.setMin(errRange.first);
	const int massPlotBins = std::max(2, obs_mass.getBins() * 2);
	const int timePlotBins = std::max(2, obs_time.getBins());
	const int timeErrPlotBins = std::max(2, obs_timeErr.getBins());

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
	RooRealVar timeErrGaus1Mean("timeErrGaus1Mean", "timeErrGaus1Mean", gaus1MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus1Sigma("timeErrGaus1Sigma", "timeErrGaus1Sigma", gaus1SigmaInit, 5e-4, std::max(0.15 * errSpan, 2e-3));
	RooRealVar timeErrGaus2Mean("timeErrGaus2Mean", "timeErrGaus2Mean", gaus2MeanInit, errRange.first, errRange.second);
	RooRealVar timeErrGaus2Sigma("timeErrGaus2Sigma", "timeErrGaus2Sigma", gaus2SigmaInit, 1e-3, std::max(0.30 * errSpan, 4e-3));
	RooRealVar timeErrTailMpv("timeErrTailMpv", "timeErrTailMpv", tailMpvInit, errRange.first, errRange.second);
	RooRealVar timeErrTailWidth("timeErrTailWidth", "timeErrTailWidth", tailWidthInit, 1e-3, std::max(0.50 * errSpan, 5e-3));
	RooRealVar timeErrTail2Mpv("timeErrTail2Mpv", "timeErrTail2Mpv", tail2MpvInit, errRange.first, errRange.second);
	RooRealVar timeErrTail2Width("timeErrTail2Width", "timeErrTail2Width", tail2WidthInit, 1e-3, std::max(0.80 * errSpan, 8e-3));
	RooRealVar timeErrLognM0("timeErrLognM0", "timeErrLognM0", lognM0Init, std::max(errRange.first, 1e-6), errRange.second);
	RooRealVar timeErrLognK("timeErrLognK", "timeErrLognK", lognKInit, 0.05, 2.0);
	RooRealVar timeErrCore1Frac("timeErrCore1Frac", "timeErrCore1Frac", 0.60, 0.05, 0.95);
	RooRealVar timeErrTailFrac("timeErrTailFrac", "timeErrTailFrac", 0.08, 0.001, 0.30);
	RooRealVar timeErrTail2Frac("timeErrTail2Frac", "timeErrTail2Frac", 0.05, 0.001, 0.30);
	RooRealVar timeErrLognFrac("timeErrLognFrac", "timeErrLognFrac", 0.03, 0.0005, 0.20);
	if (useStagedTimeErrFit)
	{
		// In the forward 14-20 bin, seed the Landau terms in the genuine high-error tail
		// so they do not collapse into narrow core-like components.
		timeErrTailMpv.setVal(errRange.first + 0.18 * errSpan);
		timeErrTailMpv.setRange(errRange.first + 0.10 * errSpan, errRange.first + 0.55 * errSpan);
		timeErrTailWidth.setVal(std::max(0.08 * errSpan, 1.5e-3));
		timeErrTailWidth.setRange(std::max(0.03 * errSpan, 1.2e-3), std::max(0.35 * errSpan, 4e-3));
		timeErrTail2Mpv.setVal(errRange.first + 0.35 * errSpan);
		timeErrTail2Mpv.setRange(errRange.first + 0.18 * errSpan, errRange.second);
		timeErrTail2Width.setVal(std::max(0.12 * errSpan, 2.0e-3));
		timeErrTail2Width.setRange(std::max(0.05 * errSpan, 1.5e-3), std::max(0.60 * errSpan, 6e-3));
		timeErrLognM0.setVal(errRange.first + 0.78 * errSpan);
		timeErrLognM0.setRange(errRange.first + 0.45 * errSpan, errRange.second);
		timeErrLognK.setVal(0.22);
		timeErrLognK.setRange(0.05, 0.80);
		timeErrTailFrac.setVal(0.03);
		timeErrTailFrac.setRange(0.0005, 0.15);
		timeErrTail2Frac.setVal(0.01);
		timeErrTail2Frac.setRange(0.0001, 0.08);
		timeErrLognFrac.setVal(0.01);
		timeErrLognFrac.setRange(0.0001, 0.08);
	}
	RooGaussian timeErrGaus1("timeErrGaus1", "timeErrGaus1", obs_timeErr, timeErrGaus1Mean, timeErrGaus1Sigma);
	RooGaussian timeErrGaus2("timeErrGaus2", "timeErrGaus2", obs_timeErr, timeErrGaus2Mean, timeErrGaus2Sigma);
	RooLandau timeErrTail("timeErrTail", "timeErrTail", obs_timeErr, timeErrTailMpv, timeErrTailWidth);
	RooLandau timeErrTail2("timeErrTail2", "timeErrTail2", obs_timeErr, timeErrTail2Mpv, timeErrTail2Width);
	RooLognormal timeErrLogn("timeErrLogn", "timeErrLogn", obs_timeErr, timeErrLognM0, timeErrLognK);
	std::unique_ptr<RooAddPdf> timeErrPdfAnalytic = std::make_unique<RooAddPdf>(
			"timeErrPdf", "timeErrPdf",
			RooArgList(timeErrTail, timeErrTail2, timeErrLogn, timeErrGaus1, timeErrGaus2),
			RooArgList(timeErrTailFrac, timeErrTail2Frac, timeErrLognFrac, timeErrCore1Frac),
			true);
	std::unique_ptr<RooDataHist> timeErrHistData;
	std::unique_ptr<RooHistPdf> timeErrPdfHist;
	RooAbsPdf *timeErrPdf = timeErrPdfAnalytic.get();
	// Fit the standalone ctau-error distribution first, then freeze it for the lifetime fit.
	if (errPdfOpt == kErrPdfAnalytic && useStagedTimeErrFit)
	{
		// First determine the core with a 2-Gaussian fit, then release the Landau tails.
		RooAddPdf timeErrPdfCoreOnly(
			"timeErrPdfCoreOnly", "timeErrPdfCoreOnly",
			RooArgList(timeErrGaus1, timeErrGaus2),
			RooArgList(timeErrCore1Frac),
			true);
		timeErrPdfCoreOnly.fitTo(*timeErrData, PrintLevel(-1), SumW2Error(isWeight));
	}
	RooFitResult *timeErrResult = nullptr;
	if (errPdfOpt == kErrPdfHist)
	{
		timeErrHistData = std::make_unique<RooDataHist>(
			"timeErrHistData", "timeErrHistData",
			RooArgSet(obs_timeErr), *timeErrData);
		timeErrPdfHist = std::make_unique<RooHistPdf>(
			"timeErrHistPdf", "timeErrHistPdf",
			RooArgSet(obs_timeErr), *timeErrHistData, histPdfInterpolationOrder);
		timeErrPdf = timeErrPdfHist.get();
	}
	else
	{
		timeErrResult = timeErrPdfAnalytic->fitTo(*timeErrData, Save(true), PrintLevel(-1), SumW2Error(isWeight));
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
		timeErrCore1Frac.setConstant(true);
		timeErrTailFrac.setConstant(true);
		timeErrTail2Frac.setConstant(true);
		timeErrLognFrac.setConstant(true);
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
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), DataError(RooAbsData::SumW2), Name("data"));
	else
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), Name("data"));
	timeErrPdf->plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (errPdfOpt == kErrPdfAnalytic)
	{
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrTail2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("tail2"));
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrLogn), LineColor(kMagenta + 1), LineStyle(kDashed), LineWidth(2), Name("logn"));
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
		timeErrPdf->plotOn(timeErrPlot, Components(timeErrGaus2), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("gaus2"));
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
		timeErrLeg.AddEntry(o, errPdfOpt == kErrPdfHist ? "RooHistPdf" : "Fit", "l");
	if (errPdfOpt == kErrPdfAnalytic)
	{
		if (auto *o = findObj(timeErrPlot, "tail"))
			timeErrLeg.AddEntry(o, "Landau tail", "l");
		if (auto *o = findObj(timeErrPlot, "tail2"))
			timeErrLeg.AddEntry(o, "Landau tail 2", "l");
		if (auto *o = findObj(timeErrPlot, "logn"))
			timeErrLeg.AddEntry(o, "Log-normal tail", "l");
		if (auto *o = findObj(timeErrPlot, "gaus1"))
			timeErrLeg.AddEntry(o, "Gauss 1", "l");
		if (auto *o = findObj(timeErrPlot, "gaus2"))
			timeErrLeg.AddEntry(o, "Gauss 2", "l");
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
			int minimize = -1;
			int hesse = -1;
			if (timeErrResult)
			{
				for (UInt_t i = 0, n = timeErrResult->numStatusHistory(); i < n; ++i)
				{
					const char *lab = timeErrResult->statusLabelHistory(i);
					if (lab)
					{
						if (!strcmp(lab, "MINIMIZE") || !strcmp(lab, "MIGRAD"))
							minimize = timeErrResult->statusCodeHistory(i);
						else if (!strcmp(lab, "HESSE"))
							hesse = timeErrResult->statusCodeHistory(i);
					}
			}
		}
		if (errPdfOpt == kErrPdfHist)
			tx.DrawLatex(0.19, 0.765, Form("Status : RooHistPdf template (%d)", histPdfInterpolationOrder));
		else
			tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
	}
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
			tp.DrawLatex(xtext, y0 + dy * k++, Form("bins = %d", obs_timeErr.getBins()));
		}
		else
		{
					printVar("mean_{1}", timeErrGaus1Mean);
					printVar("#sigma_{1}", timeErrGaus1Sigma);
					printVar("mean_{2}", timeErrGaus2Mean);
					printVar("#sigma_{2}", timeErrGaus2Sigma);
					printVar("mpv_{L}", timeErrTailMpv);
					printVar("#sigma_{L}", timeErrTailWidth);
					printVar("mpv_{L2}", timeErrTail2Mpv);
					printVar("#sigma_{L2}", timeErrTail2Width);
					printVar("m0_{LN}", timeErrLognM0);
					printVar("k_{LN}", timeErrLognK);
					printVar("f_{tail}", timeErrTailFrac);
					printVar("f_{tail2}", timeErrTail2Frac);
					printVar("f_{LN}", timeErrLognFrac);
					printVar("f_{core1}", timeErrCore1Frac);
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", timeErrChi.first, timeErrNdf, timeErrPvalue));
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
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale", 1.0, 0.5, 3.0);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", 0.0, -0.5, 0.5);
	ctauTime2Mean.setConstant(true);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", 0.8, 0.05, 5.0);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", 0.0, -0.5, 0.5);
	ctauTime3Mean.setConstant(true);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", 1.7, 0.05, 10.0);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", 0.0, -0.5, 0.5);
	ctauTime4Mean.setConstant(true);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", 2.5, 0.05, 15.0);
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
		Nctau1.setVal(data->numEntries() * 0.20);
		Nctau2.setVal(data->numEntries() * 0.60);
		Nctau3.setVal(data->numEntries() * 0.20);
		Nctau4.setVal(data->numEntries() * 0.05);

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
	RooFitResult *time_result = time_pdf.fitTo(*data, Extended(), Save(), SumW2Error(isWeight));
	time_result->Print();

		const double nCtauTotal = std::max(1e-12,
			Nctau1.getVal() +
			(nResolutionComponents >= 2 ? Nctau2.getVal() : 0.0) +
			(nResolutionComponents >= 3 ? Nctau3.getVal() : 0.0) +
			(nResolutionComponents >= 4 ? Nctau4.getVal() : 0.0));
		const double ctauFrac1 = Nctau1.getVal() / nCtauTotal;
		const double ctauFrac2 = Nctau2.getVal() / nCtauTotal;
		const double ctauFrac3 = Nctau3.getVal() / nCtauTotal;
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", resultDir.Data(), figTag.Data());
	{
		TFile resolutionFile(resolutionFileName, "RECREATE");
		if (!resolutionFile.IsZombie())
		{
			TParameter<int>("errPdfOpt", errPdfOpt).Write();
			TParameter<int>("histPdfInterpolationOrder", histPdfInterpolationOrder).Write();
			TParameter<int>("timeErrPlotBins", timeErrPlotBins).Write();
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
				printVar("mean_{1}", ctauTime1Mean);
				printVar("#sigma_{1}", ctauTime1Scale);
				if (nResolutionComponents >= 2)
				{
					printVar("mean_{2}", ctauTime2Mean);
					printVar("#Delta #sigma_{21}", ctauTime2Delta);
				}
				if (nResolutionComponents >= 3)
				{
					printVar("N_{3}", Nctau3);
					printVar("mean_{3}", ctauTime3Mean);
					printVar("#Delta #sigma_{32}", ctauTime3Delta);
				}
				if (nResolutionComponents >= 4)
				{
					printVar("N_{4}", Nctau4);
					printVar("mean_{4}", ctauTime4Mean);
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", chiM.first, ndf, pvalue));
	}

	TLine line(ctRange.first, 0.0, ctRange.second, 0.0);
	line.SetLineStyle(2);
	line.Draw("same");

	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

	// ------------------------------------------------------------------
	// print final result
	// ------------------------------------------------------------------
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	cout << "Prompt resolution components used in PR fit: " << nResolutionComponents << endl;
	cout << "Time-error model placeholder in PR fit: " << static_cast<int>(timeErrModel) << endl;
	cout << "errPdfOpt in PR fit: " << errPdfOpt
			 << (errPdfOpt == kErrPdfHist ? " (RooHistPdf)" : " (analytic)") << endl;
	if (timeErrResult)
	{
		cout << "------------------ FIT RESULT FOR TIME ERR ---------------" << endl;
		timeErrResult->Print("v");
	}
	time_result->Print();
	const TString figTimeErr = figName("timeerr_model");
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_pr err fit : " << figTimeErr << std::endl;
	std::cout << "[FIG] ctau_pr lifetime fit : " << figLifetime << std::endl;
	// time_pdf.Print("V");
}

// //
// 	// THE LIFETIME FIT: decompose into two resolved components for plotting
// 	RooProdPdf time_pdf_component1(
// 		"time_pdf_component1",
// 		"time_pdf_component1",
// 		RooArgSet(timeErrPdf),
// 		Conditional(RooArgSet(ctauTime1), RooArgSet(obs_time))
// 	);
// 	RooProdPdf time_pdf_component2(
// 		"time_pdf_component2",
// 		"time_pdf_component2",
// 		RooArgSet(timeErrPdf),
// 		Conditional(RooArgSet(ctauTime2), RooArgSet(obs_time))
// 	);
// 	RooAddPdf time_pdf(
// 		"time_pdf",
// 		"time_pdf",
// 		RooArgList(time_pdf_component1, time_pdf_component2),
// 		RooArgList(Nctau1, Nctau2)
// 	);
// 	RooFitResult *time_result = time_pdf.fitTo(*data,Extended(),Save());
// 	time_result->Print();
