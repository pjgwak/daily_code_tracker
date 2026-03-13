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

void ctau_np(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4)
{
	bool isWeight = false;
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	// Choose how many SingleSided signal components to use in each (pt, y) bin.
	int nSignalSSComponents = 2;
	if (yLow == 0.0f)
	{
		if (ptLow == 200.0f && ptHigh == 350.0f)
		{
			nSignalSSComponents = 2; // 1~3
		}
	}
	nSignalSSComponents = std::clamp(nSignalSSComponents, 1, 3);

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	// First we open the actual RooDataSet used in the prompt ctau analysis.
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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
	const int timePlotBins = std::max(2, obs_time.getBins() * 2);
	const int timeErrPlotBins = std::max(2, obs_timeErr.getBins() * 2);
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", prResultDir.Data(), figTag.Data());

	// Keep the decay constants away from the numerically unstable tau -> 0 limit.
	const double signalLifetimeFloor = std::max(0.5 * errRange.first, 5e-3);
	const double bkgLifetimeFloor = std::max(1.0 * errRange.first, 1e-2);
	const double maxStableLifetime = std::max(ctRange.second - ctRange.first, 20.0 * bkgLifetimeFloor);

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
	RooAddPdf timeErrPdf(
			"timeErrPdf", "timeErrPdf",
			RooArgList(timeErrTail, timeErrGaus1, timeErrGaus2),
			RooArgList(timeErrTailFrac, timeErrCore1Frac),
			true);
	std::unique_ptr<RooFitResult> timeErrResult(timeErrPdf.fitTo(*timeErrData, Save(), PrintLevel(-1), SumW2Error(isWeight)));
	timeErrGaus1Mean.setConstant(true);
	timeErrGaus1Sigma.setConstant(true);
	timeErrGaus2Mean.setConstant(true);
	timeErrGaus2Sigma.setConstant(true);
	timeErrTailMpv.setConstant(true);
	timeErrTailWidth.setConstant(true);
	timeErrCore1Frac.setConstant(true);
	timeErrTailFrac.setConstant(true);

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
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), DataError(RooAbsData::SumW2), Name("data"));
	else
		timeErrData->plotOn(timeErrPlot, Binning(timeErrPlotBins), Name("data"));
	timeErrPdf.plotOn(timeErrPlot, LineColor(kBlack), LineWidth(2), Name("model"));
	timeErrPdf.plotOn(timeErrPlot, Components(timeErrTail), LineColor(kBlue + 2), LineStyle(kDashed), LineWidth(2), Name("tail"));
	timeErrPdf.plotOn(timeErrPlot, Components(timeErrGaus1), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("gaus1"));
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
	if (auto *o = findObj(timeErrPlot, "tail"))
		timeErrLeg.AddEntry(o, "Landau tail", "l");
	if (auto *o = findObj(timeErrPlot, "gaus1"))
		timeErrLeg.AddEntry(o, "Gauss 1", "l");
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
		printVar("#mu_{1}", timeErrGaus1Mean);
		printVar("#sigma_{1}", timeErrGaus1Sigma);
		printVar("#mu_{2}", timeErrGaus2Mean);
		printVar("#sigma_{2}", timeErrGaus2Sigma);
		printVar("mpv_{L}", timeErrTailMpv);
		printVar("#sigma_{L}", timeErrTailWidth);
		printVar("f_{tail}", timeErrTailFrac);
		printVar("f_{G1}", timeErrCore1Frac);
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
	const int timeErrNPar = timeErrResult ? timeErrResult->floatParsFinal().getSize() : 0;
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
	// build signal ctau model
	// ------------------------------------------------------------------
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
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", ctauTime2MeanVal);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", std::max(0.05, ctauTime2ScaleVal - ctauTime1ScaleVal));
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac1("ctauFrac1", "ctauFrac1", ctauFrac1Val);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", ctauTime3MeanVal);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", std::max(0.05, ctauTime3ScaleVal - ctauTime2ScaleVal));
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac2("ctauFrac2", "ctauFrac2", ctauFrac2Val);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", ctauTime4MeanVal);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal));
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale, ctauMeanScale, obs_timeErr);
	RooRealVar ctauFrac3("ctauFrac3", "ctauFrac3", ctauFrac3Val);
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

	// Build up to three positive lifetime components that share the prompt resolution kernel.
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
	// fit the signal ctau model
	// ------------------------------------------------------------------
	RooProdPdf time_pdf("time_pdf", "time_pdf", RooArgSet(timeErrPdf), Conditional(RooArgSet(*signal_time), RooArgSet(obs_time)));
	RooFitResult *time_result = time_pdf.fitTo(*data, Extended(), Save(), SumW2Error(isWeight));
	time_result->Print();

	const TString npModelFileName = TString::Format("%s/ctau_np_model_%s.root", resultDir.Data(), figTag.Data());
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
			modelFile.Write();
		}
		else
		{
			std::cerr << "ERROR: cannot create NP model file: " << npModelFileName << std::endl;
		}
	}
	std::cout << "Saved ctau NP model parameters to " << npModelFileName << std::endl;

	// Build component-only PDFs for plotting with the same conditional structure.
	std::unique_ptr<RooProdPdf> time_pdf_ss1;
	std::unique_ptr<RooProdPdf> time_pdf_ss2;
	std::unique_ptr<RooProdPdf> time_pdf_ss3;
	time_pdf_ss1 = std::make_unique<RooProdPdf>(
			"time_pdf_ss1",
			"time_pdf_ss1",
			RooArgSet(timeErrPdf),
			Conditional(RooArgSet(signal_ss1_time), RooArgSet(obs_time)));
	if (nSignalSSComponents >= 2)
	{
		time_pdf_ss2 = std::make_unique<RooProdPdf>(
				"time_pdf_ss2",
				"time_pdf_ss2",
				RooArgSet(timeErrPdf),
				Conditional(RooArgSet(signal_ss2_time), RooArgSet(obs_time)));
	}
	if (nSignalSSComponents == 3)
	{
		time_pdf_ss3 = std::make_unique<RooProdPdf>(
				"time_pdf_ss3",
				"time_pdf_ss3",
				RooArgSet(timeErrPdf),
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
		auto printVar = [&](const char *title, const RooAbsReal &var)
		{
			const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
			if (rrv && rrv->isConstant())
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
			else if (rrv)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
			else
				tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g", title, var.getVal()));
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
			if (nResolutionComponents >= 2)
				printVar("f_{res1}", ctauFrac1);
			if (nResolutionComponents >= 3)
				printVar("f_{res2}", ctauFrac2);
			if (nResolutionComponents >= 4)
				printVar("f_{res3}", ctauFrac3);
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
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", chiM.first, ndf, pvalue));
	}

	TLine line(ctRange.first, 0.0, ctRange.second, 0.0);
	line.SetLineStyle(2);
	line.Draw("same");

	cLifetime->Print(figName("lifetime_fit"));
	delete cLifetime;

	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	if (timeErrResult)
	{
		cout << "------------------ FIT RESULT FOR TIME ERR ---------------" << endl;
		timeErrResult->Print("v");
	}
	time_result->Print("v");
	cout << "Prompt resolution components used in NP fit: " << nResolutionComponents << endl;
	cout << "SingleSided components used in NP fit: " << nSignalSSComponents << endl;
}
