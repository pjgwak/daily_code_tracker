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
	int nResolutionComponents = 3; // 1~4 components
	
	if (yLow == 1.6f)
	{
		if (ptLow == 1.0f && ptHigh == 2.0f)
		{
			nResolutionComponents = 3;
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
	if (!massTmp || !timeTmp)
	{
		std::cerr << "ERROR: required variables mass/ctau3D are missing in dataset." << std::endl;
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

	TString cutAll = Form(
			"(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
			"(ctau3D >= %g && ctau3D <= %g)",
			ptLow, ptHigh, yLow, yHigh, ctRange.first, ctRange.second);
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
	// Use only lifetime fit (mass fit skipped on purpose).
	// THE CTAU TIME PDF
	// Modeled as a sum of up to four Gaussian resolution components.
	// Use per-event ctau errors as sigma scale factors while keeping the mean scale fixed to 1.
	RooConstVar ctauMeanScale("ctauMeanScale", "ctauMeanScale", 1.0);
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", 0.0, -0.01, 0.01);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale", 0.045, 0.005, 0.20);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", obs_time, ctauTime1Mean, ctauTime1Scale);
	RooFormulaVar ctauTime2Mean("ctauTime2Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta", 0.020, 0.001, 0.30);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", obs_time, ctauTime2Mean, ctauTime2Scale);
	RooFormulaVar ctauTime3Mean("ctauTime3Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta", 0.060, 0.001, 0.60);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", obs_time, ctauTime3Mean, ctauTime3Scale);
	RooFormulaVar ctauTime4Mean("ctauTime4Mean", "@0", RooArgList(ctauTime1Mean));
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta", 1.5, 0.05, 8.0);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", obs_time, ctauTime4Mean, ctauTime4Scale);

	// ------------------------------------------------------------------
	// build ctau fit model
	// ------------------------------------------------------------------
		RooRealVar ctauFrac1Var("ctauFrac1Var", "ctauFrac1Var", 0.75, 0.05, 0.98);
		RooRealVar ctauFrac2Var("ctauFrac2Var", "ctauFrac2Var", 0.60, 0.01, 0.98);
		RooRealVar ctauFrac3Var("ctauFrac3Var", "ctauFrac3Var", 0.30, 0.01, 0.95);

		// THE LIFETIME FIT
		// Fit the prompt resolution as a normalized mixture. Fractions are more stable here than
		// extended component yields, which were collapsing to zero and spoiling the pull.
		RooArgList promptTimePdfList;
		RooArgList promptTimeFracList;
		promptTimePdfList.add(ctauTime1);
		if (nResolutionComponents >= 2)
		{
			promptTimePdfList.add(ctauTime2);
			promptTimeFracList.add(ctauFrac1Var);
		}
		if (nResolutionComponents >= 3)
		{
			promptTimePdfList.add(ctauTime3);
			promptTimeFracList.add(ctauFrac2Var);
		}
		else
		{
			ctauFrac2Var.setVal(0.0);
			ctauFrac2Var.setConstant(true);
			ctauTime3Delta.setConstant(true);
		}
		if (nResolutionComponents >= 4)
		{
			promptTimePdfList.add(ctauTime4);
			promptTimeFracList.add(ctauFrac3Var);
		}
		else
		{
			ctauFrac3Var.setVal(0.0);
			ctauFrac3Var.setConstant(true);
			ctauTime4Delta.setConstant(true);
		}
		if (nResolutionComponents < 2)
		{
			ctauFrac1Var.setVal(0.0);
			ctauFrac1Var.setConstant(true);
			ctauTime2Delta.setConstant(true);
		}
	RooAddPdf time_pdf("time_pdf", "time_pdf", promptTimePdfList, promptTimeFracList, true);
	RooFitResult *time_result = time_pdf.fitTo(*data, Save(), SumW2Error(isWeight));
	time_result->Print();

	const double ctauFrac1 = (nResolutionComponents >= 2) ? ctauFrac1Var.getVal() : 1.0;
	const double ctauFrac2 = (nResolutionComponents >= 3) ? ctauFrac2Var.getVal() : 0.0;
	const double ctauFrac3 = (nResolutionComponents >= 4) ? ctauFrac3Var.getVal() : 0.0;
	const double compFrac1 = (nResolutionComponents == 1) ? 1.0 : ctauFrac1;
	const double compFrac2 = (nResolutionComponents >= 2) ? ((nResolutionComponents == 2) ? (1.0 - ctauFrac1) : (1.0 - ctauFrac1) * ctauFrac2) : 0.0;
	const double compFrac3 = (nResolutionComponents >= 3) ? ((nResolutionComponents == 3) ? (1.0 - ctauFrac1) * (1.0 - ctauFrac2) : (1.0 - ctauFrac1) * (1.0 - ctauFrac2) * ctauFrac3) : 0.0;
	const double compFrac4 = (nResolutionComponents >= 4) ? ((1.0 - ctauFrac1) * (1.0 - ctauFrac2) * (1.0 - ctauFrac3)) : 0.0;
	const double nEntries = data->numEntries();
	const TString resolutionFileName = TString::Format("%s/ctau_resolution_%s.root", resultDir.Data(), figTag.Data());
	{
		TFile resolutionFile(resolutionFileName, "RECREATE");
		if (!resolutionFile.IsZombie())
		{
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
	time_pdf.plotOn(timefitplot, LineColor(kBlack), LineWidth(2), Name("model"));
	time_pdf_component1.plotOn(
			timefitplot,
			LineColor(kBlue + 1),
			LineStyle(kDashed),
			LineWidth(2),
			Normalization(compFrac1 * nEntries, RooAbsReal::NumEvent),
			Name("ctau_component_1"));
	if (nResolutionComponents >= 2)
	{
		time_pdf_component2.plotOn(
				timefitplot,
				LineColor(kRed + 1),
				LineStyle(kDashed),
				LineWidth(2),
				Normalization(compFrac2 * nEntries, RooAbsReal::NumEvent),
				Name("ctau_component_2"));
	}
		if (nResolutionComponents >= 3)
		{
			time_pdf_component3.plotOn(
					timefitplot,
					LineColor(kMagenta + 1),
					LineStyle(kDashed),
					LineWidth(2),
					Normalization(compFrac3 * nEntries, RooAbsReal::NumEvent),
					Name("ctau_component_3"));
		}
		if (nResolutionComponents >= 4)
		{
			time_pdf_component4.plotOn(
					timefitplot,
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
			printVar("f_{1}", ctauFrac1Var);
			if (nResolutionComponents >= 2)
				printVar("f_{2}", ctauFrac2Var);
				printVar("mean_{1}", ctauTime1Mean);
				printVar("#sigma_{1}", ctauTime1Scale);
				if (nResolutionComponents >= 2)
				{
					printVar("mean_{2}", ctauTime2Mean);
					printVar("#Delta #sigma_{21}", ctauTime2Delta);
				}
				if (nResolutionComponents >= 3)
				{
					printVar("f_{3}", ctauFrac3Var);
					printVar("mean_{3}", ctauTime3Mean);
					printVar("#Delta #sigma_{32}", ctauTime3Delta);
				}
				if (nResolutionComponents >= 4)
				{
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
	cout << "------------------ FIT RESULT SUMMARY --------------------" << endl;
	cout << "Prompt resolution components used in PR fit: " << nResolutionComponents << endl;
	cout << Form("Prompt ctau chi2/ndf: %.1f/%d = %.3f, p=%.3g", chiM.first, ndf, ndf > 0 ? chiM.first / ndf : 0.0, pvalue) << endl;
	cout << "------------------ FIT RESULT FOR TIME ONLY --------------" << endl;
	time_result->Print();
	const TString figLifetime = figName("lifetime_fit");
	std::cout << "[FIG] ctau_pr lifetime fit : " << figLifetime << std::endl;
	// time_pdf.Print("V");
}

// //
// 	// THE LIFETIME FIT: decompose into two resolved components for plotting
// 	RooProdPdf time_pdf_component1(
// 		"time_pdf_component1",
// 		"time_pdf_component1",
// 		Conditional(RooArgSet(ctauTime1), RooArgSet(obs_time))
// 	);
// 	RooProdPdf time_pdf_component2(
// 		"time_pdf_component2",
// 		"time_pdf_component2",
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
