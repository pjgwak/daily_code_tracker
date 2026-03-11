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
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
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
#include "TMath.h"
#include "TParameter.h"
#include "TString.h"
#include "RooHist.h"
#include <algorithm>
#include <cmath>
#include <iostream>
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

void mass(float ptLow = 1, float ptHigh = 2, float yLow = 1.6, float yHigh = 2.4)
{
	// ------------------------------------------------------------------
	// model control
	// ------------------------------------------------------------------
	// Choose how many mass components to use in each (pt, y) bin.
	int nSignalCBComponents = 2;
	int nSignalGaussComponents = 1;
	int nBkgExpComponents = 0;
	int nBkgChebyOrder = 3;
	if (yLow == 0.0f)
	{
		if (ptLow == 200.0f && ptHigh == 350.0f) //dummy
		{
			nSignalCBComponents = 1;    // 0~2, automatic
			nSignalGaussComponents = 1; // 0~2, automatic
			nBkgExpComponents = 0;      // 0~1
			nBkgChebyOrder = 2;         // 0~6
		}
	}
	nSignalCBComponents = std::clamp(nSignalCBComponents, 0, 2);
	nSignalGaussComponents = std::clamp(nSignalGaussComponents, 0, 2);
	nBkgExpComponents = std::clamp(nBkgExpComponents, 0, 1);
	nBkgChebyOrder = std::clamp(nBkgChebyOrder, 0, 6);
	if (nSignalCBComponents + nSignalGaussComponents <= 0)
	{
		std::cerr << "ERROR: at least one signal component is required." << std::endl;
		return;
	}
	if (nBkgExpComponents > 0 && nBkgChebyOrder > 0)
	{
		std::cerr << "ERROR: choose either exponential or Chebychev background." << std::endl;
		return;
	}
	if (nBkgExpComponents + nBkgChebyOrder <= 0)
	{
		std::cerr << "ERROR: no background model is configured." << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
	bool isWeight = false;
	
	const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
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
	if (!dataBasic || dataBasic->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	TString finalCut = cutBasic;
	auto *ctau3DTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
	auto *ctau3DErrTmp = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DErr"));
	if (!ctau3DTmp || !ctau3DErrTmp) {
		std::cerr << "ERROR: required variables ctau3D/ctau3DErr are missing in dataset." << std::endl;
		return;
	}

	auto quantileRange = [&](RooDataSet &ds, RooRealVar &var) {
		std::vector<double> vals;
		vals.reserve(ds.numEntries());
		for (int i = 0; i < ds.numEntries(); ++i) {
			ds.get(i);
			const double v = var.getVal();
			if (!std::isfinite(v)) continue;
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
		return std::make_pair(qAt(0.001), qAt(0.999));
	};

	const auto ctau3DRange = quantileRange(*dataBasic, *ctau3DTmp);
	const auto ctau3DErrRange = quantileRange(*dataBasic, *ctau3DErrTmp);
	finalCut = Form(
		"%s && (ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)",
		cutBasic.Data(), ctau3DRange.first, ctau3DRange.second, ctau3DErrRange.first, ctau3DErrRange.second
	);

	auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(finalCut)));
	if (!dataSel || dataSel->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after selection: " << finalCut << std::endl;
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
	const TString figDir = TString::Format("figs/%s/mass", yTag.Data());
	const TString resultDir = TString::Format("roots/%s/mass", yTag.Data());
	const TString figTag = yTag + "_" + ptTag;
	const TString mcResultDir = TString::Format("roots/%s/mc_mass", yTag.Data());
	const TString mcModelFileName = TString::Format("%s/mc_mass_model_%s.root", mcResultDir.Data(), figTag.Data());
	const TString modelFileName = TString::Format("%s/mass_model_%s.root", resultDir.Data(), figTag.Data());
	auto figName = [&](const char *name)
	{
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(resultDir, true);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	TFile *mcModelFile = TFile::Open(mcModelFileName);
	if (!mcModelFile || mcModelFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open MC mass model file: " << mcModelFileName << std::endl;
		return;
	}
	auto readIntParam = [&](const char *name, int fallback)
	{
		auto *param = dynamic_cast<TParameter<int> *>(mcModelFile->Get(name));
		return param ? param->GetVal() : fallback;
	};
	auto readDoubleParam = [&](const char *name, double fallback)
	{
		auto *param = dynamic_cast<TParameter<double> *>(mcModelFile->Get(name));
		return param ? param->GetVal() : fallback;
	};
	nSignalGaussComponents = std::clamp(readIntParam("nSignalGaussComponents", nSignalGaussComponents), 0, 2);
	nSignalCBComponents = std::clamp(readIntParam("nSignalCBComponents", nSignalCBComponents), 0, 2);
	if (nSignalCBComponents + nSignalGaussComponents <= 0)
	{
		std::cerr << "ERROR: invalid signal component configuration in " << mcModelFileName << std::endl;
		return;
	}

	RooRealVar &obs_mass = *static_cast<RooRealVar *>(dataSel->get()->find("mass"));
	obs_mass.SetTitle("mass");
	obs_mass.setUnit("GeV/c^{2}");
	const int massPlotBins = 200;

	RooDataSet *data = dataSel.get();
	auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
	{
		return fr ? fr->findObject(n) : nullptr;
	};

	/**********************THE MASS FIT********************************
	 * Some hints: The signal mass is a simple gaussian, you can use a
	 * RooGaussian PDF for the signal model. The background is an O(1)
	 * polynomial with a negative gradient. Polynomials are
	 * surprisingly hard to fit, so to help you out I'll tell you that
	 * the zero order term is 5500 and should be _fixed_, allowing only
	 * the gradient to float. Limits are important: Is it worth trying
	 * to fit the first order term on a range >0? Try the range -2.0 -> 0.0
	 */

	RooRealVar signal_mass_mean("signal_mass_mean", "signal_mass_mean", 3.096, 3.05, 3.15, "GeV/c^{2}");
	RooRealVar signal_mass_sigma("signal_mass_sigma", "signal_mass_sigma", 0.03, 0.005, 0.080, "GeV/c^{2}");
	RooRealVar signal_mass_sigma2("signal_mass_sigma2", "signal_mass_sigma2", 0.06, 0.010, 0.120, "GeV/c^{2}");
	RooRealVar signal_mass_cb_sigma("signal_mass_cb_sigma", "signal_mass_cb_sigma", 0.035, 0.008, 0.080, "GeV/c^{2}");
	RooRealVar signal_mass_cb_sigma2("signal_mass_cb_sigma2", "signal_mass_cb_sigma2", 0.055, 0.010, 0.120, "GeV/c^{2}");
	RooRealVar signal_mass_cb_alpha("signal_mass_cb_alpha", "signal_mass_cb_alpha", 1.5, -5.0, 5.0);
	RooRealVar signal_mass_cb_n("signal_mass_cb_n", "signal_mass_cb_n", 3.0, 1.0, 8.0);
	RooRealVar signal_mass_cb_alpha2("signal_mass_cb_alpha2", "signal_mass_cb_alpha2", 2.0, -5.0, 5.0);
	RooRealVar signal_mass_frac1("signal_mass_frac1", "signal_mass_frac1", 0.65, 0.0, 1.0);
	RooRealVar signal_mass_frac2("signal_mass_frac2", "signal_mass_frac2", 0.20, 0.0, 1.0);
	RooRealVar signal_mass_frac3("signal_mass_frac3", "signal_mass_frac3", 0.10, 0.0, 1.0);

	signal_mass_mean.setVal(readDoubleParam("signal_mass_mean", signal_mass_mean.getVal()));
	signal_mass_sigma.setVal(readDoubleParam("signal_mass_sigma", signal_mass_sigma.getVal()));
	signal_mass_sigma2.setVal(readDoubleParam("signal_mass_sigma2", signal_mass_sigma2.getVal()));
	signal_mass_cb_sigma.setVal(readDoubleParam("signal_mass_cb_sigma", signal_mass_cb_sigma.getVal()));
	signal_mass_cb_sigma2.setVal(readDoubleParam("signal_mass_cb_sigma2", signal_mass_cb_sigma2.getVal()));
	signal_mass_cb_alpha.setVal(readDoubleParam("signal_mass_cb_alpha", signal_mass_cb_alpha.getVal()));
	signal_mass_cb_alpha2.setVal(readDoubleParam("signal_mass_cb_alpha2", signal_mass_cb_alpha2.getVal()));
	signal_mass_cb_n.setVal(readDoubleParam("signal_mass_cb_n", signal_mass_cb_n.getVal()));
	signal_mass_frac1.setVal(readDoubleParam("signal_mass_frac1", signal_mass_frac1.getVal()));
	signal_mass_frac2.setVal(readDoubleParam("signal_mass_frac2", signal_mass_frac2.getVal()));
	signal_mass_frac3.setVal(readDoubleParam("signal_mass_frac3", signal_mass_frac3.getVal()));
	signal_mass_cb_alpha.setConstant(true);
	signal_mass_cb_alpha2.setConstant(true);
	signal_mass_cb_n.setConstant(true);

	RooGaussian signal_mass_gaus("signal_mass_gaus", "signal_mass_gaus", obs_mass, signal_mass_mean, signal_mass_sigma);
	RooGaussian signal_mass_gaus2("signal_mass_gaus2", "signal_mass_gaus2", obs_mass, signal_mass_mean, signal_mass_sigma2);
	RooCBShape signal_mass_cb("signal_mass_cb", "signal_mass_cb", obs_mass, signal_mass_mean, signal_mass_cb_sigma, signal_mass_cb_alpha, signal_mass_cb_n);
	RooCBShape signal_mass_cb2("signal_mass_cb2", "signal_mass_cb2", obs_mass, signal_mass_mean, signal_mass_cb_sigma2, signal_mass_cb_alpha2, signal_mass_cb_n);

	std::vector<RooAbsPdf *> signalComponents;
	if (nSignalGaussComponents >= 1)
		signalComponents.push_back(&signal_mass_gaus);
	if (nSignalCBComponents >= 1)
		signalComponents.push_back(&signal_mass_cb);
	if (nSignalGaussComponents >= 2)
		signalComponents.push_back(&signal_mass_gaus2);
	if (nSignalCBComponents >= 2)
		signalComponents.push_back(&signal_mass_cb2);

	RooArgList signalPdfList;
	for (RooAbsPdf *pdf : signalComponents)
		signalPdfList.add(*pdf);
	RooArgList signalFracList;
	if (signalComponents.size() >= 2)
		signalFracList.add(signal_mass_frac1);
	if (signalComponents.size() >= 3)
		signalFracList.add(signal_mass_frac2);
	if (signalComponents.size() >= 4)
		signalFracList.add(signal_mass_frac3);
	std::unique_ptr<RooAbsPdf> signal_mass_owned;
	RooAbsPdf *signal_mass_pdf = nullptr;
	if (signalComponents.size() == 1)
	{
		signal_mass_pdf = signalComponents.front();
	}
	else
	{
		signal_mass_owned = std::make_unique<RooAddPdf>("signal_mass", "signal_mass", signalPdfList, signalFracList, true);
		signal_mass_pdf = signal_mass_owned.get();
	}

	RooRealVar Nsig("Nsig", "Nsig", std::max(1.0, 0.8 * data->numEntries()), 0.0, data->numEntries());
	RooRealVar Nbkg("Nbkg", "Nbkg", std::max(1.0, 0.2 * data->numEntries()), 0.0, data->numEntries());
	RooRealVar bkg_mass_p1("bkg_mass_p1", "bkg_mass_p1", 0.0, -0.4, 0.4);
	RooRealVar bkg_mass_p2("bkg_mass_p2", "bkg_mass_p2", 0.0, -0.2, 0.2);
	RooRealVar bkg_mass_p3("bkg_mass_p3", "bkg_mass_p3", 0.0, -0.2, 0.2);
	RooRealVar bkg_mass_p4("bkg_mass_p4", "bkg_mass_p4", 0.0, -0.2, 0.2);
	RooRealVar bkg_mass_p5("bkg_mass_p5", "bkg_mass_p5", 0.0, -0.2, 0.2);
	RooRealVar bkg_mass_p6("bkg_mass_p6", "bkg_mass_p6", 0.0, -0.2, 0.2);
	RooRealVar bkg_mass_lambda("bkg_mass_lambda", "bkg_mass_lambda", -1.0, -10.0, -1e-4);
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
	std::unique_ptr<RooAbsPdf> bkg_mass;
	if (nBkgExpComponents == 1)
		bkg_mass = std::make_unique<RooExponential>("bkg_mass", "bkg_mass", obs_mass, bkg_mass_lambda);
	else if (nBkgChebyOrder >= 1)
		bkg_mass = std::make_unique<RooChebychev>("bkg_mass", "bkg_mass", obs_mass, chebyCoeffList);
	else
		bkg_mass = std::make_unique<RooExponential>("bkg_mass", "bkg_mass", obs_mass, bkg_mass_lambda);
	const TString bkgLegendLabel = nBkgExpComponents == 1 ? "Background (E)" : TString::Format("Background (%d)", nBkgChebyOrder);

	std::unique_ptr<RooAbsPdf> mass_pdf;
	mass_pdf = std::make_unique<RooAddPdf>("mass_pdf", "mass_pdf", RooArgList(*signal_mass_pdf, *bkg_mass), RooArgList(Nsig, Nbkg));

	std::vector<std::unique_ptr<RooRealVar>> constraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> constraintPdfs;
	RooArgSet constraints;
	auto addConstraint = [&](const char *name, RooRealVar &var, bool freeThis)
	{
		if (freeThis)
			return;
		const double central = readDoubleParam(name, var.getVal());
		const double error = readDoubleParam((TString(name) + "_err").Data(), 0.0);
		if (!std::isfinite(error) || error <= 0.0)
			return;
		var.setVal(central);
		const TString meanName = TString::Format("%s_mc_mean", name);
		const TString sigmaName = TString::Format("%s_mc_sigma", name);
		const TString pdfName = TString::Format("%s_constraint", name);
		constraintConsts.push_back(std::make_unique<RooRealVar>(meanName, meanName, central));
		constraintConsts.back()->setConstant(true);
		const double sigmaVal = std::max(error, 1e-6);
		constraintConsts.push_back(std::make_unique<RooRealVar>(sigmaName, sigmaName, sigmaVal, 1e-9, std::max(10.0 * sigmaVal, 1e-8)));
		constraintConsts.back()->setConstant(true);
		constraintPdfs.push_back(std::make_unique<RooGaussian>(pdfName, pdfName, var, *constraintConsts[constraintConsts.size() - 2], *constraintConsts.back()));
		constraints.add(*constraintPdfs.back());
	};
	bool firstSigmaFreed = false;
	addConstraint("signal_mass_mean", signal_mass_mean, true);
	if (nSignalGaussComponents >= 1)
	{
		addConstraint("signal_mass_sigma", signal_mass_sigma, !firstSigmaFreed);
		firstSigmaFreed = true;
	}
	if (nSignalCBComponents >= 1)
	{
		addConstraint("signal_mass_cb_sigma", signal_mass_cb_sigma, !firstSigmaFreed);
		firstSigmaFreed = true;
	}
	if (nSignalGaussComponents >= 2)
		addConstraint("signal_mass_sigma2", signal_mass_sigma2, false);
	if (nSignalCBComponents >= 2)
	{
		addConstraint("signal_mass_cb_sigma2", signal_mass_cb_sigma2, false);
	}
	if (signalComponents.size() >= 2)
		addConstraint("signal_mass_frac1", signal_mass_frac1, false);
	if (signalComponents.size() >= 3)
		addConstraint("signal_mass_frac2", signal_mass_frac2, false);
	if (signalComponents.size() >= 4)
		addConstraint("signal_mass_frac3", signal_mass_frac3, false);

	std::unique_ptr<RooFitResult> mass_result;
	if (constraints.getSize() > 0)
		mass_result.reset(mass_pdf->fitTo(*data, Extended(), Save(), PrintLevel(-1), SumW2Error(isWeight), ExternalConstraints(constraints)));
	else
		mass_result.reset(mass_pdf->fitTo(*data, Extended(), Save(), PrintLevel(-1), SumW2Error(isWeight)));
	mass_result->Print();

	{
		TFile outFile(modelFileName, "RECREATE");
		if (!outFile.IsZombie())
		{
			TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
			TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
			TParameter<int>("nBkgExpComponents", nBkgExpComponents).Write();
			TParameter<int>("nBkgChebyOrder", nBkgChebyOrder).Write();
			auto writeVar = [&](const char *name, const RooRealVar &var)
			{
				TParameter<double>(name, var.getVal()).Write();
				TParameter<double>(TString(name) + "_err", std::max(var.getError(), 1e-6)).Write();
			};
			writeVar("Nsig", Nsig);
			writeVar("Nbkg", Nbkg);
			writeVar("signal_mass_mean", signal_mass_mean);
			if (nSignalGaussComponents >= 1)
				writeVar("signal_mass_sigma", signal_mass_sigma);
			if (nSignalCBComponents >= 1)
			{
				writeVar("signal_mass_cb_sigma", signal_mass_cb_sigma);
				writeVar("signal_mass_cb_alpha", signal_mass_cb_alpha);
				writeVar("signal_mass_cb_n", signal_mass_cb_n);
			}
			if (nSignalGaussComponents >= 2)
				writeVar("signal_mass_sigma2", signal_mass_sigma2);
			if (nSignalCBComponents >= 2)
			{
				writeVar("signal_mass_cb_sigma2", signal_mass_cb_sigma2);
				writeVar("signal_mass_cb_alpha2", signal_mass_cb_alpha2);
			}
			if (signalComponents.size() >= 2)
				writeVar("signal_mass_frac1", signal_mass_frac1);
			if (signalComponents.size() >= 3)
				writeVar("signal_mass_frac2", signal_mass_frac2);
			if (signalComponents.size() >= 4)
				writeVar("signal_mass_frac3", signal_mass_frac3);
			if (nBkgExpComponents == 1)
				writeVar("bkg_mass_lambda", bkg_mass_lambda);
			if (nBkgChebyOrder >= 1)
				writeVar("bkg_mass_p1", bkg_mass_p1);
			if (nBkgChebyOrder >= 2)
				writeVar("bkg_mass_p2", bkg_mass_p2);
			if (nBkgChebyOrder >= 3)
				writeVar("bkg_mass_p3", bkg_mass_p3);
			if (nBkgChebyOrder >= 4)
				writeVar("bkg_mass_p4", bkg_mass_p4);
			if (nBkgChebyOrder >= 5)
				writeVar("bkg_mass_p5", bkg_mass_p5);
			if (nBkgChebyOrder >= 6)
				writeVar("bkg_mass_p6", bkg_mass_p6);
			if (mass_result)
				mass_result->Write("fit_result");
		}
	}

	// ------------------------------------------------------------------
	// draw mass fit
	// ------------------------------------------------------------------
	TCanvas *cMass = new TCanvas("cMass", "cMass", 800, 800);
	TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
	pad1->SetBottomMargin(0.00001);
	pad1->SetTopMargin(0.08);
	pad1->Draw();
	pad1->cd();

	RooPlot *massplot = obs_mass.frame(Title(""));
	if (isWeight)
		data->plotOn(massplot, Binning(massPlotBins), Name("data"), DataError(RooAbsData::SumW2));
	else
		data->plotOn(massplot, Binning(massPlotBins), Name("data"));
	mass_pdf->plotOn(massplot, LineColor(kBlack), LineWidth(2), Name("model"));
	if (nSignalGaussComponents >= 1)
		mass_pdf->plotOn(massplot, Components(signal_mass_gaus), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Name("signal_gaus_component"));
	if (nSignalCBComponents >= 1)
		mass_pdf->plotOn(massplot, Components(signal_mass_cb), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("signal_cb_component"));
	if (nSignalGaussComponents >= 2)
		mass_pdf->plotOn(massplot, Components(signal_mass_gaus2), LineColor(kAzure + 2), LineStyle(kDashed), LineWidth(2), Name("signal_gaus2_component"));
	if (nSignalCBComponents >= 2)
		mass_pdf->plotOn(massplot, Components(signal_mass_cb2), LineColor(kOrange + 7), LineStyle(kDashed), LineWidth(2), Name("signal_cb2_component"));
	mass_pdf->plotOn(massplot, Components(*bkg_mass), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Name("bkg_component"));

	double ymax = -1e300;
	if (auto *hdata = dynamic_cast<RooHist *>(massplot->getHist("data")))
	{
		for (int i = 0; i < hdata->GetN(); ++i)
		{
			double xval = 0.0, yval = 0.0;
			hdata->GetPoint(i, xval, yval);
			if (yval > ymax)
				ymax = yval;
		}
	}
	if (ymax > 0.0 && ymax < 1e300)
		massplot->SetMaximum(ymax * 1.8);
	massplot->GetYaxis()->SetTitle("Events");
	massplot->GetYaxis()->SetTitleOffset(1.6);
	massplot->GetXaxis()->SetTitle("");
	massplot->Draw("e");

	TLegend leg(0.50, 0.70, 0.74, 0.89);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextSize(0.03);
	if (auto *o = findObj(massplot, "data"))
		leg.AddEntry(o, "Data", "lep");
	if (auto *o = findObj(massplot, "model"))
		leg.AddEntry(o, "Fit", "l");
	if (auto *o = findObj(massplot, "signal_gaus_component"))
		leg.AddEntry(o, "G1", "l");
	if (auto *o = findObj(massplot, "signal_cb_component"))
		leg.AddEntry(o, "CB1", "l");
	if (nSignalGaussComponents >= 2)
		if (auto *o = findObj(massplot, "signal_gaus2_component"))
			leg.AddEntry(o, "G2", "l");
	if (nSignalCBComponents >= 2)
		if (auto *o = findObj(massplot, "signal_cb2_component"))
			leg.AddEntry(o, "CB2", "l");
	if (auto *o = findObj(massplot, "bkg_component"))
		leg.AddEntry(o, bkgLegendLabel.Data(), "l");
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
		const int status = mass_result ? mass_result->status() : -1;
		int hesse = -1;
		if (mass_result)
		{
			for (UInt_t i = 0, n = mass_result->numStatusHistory(); i < n; ++i)
			{
				const char *lab = mass_result->statusLabelHistory(i);
				if (lab && TString(lab) == "HESSE")
				{
					hesse = mass_result->statusCodeHistory(i);
					break;
				}
			}
		}
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.74, y0 = 0.87, dy = -0.045;
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
		printVar("N_{sig}", Nsig);
		printVar("N_{bkg}", Nbkg);
		printVar("#mu", signal_mass_mean);
			if (nSignalGaussComponents >= 1)
				printVar("#sigma_{G1}", signal_mass_sigma);
			if (nSignalCBComponents >= 1)
			{
				printVar("#sigma_{CB1}", signal_mass_cb_sigma);
				printVar("#alpha_{CB1}", signal_mass_cb_alpha);
				printVar("n_{CB}", signal_mass_cb_n);
			}
			if (nSignalGaussComponents >= 2)
				printVar("#sigma_{G2}", signal_mass_sigma2);
			if (nSignalCBComponents >= 2)
			{
				printVar("#sigma_{CB2}", signal_mass_cb_sigma2);
				printVar("#alpha_{CB2}", signal_mass_cb_alpha2);
			}
		if (signalComponents.size() >= 2)
			printVar("f_{1}", signal_mass_frac1);
		if (signalComponents.size() >= 3)
			printVar("f_{2}", signal_mass_frac2);
		if (signalComponents.size() >= 4)
			printVar("f_{3}", signal_mass_frac3);
		if (nBkgExpComponents == 1)
			printVar("#lambda_{E}", bkg_mass_lambda);
		if (nBkgChebyOrder >= 1)
			printVar("p_{1}", bkg_mass_p1);
		if (nBkgChebyOrder >= 2)
			printVar("p_{2}", bkg_mass_p2);
		if (nBkgChebyOrder >= 3)
			printVar("p_{3}", bkg_mass_p3);
		if (nBkgChebyOrder >= 4)
			printVar("p_{4}", bkg_mass_p4);
		if (nBkgChebyOrder >= 5)
			printVar("p_{5}", bkg_mass_p5);
		if (nBkgChebyOrder >= 6)
			printVar("p_{6}", bkg_mass_p6);
	}

	// ------------------------------------------------------------------
	// draw pull and print summary
	// ------------------------------------------------------------------
	cMass->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
	pad2->SetTopMargin(0.00001);
	pad2->SetBottomMargin(0.4);
	pad2->Draw();
	pad2->cd();

	RooPlot *fpull = obs_mass.frame(Title(""));
	RooHist *hpull = massplot->pullHist("data", "model");
	if (hpull)
		fpull->addPlotable(hpull, "P");
	fpull->GetYaxis()->SetTitle("Pull");
	fpull->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
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
	int npar = mass_result ? mass_result->floatParsFinal().getSize() : 0;
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

	TLine line(obs_mass.getMin(), 0.0, obs_mass.getMax(), 0.0);
	line.SetLineStyle(2);
	line.Draw("same");

	cMass->Print(figName("mass_fit"));
	delete cMass;

	std::cout << "------------------ FIT RESULT FOR MASS ONLY --------------" << std::endl;
	mass_result->Print("v");
}
