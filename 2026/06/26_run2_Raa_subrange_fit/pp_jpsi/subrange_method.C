#include "RooAddModel.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooCrystalBall.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGlobalFunc.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooResolutionModel.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TParameter.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1D.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

using namespace RooFit;

namespace
{
const double kCtauConstraintScale = 3;
const bool kDisableCtauConstraints = false;
const bool kFloatSignalSSFractions = false;
const double kCtauFitMinOverride = -0.2;
const char *kDataRoot = "/data/users/pjgwak/work/raa_pb18/run2_raa_pbpb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_260303_root634.root";
const char *kDatasetName = "dataset";
const char *kReferenceDir = "/data/users/pjgwak/work/daily_code_tracker/2026/06/26_run2_Raa_subrange_fit/pp_jpsi";
const char *kAccCut = "(((abs(eta1) <= 1.2) && (pt1 >= 3.5)) || ((abs(eta2) <= 1.2) && (pt2 >= 3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)))";

TString format_tag(double value)
{
	return TString::Format("%.2f", value);
}

TString bin_tag(float ptLow, float ptHigh, float yLow, float yHigh)
{
	return TString::Format("y%s_%s_pt%s_%s",
		format_tag(yLow).Data(), format_tag(yHigh).Data(),
		format_tag(ptLow).Data(), format_tag(ptHigh).Data());
}

TString y_tag(float yLow, float yHigh)
{
	return TString::Format("y%s_%s", format_tag(yLow).Data(), format_tag(yHigh).Data());
}

double read_double(TFile *file, const char *name, double fallback)
{
	if (!file)
		return fallback;
	if (auto *param = dynamic_cast<TParameter<double> *>(file->Get(name)))
		return param->GetVal();
	if (auto *var = dynamic_cast<RooRealVar *>(file->Get(name)))
		return var->getVal();
	return fallback;
}

int read_int(TFile *file, const char *name, int fallback)
{
	if (!file)
		return fallback;
	if (auto *param = dynamic_cast<TParameter<int> *>(file->Get(name)))
		return param->GetVal();
	return fallback;
}

double finite_or(double value, double fallback)
{
	return std::isfinite(value) ? value : fallback;
}

std::pair<double, double> quantile_range(RooDataSet &ds, RooRealVar &var, double qLo = 0.001, double qHi = 0.999)
{
	std::vector<double> vals;
	vals.reserve(ds.numEntries());
	for (int i = 0; i < ds.numEntries(); ++i)
	{
		ds.get(i);
		const double v = var.getVal();
		if (std::isfinite(v))
			vals.push_back(v);
	}
	if (vals.empty())
		return {var.getMin(), var.getMax()};
	std::sort(vals.begin(), vals.end());
	auto qAt = [&](double q) {
		long long idx = llround(std::clamp(q, 0.0, 1.0) * (vals.size() - 1));
		idx = std::clamp(idx, 0LL, static_cast<long long>(vals.size()) - 1);
		return vals[static_cast<size_t>(idx)];
	};
	return {qAt(qLo), qAt(qHi)};
}

std::pair<double, double> fit_range_around(double center, double hardLow, double hardHigh, double relWidth, double absWidth)
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
	return {low, high};
}

std::pair<double, int> chi2_from_pull(RooHist *hpull)
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

std::vector<double> build_strong_adaptive_ctau_edges(double ctMin, double ctMax, double ctStep)
{
	std::vector<double> edges;
	if (ctStep <= 0.0 || ctMax <= ctMin)
		return edges;

	auto addEdge = [&](double x) {
		if (x < ctMin - 1e-9 || x > ctMax + 1e-9)
			return;
		if (edges.empty() || std::abs(edges.back() - x) > 1e-9)
			edges.push_back(x);
	};
	auto appendFixedBins = [&](double start, double stop, int nBins) {
		if (stop <= start || nBins <= 0)
			return;
		addEdge(start);
		for (int i = 1; i < nBins; ++i)
			addEdge(start + (stop - start) * static_cast<double>(i) / static_cast<double>(nBins));
		addEdge(stop);
	};
	auto appendBand = [&](double start, double stop, double targetWidth) {
		start = std::max(start, ctMin);
		stop = std::min(stop, ctMax);
		if (stop <= start || targetWidth <= 0.0)
			return;
		const int nBins = std::max(1, static_cast<int>(std::ceil((stop - start) / targetWidth - 1e-9)));
		appendFixedBins(start, stop, nBins);
	};

	const double tailWide = std::max(4.0, 2.0 * ctStep);
	const double tailMid = std::max(2.0, ctStep);
	const double tailNear = std::max(1.0, 0.5 * ctStep);
	appendBand(ctMin, -20.0, tailWide);
	appendBand(-20.0, -10.0, tailMid);
	appendBand(-10.0, -5.0, tailNear);
	appendBand(-5.0, -3.0, 0.50);
	appendBand(-3.0, -2.5, 0.20);
	appendBand(-2.5, -2.0, 0.10);
	appendBand(-2.0, -1.45, 0.075);
	appendBand(-1.45, -1.05, 0.050);
	appendBand(-1.05, -0.75, 0.030);
	appendBand(-0.75, -0.50, 0.020);
	appendBand(-0.50, -0.30, 0.012);
	appendBand(-0.30, -0.15, 0.008);
	appendBand(-0.15, -0.06, 0.005);
	appendBand(-0.06, 0.06, 0.0025);
	appendBand(0.06, 0.15, 0.005);
	appendBand(0.15, 0.30, 0.008);
	appendBand(0.30, 0.50, 0.012);
	appendBand(0.50, 0.75, 0.020);
	appendBand(0.75, 1.05, 0.030);
	appendBand(1.05, 1.45, 0.050);
	appendBand(1.45, 2.00, 0.075);
	appendBand(2.00, 2.50, 0.10);
	appendBand(2.50, 3.00, 0.20);
	appendBand(3.00, 5.00, 0.50);
	appendBand(5.00, 10.00, tailNear);
	appendBand(10.00, 20.00, tailMid);
	appendBand(20.00, ctMax, tailWide);

	return edges;
}

int adaptive_mass_bins(int entries, int minBins = 18, int maxBins = 90, double targetEntriesPerBin = 35.0)
{
	if (entries <= 0)
		return minBins;
	const int bins = static_cast<int>(std::lround(static_cast<double>(entries) / targetEntriesPerBin));
	return std::clamp(bins, minBins, maxBins);
}

void apply_logy_auto_range(RooPlot *plot, const char *histName, double topScale = 20.0, double bottomScale = 1e3)
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
		double x = 0.0;
		double y = 0.0;
		hist->GetPoint(i, x, y);
		if (!std::isfinite(y) || y <= 0.0)
			continue;
		const double yHigh = y + std::max(0.0, hist->GetErrorYhigh(i));
		const double yLow = y - std::max(0.0, hist->GetErrorYlow(i));
		peak = std::max(peak, std::isfinite(yHigh) ? yHigh : y);
		minPositive = std::min(minPositive, y);
		if (std::isfinite(yLow) && yLow > 0.0)
			minPositive = std::min(minPositive, yLow);
	}
	if (peak <= 0.0)
		return;

	double ymin = peak / bottomScale;
	if (std::isfinite(minPositive) && minPositive > 0.0)
		ymin = std::min(ymin, 0.5 * minPositive);
	ymin = std::max(ymin, 1e-3);
	const double ymax = std::max(peak * topScale, ymin * 10.0);
	plot->SetMinimum(ymin);
	plot->SetMaximum(ymax);
}

bool fit_is_good(const RooFitResult *result, int minCovQual = 2)
{
	return result && result->status() == 0 && result->covQual() >= minCovQual;
}

std::unique_ptr<RooFitResult> fit_mass_with_retry(RooAbsPdf &pdf, RooDataSet &data, int &attemptUsed,
	const RooArgSet *constraints = nullptr)
{
	auto runFit = [&](int strategy, bool useMinuit2) {
		if (constraints && constraints->getSize() > 0)
		{
			if (useMinuit2)
				return pdf.fitTo(data, Extended(), Save(), Strategy(strategy), PrintLevel(-1),
					Minimizer("Minuit2", "migrad"), ExternalConstraints(*constraints));
			return pdf.fitTo(data, Extended(), Save(), Strategy(strategy), PrintLevel(-1), ExternalConstraints(*constraints));
		}
		if (useMinuit2)
			return pdf.fitTo(data, Extended(), Save(), Strategy(strategy), PrintLevel(-1), Minimizer("Minuit2", "migrad"));
		return pdf.fitTo(data, Extended(), Save(), Strategy(strategy), PrintLevel(-1));
	};

	attemptUsed = 1;
	auto result = std::unique_ptr<RooFitResult>(runFit(2, true));
	if (fit_is_good(result.get()))
		return result;

	std::cout << Form("[INFO] mass fit retry 1: status=%d covQual=%d, retrying default minimizer",
		result ? result->status() : -999, result ? result->covQual() : -999) << std::endl;
	attemptUsed = 2;
	result.reset(runFit(2, false));
	if (fit_is_good(result.get()))
		return result;

	std::cout << Form("[INFO] mass fit retry 2: status=%d covQual=%d, retrying lower strategy",
		result ? result->status() : -999, result ? result->covQual() : -999) << std::endl;
	attemptUsed = 3;
	result.reset(runFit(1, false));
	return result;
}

struct SliceResult
{
	double ctLow = 0.0;
	double ctHigh = 0.0;
	int entries = 0;
	double nsig = 0.0;
	double nsigErr = 0.0;
	double nbkg = 0.0;
	double nbkgErr = 0.0;
	double nsigSignificance = 0.0;
	int signalWindowEntries = 0;
	int shapeDrivenSignal = 0;
	int usedInCtauFit = 1;
	int status = -99;
	int covQual = -99;
	int fitAttempt = 0;
};

SliceResult fit_mass_slice(RooDataSet &inputData,
	RooRealVar &massIn,
	float ptLow,
	float ptHigh,
	float yLow,
	float yHigh,
	double ctLow,
	double ctHigh,
	TFile *referenceMassFile,
	const TString &figDir,
	const TString &rootDir,
	const TString &tag,
	bool drawMassSlices)
{
	SliceResult out;
	out.ctLow = ctLow;
	out.ctHigh = ctHigh;

	TString sliceCut = Form(
		"(mass > 2.6 && mass < 3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
		"(ctau3D >= %.12g && ctau3D < %.12g)",
		ptLow, ptHigh, yLow, yHigh, ctLow, ctHigh);
	auto dataSlice = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData.reduce(sliceCut)));
	if (!dataSlice || dataSlice->numEntries() <= 0)
		return out;
	out.entries = dataSlice->numEntries();

	RooRealVar obs_mass("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
	obs_mass.SetTitle("mass");
	obs_mass.setBins(std::max(80, massIn.getBins()));
	obs_mass.setRange("massFit", 2.6, 3.5);

	const bool hasReference = referenceMassFile && !referenceMassFile->IsZombie();
	int nSignalGaussComponents = std::clamp(read_int(referenceMassFile, "nSignalGaussComponents", 1), 0, 2);
	int nSignalCBComponents = std::clamp(read_int(referenceMassFile, "nSignalCBComponents", 2), 0, 2);
	int nBkgChebyOrder = std::clamp(read_int(referenceMassFile, "nBkgChebyOrder", yLow >= 1.6f ? 3 : 2), 1, 6);
	if (nSignalGaussComponents + nSignalCBComponents <= 0)
		nSignalGaussComponents = 1;

	const double sigmaFloor = 1e-4;
	const double signalMassMeanRef = finite_or(read_double(referenceMassFile, "signal_mass_mean", 3.096), 3.096);
	const double signalMassSigmaRef = std::max(finite_or(read_double(referenceMassFile, "signal_mass_sigma", 0.035), 0.035), sigmaFloor);
	const double signalMassSigma2Ref = std::max(finite_or(read_double(referenceMassFile, "signal_mass_sigma2", 0.065), 0.065), sigmaFloor);
	const double signalMassCbSigmaRef = std::max(finite_or(read_double(referenceMassFile, "signal_mass_cb_sigma", 0.045), 0.045), sigmaFloor);
	const double signalMassCbSigma2Ref = std::max(finite_or(read_double(referenceMassFile, "signal_mass_cb_sigma2", 0.080), 0.080), sigmaFloor);
	const double signalWindowHalfWidth = std::clamp(2.5 * std::max(signalMassSigmaRef, signalMassCbSigmaRef), 0.050, 0.120);
	const double signalWindowLow = std::max(2.6, signalMassMeanRef - signalWindowHalfWidth);
	const double signalWindowHigh = std::min(3.5, signalMassMeanRef + signalWindowHalfWidth);
	out.signalWindowEntries = static_cast<int>(std::lround(dataSlice->sumEntries(
		Form("mass >= %.12g && mass <= %.12g", signalWindowLow, signalWindowHigh))));
	RooRealVar signal_mass_mean("signal_mass_mean", "signal_mass_mean",
		std::clamp(signalMassMeanRef, 3.05, 3.15), 3.05, 3.15);
	RooRealVar signal_mass_sigma("signal_mass_sigma", "signal_mass_sigma",
		std::clamp(signalMassSigmaRef, sigmaFloor, 0.15), sigmaFloor, 0.15);
	RooRealVar signal_mass_sigma2("signal_mass_sigma2", "signal_mass_sigma2",
		std::clamp(signalMassSigma2Ref, sigmaFloor, 0.25), sigmaFloor, 0.25);
	RooRealVar signal_mass_cb_sigma("signal_mass_cb_sigma", "signal_mass_cb_sigma",
		std::clamp(signalMassCbSigmaRef, sigmaFloor, 0.20), sigmaFloor, 0.20);
	RooRealVar signal_mass_cb_sigma2("signal_mass_cb_sigma2", "signal_mass_cb_sigma2",
		std::clamp(signalMassCbSigma2Ref, sigmaFloor, 0.30), sigmaFloor, 0.30);
	RooRealVar signal_mass_cb_alpha("signal_mass_cb_alpha", "signal_mass_cb_alpha",
		read_double(referenceMassFile, "signal_mass_cb_alpha", 1.5), 0.001, 30.0);
	RooRealVar signal_mass_cb_alpha2("signal_mass_cb_alpha2", "signal_mass_cb_alpha2",
		read_double(referenceMassFile, "signal_mass_cb_alpha2", 2.0), 0.001, 30.0);
	RooRealVar signal_mass_cb_n("signal_mass_cb_n", "signal_mass_cb_n",
		read_double(referenceMassFile, "signal_mass_cb_n", 3.0), 0.001, 50.0);
	RooRealVar signal_mass_cb_n2("signal_mass_cb_n2", "signal_mass_cb_n2",
		read_double(referenceMassFile, "signal_mass_cb_n2", 4.0), 0.001, 50.0);
	RooRealVar signal_mass_frac1("signal_mass_frac1", "signal_mass_frac1",
		std::clamp(read_double(referenceMassFile, "signal_mass_frac1", 0.65), 0.001, 0.999), 0.001, 0.999);
	RooRealVar signal_mass_frac2("signal_mass_frac2", "signal_mass_frac2",
		std::clamp(read_double(referenceMassFile, "signal_mass_frac2", 0.30), 0.001, 0.999), 0.001, 0.999);

	if (hasReference)
	{
		signal_mass_cb_alpha.setConstant(true);
		signal_mass_cb_alpha2.setConstant(true);
		signal_mass_cb_n.setConstant(true);
		signal_mass_cb_n2.setConstant(true);
		signal_mass_frac1.setConstant(true);
		signal_mass_frac2.setConstant(true);
	}

	std::vector<std::unique_ptr<RooConstVar>> massConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> massConstraintPdfs;
	RooArgSet massConstraints;
	auto addMassConstraint = [&](const char *baseName, RooRealVar &var, double central, double sigma) {
		if (!hasReference)
			return;
		massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_mean", baseName), "", central));
		massConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_sigma", baseName), "", sigma));
		massConstraintPdfs.push_back(std::make_unique<RooGaussian>(
			Form("%s_constraint", baseName), Form("%s_constraint", baseName),
			var, *massConstraintConsts[massConstraintConsts.size() - 2], *massConstraintConsts.back()));
		massConstraints.add(*massConstraintPdfs.back());
	};
	addMassConstraint("signal_mass_mean", signal_mass_mean, signalMassMeanRef, 0.010);
	addMassConstraint("signal_mass_sigma", signal_mass_sigma, signalMassSigmaRef, std::max(0.005, 0.25 * signalMassSigmaRef));
	addMassConstraint("signal_mass_sigma2", signal_mass_sigma2, signalMassSigma2Ref, std::max(0.008, 0.25 * signalMassSigma2Ref));
	addMassConstraint("signal_mass_cb_sigma", signal_mass_cb_sigma, signalMassCbSigmaRef, std::max(0.006, 0.25 * signalMassCbSigmaRef));
	addMassConstraint("signal_mass_cb_sigma2", signal_mass_cb_sigma2, signalMassCbSigma2Ref, std::max(0.010, 0.25 * signalMassCbSigma2Ref));

	RooGaussian signal_mass_gaus("signal_mass_gaus", "signal_mass_gaus", obs_mass, signal_mass_mean, signal_mass_sigma);
	RooGaussian signal_mass_gaus2("signal_mass_gaus2", "signal_mass_gaus2", obs_mass, signal_mass_mean, signal_mass_sigma2);
	std::unique_ptr<RooAbsPdf> signal_mass_cb;
	if (nSignalCBComponents >= 2)
	{
		signal_mass_cb = std::make_unique<RooCrystalBall>("signal_mass_cb", "signal_mass_cb",
			obs_mass, signal_mass_mean, signal_mass_cb_sigma, signal_mass_cb_sigma2,
			signal_mass_cb_alpha, signal_mass_cb_n, signal_mass_cb_alpha2, signal_mass_cb_n2);
	}
	else if (nSignalCBComponents == 1)
	{
		signal_mass_cb = std::make_unique<RooCBShape>("signal_mass_cb", "signal_mass_cb",
			obs_mass, signal_mass_mean, signal_mass_cb_sigma, signal_mass_cb_alpha, signal_mass_cb_n);
	}

	std::vector<RooAbsPdf *> signalComponents;
	if (nSignalGaussComponents >= 1)
		signalComponents.push_back(&signal_mass_gaus);
	if (signal_mass_cb)
		signalComponents.push_back(signal_mass_cb.get());
	if (nSignalGaussComponents >= 2)
		signalComponents.push_back(&signal_mass_gaus2);

	RooArgList signalPdfList;
	for (auto *pdf : signalComponents)
		signalPdfList.add(*pdf);
	RooArgList signalFracList;
	if (signalComponents.size() >= 2)
		signalFracList.add(signal_mass_frac1);
	if (signalComponents.size() >= 3)
		signalFracList.add(signal_mass_frac2);

	std::unique_ptr<RooAbsPdf> signal_mass_owned;
	RooAbsPdf *signal_mass = signalComponents.front();
	if (signalComponents.size() > 1)
	{
		signal_mass_owned = std::make_unique<RooAddPdf>("signal_mass", "signal_mass",
			signalPdfList, signalFracList, true);
		signal_mass = signal_mass_owned.get();
	}

	RooRealVar bkg_mass_p1("bkg_mass_p1", "bkg_mass_p1", read_double(referenceMassFile, "bkg_mass_p1", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p2("bkg_mass_p2", "bkg_mass_p2", read_double(referenceMassFile, "bkg_mass_p2", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p3("bkg_mass_p3", "bkg_mass_p3", read_double(referenceMassFile, "bkg_mass_p3", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p4("bkg_mass_p4", "bkg_mass_p4", read_double(referenceMassFile, "bkg_mass_p4", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p5("bkg_mass_p5", "bkg_mass_p5", read_double(referenceMassFile, "bkg_mass_p5", 0.0), -1.0, 1.0);
	RooRealVar bkg_mass_p6("bkg_mass_p6", "bkg_mass_p6", read_double(referenceMassFile, "bkg_mass_p6", 0.0), -1.0, 1.0);
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
	RooChebychev bkg_mass_cheby("bkg_mass_cheby", "bkg_mass_cheby", obs_mass, chebyCoeffList);
	RooRealVar bkg_mass_lambda("bkg_mass_lambda", "bkg_mass_lambda", -1.0, -20.0, -1e-4);
	RooExponential bkg_mass_exp("bkg_mass_exp", "bkg_mass_exp", obs_mass, bkg_mass_lambda);

	const double entries = std::max(1, dataSlice->numEntries());
	// Keep subrange mass-slice fits positive-definite; Chebychev can go negative and stall MIGRAD in narrow ctau slices.
	RooAbsPdf *bkg_mass_pdf = static_cast<RooAbsPdf *>(&bkg_mass_exp);
	RooRealVar Nsig("Nsig", "Nsig", std::max(1.0, 0.35 * entries), 0.0, entries * 1.2);
	RooRealVar Nbkg("Nbkg", "Nbkg", std::max(1.0, 0.65 * entries), 0.0, entries * 1.2);
	RooAddPdf mass_pdf("mass_pdf", "mass_pdf", RooArgList(*signal_mass, *bkg_mass_pdf), RooArgList(Nsig, Nbkg));

	auto fitResult = fit_mass_with_retry(mass_pdf, *dataSlice, out.fitAttempt, &massConstraints);
	if (fitResult)
	{
		out.status = fitResult->status();
		out.covQual = fitResult->covQual();
	}
	out.nsig = finite_or(Nsig.getVal(), 0.0);
	out.nsigErr = std::max(0.0, finite_or(Nsig.getError(), std::sqrt(std::max(0.0, out.nsig))));
	out.nbkg = finite_or(Nbkg.getVal(), 0.0);
	out.nbkgErr = std::max(0.0, finite_or(Nbkg.getError(), std::sqrt(std::max(0.0, out.nbkg))));
	out.nsigSignificance = out.nsigErr > 0.0 ? out.nsig / out.nsigErr : 0.0;
	const double peakCount = static_cast<double>(out.signalWindowEntries);
	const double peakCompatibleSignal = 1.8 * (peakCount + std::sqrt(std::max(1.0, peakCount)));
	out.shapeDrivenSignal = (hasReference && out.nsig > 4.0 && out.nsig > std::max(6.0, peakCompatibleSignal)) ? 1 : 0;
	if (out.shapeDrivenSignal)
	{
		std::cout << Form("  [WARN] mass slice signal is shape-driven: Nsig=%.3f +/- %.3f, peak-window entries=%d",
			out.nsig, out.nsigErr, out.signalWindowEntries) << std::endl;
	}

	const TString massRootDir = TString::Format("%s/mass_slices/%s", rootDir.Data(), tag.Data());
	gSystem->mkdir(massRootDir, true);
	TString massRootName = TString::Format("%s/mass_fit_%s_ct%.3f_%.3f.root",
		massRootDir.Data(), tag.Data(), ctLow, ctHigh);
	massRootName.ReplaceAll("-", "m");
	TFile massOut(massRootName, "RECREATE");
	if (!massOut.IsZombie())
	{
		auto hMass = std::make_unique<TH1D>("hMass", "mass data in ctau slice", 80, 2.6, 3.5);
		dataSlice->fillHistogram(hMass.get(), RooArgList(obs_mass));
		hMass->Write();
		if (fitResult)
			fitResult->Write("massFitResult");
		TParameter<double>("ctLow", ctLow).Write();
		TParameter<double>("ctHigh", ctHigh).Write();
		TParameter<int>("entries", out.entries).Write();
		TParameter<double>("Nsig", out.nsig).Write();
		TParameter<double>("NsigErr", out.nsigErr).Write();
		TParameter<double>("NsigSignificance", out.nsigSignificance).Write();
		TParameter<int>("signalWindowEntries", out.signalWindowEntries).Write();
		TParameter<int>("shapeDrivenSignal", out.shapeDrivenSignal).Write();
		TParameter<double>("Nbkg", out.nbkg).Write();
		TParameter<double>("NbkgErr", out.nbkgErr).Write();
		TParameter<int>("fitStatus", out.status).Write();
		TParameter<int>("covQual", out.covQual).Write();
		TParameter<int>("fitAttempt", out.fitAttempt).Write();
		TParameter<double>("signal_mass_mean", signal_mass_mean.getVal()).Write();
		TParameter<double>("signal_mass_sigma", signal_mass_sigma.getVal()).Write();
		TParameter<double>("signal_mass_cb_sigma", signal_mass_cb_sigma.getVal()).Write();
		TParameter<int>("nMassConstraints", massConstraints.getSize()).Write();
		massOut.Write();
	}

	if (drawMassSlices)
	{
		const TString massFigDir = TString::Format("%s/mass_slices/%s", figDir.Data(), tag.Data());
		gSystem->mkdir(massFigDir, true);
		TCanvas c("c_mass_slice", "c_mass_slice", 800, 800);
		TPad padTop("pad_mass_top", "pad_mass_top", 0.0, 0.25, 1.0, 1.0);
		padTop.SetBottomMargin(0.00001);
		padTop.SetTopMargin(0.08);
		padTop.SetTicks(1, 1);
		padTop.Draw();
		padTop.cd();

		auto *frame = obs_mass.frame(Bins(80), Title(""));
		dataSlice->plotOn(frame, Binning(80), Name("dataOS"), MarkerSize(0.8));
		mass_pdf.plotOn(frame, Name("pdfMASS_tot"), LineColor(kBlack), LineWidth(2));
		mass_pdf.plotOn(frame, Components(*signal_mass), Name("sig"), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2));
		mass_pdf.plotOn(frame, Components(*bkg_mass_pdf), Name("bkg"), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2));
		frame->SetFillStyle(4000);
		frame->GetYaxis()->SetTitle("Events");
		frame->GetYaxis()->SetTitleOffset(1.6);

		double yMax = -1.0;
		if (auto *dataHistForRange = dynamic_cast<RooHist *>(frame->getHist("dataOS")))
		{
			double xPoint = 0.0;
			double yPoint = 0.0;
			for (int i = 0; i < dataHistForRange->GetN(); ++i)
			{
				dataHistForRange->GetPoint(i, xPoint, yPoint);
				const double yHighPoint = yPoint + dataHistForRange->GetErrorYhigh(i);
				if (std::isfinite(yHighPoint))
					yMax = std::max(yMax, yHighPoint);
			}
		}
		if (yMax <= 0.0)
			yMax = std::max(1.0, entries);
		frame->GetYaxis()->SetRangeUser(0.0, yMax * 1.65);
		frame->GetXaxis()->SetLabelSize(0);
		frame->GetXaxis()->SetTitleSize(0);
		frame->GetXaxis()->SetTitle("");
		frame->Draw("e");

		auto findObj = [&](const char *n) -> TObject * { return frame ? frame->findObject(n) : nullptr; };
		TLegend leg(0.50, 0.70, 0.74, 0.89);
		leg.SetBorderSize(0);
		leg.SetFillStyle(0);
		leg.SetTextSize(0.03);
		if (auto *o = findObj("dataOS"))
			leg.AddEntry(o, "Data", "lep");
		if (auto *o = findObj("pdfMASS_tot"))
			leg.AddEntry(o, "Fit", "l");
		if (auto *o = findObj("sig"))
			leg.AddEntry(o, "Signal", "l");
		if (auto *o = findObj("bkg"))
			leg.AddEntry(o, "Background", "l");
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
			tx.SetTextAlign(11);
			tx.SetTextSize(0.04);
			tx.SetTextFont(62);
			tx.DrawLatex(0.23, 0.930, "CMS");
			tx.SetTextFont(72);
			tx.DrawLatex(0.32, 0.930, "Internal");
		}
		{
			TLatex tx;
			tx.SetNDC();
			tx.SetTextSize(0.030);
			tx.SetTextFont(42);
			double xtext = 0.19;
			double y0 = 0.865;
			double dy = -0.048;
			int k = 0;
			tx.DrawLatex(xtext, y0 + dy * k++, "J/#psi data");
			if (yLow == 0)
				tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
			else
				tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
			tx.DrawLatex(xtext, y0 + dy * k++, Form("%.3f < ctau3D < %.3f mm", ctLow, ctHigh));
			tx.DrawLatex(xtext, y0 + dy * k++, Form("status/cov = %d/%d", out.status, out.covQual));
		}
		{
			TLatex tp;
			tp.SetNDC();
			tp.SetTextSize(0.024);
			tp.SetTextFont(42);
			double xtext = 0.74;
			double y0 = 0.87;
			double dy = -0.045;
			int k = 0;
			tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{sig} = %.4g #pm %.3g", out.nsig, out.nsigErr));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{bkg} = %.4g #pm %.3g", out.nbkg, out.nbkgErr));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("SR entries = %d%s", out.signalWindowEntries, out.shapeDrivenSignal ? " (veto)" : ""));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("#mu = %.5g #pm %.2g", signal_mass_mean.getVal(), signal_mass_mean.getError()));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("#sigma_{G1} = %.4g #pm %.2g", signal_mass_sigma.getVal(), signal_mass_sigma.getError()));
			if (signal_mass_cb)
				tp.DrawLatex(xtext, y0 + dy * k++, Form("#sigma_{CB1} = %.4g #pm %.2g", signal_mass_cb_sigma.getVal(), signal_mass_cb_sigma.getError()));
		}

		c.cd();
		TPad padPull("pad_mass_pull", "pad_mass_pull", 0.0, 0.0, 1.0, 0.25);
		padPull.SetTopMargin(0.00001);
		padPull.SetBottomMargin(0.4);
		padPull.SetFillStyle(4000);
		padPull.SetFrameFillStyle(4000);
		padPull.SetTicks(1, 1);
		padPull.Draw();
		padPull.cd();

		auto *frameTMP = static_cast<RooPlot *>(frame->Clone("frame_mass_pull_tmp"));
		RooHist *hpull = frameTMP->pullHist("dataOS", "pdfMASS_tot", true);
		if (hpull)
			hpull->SetMarkerSize(0.8);
		auto *pullFrame = obs_mass.frame(Title(""));
		if (hpull)
			pullFrame->addPlotable(hpull, "P");
		pullFrame->GetYaxis()->SetTitle("Pull");
		pullFrame->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
		pullFrame->GetXaxis()->CenterTitle();
		pullFrame->SetMinimum(-8);
		pullFrame->SetMaximum(8);
		pullFrame->GetYaxis()->SetNdivisions(505);
		pullFrame->GetYaxis()->SetTitleSize(0.12);
		pullFrame->GetYaxis()->SetLabelSize(0.10);
		pullFrame->GetXaxis()->SetTitleSize(0.15);
		pullFrame->GetXaxis()->SetLabelSize(0.10);
		pullFrame->Draw();

		TLine line(2.6, 0.0, 3.5, 0.0);
		line.SetLineStyle(2);
		line.Draw("same");

		auto chiM = chi2_from_pull(hpull);
		int npar = fitResult ? fitResult->floatParsFinal().getSize() : 0;
		int ndf = std::max(1, chiM.second - npar);
		{
			TLatex tc;
			tc.SetNDC();
			tc.SetTextSize(0.085);
			tc.SetTextFont(42);
			tc.SetTextAlign(33);
			tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", chiM.first, ndf));
		}

		TString name = TString::Format("%s/mass_slice_%s_ct%.3f_%.3f.pdf",
			massFigDir.Data(), tag.Data(), ctLow, ctHigh);
		name.ReplaceAll("-", "m");
		c.SaveAs(name);
		delete pullFrame;
		delete frameTMP;
		delete frame;
	}

	return out;
}
}

void draw_raw_mass_subrange_distributions(float ptLow = 3.0, float ptHigh = 6.5, float yLow = 1.6, float yHigh = 2.4,
	double ctMin = -20.0, double ctMax = 30.0, double ctStep = 2.0,
	int maxSlices = -1, bool applyQuantileCut = false, bool useLogY = false, bool saveEachSlice = true)
{
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	if (ctStep <= 0.0 || ctMax <= ctMin)
	{
		std::cerr << "ERROR: invalid ctau binning: [" << ctMin << ", " << ctMax << "] step " << ctStep << std::endl;
		return;
	}

	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
	TFile *inputFile = TFile::Open(kDataRoot);
	if (!inputFile || inputFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open input data file: " << kDataRoot << std::endl;
		return;
	}
	auto *inputData = dynamic_cast<RooDataSet *>(inputFile->Get(kDatasetName));
	if (!inputData)
	{
		std::cerr << "ERROR: RooDataSet '" << kDatasetName << "' not found in " << kDataRoot << std::endl;
		return;
	}

	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && (recoQQsign==0) && %s",
		ptLow, ptHigh, yLow, yHigh, kAccCut);
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	auto *ctau3DVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
	auto *ctau3DErrVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DErr"));
	if (!ctau3DVar || !ctau3DErrVar)
	{
		std::cerr << "ERROR: required variables ctau3D/ctau3DErr are missing in dataset." << std::endl;
		return;
	}

	TString finalCut;
	if (applyQuantileCut)
	{
		const auto ctau3DRange = quantile_range(*dataBasic, *ctau3DVar);
		const auto ctau3DErrRange = quantile_range(*dataBasic, *ctau3DErrVar);
		finalCut = Form(
			"%s && !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes) && "
			"(ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)",
			cutBasic.Data(), ctau3DRange.first, ctau3DRange.second, ctau3DErrRange.first, ctau3DErrRange.second);
	}
	else
	{
		finalCut = Form("%s && !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes)", cutBasic.Data());
	}
	auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(finalCut)));
	if (!dataSel || dataSel->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after final selection: " << finalCut << std::endl;
		return;
	}
	std::cout << Form("[INFO] selected pp data entries for raw mass diagnostic: basic=%d final=%d",
		dataBasic->numEntries(), dataSel->numEntries()) << std::endl;

	auto edges = build_strong_adaptive_ctau_edges(ctMin, ctMax, ctStep);
	if (edges.size() < 2)
	{
		std::cerr << "ERROR: no adaptive ctau edges were built." << std::endl;
		return;
	}
	if (maxSlices > 0 && maxSlices + 1 < static_cast<int>(edges.size()))
		edges.resize(maxSlices + 1);

	double minEdgeWidth = std::numeric_limits<double>::infinity();
	double maxEdgeWidth = 0.0;
	double maxNeighborWidthRatio = 1.0;
	for (size_t i = 1; i < edges.size(); ++i)
	{
		const double width = edges[i] - edges[i - 1];
		minEdgeWidth = std::min(minEdgeWidth, width);
		maxEdgeWidth = std::max(maxEdgeWidth, width);
		if (i > 1)
		{
			const double prevWidth = edges[i - 1] - edges[i - 2];
			if (prevWidth > 0.0 && width > 0.0)
				maxNeighborWidthRatio = std::max(maxNeighborWidthRatio, std::max(width / prevWidth, prevWidth / width));
		}
	}
	std::cout << Form("[INFO] strong adaptive raw-mass ctau slicing: total %zu bins, width range [%.6g, %.6g], max adjacent width ratio %.3f",
		edges.size() - 1, minEdgeWidth, maxEdgeWidth, maxNeighborWidthRatio) << std::endl;

	const TString yTag = y_tag(yLow, yHigh);
	const TString tag = bin_tag(ptLow, ptHigh, yLow, yHigh);
	TString ctRangeTag = TString::Format("_ct%.1f_%.1f", ctMin, ctMax);
	ctRangeTag.ReplaceAll("-", "m");
	ctRangeTag.ReplaceAll(".", "p");
	const TString outputTag = tag + "_rawMass_strongAdaptive" + ctRangeTag + (applyQuantileCut ? "" : "_noQuantile");
	const TString figDir = TString::Format("%s/figs/%s/subrange_mass_raw", kReferenceDir, yTag.Data());
	const TString rootDir = TString::Format("%s/roots/%s/subrange_mass_raw", kReferenceDir, yTag.Data());
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(rootDir, true);
	const TString outputName = TString::Format("%s/raw_mass_subranges_%s.root", rootDir.Data(), outputTag.Data());
	const TString multiPdfName = TString::Format("%s/raw_mass_subranges_%s.pdf", figDir.Data(), outputTag.Data());

	TFile outFile(outputName, "RECREATE");
	if (outFile.IsZombie())
	{
		std::cerr << "ERROR: cannot create output ROOT file: " << outputName << std::endl;
		return;
	}

	double ctLowOut = 0.0;
	double ctHighOut = 0.0;
	int entriesOut = 0;
	int massBinsOut = 0;
	int jpsiWindowEntriesOut = 0;
	TTree summary("rawMassSliceSummary", "raw mass subrange diagnostic summary");
	summary.Branch("ctLow", &ctLowOut, "ctLow/D");
	summary.Branch("ctHigh", &ctHighOut, "ctHigh/D");
	summary.Branch("entries", &entriesOut, "entries/I");
	summary.Branch("massBins", &massBinsOut, "massBins/I");
	summary.Branch("jpsiWindowEntries", &jpsiWindowEntriesOut, "jpsiWindowEntries/I");

	TCanvas c("c_raw_mass_slice", "c_raw_mass_slice", 800, 800);
	c.Print(multiPdfName + "[");
	RooRealVar obs_mass("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
	obs_mass.SetTitle("mass");

	int drawnSlices = 0;
	for (size_t i = 0; i + 1 < edges.size(); ++i)
	{
		const double ctLow = edges[i];
		const double ctHigh = edges[i + 1];
		TString sliceCut = Form("ctau3D >= %.12g && ctau3D < %.12g", ctLow, ctHigh);
		auto dataSlice = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(dataSel->reduce(sliceCut)));
		const int entries = dataSlice ? dataSlice->numEntries() : 0;
		const int massBins = adaptive_mass_bins(entries);
		const int jpsiWindowEntries = dataSlice ? static_cast<int>(std::lround(dataSlice->sumEntries("mass >= 2.976 && mass <= 3.216"))) : 0;

		ctLowOut = ctLow;
		ctHighOut = ctHigh;
		entriesOut = entries;
		massBinsOut = entries > 0 ? massBins : 0;
		jpsiWindowEntriesOut = jpsiWindowEntries;
		summary.Fill();

		std::cout << Form("[RAW SLICE %zu/%zu] ctau %.5g to %.5g: entries=%d, massBins=%d, J/psi-window=%d",
			i + 1, edges.size() - 1, ctLow, ctHigh, entries, massBinsOut, jpsiWindowEntries) << std::endl;
		if (!dataSlice || entries <= 0)
			continue;

		auto hMass = std::make_unique<TH1D>(Form("hMass_ct_%04zu", i), "raw mass in ctau slice", massBins, 2.6, 3.5);
		hMass->Sumw2();
		hMass->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
		hMass->GetYaxis()->SetTitle("Events");
		dataSlice->fillHistogram(hMass.get(), RooArgList(obs_mass));
		hMass->Write();

		c.Clear();
		c.SetTicks(1, 1);
		c.SetLogy(useLogY);
		hMass->SetMarkerStyle(20);
		hMass->SetMarkerSize(0.75);
		hMass->SetLineColor(kBlack);
		hMass->SetMarkerColor(kBlack);
		const double maxContent = std::max(1.0, hMass->GetMaximum());
		const double yAxisMax = useLogY ? maxContent * 24.0 : maxContent * 1.65;
		const double jpsiLineMax = useLogY ? maxContent * 8.0 : maxContent * 1.05;
		if (useLogY)
		{
			hMass->SetMinimum(0.5);
			hMass->SetMaximum(yAxisMax);
		}
		else
		{
			hMass->SetMinimum(0.0);
			hMass->SetMaximum(yAxisMax);
		}
		hMass->Draw("E");

		TLine jpsiLine(3.096, useLogY ? 0.5 : 0.0, 3.096, jpsiLineMax);
		jpsiLine.SetLineStyle(2);
		jpsiLine.SetLineColor(kBlue + 1);
		jpsiLine.Draw("same");

		TLatex tx;
		tx.SetNDC();
		tx.SetTextFont(62);
		tx.SetTextSize(0.034);
		tx.DrawLatex(0.20, 0.925, "CMS");
		tx.SetTextFont(72);
		tx.DrawLatex(0.285, 0.925, "Internal");
		tx.SetTextFont(42);
		tx.SetTextSize(0.030);
		tx.DrawLatex(0.16, 0.86, "Raw mass distribution, no mass fit");
		if (yLow == 0.0f)
			tx.DrawLatex(0.16, 0.81, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
		else
			tx.DrawLatex(0.16, 0.81, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
		tx.DrawLatex(0.16, 0.76, Form("%.5g < ctau3D < %.5g mm", ctLow, ctHigh));
		tx.DrawLatex(0.16, 0.71, Form("Entries = %d, bins = %d", entries, massBins));
		tx.DrawLatex(0.16, 0.66, Form("2.976 < m < 3.216: %d", jpsiWindowEntries));

		c.Print(multiPdfName);
		if (saveEachSlice)
		{
			TString slicePdfName = TString::Format("%s/raw_mass_%s_ct%.5f_%.5f.pdf", figDir.Data(), outputTag.Data(), ctLow, ctHigh);
			slicePdfName.ReplaceAll("-", "m");
			c.SaveAs(slicePdfName);
		}
		++drawnSlices;
	}
	c.Print(multiPdfName + "]");

	summary.Write();
	TParameter<double>("ptLow", ptLow).Write();
	TParameter<double>("ptHigh", ptHigh).Write();
	TParameter<double>("yLow", yLow).Write();
	TParameter<double>("yHigh", yHigh).Write();
	TParameter<double>("ctMin", ctMin).Write();
	TParameter<double>("ctMax", ctMax).Write();
	TParameter<double>("ctStep", ctStep).Write();
	TParameter<int>("nRawMassSlices", static_cast<int>(edges.size()) - 1).Write();
	TParameter<int>("nDrawnRawMassSlices", drawnSlices).Write();
	TParameter<int>("strongAdaptiveCtauBinning", 1).Write();
	TParameter<double>("adaptiveMinBinWidth", minEdgeWidth).Write();
	TParameter<double>("adaptiveMaxBinWidth", maxEdgeWidth).Write();
	TParameter<double>("adaptiveMaxNeighborWidthRatio", maxNeighborWidthRatio).Write();
	TParameter<int>("applyQuantileCut", applyQuantileCut ? 1 : 0).Write();
	TParameter<int>("useLogY", useLogY ? 1 : 0).Write();
	outFile.Write();

	std::cout << "[DONE] wrote " << outputName << std::endl;
	std::cout << "[DONE] wrote " << multiPdfName << std::endl;
}

void subrange_method(float ptLow = 3.0, float ptHigh = 6.5, float yLow = 1.6, float yHigh = 2.4,
	double ctMin = -8.0, double ctMax = 10.0, double ctStep = 2.0,
	int maxSlices = -1, bool drawMassSlices = true, bool trimSparseFwd = true,
	bool applyQuantileCut = true, bool reuseMassSlices = false, double centerVetoHalfWidth = 0.0)
{
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	if (ctStep <= 0.0 || ctMax <= ctMin)
	{
		std::cerr << "ERROR: invalid ctau binning: [" << ctMin << ", " << ctMax << "] step " << ctStep << std::endl;
		return;
	}

	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
	TFile *inputFile = TFile::Open(kDataRoot);
	if (!inputFile || inputFile->IsZombie())
	{
		std::cerr << "ERROR: cannot open input data file: " << kDataRoot << std::endl;
		return;
	}
	auto *inputData = dynamic_cast<RooDataSet *>(inputFile->Get(kDatasetName));
	if (!inputData)
	{
		std::cerr << "ERROR: RooDataSet '" << kDatasetName << "' not found in " << kDataRoot << std::endl;
		return;
	}
	TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && (recoQQsign==0) && %s",
		ptLow, ptHigh, yLow, yHigh, kAccCut);
	auto dataBasic = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(cutBasic)));
	if (!dataBasic || dataBasic->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after basic selection: " << cutBasic << std::endl;
		return;
	}

	auto *massVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("mass"));
	auto *ctau3DVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3D"));
	auto *ctau3DErrVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DErr"));
	auto *ctau3DResVar = dynamic_cast<RooRealVar *>(dataBasic->get()->find("ctau3DRes"));
	if (!massVar || !ctau3DVar || !ctau3DErrVar || !ctau3DResVar)
	{
		std::cerr << "ERROR: required variables mass/ctau3D/ctau3DErr/ctau3DRes are missing in dataset." << std::endl;
		return;
	}

	const auto ctau3DRange = quantile_range(*dataBasic, *ctau3DVar);
	const auto ctau3DErrRange = quantile_range(*dataBasic, *ctau3DErrVar);
	TString finalCut;
	if (applyQuantileCut)
	{
		finalCut = Form(
			"%s && !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes) && "
			"(ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)",
			cutBasic.Data(), ctau3DRange.first, ctau3DRange.second, ctau3DErrRange.first, ctau3DErrRange.second);
	}
	else
	{
		finalCut = Form("%s && !TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes)", cutBasic.Data());
		std::cout << Form("[INFO] quantile cleanup disabled: nominal ctau range would be [%.6g, %.6g] and ctauErr [%.6g, %.6g]",
			ctau3DRange.first, ctau3DRange.second, ctau3DErrRange.first, ctau3DErrRange.second) << std::endl;
	}
	auto dataSel = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData->reduce(finalCut)));
	if (!dataSel || dataSel->numEntries() <= 0)
	{
		std::cerr << "ERROR: no entries after final selection: " << finalCut << std::endl;
		return;
	}
	massVar = dynamic_cast<RooRealVar *>(dataSel->get()->find("mass"));
	std::cout << Form("[INFO] selected pp data entries: basic=%d final=%d", dataBasic->numEntries(), dataSel->numEntries()) << std::endl;

	const TString yTag = y_tag(yLow, yHigh);
	const TString tag = bin_tag(ptLow, ptHigh, yLow, yHigh);
	TString ctStepSuffix = "";
	if (std::abs(ctStep - 1.0) > 1e-9)
	{
		ctStepSuffix = "_ctStep" + format_tag(ctStep);
		ctStepSuffix.ReplaceAll(".", "p");
	}
	const TString outputTag = tag + (trimSparseFwd ? "" : "_fullrange") + (applyQuantileCut ? "" : "_noQuantile") + ctStepSuffix;
	const TString figDir = TString::Format("figs/%s/subrange", yTag.Data());
	const TString rootDir = TString::Format("roots/%s/subrange", yTag.Data());
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(rootDir, true);
	const TString outputName = TString::Format("%s/subrange_result_%s.root", rootDir.Data(), outputTag.Data());

	const TString referenceMassName = TString::Format("%s/roots/%s/mass/mass_model_%s.root",
		kReferenceDir, yTag.Data(), tag.Data());
	std::unique_ptr<TFile> referenceMass(TFile::Open(referenceMassName, "READ"));
	if (!referenceMass || referenceMass->IsZombie())
	{
		std::cout << "[INFO] reference mass fit not found, using floating default signal shape: "
				  << referenceMassName << std::endl;
		referenceMass.reset();
	}
	else
	{
		std::cout << "[INFO] using reference mass shape: " << referenceMassName << std::endl;
	}

	std::vector<double> edges;
	auto addEdge = [&](double x) {
		if (edges.empty() || std::abs(edges.back() - x) > 1e-9)
			edges.push_back(x);
	};
	auto appendFixedBins = [&](double start, double stop, int nBins) {
		if (stop <= start || nBins <= 0)
			return;
		addEdge(start);
		for (int i = 1; i < nBins; ++i)
			addEdge(start + (stop - start) * static_cast<double>(i) / static_cast<double>(nBins));
		addEdge(stop);
	};

	double fitCtMin = ctMin;
	double fitCtMax = ctMax;
	const bool isForward = yLow >= 1.6f;
	const bool isLowPtForward = isForward && ptHigh <= 2.0f;
	const bool useAdaptiveFwdTrim = trimSparseFwd && isLowPtForward;
	if (useAdaptiveFwdTrim && ctMax - ctMin > 20.0)
	{
		fitCtMin = std::max(ctMin, -2.0);
		fitCtMax = std::min(ctMax, 10.0);
		std::cout << Form("[INFO] fwd low-pT tail signal is sparse: reducing ctau range from [%.3f, %.3f] to [%.3f, %.3f]",
			ctMin, ctMax, fitCtMin, fitCtMax) << std::endl;
	}

	auto appendBand = [&](double start, double stop, double targetWidth) {
		start = std::max(start, fitCtMin);
		stop = std::min(stop, fitCtMax);
		if (stop <= start || targetWidth <= 0.0)
			return;
		const int nBins = std::max(1, static_cast<int>(std::ceil((stop - start) / targetWidth - 1e-9)));
		appendFixedBins(start, stop, nBins);
	};

	const double outerStep = std::max(0.005, std::min(ctStep, 0.50));
	appendBand(fitCtMin, -3.0, outerStep);
	appendBand(-3.0, -2.5, 0.250);
	appendBand(-2.5, -2.0, 0.125);
	appendBand(-2.0, -1.45, 0.090);
	appendBand(-1.45, -1.05, 0.060);
	appendBand(-1.05, -0.75, 0.040);
	appendBand(-0.75, -0.50, 0.025);
	appendBand(-0.50, -0.30, 0.015);
	appendBand(-0.30, -0.15, 0.010);
	appendBand(-0.15, -0.06, 0.0075);
	appendBand(-0.06, 0.06, 0.0050);
	appendBand(0.06, 0.15, 0.0075);
	appendBand(0.15, 0.30, 0.010);
	appendBand(0.30, 0.50, 0.015);
	appendBand(0.50, 0.75, 0.025);
	appendBand(0.75, 1.05, 0.040);
	appendBand(1.05, 1.45, 0.060);
	appendBand(1.45, 2.00, 0.090);
	appendBand(2.00, 2.50, 0.125);
	appendBand(2.50, 3.00, 0.250);
	appendBand(3.00, fitCtMax, outerStep);

	double minEdgeWidth = std::numeric_limits<double>::infinity();
	double maxEdgeWidth = 0.0;
	double maxNeighborWidthRatio = 1.0;
	for (size_t i = 1; i < edges.size(); ++i)
	{
		const double width = edges[i] - edges[i - 1];
		minEdgeWidth = std::min(minEdgeWidth, width);
		maxEdgeWidth = std::max(maxEdgeWidth, width);
		if (i > 1)
		{
			const double prevWidth = edges[i - 1] - edges[i - 2];
			if (prevWidth > 0.0 && width > 0.0)
				maxNeighborWidthRatio = std::max(maxNeighborWidthRatio, std::max(width / prevWidth, prevWidth / width));
		}
	}
	std::cout << Form("[INFO] adaptive ctau slicing: total %zu bins, width range [%.6g, %.6g], max adjacent width ratio %.3f; RooBinning + fitTo IntegrateBins enabled",
		edges.size() > 0 ? edges.size() - 1 : 0, minEdgeWidth, maxEdgeWidth, maxNeighborWidthRatio) << std::endl;
	if (maxSlices > 0 && maxSlices + 1 < static_cast<int>(edges.size()))
		edges.resize(maxSlices + 1);

	std::vector<SliceResult> sliceResults;
	if (reuseMassSlices)
	{
		std::unique_ptr<TFile> prevOut(TFile::Open(outputName, "READ"));
		if (!prevOut || prevOut->IsZombie())
		{
			std::cerr << "ERROR: cannot reuse mass slices; previous output not found: " << outputName << std::endl;
			return;
		}
		auto *summary = dynamic_cast<TTree *>(prevOut->Get("sliceSummary"));
		if (!summary)
		{
			std::cerr << "ERROR: sliceSummary not found in previous output: " << outputName << std::endl;
			return;
		}
		SliceResult r;
		summary->SetBranchAddress("ctLow", &r.ctLow);
		summary->SetBranchAddress("ctHigh", &r.ctHigh);
		summary->SetBranchAddress("entries", &r.entries);
		summary->SetBranchAddress("Nsig", &r.nsig);
		summary->SetBranchAddress("NsigErr", &r.nsigErr);
		summary->SetBranchAddress("Nbkg", &r.nbkg);
		summary->SetBranchAddress("NbkgErr", &r.nbkgErr);
		summary->SetBranchAddress("NsigSignificance", &r.nsigSignificance);
		if (summary->GetBranch("signalWindowEntries"))
			summary->SetBranchAddress("signalWindowEntries", &r.signalWindowEntries);
		if (summary->GetBranch("shapeDrivenSignal"))
			summary->SetBranchAddress("shapeDrivenSignal", &r.shapeDrivenSignal);
		if (summary->GetBranch("usedInCtauFit"))
			summary->SetBranchAddress("usedInCtauFit", &r.usedInCtauFit);
		summary->SetBranchAddress("fitStatus", &r.status);
		summary->SetBranchAddress("covQual", &r.covQual);
		const Long64_t nRows = summary->GetEntries();
		sliceResults.reserve(nRows);
		for (Long64_t i = 0; i < nRows; ++i)
		{
			r = SliceResult();
			summary->GetEntry(i);
			sliceResults.push_back(r);
		}
		std::cout << Form("[INFO] reusing %zu mass subrange slices from %s; skipping mass fits.",
			sliceResults.size(), outputName.Data()) << std::endl;
	}
	else
	{
		sliceResults.reserve(edges.size() > 0 ? edges.size() - 1 : 0);
		for (size_t i = 0; i + 1 < edges.size(); ++i)
		{
			std::cout << Form("[SLICE %zu/%zu] ctau %.3f to %.3f", i + 1, edges.size() - 1, edges[i], edges[i + 1]) << std::endl;
			sliceResults.push_back(fit_mass_slice(*dataSel, *massVar, ptLow, ptHigh, yLow, yHigh,
				edges[i], edges[i + 1], referenceMass.get(), figDir, rootDir, outputTag, drawMassSlices));
			const auto &r = sliceResults.back();
			std::cout << Form("  entries=%d peakWin=%d Nsig=%.3f +/- %.3f status=%d covQual=%d%s",
				r.entries, r.signalWindowEntries, r.nsig, r.nsigErr, r.status, r.covQual,
				r.shapeDrivenSignal ? " shape-veto" : "") << std::endl;
		}
	}

	if (sliceResults.empty())
	{
		std::cerr << "ERROR: no ctau slices were fitted." << std::endl;
		return;
	}

	double maxSliceSignal = 0.0;
	for (const auto &r : sliceResults)
		maxSliceSignal = std::max(maxSliceSignal, r.nsig);
	const double edgeSignalThreshold = std::max(3.0, 0.010 * maxSliceSignal);
	std::cout << Form("[INFO] mass-signal edge threshold for ctau fit: Nsig >= %.3f (1%% of max slice %.3f)",
		edgeSignalThreshold, maxSliceSignal) << std::endl;
	auto hasMassSignal = [edgeSignalThreshold](const SliceResult &r) {
		const double relErr = r.nsig > 0.0 ? r.nsigErr / r.nsig : 999.0;
		return r.entries >= 8 && r.status == 0 && r.covQual >= 2 &&
			r.nsig > 0.0 && r.nsigErr > 0.0 && r.nsigSignificance >= 2.5 &&
			r.nsig >= edgeSignalThreshold && relErr <= 0.75 && r.shapeDrivenSignal == 0;
	};
	int firstSignalSlice = -1;
	int lastSignalSlice = -1;
	for (size_t i = 0; i < sliceResults.size(); ++i)
	{
		if (!hasMassSignal(sliceResults[i]))
			continue;
		if (firstSignalSlice < 0)
			firstSignalSlice = static_cast<int>(i);
		lastSignalSlice = static_cast<int>(i);
	}
	std::vector<SliceResult> ctauFitSlices;
	if (firstSignalSlice >= 0 && lastSignalSlice >= firstSignalSlice)
	{
		const int originalFirstSignalSlice = firstSignalSlice;
		while (firstSignalSlice <= lastSignalSlice && sliceResults[firstSignalSlice].ctHigh <= kCtauFitMinOverride + 1e-9)
			++firstSignalSlice;
		if (firstSignalSlice != originalFirstSignalSlice)
		{
			std::cout << Form("[INFO] forcing ctau fit lower edge to %.3f mm: first slice %d -> %d, range starts at %.3f",
				kCtauFitMinOverride, originalFirstSignalSlice, firstSignalSlice,
				(firstSignalSlice <= lastSignalSlice ? sliceResults[firstSignalSlice].ctLow : 999.0)) << std::endl;
		}
		for (size_t i = 0; i < sliceResults.size(); ++i)
			sliceResults[i].usedInCtauFit = (static_cast<int>(i) >= firstSignalSlice && static_cast<int>(i) <= lastSignalSlice) ? 1 : 0;
		ctauFitSlices.assign(sliceResults.begin() + firstSignalSlice, sliceResults.begin() + lastSignalSlice + 1);
		if (ctauFitSlices.size() != sliceResults.size())
		{
			std::cout << Form("[INFO] excluding no-mass-signal edge slices from ctau fit: kept %zu/%zu, range [%.3f, %.3f]",
				ctauFitSlices.size(), sliceResults.size(), ctauFitSlices.front().ctLow, ctauFitSlices.back().ctHigh) << std::endl;
		}
	}
	else
	{
		std::cout << "[INFO] no significant mass-signal slice found; keeping all slices for diagnostic ctau fit." << std::endl;
		ctauFitSlices = sliceResults;
	}

	std::vector<double> histEdges;
	histEdges.reserve(ctauFitSlices.size() + 1);
	for (const auto &r : ctauFitSlices)
		histEdges.push_back(r.ctLow);
	histEdges.push_back(ctauFitSlices.back().ctHigh);
	auto hSignalYield = std::make_unique<TH1D>("hSignalYield", "signal yield from mass subranges",
		static_cast<int>(ctauFitSlices.size()), histEdges.data());
	hSignalYield->Sumw2();
	hSignalYield->GetXaxis()->SetTitle("ctau3D (mm)");
	hSignalYield->GetYaxis()->SetTitle("Signal yield from mass fit");

	for (size_t i = 0; i < ctauFitSlices.size(); ++i)
	{
		const auto &r = ctauFitSlices[i];
		hSignalYield->SetBinContent(i + 1, std::max(0.0, r.nsig));
		const double err = r.nsigErr > 0.0 ? r.nsigErr : std::sqrt(std::max(0.0, r.nsig));
		hSignalYield->SetBinError(i + 1, err);
	}

	RooRealVar ctau("ctau3D", "ctau3D", histEdges.front(), histEdges.back(), "mm");
	ctau.setBins(static_cast<int>(ctauFitSlices.size()));
	RooBinning ctauBinning(static_cast<int>(histEdges.size()) - 1, histEdges.data(), "ctauSubrangeBinning");
	ctau.setBinning(ctauBinning);
	RooDataHist signalYieldData("signalYieldData", "signalYieldData", RooArgList(ctau), hSignalYield.get());

	const double centerVeto = std::max(0.0, centerVetoHalfWidth);
	const bool useCenterVeto = centerVeto > 0.0 && histEdges.front() < -centerVeto && histEdges.back() > centerVeto;
	TString ctauFitRangeName = "";
	if (useCenterVeto)
	{
		ctau.setRange("ctauFitLeft", histEdges.front(), -centerVeto);
		ctau.setRange("ctauFitRight", centerVeto, histEdges.back());
		ctauFitRangeName = "ctauFitLeft,ctauFitRight";
		std::cout << Form("[INFO] excluding central ctau window from final fit and pull: |ctau| < %.5g mm", centerVeto) << std::endl;
	}

	const double ctRangeWidth = histEdges.back() - histEdges.front();
	double minBinWidth = ctRangeWidth;
	for (size_t i = 1; i < histEdges.size(); ++i)
		minBinWidth = std::min(minBinWidth, histEdges[i] - histEdges[i - 1]);
	const double promptSigmaMax = std::min(1.0, std::max(0.10, 0.50 * ctRangeWidth));
	const double promptSigmaSeed = std::clamp(std::max(0.15, 0.5 * minBinWidth), 0.03, promptSigmaMax);

	const TString referencePrName = TString::Format("%s/roots/%s/ctau_pr/ctau_resolution_%s.root",
		kReferenceDir, yTag.Data(), tag.Data());
	const TString referenceNpName = TString::Format("%s/roots/%s/ctau_np/ctau_np_model_%s.root",
		kReferenceDir, yTag.Data(), tag.Data());
	std::unique_ptr<TFile> referencePr(TFile::Open(referencePrName, "READ"));
	std::unique_ptr<TFile> referenceNp(TFile::Open(referenceNpName, "READ"));
	const bool useReferenceCtauMC = referencePr && !referencePr->IsZombie() && referenceNp && !referenceNp->IsZombie();
	if (useReferenceCtauMC)
	{
		std::cout << "[INFO] using MC-constrained ctau PR model: " << referencePrName << std::endl;
		std::cout << "[INFO] using MC-constrained ctau NP model: " << referenceNpName << std::endl;
	}
	else
	{
		std::cout << "[INFO] reference ctau MC files not found; using fallback Gaussian + single-sided decay model." << std::endl;
	}

	auto fitRangeAround = [](double center, double hardLow, double hardHigh, double relWidth, double absWidth) {
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

	std::vector<std::unique_ptr<RooConstVar>> ctauConstraintConsts;
	std::vector<std::unique_ptr<RooGaussian>> ctauConstraintPdfs;
	RooArgSet ctauConstraints;
	auto addCtauConstraint = [&](const char *baseName, RooRealVar &var, double central, double sigma) {
		if (kDisableCtauConstraints)
			return;
		sigma *= kCtauConstraintScale;
		ctauConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_mean", baseName), "", central));
		ctauConstraintConsts.push_back(std::make_unique<RooConstVar>(Form("%s_sigma", baseName), "", sigma));
		ctauConstraintPdfs.push_back(std::make_unique<RooGaussian>(
			Form("%s_constraint", baseName), Form("%s_constraint", baseName),
			var, *ctauConstraintConsts[ctauConstraintConsts.size() - 2], *ctauConstraintConsts.back()));
		ctauConstraints.add(*ctauConstraintPdfs.back());
	};

	RooRealVar fallbackPromptMean("fallbackPromptMean", "fallbackPromptMean", 0.0, -0.30, 0.30);
	RooRealVar fallbackPromptSigma("fallbackPromptSigma", "fallbackPromptSigma", promptSigmaSeed, 0.02, promptSigmaMax);
	RooGaussian fallbackPromptPdf("fallbackPromptPdf", "fallbackPromptPdf", ctau, fallbackPromptMean, fallbackPromptSigma);
	RooGaussModel fallbackPromptResolution("fallbackPromptResolution", "fallbackPromptResolution", ctau, fallbackPromptMean, fallbackPromptSigma);
	RooRealVar fallbackNonpromptTau("fallbackNonpromptTau", "fallbackNonpromptTau", 0.7, 0.05, std::max(10.0, ctRangeWidth));
	RooDecay fallbackNonpromptPdf("fallbackNonpromptPdf", "fallbackNonpromptPdf", ctau, fallbackNonpromptTau, fallbackPromptResolution, RooDecay::SingleSided);

	const int nResolutionComponents = useReferenceCtauMC ? std::clamp(read_int(referencePr.get(), "nResolutionComponents", 1), 1, 4) : 1;
	const double ctauTime1MeanVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime1Mean", 0.0) : 0.0;
	const double ctauTime1ScaleVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime1Scale", promptSigmaSeed) : promptSigmaSeed;
	const auto ctauTime1ScaleRange = fitRangeAround(ctauTime1ScaleVal, 0.005, 0.20, 1.0, 0.015);
	RooRealVar ctauTime1Mean("ctauTime1Mean", "ctauTime1Mean", ctauTime1MeanVal, ctauTime1MeanVal - 0.05, ctauTime1MeanVal + 0.05);
	RooRealVar ctauTime1Scale("ctauTime1Scale", "ctauTime1Scale",
		std::clamp(ctauTime1ScaleVal, ctauTime1ScaleRange.first, ctauTime1ScaleRange.second),
		ctauTime1ScaleRange.first, ctauTime1ScaleRange.second);
	RooGaussModel ctauTime1("ctauTime1", "ctauTime1", ctau, ctauTime1Mean, ctauTime1Scale);

	const double ctauTime2MeanVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime2Mean", 0.0) : 0.0;
	const double ctauTime2ScaleVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime2Scale", 2.0 * ctauTime1ScaleVal) : 2.0 * ctauTime1ScaleVal;
	const double ctauTime2DeltaVal = std::max(0.001, ctauTime2ScaleVal - ctauTime1ScaleVal);
	const auto ctauTime2DeltaRange = fitRangeAround(ctauTime2DeltaVal, 0.001, 0.30, 2.0, 0.040);
	RooRealVar ctauTime2Mean("ctauTime2Mean", "ctauTime2Mean", ctauTime2MeanVal, ctauTime2MeanVal - 0.05, ctauTime2MeanVal + 0.05);
	RooRealVar ctauTime2Delta("ctauTime2Delta", "ctauTime2Delta",
		std::clamp(ctauTime2DeltaVal, ctauTime2DeltaRange.first, ctauTime2DeltaRange.second),
		ctauTime2DeltaRange.first, ctauTime2DeltaRange.second);
	RooFormulaVar ctauTime2Scale("ctauTime2Scale", "@0+@1", RooArgList(ctauTime1Scale, ctauTime2Delta));
	RooGaussModel ctauTime2("ctauTime2", "ctauTime2", ctau, ctauTime2Mean, ctauTime2Scale);

	const double ctauTime3MeanVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime3Mean", 0.0) : 0.0;
	const double ctauTime3ScaleVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime3Scale", ctauTime2ScaleVal + 0.05) : ctauTime2ScaleVal + 0.05;
	const double ctauTime3DeltaVal = std::max(0.001, ctauTime3ScaleVal - ctauTime2ScaleVal);
	const auto ctauTime3DeltaRange = fitRangeAround(ctauTime3DeltaVal, 0.001, 0.60, 2.0, 0.080);
	RooRealVar ctauTime3Mean("ctauTime3Mean", "ctauTime3Mean", ctauTime3MeanVal, ctauTime3MeanVal - 0.05, ctauTime3MeanVal + 0.05);
	RooRealVar ctauTime3Delta("ctauTime3Delta", "ctauTime3Delta",
		std::clamp(ctauTime3DeltaVal, ctauTime3DeltaRange.first, ctauTime3DeltaRange.second),
		ctauTime3DeltaRange.first, ctauTime3DeltaRange.second);
	RooFormulaVar ctauTime3Scale("ctauTime3Scale", "@0+@1", RooArgList(ctauTime2Scale, ctauTime3Delta));
	RooGaussModel ctauTime3("ctauTime3", "ctauTime3", ctau, ctauTime3Mean, ctauTime3Scale);

	const double ctauTime4MeanVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime4Mean", 0.0) : 0.0;
	const double ctauTime4ScaleVal = useReferenceCtauMC ? read_double(referencePr.get(), "ctauTime4Scale", ctauTime3ScaleVal + 0.3) : ctauTime3ScaleVal + 0.3;
	const double ctauTime4DeltaVal = std::max(0.05, ctauTime4ScaleVal - ctauTime3ScaleVal);
	const auto ctauTime4DeltaRange = fitRangeAround(ctauTime4DeltaVal, 0.05, 1.20, 1.5, 0.200);
	RooRealVar ctauTime4Mean("ctauTime4Mean", "ctauTime4Mean", ctauTime4MeanVal, ctauTime4MeanVal - 0.05, ctauTime4MeanVal + 0.05);
	RooRealVar ctauTime4Delta("ctauTime4Delta", "ctauTime4Delta",
		std::clamp(ctauTime4DeltaVal, ctauTime4DeltaRange.first, ctauTime4DeltaRange.second),
		ctauTime4DeltaRange.first, ctauTime4DeltaRange.second);
	RooFormulaVar ctauTime4Scale("ctauTime4Scale", "@0+@1", RooArgList(ctauTime3Scale, ctauTime4Delta));
	RooGaussModel ctauTime4("ctauTime4", "ctauTime4", ctau, ctauTime4Mean, ctauTime4Scale);

	const double ctauFrac1Val = useReferenceCtauMC ? read_double(referencePr.get(), "ctauFrac1", 0.5) : 0.5;
	const double ctauFrac2Val = useReferenceCtauMC ? read_double(referencePr.get(), "ctauFrac2", 0.2) : 0.2;
	const double ctauFrac3Val = useReferenceCtauMC ? read_double(referencePr.get(), "ctauFrac3", 0.1) : 0.1;
	RooRealVar ctauFrac1("ctauFrac1", "ctauFrac1", std::clamp(ctauFrac1Val, 1e-6, 0.999999), 1e-6, 0.999999);
	RooRealVar ctauFrac2("ctauFrac2", "ctauFrac2", std::clamp(ctauFrac2Val, 1e-6, 0.999999), 1e-6, 0.999999);
	RooRealVar ctauFrac3("ctauFrac3", "ctauFrac3", std::clamp(ctauFrac3Val, 1e-6, 0.999999), 1e-6, 0.999999);
	RooFormulaVar ctauFrac2Model("ctauFrac2Model", "(1-@0)*@1", RooArgList(ctauFrac1, ctauFrac2));
	RooFormulaVar ctauFrac3Model("ctauFrac3Model", "(1-@0)*(1-@1)*@2", RooArgList(ctauFrac1, ctauFrac2, ctauFrac3));

	std::unique_ptr<RooAddModel> promptResolutionModel;
	RooResolutionModel *promptResolutionPtr = &ctauTime1;
	if (useReferenceCtauMC && nResolutionComponents >= 2)
	{
		RooArgList resolutionPdfList;
		RooArgList resolutionFracList;
		resolutionPdfList.add(ctauTime1);
		resolutionPdfList.add(ctauTime2);
		resolutionFracList.add(ctauFrac1);
		if (nResolutionComponents >= 3)
		{
			resolutionPdfList.add(ctauTime3);
			resolutionFracList.add(ctauFrac2Model);
		}
		if (nResolutionComponents >= 4)
		{
			resolutionPdfList.add(ctauTime4);
			resolutionFracList.add(ctauFrac3Model);
		}
		promptResolutionModel = std::make_unique<RooAddModel>("promptResolutionModel", "promptResolutionModel", resolutionPdfList, resolutionFracList);
		promptResolutionPtr = promptResolutionModel.get();
	}

	RooAbsPdf *promptTimePdf = useReferenceCtauMC ? static_cast<RooAbsPdf *>(promptResolutionPtr) : static_cast<RooAbsPdf *>(&fallbackPromptPdf);
	RooResolutionModel &promptResolution = useReferenceCtauMC ? *promptResolutionPtr : static_cast<RooResolutionModel &>(fallbackPromptResolution);
	if (useReferenceCtauMC)
	{
		addCtauConstraint("ctauTime1Mean", ctauTime1Mean, ctauTime1MeanVal, 0.005);
		addCtauConstraint("ctauTime1Scale", ctauTime1Scale, ctauTime1ScaleVal, 0.010);
		if (nResolutionComponents >= 2)
		{
			addCtauConstraint("ctauTime2Mean", ctauTime2Mean, ctauTime2MeanVal, 0.005);
			addCtauConstraint("ctauTime2Delta", ctauTime2Delta, ctauTime2DeltaVal, 0.020);
			addCtauConstraint("ctauFrac1", ctauFrac1, std::clamp(ctauFrac1Val, 1e-6, 0.999999), 0.05);
		}
		if (nResolutionComponents >= 3)
		{
			addCtauConstraint("ctauTime3Mean", ctauTime3Mean, ctauTime3MeanVal, 0.005);
			addCtauConstraint("ctauTime3Delta", ctauTime3Delta, ctauTime3DeltaVal, 0.040);
			addCtauConstraint("ctauFrac2", ctauFrac2, std::clamp(ctauFrac2Val, 1e-6, 0.999999), 0.05);
		}
		if (nResolutionComponents >= 4)
		{
			addCtauConstraint("ctauTime4Mean", ctauTime4Mean, ctauTime4MeanVal, 0.005);
			addCtauConstraint("ctauTime4Delta", ctauTime4Delta, ctauTime4DeltaVal, 0.080);
			addCtauConstraint("ctauFrac3", ctauFrac3, std::clamp(ctauFrac3Val, 1e-6, 0.999999), 0.05);
		}
	}

	const int nSignalSSComponents = useReferenceCtauMC ? std::clamp(read_int(referenceNp.get(), "nSignalSSComponents", 1), 1, 3) : 1;
	const double signalLifetimeFloor = useReferenceCtauMC ? read_double(referenceNp.get(), "signalLifetimeFloor", 2e-3) : 2e-3;
	const double signalLifetime1Val = useReferenceCtauMC ? read_double(referenceNp.get(), "signal_lifetime", 0.1) : 0.7;
	const double signalLifetime2Val = useReferenceCtauMC ? read_double(referenceNp.get(), "signal2_lifetime", std::max(0.2, 2.0 * signalLifetime1Val)) : 1.4;
	const double signalLifetime3Val = useReferenceCtauMC ? read_double(referenceNp.get(), "signal3_lifetime", std::max(0.4, 2.0 * signalLifetime2Val)) : 2.8;
	const double signalLifetime2Floor = useReferenceCtauMC ? read_double(referenceNp.get(), "signalLifetime2Floor", std::max(0.75 * signalLifetimeFloor, 2e-3)) : 2e-3;
	const double signalLifetime3Floor = useReferenceCtauMC ? read_double(referenceNp.get(), "signalLifetime3Floor", std::max(1.00 * signalLifetimeFloor, 2e-3)) : 2e-3;
	const double signalLogLifetimeOffsetVal = useReferenceCtauMC
		? read_double(referenceNp.get(), "signal_log_lifetime_offset", std::log(std::max(signalLifetime1Val - signalLifetimeFloor, 1e-4)))
		: std::log(std::max(signalLifetime1Val - signalLifetimeFloor, 1e-4));
	const double signal2LogLifetimeOffsetVal = useReferenceCtauMC
		? read_double(referenceNp.get(), "signal2_log_lifetime_offset", std::log(std::max(signalLifetime2Val - signalLifetime2Floor, 1e-4)))
		: std::log(std::max(signalLifetime2Val - signalLifetime2Floor, 1e-4));
	const double signal3LogLifetimeOffsetVal = useReferenceCtauMC
		? read_double(referenceNp.get(), "signal3_log_lifetime_offset", std::log(std::max(signalLifetime3Val - signalLifetime3Floor, 1e-4)))
		: std::log(std::max(signalLifetime3Val - signalLifetime3Floor, 1e-4));
	const double signalLogLifetimeMin = useReferenceCtauMC ? signalLogLifetimeOffsetVal - 0.70 : std::log(1e-4);
	const double signalLogLifetimeMax = useReferenceCtauMC ? signalLogLifetimeOffsetVal + 0.40 : std::log(std::max(20.0, ctRangeWidth + 1.0));
	const double signal2LogLifetimeMin = useReferenceCtauMC ? signal2LogLifetimeOffsetVal - 0.60 : std::log(1e-4);
	const double signal2LogLifetimeMax = useReferenceCtauMC ? signal2LogLifetimeOffsetVal + 0.60 : std::log(std::max(20.0, ctRangeWidth + 1.0));
	const double signal3LogLifetimeMin = useReferenceCtauMC ? signal3LogLifetimeOffsetVal - 0.80 : std::log(1e-4);
	const double signal3LogLifetimeMax = useReferenceCtauMC ? signal3LogLifetimeOffsetVal + 0.80 : std::log(std::max(20.0, ctRangeWidth + 1.0));
	RooRealVar signal_log_lifetime_offset(
		"signal_log_lifetime_offset", "log(signal_lifetime - floor)",
		std::clamp(signalLogLifetimeOffsetVal, signalLogLifetimeMin, signalLogLifetimeMax), signalLogLifetimeMin, signalLogLifetimeMax);
	RooFormulaVar signal_lifetime("signal_lifetime", Form("%g + exp(@0)", signalLifetimeFloor), RooArgList(signal_log_lifetime_offset));
	RooRealVar signal2_log_lifetime_offset(
		"signal2_log_lifetime_offset", "log(signal2_lifetime - floor)",
		std::clamp(signal2LogLifetimeOffsetVal, signal2LogLifetimeMin, signal2LogLifetimeMax), signal2LogLifetimeMin, signal2LogLifetimeMax);
	RooFormulaVar signal2_lifetime("signal2_lifetime", Form("%g + exp(@0)", signalLifetime2Floor), RooArgList(signal2_log_lifetime_offset));
	RooRealVar signal3_log_lifetime_offset(
		"signal3_log_lifetime_offset", "log(signal3_lifetime - floor)",
		std::clamp(signal3LogLifetimeOffsetVal, signal3LogLifetimeMin, signal3LogLifetimeMax), signal3LogLifetimeMin, signal3LogLifetimeMax);
	RooFormulaVar signal3_lifetime("signal3_lifetime", Form("%g + exp(@0)", signalLifetime3Floor), RooArgList(signal3_log_lifetime_offset));
	if (useReferenceCtauMC)
	{
		addCtauConstraint("signal_log_lifetime_offset", signal_log_lifetime_offset, signalLogLifetimeOffsetVal, 0.12);
		if (nSignalSSComponents >= 2)
			addCtauConstraint("signal2_log_lifetime_offset", signal2_log_lifetime_offset, signal2LogLifetimeOffsetVal, 0.16);
		if (nSignalSSComponents >= 3)
			addCtauConstraint("signal3_log_lifetime_offset", signal3_log_lifetime_offset, signal3LogLifetimeOffsetVal, 0.20);
	}

	RooDecay signal_ss1_time("signal_ss1_time", "signal_ss1_time", ctau, signal_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss2_time("signal_ss2_time", "signal_ss2_time", ctau, signal2_lifetime, promptResolution, RooDecay::SingleSided);
	RooDecay signal_ss3_time("signal_ss3_time", "signal_ss3_time", ctau, signal3_lifetime, promptResolution, RooDecay::SingleSided);
	const double npYield1 = useReferenceCtauMC ? std::max(0.0, read_double(referenceNp.get(), "NsignalSS1", 1.0)) : 1.0;
	const double npYield2 = useReferenceCtauMC ? std::max(0.0, read_double(referenceNp.get(), "NsignalSS2", 0.0)) : 0.0;
	const double npYield3 = useReferenceCtauMC ? std::max(0.0, read_double(referenceNp.get(), "NsignalSS3", 0.0)) : 0.0;
	const double npYieldSum = std::max(1e-9, npYield1 + npYield2 + npYield3);
	RooRealVar signalSSFrac1("signalSSFrac1", "signalSSFrac1",
		std::clamp(npYield1 / npYieldSum, 1e-6, 0.999999), 1e-6, 0.999999);
	RooRealVar signalSSFrac2("signalSSFrac2", "signalSSFrac2",
		std::clamp(npYield2 / std::max(1e-9, npYield2 + npYield3), 1e-6, 0.999999), 1e-6, 0.999999);
	signalSSFrac1.setConstant(!kFloatSignalSSFractions);
	signalSSFrac2.setConstant(!kFloatSignalSSFractions);
	std::unique_ptr<RooAddPdf> signalNpOwned;
	RooAbsPdf *signalNpTime = useReferenceCtauMC ? static_cast<RooAbsPdf *>(&signal_ss1_time) : static_cast<RooAbsPdf *>(&fallbackNonpromptPdf);
	if (useReferenceCtauMC && nSignalSSComponents == 2)
	{
		signalNpOwned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time",
			RooArgList(signal_ss1_time, signal_ss2_time), RooArgList(signalSSFrac1), true);
		signalNpTime = signalNpOwned.get();
	}
	else if (useReferenceCtauMC && nSignalSSComponents == 3)
	{
		signalNpOwned = std::make_unique<RooAddPdf>("signal_np_time", "signal_np_time",
			RooArgList(signal_ss1_time, signal_ss2_time, signal_ss3_time), RooArgList(signalSSFrac1, signalSSFrac2), true);
		signalNpTime = signalNpOwned.get();
	}

	const double totalYield = std::max(1.0, hSignalYield->Integral());
	RooRealVar Nprompt("Nprompt", "Nprompt", 0.80 * totalYield, 0.0, 2.0 * totalYield);
	RooRealVar Nnonprompt("Nnonprompt", "Nnonprompt", 0.20 * totalYield, 0.0, 2.0 * totalYield);
	RooAddPdf ctauModel("ctauModel", "ctauModel", RooArgList(*promptTimePdf, *signalNpTime), RooArgList(Nprompt, Nnonprompt));

	std::unique_ptr<RooFitResult> ctauFit;
	if (ctauFitSlices.size() >= 2 && totalYield > 0.0)
	{
		if (useReferenceCtauMC && ctauConstraints.getSize() > 0)
		{
			std::cout << Form("[INFO] final ctau fit uses fitTo with IntegrateBins, SumW2Error(false), and %d ExternalConstraints.",
				ctauConstraints.getSize()) << std::endl;
			if (useCenterVeto)
			{
				ctauFit.reset(ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					ExternalConstraints(ctauConstraints), Range(ctauFitRangeName.Data()),
					Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(1)));
			}
			else
			{
				ctauFit.reset(ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					ExternalConstraints(ctauConstraints), Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(1)));
			}
		}
		else
		{
			std::cout << "[INFO] final ctau fit uses fitTo with IntegrateBins and SumW2Error(false)." << std::endl;
			if (useCenterVeto)
			{
				ctauFit.reset(ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					Range(ctauFitRangeName.Data()), Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(1)));
			}
			else
			{
				ctauFit.reset(ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(1)));
			}
		}
		if (!ctauFit || ctauFit->status() != 0 || ctauFit->covQual() < 2)
		{
			std::cout << Form("[INFO] final ctau constrained fitTo status=%d covQual=%d; retrying strategy 2.",
				ctauFit ? ctauFit->status() : -999, ctauFit ? ctauFit->covQual() : -999) << std::endl;
			if (useReferenceCtauMC && ctauConstraints.getSize() > 0)
			{
				if (useCenterVeto)
				{
					ctauFit.reset(ctauModel.fitTo(signalYieldData,
						Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
						ExternalConstraints(ctauConstraints), Range(ctauFitRangeName.Data()),
						Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(2)));
				}
				else
				{
					ctauFit.reset(ctauModel.fitTo(signalYieldData,
						Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
						ExternalConstraints(ctauConstraints), Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(2)));
				}
			}
			else
			{
				if (useCenterVeto)
				{
					ctauFit.reset(ctauModel.fitTo(signalYieldData,
						Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
						Range(ctauFitRangeName.Data()), Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(2)));
				}
				else
				{
					ctauFit.reset(ctauModel.fitTo(signalYieldData,
						Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
						Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(2)));
				}
			}
		}
	}
	else
	{
		std::cout << "[INFO] skipping final ctau fit: need at least two non-empty subrange bins." << std::endl;
	}
	if (ctauFit)
		ctauFit->Print();

	TCanvas cCtau("c_ctau_subrange", "c_ctau_subrange", 800, 800);
	TPad *ctPad1 = new TPad("ctPad1", "ctPad1", 0.0, 0.25, 1.0, 1.0);
	ctPad1->SetBottomMargin(0.00001);
	ctPad1->SetTopMargin(0.08);
	ctPad1->SetFillStyle(4000);
	ctPad1->SetFrameFillStyle(4000);
	ctPad1->SetTicks(1, 1);
	ctPad1->SetLogy();
	ctPad1->Draw();
	ctPad1->cd();

	auto *ctFrame = ctau.frame(Title(""));
	signalYieldData.plotOn(ctFrame, Binning(ctauBinning), Name("yieldData"));
	ctauModel.plotOn(ctFrame, Name("ctauTotal"), LineColor(kBlack), LineWidth(2), Precision(1e-5));
	ctauModel.plotOn(ctFrame, Components(*promptTimePdf), Name("ctauPrompt"), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Precision(1e-5));
	ctauModel.plotOn(ctFrame, Components(*signalNpTime), Name("ctauNonPrompt"), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Precision(1e-5));
	ctFrame->SetFillStyle(4000);
	ctFrame->GetXaxis()->SetTitle("");
	ctFrame->GetXaxis()->SetTitleSize(0);
	ctFrame->GetXaxis()->SetLabelSize(0);
	ctFrame->GetYaxis()->SetTitle("Signal yield from mass fits");
	ctFrame->GetYaxis()->SetTitleOffset(1.6);
	apply_logy_auto_range(ctFrame, "yieldData");
	ctFrame->Draw("e");
	TLegend ctLeg(0.50, 0.70, 0.72, 0.89);
	ctLeg.SetBorderSize(0);
	ctLeg.SetFillStyle(0);
	ctLeg.SetTextSize(0.03);
	auto findCtObj = [&](const char *n) -> TObject * { return ctFrame ? ctFrame->findObject(n) : nullptr; };
	if (auto *o = findCtObj("yieldData"))
		ctLeg.AddEntry(o, "Mass-fit signal yield", "lep");
	if (auto *o = findCtObj("ctauTotal"))
		ctLeg.AddEntry(o, "Fit", "l");
	if (auto *o = findCtObj("ctauPrompt"))
		ctLeg.AddEntry(o, "PR", "l");
	if (auto *o = findCtObj("ctauNonPrompt"))
		ctLeg.AddEntry(o, "NP", "l");
	ctLeg.Draw("same");

	const double fitYieldSum = std::max(1e-9, Nprompt.getVal() + Nnonprompt.getVal());
	const double bFractionVal = Nnonprompt.getVal() / fitYieldSum;
	const double bFractionErr = std::sqrt(std::pow(Nprompt.getVal() * Nnonprompt.getError(), 2) +
		std::pow(Nnonprompt.getVal() * Nprompt.getError(), 2)) / std::pow(fitYieldSum, 2);
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
		tx.SetTextAlign(11);
		tx.SetTextSize(0.04);
		tx.SetTextFont(62);
		tx.DrawLatex(0.23, 0.930, "CMS");
		tx.SetTextFont(72);
		tx.DrawLatex(0.32, 0.930, "Internal");
	}
	{
		TLatex tx;
		tx.SetNDC();
		tx.SetTextSize(0.03);
		tx.SetTextFont(42);
		double xtext = 0.19, y0 = 0.865, dy = -0.05;
		int k = 0;
		tx.DrawLatex(xtext, y0 + dy * k++, "Data J/#psi #rightarrow #mu^{+}#mu^{-}");
		if (yLow == 0.0f)
			tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
		else
			tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
	}
	{
		TLatex tp;
		tp.SetNDC();
		tp.SetTextSize(0.026);
		tp.SetTextFont(42);
		double xtext = 0.74, y0 = 0.87, dy = -0.045;
		int k = 0;
		tp.DrawLatex(xtext, y0 + dy * k++, Form("f_{B} = %.4g #pm %.3g", bFractionVal, bFractionErr));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{sig} = %.4g", hSignalYield->Integral()));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{PR} = %.4g", Nprompt.getVal()));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{NP} = %.4g", Nnonprompt.getVal()));
		if (ctauFit)
			tp.DrawLatex(xtext, y0 + dy * k++, Form("status/cov = %d/%d", ctauFit->status(), ctauFit->covQual()));
	}

	cCtau.cd();
	TPad *ctPad2 = new TPad("ctPad2", "ctPad2", 0.0, 0.0, 1.0, 0.25);
	ctPad2->SetTopMargin(0.00001);
	ctPad2->SetBottomMargin(0.4);
	ctPad2->SetFillStyle(4000);
	ctPad2->SetFrameFillStyle(4000);
	ctPad2->SetTicks(1, 1);
	ctPad2->Draw();
	ctPad2->cd();

	RooPlot *ctPullPlot = ctau.frame(Title(""));
	RooHist *ctPull = ctFrame->pullHist("yieldData", "ctauTotal", true);
	if (ctPull)
		ctPull->SetMarkerSize(0.8);
	if (ctPull)
		ctPullPlot->addPlotable(ctPull, "P");
	double maxAbsPull = 0.0;
	double maxPullX = 0.0;
	double maxPullY = 0.0;
	if (ctPull)
	{
		double px = 0.0, py = 0.0;
		for (int ip = 0; ip < ctPull->GetN(); ++ip)
		{
			ctPull->GetPoint(ip, px, py);
			if (std::isfinite(py) && std::abs(py) > maxAbsPull)
			{
				maxAbsPull = std::abs(py);
				maxPullX = px;
				maxPullY = py;
			}
		}
	}
	const double pullAxisMax = std::max(8.0, std::min(80.0, std::ceil(maxAbsPull + 1.0)));
	std::cout << Form("[FIT] subrange ctau max |pull| = %.3f at ctau %.5g (pull %.3f), pull axis = +/-%.1f",
		maxAbsPull, maxPullX, maxPullY, pullAxisMax) << std::endl;
	ctPullPlot->GetYaxis()->SetTitle("Pull");
	ctPullPlot->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} [mm]");
	ctPullPlot->GetXaxis()->CenterTitle();
	ctPullPlot->SetMinimum(-pullAxisMax);
	ctPullPlot->SetMaximum(pullAxisMax);
	ctPullPlot->GetYaxis()->SetNdivisions(505);
	ctPullPlot->GetYaxis()->SetTitleSize(0.12);
	ctPullPlot->GetYaxis()->SetLabelSize(0.10);
	ctPullPlot->GetXaxis()->SetTitleSize(0.15);
	ctPullPlot->GetXaxis()->SetLabelSize(0.10);
	ctPullPlot->Draw();

	auto ctChi = chi2_from_pull(ctPull);
	int ctNpar = ctauFit ? ctauFit->floatParsFinal().getSize() : 0;
	int ctNdf = std::max(1, ctChi.second - ctNpar);
	double ctPvalue = TMath::Prob(ctChi.first, ctNdf);
	std::cout << Form("[FIT] subrange ctau chi2/ndf = %.1f/%d = %.3f, p = %.3g",
		ctChi.first, ctNdf, ctChi.first / ctNdf, ctPvalue) << std::endl;
	{
		TLatex tc;
		tc.SetNDC();
		tc.SetTextSize(0.085);
		tc.SetTextFont(42);
		tc.SetTextAlign(33);
		tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", ctChi.first, ctNdf));
	}
	TLine ctLine(histEdges.front(), 0.0, histEdges.back(), 0.0);
	ctLine.SetLineStyle(2);
	ctLine.Draw("same");

	const TString ctauFigName = TString::Format("%s/ctau_signal_yield_%s.pdf", figDir.Data(), outputTag.Data());
	cCtau.SaveAs(ctauFigName);
	delete ctPullPlot;

	TFile outFile(outputName, "RECREATE");
	if (!outFile.IsZombie())
	{
		double ctLowOut = 0.0;
		double ctHighOut = 0.0;
		int entriesOut = 0;
		double nsigOut = 0.0;
		double nsigErrOut = 0.0;
		double nbkgOut = 0.0;
		double nbkgErrOut = 0.0;
		double nsigSignificanceOut = 0.0;
		int signalWindowEntriesOut = 0;
		int shapeDrivenSignalOut = 0;
		int usedInCtauFitOut = 0;
		int statusOut = 0;
		int covQualOut = 0;
		TTree summary("sliceSummary", "mass subrange fit summary");
		summary.Branch("ctLow", &ctLowOut, "ctLow/D");
		summary.Branch("ctHigh", &ctHighOut, "ctHigh/D");
		summary.Branch("entries", &entriesOut, "entries/I");
		summary.Branch("Nsig", &nsigOut, "Nsig/D");
		summary.Branch("NsigErr", &nsigErrOut, "NsigErr/D");
		summary.Branch("Nbkg", &nbkgOut, "Nbkg/D");
		summary.Branch("NbkgErr", &nbkgErrOut, "NbkgErr/D");
		summary.Branch("NsigSignificance", &nsigSignificanceOut, "NsigSignificance/D");
		summary.Branch("signalWindowEntries", &signalWindowEntriesOut, "signalWindowEntries/I");
		summary.Branch("shapeDrivenSignal", &shapeDrivenSignalOut, "shapeDrivenSignal/I");
		summary.Branch("usedInCtauFit", &usedInCtauFitOut, "usedInCtauFit/I");
		summary.Branch("fitStatus", &statusOut, "fitStatus/I");
		summary.Branch("covQual", &covQualOut, "covQual/I");
		for (const auto &r : sliceResults)
		{
			ctLowOut = r.ctLow;
			ctHighOut = r.ctHigh;
			entriesOut = r.entries;
			nsigOut = r.nsig;
			nsigErrOut = r.nsigErr;
			nbkgOut = r.nbkg;
			nbkgErrOut = r.nbkgErr;
			nsigSignificanceOut = r.nsigSignificance;
			signalWindowEntriesOut = r.signalWindowEntries;
			shapeDrivenSignalOut = r.shapeDrivenSignal;
			usedInCtauFitOut = r.usedInCtauFit;
			statusOut = r.status;
			covQualOut = r.covQual;
			summary.Fill();
		}
		hSignalYield->Write();
		summary.Write();
		if (ctauFit)
			ctauFit->Write("ctauFitResult");
		TParameter<double>("ptLow", ptLow).Write();
		TParameter<double>("ptHigh", ptHigh).Write();
		TParameter<double>("yLow", yLow).Write();
		TParameter<double>("yHigh", yHigh).Write();
		TParameter<double>("ctMin", histEdges.front()).Write();
		TParameter<double>("ctMax", histEdges.back()).Write();
		TParameter<double>("ctStep", ctStep).Write();
		TParameter<int>("nMassSlices", static_cast<int>(sliceResults.size())).Write();
		TParameter<int>("nCtauFitSlices", static_cast<int>(ctauFitSlices.size())).Write();
		TParameter<int>("adaptiveCtauBinning", 1).Write();
		TParameter<double>("adaptiveMinBinWidth", minEdgeWidth).Write();
		TParameter<double>("adaptiveMaxBinWidth", maxEdgeWidth).Write();
		TParameter<double>("adaptiveMaxNeighborWidthRatio", maxNeighborWidthRatio).Write();
		TParameter<double>("edgeSignalThreshold", edgeSignalThreshold).Write();
		TParameter<double>("ctauFitMinOverride", kCtauFitMinOverride).Write();
		TParameter<double>("centerVetoHalfWidth", centerVeto).Write();
		TParameter<int>("useCenterVeto", useCenterVeto ? 1 : 0).Write();
		TParameter<double>("integrateBinsPrecision", 1e-6).Write();
		TParameter<int>("ctauFitUsesChi2", 0).Write();
		TParameter<int>("ctauFitSumW2Error", 0).Write();
		TParameter<int>("nCtauConstraintsBuilt", ctauConstraints.getSize()).Write();
		TParameter<int>("nCtauConstraintsApplied", (useReferenceCtauMC && ctauConstraints.getSize() > 0) ? ctauConstraints.getSize() : 0).Write();
		TParameter<double>("ctauConstraintScale", kCtauConstraintScale).Write();
		TParameter<int>("disableCtauConstraints", kDisableCtauConstraints ? 1 : 0).Write();
		TParameter<int>("floatSignalSSFractions", kFloatSignalSSFractions ? 1 : 0).Write();
		TParameter<int>("trimSparseFwd", trimSparseFwd ? 1 : 0).Write();
		TParameter<int>("applyQuantileCut", applyQuantileCut ? 1 : 0).Write();
		TParameter<int>("reuseMassSlices", reuseMassSlices ? 1 : 0).Write();
		TParameter<double>("Nprompt", Nprompt.getVal()).Write();
		TParameter<double>("NpromptErr", Nprompt.getError()).Write();
		TParameter<double>("Nnonprompt", Nnonprompt.getVal()).Write();
		TParameter<double>("NnonpromptErr", Nnonprompt.getError()).Write();
		TParameter<double>("bFraction", bFractionVal).Write();
		TParameter<double>("bFractionErr", bFractionErr).Write();
		TParameter<double>("ctauChi2", ctChi.first).Write();
		TParameter<int>("ctauNdf", ctNdf).Write();
		TParameter<double>("ctauPvalue", ctPvalue).Write();
		TParameter<double>("promptSigma", useReferenceCtauMC ? ctauTime1Scale.getVal() : fallbackPromptSigma.getVal()).Write();
		TParameter<double>("nonpromptTau", useReferenceCtauMC ? signal_lifetime.getVal() : fallbackNonpromptTau.getVal()).Write();
		TParameter<int>("useReferenceCtauMC", useReferenceCtauMC ? 1 : 0).Write();
		TParameter<int>("nResolutionComponents", nResolutionComponents).Write();
		TParameter<int>("nSignalSSComponents", nSignalSSComponents).Write();
		TParameter<double>("signalSSFrac1", signalSSFrac1.getVal()).Write();
		TParameter<double>("signalSSFrac2", signalSSFrac2.getVal()).Write();
		outFile.Write();
	}

	std::cout << "[DONE] wrote " << outputName << std::endl;
	std::cout << "[DONE] wrote " << ctauFigName << std::endl;
	delete ctFrame;
}
