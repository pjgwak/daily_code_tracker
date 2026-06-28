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
#include "subrange_ctau_config.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
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

using pp_jpsi_subrange::ManualCtauBinningOverride;
using pp_jpsi_subrange::default_ctau_binning_note;
using pp_jpsi_subrange::has_default_ctau_binning;
using pp_jpsi_subrange::find_manual_ctau_binning_override;
using pp_jpsi_subrange::find_manual_ctau_fit_range;
using pp_jpsi_subrange::target_ctau_bin_width;

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

// Compact record for one ctau subrange mass fit. The same fields are written
// into each per-slice ROOT file so failed or suspicious slices can be rerun
// independently and then collected by the parent job.
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

// Stable file naming for per-slice mass fits. The name is based on the ctau
// edges, not the loop index, so a single-slice rerun overwrites exactly the
// artifact that the full collection step will later read.
TString mass_slice_root_name(const TString &rootDir, const TString &tag, double ctLow, double ctHigh)
{
	const TString massRootDir = TString::Format("%s/mass_slices/%s", rootDir.Data(), tag.Data());
	TString name = TString::Format("%s/mass_fit_%s_ct%.3f_%.3f.root",
		massRootDir.Data(), tag.Data(), ctLow, ctHigh);
	name.ReplaceAll("-", "m");
	return name;
}

TString mass_slice_fig_name(const TString &figDir, const TString &tag, double ctLow, double ctHigh)
{
	const TString massFigDir = TString::Format("%s/mass_slices/%s", figDir.Data(), tag.Data());
	TString name = TString::Format("%s/%s_ct%.3f_%.3f.pdf",
		massFigDir.Data(), tag.Data(), ctLow, ctHigh);
	name.ReplaceAll("-", "m");
	return name;
}

void write_empty_slice_plot(const TString &figDir, const TString &tag, float ptLow, float ptHigh, float yLow, float yHigh, double ctLow, double ctHigh)
{
	const TString massFigDir = TString::Format("%s/mass_slices/%s", figDir.Data(), tag.Data());
	gSystem->mkdir(massFigDir, true);
	TCanvas c("c_empty_mass_slice", "c_empty_mass_slice", 800, 800);
	c.SetTicks(1, 1);
	c.SetTopMargin(0.08);
	c.SetRightMargin(0.05);
	c.SetLeftMargin(0.13);
	c.SetBottomMargin(0.12);
	TH1D frame("h_empty_mass_slice", "", 80, 2.6, 3.5);
	frame.SetMinimum(0.0);
	frame.SetMaximum(1.0);
	frame.GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
	frame.GetYaxis()->SetTitle("Events");
	frame.GetYaxis()->SetTitleOffset(1.35);
	frame.Draw("hist");

	TLatex tx;
	tx.SetNDC();
	tx.SetTextFont(42);
	tx.SetTextSize(0.032);
	tx.SetTextAlign(31);
	tx.DrawLatex(0.96, 0.935, "pp #sqrt{s} = 5.02 TeV (28.0 pb^{-1})");
	tx.SetTextAlign(11);
	tx.SetTextFont(62);
	tx.SetTextSize(0.040);
	tx.DrawLatex(0.18, 0.930, "CMS");
	tx.SetTextFont(52);
	tx.DrawLatex(0.27, 0.930, "Internal");
	tx.SetTextFont(42);
	tx.SetTextSize(0.030);
	double y0 = 0.84;
	double dy = -0.050;
	int k = 0;
	tx.DrawLatex(0.18, y0 + dy * k++, "J/#psi data");
	if (yLow == 0.0f)
		tx.DrawLatex(0.18, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
	else
		tx.DrawLatex(0.18, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
	tx.DrawLatex(0.18, y0 + dy * k++, Form("%.3f < ctau3D < %.3f mm", ctLow, ctHigh));
	tx.DrawLatex(0.18, y0 + dy * k++, "empty slice: no selected entries");

	c.SaveAs(mass_slice_fig_name(figDir, tag, ctLow, ctHigh));
}

// Read back a per-slice mass fit result produced either by the serial path,
// a parallel child ROOT process, or a targeted rerun from run_pp_jpsi_chain.sh.
bool read_slice_result_from_mass_file(const TString &rootDir, const TString &tag, double ctLow, double ctHigh, SliceResult &out)
{
	out = SliceResult();
	out.ctLow = ctLow;
	out.ctHigh = ctHigh;
	const TString fileName = mass_slice_root_name(rootDir, tag, ctLow, ctHigh);
	std::unique_ptr<TFile> file(TFile::Open(fileName, "READ"));
	if (!file || file->IsZombie())
	{
		std::cerr << "ERROR: missing mass-slice output: " << fileName << std::endl;
		return false;
	}
	out.entries = read_int(file.get(), "entries", 0);
	out.nsig = read_double(file.get(), "Nsig", 0.0);
	out.nsigErr = read_double(file.get(), "NsigErr", 0.0);
	out.nsigSignificance = read_double(file.get(), "NsigSignificance", 0.0);
	out.signalWindowEntries = read_int(file.get(), "signalWindowEntries", 0);
	out.shapeDrivenSignal = read_int(file.get(), "shapeDrivenSignal", 0);
	out.nbkg = read_double(file.get(), "Nbkg", 0.0);
	out.nbkgErr = read_double(file.get(), "NbkgErr", 0.0);
	out.status = read_int(file.get(), "fitStatus", -99);
	out.covQual = read_int(file.get(), "covQual", -99);
	out.fitAttempt = read_int(file.get(), "fitAttempt", 0);
	return true;
}

// Empty ctau slices are still materialized as ROOT files. This keeps parallel
// collection deterministic: missing files mean a failed child job, while empty
// files mean a valid slice with no selected entries.
void write_empty_slice_result_file(const TString &rootDir, const TString &tag, const SliceResult &out)
{
	const TString massRootDir = TString::Format("%s/mass_slices/%s", rootDir.Data(), tag.Data());
	gSystem->mkdir(massRootDir, true);
	const TString massRootName = mass_slice_root_name(rootDir, tag, out.ctLow, out.ctHigh);
	TFile massOut(massRootName, "RECREATE");
	if (massOut.IsZombie())
		return;
	TParameter<double>("ctLow", out.ctLow).Write();
	TParameter<double>("ctHigh", out.ctHigh).Write();
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
	massOut.Write();
}

const char *bool_literal(bool value)
{
	return value ? "true" : "false";
}

// Fit one mass distribution inside a single ctau subrange. This function is
// intentionally self-contained because it is called in three different modes:
// serial full production, child ROOT processes for parallel production, and
// targeted one/two-slice reruns from run_pp_jpsi_chain.sh.
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

	// ------------------------------------------------------------------
	// event selection for this ctau slice
	// ------------------------------------------------------------------
	TString sliceCut = Form(
		"(mass > 2.6 && mass < 3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && "
		"(ctau3D >= %.12g && ctau3D < %.12g)",
		ptLow, ptHigh, yLow, yHigh, ctLow, ctHigh);
	auto dataSlice = std::unique_ptr<RooDataSet>(static_cast<RooDataSet *>(inputData.reduce(sliceCut)));
	if (!dataSlice || dataSlice->numEntries() <= 0)
	{
		// Preserve a successful, zero-entry slice result so targeted reruns and
		// parallel collection can distinguish it from an execution failure.
		write_empty_slice_result_file(rootDir, tag, out);
		if (drawMassSlices)
			write_empty_slice_plot(figDir, tag, ptLow, ptHigh, yLow, yHigh, ctLow, ctHigh);
		return out;
	}
	out.entries = dataSlice->numEntries();

	// ------------------------------------------------------------------
	// mass model setup
	// ------------------------------------------------------------------
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

	// ------------------------------------------------------------------
	// mass fit and slice-quality bookkeeping
	// ------------------------------------------------------------------
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

	// ------------------------------------------------------------------
	// persist per-slice result
	// ------------------------------------------------------------------
	const TString massRootDir = TString::Format("%s/mass_slices/%s", rootDir.Data(), tag.Data());
	gSystem->mkdir(massRootDir, true);
	const TString massRootName = mass_slice_root_name(rootDir, tag, ctLow, ctHigh);
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

	// ------------------------------------------------------------------
	// optional per-slice diagnostic plot
	// ------------------------------------------------------------------
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
			tx.SetTextFont(52);
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

		c.SaveAs(mass_slice_fig_name(figDir, tag, ctLow, ctHigh));
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

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
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

	// ------------------------------------------------------------------
	// event selection and quantile cleanup
	// ------------------------------------------------------------------
	// The first reduction applies the physics selection. The optional quantile
	// cleanup removes pathological ctau/ctauErr tails before any subrange mass fit,
	// while keeping the nominal ctMin/ctMax grid independent from outliers.
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

	// ------------------------------------------------------------------
	// output naming and run-control tags
	// ------------------------------------------------------------------
	// Keep naming consistent with mc_mass/mass/ctau_pr/ctau_np outputs so the
	// reference files can be found by bin tag and the chain script can rerun only
	// selected pieces without changing downstream paths.
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

	// Raw diagnostic output stores only histograms and counting metadata; it does
	// not feed the production subrange fit.
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
		tx.SetTextFont(52);
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

// Mass-subrange workflow. The first arguments define the physics bin and ctau
// binning. This macro fits signal yields in ctau slices and writes the
// sliceSummary consumed by subrange_ctau.C.
//   firstSlice/nSlicesToFit : run a contiguous subset of adaptive ctau slices.
//                      With massSlicesOnly=true this is the targeted rerun path
//                      for suspicious subrange mass fits.
//   parallelMassSlices/maxParallelJobs : fan out full mass-slice fitting to child
//                      ROOT processes, then collect their per-slice ROOT files.
//   massSlicesOnly   : stop after writing per-slice mass fit ROOT files and skip
//                      writing a partial sliceSummary.
void subrange_mass(float ptLow = 3.0, float ptHigh = 6.5, float yLow = 1.6, float yHigh = 2.4,
	double ctMin = -8.0, double ctMax = 10.0, double ctStep = 2.0,
	int maxSlices = -1, bool drawMassSlices = true, bool trimSparseFwd = true,
	bool applyQuantileCut = true, bool reuseMassSlices = false, double centerVetoHalfWidth = 0.0,
	int firstSlice = 0, int nSlicesToFit = -1, bool parallelMassSlices = false, int maxParallelJobs = 4,
	bool massSlicesOnly = false)
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
	const TString baseOutputTag = tag + (trimSparseFwd ? "" : "_fullrange") + (applyQuantileCut ? "" : "_noQuantile");
	// Per-slice mass files always use the base tag so a targeted rerun writes
	// back into the same mass_slices directory used by the full production run.
	// Only diagnostic final-ctau subset runs get an extra outputTag suffix.
	const TString massSliceTag = baseOutputTag;
	TString outputTag = baseOutputTag;
	if (!massSlicesOnly && (firstSlice > 0 || nSlicesToFit > 0))
	{
		outputTag += TString::Format("_slice%d", std::max(0, firstSlice));
		if (nSlicesToFit > 0)
			outputTag += TString::Format("_n%d", nSlicesToFit);
	}
	const TString figDir = TString::Format("figs/%s/subrange", yTag.Data());
	const TString rootDir = TString::Format("roots/%s/subrange", yTag.Data());
	gSystem->mkdir(figDir, true);
	gSystem->mkdir(rootDir, true);
	const TString outputName = TString::Format("%s/subrange_result_%s.root", rootDir.Data(), outputTag.Data());

	// ------------------------------------------------------------------
	// reference mass-shape input
	// ------------------------------------------------------------------
	// The inclusive mass fit provides the nominal signal-shape parameters. Each
	// narrow ctau slice can then focus on yield extraction instead of freely
	// rediscovering unstable shape parameters from low statistics.
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

	// ------------------------------------------------------------------
	// adaptive ctau subrange binning
	// ------------------------------------------------------------------
	// The grid is fine near ctau=0, where prompt resolution changes rapidly, and
	// coarser in the tails, where statistics are sparse. The outer ctStep controls
	// only the far-tail band; central bins remain fixed to preserve resolution.
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
	if (const auto *manualRange = find_manual_ctau_fit_range(ptLow, ptHigh, yLow, yHigh))
	{
		const double requestedLow = manualRange->ctLow;
		const double requestedHigh = manualRange->ctHigh;
		const double manualFitCtMin = std::max(fitCtMin, requestedLow);
		const double manualFitCtMax = std::min(fitCtMax, requestedHigh);
		if (!(manualFitCtMin < manualFitCtMax))
		{
			std::cerr << Form("ERROR: enabled manual ctau mass-slice range [%.6g, %.6g] does not overlap nominal range [%.6g, %.6g] for this bin.",
				requestedLow, requestedHigh, fitCtMin, fitCtMax) << std::endl;
			return;
		}
		std::cout << Form("[INFO] manual ctau mass-slice range override applied: requested [%.6g, %.6g], using [%.6g, %.6g]",
			requestedLow, requestedHigh, manualFitCtMin, manualFitCtMax) << std::endl;
		fitCtMin = manualFitCtMin;
		fitCtMax = manualFitCtMax;
	}

	auto appendBand = [&](double start, double stop, double targetWidth) {
		start = std::max(start, fitCtMin);
		stop = std::min(stop, fitCtMax);
		if (stop <= start || targetWidth <= 0.0)
			return;
		const int nBins = std::max(1, static_cast<int>(std::ceil((stop - start) / targetWidth - 1e-9)));
		appendFixedBins(start, stop, nBins);
	};

	const double outerStep = std::max(0.005, std::min(ctStep, 0.75));
	appendBand(fitCtMin, -3.0, outerStep);
	appendBand(-3.0, -2.5, 0.250);
	appendBand(-2.5, -2.0, 0.250);
	appendBand(-2.0, -1.45, 0.090);
	appendBand(-1.45, -1.05, 0.060);
	appendBand(-1.05, -0.75, 0.060);
	appendBand(-0.75, -0.50, 0.050);
	appendBand(-0.50, -0.30, 0.015);
	appendBand(-0.30, -0.15, 0.010);
	appendBand(-0.15, -0.06, 0.0075);
	appendBand(-0.06, 0.06, 0.0050);
	appendBand(0.06, 0.15, 0.0075);
	appendBand(0.15, 0.30, 0.010);
	appendBand(0.30, 0.50, 0.015);
	appendBand(0.50, 0.75, 0.040);
	appendBand(0.75, 1.05, 0.060);
	appendBand(1.05, 1.45, 0.060);
	appendBand(1.45, 2.00, 0.090);
	appendBand(2.00, 2.50, 0.250);
	appendBand(2.50, 3.00, 0.500);
	appendBand(3.00, fitCtMax, outerStep);

	const ManualCtauBinningOverride *ctauBinningOverride = find_manual_ctau_binning_override(ptLow, ptHigh, yLow, yHigh);
	const bool optimizedCtauBinningApplied = ctauBinningOverride != nullptr || has_default_ctau_binning();
	const char *ctauBinningNote = ctauBinningOverride ? ctauBinningOverride->note : default_ctau_binning_note();
	if (optimizedCtauBinningApplied && edges.size() >= 2)
	{
		const size_t beforeMerge = edges.size() - 1;
		std::vector<double> mergedEdges;
		mergedEdges.reserve(edges.size());
		mergedEdges.push_back(edges.front());
		for (size_t i = 0; i + 1 < edges.size(); )
		{
			const double start = edges[i];
			const double width0 = edges[i + 1] - edges[i];
			const double center = 0.5 * (edges[i] + edges[i + 1]);
			const double targetWidth = target_ctau_bin_width(ctauBinningOverride, center, width0);
			size_t j = i + 1;
			while (j + 1 < edges.size() && edges[j] - start < targetWidth - 1e-9)
				++j;
			mergedEdges.push_back(edges[j]);
			i = j;
		}
		edges.swap(mergedEdges);
		std::cout << Form("[INFO] optimized ctau mass-slice binning for %s: merged %zu -> %zu bins; note: %s",
			bin_tag(ptLow, ptHigh, yLow, yHigh).Data(), beforeMerge, edges.size() - 1, ctauBinningNote) << std::endl;
	}

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

	// sliceIndices preserves the original adaptive-slice index after edges are
	// narrowed. This matters for child-process logs and for targeted reruns such
	// as firstSlice=100,nSlicesToFit=1.
	std::vector<int> sliceIndices;
	sliceIndices.reserve(edges.size() > 0 ? edges.size() - 1 : 0);
	for (size_t i = 0; i + 1 < edges.size(); ++i)
		sliceIndices.push_back(static_cast<int>(i));
	if (maxSlices > 0 && maxSlices + 1 < static_cast<int>(edges.size()))
	{
		edges.resize(maxSlices + 1);
		sliceIndices.resize(maxSlices);
	}
	const std::vector<double> summaryEdges = edges;
	if ((firstSlice > 0 || nSlicesToFit > 0) && edges.size() >= 2)
	{
		// Restrict the adaptive ctau grid to the requested contiguous slice range.
		// This is the low-level primitive behind --subrange-slices in the run script.
		const int totalSlices = static_cast<int>(edges.size()) - 1;
		const int selectedFirst = std::clamp(firstSlice, 0, totalSlices);
		const int selectedLastExclusive = nSlicesToFit > 0 ? std::min(totalSlices, selectedFirst + nSlicesToFit) : totalSlices;
		if (selectedFirst >= selectedLastExclusive)
		{
			std::cerr << Form("ERROR: requested slice selection first=%d n=%d is empty for %d available slices",
				firstSlice, nSlicesToFit, totalSlices) << std::endl;
			return;
		}
		std::vector<double> selectedEdges(edges.begin() + selectedFirst, edges.begin() + selectedLastExclusive + 1);
		std::vector<int> selectedIndices(sliceIndices.begin() + selectedFirst, sliceIndices.begin() + selectedLastExclusive);
		edges.swap(selectedEdges);
		sliceIndices.swap(selectedIndices);
		std::cout << Form("[INFO] selected ctau slice range: first=%d n=%zu, ctau [%.6g, %.6g]",
			sliceIndices.front(), sliceIndices.size(), edges.front(), edges.back()) << std::endl;
	}

	// ------------------------------------------------------------------
	// mass-slice execution mode
	// ------------------------------------------------------------------
	// Exactly one of the paths below fills sliceResults:
	//   1) parallelMassSlices: spawn child ROOT jobs, then collect per-slice files.
	//   2) serial/targeted: fit the selected slices directly in this process.
	std::vector<SliceResult> sliceResults;
	if (reuseMassSlices)
		std::cout << "[INFO] subrange_mass ignores reuseMassSlices=true; running mass-slice fits." << std::endl;
	sliceResults.reserve(edges.size() > 0 ? edges.size() - 1 : 0);
	if (parallelMassSlices && !massSlicesOnly)
		{
			// Full production acceleration path. Each child ROOT process runs this same
			// macro in massSlicesOnly mode for exactly one original slice, writes a
			// per-slice ROOT file, and exits. The parent then reads those files back in
			// ctau-edge order and writes the mass-slice summary.
			const int jobs = std::max(1, maxParallelJobs);
			const TString logDir = TString::Format("%s/logs/subrange_parallel_%s", kReferenceDir, massSliceTag.Data());
			gSystem->mkdir(logDir, true);
			std::cout << Form("[INFO] running %zu mass-slice fits in parallel ROOT subprocesses, maxParallelJobs=%d",
				edges.size() > 0 ? edges.size() - 1 : 0, jobs) << std::endl;
			int failures = 0;
			for (size_t batchStart = 0; batchStart + 1 < edges.size(); batchStart += static_cast<size_t>(jobs))
			{
				std::ostringstream script;
				script << "pids=; status=0; ";
				const size_t batchEnd = std::min(edges.size() - 1, batchStart + static_cast<size_t>(jobs));
				for (size_t i = batchStart; i < batchEnd; ++i)
				{
					const int originalSlice = sliceIndices[i];
					const TString logName = TString::Format("%s/slice_%04d.log", logDir.Data(), originalSlice);
					const TString expr = TString::Format(
						"subrange_mass(%g,%g,%g,%g,%.17g,%.17g,%.17g,-1,%s,%s,%s,false,%.17g,%d,1,false,%d,true)",
						ptLow, ptHigh, yLow, yHigh, ctMin, ctMax, ctStep,
						bool_literal(drawMassSlices), bool_literal(trimSparseFwd), bool_literal(applyQuantileCut),
						centerVetoHalfWidth, originalSlice, jobs);
					script << "root -l -b -q -e \".L " << kReferenceDir << "/subrange_mass.C\" -e \""
						   << expr.Data() << "\" > " << logName.Data() << " 2>&1 & pids=\"$pids $!\"; ";
				}
				script << "for pid in $pids; do wait $pid || status=1; done; exit $status";
				TString batchScript(script.str().c_str());
				batchScript.ReplaceAll("'", "'\\''");
				const TString cmd = TString::Format("bash -lc '%s'", batchScript.Data());
				const int rc = gSystem->Exec(cmd);
				if (rc != 0)
				{
					std::cerr << Form("ERROR: parallel mass-slice batch failed with status %d; see %s", rc, logDir.Data()) << std::endl;
					++failures;
				}
			}
			if (failures > 0)
				return;
			for (size_t i = 0; i + 1 < edges.size(); ++i)
			{
				SliceResult r;
				if (!read_slice_result_from_mass_file(rootDir, massSliceTag, edges[i], edges[i + 1], r))
					return;
				sliceResults.push_back(r);
				std::cout << Form("[SLICE %zu/%zu] ctau %.3f to %.3f: entries=%d peakWin=%d Nsig=%.3f +/- %.3f status=%d covQual=%d%s",
					i + 1, edges.size() - 1, r.ctLow, r.ctHigh, r.entries, r.signalWindowEntries,
					r.nsig, r.nsigErr, r.status, r.covQual, r.shapeDrivenSignal ? " shape-veto" : "") << std::endl;
			}
		}
	else
	{
		for (size_t i = 0; i + 1 < edges.size(); ++i)
			{
				std::cout << Form("[SLICE %zu/%zu] ctau %.3f to %.3f", i + 1, edges.size() - 1, edges[i], edges[i + 1]) << std::endl;
				sliceResults.push_back(fit_mass_slice(*dataSel, *massVar, ptLow, ptHigh, yLow, yHigh,
					edges[i], edges[i + 1], referenceMass.get(), figDir, rootDir, massSliceTag, drawMassSlices));
				const auto &r = sliceResults.back();
				std::cout << Form("  entries=%d peakWin=%d Nsig=%.3f +/- %.3f status=%d covQual=%d%s",
					r.entries, r.signalWindowEntries, r.nsig, r.nsigErr, r.status, r.covQual,
					r.shapeDrivenSignal ? " shape-veto" : "") << std::endl;
			}
	}

	if (sliceResults.empty())
	{
		std::cerr << "ERROR: no mass subrange slices were fitted." << std::endl;
		return;
	}

	auto writeSliceSummary = [&](const std::vector<SliceResult> &summaryResults, int summaryFirstSlice, int summaryNSlicesToFit,
		bool collectedFromMassSliceFiles) -> bool {
		if (summaryResults.empty())
		{
			std::cerr << "ERROR: no mass subrange slices available for summary output." << std::endl;
			return false;
		}
		TFile outFile(outputName, "RECREATE");
		if (outFile.IsZombie())
		{
			std::cerr << "ERROR: cannot create mass-slice summary: " << outputName << std::endl;
			return false;
		}

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
		int usedInCtauFitOut = 1;
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
		for (const auto &r : summaryResults)
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
		summary.Write();
		TParameter<double>("ptLow", ptLow).Write();
		TParameter<double>("ptHigh", ptHigh).Write();
		TParameter<double>("yLow", yLow).Write();
		TParameter<double>("yHigh", yHigh).Write();
		TParameter<double>("ctMin", fitCtMin).Write();
		TParameter<double>("ctMax", fitCtMax).Write();
		TParameter<double>("ctStep", ctStep).Write();
		TParameter<int>("nMassSlices", static_cast<int>(summaryResults.size())).Write();
		TParameter<int>("adaptiveCtauBinning", 1).Write();
		TParameter<int>("optimizedCtauBinningApplied", optimizedCtauBinningApplied ? 1 : 0).Write();
		TParameter<double>("adaptiveMinBinWidth", minEdgeWidth).Write();
		TParameter<double>("adaptiveMaxBinWidth", maxEdgeWidth).Write();
		TParameter<double>("adaptiveMaxNeighborWidthRatio", maxNeighborWidthRatio).Write();
		TParameter<int>("trimSparseFwd", trimSparseFwd ? 1 : 0).Write();
		TParameter<int>("applyQuantileCut", applyQuantileCut ? 1 : 0).Write();
		TParameter<int>("firstSlice", summaryFirstSlice).Write();
		TParameter<int>("nSlicesToFit", summaryNSlicesToFit).Write();
		TParameter<int>("parallelMassSlices", parallelMassSlices ? 1 : 0).Write();
		TParameter<int>("maxParallelJobs", maxParallelJobs).Write();
		TParameter<int>("collectedFromMassSliceFiles", collectedFromMassSliceFiles ? 1 : 0).Write();
		outFile.Write();
		std::cout << "[DONE] wrote mass-slice summary " << outputName << std::endl;
		return true;
	};

	if (massSlicesOnly)
	{
		bool wroteSummary = false;
		const bool selectedSubset = firstSlice > 0 || nSlicesToFit > 0;
		if (selectedSubset)
		{
			std::vector<SliceResult> collectedResults;
			collectedResults.reserve(summaryEdges.size() > 0 ? summaryEdges.size() - 1 : 0);
			bool collectOk = true;
			for (size_t i = 0; i + 1 < summaryEdges.size(); ++i)
			{
				SliceResult r;
				if (!read_slice_result_from_mass_file(rootDir, massSliceTag, summaryEdges[i], summaryEdges[i + 1], r))
				{
					collectOk = false;
					break;
				}
				collectedResults.push_back(r);
			}
			if (collectOk)
			{
				std::cout << Form("[INFO] collected %zu saved mass-slice files after targeted rerun; refreshing sliceSummary.",
					collectedResults.size()) << std::endl;
				wroteSummary = writeSliceSummary(collectedResults, 0, -1, true);
			}
			else
			{
				std::cerr << "[WARN] targeted mass-slice rerun updated per-slice ROOT files, but not all slices exist yet; sliceSummary was not refreshed." << std::endl;
				std::cerr << "[WARN] Run the full subrange_mass step once, then targeted reruns can refresh sliceSummary for subrange_ctau." << std::endl;
			}
		}
		else
		{
			wroteSummary = writeSliceSummary(sliceResults, firstSlice, nSlicesToFit, false);
		}

		std::cout << Form("[DONE] mass-slices-only mode wrote %zu slice fit(s); %s final ctau fit.",
			sliceResults.size(), wroteSummary ? "refreshed sliceSummary and skipped" : "skipped") << std::endl;
		return;
	}

	writeSliceSummary(sliceResults, firstSlice, nSlicesToFit, false);
}
