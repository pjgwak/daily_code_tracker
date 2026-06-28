#include "RooAddModel.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
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
#include "RooLinkedList.h"
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


SliceResult merge_ctau_slice_span(const std::vector<SliceResult> &slices, size_t first, size_t last)
{
	SliceResult merged = slices[first];
	merged.ctLow = slices[first].ctLow;
	merged.ctHigh = slices[last].ctHigh;
	merged.entries = 0;
	merged.nsig = 0.0;
	merged.nsigErr = 0.0;
	merged.nbkg = 0.0;
	merged.nbkgErr = 0.0;
	merged.signalWindowEntries = 0;
	merged.shapeDrivenSignal = 0;
	merged.usedInCtauFit = 1;
	merged.status = 0;
	merged.covQual = 3;
	merged.fitAttempt = 0;
	double err2 = 0.0;
	double bkgErr2 = 0.0;
	for (size_t i = first; i <= last; ++i)
	{
		const auto &r = slices[i];
		merged.entries += r.entries;
		merged.nsig += r.nsig;
		err2 += r.nsigErr * r.nsigErr;
		merged.nbkg += r.nbkg;
		bkgErr2 += r.nbkgErr * r.nbkgErr;
		merged.signalWindowEntries += r.signalWindowEntries;
		merged.shapeDrivenSignal = merged.shapeDrivenSignal || r.shapeDrivenSignal;
		if (r.status != 0)
			merged.status = r.status;
		merged.covQual = std::min(merged.covQual, r.covQual);
		merged.fitAttempt = std::max(merged.fitAttempt, r.fitAttempt);
	}
	merged.nsigErr = std::sqrt(err2);
	merged.nbkgErr = std::sqrt(bkgErr2);
	merged.nsigSignificance = merged.nsigErr > 0.0 ? merged.nsig / merged.nsigErr : 0.0;
	return merged;
}

void merge_ctau_slice_span_in_place(std::vector<SliceResult> &slices, size_t first, size_t last)
{
	if (slices.empty() || first >= slices.size() || last >= slices.size() || last < first)
		return;
	SliceResult merged = merge_ctau_slice_span(slices, first, last);
	slices.erase(slices.begin() + first, slices.begin() + last + 1);
	slices.insert(slices.begin() + first, merged);
}

size_t enforce_non_decreasing_tail_widths(std::vector<SliceResult> &slices)
{
	const double eps = 1e-9;
	auto width = [](const SliceResult &r) { return r.ctHigh - r.ctLow; };
	size_t nMerged = 0;

	while (true)
	{
		size_t first = slices.size();
		size_t last = slices.size();
		for (size_t i = 0; i < slices.size(); ++i)
		{
			if (slices[i].ctLow >= -eps)
			{
				if (first == slices.size())
					first = i;
				last = i;
			}
		}
		if (first == slices.size() || last <= first)
			break;
		bool changed = false;
		for (size_t i = first + 1; i <= last; ++i)
		{
			if (width(slices[i]) + eps >= width(slices[i - 1]))
				continue;
			if (i < last)
				merge_ctau_slice_span_in_place(slices, i, i + 1);
			else
				merge_ctau_slice_span_in_place(slices, i - 1, i);
			++nMerged;
			changed = true;
			break;
		}
		if (!changed)
			break;
	}

	while (true)
	{
		size_t first = slices.size();
		size_t last = slices.size();
		for (size_t i = 0; i < slices.size(); ++i)
		{
			if (slices[i].ctHigh <= eps)
			{
				if (first == slices.size())
					first = i;
				last = i;
			}
		}
		if (first == slices.size() || last <= first)
			break;
		bool changed = false;
		for (size_t i = last; i-- > first; )
		{
			if (width(slices[i]) + eps >= width(slices[i + 1]))
				continue;
			if (i > first)
				merge_ctau_slice_span_in_place(slices, i - 1, i);
			else
				merge_ctau_slice_span_in_place(slices, i, i + 1);
			++nMerged;
			changed = true;
			break;
		}
		if (!changed)
			break;
	}

	return nMerged;
}

}

void subrange_ctau(float ptLow = 3.0, float ptHigh = 6.5, float yLow = 1.6, float yHigh = 2.4,
	double ctMin = -8.0, double ctMax = 10.0, double ctStep = 2.0,
	int maxSlices = -1, bool drawMassSlices = true, bool trimSparseFwd = true,
	bool applyQuantileCut = true, bool reuseMassSlices = true, double centerVetoHalfWidth = 0.0,
	int firstSlice = 0, int nSlicesToFit = -1, bool parallelMassSlices = false, int maxParallelJobs = 4,
	bool massSlicesOnly = false)
{
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
	gStyle->SetOptTitle(0);
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	(void)ctMin;
	(void)ctMax;
	(void)maxSlices;
	(void)drawMassSlices;

	const TString yTag = y_tag(yLow, yHigh);
	const TString tag = bin_tag(ptLow, ptHigh, yLow, yHigh);
	const TString baseOutputTag = tag + (trimSparseFwd ? "" : "_fullrange") + (applyQuantileCut ? "" : "_noQuantile");
	TString outputTag = baseOutputTag;
	if (firstSlice > 0 || nSlicesToFit > 0)
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
	TString inputName = outputName;

	std::unique_ptr<TFile> prevOut(TFile::Open(inputName, "READ"));
	if (!prevOut || prevOut->IsZombie())
	{
		if (prevOut)
			prevOut.reset();
		TString ctStepTag = TString::Format("_ctStep%.2f", ctStep);
		ctStepTag.ReplaceAll(".", "p");
		const TString legacyOutputTag = baseOutputTag + ctStepTag + (outputTag(baseOutputTag.Length(), outputTag.Length() - baseOutputTag.Length()));
		const TString legacyInputName = TString::Format("%s/subrange_result_%s.root", rootDir.Data(), legacyOutputTag.Data());
		prevOut.reset(TFile::Open(legacyInputName, "READ"));
		if (prevOut && !prevOut->IsZombie())
		{
			inputName = legacyInputName;
			std::cout << "[INFO] using legacy ctStep-tagged mass-slice summary: " << inputName << std::endl;
		}
	}
	if (!prevOut || prevOut->IsZombie())
	{
		std::cerr << "ERROR: cannot run subrange ctau fit; mass-slice summary not found: " << outputName << std::endl;
		return;
	}
	auto *summary = dynamic_cast<TTree *>(prevOut->Get("sliceSummary"));
	if (!summary)
	{
		std::cerr << "ERROR: sliceSummary not found in mass-slice summary: " << outputName << std::endl;
		return;
	}

	double minEdgeWidth = read_double(prevOut.get(), "adaptiveMinBinWidth", 0.0);
	double maxEdgeWidth = read_double(prevOut.get(), "adaptiveMaxBinWidth", 0.0);
	double maxNeighborWidthRatio = read_double(prevOut.get(), "adaptiveMaxNeighborWidthRatio", 1.0);
	std::vector<SliceResult> sliceResults;
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
	prevOut.reset();
	std::cout << Form("[INFO] loaded %zu mass subrange slices from %s", sliceResults.size(), inputName.Data()) << std::endl;
	if (sliceResults.empty())
	{
		std::cerr << "ERROR: no mass subrange slices available for ctau fit." << std::endl;
		return;
	}
	// Automatic final-ctau range trimming: find the first and last reliable
	// mass-signal slices and drop empty/noisy tails before fitting the ctau model.
	// Manual overrides below can further replace this automatic range per bin.
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

	// Optional human override of the final ctau fit window. This does not refit
	// mass slices; it only changes which already-fitted slices are used to build
	// hSignalYield and the RooDataHist for the final ctau fit.
	int manualCtauFitRangeApplied = 0;
	double manualCtauFitLow = 0.0;
	double manualCtauFitHigh = 0.0;
	double actualManualCtauFitLow = 0.0;
	double actualManualCtauFitHigh = 0.0;
	if (const auto *manualRange = find_manual_ctau_fit_range(ptLow, ptHigh, yLow, yHigh))
	{
		manualCtauFitLow = manualRange->ctLow;
		manualCtauFitHigh = manualRange->ctHigh;
		std::vector<SliceResult> manualSlices;
		manualSlices.reserve(sliceResults.size());
		for (auto &r : sliceResults)
		{
			// Keep any adaptive slice that overlaps the requested manual interval. This
			// preserves the native bin edges instead of inventing a partial-bin yield.
			const bool use = r.ctHigh > manualRange->ctLow + 1e-9 && r.ctLow < manualRange->ctHigh - 1e-9;
			r.usedInCtauFit = use ? 1 : 0;
			if (use)
				manualSlices.push_back(r);
		}
		if (manualSlices.empty())
		{
			std::cerr << Form("ERROR: enabled manual ctau fit range [%.6g, %.6g] has no overlapping fitted slices for this bin.",
				manualRange->ctLow, manualRange->ctHigh) << std::endl;
			return;
		}
		ctauFitSlices.swap(manualSlices);
		manualCtauFitRangeApplied = 1;
		actualManualCtauFitLow = ctauFitSlices.front().ctLow;
		actualManualCtauFitHigh = ctauFitSlices.back().ctHigh;
		std::cout << Form("[INFO] manual ctau fit range override applied: requested [%.6g, %.6g], using binned range [%.6g, %.6g] with %zu slices%s%s",
			manualRange->ctLow, manualRange->ctHigh, actualManualCtauFitLow, actualManualCtauFitHigh,
			ctauFitSlices.size(), manualRange->note && manualRange->note[0] ? "; note: " : "",
			manualRange->note && manualRange->note[0] ? manualRange->note : "") << std::endl;
	}

	// Standard RooPlot::pullHist cannot assign a meaningful pull to zero-error
	// adaptive bins. Keep all mass-slice outputs, but trim zero-error edge bins
	// from the final ctau histogram so the pull panel reflects fitted slices only.
	auto hasCtauFitError = [](const SliceResult &r) {
		return r.entries > 0 && r.nsigErr > 0.0 && std::isfinite(r.nsigErr);
	};
	const size_t beforeEdgeTrim = ctauFitSlices.size();
	while (!ctauFitSlices.empty() && !hasCtauFitError(ctauFitSlices.front()))
		ctauFitSlices.erase(ctauFitSlices.begin());
	while (!ctauFitSlices.empty() && !hasCtauFitError(ctauFitSlices.back()))
		ctauFitSlices.pop_back();
	if (ctauFitSlices.empty())
	{
		std::cerr << "ERROR: no non-zero-error ctau bins remain after adaptive edge-bin trimming." << std::endl;
		return;
	}
	for (auto &r : sliceResults)
	{
		const bool use = r.ctHigh > ctauFitSlices.front().ctLow + 1e-9 &&
			r.ctLow < ctauFitSlices.back().ctHigh - 1e-9;
		r.usedInCtauFit = use ? 1 : 0;
	}
	if (ctauFitSlices.size() != beforeEdgeTrim)
	{
		std::cout << Form("[INFO] adaptive ctau pull binning: removed %zu zero-error edge bins, using %zu/%zu slices, range [%.6g, %.6g]",
			beforeEdgeTrim - ctauFitSlices.size(), ctauFitSlices.size(), beforeEdgeTrim,
			ctauFitSlices.front().ctLow, ctauFitSlices.back().ctHigh) << std::endl;
	}
	if (manualCtauFitRangeApplied)
	{
		actualManualCtauFitLow = ctauFitSlices.front().ctLow;
		actualManualCtauFitHigh = ctauFitSlices.back().ctHigh;
	}


	// Apply default ctau binning bands after automatic range trimming. Exact-bin
	// overrides from the shared config replace the default bands when needed.
	const ManualCtauBinningOverride *ctauBinningOverride = find_manual_ctau_binning_override(ptLow, ptHigh, yLow, yHigh);
	const bool optimizedCtauBinningApplied = ctauBinningOverride != nullptr || has_default_ctau_binning();
	const char *ctauBinningNote = ctauBinningOverride ? ctauBinningOverride->note : default_ctau_binning_note();
	if (optimizedCtauBinningApplied)
	{
		const size_t beforeMerge = ctauFitSlices.size();
		std::vector<SliceResult> mergedSlices;
		for (size_t i = 0; i < ctauFitSlices.size(); )
		{
			const double start = ctauFitSlices[i].ctLow;
			const double width0 = ctauFitSlices[i].ctHigh - ctauFitSlices[i].ctLow;
			const double center = 0.5 * (ctauFitSlices[i].ctLow + ctauFitSlices[i].ctHigh);
			const double targetWidth = target_ctau_bin_width(ctauBinningOverride, center, width0);
			SliceResult merged = ctauFitSlices[i];
			merged.entries = 0;
			merged.nsig = 0.0;
			merged.nsigErr = 0.0;
			merged.nbkg = 0.0;
			merged.nbkgErr = 0.0;
			merged.signalWindowEntries = 0;
			merged.shapeDrivenSignal = 0;
			merged.status = 0;
			merged.covQual = 3;
			double err2 = 0.0;
			double bkgErr2 = 0.0;
			size_t j = i;
			while (j < ctauFitSlices.size())
			{
				const auto &r = ctauFitSlices[j];
				merged.ctHigh = r.ctHigh;
				merged.entries += r.entries;
				merged.nsig += r.nsig;
				err2 += r.nsigErr * r.nsigErr;
				merged.nbkg += r.nbkg;
				bkgErr2 += r.nbkgErr * r.nbkgErr;
				merged.signalWindowEntries += r.signalWindowEntries;
				merged.shapeDrivenSignal = merged.shapeDrivenSignal || r.shapeDrivenSignal;
				if (r.status != 0)
					merged.status = r.status;
				merged.covQual = std::min(merged.covQual, r.covQual);
				++j;
				if (merged.ctHigh - start >= targetWidth - 1e-9)
					break;
			}
			merged.nsigErr = std::sqrt(err2);
			merged.nbkgErr = std::sqrt(bkgErr2);
			merged.nsigSignificance = merged.nsigErr > 0.0 ? merged.nsig / merged.nsigErr : 0.0;
			merged.usedInCtauFit = 1;
			mergedSlices.push_back(merged);
			i = j;
		}
		ctauFitSlices.swap(mergedSlices);
		const size_t monotonicTailMerges = enforce_non_decreasing_tail_widths(ctauFitSlices);
		if (monotonicTailMerges > 0)
		{
			std::cout << Form("[INFO] monotonic tail binning for %s: applied %zu outward-width merges",
				bin_tag(ptLow, ptHigh, yLow, yHigh).Data(), monotonicTailMerges) << std::endl;
		}
		if (manualCtauFitRangeApplied)
		{
			actualManualCtauFitLow = ctauFitSlices.front().ctLow;
			actualManualCtauFitHigh = ctauFitSlices.back().ctHigh;
		}
		std::cout << Form("[INFO] optimized final ctau binning for %s: merged %zu -> %zu bins; note: %s",
			bin_tag(ptLow, ptHigh, yLow, yHigh).Data(), beforeMerge, ctauFitSlices.size(), ctauBinningNote) << std::endl;
	}

	// ------------------------------------------------------------------
	// build signal-yield histogram for the final ctau fit
	// ------------------------------------------------------------------
	// The mass fits provide one signal yield per ctau slice. Convert those yields
	// into a RooDataHist with the exact adaptive bin edges so fitTo can integrate
	// the prompt/nonprompt ctau PDFs over each nonuniform bin.
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

	// ------------------------------------------------------------------
	// final ctau fit range controls
	// ------------------------------------------------------------------
	// centerVetoHalfWidth is a diagnostic switch: it masks the central prompt peak
	// from the final ctau fit without changing the mass-slice fits themselves.
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

	// ------------------------------------------------------------------
	// ctau model setup from MC-constrained references
	// ------------------------------------------------------------------
	// Prompt resolution and nonprompt lifetime seeds come from the dedicated
	// ctau_pr/ctau_np steps. If those files are absent, the fallback model keeps
	// the macro usable for diagnostics, but production should use the references.
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

	// ------------------------------------------------------------------
	// external constraints for final ctau fit
	// ------------------------------------------------------------------
	// The subrange method fits yields, not individual events. Loose Gaussian
	// constraints keep the MC-derived prompt-resolution hierarchy and nonprompt
	// lifetime scale from drifting into unphysical solutions in weak bins.
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

	// ------------------------------------------------------------------
	// final ctau fit model and yield parameters
	// ------------------------------------------------------------------
	// Nprompt and Nnonprompt float against the binned signal-yield histogram. The
	// prompt component is the resolution model; the nonprompt component is a sum
	// of single-sided decays convolved with the same prompt resolution.
	const double totalYield = std::max(1.0, hSignalYield->Integral());
	RooRealVar Nprompt("Nprompt", "Nprompt", 0.80 * totalYield, 0.0, 2.0 * totalYield);
	RooRealVar Nnonprompt("Nnonprompt", "Nnonprompt", 0.20 * totalYield, 0.0, 2.0 * totalYield);
	RooAddPdf ctauModel("ctauModel", "ctauModel", RooArgList(*promptTimePdf, *signalNpTime), RooArgList(Nprompt, Nnonprompt));

	// ------------------------------------------------------------------
	// run final ctau fit
	// ------------------------------------------------------------------
	// The final observable is a binned mass-fit signal-yield histogram with
	// nonuniform adaptive bins. IntegrateBins makes fitTo compare bin integrals,
	// not point values; SumW2Error(false) keeps the HESSE errors tied to the
	// summed yield scale of the histogram.
	std::unique_ptr<RooFitResult> ctauFit;
	if (ctauFitSlices.size() >= 2 && totalYield > 0.0)
	{
		const bool applyCtauConstraints = useReferenceCtauMC && ctauConstraints.getSize() > 0;
		if (applyCtauConstraints)
		{
			std::cout << Form("[INFO] final ctau fit uses fitTo with IntegrateBins, SumW2Error(false), and %d ExternalConstraints.",
				ctauConstraints.getSize()) << std::endl;
		}
		else
		{
			std::cout << "[INFO] final ctau fit uses fitTo with IntegrateBins and SumW2Error(false)." << std::endl;
		}
		auto runCtauFit = [&](int strategy) -> RooFitResult * {
			if (useCenterVeto)
			{
				if (applyCtauConstraints)
				{
					return ctauModel.fitTo(signalYieldData,
						Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
						ExternalConstraints(ctauConstraints), Range(ctauFitRangeName.Data()),
						Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(strategy));
				}
				return ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					Range(ctauFitRangeName.Data()), Minimizer("Minuit2", "migrad"),
					PrintLevel(-1), Strategy(strategy));
			}
			if (applyCtauConstraints)
			{
				return ctauModel.fitTo(signalYieldData,
					Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
					ExternalConstraints(ctauConstraints), Minimizer("Minuit2", "migrad"),
					PrintLevel(-1), Strategy(strategy));
			}
			return ctauModel.fitTo(signalYieldData,
				Extended(), Save(), IntegrateBins(1e-6), SumW2Error(false), Offset(true),
				Minimizer("Minuit2", "migrad"), PrintLevel(-1), Strategy(strategy));
		};
		ctauFit.reset(runCtauFit(1));
		if (!ctauFit || ctauFit->status() != 0 || ctauFit->covQual() < 2)
		{
			std::cout << Form("[INFO] final ctau fitTo status=%d covQual=%d; retrying strategy 2.",
				ctauFit ? ctauFit->status() : -999, ctauFit ? ctauFit->covQual() : -999) << std::endl;
			ctauFit.reset(runCtauFit(2));
		}
	}
	else
	{
		std::cout << "[INFO] skipping final ctau fit: need at least two non-empty subrange bins." << std::endl;
	}
	if (ctauFit)
		ctauFit->Print();

	// ------------------------------------------------------------------
	// draw final ctau projection and pull
	// ------------------------------------------------------------------
	TCanvas cCtau("c_ctau_subrange", "c_ctau_subrange", 800, 800);
	cCtau.cd();
	TPad *ctPad1 = new TPad("ctPad1", "ctPad1", 0.0, 0.25, 1.0, 1.0);
	ctPad1->SetBottomMargin(0.00001);
	ctPad1->SetTopMargin(0.08);
	ctPad1->SetLogy();
	ctPad1->Draw();
	ctPad1->cd();

	auto *ctFrame = ctau.frame(Title(""));
	ctFrame->SetTitle("");
	signalYieldData.plotOn(ctFrame, Binning(ctauBinning), Name("yieldData"));
	ctauModel.plotOn(ctFrame, Name("ctauTotal"), LineColor(kBlack), LineWidth(2), Precision(1e-5));
	ctauModel.plotOn(ctFrame, Components(*promptTimePdf), Name("ctauPrompt"), LineColor(kRed + 1), LineStyle(kDashed), LineWidth(2), Precision(1e-5));
	ctauModel.plotOn(ctFrame, Components(*signalNpTime), Name("ctauNonPrompt"), LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Precision(1e-5));
	ctFrame->SetTitle("");
	ctFrame->SetFillStyle(4000);
	ctFrame->GetXaxis()->SetTitle("");
	ctFrame->GetXaxis()->SetTitleSize(0);
	ctFrame->GetXaxis()->SetLabelSize(0);
	ctFrame->GetYaxis()->SetTitle("Events");
	ctFrame->GetYaxis()->SetTitleOffset(1.6);
	apply_logy_auto_range(ctFrame, "yieldData");
	ctFrame->Draw("e");
	TLegend ctLeg(0.50, 0.70, 0.72, 0.89);
	ctLeg.SetBorderSize(0);
	ctLeg.SetFillStyle(0);
	ctLeg.SetTextSize(0.03);
	auto findCtObj = [&](const char *n) -> TObject * { return ctFrame ? ctFrame->findObject(n) : nullptr; };
	if (auto *o = findCtObj("yieldData"))
		ctLeg.AddEntry(o, "Data", "lep");
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
		TLatex cms;
		cms.SetNDC();
		cms.SetTextAlign(11);
		cms.SetTextSize(0.040);
		cms.SetTextFont(62);
		cms.DrawLatex(0.21, 0.935, "CMS");

		TLatex extra;
		extra.SetNDC();
		extra.SetTextAlign(11);
		extra.SetTextSize(0.040);
		extra.SetTextFont(52);
		extra.DrawLatex(0.315, 0.935, "Internal");
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
		tp.SetTextSize(0.024);
		tp.SetTextFont(42);
		double xtext = 0.73, y0 = 0.87, dy = -0.045;
		int k = 0;
		tp.DrawLatex(xtext, y0 + dy * k++, Form("f_{B} = %.4g #pm %.3g", bFractionVal, bFractionErr));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{sig} = %.4g", hSignalYield->Integral()));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{PR} = %.4g", Nprompt.getVal()));
		tp.DrawLatex(xtext, y0 + dy * k++, Form("N_{NP} = %.4g", Nnonprompt.getVal()));
	}
	{
		TLatex tx;
		tx.SetNDC();
		tx.SetTextSize(0.03);
		tx.SetTextFont(42);
		int minimize = ctauFit ? ctauFit->status() : -1;
		int hesse = -1;
		if (ctauFit)
		{
			for (UInt_t i = 0, n = ctauFit->numStatusHistory(); i < n; ++i)
			{
				const char *lab = ctauFit->statusLabelHistory(i);
				if (!lab)
					continue;
				if (TString(lab) == "MINIMIZE" || TString(lab) == "MIGRAD")
					minimize = ctauFit->statusCodeHistory(i);
				else if (TString(lab) == "HESSE")
					hesse = ctauFit->statusCodeHistory(i);
			}
		}
		tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", minimize, hesse));
	}

	cCtau.cd();
	TPad *ctPad2 = new TPad("ctPad2", "ctPad2", 0.0, 0.0, 1.0, 0.25);
	ctPad2->SetTopMargin(0.00001);
	ctPad2->SetBottomMargin(0.4);
	ctPad2->Draw();
	ctPad2->cd();

	RooPlot *ctPullPlot = ctau.frame(Title(""));
	ctPullPlot->SetTitle("");
	RooHist *ctPull = ctFrame->pullHist("yieldData", "ctauTotal", true);
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
	ctPullPlot->SetMinimum(-8);
	ctPullPlot->SetMaximum(8);
	ctPullPlot->GetYaxis()->SetNdivisions(505);
	ctPullPlot->GetYaxis()->SetTitleSize(0.12);
	ctPullPlot->GetYaxis()->SetLabelSize(0.10);
	ctPullPlot->GetXaxis()->SetTitleSize(0.15);
	ctPullPlot->GetXaxis()->SetLabelSize(0.10);
	ctPullPlot->SetTitle("");
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

	// ------------------------------------------------------------------
	// write final ROOT output
	// ------------------------------------------------------------------
	// sliceSummary contains every mass slice, including slices that were excluded
	// from the final ctau fit. The usedInCtauFit branch marks the actual fit range.
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
		TParameter<int>("optimizedCtauBinningApplied", optimizedCtauBinningApplied ? 1 : 0).Write();
		TParameter<double>("adaptiveMinBinWidth", minEdgeWidth).Write();
		TParameter<double>("adaptiveMaxBinWidth", maxEdgeWidth).Write();
		TParameter<double>("adaptiveMaxNeighborWidthRatio", maxNeighborWidthRatio).Write();
		TParameter<double>("edgeSignalThreshold", edgeSignalThreshold).Write();
		TParameter<double>("ctauFitMinOverride", kCtauFitMinOverride).Write();
		// Persist run controls and manual-range bookkeeping so downstream checks can
		// tell whether a result came from automatic trimming, ctau-only reuse, a
		// parallel mass-slice run, or a manually forced ctau range.
		TParameter<int>("manualCtauFitRangeApplied", manualCtauFitRangeApplied).Write();
		TParameter<double>("manualCtauFitLow", manualCtauFitLow).Write();
		TParameter<double>("manualCtauFitHigh", manualCtauFitHigh).Write();
		TParameter<double>("actualManualCtauFitLow", actualManualCtauFitLow).Write();
		TParameter<double>("actualManualCtauFitHigh", actualManualCtauFitHigh).Write();
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
		TParameter<int>("reuseMassSlices", 1).Write();
		TParameter<int>("firstSlice", firstSlice).Write();
		TParameter<int>("nSlicesToFit", nSlicesToFit).Write();
		TParameter<int>("parallelMassSlices", parallelMassSlices ? 1 : 0).Write();
		TParameter<int>("maxParallelJobs", maxParallelJobs).Write();
		TParameter<int>("massSlicesOnly", 0).Write();
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
