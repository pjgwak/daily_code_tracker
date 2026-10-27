// Raw yield comparison with the September OO folder (old acceptance):
// dimuon candidate counts from both DATA inputs and fitted Nsig/Nbkg from both fit2d results.
// Writes drawing/raw_yield_vs_sept/raw_yield.txt (tab-separated).
// run: root -l -b -q drawing/dump_raw_yield_vs_sept.C
#include <TFile.h>
#include <TParameter.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <cstdio>
#include <vector>

namespace
{
const char *SEPT_DIR = "/data/users/pjgwak/work/daily_code_tracker/2026/09/19_OO_2dfit_tuning_Hijing_tuning";
const char *SEPT_DATA = "/data/users/pjgwak/work/daily_code_tracker/2026/04/00_OO_skims_updated/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
const char *NEW_DATA = "inputs/roodataset_roots/RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";

struct FitYield
{
	double nsig = -1, nsigErr = -1, nbkg = -1, nbkgErr = -1, ctLow = 0, ctHigh = 0;
};

FitYield readFit(const TString &dir, double yLow, double yHigh, double ptLow, double ptHigh)
{
	FitYield out;
	TString tag = Form("y%.2f_%.2f", yLow, yHigh);
	TString path = Form("%s/roots/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root", dir.Data(), tag.Data(), tag.Data(), ptLow, ptHigh);
	TFile f(path);
	auto *r = f.IsOpen() ? dynamic_cast<RooFitResult *>(f.Get("modelResult")) : nullptr;
	if (!r)
	{
		std::printf("WARNING: no fit result %s\n", path.Data());
		return out;
	}
	if (auto *v = dynamic_cast<RooRealVar *>(r->floatParsFinal().find("Nsig"))) out.nsig = v->getVal(), out.nsigErr = v->getError();
	if (auto *v = dynamic_cast<RooRealVar *>(r->floatParsFinal().find("Nbkg"))) out.nbkg = v->getVal(), out.nbkgErr = v->getError();
	if (auto *p = dynamic_cast<TParameter<double> *>(f.Get("appliedCtauLow"))) out.ctLow = p->GetVal();
	if (auto *p = dynamic_cast<TParameter<double> *>(f.Get("appliedCtauHigh"))) out.ctHigh = p->GetVal();
	return out;
}

double count(RooDataSet *ds, const TString &cut)
{
	std::unique_ptr<RooAbsData> sel(ds->reduce(cut));
	return sel->sumEntries();
}
} // namespace

void dump_raw_yield_vs_sept()
{
	TFile fOld(SEPT_DATA), fNew(NEW_DATA);
	auto *dOld = static_cast<RooDataSet *>(fOld.Get("dataset"));
	auto *dNew = static_cast<RooDataSet *>(fNew.Get("dataset"));

	struct Bin { const char *region; double yLow, yHigh, ptLow, ptHigh; };
	std::vector<Bin> bins;
	for (auto [lo, hi] : std::vector<std::pair<double, double>>{{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}, {10, 12}, {12, 14}, {14, 20}})
		bins.push_back({"fwd", 1.6, 2.4, lo, hi});
	for (auto [lo, hi] : std::vector<std::pair<double, double>>{{7, 8}, {8, 9}, {9, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 20}, {20, 30}})
		bins.push_back({"mid", 0.0, 1.6, lo, hi});

	gSystem->mkdir("drawing/raw_yield_vs_sept", true);
	FILE *out = std::fopen("drawing/raw_yield_vs_sept/raw_yield.txt", "w");
	std::fprintf(out, "region\tptLow\tptHigh\tcandOld\tcandNew\tcandOldSeptWin\tcandNewSeptWin\t"
	                  "nsigOld\tnsigOldErr\tnbkgOld\tctLowOld\tctHighOld\tnsigNew\tnsigNewErr\tnbkgNew\tctLowNew\tctHighNew\n");
	for (const auto &b : bins)
	{
		TString base = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)", b.ptLow, b.ptHigh, b.yLow, b.yHigh);
		FitYield o = readFit(SEPT_DIR, b.yLow, b.yHigh, b.ptLow, b.ptHigh);
		FitYield n = readFit(".", b.yLow, b.yHigh, b.ptLow, b.ptHigh);
		TString win = base + Form(" && (ctau >= %g && ctau <= %g)", o.ctLow, o.ctHigh);
		std::fprintf(out, "%s\t%g\t%g\t%.0f\t%.0f\t%.0f\t%.0f\t%.3f\t%.3f\t%.3f\t%g\t%g\t%.3f\t%.3f\t%.3f\t%g\t%g\n",
		             b.region, b.ptLow, b.ptHigh, count(dOld, base), count(dNew, base), count(dOld, win), count(dNew, win),
		             o.nsig, o.nsigErr, o.nbkg, o.ctLow, o.ctHigh, n.nsig, n.nsigErr, n.nbkg, n.ctLow, n.ctHigh);
	}
	std::fclose(out);
	std::printf("wrote drawing/raw_yield_vs_sept/raw_yield.txt\n");
}
