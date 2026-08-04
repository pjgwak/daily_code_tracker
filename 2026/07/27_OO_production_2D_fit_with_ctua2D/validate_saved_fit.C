// 2026-07-31 11:37 KST
// Fail the shell job when a saved per-bin RooFitResult is not usable downstream.

#include "RooFitResult.h"
#include "TFile.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>

namespace
{
bool allHistoryCodesZero(const RooFitResult &result)
{
	for (UInt_t i = 0; i < result.numStatusHistory(); ++i)
		if (result.statusCodeHistory(i) != 0)
			return false;
	return true;
}

template <typename T>
bool readParameter(TFile &file, const char *name, T &value)
{
	auto *parameter = dynamic_cast<TParameter<T> *>(file.Get(name));
	if (!parameter)
		return false;
	value = parameter->GetVal();
	return true;
}
} // namespace

void validate_saved_fit(const char *stage, double ptLow, double ptHigh,
	double yLow, double yHigh, double edmMax = 1e-3,
	double maxChi2Ndf = 2.0)
{
	const TString tag = Form("y%.2f_%.2f_pt%.2f_%.2f", yLow, yHigh, ptLow, ptHigh);
	const TString binDir = Form("roots/y%.2f_%.2f", yLow, yHigh);
	TString fileName;
	TString resultKey;

	if (TString(stage) == "mc_mass")
	{
		fileName = Form("%s/mc_mass/mc_mass_model_%s.root", binDir.Data(), tag.Data());
		resultKey = "fit_result";
	}
	else if (TString(stage) == "mass")
	{
		fileName = Form("%s/mass/mass_model_%s.root", binDir.Data(), tag.Data());
		resultKey = "fit_result";
	}
	else if (TString(stage) == "ctau_pr")
	{
		fileName = Form("%s/ctau_pr/ctau_resolution_%s.root", binDir.Data(), tag.Data());
		resultKey = "timeResult";
	}
	else if (TString(stage) == "ctau_np")
	{
		fileName = Form("%s/ctau_np/ctau_np_model_%s.root", binDir.Data(), tag.Data());
		resultKey = "timeResult";
	}
	else if (TString(stage) == "ctau_bkg")
	{
		fileName = Form("%s/ctau_bkg/ctau_bkg_fitresult_%s.root", binDir.Data(), tag.Data());
		resultKey = "timeResult";
	}
	else if (TString(stage) == "fit2d")
	{
		fileName = Form("%s/fit2d/fit2d_result_%s.root", binDir.Data(), tag.Data());
		resultKey = "modelResult";
	}
	else
	{
		std::cerr << "[VALIDATE] unsupported stage: " << stage << std::endl;
		gSystem->Exit(2);
		return;
	}

	TFile file(fileName, "READ");
	auto *result = file.IsZombie() ? nullptr
		: dynamic_cast<RooFitResult *>(file.Get(resultKey));
	bool pass = result && result->status() == 0 && allHistoryCodesZero(*result)
		&& result->covQual() == 3 && std::isfinite(result->edm())
		&& result->edm() <= edmMax;

	if (result)
	{
		std::cout << "[VALIDATE] " << stage << ' ' << tag
		          << " status=" << result->status()
		          << " covQual=" << result->covQual()
		          << " edm=" << result->edm();
	}
	else
		std::cout << "[VALIDATE] " << stage << ' ' << tag << " missing RooFitResult";

	if (pass && TString(stage) == "fit2d")
	{
		// 2026-07-31 22:06 KST: require and cross-check projection-specific NDF schema version 2.
		double massChi2 = 0.0;
		double timeChi2 = 0.0;
		int projectionNdfVersion = 0;
		int totalFloatNpar = 0;
		int massPullBins = 0;
		int timePullBins = 0;
		int massNpar = 0;
		int timeNpar = 0;
		int massNdf = 0;
		int timeNdf = 0;
		const bool hasQuality = readParameter(file, "projectionNdfVersion", projectionNdfVersion)
			&& readParameter(file, "totalFloatNpar", totalFloatNpar)
			&& readParameter(file, "massPullBins", massPullBins)
			&& readParameter(file, "timePullBins", timePullBins)
			&& readParameter(file, "massNpar", massNpar)
			&& readParameter(file, "timeNpar", timeNpar)
			&& readParameter(file, "massChi2", massChi2)
			&& readParameter(file, "massNdf", massNdf)
			&& readParameter(file, "timeChi2", timeChi2)
			&& readParameter(file, "timeNdf", timeNdf)
			&& projectionNdfVersion == 2
			&& totalFloatNpar >= massNpar && totalFloatNpar >= timeNpar
			&& massNpar >= 0 && timeNpar >= 0
			&& massNdf == massPullBins - massNpar
			&& timeNdf == timePullBins - timeNpar
			&& massNdf > 0 && timeNdf > 0;
		const double massChi2Ndf = hasQuality ? massChi2 / massNdf : -1.0;
		const double timeChi2Ndf = hasQuality ? timeChi2 / timeNdf : -1.0;
		pass = hasQuality && std::isfinite(massChi2Ndf)
			&& std::isfinite(timeChi2Ndf)
			&& massChi2Ndf <= maxChi2Ndf && timeChi2Ndf <= maxChi2Ndf;
		std::cout << " massChi2Ndf=" << massChi2Ndf
		          << " timeChi2Ndf=" << timeChi2Ndf
		          << " massNpar=" << massNpar
		          << " timeNpar=" << timeNpar;
	}

	std::cout << (pass ? " PASS" : " FAIL") << std::endl;
	gSystem->Exit(pass ? 0 : 2);
}
