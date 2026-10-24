#include "TNamed.h"
#include "RooFitResult.h"
#include "nominal_bfraction.h"
#include "RooRealVar.h"
#include <algorithm>
#include "TFile.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <limits>
#include <iostream>

namespace
{
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
												double maxChi2Ndf = 2.0, double agreementPercent = 0.0, const char *outputDir = "")
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

	// -1 selects the current per-bin target without duplicating nominal errors in shell.
	const bool enforceAgreement = agreementPercent != 0.0;
	const bool automaticAgreement = agreementPercent == -1.0;
	if (automaticAgreement)
		agreementPercent = TString(stage) == "fit2d"
													 ? nominal_bfraction::allowedTargetForBin(nominal_bfraction::find(ptLow, ptHigh, yLow, yHigh))
													 : 0.0;
	if (!std::isfinite(agreementPercent) || agreementPercent < 0)
		throw std::runtime_error("Agreement percentage must be nonnegative, or -1 for the per-bin target");
	if (outputDir[0])
		fileName = TString(outputDir) + "/" + fileName;
	TFile file(fileName, "READ");
	auto *result = file.IsZombie() ? nullptr
																 : dynamic_cast<RooFitResult *>(file.Get(resultKey));

	// Judge the saved final result; finite elevated EDM is reported as a warning.
	bool pass = result && result->status() == 0 && result->covQual() == 3 && std::isfinite(result->edm()) && result->edm() >= 0 && std::isfinite(result->minNll());
	if (result && result->edm() > edmMax)
		std::cout << "[VALIDATE] elevated EDM warning: " << result->edm() << std::endl;

	if (result)
	{
		std::cout << "[VALIDATE] " << stage << ' ' << tag
							<< " status=" << result->status()
							<< " covQual=" << result->covQual()
							<< " edm=" << result->edm();
	}
	else
		std::cout << "[VALIDATE] " << stage << ' ' << tag << " missing RooFitResult";


    if (TString(stage) != "fit2d")
    {
        int version=0, bins=0, npar=0, ndf=0;
        double chi2Ndf=-1;
        const bool quality = readParameter(file, "prefitQualityVersion", version) && version == 1
            && readParameter(file, "prefitChi2Bins", bins)
            && readParameter(file, "prefitChi2Npar", npar)
            && readParameter(file, "prefitChi2Ndf", ndf)
            && readParameter(file, "prefitChi2NdfValue", chi2Ndf)
            && result && npar == result->floatParsFinal().getSize() && npar >= 0 && ndf == bins-npar && ndf > 0
            && std::isfinite(chi2Ndf) && chi2Ndf >= 0 && chi2Ndf < maxChi2Ndf;
        pass = pass && quality;
        std::cout << " prefitChi2Ndf=" << chi2Ndf << " ndf=" << ndf
                  << " prefitQuality=" << (quality ? "PASS" : "FAIL_OR_MISSING");
    }

	if (result && TString(stage) == "fit2d")
	{
		auto *backend = file.Get<TNamed>("configuredMinimizer");
		const bool commonMinimizer = backend && std::string(backend->GetTitle()) == "Minuit2";
		pass = pass && commonMinimizer;
		std::cout << " minimizer=" << (backend ? backend->GetTitle() : "MISSING")
		          << " commonMinimizer=" << (commonMinimizer ? "PASS" : "FAIL");
		// 2026-09-10 KST: Validate the direct ROOT result without a Python dispatcher.
		auto *fraction = dynamic_cast<const RooRealVar *>(result->floatParsFinal().find("bFraction"));
		// KST - Never accept a result produced with the prohibited custom convolution.
		int convolution = -1;
		pass = pass && readParameter(file, "stableGaussianConvolution", convolution) && convolution == 0;
		int hesse = -1, minimize = -1;
		double savedB = 0., savedError = 0.;
		const bool hasMinimize = readParameter(file, "minimizeStatus", minimize);
		const bool hasHesse = readParameter(file, "hesseStatus", hesse);
		bool validFraction = fraction && std::isfinite(fraction->getVal()) && std::isfinite(fraction->getError()) && fraction->getError() > 0 && std::min(fraction->getVal() - fraction->getMin(), fraction->getMax() - fraction->getVal()) > std::max(1e-8, 0.01 * std::abs(fraction->getError())) && readParameter(file, "bFraction", savedB) && readParameter(file, "bFractionErr", savedError) && std::abs(savedB - fraction->getVal()) < 1e-10 && std::abs(savedError - fraction->getError()) < 1e-8;
		pass = pass && validFraction && hasHesse && hesse == 0 && (!hasMinimize || minimize == 0);
		for (const char *name : {"rawFinalNll", "unoffsetFullNll", "unoffsetNll"})
		{
			double value = 0.;
			if (readParameter(file, name, value))
				pass = pass && std::isfinite(value);
		}
		std::cout << " hesse=" << hesse << " minimize=" << minimize
							<< " bFree=" << (fraction != nullptr);

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
		const bool hasQuality = readParameter(file, "projectionNdfVersion", projectionNdfVersion) && readParameter(file, "totalFloatNpar", totalFloatNpar) && readParameter(file, "massPullBins", massPullBins) && readParameter(file, "timePullBins", timePullBins) && readParameter(file, "massNpar", massNpar) && readParameter(file, "timeNpar", timeNpar) && readParameter(file, "massChi2", massChi2) && readParameter(file, "massNdf", massNdf) && readParameter(file, "timeChi2", timeChi2) && readParameter(file, "timeNdf", timeNdf) && (projectionNdfVersion == 2 || projectionNdfVersion == 3) && totalFloatNpar >= massNpar && totalFloatNpar >= timeNpar && massNpar >= 0 && timeNpar >= 0 && massNdf == massPullBins - massNpar && timeNdf == timePullBins - timeNpar && massNdf > 0 && timeNdf > 0;
		const double massChi2Ndf = hasQuality ? massChi2 / massNdf : -1.0;
		const double timeChi2Ndf = hasQuality ? timeChi2 / timeNdf : -1.0;
		pass = pass && hasQuality && std::isfinite(massChi2Ndf) && std::isfinite(timeChi2Ndf) && massChi2Ndf >= 0 && timeChi2Ndf >= 0 && massChi2Ndf <= maxChi2Ndf && timeChi2Ndf <= maxChi2Ndf;
		std::cout << " massChi2Ndf=" << massChi2Ndf
							<< " timeChi2Ndf=" << timeChi2Ndf
							<< " massNpar=" << massNpar
							<< " timeNpar=" << timeNpar;
	}

	if (TString(stage) == "fit2d")
	{
		const auto &reference = nominal_bfraction::find(ptLow, ptHigh, yLow, yHigh);
		double fitted = std::numeric_limits<double>::quiet_NaN();
		const double reportTarget = agreementPercent > 0 ? agreementPercent : nominal_bfraction::allowedTargetForBin(reference);
		const bool agreement = readParameter(file, "bFraction", fitted) && nominal_bfraction::agrees(fitted, reference, reportTarget);
		std::cout << " quality=" << (pass ? "PASS" : "FAIL")
							<< " difference_percent=" << nominal_bfraction::differencePercent(fitted, reference)
							<< " target_percent=" << reportTarget
							<< " agreement_policy=" << (!enforceAgreement ? "diagnostic_only" : automaticAgreement ? "allowed_default" : "explicit_required")
							<< " agreement=" << (agreement ? "PASS" : "FAIL");
        const double sigma = nominal_bfraction::differenceSigma(fitted, reference);
        const double percent = std::abs(nominal_bfraction::differencePercent(fitted, reference));
        std::cout << " syst_sigma_nominal=" << sigma
                  << " best_0p138sigma=" << (sigma <= nominal_bfraction::bestSigma ? "PASS" : "FAIL")
                  << " maginot_0p414sigma=" << (sigma <= nominal_bfraction::maginotSigma ? "PASS" : "FAIL")
                  << " percent_lt2=" << (percent < 2.0 ? "PASS" : "FAIL")
                  << " percent_le2p5=" << (percent <= 2.5 ? "PASS" : "FAIL")
                  << " forward_percent_le3=" << (reference.mid ? "NA" : (percent <= nominal_bfraction::forwardAllowedPercent ? "PASS" : "FAIL"))
                  << " last_1sigma=" << (sigma <= 1.0 ? "PASS" : "FAIL")
                  << " acceptance_tier=" << nominal_bfraction::tier(fitted, reference);

		pass = pass && (!enforceAgreement || agreement);
	}
	std::cout << (pass ? " PASS" : " FAIL") << std::endl;
	gSystem->Exit(pass ? 0 : 2);
}