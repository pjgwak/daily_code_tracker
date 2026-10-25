#ifndef SAVED_FIT_HELPERS_H
#define SAVED_FIT_HELPERS_H

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "TSystem.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "TFile.h"
#include "TIterator.h"
#include "TParameter.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

static std::unique_ptr<RooFitResult> clone_saved_fit_result(TFile *file, const char *name)
{
	if (!file)
		return nullptr;
	auto *result = dynamic_cast<RooFitResult *>(file->Get(name));
	if (!result)
		return nullptr;
	return std::unique_ptr<RooFitResult>(static_cast<RooFitResult *>(result->Clone()));
}

static bool apply_saved_fit_result(RooFitResult *fitResult, RooArgSet *params)
{
	if (!fitResult || !params)
		return false;

	bool applied = false;
	for (auto *obj : *params)
	{
		auto *target = dynamic_cast<RooRealVar *>(obj);
		if (!target)
			continue;

		if (auto *src = dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find(target->GetName())))
		{
			target->setVal(src->getVal());
			const double err = src->getError();
			if (std::isfinite(err) && err >= 0.0)
				target->setError(err);
			target->setConstant(false);
			applied = true;
			continue;
		}
		if (auto *src = dynamic_cast<RooRealVar *>(fitResult->constPars().find(target->GetName())))
		{
			target->setVal(src->getVal());
			target->setConstant(true);
			applied = true;
		}
	}
	return applied;
}

static bool apply_saved_fit_result(RooFitResult *fitResult, RooAbsPdf &pdf, const RooArgSet &observables)
{
	std::unique_ptr<RooArgSet> params(pdf.getParameters(observables));
	return apply_saved_fit_result(fitResult, params.get());
}

static int read_saved_int_param(TFile *file, const char *name, int fallback)
{
	if (!file)
		return fallback;
	auto *param = dynamic_cast<TParameter<int> *>(file->Get(name));
	return param ? param->GetVal() : fallback;
}

static double read_saved_double_param(TFile *file, const char *name, double fallback = std::numeric_limits<double>::quiet_NaN())
{
	if (!file)
		return fallback;
	if (auto *param = dynamic_cast<TParameter<double> *>(file->Get(name)))
		return param->GetVal();
	if (auto *var = dynamic_cast<RooRealVar *>(file->Get(name)))
		return var->getVal();
	if (auto *abs = dynamic_cast<RooAbsReal *>(file->Get(name)))
		return abs->getVal();
	return fallback;
}

static bool load_saved_fit_file(std::unique_ptr<TFile> &file, const TString &path, const char *label)
{
	file.reset(TFile::Open(path, "READ"));
	if (!file || file->IsZombie())
	{
		std::cerr << "ERROR: cannot open saved " << (label ? label : "fit") << " file: " << path << std::endl;
		return false;
	}
	return true;
}

static int saved_prefit_npar(const TString &path, const char *key)
{
    TFile file(path, "READ");
    auto *result = file.Get<RooFitResult>(key);
    if (!result) throw std::runtime_error("Missing saved prefit result for NDF");
    return result->floatParsFinal().getSize();
}

// A failed or unaudited prefit must never feed a fresh final likelihood.
static void require_prefit_quality(TFile &file, const char *resultKey, double maximum = 2.0)
{
    // Diagnostic-only override (study/prefit_variant.py --env=OO_PREFIT_CHI2_MAX=...): never set in production.
    if (const char *m = gSystem->Getenv("OO_PREFIT_CHI2_MAX")) {
        maximum = std::atof(m);
        std::cout << "[PREFIT_QUALITY] WARNING: chi2/ndf gate raised to " << maximum << " for " << file.GetName() << std::endl;
    }
    auto *r = file.Get<RooFitResult>(resultKey);
    const double chi = read_saved_double_param(&file, "prefitChi2NdfValue");
    const int bins = read_saved_int_param(&file, "prefitChi2Bins", -1);
    const int npar = read_saved_int_param(&file, "prefitChi2Npar", -1);
    const int ndf = read_saved_int_param(&file, "prefitChi2Ndf", -1);
    if (!r || r->status() != 0 || r->covQual() != 3 || !std::isfinite(r->edm())
        || r->edm() < 0 || !std::isfinite(r->minNll())
        || read_saved_int_param(&file, "prefitQualityVersion", 0) != 1
        || npar != r->floatParsFinal().getSize() || ndf != bins-npar || ndf <= 0
        || !std::isfinite(chi) || chi < 0 || chi >= maximum)
        throw std::runtime_error(std::string("Invalid or unaudited prefit: ") + file.GetName());
}

// Match RooCurve::chiSquare: nonzero data points inside the plotted curve.
static int prefit_chi2_bins(RooPlot &plot)
{
    auto *hist = plot.getHist("data");
    auto *curve = plot.getCurve("model");
    if (!hist || !curve || curve->GetN() < 2) return 0;
    int count = 0;
    for (int i = 0; i < hist->GetN(); ++i)
        if (hist->GetPointY(i) != 0 && hist->GetPointX(i) >= curve->GetPointX(0)
            && hist->GetPointX(i) <= curve->GetPointX(curve->GetN()-1)) ++count;
    return count;
}

static void record_prefit_quality(const TString &path, double chi2Ndf, int bins, int npar, bool drawOnly)
{
    const int ndf = bins - npar;
    std::cout << std::setprecision(17) << "[PREFIT_QUALITY] file=" << path
              << " chi2Ndf=" << chi2Ndf << " bins=" << bins << " npar=" << npar
              << " ndf=" << ndf << std::endl;
    if (drawOnly) return;
    TFile file(path, "UPDATE");
    if (file.IsZombie()) throw std::runtime_error("Cannot save prefit quality");
    TParameter<int>("prefitQualityVersion", 1).Write("", TObject::kOverwrite);
    TParameter<int>("prefitChi2Bins", bins).Write("", TObject::kOverwrite);
    TParameter<int>("prefitChi2Npar", npar).Write("", TObject::kOverwrite);
    TParameter<int>("prefitChi2Ndf", ndf).Write("", TObject::kOverwrite);
    TParameter<double>("prefitChi2NdfValue", chi2Ndf).Write("", TObject::kOverwrite);
}

#endif
