#ifndef SAVED_FIT_HELPERS_H
#define SAVED_FIT_HELPERS_H

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
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

#endif
