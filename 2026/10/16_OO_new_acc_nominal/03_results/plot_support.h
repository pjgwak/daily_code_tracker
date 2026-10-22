#ifndef OO_DRAWING_LATEST_FIT_RESULTS_H
#define OO_DRAWING_LATEST_FIT_RESULTS_H

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace oo_drawing {
struct FitValues { double b, bError, yield, yieldError; };
inline FitValues readFit(const TString &path) {
    TFile file(path, "READ");
    auto *result=file.Get<RooFitResult>("modelResult");
    auto *backend=file.Get<TNamed>("configuredMinimizer");
    auto *migrad=file.Get<TParameter<int>>("minimizeStatus");
    auto *hesse=file.Get<TParameter<int>>("hesseStatus");
    if(file.IsZombie() || !result || result->status()!=0 || result->covQual()!=3 ||
       !std::isfinite(result->edm()) || !backend || std::string(backend->GetTitle())!="Minuit2" ||
       !migrad || migrad->GetVal()!=0 || !hesse || hesse->GetVal()!=0)
        throw std::runtime_error("Invalid or incompatible final fit: "+std::string(path.Data()));
    auto *b=dynamic_cast<const RooRealVar*>(result->floatParsFinal().find("bFraction"));
    auto *n=dynamic_cast<const RooRealVar*>(result->floatParsFinal().find("Nsig"));
    if(!b || !n || !std::isfinite(b->getVal()) || !std::isfinite(b->getError()) ||
       !std::isfinite(n->getVal()) || !std::isfinite(n->getError()) || b->getVal()<0 || b->getVal()>1 ||
       b->getError()<=0 || n->getVal()<=0 || n->getError()<=0)
        throw std::runtime_error("Missing or invalid B/Nsig: "+std::string(path.Data()));
    return {b->getVal(),b->getError(),n->getVal(),n->getError()};
}
}
namespace nominal_bfraction {
constexpr double targetPercent=2.5, forwardAllowedPercent=3., preferredPercent=2., bestSigma=.138, maginotSigma=.414, lastSigma=1.;
struct Bin { double low, high, value, error; bool mid; };
inline double differencePercent(double fitted, const Bin &bin) {
    return 100.0 * (fitted/bin.value - 1.0);
}
inline double differenceSigma(double fitted, const Bin &bin) {
    return std::abs(fitted-bin.value)/bin.error;
}
inline double targetForBin(const Bin &bin) {
    return 100.0 * bestSigma * bin.error / bin.value;
}
// The last tier is explicit: validation prints both percent and sigma criteria.
inline double allowedTargetForBin(const Bin &bin) {
    return std::max(bin.mid ? targetPercent : forwardAllowedPercent, 100.0 * lastSigma * bin.error / bin.value);
}
inline const char *tier(double fitted, const Bin &bin) {
    const double sigma = differenceSigma(fitted, bin);
    const double percent = std::abs(differencePercent(fitted, bin));
    if (!std::isfinite(sigma) || !std::isfinite(percent)) return "OUTSIDE";
    if (sigma <= bestSigma) return "BEST";
    if (sigma <= maginotSigma) return "MAGINOT";
    if (percent < preferredPercent) return "PERCENT_2";
    if (percent <= targetPercent) return "PERCENT_2P5";
    if (!bin.mid && percent <= forwardAllowedPercent) return "PERCENT_3";
    if (sigma <= lastSigma) return "LAST_1SIGMA";
    return "OUTSIDE";
}
}
#endif
