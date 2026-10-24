#ifndef OO_ROBUST_CHI2_H
#define OO_ROBUST_CHI2_H

// Binning-robust goodness of fit for unbinned prefits (2026-10-15).
// RooPlot::chiSquare on a fixed uniform binning breaks down when the bin width is far
// from the scale of the distribution or the window reaches regions without events
// (e.g. NP truth ctau in [0,6] mm with tau ~ 0.4 mm). Here the bin edges are data
// quantiles (equal population, no empty bins) inside the populated range, and the
// expectation in each bin is the normalized PDF integral times the event count.

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TMath.h"
#include "TParameter.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <cstdlib>
#include <string>

struct OORobustChi2 { double chi2 = NAN; int bins = 0, ndf = 0; double lo = NAN, hi = NAN; };

inline OORobustChi2 oo_robust_chi2(RooAbsPdf &pdf, RooDataSet &data, RooRealVar &obs, int npar, int minPerBin = 50, int maxBins = 100)
{
    OORobustChi2 out;
    std::vector<double> v;
    v.reserve(data.numEntries());
    for (int i = 0; i < data.numEntries(); ++i) v.push_back(data.get(i)->getRealValue(obs.GetName()));
    std::sort(v.begin(), v.end());
    const int n = v.size();
    const int nb = std::clamp(n / minPerBin, 5, maxBins);
    if (n < 5 * minPerBin) return out;
    std::vector<double> edges{obs.getMin()};
    for (int k = 1; k < nb; ++k) edges.push_back(v[size_t(double(k) * n / nb)]);
    edges.push_back(obs.getMax());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    const double saveMin = obs.getMin(), saveMax = obs.getMax();
    std::unique_ptr<RooAbsReal> total(pdf.createIntegral(obs, RooFit::NormSet(obs)));
    const double norm = total->getVal();
    double chi2 = 0;
    int used = 0;
    for (size_t k = 0; k + 1 < edges.size(); ++k) {
        const double a = edges[k], b = edges[k + 1];
        const long observed = std::lower_bound(v.begin(), v.end(), b) - std::lower_bound(v.begin(), v.end(), a)
                              + (k + 2 == edges.size() ? std::count(v.begin(), v.end(), b) : 0);
        obs.setRange("oo_robust_bin", a, b);
        std::unique_ptr<RooAbsReal> part(pdf.createIntegral(obs, RooFit::NormSet(obs), RooFit::Range("oo_robust_bin")));
        const double expected = n * part->getVal() / norm;
        if (!(expected > 0)) continue;
        chi2 += (observed - expected) * (observed - expected) / expected;
        ++used;
    }
    obs.setRange(saveMin, saveMax);
    out.chi2 = chi2; out.bins = used; out.ndf = used - npar; out.lo = edges.front(); out.hi = edges.back();
    std::cout << "[ROBUST_CHI2] range=[" << edges.front() << "," << edges.back() << "] n=" << n << '\n';
    std::cout << "[ROBUST_CHI2] equal-population bins=" << used << " chi2/ndf=" << chi2 << "/" << out.ndf << " = "
              << (out.ndf > 0 ? chi2 / out.ndf : NAN) << " p=" << (out.ndf > 0 ? TMath::Prob(chi2, out.ndf) : NAN) << std::endl;
    return out;
}

// Closure check (env OO_ROBUST_TOY=1): the same statistic on a toy drawn from the fitted PDF should give ~1.
inline void oo_robust_chi2_toy(RooAbsPdf &pdf, RooDataSet &data, RooRealVar &obs, int npar)
{
    const char *e = std::getenv("OO_ROBUST_TOY");
    if (!e || std::string(e) != "1") return;
    std::unique_ptr<RooDataSet> toy(pdf.generate(RooArgSet(obs), data.numEntries()));
    std::cout << "[ROBUST_TOY] ";
    oo_robust_chi2(pdf, *toy, obs, npar);
}

inline void oo_write_robust_chi2(const char *path, const OORobustChi2 &r)
{
    TFile f(path, "UPDATE");
    if (f.IsZombie()) return;
    TParameter<double>("prefitChi2NdfRobust", r.ndf > 0 ? r.chi2 / r.ndf : NAN).Write("", TObject::kOverwrite);
    TParameter<int>("prefitChi2RobustBins", r.bins).Write("", TObject::kOverwrite);
    TParameter<int>("prefitChi2RobustNdf", r.ndf).Write("", TObject::kOverwrite);
}

#endif
