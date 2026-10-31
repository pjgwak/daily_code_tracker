#ifndef OO_PREFIT_JITTER_H
#define OO_PREFIT_JITTER_H

// 24 study hook: OO_PREFIT_JITTER=<seed> moves every floating prefit parameter to a random start
// before the first fit (initial-value study). OO_PREFIT_JITTER_FRAC (default 0.3) sets the size:
// shift = frac * U(-1,1) * min(range, max(|value|, 0.1 * range)), kept inside the range.
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "TRandom3.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>

inline void oo_prefit_jitter(RooAbsPdf &pdf, const RooArgSet &observables) {
    const char *seedText = gSystem->Getenv("OO_PREFIT_JITTER");
    if (!seedText) return;
    const double frac = gSystem->Getenv("OO_PREFIT_JITTER_FRAC") ? std::atof(gSystem->Getenv("OO_PREFIT_JITTER_FRAC")) : 0.3;
    TRandom3 rng(std::atoi(seedText));
    std::unique_ptr<RooArgSet> params(pdf.getParameters(observables));
    for (auto *a : *params) {
        auto *v = dynamic_cast<RooRealVar *>(a);
        if (!v || v->isConstant() || !v->hasMin() || !v->hasMax()) continue;
        const double lo = v->getMin(), hi = v->getMax(), range = hi - lo, before = v->getVal();
        const double step = frac * rng.Uniform(-1, 1) * std::min(range, std::max(std::abs(before), 0.1 * range));
        v->setVal(std::clamp(before + step, lo + 1e-6 * range, hi - 1e-6 * range));
        std::printf("[PREFIT_JITTER] seed=%s %s %.6g -> %.6g\n", seedText, v->GetName(), before, v->getVal());
    }
}

#endif
