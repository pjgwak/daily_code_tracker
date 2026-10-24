// Collect B and fit quality for all 20 bins under one results directory into CSV.
#include "../utils/nominal_bfraction.h"
#include "TFile.h"
#include "TParameter.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include <cstdio>
template <class T> static double par(TFile &f, const char *n) { auto *p = f.Get<TParameter<T>>(n); return p ? double(p->GetVal()) : NAN; }
void collect(const char *base, const char *out) {
  FILE *o = fopen(out, "w");
  fprintf(o, "region,ptLow,ptHigh,B,Berr,Bnom,BnomErr,pull,percent,status,covQual,edm,minimize,hesse,massChi2,massNdf,timeChi2,timeNdf,rawNll,massPvalue,timePvalue\n");
  for (int mid = 0; mid < 2; ++mid)
    for (int i = 0; i < (mid ? 8 : 12); ++i) {
      const auto &ref = mid ? nominal_bfraction::midBins()[i] : nominal_bfraction::bins()[i];
      const char *tag = mid ? "y0.00_1.60" : "y1.60_2.40";
      TFile f(Form("%s/roots/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root", base, tag, tag, ref.low, ref.high));
      auto *r = f.IsZombie() ? nullptr : f.Get<RooFitResult>("modelResult");
      double B = NAN, Be = NAN;
      if (r) if (auto *b = dynamic_cast<RooRealVar *>(r->floatParsFinal().find("bFraction"))) { B = b->getVal(); Be = b->getError(); }
      fprintf(o, "%s,%g,%g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%g,%g,%.17g,%g,%.17g,%g,%.17g,%.17g,%.17g\n", mid ? "mid" : "fwd", ref.low, ref.high, B, Be, ref.value, ref.error,
              (B - ref.value) / ref.error, 100 * (B / ref.value - 1), r ? r->status() : -9, r ? r->covQual() : -9, r ? r->edm() : NAN,
              par<int>(f, "minimizeStatus"), par<int>(f, "hesseStatus"), par<double>(f, "massChi2"), par<int>(f, "massNdf"), par<double>(f, "timeChi2"), par<int>(f, "timeNdf"), par<double>(f, "rawFinalNll"), par<double>(f, "massPvalue"), par<double>(f, "timePvalue"));
    }
  fclose(o);
}
