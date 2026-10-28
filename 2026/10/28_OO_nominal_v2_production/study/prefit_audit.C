// Prefit quality per bin: status/covQual/chi2 of every prefit stage, PR resolution widths and fractions,
// and floating parameters stuck at a limit or with error > |value| (unused components).
//   root -l -b -q 'study/prefit_audit.C("roots","study/prefit_audit.txt")'
#include "TFile.h"
#include "TParameter.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include <cmath>
#include <cstdio>
#include <string>

static double val(TFile &f, const char *n) {
  if (auto *p = f.Get<TParameter<double>>(n)) return p->GetVal();
  if (auto *p = f.Get<TParameter<int>>(n)) return p->GetVal();
  return NAN;
}

static void stage(FILE *o, const char *base, const char *y, const char *st, const char *file, const char *tag, const char *resName) {
  TFile f(Form("%s/%s/%s/%s_%s.root", base, y, st, file, tag));
  if (f.IsZombie()) { fprintf(o, "  %-9s missing\n", st); return; }
  auto *r = f.Get<RooFitResult>(resName);
  fprintf(o, "  %-9s st%d cq%d chi2/ndf %.3f", st, r ? r->status() : -9, r ? r->covQual() : -9, val(f, "prefitChi2NdfValue"));
  const double rob = val(f, "prefitChi2NdfRobust");
  if (std::isfinite(rob)) fprintf(o, " (robust %.3f)", rob);
  if (!strcmp(st, "ctau_pr") || !strcmp(st, "ctau_bkg") || !strcmp(st, "ctau_np")) {
    const int n = int(val(f, "nResolutionComponents"));
    fprintf(o, "  res n=%d sigma", n);
    for (int i = 1; i <= n; ++i) fprintf(o, " %.4f", val(f, Form("ctauTime%dScale", i)));
    fprintf(o, " frac");
    for (int i = 1; i < n; ++i) fprintf(o, " %.4f", val(f, Form("ctauFrac%d", i)));
  }
  fprintf(o, "\n");
  if (!r) return;
  for (auto *a : r->floatParsFinal()) {
    auto *v = dynamic_cast<RooRealVar *>(a);
    if (!v) continue;
    const double lo = v->getMin(), hi = v->getMax(), x = v->getVal(), e = v->getError(), w = hi - lo;
    const bool atLimit = std::isfinite(w) && w > 0 && (x - lo < 1e-3 * w || hi - x < 1e-3 * w);
    const bool loose = e > std::fabs(x) && std::fabs(x) > 0;
    if (atLimit || loose)
      fprintf(o, "      ! %-32s %.5g +- %.3g  [%g, %g]%s%s\n", v->GetName(), x, e, lo, hi, atLimit ? " AT-LIMIT" : "", loose ? " err>|val|" : "");
  }
}

void prefit_audit(const char *base = "roots", const char *out = "study/prefit_audit.txt") {
  FILE *o = fopen(out, "w");
  const double fwd[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20}, mid[] = {7, 8, 9, 10, 12, 14, 16, 20, 30};
  for (int m = 0; m < 2; ++m) {
    const char *y = m ? "y0.00_1.60" : "y1.60_2.40";
    const double *e = m ? mid : fwd; const int n = m ? 8 : 12;
    for (int i = 0; i < n; ++i) {
      const std::string tag = Form("%s_pt%.2f_%.2f", y, e[i], e[i + 1]);
      fprintf(o, "%s %g-%g\n", m ? "mid" : "fwd", e[i], e[i + 1]);
      stage(o, base, y, "mc_mass", "mc_mass_model", tag.c_str(), "fit_result");
      stage(o, base, y, "ctau_pr", "ctau_resolution", tag.c_str(), "timeResult");
      stage(o, base, y, "ctau_np", "ctau_np_model", tag.c_str(), "timeResult");
      stage(o, base, y, "mass", "mass_model", tag.c_str(), "fit_result");
      stage(o, base, y, "ctau_bkg", "ctau_bkg_fitresult", tag.c_str(), "timeResult");
    }
  }
  fclose(o);
}
