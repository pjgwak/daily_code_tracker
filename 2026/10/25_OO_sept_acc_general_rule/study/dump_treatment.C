// Dump every final-fit parameter with its treatment (fixed / constrained / free) to a tab-separated table.
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include <cstdio>
#include <set>
#include <string>
void dump_treatment(const char *base, const char *out) {
  FILE *o = fopen(out, "w");
  fprintf(o, "bin\tparameter\ttreatment\tvalue\terror\n");
  const double fe[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20}, me[] = {7, 8, 9, 10, 12, 14, 16, 20, 30};
  for (int mid = 0; mid < 2; ++mid)
    for (int i = 0; i < (mid ? 8 : 12); ++i) {
      const double lo = mid ? me[i] : fe[i], hi = mid ? me[i + 1] : fe[i + 1];
      const char *tag = mid ? "y0.00_1.60" : "y1.60_2.40";
      const char *path = Form("%s/roots/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root", base, tag, tag, lo, hi);
      if (gSystem->AccessPathName(path)) continue;
      TFile f(path);
      auto *r = f.IsZombie() ? nullptr : f.Get<RooFitResult>("modelResult");
      if (!r) continue;
      std::set<std::string> con;
      if (auto *n = f.Get<TNamed>("constrainedParameters")) {
        std::string s = n->GetTitle(), cur; bool inList = false;
        for (char c : s) {
          if (c == ':') { inList = true; cur.clear(); }
          else if (c == ',' || c == ';') { if (inList && !cur.empty()) con.insert(cur); cur.clear(); if (c == ';') inList = false; }
          else cur += c;
        }
      }
      const char *bin = Form("%s_%g_%g", mid ? "mid" : "fwd", lo, hi);
      for (auto *a : r->floatParsFinal()) { auto *v = (RooRealVar *)a; fprintf(o, "%s\t%s\t%s\t%.8g\t%.8g\n", bin, v->GetName(), con.count(v->GetName()) ? "constrained" : "free", v->getVal(), v->getError()); }
      for (auto *a : r->constPars()) { auto *v = (RooRealVar *)a; fprintf(o, "%s\t%s\tfixed\t%.8g\t0\n", bin, v->GetName(), v->getVal()); }
    }
  fclose(o);
}
