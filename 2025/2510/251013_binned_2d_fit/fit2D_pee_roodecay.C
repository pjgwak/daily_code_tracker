#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include <memory>
using namespace RooFit;

void fit2D_pee_roodecay()
{
  // --- Observables ---
  RooRealVar mass("mass", "mass", 2.6, 3.5);
  RooRealVar ct("ct", "ct", -0.5, 2.0);
  RooRealVar dt("dt", "per-event error", 0.005, 0.2);

  // --- Mass PDF ---
  RooRealVar m0("m0", "mean", 3.096, 2.9, 3.2);
  RooRealVar ms("ms", "sigma", 0.020, 0.005, 0.050);
  RooGaussian massPdf("massPdf", "mass pdf", mass, m0, ms);

  // --- Resolution: sigma = sf * dt  (PEE) ---
  RooRealVar bias("bias", "bias", 0.0);
  bias.setConstant(kTRUE);
  RooRealVar sf("sf", "scale", 1.40, 0.6, 3.0);
  dt.setMin(1e-6);
  sf.setMin(1e-3);
  RooGaussModel res("res", "res", ct, bias, sf, dt);

  // --- 3 RooDecay ---
  RooRealVar tau1("tau1", "tau1", 0.25, 0.05, 1.5);
  RooRealVar tau2("tau2", "tau2", 0.60, 0.05, 3.0);
  RooRealVar tau3("tau3", "tau3", 0.12, 0.02, 1.0);
  RooDecay d1("d1", "NP1", ct, tau1, res, RooDecay::SingleSided);
  RooDecay d2("d2", "NP2", ct, tau2, res, RooDecay::SingleSided);
  RooDecay d3("d3", "Sym", ct, tau3, res, RooDecay::DoubleSided);

  // --- p(dt) (Punzi term) : RooHistPdf -> RooHistFunc wrapped by RooGenericPdf ---
  RooDataSet toy_dt("toy_dt", "toy_dt", RooArgSet(dt));
  for (int i = 0; i < 50000; ++i)
  {
    double v = gRandom->Gaus(0.04, 0.01);
    if (v < 0.005)
      v = 0.005;
    if (v > 0.2)
      v = 0.2;
    dt.setVal(v);
    toy_dt.add(RooArgSet(dt));
  }
  RooDataHist hdt("hdt", "hdt", RooArgList(dt), toy_dt);
  // RooHistPdf histPdf("histPdf", "p(dt)", RooArgList(dt), hdt, 0);

  RooHistFunc fdt("fdt", "f(dt) raw", RooArgList(dt), hdt);
  std::unique_ptr<RooAbsReal> I(fdt.createIntegral(RooArgSet(dt)));
  RooFormulaVar fdt_unit("fdt_unit", "@0/@1", RooArgList(fdt, *I)); // make norm become 1
  RooGenericPdf pdt("pdt", "p(dt) unit", "@0", RooArgList(fdt_unit));

  RooHistPdf pdt_hist("pdt_hist", "p(dt) hist", RooArgList(dt), hdt, 0); // to compare

  // --- plotting 1 - compare RooHistPdf and RooHistFunc ---
  const double dtMin = dt.getMin();
  const double dtMax = dt.getMax();
  dt.setRange("dtPlot", dtMin, dtMax);

  TCanvas c_dt("c_dt", "p(dt) comparison", 800, 600);
  RooPlot *f_dt = dt.frame(Title("p(dt) : HistPdf vs HistFunc(GenericPdf)"), Range(dtMin, dtMax));

  hdt.plotOn(f_dt, LineColor(kBlack), MarkerColor(kBlack), DataError(RooAbsData::SumW2), Name("h_data_dt"));

  pdt_hist.plotOn(f_dt,
                  // Range("dtPlot"), NormRange("dtPlot"),
                  LineColor(kRed), LineWidth(2), Name("p_hist"));
  
  // Don't use NormRange for the RooGenericPdf -> It tries numeric integral in that range -> Too slow
  pdt.plotOn(f_dt,
            //  Range("dtPlot"), NormRange("dtPlot"),
             LineColor(kBlue + 1), LineStyle(kDashed), LineWidth(2), Name("p_func"));

  TLegend leg_dt(0.58, 0.70, 0.88, 0.88);
  leg_dt.SetBorderSize(0);
  leg_dt.SetFillStyle(0);
  leg_dt.AddEntry(f_dt->findObject("h_data_dt"), "Data (hdt)", "lep");
  leg_dt.AddEntry(f_dt->findObject("p_hist"), "RooHistPdf p(dt)", "l");
  leg_dt.AddEntry(f_dt->findObject("p_func"), "RooHistFunc -> GenericPdf", "l");
  leg_dt.Draw();

  f_dt->Draw();
  c_dt.SaveAs("pdt_compare.png");
  // --- end of plotting 1 ---

  TH1 *hdt_h1 = hdt.createHistogram("hdt_h1", dt);
  dt.setRange(hdt_h1->GetXaxis()->GetXmin() + 1e-6, hdt_h1->GetXaxis()->GetXmax());

  // --- combine components: comp_i = massPdf × (d_i × pdt) ---
  // Conditional(RooArgSet(d_i), RooArgSet(dt)): "d_i" is conditional to dt
  RooProdPdf d1_punzi("d1_punzi", "d1×pdt", RooArgSet(d1, pdt),
                      Conditional(RooArgSet(d1), RooArgSet(dt)));
  RooProdPdf d2_punzi("d2_punzi", "d2×pdt", RooArgSet(d2, pdt),
                      Conditional(RooArgSet(d2), RooArgSet(dt)));
  RooProdPdf d3_punzi("d3_punzi", "d3×pdt", RooArgSet(d3, pdt),
                      Conditional(RooArgSet(d3), RooArgSet(dt)));

  RooProdPdf comp1("comp1", "mass×(d1×pdt)", RooArgSet(massPdf, d1_punzi));
  RooProdPdf comp2("comp2", "mass×(d2×pdt)", RooArgSet(massPdf, d2_punzi));
  RooProdPdf comp3("comp3", "mass×(d3×pdt)", RooArgSet(massPdf, d3_punzi));

  // --- Extended yield -> To avoid RooGaussian out of range issue ---
  RooRealVar N1("N1", "yield1", 15000, 0, 1e8);
  RooRealVar N2("N2", "yield2", 12000, 0, 1e8);
  RooRealVar N3("N3", "yield3", 10000, 0, 1e8);
  RooExtendPdf e1("e1", "e1", comp1, N1);
  RooExtendPdf e2("e2", "e2", comp2, N2);
  RooExtendPdf e3("e3", "e3", comp3, N3);

  // --- final pdf ---
  RooAddPdf model("model", "sum of components", RooArgList(e1, e2, e3));

  // --- make toyMC ---
  std::unique_ptr<RooDataSet> toy(
      model.generate(RooArgSet(mass, ct, dt), 50000) //
  );

  // --- Fit (unbinned + extended + PEE) ---
  std::unique_ptr<RooFitResult> fitres(
      model.fitTo(*toy,
                  Save(true),
                  SumW2Error(true),
                  ConditionalObservables(RooArgSet(dt)), // PEE
                  PrintEvalErrors(-1)
                  // NumCPU(32)
                  ));

  // --- plot ctau3D ---
  // Full conditional PDF -> Don't use ProjWData
  TCanvas c1("c1", "ct proj (extended punzi)", 800, 600);
  RooPlot *fct = ct.frame(Title("Projection on ct"));
  toy->plotOn(fct, Binning(120), DataError(RooAbsData::SumW2));
  model.plotOn(fct, LineColor(kBlue), NumCPU(16));
  model.plotOn(fct, Components("e1"), LineStyle(kDashed), LineColor(kRed), NumCPU(16));
  model.plotOn(fct, Components("e2"), LineStyle(kDashed), LineColor(kGreen + 2), NumCPU(16));
  model.plotOn(fct, Components("e3"), LineStyle(kDashed), LineColor(kMagenta + 2), NumCPU(16));
  fct->Draw();
  c1.SaveAs("punzi_ext_ct.png");

  // --- plot mass ---
  TCanvas c2("c2", "mass proj (extended punzi)", 800, 600);
  RooPlot *fm = mass.frame(Title("Projection on mass"));
  toy->plotOn(fm, Binning(80), DataError(RooAbsData::SumW2));
  model.plotOn(fm, LineColor(kBlue), NumCPU(16));
  model.plotOn(fm, Components("e1"), LineStyle(kDashed), LineColor(kRed), NumCPU(16));
  model.plotOn(fm, Components("e2"), LineStyle(kDashed), LineColor(kGreen + 2), NumCPU(16));
  model.plotOn(fm, Components("e3"), LineStyle(kDashed), LineColor(kMagenta + 2), NumCPU(16));
  fm->Draw();
  c2.SaveAs("punzi_ext_mass.png");
}
