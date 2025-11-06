#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooAddModel.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooHist.h"

using namespace RooFit;

void toyMC2()
{
  // --- Observable ---
  RooRealVar t("t", "decay time", -1.0, 5.0);

  // --- Resolution: 3-Gaussian mixture ---
  RooRealVar mean("mean", "mean", 0.0);
  RooRealVar sigma1("sigma1", "sigma1", 0.02, 0.001, 0.1);
  RooRealVar sigma2("sigma2", "sigma2", 0.05, 0.001, 0.2);
  RooRealVar sigma3("sigma3", "sigma3", 0.15, 0.01, 0.5);

  RooGaussModel g1("g1", "res1", t, mean, sigma1);
  RooGaussModel g2("g2", "res2", t, mean, sigma2);
  RooGaussModel g3("g3", "res3", t, mean, sigma3);

  RooRealVar f1("f1", "frac1", 0.5, 0., 1.);
  RooRealVar f2("f2", "frac2", 0.3, 0., 1.);

  RooAddModel resModel("resModel", "3-Gaussian resolution",
                       RooArgList(g1, g2, g3),
                       RooArgList(f1, f2));

  // --- Prompt = delta ⊗ resolution ---
  RooAbsPdf *prompt = &resModel;

  // --- Nonprompt = exponential ⊗ resolution ---
  // RooRealVar tau("tau", "lifetime", 1.5, 5.0, 0.1); // 음수 slope
  // RooDecay decay("decay", "nonprompt", t, tau, resModel, RooDecay::SingleSided);
  RooRealVar tau1("tau1", "short lifetime", 0.3, 0.05, 2.0);
  RooRealVar tau2("tau2", "long lifetime", 1.5, 0.5, 20.0);
  RooRealVar fNP("fNP", "frac short", 0.5, 0., 1.);
  RooDecay decay1("decay1", "NP short", t, tau1, resModel, RooDecay::SingleSided);
  RooDecay decay2("decay2", "NP long", t, tau2, resModel, RooDecay::SingleSided);
  RooAddPdf decay("decay", "NP sum",
                      RooArgList(decay1, decay2), RooArgList(fNP));

  // --- Background = exponential ⊗ resolution ---
  // RooRealVar lambda("lambda", "bkg slope", 1.5, 0, 3); // 음수 slope
  // RooDecay bkg("bkg", "background", t, lambda, resModel, RooDecay::SingleSided);

  RooRealVar lambda1("lambda1", "short lifetime", 1.5, 0, 3); // 음수 slope
  RooRealVar lambda2("lambda2", "long lifetime", 0.5, 0, 1);  // 음수 slope
  RooRealVar fBkg("fBkg", "frac short", 0.6, 0.3, 1);           // 음수 slope
  RooDecay bkg1("bkg1", "NP short", t, lambda1, resModel, RooDecay::SingleSided);
  RooDecay bkg2("bkg2", "NP long", t, lambda2, resModel, RooDecay::SingleSided);
  RooAddPdf bkg("bkg", "NP sum",
                  RooArgList(bkg1, bkg2), RooArgList(fBkg));

  // --- Yields ---
  RooRealVar Npr("Npr", "yield prompt", 3000, 0, 1e6);
  RooRealVar Nnp("Nnp", "yield nonprompt", 2000, 0, 1e6);
  RooFormulaVar Nbkg("Nbkg", "@0 * 0.18", RooArgList(Nnp));
  // RooRealVar Nbkg("Nbkg", "yield bkg", 100, 0, 800);

  RooExtendPdf extPR("extPR", "extended prompt", *prompt, Npr);
  RooExtendPdf extNP("extNP", "extended nonprompt", decay, Nnp);
  RooExtendPdf extBkg("extBkg", "extended background", bkg, Nbkg);

  RooAddPdf model("model", "PR+NP+Bkg",
                  RooArgList(extPR, extNP, extBkg));

  // --- Generate toy dataset ---
  auto data = model.generate(t, 10000);

  // --- Fit ---
  model.fitTo(*data, Extended(kTRUE), PrintLevel(-1));

  // --- Plot main frame ---
  RooPlot *frame = t.frame();
  data->plotOn(frame, Name("data"));
  model.plotOn(frame, Name("curveAll")); // 전체 모델에 이름 지정
  model.plotOn(frame, Components("extPR"), LineStyle(kDashed), LineColor(kRed));
  model.plotOn(frame, Components("extNP"), LineStyle(kDashed), LineColor(kBlue));
  model.plotOn(frame, Components("extBkg"), LineStyle(kDashed), LineColor(kGreen));

  // --- Pull distribution (curveAll과 비교) ---
  RooHist *hpull = frame->pullHist("data", "curveAll");
  RooPlot *frame_pull = t.frame(Title("Pull Distribution"));
  frame_pull->addPlotable(hpull, "P");
  frame_pull->SetMinimum(-3);
  frame_pull->SetMaximum(3);

  // --- Canvas with two pads ---
  TCanvas *c = new TCanvas("c", "c", 800, 800);
  c->Divide(1, 2);

  // upper pad: logY fit result
  TPad *pad1 = (TPad *)c->cd(1);
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);
  pad1->SetLogy();
  frame->Draw();

  // --- Chi2/ndf 표시 ---
  double chi2 = frame->chiSquare("curveAll", "data");
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  // lower pad: pull
  TPad *pad2 = (TPad *)c->cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  frame_pull->Draw();

  // --- Save result ---
  c->SaveAs("toyMC2.png");
}
