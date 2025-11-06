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

void toyMC1()
{
  // --- Observable ---
  RooRealVar t("t", "decay time", -1.0, 5.0);

  // --- Resolution: 3-Gaussian mixture ---
  RooRealVar mean("mean", "mean", 0.0);
  RooRealVar sigma1("sigma1", "sigma1", 0.2, 0.01, 0.2);
  RooRealVar sigma2("sigma2", "sigma2", 0.20, 0.05, 0.5);
  RooRealVar sigma3("sigma3", "sigma3", 0.50, 0.1, 1.0);

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

  // --- Nonprompt lifetime ---
  RooRealVar tau("tau", "nonprompt lifetime", 2.8, 1.0, 6.0);
  RooDecay decay("decay", "nonprompt", t, tau, resModel, RooDecay::SingleSided);

  // --- Background lifetime ---
  RooRealVar tau_bkg("tau_bkg", "bkg lifetime", 1.0, 0.5, 3.0);
  RooDecay bkg("bkg", "background", t, tau_bkg, resModel, RooDecay::SingleSided);

  // --- Yields ---
  RooRealVar Npr("Npr", "yield prompt", 15000, 0, 1e7);
  RooRealVar Nnp("Nnp", "yield nonprompt", 15000, 0, 1e7);
  RooRealVar Nbkg("Nbkg", "yield bkg", 3000, 0, 1e7);

  RooExtendPdf extPR("extPR", "extended prompt", *prompt, Npr);
  RooExtendPdf extNP("extNP", "extended nonprompt", decay, Nnp);
  RooExtendPdf extBkg("extBkg", "extended background", bkg, Nbkg);

  RooAddPdf model("model", "PR+NP+Bkg",
                  RooArgList(extPR, extNP, extBkg));

  // --- Generate toy dataset ---
  auto data = model.generate(t, 50000);

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
  c->SaveAs("toyMC1.png");
}
