#include <TStopwatch.h>
#include <RooRealVar.h>
#include "RooFit.h"
#include "RooMinimizer.h"

using namespace RooFit;

void toy_pr_3expo_bkg()
{
  cout << "=== start toy_pr_3expo_bkg()";
  TStopwatch t;
  t.Start();

  gSystem->mkdir("figs", true);

  // === mass ===
  RooRealVar *mass = new RooRealVar("mass", "mass", 5.15, 5.6, "GeV/c^{2}");

  // --- mass signal - two gaussians ---
  RooRealVar *Mean0_mass = new RooRealVar("M_{0}", "mass mean gauss 0", 5.28, 5.27, 5.29);
  RooRealVar *Sigma0_mass = new RooRealVar("#sigma_{M0}", "mass sigma gauss0", 0.001, 0.00001, 0.5);
  RooResolutionModel *g0_mass = new RooGaussModel("g0_mass", "mass resolution", *mass, *Mean0_mass, *Sigma0_mass);

  // RooRealVar *Mean1_mass = new RooRealVar("M_{1}", "mass mean gauss 1", 5.28, 5.27, 5.29);
  RooRealVar *Sigma1_mass = new RooRealVar("#sigma_{M1}", "mass sigma gauss0", 0.005, 0.00001, 0.5);
  RooResolutionModel *g1_mass = new RooGaussModel("g1_mass", "mass resolution", *mass, *Mean0_mass, *Sigma1_mass);

  // fractions of gauss0, 1
  RooRealVar *frac_g0_mass = new RooRealVar("f_g0_M", "gauss resol.fracion 0", 0.05, 0.01, 1);
  RooRealVar *frac_g1_mass = new RooRealVar("f_g1_M", "gauss resol.fracion 1", 0.3, 0.01, 1);

  // combine gauss
  RooAbsPdf *mass_b_pdf = new RooAddPdf("mass_b_pdf", "dual gaussian Pdf", RooArgList(*g0_mass, *g1_mass), RooArgList(*frac_g0_mass));

  // extend
  RooRealVar *N_pr = new RooRealVar("N_{pr}", "signal number B Jpsi K", 6000, 1, 10000);
  RooAbsPdf *mass_b_epdf = new RooExtendPdf("mass_b_epdf", "mass sig B 2 Jpsi K", *mass_b_pdf, *N_pr);

  // --- ctau signal ---
  RooRealVar *tom = new RooRealVar("tom", "proper decay length", -0.2, 0.4, "cm");
  RooRealVar *b_tau = new RooRealVar("c#tau_{b}", "B ctau", 0.05, 0.04, 0.055, "cm");
  RooRealVar *Mean_tom = new RooRealVar("tomMean", "proper time mean", 0);
  RooRealVar *ratio_tom = new RooRealVar("s_{to}", "ratio", 1, 0.01, 50); // ratio of what?
  RooRealVar *tom_err = new RooRealVar("tom_err", "tom_err", 0, 1.0, "cm");

  RooRealVar *Sigma0_tom = new RooRealVar("#sigma_{t0}", "proper time sigma", 0.00265, 0, 0.005);
  // per event err??
  // RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *tom_err));
  RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *Sigma0_tom));

  RooRealVar *frac_g_tom = new RooRealVar("f_g_#sigma_{t0}", "gauss resol. fracion", 0.8048, 0.01, 1);

  RooResolutionModel *g0_tom = new RooGaussModel("g0_tom", "proper time resolution", *tom, *Mean_tom, *Sigma0_tom);
  RooResolutionModel *g1_tom = new RooGaussModel("g1_tom", "proper time resolution", *tom, *Mean_tom, *Sigma1_tom);

  // one gauss signal
  RooResolutionModel *properTimeRes = new RooGaussModel("properTimeRes", "proper time resolution", *tom, *Mean_tom, *Sigma1_tom);

  // two gauss signal
  // RooResolutionModel *properTimeRes = new RooAddModel("properTimeRes", "dual gaussain PDF", RooArgList(*g0_tom, *g1_tom), RooArgList(*frac_g_tom));

  RooAbsPdf *time_b_pdf = new RooDecay("time_b_pdf", "time pdf B2 Jpsi K", *tom, *b_tau, *properTimeRes, RooDecay::SingleSided);
  // OK, signal time PDF is ready

  // --- combine signal PDFs (mass + time)
  RooAbsPdf *masstimeSigpdf = new RooProdPdf("masstimeSigpdf", "signal m* t", RooArgSet(*time_b_pdf, *mass_b_epdf));

  // --- NP ---
  RooRealVar *N_np = new RooRealVar("N_{NP}", "signal number B Jpsi K", 4000, 1, 100000);
  RooAbsPdf *mass_np_epdf = new RooExtendPdf("mass_np_epdf", "mass sig B 2 Jpsi K", *mass_b_pdf, *N_np);
  auto masstimeSigpdf2 = new RooProdPdf("masstimeSigpdf2", "signal m*t", RooArgSet(*g0_tom, *mass_np_epdf));

  // === Jpsi bg ===
  // BB mass -> mass와 slope는 공유
  RooRealVar *slope = new RooRealVar("slope", "slope", 0.01, -1, 1);
  // RooAbsPdf *mass_BB_pdf = new RooPolynomial("mass_QCD_pdf", "mass pdf QCD", *mass, *slope);
  RooAbsPdf *mass_BB_pdf = new RooChebychev("mass_QCD_pdf", "mass pdf QCD",
                                            *mass, *slope);

  // prompt BB time - NP인데 Res과 함쳤다는 뜻?
  // --- base tau
  RooRealVar tau0("tau0", "base B ctau", 0.01, 1e-4, 0.2, "cm");

  // ratio parameter
  RooRealVar r_left("r_left", "scale factor left", 1.0, 0.1, 3.0);
  RooRealVar r_mid("r_mid", "scale factor middle", 1.0, 0.1, 3.0);
  RooRealVar r_right("r_right", "scale factor right", 3.0, 0.1, 10.0);

  // tau0.setMin(1e-4);
  // r_left.setMin(0.1);
  // r_mid.setMin(0.1);
  // r_right.setMin(0.1);

  // tau for bkg pdfs
  RooFormulaVar tau_left("tau_left", "@0*@1", RooArgList(tau0, r_left));
  RooFormulaVar tau_mid("tau_mid", "@0*@1", RooArgList(tau0, r_mid));
  RooFormulaVar tau_right("tau_right", "@0*@1", RooArgList(tau0, r_right));

  // --- 3 decay PDF
  RooDecay time_BB_left("time_BB_left", "Left decay",
                        *tom, tau_left, *properTimeRes, RooDecay::Flipped);

  RooDecay time_BB_mid("time_BB_mid", "Center decay",
                       *tom, tau_mid, *properTimeRes, RooDecay::DoubleSided);

  RooDecay time_BB_right("time_BB_right", "Right decay",
                         *tom, tau_right, *properTimeRes, RooDecay::SingleSided);

  RooRealVar f_left("f_left", "fraction left", 0.3, 0.0, 1.0);
  RooRealVar f_mid("f_mid", "fraction middle", 0.3, 0.0, 1.0);

  auto time_BB_pdf = new RooAddPdf("time_BB_pdf", "sum of 3 decays",
                                   RooArgList(time_BB_left, time_BB_mid, time_BB_right),
                                   RooArgList(f_left, f_mid));

  RooRealVar *N_BB = new RooRealVar("N_{BB}", "Signal number BB", 10000, 1, 10000);
  RooAbsPdf *mass_BB_epdf = new RooExtendPdf("mass_BB_epdf", "mass BB", *mass_BB_pdf, *N_BB);

  // bkg BB mass + timie
  RooAbsPdf *masstimeBkgpdf = new RooProdPdf("masstimeBkgpdf", "background2 m * t", RooArgSet(*time_BB_pdf, *mass_BB_epdf));
  /// OK bkg2 BB is readg

  // === Sig + Bkg PDF ===
  RooAddPdf *pdf_fit = new RooAddPdf(
      "pdf_fit", "Total(B+S) PDF",
      RooArgSet(*masstimeSigpdf, *masstimeSigpdf2, *masstimeBkgpdf));

  // === Fit ===
  /// Real data generate -> Toy MC로 테스트
  RooArgSet *pdf_obs = new RooArgSet(*mass, *tom);
  // RooDataSet *dataB = RooDataSet::read("fitsample.txt", RooArgList(*pdf_obs));
  std::unique_ptr<RooDataSet> data(pdf_fit->generate(RooArgSet(*mass, *tom), 20000));

  // fit
  // Define difference -log(L) for B->ff, B->ff fbar and B->fbar f
  auto nllB = pdf_fit->createNLL(*data, Extended(true), NumCPU(32), EvalBackend("legacy"));

  RooMinimizer m1(*nllB);
  m1.setVerbose(kFALSE);
  // m1.setLogFile(logFile); - 이게 나을지도??
  m1.setProfile(0);
  m1.setPrintLevel(-1);
  m1.setStrategy(2);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(2).removeTopic(RooFit::Eval);
  gErrorIgnoreLevel = kFatal;
  m1.migrad();
  m1.hesse();
  m1.migrad();
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  std::unique_ptr<RooFitResult> result(m1.save());

  {
  TCanvas c1("c1", "", 800, 800);
  c1.Divide(1, 2);
  TPad *pad1 = (TPad *)c1.cd(1);
  gPad->SetLogy();
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);

  TPad *pad2 = (TPad *)c1.cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);

  // main plot
  pad1->cd();
  RooPlot *frame = tom->frame();
  data->plotOn(frame, Name("data"));
  pdf_fit->plotOn(frame, Name("model"));
  pdf_fit->plotOn(frame,
                  Components("masstimeSigpdf"),
                  Name("np"),
                  LineColor(kRed),
                  LineStyle(kDashed),
                  LineWidth(4));
  pdf_fit->plotOn(frame,
                  Components("masstimeSigpdf2"),
                  Name("pr"),
                  LineColor(kOrange),
                  LineStyle(kDashed),
                  LineWidth(4));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf"),
                  Name("bkg"),
                  LineColor(kGreen+2),
                  LineStyle(kDashed),
                  LineWidth(4));
  
  // pdf_fit->plotOn(frame,
  //                 Components("time_BB_left"),
  //                 Name("time_BB_left"),
  //                 LineColor(kMagenta + 2),
  //                 LineStyle(kDashDotted),
  //                 LineWidth(2));
  // pdf_fit->plotOn(frame,
  //                 Components("time_BB_mid"),
  //                 Name("time_BB_mid"),
  //                 LineColor(kViolet),
  //                 LineStyle(kDashDotted),
  //                 LineWidth(2));
  // pdf_fit->plotOn(frame,
  //                 Components("time_BB_right"),
  //                 Name("time_BB_right"),
  //                 LineColor(kViolet + 4),
  //                 LineStyle(kDashDotted),
  //                 LineWidth(2));
  frame->Draw();

  // --- legend ---
  auto leg = new TLegend(0.6, 0.65, 0.95, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(frame->findObject("data"), "data", "pe");
  leg->AddEntry(frame->findObject("model"), "Total fit", "e");
  leg->AddEntry(frame->findObject("pr"), "Prompt J/#psi", "e");
  leg->AddEntry(frame->findObject("np"), "NP J/#psi", "e");
  leg->AddEntry(frame->findObject("bkg"), "Bkg", "e");
  leg->Draw();

  // pull
  pad2->cd();
  RooHist *hpull = frame->pullHist("data", "model"); // data - fit / error
  RooPlot *frame_pull = tom->frame();
  frame_pull->addPlotable(hpull, "P");
  frame_pull->SetTitle("");
  frame_pull->GetYaxis()->SetTitle("Pull");
  frame_pull->GetYaxis()->SetNdivisions(505);
  frame_pull->GetYaxis()->SetTitleSize(0.08);
  frame_pull->GetYaxis()->SetLabelSize(0.08);
  frame_pull->GetXaxis()->SetTitleSize(0.1);
  frame_pull->GetXaxis()->SetLabelSize(0.08);
  frame_pull->Draw();

  // --- chi2/ndf ---
  int nFitParam = pdf_fit->getParameters(*data)->selectByAttrib("Constant", kFALSE)->getSize();
  double chi2 = frame->chiSquare("model", "data", nFitParam); // (chi2/ndf)
  std::cout << "Chi2/ndf = " << chi2 << std::endl;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.DrawLatex(0.85, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  c1.SaveAs("figs/toy_pr_3expo_ctau.png");
  // c1.SaveAs("figs/toy_pr_3expo_ctau.pdf");
  }

  // === mass plot ===
  {
  TCanvas c2("c2", "", 800, 800);
  c2.Divide(1, 2);
  TPad *pad1 = (TPad *)c2.cd(1);
  // gPad->SetLogy();
  pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.02);

  TPad *pad2 = (TPad *)c2.cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.25);

  // main plot
  pad1->cd();
  RooPlot *frame = mass->frame();
  data->plotOn(frame, Name("data"));
  pdf_fit->plotOn(frame, Name("model"));
  pdf_fit->plotOn(frame,
                  Components("masstimeSigpdf"),
                  Name("np"),
                  LineColor(kRed),
                  LineStyle(kDotted),
                  LineWidth(3));
  pdf_fit->plotOn(frame,
                  Components("masstimeSigpdf2"),
                  Name("pr"),
                  LineColor(kOrange),
                  LineStyle(kDashed),
                  LineWidth(2));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf"),
                  Name("bkg"),
                  LineColor(kGreen + 2),
                  LineStyle(kDashed),
                  LineWidth(2));
  frame->Draw();

  // --- legend ---
  auto leg = new TLegend(0.6, 0.65, 0.95, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(frame->findObject("data"), "data", "pe");
  leg->AddEntry(frame->findObject("model"), "Total fit", "e");
  leg->AddEntry(frame->findObject("pr"), "Prompt J/#psi", "e");
  leg->AddEntry(frame->findObject("np"), "NP J/#psi", "e");
  leg->AddEntry(frame->findObject("bkg"), "Bkg", "e");
  leg->Draw();

  // pull
  pad2->cd();
  RooHist *hpull = frame->pullHist("data", "model"); // data - fit / error
  RooPlot *frame_pull = mass->frame();
  frame_pull->addPlotable(hpull, "P");
  frame_pull->SetTitle("");
  frame_pull->GetYaxis()->SetTitle("Pull");
  frame_pull->GetYaxis()->SetNdivisions(505);
  frame_pull->GetYaxis()->SetTitleSize(0.08);
  frame_pull->GetYaxis()->SetLabelSize(0.08);
  frame_pull->GetXaxis()->SetTitleSize(0.1);
  frame_pull->GetXaxis()->SetLabelSize(0.08);
  frame_pull->Draw();

  // --- chi2/ndf ---
  int nFitParam = pdf_fit->getParameters(*data)->selectByAttrib("Constant", kFALSE)->getSize();
  double chi2 = frame->chiSquare("model", "data", nFitParam); // (chi2/ndf)
  std::cout << "Chi2/ndf = " << chi2 << std::endl;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.08);
  latex.DrawLatex(0.85, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  c2.SaveAs("figs/toy_pr_3expo_mass.png");
  // c2.SaveAs("figs/toy_pr_3expo_mass.pdf");
  }

  result->Print("v");

  cout << "b fraction = N_{NP} / (N_NP + N_PR) = " << N_np->getVal() / (N_np->getVal() + N_pr->getVal()) << "\n";

  cout
      << "=== finish toy_pr_3expo_bkg()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}