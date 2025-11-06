#include <TStopwatch.h>
#include <RooRealVar.h>
#include "RooFit.h"
#include "RooMinimizer.h"

using namespace RooFit;

void toy_2009_model_test()
{
  cout << "=== start toy_2009_model_test()";
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
  RooRealVar *N_b = new RooRealVar("N_{b}", "signal number B Jpsi K", 1500, 1, 3000);
  RooAbsPdf *mass_b_epdf = new RooExtendPdf("mass_b_epdf", "mass sig B 2 Jpsi K", *mass_b_pdf, *N_b);

  // --- ctau signal ---
  RooRealVar *tom = new RooRealVar("tom", "proper decay length", -0.05, 0.4, "cm");
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

  // === Jpsi bkg ===
  // ctau에서 Prompt에 해당하는 것과 묶었으니까 prompt J/psi mass bkg라고 말할 수 있다. -> vs BB bkg
  // Prompt -> 진짜 관심 있는 Jpsi
  // 총 Sig + (Jpsi bkg + BB bkg)로 구성
  // --- prompt J/psi mass ---
  RooRealVar *slope = new RooRealVar("slope", "slope", 0.0, -0.1, 0.01);
  RooAbsPdf *mass_pJpsi_pdf = new RooPolynomial("mass_pJpsi_pdf", "mass pdf pJpsi", *mass, *slope); // 이게 prompt인지 어떻게 알고?
  RooRealVar *N_pJpsi = new RooRealVar("N_{p-J/#psi}", "signal number prompt Jpsi", 10000, 1, 20000);
  RooAbsPdf *mass_pJpsi_epdf = new RooExtendPdf("mass_pJpsi_epdf", "mass sig pJpsi", *mass_pJpsi_pdf, *N_pJpsi);

  // --- prompt J/psi ctau ---
  // 이름이 몇 번 중복되는데 ctua_gm1과 같다.
  RooAbsPdf *time_pJpsi_pdf = new RooGaussModel("time_pJpsi_pdf", "proper time resolution", *tom, *Mean_tom, *Sigma1_tom);

  // --- bkg Jpsi mass + time ---
  RooAbsPdf *masstimeBkgpdf = new RooProdPdf("masstimeBkgpdf", "background m * t", RooArgSet(*time_pJpsi_pdf, *mass_pJpsi_epdf));
  //// OK background jpsi is ready

  // === Jpsi bg ===
  // BB mass -> mass와 slope는 공유
  RooAbsPdf *mass_BB_pdf = new RooPolynomial("mass_QCD_pdf", "mass pdf QCD", *mass, *slope);
  RooRealVar *N_BB = new RooRealVar("N_{BB}", "Signal number BB", 500, 1, 2000);
  RooAbsPdf *mass_BB_epdf = new RooExtendPdf("mass_BB_epdf", "mass BB", *mass_BB_pdf, *N_BB);

  // prompt BB time - NP인데 Res과 함쳤다는 뜻?
  RooAbsPdf *time_BB_pdf = new RooDecay("time_BB_pdf", "time pdf BB", *tom, *b_tau, *properTimeRes, RooDecay::SingleSided);

  // bkg BB mass + timie
  RooAbsPdf *masstimeBkgpdf2 = new RooProdPdf("masstimeBkgpdf2", "background2 m * t", RooArgSet(*time_BB_pdf, *mass_BB_epdf));
  /// OK bkg2 BB is readg


  // === Sig + Bkg PDF ===
  RooAddPdf *pdf_fit = new RooAddPdf(
      "pdf_fit", "Total(B+S) PDF",
      RooArgSet(*masstimeSigpdf, *masstimeBkgpdf, *masstimeBkgpdf2));

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
  m1.migrad();
  m1.hesse();
  m1.migrad();

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
                  Name("signal"),
                  LineColor(kRed),
                  LineStyle(kDotted),
                  LineWidth(3));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf"),
                  Name("pr_bkg"),
                  LineColor(kOrange),
                  LineStyle(kDashed),
                  LineWidth(2));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf2"),
                  Name("np_bkg"),
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
  leg->AddEntry(frame->findObject("signal"), "Signal", "e");
  leg->AddEntry(frame->findObject("pr_bkg"), "Prompt bkg", "e");
  leg->AddEntry(frame->findObject("np_bkg"), "Nonprompt bkg", "e");
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

  c1.SaveAs("figs/toy_2009_model_ctau.png");
  // c1.SaveAs("figs/toy_2009_model_ctau.pdf");
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
                  Name("signal"),
                  LineColor(kRed),
                  LineStyle(kDotted),
                  LineWidth(3));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf"),
                  Name("pr_bkg"),
                  LineColor(kOrange),
                  LineStyle(kDashed),
                  LineWidth(2));
  pdf_fit->plotOn(frame,
                  Components("masstimeBkgpdf2"),
                  Name("np_bkg"),
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
  leg->AddEntry(frame->findObject("signal"), "Signal", "e");
  leg->AddEntry(frame->findObject("pr_bkg"), "Prompt bkg", "e");
  leg->AddEntry(frame->findObject("np_bkg"), "Nonprompt bkg", "e");
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

  c2.SaveAs("figs/toy_2009_model_mass.png");
  // c2.SaveAs("figs/toy_2009_model_mass.pdf");
  }

  result->Print("v");

  cout << "=== finish toy_2009_model_test()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());

  // fitResult->Print("V");
}