#include <gsl/gsl_errno.h>
#include <TStopwatch.h>
#include <RooRealVar.h>
#include "RooFit.h"
#include "RooMinimizer.h"

using namespace RooFit;

void final_2d_fit()
{
  cout << "=== start final_2d_fit() ===\n";
  TStopwatch t;
  t.Start();

  double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
  gSystem->mkdir("figs", true);

  // === mass ===
  int nMassBins = 200;
  RooRealVar mass("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
  mass.setRange("massRange", massLow, massHigh);
  RooBinning binningMass(massLow, massHigh);
  binningMass.addUniform(nMassBins, massLow, massHigh);
  mass.setBinning(binningMass);
  
  // RooBinning binningMass(2.6, 3.5);
  // binningMass.addUniform(nMassBins, 2.6, 3.5);
  // mass.setBinning(binningMass);

  // keep ctau bin width about 83 micro-meter
  int nCtauBins = 200;
  double ctMin = -0.5, ctMax = 2.5;
  RooRealVar ctau3D("ctau3D", "ctau3D", ctMin, ctMax, "mm");

  RooBinning binningCtau(ctMin, ctMax);
  binningCtau.addUniform(nCtauBins, ctMin, ctMax);
  ctau3D.setBinning(binningCtau);

  // RooBinning binningCtau(ctMin, ctMax);
  // binningCtau.addUniform(nCtauBins, ctMin, ctMax);
  // ctau3D.setBinning(binningCtau);

  // --- 공통 평균 ---
  RooRealVar mean("mean", "mean of CB", 3.096, 3.09, 3.12);

  // === CB1 ===
  RooRealVar sigma1("sigma1", "CB1 sigma", 0.01, 0.005, 0.1);
  RooRealVar alpha1("alpha1", "CB1 alpha", 1.5, 0.1, 10.0);
  RooRealVar n1("n1", "CB1 n", 3.0, 1.0, 10.0);

  RooCBShape cb1("cb1", "Crystal Ball 1",
                 mass, mean, sigma1, alpha1, n1);

  // === CB2 ===
  RooRealVar sigma2("sigma2", "CB2 sigma", 0.01, 0.005, 0.1);
  RooRealVar alpha2("alpha2", "CB2 alpha", -1.5, -10.0, -0.01);
  RooRealVar n2("n2", "CB2 n", 3.0, 1.0, 10.0);

  RooCBShape cb2("cb2", "Crystal Ball 2",
                 mass, mean, sigma2, alpha2, n1);

  // === Gaussian ===
  RooRealVar sigmaG("sigmaG", "Gaussian sigma", 0.1, 0.01, 0.1);
  RooGaussian gaus("gaus", "extra Gaussian",
                   mass, mean, sigmaG);

  // === fraction 파라미터 ===
  // f1 = CB1 비율, f2 = CB2 비율, 나머지 = Gauss
  RooRealVar f1("f1", "frac CB1", 0.5, 0.0, 1.0);
  RooRealVar f2("f2", "frac CB2", 0.3, 0.0, 1.0);

  // === 최종 Signal PDF ===
  RooAddPdf signal("signal", "CB1 + CB2 + Gauss",
                   RooArgList(cb1, cb2, gaus),
                   RooArgList(f1, f2), kTRUE);

  // === mass bkg ===
  // --- Chebychev 계수 (2차이므로 c0, c1, c2 필요) ---
  RooRealVar s1("s1", "Chebychev coefficient 0", 0.1, -1.0, 1.0);
  RooRealVar s2("s2", "Chebychev coefficient 1", 0.0, -1.0, 1.0);

  // --- Background PDF ---
  RooChebychev bkg("bkg", "2nd order Chebychev background",
                   mass, RooArgList(s1, s2));



  // === import dataset ===
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f) && (ctau3D>%.3f && ctau3D < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);
  TString osCut = "(recoQQsign == 0)";
  TString fullCut = Form("%s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data());

  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  auto data_tmp = (RooDataSet *)ds->reduce(Cut(fullCut));
  RooArgSet obs(mass, ctau3D);
  auto data1 = new RooDataSet("data1", "dataset with local vars", obs, Import(*data_tmp));

  // =======================
  // ===     ctau3D      ===
  // =======================
  // ================= Resolution (3-Gauss) =================
  // RooRealVar mu0_ctau("mu0_ctau", "mean ctau", 0.0, -0.01, 0.01);
  RooRealVar mu0_ctau("mu0_ctau", "mean ctau", 0.0);

  RooRealVar sigma1_ctau("sigma1_ctau", "smallest sigma", 0.03, 0.001, 0.1);

  // sigma2 = sigma1 * r2, sigma3 = sigma2 * r3
  // r2, r3 > 1 이므로 항상 sigma1 < sigma2 < sigma3
  RooRealVar r2_ctau("r2_ctau", "scale factor 2", 2.0, 1.01, 5.0);
  RooRealVar r3_ctau("r3_ctau", "scale factor 3", 1.5, 1.01, 5.0);

  RooFormulaVar sigma2_ctau("sigma2_ctau", "@0*@1", RooArgList(sigma1_ctau, r2_ctau));
  RooFormulaVar sigma3_ctau("sigma3_ctau", "@0*@1", RooArgList(sigma2_ctau, r3_ctau));

  RooRealVar f1_ctau("f1_ctau", "frac G1 ctau", 0.6, 0.2, 1.0);
  RooRealVar f2_ctau("f2_ctau", "frac G2 ctau", 0.01, 0.0001, 0.1);

  RooGaussModel g1_ctau("g1_ctau", "Gauss1 ctau", ctau3D, mu0_ctau, sigma1_ctau);
  RooGaussModel g2_ctau("g2_ctau", "Gauss2 ctau", ctau3D, mu0_ctau, sigma2_ctau);
  RooGaussModel g3_ctau("g3_ctau", "Gauss3 ctau", ctau3D, mu0_ctau, sigma3_ctau);

  RooAddModel resModel_ctau("resModel_ctau", "3-Gauss Resolution ctau",
                            RooArgList(g1_ctau, g2_ctau, g3_ctau),
                            RooArgList(f1_ctau, f2_ctau));

  // ================= Nonprompt (Decay ⊗ Res) =================
  RooRealVar lambdaNP_ctau("lambdaNP_ctau", "nonprompt decay const", 1.0, 0.01, 5.0);

  RooDecay decayNP_ctau("decayNP_ctau", "nonprompt decay",
                        ctau3D, lambdaNP_ctau, resModel_ctau, RooDecay::SingleSided);

  // ================= Signal = Prompt + Nonprompt =================
  RooRealVar bFraction("bFraction", "fraction prompt", 0.5, 0.4, 0.8);
  RooFormulaVar fPrSig("fPrSig", "1.0-@0", RooArgList(bFraction));

  // 항상 "res 먼저"로 고정
  RooAddPdf sigModel_ctau("sigModel_ctau", "Signal (Prompt+Nonprompt)",
                          RooArgList(resModel_ctau, decayNP_ctau),
                          RooArgList(fPrSig));

  // RooRealVar bFracion("bFracion", "fraction prompt", 0.6, 0.0, 1.0);
  // RooAddPdf sigModel_ctau("sigModel_ctau", "Signal (Prompt+Nonprompt)",
  //                         RooArgList(resModel_ctau, decayNP_ctau),
  //                         RooArgList(fPrompt_ctau), kTRUE);



  // ================= Background =================
  RooRealVar lambdaL_ctau("lambdaL_ctau", "lambda left", 0.1, 0.01, 10.0);
  RooRealVar lambdaR_ctau("lambdaR_ctau", "lambda right", 1.0, 0.1, 10.0);
  RooRealVar lambdaDS_ctau("lambdaDS_ctau", "lambda double-sided", 1.0, 0.01, 10.0);

  // === Background components (⊗ resolution) ===
  RooDecay bkgL_ctau("bkgL_ctau", "left decay",
                     ctau3D, lambdaL_ctau, resModel_ctau, RooDecay::Flipped);

  RooDecay bkgR_ctau("bkgR_ctau", "right decay",
                     ctau3D, lambdaR_ctau, resModel_ctau, RooDecay::SingleSided);

  RooDecay bkgDS_ctau("bkgDS_ctau", "double sided decay",
                      ctau3D, lambdaDS_ctau, resModel_ctau, RooDecay::DoubleSided);

  // === Fractions ===
  RooRealVar fracRes("fracRes", "fraction Res", 0.1, 0.01, 1.0);
  RooRealVar fracMid("fracMid", "fraction Res", 0.2, 0.01, 1.0);
  RooRealVar fracRight("fracRight", "fraction Res", 0.5, 0.01, 1.0);
  RooRealVar fracLeft("fracLeft", "fraction Res", 0.1, 0.01, 1.0);
  // RooRealVar fRes("fRes", "fraction Res", 0.5, 0.01, 1.0);

  // === background model ===
  RooAddPdf bkgModel_ctau("bkgModel_ctau", "Background ctau",
                          RooArgList(resModel_ctau, bkgDS_ctau, bkgL_ctau, bkgR_ctau),
                          RooArgList(fracRes, fracMid, fracLeft));

  // ======================================
  // 2D fit
  // ======================================
  // ===== 2D observables =====
  RooArgSet obs2D(mass, ctau3D);

  // ===== 2D Signal (factorised) =====
  RooProdPdf sig2D("sig2D", "2D Signal",
                   RooArgList(signal, sigModel_ctau));
  // RooProdPdf sig2D("sig2D", "2D Signal",
  //                  RooArgList(signal_binned, sig_ctau_binned),
  //                  RooFit::Conditional(RooArgSet(signal_binned), RooArgSet(mass)),
  //                  RooFit::Conditional(RooArgSet(sig_ctau_binned), RooArgSet(ctau3D)));

  // ===== 2D Background (factorised) =====
  RooProdPdf bkg2D("bkg2D", "2D Background",
                   RooArgList(bkg, bkgModel_ctau));
  // RooProdPdf bkg2D("bkg2D", "2D Background",
  //                  RooArgList(bkg_binned, bkg_ctau_binned),
  //                  RooFit::Conditional(RooArgSet(bkg_binned), RooArgSet(mass)),
  //                  RooFit::Conditional(RooArgSet(bkg_ctau_binned), RooArgSet(ctau3D)));

  // ===== Total 2D model =====
  RooRealVar Nsig2D("Nsig2D", "yield signal 2D", 1000000, 800000, 1500000);
  RooRealVar Nbkg2D("Nbkg2D", "yield background 2D", 150000, 100000, 300000);

  RooAddPdf model2D("model2D", "Signal + Background (2D)",
                    RooArgList(sig2D, bkg2D),
                    RooArgList(Nsig2D, Nbkg2D));

  // ===== 2D dataset =====
  RooDataHist data2D("data2D", "2D dataset",
                     RooArgList(mass, ctau3D), *data1);

  // ===== Fit =====
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // setFitRangeFromData(data2D, ctau3D, "fitRange");

  auto fitResult2D = model2D.fitTo(data2D,
                                   Save(), Hesse(), Strategy(2),
                                   IntegrateBins(1e-6),
                                   Offset("bins"),
                                   Extended(kTRUE),
                                   PrintLevel(-1), NumCPU(32),
                                   EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(1.),
                                   RooFit::PrintEvalErrors(-1),
                                   RooFit::PrintLevel(-1));

  // auto fitResult2D = model2D.fitTo(*data2D,
  //                                  Save(),
  //                                  IntegrateBins(1e-3),
  //                                  Offset("bin"),
  //                                  Extended(kTRUE),
  //                                  PrintLevel(-1), NumCPU(32),
  //                                  EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(2.),
  //                                  RooFit::PrintEvalErrors(-1),
  //                                  RooFit::PrintLevel(-1));

  // ===== Projection plots =====
  // mass projection
  // === Mass projection with pull ===
  RooPlot *frame_mass2D = mass.frame();
  data2D.plotOn(frame_mass2D, Name("data_mass"));
  model2D.plotOn(frame_mass2D, Name("model_mass"));
  model2D.plotOn(frame_mass2D, Components("sig2D"), LineStyle(kDotted), LineColor(kRed));
  model2D.plotOn(frame_mass2D, Components("bkg2D"), LineStyle(kDotted), LineColor(kGreen));

  int nParMass = fitResult2D->floatParsFinal().getSize();
  double chi2_mass = frame_mass2D->chiSquare("model_mass", "data_mass", nParMass);

  // --- pull hist ---
  RooHist *hpull_mass = frame_mass2D->pullHist("data_mass", "model_mass");
  RooPlot *frame_pull_mass = mass.frame();
  frame_pull_mass->addPlotable(hpull_mass, "P");

  // === Ctau projection with pull ===
  RooPlot *frame_ctau2D = ctau3D.frame();
  data2D.plotOn(frame_ctau2D, Name("data_ctau"));
  model2D.plotOn(frame_ctau2D, Name("model_ctau"));
  model2D.plotOn(frame_ctau2D, Components("sig2D"), LineStyle(kDotted) , LineColor(kRed));
  model2D.plotOn(frame_ctau2D, Components("bkg2D"), LineStyle(kDotted), LineColor(kGreen));

  int nParCtau = fitResult2D->floatParsFinal().getSize();
  double chi2_ctau2d = frame_ctau2D->chiSquare("model_ctau", "data_ctau", nParCtau);

  // --- pull hist ---
  RooHist *hpull_ctau2d = frame_ctau2D->pullHist("data_ctau", "model_ctau");
  RooPlot *frame_pull_ctau2d = ctau3D.frame();
  frame_pull_ctau2d->addPlotable(hpull_ctau2d, "P");
  frame_pull_ctau2d->GetYaxis()->SetRangeUser(-5, 5);

  // --- mass plot ---
  TCanvas c2D("c2D", "2D fit projections with pulls", 1200, 800);
  c2D.Divide(2, 2);

  c2D.cd(1);
  gPad->SetPad(0, 0.5, 0.5, 1.0);
  frame_mass2D->Draw();
  TLatex latex_mass;
  latex_mass.SetNDC();
  latex_mass.SetTextSize(0.06);
  latex_mass.DrawLatex(0.6, 0.85, Form("#chi^{2}/ndf = %.2f", chi2_mass));

  // --- mass pull ---
  c2D.cd(3);
  gPad->SetPad(0, 0, 0.5, 0.5);
  frame_pull_mass->Draw();
  TLine *line0_mass = new TLine(mass.getMin(), 0, mass.getMax(), 0);
  line0_mass->SetLineStyle(2);
  line0_mass->Draw("same");

  // --- ctau plot ---
  c2D.cd(2);
  gPad->SetPad(0.5, 0.5, 1.0, 1.0);
  gPad->SetLogy();
  frame_ctau2D->Draw();

  TLatex latex_ctau2d;
  latex_ctau2d.SetNDC();
  latex_ctau2d.SetTextSize(0.06);
  latex_ctau2d.DrawLatex(0.6, 0.85, Form("#chi^{2}/ndf = %.2f", chi2_ctau2d));

  double nsigVal = Nsig2D.getVal();
  double nsigErr = Nsig2D.getError();
  double nbkgVal = Nbkg2D.getVal();
  double nbkgErr = Nbkg2D.getError();
  double bFracVal = bFraction.getVal();
  double bFracErr = bFraction.getError();

  latex_ctau2d.DrawLatex(0.60, 0.70, Form("N_{sig} = %.0f #pm %.0f", nsigVal, nsigErr));
  latex_ctau2d.DrawLatex(0.60, 0.65, Form("N_{bkg} = %.0f #pm %.0f", nbkgVal, nbkgErr));
  latex_ctau2d.DrawLatex(0.60, 0.60, Form("f_{NP} = %.3f #pm %.3f", bFracVal, bFracErr));

  // --- ctau pull ---
  c2D.cd(4);
  gPad->SetPad(0.5, 0, 1.0, 0.5);
  frame_pull_ctau2d->Draw();
  TLine *line0_ctau = new TLine(ctau3D.getMin(), 0, ctau3D.getMax(), 0);
  line0_ctau->SetLineStyle(2);
  line0_ctau->Draw("same");

  c2D.SaveAs(Form("fit2D_projections_with_pulls_pT%.1f.png", ptLow));
  fitResult2D->Print("V");


  cout
      << "=== finish final_2d_fit()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}