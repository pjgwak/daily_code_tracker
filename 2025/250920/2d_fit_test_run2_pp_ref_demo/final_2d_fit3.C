#include <gsl/gsl_errno.h>
#include <TStopwatch.h>
#include <RooRealVar.h>
#include "RooFit.h"
#include "RooMinimizer.h"

using namespace RooFit;

void final_2d_fit2()
{
  cout << "=== start final_2d_fit() ===\n";
  TStopwatch t;
  t.Start();

  double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
  gSystem->mkdir("figs", true);

  // === mass ===
  RooRealVar mass("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
  int nMassBins = 90;

  // keep ctau bin width about 83 micro-meter
  double ctMin = -1, ctMax = 4;
  RooRealVar ctau3D("ctau3D", "ctau3D", ctMin, ctMax, "mm");
  int nCtauBins = 180;

  // --- 공통 평균 ---
  RooRealVar mean("mean", "mean of CB", 3.1, 3.0, 3.2);

  // === CB1 ===
  RooRealVar sigma1("sigma1", "CB1 sigma", 0.03, 0.005, 0.1);
  RooRealVar alpha1("alpha1", "CB1 alpha", 1.5, 0.1, 10.0);
  RooRealVar n1("n1", "CB1 n", 3.0, 1.0, 10.0);

  RooCBShape cb1("cb1", "Crystal Ball 1",
                 mass, mean, sigma1, alpha1, n1);

  // === CB2 ===
  RooRealVar sigma2("sigma2", "CB2 sigma", 0.02, 0.005, 0.1);
  RooRealVar alpha2("alpha2", "CB2 alpha", -1.5, -10.0, -0.1);
  RooRealVar n2("n2", "CB2 n", 3.0, 1.0, 10.0);

  RooCBShape cb2("cb2", "Crystal Ball 2",
                 mass, mean, sigma2, alpha2, n2);

  // === Gaussian ===
  RooRealVar sigmaG("sigmaG", "Gaussian sigma", 0.015, 0.005, 0.4);
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

  // === mass - Binned PDF ===
  RooBinSamplingPdf signal_binned(
      "signal_binned", "Binned Double Crystal Ball",
      mass,
      signal,
      nMassBins);

  RooBinSamplingPdf bkg_binned(
      "bkg_binned", "Binned Chebychev2",
      mass,
      bkg,
      nMassBins);

  // combine for test mass fit
  RooRealVar nsig("nsig", "signal yield", 800000, 600000, 1e6);
  RooRealVar nbkg("nbkg", "background yield", 160000, 100000, 1e6);
  RooAddPdf model_binned(
      "model_binned", "sig+bkg binned model",
      RooArgList(signal_binned, bkg_binned),
      RooArgList(nsig, nbkg));

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

  // RooBinning massBinning(nMassBins, massLow, massHigh); // 90 bins, 2.6~3.5 GeV
  auto dataHist = new RooDataHist("dataHist", "binned dataset",
                                  RooArgSet(mass), *data1);

  // binned model fit
  auto fitResult = model_binned.fitTo(
      *dataHist,
      Save(),
      Extended(kTRUE),
      PrintLevel(-1));

  fitResult->Print("v");

  // ===== 3. 메인 plot =====
  RooPlot *frame = mass.frame();
  dataHist->plotOn(frame, Name("dataHist"));
  model_binned.plotOn(frame, Name("model_binned"));
  model_binned.plotOn(frame, Components("signal_binned"), LineColor(kRed));
  model_binned.plotOn(frame, Components("bkg_binned"), LineColor(kBlue));

  // ===== 4. Pull plot =====
  RooHist *hpull = frame->pullHist("dataHist", "model_binned"); // (data - fit)/error
  RooPlot *frame_pull = mass.frame();
  frame_pull->addPlotable(hpull, "P"); // pull 점 찍기

  // ===== 5. 캔버스 구성 =====
  TCanvas c("c", "fit with pull", 800, 800);
  c.Divide(1, 2);

  // 상단: 메인 fit
  c.cd(1);
  gPad->SetPad(0, 0.3, 1, 1.0); // 상단 70%
  frame->Draw();

  // 하단: pull
  c.cd(2);
  gPad->SetPad(0, 0, 1, 0.3);                  // 하단 30%
  frame_pull->GetYaxis()->SetRangeUser(-5, 5); // pull 축
  frame_pull->Draw();

  double chi2 = frame->chiSquare("model_binned", "dataHist",
                                      fitResult->floatParsFinal().getSize());
  TLatex latex;
  latex.SetNDC(); // Normalized coordinates
  latex.SetTextSize(0.08);
  latex.DrawLatex(0.72, 0.81, Form("#chi^{2}/ndf = %.2f", chi2));

  // ===== 6. 결과 저장 =====
  // c.SaveAs("fit_mass_with_pull.pdf");
  c.SaveAs("fit_mass_with_pull.png");

  // =======================
  // ===     ctau3D      ===
  // =======================
  // ================= Resolution (3-Gauss) =================
  // RooRealVar mu0_ctau("mu0_ctau", "mean ctau", 0.0, -0.01, 0.01);
  RooRealVar mu0_ctau("mu0_ctau", "mean ctau", 0.0);

  RooRealVar sigma1_ctau("sigma1_ctau", "smallest sigma", 0.03, 0.005, 0.1);

  // sigma2 = sigma1 * r2, sigma3 = sigma2 * r3
  // r2, r3 > 1 이므로 항상 sigma1 < sigma2 < sigma3
  RooRealVar r2_ctau("r2_ctau", "scale factor 2", 2.0, 1.01, 5.0);
  RooRealVar r3_ctau("r3_ctau", "scale factor 3", 1.5, 1.01, 5.0);

  RooFormulaVar sigma2_ctau("sigma2_ctau", "@0*@1", RooArgList(sigma1_ctau, r2_ctau));
  RooFormulaVar sigma3_ctau("sigma3_ctau", "@0*@1", RooArgList(sigma2_ctau, r3_ctau));

  RooRealVar f1_ctau("f1_ctau", "frac G1 ctau", 0.6, 0.5, 1.0);
  RooRealVar f2_ctau("f2_ctau", "frac G2 ctau", 0.3, 0.0, 1.0);

  RooGaussModel g1_ctau("g1_ctau", "Gauss1 ctau", ctau3D, mu0_ctau, sigma1_ctau);
  RooGaussModel g2_ctau("g2_ctau", "Gauss2 ctau", ctau3D, mu0_ctau, sigma2_ctau);
  RooGaussModel g3_ctau("g3_ctau", "Gauss3 ctau", ctau3D, mu0_ctau, sigma3_ctau);

  RooAddModel resModel_ctau("resModel_ctau", "3-Gauss Resolution ctau",
                            RooArgList(g1_ctau, g2_ctau, g3_ctau),
                            RooArgList(f1_ctau, f2_ctau));

  // ================= Nonprompt (Decay ⊗ Res) =================
  RooRealVar lambdaNP_ctau("lambdaNP_ctau", "nonprompt decay const", 1.0, 0.2, 5.0);

  RooDecay decayNP_ctau("decayNP_ctau", "nonprompt decay",
                        ctau3D, lambdaNP_ctau, resModel_ctau, RooDecay::SingleSided);

  // ================= Signal = Prompt + Nonprompt =================
  RooRealVar bFraction("bFraction", "fraction prompt", 0.6, 0.0, 1.0);
  RooFormulaVar fPrSig("fPrSig", "1.0-@0", RooArgList(bFraction));

  // 항상 "res 먼저"로 고정
  RooAddPdf sigModel_ctau("sigModel_ctau", "Signal (Prompt+Nonprompt)",
                          RooArgList(resModel_ctau, decayNP_ctau),
                          RooArgList(fPrSig, bFraction));

  // RooRealVar bFracion("bFracion", "fraction prompt", 0.6, 0.0, 1.0);
  // RooAddPdf sigModel_ctau("sigModel_ctau", "Signal (Prompt+Nonprompt)",
  //                         RooArgList(resModel_ctau, decayNP_ctau),
  //                         RooArgList(fPrompt_ctau), kTRUE);

  // --- wrap with RooBinSamplingPdf ---
  RooBinSamplingPdf sig_ctau_binned("sig_ctau_binned", "binned signal ctau",
                                    ctau3D, sigModel_ctau, nCtauBins);

  // ================= Background =================
  RooRealVar lambdaL_ctau("lambdaL_ctau", "lambda left", 0.01, 0.0001, 100.0);
  RooRealVar lambdaR_ctau("lambdaR_ctau", "lambda right", 1.0, 0.1, 5.0);
  RooRealVar lambdaDS_ctau("lambdaDS_ctau", "lambda double-sided", 1.0, 0.1, 5.0);

  // === Background components (⊗ resolution) ===
  RooDecay bkgL_ctau("bkgL_ctau", "left decay",
                     ctau3D, lambdaL_ctau, resModel_ctau, RooDecay::Flipped);

  RooDecay bkgR_ctau("bkgR_ctau", "right decay",
                     ctau3D, lambdaR_ctau, resModel_ctau, RooDecay::SingleSided);

  RooDecay bkgDS_ctau("bkgDS_ctau", "double sided decay",
                      ctau3D, lambdaDS_ctau, resModel_ctau, RooDecay::DoubleSided);

  // === Fractions ===
  RooRealVar fRes("fRes", "fraction Res", 0.2, 0.0, 1.0); // resolution fraction
  RooRealVar alpha("alpha", "fraction Mid within non-Res", 0.3, 0.0, 1.0);
  RooRealVar beta("beta", "fraction Right within non-Res & non-Mid", 0.3, 0.0, 1.0);

  // === derived fractions ===
  // Res fraction
  RooFormulaVar fracRes("fracRes", "@0", RooArgList(fRes));

  // Mid fraction = (1-fRes) * alpha
  RooFormulaVar fracMid("fracMid", "(1-@0)*@1", RooArgList(fRes, alpha));

  // Right fraction = (1-fRes) * (1-alpha) * beta
  RooFormulaVar fracRight("fracRight", "(1-@0)*(1-@1)*@2",
                          RooArgList(fRes, alpha, beta));

  // Left fraction = (1-fRes) * (1-alpha) * (1-beta)
  RooFormulaVar fracLeft("fracLeft", "(1-@0)*(1-@1)*(1-@2)",
                         RooArgList(fRes, alpha, beta));

  // === background model ===
  RooAddPdf bkgModel_ctau("bkgModel_ctau", "Background ctau",
                          RooArgList(resModel_ctau, bkgDS_ctau, bkgR_ctau, bkgL_ctau),
                          RooArgList(fracRes, fracMid, fracRight, fracLeft));

  // --- wrap with RooBinSamplingPdf ---
  RooBinSamplingPdf bkg_ctau_binned("bkg_ctau_binned", "binned bkg ctau",
                                    ctau3D, bkgModel_ctau, nCtauBins);

  // ================= Total model =================
  RooRealVar nsig_ctau("nsig_ctau", "signal yield ctau", 5000, 0, 1e6);
  RooRealVar nbkg_ctau("nbkg_ctau", "bkg yield ctau", 2000, 0, 1e6);

  RooAddPdf totalModel_ctau("totalModel_ctau", "Signal+Bkg ctau",
                            RooArgList(sig_ctau_binned, bkg_ctau_binned),
                            RooArgList(nsig_ctau, nbkg_ctau));

  // ===== 6. Binned data =====
  // ctau3D.setBins(80, "fitRange");
  RooDataHist dataHist_ctau("dataHist_ctau", "binned dataset (ctau3D)",
                            RooArgSet(ctau3D), *data1);
  // ===== 7. Fit =====
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  gsl_set_error_handler_off();
  auto fitResult_ctau = totalModel_ctau.fitTo(dataHist_ctau,
                                                Range(ctMin,ctMax),
                                                Save(),
                                                Extended(kTRUE),
                                                RooFit::RecoverFromUndefinedRegions(1.),
                                                RooFit::PrintEvalErrors(-1),
                                                RooFit::PrintLevel(-1));

  // ===== 8. Plot + Pull =====
  RooPlot *frame_ctau = ctau3D.frame();
  dataHist_ctau.plotOn(frame_ctau, Name("data_ctau"));
  totalModel_ctau.plotOn(frame_ctau, Name("model_ctau"));
  totalModel_ctau.plotOn(frame_ctau, Components("sig_ctau_binned"), LineColor(kRed));
  totalModel_ctau.plotOn(frame_ctau, Components("bkg_ctau_binned"), LineColor(kGreen));

  // pull
  RooHist *hpull_ctau = frame_ctau->pullHist("data_ctau", "model_ctau");
  RooPlot *frame_pull_ctau = ctau3D.frame();
  frame_pull_ctau->addPlotable(hpull_ctau, "P");

  // χ² 계산
  int npar_ctau = fitResult_ctau->floatParsFinal().getSize();
  double chi2_ctau = frame_ctau->chiSquare("model_ctau", "data_ctau", npar_ctau);

  // ===== 9. Canvas =====
  TCanvas c_ctau("c_ctau", "ctau3D fit with pull", 800, 800);
  c_ctau.Divide(1, 2);

  c_ctau.cd(1);
  gPad->SetPad(0, 0.3, 1, 1.0);
  gPad->SetLogy();
  frame_ctau->Draw();
  TLatex latex_ctau;
  latex_ctau.SetNDC();
  latex_ctau.SetTextSize(0.04);
  latex_ctau.DrawLatex(0.6, 0.85, Form("#chi^{2}/ndf = %.2f", chi2_ctau));

  c_ctau.cd(2);
  gPad->SetPad(0, 0, 1, 0.3);
  frame_pull_ctau->GetYaxis()->SetRangeUser(-5, 5);
  frame_pull_ctau->Draw();

  c_ctau.SaveAs("fit_ctau3D_with_pull.png");
  fitResult_ctau->Print("V");

  // ======================================
  // 2D fit
  // ======================================
  // ===== 2D observables =====
  RooArgSet obs2D(mass, ctau3D);

  // ===== 2D Signal (factorised) =====
  RooProdPdf sig2D("sig2D", "2D Signal",
                   RooArgList(signal, sig_ctau_binned));
  // RooProdPdf sig2D("sig2D", "2D Signal",
  //                  RooArgList(signal, sig_ctau_binned),
  //                  RooFit::Conditional(RooArgSet(signal), RooArgSet(mass)),
  //                  RooFit::Conditional(RooArgSet(sig_ctau_binned), RooArgSet(ctau3D)));

  // ===== 2D Background (factorised) =====
  RooProdPdf bkg2D("bkg2D", "2D Background",
                   RooArgList(bkg, bkg_ctau_binned));
  // RooProdPdf bkg2D("bkg2D", "2D Background",
  //                  RooArgList(bkg, bkg_ctau_binned),
  //                  RooFit::Conditional(RooArgSet(bkg), RooArgSet(mass)),
  //                  RooFit::Conditional(RooArgSet(bkg_ctau_binned), RooArgSet(ctau3D)));

  // ===== Total 2D model =====
  RooRealVar Nsig2D("Nsig2D", "yield signal 2D", 100000, 0, 1e6);
  RooRealVar Nbkg2D("Nbkg2D", "yield background 2D", 50000, 0, 1e6);

  RooAddPdf model2D("model2D", "Signal + Background (2D)",
                    RooArgList(sig2D, bkg2D),
                    RooArgList(Nsig2D, Nbkg2D));

  // ===== 2D dataset =====
  RooDataHist data2D("data2D", "2D dataset",
                     RooArgList(mass, ctau3D), *data1);

  // ===== Fit =====
  auto fitResult2D = model2D.fitTo(data2D,
                                   Save(),
                                   Extended(kTRUE),
                                   PrintLevel(-1), NumCPU(32),
                                   EvalBackend("legacy"),RooFit::RecoverFromUndefinedRegions(1.),
                                   RooFit::PrintEvalErrors(-1),
                                   RooFit::PrintLevel(-1));

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

  c2D.SaveAs("fit2D_projections_with_pulls.png");
  fitResult2D->Print("V");

  // // --- mass signal - two gaussians ---
  // double mean0 = 3.0969;
  // double sigma1_0 = 0.020;
  // double sigma2_0 = 0.040;
  // double sigmaG_0 = 0.015;

  // // === bring MC fit results ===
  // TFile fin(Form("roots/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  // RooFitResult *fr = (RooFitResult *)fin.Get("fitResult");
  // if (!fr)
  // {
  //   std::cerr << "fitResult not found!" << std::endl;
  //   return;
  // }

  // const RooArgList &params = fr->floatParsFinal();
  // RooRealVar *sigmaGVar = (RooRealVar *)params.find("sigmaG");
  // RooRealVar *alphaLVar = (RooRealVar *)params.find("alphaL");
  // RooRealVar *nLVar = (RooRealVar *)params.find("nL");
  // RooRealVar *sigmaRVar = (RooRealVar *)params.find("sigmaR");
  // RooRealVar *alphaRatioVar = (RooRealVar *)params.find("alphaRatio");
  // RooRealVar *nRatioVar = (RooRealVar *)params.find("nRatio");
  // RooRealVar *sigmaRatioVar = (RooRealVar *)params.find("sigmaRatio");
  // RooRealVar *sigmaGRatioVar = (RooRealVar *)params.find("sigmaGRatio");

  // // if (sigmaL)
  // // {
  // //   sigmaL->setConstant(kTRUE);
  // //   std::cout << "Fixed sigmaL = " << sigmaL->getVal() << std::endl;
  // // }

  // // if (alphaLVar)
  // // {
  // //   alphaLVar->setConstant(kTRUE);
  // //   std::cout << "\nFixed alphaL = " << alphaLVar->getVal() << std::endl;
  // // }
  // // if (nLVar)
  // // {
  // //   nLVar->setConstant(kTRUE);
  // //   std::cout << "Fixed nL = " << nLVar->getVal() << std::endl;
  // // }
  // // if (sigmaGVar)
  // // {
  // //   sigmaGVar->setConstant(kTRUE);
  // //   std::cout << "Fixed sigmaG = " << sigmaGVar->getVal() << std::endl;
  // // }
  // // if (alphaRatioVar)
  // // {
  // //   alphaRatioVar->setConstant(kTRUE);
  // //   std::cout << "Fixed alphaRatio = " << alphaRatioVar->getVal() << std::endl;
  // // }
  // // if (nRatioVar)
  // // {
  // //   nRatioVar->setConstant(kTRUE);
  // //   std::cout << "Fixed nRatio = " << nRatioVar->getVal() << std::endl;
  // // }
  // // if (sigmaRatioVar)
  // // {
  // //   sigmaRatioVar->setConstant(kTRUE);
  // //   std::cout << "Fixed sigmaRatio = " << sigmaRatioVar->getVal() << std::endl;
  // // }
  // // if (sigmaGRatioVar)
  // // {
  // //   sigmaGRatioVar->setConstant(kTRUE);
  // //   std::cout << "Fixed sigmaGRatio = " << sigmaGRatioVar->getVal() << std::endl;
  // // }

  // // --- mass fit result ---
  // TFile fInMass(Form("roots/mass_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  // RooFitResult *fRMass = (RooFitResult *)fInMass.Get("fitResult");
  // if (!fRMass)
  // {
  //   std::cerr << "fitResult not found!" << std::endl;
  //   return;
  // }

  // const RooArgList &paraMass = fRMass->floatParsFinal();
  // RooRealVar *s1 = (RooRealVar *)paraMass.find("c1");
  // RooRealVar *f_cbR = (RooRealVar *)paraMass.find("f_cbR");
  // RooRealVar *f_gaus = (RooRealVar *)paraMass.find("f_gaus");
  // RooRealVar *meanMass = (RooRealVar *)paraMass.find("mean");
  // RooRealVar *sigmaL = (RooRealVar *)paraMass.find("sigmaL");

  // // if (s1)
  // // {
  // //   s1->setConstant(kTRUE);
  // //   std::cout << "Fixed s1 = " << s1->getVal() << std::endl;
  // // }
  // // if (f_cbR)
  // // {
  // //   f_cbR->setConstant(kTRUE);
  // //   std::cout << "Fixed f_cbR = " << f_cbR->getVal() << std::endl;
  // // }

  // // if (f_gaus)
  // // {
  // //   f_gaus->setConstant(kTRUE);
  // //   std::cout << "Fixed f_gaus = " << f_gaus->getVal() << std::endl;
  // // }

  // // if (meanMass)
  // // {
  // //   meanMass->setConstant(kTRUE);
  // //   std::cout << "Fixed meanMass = " << meanMass->getVal() << std::endl;
  // // }
  // // if (sigmaL)
  // // {
  // //   sigmaL->setConstant(kTRUE);
  // //   std::cout << "Fixed sigmaL = " << sigmaL->getVal() << std::endl;
  // // }

  // // === build model ===
  // RooCBShape cbLeft("cbLeft", "CB left tail", *mass, *meanMass, *sigmaL, *alphaLVar, *nLVar);

  // // --- CB right ---
  // RooFormulaVar sigmaR("sigmaR", "sigma right", "@0*@1", RooArgList(*sigmaL, *sigmaRatioVar));
  // RooFormulaVar alphaR("alphaR", "alpha right", "-@0*@1", RooArgList(*alphaLVar, *alphaRatioVar));
  // RooFormulaVar nR("nR", "n right", "@0*@1", RooArgList(*nLVar, *nRatioVar));
  // RooCBShape cbRight("cbRight", "CB right tail", *mass, *meanMass, sigmaR, alphaR, nR);

  // // --- gauss ---
  // // RooRealVar sigmaGRatio("sigmaGRatio", "", 1.0, 0.5, 10.0);
  // RooFormulaVar sigmaG("sigmaG", "sigma gauss", "@0*@1", RooArgList(*sigmaL, *sigmaGRatioVar));
  // RooGaussian gaus("gaus", "gaus core", *mass, *meanMass, sigmaG);

  // // --- combine sig components ---
  // // RooRealVar f_cbR("f_cbR", "frac of right CB", 0.17, 0, 1.0);
  // // RooRealVar f_gaus("f_gaus", "frac of gaus", 0.45, 0, 1);
  // auto mass_b_pdf = new RooAddPdf("mass_b_pdf", "", RooArgList(cbLeft, cbRight, gaus), RooArgList(*f_cbR, *f_gaus), kTRUE);

  // // extend
  // RooRealVar *N_np = new RooRealVar("N_{np}", "signal number B Jpsi K", 400000, 1000, 600000);
  // RooAbsPdf *mass_b_epdf = new RooExtendPdf("mass_b_epdf", "mass sig B 2 Jpsi K", *mass_b_pdf, *N_np);

  // // --- ctau signal ---
  // RooRealVar *ctau3D = new RooRealVar("ctau3D", "proper decay length", -1, 4, "mm");
  // RooRealVar *b_tau = new RooRealVar("c#tau_{b}", "B ctau", 0.5, 0.001, 1, "mm");
  // RooRealVar *Mean_tom = new RooRealVar("tomMean", "proper time mean", 0);
  // RooRealVar *ratio_tom = new RooRealVar("s_{to}", "ratio", 1.5, 1, 50);  // ratio of what?
  // RooRealVar *ratio_tom2 = new RooRealVar("s_{to}", "ratio", 1.5, 1, 50); // ratio of what?
  // RooRealVar *ctau3D_err = new RooRealVar("tom_err", "tom_err", 0, 1.0, "mm");

  // RooRealVar *Sigma0_tom = new RooRealVar("#sigma_{t0}", "proper time sigma", 0.05, 0.005, 0.2);
  // // per event err??
  // // RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *ctau3D_err));
  // RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *Sigma0_tom));
  // RooFormulaVar *Sigma2_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom2, *Sigma1_tom));

  // RooRealVar *frac_g_tom = new RooRealVar("f_g_#sigma_{t0}", "gauss resol. fracion", 0.8048, 0.01, 1);
  // RooRealVar *frac_g_tom2 = new RooRealVar("f_g_#sigma2", "gauss resol. fracion", 0.8048, 0.01, 1);

  // RooResolutionModel *g0_tom = new RooGaussModel("g0_tom", "proper time resolution", *ctau3D, *Mean_tom, *Sigma0_tom);
  // RooResolutionModel *g1_tom = new RooGaussModel("g1_tom", "proper time resolution", *ctau3D, *Mean_tom, *Sigma1_tom);
  // RooResolutionModel *g2_tom = new RooGaussModel("g2_tom", "proper time resolution", *ctau3D, *Mean_tom, *Sigma2_tom);

  // // one gauss signal
  // // RooResolutionModel *properTimeRes = new RooGaussModel("properTimeRes", "proper time resolution", *ctau3D, *Mean_tom, *Sigma1_tom);

  // // two gauss signal
  // // RooResolutionModel *properTimeRes = new RooAddModel("properTimeRes", "dual gaussain PDF", RooArgList(*g0_tom, *g1_tom), RooArgList(*frac_g_tom));

  // // three gauss signal
  // RooResolutionModel *properTimeRes = new RooAddModel("properTimeRes", "dual gaussain PDF", RooArgList(*g0_tom, *g1_tom, *g2_tom), RooArgList(*frac_g_tom, *frac_g_tom2));

  // RooAbsPdf *time_b_pdf = new RooDecay("time_b_pdf", "time pdf B2 Jpsi K", *ctau3D, *b_tau, *properTimeRes, RooDecay::SingleSided);
  // // OK, signal time PDF is ready

  // // --- combine signal PDFs (mass + time)
  // RooAbsPdf *masstimeSigpdf = new RooProdPdf("masstimeSigpdf", "signal m* t", RooArgSet(*time_b_pdf, *mass_b_epdf));

  // // --- NP ---
  // RooRealVar *N_pr = new RooRealVar("N_{PR}", "signal number B Jpsi K", 800000, 500000, 1000000);
  // RooAbsPdf *mass_np_epdf = new RooExtendPdf("mass_np_epdf", "mass sig B 2 Jpsi K", *mass_b_pdf, *N_pr);
  // RooAbsPdf *masstimeSigpdf2 = new RooProdPdf("masstimeSigpdf2", "signal m*t", RooArgSet(*properTimeRes, *mass_np_epdf));

  // // === Jpsi bg ===
  // // BB mass -> mass와 slope는 공유
  // // RooRealVar *slope = new RooRealVar("slope", "slope", 0.01, -1, 1);
  // // RooAbsPdf *mass_BB_pdf = new RooPolynomial("mass_QCD_pdf", "mass pdf QCD", *mass, *slope);
  // RooAbsPdf *mass_BB_pdf = new RooChebychev("mass_QCD_pdf", "mass pdf QCD",
  //                                           *mass, *s1);

  // // prompt BB time - NP인데 Res과 함쳤다는 뜻?
  // // --- base tau
  // // RooRealVar tau0("tau0", "base B ctau", 0.01, 0.001, 0.2, "mm");

  // // RooRealVar tau_left("tau_left", "scale factor left", 1.0, 0.001, 100.0);
  // // RooRealVar tau_mid("tau_mid", "scale factor middle", 3, 0.01, 10.0);
  // // RooRealVar tau_right("tau_right", "scale factor right", 3.0, 0.1, 100.0);

  // // --- 3 decay PDF ---
  // TFile fCBkg(Form("roots/ctau_CtauFull_RSB_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  // RooFitResult *fRCBkg = (RooFitResult *)fCBkg.Get("fitResult");
  // if (!fRCBkg)
  // {
  //   std::cerr << "fitResult not found!" << std::endl;
  //   return;
  // }

  // const RooArgList &paraCBkg = fRCBkg->floatParsFinal();
  // RooRealVar *f_left = (RooRealVar *)paraCBkg.find("f_left");
  // RooRealVar *f_mid = (RooRealVar *)paraCBkg.find("f_mid");
  // RooRealVar *f_right = (RooRealVar *)paraCBkg.find("f_right");
  // RooRealVar *tau_left = (RooRealVar *)paraCBkg.find("tau_left");
  // RooRealVar *tau_mid = (RooRealVar *)paraCBkg.find("tau_mid");
  // RooRealVar *tau_right = (RooRealVar *)paraCBkg.find("tau_right");

  // // if (f_mid)
  // // {
  // //   f_mid->setConstant(kTRUE);
  // //   std::cout << "Fixed f_mid = " << f_mid->getVal() << std::endl;
  // // }
  // // if (f_right)
  // // {
  // //   f_right->setConstant(kTRUE);
  // //   std::cout << "Fixed f_right = " << f_right->getVal() << std::endl;
  // // }
  // // if (f_left)
  // // {
  // //   f_left->setConstant(kTRUE);
  // //   std::cout << "Fixed f_left = " << f_left->getVal() << std::endl;
  // // }
  // // if (tau_left)
  // // {
  // //   tau_left->setConstant(kTRUE);
  // //   std::cout << "Fixed tau_left = " << tau_left->getVal() << std::endl;
  // // }
  // // if (tau_mid)
  // // {
  // //   tau_mid->setConstant(kTRUE);
  // //   std::cout << "Fixed tau_mid = " << tau_mid->getVal() << std::endl;
  // // }
  // // if (tau_right)
  // // {
  // //   tau_right->setConstant(kTRUE);
  // //   std::cout << "Fixed tau_right = " << tau_right->getVal() << std::endl;
  // // }

  // RooDecay time_BB_left("time_BB_left", "Left decay",
  //                       *ctau3D, *tau_left, *properTimeRes, RooDecay::Flipped);

  // RooDecay time_BB_mid("time_BB_mid", "Center decay",
  //                      *ctau3D, *tau_mid, *properTimeRes, RooDecay::DoubleSided);

  // RooDecay time_BB_right("time_BB_right", "Right decay",
  //                        *ctau3D, *tau_right, *properTimeRes, RooDecay::SingleSided);

  // // RooRealVar f_left("f_left", "fraction left", 0.3, 0.0, 1.0);
  // // RooRealVar f_mid("f_mid", "fraction middle", 0.3, 0.0, 1.0);
  // auto time_BB_pdf = new RooAddPdf("time_BB_pdf", "sum of 2 decays",
  //                                  RooArgList(time_BB_left, time_BB_mid, time_BB_right, *properTimeRes),
  //                                  RooArgList(*f_left, *f_mid, *f_right), kTRUE);

  // // RooRealVar f_resol_bkg("f_resol_bkg", "fraction of resolution", 0.1, 0.0, 1.0);
  // // auto time_BB_pdf = new RooAddPdf("time_BB_pdf", "decay+resol combined",
  // //                                  RooArgList(*properTimeRes, *time_BB_pdf_old),
  // //                                  RooArgList(f_resol_bkg));

  // RooRealVar *N_BB = new RooRealVar("N_{BB}", "Signal number BB", 200000, 100000, 300000);
  // RooAbsPdf *mass_BB_epdf = new RooExtendPdf("mass_BB_epdf", "mass BB", *mass_BB_pdf, *N_BB);

  // // bkg BB mass + timie
  // RooAbsPdf *masstimeBkgpdf = new RooProdPdf("masstimeBkgpdf", "background2 m * t", RooArgSet(*time_BB_pdf, *mass_BB_epdf));
  // /// OK bkg2 BB is readg

  // // === Sig + Bkg PDF ===
  // RooAddPdf *pdf_fit = new RooAddPdf(
  //     "pdf_fit", "Total(B+S) PDF",
  //     RooArgSet(*masstimeSigpdf, *masstimeSigpdf2, *masstimeBkgpdf));

  // // === Fit ===
  // /// Real data generate -> Toy MC로 테스트
  // RooArgSet *pdf_obs = new RooArgSet(*mass, *ctau3D);

  // TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  // if (!fInput || fInput->IsZombie())
  // {
  //   cout << "Error: cannot open input file\n";
  //   return;
  // }
  // // read dataset
  // RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  // if (!ds)
  // {
  //   cout << "Error: cannot find RooDataSet\n";
  //   return;
  // }

  // // RooMsgService::instance().getStream(RooFit::ERROR).removeTopic(RooFit::Tracing);

  // // to make sure that 2 variables has proper ranges
  // // origianl RooDataSet could include variables having very wide range
  // RooArgSet obs(*mass, *ctau3D);

  // double ctMin = -2, ctMax = 10;
  // mass->setRange(2.6, 3.5);
  // mass->setRange("fitRange", 2.6, 3.5);
  // ctau3D->setRange(ctMin, ctMax);
  // ctau3D->setRange("fitRange", ctMin, ctMax);

  // mass->setBins(40);
  // ctau3D->setBins(80);

  // TString kineCut = Form( // tmp: no cbin
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f) && (ctau3D>%.3f && ctau3D < %.3f)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);

  // auto data_tmp = (RooDataSet *)ds->reduce(Cut(kineCut));
  // auto data1 = new RooDataSet("data1", "dataset with local vars", obs, Import(*data_tmp));
  // // std::unique_ptr<RooDataSet> data(pdf_fit->generate(RooArgSet(*mass, *ctau3D), 20000));

  // // RooBinSamplingPdf binSampler("binSampler", "", *ctau3D, *fit_pdf);
  // // binSampler.fitTo(data);

  // auto data = new RooDataHist("data", "binned dataset", obs, *data1);
  // auto gausBinned = new RooBinSamplingPdf("gausBinned", "Gauss with fine sampling in each bin", obs, *pdf_fit);

  // // fit

  // // --- remove message ---
  // // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  // // "Eval" (NaN)
  // for (int i = 0; i < 3; i++)
  // {
  //   RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  // }
  // // Tracing
  // for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  // {
  //   RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  // }
  // gErrorIgnoreLevel = kFatal;

  // std::unique_ptr<RooFitResult> result{pdf_fit->fitTo(*data,
  //                                                     RooFit::Save(),
  //                                                     Range("fitRange"),
  //                                                     Extended(true),
  //                                                     NumCPU(32),
  //                                                     EvalBackend("legacy"),
  //                                                     RooFit::RecoverFromUndefinedRegions(1.5),
  //                                                     RooFit::PrintEvalErrors(-1),
  //                                                     RooFit::PrintLevel(-1))};
  // // Define difference -log(L) for B->ff, B->ff fbar and B->fbar f
  // // auto nllB = pdf_fit->createNLL(*data,
  // //                                Extended(true),
  // //                                NumCPU(32),
  // //                                EvalBackend("legacy"),
  // //                                SumW2Error(true), Offset(true), Verbose(kFALSE));

  // // RooMinimizer m1(*nllB);
  // // // m1.RecoverFromUndefinedRegions(1);
  // // m1.setVerbose(kFALSE);
  // // m1.setProfile(0);
  // // m1.setPrintLevel(-1);
  // // m1.setStrategy(2);
  // // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // // // fit
  // // m1.migrad();
  // // m1.hesse();
  // // m1.migrad();
  // // std::unique_ptr<RooFitResult> result(m1.save());

  // {
  //   TCanvas c1("c1", "", 800, 800);
  //   c1.Divide(1, 2);
  //   TPad *pad1 = (TPad *)c1.cd(1);
  //   gPad->SetLogy();
  //   pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  //   pad1->SetBottomMargin(0.02);

  //   TPad *pad2 = (TPad *)c1.cd(2);
  //   pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  //   pad2->SetTopMargin(0.05);
  //   pad2->SetBottomMargin(0.25);

  //   // main plot
  //   pad1->cd();
  //   RooPlot *frame = ctau3D->frame();
    
  //   data->plotOn(frame, Name("data"));
  //   pdf_fit->plotOn(frame, Name("model"));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeSigpdf"),
  //                   Name("np"),
  //                   LineColor(kRed),
  //                   LineStyle(kDashed),
  //                   LineWidth(4));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeSigpdf2"),
  //                   Name("pr"),
  //                   LineColor(kOrange),
  //                   LineStyle(kDashed),
  //                   LineWidth(4));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeBkgpdf"),
  //                   Name("bkg"),
  //                   LineColor(kGreen + 2),
  //                   LineStyle(kDashed),
  //                   LineWidth(4));
  //   frame->SetMinimum(0.01);
  //   // frame->SetMaximum(ymax * 10.0);
  //   frame->Draw();

  //   // --- legend ---
  //   auto leg = new TLegend(0.6, 0.65, 0.95, 0.88);
  //   leg->SetBorderSize(0);
  //   leg->SetFillStyle(0);
  //   leg->AddEntry(frame->findObject("data"), "data", "pe");
  //   leg->AddEntry(frame->findObject("model"), "Total fit", "e");
  //   leg->AddEntry(frame->findObject("pr"), "Prompt J/#psi", "e");
  //   leg->AddEntry(frame->findObject("np"), "NP J/#psi", "e");
  //   leg->AddEntry(frame->findObject("bkg"), "Bkg", "e");
  //   leg->Draw();

  //   // pull
  //   pad2->cd();
  //   RooHist *hpull = frame->pullHist("data", "model"); // data - fit / error
  //   RooPlot *frame_pull = ctau3D->frame();
  //   frame_pull->addPlotable(hpull, "P");
  //   frame_pull->SetTitle("");
  //   frame_pull->GetYaxis()->SetTitle("Pull");
  //   frame_pull->GetYaxis()->SetNdivisions(505);
  //   frame_pull->GetYaxis()->SetTitleSize(0.08);
  //   frame_pull->GetYaxis()->SetLabelSize(0.08);
  //   frame_pull->GetXaxis()->SetTitleSize(0.1);
  //   frame_pull->GetXaxis()->SetLabelSize(0.08);
  //   frame_pull->Draw();

  //   // --- chi2/ndf ---
  //   int nFitParam = pdf_fit->getParameters(*data)->selectByAttrib("Constant", kFALSE)->getSize();
  //   double chi2 = frame->chiSquare("model", "data", nFitParam); // (chi2/ndf)
  //   std::cout << "Chi2/ndf = " << chi2 << std::endl;
  //   TLatex latex;
  //   latex.SetNDC();
  //   latex.SetTextSize(0.08);
  //   latex.DrawLatex(0.85, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  //   c1.SaveAs("figs/toy_pr_3expo_ctau.png");
  //   // c1.SaveAs("figs/toy_pr_3expo_ctau.pdf");
  // }

  // // === mass plot ===
  // {
  //   TCanvas c2("c2", "", 800, 800);
  //   c2.Divide(1, 2);
  //   TPad *pad1 = (TPad *)c2.cd(1);
  //   // gPad->SetLogy();
  //   pad1->SetPad(0.0, 0.3, 1.0, 1.0);
  //   pad1->SetBottomMargin(0.02);

  //   TPad *pad2 = (TPad *)c2.cd(2);
  //   pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  //   pad2->SetTopMargin(0.05);
  //   pad2->SetBottomMargin(0.25);

  //   // main plot
  //   pad1->cd();
  //   RooPlot *frame = mass->frame();
  //   data->plotOn(frame, Name("data"));
  //   pdf_fit->plotOn(frame, Name("model"));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeSigpdf"),
  //                   Name("np"),
  //                   LineColor(kRed),
  //                   LineStyle(kDotted),
  //                   LineWidth(3));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeSigpdf2"),
  //                   Name("pr"),
  //                   LineColor(kOrange),
  //                   LineStyle(kDashed),
  //                   LineWidth(2));
  //   pdf_fit->plotOn(frame,
  //                   Components("masstimeBkgpdf"),
  //                   Name("bkg"),
  //                   LineColor(kGreen + 2),
  //                   LineStyle(kDashed),
  //                   LineWidth(2));
  //   frame->Draw();

  //   // --- legend ---
  //   auto leg = new TLegend(0.6, 0.65, 0.95, 0.88);
  //   leg->SetBorderSize(0);
  //   leg->SetFillStyle(0);
  //   leg->AddEntry(frame->findObject("data"), "data", "pe");
  //   leg->AddEntry(frame->findObject("model"), "Total fit", "e");
  //   leg->AddEntry(frame->findObject("pr"), "Prompt J/#psi", "e");
  //   leg->AddEntry(frame->findObject("np"), "NP J/#psi", "e");
  //   leg->AddEntry(frame->findObject("bkg"), "Bkg", "e");
  //   leg->Draw();

  //   // pull
  //   pad2->cd();
  //   RooHist *hpull = frame->pullHist("data", "model"); // data - fit / error
  //   RooPlot *frame_pull = mass->frame();
  //   frame_pull->addPlotable(hpull, "P");
  //   frame_pull->SetTitle("");
  //   frame_pull->GetYaxis()->SetTitle("Pull");
  //   frame_pull->GetYaxis()->SetNdivisions(505);
  //   frame_pull->GetYaxis()->SetTitleSize(0.08);
  //   frame_pull->GetYaxis()->SetLabelSize(0.08);
  //   frame_pull->GetXaxis()->SetTitleSize(0.1);
  //   frame_pull->GetXaxis()->SetLabelSize(0.08);
  //   frame_pull->Draw();

  //   // --- chi2/ndf ---
  //   int nFitParam = pdf_fit->getParameters(*data)->selectByAttrib("Constant", kFALSE)->getSize();
  //   double chi2 = frame->chiSquare("model", "data", nFitParam); // (chi2/ndf)
  //   std::cout << "Chi2/ndf = " << chi2 << std::endl;
  //   TLatex latex;
  //   latex.SetNDC();
  //   latex.SetTextSize(0.08);
  //   latex.DrawLatex(0.85, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  //   c2.SaveAs("figs/toy_pr_3expo_mass.png");
  //   // c2.SaveAs("figs/toy_pr_3expo_mass.pdf");
  // }

  // result->Print("v");

  // cout << "b fraction = N_{NP} / (N_NP + N_PR) = " << N_np->getVal() / (N_np->getVal() + N_pr->getVal()) << "\n";

  cout
      << "=== finish final_2d_fit()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}