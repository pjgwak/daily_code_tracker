#include <gsl/gsl_errno.h>
#include <TStopwatch.h>
#include <RooRealVar.h>
#include "RooFit.h"
#include "RooMinimizer.h"

using namespace RooFit;

// === helper function ===
void fixAndReport(const RooArgList &params, const char *parName)
{
  RooRealVar *par = (RooRealVar *)params.find(parName);
  if (par)
  {
    par->setConstant(kTRUE);
    cout << "Fixed " << parName << " = " << par->getVal() << "\n";
  }
  else
  {
    cerr << "Parameter " << parName << " not found in fitResult" << "\n";
  }
}

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
  double ctMin = -0.5, ctMax = 2;
  RooRealVar ctau3D("ctau3D", "ctau3D", ctMin, ctMax, "mm");

  RooBinning binningCtau(ctMin, ctMax);
  binningCtau.addUniform(nCtauBins, ctMin, ctMax);
  ctau3D.setBinning(binningCtau);

  // RooBinning binningCtau(ctMin, ctMax);
  // binningCtau.addUniform(nCtauBins, ctMin, ctMax);
  // ctau3D.setBinning(binningCtau);

  // === import dataset ===
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form(
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f) && (ctau3D>%.3f && ctau3D < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);
  TString osCut = "(recoQQsign == 0)";
  TString fullCut = Form("%s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data());

  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");

  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!dsRaw)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }
  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight",
                                        dsRaw, *dsRaw->get(), 0, "weight");

  // === new dataset with cuts ===
  RooDataSet *dsReduced = (RooDataSet *)dsWeight->reduce(Cut(fullCut));
  if (!dsReduced || dsReduced->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced->Print();

  // check weight
  const bool hasWeight = (dsReduced->isWeighted() || dsReduced->weightVar() != nullptr);
  // const bool hasWeight = false;
  if (hasWeight)
  {
    // MC's weight value = 1 -> same with no-weight case
    cout << "[Info] Using weighted dataset \n";
  }
  else
  {
    cout << "[Info] Using UN-weighted dataset \n";
  }


  // Build mass model
  // ----------------

  // --- Bring mass fit result ---
  TFile fMassFit(Form("roots/mass_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *resultMass = (RooFitResult *)fMassFit.Get("fitResult");
  if (!resultMass)
  {
    std::cerr << "fitResult not found!" << std::endl;
    return;
  }
  RooArgList paraMass;
  paraMass.add(resultMass->floatParsFinal());
  paraMass.add(resultMass->constPars());

  RooRealVar *NsigMass = (RooRealVar *)paraMass.find("Nsig");
  RooRealVar *NbkgMass = (RooRealVar *)paraMass.find("Nbkg");
  RooRealVar *alphaLMass = (RooRealVar *)paraMass.find("alphaL");
  RooRealVar *alphaRatioMass = (RooRealVar *)paraMass.find("alphaRatio");
  RooRealVar *nLMass = (RooRealVar *)paraMass.find("nL");
  RooRealVar *nRatioMass = (RooRealVar *)paraMass.find("nRatio");
  RooRealVar *sigmaGRatioMass = (RooRealVar *)paraMass.find("sigmaGRatio");
  RooRealVar *sigmaRatioMass = (RooRealVar *)paraMass.find("sigmaRatio");

  RooRealVar *s1Mass = (RooRealVar *)paraMass.find("c1");
  RooRealVar *f_cbRMass = (RooRealVar *)paraMass.find("f_cbR");
  RooRealVar *f_gausMass = (RooRealVar *)paraMass.find("f_gaus");
  RooRealVar *meanMass = (RooRealVar *)paraMass.find("mean");
  RooRealVar *sigmaLMass = (RooRealVar *)paraMass.find("sigmaL");

  fixAndReport(paraMass, "Nsig");
  fixAndReport(paraMass, "Nbkg");
  fixAndReport(paraMass, "alphaL");
  fixAndReport(paraMass, "alphaRatio");
  fixAndReport(paraMass, "nL");
  fixAndReport(paraMass, "nRatio");
  fixAndReport(paraMass, "sigmaGRatio");
  fixAndReport(paraMass, "sigmaRatio");

  fixAndReport(paraMass, "c1");
  fixAndReport(paraMass, "f_cbR");
  fixAndReport(paraMass, "f_gaus");
  fixAndReport(paraMass, "mean");
  fixAndReport(paraMass, "sigmaL");

  // --- Build mass model ---

  // mass sig
  RooFormulaVar sigmaR("sigmaR", "sigma right", "@0*@1", RooArgList(*sigmaLMass, *sigmaRatioMass));
  RooFormulaVar alphaR("alphaR", "alpha right", "-@0*@1", RooArgList(*alphaLMass, *alphaRatioMass));
  RooFormulaVar nR("nR", "n right", "@0*@1", RooArgList(*nLMass, *nRatioMass));
  RooFormulaVar sigmaG("sigmaG", "sigma gauss", "@0*@1", RooArgList(*sigmaLMass, *sigmaGRatioMass));

  RooCrystalBall cbLeft("cbLeft", "CB left tail", mass, *meanMass, *sigmaLMass, *alphaLMass, *nLMass);
  RooCrystalBall cbRight("cbRight", "CB right tail", mass, *meanMass, sigmaR, alphaR, nR);
  RooGaussian gaus("gaus", "gaus core", mass, *meanMass, sigmaG);

  RooAddPdf signal("signal", "", RooArgList(cbLeft, cbRight, gaus), RooArgList(*f_cbRMass, *f_gausMass));

  // mass bkg
  RooChebychev bkg("bkg", "Chebychev", mass, RooArgList(*s1Mass));


  // Build ctau Model
  // ----------------
  
  // --- Bring bkg fit result ---
  TFile fCBkgFit(Form("roots/ctau_bkg_RSB_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *resultCBkg = (RooFitResult *)fCBkgFit.Get("fitResult");
  if (!resultCBkg)
  {
    std::cerr << "fitResult not found!" << std::endl;
    return;
  }
  RooArgList paraCBkg;
  paraCBkg.add(resultCBkg->floatParsFinal());
  paraCBkg.add(resultCBkg->constPars());

  RooRealVar *mu0CBkg = (RooRealVar *)paraCBkg.find("mu0");
  RooRealVar *f12CBkg = (RooRealVar *)paraCBkg.find("f12");
  RooRealVar *f_bkg_leftCBkg = (RooRealVar *)paraCBkg.find("f_bkg_left");
  RooRealVar *f_bkg_midCBkg = (RooRealVar *)paraCBkg.find("f_bkg_mid");
  // RooRealVar *f_bkg_resCBkg = (RooRealVar *)paraCBkg.find("f_bkg_res");
  RooRealVar *frac_g1_in_g12CBkg = (RooRealVar *)paraCBkg.find("frac_g1_in_g12");
  RooRealVar *r21CBkg = (RooRealVar *)paraCBkg.find("r21");
  RooRealVar *r32CBkg = (RooRealVar *)paraCBkg.find("r32");
  RooRealVar *r_leftCBkg = (RooRealVar *)paraCBkg.find("r_left");
  RooRealVar *r_rightCBkg = (RooRealVar *)paraCBkg.find("r_right");
  RooRealVar *sigma1CBkg = (RooRealVar *)paraCBkg.find("sigma1");
  RooRealVar *tau_midCBkg = (RooRealVar *)paraCBkg.find("tau_mid");

  fixAndReport(paraCBkg, "mu0");
  // fixAndReport(paraCBkg, "f12");
  // fixAndReport(paraCBkg, "frac_g1_in_g12");
  // fixAndReport(paraCBkg, "sigma1");
  fixAndReport(paraCBkg, "r21");
  fixAndReport(paraCBkg, "r32");
  fixAndReport(paraCBkg, "r_left");
  fixAndReport(paraCBkg, "r_right");
  fixAndReport(paraCBkg, "tau_mid");
  fixAndReport(paraCBkg, "f_bkg_left");
  // fixAndReport(paraCBkg, "f_bkg_mid");
  // fixAndReport(paraCBkg, "f_bkg_res");

  // --- gauss resolution ---
  RooFormulaVar sigma2_ctau("sigma2_ctau", "@0*@1", RooArgList(*sigma1CBkg, *r21CBkg));
  RooFormulaVar sigma3_ctau("sigma3_ctau", "@0*@1", RooArgList(sigma2_ctau, *r32CBkg));
  
  // fg1 = f12 * frac_g1_in_g12
  // fg2 = f12 * (1-frac_g1_in_g12)
  // fg3 = 1 - f12
  RooFormulaVar fg1_ctau("fg1_ctau", "@0*@1", RooArgList(*f12CBkg, *frac_g1_in_g12CBkg));
  RooFormulaVar fg2_ctau("fg2_ctau", "@0*(1-@1)", RooArgList(*f12CBkg, *frac_g1_in_g12CBkg));
  
  // RooRealVar f1_ctau("f1_ctau", "frac G1 ctau", 0.6, 0.2, 1.0);
  // RooRealVar f2_ctau("f2_ctau", "frac G2 ctau", 0.01, 0.1, 1.0);

  RooGaussModel g1_ctau("g1_ctau", "Gauss1 ctau", ctau3D, *mu0CBkg, *sigma1CBkg);
  RooGaussModel g2_ctau("g2_ctau", "Gauss2 ctau", ctau3D, *mu0CBkg, sigma2_ctau);
  RooGaussModel g3_ctau("g3_ctau", "Gauss3 ctau", ctau3D, *mu0CBkg, sigma3_ctau);

  RooAddModel resModel_ctau("resModel_ctau", "3-Gauss Resolution ctau",
                            RooArgList(g1_ctau, g2_ctau, g3_ctau),
                            RooArgList(fg1_ctau, fg2_ctau));
  
  
  // --- ctau continuum bkg model ---
  RooFormulaVar tau_left("tau_left", "@0*@1", RooArgList(*tau_midCBkg, *r_leftCBkg));
  RooFormulaVar tau_right("tau_right", "@0*@1", RooArgList(*tau_midCBkg, *r_rightCBkg));

  // deacy models
  RooDecay decayBkgL("decayBkgL", "Left decay",
                     ctau3D, tau_left, resModel_ctau, RooDecay::Flipped);
  RooDecay decayBkgMid("decayBkgMid", "Center decay",
                       ctau3D, *tau_midCBkg, resModel_ctau, RooDecay::DoubleSided);
  RooDecay decayBkgR("decayBkgR", "Right decay",
                     ctau3D, tau_right, resModel_ctau, RooDecay::SingleSided);

  // --- total ctau bkg model (PR-like bkg + Continuum bkg) ---
  RooRealVar f_bkg_resCBkg2("f_bkg_resCBkg2", "", 0.8, 0.01, 1);
  RooAddPdf bkgModel_ctau("bkgModel_ctau", "",
                          RooArgList(resModel_ctau, decayBkgMid, decayBkgL, decayBkgR),
                          RooArgList(f_bkg_resCBkg2, *f_bkg_midCBkg, *f_bkg_leftCBkg));

  // --- Bring NP fit result ---
  TFile fCTrueFit(Form("roots/ctau_true_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *resultCTrue = (RooFitResult *)fCTrueFit.Get("fitResult");
  if (!resultCTrue)
  {
    std::cerr << "fitResult not found!" << std::endl;
    return;
  }
  RooArgList paraCTrue;
  paraCTrue.add(resultCTrue->floatParsFinal());
  paraCTrue.add(resultCTrue->constPars());

  RooRealVar *f_right1CTrue = (RooRealVar *)paraCTrue.find("f_right1");
  RooRealVar *tauR1CTrue = (RooRealVar *)paraCTrue.find("tauR1");
  RooRealVar *tauRatioCTrue = (RooRealVar *)paraCTrue.find("tauRatio");

  fixAndReport(paraCTrue, "f_right1");
  // fixAndReport(paraCTrue, "tauR1");
  fixAndReport(paraCTrue, "tauRatio");

  // --- NP model ---
  RooFormulaVar tauR2("tauR2", "@0*@1", RooArgList(*tauR1CTrue, *tauRatioCTrue));
  RooDecay decayTrueR1("decayTrueR1", "",
                       ctau3D, *tauR1CTrue, resModel_ctau, RooDecay::SingleSided);
  RooDecay decayTrueR2("decayTrueR2", "",
                       ctau3D, tauR2, resModel_ctau, RooDecay::SingleSided);

  // combine models
  RooAddPdf decayNP_ctau("decayNP_ctau", "",
                         RooArgList(decayTrueR1, decayTrueR2),
                         RooArgList(*f_right1CTrue), kTRUE);

  // --- ctau signal = Prompt + NP ---
  RooRealVar bFraction("bFraction", "fraction prompt", 0.5, 0.1, 0.8);
  RooFormulaVar fPrSig("fPrSig", "1.0-@0", RooArgList(bFraction));
  RooAddPdf sigModel_ctau("sigModel_ctau", "Signal (Prompt+Nonprompt)",
                          RooArgList(resModel_ctau, decayNP_ctau),
                          RooArgList(fPrSig));



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
  RooRealVar Nsig2D("Nsig2D", "yield signal 2D", 800000, 700000, 1000000);
  RooRealVar Nbkg2D("Nbkg2D", "yield background 2D", 150000, 100000, 200000);

  RooAddPdf model2D("model2D", "Signal + Background (2D)",
                    RooArgList(sig2D, bkg2D),
                    RooArgList(Nsig2D, Nbkg2D));

  // ===== 2D dataset =====
  RooDataHist data2D("data2D", "2D dataset",
                     RooArgList(mass, ctau3D), *dsReduced);

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
                                   Save(), SumW2Error(hasWeight),
                                  Strategy(2),
                                   IntegrateBins(1e-4),
                                   Offset("bins"),
                                   Extended(kTRUE),
                                   NumCPU(32),
                                   EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(1.),
                                   RooFit::PrintEvalErrors(-1),
                                   RooFit::PrintLevel(-1));

  // auto fitResult2D = model2D.fitTo(*dsReduced,
  //                                  Save(), SumW2Error(hasWeight),
  //                                 //  Strategy(2),
  //                                  Offset(),
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


  cout << "=== finish final_2d_fit()";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}