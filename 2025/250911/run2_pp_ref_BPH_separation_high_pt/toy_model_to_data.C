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

void toy_model_to_data()
{
  // ROOT::EnableImplicitMT(24);
  float ptLow = 20, ptHigh = 40;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "cutCtauFull", region = "SR";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_shape_test() ===\n";

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs_region6_test", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!ds)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }

  // === declare cuts ===
  // --- basic cuts ---
  // acceptance
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut

  // kinematics cuts
  //  - correct? (<= cBin <) -> maybe (< cBin <=) ??
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  TString osCut = "(recoQQsign == 0)";

  // --- region6 cuts ---
  const TString cutPR = "(abs(ctau3D) < 0.05)";
  const TString cutNP = "(ctau3D >= 0.10 && ctau3D <= 0.80)";
  const TString cutCtauFull = "(ctau3D >= -0.1 && ctau3D <= 0.5)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";

  // --- combine cuts ---
  TString compCut;
  if (comp == "PR")
    compCut = cutPR;
  else if (comp == "NP")
    compCut = cutNP;
  else
    compCut = cutCtauFull;

  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else
    regionCut = "(1)";

  TString fullCut = Form("%s && %s && %s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data(),
                         compCut.Data(),
                         regionCut.Data());

  // === new dataset with cuts ===
  RooDataSet *data = (RooDataSet *)ds->reduce(Cut(fullCut));
  if (!data || data->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  data->Print();

  // check weight
  // const bool hasWeight = (data->isWeighted() || data->weightVar() || data->get()->find("weight"));
  const bool hasWeight = false;

  // variables
  // === observables ===
  auto *ctau = dynamic_cast<RooRealVar *>(data->get()->find("ctau3D"));
  auto *ctRes = dynamic_cast<RooRealVar *>(data->get()->find("ctau3DRes")); // per-event res
  // ===== 범위: 데이터는 전체, NP는 ct>0만 기여/정규화 =====
  ctau->setMin(-0.1);
  ctau->setMax(0.5);
  ctau->setRange(-0.1, 0.5);
  ctau->setRange("Full", -0.1, 0.5);
  ctau->setRange("Res",-0.05, 0.05);
  // ctau->setRange("pos", 0, 0.05);

  // --- Resolution: 3-Gaussian mixture ---
  // RooRealVar mean("mean", "mean", 9.7162e-04, -0.01, 0.01);
  RooRealVar mean("mean", "mean", 0);
  RooRealVar sigma1("sigma1", "sigma1", 0.3, 0.001, 1);
  // sigma1.setConstant();
  // RooRealVar sigma2("sigma2", "sigma2", 3.5139e-02, 0.01, 1);
  // RooRealVar sigma3("sigma3", "sigma3", 5.7212e-02, 0.01, 0.5);

  // RooRealVar sigma1("sigma1", "sigma1", 0.2, 0.001, 1); // 가장 큼
  // RooRealVar r21("r21", "sigma2/sigma1", 1.1, 1, 3);    // 0<r21<1
  // RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, r21));

  // RooRealVar r32("r32", "sigma3/sigma2", 1.1, 1, 10); // 0<r32<1
  // RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, r32));

  //     Floating Parameter    FinalValue +/-  Error
  // --------------------  --------------------------
                //    fG1    7.0639e-01 +/-  3.07e-02
                //    fG2    2.7884e-01 +/-  2.87e-02
                //    mu0    9.7162e-04 +/-  2.43e-03
                //    r21    1.8848e+00 +/-  3.83e-02
                //    r32    2.6312e+00 +/-  9.12e-02
                // sigma1    7.2705e-01 +/-  1.26e-02

  RooRealVar r2("r2", "ratio sigma2/sigma1", 1.4413e+00, 1.0, 10.0);
  RooRealVar r3("r3", "ratio sigma3/sigma1", 2.2415e+00, 1.0, 10.0);
  r2.setConstant();
  r3.setConstant();

  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, r2));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma1, r3));

  RooGaussModel g1("g1", "res1", *ctau, mean, sigma1);
  RooGaussModel g2("g2", "res2", *ctau, mean, sigma2);
  RooGaussModel g3("g3", "res3", *ctau, mean, sigma3);

  RooRealVar f1("f1", "frac1", 5.3087e-01, 0, 1);
  RooRealVar f2("f2", "frac2", 9.2618e-01, 0, 1);
  // f1.setConstant(true);
  // f2.setConstant(true);

  RooAddModel resModel("resModel", "3-Gaussian resolution",
                       RooArgList(g1, g2, g3),
                       RooArgList(f1, f2));

  // --- Prompt = delta ⊗ resolution ---
  RooAbsPdf *prompt = &resModel;

  // --- Nonprompt lifetime ---
  RooRealVar tau("tau", "nonprompt lifetime", 0.45, 0.3, 1);
  RooDecay decay("decay", "nonprompt", *ctau, tau, resModel, RooDecay::SingleSided);

  // RooRealVar tau1("tau1", "short lifetime", 0.3, 0.1, 1.0);
  // RooRealVar tau1Ratio("tau1Ratio", "", 1.14, 0.5, 2);
  // tau1Ratio.setConstant();
  // RooFormulaVar tau2("tau2", "@0 * @1", RooArgList(tau1, tau1Ratio));
  // // RooRealVar tau2("tau2", "long lifetime", 0.5, 0.1, 1.0)
  
  // RooRealVar fNP("fNP", "frac short", 0.571, 0., 1.);
  // fNP.setConstant();

  // RooDecay decay1("decay1", "NP short", *ctau, tau1, resModel, RooDecay::SingleSided);
  // RooDecay decay2("decay2", "NP long", *ctau, tau2, resModel, RooDecay::SingleSided);
  // RooAddPdf decay("decay", "NP sum",
  //                 RooArgList(decay1, decay2), RooArgList(fNP));

  // --- Background lifetime ---
  // RooRealVar tau_bkg("tau_bkg", "bkg lifetime", 4.4669e-01, 0.1, 1);
  // tau_bkg.setConstant(true);
  // RooDecay bkg("bkg", "background", *ctau, tau_bkg, resModel, RooDecay::SingleSided);

  RooRealVar tau_bkg1("tau_bkg1", "short lifetime", 4.2119e-01, 0.1, 1.0);
  // RooRealVar tau_bkg1Ratio("tau_bkg1Ratio", "", 1.14, 0.5, 2);
  // tau_bkg1Ratio.setConstant();
  // RooFormulaVar tau_bkg2("tau_bkg2", "@0 * @1", RooArgList(tau_bkg1, tau_bkg1Ratio));
  RooRealVar tau_bkg2("tau_bkg2", "long lifetime", 4.2422e-01, 0.01, 1.0);

  tau_bkg1.setConstant();
  tau_bkg2.setConstant();

  RooRealVar fNP_bkg("fNP_bkg", "frac short", 4.6362e-03, 0., 1.);
  fNP_bkg.setConstant();

  RooDecay decay_bkg1("decay_bkg1", "NP short", *ctau, tau_bkg1, resModel, RooDecay::SingleSided);
  RooDecay decay_bkg2("decay_bkg2", "NP long", *ctau, tau_bkg2, resModel, RooDecay::SingleSided);
  RooAddPdf bkg("bkg", "NP sum",
                  RooArgList(decay_bkg1, decay_bkg2), RooArgList(fNP_bkg));

  // --- Yields ---
  RooRealVar Npr("Npr", "yield prompt", 1000, 100, 1000000);
  RooRealVar Nnp("Nnp", "yield nonprompt", 500, 100, 50000);
  RooRealVar NbkgRatio("NbkgRatio", "yield nonprompt", 0.26, 0.01, 0.4);
  NbkgRatio.setConstant();
  RooFormulaVar Nbkg("Nbkg", "@0 * @1", RooArgList(Nnp, NbkgRatio));
  // RooRealVar Nbkg("Nbkg", "yield bkg", 30000, 0, 1e7);

  RooExtendPdf extPR("extPR", "extended prompt", *prompt, Npr);
  RooExtendPdf extNP("extNP", "extended nonprompt", decay, Nnp);
  RooExtendPdf extBkg("extBkg", "extended background", bkg, Nbkg);

  RooAddPdf model("model", "PR+NP+Bkg",
                  RooArgList(extPR, extNP, extBkg));

  // --- Generate toy dataset ---
  // auto data = model.generate(t, 50000);

  // --- Fit ---
  auto dh = new RooDataHist("dh", "binned dataset", *ctau, *data);

  auto result = model.fitTo(*dh, Extended(kTRUE), PrintLevel(-1), Save(), Offset(true), NumCPU(32), EvalBackend("legacy"), Strategy(2)); //, BatchMode(true)

  // --- Plot main frame ---
  RooPlot *frame = ctau->frame();
  dh->plotOn(frame, Name("data"));
  model.plotOn(frame, Name("curveAll"), Range("Full"), NormRange("Full")); // 전체 모델에 이름 지정
  model.plotOn(frame, Components("extPR"), LineStyle(kDashed), LineColor(kRed));
  model.plotOn(frame, Components("g1"), LineStyle(kDotted), LineColor(kSpring));
  model.plotOn(frame, Components("g2"), LineStyle(kDotted), LineColor(kAzure));
  model.plotOn(frame, Components("g3"), LineStyle(kDotted), LineColor(kBlack));
  model.plotOn(frame, Components("extNP"), LineStyle(kDashed), LineColor(kMagenta));
  model.plotOn(frame, Components("extBkg"), LineStyle(kDashed), LineColor(kViolet-9));
  model.plotOn(frame, Components("extBkg,extNP"), LineStyle(kDashed), LineColor(kViolet), Name("NPContribution"));

  // --- Pull distribution (curveAll과 비교) ---
  RooHist *hpull = frame->pullHist("data", "curveAll");
  RooPlot *frame_pull = ctau->frame(Title("Pull Distribution"));
  frame_pull->addPlotable(hpull, "P");
  frame_pull->SetMinimum(-8);
  frame_pull->SetMaximum(8);

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
  int nFitParam = result->floatParsFinal().getSize();
  double chi2ndf = frame->chiSquare("curveAll", "data", nFitParam);
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  // lower pad: pull
  TPad *pad2 = (TPad *)c->cd(2);
  pad2->SetPad(0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  frame_pull->Draw();

  // --- Save result ---
  c->SaveAs("fit_result_with_pull.png");

  result->Print("V");

  cout << "Chi2/ndf: " << chi2ndf << "\n";

  ctau->setRange("sig", -0.05, 0.05);

  // 2) 관측 변수 집합
  RooArgSet obs(*ctau);

  // 3) 컴포넌트별 구간 확률 (NormSet 꼭!)
  auto fPrompt = prompt->createIntegral(obs, NormSet(obs), Range("sig"))->getVal();
  auto fNonprompt = decay.createIntegral(obs, NormSet(obs), Range("sig"))->getVal();
  auto fBkg = bkg.createIntegral(obs, NormSet(obs), Range("sig"))->getVal();

  // 4) 구간 기대 이벤트 수 = (구간 확률) × (전체 yield 값)
  double Npr_in = Npr.getVal() * fPrompt;
  double Nnp_in = Nnp.getVal() * fNonprompt;
  double Nbkg_in = Nbkg.getVal() * fBkg;

  std::cout << "[-0.05,0.05] yields:\n"
            << "  Prompt     = " << Npr_in << "\n"
            << "  Nonprompt  = " << Nnp_in << "\n"
            << "  Background = " << Nbkg_in << "\n"
            << "  Total      = " << (Npr_in + Nnp_in + Nbkg_in) << std::endl;

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
