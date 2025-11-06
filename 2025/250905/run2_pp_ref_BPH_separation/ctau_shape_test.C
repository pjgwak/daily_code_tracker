#include <iostream> // std::cout
#include <cstdio>   // printf
#include <string>
#include <TStopwatch.h>
#include <TSystem.h> // gSystem
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <RooGlobalFunc.h> // using namespace RooFit;
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>     // RooFormulaVar parameter list
#include <RooExponential.h>
#include <RooFitResult.h>

using namespace RooFit;
using std::cout;
using std::string;

// --- next macro: ctau_shape_test() ---
// apply cuts
// basic cuts
// region6

void ctau_shape_test()
{
  float ptLow = 6.5, ptHigh = 9;
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
  RooDataSet *ds_red = (RooDataSet *)ds->reduce(Cut(fullCut));
  if (!ds_red || ds_red->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  ds_red->Print();

  // check weight
  // const bool hasWeight = (ds_red->isWeighted() || ds_red->weightVar() || ds_red->get()->find("weight"));
  const bool hasWeight = false;

  // variables
  // === observables ===
  auto *ctau = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3D"));
  auto *ctRes = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3DRes")); // per-event res
  // ===== 범위: 데이터는 전체, NP는 ct>0만 기여/정규화 =====
  ctau->setRange(-0.1, 0.5);
  ctau->setRange("fullRange", -0.1, 0.5);
  ctau->setRange("posRange", 0.0, ctau->getMax());

  // ===== 공통: 분해능 평균 =====
  RooRealVar muRes("muRes", "resolution mean", 0.0, -0.02, 0.02);

  // ===== 소프트-플로어 σ (항상 > 0): sigma = sqrt(floor^2 + a^2) =====
  RooRealVar sFloor("sFloor", "sigma floor", 0.005); // 고정 바닥
  sFloor.setConstant(kTRUE);

  // --- NP용 폭(컨볼루션용 RooGaussModel) ---
  RooRealVar a1("a1", "raw width 1", 0.010, 0.001, 0.30);
  RooRealVar a2("a2", "raw width 2", 0.020, 0.001, 0.60);
  RooRealVar a3("a3", "raw width 3", 0.050, 0.001, 1.00);

  RooFormulaVar sig1("sig1", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a1));
  RooFormulaVar sig2("sig2", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a2));
  RooFormulaVar sig3("sig3", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a3));

  // --- NP 분해능: 3-Gauss (항상 양수, 합=1) ---
  RooGaussModel G1("G1", "G1", *ctau, muRes, sig1);
  RooGaussModel G2("G2", "G2", *ctau, muRes, sig2);
  RooGaussModel G3("G3", "G3", *ctau, muRes, sig3);

  RooRealVar f1("f1", "core fraction", 0.70, 0.0, 1.0);
  RooRealVar f2p("f2p", "tail2 in remainder", 0.50, 0.0, 1.0);
  RooFormulaVar f2("f2", "(1.0-@0)*@1", RooArgList(f1, f2p));
  // RooAddModel의 마지막 항은 자동으로 1-f1-f2 → 음수 불가
  RooAddModel Res3G("Res3G", "triple Gaussian res", RooArgList(G1, G2, G3), RooArgList(f1, f2));

  // --- Nonprompt: [exp(-ct/t)×H(ct)] ⊗ Res3G  (ct<0에서 정확히 0) ---
  RooRealVar tauFloor("tauFloor", "lifetime soft floor", 0.05); // 단위: ct, 고정
  tauFloor.setConstant(kTRUE);

  // 원하는 초기 τ≈0.20에 맞춰 aTau 초기값 지정:
  // aTau0 = sqrt(0.20^2 - 0.05^2) ≈ 0.193649
  RooRealVar aTau("aTau", "lifetime free term", 0.193649, 0.0, 10.0);

  RooFormulaVar tau("tau", "sqrt(@0*@0 + @1*@1)", RooArgList(tauFloor, aTau));

  RooDecay np_raw("np_raw", "NP raw", *ctau, tau, Res3G, RooDecay::SingleSided);

  // ===== Prompt: 3-Gauss (pdf 자체, 항상 ≥0, 합=1) =====
  RooRealVar a1p("a1p", "raw width 1P", 0.006, 0.001, 0.30);
  RooRealVar a2p("a2p", "raw width 2P", 0.015, 0.001, 0.60);
  RooRealVar a3p("a3p", "raw width 3P", 0.040, 0.001, 1.00);

  RooFormulaVar sp1("sp1", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a1p));
  RooFormulaVar sp2("sp2", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a2p));
  RooFormulaVar sp3("sp3", "sqrt(@0*@0 + @1*@1)", RooArgList(sFloor, a3p));

  RooRealVar muP("muP", "prompt mean", 0.0, -0.02, 0.02);
  RooGaussian P1("P1", "P1", *ctau, muP, sp1);
  RooGaussian P2("P2", "P2", *ctau, muP, sp2);
  RooGaussian P3("P3", "P3", *ctau, muP, sp3);

  RooRealVar g1("g1", "prompt core frac", 0.80, 0.0, 1.0);
  RooRealVar g2p("g2p", "prompt tail2 in rem", 0.50, 0.0, 1.0);
  RooFormulaVar g2("g2", "(1.0-@0)*@1", RooArgList(g1, g2p));
  RooAddPdf prPdf("prPdf", "prompt 3G", RooArgList(P1, P2, P3), RooArgList(g1, g2));

  // ===== 확장 파트: 성분별 범위 해석(중요) =====
  double nData = ds_red->sumEntries();

  RooRealVar Npr("Npr", "yield prompt (full)", 0.85 * nData, 0, 10. * nData);
  RooRealVar Nnp("Nnp", "yield NP (ct>0 only)", 0.15 * nData, 0, 10. * nData);

  // Prompt는 fullRange 해석(기본). NP는 posRange에서의 이벤트 수로 해석.
  RooExtendPdf pr_ext("pr_ext", "prompt ext", prPdf, Npr /* no rangeName */);
  RooExtendPdf np_ext("np_ext", "NP ext (posRange)", np_raw, Nnp, "posRange");

  // ===== 최종 모델 (항상 ≥0) =====
  RooAddPdf model("model", "pr+np", RooArgList(pr_ext, np_ext));

  // ===== combined model =====
  // RooAddPdf model("model", "prompt+nonprompt",
  //                 RooArgList(pr, np), RooArgList(Npr, Nnp));

  // === full fit (여기서만 Extended) ===
  auto r = model.fitTo(*ds_red, Save(), Range("fitAll"), SumW2Error(true), Extended(kTRUE), PrintLevel(-1));

  // === plot: Range/NormRange + ConditionalObservables를 명시 ===
  TCanvas c1("c1", "c1", 800, 800);
  auto fr = ctau->frame(Range("fitAll"));
  ds_red->plotOn(fr, DataError(RooAbsData::SumW2));
  model.plotOn(fr, Range("fitAll"), NormRange("fitAll"));
  fr->Draw();
  c1.SaveAs("shape_test.png");

  r->Print("V");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}