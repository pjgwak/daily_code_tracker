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
#include <RooArgList.h> // RooFormulaVar parameter list
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooCrystalBall.h>

using namespace RooFit;
using std::cout;
using std::string;

// --- next macro: mass_center_fit() ---
// apply cuts
// basic cuts
// region6

void mass_center_fit()
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "NP", region = "MassFull";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start mass_center_fit() ===\n";

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs_region6_test", true);

  // read input
  // TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root");

  // OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root
  // OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root
  // OniaRooDataSet_JPsi_pp_GENONLY_NonPrompt_230215.root

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
  const TString cutCtauFull = "(ctau3D >= -.10 && ctau3D <= 0.80)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";
  const TString cutMassFull = "(mass >= 2.6 && mass <= 3.5)";

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
  else if (region == "MassFull")
    regionCut = cutMassFull;
  else regionCut = "(1)";

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
  const bool hasWeight = (ds_red->isWeighted() || ds_red->weightVar() || ds_red->get()->find("weight"));
  // const bool hasWeight = false;

  // variables
  RooRealVar *massVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("mass"));
  if (!massVar)
    cout << "Warn: There is no variable 'mass'\n";
  massVar->setRange(2.6, 3.5); // Change PDF's range - Need.
  massVar->setRange("fitRange", 2.6, 3.5);

  // === fitting ===
  double mean0 = 3.0969;
  double sigma1_0 = 0.020;
  double sigma2_0 = 0.040;
  double sigmaG_0 = 0.015;

  // ===== 공통: 시그널 평균(공유) =====
  RooRealVar mean("mean", "signal mean", mean0, mean0 - 0.050, mean0 + 0.050);

  // ----- 왼쪽 테일 Crystal Ball -----
  RooRealVar sigmaL("sigmaL", "sigma left", 0.025, 0.001, 0.100);

  // RooRealVar alphaL("alphaL", "alpha left", 1.5, 0.5, 2.0); // >0 → 왼쪽 테일
  // RooRealVar nL("nL", "n left", 3.0, 1, 100.0);

  RooRealVar aLprime("aLprime", "log(alphaL)", std::log(1.5), std::log(0.5), std::log(2.0));
  RooFormulaVar alphaL("alphaL", "exp(@0)", RooArgList(aLprime));
  RooRealVar nLprime("nLprime", "log(nL-1)", std::log(3.0 - 1.0), std::log(0.3), std::log(20.0));
  RooFormulaVar nL("nL", "1+exp(@0)", RooArgList(nLprime));

  RooCBShape cbLeft("cbLeft", "CB left tail", *massVar, mean, sigmaL, alphaL, nL);

  // ----- 오른쪽 테일 Crystal Ball -----
  RooRealVar sigmaR("sigmaR", "sigma right", 0.025, 0.01, 0.05);
  // RooRealVar alphaR("alphaR", "alpha right", -1.5, -5.0, -0.01);
  RooRealVar alphaRatio("alphaRatio", "alphaR/alphaL ratio", 1.0, 0.5, 5.0);
  RooFormulaVar alphaR("alphaR", "alpha right", "-@0*@1", RooArgList(alphaL, alphaRatio));
  RooRealVar nRatio("nRatio", "nR/nL ratio", 1.0, 0.5, 100.0);
  RooFormulaVar nR("nR", "n right", "@0*@1", RooArgList(nL, nRatio));
  RooCBShape cbRight("cbRight", "CB right tail", *massVar, mean, sigmaR, alphaR, nR);

  // ----- 추가 Gaussian -----
  RooRealVar sigmaG("sigmaG", "sigma gaussian", 0.015, 0.001, 0.080);
  RooGaussian gaus("gaus", "gaus core", *massVar, mean, sigmaG);

  // ----- 혼합 -----
  RooRealVar f_cbR("f_cbR", "frac of right CB", 0.5, 0, 1.0);
  RooRealVar f_gaus("f_gaus", "frac of gaus", 0.05, 0, 1);
  RooAddPdf sig("sig", "signal = CBleft ⊕ CBright ⊕ Gauss",
                RooArgList(cbLeft, cbRight, gaus),
                RooArgList(f_cbR, f_gaus));
  // RooAddPdf sig("sig", "signal = CBleft ⊕ CBright ⊕ Gauss",
  //             RooArgList(cbLeft, cbRight),
  //               RooArgList(f_cbR));

  // ===== 배경: 2차 Chebyshev =====
  // RooChebychev의 차수 2 → 계수 3개 필요 (c0, c1, c2)
  RooRealVar c1("c1", "c1", 0.0, -1, 1);
  RooRealVar c2("c2", "c2", 0.0, -1.0, 1.0);
  RooRealVar c3("c3", "c3", 0.0, -1.0, 1.0);
  RooChebychev bkg("pdf_mass_bkg", "2nd-order Chebychev", *massVar, RooArgList(c1, c2, c3));

  // ===== Extended yields =====
  RooRealVar Nsig("Nsig", "signal yield", 6.7089e+05, 100000, 800000);
  RooRealVar Nbkg("Nbkg", "background yield", 5.9070e+04, 10000, 90000);

  RooAddPdf model("pdf_mass_tot", "signal ⊕ background (extended)",
                RooArgList(sig, bkg), RooArgList(Nsig, Nbkg));

  // --- perform fit ---
  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"), SumW2Error(hasWeight), Offset(true), Extended(kTRUE), PrintLevel(-1), Warnings(kFALSE), Verbose(kFALSE), NumCPU(32), EvalBackend("legacy"), Strategy(2));

  // RooFitResult *fitResult = nullptr;

  // // Strategy 0 → 1 → 2 순차 시도
  // for (int strat = 0; strat <= 2; ++strat)
  // {
  //   fitResult = model.fitTo(
  //       *ds_red,
  //       Save(),
  //       Range("fitRange"),
  //       SumW2Error(hasWeight),
  //       Offset(true),
  //       Extended(kTRUE),
  //       PrintLevel(-1),
  //       NumCPU(30),
  //       Warnings(kFALSE),
  //       Verbose(kFALSE),
  //       Strategy(strat));

  //   // 성공 판정: status==0 && edm < 기준값
  //   if (fitResult && fitResult->status() == 0 && fitResult->edm() < 1e-4)
  //   {
  //     std::cout << "Converged with strategy " << strat << std::endl;
  //     break;
  //   }
  // }

  // RooFitResult *fitResult = nullptr;
  // struct Step
  // {
  //   const char *algo;
  //   int strat;
  //   bool hesse;
  // };
  // std::vector<Step> steps = {
  //     {"Simplex", 0, false}, // 자리잡기
  //     {"Migrad", 0, true},   // 기본
  //     {"Migrad", 1, true},   // 강화
  //     {"Migrad", 2, true}    // 최강
  // };

  // const double EDM_TARGET = 1e-4;

  // for (const auto &s : steps)
  // {
  //   fitResult = model.fitTo(
  //       *ds_red,
  //       Save(),
  //       Range("fitRange"),
  //       SumW2Error(hasWeight),
  //       Offset(true),
  //       Extended(kTRUE),
  //       Minimizer("Minuit2", s.algo), // ← SIMPLEX / MIGRAD 전환
  //       Strategy(s.strat),
  //       Hesse(s.hesse), // MIGRAD 단계에서 HESSE 함께
  //       PrintLevel(-1),
  //       NumCPU(30),
  //       Warnings(kFALSE),
  //       Verbose(kFALSE));

  //   // 수렴 판정
  //   if (fitResult &&
  //       fitResult->status() == 0 &&
  //       fitResult->covQual() >= 2 && // 2~3 권장
  //       fitResult->edm() < EDM_TARGET)
  //   {
  //     std::cout << "[OK] algo=" << s.algo
  //               << " strat=" << s.strat
  //               << " edm=" << fitResult->edm()
  //               << " covQual=" << fitResult->covQual()
  //               << std::endl;
  //     break;
  //   }
  //   else
  //   {
  //     std::cout << "[Retry] algo=" << s.algo
  //               << " strat=" << s.strat
  //               << " edm=" << (fitResult ? fitResult->edm() : -1)
  //               << " covQual=" << (fitResult ? fitResult->covQual() : -1)
  //               << " status=" << (fitResult ? fitResult->status() : -1)
  //               << std::endl;
  //   }
  // }

  // === draw ctau3D ===
  double chi2ndf = 0;

  if (massVar)
  {
    // --- divided canvas ---
    TCanvas c_mass("c_mass", "c_mass", 800, 800);

    // --- main plot ---
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1->SetBottomMargin(0.00001);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    double massMin = 2.6, massMax = 3.5;
    RooPlot *f_ctau = massVar->frame(Range(massMin, massMax), Title("")); // Bins(80)
    if (hasWeight)
      ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2), WeightVar("weight"), Name("ds_red"));
    else
      ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("ds_red"));
    model.plotOn((f_ctau), NormRange("fitRange"), Range("fitRange"), Name("model"));

    // y axis: logY style
    double ymin = 1e300, ymax = -1e300;
    RooHist *hdata = (RooHist *)f_ctau->getHist("ds_red"); // use first dataset on f_ctau
    if (hdata)
    {
      for (int i = 0; i < hdata->GetN(); i++)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }

    double floor = 1e-3;
    if (ymin <= 0 || ymin == 1e300)
      ymin = floor;

    f_ctau->SetMinimum(ymin * 0.5);
    f_ctau->SetMaximum(ymax * 10.0);

    // title
    f_ctau->GetYaxis()->SetTitle("Events");
    f_ctau->GetXaxis()->SetTitle("");
    f_ctau->Draw("e");

    // --- pull pad ---
    c_mass.cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();

    RooHist *hpull = f_ctau->pullHist();
    RooPlot *f_pull = massVar->frame(Range(massMin, massMax), Title(""));
    f_pull->addPlotable(hpull, "P"); // P: points only

    f_pull->GetYaxis()->SetTitle("Pull");
    f_pull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^2]");
    f_pull->GetXaxis()->CenterTitle();
    f_pull->SetMinimum(-8);
    f_pull->SetMaximum(8);
    f_pull->GetYaxis()->SetNdivisions(505);
    f_pull->GetYaxis()->SetTitleSize(0.12);
    f_pull->GetYaxis()->SetLabelSize(0.10);
    f_pull->GetXaxis()->SetTitleSize(0.15);
    f_pull->GetXaxis()->SetLabelSize(0.10);
    f_pull->Draw();

    // --- draw pull = 0 line ---
    double xmin = massMin;
    double xmax = massMax;
    TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
    // line->SetLineColor();
    line->SetLineStyle(2);
    line->Draw("same");

    // --- compute and draw chi square ---
    int nFitParam = fitResult->floatParsFinal().getSize();
    chi2ndf = f_ctau->chiSquare("model", "ds_red", nFitParam);

    TLatex latex;
    latex.SetNDC(); // use pad coordinates (0~1)
    latex.SetTextSize(0.1);
    latex.DrawLatex(0.82, 0.88, Form("#chi^{2}/ndf = %.2f", chi2ndf));

    c_mass.SaveAs(Form("figs_region6_test/mass_center_%s.png", comp.c_str()));
  }

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}