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

// --- next macro: ctau_main_fit() ---
// apply cuts
// basic cuts
// region6

void ctau_main_fit()
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "cutCtauFull", region = "SR";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_main_fit() ===\n";
  
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
  const TString cutCtauFull = "(ctau3D >= -0.1 && ctau3D <= 0.1)";

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
  RooRealVar *ctau3DVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3D"));
  if (!ctau3DVar)
    cout << "Warn: There is no variable 'ctau3D'\n";
  ctau3DVar->setRange(-0.1, 0.1); // Change PDF's range - Need.
  ctau3DVar->setRange("fitRange", -0.1, 0.1);

  // === bring parameters from previous steps ===
  // --- ctau MC Res ---
  TFile fMcRes("fitResult_mc_ctau.root");
  RooFitResult *fr_McRes = (RooFitResult *)fMcRes.Get("fitResult");
  if (!fr_McRes)
  {
    std::cerr << "fitResult MC Res not found!" << std::endl;
    return;
  }

  // // 최종 파라미터 리스트
  const RooArgList &params_McRes = fr_McRes->floatParsFinal();

  // // 예시: mean, sigmaL, sigmaR 만 가져와서 고정
  RooRealVar *sigma1Var = (RooRealVar *)params_McRes.find("sigma1");
  RooRealVar *mu0Var = (RooRealVar *)params_McRes.find("mu0");
  RooRealVar *r21Var = (RooRealVar *)params_McRes.find("r21");
  RooRealVar *r32Var = (RooRealVar *)params_McRes.find("r32");
  RooRealVar *fG1Var = (RooRealVar *)params_McRes.find("fG1");
  RooRealVar *fG2Var = (RooRealVar *)params_McRes.find("fG2");

  // if (sigma1Var)
  // {
  //   sigma1Var->setConstant(kTRUE);
  //   std::cout << "Fixed sigma1 = " << sigma1Var->getVal() << std::endl;
  // }
  if (mu0Var)
  {
    mu0Var->setConstant(kTRUE);
    std::cout << "Fixed m0 = " << mu0Var->getVal() << std::endl;
  }
  if (r21Var)
  {
    r21Var->setConstant(kTRUE);
    std::cout << "Fixed r21 = " << r21Var->getVal() << std::endl;
  }
  if (r32Var)
  {
    r32Var->setConstant(kTRUE);
    std::cout << "Fixed r32 = " << r32Var->getVal() << std::endl;
  }
  if (fG1Var)
  {
    fG1Var->setConstant(kTRUE);
    std::cout << "Fixed fG1 = " << fG1Var->getVal() << std::endl;
  }
  if (fG2Var)
  {
    fG2Var->setConstant(kTRUE);
    std::cout << "Fixed fG2 = " << fG2Var->getVal() << std::endl;
  }

  // --- ctau Data Side SR ---
  TFile fSideSR("fitResult_ctau_side_SR.root");
  RooFitResult *fr_SideSR = (RooFitResult *)fSideSR.Get("fitResult");
  if (!fr_SideSR)
  {
    std::cerr << "fitResult Side SR not found!" << std::endl;
    return;
  }

  const RooArgList &params_SideSR = fr_SideSR->floatParsFinal();
  RooRealVar *N1Var = (RooRealVar *)params_SideSR.find("N");

  // --- ctau Data Side LSB ---
  TFile fSideLSB("fitResult_ctau_side_LSB.root");
  RooFitResult *fr_SideLSB = (RooFitResult *)fSideLSB.Get("fitResult");
  if (!fr_SideLSB)
  {
    std::cerr << "fitResult Side LSB not found!" << std::endl;
    return;
  }

  const RooArgList &params_SideLSB = fr_SideLSB->floatParsFinal();
  RooRealVar *N2Var = (RooRealVar *)params_SideLSB.find("N");
  RooRealVar *tau1Var = (RooRealVar *)params_SideLSB.find("tau1");
  RooRealVar *tau2Var = (RooRealVar *)params_SideLSB.find("tau2");
  RooRealVar *f1Var = (RooRealVar *)params_SideLSB.find("f1");

  RooFormulaVar c1("c1", "-1.0/@0", RooArgList(*tau1Var));
  RooFormulaVar c2("c2", "-1.0/@0", RooArgList(*tau2Var));
  RooExponential e1("e1", "e1", *ctau3DVar, c1);
  RooExponential e2("e2", "e2", *ctau3DVar, c2);
  RooAddPdf expo2("expo2", "f1*e1 + (1-f1)*e2", RooArgList(e1, e2), RooArgList(*f1Var));

  // if (XXX)
  // {
  //   XXX->setConstant(kTRUE);
  //   std::cout << "Fixed alphaL = " << XXX->getVal() << std::endl;
  // }

  RooRealVar sigma1("sigma1", "sigma1", 0.001, 0.00001, 0.01); // 가장 큼
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(*sigma1Var, *r21Var));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, *r32Var));
  RooGaussian G1("G1", "G1", *ctau3DVar, *mu0Var, *sigma1Var);
  RooGaussian G2("G2", "G2", *ctau3DVar, *mu0Var, sigma2);
  RooGaussian G3("G3", "G3", *ctau3DVar, *mu0Var, sigma3);
  RooAddPdf Res3("Res3", "3-Gauss resolution",
                 RooArgList(G1, G2, G3), RooArgList(*fG1Var, *fG2Var));

  // RooDecay에 넣기 위한 "resolution model" (RooResolutionModel 계열)
  RooGaussModel GM1("GM1", "GM1", *ctau3DVar, *mu0Var, *sigma1Var);
  RooGaussModel GM2("GM2", "GM2", *ctau3DVar, *mu0Var, sigma2);
  RooGaussModel GM3("GM3", "GM3", *ctau3DVar, *mu0Var, sigma3);
  RooAddModel ResModel("ResModel", "3-Gauss resolution model",
                       RooArgList(GM1, GM2, GM3), RooArgList(*fG1Var, *fG2Var));

  // 3) PR = Res
  //  - Prompt 성분은 delta(0) ⊗ Resolution = Resolution 자체
  RooAbsPdf &PR = Res3; // 이름별도로 원하시면 RooAddPdf PR("PR","PR", ... )로 별칭 가능

  // 4) Background = NP + NP continuum
  //  - 둘 다 ctau>0 단측 지수 ⊗ Resolution
  RooRealVar tau_np("tau_np", "tau_np [mm]", 0.20, 0.00001, 10.0);
  RooRealVar tau_cont("tau_cont", "tau_cont [mm]", 1.00, 0.00000, 20.0);

  // decreasing exponential (ctau>0): RooDecay::SingleSided
  RooDecay NP_core("NP_core", "NP core (exp ⊗ Res)",
                   *ctau3DVar, tau_np, ResModel, RooDecay::SingleSided);
  RooDecay NPC_core("NPC_core", "NP continuum (exp ⊗ Res)",
                    *ctau3DVar, tau_cont, ResModel, RooDecay::SingleSided);

  // Bkg 내부 비율(이전 결과에서 고정한다고 하셨으니 setConstant(true) 예시)
  RooRealVar fNP_in_Bkg("fNP_in_Bkg", "frac of NP in Bkg", 0.2, 0.0, 1.0);
  fNP_in_Bkg.setVal(N2Var->getVal() / N1Var->getVal());
  fNP_in_Bkg.setConstant(true); // 이전 결과 고정

  RooAddPdf Bkg("Bkg", "Bkg = NP + NP_cont",
                RooArgList(NPC_core, NP_core), RooArgList(fNP_in_Bkg));

  // (2) Extended 방식 (권장: 수치적으로 안정적)
  //   N_PR, N_BKG 각각의 절대 yield 사용
  RooRealVar N_PR("N_PR", "prompt yield", 6.4075e+04, 1, 1e6);
  RooRealVar N_BK("N_BK", "background yield", 6.9891e+07, 5e+07, 1e8);
  RooAddPdf model("model", "N_PR*PR + N_BK*Bkg",
                      RooArgList(PR, Bkg), RooArgList(N_PR, N_BK));

  // --- perform fit ---
  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(1).removeTopic(Tracing);
  // RooMsgService::instance().getStream(0).removeTopic(Caching);
  // RooMsgService::instance().getStream(1).removeTopic(Caching);
  // RooMsgService::instance().getStream(0).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Integration);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  // RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"), Extended(), SumW2Error(hasWeight), Offset(true), PrintLevel(-1), NumCPU(30), Warnings(kFALSE), Verbose(kFALSE), RecoverFromUndefinedRegions(1.5), PrintEvalErrors(-1));

  // === draw ctau3D ===

  double chi2ndf = 0;

  if (ctau3DVar)
  {
    // --- divided canvas ---
    TCanvas c_ctau("c_ctau", "c_ctau", 800, 800);

    // --- main plot ---
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1->SetBottomMargin(0.00001);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    double ctMin = -0.1, ctMax = 0.1;
    RooPlot *f_ctau = ctau3DVar->frame(Range(ctMin, ctMax), Title("")); // Bins(80)
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
    c_ctau.cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();

    RooHist *hpull = f_ctau->pullHist();
    RooPlot *f_pull = ctau3DVar->frame(Range(ctMin, ctMax), Title(""));
    f_pull->addPlotable(hpull, "P"); // P: points only

    f_pull->GetYaxis()->SetTitle("Pull");
    f_pull->GetXaxis()->SetTitle("c#tau_{3D} [mm]");
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
    double xmin = ctMin;
    double xmax = ctMax;
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

    c_ctau.SaveAs("figs_region6_test/ctau3D_cent.png");
  }

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}