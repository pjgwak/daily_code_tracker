#include "utils/mass_fit_utils.h"
#include "cfgs/pp_run2_pt40p50_y0_1p6.C"
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFormulaVar.h>
#include <RooFitResult.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>

using namespace RooFit;
using std::cout;

void defineMassSig(RooWorkspace *ws, FitConfigMcMass &cfg);
void draw_mass_fit(RooWorkspace *ws, RooRealVar &mass, RooAbsPdf &model, RooAbsData &data, const FitConfigMcMass &cfg, const RooFitResult *fitResult, int nBins = 80, const char *outPrefix = "mc_mass");

void mc_mass_fit(const char *configFile = "cfgs/pp_run2_pt6p5_9_y0_1p6.C")
{
  gROOT->ProcessLineSync(Form(".x %s", configFile));
  FitConfigMcMass cfg = getConfigMcMass();

  cout << "\n === Start mc_mass_fit() ===\n";
  cout << "Config: " << configFile << "\n"
       << "System: " << cfg.system << "\n"
       << "Mode: " << cfg.mode << "\n"
       << "Prompt: " << (cfg.isPrompt ? "true" : "false") << "\n";


  TStopwatch t; t.Start();

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(WARNING);

  // use rootlogon
  gROOT->Macro(cfg.rootlogon.c_str());

  // make output dirs
  gSystem->mkdir(cfg.figDir.c_str(), true);
  gSystem->mkdir(cfg.rootDir.c_str(), true);
  std::cout << "[Info] Output dirs prepared: " << cfg.figDir << ", " << cfg.rootDir << "\n";

  // read input
  TFile *fInput = TFile::Open(cfg.inputFile.c_str());
  if (!fInput || fInput->IsZombie())
  {
    cout << "[Error] Cannot open input file\n";
    return;
  }
  // OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root
  // OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root
  // OniaRooDataSet_JPsi_pp_GENONLY_NonPrompt_230215.root


  // read dataset
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!dsRaw)
  {
    cout << "[Error]: Cannot find RooDataSet\n";
    return;
  }

  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0, "weight");

  // define cut
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form( "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f)", cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh, cfg.massLow, cfg.massHigh);

  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim later
  TString osCut = "(recoQQsign == 0)";

  TString fullCut = Form("%s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data());

  // === new dataset with cuts ===
  RooDataSet *dsReduced_tmp = (RooDataSet *)dsWeight->reduce(Cut(fullCut));
  if (!dsReduced_tmp || dsReduced_tmp->numEntries() == 0)
  {
    cout << "[Error] Reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced_tmp->Print();

  // check weight
  const bool hasWeight = (dsReduced_tmp->isWeighted() || dsReduced_tmp->weightVar() != nullptr);
  if (hasWeight) {
    // MC's weight value = 1 -> same with no-weight case
    cout << "[Info] Using weighted dataset \n";
  }
  else {
    cout << "[Info] Using UN-weighted dataset \n";
  }


  // observable
  auto mass = new RooRealVar("mass", "invariant mass", cfg.massLow , cfg.massHigh, "GeV/c^{2}");
  RooArgSet obs(*mass);
  auto dsReduced = new RooDataSet("dsReduced", "dataset with local vars", obs, Import(*dsReduced_tmp));


  // === Create workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*dsReduced);


  // === build model ===
  defineMassSig(ws, cfg);


  // === perform fit ===
  cout << "\n--- start fitting --- \n";
  auto fitResult = ws->pdf("model")->fitTo(*ws->data("dsReduced"), Save(), SumW2Error(hasWeight), Offset(true), PrintLevel(-1), NumCPU(cfg.numCPU), EvalBackend("legacy"), PrintEvalErrors(-1), Verbose(false));

  // === draw plots ===
  cout << "\n--- draw plots --- \n";
  draw_mass_fit(ws, *ws->var("mass"), *ws->pdf("model"), *dsReduced, cfg, fitResult, 80, "mc_mass");


  // === save results ===
  cout << "\n--- save results ---\n";
  TFile fout(Form("%s/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "RECREATE");
  fitResult->Write("fitResult");
  ws->pdf("model")->Write();
  fout.Close();

  fitResult->Print("V");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}


// === Helper functions ===
// ------------------------
void defineMassSig(RooWorkspace *ws, FitConfigMcMass &cfg)
{
  cout << "\n--- build mass signal  if (cfg.massSigPdf == "DCBGauss") {
    // cb1
    ws->factory(Form(
        "CBShape::cbLeft(mass,"
        " mean[%g,%g,%g],"
        " sigmaL[%g,%g,%g],"
        " alphaL[%g,%g,%g],"
        " nL[%g,%g,%g])",
        cfg.mean.init, cfg.mean.min, cfg.mean.max,
        cfg.sigmaL.init, cfg.sigmaL.min, cfg.sigmaL.max,
        cfg.alphaL.init, cfg.alphaL.min, cfg.alphaL.max,
        cfg.nL.init, cfg.nL.min, cfg.nL.max));

    // cb2
    ws->factory(Form("sigmaRatio[%g,%g,%g]", cfg.sigmaRatio.init, cfg.sigmaRatio.min, cfg.sigmaRatio.max));
    ws->factory(Form("alphaRatio[%g,%g,%g]", cfg.alphaRatio.init, cfg.alphaRatio.min, cfg.alphaRatio.max));
    ws->factory(Form("nRatio[%g,%g,%g]", cfg.nRatio.init, cfg.nRatio.min, cfg.nRatio.max));

    ws->factory("expr::sigmaR('@0*@1',{sigmaL,sigmaRatio})");
    ws->factory("expr::alphaR('-@0*@1',{alphaL,alphaRatio})");
    ws->factory("expr::nR('@0*@1',{nL,nRatio})");

    ws->factory("CBShape::cbRight(mass, mean, sigmaR, alphaR, nR)");

    // gauss
    ws->factory(Form("sigmaGRatio[%g,%g,%g]", cfg.sigmaGRatio.init, cfg.sigmaGRatio.min, cfg.sigmaGRatio.max));
    ws->factory("expr::sigmaG('@0*@1',{sigmaL,sigmaGRatio})");
    ws->factory("Gaussian::gaus(mass, mean, sigmaG)");

    // sum
    ws->factory(Form("f_cbR[%g,%g,%g]", cfg.f_cbR.init, cfg.f_cbR.min, cfg.f_cbR.max));
    ws->factory(Form("f_gaus[%g,%g,%g]", cfg.f_gaus.init, cfg.f_gaus.min, cfg.f_gaus.max));
    ws->factory("SUM::model(f_gaus*gaus, f_cbR*cbRight, cbLeft)");

    cout << "mass sig pdf: DCBGauss\n"; model ---\n";
  // select function according to cfg setting

  }
  else if (cfg.massSigPdf == "DCBGauss") {
    // cb1
    ws->factory(Form(
        "CBShape::cbLeft(mass,"
        " mean[%g,%g,%g],"
        " sigmaL[%g,%g,%g],"
        " alphaL[%g,%g,%g],"
        " nL[%g,%g,%g])",
        cfg.mean.init, cfg.mean.min, cfg.mean.max,
        cfg.sigmaL.init, cfg.sigmaL.min, cfg.sigmaL.max,
        cfg.alphaL.init, cfg.alphaL.min, cfg.alphaL.max,
        cfg.nL.init, cfg.nL.min, cfg.nL.max));

    // cb2
    ws->factory(Form("sigmaRatio[%g,%g,%g]", cfg.sigmaRatio.init, cfg.sigmaRatio.min, cfg.sigmaRatio.max));
    ws->factory(Form("alphaRatio[%g,%g,%g]", cfg.alphaRatio.init, cfg.alphaRatio.min, cfg.alphaRatio.max));
    ws->factory(Form("nRatio[%g,%g,%g]", cfg.nRatio.init, cfg.nRatio.min, cfg.nRatio.max));

    ws->factory("expr::sigmaR('@0*@1',{sigmaL,sigmaRatio})");
    ws->factory("expr::alphaR('-@0*@1',{alphaL,alphaRatio})");
    ws->factory("expr::nR('@0*@1',{nL,nRatio})");

    ws->factory("CBShape::cbRight(mass, mean, sigmaR, alphaR, nR)");

    // sum
    ws->factory(Form("f_cbR[%g,%g,%g]", cfg.f_cbR.init, cfg.f_cbR.min, cfg.f_cbR.max));
    ws->factory("SUM::model(f_gaus*gaus, f_cbR*cbRight, cbLeft)");
    
    cout << "mass sig pdf: DCB\n";
  }
  else {
    cerr << "[Error]: Mass signal model is wrong: " << cfg.massSigPdf << "\n";
    cerr << "[Error]: You can choose: DCBGauss, DCB\n";
    return;
  }
}

void draw_mass_fit(RooWorkspace *ws, RooRealVar &mass, RooAbsPdf &model, RooAbsData &data, const FitConfigMcMass &cfg, const RooFitResult *fitResult, int nBins, const char *outPrefix)
{
  // --- small local helpers ---
  auto hasPdf = [&](const char *n) -> bool
  { return ws && ws->pdf(n); };
  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  { return fr ? fr->findObject(n) : nullptr; };
  const bool wantGaus = (cfg.massSigPdf == "DCBGauss");

  // --- canvas & top pad ---
  TCanvas c("c_mass", "c_mass", 800, 800);
  TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1.SetBottomMargin(0.00001);
  pad1.SetLogy();
  pad1.Draw();
  pad1.cd();

  // --- frame & plot ---
  const double massMin = cfg.massLow;
  const double massMax = cfg.massHigh;
  mass.setRange(massMin, massMax);

  RooPlot *fr = mass.frame(Range(massMin, massMax), Bins(nBins), Title(""));
  data.plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"));

  // component
  model.plotOn(fr, Components("cbLeft"), LineStyle(kDotted), LineColor(kRed), Name("cbLeft"));
  model.plotOn(fr, Components("cbRight"), LineStyle(kDotted), LineColor(kAzure), Name("cbRight"));
  if (wantGaus && hasPdf("gaus"))
  {
    model.plotOn(fr, Components("gaus"), LineStyle(kDotted), LineColor(kViolet), Name("gaus"));
  }

  // --- dynamic y-range for log scale ---
  double ymin = 1e300, ymax = -1e300;
  if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
  {
    for (int i = 0; i < hdata->GetN(); ++i)
    {
      double x, y;
      hdata->GetPoint(i, x, y);
      if (y > 0 && y < ymin)
        ymin = y;
      if (y > ymax)
        ymax = y;
    }
  }
  if (ymin <= 0 || ymin == 1e300)
    ymin = 1e-3;
  fr->SetMinimum(ymin * 0.5);
  fr->SetMaximum(std::max(ymax, ymin) * 1e6);

  fr->GetYaxis()->SetTitle("Events");
  fr->GetXaxis()->SetTitle("");
  fr->Draw("e");

  // --- legend ---
  TLegend leg(0.49, 0.65, 0.70, 0.94);
  {
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    if (auto *o = findObj(fr, "data"))
      leg.AddEntry(o, "Data", "lep");
    if (auto *o = findObj(fr, "model"))
      leg.AddEntry(o, "Fit model", "pe");
    if (auto *o = findObj(fr, "cbLeft"))
      leg.AddEntry(o, "CB1", "pe");
    if (auto *o = findObj(fr, "cbRight"))
      leg.AddEntry(o, "CB2", "pe");
    if (wantGaus)
    {
      if (auto *o = findObj(fr, "gaus"))
        leg.AddEntry(o, "Gauss", "pe");
    }
    leg.Draw("same");
  }

  // --- CMS/info latex ---
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    double x = 0.19, y0 = 0.90, dy = -0.06;
    int k = 0;
    tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
    tx.DrawLatex(x, y0 + dy * k++, "Prompt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
    if (cfg.yLow == 0)
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", cfg.ptLow, cfg.ptHigh, cfg.yHigh));
    else
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh));
  }

  // --- parameter latex ---
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.025);
    tp.SetTextFont(42);
    double x = 0.71, y0 = 0.91, dy = -0.045;
    int k = 0;
    auto print = [&](const char *title, const char *vname)
    {
      if (auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname)))
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, v->getVal(), v->getError()));
    };
    print("mean", "mean");
    print("#alpha_{L}", "alphaL");
    print("n_{L}", "nL");
    print("#sigma_{L}", "sigmaL");
    print("#alpha_{ratio}", "alphaRatio");
    print("n_{ratio}", "nRatio");
    print("#sigma_{ratio,R}", "sigmaRatio");
    if (wantGaus)
    {
      print("#sigma_{ratio,Gauss}", "sigmaGRatio");
      print("f_{Gauss}", "f_gaus");
    }
    print("f_{CB,R}", "f_cbR");
  }

  // --- pull pad ---
  c.cd();
  TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2.SetTopMargin(0.00001);
  pad2.SetBottomMargin(0.4);
  pad2.Draw();
  pad2.cd();

  RooHist *hpull = fr->pullHist("data", "model");
  RooPlot *fpull = mass.frame(Range(massMin, massMax), Title(""));
  fpull->addPlotable(hpull, "P");
  fpull->GetYaxis()->SetTitle("Pull");
  fpull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
  fpull->GetXaxis()->CenterTitle();
  fpull->SetMinimum(-8);
  fpull->SetMaximum(8);
  fpull->GetYaxis()->SetNdivisions(505);
  fpull->GetYaxis()->SetTitleSize(0.12);
  fpull->GetYaxis()->SetLabelSize(0.10);
  fpull->GetXaxis()->SetTitleSize(0.15);
  fpull->GetXaxis()->SetLabelSize(0.10);
  fpull->Draw();

  TLine line(massMin, 0.0, massMax, 0.0);
  line.SetLineStyle(2);
  line.Draw("same");

  // --- chi2/ndf ---
  if (fitResult)
  {
    int npar = fitResult->floatParsFinal().getSize();
    double chi2ndf = fr->chiSquare("model", "data", npar);
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.10);
    tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
    cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
  }

  // --- save ---
  // gSystem->mkdir(cfg.figDir.c_str(), true);
  TString out = Form("%s/%s_pT%.1f_%.1f_y%.1f_%.1f",
                     cfg.figDir.c_str(), outPrefix,
                     cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh);
  c.SaveAs(out + ".png");
  c.SaveAs(out + ".pdf");
}