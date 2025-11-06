#include "utils/ctau_fit_utils.h"
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

void importErrPdf(RooWorkspace *ws, FitConfigCtauPrMc &cfg);
void defineCtauPrMc(RooWorkspace *ws, FitConfigCtauPrMc &cfg);
void drawCtauFit(RooWorkspace *ws, RooRealVar &ctau3D, RooAbsPdf &model, RooAbsData &data, const FitConfigCtauPrMc &cfg, const RooFitResult *fitResult, int nBins = 80, const char *outPrefix = "ctau_pr");

void ctau_mc_fit(const char *configFile = "cfgs/pp_run2_pt6p5_9_y0_1p6.C")
{
  gROOT->ProcessLineSync(Form(".x %s", configFile));
  FitConfigCtauPrMc cfg = getConfigCtauPrMc();
  FitConfigSB cfgSB = getConfigSB();

  cout << "\n === Start ctau_mc_fit() ===\n";
  cout << "Config: " << configFile << "\n"
       << "System: " << cfg.system << "\n"
       << "Mode: " << cfg.mode << "\n"
       << "Prompt: " << (cfg.isPrompt ? "true" : "false") << "\n";

  TStopwatch t;
  t.Start();

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
  // RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0);

  // define cut
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form("(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
                         " mass >= %.3f && mass < %.3f)",
                         cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh, cfg.massLow, cfg.massHigh);

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

  // check weight
  const bool hasWeight = (dsReduced_tmp->isWeighted() || dsReduced_tmp->weightVar() != nullptr);
  if (hasWeight)
  {
    // MC's weight value = 1 -> same with no-weight case
    cout << "[Info] Using weighted dataset \n";
  }
  else
  {
    cout << "[Info] Using UN-weighted dataset \n";
  }

  auto ctau3D = new RooRealVar("ctau3D", "", cfg.ctLow, cfg.ctHigh, "L_{J/#psi} [mm]");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", cfgSB.ctErrLow, cfgSB.ctErrHigh, "#sigma_{L_{J/#psi}}");
  auto weight = new RooRealVar("weight", "", 0, 10000, "");
  RooArgSet obs(*ctau3D, *ctau3DErr, *weight);
  auto dsReduced = new RooDataSet("dsReduced", "dataset with local vars", obs, Import(*dsReduced_tmp));

  // === bring ctauErrPdf ===
  TFile inputErr(Form("%s/ctau_err_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "read");
  if (inputErr.IsZombie())
  {
    cerr << "[Error] Can't open MC fit file\n";
    return;
  }

  auto errPdfSig = dynamic_cast<RooHistPdf *>(inputErr.Get("errPdfSig")->Clone("errPdfSig"));
  if (!errPdfSig)
  {
    cerr << "[Error] Missing errPdfSig or errPdfSigInter2 in file\n";
    return;
  }

  // === create workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*dsReduced);
  ws->import(*errPdfSig);

  // set ctau3D range
  ws->var("ctau3D")->setRange(cfg.ctLow, cfg.ctHigh);
  ws->var("ctau3D")->setRange("fitRange", cfg.ctLow, cfg.ctHigh);
  // cte_in->setRange(0.00, 0.3); // no ctau3DErr cut for MC

  // === buld model ===
  defineCtauPrMc(ws, cfg);

  RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE",
                      *ws->pdf("errPdfSig"),
                      Conditional(*ws->pdf("prResPdf"), *ws->var("ctau3D")));
  ws->import(CtPR_PEE);

  // === perform fit ===
  cout << "\n--- start fitting ---\n";
  auto fitResult = ws->pdf("CtPR_PEE")->fitTo(*dsReduced, Save(), Range("fitRange"), Offset(true), PrintLevel(-1), PrintEvalErrors(-1), SumW2Error(hasWeight), ConditionalObservables(RooArgSet(*ctau3DErr)), NumCPU(cfg.numCPU), EvalBackend("legacy"));

  // === draw plot ===
  drawCtauFit(ws, *ws->var("ctau3D"), *ws->pdf("CtPR_PEE"), *dsReduced, cfg, fitResult, 80, "ctau_pr");
  
  // === save results ===
  cout << "\n--- save results ---\n";
  TFile fout(Form("%s/ctau_pr_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "RECREATE");
  // CtPR_PEE.Write();
  fitResult->Write("fitResult");
  ws->pdf("CtPR_PEE")->Write();
  fout.Close();

  if (fitResult)
    fitResult->Print("V");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}

// 
// 
// 
// 
// 
// 

// === Helper functions ===
// ------------------------
void defineCtauPrMc(RooWorkspace *ws, FitConfigCtauPrMc &cfg)
{
  cout << "\n--- build ctau Res model ---\n";
  // select function according to cfg setting
  if (cfg.ctauResPdf == "Gauss1")
  {
    ws->factory(Form("mu0[%g,%g,%g]", cfg.mu0.init, cfg.mu0.min, cfg.mu0.max));
    ws->factory(Form("sigma1[%g,%g,%g]", cfg.sigma1.init, cfg.sigma1.min, cfg.sigma1.max));

    // --- components ---
    ws->factory("GaussModel::g1(ctau3D,mu0,sigma1,ctau3DErr)");

    // --- 3-Gauss resolution model ---
    // to keep the name consistently
    ws->factory("GaussModel::prResPdf(ctau3D,mu0,sigma1)");

    cout << "ctau Res pdf: Gauss1\n";
  }
  else if (cfg.ctauResPdf == "Gauss2")
  {
    ws->factory(Form("mu0[%g,%g,%g]", cfg.mu0.init, cfg.mu0.min, cfg.mu0.max));
    ws->factory(Form("sigma1[%g,%g,%g]", cfg.sigma1.init, cfg.sigma1.min, cfg.sigma1.max));
    ws->factory(Form("r21[%g,%g,%g]", cfg.r21.init, cfg.r21.min, cfg.r21.max)); // sigma2/sigma1

    ws->factory("expr::sigma2('@0*@1',{sigma1,r21})");

    ws->factory("GaussModel::g1(ctau3D,mu0,sigma1,one[1],ctau3DErr)");
    ws->factory("GaussModel::g2(ctau3D,mu0,sigma2,one,ctau3DErr)");

    ws->factory(Form("fg1[%g,%g,%g]", cfg.fg1.init, cfg.fg1.min, cfg.fg1.max));

    ws->factory("AddPdf::prResPdf({g1,g2},{fg1})"); // AddPdf is better if Resolution model is the fit model
    // ws->factory("RooAddModel::prResPdf({g1,g2},{fg1})"); // Must use AddModel when insert res model into the decay model

    cout << "ctau Res pdf: Gauss2\n";
  }
  else if (cfg.ctauResPdf == "Gauss3")
  {
    ws->factory(Form("mu0[%g,%g,%g]", cfg.mu0.init, cfg.mu0.min, cfg.mu0.max));
    ws->factory(Form("sigma1[%g,%g,%g]", cfg.sigma1.init, cfg.sigma1.min, cfg.sigma1.max));
    ws->factory(Form("r21[%g,%g,%g]", cfg.r21.init, cfg.r21.min, cfg.r21.max)); // sigma2/sigma1
    ws->factory(Form("r32[%g,%g,%g]", cfg.r32.init, cfg.r32.min, cfg.r32.max)); // sigma3/sigma2

    ws->factory("expr::sigma2('@0*@1',{sigma1,r21})");
    ws->factory("expr::sigma3('@0*@1',{sigma2,r32})");

    ws->factory("GaussModel::g1(ctau3D,mu0,sigma1,one[1],ctau3DErr)");
    ws->factory("GaussModel::g2(ctau3D,mu0,sigma2,one,ctau3DErr)");
    ws->factory("GaussModel::g3(ctau3D,mu0,sigma3,one,ctau3DErr)");

    ws->factory(Form("fg1[%g,0,1]", cfg.fg1.init));
    ws->factory(Form("fg23[%g,0,1]", 0.1)); // fraction of G2 in (1 - fG1)

    ws->factory("expr::fg2('(1-@0)*@1',{fg1,fg23})");     // fg2 = (1-a)*b
    ws->factory("expr::fg3('(1-@0)*(1-@1)',{fg1,fg23})"); // fg3 = (1-a)*(1-b)

    ws->factory("AddPdf::prResPdf({g1,g2,g3},{fg1,fg2})"); // AddPdf is better if Resolution model is the fit model
    // ws->factory("RooAddModel::prResPdf({g1,g2,g3},{fg1,fg2})"); // Must use AddModel when insert res model into the decay model

    cout << "ctau Res pdf: Gauss3\n";
  }
  else
  {
    cerr << "[Error]: ctau Res model is wrong: " << cfg.ctauResPdf << "\n";
    cerr << "[Error]: You can choose: Gauss1 ~3\n";
    return;
  }
}

void drawCtauFit(RooWorkspace *ws, RooRealVar &mass, RooAbsPdf &model, RooAbsData &data, const FitConfigCtauPrMc &cfg, const RooFitResult *fitResult, int nBins, const char *outPrefix)
{
  // --- small local helpers ---
  auto hasPdf = [&](const char *n) -> bool
  { return ws && ws->pdf(n); };
  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  { return fr ? fr->findObject(n) : nullptr; };

  // --- canvas & top pad ---
  TCanvas c("c_mass", "c_mass", 800, 800);
  TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1.SetBottomMargin(0.00001);
  pad1.SetLogy();
  pad1.Draw();
  pad1.cd();

  // --- frame & plot ---
  const double massMin = cfg.ctLow;
  const double massMax = cfg.ctHigh;
  mass.setRange(massMin, massMax);

  RooPlot *fr = mass.frame(Range(massMin, massMax), Bins(nBins), Title(""));
  data.plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsReduced")));

  // component
  model.plotOn(fr, Components("g1"), LineStyle(kDotted), LineColor(kRed), Name("g1"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsReduced")), Normalization(ws->var("fg1")->getVal(), RooAbsReal::Relative));
  if (cfg.ctauResPdf == "Gauss2" || cfg.ctauResPdf == "Gauss3")
    model.plotOn(fr, Components("g2"), LineStyle(kDotted), LineColor(kOrange), Name("g2"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsReduced")), Normalization((1-ws->var("fg1")->getVal())*ws->var("fg23")->getVal(), RooAbsReal::Relative));
  if (cfg.ctauResPdf == "Gauss3")
    model.plotOn(fr, Components("g3"), LineStyle(kDotted), LineColor(kGreen), Name("g3"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsReduced")), Normalization((1 - ws->var("fg1")->getVal()) * (1 - ws->var("fg23")->getVal()), RooAbsReal::Relative));

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
    if (auto *o = findObj(fr, "g1"))
      leg.AddEntry(o, "Gauss1", "pe");
    if (auto *o = findObj(fr, "g2"))
      leg.AddEntry(o, "Gauss2", "pe");
    if (auto *o = findObj(fr, "g3"))
      leg.AddEntry(o, "Gauss3", "pe");

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
    print("mean", "mu0");
    print("#sigma_{G1}", "sigma1");
    print("#sigma_{21}", "r21");
    print("#sigma_{32}", "r32");
    print("f_{G1}", "fg1");
    print("f_{G2}/(f_{G2}+f_{G3})", "fg23"); // fraction of G2 in (1 - fG1)
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