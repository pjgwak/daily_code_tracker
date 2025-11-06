#include "utils/ctau_fit_utils.h"
#include "cfgs/pp_run2_pt6p5_9_y0_1p6.C"
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

  // RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0, "weight");
  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0);

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

  // print out object list
  const char *wname = ""; // 가중치 쓰면 "weight"

  // 2) dsReduced_tmp 안에서 ct 변수 '그 자체'를 꺼내옵니다.
  RooRealVar *ct_in = dynamic_cast<RooRealVar *>(dsReduced_tmp->get()->find("ctau3D"));
  RooRealVar *cte_in = dynamic_cast<RooRealVar *>(dsReduced_tmp->get()->find("ctau3DErr"));

  // 4) 범위는 '데이터가 실제로 가진' 값으로 설정
  double xmin, xmax;
  dsReduced_tmp->getRange(*ct_in, xmin, xmax); // 데이터 기반
  ct_in->setRange(-0.15, 0.15);
  ct_in->setRange("fitRange", -0.15, 0.15);
  cte_in->setRange(0.01, 0.3);
  cte_in->setRange(0.01, 0.3);

  // 5) 2-가우스 모델 (변수는 '새로 만들지 않고' ct_in을 그대로 사용)
  RooRealVar mu("mu", "mean", 0.0, -0.2, 0.2);
  RooRealVar s1("sigma1", "sigma1", 0.08, 0.002, 5);
  RooRealVar s2("sigma2", "sigma2", 0.25, 0.01, 2.0);
  RooRealVar f1("f1", "frac g1", 0.6, 0.0, 1.0);
  RooRealVar one("one", "frac g1", 1.0);

  RooGaussModel g1("g1", "g1", *ct_in, mu, s1, one, *cte_in);
  RooGaussModel g2("g2", "g2", *ct_in, mu, s2, one, *cte_in);
  RooAddPdf model("model", "g1+g2", RooArgList(g1, g2), RooArgList(f1));

  cout << "\npin 1\n";

  // 6) 피팅 (옵션 최소화 → 먼저 움직이는지 확인)
  auto fr = model.fitTo(*dsReduced_tmp, Save(), Range("fitRange"), Offset(true), PrintLevel(1), ConditionalObservables(RooArgSet(*cte_in)));
  
  cout << "\npin 2\n";

  // 7) 플롯 (불필요한 ProjWData 절대 X)
  TCanvas c1("c1", "", 800, 800);
  RooPlot *frame = ct_in->frame(Title(""));
  dsReduced_tmp->plotOn(frame, DataError(RooAbsData::SumW2));
  cout << "\npin 3\n";
  model.plotOn(frame, Range("fitRange"), NormRange("fitRange"), ProjWData(*cte_in, *dsReduced_tmp));
  cout << "\npin 4\n";
  // model.plotOn(frame, Components(g1), LineStyle(kDashed), LineColor(kRed));
  // model.plotOn(frame, Components(g2), LineStyle(kDashed), LineColor(kBlue));

  frame->Draw();
  c1.Draw();
  c1.SaveAs("test.png");

  // 8) 결과 확인
  if (fr)
    fr->Print("v");

  // auto ws = new RooWorkspace("ws");
  // ws->import(*dsReduced_tmp);

  // importErrPdf(ws, cfg);
  // defineCtauPrMc(ws, cfg);
  // auto fr = ws->pdf("resModel")->fitTo(*dsReduced_tmp, Save(), Range("fitRange"), Offset(true), PrintLevel(0), NumCPU(cfg.numCPU), EvalBackend("legacy"), PrintEvalErrors(-1), Verbose(false));

  // fr->Print("V");

  // ------------------

  // // === draw plots ===
  // cout << "\n--- draw plots --- \n";
  // // drawCtauFit(ws, *ws->var("ctau3D"), *ws->pdf("CtPR_PEE"), *dsReduced, cfg, fitResult, 80, "ctau_pr");

  // // === save results ===
  // cout << "\n--- save results ---\n";
  // // TFile fout(Form("%s/ctau_pr_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "RECREATE");
  // // fitResult->Write("fitResult");
  // // ws->pdf("resModel")->Write();
  // // fout.Close();

  // fitResult->Print("V");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}

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
    ws->factory("GaussModel::resModel(ctau3D,mu0,sigma1)");

    cout << "ctau Res pdf: Gauss1\n";
  }
  else if (cfg.ctauResPdf == "Gauss2")
  {
    ws->factory(Form("mu0[%g,%g,%g]", cfg.mu0.init, cfg.mu0.min, cfg.mu0.max));
    ws->factory(Form("sigma1[%g,%g,%g]", cfg.sigma1.init, cfg.sigma1.min, cfg.sigma1.max));
    ws->factory(Form("r21[%g,%g,%g]", cfg.r21.init, cfg.r21.min, cfg.r21.max)); // sigma2/sigma1
    ws->factory(Form("fg1[%g,%g,%g]", cfg.fg1.init, cfg.fg1.min, cfg.fg1.max));

    // --- derived sigmas ---
    ws->factory("expr::sigma2('@0*@1',{sigma1,r21})");

    // --- components
    ws->factory("GaussModel::g1(ctau3D,mu0,sigma1,ctau3DErr)");
    ws->factory("GaussModel::g2(ctau3D,mu0,sigma2,ctau3DErr)");

    // --- 3-Gauss resolution model ---
    ws->factory("RooAddModel::resModel({g1,g2},{fg1})");

    cout << "ctau Res pdf: Gauss2\n";
  }
  else if (cfg.ctauResPdf == "Gauss3")
  {
    ws->factory(Form("mu0[%g,%g,%g]", cfg.mu0.init, cfg.mu0.min, cfg.mu0.max));
    ws->factory(Form("sigma1[%g,%g,%g]", cfg.sigma1.init, cfg.sigma1.min, cfg.sigma1.max));
    ws->factory(Form("r21[%g,%g,%g]", cfg.r21.init, cfg.r21.min, cfg.r21.max)); // sigma2/sigma1
    ws->factory(Form("r32[%g,%g,%g]", cfg.r32.init, cfg.r32.min, cfg.r32.max)); // sigma3/sigma2
    // ws->factory(Form("fg1[%g,%g,%g]", cfg.fg1.init, cfg.fg1.min, cfg.fg1.max));
    // ws->factory(Form("fg2[%g,%g,%g]", cfg.fg2.init, cfg.fg2.min, cfg.fg2.max));

    // // --- derived sigmas ---
    ws->factory("expr::sigma2('@0*@1',{sigma1,r21})");
    ws->factory("expr::sigma3('@0*@1',{sigma2,r32})");

    // // --- components
    ws->factory("GaussModel::g1(ctau3D,mu0,sigma1,ctau3DErr)");
    ws->factory("GaussModel::g2(ctau3D,mu0,sigma2,ctau3DErr)");
    ws->factory("GaussModel::g3(ctau3D,mu0,sigma3,ctau3DErr)");

    // // --- 3-Gauss resolution model ---
    // ws->factory("RooAddModel::resModel({g1,g2,g3},{fg1,fg2})");

    // -------------------
    ws->factory(Form("g1Frac[%g,0,1]", cfg.fg1.init));
    ws->factory(Form("g2FracInRest[%g,0,1]", 0.1)); // fraction of G2 in (1 - fG1)

    // for printing
    ws->factory("expr::fg1('@0',{g1Frac})");                         // fg1 = a
    ws->factory("expr::fg2('(1-@0)*@1',{g1Frac,g2FracInRest})");     // fg2 = (1-a)*b
    ws->factory("expr::fg3('(1-@0)*(1-@1)',{g1Frac,g2FracInRest})"); // fg3 = (1-a)*(1-b)

    ws->factory("RooAddModel::resModel({g1,g2,g3},{fg1,fg2})");
    // ws->factory("RooAddModel::resModel({g1,g2,g3},{fg1,fg2})");

    cout << "ctau Res pdf: Gauss3\n";
  }
  else
  {
    cerr << "[Error]: ctau Res model is wrong: " << cfg.ctauResPdf << "\n";
    cerr << "[Error]: You can choose: Gauss1 ~3\n";
    return;
  }

  RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE",
                      *(ws->pdf("errPdfSig")),
                      Conditional(*(ws->pdf("resModel")), RooArgList(*(ws->var("ctau3D")))));
  ws->import(CtPR_PEE);
}

void importErrPdf(RooWorkspace *ws, FitConfigCtauPrMc &cfg)
{
  TFile fin(Form("%s/ctau_err_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "read");
  if (fin.IsZombie())
  {
    cerr << "[Error] Can't open MC fit file\n";
    return;
  }

  auto h1 = dynamic_cast<RooHistPdf*>(fin.Get("errPdfSig")->Clone("errPdfSig"));
  auto h2 = dynamic_cast<RooHistPdf*>(fin.Get("errPdfSigInter2")->Clone("errPdfSigInter2"));

  if (!h1 || !h2)
  {
    cerr << "[Error] Missing errPdfSig or errPdfSigInter2 in file\n";
    return;
  }

  ws->import(*h1);
  ws->import(*h2);
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
  const double massMin = cfg.ctResLow;
  const double massMax = cfg.ctResHigh;
  mass.setRange(massMin, massMax);

  RooPlot *fr = mass.frame(Range(massMin, massMax), Bins(nBins), Title(""));
  data.plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"));

  // component
  // model.plotOn(fr, Components("g1"), LineStyle(kDotted), LineColor(kRed), Name("g1"));
  // if (cfg.ctauResPdf == "Gauss2" || cfg.ctauResPdf == "Gauss3")
  //   model.plotOn(fr, Components("g2"), LineStyle(kDotted), LineColor(kAzure), Name("g2"));
  // if (cfg.ctauResPdf == "Gauss3")
  //   model.plotOn(fr, Components("g3"), LineStyle(kDotted), LineColor(kViolet), Name("g3"));

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
    print("f_{G1}", "g1Frac");
    print("f_{G2}/(f_{G2}+f_{G3})", "g2FracInRest"); // fraction of G2 in (1 - fG1)
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