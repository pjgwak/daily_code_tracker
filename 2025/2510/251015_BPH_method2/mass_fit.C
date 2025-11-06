#include <iostream>
#include <string>
#include <RooRealVar.h>
#include <TFile.h>
#include <RooDataSet.h>
#include <TStopwatch.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <RooChebychev.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TLatex.h>
#include <TAxis.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>

using std::cout; using std::string;
using namespace RooFit;

void mass_fit(double ptLow=6.5, double ptHigh=7.5, double ctMin=0.5, double ctMax=0.6)
{
  cout << "=== start mass_fit() ===\n";
  TStopwatch time;
  time.Start();

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(WARNING);

  // --- set variables ---
  // kinematics
  float massLow = 2.6, massHigh = 3.5;
  float yLow = 0, yHigh = 2.4;
  int cLow = 0, cHigh = 180;

  TString figDir = Form("figs/pT%.1f_%.1f_y%.1f_%.1f", ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots/pT%.1f_%.1f_y%.1f_%.1f", ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);

  // --- read input ---
  cout << "\n=== Import inputs ===\n";
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << "\n";
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");
  
  // --- set kinematics ---
  cout << "\n=== Set kienematic cuts ===\n";
  char reduceDS_woCtErr[3000];

  sprintf(reduceDS_woCtErr, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f)"

          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);

  // reduce the dataset
  cout << "\n=== Reduce datasets ===\n";
  RooDataSet *redData;

  
  redData = (RooDataSet *)data->reduce(reduceDS_woCtErr);

  // --- build model ---
  cout << "\n=== Build model ===\n";
  
  // signal
  const RooArgSet *row = redData->get();
  RooRealVar *mass = const_cast<RooRealVar *>(static_cast<const RooRealVar *>(row->find("mass")));
  mass->setRange(2.6, 3.5);
  mass->setRange("massRange", 2.6, 3.5);

  RooRealVar mean("mean", "mean", 3.096, 3.04, 3.13);
  RooRealVar sigma1("sigma1", "CB sigma", 0.030, 0.005, 0.100);
  RooRealVar sigma12("sigma12", "CB sigma", 1.1, 1, 10);
  RooFormulaVar sigma2("sigma2", "@0*@1", {sigma1, sigma12});
  RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.2, 5.0);

  // RooRealVar nL("nL", "nL", 1.3540);
  RooRealVar nL("nL", "nL", 3.0, 1.0, 100.0);
  RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2);
  RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);

  RooRealVar sigma1G("sigma1G", "CB sigma", 1.1, 1, 10);
  RooFormulaVar sigmaG("sigmaG", "@0*@1", {sigma1, sigma1G});
  ;
  RooRealVar fRight("fRight", "frac(CBR in DCB)", 0.50, 0.00, 1.00);
  RooRealVar fG("fG", "frac(Gauss in total)", 0.20, 0.00, 1.00);

  RooCBShape CBL("CBL", "left-tail CB", *mass, mean, sigma1, alphaL, nL);
  RooCBShape CBR("CBR", "right-tail CB", *mass, mean, sigma2, alphaL, nL);
  RooAddPdf DCB("DCB", "Double-CB approx",
                RooArgList(CBR, CBL),
                RooArgList(fRight), // fRight * CBR + (1-fRight) * CBL
                true);

  RooGaussian G("G", "core Gauss", *mass, mean, sigmaG);
  RooAddPdf sigMass("sigMass", "Signal DCB+G",
                RooArgList(CBR, CBL, G),
                RooArgList(fRight, fG), // fG * G + (1-fG) * DCB
                true);

  // bkg
  RooRealVar sl1("sl1", "sl1", 0.01, -1.0, 1.0);
  RooRealVar sl2("sl2", "sl2", 0.01, -1.0, 1.0);
  RooRealVar sl3("sl3", "sl3", 0.01, -1.0, 1.0);
  RooChebychev bkgMass("bkgMass", "Chebyshev order 3", *mass, RooArgList(sl1, sl2));

  // total model
  double sigInit, sigLow, sigHigh;
  double bkgInit, bkgLow, bkgHigh;

  if (ctMin >= -0.1 && ctMax < 0.1) {
    sigInit = 10000;
    sigLow = 1;
    sigHigh = 1e8;
    bkgInit = 1000;
    bkgLow = 1;
    bkgHigh = 1e6;
  } else {
    sigInit = 1000;
    sigLow = 1;
    sigHigh = 1e6;
    bkgInit = 100;
    bkgLow = 1;
    bkgHigh = 1e4;
  }
  RooRealVar nSig("nSig", "signal yield", sigInit, sigLow, sigHigh);
  RooRealVar nBkg("nBkg", "background yield", bkgInit, bkgLow, bkgHigh);

  // 신호+배경 합
  RooAddPdf model("model", "sigMass + Cheby3 (extended)",
                  RooArgList(sigMass, bkgMass),
                  RooArgList(nSig, nBkg));


  // --- fix from MC ---
  TFile inputMcFit(Form("%s/mc_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
  RooFitResult *fitMcMass = dynamic_cast<RooFitResult *>(inputMcFit.Get("fitMass"));
  RooArgSet targets;
  targets.add(nL);
  targets.add(fG);
  targets.add(sigma1G);
  targets.add(sigma12);
  // set MC value and fix

  auto apply_from_list = [&](const RooArgList &lst)
  {
    for (int i = 0; i < lst.getSize(); ++i)
    {
      const RooAbsArg *a = lst.at(i);
      const RooRealVar *src = dynamic_cast<const RooRealVar *>(a);
      if (!src)
        continue;
      if (RooAbsArg *destA = targets.find(src->GetName()))
      {
        if (auto *dest = dynamic_cast<RooRealVar *>(destA))
        {
          dest->setVal(src->getVal());
          dest->setConstant(true); //
          // std::cout << "[MC->DATA] " << src->GetName() << " = " << src->getVal() << " (frozen)\n";
        }
      }
    }
  };
  apply_from_list(fitMcMass->floatParsFinal());
  apply_from_list(fitMcMass->constPars());

  // --- fit ---
  RooFitResult *fitMass = model.fitTo(*redData, Range("massRange"), Save(), Extended(), PrintLevel(0), NumCPU(32), EvalBackend("legacy"));

  // === start plot ===
  // --- canvas & top pad ---
  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  { return fr ? fr->findObject(n) : nullptr; };

  TCanvas c("c_mass", "c_mass", 800, 800);
  TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1.SetBottomMargin(0.00001);
  pad1.SetLogy();
  pad1.Draw();
  pad1.cd();

  // --- frame & plot ---
  RooPlot *fr = mass->frame(Range(massLow, massHigh), Title("")); // , Bins(nBins)
  redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"));
  model.plotOn(fr, Components("sigMass"), LineStyle(kDotted), LineColor(kRed), Name("sigMass"));
  model.plotOn(fr, Components("bkgMass"), LineStyle(kDotted), LineColor(kAzure), Name("bkgMass"));

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
  fr->SetMaximum(std::max(ymax, ymin) * 1e2);

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
    if (auto *o = findObj(fr, "sigMass"))
      leg.AddEntry(o, "Signal", "pe");
    if (auto *o = findObj(fr, "bkgMass"))
      leg.AddEntry(o, "Bkg", "pe");

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
    tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
    if (yLow == 0)
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
    tx.DrawLatex(x, y0 + dy * k++, Form("%.2f < c#tau < %.2f", ctMin, ctMax));
    
    // fit status
    int st = fitMass->status();  // 0 = success
    int cq = fitMass->covQual(); // 3 = best quality
    if (st != 0 || cq != 3) {
      tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d  CovQual=%d", st, cq));
      std::ofstream flog(Form("logs/fit_status_pT%.1f_%.1f_ctau%.2f_%.2f.txt", ptLow, ptHigh, ctMin, ctMax), std::ios::app);
      flog.close();
    }
    
  }

  // --- parameter latex ---
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.024);
    tp.SetTextFont(42);
    double x = 0.71, y0 = 0.91, dy = -0.04;
    int k = 0;
    auto print = [&](const char *title, const char *vname)
    {
      if (auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname)))
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, v->getVal(), v->getError()));
    };

    // print
    print("N_{Sig}", "nSig");
    print("N_{Bkg}", "nBkg");

    print("mean", "mean");
    print("#alpha_{L}", "alphaL");
    print("n_{L}", "nL");
    print("#sigma_{1}", "sigma1");
    print("#sigma_{2/1}", "sigma12");
    print("#sigma_{G/1}", "sigma1G");
    print("f_{G}", "fG");
    print("f_{CB,R}", "fRight");

    print("#tau", "tauMass");
    print("sl_{1}", "sl1");
    print("sl_{2}", "sl2");
    print("sl_{3}", "sl3");
    print("sl_{4}", "sl4");
    print("sl_{5}", "sl5");
    print("sl_{6}", "sl6");
  }

  // --- pull pad ---
  c.cd();
  TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2.SetTopMargin(0.00001);
  pad2.SetBottomMargin(0.4);
  pad2.Draw();
  pad2.cd();

  RooHist *hpull = fr->pullHist("data", "model");
  RooPlot *fpull = mass->frame(Range(massLow, massHigh), Title(""));
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

  TLine line(massLow, 0.0, massHigh, 0.0);
  line.SetLineStyle(2);
  line.Draw("same");

  // --- chi2/ndf ---
  if (fitMass)
  {
    int npar = fitMass->floatParsFinal().getSize();
    double chi2ndf = fr->chiSquare("model", "data", npar);
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.10);
    tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
    cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
  }
  c.SaveAs(Form("%s/mass_fit_pT%.1f_%.1f_ct%.2f_%.2f.png", figDir.Data(), ptLow, ptHigh, ctMin, ctMax));
  // === finish plot ===

  // print fit result
  fitMass->Print("V");

  // --- save ---
  TFile f(Form("%s/mass_fit_pT%.1f_%.1f_ct%.2f_%.2f.root", rootDir.Data(), ptLow, ptHigh, ctMin, ctMax), "RECREATE");
  fitMass->Write("fitMass");
  f.Close();

  cout << "\n=== finish mass_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}