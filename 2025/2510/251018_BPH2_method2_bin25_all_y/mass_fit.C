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
#include <TROOT.h>   // gROOT
#include <TSystem.h> // gSystem
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <TLine.h>
#include <fstream> // ofstream
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>

using std::cout; using std::string;
using namespace RooFit;

void mass_fit(double ptLow = 7.5, double ptHigh = 8.5, double ctMin = -0.5, double ctMax=-0.4)
{
  // double binWidth = 0.05;
  //   = ctMin + binWidth;
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

  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
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
  RooRealVar sigma12("sigma12", "ratio of sigma 1 vs 2", 1.1, 0.1, 10);
  RooFormulaVar sigma2("sigma2", "@0*@1", {sigma1, sigma12});

  // CBs share nL and nR values
  RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.01, 10.0);
  RooRealVar nL("nL", "nL", 3.0, 1.0, 10.0);
  // RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2)
  // RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);

  RooRealVar sigma1G("sigma1G", "ratio of sigma1 vs sigmaG", 1.1, 1, 105);
  RooFormulaVar sigmaG("sigmaG", "@0*@1", {sigma1, sigma1G});

  RooRealVar fCB2("fCB2", "frac(CB2 in DCB)", 0.50, 0.00, 1.00);
  RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.20, 0.00, 1.00);

  RooCBShape CB1("CB1", "left-tail CB", *mass, mean, sigma1, alphaL, nL);
  RooCBShape CB2("CB2", "right-tail CB", *mass, mean, sigma2, alphaL, nL);

  // Gauss
  RooGaussian G("G", "core Gauss", *mass, mean, sigmaG);

  RooAddPdf sigMass("sigMass", "Signal DCB+G",
                    RooArgList(CB1, CB2, G),
                    RooArgList(fCB1, fCB2),
                    true); // Recursive fracton

  // bkg
  RooRealVar sl1("sl1", "sl1", 0.05, -1.0, 1.0);
  RooRealVar sl2("sl2", "sl2", 0.05, -1.0, 1.0);
  RooRealVar sl3("sl3", "sl3", 0.01, -1.0, 1.0);
  RooChebychev bkgMass("bkgMass", "", *mass, RooArgList(sl1, sl2));

  // total model
  double sigInit, sigLow, sigHigh;
  double bkgInit, bkgLow, bkgHigh;

  double totalEvent = redData->sumEntries();
  sigInit = totalEvent*0.7;
  sigLow = 1;
  sigHigh = totalEvent;
  bkgInit = totalEvent*0.3;
  bkgLow = 1;
  bkgHigh = totalEvent;

  // if (ctMin >= -0.1 && ctMax < 0.1) {
  //   sigInit = 10000;
  //   sigLow = 1;
  //   sigHigh = 1e8;
  //   bkgInit = 1000;
  //   bkgLow = 1;
  //   bkgHigh = 1e6;
  // } else {
  //   sigInit = 1000;
  //   sigLow = 1;
  //   sigHigh = 1e6;
  //   bkgInit = 100;
  //   bkgLow = 1;
  //   bkgHigh = 1e4;
  // }

  // sigInit = 100;
  // sigLow = 1;
  // sigHigh = 1000;

  // bkgInit = 100;
  // bkgLow = 1;
  // bkgHigh = 1000;

  RooRealVar nSig("nSig", "signal yield", sigInit, sigLow, sigHigh);
  RooRealVar nBkg("nBkg", "background yield", bkgInit, bkgLow, bkgHigh);

  RooAddPdf model("model", "",
                  RooArgList(sigMass, bkgMass),
                  RooArgList(nSig, nBkg));

  // --- fix from MC ---
  TFile inputMcFit(Form("%s/mc_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
  RooFitResult *fitMcMass = dynamic_cast<RooFitResult *>(inputMcFit.Get("fitMass"));
  RooArgSet targets; // fixed parameter list
  targets.add(nL);
  targets.add(fCB2);
  targets.add(sigma1G);
  targets.add(sigma12);

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
  RooFitResult *fitMass = model.fitTo(*redData, Range("massRange"), Save(), Extended(), PrintLevel(-1), Verbose(false), PrintEvalErrors(-1), NumCPU(32), EvalBackend("legacy"));

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
  fr->SetMaximum(std::max(ymax, ymin) * 1e4);

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
    tx.DrawLatex(x, y0 + dy * k++, Form("%.3f < c#tau < %.3f", ctMin, ctMax));
    
    // fit status
    int st = fitMass->status();  // 0 = success
    if (st != 0) {
      tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
      std::ofstream flog(Form("logs%s/mass_status_pT%.1f_%.1f_ctau%.3f_%.3f.txt", userLabel.Data(),  ptLow, ptHigh, ctMin, ctMax), std::ios::app);
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

    // lambda function for printing
    auto print = [&](const char *title, const char *vname)
    {
      auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname));
      if (!v)
        return;

      const double val = v->getVal(), err = v->getError();
      const bool fixed = v->isConstant();

      const bool hasMin = v->hasMin(), hasMax = v->hasMax();
      const double vmin = hasMin ? v->getMin() : 0.0;
      const double vmax = hasMax ? v->getMax() : 0.0;
      const double span = hasMin && hasMax ? std::fabs(vmax - vmin) : 0.0;

      // tol: 1e-12, scale: 1e-9, range scale: 1e-6
      const double tol = std::max({1e-12,
                                   1e-9 * (1.0 + std::max({std::fabs(val), std::fabs(vmin), std::fabs(vmax)})),
                                   (span > 0 ? 1e-6 * span : 0.0)});

      const bool atMin = hasMin && (val - vmin) <= tol;
      const bool atMax = hasMax && (vmax - val) <= tol;
      const bool atBound = atMin || atMax;

      TString note;
      if (fixed)
        note += "(fixed)";
      if (atBound)
        note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                              : (atMin ? ", (at min)" : ", (at max)");
      if (!fixed && (err <= 0 || std::isnan(err)))
        note += note.IsNull() ? "(no error)" : ", (no error)";

      const Int_t oldColor = tp.GetTextColor();
      if (atBound)
        tp.SetTextColor(kRed + 1);

      if (fixed)
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.6g %s", title, val, note.Data()));
      else
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.6g #pm %.3g %s", title, val, err, note.Data()));

      tp.SetTextColor(oldColor);
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
    print("f_{CB1}", "fCB1");
    print("f_{CB2}", "fCB2");

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
    tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.3f", chi2ndf));
    cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
  }
  c.SaveAs(Form("%s/mass_fit_pT%.1f_%.1f_ct%.3f_%.3f.png", figDir.Data(), ptLow, ptHigh, ctMin, ctMax));
  // === finish plot ===

  // print fit result
  fitMass->Print("V");

  // --- save ---
  TFile f(Form("%s/mass_fit_pT%.1f_%.1f_ct%.3f_%.3f.root", rootDir.Data(), ptLow, ptHigh, ctMin, ctMax), "RECREATE");
  fitMass->Write("fitMass");
  f.Close();

  cout << "\n=== finish mass_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}