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
#include <TROOT.h> // gROOT
#include <TSystem.h> // gSystem
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <TLine.h>
#include <fstream> // ofstream


using std::cout;
using std::string;
using namespace RooFit;

void ctau_res(double ptLow = 6.5, double ptHigh = 7.5)
{
  cout << "=== start ctau_res() ===\n";
  TStopwatch time;
  time.Start();

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // silent mode
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

  // --- make output folders ---
  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);
  gSystem->mkdir(Form("logs%s", userLabel.Data()), kTRUE);

  // --- read input ---
  cout << "\n=== Import inputs ===\n";
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << "\n";
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");

  // --- cuts ---
  cout << "\n=== Set kienematic cuts ===\n";
  char reduceDS_woCtErr[3000];

  sprintf(reduceDS_woCtErr,
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f)"

          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

  // reduce the dataset
  cout << "\n=== Reduce datasets ===\n";
  RooDataSet *redData = (RooDataSet *)data->reduce(reduceDS_woCtErr);


  // --- build model ---
  cout << "\n=== Build model ===\n";
  const RooArgSet *row = redData->get();
  RooRealVar *ctau3DRes = const_cast<RooRealVar *>(static_cast<const RooRealVar *>(row->find("ctau3DRes")));
  ctau3DRes->setRange(-10, 10);
  ctau3DRes->setRange("resRange", -10, 10);

  // --- signal ---
  RooRealVar mu("mu", "mean", 0.0, -0.05, 0.05);
  RooRealVar sigma1("sigma1","scale1",1.00,0.05, 5.0);

  RooRealVar sigma12("sigma12", "", 1.1, 1, 10);
  RooRealVar sigma23("sigma23", "", 1.1, 1, 10);
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, sigma12));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, sigma23));

  RooGaussian G1("G1","G1", *ctau3DRes, mu, sigma1);
  RooGaussian G2("G2","G2", *ctau3DRes, mu, sigma2);
  RooGaussian G3("G3","G3", *ctau3DRes, mu, sigma3);

  RooRealVar f1("f1","frac1",0.70,0.0,1.0);
  RooRealVar f2("f2","frac2",0.25,0.0,1.0);

  RooAddPdf model("model", "", RooArgList(G1, G2, G3), RooArgList(f1, f2), true);

  // --- fit ---
  RooFitResult *fitRes = model.fitTo(*redData, Range("resRange"), Save(), PrintLevel(0), NumCPU(32), EvalBackend("legacy"));

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
  RooPlot *fr = ctau3DRes->frame(Title("")); // , Bins(nBins)
  redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"));
  model.plotOn(fr, Components("G1"), Name("G1"), LineStyle(kDotted), LineColor(kRed));
  model.plotOn(fr, Components("G2"), Name("G2"), LineStyle(kDotted), LineColor(kOrange));
  model.plotOn(fr, Components("G3"), Name("G3"), LineStyle(kDotted), LineColor(kGreen));

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
  TLegend leg(0.49, 0.66, 0.70, 0.94);
  {
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    if (auto *o = findObj(fr, "data"))
      leg.AddEntry(o, "Data", "lep");
    if (auto *o = findObj(fr, "model"))
      leg.AddEntry(o, "Model", "pe");
    if (auto *o = findObj(fr, "G1"))
      leg.AddEntry(o, "Gauss1", "pe");
    if (auto *o = findObj(fr, "G2"))
      leg.AddEntry(o, "Gauss2", "pe");
    if (auto *o = findObj(fr, "G3"))
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
    if (yLow == 0)
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
    else
      tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
    
    // fit status
    int st = fitRes->status(); // 0 = success
    if (st != 0)
    {
      tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
      std::ofstream flog(Form("logs%s/ctau_res_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
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
      const double eps = 1e-9 * (1.0 + std::fabs(val));
      const bool fixed = v->isConstant();
      const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
      const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
      const bool atBound = atMin || atMax;

      TString note;
      if (fixed)
        note += "(fixed)";
      if (atBound)
        note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                              : (atMin ? ", (at min)" : ", (at max)");

      const Int_t oldColor = tp.GetTextColor();
      if (atBound)
        tp.SetTextColor(kRed + 1);

      if (fixed)
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
      else
        tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

      tp.SetTextColor(oldColor);
    };

    // print

    print("mean", "mu");
    print("#sigma1", "sigma1");
    print("#sigma_{2/1}", "sigma12");
    print("#sigma_{3/2}", "sigma23");
    print("f_{1}", "f1");
    print("f_{2}", "f2");
  }

  // --- pull pad ---
  c.cd();
  TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2.SetTopMargin(0.00001);
  pad2.SetBottomMargin(0.4);
  pad2.Draw();
  pad2.cd();

  RooHist *hpull = fr->pullHist("data", "model");
  RooPlot *fpull = ctau3DRes->frame(Title(""));
  fpull->addPlotable(hpull, "P");
  fpull->GetYaxis()->SetTitle("Pull");
  fpull->GetXaxis()->SetTitle("ctau3DRes");
  fpull->GetXaxis()->CenterTitle();
  fpull->SetMinimum(-8);
  fpull->SetMaximum(8);
  fpull->GetYaxis()->SetNdivisions(505);
  fpull->GetYaxis()->SetTitleSize(0.12);
  fpull->GetYaxis()->SetLabelSize(0.10);
  fpull->GetXaxis()->SetTitleSize(0.15);
  fpull->GetXaxis()->SetLabelSize(0.10);
  fpull->Draw();

  // TLine line(massLow, 0.0, massHigh, 0.0);
  // line.SetLineStyle(2);
  // line.Draw("same");

  // --- chi2/ndf ---
  if (fitRes)
  {
    int npar = fitRes->floatParsFinal().getSize();
    double chi2ndf = fr->chiSquare("model", "data", npar);
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.10);
    tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
    cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
  }
  c.SaveAs(Form("%s/ctau_res_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  // === finish plot ===

  // print fit result
  fitRes->Print("V");

  // --- save ---
  TFile f(Form("%s/ctau_res_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitRes->Write("fitRes");
  f.Close();

  cout << "\n=== finish ctau_res() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}