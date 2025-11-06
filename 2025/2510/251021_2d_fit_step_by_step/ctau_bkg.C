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
#include <RooGaussModel.h>
#include <RooAddModel.h>
#include <RooDecay.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooArgSet.h>

using std::cout;
using std::string;
using namespace RooFit;

void ctau_bkg(double ptLow = 6.5, double ptHigh = 7.5)
{
  cout << "=== start ctau_bkg() ===\n";
  TStopwatch time;
  time.Start();

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // silent mode
  // RooMsgService::instance().getStream(0).removeTopic(Eval);
  // RooMsgService::instance().getStream(1).removeTopic(Eval);
  // RooMsgService::instance().getStream(0).removeTopic(Caching);
  // RooMsgService::instance().getStream(1).removeTopic(Caching);
  // RooMsgService::instance().getStream(0).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Integration);
  // RooMsgService::instance().setGlobalKillBelow(WARNING);

  // --- set variables ---
  // kinematics
  float massLow = 2.6, massHigh = 3.5;
  float yLow = 0, yHigh = 2.4;
  int cLow = 0, cHigh = 180;
  double ctLow = -1, ctHigh = 4;
  double errLow = 0, errHigh = 0.999;

  // --- make output folders ---
  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);
  gSystem->mkdir(Form("logs%s", userLabel.Data()), kTRUE);

  // --- read input ---
  cout << "\n=== Import inputs ===\n";
  TString ctErrPath = Form("%s/ctau_err_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh);
  cout << ctErrPath.Data() << "\n";
  TFile fInData(ctErrPath.Data());
  RooDataSet *data = (RooDataSet *)fInData.Get("redDataSB");
  data->SetName("data");

  auto getRooVal = [](TFile& f, const char* name)->double {
    if (auto v = (RooRealVar*)f.Get(name)) {
      cout << name << ", " << v->getVal() << "\n";
      return v->getVal();
    } 
    cout << name << ", 0\n";
    return 0.0;
  };

  ctLow  = getRooVal(fInData, "rooCtLow");
  ctHigh = getRooVal(fInData, "rooCtHigh");
  errLow = getRooVal(fInData, "rooCtErrLow");
  errHigh= getRooVal(fInData, "rooCtErrHigh");

  // --- cuts ---
  cout << "\n=== Set kienematic cuts ===\n";
  char reduceDS_woCtErr[3000];

  sprintf(reduceDS_woCtErr,
          "(ctau3D > %.3f && ctau3D < %.3f && ctau3DErr > %.3f && ctau3DErr < %.3f)"
          ,ctLow, ctHigh, errLow, errHigh);

  // reduce the dataset
  cout << "\n=== Reduce datasets ===\n";
  RooDataSet *redData = (RooDataSet *)data->reduce(reduceDS_woCtErr);


  // --- build model ---
  cout << "\n=== Build model ===\n";
  const RooArgSet *row = redData->get();
  RooRealVar *ctau3D = const_cast<RooRealVar *>(static_cast<const RooRealVar *>(row->find("ctau3D")));
  RooRealVar *ctau3DErr = const_cast<RooRealVar *>(static_cast<const RooRealVar *>(row->find("ctau3DErr")));

  // read from err file
  ctau3D->setRange(ctLow, ctHigh);
  ctau3D->setRange("ctauRange", ctLow, ctHigh);
  ctau3DErr->setRange(errLow, errHigh);
  ctau3DErr->setRange("ctauRange", errLow, errHigh);

  // -- bkg punzi term ---
  auto bkgPdf = dynamic_cast<RooHistPdf *>(fInData.Get("bkgPdf")->Clone("bkgPdf"));

  // --- signal ---
  RooRealVar mu("mu", "mean", 0.0, -0.05, 0.05);
  RooRealVar sigma1("sigma1","scale1", 0.01, 0.0001, 50.1);

  RooRealVar sigma12("sigma12", "", 1.1, 1.00000001, 10);
  RooRealVar sigma23("sigma23", "", 1.1, 1.00000001, 10);
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, sigma12));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, sigma23));

  RooGaussModel R1("R1", "R1", *ctau3D, mu, sigma1, *ctau3DErr);
  RooGaussModel R2("R2", "R2", *ctau3D, mu, sigma2, *ctau3DErr);
  RooGaussModel R3("R3", "R3", *ctau3D, mu, sigma3, *ctau3DErr);

  RooRealVar fRes1("fRes1","",0.70,0.0,1.0);
  RooRealVar fRes2("fRes2","",0.25,0.0,1.0);

  RooAddModel ctRes("ctRes", "", RooArgList(R1, R2, R3), RooArgList(fRes1, fRes2));


  // --- deacy model ---
  RooRealVar lambdap1("lambdap1", "lambda + #1", 0.01, 0.0004, 1);
  RooRealVar lambdap2("lambdap2", "lambda + #2", 0.01, 0.0004, 1);
  RooRealVar lambdam1("lambdam1", "lambda - #1", 0.025, 0.0004, 1);
  RooRealVar lambdam2("lambdam2", "lambda - #2", 0.06, 0.0004, 1);
  RooRealVar lambdas1("lambdas1", "lambda 0 #1", 0.20, 0.0004, 1);
  RooRealVar lambdas2("lambdas2", "lambda 0 #2", 0.04, 0.0004, 1);

  RooDecay CtPos1("CtPos1", "right-1", *ctau3D, lambdap1, ctRes, RooDecay::SingleSided);
  RooDecay CtPos2("CtPos2", "right-2", *ctau3D, lambdap2, ctRes, RooDecay::SingleSided);

  RooDecay CtNeg1("CtNeg1", "left-1", *ctau3D, lambdam1, ctRes, RooDecay::Flipped);
  RooDecay CtNeg2("CtNeg2", "left-2", *ctau3D, lambdam2, ctRes, RooDecay::Flipped);

  RooDecay ctMid1("ctMid1", "mid-1", *ctau3D, lambdas1, ctRes, RooDecay::DoubleSided);
  RooDecay ctMid2("ctMid2", "mid-2", *ctau3D, lambdas2, ctRes, RooDecay::DoubleSided);

  RooProdPdf CtPos1_PEE("CtPos1_PEE", "CtPos1 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(CtPos1, RooArgSet(*ctau3D)));

  RooProdPdf CtPos2_PEE("CtPos2_PEE", "CtPos2 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(CtPos2, RooArgSet(*ctau3D)));

  RooProdPdf CtNeg1_PEE("CtNeg1_PEE", "CtNeg1 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(CtNeg1, RooArgSet(*ctau3D)));

  RooProdPdf CtNeg2_PEE("CtNeg2_PEE", "CtNeg2 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(CtNeg2, RooArgSet(*ctau3D)));

  RooProdPdf CtMid1_PEE("CtMid1_PEE", "CtMid1 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(ctMid1, RooArgSet(*ctau3D)));

  RooProdPdf CtMid2_PEE("CtMid2_PEE", "CtMid2 * P(err)",
                        RooArgSet(*bkgPdf),
                        RooFit::Conditional(ctMid2, RooArgSet(*ctau3D)));

  RooProdPdf ctRes_PEE("ctRes_PEE", "ctRes * P(err)",
                       RooArgSet(*bkgPdf),
                       RooFit::Conditional(ctRes, RooArgSet(*ctau3D)));

  RooRealVar fR12("fR12", "frac right 1/(1+2)", 0.1, 0.0, 1.0);
  RooRealVar fL12("fL12", "frac left 1/(1+2)", 0.1, 0.0, 1.0);
  RooRealVar fM12("fM12", "frac mid 1/(1+2)", 0.1, 0.0, 1.0);
  RooRealVar fTotRes("fTotRes", "Resolution weight", 0.1, 0.0, 1.0);
  RooRealVar fTotR("fTotR", "right  weight", 0.1, 0.0, 1.0);
  RooRealVar fTotL("fTotL", "Left  weight", 0.1, 0.0, 1.0);

  // RooAddPdf model("model", "Total Decay (L/R/M with ctRes)",
  //                 {ctRes_PEE,CtPos1_PEE, CtNeg1_PEE, CtMid2_PEE},
  //                 {fR12, fL12, fM12},
  //                 true);
  // RooAddPdf model("model", "Total Decay (L/R/M with ctRes)",
  //               {ctRes_PEE, CtPos1_PEE, CtNeg1_PEE, CtMid1_PEE, CtPos2_PEE, CtNeg2_PEE, CtMid2_PEE},
  //               {fR12, fL12, fM12, fTotRes, fTotR, fTotL}, true);

  RooAddPdf model("model", "Total Decay (L/R/M with ctRes)",
                  {ctRes_PEE, CtPos1_PEE, CtNeg1_PEE, CtMid1_PEE,CtNeg2_PEE},
                  {fR12, fL12, fM12, fTotRes}, true);

  // --- fit ---
  cout << "\n=== Fit ===\n";
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  RooFitResult *fitRes = model.fitTo(*redData, Range("ctauRange"), Save(), PrintLevel(0), PrintEvalErrors(-1), Strategy(2), NumCPU(32), EvalBackend("legacy"), ConditionalObservables(RooArgSet(*ctau3DErr)), RecoverFromUndefinedRegions(2));
  //  

  // === start plot ===
  cout << "\n=== Draw plots ===\n";
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
  RooPlot *fr = ctau3D->frame(Title("")); // , Bins(nBins)
  redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"), ProjWData(RooArgList(*ctau3DErr), *redData), NumCPU(16), Normalization(redData->sumEntries(), RooAbsReal::NumEvent));
  // model.plotOn(fr, Components("ctRes"), Name("ctRes"), LineStyle(kDotted), LineColor(kRed));
  // model.plotOn(fr, Components("ctRight"), Name("ctRight"), LineStyle(kDotted), LineColor(kOrange));
  // model.plotOn(fr, Components("ctLeft"), Name("ctLeft"), LineStyle(kDotted), LineColor(kGreen));
  // model.plotOn(fr, Components("ctMid"), Name("ctMid"), LineStyle(kDotted), LineColor(kYellow));

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
    // if (auto *o = findObj(fr, "ctRes"))
    //   leg.AddEntry(o, "Peak", "pe");
    // if (auto *o = findObj(fr, "ctMid"))
    //   leg.AddEntry(o, "DoulbeSided", "pe");
    // if (auto *o = findObj(fr, "ctRight"))
    //   leg.AddEntry(o, "SingleSided", "pe");
    // if (auto *o = findObj(fr, "ctLeft"))
    //   leg.AddEntry(o, "Flipped", "pe");

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
    
    // fit status
    int st = fitRes->status(); // 0 = success
    if (st != 0)
    {
      tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
      std::ofstream flog(Form("logs%s/ctau_bkg_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
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
    print("f_{Res1}", "fRes1");
    print("f_{Res2}", "fRes2");
  }

  // --- pull pad ---
  c.cd();
  TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2.SetTopMargin(0.00001);
  pad2.SetBottomMargin(0.4);
  pad2.Draw();
  pad2.cd();

  RooHist *hpull = fr->pullHist("data", "model");
  RooPlot *fpull = ctau3D->frame(Title(""));
  fpull->addPlotable(hpull, "P");
  fpull->GetYaxis()->SetTitle("Pull");
  fpull->GetXaxis()->SetTitle("ctau3D");
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
  c.SaveAs(Form("%s/ctau_bkg_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  // === finish plot ===

  // print fit result
  fitRes->Print("V");

  // --- save ---
  TFile f(Form("%s/ctau_bkg_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitRes->Write("fitRes");
  f.Close();

  cout << "\n=== finish ctau_bkg() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}