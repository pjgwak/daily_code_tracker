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
#include <TROOT.h>   // gROOT
#include <TSystem.h> // gSystem
#include <TCanvas.h>
#include <TLegend.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TLatex.h>
#include <TAxis.h>
#include <RooGaussModel.h>
#include <TH1D.h>
#include <RooAddModel.h>
#include <RooDecay.h>
#include <fstream>
#include <TLine.h>
#include <RooGaussian.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>

using std::cout;
using std::string;
using namespace RooFit;

static const double CT_EDGES[] = {-0.500, -0.400, -0.300, -0.200, -0.100, 0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900, 2.000}; // widht 0.125 - 25 bins
static const int N_CT_BINS = (int)(sizeof(CT_EDGES) / sizeof(CT_EDGES[0])) - 1;

void ctau_fit_pT6p5_7p5(double ptLow = 17.5, double ptHigh = 20.0)
{
  cout << "=== start ctau_fit_pT6p5_7p5() ===\n";
  TStopwatch time;
  time.Start();

  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // --- set variables ---
  // kinematics
  double ctMin = 2.6, ctMax = 3.5;
  double yLow = 0, yHigh = 2.4;
  int cLow = 0, cHigh = 180;

  TString userLabel = "";
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);

  // --- prepare histogram ---
  TH1D *hNsig = new TH1D("hNsig", "", N_CT_BINS, CT_EDGES);
  hNsig->Sumw2();

  // --- read inputs and fill the hist ---
  for (int i = 0; i < N_CT_BINS; ++i)
  {
    double ctMin = CT_EDGES[i];
    double ctMax = CT_EDGES[i + 1];

    char fname[512];
    std::snprintf(fname, sizeof(fname),
                  "%s/mass_fit_pT%.1f_%.1f_ct%.3f_%.3f.root", rootDir.Data(), ptLow, ptHigh, ctMin, ctMax);

    TFile f(fname, "READ");
    if (f.IsZombie())
    {
      std::cerr << "[Error] Cannot open: " << fname << "\n";
      continue;
    }

    RooFitResult *fr = nullptr;
    fr = dynamic_cast<RooFitResult *>(f.Get("fitMass"));
    if (!fr)
    {
      std::cerr << "[Error] RooFitResult 'fitMass' not found in " << fname << "\n";
      continue;
    }

    const RooArgList &fl = fr->floatParsFinal();
    auto *a = fl.find("nSig");
    auto *v = dynamic_cast<RooRealVar *>(a);
    if (!v)
    {
      std::cerr << "[Error] nSig not found in " << fname << "\n";
      continue;
    }
    
    const double val = v->getVal();
    const double err = v->getError();

    // cout << "\nctau3D [" << ctMin << ", " << ctMax << "]\n";
    // cout << "val: " << val << ", err: " << err << "\n";

    hNsig->SetBinContent(i + 1, val);
    hNsig->SetBinError(i + 1, err);
  }

  // --- check histogram ---
  TCanvas c_hist("c_hist", "", 900, 600);
  hNsig->SetMarkerStyle(20);
  hNsig->SetMarkerSize(0.9);
  hNsig->SetLineWidth(2);
  hNsig->Draw("e");
  c_hist.SaveAs(Form("figs%s/ctau_hist_pT%.1f_%.1f.png",userLabel.Data(), ptLow, ptHigh));

  // --- prepare dataset ---
  ctMin = hNsig->GetXaxis()->GetXmin();
  ctMax = hNsig->GetXaxis()->GetXmax();
  RooRealVar ct("ct", "ct", ctMin, ctMax);

  RooDataHist dh("dh", "", RooArgList(ct), hNsig);

  // --- resolution model model ---
  RooRealVar mean("mean", "mean", 0.0, -0.02, 0.02);
  mean.setVal(0);
  mean.setConstant(true);

  RooRealVar sigma1("sigma1", "core width", 0.02, 0.001, 0.05);
  RooRealVar sigma12("sigma12", "sigma 1 vs 2", 1.5, 1.0, 10);
  RooRealVar sigma23("sigma23", "sigma 2 vs 3", 1.1, 1.0, 5);

  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, sigma12));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, sigma23));

  RooGaussModel g1("g1", "G1 res", ct, mean, sigma1);
  RooGaussModel g2("g2", "G2 res", ct, mean, sigma2);
  RooGaussModel g3("g3", "G2 res", ct, mean, sigma3);

  // recursive fractions
  RooRealVar fG1("fG1", "frac of g1", 0.4, 0.0, 1.0);
  RooRealVar fG2("fG2", "frac of g2", 0.2, 0.0, 1.0);

  RooFormulaVar w1("w1", "@0", RooArgList(fG1));                // = fG1
  RooFormulaVar w2("w2", "(1-@0)*@1", RooArgList(fG1, fG2));     // = (1-fG1)*fG2
  RooFormulaVar w3("w3", "(1-@0)*(1-@1)", RooArgList(fG1, fG2)); // = (1-fG1)*(1-fG2)

  RooAddModel Res("Res", "", RooArgList(g1, g2), RooArgList(fG1));
  // RooAddModel Res("Res", "", RooArgList(g1, g2, g3), RooArgList(w1, w2, w3));

  // RooGaussian Ga1("Ga1", "Gaussian (sigma1)", ct, mean, sigma1);
  // RooGaussian Ga2("Ga2", "Gaussian (sigma2)", ct, mean, sigma2);
  // RooAddPdf GaussModel("GaussModel", "Gauss1+Gauss2 PDF", RooArgList(Ga1, Ga2), RooArgList(fG1));

  // --- decay model ---
  // RooTruthModel TR("TR", "truth", ct);
  RooRealVar lamP1("lamP1", "", 0.02, 0.001, 0.1);
  RooRealVar lamP12("lamP12", "l", 1.1, 0.1, 10);
  RooFormulaVar lamP2("lamP2", "@0*@1", {lamP1, lamP12});

  RooRealVar lamP1N1("lamP1N1", "", 1.1, 0.001, 10);
  RooFormulaVar lamN1("lamN1", "@0*@1", {lamP1N1, lamP1});

  RooDecay decayP1("decayP1", "right exp", ct, lamP1, Res, RooDecay::SingleSided);
  RooDecay decayP2("decayP2", "right exp", ct, lamP2, Res, RooDecay::SingleSided);
  RooDecay decayN("decayN", "left  exp", ct, lamN1, Res, RooDecay::Flipped);

  RooRealVar fP1("fP1", "", 0.1, 0, 1);
  RooRealVar fN1("fN1", "", 0.1, 0, 1);

  // RooFormulaVar wD2("wD2", "(1-@0)*@1", RooArgList(fG1, fG2));     // = (1-fG1)*fG2
  // RooFormulaVar wD3("wD3", "(1-@0)*(1-@1)", RooArgList(fG1, fG2)); // = (1-fG1)*(1-fG2)
  RooAddPdf npModel("npModel", "", {decayN, decayP1, decayP2}, {fN1, fP1}, true);

  // const double Ntot = hNsig->Integral(); // get total number of events in the hist. not used now

  // --- final model ---
  RooRealVar bFrac("bFrac", "The b-fraction", 0.18, 0.10, 0.25);
  RooFormulaVar fPrompt("fPrompt", "1-@0", bFrac);

  RooAddPdf model("model", "", RooArgList(Res, npModel), RooArgList(fPrompt));
  // RooAddPdf model("model", "", RooArgList(Res, decayP1, decayP2, decayN), {bFrac, fP1, fN1}, true);

  // RooRealVar nRes("nRes", "The b-fraction", 100000, 1, 1e8);
  // RooRealVar nP1("nP1", "The b-fraction", 100000, 1, 1e8);
  // RooRealVar nP2("nP2", "The b-fraction", 100000, 1, 1e8);
  // RooRealVar nN("nN", "The b-fraction", 100000, 1, 1e8);
  // RooAddPdf model("model", "", RooArgList(Res, decayP1, decayP2, decayN), {nRes, nP1, nP2, nN});

  // --- fit ---
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  // RooNumIntConfig &cfG1 = *RooAbsReal::defaultIntegratorConfig();
  // cfG1.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  // cfG1.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D")
  //     .setRealValue("epsRel", 1e-5);
  // cfG1.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 10000);

  // RooFitResult *fitCtau = model.chi2FitTo(dh, Save(), PrintLevel(-1), PrintEvalErrors(-1), IntegrateBins(1e-8), Offset(true), Verbose(false), Strategy(2));

  RooFitResult *fitCtau = model.chi2FitTo(dh, Save(true), SumW2Error(true), PrintLevel(-1), PrintEvalErrors(-1), Strategy(2), Optimize(0), IntegrateBins(1e-8), Verbose(false));
  // model.chi2FitTo(dh, Save(true), DataError(RooAbsData::SumW2), PrintLevel(-1));
  //  Extended()

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
  RooPlot *fr = ct.frame(Range(ctMin, ctMax), Title("")); // , Bins(nBins)
  dh.plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(fr, Name("model"));
  model.plotOn(fr, Components("Res"), LineStyle(kDotted), LineColor(kRed), Name("Res"));
  model.plotOn(fr, Components("npModel"), LineStyle(kDotted), LineColor(kOrange), Name("npModel"));

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
      leg.AddEntry(o, "Model", "pe");
    if (auto *o = findObj(fr, "Res"))
      leg.AddEntry(o, "Prompt", "pe");
    if (auto *o = findObj(fr, "npModel"))
      leg.AddEntry(o, "Nonprompt", "pe");

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
    int st = fitCtau->status(); // 0 = success
    if (st != 0)
    {
      tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
      std::ofstream flog(Form("logs%s/ctau_status_pT%.1f_%.1f_ctau%.3f_%.3f.txt", userLabel.Data(), ptLow, ptHigh, ctMin, ctMax), std::ios::app);
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
    print("b-fraction", "bFrac");

    print("mean_{G}", "mean");
    print("#sigma_{G1}", "sigma1");
    print("#sigma_{1/2}", "sigma12");
    print("f_{Res, G1}", "fG1");
    // print("f_{Res, G2}", "fG2");

    print("f_{Decay, P1}", "fP1");
    print("f_{Decay, N}", "fN1");

    print("#lambda_{P1}", "lamP1");
    print("#lambda_{P1/P2}", "lamP12");
    print("#lambda_{P1/N1}", "lamP1N1");
  }

  // --- pull pad ---
  c.cd();
  TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2.SetTopMargin(0.00001);
  pad2.SetBottomMargin(0.4);
  pad2.Draw();
  pad2.cd();

  RooHist *hpull = fr->pullHist("data", "model");
  RooPlot *fpull = ct.frame(Range(ctMin, ctMax), Title(""));
  fpull->addPlotable(hpull, "P");
  fpull->GetYaxis()->SetTitle("Pull");
  fpull->GetXaxis()->SetTitle("c#tau [mm]");
  fpull->GetXaxis()->CenterTitle();
  // fpull->SetMinimum(-8);
  // fpull->SetMaximum(8);
  fpull->GetYaxis()->SetNdivisions(505);
  fpull->GetYaxis()->SetTitleSize(0.12);
  fpull->GetYaxis()->SetLabelSize(0.10);
  fpull->GetXaxis()->SetTitleSize(0.15);
  fpull->GetXaxis()->SetLabelSize(0.10);
  fpull->Draw();

  TLine line(ctMin, 0.0, ctMax, 0.0);
  line.SetLineStyle(2);
  line.Draw("same");

  // --- chi2/ndf ---
  if (fitCtau)
  {
    int npar = fitCtau->floatParsFinal().getSize();
    double chi2ndf = fr->chiSquare("model", "data", npar);
    TLatex tc;
    tc.SetNDC();
    tc.SetTextSize(0.10);
    tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.1f", chi2ndf));
    cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
  }
  c.SaveAs(Form("%s/ctau_fit_pT%.1f_%.1f_ct%.3f_%.3f.png", figDir.Data(), ptLow, ptHigh, ctMin, ctMax));
  // === finish plot ===

  // print fit result
  fitCtau->Print("V");

  // cout << "\nb-fraction: " << (nN.getVal() + nP1.getVal() + nP2.getVal()) / (nN.getVal() + nP1.getVal() + nP2.getVal() + nRes.getVal()) << "\n";

  cout << "\n=== finish ctau_fit_pT6p5_7p5() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}