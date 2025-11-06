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

// structure

void mass_fit(double ptLow=6.5, double ptHigh=7.5, double ctMin=0.5, double ctMax=1.0)
{
  cout << "=== start mass_fit() ===\n";
  TStopwatch time;
  time.Start();

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

  RooRealVar mean("mean", "peak mean", 3.096, 3.0, 3.2);
  RooRealVar sigmaL("sigmaL", "sigma left", 0.04, 0.001, 0.1);
  RooRealVar sigmaR("sigmaR", "sigma right", 0.04, 0.001, 0.1);

  RooRealVar alphaL("alphaL", "alpha left", 1.5, 0.2, 5.0);
  RooRealVar nL("nL", "n left", 3.0, 1, 200.0);
  RooRealVar alphaR("alphaR", "alpha right", 2.0, 0.2, 5.0);
  RooRealVar nR("nR", "n right", 3.5, 1, 200.0);

  // double sided CB
  RooCrystalBall dcb("dcb", "Double-sided CrystalBall",
                       *mass, mean, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

  // bkg
  RooRealVar c1("c1", "c1", 0.01, -1.0, 1.0);
  RooRealVar c2("c2", "c2", 0.01, -1.0, 1.0);
  RooRealVar c3("c3", "c3", 0.01, -1.0, 1.0);
  RooChebychev cheb3("cheb3", "Chebyshev order 3", *mass, RooArgList(c1, c2));

  // total model
  RooRealVar nSig("nSig", "signal yield", 1000, 1, 1e6);
  RooRealVar nBkg("nBkg", "background yield", 1000, 1, 1e5);

  // 신호+배경 합
  RooAddPdf model("model", "DCB + Cheby3 (extended)",
                  RooArgList(dcb, cheb3),
                  RooArgList(nSig, nBkg));

  // --- fit ---
  RooFitResult *fitMass = model.fitTo(*redData, Range("massRange"), Save(), Extended(), PrintLevel(0), NumCPU(32), EvalBackend("legacy"));

  // --- draw ---
  TCanvas *c = new TCanvas("c_massfit", "mass fit", 900, 800);
  RooPlot *fMass = mass->frame(Title("Mass fit: DCB + Cheby3"));
  redData->plotOn(fMass, Name("h_data"), DataError(RooAbsData::SumW2));
  model.plotOn(fMass, Name("c_model"), Range("massRange"), LineColor(kBlue), LineWidth(2));
  fMass->Draw();
  c->SaveAs(Form("figs/mass_fit_pT%.1f_%.1f_ct%.2f_%.2f.png", ptLow, ptHigh, ctMin, ctMax));

  fitMass->Print("V");

  // --- save ---
  TFile f(Form("roots/mass_fit_pT%.1f_%.1f_ct%.2f_%.2f.root", ptLow, ptHigh, ctMin, ctMax), "RECREATE");
  fitMass->Write("fitMass");
  f.Close();

  cout << "\n=== finish mass_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}