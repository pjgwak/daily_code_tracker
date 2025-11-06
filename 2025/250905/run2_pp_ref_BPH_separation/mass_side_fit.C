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
#include <RooArgList.h> // RooFormulaVar parameter list
#include <RooExponential.h>
#include <RooFitResult.h>

using namespace RooFit;
using std::cout;
using std::string;

// --- next macro: mass_side_fit() ---
// apply cuts
// basic cuts
// region6

void mass_side_fit()
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "NP", region = "LSB";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start mass_side_fit() ===\n";

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
  const TString cutCtauFull = "(ctau3D >= -.10 && ctau3D <= 0.80)";

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
  RooRealVar *massVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("mass"));
  if (!massVar)
    cout << "Warn: There is no variable 'mass'\n";
  massVar->setRange(2.6, 2.95); // Change PDF's range - Need.
  massVar->setRange("fitRange", 2.6, 2.95);

  // === fit by Expo ===
  // --- 1-Expo ---
  RooRealVar tau("tau", "lifetiem", 0.1, 0.01, 10);
  RooFormulaVar c("c", "-1.0/@0", RooArgList(tau)); // c = -1/tau
  RooExponential expo1("expo1", "exp(c*t)", *massVar, c);

  // --- 2-Expo ---
  RooRealVar tau1("tau1", "lifetime1", 0.05, 0.001, 10.0);
  RooRealVar tau2("tau2", "lifetime2", 0.20, 0.001, 10.0);
  RooFormulaVar c1("c1", "-1.0/@0", RooArgList(tau1));
  RooFormulaVar c2("c2", "-1.0/@0", RooArgList(tau2));
  RooExponential e1("e1", "e1", *massVar, c1);
  RooExponential e2("e2", "e2", *massVar, c2);
  RooRealVar f1("f1", "frac fast", 0.5, 0.0, 1.0);
  RooAddPdf expo2("expo2", "f1*e1 + (1-f1)*e2", RooArgList(e1, e2), RooArgList(f1));

  // --- extended fit ---
  RooRealVar N("N", "yield", ds_red->numEntries(), 1, 1e9);
  RooExtendPdf model("model", "extended 1-Expo", expo2, N);

  // --- perform fit ---
  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"));

  // === draw ctau3D ===

  double chi2ndf = 0;

  if (massVar)
  {
    // --- divided canvas ---
    TCanvas c_mass("c_mass", "c_mass", 800, 800);

    // --- main plot ---
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1->SetBottomMargin(0.00001);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    double massMin = 2.6, massMax = 3.5;
    RooPlot *f_ctau = massVar->frame(Range(massMin, massMax), Title("")); // Bins(80)
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
    c_mass.cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();

    RooHist *hpull = f_ctau->pullHist();
    RooPlot *f_pull = massVar->frame(Range(massMin, massMax), Title(""));
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
    double xmin = massMin;
    double xmax = massMax;
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

    c_mass.SaveAs("figs_region6_test/mass_side.png");
  }

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}