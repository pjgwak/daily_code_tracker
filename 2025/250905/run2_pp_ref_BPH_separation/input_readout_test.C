#include <iostream> // std::cout
#include <cstdio>   // printf
#include "TStopwatch.h"
#include "TSystem.h" // gSystem
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"

#include "RooGlobalFunc.h" // using namespace RooFit;
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"

using namespace RooFit;
using std::cout;

void input_readout_test()
{
  
  TStopwatch t;
  t.Start();
  cout << "\n=== Start input_readout_test() ===\n";

  // make output folder
  gSystem->mkdir("figs_readout", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!ds) {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in dataset ---\n";
  ds->Print();

  // draw simplt plots
  // check weight
  // const bool hasWeight = (ds->isWeighted() || ds->weightVar() || ds->get()->find("weight"));
  const bool hasWeight = false;

  // variables
  RooRealVar *massVar = dynamic_cast<RooRealVar *>(ds->get()->find("mass"));
  RooRealVar *ctau3DVar = dynamic_cast<RooRealVar *>(ds->get()->find("ctau3D"));

  if (!massVar) cout << "Warn: There is no variable 'mass'\n";
  if (!ctau3DVar) cout << "Warn: There is no variable 'ctau3D'\n";

  // draw mass
  if (massVar)
  {
    double mMin = 2.6, mMax = 3.5;
    RooPlot *f_mass = massVar->frame(Range(mMin, mMax), Bins(60), Title(""));
    if (hasWeight)
      ds->plotOn(f_mass, DataError(RooAbsData::SumW2), WeightVar("weight"));
    else
      ds->plotOn(f_mass, DataError(RooAbsData::SumW2));

    TCanvas c_mass("c_mass", "c_mass", 800, 600);
    c_mass.SetLogy();
    f_mass->GetXaxis()->SetTitle("m [GeV/c^2]");
    f_mass->GetYaxis()->SetTitle("Events");
    f_mass->Draw("e");
    c_mass.SaveAs("figs_readout/mass_data.png");
  }

  // draw ctau3D
  if (ctau3DVar)
  {
    double ctMin = -0.1, ctMax = 0.8;
    RooPlot *f_ctau = ctau3DVar->frame(Range(ctMin, ctMax), Bins(80), Title(""));
    if (hasWeight)
      ds->plotOn(f_ctau, DataError(RooAbsData::SumW2), WeightVar("weight"));
    else
      ds->plotOn(f_ctau, DataError(RooAbsData::SumW2));

    TCanvas c_ctau("c_ctau", "c_ctau", 800, 600);
    c_ctau.SetLogy();
    f_ctau->GetXaxis()->SetTitle("c#tau_{3D} [mm]");
    f_ctau->GetYaxis()->SetTitle("Events");
    f_ctau->Draw("e");
    c_ctau.SaveAs("figs_readout/ctau3D_data.png");
  }

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}