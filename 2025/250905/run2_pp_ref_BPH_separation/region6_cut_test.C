

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

using namespace RooFit;
using std::cout; using std::string;

// --- next macro: region6_cut_test() ---
// apply cuts
// basic cuts
// region6

void region6_cut_test()
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "NP", region = "LSB";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start region6_cut_test() ===\n";

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
  RooRealVar *ctau3DVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3D"));

  if (!massVar) cout << "Warn: There is no variable 'mass'\n";
  if (!ctau3DVar) cout << "Warn: There is no variable 'ctau3D'\n";

  // // draw mass
  if (massVar)
  {
    double mMin = 2.6, mMax = 3.5;
    RooPlot *f_mass = massVar->frame(Range(mMin, mMax), Title("")); // Bins(60)
    if (hasWeight)
      ds_red->plotOn(f_mass, DataError(RooAbsData::SumW2), WeightVar("weight"));
    else
      ds_red->plotOn(f_mass, DataError(RooAbsData::SumW2));

    TCanvas c_mass("c_mass", "c_mass", 800, 600);
    c_mass.SetLogy();
    f_mass->GetXaxis()->SetTitle("m [GeV/c^2]");
    f_mass->GetYaxis()->SetTitle("Events");
    f_mass->Draw("e");
    c_mass.SaveAs("figs_region6_test/mass_data.png");
  }

  // draw ctau3D
  if (ctau3DVar)
  {
    double ctMin = -0.1, ctMax = 0.8;
    RooPlot *f_ctau = ctau3DVar->frame(Range(ctMin, ctMax), Title("")); // Bins(80)

    if (hasWeight)
      ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2), WeightVar("weight"));
    else
      ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2));

    TCanvas c_ctau("c_ctau", "c_ctau", 800, 600);
    c_ctau.SetLogy();
    f_ctau->GetXaxis()->SetTitle("c#tau_{3D} [mm]");
    f_ctau->GetYaxis()->SetTitle("Events");
    f_ctau->Draw("e");
    c_ctau.SaveAs("figs_region6_test/ctau3D_data.png");
  }

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}