#include "TFile.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "TCanvas.h"

void draw_mass_from_roodataset()
{
  // open input
  TFile *f = TFile::Open("/work/pjgwak/pol24/input_roodataset/roots/OniaRooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw0_Accw0_PtW0_TnP0_Run3_PbPb_ptWeightFit.root");
  if (!f || f->IsZombie())
  {
    std::cerr << "Error: cannot open file." << std::endl;
    return 1;
  }

  // load dataset
  RooDataSet *data = (RooDataSet *)f->Get("dataset");
  if (!data)
  {
    std::cerr << "Error: cannot find RooDataSet in file." << std::endl;
    return 1;
  }

  // load mass variable
  RooRealVar *mass = (RooRealVar *)data->get()->find("mass");
  if (!mass)
  {
    std::cerr << "Error: cannot find variable 'mass'." << std::endl;
    return 1;
  }
  mass->setRange(2.6, 3.5);

  // draw mass
  RooPlot *frame = mass->frame();
  data->plotOn(frame);

  TCanvas *c1 = new TCanvas("c1", "Mass distribution", 800, 600);
  frame->Draw();
  c1->SaveAs("mass_distribution.png");


  delete c1;
  f->Close();
  delete f;

  // count
  double nInRegion = data->sumEntries("", "signalRegion");
  std::cout << "Number of events in 3.0 ~ 3.2 GeV region: " << nInRegion << std::endl;
}
