#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

void draw_flowSkim_mass()
{
  // open file
  TFile *f = TFile::Open("/data/users/pjgwak/work/daily_code_tracker/2025/250826/Pb23_flowSkim/FlowSkim_2023PbPbPromptRecoData_132X_miniAOD_noTrack_250826.root");

  // get tree
  TTree *tree = (TTree *)f->Get("myTree");

  // make canvas
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);

  // draw branch
  tree->Draw("mass >> h(100, 2.6, 3.5)", "pt>6.5 && fabs(y)<2.4");

  // double minVal = tree->GetMinimum("mass");
  // double maxVal = tree->GetMaximum("mass");
  // cout << "mass: [" << minVal << ", " << maxVal << "]\n";

  // save results
  c1->SaveAs("mass.png");
}