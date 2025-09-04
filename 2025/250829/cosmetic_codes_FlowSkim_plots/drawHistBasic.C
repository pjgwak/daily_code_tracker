#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TStyle.h>
#include <TAxis.h>
#include <iostream>

void drawHistBasic()
{
  // read rootlogon.C
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0); // turn off histogram stat box

  // open file
  TFile *f = TFile::Open("/data/users/pjgwak/work/daily_code_tracker/2025/250827/Pb23_draw_plot_from_FlowSkim/figs/hists.root");
  if (!f)
    std::cerr << "[ERROR] can't open input file\n";

  // get histogram
  // h1
  //  
  // h2
  //  
  //
  TH1D *h = (TH1D *)f->Get("h1/ctau3D_region6_NP_SR_fwd_pt30to50");
  if (!h)
    std::cerr << "[ERROR] can't find histogram\n";

  // draw
  TCanvas *c = new TCanvas("c", "c", 800, 600);

  // cosmetcis - general
  // h->SetLineColor(kBlue + 2);
  // h->SetLineWidth(2);
  // h->SetMarkerStyle(20);
  // h->SetMarkerColor(kBlue + 2);
  h->SetMarkerSize(0.9);

  // cosmetics - specific
  // h->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  h->GetXaxis()->SetTitle("ctau3D (mm)");
  // h->GetYaxis()->SetTitle("Events");
  h->GetXaxis()->CenterTitle();
  // h->GetYaxis()->CenterTitle();

  h->SetStats(0); // double check to remove the stat box
  h->Draw("e");

  c->SaveAs("hist_basic.png");
  // c->SaveAs("hist_basic.pdf");
}