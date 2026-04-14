#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <string>

void plot_muon_distribution(bool isMC = false, bool isPr = false, bool useGlobal = true)
{
  std::cout << "Start plot_muon_distributon()\n";

  // ------------------------------------------------------------------
  // style
  // ------------------------------------------------------------------
  // CMS plot style
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gSystem->mkdir("figs", true);

  // ------------------------------------------------------------------
  // input
  // ------------------------------------------------------------------
  std::string mcLabel = "";
  if (isMC)
    mcLabel = isPr ? "PR" : "NP";
  const std::string globalLabel = useGlobal ? "globalOn" : "globalOff";
  const std::string inputFile =
      Form("skim_roots/skim_OO2025_isMC%d%s_%s_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root",
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), globalLabel.c_str());
  std::cout << "[INFO] input file: " << inputFile << "\n";

  TFile *testFile = TFile::Open(inputFile.c_str(), "READ");
  if (!testFile || testFile->IsZombie())
  {
    std::cout << "[ERROR] failed to open input file: " << inputFile << "\n";
    return;
  }
  testFile->Close();
  delete testFile;

  TChain *myTree = new TChain("myTree");
  myTree->Add(inputFile.c_str());
  if (myTree->GetEntries() <= 0)
  {
    std::cout << "[ERROR] no entries found in: " << inputFile << "\n";
    return;
  }

  // ------------------------------------------------------------------
  // histogram
  // ------------------------------------------------------------------
  TCanvas c1("c1", "", 900, 800);
  gStyle->SetOptStat(0);
  // Single 2D histogram with both muons (mu1 + mu2)
  myTree->Draw("pt1:eta1>>h(40,-2.5,2.5,40,0,10)", "", "goff");
  myTree->Draw("pt2:eta2>>+h", "", "goff");

  // ------------------------------------------------------------------
  // plotting
  // ------------------------------------------------------------------
  c1.cd();
  TH2F *h = (TH2F *)gDirectory->Get("h");
  if (h) {
    h->SetTitle("Soft muon p_{T} vs #eta; #eta; p_{T} (GeV)");
    h->SetMinimum(1.0); // show only bins with content >= 1
    h->Draw("colz");
  }
  else
  {
    std::cout << "[ERROR] histogram h was not created.\n";
    return;
  }
  gPad->SetLogz();
  gPad->SetRightMargin(0.14);
  gPad->SetTopMargin(0.06);

  // ------------------------------------------------------------------
  // labels
  // ------------------------------------------------------------------
  // Latex labels (match CMS Performance plot)
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.045);
    tx.SetTextFont(62); // bold
    tx.SetTextAlign(13);
    tx.DrawLatex(0.15, 0.985, "CMS");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.045);
    tx.SetTextFont(52); // italic
    tx.SetTextAlign(13);
    tx.DrawLatex(0.25, 0.985, "Performance");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.04);
    tx.SetTextFont(42);
    tx.SetTextAlign(33);
    tx.DrawLatex(0.9, 0.985, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.04);
    tx.SetTextFont(42);
    tx.SetTextAlign(13);
    tx.SetTextColor(kRed + 1);
    tx.DrawLatex(0.18, 0.92,
                 Form("Soft Muon (%s, %s)",
                      isMC ? (isPr ? "MC PR" : "MC NP") : "Data",
                      globalLabel.c_str()));
  }

  // ------------------------------------------------------------------
  // acceptance boundary
  // ------------------------------------------------------------------
  // IsAcceptanceQQ pt cut curves (piecewise)
  const double y1 = 3.3;
  const double y2 = 2.1;
  const double y3 = 1.0;
  const double eta_max = 2.4;
  // |eta| <= 1.0
  TLine *h1 = new TLine(0.0, y1, 1.0, y1);
  h1->SetLineStyle(2);
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(4);
  h1->Draw();
  TLine *h1n = new TLine(-1.0, y1, 0.0, y1);
  h1n->SetLineStyle(2);
  h1n->SetLineColor(kBlack);
  h1n->SetLineWidth(4);
  h1n->Draw();
  // 1.0 < |eta| <= 1.3
  TLine *d1 = new TLine(1.0, y1, 1.3, y2);
  d1->SetLineStyle(2);
  d1->SetLineColor(kBlack);
  d1->SetLineWidth(4);
  d1->Draw();
  TLine *d1n = new TLine(-1.3, y2, -1.0, y1);
  d1n->SetLineStyle(2);
  d1n->SetLineColor(kBlack);
  d1n->SetLineWidth(4);
  d1n->Draw();
  // 1.3 < |eta| <= 1.7
  TLine *d2 = new TLine(1.3, y2, 1.7, y3);
  d2->SetLineStyle(2);
  d2->SetLineColor(kBlack);
  d2->SetLineWidth(4);
  d2->Draw();
  TLine *d2n = new TLine(-1.7, y3, -1.3, y2);
  d2n->SetLineStyle(2);
  d2n->SetLineColor(kBlack);
  d2n->SetLineWidth(4);
  d2n->Draw();
  // 1.7 < |eta| <= 2.4
  TLine *h2 = new TLine(1.7, y3, 2.4, y3);
  h2->SetLineStyle(2);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(4);
  h2->Draw();
  TLine *h2n = new TLine(-2.4, y3, -1.7, y3);
  h2n->SetLineStyle(2);
  h2n->SetLineColor(kBlack);
  h2n->SetLineWidth(4);
  h2n->Draw();

  // vertical edges at |eta| = 2.4 for pt >= 1.0
  TLine *vmaxp = new TLine(eta_max, y3, eta_max, 10.0);
  vmaxp->SetLineStyle(2);
  vmaxp->SetLineColor(kBlack);
  vmaxp->SetLineWidth(4);
  vmaxp->Draw();
  TLine *vmaxn = new TLine(-eta_max, y3, -eta_max, 10.0);
  vmaxn->SetLineStyle(2);
  vmaxn->SetLineColor(kBlack);
  vmaxn->SetLineWidth(4);
  vmaxn->Draw();

  // ------------------------------------------------------------------
  // output
  // ------------------------------------------------------------------
  const std::string outputName =
      Form("figs/soft_muon_pT_eta_both_isMC%d%s_%s.png",
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), globalLabel.c_str());
  const std::string outputPdfName =
      Form("figs/soft_muon_pT_eta_both_isMC%d%s_%s.pdf",
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), globalLabel.c_str());
  std::cout << "[INFO] output file: " << outputName << "\n";
  c1.Update();
  c1.SaveAs(outputName.c_str());
  c1.SaveAs(outputPdfName.c_str());
}
