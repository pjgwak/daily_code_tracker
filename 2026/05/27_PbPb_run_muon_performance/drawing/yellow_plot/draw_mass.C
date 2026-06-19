#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>

#include <iostream>
#include <string>
#include <vector>

void draw_mass(const char *inputFile =
                   "/data/users/pjgwak/work/daily_code_tracker/2026/05/"
                   "27_PbPb_run_muon_performance/skimming/skim_roots/PromptReco/"
                   "skim_PbPb2026_isMC0_globalOn_Dimuon_MiniAOD_PromptReco_404337_404474_L1DoubleMu0_MaxDr3p5_Open_yellow.root",
                  //  404337_404358
                  // 404402_404468
               double massLow = 0.2,
               double massHigh = 300,
               int nBins = 1500)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gSystem->mkdir("figs", true);

  TFile *testFile = TFile::Open(inputFile, "READ");
  if (!testFile || testFile->IsZombie()) {
    std::cerr << "[ERROR] failed to open input file: " << inputFile << "\n";
    return;
  }
  testFile->Close();
  delete testFile;

  TChain tree("myTree");
  tree.Add(inputFile);
  if (tree.GetEntries() <= 0) {
    std::cerr << "[ERROR] no entries found in myTree: " << inputFile << "\n";
    return;
  }

  std::string inputName(inputFile);
  std::string baseName = inputName.substr(inputName.find_last_of("/") + 1);
  std::string runLabel = "unknown";
  std::string triggerLabel = "unknown";
  const std::string recoTag = "PromptReco_";
  const size_t runStartTag = baseName.find(recoTag);
  if (runStartTag != std::string::npos) {
    const size_t runStart = runStartTag + recoTag.size();
    const size_t firstRunSep = baseName.find("_", runStart);
    const size_t secondRunSep = firstRunSep == std::string::npos ? std::string::npos : baseName.find("_", firstRunSep + 1);
    if (secondRunSep != std::string::npos) {
      runLabel = baseName.substr(runStart, secondRunSep - runStart);
      std::replace(runLabel.begin(), runLabel.end(), '_', '-');
      const size_t triggerStart = secondRunSep + 1;
      size_t triggerEnd = baseName.rfind("_yellow.root");
      if (triggerEnd == std::string::npos) triggerEnd = baseName.rfind(".root");
      if (triggerEnd != std::string::npos && triggerEnd > triggerStart) {
        triggerLabel = baseName.substr(triggerStart, triggerEnd - triggerStart);
      }
    }
  }

  TCanvas canvas("canvas", "", 900, 800);
  canvas.SetLeftMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.08);
  canvas.SetBottomMargin(0.12);
  canvas.SetTickx(1);
  canvas.SetTicky(1);
  canvas.SetLogy();
  canvas.SetLogx();

  const double xMin = massLow > 0.0 ? massLow : 1e-3;
  std::vector<double> xBins(nBins + 1);
  const double logLow = std::log10(xMin);
  const double logHigh = std::log10(massHigh);
  for (int i = 0; i <= nBins; ++i) {
    xBins[i] = std::pow(10.0, logLow + (logHigh - logLow) * i / nBins);
  }

  TH1F hist("hist", "", nBins, xBins.data());
  hist.Sumw2();
  hist.SetFillColor(kYellow);
  hist.SetFillStyle(1001);
  hist.SetLineColor(kBlack);
  hist.SetLineWidth(2);
  hist.GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
  hist.GetXaxis()->SetTitleFont(132);
  hist.GetXaxis()->SetLabelFont(132);
  hist.GetXaxis()->SetLabelSize(0);
  hist.GetYaxis()->SetTitleFont(132);
  hist.GetYaxis()->SetLabelFont(132);
  hist.GetYaxis()->SetTitle("Events");
  hist.GetXaxis()->CenterTitle();
  hist.GetYaxis()->CenterTitle();
  hist.GetXaxis()->SetTitleOffset(1.1);
  hist.GetYaxis()->SetTitleOffset(1.1);
  hist.GetXaxis()->SetMoreLogLabels(kFALSE);
  hist.GetXaxis()->SetNoExponent(kTRUE);

  const std::string drawExpr = "mass>>hist";
  // const std::string cut = Form("mass > 0.2");
  const std::string cut = Form("mass > 0.2 && mass < %g && pt1 > 4 && pt2 > 4 && abs(eta1) < 2.4 && abs(eta2) < 2.4", massHigh);
  const Long64_t nDraw = tree.Draw(drawExpr.c_str(), cut.c_str(), "goff");
  if (nDraw <= 0 || hist.GetEntries() <= 0) {
    std::cerr << "[ERROR] histogram is empty. Cut: " << cut << "\n";
    return;
  }

  hist.SetMinimum(0.5);
  hist.SetMaximum(hist.GetMaximum() * 15.0);
  hist.Draw("HIST");

  {
    auto xToNdc = [&](double x) {
      return canvas.GetLeftMargin() +
             (std::log10(x) - logLow) / (logHigh - logLow) *
                 (1.0 - canvas.GetLeftMargin() - canvas.GetRightMargin());
    };

    TLatex xLabel;
    xLabel.SetNDC();
    xLabel.SetTextFont(132);
    xLabel.SetTextSize(0.045);
    xLabel.SetTextAlign(23);
    xLabel.DrawLatex(xToNdc(1.0), 0.105, "1");
    xLabel.DrawLatex(xToNdc(3.0), 0.105, "3");
    xLabel.DrawLatex(xToNdc(10.0), 0.105, "10");
    xLabel.DrawLatex(xToNdc(100.0), 0.105, "100");
  }

  {
    auto peakY = [&](double low, double high, double scale) {
      int b1 = hist.FindBin(low);
      int b2 = hist.FindBin(high);
      double y = 1.0;
      for (int b = b1; b <= b2; ++b) {
        y = std::max(y, hist.GetBinContent(b));
      }
      return y * scale;
    };

    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextSize(0.040);
    latex.SetTextAlign(23);
    latex.DrawLatex(0.79, peakY(0.65, 0.90, 1.6), "#rho/#omega");
    latex.DrawLatex(1.02, peakY(0.98, 1.10, 2.1), "#phi");
    latex.DrawLatex(3.20, peakY(2.95, 3.25, 1.6), "J/#psi");
    {
      TLatex psi2S;
      psi2S.SetTextFont(42);
      psi2S.SetTextSize(0.030);
      psi2S.SetTextAlign(21);
      psi2S.DrawLatex(4.35, peakY(3.55, 3.90, 1.35), "#psi(2S)");
    }
    latex.DrawLatex(9.50, peakY(8.5, 10.8, 1.8), "#Upsilon(1,2,3S)");
    latex.DrawLatex(91.2, peakY(75.0, 100.0, 2.2), "Z");
  }

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, Form("PbPb PromptReco run %s", runLabel.c_str()));
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.04);
    tx.SetTextFont(72);
    tx.DrawLatex(0.16, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    const double xtext = 0.16;
    const double y0 = 0.865;
    const double dy = -0.05;
    int k = 0;
    tx.DrawLatex(xtext, y0 + dy * k++, "Data #mu^{+}#mu^{-}");
    tx.DrawLatex(xtext, y0 + dy * k++, triggerLabel.c_str());
    tx.DrawLatex(xtext, y0 + dy * k++, "p_{T}^{#mu} > 4 GeV, |#eta^{#mu}| < 2.4");
  }

  canvas.SaveAs("figs/mass_yellow.pdf");

  TFile outFile("figs/mass_yellow.root", "RECREATE");
  hist.Write("mass");
  outFile.Close();

  std::cout << "[INFO] tree entries: " << tree.GetEntries() << "\n";
  std::cout << "[INFO] selected dimuon candidates: " << nDraw << "\n";
  std::cout << "[INFO] saved figs/mass_yellow.{pdf,root}\n";
}
