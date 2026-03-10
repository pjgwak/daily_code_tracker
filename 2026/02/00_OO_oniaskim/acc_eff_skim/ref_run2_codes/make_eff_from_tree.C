#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>

using namespace std;

static inline int ybin_abs_eff(double yAbs)
{
  if (yAbs >= 0.0 && yAbs < 1.6)
    return 0;
  if (yAbs >= 1.6 && yAbs < 2.4)
    return 1;
  return -1;
}

void make_eff_from_tree(const char *inFileName = "roots/efficiency_Prompt_jpsi_ptWgtNoPbPb2023_ExplicitHLT.root",
                        const char *outFileName = "roots/eff_maps_ExplicitHLT.root")
{
  TFile *fin = TFile::Open(inFileName, "READ");
  if (!fin || fin->IsZombie())
  {
    cout << "[ERROR] Cannot open input file: " << inFileName << endl;
    return;
  }

  TTree *tDen = (TTree *)fin->Get("tDen");
  TTree *tNum = (TTree *)fin->Get("tNum");
  if (!tDen || !tNum)
  {
    cout << "[ERROR] tDen or tNum not found in file " << inFileName << endl;
    fin->Close();
    return;
  }

  cout << "tDen entries = " << tDen->GetEntries() << endl;
  cout << "tNum entries = " << tNum->GetEntries() << endl;

  Float_t den_pt, den_y, den_cosHX, den_phiHX;
  Float_t num_pt, num_y, num_cosHX, num_phiHX;

  tDen->SetBranchAddress("pt", &den_pt);
  tDen->SetBranchAddress("y", &den_y);
  tDen->SetBranchAddress("cosHX", &den_cosHX);
  tDen->SetBranchAddress("phiHX", &den_phiHX);

  tNum->SetBranchAddress("pt", &num_pt);
  tNum->SetBranchAddress("y", &num_y);
  tNum->SetBranchAddress("cosHX", &num_cosHX);
  tNum->SetBranchAddress("phiHX", &num_phiHX);

  // ================================
  // === (y, pT) 2D Efficiency (unchanged)
  // ================================
  const int nY = 12;
  double yBins[nY + 1] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
                          0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

  const int nPt = 9;
  double ptBins[nPt + 1] = {0, 3, 6.5, 9, 12, 15, 20, 25, 30, 50};

  TH2D *hDen_y_pt = new TH2D("hDen_y_pt", "Denominator; y; p_{T} (GeV)",
                             nY, yBins, nPt, ptBins);
  TH2D *hNum_y_pt = new TH2D("hNum_y_pt", "Numerator; y; p_{T} (GeV)",
                             nY, yBins, nPt, ptBins);
  TH2D *hAcc_y_pt = (TH2D *)hNum_y_pt->Clone("hAcc_y_pt");
  hAcc_y_pt->SetTitle("Efficiency; y; p_{T} (GeV)");

  Long64_t nDen = tDen->GetEntries();
  for (Long64_t i = 0; i < nDen; ++i)
  {
    tDen->GetEntry(i);
    hDen_y_pt->Fill(den_y, den_pt);
  }

  Long64_t nNum = tNum->GetEntries();
  for (Long64_t i = 0; i < nNum; ++i)
  {
    tNum->GetEntry(i);
    hNum_y_pt->Fill(num_y, num_pt);
  }

  hAcc_y_pt->Divide(hNum_y_pt, hDen_y_pt, 1.0, 1.0, "B");

  // ================================
  // === (pT, |cosHX|) 2D Efficiency split by |y|
  // ================================
  const int nCos = 10;
  const char *ytag[2] = {"y0_1p6", "y1p6_2p4"};

  TH2D *hDen_cos_pt[2], *hNum_cos_pt[2], *hAcc_cos_pt[2];
  for (int iy = 0; iy < 2; ++iy)
  {
    hDen_cos_pt[iy] = new TH2D(Form("hDen_cos_pt_%s", ytag[iy]),
                               Form("Denominator (%s); |cos#theta_{HX}|; p_{T} (GeV)", ytag[iy]),
                               nCos, 0, 1, nPt, ptBins);

    hNum_cos_pt[iy] = new TH2D(Form("hNum_cos_pt_%s", ytag[iy]),
                               Form("Numerator (%s); |cos#theta_{HX}|; p_{T} (GeV)", ytag[iy]),
                               nCos, 0, 1, nPt, ptBins);

    hAcc_cos_pt[iy] = (TH2D *)hNum_cos_pt[iy]->Clone(Form("hAcc_cos_pt_%s", ytag[iy]));
    hAcc_cos_pt[iy]->SetTitle(Form("Efficiency (%s); |cos#theta_{HX}|; p_{T} (GeV)", ytag[iy]));
  }

  for (Long64_t i = 0; i < nDen; ++i)
  {
    tDen->GetEntry(i);
    int iy = ybin_abs_eff(fabs(den_y));
    if (iy < 0)
      continue;
    hDen_cos_pt[iy]->Fill(fabs(den_cosHX), fabs(den_pt));
  }

  for (Long64_t i = 0; i < nNum; ++i)
  {
    tNum->GetEntry(i);
    int iy = ybin_abs_eff(fabs(num_y));
    if (iy < 0)
      continue;
    hNum_cos_pt[iy]->Fill(fabs(num_cosHX), fabs(num_pt));
  }

  for (int iy = 0; iy < 2; ++iy)
  {
    hAcc_cos_pt[iy]->Divide(hNum_cos_pt[iy], hDen_cos_pt[iy], 1.0, 1.0, "B");
  }

  // ================================
  // === (cosHX, phiHX) 2D Efficiency split by |y|
  // ================================
  const int nPhi = 24;

  TH2D *hDen_cos_phi_HX[2], *hNum_cos_phi_HX[2], *hAcc_cos_phi_HX[2];
  for (int iy = 0; iy < 2; ++iy)
  {
    hDen_cos_phi_HX[iy] = new TH2D(Form("hDen_cos_phi_HX_%s", ytag[iy]),
                                   Form("Denominator (%s); cos#theta_{HX}; #phi_{HX} (rad)", ytag[iy]),
                                   nCos, -1., 1., nPhi, -TMath::Pi(), TMath::Pi());

    hNum_cos_phi_HX[iy] = new TH2D(Form("hNum_cos_phi_HX_%s", ytag[iy]),
                                   Form("Numerator (%s); cos#theta_{HX}; #phi_{HX} (rad)", ytag[iy]),
                                   nCos, -1., 1., nPhi, -TMath::Pi(), TMath::Pi());

    hAcc_cos_phi_HX[iy] = (TH2D *)hNum_cos_phi_HX[iy]->Clone(Form("hAcc_cos_phi_HX_%s", ytag[iy]));
    hAcc_cos_phi_HX[iy]->SetTitle(Form("Efficiency (%s); cos#theta_{HX}; #phi_{HX} (rad)", ytag[iy]));
  }

  for (Long64_t i = 0; i < nDen; ++i)
  {
    tDen->GetEntry(i);
    int iy = ybin_abs_eff(fabs(den_y));
    if (iy < 0)
      continue;
    hDen_cos_phi_HX[iy]->Fill(den_cosHX, den_phiHX);
  }

  for (Long64_t i = 0; i < nNum; ++i)
  {
    tNum->GetEntry(i);
    int iy = ybin_abs_eff(fabs(num_y));
    if (iy < 0)
      continue;
    hNum_cos_phi_HX[iy]->Fill(num_cosHX, num_phiHX);
  }

  for (int iy = 0; iy < 2; ++iy)
  {
    hAcc_cos_phi_HX[iy]->Divide(hNum_cos_phi_HX[iy], hDen_cos_phi_HX[iy], 1.0, 1.0, "B");
  }

  // === draw (split panels) + save pdf
  {
    TCanvas *c = new TCanvas("cEff_cos_pt_byY_eff", "Efficiency: |cosHX| vs pT by |y|", 1200, 500);
    c->Divide(2, 1);
    for (int iy = 0; iy < 2; ++iy)
    {
      c->cd(iy + 1);
      gPad->SetRightMargin(0.12);
      hAcc_cos_pt[iy]->Draw("COLZ");
    }
    c->SaveAs("figs/eff_cospt_byY.pdf");
    delete c;
  }

  {
    TCanvas *c = new TCanvas("cEff_cos_phi_byY_eff", "Efficiency: cosHX vs phiHX by |y|", 1200, 500);
    c->Divide(2, 1);
    for (int iy = 0; iy < 2; ++iy)
    {
      c->cd(iy + 1);
      gPad->SetRightMargin(0.12);
      hAcc_cos_phi_HX[iy]->Draw("COLZ");
    }
    c->SaveAs("figs/eff_cosphi_byY.pdf");
    delete c;
  }

  // === save result ===
  TFile *fout = new TFile(outFileName, "RECREATE");
  fout->cd();

  hDen_y_pt->Write();
  hNum_y_pt->Write();
  hAcc_y_pt->Write();

  for (int iy = 0; iy < 2; ++iy)
  {
    hDen_cos_pt[iy]->Write();
    hNum_cos_pt[iy]->Write();
    hAcc_cos_pt[iy]->Write();

    hDen_cos_phi_HX[iy]->Write();
    hNum_cos_phi_HX[iy]->Write();
    hAcc_cos_phi_HX[iy]->Write();
  }

  fout->Close();
  fin->Close();

  cout << "[INFO] Efficiency histograms written to: " << outFileName << endl;
}
