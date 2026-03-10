#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>
#include "make_acc_from_tree.C"
#include "make_eff_from_tree.C"

using namespace std;

// external macros (same directory or in include path)
void make_eff_from_tree(const char *inFileName, const char *outFileName);
void make_acc_from_tree(const char *inFileName, const char *outFileName);

using namespace std;

TH2D *makeAccEff(TH2D *hAcc, TH2D *hEff,
                 const char *name, const char *title)
{
  if (!hAcc || !hEff)
    return nullptr;

  TH2D *h = (TH2D *)hAcc->Clone(name);
  h->SetDirectory(nullptr);
  h->SetTitle(title);
  h->Sumw2();
  h->Multiply(hEff); // bin-by-bin acc * eff (error propagated)
  return h;
}

void make_acc_eff_maps(const char *accFileName = "roots/acc_maps.root",
                       const char *effFileName = "roots/eff_maps_ExplicitHLT.root",
                       const char *outFileName = "roots/acc_eff_maps.root")
{
  // --- open input files ---
  TFile *fAcc = TFile::Open(accFileName, "READ");
  TFile *fEff = TFile::Open(effFileName, "READ");

  if (!fAcc || fAcc->IsZombie() || !fEff || fEff->IsZombie())
  {
    cout << "[ERROR] Cannot open acc/eff input files" << endl;
    return;
  }

  TFile *fout = new TFile(outFileName, "RECREATE");

  // =========================
  // (y, pT)
  // =========================
  TH2D *hAcc_y_pt = (TH2D *)fAcc->Get("hAcc_y_pt");
  TH2D *hEff_y_pt = (TH2D *)fEff->Get("hAcc_y_pt");

  TH2D *hAccEff_y_pt =
      makeAccEff(hAcc_y_pt, hEff_y_pt,
                 "hAccEff_y_pt",
                 "Acceptance#timesEfficiency; y; p_{T} (GeV)");

  // =========================
  // split by |y|
  // =========================
  const char *ytag[2] = {"y0_1p6", "y1p6_2p4"};

  TH2D *hAccEff_cos_pt[2];
  TH2D *hAccEff_cos_phi_HX[2];

  for (int iy = 0; iy < 2; ++iy)
  {

    // (|cosHX|, pT)
    TH2D *hAcc_cos_pt =
        (TH2D *)fAcc->Get(Form("hAcc_cos_pt_%s", ytag[iy]));
    TH2D *hEff_cos_pt =
        (TH2D *)fEff->Get(Form("hAcc_cos_pt_%s", ytag[iy]));

    hAccEff_cos_pt[iy] =
        makeAccEff(hAcc_cos_pt, hEff_cos_pt,
                   Form("hAccEff_cos_pt_%s", ytag[iy]),
                   Form("Acc#timesEff (%s); |cos#theta_{HX}|; p_{T} (GeV)", ytag[iy]));

    // (cosHX, phiHX)
    TH2D *hAcc_cos_phi =
        (TH2D *)fAcc->Get(Form("hAcc_cos_phi_HX_%s", ytag[iy]));
    TH2D *hEff_cos_phi =
        (TH2D *)fEff->Get(Form("hAcc_cos_phi_HX_%s", ytag[iy]));

    hAccEff_cos_phi_HX[iy] =
        makeAccEff(hAcc_cos_phi, hEff_cos_phi,
                   Form("hAccEff_cos_phi_HX_%s", ytag[iy]),
                   Form("Acc#timesEff (%s); cos#theta_{HX}; #phi_{HX}", ytag[iy]));
  }

  // =========================
  // write
  // =========================
  fout->cd();

  if (hAccEff_y_pt)
    hAccEff_y_pt->Write();

  for (int iy = 0; iy < 2; ++iy)
  {
    if (hAccEff_cos_pt[iy])
      hAccEff_cos_pt[iy]->Write();
    if (hAccEff_cos_phi_HX[iy])
      hAccEff_cos_phi_HX[iy]->Write();
  }

  fout->Close();
  fAcc->Close();
  fEff->Close();

  cout << "[INFO] acc*eff histograms written to: "
       << outFileName << endl;
}