#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace RooFit;

void toy2D_mass_ctau3D_binned()
{
  // -------------------------------
  // 1. Variables
  // -------------------------------
  RooRealVar mass("mass", "mass", 2.6, 3.5);
  RooRealVar ctau3D("ctau3D", "ctau3D", -0.5, 4.0);
  RooRealVar ctau3DErr("ctau3DErr", "ctau3DErr", 0.005, 0.08);

  // -------------------------------
  // 2. Generate Toy Data
  // -------------------------------
  TRandom3 rnd(12345);
  const int nEvents = 1000;
  RooArgSet obs(mass, ctau3D, ctau3DErr);
  RooDataSet data("data", "ToyMC", obs);

  for (int i = 0; i < nEvents; i++)
  {
    double sig = (rnd.Rndm() < 0.7);
    double m = sig ? rnd.Gaus(3.1, 0.03) : 2.6 + rnd.Exp(0.25);
    while (m < 2.6 || m > 3.5)
      m = 2.6 + rnd.Exp(0.25);

    double err = rnd.Gaus(0.02, 0.005);
    if (err < 0.005)
      err = 0.005;

    double u = rnd.Rndm();
    double ct = (u < 0.7) ? rnd.Gaus(0.0, 0.008) : rnd.Exp(0.3);
    double smeared_ct = rnd.Gaus(ct, err);

    mass.setVal(m);
    ctau3DErr.setVal(err);
    ctau3D.setVal(smeared_ct);
    data.add(obs);
  }

  // -------------------------------
  // 3. PDFs
  // -------------------------------
  // --- Mass ---
  RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  RooRealVar sigma("sigma", "sigma", 0.03, 0.005, 0.1);
  RooGaussian massSig("massSig", "mass signal", mass, mean, sigma);

  RooRealVar lambda("lambda", "lambda", -1.5, -10.0, 0.0);
  RooExponential massBkg("massBkg", "mass background", mass, lambda);

  RooRealVar fracSig("fracSig", "signal fraction", 0.7, 0.0, 1.0);
  RooAddPdf massPDF("massPDF", "mass total", RooArgList(massSig, massBkg), fracSig);

  // --- ctau3DErr PDF (데이터 분포 기반) ---
  TH1 *hErr = data.createHistogram("hErr", ctau3DErr, Binning(60));
  RooDataHist hErrDH("hErrDH", "binned err", RooArgSet(ctau3DErr), hErr);
  RooHistPdf errPdf("errPdf", "P(err)", RooArgSet(ctau3DErr), hErrDH);

  // --- Resolution model ---
  RooRealVar mean_ct("mean_ct", "mean_ct", 0.0);
  RooGaussModel resModel("resModel", "resolution model", ctau3D, mean_ct, ctau3DErr);

  // --- Prompt & NonPrompt ---
  RooRealVar tauP("tauP", "prompt tau", 0.005, 0.0001, 0.02);
  RooDecay ctPrompt("ctPrompt", "Prompt", ctau3D, tauP, resModel, RooDecay::SingleSided);

  RooRealVar tauNP("tauNP", "nonprompt tau", 0.3, 0.05, 2.0);
  RooDecay ctNP("ctNP", "NonPrompt", ctau3D, tauNP, resModel, RooDecay::SingleSided);

  RooRealVar fPrompt("fPrompt", "prompt fraction", 0.7, 0.0, 1.0);
  RooAddPdf ctCore("ctCore", "ct total", RooArgList(ctPrompt, ctNP), fPrompt);

  // --- errPdf를 곱해서 marginalize ---
  RooProdPdf ctauFull("ctauFull", "ct * P(err)", RooArgSet(ctCore, errPdf));

  // --- 최종 2D 모델 ---
  RooProdPdf totPDF("totPDF", "mass*(ct*P(err))", RooArgSet(massPDF, ctauFull));

  // -------------------------------
  // 4. Fit
  // -------------------------------
  RooFitResult *fitres = totPDF.fitTo(data, Save(), NumCPU(8), PrintLevel(-1), Verbose(kFALSE));

  // -------------------------------
  // 5. Histogram Comparison
  // -------------------------------
  const int NB_MASS = 120;
  const int NB_CTAU = 150;

  TH1D *hMassData = (TH1D *)data.createHistogram("hMassData", mass, Binning(NB_MASS));
  TH1D *hCtauData = (TH1D *)data.createHistogram("hCtauData", ctau3D, Binning(NB_CTAU));
  hMassData->Sumw2();
  hCtauData->Sumw2();

  TH1D *hMassPdf = (TH1D *)totPDF.createHistogram("hMassPdf", mass, Binning(NB_MASS));
  TH1D *hCtauPdf = (TH1D *)totPDF.createHistogram("hCtauPdf", ctau3D, Binning(NB_CTAU));

  hMassPdf->Scale(hMassData->Integral() / hMassPdf->Integral());
  hCtauPdf->Scale(hCtauData->Integral() / hCtauPdf->Integral());

  // -------------------------------
  // 6. Pull & chi2/NDF
  // -------------------------------
  auto make_pull = [](TH1D *hData, TH1D *hPdf, double &chi2, int &ndf)
  {
    TH1D *hPull = (TH1D *)hData->Clone("pull");
    hPull->Reset();
    chi2 = 0;
    ndf = 0;
    for (int i = 1; i <= hData->GetNbinsX(); ++i)
    {
      double d = hData->GetBinContent(i), e = hData->GetBinError(i);
      double p = hPdf->GetBinContent(i);
      if (e <= 0)
        e = 1.0;
      double pull = (d - p) / e;
      hPull->SetBinContent(i, pull);
      hPull->SetBinError(i, 1.0);
      if (std::isfinite(pull))
      {
        chi2 += pull * pull;
        ndf++;
      }
    }
    return hPull;
  };

  double chi2_m, chi2_c;
  int ndf_m, ndf_c;
  TH1D *hMassPull = make_pull(hMassData, hMassPdf, chi2_m, ndf_m);
  TH1D *hCtauPull = make_pull(hCtauData, hCtauPdf, chi2_c, ndf_c);

  double chi2ndf_m = chi2_m / std::max(1, ndf_m - (int)fitres->floatParsFinal().getSize());
  double chi2ndf_c = chi2_c / std::max(1, ndf_c - (int)fitres->floatParsFinal().getSize());

  // -------------------------------
  // 7. Draw
  // -------------------------------
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "Projections + Pull", 1600, 800);
  c1->Divide(2, 2);

  // --- MASS ---
  c1->cd(1);
  hMassData->SetMarkerStyle(20);
  hMassData->Draw("E");
  hMassPdf->SetLineColor(kBlue);
  hMassPdf->SetLineWidth(2);
  hMassPdf->Draw("HISTSAME");

  TLatex t1;
  t1.SetNDC();
  t1.SetTextSize(0.05);
  t1.SetTextColor(kBlue + 2);
  t1.DrawLatex(0.7, 0.9, Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_m, ndf_m, chi2ndf_m));

  c1->cd(3);
  hMassPull->GetYaxis()->SetRangeUser(-5, 5);
  hMassPull->Draw("E");

  // --- CTAU ---
  c1->cd(2);
  hCtauData->SetMarkerStyle(20);
  hCtauData->Draw("E");
  hCtauPdf->SetLineColor(kRed);
  hCtauPdf->SetLineWidth(2);
  hCtauPdf->Draw("HISTSAME");

  TLatex t2;
  t2.SetNDC();
  t2.SetTextSize(0.05);
  t2.SetTextColor(kRed + 1);
  t2.DrawLatex(0.7, 0.9, Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_c, ndf_c, chi2ndf_c));

  c1->cd(4);
  hCtauPull->GetYaxis()->SetRangeUser(-5, 5);
  hCtauPull->Draw("E");

  c1->SaveAs("toy2D_mass_ctau3D_doubleModel_v3.pdf");

  std::cout << "\n=== chi2/NDF (mass) = " << chi2ndf_m
            << ", (ctau3D) = " << chi2ndf_c << " ===" << std::endl;
}
