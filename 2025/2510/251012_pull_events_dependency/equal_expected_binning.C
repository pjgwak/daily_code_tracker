#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TRandom3.h"
#include "RooFitResult.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace RooFit;

void equal_expected_binning()
{
  // 1. Observables
  RooRealVar mass("mass", "mass", 2.6, 3.5);
  RooRealVar ctau3D("ctau3D", "ctau3D", -0.5, 4.0);
  RooRealVar ctau3DErr("ctau3DErr", "ctau3DErr", 0.005, 0.08);

  // 2. Generate ToyMC
  TRandom3 rnd(12345);
  const int nEvents = 100000;
  RooArgSet obs(mass, ctau3D, ctau3DErr);
  RooDataSet data("data", "ToyMC data", obs);

  for (int i = 0; i < nEvents; i++)
  {
    double sig = (rnd.Rndm() < 0.7);
    double m = sig ? rnd.Gaus(3.1, 0.03) : 2.6 + rnd.Exp(0.25);
    while (m < 2.6 || m > 3.5)
      m = 2.6 + rnd.Exp(0.25); // 배경만 재샘플 or 공통 적용

    double err = rnd.Gaus(0.02, 0.005);
    if (err < 0.005)
      err = 0.005;

    // --- Ctau3D distribution (extremely peaked) ---
    double ct;
    double u = rnd.Rndm();
    if (u < 0.65)
    {
      // 65% of events: ultra-narrow Gaussian around 0
      ct = rnd.Gaus(0.0001, 0.01);
    }
    else if (u < 0.9)
    {
      // 25% of events: moderate exponential tail
      ct = rnd.Exp(0.3);
    }
    else
    {
      // 10% of events: long tail (background-like)
      ct = rnd.Exp(1.0);
    }
    double smeared_ct = rnd.Gaus(ct, err);

    mass.setVal(m);
    ctau3DErr.setVal(err);
    ctau3D.setVal(smeared_ct);
    data.add(obs);
  }

  // 3. Output dir
  TString baseDir = "figs";
  TString subDir = Form("%d", nEvents);
  TString outDir = baseDir + "/" + subDir;
  gSystem->mkdir(baseDir, kTRUE);
  gSystem->mkdir(outDir, kTRUE);

  // 4. PDFs
  RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  RooRealVar sigma("sigma", "sigma", 0.03, 0.005, 0.1);
  RooGaussian massSig("massSig", "mass signal", mass, mean, sigma);

  RooRealVar lambda("lambda", "lambda", -1.5, -10.0, 0.0);
  RooExponential massBkg("massBkg", "mass background", mass, lambda);

  RooRealVar fracSig("fracSig", "signal fraction", 0.6, 0.0, 1.0);
  RooAddPdf massPDF("massPDF", "mass total", RooArgList(massSig, massBkg), fracSig);

  // --- Prompt Gaussian core ---
  RooRealVar mean_ct("mean_ct", "mean_ct", 0.0);
  RooRealVar sigma_ct("sigma_ct", "sigma_ct", 0.01, 0.001, 0.05);
  RooGaussian ctPrompt("ctPrompt", "Prompt core", ctau3D, mean_ct, sigma_ct);

  // --- Tail 1: short-lived exponential ---
  RooRealVar tau1("tau1", "tau1", 0.3, 0.05, 2.0);
  RooGaussModel resModel("resModel", "Resolution model", ctau3D, mean_ct, ctau3DErr);
  RooDecay ctTail1("ctTail1", "Short tail", ctau3D, tau1, resModel, RooDecay::SingleSided);

  // --- Tail 2: long-lived exponential ---
  RooRealVar tau2("tau2", "tau2", 1.0, 0.3, 5.0);
  RooDecay ctTail2("ctTail2", "Long tail", ctau3D, tau2, resModel, RooDecay::SingleSided);

  // --- Combine tails ---
  RooRealVar fracShort("fracShort", "fracShort", 0.7, 0.0, 1.0);
  RooAddPdf ctTailSum("ctTailSum", "Tail sum",
                      RooArgList(ctTail1, ctTail2),
                      RooArgList(fracShort));

  // --- Combine with prompt ---
  RooRealVar fracPrompt("fracPrompt", "Prompt fraction", 0.65, 0.0, 1.0);
  RooAddPdf ctauPDF("ctauPDF", "Prompt + 2-tail",
                    RooArgList(ctPrompt, ctTailSum),
                    RooArgList(fracPrompt));
  RooProdPdf totPDF("totPDF", "mass*ctau model",
                    RooArgSet(massPDF, ctauPDF),
                    Conditional(RooArgSet(ctauPDF), RooArgSet(ctau3DErr)));

  // 5. Fit
  RooFitResult *fitres = totPDF.fitTo(data, Save(),
                                      ConditionalObservables(RooArgSet(ctau3DErr)),
                                      NumCPU(8), PrintLevel(-1), Verbose(kFALSE));
  fitres->Print();

  // 6. 1D Projection + Pull, with real-event scaling
  TCanvas *c1D = new TCanvas("c1D", "1D Projections", 1600, 800);
  c1D->Divide(2, 2);

  // helper χ² from pull
  auto chi2_from_pull = [](RooHist *hpull) -> pair<double, int>
  {
    double chi2 = 0;
    int n = 0;
    for (int i = 0; i < hpull->GetN(); i++)
    {
      double x, y;
      hpull->GetPoint(i, x, y);
      if (std::isfinite(y))
      {
        chi2 += y * y;
        n++;
      }
    }
    return {chi2, n};
  };

  // --- MASS ---
  RooPlot *fMass = mass.frame(Title("Mass projection"));
  // fMass->updateNormVars(RooArgSet(mass, ctau3D, ctau3DErr));
  data.plotOn(fMass, Binning(300), Name("massData"), DataError(RooAbsData::SumW2));
  totPDF.plotOn(fMass,
                LineColor(kBlue),
                // ProjWData(RooArgSet(ctau3DErr), data),
                // Normalization(data.sumEntries(), RooAbsReal::NumEvent),
                Name("massModel"));
  // 여기에서 y축 단위 맞추기 위해
  // frame 내부에서 bin width 고려된 스케일이 적용된다는 것을 이용

  RooHist *hMassPull = fMass->pullHist("massData", "massModel");
  auto [chi2_m, nbin_m] = chi2_from_pull(hMassPull);
  int npar = fitres->floatParsFinal().getSize();
  int ndf_m = std::max(1, nbin_m - npar);
  double chi2ndf_m = chi2_m / ndf_m;

  c1D->cd(1);
  fMass->Draw();
  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.06);
  latex1.SetTextFont(62);
  latex1.SetTextColor(kBlue + 2);
  latex1.DrawLatex(0.6, 0.93, Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_m, ndf_m, chi2ndf_m));

  c1D->cd(3);
  RooPlot *fMassPull = mass.frame(Title("Mass Pull"));
  fMassPull->addPlotable(hMassPull, "P");
  fMassPull->SetMinimum(-5);
  fMassPull->SetMaximum(5);
  fMassPull->Draw();

  // --- CTAU3D ---
  RooPlot *fCtau = ctau3D.frame(Title("Ctau3D projection"));
  fCtau->updateNormVars(RooArgSet(mass,ctau3D,ctau3DErr));
  data.plotOn(fCtau, Binning(300), Name("ctauData"), DataError(RooAbsData::SumW2));
  ctauPDF.plotOn(fCtau,
                 LineColor(kRed),
                 ProjWData(RooArgSet(ctau3DErr), data),
                 Normalization(1, RooAbsReal::NumEvent),
                //  Normalization(data.sumEntries(), RooAbsReal::NumEvent),
                 Name("ctauModel"), NumCPU(32));
  RooHist *hCtauPull = fCtau->pullHist("ctauData", "ctauModel");
  auto [chi2_ct, nbin_ct] = chi2_from_pull(hCtauPull);
  int ndf_ct = std::max(1, nbin_ct - npar);
  double chi2ndf_ct = chi2_ct / ndf_ct;

  c1D->cd(2);
  fCtau->Draw();
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.06);
  latex2.SetTextFont(62);
  latex2.SetTextColor(kRed + 1);
  latex2.DrawLatex(0.6, 0.93, Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_ct, ndf_ct, chi2ndf_ct));

  c1D->cd(4);
  RooPlot *fCtauPull = ctau3D.frame(Title("Ctau3D Pull"));
  fCtauPull->addPlotable(hCtauPull, "P");
  fCtauPull->SetMinimum(-5);
  fCtauPull->SetMaximum(5);
  fCtauPull->Draw();

  c1D->SaveAs(outDir + "/toy2D_1D_withRealCounts.pdf");

  std::cout << "\n=== Saved to " << outDir << " ===\n";
}
