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

using namespace RooFit;

void toy2D_mass_ctau3D_fit()
{
  // ==============================
  // 1. Observables
  // ==============================
  RooRealVar mass("mass", "mass", 2.6, 3.5);
  RooRealVar ctau3D("ctau3D", "ctau3D", -0.5, 4.0);
  RooRealVar ctau3DErr("ctau3DErr", "ctau3DErr", 0.005, 0.08);

  // ==============================
  // 2. Generate ToyMC
  // ==============================
  TRandom3 rnd(12345);
  const int nEvents = 240000;
  RooArgSet obs(mass, ctau3D, ctau3DErr);
  RooDataSet data("data", "ToyMC data", obs);

  for (int i = 0; i < nEvents; i++)
  {
    double sig = (rnd.Rndm() < 0.7);
    double m = sig ? rnd.Gaus(3.1, 0.03) : 2.6 + rnd.Exp(0.25);
    double err = rnd.Gaus(0.02, 0.005);
    if (err < 0.005)
      err = 0.005;
    double ct = rnd.Exp(0.3);
    double smeared_ct = rnd.Gaus(ct, err);

    mass.setVal(m);
    ctau3DErr.setVal(err);
    ctau3D.setVal(smeared_ct);
    data.add(obs);
  }

  // ==============================
  // 3. Output directory setup
  // ==============================
  TString baseDir = "figs";
  TString subDir = Form("%d", nEvents);
  TString outDir = baseDir + "/" + subDir;

  gSystem->mkdir(baseDir, kTRUE);
  gSystem->mkdir(outDir, kTRUE);

  std::cout << "Figures will be saved in: " << outDir << std::endl;

  // ==============================
  // 4. PDFs
  // ==============================
  RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  RooRealVar sigma("sigma", "sigma", 0.03, 0.005, 0.1);
  RooGaussian massSig("massSig", "mass signal", mass, mean, sigma);

  RooRealVar lambda("lambda", "lambda", -1.5, -10.0, 0.0);
  RooExponential massBkg("massBkg", "mass background", mass, lambda);

  RooRealVar fracSig("fracSig", "signal fraction", 0.6, 0.0, 1.0);
  RooAddPdf massPDF("massPDF", "mass total", RooArgList(massSig, massBkg), fracSig);

  RooRealVar tau("tau", "lifetime", 0.3, 0.05, 2.0);
  RooRealVar mean_ct("mean_ct", "mean_ct", 0.0);
  RooFormulaVar sigma_ct("sigma_ct", "@0", RooArgList(ctau3DErr));

  RooGaussModel resModel("resModel", "Resolution Model", ctau3D, mean_ct, sigma_ct);
  RooDecay ctauPDF("ctauPDF", "Lifetime model", ctau3D, tau, resModel, RooDecay::SingleSided);

  RooProdPdf totPDF("totPDF", "mass*ctau model",
                    RooArgSet(massPDF, ctauPDF),
                    Conditional(RooArgSet(ctauPDF), RooArgSet(ctau3DErr)));

  // ==============================
  // 5. Fit
  // ==============================
  RooFitResult *fitres = totPDF.fitTo(data, Save(),
                                      ConditionalObservables(RooArgSet(ctau3DErr)),
                                      NumCPU(8), PrintLevel(-1), Verbose(kFALSE));
  fitres->Print();

  // ==============================
  // 6. 1D Projection + Pull
  // ==============================
  TCanvas *c1D = new TCanvas("c1D", "1D Projections with Pulls", 1600, 800);
  c1D->Divide(2, 2);

  // ------------------------------------
  // Helper: pull 기반 χ² 계산 람다
  // ------------------------------------
  auto chi2_from_pull = [](RooHist *hpull) -> std::pair<double, int>
  {
    double chi2 = 0.0;
    int n = 0;
    for (int i = 0; i < hpull->GetN(); ++i)
    {
      double x, y;
      hpull->GetPoint(i, x, y);
      if (std::isfinite(y))
      {
        chi2 += y * y;
        n++;
      }
    }
    return std::make_pair(chi2, n);
  };

  // ------------------------------------
  // MASS
  // ------------------------------------
  RooPlot *fMass = mass.frame(Title("Mass projection"));
  data.plotOn(fMass, Binning(300), Name("massData"));
  // totPDF.plotOn(fMass, LineColor(kBlue));
  totPDF.plotOn(
      fMass,
      LineColor(kBlue),
      Normalization(data.sumEntries(), RooAbsReal::NumEvent), Name("massPlot"));
  RooHist *hMassPull = fMass->pullHist("massData", "massPlot");

  // pull 기반 chi2 계산
  auto [chi2_m, nbin_m] = chi2_from_pull(hMassPull);
  int npar_mass = fitres->floatParsFinal().getSize(); // 유효 파라미터 자동 카운트
  double ndf_m = std::max(1, nbin_m - npar_mass);
  double chi2ndf_m = chi2_m / ndf_m;

  // plot + chi2 표시
  c1D->cd(1);
  fMass->Draw();
  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.06);
  latex1.SetTextFont(62);
  latex1.SetTextColor(kBlue + 2);
  latex1.SetTextAlign(23);
  latex1.DrawLatex(0.75, 0.88,
                   Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_m, (int)ndf_m, chi2ndf_m));

  // Pull plot
  RooPlot *fMassPull = mass.frame(Title("Mass Pull"));
  fMassPull->addPlotable(hMassPull, "P");
  fMassPull->SetMinimum(-5);
  fMassPull->SetMaximum(5);
  c1D->cd(3);
  fMassPull->Draw();

  // ------------------------------------
  // CTAU3D
  // ------------------------------------
  RooPlot *fCtau = ctau3D.frame(Title("Ctau3D projection"));
  data.plotOn(fCtau, Binning(300));
  // totPDF.plotOn(fCtau, LineColor(kRed), ProjWData(RooArgSet(ctau3DErr), data));
  totPDF.plotOn(
      fCtau,
      LineColor(kRed),
      ProjWData(RooArgSet(ctau3DErr), data),
      Normalization(data.sumEntries(), RooAbsReal::NumEvent));
  RooHist *hCtauPull = fCtau->pullHist();

  // pull 기반 chi2 계산
  auto [chi2_ct, nbin_ct] = chi2_from_pull(hCtauPull);
  int npar_ctau = fitres->floatParsFinal().getSize();
  double ndf_ct = std::max(1, nbin_ct - npar_ctau);
  double chi2ndf_ct = chi2_ct / ndf_ct;

  // plot + chi2 표시
  c1D->cd(2);
  fCtau->Draw();
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.06);
  latex2.SetTextFont(62);
  latex2.SetTextColor(kRed + 1);
  latex2.SetTextAlign(23);
  latex2.DrawLatex(0.75, 0.88,
                   Form("#chi^{2}/NDF = %.1f / %d = %.3f", chi2_ct, (int)ndf_ct, chi2ndf_ct));

  // Pull plot
  RooPlot *fCtauPull = ctau3D.frame(Title("Ctau3D Pull"));
  fCtauPull->addPlotable(hCtauPull, "P");
  fCtauPull->SetMinimum(-5);
  fCtauPull->SetMaximum(5);
  c1D->cd(4);
  fCtauPull->Draw();

  c1D->SaveAs(outDir + "/toy2D_mass_ctau3D_1D_withPull.pdf");

  // ==============================
  // 7. Punzi term (ctau3DErr)
  // ==============================
  TH1 *hErr = data.createHistogram("hErr", ctau3DErr, Binning(50));
  RooDataHist hErrData("hErrData", "binned err", RooArgSet(ctau3DErr), hErr);
  RooHistPdf punziPdf("punziPdf", "Punzi term", RooArgSet(ctau3DErr), hErrData);

  TCanvas *cPunzi = new TCanvas("cPunzi", "ctau3DErr (Punzi term)", 600, 500);
  RooPlot *fErr = ctau3DErr.frame(Title("ctau3DErr (Punzi term)"));
  hErrData.plotOn(fErr);
  punziPdf.plotOn(fErr, LineColor(kGreen + 2));
  fErr->Draw();
  cPunzi->SaveAs(outDir + "/toy2D_mass_ctau3DErr_Punzi.pdf");

  // ==============================
  // 8. 2D Visualization (Surface only)
  // ==============================
  TH1 *h2_pdf = totPDF.createHistogram("h2_pdf", mass, Binning(120), YVar(ctau3D, Binning(120)));
  TH1 *h2_data = data.createHistogram("h2_data", mass, Binning(120), YVar(ctau3D, Binning(120)));

  // Scale the PDF according to the number of data events
  h2_pdf->Scale(h2_data->Integral() / h2_pdf->Integral()); // w/o bin widths
  // double dx = (mass.getMax() - mass.getMin()) / h2_pdf->GetNbinsX();
  // double dy = (ctau3D.getMax() - ctau3D.getMin()) / h2_pdf->GetNbinsY();
  // double scaleFactor = h2_data->Integral() / (h2_pdf->Integral("width") / (dx * dy));
  // h2_pdf->Scale(scaleFactor);

  // ---- chi2/NDF ----
  double chi2 = 0.0;
  int nBins = 0;
  for (int ix = 1; ix <= h2_data->GetNbinsX(); ix++)
  {
    for (int iy = 1; iy <= h2_data->GetNbinsY(); iy++)
    {
      double d = h2_data->GetBinContent(ix, iy);
      double p = h2_pdf->GetBinContent(ix, iy);
      if (p > 0)
      {
        chi2 += (d - p) * (d - p) / p;
        nBins++;
      }
    }
  }
  int nParams = fitres->floatParsFinal().getSize();
  double ndf = nBins - nParams;
  double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

  // ---- 2D surface plot ----
  gStyle->SetOptStat(0);
  TCanvas *c2D = new TCanvas("c2D", "2D Fit Visualization", 800, 600);
  c2D->SetRightMargin(0.15);
  h2_pdf->SetTitle("PDF 2D Surface (lego/surf)");
  h2_pdf->SetLineColor(kRed);
  h2_pdf->SetMaximum(h2_pdf->GetMaximum() * 1.2);
  h2_pdf->Draw("SURF1");

  TLatex latex;
  latex.SetNDC(); // Normalize coordinates to canvas (0–1)
  latex.SetTextSize(0.045);
  latex.SetTextAlign(33); // right-top alignment
  latex.DrawLatex(0.95, 0.92,
                  Form("#chi^{2}/NDF = %.1f / %d = %.4f", chi2, nBins, chi2ndf));

  c2D->SaveAs(outDir + "/toy2D_mass_ctau3D_2Dsurface.png");

  std::cout << "\n=== All figures saved under " << outDir << " ===\n";
}
