#include <iostream>
#include <TStopwatch.h>
#include <string.h>
#include <RooDataSet.h>
#include <TFile.h>
#include <RooRealVar.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <RooChebychev.h>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>

using namespace RooFit;

using std::cout;
using std::string;

void s04_mass_fit()
{
  cout << "=== start s04_mass_fit() ===\n";
  TStopwatch time;
  time.Start();

  // set variables
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // observable ranges
  double ctMin = -0.5, ctMax = 2; // lmin, lmax: 1 for lowpT, 2 for higpT
  float errmin = 0.008, errmax = 0.3;

  // read inputs
  string fileNameCErr = "roots/ctau_err.root";
  TFile fInPRMC(fileNameCErr.c_str());
  RooDataSet *dataPRMC = (RooDataSet *)fInPRMC.Get("dsReduced");
  dataPRMC->SetName("dataPRMC");
  const RooArgSet *vars = dataPRMC->get();
  RooRealVar *mass = const_cast<RooRealVar *>(static_cast<const RooRealVar *>(vars->find("mass")));

  // build model
  RooRealVar mean("mean", "mean", 3.096, 3.04, 3.13);
  RooRealVar sigma("sigma", "CB sigma", 0.030, 0.005, 0.100);
  RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.2, 5.0);
  RooRealVar nL("nL", "nL", 3.0, 1.0, 20.0);
  RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2);
  RooRealVar nR("nR", "nR", 3.0, 1.0, 100.0);
  RooRealVar sigmaG("sigmaG", "Gauss sigma", 0.020, 0.002, 0.100);
  RooRealVar fRight("fRight", "frac(CBR in DCB)", 0.50, 0.00, 1.00);
  RooRealVar fG("fG", "frac(Gauss in total)", 0.20, 0.00, 1.00);

  RooCBShape CBL("CBL", "left-tail CB", *mass, mean, sigma, alphaL, nL);
  RooCBShape CBR("CBR", "right-tail CB", *mass, mean, sigma, alphaR, nR);
  RooAddPdf DCB("DCB", "Double-CB approx",
                RooArgList(CBR, CBL),
                RooArgList(fRight), // fRight * CBR + (1-fRight) * CBL
                true);

  RooGaussian G("G", "core Gauss", *mass, mean, sigmaG);
  RooAddPdf Sig("Sig", "Signal DCB+G",
                RooArgList(G, DCB),
                RooArgList(fG), // fG * G + (1-fG) * DCB
                true);

  // --- chebychev ---
  RooRealVar sl1("sl1", "Cheb c1", 0.0, -1.0, 1.0);
  RooRealVar sl2("sl2", "Cheb c2", 0.0, -1.0, 1.0);
  RooChebychev Bkg("Bkg", "Chebyshev(2)", *mass, RooArgList(sl1, sl2));

  RooRealVar nSig("nSig", "signal fraction", 100000, 1.0, 1e8);
  RooRealVar nBkg("nBkg", "signal fraction", 10000, 1.0, 1e8);
  RooAddPdf Tot("Tot", "Total = Sig + Bkg",
                RooArgList(Sig, Bkg), RooArgList(nSig, nBkg));

  // --- use MC result ---
  TFile inputMcFit("roots/mc_mass_fit.root", "READ");
  if (inputMcFit.IsZombie())
  {
    ::Error("fit", "MC fit file cannot be opened: roots/mc_mass_fit.root");
  }
  else
  {
    RooFitResult *fitMcMass = dynamic_cast<RooFitResult *>(inputMcFit.Get("fitMcMass"));
    if (!fitMcMass)
    {
      ::Warning("fit", "RooFitResult 'fitMcMass' not found in roots/mc_mass_fit.root");
    }
    else
    {
      RooArgSet targets;
      // targets.add(alphaL);
      targets.add(nL);
      // targets.add(alphaR);
      targets.add(nR);
      // targets.add(fRight);
      // targets.add(fG);
      // targets.add(sl1);
      // targets.add(sl2);

      // set MC value and fix
      auto apply_from_list = [&](const RooArgList &lst)
      {
        for (int i = 0; i < lst.getSize(); ++i)
        {
          const RooAbsArg *a = lst.at(i);
          const RooRealVar *src = dynamic_cast<const RooRealVar *>(a);
          if (!src)
            continue;
          if (RooAbsArg *destA = targets.find(src->GetName()))
          {
            if (auto *dest = dynamic_cast<RooRealVar *>(destA))
            {
              dest->setVal(src->getVal());
              dest->setConstant(true); //
              // std::cout << "[MC->DATA] " << src->GetName() << " = " << src->getVal() << " (frozen)\n";
            }
          }
        }
      };

      apply_from_list(fitMcMass->floatParsFinal());
      apply_from_list(fitMcMass->constPars());
    }
  }

  // fit
  RooFitResult *fr = Tot.fitTo(*dataPRMC, Extended(), Save(), PrintLevel(-1), PrintEvalErrors(-1), NumCPU(32), EvalBackend("legacy"));

  // draw
  RooPlot *f = mass->frame();
  dataPRMC->plotOn(f, Name("h_data"));
  Tot.plotOn(f, Name("curve_all"), LineColor(kBlue));
  Tot.plotOn(f, Components(DCB), LineStyle(kDashed), LineColor(kRed), Name("curve_dcb"));
  Tot.plotOn(f, Components(G), LineStyle(kDotted), LineColor(kGreen + 2), Name("curve_g"));

  RooHist *pull = f->pullHist("h_data", "curve_all");
  RooPlot *fp = mass->frame();
  fp->addPlotable(pull, "P");
  fp->SetTitle("Pull");
  fp->GetYaxis()->SetTitle("(Data-Model)/#sigma");

  TCanvas c("c", "DCB+Gauss (no WS)", 900, 800);
  c.Divide(1, 2);
  c.cd(1);
  f->Draw();
  {
    TLegend leg(0.62, 0.65, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(f->findObject("h_data"), "Data", "p");
    leg.AddEntry(f->findObject("curve_all"), "Total (DCB+G)", "l");
    leg.AddEntry(f->findObject("curve_dcb"), "DCB", "l");
    leg.AddEntry(f->findObject("curve_g"), "Gauss", "l");
    leg.Draw();
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.035);
    const int nPar = fr ? fr->floatParsFinal().getSize() : 8;
    tx.DrawLatex(0.60, 0.60, Form("#chi^{2}/ndf(vis.) = %.2f", f->chiSquare("curve_all", "h_data", nPar)));
  }
  c.cd(2);
  fp->Draw();
  c.SaveAs("figs/mass_fit.pdf");

  // print out fit result
  fr->Print("V");

  TFile output("roots/mass_fit.root", "recreate");
  fr->Write("fitMass");
  output.Close();

  cout << "=== finish s04_mass_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}