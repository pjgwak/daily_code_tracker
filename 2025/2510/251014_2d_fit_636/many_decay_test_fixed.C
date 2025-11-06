/*
 * Corrected version of the many_decay_test example.
 *
 * Changes from the original snippet:
 *  - Use setVal() instead of assignment when updating RooRealVar objects.
 *  - Clamp the generated per‑event error to be within the defined range of dtErr.
 *  - Ensure the random error is positive by taking its absolute value.
 *
 * Note: This code requires ROOT and RooFit to be available in the build
 * environment. To compile, use a command like:
 *   g++ -O2 -o many_decay_test_fixed many_decay_test_fixed.cpp \
 *     `root-config --cflags --libs` -lRooFit -lRooFitCore -lRooStats
 */

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "TCanvas.h"

using namespace RooFit;

void many_decay_test_fixed()
{
  // ======================================================================
  // STEP 1: Define variables and dataset
  // ======================================================================

  // Observables: decay time (dt) and per-event error (dtErr)
  RooRealVar dt("dt", "decay time", -10, 10);
  RooRealVar dtErr("dtErr", "per-event error", 0.01, 1.0);

  // Dataset with both observables
  std::unique_ptr<RooDataSet> data(
      new RooDataSet("data", "data", RooArgSet(dt, dtErr)));

  // Random number generator seed
  RooRandom::randomGenerator()->SetSeed(42);
  const int n_events = 5000;

  for (int i = 0; i < n_events; ++i)
  {
    // Generate a true decay time for a particle with tau=1.5
    double trueDt = RooRandom::randomGenerator()->Exp(1.5);

    // Generate a positive per-event error from a Gaussian distribution.
    // Use absolute value to avoid negative errors and clamp to the dtErr range.
    double error = std::abs(RooRandom::randomGenerator()->Gaus(0.2, 0.05));
    if (error < dtErr.getMin())
      error = dtErr.getMin();
    if (error > dtErr.getMax())
      error = dtErr.getMax();

    // Smear the true decay time by the error and set the values
    dt.setVal(trueDt + RooRandom::randomGenerator()->Gaus(0.0, error));
    dtErr.setVal(error);

    // Add event to the dataset
    data->add(RooArgSet(dt, dtErr));
  }

  // ======================================================================
  // STEP 2: Build the decay models with per-event resolution
  // ======================================================================

  // Common resolution scale factor for all components
  RooRealVar resScale("resScale", "resolution scale factor", 1.0, 0.5, 2.0);

  // --- Component 1: Signal decay ---
  RooRealVar tau1("tau1", "tau1", 1.5, 1.0, 2.0);
  RooDecay decay1_pdf("decay1_pdf", "decay1_pdf", dt, tau1,
                      RooDecay::SingleSided);
  RooProduct sigma1("sigma1", "per-event error scaled",
                    RooArgList(resScale, dtErr));
  RooGaussModel res1_model("res1_model", "resolution model 1", dt,
                           RooConst(0.0), sigma1);
  RooFFTConvPdf conv1_pdf("conv1_pdf", "convolved PDF 1", dt,
                          decay1_pdf, res1_model);

  // --- Component 2: Background decay 1 ---
  RooRealVar tau2("tau2", "tau2", 3.0, 2.0, 4.0);
  RooDecay decay2_pdf("decay2_pdf", "decay2_pdf", dt, tau2,
                      RooDecay::SingleSided);
  RooProduct sigma2("sigma2", "per-event error scaled",
                    RooArgList(resScale, dtErr));
  RooGaussModel res2_model("res2_model", "resolution model 2", dt,
                           RooConst(0.0), sigma2);
  RooFFTConvPdf conv2_pdf("conv2_pdf", "convolved PDF 2", dt,
                          decay2_pdf, res2_model);

  // --- Component 3: Background decay 2 ---
  RooRealVar tau3("tau3", "tau3", 5.0, 4.0, 6.0);
  RooDecay decay3_pdf("decay3_pdf", "decay3_pdf", dt, tau3,
                      RooDecay::SingleSided);
  RooProduct sigma3("sigma3", "per-event error scaled",
                    RooArgList(resScale, dtErr));
  RooGaussModel res3_model("res3_model", "resolution model 3", dt,
                           RooConst(0.0), sigma3);
  RooFFTConvPdf conv3_pdf("conv3_pdf", "convolved PDF 3", dt,
                          decay3_pdf, res3_model);

  // ======================================================================
  // STEP 3: Combine components and fit to data
  // ======================================================================

  // Extended yields for each component
  RooRealVar n_sig("n_sig", "number of signal events", 2000, 0, 5000);
  RooRealVar n_bkg1("n_bkg1", "number of background events 1", 1000, 0, 5000);
  RooRealVar n_bkg2("n_bkg2", "number of background events 2", 2000, 0, 5000);

  // Composite model: sum of convolved PDFs with yields
  RooAddPdf composite_model("composite_model", "3-component decay model",
                            RooArgList(conv1_pdf, conv2_pdf, conv3_pdf),
                            RooArgList(n_sig, n_bkg1, n_bkg2));

  // Fit the model to the dataset; specify dtErr as conditional observable
  composite_model.fitTo(*data, ConditionalObservables(RooArgSet(dtErr)));

  // ======================================================================
  // STEP 4: Plot the model and data
  // ======================================================================

  // Create a canvas and frame for plotting
  TCanvas *c = new TCanvas("c", "Decay Fit Plot", 800, 600);
  std::unique_ptr<RooPlot> frame(dt.frame(Title("Three Decay Models with Per-Event Resolution")));

  // Plot data with error bars
  data->plotOn(frame.get(), DataError(RooAbsData::SumW2), MarkerStyle(20));

  // Plot full model
  composite_model.plotOn(frame.get(), LineColor(kBlue));

  // Plot individual components
  composite_model.plotOn(frame.get(), Components(conv1_pdf), LineStyle(kDashed), LineColor(kRed));
  composite_model.plotOn(frame.get(), Components(conv2_pdf), LineStyle(kDashed), LineColor(kGreen));
  composite_model.plotOn(frame.get(), Components(conv3_pdf), LineStyle(kDashed), LineColor(kOrange));

  // Draw the frame and save the canvas
  frame->Draw();
  c->SaveAs("combined_decay_fit_fixed.png");

  // Clean up: canvas and dataset will be deleted automatically when going out of scope
}