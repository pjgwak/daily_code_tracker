#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooDecay.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooCurve.h"
#include "TLine.h"
#include "TLegend.h"

using namespace RooFit;

void plotCombinedSignals()
{
  // 1. Setup RooFit and define variables
  RooRealVar x("x", "x", 0, 10);
  RooRealVar dt("dt", "per-event error", 0, 1);
  RooRealVar tau3("tau3", "decay time 3", 1.5);
  RooRealVar tau2("tau2", "decay time 2", 2.0);

  // Create a mock dataset for ProjWData
  // In a real application, this would be your actual data
  RooDataSet *data1 = new RooDataSet("data1", "dataset with per-event error", RooArgSet(x, dt));
  for (int i = 0; i < 5000; ++i)
  {
    x.setVal(gRandom->Gaus(5, gRandom->Gaus(0.5, 0.1)));
    dt.setVal(gRandom->Uniform(0.1, 1.0));
    data1->add(RooArgSet(x, dt));
  }

  // Create a RooHistPdf for the per-event variable
  // This is necessary for the ProjWData to work correctly
  TH1 *h_dt = data1->createHistogram("h_dt", dt, Binning(50));
  RooDataHist dataHist_dt("dataHist_dt", "dataHist_dt", dt, h_dt);
  RooHistPdf histPdf_dt("histPdf_dt", "histPdf_dt", RooArgSet(dt), dataHist_dt);

  // 2. Create placeholder models (sig_3dV3 and sig_3dV2)
  // These models include a per-event error, mimicking your setup
  RooDecay decay3("decay3", "decay3", x, tau3, dt, RooDecay::SingleSided);
  RooDecay decay2("decay2", "decay2", x, tau2, dt, RooDecay::SingleSided);

  RooProdPdf sig_3dV3("sig_3dV3", "Signal 3 (Prod)", RooArgSet(decay3, histPdf_dt));
  RooProdPdf sig_3dV2("sig_3dV2", "Signal 2 (Prod)", RooArgSet(decay2, histPdf_dt));

  // 3. Create and Plot Individual Curves (to generate RooCurve objects)
  // We plot on temporary frames to retrieve the RooCurve objects
  RooPlot *temp_frame1 = x.frame();
  sig_3dV3.plotOn(temp_frame1, ProjWData(RooArgSet(dt), *data1));
  RooCurve *curve1 = (RooCurve *)temp_frame1->getObject(0);

  RooPlot *temp_frame2 = x.frame();
  sig_3dV2.plotOn(temp_frame2, ProjWData(RooArgSet(dt), *data1));
  RooCurve *curve2 = (RooCurve *)temp_frame2->getObject(0);

  // 4. Manually Sum the Curves
  RooCurve *sum_curve = new RooCurve(*curve1);
  sum_curve->add(*curve2);

  // 5. Draw Final Plot
  TCanvas *canvas = new TCanvas("canvas", "Combined Signals", 800, 600);
  RooPlot *frame = x.frame(Title("Merged Signal Plots (Manual Sum)"));

  // Plot individual components (for comparison)
  curve1->SetLineColor(kBlue);
  curve1->SetLineStyle(kDashed);
  frame->addObject(curve1);

  curve2->SetLineColor(kRed);
  curve2->SetLineStyle(kDashed);
  frame->addObject(curve2);

  // Plot the manually summed curve
  sum_curve->SetLineColor(kBlack);
  sum_curve->SetLineWidth(1);
  frame->addObject(sum_curve);

  // Draw the data points from a potential total fit (for reference)
  // Assuming a total model and fit, then plotting the dataset
  // For this example, we just plot the mock data we created
  data1->plotOn(frame, MarkerStyle(kPlus), MarkerColor(kGray));

  // Create a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(curve1, "Signal 3dV3 (Component)", "l");
  legend->AddEntry(curve2, "Signal 3dV2 (Component)", "l");
  legend->AddEntry(sum_curve, "Manual Sum", "l");
  legend->AddEntry(frame->getHist(), "Data Points", "ep");

  frame->Draw();
  legend->Draw("same");

  canvas->SaveAs("combined_signals.png");
}
