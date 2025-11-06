#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"

using namespace RooFit;

void roofitPullDistribution() {
    // 1. Define variables and a sharp peak model
    RooRealVar x("x", "x", 0, 10);
    RooRealVar mean("mean", "mean of gaussian", 5, 4, 6);
    RooRealVar sigma("sigma", "width of gaussian", 0.1, 0.01, 0.5);
    RooRealVar c0("c0", "coefficient for linear bg", 0.5, 0, 1);
    
    RooGaussian sig("sig", "signal PDF", x, mean, sigma);
    RooPolynomial bkg("bkg", "background PDF", x, RooArgList(c0));
    
    RooRealVar nsig("nsig", "number of signal events", 1000, 0, 2000);
    RooRealVar nbkg("nbkg", "number of background events", 5000, 0, 10000);
    RooAddPdf model("model", "signal+background", RooArgList(sig, bkg), RooArgList(nsig, nbkg));

    // 2. Generate pseudo-data for the fit
    RooDataSet* data = model.generate(x, 600000);
    
    // 3. Perform the fit
    model.fitTo(*data, Save());
    
    // 4. Create a RooPlot and plot the data and model
    RooPlot* frame = x.frame(Title("RooFit Pull Plot Example"));
    data->plotOn(frame);
    model.plotOn(frame);
    
    // 5. Add the pull distribution to the plot
    RooHist* hpull = frame->pullHist();
    
    // 6. Create the canvas and pads
    TCanvas* c = new TCanvas("rf109_pull", "rf109_pull", 800, 800);
    c->Divide(1, 2, 0, 0, 0);

    // Top pad for the main plot
    c->cd(1);
    gPad->SetPad(0.01, 0.33, 0.99, 0.99);
    gPad->SetBottomMargin(0);
    frame->GetYaxis()->SetTitleOffset(1.6);
    frame->Draw();
    
    // Bottom pad for the pull plot
    c->cd(2);
    gPad->SetPad(0.01, 0.01, 0.99, 0.32);
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3);
    RooPlot* pullFrame = x.frame(Title("Pull Distribution"));
    pullFrame->addPlotable(hpull, "P");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetTitleSize(0.1);
    pullFrame->GetYaxis()->SetLabelSize(0.1);
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetXaxis()->SetTitleSize(0.12);
    pullFrame->GetXaxis()->SetLabelSize(0.1);
    pullFrame->Draw();

    c->SaveAs("peak_test.png");
}
