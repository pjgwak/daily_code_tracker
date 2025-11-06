#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit ;

void gen2(const char* input, const char* output, const char* suffix) {
   TFile *f = new TFile(input);
   RooWorkspace *ws_in = (RooWorkspace*) f->Get("workspace");
   RooAbsPdf *thepdf = ws_in->pdf(Form("tot_%s",suffix));
   RooDataSet* data_out = (RooDataSet*) ws_in->data(Form("data_%s", suffix))->reduce("2.2<invMass && invMass<4.5");
   RooWorkspace *ws_out = new RooWorkspace("workspace");
   ws_out->import(*data_out);
   ws_out->import(*thepdf);

   // now fit
   RooRealVar *invMass = ws_out->var("invMass");
   invMass->setMin(2.2);
   invMass->setMax(4.5);
   invMass->setRange("MassWindow",2.2,4.5);

   thepdf = ws_out->pdf(Form("tot_%s",suffix));
   thepdf->fitTo(*data_out,Offset(kTRUE), Extended(kTRUE), NumCPU(3), Range("MassWindow"), Save(), Minimizer("Minuit2","Migrad"));

   ws_out->writeToFile(output);
}
