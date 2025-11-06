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

void gen1(const char* output, const char* suffix) {
  // S e t u p   m o d e l 
  // ---------------------

   gRandom = new TRandom();
   gRandom->SetSeed(0);
   RooWorkspace *ws = new RooWorkspace("workspace");
  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  ws->factory("invMass[2,5]");
  RooRealVar *invMass = ws->var("invMass");
  RooRealVar *weight = new RooRealVar("weight","weight",1,0,1e6);
  ws->factory(Form("Gaussian::siga_nonorm_%s(invMass,m_%s[3,2.5,3.5],sigmaa_%s[0.02,0.,0.5])",suffix,suffix,suffix));
  ws->factory(Form("Gaussian::sigb_nonorm_%s(invMass,m_%s,sigmab_%s[0.1,0.,0.5])",suffix,suffix,suffix));
  ws->factory(Form("SUM::sig_nonorm_%s(f_%s[0.8,0,1]*siga_nonorm_%s,sigb_nonorm_%s)",suffix,suffix,suffix,suffix));
  ws->factory(Form("RooExtendPdf::sig_%s(sig_nonorm_%s,Nsig_%s[1e4,-1e5,1e5])",suffix,suffix,suffix));
  ws->factory(Form("Gaussian::sig2_nonorm_%s(invMass,m2_%s[3.5,3.,4.],sigma2_%s[0.1,0.,0.5])",suffix,suffix,suffix));
  ws->factory(Form("RooFormulaVar::Nsig2_%s('@0*@1',{Nsig_%s,frac_%s[0.5,-10,10]})",suffix,suffix,suffix));
  ws->factory(Form("RooExtendPdf::sig2_%s(sig2_nonorm_%s,Nsig2_%s)",suffix,suffix,suffix));
  ws->factory(Form("Chebychev::bkg_nonorm_%s(invMass,{lambda0_%s[0.1,-1.5,1.5],lambda1_%s[0.1,-1.5,1.5]})",suffix,suffix,suffix));
  ws->factory(Form("RooExtendPdf::bkg_%s(bkg_nonorm_%s,Nbkg_%s[1e4,-1e5,1e5])",suffix,suffix,suffix));
  // ws->factory(Form("SUM::tot_%s( sig_%s,  bkg_%s)",suffix,suffix,suffix));
  // RooAbsPdf *thepdf = ws->pdf(Form("tot_%s",suffix));
  RooAbsPdf *thepdf = new RooAddPdf(Form("tot_%s",suffix), Form("tot_%s",suffix), 
        RooArgList(*ws->pdf(Form("sig_%s",suffix)), *ws->pdf(Form("sig2_%s",suffix)), *ws->pdf(Form("bkg_%s",suffix))));
  ws->import(*thepdf);

  // G e n e r a t e   e v e n t s 
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  RooDataSet* data = thepdf->generate(*invMass,2e4,Name(Form("data_%s",suffix))) ;  
  ws->import(*data);

  ws->writeToFile(output);
}
