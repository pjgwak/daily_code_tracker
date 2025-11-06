#ifndef bfrac_h
#define bfrac_h

#include "bin.h"
#include "resultUtils.h"
#include "../../../Efficiency/plotEffs.C"

#include <vector>
#include "TString.h"
#include "RooStats/SamplingDistribution.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooMultiVarGaussian.h"
#include "RooUniform.h"
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

SamplingDistribution combDist(RooFormulaVar formula, RooArgSet pdfs, int nevts=1e4);
double mean(SamplingDistribution* s);
void testCombDist();
RooRealVar Bval(const char* pr_fits, const char* npr_fits, anabin bin);
RooRealVar alphaval(const char* pr_fits, const char* npr_fits, anabin bin, bool is2S=false);

RooRealVar bfrac (const char* pr_fits, // name of the prompt fits directory
      const char* npr_fits,         // name of the non-prompt fits directory
      anabin bin,                  // bin to be considered
      TString obs                  // name of the observable (eg Bfrac_PbPb or N_Psi2S_NPrompt_PP)
      ) {
   bool isPbPb = (obs.Index("_PbPb") != kNPOS);
   bool is2S = (obs.Index("_Psi2S") != kNPOS);
   bool doN = (obs.Index("N_") != kNPOS);
   bool doRfrac = (obs.Index("RFrac2Svs1S_") != kNPOS); 
   bool doNonPrompt = (obs.Index("_NPrompt") != kNPOS);

   RooRealVar bfrac("bfrac","bfrac",0);

   const char* token = isPbPb ? "PbPb" : "PP";
   vector<TString> files_pr = fileList(pr_fits, token);
   vector<TString> files_npr = fileList(npr_fits, token);

   // fix the centrality if needed
   bin.setcentbin(binI(bin.centbin().low(),abs(bin.centbin().high())));
   if (!isPbPb) bin.setcentbin(binI(0,200));

   vector<TString>::const_iterator it;
   TString file_pr, file_npr;
   int nok=0;
   TString thebin_str(Form("pt%i%i_rap%i%i_cent%i%i.root",
            (int) (10*bin.ptbin().low()), (int) (10*bin.ptbin().high()),
            (int) (10*bin.rapbin().low()), (int) (10*bin.rapbin().high()),
            bin.centbin().low(),bin.centbin().high()));
   for (vector<TString>::const_iterator it=files_pr.begin(); it!=files_pr.end(); it++) {
      if (it->Index(thebin_str) != TString::kNPOS) {
         file_pr = *it;
         nok++;
      }
   }
   if (nok != 1) {
      cout << "Error, found " << nok << " files corresponding to the requested bin " << thebin_str << " in " << pr_fits << ". Exiting." << endl;
      return bfrac;
   }
   nok=0;
   for (vector<TString>::const_iterator it=files_npr.begin(); it!=files_npr.end(); it++) {
      if (it->Index(thebin_str) != TString::kNPOS) {
         file_npr = *it;
         nok++;
      }
   }
   if (nok != 1) {
      cout << "Error, found " << nok << " files corresponding to the requested bin " << thebin_str << " in " << npr_fits << ". Exiting." << endl;
      return bfrac;
   }

   // get the different ingredients of the formula
   // first, single ratios
   TFile *fpr = TFile::Open(file_pr); if (!fpr || !fpr->IsOpen()) return bfrac;
   TFile *fnpr = TFile::Open(file_npr); if (!fnpr || !fnpr->IsOpen()) return bfrac;
   RooWorkspace *wspr = (RooWorkspace*) fpr->Get("workspace"); if (!wspr) return bfrac;
   RooWorkspace *wsnpr = (RooWorkspace*) fnpr->Get("workspace"); if (!wsnpr) return bfrac;
   RooRealVar Njpsipr_var(*wspr->var(Form("N_Jpsi_%s", isPbPb ? "PbPb" : "PP"))); Njpsipr_var.SetName("Njpsipr_var");
   RooRealVar Njpsipr("Njpsipr","",Njpsipr_var.getVal());
   RooRealVar Njpsinpr_var(*wsnpr->var(Form("N_Jpsi_%s", isPbPb ? "PbPb" : "PP"))); Njpsinpr_var.SetName("Njpsinpr_var");
   RooRealVar Njpsinpr("Njpsinpr","",Njpsinpr_var.getVal());
   // the number of psi' is a RooFormula, so we have to build it and compute its uncertainty
   RooRealVar rfracpr_var(*wspr->var(Form("RFrac2Svs1S_%s", isPbPb ? "PbPb" : "PP"))); rfracpr_var.SetName("rfracpr_var");
   RooRealVar rfracpr("rfracpr","",rfracpr_var.getVal());
   RooRealVar rfracnpr_var(*wsnpr->var(Form("RFrac2Svs1S_%s", isPbPb ? "PbPb" : "PP"))); rfracnpr_var.SetName("rfracnpr_var");
   RooRealVar rfracnpr("rfracnpr","",rfracnpr_var.getVal());
   RooFitResult *frpr = (RooFitResult*) wspr->obj(Form("fitResult_pdfMASS_Tot_%s", isPbPb ? "PbPb" : "PP"));
   RooFitResult *frnpr = (RooFitResult*) wsnpr->obj(Form("fitResult_pdfMASS_Tot_%s", isPbPb ? "PbPb" : "PP"));
   double corrpr = frpr->correlation(Form("N_Jpsi_%s", isPbPb ? "PbPb" : "PP"), Form("RFrac2Svs1S_%s", isPbPb ? "PbPb" : "PP"));
   double corrnpr = frnpr->correlation(Form("N_Jpsi_%s", isPbPb ? "PbPb" : "PP"), Form("RFrac2Svs1S_%s", isPbPb ? "PbPb" : "PP"));
   double valpr = Njpsipr.getVal() * rfracpr.getVal();
   double valnpr = Njpsinpr.getVal() * rfracnpr.getVal();
   // cf https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
   double errpr = fabs(valpr) * sqrt(pow(Njpsipr_var.getError()/Njpsipr.getVal(),2)
         +pow(rfracpr_var.getError()/rfracpr.getVal(),2)
         +2.*fabs(corrpr*Njpsipr_var.getError()*rfracpr_var.getError()/Njpsipr.getVal()/rfracpr.getVal()));
   double errnpr = fabs(valnpr) * sqrt(pow(Njpsinpr_var.getError()/Njpsinpr.getVal(),2)
         +pow(rfracnpr_var.getError()/rfracnpr.getVal(),2)
         +2.*fabs(corrnpr*Njpsinpr_var.getError()*rfracnpr_var.getError()/Njpsinpr.getVal()/rfracnpr.getVal()));
   RooRealVar Npsi2Spr_var("Npsi2Spr_var","",valpr); Npsi2Spr_var.setError(errpr); RooRealVar Npsi2Spr("Npsi2Spr","",Npsi2Spr_var.getVal());
   RooRealVar Npsi2Snpr_var("Npsi2Snpr_var","",valnpr); Npsi2Snpr_var.setError(errnpr); RooRealVar Npsi2Snpr("Npsi2Snpr","",Npsi2Snpr_var.getVal());

   RooRealVar Njpsipr_err("Njpsipr_err","Njpsipr_err",Njpsipr_var.getError());
   RooRealVar Njpsinpr_err("Nijpsinpr_err","Njpsinpr_err",Njpsinpr_var.getError());
   RooRealVar Npsi2Spr_err("Npsi2Spr_err","Npsi2Spr_err",Npsi2Spr_var.getError());
   RooRealVar Npsi2Snpr_err("Npsi2Snpr_err","Npsi2Snpr_err",Npsi2Snpr_var.getError());
   RooRealVar rfracpr_err("rfracpr_err","rfracpr_err",rfracpr_var.getError());
   RooRealVar rfracnpr_err("Nrfracnpr_err","rfracnpr_err",rfracnpr_var.getError());

   // then efficiencies
   TString feffjpsipr("../Efficiency/files/histos_"); feffjpsipr += "jpsi_"; feffjpsipr += isPbPb ? "pbpb.root" : "pp.root";
   TString feffpsippr("../Efficiency/files/histos_"); feffpsippr += "psi2s_"; feffpsippr += isPbPb ? "pbpb.root" : "pp.root";
   TString feffnpr("../Efficiency/files/histos_npjpsi_"); feffnpr += isPbPb ? "pbpb.root" : "pp.root";
   TString numname("hnumptdepcut_"); numname += (bin.centbin() == binI(0,200)) ? "pt" : "cent"; numname += (bin.rapbin() == binF(0,1.6)) ? "mid" : "fwd";
   TString denname("hnum_"); denname += (bin.centbin() == binI(0,200)) ? "pt" : "cent"; denname += (bin.rapbin() == binF(0,1.6)) ? "mid" : "fwd";
   TFile *tfeffjpsipr = TFile::Open(feffjpsipr);
   TFile *tfeffpsippr = TFile::Open(feffpsippr);
   TFile *tfeffnpr = TFile::Open(feffnpr);

   TH1F *hnumjpsipr = (TH1F*) tfeffjpsipr->Get(numname);
   TH1F *hnumpsippr = (TH1F*) tfeffpsippr->Get(numname);
   TH1F *hnumnpr = (TH1F*) tfeffnpr->Get(numname);
   TH1F *hdenjpsipr = (TH1F*) tfeffjpsipr->Get(denname);
   TH1F *hdenpsippr = (TH1F*) tfeffpsippr->Get(denname);
   TH1F *hdennpr = (TH1F*) tfeffnpr->Get(denname);

   RooRealVar numjpsipr("numjpsipr","",0), numjpsipr_err("numjpsipr_err","",0), numnjpsipr("numnjpsipr","",0), numnjpsipr_err("numnjpsipr_err","",0);
   RooRealVar numpsippr("numpsippr","",0), numpsippr_err("numpsippr_err","",0), numnpsippr("numnpsippr","",0), numnpsippr_err("numnpsippr_err","",0);
   RooRealVar numnpr("numnpr","",0), numnpr_err("numnpr_err","",0), numnnpr("numnnpr","",0), numnnpr_err("numnnpr_err","",0);
   RooRealVar failjpsipr("failjpsipr","",0), failjpsipr_err("failjpsipr_err","",0), failnjpsipr("failnjpsipr","",0), failnjpsipr_err("failnjpsipr_err","",0);
   RooRealVar failpsippr("failpsippr","",0), failpsippr_err("failpsippr_err","",0), failnpsippr("failnpsippr","",0), failnpsippr_err("failnpsippr_err","",0);
   RooRealVar failnpr("failnpr","",0), failnpr_err("failnpr_err","",0), failnnpr("failnnpr","",0), failnnpr_err("failnnpr_err","",0);
   if (bin == anabin(0,1.6,6.5,30,0,200) || bin == anabin(1.6,2.4,3,30,0,200)) { // Min.Bias bin
      hnumjpsipr = integrateHist(hnumjpsipr); numjpsipr.setVal(hnumjpsipr->GetBinContent(1)); numjpsipr_err.setVal(hnumjpsipr->GetBinError(1));
      hnumpsippr = integrateHist(hnumpsippr); numpsippr.setVal(hnumpsippr->GetBinContent(1)); numpsippr_err.setVal(hnumpsippr->GetBinError(1));
      hnumnpr = integrateHist(hnumnpr); numnpr.setVal(hnumnpr->GetBinContent(1)); numnpr_err.setVal(hnumnpr->GetBinError(1));
      hdenjpsipr = integrateHist(hdenjpsipr); failjpsipr.setVal(hdenjpsipr->GetBinContent(1)-numjpsipr.getVal()); failjpsipr_err.setVal(sqrt(pow(hdenjpsipr->GetBinError(1),2)-pow(numjpsipr_err.getVal(),2)));
      hdenpsippr = integrateHist(hdenpsippr); failpsippr.setVal(hdenpsippr->GetBinContent(1)-numpsippr.getVal()); failpsippr_err.setVal(sqrt(pow(hdenpsippr->GetBinError(1),2)-pow(numpsippr_err.getVal(),2)));
      hdennpr = integrateHist(hdennpr); failnpr.setVal(hdennpr->GetBinContent(1)-numnpr.getVal()); failnpr_err.setVal(sqrt(pow(hdennpr->GetBinError(1),2)-pow(numnpr_err.getVal(),2)));
   } else {
      for (int i=1; i<=hnumjpsipr->GetNbinsX(); i++) { // loop over the histogram bins to find the one we want
         if ((bin.centbin()==binI(0,200) && bin.ptbin()==binF(hnumjpsipr->GetXaxis()->GetBinLowEdge(i),hnumjpsipr->GetXaxis()->GetBinUpEdge(i))) // this is the bin we're looking for (pt)
               || bin.centbin()==binI(2*hnumjpsipr->GetXaxis()->GetBinLowEdge(i),2*hnumjpsipr->GetXaxis()->GetBinUpEdge(i))) { // this is the bin we're looking for (centrality)
            numjpsipr.setVal(hnumjpsipr->GetBinContent(i));
            numpsippr.setVal(hnumpsippr->GetBinContent(i));
            numjpsipr_err.setVal(hnumjpsipr->GetBinError(i));
            numpsippr_err.setVal(hnumpsippr->GetBinError(i));
            numnpr.setVal(hnumnpr->GetBinContent(i));
            numnpr_err.setVal(hnumnpr->GetBinError(i));
            failjpsipr.setVal(hdenjpsipr->GetBinContent(i)-numjpsipr.getVal());
            failpsippr.setVal(hdenpsippr->GetBinContent(i)-numpsippr.getVal());
            failnpr.setVal(hdennpr->GetBinContent(i)-numnpr.getVal());
            failjpsipr_err.setVal(sqrt(pow(hdenjpsipr->GetBinError(i),2)-pow(numjpsipr_err.getVal(),2)));
            failpsippr_err.setVal(sqrt(pow(hdenpsippr->GetBinError(i),2)-pow(numpsippr_err.getVal(),2)));
            failnpr_err.setVal(sqrt(pow(hdennpr->GetBinError(i),2)-pow(numnpr_err.getVal(),2)));
            break;
         }
      }
   }

   delete fpr; fpr=0;
   delete fnpr; fnpr=0;
   delete tfeffjpsipr; tfeffjpsipr=0;
   delete tfeffpsippr; tfeffpsippr=0;
   delete tfeffnpr; tfeffnpr=0;

   // build the RooRealVar for each ingredient
   RooRealVar numjpsipr_var("numjpsipr_var","numjpsipr_var",numjpsipr.getVal(),0,numjpsipr.getVal()+100.*numjpsipr_err.getVal()); numjpsipr_var.setError(numjpsipr_err.getVal());
   RooRealVar numpsippr_var("numpsippr_var","numpsippr_var",numpsippr.getVal(),0,numpsippr.getVal()+100.*numpsippr_err.getVal()); numpsippr_var.setError(numpsippr_err.getVal());
   RooRealVar numnpr_var("numnpr_var","numnpr_var",numnpr.getVal(),0,numnpr.getVal()+100.*numnpr_err.getVal()); numnpr_var.setError(numnpr_err.getVal());
   RooRealVar failjpsipr_var("failjpsipr_var","failjpsipr_var",failjpsipr.getVal(),0,failjpsipr.getVal()+100.*failjpsipr_err.getVal()); failjpsipr_var.setError(failjpsipr_err.getVal());
   RooRealVar failpsippr_var("failpsippr_var","failpsippr_var",failpsippr.getVal(),0,failpsippr.getVal()+100.*failpsippr_err.getVal()); failpsippr_var.setError(failpsippr_err.getVal());
   RooRealVar failnpr_var("failnpr_var","failnpr_var",failnpr.getVal(),0,failnpr.getVal()+100.*failnpr_err.getVal()); failnpr_var.setError(failnpr_err.getVal());
   RooRealVar *Npr_var, *Npr, *Npr_err, *Nnpr_var, *Nnpr, *Nnpr_err;
   if (!is2S) {
      Npr_var = &Njpsipr_var; Npr = &Njpsipr; Npr_err = &Njpsipr_err;
      Nnpr_var = &Njpsinpr_var; Nnpr = &Njpsinpr; Nnpr_err = &Njpsinpr_err;
   } else {
      Npr_var = &Npsi2Spr_var; Npr = &Npsi2Spr; Npr_err = &Npsi2Spr_err;
      Nnpr_var = &Npsi2Snpr_var; Nnpr = &Npsi2Snpr; Nnpr_err = &Npsi2Snpr_err;
   }

   // build the PDFs: Gaussian for numbers of events, Gaussian (PbPb) or Poissonian (pp) for MC numbers of events
   RooGaussian Npr_pdf("Npr_pdf","",*Npr_var,*Npr,*Npr_err);
   RooGaussian Nnpr_pdf("Nnpr_pdf","",*Nnpr_var,*Nnpr,*Nnpr_err);
   double elempr[4] = {pow(Njpsipr_var.getError(),2),corrpr*Njpsipr_var.getError()*rfracpr_var.getError(),corrpr*Njpsipr_var.getError()*rfracpr_var.getError(),pow(rfracpr_var.getError(),2)};
   double elemnpr[4] = {pow(Njpsinpr_var.getError(),2),corrnpr*Njpsinpr_var.getError()*rfracnpr_var.getError(),corrnpr*Njpsinpr_var.getError()*rfracnpr_var.getError(),pow(rfracnpr_var.getError(),2)};
   RooMultiVarGaussian Nmultipr_pdf("Nmultipr_pdf","",RooArgList(Njpsipr_var,rfracpr_var),RooArgList(Njpsipr,rfracpr),TMatrixDSym(2,elempr));
   RooMultiVarGaussian Nmultinpr_pdf("Nmultinpr_pdf","",RooArgList(Njpsinpr_var,rfracnpr_var),RooArgList(Njpsinpr,rfracnpr),TMatrixDSym(2,elemnpr));
   RooAbsPdf *numjpsipr_pdf=NULL; 
   numjpsipr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("numjpsipr_pdf","",numjpsipr_var,numjpsipr,numjpsipr_err) : 
      (RooAbsPdf*) new RooPoisson("numjpsipr_pdf","",numjpsipr_var,numjpsipr);
   RooAbsPdf *numpsippr_pdf=NULL; 
   numpsippr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("numpsippr_pdf","",numpsippr_var,numpsippr,numpsippr_err) : 
      (RooAbsPdf*) new RooPoisson("numpsippr_pdf","",numpsippr_var,numpsippr);
   RooAbsPdf *numnpr_pdf=NULL; 
   numnpr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("numnpr_pdf","",numnpr_var,numnpr,numnpr_err) : 
      (RooAbsPdf*) new RooPoisson("numnpr_pdf","",numnpr_var,numnpr);
   RooAbsPdf *failjpsipr_pdf=NULL; 
   failjpsipr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("failjpsipr_pdf","",failjpsipr_var,failjpsipr,failjpsipr_err) : 
      (RooAbsPdf*) new RooPoisson("failjpsipr_pdf","",failjpsipr_var,failjpsipr);
   RooAbsPdf *failpsippr_pdf=NULL; 
   failpsippr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("failpsippr_pdf","",failpsippr_var,failpsippr,failpsippr_err) : 
      (RooAbsPdf*) new RooPoisson("failpsippr_pdf","",failpsippr_var,failpsippr);
   RooAbsPdf *failnpr_pdf=NULL; 
   failnpr_pdf = isPbPb ? 
      (RooAbsPdf*) new RooGaussian("failnpr_pdf","",failnpr_var,failnpr,failnpr_err) : 
      (RooAbsPdf*) new RooPoisson("failnpr_pdf","",failnpr_var,failnpr);

   // variable list
   RooArgList *lvar = new RooArgList(*Npr_var,*Nnpr_var,numjpsipr_var,numnpr_var,failjpsipr_var,failnpr_var);
   if (is2S) lvar = new RooArgList(*Npr_var,*Nnpr_var,numpsippr_var,numnpr_var,failpsippr_var,failnpr_var);
   if (doRfrac) lvar->add(RooArgList(rfracpr_var,rfracnpr_var,numpsippr_var,failpsippr_var));

   // pdf list
   RooArgSet *pdfs = new RooArgSet(Npr_pdf,Nnpr_pdf,*numjpsipr_pdf,*numnpr_pdf,*failjpsipr_pdf,*failnpr_pdf);
   if (is2S) pdfs = new RooArgSet(Npr_pdf,Nnpr_pdf,*numpsippr_pdf,*numnpr_pdf,*failpsippr_pdf,*failnpr_pdf);
   if (doRfrac) pdfs = new RooArgSet(Nmultipr_pdf,Nmultinpr_pdf,*numjpsipr_pdf,*numpsippr_pdf,*numnpr_pdf,*failjpsipr_pdf,*failpsippr_pdf,*failnpr_pdf);

   // formula
   TString formulastr;
   if (!doN) formulastr = "((@2/(@2+@4)) - (@0/(@0+@1))) / ((@2/(@2+@4)) - (@3/(@3+@5)))"; // B-fraction
   else if (doNonPrompt) formulastr = "(@0+@1) * (((@2/(@2+@4)) - (@0/(@0+@1))) / ((@2/(@2+@4)) - (@3/(@3+@5))))"; // Nb of non-prompt onia
   else formulastr = "(@0+@1) * (1-(((@2/(@2+@4)) - (@0/(@0+@1))) / ((@2/(@2+@4)) - (@3/(@3+@5)))))"; // Nb of prompt onia
   if (doRfrac) {
      if (doNonPrompt) formulastr = "((@0*@6+@1*@7)/(@0+@1))" // (Njpsi,prompt*Rprompt+Njspi,nonprompt*Rnonprompt) / (Njpsi,prompt+Njpsi,nonprompt)
         " * (((@8/(@8+@9)) - (@0*@6/(@0*@6+@1*@7))) / ((@8/(@8+@9)) - (@3/(@3+@5))))" // B fraction for psi prime
         " / (((@2/(@2+@4)) - (@0/(@0+@1))) / ((@2/(@2+@4)) - (@3/(@3+@5))))"; // B fraction for jpsi
      else formulastr = "((@0*@6+@1*@7)/(@0+@1))" // (Njpsi,prompt*Rprompt+Njspi,nonprompt*Rnonprompt) / (Njpsi,prompt+Njpsi,nonprompt)
         " * (1-(((@8/(@8+@9)) - (@0*@6/(@0*@6+@1*@7))) / ((@8/(@8+@9)) - (@3/(@3+@5)))))" // 1- (B fraction for psi prime)
         " / (1-(((@2/(@2+@4)) - (@0/(@0+@1))) / ((@2/(@2+@4)) - (@3/(@3+@5)))))"; // 1 - (B fraction for jpsi)
   }
   RooFormulaVar formula("formula","formula",formulastr.Data(),*lvar);

   // compute the final variable
   double bfracval = formula.getVal();
   // for (int i=0; i<lvar->getSize(); i++) cout << ((RooRealVar*) &((*lvar)[i]))->getVal() << " +- " << ((RooRealVar*) &((*lvar)[i]))->getError() << ", ";
   // cout << endl;
   // cout << formulastr.Data() << endl;
   // cout << bfracval << endl;
   SamplingDistribution s = combDist(formula, *pdfs, 1e4);

   // bfrac.setVal(s.InverseCDF(0.5)); // median
   // bfrac.setVal(mean(&s)); // mean
   bfrac.setVal(bfracval); // mode?
   bfrac.setAsymError(s.InverseCDF(ROOT::Math::normal_cdf(-1))-bfrac.getVal(),s.InverseCDF(ROOT::Math::normal_cdf(1))-bfrac.getVal());
   bfrac.setError((fabs(bfrac.getErrorLo())+fabs(bfrac.getErrorHi()))/2.);
   // cout << bfracval << " " << s.InverseCDF(0.5) << " " << mean(&s) << ", " << s.InverseCDF(ROOT::Math::normal_cdf(-1)) << " " << s.InverseCDF(ROOT::Math::normal_cdf(1)) << endl;

   delete numjpsipr_pdf; numjpsipr_pdf=0;
   delete numpsippr_pdf; numpsippr_pdf=0;
   delete numnpr_pdf; numnpr_pdf=0;
   delete failjpsipr_pdf; failjpsipr_pdf=0;
   delete failpsippr_pdf; failpsippr_pdf=0;
   delete failnpr_pdf; failnpr_pdf=0;

   return bfrac;
}

SamplingDistribution combDist(RooFormulaVar formula, RooArgSet pdfs, int nevts) {
   RooRealVar formulaRRV("formulaRRV","formulaRRV",-1e99,1e99);
   RooDataSet data("sampledist","sampledist",RooArgSet(formulaRRV));
   double themin=1e99, themax=-1e99;
   for (int i=0; i<nevts; i++) {
      TIterator* it = pdfs.createIterator();
      RooAbsPdf *thePdf = (RooAbsPdf*) it->Next();
      while (thePdf) {
         thePdf->generateEvent((TString(thePdf->ClassName())=="RooMultiVarGaussian") ? -1 : 1); // generateEvent(-1) for MultiVarGaussian, (1) for others
         thePdf = (RooAbsPdf*) it->Next();
      }
      // cout << formula->getVal() << endl;
      double theval = formula.getVal();
      if (theval<themin) themin = theval;
      if (theval>themax) themax = theval;
      formulaRRV.setVal(theval);
      data.add(RooArgSet(formulaRRV));
   }

   // // plot the distribution
   // RooPlot *frame = formulaRRV.frame(Bins(100),Range(themin, themax));
   // // RooPlot *frame = formulaRRV.frame(Bins(100),Range(-0.5, 1.5));
   // TCanvas *c1 = new TCanvas();
   // data.plotOn(frame);
   // frame->Draw();
   // c1->SaveAs("test.pdf");
   // delete c1; c1=0;
   // delete frame; frame=0;

   return SamplingDistribution("s","s",data);
}

void testCombDist() {
   RooRealVar x("x","x",-5,5);
   RooRealVar y("y","y",-5,5);
   RooFormulaVar f("f","f","@0+@1",RooArgList(x,y));
   RooUniform ux("ux","ux",RooArgSet(x));
   RooUniform uy("uy","uy",RooArgSet(y));
   RooArgSet s(ux,uy);
   combDist(f,s,100000);
}

double mean(SamplingDistribution *s) {
   vector<double> v = s->GetSamplingDistribution();
   if (v.size()==0) return 0;
   double ans=0;
   for (vector<double>::const_iterator it=v.begin(); it!=v.end(); it++) ans+=*it;

   return ans / v.size();
}

RooRealVar Bval(const char* pr_fits, const char* npr_fits, anabin bin) {
   RooRealVar alpha = alphaval(pr_fits, npr_fits, bin);
   RooRealVar Nnp_PbPb = bfrac(pr_fits, npr_fits, bin, "N_Jpsi_NPrompt_PbPb");
   RooRealVar Np_PbPb = bfrac(pr_fits, npr_fits, bin, "N_Jpsi_Prompt_PbPb");
   RooRealVar Nnp_PP = bfrac(pr_fits, npr_fits, bin, "N_Jpsi_NPrompt_PP");
   RooRealVar Np_PP = bfrac(pr_fits, npr_fits, bin, "N_Jpsi_Prompt_PP");

   double val = alpha.getVal() * (Nnp_PbPb.getVal()/Nnp_PP.getVal()) / (Np_PbPb.getVal()/Np_PP.getVal());
   double error = val * sqrt(pow(alpha.getError()/alpha.getVal(),2) 
         + pow(Nnp_PbPb.getError()/Nnp_PbPb.getVal(),2)
         + pow(Np_PbPb.getError()/Np_PbPb.getVal(),2)
         + pow(Nnp_PP.getError()/Nnp_PP.getVal(),2)
         + pow(Np_PP.getError()/Np_PP.getVal(),2)
         );
   RooRealVar ans("ans","",val);
   ans.setError(error);
   return ans;
}

RooRealVar alphaval(const char* pr_fits, const char* npr_fits, anabin bin, bool is2S) {
   bin.setcentbin(binI(0,200));
   RooRealVar f = bfrac(pr_fits, npr_fits, bin, Form("Bfrac%s_PP",is2S ? "_Psi2S" : ""));
   // then efficiencies
   TString feffjpsipr(Form("../Efficiency/files/histos_%s_pp.root",is2S ? "psi2s" : "jpsi"));
   TString feffnpr("../Efficiency/files/histos_npjpsi_pp.root");
   TString numname("hnumptdepcut_"); numname += (bin.centbin() == binI(0,200)) ? "pt" : "cent"; numname += (bin.rapbin() == binF(0,1.6)) ? "mid" : "fwd";
   TString denname("hnum_"); denname += (bin.centbin() == binI(0,200)) ? "pt" : "cent"; denname += (bin.rapbin() == binF(0,1.6)) ? "mid" : "fwd";
   TFile *tfeffjpsipr = TFile::Open(feffjpsipr);
   TFile *tfeffnpr = TFile::Open(feffnpr);

   TH1F *hnumjpsipr = (TH1F*) tfeffjpsipr->Get(numname);
   TH1F *hnumnpr = (TH1F*) tfeffnpr->Get(numname);
   TH1F *hdenjpsipr = (TH1F*) tfeffjpsipr->Get(denname);
   TH1F *hdennpr = (TH1F*) tfeffnpr->Get(denname);

   double numjpsipr, numnpr, denjpsipr, dennpr;

   if (bin == anabin(0,1.6,6.5,30,0,200) || bin == anabin(1.6,2.4,3,30,0,200)) { // Min.Bias bin
      hnumjpsipr = integrateHist(hnumjpsipr); numjpsipr = hnumjpsipr->GetBinContent(1);
      hnumnpr = integrateHist(hnumnpr); numnpr = hnumnpr->GetBinContent(1);
      hdenjpsipr = integrateHist(hdenjpsipr); denjpsipr = hdenjpsipr->GetBinContent(1);
      hdennpr = integrateHist(hdennpr); dennpr = hdennpr->GetBinContent(1);
   } else {
      for (int i=1; i<=hnumjpsipr->GetNbinsX(); i++) { // loop over the histogram bins to find the one we want
         if ((bin.centbin()==binI(0,200) && bin.ptbin()==binF(hnumjpsipr->GetXaxis()->GetBinLowEdge(i),hnumjpsipr->GetXaxis()->GetBinUpEdge(i))) // this is the bin we're looking for (pt)
               || bin.centbin()==binI(2*hnumjpsipr->GetXaxis()->GetBinLowEdge(i),2*hnumjpsipr->GetXaxis()->GetBinUpEdge(i))) { // this is the bin we're looking for (centrality)
            numjpsipr = hnumjpsipr->GetBinContent(i);
            numnpr = hnumnpr->GetBinContent(i);
            denjpsipr = hdenjpsipr->GetBinContent(i);
            dennpr = hdennpr->GetBinContent(i);
            break;
         }
      }
   }

   delete tfeffjpsipr; tfeffjpsipr=0;
   delete tfeffnpr; tfeffnpr=0;

   double effjpsipr = numjpsipr/denjpsipr;
   double effnpr = numnpr/dennpr;

   // bin.print();
   // cout << "effnpr = " << numnpr << "/" << dennpr << ", effpr = " << numjpsipr << "/" << denjpsipr << ", f = " << f.getVal() << endl;
   double val = (effnpr/effjpsipr) * f.getVal() / (1-f.getVal());
   double error = val*sqrt(2.)*f.getError()/f.getVal();
   RooRealVar ans("ans", "", val);
   ans.setError(error);
   return ans;
}

#endif // ifndef bfrac_h
