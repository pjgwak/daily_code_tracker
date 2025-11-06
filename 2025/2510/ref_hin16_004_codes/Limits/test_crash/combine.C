#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooExtendPdf.h"

#include <list>
#include <string>
#include <iostream>

using namespace RooFit;
using namespace std;

RooWorkspace* combine(const char* name_pbpb, const char* name_pp) {
   RooFit::PrintLevel(0);
   RooTrace::active(kTRUE);
   RooTrace::verbose(kTRUE);
   RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG); 
   // const char *poiname="N_{#Upsilon(3S)}";
   TFile *f = new TFile(name_pbpb) ;
   TFile *f_pp = new TFile(name_pp) ;

   // Retrieve workspace from file
   RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
   RooWorkspace* ws_pp = (RooWorkspace*) f_pp->Get("workspace");

   // RooRealVar *theVar; 
   RooDataSet *data; RooAbsPdf *pdf;
   RooRealVar *theVar_pp; RooDataSet *data_pp; RooAbsPdf *pdf_pp, *pdf_pbpb;

   // theVar = ws->var(poiname);
   // pdf = ws->pdf("pdf");
   data =(RooDataSet *) ws->data("data_PbPb");
   data_pp =(RooDataSet *) ws_pp->data("data_PP");

   // make the combined dataset
	RooCategory sample("sample", "sample");
	sample.defineType("PbPb");
	sample.defineType("PP");

   RooRealVar *invMass = ws->var("invMass");
   // RooRealVar *invMass = (RooRealVar*) ws->var("invMass")->Clone("invMass");
   // RooRealVar *invMass_pp = ws_pp->var("invMass");
   // add missing binnings
   // list<string> binningNames = invMass_pp->getBinningNames();
   // for (list<string>::const_iterator it=binningNames.begin(); it!=binningNames.end(); it++) {
   //    if (invMass->hasBinning(it->c_str())) continue;
   //    invMass->setBinning(invMass_pp->getBinning(it->c_str()));
   // }

   RooArgSet cols(*invMass);

	RooDataSet data_combo("data", "data", cols, RooFit::Index(sample),
/*Only for track rotation*/
		RooFit::Import("PbPb", *data), RooFit::Import("PP", *data_pp));

   RooWorkspace *wcombo = new RooWorkspace("wcombo","workspace for PbPb + pp");
   wcombo->import(data_combo);

   // build the pp pdf
   // pdf_pp = ws_pp->pdf("pdfMASS_Tot_PP");
   RooAbsPdf *sig1S_PP = ws_pp->pdf("sig_nonorm_PP");
   RooAbsPdf *sig2S_PP = ws_pp->pdf("sig2_nonorm_PP");
   RooAbsPdf *pdf_combinedbkgd_PP = ws_pp->pdf("bkg_nonorm_PP");
   RooRealVar *nsig1f_PP = ws_pp->var("Nsig_PP");
   RooRealVar *frac_PP = ws_pp->var("frac_PP");
   RooFormulaVar *nsig2f_PP = new RooFormulaVar("Nsig2_PP","@0*@1",RooArgList(*frac_PP,*nsig1f_PP));
   RooRealVar *nbkgd_PP = ws_pp->var("Nbkg_PP");

   pdf_pp = new RooAddPdf ("tot_PP","new total p.d.f. PP",
         RooArgList(*sig1S_PP,*sig2S_PP,*pdf_combinedbkgd_PP),
         RooArgList(*nsig1f_PP,*nsig2f_PP,*nbkgd_PP));

   // build the pbpb pdf
   RooAbsPdf *sig1S_PbPb = ws->pdf("sig_nonorm_PbPb");
   RooAbsPdf *sig2S_PbPb = ws->pdf("sig2_nonorm_PbPb");
   RooAbsPdf *pdf_combinedbkgd_PbPb = ws->pdf("bkg_nonorm_PbPb");
   RooRealVar *nsig1f_PbPb = ws->var("Nsig_PbPb");
   RooRealVar *frac_PbPb = ws->var("frac_PbPb");
   RooRealVar *frac = new RooRealVar("frac","frac",1,-10,10);
   RooFormulaVar *nsig2f_PbPb = new RooFormulaVar("Nsig2_PbPb","@0*@1*@2",RooArgList(*frac,*frac_PbPb,*nsig1f_PP));
   RooRealVar *nbkgd_PbPb = ws->var("Nbkg_PbPb");

   pdf_pbpb = new RooAddPdf ("tot_PbPb","new total p.d.f. PbPb",
         RooArgList(*sig1S_PbPb,*sig2S_PbPb,*pdf_combinedbkgd_PbPb),
         RooArgList(*nsig1f_PbPb,*nsig2f_PbPb,*nbkgd_PbPb));

	RooSimultaneous simPdf("simPdf", "simPdf", sample);
	simPdf.addPdf(*pdf_pbpb, "PbPb");
	simPdf.addPdf(*pdf_pp, "PP");

   // import PDFs to the workspace
   wcombo->import(*pdf_pp,	RecycleConflictNodes());
   wcombo->import(*pdf_pbpb, RecycleConflictNodes());
   wcombo->import(simPdf, RecycleConflictNodes());

   simPdf.fitTo(data_combo); // crashes sometimes but not always?? adding Range("MassWindow") or NumCPU(2) improves stability
   // simPdf.fitTo(data_combo,NumCPU(3),Range("MassWindow"));

   RooPlot* xframeSB = invMass->frame();
   data_combo.plotOn(xframeSB,Cut("sample==sample::PP"));
   simPdf.plotOn(xframeSB,Slice(sample,"PP"),ProjWData(sample,data_combo));
   TCanvas *cPP = new TCanvas();
   cPP->cd(); xframeSB->Draw();
   cPP->SaveAs("cPP.pdf");
   xframeSB = invMass->frame();
   data_combo.plotOn(xframeSB,Cut("sample==sample::PbPb"));
   simPdf.plotOn(xframeSB,Slice(sample,"PbPb"),ProjWData(sample,data_combo));
   TCanvas *cPbPb = new TCanvas();
   cPbPb->cd(); xframeSB->Draw();
   cPbPb->SaveAs("cPbPb.pdf");

   wcombo->writeToFile("fitresult_combo.root");
   return wcombo;
}
