#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"

#include <list>
#include <string>
#include <iostream>

using namespace RooFit;
using namespace std;

RooWorkspace* test_combine(const char* name_pbpb="fitresult.root", const char* name_pp="fitresult_pp.root", bool commonParams=false, int nCPU=2)
{
   // RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG); 
   TFile *f = TFile::Open(name_pbpb) ;
   TFile *f_pp = TFile::Open(name_pp) ;

   // Retrieve workspace from file
   RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
   RooWorkspace* ws_pp = (RooWorkspace*) f_pp->Get("workspace");

   RooDataSet *data=NULL;
   RooDataSet *data_pp=NULL;

   // let's massively clone everything into the new workspace
   RooWorkspace *wcombo = new RooWorkspace("workspace","workspace for PbPb + pp");

   // first, variables
   RooArgSet allVars = ws->allVars(); allVars.add(ws_pp->allVars());
   TIterator* it = allVars.createIterator();
   RooRealVar *theVar = (RooRealVar*) it->Next();
   while (theVar) {
      // restrict the range to [val-5sigma, val+5sigma]
      double newmin=theVar->getMin(), newmax=theVar->getMax();
      double val=theVar->getVal();
      double err=theVar->getError();
      TString varName(theVar->GetName());
      if ( (varName != "invMass") && (varName != "pt") && (varName != "rap") && (varName != "cent")) {
         newmin = max(val-5.*fabs(err), newmin);
         newmax = min(val+5.*fabs(err), newmax);
      }

      RooRealVar *theVarCopy = new RooRealVar(varName.Data(),theVar->GetTitle(),val,newmin,newmax);
      theVarCopy->setConstant(theVar->isConstant());
      theVarCopy->setError(err);
      if ( (theVar->getMin() == theVar->getMax()) || theVar->getError() == 0. ) theVarCopy->setConstant();
      wcombo->import(*theVarCopy);
      theVar = (RooRealVar*) it->Next();
   }


   // then, the formulae. Here we assume that they are all products.
   RooArgSet allFunctions = ws_pp->allFunctions(); allFunctions.add(ws->allFunctions());
   it = allFunctions.createIterator();
   RooFormulaVar *theFunc = (RooFormulaVar*) it->Next();
   while (theFunc) {
      // skip the number of Psi2S in PbPb and PP as we'll build them later (using the double ratio)
      if (TString(theFunc->GetName())=="N_Psi2S_PbPb") {theFunc = (RooFormulaVar*) it->Next(); continue;}
      RooArgSet* theVars = theFunc->getVariables();
      RooArgSet* theVars_copy = new RooArgSet();
      TIterator* it2 = theVars->createIterator();
      theVar = (RooRealVar*) it2->Next();
      while (theVar) {
         RooRealVar *theVar_copy = wcombo->var(theVar->GetName());
         if (!theVar_copy) {cout << "Error, missing variable " << theVar->GetName() << endl; theVar_copy = theVar;}
         theVars_copy->add(*theVar_copy);
         theVar = (RooRealVar*) it2->Next();
      }
      TString args("@0"); for (int i=1; i<theVars->getSize(); i++) args+=Form("*@%i",i);
      RooFormulaVar *theFuncCopy = new RooFormulaVar(theFunc->GetName(),args.Data(),RooArgList(*theVars_copy));
      wcombo->import(*theFuncCopy);
      theFunc = (RooFormulaVar*) it->Next();
   }

   // at last, import all pdf's
   RooArgSet allPdfs = ws_pp->allPdfs(); allPdfs.add(ws->allPdfs());
   int nok=0;
   while (allPdfs.getSize() != nok) {
      cout << nok << endl;
      nok=0;
      it = allPdfs.createIterator();
      RooAbsPdf *thePdf = (RooAbsPdf*) it->Next();
      while (thePdf) {
         if (wcombo->pdf(thePdf->GetName())) {thePdf = (RooAbsPdf*) it->Next(); nok++; continue;}
         TString className(thePdf->ClassName());
         ostringstream oss; thePdf->printArgs(oss); TString theargs(oss.str().c_str());
         if (className == "RooCBShape") {
            TString m,m0,sigma,alpha,n;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("m=")!=kNPOS) m=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("m0=")!=kNPOS) m0=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("sigma=")!=kNPOS) sigma=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("alpha=")!=kNPOS) alpha=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String();
            }
            cout << Form("CBShape::%s(%s, %s, %s, %s, %s)",thePdf->GetName(),m.Data(),m0.Data(),sigma.Data(),alpha.Data(),n.Data()) << endl;
            wcombo->factory(Form("CBShape::%s(%s, %s, %s, %s, %s)",thePdf->GetName(),m.Data(),m0.Data(),sigma.Data(),alpha.Data(),n.Data()));
         } else if (className == "RooExtendPdf") {
            TString pdf,n;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("pdf=")!=kNPOS) pdf=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String();
            }
            cout << Form("RooExtendPdf::%s(%s, %s)",thePdf->GetName(),pdf.Data(),n.Data()) << endl;
            wcombo->factory(Form("RooExtendPdf::%s(%s, %s)",thePdf->GetName(),pdf.Data(),n.Data()));
         } else if (className == "RooChebychev") {
            TString obs,args;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("coefList=")!=kNPOS) args=((TObjString*) t.Tokenize("=")->At(1))->String();
            }
            args.ReplaceAll("(","{"); args.ReplaceAll(")","}");
            cout << Form("Chebychev::%s(%s, %s)",thePdf->GetName(),obs.Data(),args.Data()) << endl;
            wcombo->factory(Form("Chebychev::%s(%s, %s)",thePdf->GetName(),obs.Data(),args.Data()));
         } else if (className == "RooUniform") {
           TString obs;
           TString t; Int_t from = 0;
           while (theargs.Tokenize(t, from , " ")) {
             if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
           }
           obs.ReplaceAll("(",""); obs.ReplaceAll(")","");
           cout << Form("Uniform::%s(%s)",thePdf->GetName(),obs.Data()) << endl;
           wcombo->factory(Form("Uniform::%s(%s)",thePdf->GetName(),obs.Data()));
         } else if (className == "RooAddPdf") {
            // skip the total PDFs as we will rebuild them
            if (TString(thePdf->GetName()).Index("pdfMASS_Tot_")!=kNPOS) {thePdf = (RooAbsPdf*) it->Next(); nok++; continue;}
            RooArgList pdfList, coefList, pdfList_orig=((RooAddPdf*)thePdf)->pdfList(), coefList_orig=((RooAddPdf*)thePdf)->coefList();
            for (int i=0; i<pdfList_orig.getSize(); i++) {
               RooAbsPdf *toadd = wcombo->pdf(pdfList_orig[i].GetName());
               if (toadd) pdfList.add(*toadd);
            }
            for (int i=0; i<coefList_orig.getSize(); i++) {
               RooAbsReal *coef = wcombo->var(coefList_orig[i].GetName());
               if (!coef) coef = wcombo->function(coefList_orig[i].GetName());
               if (!coef) {cout << "ERROR! Cannot find " << coefList_orig[i].GetName() << endl; return NULL;}
               coefList.add(*coef);
            }
            RooAddPdf *newpdf=NULL;
            if (pdfList.getSize()<coefList.getSize()) {thePdf = (RooAbsPdf*) it->Next(); continue;}
            if (coefList.getSize()>0) newpdf = new RooAddPdf(thePdf->GetName(),thePdf->GetTitle(),pdfList,coefList);
            else newpdf = new RooAddPdf(thePdf->GetName(),thePdf->GetTitle(),pdfList);
            wcombo->import(*newpdf);
         } else {
            cout << "Error, I don't know this type of PDF: " << className << "::" << thePdf->GetName() << ". I will ignore it." << endl;
         }
         nok++;
         thePdf = (RooAbsPdf*) it->Next();
      }
   }



   // D A T A
   data =(RooDataSet *) ws->data("dOS_DATA_PbPb");
   data_pp =(RooDataSet *) ws_pp->data("dOS_DATA_PP");

	RooCategory sample("sample", "sample");
	sample.defineType("PbPb");
	sample.defineType("PP");
   RooRealVar *invMass = wcombo->var("invMass");
   RooArgSet cols(*invMass);

   // hand-copy the pp and pbpb datasets
   // RooDataSet data_pp_copy("data_pp_copy","data_pp_copy",cols);

   RooDataSet data_combo("dOS_DATA", "dOS_DATA", cols, RooFit::Index(sample),
         /*Only for track rotation*/
         RooFit::Import("PbPb", *data), RooFit::Import("PP", *data_pp));

   wcombo->import(data_combo);

   // build the double ratio
   RooRealVar *RFrac2Svs1S_PbPbvsPP = new RooRealVar("RFrac2Svs1S_PbPbvsPP","RFrac2Svs1S_PbPbvsPP",0.5,-10,10);
   // help RooFit by telling it what to expect for the double ratio
   RooRealVar *RFrac2Svs1S_PbPb = wcombo->var("RFrac2Svs1S_PbPb");
   RooRealVar *RFrac2Svs1S_PP = wcombo->var("RFrac2Svs1S_PP");
   RFrac2Svs1S_PbPbvsPP->setVal(RFrac2Svs1S_PbPb->getVal()/RFrac2Svs1S_PP->getVal());
   RFrac2Svs1S_PbPbvsPP->setError(
         sqrt(pow(RFrac2Svs1S_PbPb->getError()/RFrac2Svs1S_PbPb->getVal(),2) + 
            pow(RFrac2Svs1S_PP->getError()/RFrac2Svs1S_PP->getVal(),2))
         *RFrac2Svs1S_PbPbvsPP->getVal());
  
   // build the pp pdf
  RooAbsPdf *sig1S_PP = wcombo->pdf("pdfMASS_Jpsi_PP");
  RooAbsPdf *sig2S_PP = wcombo->pdf("pdfMASS_Psi2S_PP");
  RooAbsPdf *pdf_combinedbkgd_PP = wcombo->pdf("pdfMASS_Bkg_PP");
  RooRealVar *nsig1f_PP = wcombo->var("N_Jpsi_PP");
  RooRealVar *nbkgd_PP = wcombo->var("N_Bkg_PP");
  RooFormulaVar *nsig2f_PP = (RooFormulaVar*) wcombo->function("N_Psi2S_PP");
  RooAddPdf *pdf_pp = new RooAddPdf ("pdfMASS_Tot_PP","new total p.d.f. PP",
                                     RooArgList(*sig1S_PP,*sig2S_PP,*pdf_combinedbkgd_PP),
                                     RooArgList(*nsig1f_PP,*nsig2f_PP,*nbkgd_PP));

   // build the pbpb pdf
   RooAbsPdf *sig1S = wcombo->pdf("pdfMASS_Jpsi_PbPb");
   RooAbsPdf *sig2S = wcombo->pdf("pdfMASS_Psi2S_PbPb");
   RooAbsPdf *pdf_combinedbkgd = wcombo->pdf("pdfMASS_Bkg_PbPb");
   RooRealVar *nsig1f = wcombo->var("N_Jpsi_PbPb");
   RooRealVar *nbkgd = wcombo->var("N_Bkg_PbPb");
   RooFormulaVar *nsig2f = new RooFormulaVar("N_Psi2S_PbPb","@0*@1*@2",RooArgList(*RFrac2Svs1S_PbPbvsPP,*RFrac2Svs1S_PP,*nsig1f));
   RooAbsPdf *pdf_pbpb = new RooAddPdf ("pdfMASS_Tot_PbPb","new total p.d.f.",
         RooArgList(*sig1S,*sig2S,*pdf_combinedbkgd),
         RooArgList(*nsig1f,*nsig2f,*nbkgd));

   // build the simunltaneous pdf
	 RooSimultaneous simPdf("simPdf", "simPdf", sample);
	 simPdf.addPdf(*pdf_pbpb, "PbPb");
	 simPdf.addPdf(*pdf_pp, "PP");

   // import PDFs to the workspace
   if (!commonParams) {
      wcombo->import(*pdf_pp,	RecycleConflictNodes());
   } else {
      wcombo->import(*pdf_pp,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PP","f_Jpsi"), RenameVariable("m_Jpsi_PP","m_Jpsi"), RenameVariable("sigma1_Jpsi_PP","sigma1_Jpsi"), RenameVariable("rSigma21_Jpsi_PP","rSigma21_Jpsi"));
   }
   if (!commonParams) {
      wcombo->import(*pdf_pbpb, RecycleConflictNodes());
   } else {
      wcombo->import(*pdf_pbpb,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PbPb","f_Jpsi"), RenameVariable("m_Jpsi_PbPb","m_Jpsi"), RenameVariable("sigma1_Jpsi_PbPb","sigma1_Jpsi"), RenameVariable("rSigma21_Jpsi_PbPb","rSigma21_Jpsi"));
   }
   wcombo->import(simPdf, RecycleConflictNodes());

   wcombo->Print();
   // return NULL;

   ws->Delete(); ws_pp->Delete(); f->Delete(); f_pp->Delete();

   simPdf.fitTo(data_combo,NumCPU(nCPU), Extended(kTRUE), Minimizer("Minuit2","Migrad"));

   // wcombo->writeToFile("fitresult_combo.root");
   return wcombo;
}
