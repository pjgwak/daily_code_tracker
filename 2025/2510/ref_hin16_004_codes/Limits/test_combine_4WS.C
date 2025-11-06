#include "RooGlobalFunc.h"
#include "RooWorkspace.h"
#include "RooExtendPdf.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"

#include "../Fitter/Macros/Utilities/monster.h"
#include "../Fitter/Systematics/syst.h"
#include "../Fitter/Systematics/syst.C"
#include "../Efficiency/eff.h"

#include <list>
#include <string>
#include <iostream>

using namespace RooFit;
using namespace std;


RooWorkspace* test_combine_4WS(const char* name_pbpb_pr, const char* name_pp_pr, const char* name_pbpb_npr, const char* name_pp_npr, bool commonParams=false, int nCPU=2)
{
   // RooMsgService::instance().addStream(RooFit::MsgLevel::DEBUG); 
   TFile *f_pbpb_pr = TFile::Open(name_pbpb_pr) ;
   TFile *f_pp_pr = TFile::Open(name_pp_pr) ;
   TFile *f_pbpb_npr = TFile::Open(name_pbpb_npr) ;
   TFile *f_pp_npr = TFile::Open(name_pp_npr) ;

   // Retrieve workspace from file
   RooWorkspace* ws_pbpb_pr = (RooWorkspace*) f_pbpb_pr->Get("workspace");
   RooWorkspace* ws_pp_pr = (RooWorkspace*) f_pp_pr->Get("workspace");
   RooWorkspace* ws_pbpb_npr = (RooWorkspace*) f_pbpb_npr->Get("workspace");
   RooWorkspace* ws_pp_npr = (RooWorkspace*) f_pp_npr->Get("workspace");

   RooDataSet *data_pbpb_pr=NULL;
   RooDataSet *data_pp_pr=NULL;
   RooDataSet *data_pbpb_npr=NULL;
   RooDataSet *data_pp_npr=NULL;

   // let's massively clone everything into the new workspace
   RooWorkspace *wcombo = new RooWorkspace("workspace","workspace for PbPb + pp");

   // first, variables
   // prompt
   RooArgSet allVars = ws_pbpb_pr->allVars(); 
   allVars.add(ws_pp_pr->allVars());
   TIterator* it = allVars.createIterator();
   RooRealVar *theVar = (RooRealVar*) it->Next();
   while (theVar) {
      // restrict the range to [val-5sigma, val+5sigma]
      double newmin=theVar->getMin(), newmax=theVar->getMax();
      double val=theVar->getVal();
      double err=theVar->getError();
      TString varName(theVar->GetName());
      TString varTitle(theVar->GetTitle());
      if ( (varName != "invMass") && (varName != "pt") && (varName != "rap") && (varName != "cent")) {
         newmin = max(val-5.*fabs(err), newmin);
         newmax = min(val+5.*fabs(err), newmax);
         varName = Form("%s_pass",varName.Data());
         varTitle = Form("%s (pass)",theVar->GetTitle());
      }

      RooRealVar *theVarCopy = new RooRealVar(varName.Data(),varTitle.Data(),val,newmin,newmax);
      theVarCopy->setConstant(theVar->isConstant());
      theVarCopy->setError(err);
      if ( (theVar->getMin() == theVar->getMax()) || theVar->getError() == 0. ) theVarCopy->setConstant();
      wcombo->import(*theVarCopy);
      theVar = (RooRealVar*) it->Next();
   }

   // non-prompt
   RooArgSet allVars_np = ws_pbpb_npr->allVars();
   allVars_np.add(ws_pp_npr->allVars());
   it = allVars_np.createIterator();
   theVar = (RooRealVar*) it->Next();
   while (theVar) {
      // restrict the range to [val-5sigma, val+5sigma]
      double newmin=theVar->getMin(), newmax=theVar->getMax();
      double val=theVar->getVal();
      double err=theVar->getError();
      TString varName(theVar->GetName());
      TString varTitle(theVar->GetTitle());
      if ( (varName != "invMass") && (varName != "pt") && (varName != "rap") && (varName != "cent")) {
         newmin = max(val-5.*fabs(err), newmin);
         newmax = min(val+5.*fabs(err), newmax);
         varName = Form("%s_fail",varName.Data());
         varTitle = Form("%s (fail)",theVar->GetTitle());
      }

      RooRealVar *theVarCopy = new RooRealVar(varName.Data(),varTitle.Data(),val,newmin,newmax);
      theVarCopy->setConstant(theVar->isConstant());
      theVarCopy->setError(err);
      if ( (theVar->getMin() == theVar->getMax()) || theVar->getError() == 0. ) theVarCopy->setConstant();
      wcombo->import(*theVarCopy);
      theVar = (RooRealVar*) it->Next();
   }


   // then, the formulae. Here we assume that they are all products.
   // prompt
   RooArgSet allFunctions = ws_pp_pr->allFunctions(); allFunctions.add(ws_pbpb_pr->allFunctions());
   it = allFunctions.createIterator();
   RooFormulaVar *theFunc = (RooFormulaVar*) it->Next();
   while (theFunc) {
      TString funcName(theFunc->GetName());
      // skip the number of Psi2S in PbPb and PP as we'll build them later (using the double ratio)
      if (funcName=="N_Psi2S_PbPb") {theFunc = (RooFormulaVar*) it->Next(); continue;}
      RooArgSet* theVars = theFunc->getVariables();
      RooArgSet* theVars_copy = new RooArgSet();
      TIterator* it2 = theVars->createIterator();
      theVar = (RooRealVar*) it2->Next();
      while (theVar) {
         const char* varname_pr = Form("%s_pass",theVar->GetName());
         RooRealVar *theVar_copy = wcombo->var(varname_pr);
         if (!theVar_copy) {cout << "Error, missing variable " << varname_pr << endl; theVar_copy = theVar; theVar->SetName(varname_pr);}
         theVars_copy->add(*theVar_copy);
         theVar = (RooRealVar*) it2->Next();
      }
      TString args("@0"); for (int i=1; i<theVars->getSize(); i++) args+=Form("*@%i",i);
      RooFormulaVar *theFuncCopy = new RooFormulaVar(Form("%s_pass",theFunc->GetName()),args.Data(),RooArgList(*theVars_copy));
      wcombo->import(*theFuncCopy);
      theFunc = (RooFormulaVar*) it->Next();
   }

   // non-prompt
   RooArgSet allFunctions_np = ws_pp_npr->allFunctions(); allFunctions_np.add(ws_pbpb_npr->allFunctions());
   it = allFunctions_np.createIterator();
   theFunc = (RooFormulaVar*) it->Next();
   while (theFunc) {
      RooArgSet* theVars = theFunc->getVariables();
      RooArgSet* theVars_copy = new RooArgSet();
      TIterator* it2 = theVars->createIterator();
      theVar = (RooRealVar*) it2->Next();
      while (theVar) {
         const char* varname_npr = Form("%s_fail",theVar->GetName());
         RooRealVar *theVar_copy = wcombo->var(varname_npr);
         if (!theVar_copy) {cout << "Error, missing variable " << varname_npr << endl; theVar_copy = theVar; theVar->SetName(varname_npr);}
         theVars_copy->add(*theVar_copy);
         theVar = (RooRealVar*) it2->Next();
      }
      TString args("@0"); for (int i=1; i<theVars->getSize(); i++) args+=Form("*@%i",i);
      RooFormulaVar *theFuncCopy = new RooFormulaVar(Form("%s_fail",theFunc->GetName()),args.Data(),RooArgList(*theVars_copy));
      wcombo->import(*theFuncCopy);
      theFunc = (RooFormulaVar*) it->Next();
   }

   // at last, import all pdf's
   // prompt
   RooArgSet allPdfs = ws_pp_pr->allPdfs(); allPdfs.add(ws_pbpb_pr->allPdfs());
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
               if (t.Index("m0=")!=kNPOS) m0=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
               if (t.Index("sigma=")!=kNPOS) sigma=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
               if (t.Index("alpha=")!=kNPOS) alpha=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
            }
            cout << Form("CBShape::%s_pass(%s, %s, %s, %s, %s)",thePdf->GetName(),m.Data(),m0.Data(),sigma.Data(),alpha.Data(),n.Data()) << endl;
            wcombo->factory(Form("CBShape::%s_pass(%s, %s, %s, %s, %s)",thePdf->GetName(),m.Data(),m0.Data(),sigma.Data(),alpha.Data(),n.Data()));
         } else if (className == "RooExtendPdf") {
            TString pdf,n;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("pdf=")!=kNPOS) pdf=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String() + "_pass";
            }
            cout << Form("RooExtendPdf::%s_pass(%s, %s)",thePdf->GetName(),pdf.Data(),n.Data()) << endl;
            wcombo->factory(Form("RooExtendPdf::%s_pass(%s, %s)",thePdf->GetName(),pdf.Data(),n.Data()));
         } else if (className == "RooChebychev") {
            TString obs,args;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("coefList=")!=kNPOS) args=((TObjString*) t.Tokenize("=")->At(1))->String();
            }
            args.ReplaceAll("(","{"); args.ReplaceAll(")","}");
            args.ReplaceAll("_PP",TString("_PP")+ "_pass"); args.ReplaceAll("_PbPb",TString("_PbPb")+"_pass");
            cout << Form("Chebychev::%s_pass(%s, %s)",thePdf->GetName(),obs.Data(),args.Data()) << endl;
            wcombo->factory(Form("Chebychev::%s_pass(%s, %s)",thePdf->GetName(),obs.Data(),args.Data()));
         } else if (className == "RooUniform") {
           TString obs;
           TString t; Int_t from = 0;
           while (theargs.Tokenize(t, from , " ")) {
             if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
           }
           obs.ReplaceAll("(",""); obs.ReplaceAll(")","");
           cout << Form("Uniform::%s_pass(%s)",thePdf->GetName(),obs.Data()) << endl;
           wcombo->factory(Form("Uniform::%s_pass(%s)",thePdf->GetName(),obs.Data()));
         } else if (className == "RooAddPdf") {
            // skip the total PDFs as we will rebuild them
            if (TString(thePdf->GetName()).Index("pdfMASS_Tot_")!=kNPOS) {thePdf = (RooAbsPdf*) it->Next(); nok++; continue;}
            RooArgList pdfList, coefList, pdfList_orig=((RooAddPdf*)thePdf)->pdfList(), coefList_orig=((RooAddPdf*)thePdf)->coefList();
            for (int i=0; i<pdfList_orig.getSize(); i++) {
               RooAbsPdf *toadd = wcombo->pdf(Form("%s_pass",pdfList_orig[i].GetName()));
               if (toadd) pdfList.add(*toadd);
            }
            for (int i=0; i<coefList_orig.getSize(); i++) {
               RooAbsReal *coef = wcombo->var(Form("%s_pass",coefList_orig[i].GetName()));
               if (!coef) coef = wcombo->function(Form("%s_pass",coefList_orig[i].GetName()));
               if (!coef) {cout << "ERROR! Cannot find " << Form("%s_pass",coefList_orig[i].GetName()) << endl; return NULL;}
               coefList.add(*coef);
            }
            RooAddPdf *newpdf=NULL;
            if (pdfList.getSize()<coefList.getSize()) {thePdf = (RooAbsPdf*) it->Next(); continue;}
            if (coefList.getSize()>0) newpdf = new RooAddPdf(Form("%s_pass",thePdf->GetName()),Form("%s (pass)",thePdf->GetTitle()),pdfList,coefList);
            else newpdf = new RooAddPdf(Form("%s_pass",thePdf->GetName()),Form("%s (pass)",thePdf->GetTitle()),pdfList);
            wcombo->import(*newpdf);
         } else {
            cout << "Error, I don't know this type of PDF: " << className << "::" << thePdf->GetName() << ". I will ignore it." << endl;
         }
         nok++;
         thePdf = (RooAbsPdf*) it->Next();
      }
   }

   // non-prompt
   RooArgSet allPdfs_np = ws_pp_npr->allPdfs(); allPdfs_np.add(ws_pbpb_npr->allPdfs());
   nok=0;
   while (allPdfs_np.getSize() != nok) {
      cout << nok << endl;
      nok=0;
      it = allPdfs_np.createIterator();
      RooAbsPdf *thePdf = (RooAbsPdf*) it->Next();
      while (thePdf) {
         if (wcombo->pdf(Form("%s_fail",thePdf->GetName()))) {thePdf = (RooAbsPdf*) it->Next(); nok++; continue;}
         TString className(thePdf->ClassName());
         ostringstream oss; thePdf->printArgs(oss); TString theargs(oss.str().c_str());
         if (className == "RooCBShape") {
            TString m,m0,sigma,alpha,n;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("m=")!=kNPOS) m=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("m0=")!=kNPOS) m0=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
               if (t.Index("sigma=")!=kNPOS) sigma=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
               if (t.Index("alpha=")!=kNPOS) alpha=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
            }
            TString factorystring = Form("CBShape::%s_fail(%s, %s, %s, %s, %s)",thePdf->GetName(),m.Data(),m0.Data(),sigma.Data(),alpha.Data(),n.Data());
            cout << factorystring.Data() << endl;
            wcombo->factory(factorystring.Data());
         } else if (className == "RooExtendPdf") {
            TString pdf,n;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("pdf=")!=kNPOS) pdf=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
               if (t.Index("n=")!=kNPOS) n=((TObjString*) t.Tokenize("=")->At(1))->String() + "_fail";
            }
            TString factorystring = Form("RooExtendPdf::%s_fail(%s, %s)",thePdf->GetName(),pdf.Data(),n.Data());
            cout << factorystring.Data() << endl;
            wcombo->factory(factorystring.Data());
         } else if (className == "RooChebychev") {
            TString obs,args;
            TString t; Int_t from = 0;
            while (theargs.Tokenize(t, from , " ")) {
               if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
               if (t.Index("coefList=")!=kNPOS) args=((TObjString*) t.Tokenize("=")->At(1))->String();
            }
            args.ReplaceAll("(","{"); args.ReplaceAll(")","}");
            args.ReplaceAll("_PP",TString("_PP")+ "_fail"); args.ReplaceAll("_PbPb",TString("_PbPb")+"_fail");
            TString factorystring = Form("Chebychev::%s_fail(%s, %s)",thePdf->GetName(),obs.Data(),args.Data());
            cout << factorystring.Data() << endl;
            wcombo->factory(factorystring.Data());
         } else if (className == "RooUniform") {
           TString obs;
           TString t; Int_t from = 0;
           while (theargs.Tokenize(t, from , " ")) {
             if (t.Index("x=")!=kNPOS) obs=((TObjString*) t.Tokenize("=")->At(1))->String();
           }
           obs.ReplaceAll("(",""); obs.ReplaceAll(")","");
           cout << Form("Uniform::%s(%s)",thePdf->GetName(),obs.Data()) << endl;
           wcombo->factory(Form("Uniform::%s_fail(%s)",thePdf->GetName(),obs.Data()));
         } else if (className == "RooAddPdf") {
            // skip the total PDFs as we will rebuild them
            if (TString(thePdf->GetName()).Index("pdfMASS_Tot_")!=kNPOS) {thePdf = (RooAbsPdf*) it->Next(); nok++; continue;}
            RooArgList pdfList, coefList, pdfList_orig=((RooAddPdf*)thePdf)->pdfList(), coefList_orig=((RooAddPdf*)thePdf)->coefList();
            for (int i=0; i<pdfList_orig.getSize(); i++) {
               RooAbsPdf *toadd = wcombo->pdf(Form("%s_fail",pdfList_orig[i].GetName()));
               if (toadd) pdfList.add(*toadd);
            }
            for (int i=0; i<coefList_orig.getSize(); i++) {
               RooAbsReal *coef = wcombo->var(Form("%s_fail",coefList_orig[i].GetName()));
               if (!coef) coef = wcombo->function(Form("%s_fail",coefList_orig[i].GetName()));
               if (!coef) {cout << "ERROR! Cannot find " << Form("%s_fail",coefList_orig[i].GetName()) << endl; return NULL;}
               coefList.add(*coef);
            }
            RooAddPdf *newpdf=NULL;
            if (pdfList.getSize()<coefList.getSize()) {thePdf = (RooAbsPdf*) it->Next(); continue;}
            if (coefList.getSize()>0) newpdf = new RooAddPdf(Form("%s_fail",thePdf->GetName()),Form("%s (fail)",thePdf->GetTitle()),pdfList,coefList);
            else newpdf = new RooAddPdf(Form("%s_fail",thePdf->GetName()),Form("%s (fail)",thePdf->GetTitle()),pdfList);
            wcombo->import(*newpdf);
         } else {
            cout << "Error, I don't know this type of PDF: " << className << "::" << thePdf->GetName() << ". I will ignore it." << endl;
         }
         nok++;
         thePdf = (RooAbsPdf*) it->Next();
      }
   }


   // D A T A
   data_pbpb_pr =(RooDataSet *) ws_pbpb_pr->data("dOS_DATA_PbPb");
   data_pp_pr =(RooDataSet *) ws_pp_pr->data("dOS_DATA_PP");
   data_pbpb_npr =(RooDataSet *) ws_pbpb_npr->data("dOS_DATA_PbPb");
   data_pp_npr =(RooDataSet *) ws_pp_npr->data("dOS_DATA_PP");

   RooCategory sample("sample", "sample");
	 sample.defineType("PbPb_pass");
	 sample.defineType("PP_pass");
	 sample.defineType("PbPb_fail");
	 sample.defineType("PP_fail");
   RooRealVar *invMass = wcombo->var("invMass");
   RooArgSet cols(*invMass);

   // hand-copy the pp and pbpb datasets
   // RooDataSet data_pp_copy("data_pp_copy","data_pp_copy",cols);

   RooDataSet data_combo("dOS_DATA", "dOS_DATA", cols, RooFit::Index(sample),
         /*Only for track rotation*/
         RooFit::Import("PbPb_pass", *data_pbpb_pr), RooFit::Import("PP_pass", *data_pp_pr),
         RooFit::Import("PbPb_fail", *data_pbpb_npr), RooFit::Import("PP_fail", *data_pp_npr));

   wcombo->import(data_combo);


   // build the double ratio prompt
   RooRealVar *RFrac2Svs1S_PbPbvsPP_P = new RooRealVar("RFrac2Svs1S_PbPbvsPP_P","RFrac2Svs1S_PbPbvsPP_P",0.5,-10,10);
  
   // help RooFit by telling it what to expect for the double ratio prompt
   anabin thebin = binFromFile(name_pbpb_pr);
   double effjpsi_pp_P_val = ljpsieff("jpsi","pp",thebin,"../Efficiency");
   double effpsip_pp_P_val = ljpsieff("psi2s","pp",thebin,"../Efficiency");
   double effjpsi_pp_NP_val = ljpsieff("npjpsi","pp",thebin,"../Efficiency");
   const double vars[18] = {wcombo->var("N_Jpsi_PbPb_pass")->getVal(),wcombo->var("RFrac2Svs1S_PbPb_pass")->getVal(),wcombo->var("N_Jpsi_PP_pass")->getVal(),wcombo ->var("RFrac2Svs1S_PP_pass")->getVal(),wcombo->var("N_Jpsi_PbPb_fail")->getVal(),wcombo->var("RFrac2Svs1S_PbPb_fail")->getVal(),wcombo->var("N_Jpsi_PP_fail")->getVal() ,wcombo->var("RFrac2Svs1S_PP_fail")->getVal(),effjpsi_pp_P_val,effpsip_pp_P_val,effjpsi_pp_NP_val,0.,0.,0.,0.,0.,0.,0.};
    RFrac2Svs1S_PbPbvsPP_P->setVal(doubleratio_prompt(vars));
   wcombo->import(*RFrac2Svs1S_PbPbvsPP_P);
  
   // build the prompt pp pdf
   RooAbsPdf *sig1S_PP_pass = wcombo->pdf("pdfMASS_Jpsi_PP_pass");
   RooAbsPdf *sig2S_PP_pass = wcombo->pdf("pdfMASS_Psi2S_PP_pass");
   RooAbsPdf *pdf_combinedbkgd_PP_pass = wcombo->pdf("pdfMASS_Bkg_PP_pass");
   RooRealVar *nsig1f_PP_pass = wcombo->var("N_Jpsi_PP_pass");
   RooRealVar *nbkgd_PP_pass = wcombo->var("N_Bkg_PP_pass");
   RooFormulaVar *nsig2f_PP_pass = (RooFormulaVar*) wcombo->function("N_Psi2S_PP_pass");
   RooAddPdf *pdf_pp_pass = new RooAddPdf ("pdfMASS_Tot_PP_pass","new total pass p.d.f. PP",
                                           RooArgList(*sig1S_PP_pass,*sig2S_PP_pass,*pdf_combinedbkgd_PP_pass),
                                           RooArgList(*nsig1f_PP_pass,*nsig2f_PP_pass,*nbkgd_PP_pass));
  
   // build the non-prompt pp pdf
   RooAbsPdf *sig1S_PP_fail = wcombo->pdf("pdfMASS_Jpsi_PP_fail");
   RooAbsPdf *sig2S_PP_fail = wcombo->pdf("pdfMASS_Psi2S_PP_fail");
   RooAbsPdf *pdf_combinedbkgd_PP_fail = wcombo->pdf("pdfMASS_Bkg_PP_fail");
   RooRealVar *nsig1f_PP_fail = wcombo->var("N_Jpsi_PP_fail");
   RooRealVar *nbkgd_PP_fail = wcombo->var("N_Bkg_PP_fail");
   RooFormulaVar *nsig2f_PP_fail = (RooFormulaVar*) wcombo->function("N_Psi2S_PP_fail");
   RooAddPdf *pdf_pp_fail = new RooAddPdf ("pdfMASS_Tot_PP_fail","new total fail p.d.f. PP",
                                      RooArgList(*sig1S_PP_fail,*sig2S_PP_fail,*pdf_combinedbkgd_PP_fail),
                                      RooArgList(*nsig1f_PP_fail,*nsig2f_PP_fail,*nbkgd_PP_fail));
  
   // build the prompt pbpb pdf
   RooAbsPdf *sig1S_pass = wcombo->pdf("pdfMASS_Jpsi_PbPb_pass");
   RooAbsPdf *sig2S_pass = wcombo->pdf("pdfMASS_Psi2S_PbPb_pass");
   RooAbsPdf *pdf_combinedbkgd_pass = wcombo->pdf("pdfMASS_Bkg_PbPb_pass");
   RooRealVar *nsig1f_pass = wcombo->var("N_Jpsi_PbPb_pass");
   RooRealVar *nbkgd_pass = wcombo->var("N_Bkg_PbPb_pass");
  
   RooRealVar* effjpsi_pp_P = new RooRealVar("effjpsi_pp_P","effjpsi_pp_P",effjpsi_pp_P_val);
   RooRealVar* effpsip_pp_P = new RooRealVar("effpsip_pp_P","effpsip_pp_P",effpsip_pp_P_val);
   RooRealVar* effjpsi_pp_NP = new RooRealVar("effjpsi_pp_NP","effjpsi_pp_NP",effjpsi_pp_NP_val);
   Double_t Npsi2SPbPbPass = npsip_pbpb_pass_from_doubleratio_prompt(wcombo, RooArgList(*effjpsi_pp_P,*effpsip_pp_P,*effjpsi_pp_NP)); // Create and import the needed variables
  
   RooFormulaVar* nsig2f_pass = (RooFormulaVar*) wcombo->function("N_Psi2S_PbPb_pass");
   RooAbsPdf *pdf_pbpb_pass = new RooAddPdf ("pdfMASS_Tot_PbPb_pass","new total pass p.d.f. PbPb",
                                        RooArgList(*sig1S_pass,*sig2S_pass,*pdf_combinedbkgd_pass),
                                        RooArgList(*nsig1f_pass,*nsig2f_pass,*nbkgd_pass));
  
  // build the non-prompt pbpb pdf
  RooAbsPdf *sig1S_fail = wcombo->pdf("pdfMASS_Jpsi_PbPb_fail");
  RooAbsPdf *sig2S_fail = wcombo->pdf("pdfMASS_Psi2S_PbPb_fail");
  RooAbsPdf *pdf_combinedbkgd_fail = wcombo->pdf("pdfMASS_Bkg_PbPb_fail");
  RooRealVar *nsig1f_fail = wcombo->var("N_Jpsi_PbPb_fail");
  RooRealVar *nbkgd_fail = wcombo->var("N_Bkg_PbPb_fail");
  RooFormulaVar *nsig2f_fail = (RooFormulaVar*) wcombo->function("N_Psi2S_PbPb_fail");
  RooAddPdf *pdf_pbpb_fail = new RooAddPdf ("pdfMASS_Tot_PbPb_fail","new total fail p.d.f. PbPb",
                                          RooArgList(*sig1S_fail,*sig2S_fail,*pdf_combinedbkgd_fail),
                                          RooArgList(*nsig1f_fail,*nsig2f_fail,*nbkgd_fail));
  
   // build the simunltaneous pdf
   RooSimultaneous simPdf("simPdf", "simPdf", sample);
   simPdf.addPdf(*pdf_pbpb_pass, "PbPb_pass");
   simPdf.addPdf(*pdf_pbpb_fail, "PbPb_fail");
	 simPdf.addPdf(*pdf_pp_pass, "PP_pass");
   simPdf.addPdf(*pdf_pp_fail, "PP_fail");
  
   // import PDFs to the workspace
   if (!commonParams) {
     wcombo->import(*pdf_pp_pass,	RecycleConflictNodes());
     wcombo->import(*pdf_pp_fail,	RecycleConflictNodes());
   } else {
     wcombo->import(*pdf_pp_pass,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PP_pass","f_Jpsi_pass"), RenameVariable("m_Jpsi_PP_pass","m_Jpsi_pass"), RenameVariable("sigma1_Jpsi_PP_pass","sigma1_Jpsi_pass"), RenameVariable("rSigma21_Jpsi_PP_pass","rSigma21_Jpsi_pass"));
     wcombo->import(*pdf_pp_fail,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PP_fail","f_Jpsi_fail"), RenameVariable("m_Jpsi_PP_fail","m_Jpsi_fail"), RenameVariable("sigma1_Jpsi_PP_fail","sigma1_Jpsi_fail"), RenameVariable("rSigma21_Jpsi_PP_fail","rSigma21_Jpsi_fail"));
   }
   if (!commonParams) {
     wcombo->import(*pdf_pbpb_pass, RecycleConflictNodes());
     wcombo->import(*pdf_pbpb_fail, RecycleConflictNodes());
   } else {
     wcombo->import(*pdf_pbpb_pass,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PbPb_pass","f_Jpsi_pass"), RenameVariable("m_Jpsi_PbPb_pass","m_Jpsi_pass"), RenameVariable("sigma1_Jpsi_PbPb_pass","sigma1_Jpsi_pass"),  RenameVariable("rSigma21_Jpsi_PbPb_pass","rSigma21_Jpsi_pass"));
     wcombo->import(*pdf_pbpb_fail,	RecycleConflictNodes(), RenameVariable("f_Jpsi_PbPb_fail","f_Jpsi_fail"), RenameVariable("m_Jpsi_PbPb_fail","m_Jpsi_fail"), RenameVariable("sigma1_Jpsi_PbPb_fail","sigma1_Jpsi_fail"),  RenameVariable("rSigma21_Jpsi_PbPb_fail","rSigma21_Jpsi_fail"));
   }
   wcombo->import(simPdf, RecycleConflictNodes());
  
   wcombo->Print();

   ws_pbpb_pr->Delete(); ws_pp_pr->Delete(); ws_pbpb_npr->Delete(); ws_pp_npr->Delete(); f_pbpb_pr->Delete(); f_pp_pr->Delete(); f_pbpb_npr->Delete(); f_pp_npr->Delete();

   // simPdf.fitTo(data_combo,NumCPU(nCPU), Extended(kTRUE), Minimizer("Minuit2","Migrad"));

   // wcombo->writeToFile("fitresult_combo.root");
   return wcombo;
}
