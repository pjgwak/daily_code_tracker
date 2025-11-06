#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooMinuit.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SimpleInterval.h"
#include "TAxis.h"
#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TROOT.h>
#include <fstream>
#include <TGraph.h>
#include "TMath.h"
#include "TF1.h"

using namespace RooFit;
using namespace RooStats;

const bool dosyst = false;

void combinedWorkspace_fromSimFit(const char* name_in, const char* name_out){

   TFile *fin = TFile::Open(name_in);

   RooWorkspace * ws_in = (RooWorkspace*) fin->Get("workspace");
   if (!ws_in) return;
   // File.GetObject("wcombo", ws);
   // ws->Print();
   RooWorkspace *ws = new RooWorkspace("wcombo","workspace for PbPb + pp");

   // create the combined dataset
   RooDataSet * data_PP = (RooDataSet*) ws_in->data("dOS_DATA_PP");
   RooDataSet * data_PbPb = (RooDataSet*) ws_in->data("dOS_DATA_PbPb");
   RooCategory *sample =  ws_in->cat("sample");
   RooRealVar *invMass = ws_in->var("invMass");
	RooArgSet cols(*invMass);
   RooDataSet *data = new RooDataSet("dOS_DATA", "dOS_DATA", cols, RooFit::Index(*sample),
         RooFit::Import("PbPb", *data_PbPb), RooFit::Import("PP", *data_PP));
   ws->import(*data);
   ws->import(*(ws_in->pdf("simPdf")));


   // systematics
   ws->factory( "syst_PP_kappa[1.10]" );
   ws->factory( "expr::alpha_syst_PP('pow(syst_PP_kappa,beta_syst_PP)',syst_PP_kappa,beta_syst_PP[0,-5,5])" );
   // ws->factory( "SUM::syst_PP_nom(alpha_syst_PP*pdfMASS_Bkg_PP)" );
   ws->factory( "prod::N_Bkg_PP_syst(alpha_syst_PP,N_Bkg_PP)" );
   ws->factory( "Gaussian::constr_syst_PP(beta_syst_PP,glob_syst_PP[0,-5,5],1)" );
   ws->factory( "EDIT::pdfMASSTot_Bkg_PP_syst(pdfMASSTot_Bkg_PP, N_Bkg_PP=N_Bkg_PP_syst)");
   ws->factory( "EDIT::pdfMASSTot_PP_syst(pdfMASS_Tot_PP, pdfMASSTot_Bkg_PP = pdfMASSTot_Bkg_PP_syst)");
   ws->factory( "PROD::pdfMASSTot_Bkg_PP_norm(pdfMASSTot_Bkg_PP_syst, constr_syst_PP)");
   ws->factory( "EDIT::simPdf_syst(simPdf, pdfMASS_Tot_PP=pdfMASSTot_PP_syst)");

   /////////////////////////////////////////////////////////////////////
   RooRealVar * pObs = ws->var("invMass"); // get the pointer to the observable
   RooArgSet obs("observables");
   obs.add(*pObs);
   obs.add( *ws->cat("sample"));    
   //  /////////////////////////////////////////////////////////////////////
   ws->var("glob_syst_PP")->setConstant(true);
   // ws->var("glob_syst_pbpb")->setConstant(true);
   RooArgSet globalObs("global_obs");
   if (dosyst) {
      // globalObs.add( *ws->var("glob_syst_pbpb") );
      globalObs.add( *ws->var("glob_syst_PP") );
   }

   // ws->Print();

   RooArgSet poi("poi");
   poi.add( *ws->var("RFrac2Svs1S_PbPbvsPP") );

   RooRealVar* RFrac2Svs1S_PbPbvsPP = ws->var("RFrac2Svs1S_PbPbvsPP");
   RooRealVar* leftEdge = new RooRealVar("leftEdge","leftEdge",-10);
   RooRealVar* rightEdge = new RooRealVar("rightEdge","rightEdge",10);
   RooGenericPdf step("step", "step", "(@0 >= @1) && (@0 < @2)", RooArgList(*RFrac2Svs1S_PbPbvsPP, *leftEdge, *rightEdge));
   ws->import(step);
   ws->factory( "Uniform::flat(RFrac2Svs1S_PbPbvsPP)" );


   // create set of nuisance parameters
   RooArgSet nuis("nuis");
   // nuis.add( *ws->var("beta_syst_pbpb") );
   nuis.add( *ws->var("beta_syst_PP") );

   // set parameters constant
   RooArgSet allVars = ws->allVars();
   TIterator* it = allVars.createIterator();
   RooRealVar *theVar = (RooRealVar*) it->Next();
   while (theVar) {
      TString varname(theVar->GetName());
      if (varname != "RFrac2Svs1S_PbPbvsPP"
            && varname != "invMass"
            && varname != "sample"
            ) 
         theVar->setConstant();
      if (varname=="glob_syst_PP"
            || varname=="beta_syst_PP"
            ) {
         cout << varname << endl;
         theVar->setConstant(!dosyst);
      }
      theVar = (RooRealVar*) it->Next();
   }

   // create signal+background Model Config
   RooStats::ModelConfig sbHypo("SbHypo");
   sbHypo.SetWorkspace( *ws );
   sbHypo.SetPdf( *ws->pdf("simPdf_syst") );
   sbHypo.SetObservables( obs );
   sbHypo.SetGlobalObservables( globalObs );
   sbHypo.SetParametersOfInterest( poi );
   sbHypo.SetNuisanceParameters( nuis );
   sbHypo.SetPriorPdf( *ws->pdf("step") ); // this is optional

   // ws->Print();
   /////////////////////////////////////////////////////////////////////
   RooAbsReal * pNll = sbHypo.GetPdf()->createNLL( *data,NumCPU(2),Offset(kTRUE), Extended(kTRUE));
   RooMinuit(*pNll).migrad(); // minimize likelihood wrt all parameters before making plots
   RooPlot *framepoi = ((RooRealVar *)poi.first())->frame(Bins(10),Range(0.,1.),Title("LL and profileLL in RFrac2Svs1S_PbPbvsPP"));
   pNll->plotOn(framepoi,ShiftToZero());
   
   RooAbsReal * pProfile = pNll->createProfile( globalObs ); // do not profile global observables
   pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
   pProfile->plotOn(framepoi,LineColor(kRed));
   framepoi->SetMinimum(0);
   framepoi->SetMaximum(10);
   TCanvas *cpoi = new TCanvas();
   cpoi->cd(); framepoi->Draw();
   cpoi->SaveAs("cpoi.pdf");

   ((RooRealVar *)poi.first())->setMin(0.);
   RooArgSet * pPoiAndNuisance = new RooArgSet("poiAndNuisance");
   // pPoiAndNuisance->add(*sbHypo.GetNuisanceParameters());
   // pPoiAndNuisance->add(*sbHypo.GetParametersOfInterest());
   pPoiAndNuisance->add( nuis );
   pPoiAndNuisance->add( poi );
   sbHypo.SetSnapshot(*pPoiAndNuisance);

   RooPlot* xframeSB = pObs->frame(Title("SBhypo"));
   data->plotOn(xframeSB,Cut("sample==sample::PP"));
   RooAbsPdf *pdfSB = sbHypo.GetPdf();
   sample = ws->cat("sample");
   pdfSB->plotOn(xframeSB,Slice(*sample,"PP"),ProjWData(*sample,*data));
   TCanvas *c1 = new TCanvas();
   c1->cd(); xframeSB->Draw();
   c1->SaveAs("c1.pdf");
   RooPlot* xframeB = pObs->frame(Title("SBhypo_PbPb"));
   data->plotOn(xframeB,Cut("sample==sample::PbPb"));
   RooAbsPdf *pdfB = sbHypo.GetPdf();
   pdfB->plotOn(xframeB,Slice(*sample,"PbPb"),ProjWData(*sample,*data));
   TCanvas *c2 = new TCanvas();
   c2->cd(); xframeB->Draw();
   c2->SetLogy();
   c2->SaveAs("c2.pdf");
   c2->SaveAs("c2.root");

   delete pProfile;
   delete pNll;
   delete pPoiAndNuisance;
   ws->import( sbHypo );
   /////////////////////////////////////////////////////////////////////
   RooStats::ModelConfig bHypo = sbHypo;
   bHypo.SetName("BHypo");
   bHypo.SetWorkspace(*ws);
   pNll = bHypo.GetPdf()->createNLL( *data,NumCPU(2) );
   RooArgSet poiAndGlobalObs("poiAndGlobalObs");
   poiAndGlobalObs.add( poi );
   poiAndGlobalObs.add( globalObs );
   pProfile = pNll->createProfile( poiAndGlobalObs ); // do not profile POI and global observables
   ((RooRealVar *)poi.first())->setVal( 0 );  // set RFrac2Svs1S_PbPbvsPP=0 here
   pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
   pPoiAndNuisance = new RooArgSet( "poiAndNuisance" );
   pPoiAndNuisance->add( nuis );
   pPoiAndNuisance->add( poi );
   bHypo.SetSnapshot(*pPoiAndNuisance);

   delete pProfile;
   delete pNll;
   delete pPoiAndNuisance;

   // import model config into workspace
   bHypo.SetWorkspace(*ws);
   ws->import( bHypo );
   /////////////////////////////////////////////////////////////////////
   ws->Print();
   bHypo.Print();
   sbHypo.Print();

   // save workspace to file
   ws -> SaveAs(name_out);

   return;
}




