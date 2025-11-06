#include "Macros/Utilities/monster.h"

#include "Math/Functor.h"
#include "RooFunctorBinding.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFunctorBinding.h" 
#include "TF1.h"
#include "RooTFnBinding.h"

#include <iostream>

using namespace std;

void testMonster() {
   RooWorkspace *ws = new RooWorkspace("ws","ws");

   RooArgList vars;
   RooRealVar njpsi_pbpb_pr("njpsi_pbpb_pr","",50,0,200); vars.add(njpsi_pbpb_pr);
   RooRealVar rfrac_pbpb_pr("rfrac_pbpb_pr","",1,0,2); vars.add(rfrac_pbpb_pr);
   RooRealVar njpsi_pp_pr("njpsi_pp_pr","",100,0,200); vars.add(njpsi_pp_pr);
   RooRealVar rfrac_pp_pr("rfrac_pp_pr","",1,0,2); vars.add(rfrac_pp_pr);
   RooRealVar njpsi_pbpb_npr("njpsi_pbpb_npr","",1,0,2); vars.add(njpsi_pbpb_npr);
   RooRealVar rfrac_pbpb_npr("rfrac_pbpb_npr","",1,0,2); vars.add(rfrac_pbpb_npr);
   RooRealVar njpsi_pp_npr("njpsi_pp_npr","",1,0,2); vars.add(njpsi_pp_npr);
   RooRealVar rfrac_pp_npr("rfrac_pp_npr","",1,0,2); vars.add(rfrac_pp_npr);
   RooRealVar effjpsi_pp_pr("effjpsi_pp_pr","",0.9,0,1); vars.add(effjpsi_pp_pr);
   RooRealVar effpsip_pp_pr("effpsip_pp_pr","",0.9,0,1); vars.add(effpsip_pp_pr);
   RooRealVar effjpsi_pp_npr("effjpsi_pp_npr","",0.2,0,1); vars.add(effjpsi_pp_npr);
   RooRealVar alpha_eff("alpha_eff","",1,0,2); vars.add(alpha_eff);
   RooRealVar alpha_fitpp("alpha_fitpp","",1,0,2); vars.add(alpha_fitpp);
   RooRealVar alpha_fitpbpb("alpha_fitpbpb","",1,0,2); vars.add(alpha_fitpbpb);
   RooRealVar alpha_effjpsi_pp_pr("alpha_effjpsi_pp_pr","",1,0,2); vars.add(alpha_effjpsi_pp_pr);
   RooRealVar alpha_effpsip_pp_pr("alpha_effpsip_pp_pr","",1,0,2); vars.add(alpha_effpsip_pp_pr);
   RooRealVar alpha_effjpsi_pp_npr("alpha_effjpsi_pp_npr","",1,0,2); vars.add(alpha_effjpsi_pp_npr);
   RooRealVar alpha_effpppbpb("alpha_effpppbpb","",1,0,2); vars.add(alpha_effpppbpb);

   double xx[npar_monster];
   for (int i=0; i<npar_monster; i++) {xx[i] = ((RooRealVar*) &(vars[i]))->getVal(); cout << xx[i] << endl;}
   cout << monsterFunction(xx) << endl;

   // ROOT::Math::Functor func(&monsterFunction,npar_monster);
   // RooAbsReal *monster = RooFit::bindFunction("monster",func,vars);
   // RooFunctorBinding *monster = new RooFunctorBinding("monster","",func,vars);
   // TF1 *fmonster = new TF1("tfmonster",monsterFunction2,0,1,npar_monster);
   // RooAbsReal *monster = RooFit::bindFunction(fmonster,njpsi_pbpb_pr,vars);
   // ws->import(*fmonster);
   // ws->import(*monster);

   cout << monsterFunction(ws,vars);
   ws->writeToFile("test.root");
}
