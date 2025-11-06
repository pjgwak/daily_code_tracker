#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void checkFit3sigma(const char* filename) {
   TFile *f = new TFile(filename);
   if (!f) return;

   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) {delete f; return;}

   TString thePoiNames="N_Jpsi,b_Jpsi,f_Jpsi,m_Jpsi,sigma1_Jpsi,alpha_Jpsi,n_Jpsi,sigma2_Jpsi,MassRatio,rSigma21_Jpsi,"
      "lambda1_Bkg,lambda2_Bkg,lambda3_Bkg,lambda4_Bkg,lambda5_Bkg,N_Bkg,b_Bkg,"
      "ctau1_CtauRes,ctau2_CtauRes,f_CtauRes,rSigma21_CtauRes,s1_CtauRes,"
      "fDFSS_BkgNoPR,fDLIV_BkgNoPR,lambdaDDS_BkgNoPR,lambdaDF_BkgNoPR,lambdaDSS_BkgNoPR,lambdaDSS_JpsiNoPR,"
      "eff,effnp,lumi,taa,ncoll,npart,correl_N_Jpsi_vs_b_Jpsi";

   TString t; Int_t from = 0;
   while (thePoiNames.Tokenize(t, from , ",")) {
      TString cutstr = t + "_err>0 && " + t + "_min!=" + t + "_max && ((" + t + "_val - " + t + "_min)/" + t + "_err<3 || (" + t + "_max - " + t + "_val)/" + t + "_err<3)";
      TString colstr = "collSystem:ymin:ymax:ptmin:ptmax:centmin:centmax:" + t + "_val:" + t + "_err:" + t + "_min:" + t + "_max";
      cout << cutstr.Data() << endl;
      tr->Scan(colstr,cutstr);
   }

   f->Close();
   delete f;
}
