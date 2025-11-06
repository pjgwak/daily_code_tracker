#ifndef resultUtils_h
#define resultUtils_h

#include "bin.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "TString.h"
#include "TFile.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "RooFitResult.h"
#include "TFitResult.h"
#include "TTree.h"

#include <string>
#include <vector>
#include <iostream>

///////////////
// CONSTANTS //
///////////////

// NB: luminosities in mub-1
const double lumipp = 28.0e6;
const double lumipbpb_ABCD = 351*1.049;
const double lumipbpb_peri = 464;
const double NMB = 2.366003e9*1.049;

using namespace std;

RooRealVar* poiFromFile(const char* filename, const char* token, const char* thepoiname);
RooRealVar* poiFromWS(RooWorkspace* ws, const char* token, const char* thepoiname);
RooAbsPdf* pdfFromWS(RooWorkspace* ws, const char* token, const char* thepdfname);
RooAbsData* dataFromWS(RooWorkspace* ws, const char* token, const char* thedataname);
vector<TString> fileList(const char* input, const char* token="", const char* DSTag="DATA", const char* prependPath="", const char* fitType="");
vector<TString> combFileList(const char* input, const char* token="", const char* prependPath="");
vector<TString> limitsFileList(const char* input, const char* token="", const char* prependPath="");
RooRealVar* ratioVar(RooRealVar *num, RooRealVar *den, bool usedenerror=true);
anabin binFromFile(const char* filename);
bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck, bool override=true);
bool binok(anabin thecat, string xaxis, anabin &tocheck, bool override=true);
bool isSameBinPPPbPb(const char* filenamePbPb, const char* filenamePP);
TString treeFileName(const char* workDirName, const char* DSTag="DATA", const char* prependPath="", const char* fitType = "");
double poiFromBin(const char* workDirName, const char* collSystem, const char* thepoiname, anabin thebin, const char* DSTag="DATA", const char* prependPath="");
double poiErrFromBin(const char* workDirName, const char* collSystem, const char* thepoiname, anabin thebin, const char* DSTag="DATA", const char* prependPath="");
void prune(vector<anabin> &v, bool keepshortest=true);
void prune(TGraphAsymmErrors *g, TGraphAsymmErrors *gsyst=NULL, bool keepshortest=true);
void prune(TGraphErrors *g, bool keepshortest=true);
TGraphAsymmErrors* result12007_mid_cent();
TGraphAsymmErrors* result12007_fwd_cent();
TGraphAsymmErrors* result12007_mid_cent_syst();
TGraphAsymmErrors* result12007_fwd_cent_syst();
TGraphAsymmErrors* result12007_mid();
TGraphAsymmErrors* result12007_fwd();
TGraphAsymmErrors* result12007_mid_syst();
TGraphAsymmErrors* result12007_fwd_syst();

void results2tree(
      const char* workDirName, 
      const char* DSTag="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
      const char* prependPath="",
      const char* fitType = "",
      const char* thePoiNames="N_Jpsi,b_Jpsi,f_Jpsi,m_Jpsi,sigma1_Jpsi,alpha_Jpsi,n_Jpsi,sigma2_Jpsi,rSigma21_Jpsi,"
      "lambda1_Bkg,lambda2_Bkg,lambda3_Bkg,lambda4_Bkg,lambda5_Bkg,N_Bkg,b_Bkg,"
      "ctau1_CtauRes,ctau2_CtauRes,f_CtauRes,rSigma21_CtauRes,s1_CtauRes,"
      "fDFSS_BkgNoPR,fDLIV_BkgNoPR,lambdaDDS_BkgNoPR,lambdaDF_BkgNoPR,lambdaDSS_BkgNoPR,lambdaDSS_JpsiNoPR,"
      "eff,effnp,acc,accnp,lumi,taa,ncoll,npart,correl_N_Jpsi_vs_b_Jpsi",
      bool wantPureSMC=false
      );
#include "../../results2tree.C"

RooRealVar* poiFromFile(const char* filename, const char* token, const char* thepoiname) {
   TFile *f = TFile::Open(filename);
   if (!f) {
      cout << "Error, file " << filename << " does not exist." << endl;
      return NULL;
   }
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   if (!ws) {
      cout << "Error, file " << filename << " is bad." << endl;
      return NULL;
   }
   RooRealVar *ans = poiFromWS(ws, token, thepoiname);
   if (!ans) return NULL;
   TString poiname_and_token = TString(thepoiname) + TString(token);
   RooRealVar* ansc = new RooRealVar(*ans,poiname_and_token + Form("_from_%s",filename));
   delete ws;
   delete f;
   return ansc;
}

RooRealVar* poiFromWS(RooWorkspace* ws, const char* token, const char* thepoiname) {
   TString poiname_and_token = TString(thepoiname) + TString(token);
   RooRealVar *ans = (RooRealVar*) ws->var(poiname_and_token);
   if (!ans && TString(thepoiname) == "N_Psi2S") { // we want the number of psi2S but it's not in the ws, let's compute it
      RooRealVar *njpsi = (RooRealVar*) ws->var(Form("N_Jpsi%s",token));
      // cout << Form("N_Jpsi%s",token) << endl;
      if (!njpsi) return ans;
      RooRealVar *rfrac = (RooRealVar*) ws->var(Form("RFrac2Svs1S%s",token));
      // cout << Form("RFrac2Svs1S%s",token) << endl;
      if (!rfrac) return ans;
      RooFitResult *fr = (RooFitResult*) ws->obj(Form("fitResult_pdfMASS_Tot%s", token));
      double corr=0;
      if (fr) {
         corr = fr->correlation(Form("N_Jpsi%s",token), Form("RFrac2Svs1S%s",token));
      } 
      double npsip = njpsi->getVal()*rfrac->getVal();
      double dnpsip = fabs(npsip) * sqrt(pow(njpsi->getError()/njpsi->getVal(),2)
            +pow(rfrac->getError()/rfrac->getVal(),2)
            +2.*fabs(corr*njpsi->getError()*rfrac->getError()/njpsi->getVal()/rfrac->getVal()));
      ans = new RooRealVar(Form("N_Psi2S%s",token),"",npsip);
      ans->setError(dnpsip);
   }
   return ans;
}

RooAbsPdf* pdfFromWS(RooWorkspace* ws, const char* token, const char* thepdfname) {
   TString pdfname_and_token = TString(thepdfname) + TString(token);
   RooAbsPdf *ans = (RooAbsPdf*) ws->pdf(pdfname_and_token);
   return ans;
}

RooAbsData* dataFromWS(RooWorkspace* ws, const char* token, const char* thedataname) {
   TString dataname_and_token = TString(thedataname) + TString(token);
   RooAbsData *ans = (RooAbsData*) ws->data(dataname_and_token);
   return ans;
}

vector<TString> fileList(const char* input, const char* token, const char* DSTag, const char* prependPath, const char* fitType) {
   vector<TString> ans;
  
   TString basedir("");
   if (!strcmp(fitType,"")) basedir = Form("Output/%s/result/%s/",input,DSTag);
   else basedir = Form("Output/%s/%s/%s/result/",input,fitType,DSTag);
  
   if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
   TSystemDirectory dir(input,basedir);

   TList *files = dir.GetListOfFiles();

   if (files) {
      TIter next(files);
      TSystemFile *file;
      TString fname;

      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (fname.EndsWith(".root") && fname.Index("FIT") != kNPOS
               && (TString(token) == "" || fname.Index(token) != kNPOS)) {
            ans.push_back(basedir+fname);
         }
      }
   }

   delete files;
   return ans;
}

vector<TString> combFileList(const char* input, const char* token, const char* prependPath) {
  vector<TString> ans;
  
  TString basedir(Form("CombinedWorkspaces/%s/",input));
  if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
  TSystemDirectory dir(input,basedir);
  
  TList *files = dir.GetListOfFiles();
  
  if (files) {
    TIter next(files);
    TSystemFile *file;
    TString fname;
    
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (fname.EndsWith(".root") && fname.Index("combined_PbPbPP_") != kNPOS
          && (fname.Index(token) != kNPOS)) {
        ans.push_back(basedir+fname);
      }
    }
  }
  
  delete files;
  return ans;
}

vector<TString> limitsFileList(const char* input, const char* token, const char* prependPath) {
  vector<TString> ans;
  
  TString basedir(Form("Output/%s/",input));
  if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
  TSystemDirectory dir(input,basedir);
  
  TList *files = dir.GetListOfFiles();
  
  if (files) {
    TIter next(files);
    TSystemFile *file;
    TString fname;
    
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (fname.EndsWith(".root") && fname.Index("combined_PbPbPP_") != kNPOS
          && (fname.Index("_Scan") == kNPOS) && (fname.Index(token) != kNPOS)) {
        ans.push_back(basedir+fname);
      }
    }
  }
  
  delete files;
  return ans;
}

RooRealVar* ratioVar(RooRealVar *num, RooRealVar *den, bool usedenerror) {
   double n = num->getVal();
   double d = den->getVal();
   double dn = num->getError();
   double dd = den->getError();

   double r = d!=0 ? n/d : 0;
   double dr = n!=0 && d!=0 ? fabs(r * sqrt(pow(dn/n,2) + pow(dd/d,2))) : 0;
   if (!usedenerror && n!=0) dr = fabs((dn/n)*r);
   RooRealVar *ans = new RooRealVar(Form("%s_over_%s",num->GetName(),den->GetName()), Form("%s / %s",num->GetTitle(),den->GetTitle()), r);
   ans->setError(dr);

   return ans;
}

anabin binFromFile(const char* filename) {
   TFile *f = TFile::Open(filename);
   if (!f) {
      cout << "Error, file " << filename << " does not exist." << endl;
      return anabin(0,0,0,0,0,0);
   }
   RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
   if (!ws) {
      cout << "Error, file " << filename << " is bad." << endl;
      return anabin(0,0,0,0,0,0);
   }
   RooRealVar *pt = (RooRealVar*) ws->var("pt");
   RooRealVar *rap = (RooRealVar*) ws->var("rap");
   RooRealVar *cent = (RooRealVar*) ws->var("cent");
   if (!pt || !rap || !cent) {
      cout << "Error, file " << filename << " is bad." << endl;
      return anabin(0,0,0,0,0,0);
   }
   anabin ans(rap->getMin(),rap->getMax(),pt->getMin(),pt->getMax(),cent->getMin(),cent->getMax());
   delete ws;
   delete f;
   return ans;
}

bool binok(vector<anabin> thecats, string xaxis, anabin &tocheck, bool override) {
   bool ok=false;

   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      if (xaxis=="pt" && it->rapbin()==tocheck.rapbin() && it->centbin()==tocheck.centbin()
            && ! (it->ptbin()==tocheck.ptbin())) {
         ok=true;
         if (override) tocheck.setptbin(it->ptbin());
         break;
      } else if (xaxis=="cent" && it->rapbin()==tocheck.rapbin() && it->ptbin()==tocheck.ptbin()
            && ! (it->centbin()==tocheck.centbin())) {
         ok=true;
         if (override) tocheck.setcentbin(it->centbin());
         break;
      } else if (xaxis=="rap" && it->centbin()==tocheck.centbin() && it->ptbin()==tocheck.ptbin()
            && ! (it->rapbin()==tocheck.rapbin())) {
         ok=true;
         if (override) tocheck.setrapbin(it->rapbin());
         break;
      } else if ((it->centbin().low()<=0 && it->centbin().high()<=0)
            && it->rapbin()==tocheck.rapbin() && it->ptbin()==tocheck.ptbin()
            &&  (abs(it->centbin().low())==abs(tocheck.centbin().low()) && abs(it->centbin().high())==abs(tocheck.centbin().high()))) {
         ok=true;
         break;
      }
   }

   return ok;
}

bool binok(anabin thecat, string xaxis, anabin &tocheck, bool override) {
   vector<anabin> thecats; thecats.push_back(thecat);
   return binok(thecats, xaxis, tocheck, override);
}

bool isSameBinPPPbPb(const char* filenamePbPb, const char* filenamePP)
{
  TString fPbPb(filenamePbPb);
  TString fPP(filenamePP);
  
  anabin thebin_PbPb = binFromFile(fPbPb.Data());
  anabin thebin_PP = binFromFile(fPP.Data());
  
  if ( fPbPb.Contains("_PbPb_") && fPP.Contains("_PP_"))
  {
    if ( (thebin_PbPb.ptbin().low() == thebin_PP.ptbin().low()) && (thebin_PbPb.ptbin().high() == thebin_PP.ptbin().high())  && (thebin_PbPb.rapbin().low() == thebin_PP.rapbin().low()) && (thebin_PbPb.rapbin().high() == thebin_PP.rapbin().high()) ) return true;
    else return false;
  }
  else
  {
    cout << "[ERROR]: Comparing bins of same collision system" << endl;
    return false;
  }
}

TString treeFileName(const char* workDirName, const char* DSTag, const char* prependPath, const char* fitType) {
   TString outputFileName("");
   if (!strcmp(fitType,"")) outputFileName = Form("Output/%s/result/%s/tree_allvars.root",workDirName,DSTag);
   else outputFileName = Form("Output/%s/%s/%s/result/tree_allvars.root",workDirName,fitType,DSTag);
   if ( strcmp(prependPath,"") ) outputFileName.Prepend(Form("%s/",prependPath));
   return outputFileName;
}

double poiFromBin(const char* workDirName, const char* theCollSystem, const char* thepoiname, anabin thebin, const char* DSTag, const char* prependPath) {
   TString tfname = treeFileName(workDirName, DSTag, prependPath);
   TFile *f = TFile::Open(tfname);
   if (!f || !f->IsOpen()) {
      results2tree(workDirName,DSTag,prependPath);
      f = new TFile(tfname);
      if (!f) return -1e99;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return -1e99;

   // fix centrality for pp
   if (TString(theCollSystem) == "PP") thebin.setcentbin(binI(0,200));

   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   float val=-999.;
   int valI=-999.;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);
   if (!TString(thepoiname).Contains("chi2") && !TString(thepoiname).Contains("ndof")) tr->SetBranchAddress(Form("%s_val",thepoiname),&val);
   else if (!TString(thepoiname).Contains("ndof")) tr->SetBranchAddress(thepoiname,&val);
   else tr->SetBranchAddress(thepoiname,&valI);
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   bool found=false;
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);
      if ((anabin(ymin, ymax, ptmin, ptmax, centmin, centmax) == thebin) && (TString(collSystem) == TString(theCollSystem))) {
         found=true;
         break;
      }
   }
   f->Close(); delete f;
   if (!found) val = -999;
   if (TString(thepoiname).Contains("ndof")) val = valI;
   return val;
}

double poiErrFromBin(const char* workDirName, const char* theCollSystem, const char* thepoiname, anabin thebin, const char* DSTag, const char* prependPath) {
   TString tfname = treeFileName(workDirName, DSTag, prependPath);
   TFile *f = TFile::Open(tfname);
   if (!f || !f->IsOpen()) {
      results2tree(workDirName,DSTag,prependPath);
      f = new TFile(tfname);
      if (!f) return -1e99;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return -1e99;

   // fix centrality for pp
   if (TString(theCollSystem) == "PP") thebin.setcentbin(binI(0,200));

   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   float err=0;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);
   tr->SetBranchAddress(Form("%s_errL",thepoiname),&err); //Only take low error (for most of pois, except Taa, Npart or Ncoll the error is symmetric)
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   bool found=false;
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);
      if ((anabin(ymin, ymax, ptmin, ptmax, centmin, centmax) == thebin) && (TString(collSystem) == TString(theCollSystem))) {
         found=true;
         break;
      }
   }
   f->Close(); delete f;
   if (!found) err=-999;
   return err;
}

void prune(vector<anabin> &v, bool keepshort) {
   vector<anabin> ans;

   vector<anabin>::const_iterator it1, it2;
   for (it1=v.begin(); it1!=v.end(); it1++) {
      bool binok=true;
      for (it2=v.begin(); it2!=v.end(); it2++) {
         if (*it1==*it2) continue;
         if (it1->rapbin()==it2->rapbin() && it1->ptbin()==it2->ptbin()) {
            binI cb1 = it1->centbin();
            binI cb2 = it2->centbin();
            if (!(cb1==binI(0,200)) && cb1.low()==cb2.low()) { // the bin is not MB and there is another bin with the same lower edge
               if (keepshort && cb1.high()>cb2.high()) binok=false;
               if (!keepshort && cb1.high()<cb2.high()) binok=false;
            }
         } // same pt and rap bins
      } // for it2
      if (binok) ans.push_back(*it1);
   } // for it1

   v = ans;
}

void prune(TGraphAsymmErrors *g, TGraphAsymmErrors *gsyst, bool keepshort) {
   int n = g->GetN();
   for (int i1=0; i1<n; i1++) {
      double xl1 = g->GetX()[i1]-g->GetErrorXlow(i1);
      double xh1 = g->GetX()[i1]+g->GetErrorXhigh(i1);
      bool binok=true;
      for (int i2=0; i2<n; i2++) {
         if (i2==i1) continue;
         double xl2 = g->GetX()[i2]-g->GetErrorXlow(i2);
         double xh2 = g->GetX()[i2]+g->GetErrorXhigh(i2);
         if (fabs(xl1-xl2)<1e-3 && keepshort && xh1>xh2) binok=false;
         if (fabs(xl1-xl2)<1e-3 && !keepshort && xh1<xh2) binok=false;
      } // for i2
      if (!binok) {
         g->SetPoint(i1,g->GetX()[i1],-g->GetY()[i1]);
         if (gsyst) gsyst->SetPoint(i1,gsyst->GetX()[i1],-gsyst->GetY()[i1]);
      }
   } // for i1
}

void prune(TGraphErrors *g, bool keepshort) {
   int n = g->GetN();
   for (int i1=0; i1<n; i1++) {
      double xl1 = g->GetX()[i1]-g->GetErrorX(i1);
      double xh1 = g->GetX()[i1]+g->GetErrorX(i1);
      bool binok=true;
      for (int i2=0; i2<n; i2++) {
         if (i2==i1) continue;
         double xl2 = g->GetX()[i2]-g->GetErrorX(i2);
         double xh2 = g->GetX()[i2]+g->GetErrorX(i2);
         if (fabs(xl1-xl2)<1e-3 && keepshort && xh1>xh2) binok=false;
         if (fabs(xl1-xl2)<1e-3 && !keepshort && xh1<xh2) binok=false;
      } // for i2
      if (!binok) g->SetPoint(i1,g->GetX()[i1],-g->GetY()[i1]);
   } // for i1
}

TGraphAsymmErrors* result12007_mid_cent() {
   double x[3] = {308.4-5, 158.6-5, 32.8-5};
   double y[3] = {0.43, 0.67, 0.};
   double ex[3] = {0, 0, 0};
   double ey[3] = {0.2, 0.2, 0.47};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(3,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_mid_cent");
   return ans;
}
TGraphAsymmErrors* result12007_fwd_cent() {
   double x[3] = {308.4+5, 158.6+5, 32.8+5};
   double y[3] = {2.31, 0.93, 0.89};
   double ex[3] = {0, 0, 0};
   double ey[3] = {0.53, 0.47, 0.39};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(3,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_fwd_cent");
   return ans;
}
TGraphAsymmErrors* result12007_mid_cent_syst() {
   double x[3] = {308.4-5, 158.6-5, 32.8-5};
   double y[3] = {0.43, 0.67, 0.};
   double ex[3] = {5, 5, 5};
   double ey[3] = {0.07, 0.09, 0.47};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(3,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_mid_cent_syst");
   return ans;
}
TGraphAsymmErrors* result12007_fwd_cent_syst() {
   double x[3] = {308.4+5, 158.6+5, 32.8+5};
   double y[3] = {2.31, 0.93, 0.89};
   double ex[3] = {5, 5, 5};
   double ey[3] = {0.37, 0.27, 0.17};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(3,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_fwd_cent_syst");
   return ans;
}
TGraphAsymmErrors* result12007_mid() {
   double x[1] = {200};
   double y[1] = {0.45};
   double ex[1] = {0};
   double ey[1] = {0.13};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(1,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_mid");
   return ans;
}
TGraphAsymmErrors* result12007_fwd() {
   double x[1] = {200};
   double y[1] = {1.67};
   double ex[1] = {0};
   double ey[1] = {0.34};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(1,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_fwd");
   return ans;
}
TGraphAsymmErrors* result12007_mid_syst() {
   double x[1] = {200};
   double y[1] = {0.45};
   double ex[1] = {25};
   double ey[1] = {0.07};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(1,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_mid_syst");
   return ans;
}
TGraphAsymmErrors* result12007_fwd_syst() {
   double x[1] = {200};
   double y[1] = {1.67};
   double ex[1] = {25};
   double ey[1] = {0.27};
   TGraphAsymmErrors *ans = new TGraphAsymmErrors(1,x,y,ex,ex,ey,ey);
   ans->SetName("graph_12007_fwd_syst");
   return ans;
}
#endif // ifndef resultUtils_h
