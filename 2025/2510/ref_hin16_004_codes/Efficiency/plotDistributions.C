//#include "Macros/CMS/CMS_lumi.C"
//#include "Macros/CMS/tdrstyle.C"
#include "../Fitter/Macros/Utilities/bin.h"
#include "../Fitter/Macros/Utilities/bfrac.h"
#include "../Fitter/Macros/Utilities/EVENTUTILS.h"
#include "../Fitter/Macros/Utilities/initClasses.h"
#include "../Fitter/Macros/Utilities/resultUtils.h"
#include "../Fitter/Macros/Utilities/texUtils.h"
#include "../Fitter/Macros/Utilities/monster.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TArrow.h"

using namespace std;

////////////////
// PARAMETERS //
////////////////

#ifndef poiname_check
#define poiname_check
const char* poiname       = "N_Jpsi";
#endif
bool  fiterrors     = true;  // statistical errors are from the fit
bool  doprompt      = true;  // prompt Jpsi
bool  dononprompt   = false;  // nonprompt Jpsi
string nameTag_base = "_mc";    // can put here e.g. "_data", "_mc", ...

//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void setOptions(bool afiterrors, bool adoprompt, bool adononprompt, string anameTag_base);
void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, bool isPP=false);
void plot(vector<anabin> thecats, string xaxis, string workDirName);
int color(int i);
int markerstyle(int i);
string nameTag;
bool isData = nameTag_base.find("_data")!=std::string::npos;

class distr_input {
public:
  double npp;
  double dnpp_stat;
  double naa;
  double dnaa_stat;
  double bfracpp;
  double dbfracpp;
  double bfracaa;
  double dbfracaa;
  double correlaa;
  double correlpp;
};


/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

void plotPt(string workDirName, int iplot) {
  string xaxis = "pt";
  vector<anabin> theCats;
  
  // 4 rapidity intervals
  if (iplot==0) {
    theCats.push_back(anabin(0,0.6,6.5,50,0,200));
    theCats.push_back(anabin(0.6,1.2,6.5,50,0,200));
    theCats.push_back(anabin(1.2,1.8,6.5,50,0,200));
    theCats.push_back(anabin(1.8,2.4,6.5,50,0,200));
  }
  
  // 3 centrality intervals
  if (iplot==1) {
    
    theCats.push_back(anabin(0,2.4,6.5,50,0,20));
    theCats.push_back(anabin(0,2.4,6.5,50,20,60));
    theCats.push_back(anabin(0,2.4,6.5,50,60,200));
  }
  
  // 1 rapidity interval
  if (iplot==2) {
    theCats.push_back(anabin(0,2.4,6.5,50,0,200));
  }
  
  plot(theCats,xaxis,workDirName);
};

void plotRap(string workDirName) {
  string xaxis = "rap";
  vector<anabin> theCats;
  
  theCats.push_back(anabin(0,2.4,6.5,50,0,200));
  
  plot(theCats,xaxis,workDirName);
};

void plotAll(string workDirName) {
  
  if (doprompt) nameTag_base += "_prompt";
  if (dononprompt) nameTag_base += "_nonprompt";
  nameTag = nameTag_base;
  
  plotPt(workDirName,0);
//  plotPt(workDirName,1);
  plotPt(workDirName,2);
//  if (isData) plotRap(workDirName);
};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plot(vector<anabin> thecats, string xaxis, string outputDir) {
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent"
  
  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  
  TString sampleName(isData ? "DATA" : (doprompt? "MCJPSIPR" : "MCJPSINOPR"));
  TFile *f = new TFile(treeFileName(outputDir.c_str(),sampleName.Data(),"../Fitter"));
  if (f) cout << "Using file : " << f->GetName() << endl;
  else{
    cout << "[ERROR] file not found" << endl;
    return;
  }
  /*if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),sampleName.Data(),"../Fitter");
    f = new TFile(treeFileName(outputDir.c_str(),sampleName.Data(),"../Fitter"));
    if (!f) return;
  }*/
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) {cout << "{ERROR] Tree not dound in file " << f->GetName() << endl; return;}
  
  map<anabin, distr_input> theVars_inputs;
  
  vector<double> x, ex, y, ey;
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  float val, err=0;
  float bfrac, bfrac_err, correl=0.;
  int ival=-999;
  char collSystem[5];
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("collSystem",collSystem);
  if (isData)
  {
    if (TString(poiname).Contains("N_Jpsi"))
    {
      tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
      tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&err);
    }
    else
    {
      tr->SetBranchAddress(Form("%s_val",poiname),&val);
      tr->SetBranchAddress(Form("%s_errL",poiname),&err);
    }
    tr->SetBranchAddress("b_Jpsi_val",&bfrac);
    tr->SetBranchAddress("b_Jpsi_errL",&bfrac_err);
    tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  }
  else
  {
    tr->SetBranchAddress(Form("%s_val",poiname),&val);
    tr->SetBranchAddress(Form("%s_err",poiname),&err);
  }
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);
    
    bool ispp = (TString(collSystem)=="PP");
    
    if (ispp) {
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp_stat = err;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_err;
      theVars_inputs[thebin].correlpp = correl;
    } else {
      theVars_inputs[thebin].naa = val;
      theVars_inputs[thebin].dnaa_stat = err;
      theVars_inputs[thebin].bfracaa = bfrac;
      theVars_inputs[thebin].dbfracaa = bfrac_err;
      theVars_inputs[thebin].correlaa = correl;
    }
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinnedPbPb;
  map<anabin, vector<double> > theVarsBinnedPbPb_stat;
  map<anabin, vector<double> > theVarsBinnedPP;
  map<anabin, vector<double> > theVarsBinnedPP_stat;
  map<anabin, TGraphAsymmErrors* > theGraphsPbPb;
  map<anabin, TGraphAsymmErrors* > theGraphsPP;
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinnedPbPb[*it] = vector<double>();
    theVarsBinnedPP[*it] = vector<double>();
    theVarsBinnedPbPb_stat[*it] = vector<double>();
    theVarsBinnedPP_stat[*it] = vector<double>();
  }
  
  for (map<anabin, distr_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    distr_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    distr_input spp = theVars_inputs[thebinPP];
    
    if (s.naa <= 0 || spp.npp <= 0) continue;
    if ((doprompt || dononprompt) && isData && (spp.bfracpp<=0 || s.bfracaa<=0)) continue;
    if (!isData && s.naa <= 0) continue;
    
    theBins[thebin].push_back(it->first);
    
    double naa = s.naa;
    double npp = spp.npp;
    double dnaa = s.dnaa_stat;
    double dnpp = spp.dnpp_stat;
    
    if (isData)
    {
      if (doprompt) {
        naa = s.naa*(1.-s.bfracaa);
        npp = spp.npp*(1.-spp.bfracpp);
        
        dnaa = naa*sqrt(pow(s.dnaa_stat/s.naa,2)
                        - 2.*s.correlaa*s.dnaa_stat*s.dbfracaa/(s.naa*s.bfracaa)
                        + pow(s.dbfracaa/s.bfracaa,2));
        dnpp = npp*sqrt(pow(spp.dnpp_stat/spp.npp,2)
                        - 2.*spp.correlpp*spp.dnpp_stat*spp.dbfracpp/(spp.npp*spp.bfracpp)
                        + pow(spp.dbfracpp/spp.bfracpp,2));
      }
      if (dononprompt) {
        naa = s.naa*s.bfracaa;
        npp = spp.npp*spp.bfracpp;
        
        dnaa = naa*sqrt(pow(s.dnaa_stat/s.naa,2)
                        + 2.*s.correlaa*s.dnaa_stat*s.dbfracaa/(s.naa*s.bfracaa)
                        + pow(s.dbfracaa/s.bfracaa,2));
        dnpp = npp*sqrt(pow(spp.dnpp_stat/spp.npp,2)
                        + 2.*spp.correlpp*spp.dnpp_stat*spp.dbfracpp/(spp.npp*spp.bfracpp)
                        + pow(spp.dbfracpp/spp.bfracpp,2));
      }
    }
    theVarsBinnedPbPb[thebin].push_back(naa);
    theVarsBinnedPbPb_stat[thebin].push_back(dnaa);
    
    theVarsBinnedPP[thebin].push_back(npp);
    theVarsBinnedPP_stat[thebin].push_back(dnpp);
  }
  
  // make TGraphAsymmErrors
  int cnt=0;
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    int n = theBins[*it].size();
    if(n==0) {
      cout << "Error, nothing found for category" << endl;
      theGraphsPbPb[*it] = NULL;
      theGraphsPP[*it] = NULL;
      continue;
    }
    
    theGraphsPbPb[*it] = new TGraphAsymmErrors(n);
    theGraphsPbPb[*it]->SetName("PbPbDistribution");
    
    theGraphsPP[*it] = new TGraphAsymmErrors(n);
    theGraphsPP[*it]->SetName("PPDistribution");
    
    for (int i=0; i<n; i++) {
      double x=0, exl=0, exh=0, yPbPb=0, eylPbPb=0, eyhPbPb=0;
      double yPP=0, eylPP=0, eyhPP=0;
      double low=0, high=0;
      anabin thebin = theBins[*it][i];
      yPbPb = theVarsBinnedPbPb[*it][i];
      yPP = theVarsBinnedPP[*it][i];
      if (xaxis=="pt" || xaxis=="rap") {
        if (xaxis=="pt") {
          low= thebin.ptbin().low();
          high = thebin.ptbin().high();
        } else {
          low= thebin.rapbin().low();
          high = thebin.rapbin().high();
        }
        x = (low+high)/2.;
        exh = (high-low)/2.;
        exl = (high-low)/2.;
      }
      if (xaxis=="cent") {
        low= thebin.centbin().low();
        high = thebin.centbin().high();

        x = (low+high)/2./2.;
        exh = (high-low)/2./2.;
        exl = (high-low)/2./2.;
      }
      eylPbPb = fabs(theVarsBinnedPbPb_stat[*it][i]);
      eyhPbPb = eylPbPb;
      
      eylPP = fabs(theVarsBinnedPP_stat[*it][i]);
      eyhPP = eylPP;
      
      theGraphsPbPb[*it]->SetPoint(i,x,yPbPb);
      theGraphsPbPb[*it]->SetPointError(i,exl,exh,eylPbPb,eyhPbPb);
      
      theGraphsPP[*it]->SetPoint(i,x,yPP);
      theGraphsPP[*it]->SetPointError(i,exl,exh,eylPP,eyhPP);

    }
    cnt++;
  }
  
  // plot
  plotGraph(theGraphsPbPb, xaxis);
  plotGraph(theGraphsPP, xaxis,true);
}


void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, bool isPP) {
  
  vector<anabin> theCats;
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    if (!tg) {cout << "TGraph not found for bin" << endl; continue;}
    
    theCats.push_back(thebin);
    
    TString centLow(Form("%d",(int)(it->first.centbin().low()/2.)));
    TString centHigh(Form("%d",(int)(it->first.centbin().high()/2.)));
    TString rapLow(Form("%.1f",it->first.rapbin().low())); rapLow.Remove(rapLow.First('.'),1);
    TString rapHigh(Form("%.1f",it->first.rapbin().high()));rapHigh.Remove(rapHigh.First('.'),1);
    TString ptLow(Form("%.1f",it->first.ptbin().low()));ptLow.Remove(ptLow.First('.'),1);
    TString ptHigh(Form("%.1f",it->first.ptbin().high()));ptHigh.Remove(ptHigh.First('.'),1);
    TString distrName(Form("c%s_%s_",isPP ? "PP" : "PbPb", xaxis.c_str()));
    if (xaxis == "cent") distrName += Form("rap%s%s_pt%s%s",rapLow.Data(),rapHigh.Data(),ptLow.Data(),ptHigh.Data());
    else if (xaxis == "rap") distrName += Form("pt%s%s_cent%s%s",ptLow.Data(),ptHigh.Data(),centLow.Data(),centHigh.Data());
    else if (xaxis == "pt") distrName += Form("rap%s%s_cent%s%s",rapLow.Data(),rapHigh.Data(),centLow.Data(),centHigh.Data());
    else {cout << "Unkown axis" << endl; return;};
    
    tg->SetMarkerStyle(markerstyle(cnt));
    tg->SetMarkerColor(color(cnt));
    tg->SetLineColor(color(cnt));
    tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);

    gStyle->SetEndErrorSize(5);

    TString baseDir(Form("Distributions%s/", nameTag.c_str()));
    gSystem->mkdir(baseDir.Data(), kTRUE);
    TFile* fd = new TFile(Form("%s/%s.root", baseDir.Data(), distrName.Data()),"RECREATE");
    tg->Write();
    fd->Close();
    delete fd;
    
    cnt++;
  }
}

int color(int i) {
  if (i==0) return kRed+2;
  else if (i==1) return kBlue+2;
  else if (i==2) return kGreen+2;
  else if (i==3) return kCyan+2;
  else if (i==4) return kMagenta+2;
  else if (i==5) return kOrange+2;
  else return kBlack;
}

int markerstyle(int i) {
  if (i==0) return kFullSquare;
  else if (i==1) return kFullCircle;
  else if (i==2) return kFullStar;
  else if (i==3) return kFullCross;
  else if (i==4) return kOpenSquare;
  else if (i==5) return kOpenCircle;
  else return kOpenStar;
}

void setOptions(bool afiterrors, bool adoprompt, bool adononprompt, string anameTag_base) {
  fiterrors = afiterrors;
  doprompt = adoprompt;
  dononprompt = adononprompt;
  nameTag_base = anameTag_base;
  if (doprompt) nameTag_base += "_prompt";
  if (dononprompt) nameTag_base += "_nonprompt";
}

void printOptions() {
  cout <<
  "fiterrors = " << fiterrors << ", "
  "doprompt = " << doprompt << ", "
  "dononprompt = " << dononprompt << ", "
  "nameTag_base = \"" << nameTag_base << "\"" <<
  endl;
}
