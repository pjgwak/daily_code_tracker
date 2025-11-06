#ifndef plotVars_C
#define plotVars_C

#include "../../Fitter/Macros/Utilities/resultUtils.h"
#include "../../Fitter/Macros/Utilities/texUtils.h"
#include "../../Fitter/Macros/Utilities/bin.h"
#include "../../Fitter/Macros/Utilities/bfrac.h"
#include "../../Fitter/Macros/CMS/CMS_lumi.C"
#include "../../Fitter/Macros/CMS/tdrstyle.C"
#include "syst.h"
#include "../../Fitter/results2tree.C"

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "RooFitResult.h" 
#include "TFitResult.h" 
#include "TLine.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

//////////////
// SETTINGS //
//////////////
#define normalCutsDir "nominal"
#define invCutsDir "nonprompt"
string addTagBase="";

/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

// will plot the dependence of varname in workDirName, as a function of pt, centrality or rapidity. 
// collTag should be PP or PbPb. If collTag="", then both PP and PbPb are plotted.
void plotPt(const char* workDirName, const char* varname, int iplot, const char* collTag="", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);
void plotCent(const char* workDirName, const char* varname, int iplot, const char* collTag="", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);
void plotRap(const char* workDirName, const char* varname, const char* collTag="", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);
void plotAll(const char* workDirName, const char* varname, const char* collTag="", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);

// will plot the dependence of varname as a function to xaxis (=pt, cent or rap) for the file in workDirNames (of the form "dir1,dir2,dir3,...")
void plotFiles(const char* workDirNames, const char* varname, const char* xaxis, float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax, 
      const char* collTag="PP", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);
void plotFilesAll(const char* workDirNames, const char* varname, const char* collTag="PP", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);

// will plot the dependence of varname as a function of the analysis bin (scanning all analysis bins)
// if xaxis=="", all bins are plotted. Otherwise only bins differential in xaxis are shown.
void plotFiles(const char* workDirNames, const char* varname, const char* xaxis="", const char* collTag="PP", bool plotErr=true, const char* DSTag="DATA", bool doSig=false);

// will plot the dependence of varname1 and varname2 (from the same file) as a function of the analysis bin (scanning all analysis bins)
// if xaxis=="", all bins are plotted. Otherwise only bins differential in xaxis are shown.
void plotComparisonVars(const char* workDirName, const char* varname1, const char* varname2, const char* xaxis="", const char* collTag="", bool plotErr=true, const char* DSTag="DATA");

////////////////////////
// OTHER DECLARATIONS //
////////////////////////

TGraphErrors* plotVar(TTree *tr, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true, bool doSig=false);
vector<TGraphErrors*> plotVar(TTree *tr, const char* varname, vector<anabin> theBin, string xaxis, string collTag, bool plotErr=true, bool doSig=false);
TGraphErrors* plotVar(const char* filename, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr=true, bool doSig=false);
void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags, const char* workDirName, string collTag="", const char* basename="", const char* label="");
bool binok(anabin thebin, const char* xaxis);
int color(int i);
int markerstyle(int i);
string addTag="";


////////////////////
// IMPLEMENTATION //
////////////////////

void plotPt(const char* workDirName, const char* varname, int iplot, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
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

   addTag = addTagBase + Form("_%i",iplot);

   TFile *f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
   if (!f || !f->IsOpen()) {
      results2tree(workDirName,DSTag);
      f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg;
   if (string(collTag) != "") tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr, doSig);
   else {
      // plot for pp and pbpb and concatenate the results to tg
     vector<TGraphErrors*> tgpp =  plotVar(tr, varname, theCats, xaxis, "PP", plotErr, doSig);
     vector<TGraphErrors*> tgpbpb =  plotVar(tr, varname, theCats, xaxis, "PbPb", plotErr, doSig);
      tg.insert(tg.end(), tgpp.begin(), tgpp.end());
      tg.insert(tg.end(), tgpbpb.begin(), tgpbpb.end());
   }
   vector<string> tags;
   vector<string> collTags;
   if (string(collTag) != "") collTags.push_back(collTag);
   else {
      collTags.push_back("PP");
      collTags.push_back("PbPb");
   }
   vector<anabin>::const_iterator it;
   vector<string>::const_iterator itc;
   for (itc=collTags.begin(); itc!=collTags.end(); itc++) {
      for (it=theCats.begin(); it!=theCats.end(); it++) {
         ostringstream oss;
         oss.precision(1);
         oss.setf(ios::fixed);
         oss << *itc << ": "
            << it->rapbin().low() << " < |y| < " << it->rapbin().high() << ", " 
            << it->centbin().low()/2. << "%-" << it->centbin().high()/2. << "%";
         tags.push_back(oss.str());
      }
   }
   plotGraphs(tg, tags, workDirName, collTag);
}

void plotCent(const char* workDirName, const char* varname, int iplot, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
   string xaxis = "cent";
   vector<anabin> theCats;

   // 4 rapidity intervals
   if (iplot==0) {
      theCats.push_back(anabin(0,0.6,6.5,50,0,200));
      theCats.push_back(anabin(0.6,1.2,6.5,50,0,200));
      theCats.push_back(anabin(1.2,1.8,6.5,50,0,200));
      theCats.push_back(anabin(1.8,2.4,6.5,50,0,200));
   }

   // 1 rapidity interval
   if (iplot==1) {
      theCats.push_back(anabin(0,2.4,6.5,50,0,200));
   }

   // fwd and low pt
   if (iplot==2) {
      theCats.push_back(anabin(1.8,2.4,3,6.5,0,200));
      theCats.push_back(anabin(1.8,2.4,6.5,50,0,200));
   }

   addTag = addTagBase + Form("_%i",iplot);

   TFile *f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
   if (!f || !f->IsOpen()) {
      results2tree(workDirName,DSTag);
      f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg;
   if (string(collTag) != "") tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr, doSig);
   else {
      // plot for pp and pbpb and concatenate the results to tg
     vector<TGraphErrors*> tgpp =  plotVar(tr, varname, theCats, xaxis, "PP", plotErr, doSig);
     vector<TGraphErrors*> tgpbpb =  plotVar(tr, varname, theCats, xaxis, "PbPb", plotErr, doSig);
      tg.insert(tg.end(), tgpp.begin(), tgpp.end());
      tg.insert(tg.end(), tgpbpb.begin(), tgpbpb.end());
   }
   vector<string> tags;
   vector<string> collTags;
   if (string(collTag) != "") collTags.push_back(collTag);
   else {
      collTags.push_back("PP");
      collTags.push_back("PbPb");
   }
   vector<anabin>::const_iterator it;
   vector<string>::const_iterator itc;
   for (itc=collTags.begin(); itc!=collTags.end(); itc++) {
      for (it=theCats.begin(); it!=theCats.end(); it++) {
         ostringstream oss;
         oss.precision(1);
         oss.setf(ios::fixed);
         oss << *itc << ": "
            << it->rapbin().low() << " < |y| < " << it->rapbin().high() << ", " 
            << it->ptbin().low() << " < p_{T} < " << it->ptbin().high() << " GeV/c";
         tags.push_back(oss.str());
      }
   }
   plotGraphs(tg, tags, workDirName, collTag);
}

void plotRap(const char* workDirName, const char* varname, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
   string xaxis = "rap";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,2.4,6.5,50,0,200));
   addTag = addTagBase;

   TFile *f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
   if (!f || !f->IsOpen()) {
      results2tree(workDirName,DSTag);
      f = new TFile(treeFileName(workDirName,DSTag,"","mass"));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   vector<TGraphErrors*> tg;
   if (string(collTag) != "") tg = plotVar(tr, varname, theCats, xaxis, collTag, plotErr, doSig);
   else {
      // plot for pp and pbpb and concatenate the results to tg
     vector<TGraphErrors*> tgpp =  plotVar(tr, varname, theCats, xaxis, "PP", plotErr, doSig);
     vector<TGraphErrors*> tgpbpb =  plotVar(tr, varname, theCats, xaxis, "PbPb", plotErr, doSig);
      tg.insert(tg.end(), tgpp.begin(), tgpp.end());
      tg.insert(tg.end(), tgpbpb.begin(), tgpbpb.end());
   }
   vector<string> tags;
   vector<string> collTags;
   if (string(collTag) != "") collTags.push_back(collTag);
   else {
      collTags.push_back("PP");
      collTags.push_back("PbPb");
   }
   vector<anabin>::const_iterator it;
   vector<string>::const_iterator itc;
   for (itc=collTags.begin(); itc!=collTags.end(); itc++) {
      for (it=theCats.begin(); it!=theCats.end(); it++) {
         ostringstream oss;
         oss.precision(1);
         oss.setf(ios::fixed);
         oss << *itc << ": "
            << it->ptbin().low() << " < p_{T} < " << it->ptbin().high() << " GeV/c, " 
            << it->centbin().low() << "-" << it->centbin().high() << "\%";
         tags.push_back(oss.str());
      }
   }
   plotGraphs(tg, tags, workDirName, collTag);
}

void plotAll(const char* workDirName, const char* varname, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
   plotPt(workDirName, varname, 0, collTag, plotErr, DSTag, doSig);
   plotPt(workDirName, varname, 1, "PbPb", plotErr, DSTag, doSig);
   plotPt(workDirName, varname, 2, collTag, plotErr, DSTag, doSig);
   plotCent(workDirName, varname, 0, "PbPb", plotErr, DSTag, doSig);
   plotCent(workDirName, varname, 1, "PbPb", plotErr, DSTag, doSig);
   plotCent(workDirName, varname, 2, "PbPb", plotErr, DSTag, doSig);
   plotRap(workDirName, varname, collTag, plotErr, DSTag, doSig);
}

void plotFilesAll(const char* workDirNames, const char* varname, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
   plotFiles(workDirNames, varname, "rap", 0, 2.4, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   plotFiles(workDirNames, varname, "pt", 0, 2.4, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   plotFiles(workDirNames, varname, "pt", 0, 0.6, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   plotFiles(workDirNames, varname, "pt", 0.6, 1.2, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   plotFiles(workDirNames, varname, "pt", 1.2, 1.8, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   plotFiles(workDirNames, varname, "pt", 1.8, 2.4, 3, 50, 0, 200, collTag, plotErr, DSTag, doSig);
   if (TString(collTag)=="PbPb") {
      plotFiles(workDirNames, varname, "pt", 0, 2.4, 6.5, 50, 0, 20, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "pt", 0, 2.4, 6.5, 50, 20, 60, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "pt", 0, 2.4, 6.5, 50, 60, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 0, 2.4, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 0, 0.6, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 0.6, 1.2, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 1.2, 1.8, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 1.8, 2.4, 6.5, 50, 0, 200, collTag, plotErr, DSTag, doSig);
      plotFiles(workDirNames, varname, "cent", 1.8, 2.4, 3, 6.5, 0, 200, collTag, plotErr, DSTag, doSig);
   }
}

void plotFiles(const char* workDirNames, const char* varname, const char* xaxis, float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax, 
               const char* collTag, bool plotErr, const char* DSTag, bool doSig) {

   vector<TGraphErrors*> tg;
   vector<string> tags;
   anabin theBin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);

   TString workDirNamesStr(workDirNames);
   TString workDirName; Int_t from = 0;
   while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
     TGraphErrors *tgg = plotVar(treeFileName(workDirName,DSTag), varname, theBin, xaxis, collTag, plotErr, doSig);
      if (!tgg) {
         results2tree(workDirName,DSTag);
         tgg = plotVar(treeFileName(workDirName,DSTag), varname, theBin, xaxis, collTag, plotErr, doSig);
      }
      if (tgg) {
         tg.push_back(tgg);
         tags.push_back(workDirName.Data());
      }
   }

   if (tg.size()>0) plotGraphs(tg, tags, tags[0].c_str(), collTag, 
         Form("%s_pt%i%i_rap%i%i_cent%i%i", collTag,(int)ptmin*10,(int)ptmax*10,(int)rapmin*10,(int)rapmax*10,centmin,centmax),
         Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f, %i-%i%s",rapmin,rapmax,ptmin,ptmax,centmin/2,centmax/2,"%"));
}

TGraphErrors* plotVar(TTree *tr, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr, bool doSig) {
   if (!tr) return NULL;
   vector<double> x, ex, y, ey;
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   float val, val_err=0;
   int ival=-999;
   float valmax=0, valmin=0;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);

   if (string(varname)=="nll" || string(varname)=="chi2" || string(varname)=="normchi2" || string(varname) == "chi2prob") {
      tr->SetBranchAddress(varname,&val);
   } else if (string(varname).find("npar") != string::npos) {
      tr->SetBranchAddress(varname,&ival);
   } else {
      tr->SetBranchAddress(Form("%s_val",varname),&val);
   }
   if (plotErr) tr->SetBranchAddress(Form("%s_err",varname),&val_err);
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);

      // special case of npar
      if (ival>=0) val=ival;

      anabin trbin(ymin, ymax, ptmin, ptmax, centmin, centmax);
      // special case of PP and centrality
      if (collTag=="PP" /*&& xaxis=="cent"*/) trbin = anabin(ymin, ymax, ptmin, ptmax, 0, 200);
      // general case
      if (!binok(theBin, xaxis, trbin, false)) continue;
      if (string(collSystem) != collTag) continue;

      // special cases of Bfrac and Nprompt
      if (TString(varname).Index("Bfrac") != kNPOS || TString(varname).Index("Prompt") != kNPOS) {
         RooRealVar bfracvar = bfrac(normalCutsDir,invCutsDir,trbin,Form("%s_%s",varname,collTag.c_str()));
         val = bfracvar.getVal();
         val_err = bfracvar.getError();
      }
      if (doSig) {
        val = val / val_err;
        val_err = 0.0;
      }


      if (xaxis=="pt") {
         x.push_back((ptmin+ptmax)/2.);
         ex.push_back((ptmax-ptmin)/2.);
      } else if (xaxis=="cent") {
         x.push_back((centmin+centmax)/4.);
         ex.push_back((centmax-centmin)/4.);
      } else { // if (xaxis=="rap") 
         x.push_back((ymin+ymax)/2.);
         ex.push_back((ymax-ymin)/2.);
      }
      y.push_back(val);
      ey.push_back(val_err);

      // min and max
      valmax = max(valmax, (float) 1.6*(val+val_err));
      valmin = min(valmin, (float) (val>0 ? 0 : 1.1*(val-val_err)));
   }

   int n = x.size();

   TString name = Form("%s_%.1f_%.1f_%.1f_%.1f_%i_%i",collTag.c_str(),
         theBin.rapbin().low(),theBin.rapbin().high(),
         theBin.ptbin().low(),theBin.ptbin().high(),
         theBin.centbin().low(),theBin.centbin().high());
   TString hname = "haxes_" + name;
   TString gname = "gr_" + name;
   TGraphErrors *ans = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
   ans->SetName(gname);
   
   // if plotting vs. centrality, do not plot wide bins
   prune(ans);

   TH1F *haxes=NULL;
   if (xaxis=="pt") {
      haxes = new TH1F(hname,Form(";p_{T} (GeV/c);%s",varname),1,0,50);
   } else if (xaxis=="cent") {
      haxes = new TH1F(hname,Form(";Centrality bin;%s",varname),1,0,100);
   } else { // if (xaxis=="rap")
      haxes = new TH1F(hname,Form(";|y|;%s",varname),1,0,2.4);
   }
   haxes->GetYaxis()->SetLimits(valmin, valmax);
   haxes->GetYaxis()->SetRangeUser(valmin, valmax);
   haxes->GetYaxis()->SetTitleOffset(1.4);
   if (doSig) haxes->GetYaxis()->SetTitle((string(varname)+ " Significance").c_str());
   ans->SetHistogram(haxes);

   ans->Sort();

   return ans;
}

vector<TGraphErrors*> plotVar(TTree *tr, const char* varname, vector<anabin> theBin, string xaxis, string collTag, bool plotErr, bool doSig) {
   vector<TGraphErrors*> ans;

   vector<anabin>::const_iterator it;
   for (it=theBin.begin(); it!=theBin.end(); it++) {
     ans.push_back(plotVar(tr, varname, *it, xaxis, collTag, plotErr, doSig));
   }

   return ans;
}

TGraphErrors* plotVar(const char* filename, const char* varname, anabin theBin, string xaxis, string collTag, bool plotErr, bool doSig) {
   TFile *f = new TFile(filename);
   if (!f) return NULL;
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return NULL;
   return plotVar(tr, varname, theBin, xaxis, collTag, plotErr, doSig);
}

void plotGraphs(vector<TGraphErrors*> graphs, vector<string> tags, const char* workDirName, string collTag, const char* basename, const char* label) {
   if (graphs.size() != tags.size()) {
      cout << "Different number of graphs and legends" << endl;
      return;
   }
   if (graphs.size() == 0) return;

   setTDRStyle();
   TCanvas *c1 = new TCanvas("c1","c1",600,600);

   TLegend *tleg = new TLegend(0.18,0.69,0.52,0.89);
   tleg->SetBorderSize(0);
   tleg->SetTextSize(0.03);

   TH1 *haxes = NULL;
   double ymin=0, ymax=0;

   bool issysts = string(workDirName)=="systematics";

   for (unsigned int i=0; i<graphs.size(); i++) {
      int thecolor = (issysts && tags[i] == "Total") ? kBlack : color(i);
      graphs[i]->SetLineColor(thecolor);
      graphs[i]->SetMarkerColor(thecolor);
      graphs[i]->SetMarkerStyle(markerstyle(i));
      graphs[i]->SetMarkerSize(issysts ? 0 : 1.5);
      if (issysts) {
         graphs[i]->SetLineWidth(2);
         graphs[i]->SetFillStyle(0);
      }
      if (i==0) {
         graphs[i]->Draw(!issysts ? "AP" : "A2");
         haxes = graphs[i]->GetHistogram();
         ymin = haxes->GetYaxis()->GetXmin();
         ymax = haxes->GetYaxis()->GetXmax();
      }
      else {
         graphs[i]->Draw(!issysts ? "P" : "2");
         ymin = min(ymin, graphs[i]->GetYaxis()->GetXmin());
         ymax = max(ymax, graphs[i]->GetYaxis()->GetXmax());
      }
      tleg->AddEntry(graphs[i],tags[i].c_str(),"LP");
   }

   if (issysts) graphs[0]->Draw("2"); // in the case of systs, re-draw the total on top

   if (haxes) haxes->GetYaxis()->SetRangeUser(ymin, ymax);
   // c1->Update();

   tleg->Draw();

   int iPos = 33;
   int ilumi = 106;
   if (collTag=="PP") ilumi = 107;
   if (collTag=="PbPb") ilumi = 108;
   CMS_lumi( (TPad*) gPad, ilumi, iPos, "" );

   string yaxis = haxes->GetYaxis()->GetTitle();
   string xaxis = "rap";
   TString txaxis(haxes->GetXaxis()->GetTitle());
   if (txaxis.Index("Centrality") != kNPOS) xaxis = "cent";
   if (txaxis.Index("p_{T}") != kNPOS) xaxis = "pt";

   TLatex tl;
   tl.SetTextAlign(32); // right align
   tl.DrawLatexNDC(0.93,0.17,label);

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", workDirName), kTRUE); 
   c1->SaveAs(Form("Output/%s/plot/RESULT/root/plot_%s_%s_%s_vs_%s%s.root",
            workDirName, basename, yaxis.c_str(), collTag.c_str(), xaxis.c_str(),addTag.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", workDirName), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/png/plot_%s_%s_%s_vs_%s%s.png",
            workDirName, basename, yaxis.c_str(), collTag.c_str(), xaxis.c_str(),addTag.c_str()));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", workDirName), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/plot_%s_%s_%s_vs_%s%s.pdf",
            workDirName, basename, yaxis.c_str(), collTag.c_str(), xaxis.c_str(),addTag.c_str()));

   // print the result tables
   string xname = "$|y|$";
   if (xaxis=="cent") xname = "Centrality";
   if (xaxis=="pt") xname = "\\pt";
   string yname = latexSafe(yaxis); // make the name safe for LaTeX
   gSystem->mkdir(Form("Output/%s/tex/", workDirName), kTRUE); 
   char texname[2048]; 
   sprintf(texname, "Output/%s/tex/result_%s_%s_vs_%s.tex",workDirName, basename, yaxis.c_str(), xaxis.c_str());
   vector<string> tags_fixed;
   for (vector<string>::const_iterator it=tags.begin(); it!=tags.end(); it++) {
      string tagfixed=latexSafe(*it);
      tags_fixed.push_back(tagfixed);
   }
   bool samesize=true;
   for (unsigned int i=0; i<graphs.size(); i++) {
      if (graphs[i]->GetN() != graphs[0]->GetN()) samesize=false;
   }
   if (samesize) {
      inittex(texname, xname.c_str(), tags_fixed);
      addline(texname,yname,graphs.size());
      printGraph(graphs, texname);
   } else {
      inittex(texname, xname.c_str(), yname);
      for (unsigned int i=0; i<graphs.size(); i++) {
         addline(texname,tags_fixed[i]);
         printGraph(graphs[i], texname);
      }
   }
   closetex(texname);
   cout << "Closed " << texname << endl;
   cout << "It is advised that you check the contents of " << texname << " as it may not compile nicely as is." << endl;
}

void plotFiles(const char* workDirNames, const char* varname, const char* xaxis, const char* collTag, bool plotErr, const char* DSTag, bool doSig) {
   setTDRStyle();

   typedef pair<anabin,string> anabinc;

   vector<string> tags;
   vector<map<anabinc, float> > vals;
   vector<map<anabinc, float> > errs;
   map<anabinc,int> binmap;

   float pull=0, pullmin=0, pullmax=0;
   TTree *tpull = new TTree("tpull","tpull");
   tpull->Branch("pull",&pull,"pull/F");

   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   char collSystem[5];
   float val, val_err=0;

   TString workDirNamesStr(workDirNames);
   TString workDirName; Int_t from = 0;
   int ibin=1;
   while (workDirNamesStr.Tokenize(workDirName, from , ",")) {
      TFile *f = new TFile(treeFileName(workDirName,DSTag));
      if (!f || !f->IsOpen()) {
         results2tree(workDirName,DSTag);
         f = new TFile(treeFileName(workDirName,DSTag));
         if (!f) return;
      }
      TTree *tr = (TTree*) f->Get("fitresults");
      if (!tr) continue;
      tr->SetBranchAddress("ptmin",&ptmin);
      tr->SetBranchAddress("ptmax",&ptmax);
      tr->SetBranchAddress("ymin",&ymin);
      tr->SetBranchAddress("ymax",&ymax);
      tr->SetBranchAddress("centmin",&centmin);
      tr->SetBranchAddress("centmax",&centmax);
      tr->SetBranchAddress("collSystem",collSystem);
      if (string(varname)=="nll" || string(varname)=="chi2" || string(varname)=="normchi2" || string(varname) == "chi2prob") {
         tr->SetBranchAddress(varname,&val);
      } else if (string(varname).find("npar") != string::npos) {
         tr->SetBranchAddress(varname,&val);
      } else {
         tr->SetBranchAddress(Form("%s_val",varname),&val);
      }
      if (plotErr) tr->SetBranchAddress(Form("%s_err",varname),&val_err);

      map<anabinc,float> vals0;
      map<anabinc,float> errs0;

      for (int i=0; i<tr->GetEntries(); i++) {
         tr->GetEntry(i);
         anabinc thebin(anabin(ymin,ymax,ptmin,ptmax,centmin,centmax), collSystem);
         if ((string(collTag) == "" || string(collTag) == collSystem) && binok(thebin.first,xaxis)) {
            if (doSig) {
               val = val/val_err;
               val_err=0;
            }
            vals0[thebin] = val;
            errs0[thebin] = val_err;

            if (binmap.find(thebin) == binmap.end()) {
               binmap[thebin] = ibin;
               ibin++;
            }
         }
      }

      vals.push_back(vals0);
      errs.push_back(errs0);
      tags.push_back(string(workDirName.Data()));

      delete f;
   }

   int nbins = binmap.size();
   TCanvas *c1 = new TCanvas("c1", "",139,267,1265,537);
   c1->Range(-2.920651,-1.349463,34.28077,6.101919);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.07850912);
   c1->SetRightMargin(0.1419508);
   c1->SetTopMargin(0.03543307);
   c1->SetBottomMargin(0.1811024);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   TLegend *tleg = new TLegend(0.18,0.73,0.52,0.89);
   tleg->SetBorderSize(0);

   // for the ratios
   TCanvas *cratio = new TCanvas("cratio", "",139,267,1265,537);
   cratio->Range(-2.920651,-1.349463,34.28077,6.101919);
   cratio->SetFillColor(0);
   cratio->SetBorderMode(0);
   cratio->SetBorderSize(2);
   cratio->SetTickx(1);
   cratio->SetTicky(1);
   cratio->SetLeftMargin(0.07850912);
   cratio->SetRightMargin(0.1419508);
   cratio->SetTopMargin(0.03543307);
   cratio->SetBottomMargin(0.1811024);
   cratio->SetFrameFillStyle(0);
   cratio->SetFrameBorderMode(0);
   cratio->SetFrameFillStyle(0);
   cratio->SetFrameBorderMode(0);
   cratio->SetGridy();
   TH1F *hden = new TH1F("hden",";;ratio",nbins,0,nbins);
   TLegend *tleg_ratio = new TLegend(0.18,0.73,0.52,0.89);
   tleg_ratio->SetBorderSize(0);

   c1->cd();

   map<anabinc,float> vals0 = vals[0];
   map<anabinc,float> errs0 = errs[0];

   double vmin=0,vmax=0;
   TH1F *haxes=NULL;
   for (unsigned int i=0; i<vals.size(); i++) {
      map<anabinc,float> valsi = vals[i];
      map<anabinc,float> errsi = errs[i];
      string tagi = tags[i];
      TH1F *h = new TH1F(Form("h%i",i),Form(";;%s",varname),nbins,0,nbins);
      if (i==0) haxes=h;
      h->SetLineColor((i<4) ? 1+i : 2+i); // skip yellow
      h->SetMarkerColor((i<4) ? 1+i : 2+i);
      h->SetMarkerStyle(20+i);
      h->SetMarkerSize(1.5);

      map<anabinc,float>::const_iterator itv=valsi.begin(),ite=errsi.begin();
      for (;itv!=valsi.end();itv++,ite++) {
         anabinc thebinc = itv->first;
         anabin thebin = thebinc.first;
         h->SetBinContent(binmap[thebinc],itv->second);
         h->SetBinError(binmap[thebinc],ite->second);
         h->GetXaxis()->SetBinLabel(binmap[thebinc],Form("%s-Pt[%.1f-%.1f]-Y[%.1f-%.1f]-C[%i-%i]",thebinc.second.c_str(),
                  thebin.ptbin().low(),thebin.ptbin().high(),
                  thebin.rapbin().low(),thebin.rapbin().high(),
                  thebin.centbin().low(),thebin.centbin().high()));
         // h->GetXaxis()->LabelsOption("v");

         // pulls
         if (i>1 && i%2==0 && plotErr && ite->second != 0) {
            pull = (itv->second-vals[i-1][thebinc]) / sqrt(pow(ite->second,2) + pow(errs[i-1][thebinc],2));
            // hpull->Fill(pull);
            tpull->Fill(); 
            pullmin = min(pullmin,pull); pullmax = max(pullmax,pull);
         }
         if (i==1 && vals.size()==2) {
            pull = (itv->second-vals[0][thebinc]) / sqrt(pow(ite->second,2) + pow(errs[0][thebinc],2));
            // hpull->Fill(pull);
            tpull->Fill(); 
            pullmin = min(pullmin,pull); pullmax = max(pullmax,pull);
         }
      }

      int binmin = h->GetMinimumBin();
      int binmax = h->GetMaximumBin();
      vmin = min(vmin,h->GetBinContent(binmin)-1.2*h->GetBinError(binmin));
      vmax = max(vmax,h->GetBinContent(binmax)+1.2*h->GetBinError(binmax));

      c1->cd();
      h->Draw(i==0 ? "" : "same");
      tleg->AddEntry(h,tagi.c_str(),"l");

      // ratios
      if (i==0) {
         hden = h;
      } else {
         cratio->cd();
         TH1F *hratio = (TH1F*) h->Clone(Form("%s_ratio",h->GetName()));
         if (doSig) hratio->GetYaxis()->SetTitle((string(varname)+ "_Significance").c_str());
         hratio->GetYaxis()->SetTitleOffset(0.5);
         hratio->GetXaxis()->SetLabelSize(0.04);
         hratio->Divide(hden);
         hratio->Draw(i == 1 ? "" : "same");
         tleg_ratio->AddEntry(hratio,Form("%s / %s",tagi.c_str(),tags[0].c_str()),"l");
         c1->cd();
      }
   }

   if (haxes) haxes->GetYaxis()->SetRangeUser(vmin,vmax);
   c1->cd();
   tleg->Draw();
   int iPos = 33;
   int ilumi = 106;
   if (string(collTag)=="PP") ilumi = 107;
   if (string(collTag)=="PbPb") ilumi = 108;
   CMS_lumi( (TPad*) gPad, ilumi, iPos, "" );

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", tags[0].c_str()), kTRUE); 
   c1->SaveAs(Form("Output/%s/plot/RESULT/root/plot_%s_%s_vs_%s.root",tags[0].c_str(), collTag, varname, xaxis));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", tags[0].c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/png/plot_%s_%s_vs_%s.png",tags[0].c_str(), collTag, varname, xaxis));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", tags[0].c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/plot_%s_%s_vs_%s.pdf",tags[0].c_str(), collTag, varname, xaxis));

   // ratios
   cratio->cd();
   tleg_ratio->Draw();
   CMS_lumi( (TPad*) gPad, ilumi, iPos, "" );

   cratio->cd();
   cratio->Update();
   cratio->RedrawAxis();
   cratio->SaveAs(Form("Output/%s/plot/RESULT/root/ratio_%s_%s_vs_%s.root",tags[0].c_str(), collTag, varname, xaxis));
   cratio->SaveAs(Form("Output/%s/plot/RESULT/png/ratio_%s_%s_vs_%s.png",tags[0].c_str(), collTag, varname, xaxis));
   cratio->SaveAs(Form("Output/%s/plot/RESULT/pdf/ratio_%s_%s_vs_%s.pdf",tags[0].c_str(), collTag, varname, xaxis));

   // pulls
   TCanvas *cpull = new TCanvas("cpull","cpull",600,600); cpull->cd();

   // find plotting range
   float pullmaxplot = ((int) (max(fabs(pullmin),fabs(pullmax)))) + 1;
   TH1F *hpull = new TH1F("hpull","Pull distribution;Pull;Entries",10,-pullmaxplot,pullmaxplot);
   hpull->Sumw2(kFALSE); hpull->SetBinErrorOption(TH1::kPoisson);
   for (int i=0; i<tpull->GetEntries(); i++) {
   tpull->GetEntry(i);
   hpull->Fill(pull);
   }

   hpull->Draw("E");
   auto r = hpull->Fit("gaus","LES");
   double chi2_BC =  2.* r->MinFcnValue(); 
   cout << "chi2 (Baker-Cousins) = " << chi2_BC << endl;
   CMS_lumi( (TPad*) gPad, ilumi, iPos, "" );

   cpull->cd();
   cpull->Update();
   cpull->RedrawAxis();
   cpull->SaveAs(Form("Output/%s/plot/RESULT/root/pull_%s_%s_vs_%s.root",tags[0].c_str(), collTag, varname, xaxis));
   cpull->SaveAs(Form("Output/%s/plot/RESULT/png/pull_%s_%s_vs_%s.png",tags[0].c_str(), collTag, varname, xaxis));
   cpull->SaveAs(Form("Output/%s/plot/RESULT/pdf/pull_%s_%s_vs_%s.pdf",tags[0].c_str(), collTag, varname, xaxis));

}

void plotComparisonVars(const char* workDirName, const char* varname1, const char* varname2, const char* xaxis, const char* collTag, bool plotErr, const char* DSTag) {
  setTDRStyle();
  
  typedef pair<anabin,string> anabinc;
  
  map<anabinc,int> binmap;
  
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  char collSystem[5];
  float val1, val_err1=0;
  float val2, val_err2=0;
  
  int ibin=1;
  TFile *f = new TFile(treeFileName(workDirName,DSTag));
  if (!f || !f->IsOpen()) {
    results2tree(workDirName,DSTag);
    f = new TFile(treeFileName(workDirName,DSTag));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr)
  {
    cout << "[Error]: Results tree could not be retrieved" << endl;
    return;
  }
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("collSystem",collSystem);
  if (string(varname1)=="nll" || string(varname1)=="chi2" || string(varname1)=="normchi2" || string(varname1) == "chi2prob") {
    tr->SetBranchAddress(varname1,&val1);
    tr->SetBranchAddress(varname2,&val2);
  } else if (string(varname1).find("npar") != string::npos) {
    tr->SetBranchAddress(varname1,&val1);
    tr->SetBranchAddress(varname2,&val2);
  } else {
    tr->SetBranchAddress(Form("%s_val",varname1),&val1);
    tr->SetBranchAddress(Form("%s_val",varname2),&val2);
  }
  if (plotErr)
  {
    tr->SetBranchAddress(Form("%s_err",varname1),&val_err1);
    tr->SetBranchAddress(Form("%s_err",varname2),&val_err2);
  }
  
  map<anabinc,float> vals1;
  map<anabinc,float> errs1;
  
  map<anabinc,float> vals2;
  map<anabinc,float> errs2;
  
  for (int i=0; i<tr->GetEntries(); i++) {
    tr->GetEntry(i);
    anabinc thebin(anabin(ymin,ymax,ptmin,ptmax,centmin,centmax), collSystem);
    if ((string(collTag) == "" || string(collTag) == collSystem) && binok(thebin.first,xaxis)) {
      vals1[thebin] = val1;
      errs1[thebin] = val_err1;
      
      vals2[thebin] = val2;
      errs2[thebin] = val_err2;
      
      if (binmap.find(thebin) == binmap.end()) {
        binmap[thebin] = ibin;
        ibin++;
      }
    }
  }
  
  delete f;
  
  int nbins = binmap.size();
  TCanvas *c1 = new TCanvas("c1", "",139,267,1265,537);
  c1->Range(-2.920651,-1.349463,34.28077,6.101919);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.07850912);
  c1->SetRightMargin(0.1419508);
  c1->SetTopMargin(0.03543307);
  c1->SetBottomMargin(0.1811024);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  
//  TLegend *tleg = new TLegend(0.18,0.73,0.52,0.89);
//  tleg->SetBorderSize(0);
  
//  TH1F *h1 = new TH1F("h1",Form(";;%s",varname1),nbins,0,nbins);
  TH1F *h1 = new TH1F("h1",Form(";;%s/%s",varname1,varname2),nbins,0,nbins);
  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.5);
  h1->GetYaxis()->SetTitleOffset(0.5);
  h1->GetXaxis()->SetLabelSize(0.04);
  
  TH1F *h2 = new TH1F("h2",Form(";;%s",varname2),nbins,0,nbins);
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(1.5);
  
  
  map<anabinc,float>::const_iterator itv1=vals1.begin(),ite1=errs1.begin();
  for (;itv1!=vals1.end();itv1++,ite1++) {
    anabinc thebinc = itv1->first;
    anabin thebin = thebinc.first;
    h1->SetBinContent(binmap[thebinc],itv1->second);
    h1->SetBinError(binmap[thebinc],ite1->second);
    h1->GetXaxis()->SetBinLabel(binmap[thebinc],Form("%s-Pt[%.1f-%.1f]-Y[%.1f-%.1f]-C[%i-%i]",thebinc.second.c_str(),
                                                     thebin.ptbin().low(),thebin.ptbin().high(),
                                                     thebin.rapbin().low(),thebin.rapbin().high(),
                                                     thebin.centbin().low(),thebin.centbin().high()));
  }
  
  map<anabinc,float>::const_iterator itv2=vals2.begin(),ite2=errs2.begin();
  for (;itv2!=vals2.end();itv2++,ite2++) {
    anabinc thebinc = itv2->first;
    anabin thebin = thebinc.first;
    h2->SetBinContent(binmap[thebinc],itv2->second);
//    h2->SetBinError(binmap[thebinc],ite2->second);
    h2->SetBinError(binmap[thebinc],0.);
    h2->GetXaxis()->SetBinLabel(binmap[thebinc],Form("%s-Pt[%.1f-%.1f]-Y[%.1f-%.1f]-C[%i-%i]",thebinc.second.c_str(),
                                                     thebin.ptbin().low(),thebin.ptbin().high(),
                                                     thebin.rapbin().low(),thebin.rapbin().high(),
                                                     thebin.centbin().low(),thebin.centbin().high()));
  }
  
  
  h1->Divide(h2);
  
  double vmin = h1->GetMinimum();
  double vmax = h1->GetMaximum();
  h1->GetYaxis()->SetRangeUser(0.,vmax+vmax*0.1);
  
//  tleg->AddEntry(h1,varname1,"l");
//  tleg->AddEntry(h2,varname2,"l");
  TLine* l1 = new TLine(0.,1.,h1->GetXaxis()->GetXmax(),1.);
  l1->SetLineWidth(3);
  
  h1->GetListOfFunctions()->Add(l1);
  
  h1->Draw();
//  h2->Draw("same");
//  tleg->Draw("same");
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", workDirName), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/comparison_%s_%s.root",workDirName, varname1, varname2));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", workDirName), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/comparison_%s_%s.pdf",workDirName, varname1, varname2));
}

bool binok(anabin thebin, const char* xaxis) {
   if (string(xaxis) == "") return true;

   if (string(xaxis) == "rap") {
      if (thebin.ptbin()==binF(6.5,50)&&thebin.centbin()==binI(0,200)) return true;
      return false;
   }

   // if the xaxis is not "" not "rap", discard pt- and centrality-integrated bins
   if (thebin.centbin()==binI(0,200) && (thebin.ptbin()==binF(6.5,50) || thebin.ptbin()==binF(3,50))) return false;

   if (string(xaxis) == "pt" && thebin.centbin() == binI(0,200)) return true; 

   if (string(xaxis) == "cent" && thebin.rapbin() == binF(0,2.4) && thebin.ptbin() == binF(6.5,50)) return true; 
   if (string(xaxis) == "cent" && thebin.rapbin() == binF(0,0.6) && thebin.ptbin() == binF(6.5,50)) return true; 
   if (string(xaxis) == "cent" && thebin.rapbin() == binF(0.6,1.2) && thebin.ptbin() == binF(6.5,50)) return true; 
   if (string(xaxis) == "cent" && thebin.rapbin() == binF(1.2,1.8) && thebin.ptbin() == binF(6.5,50)) return true; 
   if (string(xaxis) == "cent" && thebin.rapbin() == binF(1.8,2.4) && thebin.ptbin() == binF(6.5,50)) return true; 
   if (string(xaxis) == "cent" && thebin.rapbin() == binF(1.8,2.4) && thebin.ptbin() == binF(3,6.5)) return true; 

   return false;
}

int color(int i) {
   if (i==0) return kRed+2;
   else if (i==1) return kBlue+2;
   else if (i==2) return kGreen+2;
   else if (i==3) return kCyan+2;
   else if (i==4) return kMagenta+2;
   else if (i==5) return kOrange+2;
   else if (i==6) return kViolet+2;
   else if (i==7) return kTeal+2;
   else if (i==8) return kSpring+2;
   else if (i==9) return kYellow+2;
   else if (i==10) return kPink+2;
   else if (i==11) return kAzure+2;
   else return color(i%12)+(i/12);
}

int markerstyle(int i) {
   if (i==0) return kFullSquare;
   else if (i==1) return kFullCircle;
   else if (i==2) return kFullStar;
   else if (i==3) return kFullCross;
   else if (i==4) return kOpenSquare;
   else if (i==5) return kOpenCircle;
   else if (i==6) return kOpenStar;
   else if (i==7) return kFullTriangleUp;
   else if (i==8) return kFullTriangleDown;
   else if (i==9) return kOpenTriangleUp;
   else if (i==10) return kOpenTriangleDown;
   else return markerstyle(i%11);
}

#endif // #define plotvars_C
