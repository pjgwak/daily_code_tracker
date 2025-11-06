#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/bfrac.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/monster.h"
#include "Macros/Utilities/rappGraph.h"
#include "Systematics/syst.h"
#include "../Limits/limits.h"

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
const char* poiname = "RFrac2Svs1S"; // for double ratios
// const char* poiname = "N_Jpsi"; // for RAA (will correct automatically for efficiency)
// const char* poiname = "N_Psi2S"; // for RAA (will correct automatically for efficiency)
#endif
const char* ylabel = "(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}";
const bool  plot12007_mid = false; // plot 12-007, midrapidity
const bool  plot12007_fwd = false; // plot 12-007, fwdrapidity
const bool  fiterrors     = true;  // statistical errors are from the fit
const bool  FCerrors      = false; // statistical errors are from the Feldman-Cousins intervals ("limits")
const bool  promptonly    = true;  // plot the prompt only double ratio
const bool  nonpromptonly = false; // plot the non-prompt only double ratio
const bool  plotlimits95  = true;  // display 95% CL limits (when the lower limit is 0)
const bool  plotsysts     = true;  // display systematics
const bool  plotrapp      = false; // plot Rapp and Du's predictions
const char* nameTag="";            // can put here e.g. "_prompt", "_nonprompt", ...

const bool plot12007 = plot12007_mid || plot12007_fwd;

///////////////
// CONSTANTS //
///////////////

const double xfrac = 0.8; // for the 2-panel plot for the centrality dependence.
const double lumipp = 27.7e6;
const double lumipbpb_ABCD = 351;
const double lumipbpb_peri = 464;

//////////////////
// DECLARATIONS //
//////////////////

RooRealVar* poiFromFile(const char* filename, const char* token="");
// plot
void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst);
void plot(vector<anabin> thecats, string xaxis, string workDirName, string workDirNameFail="");
void centrality2npart(TGraphAsymmErrors* tg, bool issyst=false, bool isMB=false, double xshift=0.);
void plotLimits(vector<anabin> theCats, string xaxis, const char* filename="../Limits/csv/Limits_95.csv", double xshift=0, bool ULonly=true, bool isInclusive=false);
void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color);



/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

void plotPt(string workDirName, string workDirNameFail="") {
   string xaxis = "pt";
   vector<anabin> theCats;
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   plot(theCats,xaxis,workDirName, workDirNameFail);
};

void plotCent(string workDirName, string workDirNameFail="") {
   string xaxis = "cent";
   vector<anabin> theCats;

   // centrality dependence
   theCats.push_back(anabin(0,1.6,6.5,30,0,200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,200));

   // for minimum bias
   theCats.push_back(anabin(0,1.6,6.5,30,0,-200));
   theCats.push_back(anabin(1.6,2.4,3,30,0,-200));

   plot(theCats,xaxis,workDirName, workDirNameFail);
};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plot(vector<anabin> thecats, string xaxis, string outputDir, string outputDirFail) {
   // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
   // xaxis is the variable to be plotted. "pt", "rap" or "cent"

   TFile *f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
   if (!f || !f->IsOpen()) {
      results2tree(outputDir.c_str(),"DATA");
      f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;

   map<anabin, double> theVars_val;
   map<anabin, double> theVars_stat;
   map<anabin, double> theVars_syst;
   map<anabin,syst> stat_PP; // statistical uncertainty on PP fits (useful for the case of the centrality dependence)

   vector<double> x, ex, y, ey;
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   float val, err=0;
   int ival=-999;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);
   tr->SetBranchAddress(Form("%s_val",poiname),&val);
   tr->SetBranchAddress(Form("%s_err",poiname),&err);
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);
      anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);

      if (TString(collSystem)=="PP") continue;

      double dR, dR_stat, dR_statPP, dR_syst;
      const char* opt = (xaxis=="cent") ? "PbPb" : ""; // for the centrality dependence, exclude PbPb errors
      if (!promptonly && !nonpromptonly) {
         dR = doubleratio_pass_nominal(outputDir.c_str(), thebin);
         dR_stat = doubleratio_pass_stat(outputDir.c_str(), thebin, "", opt);
         dR_statPP = doubleratio_pass_stat(outputDir.c_str(), thebin, "", "PP");
         dR_syst = doubleratio_pass_syst(outputDir.c_str(), thebin, "", opt);
      } else if (promptonly) {
         dR = doubleratio_prompt_nominal(outputDir.c_str(), outputDirFail.c_str(), thebin);
         dR_stat = doubleratio_prompt_stat(outputDir.c_str(), outputDirFail.c_str(), thebin, "", opt);
         dR_statPP = doubleratio_prompt_stat(outputDir.c_str(), outputDirFail.c_str(), thebin, "", "PP");
         dR_syst = doubleratio_prompt_syst(outputDir.c_str(), outputDirFail.c_str(), thebin, "", opt);
      } else {
         dR = doubleratio_nonprompt_nominal(outputDir.c_str(), outputDirFail.c_str(), thebin);
         dR_stat = doubleratio_nonprompt_stat(outputDir.c_str(), outputDirFail.c_str(), thebin, "", opt);
         dR_statPP = doubleratio_nonprompt_stat(outputDir.c_str(), outputDirFail.c_str(), thebin, "", "PP");
         dR_syst = doubleratio_nonprompt_syst(outputDir.c_str(), outputDirFail.c_str(), thebin, "", opt);
      }

      theVars_val[thebin] = dR;
      theVars_stat[thebin] = dR_stat;
      theVars_syst[thebin] = dR_syst;

      syst thestat_PP;
      thestat_PP.name = "stat_PP";
      thestat_PP.value_dR = dR_statPP;
      thestat_PP.value = thestat_PP.value_dR / dR; 
      thestat_PP.value_dR_rel = thestat_PP.value; 
      stat_PP[thebin] = thestat_PP;
   }

   map<anabin, vector<anabin> > theBins;
   map<anabin, vector<double> > theVarsBinned;
   map<anabin, vector<double> > theVarsBinned_stat;
   map<anabin, vector<double> > theVarsBinned_syst;
   map<anabin, TGraphAsymmErrors* > theGraphs;
   map<anabin, TGraphAsymmErrors* > theGraphs_syst;

   // initialize the maps
   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      theBins[*it] = vector<anabin>();
      theVarsBinned[*it] = vector<double>();
   }

   for (map<anabin, double>::const_iterator it=theVars_val.begin(); it!=theVars_val.end(); it++) {
      anabin thebin = it->first;
      if (!binok(thecats,xaxis,thebin)) continue;
      theBins[thebin].push_back(it->first);
      theVarsBinned[thebin].push_back(it->second);
      theVarsBinned_stat[thebin].push_back(theVars_stat[it->first]);
      theVarsBinned_syst[thebin].push_back(theVars_stat[it->first]);
   }

   // systematics
   map<anabin, syst> syst_PP; 
   if (plotsysts) {
      if (!promptonly && !nonpromptonly) syst_PP = readSyst_all_pass("PP","",outputDir.c_str());
      else if (promptonly) syst_PP = readSyst_all_prompt("PP","",outputDir.c_str(),outputDirFail.c_str());
      else syst_PP = readSyst_all_nonprompt("PP","",outputDir.c_str(),outputDirFail.c_str());
   }
   if (xaxis=="cent") { // put the PP stat error into the PP syst, that will go into a box at 1
      vector< map<anabin, syst> > all_PP;
      all_PP.push_back(syst_PP);
      all_PP.push_back(stat_PP);
      syst_PP = combineSyst(all_PP,"statsyst_PP");
   }
   map<anabin, syst> syst_PbPb;
   if (plotsysts) {
      if (!promptonly && !nonpromptonly) syst_PbPb = readSyst_all_pass("PbPb","",outputDir.c_str());
      else if (promptonly) syst_PbPb = readSyst_all_prompt("PbPb","",outputDir.c_str(),outputDirFail.c_str());
      else syst_PbPb = readSyst_all_nonprompt("PbPb","",outputDir.c_str(),outputDirFail.c_str());
   }
   // map<anabin, syst> syst_PbPb_NP_add = readSyst("Systematics/csv/syst_PbPb_bhad_add.csv");

   // make TGraphAsymmErrors
   int cnt=0;
   for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
      int n = theBins[*it].size();
      if(n==0) {
         cout << "Error, nothing found for category" << endl;
         theGraphs[*it] = NULL;
         continue;
      }
      bool isMB = (it->centbin().low()<=0 && it->centbin().high()<=0);
      if (isMB && n!=1) cout << "Warning, I have " << n << " bins for a MinBias category, I expected 1." << endl;

      theGraphs[*it] = new TGraphAsymmErrors(n);
      theGraphs[*it]->SetName(Form("bin_%i",cnt));
      theGraphs_syst[*it] = new TGraphAsymmErrors(n);
      theGraphs_syst[*it]->SetName(Form("bin_%i_syst",cnt));

      for (int i=0; i<n; i++) {
         double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
         double exsyst=0, eysyst=0;
         double low=0, high=0; 
         anabin thebin = theBins[*it][i];
         if (xaxis=="pt") {
            low= thebin.ptbin().low();
            high = thebin.ptbin().high();
            x = (low+high)/2.;
            exh = (high-low)/2.;
            exl = (high-low)/2.;
            exsyst = 0.5;
            eysyst = plotsysts ? sqrt(pow(syst_PP[thebin].value_dR,2) + pow(syst_PbPb[thebin].value_dR,2)) : 0; // quadratic sum of PP and PbPb systs
         }
         if (xaxis=="cent") {
            low= thebin.centbin().low();
            high = thebin.centbin().high();
            // x = isMB ? 150 + (150./1.6)*it->rapbin().low() : HI::findNpartAverage(low,high);
            // exl = 0.;
            // exh = 0.;
            x = (low+high)/2./2.;
            exh = (high-low)/2./2.;
            exl = (high-low)/2./2.;
            // exsyst = !isMB ? 5 : 5./(1.-xfrac);
            exsyst = exl;
            eysyst = plotsysts ? syst_PbPb[thebin].value_dR : 0; // only PbPb syst: the PP one will go to a dedicated box
            if (plotsysts && isMB) eysyst = sqrt(pow(syst_PP[thebin].value_dR,2) + pow(syst_PbPb[thebin].value_dR,2));
         }
         y = theVarsBinned[*it][i];
         if (fiterrors) {
            eyl = fabs(theVarsBinned_stat[*it][i]);
            if (isMB) eyl = sqrt(pow(theVarsBinned_stat[*it][i],2) + pow(stat_PP[*it].value_dR,2));
            eyh = eyl;
         } else {
            map<anabin, limits> maplim = readLimits("../Limits/csv/Limits_68.csv");
            if (maplim.find(*it) != maplim.end()) {
               limits lim = maplim[*it];
               eyl = lim.val.first;
               eyh = lim.val.second;
            } else {
               eyl = fabs(theVarsBinned_stat[*it][i]);
               eyh = theVarsBinned_stat[*it][i];
            }
         }

         // eysyst = y*eysyst;
         // add the additive part of the NP contamination syst
         // if (plotsysts && syst_PbPb_NP_add.find(thebin) != syst_PbPb_NP_add.end()) eysyst = sqrt(pow(syst_PbPb_NP_add[thebin].value,2) + pow(eysyst,2)); 

         theGraphs[*it]->SetPoint(i,x,y);
         theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
         theGraphs_syst[*it]->SetPoint(i,x,y);
         // theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst,eysyst);
         theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst,eysyst);
         // cout << x << " " << y << " " << eyl << " " << eyh << " " << eysyst << endl;

         // theGraphs[*it]->Sort();
         // theGraphs_syst[*it]->Sort();
      }
      cnt++;
   }

   // plot
   plotGraph(theGraphs, theGraphs_syst, xaxis, outputDir, syst_PP);
}

RooRealVar* poiFromFile(const char* filename, const char* token) {
   RooRealVar *ans = poiFromFile(filename,token,poiname);
   
   // case of N_Jpsi and N_Psi2S: we want the correct normalization because we'll make a RAA
   bool isjpsi = (TString(poiname) == "N_Jpsi");
   bool ispsip = (TString(poiname) == "N_Psi2S");
   if (isjpsi || ispsip) {
      bool isPP = (TString(token).Index("PP") != kNPOS);
      double normfactor = 1.;
      anabin thebin = binFromFile(filename);

      // luminosity and Ncoll
      if (isPP) {
         normfactor = 1./lumipp;
      } else {
         if (thebin.centbin().low()>=60) normfactor = 1./lumipbpb_peri;
         else normfactor = 1./lumipbpb_ABCD;
         normfactor *= 1./(208.*208.*(HI::findNcollAverage(thebin.centbin().low(),thebin.centbin().high())/HI::findNcollAverage(0,200)));
         normfactor *= 200./(thebin.centbin().high()-thebin.centbin().low());
      }

      // efficiency
      TFile *f = TFile::Open(Form("../Efficiency/files/histos_%s_%s.root", isjpsi ? "jpsi" : "psi2s", isPP ? "pp" : "pbpb"));
      bool ismid = (thebin.rapbin().low() < 1.);
      bool isptdep = (thebin.centbin() == binI(0,200));
      TH1F *hnum = (TH1F*) f->Get(Form("hnumptdepcut_%s%s", isptdep ? "pt" : "cent", ismid ? "mid" : "fwd"));
      TH1F *hden = (TH1F*) f->Get(Form("hden_%s%s", isptdep ? "pt" : "cent", ismid ? "mid" : "fwd"));
      double numval, numerr, denval, denerr;
      bool isintegrated = (isptdep && ((thebin.rapbin()==binF(0,1.6)&&thebin.ptbin()==binF(6.5,30)) || (thebin.rapbin()==binF(1.6,2.4)&&thebin.ptbin()==binF(3,30))));
      if (isintegrated) {
         numval = hnum->IntegralAndError(1,hnum->GetNbinsX(),numerr);
         denval = hden->IntegralAndError(1,hden->GetNbinsX(),denerr);
      } else {
         int ibin = hnum->FindBin((thebin.centbin().low()+thebin.centbin().high())/4.);
         if (isptdep) ibin = hnum->FindBin((thebin.ptbin().low()+thebin.ptbin().high())/2.);
         numval = hnum->GetBinContent(ibin);
         numerr = hnum->GetBinError(ibin);
         denval = hden->GetBinContent(ibin);
         denerr = hden->GetBinError(ibin);
      }
      double efficiency = numval / denval;
      normfactor = normfactor / efficiency;
      delete f;

      ans->setVal(ans->getVal()*normfactor);
      ans->setError(ans->getError()*normfactor);
   }

   return ans;
}


void plotGraph(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst) {
   setTDRStyle();

   vector<anabin> theCats;

   TCanvas *c1 = NULL;
   if (xaxis=="cent") c1 = new TCanvas("c1","c1",600/xfrac,600);
   else c1 = new TCanvas("c1","c1",600,600);

   // in the case of the centrality dependence, we need the minimum bias panel on the right
   TPad *padl=NULL, *padr=NULL;
   if (xaxis=="cent") {
         padl = new TPad("padl_left","",0,0,xfrac,1);
         padl->SetLeftMargin(gStyle->GetPadLeftMargin()/xfrac);
         padl->SetRightMargin(0);
         // padl->SetBottomMargin(gStyle->GetPadBottomMargin());
         // padl->SetTopMargin(gStyle->GetPadTopMargin());
         padl->SetFrameBorderMode(0);
         padl->SetBorderMode(0);
         padl->SetBorderSize(0);
         padl->Draw();

         padr = new TPad("padr_right","",xfrac,0,1,1);
         padr->SetRightMargin(gStyle->GetPadRightMargin()/(1.-xfrac));
         padr->SetLeftMargin(0);
         // padr->SetBottomMargin(gStyle->GetPadBottomMargin());
         // padr->SetTopMargin(gStyle->GetPadTopMargin());
         padr->SetFrameBorderMode(0);
         padr->SetBorderMode(0);
         padr->SetBorderSize(0);
         padr->Draw();

         padl->cd();
   }

   // the axes
   TH1F *haxes=NULL, *haxesr=NULL; TLine line, liner;
   if (xaxis=="pt") {
      haxes = new TH1F("haxes","haxes",1,0,30);
      line = TLine(0,1,30,1);
   }
   if (xaxis=="cent") {
      haxes = new TH1F("haxesl","haxesl",1,0,420);
      haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X")/xfrac);
      // haxes->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")*xfrac);
      line = TLine(0,1,420,1);
      haxesr = new TH1F("haxesr","haxesr",1,0,420);
      haxesr->GetXaxis()->SetTickLength(0);
      haxesr->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")/(1.-xfrac));
      haxesr->GetYaxis()->SetRangeUser(0,1.5);
      if (nonpromptonly) haxesr->GetYaxis()->SetRangeUser(0,2.5);
      if (plot12007) haxesr->GetYaxis()->SetRangeUser(0,3.2);
      haxesr->GetXaxis()->SetTitleSize(0);
      haxesr->GetXaxis()->SetLabelSize(0);
      liner = TLine(0,1,420,1);
   }
   haxes->GetYaxis()->SetRangeUser(0,1.5);
   if (nonpromptonly) haxes->GetYaxis()->SetRangeUser(0,2.5);
   if (plot12007) haxes->GetYaxis()->SetRangeUser(0,3.2);
   haxes->GetYaxis()->SetTitle(ylabel);
   const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "N_{part}";
   haxes->GetXaxis()->SetTitle(xlabel);
   haxes->Draw();
   line.Draw();
   if (xaxis=="cent") {
      padr->cd();
      haxesr->Draw();
      liner.Draw();
      padl->cd();
   }

   double xshift=0.;
   if (xaxis=="cent") xshift=0.05;
   TLegend *tleg = new TLegend(0.16+xshift,0.73,0.50,0.89);
   tleg->SetBorderSize(0);
   tleg->SetTextSize(0.03);

   // plot Rapp's theory first if we have to
   if (plotrapp) {
      if (plot12007_mid) {
         TGraphErrors *gmid_12007 = rapp_cent_mid_276();
         gmid_12007->SetFillColorAlpha(kCyan-7,0.5); 
         gmid_12007->Draw("3");
      }
      if (plot12007_fwd) {
         TGraphErrors *gfwd_12007 = rapp_cent_fwd_276();
         gfwd_12007->SetFillColorAlpha(kMagenta-7,0.5); 
         gfwd_12007->Draw("3");
      }

      TGraphErrors *gmid=NULL, *gfwd=NULL;
      if (xaxis=="cent") {
         gmid = rapp_cent_mid_502();
         gfwd = rapp_cent_fwd_502();
      }
      if (xaxis=="pt") {
         gmid = rapp_pt_mid_502();
         gfwd = rapp_pt_fwd_502();
      }
      if (gmid && !(plot12007_fwd && !plot12007_mid)) { // if asking for 2.76TeV fwd only, do not plot 5TeV mid
         gmid->SetFillColorAlpha(kBlue-7,0.5); 
         // gmid->SetFillStyle(3440);
         gmid->Draw("3");
      }
      if (gfwd && !(plot12007_mid && !plot12007_fwd)) { // if asking for 2.76TeV mid only, do not plot 5TeV fwd
         gfwd->SetFillColorAlpha(kRed-7,0.5);
         // gfwd->SetFillStyle(3404);
         gfwd->Draw("3");
      }
   }

   // prepare for the printing of the result tables
   const char* xname = (xaxis=="cent") ? "Centrality" : "\\pt";
   gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE); 
   char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s%s.tex",outputDir.c_str(),xaxis.c_str(),nameTag);
   string yname("\\doubleRatio");
   inittex(texname, xname, yname);

   int cnt=0;
   map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
   map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
   for (; it!=theGraphs.end(); it++) {
      anabin thebin = it->first;
      bool isMB = (thebin.centbin().low()<=0 && thebin.centbin().high()<=0);
      TGraphAsymmErrors* tg = it->second;
      TGraphAsymmErrors* tg_syst = it_syst->second;
      if (!tg || !tg_syst) continue;

      theCats.push_back(thebin);

      if (thebin.rapbin() == binF(0.,1.6)) {
         tg->SetMarkerStyle(kFullSquare);
         tg->SetMarkerColor(kBlue);
         tg->SetLineColor(kBlue);
         tg_syst->SetFillColorAlpha(kBlue, 0.5);
      } else if (thebin.rapbin() == binF(1.6,2.4)) {
         tg->SetMarkerStyle(kFullCircle);
         tg->SetMarkerColor(kRed);
         tg->SetLineColor(kRed);
         tg_syst->SetFillColorAlpha(kRed, 0.5);
      } else {
         tg->SetMarkerStyle(kFullTriangleUp);
         tg->SetMarkerColor(kGreen);
         tg->SetLineColor(kGreen);
         tg_syst->SetFillColorAlpha(kGreen, 0.5);
      }
      tg->SetMarkerSize(1.5);
      tg->SetLineWidth(tg->GetLineWidth()*2);

      if (xaxis=="cent") {
         if (thebin.centbin().low()<=0 && thebin.centbin().high()<=0) padr->cd();
         else padl->cd();

         // do not plot wide centrality bins
         prune(tg, tg_syst);
      }
      tg_syst->Draw("2");      
      gStyle->SetEndErrorSize(5);
      tg->Draw("P");
      // tg->Draw("[]");

      TString raplabel = Form("%.1f < |y| < %.1f, ",it->first.rapbin().low(),it->first.rapbin().high());
      if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f, ",it->first.rapbin().high());
      TString otherlabel = "BWAA";
      if (xaxis == "pt") otherlabel.Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
      if (xaxis == "cent") otherlabel.Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
      if (!isMB) {
         tleg->AddEntry(tg, (raplabel + otherlabel), "p");
      }

      // print tex
      ostringstream oss;
      oss.precision(1); oss.setf(ios::fixed);
      oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
      if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
      if (xaxis == "cent") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";

      addline(texname,oss.str());
      printGraph(tg, tg_syst, texname);

      // for the centrality dependence: we want Npart plotted, not the centrality
      if (xaxis == "cent") {
         centrality2npart(tg, false, isMB, (150./1.6)*it->first.rapbin().low());
         centrality2npart(tg_syst, true, isMB, (150./1.6)*it->first.rapbin().low());
      }

      // in the case where the centrality dependence is plotted: treat the PP uncertainties as global systematics
      if (xaxis == "cent" && !isMB) {
         double x, dx, y, dy;
         dx = 10;
         x = 2*dx*(thebin.rapbin().low()/1.6) + dx;
         y = 1;
         anabin thebin(it->first.rapbin().low(),
               it->first.rapbin().high(),
               it->first.ptbin().low(),
               it->first.ptbin().high(),
               0,200);
         dy = gsyst[thebin].value_dR_rel;
         cout << "global syst: " << dy << endl;
         TBox *tbox = new TBox(x-dx,y-dy,x+dx,y+dy);
         if (thebin.rapbin() == binF((float) 0.,1.6)) tbox->SetFillColorAlpha(kBlue, 0.5);
         else if (thebin.rapbin() == binF(1.6,2.4)) tbox->SetFillColorAlpha(kRed, 0.5);
         else tbox->SetFillColorAlpha(kGreen, 0.5);
         tbox->Draw();
      }

      cnt++;
      it_syst++;
   }

   // now plot the 12007 if we need to
   if (plot12007 && xaxis == "cent") {
      padl->cd();
      if (plot12007_mid) {
         TGraphAsymmErrors *g_12007_mid_cent = result12007_mid_cent();
         TGraphAsymmErrors *g_12007_mid_cent_syst = result12007_mid_cent_syst();
         g_12007_mid_cent->SetMarkerStyle(kFullTriangleUp);
         g_12007_mid_cent->SetMarkerSize(1.5);
         g_12007_mid_cent->SetMarkerColor(kCyan);
         g_12007_mid_cent->SetLineColor(kCyan);
         g_12007_mid_cent_syst->SetFillColorAlpha(kCyan, 0.5);
         g_12007_mid_cent_syst->Draw("2");
         g_12007_mid_cent->Draw("P");
         tleg->AddEntry(g_12007_mid_cent,"0 < |y| < 1.6, 6.5 < p_{T} < 30 GeV/c, #sqrt{s_{NN}} = 2.76 TeV", "lp");
      }
      if (plot12007_fwd) {
         TGraphAsymmErrors *g_12007_fwd_cent = result12007_fwd_cent();
         TGraphAsymmErrors *g_12007_fwd_cent_syst = result12007_fwd_cent_syst();
         g_12007_fwd_cent->SetMarkerStyle(kFullTriangleDown);
         g_12007_fwd_cent->SetMarkerSize(1.5);
         g_12007_fwd_cent->SetMarkerColor(kMagenta);
         g_12007_fwd_cent->SetLineColor(kMagenta);
         g_12007_fwd_cent_syst->SetFillColorAlpha(kMagenta, 0.5);
         g_12007_fwd_cent_syst->Draw("2");
         g_12007_fwd_cent->Draw("P");
         tleg->AddEntry(g_12007_fwd_cent,"1.6 < |y| < 2.4, 3 < p_{T} < 30 GeV/c, #sqrt{s_{NN}} = 2.76 TeV", "lp");
      }

      padr->cd();
      if (plot12007_mid) {
         TGraphAsymmErrors *g_12007_mid = result12007_mid();
         TGraphAsymmErrors *g_12007_mid_syst = result12007_mid_syst();
         g_12007_mid->SetMarkerStyle(kFullTriangleUp);
         g_12007_mid->SetMarkerSize(1.5);
         g_12007_mid->SetMarkerColor(kCyan);
         g_12007_mid->SetLineColor(kCyan);
         g_12007_mid_syst->SetFillColorAlpha(kCyan, 0.5);
         g_12007_mid_syst->Draw("2");
         g_12007_mid->Draw("P");
      }
      if (plot12007_fwd) {
         TGraphAsymmErrors *g_12007_fwd = result12007_fwd();
         TGraphAsymmErrors *g_12007_fwd_syst = result12007_fwd_syst();
         g_12007_fwd->SetMarkerStyle(kFullTriangleDown);
         g_12007_fwd->SetMarkerSize(1.5);
         g_12007_fwd->SetMarkerColor(kMagenta);
         g_12007_fwd->SetLineColor(kMagenta);
         g_12007_fwd_syst->SetFillColorAlpha(kMagenta, 0.5);
         g_12007_fwd_syst->Draw("2");
         g_12007_fwd->Draw("P");
      }

      padl->cd();
   }

   if (plotlimits95) plotLimits(theCats,xaxis,"../Limits/csv/Limits_95.csv",0);
   if (fiterrors && FCerrors) plotLimits(theCats,xaxis,"../Limits/csv/Limits_68.csv",xaxis=="cent" ? 5 : 1, false);

   // limits for inclusive
   if (xaxis=="cent") {
      padr->cd();
      if (plotlimits95) plotLimits(theCats,xaxis,"../Limits/csv/Limits_95.csv", 0 , true, true);
      padl->cd();
   }

   if (xaxis=="cent") padl->cd();
   tleg->Draw();

   TLatex tl;
   double tlx = 0.2+xshift;
   double tly = 0.69;
   if (!promptonly && !nonpromptonly) tl.DrawLatexNDC(tlx,tly,"Passing #font[12]{l}_{J/#psi}^{3D} cut");
   else if (promptonly) tl.DrawLatexNDC(tlx,tly,"Prompt only");
   else tl.DrawLatexNDC(tlx,tly,"Non-prompt only");

   int iPos = 33;
   CMS_lumi( (TPad*) gPad, 106, iPos, "" );
   // CMS_lumi( (TPad*) gPad, 103, iPos, "" );

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE); 
   c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s%s%s.root",outputDir.c_str(), xaxis.c_str(), nameTag, plot12007 ? "_12007" : ""));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s%s%s.png",outputDir.c_str(), xaxis.c_str(), nameTag, plot12007 ? "_12007" : ""));
   gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
   c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s%s%s.pdf",outputDir.c_str(), xaxis.c_str(), nameTag, plot12007 ? "_12007" : ""));

   delete tleg;
   delete haxes; delete haxesr;
   delete padl; delete padr;
   delete c1;

   // close tex
   closetex(texname);
   cout << "Closed " << texname << endl;
}

void centrality2npart(TGraphAsymmErrors* tg, bool issyst, bool isMB, double xshift) {
   int n = tg->GetN();
   for (int i=0; i<n; i++) {
      double x, y, exl, exh, eyl, eyh;
      x = tg->GetX()[i];
      y = tg->GetY()[i];
      exl = tg->GetErrorXlow(i);
      exh = tg->GetErrorXhigh(i);
      eyl = tg->GetErrorYlow(i);
      eyh = tg->GetErrorYhigh(i);
      x = isMB ? 150 + xshift : HI::findNpartAverage(2.*(x-exl),2.*(x+exh));
      if (!issyst) {
         if (isMB) {exl=0; exh=0;} else {
            exl = HI::findNpartSyst(2.*(x-exl),2.*(x+exh));//0.;
            exh = HI::findNpartSyst(2.*(x-exl),2.*(x+exh));//0.;
         }
      } else {
         exl = !isMB ? 5 : 5./(1.-xfrac);
         exh = exl;
      }
      tg->SetPoint(i,x,y);
      tg->SetPointError(i,exl,exh,eyl,eyh);
   }
}

void plotLimits(vector<anabin> theCats, string xaxis, const char* filename, double xshift, bool ULonly, bool isInclusive) {
   map<anabin,limits> maplim = readLimits(filename);

   map<anabin,limits>::const_iterator it;
   for (it=maplim.begin(); it!=maplim.end(); it++) {
      anabin thebin = it->first;
      if (!binok(theCats,xaxis,thebin,false)) continue;
      limits lim = it->second;
      if (ULonly && lim.val.first>0) continue; // only draw upper limits, ie interval which lower limit is 0
      bool isInclusiveBin = (xaxis=="cent" && thebin.centbin()==binI(0,200));
      if (isInclusiveBin && !isInclusive) continue;
      if (!isInclusiveBin && isInclusive) continue;
      // draw arrow in the right place and with the right color...
      Color_t color=kBlack;
      double x=0, y=0, dx=0;
      double ylow = lim.val.first;
      y = lim.val.second;
      if (xaxis=="pt") {
         double low= thebin.ptbin().low();
         double high = thebin.ptbin().high();
         x = (low+high)/2.;
         dx = 0.5;
      } else if (xaxis=="cent") {
         double low= thebin.centbin().low();
         double high = thebin.centbin().high();
         x = isInclusive ? 150 + (150./1.6)*thebin.rapbin().low(): HI::findNpartAverage(low,high);
         dx = isInclusive ? 10./(1.-xfrac) : 10;
      }
      if (thebin.rapbin() == binF(0.,1.6)) {
         color = kBlue;
      } else if (thebin.rapbin() == binF(1.6,2.4)) {
         color = kRed;
      }
      drawArrow(x+xshift, ylow, y, dx, color);
   }
}

void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color) {
   TArrow *arrow = new TArrow(x,yhigh,x,ylow<=0. ? 0.05 : ylow,0.03,ylow<=0. ? ">" : "<>");
   arrow->SetLineColor(color);
   arrow->Draw();
   if (ylow<=0.) {
      TLine *line = new TLine(x-dx,yhigh,x+dx,yhigh);
      line->SetLineColor(color);
      line->Draw();
   }
}
