#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/monster.h"
#include "Systematics/syst.h"

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
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"

using namespace std;

////////////////
// PARAMETERS //
////////////////

//#ifndef poiname_check
//#define poiname_check
//const char* poiname       = "N_Jpsi"; // for RAA (will correct automatically for efficiency)
//#endif
//const char* ylabel        = "R_{AA}";

bool  doprompt      = false;  // prompt Jpsi
bool  dononprompt   = true;  // nonprompt Jpsi
bool  is16004       = false; // plot results in 16-004 bins. If false use 16-025 bins
bool  plotPsi2S     = false; // plot Psi2S
bool  plot14005     = false; // plot results at 2.76 TeV
bool  applyEff      = true;
bool  applyAcc      = true;
bool  doLogPt       = false;
bool  includeEffSyst = true;
bool  excludeNonFitSyst = false;
bool  plotFwdMid    = false;
bool  isPreliminary = false;
string nameTag_base = "_prompt";    // can put here e.g. "_prompt", "_nonprompt", ...

const bool useNcoll = false; // false -> use TAA / NMB, true -> use Ncoll / lumiPbPb

//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void setOptions(bool adoprompt, bool adononprompt, bool aplotFwdMid, bool ais16004, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, bool aincludeEffSyst, bool  aexcludeNonFitSyst, string anameTag_base="");
void plotGraphRAA(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high);
void plotGraphXS(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high, const char* collName);
void plotCombGraphXS(map<anabin, TGraphAsymmErrors*> theGraphs_pp, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pp, map<anabin, TGraphAsymmErrors*> theGraphs_pbpb, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pbpb, string xaxis, string outputDir, map<anabin, syst> gsyst_low_pp, map<anabin, syst> gsyst_high_pp, map<anabin, syst> gsyst_low_pbpb, map<anabin, syst> gsyst_high_pbpb, const char* collName="PPandPbPb");
void plotGraphBF(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, const char* collName);
void plotCombGraphBF(map<anabin, TGraphAsymmErrors*> theGraphs_pp, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pp, map<anabin, TGraphAsymmErrors*> theGraphs_pbpb, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pbpb, string xaxis, string outputDir, const char* collName="PPandPbPb");
void plotRAA(vector<anabin> thecats, string xaxis, string workDirName);
void plotXS(vector<anabin> thecats, string xaxis, string workDirName);
void plotBF(vector<anabin> thecats, string xaxis, string workDirName);
void centrality2npart(TGraphAsymmErrors* tg, bool issyst=false, double xshift=0.);
int color(int i);
int markerstyle(int i);
string nameTag;

class raa_input {
public:
  double npp;
  double dnpp_stat;
  double systpp;
  double naa;
  double dnaa_stat;
  double systaa;
  double effpp;
  double effaa;
  double accpp;
  double accaa;
  double effppP; // Only used for the b fraction correction
  double effaaP; // Only used for the b fraction correction
  double accppP; // Only used for the b fraction correction
  double accaaP; // Only used for the b fraction correction
  double effppNP; // Only used for the b fraction correction
  double effaaNP; // Only used for the b fraction correction
  double accppNP; // Only used for the b fraction correction
  double accaaNP; // Only used for the b fraction correction
  double systeffppP; // Only used for the b fraction correction
  double systeffaaP; // Only used for the b fraction correction
  double systeffppNP; // Only used for the b fraction correction
  double systeffaaNP; // Only used for the b fraction correction
  double taa;
  double lumipp;
  double lumiaa;
  double ncoll;
  double bfracpp;
  double dbfracpp;
  double systbfracpp;
  double bfracaa;
  double dbfracaa;
  double systbfracaa;
  double correlaa;
  double correlpp;
  syst   statpp;
};

map<anabin, raa_input > readResults(const char* resultsFile);
void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color);

/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////

void plotPt(string workDirName, string poiname, int iplot) {
  // poiname is the observable to be plotted = "RAA", "BF" and "XS" (For RAA, b-fraction and cross-section)
  
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
  
  // 2 rapidity intervals for 16-004 bins
  if (iplot==3) {
    theCats.push_back(anabin(0,1.6,6.5,30,0,200));
    theCats.push_back(anabin(1.6,2.4,3.0,30,0,200));
  }
  
  
  nameTag = nameTag_base + Form("_%i",iplot);
  
  if (poiname.find("RAA")!=std::string::npos) plotRAA(theCats,xaxis,workDirName);
  else if (poiname.find("BF")!=std::string::npos && !plotPsi2S) plotBF(theCats,xaxis,workDirName);
  else if (poiname.find("XS")!=std::string::npos && !plotPsi2S) plotXS(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
  
};

void plotCent(string workDirName, string poiname, int iplot) {
  // poiname is the observable to be plotted = "RAA", "BF" and "XS" (For RAA, b-fraction and cross-section)
  
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
  
  // 2 rapidity intervals for 16-004 bins
  if (iplot==3) {
    theCats.push_back(anabin(0,1.6,6.5,30,0,200));
    theCats.push_back(anabin(1.6,2.4,3.0,30,0,200));
  }
  
  nameTag = nameTag_base + Form("_%i",iplot);
  if (poiname.find("RAA")!=std::string::npos) plotRAA(theCats,xaxis,workDirName);
  else if (poiname.find("BF")!=std::string::npos && !plotPsi2S) plotBF(theCats,xaxis,workDirName);
  else if (poiname.find("XS")!=std::string::npos && !plotPsi2S) plotXS(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};

void plotRap(string workDirName, string poiname) {
  // poiname is the observable to be plotted = "RAA", "BF" and "XS" (For RAA, b-fraction and cross-section)
  
  string xaxis = "rap";
  vector<anabin> theCats;
  
  theCats.push_back(anabin(0,2.4,6.5,50,0,200));
  
  nameTag = nameTag_base;
  if (poiname.find("RAA")!=std::string::npos) plotRAA(theCats,xaxis,workDirName);
  else if (poiname.find("BF")!=std::string::npos && !plotPsi2S) plotBF(theCats,xaxis,workDirName);
  else if (poiname.find("XS")!=std::string::npos && !plotPsi2S) plotXS(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};

void plotAll(string workDirName, string poiname) {
  // poiname is the observable to be plotted = "RAA", "BF" and "XS" (For RAA, b-fraction and cross-section)
//  if (dononprompt) nameTag_base = "_nonprompt";
  if (!doprompt && !dononprompt) nameTag_base = "";
  
  if (is16004)
  {
    plotPt(workDirName,poiname,3);
    plotCent(workDirName,poiname,3);
  }
  else{
    plotPt(workDirName,poiname,0);
    plotPt(workDirName,poiname,1);
    plotPt(workDirName,poiname,2);
    if ((poiname.find("XS")==std::string::npos) && (poiname.find("BF")==std::string::npos))
    {
      plotCent(workDirName,poiname,0);
      plotCent(workDirName,poiname,1);
      plotCent(workDirName,poiname,2);
    }
    plotRap(workDirName,poiname);
  }
};

void doAllplots(bool is16004=false) {
  //bool adoprompt, bool adononprompt, bool aplotFwdMid, bool ais16004, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, bool aincludeEffSyst, bool  aexcludeNonFitSyst, string anameTag_base=""
  
  if (!is16004)
  {
    // pr J/psi
    // noCorr with fitSyst
    setOptions(true,false,false,false,false,false,false,false,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","XS");
    
    // accEffCorr with all syst
    setOptions(true,false,false,false,false,false,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","XS");
    
    setOptions(true,false,false,false,false,true,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA"); // Comparison with 2.76
    
    // accEffCorr only 2 rap ranges
    setOptions(true,false,true,false,false,false,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");
    
    // np J/psi
    // noCorr with fitSyst
    setOptions(false,true,false,false,false,false,false,false,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","XS");
    
    // accEffCorr with all syst -> Systematics on POINTS are not working for RAA vs. pT and |y|
    setOptions(false,true,false,false,false,false,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","XS");
   
    setOptions(false,true,false,false,false,true,true,true,false,true,false);
    printOptions();
      plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA"); // Comparison with 2.76
    
    // accEffCorr only 2 rap ranges
    setOptions(false,true,true,false,false,false,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","RAA");

    
    // noCorr with fitSyst
    setOptions(false,false,false,false,false,false,false,false,false,true,false,"");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","BF");
    
    // accEffCorr with all syst
    setOptions(false,false,false,false,false,false,true,true,false,true,false,"");
    plotAll("DataFits_16025_2D_2CB_polBkg_nominal","BF");
    
    // Tables of systematics
//    map<anabin, syst> syst_PP = readSyst_all("PP","NJpsi_prompt","16025",true,true,"",true,"Systematics/systs_PP_prompt_all.tex");
//    syst_PP = readSyst_all("PP","NJpsi_nonprompt","16025",true,true,"",true,"Systematics/systs_PP_nonprompt_all.tex");
//    
//    map<anabin, syst> syst_PbPb = readSyst_all("PbPb","NJpsi_prompt","16025",true,true,"",true,"Systematics/systs_PbPb_prompt_all.tex");
//    syst_PbPb = readSyst_all("PbPb","NJpsi_nonprompt","16025",true,true,"",true,"Systematics/systs_PbPb_nonprompt_all.tex");
  }
  else
  {
//    setOptions(true,false,true,false,false,false,false,false,false,true);
//    printOptions();
//    plotAll("DataFits_16004_2D_2CB_polBkg_nominal","RAA");
//    
//    
//    setOptions(true,false,true,false,false,true,false,false,false,true);
//    printOptions();
//    plotAll("DataFits_16004_2D_2CB_polBkg_nominal","RAA");
    
//    setOptions(true,false,true,false,false,true,true,false,true,false);
//    printOptions();
//    plotAll("DataFits_16004_2D_2CB_polBkg_nominal","RAA");
    
  //bool adoprompt, bool adononprompt, bool ais16004, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, bool aincludeEffSyst, bool  aexcludeNonFitSyst, string anameTag_base=""
    setOptions(true,false,false,true,true,false,true,true,false,true,false);
    printOptions();
    plotAll("DataFits_16004_2D_2CB_polBkg_nominal","RAA");
  }

};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plotRAA(vector<anabin> thecats, string xaxis, string outputDir) {
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent"
  // outputDir is the directory to save the plots
  
  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  
  TString poi("NJpsi");
  if (doprompt) poi = "NJpsi_prompt";
  if (dononprompt) poi = "NJpsi_nonprompt";
  
  TFile *f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
  if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),"DATA");
    f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  
  TString sTag("16025");
  if (is16004) sTag = "16004";
  
  map<anabin, raa_input> theVars_inputs;
  
  map<anabin, syst> syst_PP = readSyst_all("PP",poi.Data(),sTag.Data(),includeEffSyst,false);
  map<anabin, syst> syst_PbPb = readSyst_all("PbPb",poi.Data(),sTag.Data(),includeEffSyst,false);
  map<anabin, syst> syst_taa_low = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_low.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_taa_high = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_high.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_Nmb = readSyst(Form("Systematics/csv/syst_%s_PbPb_Nmb.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_lumipp = readSyst(Form("Systematics/csv/syst_%s_PP_lumi.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> stat_PP; // for PP statistics
  map<anabin, syst> syst_glb_low; // for the boxes at 1
  map<anabin, syst> syst_glb_high; // for the boxes at 1
  
  vector<double> x, ex, y, ey;
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  float eff, acc, lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL,bfrac_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
//  tr->SetBranchAddress("N_Jpsi_val",&val);
//  tr->SetBranchAddress("N_Jpsi_errL",&errL);
//  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  if (!dononprompt)
  {
    tr->SetBranchAddress("eff_val",&eff);
    tr->SetBranchAddress("acc_val",&acc);
  }
  else
  {
    tr->SetBranchAddress("effnp_val",&eff);
    tr->SetBranchAddress("accnp_val",&acc);
  }
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    
    if (xaxis=="rap" && ((ymin==0 && ymax<=0.61 && ymax>=0.59 ) || (ymin>=0.59 && ymin<=0.61 && ymax>=1.19 && ymax <=1.21) || (ymin>=1.19 && ymin<=1.21 && ymax>=1.79 && ymax <=1.81) || (ymin>=1.79 && ymin<=1.81 && ymax>=2.39 && ymax <=2.41))) continue;
    
    anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);
    
    bool ispp = (TString(collSystem)=="PP");
    
    if (ispp) {
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp_stat = errL;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_errL;
      theVars_inputs[thebin].systpp = syst_PP[thebin].value;
      syst thestat_PP; thestat_PP.name = "stat_PP"; thestat_PP.value = errL/val;
      stat_PP[thebin] = thestat_PP;
      theVars_inputs[thebin].statpp = thestat_PP;
      theVars_inputs[thebin].lumipp = lumi;
      theVars_inputs[thebin].effpp = eff;
      theVars_inputs[thebin].accpp = acc;
      theVars_inputs[thebin].correlpp = correl;
    } else {
      theVars_inputs[thebin].naa = val;
      theVars_inputs[thebin].dnaa_stat = errL;
      theVars_inputs[thebin].bfracaa = bfrac;
      theVars_inputs[thebin].dbfracaa = bfrac_errL;
      theVars_inputs[thebin].systaa = syst_PbPb[thebin].value;
      theVars_inputs[thebin].lumiaa = lumi;
      theVars_inputs[thebin].effaa = eff;
      theVars_inputs[thebin].accaa = acc;
      theVars_inputs[thebin].taa = taa;
      theVars_inputs[thebin].ncoll = ncoll;
      theVars_inputs[thebin].correlaa = correl;
    }
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned;
  map<anabin, vector<double> > theVarsBinned_stat;
  map<anabin, vector<double> > theVarsBinned_syst_low;
  map<anabin, vector<double> > theVarsBinned_syst_high;
  map<anabin, TGraphAsymmErrors* > theGraphs;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst;
  
  map<anabin, raa_input > theResults16004; //In case they are needed to compute the Psi2S Raa
  if (is16004 && plotPsi2S) theResults16004 = readResults("raa_16004.csv");
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned[*it] = vector<double>();
  }
  
  for (map<anabin, raa_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    raa_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    raa_input spp = theVars_inputs[thebinPP];
//    cout << "Result: " << endl;
//    thebinOrig.print();
//    cout << s.naa  << " +- " << s.dnaa_stat << " ; " << spp.npp << " +- " << spp.dnpp_stat << " ; " << s.bfracaa << " ; " << spp.bfracpp << endl;
    if (s.naa <= 0 || spp.npp <= 0) continue;
    if (applyEff && (s.effaa <= 0 || spp.effpp <= 0))
    {
      thebinOrig.print();
      cout << s.effaa << " " << spp.effpp << endl;
      cout << "[WARNING] : You want to correct by efficiency but the value is <=0 for this bin " << endl;
      continue;
    }
    if (applyAcc && (s.accaa <= 0 || spp.accpp <= 0))
    {
      thebinOrig.print();
      cout << "[WARNING] : You want to correct by acceptance but the value is <=0 for this bin " << endl;
      continue;
    }
    if ((doprompt || dononprompt) && (spp.bfracpp<=0 || s.bfracaa<=0)) continue;
    
    theBins[thebin].push_back(it->first);
    
    double normfactorpp = 1., normfactoraa = 1.;
    normfactorpp = 1./spp.lumipp;
    
    if (useNcoll) {
      normfactoraa = 1./s.lumiaa;
      normfactoraa *= 1./(208.*208.*(HI::findNcollAverage(it->first.centbin().low(),it->first.centbin().high())/HI::findNcollAverage(0,200)));
    } else {
      double myNmb = NMB * s.lumiaa / lumipbpb_ABCD;
      normfactoraa = 1./(myNmb*s.taa*1e-3); // the 1e-3 factor is because taa is in mb-1 while lumis are in mub-1
//      thebinOrig.print();
//      cout << "Taa = " << s.taa << endl;
    }
    normfactoraa *= 200./(it->first.centbin().high()-it->first.centbin().low());
    
    if (applyEff && applyAcc)
    {
      normfactorpp = normfactorpp / (spp.effpp*spp.accpp);
      normfactoraa = normfactoraa / (s.effaa*s.accaa);
    }
    else if (applyEff)
    {
      normfactorpp = normfactorpp / (spp.effpp);
      normfactoraa = normfactoraa / (s.effaa);
    }

    double naa = s.naa;
    double npp = spp.npp;
    double dnaa = s.dnaa_stat;
    double dnpp = spp.dnpp_stat;
//    cout << "correl b-NJpsi = " << spp.bfracpp << " ; " << s.correlaa << endl;
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
//    cout << "after prompt = " << naa << " ; " << npp << " ; " << dnaa << " ; " << dnpp << endl;
    naa *= normfactoraa;
    npp *= normfactorpp;
    dnaa *= normfactoraa;
    dnpp *= normfactorpp;
    
    double raa = npp>0 ? naa / npp : 0;
    double draa = raa>0 ? raa*sqrt(pow(dnaa/naa,2) + pow(dnpp/npp,2)) : 0;
//    cout << it->first.ptbin().low() << " " << it->first.ptbin().high() << " -> " << raa << " " << s.naa << " " << spp.npp << " " << s.effaa << " "  << spp.effpp << " " << s.accaa << " " << s.accpp << ", " << s.bfracaa << " " << spp.bfracpp <<  " " << s.taa << endl;
    double syst_low = raa*sqrt(pow(spp.systpp,2)+pow(s.systaa,2));
    double syst_high = syst_low;
    
    // case of the centrality dependence: factor out pp uncertainties, but include taa
    if (xaxis=="cent") {
      draa = raa*dnaa/naa;
      syst_low = raa*sqrt(pow(s.systaa,2)+pow(syst_taa_low[it->first].value,2));
      syst_high = raa*sqrt(pow(s.systaa,2)+pow(syst_taa_high[it->first].value,2));
    }
    
    if (is16004 && plotPsi2S) // Compute Psi2S Raa from the product of Jpsi Raa and 16-004 double ratio
    {
//      thebin.print();
      raa_input doubleR16004 = theResults16004[thebinOrig];
      
      if (doubleR16004.naa == 0) // It means that the bin has confidence intervals
      {
//        thebinOrig.print();
//        cout << "raa = " << raa << " ; dnaa_stat = " << doubleR16004.dnaa_stat << " ; systaa = " << doubleR16004.systaa << endl;
        draa = raa*doubleR16004.dnaa_stat; // Lower limit
        syst_low =0;
        syst_high = raa*doubleR16004.systaa; // Upper limit
        raa = 0.;
        
//        cout << raa << " ; " << draa << " ; " << syst << endl;
      }
      else
      {
        draa = TMath::Sqrt(TMath::Power(doubleR16004.dnaa_stat/doubleR16004.naa,2.) + TMath::Power(draa/raa,2.));
        syst_low = TMath::Sqrt(TMath::Power(doubleR16004.systaa/doubleR16004.naa,2.) + TMath::Power(syst_low/raa,2.));
        syst_high = TMath::Sqrt(TMath::Power(doubleR16004.systaa/doubleR16004.naa,2.) + TMath::Power(syst_high/raa,2.));
        
        raa = raa*doubleR16004.naa;
        draa = raa*draa;
        syst_low = raa*syst_low;
        syst_high = raa*syst_high;
      }
    }
//    thebinOrig.print();
//    cout << "Raa = " << raa << " +- " << draa << " +- " << syst << endl;
    theVarsBinned[thebin].push_back(raa);
    theVarsBinned_stat[thebin].push_back(draa);
    theVarsBinned_syst_low[thebin].push_back(syst_low);
    theVarsBinned_syst_high[thebin].push_back(syst_high);
  }
  
  // systematics
  if (xaxis=="cent") { // put the PP stat error into the PP syst, that will go into a box at 1
    vector< map<anabin, syst> > all_glb;
    all_glb.push_back(syst_PP);
    all_glb.push_back(stat_PP);
    all_glb.push_back(syst_lumipp);
    all_glb.push_back(syst_Nmb);
    syst_glb_low = combineSyst(all_glb,"global");
    syst_glb_high = syst_glb_low;
  } else {
    vector< map<anabin, syst> > all_glb_low;
    all_glb_low.push_back(syst_taa_low);
    all_glb_low.push_back(syst_lumipp);
    all_glb_low.push_back(syst_Nmb);
    syst_glb_low = combineSyst(all_glb_low,"global");
    
    vector< map<anabin, syst> > all_glb_high;
    all_glb_high.push_back(syst_taa_high);
    all_glb_high.push_back(syst_lumipp);
    all_glb_high.push_back(syst_Nmb);
    syst_glb_high = combineSyst(all_glb_high,"global");
  }
  
  // make TGraphAsymmErrors
  int cnt=0;
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    int n = theBins[*it].size();
    if(n==0) {
      cout << "Error, nothing found for category" << endl;
      theGraphs[*it] = NULL;
      continue;
    }
    
    theGraphs[*it] = new TGraphAsymmErrors(n);
    theGraphs[*it]->SetName(Form("bin_%i",cnt));
    theGraphs_syst[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst[*it]->SetName(Form("bin_%i_syst",cnt));
    
    for (int i=0; i<n; i++) {
      double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
      double exsyst=0, eysyst_low=0,eysyst_high=0;
      double low=0, high=0;
      anabin thebin = theBins[*it][i];
      y = theVarsBinned[*it][i];
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
        exsyst = (xaxis=="pt") ? 0.5 : 0.05;
        eysyst_low = theVarsBinned_syst_low[*it][i];
        eysyst_high = theVarsBinned_syst_high[*it][i];
      }
      if (xaxis=="cent") {
        low= thebin.centbin().low();
        high = thebin.centbin().high();
        // exl = 0.;
        // exh = 0.;
        x = (low+high)/2./2.;
        exh = (high-low)/2./2.;
        exl = (high-low)/2./2.;
        exsyst = exl;
        eysyst_low = theVarsBinned_syst_low[*it][i];
        eysyst_high = theVarsBinned_syst_high[*it][i];
      }
      eyl = fabs(theVarsBinned_stat[*it][i]);
      eyh = eyl;
      
      // eysyst = y*eysyst;
      
      theGraphs[*it]->SetPoint(i,x,y);
      theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
      theGraphs_syst[*it]->SetPoint(i,x,y);
      theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst_low,eysyst_high);
//       cout << "final = " << x << " " << y << " " << eyl << " " << eyh << " " << eysyst << endl;
      
      // theGraphs[*it]->Sort();
      // theGraphs_syst[*it]->Sort();
    }
    cnt++;
  }
  
  // plot
  plotGraphRAA(theGraphs, theGraphs_syst, xaxis, outputDir, syst_glb_low, syst_glb_high);
}

void plotXS(vector<anabin> thecats, string xaxis, string outputDir) {
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent"
  // outputDir is the directory to save the plots
  
  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  
  TString poi("NJpsi");
  if (doprompt) poi = "NJpsi_prompt";
  if (dononprompt) poi = "NJpsi_nonprompt";
  
  if (xaxis=="cent") {
    cout << "ERROR cross-sections vs centrality not implemented" << endl; // Important note: if you remove this, remeber to include Taa uncertainty
    return;
  }
  
  TFile *f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
  if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),"DATA");
    f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  
  TString sTag("16025");
  if (is16004) sTag = "16004";
  
  map<anabin, raa_input> theVars_inputs;
  
  map<anabin, syst> syst_PP = readSyst_all("PP",poi.Data(),sTag.Data(),includeEffSyst,includeEffSyst);
  map<anabin, syst> syst_PbPb = readSyst_all("PbPb",poi.Data(),sTag.Data(),includeEffSyst,includeEffSyst);
  map<anabin, syst> syst_taa_low = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_low.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_taa_high = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_high.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_Nmb = readSyst(Form("Systematics/csv/syst_%s_PbPb_Nmb.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_lumipp = readSyst(Form("Systematics/csv/syst_%s_PP_lumi.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_glb_aa_low; // for the boxes at 1
  map<anabin, syst> syst_glb_aa_high; // for the boxes at 1
  map<anabin, syst> syst_glb_pp; // for the boxes at 1
  
  vector<double> x, ex, y, ey;
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  float eff, acc, lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL,bfrac_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
//  tr->SetBranchAddress("N_Jpsi_val",&val);
//  tr->SetBranchAddress("N_Jpsi_errL",&errL);
//  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  if (!dononprompt)
  {
    tr->SetBranchAddress("eff_val",&eff);
    tr->SetBranchAddress("acc_val",&acc);
  }
  else
  {
    tr->SetBranchAddress("effnp_val",&eff);
    tr->SetBranchAddress("accnp_val",&acc);
  }
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {

    tr->GetEntry(i);
    
    if (xaxis=="rap" && ((ymin==0 && ymax<=0.61 && ymax>=0.59 ) || (ymin>=0.59 && ymin<=0.61 && ymax>=1.19 && ymax <=1.21) || (ymin>=1.19 && ymin<=1.21 && ymax>=1.79 && ymax <=1.81) || (ymin>=1.79 && ymin<=1.81 && ymax>=2.39 && ymax <=2.41))) continue;
    
    anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);
    
    bool ispp = (TString(collSystem)=="PP");
    
    if (ispp) {
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp_stat = errL;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_errL;
      theVars_inputs[thebin].systpp = syst_PP[thebin].value;
      theVars_inputs[thebin].lumipp = lumi;
      theVars_inputs[thebin].effpp = eff;
      theVars_inputs[thebin].accpp = acc;
      theVars_inputs[thebin].correlpp = correl;
    } else {
      theVars_inputs[thebin].naa = val;
      theVars_inputs[thebin].dnaa_stat = errL;
      theVars_inputs[thebin].bfracaa = bfrac;
      theVars_inputs[thebin].dbfracaa = bfrac_errL;
      theVars_inputs[thebin].systaa = syst_PbPb[thebin].value;
      theVars_inputs[thebin].lumiaa = lumi;
      theVars_inputs[thebin].effaa = eff;
      theVars_inputs[thebin].accaa = acc;
      theVars_inputs[thebin].taa = taa;
      theVars_inputs[thebin].ncoll = ncoll;
      theVars_inputs[thebin].correlaa = correl;
    }
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned_aa;
  map<anabin, vector<double> > theVarsBinned_pp;
  map<anabin, vector<double> > theVarsBinned_stat_aa;
  map<anabin, vector<double> > theVarsBinned_stat_pp;
  map<anabin, vector<double> > theVarsBinned_syst_aa;
  map<anabin, vector<double> > theVarsBinned_syst_pp;
  map<anabin, TGraphAsymmErrors* > theGraphs_aa;
  map<anabin, TGraphAsymmErrors* > theGraphs_pp;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst_aa;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst_pp;
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned_aa[*it] = vector<double>();
    theVarsBinned_pp[*it] = vector<double>();
  }
  
  for (map<anabin, raa_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    raa_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    raa_input spp = theVars_inputs[thebinPP];
    
    if (s.naa <= 0 || spp.npp <= 0) continue;
    if (applyEff && (s.effaa <= 0 || spp.effpp <= 0))
    {
      thebinOrig.print();
      cout << "[WARNING] : You want to correct by efficiency but the value is <=0 for this bin " << endl;
      continue;
    }
    if (applyAcc && (s.accaa <= 0 || spp.accpp <= 0))
    {
      thebinOrig.print();
      cout << "[WARNING] : You want to correct by acceptance but the value is <=0 for this bin " << endl;
      continue;
    }
    if ((doprompt || dononprompt) && (spp.bfracpp<=0 || s.bfracaa<=0)) continue;
    
    theBins[thebin].push_back(it->first);
    
    double deltaPt = it->first.ptbin().high() - it->first.ptbin().low();
    double deltaRap = it->first.rapbin().high() - it->first.rapbin().low();
    
//    cout << "deltaPt = " << deltaPt << " ; deltaY = " << deltaRap << endl;
    
    double normfactorpp = 1., normfactoraa = 1.;
    normfactorpp = 1./(spp.lumipp*deltaPt*deltaRap*1e-3); // the 1e-3 factor is because lumis are in mub-1 and we give xsection in nb-1
    normfactoraa = 1./(deltaPt*deltaRap*1e-3);
    
    if (useNcoll) {
      normfactoraa = 1./(s.lumiaa*deltaPt*deltaRap*1e-3);
      normfactoraa *= 1./(208.*208.*(HI::findNcollAverage(it->first.centbin().low(),it->first.centbin().high())/HI::findNcollAverage(0,200)));
    } else {
      double myNmb = NMB * s.lumiaa / lumipbpb_ABCD;
      normfactoraa *= 1./(myNmb*s.taa*1e-3); // the 1e-3 factor is because taa is in mb-1 while lumis are in mub-1
    }
    
    normfactoraa *= 200./(it->first.centbin().high()-it->first.centbin().low());
    
    if (applyEff && applyAcc)
    {
      normfactorpp = normfactorpp / (spp.effpp*spp.accpp);
      normfactoraa = normfactoraa / (s.effaa*s.accaa);
    }
    else if (applyEff)
    {
      normfactorpp = normfactorpp / (spp.effpp);
      normfactoraa = normfactoraa / (s.effaa);
    }
    
    double naa = s.naa;
    double npp = spp.npp;
    double dnaa = s.dnaa_stat;
    double dnpp = spp.dnpp_stat;
    
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
    
    
    double xsec_aa = naa*normfactoraa;
    double delta_xsec_aa = dnaa*normfactoraa;
    
    double xsec_pp = npp*normfactorpp;
    double delta_xsec_pp = dnpp*normfactorpp;
    
    double syst_xsec_aa = xsec_aa*sqrt(pow(s.systaa,2));
    double syst_xsec_pp = xsec_pp*sqrt(pow(spp.systpp,2));
    
//    cout << "PbPb: " << xsec_aa << " pp : " << xsec_pp << endl;
    
    theVarsBinned_aa[thebin].push_back(xsec_aa);
    theVarsBinned_pp[thebin].push_back(xsec_pp);
    theVarsBinned_stat_aa[thebin].push_back(delta_xsec_aa);
    theVarsBinned_stat_pp[thebin].push_back(delta_xsec_pp);
    theVarsBinned_syst_aa[thebin].push_back(syst_xsec_aa);
    theVarsBinned_syst_pp[thebin].push_back(syst_xsec_pp);
  }
  
  // systematics
  vector< map<anabin, syst> > all_glb_aa_low;
  all_glb_aa_low.push_back(syst_taa_low);
  all_glb_aa_low.push_back(syst_Nmb);
  syst_glb_aa_low = combineSyst(all_glb_aa_low,"global");
  
  vector< map<anabin, syst> > all_glb_aa_high;
  all_glb_aa_high.push_back(syst_taa_high);
  all_glb_aa_high.push_back(syst_Nmb);
  syst_glb_aa_high = combineSyst(all_glb_aa_high,"global");
  
  vector< map<anabin, syst> > all_glb_pp;
  all_glb_pp.push_back(syst_lumipp);
  syst_glb_pp = combineSyst(all_glb_pp,"global");
  
  // make TGraphAsymmErrors
  int cnt=0;
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    int n = theBins[*it].size();
    if(n==0) {
      cout << "Error, nothing found for category" << endl;
      theGraphs_aa[*it] = NULL;
      theGraphs_pp[*it] = NULL;
      continue;
    }
    
    theGraphs_aa[*it] = new TGraphAsymmErrors(n);
    theGraphs_aa[*it]->SetName(Form("bin_%i_aa",cnt));
    theGraphs_syst_aa[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst_aa[*it]->SetName(Form("bin_%i_syst_aa",cnt));
    
    theGraphs_pp[*it] = new TGraphAsymmErrors(n);
    theGraphs_pp[*it]->SetName(Form("bin_%i_pp",cnt));
    theGraphs_syst_pp[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst_pp[*it]->SetName(Form("bin_%i_syst_pp",cnt));
    
    for (int i=0; i<n; i++) {
      double x=0, exl=0, exh=0, yaa=0, eylaa=0, eyhaa=0, ypp=0, eylpp=0, eyhpp=0;
      double exsyst=0;
      double eysyst_aa=0, eysyst_pp=0;
      double low=0, high=0;
      anabin thebin = theBins[*it][i];
      yaa = theVarsBinned_aa[*it][i];
      ypp = theVarsBinned_pp[*it][i];
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
        exsyst = (xaxis=="pt") ? 0.5 : 0.05;
        eysyst_aa = theVarsBinned_syst_aa[*it][i];
        eysyst_pp = theVarsBinned_syst_pp[*it][i];
      }
      eylaa = fabs(theVarsBinned_stat_aa[*it][i]);
      eyhaa = eylaa;
      eylpp = fabs(theVarsBinned_stat_pp[*it][i]);
      eyhpp = eylpp;
      
      theGraphs_aa[*it]->SetPoint(i,x,yaa);
      theGraphs_aa[*it]->SetPointError(i,exl,exh,eylaa,eyhaa);
      theGraphs_syst_aa[*it]->SetPoint(i,x,yaa);
      theGraphs_syst_aa[*it]->SetPointError(i,exsyst,exsyst,eysyst_aa,eysyst_aa);
      // cout << x << " " << yaa << " " << eylaa << " " << eyhaa << " " << eysyst_aa << endl;
      
      theGraphs_pp[*it]->SetPoint(i,x,ypp);
      theGraphs_pp[*it]->SetPointError(i,exl,exh,eylpp,eyhpp);
      theGraphs_syst_pp[*it]->SetPoint(i,x,ypp);
      theGraphs_syst_pp[*it]->SetPointError(i,exsyst,exsyst,eysyst_pp,eysyst_pp);
      // cout << x << " " << ypp << " " << eylpp << " " << eyhpp << " " << eysyst_pp << endl;
      
      // theGraphs[*it]->Sort();
      // theGraphs_syst[*it]->Sort();
    }
    cnt++;
  }
  
  // plot
  plotGraphXS(theGraphs_aa, theGraphs_syst_aa, xaxis, outputDir, syst_glb_aa_low, syst_glb_aa_high,"pbpb");
  plotGraphXS(theGraphs_pp, theGraphs_syst_pp, xaxis, outputDir, syst_glb_pp, syst_glb_pp, "pp");
  plotCombGraphXS(theGraphs_pp, theGraphs_syst_pp, theGraphs_aa, theGraphs_syst_aa, xaxis, outputDir, syst_glb_pp, syst_glb_pp, syst_glb_aa_low, syst_glb_aa_high);
}

void plotBF(vector<anabin> thecats, string xaxis, string outputDir) {
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent"
  // outputDir is the directory to save the plots
  
  if (doprompt || dononprompt) {
    cout << "ERROR for b fractions both doprompt and dononprompt have to be set to false." << endl;
    return;
  }
  
  TFile *f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
  if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),"DATA");
    f = new TFile(treeFileName(outputDir.c_str(),"DATA"));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  
  TString sTag("16025");
  if (is16004) sTag = "16004";
  
  map<anabin, raa_input> theVars_inputs;
  
  map<anabin, syst> syst_PP = readSyst_all("PP","BJpsi",sTag.Data(),false,false);
  map<anabin, syst> syst_PbPb = readSyst_all("PbPb","BJpsi",sTag.Data(),false,false);
  map<anabin, syst> systEffP_PP = readSyst_all("PP","NJpsi_prompt",sTag.Data(),true,true);
  map<anabin, syst> systEffP_PbPb = readSyst_all("PbPb","NJpsi_prompt",sTag.Data(),true,true);
  map<anabin, syst> systEffNP_PP = readSyst_all("PP","NJpsi_nonprompt",sTag.Data(),true,true);
  map<anabin, syst> systEffNP_PbPb = readSyst_all("PbPb","NJpsi_nonprompt",sTag.Data(),true,true);
  
  vector<double> x, ex, y, ey;
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  float effp, accp,effnp, accnp;
  float val, err=0;
  float bfrac, bfrac_err;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_err);
  tr->SetBranchAddress("eff_val",&effp);
  tr->SetBranchAddress("acc_val",&accp);
  tr->SetBranchAddress("effnp_val",&effnp);
  tr->SetBranchAddress("accnp_val",&accnp);
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    
    tr->GetEntry(i);
    
    bool ispp = (TString(collSystem)=="PP");
    
    if (xaxis=="rap" && ((ymin==0 && ymax<=0.61 && ymax>=0.59 ) || (ymin>=0.59 && ymin<=0.61 && ymax>=1.19 && ymax <=1.21) || (ymin>=1.19 && ymin<=1.21 && ymax>=1.79 && ymax <=1.81) || (ymin>=1.79 && ymin<=1.81 && ymax>=2.39 && ymax <=2.41))) continue;
    if (xaxis=="pt" && (ptmin>=2.99 && ptmin<=3.01 && ptmax<=6.51 && ptmax>=6.49)) continue;
    //fixing the error bin problems?
    //if (xaxis=="pt" && ((centmin == 0) && (centmax == 200) && ymin==0 && ymax>=2.39 && ymax <=2.41 && ((ptmin>=14.99 && ptmin<15.01 && ptmax <= 20.01 && ptmax > 19.99) || (ptmin>=19.99 && ptmin<20.01 && ptmax <= 30.01 && ptmax > 29.99) || (ptmin>=29.99 && ptmin<30.01 && ptmax <= 50.01 && ptmax > 49.99)))) continue;
    
    anabin thebin(ymin, ymax, ptmin, ptmax, centmin, centmax);

    if (ispp) {
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_err;
      theVars_inputs[thebin].effppP = effp;
      theVars_inputs[thebin].systeffppP = systEffP_PP[thebin].value;
      theVars_inputs[thebin].accppP = accp;
      theVars_inputs[thebin].effppNP = effnp;
      theVars_inputs[thebin].systeffppNP = systEffNP_PP[thebin].value;
      theVars_inputs[thebin].accppNP = accnp;
      theVars_inputs[thebin].systbfracpp = syst_PP[thebin].value;
    } else {
      theVars_inputs[thebin].bfracaa = bfrac;
      theVars_inputs[thebin].dbfracaa = bfrac_err;
      theVars_inputs[thebin].effaaP = effp;
      theVars_inputs[thebin].systeffaaP = systEffP_PbPb[thebin].value;
      theVars_inputs[thebin].accaaP = accp;
      theVars_inputs[thebin].effaaNP = effnp;
      theVars_inputs[thebin].systeffaaNP = systEffNP_PbPb[thebin].value;
      theVars_inputs[thebin].accaaNP = accnp;
      theVars_inputs[thebin].systbfracaa = syst_PbPb[thebin].value;
    }
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned_aa;
  map<anabin, vector<double> > theVarsBinned_pp;
  map<anabin, vector<double> > theVarsBinned_stat_aa;
  map<anabin, vector<double> > theVarsBinned_stat_pp;
  map<anabin, vector<double> > theVarsBinned_syst_aa;
  map<anabin, vector<double> > theVarsBinned_syst_pp;
  map<anabin, TGraphAsymmErrors* > theGraphs_aa;
  map<anabin, TGraphAsymmErrors* > theGraphs_pp;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst_aa;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst_pp;
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned_aa[*it] = vector<double>();
    theVarsBinned_pp[*it] = vector<double>();
  }
  
  for (map<anabin, raa_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    raa_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    raa_input spp = theVars_inputs[thebinPP];
    
    if ((doprompt || dononprompt) && (spp.effppP<=0 || s.bfracaa<=0)) continue;
    if (applyEff && (s.effaaP <= 0 || spp.effppP <= 0 || s.effaaNP <= 0 || spp.effppNP <= 0))
    {
      thebinOrig.print();
      cout << "[WARNING] : You want to correct by efficiency but the value is <=0 for this bin " << endl;
      continue;
    }
    if (applyAcc && (s.accaaP <= 0 || spp.accppP <= 0 || s.accaaNP <= 0 || spp.accppNP <= 0))
    {
      thebinOrig.print();
      cout << "[WARNING] : You want to correct by acceptance but the value is <=0 for this bin " << endl;
      continue;
    }
    
    theBins[thebin].push_back(it->first);
    
    double bf_aa = s.bfracaa;
    double delta_bf_aa = s.dbfracaa;
    
    double bf_pp = spp.bfracpp;
    double delta_bf_pp = spp.dbfracpp;
    
    double accxEffaa_p(1.);
    double accxEffaa_np(1.);
    
    double accxEffpp_p(1.);
    double accxEffpp_np(1.);
    
    double systEff_pp(0.);
    double systEff_aa(0.);
    double corrFactor_pp(1.);
    double corrFactor_aa(1.);
    if (applyEff && applyAcc)
    {
      accxEffaa_p = (s.effaaP*s.accaaP);
      accxEffaa_np = (s.effaaNP*s.accaaNP);
      
      accxEffpp_p = (spp.effppP*spp.accppP);
      accxEffpp_np = (spp.effppNP*spp.accppNP);
    
      
//      systEff_aa = abs(bf_aa -(bf_aa*(accxEffaa_p/((bf_aa*accxEffaa_p) + (1.-bf_aa)*accxEffaa_np))))/bf_aa;
//      systEff_pp = abs(bf_pp -(bf_pp*(accxEffpp_p/((bf_pp*accxEffpp_p) + (1.-bf_pp)*accxEffpp_np))))/bf_pp;
    }
    else if (applyEff)
    {
      accxEffaa_p = s.effaaP;
      accxEffaa_np = s.effaaNP;
      
      accxEffpp_p = spp.effppP;
      accxEffpp_np = spp.effppNP;
      
//      systEff_aa = abs(bf_aa -(bf_aa*(effaa_p/((bf_aa*effaa_p) + (1.-bf_aa)*effaa_np))))/bf_aa;
//      systEff_pp = abs(bf_pp -(bf_pp*(effpp_p/((bf_pp*effpp_p) + (1.-bf_pp)*effpp_np))))/bf_pp;
    }
    
    corrFactor_pp = (accxEffpp_p/((bf_pp*accxEffpp_p) + (1.-bf_pp)*accxEffpp_np));
    bf_pp *= corrFactor_pp;
    delta_bf_pp *= corrFactor_pp;
    
//    systEff_pp =(accxEffpp_np/accxEffpp_p)*corrFactor_pp*sqrt(pow(spp.systbfracpp/spp.bfracpp,2) + pow(1-spp.bfracpp,2)*(pow(spp.systeffaaP/accxEffpp_p,2)+pow(spp.systeffppNP/accxEffpp_np,2)));
    systEff_pp =(accxEffpp_np/accxEffpp_p)*corrFactor_pp*sqrt(pow(spp.systbfracpp,2) + pow(1-spp.bfracpp,2)*(pow(spp.systeffaaP,2)+pow(spp.systeffppNP,2)));
    
    corrFactor_aa = (accxEffaa_p/((bf_aa*accxEffaa_p) + (1.-bf_aa)*accxEffaa_np));
    bf_aa *= corrFactor_aa;
    delta_bf_aa *= corrFactor_aa;
    
//    systEff_aa = (accxEffaa_np/accxEffaa_p)*corrFactor_aa*sqrt(pow(s.systbfracaa/s.bfracaa,2) + pow(1-s.bfracaa,2)*(pow(s.systeffaaP/accxEffaa_p,2)+pow(s.systeffaaNP/accxEffaa_np,2)));
    systEff_aa = (accxEffaa_np/accxEffaa_p)*corrFactor_aa*sqrt(pow(s.systbfracaa,2) + pow(1-s.bfracaa,2)*(pow(s.systeffaaP,2)+pow(s.systeffaaNP,2)));
    
//    double syst_bf_aa = bf_aa*sqrt(pow(s.systbfracaa,2) + pow(systEff_aa,2));
//    double syst_bf_pp = bf_pp*sqrt(pow(spp.systbfracpp,2) + pow(systEff_pp,2));

    double syst_bf_aa = bf_aa*sqrt(pow(s.systbfracaa,2) + pow(systEff_aa,2));
    double syst_bf_pp = bf_pp*sqrt(pow(spp.systbfracpp,2) + pow(systEff_pp,2));
    
    theVarsBinned_aa[thebin].push_back(bf_aa);
    theVarsBinned_pp[thebin].push_back(bf_pp);
    theVarsBinned_stat_aa[thebin].push_back(delta_bf_aa);
    theVarsBinned_stat_pp[thebin].push_back(delta_bf_pp);
    theVarsBinned_syst_aa[thebin].push_back(syst_bf_aa);
    theVarsBinned_syst_pp[thebin].push_back(syst_bf_pp);
  }
  
  // make TGraphAsymmErrors
  int cnt=0;
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    int n = theBins[*it].size();
    if(n==0) {
      cout << "Error, nothing found for category" << endl;
      theGraphs_aa[*it] = NULL;
      theGraphs_pp[*it] = NULL;
      continue;
    }
    
    theGraphs_aa[*it] = new TGraphAsymmErrors(n);
    theGraphs_aa[*it]->SetName(Form("bin_%i_aa",cnt));
    theGraphs_syst_aa[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst_aa[*it]->SetName(Form("bin_%i_syst_aa",cnt));
    
    theGraphs_pp[*it] = new TGraphAsymmErrors(n);
    theGraphs_pp[*it]->SetName(Form("bin_%i_pp",cnt));
    theGraphs_syst_pp[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst_pp[*it]->SetName(Form("bin_%i_syst_pp",cnt));
    
    for (int i=0; i<n; i++) {
      double x=0, exl=0, exh=0, yaa=0, eylaa=0, eyhaa=0, ypp=0, eylpp=0, eyhpp=0;
      double exsyst=0;
      double eysyst_aa=0, eysyst_pp=0;
      double low=0, high=0;
      anabin thebin = theBins[*it][i];
      yaa = theVarsBinned_aa[*it][i];
      ypp = theVarsBinned_pp[*it][i];
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
        exsyst = (xaxis=="pt") ? 0.5 : 0.05;
        eysyst_aa = theVarsBinned_syst_aa[*it][i];
        eysyst_pp = theVarsBinned_syst_pp[*it][i];
      }
      eylaa = fabs(theVarsBinned_stat_aa[*it][i]);
      eyhaa = eylaa;
      eylpp = fabs(theVarsBinned_stat_pp[*it][i]);
      eyhpp = eylpp;
      
      theGraphs_aa[*it]->SetPoint(i,x,yaa);
      theGraphs_aa[*it]->SetPointError(i,exl,exh,eylaa,eyhaa);
      theGraphs_syst_aa[*it]->SetPoint(i,x,yaa);
      theGraphs_syst_aa[*it]->SetPointError(i,exsyst,exsyst,eysyst_aa,eysyst_aa);
      // cout << x << " " << yaa << " " << eylaa << " " << eyhaa << " " << eysyst_aa << endl;
      
      theGraphs_pp[*it]->SetPoint(i,x,ypp);
      theGraphs_pp[*it]->SetPointError(i,exl,exh,eylpp,eyhpp);
      theGraphs_syst_pp[*it]->SetPoint(i,x,ypp);
      theGraphs_syst_pp[*it]->SetPointError(i,exsyst,exsyst,eysyst_pp,eysyst_pp);
      // cout << x << " " << ypp << " " << eylpp << " " << eyhpp << " " << eysyst_pp << endl;
      
      // theGraphs[*it]->Sort();
      // theGraphs_syst[*it]->Sort();
    }
    cnt++;
  }
  
  // plot
  plotGraphBF(theGraphs_aa, theGraphs_syst_aa, xaxis, outputDir,"pbpb");
  plotGraphBF(theGraphs_pp, theGraphs_syst_pp, xaxis, outputDir,"pp");
  plotCombGraphBF(theGraphs_pp, theGraphs_syst_pp,theGraphs_aa, theGraphs_syst_aa, xaxis, outputDir);
}

void plotGraphRAA(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high) {
  setTDRStyle();
  
  const char* ylabel = "R_{AA}";
  
  int intervals2Plot = theGraphs.size();
  
  if (intervals2Plot>1 && plot14005) return;
  
  if (plotFwdMid && intervals2Plot<=3) return;
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
    if (intervals2Plot != 1) {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0, 30);
      line = TLine(doLogPt ? 3 : 0, 1, 30, 1);
    } else {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is16004 ? 30 : 50);
      line = TLine(doLogPt ? 3 : 0,1,is16004 ? 30 : 50,1);
    }
    if (doLogPt) c1->SetLogx();
  }
  if (xaxis=="rap") {
    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
    line = TLine(0,1,2.4,1);
  }
  if (xaxis=="cent") {
    haxes = new TH1F("haxes","haxes",1,0,420);
    haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X"));
    line = TLine(0,1,420,1);
  }
  haxes->GetYaxis()->SetRangeUser(0,1.5);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : ((xaxis=="rap") ? "|y|" : "N_{part}");
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
//  line.Draw();
  
  double xshift=0.025;
//  if (xaxis=="cent") xshift=0.05;
//  TLegend *tleg = new TLegend(0.16+xshift,0.67,0.50, ((xaxis == "cent") && is16004) ? 0.89 : 0.83);
  TLegend *tleg(0x0);
  if (xaxis!="cent" && intervals2Plot == 2) tleg = new TLegend(0.44,0.50,0.76,0.62);
  else if (xaxis=="cent" && intervals2Plot == 2) tleg = new TLegend(0.19,0.16,0.51,0.28);
  else if ((xaxis=="cent" || xaxis=="rap")  && intervals2Plot == 1 && plot14005) tleg = new TLegend(0.51,0.47,0.83,0.62);
  else if (xaxis=="pt" && intervals2Plot == 1 && plot14005) tleg = new TLegend(0.19,0.49,0.51,0.64);
  else if (dononprompt && intervals2Plot == 3) tleg = new TLegend(0.56,0.47,0.88,0.62);
  else if (doprompt && intervals2Plot == 3) tleg = new TLegend(0.19,0.49,0.51,0.64);
  else if (plotFwdMid && intervals2Plot == 4) tleg = new TLegend(0.56,0.47,0.88,0.62);
  else tleg = new TLegend(0.56,0.42,0.88,0.62);
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  tleg->SetTextFont(42);
  tleg->SetTextSize(0.04);
  bool drawLegend(true);
  
  // prepare for the printing of the result tables
  const char* xname = (xaxis=="cent") ? "Centrality" : (xaxis=="pt" ? "\\pt" : "$|y|$");
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s_RAA_%s%s%s_%s.tex",outputDir.c_str(),plotPsi2S ? "Psi2S" : "JPsi",xaxis.c_str(),nameTag.c_str(), (xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr");
  string yname(plotPsi2S ? "$\\raa (\\psiP)$" : "$\\raa (\\Jpsi)$");
  inittex(texname, xname, yname);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    // if needed, draw 14-005 first
    TGraphErrors *tg_16004(0x0);
    if (plot14005) {
      TString name14005 = Form("/afs/cern.ch/work/m/miheejo/public/2015JpsiRAA5TeV/276TeVRaa/histRaa/makeRaa_%s.root",xaxis=="rap" ? "y" : xaxis.c_str());
      TString name14005_syst = Form("/afs/cern.ch/work/m/miheejo/public/2015JpsiRAA5TeV/276TeVRaa/histSyst/raaSystUncert_%s.root",xaxis=="rap" ? "y" : xaxis.c_str());
      TFile *fraa = new TFile(name14005);
      TFile *fraa_syst = new TFile(name14005_syst);
      if (fraa && fraa_syst) {
        // the main graph
        name14005 = Form("g%sJpsi",dononprompt ? "NonPr" : "Pr");
        if (thebin != anabin(0,2.4,6.5,50,0,200) && thebin != anabin(0,2.4,6.5,30,0,200)) {
          int ptlowi = thebin.ptbin().low()*10;
          if (ptlowi == ((int) thebin.ptbin().low())*10) ptlowi = thebin.ptbin().low();
          int pthighi = (thebin.ptbin().high()<30) ? thebin.ptbin().high() : 30;
          int raplowi = thebin.rapbin().low()*10;
          int raphighi = thebin.rapbin().high()*10;
          name14005 += Form("_pt%i%iy%i%i", ptlowi, pthighi, raplowi, raphighi);
        }
        TGraphErrors *tg = (TGraphErrors*) fraa->Get(name14005);
        tg_16004 = tg;
        
        // the syst graph
        name14005 = Form("g%sJpsiSyst",dononprompt ? "NonPr" : "Pr");
        if (thebin != anabin(0,2.4,6.5,50,0,200) && thebin != anabin(0,2.4,6.5,30,0,200)) {
          int ptlowi = thebin.ptbin().low()*10;
          if (ptlowi == ((int) thebin.ptbin().low())*10) ptlowi = thebin.ptbin().low();
          int pthighi = (thebin.ptbin().high()<30) ? thebin.ptbin().high() : 30;
          int raplowi = thebin.rapbin().low()*10;
          int raphighi = thebin.rapbin().high()*10;
          name14005 += Form("_pt%i%iy%i%i", ptlowi, pthighi, raplowi, raphighi);
        }
        TGraphErrors *tg_syst = (TGraphErrors*) fraa_syst->Get(name14005);
        
        if (tg && tg_syst) {
          
          // the syst box
          name14005 = Form("lumi_%s_y024_pt6530",dononprompt ? "npr" : "pr");
          if (xaxis == "pt" || xaxis == "rap") name14005 = "lumi";
          TBox *b14005 = (TBox*) fraa_syst->Get(name14005);
          double glob_14005;
          if (!b14005) glob_14005 = 0;
          else glob_14005 = b14005->GetY2() - 1;
          
          // now we have all the ingredients! Let's draw them.
          tg->SetMarkerStyle(markerstyle(cnt));
          tg->SetMarkerColor(13+cnt);
          tg->SetLineColor(13+cnt);
          tg_syst->SetLineColor(13+cnt);
          tg_syst->SetFillColorAlpha(13+cnt, 0.5);
//          if (markerstyle(cnt+4) == kFullStar) tg->SetMarkerSize(2.3);
//          else if (markerstyle(cnt+4) == kFullDiamond) tg->SetMarkerSize(2.2);
//          else if (markerstyle(cnt+4) == kFullCross) tg->SetMarkerSize(2.0);
//          else tg->SetMarkerSize(1.5);
          tg->SetMarkerSize(1.5);
          tg->SetLineWidth(tg->GetLineWidth()*2);
          
          if (xaxis=="cent") {
            // do not plot wide centrality bins
            prune(tg, tg_syst);
          }
          if (xaxis=="rap") tg_syst->GetXaxis()->SetNdivisions(306,false);
          
          tg_syst->Draw("5");
          gStyle->SetEndErrorSize(5);
          tg->Draw("P");
          
//          tleg->AddEntry(tg, "#splitline{#sqrt[]{s_{NN}} = 2.76 TeV}{#scale[0.6]{(Eur.Phys.J. C77 (2017) no.4, 252)}}", "p");
          
          // in the case where the centrality dependence is plotted: treat the PP uncertainties as global systematics
          double x(0.), dx(0.), y(0.), dy(0.);
          double rightA = 0.;
          int centminGlob(0),centmaxGlob(200);
          if (xaxis=="cent") {
            dx = 10;
            rightA = 420.;
          } else if (xaxis=="pt") {
            dx = 0.625;
            if (intervals2Plot == 1)
            {
              rightA = 50.;
              dx = 1.25;
            }
            else
            {
              rightA = 30.;
              dx = 0.65;
            }
            if (intervals2Plot == 3)
            {
              centminGlob = it->first.centbin().low();
              centmaxGlob = it->first.centbin().high();
            }
          } else if (xaxis=="rap") {
            dx = 0.06;
            rightA = 2.4;
          }
          x = rightA - (2*dx*cnt + dx);
          y = 1;

          dy = glob_14005;
          TBox *tbox = new TBox(x-dx,y-dy,x+dx,y+dy);
          tbox->SetFillColorAlpha(13+cnt, 1);
          tbox->SetLineColor(kBlack);
          tbox->Draw("l");
          TBox *tboxl = (TBox*) tbox->Clone("tboxl");
          tboxl->SetFillStyle(0);
          tboxl->Draw("l");
          
          cnt++;
        } // if (tg && tg_syst)
      } // if (fraa && fraa_syst)
      
      if (fraa) delete fraa;
      if (fraa_syst) delete fraa_syst;
    } // if (plot14005)
    
    int style = plot14005 ? cnt -1 : cnt;
    int colorI = plot14005 ? cnt -1 : cnt;
    int colorF = color(colorI)-11;
    if (intervals2Plot==2)
    {
      style = 4;
      cnt==0 ? colorI = 9 : colorI =4;
      colorF = color(colorI)-11;
    }
    if (intervals2Plot==3)
    {
      style = 0;
      colorI = cnt+6;
      colorF = color(colorI)-10;
    }
    if (intervals2Plot>3)
    {
      style = cnt+1;
      colorI = cnt+1;
      colorF = color(colorI)-11;
    }
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    if (xaxis=="cent") {
      // do not plot wide centrality bins
      prune(tg, tg_syst);
    }
    if (is16004 && plotPsi2S)
    {
      Int_t nPoints = tg->GetN();
      double x, y, exl, exh;
//      double x_sys, y_sys;
      double liml, limh;
      for (Int_t i = 0 ; i < nPoints ; i++)
      {
        tg->GetPoint(i,x,y);
        if (y == 0)
        {
          exl = tg->GetErrorXlow(i);
          exh = tg->GetErrorXhigh(i);
          
          liml = tg->GetErrorYhigh(i);
//          tg_syst->GetPoint(i,x_sys,y_sys);
          limh = tg_syst->GetErrorYhigh(i);;
          
          double dx(0.);
          if (xaxis=="cent")
          {
            x = HI::findNpartAverage(2.*(x-exl),2.*(x+exh));
            dx =10;
          }
          if (xaxis=="pt") dx = 0.5;
          
          drawArrow(x, liml, limh, dx, color(cnt));
          tg->SetPoint(i,-1000.,0);
          tg->SetPointError(i,0,0,0,0);
          tg_syst->SetPoint(i,-1000.,0);
          tg_syst->SetPointError(i,0,0,0,0);
            
          //delete point in res and syst
        }
      }
    }
    
    bool plot = true;
    if (plotFwdMid && (cnt!=0 && cnt!=3)) plot = false;
    
    // in the case where the centrality dependence is plotted: treat the PP uncertainties as global systematics
    // if (xaxis == "cent") {
    if (plot)
    {
      double x(0.), dx(0.), y(0.), dy_low(0.), dy_high(0.);
      double rightA = 0.;
      int centminGlob(0),centmaxGlob(200);
      if (xaxis=="cent") {
        dx = 10;
        rightA = 420.;
      } else if (xaxis=="pt") {
        dx = 0.625;
        if (intervals2Plot == 1)
        {
          rightA = 50.;
          dx = 1.25;
        }
        else
        {
          rightA = 30.;
          dx = 0.65;
        }
        if (intervals2Plot == 3)
        {
          centminGlob = it->first.centbin().low();
          centmaxGlob = it->first.centbin().high();
        }
      } else if (xaxis=="rap") {
        dx = 0.06;
        rightA = 2.4;
      }
      x = rightA - (2*dx*cnt + dx);
      if (plotFwdMid && cnt>0) x = rightA - (2*dx + dx);
      y = 1;
      
      anabin thebinglb(it->first.rapbin().low(),
                       it->first.rapbin().high(),
                       it->first.ptbin().low(),
                       it->first.ptbin().high(),
                       centminGlob,centmaxGlob);
      thebinglb.print();
      dy_low = gsyst_low[thebinglb].value;
      dy_high = gsyst_high[thebinglb].value;
      cout << "global syst: " << "+" << dy_high << " -" << dy_low << endl;
      TBox *tbox = new TBox(x-dx,y-dy_low,x+dx,y+dy_high);
      tbox->SetFillColorAlpha(colorF, 1);
      tbox->SetLineColor(color(colorI));
      tbox->Draw("l");
      TBox *tboxl = (TBox*) tbox->Clone("tboxl");
      tboxl->SetFillStyle(0);
      tboxl->Draw("l");
    }
    // }
    
    // Plot graphs after uncertainties to avoid overlap

    if (plot)
    {
      tg_syst->Draw("5");
      gStyle->SetEndErrorSize(5);
      tg->Draw("P");
      // tg->Draw("[]");
    }
  
    TLatex tl;
    double tlx = 0.25; //0.92;
    double tly = 0.80; //0.69;
    tl.SetTextFont(42); // consistent font for symbol and plain text

    if (plot)
    {
      if (plot14005)
      {
        tleg->AddEntry(tg, "#sqrt[]{s_{NN}} = 5.02 TeV", "p");
//        tleg->AddEntry(tg_16004, "#splitline{#sqrt[]{s_{NN}} = 2.76 TeV}{#scale[0.6]{Eur.Phys.J. C77 (2017) 252}}", "p");
        tleg->AddEntry(tg_16004, "#sqrt[]{s_{NN}} = 2.76 TeV", "p");
        
        
        double xp = 0.59;
        double yp = 0.45;
        if (xaxis == "pt")
        {
          xp = 0.27;
          yp = 0.47;
        }
        
        TLatex *tex = new TLatex(xp,yp,"#scale[0.6]{Eur.Phys.J. C77 (2017) 252}");
        tex->SetNDC();
        tex->SetTextSize(0.04);
        tex->SetTextFont(42);
        tex->SetLineWidth(2);
        tex->Draw();
      }
      else
      {
        TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
        if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
        TString ptlabel = Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
        TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
        
        if (xaxis == "pt")
        {
          if (is16004) tleg->AddEntry(tg, raplabel, "p");
          else
          {
            if (intervals2Plot > 3) tleg->AddEntry(tg, raplabel, "p");
            else if (intervals2Plot == 3) tleg->AddEntry(tg, Form("Cent. %s",centlabel.Data()), "p");
            else if (intervals2Plot == 1) drawLegend = false;
            else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",raplabel.Data(),centlabel.Data()));
            //else tleg->AddEntry(tg, Form("%s, Cent. %s",raplabel.Data(),centlabel.Data()), "p");
          }
        }
        if (xaxis == "cent")
        {
          if (is16004) tleg->AddEntry(tg, Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()), "p");
          else
          {
            if (intervals2Plot > 3 ) tleg->AddEntry(tg, raplabel, "p");
            else if (intervals2Plot == 2) tleg->AddEntry(tg, ptlabel, "p");
            else if (intervals2Plot == 1) drawLegend = false;
            else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()));
            //else tleg->AddEntry(tg, Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()), "p");
          }
        }
        if (xaxis == "rap")
        {
          if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
          else if (intervals2Plot == 1) drawLegend = false;
          else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
          //else tleg->AddEntry(tg, Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()), "p");
        }
      }
    }
    
    // print tex
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
    if (xaxis == "cent" || xaxis == "rap") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";
    
    addline(texname,oss.str());
    printGraph(tg, tg_syst, texname);
    
    // for the centrality dependence: we want Npart plotted, not the centrality
    if (xaxis == "cent") {
      centrality2npart(tg, false, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
      centrality2npart(tg_syst, true, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
    }
    
    cnt++;
    it_syst++;
  }
  
  
  if (drawLegend) tleg->Draw();
  line.Draw();
  
  if(xaxis!="cent" && intervals2Plot > 3)
  {
    //      TLatex *tex = new TLatex(0.21,0.86,"Cent. 0-100%");
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 2 && !is16004)
  {
    TLatex *tex = new TLatex(0.2,0.78,"1.8 < |y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot > 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"6.5 < p_{T} < 50 GeV/c");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="pt" && is16004)
  {
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="rap" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  
  TLatex tl;
//  double tlx = 0.54; //0.92;
  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;
//  tl.SetTextAlign(32); // right adjusted
  tl.SetTextFont(42); // consistent font for symbol and plain text
  tl.SetTextSize(0.057); 
  if (doprompt) tl.DrawLatexNDC(tlx,tly,plotPsi2S ? "Prompt #psi(2S)" : "Prompt J/#psi");
  //if (doprompt) tl.DrawLatexNDC(tlx-0.08,tly,plotPsi2S ? "Prompt #psi(2S)" : "Prompt J/#psi");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi from b hadrons");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi (b hadrons)");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi #leftarrow B");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"B #rightarrow J/#psi");
  tl.SetTextSize(0.046);
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  // CMS_lumi( (TPad*) gPad, 103, iPos, "" );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_RAA_%s%s%s_%s.root",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_RAA_%s%s%s_%s.png",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_RAA_%s%s%s_%s.pdf",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  if (plotFwdMid)
  {
    gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_RAA_%s%s%s_2Rapranges_%s.root",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_RAA_%s%s%s_2Rapranges_%s.png",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_RAA_%s%s%s_2Rapranges_%s.pdf",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  }
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}

void plotGraphXS(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high, const char* collName) {
  setTDRStyle();
  
  const char* ylabel = "";
  TString scollName(collName);
  bool isPP = scollName.Contains("pp");
  
  int intervals2Plot = theGraphs.size();
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
    if (scollName.Contains("pp")) ylabel = Form("#bf{#it{#Beta}} #times d%s#sigma/dp_{T}%s (nb / GeV/#it{c})",intervals2Plot>1 ? "^{2}" : "",intervals2Plot>1 ? "dy" : "");
    else ylabel = Form("#bf{#it{#Beta}} #times (1/T_{AA}) #times d%sN/dp_{T}%s (nb / GeV/#it{c})",intervals2Plot>1 ? "^{2}" : "",intervals2Plot>1 ? "dy" : "");
//    if (scollName.Contains("pp")) ylabel = "#bf{#it{#Beta}} #times d#sigma/dp_{T} (nb / GeV/#it{c})";
//    else ylabel = "#bf{#it{#Beta}} #times (1/T_{AA}) #times dN/dp_{T} (nb / GeV/#it{c})";

    if (intervals2Plot != 1) {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0, 30);
      line = TLine(doLogPt ? 3 : 0, 1, 30, 1);
    } else {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is16004 ? 30 : 50);
      line = TLine(doLogPt ? 3 : 0,1,is16004 ? 30 : 50,1);
    }

    if (doLogPt) c1->SetLogx();
    haxes->GetYaxis()->SetRangeUser(0.0001,5E2);
    
    c1->SetLogy();
  }
  if (xaxis=="rap") {
    //    ylabel = "B #cdot #frac{d#sigma}{d#it{y}} (#mu b)";
    if (scollName.Contains("pp")) ylabel = "#bf{#it{#Beta}} #times d#sigma/d#it{y} (nb)";
    else ylabel = "#bf{#it{#Beta}} #times (1/T_{AA}) #times d#sigma/d#it{y} (nb)";
//    if (scollName.Contains("pp")) ylabel = "#bf{#it{#Beta}} #times d#sigma/d#it{y} (nb)";
//    else ylabel = "#bf{#it{#Beta}} #times (1/T_{AA}) #times d#sigma/d#it{y} (nb)";
    
//    c1->SetLogy();

    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetYaxis()->SetRangeUser(0.001, 1);
//    haxes->GetYaxis()->SetRangeUser(0.001,scollName.Contains("pp") ? 1 : 0.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
  }
  if (xaxis=="cent") {
    cout << "ERROR centrality plots are not implemented for cross section" << endl;
    return;
  }
  
  haxes->GetYaxis()->SetTitle(ylabel);
  haxes->GetYaxis()->SetTitleSize(0.055);
  haxes->GetYaxis()->SetTitleOffset(1.35);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "|y|";
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
  
  double xshift=0.025;
  //  if (xaxis=="cent") xshift=0.05;
  TLegend *tleg = new TLegend(0.20,0.17,0.50,0.33);
  //TLegend *tleg = new TLegend(0.16+xshift,0.67,0.50,0.83);
  tleg->SetFillStyle(0);
  tleg->SetBorderSize(0);
  tleg->SetTextSize(0.037);
  tleg->SetTextFont(42);
  
  // prepare for the printing of the result tables
  const char* xname = xaxis=="pt" ? "\\pt" : "$|y|$";
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s_XS_%s_%s%s%s_%s.tex",outputDir.c_str(),"JPsi", scollName.Data(),xaxis.c_str(),nameTag.c_str(), (xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr");
  string yname(Form("$\\dsd%s_%s (\\Jpsi)$",xaxis.c_str(),scollName.Data()));
  inittex(texname, xname, yname);
  
  int cnt=0;
  double dy_low(0.), dy_high(0.); //For global uncertainties
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;
    if (intervals2Plot==2)
    {
      style = 4;
      cnt==0 ? colorI = 9 : colorI =4;
      colorF = color(colorI)-11;
    }
    if (intervals2Plot==3)
    {
      style = 0;
      colorI = cnt+6;
      colorF = color(colorI)-10;
    }
    if (intervals2Plot>3)
    {
      style = cnt+1;
      colorI = cnt+1;
      colorF = color(colorI)-11;
    }
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");
    
    TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
    if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
    TString ptlabel = Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
    TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
    
    TLatex tl;
    double tlx = 0.92;
    double tly = isPP ? 0.59 : 0.57; //0.69;
    tl.SetTextAlign(32); // right adjusted
    tl.SetTextFont(42); // consistent font for symbol and plain text
    tl.SetTextSize(0.037); 

    if (xaxis == "pt")
    {
      if (intervals2Plot > 3) tleg->AddEntry(tg, raplabel, "p");
      else if (intervals2Plot == 3) tleg->AddEntry(tg, Form("Cent. %s",centlabel.Data()), "p");
      else if (isPP) tl.DrawLatexNDC(tlx,tly,raplabel.Data());
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",raplabel.Data(),centlabel.Data()));
    }
    if (xaxis == "cent")
    {
      if (intervals2Plot > 3 ) tleg->AddEntry(tg, raplabel, "p");
      else if (intervals2Plot == 2) tleg->AddEntry(tg, ptlabel, "p");
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()));
    }
    if (xaxis == "rap")
    {
      if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
      else if (isPP) tl.DrawLatexNDC(tlx,tly,ptlabel.Data());
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
    }
    
    // print tex
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
    
    addline(texname,oss.str());
    printGraph(tg, tg_syst, texname);
    
    //  global systematics
    int centminGlob(0),centmaxGlob(200);
    if (xaxis=="pt" && intervals2Plot == 3)
    {
      centminGlob = it->first.centbin().low();
      centmaxGlob = it->first.centbin().high();
    }
    
    anabin thebinglb(it->first.rapbin().low(),
                     it->first.rapbin().high(),
                     it->first.ptbin().low(),
                     it->first.ptbin().high(),
                     centminGlob,centmaxGlob);
    thebinglb.print();
    dy_low = gsyst_low[thebinglb].value;
    dy_high = gsyst_high[thebinglb].value;
    cout << "global syst: " << "+" << dy_high << " -" << dy_low << endl;
    
    
    cnt++;
    it_syst++;
  }
  
  tleg->Draw();
  
  // Draw global syst
  TLatex *texG(0x0);
  double xG = 0.92;
  double yG = 0.59;
  if (intervals2Plot == 1) yG = 0.67;
  if (isPP) texG = new TLatex(xG,yG,Form("Global unc. %1.1f %%",dy_high*100.));
  else texG = new TLatex(xG,yG,Form("Global unc. #splitline{+%1.1f %%}{-%1.1f %%}",dy_high*100.,dy_low*100.));
  texG->SetTextAlign(32);
  texG->SetNDC();
  texG->SetTextSize(0.037);
  texG->SetTextFont(42);
  texG->SetLineWidth(2);
  texG->Draw();
  
  if(xaxis=="pt" && intervals2Plot == 3)
  {
    TLatex *tex = new TLatex(0.22,0.35,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }

  
  
  TLatex tl;
//  double tlx = 0.92; // position for the default font + right adjust
//  double tly = 0.69;
//  tl.SetTextAlign(32); // right adjusted
//  tl.SetTextFont(42); // consistent font for symbol and plain text
//  if (doprompt) tl.DrawLatexNDC(tlx,tly,"Prompt J/#psi");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  
  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;
  //  tl.SetTextAlign(32); // right adjusted
  tl.SetTextFont(42); // consistent font for symbol and plain text
  tl.SetTextSize(0.055);
  if (doprompt) tl.DrawLatexNDC(tlx,tly,"Prompt J/#psi");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi from b hadrons");
  tl.SetTextSize(0.046);
  
  
//  TLatex t2;
//  t2.SetTextAlign(32); // right adjusted
//  tl.DrawLatexNDC(tlx,tly+0.1,Form("Global unc. : +%.1%f-%.1%f",dy_high,dy_low));
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  // CMS_lumi( (TPad*) gPad, 103, iPos, "" );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_XS_%s_%s%s%s_%s.root",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_XS_%s_%s%s%s_%s.png",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_XS_%s_%s%s%s_%s.pdf",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}

void plotCombGraphXS(map<anabin, TGraphAsymmErrors*> theGraphs_pp, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pp, map<anabin, TGraphAsymmErrors*> theGraphs_pbpb, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pbpb, string xaxis, string outputDir, map<anabin, syst> gsyst_low_pp, map<anabin, syst> gsyst_high_pp, map<anabin, syst> gsyst_low_pbpb, map<anabin, syst> gsyst_high_pbpb, const char* collName) {
  setTDRStyle();
  
  const char* ylabel = "";
  TString scollName(collName);
  
  int intervals2Plot = theGraphs_pp.size();
  
  if (intervals2Plot > 1)
  {
    cout << "[WARNING] combined plotting of XS not supported for more than one range" << endl;
    return;
  }
  if (xaxis=="cent") {
    cout << "[WARNING] combined plotting of XS not supported vs centrality" << endl;
    return;
  }
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
//    if (scollName.Contains("pp")) ylabel = Form("#bf{#it{#Beta}} #times d%s#sigma/dp_{T}%s (nb / GeV/#it{c})",intervals2Plot>1 ? "^{2}" : "",intervals2Plot>1 ? "dy" : "");
//    else ylabel = Form("#bf{#it{#Beta}} #times (1/T_{AA}) #times d%sN/dp_{T}%s (nb / GeV/#it{c})",intervals2Plot>1 ? "^{2}" : "",intervals2Plot>1 ? "dy" : "");
//    
    ylabel = "d#sigma/dp_{T} (nb / GeV/#it{c})";
    
    haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is16004 ? 30 : 50);
    line = TLine(doLogPt ? 3 : 0,1,is16004 ? 30 : 50,1);
    
    if (doLogPt) c1->SetLogx();
    haxes->GetYaxis()->SetRangeUser(0.0001,5E2);
    
    c1->SetLogy();
  }
  if (xaxis=="rap") {
    //    ylabel = "B #cdot #frac{d#sigma}{d#it{y}} (#mu b)";
//    if (scollName.Contains("pp")) ylabel = "#bf{#it{#Beta}} #times d#sigma/d#it{y} (nb)";
//    else ylabel = "#bf{#it{#Beta}} #times (1/T_{AA}) #times d#sigma/d#it{y} (nb)";
    ylabel = "d#sigma/dy (nb / GeV/#it{c})";
    //    if (scollName.Contains("pp")) ylabel = "#bf{#it{#Beta}} #times d#sigma/d#it{y} (nb)";
    //    else ylabel = "#bf{#it{#Beta}} #times (1/T_{AA}) #times d#sigma/d#it{y} (nb)";
    
    //    c1->SetLogy();
    
    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetYaxis()->SetRangeUser(0.001, 1.2);
    //    haxes->GetYaxis()->SetRangeUser(0.001,scollName.Contains("pp") ? 1 : 0.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
  }
  
  haxes->GetYaxis()->SetTitle(ylabel);
  haxes->GetYaxis()->SetTitleSize(0.055);
  haxes->GetYaxis()->SetTitleOffset(1.35);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "|y|";
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
  
  double xshift=0.025;
  //  if (xaxis=="cent") xshift=0.05;
  TLegend *tleg(0x0);
  if (xaxis=="rap")
  {
//    if (doprompt) tleg = new TLegend(0.16,0.55,0.46,0.71);
//    if (dononprompt) tleg = new TLegend(0.16,0.45,0.46,0.61);
    tleg = new TLegend(0.16,0.55,0.46,0.71);
  }
  else tleg = new TLegend(0.20,0.17,0.50,0.33);
  
  //TLegend *tleg = new TLegend(0.16+xshift,0.67,0.50,0.83);
  tleg->SetFillStyle(0);
  tleg->SetBorderSize(0);
  tleg->SetTextSize(0.037);
  tleg->SetTextFont(42);
  
  int cnt=0;
  double dy_low_pp(0.), dy_high_pp(0.); //For global uncertainties
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs_pp.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst_pp.begin();
  for (; it!=theGraphs_pp.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = kOpenCircle;
    int colorI = kBlue;
    int colorF = kBlue-11;//kWhite;
    
    tg->SetMarkerStyle(style);
    tg->SetMarkerColor(colorI);
    tg->SetLineColor(colorI);
    tg_syst->SetLineColor(colorI);
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");
    
    tleg->AddEntry(tg, "pp", "p");
    
    //  global systematics
    anabin thebinglb(it->first.rapbin().low(),
                     it->first.rapbin().high(),
                     it->first.ptbin().low(),
                     it->first.ptbin().high(),
                     0,200);
    thebinglb.print();
    dy_low_pp = gsyst_low_pp[thebinglb].value;
    dy_high_pp = gsyst_high_pp[thebinglb].value;
    cout << "pp global syst: " << "+" << dy_high_pp << " -" << dy_low_pp << endl;
    cnt++;
    it_syst++;
  }
  
  cnt=0;
  double dy_low_pbpb(0.), dy_high_pbpb(0.); //For global uncertainties
  it=theGraphs_pbpb.begin();
  it_syst=theGraphs_syst_pbpb.begin();
  for (; it!=theGraphs_pbpb.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");
    
    tleg->AddEntry(tg, "PbPb, Cent. 0-100%", "p");
    
    //  global systematics
    anabin thebinglb(it->first.rapbin().low(),
                     it->first.rapbin().high(),
                     it->first.ptbin().low(),
                     it->first.ptbin().high(),
                     0,200);
    thebinglb.print();
    dy_low_pbpb = gsyst_low_pbpb[thebinglb].value;
    dy_high_pbpb = gsyst_high_pbpb[thebinglb].value;
    cout << "pbpb global syst: " << "+" << dy_high_pbpb << " -" << dy_low_pbpb << endl;
    cnt++;
    it_syst++;
  }
  
  tleg->Draw();
  
  // Draw global syst
  TLatex *texG(0x0);
  double xG = 0.92;
  double yG = 0.59;
  if (intervals2Plot == 1) yG = 0.67;
  texG = new TLatex(xG-0.055,yG,Form("pp global unc. %1.1f %%",dy_high_pp*100.));
  texG->SetTextAlign(32);
  texG->SetNDC();
  texG->SetTextSize(0.034);
  texG->SetTextFont(42);
  texG->SetLineWidth(2);
  texG->Draw();
  
  TLatex *texGpbpb(0x0);
  texGpbpb = new TLatex(xG,yG-0.07,Form("PbPb global unc. #splitline{+%1.1f %%}{#minus%1.1f %%}",dy_high_pbpb*100.,dy_low_pbpb*100.));
  texGpbpb->SetTextAlign(32);
  texGpbpb->SetNDC();
  texGpbpb->SetTextSize(0.034);
  texGpbpb->SetTextFont(42);
  texGpbpb->SetLineWidth(2);
  texGpbpb->Draw();
  
  
  TLatex tl;
  TLatex t2;
  //  double tlx = 0.92; // position for the default font + right adjust
  //  double tly = 0.69;
  //  tl.SetTextAlign(32); // right adjusted
  //  tl.SetTextFont(42); // consistent font for symbol and plain text
  //  if (doprompt) tl.DrawLatexNDC(tlx,tly,"Prompt J/#psi");
  //  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  
  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;
  //  tl.SetTextAlign(32); // right adjusted
  tl.SetTextFont(42); // consistent font for symbol and plain text
  tl.SetTextSize(0.055);
  t2.SetTextFont(42); // consistent font for symbol and plain text
  t2.SetTextSize(0.055);
  if (doprompt) tl.DrawLatexNDC(tlx,tly,"Prompt J/#psi");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi from b hadrons");
  t2.DrawLatexNDC(tlx,tly-0.05,"J/#psi #rightarrow #mu^{+}#mu^{-}");
  tl.SetTextSize(0.046);
  t2.SetTextSize(0.046);
  
  if(xaxis=="pt")
  {
    TLatex *tex = new TLatex(tlx,tly-0.1,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="rap")
  {
    TLatex *tex = new TLatex(tlx,tly-0.1,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_XS_%s_%s%s%s_%s.root",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_XS_%s_%s%s%s_%s.png",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_XS_%s_%s%s%s_%s.pdf",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  delete tleg;
  delete haxes;
  delete c1;
}

void plotGraphBF(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, const char* collName) {
  setTDRStyle();
  
  const char* ylabel = "Nonprompt J/#psi fraction";
  TString scollName(collName);
  bool isPP = scollName.Contains("pp");
  
  int intervals2Plot = theGraphs.size();
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
    if (intervals2Plot != 1) {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0, 30);
      line = TLine(doLogPt ? 3 : 0, 1, 30, 1);
    } else {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is16004 ? 30 : 50);
      line = TLine(doLogPt ? 3 : 0,1,is16004 ? 30 : 50,1);
    }
    if (doLogPt) c1->SetLogx();
  }
  if (xaxis=="rap") {
    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
  }
  if (xaxis=="cent") {
    haxes = new TH1F("haxes","haxes",1,0,420);
    haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X"));
  }
  haxes->GetYaxis()->SetRangeUser(0,0.9);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : ((xaxis=="rap") ? "|y|" : "N_{part}");
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
  
  double xshift=0.025;
  //  if (xaxis=="cent") xshift=0.05;
  TLegend *tleg;
  if(xaxis=="cent" && intervals2Plot == 2) tleg = new TLegend(0.21,0.20,0.50,0.36);
  else tleg = new TLegend(0.16+xshift,0.67,0.50,0.83);
  tleg->SetFillStyle(0);
  tleg->SetBorderSize(0);
  tleg->SetTextSize(0.037);
  tleg->SetTextFont(42);
  
  // prepare for the printing of the result tables
  const char* xname = (xaxis=="cent") ? "Centrality" : (xaxis=="pt" ? "\\pt" : "$|y|$");
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s_BF_%s_%s%s%s_%s.tex",outputDir.c_str(),"JPsi", scollName.Data(),xaxis.c_str(),nameTag.c_str(), (xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr");
  string yname("$\\raa (\\Jpsi)$");
  inittex(texname, xname, yname);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;
    if (intervals2Plot==2)
    {
      style = 4;
      cnt==0 ? colorI = 9 : colorI =4;
      colorF = color(colorI)-11;
    }
    if (intervals2Plot==3)
    {
      style = 0;
      colorI = cnt+6;
      colorF = color(colorI)-10;
    }
    if (intervals2Plot>3)
    {
      style = cnt+1;
      colorI = cnt+1;
      colorF = color(colorI)-11;
    }
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    if (xaxis=="cent") {
      // do not plot wide centrality bins
      prune(tg, tg_syst);
    }
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");

    TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
    if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
    TString ptlabel = Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
    TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
    
    TLatex tl;
    double tlx = 0.20;
    double tly = 0.80; //0.69;
    tl.SetTextSize(0.037);
    tl.SetTextFont(42); // consistent font for symbol and plain text

    if (xaxis == "pt")
    {
      if (intervals2Plot > 3) tleg->AddEntry(tg, raplabel, "p");
      else if (intervals2Plot == 3) tleg->AddEntry(tg, Form("Cent. %s",centlabel.Data()), "p");
      else if (isPP) tl.DrawLatexNDC(tlx,tly,raplabel.Data());
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",raplabel.Data(),centlabel.Data()));
    }
    if (xaxis == "cent")
    {
      if (intervals2Plot > 3 ) tleg->AddEntry(tg, raplabel, "p");
      else if (intervals2Plot == 2) tleg->AddEntry(tg, ptlabel, "p");
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()));
    }
    if (xaxis == "rap")
    {
      if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
      else if (isPP) tl.DrawLatexNDC(tlx,tly, ptlabel.Data());
      else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
    }
    
    // print tex
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
    if (xaxis == "cent" || xaxis == "rap") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";
    
    addline(texname,oss.str());
    printGraph(tg, tg_syst, texname);
    
    // for the centrality dependence: we want Npart plotted, not the centrality
    if (xaxis == "cent") {
      centrality2npart(tg, false, (150./1.6)*it->first.rapbin().low());
      centrality2npart(tg_syst, true, (150./1.6)*it->first.rapbin().low());
    }
    
    cnt++;
    it_syst++;
  }
  
  
  tleg->Draw();
  
  if(xaxis!="cent" && intervals2Plot > 3 && !isPP)
  {
    TLatex *tex = new TLatex(0.20,0.86,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 2)
  {
    TLatex *tex = new TLatex(0.20,0.86,"1.8 < |y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot > 3)
  {
    TLatex *tex = new TLatex(0.20,0.86,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 3)
  {
    TLatex *tex = new TLatex(0.20,0.86,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_BF_%s_%s%s%s_%s.root",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_BF_%s_%s%s%s_%s.png",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_BF_%s_%s%s%s_%s.pdf",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}

void plotCombGraphBF(map<anabin, TGraphAsymmErrors*> theGraphs_pp, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pp, map<anabin, TGraphAsymmErrors*> theGraphs_pbpb, map<anabin, TGraphAsymmErrors*> theGraphs_syst_pbpb, string xaxis, string outputDir, const char* collName)
{
  setTDRStyle();
  
  const char* ylabel = "Nonprompt J/#psi fraction";
  TString scollName(collName);
  
  int intervals2Plot = theGraphs_pp.size();
  
  if (intervals2Plot > 1)
  {
    cout << "[WARNING] combined plotting of BF not supported for more than one range" << endl;
    return;
  }
  if (xaxis=="cent") {
    cout << "[WARNING] combined plotting of BF not supported vs centrality" << endl;
    return;
  }
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
    if (intervals2Plot != 1) {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0, 30);
      line = TLine(doLogPt ? 3 : 0, 1, 30, 1);
    } else {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is16004 ? 30 : 50);
      line = TLine(doLogPt ? 3 : 0,1,is16004 ? 30 : 50,1);
    }
    if (doLogPt) c1->SetLogx();
  }
  if (xaxis=="rap") {
    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
  }
  
  haxes->GetYaxis()->SetRangeUser(0,0.9);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : "|y|";
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
  
  double xshift=0.025;
  //  if (xaxis=="cent") xshift=0.05;
  TLegend *tleg;
  tleg = new TLegend(0.16+xshift,0.67,0.50,0.83);
  tleg->SetFillStyle(0);
  tleg->SetBorderSize(0);
  tleg->SetTextSize(0.037);
  tleg->SetTextFont(42);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs_pp.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst_pp.begin();
  for (; it!=theGraphs_pp.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = kOpenCircle;
    int colorI = kBlue;
    int colorF = kBlue-11;//kWhite;
    
    tg->SetMarkerStyle(style);
    tg->SetMarkerColor(colorI);
    tg->SetLineColor(colorI);
    tg_syst->SetLineColor(colorI);
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");
    
    tleg->AddEntry(tg, "pp", "p");
    
    cnt++;
    it_syst++;
  }
  
  cnt=0;
  it=theGraphs_pbpb.begin();
  it_syst=theGraphs_syst_pbpb.begin();
  for (; it!=theGraphs_pbpb.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    tg_syst->Draw("5");
    gStyle->SetEndErrorSize(5);
    tg->Draw("P");
    
    tleg->AddEntry(tg, "PbPb, Cent. 0-100%", "p");
    
    cnt++;
    it_syst++;
  }
  
  
  tleg->Draw();
  
  if(xaxis=="pt")
  {
    TLatex *tex = new TLatex(0.20,0.86,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="rap")
  {
    TLatex *tex = new TLatex(0.20,0.86,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.037);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_BF_%s_%s%s%s_%s.root",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_BF_%s_%s%s%s_%s.png",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_BF_%s_%s%s%s_%s.pdf",outputDir.c_str(), "JPsi", scollName.Data(), xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  delete tleg;
  delete haxes;
  delete c1;
  
}

void centrality2npart(TGraphAsymmErrors* tg, bool issyst, double xshift) {
  int n = tg->GetN();
  for (int i=0; i<n; i++) {
    double x, y, exl, exh, exlC, exhC, eyl, eyh;
    x = tg->GetX()[i];
    if (x>0)
    {
      y = tg->GetY()[i];
      exlC = tg->GetErrorXlow(i);
      exhC = tg->GetErrorXhigh(i);
      eyl = tg->GetErrorYlow(i);
      eyh = tg->GetErrorYhigh(i);
      if (!issyst) {
        exl = HI::findNpartSyst_low(2.*(x-exlC),2.*(x+exhC));//0.;
        exh = HI::findNpartSyst_high(2.*(x-exlC),2.*(x+exhC));//0.;
      } else {
        exl = 5;
        exh = exl;
      }
      x = HI::findNpartAverage(2.*(x-exlC),2.*(x+exhC));
      tg->SetPoint(i,x+xshift,y);
      tg->SetPointError(i,exl,exh,eyl,eyh);
    }
    else
    {
      tg->SetPoint(i,-1000,0);
      tg->SetPointError(i,0,0,0,0);
    }
  }
}

int color(int i) {
  if (i==0) return kRed+2;
  else if (i==1) return kBlue+2;
  else if (i==2) return kMagenta+2;
  else if (i==3) return kCyan+2;
  else if (i==4) return kGreen+2;
  else if (i==5) return kOrange+2;
  else if (i==6) return kRed+1;
  else if (i==7) return kYellow+1;
  else if (i==8) return kAzure+1;
  else if (i==9) return kBlack;
  else return kBlack;
}

int markerstyle(int i) {
  if (i==0) return kFullSquare;
  else if (i==1) return kFullCircle;
  else if (i==2) return kFullStar;
  else if (i==3) return kFullCross;
  else if (i==4) return kFullDiamond;
  else if (i==5) return kOpenSquare;
  else if (i==6) return kOpenCircle;
  else if (i==7) return kOpenStar;
  else if (i==8) return kFullTriangleDown;
  else return kOpenCross;
}

void setOptions(bool adoprompt, bool adononprompt, bool aplotFwdMid, bool ais16004, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, bool aincludeEffSyst, bool  aexcludeNonFitSyst, string anameTag_base) {
  doprompt = adoprompt;
  dononprompt = adononprompt;
  plotFwdMid = aplotFwdMid;
  is16004 = ais16004;
  plotPsi2S = aplotPsi2S;
  plot14005 = aplot14005;
  applyEff = aapplyEff;
  applyAcc = aapplyAcc;
  doLogPt = adoLogPt;
  includeEffSyst = aincludeEffSyst;
  excludeNonFitSyst = aexcludeNonFitSyst;
  
  nameTag_base = anameTag_base;
  
  if (plotPsi2S && !is16004)
  {
    cout << "[ERROR] Options not set: you want to plot Psi2 but bins are not 16004 ones" << endl;
    printOptions();
    return;
  }
  if (doprompt) nameTag_base += "_prompt";
  if (dononprompt) nameTag_base += "_nonprompt";
  
//  if (aplotPsi2S) nameTag_base += "_Psi2S";
//  else nameTag_base += "_JPsi";
  
  if (plotFwdMid) nameTag_base += "_2RapRanges";
    
  if (is16004)  nameTag_base += "_16004";
  if (plot14005) nameTag_base += "_14005";
}

void printOptions() {
  cout <<
  "doprompt = " << doprompt << ", " <<
  "dononprompt = " << dononprompt << ", " <<
  "is16004 = " << is16004 << ", " <<
  "plotPsi2S = " << plotPsi2S << ", " <<
  "plot14005 = " << plot14005 << ", " <<
  "applyEff = " << applyEff << ", " <<
  "applyAcc = " << applyAcc << ", " <<
  "doLogPt = " << doLogPt << ", " <<
  "includeEffSyst = " << includeEffSyst << ", " <<
  "excludeNonFitSyst = " << excludeNonFitSyst << ", " <<
  "nameTag_base = \"" << nameTag_base << "\"" <<
  endl;
}

map<anabin, raa_input > readResults(const char* resultsFile)
{
  map<anabin, raa_input> ans;
  raa_input theresult;
  
  ifstream file(resultsFile);
  if (!(file.good())) return ans;
  
  string resultname; getline(file,resultname);
  
  cout << "[INFO] Reading results : " << resultname.c_str() << endl;
  
  string line;
  double rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;
  
  while (file.good()) {
    getline(file,line);
    if (line.size()==0) break;
    TString tline(line.c_str());
    TString t; Int_t from = 0, cnt=0;
    while (tline.Tokenize(t, from , ",")) {
      t.Strip(TString::kBoth,' ');
      value = atof(t.Data());
      if (cnt==0) rapmin = atof(t.Data());
      else if (cnt==1) rapmax = value;
      else if (cnt==2) ptmin = value;
      else if (cnt==3) ptmax = value;
      else if (cnt==4) centmin = value;
      else if (cnt==5) centmax = value;
      else if (cnt==6) theresult.naa = value;
      else if (cnt==7) theresult.dnaa_stat = value;
      else if (cnt==8) theresult.systaa = value;
      else if (cnt>8) {
        cout << "Warning, too many fields, I'll take the last one." << endl;
        continue;
      }
      cnt++;
    }
    anabin thebin(rapmin, rapmax, ptmin, ptmax, centmin, centmax);
    ans[thebin] = theresult;
  }
  
  file.close();
  
  return ans;
}

void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color) {
  TArrow *arrow = new TArrow(x,yhigh,x,ylow<=0. ? 0.01 : ylow,0.03,ylow<=0. ? ">" : "<>");
  arrow->SetLineColor(color);
  arrow->Draw();
  if (ylow<=0.) {
    TLine *line = new TLine(x-dx,yhigh,x+dx,yhigh);
    line->SetLineColor(color);
    line->Draw();
  }
}
