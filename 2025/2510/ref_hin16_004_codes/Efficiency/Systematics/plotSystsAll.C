///// IMPORTANT NOTE /////
// 1) Copy the same file from Fitter/plotVars.C, Fitter/Systematics/syst.C and Fitter/Systematics/syst.h.
// 2) Change all headers' relative location
// 3) Change input file names to the efficiency csv files in plotSystsAll.C and syst.C

#include "plotVars.C"
#include "syst.C"

#include "TString.h"
#include <sstream>

using namespace std;

TString spoiname(""); // It can be NJpsi, BJpsi, NJpsi_prompt or NJpsi_nonprompt
bool bins16004=false; // If false use 16-025 bins, if true use 16-004 bins
TString nameTag("16025"); // It can be 16025 or 16004 or...

// DECLARATIONS
void plotSysts(anabin thebin, string xaxis, string collTag, bool plotEffSyst=false, bool plotSigSyst=true, bool plotGlobalSysts=false);
void plotSystsAll(const char* apoiname="NJpsi", bool plotEffSyst=false, bool plotSigSyst=true, bool plotGlobalSysts=false, bool is16004=false);


void plotSystsAll(const char* apoiname, bool plotEffSyst, bool plotSigSyst, bool plotGlobalSysts, bool is16004) {
  spoiname = apoiname;
  if (!spoiname.CompareTo("NJpsi") && !spoiname.CompareTo("BJpsi") && !spoiname.CompareTo("NJpsi_prompt") && !spoiname.CompareTo("NJpsi_nonprompt"))
  {
    cout << "[ERROR] : unknown systematic" << endl;
    return;
  }
  
  if (!is16004)
  {
    nameTag = "16025";
    //pt dependence in rapidity bins
    plotSysts(anabin(0,0.6,6.5,50,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0.6,1.2,6.5,50,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.2,1.8,6.5,50,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.8,2.4,6.5,50,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    plotSysts(anabin(0,0.6,6.5,50,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0.6,1.2,6.5,50,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.2,1.8,6.5,50,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.8,2.4,6.5,50,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    //pt dependence in centrality bins
    plotSysts(anabin(0,2.4,6.5,50,0,20),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,20,60),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,60,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    //centrality dependence in rapidity bins
    plotSysts(anabin(0,0.6,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0.6,1.2,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.2,1.8,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.8,2.4,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    //centrality dependence at fwd rapidity and low pt
    plotSysts(anabin(1.8,2.4,3,6.5,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    //rapidity dependence
    plotSysts(anabin(0,2.4,6.5,50,0,200),"rap","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(0,2.4,6.5,50,0,200),"rap","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
  }
  else
  {
    bins16004 = true;
    nameTag = "16004";
    
    //pt dependence in rapidity bins
    plotSysts(anabin(0,1.6,6.5,30,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.6,2.4,6.5,30,0,200),"pt","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    plotSysts(anabin(0,1.6,6.5,30,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.6,2.4,6.5,30,0,200),"pt","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
    
    //centrality dependence in rapidity bins
    plotSysts(anabin(0,1.6,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
    plotSysts(anabin(1.6,2.4,6.5,50,0,200),"cent","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
  }
  
}

void plotSysts(anabin thebin, string xaxis, string collTag, bool plotEffSyst, bool plotSigSyst, bool plotGlobalSysts) {
  float rapmin = thebin.rapbin().low();
  float rapmax = thebin.rapbin().high();
  float ptmin = thebin.ptbin().low();
  float ptmax = thebin.ptbin().high();
  int centmin = thebin.centbin().low();
  int centmax = thebin.centbin().high();
  
  vector<string> tags;
  vector<TGraphErrors*> graphs;
  vector<map<anabin, syst> > systs;
  
  // total
  systs.push_back(readSyst_all(collTag.c_str(),spoiname.Data(),nameTag.Data(),plotEffSyst,plotSigSyst,"../"));
  tags.push_back(systs.back().begin()->second.name);
  
  // detail
  if (collTag=="PbPb") {
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_trg_binned.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_muid_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_sta_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_trg_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_muid_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_sta_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_trg_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_ptWeighting.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_acc.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PbPb_eff.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
  } else {
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_trg_binned.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_muid_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_sta_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_trg_syst.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_muid_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_sta_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_trg_stat.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_ptWeighting.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_acc.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst(Form("csv/stat_%s_%s_PP_eff.csv",nameTag.Data(),spoiname.Data())));
    tags.push_back(systs.back().begin()->second.name);
  }
  
  // global systs
  if (plotGlobalSysts) {
    systs.push_back(readSyst("csv/syst_PbPb_taa.csv"));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst("csv/syst_PbPb_Nmb.csv"));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst("csv/syst_PP_lumi.csv"));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst("csv/syst_PP_trk.csv"));
    tags.push_back(systs.back().begin()->second.name);
    systs.push_back(readSyst("csv/syst_PbPb_trk.csv"));
    tags.push_back(systs.back().begin()->second.name);
  }
  
  set<anabin> sb;
  if (bins16004) sb = allbins16004();
  else sb = allbins();
  
  for (unsigned int i=0; i<systs.size(); i++) {
    map<anabin, syst> thesyst = systs[i];
    vector<double> x, y, dx, dy;
    
    double valmax=-1e99;
    
    for (set<anabin>::const_iterator it=sb.begin(); it!=sb.end(); it++) {
      anabin it2 = *it;
      
      if (!binok(thebin,xaxis,it2,false)) continue;
      if (thesyst.find(*it)==thesyst.end()) continue;
      
      double low, high;
      if (xaxis=="pt") {
        low= it->ptbin().low();
        high = it->ptbin().high();
      } else if (xaxis=="rap") {
        low= it->rapbin().low();
        high = it->rapbin().high();
        if ((low==0 && high<=0.61 && high>=0.59 ) || (low>=0.59 && low<=0.61 && high>=1.19 && high <=1.21) || (low>=1.19 && low<=1.21 && high>=1.79 && high <=1.81) || (low>=1.79 && low<=1.81 && high>=2.39 && high <=2.41)) continue;
      } else {
        low= it->centbin().low()/2;
        high = it->centbin().high()/2;
      }
      x.push_back((low+high)/2.);
      dx.push_back((high-low)/2.);
      y.push_back(0);
      dy.push_back(thesyst[*it].value);
      valmax = max(valmax,dy.back());
    }
    
    TGraphErrors *thegraph = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
    TH1F *haxes=NULL;
    if (xaxis=="pt") {
      haxes = new TH1F(tags[i].c_str(),Form(";p_{T} (GeV/c);Syst_%s",spoiname.Data()),1,0,50);
    } else if (xaxis=="cent") {
      haxes = new TH1F(tags[i].c_str(),Form(";Centrality percentile;Syst_%s",spoiname.Data()),1,0,100);
    } else { // if (xaxis=="rap")
      haxes = new TH1F(tags[i].c_str(),Form(";|y|;Syst_%s",spoiname.Data()),1,0,2.4);
    }
    haxes->GetYaxis()->SetLimits(-1.3*valmax, 2.*valmax);
    haxes->GetYaxis()->SetRangeUser(-1.3*valmax, 2.*valmax);
    haxes->GetYaxis()->SetTitleOffset(1.4);
    thegraph->SetHistogram(haxes);

    thegraph->Sort();
    graphs.push_back(thegraph);
  }
  
  
  plotGraphs(graphs, tags, "systematics", collTag,
             Form("pt%i%i_rap%i%i_cent%i%i_%s",(int)ptmin*10,(int)ptmax*10,(int)rapmin*10,(int)rapmax*10,centmin,centmax,nameTag.Data()),
             Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f, %i-%i%s",rapmin,rapmax,ptmin,ptmax,centmin/2,centmax/2,"%"));
}
