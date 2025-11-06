#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"

#include "../Fitter/Macros/CMS/CMS_lumi.C"
#include "../Fitter/Macros/CMS/tdrstyle.C"
#include "../Fitter/Macros/Utilities/texUtils.h"
#include "JpsiDrawing.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


TH1F* integrateHist(TH1F *hist) {
   TString name = hist->GetName(); name = name + "_int";
   TString title = hist->GetTitle(); title = title + " integrated";
   double integral, integralerror;
   int nbins = hist->GetNbinsX();
   integral = hist->IntegralAndError(1,nbins,integralerror);
   TH1F *ans = new TH1F(name, title, 1, hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinUpEdge(nbins));
   ans->SetBinContent(1,integral);
   ans->SetBinError(1,integralerror);
   return ans;
}


class drawingEff {
public:

  // Bin boundaries for efficiency histograms
  vector<float> bins_4rap;
  vector<float> bins_3cent;
  int nbins_4rap;
  int nbins_3cent;
  bool ispbpb, isacc;

  // Input file
  TFile *finput;

  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  vector<TH1F*> hnum_cent_rap;
  vector<TH1F*> hden_cent_rap;
  vector<TGraphAsymmErrors*> heff_cent_rap;
  // Eff vs pT in 4+1 |y| regions
  vector<TH1F*> hnum_pt_rap;
  vector<TH1F*> hden_pt_rap;
  vector<TGraphAsymmErrors*> heff_pt_rap;
  // Eff vs pT in 3 centrality regions
  vector<TH1F*> hnum_pt_cent;
  vector<TH1F*> hden_pt_cent;
  vector<TGraphAsymmErrors*> heff_pt_cent;
  // Eff vs rap integrated
  TH1F *hnum_rap;
  TH1F *hden_rap;
  TGraphAsymmErrors *heff_rap;

  drawingEff(string fname, bool ispbpb, bool isacc);
  ~drawingEff();
  void loadHisto();
  void getEfficiency();
  void drawHisto(vector<string> *outname);
  void fixCentPp(TH1F*);
  void checkUnderFlow(TH1 *hnum, TH1 *hden);
};

void drawingEff::fixCentPp(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   float y = hist->GetBinContent(1);
   float dy = hist->GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist->SetBinContent(i,y);
      hist->SetBinError(i,dy);
   }
}

drawingEff::~drawingEff(){
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  for (size_t i=0; i<hnum_cent_rap.size(); i++) {
    delete hnum_cent_rap[i];
    delete hden_cent_rap[i];
    delete heff_cent_rap[i];
  }
  // Eff vs pT in 4+1 |y| regions
  for (size_t i=0; i<hnum_pt_rap.size(); i++) {
    delete hnum_pt_rap[i];
    delete hden_pt_rap[i];
    delete heff_pt_rap[i];
  }
  // Eff vs pT in 3 centrality regions
  for (size_t i=0; i<hnum_pt_cent.size(); i++) {
    delete hnum_pt_cent[i];
    delete hden_pt_cent[i];
    delete heff_pt_cent[i];
  }
  // Eff vs rap integrated
  delete hnum_rap;
  delete hden_rap;
  delete heff_rap;

  delete finput;
}

drawingEff::drawingEff(string fname, bool ispbpb_, bool isacc_){
  finput = new TFile(fname.c_str());
  ispbpb = ispbpb_;
  isacc = isacc_;

  float bins_4rap_[] = {0, 0.6, 1.2, 1.8, 2.4};
  float bins_3cent_[] = {0, 10, 30, 100};
  
  nbins_4rap = sizeof(bins_4rap_)/sizeof(float) -1;
  for (int i=0; i<=nbins_4rap; i++)
    bins_4rap.push_back(bins_4rap_[i]);

  nbins_3cent = sizeof(bins_3cent_)/sizeof(float) -1;
  for (int i=0; i<=nbins_3cent; i++) 
    bins_3cent.push_back(bins_3cent_[i]);
}

void drawingEff::checkUnderFlow(TH1 *hnum, TH1 *hden){
  for (int j=0; j<=hnum->GetNbinsX()+1; j++) {
    double num0 = hnum->GetBinContent(j);
    double den0 = hden->GetBinContent(j);
    
    // If underflow bin has more entries in numerator than denominator, set it to 0
    if ((j==0 && num0>den0) || (j==hnum->GetNbinsX()+1 && num0>den0)) {
      hnum->SetBinContent(j,0);
      hden->SetBinContent(j,0);
      hnum->SetBinError(j,0);
      hden->SetBinError(j,0);
    } else if (num0>den0) {
      hnum->SetBinContent(j,0);
      hden->SetBinContent(j,0);
      hnum->SetBinError(j,0);
      hden->SetBinError(j,0);
    }
    
  }
}

void drawingEff::loadHisto(){
  if (!finput) cout << "Cannot open input file!\n" << endl;
  
  // Load denominator, numerators
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  TH1F *hnum = (TH1F*)finput->Get("hnum_cent_rap0024");
  TH1F *hden = (TH1F*)finput->Get("hden_cent_rap0024");
  if (!ispbpb) { // set all bin contents to same as the 1st bin
    fixCentPp(hnum);
    fixCentPp(hden);
  }
  checkUnderFlow(hnum,hden);
  hnum_cent_rap.push_back(hnum);
  hden_cent_rap.push_back(hden);
  for (int i=0; i<nbins_4rap; i++) {
    TH1F *hnum1 = (TH1F*)finput->Get( Form("hnum_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TH1F *hden1 = (TH1F*)finput->Get( Form("hden_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    if (!ispbpb) { // set all bin contents to same as the 1st bin
      fixCentPp(hnum1);
      fixCentPp(hden1);
    }
    checkUnderFlow(hnum1,hden1);
    hnum_cent_rap.push_back(hnum1);
    hden_cent_rap.push_back(hden1);
  }
  hnum = (TH1F*)finput->Get("hnum_cent_rap1824_pt3065");
  hden = (TH1F*)finput->Get("hden_cent_rap1824_pt3065");
  if (!ispbpb) { // set all bin contents to same as the 1st bin
    fixCentPp(hnum);
    fixCentPp(hden);
  }
  checkUnderFlow(hnum,hden);
  hnum_cent_rap.push_back(hnum);
  hden_cent_rap.push_back(hden);

  // Eff vs pT in 4+1 |y| regions
  hnum = (TH1F*)finput->Get("hnum_pt_rap0024");
  hden = (TH1F*)finput->Get("hden_pt_rap0024");
  checkUnderFlow(hnum,hden);
  hnum_pt_rap.push_back(hnum);
  hden_pt_rap.push_back(hden);
  for (int i=0; i<nbins_4rap; i++) {
    TH1F *hnum1 = (TH1F*)finput->Get( Form("hnum_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TH1F *hden1 = (TH1F*)finput->Get( Form("hden_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    checkUnderFlow(hnum1,hden1);
    hnum_pt_rap.push_back(hnum1);
    hden_pt_rap.push_back(hden1);
  }

  // Eff vs pT in 3 centrality regions
  for (int i=0; i<nbins_3cent-1; i++) { // for pp, only the 1st bin is necessary
    TH1F *hnum1 = (TH1F*)finput->Get( Form("hnum_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    TH1F *hden1 = (TH1F*)finput->Get( Form("hden_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    checkUnderFlow(hnum1,hden1);
    hnum_pt_cent.push_back(hnum1);
    hden_pt_cent.push_back(hden1);
  }
  hnum = (TH1F*)finput->Get( Form("hnum_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
  hden = (TH1F*)finput->Get( Form("hden_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
  checkUnderFlow(hnum,hden);
  hnum_pt_cent.push_back(hnum);
  hden_pt_cent.push_back(hden);

  // Eff vs rap integrated
  hnum_rap = (TH1F*)finput->Get("hnum_rap");
  hden_rap = (TH1F*)finput->Get("hden_rap");
  checkUnderFlow(hnum_rap,hden_rap);
}

void drawingEff::getEfficiency(){

  setTDRStyle();
  gStyle->SetEndErrorSize(3);

  for (int i=0; i<nbins_4rap+2; i++) {
    TGraphAsymmErrors *geff = new TGraphAsymmErrors(hnum_cent_rap[i],hden_cent_rap[i],"n");
    const int np = geff->GetN();
    cout << hnum_cent_rap[i]->GetName() << endl;
    for (int a=0; a<np; a++) {
      double gey = geff->GetErrorYhigh(a);
      double gex = geff->GetErrorXhigh(a);
      if (gey!=gey) { // when error is nan, error will be calculated by hand!
        double eyd = hden_cent_rap[i]->GetBinError(a+1);
        double eyn = hnum_cent_rap[i]->GetBinError(a+1);
        double yd = hden_cent_rap[i]->GetBinContent(a+1);
        double yn = hnum_cent_rap[i]->GetBinContent(a+1);
        double uncert = yn/yd * sqrt( pow(eyd/yd,2) + pow(eyn/yn,2) ) / 2.0;
        cout << "   getEfficiency("<<a<<")\t" << yn/yd << "  " << uncert << endl;
        geff->SetPointError(a,gex,gex,uncert,uncert);
      }
    }
    heff_cent_rap.push_back(geff);

    string gname = hnum_cent_rap[i]->GetName();
    gname.replace(gname.begin(),gname.begin()+4,"heff");
    heff_cent_rap[i]->SetName(gname.c_str());
    
    heff_cent_rap[i]->GetXaxis()->SetTitle(hnum_cent_rap[i]->GetXaxis()->GetTitle());
    heff_cent_rap[i]->GetYaxis()->SetTitle(hnum_cent_rap[i]->GetYaxis()->GetTitle());
    heff_cent_rap[i]->GetXaxis()->SetRangeUser(0,100);
    SetHistStyle(heff_cent_rap[i],i,i,0,1.3);
  }
  
  for (int i=0; i<nbins_4rap+1; i++) {
    TGraphAsymmErrors *geff = new TGraphAsymmErrors(hnum_pt_rap[i],hden_pt_rap[i],"n");
    const int np = geff->GetN();
    cout << hnum_pt_rap[i]->GetName() << endl;
    for (int a=0; a<np; a++) {
      double gey = geff->GetErrorYhigh(a);
      double gex = geff->GetErrorXhigh(a);
      if (gey!=gey) { // when error is nan, error will be calculated by hand!
        double eyd = hden_pt_rap[i]->GetBinError(a+1);
        double eyn = hnum_pt_rap[i]->GetBinError(a+1);
        double yd = hden_pt_rap[i]->GetBinContent(a+1);
        double yn = hnum_pt_rap[i]->GetBinContent(a+1);
        double uncert = yn/yd * sqrt( pow(eyd/yd,2) + pow(eyn/yn,2) ) / 2.0;
        cout << "   getEfficiency("<<a<<")\t" << yn/yd << "  " << uncert << endl;
        geff->SetPointError(a,gex,gex,uncert,uncert);
      }
    }
    heff_pt_rap.push_back(geff);
   
    string gname = hnum_pt_rap[i]->GetName();
    gname.replace(gname.begin(),gname.begin()+4,"heff");
    heff_pt_rap[i]->SetName(gname.c_str());

    heff_pt_rap[i]->GetXaxis()->SetTitle(hnum_pt_rap[i]->GetXaxis()->GetTitle());
    heff_pt_rap[i]->GetYaxis()->SetTitle(hnum_pt_rap[i]->GetYaxis()->GetTitle());
    heff_pt_rap[i]->GetXaxis()->SetRangeUser(0,50);
    SetHistStyle(heff_pt_rap[i],i,i,0,1.3);
  }
 
  for (int i=0; i<nbins_3cent; i++) {
    if (!ispbpb && i>0) {
      TGraphAsymmErrors *geff = (TGraphAsymmErrors *)heff_pt_cent[0]->Clone();
      string gname = hnum_pt_cent[i]->GetName();
      cout << hnum_pt_cent[i]->GetName() << endl;
      gname.replace(gname.begin(),gname.begin()+4,"heff");
      geff->SetName(gname.c_str());
      heff_pt_cent.push_back(geff);
    } else {
      TGraphAsymmErrors *geff = new TGraphAsymmErrors(hnum_pt_cent[i],hden_pt_cent[i],"n");
      const int np = geff->GetN();
      cout << hnum_pt_cent[i]->GetName() << endl;
      for (int a=0; a<np; a++) {
        double gey = geff->GetErrorYhigh(a);
        double gex = geff->GetErrorXhigh(a);
        if (gey!=gey) { // when error is nan, error will be calculated by hand!
          double eyd = hden_pt_cent[i]->GetBinError(a+1);
          double eyn = hnum_pt_cent[i]->GetBinError(a+1);
          double yd = hden_pt_cent[i]->GetBinContent(a+1);
          double yn = hnum_pt_cent[i]->GetBinContent(a+1);
          double uncert = yn/yd * sqrt( pow(eyd/yd,2) + pow(eyn/yn,2) ) / 2.0;
          cout << "   getEfficiency("<<a<<")\t" << yn/yd << "  " << uncert << endl;
          geff->SetPointError(a,gex,gex,uncert,uncert);
        }
      }
      heff_pt_cent.push_back(geff);
      
      string gname = hnum_pt_cent[i]->GetName();
      gname.replace(gname.begin(),gname.begin()+4,"heff");
      heff_pt_cent[i]->SetName(gname.c_str());

      heff_pt_cent[i]->GetXaxis()->SetTitle(hnum_pt_cent[i]->GetXaxis()->GetTitle());
      heff_pt_cent[i]->GetYaxis()->SetTitle(hnum_pt_cent[i]->GetYaxis()->GetTitle());
      heff_pt_cent[i]->GetXaxis()->SetRangeUser(0,50);
      SetHistStyle(heff_pt_cent[i],i,i,0,1.3);
    }
  }

  heff_rap = new TGraphAsymmErrors(hnum_rap,hden_rap,"n");
  const int np = heff_rap->GetN();
  cout << hnum_rap->GetName() << endl;
  for (int a=0; a<np; a++) {
    double gey = heff_rap->GetErrorYhigh(a);
    double gex = heff_rap->GetErrorXhigh(a);
    if (gey!=gey) { // when error is nan, error will be calculated by hand!
      double eyd = hden_rap->GetBinError(a+1);
      double eyn = hnum_rap->GetBinError(a+1);
      double yd = hden_rap->GetBinContent(a+1);
      double yn = hnum_rap->GetBinContent(a+1);
      double uncert = yn/yd * sqrt( pow(eyd/yd,2) + pow(eyn/yn,2) ) / 2.0;
      cout << "   getEfficiency("<<a<<")\t" << yn/yd << "  " << uncert << endl;
      heff_rap->SetPointError(a,gex,gex,uncert,uncert);
    }
  }
  string gname = hnum_rap->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  heff_rap->SetName(gname.c_str());

  heff_rap->GetXaxis()->SetTitle(hnum_rap->GetXaxis()->GetTitle());
  heff_rap->GetYaxis()->SetTitle(hnum_rap->GetYaxis()->GetTitle());
  heff_rap->GetXaxis()->SetRangeUser(0,2.4);
  SetHistStyle(heff_rap,0,0,0,1.3);
}

void drawingEff::drawHisto(vector<string> *outname){

  gStyle->SetEndErrorSize(3);

  TCanvas *can = new TCanvas("can","can",600,600);
  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.035);
  TLegend *leg = new TLegend(0.55,0.60,0.94,0.88);
  SetLegendStyle(leg);

  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  if (isacc) heff_cent_rap[0]->GetYaxis()->SetTitle("Acceptance");
  heff_cent_rap[0]->Draw("ap");
  leg->AddEntry(heff_cent_rap[0],"|y|: 0-2.4, 6.5-50 GeV/c","p");
  for (int i=0; i<nbins_4rap; i++) {
    if (isacc) heff_cent_rap[i+1]->GetYaxis()->SetTitle("Acceptance");
    heff_cent_rap[i+1]->Draw("p");
    leg->AddEntry(heff_cent_rap[i+1],Form("|y|: %.1f-%.1f, 6.5-50 GeV/c",bins_4rap[i],bins_4rap[i+1]),"p");
  }
  if (isacc) heff_cent_rap[nbins_4rap+1]->GetYaxis()->SetTitle("Acceptance");
  heff_cent_rap[nbins_4rap+1]->Draw("p");
  leg->AddEntry(heff_cent_rap[nbins_4rap+1],Form("|y|: %.1f-%.1f, 3-6.5 GeV/c",bins_4rap[nbins_4rap-1],bins_4rap[nbins_4rap]),"p");
  leg->Draw();
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  can->SaveAs(Form("%s_cent_rap.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_cent_rap.pdf",outname->at(1).c_str()));
  can->SaveAs(Form("%s_cent_rap.root",outname->at(1).c_str()));

  delete can;
  delete leg;
  

  // Eff vs pT in 4+1 |y| regions
  can = new TCanvas("can","can",600,600);
  leg = new TLegend(0.6,0.68,0.9,0.88);
  SetLegendStyle(leg);

  if (isacc) heff_pt_rap[0]->GetYaxis()->SetTitle("Acceptance");
  heff_pt_rap[0]->Draw("ap");
  if (ispbpb) leg->AddEntry(heff_pt_rap[0],"|y|: 0-2.4, 0-100%","p");
  else leg->AddEntry(heff_pt_rap[0],"|y|: 0-2.4","p");
  for (int i=0; i<nbins_4rap; i++) {
    if (isacc) heff_pt_rap[i+1]->GetYaxis()->SetTitle("Acceptance");
    heff_pt_rap[i+1]->Draw("p");
    if (ispbpb) leg->AddEntry(heff_pt_rap[i+1],Form("|y|: %.1f-%.1f, 0-100%%",bins_4rap[i],bins_4rap[i+1]),"p");
    else leg->AddEntry(heff_pt_rap[i+1],Form("|y|: %.1f-%.1f",bins_4rap[i],bins_4rap[i+1]),"p");
  }
  leg->Draw();
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  can->SaveAs(Form("%s_pt_rap.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_pt_rap.pdf",outname->at(1).c_str()));
  can->SaveAs(Form("%s_pt_rap.root",outname->at(1).c_str()));

  delete can;
  delete leg;
  
 
  // Eff vs pT in 3 centrality regions
  can = new TCanvas("can","can",600,600);
  leg = new TLegend(0.67,0.68,0.9,0.88);
  SetLegendStyle(leg);

  if (isacc) heff_pt_cent[0]->GetYaxis()->SetTitle("Acceptance");
  heff_pt_cent[0]->Draw("ap");
  leg->AddEntry(heff_pt_cent[0],Form("%.0f-%.0f%%",bins_3cent[0],bins_3cent[1]),"p");
  if (ispbpb) {
    for (int i=1; i<nbins_3cent; i++) {
      if (isacc) heff_pt_cent[i]->GetYaxis()->SetTitle("Acceptance");
      heff_pt_cent[i]->Draw("p");
      leg->AddEntry(heff_pt_cent[i],Form("%.0f-%.0f%%",bins_3cent[i],bins_3cent[i+1]),"p");
    }
    leg->Draw();
  }
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  lat->DrawLatex(0.2,0.80,"|y| < 2.4");
  can->SaveAs(Form("%s_pt_cent.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_pt_cent.pdf",outname->at(1).c_str()));
  can->SaveAs(Form("%s_pt_cent.root",outname->at(1).c_str()));

  delete can;
  delete leg;
  

  // Eff vs rap integrated
  can = new TCanvas("can","can",600,600);

  if (isacc) heff_rap->GetYaxis()->SetTitle("Acceptance");
  heff_rap->Draw("ap");
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  lat->DrawLatex(0.2,0.80,"6.5-50 GeV/c");
  if (ispbpb) lat->DrawLatex(0.2,0.75,"0-100%");
  can->SaveAs(Form("%s_rap.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_rap.pdf",outname->at(1).c_str()));
  can->SaveAs(Form("%s_rap.root",outname->at(1).c_str()));

  delete can;
  delete lat;
  
}










////////////////////////////////////////////////////////////////
////// Draw eff plots in different kinematic regions
////////////////////////////////////////////////////////////////
void plotEffs(string subdir="eff/nominal") {
  string dir = "figs";
  bool ispbpb=false, isacc=false;

  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  vector<string> latex;
  latex.push_back("prompt J/#psi (pp)");
  latex.push_back(Form("%s/%s/jpsi_pp",dir.c_str(),subdir.c_str()));
  drawingEff fjpsi_pp(Form("files/%s/histos_jpsi_pp.root",subdir.c_str()),ispbpb,isacc);
  fjpsi_pp.loadHisto();
  fjpsi_pp.getEfficiency();
  fjpsi_pp.drawHisto(&latex);

  latex.clear();
  latex.push_back("nonprompt J/#psi (pp)");
  latex.push_back(Form("%s/%s/npjpsi_pp",dir.c_str(),subdir.c_str()));
  drawingEff fnpjpsi_pp(Form("files/%s/histos_npjpsi_pp.root",subdir.c_str()),ispbpb,isacc);
  fnpjpsi_pp.loadHisto();
  fnpjpsi_pp.getEfficiency();
  fnpjpsi_pp.drawHisto(&latex);

  ispbpb=true; 

  latex.clear();
  latex.push_back("prompt J/#psi (PbPb)");
  latex.push_back(Form("%s/%s/jpsi_pbpb",dir.c_str(),subdir.c_str()));
  drawingEff fjpsi_pbpb(Form("files/%s/histos_jpsi_pbpb.root",subdir.c_str()),ispbpb&&!isacc,isacc);
  fjpsi_pbpb.loadHisto();
  fjpsi_pbpb.getEfficiency();
  fjpsi_pbpb.drawHisto(&latex);

  latex.clear();
  latex.push_back("nonprompt J/#psi (PbPb)");
  latex.push_back(Form("%s/%s/npjpsi_pbpb",dir.c_str(),subdir.c_str()));
  drawingEff fnpjpsi_pbpb(Form("files/%s/histos_npjpsi_pbpb.root",subdir.c_str()),ispbpb&&!isacc,isacc);
  fnpjpsi_pbpb.loadHisto();
  fnpjpsi_pbpb.getEfficiency();
  fnpjpsi_pbpb.drawHisto(&latex);
}









//////////////////////////////////////////////////////////////////
// For 1 kinematic region, draw different SFs on a same canvas
// -> Call "drawMultiples()" function to execute this macro
//////////////////////////////////////////////////////////////////
void plotMultipleSamples(vector<TGraphAsymmErrors*> heff, vector<string> *histname, string *outname, int xVar, const double *xaxis, string rap, string pt, string cent, bool ispbpb, bool isacc=false) {

  setTDRStyle();
  gStyle->SetEndErrorSize(3);

  const int NGraph = heff.size();
  // Draw ratios
  vector<TGraphAsymmErrors*> ratio;
  vector<TH1F*> geff;
  for (int i=0; i<NGraph; i++) {
    int Ntotal = heff[i]->GetN();
    
    // create ratio histograms
    double *y = heff[i]->GetY();
    double *eylow = heff[i]->GetEYlow();
    double *eyhigh = heff[i]->GetEYhigh();
    
    TH1F *htmp = new TH1F(Form("h%d_%s",i,heff[i]->GetName()),"htmp",Ntotal,xaxis);
    cout << htmp->GetName() << endl;
    for (int a=0; a<Ntotal; a++) {
      htmp->SetBinContent(a+1,y[a]);
      double ey = (eylow[a]>eyhigh[a] ? eylow[a] : eyhigh[a]);
      htmp->SetBinError(a+1,ey);
      cout << "\t" << a << "  " << xaxis[a] << " " << xaxis[a+1] << "  " << y[a] << "  " << ey << endl;
    }
    geff.push_back(htmp);
  }
    
  // except the 0th item, other histograms will be divided by 0th item
  double max = 0; // y-axis maximum for ratio plots
  for (int i=1; i<NGraph; i++) {
    TGraphAsymmErrors *rtmp = new TGraphAsymmErrors(geff[i],geff[0],"pois");
    ratio.push_back(rtmp);
    max = TMath::MaxElement(rtmp->GetN(), rtmp->GetY());
  }
  for (int i=1; i<NGraph; i++) {
    SetHistStyle(ratio[i-1],i,i,0.90,max+0.17);
    if (rap=="1.8-2.4" && pt=="3-6.5") {
      SetHistStyle(ratio[i-1],i,i,0.90,max+0.17);
    }
    if (cent=="") {
      ratio[i-1]->GetXaxis()->SetRangeUser(0,100);
    } else if (pt=="" && cent=="0-100") {
      ratio[i-1]->GetXaxis()->SetRangeUser(0,50);
    } else if (pt=="" && rap=="0-2.4") {
      ratio[i-1]->GetXaxis()->SetRangeUser(0,50);
    } else if (rap=="" && pt=="6.5-30" && cent=="0-100") {
      ratio[i-1]->GetXaxis()->SetRangeUser(0,2.4);
    }
    ratio[i-1]->GetXaxis()->SetTitleOffset(1.00);
    ratio[i-1]->GetYaxis()->SetTitleSize(0.10);
    ratio[i-1]->GetYaxis()->SetTitleOffset(0.60);
    ratio[i-1]->GetXaxis()->SetTitleSize(0.10);
    ratio[i-1]->GetXaxis()->SetLabelSize(0.10);
    ratio[i-1]->GetYaxis()->SetLabelSize(0.10);
    string title = heff[0]->GetXaxis()->GetTitle();
    ratio[i-1]->GetXaxis()->SetTitle(title.c_str());
    // setting for pp/pbpb ratio
    ratio[i-1]->GetYaxis()->SetTitle("Ratio(PbPb/pp)");
    if (isacc) ratio[i-1]->GetYaxis()->SetRangeUser(0.9,1.5);
    else ratio[i-1]->GetYaxis()->SetRangeUser(0.4,1.1);
    // setting for tnp scale factors ratio
//    ratio[i-1]->GetYaxis()->SetTitle("Ratio(SF/no SF)");
//    ratio[i-1]->GetYaxis()->SetRangeUser(0.9,1.5);
  }

  TCanvas *can = new TCanvas("can","can",600,700);
  can->Draw();
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.35,1.0,1);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.35);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.30);
  pad2->Draw();
  
  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.05);
  TLegend *leg = new TLegend(0.60,0.65,0.94,0.88);
  SetLegendStyle(leg);

  for (int i=0; i<NGraph; i++) {
    pad1->cd();
    SetHistStyle(heff[i],i,i,0,1.3);
    leg->AddEntry(heff[i],histname->at(i).c_str(),"p");
    heff[i]->GetYaxis()->SetTitle(isacc ? "Acceptance" : "Efficiency");
    
    if (i==0) heff[i]->Draw("ap");
    else heff[i]->Draw("p");
  }
  leg->Draw();
  can->Update();


  pad2->cd();
  TLine *line = new TLine(); line->SetLineColor(kGray+2);
  for (int i=1; i<NGraph; i++) {
    if (cent=="" && ispbpb) {
      line->DrawLine(0,1,100,1);
    } else if (pt=="" && cent=="0-100") {
      line->DrawLine(0,1,50,1);
    } else if (pt=="" && rap=="0-2.4") {
      line->DrawLine(0,1,50,1);
    } else if (rap=="" && pt=="6.5-30" && cent=="0-100") {
      line->DrawLine(0,1,2.4,1);
    }
    can->Update();
    if (i==1) ratio[i-1]->Draw("ap");
    else ratio[i-1]->Draw("p");
    can->Update();
  }

  gSystem->mkdir(Form("%s",outname->c_str()),kTRUE);

  pad1->cd();
  if (xVar == 0) { // vs. Cent
    lat->DrawLatex(0.2,0.85,Form("|y|: %s, %s GeV/c",rap.c_str(),pt.c_str()));
    can->Update();
    
    can->SaveAs(Form("%s/cent_rap%s_pt%s.png",outname->c_str(),rap.c_str(),pt.c_str()));
    can->SaveAs(Form("%s/cent_rap%s_pt%s.pdf",outname->c_str(),rap.c_str(),pt.c_str()));

  } else if (xVar == 1) { // vs. pt in 4+1 |y| regions
    if (ispbpb) lat->DrawLatex(0.2,0.85,Form("|y|: %s, 0-100%%",rap.c_str()));
    else lat->DrawLatex(0.2,0.85,Form("|y|: %s",rap.c_str()));
    can->Update();
    
    can->SaveAs(Form("%s/pt_rap%s.png",outname->c_str(),rap.c_str()));
    can->SaveAs(Form("%s/pt_rap%s.pdf",outname->c_str(),rap.c_str()));
  
  } else if (xVar == 2) { // vs. pt in 3 cent regions
    if (ispbpb) lat->DrawLatex(0.2,0.85,Form("%s%%",cent.c_str()));
    can->Update();

    can->SaveAs(Form("%s/pt_cent%s.png",outname->c_str(),cent.c_str()));
    can->SaveAs(Form("%s/pt_cent%s.pdf",outname->c_str(),cent.c_str()));

  } else if (xVar == 3) { // vs. rap integrated over all regions
    lat->DrawLatex(0.2,0.80,"6.5-50 GeV/c");
    if (ispbpb) lat->DrawLatex(0.2,0.75,"0-100%");
    can->Update();

    can->SaveAs(Form("%s/rap.png",outname->c_str()));
    can->SaveAs(Form("%s/rap.pdf",outname->c_str()));
  }
   
  for (int i=0; i<NGraph; i++) {
    delete geff[i];
    if (i<NGraph-1) { // except the 0th item, other histograms will be divided by 0th item
      delete ratio[i];
    }
  }
 
  delete leg;
  delete lat;
  delete pad1;
  delete pad2;
  delete can;
}

void callMultiples(vector<string> *inputfile, bool ispbpb, string *outname, vector<string> *histname) {
  
  const unsigned int isize = inputfile->size();
  vector<drawingEff *> fsource;
  for (unsigned int i=0; i<isize; i++) {
    cout << "creating drawingEff: " << inputfile->at(i) << endl;
    fsource.push_back(new drawingEff(inputfile->at(i).c_str(),ispbpb,false));
    fsource[i]->loadHisto();
    fsource[i]->getEfficiency();
  }

  string raparr[] = {"0-2.4","0-0.6","0.6-1.2","1.2-1.8","1.8-2.4","1.8-2.4"};
  string ptarr[] = {"6.5-30","3-6.5"};
  string centarr[] = {"0-10","10-30","30-100"};

  vector<TGraphAsymmErrors*> heff_cent_rap;
  vector<TGraphAsymmErrors*> heff_pt_rap;
  vector<TGraphAsymmErrors*> heff_pt_cent;
  vector<TGraphAsymmErrors*> heff_rap;
  
  //xVar is for drawing latex on plots
  //xVar == 0 : vs. Cent
  //xVar == 1 : vs. pt in 4+1 |y| regions
  //xVar == 2 : vs. pt in 3 cent regions
  //xVar == 3 : vs. rap integrated over all regions
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  unsigned int heff_size = fsource[0]->heff_cent_rap.size();
  cout << "heff_size: " << heff_size << endl;
  for (unsigned int i=0; i<heff_size; i++) { // heff_size == number of loops for each 4+1 |y| region
    int xVar = 0;
    for (unsigned int fidx=0; fidx<isize; fidx++) {
      heff_cent_rap.push_back(fsource[fidx]->heff_cent_rap[i]); // Compare fdix input files
    }
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_cent_rap[i];
    const int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = raparr[i];
    string pt = "6.5-30";
    string cent = "";
    if (i+1 == heff_size) {
      rap = "1.8-2.4";
      pt = "3-6.5";
    }
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_cent_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb);
    
    heff_cent_rap.clear();
    delete[] xaxis;
  }
  // Eff vs pT in 4+1 |y| regions
  heff_size = fsource[0]->heff_pt_rap.size();
  for (unsigned int i=0; i<heff_size; i++) {
    int xVar = 1;
    for (unsigned int fidx=0; fidx<isize; fidx++) {
      heff_pt_rap.push_back(fsource[fidx]->heff_pt_rap[i]);
    }
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_pt_rap[i];
    int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = raparr[i];
    string pt = "";
    string cent = "0-100";
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_pt_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb);

    heff_pt_rap.clear();
    delete[] xaxis;
  }
  // Eff vs pT in 3 centrality regions
  heff_size = fsource[0]->heff_pt_cent.size();
  for (unsigned int i=0; i<heff_size; i++) {
    int xVar = 2;
    for (unsigned int fidx=0; fidx<isize; fidx++) {
      heff_pt_cent.push_back(fsource[fidx]->heff_pt_cent[i]);
    }
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_pt_cent[i];
    int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = "0-2.4";
    string pt = "";
    string cent = centarr[i];
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_pt_cent, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb);

    heff_pt_cent.clear();
    delete[] xaxis;
  }
  // Eff vs rap integrated
  int xVar = 3;
  for (unsigned int fidx=0; fidx<isize; fidx++) {
    heff_rap.push_back(fsource[fidx]->heff_rap);
  }
  // Defining x-axis to re-making ratio histogram
  TH1F *hden = fsource[0]->hden_rap;
  int nbins = hden->GetNbinsX();
  double *xaxis = new double[nbins+1];
  for (int a=0; a<nbins; a++) {
    xaxis[a] = hden->GetBinLowEdge(a+1);
  }
  xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
  for (int a=0; a<nbins+1; a++) {
    cout << "xaxis[" << a << "] " << xaxis[a] << endl;
  }
  // done with definining x-axis

  string rap = "";
  string pt = "6.5-30";
  string cent = "0-100";
  cout << rap << " " << pt << " " << cent << endl;
  plotMultipleSamples(heff_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb);

  heff_rap.clear();
  delete[] xaxis;
  for (unsigned int i=0; i<isize; i++) {
    delete fsource[i];
  }
}




void ppPbPbRatio(vector<string> *inputfile, string *outname, vector<string> *histname, bool isacc) {
  
  const unsigned int isize = inputfile->size();
  vector<drawingEff *> fsource;
  for (unsigned int i=0; i<isize; i++) {
    cout << "creating drawingEff: " << inputfile->at(i) << endl;
    //i: 0th : pp, 1st: pbpb
    fsource.push_back(new drawingEff(inputfile->at(i).c_str(),i,isacc));
    fsource[i]->loadHisto();
    fsource[i]->getEfficiency();
  }

  string raparr[] = {"0-2.4","0-0.6","0.6-1.2","1.2-1.8","1.8-2.4","1.8-2.4"};
  string ptarr[] = {"6.5-30","3-6.5"};
  string centarr[] = {"0-10","10-30","30-100"};

  vector<TGraphAsymmErrors*> heff_cent_rap;
  vector<TGraphAsymmErrors*> heff_pt_rap;
  vector<TGraphAsymmErrors*> heff_pt_cent;
  vector<TGraphAsymmErrors*> heff_rap;
  
  bool ispbpb=true;
  //xVar is for drawing latex on plots
  //xVar == 0 : vs. Cent
  //xVar == 1 : vs. pt in 4+1 |y| regions
  //xVar == 2 : vs. pt in 3 cent regions
  //xVar == 3 : vs. rap integrated over all regions
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  unsigned int heff_size = fsource[0]->heff_cent_rap.size();
  cout << "heff_size: " << heff_size << endl;
  for (unsigned int i=0; i<heff_size; i++) { // heff_size == number of loops for each 4+1 |y| region
    int xVar = 0;
    for (unsigned int fidx=0; fidx<isize; fidx++) {
      heff_cent_rap.push_back(fsource[fidx]->heff_cent_rap[i]); // Compare fdix input files
    }
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_cent_rap[i];
    const int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = raparr[i];
    string pt = "6.5-30";
    string cent = "";
    if (i+1 == heff_size) {
      rap = "1.8-2.4";
      pt = "3-6.5";
    }
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_cent_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb, isacc);
    
    heff_cent_rap.clear();
    delete[] xaxis;
  }
  // Eff vs pT in 4+1 |y| regions
  heff_size = fsource[0]->heff_pt_rap.size();
  for (unsigned int i=0; i<heff_size; i++) {
    int xVar = 1;
    for (unsigned int fidx=0; fidx<isize; fidx++) {
      heff_pt_rap.push_back(fsource[fidx]->heff_pt_rap[i]);
    }
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_pt_rap[i];
    int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = raparr[i];
    string pt = "";
    string cent = "0-100";
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_pt_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb, isacc);

    heff_pt_rap.clear();
    delete[] xaxis;
  }
  // Eff vs pT in 3 centrality regions
  heff_size = fsource[0]->heff_pt_cent.size();
  cout << "heff_size " << heff_size << endl;
  for (unsigned int i=0; i<heff_size; i++) {
    int xVar = 2;
    heff_pt_cent.push_back(fsource[0]->heff_pt_cent[0]); // pp doesn't have centrality differential histos
    heff_pt_cent.push_back(fsource[1]->heff_pt_cent[i]);
    
    // Defining x-axis to re-making ratio histogram
    TH1F *hden = fsource[0]->hden_pt_cent[i];
    int nbins = hden->GetNbinsX();
    double *xaxis = new double[nbins+1];
    for (int a=0; a<nbins; a++) {
      xaxis[a] = hden->GetBinLowEdge(a+1);
    }
    xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
    for (int a=0; a<nbins+1; a++) {
      cout << "xaxis[" << a << "] " << xaxis[a] << endl;
    }
    // done with definining x-axis

    string rap = "0-2.4";
    string pt = "";
    string cent = centarr[i];
    cout << rap << " " << pt << " " << cent << endl;
    plotMultipleSamples(heff_pt_cent, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb, isacc);

    heff_pt_cent.clear();
    delete[] xaxis;
  }
  // Eff vs rap integrated
  int xVar = 3;
  for (unsigned int fidx=0; fidx<isize; fidx++) {
    heff_rap.push_back(fsource[fidx]->heff_rap);
  }
  // Defining x-axis to re-making ratio histogram
  TH1F *hden = fsource[0]->hden_rap;
  int nbins = hden->GetNbinsX();
  double *xaxis = new double[nbins+1];
  for (int a=0; a<nbins; a++) {
    xaxis[a] = hden->GetBinLowEdge(a+1);
  }
  xaxis[nbins] = xaxis[nbins-1] + hden->GetBinWidth(nbins);
  for (int a=0; a<nbins+1; a++) {
    cout << "xaxis[" << a << "] " << xaxis[a] << endl;
  }
  // done with definining x-axis

  string rap = "";
  string pt = "6.5-30";
  string cent = "0-100";
  cout << rap << " " << pt << " " << cent << endl;
  plotMultipleSamples(heff_rap, histname, outname, xVar, xaxis, rap, pt, cent, ispbpb, isacc);

  heff_rap.clear();
  delete[] xaxis;
  for (unsigned int i=0; i<isize; i++) {
    delete fsource[i];
  }
}

void drawMultiples() {
  string dir = "figs";
  string subdir = "ratios";

  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  vector<string> sourcedir;
  sourcedir.push_back("noSF");
  sourcedir.push_back("nominal");
  sourcedir.push_back("trg__muid");
  sourcedir.push_back("trg__sta");

  vector<string> histname;
  histname.push_back("No TnP SF");
  histname.push_back("Nominal (Trg)");
  histname.push_back("Trg + MuID");
  histname.push_back("Trg + STA");

  vector<string> inputfile;
  string outname;
  bool ispbpb;
  
  for (size_t i=0; i<sourcedir.size(); i++) {
    inputfile.push_back(Form("files/eff/%s/histos_jpsi_pp.root",sourcedir[i].c_str()));
    outname = dir + "/" + subdir + "/jpsi_pp";
    ispbpb = false;
  }
  callMultiples(&inputfile, ispbpb, &outname, &histname);
  inputfile.clear();

  for (size_t i=0; i<sourcedir.size(); i++) {
    inputfile.push_back(Form("files/eff/%s/histos_npjpsi_pp.root",sourcedir[i].c_str()));
    outname = dir + "/" + subdir + "/npjpsi_pp";
    ispbpb = false;
  }
  callMultiples(&inputfile, ispbpb, &outname, &histname);
  inputfile.clear();

  for (size_t i=0; i<sourcedir.size(); i++) {
    inputfile.push_back(Form("files/eff/%s/histos_jpsi_pbpb.root",sourcedir[i].c_str()));
    outname = dir + "/" + subdir + "/jpsi_pbpb";
    ispbpb = true;
  }
  callMultiples(&inputfile, ispbpb, &outname, &histname);
  inputfile.clear();

  for (size_t i=0; i<sourcedir.size(); i++) {
    inputfile.push_back(Form("files/eff/%s/histos_npjpsi_pbpb.root",sourcedir[i].c_str()));
    outname = dir + "/" + subdir + "/npjpsi_pbpb";
    ispbpb = true;
  }
  callMultiples(&inputfile, ispbpb, &outname, &histname);
  inputfile.clear();


}

void drawRatios() {
  string dir = "figs";
  string subdir = "ratios_pr_eff";
  bool isacc=false;

  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  vector<string> inputfile;
  vector<string> histname;
  string outname;
 
  // input files should be given in the pp -> pbpb order 
  inputfile.push_back("files/eff/nominal/histos_jpsi_pp.root");
  inputfile.push_back("files/eff/nominal/histos_jpsi_pbpb.root");
  //inputfile.push_back("files/acc/nominal/histos_jpsi_pp.root");
  //inputfile.push_back("files/acc/nominal/histos_jpsi_pbpb.root");
  histname.push_back("pp");
  histname.push_back("PbPb");

  outname = dir + "/" + subdir ;
  
  ppPbPbRatio(&inputfile, &outname, &histname, isacc);
  inputfile.clear();

}

