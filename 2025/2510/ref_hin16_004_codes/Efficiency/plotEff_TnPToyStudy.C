#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include <TProfile.h>
#include <TObjArray.h>

#include "../Fitter/Macros/CMS/CMS_lumi.C"
#include "../Fitter/Macros/CMS/tdrstyle.C"
#include "../Fitter/Macros/Utilities/texUtils.h"
#include "JpsiDrawing.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


class drawingEff_TnPToyStudy {
public:

  // Bin boundaries for efficiency histograms
  vector<float> bins_4rap;
  vector<float> bins_3cent;
  int nbins_4rap;
  int nbins_3cent;
  bool ispbpb;

  // Input file
  TFile *finput;

  drawingEff_TnPToyStudy(string fname, bool ispbpb);
  ~drawingEff_TnPToyStudy();
  void loadHisto(vector<string> *outname);
  TProfile* simpleEff(TObjArray *anum, TObjArray *aden, bool dofixcent);
  TProfile* newTProfile(TH1F *hist);
  void fillTProf(TProfile *tp, TH1F *th);
  void fixCentPp(TH1F*);
  void checkUnderFlow(TH1 *hnum, TH1 *hden);
};

void drawingEff_TnPToyStudy::fixCentPp(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   float y = hist->GetBinContent(1);
   float dy = hist->GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist->SetBinContent(i,y);
      hist->SetBinError(i,dy);
   }
}

TProfile* drawingEff_TnPToyStudy::newTProfile(TH1F *hist) {
  int nbins = hist->GetNbinsX();
  double *xbins = new double[nbins+1];
  xbins[0] = hist->GetXaxis()->GetBinLowEdge(1);
  for (int i=1; i<nbins+1; i++) xbins[i] = hist->GetXaxis()->GetBinUpEdge(i);
  return new TProfile(Form("tprof_%s",hist->GetName()),"",nbins,xbins,"S");
}

void drawingEff_TnPToyStudy::fillTProf(TProfile *tp, TH1F *th) {
  if (tp->GetNbinsX() != th->GetNbinsX()) {
    cout << "Error, tp and th have different number of bins" << endl;
    return;
  }
  for (int i=0; i<=tp->GetNbinsX()+1; i++) {
    tp->Fill(th->GetBinCenter(i),th->GetBinContent(i));
  }
}

TProfile* drawingEff_TnPToyStudy::simpleEff(TObjArray *anum, TObjArray *aden, bool dofixcent) {
  setTDRStyle();
  gStyle->SetEndErrorSize(3);

  TProfile *tprof = newTProfile((TH1F*) anum->At(0));

  for (int i=0; i<anum->GetSize(); i++) {
    TH1F *hnum = (TH1F*) anum->At(i);
    TH1F *hden = (TH1F*) aden->At(i);
    if (dofixcent && !ispbpb) { // set all bin contents to same as the 1st bin
      fixCentPp(hnum);
      fixCentPp(hden);
    }
    checkUnderFlow(hnum,hden);
    TH1F *hratio = (TH1F*) hnum->Clone("hratio");
    hratio->Divide(hnum,hden,1,1,"B");
    fillTProf(tprof,hratio);
    delete hnum;
    delete hden;
    delete hratio;
  }

  return tprof;
}


drawingEff_TnPToyStudy::~drawingEff_TnPToyStudy(){
  delete finput;
}

drawingEff_TnPToyStudy::drawingEff_TnPToyStudy(string fname, bool ispbpb_){
  finput = new TFile(fname.c_str());
  ispbpb = ispbpb_;

  float bins_4rap_[] = {0, 0.6, 1.2, 1.8, 2.4};
  float bins_3cent_[] = {0, 10, 30, 100};
  
  nbins_4rap = sizeof(bins_4rap_)/sizeof(float) -1;
  for (int i=0; i<=nbins_4rap; i++)
    bins_4rap.push_back(bins_4rap_[i]);

  nbins_3cent = sizeof(bins_3cent_)/sizeof(float) -1;
  for (int i=0; i<=nbins_3cent; i++) 
    bins_3cent.push_back(bins_3cent_[i]);
}

void drawingEff_TnPToyStudy::checkUnderFlow(TH1 *hnum, TH1 *hden){

//  cout << "\ncheckUnderFlow: " << hnum->GetName() << " " << hden->GetName() ;
  for (int j=0; j<=hnum->GetNbinsX()+1; j++) {
    double num0 = hnum->GetBinContent(j);
    double den0 = hden->GetBinContent(j);
    
    // If underflow bin has more entries in numerator than denominator, set it to 0
    if ((j==0 && num0>den0) || (j==hnum->GetNbinsX()+1 && num0>den0)) {
      hnum->SetBinContent(j,0);
      hden->SetBinContent(j,0);
      hnum->SetBinError(j,0);
      hden->SetBinError(j,0);
    }
    else if (num0>den0) {
      cout << "checkUnderFlow: " << hnum->GetBinContent(j) << "/" << hden->GetBinContent(j) << endl;
    }
  }

}

void drawingEff_TnPToyStudy::loadHisto(vector<string> *outname) {
  if (!finput) cout << "Cannot open input file!\n" << endl;

  TFile *foutput = new TFile(Form("%s_toySFToy.root",outname->at(1).c_str()),"recreate");

  setTDRStyle();
  gStyle->SetEndErrorSize(3);

  TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.035);
  TCanvas *can = new TCanvas("can","can",600,600);
  TLegend *leg = new TLegend(0.55,0.60,0.94,0.88);
  SetLegendStyle(leg);
 
  cout << "\n\t" << outname->at(0) << " " << outname->at(1) << endl;

  // Load denominator, numerators
  ///// Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  TObjArray *hnum = (TObjArray*)finput->Get("hnum_cent_rap0024");
  TObjArray *hden = (TObjArray*)finput->Get("hden_cent_rap0024");
  bool dofixcent = true; // True: Only if it is centrality drawing (no matter if it's pp or pbpb)
  TProfile *tp_cent_rap[5];
  tp_cent_rap[0] = simpleEff(hnum, hden, dofixcent);
  tp_cent_rap[0]->SetName("tnpSFToy_cent_rap0024");
  
  TH1F *haxes = new TH1F("haxes",";;Efficiency",1,0,100);
  haxes->GetXaxis()->SetTitle("Centrality [%]");
  haxes->Draw();

  tp_cent_rap[0]->Draw("same");
  leg->AddEntry(tp_cent_rap[0],"|y|: 0-2.4, 6.5-50 GeV/c","p");
  leg->Draw();
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  can->SaveAs(Form("%s_cent_rap_0.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_cent_rap_0.pdf",outname->at(1).c_str()));

  vector<TH1*> tps;
  vector<string> yname;
 
  yname.push_back("$|y|<2.4$");
  tps.push_back(tp_cent_rap[0]);
  TString texname = outname->at(1) + "_cent.tex";
  inittex(texname.Data(), "Centrality [\\%]", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  delete can;
  delete leg;

  for (int i=0; i<nbins_4rap; i++) {
    TObjArray *hnum1 = (TObjArray*)finput->Get( Form("hnum_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TObjArray *hden1 = (TObjArray*)finput->Get( Form("hden_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    tp_cent_rap[i+1] = simpleEff(hnum1, hden1, dofixcent);
    tp_cent_rap[i+1]->SetName( Form("tnpSFToy_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    
    can = new TCanvas("can","can",600,600);
    leg = new TLegend(0.55,0.60,0.94,0.88);
    SetLegendStyle(leg);
    haxes->Draw();
    tp_cent_rap[i+1]->Draw("same");
    leg->AddEntry(tp_cent_rap[i+1],Form("|y|: %.1f-%.1f, 6.5-50 GeV/c",bins_4rap[i],bins_4rap[i+1]),"p");
    leg->Draw();
    lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
    can->SaveAs(Form("%s_cent_rap_%d.png",outname->at(1).c_str(),i+1));
    can->SaveAs(Form("%s_cent_rap_%d.pdf",outname->at(1).c_str(),i+1));
    
    yname.push_back(Form("$%0.1f<|y|<%0.1f$",bins_4rap[i],bins_4rap[i+1]));
    tps.push_back(tp_cent_rap[i+1]);

    delete can;
    delete leg;
  }
  texname = outname->at(1) + "_cent_pt65500.tex";
  inittex(texname.Data(), "Centrality [\\%]", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();


  hnum = (TObjArray*)finput->Get("hnum_cent_rap1824_pt3065");
  hden = (TObjArray*)finput->Get("hden_cent_rap1824_pt3065");
  tp_cent_rap[4] = simpleEff(hnum, hden, dofixcent);
  tp_cent_rap[4]->SetName("tnpSFToy_cent_rap1824_pt3065");
  
  can = new TCanvas("can","can",600,600);
  leg = new TLegend(0.55,0.60,0.94,0.88);
  SetLegendStyle(leg);
  haxes->Draw();
  tp_cent_rap[4]->Draw("same");
  leg->AddEntry(tp_cent_rap[4],Form("|y|: %.1f-%.1f, 3-6.5 GeV/c",bins_4rap[nbins_4rap-1],bins_4rap[nbins_4rap]),"p");
  leg->Draw();
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  can->SaveAs(Form("%s_cent_rap_%d.png",outname->at(1).c_str(),nbins_4rap));
  can->SaveAs(Form("%s_cent_rap_%d.pdf",outname->at(1).c_str(),nbins_4rap));
  
  yname.push_back(Form("$%0.1f<|y|<%0.1f$, $3<p_{T}<6.5$ GeV/c",bins_4rap[nbins_4rap-1],bins_4rap[nbins_4rap]));
  tps.push_back(tp_cent_rap[4]);

  texname = outname->at(1) + "_cent_pt30500.tex";
  inittex(texname.Data(), "Centrality [\\%]", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  delete can;
  delete leg;  

  ///// Eff vs pT in 4+1 |y| regions
  hnum = (TObjArray*)finput->Get("hnum_pt_rap0024");
  hden = (TObjArray*)finput->Get("hden_pt_rap0024");
  dofixcent = false;
  TProfile *tp_pt_rap[5];
  tp_pt_rap[0] = simpleEff(hnum, hden, dofixcent);
  tp_pt_rap[0]->SetName("tnpSFToy_pt_rap0024");

  haxes->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  haxes->SetBins(1,0,50);
  
  can = new TCanvas("can","can",600,600);
  leg = new TLegend(0.6,0.68,0.9,0.88);
  SetLegendStyle(leg);
  haxes->Draw();
  tp_pt_rap[0]->Draw("same");
  if (ispbpb) leg->AddEntry(tp_pt_rap[0],"|y|: 0-2.4, 0-100%","p");
  else leg->AddEntry(tp_pt_rap[0],"|y|: 0-2.4","p");
  leg->Draw();
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  can->SaveAs(Form("%s_pt_rap_0.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_pt_rap_0.pdf",outname->at(1).c_str()));

  yname.push_back("$|y|<2.4$");
  tps.push_back(tp_pt_rap[0]);
  
  // |y| ranges have different number of pT bins -> create separate tex files
  texname = outname->at(1) + "_pt_rap__rap0024.tex";
  inittex(texname.Data(), "$p_{T}$ $\\mathrm{[GeV/c]}$", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  delete can;
  delete leg;

  for (int i=0; i<nbins_4rap-1; i++) {
    TObjArray *hnum1 = (TObjArray*)finput->Get( Form("hnum_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TObjArray *hden1 = (TObjArray*)finput->Get( Form("hden_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    tp_pt_rap[i+1] = simpleEff(hnum1, hden1, dofixcent);
    tp_pt_rap[i+1]->SetName( Form("tnpSFToy_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    
    can = new TCanvas("can","can",600,600);
    leg = new TLegend(0.6,0.68,0.9,0.88);
    SetLegendStyle(leg);
    haxes->Draw();
    tp_pt_rap[i+1]->Draw("same");
    if (ispbpb) leg->AddEntry(tp_pt_rap[i+1],Form("|y|: %.1f-%.1f, 0-100%%",bins_4rap[i],bins_4rap[i+1]),"p");
    else leg->AddEntry(tp_pt_rap[i+1],Form("|y|: %.1f-%.1f",bins_4rap[i],bins_4rap[i+1]),"p");
    leg->Draw();
    lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
    can->SaveAs(Form("%s_pt_rap_%d.png",outname->at(1).c_str(),i+1));
    can->SaveAs(Form("%s_pt_rap_%d.pdf",outname->at(1).c_str(),i+1));

    yname.push_back(Form("$%0.1f<|y|<%0.1f$",bins_4rap[i],bins_4rap[i+1]));
    tps.push_back(tp_pt_rap[i+1]);

    delete can;
    delete leg;
  }
  texname = outname->at(1) + "_pt_rap__pt65500.tex";
  inittex(texname.Data(), "$p_{T}$ $\\mathrm{[GeV/c]}$", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  for (int i=nbins_4rap-1; i<nbins_4rap; i++) {
    TObjArray *hnum1 = (TObjArray*)finput->Get( Form("hnum_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TObjArray *hden1 = (TObjArray*)finput->Get( Form("hden_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    tp_pt_rap[i+1] = simpleEff(hnum1, hden1, dofixcent);
    tp_pt_rap[i+1]->SetName( Form("tnpSFToy_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    
    can = new TCanvas("can","can",600,600);
    leg = new TLegend(0.6,0.68,0.9,0.88);
    SetLegendStyle(leg);
    haxes->Draw();
    tp_pt_rap[i+1]->Draw("same");
    if (ispbpb) leg->AddEntry(tp_pt_rap[i+1],Form("|y|: %.1f-%.1f, 0-100%%",bins_4rap[i],bins_4rap[i+1]),"p");
    else leg->AddEntry(tp_pt_rap[i+1],Form("|y|: %.1f-%.1f",bins_4rap[i],bins_4rap[i+1]),"p");
    leg->Draw();
    lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
    can->SaveAs(Form("%s_pt_rap_%d.png",outname->at(1).c_str(),i+1));
    can->SaveAs(Form("%s_pt_rap_%d.pdf",outname->at(1).c_str(),i+1));

    yname.push_back(Form("$%0.1f<|y|<%0.1f$",bins_4rap[i],bins_4rap[i+1]));
    tps.push_back(tp_pt_rap[i+1]);

    delete can;
    delete leg;
  }
  texname = outname->at(1) + "_pt_rap__pt30500.tex";
  inittex(texname.Data(), "$p_{T}$ $\\mathrm{[GeV/c]}$", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  ///// Eff vs pT in 3 centrality regions
  TProfile *tp_pt_cent[nbins_3cent];
  for (int i=0; i<(ispbpb?nbins_3cent-1:1); i++) { // for pp, only the 1st bin is necessary
    TObjArray *hnum1 = (TObjArray*)finput->Get( Form("hnum_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    TObjArray *hden1 = (TObjArray*)finput->Get( Form("hden_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    tp_pt_cent[i] = simpleEff(hnum1, hden1, dofixcent);
    tp_pt_cent[i]->SetName( Form("tnpSFToy_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    
    haxes->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    haxes->SetBins(1,0,50);
    
    can = new TCanvas("can","can",600,600);
    leg = new TLegend(0.67,0.68,0.9,0.88);
    SetLegendStyle(leg);
    haxes->Draw();
    tp_pt_cent[i]->Draw("same");
    if (ispbpb) {
      leg->AddEntry(tp_pt_cent[i],Form("%.0f-%.0f%%",bins_3cent[i],bins_3cent[i+1]),"p");
      leg->Draw();
    }
    lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
    lat->DrawLatex(0.2,0.80,"|y| < 2.4");
    can->SaveAs(Form("%s_pt_cent_%d.png",outname->at(1).c_str(),i));
    can->SaveAs(Form("%s_pt_cent_%d.pdf",outname->at(1).c_str(),i));

    yname.push_back(Form("$%0.1f-%0.1f$\\%%",bins_3cent[i],bins_3cent[i+1]));
    tps.push_back(tp_pt_cent[i]);

    delete can;
    delete leg;
  }
  if (ispbpb) {
    hnum = (TObjArray*)finput->Get( Form("hnum_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    hden = (TObjArray*)finput->Get( Form("hden_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    tp_pt_cent[nbins_3cent-1] = simpleEff(hnum, hden, dofixcent);
    tp_pt_cent[nbins_3cent-1]->SetName( Form("tnpSFToy_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    
    haxes->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    haxes->SetBins(1,0,50);
    
    can = new TCanvas("can","can",600,600);
    leg = new TLegend(0.67,0.68,0.9,0.88);
    SetLegendStyle(leg);
    haxes->Draw();
    tp_pt_cent[nbins_3cent-1]->Draw("same");
    leg->AddEntry(tp_pt_cent[nbins_3cent-1],Form("%.0f-%.0f%%",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]),"p");
    leg->Draw();
    lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
    lat->DrawLatex(0.2,0.80,"|y| < 2.4");
    can->SaveAs(Form("%s_pt_cent_%d.png",outname->at(1).c_str(),nbins_3cent));
    can->SaveAs(Form("%s_pt_cent_%d.pdf",outname->at(1).c_str(),nbins_3cent));

    yname.push_back(Form("$%0.1f-%0.1f$\\%%",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]));
    tps.push_back(tp_pt_cent[nbins_3cent-1]);
    
    delete can;
    delete leg;
  }

  texname = outname->at(1) + "_pt_cent.tex";
  inittex(texname.Data(), "$p_{T}$ $\\mathrm{[GeV/c]}$", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  ///// Eff vs rap integrated
  hnum = (TObjArray*)finput->Get("hnum_rap");
  hden = (TObjArray*)finput->Get("hden_rap");
  TProfile *tp_rap = simpleEff(hnum, hden, dofixcent);
  tp_rap->SetName("tnpSFToy_rap");

  haxes->GetXaxis()->SetTitle("|y|");
  haxes->SetBins(1,0,2.4);
    
  can = new TCanvas("can","can",600,600);
  haxes->Draw();
  tp_rap->Draw("same");
  lat->DrawLatex(0.2,0.85,outname->at(0).c_str());
  lat->DrawLatex(0.2,0.80,"6.5-50 GeV/c");
  if (ispbpb) lat->DrawLatex(0.2,0.75,"0-100%");
  can->SaveAs(Form("%s_rap.png",outname->at(1).c_str()));
  can->SaveAs(Form("%s_rap.pdf",outname->at(1).c_str()));

  yname.push_back("$6.5<p_{T}<50$ GeV/c, 0-100%");
  tps.push_back(tp_rap);

  texname = outname->at(1) + "_rap.tex";
  inittex(texname.Data(), "|y|", yname);
  printHist(tps, texname.Data());
  closetex(texname.Data());
  tps.clear();
  yname.clear();

  delete can;
  delete lat;


  ///// Write TProfiles into .root file
  foutput->cd();
  for (int i=0; i<5; i++)
    tp_cent_rap[i]->Write();
  for (int i=0; i<nbins_4rap+1; i++) 
    tp_pt_rap[i]->Write();
  for (int i=0; i<(ispbpb?nbins_3cent:1); i++) 
    tp_pt_cent[i]->Write();
  tp_rap->Write(); 
  

  // Clean up memory
  for (int i=0; i<5; i++)
    delete tp_cent_rap[i];
  for (int i=0; i<nbins_4rap+1; i++) 
    delete tp_pt_rap[i];
  for (int i=0; i<(ispbpb?nbins_3cent:1); i++) 
    delete tp_pt_cent[i];
  delete tp_rap;
}









////////////////////////////////////////////////////////////////
////// Draw eff plots in different kinematic regions
////////////////////////////////////////////////////////////////
void plotEff_TnPToyStudy(string subdir="eff/trg_toy") {

  string dir = "figs";

  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  vector<string> latex;
  latex.push_back("prompt J/#psi (pp)");
  latex.push_back(Form("%s/%s/jpsi_pp",dir.c_str(),subdir.c_str()));
  drawingEff_TnPToyStudy fjpsi_pp(Form("files/%s/histos_jpsi_pp.root",subdir.c_str()),false);
  fjpsi_pp.loadHisto(&latex);

  latex.clear();
  latex.push_back("nonprompt J/#psi (pp)");
  latex.push_back(Form("%s/%s/npjpsi_pp",dir.c_str(),subdir.c_str()));
  drawingEff_TnPToyStudy fnpjpsi_pp(Form("files/%s/histos_npjpsi_pp.root",subdir.c_str()),false);
  fnpjpsi_pp.loadHisto(&latex);

  latex.clear();
  latex.push_back("prompt J/#psi (PbPb)");
  latex.push_back(Form("%s/%s/jpsi_pbpb",dir.c_str(),subdir.c_str()));
  drawingEff_TnPToyStudy fjpsi_pbpb(Form("files/%s/histos_jpsi_pbpb.root",subdir.c_str()),true);
  fjpsi_pbpb.loadHisto(&latex);

  latex.clear();
  latex.push_back("nonprompt J/#psi (PbPb)");
  latex.push_back(Form("%s/%s/npjpsi_pbpb",dir.c_str(),subdir.c_str()));
  drawingEff_TnPToyStudy fnpjpsi_pbpb(Form("files/%s/histos_npjpsi_pbpb.root",subdir.c_str()),true);
  fnpjpsi_pbpb.loadHisto(&latex);
}


