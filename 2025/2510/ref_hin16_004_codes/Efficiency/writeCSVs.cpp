#include "oniaEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <fstream>

// Bin boundaries for efficiency histograms
const double bins_4rap[] = {0, 0.6, 1.2, 1.8, 2.4};
const int nbins_4rap = sizeof(bins_4rap)/sizeof(double) -1;
const double bins_3cent[] = {0, 10, 30, 100};
const int nbins_3cent = sizeof(bins_3cent)/sizeof(double) -1;

using namespace HI;
using namespace std;

class getEffSyst{
  public:
    // HIN-16-004 binning
    struct h16004 {
      TH1F *hnum_centfwd;
      TH1F *hden_centfwd;
      TGraphAsymmErrors *heff_centfwd;
      TProfile *peff_centfwd;
      TH1F *hnum_ptfwd;
      TH1F *hden_ptfwd;
      TGraphAsymmErrors *heff_ptfwd;
      TProfile *peff_ptfwd;
      TH1F *hnum_centmid;
      TH1F *hden_centmid;
      TGraphAsymmErrors *heff_centmid;
      TProfile *peff_centmid;
      TH1F *hnum_ptmid;
      TH1F *hden_ptmid;
      TGraphAsymmErrors *heff_ptmid;
      TProfile *peff_ptmid;
    };

    struct h16025 {   
      // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
      vector<TH1F*> hnum_cent_rap;
      vector<TH1F*> hden_cent_rap;
      vector<TGraphAsymmErrors*> heff_cent_rap;
      vector<TProfile*> peff_cent_rap;
      // Eff vs pT in 4+1 |y| regions
      vector<TH1F*> hnum_pt_rap;
      vector<TH1F*> hden_pt_rap;
      vector<TGraphAsymmErrors*> heff_pt_rap;
      vector<TProfile*> peff_pt_rap;
      // Eff vs pT in 3 centrality regions
      vector<TH1F*> hnum_pt_cent;
      vector<TH1F*> hden_pt_cent;
      vector<TGraphAsymmErrors*> heff_pt_cent;
      vector<TProfile*> peff_pt_cent;
      // Eff vs rap integrated
      TH1F *hnum_rap;
      TH1F *hden_rap;
      TGraphAsymmErrors *heff_rap;
      TProfile *peff_rap;
      // Eff all integrated (calculated by hnum_rap, hden_rap)
      TH1F *hnum_all;
      TH1F *hden_all;
      TGraphAsymmErrors *heff_all;
      TProfile *peff_all;
    };

    getEffSyst(string fnameSyst, bool ispbpb, bool takemax, bool is16025, bool isacc=false);
    getEffSyst(string fname, string fnameSyst, bool ispbpb, bool takemax, bool is16025, bool isacc=false);
    getEffSyst(string fname, string fnameSyst, string fnameSyst2, bool ispbpb, bool takemax, bool is16025, bool isacc=false);
    virtual ~getEffSyst();
    void runStat(bool isToy=true);
    void runSyst();
    virtual void writeStat(vector<string> *outname);
    virtual void writeStatToy(vector<string> *outname);
    virtual void writeSyst(vector<string> *outname);


  private:
    // Input file
    TFile *finput, *finputSyst1, *finputSyst2;

    bool ispbpb, takemax, is16025, isacc, istoy;

    // nom is nominal
    // syst1 is used for syst with 100 toys or stat
    // syst2 is for the 2nd systematic histograms when get the max(syst1,syst2)
    h16004 nom_04, syst1_04, syst2_04;
    h16025 nom_25, syst1_25, syst2_25;
    
    void clean16004(h16004 &syst);
    void clean16025(h16025 &syst);
    void loadHisto16004(TFile *f, h16004 &h);
    void loadHisto16025(TFile *f, h16025 &h);
    void getEfficiency16004(h16004 &h);
    void getEfficiency16025(h16025 &h);
    void loadProfile16004(TFile *f, h16004 &h);
    void loadProfile16025(TFile *f, h16025 &h);
    
    void fixCentPp(TH1F* h);
    void fixCentPp(TH1F& h);
    void fixErrorNan(TGraphAsymmErrors *geff, TH1F *hden, TH1F *hnum);
    void checkUnderFlow(TH1 *hnum, TH1 *hden);
    TProfile* allIntegratedEff(TObjArray *anum, TObjArray *aden);
    TProfile* simpleEff(TObjArray *anum, TObjArray *aden, bool iscent);
    TProfile* newTProfile(TH1F *hist);
    void fillTProf(TProfile *tp, TH1F *th);


};

void getEffSyst::runStat(bool isToy_) {
  istoy = isToy_;
  if (istoy) {
    if (is16025) loadProfile16025(finputSyst1, syst1_25);
    else loadProfile16004(finputSyst1, syst1_04);
  } else {
    if (is16025) {
      loadHisto16025(finputSyst1, syst1_25);
      getEfficiency16025(syst1_25);
    } else {
      loadHisto16004(finputSyst1, syst1_04);
      getEfficiency16004(syst1_04);
    }
  }
}

void getEffSyst::runSyst() {
  istoy = false;
  if (is16025) {
    loadHisto16025(finput, nom_25);
    loadHisto16025(finputSyst1, syst1_25);
    getEfficiency16025(nom_25);
    getEfficiency16025(syst1_25);
    if (finputSyst2) {
      loadHisto16025(finputSyst2, syst2_25);
      getEfficiency16025(syst2_25);
    }
  } else {
    loadHisto16004(finput, nom_04);
    loadHisto16004(finputSyst1, syst1_04);
    getEfficiency16004(nom_04);
    getEfficiency16004(syst1_04);
    if (finputSyst2) {
      loadHisto16004(finputSyst2, syst2_04);
      getEfficiency16004(syst2_04);
    }
  }

}

void getEffSyst::fixCentPp(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   float y = hist->GetBinContent(1);
   float dy = hist->GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist->SetBinContent(i,y);
      hist->SetBinError(i,dy);
   }
}
void getEffSyst::fixCentPp(TH1F &hist) {
   int nbins = hist.GetNbinsX();
   float y = hist.GetBinContent(1);
   float dy = hist.GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist.SetBinContent(i,y);
      hist.SetBinError(i,dy);
   }
}

TProfile* getEffSyst::newTProfile(TH1F *hist) {
  int nbins = hist->GetNbinsX();
  double *xbins = new double[nbins+1];
  xbins[0] = hist->GetXaxis()->GetBinLowEdge(1);
  for (int i=1; i<nbins+1; i++) xbins[i] = hist->GetXaxis()->GetBinUpEdge(i);
  return new TProfile(Form("tprof_%s",hist->GetName()),"",nbins,xbins,"S");
}

void getEffSyst::fillTProf(TProfile *tp, TH1F *th) {
  if (tp->GetNbinsX() != th->GetNbinsX()) {
    cout << "Error, tp and th have different number of bins" << endl;
    return;
  }
  for (int i=0; i<=tp->GetNbinsX()+1; i++) {
    tp->Fill(th->GetBinCenter(i),th->GetBinContent(i));
  }
}

TProfile* getEffSyst::allIntegratedEff(TObjArray *anum, TObjArray *aden) {
  // get all integrated bin's efficiency from _rap histograms
  TProfile *tprof = new TProfile("tprof_peff_all","",1,0,1,"S");

  for (int i=0; i<anum->GetSize(); i++) {
    TH1F *hnum = (TH1F*) anum->At(i);
    TH1F *hden = (TH1F*) aden->At(i);
    checkUnderFlow(hnum,hden);

    double numerr, denerr;
    double numval = hnum->IntegralAndError(1,hnum->GetNbinsX(),numerr);
    double denval = hden->IntegralAndError(1,hden->GetNbinsX(),denerr);
    
    TH1F *hnum_all = new TH1F("hnum_all","",1,0,1);
    TH1F *hden_all = new TH1F("hden_all","",1,0,1);
    hnum_all->SetBinContent(1,numval);
    hden_all->SetBinContent(1,denval);
    hnum_all->SetBinError(1,numerr);
    hden_all->SetBinError(1,denerr);

    TH1F *hratio = new TH1F("peff_all","",1,0,1);
    hratio->Divide(hnum_all,hden_all);
    fillTProf(tprof,hratio);
    delete hnum;
    delete hden;
    delete hnum_all;
    delete hden_all;
    delete hratio;
  }

  return tprof;
}

TProfile* getEffSyst::simpleEff(TObjArray *anum, TObjArray *aden, bool iscent) {
  TProfile *tprof = newTProfile((TH1F*) anum->At(0));

  for (int i=0; i<anum->GetSize(); i++) {
    TH1F *hnum = (TH1F*) anum->At(i);
    TH1F *hden = (TH1F*) aden->At(i);
    if (!ispbpb && iscent) { // set all bin contents to same as the 1st bin
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

getEffSyst::getEffSyst(string fnameSyst, bool ispbpb_, bool takemax_, bool is16025_, bool isacc_) {
  finputSyst1 = new TFile(fnameSyst.c_str());
  finput = 0;
  finputSyst2 = 0;
  ispbpb = ispbpb_;
  takemax = takemax_;
  is16025 = is16025_;
  isacc = isacc_;
}

getEffSyst::getEffSyst(string fname, string fnameSyst, bool ispbpb_, bool takemax_, bool is16025_, bool isacc_) {
  finput = new TFile(fname.c_str());
  finputSyst1 = new TFile(fnameSyst.c_str());
  finputSyst2 = 0;
  ispbpb = ispbpb_;
  takemax = takemax_;
  is16025 = is16025_;
  isacc = isacc_;
}

getEffSyst::getEffSyst(string fname, string fnameSyst, string fnameSyst2, bool ispbpb_, bool takemax_, bool is16025_, bool isacc_) {
  finput = new TFile(fname.c_str());
  finputSyst1 = new TFile(fnameSyst.c_str());
  finputSyst2 = new TFile(fnameSyst2.c_str());
  ispbpb = ispbpb_;
  takemax = takemax_;
  is16025 = is16025_;
  isacc = isacc_;
}

getEffSyst::~getEffSyst() {
  if (finput) {
    if (is16025) clean16025(nom_25);
    else clean16004(nom_04);
    delete finput;
  }
  if (finputSyst1) {
    if (is16025) clean16025(syst1_25);
    else clean16004(syst1_04);
    delete finputSyst1;
  }
  if (finputSyst2) {
    if (is16025) clean16025(syst2_25);
    else clean16004(syst2_04);
    delete finputSyst2;
  }
}

void getEffSyst::clean16004(h16004 &syst) {
  // hin-16-004 binning 
  if (istoy) {
    if (syst.peff_centfwd) delete syst.peff_centfwd;
    if (syst.peff_ptfwd) delete syst.peff_ptfwd;
    if (syst.peff_centmid) delete syst.peff_centmid;
    if (syst.peff_ptmid) delete syst.peff_ptmid;
  } else {
    if (syst.hnum_centfwd) delete syst.hnum_centfwd;
    if (syst.hden_centfwd) delete syst.hden_centfwd;
    if (syst.heff_centfwd) delete syst.heff_centfwd;
    if (syst.hnum_ptfwd) delete syst.hnum_ptfwd;
    if (syst.hden_ptfwd) delete syst.hden_ptfwd;
    if (syst.heff_ptfwd) delete syst.heff_ptfwd;
    if (syst.hnum_centmid) delete syst.hnum_centmid;
    if (syst.hden_centmid) delete syst.hden_centmid;
    if (syst.heff_centmid) delete syst.heff_centmid;
    if (syst.hnum_ptmid) delete syst.hnum_ptmid;
    if (syst.hden_ptmid) delete syst.hden_ptmid;
    if (syst.heff_ptmid) delete syst.heff_ptmid;
  }
}

void getEffSyst::clean16025(h16025 &syst) {
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  for (size_t i=0; i<syst.hnum_cent_rap.size(); i++) 
    if (syst.hnum_cent_rap[i]) delete syst.hnum_cent_rap[i];
  for (size_t i=0; i<syst.hden_cent_rap.size(); i++) 
    if (syst.hden_cent_rap[i]) delete syst.hden_cent_rap[i];
  for (size_t i=0; i<syst.heff_cent_rap.size(); i++) 
    if (syst.heff_cent_rap[i]) delete syst.heff_cent_rap[i];
  for (size_t i=0; i<syst.peff_cent_rap.size(); i++) 
    if (syst.peff_cent_rap[i]) delete syst.peff_cent_rap[i];
  
  // Eff vs pT in 4+1 |y| regions
  for (size_t i=0; i<syst.hnum_pt_rap.size(); i++) 
    if (syst.hnum_pt_rap[i]) delete syst.hnum_pt_rap[i];
  for (size_t i=0; i<syst.hden_pt_rap.size(); i++) 
    if (syst.hden_pt_rap[i]) delete syst.hden_pt_rap[i];
  for (size_t i=0; i<syst.heff_pt_rap.size(); i++) 
    if (syst.heff_pt_rap[i]) delete syst.heff_pt_rap[i];
  for (size_t i=0; i<syst.peff_pt_rap.size(); i++) 
    if (syst.peff_pt_rap[i]) delete syst.peff_pt_rap[i];

  // Eff vs pT in 3 centrality regions
  for (size_t i=0; i<syst.hnum_pt_cent.size(); i++) 
    if (syst.hnum_pt_cent[i]) delete syst.hnum_pt_cent[i];
  for (size_t i=0; i<syst.hden_pt_cent.size(); i++) 
    if (syst.hden_pt_cent[i]) delete syst.hden_pt_cent[i];
  for (size_t i=0; i<syst.heff_pt_cent.size(); i++) 
    if (syst.heff_pt_cent[i]) delete syst.heff_pt_cent[i];
  for (size_t i=0; i<syst.peff_pt_cent.size(); i++) 
    if (syst.peff_pt_cent[i]) delete syst.peff_pt_cent[i];

  // Eff vs rap integrated
  if (syst.hnum_pt_cent.size() != 0) {
    if (syst.hnum_rap) delete syst.hnum_rap;
    if (syst.hden_rap) delete syst.hden_rap;
    if (syst.heff_rap) delete syst.heff_rap;
  }
  if (syst.peff_pt_cent.size() != 0) 
    if (syst.peff_rap) delete syst.peff_rap;
  
  // Eff all integrated (calculated by hnum_rap, hden_rap)
  if (syst.hnum_pt_cent.size() != 0) {
    if (syst.hnum_all) delete syst.hnum_all;
    if (syst.hden_all) delete syst.hden_all;
    if (syst.heff_all) delete syst.heff_all;
  }
  if (syst.peff_pt_cent.size() != 0) 
    if (syst.peff_all) delete syst.peff_all;
}


void getEffSyst::checkUnderFlow(TH1 *hnum, TH1 *hden){
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


void getEffSyst::loadProfile16025(TFile *f, h16025 &h) {
  if (!f) cout << "Cannot open input file!\n" << endl;

  // Load denominator, numerators
  TObjArray *hnum, *hden;
  TProfile *peff;

  bool iscent = true;
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  hnum = (TObjArray*)f->Get("hnum_cent_rap0024");
  hden = (TObjArray*)f->Get("hden_cent_rap0024");
  peff = simpleEff(hnum, hden, iscent);
  peff->SetName("peff_cent_rap0024");
  h.peff_cent_rap.push_back(peff);
 
  for (int i=0; i<nbins_4rap; i++) {
    hnum = (TObjArray*)f->Get( Form("hnum_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    hden = (TObjArray*)f->Get( Form("hden_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    peff = simpleEff(hnum, hden, iscent);
    peff->SetName(Form("peff_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10));
    h.peff_cent_rap.push_back(peff);
  }
 
  hnum = (TObjArray*)f->Get("hnum_cent_rap1824_pt3065");
  hden = (TObjArray*)f->Get("hden_cent_rap1824_pt3065");
  peff = simpleEff(hnum, hden, iscent);
  peff->SetName("peff_cent_rap1824_pt3065");
  h.peff_cent_rap.push_back(peff);

  iscent = false;
  // Eff vs pT in 4+1 |y| regions
  hnum = (TObjArray*)f->Get("hnum_pt_rap0024");
  hden = (TObjArray*)f->Get("hden_pt_rap0024");
  peff = simpleEff(hnum, hden, iscent);
  peff->SetName("peff_pt_rap0024");
  h.peff_pt_rap.push_back(peff);
  
  for (int i=0; i<nbins_4rap; i++) {
    hnum = (TObjArray*)f->Get( Form("hnum_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    hden = (TObjArray*)f->Get( Form("hden_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    peff = simpleEff(hnum, hden, iscent);
    peff->SetName(Form("peff_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10));
    h.peff_pt_rap.push_back(peff);
  }

  // Eff vs pT in 3 centrality regions
  for (int i=0; i<(ispbpb?nbins_3cent-1:1); i++) { // for pp, only the 1st bin is necessary
    hnum = (TObjArray*)f->Get( Form("hnum_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    hden = (TObjArray*)f->Get( Form("hden_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    peff = simpleEff(hnum, hden, iscent);
    peff->SetName(Form("peff_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]));
    h.peff_pt_cent.push_back(peff);
  }
  if (ispbpb) {
    hnum = (TObjArray*)f->Get( Form("hnum_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    hden = (TObjArray*)f->Get( Form("hden_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    peff = simpleEff(hnum, hden, iscent);
    peff->SetName(Form("peff_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]));
    h.peff_pt_cent.push_back(peff);
  }

  // Eff vs rap integrated
  hnum = (TObjArray*)f->Get("hnum_rap");
  hden = (TObjArray*)f->Get("hden_rap");
  h.peff_rap = simpleEff(hnum, hden, iscent);
  h.peff_rap->SetName("peff_rap");
  
  // Eff all integrated (calculated by hnum_rap, hden_rap)
  hnum = (TObjArray*)f->Get("hnum_rap");
  hden = (TObjArray*)f->Get("hden_rap");
  h.peff_all = allIntegratedEff(hnum, hden);
  h.peff_all->SetName("peff_all");

}

void getEffSyst::loadProfile16004(TFile *f, h16004 &h) {
  TObjArray *hnum, *hden;
  TProfile *peff;
  string gname;

  bool iscent = true;
  hnum = (TObjArray*)f->Get("hnum_cent_rap0016"); 
  hden = (TObjArray*)f->Get("hden_cent_rap0016"); 
  h.peff_centmid = simpleEff(hnum, hden, iscent);
  h.peff_centmid->SetName("peff_cent_rap0016");

  hnum = (TObjArray*)f->Get("hnum_cent_rap1624"); 
  hden = (TObjArray*)f->Get("hden_cent_rap1624"); 
  h.peff_centfwd = simpleEff(hnum, hden, iscent);
  h.peff_centfwd->SetName("peff_cent_rap1624");
  
  iscent = false;
  hnum = (TObjArray*)f->Get("hnum_pt_rap0016"); 
  hden = (TObjArray*)f->Get("hden_pt_rap0016"); 
  h.peff_ptmid = simpleEff(hnum, hden, iscent);
  h.peff_ptmid->SetName("peff_pt_rap0016");

  hnum = (TObjArray*)f->Get("hnum_pt_rap1624"); 
  hden = (TObjArray*)f->Get("hden_pt_rap1624"); 
  h.peff_ptfwd = simpleEff(hnum, hden, iscent);
  h.peff_ptfwd->SetName("peff_pt_rap1624");
}


void getEffSyst::loadHisto16025(TFile *f, h16025 &h) {
  if (!f) cout << "Cannot open input file!\n" << endl;

  // Load denominator, numerators
  // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
  TH1F *hnum = (TH1F*)f->Get("hnum_cent_rap0024");
  TH1F *hden = (TH1F*)f->Get("hden_cent_rap0024");
  if (!ispbpb || isacc) { // set all bin contents to same as the 1st bin
    fixCentPp(hnum);
    fixCentPp(hden);
  }
  checkUnderFlow(hnum,hden);
  h.hnum_cent_rap.push_back(hnum);
  h.hden_cent_rap.push_back(hden);
  for (int i=0; i<nbins_4rap; i++) {
    TH1F *hnum1 = (TH1F*)f->Get( Form("hnum_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TH1F *hden1 = (TH1F*)f->Get( Form("hden_cent_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    if (!ispbpb || isacc) { // set all bin contents to same as the 1st bin
      fixCentPp(hnum1);
      fixCentPp(hden1);
    }
    checkUnderFlow(hnum1,hden1);
    h.hnum_cent_rap.push_back(hnum1);
    h.hden_cent_rap.push_back(hden1);
  }
  hnum = (TH1F*)f->Get("hnum_cent_rap1824_pt3065");
  hden = (TH1F*)f->Get("hden_cent_rap1824_pt3065");
  if (!ispbpb || isacc) { // set all bin contents to same as the 1st bin
    fixCentPp(hnum);
    fixCentPp(hden);
  }
  checkUnderFlow(hnum,hden);
  h.hnum_cent_rap.push_back(hnum);
  h.hden_cent_rap.push_back(hden);

  // Eff vs pT in 4+1 |y| regions
  hnum = (TH1F*)f->Get("hnum_pt_rap0024");
  hden = (TH1F*)f->Get("hden_pt_rap0024");
  checkUnderFlow(hnum,hden);
  h.hnum_pt_rap.push_back(hnum);
  h.hden_pt_rap.push_back(hden);
  for (int i=0; i<nbins_4rap; i++) {
    TH1F *hnum1 = (TH1F*)f->Get( Form("hnum_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    TH1F *hden1 = (TH1F*)f->Get( Form("hden_pt_rap%02.0f%02.0f",bins_4rap[i]*10,bins_4rap[i+1]*10) );
    checkUnderFlow(hnum1,hden1);
    h.hnum_pt_rap.push_back(hnum1);
    h.hden_pt_rap.push_back(hden1);
  }

  // Eff vs pT in 3 centrality regions
  for (int i=0; i<(ispbpb?nbins_3cent-1:1); i++) { // for pp, only the 1st bin is necessary
    TH1F *hnum1 = (TH1F*)f->Get( Form("hnum_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    TH1F *hden1 = (TH1F*)f->Get( Form("hden_pt_cent%02.0f%02.0f",bins_3cent[i],bins_3cent[i+1]) );
    checkUnderFlow(hnum1,hden1);
    h.hnum_pt_cent.push_back(hnum1);
    h.hden_pt_cent.push_back(hden1);
  }
  if (ispbpb) {
    hnum = (TH1F*)f->Get( Form("hnum_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    hden = (TH1F*)f->Get( Form("hden_pt_cent%02.0f%03.0f",bins_3cent[nbins_3cent-1],bins_3cent[nbins_3cent]) );
    checkUnderFlow(hnum,hden);
    h.hnum_pt_cent.push_back(hnum);
    h.hden_pt_cent.push_back(hden);
  }

  // Eff vs rap integrated
  h.hnum_rap = (TH1F*)f->Get("hnum_rap");
  h.hden_rap = (TH1F*)f->Get("hden_rap");
  checkUnderFlow(h.hnum_rap,h.hden_rap);
  
  // Eff all integrated (calculated by hnum_rap, hden_rap)
  double numerr, denerr;
  double numval = h.hnum_rap->IntegralAndError(1,h.hnum_rap->GetNbinsX(),numerr);
  double denval = h.hden_rap->IntegralAndError(1,h.hden_rap->GetNbinsX(),denerr);
  int i=0;
  if (finput==f) i=0;
  else if (finputSyst1==f) i=1;
  else if (finputSyst2==f) i=2;
  h.hnum_all = new TH1F(Form("hnum_all_%d",i),"",1,0,1);
  h.hden_all = new TH1F(Form("hden_all_%d",i),"",1,0,1);
  h.hnum_all->SetBinContent(1,numval);
  h.hden_all->SetBinContent(1,denval);
  h.hnum_all->SetBinError(1,numerr);
  h.hden_all->SetBinError(1,denerr);
}

void getEffSyst::loadHisto16004(TFile *f, h16004 &h) {
  h.hnum_centmid = (TH1F*)f->Get("hnum_cent_rap0016"); 
  h.hnum_centfwd = (TH1F*)f->Get("hnum_cent_rap1624"); 
  h.hden_centmid = (TH1F*)f->Get("hden_cent_rap0016"); 
  h.hden_centfwd = (TH1F*)f->Get("hden_cent_rap1624"); 
  if (!ispbpb || isacc) { // set all bin contents to same as the 1st bin
    fixCentPp(h.hnum_centmid);
    fixCentPp(h.hnum_centfwd);
    fixCentPp(h.hden_centmid);
    fixCentPp(h.hden_centfwd);
  }
  checkUnderFlow(h.hnum_centmid,h.hden_centmid);
  checkUnderFlow(h.hnum_centfwd,h.hden_centfwd);
  h.hnum_ptmid = (TH1F*)f->Get("hnum_pt_rap0016"); 
  h.hnum_ptfwd = (TH1F*)f->Get("hnum_pt_rap1624"); 
  h.hden_ptmid = (TH1F*)f->Get("hden_pt_rap0016"); 
  h.hden_ptfwd = (TH1F*)f->Get("hden_pt_rap1624"); 
  checkUnderFlow(h.hnum_ptmid,h.hden_ptmid);
  checkUnderFlow(h.hnum_ptfwd,h.hden_ptfwd);
}

void getEffSyst::fixErrorNan(TGraphAsymmErrors *geff, TH1F *hden, TH1F *hnum) {
  const int np = geff->GetN();
  for (int a=0; a<np; a++) {
    double gey = geff->GetErrorYhigh(a);
    double gex = geff->GetErrorXhigh(a);
    if (gey!=gey) { // when error is nan, error will be calculated by hand!
      double eyd = hden->GetBinError(a+1);
      double eyn = hnum->GetBinError(a+1);
      double yd = hden->GetBinContent(a+1);
      double yn = hnum->GetBinContent(a+1);
      double uncert = yn/yd * sqrt( pow(eyd/yd,2) + pow(eyn/yn,2) ) / 2.0;
      geff->SetPointError(a,gex,gex,uncert,uncert);
    }
  }
}

void getEffSyst::getEfficiency16025(h16025 &h) {
  // Calculate efficiency
  for (int i=0; i<nbins_4rap+2; i++) {
    TGraphAsymmErrors *geff = new TGraphAsymmErrors(h.hnum_cent_rap[i],h.hden_cent_rap[i],"n");
    fixErrorNan(geff,h.hnum_cent_rap[i],h.hden_cent_rap[i]);
    h.heff_cent_rap.push_back(geff);
    string gname = h.hnum_cent_rap[i]->GetName();
    gname.replace(gname.begin(),gname.begin()+4,"heff");
    h.heff_cent_rap[i]->SetName(gname.c_str());
  }
  for (int i=0; i<nbins_4rap+1; i++) {
    TGraphAsymmErrors *geff = new TGraphAsymmErrors(h.hnum_pt_rap[i],h.hden_pt_rap[i],"n");
    fixErrorNan(geff,h.hnum_pt_rap[i],h.hden_pt_rap[i]);
    h.heff_pt_rap.push_back(geff);
    string gname = h.hnum_pt_rap[i]->GetName();
    gname.replace(gname.begin(),gname.begin()+4,"heff");
    h.heff_pt_rap[i]->SetName(gname.c_str());
  }
  for (int i=0; i<(ispbpb?nbins_3cent:1); i++) {
    TGraphAsymmErrors *geff = new TGraphAsymmErrors(h.hnum_pt_cent[i],h.hden_pt_cent[i],"n");
    fixErrorNan(geff,h.hnum_pt_cent[i],h.hden_pt_cent[i]);
    h.heff_pt_cent.push_back(geff);
    string gname = h.hnum_pt_cent[i]->GetName();
    gname.replace(gname.begin(),gname.begin()+4,"heff");
    h.heff_pt_cent[i]->SetName(gname.c_str());
  }
  h.heff_rap = new TGraphAsymmErrors(h.hnum_rap,h.hden_rap,"n");
  fixErrorNan(h.heff_rap,h.hnum_rap,h.hden_rap);
  string gname = h.hnum_rap->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_rap->SetName(gname.c_str());
  
  h.heff_all = new TGraphAsymmErrors(h.hnum_all,h.hden_all,"n");
  fixErrorNan(h.heff_all,h.hnum_all,h.hden_all);
  gname = h.hnum_all->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_all->SetName(gname.c_str());
}

void getEffSyst::getEfficiency16004(h16004 &h) {
  string gname;

  h.heff_centmid = new TGraphAsymmErrors(h.hnum_centmid,h.hden_centmid,"n");
  fixErrorNan(h.heff_centmid,h.hnum_centmid,h.hden_centmid);
  gname = h.hnum_centmid->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_centmid->SetName(gname.c_str());

  h.heff_ptmid = new TGraphAsymmErrors(h.hnum_ptmid,h.hden_ptmid,"n");
  fixErrorNan(h.heff_ptmid,h.hnum_ptmid,h.hden_ptmid);
  gname = h.hnum_ptmid->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_ptmid->SetName(gname.c_str());

  h.heff_centfwd = new TGraphAsymmErrors(h.hnum_centfwd,h.hden_centfwd,"n");
  fixErrorNan(h.heff_centfwd,h.hnum_centfwd,h.hden_centfwd);
  gname = h.hnum_centfwd->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_centfwd->SetName(gname.c_str());

  h.heff_ptfwd = new TGraphAsymmErrors(h.hnum_ptfwd,h.hden_ptfwd,"n");
  fixErrorNan(h.heff_ptfwd,h.hnum_ptfwd,h.hden_ptfwd);
  gname = h.hnum_ptfwd->GetName();
  gname.replace(gname.begin(),gname.begin()+4,"heff");
  h.heff_ptfwd->SetName(gname.c_str());
}



void getEffSyst::writeSyst(vector<string> *outname) {
  string foutdir = outname->at(1);
  string foutname = outname->at(2);
  cout << foutdir << "/" <<foutname << endl;
  gSystem->mkdir(Form("%s",foutdir.c_str()),kTRUE);
  ofstream foSyst(Form("%s/%s",foutdir.c_str(),foutname.c_str()), ios::out);

  foSyst << outname->at(0) << endl;
 
  /////// relative uncertainties bin by bin
  double *effSystY, *effSystY2, *effY;
  
  if (is16025) {
    // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
    for (int i=0; i<=nbins_4rap+1; i++) {
      effSystY = syst1_25.heff_cent_rap[i]->GetY();
      effY = nom_25.heff_cent_rap[i]->GetY();
      effSystY2 = (takemax ? syst2_25.heff_cent_rap[i]->GetY() : effY);

      for (int j=0; j<nom_25.hnum_cent_rap[i]->GetNbinsX(); j++) {
        // centrality bin index is 0-200, so multiply 2
        double lowEdge = nom_25.hnum_cent_rap[i]->GetBinLowEdge(j+1) *2;
        double binWidth = nom_25.hnum_cent_rap[i]->GetBinWidth(j+1) *2;
        double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
        if (i==0) foSyst << "0, 2.4, 6.5, 50, ";
        else if (i==nbins_4rap+1) foSyst << "1.8, 2.4, 3, 6.5, ";
        else foSyst << bins_4rap[i-1] << ", " << bins_4rap[i] << ", 6.5, 50, ";
        if (ispbpb) foSyst << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
        else { foSyst << "0, 200, " << relSyst << endl; break; }
      }
    } // end of 1 line
   
    // Eff vs pT in 4+1 |y| regions
    for (int i=0; i<=nbins_4rap; i++) {
      effSystY = syst1_25.heff_pt_rap[i]->GetY();
      effY = nom_25.heff_pt_rap[i]->GetY();
      effSystY2 = (takemax ? syst2_25.heff_pt_rap[i]->GetY() : effY);
      
      for (int j=0; j<nom_25.hnum_pt_rap[i]->GetNbinsX(); j++) {
        double lowEdge = nom_25.hnum_pt_rap[i]->GetBinLowEdge(j+1);
        double binWidth = nom_25.hnum_pt_rap[i]->GetBinWidth(j+1);
        double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
        if (i==0) foSyst << "0, 2.4, " ;
        else foSyst << bins_4rap[i-1] << ", " << bins_4rap[i] << ", " ;
        foSyst << lowEdge << ", " << lowEdge+binWidth << ", 0, 200, " << relSyst << endl;
      }
    } // end of 1 line

    // Eff vs pT in 3 centrality regions
    for (int i=0; i<(ispbpb?nbins_3cent:1); i++) {
      effSystY = syst1_25.heff_pt_cent[i]->GetY();
      effY = nom_25.heff_pt_cent[i]->GetY();
      effSystY2 = (takemax ? syst2_25.heff_pt_cent[i]->GetY() : effY);
      
      for (int j=0; j<nom_25.hnum_pt_cent[i]->GetNbinsX(); j++) {
        double lowEdge = nom_25.hnum_pt_cent[i]->GetBinLowEdge(j+1);
        double binWidth = nom_25.hnum_pt_cent[i]->GetBinWidth(j+1);
        double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
        foSyst << "0, 2.4, " ;
        foSyst << lowEdge << ", " << lowEdge+binWidth << ", ";
        if (ispbpb) foSyst << bins_3cent[i]*2 << ", " << bins_3cent[i+1]*2 << ", " << relSyst << endl;
        else foSyst << "0, 200, " << relSyst << endl;
      }
    } // end of 1 line

    // Eff vs rap integrated
    effSystY = syst1_25.heff_rap->GetY();
    effY = nom_25.heff_rap->GetY();
    effSystY2 = (takemax ? syst2_25.heff_rap->GetY() : effY);
    for (int j=0; j<nom_25.hnum_rap->GetNbinsX(); j++) {
      double lowEdge = nom_25.hnum_rap->GetBinLowEdge(j+1);
      double binWidth = nom_25.hnum_rap->GetBinWidth(j+1);
      double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
      foSyst << lowEdge << ", " << lowEdge+binWidth << ", ";
      foSyst << "6.5, 50, 0, 200, " << relSyst << endl;
    }
    
    if (ispbpb) {
      // Eff all integrated (calculated by hnum_rap, hden_rap)
      effSystY = syst1_25.heff_all->GetY();
      effY = nom_25.heff_all->GetY();
      effSystY2 = (takemax ? syst2_25.heff_all->GetY() : effY);
      for (int j=0; j<nom_25.hnum_all->GetNbinsX(); j++) {
        double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
        foSyst << "0, 2.4, 6.5, 50, 0, 200, " << relSyst << endl;
      }
    }

  } else {
    // HIN-16-004 binning
    effSystY = syst1_04.heff_centmid->GetY();
    effY = nom_04.heff_centmid->GetY();
    effSystY2 = (takemax ? syst2_04.heff_centmid->GetY() : effY);
    for (int j=0; j<nom_04.hnum_centmid->GetNbinsX(); j++) {
      double lowEdge = nom_04.hnum_centmid->GetBinLowEdge(j+1) *2;
      double binWidth = nom_04.hnum_centmid->GetBinWidth(j+1) *2;
      double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
      foSyst << "0, 1.6, 6.5, 30, ";
      if (ispbpb) foSyst << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
      else { foSyst << "0, 200, " << relSyst << endl; break; }
    }

    effSystY = syst1_04.heff_centfwd->GetY();
    effY = nom_04.heff_centfwd->GetY();
    effSystY2 = (takemax ? syst2_04.heff_centfwd->GetY() : effY);
    for (int j=0; j<nom_04.hnum_centfwd->GetNbinsX(); j++) {
      double lowEdge = nom_04.hnum_centfwd->GetBinLowEdge(j+1) *2;
      double binWidth = nom_04.hnum_centfwd->GetBinWidth(j+1) *2;
      double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
      foSyst << "1.6, 2.4, 3, 30, ";
      if (ispbpb) foSyst << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
      else { foSyst << "0, 200, " << relSyst << endl; break; }
    }

    effSystY = syst1_04.heff_ptmid->GetY();
    effY = nom_04.heff_ptmid->GetY();
    effSystY2 = (takemax ? syst2_04.heff_ptmid->GetY() : effY);
    for (int j=0; j<nom_04.hnum_ptmid->GetNbinsX(); j++) {
      double lowEdge = nom_04.hnum_ptmid->GetBinLowEdge(j+1);
      double binWidth = nom_04.hnum_ptmid->GetBinWidth(j+1);
      double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
      foSyst << "0, 1.6, ";
      foSyst << lowEdge << ", " << lowEdge+binWidth << ", ";
      foSyst << "0, 200, " << relSyst << endl;
    }

    effSystY = syst1_04.heff_ptfwd->GetY();
    effY = nom_04.heff_ptfwd->GetY();
    effSystY2 = (takemax ? syst2_04.heff_ptfwd->GetY() : effY);
    for (int j=0; j<nom_04.hnum_ptfwd->GetNbinsX(); j++) {
      double lowEdge = nom_04.hnum_ptfwd->GetBinLowEdge(j+1);
      double binWidth = nom_04.hnum_ptfwd->GetBinWidth(j+1);
      double relSyst = max( TMath::Abs(effSystY[j]-effY[j]),TMath::Abs(effSystY2[j]-effY[j]) )/effY[j];
      foSyst << "1.6, 2.4, ";
      foSyst << lowEdge << ", " << lowEdge+binWidth << ", ";
      foSyst << "0, 200, " << relSyst << endl;
    }
  }

  foSyst.close();
}

void getEffSyst::writeStat(vector<string> *outname) {
  string foutdir = outname->at(1);
  string foutname = outname->at(2);
  cout << foutdir << "/" <<foutname << endl;
  gSystem->mkdir(Form("%s",foutdir.c_str()),kTRUE);
  ofstream foStat(Form("%s/%s",foutdir.c_str(),foutname.c_str()), ios::out);

  foStat << outname->at(0) << endl;
  
  /////// relative uncertainties bin by bin
  double *effY, *effEY;
  
  if (is16025) {
    // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
    for (int i=0; i<=nbins_4rap+1; i++) {
      effY = syst1_25.heff_cent_rap[i]->GetY();
      effEY = syst1_25.heff_cent_rap[i]->GetEYhigh();

      for (int j=0; j<syst1_25.hnum_cent_rap[i]->GetNbinsX(); j++) {
        // centrality bin index is 0-200, so multiply 2
        double lowEdge = syst1_25.hnum_cent_rap[i]->GetBinLowEdge(j+1) *2;
        double binWidth = syst1_25.hnum_cent_rap[i]->GetBinWidth(j+1) *2;
        double relSyst = TMath::Abs(effEY[j]/effY[j]);
        if (i==0) foStat << "0, 2.4, 6.5, 50, ";
        else if (i==nbins_4rap+1) foStat << "1.8, 2.4, 3, 6.5, ";
        else foStat << bins_4rap[i-1] << ", " << bins_4rap[i] << ", 6.5, 50, ";
        if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
        else { foStat << "0, 200, " << relSyst << endl; break; }
      }
    } // end of 1 line
   
    // Eff vs pT in 4+1 |y| regions
    for (int i=0; i<=nbins_4rap; i++) {
      effY = syst1_25.heff_pt_rap[i]->GetY();
      effEY = syst1_25.heff_pt_rap[i]->GetEYhigh();
      
      for (int j=0; j<syst1_25.hnum_pt_rap[i]->GetNbinsX(); j++) {
        double lowEdge = syst1_25.hnum_pt_rap[i]->GetBinLowEdge(j+1);
        double binWidth = syst1_25.hnum_pt_rap[i]->GetBinWidth(j+1);
        double relSyst = TMath::Abs(effEY[j]/effY[j]);
        if (i==0) foStat << "0, 2.4, " ;
        else foStat << bins_4rap[i-1] << ", " << bins_4rap[i] << ", " ;
        foStat << lowEdge << ", " << lowEdge+binWidth << ", 0, 200, " << relSyst << endl;
      }
    } // end of 1 line

    // Eff vs pT in 3 centrality regions
    for (int i=0; i<(ispbpb?nbins_3cent:1); i++) {
      int idx = i;
      if (isacc) idx = 0; //force the 0th histogram to read acceptance 

      effY = syst1_25.heff_pt_cent[idx]->GetY();
      effEY = syst1_25.heff_pt_cent[idx]->GetEYhigh();
      
      for (int j=0; j<syst1_25.hnum_pt_cent[idx]->GetNbinsX(); j++) {
        double lowEdge = syst1_25.hnum_pt_cent[idx]->GetBinLowEdge(j+1);
        double binWidth = syst1_25.hnum_pt_cent[idx]->GetBinWidth(j+1);
        double relSyst = TMath::Abs(effEY[j]/effY[j]);
        foStat << "0, 2.4, " ;
        foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
        if (ispbpb) foStat << bins_3cent[idx]*2 << ", " << bins_3cent[idx+1]*2 << ", " << relSyst << endl;
        else foStat << "0, 200, " << relSyst << endl;
      }
    } // end of 1 line

    // Eff vs rap integrated
    effY = syst1_25.heff_rap->GetY();
    effEY = syst1_25.heff_rap->GetEYhigh();
    
    for (int j=0; j<syst1_25.hnum_rap->GetNbinsX(); j++) {
      double lowEdge = syst1_25.hnum_rap->GetBinLowEdge(j+1);
      double binWidth = syst1_25.hnum_rap->GetBinWidth(j+1);
      double relSyst = TMath::Abs(effEY[j]/effY[j]);
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "6.5, 50, 0, 200, " << relSyst << endl;
    }
    
    if (ispbpb) {
      // Eff all integrated (calculated by hnum_rap, hden_rap)
      effY = syst1_25.heff_all->GetY();
      effEY = syst1_25.heff_all->GetEYhigh();

      for (int j=0; j<syst1_25.hnum_all->GetNbinsX(); j++) {
        double relSyst = TMath::Abs(effEY[j]/effY[j]);
        foStat << "0, 2.4, 6.5, 50, 0, 200, " << relSyst << endl;
      }
    }

  } else {
    // HIN-16-004 binning
    effY = syst1_04.heff_centmid->GetY();
    effEY = syst1_04.heff_centmid->GetEYhigh();

    for (int j=0; j<syst1_04.hnum_centmid->GetNbinsX(); j++) {
      double lowEdge = syst1_04.hnum_centmid->GetBinLowEdge(j+1) *2;
      double binWidth = syst1_04.hnum_centmid->GetBinWidth(j+1) *2;
      double relSyst = TMath::Abs(effEY[j]/effY[j]);
      foStat << "0, 1.6, 6.5, 30, ";
      if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
      else { foStat << "0, 200, " << relSyst << endl; break; }
    }

    effY = syst1_04.heff_centfwd->GetY();
    effEY = syst1_04.heff_centfwd->GetEYhigh();

    for (int j=0; j<syst1_04.hnum_centfwd->GetNbinsX(); j++) {
      double lowEdge = syst1_04.hnum_centfwd->GetBinLowEdge(j+1) *2;
      double binWidth = syst1_04.hnum_centfwd->GetBinWidth(j+1) *2;
      double relSyst = TMath::Abs(effEY[j]/effY[j]);
      foStat << "1.6, 2.4, 3, 30, ";
      if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relSyst << endl;
      else { foStat << "0, 200, " << relSyst << endl; break; }
    }

    effY = syst1_04.heff_ptmid->GetY();
    effEY = syst1_04.heff_ptmid->GetEYhigh();

    for (int j=0; j<syst1_04.hnum_ptmid->GetNbinsX(); j++) {
      double lowEdge = syst1_04.hnum_ptmid->GetBinLowEdge(j+1);
      double binWidth = syst1_04.hnum_ptmid->GetBinWidth(j+1);
      double relSyst = TMath::Abs(effEY[j]/effY[j]);
      foStat << "0, 1.6, ";
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "0, 200, " << relSyst << endl;
    }

    effY = syst1_04.heff_ptfwd->GetY();
    effEY = syst1_04.heff_ptfwd->GetEYhigh();

    for (int j=0; j<syst1_04.hnum_ptfwd->GetNbinsX(); j++) {
      double lowEdge = syst1_04.hnum_ptfwd->GetBinLowEdge(j+1);
      double binWidth = syst1_04.hnum_ptfwd->GetBinWidth(j+1);
      double relSyst = TMath::Abs(effEY[j]/effY[j]);
      foStat << "1.6, 2.4, ";
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "0, 200, " << relSyst << endl;
    }
  }

  foStat.close();
}

void getEffSyst::writeStatToy(vector<string> *outname) {
  string foutdir = outname->at(1);
  string foutname = outname->at(2);
  cout << foutdir << "/" <<foutname << endl;
  gSystem->mkdir(Form("%s",foutdir.c_str()),kTRUE);
  ofstream foStat(Form("%s/%s",foutdir.c_str(),foutname.c_str()), ios::out);

  foStat << outname->at(0) << endl;
 
  /////// relative uncertainties bin by bin
  if (is16025) {
    // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
    for (int i=0; i<=nbins_4rap+1; i++) {
      for (int j=0; j<syst1_25.peff_cent_rap[i]->GetNbinsX(); j++) {
        double effStatY = syst1_25.peff_cent_rap[i]->GetBinError(j+1);
        double effY = syst1_25.peff_cent_rap[i]->GetBinContent(j+1);
        // centrality bin index is 0-200, so multiply 2
        double lowEdge = syst1_25.peff_cent_rap[i]->GetBinLowEdge(j+1) *2;
        double binWidth = syst1_25.peff_cent_rap[i]->GetBinWidth(j+1) *2;
        double relStat = TMath::Abs(effStatY/effY);
        if (i==0) foStat << "0, 2.4, 6.5, 50, ";
        else if (i==nbins_4rap+1) foStat << "1.8, 2.4, 3, 6.5, ";
        else foStat << bins_4rap[i-1] << ", " << bins_4rap[i] << ", 6.5, 50, ";
        // pp has bin contents in the 1st bin only
        if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relStat << endl;
        else { foStat << "0, 200, " << relStat << endl; break; }
      }
    } // end of 1 line
   
    // Eff vs pT in 4+1 |y| regions
    for (int i=0; i<=nbins_4rap; i++) {
      for (int j=0; j<syst1_25.peff_pt_rap[i]->GetNbinsX(); j++) {
        double effStatY = syst1_25.peff_pt_rap[i]->GetBinError(j+1);
        double effY = syst1_25.peff_pt_rap[i]->GetBinContent(j+1);
        double lowEdge = syst1_25.peff_pt_rap[i]->GetBinLowEdge(j+1);
        double binWidth = syst1_25.peff_pt_rap[i]->GetBinWidth(j+1);
        double relStat = TMath::Abs(effStatY/effY);
        if (i==0) foStat << "0, 2.4, " ;
        else foStat << bins_4rap[i-1] << ", " << bins_4rap[i] << ", " ;
        foStat << lowEdge << ", " << lowEdge+binWidth << ", 0, 200, " << relStat << endl;
      }
    } // end of 1 line

    // Eff vs pT in 3 centrality regions
    for (int i=0; i<(ispbpb?nbins_3cent:1); i++) {
      for (int j=0; j<syst1_25.peff_pt_cent[i]->GetNbinsX(); j++) {
        double effStatY = syst1_25.peff_pt_cent[i]->GetBinError(j+1);
        double effY = syst1_25.peff_pt_cent[i]->GetBinContent(j+1);
        double lowEdge = syst1_25.peff_pt_cent[i]->GetBinLowEdge(j+1);
        double binWidth = syst1_25.peff_pt_cent[i]->GetBinWidth(j+1);
        double relStat = TMath::Abs(effStatY/effY);
        foStat << "0, 2.4, " ;
        foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
        if (ispbpb) foStat << bins_3cent[i]*2 << ", " << bins_3cent[i+1]*2 << ", " << relStat << endl;
        else foStat << "0, 200, " << relStat << endl;
      }
    } // end of 1 line

    // Eff vs rap integrated
    for (int j=0; j<syst1_25.peff_rap->GetNbinsX(); j++) {
      double effStatY = syst1_25.peff_rap->GetBinError(j+1);
      double effY = syst1_25.peff_rap->GetBinContent(j+1);
      double lowEdge = syst1_25.peff_rap->GetBinLowEdge(j+1);
      double binWidth = syst1_25.peff_rap->GetBinWidth(j+1);
      double relStat = TMath::Abs(effStatY/effY);
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "6.5, 50, 0, 200, " << relStat << endl;
    }
    
    if (ispbpb) {
      // Eff all integrated (calculated by hnum_rap, hden_rap)
      for (int j=0; j<syst1_25.peff_all->GetNbinsX(); j++) {
        double effStatY = syst1_25.peff_all->GetBinError(j+1);
        double effY = syst1_25.peff_all->GetBinContent(j+1);
        double relStat = TMath::Abs(effStatY/effY);
        foStat << "0, 2.4, 6.5, 50, 0, 200, " << relStat << endl;
      }
    }

  } else {
    // HIN-16-004 binning
    for (int j=0; j<syst1_04.peff_centmid->GetNbinsX(); j++) {
      double effStatY = syst1_04.peff_centmid->GetBinError(j+1);
      double effY = syst1_04.peff_centmid->GetBinContent(j+1);
      double lowEdge = syst1_04.peff_centmid->GetBinLowEdge(j+1) *2;
      double binWidth = syst1_04.peff_centmid->GetBinWidth(j+1) *2;
      double relStat = TMath::Abs(effStatY/effY);
      foStat << "0, 1.6, 6.5, 30, ";
      if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relStat << endl;
      else { foStat << "0, 200, " << relStat << endl; break; }
    }

    for (int j=0; j<syst1_04.peff_centfwd->GetNbinsX(); j++) {
      double effStatY = syst1_04.peff_centfwd->GetBinError(j+1);
      double effY = syst1_04.peff_centfwd->GetBinContent(j+1);
      double lowEdge = syst1_04.peff_centfwd->GetBinLowEdge(j+1) *2;
      double binWidth = syst1_04.peff_centfwd->GetBinWidth(j+1) *2;
      double relStat = TMath::Abs(effStatY/effY);
      foStat << "1.6, 2.4, 3, 30, ";
      if (ispbpb) foStat << lowEdge << ", " << lowEdge+binWidth << ", " << relStat << endl;
      else { foStat << "0, 200, " << relStat << endl; break; }
    }

    for (int j=0; j<syst1_04.peff_ptmid->GetNbinsX(); j++) {
      double effStatY = syst1_04.peff_ptmid->GetBinError(j+1);
      double effY = syst1_04.peff_ptmid->GetBinContent(j+1);
      double lowEdge = syst1_04.peff_ptmid->GetBinLowEdge(j+1);
      double binWidth = syst1_04.peff_ptmid->GetBinWidth(j+1);
      double relStat = TMath::Abs(effStatY/effY);
      foStat << "0, 1.6, ";
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "0, 200, " << relStat << endl;
    }

    for (int j=0; j<syst1_04.peff_ptfwd->GetNbinsX(); j++) {
      double effStatY = syst1_04.peff_ptfwd->GetBinError(j+1);
      double effY = syst1_04.peff_ptfwd->GetBinContent(j+1);
      double lowEdge = syst1_04.peff_ptfwd->GetBinLowEdge(j+1);
      double binWidth = syst1_04.peff_ptfwd->GetBinWidth(j+1);
      double relStat = TMath::Abs(effStatY/effY);
      foStat << "1.6, 2.4, ";
      foStat << lowEdge << ", " << lowEdge+binWidth << ", ";
      foStat << "0, 200, " << relStat << endl;
    }
  }

  foStat.close();
}

////////////////////////////////////////////
///////////////// M A I N  //////////////////
/////////////////////////////////////////////
void writeCSVs(bool is16025 = true) {

  string dir = "syst";
  string subdir = "eff";

  string statdirsToy[] = {"trg_ptWeighting","trg_toy","trg__muid_toy","trg__sta_toy"};
  string systdirs[] = {"trg_binned"};
  string systdirsMax[] = {
  "trg_minus1sigma","trg_plus1sigma",
  "trg_trk_plus1sigma","trg_trk_minus1sigma",
  "trg__muid_plus1sigma","trg__muid_minus1sigma",
  "trg__sta_plus1sigma","trg__sta_minus1sigma"};

  string statdirsToy_name[] = {"ptWeighting","trg_stat","muid_stat","sta_stat"};
  string systdirsMax_name[] = {"trg_syst","trk_syst","muid_syst","sta_syst"};

  const int nstatdirToy = sizeof(statdirsToy)/sizeof(string);
  const int nsystdir = sizeof(systdirs)/sizeof(string);
  const int nsystdirMax = sizeof(systdirsMax)/sizeof(string);
  
  string nameTag = "16004";
  if (is16025) nameTag = "16025";
  
  // treat statistical efficiencies 
  dir = "stat";
  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);
/* 
  for (int i=0; i<2; i++) {
    bool ispbpb=false;
    bool takemax=false;
    bool isToy = false;
    bool isacc = i==0 ? true : false;
    subdir = i==0? "acc" : "eff";
    
    vector<string> latex;
    latex.push_back("nominal");
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_prompt_PP_%s.csv",nameTag.c_str(),subdir.c_str()));
    getEffSyst pp_pr(
        Form("files/%s/nominal/histos_jpsi_pp.root",subdir.c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_pr.runStat(isToy);
    pp_pr.writeStat(&latex);
    latex.clear();

    latex.push_back("nominal");
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_nonprompt_PP_%s.csv",nameTag.c_str(),subdir.c_str()));
    getEffSyst pp_np(
        Form("files/%s/nominal/histos_npjpsi_pp.root",subdir.c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_np.runStat(isToy);
    pp_np.writeStat(&latex);
    latex.clear();

    ispbpb = true;

    latex.push_back("nominal");
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_prompt_PbPb_%s.csv",nameTag.c_str(),subdir.c_str()));
    getEffSyst pbpb_pr(
        Form("files/%s/nominal/histos_jpsi_pbpb.root",subdir.c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_pr.runStat(isToy);
    pbpb_pr.writeStat(&latex);
    latex.clear();

    latex.push_back("nominal");
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_nonprompt_PbPb_%s.csv",nameTag.c_str(),subdir.c_str()));
    getEffSyst pbpb_np(
        Form("files/%s/nominal/histos_npjpsi_pbpb.root",subdir.c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_np.runStat(isToy);
    pbpb_np.writeStat(&latex);
  }
*/
  // treat statistical efficiencies with 100 Toys
  dir = "stat";
  subdir = "eff";
  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);
  
  for (int i=0; i<nstatdirToy; i++) {
    bool ispbpb=false;
    bool takemax=false;
    bool isToy = true;
    bool isacc = false;
    
    vector<string> latex;
    latex.push_back(statdirsToy_name[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_prompt_PP_%s.csv",nameTag.c_str(),statdirsToy_name[i].c_str()));
    getEffSyst pp_pr(
        Form("files/%s/%s/histos_jpsi_pp.root",subdir.c_str(),statdirsToy[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_pr.runStat(isToy);
    pp_pr.writeStatToy(&latex);
    latex.clear();

    latex.push_back(statdirsToy_name[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_nonprompt_PP_%s.csv",nameTag.c_str(),statdirsToy_name[i].c_str()));
    getEffSyst pp_np(
        Form("files/%s/%s/histos_npjpsi_pp.root",subdir.c_str(),statdirsToy[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_np.runStat(isToy);
    pp_np.writeStatToy(&latex);
    latex.clear();

    ispbpb=true;

    latex.push_back(statdirsToy_name[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_prompt_PbPb_%s.csv",nameTag.c_str(),statdirsToy_name[i].c_str()));
    getEffSyst pbpb_pr(
        Form("files/%s/%s/histos_jpsi_pbpb.root",subdir.c_str(),statdirsToy[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_pr.runStat(isToy);
    pbpb_pr.writeStatToy(&latex);
    latex.clear();

    latex.push_back(statdirsToy_name[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("stat_%s_NJpsi_nonprompt_PbPb_%s.csv",nameTag.c_str(),statdirsToy_name[i].c_str()));
    getEffSyst pbpb_np(
        Form("files/%s/%s/histos_npjpsi_pbpb.root",subdir.c_str(),statdirsToy[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_np.runStat(isToy);
    pbpb_np.writeStatToy(&latex);
  }
  
  
  // treat systematical efficiencies
  dir = "syst";
  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  for (int i=0; i<nsystdir; i++) {
    bool ispbpb=false;
    bool takemax=false;
    bool isacc = false;

    vector<string> latex;
    latex.push_back(systdirs[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_prompt_PP_%s.csv",nameTag.c_str(),systdirs[i].c_str()));
    getEffSyst pp_pr(
        Form("files/%s/nominal/histos_jpsi_pp.root",subdir.c_str()),
        Form("files/%s/%s/histos_jpsi_pp.root",subdir.c_str(),systdirs[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_pr.runSyst();
    pp_pr.writeSyst(&latex);
    latex.clear();

    latex.push_back(systdirs[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_nonprompt_PP_%s.csv",nameTag.c_str(),systdirs[i].c_str()));
    getEffSyst pp_np(
        Form("files/%s/nominal/histos_npjpsi_pp.root",subdir.c_str()),
        Form("files/%s/%s/histos_npjpsi_pp.root",subdir.c_str(),systdirs[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_np.runSyst();
    pp_np.writeSyst(&latex);
    latex.clear();

    ispbpb = true;

    latex.push_back(systdirs[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_prompt_PbPb_%s.csv",nameTag.c_str(),systdirs[i].c_str()));
    getEffSyst pbpb_pr(
        Form("files/%s/nominal/histos_jpsi_pbpb.root",subdir.c_str()),
        Form("files/%s/%s/histos_jpsi_pbpb.root",subdir.c_str(),systdirs[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_pr.runSyst();
    pbpb_pr.writeSyst(&latex);
    latex.clear();

    latex.push_back(systdirs[i].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_nonprompt_PbPb_%s.csv",nameTag.c_str(),systdirs[i].c_str()));
    getEffSyst pbpb_np(
        Form("files/%s/nominal/histos_npjpsi_pbpb.root",subdir.c_str()),
        Form("files/%s/%s/histos_npjpsi_pbpb.root",subdir.c_str(),systdirs[i].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_np.runSyst();
    pbpb_np.writeSyst(&latex);
  }
 

  // treat systematical efficiencies (maximum among 2 sets will be taken)
  dir = "syst";
  gSystem->mkdir(dir.c_str(),kTRUE);
  gSystem->mkdir(Form("%s/%s",dir.c_str(),subdir.c_str()),kTRUE);

  for (int i=0; i<nsystdirMax; i=i+2) {
    bool ispbpb=false;
    bool takemax=true;
    bool isacc = false;

    vector<string> latex;
    latex.push_back(systdirsMax_name[i/2].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_prompt_PP_%s.csv",nameTag.c_str(),systdirsMax_name[i/2].c_str()));
    getEffSyst pp_pr(
        Form("files/%s/nominal/histos_jpsi_pp.root",subdir.c_str()),
        Form("files/%s/%s/histos_jpsi_pp.root",subdir.c_str(),systdirsMax[i].c_str()),
        Form("files/%s/%s/histos_jpsi_pp.root",subdir.c_str(),systdirsMax[i+1].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_pr.runSyst();
    pp_pr.writeSyst(&latex);
    latex.clear();

    latex.push_back(systdirsMax_name[i/2].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_nonprompt_PP_%s.csv",nameTag.c_str(),systdirsMax_name[i/2].c_str()));
    getEffSyst pp_np(
        Form("files/%s/nominal/histos_npjpsi_pp.root",subdir.c_str()),
        Form("files/%s/%s/histos_npjpsi_pp.root",subdir.c_str(),systdirsMax[i].c_str()),
        Form("files/%s/%s/histos_npjpsi_pp.root",subdir.c_str(),systdirsMax[i+1].c_str()),
        ispbpb,takemax,is16025,isacc);
    pp_np.runSyst();
    pp_np.writeSyst(&latex);
    latex.clear();
    
    ispbpb=true;

    latex.push_back(systdirsMax_name[i/2].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_prompt_PbPb_%s.csv",nameTag.c_str(),systdirsMax_name[i/2].c_str()));
    getEffSyst pbpb_pr(
        Form("files/%s/nominal/histos_jpsi_pbpb.root",subdir.c_str()),
        Form("files/%s/%s/histos_jpsi_pbpb.root",subdir.c_str(),systdirsMax[i].c_str()),
        Form("files/%s/%s/histos_jpsi_pbpb.root",subdir.c_str(),systdirsMax[i+1].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_pr.runSyst();
    pbpb_pr.writeSyst(&latex);
    latex.clear();
    
    latex.push_back(systdirsMax_name[i/2].c_str());
    latex.push_back(Form("%s/%s/",dir.c_str(),subdir.c_str()));
    latex.push_back(Form("syst_%s_NJpsi_nonprompt_PbPb_%s.csv",nameTag.c_str(),systdirsMax_name[i/2].c_str()));
    getEffSyst pbpb_np(
        Form("files/%s/nominal/histos_npjpsi_pbpb.root",subdir.c_str()),
        Form("files/%s/%s/histos_npjpsi_pbpb.root",subdir.c_str(),systdirsMax[i].c_str()),
        Form("files/%s/%s/histos_npjpsi_pbpb.root",subdir.c_str(),systdirsMax[i+1].c_str()),
        ispbpb,takemax,is16025,isacc);
    pbpb_np.runSyst();
    pbpb_np.writeSyst(&latex);
  }


}
