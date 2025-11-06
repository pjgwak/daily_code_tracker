#ifndef results2tree_C
#define results2tree_C

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"

#include <vector>
#include <cstring>

#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/EVENTUTILS.h"

using namespace std;
using namespace RooFit;

struct poi {
  Char_t name[64];
  float val;
  float errH;
  float errL;
  float parIni_mass;
  float parIni_mass_err;
  float parLoad_mass;
  float parLoad_mass_err;
  float parIni_ctaumass;
  float parIni_ctaumass_err;
  float parFit_ctaumass;
  float parFit_ctaumass_err;
  float min;
  float max;
};

const int nBins = 46;

void results2tree(
                  const char* workDirName,
                  const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
                  const char* prependPath, //="",
                  const char* fitType, // "mass", "ctau"...
                  const char* thePoiNames, //="N_Jpsi,b_Jpsi,f_Jpsi,m_Jpsi,sigma1_Jpsi,alpha_Jpsi,n_Jpsi,sigma2_Jpsi,rSigma21_Jpsi,
                  // lambda1_Bkg,lambda2_Bkg,lambda3_Bkg,lambda4_Bkg,lambda5_Bkg,N_Bkg,b_Bkg,
                  // ctau1_CtauRes,ctau2_CtauRes,f_CtauRes,rSigma21_CtauRes,sigma1_CtauRes,
                  // fDFSS_BkgNoPR,fDLIV_BkgNoPR,lambdaDDS_BkgNoPR,lambdaDF_BkgNoPR,lambdaDSS_BkgNoPR,lambdaDSS_JpsiNoPR,
                  // eff,effnp,acc,accnp,lumi,taa,ncoll,npart,correl_N_Jpsi_vs_b_Jpsi",
                  bool wantPureSMC //=false
) {
  // workDirName: usual tag where to look for files in Output
  // thePoiNames: comma-separated list of parameters to store ("par1,par2,par3"). Default: all
  
  TFile *f = new TFile(treeFileName(workDirName,DSTag,prependPath,fitType),"RECREATE");
  TTree *tr = new TTree("fitresults","fit results");
  
  bin<int> centbin_int(0,200); // Integrated bin to use in pp (to avoid centrality [-1,1000])
  
  // bin edges
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  // model names
  Char_t jpsiName[128], bkgName[128];
  // collision system
  Char_t collSystem[8];
  // goodness of fit
  float nll, chi2_invMass, chi2_ctau, chi2prob, chi2prob_invMass, chi2prob_ctau, normchi2_invMass,normchi2_ctau;
  int npar, nparbkg,ndof_invMass, ndof_ctau;
  // parameters to store: make it a vector
  vector<poi> thePois;
  TString thePoiNamesStr(thePoiNames);
  TString t; Int_t from = 0;
  cout << ">>>>>> Parameters to be extracted from workspaces:" << endl;
  while (thePoiNamesStr.Tokenize(t, from , ",")) {
    poi p; strcpy(p.name, t.Data());
    cout << p.name << endl;
    thePois.push_back(p);
  }
  cout << endl;
  
  // create tree branches
  tr->Branch("ptmin",&ptmin,"ptmin/F");
  tr->Branch("ptmax",&ptmax,"ptmax/F");
  tr->Branch("ymin",&ymin,"ymin/F");
  tr->Branch("ymax",&ymax,"ymax/F");
  tr->Branch("centmin",&centmin,"centmin/F");
  tr->Branch("centmax",&centmax,"centmax/F");
  tr->Branch("jpsiName",jpsiName,"jpsiName/C");
  tr->Branch("bkgName",bkgName,"bkgName/C");
  tr->Branch("collSystem",collSystem,"collSystem/C");
  tr->Branch("nll",&nll,"nll/F");
  tr->Branch("chi2_invMass",&chi2_invMass,"chi2_invMass/F");
  tr->Branch("chi2prob_invMass",&chi2prob_invMass,"chi2prob_invMass/F");
  tr->Branch("normchi2_invMass",&normchi2_invMass,"normchi2_invMass/F");
  tr->Branch("chi2prob_invMass",&chi2prob_invMass,"chi2prob_invMass/F");
  tr->Branch("ndof_invMass",&ndof_invMass,"ndof_invMass/I");
  tr->Branch("chi2_ctau",&chi2_ctau,"chi2_ctau/F");
  tr->Branch("chi2prob_ctau",&chi2prob_ctau,"chi2prob_ctau/F");
  tr->Branch("normchi2_ctau",&normchi2_ctau,"normchi2_ctau/F");
  tr->Branch("chi2prob_ctau",&chi2prob_ctau,"chi2prob_ctau/F");
  tr->Branch("ndof_ctau",&ndof_ctau,"ndof_ctau/I");
  tr->Branch("npar",&npar,"npar/I");
  tr->Branch("nparbkg",&nparbkg,"nparbkg/I");
  
  for (vector<poi>::iterator it=thePois.begin(); it!=thePois.end(); it++) {
    tr->Branch(Form("%s_val",it->name),&(it->val),Form("%s_val/F",it->name));
    tr->Branch(Form("%s_errL",it->name),&(it->errL),Form("%s_errL/F",it->name));
    tr->Branch(Form("%s_errH",it->name),&(it->errH),Form("%s_errH/F",it->name));
    tr->Branch(Form("%s_min",it->name),&(it->min),Form("%s_min/F",it->name));
    tr->Branch(Form("%s_max",it->name),&(it->max),Form("%s_max/F",it->name));
    tr->Branch(Form("%s_parIni_mass",it->name),&(it->parIni_mass),Form("%s_parIni_mass/F",it->name));
    tr->Branch(Form("%s_parIni_mass_err",it->name),&(it->parIni_mass_err),Form("%s_parIni_mass_err/F",it->name));
    tr->Branch(Form("%s_parLoad_mass",it->name),&(it->parLoad_mass),Form("%s_parLoad_mass/F",it->name));
    tr->Branch(Form("%s_parLoad_mass_err",it->name),&(it->parLoad_mass_err),Form("%s_parLoad_mass_err/F",it->name));
    tr->Branch(Form("%s_parIni_ctaumass",it->name),&(it->parIni_ctaumass),Form("%s_parIni_ctaumass/F",it->name));
    tr->Branch(Form("%s_parIni_ctaumass_err",it->name),&(it->parIni_ctaumass_err),Form("%s_parIni_ctaumass_err/F",it->name));
    tr->Branch(Form("%s_parFit_ctaumass",it->name),&(it->parFit_ctaumass),Form("%s_parFit_ctaumass/F",it->name));
    tr->Branch(Form("%s_parFit_ctaumass_err",it->name),&(it->parFit_ctaumass_err),Form("%s_parFit_ctaumass_err/F",it->name));
  }
  
  // list of files
  vector<TString> theFiles = fileList(workDirName,"",DSTag,"",fitType);
  
  cout << ">>>>>> Number of files to be parsed = " << theFiles.size() << endl;
  cout << endl;
  
  int cnt=0;
  for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
    cout << "Parsing file " << cnt << " / " << theFiles.size() << ": " << *it << endl;
    
    // parse the file name to get info
    anabin thebin = binFromFile(*it);
    ptmin = thebin.ptbin().low();
    ptmax = thebin.ptbin().high();
    ymin = thebin.rapbin().low();
    ymax = thebin.rapbin().high();
    strcpy(collSystem, (it->Contains("PbPb")) ? "PbPb" : "PP");
    bool isPP = !(it->Contains("PbPb"));
    if (isPP) thebin.setcentbin(centbin_int);
    centmin = thebin.centbin().low();
    centmax = thebin.centbin().high();
    
    // get the model names
    from = 0;
    bool catchjpsi=false, catchbkg=false, catchtype=false;
    TString modelType;
    while (it->Tokenize(t, from, "_")) {
      if (catchjpsi) {strcpy(jpsiName, t.Data()); catchjpsi=false;}
      if (catchbkg) {strcpy(bkgName, t.Data()); catchbkg=false;}
      if (catchtype) {modelType = t; catchtype=false;}
      if ((t=="Jpsi") || t=="CtauRes") catchjpsi=true;
      if (t=="Bkg") catchbkg=true;
      if (t.EndsWith("FIT")) catchtype=true;
    }
    
    TFile *f = TFile::Open(*it); RooWorkspace *ws = NULL;
    if (!f) {
      cout << "[ERROR] file " << *it << " does not exist." << endl;
    } else {
      ws = (RooWorkspace*) f->Get("workspace");
      if (!ws) {
        cout << "[ERROR] workspace not found in " << *it << "." << endl;
      }
    }
    
    nll=0; /*chi2=0;*/ npar=0; ndof_ctau=0,ndof_invMass=0;
    chi2_invMass=0; //ndof_invMass=0;
    chi2_ctau=0; //ndof_ctau=0;
    if (f && ws) {
      // get the model for nll and npar
      RooAbsPdf *model = (strcmp(modelType.Data(),"CTAURES")) ? pdfFromWS(ws, Form("_%s",collSystem), Form("pdf%s_Tot",modelType.Data())) : pdfFromWS(ws, Form("_%s",collSystem), "pdfCTAUNRES_Tot");
      cout << "The parameters will be read from the PDF " << model->GetName() << endl;
      
      RooAbsPdf *model_bkg = pdfFromWS(ws, Form("_%s",collSystem), Form("pdf%s_Bkg",modelType.Data()));
      
      // get the dataset
      const char* token = (strcmp(DSTag,"DATA") && wantPureSMC) ? Form("_%s_NoBkg",collSystem) : Form("_%s",collSystem);
      RooAbsData *dat = dataFromWS(ws, token, Form("dOS_%s", DSTag));
      
      if (dat) {
        if (model) {
          RooAbsReal *NLL = model->createNLL(*dat);
          if (NLL) nll = NLL->getVal();
          npar = model->getParameters(dat)->selectByAttrib("Constant",kFALSE)->getSize();
          
          RooRealVar *chi2var_invMass = ws->var("chi2_invMass");
          RooRealVar *ndofvar_invMass = ws->var("ndof_invMass");
          RooRealVar *chi2var_ctau = ws->var("chi2_ctau");
          RooRealVar *ndofvar_ctau = ws->var("ndof_ctau");
          if (chi2var_invMass && ndofvar_invMass) {
            chi2_invMass = chi2var_invMass->getVal();
            ndof_invMass = ndofvar_invMass->getVal();
            normchi2_invMass = chi2_invMass/ndof_invMass;
            chi2prob_invMass = TMath::Prob(chi2_invMass,ndof_invMass);
          }
          else {
            cout << "[WARNING]: chi2 or ndof for mass not found" << endl;
          }
          if (chi2var_ctau && ndofvar_ctau) {
            chi2_ctau = chi2var_ctau->getVal();
            ndof_ctau = ndofvar_ctau->getVal();
            normchi2_ctau = chi2_ctau/ndof_ctau;
            chi2prob_ctau = TMath::Prob(chi2_ctau,ndof_ctau);
          }
          else {
            cout << "[WARNING]: chi2 or ndof for ctau not found" << endl;
          }
        }
        else {
          cout << "[ERROR] << the PDF was not found!" << endl;
          return;
        }
        if (model_bkg) {
          nparbkg = model_bkg->getParameters(dat)->selectByAttrib("Constant",kFALSE)->getSize();
        }
      }
      else{
        cout << "[ERROR] DataSet not found!" << endl;
        return;
      }
      
      // get the snapshots
      const RooArgSet *pdfMASS_Tot_parIni = ws->getSnapshot(Form("pdfMASS_Tot_%s_parIni",collSystem));
      const RooArgSet *pdfMASS_Tot_parLoad = ws->getSnapshot(Form("pdfMASS_Tot_%s_parLoad",collSystem));
      const RooArgSet *pdfCTAUMASS_Tot_parIni = ws->getSnapshot(Form("pdfCTAUMASS_Tot_%s_parIni",collSystem));
      const RooArgSet *pdfCTAUMASS_Tot_parFit = ws->getSnapshot(Form("pdfCTAUMASS_Tot_%s_parFit",collSystem));
      
      // get the POIs
      for (vector<poi>::iterator itpoi=thePois.begin(); itpoi!=thePois.end(); itpoi++) {
        RooRealVar *thevar = poiFromWS(ws, Form("_%s",collSystem), itpoi->name);
        RooRealVar *thevar_parIni_mass = pdfMASS_Tot_parIni ? (RooRealVar*) pdfMASS_Tot_parIni->find(Form("%s_%s",itpoi->name,collSystem)) : 0;
        RooRealVar *thevar_parLoad_mass = pdfMASS_Tot_parLoad ? (RooRealVar*) pdfMASS_Tot_parLoad->find(Form("%s_%s",itpoi->name,collSystem)) : 0;
        RooRealVar *thevar_parIni_ctaumass = pdfCTAUMASS_Tot_parIni ? (RooRealVar*) pdfCTAUMASS_Tot_parIni->find(Form("%s_%s",itpoi->name,collSystem)) : 0;
        RooRealVar *thevar_parFit_ctaumass = pdfCTAUMASS_Tot_parFit ? (RooRealVar*) pdfCTAUMASS_Tot_parFit->find(Form("%s_%s",itpoi->name,collSystem)) : 0;
        itpoi->val = thevar ? thevar->getVal() : 0;
        itpoi->errL = thevar ? thevar->getError() : 0;
        itpoi->errH = itpoi->errL;
        itpoi->min = thevar ? thevar->getMin() : 0;
        itpoi->max = thevar ? thevar->getMax() : 0;
        itpoi->parIni_mass = thevar_parIni_mass ? thevar_parIni_mass->getVal() : 0;
        itpoi->parIni_mass_err = thevar_parIni_mass ? thevar_parIni_mass->getError() : 0;
        itpoi->parLoad_mass = thevar_parLoad_mass ? thevar_parLoad_mass->getVal() : 0;
        itpoi->parLoad_mass_err = thevar_parLoad_mass ? thevar_parLoad_mass->getError() : 0;
        itpoi->parIni_ctaumass = thevar_parIni_ctaumass ? thevar_parIni_ctaumass->getVal() : 0;
        itpoi->parIni_ctaumass_err = thevar_parIni_ctaumass ? thevar_parIni_ctaumass->getError() : 0;
        itpoi->parFit_ctaumass = thevar_parFit_ctaumass ? thevar_parFit_ctaumass->getVal() : 0;
        itpoi->parFit_ctaumass_err = thevar_parFit_ctaumass ? thevar_parFit_ctaumass->getError() : 0;
        
        if (TString(itpoi->name).Contains("correl_")) {
          // correlation between two variables
          TString toparse = TString(itpoi->name).ReplaceAll("correl_","");
          TString var1name, var2name;
          from = 0; int i=0;
          while (toparse.Tokenize(t, from, "_vs_")) {
            if (i==0) var1name = t;
            if (i==1) var2name = t;
            i++;
          }
          
          RooFitResult *fr = (RooFitResult*) ws->obj(Form("fitResult_pdf%s_Tot_%s",modelType.Data(),collSystem));
          RooRealVar *thevar1 = poiFromWS(ws, Form("_%s",collSystem), var1name);
          RooRealVar *thevar2 = poiFromWS(ws, Form("_%s",collSystem), var2name);
          if (fr) {
            itpoi->val = fr->correlation(Form("%s_%s",var1name.Data(),collSystem),Form("%s_%s",var2name.Data(),collSystem));
            itpoi->errL = thevar1 && thevar2 ? itpoi->val * thevar1->getError() * thevar2->getError() : 0;
            itpoi->errH = itpoi->errL;
          }
        }
        else if (TString(itpoi->name).Contains("eff") || TString(itpoi->name).Contains("acc")) {
          // efficiency
          TFile *feff(0x0);
          if (TString(itpoi->name).Contains("eff")) feff = TFile::Open(Form("../Efficiency/files_eff/nominal/histos_%s_%s.root",
                                         TString(itpoi->name)=="effnp" ? "npjpsi" : "jpsi",
                                         isPP ? "pp" : "pbpb"));
          else feff = TFile::Open(Form("../Efficiency/files_acc/nominal/histos_%s_%s.root",
                                  TString(itpoi->name)=="accnp" ? "npjpsi" : "jpsi",
                                  isPP ? "pp" : "pbpb"));
          
          int catmin, catmax;
          bool isallcent = (thebin.centbin() == binI(0,200));
          bool isallrap = (thebin.rapbin() == binF(0,2.4));
          bool islowpt = (thebin.ptbin() == binF(3,6.5));
          bool israpFwd = (thebin.rapbin() == binF(1.8,2.4));
          bool ishighpt = ((thebin.ptbin() == binF(6.5,50)) || (thebin.ptbin() == binF(6.5,30)) || (thebin.ptbin() == binF(3.0,30)) );
          bool israpBin = (thebin.rapbin() == binF(0.,0.4)) || (thebin.rapbin() == binF(0.4,0.8)) || (thebin.rapbin() == binF(0.8,1.2)) || (thebin.rapbin() == binF(1.2,1.6)) || (thebin.rapbin() == binF(1.6,2.0)) || (thebin.rapbin() == binF(2.0,2.4));
          bool israpBin4PtCent = ((thebin.rapbin() == binF(0.,0.6)) || (thebin.rapbin() == binF(0.6,1.2)) || (thebin.rapbin() == binF(1.2,1.8)) || (thebin.rapbin() == binF(1.8,2.4)) || (thebin.rapbin() == binF(0.,1.6)) || (thebin.rapbin() == binF(1.6,2.4)));
          
          bool isPtCentPbPb = isallrap && !ishighpt && !islowpt && ((thebin.centbin() == binI(0,20)) || (thebin.centbin() == binI(20,60)) || (thebin.centbin() == binI(60,200)));
          bool isPtCentPP = !isPtCentPbPb && isPP && isallrap && !ishighpt && !islowpt && ((thebin.centbin() == binI(0,200)) && (thebin.ptbin() == binF(15,20) || thebin.ptbin() == binF(20,30) || thebin.ptbin() == binF(30,50)));
          bool isPtRap = !isPtCentPbPb && !isPtCentPP && isallcent && !ishighpt && ((thebin.rapbin() == binF(1.6,2.4) && islowpt) || !islowpt) && israpBin4PtCent;
          bool isPt = !isPtCentPbPb && !isPtCentPP && !isPtRap && isallcent && isallrap && !ishighpt && !islowpt;
          bool isRap = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && isallcent && ishighpt && !isallrap && israpBin;
          bool isCent = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && isallrap && !isallcent && ishighpt;
          bool isCentRap = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && !isCent && ((!isPP && !isallcent) || isPP) && israpBin4PtCent && ishighpt;
          bool isCentRapIntPbPb = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && !isCent && !isCentRap && (!isPP && isallcent) && ishighpt && israpBin4PtCent;
          bool isFwdLowPtCentPbPb = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && !isCent && !isCentRap && !isCentRapIntPbPb && !ishighpt && islowpt && israpFwd && !isPP && !isallcent;
          bool isFwdLowPtCentPP = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && !isCent && !isCentRap && !isCentRapIntPbPb && !isFwdLowPtCentPbPb && !ishighpt && islowpt && israpFwd && isPP && isallcent;
          bool isallInt = !isPtCentPbPb && !isPtCentPP && !isPtRap && !isPt && !isRap && !isCent && !isCentRap && !isFwdLowPtCentPbPb && !isFwdLowPtCentPP && isallcent && isallrap && ishighpt;
          
          TString hnumname, hdenname;
          TString tag1, tag2;
          if (isPtCentPbPb || isPtCentPP || isPtRap || isPt || isCent || isCentRap || isCentRapIntPbPb)
          {
            if (isPtCentPbPb || isPtCentPP) {tag2 = "cent"; tag1 = "pt";}
            if (isPtRap || isPt) {tag2 = "rap"; tag1 = "pt";}
            if (isCent || isCentRap || isCentRapIntPbPb) {tag2 = "rap"; tag1 = "cent";}
            
            if (tag2 == "cent")
            {
              catmin = (TString(itpoi->name).Contains("acc") || isPtCentPP) ? 0 : thebin.centbin().low()/2;
              catmax = (TString(itpoi->name).Contains("acc") || isPtCentPP) ? 10 : thebin.centbin().high()/2;
            }
            else
            {
              catmin = thebin.rapbin().low()*10;
              catmax = thebin.rapbin().high()*10;
            }
            
            hnumname = Form("hnum_%s_%s%02i%02i", tag1.Data(), tag2.Data(), catmin, catmax);
            hdenname = Form("hden_%s_%s%02i%02i", tag1.Data(), tag2.Data(), catmin, catmax);
          }
          else if (isRap || isallInt)
          {
            tag1 = "rap";
            hnumname = "hnum_rap";
            hdenname = "hden_rap";
          }
          else if (isFwdLowPtCentPbPb || isFwdLowPtCentPP)
          {
            tag1 = "cent";
            hnumname = "hnum_cent_rap1824_pt3065";
            hdenname = "hden_cent_rap1824_pt3065";
          }
          else
          {
            cout << "ERROR: no histo defined for this bin:" << endl;
            thebin.print();
          }

          TH1F *hnum = (TH1F*) feff->Get(hnumname);
          TH1F *hden = (TH1F*) feff->Get(hdenname);
          if (!hnum || !hden) {
            thebin.print();
            cout << hnumname << " not found!" << endl;
            itpoi->val = 0;
            itpoi->errL = 0;
            itpoi->errH = itpoi->errL;
            continue;
          }
          
          double numval, numerr, denval, denerr;
          int ibin = hnum->FindBin((thebin.centbin().low()+thebin.centbin().high())/4.);
          if (tag1 == "cent" && (TString(itpoi->name).Contains("acc") || isFwdLowPtCentPP || (isCentRap && isPP))) ibin = 1;
          else if (tag1 == "pt") ibin = hnum->FindBin((thebin.ptbin().low()+thebin.ptbin().high())/2.);
          else if (tag1 == "rap") ibin = hnum->FindBin((thebin.rapbin().low()+thebin.rapbin().high())/2.);
          cout << "bin: " << ibin << " ; " << hnum->GetBinCenter(ibin) << endl;
          numval = hnum->GetBinContent(ibin);
          numerr = hnum->GetBinError(ibin);
          denval = hden->GetBinContent(ibin);
          denerr = hden->GetBinError(ibin);
          if (isallInt || isCentRapIntPbPb) {
            cout << "Will use integrated efficiency" << endl;
            numval = hnum->IntegralAndError(1,hnum->GetNbinsX(),numerr);
            denval = hden->IntegralAndError(1,hden->GetNbinsX(),denerr);
          }
          double efficiency = (denval>0) ? numval / denval : 0;
          itpoi->val = efficiency;
          itpoi->errL = (numval>0 && denval>0) ? efficiency*sqrt(pow(numerr/numval,2)+pow(denerr/denval,2)) : 0;
          itpoi->errH = itpoi->errL;
          
          thebin.print();
          cout << hnumname << " ; " << hdenname << (TString(itpoi->name).Contains("eff") ? " ; eff = " : " ; acc = ") << efficiency << endl;
          
          delete feff;
        }
        else if (TString(itpoi->name)=="lumi") {
          // luminosity and Ncoll
          if (isPP) {
            itpoi->val = lumipp;
            itpoi->errL = 0.023; // from LUM-16-001
            itpoi->errH = itpoi->errL;
          } else {
            if (thebin.centbin().low()>=60) itpoi->val = lumipbpb_peri;
            else itpoi->val = lumipbpb_ABCD;
            itpoi->errL = 0; // FIXME
            itpoi->errH = itpoi->errL;
          }
        }
        else if (TString(itpoi->name)=="ncoll") {
          if (isPP) {
            itpoi->val=1;
            itpoi->errL=0;
            itpoi->errH = itpoi->errL;
          } else {
            itpoi->val = HI::findNcollAverage(thebin.centbin().low(),thebin.centbin().high());
            itpoi->errL = 0; //FIXME
            itpoi->errH = itpoi->errL;
          }
        }
        else if (TString(itpoi->name)=="npart") {
          if (isPP) {
            itpoi->val=2;
            itpoi->errL=0;
            itpoi->errH = itpoi->errL;
          } else {
            itpoi->val = HI::findNpartAverage(thebin.centbin().low(),thebin.centbin().high());
            itpoi->errL = 0; //FIXME
            itpoi->errH = itpoi->errL;
          }
        }
        if (TString(itpoi->name)=="taa") {
          itpoi->val=HI::findTaaAverage(thebin.centbin().low(),thebin.centbin().high());
//          cout << "Taa = " << HI::findTaaAverage(thebin.centbin().low(),thebin.centbin().high()) << endl;
          itpoi->errL=HI::findTaaAverage_err_low(thebin.centbin().low(),thebin.centbin().high());
          itpoi->errH = HI::findTaaAverage_err_high(thebin.centbin().low(),thebin.centbin().high());
        }
        
      }
      
      // delete model;
      // delete model_bkg;
      // delete dat;
      delete ws;
      delete f;
    }
    else {
      for (vector<poi>::iterator itpoi=thePois.begin(); itpoi!=thePois.end(); itpoi++) {
        itpoi->val = 0;
        itpoi->errL = 0;
        itpoi->errH = itpoi->errL;
      }
    }
    
    // fill the tree
    tr->Fill();
    cnt++;
    cout << endl;
  } // loop on the files
  
  f->Write();
  f->Close();
  cout << "Tree with results created" << endl;
  cout << endl;
}

#endif // #ifndef results2tree_C
