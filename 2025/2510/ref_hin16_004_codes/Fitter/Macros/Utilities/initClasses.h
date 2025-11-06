#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "RooBinning.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussModel.h"
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"

#include "../CMS/tdrstyle.C"
#include "../CMS/CMS_lumi.C"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>


using namespace std;
using namespace RooFit;

typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct EvtPar {
  StartEnd RunNb;
  int      TriggerBit;
} EvtPar;

typedef struct DiMuonPar {
  MinMax ctau, ctauN, ctauNRes, ctauRes, ctauErr, ctauTrue, M, Pt, AbsRap;
  string ctauCut;
} DiMuonPar;

typedef struct SiMuonPar {
  MinMax Pt, Eta;
} SiMuonPar;

typedef struct InputOpt {
  int        oniaMode;
  EvtPar     PbPb, pp;
} InputOpt;

typedef struct KinCuts {
  StartEnd   Centrality;
  SiMuonPar  sMuon;
  DiMuonPar  dMuon;
} KinCuts;

bool isEqualKinCuts(struct KinCuts cutA, struct KinCuts cutB, bool isPbPb) 
{
  bool cond = true;

  if (isPbPb) {
    cond = cond && (cutA.Centrality.Start    == cutB.Centrality.Start);
    cond = cond && (cutA.Centrality.End      == cutB.Centrality.End);
  }
  cond = cond && (cutA.sMuon.Pt.Min        == cutB.sMuon.Pt.Min);
  cond = cond && (cutA.sMuon.Pt.Max        == cutB.sMuon.Pt.Max);
  cond = cond && (cutA.sMuon.Eta.Min       == cutB.sMuon.Eta.Min);
  cond = cond && (cutA.sMuon.Eta.Max       == cutB.sMuon.Eta.Max);

  cond = cond && (cutA.dMuon.ctau.Min      == cutB.dMuon.ctau.Min);
  cond = cond && (cutA.dMuon.ctau.Max      == cutB.dMuon.ctau.Max);
  cond = cond && (cutA.dMuon.ctauErr.Min   == cutB.dMuon.ctauErr.Min);
  cond = cond && (cutA.dMuon.ctauErr.Max   == cutB.dMuon.ctauErr.Max);
  cond = cond && (cutA.dMuon.ctauTrue.Min  == cutB.dMuon.ctauTrue.Min);
  cond = cond && (cutA.dMuon.ctauTrue.Max  == cutB.dMuon.ctauTrue.Max);
  cond = cond && (cutA.dMuon.ctauNRes.Min  == cutB.dMuon.ctauNRes.Min);
  cond = cond && (cutA.dMuon.ctauNRes.Max  == cutB.dMuon.ctauNRes.Max);
  cond = cond && (cutA.dMuon.ctauRes.Min   == cutB.dMuon.ctauRes.Min);
  cond = cond && (cutA.dMuon.ctauRes.Max   == cutB.dMuon.ctauRes.Max);
  cond = cond && (cutA.dMuon.ctauN.Min     == cutB.dMuon.ctauN.Min);
  cond = cond && (cutA.dMuon.ctauN.Max     == cutB.dMuon.ctauN.Max);
  cond = cond && (cutA.dMuon.ctauCut       == cutB.dMuon.ctauCut);
  cond = cond && (cutA.dMuon.M.Min         == cutB.dMuon.M.Min);
  cond = cond && (cutA.dMuon.M.Max         == cutB.dMuon.M.Max);
  cond = cond && (cutA.dMuon.Pt.Min        == cutB.dMuon.Pt.Min);
  cond = cond && (cutA.dMuon.Pt.Max        == cutB.dMuon.Pt.Max);
  cond = cond && (cutA.dMuon.AbsRap.Min    == cutB.dMuon.AbsRap.Min);
  cond = cond && (cutA.dMuon.AbsRap.Max    == cutB.dMuon.AbsRap.Max);

  return cond;
}


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z; };
ParticleMass Mass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188};


enum class MassModel 
{
    InvalidModel =0,
    SingleGaussian=1,
    DoubleGaussian=2,
    SingleCrystalBall=3,
    ExtendedCrystalBall=4,
    DoubleCrystalBall=5,
    GaussianAndCrystalBall=6,
    Uniform=7,
    Chebychev1=8,
    Chebychev2=9,
    Chebychev3=10, 
    Chebychev4=11,
    Chebychev5=12,
    Chebychev6=13,
    ExpChebychev1=14,
    ExpChebychev2=15,
    ExpChebychev3=16,
    ExpChebychev4=17,
    ExpChebychev5=18,
    ExpChebychev6=19,
    Exponential=20
};
map< string , MassModel > MassModelDictionary = {
  {"InvalidModel",            MassModel::InvalidModel},
  {"SingleGaussian",          MassModel::SingleGaussian},
  {"DoubleGaussian",          MassModel::DoubleGaussian},
  {"SingleCrystalBall",       MassModel::SingleCrystalBall},
  {"ExtendedCrystalBall",     MassModel::ExtendedCrystalBall},
  {"DoubleCrystalBall",       MassModel::DoubleCrystalBall},
  {"GaussianAndCrystalBall",  MassModel::GaussianAndCrystalBall},
  {"Uniform",                 MassModel::Uniform},
  {"Chebychev1",              MassModel::Chebychev1},
  {"Chebychev2",              MassModel::Chebychev2},
  {"Chebychev3",              MassModel::Chebychev3},
  {"Chebychev4",              MassModel::Chebychev4},
  {"Chebychev5",              MassModel::Chebychev5},
  {"Chebychev6",              MassModel::Chebychev6},
  {"ExpChebychev1",           MassModel::ExpChebychev1},
  {"ExpChebychev2",           MassModel::ExpChebychev2},
  {"ExpChebychev3",           MassModel::ExpChebychev3},
  {"ExpChebychev4",           MassModel::ExpChebychev4},
  {"ExpChebychev5",           MassModel::ExpChebychev5},
  {"ExpChebychev6",           MassModel::ExpChebychev6},
  {"Exponential",             MassModel::Exponential}
};


enum class CtauModel 
{     
    InvalidModel=0,
    QuadrupleGaussianResolution=1,
    TripleGaussianResolution=2,
    DoubleGaussianResolution=3, 
    SingleGaussianResolution=4,
    TripleDecay=5,
    DoubleSingleSidedDecay=6,
    SingleSidedDecay=7,
    Delta=8,
    QuadrupleDecay=9
};
map< string , CtauModel > CtauModelDictionary = {
  {"InvalidModel",             CtauModel::InvalidModel},
  {"QuadrupleGaussianResolution", CtauModel::QuadrupleGaussianResolution},
  {"TripleGaussianResolution", CtauModel::TripleGaussianResolution},
  {"DoubleGaussianResolution", CtauModel::DoubleGaussianResolution},
  {"SingleGaussianResolution", CtauModel::SingleGaussianResolution},
  {"TripleDecay",              CtauModel::TripleDecay},
  {"QuadrupleDecay",           CtauModel::QuadrupleDecay},
  {"DoubleSingleSidedDecay",   CtauModel::DoubleSingleSidedDecay},
  {"SingleSidedDecay",         CtauModel::SingleSidedDecay},
  {"Delta",                    CtauModel::Delta}
};

typedef struct CtauPNP {
  CtauModel    Prompt, NonPrompt;
} CtauPNP;

typedef struct CtauMassModel {
  MassModel  Mass;
  CtauPNP    Ctau;
} FitModel;

typedef struct CharmModel {
  CtauMassModel  Jpsi, Psi2S, Bkg;
  CtauModel CtauRes;
  CtauModel CtauTrue, CtauTrueRes;
} CharmModel;

typedef struct OniaModel {
  CharmModel  PbPb, PP;
} OniaModel;




#include <boost/algorithm/string/replace.hpp>

void setFixedVarsToContantVars(RooWorkspace& ws)
{
  RooArgSet listVar = ws.allVars();
  TIterator* parIt = listVar.createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if ( it->getMin()==it->getMax() && !it->isConstant() ) {
      cout << "[INFO] Setting " << it->GetName() << " constant!" << endl;
      it->setConstant(kTRUE);
    }
  }
};


template<typename T>
T* clone(const T& s)
{
  RooArgSet* cloneSet = (RooArgSet*)RooArgSet(s,s.GetName()).snapshot(kTRUE);
  if (!cloneSet) { cout << "[ERROR] Couldn't deep-clone " << s.GetName() << endl; return NULL; }
  T* obj = (RooAbsPdf*)cloneSet->find(s.GetName());
  if (!obj) { cout << "[ERROR] Couldn't deep-clone " << s.GetName() << endl; delete cloneSet; return NULL; }
  obj->setOperMode(RooAbsArg::ADirty, kTRUE);
  return obj;
};


bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2) {
  TIterator* parIt = pars1->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2->getRealValue(it->GetName(),-1e99);
    if ( strcmp(it->GetName(),"ctauErr")==0 || strcmp(it->GetName(),"ctau")==0 || strcmp(it->GetName(),"ctauTrue")==0 || strcmp(it->GetName(),"ctauRes")==0 || strcmp(it->GetName(),"ctauNRes")==0 || strcmp(it->GetName(),"ctauN")==0 ) continue;
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false;  // the parameter was found, but with a different value!
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMin() != it->getMin() ) return false;  // the parameter has different lower limit
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMax() != it->getMax() ) return false;  // the parameter has different upper limit
    if (string(it->GetName()).find("N_")!=std::string::npos) {
      if ( (((RooRealVar&)(*pars2)[it->GetName()]).getMin() != it->getMin()) && 
           (((RooRealVar&)(*pars2)[it->GetName()]).getMin() != -2000000.0  ) && 
           (((RooRealVar&)(*pars2)[it->GetName()]).getMin() != 0.0         ) ) return false;  // the parameter has different lower limit
      if ( (((RooRealVar&)(*pars2)[it->GetName()]).getMax() != it->getMax()) &&
           (((RooRealVar&)(*pars2)[it->GetName()]).getMax() != 2000000.0   ) ) return false;  // the parameter has different upper limit
    }
  }
  return true;
};


bool setConstant( RooWorkspace& myws, string parName, bool CONST)
{
  if (myws.var(parName.c_str())) { 
    myws.var(parName.c_str())->setConstant(CONST);
    if (CONST) { cout << "[INFO] Setting parameter " << parName << " : " << myws.var(parName.c_str())->getVal() << " to constant value!" << endl; }
  }
  else if (!myws.function(parName.c_str())) { 
    cout << "[ERROR] Parameter " << parName << " was not found!" << endl;
    return false;
  }

  return true;
};


bool saveWorkSpace(RooWorkspace& myws, string outputDir, string FileName)
{
  // Save the workspace
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  cout << FileName << endl;
  TFile *file =  new TFile(FileName.c_str(), "RECREATE");
  if (!file) {
    file->Close(); delete file;
    cout << "[ERROR] Output root file with fit results could not be created!" << endl; return false; 
  } else {
    file->cd();    
    myws.Write("workspace"); 
    file->Write(); file->Close(); delete file;
  }
  return true;
};


bool isCompatibleDataset(const RooDataSet& ds, const RooDataSet& ref, bool checkRange=true)
{
  // Check that the DataSets have the same number of events
  if (ds.numEntries()!=ref.numEntries()){ cout << "[ERROR] DataSets " << ds.GetName() << " and " << ref.GetName() << " don't have the same number of events!" << endl; return false; }
  if (ds.sumEntries()!=ref.sumEntries()){ cout << "[ERROR] DataSets " << ds.GetName() << " and " << ref.GetName() << " don't have the same sum of weights!" << endl; return false; }
  // Check that the input DataSet have the same variables and distributions as the reference
  const RooArgSet* listVar = ref.get();
  TIterator* parIt = listVar->createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    string name = it->GetName();
    if ( ds.get()->find(it->GetName()) == NULL ) { cout << "[ERROR] DataSet " << ds.GetName() << " does not contain the variable " << it->GetName() << " !" << endl; return false; }
    if (checkRange) {
      if ( it->getMin() != ((RooRealVar*)ds.get()->find(it->GetName()))->getMin() ) { cout << "[ERROR] " << it->GetName() << " Min Range disaggrement : "
                                                                                           << ds.GetName() << " ( " << ((RooRealVar*)ds.get()->find(it->GetName()))->getMin() << " ) " << " and "
                                                                                           << ref.GetName() << " ( " << it->getMin() << " ) " << " ! " << endl; return false; }
      if ( it->getMax() != ((RooRealVar*)ds.get()->find(it->GetName()))->getMax() ) { cout << "[ERROR] " << it->GetName() << " Max Range disaggrement : "
                                                                                           << ds.GetName() << " ( " << ((RooRealVar*)ds.get()->find(it->GetName()))->getMax() << " ) " << " and "
                                                                                           << ref.GetName() << " ( " << it->getMax() << " ) " << " ! " << endl; return false; }
    }
    if ( ref.mean(*it) != ds.mean(*it) ) { cout << "[ERROR] " << it->GetName() << " Mean Value disaggrement : "
                                                << ds.GetName() << " ( " << ds.mean(*it) << " ) " << " and "
                                                << ref.GetName() << " ( " << ref.mean(*it) << " ) " << " ! " << endl; return false; }
    if ( ref.moment(*it, 2) != ds.moment(*it, 2) ) { cout << "[ERROR] " << it->GetName() << " Variance Value disaggrement : "
                                                          << ds.GetName() << " ( " << ds.moment(*it, 2) << " ) " << " and "
                                                          << ref.GetName() << " ( " << ref.moment(*it, 2) << " ) " << " ! " << endl; return false; }
    if ( ref.moment(*it, 3) != ds.moment(*it, 3) ) { cout << "[ERROR] " << it->GetName() << " 3rd Moment Value disaggrement : "
                                                          << ds.GetName() << " ( " << ds.moment(*it, 3) << " ) " << " and "
                                                          << ref.GetName() << " ( " << ref.moment(*it, 3) << " ) " << " ! " << endl; return false; }
    if ( ref.moment(*it, 4) != ds.moment(*it, 4) ) { cout << "[ERROR] " << it->GetName() << " 4rd Moment Value disaggrement : "
                                                          << ds.GetName() << " ( " << ds.moment(*it, 4) << " ) " << " and "
                                                          << ref.GetName() << " ( " << ref.moment(*it, 4) << " ) " << " ! " << endl; return false; }
  }
  // DataSets are compatible if they passed all tests
  cout << "[INFO] DataSets " << ds.GetName() << " and " << ref.GetName() << " are compatible!" << endl; return true;
};


bool isFitAlreadyFound(RooArgSet *newpars, string FileName, string pdfName) 
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] FileName: " << FileName << " was not found" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found" << endl;
    file->Close(); delete file;
    return false;
  }
  string snapShotName = Form("%s_parIni", pdfName.c_str());
  const RooArgSet *params = ws->getSnapshot(snapShotName.c_str());
  if (!params) {
    cout << "[INFO] Snapshot parIni was not found" << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }
  bool result = compareSnapshots(newpars, params);
  delete ws;
  file->Close(); delete file; 

  return result;
};


bool loadPreviousFitResult(RooWorkspace& myws, string FileName, string DSTAG, bool isPbPb, bool loadNumberOfEvents, bool updateN)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  cout <<  "[INFO] Loading variables and functions from: " << FileName << endl;
  RooArgSet listVar = ws->allVars();
  TIterator* parIt = listVar.createIterator();
  string print = "[INFO] Variables loaded: ";
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    string name = it->GetName();
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || name=="ctauRes" || name=="ctauNRes" || name=="ctauN" ||
         name=="ctauTrue" || name=="pt" || name=="cent" || 
         name=="rap" || name=="One" ) continue;
    if ( !loadNumberOfEvents && name.find("N_")!=std::string::npos ) continue;
    if (myws.var(name.c_str())) {
      print = print + Form("  %s: %.5f->%.5f ", name.c_str(), myws.var(name.c_str())->getValV(), ws->var(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->var(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( ws->var(name.c_str())->getError() );
    } else {
      if ( (name==Form("lambdaDSS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))) && myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))) ) {
        print = print + Form("  %s: %.5f->%.5f", Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
      if ( (name==Form("sigmaMC_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))) && myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))) {
        print = print + Form("  %s: %.5f->%.5f  ", Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->getValV(), ws->var(name.c_str())->getValV()) ;
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->setVal  ( ws->var(name.c_str())->getValV()  );
        myws.var(Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))->setError( ws->var(name.c_str())->getError() ); 
      }
    }
  }
  cout << print << endl;
  RooArgSet listFun = ws->allFunctions();
  TIterator* parFunIt = listFun.createIterator();
  string printFun = "[INFO] Functions loaded: ";
  for (RooRealVar* it = (RooRealVar*)parFunIt->Next(); it!=NULL; it = (RooRealVar*)parFunIt->Next() ) {
    string name = it->GetName();
    if ( name=="invMass" || name=="ctau" || name=="ctauErr" || name=="ctauRes" || name=="ctauNRes" || name=="ctauN" ||
         name=="ctauTrue" || name=="pt" || name=="cent" || 
         name=="rap" || name=="One" ) continue;
    if ( !loadNumberOfEvents && name.find("N_")!=std::string::npos ) continue;
    if (myws.var(name.c_str())) { 
      printFun = printFun + Form("  %s: %.5f->%.5f  ", name.c_str(), myws.var(name.c_str())->getValV(), ws->function(name.c_str())->getValV()) ;
      myws.var(name.c_str())->setVal  ( ws->function(name.c_str())->getValV()  );
      myws.var(name.c_str())->setError( 0.0 );
    } else {
      Double_t MassRatio = (Mass.Psi2S/Mass.JPsi);
      string reName = name.c_str();
      boost::replace_all(reName, "Jpsi", "Psi2S");
      if (myws.var(reName.c_str())) {
        Double_t value = 0.0;
        if ( (reName==Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) { value = ws->var(name.c_str())->getValV() * MassRatio; }
        else if ( (reName==Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) { value = ws->var(name.c_str())->getValV() * MassRatio; }
        else { value = ws->var(name.c_str())->getValV(); }
        printFun = printFun + Form("  %s: %.5f->%.5f  ", reName.c_str(), myws.var(reName.c_str())->getValV(), value) ;
        myws.var(reName.c_str())->setVal  ( value );
        myws.var(reName.c_str())->setError( 0.0 );
      }
    }
  }
  setFixedVarsToContantVars(myws);
  cout << printFun << endl;

  if (updateN && myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))) {
    string dsName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbPb ?"PbPb":"PP"));
    string dsNameCut = dsName+"_CTAUCUT"; if (!myws.data(dsNameCut.c_str())) dsNameCut = dsName;
    cout << "[INFO] Checking if local dataset " << dsNameCut << " is compatible with Mass fit dataset " << dsName << " !" << endl;
    if (myws.data(dsNameCut.c_str()) && ws->data(dsName.c_str()) && !isCompatibleDataset(*(RooDataSet*)myws.data(dsNameCut.c_str()), *(RooDataSet*)ws->data(dsName.c_str()), false)) {
      // Let's fit again the mass with only the N parameters free, to account for possible ctau or ctauErr cuts in the input datasets
      myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE);
      std::vector< std::string > objs = {"Bkg", "Jpsi", "Psi2S"}; std::map< std::string , double > dN_Old , dN_New, dNErr_Old;
      for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP"))))  setConstant( myws, Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")), false); }
      std::cout << "[INFO] Fitting the mass distribution to update the number of events!" << std::endl;
      for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) {
          dN_Old[obj]    = myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->getVal();
          dNErr_Old[obj] = myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->getError();
        }
      }
      ((RooFitResult*)myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->fitTo(*myws.data(dsNameCut.c_str()), Extended(kTRUE), Range("MassWindow"), NumCPU(32), PrintLevel(-1), Save()))->Print("v");
      for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) { dN_New[obj]    = myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->getVal(); } }
      bool update = true; for (auto obj : objs) {  if (dN_Old.count(obj)>0 && ((dN_New.at(obj) - dN_Old.at(obj))>dNErr_Old.at(obj))) { update = false; } }
      if (update) {
        for (auto obj : objs) {
          if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) cout << "[INFO] Change in " << Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP"))
                                                                                 << " : " << dN_Old.at(obj) << " -> " << dN_New.at(obj) << endl;
        }
      }
      else {
        for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) { myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->setVal(dN_Old.at(obj)); } }
        cout << "[WARNING] Number of Re-Fitted events is larger than in mass fits, setting the N parameter values to its results from mass fits!" << endl;
      }
      // Set the re-fitted N errors to their values from the mass fits
      for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) { myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->setError(dNErr_Old.at(obj)); } }
      myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kFALSE);
    }
  }

  delete ws;
  file->Close(); delete file;
  return true;
};


bool loadCtauErrRange(string FileName, struct KinCuts& cut)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[ERROR] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[ERROR] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (ws->var("ctauErr")) {
    cut.dMuon.ctauErr.Min = ws->var("ctauErr")->getMin();
    cut.dMuon.ctauErr.Max = ws->var("ctauErr")->getMax();
  } else {
    cout << Form("[ERROR] ctauErr was not found!") << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }

  delete ws;
  file->Close(); delete file;
  return true;
};


bool loadYields(RooWorkspace& myws, string FileName, string dsName, string pdfName)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  bool compDS = true;
  if (ws->data(dsName.c_str()) && myws.data(dsName.c_str())) {
    if (isCompatibleDataset(*(RooDataSet*)myws.data(dsName.c_str()), *(RooDataSet*)ws->data(dsName.c_str()))) {
      bool isPbPb = true; if (FileName.find("_PP_")!=std::string::npos) isPbPb = false;
      const RooArgSet *params = ws->getSnapshot(Form("%s_parIni", pdfName.c_str()));
      std::vector< std::string > objs = {"Bkg", "Jpsi", "Psi2S"};
      for (auto obj : objs) {
        string name = Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP"));
        if (params->find(name.c_str()))  {
          string par = Form("%s[ %.1f, %.1f, %.1f ]", name.c_str(), 
                            ((RooRealVar*)params->find(name.c_str()))->getVal(), 
                            ((RooRealVar*)params->find(name.c_str()))->getMin(), 
                            ((RooRealVar*)params->find(name.c_str()))->getMax());
          myws.factory(par.c_str());
          myws.var(name.c_str())->setVal(((RooRealVar*)params->find(name.c_str()))->getVal());
          myws.var(name.c_str())->setMin(((RooRealVar*)params->find(name.c_str()))->getMin());
          myws.var(name.c_str())->setMax(((RooRealVar*)params->find(name.c_str()))->getMax());
          cout << "[INFO] Yield loaded : " << par << endl;
        }
      }
    }
    else { cout << "[INFO] RooDatasets used to extract the Yields are not compatible!" << endl; compDS = false; }
  }
  else { cout << "[INFO] RooDatasets used to extract the Yields were not found!" << endl; compDS = false; }

  delete ws;
  file->Close(); delete file;
  return compDS;
};


bool loadSPlotDS(RooWorkspace& myws, string FileName, string dsName)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[ERROR] File " << FileName << " was not found!" << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[ERROR] Workspace was not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  if (ws->data(dsName.c_str())) {
    myws.import(*ws->data(dsName.c_str()), Rename((dsName+"_INPUT").c_str()));
    if (myws.data((dsName+"_INPUT").c_str())) { cout << "[INFO] RooDataset " << (dsName+"_INPUT") << " was imported!" << endl; }
    else { cout << "[ERROR] Importing RooDataset " << (dsName+"_INPUT") << " failed!" << endl; }
  } else {
    cout << "[ERROR] RooDataset " << dsName << " was not found!" << endl;
    delete ws;
    file->Close(); delete file;
    return false;
  }

  delete ws;
  file->Close(); delete file;
  return true;
};


int importDataset(RooWorkspace& myws, const RooWorkspace& inputWS, struct KinCuts cut, string label, bool cutSideBand=false)
{
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  if (cutSideBand) {
    indMuonMass =  indMuonMass + "&&" + "((2.0 < invMass && invMass < 2.8) || (3.3 < invMass && invMass < 3.5) || (3.9 < invMass && invMass < 5.0))";
  }
  string indMuonRap     = Form("(%.6f <= abs(rap) && abs(rap) < %.6f)",    cut.dMuon.AbsRap.Min,   cut.dMuon.AbsRap.Max);
  string indMuonPt      = Form("(%.6f <= pt && pt < %.6f)",                cut.dMuon.Pt.Min,       cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.6f < ctau && ctau <= %.6f)",            cut.dMuon.ctau.Min,     cut.dMuon.ctau.Max); 
  if(cut.dMuon.ctauCut!=""){ indMuonCtau = cut.dMuon.ctauCut; }
  string indMuonCtauErr = Form("(%.12f < ctauErr && ctauErr < %.12f)",     cut.dMuon.ctauErr.Min,  cut.dMuon.ctauErr.Max);
  string inCentrality   = Form("(%d <= cent && cent < %d)",                cut.Centrality.Start,   cut.Centrality.End);
  string indMuonCtauTrue = Form("(%.12f < ctauTrue && ctauTrue < %.12f)",  cut.dMuon.ctauTrue.Min, cut.dMuon.ctauTrue.Max);
  string indMuonCtauRes = Form("(%.12f < ctauRes && ctauRes < %.12f)",     cut.dMuon.ctauRes.Min,  cut.dMuon.ctauRes.Max);
  string indMuonCtauNRes = Form("(%.12f < ctauNRes && ctauNRes < %.12f)",  cut.dMuon.ctauNRes.Min, cut.dMuon.ctauNRes.Max);
  string indMuonCtauN   = Form("(%.12f < ctauN && ctauN < %.12f)",        cut.dMuon.ctauN.Min, cut.dMuon.ctauN.Max);
  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr;
  if (label.find("PbPb")!=std::string::npos){ strCut = strCut +"&&"+ inCentrality; }
  if (label.find("MC")!=std::string::npos){ strCut = strCut +"&&"+ indMuonCtauTrue +"&&"+ indMuonCtauNRes +"&&"+ indMuonCtauRes;  }
  else { strCut = strCut +"&&"+ indMuonCtauN;  }

  // Reduce and import the datasets
  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return -1;
  }
  cout << "[INFO] Importing local RooDataSet with cuts: " << strCut << endl;
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataOS->sumEntries()==0){ 
    cout << "[ERROR] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return -1;
  }
  myws.import(*dataOS);
  delete dataOS;
  
  if (label.find("NoBkg")==std::string::npos && label.find("AccEff")==std::string::npos && label.find("lJpsiEff")==std::string::npos) // Don't try to find SS dataset if label contais NoBkg or correction
  {
    if (!(inputWS.data(Form("dSS_%s", label.c_str())))){
      cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
      return -1;
    }
    RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
    if (dataSS->sumEntries()==0){
      cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    }
    myws.import(*dataSS);
    delete dataSS;
  }

  // Set the range of each global parameter in the local roodataset
  const RooArgSet* rowOS = myws.data(Form("dOS_%s", label.c_str()))->get();
  ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);        
  ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
  ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);            
  ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
  ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);        
  ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);
  ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
  if (label.find("PbPb")!=std::string::npos){
    ((RooRealVar*)rowOS->find("cent"))->setMin(cut.Centrality.Start);      
    ((RooRealVar*)rowOS->find("cent"))->setMax(cut.Centrality.End);
  }
  if (label.find("MC")!=std::string::npos){
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);      
    ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
    ((RooRealVar*)rowOS->find("ctauRes"))->setMin(cut.dMuon.ctauRes.Min);      
    ((RooRealVar*)rowOS->find("ctauRes"))->setMax(cut.dMuon.ctauRes.Max);
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMin(cut.dMuon.ctauNRes.Min);      
    ((RooRealVar*)rowOS->find("ctauNRes"))->setMax(cut.dMuon.ctauNRes.Max);
  }
  else {
    ((RooRealVar*)rowOS->find("ctauN"))->setMin(cut.dMuon.ctauN.Min);      
    ((RooRealVar*)rowOS->find("ctauN"))->setMax(cut.dMuon.ctauN.Max);
  }

  if (myws.data(Form("dOS_%s_SPLOT_INPUT", label.c_str()))){
    RooDataSet* dataOS = (RooDataSet*)myws.data(Form("dOS_%s_SPLOT_INPUT", label.c_str()))->reduce(strCut.c_str()); 
    if (dataOS) {
      // Set the range of each global parameter in the local roodataset
      const RooArgSet* rowOS = dataOS->get();
      ((RooRealVar*)rowOS->find("invMass"))->setMin(cut.dMuon.M.Min);        
      ((RooRealVar*)rowOS->find("invMass"))->setMax(cut.dMuon.M.Max);
      ((RooRealVar*)rowOS->find("pt"))->setMin(cut.dMuon.Pt.Min);            
      ((RooRealVar*)rowOS->find("pt"))->setMax(cut.dMuon.Pt.Max);
      ((RooRealVar*)rowOS->find("ctau"))->setMin(cut.dMuon.ctau.Min);        
      ((RooRealVar*)rowOS->find("ctau"))->setMax(cut.dMuon.ctau.Max);
      ((RooRealVar*)rowOS->find("ctauErr"))->setMin(cut.dMuon.ctauErr.Min);
      ((RooRealVar*)rowOS->find("ctauErr"))->setMax(cut.dMuon.ctauErr.Max);
      if (label.find("PbPb")!=std::string::npos){
        ((RooRealVar*)rowOS->find("cent"))->setMin(cut.Centrality.Start);      
        ((RooRealVar*)rowOS->find("cent"))->setMax(cut.Centrality.End);
      }
      if (label.find("MC")!=std::string::npos){
        ((RooRealVar*)rowOS->find("ctauTrue"))->setMin(cut.dMuon.ctauTrue.Min);      
        ((RooRealVar*)rowOS->find("ctauTrue"))->setMax(cut.dMuon.ctauTrue.Max);
        ((RooRealVar*)rowOS->find("ctauRes"))->setMin(cut.dMuon.ctauRes.Min);      
        ((RooRealVar*)rowOS->find("ctauRes"))->setMax(cut.dMuon.ctauRes.Max);
        ((RooRealVar*)rowOS->find("ctauNRes"))->setMin(cut.dMuon.ctauNRes.Min);      
        ((RooRealVar*)rowOS->find("ctauNRes"))->setMax(cut.dMuon.ctauNRes.Max);
      }
      else {
        ((RooRealVar*)rowOS->find("ctauN"))->setMin(cut.dMuon.ctauN.Min);      
        ((RooRealVar*)rowOS->find("ctauN"))->setMax(cut.dMuon.ctauN.Max);
      }
      if (dataOS->sumEntries()==0){ cout << "[ERROR] No events from dataset " <<  Form("dOS_%s_SPLOT_INPUT", label.c_str()) << " passed the kinematic cuts!" << endl; }
      else if (!isCompatibleDataset(*dataOS, *(RooDataSet*)myws.data(Form("dOS_%s", label.c_str())))){ cout << "[ERROR] sPlot and Original Datasets are inconsistent!" << endl; delete dataOS; return -1; }
      else {
        myws.import(*dataOS, Rename(Form("dOS_%s_SPLOT", label.c_str())));
        if (myws.data(Form("dOS_%s_SPLOT", label.c_str()))) { cout << "[INFO] RooDataset " << Form("dOS_%s_SPLOT", label.c_str()) << " was imported!" << endl; }
        else { cout << "[ERROR] Importing RooDataset " << Form("dOS_%s_SPLOT", label.c_str()) << " failed!" << endl; delete dataOS; return -1; }
        cout << "[INFO] SPlotDS Events: " << dataOS->sumEntries() << " , origDS Events: " << myws.data(Form("dOS_%s", label.c_str()))->sumEntries() << std::endl;
      }
      delete dataOS;
    }
  }

  // Set the range of each global parameter in the local workspace
  myws.var("invMass")->setMin(cut.dMuon.M.Min);        
  myws.var("invMass")->setMax(cut.dMuon.M.Max);
  myws.var("pt")->setMin(cut.dMuon.Pt.Min);            
  myws.var("pt")->setMax(cut.dMuon.Pt.Max);
  myws.var("rap")->setMin(cut.dMuon.AbsRap.Min);       
  myws.var("rap")->setMax(cut.dMuon.AbsRap.Max);
  myws.var("ctau")->setMin(cut.dMuon.ctau.Min);        
  myws.var("ctau")->setMax(cut.dMuon.ctau.Max);
  myws.var("ctauErr")->setMin(cut.dMuon.ctauErr.Min);  
  myws.var("ctauErr")->setMax(cut.dMuon.ctauErr.Max);
  if (label.find("PbPb")!=std::string::npos){
    myws.var("cent")->setMin(cut.Centrality.Start);      
    myws.var("cent")->setMax(cut.Centrality.End);
  }
  if (label.find("MC")!=std::string::npos){
    myws.var("ctauTrue")->setMin(cut.dMuon.ctauTrue.Min);      
    myws.var("ctauTrue")->setMax(cut.dMuon.ctauTrue.Max);
    myws.var("ctauRes")->setMin(cut.dMuon.ctauRes.Min);      
    myws.var("ctauRes")->setMax(cut.dMuon.ctauRes.Max);
    myws.var("ctauNRes")->setMin(cut.dMuon.ctauNRes.Min);      
    myws.var("ctauNRes")->setMax(cut.dMuon.ctauNRes.Max);
  }
  else {
    myws.var("ctauN")->setMin(cut.dMuon.ctauN.Min);      
    myws.var("ctauN")->setMax(cut.dMuon.ctauN.Max);
  }
  cout << "[INFO] Analyzing bin: " << Form(
                                           "%.3f < pt < %.3f, %.3f < rap < %.3f, %d < cent < %d", 
                                           cut.dMuon.Pt.Min,
                                           cut.dMuon.Pt.Max,
                                           cut.dMuon.AbsRap.Min,
                                           cut.dMuon.AbsRap.Max,
                                           cut.Centrality.Start,
                                           cut.Centrality.End
                                           ) << endl;

  return 1;
};


void printChi2(RooWorkspace& myws, TPad* Pad, RooPlot* frame, string varLabel, string dataLabel, string pdfLabel, int nBins, bool useDefaultName=true)
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1); 
  unsigned int nFitPar = myws.pdf(pdfLabel.c_str())->getParameters(*myws.data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize();
  TH1* hdatact = myws.data(dataLabel.c_str())->createHistogram("hdatact", *myws.var(varLabel.c_str()), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
//  RooHist *hpull = frame->pullHist("hdatact",pdfLabel.c_str(), true);
  RooHist *hpull = frame->pullHist(0,0, true);
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i+1) > 0.0) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  if (useDefaultName) {
    RooRealVar chi2Var("chi2","chi2",chi2);
    RooRealVar ndofVar("ndof","ndof",ndof);
    myws.import(chi2Var); myws.import(ndofVar);
  } else {
    RooRealVar chi2Var((string("chi2_")+varLabel).c_str(),(string("chi2_")+varLabel).c_str(),chi2);
    RooRealVar ndofVar((string("ndof_")+varLabel).c_str(),(string("ndof_")+varLabel).c_str(),ndof);
    myws.import(chi2Var); myws.import(ndofVar);
  }
  delete hdatact; 
  delete hpull;
};


void getRange(TH1* hist, int nMaxBins, vector<double>& rangeErr)
{
  // 1) Find the bin with the maximum Y value
  int binMaximum = hist->GetMaximumBin();
  // 2) Loop backward and find the first bin
  int binWithContent = -1;
  int firstBin = 1;
  for (int i=binMaximum; i>0; i--) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((binWithContent-i) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { firstBin = binWithContent; break; }
      else { binWithContent = i; }
    }
  }
  // 3) Loop forward and find the last bin
  binWithContent = -1;
  int lastBin = hist->GetNbinsX();
  for (int i=binMaximum; i<hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i)>0.0) {
      if ( binWithContent>0 && ((i - binWithContent) > nMaxBins) && hist->GetBinContent(i)<5.0 ) { lastBin = binWithContent+1; break; }
      else { binWithContent = i; }
    }
  }
  // 4) Build the set of bins
  int startBin = ( (firstBin>1) ? (firstBin-1) : firstBin );
  const int nNewBins = lastBin - startBin + 1;
  double binning[nNewBins+1];
  binning[0] = hist->GetXaxis()->GetXmin();
  binning[nNewBins] = hist->GetXaxis()->GetXmax();
  for (int i=1; i<nNewBins; i++) {
    int iBin = startBin + i;
    binning[i] = hist->GetBinLowEdge(iBin); 
  }
  rangeErr.push_back(binning[(firstBin>1)?1:0]);
  rangeErr.push_back(binning[nNewBins-1]);

  return;
};


bool isSPlotDSAlreadyFound(RooWorkspace& myws, string FileName, vector<string> dsNames, bool loadDS)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  bool found = true;
  for (unsigned int i=0; i<dsNames.size(); i++) {
    string dsName = dsNames.at(i);
    if ( !(ws->data(dsName.c_str())) ) {
      cout << "[INFO] " << dsName << " was not found in: " << FileName << endl; found = false;
    }
    if (loadDS && found) {
      myws.import(*(ws->data(dsName.c_str())));
      if (myws.data(dsName.c_str())) { cout << "[INFO] sPlot DataSet " << dsName << " succesfully imported!" << endl; }
      else {  cout << "[ERROR] sPlot DataSet " << dsName << " import failed!" << endl; found = false; }
    }
  }
  
  delete ws;
  file->Close(); delete file;

  return found;
};

bool isPdfAlreadyFound(RooWorkspace& myws, string FileName, vector<string> pdfNames, bool loadCtauErrPdf)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }

  bool found = true;
  for (unsigned int i=0; i<pdfNames.size(); i++) {
    string pdfName = pdfNames.at(i);
    string dataName = pdfName;
    dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
    if ( !(ws->pdf(pdfName.c_str())) ) {
      cout << "[INFO] " << pdfName << " was not found in: " << FileName << endl; found = false;
    }
    if (loadCtauErrPdf && found) {
      myws.import(*(ws->pdf(pdfName.c_str())));
      if (myws.pdf(pdfName.c_str()))   { cout << "[INFO] Pdf " << pdfName << " succesfully imported!" << endl;       }
      else { cout << "[ERROR] Pdf " << pdfName << " import failed!" << endl; found = false; }
      if (ws->data(dataName.c_str())) { myws.import(*(ws->data(dataName.c_str()))); }
      if (myws.data(dataName.c_str())) { cout << "[INFO] DataHist " << dataName << " succesfully imported!" << endl; }
    }
  }
  
  delete ws;
  file->Close(); delete file;

  return found;
};


#endif // #ifndef initClasses_h
