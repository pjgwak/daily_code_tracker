// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program creates two RooDataSets (opposite-sign and same-sign dimuons) from an Onia Tree.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 26 19:08 CET 2016
//
//

#include "Utilities/initOniaTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"

#include "TObjArray.h"

map<int, double>   fCentMap; // map for centrality-Ncoll mapping
double             fCentBinning[200];
int                fCentBins;
TObjArray*         fcorrArray = NULL; // Array with the 2D correction for weighting

string  findMyTree(string FileName);
bool    getTChain(TChain* fChain, vector<string> FileNames);
void    iniBranch(TChain* fChain,bool isMC=false);
bool    checkDS(RooDataSet* DS, string DSName);
double  deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
bool    isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR=0.03);
double  getNColl(int centr, bool isPP);
double  getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP);
bool    readCorrection(const char* file);
void    setCentralityMap(const char* file);

bool tree2DataSet(RooWorkspace& Workspace, vector<string> InputFileNames, string DSName, string OutputFileName, bool UpdateDS=false)
{
  RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL; RooDataSet* dataOSNoBkg = NULL;
  
  bool isMC = false;
  if (DSName.find("MC")!=std::string::npos) isMC =true;
  
  bool isPbPb = false;
  if (DSName.find("PbPb")!=std::string::npos) isPbPb =true;
  int triggerIndex_PP   = PP::HLT_HIL1DoubleMu0_v1;
  int triggerIndex_PbPb = HI::HLT_HIL1DoubleMu0_v1;
  int CentFactor = 1;
  
  bool usePeriPD = false;
  if (InputFileNames[0].find("HIOniaPeripheral30100")!=std::string::npos) {
    cout << "[INFO] Working with Peripheral PbPb PD" << endl;
    usePeriPD = true;
    triggerIndex_PbPb = HI::HLT_HIL1DoubleMu0_2HF_Cent30100_v1;
  }
  
  bool applyWeight = false;
  if (isMC && isPbPb) applyWeight = true;
  
  bool isPureSDataset = false;
  if (OutputFileName.find("_PureS")!=std::string::npos) isPureSDataset = true;

  bool applyWeight_Corr = false;
  if ( (OutputFileName.find("_AccEff")!=std::string::npos) || (OutputFileName.find("_lJpsiEff")!=std::string::npos) ) applyWeight_Corr = true;
  if(applyWeight == true) applyWeight_Corr = false;
  
  TString corrName = "";
  TString corrFileName = "";
  if (OutputFileName.find("_AccEff")!=std::string::npos)
  {
    corrFileName = "correction_AccEff.root";
    corrName = "AccEff";
  }
  else if (OutputFileName.find("_lJpsiEff")!=std::string::npos)
  {
    corrFileName = "correction_lJpsiEff.root";
    corrName = "lJpsiEff";
  }
  
  bool createDS = ( gSystem->AccessPathName(OutputFileName.c_str()) || UpdateDS );
  if ( !gSystem->AccessPathName(OutputFileName.c_str()) ) {
    cout << "[INFO] Loading RooDataSet from " << OutputFileName << endl;
    
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"READ");
    if (isMC && isPureSDataset) {
      dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg", DSName.c_str()));
      if (checkDS(dataOSNoBkg, DSName)==false) { createDS = true; }
    } 
    else if (applyWeight_Corr) {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
    }
    else {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s", DSName.c_str()));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
      dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s", DSName.c_str()));
      if (checkDS(dataSS, DSName)==false) { createDS = true; }
    }
    DBFile->Close(); delete DBFile;
  }

  if (createDS) {
    cout << "[INFO] Creating " << (isPureSDataset ? "pure signal " : "") << "RooDataSet for " << DSName << endl;
    TreeName = findMyTree(InputFileNames[0]); if(TreeName==""){return false;}
    
    TChain* theTree = new TChain(TreeName.c_str(),"");
    if(!getTChain(theTree, InputFileNames)){ return false; }     // Import files to TChain
    initOniaTree(theTree);                                       // Initialize the Onia Tree
    iniBranch(theTree,isMC);                                     // Initialize the Branches
    
    RooRealVar* mass    = new RooRealVar("invMass","#mu#mu mass", 1.0, 6.0, "GeV/c^{2}");
    RooRealVar* ctau    = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauN    = new RooRealVar("ctauN","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauTrue = new RooRealVar("ctauTrue","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauNRes = new RooRealVar("ctauNRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauRes = new RooRealVar("ctauRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauErr = new RooRealVar("ctauErr","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ptQQ    = new RooRealVar("pt","#mu#mu p_{T}", -1.0, 10000.0, "GeV/c");
    RooRealVar* rapQQ   = new RooRealVar("rap","#mu#mu y", -2.5, 2.5, "");
    RooRealVar* cent    = new RooRealVar("cent","centrality", -1.0, 1000.0, "");
    RooRealVar* weight  = new RooRealVar("weight","MC weight", 0.0, 10000000.0, "");
    RooRealVar* weightCorr   = new RooRealVar("weightCorr","Data correction weight", 0.0, 10000000.0, "");
    RooArgSet*  cols    = NULL;
    
    if (applyWeight)
    {
      setCentralityMap(Form("%s/Input/CentralityMap_PbPb2015.txt",gSystem->ExpandPathName(gSystem->pwd())));
      if (isMC) {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weight);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
      } else {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weight);
        cols->add(*ctauN);
      }
      dataOS = new RooDataSet(Form("dOS_%s", DSName.c_str()), "dOS", *cols, WeightVar(*weight), StoreAsymError(*mass));
      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, WeightVar(*weight), StoreAsymError(*mass));
      if (isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg", DSName.c_str()), "dOSNoBkg", *cols, WeightVar(*weight), StoreAsymError(*mass));
    }
    else if (applyWeight_Corr)
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weightCorr);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
      } else {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weightCorr);
        cols->add(*ctauN);
      }
      if (!readCorrection(Form("%s/Input/%s",gSystem->ExpandPathName(gSystem->pwd()),corrFileName.Data()))){ return false; }
      dataOS = new RooDataSet(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()), "dOS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      //      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      cout<<"--- 1./applyWeight_Corr applied---"<<endl;
    }
    else
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
      } else {
        cols = new RooArgSet(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent);
        cols->add(*ctauN);
      }  
      dataOS = new RooDataSet(Form("dOS_%s", DSName.c_str()), "dOS", *cols, StoreAsymError(*mass));
      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, StoreAsymError(*mass));
      if (isMC && isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg", DSName.c_str()), "dOSNoBkg", *cols, StoreAsymError(*mass));
    }
    
    Long64_t nentries = theTree->GetEntries();
    //nentries = 50000;
    
    float normF = 0.;
    if (isMC && isPbPb)
    {
      cout << "[INFO] Computing sum of weights for " << nentries << " nentries" << endl;
      
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries << endl;
        
        if (theTree->LoadTree(jentry)<0) break;
        if (theTree->GetTreeNumber()!=fCurrent) {
          fCurrent = theTree->GetTreeNumber();
          cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
        }
        
        theTree->GetEntry(jentry);
       // normF += theTree->GetWeight()*getNColl(Centrality,!isPbPb);
        normF += theTree->GetWeight();  
    }
      
      normF = nentries/normF;
    }
    
    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries << endl;
      
      if (theTree->LoadTree(jentry)<0) break;
      if (theTree->GetTreeNumber()!=fCurrent) {
        fCurrent = theTree->GetTreeNumber();
        cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
      }
      Reco_QQ_4mom->Clear();
      Reco_QQ_mumi_4mom->Clear();
      Reco_QQ_mupl_4mom->Clear();
      if (isMC) {
        Gen_QQ_mumi_4mom->Clear();
        Gen_QQ_mupl_4mom->Clear();
      }
      theTree->GetEntry(jentry);
      
      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
        TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
        mass->setVal(RecoQQ4mom->M());
        if (theTree->GetBranch("Reco_QQ_ctau3D")) { ctau->setVal(Reco_QQ_ctau3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctau")) { ctau->setVal(Reco_QQ_ctau[iQQ]); }
        else { cout << "[ERROR] No ctau information found in the Onia Tree" << endl; }
        if (theTree->GetBranch("Reco_QQ_ctauErr3D")) { ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctauErr")) { ctauErr->setVal(Reco_QQ_ctauErr[iQQ]); }
        else { cout << "[ERROR] No ctauErr information found in the Onia Tree" << endl; }
        
        ctauN->setVal(ctau->getVal()/ctauErr->getVal());
        
        ptQQ->setVal(RecoQQ4mom->Pt());
        rapQQ->setVal(RecoQQ4mom->Rapidity());
        cent->setVal(Centrality*CentFactor);
        if (isMC) {
          if (theTree->GetBranch("Reco_QQ_ctauTrue3D")) { ctauTrue->setVal(Reco_QQ_ctauTrue3D[iQQ]); }
          else if (theTree->GetBranch("Reco_QQ_ctauTrue")) { ctauTrue->setVal(Reco_QQ_ctauTrue[iQQ]); }
          else { cout << "[ERROR] No ctauTrue information found in the Onia Tree" << endl; }
          ctauNRes->setVal( (ctau->getValV() - ctauTrue->getValV())/(ctauErr->getValV()) );
          ctauRes->setVal( (ctau->getValV() - ctauTrue->getValV()) );
        }

        if (applyWeight){
          double w = theTree->GetWeight();
          if (isMC && isPbPb) w = w*normF;//*getNColl(Centrality,!isPbPb)*normF;
          weight->setVal(w);
        }
        else if (applyWeight_Corr)
        {
          double Corr = getCorr(RecoQQ4mom->Rapidity(),RecoQQ4mom->Pt(),RecoQQ4mom->M(),!isPbPb);
          double wCorr = 1/Corr;
          weightCorr->setVal(wCorr);
        }
        
        if (
            ( RecoQQ::areMuonsInAcceptance2015(iQQ) ) &&  // 2015 Global Muon Acceptance Cuts
            ( RecoQQ::passQualityCuts2015(iQQ)      ) &&  // 2015 Soft Global Muon Quality Cuts
            ( isPbPb ? (RecoQQ::isTriggerMatch(iQQ,triggerIndex_PbPb) || (usePeriPD ? RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF0_Cent30100_v1) : (RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF_v1) || RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF0_v1)))) :
              RecoQQ::isTriggerMatch(iQQ, triggerIndex_PP) )     // if PbPb && !periPD then (HLT_HIL1DoubleMu0_v1 || HLT_HIL1DoubleMu0_2HF_v1)
            )
        {
          if (Reco_QQ_sign[iQQ]==0) { // Opposite-Sign dimuons
            if (isMC && isPureSDataset && isMatchedRecoDiMuon(iQQ)) dataOSNoBkg->add(*cols, (applyWeight ? weight->getVal() : 1.0)); // Signal-only dimuons
            else if (applyWeight_Corr) dataOS->add(*cols,weightCorr->getVal()); //Signal and background dimuons
            else dataOS->add(*cols, ( applyWeight ? weight->getVal() : 1.0)); //Signal and background dimuons            
          }
          else { // Like-Sign dimuons
            if (!isPureSDataset && !applyWeight_Corr ) dataSS->add(*cols, ( applyWeight  ? weight->getVal() : 1.0));
          }
        }
      }
    }
    // Close the TChain and all its pointers
    delete Reco_QQ_4mom; delete Reco_QQ_mumi_4mom; delete Reco_QQ_mupl_4mom; delete Gen_QQ_mumi_4mom; delete Gen_QQ_mupl_4mom;
    theTree->Reset(); delete theTree;
    
    // Save all the datasets
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"RECREATE");
    DBFile->cd();
    if (isMC && isPureSDataset) dataOSNoBkg->Write(Form("dOS_%s_NoBkg", DSName.c_str()));
    else if (applyWeight_Corr) dataOS->Write(Form("dOS_%s_%s", DSName.c_str(),corrName.Data()));
    else
    {
      dataOS->Write(Form("dOS_%s", DSName.c_str()));
      dataSS->Write(Form("dSS_%s", DSName.c_str()));
    }
    DBFile->Write(); DBFile->Close(); delete DBFile;
  }
  
  // Import datasets to workspace
  if (isMC && isPureSDataset)
  {
    if (!dataOSNoBkg) { cout << "[ERROR] " << DSName << "_NoBkg was not found" << endl; return false; }
    Workspace.import(*dataOSNoBkg);
  }
  else if (applyWeight_Corr)
  {
    if(!dataOS) { cout << "[ERROR] " << DSName << "_" << corrName.Data() << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
  }
  else
  {
    if(!dataOS || !dataSS) { cout << "[ERROR] " << DSName << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
    Workspace.import(*dataSS);
  }
  
  // delete the local datasets
  delete dataSS; delete dataOS; delete dataOSNoBkg;
  
  // delete the correction array
  if (fcorrArray) delete fcorrArray;
  
  return true;
};

string findMyTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  string name = "";
  if(f->GetListOfKeys()->Contains("hionia")){ name = "hionia/myTree"; }
  else if(f->GetListOfKeys()->Contains("myTree")){ name = "myTree"; }
  else { cout << "[ERROR] myTree was not found in: " << FileName << endl; }
  f->Close(); delete f;
  return name;
};

bool getTChain(TChain *fChain, vector<string> FileNames)
{
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;
  for (vector<string>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
    cout << "[INFO] Adding TFile " << FileName->c_str() << endl;
    fChain->Add(Form("%s/%s", FileName->c_str(),  TreeName.c_str()));
  }
  if (!fChain) { cout << "[ERROR] fChain was not created, some input files are missing" << endl; return false; }
  return true;
};

void iniBranch(TChain* fChain, bool isMC)
{
  cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(false); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(false); }
  if (isMC) {
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->GetBranch("Gen_QQ_mupl_4mom")->SetAutoDelete(false); }
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->GetBranch("Gen_QQ_mumi_4mom")->SetAutoDelete(false); }
  }
  fChain->SetBranchStatus("*",0);
  RecoQQ::iniBranches(fChain);
  if (fChain->GetBranch("Centrality"))        { fChain->SetBranchStatus("Centrality",1);        }
  if (fChain->GetBranch("Reco_QQ_size"))      { fChain->SetBranchStatus("Reco_QQ_size",1);      }
  if (fChain->GetBranch("Reco_QQ_sign"))      { fChain->SetBranchStatus("Reco_QQ_sign",1);      }
  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->SetBranchStatus("Reco_QQ_4mom",1);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_ctau3D"))    { fChain->SetBranchStatus("Reco_QQ_ctau3D",1);    }
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) { fChain->SetBranchStatus("Reco_QQ_ctauErr3D",1); }
  if (fChain->GetBranch("Reco_QQ_ctau"))      { fChain->SetBranchStatus("Reco_QQ_ctau",1);      }
  if (fChain->GetBranch("Reco_QQ_ctauErr"))   { fChain->SetBranchStatus("Reco_QQ_ctauErr",1);   }
  if (isMC)
  {
    if (fChain->GetBranch("Gen_QQ_size"))      { fChain->SetBranchStatus("Gen_QQ_size",1);      }
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1); }
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue3D")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue3D",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue",1); }
  }
};

bool checkDS(RooDataSet* DS, string DSName)
{
  bool incCent     = (DSName.find("PbPb")!=std::string::npos);
  bool incCtauTrue = (DSName.find("MC")!=std::string::npos);
  const RooArgSet* row = DS->get();
  if (
      (row->find("invMass")!=0) &&
      (row->find("pt")!=0)      &&
      (row->find("ctau")!=0)    &&
      (row->find("ctauErr")!=0) &&
      (incCent     ? row->find("cent")!=0     : true) &&
      (incCtauTrue ? row->find("ctauTrue")!=0 : true) &&
      (incCtauTrue ? row->find("ctauRes")!=0 : true) &&
      (incCtauTrue ? row->find("ctauNRes")!=0 : true) &&
      (incCtauTrue ? true : row->find("ctauN")!=0)
      ) 
    { return true; }
  else 
    { cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << endl; }

  return false;
};

double deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

bool isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
{
  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);
  
  bool isMatched(false);
  int iGenMuon(0);
  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
  {
    TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
    TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) isMatched = true;
    iGenMuon++;
  }
  
  return isMatched;
};

double getNColl(int centr, bool isPP)
{
  // Returns the corresponding Ncoll value to the "centr" centrality bin
  
  if ( isPP ) return 1.;
  
  int normCent = TMath::Nint(centr/2.);
  
  int lcent = 0;
  int ucent = 0;
  for ( int i = 0 ; i < fCentBins ; i++ )
  {
    ucent = fCentBinning[i];
    if ( (normCent >= lcent) && (normCent < ucent) ) return fCentMap[ucent];
    else lcent = ucent;
  }
  return 1.;
};

void setCentralityMap(const char* file)
{
  // Creates a mapping between centrality and Ncoll, based on a text file (taken from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHeavyIonCentrality)
  
  if ( strlen(file) > 0 )
  {
    char line[1024];
    ifstream in(file);
    float lcent;
    float ucent;
    float Ncoll;
    
    fCentBins = 0;
    while ( in.getline(line,1024,'\n'))
    {
      sscanf(line,"%f %f %f",&lcent,&ucent,&Ncoll);
      
      fCentMap[ucent] = Ncoll;
      fCentBinning[fCentBins++] = ucent;
    }
    if ( fCentBins == 0 ) std::cout << "[INFO] No centrality map could be defined: The file provided is empty" << std::endl;
    else std::cout << "[INFO] Defining centrality map" << std::endl;
  }
  else
  {
    fCentBins = 0;
    std::cout << "[INFO] No centrality map could be defined: No file provided" << std::endl;
  }
};

bool readCorrection(const char* file)
{
  TFile *froot = new TFile(file,"READ");
  if (!froot)
  {
    cout << "[ERROR] File "<< file << " for correction of events not found" << endl;
    return false;
  }
  
  TList* lcorr = froot->GetListOfKeys();
  TIter nextCorr(lcorr);
  
  fcorrArray = new TObjArray();
  fcorrArray->SetOwner(kTRUE);
  
  TObjString* fname(0x0);
  while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
  {
    TH2* h = static_cast<TH2*>(froot->FindObjectAny(fname->GetString().Data()));

    TString sName(h->GetName());
    if ( sName.Contains("hcorr") ) fcorrArray->Add(h->Clone());
    else cout << "[WARNING] Correction histo " << sName.Data() << " not according to naming convention. Not included in correction array" << endl;
  }
  
  if (!(fcorrArray->GetSize()>0))
  {
    cout << "[ERROR] Correction array empty: No corrections found." << endl;
    return false;
  }
  delete lcorr;
//  froot->Close(); delete froot;
  
  return true;
};

double getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP)
{
  const char* collName = "PbPb";
  const char* massName = "Interp";
  if (isPP) collName = "PP";
  if (mass>3.5) massName = "Psi2S";
  else if (mass<3.3) massName = "Jpsi";
  
  if (!fcorrArray)
  {
    cout << "[ERROR] No correction array exist" << endl;
    return 0;
  }

  Double_t corr = 1.;
  if (!strcmp(massName,"Interp"))
  {
    TH2* corrHistoJpsi = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s",collName)));
    TH2* corrHistoPsi2S = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_Psi2S_%s",collName)));
    if (!corrHistoJpsi || !corrHistoPsi2S)
    {
      std::cout << "[Error] No histogram provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
      return 1.;
    }
    
    Int_t binJpsi = corrHistoJpsi->FindBin(fabs(rapidity), pt);
    Double_t corrJpsi = corrHistoJpsi->GetBinContent(binJpsi);
    
    Int_t binPsi2S = corrHistoPsi2S->FindBin(fabs(rapidity), pt);
    Double_t corrPsi2S = corrHistoPsi2S->GetBinContent(binPsi2S);
    
    corr = ((corrPsi2S - corrJpsi)/(3.5-3.3))*(mass-3.3) + corrJpsi;
    
  }
  else
  {
    TH2* corrHisto = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_%s_%s",massName,collName)));
    if (!corrHisto)
    {
      std::cout << "[Error] No histogram provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
      return 1.;
    }
    
    Int_t bin = corrHisto->FindBin(fabs(rapidity), pt);
    corr = corrHisto->GetBinContent(bin);
  }
  if(corr<0.00001) corr=1.0;
  
  return corr;
};
