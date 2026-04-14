// For PbPb2023 data

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TVector3.h>
#include <TLorentzVector.h>
// #include "cutsAndBins.h"
using std::cout; using std::string;

static const long MAXTREESIZE = 1000000000000;

// Run3 - Man
// bool IsAcceptanceQQ(double pt, double eta)
// {
//   const double aeta = std::abs(eta);

//   if (aeta > 2.4)
//     return false;
//   if (aeta <= 1.0)
//     return pt >= 3.3;
//   if (aeta <= 1.3)
//     return pt >= ((2.1 - 3.3) / (1.3 - 1.0)) * (aeta - 1.0) + 3.3;
//   if (aeta <= 1.7)
//     return pt >= ((1.0 - 2.1) / (1.7 - 1.3)) * (aeta - 1.3) + 2.1;

//   return pt >= 1.0;
// }

// Run3 - JunChen
double PtThreshold(double eta)
{
  const double x = eta;

  const int nPts = 10;
  const double mEta[nPts] = {-2.4, -1.7, -1.3, -1.3, -1.0,
                             1.0, 1.3, 1.3, 1.7, 2.4};

  const double mPt[nPts] = {1.0, 1.0, 1.53, 2.1, 3.3,
                            3.3, 2.1, 1.53, 1.0, 1.0};

  int iseg = -1;

  for (int i = 0; i < nPts - 1; i++)
  {
    if (x >= mEta[i] && x < mEta[i + 1])
    {
      iseg = i;
      break;
    }
  }

  if (x == mEta[nPts - 1])
    iseg = nPts - 2;

  if (iseg < 0)
    return 999999.;

  double mSlope = (mPt[iseg + 1] - mPt[iseg]) /
                  (mEta[iseg + 1] - mEta[iseg]);

  double mPtTh = mSlope * (x - mEta[iseg]) + mPt[iseg];

  return mPtTh;
}
bool IsAcceptanceQQ(double pt, double eta)
{
  double th = PtThreshold(eta);
  return pt >= th;
}

void onia_to_skim_vector(long nEvt = -1, bool isMC = false, bool isPr = false, bool useGlobal = true)
{
  cout << "Start onia_to_skim_data()\n";
  TStopwatch timer;
  timer.Start();

  // ===== read input =====
  TChain *oniaTree = new TChain("hionia/myTree");
  std::string inputFile = "/data/Oniatree/light_ions_Raa/OO_Data/OniaTree_OODimuon_MINIAOD_Run2025OO_PromptReco_v1_Jul12_merged.root";
  if (isMC)
  {
    inputFile = isPr
        ? "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcPr.root"
        : "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcNp.root";
  }
  oniaTree->Add(inputFile.c_str());

  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 2.6, massHigh = 3.5;

  // ===== labeling =====
  std::string mcLabel = "";
  if (isMC)
    mcLabel = isPr ? "PR" : "NP";
  const std::string globalLabel = useGlobal ? "globalOn" : "globalOff";

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 1000;

  std::vector<float> *Reco_QQ_4mom_m = nullptr;
  std::vector<float> *Reco_QQ_4mom_pt = nullptr;
  std::vector<float> *Reco_QQ_4mom_phi = nullptr;
  std::vector<float> *Reco_QQ_4mom_eta = nullptr;
  
  std::vector<float> *Reco_mu_4mom_pt = nullptr;
  std::vector<float> *Reco_mu_4mom_phi = nullptr;
  std::vector<float> *Reco_mu_4mom_eta = nullptr;

  oniaTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);

  oniaTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
  oniaTree->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
  oniaTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
  

  // event-level scalars
  UInt_t runNb = 0;
  UInt_t eventNb = 0;
  UInt_t LS = 0;
  float zVtx = 0.f;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Float_t SumET_HF;
  Float_t Gen_weight; // MC

  // // collection sizes
  Short_t Reco_QQ_size;
  // Short_t Reco_mu_size;

  // // TClonesArray pointer
  // TClonesArray *Reco_QQ_4mom = nullptr;
  // TClonesArray *Reco_mu_4mom = nullptr;

  // // per-dimuon (size = Reco_QQ_size)
  ULong64_t Reco_QQ_trig[maxBranchSize];
  Float_t Reco_QQ_VtxProb[maxBranchSize];
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_sign[maxBranchSize];
  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctauErr3D[maxBranchSize];
  Float_t Reco_QQ_ctau[maxBranchSize];
  Float_t Reco_QQ_ctauErr[maxBranchSize];

  // // per-muon (size = Reco_mu_size)
  ULong64_t Reco_mu_trig[maxBranchSize];
  Bool_t Reco_mu_highPurity[maxBranchSize];
  Int_t Reco_mu_nTrkHits[maxBranchSize];
  // Float_t Reco_mu_normChi2_global[maxBranchSize];
  Int_t Reco_mu_nMuValHits[maxBranchSize];
  // Int_t Reco_mu_StationsMatched[maxBranchSize];
  // Float_t Reco_mu_dxy[maxBranchSize];
  // Float_t Reco_mu_dxyErr[maxBranchSize];
  Float_t Reco_mu_dz[maxBranchSize];
  Float_t Reco_mu_dzErr[maxBranchSize];
  Int_t Reco_mu_nTrkWMea[maxBranchSize];
  Int_t Reco_mu_nPixWMea[maxBranchSize];
  Int_t Reco_mu_nPixValHits[maxBranchSize];
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Int_t Reco_mu_whichGen[maxBranchSize]; // MC
  Bool_t Reco_mu_isSoftCutBased[maxBranchSize];
  Bool_t Reco_mu_isGlobal[maxBranchSize];
  Float_t Reco_QQ_dca[maxBranchSize];

  // // ----- SetBranchAddress -----
  if (oniaTree->GetBranch("runNb"))
    oniaTree->SetBranchAddress("runNb", &runNb);
  if (oniaTree->GetBranch("eventNb"))
    oniaTree->SetBranchAddress("eventNb", &eventNb);
  if (oniaTree->GetBranch("LS"))
    oniaTree->SetBranchAddress("LS", &LS);
  if (oniaTree->GetBranch("zVtx"))
    oniaTree->SetBranchAddress("zVtx", &zVtx);
  // oniaTree->SetBranchAddress("Centrality", &Centrality);
  oniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
  // oniaTree->SetBranchAddress("SumET_HF", &SumET_HF);

  // // sizes
  oniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  // oniaTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

  // TLorentzVecotr
  // oniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  // oniaTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // // dimuon
  oniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  oniaTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  oniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  oniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  oniaTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  oniaTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr);

  // single-muon
  oniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  oniaTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  oniaTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  // oniaTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  oniaTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  // oniaTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  // oniaTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  // oniaTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  // oniaTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  // oniaTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  oniaTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  oniaTree->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased);
  oniaTree->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal);
  oniaTree->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca);

  // MC only
  if (isMC)
  {
    if (oniaTree->GetBranch("Gen_weight"))
      oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
    if (oniaTree->GetBranch("Reco_mu_whichGen"))
      oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }

  // ----- declare local variables ----
  const static long long int nMaxDimu = 1000;
  
  // event-level scalar
  Int_t event = 0;
  Int_t runN = 0;
  Int_t lumi = 0;
  Float_t vz = 0;

  // dimuon size
  Int_t nDimu = 0;

  // dimuon candidate-level variables (size = nDimu)
  Float_t mass[nMaxDimu];
  Float_t y[nMaxDimu];
  Float_t pt[nMaxDimu], pt1[nMaxDimu], pt2[nMaxDimu];
  Float_t eta[nMaxDimu], eta1[nMaxDimu], eta2[nMaxDimu];
  Float_t phi[nMaxDimu], phi1[nMaxDimu], phi2[nMaxDimu];
  Int_t recoQQsign[nMaxDimu], cBin[nMaxDimu];
  Float_t ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];
  Float_t ctau[nMaxDimu], ctauErr[nMaxDimu], ctauRes[nMaxDimu];
  Float_t recoQQdca[nMaxDimu];

  // weights
  Double_t weight = 1.0; // event-level weight
  Double_t TnPweight[nMaxDimu]; // per-candidate level weight

  TFile *fFlowSkim = new TFile(
      Form("skim_roots/skim_OO2025_isMC%d%s_%s_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root",
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), globalLabel.c_str()),
      "recreate");

  cout << "useGlobal: " << useGlobal << "\n";
  cout << "Output file label: " << globalLabel << "\n";

  // fFlowSkim->SetCompressionSettings(207); // LZMA: 207, ZLIB: 1xx
  TTree *flowTree = new TTree("myTree", "");
  flowTree->SetMaxTreeSize(MAXTREESIZE);
  // flowTree->SetAutoSave(50 * 1024 * 1024); // save every 50 MB

  // ----- SetBranchAddress -----
  // event values
  flowTree->Branch("eventNb", &event, "eventNb/I");
  flowTree->Branch("runNb", &runN, "runNb/I");
  flowTree->Branch("LS", &lumi, "LS/I");
  flowTree->Branch("zVtx", &vz, "zVtx/F");

  // dimuon size
  flowTree->Branch("nDimu", &nDimu, "nDimu/I");

  // per-candidate (size = [nDimu])
  flowTree->Branch("mass", mass, "mass[nDimu]/F");
  flowTree->Branch("y", y, "y[nDimu]/F");
  flowTree->Branch("cBin", cBin, "cBin[nDimu]/I");
  flowTree->Branch("pt", pt, "pt[nDimu]/F");
  flowTree->Branch("pt1", pt1, "pt1[nDimu]/F");
  flowTree->Branch("pt2", pt2, "pt2[nDimu]/F");
  flowTree->Branch("eta", eta, "eta[nDimu]/F");
  flowTree->Branch("eta1", eta1, "eta1[nDimu]/F");
  flowTree->Branch("eta2", eta2, "eta2[nDimu]/F");
  flowTree->Branch("phi", phi, "phi[nDimu]/F");
  flowTree->Branch("phi1", phi1, "phi1[nDimu]/F");
  flowTree->Branch("phi2", phi2, "phi2[nDimu]/F");
  flowTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  flowTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  flowTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  flowTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
  flowTree->Branch("ctau", ctau, "ctau[nDimu]/F");
  flowTree->Branch("ctauErr", ctauErr, "ctauErr[nDimu]/F");
  flowTree->Branch("ctauRes", ctauRes, "ctauRes[nDimu]/F");
  flowTree->Branch("recoQQdca", recoQQdca, "recoQQdca[nDimu]/F");

  // weights
  flowTree->Branch("weight", &weight, "weight/D");
  flowTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");

  // ===== event loop =====
  // counters
  long count_dimuon = 0;
  long recoCandTotal = 0;
  long recoPassMass = 0;
  long recoPassGenMatch = 0;
  long recoPassSoftMuon = 0;
  long recoPassMuonTrig = 0;
  long recoPassAcceptance = 0;
  long recoPassGlobal = 0;
  long recoPassFiniteCtau = 0;
  long recoPassSign = 0;
  long count_eta_mid = 0;
  long count_eta_fwd = 0;

  // ---- preaparation before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);
  const Long64_t reportEvery = std::clamp<Long64_t>(nEvt / 100, 1000, 200000);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level
  const ULong64_t trigMask = (1ULL << 2);

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    if ((iev % reportEvery) == 0 || (iev + 1) == nEvt)
    {
      const double elapsed = timer.RealTime();
      timer.Continue();
      const double frac = (nEvt > 0) ? double(iev + 1) / double(nEvt) : 0.0;
      const double estTotal = (frac > 0.0) ? (elapsed / frac) : 0.0;
      const double etaSec = estTotal - elapsed;

      auto toHms = [](double sec) {
        int h = int(sec / 3600);
        sec -= 3600 * h;
        int m = int(sec / 60);
        sec -= 60 * m;
        int s = int(sec + 0.5);
        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << h << ":"
           << std::setw(2) << m << ":" << std::setw(2) << s;
        return os.str();
      };

      cout << "[" << std::fixed << std::setprecision(1) << (frac * 100.0) << "%] "
           << (iev + 1) << " / " << nEvt
           << "  | elapsed " << toHms(elapsed)
           << "  | ETA " << toHms(std::max(0.0, etaSec)) << "\n";
      cout << "Saved dimuon: " << count_dimuon << "\n";
      cout << "Selected dimuon by absRecoEta: mid=" << count_eta_mid
           << ", fwd=" << count_eta_fwd << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total saved dimuon: " << count_dimuon << "\n";
      cout << "Selected dimuon by absRecoEta: mid=" << count_eta_mid
           << ", fwd=" << count_eta_fwd << "\n";
    }

    // read event in Oniatree
    oniaTree->GetEntry(iev);

    // cout << "is Reco_QQ_4mom_pt->size() == Reco_QQ_size: " << (n == qq_size) << endl;
    
    // event-level values
    runN = runNb;
    event = eventNb;
    lumi = LS;
    vz = zVtx;

    // centrality
    int cBin_ = -999; // no centrality
    // if(isMC) cBin_ = Centrality;
    // else if(!isMC) {
    //   cBin_ = getHiBinFromhiHF(SumET_HF);
    //   // if(hiHFBinEdge ==0) cBin_ = getHiBinFromhiHF(SumET_HF);
    //   // else if(hiHFBinEdge == 1) cBin_ = getHiBinFromhiHF_Up(SumET_HF);
    //   // else if(hiHFBinEdge == -1) cBin_ = getHiBinFromhiHF_Down(SumET_HF);
    // } 
    // if(cBin_==-999){ cout << "ERROR!!! No HF Centrality Matching!!" << endl; return;}
    // if(cBin_ < cLow || cBin_ > cHigh) continue;
    
    // NColl - MC only, Galuber model
    double Ncoll_weight = 0, Gen_weight_ = 0;
    // if(isMC){
    //   weight = findNcoll(Centrality) * Gen_weight;
    //   Ncoll_weight = findNcoll(Centrality);
    //   Gen_weight_ = Gen_weight;
    // }

    // check dimoun number
    if (Reco_QQ_size<0) continue;

    // apply HLT trigger
    if ((HLTriggers & trigMask) == 0)
    {
      continue;
    }

    // if (TMath::Abs(vz) > 15) continue; // not included in the Jun Chen's slide
    
    // ----- dimuon loop -----
    // reset number of dimuon
    nDimu = 0;

    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      ++recoCandTotal;

      if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      const int iMuPl = Reco_QQ_mupl_idx[irqq], iMuMi = Reco_QQ_mumi_idx[irqq]; // index guard
      if (iMuPl < 0 || iMuMi < 0 || iMuPl >= maxBranchSize || iMuMi >= maxBranchSize) continue;

      if (!Reco_QQ_4mom_m || !Reco_QQ_4mom_pt || !Reco_QQ_4mom_phi || !Reco_QQ_4mom_eta ||
          !Reco_mu_4mom_pt || !Reco_mu_4mom_phi || !Reco_mu_4mom_eta) continue;

      if (irqq >= static_cast<Int_t>(Reco_QQ_4mom_m->size()) ||
          irqq >= static_cast<Int_t>(Reco_QQ_4mom_pt->size()) ||
          irqq >= static_cast<Int_t>(Reco_QQ_4mom_phi->size()) ||
          irqq >= static_cast<Int_t>(Reco_QQ_4mom_eta->size()) ||
          iMuPl >= static_cast<Int_t>(Reco_mu_4mom_pt->size()) ||
          iMuPl >= static_cast<Int_t>(Reco_mu_4mom_phi->size()) ||
          iMuPl >= static_cast<Int_t>(Reco_mu_4mom_eta->size()) ||
          iMuMi >= static_cast<Int_t>(Reco_mu_4mom_pt->size()) ||
          iMuMi >= static_cast<Int_t>(Reco_mu_4mom_phi->size()) ||
          iMuMi >= static_cast<Int_t>(Reco_mu_4mom_eta->size())) continue;

      // user-defined kinematic cuts
      float mass_ = Reco_QQ_4mom_m->at(irqq);

      float eta_ = Reco_QQ_4mom_eta->at(irqq);
      float eta1_ = Reco_mu_4mom_eta->at(iMuPl);
      float eta2_ = Reco_mu_4mom_eta->at(iMuMi);
      
      float pt_ = Reco_QQ_4mom_pt->at(irqq);
      float pt1_ = Reco_mu_4mom_pt->at(iMuPl);
      float pt2_ = Reco_mu_4mom_pt->at(iMuMi);
      
      float phi_ = Reco_QQ_4mom_phi->at(irqq);
      float phi1_ = Reco_mu_4mom_phi->at(iMuPl);
      float phi2_ = Reco_mu_4mom_phi->at(iMuMi);


      // reconstrcut TLorentz vector to get the rapidityu
      TLorentzVector JP;
      JP.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
      float y_ = JP.Rapidity();
      const float absEta_ = std::fabs(eta_);
      
      if(mass_ < massLow || mass_ > massHigh) continue;
      ++recoPassMass;

      // check MC Reco-Gen match 
      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      ++recoPassGenMatch;
      
      // isSoftCutBased
      if (!Reco_mu_isSoftCutBased[iMuPl] || !Reco_mu_isSoftCutBased[iMuMi]) continue;
      ++recoPassSoftMuon;

      if (((Reco_mu_trig[iMuPl] & trigMask) == 0) &&
          ((Reco_mu_trig[iMuMi] & trigMask) == 0)) continue;
      ++recoPassMuonTrig;

      // dimuon y < 2.4 and single muon acceptance cut
      if ( !(TMath::Abs(y_) < 2.4) ||
           !IsAcceptanceQQ(pt1_, eta1_) ||
           !IsAcceptanceQQ(pt2_, eta2_)) continue;
      ++recoPassAcceptance;

      if (useGlobal)
      {
        if (!Reco_mu_isGlobal[iMuPl] || !Reco_mu_isGlobal[iMuMi]) continue;
      }
      ++recoPassGlobal;

      // vertex probability cut
      // if (Reco_QQ_VtxProb[irqq] < 0.01f) continue; // not included in JunChen's slide

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      if (!std::isfinite(Reco_QQ_ctau3D[irqq]) || !std::isfinite(Reco_QQ_ctauErr3D[irqq])) continue;
      ++recoPassFiniteCtau;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];
      if (Reco_QQ_sign[irqq] != 0) continue; // opposite sign (basic)
      ++recoPassSign;

      // ----- fill per-candidate outputs -----
      if (isMC) TnPweight[nDimu] = tnp_weight;
      else TnPweight[nDimu] = 1.0;

      mass[nDimu] = mass_;
      
      y[nDimu] = y_;

      cBin[nDimu] = cBin_;
      
      pt[nDimu] = pt_;
      pt1[nDimu] = pt1_;
      pt2[nDimu] = pt2_;

      eta[nDimu] = eta_;
      eta1[nDimu] = eta1_;
      eta2[nDimu] = eta2_;

      phi[nDimu] = phi_;
      phi1[nDimu] = phi1_;
      phi2[nDimu] = phi2_;

      recoQQsign[nDimu] = Reco_QQ_sign[irqq];
      recoQQdca[nDimu] = Reco_QQ_dca[irqq];

      // ctau3D
      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      ctau3DRes[nDimu] = (std::fabs(ctau3DErr[nDimu]) > 0.f)
                             ? (ctau3D[nDimu] / ctau3DErr[nDimu])
                             : 0.f;

      // ctau (2D)
      ctau[nDimu] = Reco_QQ_ctau[irqq];
      ctauErr[nDimu] = Reco_QQ_ctauErr[irqq];
      ctauRes[nDimu] = (std::fabs(ctauErr[nDimu]) > 0.f)
                             ? (ctau[nDimu] / ctauErr[nDimu])
                             : 0.f;

      if (absEta_ < 1.6f) ++count_eta_mid;
      else if (absEta_ < 2.4f) ++count_eta_fwd;

      ++nDimu;
    } // end of dimuon loop

    if (nDimu > 0)
    {
      count_dimuon += nDimu;
      flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== save results =====
  fFlowSkim->cd();
  flowTree->Write("myTree");
  fFlowSkim->Close();

  // ----- print -----
  cout << "Total saved dimuon: " << count_dimuon << "\n";
  cout << "Selected dimuon by absRecoEta: mid=" << count_eta_mid
       << ", fwd=" << count_eta_fwd << "\n";
  cout << "Reco candidate cutflow: total=" << recoCandTotal
       << ", mass=" << recoPassMass
       << ", genMatch=" << recoPassGenMatch
       << ", softMuon=" << recoPassSoftMuon
       << ", muTrig=" << recoPassMuonTrig
       << ", acceptance=" << recoPassAcceptance
       << ", global=" << recoPassGlobal
       << ", finiteCtau=" << recoPassFiniteCtau
       << ", sign=" << recoPassSign << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
