// For PbPb2023 data

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
// #include "cutsAndBins.h"
using std::cout; using std::string;

static const long MAXTREESIZE = 1000000000000;

// 2018
// bool IsAcceptanceQQ(double pt, double eta)
// {
//   return ((fabs(eta) < 0.3 && pt > 3.4) ||
//           (0.3 < fabs(eta) && fabs(eta) < 1.1 && pt > 3.3) ||
//           (1.1 < fabs(eta) && fabs(eta) < 1.5 && pt > 9.08 - 5.25*fabs(eta)) ||
//           (1.5 < fabs(eta) && fabs(eta) < 2.4 && pt > 0.8 && pt>2.4-0.8*fabs(eta)));
// }

// Run3
bool IsAcceptanceQQ(double pt, double eta)
{
  const double aeta = std::abs(eta);

  if (aeta > 2.4)
    return false;
  if (aeta <= 1.0)
    return pt >= 3.3;
  if (aeta <= 1.3)
    return pt >= ((2.1 - 3.3) / (1.3 - 1.0)) * (aeta - 1.0) + 3.3;
  if (aeta <= 1.7)
    return pt >= ((1.0 - 2.1) / (1.7 - 1.3)) * (aeta - 1.3) + 2.1;

  return pt >= 1.0;
}

double GetPtWeight(TF1 *fW, double pt, double y)
{
  (void)y;
  if (!fW)
    return 1.0;
  return fW->Eval(pt);
}


void efficiency_1d(long nEvt = 1000, bool isPr = true, bool isPtWeight = false, bool isMC = true)
{
  cout << "Start onia_to_skim_data()\n";

  // ===== read input =====
  TChain *oniaTree = new TChain("hionia/myTree");
  
  // input is always MC
  std::string inputFile = isPr ? "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcPr.root" 
                          : "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcNp.root";
  oniaTree->Add(inputFile.c_str());


  // ===== read pT reweight input =====
  TF1 *fWpt = nullptr;
  if (isPtWeight)
  {
    // Placeholder switch: set true when pt-weight TF1 is prepared in a ROOT file.
    const bool usePtWeightFromRoot = false;
    const std::string ptWeightFilePath = "/data/Oniatree/light_ions_Raa/OO_PrivateMc/pt_weight_input.root"; // TODO: set real file name
    const std::string ptWeightFuncName = "fWpt2"; // TODO: set real object name
    if (usePtWeightFromRoot)
    {
      TFile *fPtWeightIn = TFile::Open(ptWeightFilePath.c_str(), "READ");
      if (fPtWeightIn && !fPtWeightIn->IsZombie())
      {
        TF1 *fWptIn = dynamic_cast<TF1 *>(fPtWeightIn->Get(ptWeightFuncName.c_str()));
        if (fWptIn)
        {
          fWpt = (TF1 *)fWptIn->Clone("fWpt_loaded");
        }
        else
        {
          cout << "[WARN] pT weight TF1 not found: " << ptWeightFuncName << "\n";
        }
      }
      else
      {
        cout << "[WARN] failed to open pT weight file: " << ptWeightFilePath << "\n";
      }
      if (fPtWeightIn)
      {
        fPtWeightIn->Close();
        delete fPtWeightIn;
      }
    }

    // Fallback/default: unit weight (tune parameters after fitting if needed).
    if (!fWpt)
    {
      fWpt = new TF1("fWpt", "exp([0]+[1]*x+[2]*x*x)", 0, 100);
      fWpt->SetParameters(0.0, 0.0, 0.0); // default weight = 1
    }
  }

  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 0, massHigh = 200;

  // ===== labeling =====
  std::string mcLabel = "";
  if (isMC)
    mcLabel = isPr ? "PR" : "NP";

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 1000;

  std::vector<float> *Reco_QQ_4mom_m = nullptr;
  std::vector<float> *Reco_QQ_4mom_pt = nullptr;
  std::vector<float> *Reco_QQ_4mom_phi = nullptr;
  std::vector<float> *Reco_QQ_4mom_eta = nullptr;
  
  std::vector<float> *Reco_mu_4mom_pt = nullptr;
  std::vector<float> *Reco_mu_4mom_phi = nullptr;
  std::vector<float> *Reco_mu_4mom_eta = nullptr;
  std::vector<float> *Gen_QQ_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_pt = nullptr;
  std::vector<float> *Gen_QQ_4mom_phi = nullptr;
  std::vector<float> *Gen_QQ_4mom_eta = nullptr;
  std::vector<float> *Gen_mu_4mom_pt = nullptr;
  std::vector<float> *Gen_mu_4mom_phi = nullptr;
  std::vector<float> *Gen_mu_4mom_eta = nullptr;

  oniaTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
  oniaTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);

  oniaTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
  oniaTree->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
  oniaTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
  oniaTree->SetBranchAddress("Gen_QQ_4mom_m", &Gen_QQ_4mom_m);
  oniaTree->SetBranchAddress("Gen_QQ_4mom_pt", &Gen_QQ_4mom_pt);
  oniaTree->SetBranchAddress("Gen_QQ_4mom_phi", &Gen_QQ_4mom_phi);
  oniaTree->SetBranchAddress("Gen_QQ_4mom_eta", &Gen_QQ_4mom_eta);
  oniaTree->SetBranchAddress("Gen_mu_4mom_pt", &Gen_mu_4mom_pt);
  oniaTree->SetBranchAddress("Gen_mu_4mom_phi", &Gen_mu_4mom_phi);
  oniaTree->SetBranchAddress("Gen_mu_4mom_eta", &Gen_mu_4mom_eta);
  

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
  Short_t Gen_QQ_size;
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
  Short_t Gen_QQ_mupl_idx[maxBranchSize];
  Short_t Gen_QQ_mumi_idx[maxBranchSize];

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
  Short_t Reco_mu_whichGen[maxBranchSize]; // MC
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
  oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
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
  oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
  oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

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
    // oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
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
      Form("skim_roots/eff_OO2025_isMC%d%s_ptW%d_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root",
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), isPtWeight),
      "recreate");

  const double ptBinsMid[] = {7, 8, 9, 10, 12, 14, 16, 20, 35};
  const double ptBinsFwd[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20};
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;

  TH1D *hist_eff_den_mid = new TH1D("hist_eff_den_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_den_fwd = new TH1D("hist_eff_den_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_num_mid = new TH1D("hist_eff_num_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_num_fwd = new TH1D("hist_eff_num_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  

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
  long count_gen_dimuon = 0;
  long count_reco_dimuon = 0;

  // ---- preaparation before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if ((iev % 100000 == 0) && (iev != 0))
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "GEN accepted dimuon: " << count_gen_dimuon << "\n";
      cout << "Reco selected dimuon: " << count_reco_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total GEN accepted dimuon: " << count_gen_dimuon << "\n";
      cout << "Total Reco selected dimuon: " << count_reco_dimuon << "\n";
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
    if (Gen_QQ_size < 0) continue;
    if (Reco_QQ_size<0) continue;

    // apply HLT trigger
    int kTrigSel = 2; // HLT_OxyL1SingleMuOpen_v1
    if (!((HLTriggers & ((int)std::pow(2, kTrigSel))) == ((int)std::pow(2, kTrigSel))))
    {
      continue;
    }

    // if (TMath::Abs(vz) > 15) continue; // not included in the Jun Chen's slide

    // ===== start GEN dimuon loop =====
    // reset number of GEN-dimuon
    int nGenDimu = 0;

    for (Int_t irqq = 0; irqq < Gen_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID)
        cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nGenDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nGenDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      // check TLorentzVector pointers
      const int iMuPl = Gen_QQ_mupl_idx[irqq], iMuMi = Gen_QQ_mumi_idx[irqq]; // index guard

      // user-defined kinematic cuts
      float gen_mass = Gen_QQ_4mom_m->at(irqq);
      float gen_eta = Gen_QQ_4mom_eta->at(irqq);
      float gen_eta1 = Gen_mu_4mom_eta->at(iMuPl);
      float gen_eta2 = Gen_mu_4mom_eta->at(iMuMi);
      float gen_pt = Gen_QQ_4mom_pt->at(irqq);
      float gen_pt1 = Gen_mu_4mom_pt->at(iMuPl);
      float gen_pt2 = Gen_mu_4mom_pt->at(iMuMi);
      float gen_phi = Gen_QQ_4mom_phi->at(irqq);

      // reconstrcut TLorentz vector to get the rapidityu
      TLorentzVector JP;
      JP.SetPtEtaPhiM(gen_pt, gen_eta, gen_phi, gen_mass);
      float gen_y = JP.Rapidity();

      if (gen_mass < massLow || gen_mass > massHigh)
        continue;

      // ===== Efficiency denominator: same as acceptance numerator =====
      if (!(TMath::Abs(gen_y) < 2.4) ||
          !IsAcceptanceQQ(gen_pt1, gen_eta1) ||
          !IsAcceptanceQQ(gen_pt2, gen_eta2))
        continue;

      if (std::abs(gen_y) < 1.6)
      {
        const double wpt = isPtWeight ? GetPtWeight(fWpt, gen_pt, gen_y) : 1.0;
        const int ptBin = hist_eff_den_mid->FindBin(gen_pt);
        if (ptBin > 0 && ptBin <= hist_eff_den_mid->GetNbinsX())
          hist_eff_den_mid->Fill(gen_pt, weight * wpt);
      }
      else
      {
        const double wpt = isPtWeight ? GetPtWeight(fWpt, gen_pt, gen_y) : 1.0;
        const int ptBin = hist_eff_den_fwd->FindBin(gen_pt);
        if (ptBin > 0 && ptBin <= hist_eff_den_fwd->GetNbinsX())
          hist_eff_den_fwd->Fill(gen_pt, weight * wpt);
      }

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      ++nGenDimu;
    } // end of dimuon loop

    if (nGenDimu > 0)
    {
      count_gen_dimuon += nGenDimu;
      // flowTree->Fill();
      if (VERBOSE_EVENT)
        cout << "  -> GEN accepted (" << nGenDimu << " dimuons)\n";
    }

    // ===== Reco dimuon loop =====
    // reset number of dimuon
    nDimu = 0;

    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      // check TLorentzVector pointers
      const int iMuPl = Reco_QQ_mupl_idx[irqq], iMuMi = Reco_QQ_mumi_idx[irqq]; // index guard

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
      
      if(mass_ < massLow || mass_ > massHigh) continue;

      // check MC Reco-Gen match 
      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      
      // isSoftCutBased
      if (!Reco_mu_isSoftCutBased[iMuPl] || !Reco_mu_isSoftCutBased[iMuMi]) continue;

      // dimuon y < 2.4 and single muon acceptance cut - 2018 soft
      if ( !(TMath::Abs(y_) < 2.4) ||
           !IsAcceptanceQQ(pt1_, eta1_) ||
           !IsAcceptanceQQ(pt2_, eta2_)) continue;

      // isGlobal - all muons in the test file pass isTracker
      if (!Reco_mu_isGlobal[iMuPl] || !Reco_mu_isGlobal[iMuMi]) continue;

      // vertex probability cut
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue; // not included in JunChen's slide

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      if (Reco_QQ_sign[irqq] != 0) continue; // basic

      // recoQQsign[nDimu] = Reco_QQ_sign[irqq];

      // ===== Efficiency numerator =====
      if (std::abs(y_) < 1.6)
      {
        const double wpt = isPtWeight ? GetPtWeight(fWpt, pt_, y_) : 1.0;
        const int ptBin = hist_eff_num_mid->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_eff_num_mid->GetNbinsX())
          hist_eff_num_mid->Fill(pt_, weight * wpt);
      }
      else
      {
        const double wpt = isPtWeight ? GetPtWeight(fWpt, pt_, y_) : 1.0;
        const int ptBin = hist_eff_num_fwd->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_eff_num_fwd->GetNbinsX())
          hist_eff_num_fwd->Fill(pt_, weight * wpt);
      }

      ++nDimu;
    } // end of dimuon loop

    if (nDimu > 0)
    {
      count_reco_dimuon += nDimu;
      // flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== compute efficiency =====
  hist_eff_den_mid->Scale(1.0, "width");
  hist_eff_den_fwd->Scale(1.0, "width");
  hist_eff_num_mid->Scale(1.0, "width");
  hist_eff_num_fwd->Scale(1.0, "width");

  TH1D *hist_eff_mid = (TH1D *)hist_eff_num_mid->Clone("hist_eff_mid");
  hist_eff_mid->SetTitle(";p_{T} (GeV/c);Efficiency");
  hist_eff_mid->Divide(hist_eff_num_mid, hist_eff_den_mid, 1.0, 1.0, "B");

  TH1D *hist_eff_fwd = (TH1D *)hist_eff_num_fwd->Clone("hist_eff_fwd");
  hist_eff_fwd->SetTitle(";p_{T} (GeV/c);Efficiency");
  hist_eff_fwd->Divide(hist_eff_num_fwd, hist_eff_den_fwd, 1.0, 1.0, "B");

  // ===== save results =====
  fFlowSkim->cd();
  hist_eff_den_mid->Write();
  hist_eff_den_fwd->Write();
  hist_eff_num_mid->Write();
  hist_eff_num_fwd->Write();
  hist_eff_mid->Write();
  hist_eff_fwd->Write();
  // flowTree->Write("myTree");
  fFlowSkim->Close();

  // ----- print -----
  cout << "Total GEN accepted dimuon: " << count_gen_dimuon << "\n";
  cout << "Total Reco selected dimuon: " << count_reco_dimuon << "\n";
  // cout << "Total dimuon (GEN+Reco): " << (count_gen_dimuon + count_reco_dimuon) << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
