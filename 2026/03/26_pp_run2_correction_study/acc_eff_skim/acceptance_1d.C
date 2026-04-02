// For PbPb2023 data

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
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

// Run3 OO MuAccCut - from Jun Chen
bool IsAcceptanceQQJunChenOO(double pt, double eta)
{
  const double aeta = std::abs(eta);
  const double p = pt * std::cosh(eta);

  return (pt >= 0.5 && p >= 2.5 && aeta <= 2.5);
}

double GetPtWeight(TF1 *fW, double pt, double y)
{
  (void)y;
  if (!fW)
    return 1.0;
  return fW->Eval(pt);
}

void acceptance_1d(long nEvt = 10000, bool isPr = true,
                   bool isNCollW = false, bool isGenW = true, bool isPtW = false)
{
  cout << "Start onia_to_skim_data()\n";

  // ===== read input =====
  TChain *oniaTree = new TChain("hionia/myTree");
  std::string inputFile = isPr ? "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcPr.root"
                          : "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcNp.root";
  const int nFilesAdded = oniaTree->Add(inputFile.c_str());

  TF1 *fPtW_mid = nullptr;
  TF1 *fPtW_fwd = nullptr;
  if (isPtW)
  {
    const std::string ptWeightFilePath =
        "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/pT_reweight/fit_ratio_outputs/roots/fit_ratio_exp2.root";
    const std::string ptWeightMidName = isPr ? "exp2_pr_mid" : "exp2_np_mid";
    const std::string ptWeightFwdName = isPr ? "exp2_pr_fwd" : "exp2_np_fwd";

    TFile *fPtWeightIn = TFile::Open(ptWeightFilePath.c_str(), "READ");
    if (!fPtWeightIn || fPtWeightIn->IsZombie())
    {
      cout << "[ERROR] failed to open pT weight file: " << ptWeightFilePath << "\n";
      return;
    }

    TF1 *fPtWMidIn = dynamic_cast<TF1 *>(fPtWeightIn->Get(ptWeightMidName.c_str()));
    TF1 *fPtWFwdIn = dynamic_cast<TF1 *>(fPtWeightIn->Get(ptWeightFwdName.c_str()));
    if (!fPtWMidIn || !fPtWFwdIn)
    {
      cout << "[ERROR] failed to load pT weight functions: "
           << ptWeightMidName << ", " << ptWeightFwdName << "\n";
      fPtWeightIn->Close();
      delete fPtWeightIn;
      return;
    }

    fPtW_mid = static_cast<TF1 *>(fPtWMidIn->Clone(Form("%s_clone", ptWeightMidName.c_str())));
    fPtW_fwd = static_cast<TF1 *>(fPtWFwdIn->Clone(Form("%s_clone", ptWeightFwdName.c_str())));
    fPtWeightIn->Close();
    delete fPtWeightIn;
  }


  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 2.6, massHigh = 3.5;

  // ===== labeling =====
  const std::string mcLabel = isPr ? "PR" : "NP";
  const std::string weightLabel = Form("_ncollW%d_genW%d_ptW%d",
                                       isNCollW, isGenW, isPtW);

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 5000;

  std::vector<float> *Gen_QQ_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_pt = nullptr;
  std::vector<float> *Gen_QQ_4mom_phi = nullptr;
  std::vector<float> *Gen_QQ_4mom_eta = nullptr;
  
  std::vector<float> *Gen_mu_4mom_pt = nullptr;
  std::vector<float> *Gen_mu_4mom_phi = nullptr;
  std::vector<float> *Gen_mu_4mom_eta = nullptr;

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
  Float_t Gen_weight = 1.f; // MC

  // // collection sizes
  Short_t Gen_QQ_size;
  // Short_t Gen_mu_size;

  // // TClonesArray pointer
  // TClonesArray *Gen_QQ_4mom = nullptr;
  // TClonesArray *Gen_mu_4mom = nullptr;

  // // per-dimuon (size = Gen_QQ_size)
  Short_t Gen_QQ_mupl_idx[maxBranchSize];
  Short_t Gen_QQ_mumi_idx[maxBranchSize];
  Short_t Gen_mu_charge[maxBranchSize];
  Float_t Gen_QQ_ctau3D[maxBranchSize];
  Float_t Gen_QQ_ctau[maxBranchSize];

  // // per-muon (size = Gen_mu_size)
  // Float_t Gen_mu_normChi2_global[maxBranchSize];
  // Int_t Gen_mu_StationsMatched[maxBranchSize];
  // Float_t Gen_mu_dxy[maxBranchSize];
  // Float_t Gen_mu_dxyErr[maxBranchSize];
  Float_t Gen_mu_dz[maxBranchSize];
  Float_t Gen_mu_dzErr[maxBranchSize];
  Int_t Gen_mu_SelectionType[maxBranchSize];
  Int_t Gen_mu_whichGen[maxBranchSize]; // MC

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
  // oniaTree->SetBranchAddress("Gen_mu_size", &Gen_mu_size);

  // TLorentzVecotr
  // oniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
  // oniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);

  // // dimuon
  oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
  oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);
  oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge);
  oniaTree->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D);
  oniaTree->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau);

  // single-muon
  // oniaTree->SetBranchAddress("Gen_mu_normChi2_global", Gen_mu_normChi2_global);
  // oniaTree->SetBranchAddress("Gen_mu_StationsMatched", Gen_mu_StationsMatched);
  // oniaTree->SetBranchAddress("Gen_mu_dxy", Gen_mu_dxy);
  // oniaTree->SetBranchAddress("Gen_mu_dxyErr", Gen_mu_dxyErr);
  // oniaTree->SetBranchAddress("Gen_mu_dz", Gen_mu_dz);
  // oniaTree->SetBranchAddress("Gen_mu_dzErr", Gen_mu_dzErr);

  // MC only
  if (isGenW)
  {
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
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

  const char *outDir = "skim_roots";
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] failed to create output directory: " << outDir << "\n";
    return;
  }

  std::string outFilePath = Form(
      "%s/acc_OO2025_isMC1_%s%s_Dimuon_MiniAOD_Private_MC.root",
      outDir, mcLabel.c_str(), weightLabel.c_str());

  TFile *fFlowSkim = TFile::Open(outFilePath.c_str(), "RECREATE");
  if (!fFlowSkim || fFlowSkim->IsZombie() || !fFlowSkim->IsWritable())
  {
    cout << "[ERROR] cannot open output file for writing: " << outFilePath << "\n";
    return;
  }

  const double ptBinsMid[] = {7, 8, 9, 10, 12, 14, 16, 20, 35};
  const double ptBinsFwd[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20};
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;

  TH1D *hist_den_mid = new TH1D("hist_den_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_den_fwd = new TH1D("hist_den_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_num_mid = new TH1D("hist_num_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_num_fwd = new TH1D("hist_num_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);

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

  // ---- preaparation before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  if (isNCollW)
  {
    cout << "[ERROR] OO NCollWeight Not Ready\n";
    return;
  }

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if ((iev % 100000 == 0) && (iev != 0))
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Saved dimuon: " << count_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total saved dimuon: " << count_dimuon << "\n";
    }

    // read event in Oniatree
    oniaTree->GetEntry(iev);

    // cout << "is Gen_QQ_4mom_pt->size() == Gen_QQ_size: " << (n == qq_size) << endl;
    
    // event-level values
    runN = runNb;
    event = eventNb;
    lumi = LS;
    vz = zVtx;

    double Ncoll_weight = 1.0;
    if (isNCollW)
    {
      Ncoll_weight = 1.0;
    }
    double Gen_weight_ = 1.0;
    if (isGenW)
    {
      Gen_weight_ = Gen_weight;
    }
    const double eventWeightBase = Ncoll_weight * Gen_weight_;

    // check dimoun number
    if (Gen_QQ_size<0) continue;

    // if (TMath::Abs(vz) > 15) continue; // not included in the Jun Chen's slide
    
    // ----- dimuon loop -----
    // reset number of dimuon
    nDimu = 0;

    for (Int_t irqq = 0; irqq < Gen_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      // check TLorentzVector pointers
      const int iMuPl = Gen_QQ_mupl_idx[irqq], iMuMi = Gen_QQ_mumi_idx[irqq]; // index guard

      // user-defined kinematic cuts
      float mass_ = Gen_QQ_4mom_m->at(irqq);

      float eta_ = Gen_QQ_4mom_eta->at(irqq);
      float eta1_ = Gen_mu_4mom_eta->at(iMuPl);
      float eta2_ = Gen_mu_4mom_eta->at(iMuMi);
      
      float pt_ = Gen_QQ_4mom_pt->at(irqq);
      float pt1_ = Gen_mu_4mom_pt->at(iMuPl);
      float pt2_ = Gen_mu_4mom_pt->at(iMuMi);
      
      float phi_ = Gen_QQ_4mom_phi->at(irqq);
      float phi1_ = Gen_mu_4mom_phi->at(iMuPl);
      float phi2_ = Gen_mu_4mom_phi->at(iMuMi);

      // reconstrcut TLorentz vector to get the rapidityu
      TLorentzVector JP;
      JP.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
      float y_ = JP.Rapidity();
      const bool isMidRap = std::abs(y_) < 1.6;
      TF1 *fPtW = isMidRap ? fPtW_mid : fPtW_fwd;
      const double ptWeight = isPtW ? GetPtWeight(fPtW, pt_, y_) : 1.0;
      
      if (mass_ < massLow || mass_ > massHigh) continue;
      if (Gen_mu_charge[iMuPl] * Gen_mu_charge[iMuMi] > 0) continue;

      // ===== Acceptance denominator cut =====
      if (!(TMath::Abs(y_) < 2.4))
        continue;
      
      // fill the denominator
      if (isMidRap) {
        const int ptBin = hist_den_mid->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_den_mid->GetNbinsX())
          hist_den_mid->Fill(pt_, eventWeightBase * ptWeight);
      } else {
        const int ptBin = hist_den_fwd->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_den_fwd->GetNbinsX())
          hist_den_fwd->Fill(pt_, eventWeightBase * ptWeight);
      }
      
      // ===== Acceptance numerator =====
      if ( !IsAcceptanceQQJunChenOO(pt1_, eta1_) ||
           !IsAcceptanceQQJunChenOO(pt2_, eta2_)) continue;
      
      // fill the numerator
      if (isMidRap) {
        const int ptBin = hist_num_mid->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_num_mid->GetNbinsX())
          hist_num_mid->Fill(pt_, eventWeightBase * ptWeight);
      } else {
        const int ptBin = hist_num_fwd->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_num_fwd->GetNbinsX())
          hist_num_fwd->Fill(pt_, eventWeightBase * ptWeight);
      }


      // init weights
      weight = eventWeightBase;
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      ++nDimu;
    } // end of dimuon loop

    if (nDimu > 0)
    {
      count_dimuon += nDimu;
      // flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== compute acceptance =====
  hist_den_mid->Scale(1.0, "width");
  hist_den_fwd->Scale(1.0, "width");
  hist_num_mid->Scale(1.0, "width");
  hist_num_fwd->Scale(1.0, "width");

  // h_acc_fwd = h_num_fwd / h_den_fwd
  TH1D *hist_acc_fwd = (TH1D *)hist_num_fwd->Clone("hist_acc_fwd");
  hist_acc_fwd->SetTitle(";p_{T} (GeV/c);Acceptance");
  hist_acc_fwd->Divide(hist_num_fwd, hist_den_fwd, 1.0, 1.0, "B");

  // h_acc_mid = h_num_mid / h_den_mid
  TH1D *hist_acc_mid = (TH1D *)hist_num_mid->Clone("hist_acc_mid");
  hist_acc_mid->SetTitle(";p_{T} (GeV/c);Acceptance");
  hist_acc_mid->Divide(hist_num_mid, hist_den_mid, 1.0, 1.0, "B");

  // ===== save results =====
  fFlowSkim->cd();
  hist_den_mid->Write();
  hist_den_fwd->Write();
  hist_num_mid->Write();
  hist_num_fwd->Write();
  hist_acc_mid->Write();
  hist_acc_fwd->Write();
  // flowTree->Write("myTree");
  fFlowSkim->Close();

  // ----- print -----
  cout << "Total saved dimuon: " << count_dimuon << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
