#include "TBranch.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TLeaf.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using std::cout;

namespace
{
  TString resolveBaseDir()
  {
    return gSystem->DirName(__FILE__);
  }

  bool IsAcceptanceQQRun2018(double pt, double eta)
  {
    const double aeta = std::abs(eta);

    return ((aeta < 1.2 && pt >= 3.5) ||
            (1.2 <= aeta && aeta < 2.1 && pt >= 5.47 - 1.89 * aeta) ||
            (2.1 <= aeta && aeta < 2.4 && pt >= 1.5));
  }

  double GetPtWeight(TF1 *fW, double pt)
  {
    return fW ? fW->Eval(pt) : 1.0;
  }

  bool BranchUsesShort(TTree *tree, const char *branchName)
  {
    if (!tree)
      return false;
    TBranch *branch = tree->GetBranch(branchName);
    if (!branch)
      return false;
    TLeaf *leaf = branch->GetLeaf(branchName);
    if (!leaf && branch->GetListOfLeaves() && branch->GetListOfLeaves()->GetEntries() > 0)
      leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->At(0));
    if (!leaf)
      return false;
    const std::string typeName = leaf->GetTypeName();
    return typeName == "Short_t" || typeName == "UShort_t";
  }

  TH1D *BuildAcceptanceHist(TH1D *num, TH1D *den, const char *name, const char *title)
  {
    if (!num || !den)
      return nullptr;

    TH1D *acc = static_cast<TH1D *>(num->Clone(name));
    acc->SetTitle(title);
    acc->Divide(num, den, 1.0, 1.0, "B");
    return acc;
  }
} // namespace

void acceptance_1d(long nEvt = 100000, bool isPr = true, bool isPbPb = false, bool isGenW = true, bool isPtW = true)
{
  cout << "Start acceptance_1d()\n";
  TH1::SetDefaultSumw2(kTRUE);

  TChain *oniaTree = new TChain("hionia/myTree");
  const std::string inputFile = isPr
                                    ? "/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root"
                                    : "/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root";
  cout << "[INFO] input MC is fixed to pp source: " << inputFile << "\n";
  if (oniaTree->Add(inputFile.c_str()) <= 0)
  {
    cout << "[ERROR] failed to add input file: " << inputFile << "\n";
    return;
  }

  TF1 *fPtW_mid = nullptr;
  TF1 *fPtW_fwd = nullptr;
  if (isPtW)
  {
    const std::string ptWeightPath =
        isPbPb
            ? (isPr
                   ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310.root"
                   : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310.root")
            : (isPr
                   ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310.root"
                   : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root");

    TFile *fMid = TFile::Open(ptWeightPath.c_str(), "READ");
    TFile *fFwd = TFile::Open(ptWeightPath.c_str(), "READ");
    if (!fMid || fMid->IsZombie() || !fFwd || fFwd->IsZombie())
    {
      cout << "[ERROR] failed to open pT weight file: " << ptWeightPath << "\n";
      return;
    }

    TF1 *fMidIn = dynamic_cast<TF1 *>(fMid->Get("dataMC_Ratio1"));
    TF1 *fFwdIn = dynamic_cast<TF1 *>(fFwd->Get("dataMC_Ratio1"));
    if (!fMidIn || !fFwdIn)
    {
      cout << "[ERROR] failed to load pT weight function dataMC_Ratio1 from " << ptWeightPath << "\n";
      return;
    }

    fPtW_mid = static_cast<TF1 *>(fMidIn->Clone("acc_run2_ptw_mid"));
    fPtW_fwd = static_cast<TF1 *>(fFwdIn->Clone("acc_run2_ptw_fwd"));
    fMid->Close();
    fFwd->Close();
    delete fMid;
    delete fFwd;
  }

  const bool useShortSizeBranches = BranchUsesShort(oniaTree, "Gen_QQ_size");
  const bool useShortIndexBranches = BranchUsesShort(oniaTree, "Gen_QQ_mupl_idx");
  const bool usesSplitMomentumBranches = (oniaTree->GetBranch("Gen_QQ_4mom_pt") != nullptr);
  const bool hasGenWeightBranch = (oniaTree->GetBranch("Gen_weight") != nullptr);

  TClonesArray *Gen_QQ_4mom = nullptr;
  TClonesArray *Gen_mu_4mom = nullptr;
  std::vector<float> *Gen_QQ_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_pt = nullptr;
  std::vector<float> *Gen_QQ_4mom_phi = nullptr;
  std::vector<float> *Gen_QQ_4mom_eta = nullptr;
  std::vector<float> *Gen_mu_4mom_pt = nullptr;
  std::vector<float> *Gen_mu_4mom_phi = nullptr;
  std::vector<float> *Gen_mu_4mom_eta = nullptr;

  if (usesSplitMomentumBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_4mom_m", &Gen_QQ_4mom_m);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_pt", &Gen_QQ_4mom_pt);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_phi", &Gen_QQ_4mom_phi);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_eta", &Gen_QQ_4mom_eta);
    oniaTree->SetBranchAddress("Gen_mu_4mom_pt", &Gen_mu_4mom_pt);
    oniaTree->SetBranchAddress("Gen_mu_4mom_phi", &Gen_mu_4mom_phi);
    oniaTree->SetBranchAddress("Gen_mu_4mom_eta", &Gen_mu_4mom_eta);
  }
  else
  {
    oniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
    oniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
  }

  Short_t Gen_QQ_size_short = 0;
  Int_t Gen_QQ_size_int = 0;
  Short_t Gen_QQ_mupl_idx_short[1000];
  Short_t Gen_QQ_mumi_idx_short[1000];
  Short_t Gen_mu_charge_short[1000];
  Int_t Gen_QQ_mupl_idx_int[1000];
  Int_t Gen_QQ_mumi_idx_int[1000];
  Int_t Gen_mu_charge_int[1000];
  Float_t Gen_weight = 1.f;

  if (useShortSizeBranches)
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_short);
  else
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_int);

  if (useShortIndexBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_short);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_short);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_short);
  }
  else
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_int);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_int);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_int);
  }

  if (isGenW && hasGenWeightBranch)
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);

  auto GetGenQQSize = [&]() -> Int_t
  {
    return useShortSizeBranches ? static_cast<Int_t>(Gen_QQ_size_short) : Gen_QQ_size_int;
  };
  auto GetGenQQMuplIdx = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_QQ_mupl_idx_short[idx]) : Gen_QQ_mupl_idx_int[idx];
  };
  auto GetGenQQMumiIdx = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_QQ_mumi_idx_short[idx]) : Gen_QQ_mumi_idx_int[idx];
  };
  auto GetGenMuCharge = [&](Int_t idx) -> Int_t
  {
    return useShortIndexBranches ? static_cast<Int_t>(Gen_mu_charge_short[idx]) : Gen_mu_charge_int[idx];
  };

  const TString outDir = resolveBaseDir() + "/skim_roots";
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] failed to create output directory: " << outDir << "\n";
    return;
  }

  const std::string collLabel = isPbPb ? "PbPb2018" : "pp2018";
  const std::string mcLabel = isPr ? "PR" : "NP";
  const std::string outFilePath = Form("%s/acc_%s_ppInput_isMC1_%s_ncollW0_genW%d_ptW%d.root",
                                       outDir.Data(), collLabel.c_str(), mcLabel.c_str(), isGenW, isPtW);
  TFile *fout = TFile::Open(outFilePath.c_str(), "RECREATE");
  if (!fout || fout->IsZombie() || !fout->IsWritable())
  {
    cout << "[ERROR] cannot open output file for writing: " << outFilePath << "\n";
    return;
  }

  const double ptBinsMid[] = {6.5, 9, 12, 15, 20, 25, 40};
  const double ptBinsFwd[] = {3.5, 6.5, 9, 12, 40};
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;

  TH1D *hist_acc_den_mid = new TH1D("hist_acc_den_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_acc_den_fwd = new TH1D("hist_acc_den_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_acc_num_mid = new TH1D("hist_acc_num_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_acc_num_fwd = new TH1D("hist_acc_num_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);

  const double massLow = 2.6;
  const double massHigh = 3.5;
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt < 0 || nEvt > nTot)
    nEvt = nTot;

  long denCount = 0;
  long numCount = 0;
  for (Long64_t iev = 0; iev < nEvt; ++iev)
  {
    if ((iev % 100000 == 0) && iev != 0)
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
           << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Acceptance denominator: " << denCount << "\n";
      cout << "Acceptance numerator: " << numCount << "\n";
    }

    oniaTree->GetEntry(iev);

    const Int_t Gen_QQ_size = GetGenQQSize();
    if (Gen_QQ_size < 0)
      continue;

    double eventWeightBase = 1.0;
    if (isGenW && hasGenWeightBranch)
      eventWeightBase *= Gen_weight;

    for (Int_t irqq = 0; irqq < Gen_QQ_size; ++irqq)
    {
      const Int_t iMuPl = GetGenQQMuplIdx(irqq);
      const Int_t iMuMi = GetGenQQMumiIdx(irqq);
      if (iMuPl < 0 || iMuMi < 0)
        continue;

      float mass = 0.f;
      float pt = 0.f;
      float eta = 0.f;
      float phi = 0.f;
      float pt1 = 0.f;
      float pt2 = 0.f;
      float eta1 = 0.f;
      float eta2 = 0.f;

      if (usesSplitMomentumBranches)
      {
        if (!Gen_QQ_4mom_m || !Gen_QQ_4mom_pt || !Gen_QQ_4mom_eta || !Gen_QQ_4mom_phi ||
            !Gen_mu_4mom_pt || !Gen_mu_4mom_eta)
          continue;
        mass = Gen_QQ_4mom_m->at(irqq);
        pt = Gen_QQ_4mom_pt->at(irqq);
        eta = Gen_QQ_4mom_eta->at(irqq);
        phi = Gen_QQ_4mom_phi->at(irqq);
        pt1 = Gen_mu_4mom_pt->at(iMuPl);
        pt2 = Gen_mu_4mom_pt->at(iMuMi);
        eta1 = Gen_mu_4mom_eta->at(iMuPl);
        eta2 = Gen_mu_4mom_eta->at(iMuMi);
      }
      else
      {
        auto *qq = static_cast<TLorentzVector *>(Gen_QQ_4mom->At(irqq));
        auto *muPl = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuPl));
        auto *muMi = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuMi));
        if (!qq || !muPl || !muMi)
          continue;
        mass = qq->M();
        pt = qq->Pt();
        eta = qq->Eta();
        phi = qq->Phi();
        pt1 = muPl->Pt();
        pt2 = muMi->Pt();
        eta1 = muPl->Eta();
        eta2 = muMi->Eta();
      }

      if (mass < massLow || mass > massHigh)
        continue;
      if (GetGenMuCharge(iMuPl) * GetGenMuCharge(iMuMi) > 0)
        continue;

      TLorentzVector qq;
      qq.SetPtEtaPhiM(pt, eta, phi, mass);
      const float absY = std::abs(qq.Rapidity());
      if (absY >= 2.4f)
        continue;

      bool isMid = false;
      bool isFwd = false;
      if (absY < 1.6f && pt >= 6.5f && pt < 40.f)
        isMid = true;
      if (absY >= 1.6f && absY < 2.4f && pt >= 3.5f && pt < 40.f)
        isFwd = true;
      if (!isMid && !isFwd)
        continue;

      const double ptWeight = isPtW ? GetPtWeight(isMid ? fPtW_mid : fPtW_fwd, pt) : 1.0;
      const double weight = eventWeightBase * ptWeight;

      if (isMid)
      {
        hist_acc_den_mid->Fill(pt, weight);
      }
      else
      {
        hist_acc_den_fwd->Fill(pt, weight);
      }
      ++denCount;

      if (!IsAcceptanceQQRun2018(pt1, eta1) || !IsAcceptanceQQRun2018(pt2, eta2))
        continue;

      if (isMid)
      {
        hist_acc_num_mid->Fill(pt, weight);
      }
      else
      {
        hist_acc_num_fwd->Fill(pt, weight);
      }
      ++numCount;
    }
  }

  hist_acc_den_mid->Scale(1.0, "width");
  hist_acc_den_fwd->Scale(1.0, "width");
  hist_acc_num_mid->Scale(1.0, "width");
  hist_acc_num_fwd->Scale(1.0, "width");

  TH1D *hist_acc_mid = BuildAcceptanceHist(hist_acc_num_mid, hist_acc_den_mid, "hist_acc_mid", ";p_{T} (GeV/c);Acceptance");
  TH1D *hist_acc_fwd = BuildAcceptanceHist(hist_acc_num_fwd, hist_acc_den_fwd, "hist_acc_fwd", ";p_{T} (GeV/c);Acceptance");

  fout->cd();
  hist_acc_den_mid->Write();
  hist_acc_den_fwd->Write();
  hist_acc_num_mid->Write();
  hist_acc_num_fwd->Write();
  hist_acc_mid->Write();
  hist_acc_fwd->Write();
  fout->Close();

  cout << "Acceptance denominator count: " << denCount << "\n";
  cout << "Acceptance numerator count: " << numCount << "\n";
  cout << "Output ROOT: " << outFilePath << "\n";
  cout << "Finish acceptance_1d()\n";
}
