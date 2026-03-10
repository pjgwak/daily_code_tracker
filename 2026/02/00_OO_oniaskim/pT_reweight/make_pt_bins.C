// For PbPb2023 data

#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TSystem.h"
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

namespace
{
// ----------------------------------------------------------------------
// Save one pT-binned histogram with a CMS-style single-panel layout.
// ----------------------------------------------------------------------
void savePtHistPlot(TH1D *hist, const TString &outPath, const TString &kinLabel, const TString &sampleLabel)
{
  TCanvas canv("c_fill_bins", "", 800, 800);
  canv.cd();

  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.Draw();
  pad1.cd();

  hist->SetTitle("");
  hist->SetLineColor(kAzure + 2);
  hist->SetMarkerColor(kAzure + 2);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.7);
  hist->SetLineWidth(3);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetXaxis()->SetLimits(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  double maxY = hist->GetMaximum();
  if (maxY <= 0.0)
    maxY = 1.0;
  hist->SetMinimum(0.0);
  hist->SetMaximum(maxY * 2.0);
  hist->Draw("E1");

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.040);
    tx.SetTextFont(72);
    tx.DrawLatex(0.19, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.030);
    tx.SetTextFont(42);
    tx.DrawLatex(0.19, 0.865, sampleLabel);
    tx.DrawLatex(0.19, 0.815, kinLabel);
  }

  TLegend leg(0.50, 0.78, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.024);
  leg.SetEntrySeparation(0.01);
  leg.AddEntry(hist, hist->GetYaxis()->GetTitle(), "lp");
  leg.Draw();

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

// ----------------------------------------------------------------------
// Muon acceptance helper.
// The 2018 version is kept below as a reference, while the Run-3 version
// is the one currently used in this workflow.
// ----------------------------------------------------------------------
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

void make_pt_bins(long nEvt = 1000, bool isMC = false, bool isPr = true)
{
  cout << "Start onia_to_skim_data()\n";

  // --------------------------------------------------------------------
  // Global ROOT drawing setup.
  // --------------------------------------------------------------------
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // --------------------------------------------------------------------
  // Build the input TChain and choose the source file from the run mode.
  //   data        : inclusive OO data
  //   MC + isPr   : prompt MC
  //   MC + !isPr  : non-prompt MC
  // --------------------------------------------------------------------
  TChain *oniaTree = new TChain("hionia/myTree");
  
  std::string inputFile = "/data/Oniatree/light_ions_Raa/OO_Data/OniaTree_OODimuon_MINIAOD_Run2025OO_PromptReco_v1_Jul12_merged.root";
  if (isMC)
  {
    inputFile = isPr
        ? "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcPr.root"
        : "/data/Oniatree/light_ions_Raa/OO_PrivateMc/Oniatree_OO2025PrivateMcNp.root";
  }
  oniaTree->Add(inputFile.c_str());

  // --------------------------------------------------------------------
  // Analysis-level selection parameters used in the event loop.
  // --------------------------------------------------------------------
  int cLow = 0, cHigh = 180;
  double massLow = 2.6, massHigh = 3.5;

  // --------------------------------------------------------------------
  // Output-file label for the MC flavor.
  // --------------------------------------------------------------------
  std::string mcLabel = "";
  if (isMC)
    mcLabel = isPr ? "PR" : "NP";

  // --------------------------------------------------------------------
  // Set input branch addresses from the onia tree.
  // This block defines the exact information needed for candidate selection
  // and for constructing the pT-binned output histograms.
  // --------------------------------------------------------------------
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
    // oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
    if (oniaTree->GetBranch("Reco_mu_whichGen"))
      oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }

  // --------------------------------------------------------------------
  // Declare local output-side buffers.
  // Most of the TTree branches are currently left in place for possible
  // skim production, even though only the histograms are written here.
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // Prepare output ROOT file and histograms.
  //   hist_raw_* : unnormalized candidate counts, used later for error work
  //   hist_eff_* : normalized shape histograms used for direct plotting
  // --------------------------------------------------------------------
  const TString outBaseDir = "fill_bins_outputs";
  const TString outRootsDir = outBaseDir + "/roots";
  const TString outFigsDir = outBaseDir + "/figs";
  gSystem->mkdir(outRootsDir, true);
  gSystem->mkdir(outFigsDir, true);
  TFile *fFlowSkim = new TFile(
      Form("%s/pt_bins_OO2025_isMC%d%s_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root",
           outRootsDir.Data(),
           isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str()),
      "recreate");

  const double ptBinsMid[] = {7, 8, 9, 10, 12, 14, 16, 20, 35};
  const double ptBinsFwd[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20};
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;

  TH1D *hist_eff_mid = new TH1D("hist_eff_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_fwd = new TH1D("hist_eff_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_raw_mid = new TH1D("hist_raw_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_raw_fwd = new TH1D("hist_raw_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  hist_eff_mid->SetDirectory(nullptr);
  hist_eff_fwd->SetDirectory(nullptr);
  hist_raw_mid->SetDirectory(nullptr);
  hist_raw_fwd->SetDirectory(nullptr);
  hist_eff_mid->Sumw2();
  hist_eff_fwd->Sumw2();
  hist_raw_mid->Sumw2();
  hist_raw_fwd->Sumw2();
  hist_eff_mid->GetYaxis()->SetTitle("Normalized counts");
  hist_eff_fwd->GetYaxis()->SetTitle("Normalized counts");
  hist_raw_mid->GetYaxis()->SetTitle("Counts");
  hist_raw_fwd->GetYaxis()->SetTitle("Counts");
  

  // --------------------------------------------------------------------
  // Optional skim tree definition. The branches are still declared so the
  // macro can be extended back into a skim-maker if needed.
  // --------------------------------------------------------------------
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

  // --------------------------------------------------------------------
  // Event-loop bookkeeping.
  // --------------------------------------------------------------------
  long count_reco_dimuon = 0;

  // --------------------------------------------------------------------
  // Resolve the number of entries to process.
  // --------------------------------------------------------------------
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  // --------------------------------------------------------------------
  // Main event loop.
  // --------------------------------------------------------------------
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if ((iev % 100000 == 0) && (iev != 0))
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Reco selected dimuon: " << count_reco_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
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

    // ----------------------------------------------------------------
    // Event-level quantities. Centrality-related code is kept commented
    // because this pT-binning workflow currently does not use it.
    // ----------------------------------------------------------------
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
    
    // MC-only global weights can be inserted here when needed.
    double Ncoll_weight = 0, Gen_weight_ = 0;
    // if(isMC){
    //   weight = findNcoll(Centrality) * Gen_weight;
    //   Ncoll_weight = findNcoll(Centrality);
    //   Gen_weight_ = Gen_weight;
    // }

    // check dimoun number
    if (Reco_QQ_size<0) continue;

    // apply HLT trigger
    int kTrigSel = 2; // HLT_OxyL1SingleMuOpen_v1
    // Keep the original trigger-mask expression unchanged.
    // ACLiC may report a signedness warning here, but do not modify this line
    // unless the trigger-definition convention itself is intentionally updated.
    if (!((HLTriggers & ((int)std::pow(2, kTrigSel))) == ((int)std::pow(2, kTrigSel))))
    {
      continue;
    }

    // if (TMath::Abs(vz) > 15) continue; // not included in the Jun Chen's slide

    // ----------------------------------------------------------------
    // Reconstructed dimuon loop.
    // ----------------------------------------------------------------
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

      // ----------------------------------------------------------------
      // Read dimuon and single-muon kinematics for this candidate.
      // ----------------------------------------------------------------
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


      // Reconstruct the dimuon four-vector to compute rapidity.
      TLorentzVector JP;
      JP.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
      float y_ = JP.Rapidity();
      
      if(mass_ < massLow || mass_ > massHigh) continue;

      // Require reco-gen matching for MC templates.
      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      
      // Apply the single-muon quality and acceptance selections.
      if (!Reco_mu_isSoftCutBased[iMuPl] || !Reco_mu_isSoftCutBased[iMuMi]) continue;

      // dimuon y < 2.4 and single muon acceptance cut - 2018 soft
      if ( !(TMath::Abs(y_) < 2.4) ||
           !IsAcceptanceQQ(pt1_, eta1_) ||
           !IsAcceptanceQQ(pt2_, eta2_)) continue;

      // Require both muons to be global muons.
      if (!Reco_mu_isGlobal[iMuPl] || !Reco_mu_isGlobal[iMuMi]) continue;

      // Apply the dimuon vertex quality selection.
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue; // not included in JunChen's slide

      // Candidate weight placeholder. This stays at unity for now.
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      if (Reco_QQ_sign[irqq] != 0) continue; // basic

      // recoQQsign[nDimu] = Reco_QQ_sign[irqq];

      // ----------------------------------------------------------------
      // Fill the raw-count histogram and the plotting histogram for the
      // appropriate rapidity region.
      // ----------------------------------------------------------------
      if (std::abs(y_) < 1.6)
      {
        const int ptBin = hist_eff_mid->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_eff_mid->GetNbinsX())
        {
          hist_eff_mid->Fill(pt_, weight);
          hist_raw_mid->Fill(pt_, weight);
        }
      }
      else
      {
        const int ptBin = hist_eff_fwd->FindBin(pt_);
        if (ptBin > 0 && ptBin <= hist_eff_fwd->GetNbinsX())
        {
          hist_eff_fwd->Fill(pt_, weight);
          hist_raw_fwd->Fill(pt_, weight);
        }
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

  // --------------------------------------------------------------------
  // Convert the plotting histograms into normalized pT shapes.
  // The raw histograms remain untouched and preserve count-based errors.
  // --------------------------------------------------------------------
  const double intMid = hist_eff_mid->Integral(1, hist_eff_mid->GetNbinsX());
  if (intMid > 0.0)
    hist_eff_mid->Scale(1.0 / intMid, "width");

  const double intFwd = hist_eff_fwd->Integral(1, hist_eff_fwd->GetNbinsX());
  if (intFwd > 0.0)
    hist_eff_fwd->Scale(1.0 / intFwd, "width");

  // --------------------------------------------------------------------
  // Save both raw and normalized histograms to the output ROOT file.
  // --------------------------------------------------------------------
  fFlowSkim->cd();
  hist_raw_mid->Write();
  hist_raw_fwd->Write();
  hist_eff_mid->Write();
  hist_eff_fwd->Write();
  // flowTree->Write("myTree");
  fFlowSkim->Close();

  // --------------------------------------------------------------------
  // Save quick-look figures for raw counts and normalized shapes.
  // --------------------------------------------------------------------
  const TString sampleLabel = isMC ? (isPr ? "Prompt J/#psi MC" : "Nonprompt J/#psi MC")
                                   : "Inclusive J/#psi data";
  savePtHistPlot(hist_raw_mid, outFigsDir + Form("/mid_raw_isMC%d%s", isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str()),
                 "|y| < 1.6", sampleLabel);
  savePtHistPlot(hist_raw_fwd, outFigsDir + Form("/fwd_raw_isMC%d%s", isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str()),
                 "1.6 < |y| < 2.4", sampleLabel);
  savePtHistPlot(hist_eff_mid, outFigsDir + Form("/mid_norm_isMC%d%s", isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str()),
                 "|y| < 1.6", sampleLabel);
  savePtHistPlot(hist_eff_fwd, outFigsDir + Form("/fwd_norm_isMC%d%s", isMC, mcLabel.empty() ? "" : ("_" + mcLabel).c_str()),
                 "1.6 < |y| < 2.4", sampleLabel);

  // --------------------------------------------------------------------
  // Print a minimal processing summary.
  // --------------------------------------------------------------------
  cout << "Total Reco selected dimuon: " << count_reco_dimuon << "\n";
  cout << "Saved ROOT file under: " << outRootsDir << "\n";
  cout << "Saved plots under: " << outFigsDir << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
