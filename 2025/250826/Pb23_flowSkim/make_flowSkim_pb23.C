// For PbPb2023 data

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>

using std::cout; using std::string;

static const long MAXTREESIZE = 1000000000000;

void make_flowSkim_pb23(bool isMC = false, long nEvt = -1)
{
  cout << "Start make_flowSkim_pb23()\n";

  // ===== read input =====
  TFile *fInput = TFile::Open("/data/users/pjgwak/input_files/PbPb2023_Data_noTrack/Oniatree_2023PbPbPromptRecoData_132X_miniAOD_noTrack.root", "read");
  TTree *oniaTree = (TTree *)fInput->Get("hionia/myTree");

  // ===== labeling =====
  // PD type
  // Data vs MC
  // Trigger
  // Muon Selection
  // kinematics (?)

  // TString fPD;
  // if (PDtype == 1)
  //   fPD = "DoubleMuon";
  // else if (PDtype == 2)
  //   fPD = "DBPeri";

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 1000;

  // event-level scalars
  UInt_t runNb;
  UInt_t eventNb;
  UInt_t LS;
  float zVtx;
  Int_t Centrality;
  ULong64_t HLTriggers;
  Float_t SumET_HF;
  Float_t Gen_weight; // MC

  // collection sizes
  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;

  // TClonesArray pointer
  TClonesArray *Reco_QQ_4mom = nullptr;
  TClonesArray *Reco_mu_4mom = nullptr;

  // per-dimuon (size = Reco_QQ_size)
  ULong64_t Reco_QQ_trig[maxBranchSize];
  Float_t Reco_QQ_VtxProb[maxBranchSize];
  Short_t Reco_QQ_mupl_idx[maxBranchSize];
  Short_t Reco_QQ_mumi_idx[maxBranchSize];
  Short_t Reco_QQ_sign[maxBranchSize];
  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctauErr3D[maxBranchSize];

  // per-muon (size = Reco_mu_size)
  ULong64_t Reco_mu_trig[maxBranchSize];
  Bool_t Reco_mu_highPurity[maxBranchSize];
  Int_t Reco_mu_nTrkHits[maxBranchSize];
  Float_t Reco_mu_normChi2_global[maxBranchSize];
  Int_t Reco_mu_nMuValHits[maxBranchSize];
  Int_t Reco_mu_StationsMatched[maxBranchSize];
  Float_t Reco_mu_dxy[maxBranchSize];
  Float_t Reco_mu_dxyErr[maxBranchSize];
  Float_t Reco_mu_dz[maxBranchSize];
  Float_t Reco_mu_dzErr[maxBranchSize];
  Int_t Reco_mu_nTrkWMea[maxBranchSize];
  Int_t Reco_mu_nPixWMea[maxBranchSize];
  Int_t Reco_mu_nPixValHits[maxBranchSize];
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Int_t Reco_mu_whichGen[maxBranchSize]; // MC

  // ----- SetBranchAddress -----
  oniaTree->SetBranchAddress("runNb", &runNb);
  oniaTree->SetBranchAddress("eventNb", &eventNb);
  oniaTree->SetBranchAddress("LS", &LS);
  oniaTree->SetBranchAddress("zVtx", &zVtx);
  oniaTree->SetBranchAddress("Centrality", &Centrality);
  oniaTree->SetBranchAddress("HLTriggers", &HLTriggers);
  oniaTree->SetBranchAddress("SumET_HF", &SumET_HF);

  // sizes
  oniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  oniaTree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);

  // TLorentzVecotr
  oniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  oniaTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);

  // dimuon
  oniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  oniaTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  oniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  oniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  oniaTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  oniaTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);

  // single-muon
  oniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  oniaTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  oniaTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  oniaTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  oniaTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  oniaTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  oniaTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  oniaTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  oniaTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  oniaTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  oniaTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  oniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);

  // MC only
  if (isMC)
  {
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
    oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }

  // ===== set output file and  FlowSkim branch address (output) =====
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
  Int_t recoQQsign[nMaxDimu];
  Float_t ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];

  // weights
  Double_t weight = 1.0; // event-level weight
  Double_t TnPweight[nMaxDimu]; // per-candidate level weight

  // TLorentzVector
  TLorentzVector *JP_Reco = new TLorentzVector;
  TLorentzVector *mupl_Reco = new TLorentzVector;
  TLorentzVector *mumi_Reco = new TLorentzVector;

  TFile *fFlowSkim = new TFile("output_test.root", "recreate");
  // TFile *fFlowSkim = new TFile(
  //     Form("OniaFlowSkim_%sTrig_%sPD_PbPb_isMC%d_250826.root",
  //          fTrigName[trigIndx].Data(), fPD.Data(), isMC),
  //     "RECREATE");

  // compression and autosiave options
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

  // weights
  flowTree->Branch("weight", &weight, "weight/D");
  flowTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");

  // ===== event loop =====
  // counters
  long count_dimuon = 0;

  // ---- preapare before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = min<Long64_t>(nEvt, nTot);

  // ----- bit mask - change here -----
  int trigIndx = 0;
  int dummpy_kTrigSel = 0; // for test - it's not used for MinBias
  // if (kTrigSel == kTrigJpsipp) trigIndx = 0; // PbPb : 0
  // else if (kTrigSel == kTrigUps) trigIndx = 1;
  // else if (kTrigSel == kTrigL1DBOS40100) trigIndx = 2;
  // else if (kTrigSel == kTrigL1DB50100) trigIndx = 3;

  int kL2filter = 16;
  int kL3filter = 17;

  const ULong64_t HLT_MASK = (1ULL << dummpy_kTrigSel);
  const ULong64_t L2_MASK = (1ULL << kL2filter);
  const ULong64_t L3_MASK = (1ULL << kL3filter);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if (iev % 100000 == 0)
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / max<Long64_t>(1, nTot)) << "%)\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << endl;
      cout << "Total saved dimuon: " << count_dimuon << "\n";
    }

    // read event in Oniatree
    oniaTree->GetEntry(iev);

    // event-level values
    runN = runNb;
    event = eventNb;
    lumi = LS;
    vz = zVtx;
    // cBin = -999;
    // if(hiHFBinEdge ==0) cBin = getHiBinFromhiHF(SumET_HF);
    // else if(hiHFBinEdge == 1) cBin = getHiBinFromhiHF_Up(SumET_HF);
    // else if(hiHFBinEdge == -1) cBin = getHiBinFromhiHF_Down(SumET_HF);
    // if(cBin==-999){ cout << "ERROR!!! No HF Centrality Matching!!" << endl; return;}
    

    // apply HLT - not used for PbPb2023 (minBias)
    // if ((HLTriggers & HLT_MASK) == 0ULL) continue;
    
    // check dimoun number
    if (Reco_QQ_size<0) continue;

    // check TClonesArray
    if (!Reco_QQ_4mom || !Reco_mu_4mom) continue;
    
    // ----- dimuon loop -----
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

      // apply HLT to dimuon candidate - not used??? - check later
      // if ((Reco_QQ_trig[irqq] & HLT_MASK) == 0ULL) continue;

      // check TLorentzVector pointers
      const int iMuPl = Reco_QQ_mupl_idx[irqq], iMuMi = Reco_QQ_mumi_idx[irqq]; // index guard
      if (iMuPl < 0 || iMuMi < 0 || iMuPl >= Reco_mu_size || iMuMi >= Reco_mu_size) continue;
      if (Reco_QQ_4mom->GetEntriesFast() <= irqq) continue;
      if (Reco_mu_4mom->GetEntriesFast() <= iMuPl || Reco_mu_4mom->GetEntriesFast() <= iMuMi) continue;

      // bring TLorentzVector
      const auto &JP = *static_cast<TLorentzVector *>(Reco_QQ_4mom->At(irqq));
      const auto &mupl = *static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuPl));
      const auto &mumi = *static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuMi));

      // check MC Reco-Gen match 
      if (isMC)
      {
        if (Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if (Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }

      // selection bits
      const bool passMuonTypePl =
          ((Reco_mu_SelectionType[iMuPl] & (1 << 1)) != 0) &&
          ((Reco_mu_SelectionType[iMuPl] & (1 << 3)) != 0);
      const bool passMuonTypeMi =
          ((Reco_mu_SelectionType[iMuMi] & (1 << 1)) != 0) &&
          ((Reco_mu_SelectionType[iMuMi] & (1 << 3)) != 0);
      
      // soft id cut
      const bool muplSoft =
          (Reco_mu_nTrkWMea[iMuPl] > 5) &&
          (Reco_mu_nPixWMea[iMuPl] > 0) &&
          (std::fabs(Reco_mu_dxy[iMuPl]) < 0.3f) &&
          (std::fabs(Reco_mu_dz[iMuPl]) < 20.f) &&
          passMuonTypePl;
      const bool mumiSoft =
          (Reco_mu_nTrkWMea[iMuMi] > 5) &&
          (Reco_mu_nPixWMea[iMuMi] > 0) &&
          (std::fabs(Reco_mu_dxy[iMuMi]) < 0.3f) &&
          (std::fabs(Reco_mu_dz[iMuMi]) < 20.f) &&
          passMuonTypeMi;
      if (!(muplSoft && mumiSoft)) continue;

      // vertex probability cut
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];

      // ----- TnP (Skip now) -----
      // if (isMC)
      // {
      //   // apply TnP weight 
      //   tnp_weight *= tnp_weight_muid_pbpb(mupl.Pt(), mupl.Eta(), 0) * tnp_weight_muid_pbpb(mumi.Pt(), mumi.Eta(), 0); // muon id
      //   tnp_weight *= tnp_weight_trk_pbpb(mupl.Eta(), 0) * tnp_weight_trk_pbpb(mumi.Eta(), 0); //inner tracker

      //   // L2, L3 filter bits
      //   // Trigger part
      //   if (!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]] & ((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter))))
      //   {
      //     //         cout << "irqq : " << irqq << " - iev : " << iev << endl;
      //     //         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
      //     continue;
      //   }

      //   // should pass L2 filters
      //   const bool mupl_L2 = ((Reco_mu_trig[iMuPl] & L2_MASK) != 0ULL);
      //   const bool mupl_L3 = ((Reco_mu_trig[iMuPl] & L3_MASK) != 0ULL);
      //   const bool mumi_L2 = ((Reco_mu_trig[iMuMi] & L2_MASK) != 0ULL);
      //   const bool mumi_L3 = ((Reco_mu_trig[iMuMi] & L3_MASK) != 0ULL);
      //   if (!(mupl_L2 && mumi_L2)) {
      //     cout << "[ERROR] TnP ERROR: No matched L2 filter2\n";
      //     continue;
      //   }

      //   // decide L2, L3 muons
      //   const bool mupl_isL2 = (mupl_L2 && !mupl_L3);
      //   const bool mupl_isL3 = (mupl_L2 && mupl_L3);
      //   const bool mumi_isL2 = (mumi_L2 && !mumi_L3);
      //   const bool mumi_isL3 = (mumi_L2 && mumi_L3);

      //   bool SelDone = false;
      //   if (mupl_isL2 && mumi_isL3)
      //   {
      //     tnp_trig_w_pl = tnp_weight_trg_pbpb(mupl.Pt(), mupl.Eta(), 2, 0);
      //     tnp_trig_w_mi = tnp_weight_trg_pbpb(mumi.Pt(), mumi.Eta(), 3, 0);
      //     SelDone = true;
      //   }
      //   else if (mupl_isL3 && mumi_isL2)
      //   {
      //     tnp_trig_w_pl = tnp_weight_trg_pbpb(mupl.Pt(), mupl.Eta(), 3, 0);
      //     tnp_trig_w_mi = tnp_weight_trg_pbpb(mumi.Pt(), mumi.Eta(), 2, 0);
      //     SelDone = true;
      //   }
      //   else if (mupl_isL3 && mumi_isL3)
      //   {
      //     // if both are L3 muon, consider one of them as L2 (randomly)
      //     const int t = (std::rand() % 2 == 0) ? -1 : 1;
      //     if (t == -1)
      //     {
      //       tnp_trig_w_pl = tnp_weight_trg_pbpb(mupl.Pt(), mupl.Eta(), 2, 0);
      //       tnp_trig_w_mi = tnp_weight_trg_pbpb(mumi.Pt(), mumi.Eta(), 3, 0);
      //     }
      //     else
      //     {
      //       tnp_trig_w_pl = tnp_weight_trg_pbpb(mupl.Pt(), mupl.Eta(), 3, 0);
      //       tnp_trig_w_mi = tnp_weight_trg_pbpb(mumi.Pt(), mumi.Eta(), 2, 0);
      //     }
      //     SelDone = true;
      //   }

      //   if (!SelDone || (tnp_trig_w_pl < 0 || tnp_trig_w_mi < 0)) {
      //     cout << "[ERROR] TnP ERROR: No muon L2, L3 filter combination selected\n";
      //     continue;
      //   }
          
      //   tnp_weight *= tnp_trig_w_pl * tnp_trig_w_mi;
      // }

      // ----- fill per-candidate outputs -----
      phi[nDimu] = JP.Phi();
      phi1[nDimu] = mupl.Phi();
      phi2[nDimu] = mumi.Phi();

      if (isMC) TnPweight[nDimu] = tnp_weight;
      else TnPweight[nDimu] = 1.0;

      mass[nDimu] = JP.M();
      
      y[nDimu] = JP.Rapidity();
      
      pt[nDimu] = JP.Pt();
      pt1[nDimu] = mupl.Pt();
      pt2[nDimu] = mumi.Pt();

      eta[nDimu] = JP.Eta();
      eta1[nDimu] = mupl.Eta();
      eta2[nDimu] = mumi.Eta();

      phi[nDimu] = JP.Phi();
      phi1[nDimu] = mupl.Phi();
      phi2[nDimu] = mumi.Phi();

      recoQQsign[nDimu] = Reco_QQ_sign[irqq];

      ctau3D[nDimu] = Reco_QQ_ctau3D[irqq];
      ctau3DErr[nDimu] = Reco_QQ_ctauErr3D[irqq];
      ctau3DRes[nDimu] = (std::fabs(ctau3DErr[nDimu]) > 0.f)
                             ? (ctau3D[nDimu] / ctau3DErr[nDimu])
                             : 0.f;
      
      ++nDimu;
    } // end of dimuon loop

    if (nDimu > 0)
    {
      ++count_dimuon;
      flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== save results =====
  fFlowSkim->cd();
  flowTree->Write("myTree");
  fFlowSkim->Close();

  // ===== release memory =====
  delete JP_Reco;
  delete mupl_Reco;
  delete mumi_Reco;

  // ----- print -----
  cout << "Total saved dimuon: " << count_dimuon << "\n";

  cout << "Finish make_flowSkim_pb23()\n";
}