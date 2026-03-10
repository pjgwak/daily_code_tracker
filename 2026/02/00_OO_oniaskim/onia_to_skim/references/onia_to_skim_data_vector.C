// For PbPb2023 data

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
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
  return ((fabs(eta) < 1.0 && pt > 3.3) ||
          (1.0 < fabs(eta) && fabs(eta) < 1.3 && pt > 7.3 - 4*fabs(eta))||
          (1.3 < fabs(eta) && fabs(eta) < 1.7 && pt > 3.2525 - 1.325*fabs(eta)) ||
          (1.7 < fabs(eta) && fabs(eta) < 2.4 && pt > 1.0)
        );
}

void onia_to_skim_data_vector(bool isMC = false, long nEvt = -1)
{
  cout << "Start onia_to_skim_data()\n";

  // ===== read input =====
  //TFile *fInput = TFile::Open("oniatrees/Oniatree_DiQuarkonia_PbPb2024PromptReco_141X.root", "read");
  // TTree *oniaTree = (TTree *)fInput->Get("hionia/myTree");
  TString basePath = "/eos/cms/store/group/phys_heavyions/dileptons/Data2025/PbPb/FastOniatrees/"
//                   "Run399717_399727_DCS_only/Oniatree_PbPb2025PromptRecoData_miniAOD_Run399717_399727_DCS_only";
    "Run399717_399727_DCS/Oniatree_PbPb2025PromptRecoData_miniAOD_Run399717_399727_DCS_PD";

  TChain *oniaTree = new TChain("hionia/myTree");
  for (int i = 0; i <= 59; i++) {
      TString file = Form("%s%d.root", basePath.Data(), i);
      oniaTree->Add(file);
  }
  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 0, massHigh = 200;

  // ===== labeling =====
  // bool isPr=true
  // string mcLabel = "";
  // if (isMC==true, isPr==true) mcLabel = "PR";
  // else if (isMC==true, isPr==false) mcLabel = "NP";
  // skip

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
  UInt_t runNb;
  UInt_t eventNb;
  UInt_t LS;
  float zVtx;
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
  oniaTree->SetBranchAddress("runNb", &runNb);
  oniaTree->SetBranchAddress("eventNb", &eventNb);
  oniaTree->SetBranchAddress("LS", &LS);
  oniaTree->SetBranchAddress("zVtx", &zVtx);
  oniaTree->SetBranchAddress("Centrality", &Centrality);
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
  // if (isMC)
  // {
  //   oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
  //   oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  // }

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

  TFile *fFlowSkim = new TFile(Form("skim_PbPb2025_isMC%d_Run399717_399727_DCS.root", isMC), "recreate");

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

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if (iev % 100000 == 0)
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / max<Long64_t>(1, nTot)) << "%)\n";
      cout << "Saved dimuon: " << count_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total saved dimuon: " << count_dimuon << "\n";
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

    if (TMath::Abs(vz) > 15) continue;
    
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
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;

      // init weights
      weight = 1.;
      // if(isMC) weight = findNcoll(Centrality) * Gen_weight; // need it for PbPb23
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      recoQQsign[irqq] = Reco_QQ_sign[irqq];
      
      if (Reco_QQ_sign[irqq] != 0) continue;

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

  // ----- print -----
  cout << "Total saved dimuon: " << count_dimuon << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
