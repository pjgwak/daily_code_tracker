#include <iostream>
// #include "Style.h"
// #include "tdrstyle.C"
// #include "rootFitHeaders.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "TLegend.h"
#include "TString.h"
#include "TStopwatch.h"

using namespace RooFit;


bool isMuonAcc2023(double pt, double eta);
TVector3 MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);
TVector3 MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab);

void get_eff_file()
{
  // start timer
  TStopwatch *t = new TStopwatch;
  t -> Start();

  // set labels
  string ptWgtLabel = "No";
  string userLabel = "PbPb2023_ExplicitHLT";

  // === read tree->input ===
  TFile *rf = new TFile("/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Prompt_with_Track/Oniatree_MC_miniAOD_PbPb23_Prompt_with_Track.root");
  TTree *tree = (TTree*) rf->Get("hionia/myTree");
  // tree->Print();

  // [Not] === read pT weight file ===
  // need to produce it -> Need NP separation first
  // TFile *fPtW1 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_y0_1p6_211201.root","read");
  // TFile *fPtW2 = new TFile("./roots/ratioDataMC_AA_Jpsi_DATA_Forward_y_211218.root","read");
  // TF1* fptw1 = (TF1*) fPtW1->Get("dataMC_Ratio1");
  // TF1* fptw2 = (TF1*) fPtW2->Get("dataMC_Ratio1");

  // === connect to branch ===
  // GEN branches
  Short_t Gen_QQ_size;
  Short_t Gen_QQ_mupl_idx[1000];
  Short_t Gen_QQ_mumi_idx[1000];

  Float_t Gen_QQ_ctau[1000];
  Float_t Gen_QQ_ctau3D[1000];

  TClonesArray *Gen_mu_4mom  = nullptr;
  TClonesArray *Gen_QQ_4mom  = nullptr;

  // === Connect branches (no branch pointer needed) ===
  tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size);
  tree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
  tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);

  tree->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau);
  tree->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D);

  tree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
  tree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);

  // RECO branches
  Short_t Reco_QQ_size;
  Short_t Reco_mu_size;
  Short_t Reco_QQ_mupl_idx[1000];
  Short_t Reco_QQ_mumi_idx[1000];

  Float_t Reco_QQ_ctau[1000];
  Float_t Reco_QQ_ctau3D[1000];

  TClonesArray *Reco_mu_4mom = nullptr;
  TClonesArray *Reco_QQ_4mom = nullptr;

  // --------- Added ------------
  Float_t Reco_QQ_VtxProb[1000];
  Short_t Reco_QQ_sign[1000];
  Short_t Reco_mu_whichGen[1000]; // MC

  Bool_t Reco_mu_isSoftCutBased[1000];
  Bool_t Reco_mu_isGlobal[1000];
  Bool_t Reco_mu_isTracker[1000];

  ULong64_t HLTriggers;
  // ----------------------------

  // === Connect branches (no branch pointer needed) ===
  tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  tree->SetBranchAddress("Reco_mu_size", &Reco_mu_size);
  tree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
  tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);

  tree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau);
  tree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);

  tree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  tree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);

  tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  tree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  tree->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased);
  tree->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal);
  tree->SetBranchAddress("Reco_mu_isTracker", Reco_mu_isTracker);

  tree->SetBranchAddress("HLTriggers", &HLTriggers);

  // === preapre output file and trees ===
  TFile *outfile = new TFile(Form("roots/efficiency_Prompt_jpsi_ptWgt%s%s.root", ptWgtLabel.c_str(), userLabel.c_str()),"RECREATE");
  TTree *tDen = new TTree("tDen", "Denominator dimuons");
  TTree *tNum = new TTree("tNum", "Numerator dimuons");

  // connect branch
  Float_t den_pt = 0.0, den_y = 0.0;
  Float_t den_cosHX = 0.0, den_cosCS = 0.0;
  Float_t den_phiHX = 0.0, den_phiCS = 0.0;
  Float_t num_pt = 0.0, num_y = 0.0;
  Float_t num_cosHX = 0.0, num_cosCS = 0.0;
  Float_t num_phiHX = 0.0, num_phiCS = 0.0;

  tDen->Branch("pt", &den_pt, "pt/F");
  tDen->Branch("y", &den_y, "y/F");
  tDen->Branch("cosHX", &den_cosHX, "cosHX/F");
  tDen->Branch("cosCS", &den_cosCS, "cosCS/F");
  tDen->Branch("phiHX", &den_phiHX, "phiHX/F");
  tDen->Branch("phiCS", &den_phiCS, "phiCS/F");

  tNum->Branch("pt", &num_pt, "pt/F");
  tNum->Branch("y", &num_y, "y/F");
  tNum->Branch("cosHX", &num_cosHX, "cosHX/F");
  tNum->Branch("cosCS", &num_cosCS, "cosCS/F");
  tNum->Branch("phiHX", &num_phiHX, "phiHX/F");
  tNum->Branch("phiCS", &num_phiCS, "phiCS/F");

  // === main loop start ===
  // prepare 4-vector variables
  TLorentzVector* JP= new TLorentzVector;
  TLorentzVector* Mu1= new TLorentzVector;
  TLorentzVector* Mu2= new TLorentzVector;

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;

  double pt = 0.0, y = 0.0, mu1_pt = 0.0, mu1_eta = 0.0;
  double mu2_pt = 0.0, mu2_eta = 0.0;
  double mass = 0.0, muPtCut = 0.0;

  // event loop
  for (int i = 0; i < nEvt /*nEvt*/; i++)
  {
    tree->GetEntry(i);
    if (i%100000==0) cout << ">>>>> EVENT " << i << " / " << tree->GetEntries() <<  endl;

    //cout<<"# of Gen QQ : "<<Gen_QQ_size<<endl;
    //if(Gen_QQ_size > 1) continue;
    
    // GEN dimuon loop
    for(int j = 0; j < Gen_QQ_size; j++) {
      // bring 4-vectors
      JP = (TLorentzVector*) Gen_QQ_4mom->At(j);
      Mu1 = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[j]);
      Mu2 = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[j]);
      
      // bring kinematics
      pt = JP->Pt();
      mass = JP->M();
      y = JP->Rapidity();
      mu1_pt = Mu1->Pt();
      mu1_eta = Mu1->Eta();
      mu2_pt = Mu2->Pt();
      mu2_eta = Mu2->Eta();

      // === frame transformation ===
      // HX, CS - pp doen't make EP -> dimuon level
      // HX
      TVector3 muPlus_HX = MuPlusVector_Helicity(*JP, *Mu1);
			float cosHXVar = muPlus_HX.CosTheta();
			float phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

      // CS
      TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*JP, *Mu1);
      float cosCSVar = muPlus_CS.CosTheta();
			float phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

      // === denominator cut ===
      // dimuon satisfying |y| < 2.4, each muons pass single muon acceptance cut
      if (fabs(y) > 2.4) continue;
      bool mu1pass = isMuonAcc2023(mu1_pt, mu1_eta);
      bool mu2pass = isMuonAcc2023(mu2_pt, mu2_eta);
      if (!(mu1pass && mu2pass)) continue;

      den_pt = pt;
      den_y = y;
      den_cosHX = cosHXVar;
      den_phiHX = phiHXVar;
      den_cosCS = cosCSVar;
      den_phiCS = phiCSVar;
      tDen->Fill();
    }

    // RECO dimuon loop
    for (int j = 0; j < Reco_QQ_size; j++)
    {
      // bring 4-vectors
      JP = (TLorentzVector *)Reco_QQ_4mom->At(j);
      Mu1 = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mupl_idx[j]);
      Mu2 = (TLorentzVector *)Reco_mu_4mom->At(Reco_QQ_mumi_idx[j]);

      // bring kinematics
      pt = JP->Pt();
      mass = JP->M();
      y = JP->Rapidity();
      mu1_pt = Mu1->Pt();
      mu1_eta = Mu1->Eta();
      mu2_pt = Mu2->Pt();
      mu2_eta = Mu2->Eta();

      // === frame transformation ===
      // HX, CS - pp doen't make EP -> dimuon level
      // HX
      TVector3 muPlus_HX = MuPlusVector_Helicity(*JP, *Mu1);
      float cosHXVar = muPlus_HX.CosTheta();
      float phiHXVar = muPlus_HX.Phi() * 180 / TMath::Pi();

      // if (cosHXVar < 0) {
      // 	// if phi value is smaller than -pi, add 2pi
      // 	if ((phiHXVar - 135) < -180)
      // 		phiHXVar += 225;
      // 	else
      // 		phiHXVar -= 135;
      // }

      // else if (cosHXVar > 0) {
      // 	// if phi value is smaller than -pi, add 2pi
      // 	if ((phiHXVar - 45) < -180)
      // 		phiHXVar += 315;
      // 	else
      // 		phiHXVar -= 45;
      // }

      // CS
      TVector3 muPlus_CS = MuPlusVector_CollinsSoper(*JP, *Mu1);
      float cosCSVar = muPlus_CS.CosTheta();
      float phiCSVar = muPlus_CS.Phi() * 180 / TMath::Pi();

      // if (cosCSVar < 0) {
      // 	// if phi value is smaller than -pi, add 2pi
      // 	if ((phiCSVar - 135) < -180)
      // 		phiCSVar += 225;
      // 	else
      // 		phiCSVar -= 135;
      // }

      // else if (cosCSVar > 0) {
      // 	// if phi value is smaller than -pi, add 2pi
      // 	if ((phiCSVar - 45) < -180)
      // 		phiCSVar += 315;
      // 	else
      // 		phiCSVar -= 45;
      // }

      // === numerator cut ===
      // dimuon satisfying |y| < 2.4, each muons pass single muon acceptance cut
      // + full reconstruction cuts
      if (fabs(y) > 2.4) continue;
      bool mu1pass = isMuonAcc2023(mu1_pt, mu1_eta);
      bool mu2pass = isMuonAcc2023(mu2_pt, mu2_eta);
      if (!(mu1pass && mu2pass)) continue;

      if (Reco_QQ_VtxProb[j] < 0.01f) continue;
      if (Reco_QQ_sign[j] != 0) continue;
      if (Reco_mu_whichGen[Reco_QQ_mupl_idx[j]] == -1) continue;
      if (Reco_mu_whichGen[Reco_QQ_mumi_idx[j]] == -1) continue;

      if (!Reco_mu_isSoftCutBased[Reco_QQ_mupl_idx[j]] || !Reco_mu_isSoftCutBased[Reco_QQ_mumi_idx[j]]) continue;
      if (!Reco_mu_isGlobal[Reco_QQ_mupl_idx[j]] || !Reco_mu_isGlobal[Reco_QQ_mumi_idx[j]]) continue;
      if (!Reco_mu_isTracker[Reco_QQ_mupl_idx[j]] || !Reco_mu_isTracker[Reco_QQ_mumi_idx[j]]) continue;

      // apply HLT - for PbPb2023 MC
      bool HLTPass = false;
      if (
        (HLTriggers & ((ULong64_t)pow(2, 24))) == ((ULong64_t)pow(2, 24)) ||
        (HLTriggers & ((ULong64_t)pow(2, 25))) == ((ULong64_t)pow(2, 25)) ||
        (HLTriggers & ((ULong64_t)pow(2, 26))) == ((ULong64_t)pow(2, 26))) HLTPass = true;

      // if((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))) HLTPass=true;
      if (HLTPass == false) continue;

      // === TnP Here Later ===
      // === End of TnP ===

      num_pt = pt;
      num_y = y;
      num_cosHX = cosHXVar;
      num_phiHX = phiHXVar;
      num_cosCS = cosCSVar;
      num_phiCS = phiCSVar;
      tNum->Fill();

      // bool mu1pass = isMuonAcc2023(mu1_pt,mu1_eta);
      // bool mu2pass = isMuonAcc2023(mu2_pt,mu2_eta);
      // if (mu1pass!=true || mu2pass!=true) continue;

      // [Not] === apply pT weight - numerator ===
      // if (1.6<=fabs(y) && fabs(y)<2.4){
      //   hNumPt_2021_Fory->Fill(pt,wt2);
      // }
      // if (1.6<=fabs(y) && fabs(y)<2.4 && pt > 3.5 && pt < 50) {
      //   hNumPt_2021_Fory_Int->Fill(1,wt2); hNumY_2021->Fill(y,wt2);
      // }
      // if (fabs(y)<1.6) {
      //   hNumPt_2021_midy->Fill(pt,wt1);
      // }
      // if (fabs(y)<1.6 && pt >6.5 && pt < 50) {
      //   hNumPt_2021_midy_Int->Fill(1,wt1); hNumY_2021->Fill(y,wt1);
      // }
    }
  }

  // Tree fill check
  // tDen->Print();
  // tNum->Print();

  // === save output ===
  outfile->cd();
  tDen->Write();
  tNum->Write();
  outfile->Close();

  // print run time
  t -> Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",t->RealTime(),t->CpuTime());
}

// [not] === Acc calculation ===
// Use other file for calculation


// 2023 soft single muon acceptance cut
bool isMuonAcc2023(double pt, double eta)
{
  return ((fabs(eta) < 1.0 && pt > 3.3) ||
          (1.0 < fabs(eta) && fabs(eta) < 1.3 && pt > 7.3 - 4*fabs(eta))||
          (1.3 < fabs(eta) && fabs(eta) < 1.7 && pt > 3.2525 - 1.325*fabs(eta)) ||
          (1.7 < fabs(eta) && fabs(eta) < 2.4 && pt > 1.0)
        );
}

TVector3 MuPlusVector_Helicity(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
	// ******** Transform variables of muons from the lab frame to the upsilon's rest frame ******** //
	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	TLorentzVector MuPlusLV_QQRestFrame(MuPlusLV_Lab);

	//(Note. TLorentzVector.BoostVector() gives beta(=px/E,py/E,pz/E) of the parents)
	//(TLorentzVector.Boost() boosts from the rod frame to the lab frame, so plug in -beta to get lab to rod)
	MuPlusLV_QQRestFrame.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 MuPlusVec_Boosted = MuPlusLV_QQRestFrame.Vect();

	//Note: TVector3.Rotate() rotates the vectors, not the coordinates, so should rotate -phi and -theta

	MuPlusVec_Boosted.RotateZ(-QQVector_Lab.Phi());

	MuPlusVec_Boosted.RotateY(-QQVector_Lab.Theta());

	return MuPlusVec_Boosted;
}

// Lab to Collins-Soper
// requires the beam parameters
TVector3 MuPlusVector_CollinsSoper(const TLorentzVector &QQLV_Lab, const TLorentzVector &MuPlusLV_Lab) {
	// ******** Set beam energy for the Collins-Soper reference frame ******** //
	double sqrt_S_NN = 5.32;                    //(Center of mass Energy per nucleon pair in TeV)
	double beamEnergy = sqrt_S_NN * 1000. / 2.; //(in GeV) (Note. sqrt_S_NN = sqrt(2*E1*E2+2*p1*p2) = 2E1 when two beams have the same E)

	// ******** HX to CS (rotation from HX frame to CS frame) ******** //
	// (1. Boost two beams to upsilon's rest frame)
	// (2. Rotate the coordinates)
	// (3. Get angle between two beams(b1 and -b2), and between b1 and ZHX in the upsilon's rest frame)
	// (4. Calculate delta (angle btw ZHX and ZCS))

	// ******** Transform variables of beams from the lab frame to the upsilon's rest frame ******** //
	TLorentzVector Beam1LV_Boosted(0., 0., beamEnergy, beamEnergy);
	TLorentzVector Beam2LV_Boosted(0., 0., -beamEnergy, beamEnergy); // mind the sign!!

	Beam1LV_Boosted.Boost(-QQLV_Lab.BoostVector());
	Beam2LV_Boosted.Boost(-QQLV_Lab.BoostVector());

	// ******** Rotate the coordinates ******** //
	TVector3 Beam1Vector_QQRestFrame(Beam1LV_Boosted.Vect());
	TVector3 Beam2Vector_QQRestFrame(Beam2LV_Boosted.Vect());

	TVector3 QQVector_Lab = QQLV_Lab.Vect();
	Beam1Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam1Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	Beam2Vector_QQRestFrame.RotateZ(-QQVector_Lab.Phi());
	Beam2Vector_QQRestFrame.RotateY(-QQVector_Lab.Theta());

	// ******** Calculate the angle between z_HX and z_CS ******** //
	TVector3 ZHXunitVec(0, 0, 1.); //(define z_HX unit vector)
	double Angle_B1ZHX = Beam1Vector_QQRestFrame.Angle(ZHXunitVec); //(angle between beam1 and z_HX)
	double Angle_B2ZHX = Beam2Vector_QQRestFrame.Angle(-ZHXunitVec); //(angle between beam2 and -z_HX =(-beam2 and z_HX) )
	double Angle_B1miB2 = Beam1Vector_QQRestFrame.Angle(-Beam2Vector_QQRestFrame); //(angle between beam1 and -beam2)

	double delta = 0; //(define and initialize the angle between z_HX and z_CS)

	// Maths for calculating the angle between z_HX and z_CS is different depending on the sign of the beam1's z-coordinate)
	if (Angle_B1ZHX > Angle_B2ZHX)
		delta = Angle_B2ZHX + Angle_B1miB2 / 2.;
	else if (Angle_B1ZHX < Angle_B2ZHX)
		delta = Angle_B1ZHX + Angle_B1miB2 / 2.;
	else
		std::cout << "beam1PvecBoosted.Pz() = 0?" << std::endl;

	// ******** Rotate the coordinates along the y-axis by the angle between z_HX and z_CS ******** //
	TVector3 MuPlusVec_CS(MuPlusVector_Helicity(QQLV_Lab, MuPlusLV_Lab));

	MuPlusVec_CS.RotateY(delta);

	return MuPlusVec_CS;
}