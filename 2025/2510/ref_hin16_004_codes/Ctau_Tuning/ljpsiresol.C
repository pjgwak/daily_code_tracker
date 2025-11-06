#define ljpsiresol_cxx
#include "ljpsiresol.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"

#include <iostream>

// define bins
const int nbins_centmid = 6;
const float bins_centmid[nbins_centmid+1] = {0, 10, 20, 30, 40, 50, 100};
const int nbins_centfwd = 3;
const float bins_centfwd[nbins_centfwd+1] = {0, 20, 40, 100};
const int nbins_ptmid = 5;
const float bins_ptmid[nbins_ptmid+1] = {6.5, 9, 12, 15, 20, 30};
const int nbins_ptfwd = 3;
const float bins_ptfwd[nbins_ptfwd+1] = {3, 6.5, 12, 30};


// other settings
const double maxdr = 0.03;
const double massjpsi = 3.096;
const double masspsip = 3.686;
const double massdown = 0.15;
const double massup = 0.10;

using namespace HI;
using namespace std;

void ljpsiresol::Loop(const char* fname, bool isdata)
{
//   In a ROOT session, you can do:
//      root> .L ljpsiresol.C
//      root> ljpsiresol t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("HLTriggers",1);  
   fChain->SetBranchStatus("Reco_QQ_*",1);  
   fChain->SetBranchStatus("Gen_QQ_*",!isdata);  
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *f = new TFile(fname, "RECREATE");
   f->cd();

   // define the histos
   // we need 2 numerators (with and without ctau3d cut), 1 denominator. Let's make them in arrays
   // TH1F *hnum_centmid = new TH1F("hnum_centmid","hnum_centmid",nbins_centmid,bins_centmid);
   vector<TH1F*> hists_mid, hists_fwd;
   for (int i=0; i<nbins_ptmid; i++) {
      double ptmin = bins_ptmid[i];
      double ptmax = bins_ptmid[i+1];
      TH1F *hist = new TH1F(Form("hmid_pt%.1f-%.1f",ptmin,ptmax),"",70,-0.15,0.2);
      hists_mid.push_back(hist);
   }
   for (int i=0; i<nbins_ptfwd; i++) {
      double ptmin = bins_ptfwd[i];
      double ptmax = bins_ptfwd[i+1];
      TH1F *hist = new TH1F(Form("hfwd_pt%.1f-%.1f",ptmin,ptmax),"",70,-0.15,0.2);
      hists_fwd.push_back(hist);
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // if (ientry>1000) break;

      // skip the event if Gen_QQ_size != 1 for now
      TLorentzVector *tlvgenpl=NULL, *tlvgenmi=NULL, *tlvgenqq=NULL;
      if (!isdata) {
         b_Gen_QQ_size->GetEntry(ientry); //read only this branch
         if (Gen_QQ_size !=1) continue;

         // apply gen acceptance cuts
         b_Gen_QQ_mupl_4mom->GetEntry(ientry); //read only this branch
         b_Gen_QQ_mumi_4mom->GetEntry(ientry); //read only this branch
         tlvgenpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(0);
         tlvgenmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(0);
         if (!isGlobalMuonInAccept2015(tlvgenpl) || !isGlobalMuonInAccept2015(tlvgenmi)) continue;

         // fill denominators
         b_Gen_QQ_4mom->GetEntry(ientry); //read only this branch
         TLorentzVector *tlvgenqq = (TLorentzVector*) Gen_QQ_4mom->At(0);
         double genpt = tlvgenqq->Pt();
         b_Centrality->GetEntry(ientry); //read only this branch

         // check that the dimuon is inside the analysis bins
         bool gen_inbin = false;
         if (fabs(tlvgenqq->Rapidity())<1.6 && tlvgenqq->Pt()>6.5 && tlvgenqq->Pt()<30.) gen_inbin = true;
         if (fabs(tlvgenqq->Rapidity())>1.6 && fabs(tlvgenqq->Rapidity())<2.4 && tlvgenqq->Pt()>3. && tlvgenqq->Pt()<30.) gen_inbin = true;
         if (!gen_inbin) continue;
      }

      double weight = 1.;

      // is there at least one reco dimuon?
      b_Reco_QQ_size->GetEntry(ientry);
      if (Reco_QQ_size==0) continue;

      // ok, now read all other branches
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // loop on the reconstructed dimuons to find the one matched to the gen one
      double mindr=999.;
      int ibestqq=-1;
      for (int i=0; i<Reco_QQ_size; i++) {
         // acceptance
         if (!areMuonsInAcceptance2015(i)) continue;
         // sign
         if (Reco_QQ_sign[i]!=0) continue;
         // quality cuts
         if (!passQualityCuts2015(i)) continue;
         // trigger matching
         if (!isTriggerMatch(i,0)) continue; // HLT_HIL1DoubleMu0_v1
         // mass cut
         double mass = ((TLorentzVector*) Reco_QQ_4mom->At(i))->M();
         double mass0 = massjpsi;
         if (mass<(mass0-massdown) || mass>(mass0+massup)) continue;
         // check that the dimuon is inside the analysis bins
         bool rec_inbin = false;
         TLorentzVector *tlvrecqq = (TLorentzVector*) Reco_QQ_4mom->At(i);
         if (fabs(tlvrecqq->Rapidity())<1.6 && tlvrecqq->Pt()>6.5 && tlvrecqq->Pt()<30.) rec_inbin = true;
         if (fabs(tlvrecqq->Rapidity())>1.6 && fabs(tlvrecqq->Rapidity())<2.4 && tlvrecqq->Pt()>3. && tlvrecqq->Pt()<30.) rec_inbin = true;
         if (!rec_inbin) continue;
         // gen-reco matching
         TLorentzVector *tlvrecpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
         TLorentzVector *tlvrecmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
         if (!isdata) {
            double dr = max(tlvrecpl->DeltaR(*tlvgenpl),tlvrecmi->DeltaR(*tlvgenmi));
            if (dr<mindr) {
               mindr = dr;
               ibestqq = i;
            }
         } else {
            ibestqq=0;
         }
      } // Reco_QQ loop

      if (ibestqq<0) continue;

      // fill the numerators
      TLorentzVector *tlvrecqq = (TLorentzVector*) Reco_QQ_4mom->At(ibestqq);
      double rap = fabs(tlvrecqq->Rapidity());
      double pt = tlvrecqq->Pt();

      if (rap<1.6) {
         for (int i=0; i<nbins_ptmid; i++)
            if (pt>=bins_ptmid[i]&&pt<bins_ptmid[i+1])
               hists_mid[i]->Fill(Reco_QQ_ctau3D[ibestqq]);
      } else if (rap<2.4) {
         for (int i=0; i<nbins_ptfwd; i++)
            if (pt>=bins_ptfwd[i]&&pt<bins_ptfwd[i+1])
               hists_fwd[i]->Fill(Reco_QQ_ctau3D[ibestqq]);
      }
   } // event loop

   f->Write();
   f->Close();
}
