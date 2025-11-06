#define oniaEff_cxx
#include "oniaEff.h"
#include "tnp_weight.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>

#include <algorithm>
#include <iostream>



// define bins (Jpsi Raa new binning)
const double bins_cent_rap0024[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100};
const int nbins_cent_rap0024 = sizeof(bins_cent_rap0024)/sizeof(double) -1;
const double bins_cent_rap1824_pt3065[] = {0, 10, 20, 30, 40, 50, 100}; //Low pT forward
const int nbins_cent_rap1824_pt3065 = sizeof(bins_cent_rap1824_pt3065)/sizeof(double) -1; //Low pT forward
const double bins_cent_rap0006[] = {0, 10, 20, 30, 40, 50, 100};
const int nbins_cent_rap0006 = sizeof(bins_cent_rap0006)/sizeof(double) -1;
const double bins_cent_rap0612[] = {0, 10, 20, 30, 40, 50, 100};
const int nbins_cent_rap0612 = sizeof(bins_cent_rap0612)/sizeof(double) -1;
const double bins_cent_rap1218[] = {0, 10, 20, 30, 40, 50, 100};
const int nbins_cent_rap1218 = sizeof(bins_cent_rap1218)/sizeof(double) -1;
const double bins_cent_rap1824[] = {0, 10, 20, 30, 40, 50, 100};
const int nbins_cent_rap1824 = sizeof(bins_cent_rap1824)/sizeof(double) -1;
const double bins_pt_rap0024[] = {6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50};
const int nbins_pt_rap0024 = sizeof(bins_pt_rap0024)/sizeof(double) -1;
const double bins_pt_rap0006[] = {6.5, 8.5, 9.5, 11, 15, 20, 30, 50};
const int nbins_pt_rap0006 = sizeof(bins_pt_rap0006)/sizeof(double) -1;
const double bins_pt_rap0612[] = {6.5, 8.5, 9.5, 11, 15, 20, 30, 50};
const int nbins_pt_rap0612 = sizeof(bins_pt_rap0612)/sizeof(double) -1;
const double bins_pt_rap1218[] = {6.5, 7.5, 8.5, 9.5, 11, 15, 20, 30, 50};
const int nbins_pt_rap1218 = sizeof(bins_pt_rap1218)/sizeof(double) -1;
const double bins_pt_rap1824[] = {3, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 15, 20, 30, 50};
const int nbins_pt_rap1824 = sizeof(bins_pt_rap1824)/sizeof(double) -1;
const double bins_pt_cent010[] = {6.5, 7.5, 8.5, 9.5, 11, 13, 15, 20, 30, 50};
const int nbins_pt_cent010 = sizeof(bins_pt_cent010)/sizeof(double) -1;
const double bins_pt_cent1030[] = {6.5, 7.5, 8.5, 9.5, 11, 13, 15, 20, 30, 50};
const int nbins_pt_cent1030 = sizeof(bins_pt_cent1030)/sizeof(double) -1;
const double bins_pt_cent30100[] = {6.5, 7.5, 8.5, 9.5, 11, 13, 15, 20, 30, 50};
const int nbins_pt_cent30100 = sizeof(bins_pt_cent30100)/sizeof(double) -1;
const double bins_rap[] = {0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
const int nbins_rap = sizeof(bins_rap)/sizeof(double) -1;


// Bin boundaries for efficiency histograms
const double bins_4rap[] = {0, 0.6, 1.2, 1.8, 2.4};
const int nbins_4rap = sizeof(bins_4rap)/sizeof(double) -1;
const double bins_3cent[] = {0, 10, 30, 100};
const int nbins_3cent = sizeof(bins_3cent)/sizeof(double) -1;

// define bins (Double ratio binning)
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
const double massdown = 0.60;
const double massup = 0.40;

using namespace HI;
using namespace std;

vector<TObjArray*> oniaEff::ReadFileWeight(bool ispbpb, bool isprompt) {
   string wfilePbPb_prompt[] = {"weights_JPsi_PbPb_006_prompt.root","weights_JPsi_PbPb_0612_prompt.root","weights_JPsi_PbPb_1218_prompt.root","weights_JPsi_PbPb_1824_prompt.root"};
   string wfilePP_prompt[] = {"weights_JPsi_PP_006_prompt.root","weights_JPsi_PP_0612_prompt.root","weights_JPsi_PP_1218_prompt.root","weights_JPsi_PP_1824_prompt.root"};
   string wfilePbPb_nonprompt[] = {"weights_JPsi_PbPb_006_nonprompt.root","weights_JPsi_PbPb_0612_nonprompt.root","weights_JPsi_PbPb_1218_nonprompt.root","weights_JPsi_PbPb_1824_nonprompt.root"};
   string wfilePP_nonprompt[] = {"weights_JPsi_PP_006_nonprompt.root","weights_JPsi_PP_0612_nonprompt.root","weights_JPsi_PP_1218_nonprompt.root","weights_JPsi_PP_1824_nonprompt.root"};

   const int nidxf = sizeof(wfilePbPb_prompt)/sizeof(string);
   
   TFile *fweight[nidxf];
   vector<TObjArray*> objarr;
   
   cout << "ReadFileWeight: " << ispbpb << " " << isprompt << endl;
   for (int idxf=0; idxf<nidxf; idxf++) {
     if (ispbpb && isprompt)
       fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePbPb_prompt[idxf].c_str()));
     else if (ispbpb && !isprompt)
       fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePbPb_nonprompt[idxf].c_str()));
     else if (!ispbpb && isprompt)
       fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePP_prompt[idxf].c_str()));
     else if (!ispbpb && !isprompt)
       fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePP_nonprompt[idxf].c_str()));
     else
       cout <<"Cannot load files " << endl;
     
     cout << "weighting file : " << fweight[idxf]->GetName() << endl;
     TObjArray *objtmp = (TObjArray*)fweight[idxf]->Get("DataOverMC");
     TObjArray *obj = (TObjArray*)objtmp->Clone(Form("%s_copy",objtmp->GetName()));
     objarr.push_back(obj);
   }

   for (int idxf=0; idxf<nidxf; idxf++) {
     fweight[idxf]->Close();
   }

   return objarr;
}  

void oniaEff::Loop(const char* fname, bool ispbpb, bool isprompt, const int tnptype, const bool isacc)
{
//   In a ROOT session, you can do:
//      root> .L oniaEff.C
//      root> oniaEff t
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
   fChain->SetBranchStatus("Centrality",1);  
   fChain->SetBranchStatus("HLTriggers",1);  
   fChain->SetBranchStatus("Reco_QQ_*",1);  
   fChain->SetBranchStatus("Gen_QQ_*",1);  
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   cout << "LOOP: " << ispbpb << " " << isprompt << " " << isacc << endl;
   // Load pt weighting curves from external files
   vector<TObjArray *> wHistograms = ReadFileWeight(ispbpb, isprompt);
   
   TFile *f = new TFile(fname, "RECREATE");
   f->cd();

   TH1::SetDefaultSumw2();
   // define the histos
   // HIN-16-004 binning
   TH1F *hnum_centmid = new TH1F("hnum_cent_rap0016","hnum_cent_rap0016",nbins_centmid,bins_centmid);
   TH1F *hnum_centfwd = new TH1F("hnum_cent_rap1624","hnum_cent_rap1624",nbins_centfwd,bins_centfwd);
   TH1F *hnum_ptmid = new TH1F("hnum_pt_rap0016","hnum_pt_rap0016",nbins_ptmid,bins_ptmid);
   TH1F *hnum_ptfwd = new TH1F("hnum_pt_rap1624","hnum_pt_rap1624",nbins_ptfwd,bins_ptfwd);
   TH1F *hden_centmid = new TH1F("hden_cent_rap0016","hden_cent_rap0016",nbins_centmid,bins_centmid);
   TH1F *hden_centfwd = new TH1F("hden_cent_rap1624","hden_cent_rap1624",nbins_centfwd,bins_centfwd);
   TH1F *hden_ptmid = new TH1F("hden_pt_rap0016","hden_pt_rap0016",nbins_ptmid,bins_ptmid);
   TH1F *hden_ptfwd = new TH1F("hden_pt_rap1624","hden_pt_rap1624",nbins_ptfwd,bins_ptfwd);
   
   // HIN-16-025 binning
   // Eff vs centrality in 4+1 |y| regions (6.5-50 GeV/c), forward & low pT region
   TH1F *hnum_cent_rap[nbins_4rap+2];
   TH1F *hden_cent_rap[nbins_4rap+2];
   hnum_cent_rap[0] = new TH1F("hnum_cent_rap0024",";Centrality;Efficiency",nbins_cent_rap0024,bins_cent_rap0024);
   hden_cent_rap[0] = new TH1F("hden_cent_rap0024",";Centrality;Efficiency",nbins_cent_rap0024,bins_cent_rap0024);
   hnum_cent_rap[1] = new TH1F("hnum_cent_rap0006",";Centrality;Efficiency",nbins_cent_rap0006,bins_cent_rap0006);
   hden_cent_rap[1] = new TH1F("hden_cent_rap0006",";Centrality;Efficiency",nbins_cent_rap0006,bins_cent_rap0006);
   hnum_cent_rap[2] = new TH1F("hnum_cent_rap0612",";Centrality;Efficiency",nbins_cent_rap0612,bins_cent_rap0612);
   hden_cent_rap[2] = new TH1F("hden_cent_rap0612",";Centrality;Efficiency",nbins_cent_rap0612,bins_cent_rap0612);
   hnum_cent_rap[3] = new TH1F("hnum_cent_rap1218",";Centrality;Efficiency",nbins_cent_rap1218,bins_cent_rap1218);
   hden_cent_rap[3] = new TH1F("hden_cent_rap1218",";Centrality;Efficiency",nbins_cent_rap1218,bins_cent_rap1218);
   hnum_cent_rap[4] = new TH1F("hnum_cent_rap1824",";Centrality;Efficiency",nbins_cent_rap1824,bins_cent_rap1824);
   hden_cent_rap[4] = new TH1F("hden_cent_rap1824",";Centrality;Efficiency",nbins_cent_rap1824,bins_cent_rap1824);
   hnum_cent_rap[5] = new TH1F("hnum_cent_rap1824_pt3065",";Centrality;Efficiency",nbins_cent_rap1824_pt3065,bins_cent_rap1824_pt3065);
   hden_cent_rap[5] = new TH1F("hden_cent_rap1824_pt3065",";Centrality;Efficiency",nbins_cent_rap1824_pt3065,bins_cent_rap1824_pt3065);
   // Eff vs pT in 4+1 |y| regions
   TH1F *hnum_pt_rap[nbins_4rap+1];
   TH1F *hden_pt_rap[nbins_4rap+1];
   hnum_pt_rap[0] = new TH1F("hnum_pt_rap0024",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0024,bins_pt_rap0024);
   hden_pt_rap[0] = new TH1F("hden_pt_rap0024",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0024,bins_pt_rap0024);
   hnum_pt_rap[1] = new TH1F("hnum_pt_rap0006",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0006,bins_pt_rap0006);
   hden_pt_rap[1] = new TH1F("hden_pt_rap0006",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0006,bins_pt_rap0006);
   hnum_pt_rap[2] = new TH1F("hnum_pt_rap0612",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0612,bins_pt_rap0612);
   hden_pt_rap[2] = new TH1F("hden_pt_rap0612",";p_{T} [GeV/c];Efficiency",nbins_pt_rap0612,bins_pt_rap0612);
   hnum_pt_rap[3] = new TH1F("hnum_pt_rap1218",";p_{T} [GeV/c];Efficiency",nbins_pt_rap1218,bins_pt_rap1218);
   hden_pt_rap[3] = new TH1F("hden_pt_rap1218",";p_{T} [GeV/c];Efficiency",nbins_pt_rap1218,bins_pt_rap1218);
   hnum_pt_rap[4] = new TH1F("hnum_pt_rap1824",";p_{T} [GeV/c];Efficiency",nbins_pt_rap1824,bins_pt_rap1824);
   hden_pt_rap[4] = new TH1F("hden_pt_rap1824",";p_{T} [GeV/c];Efficiency",nbins_pt_rap1824,bins_pt_rap1824);
   // Eff vs pT in 3 centrality regions
   TH1F *hnum_pt_cent[nbins_3cent];
   TH1F *hden_pt_cent[nbins_3cent];
   hnum_pt_cent[0] = new TH1F("hnum_pt_cent0010",";p_{T} [GeV/c];Efficiency",nbins_pt_cent010,bins_pt_cent010);
   hden_pt_cent[0] = new TH1F("hden_pt_cent0010",";p_{T} [GeV/c];Efficiency",nbins_pt_cent010,bins_pt_cent010);
   hnum_pt_cent[1] = new TH1F("hnum_pt_cent1030",";p_{T} [GeV/c];Efficiency",nbins_pt_cent1030,bins_pt_cent1030);
   hden_pt_cent[1] = new TH1F("hden_pt_cent1030",";p_{T} [GeV/c];Efficiency",nbins_pt_cent1030,bins_pt_cent1030);
   hnum_pt_cent[2] = new TH1F("hnum_pt_cent30100",";p_{T} [GeV/c];Efficiency",nbins_pt_cent30100,bins_pt_cent30100);
   hden_pt_cent[2] = new TH1F("hden_pt_cent30100",";p_{T} [GeV/c];Efficiency",nbins_pt_cent30100,bins_pt_cent30100);
   // Eff vs rap integrated
   TH1F *hnum_rap = new TH1F("hnum_rap",";|y|;Efficiency",nbins_rap,bins_rap);
   TH1F *hden_rap = new TH1F("hden_rap",";|y|;Efficiency",nbins_rap,bins_rap);
   
   // also store more finely binned histos
   TH1F *hcentfine = new TH1F("hcentfine","hcentfine",200,0,200);
   TH1F *hnum_pt = new TH1F("hnum_pt",";p_{T};Entries",100,0,50);
   TH1F *hden_pt = new TH1F("hden_pt",";p_{T};Entries",100,0,50);
   TH2F *hnum2d = new TH2F("hnum2d","hnum2d",4,0,2.4,nbins_pt_rap1824,bins_pt_rap1824);
   TH2F *hden2d = new TH2F("hden2d","hden2d",4,0,2.4,nbins_pt_rap1824,bins_pt_rap1824);


   Long64_t nentries = fChain->GetEntries();
   cout << "nentries: " << nentries << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //if (!ispbpb && ientry>=10000000) {
      //  break;
      //}
      if (jentry%100000==0)
        cout << "jentry/ientry: " << jentry << " " << ientry << endl;

      // skip the event if Gen_QQ_size != 1 for now
      b_Gen_QQ_size->GetEntry(ientry); //read only this branch
      if (Gen_QQ_size !=1) continue;

      // apply gen acceptance cuts
      b_Gen_QQ_mupl_4mom->GetEntry(ientry); //read only this branch
      b_Gen_QQ_mumi_4mom->GetEntry(ientry); //read only this branch
      TLorentzVector *tlvgenpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(0);
      TLorentzVector *tlvgenmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(0);
      //if (!isGlobalMuonInAccept2015(tlvgenpl) || !isGlobalMuonInAccept2015(tlvgenmi)) continue;
      bool isgenok = true;
      if (!isacc && (!isGlobalMuonInAccept2015(tlvgenpl) || !isGlobalMuonInAccept2015(tlvgenmi)) ) isgenok=false;
      
      // fill denominators
      b_Gen_QQ_4mom->GetEntry(ientry); //read only this branch
      TLorentzVector *tlvgenqq = (TLorentzVector*) Gen_QQ_4mom->At(0);
      double genpt = tlvgenqq->Pt();
      double genrap = fabs(tlvgenqq->Rapidity());
      b_Centrality->GetEntry(ientry); //read only this branch

      // check that the dimuon is inside the analysis bins
      if (!isacc) {
        bool gen_inbin = false;
        if (genrap<1.8 && genpt>=6.5 && genpt<50.) gen_inbin = true;
        if (genrap>=1.8 && genrap<2.4 && genpt>=3. && genpt<50.) gen_inbin = true;
        if (!gen_inbin) isgenok=false;
      }
      
      double weight = (ispbpb && !isacc) ? fChain->GetWeight()*findNcoll(Centrality) : 1.;
      
      // Apply Data/MC pT ratio as a weight
      TH1D *curve;
      if (genrap>=0 && genrap<0.6)        curve = (TH1D*) wHistograms[0]->At(0);
      else if (genrap>=0.6 && genrap<1.2) curve = (TH1D*) wHistograms[1]->At(0);
      else if (genrap>=1.2 && genrap<1.8) curve = (TH1D*) wHistograms[2]->At(0);
      else if (genrap>=1.8 && genrap<2.4) curve = (TH1D*) wHistograms[3]->At(0);
      int ptbin = curve->FindBin(genpt);
      double ptweight = curve->GetBinContent(ptbin);
    
      weight = weight*ptweight;

      hcentfine->Fill(Centrality,weight);
      hden2d->Fill(genrap,genpt,weight);
      hden_pt->Fill(genpt,weight);
      
      if (isgenok) {
        // HIN-16-025 binning
        // Eff vs rap integrated
        if (genrap>=0 && genrap<2.4 && genpt>=6.5 && genpt<50) {
          hden_rap->Fill(genrap,weight);
        }

        // Eff vs cent in 4 |y| regions (6.5-50 GeV/c)
        if (genrap>=bins_4rap[0] && genrap<bins_4rap[nbins_4rap] && genpt>=6.5 && genpt<50) {
          hden_cent_rap[0]->Fill(Centrality/2.0,weight); //[0] is for 0-100
        }
        
        for (int i=0; i<nbins_4rap; i++) {
          if (genrap>=bins_4rap[i] && genrap<bins_4rap[i+1] && genpt>=6.5 && genpt<50) {
            hden_cent_rap[i+1]->Fill(Centrality/2.0,weight); //[0] is for 0-100
          }
        }
        
        // Eff vs cent at 1.8-2.4, 3-6.5
        if (genrap>=1.8 && genrap<2.4 && genpt>=3 && genpt<6.5) {
          hden_cent_rap[nbins_4rap+1]->Fill(Centrality/2.0,weight); //last histogram is for low pT
        }
        
        // Eff vs pt in 4 |y| regions
        if (genrap>=bins_4rap[0] && genrap<bins_4rap[nbins_4rap] && genpt>=6.5 && genpt<50) {
          hden_pt_rap[0]->Fill(genpt,weight); //[0] is for 0024
        }
        
        for (int i=0; i<nbins_4rap; i++) {
          if ( genrap>=bins_4rap[i] && genrap<bins_4rap[i+1] && genpt<50 &&
               ((genrap<1.8 && genpt>=6.5) || (genrap>1.8 && genpt>=3))
             ) {
            hden_pt_rap[i+1]->Fill(genpt,weight); //[0] is for 0024
          }
        }
        
        // Eff vs pt in 3 cent regions
        for (int i=0; i<nbins_3cent; i++) {
          if (Centrality/2.0>=bins_3cent[i] && Centrality/2.0<bins_3cent[i+1] && genrap<2.4 &&
              genpt<50 && genpt>=6.5) {
            hden_pt_cent[i]->Fill(genpt,weight);
          }
        }

        // HIN-16-004 binning
        if (genrap < 1.6 && genpt>=6.5 && genpt<30) {
          hden_centmid->Fill(Centrality/2.0,weight);
          hden_ptmid->Fill(genpt,weight);
        } else if (genrap >= 1.6 && genrap < 2.4 && genpt>=3 && genpt<30) {
          hden_centfwd->Fill(Centrality/2.0,weight);
          hden_ptfwd->Fill(genpt,weight);
        }
      } // end of if(isgenok)

      ///// start gen after single muon selection - acceptance numerator 
      if (isacc) {
        if (isGlobalMuonInAccept2015(tlvgenpl) && isGlobalMuonInAccept2015(tlvgenmi)) {
          if (genrap>=0 && genrap<2.4 && genpt>=6.5 && genpt<50) {
            hnum_rap->Fill(genrap,weight);
          }

          // Eff vs cent in 4 |y| regions (6.5-50 GeV/c)
          if (genrap>=bins_4rap[0] && genrap<bins_4rap[nbins_4rap] && genpt>=6.5 && genpt<50) {
            hnum_cent_rap[0]->Fill(Centrality/2.0,weight); //[0] is for 0-100
          }
          
          for (int i=0; i<nbins_4rap; i++) {
            if (genrap>=bins_4rap[i] && genrap<bins_4rap[i+1] && genpt>=6.5 && genpt<50) {
              hnum_cent_rap[i+1]->Fill(Centrality/2.0,weight); //[0] is for 0-100
            }
          }
          // Eff vs cent at 1.8-2.4, 3-6.5
          if (genrap>=1.8 && genrap<2.4 && genpt>=3 && genpt<6.5) {
              hnum_cent_rap[nbins_4rap+1]->Fill(Centrality/2.0,weight); //last histogram is for low pT
          }
          
          // Eff vs pt in 4 |y| regions
          if (genrap>=bins_4rap[0] && genrap<bins_4rap[nbins_4rap] && genpt>=6.5 && genpt<50) {
            hnum_pt_rap[0]->Fill(genpt,weight); //[0] is for 0024
          }
          
          for (int i=0; i<nbins_4rap; i++) {
            if ( genrap>=bins_4rap[i] && genrap<bins_4rap[i+1] && genpt<50 &&
                 ((genrap<1.8 && genpt>=6.5) || (genrap>1.8 && genpt>=3))
               ) {
              hnum_pt_rap[i+1]->Fill(genpt,weight); //[0] is for 0024
            }
          }
          
          // Eff vs pt in 3 cent regions
          for (int i=0; i<nbins_3cent; i++) {
            if (Centrality/2.0>=bins_3cent[i] && Centrality/2.0<bins_3cent[i+1] && genrap<2.4 &&
                genpt<50 && genpt>=6.5) {
              hnum_pt_cent[i]->Fill(genpt,weight);
            }
          }
          
          // HIN-16-004 binning
          if (genrap < 1.6 && genpt>=6.5 && genpt<30) {
            hnum_centmid->Fill(Centrality/2.0,weight);
            hnum_ptmid->Fill(genpt,weight);
          } else if (genrap >= 1.6 && genrap < 2.4 && genpt>=3 && genpt<30) {
            hnum_centfwd->Fill(Centrality/2.0,weight);
            hnum_ptfwd->Fill(genpt,weight);
          }

          hnum2d->Fill(genrap,genpt,weight);
        }

      } else { ///// start reco loop - efficiency numerator
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
           
          // gen-reco matching
          TLorentzVector *tlvrecpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
          TLorentzVector *tlvrecmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
          double dr = max(tlvrecpl->DeltaR(*tlvgenpl),tlvrecmi->DeltaR(*tlvgenmi));
          if (dr<mindr) {
            mindr = dr;
            ibestqq = i;
          }

        } // Reco_QQ loop
        
        if (ibestqq<0) continue;
        
        TLorentzVector *tlvrecqq = (TLorentzVector*) Reco_QQ_4mom->At(ibestqq);
        double recopt = tlvrecqq->Pt();
        double recorap = fabs(tlvrecqq->Rapidity());

        // check that the dimuon is inside the analysis bins
        bool rec_inbin = false;
        if (recorap<1.8 && recopt>=6.5 && recopt<50.) rec_inbin = true;
        if (recorap>=1.8 && recorap<2.4 && recopt>=3. && recopt<50.) rec_inbin = true;
        if (!rec_inbin) continue;
         
        TLorentzVector *tlvrecpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(ibestqq);
        TLorentzVector *tlvrecmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(ibestqq);
        double recMuPlpt = tlvrecpl->Pt();
        double recMuPlEta = tlvrecpl->Eta();
        double recMuMipt = tlvrecmi->Pt();
        double recMuMiEta = tlvrecmi->Eta();
       
        double tnp_weight = 1.0;
        if (ispbpb) { 
          if (tnptype == trg) { // nominal
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg_binned) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,-10) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,-10) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg_plus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,-1) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg_minus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,-2) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__muid__sta) {
            tnp_weight = tnp_weight_sta_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_muid_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__muid) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_muid_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__muid_plus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pbpb(recMuPlpt,recMuPlEta,-1) * tnp_weight_muid_pbpb(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__muid_minus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pbpb(recMuPlpt,recMuPlEta,-2) * tnp_weight_muid_pbpb(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__sta) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__sta_plus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,-1) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg__sta_minus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,-2) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pbpb(0) * tnp_weight_trk_pbpb(0);
          } else if (tnptype == trg_trk_plus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(-1) * tnp_weight_trk_pbpb(-1);
          } else if (tnptype == trg_trk_minus1sigma) {
            tnp_weight = tnp_weight_trg_pbpb(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pbpb(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pbpb(-2) * tnp_weight_trk_pbpb(-2);
          }
        } else { //pp
          if (tnptype == trg) { // nominal
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg_binned) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,-10) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,-10) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg_plus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,-1) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg_minus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,-2) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__muid__sta) {
            tnp_weight = tnp_weight_sta_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_muid_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__muid) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_muid_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__muid_plus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pp(recMuPlpt,recMuPlEta,-1) * tnp_weight_muid_pp(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__muid_minus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_muid_pp(recMuPlpt,recMuPlEta,-2) * tnp_weight_muid_pp(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__sta) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__sta_plus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,-1) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,-1) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg__sta_minus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_sta_pp(recMuPlpt,recMuPlEta,-2) * tnp_weight_sta_pp(recMuMipt,recMuMiEta,-2) *
                         tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
          } else if (tnptype == trg_trk_plus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(-1) * tnp_weight_trk_pp(-1);
          } else if (tnptype == trg_trk_minus1sigma) {
            tnp_weight = tnp_weight_trg_pp(recMuPlpt,recMuPlEta,0) * tnp_weight_trg_pp(recMuMipt,recMuMiEta,0) *
                         tnp_weight_trk_pp(-2) * tnp_weight_trk_pp(-2);
          }
        }

        
        if (tnptype != noTnPSFs) {
          // tnp scale factors applied only when it's requested
          weight = weight*tnp_weight;
        }
        
        // fill the numerators
        // Eff vs rap integrated
        if (recorap>=0 && recorap<2.4 && recopt>=6.5 && recopt<50)
          hnum_rap->Fill(genrap,weight);

        // Eff vs cent in 4 |y| regions (6.5-50 GeV/c)
        if (recorap>=bins_4rap[0] && recorap<bins_4rap[nbins_4rap] && recopt>=6.5 && recopt<50)
          hnum_cent_rap[0]->Fill(Centrality/2.0,weight); //[0] is for 0-100
        
        for (int i=0; i<nbins_4rap; i++) {
          if (recorap>=bins_4rap[i] && recorap<bins_4rap[i+1] && recopt>=6.5 && recopt<50)
            hnum_cent_rap[i+1]->Fill(Centrality/2.0,weight); //[0] is for 0-100
        }
        // Eff vs cent at 1.8-2.4, 3-6.5
        if (recorap>=1.8 && recorap<2.4 && recopt>=3 && recopt<6.5)
            hnum_cent_rap[nbins_4rap+1]->Fill(Centrality/2.0,weight); //last histogram is for low pT
        
        // Eff vs pt in 4 |y| regions
        if (recorap>=bins_4rap[0] && recorap<bins_4rap[nbins_4rap] && recopt>=6.5 && recopt<50)
          hnum_pt_rap[0]->Fill(genpt,weight); //[0] is for 0024
        
        for (int i=0; i<nbins_4rap; i++) {
          if ( recorap>=bins_4rap[i] && recorap<bins_4rap[i+1] && recopt<50 &&
               ((recorap<1.8 && recopt>=6.5) || (recorap>1.8 && recopt>=3))
             ) {
            hnum_pt_rap[i+1]->Fill(genpt,weight); //[0] is for 0024
          }
        }
        
        // Eff vs pt in 3 cent regions
        for (int i=0; i<nbins_3cent; i++) {
          if (Centrality/2.0>=bins_3cent[i] && Centrality/2.0<bins_3cent[i+1] && recorap<2.4 &&
              recopt<50 && recopt>=6.5) {
            hnum_pt_cent[i]->Fill(genpt,weight);
          }
        }
        
        // HIN-16-004 binning
        if (recorap < 1.6 && recopt>=6.5 && recopt<30) {
          hnum_centmid->Fill(Centrality/2.0,weight);
          hnum_ptmid->Fill(genpt,weight);
        } else if (recorap >= 1.6 && recorap < 2.4 && recopt>=3 && recopt<30) {
          hnum_centfwd->Fill(Centrality/2.0,weight);
          hnum_ptfwd->Fill(genpt,weight);
        }

        hnum2d->Fill(recorap,recopt,weight);
        hnum_pt->Fill(recopt,weight);
      } // end of isacc condition - fill up numerators   
   } // event loop
  

   f->Write();
   f->Close();

}
