#define drawPlot_cxx
#include "drawPlot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>

void drawPlot::Loop()
{
//   In a ROOT session, you can do:
//      root> .L drawPlot.C
//      root> drawPlot t
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
   fChain->SetBranchStatus("Reco_QQ_*",1);
   fChain->SetBranchStatus("HLTriggers",1);
   fChain->SetBranchStatus("Centrality",1);
   fChain->SetBranchStatus("SumET_HF",1);
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *f = new TFile("histos.root","RECREATE");
   f->cd();
   TH1F *hVprobSig = new TH1F("hVprobSig","Vertex probability distribution for J/#psi",100,0,1);
   TH1F *hVprobBkg = new TH1F("hVprobBkg","Vertex probability distribution for background",100,0,1);
   TH1F *hCentrality = new TH1F("hCentrality","Centrality",200,0,100);
   TH1F *hCentrality_mid = new TH1F("hCentrality_mid","Centrality",200,0,100);
   TH1F *hCentrality_fwd = new TH1F("hCentrality_fwd","Centrality",200,0,100);
   const int nbins = 200;
   const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
   double *Ncollbins = new double[nbins+1]; Ncollbins[0]=0; for (int i=1; i<nbins+1; i++) Ncollbins[i] = Ncoll[nbins-i]+0.01;
   const float Npart[nbins] = {401.99, 398.783, 396.936, 392.71, 387.901, 383.593, 377.914, 374.546, 367.507, 361.252, 356.05, 352.43, 345.701, 341.584, 335.148, 330.581, 325.135, 320.777, 315.074, 310.679, 306.687, 301.189, 296.769, 291.795, 287.516, 283.163, 277.818, 274.293, 269.29, 265.911, 260.574, 256.586, 252.732, 249.194, 245.011, 241.292, 236.715, 232.55, 229.322, 225.328, 221.263, 218.604, 214.728, 210.554, 206.878, 203.924, 200.84, 196.572, 193.288, 189.969, 186.894, 183.232, 180.24, 177.36, 174.008, 171.222, 168.296, 165.319, 162.013, 158.495, 156.05, 154.218, 150.559, 148.455, 145.471, 142.496, 139.715, 137.395, 134.469, 131.926, 129.817, 127.045, 124.467, 122.427, 119.698, 117.607, 114.543, 112.662, 110.696, 108.294, 105.777, 103.544, 101.736, 99.943, 97.4951, 95.4291, 93.2148, 91.2133, 89.5108, 87.2103, 85.7498, 83.5134, 81.9687, 79.7456, 78.1684, 76.4873, 74.7635, 72.761, 71.0948, 69.6102, 67.7806, 66.2215, 64.5813, 63.0269, 61.4325, 59.8065, 58.2423, 57.2432, 55.8296, 54.2171, 52.8809, 51.3254, 49.9902, 48.6927, 47.5565, 46.136, 44.8382, 43.6345, 42.3964, 41.4211, 39.9681, 39.178, 37.9341, 36.9268, 35.5626, 34.5382, 33.6912, 32.8156, 31.6695, 30.6552, 29.7015, 28.8655, 27.9609, 27.0857, 26.105, 25.3163, 24.4872, 23.6394, 23.0484, 22.2774, 21.4877, 20.5556, 19.9736, 19.3296, 18.5628, 17.916, 17.2928, 16.6546, 16.1131, 15.4013, 14.8264, 14.3973, 13.7262, 13.2853, 12.8253, 12.2874, 11.7558, 11.2723, 10.8829, 10.4652, 9.96477, 9.6368, 9.09316, 8.84175, 8.48084, 8.05694, 7.64559, 7.29709, 7.07981, 6.70294, 6.45736, 6.10284, 5.91788, 5.5441, 5.33311, 5.06641, 4.96415, 4.6286, 4.38214, 4.2076, 4.01099, 3.81054, 3.63854, 3.43403, 3.23244, 3.08666, 2.86953, 2.74334, 2.62787, 2.48354, 2.38115, 2.26822, 2.23137, 2.1665, 2.14264, 2.10636, 2.07358, 2.05422, 2.04126, 2.00954};
   double *Npartbins = new double[nbins+1]; Npartbins[0]=0; for (int i=1; i<nbins+1; i++) Npartbins[i] = Npart[nbins-i]+0.01;
   TH1F *hNpart = new TH1F("hNpart","Npart",200,Npartbins);
   TH1F *hNpart_mid = new TH1F("hNpart_mid","Npart",200,Npartbins);
   TH1F *hNpart_fwd = new TH1F("hNpart_fwd","Npart",200,Npartbins);
   TH1F *hNcoll = new TH1F("hNcoll","Ncoll",200,Ncollbins);
   TH1F *hNcoll_mid = new TH1F("hNcoll_mid","Ncoll",200,Ncollbins);
   TH1F *hNcoll_fwd = new TH1F("hNcoll_fwd","Ncoll",200,Ncollbins);
   TH1F *hdNpart = new TH1F("hdNpart","dN/dNpart",200,Npartbins);
   TH1F *hdNpart_mid = new TH1F("hdNpart_mid","dN/dNpart",200,Npartbins);
   TH1F *hdNpart_fwd = new TH1F("hdNpart_fwd","dN/dNpart",200,Npartbins);
   TH1F *hdNcoll = new TH1F("hdNcoll","dN/dNcoll",200,Ncollbins);
   TH1F *hdNcoll_mid = new TH1F("hdNcoll_mid","dN/dNcoll",200,Ncollbins);
   TH1F *hdNcoll_fwd = new TH1F("hdNcoll_fwd","dN/dNcoll",200,Ncollbins);
   TH1F *hHF = new TH1F("hHF","E_{HF} 0-100%",100,0,10000);
   TH1F *hHF010 = new TH1F("hHF010","E_{HF} 0-10%",100,0,10000);
   TH1F *hHF1020 = new TH1F("hHF1020","E_{HF} 10-20%",100,0,10000);
   TH1F *hHF2030 = new TH1F("hHF2030","E_{HF} 20-30%",100,0,10000);
   TH1F *hHF3040 = new TH1F("hHF3040","E_{HF} 30-40%",100,0,10000);
   TH1F *hHF4050 = new TH1F("hHF4050","E_{HF} 40-50%",100,0,10000);
   TH1F *hHF50100 = new TH1F("hHF50100","E_{HF} 50-100%",100,0,10000);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (Reco_QQ_size == 0) continue;

      for (int i=0; i<Reco_QQ_size; i++) {
         if (!areMuonsInAcceptance2015(i)) continue;
         // sign
         if (Reco_QQ_sign[i]!=0) continue;
         // quality cuts
         if (!passQualityCuts2015(i)) continue;
         // trigger matching
         if (!isTriggerMatch(i,0)) continue; // HLT_HIL1DoubleMu0_v1
         TLorentzVector *tlv = (TLorentzVector*) Reco_QQ_4mom->At(i);
         double mass = tlv->M();
         if (mass>3.&&mass<3.2) {
            hVprobSig->Fill(Reco_QQ_VtxProb[i]);
         } else if (mass>2. && mass<2.5) {
            hVprobBkg->Fill(Reco_QQ_VtxProb[i]);
         }
         if (mass>2.2&&mass<4.5) {
            hCentrality->Fill(Centrality/2.);
            hNpart->Fill(Npart[Centrality]);
            hNcoll->Fill(Ncoll[Centrality]);
            if (fabs(tlv->Rapidity()<1.6)) {
               hCentrality_mid->Fill(Centrality/2.);
               hNpart_mid->Fill(Npart[Centrality]);
               hNcoll_mid->Fill(Ncoll[Centrality]);
            } else {
               hCentrality_fwd->Fill(Centrality/2.);
               hNpart_fwd->Fill(Npart[Centrality]);
               hNcoll_fwd->Fill(Ncoll[Centrality]);
            }
            hHF->Fill(SumET_HF);
            if (Centrality<20) hHF010->Fill(SumET_HF);
            else if (Centrality<40) hHF1020->Fill(SumET_HF);
            else if (Centrality<60) hHF2030->Fill(SumET_HF);
            else if (Centrality<80) hHF3040->Fill(SumET_HF);
            else if (Centrality<100) hHF4050->Fill(SumET_HF);
            else hHF50100->Fill(SumET_HF);
         }
      }
   }

   for (int i=1; i<nbins+1; i++) {
      hdNpart->SetBinContent(i,hNpart->GetBinContent(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNpart_mid->SetBinContent(i,hNpart_mid->GetBinContent(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNpart_fwd->SetBinContent(i,hNpart_fwd->GetBinContent(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNcoll->SetBinContent(i,hNcoll->GetBinContent(i)/(Ncollbins[i]-Ncollbins[i-1]));
      hdNcoll_mid->SetBinContent(i,hNcoll_mid->GetBinContent(i)/(Ncollbins[i]-Ncollbins[i-1]));
      hdNcoll_fwd->SetBinContent(i,hNcoll_fwd->GetBinContent(i)/(Ncollbins[i]-Ncollbins[i-1]));
      hdNpart->SetBinError(i,hNpart->GetBinError(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNpart_mid->SetBinError(i,hNpart_mid->GetBinError(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNpart_fwd->SetBinError(i,hNpart_fwd->GetBinError(i)/(Npartbins[i]-Npartbins[i-1]));
      hdNcoll->SetBinError(i,hNcoll->GetBinError(i)/(Ncollbins[i]-Ncollbins[i-1]));
      hdNcoll_mid->SetBinError(i,hNcoll_mid->GetBinError(i)/(Ncollbins[i]-Ncollbins[i-1]));
      hdNcoll_fwd->SetBinError(i,hNcoll_fwd->GetBinError(i)/(Ncollbins[i]-Ncollbins[i-1]));
   }

   f->Write();
   f->Close();
}
