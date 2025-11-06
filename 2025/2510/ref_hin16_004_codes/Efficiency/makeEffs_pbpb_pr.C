#include "oniaEff.C"
#include "oniaEff_TnPToyStudy.C"
#include "oniaEff_pTShapeVary.C"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TDirectory.h>

#include <iostream>

using namespace std;

void makeEffs_pbpb_pr(int opt=1, bool isacc=false){
   bool ispbpb=true;
   bool isprompt=true;
   TChain *tch_jpsi_pbpb = new TChain("hionia/myTree");
   if (isacc) {
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
   } else {
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
     tch_jpsi_pbpb->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/PbPb502TeV/TTrees/OniaTree_Pythia8_JpsiMM_ptJpsi_30_Inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
   }

   // make the efficiency histos
   string dir = "files/"; // output files will be stored under this directory
   if (!isacc) {
     dir = dir+"eff/";
     if (opt==0) dir = dir+"noSF";
     else if (opt==1) dir = dir+"nominal";
     else if (opt==2) dir = dir+"trg__muid__sta";
     else if (opt==3) dir = dir+"trg__muid";
     else if (opt==4) dir = dir+"trg__sta";
     else if (opt==5) dir = dir+"trg_ptWeighting";
     else if (opt==6) dir = dir+"trg_toy";
     else if (opt==7) dir = dir+"trg__muid_toy";
     else if (opt==8) dir = dir+"trg__sta_toy";
     else if (opt==10) dir = dir+"trg_binned";
     else if (opt==11) dir = dir+"trg_plus1sigma";
     else if (opt==12) dir = dir+"trg_minus1sigma";
     else if (opt==13) dir = dir+"trg_trk_plus1sigma";
     else if (opt==14) dir = dir+"trg_trk_minus1sigma";
     else if (opt==15) dir = dir+"trg__muid_plus1sigma";
     else if (opt==16) dir = dir+"trg__muid_minus1sigma";
     else if (opt==17) dir = dir+"trg__sta_plus1sigma";
     else if (opt==18) dir = dir+"trg__sta_minus1sigma";
   } else {
     if (opt==1) dir = dir+"acc/nominal";
     else {
       cout << "cannot do acc with this option, exit" << endl;
       return;
     }
   }
   
   gSystem->mkdir(Form("./%s",dir.c_str()), kTRUE);
   cout << "Efficiencies for pbpb prompt Jpsi" << endl;
   if (opt==5) {
     oniaEff_pTShapeVary obj_jpsi_pbpb(tch_jpsi_pbpb);
     obj_jpsi_pbpb.LoopVary(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,obj_jpsi_pbpb.trg_ptWeighting);
   } else if (opt>=6 && opt<=8) {
     oniaEff_TnPToyStudy obj_jpsi_pbpb(tch_jpsi_pbpb);
     obj_jpsi_pbpb.LoopVary(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,opt);
   } else {
     oniaEff obj_jpsi_pbpb(tch_jpsi_pbpb);
     obj_jpsi_pbpb.Loop(Form("%s/histos_jpsi_pbpb.root",dir.c_str()),ispbpb,isprompt,opt,isacc);
   }

}
