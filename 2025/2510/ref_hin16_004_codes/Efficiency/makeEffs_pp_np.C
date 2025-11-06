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

void makeEffs_pp_np(int opt=1, bool isacc=false){
   bool ispbpb=false;
   bool isprompt=false;
   TChain *tch_npjpsi_pp = new TChain("hionia/myTree");
   if (isacc)
     tch_npjpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");
   else
     tch_npjpsi_pp->Add("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_Extended.root");

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
   cout << "Efficiencies for pp non-prompt Jpsi" << endl;
   if (opt==5) {
     oniaEff_pTShapeVary obj_npjpsi_pp(tch_npjpsi_pp);
     obj_npjpsi_pp.LoopVary(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,obj_npjpsi_pp.trg_ptWeighting);
   } else if (opt>=6 && opt<=8) {
     oniaEff_TnPToyStudy obj_npjpsi_pp(tch_npjpsi_pp);
     obj_npjpsi_pp.LoopVary(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,opt);
   } else {
     oniaEff obj_npjpsi_pp(tch_npjpsi_pp);
     obj_npjpsi_pp.Loop(Form("%s/histos_npjpsi_pp.root",dir.c_str()),ispbpb,isprompt,opt,isacc);
   }

}
