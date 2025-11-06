#include "Macros/Utilities/bfrac.h"
#include "results2tree.C"

#define normalCutsDir "nominal"
#define invCutsDir "nonprompt"

void bfrac_syst(int mode, const char* systfilename="Systematics/csv/syst_PbPb_bhad.csv", double nsigma=0, double chiNP=1.) {
   // mode: 0 for the multiplicative factor (1/(1+B)), 1 for the additive (B/(1+B))
   // nsigma: vary the value by nsigma (e.g. +1)
   // chiNP: assume this value for the non-prompt double ratio (to degrade systs)

   string fullSystFileName(systfilename);
   ofstream file(fullSystFileName.c_str());
   file << "B-had" << endl;
   TFile *f = new TFile(treeFileName(normalCutsDir));
   if (!f || !f->IsOpen()) {
      results2tree(normalCutsDir);
      f = new TFile(treeFileName(normalCutsDir));
      if (!f) return;
   }
   TTree *tr = (TTree*) f->Get("fitresults");
   if (!tr) return;
   float ptmin, ptmax, ymin, ymax, centmin, centmax;
   char collSystem[5];
   tr->SetBranchAddress("ptmin",&ptmin);
   tr->SetBranchAddress("ptmax",&ptmax);
   tr->SetBranchAddress("ymin",&ymin);
   tr->SetBranchAddress("ymax",&ymax);
   tr->SetBranchAddress("centmin",&centmin);
   tr->SetBranchAddress("centmax",&centmax);
   tr->SetBranchAddress("collSystem",collSystem);

   int ntr = tr->GetEntries();
   for (int i=0; i<ntr; i++) {
      tr->GetEntry(i);
      if (string(collSystem) != "PbPb") continue;

      anabin trbin(ymin, ymax, ptmin, ptmax, centmin, centmax);

      RooRealVar B = Bval(normalCutsDir,invCutsDir,trbin);
      RooRealVar alpha1 = alphaval(normalCutsDir,invCutsDir,trbin,false);
      RooRealVar alpha2 = alphaval(normalCutsDir,invCutsDir,trbin,true);

      double Bval_syst = B.getVal() + nsigma*B.getError();
      double alpha1val_syst = alpha1.getVal() + nsigma*alpha1.getError();
      double alpha2val_syst = alpha2.getVal() + nsigma*alpha2.getError();

      double value=(1.+alpha1val_syst)/(1.+alpha2val_syst);
      if (mode==0) value = fabs(1.-value/(1.+Bval_syst));
      else value = value * (Bval_syst/(1.+Bval_syst)) * (alpha2val_syst/alpha1val_syst) * chiNP;


      file << ymin << ", " << ymax << ", " 
         << ptmin << ", " << ptmax << ", " 
         << centmin << ", " << centmax << ", " 
         << value << endl;

      // cout << ymin << ", " << ymax << ", " 
      //    << ptmin << ", " << ptmax << ", " 
      //    << centmin << ", " << centmax << ", " 
      //    << " -> B=" << B.getVal() << ", alpha1=" << alpha1.getVal() << ", alpha2=" << alpha2.getVal() << ". used "
      //    << " -> B=" << Bval_syst << ", alpha1=" << alpha1val_syst << ", alpha2=" << alpha2val_syst << endl;
   }

   f->Close();
   delete f;

   file.close();
   cout << "Closed " << fullSystFileName << endl << endl;

}
