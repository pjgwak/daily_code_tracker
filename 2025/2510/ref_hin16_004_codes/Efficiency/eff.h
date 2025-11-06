#ifndef eff_h
#define eff_h

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "../Fitter/Macros/Utilities/bin.h"
#include <iostream>

using namespace std;

double ljpsieff(const char* part, const char* collSystem, anabin thebin, const char* prependPath = "") {
   TString path = Form("files/");
   if (TString(prependPath) != "") path = TString(prependPath) + "/" + path;
   TString filename = path + "histos_" + TString(part) + "_" + TString(collSystem) + ".root";
   TFile *f = TFile::Open(filename);
   if (!f || !f->IsOpen()) {
      cout << "Error, could not find " << filename << endl;
      return 0;
   }

   bool iscent = ((thebin.centbin()==binI(0,200)) && TString(part) != "pp");

   TString hnumname("hnumptdepcut_");
   if (iscent) hnumname += "cent";
   else hnumname += "pt";
   if (thebin.rapbin()==binF(0,1.6)) hnumname += "mid";
   else hnumname += "fwd";
   TH1F *hnum = (TH1F*) f->Get(hnumname.Data());

   TString hdenname("hnum_");
   if (iscent) hdenname += "cent";
   else hdenname += "pt";
   if (thebin.rapbin()==binF(0,1.6)) hdenname += "mid";
   else hdenname += "fwd";
   TH1F *hden = (TH1F*) f->Get(hdenname.Data());

   if (!hnum || !hden) {
      cout << "Error, hists " << hnumname << ", " << hdenname << " not found in " << filename << endl;
      return 0;
   }

   double xval = 0.1 + (iscent ? thebin.centbin().low() : thebin.ptbin().low());
   double numval = hnum->GetBinContent(hnum->FindBin(xval));
   double denval = hden->GetBinContent(hden->FindBin(xval));
   double eff = numval / denval;
   f->Close(); delete f;
   return eff;
}

#endif // ifndef eff_h
