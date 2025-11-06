// merge the contents of different files containing HypoTestInverterResult for the same bin
// Compile with g++  `root-config --cflags --libs --ldflags` -lRooStats -o hadd_hypotestinvresult hadd_hypotestinvresult.C
//
// Usage: ./hadd_hypotestinvresult output.root input1.root input2.root ...


#include "TFile.h"
#include "TString.h"
#include "RooStats/HypoTestInverterResult.h"
#include <iostream>

using namespace RooStats;
using namespace std;

int main(int argc, const char** argv) {
   if (argc<3) {
      cout << "Usage: " << argv[0] << " output.root input1.root input2.root ..." << endl;
      return 1;
   }

   TFile *output = new TFile(argv[1],"RECREATE");
   HypoTestInverterResult *res = NULL;
   TString name;

   for (int i=0; i<argc-2; i++) {
      TFile *f = new TFile(argv[i+2]);
      if (!f) {
         cout << argv[i+2] << " does not exist." << endl;
         continue;
      }
      RooStats::HypoTestInverterResult* r = static_cast<RooStats::HypoTestInverterResult*>(f->FindObjectAny("result_RFrac2Svs1S_PbPbvsPP"));
      if (!r) r = static_cast<RooStats::HypoTestInverterResult*>(f->FindObjectAny("result_RFrac2Svs1S_PbPbvsPP_P"));
      if (!r) {
         cout << "[ERROR] Could not find HypoTestInverterResult in file "<< argv[i+2] << endl;
         continue;
      }
      if (!res) {
         name = r->GetName();
         output->cd();
         res = new HypoTestInverterResult(*r,"merged");
      } else {
         res->Add(*r);
      }
      f->Close();
      delete f;
   }

   res->SetName(name);
   output->cd();
   res->Write();
   output->Write();
   output->Close();
}
