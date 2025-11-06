#include <iostream>
#include <TStopwatch.h>
#include <string>
#include <TFile.h>
#include <RooDataSet.h>
#include <TString.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>
// #include <XX.h>

using std::cout;
using std::vector;
using std::string;

void make_sub_datasets()
{
  cout << "=== start make_sub_datasets() ===\n";
  TStopwatch time;
  time.Start();

  // kinematcis
  float massLow = 2.6, massHigh = 3.5;
  float yLow = 0, yHigh = 2.4;
  int cLow = 0, cHigh = 180;
  double ptLow = 6.5, ptHigh = 7.5;

  // --- read input ---
  cout << "\n=== Import inputs ===\n";
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile inputData(fileNameData.c_str());
  cout << fileNameData.c_str() << "\n";
  RooDataSet *data = (RooDataSet *)inputData.Get("dataset");
  data->SetName("data");

  // --- pass common cuts ---
  char commonCut[3000];
  sprintf(commonCut, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f)"

          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

  auto redData = (RooDataSet *)data->reduce(commonCut);

//  && ctau3D >= %.3f && ctau3D < %.3f)
// , ctMin, ctMax
  // --- main loop ---
  vector<double> ctEdge = {-1.0, -0.5, -0.2, -0.1, -0.07, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.5, 1.0, 2.0};
  const char *v = "ctau3D";

  for (size_t i = 0; i + 1 < ctEdge.size();  ++i) {
    double lo = ctEdge[i], hi = ctEdge[i + 1];
    TString cut = Form("ctau3D > %.2f && ctau3D < %.2f", lo, hi);
    auto ds = (RooDataSet *)redData->reduce(cut);

    TFile f(Form("roots/pt%.1f_%.1f_ctau%.2f_%.2f.root", ptLow, ptHigh, lo, hi), "recreate");
    ds->Write("dataset");
    f.Close();
    delete ds;
  }
  // make 3 dataset according to actau cut
  // save and cloe the output

  cout << "\n=== finish make_sub_datasets() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}