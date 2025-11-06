#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"


using namespace RooFit;

void plotVar(RooDataSet *dataset, const char *varName, const char *title, double rangeMin, double rangeMax, const char *outname);


void r0_draw_basic_plots()
{
  gSystem->mkdir("figs/basic", true);

  // apply weight? -> later

  // Read imputs
  // -----------
  // --- NP MC ---
  string fileNamePrMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";
  TFile fInPRMC(fileNamePrMc.c_str());
  cout << fileNamePrMc.c_str() << endl;
  RooDataSet *dataPRMC = (RooDataSet *)fInPRMC.Get("dataset");
  dataPRMC->SetName("dataPRMC");

  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");


  // Draw raw plots
  // ----------
  // draw raw plots here
  RooRealVar *ctau3D = (RooRealVar *)data->get()->find("ctau3D");
  RooRealVar *ctau3DErr = (RooRealVar *)data->get()->find("ctau3DErr");

  plotVar(data, "ctau3D", "", -1, 2, "ctau3D_raw_data.png");
  plotVar(data, "ctau3DErr", "", -0.1, 0.3, "ctauErr_raw_data.png");
  plotVar(data, "ctau3DRes", "", -10, 50, "ctauRes_raw_data_n10_50.png");
  plotVar(data, "ctau3DRes", "", -10, 10, "ctauRes_raw_data_n10_10.png");
  plotVar(dataPRMC, "ctau3D", "", -1, 1, "ctau3D_raw_MC.png");
  plotVar(dataPRMC, "ctau3DErr", "", -0.1, 0.3, "ctauErr_MC_raw.png");
  plotVar(dataPRMC, "ctau3DRes", "", -10, 10, "ctauRes_MC_raw.png");


  // Draw plot with kinematic cut
  // ----------
  double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
  double ctMin = -1, ctMax = 4, errmin = 0.008, errmax = 0.3;

  char reduceDS[3000], reduceDS_woCtErr[3000];
  // build cuts
  sprintf(reduceDS, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"

          "&& (recoQQsign == 0)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errmin, errmax);

  // sprintf(reduceNpMc,
  //         "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

  //         "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3Dtrue >= %.3f && ctau3Dtrue < %.3f)",
  //         ptLow, ptHigh, yLow, yHigh, massLow, massHigh, 0.0001, ctMax);

  RooDataSet *redPRMC, *redData;
  redPRMC = (RooDataSet *)dataPRMC->reduce(reduceDS); // smae cut with Data
  redData = (RooDataSet *)data->reduce(reduceDS);

  plotVar(redData, "ctau3D", "", -1, 2, "ctau3D_KineCut_data.png");
  plotVar(redData, "ctau3DErr", "", -0.1, 0.3, "ctauErr_KineCut_data.png");
  plotVar(redPRMC, "ctau3D", "", -1, 1, "ctau3D_KineCut_MC.png");
  plotVar(redPRMC, "ctau3DErr", "", -0.1, 0.3, "ctauErr_MC_KineCut.png");
}

void plotVar(RooDataSet *dataset,
             const char *varName,
             const char *title,
             double rangeMin, double rangeMax,
             const char *outname)
{
  // get observable
  RooRealVar *var = (RooRealVar *)dataset->get()->find(varName);
  if (!var)
  {
    std::cerr << "Error: variable " << varName << " not found in dataset!" << std::endl;
    return;
  }

  // set range
  var->setRange(rangeMin, rangeMax);

  // draw
  RooPlot *frame = var->frame(Title(title));
  dataset->plotOn(frame, Name("data"));

  TCanvas *c = new TCanvas(Form("c_%s", var->GetName()), "", 800, 600);
  frame->Draw();
  string outpath = "figs/basic/" + string(outname);
  c->SaveAs(outpath.c_str());

  delete c;
}