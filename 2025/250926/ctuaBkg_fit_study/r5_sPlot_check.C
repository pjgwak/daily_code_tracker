#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;

// Draw plot with kinematic cut
// ----------
double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
double ctMin = -10, ctMax = 10;
double errMin = 0, errMax = 10;

// 0.016 - 0.03: 33초

// 1. 산점도
void drawScatter(RooDataSet* data, RooRealVar& mass, RooRealVar& ctauErr) {
  TCanvas* c1 = new TCanvas("c1","Scatter mass vs ctauErr",800,600);

  TH1* h2 = data->createHistogram("h2_mass_ctauErr", mass,
                                   Binning(50),
                                   YVar(ctauErr, Binning(50)));

  h2->SetTitle("mass vs ctauErr;mass;ctauErr");
  h2->Draw("colz");

  c1->SaveAs("scatter_mass_ctauErr.png");
}

// 2. 상관계수
double computeCorrelation(RooDataSet* data, RooRealVar& mass, RooRealVar& ctauErr) {
  double corr = data->correlation(mass, ctauErr);
  std::cout << "[INFO] Correlation(mass, ctauErr) = " << corr << std::endl;
  return corr;
}

// 3. 프로파일
void drawProfile(RooDataSet* data, RooRealVar& mass, RooRealVar& ctauErr) {
  TTree* t = (TTree*) data->GetClonedTree();
  TCanvas* c2 = new TCanvas("c2","Profile ctauErr vs mass",800,600);
  t->Draw(Form("%s:%s >> prof(50)", ctauErr.GetName(), mass.GetName()), "", "prof");
  TProfile* prof = (TProfile*) gDirectory->Get("prof");
  prof->SetTitle("Profile of ctauErr vs mass;mass;mean(ctauErr)");
  prof->Draw("E1");
  c2->SaveAs("profile_mass_ctauErr.png");
}

void r5_sPlot_check()
{
  gSystem->mkdir("figs/fit", true);

  RooRealVar mass("mass", "mass", 2.6, 3.5, "GeV/c^{2}");
  RooRealVar ctau3DErr("ctau3DErr", "", errMin, errMax, "mm");
  
  // RooBinning binningMass(2.6, 3.5);
  // binningMass.addUniform(nMassBins, 2.6, 3.5);
  // mass.setBinning(binningMass);

  // keep ctau bin width about 83 micro-meter
  

   TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f) && (ctau3D>%.3f && ctau3D < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);
  TString osCut = "(recoQQsign == 0)";
  TString fullCut = Form("%s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data());

  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  auto data_tmp = (RooDataSet *)ds->reduce(Cut(fullCut));
  RooArgSet obs(mass, ctau3DErr);
  auto data1 = new RooDataSet("data1", "dataset with local vars", obs, Import(*data_tmp));

  drawScatter(data1, mass, ctau3DErr);
  computeCorrelation(data1, mass, ctau3DErr);
  drawProfile(data1, mass, ctau3DErr);
  // char reduceDS[3000], reduceDS_woCtErr[3000];
  // // build cuts
  // sprintf(reduceDS, // no multiplicty case
  //         "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

  //         "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"

  //         "&& (recoQQsign == 0)",
  //         ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errMin, errMax);

  // RooDataSet *redData;
  // redData = (RooDataSet *)data->reduce(reduceDS);

  // // sideband
  // RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  // RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  // RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  // RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");

  // sPlot_check(redDataSBL, ctMin, ctMax, "conditional_SBL.png");
}