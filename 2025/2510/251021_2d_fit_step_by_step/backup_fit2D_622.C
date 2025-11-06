#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
// #include "RooHistPdfConv.h"
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>
#include <TStopwatch.h>

using namespace std;
using namespace RooFit;

/////////////////////////////////// fit quality information ///////////////////////////////////
// 1) for mass distributions
double UnNormChi2_mass = 0;  // un normalized chi2
int nFitParam_mass = 0;      // # of floating fit parameters
int nFullBinsPull_mass = 0;  // # of bins in pull dist
int Dof_mass = 0;            // ndf = nFullBinsPull - nFitPram
double Chi2_mass = 0;        // chi2/ndf
double theChi2Prob_mass = 0; // check the chi2 fit probability
// 2) for time distributions
double UnNormChi2_time = 0;  // un normalized chi2
int nFitParam_time = 0;      // # of floating fit parameters
int nFullBinsPull_time = 0;  // # of bins in pull dist
int Dof_time = 0;            // ndf = nFullBinsPull - nFitPram
double Chi2_time = 0;        // chi2/ndf
double theChi2Prob_time = 0; // check the chi2 fit probability
// 3) for time side distributions
double UnNormChi2_side = 0;  // un normalized chi2
int nFitParam_side = 0;      // # of floating fit parameters
int nFullBinsPull_side = 0;  // # of bins in pull dist
int Dof_side = 0;            // ndf = nFullBinsPull - nFitPram
double Chi2_side = 0;        // chi2/ndf
double theChi2Prob_side = 0; // check the chi2 fit probability



////// **** Global objects for legend
// TGraphErrors *gfake1;
// TH1F hfake11, hfake21, hfake31, hfake311, hfake41;
// TLatex *t = new TLatex();
// t->SetNDC();
// t->SetTextAlign(12);
// t->SetTextFont(42);
// Double_t fx[2], fy[2], fex[2], fey[2];
// //// black points (data, etc)
// gfake1 = new TGraphErrors(2, fx, fy, fex, fey);
// gfake1->SetMarkerStyle(20);
// gfake1->SetMarkerSize(1);
// //// blue line (background)
// hfake11 = TH1F("hfake11", "hfake1", 100, 200, 300);
// // hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
// hfake11.SetLineColor(kBlue - 2);
// hfake11.SetLineWidth(4);
// hfake11.SetLineStyle(7);
// hfake11.SetFillColor(kBlue - 10);
// hfake11.SetFillStyle(1001);
// //// black line (total fit)
// hfake21 = TH1F("hfake21", "hfake2", 100, 200, 300);
// // hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
// hfake21.SetLineColor(kBlack);
// hfake21.SetLineWidth(4);
// hfake21.SetFillColor(kGray + 2);
// hfake21.SetFillStyle(3354);
// //// red line (nonprompt)
// hfake31 = TH1F("hfake31", "hfake3", 100, 200, 300);
// // hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);
// hfake31.SetLineColor(kPink - 6);
// hfake31.SetMarkerStyle(kCircle);
// hfake31.SetLineWidth(4);
// hfake31.SetMarkerColor(kPink - 6);
// hfake31.SetLineStyle(12);
// hfake31.SetFillColor(kRed - 7);
// hfake31.SetFillStyle(3345);
// hfake311 = TH1F("hfake311", "hfake311", 100, 200, 300);
// hfake311.SetLineColor(kPink - 6);
// hfake311.SetMarkerStyle(kCircle);
// hfake311.SetLineWidth(4);
// hfake311.SetMarkerColor(kPink - 6);
// hfake311.SetLineStyle(12);
// hfake311.SetFillColor(kRed - 7);
// hfake311.SetFillStyle(3444);
// //// green line (prompt)
// hfake41 = TH1F("hfake41", "hfake4", 100, 200, 300);
// hfake41.SetLineColor(kGreen + 3);
// hfake41.SetMarkerStyle(kCircle);
// hfake41.SetLineWidth(4);
// hfake41.SetMarkerColor(kGreen + 3);
// hfake41.SetLineStyle(kDashDotted);
// hfake41.SetFillColor(kGreen - 7);
// hfake41.SetFillStyle(3444);

// Helper functions
// ----------------
void drawCtauSBPlots(RooWorkspace *ws, RooDataSet *redDataSB, RooDataHist *binDataCtErrSB, RooFitResult *fitCt_Bkg, float lmin, float lmax, double *UnNormChi2_side_t, int *nFitParam_side_t, int *nFullBinsPull_side_t, int *Dof_side_t, double *Chi2_side_t);

void setWSRange(RooWorkspace *ws, float lmin, float lmax, float errmin, float errmax)
{
  float minRangeForPF = -4 * errmax;
  if (minRangeForPF < -lmin)
    minRangeForPF = -lmin;

  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("ctau3Dtrue")->setRange(-lmin, lmax);
  ws->var("ctau3Dtrue")->setRange("trueRange", -lmin, lmax);
  ws->var("ctau3D")->setRange("promptMCfit", -0.2, 0.2);
  
  // ws->var("ctau3D")->setRange("promptMCfit", minRangeForPF, 4 * errmax);
  ws->var("ctau3D")->setRange(-lmin, lmax);
  // ws->var("ctau3D")->setRange("CBkgFit", lmin, lmax);
  ws->var("ctau3DErr")->setRange(errmin, errmax);

  return;
}

RooBinning setCtBinning(float lmin, float lmax)
{
  RooBinning rbct(-lmin, lmax);
  // if (lmax + lmin > 4.9)
  // {
  //   rbct.addBoundary(-1.5);
  //   rbct.addBoundary(-1.0);
  //   rbct.addBoundary(-0.8);
  //   rbct.addBoundary(-0.6);
  //   rbct.addBoundary(-0.5);
  //   rbct.addUniform(6, -0.5, -0.2);
  //   rbct.addUniform(12, -0.2, 0.1);
  //   rbct.addUniform(8, 0.1, 0.5);
  //   rbct.addUniform(5, 0.5, 1.0);
  //   rbct.addUniform(15, 1.0, lmax);
  // }
  // else if (lmax + lmin > 4.4)
  // { // this is what we use KYO!!!
  //   rbct.addBoundary(-1.5);
  //   rbct.addBoundary(-1.0);
  //   rbct.addBoundary(-0.8);
  //   rbct.addBoundary(-0.6);
  //   rbct.addBoundary(-0.5);
  //   rbct.addUniform(9, -0.5, -0.2);
  //   rbct.addUniform(20, -0.2, 0.1);
  //   rbct.addUniform(11, 0.1, 0.5);
  //   rbct.addUniform(5, 0.5, 1.0);
  //   rbct.addUniform(5, 1.0, lmax);
  // }
  // else
  // {
  //   rbct.addBoundary(-lmin);
  //   rbct.addBoundary(-0.7);
  //   rbct.addBoundary(-0.6);
  //   rbct.addBoundary(-0.5);
  //   rbct.addUniform(9, -0.5, -0.2);
  //   rbct.addUniform(28, -0.2, 0.1);
  //   rbct.addUniform(11, 0.1, 0.5);
  //   rbct.addUniform(15, 0.5, 1.2);
  //   rbct.addUniform(8, 1.2, lmax);
  // }

  
  rbct.addUniform(20, -0.2, 0.1);
  rbct.addUniform(11, 0.1, 0.4);
  // rbct.addUniform(5, 0.5, 1.0);
  // rbct.addUniform(5, 1.0, lmax);
  return rbct;
}

void defineMassSig(RooWorkspace *ws)
{
  // gauss: meanSig, sigmaSgi1
  ws->factory("Gaussian::G1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.008,0.1])");
  
  // cb: meanSig, sigmaSig2
  ws->factory("CBShape::CB1Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha[1.9,1.2,2.8],enne[2.5,1.0,4.0])");
  
  // Sum of G1 and CB1
  ws->factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*G1Sig,CB1Sig)");
  return;
}

void defineMassBkg(RooWorkspace *ws)
{
  // 1st order polynomial
  ws->factory("Polynomial::polBkg(mass,{coefPol[-0.05,-5.,5.]})");

  // expo
  ws->factory("Exponential::expBkg(mass,coefExp[-1.,-3.,1.])");
  return;
}
/////////////////////////////////////////////////////////

void defineCtPRRes(RooWorkspace *ws)
{
  ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[0.2, 0.1, 8],one[1.0],ctau3DErr)");
  ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.05,0.001, 1],one,ctau3DErr)");
  ws->factory("AddModel::CtPRRes({GW_PRRes,GN_PRRes},{fracRes[0.9,0.01,0.99999]})");

  // ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 10],one[1.0],ctau3DErr)");
  // ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("GaussModel::GVW_PRRes(ctau3D,meanPRResVW[0],sigmaPRResVW[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("AddModel::CtPRRes({GN_PRRes, GW_PRRes, GVW_PRRes},{fracResN[0.01,0.001,0.999], fracResW[0.01,0.001,0.999]})");
  return;
}
/////////////////////////////////////////////////////////

void defineCtBkg(RooWorkspace *ws)
{
  // ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.42,0.05,1.5],CtPRRes,RooDecay::SingleSided)");
  // ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.79,0.02,1.5],CtPRRes,RooDecay::Flipped)");
  // ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.69,0.02,5.0],CtPRRes,RooDecay::DoubleSided)");

  // -- original ---
  // ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.42,0.01, 2],CtPRRes,RooDecay::SingleSided)");
  // ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.02,0.00001, 2],CtPRRes,RooDecay::Flipped)");
  // ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.2,0.01, 2],CtPRRes,RooDecay::DoubleSided)");
  // ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.3,0.01,0.99]*CtBkgPos,CtBkgNeg)");
  // ws->factory("SUM::CtBkgSum2(fracCtBkg2[0.1,0.01,0.99]*CtBkgSum1,CtBkgDbl)");
  // ws->factory("SUM::CtBkgTot(fracCtBkg3[0.1,0.01,0.99]*CtPRRes,CtBkgSum2)");

  // -- 2L, 2M, 2R ---
  ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.42,0.01, 10],CtPRRes,RooDecay::SingleSided)");
  ws->factory("Decay::CtBkgPos2(ctau3D,lambdap2[0.42,0.01, 2],CtPRRes,RooDecay::SingleSided)");

  ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.02,0.001, 2],CtPRRes,RooDecay::Flipped)");
  ws->factory("Decay::CtBkgNeg2(ctau3D,lambdam2[0.02,0.001, 2],CtPRRes,RooDecay::DoubleSided)");

  ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.2,0.01, 2],CtPRRes,RooDecay::DoubleSided)");
  ws->factory("Decay::CtBkgDbl2(ctau3D,lambdasym2[0.2,0.00001, 2],CtPRRes,RooDecay::DoubleSided)");

  ws->factory("SUM::CtBkgTot(fracCtBkg1[0.2,0.01,0.99]*CtBkgPos_PEE,\
                                 fracCtBkg2[0.2,0.01,0.99]*CtBkgPos2_PEE,\
                                 fracCtBkg3[0.2,0.01,0.99]*CtBkgNeg_PEE,\
                                 fracCtBkg4[0.2,0.01,0.99]*CtBkgNeg2_PEE,\
                                 fracCtBkg5[0.2,0.01,0.99]*CtBkgDbl_PEE,\
                                 fracCtBkg6[0.2,0.01,0.99]*CtBkgDbl2_PEE,\
                                 CtPRRes_PEE)");

  // ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.3,0.01,0.99]*CtBkgPos,CtBkgNeg)");
  // ws->factory("SUM::CtBkgSum2(fracCtBkg2[0.1,0.01,0.99]*CtBkgSum1,CtBkgDbl)");
  // ws->factory("SUM::CtBkgSum3(fracCtBkg3[0.1,0.01,0.99]*CtBkgSum2,CtBkgNeg2)");
  // ws->factory("SUM::CtBkgTot(fracCtBkg4[0.1,0.01,0.99]*CtPRRes,CtBkgSum3)");
  return;
}

RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr")
{

  if (binSIG->numEntries() != binSB->numEntries())
  {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }
  RooDataHist *binScaleBKG = new RooDataHist("binScaleBKG", "scaled SB", RooArgSet(*(ws->var(varName.c_str()))));

  //// **** bin-by-bin scaling
  const RooArgSet *argSIG;
  const RooArgSet *argSB;
  for (Int_t i = 0; i < binSIG->numEntries(); i++)
  {
    argSIG = binSIG->get(i);
    argSB = binSB->get(i);
    RooRealVar *thisVar = (RooRealVar *)argSIG->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    //// *** set minimum as 0.1 to prevent empty PDF
    float wBkg = binSB->weight(*argSB, 0, false);
    if (wBkg <= 0.1)
      wBkg = 0.1;
    binScaleBKG->add(RooArgSet(*(ws->var(varName.c_str()))), wBkg);
    float newWeight = binSIG->weight(*argSIG, 0, false) - scalefactor * binSB->weight(*argSB, 0, false);
    if (newWeight <= 0.1)
      newWeight = 0.1;
    binSubtrSIG->add(RooArgSet(*(ws->var(varName.c_str()))), newWeight);
  }
  return binScaleBKG;
}

void getCtErrRange(RooDataSet *data, const char *t_reduceDS_woCtErr, float lmin, float lmax, float *errmin, float *errmax)
{
  RooWorkspace *ws = new RooWorkspace("ctauerrorcheckWS");
  RooDataSet *redDataCut = (RooDataSet *)data->reduce(t_reduceDS_woCtErr);
  ws->import(*redDataCut);

  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setBins(45);
  // ws->var("ctau3D")->setRange(-lmin, lmax);
  ws->var("ctau3D")->setRange(-2, 4);
  ws->var("ctau3DErr")->setRange(0.0, 0.992);
  ws->var("ctau3DErr")->setBins(124);

  cout << " -*-*-*-*-*-*-*-*-*-*-*-*-*- getCtErrRange -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  " << endl;
  cout << " *** DATA :: N events to fit (woCtErrRange) : " << redDataCut->sumEntries() << endl;

  RooDataSet *redDataSB = (RooDataSet *)redDataCut->reduce("mass < 2.9 || mass > 3.3");
  RooDataSet *redDataSIG = (RooDataSet *)redDataCut->reduce("mass > 2.9 && mass < 3.3");

  //// *** RooDataHist
  RooDataHist *tbinDataCtErrSB = new RooDataHist("tbinDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *tbinDataCtErrSIG = new RooDataHist("tbinDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  //// *** mass fit to get coefExp of coefPol for scaleF
  defineMassBkg(ws);
  defineMassSig(ws);
  struct PARAM
  {
    double fracG1;
    double fracG1Err;
    double meanSig;
    double meanSigErr;
    double sigmaSig1;
    double sigmaSig1Err;
    double sigmaSig2;
    double sigmaSig2Err;
    double alpha;
    double alphaErr;
    double enne;
    double enneErr;
  };
  double cutValue_merged;

  char funct[200];
  double initBkg = redDataSB->sumEntries() * 9.0 / 5.0;
  double initSig = redDataCut->sumEntries() - initBkg;
  
  sprintf(funct, "SUM::MassPDF(NSig[%f,1.0,50000.0]*G1CB1Sig,NBkg[%f,1.0,500000.0]*expBkg)", initSig, initBkg);
  ws->factory(funct);

  RooRealVar *mass = ws->var("mass");
  mass->setBins(120);
  RooDataHist dh("dh", "binned data", RooArgSet(*mass), *redDataCut);

  RooFitResult *fitres = ws->pdf("MassPDF")->fitTo(
      dh,
      Extended(kTRUE),
      Minos(kFALSE),
      Save(kTRUE),
      SumW2Error(kFALSE),
      // IntegrateBins(1e-3),
      Offset(kTRUE),
      NumCPU(32));

  // ws->pdf("MassPDF")->fitTo(*redDataCut, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(32));

  //// ****  scaleF to scale down ct err dist in 2.9-3.3 GeV/c2
  float bc;
  // if (!inOpt.mBkgFunct.compare("expBkg"))
  bc = ws->var("coefExp")->getVal();
  // else if (!inOpt.mBkgFunct.compare("polBkg"))
  //   bc = ws->var("coefPol")->getVal();
  float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  //// *** (tbinSubtractedSIG) = (tbinDataCtErrSIG) - scaleF*(tbinDataCtErrSB)
  RooDataHist *tbinSubtractedSIG = new RooDataHist("tbinSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  RooDataHist *tbinScaledBKG = subtractSidebands(ws, tbinSubtractedSIG, tbinDataCtErrSIG, tbinDataCtErrSB, scaleF, "ctau3DErr");

  //// **** Check the minimum and maximum of the ctau error in signal and background regions
  TH1 *histDataCtErrSIG = tbinDataCtErrSIG->createHistogram("histDataCtErrSIG", *ws->var("ctau3DErr"));
  TH1 *histSubtractedSIG = tbinSubtractedSIG->createHistogram("histSubtractedSIG", *ws->var("ctau3DErr"));
  TH1 *histScaledBKG = tbinScaledBKG->createHistogram("histScaledBKG", *ws->var("ctau3DErr"));

  double minSig = 0.5, maxSig = 0.0, minBkg = 0.5, maxBkg = 0.0;
  double cutValue = 0.2;

  int maxBinSig = histSubtractedSIG->GetMaximumBin();
  int maxBinBkg = histScaledBKG->GetMaximumBin();

  minSig = histSubtractedSIG->GetBinLowEdge(maxBinSig);
  minBkg = histScaledBKG->GetBinLowEdge(maxBinBkg);
  // pick up lower bound next to other non-zero bins
  for (int xbins = maxBinSig; xbins > 0; xbins--)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      minSig = histSubtractedSIG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG low edge: " << histSubtractedSIG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins > 0; xbins--)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      minBkg = histScaledBKG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG low edge: " << histScaledBKG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }

  // pick up upper bound next to other non-zero bins
  maxSig = histSubtractedSIG->GetBinLowEdge(maxBinSig + 1);
  maxBkg = histScaledBKG->GetBinLowEdge(maxBinBkg + 1);
  for (int xbins = maxBinSig; xbins < histSubtractedSIG->GetNbinsX(); xbins++)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      maxSig = histSubtractedSIG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG upper edge: " << histSubtractedSIG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinSig; xbins < histScaledBKG->GetNbinsX(); xbins++)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      maxBkg = histScaledBKG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG upper edge: " << histScaledBKG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }

  // choose the higher lower limit, lower upper limit
  double tmpMin = 0, tmpMax = 0;
  if (minSig > minBkg)
    tmpMin = minSig;
  else
    tmpMin = minBkg;
  if (maxSig < maxBkg)
    tmpMax = maxSig;
  else
    tmpMax = maxBkg;

  // round off lower limit -> allow more entries on the lower limits
  tmpMin = TMath::Floor(tmpMin * 1000);
  tmpMin = tmpMin / (double)1000.0;

  // round up upper limit -> allow more entries on the upper limits
  tmpMax = TMath::Ceil(tmpMax * 1000);
  tmpMax = tmpMax / (double)1000.0;

  char reduceDS[1000];
  sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
  if (redDataTmp->sumEntries() < redDataCut->sumEntries() * 0.9)
  { // if ctau error range cuts off >10% events
    delete redDataTmp;
    sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
    tmpMin = minSig;
    tmpMax = maxSig;
  }
  if ((tmpMax - tmpMin) < 0.008)
  {
    cout << "getCtErrRange:: Maximum is less than minimum! Possibly there are few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  //// *** draw final ctau error plot
  TCanvas c0("ctau_err", "ctau_err", 500, 500);
  c0.Draw();
  c0.cd();
  c0.SetLogy(1);
  RooPlot *errframe2 = ws->var("ctau3DErr")->frame();
  // tbinDataCtErrSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kRed),LineColor(kRed));
  // tbinDataCtErrSB->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kGreen+2),LineColor(kGreen+2),MarkerStyle(24));
  // tbinScaledBKG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kBlue),MarkerStyle(24),LineColor(kBlue));
  // tbinSubtractedSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),LineColor(kWhite));
  const double max = errframe2->GetMaximum() * 1.3;
  errframe2->SetMaximum(max);
  errframe2->SetMinimum(0.2);
  errframe2->Draw();
  histDataCtErrSIG->SetMarkerColor(kRed);
  histDataCtErrSIG->SetLineColor(kWhite);
  histDataCtErrSIG->SetMarkerStyle(24);
  histDataCtErrSIG->GetXaxis()->CenterTitle(1);
  histDataCtErrSIG->GetYaxis()->CenterTitle(1);
  histDataCtErrSIG->Draw("pe");
  histScaledBKG->SetMarkerColor(kBlue);
  histScaledBKG->SetLineColor(kWhite);
  histScaledBKG->SetMarkerStyle(24);
  histScaledBKG->Draw("pe same");
  histSubtractedSIG->SetLineColor(kWhite);
  histSubtractedSIG->Draw("pe same");

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetTextSize(0.035);

  // t->DrawLatex(0.92, 0.84, opt.rapString);
  // t->DrawLatex(0.92, 0.78, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.92, 0.72, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.92, 0.72, opt.etString);

  char comment[400];
  sprintf(comment, "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
  t->SetTextSize(0.04);
  t->SetTextColor(kRed);
  t->DrawLatex(0.92, 0.6, comment);
  t->SetTextColor(kBlack);

  TLegend legsb(0.6, 0.19, 0.9, 0.35, NULL, "brNDC");
  legsb.SetFillStyle(0);
  legsb.SetBorderSize(0);
  legsb.SetShadowColor(0);
  legsb.SetMargin(0.2);
  legsb.SetTextFont(42);
  legsb.SetTextSize(0.035);
  legsb.AddEntry(histDataCtErrSIG, "sig cands", "p");
  legsb.AddEntry(histScaledBKG, "scaled bkg", "p");
  legsb.AddEntry(histSubtractedSIG, "sig (= cands - bkg)", "p");
  legsb.Draw("same");

  string titlestr =  "_CtErrGetRange_Log.pdf";
  c0.SaveAs(titlestr.c_str());

  *errmin = tmpMin;
  *errmax = tmpMax;
  //  cout << "getCtErrRange:: " << t_reduceDS_woCtErr << " " << lmin << " " << lmax << " " << *errmin << " " << *errmax << endl;

  delete ws;
  delete redDataCut;
  delete redDataTmp;
  // delete binData;
  // delete binDataCtErr;
  // delete binDataSB;
  delete tbinDataCtErrSB;
  delete tbinDataCtErrSIG;
  delete tbinSubtractedSIG;
  delete tbinScaledBKG;
  delete t;
}

void drawCtNPTrue(RooWorkspace *ws, RooDataSet *redNPMC, string titlestr)
{

  RooDataSet *redNPMCTrue = (RooDataSet *)redNPMC->reduce(RooArgSet(*(ws->var("ctau3Dtrue"))));
  //// **** fitting
  cout << " " << endl;
  cout << " -*-*-*-*-*-*-*-*-*-*-*-*-*- !!!!!FITTING!!!!! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  " << endl;
  ws->pdf("CtNPTrue")->fitTo(*redNPMCTrue, Minos(0), NumCPU(32), Range("trueRange")); //, 

  //// **** Draw
  // RooPlot *tframetrue = ws->var("ctau3Dtrue")->frame(Bins(150));
  RooPlot *tframetrue = ws->var("ctau3Dtrue")->frame();
  redNPMCTrue->plotOn(tframetrue);
  ws->pdf("CtNPTrue")->plotOn(tframetrue, LineColor(kOrange));
  tframetrue->GetXaxis()->CenterTitle(1);
  tframetrue->GetYaxis()->CenterTitle(1);

  char titlestr_lin[1000];
  TLatex t;
  t.SetNDC();
  t.SetTextAlign(12);
  t.SetTextSize(0.035);

  TCanvas c0f;
  c0f.cd();
  c0f.SetLogy(0);
  tframetrue->Draw();
  sprintf(titlestr_lin, "coefExpNPTrue: %.2f #pm %.1e", ws->var("coefExpNPTrue")->getVal(), ws->var("coefExpNPTrue")->getError());
  t.DrawLatex(0.43, 0.84, titlestr_lin);
  sprintf(titlestr_lin, "sigmaNPTrue: %.1e #pm %.1e", ws->var("sigmaNPTrue")->getVal(), ws->var("sigmaNPTrue")->getError());
  t.DrawLatex(0.43, 0.76, titlestr_lin);
  sprintf(titlestr_lin, "%s_CtNPTrue_Lin.pdf", titlestr.c_str());
  c0f.SaveAs(titlestr_lin);

  c0f.SetLogy(1);
  double truemax = tframetrue->GetMaximum();
  tframetrue->SetMinimum(0.5);
  tframetrue->SetMaximum(truemax * 5);
  sprintf(titlestr_lin, "%s_CtNPTrue_Log.pdf", titlestr.c_str());
  c0f.SaveAs(titlestr_lin);
  delete tframetrue;

  ws->var("sigmaNPTrue")->setConstant(kTRUE);

  return;
}

void defineCtNP(RooWorkspace *ws, RooDataSet *redNPMC, string titlestr)
{
  RooDataHist *binMCCutNP = new RooDataHist("binMCCutNP", "MC distribution for NP signal", RooArgSet(*(ws->var("ctau3Dtrue"))), *redNPMC);

  // simple ctau NP pdf
  ws->factory("GaussModel::BResTrue(ctau3Dtrue,mean[0.0],sigmaNPTrue[0.00002,0.0000001,0.1])"); // resolution from B meson (before CtPRRes)
  ws->factory("Decay::CtNPTrue(ctau3Dtrue,coefExpNPTrue[0.1,0.00001,10],BResTrue,RooDecay::SingleSided)");

  // fit to NP MC
  drawCtNPTrue(ws, redNPMC, titlestr);

  // --- build NP Pdf convoluted PRRes ---
  RooFormulaVar sigmaNPResW("sigmaNPResW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResW")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  ws->import(sigmaNPResW);
  RooFormulaVar sigmaNPResN("sigmaNPResN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResN")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  ws->import(sigmaNPResN);
  ws->factory("GaussModel::GW_NPRes(ctau3D,meanPRResW,sigmaNPResW)");
  ws->factory("GaussModel::GN_NPRes(ctau3D,meanPRResN,sigmaNPResN)");
  ws->factory("AddModel::CtNPRes({GW_NPRes,GN_NPRes},{fracRes})");
  
  // final model
  float coefExpNPTrueVal = ws->var("coefExpNPTrue")->getVal();
  ws->factory("Decay::CtNPTot(ctau3D,coefExpNPTrue,CtNPRes,RooDecay::SingleSided)");
  
  //// **** check NPMC Reco
  // 이거 왜 있지????
  // drawCtNPReco(ws, redNPMC, titlestr, opt);
  return;
}

// draw helper
void drawInclusiveMassPlots(RooWorkspace *ws, RooDataHist *redDataCut, RooFitResult *fitMass, bool isMC)
{

  //// **** Temporary variables for plotting
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  TLatex *ty = new TLatex();
  ty->SetNDC();
  ty->SetTextAlign(12);
  char reduceDS[1000];
  string titlestr;

  //// *** Mass plot
  RooBinning rb(ws->var("mass")->getBinning().numBins(), ws->var("mass")->getBinning().array());
  RooPlot *mframe_wob = ws->var("mass")->frame();
  redDataCut->plotOn(mframe_wob, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

  double avgBinWidth = rb.averageBinWidth();
  mframe_wob->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/c^{2 })", avgBinWidth * 1000));
  mframe_wob->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframe_wob->GetXaxis()->CenterTitle(1);
  mframe_wob->GetYaxis()->CenterTitle(1);
  const double max = mframe_wob->GetMaximum() * 1.3;
  mframe_wob->SetMaximum(max);
  mframe_wob->SetMinimum(0.1);

  // ws->pdf("MassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("MassPDF")->plotOn(mframe_wob, DrawOption("F"), FillColor(kGray + 2), FillStyle(3354), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // ws->pdf("MassPDF")->plotOn(mframe_wob,Components(opt.mBkgFunct.c_str()),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("MassPDF")->plotOn(mframe_wob, Components("expBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // ws->pdf("MassPDF")->plotOn(mframe_wob,Components("expBkg"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("MassPDF")->plotOn(mframe_wob, Components("expBkg"), LineColor(kBlue - 2), LineStyle(7), LineWidth(5), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), Name("bkg"));
  ws->pdf("MassPDF")->plotOn(mframe_wob, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  redDataCut->plotOn(mframe_wob, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

  TH1 *hdata = redDataCut->createHistogram("hdata", *ws->var("mass"), Binning(rb));
  // *** Calculate chi2/nDof for mass fitting
  int nBins = hdata->GetNbinsX();
  RooHist *hpullm;
  hpullm = mframe_wob->pullHist();
  hpullm->SetName("hpullM");
  double chi2 = 0;
  int nFullBinsPull = 0;
  double *ypull = hpullm->GetY();
  for (unsigned int i = 0; i < nBins; i++)
  {
    if (hdata->GetBinContent(i + 1) == 0)
      continue;
    chi2 += pow(ypull[i], 2);
    nFullBinsPull++;
  }
  double UnNormChi2 = chi2;
  int nFitParam = fitMass->floatParsFinal().getSize();
  int Dof = nFullBinsPull - nFitParam;
  chi2 /= (nFullBinsPull - nFitParam);

  TCanvas c1wop;
  c1wop.Draw();
  mframe_wob->SetTitleOffset(1.47, "Y");
  mframe_wob->Draw();
  // c1wop.SetLogy();

  ty->SetTextSize(0.040);
  sprintf(reduceDS, "N_{J/#psi} = %0.0f #pm %0.0f", ws->var("NSig")->getVal(), ws->var("NSig")->getError());
  ty->DrawLatex(0.20, 0.85, reduceDS);
  // sprintf(reduceDS, "#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
  ty->DrawLatex(0.20, 0.79, reduceDS);

  TLegend *leg11 = new TLegend(0.18, 0.51, 0.54, 0.72, NULL, "brNDC");
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0);
  leg11->SetShadowColor(0);
  leg11->SetTextSize(0.035);
  leg11->SetTextFont(42);
  leg11->SetMargin(0.2);
  // leg11->AddEntry(gfake1, "Data", "p");
  // leg11->AddEntry(&hfake21, "Total fit", "lf");
  // leg11->AddEntry(&hfake11, "Background", "lf");
  // leg11->AddEntry(mframe_wob->findObject("bkg"), "Background", "lf");
  leg11->Draw("same");

  if (isMC)
  {
    sprintf(reduceDS, "prompt MC");
    t->DrawLatex(0.91, 0.70, reduceDS);
  }
  sprintf(reduceDS, "#alpha_{CB}: %.2f #pm %.2f", ws->var("alpha")->getVal(), ws->var("alpha")->getError());
  t->DrawLatex(0.91, 0.64, reduceDS);
  sprintf(reduceDS, "n_{CB}: %.2f #pm %.2f", ws->var("enne")->getVal(), ws->var("enne")->getError());
  t->DrawLatex(0.91, 0.58, reduceDS);

  if (isMC)
    titlestr = "_massfitInclMC_wopull.pdf";
  else
    titlestr = "_massfitIncl_wopull.pdf";
  c1wop.SaveAs(titlestr.c_str());
  //  if (isMC) titlestr = opt.dirName + "_rap" + opt.yrange  + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfitInclMC_wopull.root";
  //  else titlestr = opt.dirName + "_rap" + opt.yrange  + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfitIncl_wopull.root";
  //  c1wop.SaveAs(titlestr.c_str());
}

void drawInclusiveMcMassPlots(RooWorkspace *ws, RooDataHist *redDataCut, RooFitResult *fitMass, bool isMC)
{

  //// **** Temporary variables for plotting
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  TLatex *ty = new TLatex();
  ty->SetNDC();
  ty->SetTextAlign(12);
  char reduceDS[1000];
  string titlestr;

  //// *** Mass plot
  RooBinning rb(ws->var("mass")->getBinning().numBins(), ws->var("mass")->getBinning().array());
  RooPlot *mframe_wob = ws->var("mass")->frame();
  redDataCut->plotOn(mframe_wob, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

  double avgBinWidth = rb.averageBinWidth();
  mframe_wob->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/c^{2 })", avgBinWidth * 1000));
  mframe_wob->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframe_wob->GetXaxis()->CenterTitle(1);
  mframe_wob->GetYaxis()->CenterTitle(1);
  const double max = mframe_wob->GetMaximum() * 1.3;
  mframe_wob->SetMaximum(max);
  mframe_wob->SetMinimum(0.1);

  // ws->pdf("MassPDF")->plotOn(mframe_wob,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("MassMCPDF")->plotOn(mframe_wob, DrawOption("F"), FillColor(kGray + 2), FillStyle(3354), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  ws->pdf("MassMCPDF")->plotOn(mframe_wob, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  redDataCut->plotOn(mframe_wob, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

  TH1 *hdataMC = redDataCut->createHistogram("hdataMC", *ws->var("mass"), Binning(rb));
  // *** Calculate chi2/nDof for mass fitting
  int nBins = hdataMC->GetNbinsX();
  RooHist *hpullm;
  hpullm = mframe_wob->pullHist();
  hpullm->SetName("hpullM");
  double chi2 = 0;
  int nFullBinsPull = 0;
  double *ypull = hpullm->GetY();
  for (unsigned int i = 0; i < nBins; i++)
  {
    if (hdataMC->GetBinContent(i + 1) == 0)
      continue;
    chi2 += pow(ypull[i], 2);
    nFullBinsPull++;
  }
  double UnNormChi2 = chi2;
  int nFitParam = fitMass->floatParsFinal().getSize();
  int Dof = nFullBinsPull - nFitParam;
  chi2 /= (nFullBinsPull - nFitParam);

  TCanvas c1wop;
  c1wop.Draw();
  mframe_wob->SetTitleOffset(1.47, "Y");
  mframe_wob->Draw();
  // c1wop.SetLogy();

  ty->SetTextSize(0.040);
  sprintf(reduceDS, "N_{J/#psi} = %0.0f #pm %0.0f", ws->var("NSig")->getVal(), ws->var("NSig")->getError());
  ty->DrawLatex(0.20, 0.85, reduceDS);
  // sprintf(reduceDS, "#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
  ty->DrawLatex(0.20, 0.79, reduceDS);

  TLegend *leg11 = new TLegend(0.18, 0.51, 0.54, 0.72, NULL, "brNDC");
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0);
  leg11->SetShadowColor(0);
  leg11->SetTextSize(0.035);
  leg11->SetTextFont(42);
  leg11->SetMargin(0.2);
  // leg11->AddEntry(gfake1, "Data", "p");
  // leg11->AddEntry(&hfake21, "Total fit", "lf");
  // leg11->AddEntry(&hfake11, "Background", "lf");
  leg11->Draw("same");

  if (isMC)
  {
    sprintf(reduceDS, "prompt MC");
    t->DrawLatex(0.91, 0.70, reduceDS);
  }
  sprintf(reduceDS, "#alpha_{CB}: %.2f #pm %.2f", ws->var("alpha")->getVal(), ws->var("alpha")->getError());
  t->DrawLatex(0.91, 0.64, reduceDS);
  sprintf(reduceDS, "n_{CB}: %.2f #pm %.2f", ws->var("enne")->getVal(), ws->var("enne")->getError());
  t->DrawLatex(0.91, 0.58, reduceDS);

  if (isMC)
    titlestr = "_massfitInclMC_wopull.pdf";
  else
    titlestr = "_massfitIncl_wopull.pdf";
  c1wop.SaveAs(titlestr.c_str());
  //  if (isMC) titlestr = opt.dirName + "_rap" + opt.yrange  + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfitInclMC_wopull.root";
  //  else titlestr = opt.dirName + "_rap" + opt.yrange  + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfitIncl_wopull.root";
  //  c1wop.SaveAs(titlestr.c_str());
}

void drawCtauErrPdf(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG)
{
  RooPlot *errframe = ws->var("ctau3DErr")->frame();

  binDataCtErrSB->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kBlue), MarkerColor(kBlue), MarkerStyle(kOpenCircle));
  ws->pdf("errPdfBkg")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binDataCtErrSB->sumEntries(), RooAbsReal::NumEvent));
  // ws->pdf("errPdfBkg")->plotOn(errframe);

  TCanvas c0;
  string titlestr;
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = "_CtErrPdfBkg_Log.pdf";
  c0.SaveAs(titlestr.c_str());
  delete errframe;

  errframe = ws->var("ctau3DErr")->frame();
  binSubtractedSIG->plotOn(errframe, DataError(RooAbsData::SumW2), DrawOption("F"), FillColor(kWhite), LineColor(kWhite)); // just for axis
  ws->pdf("errPdfSig")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binSubtractedSIG->sumEntries(), RooAbsReal::NumEvent));
  binDataCtErrSIG->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kRed), MarkerColor(kRed), MarkerStyle(kOpenCircle));
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = "_CtErrPdfSig_Log.pdf";
  c0.SaveAs(titlestr.c_str());
  delete errframe;
}



void drawCtauResolPlots(RooWorkspace *ws, bool fitMC, RooPlot *tframePR)
{
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  char reduceDS[1000];
  string titlestr;

  t->SetTextSize(0.035);
  TCanvas c00;
  c00.cd();
  tframePR->Draw();
  sprintf(reduceDS, "#sigma(G_{N}): %.2f", ws->var("sigmaPRResN")->getVal());
  t->DrawLatex(0.55, 0.31, reduceDS);
  sprintf(reduceDS, "#sigma(G_{W}): %.2f", ws->var("sigmaPRResW")->getVal());
  t->DrawLatex(0.55, 0.26, reduceDS);
  sprintf(reduceDS, "frac(G_{W}): %.2f", ws->var("fracRes")->getVal());
  t->DrawLatex(0.55, 0.21, reduceDS);
  tframePR->GetXaxis()->CenterTitle(1);
  tframePR->GetYaxis()->CenterTitle(1);
  if (fitMC)
    titlestr = "_CtPRResMC_Lin.pdf";
  else
    titlestr = "_CtPRResData_Lin.pdf";
  c00.SaveAs(titlestr.c_str());

  c00.SetLogy(1);
  tframePR->GetXaxis()->CenterTitle(1);
  tframePR->GetYaxis()->CenterTitle(1);
  double prmax = tframePR->GetMaximum();
  tframePR->SetMinimum(0.5);
  tframePR->SetMaximum(prmax * 5);
  if (fitMC)
    titlestr = "_CtPRResMC_Log.pdf";
  else
    titlestr = "_CtPRResData_Log.pdf";
  c00.SaveAs(titlestr.c_str());
}



void drawFinalMass(RooWorkspace *ws, RooDataHist* redDataCut, float NSigNP_fin, float NBkg_fin, RooFitResult *fitMass, double* UnNormChi2_mass_t, int* nFitParam_mass_t, int* nFullBinsPull_mass_t, int* Dof_mass_t,double* Chi2_mass_t) {
  //// **** Temporary variables for plotting
  TLatex *t = new TLatex();  t->SetNDC();  t->SetTextAlign(32);
  TLatex *ty = new TLatex();  ty->SetNDC();  ty->SetTextAlign(12);
  char reduceDS[1000];
  string titlestr;
  
  // RooBinning rb(ws->var("mass")->getBinning().numBins(), ws->var("mass")->getBinning().array()); // 14-009
  RooBinning rb(300, 2.6, 3.5);
  RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
  RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);

  //// *** Mass plot
  RooPlot *mframe = ws->var("mass")->frame();
  redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rb));
  double avgBinWidth = rb.averageBinWidth();
  mframe->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/c^{2 })",avgBinWidth*1000));
  mframe->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
//  mframe->GetXaxis()->SetTitleSize(0.048*1.2); //PAPER
//  mframe->GetYaxis()->SetTitleSize(0.048*1.2); //PAPER
  mframe->GetXaxis()->CenterTitle(1);
  mframe->GetYaxis()->CenterTitle(1);
  const double max = mframe->GetMaximum() * 1.3;
  mframe->SetMaximum(max);
  mframe->SetMinimum(0);

    //// **** Fill color
    // G1CB1Sig, expBkg
    //ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kGray+2),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent), NumCPU(32));
    RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf("G1CB1Sig")),*(ws->pdf("expBkg"))),RooArgList(tmpVar1,tmpVar2));
    //tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kWhite), FillStyle(1001), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
    //tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
    gStyle->SetHatchesLineWidth(2);
    //ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
    //// **** Line color
    //ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
    ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), LineColor(kBlue - 2), LineStyle(7), LineWidth(5), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
    //tmpPDF.plotOn(mframe,LineColor(kRed),LineStyle(9),LineWidth(5),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
    tmpPDF.plotOn(mframe, LineColor(kPink - 6), LineStyle(11), LineWidth(5), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
    ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
    redDataCut->plotOn(mframe, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

    TH1 *hdata = redDataCut->createHistogram("hdata", *ws->var("mass"), Binning(rb));
    //// *** Calculate chi2/nDof for mass fitting
    int nBins = hdata->GetNbinsX();
    RooHist *hpullm;
    hpullm = mframe->pullHist();
    hpullm->SetName("hpullM");
    double Chi2 = 0;
    int nFullBinsPull = 0;
    double *ypull = hpullm->GetY();
    for (unsigned int i = 0; i < nBins; i++)
    {
      if (hdata->GetBinContent(i + 1) == 0)
        continue;
      nFullBinsPull++;
      Chi2 = Chi2 + pow(ypull[i], 2);
    }
  double UnNormChi2 = Chi2;
  *UnNormChi2_mass_t = Chi2;
  int nFitParam = fitMass->floatParsFinal().getSize();
  int Dof = nFullBinsPull - nFitParam;
  Chi2 /= (nFullBinsPull - nFitParam);
  *nFitParam_mass_t = nFitParam;
  *nFullBinsPull_mass_t = nFullBinsPull;
  *Dof_mass_t = Dof;
  *Chi2_mass_t = Chi2;

  // PAPER
  TCanvas* c1wop = new TCanvas("c1wop","The Canvas",200,10,600,600);
  c1wop->cd(); c1wop->Draw(); 
  
  mframe->SetTitleOffset(1.47,"Y");
  mframe->Draw();
  //// **** different lumiTextOffset for massfit_wopull
  // lumiTextOffset   = 0.20;
  // CMS_lumi(c1wop, opt.isPA, iPosPaper);
  t->SetTextSize(0.035);
//  t->DrawLatex(0.91,0.85,opt.rapString);
//  t->DrawLatex(0.91,0.78,opt.ptString);
  ty->SetTextSize(0.035); // PAPER
  // ty->DrawLatex(0.20,0.86,opt.rapString); //PAPER
  // ty->DrawLatex(0.20,0.80,opt.ptString); //PAPER
  ty->SetTextSize(0.040);
  sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NSig")->getVal(),ws->var("NSig")->getError());
//  ty->DrawLatex(0.20,0.85,reduceDS);
  // sprintf(reduceDS,"#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
//  ty->DrawLatex(0.20,0.79,reduceDS);

  //TLegend * legpaper = new TLegend(0.59,0.53,0.90,0.72,NULL,"brNDC");
  //TLegend * legpaper = new TLegend(0.17,0.53,0.63,0.72,NULL,"brNDC");
  TLegend * legpaper = new TLegend(0.17,0.55,0.57,0.75,NULL,"brNDC");
  legpaper->SetFillStyle(0); legpaper->SetBorderSize(0); legpaper->SetShadowColor(0);
  legpaper->SetTextSize(0.035); legpaper->SetTextFont(42);
  legpaper->SetMargin(0.2);
  // legpaper->AddEntry(gfake1,"Data","p");
  // legpaper->AddEntry(&hfake21,"Total fit","lf");
  // legpaper->AddEntry(&hfake31,"Bkg + nonprompt","lf"); 
  // legpaper->AddEntry(&hfake11,"Background","lf");
  legpaper->Draw("same");

  titlestr = "_massfit_wopull.pdf";
  c1wop->SaveAs(titlestr.c_str());
//  titlestr = opt.dirName + "_rap" + opt.yrange + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfit_wopull.root";
//  c1wop.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////  

  TCanvas c1("c1","The mass Canvas",200,10,600,750);
  c1.cd();
  TPad *padm1 = new TPad("padm1","This is pad1",0.0,0.3,1.0,1.0);
  padm1->SetLeftMargin(0.14);
  padm1->SetRightMargin(0.03);
  padm1->SetTopMargin(0.075);
  padm1->SetBottomMargin(0);  padm1->Draw();
  TPad *padm2 = new TPad("padm2","This is pad2",0.00,0.00,1.0,0.3);
  padm2->SetLeftMargin(0.14);
  padm2->SetRightMargin(0.03);
  padm2->SetTopMargin(0);  
  padm2->SetBottomMargin(0.30);  padm2->Draw();

  padm1->cd();  mframe->Draw();
  // lumiTextOffset   = 0.45;
  // CMS_lumi(&c1, opt.isPA, iPos);
  // t->SetTextSize(0.035);
  // t->DrawLatex(0.91,0.90,opt.rapString);
  // t->DrawLatex(0.91,0.85,opt.ptString);
  // if (opt.EventActivity ==1) t->DrawLatex(0.91,0.80,opt.ntrkString);
  // else if (opt.EventActivity ==2) t->DrawLatex(0.91,0.80,opt.etString);
  
  ty->SetTextSize(0.040);
  // sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NSig")->getVal(),ws->var("NSig")->getError());
  ty->DrawLatex(0.20,0.89,reduceDS);
  // sprintf(reduceDS,"#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
  ty->DrawLatex(0.20,0.84,reduceDS);
  
  TLegend * leg11 = new TLegend(0.18,0.51,0.54,0.72,NULL,"brNDC");
  leg11->SetFillStyle(0); leg11->SetBorderSize(0); leg11->SetShadowColor(0);
  leg11->SetTextSize(0.035); leg11->SetTextFont(42);
  leg11->SetMargin(0.2);
  // leg11->AddEntry(gfake1,"Data","p");
  // leg11->AddEntry(&hfake21,"Total fit","lf");
  // leg11->AddEntry(&hfake31,"Bkg + nonprompt","lf"); 
  // leg11->AddEntry(&hfake11,"Background","lf");
  leg11->Draw("same");
  c1.Update();

  //// **** pull
  RooPlot* mframepull =  ws->var("mass")->frame(Title("Pull")) ;
  mframepull->GetYaxis()->SetTitle("Pull");
  mframepull->GetYaxis()->CenterTitle(1);
  mframepull->SetLabelSize(0.04*2.5,"XYZ");
  mframepull->SetTitleSize(0.048*2.5,"XYZ");
  mframepull->SetTitleOffset(0.47,"Y");
  mframepull->addPlotable(hpullm,"PX") ;
  double mframemax = 0;
  if (mframepull->GetMinimum()*-1 > mframepull->GetMaximum()) mframemax = mframepull->GetMinimum()*-1;
  else mframemax = mframepull->GetMaximum();
  mframepull->SetMaximum(mframemax); 
  mframepull->SetMinimum(-1*mframemax); 
  mframepull->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframepull->GetXaxis()->CenterTitle(1);

  padm2->cd(); mframepull->Draw();
  TLine* line1 = new TLine(2.6,0,3.5,0.); line1->SetLineStyle(7); line1->Draw();

  TLatex *t2 = new TLatex();
  t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.035*3);
  sprintf(reduceDS,"#chi^{2}/dof = %.1f/%d",UnNormChi2,Dof);
  t2->DrawLatex(0.78,0.86,reduceDS);
  c1.Update();

  titlestr = "_massfit.pdf";
  c1.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////  
  mframe->SetMinimum(0.5);
  mframe->SetMaximum(max*50);
  padm1->cd();
  padm1->SetLogy(1);
  titlestr = "_massfit_Log.pdf";
  c1.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////  
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void drawFinalCtau(RooWorkspace *ws, RooDataHist *redDataCut, RooDataHist* binDataCtErr, float NSigNP_fin, float NBkg_fin, float Bfrac_fin, float ErrBfrac_fin, RooFitResult *fit2D, float lmin, float lmax, double* UnNormChi2_time_t, int* nFitParam_time_t, int* nFullBinsPull_time_t, int* Dof_time_t,double* Chi2_time_t) {
  char reduceDS[1000];
  string titlestr;
  
  // RooBinning rb(ws->var("ctau3D")->getBinning().numBins(), ws->var("ctau3D")->getBinning().array()); // 14-009
  RooBinning rb(300, -0.5, 2.0);

  RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
  RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);
  
  RooPlot *tframe = ws->var("ctau3D")->frame();
  tframe->SetTitleOffset(1.47,"Y");
  double avgBinWidth = rb.averageBinWidth();
  //tframe->GetYaxis()->SetTitle(Form("Counts / %.2f (mm)",avgBinWidth));
  tframe->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)",avgBinWidth*1000));
  tframe->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframe->GetXaxis()->CenterTitle(1);
  tframe->GetYaxis()->CenterTitle(1);

  //// **** Ctau total distributions
  RooHist *hpulltot;
  redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb),MarkerSize(1));

    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent));
    hpulltot = tframe->pullHist(); hpulltot->SetName("hpulltot");
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtBkg"),LineColor(kBlue-2),LineWidth(5),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(7));
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtNP"),LineColor(kPink-6),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(12));
    ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtPR"),LineColor(kGreen+3),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(kDashDotted));
    ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent));

  TH1 *hdatact = redDataCut->createHistogram("hdatact",*ws->var("ctau3D"),Binning(rb));
  double chi2 = 0, unNormChi2 = 0;
  int dof = 0;
  double *ypulls = hpulltot->GetY();
  unsigned int nBins = ws->var("ctau3D")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++) {
    if (hdatact->GetBinContent(i+1) == 0) continue;
    chi2 += ypulls[i]*ypulls[i];
    nFullBins++;
  }
  unNormChi2 = chi2;
  *UnNormChi2_time_t = chi2;
  int nFitPar = fit2D->floatParsFinal().getSize();
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  *nFitParam_time_t = nFitPar;
  *nFullBinsPull_time_t = nFullBins;
  *Dof_time_t =dof;
  *Chi2_time_t = chi2;
  
  TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,750);
  c2->cd();
  TPad *pad1 = new TPad("pad1","This is pad1",0.0,0.3,1.0,1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.075);
  pad1->SetBottomMargin(0);  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.0,0.0,1.0,0.3);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0);  
  pad2->SetBottomMargin(0.30);  pad2->Draw();

  pad1->cd(); tframe->Draw();
  
  // lumiTextOffset   = 0.45;
  // CMS_lumi(c2, opt.isPA, iPos);
  TLatex *ty = new TLatex();  ty->SetNDC();  ty->SetTextAlign(32);
  ty->SetTextSize(0.035);
  // ty->DrawLatex(0.91,0.90,opt.rapString);
  // ty->DrawLatex(0.91,0.85,opt.ptString);
  // if (opt.EventActivity ==1) ty->DrawLatex(0.91,0.80,opt.ntrkString);
  // else if (opt.EventActivity ==2) ty->DrawLatex(0.91,0.80,opt.etString);

  //TLegend * leg = new TLegend(0.66,0.56,0.85,0.75,NULL,"brNDC");
  TLegend * leg = new TLegend(0.66,0.56,0.99,0.75,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetShadowColor(0);
  leg->SetTextSize(0.035); leg->SetTextFont(42);
  leg->SetMargin(0.2);
  // leg->AddEntry(gfake1,"Data","p");
  // leg->AddEntry(&hfake21,"Total fit","l");
  // leg->AddEntry(&hfake41,"Prompt","l"); 
  // leg->AddEntry(&hfake311,"Nonprompt","l"); 
  // leg->AddEntry(&hfake11,"Background","l");
  leg->Draw("same"); 
 
  // KYO : write down Bfrac on plot 
  TLatex *tbfrac = new TLatex();
  tbfrac->SetNDC(); tbfrac->SetTextAlign(12); tbfrac->SetTextSize(0.035);
  sprintf(reduceDS,"B frac. = %.2f #pm %.2f",Bfrac_fin,ErrBfrac_fin);
  tbfrac->DrawLatex(0.19,0.91,reduceDS);

  RooPlot* tframepull =  ws->var("ctau3D")->frame(Title("Pull")) ;
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->GetYaxis()->CenterTitle(1);
  tframepull->SetLabelSize(0.04*2.5,"XYZ");
  tframepull->SetTitleSize(0.048*2.5,"XYZ");
  tframepull->SetTitleOffset(0.47,"Y");
  tframepull->addPlotable(hpulltot,"PX") ;
  double tframemax = 0;
  if (tframepull->GetMinimum()*-1 > tframepull->GetMaximum()) tframemax = tframepull->GetMinimum()*-1;
  else tframemax = tframepull->GetMaximum();
  tframepull->SetMaximum(tframemax); 
  tframepull->SetMinimum(-1*tframemax);
  tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframepull->GetXaxis()->CenterTitle(1);

  pad2->cd(); tframepull->Draw();
  TLine* line1 = new TLine(-lmin,0,lmax,0.); line1->SetLineStyle(7); line1->Draw();

  int nDOF = ws->var("ctau3D")->getBinning().numBins() - nFitPar;

  TLatex *t2 = new TLatex();
  t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.035*3);
  sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
  t2->DrawLatex(0.78,0.90,reduceDS);
  
  c2->Update();
  titlestr = "_timefit_Lin.pdf";
  c2->SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////
  tframe->SetMaximum(tframe->GetMaximum()*9); 
  tframe->SetMinimum(0.5); 
  pad1->SetLogy(1);
  titlestr = "_timefit_Log.pdf";
  c2->SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////
}

void fit2D_622()
{
  cout << "=== start fit2D_622() ===\n";
  TStopwatch time;
  time.Start();

  // === Kinematics and observables ranges ===
  // --- kinematics ---
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // --- observable ranges
  double ctMin = -0.5, ctMax = 2; // lmin, lmax: 1 for lowpT, 2 for higpT
  float errmin = 0.008, errmax = 0.3;

  // --- silent ---
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only print from WARNING to FATAL

  // --- rootlogon ---
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // Read MC and Data files
  // ----------------------
  cout << "=== Import inputs ===\n";

  // --- Prompt MC ---
  string fileNamePrMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";
  TFile fInPRMC(fileNamePrMc.c_str());
  cout << fileNamePrMc.c_str() << endl;
  // RooDataSet *dataPRMC = (RooDataSet *)fInPRMC.Get("dataset");

  RooDataSet *dataPRMC = (RooDataSet *)fInPRMC.Get("dataset");
  dataPRMC->SetName("dataPRMC");

  // --- NP MC ---
  string fileNameNpMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_y0_2p4_230829.root";
  TFile fInNPMC(fileNameNpMc.c_str());
  cout << fileNameNpMc.c_str() << endl;
  RooDataSet *dataNPMC = (RooDataSet *)fInNPMC.Get("dataset");
  dataNPMC->SetName("dataNPMC");

  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");

  // --- Created workspace ---
  RooWorkspace *ws = new RooWorkspace("workspace");

  // === Reduce data to reduceDS ===
  // declare ds
  cout << "\n=== Make cuts ===\n";
  char reduceDS[3000], reduceDS_woCtErr[3000], reduceNpMc[3000];
  // reduceNpMc_woCtErr[3000]; -> dataest 만들 때 안 쓰이니까 삭제

  sprintf(reduceDS_woCtErr, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"
          
          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f)"
          
          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);

  getCtErrRange(data,reduceDS_woCtErr, -ctMin, ctMax, &errmin, &errmax);

  // with CtErrCut
  sprintf(reduceDS, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"
          
          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"
          
          "&& (recoQQsign == 0)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errmin, errmax);

  sprintf(reduceNpMc,
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3Dtrue >= %.3f && ctau3Dtrue < %.3f)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, 0.0001, ctMax);

  cout << "reduceDS: " << reduceDS << endl;
  cout << "reduceNpMc: " << reduceNpMc << endl;

  RooDataSet *redPRMC, *redNPMC, *redData, *redData_woCtErr; // with CtErr cut  for for isPEE==1
  redPRMC = (RooDataSet *)dataPRMC->reduce(reduceDS); // smae cut with Data
  redNPMC = (RooDataSet *)dataNPMC->reduce(reduceNpMc);
  redData = (RooDataSet *)data->reduce(reduceDS);
  redData_woCtErr = (RooDataSet *)data->reduce(reduceDS_woCtErr);

  ws->import(*redPRMC);
  ws->import(*redNPMC);
  ws->import(*redData);

  setWSRange(ws, -ctMin, ctMax, errmin, errmax);


  // === Binning ===
  // global variables
  string titlestr;
  TCanvas c0;

  // variable title
  ws->var("mass")->SetTitle("m_{#mu#mu}");
  ws->var("ctau3D")->SetTitle("#font[12]{l}_{J/#psi}");

  // mass bin
  ws->var("mass")->setBins(120); // default 45: bin/0.02 GeV

  // ct bin
  // refer to setCtBinning() -> different range and number according to range
  // RooBinning rbct = setCtBinning(-ctMin, ctMax);
  // ws->var("ctau3D")->setBinning(rbct);
  ws->var("ctau3D")->setBins(200);

  // ctTrue binning
  // RooBinning rbtrue(-0.1, 4.0);
  // rbtrue.addUniform(5, -0.1, 0.0);
  // rbtrue.addUniform(100, 0.0, 0.5);
  // rbtrue.addUniform(15, 0.5, 1.0);
  // rbtrue.addUniform(20, 1.0, 2.5);
  // rbtrue.addUniform(5, 2.5, 4.0);
  // ws->var("ctau3Dtrue")->setBinning(rbtrue);
  ws->var("ctau3Dtrue")->setBins(100); // default

  // ctErr binning
  ws->var("ctau3DErr")->setBins(50);


  // --- Make subrange Dataset ---
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");


  // --- Define PDFs with params (mass and ctau) ---
  // mass
  defineMassBkg(ws);
  defineMassSig(ws);

  // ctau
  defineCtPRRes(ws); // R(l) : resolution function
  defineCtBkg(ws);   // theta(l') convolution R(l')
  titlestr = "any_tile_is_good";


  // --- Build NP ctau PDF ---
  defineCtNP(ws, redNPMC, titlestr); // F_B(l) : X_mc(l')


  // --- Build mass PDF ---
  char funct[200];
  double initBkg = redDataSB->sumEntries() * 9.0 / 5.0;
  double initSig = redData->sumEntries() - initBkg;
  initBkg = initBkg + 10000;
  ////// *** mSigFunct + mBkgFunct
  sprintf(funct,"SUM::MassPDF(NSig[%f,1.0,5000000.0]*%s, NBkg[%f,1.0,50000000.0]*%s)", initSig, "G1CB1Sig", initBkg, "expBkg");
  ws->factory(funct);

  // make MC PR Signal only PDF
  sprintf(funct, "SUM::MassMCPDF(NSig[%f,1.0,5000000.0]*%s)", initSig, "G1CB1Sig");
  ws->factory(funct);

  // --- MC mass fit ---
  RooFitResult *fitMassMC;

  RooRealVar *mass = ws->var("mass");
  RooRealVar *ctau3D = ws->var("ctau3D");
  RooRealVar *ctau3DErr = ws->var("ctau3DErr");
  mass->setBins(120);
  ctau3D->setBins(50);
  ctau3DErr->setBins(50);
  auto dhPRMC = new RooDataHist("dhPRMC", "", RooArgSet(*mass, *ctau3D, *ctau3DErr), *redPRMC);

  fitMassMC = ws->pdf("MassMCPDF")->fitTo(*dhPRMC, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(32));
  // fitMassMC = ws->pdf("MassMCPDF")->fitTo(*redPRMC, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(32));
  fitMassMC->Print("v");

  drawInclusiveMcMassPlots(ws, dhPRMC, fitMassMC, 1);

  ws->var("alpha")->setConstant(kFALSE);
  // ws->var("enne")->setVal(2.1); // fix???
  ws->var("enne")->setConstant(kTRUE);
  ws->var("coefPol")->setRange(-5., 5.);
  ws->var("coefPol")->setVal(-0.05);
  ws->var("coefPol")->setConstant(kFALSE);

  // --- Mass fit ---
  RooFitResult *fitMass;
  auto dhData = new RooDataHist("dhData", "", RooArgSet(*mass, *ctau3D, *ctau3DErr), *redData);

  fitMass = ws->pdf("MassPDF")->fitTo(*dhData, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(32));

  // fitMass = ws->pdf("MassPDF")->fitTo(*redData, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(32));
  // , EvalBackend("legacy")
  fitMass->Print("v");
  ws->var("alpha")->setConstant(kTRUE);
  ws->var("enne")->setConstant(kTRUE);
  ws->var("fracG1")->setConstant(kTRUE);
  ws->var("sigmaSig1")->setConstant(kTRUE);
  ws->var("sigmaSig2")->setConstant(kTRUE);
  ws->var("meanSig")->setConstant(kTRUE);
  ws->var("coefExp")->setConstant(kTRUE);
  ws->var("coefPol")->setConstant(kTRUE);
  ws->var("NSig")->setConstant(kTRUE);
  ws->var("NBkg")->setConstant(kTRUE);

  // combinedWidth -> Print 할 때 sigma를 하나만 보여주는 목적? -> 둘 다 보여주면 되니까 삭제
  // inOpt.PcombinedWidth = inOpt.combinedWidth * 1000;

  drawInclusiveMassPlots(ws, dhData, fitMass, 0);

  Double_t NSig_fin = ws->var("NSig")->getVal();
  Double_t ErrNSig_fin = ws->var("NSig")->getError();
  Double_t NBkg_fin = ws->var("NBkg")->getVal();
  Double_t ErrNBkg_fin = ws->var("NBkg")->getError();


  // --- Get ctauErr PDF for PEE ---
  // scaleF to scale down ctErr distribution in 2.9-3.3 GeV/c2
  float bc;
  bc = ws->var("coefExp")->getVal(); // expBkg
  float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  // RooDataSet(unbinned) to RooDataHist (binned)
  RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr", "binDataCtErr", RooArgSet(*(ws->var("ctau3DErr"))), *redData);
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  // extract Sig and Bkg: (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  RooDataHist *binSubtractedSIG, *binScaledBKG;
  binSubtractedSIG = new RooDataHist("binSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  binScaledBKG = subtractSidebands(ws, binSubtractedSIG, binDataCtErrSIG, binDataCtErrSB, scaleF, "ctau3DErr");

  // error PDF
  RooHistPdf errPdfSig("errPdfSig", "Error PDF signal", RooArgSet(*(ws->var("ctau3DErr"))), *binSubtractedSIG);
  ws->import(errPdfSig);
  //  RooHistPdf errPdfBkgRaw("errPdfBkg","Error PDF bkg before scaling",RooArgSet(*(ws->var("ctau3DErr"))),*binDataCtErrSB);  ws->import(errPdfBkg);
  RooHistPdf errPdfBkg("errPdfBkg", "Error PDF bkg scaled", RooArgSet(*(ws->var("ctau3DErr"))), *binScaledBKG);
  ws->import(errPdfBkg);

  //// **** Draw CtErr PDF
  auto *ct = ws->var("ctau3D");
  auto *err = ws->var("ctau3DErr");
  auto *Perr = ws->pdf("errPdfBkg"); // RooHistPdf: P(err)

  // === 개별 Decay pdf ===
  auto *CtPRRes = ws->pdf("CtPRRes");
  auto *CtBkgPos = ws->pdf("CtBkgPos");
  auto *CtBkgPos2 = ws->pdf("CtBkgPos2");
  auto *CtBkgNeg = ws->pdf("CtBkgNeg");
  auto *CtBkgNeg2 = ws->pdf("CtBkgNeg2");
  auto *CtBkgDbl = ws->pdf("CtBkgDbl");
  auto *CtBkgDbl2 = ws->pdf("CtBkgDbl2");

  // === fraction 변수 자동 생성 ===
  for (int i = 1; i <= 6; ++i)
  {
    TString name = Form("fracCtBkg%d", i);
    if (!ws->var(name))
    {
      RooRealVar f(name, name, 0.2, 0.01, 0.99);
      ws->import(f);
    }
  }

  // alias
  auto *f1 = ws->var("fracCtBkg1");
  auto *f2 = ws->var("fracCtBkg2");
  auto *f3 = ws->var("fracCtBkg3");
  auto *f4 = ws->var("fracCtBkg4");
  auto *f5 = ws->var("fracCtBkg5");
  auto *f6 = ws->var("fracCtBkg6");

  // === 1) 각 성분에 P(err) 곱 (ct만 조건) ===
  auto makePEE = [&](const char *name, RooAbsPdf *src)
  {
    RooProdPdf tmp(Form("%s_PEE", name), Form("%s * P(err)", name),
                   RooArgSet(*Perr),
                   RooFit::Conditional(*src, RooArgList(*ct)));
    ws->import(tmp);
  };

  makePEE("CtPRRes", CtPRRes);
  makePEE("CtBkgPos", CtBkgPos);
  makePEE("CtBkgPos2", CtBkgPos2);
  makePEE("CtBkgNeg", CtBkgNeg);
  makePEE("CtBkgNeg2", CtBkgNeg2);
  makePEE("CtBkgDbl", CtBkgDbl);
  makePEE("CtBkgDbl2", CtBkgDbl2);

  RooAddPdf CtBkgTot_PEE("CtBkgTot_PEE", "",
                         {*ws->pdf("CtPRRes_PEE"), *ws->pdf("CtBkgNeg_PEE"), *ws->pdf("CtBkgDbl_PEE"), *ws->pdf("CtBkgPos_PEE"), *ws->pdf("CtBkgDbl2_PEE"), *ws->pdf("CtBkgNeg2_PEE"), *ws->pdf("CtBkgPos2_PEE")},
                         {*f1, *f2, *f3, *f4, *f5, *f6}, true);
  ws->import(CtBkgTot_PEE);
  // CtPRRes_PEE, CtBkgPos_PEE
  // Preare total 2D PDFs with PEE
  // e.g  e.g.) AAA_PEE : AAA x errPdf
  // -----------------------------

  // Conditional to select Ct and Mass only :: for step[[4]] and step[[5]]
  // RooProdPdf CtBkgTot_PEE("CtBkgTot_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
  //                         Conditional(*(ws->pdf("CtBkgTot")), RooArgList(*(ws->var("ctau3D")))));
  // ws->import(CtBkgTot_PEE);

  RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE", *(ws->pdf("errPdfSig")),
                      Conditional(*(ws->pdf("CtPRRes")), RooArgList(*(ws->var("ctau3D")))));
  ws->import(CtPR_PEE);


  
  // // Build 2D PDF (mass x ctau) :: for step[[6]]
  // sprintf(funct, "PROD::MassCtPR(%s,CtPRRes)", "G1CB1Sig"); ws->factory(funct);
  // sprintf(funct, "PROD::MassCtNP(%s,CtNPTot)", "G1CB1Sig"); ws->factory(funct);
  // sprintf(funct, "PROD::MassCtBkg(%s,CtBkgTot)", "expBkg"); ws->factory(funct);

  // RooProdPdf MassCtPR_PEE("MassCtPR_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
  //                         Conditional(*(ws->pdf("MassCtPR")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  // ws->import(MassCtPR_PEE);
  // RooProdPdf MassCtNP_PEE("MassCtNP_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
  //                         Conditional(*(ws->pdf("MassCtNP")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  // ws->import(MassCtNP_PEE);
  // RooProdPdf MassCtBkg_PEE("MassCtBkg_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
  //                          Conditional(*(ws->pdf("MassCtBkg")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  // ws->import(MassCtBkg_PEE);

  // Total fit PDF = promptPDF + nonpromptPDF + bkgPDF
  // 앞의 둘은 mass 시그널에서 개수를 정하고 마지막은 mass bkg에서 개수를 정한다.
  // RooFormulaVar fracBkg("fracBkg", "@0/(@0+@1)", RooArgList(*(ws->var("NBkg")), *(ws->var("NSig"))));
  // ws->import(fracBkg);
  // ws->factory("RSUM::totPDF_PEE(fracBkg*MassCtBkg_PEE,Bfrac[0.25,0.01, 0.99]*MassCtNP_PEE,MassCtPR_PEE)");
  // ws->factory("SUM::totPDF_PEE(fracBkg*MassCtBkg_PEE, Bfrac[0.25,0.0,1.]*MassCtNP_PEE, MassCtPR_PEE)");

  // fit prompt MC ctau Resolution
  // -----------------------------
  // RooFitResult *fitCt_PRMC = ws->pdf("CtPR_PEE")->fitTo(*dhPRMC, Range("promptMCfit"), SumW2Error(kTRUE), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), Save(1), NumCPU(32));
  // // EvalBackend("legacy")
  // fitCt_PRMC->Print("v");
  
  // ws->var("meanPRResW")->setConstant(kTRUE);
  
  // ws->var("meanPRResN")->setConstant(kTRUE);
  
  

  // // --- Check goodness of promptMCfit with per event error fit //// CtWeighted = lxy/(lxy_err) ---
  // // fit은 앞에서 했고 문제 없는지 그림 그려서 보기
  // // Res 변수 추가하는 과정
  // // 요즘에는 그냥 ctau3DRes라는 이름으로 따로 저장해놓는다. -> 일단 오리지널 코드 방식으로
  // RooRealVar *CtWeighted = new RooRealVar("CtWeighted", "#font[12]{l}_{J/#psi} / #sigma( #font[12]{l}_{J/#psi} )", -5., 5.);
  // ws->import(*CtWeighted);
  // const RooArgSet *thisRow = (RooArgSet *)dhPRMC->get(0); // prompt MC rooData
  // RooArgSet *newRow = new RooArgSet(*CtWeighted);
  // RooDataSet *tempSet = new RooDataSet("tempSet", "new data set with CtWeighted", *newRow);

  // // Per-event-error를 생으로 곱하고 있다.
  // // Res 함수 정의를 보면 평범한 방식으로 변수 입력 방식으로 처리되어 있다.
  // // 여기서는 빠르게 그림을 확인하려고 이렇게?
  // for (Int_t iSamp = 0; iSamp < dhPRMC->numEntries(); iSamp++)
  // {
  //   thisRow = (RooArgSet *)dhPRMC->get(iSamp);
  //   RooRealVar *myct = (RooRealVar *)thisRow->find("ctau3D");
  //   RooRealVar *mycterr = (RooRealVar *)thisRow->find("ctau3DErr");
  //   CtWeighted->setVal(myct->getVal() / mycterr->getVal());
  //   RooArgSet *tempRow = new RooArgSet(*CtWeighted);
  //   tempSet->add(*tempRow);
  // }

  // ws->factory("Gaussian::tmpGW_PRRes(CtWeighted,meanPRResW,sigmaPRResW)");
  // ws->factory("Gaussian::tmpGN_PRRes(CtWeighted,meanPRResN,sigmaPRResN)");
  // ws->factory("SUM::tmpCtPRRes(fracRes*tmpGW_PRRes,tmpGN_PRRes)");

  // RooPlot *tempFramePR = ws->var("CtWeighted")->frame();
  // tempSet->plotOn(tempFramePR, DataError(RooAbsData::SumW2));
  // ws->pdf("tmpCtPRRes")->plotOn(tempFramePR, LineColor(kGreen + 1), Normalization(tempSet->sumEntries(), RooAbsReal::NumEvent));

  // drawCtauResolPlots(ws, true, tempFramePR);


  // [5] ctau bkg fit
  // ----------------
  auto dhDataSB = new RooDataHist("dhDataSB", "", RooArgSet(*mass, *ctau3D, *ctau3DErr), *redDataSB);

  cout << " *** DATA :: N events to fit on SIDEBANDS : " << redDataSB->sumEntries() << endl;
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  RooFitResult *fitCt_Bkg = ws->pdf("CtBkgTot_PEE")->fitTo(*redDataSB, AsymptoticError(kTRUE), Minos(0), Save(1), NumCPU(32), PrintEvalErrors(-1), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), PrintLevel(0), EvalBackend("legacy")); //

  // RooFitResult *fitCt_Bkg = ws->pdf("CtBkgTot_PEE")->chi2FitTo(*dhDataSB, SumW2Error(kTRUE), Save(1),NumCPU(32), PrintEvalErrors(-1), PrintLevel(-1));
  // EvalBackend("legacy"),
  fitCt_Bkg->Print("v");
  ws->var("fracCtBkg1")->setConstant(kTRUE);
  ws->var("fracCtBkg2")->setConstant(kTRUE);
  ws->var("fracCtBkg3")->setConstant(kTRUE);
  ws->var("lambdap")->setConstant(kTRUE);
  ws->var("lambdam")->setConstant(kTRUE);
  ws->var("lambdasym")->setConstant(kTRUE);
  
  ws->var("sigmaPRResW")->setConstant(kTRUE); // test

  ws->var("fracRes")->setConstant(kTRUE); // test

  drawCtauSBPlots(ws, redDataSB, binDataCtErrSB, fitCt_Bkg, -ctMin, ctMax, &UnNormChi2_side, &nFitParam_side, &nFullBinsPull_side, &Dof_side, &Chi2_side);

  // // --- Perform final 2D fit ---
  // //// *** Get NSig, NBkg, Bfraction and their errors
  // Double_t NSigPR_fin, ErrNSigPR_fin;
  // Double_t NSigNP_fin, ErrNSigNP_fin;
  // Double_t Bfrac_fin, ErrBfrac_fin;
  // int nFitPar;
  // Double_t theNLL;
  // double resol, Errresol;
  // RooFitResult *fit2D;

  // cout << "\n-*-*-*-*-*-*-*-*-*-*-*- total fitting -*-*-*-*-*-*-*-*-*-*-*-  " << endl;
  // fit2D = ws->pdf("totPDF_PEE")->fitTo(*dhData, Minos(0), Save(1), SumW2Error(kTRUE), PrintEvalErrors(-1), NumCPU(32));
  // // fit2D = ws->pdf("totPDF_PEE")->fitTo(*dhData, Minos(0), Save(1), SumW2Error(kTRUE), PrintEvalErrors(-1), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), EvalBackend("legacy"));
  // // RecoverFromUndefinedRegions(1)
  // // NumCPU(32), EvalBackend("legacy"),
  // // 8 cpu, No EvalBAckend - 389 s
  // // 16 cpu, No EvalBAckend - 374s
  // // 8 cpu, EvalBAckend(legacy) - 403s
  // // 8 cpu, EvalBAckend(cpu) 378 ss
  // fit2D->Print("v");

  // // --- draw mass ---
  // TCanvas *c_mass = new TCanvas("c_mass", "2D fit projections", 600, 600);
  // RooPlot *fMass = mass->frame(Title("Projection on mass"));

  // dhData->plotOn(fMass, DataError(RooAbsData::SumW2)); // Binning(120),
  // ws->pdf("totPDF_PEE")->plotOn(fMass, ProjWData(RooArgSet(*ctau3D, *ctau3DErr), *dhData, true), LineColor(kBlue), Normalization(dhData->sumEntries(), RooAbsReal::NumEvent));
  // fMass->Draw();
  // c_mass->SaveAs("final_mass.png");

  // // --- draw ctau ---
  // TCanvas *c_ctau3D = new TCanvas("c_ctau3D", "2D fit projections (ctau3D)", 600, 600);
  // RooPlot *f_ctau3D = ctau3D->frame(Title("Projection on ctau3D"));
  // f_ctau3D->SetName("frame_ctau3D");

  // dhData->plotOn(f_ctau3D, DataError(RooAbsData::SumW2), Name("h_data_ctau3D"));
  // ws->pdf("totPDF_PEE")->plotOn(f_ctau3D, LineColor(kBlue), ProjWData(RooArgSet(*mass, *ctau3DErr), *dhData, true), Normalization(dhData->sumEntries(), RooAbsReal::NumEvent), Name("model_ctau3D"));
  // f_ctau3D->Draw();
  // c_ctau3D->SaveAs("final_ctau3D.png");

  // // Pull(ctau3D)
  // RooHist *hpull_ctau3D = f_ctau3D->pullHist("h_data_ctau3D", "model_ctau3D");
  // hpull_ctau3D->SetName("hpull_ctau3D");
  // TCanvas *c_pull_ctau3D = new TCanvas("c_pull_ctau3D", "Pull(ctau3D)", 600, 300);
  // RooPlot *fPull_ctau3D = ctau3D->frame(Title("Pull ctau3D"));
  // fPull_ctau3D->addPlotable(hpull_ctau3D, "P");
  
  // int npar = (fitCt_Bkg ? fit2D->floatParsFinal().getSize() : 0);
  // double chi2 = 0;
  // int n = 0;
  // for (int i = 0; i < hpull_ctau3D->GetN(); ++i)
  // {
  //   double x, y;
  //   hpull_ctau3D->GetPoint(i, x, y);
  //   if (TMath::Finite(y))
  //   {
  //     chi2 += y * y;
  //     ++n;
  //   }
  // }
  // int ndf = std::max(1, n - npar);
  // double chi2ndf = chi2 / ndf, p = TMath::Prob(chi2, ndf);

  // c_pull_ctau3D->cd();
  // fPull_ctau3D->SetMinimum(-5);
  // fPull_ctau3D->SetMaximum(5);
  // fPull_ctau3D->Draw();
  // fPull_ctau3D->addPlotable(hpull_ctau3D, "P");
  // (new TLine(-0.5, 0, 2.0, 0))->Draw("same"); // 0-라인
  // TLatex l;
  // l.SetNDC();
  // l.SetTextSize(0.045);
  // l.DrawLatex(0.60, 0.88, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));
  // l.DrawLatex(0.60, 0.82, Form("p = %.3g", p));
  // c_pull_ctau3D->SaveAs("final_ctau3D_pull.png");

  // ///////////////////////////////// ****  1) Final mass plots **** /////////////////////////////////
  // drawFinalMass(ws, dhData, NSigNP_fin, NBkg_fin, fitMass, &UnNormChi2_mass, &nFitParam_mass, &nFullBinsPull_mass, &Dof_mass, &Chi2_mass);
  // ///////////////////////////////// ****  2) Final ctau plots **** /////////////////////////////////
  // drawFinalCtau(ws, dhData, binDataCtErr, NSigNP_fin, NBkg_fin, Bfrac_fin, ErrBfrac_fin, fit2D, -ctMin, ctMax, &UnNormChi2_time, &nFitParam_time, &nFullBinsPull_time, &Dof_time, &Chi2_time);

  cout << "=== finish fit2D_622() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}

void drawCtauSBPlots(RooWorkspace *ws, RooDataSet *redDataSB, RooDataHist *binDataCtErrSB, RooFitResult *fitCt_Bkg, float lmin, float lmax, double *UnNormChi2_side_t, int *nFitParam_side_t, int *nFullBinsPull_side_t, int *Dof_side_t, double *Chi2_side_t)
{

  char reduceDS[1000];
  string titlestr;
  double unNormChi2;
  int dof;

  RooBinning rb(ws->var("ctau3D")->getBinning().numBins(), ws->var("ctau3D")->getBinning().array());

  TLegend leg11(0.64, 0.65, 0.9, 0.74, NULL, "brNDC");
  leg11.SetFillStyle(0);
  leg11.SetBorderSize(0);
  leg11.SetShadowColor(0);
  leg11.SetMargin(0.2);
  leg11.SetTextSize(0.040);
  leg11.SetTextFont(42);
  // leg11.AddEntry(gfake1, "sideband data", "p");
  // leg11.AddEntry(&hfake11, "background", "l");

  RooPlot *tframe1 = ws->var("ctau3D")->frame(Bins(100));
  double avgBinWidth = rb.averageBinWidth();
  tframe1->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)", avgBinWidth * 1000));
  tframe1->GetYaxis()->CenterTitle(1);
  redDataSB->plotOn(tframe1, DataError(RooAbsData::SumW2));

  // original
  ws->pdf("CtBkgTot_PEE")->plotOn(tframe1, ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErrSB, kTRUE), NumCPU(16), Normalization(binDataCtErrSB->sumEntries(), RooAbsReal::NumEvent), LineStyle(7));

  // ----------------------------------------------
  // // pjgwak
  // ws->pdf("CtBkgPos")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kBlue), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * ws->var("fracCtBkg1")->getVal(), RooAbsReal::NumEvent));
  // RooCurve *curvePos = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // // --- CtBkgNeg ---
  // ws->pdf("CtBkgNeg")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kGreen), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * (1 - ws->var("fracCtBkg1")->getVal()), RooAbsReal::NumEvent));
  // RooCurve *curveNeg = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // // --- CtBkgDbl ---
  // ws->pdf("CtBkgDbl")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kOrange), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * (1 - ws->var("fracCtBkg3")->getVal()) * (1 - ws->var("fracCtBkg2")->getVal()), RooAbsReal::NumEvent));
  // RooCurve *curveDbl = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // // --- CtPRRes (Resolution) ---
  // ws->pdf("CtPRRes")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kMagenta), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * ws->var("fracCtBkg3")->getVal(), RooAbsReal::NumEvent));
  // RooCurve *curveRes = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // auto sum12 = new RooCurve("sum12", "", *curvePos, *curveNeg);
  // auto sum123 = new RooCurve("sum123", "", *sum12, *curveDbl);
  // auto sumAll = new RooCurve("sumAll", "", *sum123, *curveRes);

  // sumAll->SetLineColor(kBlack);
  // sumAll->SetLineStyle(kSolid);
  // sumAll->SetLineWidth(3);
  // sumAll->SetMarkerStyle(1);
  // sumAll->SetDrawOption("L");

  // tframe1->addObject(sumAll);
  // ----------------------------------------------

  TCanvas *c3 = new TCanvas("c3", "The Canvas", 200, 10, 600, 750);
  c3->cd();
  TPad *pad1 = new TPad("pad1", "This is pad1", 0.00, 0.3, 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.075);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2", "This is pad2", 0.00, 0.00, 1.0, 0.3);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.30);
  pad2->Draw();

  pad1->cd();
  tframe1->Draw();
  // lumiTextOffset = 0.45;
  // CMS_lumi(c3, opt.isPA, iPos);
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetTextSize(0.040);
  // t->DrawLatex(0.92, 0.89, opt.rapString);
  // t->DrawLatex(0.92, 0.83, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.92, 0.78, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.92, 0.78, opt.etString);
  leg11.Draw("same");

  //// *** pull
  TH1 *hdatasb = redDataSB->createHistogram("hdatasb", *ws->var("ctau3D"), Binning(rb));
  RooHist *hpullsb = tframe1->pullHist();
  hpullsb->SetName("hpullSB");
  int nFitPar = fitCt_Bkg->floatParsFinal().getSize();
  double chi2 = 0;
  double *ypullssb = hpullsb->GetY();
  unsigned int nBins = ws->var("ctau3D")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++)
  {
    if (hdatasb->GetBinContent(i + 1) == 0)
      continue;
    chi2 += ypullssb[i] * ypullssb[i];
    nFullBins++;
  }
  unNormChi2 = chi2;
  *UnNormChi2_side_t = chi2;
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  int nDOF = ws->var("ctau3D")->getBinning().numBins() - nFitPar;
  *nFitParam_side_t = nFitPar;
  *nFullBinsPull_side_t = nFullBins;
  *Dof_side_t = dof;
  *Chi2_side_t = chi2;

  RooPlot *tframepull = ws->var("ctau3D")->frame(Title("Pull"));
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->GetYaxis()->CenterTitle(1);
  tframepull->SetLabelSize(0.04 * 2.5, "XYZ");
  tframepull->SetTitleSize(0.048 * 2.5, "XYZ");
  tframepull->SetTitleOffset(0.47, "Y");
  tframepull->addPlotable(hpullsb, "PX");
  double tframemax = 0;
  if (tframepull->GetMinimum() * -1 > tframepull->GetMaximum())
    tframemax = tframepull->GetMinimum() * -1;
  else
    tframemax = tframepull->GetMaximum();
  tframepull->SetMaximum(tframemax);
  tframepull->SetMinimum(-1 * tframemax);
  tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframepull->GetXaxis()->CenterTitle(1);

  pad2->cd();
  tframepull->Draw();
  TLine *line1 = new TLine(-lmin, 0, lmax, 0.);
  line1->SetLineStyle(7);
  line1->Draw();

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  t2->SetTextSize(0.035 * 3);
  sprintf(reduceDS, "#chi^{2}/dof = %.2f/%d", unNormChi2, dof);
  t2->DrawLatex(0.76, 0.90, reduceDS);

  titlestr = "_CtSB_Lin.png";
  c3->SaveAs(titlestr.c_str());

  pad1->SetLogy(1);
  double originalmax = tframe1->GetMaximum();
  tframe1->SetMaximum(originalmax * 10);
  tframe1->SetMinimum(0.5);
  titlestr = "_CtSB_Log.pdf";
  c3->SaveAs(titlestr.c_str());

  delete pad1;
  delete pad2;
  delete c3;
}