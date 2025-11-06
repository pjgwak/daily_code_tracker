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
#include "RooHistPdfConv.h"
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
TGraphErrors *gfake1;
TH1F hfake11, hfake21, hfake31, hfake311, hfake41;
TLatex *t = new TLatex();
t->SetNDC();
t->SetTextAlign(12);
t->SetTextFont(42);
Double_t fx[2], fy[2], fex[2], fey[2];
//// black points (data, etc)
gfake1 = new TGraphErrors(2, fx, fy, fex, fey);
gfake1->SetMarkerStyle(20);
gfake1->SetMarkerSize(1);
//// blue line (background)
hfake11 = TH1F("hfake11", "hfake1", 100, 200, 300);
// hfake11.SetLineColor(kBlue); hfake11.SetLineWidth(4); hfake11.SetLineStyle(7); hfake11.SetFillColor(kAzure-9); hfake11.SetFillStyle(1001);
hfake11.SetLineColor(kBlue - 2);
hfake11.SetLineWidth(4);
hfake11.SetLineStyle(7);
hfake11.SetFillColor(kBlue - 10);
hfake11.SetFillStyle(1001);
//// black line (total fit)
hfake21 = TH1F("hfake21", "hfake2", 100, 200, 300);
// hfake21.SetLineColor(kBlack); hfake21.SetLineWidth(4); hfake21.SetFillColor(kBlack); hfake21.SetFillStyle(3354);
hfake21.SetLineColor(kBlack);
hfake21.SetLineWidth(4);
hfake21.SetFillColor(kGray + 2);
hfake21.SetFillStyle(3354);
//// red line (nonprompt)
hfake31 = TH1F("hfake31", "hfake3", 100, 200, 300);
// hfake31.SetLineColor(kRed); hfake31.SetMarkerStyle(kCircle); hfake31.SetLineWidth(4); hfake31.SetMarkerColor(kRed); hfake31.SetLineStyle(9); hfake31.SetFillColor(kRed-7); hfake31.SetFillStyle(3444);
hfake31.SetLineColor(kPink - 6);
hfake31.SetMarkerStyle(kCircle);
hfake31.SetLineWidth(4);
hfake31.SetMarkerColor(kPink - 6);
hfake31.SetLineStyle(12);
hfake31.SetFillColor(kRed - 7);
hfake31.SetFillStyle(3345);
hfake311 = TH1F("hfake311", "hfake311", 100, 200, 300);
hfake311.SetLineColor(kPink - 6);
hfake311.SetMarkerStyle(kCircle);
hfake311.SetLineWidth(4);
hfake311.SetMarkerColor(kPink - 6);
hfake311.SetLineStyle(12);
hfake311.SetFillColor(kRed - 7);
hfake311.SetFillStyle(3444);
//// green line (prompt)
hfake41 = TH1F("hfake41", "hfake4", 100, 200, 300);
hfake41.SetLineColor(kGreen + 3);
hfake41.SetMarkerStyle(kCircle);
hfake41.SetLineWidth(4);
hfake41.SetMarkerColor(kGreen + 3);
hfake41.SetLineStyle(kDashDotted);
hfake41.SetFillColor(kGreen - 7);
hfake41.SetFillStyle(3444);

void drawCtauSBPlots(RooWorkspace *ws, RooDataSet *redDataSB, RooDataHist *binDataCtErrSB, RooFitResult *fitCt_Bkg, float lmin, float lmax, double *UnNormChi2_side_t, int *nFitParam_side_t, int *nFullBinsPull_side_t, int *Dof_side_t, double *Chi2_side_t)
{
  char reduceDS[512];
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
  leg11.AddEntry(gfake1, "sideband data", "p");
  leg11.AddEntry(&hfake11, "background", "l");

  RooPlot *tframe1 = ws->var("ctau3D")->frame(); // Bins(100)
  double avgBinWidth = rb.averageBinWidth();
  tframe1->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)", avgBinWidth * 1000));
  tframe1->GetYaxis()->CenterTitle(1);
  redDataSB->plotOn(tframe1, DataError(RooAbsData::SumW2));

  // original
  // ws->pdf("CtBkgTot_PEE")->plotOn(tframe1, ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErrSB, kTRUE), NumCPU(16), Normalization(1, RooAbsReal::NumEvent), LineStyle(7));

  // ----------------------------------------------
  redDataSB->plotOn(tframe1, DataError(RooAbsData::SumW2), Name("dataPoints"));
  RooHist *hdata = tframe1->getHist("dataPoints");
  if (hdata)
  {
    std::cout << "Npoints = " << hdata->GetN() << std::endl;
  }
  int nPoints = hdata->GetN();

  // pjgwak
  double scalePos = (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * ws->var("fracCtBkg1")->getVal();
  if (scalePos < 0.01)
    scalePos = 0.01;
  ws->pdf("CtBkgPos")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kBlue), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * scalePos, RooAbsReal::NumEvent), Name("CtBkgPos"), Precision(1e-4));
  auto curvePos = dynamic_cast<RooCurve *>(tframe1->getCurve("CtBkgPos"));
  // RooCurve *curvePos = (RooCurve *)tframe1->getObject(tframe1->numItems() - 1);

  // --- CtBkgNeg ---
  double scaleNeg = (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * (1 - ws->var("fracCtBkg1")->getVal());
  if (scaleNeg < 0.01)
    scaleNeg = 0.01;

  ws->pdf("CtBkgNeg")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kGreen), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * scaleNeg, RooAbsReal::NumEvent), Name("CtBkgNeg"), Precision(1e-4));
  RooCurve *curveNeg = dynamic_cast<RooCurve *>(tframe1->getCurve("CtBkgNeg"));

  // --- CtBkgDbl ---
  double scaleDbl = (1 - ws->var("fracCtBkg3")->getVal()) * (1 - ws->var("fracCtBkg2")->getVal());
  if (scaleDbl < 0.01)
    scaleDbl = 0.01;

  ws->pdf("CtBkgDbl")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kOrange), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * scaleDbl, RooAbsReal::NumEvent), Name("CtBkgDbl"), Precision(1e-4));
  RooCurve *curveDbl = dynamic_cast<RooCurve *>(tframe1->getCurve("CtBkgDbl"));

  // --- CtPRRes (Resolution) ---
  ws->pdf("CtPRRes")->plotOn(tframe1, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *binDataCtErrSB), LineColor(kMagenta), LineStyle(kDashed), Normalization(binDataCtErrSB->sumEntries() * ws->var("fracCtBkg3")->getVal(), RooAbsReal::NumEvent), Name("CtPRRes"), Precision(1e-4));
  RooCurve *curveRes = dynamic_cast<RooCurve *>(tframe1->getCurve("CtPRRes"));

  // std::cerr << "Npos=" << (curvePos ? curvePos->GetN() : -1)
  //           << " Nneg=" << (curveNeg ? curveNeg->GetN() : -1)
  //           << " Ndbl=" << (curveDbl ? curveDbl->GetN() : -1)
  //           << " Nres=" << (curveRes ? curveRes->GetN() : -1) << std::endl;

  auto sum12 = new RooCurve("sum12", "", *curvePos, *curveNeg);
  auto sum123 = new RooCurve("sum123", "", *sum12, *curveDbl);
  auto sumAll = new RooCurve("sumAll", "", *sum123, *curveRes);

  sumAll->SetLineColor(kBlack);
  sumAll->SetLineStyle(kSolid);
  sumAll->SetLineWidth(3);
  sumAll->SetMarkerStyle(1);
  sumAll->SetDrawOption("L");

  tframe1->addObject(sumAll);
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
  titlestr = "_CtSB_Lin.pdf";
  c3->Modified();
  c3->Update();
  c3->SaveAs(titlestr.c_str());

  pad1->SetLogy(1);
  double originalmax = tframe1->GetMaximum();
  tframe1->SetMaximum(originalmax * 10);
  tframe1->SetMinimum(0.5);
  titlestr = "_CtSB_Log.pdf";
  c3->Modified();
  c3->Update();
  c3->SaveAs(titlestr.c_str());

  delete pad1;
  delete pad2;
  delete c3;
}

// Helper functions
RooBinning setCtBinning(float lmin, float lmax)
{
  RooBinning rbct(-lmin, lmax);
  if (lmax + lmin > 4.9)
  {
    rbct.addBoundary(-1.5);
    rbct.addBoundary(-1.0);
    rbct.addBoundary(-0.8);
    rbct.addBoundary(-0.6);
    rbct.addBoundary(-0.5);
    rbct.addUniform(6, -0.5, -0.2);
    rbct.addUniform(12, -0.2, 0.1);
    rbct.addUniform(8, 0.1, 0.5);
    rbct.addUniform(5, 0.5, 1.0);
    rbct.addUniform(15, 1.0, lmax);
  }
  else if (lmax + lmin > 4.4)
  { // this is what we use KYO!!!
    rbct.addBoundary(-1.5);
    rbct.addBoundary(-1.0);
    rbct.addBoundary(-0.8);
    rbct.addBoundary(-0.6);
    rbct.addBoundary(-0.5);
    rbct.addUniform(9, -0.5, -0.2);
    rbct.addUniform(20, -0.2, 0.1);
    rbct.addUniform(11, 0.1, 0.5);
    rbct.addUniform(5, 0.5, 1.0);
    rbct.addUniform(5, 1.0, lmax);
  }
  else
  {
    rbct.addBoundary(-lmin);
    rbct.addBoundary(-0.7);
    rbct.addBoundary(-0.6);
    rbct.addBoundary(-0.5);
    rbct.addUniform(9, -0.5, -0.2);
    rbct.addUniform(28, -0.2, 0.1);
    rbct.addUniform(11, 0.1, 0.5);
    rbct.addUniform(15, 0.5, 1.2);
    rbct.addUniform(8, 1.2, lmax);
  }

  
  // rbct.addUniform(20, -0.2, 0.1);
  // rbct.addUniform(11, 0.1, 0.4);
  // rbct.addUniform(5, 0.5, 1.0);
  // rbct.addUniform(5, 1.0, lmax);
  return rbct;
}


/////////////////////////////////////////////////////////



void drawFinalMass(RooWorkspace *ws, RooDataSet *redDataCut, float NSigNP_fin, float NBkg_fin, RooFitResult *fitMass, double *UnNormChi2_mass_t, int *nFitParam_mass_t, int *nFullBinsPull_mass_t, int *Dof_mass_t, double *Chi2_mass_t)
{
  //// **** Temporary variables for plotting
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  TLatex *ty = new TLatex();
  ty->SetNDC();
  ty->SetTextAlign(12);
  char reduceDS[512];
  string titlestr;

  RooBinning rb(ws->var("mass")->getBinning().numBins(), ws->var("mass")->getBinning().array());
  // RooRealVar tmpVar1("tmpVar1", "tmpVar1", NSigNP_fin);
  // RooRealVar tmpVar2("tmpVar2", "tmpVar2", NBkg_fin);

  //// *** Mass plot
  RooPlot *mframe = ws->var("mass")->frame();

  redDataCut->plotOn(mframe, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1)); // Binning(rb)
  double avgBinWidth = rb.averageBinWidth();
  mframe->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/c^{2 })", avgBinWidth * 1000));
  mframe->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  //  mframe->GetXaxis()->SetTitleSize(0.048*1.2); //PAPER
  //  mframe->GetYaxis()->SetTitleSize(0.048*1.2); //PAPER
  mframe->GetXaxis()->CenterTitle(1);
  mframe->GetYaxis()->CenterTitle(1);
  const double max = mframe->GetMaximum() * 1.3;
  mframe->SetMaximum(max);
  mframe->SetMinimum(0);

  // ----------------
  //// **** Fill color
  // 전체 -> mass는 projWData 없이? -> 느리다.
  
  // ws->pdf("totPDF_PEE")->plotOn(mframe, DrawOption("F"), FillColor(kGray + 2), FillStyle(3354), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));

  // Signal PR (Just mass signal model) + bkg
  // Normalize to number of Jpsi in the mass window
  
  
  

  // MassCtBkg_PEE,Bfrac[0.25,0.0,1.]*MassCtNP_PEE,MassCtPR_PEE

  // double sigPrScale = ws->function("fracBkg")->getVal() + (1 - ws->function("fracBkg")->getVal()) * (1.0 - ws->var("Bfrac")->getVal());

  double fracBkg = ws->function("fracBkg")->getVal();
  double sigNpScale = (1 - ws->function("fracBkg")->getVal()) * ws->var("Bfrac")->getVal();

  RooAddPdf tmpPDF("tmpPDF", "tmpPDF", RooArgList(*(ws->pdf("G1CB1Sig")), *(ws->pdf("expBkg"))), RooArgList(sigNpScale, fracBkg));

  tmpPDF.plotOn(mframe, LineColor(kPink - 6), Normalization(redDataCut->sumEntries() * (sigNpScale+fracBkg), RooAbsReal::NumEvent));
  
  ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), LineColor(kGreen - 10), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));

  // ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kPink - 6), Normalization(redDataCut->sumEntries() * sigNpScale, RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kOrange - 2), Components(" expBkg"));
  

  // // Sig NP + bkg
  // // tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, Components("G1CB1Sig, expBkg"), LineColor(kPink - 6), DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), Normalization(redDataCut->sumEntries() * sigNpScale, RooAbsReal::NumEvent));
  // // gStyle->SetHatchesLineWidth(2);

  // bkg only
  // ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  
  // //// **** Line color
  // ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), LineColor(kBlue - 2), LineStyle(7), LineWidth(5), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));

  // tmpPDF.plotOn(mframe, LineColor(kPink - 6), LineStyle(11), LineWidth(5), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent));

  // ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // // ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  ws->pdf("totPDF_PEE")->plotOn(mframe, Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // ----------------


  // ------------ original code start -----------------
  // //// **** Fill color
  // // ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, DrawOption("F"), FillColor(kGray + 2), FillStyle(3354), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // RooAddPdf tmpPDF("tmpPDF", "tmpPDF", RooArgList(*(ws->pdf("G1CB1Sig")), *(ws->pdf("expBkg"))), RooArgList(tmpVar1, tmpVar2));
  // // tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
  // tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kWhite), FillStyle(1001), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent));
  // // tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
  // tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent));
  // gStyle->SetHatchesLineWidth(2);
  // // ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // //// **** Line color
  // // ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), LineColor(kBlue - 2), LineStyle(7), LineWidth(5), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // // tmpPDF.plotOn(mframe,LineColor(kRed),LineStyle(9),LineWidth(5),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
  // tmpPDF.plotOn(mframe, LineColor(kPink - 6), LineStyle(11), LineWidth(5), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // redDataCut->plotOn(mframe, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));
  // ------------ original code end -----------------

  TH1 *hdata = redDataCut->createHistogram("hdata", *ws->var("mass"), Binning(rb));//Binning(rb)
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
  TCanvas *c1wop = new TCanvas("c1wop", "The Canvas", 200, 10, 600, 600);
  c1wop->cd();
  c1wop->Draw();

  mframe->SetTitleOffset(1.47, "Y");
  mframe->Draw();
  // //// **** different lumiTextOffset for massfit_wopull
  // lumiTextOffset = 0.20;
  // CMS_lumi(c1wop, opt.isPA, iPosPaper);
  // t->SetTextSize(0.035);
  // //  t->DrawLatex(0.91,0.85,opt.rapString);
  // //  t->DrawLatex(0.91,0.78,opt.ptString);
  // ty->SetTextSize(0.035);                   // PAPER
  // ty->DrawLatex(0.20, 0.86, opt.rapString); // PAPER
  // ty->DrawLatex(0.20, 0.80, opt.ptString);  // PAPER
  // ty->SetTextSize(0.040);
  // sprintf(reduceDS, "N_{J/#psi} = %0.0f #pm %0.0f", ws->var("NSig")->getVal(), ws->var("NSig")->getError());
  // //  ty->DrawLatex(0.20,0.85,reduceDS);
  // sprintf(reduceDS, "#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
  // //  ty->DrawLatex(0.20,0.79,reduceDS);

  // TLegend * legpaper = new TLegend(0.59,0.53,0.90,0.72,NULL,"brNDC");
  // TLegend * legpaper = new TLegend(0.17,0.53,0.63,0.72,NULL,"brNDC");
  TLegend *legpaper = new TLegend(0.17, 0.55, 0.57, 0.75, NULL, "brNDC");
  legpaper->SetFillStyle(0);
  legpaper->SetBorderSize(0);
  legpaper->SetShadowColor(0);
  legpaper->SetTextSize(0.035);
  legpaper->SetTextFont(42);
  legpaper->SetMargin(0.2);
  legpaper->AddEntry(gfake1, "Data", "p");
  legpaper->AddEntry(&hfake21, "Total fit", "lf");
  legpaper->AddEntry(&hfake31, "Bkg + nonprompt", "lf");
  legpaper->AddEntry(&hfake11, "Background", "lf");
  legpaper->Draw("same");

  titlestr = "_massfit_wopull.pdf";
  c1wop->SaveAs(titlestr.c_str());
  //  titlestr = opt.dirName + "_rap" + opt.yrange + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfit_wopull.root";
  //  c1wop.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////

  TCanvas c1("c1", "The mass Canvas", 200, 10, 600, 750);
  c1.cd();
  TPad *padm1 = new TPad("padm1", "This is pad1", 0.0, 0.3, 1.0, 1.0);
  padm1->SetLeftMargin(0.14);
  padm1->SetRightMargin(0.03);
  padm1->SetTopMargin(0.075);
  padm1->SetBottomMargin(0);
  padm1->Draw();
  TPad *padm2 = new TPad("padm2", "This is pad2", 0.00, 0.00, 1.0, 0.3);
  padm2->SetLeftMargin(0.14);
  padm2->SetRightMargin(0.03);
  padm2->SetTopMargin(0);
  padm2->SetBottomMargin(0.30);
  padm2->Draw();

  padm1->cd();
  mframe->Draw();
  // lumiTextOffset = 0.45;
  // CMS_lumi(&c1, opt.isPA, iPos);
  // t->SetTextSize(0.035);
  // t->DrawLatex(0.91, 0.90, opt.rapString);
  // t->DrawLatex(0.91, 0.85, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.91, 0.80, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.91, 0.80, opt.etString);

  ty->SetTextSize(0.040);
  sprintf(reduceDS, "N_{J/#psi} = %0.0f #pm %0.0f", ws->var("NSig")->getVal(), ws->var("NSig")->getError());
  ty->DrawLatex(0.20, 0.89, reduceDS);
  // sprintf(reduceDS, "#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
  // ty->DrawLatex(0.20, 0.84, reduceDS);

  TLegend *leg11 = new TLegend(0.18, 0.51, 0.54, 0.72, NULL, "brNDC");
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0);
  leg11->SetShadowColor(0);
  leg11->SetTextSize(0.035);
  leg11->SetTextFont(42);
  leg11->SetMargin(0.2);
  leg11->AddEntry(gfake1, "Data", "p");
  leg11->AddEntry(&hfake21, "Total fit", "lf");
  leg11->AddEntry(&hfake31, "Bkg + nonprompt", "lf");
  leg11->AddEntry(&hfake11, "Background", "lf");
  leg11->Draw("same");
  c1.Update();

  //// **** pull
  RooPlot *mframepull = ws->var("mass")->frame(Title("Pull"));
  mframepull->GetYaxis()->SetTitle("Pull");
  mframepull->GetYaxis()->CenterTitle(1);
  mframepull->SetLabelSize(0.04 * 2.5, "XYZ");
  mframepull->SetTitleSize(0.048 * 2.5, "XYZ");
  mframepull->SetTitleOffset(0.47, "Y");
  mframepull->addPlotable(hpullm, "PX");
  double mframemax = 0;
  if (mframepull->GetMinimum() * -1 > mframepull->GetMaximum())
    mframemax = mframepull->GetMinimum() * -1;
  else
    mframemax = mframepull->GetMaximum();
  mframepull->SetMaximum(mframemax);
  mframepull->SetMinimum(-1 * mframemax);
  mframepull->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  mframepull->GetXaxis()->CenterTitle(1);

  padm2->cd();
  mframepull->Draw();
  TLine *line1 = new TLine(2.6, 0, 3.5, 0.);
  line1->SetLineStyle(7);
  line1->Draw();

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  t2->SetTextSize(0.035 * 3);
  sprintf(reduceDS, "#chi^{2}/dof = %.1f/%d", UnNormChi2, Dof);
  t2->DrawLatex(0.78, 0.86, reduceDS);
  c1.Update();

  titlestr = "_massfit.pdf";
  c1.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////
  mframe->SetMinimum(0.5);
  mframe->SetMaximum(max * 50);
  padm1->cd();
  padm1->SetLogy(1);
  titlestr = "_massfit_Log.pdf";
  c1.SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

void drawFinalCtau(RooWorkspace *ws, RooDataSet *redDataCut, RooDataHist *binDataCtErr, float NSigNP_fin, float NBkg_fin, float Bfrac_fin, float ErrBfrac_fin, RooFitResult *fit2D, float lmin, float lmax, double *UnNormChi2_time_t, int *nFitParam_time_t, int *nFullBinsPull_time_t, int *Dof_time_t, double *Chi2_time_t)
{
  // pjgwak - here!
  char reduceDS[512];
  string titlestr;

  RooBinning rb(ws->var("ctau3D")->getBinning().numBins(), ws->var("ctau3D")->getBinning().array());

  RooRealVar tmpVar1("tmpVar1", "tmpVar1", NSigNP_fin);
  RooRealVar tmpVar2("tmpVar2", "tmpVar2", NBkg_fin);

  RooPlot *tframe = ws->var("ctau3D")->frame();
  tframe->SetTitleOffset(1.47, "Y");
  double avgBinWidth = rb.averageBinWidth();
  // tframe->GetYaxis()->SetTitle(Form("Counts / %.2f (mm)",avgBinWidth));
  tframe->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)", avgBinWidth * 1000));
  tframe->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  tframe->GetXaxis()->CenterTitle(1);
  tframe->GetYaxis()->CenterTitle(1);

  //// **** Ctau total distributions
  RooHist *hpulltot;
  redDataCut->plotOn(tframe, DataError(RooAbsData::SumW2), Binning(rb), MarkerSize(1), Name("data"));



  // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtBkg"), LineColor(kBlue - 2), LineWidth(5), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(7));

  // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtNP"), LineColor(kPink - 6), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(12));

  // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtPR"), LineColor(kGreen + 3), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(kDashDotted));

  // CtPRRes *(1 - Bfrac - fracBkg)
  //   ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 5],one[1.0],ctau3DErr)");
  // ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.01, 2],one,ctau3DErr)");
  // ws->factory("AddModel::CtPRRes({GW_PRRes,GN_PRRes},{fracRes[0.5,0.01,0.999]})");

  double scalePR = (1 - ws->function("fracBkg")->getVal())*(1.0 - ws->var("Bfrac")->getVal());
  double fracRes = ws->var("fracRes")->getVal();

  // --- GW_PRRes (wide) ---
  ws->pdf("GW_PRRes")->plotOn(tframe, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut), LineColor(kBlue), LineStyle(kDashed), Normalization(redDataCut->sumEntries() * scalePR * fracRes, RooAbsReal::NumEvent));
  RooCurve *curveGW = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  // --- GN_PRRes (narrow) ---
  ws->pdf("GN_PRRes")->plotOn(tframe, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut), LineColor(kRed), LineStyle(kDashed), Normalization(redDataCut->sumEntries() * scalePR * (1 - fracRes), RooAbsReal::NumEvent));
  RooCurve *curveGN = (RooCurve *)tframe->getObject(tframe->numItems() - 1);


  // CtNPTot * Bfrac
  // ws->factory("Decay::CtNPTot(ctau3D,coefExpNPTrue,CtNPRes,RooDecay::SingleSided)");
  ws->pdf("CtNPTot")->plotOn(tframe,
                             ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut),
                             LineColor(kRed), LineStyle(kDashed),
                             Normalization(
                                 redDataCut->sumEntries() * (1 - ws->function("fracBkg")->getVal()) * ws->var("Bfrac")->getVal(),
                                 RooAbsReal::NumEvent));
  RooCurve *curveNP = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  // ----------------------------------------------
  // pjgwak

  // --------  CtBkgTot * fracBkg ------------
  double fracBkg = (ws->function("fracBkg")->getVal());

  // --- CtBkgPos ---
  ws->pdf("CtBkgPos")->plotOn(tframe, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut), LineColor(kBlue), LineStyle(kDashed), Normalization(redDataCut->sumEntries() * fracBkg * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * ws->var("fracCtBkg1")->getVal(), RooAbsReal::NumEvent));
  RooCurve *curvePos = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  // --- CtBkgNeg ---
  ws->pdf("CtBkgNeg")->plotOn(tframe, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut), LineColor(kGreen), LineStyle(kDashed), Normalization(redDataCut->sumEntries() * fracBkg * (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * (1 - ws->var("fracCtBkg1")->getVal()), RooAbsReal::NumEvent));
  RooCurve *curveNeg = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  // --- CtBkgDbl ---
  ws->pdf("CtBkgDbl")->plotOn(tframe, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut), LineColor(kOrange), LineStyle(kDashed), Normalization(redDataCut->sumEntries() * fracBkg * (1 - ws->var("fracCtBkg3")->getVal()) * (1 - ws->var("fracCtBkg2")->getVal()), RooAbsReal::NumEvent));
  RooCurve *curveDbl = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  // --- CtPRRes (Resolution) ---
  ws->pdf("CtPRRes")->plotOn(tframe,
                             ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataCut),
                             LineColor(kMagenta), LineStyle(kDashed),
                             Normalization(redDataCut->sumEntries() *
                                               fracBkg *
                                               ws->var("fracCtBkg3")->getVal(),
                                           RooAbsReal::NumEvent));
  RooCurve *curveRes = (RooCurve *)tframe->getObject(tframe->numItems() - 1);

  auto curveResTot = new RooCurve("curveResTot", "GW+GN", *curveGW, *curveGN);

  // --- Background (Pos + Neg + Dbl + Res)
  auto curveBkgPN = new RooCurve("curveBkgPN", "Pos+Neg", *curvePos, *curveNeg);
  auto curveBkgPND = new RooCurve("curveBkgPND", "Pos+Neg+Dbl", *curveBkgPN, *curveDbl);
  auto curveBkgTot = new RooCurve("curveBkgTot", "Pos+Neg+Dbl+Res", *curveBkgPND, *curveRes);

  // --- final: NonPrompt + Background + PR resolution (GW+GN)
  auto curveTmp = new RooCurve("curveTmp", "", *curveResTot, *curveNP);
  auto curveAll = new RooCurve("curveAll", "Total Model", *curveTmp, *curveBkgTot);

  curveAll->SetLineColor(kBlack);
  curveAll->SetLineStyle(kSolid);
  curveAll->SetLineWidth(3);
  curveAll->SetMarkerStyle(1);
  curveAll->SetDrawOption("L");
  tframe->addObject(curveAll);

  hpulltot = tframe->pullHist("data", "curveAll");
  hpulltot->SetName("hpulltot");
  //   TCanvas *c1 = new TCanvas("c1", "Fit with pull", 800, 800);
  //   c1->Divide(1, 2);
  //   c1->cd(1);
  //   gPad->SetPad(0, 0.3, 1, 1);
  //   tframe->Draw();

  // // Pull plot
  //   // ------------------------------
  //   c1->cd(2);
  //   gPad->SetPad(0, 0, 1, 0.3);
  //   RooPlot *pullFrame = tframe->emptyClone("pullFrame");

  //   // 데이터와 모델 비교 → pull
  //   RooHist *hpull = tframe->pullHist("data", "curveAll");

  //   pullFrame->addPlotable(hpull, "P");
  //   pullFrame->SetTitle("");
  //   pullFrame->GetYaxis()->SetTitle("Pull");
  //   pullFrame->GetYaxis()->SetNdivisions(505);
  //   pullFrame->GetYaxis()->SetLabelSize(0.1);
  //   pullFrame->GetXaxis()->SetLabelSize(0.12);

  //   pullFrame->Draw();

  //   double chi2 = tframe->chiSquare("curveAll", "data");
  //   std::cout << "Chi2/ndf = " << chi2 << std::endl;

  //   c1->SaveAs("fit_with_pull.png");

  // ----------------------------------------------

  // // ---------- original code ---------------
  // // ws->pdf("totPDF_PEE")->plotOn(tframe, LineColor(kBlack), LineWidth(2), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent));
  // hpulltot = tframe->pullHist();
  // hpulltot->SetName("hpulltot");
  // // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtBkg"), LineColor(kBlue - 2), LineWidth(5), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(7));
  // // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtNP"), LineColor(kPink - 6), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(12));
  // // ws->pdf("totPDF_PEE")->plotOn(tframe, Components("MassCtPR"), LineColor(kGreen + 3), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(kDashDotted));
  // // ws->pdf("totPDF_PEE")->plotOn(tframe, LineColor(kBlack), LineWidth(2), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent));
  // // ---------- original code end ---------------

  TH1 *hdatact = redDataCut->createHistogram("hdatact", *ws->var("ctau3D"), Binning(rb));
  double chi2 = 0, unNormChi2 = 0;
  int dof = 0;
  double *ypulls = hpulltot->GetY();
  unsigned int nBins = ws->var("ctau3D")->getBinning().numBins();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nBins; i++)
  {
    if (hdatact->GetBinContent(i + 1) == 0)
      continue;
    chi2 += ypulls[i] * ypulls[i];
    nFullBins++;
  }
  unNormChi2 = chi2;
  *UnNormChi2_time_t = chi2;
  int nFitPar = fit2D->floatParsFinal().getSize();
  dof = nFullBins - nFitPar;
  chi2 /= (nFullBins - nFitPar);
  *nFitParam_time_t = nFitPar;
  *nFullBinsPull_time_t = nFullBins;
  *Dof_time_t = dof;
  *Chi2_time_t = chi2;

  TCanvas *c2 = new TCanvas("c2", "The Canvas", 200, 10, 600, 750);
  c2->cd();
  TPad *pad1 = new TPad("pad1", "This is pad1", 0.0, 0.3, 1.0, 1.0);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.075);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2", "This is pad2", 0.0, 0.0, 1.0, 0.3);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.30);
  pad2->Draw();

  pad1->cd();
  tframe->Draw();

  // lumiTextOffset = 0.45;
  // CMS_lumi(c2, opt.isPA, iPos);
  // TLatex *ty = new TLatex();
  // ty->SetNDC();
  // ty->SetTextAlign(32);
  // ty->SetTextSize(0.035);
  // ty->DrawLatex(0.91, 0.90, opt.rapString);
  // ty->DrawLatex(0.91, 0.85, opt.ptString);
  // if (opt.EventActivity == 1)
  //   ty->DrawLatex(0.91, 0.80, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   ty->DrawLatex(0.91, 0.80, opt.etString);

  // TLegend * leg = new TLegend(0.66,0.56,0.85,0.75,NULL,"brNDC");
  TLegend *leg = new TLegend(0.66, 0.56, 0.99, 0.75, NULL, "brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  leg->SetMargin(0.2);
  leg->AddEntry(gfake1, "Data", "p");
  leg->AddEntry(&hfake21, "Total fit", "l");
  leg->AddEntry(&hfake41, "Prompt", "l");
  leg->AddEntry(&hfake311, "Nonprompt", "l");
  leg->AddEntry(&hfake11, "Background", "l");
  leg->Draw("same");

  // KYO : write down Bfrac on plot
  TLatex *tbfrac = new TLatex();
  tbfrac->SetNDC();
  tbfrac->SetTextAlign(12);
  tbfrac->SetTextSize(0.035);
  sprintf(reduceDS, "B frac. = %.2f #pm %.2f", Bfrac_fin, ErrBfrac_fin);
  tbfrac->DrawLatex(0.19, 0.91, reduceDS);

  RooPlot *tframepull = ws->var("ctau3D")->frame(Title("Pull"));
  tframepull->GetYaxis()->SetTitle("Pull");
  tframepull->GetYaxis()->CenterTitle(1);
  tframepull->SetLabelSize(0.04 * 2.5, "XYZ");
  tframepull->SetTitleSize(0.048 * 2.5, "XYZ");
  tframepull->SetTitleOffset(0.47, "Y");
  tframepull->addPlotable(hpulltot, "PX");
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

  int nDOF = ws->var("ctau3D")->getBinning().numBins() - nFitPar;

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextAlign(22);
  t2->SetTextSize(0.035 * 3);
  sprintf(reduceDS, "#chi^{2}/dof = %.2f/%d", unNormChi2, dof);
  t2->DrawLatex(0.78, 0.90, reduceDS);

  c2->Update();
  titlestr = "_timefit_Lin.pdf";
  c2->SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////
  tframe->SetMaximum(tframe->GetMaximum() * 9);
  tframe->SetMinimum(0.5);
  pad1->SetLogy(1);
  titlestr = "_timefit_Log.pdf";
  c2->SaveAs(titlestr.c_str());
  /////////////////////////////////////////////////////////////////////////

  // TCanvas *c2b = new TCanvas("c2b", "The Canvas", 200, 10, 600, 600);
  // c2b->cd();
  // c2b->Draw();
  // c2b->SetLogy(1);

  // RooPlot *tframefill = ws->var("ctau3D")->frame();
  // tframefill->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)", avgBinWidth * 1000));
  // redDataCut->plotOn(tframefill, DataError(RooAbsData::SumW2), Binning(rb), MarkerSize(1));
  // // ws->pdf("totPDF_PEE")->plotOn(tframefill,DrawOption("F"),FillColor(kBlack),FillStyle(3354),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(8),Normalization(1,RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(tframefill, DrawOption("F"), FillColor(kGray + 2), FillStyle(3354), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(tframefill, LineColor(kBlack), LineWidth(2), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent));
  // RooAddPdf tmpPDF2("tmpPDF2", "tmpPDF2", RooArgList(*(ws->pdf("MassCtNP")), *(ws->pdf("MassCtBkg"))), RooArgList(tmpVar1, tmpVar2));
  // tmpPDF2.plotOn(tframefill, DrawOption("F"), FillColor(kWhite), FillStyle(1001), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization((NSigNP_fin + NBkg_fin) / redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // // tmpPDF2.plotOn(tframefill,DrawOption("F"),FillColor(kRed),FillStyle(3444),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(8),Normalization((NSigNP_fin+NBkg_fin)/redDataCut->sumEntries(),RooAbsReal::NumEvent));
  // tmpPDF2.plotOn(tframefill, DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization((NSigNP_fin + NBkg_fin) / redDataCut->sumEntries(), RooAbsReal::NumEvent));
  // gStyle->SetHatchesLineWidth(2);
  // // ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("MassCtBkg"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(8),Normalization(1,RooAbsReal::NumEvent));
  // ws->pdf("totPDF_PEE")->plotOn(tframefill, Components("MassCtBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent));
  // // ws->pdf("totPDF_PEE")->plotOn(tframefill,Components("MassCtBkg"),LineColor(kBlue),LineWidth(5),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(8),Normalization(1,RooAbsReal::NumEvent),LineStyle(7));
  // ws->pdf("totPDF_PEE")->plotOn(tframefill, Components("MassCtBkg"), LineColor(kBlue - 2), LineWidth(5), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization(1, RooAbsReal::NumEvent), LineStyle(7));
  // // tmpPDF2.plotOn(tframefill,LineColor(kRed),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(8),Normalization((NSigNP_fin+NBkg_fin)/redDataCut->sumEntries(),RooAbsReal::NumEvent),LineWidth(5),LineStyle(9));
  // tmpPDF2.plotOn(tframefill, LineColor(kPink - 6), ProjWData(RooArgList(*(ws->var("ctau3DErr"))), *binDataCtErr, kTRUE), NumCPU(8), Normalization((NSigNP_fin + NBkg_fin) / redDataCut->sumEntries(), RooAbsReal::NumEvent), LineWidth(5), LineStyle(11));

  // redDataCut->plotOn(tframefill, DataError(RooAbsData::SumW2), Binning(rb), MarkerSize(1));

  // tframefill->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  // tframefill->GetXaxis()->CenterTitle(1);
  // //  tframefill->GetXaxis()->SetTitleSize(0.048*1.2); //PAPER
  // //  tframefill->GetYaxis()->SetTitleSize(0.048*1.2); //PAPER
  // //  tframefill->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)",avgBinWidth*1000));
  // tframefill->GetYaxis()->SetTitle(Form("Counts / <%.0f #mum>", avgBinWidth * 1000)); // PAPER
  // tframefill->GetYaxis()->CenterTitle(1);
  // tframefill->SetMaximum(tframefill->GetMaximum() * 9);
  // tframefill->SetMinimum(0.5);
  // tframefill->Draw();

  // // TLegend * legpaper = new TLegend(0.59,0.53,0.90,0.72,NULL,"brNDC");
  // // TLegend * legpaper = new TLegend(0.54,0.54,0.99,0.74,NULL,"brNDC");
  // TLegend *legpaper = new TLegend(0.59, 0.55, 0.99, 0.75, NULL, "brNDC");
  // legpaper->SetFillStyle(0);
  // legpaper->SetBorderSize(0);
  // legpaper->SetShadowColor(0);
  // legpaper->SetTextSize(0.035);
  // legpaper->SetTextFont(42);
  // legpaper->SetMargin(0.2);
  // legpaper->AddEntry(gfake1, "Data", "p");
  // legpaper->AddEntry(&hfake21, "Total fit", "lf");
  // legpaper->AddEntry(&hfake31, "Bkg + nonprompt", "lf");
  // legpaper->AddEntry(&hfake11, "Background", "lf");
  // legpaper->Draw("same");

  // // //// **** different lumiTextOffset for timefit_wopull
  // // lumiTextOffset = 0.20;
  // // CMS_lumi(c2b, opt.isPA, iPosPaper);
  // // TLatex *t = new TLatex();
  // // t->SetNDC();
  // // t->SetTextAlign(32);
  // // t->SetTextSize(0.035);
  // // //  t->DrawLatex(0.91,0.85,opt.rapString);
  // // //  t->DrawLatex(0.91,0.78,opt.ptString);
  // // t->SetTextAlign(12);                     // PAPER
  // // t->DrawLatex(0.20, 0.86, opt.rapString); // PAPER
  // // t->DrawLatex(0.20, 0.80, opt.ptString);  // PAPER

  // c2b->Update();
  // titlestr = "_timefit_Log_wopull.pdf";
  // c2b->SaveAs(titlestr.c_str());
  // //  titlestr = opt.dirName + "_rap" + opt.yrange + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_timefit_Log_wopull.root";
  // //  c2b->SaveAs(titlestr.c_str());
  // delete c2b;
}

void fit2D()
{
  cout << "=== start fit2D() ===\n";
  TStopwatch t;
  t.Start();

  cout << "=== finish fit2D() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}