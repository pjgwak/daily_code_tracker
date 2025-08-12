#include "Final2DFit.h"
#include <iostream>
#include "TStyle.h"  // gStyle
#include <algorithm> // min()
#include "TSystem.h" // gSystem (make folder)
#include "RooMsgService.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TMath.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TLine.h"
#include <RooChebychev.h>
#include <RooCBShape.h>
#include "../../headers/polarizationUtilities.h"

// #include <RooGaussian.h>
// #include <RooFormulaVar.h>
// #include <RooCBShape.h>
// #include <RooWorkspace.h>
// #include <RooChebychev.h>
// #include <RooPolynomial.h>
// #include "RooPlot.h"
// #include "TText.h"
// #include "TArrow.h"
// #include "TFile.h"
// #include "RooDataHist.h"
// #include "RooCategory.h"
// #include "RooSimultaneous.h"
// #include "RooStats/SPlot.h"
// #include "../headers/cutsAndBin.h"
// #include "../headers/CMS_lumi_v2mass.C"
// #include "../headers/tdrstyle.C"
// #include "../headers/rootFitHeaders.h"
// #include "../headers/commonUtility.h"
// #include "../headers/JpsiUtility.h"

using std::cout; using std::endl;
using namespace RooFit;

Final2DFit::Final2DFit(float ptLow, float ptHigh,
                       float yLow, float yHigh,
                       int cLow, int cHigh,
                       float cosLow, float cosHigh,
                       int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE)
    : ptLow(ptLow), ptHigh(ptHigh),
      yLow(yLow), yHigh(yHigh),
      cLow(cLow), cHigh(cHigh),
      cosLow(cosLow), cosHigh(cosHigh),
      PR(PR), PRw(PRw),
      fEffW(fEffW), fAccW(fAccW), isPtW(isPtW), isTnP(isTnP),
      DATE(DATE)
{
  gStyle->SetEndErrorSize(0); // remove x direction error bars (cosmetic)
}

Final2DFit::~Final2DFit()
{
  delete fInputData;
}

void Final2DFit::init()
{
  cout << "--- init() ---\n\n";
  setLabels();
  makeOutputFolder();
  turnOffRooFitMessage();
  openInput();
  processDataset();
  setVariableRanges();
  fixParameters();
  buildCtauNpModel();
  buildMassCtauModel();
}

void Final2DFit::run()
{
  cout << "--- run() ---\n\n";
  do2dFit();
  // rootlogon(isLogon);
  drawCtauPull();
  drawMass();
  drawCtauRatio();
  saveOutput();
}

void Final2DFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
}

void Final2DFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
}

void Final2DFit::turnOffRooFitMessage()
{
  cout << "--- turnOffRooFitMessage() ---\n\n";
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);

  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}

void Final2DFit::openInput()
{
  cout << "--- openInput() ---\n\n";
}

void Final2DFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
}

void Final2DFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
}

void Final2DFit::fixParameters()
{
  cout << "--- fixParameters() ---\n\n";
}

void Final2DFit::buildCtauNpModel()
{
  cout << "--- buildCtauNpModel() ---\n\n";
}

void Final2DFit::buildMassCtauModel()
{
  cout << "--- buildMassCtauModel() ---\n\n";
}

void Final2DFit::do2dFit()
{
  cout << "--- do2dFit() ---\n\n";
}

void Final2DFit::drawCtauPull()
{
  cout << "--- drawCtauPull() ---\n\n";
}

void Final2DFit::drawMass()
{
  cout << "--- drawMass() ---\n\n";
}

void Final2DFit::drawCtauRatio()
{
  cout << "--- drawCtauRatio() ---\n\n";
}

void Final2DFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
}

// =================================
// ===== legacy codes or notes =====
// =================================

