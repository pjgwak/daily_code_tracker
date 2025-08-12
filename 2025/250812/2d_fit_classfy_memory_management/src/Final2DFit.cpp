#include "Final2DFit.h"
#include <iostream>
#include <algorithm> // min()
#include "TStyle.h" // gStyle
#include "RooMsgService.h"
#include "TSystem.h" // gSystem (make folder)
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TMath.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TLine.h"
#include "/work/pjgwak/pol24/headers/polarizationUtilities.h"

using std::cout; using std::endl;
using std::min;
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

Final2DFit::~Final2DFit() {
  delete fInputData;
  delete fMass;
  delete fCErr;
  delete fCRes;
  delete fCBkg;
  delete fCTrue;
  delete fMcParams;
  delete fitResult;
  delete ws;
}

void Final2DFit::init() {
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

void Final2DFit::run() {
  cout << "--- run() ---\n\n";
  do2dFit();
  rootlogon(isLogon);
  drawCtauPull();
  drawMass();
  drawCtauRatio();
  saveOutput();
}

void Final2DFit::setLabels() {
  cout << "--- setLabels() ---\n\n";
  if (PRw == 1) fname = "PR";
  else if (PRw == 2) fname = "NP";

  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);
}

void Final2DFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/Final", DATE.Data()));
  gSystem->mkdir(Form("figs/2DFit_%s/Final", DATE.Data()));
}

void Final2DFit::turnOffRooFitMessage() {
  cout << "--- turnOffRooFitMessage() ---\n\n";
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling); // ignore successful RooWorkspace importing message

  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(1).removeTopic(Tracing);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}

void Final2DFit::openInput() {
  cout << "--- openInput() ---\n\n";
  fInputData = new TFile(Form("/disk1/Oniatree/miniAOD/Run2025OO/OniaRooDataSet_miniAOD_2025OORun_isMC0_Charmonia_Effw0_Accw0_PtW0_TnP0_250722.root"));
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCBkg = new TFile(Form("roots/2DFit_%s/CtauBkg/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCTrue = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_Inclusive_%s.root", DATE.Data(), kineLabel.Data()));
}

void Final2DFit::processDataset() {
  cout << "--- processDataset() ---\n\n";

  // --- set cut label ---
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  TString OS = "recoQQsign==0 &&";
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  kineCut = OS + accCut + kineCut;

  // --- load datasets and PDFs---
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet *)fCErr->Get("dataw_Bkg");
  RooHistPdf *pdfCTAUERR_Tot = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Tot");
  RooHistPdf *pdfCTAUERR_Jpsi = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Jpsi");
  RooHistPdf *pdfCTAUERR_Bkg = (RooHistPdf *)fCErr->Get("pdfCTAUERR_Bkg");
  RooAddPdf *GaussModel_Tot = (RooAddPdf *)fCRes->Get("GaussModel_Tot");
  RooAddPdf *TrueModel_Tot = (RooAddPdf *)fCTrue->Get("TrueModel_Tot");
  RooAddPdf *pdfCTAU_Bkg_Tot = (RooAddPdf *)fCBkg->Get("pdfCTAU_Bkg_Tot");

  // --- declare workspace and import datasets ---
  ws = new RooWorkspace("workspace");
  
  ws->import(*dataw_Bkg);
  ws->import(*dataset); // total
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  // ws->import(*dataw_Sig);
  // ws->import(*GaussModel_Tot);
  ws->import(*TrueModel_Tot);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  // ws->import(*pdfCTAUERR_Bkg);
  ws->import(*pdfCTAU_Bkg_Tot);

  // cut dataset and import process dataset for the final fit
  RooArgSet *argSet = new RooArgSet();
  argSet->add(*(ws->var("ctau3D")));
  argSet->add(*(ws->var("mass")));
  argSet->add(*(ws->var("pt")));
  argSet->add(*(ws->var("y")));
  argSet->add(*(ws->var("weight")));
  argSet->add(*(ws->var("ctau3DRes")));
  argSet->add(*(ws->var("ctau3DErr")));
  argSet->add(*(ws->var("pt1")));
  argSet->add(*(ws->var("pt2")));
  argSet->add(*(ws->var("eta1")));
  argSet->add(*(ws->var("eta2")));
  argSet->add(*(ws->var("recoQQsign"))); // argSet->add(*(ws->var("cBin")) );

  RooDataSet *datasetW = nullptr;
  if (isWeighted)
    datasetW = new RooDataSet("datasetW", "dataset with event weight", *argSet, Import(*dataset), WeightVar(*ws->var("weight")));
  else
    datasetW = new RooDataSet("datasetW", "datawet witout event weight", *argSet, Import(*dataset));
  ws->import(*datasetW);

  RooDataSet *dsTot = (RooDataSet *)datasetW->reduce(*argSet, kineCut.Data());
  ws->import(*dsTot, Rename("dsTot"));
}

void Final2DFit::setVariableRanges() {
  cout << "--- setVariableRanges() ---\n\n";
  double ctauErrMin;
  double ctauErrMax;
  ctauErrMin = ws->var("ctau3DErr")->getMin();
  ctauErrMax = ws->var("ctau3DErr")->getMax();
  // ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  // ws->var("ctau3DErr")->setRange("ctauRange", ctauErrMin, ctauErrMax);
  cout << "CtauErr Min: " << ctauErrMin << ", Max: " << ctauErrMax << endl;

  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setRange("ctauRange", massLow, massHigh);

  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  
  ws->var("ctau3DRes")->setRange(-10, 10);
  // ws->var("ctau3DRes")->setRange("ctauRange", -10, 10);

  // --- print to check ranges ---
  // ws->var("mass")->Print();
  // ws->var("ctau3D")->Print();
  // ws->var("ctau3DErr")->Print();
  // ws->var("ctau3DRes")->Print();
}

void Final2DFit::fixParameters() {
  cout << "--- fixParameters() ---\n\n";
  // --- fix the parameters got from previous procedures ---

  ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass"), *ws->pdf("pdfMASS_Jpsi"), *ws->pdf("pdfMASS_bkg")))->setAttribAll("Constant", kTRUE); // mass result
  ws->pdf("pdfCTAU_Bkg_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")))->setAttribAll("Constant", kTRUE); // ctauBkg result
}

void Final2DFit::buildCtauNpModel() {
  cout << "--- buildCtauNpModel() ---\n\n";
  double lambda = ws->var("lambdaDSS")->getVal();        // RooRealVar
  double lambda2 = ws->function("lambdaDSS2")->getVal(); // RooFormulaVar 는 var() 말고 function()으로!
  double fdss = ws->var("fDSS")->getVal();

  ws->factory(Form("lambdaDSS_test1[%.12f]", lambda)); //
  ws->factory(Form("lambdaDSS_test2[%.12f]", lambda2));
  ws->factory(Form("fDSS1_test[%.2f]", fdss));

  // ws->factory(Form("lambdaDSS_test1[%.2f]", lambda));
  // ws->factory(Form("lambdaDSS_test2[%.2f]", lambda1));
  // ws->factory(Form("lambdaDSS_test3[%.2f]", lambda2));

  // ws->factory(Form("fDSS2_test[%.2f]", fdss1));

  //   cout << lambda2 << endl;
  //   return;

  // NoPR{
  // 1exp
  // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUCOND_JpsiNoPR", "ctau3D", "lambdaDSS", "pdfCTAURES")); //NP
  // 3exp
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test1", "ctau3D", "lambdaDSS_test1", "pdfCTAURES")); // NP
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test2", "ctau3D", "lambdaDSS_test2", "pdfCTAURES")); // NP
  // ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUE_test3", "ctau3D", "lambdaDSS_test3", "pdfCTAURES")); // NP

  //   ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test1", "pdfCTAUTRUE_test2", "fDSS1_test")); //618
  //   ws->factory(Form("AddModel::%s({%s , %s}, %s)", "pdfCTAUCOND_JpsiNoPR", "pdfCTAUTRUE_test12", "pdfCTAUTRUE_test3", "fDSS2_test")); //618

  // auto pdfCTAUTRUE_test12 = new RooAddPdf("pdfCTAUTRUE_test12", "",
  //                                         RooArgList(*ws->pdf("pdfCTAUTRUE_test1"), *ws->pdf("pdfCTAUTRUE_test2")),
  //                                         RooArgList(*ws->var("fDSS1_test")));
  // ws->import(*pdfCTAUTRUE_test12);
  // auto pdfCTAUCOND_JpsiNoPR = new RooAddPdf("pdfCTAUCOND_JpsiNoPR", "",
  //                                           RooArgList(*ws->pdf("pdfCTAUTRUE_test12"), *ws->pdf("pdfCTAUTRUE_test3")),
  //                                           RooArgList(*ws->var("fDSS2_test")));
  // ws->import(*pdfCTAUCOND_JpsiNoPR);

  // auto pdfCTAUTRUE_test12 = new RooAddPdf("pdfCTAUTRUE_test12", "",
  //                                         RooArgList(*ws->pdf("pdfCTAUTRUE_test1"), *ws->pdf("pdfCTAUTRUE_test2")),
  //                                         RooArgList(*ws->var("fDSS1_test")));
  // ws->import(*pdfCTAUTRUE_test12);
  auto pdfCTAUCOND_JpsiNoPR = new RooAddPdf("pdfCTAUCOND_JpsiNoPR", "",
                                            RooArgList(*ws->pdf("pdfCTAUTRUE_test1"), *ws->pdf("pdfCTAUTRUE_test2")),
                                            RooArgList(*ws->var("fDSS1_test")));
  ws->import(*pdfCTAUCOND_JpsiNoPR);
}

void Final2DFit::buildMassCtauModel() {
  cout << "--- buildMassCtauModel() ---\n\n";
  // build mass x ctau models
  // PR
  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_JpsiPR", "pdfCTAURES"));
  ws->factory("b_Jpsi[0.31, 0.1, 0.9]"); // NP fraction for Sig

  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgPR",
                   "pdfCTAU_BkgPR",
                   "pdfMASS_bkg"));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_BkgNoPR",
                   "pdfCTAU_BkgNoPR",
                   "pdfMASS_bkg"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Bkg",
                   "b_Bkg",
                   "pdfCTAUMASS_BkgNoPR",
                   "pdfCTAUMASS_BkgPR"));

  // product Punzi terms (signal PR, NP). FYI: Bkg was handled in ctauBkg fit
  RooProdPdf pdfJpsiPR("pdfCTAU_JpsiPR", "", *ws->pdf("pdfCTAUERR_Jpsi"),
                       Conditional(*ws->pdf("pdfCTAUCOND_JpsiPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiPR);
  RooProdPdf pdfJpsiNoPR("pdfCTAU_JpsiNoPR", "", *ws->pdf("pdfCTAUERR_Jpsi"),
                         Conditional(*ws->pdf("pdfCTAUCOND_JpsiNoPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfJpsiNoPR);

  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiPR",
                   "pdfCTAU_JpsiPR",
                   "pdfMASS_Jpsi"));
  ws->factory(Form("PROD::%s(%s, %s)", "pdfCTAUMASS_JpsiNoPR",
                   "pdfCTAU_JpsiNoPR",
                   "pdfMASS_Jpsi"));
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUMASS_Jpsi",
                   "b_Jpsi",
                   "pdfCTAUMASS_JpsiNoPR",
                   "pdfCTAUMASS_JpsiPR"));
  RooAbsPdf *themodel = NULL;

  // double njpsi = ws->var("N_Jpsi")->getVal();
  // ws->factory(Form("N_Jpsi[%.3f, %.3f, %.3f]",njpsi, njpsi*0.9, njpsi*1.1));

  ws->var("N_Jpsi")->setConstant(kTRUE);
  ws->var("N_Bkg")->setConstant(kTRUE);
  // ws->var("sigma_1_A")->setConstant(false);

  themodel = new RooAddPdf("pdfCTAUMASS_Tot", "pdfCTAUMASS_Tot",
                           RooArgList(*ws->pdf("pdfCTAUMASS_Bkg"), *ws->pdf("pdfCTAUMASS_Jpsi")),
                           RooArgList(*ws->var("N_Bkg"), *ws->var("N_Jpsi")));
  ws->import(*themodel);
}

void Final2DFit::do2dFit() {
  cout << "--- do2dFit() ---\n\n";
  // Legacy - I thinkg it's a meaningless cloning
  RooDataSet *dsToFit = (RooDataSet *)ws->data("dsTot");
  ws->import(*dsToFit, Rename("dsToFit"));

  cout << "##############START TOTAL CTAU FIT############" << endl;
  bool isWeighted = ws->data("dsTot")->isWeighted();

  // RooFitResult* fitResult = ws->pdf("pdfCTAUMASS_Tot")->fitTo(*dsToFit, Extended(kTRUE), ExternalConstraints(*ws->set("ConstrainPdfList")), NumCPU(nCPU), SumW2Error(isWeighted), PrintLevel(3), Save());
  fitResult = ws->pdf("pdfCTAUMASS_Tot")->fitTo(*dsToFit, Extended(true), NumCPU(nCPU), Strategy(1), SumW2Error(isWeighted), PrintLevel(-1), Save(), ConditionalObservables(RooArgSet(*ws->var("ctau3DErr"))), Range("ctauRange"), RecoverFromUndefinedRegions(1.));

  // --- note for fit options ---
  // SumW2Error(isWeighted)
  // AsymptoticError
  
  // ws->import(*fitResult, "fitResult_pdfCTAUMASS_Tot");
}

void Final2DFit::drawCtauPull() {
  cout << "--- drawCtauPull() ---\n\n";
  // --- calculate ratio for normalization -> Legacy. I think we don't need it ---
  // double normDSTot = ws->data("dsToFit")->sumEntries()/ws->data("dsTot")->sumEntries();
  double normDSTot = 1; // one means no user-custom scale

  // buld canvas and pad
  TCanvas *c_G = new TCanvas("canvas_G", "My plots", 1108, 4, 550, 520);
  c_G->cd();

  TPad *pad_G_1 = new TPad("pad_G_1", "pad_G_1", 0, 0.16, 0.98, 1.0);
  pad_G_1->SetTicks(1, 1);
  pad_G_1->Draw();
  pad_G_1->cd();
  gPad->SetLogy();

  // build frame
  RooPlot *myPlot_G = ws->var("ctau3D")->frame(Bins(nCtauBins), Range("ctauRange")); // bins
  // RooPlot* myPlot_H = ws->var("mass")->frame(Bins(nMassBin), Range(2.6,3.5)); // bins
  myPlot_G->SetTitle("");

  // --- plotOn ---
  RooPlot *myPlot2_G = (RooPlot *)myPlot_G->Clone(); // clone frame to avoid memory error
  myPlot2_G->updateNormVars(RooArgSet(*ws->var("ctau3D")));

  // plot dataset
  ws->data("dsToFit")->plotOn(myPlot_G, Name("dataHist_ctau"), DataError(RooAbsData::SumW2));

  // plot pdfs
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_Tot"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), FillStyle(1001), FillColor(kViolet + 6), VLines(), DrawOption("LF"), NumCPU(nCPU), LineColor(kBlack), Precision(1e-4), Normalization(normDSTot, RooAbsReal::RelativeExpected));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_Bkg"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_Bkg"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), FillStyle(1001), FillColor(kAzure - 9), VLines(), DrawOption("F"), NumCPU(nCPU));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_JpsiPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kRed + 3), NumCPU(nCPU));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_JpsiNoPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiNoPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kGreen + 3), NumCPU(nCPU));
  ws->data("dsToFit")->plotOn(myPlot2_G, Name("data_Ctau"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));

  // --- find good y-max ---
  myPlot2_G->GetYaxis()->SetRangeUser(10e-2, 10e8);
  TH1 *hCtau = ws->data("dsToFit")->createHistogram("histCtau", *ws->var("ctau3D"), Binning(myPlot_G->GetNbinsX(), myPlot_G->GetXaxis()->GetXmin(), myPlot_G->GetXaxis()->GetXmax()));
  Double_t YMaxCtau = hCtau->GetBinContent(hCtau->GetMaximumBin());
  Double_t YMinCtau = 1e99;
  for (int i = 1; i <= hCtau->GetNbinsX(); i++)
    if (hCtau->GetBinContent(i) > 0)
      YMinCtau = min(YMinCtau, hCtau->GetBinContent(i));
  Double_t YupCtau(0.), YdownCtau(0.);
  YupCtau = YMaxCtau * TMath::Power((YMaxCtau / 0.1), 0.5);
  
  YdownCtau = 0.1; // fix y-min

  // --- other cosmetic ---
  myPlot2_G->GetYaxis()->SetRangeUser(YdownCtau, YupCtau);
  myPlot2_G->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_G->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_G->SetFillStyle(4000);
  myPlot2_G->GetXaxis()->SetLabelSize(0);
  myPlot2_G->GetXaxis()->SetTitleSize(0);
  myPlot2_G->Draw();

  // --- draw legend ---
  TLegend *leg_G = new TLegend(text_x+0.033 + 0.25, text_y+0.09 + 0.05, text_x+0.033 + 0.38, text_y+0.09 - 0.11);
  leg_G->SetTextSize(text_size);
  leg_G->SetTextFont(43);
  leg_G->SetBorderSize(0);
  leg_G->AddEntry(myPlot2_G->findObject("data_Ctau"), "Data", "pe");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Tot"), "Total fit", "fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Bkg"), "Background", "fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiPR"), "J/#psi Prompt", "l");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiNoPR"), "J/#psi Non-Prompt", "l");
  leg_G->Draw("same");

  // --- print latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x+0.033, text_y+0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x+0.033, text_y+0.09 - y_diff * 2, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff, text_color, text_size);
  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff * 2, text_color, text_size);
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f", ws->var("b_Jpsi")->getVal(), ws->var("b_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff * 3, text_color, text_size);

  // --- pull padd ---
  TPad *pad_G_2 = new TPad("pad_G_2", "pad_G_2", 0, 0.006, 0.98, 0.227);
  RooPlot *frameTMP_G = (RooPlot *)myPlot2_G->Clone("TMP_G");
  RooHist *hpull_G;

  pullDist(ws, pad_G_2, c_G, frameTMP_G, hpull_G, "data_Ctau", "Ctau_Tot", "ctau3D", nCtauBins, -4, 6.0, "#font[12]{l}_{J/#psi} (mm)");
  printChi2(ws, pad_G_2, frameTMP_G, fitResult, "ctau3D", "data_Ctau", "Ctau_Tot", nCtauBins);
  pad_G_2->Update();

  c_G->Update();
  c_G->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}

void Final2DFit::drawMass() {
  cout << "--- drawMass() ---\n\n";
  double normDSTot = 1;

  // build canvas and pad
  TCanvas *c_H = new TCanvas("c_H", "My plots", 800, 800);
  c_H->cd();

  TPad *pad_H_1 = new TPad("pad_H_1", "pad_H_1", 0.0, 0.3, 1.0, 1.0);
  pad_H_1->SetBottomMargin(0);
  pad_H_1->SetTicks(1, 1);
  pad_H_1->Draw();
  pad_H_1->cd();
  gPad->SetLogy();

  // buld frame
  RooPlot *myPlot_H = ws->var("mass")->frame(Bins(nMassBin));
  myPlot_H->SetTitle("");

  // --- plotting ---
  RooPlot *myPlot2_H = (RooPlot *)myPlot_H->Clone(); // clone the frame to avoid memory error later
  myPlot2_H->updateNormVars(RooArgSet(*ws->var("mass")));

  // dataset
  ws->data("dsToFit")->plotOn(myPlot_H, Name("dataHist_mass"));
  
  // pdfs
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H, Name("Mass_Bkg"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_Bkg"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), FillStyle(1001), FillColor(kAzure - 9), VLines(), DrawOption("LCF"), NumCPU(nCPU), LineStyle(kDashed));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H, Name("Mass_JpsiPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiPR"), *ws->pdf("pdfCTAUMASS_Bkg"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kRed + 3), LineStyle(1), Precision(1e-4), NumCPU(nCPU));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H, Name("Mass_JpsiNoPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiNoPR"), *ws->pdf("pdfCTAUMASS_Bkg"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kGreen + 3), LineStyle(1), Precision(1e-4), NumCPU(nCPU));
  ws->data("dsToFit")->plotOn(myPlot2_H, Name("data_Mass"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_H, Name("Mass_Tot"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), NumCPU(nCPU), LineColor(kBlack));

  // --- find y-max value ---
  TH1 *h = ws->data("dsToFit")->createHistogram("hist2", *ws->var("mass"), Binning(myPlot_H->GetNbinsX(), myPlot_H->GetXaxis()->GetXmin(), myPlot_H->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBiAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));

  // --- other cosmetics ---
  myPlot2_H->GetYaxis()->SetRangeUser(Ydown, Yup);
  myPlot2_H->GetXaxis()->SetRangeUser(massLow, massHigh);
  myPlot2_H->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-} (GeV/c^{2})}");
  myPlot2_H->SetFillStyle(4000);
  myPlot2_H->GetXaxis()->SetLabelSize(0);
  myPlot2_H->GetXaxis()->SetTitleSize(0);
  myPlot2_H->Draw();

  // -- draw legends ---
  TLegend *leg_H = new TLegend(text_x+0.033 + 0.25, text_y+0.09 + 0.03, text_x+0.033 + 0.38, text_y+0.09 - 0.17);
  leg_H->SetTextSize(text_size);
  leg_H->SetTextFont(43);
  leg_H->SetBorderSize(0);
  leg_H->AddEntry(myPlot2_H->findObject("data_Mass"), "Data", "pe");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_Tot"), "Total fit", "fl");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_Bkg"), "Background", "fl");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_JpsiPR"), "J/#psi Prompt", "l");
  leg_H->AddEntry(myPlot2_H->findObject("Mass_JpsiNoPR"), "J/#psi Non-Prompt", "l");
  leg_H->Draw("same");

  // --- print latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x+0.033, text_y+0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x+0.033, text_y+0.09 - y_diff * 2, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 - y_diff, text_color, text_size);
  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x+0.033 + 0.5, text_y+0.09 - y_diff * 2, text_color, text_size);
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f", ws->var("b_Jpsi")->getVal(), ws->var("b_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 - y_diff * 3, text_color, text_size);
  
  // --- pull pad ---
  TPad *pad_H_2 = new TPad("pad_H_2", "pad_H_2", 0.0, 0.0, 1.0, 0.3);
  RooPlot *frameTMP_H = (RooPlot *)myPlot2_H->Clone("TMP_H");
  RooHist *hpull_H;

  pullDist(ws, pad_H_2, c_H, frameTMP_H, hpull_H, "data_Mass", "Mass_Tot", "mass", nMassBin, massLow, massHigh, "m_{#mu^{+}#mu^{-} (GeV/c^{2})}");
  pad_H_2->Update();

  printChi2(ws, pad_H_2, frameTMP_H, fitResult, "mass", "data_Mass", "Mass_Tot", nMassBin, false);

  cout << "############################################################################" << endl;

  c_H->Update();

  // --- save ---
  c_H->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}

void Final2DFit::drawCtauRatio() {
  cout << "--- drawCtauRatio() ---\n\n";
  double normDSTot = 1;

  // ctau canvas
  TCanvas *c_G = new TCanvas("canvas_G", "My plots", 1108, 4, 550, 520);
  c_G->cd();
  TPad *pad_G_1 = new TPad("pad_G_1", "pad_G_1", 0, 0.16, 0.98, 1.0);
  pad_G_1->SetTicks(1, 1);
  pad_G_1->Draw();
  pad_G_1->cd();
  gPad->SetLogy();
  RooPlot *myPlot_G = ws->var("ctau3D")->frame(Bins(nCtauBins), Range("ctauRange")); // bins
  // RooPlot* myPlot_H = ws->var("mass")->frame(Bins(nMassBin), Range(2.6,3.5)); // bins
  myPlot_G->SetTitle("");

  c_G->cd();
  c_G->SetLogy();
  pad_G_1->cd();

  // drawing part
  RooPlot *myPlot2_G = (RooPlot *)myPlot_G->Clone();
  myPlot2_G->updateNormVars(RooArgSet(*ws->var("ctau3D")));

  ws->data("dsToFit")->plotOn(myPlot_G, Name("dataHist_ctau"), DataError(RooAbsData::SumW2));
  //   myPlot2_G->updateNormVars(RooArgSet(*ws->var("mass"), *ws->var("ctau3D"), *ws->var("ctau3DErr")));
  // ws->data("dsToFit")->plotOn(myPlot2_G,Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));

  // ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_Tot"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::NumEvent), FillStyle(1001), FillColor(kViolet + 6), VLines(), DrawOption("LF"), NumCPU(nCPU), LineColor(kBlack));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_Tot"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), FillStyle(1001), FillColor(kViolet + 6), VLines(), DrawOption("LF"), NumCPU(nCPU), LineColor(kBlack), Precision(1e-4), Normalization(normDSTot, RooAbsReal::RelativeExpected));

  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_Bkg"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_Bkg"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), FillStyle(1001), FillColor(kAzure - 9), VLines(), DrawOption("F"), NumCPU(nCPU));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_JpsiPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kRed + 3), NumCPU(nCPU));
  ws->pdf("pdfCTAUMASS_Tot")->plotOn(myPlot2_G, Name("Ctau_JpsiNoPR"), Components(RooArgSet(*ws->pdf("pdfCTAUMASS_JpsiNoPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dsToFit"), kTRUE), Normalization(normDSTot, RooAbsReal::RelativeExpected), LineColor(kGreen + 3), NumCPU(nCPU));
  ws->data("dsToFit")->plotOn(myPlot2_G, Name("data_Ctau"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack));

  myPlot2_G->GetYaxis()->SetRangeUser(10e-2, 10e8);
  TH1 *hCtau = ws->data("dsToFit")->createHistogram("histCtau", *ws->var("ctau3D"), Binning(myPlot_G->GetNbinsX(), myPlot_G->GetXaxis()->GetXmin(), myPlot_G->GetXaxis()->GetXmax()));
  Double_t YMaxCtau = hCtau->GetBinContent(hCtau->GetMaximumBin());
  Double_t YMinCtau = 1e99;
  for (int i = 1; i <= hCtau->GetNbinsX(); i++)
    if (hCtau->GetBinContent(i) > 0)
      YMinCtau = min(YMinCtau, hCtau->GetBinContent(i));
  Double_t YupCtau(0.), YdownCtau(0.);
  YupCtau = YMaxCtau * TMath::Power((YMaxCtau / 0.1), 0.5);
  YdownCtau = 0.1;
  myPlot2_G->GetYaxis()->SetRangeUser(YdownCtau, YupCtau);
  myPlot2_G->GetXaxis()->SetRangeUser(-4, 6);
  myPlot2_G->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  myPlot2_G->SetFillStyle(4000);
  myPlot2_G->GetXaxis()->SetLabelSize(0);
  myPlot2_G->GetXaxis()->SetTitleSize(0);
  myPlot2_G->Draw();
  TLegend *leg_G = new TLegend(text_x+0.033 + 0.25, text_y+0.09 + 0.05, text_x+0.033 + 0.38, text_y+0.09 - 0.11);
  leg_G->SetTextSize(text_size);
  leg_G->SetTextFont(43);
  leg_G->SetBorderSize(0);
  leg_G->AddEntry(myPlot2_G->findObject("data_Ctau"), "Data", "pe");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Tot"), "Total fit", "fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_Bkg"), "Background", "fl");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiPR"), "J/#psi Prompt", "l");
  leg_G->AddEntry(myPlot2_G->findObject("Ctau_JpsiNoPR"), "J/#psi Non-Prompt", "l");
  leg_G->Draw("same");
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x+0.033, text_y+0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x+0.033, text_y+0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x+0.033, text_y+0.09 - y_diff * 2, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff, text_color, text_size);
  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff * 2, text_color, text_size);
  drawText(Form("b_{J/#psi} = %.4f #pm %.4f", ws->var("b_Jpsi")->getVal(), ws->var("b_Jpsi")->getError()), text_x+0.033 + 0.5, text_y+0.09 + 0.05 - y_diff * 3, text_color, text_size);

  TCanvas *c_ratio = new TCanvas("c_ratio", "My plots", 800, 800);
  c_ratio->cd();
  pad_G_1->Draw(); // reuse previous ctau3D distribution
  TPad *ratioPad = new TPad("ratioPad", "", 0, 0.006, 0.98, 0.227);
  c_ratio->cd();
  // ratioPad->SetLeftMargin(0.15);
  // ratioPad->SetRightMargin(0.07);
  ratioPad->SetTopMargin(0); // Upper and lower plot are joined
  // ratioPad->SetTopMargin(0.05);
  ratioPad->SetBottomMargin(0.35);
  // ratioPad->SetTicks(1, 1);
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->Draw();
  ratioPad->cd();

  RooHist *h_dataPoints = (RooHist *)myPlot2_G->findObject("data_Ctau");
  RooCurve *curveRatio = (RooCurve *)myPlot2_G->findObject("Ctau_Tot");

  // RooPlot *ratioFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(myPlot2_G->GetXaxis()->GetXmin(), myPlot2_G->GetXaxis()->GetXmax()));
  RooPlot *ratioFrame = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(-4, 6));
  ratioFrame->SetTitle(" ");

  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  int point_idx = 0;
  double x_min = ws->var("ctau3D")->getMin("ctauRange");
  double x_max = ws->var("ctau3D")->getMax("ctauRange");
  for (int i = 0; i < h_dataPoints->GetN(); ++i)
  {
    double x, y;
    h_dataPoints->GetPoint(i, x, y);
    double model_val = curveRatio->Eval(x);
    if (x < x_min || x > x_max)
      continue;

    if (model_val > 1e-9 && y > 0)
    {
      double ratio = y / model_val;
      g_ratio->SetPoint(point_idx, x, ratio);
      double err_y_low = h_dataPoints->GetErrorYlow(i);
      double err_y_high = h_dataPoints->GetErrorYhigh(i);
      g_ratio->SetPointError(point_idx, 0, 0, err_y_low / model_val, err_y_high / model_val);
      point_idx++;
    }
  }

  ratioFrame->SetTitle("");
  ratioFrame->SetTitleSize(0);
  ratioFrame->GetYaxis()->SetTitleOffset(0.3);
  ratioFrame->GetYaxis()->SetTitle("Data / Fit");
  ratioFrame->GetYaxis()->SetTitleSize(0.08);
  ratioFrame->GetYaxis()->SetLabelSize(0.08);
  ratioFrame->GetYaxis()->SetRangeUser(0, 3); // 0, 3
  ratioFrame->GetYaxis()->CenterTitle();

  ratioFrame->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
  ratioFrame->GetXaxis()->SetTitleOffset(1.4);
  ratioFrame->GetXaxis()->SetLabelOffset(0.04);
  ratioFrame->GetXaxis()->SetLabelSize(0.08);
  ratioFrame->GetXaxis()->SetTitleSize(0.08);
  ratioFrame->GetXaxis()->CenterTitle();

  ratioFrame->GetYaxis()->SetTickSize(0.04);
  ratioFrame->GetYaxis()->SetNdivisions(404);
  ratioFrame->GetXaxis()->SetTickSize(0.03);
  ratioFrame->Draw();

  g_ratio->SetMarkerStyle(kFullCircle);
  g_ratio->SetMarkerSize(0.7);
  g_ratio->Draw("P SAME");

  double x_min_line = ratioFrame->GetXaxis()->GetXmin();
  double x_max_line = ratioFrame->GetXaxis()->GetXmax();
  TLine *line_at_1 = new TLine(x_min_line, 1.0, x_max_line, 1.0);
  line_at_1->SetLineColor(kRed);
  line_at_1->SetLineStyle(kDashed);
  line_at_1->Draw("same");
  c_ratio->SaveAs(Form("figs/2DFit_%s/Final/2DFit_Ctau_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d_ratio.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}

void Final2DFit::saveOutput() {
  cout << "--- saveOutput() ---\n\n";
  // for pT reweighting -> but fitResult is good enough. Isn't it?
  TH1D *outh = new TH1D("2DfitResults", "fit result", 20, 0, 20);

  float temp = ws->var("b_Jpsi")->getVal();
  float temperr = ws->var("b_Jpsi")->getError();

  outh->SetBinContent(1, temp);
  outh->SetBinError(1, temperr);

  // std::cout << std::setprecision(20);
  // std::cout << "b_Jpsi = " << ws->var("b_Jpsi")->getVal()
  //         << " +/- " << ws->var("b_Jpsi")->getError() << std::endl;
  fitResult->Print("v");
  // const TMatrixDSym &cor = fitResult->correlationMatrix();
  // cor.Print();
  TFile *outFile = new TFile(Form("roots/2DFit_%s/Final/2DFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");
  outFile->cd();
  
  // ws->Write();
  outh->Write();
  ws->pdf("pdfCTAUMASS_Tot")->Write();
  
  outFile->Close();

  // --- print fit result ---
  fitResult->Print("V");
  fitResult->correlationMatrix().Print();
}

// =================================
// ===== legacy codes or notes =====
// =================================

// ===== legacy headers =====
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

// ===== legacy code for Gaussian constraint + fixing the parameters =====
// {
//   // ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass")))->setAttribAll("Constant", kTRUE);
//   std::vector<std::string> objs = {"Bkg", "Jpsi"};
//   RooArgSet pdfList = RooArgSet("ConstraionPdfList");
//   for (auto obj : objs)
//   {
//     if (ws->var(Form("N_%s", obj.c_str())))
//     {
//       ws->factory(Form("Gaussian::%s_Gauss(%s,%s_Mean[%f],%s_Sigma[%f])",
//                        Form("N_%s", obj.c_str()), Form("N_%s", obj.c_str()),
//                        Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getValV(),
//                        Form("N_%s", obj.c_str()), ws->var(Form("N_%s", obj.c_str()))->getError()));

//       pdfList.add(*ws->pdf(Form("N_%s_Gauss", obj.c_str())), kFALSE);
//       std::cout << "[INFO] Constraining N_" << obj << " with Mean : " << ws->var(Form("N_%s_Mean", obj.c_str()))->getVal()
//                 << " and Sigma: " << ws->var(Form("N_%s_Sigma", obj.c_str()))->getVal() << std::endl;
//     }
//   }
//   ws->defineSet("ConstrainPdfList", pdfList);

//   ws->pdf("pdfCTAURES")->getParameters(RooArgSet(*ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
//   cout << ws->pdf("pdfCTAURES")->getVal() << endl;
//   ws->pdf("pdfCTAU_JpsiPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
//   ws->pdf("pdfCTAU_BkgPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);
//   ws->pdf("pdfCTAU_BkgNoPR")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr")))->setAttribAll("Constant", kTRUE);

//   RooArgSet *params = (RooArgSet *)ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("mass"), *ws->var("ctau3DErr")));
//   ws->saveSnapshot(("pdfCTAUMASS_Tot_parIni"), *params, kTRUE);
//   delete params;

//   RooArgSet *newpars = (RooArgSet *)ws->pdf("pdfCTAUMASS_Tot")->getParameters(RooArgSet(*ws->var("ctau3D"), *ws->var("ctau3DErr"), *ws->var("ctau3DRes"), *ws->var("mass")));
// }