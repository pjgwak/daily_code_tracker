#include "CtauBkgFit.h"
#include <iostream>
#include "TStyle.h" // gStyle
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
#include "/work/pjgwak/pol24/headers/polarizationUtilities.h"

using std::cout; using std::endl;
using std::min;
using namespace RooFit;

CtauBkgFit::CtauBkgFit(float ptLow, float ptHigh,
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

CtauBkgFit::~CtauBkgFit()
{
  delete fMass;
  delete fCErr;
  delete fCRes;
  
  delete fitResult;

  delete ws;
}

void CtauBkgFit::init()
{
  cout << "--- init() ---\n\n";
  setLabels();
  makeOutputFolder();
  turnOffRooFitMessage();
  openInput();
  processDataset();
  setVariableRanges();
  fixParameters();
  buildCtauFitModel();
  buildCtauCondModel();
}

void CtauBkgFit::run()
{
  cout << "--- run() ---\n\n";
  doFit();
  rootlogon(isLogon);
  drawCtauPull();
  saveOutput();
}

void CtauBkgFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  if (PRw == 1) fname = "PR";
  else if (PRw == 2) fname = "NP";
}

void CtauBkgFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauBkg", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauBkg", DATE.Data()), kTRUE);
}

void CtauBkgFit::turnOffRooFitMessage()
{
  cout << "--- turnOffRooFitMessage() ---\n\n";
  RooMsgService::instance().getStream(1).removeTopic(RooFit::ObjectHandling);

  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(1).removeTopic(Tracing);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);

  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
}

void CtauBkgFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
  fCRes = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}

void CtauBkgFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // load dataset and pdfs
  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Bkg = (RooDataSet *)fCErr->Get("dataw_Bkg");
  RooAddPdf *GaussModel_Tot = (RooAddPdf *)fCRes->Get("GaussModel_Tot");
  RooAddPdf *pdfCTAUERR_Bkg = (RooAddPdf *)fCErr->Get("pdfCTAUERR_Bkg");

  // declare workspace and import
  ws = new RooWorkspace("workspace");
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Bkg);
  ws->import(*GaussModel_Tot);
  ws->import(*pdfCTAUERR_Bkg);
}

void CtauBkgFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  double ctauErrMin, ctauErrMax;
  ws->data("dataw_Bkg")->getRange(*ws->var("ctau3DErr"), ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  ws->var("ctau3DErr")->setRange("ctauErrRange", ctauErrMin, ctauErrMax);

  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  
  // --- print ranges ---
  // ws->var("ctau3D")->Print();
  // ws->var("ctau3DErr")->Print();
}

void CtauBkgFit::fixParameters()
{
  cout << "--- fixParameters() ---\n\n";
  // // --- call this function in a doFit() ---

  // // fix the parameters from previous procedures

  // // ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass")))->setAttribAll("Constant", kTRUE);
  // ws->pdf("pdfMASS_Tot")->getParameters(RooArgSet(*ws->var("mass"), *ws->pdf("pdfMASS_Jpsi"), *ws->pdf("pdfMASS_bkg")))->setAttribAll("Constant", kTRUE);
  // // ws->pdf("GaussModelCOND_ctauRes")->getParameters(
  // //     RooArgSet(*ws->var("ctau1_CtauRes"),*ws->var("ctau2_CtauRes"),*ws->var("ctau3_CtauRes"),
  // //       *ws->var("s1_CtauRes"), *ws->var("rS21_CtauRes"), *ws->var("rS32_CtauRes"),
  // //       *ws->var("f_CtauRes"), *ws->var("f2_CtauRes")
  // //       ))->setAttribAll("Constant", kTRUE);
  // ws->pdf("pdfCTAU_Bkg_Tot")->getParameters(
  //                               // RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"),  *ws->pdf("pdfCTAUCOND_Bkg"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")
  //                               RooArgSet(*ws->var("ctau3D"), *ws->pdf("pdfCTAU_BkgPR"), *ws->pdf("pdfCTAU_BkgNoPR"), *ws->pdf("pdfCTAUCOND_BkgPR"), *ws->pdf("pdfCTAUCOND_BkgNoPR")))
  //     ->setAttribAll("Constant", kTRUE);
}

void CtauBkgFit::buildCtauFitModel()
{
  cout << "--- buildCtauFitModel() ---\n\n";
  // --- set parameters ---
  ws->factory("zeroMean[0.0]");
  ws->factory("b_Bkg[0.1, 0., 1.2]"); // NP fraction for bkg

  RooRealVar nBkgNp("nBkgNp", "", 4000, 3500, 5000);
  ws->import(nBkgNp);

  // producted fractions
  RooRealVar q1("q1", "", 0.4, 0, 0.99999);
  RooRealVar q2("q2", "", 0.08, 0.01, 0.12);
  RooRealVar q3("q3", "", 0.2, 0, 0.99999);
  RooRealVar q4("q4", "", 0.2, 0, 0.99999);

  RooFormulaVar frac1("frac1", "@0", RooArgList(q1));
  RooFormulaVar frac2("frac2", "(1-@0)*@1", RooArgList(q1, q2));
  RooFormulaVar frac3("frac3", "(1-@0)*(1-@1)*@2", RooArgList(q1, q2, q3));
  ws->import(frac1);
  ws->import(frac2);
  ws->import(frac3);

  // Base lambda
  RooRealVar lambdaP("lambdaP", "", 0.1, 0.01, 0.5);
  RooRealVar lambdaC("lambdaC", "", 0.2, 0.01, 5.0);
  RooRealVar lambdaN("lambdaN", "", 0.2, 0.01, 5.0);

  ws->import(lambdaP);
  ws->import(lambdaC);
  ws->import(lambdaN);

  // ==== Todo: make function to build ctauBkg Res model ====
  // --- fix Res parameters -> Todo: move to fix function ---
  ws->var("ctau1_CtauRes")->setConstant(kTRUE);
  ws->var("s1_CtauRes")->setConstant(kTRUE);
  ws->var("ctau2_CtauRes")->setConstant(kTRUE);
  ws->var("rS21_CtauRes")->setConstant(kTRUE);
  ws->var("f_CtauRes")->setConstant(kTRUE);

  // cout << ws->var("ctau1_CtauRes")->getVal() << endl;
  // cout << ws->var("f_CtauRes")->getVal() << "+/-" << ws->var("f_CtauRes")->getError() << endl;
  // cout << ws->var("s1_CtauRes")->getVal() << endl;

  // --- make res model ---
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes1", "ctau3D",
                   "ctau1_CtauRes", //"ctau1_CtauRes",
                   "s1_CtauRes",
                   "zeroMean",
                   "ctau3DErr"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s, %s, %s)", "ctauRes2", "ctau3D",
                   "ctau2_CtauRes", //"ctau2_CtauRes",
                   "s2_CtauRes",
                   "zeroMean",
                   "ctau3DErr"));

  ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "pdfCTAURES", "ctauRes1", "ctauRes2", "f_CtauRes")); // 618
  // ==== End of ctauBkg res ====

  // --- build ctauBkg model ---
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUDSS", "ctau3D", "lambdaP", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", "pdfCTAUDF", "ctau3D", "lambdaN", "pdfCTAURES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", "pdfCTAUDDS", "ctau3D", "lambdaC", "pdfCTAURES"));

  RooAddPdf *pdfCTAU1 = new RooAddPdf("pdfCTAU1", "DF+SS model",
                                      RooArgList(*ws->pdf("pdfCTAUDF"), *ws->pdf("pdfCTAUDSS")),
                                      RooArgList(*ws->function("frac2")));
  ws->import(*pdfCTAU1);
  RooAddPdf *pdfCTAUCOND_BkgNoPR = new RooAddPdf("pdfCTAUCOND_BkgNoPR", "full background model",
                                                 RooArgList(*ws->pdf("pdfCTAU1"), *ws->pdf("pdfCTAUDDS")),
                                                 RooArgList(*ws->function("frac1")));
  ws->import(*pdfCTAUCOND_BkgNoPR);

  ws->factory(Form("SUM::%s(%s)", "pdfCTAUCOND_BkgPR", "pdfCTAURES")); // PR 618

  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAUCOND_Bkg", "b_Bkg", "pdfCTAUCOND_BkgNoPR", "pdfCTAUCOND_BkgPR"));

  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "pdfTot_Bkg", "pdfCTAUCOND_Bkg", "nBkgNp")); // nBkgNp is number of bkg events in dataw_Bkg dataset
}

void CtauBkgFit::buildCtauCondModel()
{
  cout << "--- buildCtauCondModel() ---\n\n";
  // it's used for plotting in this code and will be used for fit in the final fit
  RooProdPdf pdfPR("pdfCTAU_BkgPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfPR);
  RooProdPdf pdfNoPR("pdfCTAU_BkgNoPR", "", *ws->pdf("pdfCTAUERR_Bkg"), Conditional(*ws->pdf("pdfCTAUCOND_BkgNoPR"), RooArgList(*ws->var("ctau3D"))));
  ws->import(pdfNoPR);
  ws->factory(Form("SUM::%s(%s*%s, %s)", "pdfCTAU_Bkg", "b_Bkg", "pdfCTAU_BkgNoPR", "pdfCTAU_BkgPR"));
  RooAbsPdf *pdfCTAU_Bkg_Tot = new RooAddPdf("pdfCTAU_Bkg_Tot", "pdfCTAU_Bkg_Tot", RooArgList(*ws->pdf("pdfCTAU_Bkg")), RooArgList(*ws->var("nBkgNp")));
  ws->import(*pdfCTAU_Bkg_Tot);
}

void CtauBkgFit::doFit()
{
  cout << "--- doFit() ---\n\n";
  // --- reduce dataset with new ctauRange ---
  TH1D *hTot = (TH1D *)ws->data("dataw_Bkg")->createHistogram(("hTot"), *ws->var("ctau3D"), Binning(nCtauBins, ctauLow, ctauHigh));
  ctauMin = hTot->GetBinLowEdge(hTot->FindFirstBinAbove(1, 1));
  ctauMax = hTot->GetBinLowEdge(hTot->FindLastBinAbove(2, 1)) + hTot->GetBinWidth(hTot->FindLastBinAbove(2, 1));

  // --- Todo: user custom range ---
  // ctauMin = -2.0; // default value
  // ctauMax = 4.0; // default value

  RooDataSet *dataToFit = (RooDataSet *)ws->data("dataw_Bkg")->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f", ctauMin, ctauMax));
  //  RooDataSet* dataToFit = (RooDataSet*)dataw_Bkg->reduce(Form("ctau3D>=%.f&&ctau3D<=%.f",-2.0, 4.0))->Clone("dataw_Bkg");
  ws->import(*dataToFit, Rename("dataToFit"));

  ws->var("ctau3D")->setRange("fitRange", ctauMin, ctauMax);

  // --- fit ---
  bool isWeighted = ws->data("dataw_Bkg")->isWeighted();
  fitResult = ws->pdf("pdfTot_Bkg")->fitTo(*dataToFit, Save(), Extended(kTRUE), NumCPU(16), PrintLevel(0), AsymptoticError(isWeighted), RecoverFromUndefinedRegions(1.), PrintEvalErrors(-1), Timer(true), Range("fitRange"));
  
  // --- fit options ---
  // SumW2Error(isWeighted)
  // ConditionalObservables(RooArgSet(*ws->var("ctau3DErr"))) 
}

void CtauBkgFit::drawCtauPull()
{
  cout << "--- drawCtauPull() ---\n\n";
  // --- make canvans and pad ---
  TCanvas *c_E = new TCanvas("canvas_E", "My plots", 800, 800);
  c_E->cd();

  TPad *pad_E_1 = new TPad("pad_E_1", "pad_E_1", 0, 0.16, 0.98, 1.0);
  pad_E_1->SetTicks(1, 1);
  pad_E_1->Draw();
  pad_E_1->cd();
  gPad->SetLogy();

  // --- make a frame---
  RooPlot *myPlot_E = ws->var("ctau3D")->frame(Bins(nCtauBins), Range(ctauLow, ctauHigh)); // bins
  myPlot_E->SetTitle("");
  
  RooPlot *myPlot2_E = (RooPlot *)myPlot_E->Clone(); // copy to avoid memory error
  myPlot2_E->updateNormVars(RooArgSet(*ws->var("ctau3D")));
  ws->pdf("pdfCTAU_Bkg_Tot")->setNormRange("fitRange");

  // --- plotOn ---
  // data
  ws->data("dataToFit")->plotOn(myPlot2_E, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlue), LineColor(kBlue), MarkerSize(0.7));

  // pdf
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E, Name("ctauBkg_Tot"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataToFit"), kTRUE), FillStyle(1001), FillColor(kAzure - 9), VLines(), DrawOption("LCF"), Precision(1e-4));
  if (ws->pdf("pdfCTAU_BkgPR"))
  {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E, Name("BKGPR"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataToFit"), kTRUE), LineColor(kRed + 2));
  }
  if (ws->pdf("pdfCTAU_BkgNoPR"))
  {
    ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E, Name("BKGNoPR"), Components(RooArgSet(*ws->pdf("pdfCTAU_BkgNoPR"))), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataToFit"), kTRUE), LineColor(kOrange + 10), Precision(1e-4));
  }
  ws->pdf("pdfCTAU_Bkg_Tot")->plotOn(myPlot2_E, Name("PDF"), ProjWData(RooArgSet(*ws->var("ctau3DErr")), *ws->data("dataToFit"), kTRUE), LineColor(kBlack), Precision(1e-4));

  // redraw dataset at top
  ws->data("dataToFit")->plotOn(myPlot2_E, Name("data_ctauBkg"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(0.7));

  // --- find y-max ---
  myPlot2_E->GetYaxis()->SetRangeUser(10e-2, 10e7);
  TH1 *h = ws->data("dataToFit")->createHistogram("hist", *ws->var("ctau3D"), Binning(myPlot_E->GetNbinsX(), myPlot_E->GetXaxis()->GetXmin(), myPlot_E->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  // Yup = YMax*TMath::Power((YMax/0.1), 0.5);
  Yup = YMax * TMath::Power((YMax / 0.01), 0.5);
  Ydown = 0.01;

  // --- other cosmetics ---
  myPlot2_E->GetYaxis()->SetRangeUser(Ydown, Yup);
  myPlot2_E->GetXaxis()->SetRangeUser(-4, 7);
  myPlot2_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi)} (mm)");
  myPlot2_E->SetFillStyle(4000);
  myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_E->GetXaxis()->SetLabelSize(0);
  myPlot2_E->GetXaxis()->SetTitleSize(0);

  // --- draw vertical lines to show ctauMin and max ---
  TLine *minline = new TLine(ctauMin, 0.0, ctauMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  myPlot2_E->addObject(minline);
  TLine *maxline = new TLine(ctauMax, 0.0, ctauMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  myPlot2_E->addObject(maxline);
  myPlot2_E->Draw();

  // --- draw legend ---
  TLegend *leg_E = new TLegend(text_x + 0.033 + 0.25, text_y + 0.09 + 0.04, text_x + 0.033 + 0.38, text_y + 0.09 - 0.15);
  leg_E->SetTextSize(text_size);
  leg_E->SetTextFont(43);
  leg_E->SetBorderSize(0);
  leg_E->AddEntry(myPlot2_E->findObject("data_ctauBkg"), "Data_Bkg", "pe");
  leg_E->AddEntry(myPlot2_E->findObject("ctauBkg_Tot"), "Total PDF", "fl");
  // leg_E->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
  leg_E->Draw("same");

  // --- print latex ---
  int drawCounterLeft = 0;
  int drawCounterRight = 0;

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.033, text_y + 0.09 - y_diff * drawCounterLeft++, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.033, text_y + 0.09 - y_diff * drawCounterLeft++, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.033, text_y + 0.09 - y_diff * drawCounterLeft++, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.09 - y_diff * drawCounterLeft++, text_color, text_size);

  drawText(Form("N_{Bkg} = %.f #pm %.f", ws->var("nBkgNp")->getVal(), ws->var("nBkgNp")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("b_{Bkg} = %.4f #pm %.4f", ws->var("b_Bkg")->getVal(), ws->var("b_Bkg")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("#lambda_{Core} = %.4f #pm %.4f", ws->var("lambdaC")->getVal(), ws->var("lambdaC")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("#lambda_{N} = %.4f #pm %.4f", ws->var("lambdaN")->getVal(), ws->var("lambdaN")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("#lambda_{P} = %.4f #pm %.4f", ws->var("lambdaP")->getVal(), ws->var("lambdaP")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("f_{Core} = %.4f #pm %.4f", ws->var("q1")->getVal(), ws->var("q1")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);
  drawText(Form("f^{weight}_{N} = %.4f #pm %.4f", ws->var("q2")->getVal(), ws->var("q2")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * drawCounterRight++, text_color, text_size - 0.2);


  // --- pull pad ---
  TPad *pad_E_2 = new TPad("pad_E_2", "pad_E_2", 0, 0.006, 0.98, 0.227);
  c_E->cd();
  pad_E_2->Draw();
  pad_E_2->cd();
  pad_E_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_E_2->SetBottomMargin(0.67);
  pad_E_2->SetBottomMargin(0.4);
  pad_E_2->SetFillStyle(4000);
  pad_E_2->SetFrameFillStyle(4000);
  pad_E_2->SetTicks(1, 1);

  RooPlot *frameTMP = (RooPlot *)myPlot2_E->Clone("TMP");
  RooHist *hpull_E = frameTMP->pullHist("data_ctauBkg", "ctauBkg_Tot", true);
  hpull_E->SetMarkerSize(0.8);
  RooPlot *pullFrame_E = ws->var("ctau3D")->frame(Title("Pull Distribution"), Bins(nCtauBins), Range(ctauLow, ctauHigh));
  pullFrame_E->addPlotable(hpull_E, "PX");
  pullFrame_E->SetTitle("");
  pullFrame_E->SetTitleSize(0);
  pullFrame_E->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_E->GetYaxis()->SetTitle("Pull");
  pullFrame_E->GetYaxis()->SetTitleSize(0.15);
  pullFrame_E->GetYaxis()->SetLabelSize(0.15);
  pullFrame_E->GetYaxis()->SetRangeUser(-7, 7);
  pullFrame_E->GetYaxis()->CenterTitle();
  pullFrame_E->GetYaxis()->SetTickSize(0.04);
  pullFrame_E->GetYaxis()->SetNdivisions(404);

  pullFrame_E->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi)} (mm)");
  // pullFrame_E->GetXaxis()->SetRangeUser(-1, 7);
  pullFrame_E->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_E->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_E->GetXaxis()->SetLabelSize(0.15);
  pullFrame_E->GetXaxis()->SetTitleSize(0.15);
  pullFrame_E->GetXaxis()->CenterTitle();
  pullFrame_E->GetXaxis()->SetTickSize(0.03);
  pullFrame_E->Draw();

  TLine *lD = new TLine(ctauLow, 0, ctauHigh, 0);
  lD->SetLineStyle(1);
  lD->Draw("same");

  printChi2(ws, pad_E_2, frameTMP, fitResult, "ctau3D", "data_ctauBkg", "ctauBkg_Tot", nCtauBins, false);

  pad_E_2->Update();

  c_E->Update();
  c_E->SaveAs(Form("figs/2DFit_%s/CtauBkg/Bkg_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP));
}


void CtauBkgFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauBkg/CtauBkgResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), fname.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");
  
  fitResult->Write();
  ws->pdf("pdfCTAU_Bkg_Tot")->Write();

  outFile->Close();

  fitResult->Print(); // print fit result
  fitResult->correlationMatrix().Print();
}

// =================================
// ===== legacy codes or notes =====
// =================================

// #include <iostream>
// #include "../headers/cutsAndBin.h"
// #include "../headers/CMS_lumi_v2mass.C"
// #include "../headers/tdrstyle.C"
// #include "../headers/rootFitHeaders.h"
// #include "../headers/commonUtility.h"
// #include "../headers/JpsiUtility.h"

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