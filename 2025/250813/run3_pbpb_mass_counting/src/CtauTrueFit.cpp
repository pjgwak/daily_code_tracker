#include "CtauTrueFit.h"
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
#include "/work/pjgwak/pol24/headers/polarizationUtilities.h"

using std::cout; using std::endl;
using std::min;
using namespace RooFit;

CtauTrueFit::CtauTrueFit(float ptLow, float ptHigh,
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

CtauTrueFit::~CtauTrueFit()
{
  delete fInputTrue;
  delete fitResult;
  delete ws;
}

void CtauTrueFit::init()
{
  cout << "--- init() ---\n\n";
  setLabels();
  makeOutputFolder();
  turnOffRooFitMessage();
  openInput();
  processDataset();
  setVariableRanges();
  fixParameters();
  buildCtauTrueModel();
}

void CtauTrueFit::run()
{
  cout << "--- run() ---\n\n";
  doFit();
  rootlogon(isLogon); // cms cosmetics
  drawPlot();
  saveOutput();
}

void CtauTrueFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  // bCont means fname in other classes
  if (PR == 0) bCont = "Prompt";
  else if (PR == 1) bCont = "NonPrompt";
  else if (PR == 2) bCont = "Inclusive"; 
}

void CtauTrueFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauTrue", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauTrue", DATE.Data()), kTRUE);
}

void CtauTrueFit::turnOffRooFitMessage()
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
}

void CtauTrueFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fInputTrue = new TFile("/work/pjgwak/pol24/input_roodataset/roots/OniaRooDataSet_miniAOD_isMC1_Jpsi_cent0_180_Effw0_Accw0_PtW0_TnP0_Run2_pp_NP_GenOnly.root", "read");
}

void CtauTrueFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // --- kinematic cut ---
  TString kineCutMC = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>2.6 && mass<3.5", ptLow, ptHigh, yLow, yHigh);

  // acc cut is not used
  // TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut

  // kineCutMC = accCut+kineCutMC;
  kineCutMC = kineCutMC;

  // --- load dataset ---
  RooDataSet *dataset = (RooDataSet *)fInputTrue->Get("dataset");

  // --- declare workspace ---
  ws = new RooWorkspace("workspace");
  ws->import(*dataset);

  // --- reduce and load dataset ---
  // GenOnly no-weight
  RooDataSet *datasetMC = new RooDataSet("datasetMC", "A sample", RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("pt1")), *(ws->var("pt2")), *(ws->var("eta1")), *(ws->var("eta2"))), Import(*dataset));
  ws->import(*datasetMC);

  RooDataSet *reducedDS_MC = (RooDataSet *)datasetMC->reduce(RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCutMC.Data());
  ws->import(*reducedDS_MC, Rename("reducedDS_MC"));

  // reducedDS_MC->Print();

  delete datasetMC;
  delete reducedDS_MC;
}

void CtauTrueFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  ws->var("ctau3D")->setRange(ctau3DMin, ctau3DMax);
  ws->var("ctau3D")->setRange("ctauTrueRange", ctau3DMin, ctau3DMax);
  ws->var("ctau3D")->setRange("ctauTruePlotRange", ctau3DMin, ctau3DHigh);
  // ws->var("ctau3D")->Print();
}

void CtauTrueFit::initVar(const std::string &varName, double init, double low, double high)
{
  RooRealVar *var = ws->var(varName.c_str());
  if (!var)
  {
    std::cerr << "[ERROR] there is no variable:: " << varName << "\n";
    exit(1);
  }

  if (init < low || init > high)
  {
    std::cerr << "[ERROR] init value out of bounds for variable: " << varName << "\n";
    std::cerr << "        init = " << init << ", range = [" << low << ", " << high << "]\n";
    exit(1);
  }

  var->setVal(init);
  // var->setMin(low);
  // var->setMax(high);
  var->setRange(low, high);
}

void CtauTrueFit::fixParameters()
{
  // cout << "--- fixParameters() ---\n\n";
}

void CtauTrueFit::buildCtauTrueModel()
{
  cout << "--- buildCtauTrueModel() ---\n\n";

  // --- declare parameters ---
  ws->factory(Form("N_Jpsi_MC[300000, 20000, 1000000]"));
  ws->factory("lambdaDSS[0.1, 0.001, 0.5]");
  ws->factory("r_lambda2[1.0, 1, 10.0]");
  ws->factory("RooFormulaVar::lambdaDSS2('@0*@1', {lambdaDSS, r_lambda2})");
  ws->factory("fDSS[0.8, 0.01, 1.]");
  ws->factory("lambdaDSS3[0.5, 0.01, 2.0]");
  ws->factory("fDSS1[0.9, 0.1, 1.]");

  // --- build models ---
  ws->factory(Form("TruthModel::%s(%s)", "pdfCTAUTRUERES", "ctau3D")); // GenOnly res

  // decay model with resolution
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUEDSS1",
                   "ctau3D",
                   "lambdaDSS",
                   "pdfCTAUTRUERES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUEDSS2",
                   "ctau3D",
                   "lambdaDSS2",
                   "pdfCTAUTRUERES"));
  ws->factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", "pdfCTAUTRUEDSS3",
                   "ctau3D",
                   "lambdaDSS3",
                   "pdfCTAUTRUERES"));

  // combine decay models - 2 or 3 cases
  RooAddPdf *pdfCTAUTRUE1 = nullptr;
  RooAddPdf *pdfCTAUTRUE = nullptr;

  if (nExp == 2) ws->factory("AddPdf::pdfCTAUTRUE({pdfCTAUTRUEDSS1,pdfCTAUTRUEDSS2},{fDSS})");
  else if (nExp == 3) {
    ws->factory("AddPdf::pdfCTAUTRUE1({pdfCTAUTRUEDSS1,pdfCTAUTRUEDSS2},{fDSS})");
    ws->factory("AddPdf::pdfCTAUTRUE({pdfCTAUTRUE1,pdfCTAUTRUEDSS3},{fDSS1})");
  }

  // --- the model for the fit ---
  ws->factory(Form("RooExtendPdf::%s(%s,%s)", "TrueModel_Tot", "pdfCTAUTRUE", "N_Jpsi_MC"));
}

void CtauTrueFit::doFit()
{
  cout << "--- doFit() ---\n\n";
  // --- apply user custom fit range cut ---
  // Todo: apply if condition
  RooDataSet *dataToFit = (RooDataSet *)ws->data("reducedDS_MC")->reduce(Form("ctau3D>=%.f&&ctau3D<%.f", ctau3DMin, ctau3DMax));
  ws->import(*dataToFit, Rename("dataToFit"));

  bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("TrueModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(nCPU), Range("ctauTrueRange"), PrintLevel(-1), PrintEvalErrors(-1), Timer(true));

  delete dataToFit;
}

void CtauTrueFit::drawPlot()
{
  cout << "--- drawPlot() ---\n\n";
  // set plot normRange
  ws->pdf("TrueModel_Tot")->setNormRange("CtauTrueRange");

  // --- canvas and pad ---
  TCanvas *c_D = new TCanvas("canvas_D", "My plots", 800, 800);
  c_D->cd();

  TPad *pad_D_1 = new TPad("pad_D_1", "pad_D_1", 0.0, 0.3, 1.0, 1.0);
  pad_D_1->SetBottomMargin(0);
  pad_D_1->SetTicks(1, 1);
  pad_D_1->Draw();
  pad_D_1->cd();
  gPad->SetLogy();

  // --- ctau frames ---  
  RooPlot *myPlot_D = ws->var("ctau3D")->frame(Bins(nBins), Range("ctauTruePlotRange"));
  myPlot_D->SetTitle("");

  RooPlot *myPlot2_D = (RooPlot *)myPlot_D->Clone(); // clone to avoid memory error

  // --- plotOn ---
  ws->data("dataToFit")->plotOn(myPlot2_D, Name("MCHist_Tot"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7));

  ws->pdf("TrueModel_Tot")->plotOn(myPlot2_D, Name("MCpdf_Tot"), LineColor(kRed + 2), NormRange("ctauTrueRange"), Range("ctauTruePlotRange"));
  // ws->pdf("TrueModel_Tot")->plotOn(myPlot2_D, Name("MCpdf_Tot"), Normalization(1, RooAbsReal::RelativeExpected), Precision(1e-4), LineColor(kRed + 2), NormRange("ctauTrueRange"), Range("ctauTrueRange"));

  // --- find y-max ---
  TH1 *h = ws->data("dataToFit")->createHistogram("hist", *ws->var("ctau3D"), Binning(myPlot2_D->GetNbinsX(), myPlot2_D->GetXaxis()->GetXmin(), myPlot2_D->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));

  Double_t Yup(0.), Ydown(0.);
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / 0.6)));
  Yup = YMax * TMath::Power((YMax / YMin), (0.3 / 0.6));

  // --- cosmetics ---
  myPlot2_D->GetYaxis()->SetRangeUser(Ydown, Yup);
  // myPlot2_D->GetYaxis()->SetRangeUser(0, 150000);
  // myPlot2_D->GetXaxis()->SetRangeUser(-1, 7);
  myPlot2_D->GetXaxis()->CenterTitle();
  myPlot2_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  myPlot2_D->SetFillStyle(4000);
  myPlot2_D->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_D->GetXaxis()->SetLabelSize(0);
  myPlot2_D->GetXaxis()->SetTitleSize(0);
  myPlot2_D->Draw();

  // --- legend ---
  TLegend *leg_D = new TLegend(text_x + 0.033 + 0.25, text_y + 0.09 + 0.04, text_x + 0.033 + 0.4, text_y + 0.09 - 0.1);
  leg_D->SetTextSize(text_size);
  leg_D->SetTextFont(43);
  leg_D->SetBorderSize(0);
  leg_D->AddEntry(myPlot2_D->findObject("MCHist_Tot"), "Data", "pe");
  // leg_D->AddEntry(myPlot2_D->findObject("MCHist_Tot_NoW"),"Data w/o WF","pe");
  leg_D->AddEntry(myPlot2_D->findObject("MCpdf_Tot"), "Total fit", "fl");
  // leg_D->AddEntry(myPlot2_E->findObject("test"),"?? PDF","l");
  leg_D->Draw("same");

  // --- latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.033, text_y + 0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.09 - y_diff * 2, text_color, text_size);

  int count_yRight = 0;
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi_MC")->getVal(), ws->var("N_Jpsi_MC")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * count_yRight++, text_color, text_size);
  drawText(Form("#lambda_{1} = %.4f #pm %.4f", ws->var("lambdaDSS")->getVal(), ws->var("lambdaDSS")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * count_yRight++, text_color, text_size);
  drawText(Form("#lambda_{2}/#lambda_{1} = %.4f #pm %.4f", ws->var("r_lambda2")->getVal(), ws->var("r_lambda2")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * count_yRight++, text_color, text_size);
  if (nExp >= 3)
    drawText(Form("#lambdaDSS3 = %.4f #pm %.4f", ws->var("lambdaDSS3")->getVal(), ws->var("lambdaDSS3")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * count_yRight++, text_color, text_size);

  // --- pull pad ---
  TPad *pad_D_2 = new TPad("pad_D_2", "pad_D_2", 0.0, 0.0, 1.0, 0.3);
  c_D->cd();
  pad_D_2->Draw();
  pad_D_2->cd();
  pad_D_2->SetTopMargin(0.001); // Upper and lower plot are joined
  pad_D_2->SetBottomMargin(0.67);
  pad_D_2->SetBottomMargin(0.4);
  pad_D_2->SetFillStyle(4000);
  pad_D_2->SetFrameFillStyle(4000);
  pad_D_2->SetTicks(1, 1);

  // --- pull frames ---
  RooPlot *frameTMP = (RooPlot *)myPlot2_D->Clone("TMP");
  RooHist *hpull_D = frameTMP->pullHist(0, 0, true);
  hpull_D->SetMarkerSize(0.8);
  RooPlot *pullFrame_D = ws->var("ctau3D")->frame(Title("Pull Distribution"), Bins(myPlot2_D->GetNbinsX()), Range(myPlot2_D->GetXaxis()->GetXmin(), myPlot2_D->GetXaxis()->GetXmax()));
  pullFrame_D->addPlotable(hpull_D, "PX");
  
  // --- pull cosmetics ---
  pullFrame_D->SetTitle("");
  pullFrame_D->SetTitleSize(0);
  pullFrame_D->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_D->GetYaxis()->SetTitle("Pull");
  pullFrame_D->GetYaxis()->SetTitleSize(0.15);
  pullFrame_D->GetYaxis()->SetLabelSize(0.15);
  pullFrame_D->GetYaxis()->SetRangeUser(-7, 7);
  pullFrame_D->GetYaxis()->CenterTitle();
  pullFrame_D->GetYaxis()->SetTickSize(0.04);
  pullFrame_D->GetYaxis()->SetNdivisions(404);

  pullFrame_D->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  pullFrame_D->GetXaxis()->SetRangeUser(ctau3DMin, ctau3DMax);
  pullFrame_D->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_D->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_D->GetXaxis()->SetLabelSize(0.15);
  pullFrame_D->GetXaxis()->SetTitleSize(0.15);
  pullFrame_D->GetXaxis()->CenterTitle();
  pullFrame_D->GetXaxis()->SetTickSize(0.03);
  pullFrame_D->Draw();

  printChi2(ws, pad_D_2, frameTMP, fitResult, "ctau3D", "MCHist_Tot", "MCpdf_Tot", nBins, false);

  // draw lines at zero
  TLine *lD = new TLine(ctau3DMin, 0, ctau3DHigh, 0);
  lD->SetLineStyle(1);
  lD->Draw("same");
  pad_D_2->Update();

  // draw and save a plot
  c_D->Update();
  c_D->SaveAs(Form("figs/2DFit_%s/CtauTrue/ctauTrue_%s_%s.pdf", DATE.Data(), bCont.Data(), kineLabel.Data()));

  delete h;
  delete c_D;
}

void CtauTrueFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauTrue/CtauTrueResult_%s_%s.root", DATE.Data(), bCont.Data(), kineLabel.Data()), "RECREATE");
  
  ws->pdf("TrueModel_Tot")->Write();
  fitResult->Write();

  outFile->Close();

  fitResult->Print("v");

  delete outFile;
}

// =================================
// ===== legacy codes or notes =====
// =================================