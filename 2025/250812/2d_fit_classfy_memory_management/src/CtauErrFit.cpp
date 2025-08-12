#include "CtauErrFit.h"
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
#include "RooStats/SPlot.h"
#include "/work/pjgwak/pol24/headers/polarizationUtilities.h"

using std::cout; using std::endl;
using std::min;
using namespace RooFit;

CtauErrFit::CtauErrFit(float ptLow, float ptHigh,
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

CtauErrFit::~CtauErrFit()
{
  delete fInputData;
  delete fMass;
  delete ws;
}

void CtauErrFit::init()
{
  cout << "--- init() ---\n\n";
  setLabels();
  makeOutputFolder();
  turnOffRooFitMessage();
  openInput();
  processDataset();
  setVariableRanges();
}

void CtauErrFit::run()
{
  cout << "--- run() ---\n\n";
  doSplotFit();
  buildOutputs();
  rootlogon(isLogon);
  drawPlot();
  saveOutput();
}

void CtauErrFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  if (PRw == 1) bCont = "PR";
  else if (PRw == 2) bCont = "NP";
}

void CtauErrFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauErr", DATE.Data()));
  gSystem->mkdir(Form("figs/2DFit_%s/CtauErr", DATE.Data()));
}

void CtauErrFit::turnOffRooFitMessage()
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

void CtauErrFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fInputData = new TFile(Form("/disk1/Oniatree/miniAOD/Run2025OO/OniaRooDataSet_miniAOD_2025OORun_isMC0_Charmonia_Effw0_Accw0_PtW0_TnP0_250722.root"));
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));
}

void CtauErrFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // --- kine cuts ---
  // kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f&& cBin>%d && cBin<%d",ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh); // with cBin
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh); // w/o cBin

  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";
  kineCut = OS + accCut + kineCut;

  // --- load objects ---
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  RooDataSet *datasetMass = (RooDataSet *)fMass->Get("datasetMass");
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");

  // --- declare workspace ---
  ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  ws->import(*datasetMass);
  ws->import(*pdfMASS_Tot);

  RooArgSet *argSet = new RooArgSet(*(ws->var("ctau3D")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight")), *(ws->var("ctau3DRes")), *(ws->var("ctau3DErr")));
  argSet->add(*(ws->var("pt1")));
  argSet->add(*(ws->var("pt2")));
  argSet->add(*(ws->var("eta1")));
  argSet->add(*(ws->var("eta2")));
  argSet->add(*(ws->var("recoQQsign")));
  // argSet->add(*(ws->var("cBin")));

  // --- dataset ---
  RooDataSet *datasetW = new RooDataSet("datasetW", "A sample", *argSet, Import(*dataset)); // basically sPlot use not-weighted dataset
  // WeightVar(*ws->var("weight")));

  ws->import(*datasetW);
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(*argSet, kineCut.Data());
  ws->import(*dsAB, Rename("dsAB"));

  // dsAB->Print();
  // cout << "Weight : " << ws->var("weight")->getVal() << endl;
}

void CtauErrFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  ws->var("ctau3DErr")->setRange(ctauErrLow, ctauHigh);

  nBins = min(int(round((ws->var("ctau3DErr")->getMax() - ws->var("ctau3DErr")->getMin()) / 0.0025)), 100);
  cout << "ctau3DErr nBin : " << nBins << endl;
}

void CtauErrFit::doSplotFit()
{
  cout << "--- doSplotFit() ---\n\n";
  // --- prepare mass yileds and PDFs ---
  RooRealVar *sigYield = ws->var("N_Jpsi");
  RooRealVar *bkgYield = ws->var("N_Bkg");
  sigYield->setMin(0);
  bkgYield->setMin(0);

  RooArgList yieldList;
  yieldList.add(*ws->var("N_Jpsi"));
  yieldList.add(*ws->var("N_Bkg"));

  cout << "Sig Yield: " << sigYield->getVal() << " +/- " << sigYield->getError() << endl;
  cout << "Bkg Yield: " << bkgYield->getVal() << " +/- " << bkgYield->getError() << endl;

  RooDataSet *data = (RooDataSet *)ws->data("dsAB");
  RooArgSet *cloneSet = (RooArgSet *)RooArgSet(*ws->pdf("pdfMASS_Tot"), "pdfMASS_Tot").snapshot(kTRUE); // legacy codes
  auto clone_mass_pdf = (RooAbsPdf *)cloneSet->find("pdfMASS_Tot");
  clone_mass_pdf->setOperMode(RooAbsArg::ADirty, kTRUE);

  // --- SPlot ---
  RooStats::SPlot sData = RooStats::SPlot("sData", "An SPlot", *data, clone_mass_pdf, yieldList);
  ws->import(*data, Rename("dataset_SPLOT"));

  cout << "[INFO] Jpsi yield -> Mass Fit:" << ws->var("N_Jpsi")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("N_Jpsi") << endl;
  cout << "[INFO] Bkg  yield -> Mass Fit:" << ws->var("N_Bkg")->getVal() << ", sWeights :" << sData.GetYieldFromSWeight("N_Bkg") << endl;
}

void CtauErrFit::buildOutputs()
{
  cout << "--- buildOutputs() ---\n\n";
  // --- decdie ctau3DErr min and max ---
  TH1D *hTotTmp = (TH1D *)ws->data("dsAB")->createHistogram(("hTotTmp"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));

  for (int i = 0; i < hTotTmp->GetNbinsX() / 2; i++)
  {
    // if(hSig->GetBinContent(i)<=0&&hSig->GetBinContent(i+1)<=0&&hSig->GetBinContent(i+2)>=1&&hSig->GetBinContent(i+3)>=1){
    if (hTotTmp->GetBinContent(i) > 1)
    { // pt 7-7.5
      // if(ptLow==3&&ptHigh==4.5&&cLow==20&&cHigh==120) ctauErrMin = hSig->GetBinLowEdge(i);
      // else ctauErrMin = hSig->GetBinLowEdge(i+1);
      ctauErrMin = hTotTmp->GetBinLowEdge(i + 1);
      break;
    }
  }

  for (int i = 0; i < hTotTmp->GetNbinsX(); i++)
  {
    if (hTotTmp->GetBinContent(i) >= 1 && hTotTmp->GetBinContent(i + 1) < 1 && hTotTmp->GetBinContent(i + 2) < 1)
    {
      // ctauErrMax = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i);
      ctauErrMax = hTotTmp->GetBinLowEdge(i) + hTotTmp->GetBinWidth(i);
      break;
    }
    // else { ctauErrMax = hSig->GetBinLowEdge(i); }
    else
    {
      ctauErrMax = hTotTmp->GetBinLowEdge(i);
    }
  }

  // Todo: user custom ctauErrMin and Max
  // call function at this point (true, min, max) // -999 is default
  // change member variables: ctauErrMax, ctauErrMin

  cout << "ctauErrMax : " << ctauErrMax << " ctauErrMin : " << ctauErrMin << endl;

  // --- set new nBins with new ctau3DErr min, max ---
  double BinWidth = (ctauErrHigh - ctauErrLow) / nBins;
  newBins = (ctauErrMax - ctauErrMin) / BinWidth;
  cout << " BinWidth : " << BinWidth << endl;
  cout << " newBins : " << newBins << endl;

  // --- build sWeighted dataset with default range ---
  RooDataSet *dataw_Bkg_b = new RooDataSet("dataw_Bkg_b", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"), RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  RooDataSet *dataw_Sig_b = new RooDataSet("dataw_Sig_b", "TMP_SIG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"),RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");

  

  // --- build RooHistPdfs ---
  TH1D *hSig = (TH1D *)dataw_Sig_b->createHistogram(("hSig"), *ws->var("ctau3DErr"), Binning(nBins, ctauErrLow, ctauErrHigh));

  // total: TH1D -> RooDataHist -> RooHistPdf
  TH1D *hTot = (TH1D *)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
  RooDataHist *totHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hTot);
  RooHistPdf *pdfCTAUERR_Tot = new RooHistPdf("pdfCTAUERR_Tot", "hist pdf", *ws->var("ctau3DErr"), *totHist);

  // sig: TH1D -> RooDataHist -> RooHistPdf
  TH1D *hSig_w = (TH1D *)dataw_Sig_b->createHistogram(("hSig_w"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
  RooDataHist *sigHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hSig_w);
  RooHistPdf *pdfCTAUERR_Jpsi = new RooHistPdf("pdfCTAUERR_Jpsi", "hist pdf", *ws->var("ctau3DErr"), *sigHist);

  // bkg: TH1D -> RooDataHist -> RooHistPdf
  TH1D *hBkg_w = (TH1D *)dataw_Bkg_b->createHistogram(("hBkg_w"), *ws->var("ctau3DErr"), Binning(newBins, ctauErrMin, ctauErrMax));
  RooDataHist *bkgHist = new RooDataHist("dsAB", "", *ws->var("ctau3DErr"), hBkg_w);
  RooHistPdf *pdfCTAUERR_Bkg = new RooHistPdf("pdfCTAUERR_Bkg", "hist pdf", *ws->var("ctau3DErr"), *bkgHist);


  // --- build sWeighted dataset with ctau3DErr min and max---
  RooDataSet *dataw_Bkg = new RooDataSet("dataw_Bkg", "TMP_BKG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"), RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Bkg_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Bkg_sw");
  RooDataSet *dataw_Sig = new RooDataSet("dataw_Sig", "TMP_SIG_DATA", (RooDataSet *)ws->data("dataset_SPLOT"), RooArgSet(*ws->var("ctau3DErr"), *ws->var("N_Jpsi_sw"), *ws->var("ctau3DRes"), *ws->var("ctau3D"), *ws->var("mass")), 0, "N_Jpsi_sw");

  // import
  ws->import(*dataw_Sig);
  ws->import(*dataw_Bkg);
  ws->import(*pdfCTAUERR_Tot);
  ws->import(*pdfCTAUERR_Jpsi);
  ws->import(*pdfCTAUERR_Bkg);

  delete hTotTmp;
}

void CtauErrFit::drawPlot()
{
  cout << "--- drawPlot() ---\n\n";
    
  // TH1D* test = (TH1D*)ws->data("dsAB")->createHistogram(("test"), *ws->var("ctau3DErr"),Binning(nBins,
  // double minRange = (double)(floor(ctauErrMin * 100.) / 100.);
  // double maxRange = (double)(ceil(ctauErrMax * 100.) / 100.);
  // RooPlot *myPlot_B = ws->var("ctau3DErr")->frame(Bins(nBins), Range(minRange - 0.01, maxRange + 0.01)); // modified

  ws->var("ctau3DErr")->setRange("ctauErrWindow", ctauErrMin, ctauErrMax);

  // --- check how many events were cut off ---
  cout << ws->data("dsAB")->numEntries() << endl;
  cout << ws->data("dataset_SPLOT")->numEntries() << endl;

  // draw canvas and pads
  TCanvas *c_B = new TCanvas("canvas_B", "My plots", 800, 800);
  c_B->cd();
  TPad *pad_B_1 = new TPad("pad_B_1", "pad_B_1", 0.0, 0.3, 1.0, 1.0);
  pad_B_1->SetBottomMargin(0);
  pad_B_1->SetTicks(1, 1);
  pad_B_1->Draw();
  pad_B_1->cd();
  gPad->SetLogy();
  // pad_B_1->SetLogy();

  // --- ctauErr frame ---
  RooPlot *myPlot_B = ws->var("ctau3DErr")->frame(Bins(newBins), Range(ctauErrMin, ctauErrMax));
  myPlot_B->SetTitle("");
  RooPlot *myPlot2_B = (RooPlot *)myPlot_B->Clone(); // clone to avoid memory error

  // --- plotOn ---
  // Total: data, pdf
  ws->data("dsAB")->plotOn(myPlot2_B, Name("dataCTAUERR_Tot"), MarkerSize(.7), Binning(newBins)); // Normalization(wsmc->data("reducedDS_MC")->sumEntries()
  ws->pdf("pdfCTAUERR_Tot")->plotOn(myPlot2_B, Name("pdfCTAUERR_Tot"), LineColor(kGreen + 1), Range(ctauErrMin, ctauErrMax), LineWidth(2));

  // signal: data, pdf
  ws->data("dataw_Sig")->plotOn(myPlot2_B, Name("dataHist_Sig"), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Jpsi")->plotOn(myPlot2_B, Name("pdfCTAUERR_Jpsi"), LineColor(kRed + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));

  // bkg: data, pdf
  ws->data("dataw_Bkg")->plotOn(myPlot2_B, Name("dataHist_Bkg"), MarkerSize(.7), LineColor(kBlue + 2), MarkerColor(kBlue + 2), Binning(newBins));
  ws->pdf("pdfCTAUERR_Bkg")->plotOn(myPlot2_B, Name("pdfCTAUERR_Bkg"), LineColor(kBlue + 2), LineWidth(2), Range(ctauErrMin, ctauErrMax));

  // --- find y-max ---
  TH1D *hTot = (TH1D *)ws->data("dsAB")->createHistogram(("hTot"), *ws->var("ctau3DErr"), Binning(myPlot_B->GetNbinsX(), myPlot_B->GetXaxis()->GetXmin(), myPlot_B->GetXaxis()->GetXmax())); // tot
  Double_t YMax = hTot->GetBinContent(hTot->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= hTot->GetNbinsX(); i++)
    if (hTot->GetBinContent(i) > 0)
      YMin = min(YMin, hTot->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.4 - 0.3)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.3 / (1.0 - 0.4 - 0.3))));
  myPlot2_B->GetYaxis()->SetRangeUser(Ydown, Yup);

  // --- check how many events were cut out ---
  // cout<<ctauErrLow<<", "<<ctauErrHigh<<endl;
  cout << ws->var("ctau3DErr")->getMin() << ", " << ws->var("ctau3DErr")->getMax() << endl;
  RooDataSet *ctauResCutDS = (RooDataSet *)ws->data("dataw_Sig")->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr"))), Form("ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ctauResCutDS->SetName("ctauResCutDS");
  ws->import(*ctauResCutDS);
  Double_t outTot = ws->data("dsAB")->numEntries();
  Double_t outRes = ws->data("ctauResCutDS")->numEntries();
  cout << "Tot evt: (" << outTot << ")" << endl;
  cout << "Res evt: (" << outRes << ")" << endl;
  cout << "lost evt: (" << ((outTot - outRes) * 100) / outTot << ")%, " << outRes << "evts" << endl;

  // --- build vertical lines describing ctau3DErr min and max
  TLine *minline = new TLine(ctauErrMin, 0.0, ctauErrMin, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  minline->SetLineStyle(2);
  minline->SetLineColor(1);
  minline->SetLineWidth(3);
  myPlot2_B->addObject(minline);
  TLine *maxline = new TLine(ctauErrMax, 0.0, ctauErrMax, (Ydown * TMath::Power((Yup / Ydown), 0.5)));
  maxline->SetLineStyle(2);
  maxline->SetLineColor(1);
  maxline->SetLineWidth(3);
  myPlot2_B->addObject(maxline);

  myPlot2_B->GetXaxis()->CenterTitle();
  myPlot2_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  myPlot2_B->SetFillStyle(4000);
  myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_B->GetXaxis()->SetLabelSize(0);
  myPlot2_B->GetXaxis()->SetTitleSize(0);
  myPlot2_B->Draw();
  // Double_t outTot = ws->data("dsAB")->numEntries();
  // Double_t outErr = ws->data("dsAB")->reduce(Form("(ctauErr>=%.6f || ctauErr<=%.6f)", ctauErrHigh, ctauErrLow))->numEntries();
  // cout<<(outErr*100)/outTot<<endl;

  // --- legend ---
  TLegend *leg_B = new TLegend(text_x + 0.033 + 0.5, text_y + 0.09 - 0.2, text_x + 0.033 + 0.7, text_y + 0.09);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot2_B->findObject("dataCTAUERR_Tot"), "Data", "pe");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Tot"), "Total PDF", "l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Jpsi"), "Signal", "l");
  leg_B->AddEntry(myPlot2_B->findObject("pdfCTAUERR_Bkg"), "Background", "l");
  leg_B->Draw("same");

  // --- latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.033, text_y + 0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.09 - y_diff * 2, text_color, text_size);
  // drawText(Form("n_{J/#psi} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Jpsi"),ws->var("N_Jpsi")->getError()),text_x+0.033,text_y+0.09-y_diff*3,text_color,text_size);
  // drawText(Form("n_{Bkg} = %.f #pm %.f",sData.GetYieldFromSWeight("N_Bkg"), ws->var("N_Bkg")->getError()),text_x+0.033,text_y+0.09-y_diff*4,text_color,text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x + 0.033, text_y + 0.09 - y_diff * 3, text_color, text_size);


  // --- pull pad ---
  TPad *pad_B_2 = new TPad("pad_B_2", "pad_B_2", 0.0, 0.0, 1.0, 0.3);
  c_B->cd();
  pad_B_2->Draw();
  pad_B_2->cd();
  pad_B_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_B_2->SetBottomMargin(0.67);
  pad_B_2->SetBottomMargin(0.4);
  pad_B_2->SetFillStyle(4000);
  pad_B_2->SetFrameFillStyle(4000);
  pad_B_2->SetTicks(1, 1);

  RooPlot *frameTMP_B = (RooPlot *)myPlot2_B->Clone("TMP");
  RooHist *hpull_B = frameTMP_B->pullHist("dataCTAUERR_Tot", "pdfCTAUERR_Tot");
  hpull_B->SetMarkerSize(0.8);
  // RooPlot* pullFrame_B = ws->var("ctau3DErr")->frame(Title("Pull Distribution")) ;
  RooPlot *pullFrame_B = ws->var("ctau3DErr")->frame(Bins(newBins), Range(ctauErrMin, ctauErrMax));
  pullFrame_B->addPlotable(hpull_B, "PX");
  pullFrame_B->SetTitle("");
  pullFrame_B->SetTitleSize(0);
  pullFrame_B->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_B->GetYaxis()->SetTitle("Pull");
  pullFrame_B->GetYaxis()->SetTitleSize(0.15);
  pullFrame_B->GetYaxis()->SetLabelSize(0.15);
  pullFrame_B->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_B->GetYaxis()->CenterTitle();

  pullFrame_B->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} Error (mm)");
  pullFrame_B->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_B->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_B->GetXaxis()->SetLabelSize(0.15);
  pullFrame_B->GetXaxis()->SetTitleSize(0.15);
  pullFrame_B->GetXaxis()->CenterTitle();

  pullFrame_B->GetYaxis()->SetTickSize(0.04);
  pullFrame_B->GetYaxis()->SetNdivisions(404);
  pullFrame_B->GetXaxis()->SetTickSize(0.03);
  pullFrame_B->Draw();

  TLine *lB = new TLine(ctauErrMin, 0, ctauErrMax, 0); // pull line at zero
  lB->SetLineStyle(1);
  lB->Draw("same");
  pad_B_2->Update();

  // --- draw the plot ---
  c_B->Update();
  c_B->SaveAs(Form("figs/2DFit_%s/CtauErr/ctauErr_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));
}

void CtauErrFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->data("dataw_Bkg")->Write();
  ws->data("dataw_Sig")->Write();
  
  ws->pdf("pdfCTAUERR_Tot")->Write();
  ws->pdf("pdfCTAUERR_Bkg")->Write();
  ws->pdf("pdfCTAUERR_Jpsi")->Write();
  outFile->Close();
}

// =================================
// ===== legacy codes or notes =====
// =================================
