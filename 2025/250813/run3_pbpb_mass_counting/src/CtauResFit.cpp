#include "CtauResFit.h"
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

CtauResFit::CtauResFit(float ptLow, float ptHigh,
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

CtauResFit::~CtauResFit()
{
  delete fMass;
  delete fCErr;
  delete fitResult;
  delete ws;
}

void CtauResFit::init()
{
  cout << "--- init() ---\n\n";
  setLabels();
  makeOutputFolder();
  turnOffRooFitMessage();
  openInput();
  processDataset();
  setVariableRanges();
  buildModel();
}

void CtauResFit::run()
{
  cout << "--- run() ---\n\n";
  doFit();
  rootlogon(isLogon);
  drawPlot();
  saveOutput();
}

void CtauResFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  if (PRw == 1) bCont = "PR";
  else if (PRw == 2) bCont = "NP";
}

void CtauResFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/CtauRes", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/CtauRes", DATE.Data()), kTRUE);
}

void CtauResFit::turnOffRooFitMessage()
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

void CtauResFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fMass = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));
  fCErr = new TFile(Form("roots/2DFit_%s/CtauErr/CtauErrResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));
}

void CtauResFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // --- load dataset and pdf ---
  RooAddPdf *pdfMASS_Tot = (RooAddPdf *)fMass->Get("pdfMASS_Tot");
  RooDataSet *dataw_Sig = (RooDataSet *)fCErr->Get("dataw_Sig");

  // --- declare workspace ---
  ws = new RooWorkspace("workspace");
  ws->import(*pdfMASS_Tot);
  ws->import(*dataw_Sig);

  // --- reduce dataset ---
  double ctauErrMin;
  double ctauErrMax;
  ctauErrMin = ws->var("ctau3DErr")->getMin();
  ctauErrMax = ws->var("ctau3DErr")->getMax();

  RooDataSet *ctauResCutDS = (RooDataSet *)dataw_Sig->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr"))), Form("ctau3DRes<0&&ctau3DErr>%f&&ctau3DErr<%f", ctauErrMin, ctauErrMax));
  ws->import(*ctauResCutDS, Rename("ctauResCutDS")); // ctauResCutDS->SetName("ctauResCutDS");

  delete ctauResCutDS;
}

void CtauResFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  ws->var("mass")->setRange(massLow, massHigh);
  // ws->var("mass")->setRange("ctauRange", massLow, massHigh);
  
  ws->var("ctau3D")->setRange(ctauLow, ctauHigh);
  // ws->var("ctau3D")->setRange("ctauRange", ctauLow, ctauHigh);
  
  // ws->var("ctau3DErr")->setRange(ctauErrMin, ctauErrMax);
  // ws->var("ctau3DErr")->setRange("ctauRange", ctauErrMin, ctauErrMax);
  
  ws->var("ctau3DRes")->setRange(-10, 10);
  // ws->var("ctau3DRes")->setRange("ctauRange", -10, 10);

  ws->var("mass")->Print();
  ws->var("ctau3D")->Print();
  ws->var("ctau3DErr")->Print();
  ws->var("ctau3DRes")->Print();
}

void CtauResFit::initVar(const std::string &varName, double init, double low, double high)
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

void CtauResFit::buildModel()
{
  cout << "--- buildModel() ---\n\n";
  ws->factory("One[1.0]");
  ws->factory("ctauRes_mean[0.0]");
  ws->factory("ctau1_CtauRes[0.]");
  ws->factory("ctau2_CtauRes[0.]"); // ws->factory("s2_CtauRes[2., 1e-6, 10.]");
  ws->factory("ctau3_CtauRes[0.]"); // ws->factory("s3_CtauRes[3,  1e-6, 10.]");
  ws->factory("ctau4_CtauRes[0.]"); // ws->factory("s4_CtauRes[5.37, 0., 10.]");

  ws->factory("s1_CtauRes[0.5, 0.001, 3.0]");
  ws->factory("rS21_CtauRes[0.8, 1.0, 3.0]");
  ws->factory("rS32_CtauRes[2.2, 1.0, 5.0]");
  ws->factory("f_CtauRes[0.31, 1e-4, 1.]");
  ws->factory("f2_CtauRes[0.31, 1e-6, 1.]");

  ws->factory("RooFormulaVar::s2_CtauRes('@0*@1',{rS21_CtauRes,s1_CtauRes})");
  ws->factory("RooFormulaVar::s3_CtauRes('@0*@1',{rS32_CtauRes,s2_CtauRes})");
  // nGauss=4
  ws->factory("f3_CtauRes[0.5, 1e-6, 1.]");
  ws->factory("rS43_CtauRes[1.5, 1.0, 10.0]");
  ws->factory("RooFormulaVar::s4_CtauRes('@0*@1',{rS43_CtauRes,s3_CtauRes})");

  // if(ptLow==6.5&&ptHigh==9&&cLow==0&&cHigh==180){
  //   ws->factory("f_CtauRes[0.9, 1e-6, 1.]");ws->factory("f2_CtauRes[0.2, 1e-6, 1.]");ws->factory("f3_CtauRes[0.7, 0., 1.]");}
  // else if(ptLow==25&&ptHigh==50&&cLow==0&&cHigh==180){
  //   ws->factory("f_CtauRes[0.8, 1e-6, 1.]");ws->factory("f2_CtauRes[0.2, 1e-6, 1.]");ws->factory("f3_CtauRes[0.7, 0., 1.]");}
  // else {ws->factory("f_CtauRes[0.2, 0., 1.]");ws->factory("f2_CtauRes[0.2, 0., 1.]");ws->factory("f3_CtauRes[0.5, 0., 1.]");}
  //  create the three PDFs
  TString varName = "ctau3DRes";
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel1_ctauRes", varName.Data(),
                   "ctau1_CtauRes", //"ctau1_CtauRes",
                   "s1_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel2_ctauRes", varName.Data(),
                   "ctau2_CtauRes", //"ctau2_CtauRes",
                   "s2_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel3_ctauRes", varName.Data(),
                   "ctau3_CtauRes", //"ctau3_CtauRes",
                   "s3_CtauRes"));
  ws->factory(Form("GaussModel::%s(%s, %s, %s,One,One)", "GaussModel4_ctauRes", varName.Data(),
                   "ctau4_CtauRes", //"ctau3_CtauRes",
                   "s4_CtauRes"));
  // combine the two PDFs
  if (nGauss == 4)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel43_ctauRes", "GaussModel4_ctauRes", "GaussModel3_ctauRes", "f3_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel43_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 3)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModel23_ctauRes", "GaussModel3_ctauRes", "GaussModel2_ctauRes", "f2_CtauRes"));
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel23_ctauRes", "f_CtauRes"));
  }
  else if (nGauss == 2)
  {
    ws->factory(Form("AddModel::%s({%s, %s}, {%s})", "GaussModelCOND_ctauRes", "GaussModel1_ctauRes", "GaussModel2_ctauRes", "f_CtauRes"));
  }

  ws->factory("SUM::GaussModel_Tot(N_Jpsi*GaussModelCOND_ctauRes)");

  // RooAddPdf *GaussModel_Tot = new RooAddPdf("GaussModel_Tot");
  // RooAbsPdf *ctauResModel = ctauResModel = new RooAddPdf("GaussModel_Tot", "GaussModel_Tot", *ws->pdf("GaussModelCOND_ctauRes"), *ws->var("N_Jpsi"));
  // ws->import(*ctauResModel);
}


void CtauResFit::doFit()
{
  cout << "--- doFit() ---\n\n";
  // --- appply ctauRes range ---
  TH1D *hTot = (TH1D *)ws->data("dataw_Sig")->createHistogram(("hTot"), *ws->var("ctau3DRes"), Binning(nCtauResBins, ctauResLow, ctauResHigh));
  
  ctauResMin = -10;
  // double ctauResMax = hTot->GetBinCenter(hTot->FindLastBinAbove(1,1));
  ctauResMax = 0;
  cout << "NBins: " << hTot->GetNbinsX() << endl;
  for (int i = 0; i < hTot->GetNbinsX() / 2; i++)
  {
    // cout<<"Content: "<<hTot->GetBinContent(i)<<endl;
    if (hTot->GetBinContent(i) <= 0 && hTot->GetBinContent(i + 1) <= 0)
    {
      // cout<<"#####"<<i<<": "<<hTot->GetBinLowEdge(i+2)<<endl;
      ctauResMin = hTot->GetBinLowEdge(i + 2);
    }
    // if(hTot->GetBinContent(i)>1)ctauResMax = hTot->GetBinCenter(i)+hTot->GetBinWidth(i);
  }

  //  if (v2==-0.3&&ptLow==6.5&&ptHigh==7.5) ctauResMin=-6.4;

  ws->var("ctau3DRes")->setRange("ctauResWindow", ctauResMin, 0);
  cout << "Fit Range: " << ctauResMin << " - 0" << endl;

  RooDataSet *dataToFit = (RooDataSet *)ws->data("ctauResCutDS")->reduce(Form("ctau3DRes>=%.f&&ctau3DRes<=0", ctauResMin));
  ws->import(*dataToFit, Rename("dataToFit"));

  // --- fit ---
  bool isWeighted = ws->data("dataToFit")->isWeighted();
  fitResult = ws->pdf("GaussModel_Tot")->fitTo(*dataToFit, Save(), SumW2Error(isWeighted), Extended(kTRUE), NumCPU(nCPU), PrintLevel(-1));

  delete hTot;
  delete dataToFit;
}

void CtauResFit::drawPlot()
{
  cout << "--- drawPlot() ---\n\n";
  // --- canvas and pad
  TCanvas *c_C = new TCanvas("canvas_C", "My plots", 800, 800);
  c_C->cd();

  TPad *pad_C_1 = new TPad("pad_C_1", "pad_C_1", 0.0, 0.3, 1.0, 1.0);
  pad_C_1->SetBottomMargin(0.001);
  pad_C_1->SetTicks(1, 1);
  pad_C_1->Draw();
  pad_C_1->cd();
  gPad->SetLogy();

  // --- ctauRes frame ---
  RooPlot *myPlot_C = ws->var("ctau3DRes")->frame(Bins(nCtauResBins), Range(ctauResLow, ctauResHigh));
  myPlot_C->SetTitle("");
  
  // setFixedVarsToContantVars(ws);
  // ws->data("dataToFit")->plotOn(myPlot_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0),
  //     MarkerSize(.7), LineColor(kBlack), MarkerColor(kBlack));
  ws->data("ctauResCutDS")->plotOn(myPlot_C, Name("dataHist_ctauRes"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(.7), LineColor(kRed + 2), MarkerColor(kRed + 2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot_C, Name("modelHist_ctauRes"), Precision(1e-6), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), LineColor(kBlack));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot_C, Name("modelHist_gm1"), Precision(1e-6), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel1_ctauRes")), LineColor(kGreen + 2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot_C, Name("modelHist_gm2"), Precision(1e-6), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel2_ctauRes")), LineColor(kRed + 2));
  ws->pdf("GaussModel_Tot")->plotOn(myPlot_C, Name("modelHist_gm3"), Precision(1e-6), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel3_ctauRes")), LineColor(kBlue + 2));
  if (nGauss == 4)
  {
    ws->pdf("GaussModel_Tot")->plotOn(myPlot_C, Name("modelHist_gm4"), Precision(1e-6), NormRange("ctauResWindow"), Normalization(ws->data("ctauResCutDS")->sumEntries(), RooAbsReal::NumEvent), Components(*ws->pdf("GaussModel4_ctauRes")), LineColor(kMagenta + 2));
  }

  // --- find y-max ---
  TH1 *h = ws->data("ctauResCutDS")->createHistogram("hist", *ws->var("ctau3DRes"), Binning(myPlot_C->GetNbinsX(), myPlot_C->GetXaxis()->GetXmin(), myPlot_C->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.), Ydown(0.);
  Yup = YMax * TMath::Power((YMax / YMin), (0.5 / (1.0 - 0.5 - 0.2)));
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.2 / (1.0 - 0.5 - 0.2))));
  myPlot_C->GetYaxis()->SetRangeUser(Ydown, Yup);
  
  // --- check here ---
  Double_t outTot = ws->data("ctauResCutDS")->numEntries();
  Double_t outRes = ws->data("dataToFit")->numEntries();
  cout << "Tot evt: (" << outTot << ")" << endl;
  cout << "Res evt: (" << outRes << ")" << endl;
  cout << "lost evt: (" << (outRes * 100) / outTot << ")%, " << outRes << "evts" << endl;

  // --- draw vertical lines to display the ranges ---
  if (outRes > 0.0)
  {
    // TLine   *minline = new TLine(rangeErr[0], 0.0, rangeErr[0], Ydown*TMath::Power((Yup/Ydown),0.4));
    TLine *minline = new TLine(ctauResMin, 0.0, ctauResMin, Ydown * TMath::Power((Yup / Ydown), 0.4));
    minline->SetLineStyle(2);
    minline->SetLineColor(1);
    minline->SetLineWidth(3);
    myPlot_C->addObject(minline);
    TLine *maxline = new TLine(ctauResMax, 0.0, ctauResMax, Ydown * TMath::Power((Yup / Ydown), 0.4));
    maxline->SetLineStyle(2);
    maxline->SetLineColor(1);
    maxline->SetLineWidth(3);
    myPlot_C->addObject(maxline);
  }
  myPlot_C->GetXaxis()->CenterTitle();
  myPlot_C->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  myPlot_C->SetFillStyle(4000);
  myPlot_C->GetYaxis()->SetTitleOffset(1.43);
  myPlot_C->GetXaxis()->SetLabelSize(0);
  myPlot_C->GetXaxis()->SetTitleSize(0);
  myPlot_C->Draw();

  cout << "ctauRes range: [" << ctauResMin << ", " << ctauResMax << "]";

  // --- legend ---
  TLegend *leg_C = new TLegend(text_x + 0.033 + 0.29, text_y + 0.09 + 0.03, text_x + 0.033 + 0.39, text_y + 0.09 - 0.17);
  leg_C->SetTextSize(text_size);
  leg_C->SetTextFont(43);
  leg_C->SetBorderSize(0);
  leg_C->AddEntry(myPlot_C->findObject("dataHist_ctauRes"), "Data", "pe");
  leg_C->AddEntry(myPlot_C->findObject("modelHist_ctauRes"), "Total PDF", "l");
  leg_C->AddEntry(myPlot_C->findObject("modelHist_gm1"), "Gauss 1", "l");
  leg_C->AddEntry(myPlot_C->findObject("modelHist_gm2"), "Gauss 2", "l");
  if (nGauss == 3)
    leg_C->AddEntry(myPlot_C->findObject("modelHist_gm3"), "Gauss 3", "l");
  if (nGauss == 4)
    leg_C->AddEntry(myPlot_C->findObject("modelHist_gm4"), "Gauss 4", "l");
  leg_C->Draw("same");


  // --- latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x + 0.033, text_y + 0.09, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x + 0.033, text_y + 0.09 - y_diff, text_color, text_size);
  drawText(Form("Cent. %d - %d%s", cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.09 - y_diff * 2, text_color, text_size);
  drawText(Form("Loss: (%.4f%s) %.f evts", (outTot - outRes) * 100 / outTot, "%", outTot - outRes), text_x + 0.033, text_y + 0.09 - y_diff * 3, text_color, text_size);
  // cout<<"lost evt: ("<<(outRes*100)/outTot<<")%, "<<outRes<<"evts"<<endl;

  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x + 0.033 + 0.5, text_y + 0.09, text_color, text_size);
  drawText(Form("s1_{Res} = %.4f #pm %.4f", ws->var("s1_CtauRes")->getVal(), ws->var("s1_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 1, text_color, text_size);
  drawText(Form("(s2/s1)_{Res} = %.4f #pm %.4f", ws->var("rS21_CtauRes")->getVal(), ws->var("rS21_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 2, text_color, text_size);
  if (nGauss == 3)
  {
    drawText(Form("(s3/s2)_{Res} = %.4f #pm %.4f", ws->var("rS32_CtauRes")->getVal(), ws->var("rS32_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 3, text_color, text_size);
    drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 4, text_color, text_size);
    drawText(Form("f2_{Res} = %.4f #pm %.4f", ws->var("f2_CtauRes")->getVal(), ws->var("f2_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 5, text_color, text_size);
  }
  else if (nGauss == 2)
  {
    drawText(Form("f_{Res} = %.4f #pm %.4f", ws->var("f_CtauRes")->getVal(), ws->var("f_CtauRes")->getError()), text_x + 0.033 + 0.5, text_y + 0.09 - y_diff * 3, text_color, text_size);
  }

  // --- pull pad ---
  TPad *pad_C_2 = new TPad("pad_C_2", "pad_C_2", 0.0, 0.0, 1.0, 0.3);
  c_C->cd();
  pad_C_2->Draw();
  pad_C_2->cd();
  pad_C_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_C_2->SetBottomMargin(0.67);
  pad_C_2->SetBottomMargin(0.4);
  pad_C_2->SetFillStyle(4000);
  pad_C_2->SetFrameFillStyle(4000);
  pad_C_2->SetTicks(1, 1);

  RooHist *hpull_C = myPlot_C->pullHist("dataHist_ctauRes", "modelHist_ctauRes", true);
  hpull_C->SetMarkerSize(0.8);
  // RooPlot* pullFrame_C = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(nCtauResBins), Range(ctauResLow, ctauResHigh)) ;
  RooPlot *pullFrame_C = ws->var("ctau3DRes")->frame(Title("Pull Distribution"), Bins(myPlot_C->GetNbinsX()), Range(myPlot_C->GetXaxis()->GetXmin(), myPlot_C->GetXaxis()->GetXmax()));
  pullFrame_C->addPlotable(hpull_C, "PX");
  pullFrame_C->SetTitle("");
  pullFrame_C->SetTitleSize(0);
  pullFrame_C->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_C->GetYaxis()->SetTitle("Pull");
  pullFrame_C->GetYaxis()->SetTitleSize(0.15);
  pullFrame_C->GetYaxis()->SetLabelSize(0.15);
  pullFrame_C->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_C->GetYaxis()->CenterTitle();

  pullFrame_C->GetXaxis()->SetTitle("#frac{#font[12]{l}_{J/#psi}}{#sigma_{#font[12]{l}_{J/#psi}}}");
  pullFrame_C->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_C->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_C->GetXaxis()->SetLabelSize(0.15);
  pullFrame_C->GetXaxis()->SetTitleSize(0.15);
  pullFrame_C->GetXaxis()->CenterTitle();

  pullFrame_C->GetYaxis()->SetTickSize(0.04);
  pullFrame_C->GetYaxis()->SetNdivisions(404);
  pullFrame_C->GetXaxis()->SetTickSize(0.03);
  pullFrame_C->Draw();

  TLine *lC = new TLine(ctauResLow, 0, ctauResHigh, 0);
  lC->SetLineStyle(1);
  lC->Draw("same");

  printChi2(ws, pad_C_2, myPlot_C, fitResult, "ctau3DRes", "dataHist_ctauRes", "modelHist_ctauRes", nCtauResBins, false);
  pad_C_2->Update();

  c_C->Update();
  c_C->SaveAs(Form("figs/2DFit_%s/CtauRes/CtauRes_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));

  delete h;
  delete c_C;
}

void CtauResFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";

  // RooArgSet *fitargs = new RooArgSet();
  // fitargs->add(fitResult->floatParsFinal());
  // RooDataSet *datasetRes = new RooDataSet("datasetRes", "dataset with Resolution Fit result", *fitargs);

  // --- declare output file ---
  TFile *outFile = new TFile(Form("roots/2DFit_%s/CtauRes/CtauResResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");
  
  ws->pdf("GaussModel_Tot")->Write("GaussModel_Tot");
  // GaussModel_Tot->Write();
  //	ctauResCutDS->Write();
  // datasetRes->Write();
  fitResult->Write();
  
  outFile->Close();

  // print result
  fitResult->Print("V");

  delete outFile;
}

// =================================
// ===== legacy codes or notes =====
// =================================
