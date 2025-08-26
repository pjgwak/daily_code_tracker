#include "McMassFit.h"
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
#include "/work/pjgwak/pol24/headers/polarizationUtilities.h"

using std::cout; using std::endl;
using std::min;
using namespace RooFit;

McMassFit::McMassFit(float ptLow, float ptHigh,
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

McMassFit::~McMassFit()
{
  delete fInputMc;
  delete fitResult;
  delete ws;
}

void McMassFit::init()
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

void McMassFit::run()
{
  cout << "--- run() ---\n\n";
  doFit();
  rootlogon(isLogon);
  drawPlot();
  saveOutput();
}

void McMassFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";

  bCont = "PR"; // always use PR MC
  // if (PR == 0) bCont = "PR";
  // else if (PR == 1) bCont = "NP";
  // else if (PR == 2) bCont = "Inclusive";
}

void McMassFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  gSystem->mkdir(Form("roots/2DFit_%s/mc_Mass", DATE.Data()), true);
  gSystem->mkdir(Form("figs/2DFit_%s/mc_Mass", DATE.Data()), true);
}

void McMassFit::turnOffRooFitMessage()
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

  // RooMsgService::instance().getStream(0).removeTopic(Caching);
  // RooMsgService::instance().getStream(1).removeTopic(Caching);
  // RooMsgService::instance().getStream(0).removeTopic(Plotting);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Integration);
  // RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  // RooMsgService::instance().getStream(1).removeTopic(Fitting);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Eval);
  // RooMsgService::instance().getStream(1).removeTopic(DataHandling);
  // // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  // RooMsgService::instance().setGlobalKillBelow(ERROR);
  // RooMsgService::instance().setSilentMode(true);
}

void McMassFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fInputMc = new TFile("/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root");
}

void McMassFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // --- kinematic cuts ---
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f && cBin>=%d && cBin<%d", ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";

  TString angle_cut = Form("&& cos_ep>%.2f && cos_ep<%.2f", cosLow, cosHigh);

  kineCut = OS + accCut + kineCut + angle_cut;

  // declare workspace
  ws = new RooWorkspace("workspace");

  // --- load inputs ---
  RooDataSet *dataset = (RooDataSet *)fInputMc->Get("dataset");
  ws->import(*dataset);

  // isWeight?
  RooDataSet *datasetW = nullptr;
  if (isWeighted)
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset), WeightVar(*ws->var("weight")));
  else
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset));

  // apply cuts
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data());
  ws->import(*dsAB, Rename("dsAB"));

  delete dsAB;
  delete datasetW;
}

void McMassFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setRange("mcMassPlot", massLow, massHigh);
  ws->var("mass")->Print();
}

void McMassFit::initVar(const std::string &varName, double init, double low, double high)
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

void McMassFit::buildModel()
{
  cout << "--- buildModel() ---\n\n";
  // --- signal parameters ---
  // The order is {sigma_1,  x, alpha_1, n_1,   f, m_lambda}
  double paramsupper[8] = {1.1, 1.1, 3.1, 3.1, 1., 25.0};
  double paramslower[8] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, -25.0};

  double sigma_1_init = 0.5;
  double x_init = 0.43;
  double alpha_1_init = 0.5;
  double n_1_init = 0.8;
  double f_init = 0.4;
  double sl1_mean = 0.35, sl2_mean = 0.004, sl3_mean = 0.006;
  double N_Jpsi_high = 100000; // 2500000
  double N_Bkg_high = 200000;
  
  double m_lambda_init = 5;


  ws->factory("mean[3.096, 3.086, 3.106]");
  ws->factory("x_A[1.1,1,3]");
  ws->factory("sigma_1_A[0.01,0.001,0.1]");
  ws->factory("alpha_1_A[1.5,0.8,5]");
  ws->factory("n_1_A[1.5,0.8,5]");
  ws->factory("f[0.6,0.05,0.95]");

  RooFormulaVar sigma_2_A("sigma_2_A", "@0*@1",
                          RooArgList(*ws->var("sigma_1_A"), *ws->var("x_A")));
  RooFormulaVar alpha_2_A("alpha_2_A", "1.0*@0",
                          RooArgList(*ws->var("alpha_1_A")));
  RooFormulaVar n_2_A("n_2_A", "1.0*@0",
                      RooArgList(*ws->var("n_1_A")));

  ws->import(sigma_2_A);
  ws->import(alpha_2_A);
  ws->import(n_2_A);

  // --- build signal pdfs ---
  ws->factory("CBShape::cb_1_A(mass, mean, sigma_1_A, alpha_1_A, n_1_A)");
  ws->factory("CBShape::cb_2_A(mass, mean, sigma_2_A, alpha_2_A, n_2_A)");
  ws->factory("AddPdf::pdfMASS_Jpsi({cb_1_A,cb_2_A},{f})");

  // --- build result fit model ---
  // yields
  ws->factory("N_Jpsi[2000000,1000000,5000000]");
  
  // pdf
  ws->factory("AddPdf::pdfMASS_Tot({pdfMASS_Jpsi},{N_Jpsi})");
}

void McMassFit::doFit()
{
  cout << "--- doFit() ---\n\n";
  bool isWeighted = ws->data("dsAB")->isWeighted();
  fitResult = ws->pdf("pdfMASS_Tot")->fitTo(*ws->data("dsAB"), Save(), Range(massLow, fitLimit), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
}

void McMassFit::drawPlot()
{
  cout << "--- drawPlot() ---\n\n";
  // --- canvas and pad ---
  TCanvas *c_A = new TCanvas("canvas_A", "My plots", 800, 800);
  c_A->cd();

  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0.0, 0.3, 1.0, 1.0);
  pad_A_1->SetBottomMargin(0.001);
  pad_A_1->SetTicks(1, 1);
  pad_A_1->Draw();
  pad_A_1->cd();
  gPad->SetLogy();

  // --- mass frame ---
  RooPlot *myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  myPlot_A->SetTitle("");

  // --- plotOn ---
  ws->data("dsAB")->plotOn(myPlot_A, Name("dataOS"), MarkerSize(.8));
  // bool isWeighted = ws->data("dsAB")->isWeighted();

  // Check and get fitted parameters
  const RooArgList &fitParams = fitResult->floatParsFinal();

  auto &fitN_Jpsi = (RooRealVar &)fitParams[0];
  auto &fitAlpha = (RooRealVar &)fitParams[1];
  auto &fitFraction = (RooRealVar &)fitParams[2];
  auto &fitMean = (RooRealVar &)fitParams[3];
  auto &fitN_1 = (RooRealVar &)fitParams[4];
  auto &fitSigma_1 = (RooRealVar &)fitParams[5];
  auto &fitX = (RooRealVar &)fitParams[6];
  // cout << fitN_Jpsi.getVal() << endl;

  /*
  for ( int i = 0; i < fitParams.getSize(); ++i)
  {
    auto & fitPar = (RooRealVar &) fitParams[i];
    std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
  }
  */

  double f_factor = (double)fitFraction.getVal();

  // pdfs
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A, Name("pdfMASS_tot"), LineColor(kBlack), Range(2.6, fitLimit));
  // ws->pdf("cb_1_A")->plotOn(myPlot_A,Name("cb_1_A"), LineColor(kBlue+2), Range(3.4, 3.87), Normalization(fitFraction.getVal()));
  // ws->pdf("cb_2_A")->plotOn(myPlot_A,Name("cb_2_A"), LineColor(kGreen+2), Range(3.4, 3.87), Normalization(1-fitFraction.getVal()));
  ws->pdf("cb_1_A")->plotOn(myPlot_A, Name("cb_1_A"), LineColor(kBlue + 2), Normalization(fitFraction.getVal()), Range(2.6, fitLimit));
  ws->pdf("cb_2_A")->plotOn(myPlot_A, Name("cb_2_A"), LineColor(kGreen + 2), Normalization(1 - fitFraction.getVal()), Range(2.6, fitLimit));
  // ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A,Name("Sig_A"),Components(RooArgSet(*pdfMASS_Jpsi)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));


  // --- find y-max ---
  TH1 *h = ws->data("dsAB")->createHistogram("hist", *ws->var("mass"), Binning(myPlot_A->GetNbinsX(), myPlot_A->GetXaxis()->GetXmin(), myPlot_A->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i = 1; i <= h->GetNbinsX(); i++)
    if (h->GetBinContent(i) > 0)
      YMin = min(YMin, h->GetBinContent(i));
  double Ydown;
  double Yup;
  Ydown = YMin / (TMath::Power((YMax / YMin), (0.1 / (1.0 - 0.1 - 0.4))));
  Yup = YMax * TMath::Power((YMax / YMin), (0.4 / (1.0 - 0.1 - 0.4)));
  myPlot_A->GetYaxis()->SetRangeUser(Ydown, Yup);

  // --- cosmetics ---
  myPlot_A->SetFillStyle(4000);
  myPlot_A->GetYaxis()->SetTitleOffset(1.43);
  // myPlot_A->GetYaxis()->CenterTitle();
  // myPlot_A->GetYaxis()->SetTitleSize(0.058);
  // myPlot_A->GetYaxis()->SetLabelSize(0.054);
  // myPlot_A->GetYaxis()->SetRangeUser(ws->var("N_Jpsi")->getVal()/100, ws->var("N_Jpsi")->getVal());
  myPlot_A->GetXaxis()->SetLabelSize(0);
  myPlot_A->GetXaxis()->SetTitleSize(0);
  myPlot_A->GetXaxis()->CenterTitle();
  myPlot_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot_A->Draw();

  // --- print latex ---
  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c", ptLow, ptHigh), text_x+0.033, text_y+0.08, text_color, text_size);
  if (yLow == 0)
    drawText(Form("|y^{#mu#mu}| < %.1f", yHigh), text_x+0.033, text_y+0.08 - y_diff, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < |y^{#mu#mu}| < %.1f", yLow, yHigh), text_x+0.033, text_y+0.08 - y_diff, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError()), text_x+0.033, text_y+0.08 - y_diff * 2, text_color, text_size);
  // drawText(Form("n_{Bkg} = %.f #pm %.f",ws->var("N_Bkg")->getVal(),ws->var("N_Bkg")->getError()),text_x+0.033,text_y+0.08-y_diff*4,text_color,text_size);
  drawText(Form("#alpha = %.4f #pm %.4f", ws->var("alpha_1_A")->getVal(), ws->var("alpha_1_A")->getError()), text_x+0.033, text_y+0.08 - y_diff * 3, text_color, text_size);
  drawText(Form("f = %.4f #pm %.4f", fitFraction.getVal(), fitFraction.getError()), text_x+0.033, text_y+0.08 - y_diff * 4, text_color, text_size);
  drawText(Form("n_{1} = %.4f #pm %.4f", fitN_1.getVal(), fitN_1.getError()), text_x+0.033, text_y+0.08 - y_diff * 5, text_color, text_size);
  drawText(Form("#sigma_{1} = %.4f #pm %.4f", fitSigma_1.getVal(), fitSigma_1.getError()), text_x+0.033, text_y+0.08 - y_diff * 6, text_color, text_size);
  drawText(Form("#sigma_{2} / #sigma_{1} = %.4f #pm %.4f", fitX.getVal(), fitX.getError()), text_x+0.033, text_y+0.08 - y_diff * 7, text_color, text_size);

  // --- pull pad ---
  TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0.0, 0.0, 1.0, 0.3);
  c_A->cd();
  pad_A_2->Draw();
  pad_A_2->cd();
  pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
  pad_A_2->SetBottomMargin(0.67);
  pad_A_2->SetBottomMargin(0.4);
  pad_A_2->SetFillStyle(4000);
  pad_A_2->SetFrameFillStyle(4000);
  pad_A_2->SetTicks(1, 1);

  RooHist *hpull_A = myPlot_A->pullHist("dataOS", "pdfMASS_tot", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot *pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame_A->addPlotable(hpull_A, "P");
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_A->GetYaxis()->SetTitle("Pull");
  pullFrame_A->GetYaxis()->SetTitleSize(0.15);
  pullFrame_A->GetYaxis()->SetLabelSize(0.15);
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_A->GetXaxis()->SetLabelSize(0.15);
  pullFrame_A->GetXaxis()->SetTitleSize(0.15);
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw();

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, myPlot_A, fitResult, "mass", "dataOS", "pdfMASS_tot", nMassBin, false);

  c_A->Update();
  c_A->SaveAs(Form("figs/2DFit_%s/mc_Mass/mc_Mass_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));

  delete h;
  delete c_A;
}

void McMassFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
  
  // --- hist for pT reweigh ---
  TH1D *outh = new TH1D("fitResults", "fit result", 20, 0, 20);
  outh->GetXaxis()->SetBinLabel(1, "Jpsi");

  float temp1 = ws->var("N_Jpsi")->getVal();
  float temp1err = ws->var("N_Jpsi")->getError();

  outh->SetBinContent(1, temp1);
  outh->SetBinError(1, temp1err);


  // --- dataset for mass fit ---
  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);

  // --- declare outFile ---
  TFile *outFile = new TFile(Form("roots/2DFit_%s/mc_Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");

  ws->pdf("pdfMASS_Tot")->Write();
  datasetMass->Write();
  outh->Write();

  outFile->Close();

  // --- print result ---
  fitResult->Print("V");
  // Double_t theNLL = fitResult->minNll();
  // cout << " *** NLL : " << theNLL << endl;

  delete datasetMass;
  delete fitargs;
  delete outh;
  delete outFile;
}

// =================================
// ===== legacy codes or notes =====
// =================================