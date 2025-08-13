#include "MassFit.h"
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

MassFit::MassFit(float ptLow, float ptHigh,
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

MassFit::~MassFit()
{
  delete fInputData;
  delete fMc;
  delete fitResult;
  delete ws;
}

void MassFit::init()
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

void MassFit::run()
{
  cout << "--- run() ---\n\n";
  doFit();
  rootlogon(isLogon);
  drawPlot();
  saveOutput();
}

void MassFit::setLabels()
{
  cout << "--- setLabels() ---\n\n";
  kineLabel = getKineLabel(ptLow, ptHigh, yLow, yHigh, 0.0, cLow, cHigh);

  if (PRw == 1) bCont = "PR";
  else if (PRw == 2) bCont = "NP";
}

void MassFit::makeOutputFolder()
{
  cout << "--- makeOutputFolder() ---\n\n";
  gSystem->mkdir(Form("roots/2DFit_%s/Mass", DATE.Data()), kTRUE);
  gSystem->mkdir(Form("figs/2DFit_%s/Mass", DATE.Data()), kTRUE);
}

void MassFit::turnOffRooFitMessage()
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

void MassFit::openInput()
{
  cout << "--- openInput() ---\n\n";
  fInputData = new TFile(Form("/disk1/Oniatree/miniAOD/Run2025OO/OniaRooDataSet_miniAOD_2025OORun_isMC0_Charmonia_Effw0_Accw0_PtW0_TnP0_250722.root"));
  fMc = new TFile(Form("roots/2DFit_5p36TeV_OO/mc_Mass/mc_MassFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));
}

void MassFit::processDataset()
{
  cout << "--- processDataset() ---\n\n";
  // --- kinematic cuts ---
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f", ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&"; // 2018 acceptance cut
  TString OS = "recoQQsign==0 &&";
  kineCut = OS + accCut + kineCut;

  // declare workspace
  ws = new RooWorkspace("workspace");

  // --- load datasets ---
  RooDataSet *dataset = (RooDataSet *)fInputData->Get("dataset");
  ws->import(*dataset);

  // isWeight?
  RooDataSet *datasetW = nullptr;
  if (isWeighted)
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset) ,WeightVar(*ws->var("weight")));
  else
    datasetW = new RooDataSet("datasetW", "A sample", *dataset->get(), Import(*dataset));
  
  // apply kinecuts
  RooDataSet *dsAB = (RooDataSet *)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")), *(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data());
  ws->import(*dsAB, Rename("dsAB"));

  delete dsAB; 
  delete datasetW;
}

void MassFit::setVariableRanges()
{
  cout << "--- setVariableRanges() ---\n\n";
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
}

void MassFit::initVar(const std::string &varName, double init, double low, double high)
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

void MassFit::buildModel()
{
  cout << "--- buildModel() ---\n\n";
  // --- get MC fit parameters ---
  RooDataSet *dsMc = (RooDataSet *)fMc->Get("datasetMass");
  RooWorkspace *ws_fit = new RooWorkspace("workspace_fit");
  ws_fit->import(*dsMc);

  Double_t alpha_MC_value = ws_fit->var("alpha_1_A")->getVal();
  Double_t alpha_MC_value_err = ws_fit->var("alpha_1_A")->getError();
  Double_t n_MC_value = ws_fit->var("n_1_A")->getVal();
  Double_t n_MC_value_err = ws_fit->var("n_1_A")->getError();
  Double_t xA_MC_value = ws_fit->var("x_A")->getVal();
  Double_t xA_MC_value_err = ws_fit->var("x_A")->getError();
  Double_t f_MC_value = ws_fit->var("f")->getVal();
  Double_t f_MC_value_err = ws_fit->var("f")->getError();
  Double_t sigma_MC_value = ws_fit->var("sigma_1_A")->getVal();
  Double_t sigma_MC_value_err = ws_fit->var("sigma_1_A")->getError();

  double sigma_index = 5;

  Double_t alpha_lower = alpha_MC_value - (sigma_index * alpha_MC_value_err);
  Double_t alpha_higher = alpha_MC_value + (sigma_index * alpha_MC_value_err);
  Double_t xA_lower = xA_MC_value - (sigma_index * xA_MC_value_err);
  Double_t xA_higher = xA_MC_value + (sigma_index * xA_MC_value_err);
  Double_t n_lower = n_MC_value - (sigma_index * n_MC_value_err);
  Double_t n_higher = n_MC_value + (sigma_index * n_MC_value_err);
  Double_t f_lower = f_MC_value - (sigma_index * f_MC_value_err);
  Double_t f_higher = f_MC_value + (sigma_index * f_MC_value_err);
  Double_t sigma_lower = sigma_MC_value - (sigma_index * sigma_MC_value_err);
  Double_t sigma_higher = sigma_MC_value + (sigma_index * sigma_MC_value_err);

  // if (n_lower<0.0)n_lower==0.0;
  if (f_higher > 1.0)
    f_higher = 1.0;
  if (f_lower < 0.0)
    f_lower = 0.0;

  // --- bulid signal model ---
  // set parameters
  double paramslower[6] = {alpha_lower, n_lower, 0.0, xA_lower, 0.0, -25.0};
  double paramsupper[6] = {alpha_higher, n_higher, 0.6, xA_higher, 1.0, 25.0};

  double alpha_1_init = alpha_MC_value;
  double n_1_init = n_MC_value;
  double sigma_1_init = sigma_MC_value;
  double x_init = xA_MC_value;
  double f_init = f_MC_value;

  ws->factory(Form("mean[%f,%f,%f]",
                   pdgMass.JPsi, pdgMass.JPsi - 0.05, pdgMass.JPsi + 0.05));
  ws->factory(Form("x_A[%f]", x_init));
  ws->factory(Form("sigma_1_A[%f,%f,%f]",
                   sigma_1_init, sigma_1_init - 0.01, sigma_1_init + 0.01));
  ws->factory(Form("alpha_1_A[%f]", alpha_1_init));
  ws->factory(Form("n_1_A[%f]", n_1_init));
  ws->factory(Form("f[%f]", f_init));

  // 2) (중요) RooFormulaVar로 수식 변수 생성 → ws에 import
  RooFormulaVar sigma_2_A("sigma_2_A", "@0*@1",
                          RooArgList(*ws->var("sigma_1_A"), *ws->var("x_A")));
  ws->import(sigma_2_A);

  RooFormulaVar alpha_2_A("alpha_2_A", "1.0*@0",
                          RooArgList(*ws->var("alpha_1_A")));
  ws->import(alpha_2_A);

  RooFormulaVar n_2_A("n_2_A", "1.0*@0",
                      RooArgList(*ws->var("n_1_A")));
  ws->import(n_2_A);

  // bulid PDFS
  ws->factory("CBShape::cb_1_A(mass, mean, sigma_1_A, alpha_1_A, n_1_A)");
  ws->factory("CBShape::cb_2_A(mass, mean, sigma_2_A, alpha_2_A, n_2_A)");
  ws->factory("SUM::pdfMASS_Jpsi(f*cb_1_A, cb_2_A)");

  // --- build bkg model ---
  // set parameters
  ws->factory("sl1[0.0,-1,1]");
  ws->factory("sl2[0.0,-1,1]");
  ws->factory("RooChebychev::pdfMASS_bkg(mass, {sl1, sl2})");

  // --- result fit model ---
  // yields
  Double_t NBkg_limit = 2.0e+07;
  Double_t NJpsi_limit = 2.0e+07;

  ws->factory(Form("N_Jpsi[%f, 1000, 1000000]", NJpsi_limit));
  ws->factory(Form("N_Bkg[%f, 1000, 1000000]", NBkg_limit));

  ws->factory("SUM::pdfMASS_Tot(N_Jpsi*pdfMASS_Jpsi, N_Bkg*pdfMASS_bkg)");

  delete ws_fit;
}

void MassFit::doFit()
{
  cout << "--- doFit() ---\n\n";
  bool isWeighted = ws->data("dsAB")->isWeighted();
  fitResult = ws->pdf("pdfMASS_Tot")->fitTo(*ws->data("dsAB"), Save(), Hesse(kTRUE), Range(massLow, massHigh), Timer(kTRUE), Extended(kTRUE), SumW2Error(isWeighted), NumCPU(nCPU));
}

void MassFit::drawPlot()
{
  cout << "--- drawPlot() ---\n\n";
  // --- canvas and pad ---
  TCanvas *c_A = new TCanvas("c_A", "My plots", 800, 800);
  c_A->cd();

  TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0.0, 0.3, 1.0, 1.0);
  pad_A_1->SetBottomMargin(0);
  pad_A_1->SetTicks(1, 1);
  pad_A_1->Draw();
  pad_A_1->SetLogy();

  // --- mass frames ---
  pad_A_1->cd();

  RooPlot *myPlot_A = (RooPlot *)ws->var("mass")->frame(Bins(nMassBin), Range(massLow, massHigh));
  myPlot_A->SetTitle("");

  // --- plotOn ---
  ws->data("dsAB")->plotOn(myPlot_A, Name("dataOS"), MarkerSize(.8));

  ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A, Name("pdfMASS_Tot"), LineColor(kBlack));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A, Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cb_1_A"))), LineColor(44));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A, Components(RooArgSet(*ws->pdf("pdfMASS_bkg"), *ws->pdf("cb_2_A"))), LineColor(8));
  ws->pdf("pdfMASS_Tot")->plotOn(myPlot_A, Name("pdfMASS_bkg"), Components(RooArgSet(*ws->pdf("pdfMASS_bkg"))), LineColor(kBlue), LineStyle(kDashed), LineWidth(2));
  ws->data("dsAB")->plotOn(myPlot_A, Name("dataOS"), MarkerSize(.8)); // draw points above lines

  // --- find good y-max ---
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
  
  // --- cosmetics ---
  myPlot_A->GetYaxis()->SetRangeUser(Ydown, Yup);
  // myPlot_A->SetMinimum(2*10);
  myPlot_A->SetFillStyle(4000);
  myPlot_A->GetYaxis()->SetTitleOffset(1.43);
  myPlot_A->GetXaxis()->SetLabelSize(0);
  myPlot_A->GetXaxis()->SetTitleSize(0);
  myPlot_A->GetXaxis()->CenterTitle();
  myPlot_A->GetXaxis()->SetRangeUser(massLow, massHigh);
  myPlot_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot_A->Draw();

  // --- legend ---
  TLegend *leg_B = new TLegend(text_x + 0.033 + 0.5, text_y + 0.12 - 0.2, text_x + 0.033 + 0.7, text_y + 0.12);
  leg_B->SetTextSize(text_size);
  leg_B->SetTextFont(43);
  leg_B->SetBorderSize(0);
  leg_B->AddEntry(myPlot_A->findObject("dataOS"), "Data", "pe");
  leg_B->AddEntry(myPlot_A->findObject("pdfMASS_Tot"), "Total", "l");
  leg_B->AddEntry(myPlot_A->findObject("pdfMASS_bkg"), "Background", "l");
  leg_B->Draw("same");

  // --- draw latex ---
  if (yLow == 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c, |y^{#mu#mu}| < %.1f, Cent. %d - %d%s ", ptLow, ptHigh, yHigh, cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.08, text_color, text_size);
  else if (yLow != 0)
    drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c; %.1f < |y^{#mu#mu}| < %.1f; Cent. %d - %d%s", ptLow, ptHigh, yLow, yHigh, cLow / 2, cHigh / 2, "%"), text_x + 0.033, text_y + 0.08, text_color, text_size);
  drawText(Form("N_{J/#psi} = %.f #pm %.f,  N_{Bkg} = %.f #pm %.f", ws->var("N_Jpsi")->getVal(), ws->var("N_Jpsi")->getError(), ws->var("N_Bkg")->getVal(), ws->var("N_Bkg")->getError()), text_x + 0.033, text_y + 0.08 - y_diff * 1, text_color, text_size);
  drawText(Form("m_{J/#psi} = %.4f #pm %.4f", ws->var("mean")->getVal(), ws->var("mean")->getError()), text_x + 0.033, text_y + 0.08 - y_diff * 2, text_color, text_size);
  drawText(Form("#alpha_{J/#psi} = %.4f (fixed)", ws->var("alpha_1_A")->getVal()), text_x + 0.033, text_y + 0.08 - y_diff * 3, text_color, text_size);
  drawText(Form("f_{J/#psi} = %.4f (fixed)", ws->var("f")->getVal()), text_x + 0.033, text_y + 0.08 - y_diff * 4, text_color, text_size);
  drawText(Form("n_{J/#psi} = %.4f (fixed)", ws->var("n_1_A")->getVal()), text_x + 0.033, text_y + 0.08 - y_diff * 5, text_color, text_size);
  drawText(Form("#sigma1_{J/#psi} = %.2f #pm %.2f MeV/c^{2}, (#sigma2/#sigma1)_{J/#psi} = %.3f (fixed)", (ws->var("sigma_1_A")->getVal()) * 1000, (ws->var("sigma_1_A")->getError()) * 1000, ws->var("x_A")->getVal()), text_x + 0.033, text_y + 0.08 - y_diff * 6, text_color, text_size);

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

  // RooPlot *frameTMP = (RooPlot *)myPlot_A->Clone("TMP");
  // RooHist *hpull_A = frameTMP->pullHist("dataOS", "pdfMASS_Tot", true);

  RooHist *hpull_A = myPlot_A->pullHist("dataOS", "pdfMASS_Tot", true);
  hpull_A->SetMarkerSize(0.8);
  RooPlot *pullFrame_A = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame_A->addPlotable(hpull_A, "P");
  pullFrame_A->SetTitle("");
  pullFrame_A->SetTitleSize(0);
  pullFrame_A->GetYaxis()->SetTitleOffset(0.3);
  pullFrame_A->GetYaxis()->SetTitle("Pull");
  pullFrame_A->GetYaxis()->SetTitleSize(0.08);
  pullFrame_A->GetYaxis()->SetLabelSize(0.08);
  pullFrame_A->GetYaxis()->SetRangeUser(-3.8, 3.8);
  pullFrame_A->GetYaxis()->CenterTitle();

  pullFrame_A->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame_A->GetXaxis()->SetTitleOffset(1.05);
  pullFrame_A->GetXaxis()->SetLabelOffset(0.04);
  pullFrame_A->GetXaxis()->SetLabelSize(0.08);
  pullFrame_A->GetXaxis()->SetTitleSize(0.08);
  pullFrame_A->GetXaxis()->CenterTitle();

  pullFrame_A->GetYaxis()->SetTickSize(0.04);
  pullFrame_A->GetYaxis()->SetNdivisions(404);
  pullFrame_A->GetXaxis()->SetTickSize(0.03);
  pullFrame_A->Draw();

  TLine *l1 = new TLine(massLow, 0, massHigh, 0);
  l1->SetLineStyle(1);
  l1->Draw("same");

  printChi2(ws, pad_A_2, myPlot_A, fitResult, "mass", "dataOS", "pdfMASS_Tot", nMassBin, false);

  c_A->Update();
  c_A->SaveAs(Form("figs/2DFit_%s/Mass/Mass_Fixed_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.pdf", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP));

  delete h;
  delete c_A;
}

void MassFit::saveOutput()
{
  cout << "--- saveOutput() ---\n\n";
  TFile *outFile = new TFile(Form("roots/2DFit_%s/Mass/Mass_FixedFitResult_%s_%sw_Effw%d_Accw%d_PtW%d_TnP%d.root", DATE.Data(), kineLabel.Data(), bCont.Data(), fEffW, fAccW, isPtW, isTnP), "recreate");
  outFile->cd();

  ws->pdf("pdfMASS_Tot")->Write();
  ws->pdf("pdfMASS_bkg")->Write();

  // dataset for next step - Can a fitResult replace it?
  RooArgSet *fitargs = new RooArgSet();
  fitargs->add(fitResult->floatParsFinal());
  // std::unique_ptr<RooArgSet> fitargs(fitResult->floatParsFinal().snapshot());
  RooDataSet *datasetMass = new RooDataSet("datasetMass", "dataset with Mass Fit result", *fitargs);
  datasetMass->add(*fitargs);
  datasetMass->Write();

  // legacy for pT reweighting -> Can a fitResult replace it?
  TH1D *outh = new TH1D("fitResults", "fit result", 20, 0, 20);
  outh->GetXaxis()->SetBinLabel(1, "Jpsi");

  float temp1 = ws->var("N_Jpsi")->getVal();
  float temp1err = ws->var("N_Jpsi")->getError();

  outh->SetBinContent(1, temp1);
  outh->SetBinError(1, temp1err);
  outh->Write();

  fitResult->Write();

  outFile->Close();

  // --- print fit result ---
  fitResult->Print("V");
  // Double_t theNLL = fitResult->minNll();
  // cout << " *** NLL : " << std::setprecision(15) << theNLL << endl;

  // delete outh;
  delete fitargs;
  delete datasetMass;
  delete outFile;
}

// =================================
// ===== legacy codes or notes =====
// =================================
