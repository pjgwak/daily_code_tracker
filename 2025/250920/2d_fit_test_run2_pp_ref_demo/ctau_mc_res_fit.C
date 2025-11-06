#include <iostream> // std::cout
#include <cstdio>   // printf
#include <string>
#include <TStopwatch.h>
#include <TSystem.h> // gSystem
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>

#include <RooGlobalFunc.h> // using namespace RooFit;
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>     // RooFormulaVar parameter list
#include <RooExponential.h>
#include <RooFitResult.h>

using namespace RooFit;
using std::cout;
using std::string;

void ctau_mc_res_fit()
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "", region = "";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_mc_res_fit() ===\n";
  
  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs", true);
  gSystem->mkdir("roots", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root");
  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *ds = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!ds)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }

  // === declare cuts ===
  // --- basic cuts ---
  // acceptance
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut

  // kinematics cuts
  //  - correct? (<= cBin <) -> maybe (< cBin <=) ??
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  TString osCut = "(recoQQsign == 0)";

  // --- region6 cuts ---
  const TString cutPR = "(abs(ctau3D) < 0.05)";
  const TString cutNP = "(ctau3D >= 0.10 && ctau3D <= 0.80)";
  const TString cutCtauFull = "(ctau3D >= -0.80 && ctau3D <= 0.80)";
  const TString cutSide = "(ctau3D >= 0.10 && ctau3D <= 0.50)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";

  // --- combine cuts ---
  TString compCut;
  if (comp == "PR")
    compCut = cutPR;
  else if (comp == "NP")
    compCut = cutNP;
  else if (comp == "Side")
    compCut = cutSide;
  else
    compCut = "(1)";

  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else
    regionCut = "(1)";

  TString fullCut = Form("%s && %s && %s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data(),
                         compCut.Data(),
                         regionCut.Data());

  // === new dataset with cuts ===
  RooDataSet *ds_red = (RooDataSet *)ds->reduce(Cut(fullCut));
  if (!ds_red || ds_red->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  ds_red->Print();

  // check weight
  const bool hasWeight = (ds_red->isWeighted() || ds_red->weightVar() || ds_red->get()->find("weight"));
  // const bool hasWeight = false;

  // variables
  RooRealVar *ctau3DVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3DRes"));
  if (!ctau3DVar)
    cout << "Warn: There is no variable 'ctau3D'\n";
  ctau3DVar->setRange(-10, 10); // Change PDF's range - Need.
  ctau3DVar->setRange("fitRange", -10, 10);
  
  // === fit by Expo ===
  RooRealVar mu0("mu0", "res mean", 0, -0.02, 0.02);
  // RooRealVar mu0("mu0", "res mean", 0);

  RooRealVar sigma1("sigma1", "sigma1", 0.5, 0.1, 1);
  RooRealVar r21("r21", "sigma2/sigma1", 1.1, 1, 2.1);
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, r21));

  RooRealVar r32("r32", "sigma3/sigma2", 2.6, 1, 3);
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, r32));

  RooRealVar fg1("fg1", "frac g1", 0.6, 0.5, 1.00);
  RooRealVar fg2("fg2", "frac g2", 0.2, 0.1, 1.00);

  RooGaussModel g1("g1", "g1", *ctau3DVar, mu0, sigma1);
  RooGaussModel g2("g2", "g2", *ctau3DVar, mu0, sigma2);
  RooGaussModel g3("g3", "g3", *ctau3DVar, mu0, sigma3);
  // RooAddModel is more natural for RooResolutionModel class
  // Use RooAddModel when it should be inserted to other function
  // Use RooAddPdf with recursive option if the model is used outest fit model
  RooAddPdf model("model", "3-Gauss resolution model",
                  RooArgList(g1, g2, g3), RooArgList(fg1, fg2), kTRUE); // recursive fraction

  // --- perform fit ---
  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"), Offset(true),
                                       Extended(true),
                                       NumCPU(32),
                                       EvalBackend("legacy"),
                                       RooFit::RecoverFromUndefinedRegions(2.),
                                       RooFit::PrintEvalErrors(-1),
                                       RooFit::PrintLevel(-1));

  // === draw===
  // --- divided canvas ---
  TCanvas c_ctau("c_ctau", "c_ctau", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  double ctMin = -10, ctMax = 10;
  // double ctMin = -0.05, ctMax = 0.05;
  RooPlot *f_ctau = ctau3DVar->frame(Range(ctMin, ctMax), Title("")); // Bins(80)
  ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn((f_ctau), NormRange("fitRange"), Range("fitRange"), Name("model"));
  model.plotOn(f_ctau, Components("g1"), LineStyle(kDotted), LineColor(kRed), Name("g1"));
  model.plotOn(f_ctau, Components("g2"), LineStyle(kDotted), LineColor(kAzure), Name("g2"));
  model.plotOn(f_ctau, Components("g3"), LineStyle(kDotted), LineColor(kViolet), Name("g3"));

  // y axis: logY style
  double ymin = 1e300, ymax = -1e300;
  RooHist *hdata = (RooHist *)f_ctau->getHist("data"); // use first dataset on f_ctau
  if (hdata)
  {
    for (int i = 0; i < hdata->GetN(); i++)
    {
      double x, y;
      hdata->GetPoint(i, x, y);
      if (y > 0 && y < ymin)
        ymin = y;
      if (y > ymax)
        ymax = y;
    }
  }

  double floor = 1e-3;
  if (ymin <= 0 || ymin == 1e300)
    ymin = floor;

  f_ctau->SetMinimum(ymin * 0.5);
  f_ctau->SetMaximum(ymax * 500.0);

  // title
  f_ctau->GetYaxis()->SetTitle("Events");
  f_ctau->GetXaxis()->SetTitle("");
  f_ctau->Draw("e");

  // --- object legend ---
  TLegend *leg = new TLegend(0.49, 0.65, 0.70, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(f_ctau->findObject("data"), "Data", "lep");
  leg->AddEntry(f_ctau->findObject("model"), "Fit model", "pe");
  leg->AddEntry(f_ctau->findObject("g1"), "Gauss 1", "pe");
  leg->AddEntry(f_ctau->findObject("g2"), "Gauss 2", "pe");
  leg->AddEntry(f_ctau->findObject("g3"), "Gauss 3", "pe");
  leg->Draw("same");

  // --- info latex ---
  TLatex latexInfo;
  latexInfo.SetNDC();
  latexInfo.SetTextSize(0.03);
  latexInfo.SetTextFont(42);

  double x_start = 0.19;
  double y_start = 0.95;
  double y_step = -0.06, y_stepCount = 1;
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, "Prompt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
  if (yLow == 0)
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
  else
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));

  // --- latex parameters ---
  TLatex latexParams;
  latexParams.SetNDC();
  latexParams.SetTextSize(0.025);
  latexParams.SetTextFont(42);

  x_start = 0.71;
  y_step = -0.045;
  y_stepCount = 1;
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("mean = 0 (fixed)"));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("mean = %.4f #pm %.4f", mu0.getVal(), mu0.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma1 = %.3f #pm %.3f", sigma1.getVal(), sigma1.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma2/#sigma1 = %.3f #pm %.3f", r21.getVal(), r21.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma3/#sigma2 = %.3f #pm %.3f", r32.getVal(), r32.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{G1} = %.3f #pm %.3f", fg1.getVal(), fg1.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{G2} = %.3f #pm %.3f", fg2.getVal(), fg2.getError()));

  // === pull pad ===
  c_ctau.cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  RooHist *hpull = f_ctau->pullHist("data", "model");
  RooPlot *f_pull = ctau3DVar->frame(Range(ctMin, ctMax), Title(""));
  f_pull->addPlotable(hpull, "P"); // P: points only

  f_pull->GetYaxis()->SetTitle("Pull");
  f_pull->GetXaxis()->SetTitle("c#tau_{3D} Res");
  f_pull->GetXaxis()->CenterTitle();
  f_pull->SetMinimum(-8);
  f_pull->SetMaximum(8);
  f_pull->GetYaxis()->SetNdivisions(505);
  f_pull->GetYaxis()->SetTitleSize(0.12);
  f_pull->GetYaxis()->SetLabelSize(0.10);
  f_pull->GetXaxis()->SetTitleSize(0.15);
  f_pull->GetXaxis()->SetLabelSize(0.10);
  f_pull->Draw();

  // --- draw pull = 0 line ---
  double xmin = ctMin;
  double xmax = ctMax;
  TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
  // line->SetLineColor();
  line->SetLineStyle(2);
  line->Draw("same");

  // --- compute and draw chi square ---
  int nFitParam = fitResult->floatParsFinal().getSize();
  double chi2ndf = f_ctau->chiSquare("model", "data", nFitParam);

  TLatex latex;
  latex.SetNDC(); // use pad coordinates (0~1)
  latex.SetTextSize(0.1);
  latex.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  c_ctau.SaveAs(Form("figs/ctau_mc_res_pT%.1f_%.1f_y%.1f_%.1f.png", ptLow, ptHigh, yLow, yHigh));

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  TFile fout(Form("roots/ctau_mc_res_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "RECREATE");
  fitResult->Write("fitResult");
  model.Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}