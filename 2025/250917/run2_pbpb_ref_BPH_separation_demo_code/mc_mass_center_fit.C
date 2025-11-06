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
#include <RooArgList.h> // RooFormulaVar parameter list
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooCrystalBall.h>

using namespace RooFit;
using std::cout;
using std::string;

void mc_mass_center_fit()
{
  float ptLow = 25, ptHigh = 40;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "PR", region = "MassFull";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start mc_mass_center_fit() ===\n";

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs", true);
  gSystem->mkdir("roots", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_miniAOD_isMC1_JPsi_Prompt_cent0_200_Effw0_Accw0_PtW0_TnP0_250421.root");

  // OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root
  // OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root
  // OniaRooDataSet_JPsi_pp_GENONLY_NonPrompt_230215.root

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

  // opposite sign -> Must be FLowSkim later
  TString osCut = "(recoQQsign == 0)";

  // --- region6 cuts ---
  const TString cutPR = "(abs(ctau3D) < 0.05)";
  const TString cutNP = "(ctau3D >= 0.10 && ctau3D <= 0.80)";
  const TString cutCtauFull = "(ctau3D >= -.10 && ctau3D <= 0.80)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";
  const TString cutMassFull = "(mass >= 2.6 && mass <= 3.5)";

  // --- combine cuts ---
  TString compCut;
  if (comp == "PR")
    compCut = cutPR;
  else if (comp == "NP")
    compCut = cutNP;
  else
    compCut = cutCtauFull;

  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else if (region == "MassFull")
    regionCut = cutMassFull;
  else regionCut = "(1)";

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
  RooRealVar *massVar = dynamic_cast<RooRealVar *>(ds_red->get()->find("mass"));
  if (!massVar)
    cout << "Warn: There is no variable 'mass'\n";
  massVar->setRange(2.6, 3.5); // Change PDF's range - Need.
  massVar->setRange("fitRange", 2.6, 3.5);
  

  // === build models ===
  // --- parameters ---
  double mean0 = 3.0969;
  RooRealVar mean("mean", "signal mean", mean0, mean0 - 0.05, mean0 + 0.05);

  // --- CB left tail ---
  RooRealVar sigmaL("sigmaL", "sigma left", 0.025, 0.01, 0.05);
  RooRealVar alphaL("alphaL", "alpha left", 1.4, 0.2, 2.0);
  RooRealVar nL("nL", "n left", 2.3, 1, 5.0);
  RooCBShape cbLeft("cbLeft", "CB left tail", *massVar, mean, sigmaL, alphaL, nL);

  // --- CB right taile ---
  RooRealVar sigmaRatio("sigmaRatio", "sigma right", 1.0, 0.5, 3.0);
  RooFormulaVar sigmaR("sigmaR", "sigma right", "@0*@1", RooArgList(sigmaL, sigmaRatio));
  // RooRealVar alphaR("alphaR", "alpha right", -1.5, -5.0, -0.01);
  RooRealVar alphaRatio("alphaRatio", "alphaR/alphaL ratio", 3.6, 0.5, 10.0);
  RooFormulaVar alphaR("alphaR", "alpha right", "-@0*@1", RooArgList(alphaL, alphaRatio)); // negative sign for right side
  RooRealVar nRatio("nRatio", "nR/nL ratio", 4.2, 3.7, 5);
  RooFormulaVar nR("nR", "n right", "@0*@1", RooArgList(nL, nRatio));
  RooCBShape cbRight("cbRight", "CB right tail", *massVar, mean, sigmaR, alphaR, nR);

  // --- gauss ---
  RooRealVar sigmaGRatio("sigmaGRatio", "", 1.6, 0.5, 3.0);
  RooFormulaVar sigmaG("sigmaG", "sigma gauss", "@0*@1", RooArgList(sigmaL, sigmaGRatio));
  // RooRealVar sigmaG("sigmaG", "sigma gaussian", 0.015, 0.001, 0.080);
  RooGaussian gaus("gaus", "gaus core", *massVar, mean, sigmaG);

  // --- combine models ---
  RooRealVar f_cbR("f_cbR", "frac of right CB", 0.5, 0, 1.0);
  RooRealVar f_gaus("f_gaus", "frac of gaus", 0.05, 0, 1);
  RooAddPdf model("sig", "", RooArgList(cbLeft, cbRight, gaus), RooArgList(f_cbR, f_gaus));


  // === perform fit ===
  auto dh = new RooDataHist("dh", "binned dataset", *massVar, *ds_red);

  auto fitResult = model.fitTo(*dh, Save(), Range("fitRange"), SumW2Error(hasWeight), Offset(true), PrintLevel(-1), Warnings(kFALSE), Verbose(kFALSE), Strategy(2), NumCPU(32), EvalBackend("legacy"));

  // === draw ===
  // --- divided canvas ---
  TCanvas c_mass("c_mass", "c_mass", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  double massMin = 2.6, massMax = 3.5;
  RooPlot *massFrame = massVar->frame(Range(massMin, massMax), Title("")); // Bins(80)
  dh->plotOn(massFrame, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(massFrame, NormRange("fitRange"), Range("fitRange"), Name("model"));
  model.plotOn(massFrame, Components("cbLeft"), LineStyle(kDotted), LineColor(kRed), Name("cbLeft"));
  model.plotOn(massFrame, Components("cbRight"), LineStyle(kDotted), LineColor(kAzure), Name("cbRight"));
  model.plotOn(massFrame, Components("gaus"), LineStyle(kDotted), LineColor(kViolet), Name("gaus"));
  

  // y axis: logY style
  double ymin = 1e300, ymax = -1e300;
  RooHist *hdata = (RooHist *)massFrame->getHist("data"); // use first dataset on massFrame
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

  massFrame->SetMinimum(ymin * 0.5);
  massFrame->SetMaximum(ymax * 500.0);

  // title
  massFrame->GetYaxis()->SetTitle("Events");
  massFrame->GetXaxis()->SetTitle("");
  massFrame->Draw("e");

  // --- object legend ---
  TLegend *leg = new TLegend(0.49, 0.65, 0.70, 0.93);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(massFrame->findObject("data"), "Data", "lep");
  leg->AddEntry(massFrame->findObject("model"), "Fit model", "pe");
  leg->AddEntry(massFrame->findObject("cbLeft"), "CB1", "pe");
  leg->AddEntry(massFrame->findObject("cbRight"), "CB2", "pe");
  leg->AddEntry(massFrame->findObject("gaus"), "Gauss", "pe");
  leg->Draw("same");

  // --- info latex ---
  TLatex latexInfo;
  latexInfo.SetNDC();
  latexInfo.SetTextSize(0.03);
  latexInfo.SetTextFont(42);

  double x_start = 0.19;
  double y_start = 0.95;
  double y_step = -0.06, y_stepCount = 1;
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, "CMS PbPb. #sqrt{s_{NN}} = 5.02 TeV");
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, "Prompt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
  if (yLow ==0 )
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
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("mean = %.3f #pm %.3f", mean.getVal(), mean.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#alpha_{L} = %.3f #pm %.3f", alphaL.getVal(), alphaL.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#alpha_{ratio} = %.3f #pm %.3f", alphaRatio.getVal(), alphaRatio.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{CB,R} = %.3f #pm %.3f", f_cbR.getVal(), f_cbR.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("f_{Gauss} = %.3f #pm %.3f", f_gaus.getVal(), f_gaus.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("n_{L} = %.3f #pm %.3f", nL.getVal(), nL.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("n_{ratio} = %.3f #pm %.3f", nRatio.getVal(), nRatio.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma_{L} = %.3f #pm %.3f", sigmaL.getVal(), sigmaL.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma_{ratio, R} = %.3f #pm %.3f", sigmaRatio.getVal(), sigmaRatio.getError()));
  latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma_{ratio, Gauss} = %.3f #pm %.3f", sigmaGRatio.getVal(), sigmaGRatio.getError()));
  
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("sigma_{L} = %.3f #pm %.3f", sigmaL.getVal(), sigmaL.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("sigma_{L} = %.3f #pm %.3f", sigmaL.getVal(), sigmaL.getError()));

  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#alpha_{L} = %.3f #pm %.3f, #alpha_{Ratio} = %.3f #pm %.3f", alphaL.getVal(), alphaL.getError(), alphaRatio.getVal(), alphaRatio.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#f_{CB,R} = %.3f #pm %.3f, #f_{gauss} = %.3f #pm %.3f", f_cbR.getVal(), f_cbR.getError(), f_gaus.getVal(), f_gaus.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#n_{L} = %.3f #pm %.3f, #n_{ratio} = %.3f #pm %.3f", nL.getVal(), nL.getError(), nRatio.getVal(), nRatio.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma_{ratio} = %.3f #pm %.3f, #sigma_{Gauss, ratio} = %.3f #pm %.3f", sigmaRatio.getVal(), sigmaRatio.getError(), sigmaGRatio.getVal(), sigmaGRatio.getError()));

  // === pull pad ===
  c_mass.cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  RooHist *hpull = massFrame->pullHist("data", "model");
  RooPlot *f_pull = massVar->frame(Range(massMin, massMax), Title(""));
  f_pull->addPlotable(hpull, "P"); // P: points only

  f_pull->GetYaxis()->SetTitle("Pull");
  f_pull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
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
  double xmin = massMin;
  double xmax = massMax;
  TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
  // line->SetLineColor();
  line->SetLineStyle(2);
  line->Draw("same");

  // --- compute and draw chi square ---
  int nFitParam = fitResult->floatParsFinal().getSize();
  double chi2ndf = massFrame->chiSquare("model", "data", nFitParam);

  TLatex latexChi2;
  latexChi2.SetNDC(); // use pad coordinates (0~1)
  latexChi2.SetTextSize(0.1);
  latexChi2.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  c_mass.SaveAs(Form("figs/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.png", ptLow, ptHigh, yLow, yHigh));
  c_mass.SaveAs(Form("figs/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.pdf", ptLow, ptHigh, yLow, yHigh));

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  TFile fout(Form("roots/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "RECREATE");
  fitResult->Write("fitResult");
  model.Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}