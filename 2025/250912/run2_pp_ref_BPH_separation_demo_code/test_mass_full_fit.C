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


void test_mass_full_fit()
{
  float ptLow = 25, ptHigh = 40;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  std::string comp = "", region = "MassFull";
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start test_mass_full_fit() ===\n";

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs", true);
  gSystem->mkdir("roots", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root");

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
  const TString cutCtauFull = "(ctau3D >= -.10 && ctau3D <= 0.80)";

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";
  const TString cutMassCenter = "(mass >= 2.95 && mass <= 3.25)";
  const TString cutMassFull = "(mass >= 2.6 && mass <= 3.5)";

  // --- combine cuts ---
  TString compCut;
  if (comp == "PR")
    compCut = cutPR;
  else if (comp == "NP")
    compCut = cutNP;
  else if (comp == "CtauFull")
    compCut = cutCtauFull;
  else
    compCut = "(1)";

  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else if (region == "MassFull")
    regionCut = cutMassFull;
  else if (region == "MassCenter")
    regionCut = cutMassCenter;
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

  massVar->setRange(massLow, massHigh);
  massVar->setRange("fitRange", massLow, massHigh);
  double massMin = massLow, massMax = massHigh;

  // === fitting ===
  double mean0 = 3.0969;
  double sigma1_0 = 0.020;
  double sigma2_0 = 0.040;
  double sigmaG_0 = 0.015;

  // === bring MC fit results ===
  TFile fin(Form("roots/mc_mass_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *fr = (RooFitResult *)fin.Get("fitResult");
  if (!fr)
  {
    std::cerr << "fitResult not found!" << std::endl;
    return;
  }

  const RooArgList &params = fr->floatParsFinal();
  RooRealVar *sigmaLVar = (RooRealVar *)params.find("sigmaL");
  RooRealVar *sigmaGVar = (RooRealVar *)params.find("sigmaG");
  RooRealVar *alphaLVar = (RooRealVar *)params.find("alphaL");
  RooRealVar *nLVar = (RooRealVar *)params.find("nL");
  RooRealVar *sigmaRVar = (RooRealVar *)params.find("sigmaR");
  RooRealVar *alphaRatioVar = (RooRealVar *)params.find("alphaRatio");
  RooRealVar *nRatioVar = (RooRealVar *)params.find("nRatio");
  RooRealVar *sigmaRatioVar = (RooRealVar *)params.find("sigmaRatio");
  RooRealVar *sigmaGRatioVar = (RooRealVar *)params.find("sigmaGRatio");

  // if (sigmaLVar)
  // {
  //   sigmaLVar->setConstant(kTRUE);
  //   std::cout << "Fixed sigmaL = " << sigmaLVar->getVal() << std::endl;
  // }

  if (alphaLVar)
  {
    alphaLVar->setConstant(kTRUE);
    std::cout << "Fixed alphaL = " << alphaLVar->getVal() << std::endl;
  }
  if (nLVar)
  {
    nLVar->setConstant(kTRUE);
    std::cout << "Fixed nL = " << nLVar->getVal() << std::endl;
  }
  if (sigmaGVar)
  {
    sigmaGVar->setConstant(kTRUE);
    std::cout << "Fixed sigmaG = " << sigmaGVar->getVal() << std::endl;
  }
  if (alphaRatioVar)
  {
    alphaRatioVar->setConstant(kTRUE);
    std::cout << "Fixed alphaRatio = " << alphaRatioVar->getVal() << std::endl;
  }
  if (nRatioVar)
  {
    nRatioVar->setConstant(kTRUE);
    std::cout << "Fixed nRatio = " << nRatioVar->getVal() << std::endl;
  }
  if (sigmaRatioVar)
  {
    sigmaRatioVar->setConstant(kTRUE);
    std::cout << "Fixed sigmaRatio = " << sigmaRatioVar->getVal() << std::endl;
  }
  if (sigmaGRatioVar)
  {
    sigmaGRatioVar->setConstant(kTRUE);
    std::cout << "Fixed sigmaGRatio = " << sigmaGRatioVar->getVal() << std::endl;
  }

  // === build model ===
  RooRealVar mean("mean", "signal mean", mean0, mean0 - 0.050, mean0 + 0.050);
  RooCBShape cbLeft("cbLeft", "CB left tail", *massVar, mean, *sigmaLVar, *alphaLVar, *nLVar);

  // --- CB right ---
  RooFormulaVar sigmaR("sigmaR", "sigma right", "@0*@1", RooArgList(*sigmaLVar, *sigmaRatioVar));
  RooFormulaVar alphaR("alphaR", "alpha right", "-@0*@1", RooArgList(*alphaLVar, *alphaRatioVar));
  RooFormulaVar nR("nR", "n right", "@0*@1", RooArgList(*nLVar, *nRatioVar));
  RooCBShape cbRight("cbRight", "CB right tail", *massVar, mean, sigmaR, alphaR, nR);

  // --- gauss ---
  // RooRealVar sigmaGRatio("sigmaGRatio", "", 1.0, 0.5, 10.0);
  RooFormulaVar sigmaG("sigmaG", "sigma gauss", "@0*@1", RooArgList(*sigmaLVar, *sigmaGRatioVar));
  RooGaussian gaus("gaus", "gaus core", *massVar, mean, sigmaG);

  // --- combine sig components ---
  RooRealVar f_cbR("f_cbR", "frac of right CB", 0.5, 0, 1.0);
  RooRealVar f_gaus("f_gaus", "frac of gaus", 0.05, 0, 1);
  RooAddPdf sig("sig", "", RooArgList(cbLeft, cbRight, gaus), RooArgList(f_cbR, f_gaus));
  // RooAddPdf sig("sig", "",
  //             RooArgList(cbLeft, cbRight),
  //               RooArgList(f_cbR));

  // --- bkg ---
  RooRealVar c1("c1", "c1", 0.01, -1, 1);
  RooRealVar c2("c2", "c2", 0.01, -1, 1);
  RooRealVar c3("c3", "c3", 0.01, -1, 1);
  RooChebychev bkg("bkg", "2nd-order Chebychev", *massVar, RooArgList(c1, c2));

  // --- Extended yields ---
  RooRealVar Nsig("Nsig", "signal yield", 50000, 10000, 100000);
  RooRealVar Nbkg("Nbkg", "background yield",20000, 5000, 50000);

  RooAddPdf model("pdf_mass_tot", "",
                RooArgList(sig, bkg), RooArgList(Nsig, Nbkg));

  auto dh = new RooDataHist("dh", "binned dataset", *massVar, *ds_red);

  // --- perform fit ---
  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"), Extended(), SumW2Error(hasWeight), Offset(true), PrintLevel(-1), NumCPU(32), EvalBackend("legacy"));

  // === draw ===
  // --- divided canvas ---
  TCanvas c_mass("c_mass", "c_mass", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  RooPlot *massFrame = massVar->frame(Range(massMin, massMax), Title("")); // Bins(80)
  ds_red->plotOn(massFrame, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(massFrame, NormRange("fitRange"), Range("fitRange"), Name("model"));
  model.plotOn(massFrame, Components("cbLeft"), LineStyle(kDotted), LineColor(kRed), Name("cbLeft"));
  model.plotOn(massFrame, Components("cbRight"), LineStyle(kDotted), LineColor(kAzure), Name("cbRight"));
  model.plotOn(massFrame, Components("gaus"), LineStyle(kDotted), LineColor(kViolet), Name("gaus"));
  model.plotOn(massFrame, Components("bkg"), LineStyle(kDotted), LineColor(kOrange), Name("bkg"));

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
  massFrame->SetMaximum(ymax * 10.0);

  // title
  massFrame->GetYaxis()->SetTitle("Events");
  massFrame->GetXaxis()->SetTitle("");
  massFrame->Draw("e");

  // --- object legend ---
  TLegend *leg = new TLegend(0.16, 0.6, 0.49, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(massFrame->findObject("data"), "Data", "lep");
  leg->AddEntry(massFrame->findObject("model"), "Fit model", "pe");
  leg->AddEntry(massFrame->findObject("cbLeft"), "CB1", "pe");
  leg->AddEntry(massFrame->findObject("cbRight"), "CB2", "pe");
  leg->AddEntry(massFrame->findObject("gaus"), "Gauss", "pe");
  leg->AddEntry(massFrame->findObject("bkg"), "Cheby2", "pe");
  leg->Draw("same");

  // --- info latex ---
  TLatex latexInfo;
  latexInfo.SetNDC();
  latexInfo.SetTextSize(0.03);
  latexInfo.SetTextFont(42);

  latexInfo.DrawLatex(0.66, 0.89, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
  latexInfo.DrawLatex(0.66, 0.83, Form("Data, %s %s", comp.c_str(), region.c_str()));
  latexInfo.DrawLatex(0.66, 0.77, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));

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

  c_mass.SaveAs(Form("mass_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.png", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh));

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  // TFile fout(Form("roots/mass_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.root", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh), "RECREATE");
  // fitResult->Write("fitResult");
  // model.Write();
  // fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}