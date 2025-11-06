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

void ctau_np_fit(string region = "Full")
{
  float ptLow = 40, ptHigh = 50;
  float yLow = 0, yHigh = 1.6;
  string comp = "CtauFull";
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;
  double ctMin = 0.00001, ctMax = 8;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_np_fit() ===\n";

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only printfrom WARNING to FATAL

  // use rootlogon
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // make output folder
  gSystem->mkdir("figs", true);
  gSystem->mkdir("roots", true);

  // read input
  TFile *fInput = TFile::Open("/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_y0_2p4_230829.root");
  // /data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_JPsi_pp_GENONLY_NonPrompt_y0_1p2_230214.root
  if (!fInput || fInput->IsZombie())
  {
    cout << "Error: cannot open input file\n";
    return;
  }

  // read dataset
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!dsRaw)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }
  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight",
                                        dsRaw, *dsRaw->get(), 0);
  // RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight",
  //                                       dsRaw, *dsRaw->get(), 0, "weight");

  // === declare cuts ===
  // --- basic cuts ---
  // acceptance
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut

  // kinematics cuts
  //  - correct? (<= cBin <) -> maybe (< cBin <=) ??
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f && ctau3Dtrue >= %.3f && ctau3Dtrue < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  // TString osCut = "(recoQQsign == 0)"; -> MC Gen No sign cut
  TString fullCut = Form("%s && %s",
                         kineCut.Data(), accCut.Data());

  // === new dataset with cuts ===
  RooDataSet *dsReduced = (RooDataSet *)dsWeight->reduce(Cut(fullCut));
  if (!dsReduced || dsReduced->numEntries() == 0)
  {
    cout << "[ERROR] reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // print out object list
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced->Print();

  // check weight
  const bool hasWeight = (dsReduced->isWeighted() || dsReduced->weightVar() != nullptr);
  if (hasWeight)
  {
    // MC's weight value = 1 -> same with no-weight case
    cout << "[Info] Using weighted dataset \n";
  }
  else
  {
    cout << "[Info] Using UN-weighted dataset \n";
  }

  // variables
  RooRealVar *ctau3D = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3Dtrue"));
  if (!ctau3D)
    cout << "Warn: There is no variable 'ctau3Dtrue'\n";
  ctau3D->setRange(ctMin, ctMax); // Change PDF's range - Need.
  ctau3D->setRange("plotRange", ctMin, ctMax);
  ctau3D->setRange("fitRange", 0.01, ctMax);

  // === fit by Expo ===
  // === resolution pdf===
  RooRealVar *meanTrue = new RooRealVar("meanTrue", "proper time mean", 0);
  RooRealVar *sigmaTrue = new RooRealVar("sigmaTrue", "proper time sigma", 0.1, 0.0000001, 0.5);
  // RooResolutionModel *trueRes = new RooGaussModel("trueRes", "proper time resolution", *ctau3D, *meanTrue, *sigmaTrue);

  auto trueRes = std::make_unique<RooTruthModel>("trueRes", "truth model", *ctau3D);

  // ratio parameter
  RooRealVar tauR1("tauR1", "base B ctau", 0.01, 0.001, 0.5, "mm");
  RooRealVar tauRatio("tauRatio", "scale factor middle", 1.5, 1, 10.0);
  RooFormulaVar tauR2("tauR2", "@0*@1", RooArgList(tauR1, tauRatio));  

  RooDecay decayTrueR1("decayTrueR1", "",
                         *ctau3D, tauR1, *trueRes, RooDecay::SingleSided);
  RooDecay decayTrueR2("decayTrueR2", "",
                          *ctau3D, tauR2, *trueRes, RooDecay::SingleSided);
  

  // combine models
  RooRealVar f_right1("f_right1", "fraction middle", 0.3, 0.0, 1.0);
  RooAddPdf model("model", "",
                      RooArgList(decayTrueR1, decayTrueR2),
                      RooArgList(f_right1), kTRUE);

  // RooRealVar N("N", "yield",100000, 1, 1000000);
  // RooExtendPdf model("model", "extended 1-Expo", model_tmp, N);

  // auto model = new RooAddPdf("model", "sum of 2 decays",
  //                            RooArgList(time_BB_left, decayTrueR2, decayTrueR1, *trueRes),
  //                            RooArgList(f_left, f_mid, f_right));

  // --- perform fit ---
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  auto fitResult = model.fitTo(*dsReduced,
                               Save(),
                              //  Range("fitRange"),
                               SumW2Error(hasWeight), Offset(true),
                               Strategy(2),
                               NumCPU(32), EvalBackend("legacy"),
                               RooFit::RecoverFromUndefinedRegions(1),
                               RooFit::PrintEvalErrors(-1),
                               RooFit::PrintLevel(-1));

  // === draw ctau3Dtrue ===
  // --- divided canvas ---
  TCanvas c_ctau("c_ctau", "c_ctau", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  RooPlot *f_ctau = ctau3D->frame(Title(""));
  dsReduced->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn(f_ctau, Name("model"));
  model.plotOn(f_ctau, Components("decayTrueR1"), LineStyle(kDotted), LineColor(kOrange), Name("decayTrueR1"));
  model.plotOn(f_ctau, Components("decayTrueR2"), LineStyle(kDotted), LineColor(kRed+1), Name("decayTrueR2"));

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
  f_ctau->SetMaximum(ymax * 1e5);

  // title
  f_ctau->GetYaxis()->SetTitle("Events");
  f_ctau->GetXaxis()->SetTitle("");
  f_ctau->Draw("e");

  // --- object legend ---
  TLegend *leg = new TLegend(0.49, 0.7, 0.75, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(f_ctau->findObject("data"), "Data", "lep");
  leg->AddEntry(f_ctau->findObject("model"), "Fit Model", "pe");
  leg->AddEntry(f_ctau->findObject("decayTrueR1"), "Decay 1", "pe");
  leg->AddEntry(f_ctau->findObject("decayTrueR2"), "Decay 2", "pe");
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
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("Data, %s %s", comp.c_str(), region.c_str()));
  if (yLow == 0)
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
  else
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));

  // --- latex parameters ---
  TLatex latexParams;
  latexParams.SetNDC();
  latexParams.SetTextSize(0.03);
  latexParams.SetTextFont(42);

  x_start = 0.71;
  y_step = -0.045;
  y_stepCount = 1;
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("N_{dimuon} = %.0f #pm %.0f", N.getVal(), N.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#tau1 = %.3f #pm %.3f", tau1.getVal(), tau1.getError()));
  // latexParams.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("#sigma_{L} = %.3f #pm %.3f", sigmaLVar->getVal(), sigmaLVar->getError()));

  // === pull pad ===
  c_ctau.cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  RooHist *hpull = f_ctau->pullHist("data", "model");
  RooPlot *f_pull = ctau3D->frame(Title(""));
  f_pull->addPlotable(hpull, "P"); // P: points only

  f_pull->GetYaxis()->SetTitle("Pull");
  f_pull->GetXaxis()->SetTitle("c#tau_{3D} [mm]");
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
  latex.DrawLatex(0.82, 0.88, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  c_ctau.SaveAs(Form("figs/ctau_true_pT%.1f_%.1f_y%.1f_%.1f.png", ptLow, ptHigh, yLow, yHigh));
  c_ctau.SaveAs(Form("figs/ctau_true_pT%.1f_%.1f_y%.1f_%.1f.pdf", ptLow, ptHigh, yLow, yHigh));

  fitResult->Print("V");
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  TFile fout(Form("roots/ctau_true_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "RECREATE");
  fitResult->Write("fitResult");
  model.Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}