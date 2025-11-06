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
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  string comp = "CtauFull";
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_np_fit() ===\n";
  
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
  // TString osCut = "(recoQQsign == 0)";
  TString osCut = "(1)";

  // --- region6 cuts ---
  double ctMin = 0.0001, ctMax = 8;
  const TString cutCtauFull = Form("(ctau3Dtrue >= %.2f && ctau3Dtrue <= %.2f)", ctMin, ctMax);

  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";

  // --- combine cuts ---
  TString compCut;
  compCut = cutCtauFull;

  TString regionCut;
  regionCut = "(1)";

  // TString fullCut = Form("%s && %s && %s && %s && %s",
  //                        osCut.Data(),
  //                        accCut.Data(),
  //                        kineCut.Data(),
  //                        compCut.Data(),
  //                        regionCut.Data());

  TString fullCut = Form("%s && %s && %s",
                         kineCut.Data(), accCut.Data(), cutCtauFull.Data());

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
  RooRealVar *ctau3D = dynamic_cast<RooRealVar *>(ds_red->get()->find("ctau3Dtrue"));
  if (!ctau3D)
    cout << "Warn: There is no variable 'ctau3Dtrue'\n";
  ctau3D->setRange(ctMin, ctMax); // Change PDF's range - Need.
  ctau3D->setRange("plotRange", ctMin, ctMax);
  ctau3D->setRange("fitRange", 0.0001, ctMax);

  // === fit by Expo ===
  // === resolution pdf===
  RooRealVar *Mean_tom = new RooRealVar("tomMean", "proper time mean", 0);
  RooRealVar *ratio_tom = new RooRealVar("ratio_tom", "ratio", 1.2, 1, 10); // ratio of what?
  RooRealVar *Sigma0_tom = new RooRealVar("Sigma0_tom", "proper time sigma", 0.1, 0.001, 0.2);
  // per event err??
  // RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *ctau3D_err));
  RooFormulaVar *Sigma1_tom = new RooFormulaVar("Sigma1_tom", "@0*@1", RooArgList(*ratio_tom, *Sigma0_tom));

  RooRealVar *frac_g_tom = new RooRealVar("f_g_#sigma_{t0}", "gauss resol. fracion", 0.8048, 0.01, 1);

  RooResolutionModel *g0_tom = new RooGaussModel("g0_tom", "proper time resolution", *ctau3D, *Mean_tom, *Sigma0_tom);
  RooResolutionModel *g1_tom = new RooGaussModel("g1_tom", "proper time resolution", *ctau3D, *Mean_tom, *Sigma1_tom);

  // one gauss signal
  // RooResolutionModel *properTimeRes = new RooGaussModel("properTimeRes", "proper time resolution", *ctau3D, *Mean_tom, *Sigma1_tom);

  auto properTimeRes = new RooTruthModel("properTimeRes", "delta function", *ctau3D);

  // two gauss signal
  // RooResolutionModel *properTimeRes = new RooAddModel("properTimeRes", "dual gaussain PDF", RooArgList(*g0_tom, *g1_tom), RooArgList(*frac_g_tom));

  // === three decay pdf ===
  RooRealVar tau_right("tau_right", "base B ctau", 0.01, 0.001, 0.5, "mm");

  // RooRealVar tau_left("tau_left", "scale factor left", 1.0, 0.001, 100.0);
  // RooRealVar tau_right("tau_right", "scale factor right", 3.0, 0.1, 100.0);
  // RooRealVar tau_right2("tau_right2", "scale factor middle", 3, 0.01, 10.0);
  

  // ratio parameter
  // RooRealVar r_left("r_left", "scale factor left", 1.0, 0.01, 100.0);
  RooRealVar r_right2("r_right2", "scale factor middle", 3, 1, 100.0);
  // RooRealVar r_right("r_right", "scale factor right", 3.0, 0.1, 100.0);

  // // tau for bkg pdfs
  // RooFormulaVar tau_left("tau_left", "@0*@1", RooArgList(tau0, r_left));
  RooFormulaVar tau_right2("tau_right2", "@0*@1", RooArgList(tau_right, r_right2));
  // RooFormulaVar tau_right("tau_right", "@0*@1", RooArgList(tau0, r_right));

  // --- 3 decay PDF
  // RooDecay time_BB_left("time_BB_left", "Left decay",
  //                       *ctau3D, tau_left, *properTimeRes, RooDecay::Flipped);

  

  RooDecay time_BB_right("time_BB_right", "Right decay",
                         *ctau3D, tau_right, *properTimeRes, RooDecay::SingleSided);
  RooDecay time_BB_right2("time_BB_right2", "Center decay",
                          *ctau3D, tau_right2, *properTimeRes, RooDecay::SingleSided);

  // RooDecay model_tmp("model_tmp", "Right decay",
  //                        *ctau3D, tau_right, *properTimeRes, RooDecay::SingleSided);

  RooRealVar f_left("f_left", "fraction left", 0.3, 0.0, 1.0);
  RooRealVar f_mid("f_mid", "fraction middle", 0.3, 0.0, 1.0);
  RooRealVar f_right2("f_right2", "fraction middle", 0.3, 0.0, 1.0);

  RooAddPdf model_tmp("model_tmp", "sum of 2 decays",
                                 RooArgList(time_BB_right2, time_BB_right),
                                 RooArgList(f_right2), kTRUE);

  RooRealVar N("N", "yield",100000, 1, 50000000);
  RooExtendPdf model("model", "extended 1-Expo", model_tmp, N);

  // auto model = new RooAddPdf("model", "sum of 2 decays",
  //                            RooArgList(time_BB_left, time_BB_right2, time_BB_right, *properTimeRes),
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

  auto fitResult = model.fitTo(*ds_red, Save(), Range("fitRange"), Offset(true),
                                NumCPU(32),
                                EvalBackend("legacy"),
                                RooFit::RecoverFromUndefinedRegions(2.), // This is how RooFit behaved prior to ROOT 6.24
                                RooFit::PrintEvalErrors(-1),             // We are expecting a lot of evaluation errors. -1 switches off printing.
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

  RooPlot *f_ctau = ctau3D->frame(Range("fitRange"), Title("")); // Bins(80)
  ds_red->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"));
  model.plotOn((f_ctau), NormRange("fitRange"), Range("fitRange"), Name("model"));

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
  f_ctau->SetMaximum(ymax * 10.0);

  // title
  f_ctau->GetYaxis()->SetTitle("Events");
  f_ctau->GetXaxis()->SetTitle("");
  f_ctau->Draw("e");

  // --- object legend ---
  TLegend *leg = new TLegend(0.49, 0.8, 0.70, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(f_ctau->findObject("data"), "Data", "lep");
  leg->AddEntry(f_ctau->findObject("model"), "Expo", "pe");
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
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
  else
    latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));

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
  RooPlot *f_pull = ctau3D->frame(Range("fitRange"), Title(""));
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

  c_ctau.SaveAs(Form("figs/ctau_NPGen_%s_%s_pT%.1f_%.1f_y%.1f_%.1f.png", comp.c_str(), region.c_str(), ptLow, ptHigh, yLow, yHigh));

  fitResult->Print();
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  TFile fout(Form("roots/ctau_NPGen_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "RECREATE");
  fitResult->Write("fitResult");
  model.Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}