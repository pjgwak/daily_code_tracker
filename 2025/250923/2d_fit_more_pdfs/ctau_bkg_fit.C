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
using std::cerr;
using std::cout;
using std::string;

// === helper function ===
void fixAndReport(const RooArgList &params, const char *parName)
{
  RooRealVar *par = (RooRealVar *)params.find(parName);
  if (par)
  {
    par->setConstant(kTRUE);
    cout << "Fixed " << parName << " = " << par->getVal() << "\n";
  }
  else
  {
    cerr << "Parameter " << parName << " not found in fitResult" << "\n";
  }
}

void ctau_bkg_fit(string region = "LSB")
{
  float ptLow = 6.5, ptHigh = 9;
  float yLow = 0, yHigh = 1.6;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;
  double ctMin = -2, ctMax = 4;

  TStopwatch t;
  t.Start();
  cout << "\n=== Start ctau_bkg_fit() ===\n";

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only printfrom WARNING to FATAL

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
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  if (!dsRaw)
  {
    cout << "Error: cannot find RooDataSet\n";
    return;
  }
  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight",
                                        dsRaw, *dsRaw->get(), 0, "weight");

  // === declare cuts ===
  // --- basic cuts ---
  // acceptance
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut

  // kinematics cuts
  //  - correct? (<= cBin <) -> maybe (< cBin <=) ??
  TString kineCut = Form( // tmp: no cbin
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= 0.008 && ctau3DErr < 0.20)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  TString osCut = "(recoQQsign == 0)";

  // --- region6 cuts ---
  
  const TString cutSR = "(mass >= 3.00 && mass <= 3.20)";
  const TString cutLSB = "(mass >= 2.60 && mass <= 2.95)";
  const TString cutRSB = "(mass >= 3.21 && mass <= 3.50)";

  // --- combine cuts ---
  TString regionCut;
  if (region == "SR")
    regionCut = cutSR;
  else if (region == "LSB")
    regionCut = cutLSB;
  else if (region == "RSB")
    regionCut = cutRSB;
  else
    regionCut = "(1)";

  TString fullCut = Form("%s && %s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data(),
                         regionCut.Data());

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
  RooRealVar *ctau3D = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3D"));
  if (!ctau3D)
    cout << "Warn: There is no variable 'ctau3D'\n";
  ctau3D->setRange(ctMin, ctMax); // Change PDF's range - Need.
  ctau3D->setRange("fitRange", ctMin, ctMax);

  // Bring MC Res fit result
  // -----------------------
  
  TFile fMcRes(Form("roots/ctau_mc_res_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "open");
  RooFitResult *fitResResult = (RooFitResult *)fMcRes.Get("fitResult");
  if (!fitResResult)
  {
    std::cerr << "fitResult not found!" << std::endl;
    return;
  }
  const RooArgList &params = fitResResult->floatParsFinal();
  // RooRealVar *sigma1 = (RooRealVar *)params.find("sigma1");
  RooRealVar *r21 = (RooRealVar *)params.find("r21");
  RooRealVar *r32 = (RooRealVar *)params.find("r32");
  // fixAndReport(params, "r21");
  // fixAndReport(params, "r32");
  // fixAndReport(params, "fg1");
  // fixAndReport(params, "fg2");

  // Build Res model
  // ---------------
  RooRealVar mu0("mu0", "res mean", 0);
  RooRealVar sigma1("sigma1", "", 0.03, 0.001, 0.1);
  RooFormulaVar sigma2("sigma2", "@0*@1", RooArgList(sigma1, *r21));
  RooFormulaVar sigma3("sigma3", "@0*@1", RooArgList(sigma2, *r32));

  RooGaussModel g1("g1", "g1", *ctau3D, mu0, sigma1);
  RooGaussModel g2("g2", "g2", *ctau3D, mu0, sigma2);
  RooGaussModel g3("g3", "g3", *ctau3D, mu0, sigma3);

  // combine res model
  RooRealVar f12("f12", "sum of fg1+fg2", 0.4, 0.1, 1.);
  RooRealVar frac_g1_in_g12("frac_g1_in_g12", "relative fraction", 0.5, 0.01, 1.);

  // fg1 = f12 * frac_g1_in_g12
  // fg2 = f12 * (1-frac_g1_in_g12)
  // fg3 = 1 - f12
  RooFormulaVar fg1("fg1", "@0*@1", RooArgList(f12, frac_g1_in_g12));
  RooFormulaVar fg2("fg2", "@0*(1-@1)", RooArgList(f12, frac_g1_in_g12));
  RooAddModel resModel("resModel", "3-Gauss resolution model",
                  RooArgList(g1, g2, g3), RooArgList(fg1, fg2));

  // Build Decay model
  // -----------------

  // parameters
  RooRealVar tau_mid("tau_mid", "scale factor middle", 0.1, 0.0001, 1.0);

  // ratio parameters
  RooRealVar r_left("r_left", "tau_left / tau_mid", 5, 1, 100);
  RooRealVar r_right("r_right", "tau_right / tau_mid", 5, 1, 100);
  RooFormulaVar tau_left("tau_left", "@0*@1", RooArgList(tau_mid, r_left));
  RooFormulaVar tau_right("tau_right", "@0*@1", RooArgList(tau_mid, r_right));

  // RooRealVar tau_left("tau_left", "scale factor left", 0.2, 0.01, 1.0);
  // RooRealVar tau_mid("tau_mid", "scale factor middle", 0.1, 0.001, 10.0);
  // RooRealVar tau_right("tau_right", "scale factor right", 0.4, 0.01, 1);

  // deacy models
  RooDecay decayBkgL("decayBkgL", "Left decay",
                        *ctau3D, tau_left, resModel, RooDecay::Flipped);

  RooDecay decayBkgMid("decayBkgMid", "Center decay",
                       *ctau3D, tau_mid, resModel, RooDecay::DoubleSided);

  RooDecay decayBkgR("decayBkgR", "Right decay",
                         *ctau3D, tau_right, resModel, RooDecay::SingleSided);

  // combine
  RooRealVar f_bkg_res("f_bkg_res", "fraction left", 0.2, 0.0, 1.0);
  RooRealVar f_bkg_mid("f_bkg_mid", "fraction middle", 0.2, 0.0, 1.0);
  RooRealVar f_bkg_left("f_bkg_left", "fraction middle", 0.2, 0.0, 1.0);
  auto model = new RooAddPdf("model", "sum of 2 decays",
                             RooArgList(decayBkgMid, decayBkgL, decayBkgR),
                             RooArgList(f_bkg_mid, f_bkg_left), kTRUE);

  // --- extended fit ---
  // RooRealVar N("N", "yield", dsReduced->numEntries(), 1, 500000);
  // RooExtendPdf model("model", "extended 1-Expo", expo1, N);


  // silent
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  // Perform fit
  // -----------
  auto fitResult = model->fitTo(*dsReduced,
                                Save(), Range("fitRange"),
                                SumW2Error(hasWeight), Offset(true),
                                Strategy(2),
                                NumCPU(32), EvalBackend("legacy"),
                                RooFit::RecoverFromUndefinedRegions(1),
                                RooFit::PrintEvalErrors(-1),
                                RooFit::PrintLevel(-1));

  // === draw ctau3D ===
  // --- divided canvas ---
  TCanvas c_ctau("c_ctau", "c_ctau", 800, 800);

  // --- main plot ---
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
  pad1->SetBottomMargin(0.00001);
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();

  RooPlot *f_ctau = ctau3D->frame(Range(ctMin, ctMax), Title("")); // Bins(80)
  dsReduced->plotOn(f_ctau, DataError(RooAbsData::SumW2), Name("data"));
  model->plotOn((f_ctau), NormRange("fitRange"), Range("fitRange"), Name("model"));
  model->plotOn(f_ctau, Components("resModel"), LineStyle(kSolid), LineColor(kOrange), Name("resModel"));
  model->plotOn(f_ctau, Components("decayBkgL"), LineStyle(kDashed), LineColor(kRed + 1), Name("decayBkgL"));
  model->plotOn(f_ctau, Components("decayBkgMid"), LineStyle(kDashed), LineColor(kGreen + 1), Name("decayBkgMid"));
  model->plotOn(f_ctau, Components("decayBkgR"), LineStyle(kDashed), LineColor(kMagenta), Name("decayBkgR"));

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
  TLegend *leg = new TLegend(0.49, 0.6, 0.74, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(f_ctau->findObject("data"), "Data", "lep");
  leg->AddEntry(f_ctau->findObject("model"), "Fit Model", "pe");
  leg->AddEntry(f_ctau->findObject("resModel"), "Prompt-like bkg", "pe");
  leg->AddEntry(f_ctau->findObject("decayBkgL"), "Bkg Left", "pe");
  leg->AddEntry(f_ctau->findObject("decayBkgMid"), "Bkg Mid", "pe");
  leg->AddEntry(f_ctau->findObject("decayBkgR"), "Bkg Right", "pe");
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
  latexInfo.DrawLatex(x_start, y_start + y_step * y_stepCount++, Form("Data, %s", region.c_str()));
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
  
  // === pull pad ===
  c_ctau.cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->Draw();
  pad2->cd();

  RooHist *hpull = f_ctau->pullHist("data", "model");
  RooPlot *f_pull = ctau3D->frame(Range(ctMin, ctMax), Title(""));
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

  c_ctau.SaveAs(Form("figs/ctau_bkg_%s_pT%.1f_%.1f_y%.1f_%.1f.png", region.c_str(), ptLow, ptHigh, yLow, yHigh));

  fitResult->Print("V");
  cout << "\n chi2/ndf = " << chi2ndf << endl;

  // === save results ===
  TFile fout(Form("roots/ctau_bkg_%s_pT%.1f_%.1f_y%.1f_%.1f.root", region.c_str(), ptLow, ptHigh, yLow, yHigh), "RECREATE");
  fitResult->Write("fitResult");
  model->Write();
  fout.Close();

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}