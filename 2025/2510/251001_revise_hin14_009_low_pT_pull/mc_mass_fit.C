#include <TStopwatch.h>
#include <RooArgList.h>

using namespace RooFit;

// === kinematics === -> Read from config!
float ptLow = 13, ptHigh = 15;
float yLow = 0, yHigh = 2.4;
float massLow = 2.6, massHigh = 3.5;
int cLow = 0, cHigh = 180;
double ctMin = -2, ctMax = 6; // lmin, lmax: 1 for lowpT, 2 for higpT
float errmin = 0.008, errmax = 0.3;
const bool isWeight = true;

void setSilent()
{
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only print from WARNING to FATAL
}

void defineMassSig(RooWorkspace *ws)
{
  // --- Gauss + CB ---
  // ws->factory("Gaussian::G1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.008,0.075])");
  // ws->factory("CBShape::CB1Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha[1.9,1.2,2.8],enne[2.5,1.0,4.0])");
  // ws->factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*G1Sig,CB1Sig)");

  // --- DCB ---
  ws->factory("CBShape::CB1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.0008,0.075],alpha1[1.9,1.2,2.8],enne1[2.5,1.0,4.0])");
  ws->factory("CBShape::CB2Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha2[1,0.1,5],enne2[2.5,1.0,500.0])");

  ws->factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*CB1Sig,CB2Sig)");
  return;
}

void drawInclusiveMcMassPlots(RooWorkspace *ws,
                              RooDataSet *data,
                              RooFitResult *fitRes,
                              int isMC)
{
  RooRealVar *mass = ws->var("mass");

  RooPlot *frame = mass->frame(Title("Inclusive MC Mass Fit"));

  data->plotOn(frame, DataError(RooAbsData::SumW2), Name("dataHist"));
  ws->pdf("MassMCPDF")->plotOn(frame, LineColor(kRed), Name("fitCurve"));


  // chi2
  int nPar = fitRes->floatParsFinal().getSize(); 
  int nBins = frame->GetNbinsX();
  double chi2ndf = frame->chiSquare("fitCurve", "dataHist", nPar);
  double chi2 = chi2ndf * (nBins - nPar);
  int ndf = nBins - nPar;

  // pull 
  RooHist *hpull = frame->pullHist("dataHist", "fitCurve");
  RooPlot *framePull = mass->frame(Title("Pull"));
  framePull->addPlotable(hpull, "P");

  TCanvas c("c", "c", 800, 800);
  c.Divide(1, 2);
  c.cd(1);
  frame->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13);
  char buf[128];
  sprintf(buf, "#chi^{2}/ndf = %.2f", chi2ndf);
  latex.DrawLatex(0.75, 0.85, buf);

  c.cd(2);
  framePull->Draw();

  gSystem->mkdir("figs", true);
  c.SaveAs(Form("figs/mc_mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.png", ptLow, ptHigh, yLow, yHigh, isWeight));
}

void mc_mass_fit()
{
  cout << "\n=== Start mc_mass_fit() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // read input
  string fileNamePrMc = Form("reduced_data/mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInPRMC(fileNamePrMc.c_str());
  cout << "Load: " << fileNamePrMc.c_str() << endl;
  auto redPRMC = (RooDataSet *)fInPRMC.Get("redPRMC");

  // make workspace
  auto ws = new RooWorkspace("ws");
  ws->import(*redPRMC);
  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setRange("fitRange", 2.6, 3.5);

  // build model
  defineMassSig(ws);
  char funct[100];
  sprintf(funct, "SUM::MassMCPDF(NSig[%f,1.0,100000.0]*%s)", 50000., "G1CB1Sig"); // factory 문법으로
  ws->factory(funct);

  // fitting
  RooFitResult *fitMassMC = ws->pdf("MassMCPDF")->fitTo(*redPRMC, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), PrintLevel(0), NumCPU(32));

  // draw plot
  drawInclusiveMcMassPlots(ws, redPRMC, fitMassMC, /*isMC*/1);

  // fix parameters?
  // ws->var("alpha")->setConstant(kFALSE);
  // ws->var("enne")->setVal(2.1); // fix??? -> Do we really want to fix it as 2.1 instad of MC result?
  // ws->var("enne")->setConstant(kTRUE);
  // ws->var("coefPol")->setRange(-5., 5.);
  // ws->var("coefPol")->setVal(-0.05);
  // ws->var("coefPol")->setConstant(kFALSE);
  
  // // === save results ===
  cout << "=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/mc_mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  ws->Write("ws_mcMass");
  fitMassMC->Write();
  fileOut.Close();

  fitMassMC->Print("v");

  cout << "=== Finish mc_mass_fit() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}

