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

void readInputs(RooDataSet *&data, const bool isWeight)
{
  // calculate average weight (first 100 evets)
  auto printAvgWeight = [](RooDataSet *ds, int nCheck = 100)
  {
    // if (!ds || !ds->weightVar())
    // {
    //   cout << "   [Check] No weight variable found" << endl;
    //   return;
    // }
    int nEntries = ds->numEntries();
    int nLoop = std::min(nCheck, nEntries);

    double sumW = 0.0;
    for (int i = 0; i < nLoop; i++)
    {
      ds->get(i); // read event
      sumW += ds->weight();
    }
    double avgW = (nLoop > 0) ? sumW / nLoop : 0.0;
    cout << "   [Check] Avg weight of first " << nLoop << " events = " << avgW << endl;

    if (fabs(avgW - 1.0) < 1e-6)
    {
      cout << "\033[31m" // color red
           << "   [WARNING] It seems weighting is 1 - UN-WEIGHTED dataset. Is it right?\n"
           << "\033[0m";
    }
  };

  // --- Data ---
  string fileNameData = Form("reduced_data/data_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("redData");
  if (isWeight)
  {
    cout << "[Info] data is WEIGHTED\n";
  }
  else 
    cout << "[Info] data is UN-WEIGHTED" << endl;
  data->SetName("redData");
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

void defineMassBkg(RooWorkspace *ws)
{
  // 1st order polynomial
  ws->factory("Polynomial::polBkg(mass,{coefPol[-0.05,-5.,5.]})");

  // expo
  ws->factory("Exponential::expBkg(mass,coefExp[-1.,-3.,2.])");
  return;
}

void fixParamsMC(RooWorkspace *ws, RooFitResult *fitResult, const vector<string> &fixList)
{
  if (!ws || !fitResult)
  {
    cerr << "[Error] fixParamsMC: null input" << "\n";
    return;
  }

  const RooArgList &pars = fitResult->floatParsFinal();

  for (const auto &name : fixList)
  {
    RooRealVar *varWS = ws->var(name.c_str());
    if (!varWS)
    {
      cerr << "[Warning] variable '" << name << "' not found in workspace" << "\n";
      continue;
    }

    RooRealVar *varFit = dynamic_cast<RooRealVar *>(pars.find(name.c_str()));
    if (!varFit)
    {
      cerr << "[Warning] variable '" << name << "' not found in fitResult" << "\n";
      continue;
    }

    // set value and fix
    varWS->setVal(varFit->getVal());
    varWS->setConstant(true);

    std::cout << "[Info] variable '" << name
              << "' set to " << varFit->getVal()
              << " and fixed." << "\n";
  }
}

void drawInclusiveMassPlots(RooWorkspace *ws, RooDataSet *data, RooFitResult *fitRes)
{
  RooRealVar *mass = ws->var("mass");

  RooPlot *frame = mass->frame(Title("Inclusive MC Mass Fit"));

  data->plotOn(frame, DataError(RooAbsData::SumW2), Name("dataHist"));
  ws->pdf("MassPDF")->plotOn(frame, LineColor(kRed), Name("fitCurve"));

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
  c.SaveAs(Form("figs/mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.png", ptLow, ptHigh, yLow, yHigh, isWeight));
}

void mass_fit_redData()
{
  cout << "\n=== Start mass_fit_redData() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // === read inputs ===
  cout << "=== Import inputs ===\n";
  RooDataSet *redData = nullptr;
  readInputs(redData, isWeight);

  // === get MC mass fit result ===
  cout << "\n=== read MC mass fit result ===\n";
  string fileNamePrMc = Form("roots/mc_mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInPRMC(fileNamePrMc.c_str());
  if (fInPRMC.IsZombie())
  {
    cerr << "[Error] Failed to load: " << fileNamePrMc.c_str() << "\n";
    return;
  }
  cout << "Load: " << fileNamePrMc.c_str() << "\n";
  auto fitMcMass = (RooFitResult *)fInPRMC.Get("fitresult_MassMCPDF_redPRMC");

  // make workspace
  auto ws = new RooWorkspace("ws");
  ws->import(*redData);
  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setRange("fitRange", 2.6, 3.5);

  // === build mass model ===
  cout << "\n=== build mass model ===\n";
  defineMassSig(ws);
  defineMassBkg(ws);

  char funct[100];
  sprintf(funct, "SUM::MassPDF(NSig[%f,1.0,5000000.0]*%s, NBkg[%f,1.0,50000000.0]*%s)", 1000000., "G1CB1Sig", 1000000., "expBkg"); // factory 문법으로 바꾸기
  ws->factory(funct);

  // fix MC parameters
  vector<string> fixList = {"alpha", "enne", "sigmaSig2", "alpha1", "alpha2", "enne1", "enne2"};
  fixParamsMC(ws, fitMcMass, fixList);

  // === fitting ===
  cout << "\n=== perform fit ===\n";
  RooFitResult *fitMass = ws->pdf("MassPDF")->fitTo(*redData, Extended(1), Minos(0), Save(1), PrintLevel(0), SumW2Error(kTRUE), NumCPU(32)); //, EvalBackend("legacy")

  // === draw plot ===
  cout << "\n=== draw plot ===\n";
  drawInclusiveMassPlots(ws, redData, fitMass);

  // // // === save results ===
  cout << "\n=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  fitMass->Write("fitMass");
  fileOut.Close();

  fitMass->Print("v");

  cout << "\n=== Finish mass_fit_redData() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
