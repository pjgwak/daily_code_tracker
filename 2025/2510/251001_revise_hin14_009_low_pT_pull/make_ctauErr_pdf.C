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

void readInputs(RooDataSet *&data, RooDataSet *&redDataSB, RooDataSet *&redDataSIG, const bool isWeight)
{
  // --- reducded Data ---
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

  redDataSB = (RooDataSet *)fInData.Get("redDataSB");
  redDataSB->SetName("redDataSB");

  redDataSIG = (RooDataSet *)fInData.Get("redDataSIG");
  redDataSIG->SetName("redDataSIG");
}

float getFitResultVal(RooFitResult *fitResult, const std::string &varName)
{
  if (!fitResult)
  {
    cerr << "[Error] getFitResultVal: null fitResult" << "\n";
    return 0.0;
  }

  RooRealVar *varFit = (RooRealVar *)fitResult->floatParsFinal().find(varName.c_str());
  if (!varFit)
  {
    cerr << "[Warning] variable '" << varName << "' not found in fitResult" << "\n";
    return 0.0;
  }

  return varFit->getVal();
}

RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr")
{

  if (binSIG->numEntries() != binSB->numEntries())
  {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }
  RooDataHist *binScaleBKG = new RooDataHist("binScaleBKG", "scaled SB", RooArgSet(*(ws->var(varName.c_str()))));

  //// **** bin-by-bin scaling
  const RooArgSet *argSIG;
  const RooArgSet *argSB;
  for (Int_t i = 0; i < binSIG->numEntries(); i++)
  {
    argSIG = binSIG->get(i);
    argSB = binSB->get(i);
    RooRealVar *thisVar = (RooRealVar *)argSIG->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    //// *** set minimum as 0.1 to prevent empty PDF
    float wBkg = binSB->weight(*argSB, 0, false);
    if (wBkg <= 0.1)
      wBkg = 0.1;
    binScaleBKG->add(RooArgSet(*(ws->var(varName.c_str()))), wBkg);
    float newWeight = binSIG->weight(*argSIG, 0, false) - scalefactor * binSB->weight(*argSB, 0, false);
    if (newWeight <= 0.1)
      newWeight = 0.1;
    binSubtrSIG->add(RooArgSet(*(ws->var(varName.c_str()))), newWeight);
  }
  return binScaleBKG;
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

void drawCtauErrPdf(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG)
{
  RooPlot *errframe = ws->var("ctau3DErr")->frame();

  binDataCtErrSB->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kBlue), MarkerColor(kBlue), MarkerStyle(kOpenCircle));
  ws->pdf("errPdfBkg")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binDataCtErrSB->sumEntries(), RooAbsReal::NumEvent));
  // ws->pdf("errPdfBkg")->plotOn(errframe);

  TCanvas c0;
  string titlestr;
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = Form("figs/ctau_err_bkg_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.png", ptLow, ptHigh, yLow, yHigh, isWeight);
  c0.SaveAs(titlestr.c_str());
  delete errframe;

  errframe = ws->var("ctau3DErr")->frame();
  binSubtractedSIG->plotOn(errframe, DataError(RooAbsData::SumW2), DrawOption("F"), FillColor(kWhite), LineColor(kWhite)); // just for axis
  ws->pdf("errPdfSig")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binSubtractedSIG->sumEntries(), RooAbsReal::NumEvent));
  binDataCtErrSIG->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kRed), MarkerColor(kRed), MarkerStyle(kOpenCircle));
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = Form("figs/ctau_err_sig_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.png", ptLow, ptHigh, yLow, yHigh, isWeight);
  c0.SaveAs(titlestr.c_str());
  delete errframe;
}

void make_ctauErr_pdf()
{
  cout << "\n=== Start make_ctauErr_pdf() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // === read inputs ===
  cout << "=== Import inputs ===\n";
  RooDataSet *redData = nullptr, *redDataSB = nullptr, *redDataSIG = nullptr;
  readInputs(redData, redDataSB, redDataSIG, isWeight);

  // === get mass fit result ===
  cout << "\n=== read mass fit result ===\n";
  string fileNameMass = Form("roots/mass_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInMass(fileNameMass.c_str());
  if (fInMass.IsZombie())
  {
    cerr << "[Error] Failed to load: " << fileNameMass.c_str() << "\n";
    return;
  }
  cout << "Load: " << fileNameMass.c_str() << "\n";
  auto fitMass = (RooFitResult *)fInMass.Get("fitMass");

  // make workspace
  auto ws = new RooWorkspace("ws");
  ws->import(*redData);

  // scaleF to scale down ctErr distribution in 2.9-3.3 GeV/c2
  float bc;
  bc = bc = getFitResultVal(fitMass, "coefExp");
  float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  // RooDataSet(unbinned) to RooDataHist (binned)
  RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr", "binDataCtErr", RooArgSet(*(ws->var("ctau3DErr"))), *redData);
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  // extract Sig and Bkg: (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  RooDataHist *binSubtractedSIG, *binScaledBKG;
  binSubtractedSIG = new RooDataHist("binSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  binScaledBKG = subtractSidebands(ws, binSubtractedSIG, binDataCtErrSIG, binDataCtErrSB, scaleF, "ctau3DErr");

  // error PDF
  RooHistPdf errPdfSig("errPdfSig", "Error PDF signal", RooArgSet(*(ws->var("ctau3DErr"))), *binSubtractedSIG);
  ws->import(errPdfSig);
  RooHistPdf errPdfBkg("errPdfBkg", "Error PDF bkg scaled", RooArgSet(*(ws->var("ctau3DErr"))), *binScaledBKG);
  ws->import(errPdfBkg);

  //// **** Draw CtErr PDF
  drawCtauErrPdf(ws, binDataCtErrSB, binDataCtErrSIG, binSubtractedSIG, binScaledBKG);


  // === save results ===
  cout << "\n=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/ctau_err_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  errPdfSig.Write();
  errPdfBkg.Write();
  fileOut.Close();

  // fitMass->Print("v");

  cout << "\n=== Finish make_ctauErr_pdf() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
