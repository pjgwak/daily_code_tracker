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

void setSilent() {
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

void setVarRange(RooAbsData *data, float lmin, float lmax, float errmin, float errmax)
{
  RooRealVar *mass = (RooRealVar *)data->get()->find("mass");
  RooRealVar *ctau3D = (RooRealVar *)data->get()->find("ctau3D");
  RooRealVar *ctau3DErr = (RooRealVar *)data->get()->find("ctau3DErr");

  if (mass)
    mass->setRange(2.6, 3.5);
  if (ctau3D)
    ctau3D->setRange(lmin, lmax);
  if (ctau3DErr)
    ctau3DErr->setRange(errmin, errmax);

  return;
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

    if (fabs(avgW - 1.0) < 1e-6) {
      cout << "\033[31m" // color red
           << "   [WARNING] It seems weighting is 1 - UN-WEIGHTED dataset. Is it right?\n"
           << "\033[0m";
    }
  };

  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("dataset");
  if (isWeight) {
    cout << "[Info] data is WEIGHTED\n";

    data = new RooDataSet("data", "dataset with weight",
                              data, *data->get(), 0,
                      "weight");
    printAvgWeight(data, 100);
  } else
  {
    cout << "[Info] data is UN-WEIGHTED" << endl;
    data->SetName("data");
  }
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

void getCtErrRange(RooDataSet *data, TCut &t_reduceDS_woCtErr, RooFitResult *fitMass, float lmin, float lmax, float *errmin, float *errmax)
{
  RooWorkspace *ws = new RooWorkspace("ctauerrorcheckWS");
  RooDataSet *redDataCut = (RooDataSet *)data->reduce(t_reduceDS_woCtErr);
  ws->import(*redDataCut);

  // ws->var("mass")->setRange(2.6, 3.5);
  // ws->var("mass")->setBins(45);
  // ws->var("ctau3D")->setRange(-lmin, lmax);
  // ws->var("ctau3D")->setRange(lmin, lmax);
  ws->var("ctau3DErr")->setRange(0.0, 0.992);
  ws->var("ctau3DErr")->setBins(124);

  cout << " -*-*-*-*-*-*-*-*-*-*-*-*-*- getCtErrRange -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  " << endl;
  RooDataSet *redDataSB = (RooDataSet *)redDataCut->reduce("mass < 2.9 || mass > 3.3");
  RooDataSet *redDataSIG = (RooDataSet *)redDataCut->reduce("mass > 2.9 && mass < 3.3");

  //// *** RooDataHist
  RooDataHist *tbinDataCtErrSB = new RooDataHist("tbinDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *tbinDataCtErrSIG = new RooDataHist("tbinDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  //// ****  scaleF to scale down ct err dist in 2.9-3.3 GeV/c2
  float bc = 0;
  bc = getFitResultVal(fitMass, "coefExp");

  // if (!inOpt.mBkgFunct.compare("expBkg"))
  // else if (!inOpt.mBkgFunct.compare("polBkg"))
  //   bc = ws->var("coefPol")->getVal();
  float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  //// *** (tbinSubtractedSIG) = (tbinDataCtErrSIG) - scaleF*(tbinDataCtErrSB)
  RooDataHist *tbinSubtractedSIG = new RooDataHist("tbinSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  RooDataHist *tbinScaledBKG = subtractSidebands(ws, tbinSubtractedSIG, tbinDataCtErrSIG, tbinDataCtErrSB, scaleF, "ctau3DErr");

  //// **** Check the minimum and maximum of the ctau error in signal and background regions
  TH1 *histDataCtErrSIG = tbinDataCtErrSIG->createHistogram("histDataCtErrSIG", *ws->var("ctau3DErr"));
  TH1 *histSubtractedSIG = tbinSubtractedSIG->createHistogram("histSubtractedSIG", *ws->var("ctau3DErr"));
  TH1 *histScaledBKG = tbinScaledBKG->createHistogram("histScaledBKG", *ws->var("ctau3DErr"));

  double minSig = 0.5, maxSig = 0.0, minBkg = 0.5, maxBkg = 0.0;
  double cutValue = 0.2;

  int maxBinSig = histSubtractedSIG->GetMaximumBin();
  int maxBinBkg = histScaledBKG->GetMaximumBin();

  minSig = histSubtractedSIG->GetBinLowEdge(maxBinSig);
  minBkg = histScaledBKG->GetBinLowEdge(maxBinBkg);
  // pick up lower bound next to other non-zero bins
  for (int xbins = maxBinSig; xbins > 0; xbins--)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      minSig = histSubtractedSIG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG low edge: " << histSubtractedSIG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins > 0; xbins--)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      minBkg = histScaledBKG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG low edge: " << histScaledBKG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }

  // pick up upper bound next to other non-zero bins
  maxSig = histSubtractedSIG->GetBinLowEdge(maxBinSig + 1);
  maxBkg = histScaledBKG->GetBinLowEdge(maxBinBkg + 1);
  for (int xbins = maxBinSig; xbins < histSubtractedSIG->GetNbinsX(); xbins++)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      maxSig = histSubtractedSIG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG upper edge: " << histSubtractedSIG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinSig; xbins < histScaledBKG->GetNbinsX(); xbins++)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      maxBkg = histScaledBKG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG upper edge: " << histScaledBKG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }

  // choose the higher lower limit, lower upper limit
  double tmpMin = 0, tmpMax = 0;
  if (minSig > minBkg)
    tmpMin = minSig;
  else
    tmpMin = minBkg;
  if (maxSig < maxBkg)
    tmpMax = maxSig;
  else
    tmpMax = maxBkg;

  // round off lower limit -> allow more entries on the lower limits
  tmpMin = TMath::Floor(tmpMin * 1000);
  tmpMin = tmpMin / (double)1000.0;

  // round up upper limit -> allow more entries on the upper limits
  tmpMax = TMath::Ceil(tmpMax * 1000);
  tmpMax = tmpMax / (double)1000.0;

  char reduceDS[512];
  sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
  if (redDataTmp->sumEntries() < redDataCut->sumEntries() * 0.9)
  { // if ctau error range cuts off >10% events
    delete redDataTmp;
    sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
    tmpMin = minSig;
    tmpMax = maxSig;
  }
  if ((tmpMax - tmpMin) < 0.008)
  {
    cout << "getCtErrRange:: Maximum is less than minimum! Possibly there are few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  //// *** draw final ctau error plot
  TCanvas c0("ctau_err", "ctau_err", 500, 500);
  c0.Draw();
  c0.cd();
  c0.SetLogy(1);
  RooPlot *errframe2 = ws->var("ctau3DErr")->frame();
  // tbinDataCtErrSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kRed),LineColor(kRed));
  // tbinDataCtErrSB->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kGreen+2),LineColor(kGreen+2),MarkerStyle(24));
  // tbinScaledBKG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kBlue),MarkerStyle(24),LineColor(kBlue));
  // tbinSubtractedSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),LineColor(kWhite));
  const double max = errframe2->GetMaximum() * 1.3;
  errframe2->SetMaximum(max);
  errframe2->SetMinimum(0.2);
  errframe2->Draw();
  histDataCtErrSIG->SetMarkerColor(kRed);
  histDataCtErrSIG->SetLineColor(kWhite);
  histDataCtErrSIG->SetMarkerStyle(24);
  histDataCtErrSIG->GetXaxis()->CenterTitle(1);
  histDataCtErrSIG->GetYaxis()->CenterTitle(1);
  histDataCtErrSIG->Draw("pe");
  histScaledBKG->SetMarkerColor(kBlue);
  histScaledBKG->SetLineColor(kWhite);
  histScaledBKG->SetMarkerStyle(24);
  histScaledBKG->Draw("pe same");
  histSubtractedSIG->SetLineColor(kWhite);
  histSubtractedSIG->Draw("pe same");

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetTextSize(0.035);

  // t->DrawLatex(0.92, 0.84, opt.rapString);
  // t->DrawLatex(0.92, 0.78, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.92, 0.72, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.92, 0.72, opt.etString);

  char comment[200];
  sprintf(comment, "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
  t->SetTextSize(0.04);
  t->SetTextColor(kRed);
  t->DrawLatex(0.92, 0.6, comment);
  t->SetTextColor(kBlack);

  TLegend legsb(0.6, 0.19, 0.9, 0.35, NULL, "brNDC");
  legsb.SetFillStyle(0);
  legsb.SetBorderSize(0);
  legsb.SetShadowColor(0);
  legsb.SetMargin(0.2);
  legsb.SetTextFont(42);
  legsb.SetTextSize(0.035);
  legsb.AddEntry(histDataCtErrSIG, "sig cands", "p");
  legsb.AddEntry(histScaledBKG, "scaled bkg", "p");
  legsb.AddEntry(histSubtractedSIG, "sig (= cands - bkg)", "p");
  legsb.Draw("same");

  c0.SaveAs(Form("figs/ctau_err_raw_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.png", ptLow, ptHigh, yLow, yHigh, isWeight));

  *errmin = tmpMin;
  *errmax = tmpMax;
  //  cout << "getCtErrRange:: " << t_reduceDS_woCtErr << " " << lmin << " " << lmax << " " << *errmin << " " << *errmax << endl;

  delete ws;
  delete redDataCut;
  delete redDataTmp;
  // delete binData;
  // delete binDataCtErr;
  // delete binDataSB;
  delete tbinDataCtErrSB;
  delete tbinDataCtErrSIG;
  delete tbinSubtractedSIG;
  delete tbinScaledBKG;
  delete t;
}

void reduce_data_with_err_cut()
{
  cout << "\n=== Start reduce_data_with_err_cut() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  
  // === workspace ===
  auto ws = new RooWorkspace("ws");

  // === read inputs ===
  cout << "=== Import inputs ===\n";
  RooDataSet *data = nullptr;
  const bool isWeight = true;
  readInputs(data, isWeight);

  // === get mass fit result ===
  cout << "\n=== read mass fit result ===\n";
  string fileNameMassRaw = Form("roots/mass_raw_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInMassRaw(fileNameMassRaw.c_str());
  if (fInMassRaw.IsZombie())
  {
    cerr << "[Error] Failed to load: " << fileNameMassRaw.c_str() << "\n";
    return;
  }
  cout << "Load: " << fileNameMassRaw.c_str() << "\n";
  auto fitMass = (RooFitResult *)fInMassRaw.Get("fitMass");

  // === reduce dataset ===
  cout << "\n=== reduce dataset ===\n";

  // --- make cuts ---
  // muon kinematic cut
  TCut muonCut =
      "((abs(eta1) <= 1.2 && pt1 >= 3.5) || "
      "(abs(eta2) <= 1.2 && pt2 >= 3.5) || "
      "((abs(eta1) > 1.2 && abs(eta1) <= 2.1) && pt1 >= 5.47-1.89*abs(eta1)) || "
      "((abs(eta2) > 1.2 && abs(eta2) <= 2.1) && pt2 >= 5.47-1.89*abs(eta2)) || "
      "((abs(eta1) > 2.1 && abs(eta1) <= 2.4) && pt1 >= 1.5) || "
      "((abs(eta2) > 2.1 && abs(eta2) <= 2.4) && pt2 >= 1.5))";

  // Opposite sign
  TCut pairCut = "recoQQsign == 0";

  // kinematic bins cuts
  TCut kineCut = Form("(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f)",
                      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

  // ctau3D 
  TCut ctCut = Form("(ctau3D >= %.3f && ctau3D < %.3f)", ctMin, ctMax);

  // MC true ctau
  TCut ctTrueCut = "(ctau3Dtrue >= 0.0001 && ctau3Dtrue < 9)";

  TCut reduceDS_woCtErr = muonCut && pairCut && kineCut && ctCut; // data without ctauErr cut

  // decide ctauErr range
  getCtErrRange(data, reduceDS_woCtErr, fitMass, ctMin, ctMax, &errmin, &errmax);

  // add ctau3DErr cuts - Must excute after getCtErrRange()
  TCut errCut = Form("(ctau3DErr >= %.3f && ctau3DErr < %.3f)", errmin, errmax);

  // total cuts
  TCut reduceDS = muonCut && pairCut && kineCut && ctCut && errCut; // data with ctauErr cut

  // --- reduce ---
  // 필요한 변수만 리턴해서 저장? -> 가장 큰 파일도 650 MB 정도라서 용량에 여유가 있다.
  // -> 코드 복잡도를 줄이기 위해 그냥 통째로 저장
  RooDataSet *redData;
  redData = (RooDataSet *)data->reduce(reduceDS);
  redData->SetName("redData");

  // --- sideband dataset ===
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  redDataSIG->SetName("redDataSIG");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  redDataSB->SetName("redDataSB");


  // === set variable range ===
  setVarRange(redData, ctMin, ctMax, errmin, errmax);
  setVarRange(redDataSIG, ctMin, ctMax, errmin, errmax);
  setVarRange(redDataSB, ctMin, ctMax, errmin, errmax);


  // === save results ===
  TFile fileOut(Form("reduced_data/data_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  redData->Write();
  redDataSIG->Write();
  redDataSB->Write();
  fileOut.Close();

  cout << "=== Finish reduce_data_with_err_cut() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}