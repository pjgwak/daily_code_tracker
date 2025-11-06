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

  // --- PRMC ---
  string fileNameData = Form("reduced_data/mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("redPRMC");
  if (isWeight)
  {
    cout << "[Info] data is WEIGHTED\n";
  }
  else
    cout << "[Info] data is UN-WEIGHTED" << endl;
  data->SetName("redPRMC");
}

void readErrPdf(RooHistPdf *&sigPdf, TFile *&fInData)
{
  string fileNameData = Form("roots/ctau_err_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  fInData = TFile::Open(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  sigPdf = dynamic_cast<RooHistPdf *>(fInData->Get("errPdfSig"));
}

void defineCtPRRes(RooWorkspace *ws)
{
  ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 10],one[1.0],ctau3DErr)");
  ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.01, 3],one,ctau3DErr)");
  ws->factory("AddModel::CtPRRes({GW_PRRes,GN_PRRes},{fracRes[0.5,0.01,0.999]})");

  // ws->factory("GaussModel::GW_PRRes(ctau3D,meanPRResW[0],sigmaPRResW[2.3, 0.001, 10],one[1.0],ctau3DErr)");
  // ws->factory("GaussModel::GN_PRRes(ctau3D,meanPRResN[0],sigmaPRResN[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("GaussModel::GVW_PRRes(ctau3D,meanPRResVW[0],sigmaPRResVW[0.8,0.001, 10],one,ctau3DErr)");
  // ws->factory("AddModel::CtPRRes({GN_PRRes, GW_PRRes, GVW_PRRes},{fracResN[0.01,0.001,0.999], fracResW[0.01,0.001,0.999]})");
  return;
}

void drawCtauResolPlots(RooWorkspace *ws, bool fitMC, RooPlot *tframePR, RooFitResult *fitRes)
{
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextSize(0.035);

  TCanvas c1("c1", "c1", 600, 700);
  c1.Divide(1, 2);
  c1.cd(1);
  gPad->SetLogy();
  gPad->SetPad(0.0, 0.3, 1.0, 1.0);
  gPad->SetBottomMargin(0.02);

  // ----- Main plot -----
  tframePR->Draw();

  char buf[256];
  sprintf(buf, "#sigma(G_{N}): %.2f", ws->var("sigmaPRResN")->getVal());
  t->DrawLatex(0.55, 0.31, buf);

  sprintf(buf, "#sigma(G_{W}): %.2f", ws->var("sigmaPRResW")->getVal());
  t->DrawLatex(0.55, 0.26, buf);

  sprintf(buf, "frac(G_{W}): %.2f", ws->var("fracRes")->getVal());
  t->DrawLatex(0.55, 0.21, buf);

  tframePR->GetXaxis()->CenterTitle(true);
  tframePR->GetYaxis()->CenterTitle(true);

  // ----- Pull plot -----
  c1.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.3);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(0.25);

  RooHist* hdata = (RooHist*)tframePR->findObject("h_data", RooHist::Class());
  double xmin = 1e9;
  double xmax = -1e9;

  for (int i = 0; i < hdata->GetN(); i++)
  {
    double x, y;
    hdata->GetPoint(i, x, y);
    if (y > 0)
    { 
      if (x < xmin)
        xmin = x;
      if (x > xmax)
        xmax = x;
    }
  }

  RooHist *hpull = tframePR->pullHist("h_data", "CtPR_PEE");

  RooPlot *framePull = ws->var("ctau3D")->frame(Range(xmin, xmax));
  // RooPlot *framePull = ws->var("ctau3D")->frame(Title("Pull distribution"));
  framePull->addPlotable(hpull, "P");
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetNdivisions(505);
  framePull->GetYaxis()->SetTitleSize(0.08);
  framePull->GetYaxis()->SetLabelSize(0.07);
  framePull->GetXaxis()->SetTitleSize(0.1);
  framePull->GetXaxis()->SetLabelSize(0.08);
  framePull->Draw();

  // ----- Chi2 -----
  c1.cd(1);
  double chi2 = tframePR->chiSquare("CtPR_PEE", "h_data", fitRes->floatParsFinal().getSize());
  sprintf(buf, "#chi^{2}/ndf = %.2f", chi2);
  t->DrawLatex(0.55, 0.16, buf);

  // ----- Save -----
  std::string titlestr;
  if (fitMC)
    titlestr = Form("figs/ctau_pr_mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.pdf", ptLow, ptHigh, yLow, yHigh, isWeight);
  else
    titlestr = Form("figs/ctau_pr_data_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.pdf", ptLow, ptHigh, yLow, yHigh, isWeight);

  c1.SaveAs(titlestr.c_str());
}

void ctau_pr_mc_fit()
{
  cout << "\n=== Start ctau_pr_mc_fit() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");


  // === import inputs ===
  // PR MC, errPdfSig
  cout << "=== Import inputs ===\n";
  RooDataSet *redPRMC = nullptr;
  readInputs(redPRMC, isWeight);

  // get errPdfSig
  TFile *fInErr = nullptr;
  RooHistPdf *sigPdf = nullptr;
  readErrPdf(sigPdf, fInErr);
  

  // === make workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*redPRMC);
  ws->import(*sigPdf);

  // === build model ===
  cout << "\n=== build resolution model ===\n";
  defineCtPRRes(ws); // resolution function

  // build PEE model
  RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE", *(ws->pdf("errPdfSig")),
                      Conditional(*(ws->pdf("CtPRRes")), RooArgList(*(ws->var("ctau3D")))));
  ws->import(CtPR_PEE);

  // === fit ===
  cout << "\n=== perform fit ===\n";
  RooFitResult *fitCt_PRMC = ws->pdf("CtPR_PEE")->fitTo(*redPRMC, SumW2Error(kTRUE), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), PrintLevel(0), Save(1), NumCPU(32)); // EvalBackend("legacy")

  // ws->var("meanPRResW")->setConstant(kTRUE);
  // ws->var("sigmaPRResW")->setConstant(kTRUE); // test
  // ws->var("meanPRResN")->setConstant(kTRUE);
  // ws->var("fracRes")->setConstant(kTRUE); // test

  // === draw plot ===
  cout << "\n=== draw plot ===\n";
  // --- Check goodness of promptMCfit with per event error fit //// CtWeighted = lxy/(lxy_err) ---
  // fit은 앞에서 했고 문제 없는지 그림 그려서 보기
  // Res 변수 추가하는 과정
  // 요즘에는 그냥 ctau3DRes라는 이름으로 따로 저장해놓는다. -> 일단 오리지널 코드 방식으로
  RooRealVar *CtWeighted = new RooRealVar("CtWeighted", "#font[12]{l}_{J/#psi} / #sigma( #font[12]{l}_{J/#psi} )", -20., 20.);
  ws->import(*CtWeighted);
  const RooArgSet *thisRow = (RooArgSet *)redPRMC->get(0); // prompt MC rooData
  RooArgSet *newRow = new RooArgSet(*CtWeighted);
  RooDataSet *tempSet = new RooDataSet("tempSet", "new data set with CtWeighted", *newRow);

  // Per-event-error를 생으로 곱하고 있다.
  // Res 함수 정의를 보면 평범한 방식으로 변수 입력 방식으로 처리되어 있다.
  // 여기서는 빠르게 그림을 확인하려고 이렇게?
  for (Int_t iSamp = 0; iSamp < redPRMC->numEntries(); iSamp++)
  {
    thisRow = (RooArgSet *)redPRMC->get(iSamp);
    RooRealVar *myct = (RooRealVar *)thisRow->find("ctau3D");
    RooRealVar *mycterr = (RooRealVar *)thisRow->find("ctau3DErr");
    CtWeighted->setVal(myct->getVal() / mycterr->getVal());
    RooArgSet *tempRow = new RooArgSet(*CtWeighted);
    tempSet->add(*tempRow);
  }

  ws->factory("Gaussian::tmpGW_PRRes(CtWeighted,meanPRResW,sigmaPRResW)");
  ws->factory("Gaussian::tmpGN_PRRes(CtWeighted,meanPRResN,sigmaPRResN)");
  ws->factory("SUM::tmpCtPRRes(fracRes*tmpGW_PRRes,tmpGN_PRRes)");

  RooPlot *tempFramePR = ws->var("CtWeighted")->frame();
  tempSet->plotOn(tempFramePR, DataError(RooAbsData::SumW2), Name("h_data"));
  ws->pdf("tmpCtPRRes")->plotOn(tempFramePR, LineColor(kGreen + 1), Normalization(tempSet->sumEntries(), RooAbsReal::NumEvent), Name("CtPR_PEE"));

  drawCtauResolPlots(ws, /*isMC*/true, tempFramePR, fitCt_PRMC);

  // === save results ===
  cout << "\n=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/ctau_pr_mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  fitCt_PRMC->Write("fitCt_PRMC");
  fileOut.Close();
  
  fitCt_PRMC->Print("v");

  cout << "\n=== Finish ctau_pr_mc_fit() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
