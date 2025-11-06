#include <TStopwatch.h>
#include <RooArgList.h>

using namespace RooFit;

// === kinematics === -> Read from config!
  float ptLow = 13, ptHigh = 13.5;
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
  // --- NPMC ---
  string fileNameData = Form("reduced_data/mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("redNPMC");
  if (isWeight)
  {
    cout << "[Info] data is WEIGHTED\n";
  }
  else
    cout << "[Info] data is UN-WEIGHTED" << endl;
  data->SetName("redNPMC");
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

void drawCtNPTrue(RooWorkspace *ws, RooDataSet *redNPMC, string titlestr, RooFitResult *fitCt_NP)
{
  RooPlot *tframetrue = ws->var("ctau3Dtrue")->frame(); // Range("fitRange")
  redNPMC->plotOn(tframetrue, Name("h_data"), DataError(RooAbsData::SumW2));
  ws->pdf("CtNPTrue")->plotOn(tframetrue, LineColor(kOrange), Name("CtNPTrue"));
  // Range("fitRange"), NormRange("fitRange"),

  tframetrue->GetXaxis()->CenterTitle(1);
  tframetrue->GetYaxis()->CenterTitle(1);

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(12);
  t.SetTextSize(0.035);

  TCanvas c1("c1", "c1", 650, 750);
  c1.Divide(1, 2);

  c1.cd(1);
  gPad->SetPad(0.0, 0.33, 1.0, 1.0);
  gPad->SetBottomMargin(0.02);
  gPad->SetLogy();

  tframetrue->Draw();

  double truemax = tframetrue->GetMaximum();
  tframetrue->SetMinimum(0.1);
  tframetrue->SetMaximum(truemax * 5);

  char buf[512];
  // sprintf(buf, "coefExpNPTrue: %.2f #pm %.1e", ws->var("coefExpNPTrue")->getVal(), ws->var("coefExpNPTrue")->getError());
  // t.DrawLatex(0.43, 0.84, buf);

  // sprintf(buf, "sigmaNPTrue: %.1e #pm %.1e", ws->var("sigmaNPTrue")->getVal(), ws->var("sigmaNPTrue")->getError());
  // t.DrawLatex(0.43, 0.76, buf);

  double chi2 = tframetrue->chiSquare("CtNPTrue", "h_data", fitCt_NP->floatParsFinal().getSize());
  cout << "Chi2/ndf: " << chi2 << "\n";

  sprintf(buf, "#chi^{2}/ndf = %.2f", chi2);
  t.DrawLatex(0.43, 0.68, buf);

  c1.cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.33);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.32);
  gPad->SetLogy(0);

  RooHist *hdata = (RooHist *)tframetrue->findObject("h_data", RooHist::Class());
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

  // RooHist *hpull = tframetrue->pullHist("h_data", "CtNPTrue");
  RooHist *hpull = tframetrue->residHist("h_data", "CtNPTrue", kTRUE);
  RooPlot *framePull = ws->var("ctau3Dtrue")->frame(Range(-1, 9));

  double *ypulls = hpull->GetY();
  double chi22 = 0;
  int ndof = 0;
  unsigned int nFullBins = 0;
  for (int i = 0; i < hdata->GetN(); i++)
  {
    // for (int i = 31; i < 56; i++) {
    if (ypulls[i] == 0)
      continue;
    cout << "ypulls[" << i << "] = " << ypulls[i] << "\n";
    chi22 += ypulls[i] * ypulls[i];
    nFullBins++;
  }
  unsigned int nFitPar = fitCt_NP->floatParsFinal().getSize();
  ndof = nFullBins - nFitPar;
  cout << "#chi^{2}/ndof = " << chi22 << "/" <<ndof << "\n";

  // RooPlot *framePull = ws->var("ctau3Dtrue")->frame(Title("Pull distribution"));
  framePull->addPlotable(hpull, "P");
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetNdivisions(505);
  framePull->GetYaxis()->SetTitleSize(0.11);
  framePull->GetYaxis()->SetLabelSize(0.10);
  framePull->GetXaxis()->SetTitleSize(0.13);
  framePull->GetXaxis()->SetLabelSize(0.12);
  // framePull->GetYaxis()->SetRangeUser(-5, 5);
  framePull->Draw();

  c1.SaveAs(Form("figs/ctau_true_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.pdf", ptLow, ptHigh, yLow, yHigh, isWeight));

  delete tframetrue;
}

void defineCtNP(RooWorkspace *ws)
{
  // --- simple ctau NP pdf ---
  // ws->factory("GaussModel::BResTrue(ctau3Dtrue,mean[0.0],sigmaNPTrue[0.00002,0.0000001,0.1])"); // resolution from B meson (before CtPRRes)
  // ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // ws->factory("Decay::CtNPTrue(ctau3Dtrue,coefExpNPTrue[0.1,0.00001,10],BResTrue,RooDecay::SingleSided)");

  // --- two decay model ---
  ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // ws->factory("GaussModel::BResTrue(ctau3Dtrue,mean[0.0],sigmaNPTrue[0.00002,0.0000001,0.1])"); // resolution from B meson (before CtPRRes)

  ws->factory("coefExpNPTrue1[0.10, 0.001, 1.0]");
  ws->factory("rCoefExpNPTrue[2.0, 1.01, 50.0]");
  ws->factory("expr::coefExpNPTrue2('coefExpNPTrue1/rCoefExpNPTrue',{coefExpNPTrue1,rCoefExpNPTrue})");

  ws->factory("Decay::CtNPTrue1(ctau3Dtrue, coefExpNPTrue1, BResTrue, RooDecay::SingleSided)");
  ws->factory("Decay::CtNPTrue2(ctau3Dtrue, coefExpNPTrue2, BResTrue, RooDecay::SingleSided)");

  ws->factory("SUM::CtNPTrue(fFracCtNPTot[0.5,0.0,1.0]*CtNPTrue1, CtNPTrue2)");

  // --- three decay model ---
  // ws->factory("TruthModel::BResTrue(ctau3Dtrue)");
  // ws->factory("coefExpNPTrue1[0.1,0.001,1]");

  // ws->factory("r21[3.0,1,10.0]");   // lambda2 / lambda1
  // ws->factory("r31[6.0,1,20.0]");   // lambda3 / lambda1

  // ws->factory("expr::coefExpNPTrue2('coefExpNPTrue1*r21', coefExpNPTrue1,r21)");
  // ws->factory("expr::coefExpNPTrue3('coefExpNPTrue2*r31', coefExpNPTrue2,r31)");

  // // Decay PDFs
  // ws->factory("Decay::CtNPTrue1(ctau3Dtrue, coefExpNPTrue1, BResTrue, RooDecay::SingleSided)");
  // ws->factory("Decay::CtNPTrue2(ctau3Dtrue, coefExpNPTrue2, BResTrue, RooDecay::SingleSided)");
  // ws->factory("Decay::CtNPTrue3(ctau3Dtrue, coefExpNPTrue3, BResTrue, RooDecay::SingleSided)");

  // ws->factory("SUM::CtNPTrue(fFrac1[0.3,0.0,1.0]*CtNPTrue1, fFrac2[0.3,0.0,1.0]*CtNPTrue2, CtNPTrue3)");

  // // -> 이 부분은 final fit으로 보낼 수도 있다.
  // // --- build NP Pdf convoluted PRRes ---
  // RooFormulaVar sigmaNPResW("sigmaNPResW", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResW")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  // ws->import(sigmaNPResW);
  // RooFormulaVar sigmaNPResN("sigmaNPResN", "sqrt((@0*@1)**2+(@2)**2)", RooArgList(*(ws->var("sigmaPRResN")), *(ws->var("ctau3DErr")), *(ws->var("sigmaNPTrue"))));
  // ws->import(sigmaNPResN);
  // ws->factory("GaussModel::GW_NPRes(ctau3D,meanPRResW,sigmaNPResW)");
  // ws->factory("GaussModel::GN_NPRes(ctau3D,meanPRResN,sigmaNPResN)");
  // ws->factory("AddModel::CtNPRes({GW_NPRes,GN_NPRes},{fracRes})");

  // // final model
  // float coefExpNPTrueVal = ws->var("coefExpNPTrue")->getVal();
  // ws->factory("Decay::CtNPTot(ctau3D,coefExpNPTrue,CtNPRes,RooDecay::SingleSided)");

  //// **** check NPMC Reco
  // 이거 왜 있지????
  // drawCtNPReco(ws, redNPMC, titlestr, opt);
  return;
}

void ctau_np_mc_fit()
{
  cout << "\n=== Start ctau_np_mc_fit() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");


  // // === import inputs ===
  cout << "=== Import inputs ===\n";
  RooDataSet *redNPMC = nullptr;
  readInputs(redNPMC, isWeight);

  // === make workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*redNPMC);

  // ws->var("ctau3Dtrue")->setRange("fitRange", 0.001, 6);
  // ws->var("ctau3Dtrue")->setBins(10);

  // === build model ===
  cout << "\n=== build ctau True model ===\n";
  defineCtNP(ws);

  // === fit ===
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  cout << "\n=== perform fit ===\n";
  RooFitResult *fitCt_NPMC = ws->pdf("CtNPTrue")->fitTo(*redNPMC, Save(), NumCPU(8), Strategy(2), SumW2Error(false), PrintLevel(-1), RooFit::PrintEvalErrors(-1));

  // === draw plot ===
  cout << "\n=== draw plot ===\n";
  drawCtNPTrue(ws, redNPMC, "", fitCt_NPMC);

  // ws->var("sigmaNPTrue")->setConstant(kTRUE);

  // === save results ===
  cout << "\n=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/ctau_true_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  fitCt_NPMC->Write("fitCt_NPMC");
  fileOut.Close();

  fitCt_NPMC->Print("V");

  cout << "\n=== Finish ctau_np_mc_fit() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}
