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
  // --- data ---
  string fileNameData = Form("reduced_data/data_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  TFile fInData(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  data = (RooDataSet *)fInData.Get("redDataSB");
  if (isWeight)
  {
    cout << "[Info] data is WEIGHTED\n";
  }
  else
    cout << "[Info] data is UN-WEIGHTED" << endl;
  data->SetName("redDataSB");
}

void readErrPdf(RooHistPdf *&bkgPdf, TFile *&fInData)
{
  string fileNameData = Form("roots/ctau_err_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  fInData = TFile::Open(fileNameData.c_str());
  cout << "Load: " << fileNameData.c_str() << endl;
  bkgPdf = dynamic_cast<RooHistPdf *>(fInData->Get("errPdfBkg"));
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

void defineCtBkg(RooWorkspace *ws)
{
  // Nominal
  ws->factory("Decay::CtBkgPos(ctau3D,lambdap[0.1,0.01, 1.5],CtPRRes,RooDecay::SingleSided)");
  ws->factory("Decay::CtBkgNeg(ctau3D,lambdam[0.3,0.01, 1.5],CtPRRes,RooDecay::Flipped)");
  ws->factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.02,0.01, 5],CtPRRes,RooDecay::DoubleSided)");
  ws->factory("SUM::CtBkgSum1(fracCtBkg1[0.3,0.01,1]*CtBkgPos,CtBkgNeg)");
  ws->factory("SUM::CtBkgSum2(fracCtBkg2[0.1,0.01,1]*CtBkgSum1,CtBkgDbl)");
  ws->factory("SUM::CtBkgTot(fracCtBkg3[0.1,0.01,1]*CtPRRes,CtBkgSum2)");

  // 2 Mid
  // ws->factory("Decay::CtBkgPos(ctau3D, lambdap[0.02,0.001, 0.1], CtPRRes, RooDecay::SingleSided)");
  // ws->factory("Decay::CtBkgNeg(ctau3D, lambdam[0.02,0.001, 0.1], CtPRRes, RooDecay::Flipped)");
  // ws->factory("Decay::CtBkgDbl1(ctau3D, lambdasym1[0.2,0.001, 0.1], CtPRRes, RooDecay::DoubleSided)");
  // ws->factory("Decay::CtBkgDbl2(ctau3D, lambdasym2[0.5,0.001, 0.1], CtPRRes, RooDecay::DoubleSided)");
  // ws->factory(
  //     "RSUM::CtBkgTot("
  //     " fracCtBkgPos[0.2,0.01,1]*CtBkgPos,"
  //     " fracCtBkgNeg[0.2,0.01,1]*CtBkgNeg,"
  //     " fracCtBkgDbl1[0.2,0.01,1]*CtBkgDbl1,"
  //     " fracCtBkgDbl2[0.1,0.01,1]*CtBkgDbl2,"
  //     " CtPRRes"
  //     ")");
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

void drawCtauBkg(RooWorkspace *ws, RooDataSet *redDataSB, RooFitResult *fitRes, std::string outname="fitCtBkg")
{
  RooRealVar *ctau3D = ws->var("ctau3D");
  RooRealVar *ctau3DErr = ws->var("ctau3DErr");
  if (!ctau3D) { std::cout << "[Error] no variable ctau3D" << std::endl; return; }

  RooPlot *frame = ctau3D->frame(Title("Background fit"));
  // frame->updateNormVars(RooArgSet(*(ws->var("ctau3D"))));
  // ws->pdf("CtBkgTot_PEE")->plotOn(frame, Name("fit"), LineColor(kBlue), Normalization(redDataSB->sumEntries(), RooAbsReal::NumEvent), Name("fit"));
  redDataSB->plotOn(frame, Name("data"));
  ws->pdf("CtBkgTot")->plotOn(frame,LineStyle(kDashed), NumCPU(32), LineColor(kBlue), ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB), Name("fit"), Normalization(1, RooAbsReal::NumEvent));

  // double scalePos = (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * ws->var("fracCtBkg1")->getVal();
  // ws->pdf("CtBkgPos")->plotOn(frame, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB), LineColor(kBlue), LineStyle(kDashed), Normalization(redDataSB->sumEntries() * scalePos, RooAbsReal::NumEvent), Name("CtBkgPos"), Precision(1e-4));
  // // auto curvePos = dynamic_cast<RooCurve *>(frame->getCurve("CtBkgPos"));
  // // RooCurve *curvePos = (RooCurve *)frame->getObject(frame->numItems() - 1);

  // // // --- CtBkgNeg ---
  // double scaleNeg = (1 - ws->var("fracCtBkg3")->getVal()) * ws->var("fracCtBkg2")->getVal() * (1 - ws->var("fracCtBkg1")->getVal());

  // ws->pdf("CtBkgNeg")->plotOn(frame, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB), LineColor(kGreen), LineStyle(kDashed), Normalization(redDataSB->sumEntries() * scaleNeg, RooAbsReal::NumEvent), Name("CtBkgNeg"), Precision(1e-4));
  // // RooCurve *curveNeg = dynamic_cast<RooCurve *>(frame->getCurve("CtBkgNeg"));

  // // // --- CtBkgDbl ---
  // double scaleDbl = (1 - ws->var("fracCtBkg3")->getVal()) * (1 - ws->var("fracCtBkg2")->getVal());

  // ws->pdf("CtBkgDbl")->plotOn(frame, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB), LineColor(kOrange), LineStyle(kDashed), Normalization(redDataSB->sumEntries() * scaleDbl, RooAbsReal::NumEvent), Name("CtBkgDbl"), Precision(1e-4));
  // // RooCurve *curveDbl = dynamic_cast<RooCurve *>(frame->getCurve("CtBkgDbl"));

  // // // --- CtPRRes (Resolution) ---
  // double scaleRes = ws->var("fracCtBkg3")->getVal();
  // ws->pdf("CtPRRes")->plotOn(frame, ProjWData(RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB), LineColor(kMagenta), LineStyle(kDashed), Normalization(redDataSB->sumEntries() * scaleRes, RooAbsReal::NumEvent), Name("CtPRRes"), Precision(1e-4));
  // // RooCurve *curveRes = dynamic_cast<RooCurve *>(frame->getCurve("CtPRRes"));

  // // std::cerr << "Npos=" << (curvePos ? curvePos->GetN() : -1)
  // //           << " Nneg=" << (curveNeg ? curveNeg->GetN() : -1)
  // //           << " Ndbl=" << (curveDbl ? curveDbl->GetN() : -1)
  // //           << " Nres=" << (curveRes ? curveRes->GetN() : -1) << std::endl;
  // RooCurve* curvePos = dynamic_cast<RooCurve*>(frame->getCurve("CtBkgPos"));
  // RooCurve* curveNeg = dynamic_cast<RooCurve*>(frame->getCurve("CtBkgNeg"));
  // RooCurve* curveDbl = dynamic_cast<RooCurve*>(frame->getCurve("CtBkgDbl"));
  // RooCurve* curveRes = dynamic_cast<RooCurve*>(frame->getCurve("CtPRRes"));

  // auto sum12 = new RooCurve("sum12", "", *curvePos, *curveNeg);
  // auto sum123 = new RooCurve("sum123", "", *sum12, *curveDbl);
  // auto fit = new RooCurve("fit", "", *sum123, *curveRes);

  // fit->SetLineColor(kBlack);
  // fit->SetLineStyle(kSolid);
  // fit->SetLineWidth(3);
  // fit->SetMarkerStyle(1);
  // fit->SetDrawOption("L");

  // frame->addObject(fit);

  // === pull hist ===
  RooHist *hpull = frame->pullHist("data","fit");
  RooPlot *framePull = ctau3D->frame(Title("Pull distribution"));
  framePull->addPlotable(hpull,"P");

  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(1,2);

  c->cd(1);
  gPad->SetPad(0.0,0.3,1.0,1.0);
  gPad->SetBottomMargin(0.02);
  gPad->SetLogy();
  frame->SetMinimum(0.1);
  frame->Draw();

  double chi2 = frame->chiSquare("fit","data");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.05);
  lat.DrawLatex(0.65,0.85,Form("#chi^{2}/ndf = %.2f",chi2));

  c->cd(2);
  gPad->SetPad(0.0,0.0,1.0,0.3);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.25);
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetTitleSize(0.12);
  framePull->GetYaxis()->SetLabelSize(0.10);
  framePull->GetXaxis()->SetTitleSize(0.12);
  framePull->GetXaxis()->SetLabelSize(0.10);
  framePull->Draw();

  c->SaveAs(Form("%s.png", outname.c_str()));
  c->SaveAs(Form("%s.pdf", outname.c_str()));
}

void ctau_bkg_fit()
{
  cout << "\n=== Start ctau_bkg_fit() ===\n";
  TStopwatch t;
  t.Start();

  // silent
  setSilent();

  // rootlogon -> To global config
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // === import inputs === 
  // data-redDataSB, errPdfBkg
  cout << "=== Import inputs ===\n";
  RooDataSet *redDataSB = nullptr;
  readInputs(redDataSB, isWeight);

  // get errPdfSig
  TFile *fInErr = nullptr;
  RooHistPdf *bkgPdf = nullptr;
  readErrPdf(bkgPdf, fInErr);

  // === make workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*redDataSB);
  ws->import(*bkgPdf);

  // === define model ===
  defineCtPRRes(ws);
  defineCtBkg(ws);
  RooProdPdf CtBkgTot_PEE("CtBkgTot_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
                          Conditional(*(ws->pdf("CtBkgTot")), RooArgList(*(ws->var("ctau3D")))));
  ws->import(CtBkgTot_PEE);

  // // === fix PR MC res parameters ===
  // cout << "\n=== read MC ctau fit result ===\n";
  // string fileNamePrMc = Form("roots/ctau_pr_mc_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight);
  // TFile fInPRMC(fileNamePrMc.c_str());
  // if (fInPRMC.IsZombie())
  // {
  //   cerr << "[Error] Failed to load: " << fileNamePrMc.c_str() << "\n";
  //   return;
  // }
  // cout << "Load: " << fileNamePrMc.c_str() << "\n";
  // auto fitCt_PRMC = (RooFitResult *)fInPRMC.Get("fitCt_PRMC");

  // // fix some parameters
  // vector<string> fixList = {"sigmaPRResW", "fracRes"};
  // fixParamsMC(ws, fitCt_PRMC, fixList);


  // === perform fit ===
  cout << "\n=== perform fit ===\n";
  cout << "Maybe, it takes a few minutes\n";
  cout << " *** DATA :: N events to fit on SIDEBANDS : " << redDataSB->sumEntries() << endl;
  
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

  RooFitResult *fitCt_Bkg = ws->pdf("CtBkgTot_PEE")->fitTo(*redDataSB, SumW2Error(kTRUE), Minos(0), Save(1), ConditionalObservables(RooArgSet(*(ws->var("ctau3DErr")))), Optimize(0), NumCPU(32),PrintEvalErrors(-1), PrintLevel(-1));
  // EvalBackend("legacy"),  RecoverFromUndefinedRegions(2.5), 
  // // ws->var("fracCtBkg1")->setConstant(kTRUE);
  // // ws->var("fracCtBkg2")->setConstant(kTRUE);
  // // ws->var("fracCtBkg3")->setConstant(kTRUE);
  // // ws->var("lambdap")->setConstant(kTRUE);
  // // ws->var("lambdam")->setConstant(kTRUE);
  // // ws->var("lambdasym")->setConstant(kTRUE);

  // // === draw plot ===
  cout << "\n=== draw plot ===\n";
  drawCtauBkg(ws, redDataSB, fitCt_Bkg, "fitCtBkg_result");

  // === save result ===
  cout << "\n=== save resuts ===\n";
  gSystem->mkdir("roots", true);
  TFile fileOut(Form("roots/ctau_bkg_pT%.1f_%.1f_y%.1f_%.1f_isWeight%d.root", ptLow, ptHigh, yLow, yHigh, isWeight), "recreate");
  fitCt_Bkg->Write("fitCt_Bkg");
  fileOut.Close();

  fitCt_Bkg->Print("v");

  cout << "\n=== Finish ctau_bkg_fit() ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}