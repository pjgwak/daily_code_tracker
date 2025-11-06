#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;

void plotVar(RooDataSet *dataset, const char *varName, const char *title, double rangeMin, double rangeMax, const char *outname);

void fitWithThreeDecay_ctau3D(RooDataSet *dataset,
                              double rangeMin, double rangeMax,
                              const char *outname)
{
  // 1) 변수 준비
  RooRealVar *ctau3D = (RooRealVar *)dataset->get()->find("ctau3D");
  if (!ctau3D)
  {
    std::cerr << "Error: variable ctau3D not found in dataset!" << std::endl;
    return;
  }
  ctau3D->setRange(rangeMin, rangeMax);

  // RooRealVar *ctau3DErr = (RooRealVar *)dataset->get()->find("ctau3DErr");
  // ctau3DErr->setRange(-0.1, 0.3);

  // Resolution model (가우시안 smearing)
  RooRealVar mean("mean", "mean", 0.0);
  // 기본 sigma1
  RooRealVar sigma1("sigma1", "resolution sigma1", 0.02, 0.01, 0.1);

  // sigma 비율 (>=1로 제한 가능)
  RooRealVar sigmaRatio("sigmaRatio", "sigma2/sigma1 ratio", 5.0, 1.0, 10.0);

  // sigma2 = sigma1 * sigmaRatio
  RooFormulaVar sigma2("sigma2", "@0 * @1", RooArgList(sigma1, sigmaRatio));

  RooGaussModel g1("g1", "Gauss Res1", *ctau3D, mean, sigma1);
  RooGaussModel g2("g2", "Gauss Res2", *ctau3D, mean, sigma2);

  RooRealVar fG1("fG1", "fraction Gauss1", 0.4, 0.01, 1.0);

  RooAddModel resModel("resModel", "Double Gaussian Resolution",
                       RooArgList(g1, g2), RooArgList(fG1));

  // 2) decay constants
  RooRealVar tauL("tauL", "lifetime left", 0.05, 0.01, 2.0);
  RooRealVar tauC("tauC", "lifetime center", 0.05, 0.01, 2.0);
  RooRealVar tauR("tauR", "lifetime right", 0.05, 0.01, 2.0);

  // 3) RooDecay PDFs
  RooDecay decayL("decayL", "Left decay", *ctau3D, tauL, resModel, RooDecay::Flipped);
  RooDecay decayC("decayC", "Central decay", *ctau3D, tauC, resModel, RooDecay::DoubleSided);
  RooDecay decayR("decayR", "Right decay", *ctau3D, tauR, resModel, RooDecay::SingleSided);

  // 4) fractions
  RooRealVar fL("fL", "fraction left", 0.3, 0.01, 1.0);
  RooRealVar fC("fC", "fraction center", 0.3, 0.01, 1.0);
  // RooFormulaVar fR("fR", "1-@0-@1", RooArgList(fL, fC)); // 나머지 → 오른쪽

  // 5) combined model
  RooAddPdf model("model", "Left+Center+Right decay model",
                  RooArgList(decayC, decayL, decayR),
                  RooArgList(fC, fL));

  // 6) 피팅
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooFitResult *fitRes = model.fitTo(*dataset,
                                     Range(rangeMin, rangeMax),
                                     Save(),
                                     Offset(),
                                     NumCPU(32),
                                     EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(1.5),
                                     RooFit::PrintEvalErrors(-1),
                                     RooFit::PrintLevel(-1));

  // --- Main frame ---
  RooPlot *frame = ctau3D->frame(Title("3 RooDecay with Double-Gauss resolution"));
  dataset->plotOn(frame, Name("data"));
  model.plotOn(frame, LineColor(kRed), LineWidth(2), Name("model"));
  model.plotOn(frame, Components("decayL"), LineColor(kBlue), LineStyle(2));
  model.plotOn(frame, Components("decayC"), LineColor(kGreen + 2), LineStyle(2));
  model.plotOn(frame, Components("decayR"), LineColor(kMagenta), LineStyle(2));

  // --- Chi2 계산 ---
  double chi2 = frame->chiSquare("model", "data");
  std::cout << "Chi2/ndf = " << chi2 << std::endl;

  // --- Pull histogram ---
  RooHist *pullHist = frame->pullHist("data", "model");

  RooPlot *framePull = ctau3D->frame(Title("Pull Distribution"));
  framePull->addPlotable(pullHist, "P");

  // --- Canvas: 상단=fit, 하단=pull ---
  TCanvas *c = new TCanvas("c_fit_ctau3D_doubleRes", "", 800, 800);
  c->Divide(1, 2);

  // pad1: 메인 플롯
  c->cd(1);
  gPad->SetPad(0.0, 0.3, 1.0, 1.0); // 위쪽 70%
  frame->Draw();

  // χ² 텍스트 추가
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  // pad2: pull 플롯
  c->cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.3); // 아래쪽 30%
  framePull->SetTitle("");
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetTitleSize(0.1);
  framePull->GetYaxis()->SetLabelSize(0.08);
  framePull->GetXaxis()->SetTitleSize(0.12);
  framePull->GetXaxis()->SetLabelSize(0.1);
  framePull->Draw();

  // 0 라인
  TLine *line = new TLine(rangeMin, 0, rangeMax, 0);
  line->SetLineStyle(2);
  line->Draw("same");

  // --- Save ---
  string outpath = "figs/fit/" + string(outname);
  c->SaveAs(outpath.c_str());

  fitRes->Print("V");

  delete c;
}

void r3_noConditional_fit()
{
  gSystem->mkdir("figs/fit", true);

  // apply weight? -> later


  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");

  // Draw plot with kinematic cut
  // ----------
  double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
  double ctMin = -1, ctMax = 3, errmin = 0.008, errmax = 0.3;

  char reduceDS[3000], reduceDS_woCtErr[3000];
  // build cuts
  sprintf(reduceDS, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"

          "&& (recoQQsign == 0)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errmin, errmax);

  RooDataSet *redData;
  redData = (RooDataSet *)data->reduce(reduceDS);

  // sideband
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");

  fitWithThreeDecay_ctau3D(redDataSBL, -1.0, 3.0, "ctau3D_fit_SBL.png");
}

void plotVar(RooDataSet *dataset,
             const char *varName,
             const char *title,
             double rangeMin, double rangeMax,
             const char *outname)
{
  // get observable
  RooRealVar *var = (RooRealVar *)dataset->get()->find(varName);
  if (!var)
  {
    std::cerr << "Error: variable " << varName << " not found in dataset!" << std::endl;
    return;
  }

  // set range
  var->setRange(rangeMin, rangeMax);

  // draw
  RooPlot *frame = var->frame(Title(title));
  dataset->plotOn(frame, Name("data"));

  TCanvas *c = new TCanvas(Form("c_%s", var->GetName()), "", 800, 600);
  frame->Draw();
  string outpath = "figs/basic/" + string(outname);
  c->SaveAs(outpath.c_str());

  delete c;
}