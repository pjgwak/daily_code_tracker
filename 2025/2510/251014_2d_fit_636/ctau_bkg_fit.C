// #include "utils/ctau_fit_utils.h"
// #include "cfgs/pp_run2_pt40p50_y0_1p6.C"
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFormulaVar.h>
#include <RooFitResult.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>

using namespace RooFit;
using std::cerr;
using std::cout;
using std::string;

// void importErrPdf(RooWorkspace *ws, FitConfigCtauBkg &cfg);
// void defineCtauBkg(RooWorkspace *ws, FitConfigCtauBkg &cfg);
// void drawCtauFit(RooWorkspace *ws, RooRealVar &ctau3D, RooAbsPdf &model, RooAbsData &data, const FitConfigCtauBkg &cfg, const RooFitResult *fitResult, int nBins = 80, const char *outPrefix = "ctau_pr");

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
  float ptLow = 40, ptHigh = 50;
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
      "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f)",
      ptLow, ptHigh, yLow, yHigh, massLow, massHigh);
  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim
  TString osCut = "(recoQQsign == 0) && (ctau3D > -0.2 && ctau3D < 2 &&ctau3DErr > 0.0001 && ctau3DErr < 0.04)";

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

  // 2) dsReduced 안에서 ct 변수 '그 자체'를 꺼내옵니다.
  RooRealVar *ctau3D = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3D"));
  RooRealVar *ctau3DErr = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3DErr"));

  // 4) 범위는 '데이터가 실제로 가진' 값으로 설정
  double xmin, xmax;
  dsReduced->getRange(*ctau3D, xmin, xmax); // 데이터 기반
  ctau3D->setRange(-0.2, 2);
  ctau3D->setRange("fitRange", -0.2, 2);
  ctau3DErr->setRange(0.0001, 0.04);
  ctau3DErr->setRange("errRange", 0.0001, 0.04);

  // 5) 2-가우스 모델 (변수는 '새로 만들지 않고' ctau3D을 그대로 사용)
  RooRealVar mu("mu", "mean", 0.0);
  RooRealVar s1("sigma1", "sigma1", 0.08, 0.002, 0.1);
  RooRealVar s2("sigma2", "sigma2", 0.25, 0.01, 5);
  RooRealVar s3("sigma3", "sigma3", 0.5, 0.01, 10);

  RooRealVar f1("f1", "frac g1", 0.5, 0.0, 1.0); // g1 fraction
  RooRealVar f2("f2", "frac g2", 0.3, 0.0, 1.0); // g2 fraction
  RooRealVar one("one", "scale", 1.0);           // resolution scale (고정)

  // --- 개별 GaussModel ---
  RooGaussModel g1("g1", "Gauss sigma1", *ctau3D, mu, s1, one, *ctau3DErr);
  RooGaussModel g2("g2", "Gauss sigma2", *ctau3D, mu, s2, one, *ctau3DErr);
  RooGaussModel g3("g3", "Gauss sigma3", *ctau3D, mu, s3, one, *ctau3DErr);

  // --- 3-Gauss mixture ---
  // f1 = g1 비율, f2 = g2 비율, g3 비율 = 1 - f1 - f2
  auto resModel = new RooAddModel("resModel", "g1+g2+g3", RooArgList(g1, g2, g3), RooArgList(f1, f2));

  cout << "\npin 1\n";

  TFile fin(Form("roots/pp_run2_pT12_15_y0_1p6/ctau_err_pT%.1f_%.1f_y%.1f_%.1f.root", ptLow, ptHigh, yLow, yHigh), "read");
  if (fin.IsZombie())
  {
    cerr << "[Error] Can't open MC fit file\n";
    return;
  }

  auto binScaleBKG = dynamic_cast<RooDataHist *>(fin.Get("binScaleBKG")->Clone("binScaleBKG"));

  RooRealVar tau_mid("tau_mid", "scale factor middle", 0.10, 0.01, 1.0);
  RooRealVar tau_mid2("tau_mid2", "scale factor middle", 0.20, 0.01, 2.0);

  RooRealVar r_left("r_left", "tau_left / tau_mid", 2.0, 0.20, 20.0);
  RooRealVar r_right("r_right", "tau_right / tau_mid", 4.0, 0.50, 50.0);

  RooFormulaVar tau_left("tau_left", "@0*@1", RooArgList(tau_mid, r_left));
  RooFormulaVar tau_right("tau_right", "@0*@1", RooArgList(tau_mid, r_right));

  // --- decay models (그대로) ---
  auto decayBkgL = new RooDecay("decayBkgL", "Left decay", *ctau3D, tau_left, *resModel, RooDecay::Flipped);
  auto decayBkgMid1 = new RooDecay("decayBkgMid1", "Center decay1", *ctau3D, tau_mid, *resModel, RooDecay::DoubleSided);
  auto decayBkgMid2 = new RooDecay("decayBkgMid2", "Center decay2", *ctau3D, tau_mid2, *resModel, RooDecay::DoubleSided);
  RooRealVar fmid12("fmid12", "frac of Mid1", 0.3, 0.0, 1.0);
  auto decayBkgMid = new RooAddPdf("decayBkgMid", "sum of 2 decays",
                                   RooArgList(*decayBkgMid1, *decayBkgMid2),
                                   RooArgList(fmid12));
  auto decayBkgR = new RooDecay("decayBkgR", "Right decay", *ctau3D, tau_right, *resModel, RooDecay::SingleSided);

  // --- p(dt) (그대로) ---
  RooHistFunc fdt("fdt", "f(dt) raw", RooArgList(*ctau3DErr), *binScaleBKG);
  std::unique_ptr<RooAbsReal> I(fdt.createIntegral(RooArgSet(*ctau3DErr)));
  RooFormulaVar fdt_unit("fdt_unit", "@0/@1", RooArgList(fdt, *I));
  auto pdt = new RooGenericPdf("pdt", "p(dt) unit", "@0", RooArgList(fdt_unit));

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // 핵심: Prodt를 성분별로 따로 곱한다 (각 decay × pdt)
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  RooProdPdf compL("compL", "pdt*L", RooArgSet(*pdt, *decayBkgL),
                   Conditional(RooArgSet(*decayBkgL), RooArgSet(*ctau3D)));
  RooProdPdf compMid("compMid", "pdt*Mid", RooArgSet(*pdt, *decayBkgMid),
                     Conditional(RooArgSet(*decayBkgMid), RooArgSet(*ctau3D)));
  RooProdPdf compR("compR", "pdt*R", RooArgSet(*pdt, *decayBkgR),
                   Conditional(RooArgSet(*decayBkgR), RooArgSet(*ctau3D)));

  // --- 최종 합성 (이전: pdt * model  →  변경: compMid + compL + compR)
  RooRealVar f_bkg_mid("f_bkg_mid", "fraction middle", 0.3, 0.0, 1.0);
  RooRealVar f_bkg_left("f_bkg_left", "fraction left", 0.3, 0.0, 1.0);
  // 3성분 → 2개 frac (마지막은 1 - f_mid - f_left)
  RooAddPdf modelPEE("modelPEE", "sum of pdt-weighted decays",
                     RooArgList(compMid, compL, compR),
                     RooArgList(f_bkg_mid, f_bkg_left));
  // RooProdPdf decayLPlot("decayLPlot", "",
  //                      RooArgSet(*h2, *decayBkgL),
  //                      Conditional(RooArgSet(*decayBkgL), RooArgSet(*ctau3D)));
  // RooProdPdf decayRPlot("decayRPlot", "",
  //                       RooArgSet(*h2, *decayBkgR),
  //                       Conditional(RooArgSet(*decayBkgR), RooArgSet(*ctau3D)));
  // RooProdPdf decayM1Plot("decayM1Plot", "",
  //                       RooArgSet(*h2, *decayBkgMid1),
  //                       Conditional(RooArgSet(*decayBkgMid1), RooArgSet(*ctau3D)));
  // RooProdPdf decayM2Plot("decayM2Plot", "",
  //                       RooArgSet(*h2, *decayBkgMid2),
  //                       Conditional(RooArgSet(*decayBkgMid2), RooArgSet(*ctau3D)));
  // RooProdPdf modelPEE("modelPEE", "CtPDF with PEE",
  //                         RooArgSet(*h2, *model),
  //                         Conditional(RooArgSet(*model), RooArgSet(*ctau3D)));

  // 6) 피팅 (옵션 최소화 → 먼저 움직이는지 확인)
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

  // auto fr = modelPEE.fitTo(*dsReduced, Save(), Range("fitRange"), Offset(true), PrintLevel(0), PrintEvalErrors(-1), NumCPU(32), EvalBackend("legacy"));

  auto fr = modelPEE.fitTo(*dsReduced, Save(), Offset(true), PrintLevel(-1), PrintEvalErrors(-1), ConditionalObservables(RooArgSet(*ctau3DErr)), NumCPU(32));
  // EvalBackend("legacy")

  cout << "\npin 2\n";

  // 7) 플롯 (불필요한 ProjWData 절대 X)
  TCanvas c1("c1", "", 800, 800);

  // pad 분할
  c1.Divide(1, 2);
  TPad *pad1 = (TPad *)c1.cd(1);
  TPad *pad2 = (TPad *)c1.cd(2);

  pad1->SetPad(0.0, 0.3, 1.0, 1.0); // 위쪽 큰 영역
  pad2->SetPad(0.0, 0.0, 1.0, 0.3); // 아래 pull 영역
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.25);

  // --- main frame ---
  RooPlot *frame = ctau3D->frame(Title(""));

  // 데이터 (이름 지정 중요)
  

  // 모델 (이름 지정 중요)
  // modelPEE.plotOn(frame,
  //                 // Range("fitRange"),
  //                 // NormRange("fitRange"),
  //                 ProjWData(RooArgSet(*ctau3DErr), *dsReduced, kTRUE), NumCPU(8),
  //                 Normalization(dsReduced->sumEntries(), RooAbsReal::NumEvent), Name("model"),
  //                 Precision(-1));
  // modelPEE.fixAddCoefRange("fitRange");

  // Use exactly same fitRanges

  dsReduced->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));
  modelPEE.plotOn(frame, Name("sum_curve"));

  // decayLPlot.plotOn(frame, ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kBlue), LineStyle(kDashed), Normalization(dsReduced->sumEntries() * f_bkg_left.getVal(), RooAbsReal::NumEvent));
  // RooCurve *curve1 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  // decayRPlot.plotOn(frame, ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kGreen), LineStyle(kDashed), Normalization(dsReduced->sumEntries() * (1 - f_bkg_left.getVal() - f_bkg_mid.getVal()), RooAbsReal::NumEvent));
  // RooCurve *curve2 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  // RooCurve *sum_curve = new RooCurve("sum_curve", "Sum of curves", *curve1, *curve2);

  // sum_curve->SetLineColor(kRed);
  // sum_curve->SetLineStyle(kSolid);
  // frame->addObject(sum_curve);

  // frame->Draw();

  // decayLPlot.plotOn(frame, Range("fitRange", "errRange"), ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kMagenta), LineStyle(kDashed), Normalization(dsReduced->sumEntries() * f_bkg_left.getVal(), RooAbsReal::NumEvent), Name("decayLPlot"));
  // RooCurve *curve1 = (RooCurve *)frame->findObject("decayLPlot");

  // decayRPlot.plotOn(frame, ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kGreen), LineStyle(kDashed), Normalization(dsReduced->sumEntries() * (1 - f_bkg_left.getVal() - f_bkg_mid.getVal()), RooAbsReal::NumEvent), Name("decayRPlot"));
  // RooCurve *curve2 = (RooCurve *)frame->findObject("decayRPlot");

  // decayM1Plot.plotOn(frame, ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kOrange), LineStyle(kDashed), Normalization((dsReduced->sumEntries() * f_bkg_mid.getVal()*fmid12.getVal()), RooAbsReal::NumEvent));
  // RooCurve *curve3 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  // decayM2Plot.plotOn(frame, ProjWData(RooArgSet(*ctau3DErr), *dsReduced), LineColor(kOrange), LineStyle(kDashed), Normalization((dsReduced->sumEntries() * f_bkg_mid.getVal() * (1-fmid12.getVal())), RooAbsReal::NumEvent));
  // RooCurve *curve4 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  // RooCurve *sum1 = new RooCurve("sum1", "", *curve1, *curve2);
  // RooCurve *sum2 = new RooCurve("sum1", "", *curve3, *curve4);
  // RooCurve *sum_curve = new RooCurve("sum_curve", "", *sum1, *sum2);

  // sum_curve->SetLineColor(kBlue);
  // sum_curve->SetLineStyle(kSolid);
  // sum_curve->SetLineWidth(3);
  // sum_curve->SetMarkerStyle(1);
  // sum_curve->SetDrawOption("L"); // draw only line

  // // 프레임에 추가
  // frame->addObject(sum_curve);


  // --- chi2 계산 ---
  double chi2ndf = frame->chiSquare("sum_curve", "data",
                                    /*nFitParam=*/4);
  // 마지막 인자: 자유 파라미터 개수 (mu,s1,s2,f1 → 4개로 가정)

  // 위쪽 pad: 데이터 + 모델
  pad1->cd();
  frame->Draw();

  // // χ²/ndf 텍스트 출력
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.6, 0.85, Form("#chi^{2}/ndf = %.2f", chi2ndf));

  // --- pull histogram ---
  RooHist *hpull = frame->pullHist("data", "sum_curve");
  RooPlot *frame_pull = ctau3D->frame(Title("Pull"));
  frame_pull->addPlotable(hpull, "P");

  // y=0 기준선
  TLine *line = new TLine(ctau3D->getMin(), 0, ctau3D->getMax(), 0);
  line->SetLineStyle(2);
  line->SetLineColor(kRed);

  // 아래 pad: pull plot
  pad2->cd();
  frame_pull->Draw();
  line->Draw("same");

  // 저장
  c1.SaveAs("test_with_pull_chi2.png");

  // 8) 결과 확인
  if (fr)
    fr->Print("v");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}