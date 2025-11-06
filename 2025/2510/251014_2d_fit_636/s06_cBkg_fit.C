#include <iostream>
#include <TStopwatch.h>
#include <string.h>
#include <RooDataSet.h>
#include <TFile.h>
#include <RooRealVar.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <RooChebychev.h>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>
// #include <>

using namespace RooFit;

using std::cout;
using std::string;

void s06_cBkg_fit()
{
  cout << "=== start s06_cBkg_fit() ===\n";
  TStopwatch time;
  time.Start();

  // set variables
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // observable ranges
  double ctMin = -0.5, ctMax = 2; // lmin, lmax: 1 for lowpT, 2 for higpT
  float errmin = 0.008, errmax = 0.3;

  // read inputs
  string fileNameCErr = "roots/ctau_err.root";
  TFile fInPRMC(fileNameCErr.c_str());
  RooDataSet *dsReduced = (RooDataSet *)fInPRMC.Get("dsReduced");
  dsReduced->SetName("dsReduced");

  // 2) dsReduced 안에서 ct 변수 '그 자체'를 꺼내옵니다.
  RooRealVar *ctau3D = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3D"));
  RooRealVar *ctau3DErr = dynamic_cast<RooRealVar *>(dsReduced->get()->find("ctau3DErr"));

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
  auto resModel = new RooAddModel("resModel", "g1+g2+g3", RooArgList(g1, g2, g3), RooArgList(f1, f2));

  auto pdt = dynamic_cast<RooHistPdf *>(fInPRMC.Get("bkgPdf")->Clone("pdt"));

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
  // RooHistFunc fdt("fdt", "f(dt) raw", RooArgList(*ctau3DErr), *binScaleBKG);
  // std::unique_ptr<RooAbsReal> I(fdt.createIntegral(RooArgSet(*ctau3DErr)));
  // RooFormulaVar fdt_unit("fdt_unit", "@0/@1", RooArgList(fdt, *I));
  // auto pdt = new RooGenericPdf("pdt", "p(dt) unit", "@0", RooArgList(fdt_unit));

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // 핵심: Prodt를 성분별로 따로 곱한다 (각 decay × pdt)
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  RooProdPdf compL("compL", "pdt*L", RooArgSet(*pdt, *decayBkgL),
                   Conditional(RooArgSet(*decayBkgL), RooArgSet(*ctau3D)));
  RooProdPdf compMid("compMid", "pdt*Mid", RooArgSet(*pdt, *decayBkgMid),
                     Conditional(RooArgSet(*decayBkgMid), RooArgSet(*ctau3D)));
  RooProdPdf compR("compR", "pdt*R", RooArgSet(*pdt, *decayBkgR),
                   Conditional(RooArgSet(*decayBkgR), RooArgSet(*ctau3D)));
  RooProdPdf compRes("compRes", "pdt*Res",
                     RooArgSet(*pdt, *resModel),
                     Conditional(RooArgSet(*resModel), RooArgSet(*ctau3D)));

  RooRealVar n_mid("n_mid", "n mid", 1000, 0, 1e12);
  RooRealVar n_left("n_left", "n left", 900, 0, 1e12);
  RooRealVar n_r("n_r", "n r", 1100, 0, 1e12);
  RooRealVar n_res("n_res", "n res", 300, 0, 1e12);
  RooAddPdf modelPEE("modelPEE", "...", RooArgList(compMid, compL, compR, compRes),
                     RooArgList(n_mid, n_left, n_r, n_res));

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

  auto fr = modelPEE.fitTo(*dsReduced, IntegrateBins(1e-6), Extended(), Save(), Offset(true), PrintLevel(0), PrintEvalErrors(-1), ConditionalObservables(RooArgSet(*ctau3DErr)), NumCPU(32), EvalBackend("legacy"));
  // 
  // 

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

  dsReduced->plotOn(frame, DataError(RooAbsData::SumW2), Name("data"));
  cout << "\npin 3\n";
  modelPEE.plotOn(frame, Name("sum_curve"), NumCPU(16));
  cout << "\npin 4\n";

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
  c1.SaveAs("figs/ctau_bkg.pdf");

  // 8) 결과 확인
  if (fr)
    fr->Print("v");

  // ----------------------------

  cout << "=== finish s06_cBkg_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}