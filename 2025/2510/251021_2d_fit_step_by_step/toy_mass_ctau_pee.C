// file: toy_mass_ctau_pee.C
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooGaussian.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
using namespace RooFit;

void toy_mass_ctau_pee(int nEvents=100000)
{
for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  // -------- Observables --------
  RooRealVar mass("mass","mass [GeV]", 2.6, 3.5);
  RooRealVar ctau3D("ctau3D","ctau3D [mm]", -0.5, 4.0);
  RooRealVar ctau3DErr("ctau3DErr","per-event ctau error [mm]", 0.005, 0.20);

  // -------- ctau3DErr PDF as RooHistPdf --------
  // 예시 분포: 약간 오른쪽 치우친 에러 분포를 TH1으로 만든 뒤 RooHistPdf로 래핑
  TH1D hErr("hErr","ctauErr template", 50, 0.005, 0.20);
  TRandom3 rng(12345);
  for (int i=0;i<300000;i++){
    // log-normal 풍의 샘플 (평균~0.03~0.05 근방)
    double v = std::exp(rng.Gaus(std::log(0.04), 0.35));
    if (v<0.005) v = 0.005;
    if (v>0.20)  v = 0.20;
    hErr.Fill(v);
  }
  RooDataHist hErrRDH("hErrRDH","", RooArgList(ctau3DErr), &hErr);
  RooHistPdf  errPdf("errPdf","p(ctErr)", RooArgSet(ctau3DErr), hErrRDH, 1);

  // -------- Resolution model with PEE (GaussModel) --------
  RooRealVar ctMean("ctMean","resolution mean", 0.0);
  RooRealVar ctSigmaCore("ctSigmaCore","base sigma", 0.03, 0.001, 0.2); // base scale
  RooRealVar one("one","scalefactor", 1.0);
  // GaussModel(obs, mean, sigma_base, scale=1, sigmaVar=ctau3DErr) → σ_eff = sigma_base * 1 ⊕ ctErr (구현상 sigmaVar가 직접 들어갑니다)
  RooGaussModel ctRes("ctRes","ct resolution with PEE", ctau3D, ctMean, ctSigmaCore, one, ctau3DErr);

  // -------- Two-sided lifetime model (two RooDecay) --------
  // 양의 쪽 (t>0)
  RooRealVar tauP("tauP","tau(+)", 1.5, 0.1, 10.0); // 단위는 ct 범위 스케일에 맞게
  RooDecay ctPos("ctPos","positive decay", ctau3D, tauP, ctRes, RooDecay::SingleSided);

  // 음의 쪽 (t<0) - flipped (반전 꼬리)
  RooRealVar tauN("tauN","tau(-)", 0.8, 0.05, 10.0);
  RooDecay ctNeg("ctNeg","negative decay", ctau3D, tauN, ctRes, RooDecay::Flipped);

  // 합치기: ctShape = f * ctPos + (1-f) * ctNeg

  // 합치기: ctShape = f * ctPos + (1-f) * ctNeg
  RooRealVar fPos("fPos","frac(+)", 0.20, 0.0, 1.0);
  RooRealVar fNeg("fNeg","frac(+)", 0.20, 0.0, 1.0);

  RooDecay ctDbl("ctDbl","negative decay", ctau3D, tauN, ctRes, RooDecay::DoubleSided);

  RooAddPdf ctShape("ctShape","two-sided lifetime", RooArgList(ctPos, ctNeg, ctDbl), RooArgList(fPos, fNeg));

  // -------- ct model with Conditional errPdf --------
  // p(ct, ctErr) = p(ct | ctErr) * p(ctErr)
  RooProdPdf ctModel("ctModel","ct with PEE and errPdf",
                     RooArgSet(ctShape),
                     Conditional(errPdf, RooArgSet(ctau3DErr)));

  // -------- mass PDF (simple Gaussian for demo) --------
  RooRealVar m0("m0","mass mean", 3.096, 3.05, 3.14);
  RooRealVar mSig("mSig","mass sigma", 0.030, 0.005, 0.080);
  RooGaussian massPdf("massPdf","signal mass", mass, m0, mSig);

  // -------- Final model: mass * ctModel --------
  RooProdPdf model("model","mass * ct * ctErr", RooArgSet(massPdf, ctModel));

  // -------- Generate ToyMC --------
  // 세 변수(mass, ct, ctErr)를 동시에 생성 → 생성된 ctErr가 곧바로 GaussModel의 per-event error로 사용됨
  std::unique_ptr<RooDataSet> data(model.generate(RooArgSet(mass, ctau3D, ctau3DErr), nEvents));

  // -------- Fit --------
  // 조건부 관측치(ctau3DErr)를 명시 (모델에 이미 Conditional이 있지만 fitTo에서도 명확화)
  std::unique_ptr<RooFitResult> fitRes(
    model.fitTo(*data, Save(), ConditionalObservables(RooArgSet(ctau3DErr)),
                PrintEvalErrors(-1), PrintLevel(-1))
  );
  fitRes->Print("v");

  // -------- Quick plots --------
  // 1) mass
  TCanvas c1("c1","mass",800,600);
  RooPlot* fMass = mass.frame(Title("mass"));
  data->plotOn(fMass, Binning(80));
  // mass는 ct, ctErr에 대해 적분 필요 → ProjWData로 데이터 전달
  model.plotOn(fMass, ProjWData(RooArgSet(ctau3D, ctau3DErr), *data), LineColor(kBlue));
  fMass->Draw();
  c1.SaveAs("toy_mass.png");

  // 2) ctau3D
  TCanvas c2("c2","ctau3D",800,600);
  RooPlot* fCt = ctau3D.frame(Title("ctau3D"));
  data->plotOn(fCt, Binning(80));
  // ct는 mass, ctErr에 대해 적분 → ProjWData에 둘 다 명시
  model.plotOn(fCt, ProjWData(RooArgSet(mass, ctau3DErr), *data), LineColor(kRed));
  fCt->Draw();
  c2.SaveAs("toy_ctau3D.png");

  // 3) ctau3DErr
  TCanvas c3("c3","ctau3DErr",800,600);
  RooPlot* fErr = ctau3DErr.frame(Title("ctau3DErr"));
  data->plotOn(fErr, Binning(50));
  // errPdf만 투영해서 비교
  errPdf.plotOn(fErr, LineColor(kGreen+2));
  fErr->Draw();
  c3.SaveAs("toy_ctErr.png");
}
