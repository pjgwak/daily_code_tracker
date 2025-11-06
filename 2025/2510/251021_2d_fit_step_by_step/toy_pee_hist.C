// file: toy_pee_hist.C
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

void toy_pee_hist(int nEvents=100000)
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
  RooRealVar fPos("fPos","frac(+)", 0.70, 0.0, 1.0);
  RooAddPdf ctShape("ctShape","two-sided lifetime", RooArgList(ctPos, ctNeg), RooArgList(fPos));

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
  mass.setBins(80, "plot");        // 플롯 x-축 bin
  ctau3D.setBins(80, "plot");
  ctau3DErr.setBins(50, "plot");

  // --- 0) 투영용 RooDataHist들 준비 ---
  // mass 플롯: ctau3D, ctau3DErr로 투영
  RooDataHist proj_ct_cterr("proj_ct_cterr","proj for mass",
                            RooArgSet(ctau3D, ctau3DErr), *data);

  // ct 플롯: mass, ctau3DErr로 투영
  RooDataHist proj_m_err("proj_m_err","proj for ct",
                        RooArgSet(mass, ctau3DErr), *data);

  // ctErr 플롯: 이미 errPdf의 소스 히스토그램 hErrRDH가 있음(더 빠름)
  // RooDataHist hErrRDH("hErrRDH","", RooArgList(ctau3DErr), &hErr); // 앞서 생성했다고 가정

  // --- 1) mass ---
  TCanvas c1("c1","mass",800,600);
  RooPlot* fMass = mass.frame(Title("mass"));
  data->plotOn(fMass, Binning(80), DataError(RooAbsData::SumW2));

  model.plotOn(fMass,
    ProjWData(RooArgSet(ctau3D, ctau3DErr), proj_ct_cterr), // <-- dataset → datahist 로 교체
    Range(mass.getMin(), mass.getMax()),
    Normalization(data->sumEntries(), RooAbsReal::NumEvent),
    LineColor(kBlue));

  fMass->Draw();
  c1.SaveAs("toy_mass.png");

  // --- 2) ctau3D ---
  TCanvas c2("c2","ctau3D",800,600);
  RooPlot* fCt = ctau3D.frame(Title("ctau3D"));
  data->plotOn(fCt, Binning(80), DataError(RooAbsData::SumW2));

  model.plotOn(fCt,
    ProjWData(RooArgSet(mass, ctau3DErr), proj_m_err),       // <-- dataset → datahist 로 교체
    Range(ctau3D.getMin(), ctau3D.getMax()),
    Normalization(data->sumEntries(), RooAbsReal::NumEvent),
    LineColor(kRed));

  fCt->Draw();
  c2.SaveAs("toy_ctau3D.png");

  // --- 3) ctau3DErr ---
  // 데이터는 그대로 찍고, 모델 비교는 errPdf 소스 히스토그램을 사용해 즉시
  TCanvas c3("c3","ctau3DErr",800,600);
  RooPlot* fErr = ctau3DErr.frame(Title("ctau3DErr"));
  data->plotOn(fErr, Binning(50), DataError(RooAbsData::SumW2));

  // errPdf는 hErrRDH를 그대로 그리는게 가장 빠름
  errPdf.plotOn(fErr, LineColor(kGreen+2));
  fErr->Draw();
  c3.SaveAs("toy_ctErr.png");
}
