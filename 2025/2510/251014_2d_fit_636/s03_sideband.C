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

// === define helper functions ===
// --- calculate mass bkg scale value ---
inline double computeScaleF_cheby(const RooArgList &allPars,
                                  int order, // 1~6
                                  double massMin, double massMax,
                                  double sbL_lo, double sbL_hi,
                                  double sig_lo, double sig_hi,
                                  double sbR_lo, double sbR_hi,
                                  const char *coeffPrefix = "sl") // 계수 접두: "sl"
{
  // --- helpers ---
  auto getVal = [&](const char *name) -> double
  {
    if (auto *obj = allPars.find(name))
    {
      if (auto *rr = dynamic_cast<RooAbsReal *>(obj))
        return rr->getVal();
    }
    // 찾지 못하면 NaN으로 신속 실패
    return std::numeric_limits<double>::quiet_NaN();
  };

  auto toUnitX = [&](double m) -> double
  {
    return 2.0 * (m - massMin) / (massMax - massMin) - 1.0; // [massMin,massMax] → [-1,1]
  };

  auto chebT = [&](int n, double x) -> double
  {
    if (n == 0)
      return 1.0;
    if (n == 1)
      return x;
    double Tn_2 = 1.0, Tn_1 = x, Tn = 0.0;
    for (int k = 2; k <= n; ++k)
    {
      Tn = 2.0 * x * Tn_1 - Tn_2;
      Tn_2 = Tn_1;
      Tn_1 = Tn;
    }
    return Tn;
  };

  auto chebVal_roofit = [&](double m, const std::vector<double> &c) -> double
  {
    // RooChebychev 정의: f = 1*T0 + c1*T1 + ... + cN*TN  (상수항 1 내장)
    const double x = toUnitX(m);
    double v = 1.0;
    for (int i = 1; i <= int(c.size()); ++i)
      v += c[i - 1] * chebT(i, x);
    return v;
  };

  auto integrateSimpson = [&](auto &&f, double a, double b, int nSeg = 2048) -> double
  {
    if (b <= a)
      return 0.0;
    if (nSeg % 2)
      ++nSeg;
    const double h = (b - a) / nSeg;
    double s = f(a) + f(b);
    for (int i = 1; i < nSeg; ++i)
    {
      const double x = a + i * h;
      s += f(x) * ((i % 2) ? 4.0 : 2.0);
    }
    return s * (h / 3.0);
  };

  // --- 입력 검사 ---
  if (order < 0 || order > 6)
  {
    ::Error("computeScaleF_cheby", "Unsupported Cheby order %d (allowed 0..6).", order);
    return std::numeric_limits<double>::quiet_NaN();
  }

  // --- 계수 읽기: sl1..slN만 사용 (검색/대체 없음) ---
  std::vector<double> c;
  c.reserve(order);
  for (int i = 1; i <= order; ++i)
  {
    char nm[32];
    std::snprintf(nm, sizeof(nm), "%s%d", coeffPrefix, i); // 예: "sl1"
    double v = getVal(nm);
    if (!std::isfinite(v))
    {
      ::Error("computeScaleF_cheby", "Coefficient '%s' not found or NaN.", nm);
      return std::numeric_limits<double>::quiet_NaN();
    }
    c.push_back(v);
  }

  // --- 적분 및 스케일 계산 ---
  auto f = [&](double m)
  { return chebVal_roofit(m, c); };
  const double I_sig = integrateSimpson(f, sig_lo, sig_hi);
  const double I_sb = integrateSimpson(f, sbL_lo, sbL_hi) + integrateSimpson(f, sbR_lo, sbR_hi);

  if (std::abs(I_sb) < 1e-300)
  {
    ::Error("computeScaleF_cheby", "Sideband integral is zero. Check ranges.");
    return std::numeric_limits<double>::quiet_NaN();
  }
  return I_sig / I_sb;
}

// --- subtract histograms ---
// - binScaleBKG : k*SB (가중치+오차)
// - binSubtrSIG : SIG - k*SB (가중치+오차)
// - dk          : scale factor의 표준오차(없으면 0)
// - floorW      : 음수/제로 방지용 최소 가중치(오차에는 적용 안 함)
inline void subtractSidebands_inplace_propagate_TH1(RooDataHist *binSubtrSIG,   // out
                                                    RooDataHist *binScaleBKG,   // out
                                                    const RooDataHist *binSIG,  // in
                                                    const RooDataHist *binSB,   // in
                                                    RooRealVar *obsVar,         // e.g. ctau3DErr
                                                    double k,                   // scaleF
                                                    const std::string &varName, // "ctau3DErr"
                                                    double dk = 0.0,
                                                    double floorW = 1e-6)
{
  if (!binSubtrSIG || !binScaleBKG || !binSIG || !binSB || !obsVar)
  {
    ::Error("subtractSidebands_inplace_propagate_TH1", "Null input.");
    return;
  }
  const int nSIG = binSIG->numEntries();
  const int nSB = binSB->numEntries();
  if (nSIG != nSB)
  {
    ::Error("subtractSidebands_inplace_propagate_TH1", "Different binning: %d vs %d", nSIG, nSB);
    return;
  }

  // === 입력 RooDataHist → TH1로 변환해서 bin 오차를 읽는다 (버전 독립)
  std::unique_ptr<TH1> hSIG(binSIG->createHistogram("hSIG_tmp", *obsVar));
  std::unique_ptr<TH1> hSB(binSB->createHistogram("hSB_tmp", *obsVar));
  if (!hSIG || !hSB)
  {
    ::Error("subtractSidebands_inplace_propagate_TH1", "createHistogram failed.");
    return;
  }
  const int nb = hSIG->GetNbinsX();
  if (nb != hSB->GetNbinsX() || nb != nSIG)
  {
    ::Error("subtractSidebands_inplace_propagate_TH1", "Bin mismatch: TH1(%d) vs RooDataHist(%d).", nb, nSIG);
    return;
  }

  // === 각 빈 좌표/가중치/오차 읽어서 출력에 set()
  for (int i = 0; i < nSIG; ++i)
  {
    const RooArgSet *argSIG = binSIG->get(i);
    const RooArgSet *argSB = binSB->get(i);
    if (!argSIG || !argSB)
      continue;
    if (!argSIG->find(varName.c_str()) || !argSB->find(varName.c_str()))
    {
      ::Error("subtractSidebands_inplace_propagate_TH1", "Var '%s' not found in argset.", varName.c_str());
      return;
    }

    const int j = i + 1;                        // TH1 bin index (1..nbins)
    const double wSIG = hSIG->GetBinContent(j); // Σw
    const double eSIG = hSIG->GetBinError(j);   // sqrt(Σw²) or Poisson

    const double wSB = hSB->GetBinContent(j);
    const double eSB = hSB->GetBinError(j);

    // k*SB
    const double wSBs = k * wSB;
    const double varSBs = (k * k) * (eSB * eSB) + (dk > 0.0 ? (wSB * wSB) * (dk * dk) : 0.0);
    const double eSBs = std::sqrt(std::max(0.0, varSBs));

    // SIG - k*SB
    const double wSUB = wSIG - wSBs;
    const double varSUB = (eSIG * eSIG) + varSBs;
    const double eSUB = std::sqrt(std::max(0.0, varSUB));

    // 바닥값은 "가중치"에만 적용
    const double wSBs_floor = (wSBs <= floorW ? floorW : wSBs);
    const double wSUB_floor = (wSUB <= floorW ? floorW : wSUB);

    // 좌표로 "대입" (누적 아님)
    binScaleBKG->set(*argSIG, wSBs_floor, eSBs);
    binSubtrSIG->set(*argSIG, wSUB_floor, eSUB);
  }
}
// make simple hist
// RooDataSet *redDataSB = (RooDataSet *)dsReduced->reduce("mass < 2.9 || mass > 3.3");
// RooDataSet *redDataSIG = (RooDataSet *)dsReduced->reduce("mass > 2.9 && mass < 3.3");

// subtract
// RooDataHist *tbinDataCtErrSB = new RooDataHist("tbinDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
// RooDataHist *tbinDataCtErrSIG = new RooDataHist("tbinDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

// === finish define helper functions ===

// void getCtErrRange(RooDataSet *data, const char *t_reduceDS_woCtErr, float lmin, float lmax, float *errmin, float *errmax);
// RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName);

void s03_sideband()
{
  cout << "=== start s03_sideband() ===\n";
  TStopwatch time;
  time.Start();

  // set variables
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // observable ranges
  double ctMin = -0.5, ctMax = 2; // lmin, lmax: 1 for lowpT, 2 for higpT
  float errmin = 0.0, errmax = 1.0;

  // sideband ranges
  double sbL_lo=2.6, sbL_hi=2.9, sig_lo=2.9, sig_hi=3.3, sbR_lo=3.3, sbR_hi=3.5;

  // read inputs
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile inputData(fileNameData.c_str());
  RooDataSet *dataPRMC = (RooDataSet *)inputData.Get("dataset");
  dataPRMC->SetName("dataPRMC");

  // apply basic cuts
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(inputData.Get("dataset"));
  dsRaw->Print();
  if (!dsRaw)
  {
    cout << "[Error]: Cannot find RooDataSet\n";
    return;
  }

  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0, "weight");
  dsWeight->Print();

  // define cut
  char reduceDS_woCtErr[3000];
  sprintf(reduceDS_woCtErr, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f)"

          " && (recoQQsign == 0) ",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax);

  RooDataSet *dsReduced_tmp = (RooDataSet *)dsWeight->reduce(reduceDS_woCtErr);

  // observable
  auto mass = new RooRealVar("mass", "invariant mass", massLow, massHigh, "GeV/c^{2}");
  auto ctau3D = new RooRealVar("ctau3D", "", -0.5, 2, "#sigma_{L_{J/#psi}} [mm]");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", 0, 1, "#sigma_{L_{J/#psi}} [mm]");
  ctau3DErr->setBins(200);
  auto weight = new RooRealVar("weight", "", 0, 10000, "");
  RooArgSet obs(*mass, *ctau3DErr, *weight);
  auto dsReduced = new RooDataSet("dsReduced", "dataset with local vars", obs, Import(*dsReduced_tmp));
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced->Print();

  RooDataSet *redDataSIG_tmp = (RooDataSet *)dsReduced->reduce("mass>2.9 && mass<3.3");
  RooDataSet *redDataSB_tmp = (RooDataSet *)dsReduced->reduce("mass<2.9 || mass>3.3");

  // === start getCtErrRange ===
  // get raw mass fit result -> 그러면 여기서 mass fit 할 필요 X
  RooRealVar sl1("sl1","Cheb c1", 0.0, -1.0, 1.0);
  RooRealVar sl2("sl2","Cheb c2", 0.0, -1.0, 1.0);

  // --- 이전 mass-fit 결과에서 Chebyshev 계수만 로드+고정 ---
  {
    TFile fPrev("roots/raw_mass_fit.root", "READ"); // ★ 파일 경로/이름 맞춰주세요
    if (fPrev.IsZombie())
    {
      ::Error("cheb-load", "Cannot open previous mass-fit file.");
    }
    else
    {
      // 결과 객체 이름은 환경에 맞게 바꾸세요. (예: "fitMass", "fitMassBkg", ...)
      RooFitResult *fr = dynamic_cast<RooFitResult *>(fPrev.Get("fitMass"));
      if (!fr)
      {
        ::Error("cheb-load", "RooFitResult 'fitMass' not found in file.");
      }
      else
      {
        auto apply_one = [&](RooRealVar &dest)
        {
          // float된 리스트에서 먼저 찾고, 없으면 const 리스트에서 보충
          const RooArgList &fl = fr->floatParsFinal();
          const RooArgList &cl = fr->constPars();
          if (auto *a = fl.find(dest.GetName()))
          {
            dest.setVal(static_cast<const RooRealVar *>(a)->getVal());
            dest.setConstant(true);
            return;
          }
          if (auto *a = cl.find(dest.GetName()))
          {
            dest.setVal(static_cast<const RooRealVar *>(a)->getVal());
            dest.setConstant(true);
            return;
          }
          ::Warning("cheb-load", "Parameter '%s' not found in fit result; keep as-is.", dest.GetName());
        };

        // ★ mean/sigma/sigmaG는 건드리지 않고, Cheb 계수만 세팅+고정
        apply_one(sl1);
        apply_one(sl2);
        // (차수가 더 높으면 sl3..slN도 같은 방식으로 추가)
      }
    }
  }

  RooArgList allPars_local;
  allPars_local.add(sl1);
  allPars_local.add(sl2);
  double scaleF = computeScaleF_cheby(allPars_local, /*order=*/2, massLow, massHigh, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi);
  cout << "scaleF: " << scaleF << "\n";

  // --- subtract sig - bkg ---
  // original RooDataHist
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", *ctau3DErr, *redDataSB_tmp);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", *ctau3DErr, *redDataSIG_tmp);

  RooDataHist* binSubtrSIG = new RooDataHist("binSubtrSIG","SIG-SB", RooArgSet(*ctau3DErr));
  RooDataHist* binScaleBKG = new RooDataHist("binScaleBKG","scaled SB", RooArgSet(*ctau3DErr));

  // dk가 없으면 0.0
  subtractSidebands_inplace_propagate_TH1(binSubtrSIG, binScaleBKG,
                                          binDataCtErrSIG, binDataCtErrSB,
                                          ctau3DErr,
                                          scaleF, "ctau3DErr",
                                          /*dk=*/0.0, /*floorW=*/1e-6);

  // --- find out proper bin ranges ---
  // RooDataHist -> TH1 생성 (축 변수: *ctau3DErr)
  TH1 *histDataCtErrSIG = binDataCtErrSIG->createHistogram("histDataCtErrSIG", *ctau3DErr);
  TH1 *histSubtractedSIG = binSubtrSIG->createHistogram("histSubtractedSIG", *ctau3DErr);
  TH1 *histScaledBKG = binScaleBKG->createHistogram("histScaledBKG", *ctau3DErr);

  double minSig = 0.5, maxSig = 0.0, minBkg = 0.5, maxBkg = 0.0;
  double cutValue = 0.2;

  // 최대 빈 찾기
  int maxBinSig = histSubtractedSIG->GetMaximumBin();
  int maxBinBkg = histScaledBKG->GetMaximumBin();

  // 초기 경계값(최대 빈의 하한/상한)
  minSig = histSubtractedSIG->GetBinLowEdge(maxBinSig);
  minBkg = histScaledBKG->GetBinLowEdge(maxBinBkg);

  maxSig = histSubtractedSIG->GetBinLowEdge(maxBinSig + 1);
  maxBkg = histScaledBKG->GetBinLowEdge(maxBinBkg + 1);

  // 아래쪽으로 이동하며 cutValue 초과 구간의 최솟값 찾기
  for (int xbins = maxBinSig; xbins > 0; --xbins)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      minSig = histSubtractedSIG->GetBinLowEdge(xbins);
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins > 0; --xbins)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      minBkg = histScaledBKG->GetBinLowEdge(xbins);
    }
    else
      break;
  }

  // 위쪽으로 이동하며 cutValue 초과 구간의 최댓값 찾기
  for (int xbins = maxBinSig; xbins <= histSubtractedSIG->GetNbinsX(); ++xbins)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      maxSig = histSubtractedSIG->GetBinLowEdge(xbins + 1);
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins <= histScaledBKG->GetNbinsX(); ++xbins)
  { // ★ fix: maxBinBkg 사용
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      maxBkg = histScaledBKG->GetBinLowEdge(xbins + 1);
    }
    else
      break;
  }

  // 공통 범위 선택
  double tmpMin = (minSig > minBkg) ? minSig : minBkg;
  double tmpMax = (maxSig < maxBkg) ? maxSig : maxBkg;

  // 라운딩(하한: 내림, 상한: 올림)
  tmpMin = TMath::Floor(tmpMin * 1000.0) / 1000.0;
  tmpMax = TMath::Ceil(tmpMax * 1000.0) / 1000.0;

  // 데이터셋 리듀스 (ws 없이: ctau3DErr 이름 그대로 사용)
  char reduceDS[512];
  std::snprintf(reduceDS, sizeof(reduceDS), "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)dsReduced->reduce(reduceDS);

  // 범위가 너무 타이트하면(SIG 기준) SIG 범위로 대체
  if (redDataTmp->sumEntries() < dsReduced->sumEntries() * 0.9)
  {
    delete redDataTmp;
    std::snprintf(reduceDS, sizeof(reduceDS), "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)dsReduced->reduce(reduceDS);
    tmpMin = minSig;
    tmpMax = maxSig;
  }

  // 최소 폭 보장
  if ((tmpMax - tmpMin) < 0.008)
  {
    std::cout << "getCtErrRange:: Maximum is less than minimum! Possibly few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  // --- draw raw plots ---
  // --- 드로잉 (RooPlot 축 + TH1 오버레이) ---
  TCanvas c0("ctau_err", "ctau_err", 500, 500);
  c0.cd();
  c0.SetLogy(1);

  // ws 없이 프레임 생성
  RooPlot *errframe2 = ctau3DErr->frame();
  errframe2->SetTitle("");
  errframe2->GetXaxis()->SetTitle("ctau3DErr");
  errframe2->GetYaxis()->SetTitle("Counts");
  double maxY = std::max({histDataCtErrSIG->GetMaximum(),
                          histScaledBKG->GetMaximum(),
                          histSubtractedSIG->GetMaximum()});
  errframe2->SetMaximum(maxY * 1.3);
  errframe2->SetMinimum(0.2);
  errframe2->Draw();

  // TH1 스타일 & 오버레이
  histDataCtErrSIG->SetMarkerColor(kRed);
  histDataCtErrSIG->SetLineColor(kWhite);
  histDataCtErrSIG->SetMarkerStyle(24);
  histDataCtErrSIG->GetXaxis()->CenterTitle(1);
  histDataCtErrSIG->GetYaxis()->CenterTitle(1);
  histDataCtErrSIG->Draw("pe same");

  histScaledBKG->SetMarkerColor(kBlue);
  histScaledBKG->SetLineColor(kWhite);
  histScaledBKG->SetMarkerStyle(24);
  histScaledBKG->Draw("pe same");

  histSubtractedSIG->SetLineColor(kWhite);
  histSubtractedSIG->Draw("pe same");

  // 주석/라벨
  TLatex t;
  t.SetNDC();
  t.SetTextAlign(32);
  t.SetTextSize(0.04);
  t.SetTextColor(kRed);
  char comment[200];
  std::snprintf(comment, sizeof(comment), "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
  t.DrawLatex(0.92, 0.60, comment);
  t.SetTextColor(kBlack);

  // 범례
  TLegend legsb(0.60, 0.19, 0.90, 0.35, nullptr, "brNDC");
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

  // 저장
  std::string titlestr = "figs/ctau_err_raw.pdf";
  c0.SaveAs(titlestr.c_str());

  errmin = tmpMin;
  errmax = tmpMax;

  // === finish getCtErrRange ===


  // === start building dataset and punzi-term with ctau3DErr cut ===
  // --- make RooDataSet for the next processes ---
  ctau3DErr = new RooRealVar("ctau3DErr", "", errmin, errmax, "#sigma_{L_{J/#psi}} [mm]");
  // ctau3DErr->setBins(50);
  RooArgSet obs2(*mass, *ctau3D, *ctau3DErr, *weight);
  auto dsReduced2 = new RooDataSet("dsReduced2", "dataset with local vars", obs2, Import(*dsReduced_tmp));;

  char reduceDS_wCtErr[3000];
  sprintf(reduceDS_wCtErr, // no multiplicty case
        "ctau3DErr > %.3f && ctau3DErr < %.3f", errmin, errmax);
  auto dsReduced3 = (RooDataSet *)dsReduced2->reduce(reduceDS_wCtErr);
  RooDataSet *redDataSB = (RooDataSet *)dsReduced3->reduce("mass<2.9 || mass>3.3");

  // --- build RooHistFunc---
  // --- RooDataHist -> RooHistFunc -> RooGenericPdf ---
  // 사용 대상: binSubtrSIG (SIG − scale×SB), binScaleBKG (scaled SB)
  // 관측변수: *ctau3DErr  (반드시 RooRealVar*)

  RooArgList obs_cterr(*ctau3DErr);

  // (1) 히스토그램 함수를 생성 (보간 차수는 환경에 맞게 조절; 0=nearest, 1=linear 등 ROOT 버전에 따라 다를 수 있음)
  RooHistPdf sigPdf("sigPdf", "Signal hist func",
                    obs_cterr, *binSubtrSIG, 0);

  RooHistPdf bkgPdf("bkgPdf", "Background hist func",
                    obs_cterr, *binScaleBKG, 0);

  // (2) RooGenericPdf로 감싸서 PDF화 (@0은 ArgList의 첫 항목을 의미)
  //  - RooGenericPdf는 주어진 식을 적분해 정규화하므로, 관측변수의 range가
  //    히스토그램 범위와 일치하도록 ctau3DErr의 range를 맞춰두세요.
  // RooGenericPdf sigPdf("sigPdf", "@0", RooArgList(sigFunc));
  // RooGenericPdf bkgPdf("bkgPdf", "@0", RooArgList(bkgFunc));

  // --- draw RooHistFunc plots ---
  // --- Overlay + Pull: SIG ---
  {
    // (메인) 데이터+모델
    RooPlot *fSig = ctau3DErr->frame();
    binSubtrSIG->plotOn(fSig,
                        DataError(RooAbsData::SumW2),
                        MarkerColor(kRed), LineColor(kRed),
                        MarkerStyle(20), Name("h_sig"));
    // RooGenericPdf 기반 모델
    sigPdf.plotOn(fSig, LineColor(kRed + 2), LineWidth(2),
                  /*Range("histRange"), NormRange("histRange"),*/
                  Name("curve_sig"));

    // (풀) Data-Model / sigma
    RooHist *pullSig = fSig->pullHist("h_sig", "curve_sig", true);
    RooPlot *fSigPull = ctau3DErr->frame();
    fSigPull->SetTitle("Pull (SIG)");
    fSigPull->GetYaxis()->SetTitle("(Data - Model) / #sigma");
    fSigPull->SetMinimum(-5);
    fSigPull->SetMaximum(5);
    fSigPull->addPlotable(pullSig, "P");

    // (캔버스) 1×2: 위-메인, 아래-풀
    TCanvas cSig("cSig", "ctau3DErr: SIG (data vs model + pull)", 800, 900);
    cSig.Divide(1, 2);

    cSig.cd(1);
    fSig->SetTitle("");
    fSig->GetXaxis()->SetTitle("ctau3DErr");
    fSig->GetYaxis()->SetTitle("Counts");
    fSig->Draw();
    // 범례
    {
      TLegend leg(0.60, 0.70, 0.88, 0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.AddEntry(fSig->findObject("h_sig"), "SIG data", "p");
      leg.AddEntry(fSig->findObject("curve_sig"), "SIG model", "l");
      leg.Draw();
    }

    cSig.cd(2);
    fSigPull->Draw();

    cSig.SaveAs("figs/ctErr_SIG_data_vs_histfunc_withPull.pdf");

    // (옵션) 눈대중 chi2/ndf
    // std::cout << "[SIG] chi2/ndf(vis) = " << fSig->chiSquare("curve_sig","h_sig") << "\n";
  }

  // --- Overlay + Pull: BKG ---
  {
    RooPlot *fBkg = ctau3DErr->frame();
    binScaleBKG->plotOn(fBkg,
                        DataError(RooAbsData::SumW2),
                        MarkerColor(kBlue), LineColor(kBlue),
                        MarkerStyle(24), Name("h_bkg"));
    bkgPdf.plotOn(fBkg,
                  LineColor(kBlue + 2), LineWidth(2),
                  /*Range("histRange"), NormRange("histRange"),*/
                  Name("curve_bkg"));

    RooHist *pullBkg = fBkg->pullHist("h_bkg", "curve_bkg", true);
    RooPlot *fBkgPull = ctau3DErr->frame();
    fBkgPull->SetTitle("Pull (BKG)");
    fBkgPull->GetYaxis()->SetTitle("(Data - Model) / #sigma");
    fBkgPull->SetMinimum(-5);
    fBkgPull->SetMaximum(5);
    fBkgPull->addPlotable(pullBkg, "P");

    TCanvas cBkg("cBkg", "ctau3DErr: BKG (data vs model + pull)", 800, 900);
    cBkg.Divide(1, 2);

    cBkg.cd(1);
    fBkg->SetTitle("");
    fBkg->GetXaxis()->SetTitle("ctau3DErr");
    fBkg->GetYaxis()->SetTitle("Counts");
    fBkg->Draw();
    {
      TLegend leg2(0.60, 0.70, 0.88, 0.88);
      leg2.SetBorderSize(0);
      leg2.SetFillStyle(0);
      leg2.AddEntry(fBkg->findObject("h_bkg"), "BKG data", "p");
      leg2.AddEntry(fBkg->findObject("curve_bkg"), "BKG model", "l");
      leg2.Draw();
    }

    cBkg.cd(2);
    fBkgPull->Draw();

    cBkg.SaveAs("figs/ctErr_BKG_data_vs_histfunc_withPull.pdf");

    // std::cout << "[BKG] chi2/ndf(vis) = " << fBkg->chiSquare("curve_bkg","h_bkg") << "\n";
  }

  // === Finish building dataset and punzi-term with ctau3DErr cut ===

  // save result
  // --- save observable ranges as RooValues ---
  TFile output("roots/ctau_err.root", "recreate");
  sigPdf.Write();
  bkgPdf.Write();
  binSubtrSIG->Write("binSubtrSIG");
  dsReduced3->Write("dsReduced");
  output.Close();
  // mass range, ctau range, ctauErr range

  cout << "=== finish s03_sideband() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}