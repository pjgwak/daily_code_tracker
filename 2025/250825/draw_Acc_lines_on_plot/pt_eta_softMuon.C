// example codes by ChatGPT5
// pt_eta_softMuon.C
// ROOT 6.x에서 컴파일 없이 바로 실행 가능: root -l -q 'pt_eta_softMuon.C()'
#include <TCanvas.h>
#include <TRandom3.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <cmath>

// --- Soft-muon 풍의 pT 최소 커브(대략적인 데모용) ---
// 그림처럼: 중앙(|eta|<1.2)에서는 상수(3.3 GeV) → 약간 상승 → 전진 시 전방(|eta|~2.1)으로 갈수록 1.0 GeV까지 감소
double pTmin_soft(double eta)
{
  double aeta = std::fabs(eta);
  if (aeta > 2.4)
    return 1e9; // 수용 범위 밖(그릴 때 위로 튀지 않게 큰 값)
  if (aeta < 1.2)
    return 3.3; // 중앙 평탄
  if (aeta < 1.6)
  { // 약간 상승: 3.3 → 3.5
    double t = (aeta - 1.2) / (1.6 - 1.2);
    return 3.3 + t * (3.5 - 3.3);
  }
  if (aeta < 2.1)
  { // 완만한 하강: 3.5 → 1.0
    double t = (aeta - 1.6) / (2.1 - 1.6);
    return 3.5 + t * (1.0 - 3.5);
  }
  return 1.0; // 전방(2.1~2.4)에서는 낮은 임계
}

// pT 스펙트럼(토이): f(pT) ~ pT * exp(-pT/p0), p0 ~ 1.2
double sample_pT(TRandom3 &R, double p0 = 1.2, double pTmax = 5.0)
{
  // 간단한 역변환 샘플링 대신 수용거부법(accept-reject)
  // g(pT)=C*exp(-pT/p0) 상한으로 쓰고, pT factor는 u< (pT/p0max) 로 보정
  // 여기선 보수적으로 accept-reject 두 번 루프 돌지만 평균 몇 번 내에 수렴함
  while (true)
  {
    double x = R.Exp(p0); // Exp(mean=p0), [0, ∞)
    if (x > pTmax)
      continue;
    double u = R.Uniform();
    // 상한 비: pT factor를 (x/pTmax)로 단순 포함(대충 맞춰도 토이에는 충분)
    if (u < x / pTmax)
      return x;
  }
}

// eta 분포(토이): 중앙이 두터운 가우시안 + 컷(|eta|<2.4)
double sample_eta(TRandom3 &R)
{
  while (true)
  {
    double e = R.Gaus(0.0, 1.0);
    if (std::fabs(e) <= 2.4)
      return e;
  }
}

void pt_eta_softMuon()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird); // 보기 좋은 팔레트
  gStyle->SetNumberContours(50);

  const int N = 200000; // 토이 MC 이벤트 수
  const double etamin = -2.5, etamax = 2.5;
  const double ptmin = 0.5, ptmax = 5.0;

  TH2D *h2 = new TH2D("h2", ";#eta; p_{T}  [GeV]", 100, etamin, etamax, 100, ptmin, ptmax);

  TRandom3 R(0xA713);
  // --- 토이 생성 ---
  for (int i = 0; i < N; ++i)
  {
    double eta = sample_eta(R);
    double pt = sample_pT(R, 1.2, ptmax);
    // 전역 분포(수용 전)를 채움: 이미지처럼 바탕 히트맵 느낌
    h2->Fill(eta, pt);
  }

  // --- 캔버스 & 스타일 ---
  TCanvas *c = new TCanvas("c", "pT vs eta with Soft Muon Acc.", 900, 640);
  c->SetRightMargin(0.15);
  c->SetLogz(); // Z축 로그로 색 대비 강화

  // 2D 히스토그램 드로우
  h2->SetMinimum(1); // logz에서 0 방지
  h2->Draw("COLZ");

  // --- 수용 경계선 그리기 (점선) ---
  const int Np = 400;
  TGraph *gAcc = new TGraph(Np);
  // 좌→우로 eta를 스캔하며 pTmin(eta) 그리기
  for (int i = 0; i < Np; ++i)
  {
    double eta = etamin + (etamax - etamin) * (i / (double)(Np - 1));
    double y = pTmin_soft(eta);
    // y가 y범위를 넘으면 캔버스 밖으로 튀지 않게 살짝 클리핑
    if (y < ptmin)
      y = ptmin;
    if (y > ptmax)
      y = ptmax;
    gAcc->SetPoint(i, eta, y);
  }
  gAcc->SetLineColor(kBlack);
  gAcc->SetLineStyle(2); // dotted
  gAcc->SetLineWidth(3);
  gAcc->Draw("L SAME");

  // (선택) 경계선의 외곽 느낌을 주는 솔리드 라인 하나 더
  TGraph *gAccSolid = (TGraph *)gAcc->Clone("gAccSolid");
  gAccSolid->SetLineStyle(1);
  gAccSolid->SetLineColor(kGray + 2);
  gAccSolid->SetLineWidth(2);
  gAccSolid->Draw("L SAME");

  // --- 라벨/범례 ---
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.17, 0.87, "Soft Muon's Acc.");
  // 작은 범례 대신 선 스타일 안내
  TLegend *leg = new TLegend(0.62, 0.82, 0.88, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gAcc, "Acceptance boundary", "l");
  leg->Draw();

  c->Update();
  c->SaveAs("pt_eta_softMuon.png");
}
