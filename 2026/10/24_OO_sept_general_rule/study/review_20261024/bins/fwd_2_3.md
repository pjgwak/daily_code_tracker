# Forward, 2 < pT < 3 GeV/c, 1.6 < |y| < 2.4

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md) · [적용한 정책](../FIT_POLICY.md)

B = **0.117343 ± 0.010307**. 기준값 대비 상대 차이 **+0.66127%**. 오차는 고정 입력과 선택된 모형을 조건으로 한 fit covariance에서 얻었다.

수렴 status=0, covQual=3, EDM=9.09e-08. ctau 적합 범위 [-4.0, 6.5] mm.

## 적용 모형과 매개변수

PR Gaussian 성분 수: 4. NP 지수 성분 수: 2.

PR 평균·폭·첫 재귀 구성비를 적합하고 나머지 재귀 구성비 두 개는 MC 값에 고정한다. 독립 폭 제약은 적용하지 않는다. 구성비 공분산 배율은 16이며 두 폭 간격비 제약을 적용한다. NP 평균 수명과 총 yield에 MC 모멘트 제약을 적용하고 두 수명은 모두 적합한다.

## 해석 시 확인할 사항

- 두 NP 수명은 약 0.42418 mm로 일치한다. 선택한 결과는 NP MC 초기값에서 얻었으며 수명 고정은 적용하지 않았다.
- 성공한 실행의 초기값·초기 오차·범위·고정 상태를 동일하게 복원했을 때 재현됐다. 최종 추정값을 새로운 초기값으로 사용한 실행에서는 F2의 수렴 실패와 M8의 상대 차이 +2.49559%를 관찰했다.

[최종 ROOT](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/fit2d_result_y1.60_2.40_pt2.00_3.00.root) · [초기값·범위](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt2.00_3.00/ACTUAL_INITIAL_START_0.tsv) · [최종 변수](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt2.00_3.00/FINAL_OBJECTIVE_PARAMETERS.tsv) · [실제 제약 항](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt2.00_3.00/EXTERNAL_CONSTRAINTS.txt)

## Prefit 적합도

| Prefit | χ²/ndf | status/covQual | EDM | Equal-population χ²/ndf |
|---|---:|:---:|---:|---:|
| MC_MASS | 0.7898 | 0/3 | 0.000228 | — |
| PR | 0.9512 | 0/3 | 0.000551 | — |
| NP | 0.8455 | 0/3 | 9.76e-10 | 1.2453 (ndf 96) |
| MASS | 0.8899 | 0/3 | 0.00095 | — |
| BKG | 0.6695 | 0/3 | 0.000304 | 0.9011 (ndf 84) |

원래 저장된 binning의 χ²/ndf와 equal-population 진단을 구별한다. 후자는 추가 binning 진단이며 최소값을 골라 통과 여부를 바꾸지 않았다.

## MC mass prefit

![MC mass prefit](../figures/prefits/fwd_2_3_mc_mass.png)

[PDF](../figures/prefits/fwd_2_3_mc_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/mc_mass/mc_mass_model_y1.60_2.40_pt2.00_3.00.root)

## Prompt MC decay-length prefit

![Prompt MC decay-length prefit](../figures/prefits/fwd_2_3_pr.png)

[PDF](../figures/prefits/fwd_2_3_pr.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/ctau_pr/ctau_resolution_y1.60_2.40_pt2.00_3.00.root)

## Nonprompt MC truth prefit

![Nonprompt MC truth prefit](../figures/prefits/fwd_2_3_np.png)

[PDF](../figures/prefits/fwd_2_3_np.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_np/roots/y1.60_2.40/ctau_np/ctau_np_model_y1.60_2.40_pt2.00_3.00.root)

## Data mass prefit

![Data mass prefit](../figures/prefits/fwd_2_3_mass.png)

[PDF](../figures/prefits/fwd_2_3_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/mass/mass_model_y1.60_2.40_pt2.00_3.00.root)

## Data sideband decay-length prefit

![Data sideband decay-length prefit](../figures/prefits/fwd_2_3_bkg.png)

[PDF](../figures/prefits/fwd_2_3_bkg.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/ctau_bkg/ctau_bkg_fitresult_y1.60_2.40_pt2.00_3.00.root)

## Invariant-mass projection

χ²/ndf = 79.075/98 = 0.8069; diagnostic p-value = 0.91949.

![Invariant-mass projection](../figures/projections/fwd_2_3_mass.png)

[PDF](../figures/projections/fwd_2_3_mass.pdf)

## Decay-length projection

χ²/ndf = 279.342/355 = 0.7869; diagnostic p-value = 0.99884.

![Decay-length projection](../figures/projections/fwd_2_3_lifetime.png)

[PDF](../figures/projections/fwd_2_3_lifetime.pdf)

## Data-fit PDF boundary review

회색 점은 전체 mass 범위의 데이터이며 빨강·파랑 곡선은 각각 피팅된 PR·NP 신호 성분이다. 배경이 포함된 데이터와 신호 성분만의 높이를 직접적인 전체 적합도 검사로 해석하지 않는다. 검은 점선은 실제 적합 경계, 초록 점선은 진단 코드가 제안한 경계이며 적용하지 않았다.

![Data-fit PR and NP boundaries](../figures/boundaries/fwd_2_3.png)

[PDF](../figures/boundaries/fwd_2_3.pdf) · [경계 수치](../figures/boundaries/windows.tsv)

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md)
