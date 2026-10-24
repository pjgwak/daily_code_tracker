# Mid, 8 < pT < 9 GeV/c, |y| < 1.6

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md) · [적용한 정책](../FIT_POLICY.md)

B = **0.271284 ± 0.009935**. 기준값 대비 상대 차이 **+2.09766%**. 오차는 고정 입력과 선택된 모형을 조건으로 한 fit covariance에서 얻었다.

수렴 status=0, covQual=3, EDM=9.67e-07. ctau 적합 범위 [-1.0, 4.0] mm.

## 적용 모형과 매개변수

PR Gaussian 성분 수: 4. NP 지수 성분 수: 2.

PR 평균·폭·구성비에 공동 Gaussian 제약과 폭 간격비 제약을 적용한다. NP 수명·구성비는 MC 제약 없이 적합한다. 7–8 및 9–10과 별도의 설정이다.

## 해석 시 확인할 사항

- NP 혼합모형의 한 성분 비중이 0 경계에 접근한다. 해당 성분의 수명을 독립적으로 결정한 결과로 해석하지 않는다.
- 성공한 실행의 초기값·초기 오차·범위·고정 상태를 동일하게 복원했을 때 재현됐다. 최종 추정값을 새로운 초기값으로 사용한 실행에서는 F2의 수렴 실패와 M8의 상대 차이 +2.49559%를 관찰했다.

[최종 ROOT](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/fit2d_result_y0.00_1.60_pt8.00_9.00.root) · [초기값·범위](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt8.00_9.00/ACTUAL_INITIAL_START_0.tsv) · [최종 변수](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt8.00_9.00/FINAL_OBJECTIVE_PARAMETERS.tsv) · [실제 제약 항](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y0.00_1.60/fit2d/diagnostics_y0.00_1.60_pt8.00_9.00/EXTERNAL_CONSTRAINTS.txt)

## Prefit 적합도

| Prefit | χ²/ndf | status/covQual | EDM | Equal-population χ²/ndf |
|---|---:|:---:|---:|---:|
| MC_MASS | 0.8881 | 0/3 | 0.000466 | — |
| PR | 0.8916 | 0/3 | 0.000431 | — |
| NP | 0.8050 | 0/3 | 0.000406 | 0.9226 (ndf 96) |
| MASS | 0.6354 | 0/3 | 0.000225 | — |
| BKG | 0.8904 | 0/3 | 3.46e-05 | 2.6388 (ndf 5) |

원래 저장된 binning의 χ²/ndf와 equal-population 진단을 구별한다. 후자는 추가 binning 진단이며 최소값을 골라 통과 여부를 바꾸지 않았다.

## MC mass prefit

![MC mass prefit](../figures/prefits/mid_8_9_mc_mass.png)

[PDF](../figures/prefits/mid_8_9_mc_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y0.00_1.60/mc_mass/mc_mass_model_y0.00_1.60_pt8.00_9.00.root)

## Prompt MC decay-length prefit

![Prompt MC decay-length prefit](../figures/prefits/mid_8_9_pr.png)

[PDF](../figures/prefits/mid_8_9_pr.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y0.00_1.60/ctau_pr/ctau_resolution_y0.00_1.60_pt8.00_9.00.root)

## Nonprompt MC truth prefit

![Nonprompt MC truth prefit](../figures/prefits/mid_8_9_np.png)

[PDF](../figures/prefits/mid_8_9_np.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_np/roots/y0.00_1.60/ctau_np/ctau_np_model_y0.00_1.60_pt8.00_9.00.root)

## Data mass prefit

![Data mass prefit](../figures/prefits/mid_8_9_mass.png)

[PDF](../figures/prefits/mid_8_9_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y0.00_1.60/mass/mass_model_y0.00_1.60_pt8.00_9.00.root)

## Data sideband decay-length prefit

![Data sideband decay-length prefit](../figures/prefits/mid_8_9_bkg.png)

[PDF](../figures/prefits/mid_8_9_bkg.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y0.00_1.60/ctau_bkg/ctau_bkg_fitresult_y0.00_1.60_pt8.00_9.00.root)

## Invariant-mass projection

χ²/ndf = 65.588/98 = 0.6693; diagnostic p-value = 0.99514.

![Invariant-mass projection](../figures/projections/mid_8_9_mass.png)

[PDF](../figures/projections/mid_8_9_mass.pdf)

## Decay-length projection

χ²/ndf = 136.423/157 = 0.8689; diagnostic p-value = 0.88064.

![Decay-length projection](../figures/projections/mid_8_9_lifetime.png)

[PDF](../figures/projections/mid_8_9_lifetime.pdf)

## Data-fit PDF boundary review

회색 점은 전체 mass 범위의 데이터이며 빨강·파랑 곡선은 각각 피팅된 PR·NP 신호 성분이다. 배경이 포함된 데이터와 신호 성분만의 높이를 직접적인 전체 적합도 검사로 해석하지 않는다. 검은 점선은 실제 적합 경계, 초록 점선은 진단 코드가 제안한 경계이며 적용하지 않았다.

![Data-fit PR and NP boundaries](../figures/boundaries/mid_8_9.png)

[PDF](../figures/boundaries/mid_8_9.pdf) · [경계 수치](../figures/boundaries/windows.tsv)

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md)
