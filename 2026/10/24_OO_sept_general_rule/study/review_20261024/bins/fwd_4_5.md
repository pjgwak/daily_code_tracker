# Forward, 4 < pT < 5 GeV/c, 1.6 < |y| < 2.4

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md) · [적용한 정책](../FIT_POLICY.md)

B = **0.175474 ± 0.009036**. 기준값 대비 상대 차이 **+1.91545%**. 오차는 고정 입력과 선택된 모형을 조건으로 한 fit covariance에서 얻었다.

수렴 status=0, covQual=3, EDM=1.72e-08. ctau 적합 범위 [-3.5, 4.5] mm.

## 적용 모형과 매개변수

PR Gaussian 성분 수: 4. NP 지수 성분 수: 2.

PR 평균·폭·구성비를 공동 Gaussian 제약 아래 적합한다. MC 상관행렬을 유지하며 표준편차는 정책표에 기록한다. NP에는 MC 공분산을 적용하고, 배경 수명·yield에는 sideband prefit 제약을 적용한다.

[최종 ROOT](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/fit2d_result_y1.60_2.40_pt4.00_5.00.root) · [초기값·범위](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt4.00_5.00/ACTUAL_INITIAL_START_0.txt) · [최종 변수](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt4.00_5.00/FINAL_OBJECTIVE_PARAMETERS.txt) · [실제 제약 항](../../23_OO_simple_rules_2p1/results/recovery_candidate20_exact_starts_g7/D7_Recovery_candidate20/roots/y1.60_2.40/fit2d/diagnostics_y1.60_2.40_pt4.00_5.00/EXTERNAL_CONSTRAINTS.txt)

## Prefit 적합도

| Prefit | χ²/ndf | status/covQual | EDM | Equal-population χ²/ndf |
|---|---:|:---:|---:|---:|
| MC_MASS | 0.9635 | 0/3 | 6.93e-05 | — |
| PR | 1.3523 | 0/3 | 0.000969 | — |
| NP | 1.5806 | 0/3 | 0.000455 | 0.9997 (ndf 96) |
| MASS | 1.0385 | 0/3 | 5.9e-05 | — |
| BKG | 1.3161 | 0/3 | 0.000435 | 1.8114 (ndf 88) |

원래 저장된 binning의 χ²/ndf와 equal-population 진단을 구별한다. 후자는 추가 binning 진단이며 최소값을 골라 통과 여부를 바꾸지 않았다.

## MC mass prefit

![MC mass prefit](../figures/prefits/fwd_4_5_mc_mass.png)

[PDF](../figures/prefits/fwd_4_5_mc_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/mc_mass/mc_mass_model_y1.60_2.40_pt4.00_5.00.root)

## Prompt MC decay-length prefit

![Prompt MC decay-length prefit](../figures/prefits/fwd_4_5_pr.png)

[PDF](../figures/prefits/fwd_4_5_pr.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/roots/y1.60_2.40/ctau_pr/ctau_resolution_y1.60_2.40_pt4.00_5.00.root)

## Nonprompt MC truth prefit

![Nonprompt MC truth prefit](../figures/prefits/fwd_4_5_np.png)

[PDF](../figures/prefits/fwd_4_5_np.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_np/roots/y1.60_2.40/ctau_np/ctau_np_model_y1.60_2.40_pt4.00_5.00.root)

## Data mass prefit

![Data mass prefit](../figures/prefits/fwd_4_5_mass.png)

[PDF](../figures/prefits/fwd_4_5_mass.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/mass/mass_model_y1.60_2.40_pt4.00_5.00.root)

## Data sideband decay-length prefit

![Data sideband decay-length prefit](../figures/prefits/fwd_4_5_bkg.png)

[PDF](../figures/prefits/fwd_4_5_bkg.pdf) · [Prefit ROOT](../../23_OO_simple_rules_2p1/results/calibrated_all20_data/roots/y1.60_2.40/ctau_bkg/ctau_bkg_fitresult_y1.60_2.40_pt4.00_5.00.root)

## Invariant-mass projection

χ²/ndf = 90.163/98 = 0.9200; diagnostic p-value = 0.70107.

![Invariant-mass projection](../figures/projections/fwd_4_5_mass.png)

[PDF](../figures/projections/fwd_4_5_mass.pdf)

## Decay-length projection

χ²/ndf = 215.576/229 = 0.9414; diagnostic p-value = 0.72863.

![Decay-length projection](../figures/projections/fwd_4_5_lifetime.png)

[PDF](../figures/projections/fwd_4_5_lifetime.pdf)

## Data-fit PDF boundary review

회색 점은 전체 mass 범위의 데이터이며 빨강·파랑 곡선은 각각 피팅된 PR·NP 신호 성분이다. 배경이 포함된 데이터와 신호 성분만의 높이를 직접적인 전체 적합도 검사로 해석하지 않는다. 검은 점선은 실제 적합 경계, 초록 점선은 진단 코드가 제안한 경계이며 적용하지 않았다.

![Data-fit PR and NP boundaries](../figures/boundaries/fwd_4_5.png)

[PDF](../figures/boundaries/fwd_4_5.pdf) · [경계 수치](../figures/boundaries/windows.txt)

[전체 빈 목차](../FIT_PLOTS_BY_BIN.md)
