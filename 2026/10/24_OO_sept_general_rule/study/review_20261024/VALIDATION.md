# 적합 결과 검증 기록

검토 대상은 `recovery_candidate20_exact_starts_g7 / D7_Recovery_candidate20`이다. 20개 빈의 결과는 이 실행에서 실제로 다시 적합한 값이다. [빈별 수치·그림](FIT_PLOTS_BY_BIN.md)과 [원본 결과 CSV](data/fit_results.tsv)를 함께 확인할 수 있다.

## Prefit

최종 적합에 사용한 MC mass·PR·NP 및 DATA mass·sideband prefit은 총 100개이다. 저장된 결과에서 모두 minimizer status=0, covQual=3, EDM<0.001 및 원래 binning의 χ²/ndf<2를 확인했다. 원본 ROOT의 SHA256을 기록하고 현재 사용한 파일과 대조했다.

| Sideband 구간 | Equal-population χ²/ndf | ndf |
|---|---:|---:|
| Forward 10–12 | 4.086 | 6 |
| Mid 7–8 | 4.159 | 21 |
| Mid 8–9 | 2.639 | 5 |
| Mid 9–10 | 2.079 | 4 |

이 별도 binning 진단은 원래 prefit의 χ²/ndf와 구분한다. 적은 빈 수와 비어 있는 영역의 영향을 검토할 때 두 값을 함께 사용하며, 가장 작은 χ²/ndf만 골라 판단하지 않는다. [100개 prefit 수치·해시](data/prefit_quality.tsv).

Prefit 그림 59개는 저장된 파라미터로 다시 출력한 것이며 재적합한 결과가 아니다. 저장된 χ²를 재현한 기록과 원본 그림의 출처를 [렌더링 검증 기록](audit/prefit_rendering_original.json)에 보존했다. 이 파일에 포함된 당시 경로는 [이동 기록](audit/cleanup_plan.json)으로 추적할 수 있다.

## 최종 fit과 기준값 비교

20개 결과 모두 수렴·공분산·EDM 및 mass·ctau projection 검사를 통과했다. p-value는 각각의 projection에 대한 진단값이며 전체 2D 분포의 단일 p-value가 아니다. 두 p-value 각각에 0.05 조건을 추가하지 않았다.

기준값에 대한 상대 차이는 100×(B/Bref−1)이며 전체 20빈이 ±2.1% 범위에 있다. 최대 절대 차이는 Mid 8–9의 **2.097661%**이다. 기준값 비교는 projection 적합도와 별도로 기록한다. B는 [0,1] 범위에서 적합하며 기준 B를 중심으로 한 제약이나 사후 환산을 적용하지 않는다.

PR 좌표는 104개가 자유이고 48개가 MC 값에 고정되어 있다. 사용하는 NP 좌표 71개는 모두 자유이며 NP 수명 고정은 없다. [좌표별 출처](data/fit_coordinates.tsv)와 [집계](audit/coordinate_provenance.json)를 보존했다. Mass·sideband prefit의 고정 입력은 이 PR·NP 좌표 개수와 구별한다.

## DATA PDF 경계

선택된 ROOT 결과에서 RooFit PDF를 재생하고 독립적인 Gaussian·exponential convolution 계산과 대조했다. **20/20**이 일치 검사를 통과했다. 실제 적용 범위에서 양쪽 경계의 기대 event 수와 범위 밖 꼬리 확률을 계산하고 모든 패널을 시각적으로 확인했다.

검토에 사용한 수치 보조 기준은 경계에서 0.02 mm당 기대 event 수 <0.01, 각 외부 꼬리 확률 ≤0.0005이다. [경계 CSV](figures/boundaries/windows.tsv)의 `old_low`·`old_high`가 실제 적용한 범위이며 `ct_low`·`ct_high`는 진단 코드의 제안 범위이다. 제안 범위로 다시 잘라서 결과를 바꾸지 않았다.

[독립 PDF 대조](figures/boundaries/replay_verification.json) · [전체 그림 PDF](figures/boundaries/all_bins.pdf).

## 초기값과 성분 분리

| 비교 실행 | Forward 2–3 | Mid 8–9 |
|---|---|---|
| 성공 fit의 최종값을 새로운 시작값으로 사용 | 수렴 실패 | +2.495585% |
| 성공 실행에 실제 입력했던 초기값을 복원 | +0.661265%, 수렴·projection 통과 | +2.097661%, 수렴·projection 통과 |

선택한 실행은 두 번째 경우이며 이 조건에서 20빈을 실제로 재실행했다. Forward 2–3의 NP 초기값은 MC prefit에서 가져온다. 핵심 네 빈의 초기값·초기 오차·범위·고정 상태를 원래 성공 실행과 비교한 [검산 기록](audit/initial_coordinates.json)을 보존했다. 성공점 재현과 임의의 초기값에서의 유일한 해 수렴은 구별한다.

- **Forward 1–2:** NP3의 세 번째 성분 비중은 약 2.8×10⁻¹⁰이다. 이 성분의 수명을 별도 측정값으로 해석하지 않는다.
- **Forward 2–3:** 두 NP 수명은 약 0.42418 mm로 일치한다.
- **Mid 8–9:** 한 NP 성분의 fit-window 비중은 약 1.1×10⁻¹¹이다. NP 제약을 제거했지만 개별 두 수명이 분리 측정된 것은 아니다.

[초기값 비교 원본 수치](data/start_comparison.json) · [피팅된 성분](data/pdf_components.tsv). 구성비의 `full`과 `fit_window` normalization 표기를 구분한다.

## 오차 해석

제시한 fit error는 선택된 모형·제약 및 고정된 prefit 입력을 조건으로 한 covariance 오차이다. 고정 입력, 모형 선택 및 제약 변화에 따른 불확실성은 별도 평가 대상이다. 현재 설정에서 sWeight를 사용한 신호 형상 고정은 없으며, 이를 성분별 mass–ctau 독립성 검증으로 해석하지 않는다.
