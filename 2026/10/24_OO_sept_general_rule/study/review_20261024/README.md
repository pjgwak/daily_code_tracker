# OO J/ψ 피팅 리뷰

현재 채택한 전체 실행 `recovery_candidate20_exact_starts_g7 / D7_Recovery_candidate20`을 기준으로 방법 서술, 설정표, 검증 수치, 빈별 그림을 정리했다. 원본 ROOT 결과와 수치 설정은 변경하지 않았다.

| 문서 | 내용 |
|---|---|
| [PAPER_METHODS.md](PAPER_METHODS.md) | 논문용 영문 방법 서술과 그림 설명문 |
| [FIT_POLICY.md](FIT_POLICY.md) | 구간별 고정·자유 변수, Gaussian 제약, 범위와 초기값 |
| [VALIDATION.md](VALIDATION.md) | Prefit·최종 적합도, PDF 경계 검증, 초기값 비교 |
| [FIT_PLOTS_BY_BIN.md](FIT_PLOTS_BY_BIN.md) | 20개 빈의 그림 검토 목차 |
| [CLEANUP.md](CLEANUP.md) | 폴더 정리 내역과 이력 보관 위치 |

![기준 결과와의 비교](figures/summary/agreement.png)

20빈 모두 수렴·projection 검사와 ±2.1% 기준값 비교를 통과했고, DATA fit의 PR·NP PDF 경계를 수치·시각적으로 확인했다. 이 기준값 비교와 데이터 분포의 적합도 검사는 구별한다. NP 수명을 고정한 빈은 없다.

삼각형은 단순화한 설정을 적용한 Forward 1–2, 2–3 및 Mid 7–8, 8–9이다. 보라색은 NP 개별 성분의 분리에 한계가 있는 Forward 1–2, 2–3 및 Mid 8–9이다. Forward 2–3과 Mid 8–9에는 초기값 의존성도 기록했다. 표시가 성공·실패를 대신하는 것은 아니며 세부 내용은 각 빈 문서에 있다.

각 빈 문서에는 실제 사용한 **prefit 5개 + 최종 projection 2개 + DATA PDF 경계 1개**가 포함된다. 그림은 PNG로 표시하고 원본 PDF와 ROOT 파일에 연결했다. Mass·ctau projection p-value 각각에 별도의 0.05 조건을 요구하지 않는다.

[최종 수치 CSV](data/fit_results.csv) · [Prefit 수치 CSV](data/prefit_quality.csv) · [설정 JSON](data/policy_by_bin.json) · [요약 PDF](figures/summary/agreement.pdf) · [전체 PDF 경계 그림](figures/boundaries/all_bins.pdf)

분석 코드와 재현 방법은 [23일 폴더 README](../23_OO_simple_rules_2p1/README.md)에 있다. 이 리뷰 폴더의 문서와 그림은 독립적으로 읽을 수 있으며 ROOT와 원본 진단 CSV 링크는 23일 폴더를 참조한다.
