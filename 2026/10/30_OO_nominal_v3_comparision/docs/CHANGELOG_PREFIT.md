# prefit 설정 변경 (28일 FINAL_V2 → 31일 FINAL_V3 → FINAL_V4)

`config/prefit_settings.txt`. 이전 파일은 `config/prefit_settings_FINAL_V2.txt`, 이전 prefit ROOT는 `roots_FINAL_V2_backup/`.

| 빈 | 변경 | 근거 |
|---|---|---|
| fwd 1–2 | ctau 창 [−4, 6.5] → [−6, 10] mm (mass·bkg prefit 재적합, `bkg_saved_ct` 갱신) | 규칙(PR 자유, NP 2성분 MC 제약, 배경 자유)을 고정하고 창만 바꾸면 B가 연속적으로 움직인다: [−4, 6.5] +2.3%, [−6, 10] −0.01%, [−8, 12] −1.2%. 배경 prefit χ²/ndf 1.07, 최종 ctau p 0.97. REPORT 5장 |
| fwd 2–3 | 질량 창 상한 3.5 → 3.4 GeV(오전 FINAL_V3) → **되돌림**(오후, 2.6–3.5, FINAL_V2 prefit 복원). 대신 PR MC fit 범위 분위수 `pr_quant` 0.0005/0.9995 → 0.001/0.999 (pr·np·bkg prefit 재적합) | 질량 창을 빈마다 다르게 두지 않기로 함(사용자). 분위수 0.001 + 구성비 MC 공분산 제약 → +1.91%, PR prefit χ²/ndf 1.31. 분위수 의존은 0.001에서 국소 최소(REPORT 9장) |

나머지 18빈의 prefit 설정은 28일과 같고(mid 4빈은 2D fit 규칙만 바뀜), 다시 돌린 prefit 값도 같다(최종 B가 스캔 값과 일치).
