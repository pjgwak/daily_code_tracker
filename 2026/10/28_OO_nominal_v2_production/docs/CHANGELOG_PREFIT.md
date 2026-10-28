# prefit 선택 변경 기록 (9월 단일 뮤온 acceptance 입력)

| 날짜 (파리) | 빈 | 변경 | 근거 (`docs/evidence/prefit_variants.txt`, `fwd23_effnp_variants.txt`) |
|---|---|---|---|
| 2026-10-25 00:30 | fwd 10–12, 12–14, 14–20 | NP 성분 2 → 1 (`np_count`) | 두 번째 NP 성분 퇴화: 10–12 yield 62±72, 12–14 오차>값, 14–20 수명 하한. 1성분 χ²/ndf 0.81, 1.02 |
| 2026-10-25 00:30 | fwd 5–6 | PR MC fit 범위 0.05–99.95% → 0.2–99.8% 분위수 (`pr_quant`) | PR χ²/ndf 1.27→1.03, 배경 0.74→0.59, 가장 넓은 Gauss 0.25→0.14 mm (24일 fwd 4–5와 같은 처리) |
| 2026-10-25 19:30 | fwd 2–3 | ctau 창 [−4, 6.5] → [−20, 20] mm, 배경 성분 2/1/1 → 2/0/1, `bkg_wide_tau=1`, `bkg_profile_frac=0` | Codex 10/26 탐색 설정의 재현(`docs/evidence/fwd23_effnp_variants.txt`). 이전 prefit은 짧은 배경 수명이 하한에 붙어 있었다. 새 prefit χ²/ndf 0.98(core 0.0728 mm 해). 2D knob도 바뀜(`study/rules/C23_effnp`): NP 1성분 유효 수명을 MC 평균에 제약, PR 폭 자유. final B +0.88%(간격비 1.8×MC) → +2.96%(MC 중심 제약만). Codex의 mass 배경 2차(+2.35%, mass χ²/ndf 0.87 → 1.35)는 근거가 B뿐이라 쓰지 않았다 |
| 2026-10-27 09:30 | fwd 10–12, 12–14 | 배경 ctau prefit: 10–12 성분 2/1/0 → 1/1/0, 둘 다 `bkg_wide_tau=1`, `bkg_core_max=3` | 이전 prefit은 10–12의 가장 짧은 수명이 코드 고정 구간 [0.062, 0.089] mm 상한에, 12–14의 core 폭이 상한(2 × MC)에 붙어 있었다. 바꾼 prefit은 붙는 매개변수가 없고 final B는 0.4% 이내로 같다(`docs/FWD1014_MC_CENTRED.md` 2.1절) |
| 2026-10-27 09:30 | fwd 10–12, 12–14, 14–20 | NP 성분 1 → 2 (10/25 변경을 되돌림) | 9월 fwd 10–20 정책(NP2 + 긴 수명 자유, 나머지 공분산 ×2)으로 복귀. 1성분 + MC 중심 제약으로는 fwd 12–14 −4.3%. 2번째 성분은 여전히 퇴화(12–14: 560 ± 1200)하며 그 임의성은 `docs/FWD1014_MC_CENTRED.md` 5.2절에 기록. 2성분 prefit χ²/ndf 0.84, 0.98, 0.86 |

`ctau_bkg.C`에 빈별 옵션 두 개를 추가했다(기본값은 이전과 같음): `bkg_wide_tau=1`은 가장 짧은 배경 수명 범위를 fwd 9–10·12–20처럼 넓게 연다, `bkg_core_max=k`는 배경 해상도 core 폭 상한을 MC PR core의 k배로 둔다(기본 2). fwd 5–12에서 시험했고 B를 0.1% 이상 바꾸지 않아 10/25 설정에는 쓰지 않았다; 10/27부터 fwd 10–12, 12–14에 쓴다. 10/27에 같은 꼴의 초기값·구간 훅 `bkg_par`(ctau_bkg.C), `pr_par`(ctau_pr.C), `np_par`(ctau_np.C)와 `bkg_plot_bins`, 질량 창 `mass_lo`/`mass_hi`를 추가했다(최종 설정에는 쓰지 않음).
