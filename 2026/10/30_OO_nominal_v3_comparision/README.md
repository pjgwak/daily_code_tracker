# OO J/ψ 질량–수명 2D fit: nominal ver3 기준 production (FINAL_V4)

28일 폴더(`../28_OO_nominal_v2_production`, FINAL_V2)를 복사해 만든 폴더다. 비교 기준을 이 폴더의 `bfraction_ctau_fit.xlsx`(nominal ver3, 10/30)로 바꾸고(`utils/nominal_bfraction.h`), 2D fit 규칙을 다시 정해 `study/rules/FINAL_V4`로 두었다(오전 판 `FINAL_V3`도 남김). prefit 설정은 fwd 1–2(ctau 창 [−6, 10])와 fwd 2–3(PR MC fit 범위 분위수 0.001/0.999)만 바꿨다(`docs/CHANGELOG_PREFIT.md`). fwd 2–3의 NP는 MC truth prefit과 같은 지수 2개를 두고 평균 수명만 MC에 제약한다(`np_moment_order 1`; 이전의 유효 1지수 방식과 같은 값, REPORT 9.1절). 결과는 20빈 모두 ver3와 ±2% 안(최대 1.96%), mid는 7/8이 1% 안이고, 보고서는 [docs/REPORT.md](docs/REPORT.md)다.

사용자가 준 입력: `bfraction_ctau_fit.xlsx`, `2dfit_fractions_nominal_compare.png`, `2dfit_systematics_{prompt,nonprompt}_yield.png`(FINAL_V2 기준 비교 그림).

## 실행

```bash
bash run_ctau_bins.sh --max-jobs=20   # prefit 전부 + 최종 2D fit (규칙 study/rules/FINAL_V4)
bash run_ctau_bins.sh --2d            # 최종 2D fit만
bash run_publish.sh                   # figs_publish/
python3 drawing/report_table.py; python3 drawing/nominal_agreement.py; python3 study/treatment.py
python3 drawing/v3_systematics.py docs/figs/final_v4 drawing/final_summary/summary.txt:FINAL_V4 docs/evidence/summary_FINAL_V3.txt:FINAL_V3 docs/evidence/summary_FINAL_V2_ver3.txt:FINAL_V2
```

## 스터디 도구 (28일 폴더에 더한 것)

| 경로 | 내용 |
|---|---|
| `study/rules.py` | 규칙 어휘에 `mcerr k`, `tailerr k`(폭 = MC prefit 오차 × k) 추가, V3/V3E/V3B/V3N/J2/J3/V3K 격자. ctau 창은 `config/prefit_settings.txt`에서 읽음 |
| `study/pr_width_err.C` → `config/pr_width_err.txt` | MC PR prefit 매개변수 오차 |
| `study/v3_summary.py` | 스캔 로그(`logs/scan_*.log`) 그룹별·빈별 정리 |
| `study/compose_final.py` + `config/final_v3_choice.txt`, `config/final_v4_choice.txt` | 빈별 규칙 선택을 모아 FINAL_V3/V4 knob 폴더 생성 (V4의 mid 4빈과 fwd 2–3은 선택표에 출처 명시) |
| `study/fwd23_window_study.py`, `study/constraint_width_check.py`, `docs/make_figs.sh` | ctau 창 변형 체인, joint 세기의 제약 폭 점검, 웹페이지 그림 PNG 변환 |
| `study/choice_summary.py` | 선택 파일 + 스캔 로그로 결과표(재적합 없음) |
| `study/bkg_mass_study.py` | 질량 창·Chebychev 차수 prefit 변형 체인 일괄 실행 |
| `drawing/v3_systematics.py` | B 비교 + NP/PR 수율 계통 오차 그림 |
| `docs/evidence/` | 스캔 전체(`all_candidate_fits.txt`), 그룹 순위(`rule_scan_*.txt`), 배경 질량 변형, prefit 감사, 결과표 |
| `roots_FINAL_V2_backup/`, `config/prefit_settings_FINAL_V2.txt` | 28일 prefit과 설정 |
