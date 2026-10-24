# OO J/ψ 질량–수명 2D fit: PDF별 규칙 (new acceptance)

최종 설정 `FINAL_v3`로 20개 빈 모두 nominal B fraction과 ±2% 또는 0.138σ 이내에서 일치한다. 방법과 결과는 [docs/REPORT.md](docs/REPORT.md), mid 7–9·fwd 2–3 상세는 [docs/MID79_FWD23.md](docs/MID79_FWD23.md)에 있다.

![nominal 비교](drawing/final_summary/nominal_agreement.png)

## 폴더 구성 (9월 튜닝 폴더와 같은 구조)

| 경로 | 내용 |
|---|---|
| `mc_mass.C`, `ctau_pr.C`, `ctau_np.C`, `mass_window_screen.C`, `ctau_bkg.C`, `fit2d.C` | 단계별 적합 매크로 (`mass.C`는 범위 인자 없는 이전 DATA mass 매크로) |
| `roots/y*/{mc_mass,ctau_pr,ctau_np,mass,ctau_bkg,fit2d}/` | 단계별 fit 결과 ROOT (prefit과 최종 2D fit) |
| `figs/y*/<단계>/` | 단계별 fit 그림 |
| `figs_publish/y*/<단계>/` | 발표용 그림 (`run_publish.sh`) |
| `drawing/final_summary/` | 20빈 결과표(`summary.tsv`, `report_table.tsv`), PDF별 고정·제약·자유 표(`treatment.md`), nominal 비교 그림 |
| `config/prefit_settings.tsv` | 빈별 prefit 선택: ctau 범위, 배경 성분 수, mass 배경 차수, NP·PR 성분 수 |
| `config/pr_calibration_windows.tsv` | MC PR prefit의 범위와 파일 지문 (fit2d가 확인) |
| `study/rules/FINAL_v3/` | 최종 2D fit knob 파일 (빈별 고정·제약·자유 설정) |
| `study/` | `production.py` (실행 스크립트의 본체), `rules.py` (그룹 규칙 정의), `treatment.py` + `dump_treatment.C` (PDF별 고정·제약·자유 표), `rules/` (knob), `review_20261024/` (ctau 범위 근거: 23일 window 스터디) |
| `inputs/` | DATA, PR MC, NP MC 입력 (new acceptance) |
| `docs/` | 보고서, 웹 페이지 원본, 근거 요약표(`docs/evidence/`) |
| `logs/` | 실행 로그 |

## 실행

ROOT 6.36.06 (`~/build/root-6.36.06/bin/root`)을 쓴다.

```bash
bash run_ctau_bins.sh                 # 전 단계: mc → pr → np → mass → bkg → 2d
bash run_ctau_bins.sh --2d            # 최종 2D fit만 (roots/의 prefit 사용)
bash run_ctau_bins.sh --bkg --pt=8:9 --rapidity=mid
bash run_ctau_bins.sh --redraw        # 저장된 결과로 figs/ 다시 그리기 (재적합 없음)
bash run_publish.sh                   # figs_publish/ 만들기 (재적합 없음)
python3 drawing/nominal_agreement.py  # nominal 비교 그림
python3 study/treatment.py            # PDF별 고정·제약·자유 표
```

`--2d` 실행이 끝나면 `drawing/final_summary/summary.tsv`가 갱신된다. PR prefit을 다시 만들면 `config/pr_calibration_windows.tsv`의 파일 지문도 자동으로 갱신된다.
