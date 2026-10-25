# 저장된 최종 fit(FINAL_S1) 그리기

9월 튜닝 폴더(`../../09/19_OO_2dfit_tuning_Hijing_tuning/drawing`)의 그림 코드를 가져와 이 폴더 결과로 그린다. 재적합 없이 `../roots/`의 저장 결과만 읽는다.

```bash
bash drawing/refresh_all.sh          # 전체 갱신 (MAX_JOBS=6 으로 동시 작업 수 조정)
```

| 폴더 | 내용 |
|---|---|
| `final_summary/` | 20빈 결과표, PDF별 고정·제약·자유 표, nominal 비교 그림 (`report_table.py`, `nominal_agreement.py`) |
| `b_frac/` | B fraction vs pT (`bFraction_OO_{fwd,mid}_pp.C`). 최종 `modelResult`의 B와 오차를 직접 읽고, 사용 값은 `*_inputs.csv`에 남는다 |
| `raw_yields/` | Nsig / pT 폭 (`raw_yield_{fwd,mid}.C`), nominal 대비 비율 그림 포함 |
| `b_frac_comparison/` | nominal 대비 비교(2.5% 기준 색), \|B−B_nominal\| 계통 띠 그림, `bfraction*_agreement.csv` |
| `model_count/` | prefit 단계별 모델 성분 수 |
| `fit_convergence/` | 전 단계 RooFitResult 상태표 `fit_result_summary.csv` (status 0, covQual 3, EDM<1e-3) |
| `ctau_zoom/` | fwd pT 1–2, 2–3 lifetime 확대 그림 (−2~3 mm) |
| `projection_diagnostics/` | 최종 2D fit의 mass·ctau 투영 진단 PDF/CSV, `README.md` 표, `binning_overview.pdf` |
| `all_plots_in_one/` | `../figs_publish/`의 단계별 그림을 pT 순서로 모은 패널, B fraction·yield 사본 |

`latest_fit_results.h`가 최종 fit의 status/covQual/MIGRAD/HESSE와 Minuit2 사용을 확인하고, 통과하지 못하면 과거 배열로 대체하지 않고 실패한다. Nominal은 `../utils/nominal_bfraction.h`를 쓴다.

## 9월 폴더와 다른 점

- 투영 진단은 이 폴더에서 `fit2d.C`에만 훅이 있어 2D fit(20빈 × mass/ctau)만 만든다. prefit 단계 진단은 없다. `run_projection_diagnostics.py`는 `run_publish.sh --2d`를 `OO_DIAGNOSTICS_DIR`과 함께 돌려 빈별 production 설정(`OO_RULE_DIR` 등)을 그대로 쓰고, 실행 전후 `roots/` 해시가 같은지 확인한다. 이를 위해 `study/production.py`가 `OO_DIAGNOSTICS_DIR`을 매크로로 넘기도록 했다.
- `ctau_zoom`은 `OO_RULE_DIR=study/rules/FINAL_S1`으로 재생한다.
- ctau 배경 모델 스캔(`run_ctau_bkg_model_scan.sh`)은 재적합 스크립트라 가져오지 않았다.
