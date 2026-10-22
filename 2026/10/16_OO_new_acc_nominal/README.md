# OO new-acceptance nominal agreement

pp_HX의 새 single-muon acceptance를 적용한 OO 분석이다. 현재 선택은
`NO_COMMON_SCALE_SIGMA0414`이며 **20/20 bin이 09 nominal 통계오차의 0.414배 이내**다.
최대 차이는 mid 9–10의 **0.3710676064**다. fit 후 B 환산은 적용하지 않는다.

2026-10-23에 `~/work/26_pol_jpsi/pp_HX/`의 단계별 구조를 참고하여 운영 코드,
최종 그림, 실행 산출물, 과거 테스트를 분리했다. 이 정리는 fit 정책을 새로 선택하거나
물리 모형을 바꾸는 작업이 아니다. 선택된 ROOT와 prefit의 내용은 그대로 보존했다.

## 먼저 볼 파일

- [agreement PDF](results/nominal_agreement.pdf) · [PNG](results/nominal_agreement.png)
- [bin별 B·오차·rule](results/agreement.csv)
- [bin별 mass/ctau fit 그림](figs/) — 20개 bin, PDF 40개
- [실제 bin별 실행 설정](config/selection.json) · [rule 파일](config/rules/)
- [파일 보존 검증](results/verification.json)

## 구조와 Git 범위

```text
01_skim/                    데이터·PR·NP 입력 변환과 고정된 acceptance
02_fit/
  calibration/              mass, PR, NP, background prefit 매크로
  variants/                 현재 선택에 실제 사용한 final-fit 매크로 6종
  validate_saved_fit.C      저장된 fit 품질 검사
03_results/                 저장된 ROOT에서 B 추출 및 agreement 그림
config/
  selection.json            20개 bin의 macro, rule, window, prefit, ROOT 연결
  rules/<bin>/              선택된 bin별 설정
  input_sources.json        원본 skim/OniaTree의 외부 경로
  source_manifest.json      옮긴 소스의 해시와 원본 소스 해시
include/                    공용 헤더와 ROOT style
references/                 고정 nominal, 선택 결과표, FT1L 출처
run.py                      실행 진입점
tools/                      skim, fit, results, verify 실행 도구
results/                    사람이 볼 최종 표·agreement PDF/PNG·검증 요약 (추적)
figs/<bin>/                 사람이 볼 mass/ctau PDF (추적)
output/                     실행 산출물 전체 (Git 제외)
  inputs/                   acceptance가 적용된 데이터·MC ROOT
  prefits/<bin>/            실제 소비한 5개 prefit ROOT
  fit/<bin>/                선택된 final ROOT와 원본 fit 그림
  runs/                     재실행 명령·로그·ACLiC build
  frozen/                   수정하지 않은 최종 결과 bundle
  restructure_checks/       구조 변경 컴파일·회귀 검사 산출물
tests/archive/              기존 연구 코드·scan·테스트 결과·README (Git 제외)
```

`02_fit/variants/`의 여섯 파일은 현재 선택에서 사용하는 서로 다른 수치 구현이다.
사용하지 않은 실험 variant와 수백 개 후보 결과는 `tests/archive/`로 이동했다.
기존 archive 안의 절대경로는 당시 출처 기록이며 현재 실행 진입점이 아니다.
운영 코드는 archive 없이 실행된다. 이동 목록은 로컬 `output/restructure/relocation.json`에 있다.

이 폴더는 정리 전 Git에 추적되지 않았다. 이번 커밋부터 운영 소스·설정·참고 표와
최종 그림만 추적한다. 원시 입력, ROOT, 로그, 컴파일 결과, `tests/`, `output/`은
로컬에 보존하며 업로드하지 않는다. Git clone에는 이 로컬 데이터가 포함되지 않는다.
다른 서버에서 재현할 때는 `output/inputs/`, `output/prefits/`, `output/fit/`을 함께
복사한 뒤 `verify`로 해시를 확인한다. 원본 입력 생성 경로는 `config/input_sources.json`에서 바꾼다.

## 실행

Python 3.10 이상, ROOT 6.36.06, `pdftoppm`이 필요하다. ROOT는 `ROOT_EXECUTABLE`,
`ROOT_BIN`, PATH 순으로 찾고, 이 서버에서는 `~/build/root-6.36.06/bin/root`도 찾는다.

```bash
./run.py bins                     # bin, 실제 rule, ctau window
./run.py status                   # final ROOT와 prefit 유무
./run.py verify                   # 20 final + 100 prefit + 입력·소스·rule 해시
./run.py verify --code-only       # 로컬 ROOT 없이 추적된 소스·설정 검사
./run.py results                  # 실제 저장 ROOT에서 표·그림 재생성
./run.py results --help

# 기존 결과를 건드리지 않는 실행 계획
./run.py run-bin mid_20_30 --dry-run

# 선택된 prefit을 사용해 final fit 하나 재실행. 새 출력 위치를 사용한다.
./run.py run-bin mid_20_30 --output-base output/rerun
```

`run-bin`은 **선택된 다섯 prefit을 그대로 사용하는 final fit 재실행**이다.
prefit 자체를 다시 생성하거나 후보를 다시 선택하지 않는다. calibration 소스는
`02_fit/calibration/`에 있고, 현재 final runner의 입력은 `output/prefits/<bin>/`이다.
기본 출력에 기존 ROOT가 있으면 덮어쓰지 않고 멈춘다.

재실행 결과로 별도의 agreement를 만들려면 `config/selection.json`을 복사한 설정에서
해당 bin의 `output_root`와 `figure_dir`을 새 출력 경로로 지정하고 다음처럼 실행한다.
원래 선택 config와 reference hash는 보존한다.

```bash
./run.py results --config output/rerun_selection.json --allow-new-fits \
  --output-dir output/rerun_results --figures-dir output/rerun_figs
```

`--allow-new-fits`는 변경된 ROOT 해시를 허용하는 명시적인 옵션이다. ROOT의 fit 품질과
공통 resolution scale 부재 검사는 계속 적용한다. 새로운 값은 agreement 기준을
통과하지 않을 수도 있으며, 표와 요약에 실제 통과 수를 기록한다.

입력 변환은 다음과 같다. 기본 `output/inputs/`의 기존 파일은 덮어쓰지 않는다.
부분 실행은 운영 입력 바깥의 별도 출력 경로가 필수다.

```bash
./run.py skim data --dry-run
./run.py skim prompt --dry-run
./run.py skim nonprompt --dry-run
./run.py skim data --max-events 10 --output-dir output/skim_check
```

## FT1L은 무엇인가

**FT1L은 bin별로 선택한 기존 튜닝 설정과 background prefit을 모아 둔 assembly 이름**이다.
단일 수식이나 모든 bin에 동일한 제약을 적용하는 rule을 뜻하지 않는다.
약자의 풀이는 원본 코드에 기록돼 있지 않다.

원래 묶음은 forward 12개와 mid 8개였다. [원본 source-rule 목록](references/ft1l_assembly.json)을
보존했다. 예를 들면 다음과 같다.

| 원래 bin | FT1L 안의 source rule |
|---|---|
| fwd 1–2 | `LG13_m1_f4_n175_b4` |
| fwd 3–4 / 4–5 | `LG35_v130` / `LG35_v1000` |
| fwd 5–6 | `LG510_c0p01_mf_winB` |
| fwd 10–12 | `LG1020_w0p5_mm16_n14` |
| mid 8개 | `CT3` |

역사적 `rules.py`의 `_assembled_from()`은 각 source rule의 knobs와 background prefit을
bin별로 연결했다. 당시 `ft_pick.py`에는 status 0/covQual 3인 후보 중 nominal과의
차이가 가장 작은 후보를 찾는 로직이 있다. FT1L assembly JSON을 만든 정확한 명령까지
확인된 것은 아니다. 원본 코드는 로컬 archive와 이전 10월 15일 분석에 남아 있다.

**현재는 FT1L forward 11개만 유지**한다. fwd 2–3과 mid 전체는 다른 설정을 쓴다.
`config/selection.json`이 각 bin에서 실제 실행할 macro와 rule을 지정한다.
`rule_mode=native_environment`인 다섯 bin의 rule 파일은 기본 정책을 설명하는 기록이고,
실제 설정은 해당 macro와 명시된 환경변수에서 온다. 모든 `.knobs` 항목이 모든 모드에서
활성인 것은 아니다.

## 현재 선택의 제한

수치 agreement와 저장된 fit 품질 검사를 통과한 선택이다. bin별 예외를 정당화하거나
window/model 안정성 및 coverage를 확립한 결과로 해석하지 않는다.

- FT1L fwd 5–6의 첫 두 resolution 폭 좌표 prior σ는 `4e-6`, `8e-6` mm로,
  09 설정보다 70배 좁다. 이 강한 제약은 이번 구조 정리에서 바꾸지 않았다.
- mid 7–8과 12–14는 MC 중심 mean prior σ를 각각 MC 평균 오차의 64배와 32배로
  사용한다. 이는 제약 완화이며 nominal을 보고 강도를 선택한 기록이 있다.
- mid 9–10의 window 상한은 1.35 mm다. NP MC 유지율 97.3707%와 전체 fitted signal
  유지율 99.3255%는 서로 다른 지표다.
- fwd 2–3에는 B=0.11777인 선택 해 외에 B=0.15635인 해가 ΔNLL=0.409에 남아 있다.
- fwd 4–5 background는 기존 품질 기준을 통과하지만 보조 robust χ²/ndf=2.606이다.

09 nominal 값과 오차는 고정한다. `0.414`는 nominal 통계오차 단위의 비교 기준이며,
두 측정의 독립적인 통계적 유의도를 뜻하지 않는다.

## 구조 변경 검증

- 선택된 final ROOT 20개, 소비한 prefit 100개, accepted 입력 3개의 SHA256 보존.
- 20개 bin의 실제 ROOT에서 얻은 B/오차와 기존 표가 일치하고, agreement PNG도 바이트 단위로 일치.
- 옮긴 final-fit macro 6종을 각각 독립 ACLiC 디렉터리에서 컴파일.
- data/PR/NP 입력 변환을 각각 10 이벤트로 실행하여 경로와 출력을 확인.
- mid 20–30 final fit을 별도 출력에서 실제 재실행. 품질 검사 통과,
  floating parameter 12개의 값과 오차가 모두 기존 결과와 정확히 일치.

이 기록은 구조 변경의 회귀 검증이며 새 튜닝이나 전체 분석 재실행을 의미하지 않는다.
