# `fit2d.C` 모델 구성 및 파라미터 처리 요약

대상 코드: `fit2d.C`

## 1. 입력과 관측 변수

- 입력 데이터: `RooDataSet_miniAOD_isMC0_globalOff_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root`의 `dataset`
- 관측 변수:
  - `mass`: `2.6 < mass < 3.5`, 단위 `GeV/c^{2}`
  - `ctau3D`: bin별로 자동 설정되는 `ctRange`, 단위 `mm`
- 기본 선택:
  - `ptLow < pt < ptHigh`
  - `yLow < |y| < yHigh`
  - `2.6 < mass < 3.5`
- `ctau3D` 범위:
  - 기본: `yLow == 0`이면 `[-1, 3]`, 아니면 `[-6, 6]`
  - forward low-pt 특수 bin:
    - `1.6 <= yLow`, `1 < pT < 2` bin: `[-20, 20]`
    - `1.6 <= yLow`, `2 < pT < 3` bin: `[-8, 10]`
    - `1.6 <= yLow`, `ptHigh >= 12`: `[-2, 4]`
- 보조 데이터셋:
  - signal-region: `2.9 <= mass <= 3.2`
  - sideband: `2.6 <= mass < 2.9` 또는 `3.2 < mass <= 3.5`

## 2. 외부 1D fit 결과 의존성

`fit2d.C`는 2D fit 모델을 처음부터 독립적으로 정하지 않고, 아래 ROOT 파일의 저장된 설정과 fit 결과를 seed/template로 사용한다.

| 용도 | 파일 | 주요 사용값 |
|---|---|---|
| mass 모델 | `roots/<yTag>/mass/mass_model_<figTag>.root` | mass signal/background 구성 수, shape 파라미터, `Nsig`, `Nbkg`, 에러 |
| prompt ctau resolution | `roots/<yTag>/ctau_pr/ctau_resolution_<figTag>.root` | resolution component 수, Gaussian mean/scale, fraction |
| nonprompt ctau signal | `roots/<yTag>/ctau_np/ctau_np_model_<figTag>.root` | nonprompt single-sided component 수, lifetime, component yield |
| background ctau | `roots/<yTag>/ctau_bkg/ctau_bkg_fitresult_<figTag>.root` | background resolution, lifetime offsets, component yields, 에러 |

## 3. Mass 모델

### Signal mass PDF

저장된 `mass.C` 결과에서 다음 component 수를 읽는다.

- `nSignalGaussComponents`: `0..2`
- `nSignalCBComponents`: `0..2`

구성 가능한 signal mass component:

- `signal_mass_gaus`: `RooGaussian(mass, signal_mass_mean, signal_mass_sigma)`
- `signal_mass_gaus2`: `RooGaussian(mass, signal_mass_mean, signal_mass_sigma2)`
- `signal_mass_cb`:
  - `nSignalCBComponents == 1`: `RooCBShape`
  - `nSignalCBComponents >= 2`: `RooCrystalBall`

결합 방식:

- component가 1개면 해당 PDF를 그대로 사용
- 2개 이상이면 `RooAddPdf("signal_mass", ..., recursiveFraction=true)` 사용
- fraction은 ratio 파라미터로부터 계산:
  - `signal_mass_frac1 = ratio1 / (1 + ratio1)`
  - `signal_mass_frac2 = ratio2 / (1 + ratio2)`
  - `signal_mass_frac3 = ratio3 / (1 + ratio3)`

초기값:

- `signal_mass_mean`, `signal_mass_sigma`, CB tail 파라미터, fraction ratio는 `mass_model` ROOT 파일에서 읽는다.
- 없으면 코드 내 fallback 사용:
  - mean `3.096`
  - sigma `0.03`
  - CB sigma `0.035`
  - alpha `1.5`, alpha2 `2.0`
  - n `3.0`, n2 `4.0`
  - ratio1 `1.86`, ratio2 `0.25`, ratio3 `0.111111`

양수/순서 제약 구현:

- `signal_mass_sigma2 = sigmaFloor + exp(signal_mass_sigma2_log_offset)`, `sigmaFloor = 1e-4`
- `nSignalGaussComponents >= 1`이면 CB sigma도 log-offset 형태로 만든다.
- `signal_mass_sigma_delta2`, `signal_mass_cb_sigma_delta2` 같은 delta formula가 정의되지만, 최종 fit 제약으로 직접 쓰이지는 않는다.

고정 여부:

- 최종 2D fit 직전 mass signal shape는 고정된다.
- 고정되는 항목:
  - `signal_mass_mean`
  - `signal_mass_sigma`
  - `signal_mass_sigma2_log_offset`
  - `signal_mass_cb_sigma_log_offset` 또는 `signal_mass_cb_sigma_base`
  - `signal_mass_cb_sigma2_log_offset`
  - `signal_mass_cb_alpha`, `signal_mass_cb_n`
  - `signal_mass_cb_alpha2`, `signal_mass_cb_n2`는 double-sided CB일 때
  - `signal_mass_frac_ratio1/2/3`

### Background mass PDF

저장된 `mass.C` 결과에서 다음을 읽는다.

- `nBkgExpComponents`: `0..1`
- `nBkgChebyOrder`: `0..6`

구성:

- `nBkgExpComponents == 1`: `RooExponential("bkg_mass", mass, bkg_mass_lambda)`
- 그 외: `RooChebychev("bkg_mass", mass, bkg_mass_p1..pN)`

검사 조건:

- exponential과 Chebychev가 동시에 켜져 있으면 에러로 종료
- background component가 하나도 없으면 에러로 종료

초기값:

- `bkg_mass_lambda`, `bkg_mass_p1..p6`는 `mass_model` ROOT 파일에서 읽는다.
- fallback:
  - lambda `-1.0`, 범위 `[-10, -1e-4]`
  - Cheby coefficient `0.0`, 범위 `[-1, 1]`

고정 여부:

- 최종 2D fit에서는 background mass shape도 고정된다.
- exponential이면 `bkg_mass_lambda` 고정
- Chebychev이면 `bkg_mass_p1..p6` 모두 고정

## 4. Prompt ctau resolution 모델

저장된 `ctau_pr.C` 결과에서 `nResolutionComponents`를 `1..4`로 읽는다.

구성:

- 각 component는 `RooGaussModel(ctau3D, mean_i, scale_i)`
- 1개면 `ctauTime1` 단독 사용
- 2개 이상이면 `RooAddModel("promptResolutionModel", ...)`로 결합

mean 처리:

- `ctauTime1Mean..ctauTime4Mean`은 `RooConstVar`
- 즉 prompt resolution mean은 2D fit에서 항상 고정

scale 처리:

- `ctauTime1Scale`: `RooRealVar`
- `ctauTime2Scale = ctauTime1Scale + ctauTime2Delta`
- `ctauTime3Scale = ctauTime2Scale + ctauTime3Delta`
- `ctauTime4Scale = ctauTime3Scale + ctauTime4Delta`

따라서 scale 순서는 delta 양수 범위로 보장된다.

초기값과 범위:

- scale seed는 `ctau_resolution` ROOT 파일에서 읽는다.
- `fitRangeAround(center, hardLow, hardHigh, relWidth, absWidth)`로 seed 주변 제한 범위를 만든다.
- component별 hard range:
  - `ctauTime1Scale`: `0.005..0.20`
  - `ctauTime2Delta`: `0.001..0.30`
  - `ctauTime3Delta`: `0.001..0.60`
  - `ctauTime4Delta`: `0.05..1.20`

fraction 처리:

- 저장된 `ctauFrac1..3`을 component fraction으로 변환한 뒤 `RooConstVar`로 둔다.
- 즉 resolution fraction은 2D fit에서 고정이다.

고정/제약:

- `fixResolutionParams == true`이면 `ctauTime1Scale`, `ctauTime2Delta`, `ctauTime3Delta`, `ctauTime4Delta`를 고정한다.
- 기본값 `fixResolutionParams == false`에서는 scale/delta가 float한다.
- 단, 최종 fit 전에 Gaussian external constraint가 추가된다.
  - `ctauTime1Scale`: seed 중심, sigma `0.010`
  - `ctauTime2Delta`: seed 중심, sigma `0.020`
  - `ctauTime3Delta`: seed 중심, sigma `0.040`
  - `ctauTime4Delta`: seed 중심, sigma `0.080`

## 5. Signal ctau 모델

### Prompt signal

- prompt signal time PDF는 prompt resolution 자체를 사용한다.
- 즉 `promptTimePdf = promptResolutionPtr`

### Nonprompt signal

저장된 `ctau_np.C` 결과에서 `nSignalSSComponents`를 `1..3`으로 읽는다.

구성:

- 각 component는 prompt resolution과 convolution된 `RooDecay(..., RooDecay::SingleSided)`
  - `signal_ss1_time`
  - `signal_ss2_time`
  - `signal_ss3_time`
- component가 1개면 단독 사용
- 2개 이상이면 `RooAddPdf("signal_np_time", ...)`로 결합

결합 coefficient:

- `NsignalSS1`, `NsignalSS2`, `NsignalSS3`를 `ctau_np` ROOT 파일에서 읽어 `RooConstVar`로 사용
- 즉 nonprompt 내부 component 비율은 2D fit에서 고정이다.

lifetime 처리:

- `signal_lifetime = signalLifetimeFloorSaved + exp(signal_log_lifetime_offset)`
- `signal2_lifetime = signal_lifetime + signal2_lifetime_delta`
- `signal3_lifetime = signal2_lifetime + signal3_lifetime_delta`

초기값:

- lifetime과 log-offset은 `ctau_np` ROOT 파일에서 읽는다.
- 첫 번째 lifetime seed는 resolution scale 기반 floor보다 작지 않게 조정한다.
- lifetime 상한은 `max(signalLifetimeInit * 5, maxStableLifetime)` 기반이다.

고정 여부:

- `signal_log_lifetime_offset`은 항상 고정된다.
- 2번째 component가 있으면 `signal2_lifetime_delta`도 고정된다.
- 3번째 component가 있으면 `signal3_lifetime_delta`도 고정된다.
- 따라서 nonprompt signal lifetime shape는 `ctau_np.C` template로 고정하고, 2D fit에서 직접 shape를 재조정하지 않는다.

### Prompt + nonprompt signal 결합

`signal_time`은 다음처럼 구성된다.

```text
signal_time = RooAddPdf(signal_np_time, promptTimePdf, bFraction, recursiveFraction=true)
```

- 첫 번째 component가 `signal_np_time`이므로 `bFraction`은 nonprompt signal fraction이다.
- prompt fraction은 `1 - bFraction`이다.

`bFraction` 초기값과 범위:

- `1 < pT < 2`, `1.6 < |y| < 2.4`: 초기값 `0.08`
- `2 < pT < 3`, `1.6 < |y| < 2.4`: 초기값 `0.11`
- 그 외: 초기값 `0.10`
- 범위는 항상 `[0.05, 0.12]`

고정/제약:

- `bFraction`은 float한다.
- Gaussian external constraint는 걸지 않는다.
- 범위 제한만 적용된다.

## 6. Background ctau 모델

background ctau 모델은 `ctau_bkg.C`의 sideband fit 결과를 seed로 사용한다.

### Background prompt-like resolution

구성:

- `bkgCtauTime1..4`: `RooGaussModel`
- `nBkgResolutionComponents`가 2 이상이면 `RooAddModel("bkgPromptResolutionModel", ...)`

mean/fraction:

- mean은 `RooConstVar`
- fraction `bkgCtauFrac1..3`도 `RooConstVar`

scale:

- signal prompt resolution과 유사하게 누적 delta 방식:
  - `bkgCtauTime2Scale = bkgCtauTime1Scale + bkgCtauTime2Delta`
  - `bkgCtauTime3Scale = bkgCtauTime2Scale + bkgCtauTime3Delta`
  - `bkgCtauTime4Scale = bkgCtauTime3Scale + bkgCtauTime4Delta`

고정 여부:

- `bkgCtauTime1Scale`, `bkgCtauTime2Delta`, `bkgCtauTime3Delta`, `bkgCtauTime4Delta`는 항상 고정된다.
- 즉 background prompt-like resolution shape는 2D fit에서 float하지 않는다.

### Background lifetime components

저장된 `ctau_bkg.C` 결과에서 component 수를 읽는다.

- `nBkgSSComponents`: single-sided, `0..3`
- `nBkgFlipComponents`: flipped, `0..3`
- `nBkgDSComponents`: double-sided, `0..3`

항상 포함되는 component:

- prompt-like background: `bkgPromptTimePdf`, coefficient `Nplb`

선택적으로 포함되는 component:

- single-sided:
  - `bkg_time_ss`, `bkg_time_ss2`, `bkg_time_ss3`
- double-sided:
  - `bkg_time_ds`, `bkg_time_ds2`, `bkg_time_ds3`
- flipped:
  - `bkg_time_flip`, `bkg_time_flip2`, `bkg_time_flip3`

lifetime 구현:

- 모든 lifetime은 floor plus exponential offset 형태다.
- 예:
  - `bkg_lifetime = bkgLifetimeFloor + exp(bkg_log_lifetime_offset)`
  - `bkg_sym_lifetime = bkgSymLifetimeFloor + exp(bkg_sym_log_lifetime_offset)`
  - `bkg_flip_lifetime = bkgFlipLifetimeFloor + exp(bkg_flip_log_lifetime_offset)`

초기값:

- 각 log-offset은 `bkgTimeResult`에서 읽는다.
- fallback은 floor보다 충분히 큰 lifetime을 만드는 log 값으로 설정한다.

범위 설정:

- `setLogLifetimeRange(...)` 또는 `setLogLifetimeRangeAround(...)`로 log-offset 범위를 설정한다.
- symmetric double-sided background lifetime은 `bkgTimeResult`의 error를 이용해 중심 주변 범위를 둔다.
- 그 외는 hard low/high lifetime으로부터 log-offset 범위를 만든다.

고정/제약:

- background lifetime log-offset들은 `setConstant(false)`로 float한다.
- 이후 Gaussian external constraint를 추가한다.
- 중심값은 `bkgTimeResult`의 최종값이다.
- sigma는 `bkgTimeResult`의 error와 fallback sigma 중 큰 값이다.
  - `bkg_log_lifetime_offset`: fallback `0.12`
  - `bkg_log_lifetime2_offset`: fallback `0.16`
  - `bkg_log_lifetime3_offset`: fallback `0.20`
  - `bkg_sym_log_lifetime_offset`: fallback `0.15`
  - `bkg_sym_log_lifetime2_offset`: fallback `0.20`
  - `bkg_sym_log_lifetime3_offset`: fallback `0.25`
  - `bkg_flip_log_lifetime_offset`: fallback `0.12`
  - `bkg_flip_log_lifetime2_offset`: fallback `0.16`
  - `bkg_flip_log_lifetime3_offset`: fallback `0.20`

### Background ctau component 결합

`addBkgComponent(...)`로 component를 하나씩 `bkgTimePdfList`에 넣고, 각 coefficient를 `RooRealVar`로 만든다.

coefficient 초기값:

- `Nplb`, `Nss1..3`, `Nds1..3`, `Nflip1..3`는 `bkgTimeResult`에서 읽는다.
- 최소 초기값은 `1e-6`

coefficient 범위:

- `[1e-6, max(10 * bkgYieldSum, 10 * init)]`

coefficient 제약:

- 각 coefficient마다 Gaussian constraint를 추가한다.
- 중심값은 sideband fit에서 읽은 yield-like 값
- sigma는 `max(0.25 * init, 0.05 * bkgYieldSum, 1.0)`

결합:

- component가 1개면 해당 PDF 단독 사용
- 2개 이상이면 `RooAddPdf("bkg_time", bkgTimePdfList, bkgTimeCoeffList)` 사용

## 7. 최종 2D 모델 결합

중간 product PDF:

```text
prompt_core = signal_mass_pdf * promptTimePdf
np_core     = signal_mass_pdf * signal_np_time
signal_core = signal_mass_pdf * signal_time
bkg_core    = bkg_mass_pdf    * bkg_time_pdf
```

실제 fit 모델:

```text
model = RooAddPdf(signal_core, bkg_core, Nsig, Nbkg)
```

- `Nsig`, `Nbkg`를 coefficient로 쓰는 extended 2D model이다.
- fit 호출도 `Extended(true)`를 사용한다.
- `prompt_core`, `np_core`는 최종 model의 직접 component가 아니라 lifetime plot에서 prompt/nonprompt 분해를 그릴 때 사용된다.

## 8. Yield 파라미터

### `Nsig`

- 초기값: `mass_model` ROOT 파일의 `Nsig`
- fallback: `max(1, 0.5 * data->numEntries())`
- 범위: `[0, 2 * data->numEntries()]`
- 고정 여부: float
- 제약: Gaussian external constraint
  - 중심: mass fit의 `Nsig`
  - sigma: mass fit error 또는 fallback `max(1, 0.05 * data->numEntries())`

### `Nbkg`

- 초기값: `mass_model` ROOT 파일의 `Nbkg`
- fallback: `max(1, 0.5 * data->numEntries())`
- 범위: `[0, 2 * data->numEntries()]`
- 고정 여부: float
- 제약: Gaussian external constraint
  - 중심: mass fit의 `Nbkg`
  - sigma: mass fit error 또는 fallback `max(1, 0.05 * data->numEntries())`

## 9. External constraints 전체 목록

최종 fit의 `ExternalConstraints(massConstraints)`에 들어가는 항목은 다음이다.

| 그룹 | 파라미터 | 의미 |
|---|---|---|
| mass yield | `Nsig`, `Nbkg` | mass 1D fit yield를 중심으로 제약 |
| prompt resolution | `ctauTime1Scale`, `ctauTime2Delta`, `ctauTime3Delta`, `ctauTime4Delta` | prompt resolution scale/delta를 `ctau_pr` seed 중심으로 제약 |
| background lifetime | `bkg_*_log_lifetime*_offset` | `ctau_bkg` sideband lifetime shape를 중심으로 제약 |
| background component coefficient | `Nplb`, `Nss*`, `Nds*`, `Nflip*` | `ctau_bkg` sideband component weight를 중심으로 제약 |

주의:

- mass shape 파라미터들은 Gaussian constraint를 거는 대신 2D fit에서 고정한다.
- background mass shape도 고정한다.
- signal nonprompt lifetime shape도 고정한다.
- `bFraction`은 external constraint 없이 float한다.

## 10. Fit 설정

활성 fit 경로:

```text
model.fitTo(
  data,
  Save(true),
  Extended(true),
  Offset(true),
  RecoverFromUndefinedRegions(1.0),
  Strategy(2),
  Minimizer("Minuit2"),
  Optimize(false),
  PrintLevel(1),
  ExternalConstraints(massConstraints)
)
```

`massConstraints`가 비어 있으면 `ExternalConstraints` 없이 같은 옵션으로 fit한다. 현재 코드 흐름에서는 일반적으로 yield, resolution, background lifetime, background coefficient constraint가 추가된다.

## 11. Float/fixed 요약

| 파라미터 그룹 | 2D fit 처리 |
|---|---|
| signal mass shape | fixed |
| background mass shape | fixed |
| `Nsig`, `Nbkg` | floating + Gaussian constrained |
| prompt resolution mean | fixed |
| prompt resolution fraction | fixed |
| prompt resolution scale/delta | floating + Gaussian constrained, 단 `fixResolutionParams=true`이면 fixed |
| signal nonprompt lifetime | fixed |
| signal nonprompt internal component weights | fixed |
| `bFraction` | floating, range constrained only |
| background prompt-like resolution | fixed |
| background lifetime log-offsets | floating + Gaussian constrained |
| background ctau component coefficients | floating + Gaussian constrained |

