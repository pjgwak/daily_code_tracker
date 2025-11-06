# fit2D_step_by_step.py
# PyROOT port of fit2D_step_by_step.C (core flow + models + scaleF + fits)
# Author: ported for 관리자님
import math
import sys
from array import array
from ROOT import (
    RooFit, RooRealVar, RooWorkspace, RooArgSet, RooArgList, RooDataSet, RooDataHist,
    RooFormulaVar, RooAddPdf, RooProdPdf, RooGaussian, RooCBShape, RooChebychev,
    RooExponential, RooGaussModel, RooDecay, RooBinning, RooPlot, TFile, TCanvas, gSystem
)

# -----------------------------
# 기본 설정 (원본과 동일 파라미터)
# -----------------------------
ptLow, ptHigh = 6.5, 7.5
yLow, yHigh   = 0.0, 2.4
massLow, massHigh = 2.6, 3.5
ctLow, ctHigh = -1.0, 4.0
ctErrLow, ctErrHigh = 0.008, 0.3
# 사이드밴드 / 시그널 대역
sbL_lo, sbL_hi = 2.6, 2.9
sig_lo, sig_hi = 2.9, 3.3
sbR_lo, sbR_hi = 3.3, 3.5
# Chebychev 차수 (원본 변수 bkgMassOrder)
bkgMassOrder = 2  # 0~6 지원으로 제한 확인 로직 포함됨. :contentReference[oaicite:0]{index=0}

# --------------------------------
# Workspace & Observables
# --------------------------------
ws = RooWorkspace("ws", True)

mass   = RooRealVar("mass", "m_{#mu#mu}", massLow, massHigh, "GeV")
ctau3D = RooRealVar("ctau3D", "#font[12]{l}_{J/#psi}", ctLow, ctHigh, "mm")
ctErr  = RooRealVar("ctau3DErr", "#sigma(#font[12]{l})", ctErrLow, ctErrHigh, "mm")
ws.importClassCode()  # 편의
getvar = ws.var

getattr(ws, 'import')(mass)
getattr(ws, 'import')(ctau3D)
getattr(ws, 'import')(ctErr)

# --------------------------------------------------
# Mass signal: Gaussian + CrystalBall (SUM::G1CB1Sig) :contentReference[oaicite:1]{index=1}
# --------------------------------------------------
def define_mass_sig(ws):
    ws.factory("Gaussian::G1Sig(mass,meanSig[3.0975,3.05,3.15],sigmaSig1[0.03,0.008,0.1])")
    ws.factory("CBShape::CB1Sig(mass,meanSig,sigmaSig2[0.03,0.0008,0.075],alpha[1.9,1.2,2.8],enne[2.5,1.0,4.0])")
    ws.factory("SUM::G1CB1Sig(fracG1[0.5,0.01,0.99]*G1Sig,CB1Sig)")

# ---------------------------------------------
# Mass background: Polynomial + Exponential 등 :contentReference[oaicite:2]{index=2}
# ---------------------------------------------
def define_mass_bkg(ws):
    ws.factory("Exponential::expBkg(mass,coefExp[-1.,-3.,1.])")
    # 필요시 Chebychev로 교체 가능: RooChebychev 대신 factory 문자열로도 가능

# ---------------------------------------------------------
# Prompt-ctau Resolution (3-Gauss AddModel) with PEE 구조 :contentReference[oaicite:3]{index=3}
# ---------------------------------------------------------
def define_ct_pr_res(ws):
    # N(협의: 핵심은 sigmaPRResN, rW, rW2로 계층적 폭 구성) :contentReference[oaicite:4]{index=4}
    ws.factory("sigmaPRResN[0.05,0.005,8]")
    ws.factory("rW[1.5,1.0001,1000]")
    ws.factory("rW2[1.5,1.0001,1000]")
    ws.factory("prod::sigmaPRResW(sigmaPRResN, rW)")
    ws.factory("prod::sigmaPRResW2(sigmaPRResW, rW2)")
    ws.factory("GaussModel::GN_PRRes(ctau3D, meanPRResN[0], sigmaPRResN, one[1.0], ctau3DErr)")
    ws.factory("GaussModel::GW_PRRes(ctau3D, meanPRResW[0], sigmaPRResW, one, ctau3DErr)")
    ws.factory("GaussModel::GW_PRRes2(ctau3D, meanPRResW[0], sigmaPRResW2, one, ctau3DErr)")
    ws.factory("AddModel::CtPRRes({GN_PRRes, GW_PRRes, GW_PRRes2},{fracRes1[0.1,0.01,0.99999], fracRes2[0.2,0.01,0.99999]})")

# ------------------------------------------------------
# Bkg ctau 모델 (좌/우/양쪽 + Res), 합성(예: 2L,2M,2R) :contentReference[oaicite:5]{index=5}
# ------------------------------------------------------
def define_ct_bkg(ws):
    # decay(+convolution with CtPRRes)
    ws.factory("Decay::CtBkgPos(ctau3D,lambdap[0.02,0.001,1],CtPRRes,RooDecay::SingleSided)")
    ws.factory("Decay::CtBkgPos2(ctau3D,lambdap2[0.02,0.001,1],CtPRRes,RooDecay::SingleSided)")
    ws.factory("Decay::CtBkgNeg(ctau3D,lambdam[0.02,0.001,2],CtPRRes,RooDecay::Flipped)")
    ws.factory("Decay::CtBkgNeg2(ctau3D,lambdam2[0.02,0.001,2],CtPRRes,RooDecay::DoubleSided)")
    ws.factory("Decay::CtBkgDbl(ctau3D,lambdasym[0.01,0.001,2],CtPRRes,RooDecay::DoubleSided)")
    ws.factory("Decay::CtBkgDbl2(ctau3D,lambdasym2[0.01,0.001,4],CtPRRes,RooDecay::DoubleSided)")

    # SUM of components + residual resolution part
    ws.factory(
        "SUM::CtBkgTot(fracCtBkg1[0.1,0.01,0.99]*CtBkgPos,"
        "               fracCtBkg2[0.1,0.01,0.99]*CtBkgPos2,"
        "               fracCtBkg3[0.1,0.01,0.99]*CtBkgNeg,"
        "               fracCtBkg4[0.1,0.01,0.99]*CtBkgNeg2,"
        "               fracCtBkg5[0.1,0.01,0.99]*CtBkgDbl,"
        "               fracCtBkg6[0.1,0.01,0.99]*CtBkgDbl2,"
        "               CtPRRes)"
    )

# ----------------------------------------------------------
# Mass PDF, 2D 결합, PEE 포함 구성의 핵심 조합 (요약형)
# (원본은 Mass/ctau × err 분포를 RooProdPdf로 결합) :contentReference[oaicite:6]{index=6}
# ----------------------------------------------------------
def build_pdfs(ws):
    # signal mass
    define_mass_sig(ws)
    # background mass
    define_mass_bkg(ws)
    # mass total (간단히: 신호+지수Bkg)
    ws.factory("NSig[10000,0,1e8]")
    ws.factory("NBkg[10000,0,1e8]")
    ws.factory("SUM::MassPDF(NSig*G1CB1Sig, NBkg*expBkg)")

    # prompt resolution
    define_ct_pr_res(ws)
    # ct background
    define_ct_bkg(ws)

    # ct signal(요약): prompt/NP/Res 등은 분석자 설정에 맞춰 추가 가능
    # 여기서는 시연을 위해 prompt-Res만 구성해 Mass와 결합 예시
    ws.factory("Bfrac[0.25,0.0,1.0]")  # placeholder

    # Mass × Ct 결합 PDF 예시 (배경항)
    # errPdf를 따로 RooHistPdf로 만드는 흐름은 데이터 로딩 후 적용.
    # 여기선 구조만 세움.
    # MassCtBkg = expBkg × CtBkgTot
    ws.factory("PROD::MassCtBkg(expBkg, CtBkgTot)")

    # MassCtSig 간단 예: G1CB1Sig × CtPRRes
    ws.factory("PROD::MassCtPR(G1CB1Sig, CtPRRes)")

    # 총합 모형 (RSUM/SUM 선택 가능). 이 예시는 배경 비율을 Mass에서 NBkg로 제어.
    ws.factory("RSUM::totPDF(fracBkg[0.5,0,1]*MassCtBkg, MassCtPR)")

# ----------------------------------------------------------
# Chebyshev 기반 SB→SIG 스케일 팩터 계산 (원본 동일 로직) :contentReference[oaicite:7]{index=7}
# ----------------------------------------------------------
def compute_scaleF(order, massMin, massMax, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi, allPars, coeffPrefix="sl"):
    if order < 0 or order > 6:
        raise ValueError(f"Unsupported Cheby order {order} (allowed 0..6)")

    def to_unit_x(m):
        return 2.0 * (m - massMin) / (massMax - massMin) - 1.0

    def chebT(n, x):
        if n == 0: return 1.0
        if n == 1: return x
        Tn_2, Tn_1 = 1.0, x
        for k in range(2, n+1):
            Tn = 2.0 * x * Tn_1 - Tn_2
            Tn_2, Tn_1 = Tn_1, Tn
        return Tn_1

    def chebVal(m, coeffs):
        x = to_unit_x(m)
        v = 1.0
        for i, c in enumerate(coeffs, start=1):
            v += c * chebT(i, x)
        return v

    # allPars: RooArgList에서 이름으로 찾아 값 얻기
    def getVal(name):
        obj = allPars.find(name)
        if obj:
            rr = obj
            # RooAbsReal 호환
            return rr.getVal()
        return float("nan")

    coeffs = []
    for i in range(1, order+1):
        val = getVal(f"{coeffPrefix}{i}")
        if not math.isfinite(val):
            raise RuntimeError(f"Coefficient '{coeffPrefix}{i}' not found or NaN.")
        coeffs.append(val)

    def integrate_simpson(f, a, b, nseg=2048):
        if b <= a: return 0.0
        if nseg % 2: nseg += 1
        h = (b - a) / nseg
        s = f(a) + f(b)
        for i in range(1, nseg):
            x = a + i*h
            s += f(x) * (4.0 if (i % 2) else 2.0)
        return s * (h/3.0)

    f = lambda m: chebVal(m, coeffs)
    I_sig = integrate_simpson(f, sig_lo, sig_hi)
    I_sb  = integrate_simpson(f, sbL_lo, sbL_hi) + integrate_simpson(f, sbR_lo, sbR_hi)
    if abs(I_sb) < 1e-300:
        raise RuntimeError("Sideband integral is zero. Check range.")
    return I_sig / I_sb

# ----------------------------------------------------------
# 예시 데이터 준비 (사용자 데이터에 맞게 교체)
# 실제로는 TTree로부터 RooDataSet/RooDataHist를 빌드해야 함.
# ----------------------------------------------------------
def make_dummy_data(ws):
    # 실제 분석에서는 mass, ctau3D, ctau3DErr를 채운 RooDataSet을 로드하세요.
    # 여기서는 예시로 빈 데이터셋을 만듭니다.
    return RooDataSet("data", "data", RooArgSet(ws.var("mass"), ws.var("ctau3D"), ws.var("ctau3DErr")))

# ----------------------------------------------------------
# Main flow = fit2D_step_by_step() 포트 (요약형) :contentReference[oaicite:8]{index=8}:contentReference[oaicite:9]{index=9}
# ----------------------------------------------------------
def fit2D_step_by_step_py():
    print("=== start fit2D_step_by_step_py() ===")
    # 1) 모델 구성
    build_pdfs(ws)

    # 2) 데이터 로드/축소
    data = make_dummy_data(ws)

    # 3) Mass 범위/빈닝 설정과 mass-only 빠른 점검 (원본 drawInclusiveMassPlots 참조) :contentReference[oaicite:10]{index=10}
    ws.var("mass").setRange(massLow, massHigh)
    ws.var("mass").setBins(300)

    # 노이즈 최소화
    # for i in range(3):
    #     RooFit.RooMsgService.instance().getStream(i).removeTopic(RooFit.Eval)
    # for i in range(RooFit.RooMsgService.instance().numStreams()):
    #     RooFit.RooMsgService.instance().getStream(i).removeTopic(RooFit.Tracing)

    # 4) 사이드밴드/시그널 분리, ctErr 히스토그램, 스케일 팩터 계산 로직 (원본 computeScaleF) :contentReference[oaicite:11]{index=11}
    # mass PDF가 Chebychev일 때 계수(sl1, sl2...)를 allPars에 담아 넘기는 형태. 여기서는 예시 변수만 생성.
    ws.factory("sl1[0.0,-5,5]")
    ws.factory("sl2[0.0,-5,5]")
    allPars = RooArgList()
    allPars.add(ws.var("sl1"))
    allPars.add(ws.var("sl2"))

    try:
        scaleF = compute_scaleF(bkgMassOrder, massLow, massHigh, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi, allPars, "sl")
        print("scaleF:", scaleF)
    except Exception as e:
        print("scaleF computation skipped:", e)

    # 5) ct-배경 피팅 (Per-event-error 투영 포함 예시) :contentReference[oaicite:12]{index=12}:contentReference[oaicite:13]{index=13}
    # 실제로는 redDataSB (SB에서 추출), err 분포의 RooHistPdf 생성 후 ProjWData로 결합
    # 여기서는 구조만 보여주며 fit 호출
    ctRangeName = "ctRange"
    ws.var("ctau3D").setRange(ctRangeName, ctLow, ctHigh)
    ws.pdf("CtBkgTot")  # 존재 확인

    # 더미로 mass window cut 적용한 서브셋 생성(실분석 때는 TTree cut 사용)
    # SB: mass < 2.9 or mass > 3.3
    redDataSB = data.reduce(f"mass<{sig_lo} || mass>{sig_hi}")

    if redDataSB and redDataSB.numEntries() >= 0:
        fitBkg = ws.pdf("CtBkgTot").fitTo(
            redDataSB, RooFit.Save(), RooFit.PrintLevel(0),
            RooFit.PrintEvalErrors(-1), RooFit.NumCPU(4), RooFit.Range(ctRangeName),
            RooFit.RecoverFromUndefinedRegions(3)
        )
        if fitBkg:
            fitBkg.Print("v")

    # 6) 최종 2D 토탈 피팅(요약) :contentReference[oaicite:14]{index=14}
    # 실제 분석: dhData(RooDataHist), err 히스토그램과 함께 ProjWData 사용.
    # 여기서는 예시로 단순 fit
    tot = ws.pdf("totPDF")
    if tot:
        fit2D = tot.fitTo(data, RooFit.Save(), RooFit.PrintLevel(0), RooFit.PrintEvalErrors(-1), RooFit.NumCPU(4))
        if fit2D:
            fit2D.Print("v")

    # 7) 기본 플롯(간단) — 필요한 경우 원본의 상세 그리기 루틴 사용 :contentReference[oaicite:15]{index=15}
    can = TCanvas("c", "c", 800, 600)
    fr = ws.var("mass").frame()
    data.plotOn(fr)
    ws.pdf("MassPDF").plotOn(fr, RooFit.Normalization(1.0, RooFit.RooAbsReal.RelativeExpected))
    fr.Draw()
    can.SaveAs("mass_check.png")

    print("=== finish fit2D_step_by_step_py() ===")

# ----------------------------------------------------------
if __name__ == "__main__":
    fit2D_step_by_step_py()
