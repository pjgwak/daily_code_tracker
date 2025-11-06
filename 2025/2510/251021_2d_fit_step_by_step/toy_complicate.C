// file: toy_pee_hist_v2.C
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooGaussian.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooAddPdf.h>
#include <RooAddModel.h>
#include <RooProdPdf.h>
#include <RooChebychev.h>
#include <RooCBShape.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooPlot.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <memory>
using namespace RooFit;

void toy_complicate(int nTot=150000, double sigFrac=0.55)
{
  // 메시지 억제
  for (int i = 0; i < 3; i++)
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);

  // -------- Observables --------
  RooRealVar mass("mass","mass [GeV]", 2.6, 3.5);
  RooRealVar ctau3D("ctau3D","ctau3D [mm]", -0.5, 4.0);
  RooRealVar ctau3DErr("ctau3DErr","per-event ctau error [mm]", 0.005, 0.20);

  // -------- errPdf: signal / background 분리 --------
  TRandom3 rng(12345);

  // sig용 에러 템플릿: 조금 더 작은 에러 쪽으로 편향
  TH1D hErrSig("hErrSig","ctErr (sig)", 50, 0.005, 0.20);
  for (int i=0;i<300000;i++){
    double v = std::exp(rng.Gaus(std::log(0.035), 0.30)); // 평균 더 작게
    if (v<0.005) v=0.005; if (v>0.20) v=0.20; hErrSig.Fill(v);
  }
  RooDataHist hErrSigRDH("hErrSigRDH","", RooArgList(ctau3DErr), &hErrSig);
  RooHistPdf  errPdf_sig("errPdf_sig","p(ctErr|sig)", RooArgSet(ctau3DErr), hErrSigRDH);

  // bkg용 에러 템플릿: 더 퍼지게
  TH1D hErrBkg("hErrBkg","ctErr (bkg)", 50, 0.005, 0.20);
  for (int i=0;i<300000;i++){
    double v = std::exp(rng.Gaus(std::log(0.055), 0.40)); // 평균 더 크게
    if (v<0.005) v=0.005; if (v>0.20) v=0.20; hErrBkg.Fill(v);
  }
  RooDataHist hErrBkgRDH("hErrBkgRDH","", RooArgList(ctau3DErr), &hErrBkg);
  RooHistPdf  errPdf_bkg("errPdf_bkg","p(ctErr|bkg)", RooArgSet(ctau3DErr), hErrBkgRDH);

  // -------- ct resolution (PEE): 2-Gauss 모델 --------
  RooRealVar ctMean("ctMean","resolution mean", 0.0);
  RooRealVar one("one","scale", 1.0);

  RooRealVar sig1("sig1","core sigma", 0.028, 0.005, 0.20);
  RooGaussModel g1("g1","GaussModel1", ctau3D, ctMean, sig1, one, ctau3DErr);

  RooRealVar k("k","tail/core ratio", 2.0, 1.1, 8.0);
  RooFormulaVar sig2("sig2","@0*@1", RooArgList(sig1, k));
  RooGaussModel g2("g2","GaussModel2", ctau3D, ctMean, (RooAbsReal&)sig2, one, ctau3DErr);

  RooRealVar fG1("fG1","frac core", 0.75, 0.0, 1.0);
  // 해상도 합(ResolutionModel로 동작)
  RooAddModel ctRes("ctRes","2G PEE resolution", RooArgList(g1,g2), RooArgList(fG1));

  // -------- ct signal: pure resolution (prompt-like) --------
  // RooAddModel 자체가 PDF로도 동작하므로 그대로 사용
  RooRealVar tauS("tauS","tau(sym)", 0.6, 0.05, 10.0);
  RooDecay ctSig ("ctSig","DoubleSided",  ctau3D, tauS, ctRes, RooDecay::DoubleSided);
//   RooAbsPdf& ctSig = ctRes;

  // -------- ct background: RooDecay(Flip, Double, Single) with same resolution --------
  RooRealVar tauP("tauP","tau(+)", 1.6, 0.1, 10.0);
  RooRealVar tauM("tauM","tau(-)", 0.9, 0.05, 10.0);
  

  RooDecay ctSingle ("ctSingle","SingleSided",  ctau3D, tauP, ctRes, RooDecay::SingleSided);
  RooDecay ctFlipped("ctFlipped","Flipped",     ctau3D, tauM, ctRes, RooDecay::Flipped);

  RooRealVar f1("f1","frac Single", 0.40, 0.0, 1.0);
  RooRealVar f2("f2","frac Flipped",0.30, 0.0, 1.0); // 나머지는 Double
  RooAddPdf  ctBkg("ctBkg","bkg ctau",
                   RooArgList(ctSingle, ctFlipped),
                   RooArgList(f1)); // 3성분 SUM: f1*Single + f2*Flipped + (1-f1-f2)*Double

  // -------- ct × err 결합 (성분별 별도 errPdf) --------
  RooProdPdf ctSigModel("ctSigModel","ctSig * errSig",
                        RooArgSet(ctSig), Conditional(errPdf_sig, RooArgSet(ctau3DErr)));

  RooProdPdf ctBkgModel("ctBkgModel","ctBkg * errBkg",
                        RooArgSet(ctBkg), Conditional(errPdf_bkg, RooArgSet(ctau3DErr)));

  // -------- mass signal: CB + Gaussian --------
  RooRealVar mean("mean","mass mean", 3.096, 3.07, 3.12);

  RooRealVar sigCB("sigCB","CB sigma", 0.018, 0.004, 0.080);
  RooRealVar aL("aL","alphaL", 1.8, 0.2, 10.0);
  RooRealVar nL("nL","nL", 5.0, 0.5, 100.0);
  RooCBShape  CB("CB","CrystalBall", mass, mean, sigCB, aL, nL);

  RooRealVar sigG("sigG","G sigma", 0.030, 0.004, 0.120);
  RooGaussian G("G","Gaussian", mass, mean, sigG);

  RooRealVar fCB("fCB","frac(CB)", 0.65, 0.0, 1.0);
  RooAddPdf massSig("massSig","CB + G", RooArgList(CB, G), RooArgList(fCB));

  // -------- mass background: 2nd Chebychev --------
  RooRealVar c0("c0","c0", 0.1, -1.0, 1.0);
  RooRealVar c1("c1","c1",-0.2, -1.0, 1.0);
  RooRealVar c2("c2","c2", 0.1, -1.0, 1.0);
  RooChebychev massBkg("massBkg","Cheby2", mass, RooArgList(c0,c1,c2));

  // -------- signal / background 2D 모델 --------
  RooProdPdf sig2D("sig2D","massSig * (ctSig*errSig)", RooArgSet(massSig, ctSigModel));
  RooProdPdf bkg2D("bkg2D","massBkg * (ctBkg*errBkg)", RooArgSet(massBkg, ctBkgModel));

  // -------- Extended sum --------
  double nSig0 = nTot*sigFrac;
  double nBkg0 = nTot*(1.0 - sigFrac);
  RooRealVar nSig("nSig","signal yield", nSig0, 0, 1e9);
  RooRealVar nBkg("nBkg","background yield", nBkg0, 0, 1e9);

  RooAddPdf model("model","sig+bkg", RooArgList(sig2D, bkg2D), RooArgList(nSig, nBkg));

  // -------- Generate ToyMC --------
  std::unique_ptr<RooDataSet> data(model.generate(RooArgSet(mass, ctau3D, ctau3DErr), nTot));

  // -------- Fit --------
  std::unique_ptr<RooFitResult> fitRes(
    model.fitTo(*data, Save(), NumCPU(32), Extended(),
                ConditionalObservables(RooArgSet(ctau3DErr)),
                PrintEvalErrors(-1), PrintLevel(-1))
  );
  fitRes->Print("v");

  // -------- Plots (ProjWData: RooDataHist로 가속) --------
//    // -------- Quick plots --------
//   // 1) mass
//   TCanvas can1("can1","mass",800,600);
//   RooPlot* fMass = mass.frame(Title("mass"));
//   data->plotOn(fMass, Binning(80));
//   // mass는 ct, ctErr에 대해 적분 필요 → ProjWData로 데이터 전달
//   model.plotOn(fMass, ProjWData(RooArgSet(ctau3D, ctau3DErr), *data), LineColor(kBlue));
//   fMass->Draw();
//   can1.SaveAs("toy_mass_v2.png");

//   // 2) ctau3D
//   TCanvas can2("can2","ctau3D",800,600);
//   RooPlot* fCt = ctau3D.frame(Title("ctau3D"));
//   data->plotOn(fCt, Binning(80));
//   // ct는 mass, ctErr에 대해 적분 → ProjWData에 둘 다 명시
//   model.plotOn(fCt, ProjWData(RooArgSet(mass, ctau3DErr), *data), LineColor(kRed));
//   fCt->Draw();
//   can2.SaveAs("toy_ctau3D_v2.png");

//   // 3) ctau3DErr
//   TCanvas can3("can3","ctau3DErr",800,600);
//   RooPlot* fErr = ctau3DErr.frame(Title("ctau3DErr"));
//   data->plotOn(fErr, Binning(50));
//   // errPdf만 투영해서 비교
//   errPdf_bkg.plotOn(fErr, LineColor(kGreen+2));
//   fErr->Draw();
//   can3.SaveAs("toy_ctErr_v2.png");

  // ===================== 2D Overlay (Data vs Model) =====================
const int nx = 80, ny = 80;

// 1) 데이터 2D 히스토그램 (mass × ctau3D)
TH2D* hData2D = (TH2D*) data->createHistogram(
  "hData2D", mass, Binning(nx), YVar(ctau3D, Binning(ny))
);
hData2D->SetTitle("Data;mass;ctau3D");

// 2) 모델 2D 히스토그램 (ctau3DErr는 data로 평균 → per-bin 적분)
//    fillHistogram(&h2, RooArgList(vars), projData)
TH2D hModel2D("hModel2D","Model;mass;ctau3D",
              nx, mass.getMin(), mass.getMax(),
              ny, ctau3D.getMin(), ctau3D.getMax());
model.fillHistogram(&hModel2D, RooArgList(mass, ctau3D), data.get());

// 3) 모델 스케일을 데이터 총계에 맞춤
double sData  = hData2D->Integral();
double sModel = hModel2D.Integral();
if (sModel > 0.0) hModel2D.Scale(sData / sModel);

// 4) Overlay: Data는 COLZ, Model은 등고선
TCanvas c2D("c2D","2D overlay (Data + Model)", 900, 780);
hData2D->Draw("COLZ");
hModel2D.SetLineColor(kBlack);
hModel2D.SetLineWidth(2);
hModel2D.Draw("CONT3 SAME");
c2D.SaveAs("toy_2D_overlay.png");

// 5) 잔차 맵: (Data - Model) / sqrt(Model)
TH2D hResid("hResid","Residuals;(mass);(ctau3D)", nx,
            mass.getMin(), mass.getMax(),
            ny, ctau3D.getMin(), ctau3D.getMax());
hResid.SetTitle("Residuals: (Data-Model)/#sqrt{Model};mass;ctau3D");
for (int ix=1; ix<=nx; ++ix){
  for (int iy=1; iy<=ny; ++iy){
    double m = hModel2D.GetBinContent(ix,iy);
    double d = hData2D->GetBinContent(ix,iy);
    double r = (m>0.0) ? (d - m) / std::sqrt(m) : 0.0;
    hResid.SetBinContent(ix,iy, r);
  }
}
TCanvas cRes("cRes","2D residual", 900, 780);
hResid.Draw("COLZ");
cRes.SaveAs("toy_2D_residual.png");
}
