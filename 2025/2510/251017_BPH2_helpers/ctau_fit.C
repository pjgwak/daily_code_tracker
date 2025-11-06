#include <iostream>
#include <string>
#include <RooRealVar.h>
#include <TFile.h>
#include <RooDataSet.h>
#include <TStopwatch.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooCrystalBall.h>
#include <RooFitResult.h>
#include <RooChebychev.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TLatex.h>
#include <TAxis.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>
// #include <TLegend.h>

using std::cout; using std::string;
using namespace RooFit;

static const double CT_EDGES[] = {
    -0.5, -0.2, -0.1, -0.07, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.5, 1.0, 2.0}; //
static const int N_CT_BINS = (int)(sizeof(CT_EDGES) / sizeof(CT_EDGES[0])) - 1;

static RooRealVar *findNSigVar(const RooFitResult *fr)
{
  if (!fr)
    return nullptr;
  const RooArgList &fl = fr->floatParsFinal();
  const RooArgList &cl = fr->constPars();

  const char *candidates[] = {"nSig", "Nsig", "NSig", "nsig", "NSIG", "NSIG", "NSIG_"};
  for (auto name : candidates)
  {
    if (auto *a = fl.find(name))
      return dynamic_cast<RooRealVar *>(a);
    if (auto *a = cl.find(name))
      return dynamic_cast<RooRealVar *>(a);
  }
  // 마지막으로 Nsig(정확) 재시도
  if (auto *a = fl.find("Nsig"))
    return dynamic_cast<RooRealVar *>(a);
  if (auto *a = cl.find("Nsig"))
    return dynamic_cast<RooRealVar *>(a);
  return nullptr;
}

void ctau_fit(double ptLow=6.5, double ptHigh=7.5, double ctMin=0.0, double ctMax=0.01)
{
  cout << "=== start ctau_fit() ===\n";
  TStopwatch time;
  time.Start();

  TH1::SetDefaultSumw2(kTRUE);

  // --- set variables ---
  // kinematics
  float massLow = 2.6, massHigh = 3.5;
  float yLow = 0, yHigh = 2.4;
  int cLow = 0, cHigh = 180;

  // --- prepare histogram ---
  TH1D *hNSig = new TH1D("hNSig", Form("nSig vs ct;ct;N_{sig} (pT %.1f–%.1f)", ptLow, ptHigh), N_CT_BINS, CT_EDGES);

  // --- read inputs and fill the hist ---
  for (int i = 0; i < N_CT_BINS; ++i)
  {
    double ctMin = CT_EDGES[i];
    double ctMax = CT_EDGES[i + 1];

    char fname[512];
    std::snprintf(fname, sizeof(fname),
                  "roots/mass_fit_pT%.1f_%.1f_ct%.2f_%.2f.root",
                  ptLow, ptHigh, ctMin, ctMax);

    TFile f(fname, "READ");
    if (f.IsZombie())
    {
      std::cerr << "[WARN] Cannot open: " << fname << "\n";
      continue;
    }

    // 표준 이름: fitMass (없으면 대안 이름도 시도)
    RooFitResult *fr = nullptr;
    fr = dynamic_cast<RooFitResult *>(f.Get("fitMass"));
    if (!fr)
      fr = dynamic_cast<RooFitResult *>(f.Get("fitres"));
    if (!fr)
    {
      std::cerr << "[WARN] RooFitResult 'fitMass' not found in " << fname << "\n";
      continue;
    }

    RooRealVar *vN = findNSigVar(fr);
    if (!vN)
    {
      std::cerr << "[WARN] nSig/Nsig not found in " << fname << "\n";
      continue;
    }

    const double val = vN->getVal();
    const double err = vN->getError();

    // 히스토그램 bin index = i+1
    hNSig->SetBinContent(i + 1, val);
    hNSig->SetBinError(i + 1, err);
  }

  // --- build model ---
  cout << "\n=== Build model ===\n";
  
  // // --- fit ---
  // RooFitResult *fitMass = model.fitTo(*redData, Strategy(2), Range("massRange"), Save(), Extended(), PrintLevel(0), NumCPU(32), EvalBackend("legacy"));

  // // --- draw ---
  TCanvas c_hist("c_hist", "", 900, 600);
  hNSig->SetMarkerStyle(20);
  hNSig->SetMarkerSize(0.9);
  hNSig->SetLineWidth(2);
  hNSig->Draw("e");
  c_hist.SaveAs(Form("ctau_fit_pT%.1f_%.1f.png", ptLow, ptHigh));

  // === fitting ===
  ctMin = hNSig->GetXaxis()->GetXmin();
  ctMax = hNSig->GetXaxis()->GetXmax();
  RooRealVar ct("ct", "ct", -0.5, 2);
  hNSig->Sumw2();

  // 2) 데이터로 변환 (binned)
  RooDataHist dh("dh", "nSig vs ct (binned)", RooArgList(ct), hNSig);

  // 3) 해상도(Resolution) 모델 — 단일 가우시안(필요시 바꾸세요)
  RooRealVar mean("mean", "mean", 0.0, -0.02, 0.02);
  RooRealVar sig1("sig1", "core width", 0.030, 0.001, 1);
  RooRealVar sig12("sig12", "tail width", 1.1, 1.0, 5);
  RooRealVar sig23("si23", "tail width", 1.1, 1.0, 5);

  RooFormulaVar sig2("sig2", "@0*@1", RooArgList(sig1, sig12));
  RooFormulaVar sig3("sig3", "@0*@1", RooArgList(sig2, sig23));

  // 단일 가우스 해상도 두 개
  RooGaussModel g1("g1", "G1 res", ct, mean, sig1);
  RooGaussModel g2("g2", "G2 res", ct, mean, sig2);
  RooGaussModel g3("g3", "G2 res", ct, mean, sig3);

  // 혼합 비율 (g1의 비율; g2는 1-fG)
  RooRealVar fG("fG", "frac of g1", 0.7, 0.0, 1.0);
  RooRealVar fG2("fG2", "frac of g1", 0.7, 0.0, 1.0);

  // ★ 2가우스 해상도 결합: RooAddModel (해상도 전용 합성자)
  RooAddModel Res("Res", "G1+G2 resolution", RooArgList(g1, g2, g3), RooArgList(fG, fG2));

  RooGaussian Ga1("Ga1","Gaussian (sig1)", ct, mean, sig1);
  RooGaussian Ga2("Ga2","Gaussian (sig2)", ct, mean, sig2);

  // 두 가우스 혼합 PDF (해상도 합과 동일한 가중 fG)
  RooAddPdf G2("G2","Gauss1+Gauss2 PDF", RooArgList(Ga1, Ga2), RooArgList(fG));

  // Left/Right exponentials (해상도 없이 단면 지수)
  RooTruthModel TR("TR", "truth", ct);                                      // resolution 없음
  RooRealVar lambdaPos("lambdaPos", "lambda+", 0.5, 0.01, 1.0);               // >0
  RooRealVar lambdaNeg("lambdaNeg", "lambda-", 0.5, 0.01, 1.0);               // >0
  RooRealVar lambdadPos2("lambdadPos2", "lambda+", 0.5, 0.01, 1.0);
  RooDecay pos("pos", "right exp", ct, lambdaPos, Res, RooDecay::SingleSided); // ct>0
  RooDecay pos2("pos2", "right exp", ct, lambdadPos2, Res, RooDecay::SingleSided); // ct>0
  RooDecay neg("neg", "left  exp", ct, lambdaNeg, Res, RooDecay::Flipped);     // ct<0
  

  // --- 3) Extended yields ---
  const double Ntot = hNSig->Integral();
  RooRealVar bFrac("bFrac", "yield G1", 0.18, 0.15, 0.25);
  RooRealVar Npos("Npos", "yield +", 0.1, 0, 1);
  RooRealVar Nneg("Nneg", "yield -", 0.1, 0, 1);

  RooAddPdf npModel("npModel", "", {pos, neg, pos2}, {Npos, Nneg}, true);

  RooAddPdf model("model", "2G + Left/Right Exp",
                  RooArgList(npModel, Res),
                  RooArgList(bFrac));

  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }

  RooNumIntConfig &cfg = *RooAbsReal::defaultIntegratorConfig();
  cfg.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  cfg.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D")
      .setRealValue("epsRel", 1e-5);
  cfg.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 10000);

  RooFitResult *fr = model.fitTo(dh,
                                 Save(),
                                //  Extended(),
                                Optimize(0),
                                 PrintEvalErrors(-1),
                                 IntegrateBins(1e-6),
                                 PrintLevel(0),
                                 Strategy(2)
                                 );
  fr->Print("V");

  // TCanvas *c = new TCanvas("c_massfit", "mass fit", 900, 800);
  // RooPlot *fMass = mass->frame(Title("Mass fit: DCB + Cheby3"));
  // redData->plotOn(fMass, Name("h_data"), DataError(RooAbsData::SumW2));
  // model.plotOn(fMass, Name("c_model"), Range("massRange"), LineColor(kBlue), LineWidth(2));
  // fMass->Draw();
  // c->SaveAs(Form("figs/ctau_fit_pT%.1f_%.1f_ct%.2f_%.2f.png", ptLow, ptHigh, ctMin, ctMax));

  // fitMass->Print("V");

  // // --- save ---
  // TFile f(Form("roots/ctau_fit_pT%.1f_%.1f_ct%.2f_%.2f.root", ptLow, ptHigh, ctMin, ctMax), "RECREATE");
  // fitMass->Write("fitMass");
  // f.Close();

  RooPlot *f = ct.frame(Title("ct fit: RooDecay (pos/neg/double-sided)"));
  dh.plotOn(f, Name("h_data"), DataError(RooAbsData::SumW2));
  model.plotOn(f, Name("c_tot"), LineColor(kBlue), LineWidth(2), Precision(1e-6));

  // Canvas
  TCanvas c("c_ct", "ct fit", 900, 800);
  f->Draw();
  c.SaveAs("ct_fit_components.png");

  auto hpull = f->pullHist("h_data", "c_tot");
  RooPlot *fpull = ct.frame(Title("Pull"));
  fpull->addPlotable(hpull, "P");
  TCanvas cpull("cpull", "ct fit", 900, 800);
  fpull->Draw();
  cpull.SaveAs("ct_pull.png");

  // ------------------
  int npar = fr->floatParsFinal().getSize();                // 자유 파라미터 수
  double chi2ndf = f->chiSquare("c_tot", "h_data", npar);         // χ²/ndf
  std::cout << "chi2/ndf (frame) = " << chi2ndf << "\n";

  // 3-B) χ²/ndf (정석법: RooChi2Var 직접)
  // double chi2ndf = f->chiSquare("h_data", "c_model", npar);

  // 유효 bin 수 = dh.numEntries() (빈마다 1), 제약/고정 파라미터 있으면 npar 조정
  int ndf = dh.numEntries() - npar;
  // double chi2ndf2 = chi2val / std::max(ndf, 1);
  std::cout << "chi2 = " << chi2ndf << ", ndf = " << ndf;

  // double nNpropt = Ng1.getVal() + Nneg.getVal() + Npos.getVal();

  cout << "\nb-fraction: " << bFrac.getVal() << "\n";

  cout << "\n=== finish ctau_fit() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}