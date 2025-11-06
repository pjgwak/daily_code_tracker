#include <iostream>
#include <string>
#include <vector>
#include <TStopwatch.h>
#include <TFile.h>
#include <RooDataSet.h>
#include <TString.h>
#include <RooCategory.h>

using namespace RooFit;
using namespace std;

void mass_sim3()
{
  cout << "=== start make_sub_datasets() ===\n";
  TStopwatch time;
  time.Start();

  // input 가져오기
  double ptLow = 6.5, ptHigh = 7.5;
  vector<double> ctEdge = {-0.5, -0.2, -0.1, -0.07, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.2, 0.5, 1.0, 2.0};

  auto fname = [&](double lo, double hi)
  { return Form("roots/pt%.1f_%.1f_ctau%.2f_%.2f.root", ptLow, ptHigh, lo, hi); };

  for (size_t i = 0; i + 1 < ctEdge.size(); ++i)
  {
    auto f = TFile::Open(fname(ctEdge[i], ctEdge[i + 1]), "READ");
    auto ds = (RooDataSet *)f->Get("dataset");
    printf("%zu: %lld\n", i, ds ? ds->numEntries() : -1LL);
    f->Close();
  }

  // get data and buld a category
  RooCategory cat("cat", "");
  vector<TFile *> files;
  vector<RooDataSet *> datasets;
  for (size_t i = 0; i+1 < ctEdge.size(); ++i)
  {
    double lo = ctEdge[i], hi = ctEdge[i + 1];
    files.push_back(new TFile(fname(lo, hi), "READ"));
    datasets.push_back((RooDataSet *)files.back()->Get("dataset"));
    cat.defineType(Form("b%zu", i));
  }
  auto any = datasets.front();
  auto mass = (RooRealVar *)any->get()->find("mass");
  mass->setRange(2.6, 3.5);

  // 변수 세팅 + 모델 만들기
  // DCB+Gauss+2nd Cheby -> 초기값 저장하고 입력할 최적의 방법 제시


  // 피팅

  // 그림 그리기

  // --- build models ---

  // 모델(공유 파라미터)
  // RooRealVar mean("mean", "mean", 3.096, 3.05, 3.14);
  RooRealVar sigma("sigma", "", 0.03, 0.005, 0.1);
  RooRealVar aL("aL", "", 1.5, 0.5, 5), nL("nL", "", 3, 1, 20), aR("aR", "", 2, 0.5, 5), nR("nR", "", 3, 1, 20);
  RooRealVar s2("s2", "", 0.07, 0.01, 0.2);
  RooRealVar fSig("fSig", "", 0.7, 0, 1);


  // RooRealVar ns("ns", "ns", any->numEntries() * 0.5, 0, 1e9);
  // RooRealVar nb("nb", "nb", any->numEntries() * 0.5, 0, 1e9);
  // RooAddPdf model("model", "model", RooArgList(sig, bkg), RooArgList(ns, nb)); // extended

  // 동시 모델
  RooSimultaneous sim("sim", "sim", cat);
  std::vector<std::unique_ptr<RooRealVar>> mean, c1, c2, ns, nb;
  std::vector<std::unique_ptr<RooAbsPdf>> dcb, gaus, sig, bkg, model;

  for (size_t i = 0; i < datasets.size(); ++i)
  {
    double N = datasets[i]->numEntries();

    mean.emplace_back(std::make_unique<RooRealVar>(Form("mean_%zu", i), "", 3.096, 3.05, 3.14));
    dcb.emplace_back(std::make_unique<RooCrystalBall>(Form("dcb_%zu", i), "", *mass, *mean.back(), sigma, aL, nL, aR, nR));
    gaus.emplace_back(std::make_unique<RooGaussian>(Form("gaus_%zu", i), "", *mass, *mean.back(), s2));
    sig.emplace_back(std::make_unique<RooAddPdf>(Form("sig_%zu", i), "", RooArgList(*dcb.back(), *gaus.back()), RooArgList(fSig)));

    c1.emplace_back(std::make_unique<RooRealVar>(Form("c1_%zu", i), "", 0, -1, 1));
    c2.emplace_back(std::make_unique<RooRealVar>(Form("c2_%zu", i), "", 0, -1, 1));
    bkg.emplace_back(std::make_unique<RooChebychev>(Form("bkg_%zu", i), "", *mass, RooArgList(*c1.back(), *c2.back())));

    ns.emplace_back(std::make_unique<RooRealVar>(Form("ns_%zu", i), "", 0.6 * N, 0, 1e9));
    nb.emplace_back(std::make_unique<RooRealVar>(Form("nb_%zu", i), "", 0.4 * N, 0, 1e9));
    model.emplace_back(std::make_unique<RooAddPdf>(Form("m_%zu", i), "", RooArgList(*sig.back(), *bkg.back()), RooArgList(*ns.back(), *nb.back())));

    sim.addPdf(*model.back(), Form("b%zu", i));
  }

  // 합친 데이터(Import로 간단히)
  RooDataSet comb("comb", "comb", RooArgSet(*mass, cat));
  for (size_t i = 0; i + 1 < ctEdge.size(); ++i)
  {
    auto *ds = static_cast<RooDataSet *>(datasets[i]->reduce(SelectVars(RooArgSet(*mass)))); // mass 통일
    RooDataSet tmp(Form("tmp%zu", i), "", RooArgSet(*mass), Index(cat),
                   Import(Form("b%zu", i), *ds));
    comb.append(tmp);
    delete ds;
  }

  // --- fit --
  auto *fitResult = sim.fitTo(comb, Extended(), Save(), PrintLevel(-1), NumCPU(32), EvalBackend("legacy"));

  // --- plot ---
  int nbins = (int)datasets.size(); // 없으면 원하는 값으로 지정
  RooBinning bins(90, mass->getMin(), mass->getMax());
  int nParEffBin = 4; // 이 빈에서 뜨는 파라미터 수(ns, nb, c1, c2)

  for (int i = 0; i < nbins; ++i)
  {
    TCanvas c(Form("c_b%d", i), "", 700, 700);
    c.Divide(1, 2, 0, 0);

    // 상단: 데이터+모델
    c.cd(1);
    RooPlot *f = mass->frame();
    comb.plotOn(f, Cut(Form("cat==cat::b%d", i)), Name(Form("d%d", i)));
    sim.plotOn(f, Slice(cat, Form("b%d", i)), ProjWData(cat, comb),
               LineColor(kBlue), Name(Form("m%d", i)));
    double chi2ndf = f->chiSquare(Form("m%d", i), Form("d%d", i), nParEffBin);
    f->SetTitle(Form("bin b%d  (chi2/ndf = %.3f)", i, chi2ndf));
    f->Draw();

    // 하단: pull
    c.cd(2);
    auto hp = f->pullHist(Form("d%d", i), Form("m%d", i));
    RooPlot *fp = mass->frame();
    fp->addPlotable(hp, "P");
    fp->SetYTitle("pull");
    fp->SetMinimum(-5);
    fp->SetMaximum(5);
    fp->Draw();

    c.SaveAs(Form("fig_mass3/mass_pull_b%d.png", i));
  }
  // === finish main loop ===

  fitResult->Print();

  cout << "\n=== finish make_sub_datasets() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}