#include "RooMCStudy.h"
#include "RooPlot.h"

using namespace RooFit;

void toy_mc_fit_test()
{
  // -------------------------------
  // 1. 관측 변수와 PDF 정의
  // -------------------------------
  RooRealVar x("x", "Observable", -10, 10);
  RooRealVar mean("mean", "Mean", 0, -1, 1);
  RooRealVar sigma("sigma", "Sigma", 1, 0.1, 5);
  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  // -------------------------------
  // 2. 실제 데이터셋 만들기 (예시)
  // -------------------------------
  RooDataSet *data = gauss.generate(x, 1000);

  // -------------------------------
  // 3. 피팅 (실제 데이터)
  // -------------------------------
  gauss.fitTo(*data);

  // -------------------------------
  // 4. Toy MC Study
  // -------------------------------
  RooMCStudy mcstudy(gauss, x,
                     FitModel(gauss),
                     Binned(kFALSE),    // unbinned fit
                     Silence(),         // 출력 억제
                     Extended(kFALSE)); // Extended 여부

  mcstudy.generateAndFit(200, 1000); // (nToys, eventsPerToy)

  // -------------------------------
  // 5. 결과 보기
  // -------------------------------
  // 파라미터 pull distribution
  RooPlot *frame1 = mcstudy.plotPull(mean, Bins(40), FrameRange(-5, 5));
  RooPlot *frame2 = mcstudy.plotPull(sigma, Bins(40), FrameRange(-5, 5));

  TCanvas *c = new TCanvas("c", "Toy MC study", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
  c->SaveAs("toyMC_pulls.png");
}
