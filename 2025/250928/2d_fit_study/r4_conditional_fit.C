#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;

// Draw plot with kinematic cut
// ----------
double ptLow = 20, ptHigh = 25, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
double ctMin = -0.5, ctMax = 6;
double errMin = 0.01, errMax = 0.06;

// 0.016 - 0.03: 33초

RooAbsPdf *getOrMakeErrPdf(RooRealVar *ctau3DErr,
                           RooDataSet *dataset,
                           double ptLow, double ptHigh,
                           const char *prefix = "errPdf")
{
  // 파일 이름을 pT bin에 맞게 생성
  std::string filename = Form("%s_%.1f_%.1f.root", prefix, ptLow, ptHigh);

  // 1) 파일에서 불러오기 시도
  TFile *fIn = TFile::Open(filename.c_str(), "READ");
  if (fIn && !fIn->IsZombie())
  {
    RooAbsPdf *pdf = dynamic_cast<RooAbsPdf *>(fIn->Get("errPdf"));
    if (pdf)
    {
      std::cout << ">>> Loaded errPdf from file: " << filename << std::endl;
      fIn->Close();
      return pdf; // 성공적으로 불러온 경우
    }
    fIn->Close();
  }

  // 2) 없으면 새로 생성
  std::cout << ">>> Creating new errPdf and saving to " << filename << std::endl;

  // RooDataHist (저장용)
  auto histErr = new RooDataHist("histErr", "ctau3DErr hist", RooArgSet(*ctau3DErr), *dataset);

  // RooKeysPdf 생성
  RooKeysPdf *errPdf = new RooKeysPdf("errPdf",
                                      "Error PDF (KDE)",
                                      *ctau3DErr,
                                      *dataset,
                                      RooKeysPdf::MirrorBoth,
                                      2.0);

  // 3) 파일에 저장
  TFile *fOut = TFile::Open(filename.c_str(), "RECREATE");
  histErr->Write("histErr");
  errPdf->Write("errPdf");
  fOut->Close();

  return errPdf;
}

void fitWithDecay_PEE(RooDataSet *dataset1,
                      double rangeMin, double rangeMax,
                      // int nBinsErr,
                      const char *outname)
{
  // === Variables ===
  auto ctau3D = new RooRealVar("ctau3D", "ctau3D", ctMin, ctMax, "GeV/c^{2}");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", errMin, errMax, "mm");

  // RooRealVar *ctau3D = (RooRealVar *)dataset->get()->find("ctau3D");
  // RooRealVar *ctau3DErr = (RooRealVar *)dataset->get()->find("ctau3DErr");
  if (!ctau3D || !ctau3DErr)
  {
    std::cerr << "Error: ctau3D or ctau3DErr not found in dataset!" << std::endl;
    return;
  }
  RooArgSet obs(*ctau3D, *ctau3DErr);
  auto dataset = new RooDataSet("data1", "dataset with local vars", obs, Import(*dataset1));
  // ctau3D->setRange(rangeMin, rangeMax);
  // ctau3DErr->setRange(errMin, errMax);
  // ctau3D->setBins(180);
  // ctau3DErr->setBins(200);
  
  cerr << "Ping 1" << endl;

  // === 1) errPdf 생성 ===
  // RooAbsPdf *errPdf = getOrMakeErrPdf(ctau3DErr, dataset, ptLow, ptHigh, "errPdf_SIG");

  cout << "\nstart to build errPDF\n";
  ctau3DErr->setBins(100);
  RooDataHist *histCtErrSBL = new RooDataHist("histCtErrSBL",
                                              "ctau3DErr distribution (SBL)",
                                              RooArgSet(*ctau3DErr),
                                              *dataset);
  // ===== 빈 값 확인 및 교정 =====
int nBins = histCtErrSBL->numEntries(); // RooArgSet entry 개수
std::cout << "[INFO] nBins in histCtErrSBL = " << nBins << std::endl;

for (int i = 0; i < nBins; i++) {
  const RooArgSet* binCoords = histCtErrSBL->get(i);  // bin 좌표
  double w = histCtErrSBL->weight(*binCoords);        // bin 내용 (weight)

  if (w <= 0) {
    std::cout << "[WARNING] Bin " << i
              << " (" << ctau3DErr->getVal()
              << ") has non-positive value = " << w
              << " → replaced with 0.0001" << std::endl;

    histCtErrSBL->set(*binCoords, 0.0001); // bin 값 강제 수정
  }
}
  auto *errPdf = new RooHistPdf("errPdf",
                               "PDF of ctau3DErr (SBL)",
                               RooArgSet(*ctau3DErr),
                               *histCtErrSBL);

  auto *errPdfPlot = new RooHistPdf("errPdfPlot",
                                "PDF of ctau3DErr (SBL)",
                                RooArgSet(*ctau3DErr),
                                *histCtErrSBL, 2);

  TCanvas *cOverlap = new TCanvas("cOverlap_ctauErrSBL", "ctau3DErr DataHist + HistPdf", 800, 600);

  RooPlot *fOverlap = ctau3DErr->frame(Title("ctau3DErr (SBL): DataHist vs HistPdf"));
  // DataHist (marker)
  histCtErrSBL->plotOn(fOverlap, MarkerStyle(20), MarkerColor(kBlack));
  // HistPdf (red line)
  errPdf->plotOn(fOverlap, LineColor(kRed), LineWidth(2));

  fOverlap->Draw();
  fOverlap->SetMinimum(1e-3);
  cOverlap->SetLogy();
  cOverlap->SaveAs("ctauErrSBL_datahist_vs_histpdf.png");
  cout << "\nfinish building errPDF\n";

  // === 2) Double Gaussian Resolution with per-event error ===
  // RooRealVar mean("mean", "resolution mean", 0.0, -0.1, 0.1);
  RooRealVar mean("mean", "resolution mean", 0.0);
  RooRealVar sigma1("sigma1", "resolution sigma1", 0.08, 0.001, 1);

  // sigma 비율 (>=1로 제한 가능)
  RooRealVar sigmaRatio("sigmaRatio", "sigma2/sigma1 ratio", 2.0, 1.0, 50.0);

  // sigma2 = sigma1 * sigmaRatio
  RooFormulaVar sigma2("sigma2", "@0 * @1", RooArgList(sigma1, sigmaRatio));

  RooGaussModel g1("g1", "Gauss Res1", *ctau3D, mean, sigma1, *ctau3DErr);
  RooGaussModel g2("g2", "Gauss Res2", *ctau3D, mean, sigma2, *ctau3DErr);

  RooRealVar fG1("fG1", "fraction Gauss1", 0.8, 0.001, 1);
  RooAddModel resModel("resModel", "Double Gaussian Resolution",
                       RooArgList(g1, g2), RooArgList(fG1));

  // === 3) Decay components ===
// ===== 파라미터 정의 =====
RooRealVar tauL("tauL", "lifetime left", 0.3, 0.001, 0.4);

RooRealVar tauC1("tauC1", "lifetime central1", 0.05, 0.001, 0.1);
RooRealVar tauC2("tauC2", "lifetime central2", 0.09, 0.001, 0.5);

RooRealVar tauR("tauR", "lifetime right", 1.2, 0.01, 1);

// fraction (혼합 비율)
RooRealVar fC1("fC1", "frac central1", 0.3, 0.0, 1.0);   // decayC1 vs decayC2 분할
RooRealVar fC("fC", "frac total central", 0.55, 0.0, 1.0);
RooRealVar fL("fL", "frac left", 0.07, 0.0, 1.0);
// 나머지 비율은 1 - fC - fL → decayR 로 자동

// ===== decay components =====
RooDecay decayL ("decayL",  "Left decay", *ctau3D, tauL,  resModel, RooDecay::Flipped);
RooDecay decayC1("decayC1", "Central decay 1", *ctau3D, tauC1, resModel, RooDecay::DoubleSided);
RooDecay decayC2("decayC2", "Central decay 2", *ctau3D, tauC2, resModel, RooDecay::DoubleSided);
RooDecay decayR ("decayR",  "Right decay", *ctau3D, tauR, resModel, RooDecay::SingleSided);

// ===== 중앙 decay 2개 먼저 합치기 =====
RooAddPdf decayCentral("decayCentral", "combined central decays",
                       RooArgList(decayC1, decayC2),
                       RooArgList(fC1));  // fC1 : decayC1 비중, (1-fC1): decayC2 비중

// ===== 전체 합치기 =====
RooFormulaVar fR("fR", "1.0 - @0 - @1", RooArgList(fL, fC));

auto decayModel = new RooAddPdf("decayModel", "All decay components",
                                RooArgList(decayL, decayCentral, decayR),
                                RooArgList(fL, fC));


  cerr << "Ping 4" << endl;

  // === 4) Full model: decay ⊗ errPdf ===
  // RooProdPdf fullModel("fullModel", "decay * errPdf",
  //                      RooArgList(decayModel, errPdf),
  //                      RooFit::Conditional(RooArgSet(*ctau3DErr), RooArgSet()));
  RooProdPdf fullModel("fullModel", "PDF with PEE", *errPdf, Conditional(*decayModel, *ctau3D));
  RooProdPdf fullModelPlot("fullModelPlot", "PDF with PEE", *errPdfPlot, Conditional(*decayModel, *ctau3D));
  RooRealVar nNorm("nNorm", "signal yield", 10000, 0, 1e7);

  // 3) RooExtendPdf로 덮기
  // RooExtendPdf fullModel("fullModel",
  //                        "Extended model with yield",
  //                        fullModel_tmp,
  //                        nNorm);
  cerr << "Ping 5" << endl;

  // === 5) Fit ===
  // 6) 피팅
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  ctau3D->setRange("meanRange", ctMin, ctMax);
  ctau3DErr->setRange("meanRange", errMin, errMax);

  RooFitResult *fitRes = decayModel->fitTo(*dataset,
                                           Save(),
                                           //  Range("meanRange"),
                                           //  Extended(),
                                           //  ConditionalObservables(*ctau3DErr),
                                           Optimize(0),
                                           Strategy(2),
                                           Offset(),
                                           NumCPU(32),
                                           EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(1.5),
                                           RooFit::PrintEvalErrors(-1),
                                           RooFit::PrintLevel(-1));
  RooFitResult *fitResPlot = decayModel->fitTo(*dataset,
                                               Save(),
                                               //  Range("meanRange"),
                                               //  Extended(),
                                               //  ConditionalObservables(*ctau3DErr),
                                               Optimize(0),
                                               Strategy(2),
                                               Offset(),
                                               NumCPU(32),
                                               EvalBackend("legacy"), RooFit::RecoverFromUndefinedRegions(1.5),
                                               RooFit::PrintEvalErrors(-1),
                                               RooFit::PrintLevel(-1));

  // === 6) Plot ===
  TH1 *hh_model = fullModel.createHistogram("hh_model", *ctau3D, Binning(50), YVar(*ctau3DErr, Binning(50)));
  hh_model->SetLineColor(kBlue);

  // Make projection of data an dt
  // ctau3D->setBins(50);
  RooPlot *frame2 = ctau3D->frame(Title("Projection of model(dt|dterr) on dt"));
  
  
  // useProof=true
  
  RooDataSet *data_mom = (RooDataSet*)dataset->reduce(*ctau3DErr);
  decayModel->plotOn(frame2, ProjWData(*ctau3DErr, *dataset),  Normalization(dataset->sumEntries(), RooAbsReal::NumEvent));
  dataset->plotOn(frame2, Name("data"));


  // RooPlot 안에 그려진 RooHist 찾기
  RooHist* hData = (RooHist*) frame2->findObject("data", RooHist::Class());
  if (hData) {
    int nPoints = hData->GetN();  // 보이는 포인트 개수
    double sumY = 0.0;

    // RooHist는 TGraph 형태 → GetY() 배열을 통해 값 접근 가능
    const double* yVals = hData->GetY();
    for (int i = 0; i < nPoints; i++) {
      sumY += yVals[i];  // 각 bin의 이벤트 수 합산
    }
    
    std::cout << "보이는 포인트 수 = " << nPoints << std::endl;
    std::cout << "눈에 보이는 이벤트 총합 = " << sumY << std::endl;
    std::cout << "Dataset command로 얻은 이벤트 수 = " << dataset->sumEntries() << endl;
  }
    
  // --------------------------
  // ProjWData(*ctau3DErr, *dataset), 
  // , Normalization(1, RooAbsReal::Relative)
  // Normalization(data_mom->sumEntries(), RooAbsReal::NumEvent), 

  // fullModel.plotOn(frame2, NumCPU(16), 
      // ProjWData(*data_mom), 
      // RooFit::Precision(-1));
  // ProjWData(*data_mom),  

  // // Draw all frame2s on canvas
  TCanvas *c11 = new TCanvas("rf307_fullpereventerrors", "rf307_fullperventerrors", 800, 400);
  c11->Divide(2);
  c11->cd(1);
  gPad->SetLeftMargin(0.20);
  hh_model->GetZaxis()->SetTitleOffset(2.5);
  hh_model->Draw("surf");
  c11->cd(2);
  gPad->SetLeftMargin(0.15);
  frame2->GetYaxis()->SetTitleOffset(1.6);
  frame2->Draw();
  c11->SaveAs("conditional_test.png");

  // full plot
  // ---------

  ctau3D->setRange("plot", ctMin, ctMax);
  


  RooPlot *frame = ctau3D->frame(Title("Decay with per-event error (ctau3DErr included)"), Bins(100), Range("plot"));
  
  
  // decayModel->plotOn(frame, Components("decayL"), LineColor(kBlue), LineStyle(2), NumCPU(16), ProjWData(*ctau3DErr, *dataset), Normalization(dataset->sumEntries(), RooAbsReal::NumEvent));
  // decayModel->plotOn(frame, Components("decayC"), LineColor(kGreen + 2), LineStyle(2), NumCPU(16), ProjWData(*ctau3DErr, *dataset), Normalization(dataset->sumEntries(), RooAbsReal::NumEvent));
  // decayModel->plotOn(frame, Components("decayR"), LineColor(kMagenta), LineStyle(2), NumCPU(16), ProjWData(*ctau3DErr, *dataset), Normalization(dataset->sumEntries(), RooAbsReal::NumEvent));
  dataset->plotOn(frame, Name("data"));
  fullModelPlot.plotOn(frame, LineColor(kRed), LineWidth(2), Name("model"), NumCPU(16), Precision(-1)); // Normalization(80742, RooAbsReal::NumEvent)

  // RooArgSet projVars(*ctau3DErr);
  // // Specify the number of bins for the projected histogram
  // cout << "\nCreate Projected HIST\n";
  // RooAbsPdf *projectedHist = fullModel.createProjection(projVars);
  // cout << "\nDraw Projected HIST\n";
  // projectedHist->plotOn(frame, RooFit::LineColor(kRed), Name("model"));

  // decayModel->plotOn(frame, Components("decayCentral"), LineColor(kBlue), LineStyle(2), NumCPU(16), Precision(-1));
  // decayModel->plotOn(frame, Components("decayL"), LineColor(kGreen + 2), LineStyle(2), NumCPU(16), Precision(-1));
  // decayModel->plotOn(frame, Components("decayR"), LineColor(kMagenta), LineStyle(2), NumCPU(16), Precision(-1));

  

  double fracC1 = fC.getVal() * fC1.getVal(); // C1
  double fracC2 = fC.getVal() * (1 - fC1.getVal()); // C2
  double fracR = 1.0 - fL.getVal() - fC.getVal(); // Right
  auto proj = RooArgList(*ctau3DErr);

  // fullModel.plotOn(frame, Components("decayCentral"), LineColor(kBlue), LineStyle(2), NumCPU(16), Precision(-1));
  // fullModel.plotOn(frame, Components("decayL"), LineColor(kGreen + 2), LineStyle(2), NumCPU(16), Precision(-1));
  // fullModel.plotOn(frame, Components("decayR"), LineColor(kMagenta), LineStyle(2), NumCPU(16), Precision(-1));
  
  // decayModel->plotOn(frame, LineColor(kBlue), LineWidth(2), NumCPU(16), Normalization(80742, RooAbsReal::NumEvent), Precision(-1), ProjWData(*ctau3DErr, *dataset));
  // decayModel->plotOn(frame, LineColor(kGreen), LineStyle(kDashed), LineWidth(2), NumCPU(16), Precision(-1), ProjWData(*ctau3DErr, *dataset));

  frame->SetMinimum(0.01);
  
  
  

  // === Chi2 계산 ===
  double chi2 = frame->chiSquare("model", "data");
  std::cout << "Chi2/ndf = " << chi2 << std::endl;

  // === Pull plot ===
  RooHist *pullHist = frame->pullHist("data", "model");
  RooPlot *framePull = ctau3D->frame(Title("Pull distribution"));
  framePull->addPlotable(pullHist, "P");

  // === Canvas 분할 ===
  TCanvas *c = new TCanvas("c_fit_ctau3D_perErr", "", 800, 800);
  c->SetLogy();
  c->Divide(1, 2);

  // 위쪽: 메인 피팅
  c->cd(1);
  gPad->SetPad(0.0, 0.3, 1.0, 1.0);
  frame->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.65, 0.85, Form("#chi^{2}/ndf = %.2f", chi2));

  // 아래쪽: pull
  c->cd(2);
  gPad->SetPad(0.0, 0.0, 1.0, 0.3);
  framePull->SetTitle("");
  framePull->GetYaxis()->SetTitle("Pull");
  framePull->GetYaxis()->SetTitleSize(0.1);
  framePull->GetYaxis()->SetLabelSize(0.08);
  framePull->GetXaxis()->SetTitleSize(0.12);
  framePull->GetXaxis()->SetLabelSize(0.1);
  framePull->Draw();

  TLine *line = new TLine(rangeMin, 0, rangeMax, 0);
  line->SetLineStyle(2);
  line->Draw("same");

  // === Save ===
  string outpath = "figs/fit/" + string(outname) + ".png";
  string outpath2 = "figs/fit/" + string(outname) + ".pdf";
  c->SaveAs(outpath.c_str());
  c->SaveAs(outpath2.c_str());

  cout << "\n=== Fit result - General RooHistPdf ===\n";
  fitRes->Print("V");

  cout << "\n=== Fit result - Interpolation 2 RooHistPdf ===\n";
  fitResPlot->Print("V");

  delete c;
}

void r4_conditional_fit()
{
  gSystem->mkdir("figs/fit", true);

  // apply weight? -> later


  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");

  

  char reduceDS[3000], reduceDS_woCtErr[3000];
  // build cuts
  sprintf(reduceDS, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"

          //  && ctau3DRes >=%.3f && ctau3DRes < %.3f

          "&& (recoQQsign == 0)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errMin, errMax);

  RooDataSet *redData;
  redData = (RooDataSet *)data->reduce(reduceDS);

  // sideband
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");

  fitWithDecay_PEE(redDataSBL, ctMin, ctMax, "conditional_SBL");
}