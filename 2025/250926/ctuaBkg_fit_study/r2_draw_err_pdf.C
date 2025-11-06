#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;

RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr");
void makeHistPdfPlot(RooDataSet *dataset, const char *varName, double rangeMin, double rangeMax, int nBins, const char *pdfName, const char *pdfTitle, const char *outname, int order = 4);

void r2_draw_err_pdf()
{
  gSystem->mkdir("figs/err_pdf", true);

  // apply weight? -> later

  // --- Data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");

  // Draw plot with kinematic cut
  // ----------
  double ptLow = 6.5, ptHigh = 9, yLow = 0, yHigh = 1.6, massLow = 2.6, massHigh = 3.5;
  double ctMin = -1, ctMax = 4, errmin = 0.008, errmax = 0.3;

  char reduceDS[3000], reduceDS_woCtErr[3000];
  // build cuts
  sprintf(reduceDS, // no multiplicty case
          "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"

          "&& (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f && ctau3DErr >= %.3f && ctau3DErr < %.3f)"

          "&& (recoQQsign == 0)",
          ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctMin, ctMax, errmin, errmax);

  RooDataSet *redData;
  redData = (RooDataSet *)data->reduce(reduceDS);

  // sideband
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");

  // RooRealVar *ctau3D = (RooRealVar *)data->get()->find("ctau3D");
  // RooRealVar *ctau3DErr = (RooRealVar *)data->get()->find("ctau3DErr");

  // // SB를 하려면 mass핏을 해야한다 -> 다음에
  // // LSB만 쓸거면 scaling 안 해도 되니까 그냥 그릴 수 있다.

  // signal은 background를 안 빼준 상황이라 그리지 않기로.
  // makeHistPdfPlot(redDataSIG, "ctau3DErr", -0.1, 0.3, 60, "pdfCtErrSBL", "PDF of ctau3DErr (SBL)", "ctauErrPdf_SIG.png");
  makeHistPdfPlot(redDataSB, "ctau3DErr", -0.1, 0.3, 60, "pdfCtErrSBL", "PDF of ctau3DErr (SB)", "ctauErrPdf_SB.png");
  makeHistPdfPlot(redDataSBL, "ctau3DErr", -0.1, 0.3, 60, "pdfCtErrSBL", "PDF of ctau3DErr (SBL)", "ctauErrPdf_SBL.png");
  makeHistPdfPlot(redDataSBR, "ctau3DErr", -0.1, 0.3, 60, "pdfCtErrSBL", "PDF of ctau3DErr (SBㄲ)", "ctauErrPdf_SBR.png");
}

void makeHistPdfPlot(RooDataSet *dataset,
                     const char *varName,
                     double rangeMin, double rangeMax,
                     int nBins,
                     const char *pdfName,
                     const char *pdfTitle,
                     const char *outname,
                    int order)
{
  RooRealVar *var = (RooRealVar *)dataset->get()->find(varName);
  if (!var)
  {
    std::cerr << "Error: variable " << varName << " not found in dataset!" << std::endl;
    return;
  }

  var->setRange(rangeMin, rangeMax);
  // var->setBins(nBins);

  // make RooDataHist
  RooDataHist hist(Form("hist_%s", varName),
                   Form("%s distribution", varName),
                   RooArgSet(*var),
                   *dataset);

  // amke RooHistPdf
  RooHistPdf pdf(pdfName, pdfTitle, RooArgSet(*var), hist, order);

  // draw
  RooPlot *frame = var->frame(Title(Form("%s distribution", varName)));
  dataset->plotOn(frame, Name("data"), MarkerColor(kBlack), LineColor(kBlack));
  pdf.plotOn(frame, LineColor(kRed), LineWidth(2), Name("pdf"));

  TCanvas *c = new TCanvas(Form("c_%s", varName), "", 800, 600);
  frame->Draw();

  string outpath = "figs/err_pdf/" + string(outname);
  c->SaveAs(outpath.c_str());

  delete c;
}

RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr")
{

  if (binSIG->numEntries() != binSB->numEntries())
  {
    cout << "ERROR subtractSidebands : different binning!" << endl;
    return 0;
  }
  RooDataHist *binScaleBKG = new RooDataHist("binScaleBKG", "scaled SB", RooArgSet(*(ws->var(varName.c_str()))));

  //// **** bin-by-bin scaling
  const RooArgSet *argSIG;
  const RooArgSet *argSB;
  for (Int_t i = 0; i < binSIG->numEntries(); i++)
  {
    argSIG = binSIG->get(i);
    argSB = binSB->get(i);
    RooRealVar *thisVar = (RooRealVar *)argSIG->find(varName.c_str());
    ws->var(varName.c_str())->setVal(thisVar->getVal());
    //// *** set minimum as 0.1 to prevent empty PDF
    float wBkg = binSB->weight(*argSB, 0, false);
    if (wBkg <= 0.1)
      wBkg = 0.1;
    binScaleBKG->add(RooArgSet(*(ws->var(varName.c_str()))), scalefactor * wBkg);
    float newWeight = binSIG->weight(*argSIG, 0, false) - scalefactor * binSB->weight(*argSB, 0, false);
    if (newWeight <= 0.1)
      newWeight = 0.1;
    binSubtrSIG->add(RooArgSet(*(ws->var(varName.c_str()))), newWeight);
  }
  return binScaleBKG;
}