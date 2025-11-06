#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace RooFit;

void plotVar(RooDataSet *dataset, const char *varName, const char *title, double rangeMin, double rangeMax, const char *outname);
RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr");

void r1_draw_sb_plots()
{
  gSystem->mkdir("figs/sb", true);

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

  // SG, LSB, RSB, SB - ctau3D, ctau3DErr
  plotVar(redDataSIG, "ctau3D", "", -1, 2, "ctau3D_KineCut_SIG.png");
  plotVar(redDataSIG, "ctau3DErr", "", -0.1, 0.3, "ctauErr_KineCut_SIG.png");

  plotVar(redDataSB, "ctau3D", "", -1, 2, "ctau3D_KineCut_SB.png");
  plotVar(redDataSB, "ctau3DErr", "", -0.1, 0.3, "ctauErr_KineCut_SB.png");

  plotVar(redDataSBL, "ctau3D", "", -1, 2, "ctau3D_KineCut_SBL.png");
  plotVar(redDataSBL, "ctau3DErr", "", -0.1, 0.3, "ctauErr_KineCut_SBL.png");

  plotVar(redDataSBR, "ctau3D", "", -1, 2, "ctau3D_KineCut_SBR.png");
  plotVar(redDataSBR, "ctau3DErr", "", -0.1, 0.3, "ctauErr_KineCut_SBR.png");

  // SB를 하려면 mass핏을 해야한다 -> 다음에
  // // scaleF to scale down ctErr distribution in 2.9-3.3 GeV/c2
  // float bc;
  // bc = ws->var("coefExp")->getVal(); // expBkg
  // float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  // // RooDataSet(unbinned) to RooDataHist (binned)
  // RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr", "binDataCtErr", RooArgSet(*(ws->var("ctau3DErr"))), *redData);
  // RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  // RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  // // extract Sig and Bkg: (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  // RooDataHist *binSubtractedSIG, *binScaledBKG;
  // binSubtractedSIG = new RooDataHist("binSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  // binScaledBKG = subtractSidebands(ws, binSubtractedSIG, binDataCtErrSIG, binDataCtErrSB, scaleF, "ctau3DErr");

  // // error PDF
  // RooHistPdf errPdfSig("errPdfSig", "Error PDF signal", RooArgSet(*(ws->var("ctau3DErr"))), *binSubtractedSIG);
  // ws->import(errPdfSig);
  // //  RooHistPdf errPdfBkgRaw("errPdfBkg","Error PDF bkg before scaling",RooArgSet(*(ws->var("ctau3DErr"))),*binDataCtErrSB);  ws->import(errPdfBkg);
  // RooHistPdf errPdfBkg("errPdfBkg", "Error PDF bkg scaled", RooArgSet(*(ws->var("ctau3DErr"))), *binScaledBKG);
  // ws->import(errPdfBkg);
}

void plotVar(RooDataSet *dataset,
             const char *varName,
             const char *title,
             double rangeMin, double rangeMax,
             const char *outname)
{
  // get observable
  RooRealVar *var = (RooRealVar *)dataset->get()->find(varName);
  if (!var)
  {
    std::cerr << "Error: variable " << varName << " not found in dataset!" << std::endl;
    return;
  }

  // set range
  var->setRange(rangeMin, rangeMax);

  // draw
  RooPlot *frame = var->frame(Title(title));
  dataset->plotOn(frame, Name("data"));

  TCanvas *c = new TCanvas(Form("c_%s", var->GetName()), "", 800, 600);
  frame->Draw();
  string outpath = "figs/sb/" + string(outname);
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