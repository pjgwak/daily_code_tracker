#include "utils/mass_fit_utils.h"
#include "utils/sideband_utils.h"
#include "cfgs/pp_run2_pt12p15_y0_1p6.C"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

using namespace RooFit;
using std::cout;

double computeScaleF(const std::string& massBkgPdf, const RooArgList& allPars, double massMin, double massMax, double sbL_lo, double sbL_hi, double sig_lo, double sig_hi, double sbR_lo, double sbR_hi);

RooDataHist *subtractSbSig(RooWorkspace *ws, RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName);
RooDataHist *subtractSbBkg(RooWorkspace *ws, RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName);
std::pair<double, double> getYRange(RooPlot *fr, const char *hname = "data");
void drawCtErrDists(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG, double ctErrMin, double ctErrMax, const FitConfigSB &cfg, const char *outPrefix);
void drawCtauErrPdf(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG);

void sideband_extraction(const char *configFile = "cfgs/pp_run2_pt6p5_9_y0_1p6.C")
{
  gROOT->ProcessLineSync(Form(".x %s", configFile));
  FitConfigSB cfg = getConfigSB();

  cout << "\n === Start sideband_extraction() ===\n";
  cout << "Config: " << configFile << "\n"
       << "System: " << cfg.system << "\n"
       << "Mode: " << cfg.mode << "\n"
       << "Prompt: " << (cfg.isPrompt ? "true" : "false") << "\n";


  TStopwatch t; t.Start();

  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(WARNING);

  // use rootlogon
  gROOT->Macro(cfg.rootlogon.c_str());

  // make output dirs
  gSystem->mkdir(cfg.figDir.c_str(), true);
  gSystem->mkdir(cfg.rootDir.c_str(), true);
  std::cout << "[Info] Output dirs prepared: " << cfg.figDir << ", " << cfg.rootDir << "\n";

  // read input
  TFile *fInput = TFile::Open(cfg.inputFile.c_str());
  if (!fInput || fInput->IsZombie())
  {
    cout << "[Error] Cannot open input file\n";
    return;
  }
  // OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root
  // OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root
  // OniaRooDataSet_JPsi_pp_GENONLY_NonPrompt_230215.root


  // read dataset
  RooDataSet *dsRaw = dynamic_cast<RooDataSet *>(fInput->Get("dataset"));
  dsRaw->Print();
  if (!dsRaw)
  {
    cout << "[Error]: Cannot find RooDataSet\n";
    return;
  }

  RooDataSet *dsWeight = new RooDataSet("dsWeight", "dataset with weight", dsRaw, *dsRaw->get(), 0, "weight");
  dsWeight->Print();

  // define cut
  TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )"; // 2018 acceptance cut
  TString kineCut = Form( "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
      " mass >= %.3f && mass < %.3f)", cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh, cfg.massLow, cfg.massHigh);

  // TString kineCut = Form(
  //     "(pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && "
  //     " mass >= %.3f && mass < %.3f && cBin >= %d && cBin < %d)",
  //     ptLow, ptHigh, yLow, yHigh, massLow, massHigh, cLow, cHigh);

  // opposite sign -> Must be FLowSkim later
  TString osCut = "(recoQQsign == 0)";

  TString fullCut = Form("%s && %s && %s",
                         osCut.Data(),
                         accCut.Data(),
                         kineCut.Data());

  // === new dataset with cuts ===
  RooDataSet *dsReduced_tmp = (RooDataSet *)dsWeight->reduce(Cut(fullCut));
  if (!dsReduced_tmp || dsReduced_tmp->numEntries() == 0)
  {
    cout << "[Error] Reduced dataset is empty. Cut = " << fullCut << "\n";
    return;
  }

  // observable
  auto mass = new RooRealVar("mass", "invariant mass", cfg.massLow , cfg.massHigh, "GeV/c^{2}");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", cfg.ctErrLow , cfg.ctErrHigh, "#sigma_{L_{J/#psi}} [mm]");
  auto weight = new RooRealVar("weight", "", 0, 10000, "");
  RooArgSet obs(*mass, *ctau3DErr, *weight);
  auto dsReduced = new RooDataSet("dsReduced", "dataset with local vars", obs, Import(*dsReduced_tmp));
  cout << "\n--- Objects in reduced dataset ---\n";
  dsReduced->Print();

  // check weight
  const bool hasWeight = (dsReduced->isWeighted() || dsReduced->weightVar() != nullptr);
  // const bool hasWeight = false;

  if (hasWeight)
  {
    // MC's weight value = 1 -> no-weight case
    cout << "[Info] Using weighted dataset \n";
  }
  else
  {
    cout << "[Info] Using UN-weighted dataset \n";
  }

  // === set ctau3DErr range ===
  // getCtErrRange(data,reduceDS_woCtErr, -ctMin, ctMax, &errmin, &errmax);

  // === Create workspace ===
  auto ws = new RooWorkspace("ws");
  ws->import(*dsReduced);
  // ws->var("ctau3DErr")->setBins(120);

  // Make subrange Dataset
  RooDataSet *redDataSIG = (RooDataSet *)dsReduced->reduce("mass>2.9 && mass<3.3");
  RooDataSet *redDataSB = (RooDataSet *)dsReduced->reduce("mass<2.9 || mass>3.3");
  RooDataSet *redDataSBL = (RooDataSet *)dsReduced->reduce("mass<2.9");
  RooDataSet *redDataSBR = (RooDataSet *)dsReduced->reduce("mass>3.3");
  

  // calculate sideband scale factor
  // scaleF to scale down ctErr distribution in 2.9-3.3 GeV/c2
  double scaleF = 0;
  {
    // bring mass fit result
    TFile fin(Form("%s/mass_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "read");
    if (fin.IsZombie())
    {
      cerr << "[Error] Can't open MC fit file\n";
      return;
    }
    RooFitResult *fr = (RooFitResult *)fin.Get("fitResult");
    if (!fr)
    {
      cerr << "[Error] There is no fitResult\n";
      return;
    }
    
    const RooArgList &floats = fr->floatParsFinal();
    const RooArgList &consts = fr->constPars();
    RooArgList allPars(floats);
    allPars.add(consts);

    // compute
    scaleF = computeScaleF(cfg.massBkgPdf, allPars,
                              cfg.massLow, cfg.massHigh,
                              cfg.sblLow, cfg.sblHigh,
                              cfg.sigLow, cfg.sigHigh,
                              cfg.sbrLow, cfg.sbrHigh);
  }
  cout << "scaleF: " << scaleF << "\n"; // TODO: 나중에 오리지널 코드 결과와 비교하기

  // === build Dataset and RooHIstPDf ===
  // RooDataSet(unbinned) to RooDataHist (binned)
  // RooDataHist *binDataCtErr = new RooDataHist("binDataCtErr", "binDataCtErr", RooArgSet(*(ws->var("ctau3DErr"))), *dsReduced);
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  // extract Sig and Bkg: (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  // binScaledBKG is not scaled
  RooDataHist *binSubtractedSIG, *binScaledBKG;
  binSubtractedSIG = new RooDataHist("binSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  binScaledBKG = new RooDataHist("binScaledBKG", "", RooArgSet(*(ws->var("ctau3DErr"))));
  // subtractSidebands(ws, binSubtractedSIG, binScaledBKG, binDataCtErrSIG, binDataCtErrSB, scaleF, "ctau3DErr");

  binSubtractedSIG = subtractSbSig(ws, binDataCtErrSIG, binDataCtErrSB, scaleF, "ctau3DErr");
  binScaledBKG = subtractSbBkg(ws, binDataCtErrSIG, binDataCtErrSB, 1, "ctau3DErr");

  auto errPdfSig = new RooHistPdf("errPdfSig", "Error PDF signal", RooArgSet(*(ws->var("ctau3DErr"))), *binSubtractedSIG);
  ws->import(*errPdfSig);
  //  RooHistPdf errPdfBkgRaw("errPdfBkg","Error PDF bkg before scaling",RooArgSet(*(ws->var("ctau3DErr"))),*binDataCtErrSB);  ws->import(errPdfBkg);
  auto errPdfBkg = new RooHistPdf("errPdfBkg", "Error PDF bkg scaled", RooArgSet(*(ws->var("ctau3DErr"))), *binScaledBKG);
  ws->import(*errPdfBkg);

  auto errPdfSigInter2 = new RooHistPdf("errPdfSigInter2", "Error PDF signal", RooArgSet(*(ws->var("ctau3DErr"))), *binSubtractedSIG, 2);
  auto errPdfBkgInter2 = new RooHistPdf("errPdfBkgInter2", "Error PDF bkg scaled", RooArgSet(*(ws->var("ctau3DErr"))), *binScaledBKG, 2);

  // === draw plots === 
  drawCtErrDists(ws, binDataCtErrSB, binDataCtErrSIG, binSubtractedSIG, binScaledBKG, cfg.ctErrLow, cfg.ctErrHigh, cfg, "ctau_err");
  // drawCtauErrPdf(ws, binDataCtErrSB, binDataCtErrSIG, binSubtractedSIG, binScaledBKG);

  // // === save results ===
  cout << "\n--- save results ---\n";
  TFile fout(Form("%s/ctau_err_pT%.1f_%.1f_y%.1f_%.1f.root", cfg.rootDir.c_str(), cfg.ptLow, cfg.ptHigh, cfg.yLow, cfg.yHigh), "RECREATE");
  errPdfSig->Write();
  errPdfBkg->Write();
  errPdfSigInter2->Write();
  errPdfBkgInter2->Write();
  fout.Close();

  // fitResult->Print("V");

  cout << "=== Done ===\n";
  t.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", t.RealTime(), t.CpuTime());
}


// === Helper functions ===
// ------------------------
std::pair<double,double> getYRange(RooPlot *fr, const char *hname)
{
  // find proper y ranges for log scale
  double ymin=1e300, ymax=-1e300;
  if (auto *hdata = dynamic_cast<RooHist*>(fr->getHist(hname))) {
    for (int i=0;i<hdata->GetN();++i){
      double x,y;
      hdata->GetPoint(i,x,y);
      if (y>0 && y<ymin) ymin=y;
      if (y>ymax) ymax=y;
    }
  }
  if (ymin<=0 || ymin==1e300) ymin=1e-3;
  return {ymin,ymax};
}

void drawCtauErrPdf(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG)
{
  RooPlot *errframe = ws->var("ctau3DErr")->frame();

  binDataCtErrSB->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kBlue), MarkerColor(kBlue), MarkerStyle(kOpenCircle));
  ws->pdf("errPdfBkg")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binDataCtErrSB->sumEntries(), RooAbsReal::NumEvent));
  // ws->pdf("errPdfBkg")->plotOn(errframe);

  TCanvas c0;
  string titlestr;
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = "_CtErrPdfBkg_Log.pdf";
  c0.SaveAs(titlestr.c_str());
  delete errframe;

  errframe = ws->var("ctau3DErr")->frame();
  binSubtractedSIG->plotOn(errframe, DataError(RooAbsData::SumW2), DrawOption("F"), FillColor(kWhite), LineColor(kWhite)); // just for axis
  ws->pdf("errPdfSig")->plotOn(errframe, LineColor(kViolet + 3), Normalization(binSubtractedSIG->sumEntries(), RooAbsReal::NumEvent));
  binDataCtErrSIG->plotOn(errframe, DataError(RooAbsData::SumW2), LineColor(kRed), MarkerColor(kRed), MarkerStyle(kOpenCircle));
  c0.Clear();
  c0.SetLogy(1);
  errframe->Draw();
  errframe->GetXaxis()->CenterTitle(1);
  errframe->GetYaxis()->CenterTitle(1);
  errframe->SetMinimum(0.01);
  titlestr = "_CtErrPdfSig_Log.pdf";
  c0.SaveAs(titlestr.c_str());
  delete errframe;
}

void drawCtErrDists(RooWorkspace *ws, RooDataHist *binDataCtErrSB, RooDataHist *binDataCtErrSIG, RooDataHist *binSubtractedSIG, RooDataHist *binScaledBKG, double ctErrMin, double ctErrMax, const FitConfigSB &cfg, const char *outPrefix)
{
  auto *ctErr = ws->var("ctau3DErr");
  ctErr->setRange(ctErrMin, ctErrMax);

  // SB
  // {
  //   TCanvas c("c_sb", "c_sb", 800, 600);
  //   c.SetLogy();
  //   RooPlot *fr = ctErr->frame(Range(ctErrMin, ctErrMax), Title("Sideband ctauErr"));
  //   binDataCtErrSB->plotOn(fr, DataError(RooAbsData::SumW2), Name("dataSB"));
  //   fr->GetXaxis()->SetTitle("ctau3DErr");
  //   fr->GetYaxis()->SetTitle("Events");
  //   fr->SetMinimum(0.1);
  //   fr->Draw();

  //   TString out = Form("%s/%s_SB", cfg.figDir.c_str(), outPrefix);
  //   c.SaveAs(out + ".png");
  //   c.SaveAs(out + ".pdf");
  // }

  // SIG
  // {
  //   TCanvas c("c_sig", "c_sig", 800, 600);
  //   c.SetLogy();
  //   RooPlot *fr = ctErr->frame(Range(ctErrMin, ctErrMax), Title("Signal ctauErr"));
  //   binDataCtErrSIG->plotOn(fr, DataError(RooAbsData::SumW2), Name("dataSIG"));
  //   fr->GetXaxis()->SetTitle("ctau3DErr");
  //   fr->GetYaxis()->SetTitle("Events");
  //   fr->SetMinimum(0.1);
  //   fr->Draw();

  //   TString out = Form("%s/%s_SIG", cfg.figDir.c_str(), outPrefix);
  //   c.SaveAs(out + ".png");
  //   c.SaveAs(out + ".pdf");
  // }

  // Subtracted SIG
  {
    TCanvas c("c_sub", "c_sub", 800, 800);
    TPad pad1("pad1","pad1",0.0,0.25,1.0,1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw(); pad1.cd();

    RooPlot *fr = ctErr->frame(Range(ctErrMin, ctErrMax), Title("Subtracted SIG vs Pdf"));
    binSubtractedSIG->plotOn(fr, Name("subSIG"), DataError(RooAbsData::SumW2));
    ws->pdf("errPdfSig")->plotOn(fr, LineColor(kRed), Name("pdfSig"));

    // set y range
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("subSIG")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }
    if (ymin <= 0 || ymin == 1e300)
      ymin = 1e-3;
    fr->SetMaximum(std::max(ymax, ymin) * 1e6);
    fr->SetMinimum(0.1);
    fr->Draw();

    TLegend leg(0.65,0.75,0.88,0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(fr->findObject("subSIG"), "Subtracted SIG", "lep");
    leg.AddEntry(fr->findObject("pdfSig"), "HistPdf SIG", "l");
    leg.Draw();

    // pull pad
    c.cd();
    TPad pad2("pad2","pad2",0.0,0.0,1.0,0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.35);
    pad2.Draw(); pad2.cd();

    RooHist *hpull = fr->pullHist("subSIG","pdfSig");
    RooPlot *fpull = ctErr->frame(Range(ctErrMin,ctErrMax),Title(""));
    fpull->addPlotable(hpull,"P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr");
    fpull->SetMinimum(-8); fpull->SetMaximum(8);
    fpull->Draw();

    TLine line(ctErrMin,0.0,ctErrMax,0.0);
    line.SetLineStyle(2); line.Draw("same");

    TString out = Form("%s/%s_Subtracted", cfg.figDir.c_str(), outPrefix);
    c.SaveAs(out + ".png");
    c.SaveAs(out + ".pdf");
  }

  // Bkg
  {
    TCanvas c("c_bkg", "c_bkg", 800, 800);
    TPad pad1("pad1","pad1",0.0,0.25,1.0,1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw(); pad1.cd();

    RooPlot *fr = ctErr->frame(Title("Bkg vs Pdf"));
    binScaledBKG->plotOn(fr, Name("bkg"));
    ws->pdf("errPdfBkg")->plotOn(fr, LineColor(kBlue), Name("pdfBkg"));
    // set y range
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("bkg")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double x, y;
        hdata->GetPoint(i, x, y);
        if (y > 0 && y < ymin)
          ymin = y;
        if (y > ymax)
          ymax = y;
      }
    }
    if (ymin <= 0 || ymin == 1e300)
      ymin = 1e-3;
    fr->SetMaximum(std::max(ymax, ymin) * 1e6);
    fr->SetMinimum(0.1);
    fr->Draw();

    TLegend leg(0.65,0.75,0.88,0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(fr->findObject("bkg"), "Bkg (raw SB)", "lep");
    leg.AddEntry(fr->findObject("pdfBkg"), "HistPdf Bkg", "l");
    leg.Draw();

    // pull pad
    c.cd();
    TPad pad2("pad2","pad2",0.0,0.0,1.0,0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.35);
    pad2.Draw(); pad2.cd();

    RooHist *hpull = fr->pullHist("bkg","pdfBkg");
    RooPlot *fpull = ctErr->frame(Range(ctErrMin,ctErrMax),Title(""));
    fpull->addPlotable(hpull,"P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr");
    fpull->SetMinimum(-8); fpull->SetMaximum(8);
    fpull->Draw();

    TLine line(ctErrMin,0.0,ctErrMax,0.0);
    line.SetLineStyle(2); line.Draw("same");

    TString out = Form("%s/%s_Bkg", cfg.figDir.c_str(), outPrefix);
    c.SaveAs(out + ".png");
    c.SaveAs(out + ".pdf");
  }
}

double computeScaleF(const std::string& massBkgPdf, const RooArgList& allPars, double massMin, double massMax, double sbL_lo, double sbL_hi, double sig_lo, double sig_hi, double sbR_lo, double sbR_hi)
{
  // === inline helper functions ===
  auto getVal = [&](const char* name)->double {
    if (auto* obj = allPars.find(name)) return static_cast<RooAbsReal*>(obj)->getVal();
    return 0.0;
  };
  auto parseChebyOrder = [&](const std::string& s)->int {
    if (s.rfind("cheby", 0) != 0) return 0;
    return std::atoi(s.c_str() + 5); // extract order from "chebyN" → N
  };
  auto Iexp = [&](double bc, double a, double b)->double {
    // fore expo
    return (std::exp(bc*b) - std::exp(bc*a)) / bc; // ∫ e^{bc m} dm
  };
  auto toUnitX = [&](double m)->double {
    // for chebychev
    return 2.0*(m - massMin)/(massMax - massMin) - 1.0; // [massMin,massMax]→[-1,1]
  };
  auto chebT = [&](int n, double x)->double {
    if (n==0) return 1.0; if (n==1) return x;
    double Tn_2=1.0, Tn_1=x, Tn=0.0;
    for (int k=2;k<=n;++k){ Tn = 2.0*x*Tn_1 - Tn_2; Tn_2=Tn_1; Tn_1=Tn; }
    return Tn;
  };
  auto chebVal_roofit = [&](double m, const std::vector<double>& c)->double {
    // f(m) = 1*T0 + c1*T1 + ... + cN*TN  (a0=1)
    const double x = toUnitX(m);
    double v = 1.0;
    for (int i=1;i<=static_cast<int>(c.size());++i) v += c[i-1]*chebT(i,x);
    return v;
  };
  auto integrateSimpson = [&](auto&& f, double a, double b, int nSeg=2048)->double {
    if (nSeg%2) ++nSeg;
    const double h=(b-a)/nSeg;
    double s=f(a)+f(b);
    for (int i=1;i<nSeg;++i){
      const double x=a+i*h; s += f(x)*((i%2)?4.0:2.0);
    }
    return s*(h/3.0);
  };

  // === compute ===
  if (massBkgPdf == "expo") {
    const double bc    = getVal("tauMass");
    const double I_sig = Iexp(bc, sig_lo, sig_hi);
    const double I_sb  = Iexp(bc, sbL_lo, sbL_hi) + Iexp(bc, sbR_lo, sbR_hi);
    return I_sig / I_sb;
  } else {
    const int order = parseChebyOrder(massBkgPdf); // cheby 1 ~ 6
    std::vector<double> c; c.reserve(order);
    for (int i=1;i<=order;++i){
      char nm[8]; std::snprintf(nm,sizeof(nm),"s%d",i);
      c.push_back(getVal(nm));
    }
    auto f = [&](double m){ return chebVal_roofit(m, c); };
    const double I_sig = integrateSimpson(f, sig_lo, sig_hi);
    const double I_sb  = integrateSimpson(f, sbL_lo, sbL_hi) + integrateSimpson(f, sbR_lo, sbR_hi);
    return I_sig / I_sb;
  }
}

// orginal subtraction code
RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName)
{
  if (binSIG->numEntries() != binSB->numEntries())
  {
    cout << "[Error] subtractSidebands : different binning!\n";
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
    binScaleBKG->add(RooArgSet(*(ws->var(varName.c_str()))), wBkg);

    float newWeight = binSIG->weight(*argSIG, 0, false) - scalefactor * binSB->weight(*argSB, 0, false);
    if (newWeight <= 0.1)
      newWeight = 0.1;
    binSubtrSIG->add(RooArgSet(*(ws->var(varName.c_str()))), newWeight);
  }

  return binScaleBKG;
}


RooDataHist *subtractSbSig(RooWorkspace *ws,RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName)
{
  if (!ws || !ws->var(varName.c_str()))
  {
    std::cout << "[Error] subtractSbSig: variable " << varName << " not found\n";
    return nullptr;
  }
  if (!binSIG || !binSB)
  {
    std::cout << "[Error] subtractSbSig: null input RooDataHist\n";
    return nullptr;
  }

  TH1 *hSIG = binSIG->createHistogram("hSIG_tmp", *ws->var(varName.c_str()));
  TH1 *hSB = binSB->createHistogram("hSB_tmp", *ws->var(varName.c_str()));

  hSIG->Sumw2();
  hSB->Sumw2();

  int nBins = hSIG->GetNbinsX();
  double xMin = hSIG->GetXaxis()->GetXmin();
  double xMax = hSIG->GetXaxis()->GetXmax();
  TH1D *hSub = new TH1D("hSubtracted", "Subtracted SIG", nBins, xMin, xMax);

  for (int i = 1; i <= nBins; i++)
  {
    double Nsig = hSIG->GetBinContent(i);
    double Esig = hSIG->GetBinError(i);

    double Nsb = hSB->GetBinContent(i);
    double Esb = hSB->GetBinError(i);

    // subtraction
    double content = Nsig - scalefactor * Nsb;
    // prevent empty bin
    if (content < 0.01)
      content = 0.01;
    double error2 = Esig * Esig + (scalefactor * Esb) * (scalefactor * Esb);

    hSub->SetBinContent(i, content);
    hSub->SetBinError(i, std::sqrt(error2));
  }

  auto binSubtrSIG = new RooDataHist("binSubtrSIG", "Subtracted SIG",
                                             RooArgList(*ws->var(varName.c_str())), hSub);

  // cleanup
  delete hSIG;
  delete hSB;
  delete hSub;

  return binSubtrSIG;
}

RooDataHist *subtractSbBkg(RooWorkspace *ws, RooDataHist *binSIG, RooDataHist *binSB, double scalefactor, string varName)
{
  if (!ws || !ws->var(varName.c_str()))
  {
    std::cout << "[Error] subtractSbBkg: variable " << varName << " not found\n";
    return nullptr;
  }
  if (!binSIG || !binSB)
  {
    std::cout << "[Error] subtractSbBkg: null input RooDataHist\n";
    return nullptr;
  }

  TH1 *hSIG = binSIG->createHistogram("hSIG_tmp", *ws->var(varName.c_str()));
  TH1 *hSB = binSB->createHistogram("hSB_tmp", *ws->var(varName.c_str()));

  hSIG->Sumw2();
  hSB->Sumw2();

  int nBins = hSB->GetNbinsX();
  double xMin = hSB->GetXaxis()->GetXmin();
  double xMax = hSB->GetXaxis()->GetXmax();
  TH1D *hBkg = new TH1D("hSubtracted", "scaled bkg", nBins, xMin, xMax);

  for (int i = 1; i <= hSB->GetNbinsX(); i++)
  {
    double Nsb = hSB->GetBinContent(i);
    double Esb = hSB->GetBinError(i);

    // background scaling
    double Nbkg = scalefactor * Nsb;
    // prevent empty bin
    if (Nbkg < 0.01)
      Nbkg = 0.01;
    double Ebkg = scalefactor * Esb;

    // 결과 저장
    hBkg->SetBinContent(i, Nbkg);
    hBkg->SetBinError(i, Ebkg);
  }

  auto binScaleBKG = new RooDataHist("binScaleBKG", "scaled bkg",
                                     RooArgList(*ws->var(varName.c_str())), hBkg);

  // cleanup
  delete hSIG;
  delete hSB;
  delete hBkg;

  return binScaleBKG;
}

void getCtErrRange(RooDataSet *data, const char *t_reduceDS_woCtErr, float lmin, float lmax, float *errmin, float *errmax)
{
  RooWorkspace *ws = new RooWorkspace("ctauerrorcheckWS");
  RooDataSet *redDataCut = (RooDataSet *)data->reduce(t_reduceDS_woCtErr);
  ws->import(*redDataCut);

  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setBins(45);
  // ws->var("ctau3D")->setRange(-lmin, lmax);
  ws->var("ctau3D")->setRange(-2, 4);
  ws->var("ctau3DErr")->setRange(0.0, 0.992);
  ws->var("ctau3DErr")->setBins(124);

  cout << " -*-*-*-*-*-*-*-*-*-*-*-*-*- getCtErrRange -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-  " << endl;
  cout << " *** DATA :: N events to fit (woCtErrRange) : " << redDataCut->sumEntries() << endl;

  RooDataSet *redDataSB = (RooDataSet *)redDataCut->reduce("mass < 2.9 || mass > 3.3");
  RooDataSet *redDataSIG = (RooDataSet *)redDataCut->reduce("mass > 2.9 && mass < 3.3");

  //// *** RooDataHist
  RooDataHist *tbinDataCtErrSB = new RooDataHist("tbinDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSB);
  RooDataHist *tbinDataCtErrSIG = new RooDataHist("tbinDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIG);

  //// *** mass fit to get coefExp of coefPol for scaleF
  // defineMassBkg(ws);
  // defineMassSig(ws);
  struct PARAM
  {
    double fracG1;
    double fracG1Err;
    double meanSig;
    double meanSigErr;
    double sigmaSig1;
    double sigmaSig1Err;
    double sigmaSig2;
    double sigmaSig2Err;
    double alpha;
    double alphaErr;
    double enne;
    double enneErr;
  };
  double cutValue_merged;

  char funct[100];
  double initBkg = redDataSB->sumEntries() * 9.0 / 5.0;
  double initSig = redDataCut->sumEntries() - initBkg;
  sprintf(funct, "SUM::MassPDF(NSig[%f,1.0,50000.0]*G1CB1Sig,NBkg[%f,1.0,500000.0]*expBkg)", initSig, initBkg);
  ws->factory(funct);
  ws->pdf("MassPDF")->fitTo(*redDataCut, Extended(1), Minos(0), Save(1), SumW2Error(kTRUE), NumCPU(8));

  //// ****  scaleF to scale down ct err dist in 2.9-3.3 GeV/c2
  float bc;
  // if (!inOpt.mBkgFunct.compare("expBkg"))
  bc = ws->var("coefExp")->getVal();
  // else if (!inOpt.mBkgFunct.compare("polBkg"))
  //   bc = ws->var("coefPol")->getVal();
  float scaleF = (exp(2.9 * bc) - exp(3.3 * bc)) / (exp(2.6 * bc) - exp(2.9 * bc) + exp(3.3 * bc) - exp(3.5 * bc));

  //// *** (tbinSubtractedSIG) = (tbinDataCtErrSIG) - scaleF*(tbinDataCtErrSB)
  RooDataHist *tbinSubtractedSIG = new RooDataHist("tbinSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  RooDataHist *tbinScaledBKG = subtractSidebands(ws, tbinSubtractedSIG, tbinDataCtErrSIG, tbinDataCtErrSB, scaleF, "ctau3DErr");

  //// **** Check the minimum and maximum of the ctau error in signal and background regions
  TH1 *histDataCtErrSIG = tbinDataCtErrSIG->createHistogram("histDataCtErrSIG", *ws->var("ctau3DErr"));
  TH1 *histSubtractedSIG = tbinSubtractedSIG->createHistogram("histSubtractedSIG", *ws->var("ctau3DErr"));
  TH1 *histScaledBKG = tbinScaledBKG->createHistogram("histScaledBKG", *ws->var("ctau3DErr"));

  double minSig = 0.5, maxSig = 0.0, minBkg = 0.5, maxBkg = 0.0;
  double cutValue = 0.2;

  int maxBinSig = histSubtractedSIG->GetMaximumBin();
  int maxBinBkg = histScaledBKG->GetMaximumBin();

  minSig = histSubtractedSIG->GetBinLowEdge(maxBinSig);
  minBkg = histScaledBKG->GetBinLowEdge(maxBinBkg);
  // pick up lower bound next to other non-zero bins
  for (int xbins = maxBinSig; xbins > 0; xbins--)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      minSig = histSubtractedSIG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG low edge: " << histSubtractedSIG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinBkg; xbins > 0; xbins--)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      minBkg = histScaledBKG->GetBinLowEdge(xbins);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG low edge: " << histScaledBKG->GetBinLowEdge(xbins) << endl;
    }
    else
      break;
  }

  // pick up upper bound next to other non-zero bins
  maxSig = histSubtractedSIG->GetBinLowEdge(maxBinSig + 1);
  maxBkg = histScaledBKG->GetBinLowEdge(maxBinBkg + 1);
  for (int xbins = maxBinSig; xbins < histSubtractedSIG->GetNbinsX(); xbins++)
  {
    if (histSubtractedSIG->GetBinContent(xbins) > cutValue)
    {
      maxSig = histSubtractedSIG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: SIG binContent: " << histSubtractedSIG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: SIG upper edge: " << histSubtractedSIG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }
  for (int xbins = maxBinSig; xbins < histScaledBKG->GetNbinsX(); xbins++)
  {
    if (histScaledBKG->GetBinContent(xbins) > cutValue)
    {
      maxBkg = histScaledBKG->GetBinLowEdge(xbins + 1);
      //          cout << "getCtErrRange:: BKG binContent: " << histScaledBKG->GetBinContent(xbins) << endl;
      //          cout << "getCtErrRange:: BKG upper edge: " << histScaledBKG->GetBinLowEdge(xbins+1) << endl;
    }
    else
      break;
  }

  // choose the higher lower limit, lower upper limit
  double tmpMin = 0, tmpMax = 0;
  if (minSig > minBkg)
    tmpMin = minSig;
  else
    tmpMin = minBkg;
  if (maxSig < maxBkg)
    tmpMax = maxSig;
  else
    tmpMax = maxBkg;

  // round off lower limit -> allow more entries on the lower limits
  tmpMin = TMath::Floor(tmpMin * 1000);
  tmpMin = tmpMin / (double)1000.0;

  // round up upper limit -> allow more entries on the upper limits
  tmpMax = TMath::Ceil(tmpMax * 1000);
  tmpMax = tmpMax / (double)1000.0;

  char reduceDS[512];
  sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
  if (redDataTmp->sumEntries() < redDataCut->sumEntries() * 0.9)
  { // if ctau error range cuts off >10% events
    delete redDataTmp;
    sprintf(reduceDS, "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)redDataCut->reduce(reduceDS);
    tmpMin = minSig;
    tmpMax = maxSig;
  }
  if ((tmpMax - tmpMin) < 0.008)
  {
    cout << "getCtErrRange:: Maximum is less than minimum! Possibly there are few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  //// *** draw final ctau error plot
  TCanvas c0("ctau_err", "ctau_err", 500, 500);
  c0.Draw();
  c0.cd();
  c0.SetLogy(1);
  RooPlot *errframe2 = ws->var("ctau3DErr")->frame();
  // tbinDataCtErrSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kRed),LineColor(kRed));
  // tbinDataCtErrSB->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kGreen+2),LineColor(kGreen+2),MarkerStyle(24));
  // tbinScaledBKG->plotOn(errframe2,DataError(RooAbsData::SumW2),MarkerColor(kBlue),MarkerStyle(24),LineColor(kBlue));
  // tbinSubtractedSIG->plotOn(errframe2,DataError(RooAbsData::SumW2),LineColor(kWhite));
  const double max = errframe2->GetMaximum() * 1.3;
  errframe2->SetMaximum(max);
  errframe2->SetMinimum(0.2);
  errframe2->Draw();
  histDataCtErrSIG->SetMarkerColor(kRed);
  histDataCtErrSIG->SetLineColor(kWhite);
  histDataCtErrSIG->SetMarkerStyle(24);
  histDataCtErrSIG->GetXaxis()->CenterTitle(1);
  histDataCtErrSIG->GetYaxis()->CenterTitle(1);
  histDataCtErrSIG->Draw("pe");
  histScaledBKG->SetMarkerColor(kBlue);
  histScaledBKG->SetLineColor(kWhite);
  histScaledBKG->SetMarkerStyle(24);
  histScaledBKG->Draw("pe same");
  histSubtractedSIG->SetLineColor(kWhite);
  histSubtractedSIG->Draw("pe same");

  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(32);
  t->SetTextSize(0.035);

  // t->DrawLatex(0.92, 0.84, opt.rapString);
  // t->DrawLatex(0.92, 0.78, opt.ptString);
  // if (opt.EventActivity == 1)
  //   t->DrawLatex(0.92, 0.72, opt.ntrkString);
  // else if (opt.EventActivity == 2)
  //   t->DrawLatex(0.92, 0.72, opt.etString);

  char comment[200];
  sprintf(comment, "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
  t->SetTextSize(0.04);
  t->SetTextColor(kRed);
  t->DrawLatex(0.92, 0.6, comment);
  t->SetTextColor(kBlack);

  TLegend legsb(0.6, 0.19, 0.9, 0.35, NULL, "brNDC");
  legsb.SetFillStyle(0);
  legsb.SetBorderSize(0);
  legsb.SetShadowColor(0);
  legsb.SetMargin(0.2);
  legsb.SetTextFont(42);
  legsb.SetTextSize(0.035);
  legsb.AddEntry(histDataCtErrSIG, "sig cands", "p");
  legsb.AddEntry(histScaledBKG, "scaled bkg", "p");
  legsb.AddEntry(histSubtractedSIG, "sig (= cands - bkg)", "p");
  legsb.Draw("same");

  string titlestr =  "_CtErrGetRange_Log.pdf";
  c0.SaveAs(titlestr.c_str());

  *errmin = tmpMin;
  *errmax = tmpMax;
  //  cout << "getCtErrRange:: " << t_reduceDS_woCtErr << " " << lmin << " " << lmax << " " << *errmin << " " << *errmax << endl;

  delete ws;
  delete redDataCut;
  delete redDataTmp;
  // delete binData;
  // delete binDataCtErr;
  // delete binDataSB;
  delete tbinDataCtErrSB;
  delete tbinDataCtErrSIG;
  delete tbinSubtractedSIG;
  delete tbinScaledBKG;
  delete t;
}