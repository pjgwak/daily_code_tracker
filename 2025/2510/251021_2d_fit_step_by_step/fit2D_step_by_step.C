#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include <TROOT.h>
#include <RooCrystalBall.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <TLine.h>
#include <TROOT.h> // gROOT
#include <TSystem.h> // gSystem
#include <RooChebychev.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
// #include "RooHistPdfConv.h"
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooKeysPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooAddModel.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>
#include <TStopwatch.h>
#include <RooFormulaVar.h>
#include <RooDataHist.h>

using namespace std;
using namespace RooFit;

TGraph* drawPdf1D_marginal(RooAbsPdf* pdf,
                           RooRealVar* obs, RooRealVar* intObs1, RooRealVar* intObs2,
                           int npts=600)
{
  intObs1->setRange("i1", intObs1->getMin(), intObs1->getMax());
  intObs2->setRange("i2", intObs2->getMin(), intObs2->getMax());
  std::unique_ptr<RooAbsReal> f(
    pdf->createIntegral(RooArgSet(*intObs1,*intObs2),
                        RooArgSet(*obs,*intObs1,*intObs2),
                        "i1,i2")
  );

  double xmin = obs->getMin(), xmax = obs->getMax();
  auto* g = new TGraph(npts);
  for (int i=0;i<npts;++i){
    double x = xmin + (xmax-xmin)*i/(npts-1.0);
    obs->setVal(x);
    g->SetPoint(i, x, f->getVal());
  }

  // 정규화
  double area=0;
  for(int i=0;i<npts-1;++i){
    double x1,y1,x2,y2; g->GetPoint(i,x1,y1); g->GetPoint(i+1,x2,y2);
    area += 0.5*(y1+y2)*(x2-x1);
  }
  if(area>0) for(int i=0;i<npts;++i){
    double x,y; g->GetPoint(i,x,y); g->SetPoint(i,x,y/area);
  }
  return g;
}

// RooDataHist *subtractSidebands(RooWorkspace *ws, RooDataHist *binSubtrSIG, RooDataHist *binSIG, RooDataHist *binSB, float scalefactor, string varName = "ctau3DErr")
// {

//   if (binSIG->numEntries() != binSB->numEntries())
//   {
//     cout << "ERROR subtractSidebands : different binning!" << endl;
//     return 0;
//   }
//   RooDataHist *binScaleBKG = new RooDataHist("binScaleBKG", "scaled SB", RooArgSet(*(ws->var(varName.c_str()))));

//   //// **** bin-by-bin scaling
//   const RooArgSet *argSIG;
//   const RooArgSet *argSB;
//   for (Int_t i = 0; i < binSIG->numEntries(); i++)
//   {
//     argSIG = binSIG->get(i);
//     argSB = binSB->get(i);
//     RooRealVar *thisVar = (RooRealVar *)argSIG->find(varName.c_str());
//     ws->var(varName.c_str())->setVal(thisVar->getVal());
//     //// *** set minimum as 0.1 to prevent empty PDF
//     float wBkg = binSB->weight(*argSB, 0, false);
//     if (wBkg <= 0.1)
//       wBkg = 0.1;
//     binScaleBKG->add(RooArgSet(*(ws->var(varName.c_str()))), wBkg);
//     float newWeight = binSIG->weight(*argSIG, 0, false) - scalefactor * binSB->weight(*argSB, 0, false);
//     if (newWeight <= 0.1)
//       newWeight = 0.1;
//     binSubtrSIG->add(RooArgSet(*(ws->var(varName.c_str()))), newWeight);
//   }
//   return binScaleBKG;
// }

// void drawCtauResolPlots(RooWorkspace *ws, bool fitMC, RooPlot *tframePR)
// {
//   TLatex *t = new TLatex();
//   t->SetNDC();
//   t->SetTextAlign(22);
//   char reduceDS[1000];
//   string titlestr;

//   t->SetTextSize(0.035);
//   TCanvas c00;
//   c00.cd();
//   tframePR->Draw();
//   sprintf(reduceDS, "#sigma(G_{N}): %.2f", ws->var("sigmaPRResN")->getVal());
//   t->DrawLatex(0.55, 0.31, reduceDS);
//   sprintf(reduceDS, "#sigma(G_{W}): %.2f", ws->function("sigmaPRResW")->getVal());
//   t->DrawLatex(0.55, 0.26, reduceDS);
//   sprintf(reduceDS, "frac(G_{W}): %.2f", ws->var("fracRes")->getVal());
//   t->DrawLatex(0.55, 0.21, reduceDS);
//   tframePR->GetXaxis()->CenterTitle(1);
//   tframePR->GetYaxis()->CenterTitle(1);
//   if (fitMC)
//     titlestr = "_CtPRResMC_Lin.pdf";
//   else
//     titlestr = "_CtPRResData_Lin.pdf";
//   c00.SaveAs(titlestr.c_str());

//   c00.SetLogy(1);
//   tframePR->GetXaxis()->CenterTitle(1);
//   tframePR->GetYaxis()->CenterTitle(1);
//   double prmax = tframePR->GetMaximum();
//   tframePR->SetMinimum(0.5);
//   tframePR->SetMaximum(prmax * 5);
//   if (fitMC)
//     titlestr = "_CtPRResMC_Log.pdf";
//   else
//     titlestr = "_CtPRResData_Log.pdf";
//   c00.SaveAs(titlestr.c_str());
// }

// void drawFinalMass(RooWorkspace *ws, RooDataHist* redDataCut, float NSigNP_fin, float NBkg_fin, RooFitResult *fitMass, double* UnNormChi2_mass_t, int* nFitParam_mass_t, int* nFullBinsPull_mass_t, int* Dof_mass_t,double* Chi2_mass_t) {
//   //// **** Temporary variables for plotting
//   TLatex *t = new TLatex();  t->SetNDC();  t->SetTextAlign(32);
//   TLatex *ty = new TLatex();  ty->SetNDC();  ty->SetTextAlign(12);
//   char reduceDS[1000];
//   string titlestr;
  
//   // RooBinning rb(ws->var("mass")->getBinning().numBins(), ws->var("mass")->getBinning().array()); // 14-009
//   RooBinning rb(300, 2.6, 3.5);
//   RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
//   RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);

//   //// *** Mass plot
//   RooPlot *mframe = ws->var("mass")->frame();
//   redDataCut->plotOn(mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1),Binning(rb));
//   double avgBinWidth = rb.averageBinWidth();
//   mframe->GetYaxis()->SetTitle(Form("Counts / (%.0f MeV/c^{2 })",avgBinWidth*1000));
//   mframe->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
// //  mframe->GetXaxis()->SetTitleSize(0.048*1.2); //PAPER
// //  mframe->GetYaxis()->SetTitleSize(0.048*1.2); //PAPER
//   mframe->GetXaxis()->CenterTitle(1);
//   mframe->GetYaxis()->CenterTitle(1);
//   const double max = mframe->GetMaximum() * 1.3;
//   mframe->SetMaximum(max);
//   mframe->SetMinimum(0);

//     //// **** Fill color
//     // G1CB1Sig, expBkg
//     //ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kBlack),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//     ws->pdf("totPDF_PEE")->plotOn(mframe,DrawOption("F"),FillColor(kGray+2),FillStyle(3354),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent), NumCPU(32));
//     RooAddPdf tmpPDF("tmpPDF","tmpPDF",RooArgList(*(ws->pdf("G1CB1Sig")),*(ws->pdf("expBkg"))),RooArgList(tmpVar1,tmpVar2));
//     //tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kWhite),FillStyle(1001),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
//     tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kWhite), FillStyle(1001), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
//     //tmpPDF.plotOn(mframe,LineColor(kRed),DrawOption("F"),FillColor(kRed),FillStyle(3444),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
//     tmpPDF.plotOn(mframe, LineColor(kPink - 6), DrawOption("F"), FillColor(kRed - 7), FillStyle(3345), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
//     gStyle->SetHatchesLineWidth(2);
//     //ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),DrawOption("F"),FillColor(kAzure-9),FillStyle(1001),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//     ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), DrawOption("F"), FillColor(kBlue - 10), FillStyle(1001), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
//     //// **** Line color
//     //ws->pdf("totPDF_PEE")->plotOn(mframe,Components("expBkg"),LineColor(kBlue),LineStyle(7),LineWidth(5),Normalization(redDataCut->sumEntries(),RooAbsReal::NumEvent));
//     ws->pdf("totPDF_PEE")->plotOn(mframe, Components("expBkg"), LineColor(kBlue - 2), LineStyle(7), LineWidth(5), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
//     //tmpPDF.plotOn(mframe,LineColor(kRed),LineStyle(9),LineWidth(5),Normalization(NSigNP_fin+NBkg_fin,RooAbsReal::NumEvent));
//     tmpPDF.plotOn(mframe, LineColor(kPink - 6), LineStyle(11), LineWidth(5), Normalization(NSigNP_fin + NBkg_fin, RooAbsReal::NumEvent), NumCPU(32));
//     ws->pdf("totPDF_PEE")->plotOn(mframe, LineColor(kBlack), LineWidth(2), Normalization(redDataCut->sumEntries(), RooAbsReal::NumEvent), NumCPU(32));
//     redDataCut->plotOn(mframe, DataError(RooAbsData::SumW2), XErrorSize(0), MarkerSize(1), Binning(rb));

//     TH1 *hdata = redDataCut->createHistogram("hdata", *ws->var("mass"), Binning(rb));
//     //// *** Calculate chi2/nDof for mass fitting
//     int nBins = hdata->GetNbinsX();
//     RooHist *hpullm;
//     hpullm = mframe->pullHist();
//     hpullm->SetName("hpullM");
//     double Chi2 = 0;
//     int nFullBinsPull = 0;
//     double *ypull = hpullm->GetY();
//     for (unsigned int i = 0; i < nBins; i++)
//     {
//       if (hdata->GetBinContent(i + 1) == 0)
//         continue;
//       nFullBinsPull++;
//       Chi2 = Chi2 + pow(ypull[i], 2);
//     }
//   double UnNormChi2 = Chi2;
//   *UnNormChi2_mass_t = Chi2;
//   int nFitParam = fitMass->floatParsFinal().getSize();
//   int Dof = nFullBinsPull - nFitParam;
//   Chi2 /= (nFullBinsPull - nFitParam);
//   *nFitParam_mass_t = nFitParam;
//   *nFullBinsPull_mass_t = nFullBinsPull;
//   *Dof_mass_t = Dof;
//   *Chi2_mass_t = Chi2;

//   // PAPER
//   TCanvas* c1wop = new TCanvas("c1wop","The Canvas",200,10,600,600);
//   c1wop->cd(); c1wop->Draw(); 
  
//   mframe->SetTitleOffset(1.47,"Y");
//   mframe->Draw();
//   //// **** different lumiTextOffset for massfit_wopull
//   // lumiTextOffset   = 0.20;
//   // CMS_lumi(c1wop, opt.isPA, iPosPaper);
//   t->SetTextSize(0.035);
// //  t->DrawLatex(0.91,0.85,opt.rapString);
// //  t->DrawLatex(0.91,0.78,opt.ptString);
//   ty->SetTextSize(0.035); // PAPER
//   // ty->DrawLatex(0.20,0.86,opt.rapString); //PAPER
//   // ty->DrawLatex(0.20,0.80,opt.ptString); //PAPER
//   ty->SetTextSize(0.040);
//   sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NSig")->getVal(),ws->var("NSig")->getError());
// //  ty->DrawLatex(0.20,0.85,reduceDS);
//   // sprintf(reduceDS,"#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
// //  ty->DrawLatex(0.20,0.79,reduceDS);

//   //TLegend * legpaper = new TLegend(0.59,0.53,0.90,0.72,NULL,"brNDC");
//   //TLegend * legpaper = new TLegend(0.17,0.53,0.63,0.72,NULL,"brNDC");
//   TLegend * legpaper = new TLegend(0.17,0.55,0.57,0.75,NULL,"brNDC");
//   legpaper->SetFillStyle(0); legpaper->SetBorderSize(0); legpaper->SetShadowColor(0);
//   legpaper->SetTextSize(0.035); legpaper->SetTextFont(42);
//   legpaper->SetMargin(0.2);
//   // legpaper->AddEntry(gfake1,"Data","p");
//   // legpaper->AddEntry(&hfake21,"Total fit","lf");
//   // legpaper->AddEntry(&hfake31,"Bkg + nonprompt","lf"); 
//   // legpaper->AddEntry(&hfake11,"Background","lf");
//   legpaper->Draw("same");

//   titlestr = "_massfit_wopull.pdf";
//   c1wop->SaveAs(titlestr.c_str());
// //  titlestr = opt.dirName + "_rap" + opt.yrange + "_pT" + opt.ptrange + "_ntrk" + inOpt.ntrrange + "_ET" + inOpt.etrange + "_massfit_wopull.root";
// //  c1wop.SaveAs(titlestr.c_str());
//   /////////////////////////////////////////////////////////////////////////////////  

//   TCanvas c1("c1","The mass Canvas",200,10,600,750);
//   c1.cd();
//   TPad *padm1 = new TPad("padm1","This is pad1",0.0,0.3,1.0,1.0);
//   padm1->SetLeftMargin(0.14);
//   padm1->SetRightMargin(0.03);
//   padm1->SetTopMargin(0.075);
//   padm1->SetBottomMargin(0);  padm1->Draw();
//   TPad *padm2 = new TPad("padm2","This is pad2",0.00,0.00,1.0,0.3);
//   padm2->SetLeftMargin(0.14);
//   padm2->SetRightMargin(0.03);
//   padm2->SetTopMargin(0);  
//   padm2->SetBottomMargin(0.30);  padm2->Draw();

//   padm1->cd();  mframe->Draw();
//   // lumiTextOffset   = 0.45;
//   // CMS_lumi(&c1, opt.isPA, iPos);
//   // t->SetTextSize(0.035);
//   // t->DrawLatex(0.91,0.90,opt.rapString);
//   // t->DrawLatex(0.91,0.85,opt.ptString);
//   // if (opt.EventActivity ==1) t->DrawLatex(0.91,0.80,opt.ntrkString);
//   // else if (opt.EventActivity ==2) t->DrawLatex(0.91,0.80,opt.etString);
  
//   ty->SetTextSize(0.040);
//   // sprintf(reduceDS,"N_{J/#psi} = %0.0f #pm %0.0f",ws->var("NSig")->getVal(),ws->var("NSig")->getError());
//   ty->DrawLatex(0.20,0.89,reduceDS);
//   // sprintf(reduceDS,"#sigma = %0.0f #pm %0.0f MeV/c^{2}", opt.PcombinedWidth, opt.PcombinedWidthErr);
//   ty->DrawLatex(0.20,0.84,reduceDS);
  
//   TLegend * leg11 = new TLegend(0.18,0.51,0.54,0.72,NULL,"brNDC");
//   leg11->SetFillStyle(0); leg11->SetBorderSize(0); leg11->SetShadowColor(0);
//   leg11->SetTextSize(0.035); leg11->SetTextFont(42);
//   leg11->SetMargin(0.2);
//   // leg11->AddEntry(gfake1,"Data","p");
//   // leg11->AddEntry(&hfake21,"Total fit","lf");
//   // leg11->AddEntry(&hfake31,"Bkg + nonprompt","lf"); 
//   // leg11->AddEntry(&hfake11,"Background","lf");
//   leg11->Draw("same");
//   c1.Update();

//   //// **** pull
//   RooPlot* mframepull =  ws->var("mass")->frame(Title("Pull")) ;
//   mframepull->GetYaxis()->SetTitle("Pull");
//   mframepull->GetYaxis()->CenterTitle(1);
//   mframepull->SetLabelSize(0.04*2.5,"XYZ");
//   mframepull->SetTitleSize(0.048*2.5,"XYZ");
//   mframepull->SetTitleOffset(0.47,"Y");
//   mframepull->addPlotable(hpullm,"PX") ;
//   double mframemax = 0;
//   if (mframepull->GetMinimum()*-1 > mframepull->GetMaximum()) mframemax = mframepull->GetMinimum()*-1;
//   else mframemax = mframepull->GetMaximum();
//   mframepull->SetMaximum(mframemax); 
//   mframepull->SetMinimum(-1*mframemax); 
//   mframepull->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
//   mframepull->GetXaxis()->CenterTitle(1);

//   padm2->cd(); mframepull->Draw();
//   TLine* line1 = new TLine(2.6,0,3.5,0.); line1->SetLineStyle(7); line1->Draw();

//   TLatex *t2 = new TLatex();
//   t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.035*3);
//   sprintf(reduceDS,"#chi^{2}/dof = %.1f/%d",UnNormChi2,Dof);
//   t2->DrawLatex(0.78,0.86,reduceDS);
//   c1.Update();

//   titlestr = "_massfit.pdf";
//   c1.SaveAs(titlestr.c_str());
//   /////////////////////////////////////////////////////////////////////////////////  
//   mframe->SetMinimum(0.5);
//   mframe->SetMaximum(max*50);
//   padm1->cd();
//   padm1->SetLogy(1);
//   titlestr = "_massfit_Log.pdf";
//   c1.SaveAs(titlestr.c_str());
//   /////////////////////////////////////////////////////////////////////////////////  
// }

// ///////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////

// void drawFinalCtau(RooWorkspace *ws, RooDataHist *redDataCut, RooDataHist* binDataCtErr, float NSigNP_fin, float NBkg_fin, float Bfrac_fin, float ErrBfrac_fin, RooFitResult *fit2D, float lmin, float lmax, double* UnNormChi2_time_t, int* nFitParam_time_t, int* nFullBinsPull_time_t, int* Dof_time_t,double* Chi2_time_t) {
//   char reduceDS[1000];
//   string titlestr;
  
//   // RooBinning rb(ws->var("ctau3D")->getBinning().numBins(), ws->var("ctau3D")->getBinning().array()); // 14-009
//   RooBinning rb(300, -0.5, 2.0);

//   RooRealVar tmpVar1("tmpVar1","tmpVar1",NSigNP_fin);
//   RooRealVar tmpVar2("tmpVar2","tmpVar2",NBkg_fin);
  
//   RooPlot *tframe = ws->var("ctau3D")->frame();
//   tframe->SetTitleOffset(1.47,"Y");
//   double avgBinWidth = rb.averageBinWidth();
//   //tframe->GetYaxis()->SetTitle(Form("Counts / %.2f (mm)",avgBinWidth));
//   tframe->GetYaxis()->SetTitle(Form("Counts / (%.0f #mum)",avgBinWidth*1000));
//   tframe->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
//   tframe->GetXaxis()->CenterTitle(1);
//   tframe->GetYaxis()->CenterTitle(1);

//   //// **** Ctau total distributions
//   RooHist *hpulltot;
//   redDataCut->plotOn(tframe,DataError(RooAbsData::SumW2),Binning(rb),MarkerSize(1));

//     ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent));
//     hpulltot = tframe->pullHist(); hpulltot->SetName("hpulltot");
//     ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtBkg"),LineColor(kBlue-2),LineWidth(5),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(7));
//     ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtNP"),LineColor(kPink-6),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(12));
//     ws->pdf("totPDF_PEE")->plotOn(tframe,Components("MassCtPR"),LineColor(kGreen+3),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent),LineStyle(kDashDotted));
//     ws->pdf("totPDF_PEE")->plotOn(tframe,LineColor(kBlack),LineWidth(2),ProjWData(RooArgList(*(ws->var("ctau3DErr"))),*binDataCtErr,kTRUE),NumCPU(32),Normalization(1,RooAbsReal::NumEvent));

//   TH1 *hdatact = redDataCut->createHistogram("hdatact",*ws->var("ctau3D"),Binning(rb));
//   double chi2 = 0, unNormChi2 = 0;
//   int dof = 0;
//   double *ypulls = hpulltot->GetY();
//   unsigned int nBins = ws->var("ctau3D")->getBinning().numBins();
//   unsigned int nFullBins = 0;
//   for (unsigned int i = 0; i < nBins; i++) {
//     if (hdatact->GetBinContent(i+1) == 0) continue;
//     chi2 += ypulls[i]*ypulls[i];
//     nFullBins++;
//   }
//   unNormChi2 = chi2;
//   *UnNormChi2_time_t = chi2;
//   int nFitPar = fit2D->floatParsFinal().getSize();
//   dof = nFullBins - nFitPar;
//   chi2 /= (nFullBins - nFitPar);
//   *nFitParam_time_t = nFitPar;
//   *nFullBinsPull_time_t = nFullBins;
//   *Dof_time_t =dof;
//   *Chi2_time_t = chi2;
  
//   TCanvas* c2 = new TCanvas("c2","The Canvas",200,10,600,750);
//   c2->cd();
//   TPad *pad1 = new TPad("pad1","This is pad1",0.0,0.3,1.0,1.0);
//   pad1->SetLeftMargin(0.14);
//   pad1->SetRightMargin(0.03);
//   pad1->SetTopMargin(0.075);
//   pad1->SetBottomMargin(0);  pad1->Draw();
//   TPad *pad2 = new TPad("pad2","This is pad2",0.0,0.0,1.0,0.3);
//   pad2->SetLeftMargin(0.14);
//   pad2->SetRightMargin(0.03);
//   pad2->SetTopMargin(0);  
//   pad2->SetBottomMargin(0.30);  pad2->Draw();

//   pad1->cd(); tframe->Draw();
  
//   // lumiTextOffset   = 0.45;
//   // CMS_lumi(c2, opt.isPA, iPos);
//   TLatex *ty = new TLatex();  ty->SetNDC();  ty->SetTextAlign(32);
//   ty->SetTextSize(0.035);
//   // ty->DrawLatex(0.91,0.90,opt.rapString);
//   // ty->DrawLatex(0.91,0.85,opt.ptString);
//   // if (opt.EventActivity ==1) ty->DrawLatex(0.91,0.80,opt.ntrkString);
//   // else if (opt.EventActivity ==2) ty->DrawLatex(0.91,0.80,opt.etString);

//   //TLegend * leg = new TLegend(0.66,0.56,0.85,0.75,NULL,"brNDC");
//   TLegend * leg = new TLegend(0.66,0.56,0.99,0.75,NULL,"brNDC");
//   leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetShadowColor(0);
//   leg->SetTextSize(0.035); leg->SetTextFont(42);
//   leg->SetMargin(0.2);
//   // leg->AddEntry(gfake1,"Data","p");
//   // leg->AddEntry(&hfake21,"Total fit","l");
//   // leg->AddEntry(&hfake41,"Prompt","l"); 
//   // leg->AddEntry(&hfake311,"Nonprompt","l"); 
//   // leg->AddEntry(&hfake11,"Background","l");
//   leg->Draw("same"); 
 
//   // KYO : write down Bfrac on plot 
//   TLatex *tbfrac = new TLatex();
//   tbfrac->SetNDC(); tbfrac->SetTextAlign(12); tbfrac->SetTextSize(0.035);
//   sprintf(reduceDS,"B frac. = %.2f #pm %.2f",Bfrac_fin,ErrBfrac_fin);
//   tbfrac->DrawLatex(0.19,0.91,reduceDS);

//   RooPlot* tframepull =  ws->var("ctau3D")->frame(Title("Pull")) ;
//   tframepull->GetYaxis()->SetTitle("Pull");
//   tframepull->GetYaxis()->CenterTitle(1);
//   tframepull->SetLabelSize(0.04*2.5,"XYZ");
//   tframepull->SetTitleSize(0.048*2.5,"XYZ");
//   tframepull->SetTitleOffset(0.47,"Y");
//   tframepull->addPlotable(hpulltot,"PX") ;
//   double tframemax = 0;
//   if (tframepull->GetMinimum()*-1 > tframepull->GetMaximum()) tframemax = tframepull->GetMinimum()*-1;
//   else tframemax = tframepull->GetMaximum();
//   tframepull->SetMaximum(tframemax); 
//   tframepull->SetMinimum(-1*tframemax);
//   tframepull->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} (mm)");
//   tframepull->GetXaxis()->CenterTitle(1);

//   pad2->cd(); tframepull->Draw();
//   TLine* line1 = new TLine(-lmin,0,lmax,0.); line1->SetLineStyle(7); line1->Draw();

//   int nDOF = ws->var("ctau3D")->getBinning().numBins() - nFitPar;

//   TLatex *t2 = new TLatex();
//   t2->SetNDC(); t2->SetTextAlign(22); t2->SetTextSize(0.035*3);
//   sprintf(reduceDS,"#chi^{2}/dof = %.2f/%d",unNormChi2,dof);
//   t2->DrawLatex(0.78,0.90,reduceDS);
  
//   c2->Update();
//   titlestr = "_timefit_Lin.pdf";
//   c2->SaveAs(titlestr.c_str());
//   /////////////////////////////////////////////////////////////////////////
//   tframe->SetMaximum(tframe->GetMaximum()*9); 
//   tframe->SetMinimum(0.5); 
//   pad1->SetLogy(1);
//   titlestr = "_timefit_Log.pdf";
//   c2->SaveAs(titlestr.c_str());
//   /////////////////////////////////////////////////////////////////////////
// }

void fit2D_step_by_step()
{
  cout << "=== start fit2D_step_by_step() ===\n";
  TStopwatch time;
  time.Start();

  // time check helper
  auto output_current_time = []() {
    auto now = chrono::system_clock::now();
    time_t current_time = chrono::system_clock::to_time_t(now);
    tm* local_tm = localtime(&current_time);
    cout << "Current time: " << put_time(local_tm, "%c") << "\n";
  };

  // === set kinematics ===
  // --- kinematics ---
  float ptLow = 6.5, ptHigh = 7.5;
  float yLow = 0, yHigh = 2.4;
  float massLow = 2.6, massHigh = 3.5;
  int cLow = 0, cHigh = 180;

  // --- initial observable range ---
  double ctLow = -1, ctHigh = 4;
  float ctErrLow = 0, ctErrHigh = 0.99;

  // --- fit skip flag ---
  bool isSkipMcMass = true, isSkipRawMass = true;
  bool isSkipMass = true, isSkipRes = true;
  bool isSkipBkg = true, isSkipTrue = true;
  bool isSkipFinal = true;

  // --- mass bkg parameters ---
  int bkgMassOrder = 2; // must be same with the order or Chebychev (check codes defining mass bkg model)
  double sbL_lo = 2.6, sbL_hi = 2.9;
  double sbR_lo = 3.3, sbR_hi = 3.5;
  double sig_lo = 2.9, sig_hi = 3.3;


  // === silence ===
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  // RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  // RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  // RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // only print from WARNING to FATAL


  // === set output properies ===
  // --- make output folders ---
  TString userLabel = "";
  TString figDir = Form("figs%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(figDir, kTRUE);
  TString rootDir = Form("roots%s/pT%.1f_%.1f_y%.1f_%.1f", userLabel.Data(), ptLow, ptHigh, yLow, yHigh);
  gSystem->mkdir(rootDir, kTRUE);
  gSystem->mkdir(Form("logs%s", userLabel.Data()), kTRUE);


  // === set cosmetics - rootlogon ===
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");


  // === read inputs ===
  cout << "\n=== import inputs ===\n";

  // --- prompt mc ---
  string fileNamePrMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";
  TFile fInPRMC(fileNamePrMc.c_str());
  // cout << fileNamePrMc.c_str() << endl;
  RooDataSet *dataPRMC = (RooDataSet *)fInPRMC.Get("dataset");
  dataPRMC->SetName("dataPRMC");

  // --- nonprompt mc ---
  string fileNameNpMc = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_JPsi_GENONLY_NonPrompt_y0_2p4_230829.root";
  TFile fInNPMC(fileNameNpMc.c_str());
  // cout << fileNameNpMc.c_str() << endl;
  RooDataSet *dataNPMC = (RooDataSet *)fInNPMC.Get("dataset");
  dataNPMC->SetName("dataNPMC");

  // --- data ---
  string fileNameData = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";
  TFile fInData(fileNameData.c_str());
  // cout << fileNameData.c_str() << endl;
  RooDataSet *data = (RooDataSet *)fInData.Get("dataset");
  data->SetName("data");


  // === created workspace ===
  RooWorkspace *ws = new RooWorkspace("ws");


  // === set cuts 1 ===
  cout << "\n=== make cuts 1 ===\n";
  string reduceDS_woCtErr = Form("( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) && (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f && mass >= %.3f && mass < %.3f && ctau3D >= %.3f && ctau3D < %.3f) && (recoQQsign == 0) ",  ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctLow, ctHigh);


  // === reduce dataset 1 ===
  cout << "\n=== reduce dataset 1 ===\n";
  RooDataSet *redPRMC = (RooDataSet *)dataPRMC->reduce(reduceDS_woCtErr.c_str()); // smae cut with Data
  redPRMC->SetName("redPRMC");
  RooDataSet *redData_woCtErr = (RooDataSet *)data->reduce(reduceDS_woCtErr.c_str());
  redData_woCtErr->SetName("redData_woCtErr");
  
  ws->import(*redPRMC);


  // === mc mass fit ===
  cout << "\n=== mc mass fit ===\n";
  // --- set observables ---
  ws->var("mass")->setRange(2.6, 3.5);
  ws->var("mass")->setRange("massRange", 2.6, 3.5);

  // --- build mass signal model ---
  RooRealVar massMean("massMean", "", 3.096, 3.085, 3.102);
  RooRealVar massSigma1("massSigma1", "CB sigma", 0.030, 0.005, 0.100);
  RooRealVar massSigma12("massSigma12", "ratio of sigma 1 vs 2", 1.5, 1, 10);
  RooFormulaVar massSigma2("massSigma2", "@0*@1", {massSigma1, massSigma12});

  RooRealVar alphaL("alphaL", "", 1.5, 0, 100.0);
  RooRealVar nL("nL", "nL", 3.0, 1, 5.0);
  
  RooRealVar alphaLR("alphaLR", "", 1.0, 0.001, 10.0);
  RooRealVar nLR("nLR", "", 3.0, 1, 100.0);

  RooFormulaVar alphaR("alphaR", "@0*@1", {alphaL, alphaLR});
  RooFormulaVar nR("nR", "@0*@1", {nL, nLR});

  // RooFormulaVar alphaR("alphaR", "-@0", {alphaL});
  // RooFormulaVar nR("nR", "@0", {nL});
  // RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2)
  // RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);

  RooRealVar fCB1("fCB1", "frac(Gauss in total)", 0.20, 0.00, 1.00);
  RooRealVar fCB2("fCB2", "frac(CB2 in DCB)", 0.50, 0.00, 1.00);
  

  RooCBShape CB1("CB1", "left-tail CB", *ws->var("mass"), massMean, massSigma1, alphaL, nL);
  RooCBShape CB2("CB2", "right-tail CB", *ws->var("mass"), massMean, massSigma2, alphaR, nR);

  // Gauss
  RooRealVar massSigma1G("massSigma1G", "ratio of sigma1 vs sigmaG", 1.1, 1, 10);
  RooFormulaVar massSigmaG("massSigmaG", "@0*@1", {massSigma1, massSigma1G});
  RooGaussian massG("massG", "core Gauss", *ws->var("mass"), massMean, massSigmaG);
  
  RooCrystalBall DCB("DCB", "", *ws->var("mass"), massMean, massSigma1, alphaL, nL, alphaR, nR);

  // RooAddPdf massSigPdf("massSigPdf", "Signal DCB+G",
  //                   RooArgList(CB1, CB2, massG),
  //                   RooArgList(fCB1, fCB2),
  //                   true); // Recursive fracton
  RooAddPdf massSigPdf("massSigPdf", "",
                    RooArgList(DCB, massG),
                    RooArgList(fCB1)); // Recursive fracton

  RooRealVar nSigMassMc("nSigMassMc", "", 1e5, 1, 1e8);
  RooAddPdf mcMassFitModel("mcMassFitModel", "",
                  RooArgList(massSigPdf),
                  RooArgList(nSigMassMc));

  // --- fit ---
  // helper function
  auto syncVar = [](RooRealVar& v, const RooFitResult& fr) {
    auto findIn = [&](const RooArgList& lst) -> const RooRealVar* {
      if (auto* a = lst.find(v.GetName())) return dynamic_cast<const RooRealVar*>(a);
      return nullptr;
    };
    const RooRealVar* src = nullptr;
    if (!(src = findIn(fr.floatParsFinal())))
        src = findIn(fr.constPars());
    if (!src) return false;
    v.setVal(src->getVal());
    v.setError(src->getError());
    return true;
  };

  RooFitResult *fitMcMass;
  if (isSkipMcMass && !gSystem->AccessPathName(Form("%s/mc_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/mc_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitMcMass", tmp);
    if (tmp) fitMcMass = (RooFitResult*) tmp->Clone("fitMcMass");

    syncVar(alphaL,  *fitMcMass);
    syncVar(alphaLR,  *fitMcMass);
    syncVar(fCB1,  *fitMcMass);
    syncVar(massMean,  *fitMcMass);
    syncVar(massSigma1,  *fitMcMass);
    syncVar(massSigma1G,  *fitMcMass);
    syncVar(nL,  *fitMcMass);
    syncVar(nLR,  *fitMcMass);
    syncVar(nSigMassMc,  *fitMcMass);
  } 
  else
    fitMcMass = mcMassFitModel.fitTo(*redPRMC, Range("massRange"), Save(), Extended(), PrintLevel(0), NumCPU(32)); //, EvalBackend("legacy")
  

  // --- draw plots ---
  // set helper functions
  auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
  { return fr ? fr->findObject(n) : nullptr; };

  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *fr = ws->var("mass")->frame(Range(massLow, massHigh), Title("")); // , Bins(nBins)
    redPRMC->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    mcMassFitModel.plotOn(fr, Range("massRange"), Name("model"));
    mcMassFitModel.plotOn(fr, Range("massRange"), Components("DCB"), Name("DCB"), LineStyle(kDotted), LineColor(kRed));
    mcMassFitModel.plotOn(fr, Range("massRange"), Components("massG"), Name("G"), LineStyle(kDotted), LineColor(kGreen));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "DCB"))
        leg.AddEntry(o, "DCB", "pe");
      if (auto *o = findObj(fr, "G"))
        leg.AddEntry(o, "Gauss", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Prompt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitMcMass->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/mc_mass_status_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // --- parameter latex ---
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;

      // lambda function for printing
      auto print = [&](const char *title, const char *vname, RooAddPdf &model)
      {
        auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)  
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("N_{Sig}", "nSigMassMc", mcMassFitModel);
      print("mean", "massMean", mcMassFitModel);
      print("#alpha_{L}", "alphaL", mcMassFitModel);
      print("n_{L}", "nL", mcMassFitModel);
      print("#sigma_{1}", "massSigma1", mcMassFitModel);
      print("#sigma_{2/1}", "massSigma12", mcMassFitModel);
      print("#sigma_{G/1}", "massSigma1G", mcMassFitModel);
      print("f_{CB1}", "fCB1", mcMassFitModel);
      print("f_{CB2}", "fCB2", mcMassFitModel);
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ws->var("mass")->frame(Range(massLow, massHigh), Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(massLow, 0.0, massHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    // --- chi2/ndf ---
    if (fitMcMass)
    {
      int npar = fitMcMass->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/mc_mass_fit_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh)); 
  }

  // --- fix parameters ---
  // // // // massMean.setConstant();
  // massSigma1.setConstant();
  massSigma1G.setConstant();
  // alphaL.setConstant();
  nL.setConstant();
  alphaLR.setConstant();
  nLR.setConstant();
  // fCB1.setConstant();

  // --- save result for later---
  fitMcMass->Print("V");
  TFile f(Form("%s/mc_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitMcMass->Write("fitMcMass");
  f.Close();


  // === mass fit 1 (for sideband) ===
  cout << "\n=== mass raw fit ===\n";
  // --- build mass bkg model ---
  RooRealVar massSl1("massSl1", "massSl1", 0.05, -1.0, 1.0);
  RooRealVar massSl2("massSl2", "massSl2", 0.05, -1.0, 1.0);
  RooRealVar massSl3("massSl3", "massSl3", 0.01, -1.0, 1.0);
  RooChebychev massBkgPdf("massBkgPdf", "", *ws->var("mass"), RooArgList(massSl1, massSl2));

  // --- fit ---
  RooRealVar nSigMass("nSigMass", "signal yield", 5e5, 1e5, 1e6);
  RooRealVar nBkgMass("nBkgMass", "signal yield", 1e5, 9e4, 2e5);
  RooAddPdf massFitModel("massFitModel", "massSigPdf + Cheby3 (extended)",
                  RooArgList(massSigPdf, massBkgPdf),
                  RooArgList(nSigMass, nBkgMass));

  // --- fit ---
  RooFitResult *fitRawMass;
  if (isSkipRawMass && !gSystem->AccessPathName(Form("%s/raw_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/raw_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitRawMass", tmp);
    if (tmp) fitRawMass = (RooFitResult*) tmp->Clone("fitRawMass");

    syncVar(alphaL, *fitRawMass); 
    syncVar(alphaLR, *fitRawMass);
    syncVar(fCB1, *fitRawMass);
    syncVar(massMean, *fitRawMass);
    syncVar(massSigma1, *fitRawMass);
    syncVar(massSigma1G, *fitRawMass);
    syncVar(nL, *fitRawMass);
    syncVar(nLR, *fitRawMass);
    syncVar(nSigMass, *fitRawMass);
    syncVar(nBkgMass, *fitRawMass);
    syncVar(massSl1, *fitRawMass);
    syncVar(massSl2, *fitRawMass);
    syncVar(massSl3, *fitRawMass);
  } 
  else
    fitRawMass = massFitModel.fitTo(*redData_woCtErr, Offset(true), Range("massRange"), Save(), Extended(), PrintLevel(-1), Strategy(2), PrintEvalErrors(-1), RecoverFromUndefinedRegions(2), NumCPU(32)); //, EvalBackend("legacy")

  // --- draw plots ---
  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *fr = ws->var("mass")->frame(Range(massLow, massHigh), Title("")); // , Bins(nBins)
    redData_woCtErr->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    massFitModel.plotOn(fr, Range("massRange"), Name("model"));
    massFitModel.plotOn(fr, Range("massRange"), Components("DCB"), Name("DCB"), LineStyle(kDotted), LineColor(kRed));
    massFitModel.plotOn(fr, Range("massRange"), Components("massG"), Name("G"), LineStyle(kDotted), LineColor(kGreen));
    massFitModel.plotOn(fr, Range("massRange"), Components("massBkgPdf"), Name("bkg"), LineStyle(kDotted), LineColor(kAzure));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "DCB"))
        leg.AddEntry(o, "DCB", "pe");
      if (auto *o = findObj(fr, "G"))
        leg.AddEntry(o, "Gauss", "pe");
      if (auto *o = findObj(fr, "bkg"))
        leg.AddEntry(o, "Bkg", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitMcMass->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/raw_mass_status_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // --- parameter latex ---
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;

      // lambda function for printing
      auto print = [&](const char *title, const char *vname, RooAddPdf &model)
      {
        auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("N_{Sig}", "nSigMass", massFitModel);
      print("N_{Bkg}", "nBkgMass", massFitModel);

      print("mean", "massMean", massFitModel);
      print("#alpha_{L}", "alphaL", massFitModel);
      print("n_{L}", "nL", massFitModel);
      print("#alpha_{L/R}", "alphaLR", massFitModel);
      print("n_{L/R}", "nLR", massFitModel);

      print("#sigma_{1}", "massSigma1", massFitModel);
      // print("#sigma_{2/1}", "massSigma12", massFitModel);
      print("#sigma_{G/1}", "massSigma1G", massFitModel);
      
      print("f_{CB1}", "fCB1", massFitModel);
      // print("f_{CB2}", "fCB2", massFitModel);
      
      print("sl1", "massSl1", massFitModel);
      print("sl2", "massSl2", massFitModel);
      print("sl3", "massSl3", massFitModel);
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ws->var("mass")->frame(Range(massLow, massHigh), Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(massLow, 0.0, massHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    // --- chi2/ndf ---
    if (fitRawMass)
    {
      int npar = fitRawMass->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/raw_mass_fit_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh)); 
  }

  // --- save results ---
  fitRawMass->Print("V");
  TFile outMass(Form("%s/raw_mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitRawMass->Write("fitRawMass");
  outMass.Close();


  // === sideband subtraction ===
  cout << "\n=== sideband extraction ===\n";
  // --- set initial range of ctau3DErr ---
  ws->var("ctau3DErr")->setRange(0, 0.992);
  ws->var("ctau3DErr")->setBins(100);

  RooDataSet *redDataSBTemp = (RooDataSet *)redData_woCtErr->reduce("mass < 2.9 || mass > 3.3");
  RooDataSet *redDataSIGTmp = (RooDataSet *)redData_woCtErr->reduce("mass > 2.9 && mass < 3.3");

  // --- build temporal RooHistPdf ---
  // it seems "t" stands for temporal
  RooDataHist *tbinDataCtErrSB = new RooDataHist("tbinDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSBTemp);
  RooDataHist *tbinDataCtErrSIG = new RooDataHist("tbinDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*(ws->var("ctau3DErr"))), *redDataSIGTmp);

  // --- compute bkg scale factor ---
  // helper function
  auto computeScaleF = [&](int order_, double massMin_, double massMax_, double sbL_lo_, double sbL_hi_, double sig_lo_, double sig_hi_, double sbR_lo_, double sbR_hi_, const RooArgList& allPars_, const char* coeffPrefix_="sl") -> double
  {
    // reald parameters from allPars
    auto getVal = [&](const char *name) -> double {
      if (auto *obj = allPars_.find(name))
        if (auto *rr = dynamic_cast<RooAbsReal*>(obj)) return rr->getVal();
      return std::numeric_limits<double>::quiet_NaN();
    };

    // mass -> [-1, 1] transformation. RooFit internally does it for the Checbychev PDF
    auto toUnitX = [&](double m) -> double {
      return 2.0 * (m - massMin_) / (massMax_ - massMin_) - 1.0;
    };

    // T_n(x)
    auto chebT = [&](int n, double x)->double {
      if (n==0) return 1.0; // RooFit c0 = 1
      if (n==1) return x;
      double Tn_2 = 1.0, Tn_1 = x, Tn = 0.0;
      for (int k = 2; k <= n; ++k) {
        Tn = 2.0 * x * Tn_1 - Tn_2;
        Tn_2 = Tn_1;
        Tn_1 = Tn;
      }
      return Tn;
    };

    // RooChebychev
    auto chebVal_roofit = [&](double m, const std::vector<double>& c) -> double {
      const double x = toUnitX(m);
      double v = 1.0; //T0 = 1.0
      for (int i = 1; i <= (int)c.size(); ++i)
        v += c[i-1] * chebT(i ,x);
      return v;
    };

    // Simpson integration
    auto integrateSimpson = [&](auto&& f, double a, double b, int nSeg = 2048) {
      if (b <= a) return 0.0;
      if (nSeg % 2) ++nSeg;
      const double h = (b - a) / nSeg;
      double s = f(a) + f(b);
      for (int i = 1; i < nSeg; ++i) {
        const double x = a + i*h;
        s += f(x) * ((i%2)?4.0:2.0);
      }
      return s * (h/3.0);
    };

    // check order of chebychev
    if (order_ < 0 || order_ > 6) {
      ::Error("computeScaleF", "Unsupported Cheby order %d (allowed 0..6)", order_);
      return std::numeric_limits<double>::quiet_NaN();
    }

    // make coefficient vector
    std::vector<double> c; c.reserve(order_);

    // for (int i = 1; i <= order_; ++i) {
    //   double v = std::numeric_limits<double>::quiet_NaN();

    //   if (i == 1 && std::isfinite(massSl1)) v = massSl1;
    //   else if (i == 2 && std::isfinite(massSl2)) v = massSl2;
    //   else if (i == 3 && std::isfinite(massSl3)) v = massSl3;
    // }

  auto getLocalSl = [&](int i) -> double {
    if (std::strcmp(coeffPrefix_, "sl") != 0) 
      return std::numeric_limits<double>::quiet_NaN();

    // massSl1/2/3 이 RooRealVar (레퍼런스)라고 가정
    if (i == 1) return massSl1.getVal();
    if (i == 2) return massSl2.getVal();
    if (i == 3) return massSl3.getVal();
    return std::numeric_limits<double>::quiet_NaN();
  };

  for (int i = 1; i <= order_; ++i) {
    double v = getLocalSl(i);  // 1차~3차는 로컬 우선
    if (!std::isfinite(v)) {
      char nm[32]; 
      std::snprintf(nm, sizeof(nm), "%s%d", coeffPrefix_, i);
      v = getVal(nm);           // 로컬에 없거나 NaN이면 allPars에서 sl1, sl2...
    }
    if (!std::isfinite(v)) {
      char nm[32]; 
      std::snprintf(nm, sizeof(nm), "%s%d", coeffPrefix_, i);
      ::Error("computeScaleF_cheby", "Coefficient '%s' not found or NaN.", nm);
      return std::numeric_limits<double>::quiet_NaN();
    }
    c.push_back(v);
  }

    // compute integration and scaleF
    auto f = [&](double m){ return chebVal_roofit(m, c); };
    const double I_sig = integrateSimpson(f, sig_lo_, sig_hi_);
    const double I_sb = integrateSimpson(f, sbL_lo_, sbL_hi_) + integrateSimpson(f, sbR_lo_, sbR_hi_);

    if (std::abs(I_sb) < 1e-300) {
      ::Error("computeScaleF", "Sideband integral is zero. Check range.");
      return std::numeric_limits<double>::quiet_NaN();
    }

    return I_sig / I_sb;
  };

  RooArgList allPars_local;
  allPars_local.add(massSl1);
  allPars_local.add(massSl2);

  double scaleF = computeScaleF(bkgMassOrder, massLow, massHigh, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi, allPars_local, "massSl");
  cout << "scaleF: " << scaleF << "\n";

  // --- subtract sigHist - bkgHist*scaleF ---
  // (tbinSubtractedSIG) = (tbinDataCtErrSIG) - scaleF*(tbinDataCtErrSB)
  // helper function
  auto subtractSidebands = [](RooDataHist* outSub, RooDataHist *outScal, const RooDataHist *inSIG, const RooDataHist *inSB, RooRealVar *obsVar, const std::string& name, double scaleF, double floorW = 1e-1)->void {
    if (!outSub||!outScal||!inSIG||!inSB||!obsVar) {
      ::Error("subtractSidebands", "Null input");
      return;
    }
    const int n = inSIG->numEntries();
    if (n!=inSB->numEntries()) {
      ::Error("subtractSidebands", "Different binning");
      return;
    }

    TH1 *hSIG = inSIG->createHistogram("hSIG_tmp", *obsVar);
    TH1 *hSB = inSB->createHistogram("hSB_tmp", *obsVar);
    if (!hSIG||!hSB) {
      ::Error("subtractSideband", "createHistogram failed");
      delete hSIG; delete hSB;
      return;
    }
    hSIG->Sumw2(kTRUE); hSB->Sumw2(kTRUE);
    if (hSIG->GetNbinsX()!=hSB->GetNbinsX() || hSIG->GetNbinsX()!=n) {
      ::Error("subtractSidebands","Histogram mismatch");
      delete hSIG; delete hSB; return;
    }

    for (int i=0; i<n; ++i) {
      const RooArgSet *aSIG = inSIG->get(i);
      const RooArgSet *aSB = inSB->get(i);
      if (!aSIG||!aSB||!aSIG->find(name.c_str())||!aSB->find(name.c_str())) {
         ::Error("subtractSidebands","Var missing");
        delete hSIG; delete hSB; return;
      }

      const int j = i+1; // TH1D start from 1
      const double wSIG = hSIG->GetBinContent(j), eSIG = hSIG->GetBinError(j);
      const double wSB  = hSB ->GetBinContent(j), eSB = hSB->GetBinError(j);

      const double wSBs = scaleF * wSB;
      const double eSBs = std::fabs(scaleF) * eSB;
      const double wSUB = wSIG - wSBs;
      const double eSUB = std::hypot(eSIG, eSBs); // compute sqrt(x^2+y^2)

      outScal->set(*aSIG, (wSBs<=floorW?floorW:wSBs), eSBs);
      outSub->set(*aSIG, (wSUB<=floorW?floorW:wSUB), eSUB);
    }

    delete hSIG; delete hSB;
  };

  RooDataHist *tbinSubtractedSIG = new RooDataHist("tbinSubtractedSIG", "Subtracted data", RooArgSet(*(ws->var("ctau3DErr"))));
  RooDataHist *tbinScaledBKG = new RooDataHist("tbinScaledBKG", "Bkg data", RooArgSet(*(ws->var("ctau3DErr"))));

  subtractSidebands(tbinSubtractedSIG, tbinScaledBKG, tbinDataCtErrSIG, tbinDataCtErrSB, ws->var("ctau3DErr"), "ctau3DErr", scaleF);

  // --- find proper min and max ---
  // Just copied and pasted codes of HIN14-009
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

  // --- final ctau3DErr cutting ---
  char reduceDS_tmp[1000];
  sprintf(reduceDS_tmp, "ctau3DErr > %.3f && ctau3DErr < %.3f", tmpMin, tmpMax);
  RooDataSet *redDataTmp = (RooDataSet *)redData_woCtErr->reduce(reduceDS_tmp);
  if (redDataTmp->sumEntries() < redData_woCtErr->sumEntries() * 0.9)
  { // if ctau error range cuts off >10% events
    delete redDataTmp;
    sprintf(reduceDS_tmp, "ctau3DErr > %.3f && ctau3DErr < %.3f", minSig, maxSig);
    redDataTmp = (RooDataSet *)redData_woCtErr->reduce(reduceDS_tmp);
    tmpMin = minSig;
    tmpMax = maxSig;
  }
  if ((tmpMax - tmpMin) < 0.008)
  {
    cout << "getCtErrRange:: Maximum is less than minimum! Possibly there are few events in this bin.\n";
    tmpMax = tmpMin + 0.008;
  }

  // --- set new ctau3DErr min and max ---
  ctErrLow = tmpMin;
  if (ctErrLow < 0.001) ctErrLow = 0.001;
  ctErrHigh = tmpMax;
  cout << "ctErrRange: [" << ctErrLow << ", " << ctErrHigh << "]\n";

  // --- deallocate memory ---
  delete redDataTmp;
  delete redDataSBTemp;
  delete redDataSIGTmp;
  delete tbinSubtractedSIG;
  delete tbinScaledBKG;
  delete tbinDataCtErrSB;
  delete tbinDataCtErrSIG;

  // --- draw ctau3DErr plots ---
  {
    TCanvas c0("ctau_err", "ctau_err", 800, 800);
    c0.cd();
    c0.SetLogy(1);

    RooPlot *errFrame = ws->var("ctau3DErr")->frame();
    errFrame->SetTitle("");
    errFrame->GetXaxis()->SetTitle("ctau3DErr [mm]");
    errFrame->GetXaxis()->CenterTitle();
    errFrame->GetYaxis()->SetTitle("Counts");
    double maxY = std::max({histDataCtErrSIG->GetMaximum(),
                            histScaledBKG->GetMaximum(),
                            histSubtractedSIG->GetMaximum()});
    errFrame->SetMaximum(maxY * 1.3);
    errFrame->SetMinimum(0.2);
    errFrame->Draw();

    histDataCtErrSIG->SetMarkerColor(kRed);
    histDataCtErrSIG->SetLineColor(kWhite);
    histDataCtErrSIG->SetMarkerStyle(24);
    histDataCtErrSIG->GetXaxis()->CenterTitle(1);
    histDataCtErrSIG->GetYaxis()->CenterTitle(1);
    histDataCtErrSIG->Draw("pe same");

    histScaledBKG->SetMarkerColor(kBlue);
    histScaledBKG->SetLineColor(kWhite);
    histScaledBKG->SetMarkerStyle(24);
    histScaledBKG->Draw("pe same");

    histSubtractedSIG->SetLineColor(kWhite);
    histSubtractedSIG->SetMarkerSize(0.3);
    histSubtractedSIG->Draw("pe same");

    TLatex t;
    t.SetNDC();
    t.SetTextAlign(32);
    t.SetTextSize(0.04);
    t.SetTextColor(kRed);
    char comment[200];
    std::snprintf(comment, sizeof(comment), "Range: %.3f-%.3f (mm)", tmpMin, tmpMax);
    t.DrawLatex(0.92, 0.60, comment);
    t.SetTextColor(kBlack);

    TLegend legsb(0.60, 0.39, 0.90, 0.55, nullptr, "brNDC");
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


    c0.SaveAs(Form("%s/ctau_err_raw_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // === set cuts 2 ===
  cout << "\n=== set cuts 2 ===\n";
  const string kMuonSel =
    "( ((abs(eta1) <= 1.2) && (pt1 >= 3.5))"
    " || ((abs(eta2) <= 1.2) && (pt2 >= 3.5))"
    " || ((abs(eta1) >  1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1))))"
    " || ((abs(eta2) >  1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2))))"
    " || ((abs(eta1) >  2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5))"
    " || ((abs(eta2) >  2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )";

  // Data: ctau3DErr 포함
  string reduceDS = Form(
    "%s && (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f"
    " && mass >= %.3f && mass < %.3f"
    " && ctau3D >= %.3f && ctau3D < %.3f"
    " && ctau3DErr >= %.3f && ctau3DErr < %.3f)"
    " && (recoQQsign == 0)",
    kMuonSel.c_str(),
    ptLow, ptHigh, yLow, yHigh, massLow, massHigh, ctLow, ctHigh, ctErrLow, ctErrHigh);

  // Nonprompt MC: true-ctau + ctau3DErr(브랜치가 있을 때만 사용)
  string reduceNpMc = Form(
    "%s && (pt >= %.3f && pt < %.3f && abs(y) >= %.3f && abs(y) < %.3f"
    " && mass >= %.3f && mass < %.3f"
    " && ctau3Dtrue >= %.5f && ctau3Dtrue < %.3f)",
    kMuonSel.c_str(),
    ptLow, ptHigh, yLow, yHigh, massLow, massHigh, 1e-5, ctHigh);
  
  
  // === redcue dataset2 ===
  cout << "\n=== reduce dataset 2 ===\n";

  // --- set observables with new ranges ---
  auto mass = new RooRealVar("mass", "", massLow, massHigh, "GeV/c^{2} [mm]");
  auto ctau3D = new RooRealVar("ctau3D", "", ctLow, ctHigh, "ctau3D [mm]");
  auto ctau3Dtrue = new RooRealVar("ctau3Dtrue", "", ctLow, ctHigh, "ctau3Dtrue [mm]");
  auto ctau3DErr = new RooRealVar("ctau3DErr", "", ctErrLow, ctErrHigh, "ctau3DErr [mm]");

  RooArgSet obsData(*mass, *ctau3D, *ctau3DErr);
  RooArgSet obsNpMc(*mass, *ctau3Dtrue);
  

  RooDataSet* redNPMC_tmp = (RooDataSet *)dataNPMC->reduce(reduceNpMc.c_str());
  RooDataSet* redData_tmp = (RooDataSet *)data->reduce(reduceDS.c_str());

  auto redNPMC = new RooDataSet("redNPMC", "", obsNpMc, Import(*redNPMC_tmp));
  auto redData = new RooDataSet("redData", "", obsData, Import(*redData_tmp));
  
  
  redNPMC->SetName("redNPMC");
  redData->SetName("redData");
  ws->import(*redNPMC);
  ws->import(*redData);

  // === set variable properties ===
  cout << "\n=== set variable properties ===\n";
  // --- set ranges 2 ---
  ws->var("mass")->setRange(2.6, 3.5);
  
  ws->var("ctau3Dtrue")->setRange(1e-5, ctHigh);
  ws->var("ctau3Dtrue")->setRange("trueRange", 1e-5, ctHigh);

  ws->var("ctau3D")->setRange(ctLow, ctHigh);
  ws->var("ctau3D")->setRange("ctRange", ctLow, ctHigh);
  // ws->var("ctau3D")->setRange("promptMCfit", -0.2, 0.2);

  ws->var("ctau3DRes")->setRange(-10, 10);
  ws->var("ctau3DRes")->setRange("ctResRange", -10, 10);

  ctau3DErr->setRange(ctErrLow, ctErrHigh);
  ctau3DErr->setRange("errRange", ctErrLow, ctErrHigh);
  ctau3DErr->setBins(100);

  // --- set labels --- 
  ws->var("mass")->SetTitle("m_{#mu#mu}");
  ws->var("ctau3D")->SetTitle("#font[12]{l}_{J/#psi}");
  
  // // -- set bins ---
  // // skip
  // ws->var("mass")->setBins(120); // default 45: bin/0.02 GeV
  // ws->var("ctau3D")->setBins(200);
  // ws->var("ctau3Dtrue")->setBins(50); // default
  // ws->var("ctau3DErr")->setBins(50);

  // === cut mass subrange dataset ===
  cout << "\n=== cut mass subrange dataset ===\n";
  RooDataSet *redDataSIG = (RooDataSet *)redData->reduce("mass > 2.9 && mass < 3.3");
  RooDataSet *redDataSB = (RooDataSet *)redData->reduce("mass<2.9 || mass>3.3");
  // RooDataSet *redDataSBL = (RooDataSet *)redData->reduce("mass<2.9");
  // RooDataSet *redDataSBR = (RooDataSet *)redData->reduce("mass>3.3");


  // === mass fit 2 ===
  cout << "\n=== mass fit 2 ===\n";
  // --- fit ---
  RooFitResult *fitMass;
  if (isSkipMass && !gSystem->AccessPathName(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitMass", tmp);
    if (tmp) fitMass = (RooFitResult*) tmp->Clone("fitMass");

    syncVar(alphaL, *fitMass); 
    syncVar(alphaLR, *fitMass);
    syncVar(fCB1, *fitMass);
    syncVar(massMean, *fitMass);
    syncVar(massSigma1, *fitMass);
    syncVar(massSigma1G, *fitMass);
    syncVar(nL, *fitMass);
    syncVar(nLR, *fitMass);
    syncVar(nSigMass, *fitMass);
    syncVar(nBkgMass, *fitMass);
    syncVar(massSl1, *fitMass);
    syncVar(massSl2, *fitMass);
    syncVar(massSl3, *fitMass);
  } 
  else
    fitMass = massFitModel.fitTo(*redData_woCtErr, Offset(true), Range("massRange"), Save(), Extended(), PrintLevel(-1), Strategy(2), PrintEvalErrors(-1), RecoverFromUndefinedRegions(2), NumCPU(32)); //, EvalBackend("legacy")

  // --- draw plots ---
  {
    TCanvas c("c_mass", "c_mass", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *fr = ws->var("mass")->frame(Range(massLow, massHigh), Title("")); // , Bins(nBins)
    redData_woCtErr->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    massFitModel.plotOn(fr, Range("massRange"), Name("model"));
    massFitModel.plotOn(fr, Range("massRange"), Components("DCB"), Name("DCB"), LineStyle(kDotted), LineColor(kRed));
    massFitModel.plotOn(fr, Range("massRange"), Components("massG"), Name("G"), LineStyle(kDotted), LineColor(kGreen));
    massFitModel.plotOn(fr, Range("massRange"), Components("massBkgPdf"), Name("bkg"), LineStyle(kDotted), LineColor(kAzure));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "DCB"))
        leg.AddEntry(o, "DCB", "pe");
      if (auto *o = findObj(fr, "G"))
        leg.AddEntry(o, "Gauss", "pe");
      if (auto *o = findObj(fr, "bkg"))
        leg.AddEntry(o, "Bkg", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitMcMass->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/raw_mass_status_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // --- parameter latex ---
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;

      // lambda function for printing
      auto print = [&](const char *title, const char *vname, RooAddPdf &model)
      {
        auto *v = dynamic_cast<RooRealVar *>(model.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("N_{Sig}", "nSigMass", massFitModel);
      print("N_{Bkg}", "nBkgMass", massFitModel);

      print("mean", "massMean", massFitModel);
      print("#alpha_{L}", "alphaL", massFitModel);
      print("n_{L}", "nL", massFitModel);
      print("#alpha_{L/R}", "alphaLR", massFitModel);
      print("n_{L/R}", "nLR", massFitModel);

      print("#sigma_{1}", "massSigma1", massFitModel);
      // print("#sigma_{2/1}", "massSigma12", massFitModel);
      print("#sigma_{G/1}", "massSigma1G", massFitModel);
      
      print("f_{CB1}", "fCB1", massFitModel);
      // print("f_{CB2}", "fCB2", massFitModel);
      
      print("sl1", "massSl1", massFitModel);
      print("sl2", "massSl2", massFitModel);
      print("sl3", "massSl3", massFitModel);
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ws->var("mass")->frame(Range(massLow, massHigh), Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("mass^{inv}_{#mu#mu} [GeV/c^{2}]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(massLow, 0.0, massHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    // --- chi2/ndf ---
    if (fitMass)
    {
      int npar = fitMass->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/mass_fit_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh)); 
  }

  // --- fix parameters ---
  alphaL.setConstant();
  fCB1.setConstant();
  massMean.setConstant();
  massSigma1.setConstant();
  massSl1.setConstant();
  massSl2.setConstant();
  nBkgMass.setConstant();
  nSigMass.setConstant();

  // --- save results ---
  fitMass->Print("V");
  TFile outMass2(Form("%s/mass_fit_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitMass->Write("fitMass");
  outMass2.Close();

  // === make RooHistPdf ===
  cout << "\n=== make RooHistPdf ===\n";
  // very similar with "sideband extraction" but we use the dataset applied ctau3DErr cut

  // --- build RooDataSet and RooHistPdf ---
  RooDataHist *binDataCtErrSB = new RooDataHist("binDataCtErrSB", "Data ct error distribution for bkg", RooArgSet(*ctau3DErr), *redDataSB);
  RooDataHist *binDataCtErrSIG = new RooDataHist("binDataCtErrSIG", "Data ct error distribution for sig", RooArgSet(*ctau3DErr), *redDataSIG);

  // --- compute bkg scale factor ---
  // RooArgList allPars_local; // use same variable
  scaleF = computeScaleF(bkgMassOrder, massLow, massHigh, sbL_lo, sbL_hi, sig_lo, sig_hi, sbR_lo, sbR_hi, allPars_local, "massSl");
  cout << "scaleF: " << scaleF << "\n";

  // --- subtract sigHist - bkgHist*scaleF ---
  // (tbinSubtractedSIG) = (binDataCtErrSIG) - scaleF*(binDataCtErrSB)
  RooDataHist *binSubtractedSIG = new RooDataHist("binSubtractedSIG", "Subtracted data", RooArgSet(*ctau3DErr));
  RooDataHist *binScaledBKG = new RooDataHist("binScaledBKG", "Bkg data", RooArgSet(*ctau3DErr));

  subtractSidebands(binSubtractedSIG, binScaledBKG, binDataCtErrSIG, binDataCtErrSB, ctau3DErr, "ctau3DErr", scaleF, 0.001);

  // --- buld RooHistPdf ---
  auto errSigPdf = new RooHistPdf("errSigPdf", "Signal hist",
                    *ctau3DErr, *binSubtractedSIG, 0);

  auto errBkgPdf = new RooHistPdf("errBkgPdf", "Background hist",
                    *ctau3DErr, *binScaledBKG, 0);
  ws->import(*errSigPdf);
  ws->import(*errBkgPdf);

  // --- draw plot - Sig ---
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *frSig = ctau3DErr->frame();
    binSubtractedSIG->plotOn(frSig, DataError(RooAbsData::SumW2), Name("data"));
    errSigPdf->plotOn(frSig, Name("model"), LineColor(kRed));

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frSig->getHist("data")))
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
    frSig->SetMinimum(ymin * 0.5);
    frSig->SetMaximum(std::max(ymax, ymin) * 1e4);

    frSig->GetYaxis()->SetTitle("Events");
    frSig->GetXaxis()->SetTitle("");
    frSig->Draw("e");

    // legend
    TLegend leg(0.49, 0.65, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(frSig, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frSig, "model"))
        leg.AddEntry(o, "Sig", "pe");

      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      // tx.DrawLatex(x, y0 + dy * k++, Form("%.3f < c#tau < %.3f", ctMin, ctMax));
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = frSig->pullHist("data", "model");
    RooPlot *fpull = ctau3DErr->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr [mm]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(ctErrLow, 0.0, ctErrHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    c.SaveAs(Form("%s/ctau_err_sig_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // draw plot - Bkg
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *frBkg = ctau3DErr->frame();
    binScaledBKG->plotOn(frBkg, DataError(RooAbsData::SumW2), Name("data"));
    errBkgPdf->plotOn(frBkg, Name("model"), LineColor(kBlue));

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frBkg->getHist("data")))
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
    frBkg->SetMinimum(ymin * 0.5);
    frBkg->SetMaximum(std::max(ymax, ymin) * 1e4);

    frBkg->GetYaxis()->SetTitle("Events");
    frBkg->GetXaxis()->SetTitle("");
    frBkg->Draw("e");

    // legend
    TLegend leg(0.49, 0.65, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(frBkg, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frBkg, "model"))
        leg.AddEntry(o, "Bkg", "pe");

      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      // tx.DrawLatex(x, y0 + dy * k++, Form("%.3f < c#tau < %.3f", ctMin, ctMax));
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = frBkg->pullHist("data", "model");
    RooPlot *fpull = ctau3DErr->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DErr [mm]");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    TLine line(ctErrLow, 0.0, ctErrHigh, 0.0);
    line.SetLineStyle(2);
    line.Draw("same");

    c.SaveAs(Form("%s/ctau_err_bkg_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // === ctau res ===
  cout << "\n=== ctau res ===\n";
  // --- build ctau resolution model (3 gauss) with PEE ---
  RooRealVar ctMean("ctMean", "mean", 0.0, -0.05, 0.05);
  ctMean.setConstant(kTRUE); // fix for MC, free for Data fit
  RooRealVar ctSigma1("ctSigma1","",1.00,0.001, 5.0);

  RooRealVar ctSigma12("ctSigma12", "", 1.1, 1, 10);
  RooRealVar ctSigma23("ctSigma23", "", 1.1, 1, 10);
  RooFormulaVar ctSigma2("ctSigma2", "@0*@1", RooArgList(ctSigma1, ctSigma12));
  RooFormulaVar ctSigma3("ctSigma3", "@0*@1", RooArgList(ctSigma2, ctSigma23));

  RooGaussian ctG1("ctG1","ctau gauss", *ws->var("ctau3DRes"), ctMean, ctSigma1);
  // *ws->var("ctau3DRes")
  RooGaussian ctG2("ctG2","", *ws->var("ctau3DRes"), ctMean, ctSigma2);
  RooGaussian ctG3("ctG3","", *ws->var("ctau3DRes"), ctMean, ctSigma3);

  RooRealVar   ctF1 ("ctF1" ,"frac1" , 0.70, 0.0, 1.0);
  RooRealVar   ctF2p("ctF2p","The p stands for recursive Parameter'", 0.83, 0.0, 1.0);
  RooFormulaVar ctF2("ctF2","(1-@0)*@1", RooArgList(ctF1, ctF2p));
  // RooFormulaVar ctF3("ctF3","(1-@0)*(1-@1)", RooArgList(ctF1, ctF2p));

  RooAddPdf ctResFitModel("ctResFitModel", "", RooArgList(ctG1, ctG2, ctG3), RooArgList(ctF1, ctF2));

  // --- combine with signal PEE ---
  // skip when using the ctau3DRes

  // --- fit ---
  // RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-7);
  // RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-9);

  RooFitResult *fitRes;
  if (isSkipRes && !gSystem->AccessPathName(Form("%s/ctau_res_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/ctau_res_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitRes", tmp);
    if (tmp) fitRes = (RooFitResult*) tmp->Clone("fitRes");

    syncVar(ctSigma1, *fitRes); 
    syncVar(ctSigma12, *fitRes);
    syncVar(ctSigma23, *fitRes);
    syncVar(ctF1, *fitRes);
    syncVar(ctF2p, *fitRes);
  } 
  else {
    fitRes = ctResFitModel.fitTo(*redPRMC, Range("ctResRange"), Offset(true), Save(), PrintLevel(-1), Strategy(2), PrintEvalErrors(-1), RecoverFromUndefinedRegions(2), NumCPU(32)); //, EvalBackend("legacy"),  Range("massRange"),
  }

  // --- draw plots ---
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // --- frame & plot ---
    RooPlot *fr = ws->var("ctau3DRes")->frame(Title("")); // , Bins(nBins)
    redPRMC->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    ctResFitModel.plotOn(fr, Name("model"));
    ctResFitModel.plotOn(fr, Components("ctG1"), Name("G1"), LineStyle(kDotted), LineColor(kRed));
    ctResFitModel.plotOn(fr, Components("ctG2"), Name("G2"), LineStyle(kDotted), LineColor(kOrange));
    ctResFitModel.plotOn(fr, Components("ctG3"), Name("G3"), LineStyle(kDotted), LineColor(kGreen));

    // --- dynamic y-range for log scale ---
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // --- legend ---
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "G1"))
        leg.AddEntry(o, "Gauss1", "pe");
      if (auto *o = findObj(fr, "G2"))
        leg.AddEntry(o, "Gauss2", "pe");
      if (auto *o = findObj(fr, "G3"))
        leg.AddEntry(o, "Gauss3", "pe");

      leg.Draw("same");
    }

    // --- CMS/info latex ---
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Prompt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitRes->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/ctau_res_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // --- parameter latex ---
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;
      
      // lambda function for printing
      auto print = [&](const char *title, const char *vname)
      {
        auto *v = dynamic_cast<RooRealVar *>(ctResFitModel.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("mean", "ctMean");
      print("#sigma1", "ctSigma1");
      print("#sigma_{2/1}", "ctSigma12");
      print("#sigma_{3/2}", "ctSigma23");
      print("f_{1}", "ctF1");
      print("f_{2p}", "ctF2p");
    }

    // --- pull pad ---
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ws->var("ctau3DRes")->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3DRes");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    // TLine line(massLow, 0.0, massHigh, 0.0);
    // line.SetLineStyle(2);
    // line.Draw("same");

    // --- chi2/ndf ---
    if (fitRes)
    {
      int npar = fitRes->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/ctau_res_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // --- fix parameters ---
  // ctMean.setConstant(false);
  ctSigma12.setConstant();
  ctSigma23.setConstant();
  ctF1.setConstant();
  ctF2p.setConstant();

  // --- save results ---
  fitRes->Print("V");
  TFile outRes(Form("%s/ctau_res_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitRes->Write("fitRes");
  outRes.Close();


  // === ctau bkg ===
  cout << "\n=== ctau bkg ===\n";
  // RooRealVar mean("mean", "resolution mean", 0.0, -0.1, 0.1);
  RooRealVar resMean("resMean", "resolution mean", 0.0);
  RooRealVar resSigma1("resSigma1", "resolution sigma1", 0.05, 0.001, 0.5);

  RooRealVar resSigma1_meas("resSigma1_meas", "Measured resSigma1", 0.1);
  RooRealVar resSigma1_sigma("resSigma1_sigma", "resSigma1 uncertainty", 0.01);
  RooGaussian resSigma1_constrain(
      "resSigma1_constrain",
      "Gaussian constraint on resSigma1",
      resSigma1, // variable
      resSigma1_meas, // expected mean
      resSigma1_sigma //expected sigma
  );

  RooRealVar resSigma12("resSigma12", "sigma2/sigma1", 2, 1.0, 3.0);
  RooRealVar resSigma23("resSigma23", "sigma3/sigma2", 2.0, 1.0, 5.0);

  // sigma2 = sigma1 * sigmaRatio
  RooFormulaVar resSigma2("resSigma2", "@0 * @1", RooArgList(resSigma1, ctSigma12));
  RooFormulaVar resSigma3("resSigma3", "@0 * @1", RooArgList(resSigma2, ctSigma23));

  RooGaussModel resG1("resG1", "Gauss Res1", *ctau3D, resMean, resSigma1, *ctau3DErr);
  RooGaussModel resG2("resG2", "Gauss Res2", *ctau3D, resMean, resSigma2, *ctau3DErr);
  RooGaussModel resG3("resG3", "Gauss Res2", *ctau3D, resMean, resSigma3, *ctau3DErr);

  RooRealVar   resF1 ("resF1" ,"frac1" , 0.8, 0, 1.0);
  RooRealVar   resF2p("resF2p","The p stands for recursive Parameter'", 0.83, 0.0, 1.0);
  RooFormulaVar resF2("resF2","(1-@0)*@1", RooArgList(resF1, resF2p));

  RooAddModel bkgResModel("bkgResModel", "", RooArgList(resG1, resG2, resG3), RooArgList(ctF1, ctF2));

  RooRealVar tauL("tauL", "lifetime left", 7.5556e-03, 0.001, 0.1);

  RooRealVar tauC1("tauC1", "lifetime central1", 0.05, 0.0001, 0.1);
  RooRealVar tauC2("tauC2", "lifetime central2", 0.09, 0.0001, 0.5);

  RooRealVar tauR("tauR", "lifetime right", 4.4341e-01, 0.01, 0.5);
  RooRealVar tauSigNp("tauSigNp", "lifetime right", 0.1, 0.01, 0.5);

  RooRealVar fBkgC1("fBkgC1", "frac central1", 0.3, 0.0, 1.0);   // decayC1 vs decayC2
  RooRealVar fBkgC("fBkgC", "frac total central", 0.55, 0.0, 1.0);
  RooRealVar fBkgL("fBkgL", "frac left", 0.07, 0.0, 1.0);

  RooDecay decayL ("decayL",  "Left decay", *ctau3D, tauL,  bkgResModel, RooDecay::Flipped);
  RooDecay decayC1("decayC1", "Central decay 1", *ctau3D, tauC1, bkgResModel, RooDecay::DoubleSided);
  RooDecay decayC2("decayC2", "Central decay 2", *ctau3D, tauC2, bkgResModel, RooDecay::DoubleSided);
  RooDecay decayR ("decayR",  "Right decay", *ctau3D, tauR, bkgResModel, RooDecay::SingleSided);

  RooProdPdf decayL_PEE ("decayL_PEE",  "L * err",  RooArgSet(*errBkgPdf),
                        RooFit::Conditional(RooArgSet(decayL),  RooArgSet(*ctau3D)));
  RooProdPdf decayC1_PEE("decayC1_PEE", "C1 * err", RooArgSet(*errBkgPdf),
                        RooFit::Conditional(RooArgSet(decayC1), RooArgSet(*ctau3D)));
  RooProdPdf decayC2_PEE("decayC2_PEE", "C2 * err", RooArgSet(*errBkgPdf),
                        RooFit::Conditional(RooArgSet(decayC2), RooArgSet(*ctau3D)));
  RooProdPdf decayR_PEE("decayR_PEE",  "R * err",  RooArgSet(*errBkgPdf),
                        RooFit::Conditional(RooArgSet(decayR),  RooArgSet(*ctau3D)));

  RooAddPdf decayCentral_PEE("decayCentral_PEE", "central with err",
                        RooArgList(decayC1_PEE, decayC2_PEE),
                        RooArgList(fBkgC1)); 

  RooAddPdf ctBkgFitModel("ctBkgFitModel", "", RooArgList(decayCentral_PEE, decayL_PEE, decayR_PEE), RooArgList(fBkgC, fBkgL));

  // --- fit ---
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  RooFitResult *fitBkg;
  if (isSkipBkg && !gSystem->AccessPathName(Form("%s/ctau_bkg_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/ctau_bkg_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitBkg", tmp);
    if (tmp) fitBkg = (RooFitResult*) tmp->Clone("fitBkg");

    syncVar(ctF1,*fitBkg);
    syncVar(ctF2p,*fitBkg);
    syncVar(ctSigma12,*fitBkg);
    syncVar(ctSigma23,*fitBkg);
    syncVar(resMean,*fitBkg);
    syncVar(fBkgC,*fitBkg);
    syncVar(fBkgC1,*fitBkg);
    syncVar(fBkgL,*fitBkg);
    syncVar(resSigma1,*fitBkg);
    syncVar(tauC1,*fitBkg);
    syncVar(tauC2,*fitBkg);
    syncVar(tauL,*fitBkg);
    syncVar(tauR,*fitBkg);
  } 
  else
    fitBkg = ctBkgFitModel.fitTo(*redDataSB, Save(), Optimize(0), Offset(), NumCPU(32), EvalBackend("legacy"),RecoverFromUndefinedRegions(1.5), PrintEvalErrors(-1), PrintLevel(-1), RooFit::ExternalConstraints(resSigma1_constrain));
  fitBkg->Print("V");
  // EvalBackend("legacy"), ExternalConstraints(RooArgSet(c_resSigma1))

  // --- draw plots ---
  // plotting is very slow. Draw when you need it.
  if (!isSkipBkg) {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *fr = ctau3D->frame(Title("")); // , Bins(nBins)
    redDataSB->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    ctBkgFitModel.plotOn(fr, Name("model"), NumCPU(16), ProjWData(RooArgList(*ctau3DErr), *redDataSB), Normalization(redDataSB->sumEntries(), RooAbsReal::NumEvent), Precision(1e-4));
    // ProjWData(RooArgList(*ctau3DErr), *redDataSB)
    
    // components - precision -1 -> to save time
    // ProjWData makes error for components - result shold be same though drawing speed is slow. 
    // -> draw components for final plot

    // ctBkgFitModel.plotOn(fr, Name("bkgC"), Components("decayCentral_PEE"), NumCPU(4), Normalization(redDataSB->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kRed), ConditionalObservables(*ctau3DErr));
    // ctBkgFitModel.plotOn(fr, Name("bkgL"), Components("decayL_PEE"), NumCPU(4), Normalization(redDataSB->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kMagenta));
    // ctBkgFitModel.plotOn(fr, Name("bkgR"), Components("decayR_PEE"), NumCPU(4), Normalization(redDataSB->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kOrange));

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    if (ymin <= 0 || ymin >= 1e300)
      ymin = 1e-3;
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);
    // fr->SetMaximum(std::max(ymax, ymin));

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // legend
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "bkgC"))
        leg.AddEntry(o, "2 DecayC", "pe");
      if (auto *o = findObj(fr, "bkgL"))
        leg.AddEntry(o, "DecayL", "pe");
      if (auto *o = findObj(fr, "bkgR"))
        leg.AddEntry(o, "DecayR", "pe");
      
      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitRes->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/ctau_bkg_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // parameter latex
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;
      
      // lambda function for printing
      auto print = [&](const char *title, const char *vname)
      {
        auto *v = dynamic_cast<RooRealVar *>(ctBkgFitModel.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("mean", "resMean");
      print("#sigma1", "resSigma1");
      print("#sigma_{2/1}", "ctSigma12");
      print("#sigma_{3/2}", "ctSigma23");

      print("f_{Res1}", "ctF1");
      print("f_{Res2}", "ctF2p");

      print("f_{BkgL}", "fBkgL");
      print("f_{BkgC}", "fBkgC");
      print("f_{BkgC1}", "fBkgC1");

      print("#tau_{L}", "tauL");
      print("#tau_{R}", "tauR");
      print("#tau_{C1}", "tauC1");
      print("#tau_{C2}", "tauC2");
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ctau3D->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3D");
    fpull->GetXaxis()->CenterTitle();
    // fpull->SetMinimum(-8);
    // fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    // TLine line(massLow, 0.0, massHigh, 0.0);
    // line.SetLineStyle(2);
    // line.Draw("same");

    // chi2/ndf
    if (fitBkg)
    {
      int npar = fitBkg->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/ctau_bkg_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }
  
  // --- fix parameters ---
  fBkgC.setConstant();
  fBkgC1.setConstant();
  fBkgL.setConstant();
  // resSigma1.setConstant();
  tauC1.setConstant();
  tauC2.setConstant();
  tauL.setConstant();
  tauR.setConstant();

  // --- save results ---
  // fitBkg->Print("V");
  TFile outBkg(Form("%s/ctau_bkg_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitBkg->Write("fitBkg");
  outBkg.Close();


  // === ctau true ===
  cout << "\n=== ctau true ===\n";
  // --- build mass signal model ---
  // --- build bkg signal model ---

  // RooRealVar mean("mean", "resolution mean", 0.0, -0.1, 0.1);
  RooRealVar trueMean("trueMean", "resolution mean", 0.0);
  RooRealVar trueSigma1("trueSigma1", "resolution sigma1", 0.05, 0.0000001, 0.5);

  RooTruthModel deltaFcn("deltaFcn", "", *ctau3Dtrue);

  RooRealVar trTau1("trTau1", "lifetime left", 7.5556e-03, 0.001, 10);
  RooRealVar trTau2("trTau2", "lifetime left", 1, 0.001, 10);

  // RooRealVar tauR("tauR", "lifetime right", 4.4341e-01, 0.01, 0.5);
  // RooRealVar tauSigNp("tauSigNp", "lifetime right", 0.1, 0.01, 0.5);

  // RooRealVar fBkgC1("fBkgC1", "frac central1", 0.3, 0.0, 1.0);   // decayC1 vs decayC2
  // RooRealVar fBkgC("fBkgC", "frac total central", 0.55, 0.0, 1.0);
  RooRealVar fTrue("fTrue", "", 0.07, 0.0, 1.0);

  // RooDecay decayL ("decayL",  "Left decay", *ctau3D, tauL,  bkgResModel, RooDecay::Flipped);
  // RooDecay decayC1("decayC1", "Central decay 1", *ctau3D, tauC1, bkgResModel, RooDecay::DoubleSided);
  RooDecay trueR1 ("trueR1",  "Right decay", *ctau3Dtrue, trTau1, deltaFcn, RooDecay::SingleSided);
  RooDecay trueR2 ("trueR2",  "Right decay", *ctau3Dtrue, trTau2, deltaFcn, RooDecay::SingleSided);

  RooAddPdf ctTrueFitModel("ctTrueFitModel", "", RooArgList(trueR1, trueR2), RooArgList(fTrue));

  // --- fit ---
  for (int i = 0; i < 3; i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Eval);
  }
  // Tracing
  for (int i = 0; i < RooMsgService::instance().numStreams(); i++)
  {
    RooMsgService::instance().getStream(i).removeTopic(RooFit::Tracing);
  }
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  RooFitResult *fitTrue;
  if (isSkipTrue && !gSystem->AccessPathName(Form("%s/ctau_true_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/ctau_true_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitTrue", tmp);
    if (tmp) fitTrue = (RooFitResult*) tmp->Clone("fitTrue");

    syncVar(trTau1,*fitTrue);
    syncVar(trTau2,*fitTrue);
    syncVar(fTrue,*fitTrue);
  } 
  else
    fitTrue = ctTrueFitModel.fitTo(*redNPMC, Save(), Optimize(0), Offset(), NumCPU(32), EvalBackend("legacy"),RecoverFromUndefinedRegions(1.5), PrintEvalErrors(-1), PrintLevel(-1));
  fitTrue->Print("V");
  
  // --- draw plots ---
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *fr = ctau3Dtrue->frame(Title("")); // , Bins(nBins)
    redNPMC->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    ctTrueFitModel.plotOn(fr, Name("model"), NumCPU(16));

    
    ctTrueFitModel.plotOn(fr, NumCPU(4), Name("trueR1"), Components("trueR1"), LineStyle(kDotted), LineColor(kRed));
    ctTrueFitModel.plotOn(fr, NumCPU(4), Name("trueR2"), Components("trueR2"), LineStyle(kDotted), LineColor(kMagenta));
    // ProjWData(RooArgList(*ctau3DErr), *redNPMC)
    
    // components - precision -1 -> to save time
    // ProjWData makes error for components - result shold be same though drawing speed is slow. 
    // -> draw components for final plot

    // ctBkgFitModel.plotOn(fr, Name("bkgC"), Components("decayCentral_PEE"), NumCPU(4), Normalization(redNPMC->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kRed), ConditionalObservables(*ctau3DErr));
    // ctBkgFitModel.plotOn(fr, Name("bkgL"), Components("decayL_PEE"), NumCPU(4), Normalization(redNPMC->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kMagenta));
    // ctBkgFitModel.plotOn(fr, Name("bkgR"), Components("decayR_PEE"), NumCPU(4), Normalization(redNPMC->sumEntries(), RooAbsReal::NumEvent), Precision(-1), LineStyle(kDotted), LineColor(kOrange));

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    if (ymin <= 0 || ymin >= 1e300)
      ymin = 1e-3;
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);
    // fr->SetMaximum(std::max(ymax, ymin));

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // legend
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "trueR1"))
        leg.AddEntry(o, "Decay1", "pe");
      if (auto *o = findObj(fr, "trueR2"))
        leg.AddEntry(o, "Decay2", "pe");
      // if (auto *o = findObj(fr, "bkgR"))
      //   leg.AddEntry(o, "DecayR", "pe");
      
      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Nonpromt MC, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitRes->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/ctau_bkg_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // parameter latex
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;
      
      // lambda function for printing
      auto print = [&](const char *title, const char *vname)
      {
        auto *v = dynamic_cast<RooRealVar *>(ctTrueFitModel.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("f_{1}", "fTrue");
      print("#tau_{1}", "trTau1");
      print("#tau_{2}", "trTau2");
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ctau3Dtrue->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3Dtrue");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    // TLine line(massLow, 0.0, massHigh, 0.0);
    // line.SetLineStyle(2);
    // line.Draw("same");

    // chi2/ndf
    if (fitTrue)
    {
      int npar = fitTrue->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/ctau_true_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }
  // --- fix parameters ---
  trTau1.setConstant();
  trTau2.setConstant();
  fTrue.setConstant();

  // --- save results ---
  TFile outTrue(Form("%s/ctau_true_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitTrue->Write("fitTrue");
  outTrue.Close();


  // === final fit ===
  cout << "\n=== final fit ===\n";

  // --- sig prompt model ---
  // --------------------------------------------
  // // massSig * (res * sigPee)
  // RooProdPdf sigRes_PEE ("sigRes_PEE",  "",  RooArgSet(*errSigPdf),
  //                       RooFit::Conditional(RooArgSet(bkgResModel),  RooArgSet(*ctau3D))); // name is not a matter
  // RooProdPdf ctMaSigPr ("ctMaSigPr",  "",  RooArgSet(massSigPdf, sigRes_PEE));
  
  // // --- sig nonprompt model ---
  // // np decay models with res
  // RooDecay sigNp1 ("sigNp1",  "Right decay", *ctau3D, trTau1, bkgResModel, RooDecay::SingleSided);
  // RooDecay sigNp2 ("sigNp2",  "Right decay", *ctau3D, trTau2, bkgResModel, RooDecay::SingleSided);

  // RooProdPdf sigNp1_PEE ("sigNp1_PEE",  "",  RooArgSet(*errSigPdf),
  //                       RooFit::Conditional(RooArgSet(sigNp1), RooArgSet(*ctau3D)));
  // RooProdPdf sigNp2_PEE ("sigNp2_PEE",  "",  RooArgSet(*errSigPdf),
  //                       RooFit::Conditional(RooArgSet(sigNp2), RooArgSet(*ctau3D)));

  // // mass * ctau
  // RooProdPdf ctMaSigNp1 ("ctMaSigNp1",  "",  RooArgSet(massSigPdf, sigNp1_PEE));
  // RooProdPdf ctMaSigNp2 ("ctMaSigNp2",  "",  RooArgSet(massSigPdf, sigNp2_PEE));
  
  // // sum
  // RooAddPdf ctMaNpTot("ctMaNpTot", "sigNp1 + sigNp2",
  //                  RooArgList(ctMaSigNp1, ctMaSigNp2),
  //                  RooArgList(fTrue));
  
  // // --- bkg model ---
  // RooProdPdf ctMaSigBkgL("ctMaSigBkgL",  "",  RooArgSet(massBkgPdf, decayL_PEE));
  // RooProdPdf ctMaSigBkgC1("ctMaSigBkgC1",  "",  RooArgSet(massBkgPdf, decayC1_PEE));
  // RooProdPdf ctMaSigBkgC2("ctMaSigBkgC2",  "",  RooArgSet(massBkgPdf, decayC2_PEE));
  // RooProdPdf ctMaSigBkgR("ctMaSigBkgR",  "",  RooArgSet(massBkgPdf, decayR_PEE));

  // // sum
  // RooAddPdf ctMaSigBkgC("ctMaSigBkgC", "central with err",
  //                       RooArgList(ctMaSigBkgC1, ctMaSigBkgC2),
  //                       RooArgList(fBkgC1)); 

  // RooAddPdf ctMaBkgTot("ctMaBkgTot", "", RooArgList(ctMaSigBkgC, ctMaSigBkgL, ctMaSigBkgR), RooArgList(fBkgC, fBkgL));

  // // --- final model ---
  // RooRealVar bFrac("bFrac", "", 0.5, 0, 1);
  // RooAddPdf ctMassSig("ctMassSig", "", RooArgList(ctMaSigPr, ctMaNpTot), RooArgList(bFrac));

  // RooFormulaVar fSignal("fSignal", "@0/(@0+@1)", {nSigMass, nBkgMass});
  // RooAddPdf finalFitModel("finalFitModel", "", RooArgList(ctMassSig, ctMaBkgTot), RooArgList(fSignal));
  // --------------------------------------------


  // --------------------------------------------
  RooProdPdf massCtPR("massCtPR","", RooArgSet(massSigPdf, bkgResModel));
  // rf307식 PEE:  err × Conditional(본체, 관측들)
  RooProdPdf MassCtPR_PEE("MassCtPR_PEE","PDF with PEE",
                          *errSigPdf,
                          RooFit::Conditional(massCtPR, RooArgSet(*ctau3D, *mass)));

  // =============== Signal: Non-prompt(2항) ===============
  // ctau part (with resolution)
  RooDecay sigNp1("sigNp1","NP1", *ctau3D, trTau1, bkgResModel, RooDecay::SingleSided);
  RooDecay sigNp2("sigNp2","NP2", *ctau3D, trTau2, bkgResModel, RooDecay::SingleSided);

  // mass × ctau
  RooProdPdf massCtNP1("massCtNP1","", RooArgSet(massSigPdf, sigNp1));
  RooProdPdf massCtNP2("massCtNP2","", RooArgSet(massSigPdf, sigNp2));

  // rf307식 PEE
  RooProdPdf massCtNP1_PEE("massCtNP1_PEE","", *errSigPdf,
                          RooFit::Conditional(massCtNP1, RooArgSet(*ctau3D, *mass)));
  RooProdPdf massCtNP2_PEE("massCtNP2_PEE","", *errSigPdf,
                          RooFit::Conditional(massCtNP2, RooArgSet(*ctau3D, *mass)));

  // sum (NP 내부 혼합)
  RooAddPdf massCtNP_TOT_PEE("massCtNP_TOT_PEE","NP1+NP2 (PEE)",
                            RooArgList(massCtNP1_PEE, massCtNP2_PEE),
                            RooArgList(fTrue));

  // Sig 합 (Prompt vs Non-prompt)
  RooRealVar bFrac("bFrac", "", 0.5, 0, 1);
  RooAddPdf massCtSig_PEE("massCtSig_PEE","Sig PR+NP (PEE)",
                          RooArgList(massCtNP_TOT_PEE, MassCtPR_PEE),
                          RooArgList(bFrac)); // bFrac = NP 비중

  // =============== Background ===============
  RooProdPdf massCtBkgL ("massCtBkgL","",  RooArgSet(massBkgPdf, decayL));
  RooProdPdf massCtBkgC1("massCtBkgC1","", RooArgSet(massBkgPdf, decayC1));
  RooProdPdf massCtBkgC2("massCtBkgC2","", RooArgSet(massBkgPdf, decayC2));
  RooProdPdf massCtBkgR ("massCtBkgR","",  RooArgSet(massBkgPdf, decayR));

  RooProdPdf massCtBkgL_PEE ("massCtBkgL_PEE","",  *errBkgPdf,
                            RooFit::Conditional(massCtBkgL,  RooArgSet(*ctau3D, *mass)));
  RooProdPdf massCtBkgC1_PEE("massCtBkgC1_PEE","", *errBkgPdf,
                            RooFit::Conditional(massCtBkgC1, RooArgSet(*ctau3D, *mass)));
  RooProdPdf massCtBkgC2_PEE("massCtBkgC2_PEE","", *errBkgPdf,
                            RooFit::Conditional(massCtBkgC2, RooArgSet(*ctau3D, *mass)));
  RooProdPdf massCtBkgR_PEE ("massCtBkgR_PEE","",  *errBkgPdf,
                            RooFit::Conditional(massCtBkgR,  RooArgSet(*ctau3D, *mass)));

  RooAddPdf massCtBkgC_PEE("massCtBkgC_PEE","central (PEE)",
                          RooArgList(massCtBkgC1_PEE, massCtBkgC2_PEE),
                          RooArgList(fBkgC1));

  RooAddPdf massCtBkgTot_PEE("massCtBkgTot_PEE","Bkg TOT (PEE)",
                            RooArgList(massCtBkgC_PEE, massCtBkgL_PEE, massCtBkgR_PEE),
                            RooArgList(fBkgC, fBkgL));

  // =============== Final ===============
  RooFormulaVar fSignal("fSignal","@0/(@0+@1)", RooArgList(nSigMass, nBkgMass));
  RooAddPdf finalFitModel("finalFitModel","Sig+Bkg (PEE)",
                          RooArgList(massCtSig_PEE, massCtBkgTot_PEE),
                          RooArgList(nSigMass, nBkgMass));
  // --------------------------------------------



  // --------------------------------------------
  // --- signal PR ---
  // massCtPR
  // RooProdPdf massCtPR("massCtPR",  "",  RooArgSet(massSigPdf, bkgResModel));

  // // massCtPR_PEE
  // RooProdPdf MassCtPR_PEE("MassCtPR_PEE", "PDF with PEE", *errSigPdf, Conditional(massCtPR, RooArgList(*ctau3D, *mass)));

  // // --- signal NP ---
  // // massCtNP
  // RooDecay sigNp1 ("sigNp1",  "Right decay", *ctau3D, trTau1, bkgResModel, RooDecay::SingleSided);
  // RooProdPdf massCtNP("massCtNP",  "",  RooArgSet(massSigPdf, sigNp1));
  
  // // massCtNP_PEE
  // RooProdPdf massCtNP_PEE("massCtNP_PEE", "PDF with PEE", *errSigPdf, Conditional(massCtNP, RooArgList(*ctau3D, *mass)));

  // // --- bkg ---
  // // massCtBkg
  // RooProdPdf massCtBkg("massCtBkg",  "",  RooArgSet(massBkgPdf, decayC1));

  // // massCtBkg_PEE
  // RooProdPdf massCtBkg_PEE("massCtBkg_PEE", "PDF with PEE", *errBkgPdf, Conditional(massCtBkg, RooArgList(*ctau3D, *mass)));


  // // final model
  // RooRealVar bFrac("bFrac", "", 0.5, 0, 1);
  // RooAddPdf massCtSig_PEE("massCtSig_PEE", "", RooArgList(massCtNP_PEE, MassCtPR_PEE), RooArgList(bFrac));

  // RooFormulaVar fSignal("fSignal", "@0/(@0+@1)", {nSigMass, nBkgMass});
  // RooAddPdf finalFitModel("finalFitModel", "", RooArgList(massCtSig_PEE, massCtBkg_PEE), RooArgList(fSignal));
  // --------------------------------------------


  //   RooProdPdf CtBkgTot_PEE("CtBkgTot_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
  //                           Conditional(*(ws->pdf("CtBkgTot")), RooArgList(*(ws->var("ctau3D")))));
  //   ws->import(CtBkgTot_PEE);

  //   RooProdPdf CtPR_PEE("CtPR_PEE", "CtPDF with PEE", *(ws->pdf("errPdfSig")),
  //                       Conditional(*(ws->pdf("CtPRRes")), RooArgList(*(ws->var("ctau3D")))));
  //   ws->import(CtPR_PEE);

  //   // // Build 2D PDF (mass x ctau) :: for step[[6]]
  //   // sprintf(funct, "PROD::MassCtPR(%s,CtPRRes)", "G1CB1Sig"); ws->factory(funct);
  //   // sprintf(funct, "PROD::MassCtNP(%s,CtNPTot)", "G1CB1Sig"); ws->factory(funct);
  //   // sprintf(funct, "PROD::MassCtBkg(%s,CtBkgTot)", "expBkg"); ws->factory(funct);

  //   // RooProdPdf MassCtPR_PEE("MassCtPR_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
  //   //                         Conditional(*(ws->pdf("MassCtPR")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  //   // ws->import(MassCtPR_PEE);
  //   // RooProdPdf MassCtNP_PEE("MassCtNP_PEE", "PDF with PEE", *(ws->pdf("errPdfSig")),
  //   //                         Conditional(*(ws->pdf("MassCtNP")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  //   // ws->import(MassCtNP_PEE);
  //   // RooProdPdf MassCtBkg_PEE("MassCtBkg_PEE", "PDF with PEE", *(ws->pdf("errPdfBkg")),
  //   //                          Conditional(*(ws->pdf("MassCtBkg")), RooArgList(*(ws->var("ctau3D")), *(ws->var("mass")))));
  //   // ws->import(MassCtBkg_PEE);

  // --- fit ---
  cout << "start final fit\n";
  resSigma1.setRange(0.001, 1);

  // RooRealVar resSigma1_meas2("resSigma1_meas2", "Measured resSigma1", 0.05);
  // RooRealVar resSigma1_sigma2("resSigma1_sigma2", "resSigma1 uncertainty", 0.01);
  // RooGaussian resSigma1_constrain2(
  //     "resSigma1_constrain2",
  //     "Gaussian constraint on resSigma1",
  //     resSigma1, // variable
  //     resSigma1_meas2, // expected mean
  //     resSigma1_sigma2 //expected sigma
  // );

  output_current_time();
  RooFitResult *fitFinal;
  if (isSkipFinal && !gSystem->AccessPathName(Form("%s/final_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), kReadPermission)) {
    TFile fin(Form("%s/final_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "READ");
    RooFitResult* tmp = nullptr; fin.GetObject("fitFinal", tmp);
    if (tmp) fitFinal = (RooFitResult*) tmp->Clone("fitFinal");

    syncVar(bFrac,*fitFinal);
    syncVar(resSigma1,*fitFinal);
    syncVar(trTau1,*fitFinal);
  } 
  else
    fitFinal = finalFitModel.fitTo(*redData, Save(), Offset(), NumCPU(32), EvalBackend("legacy"),RecoverFromUndefinedRegions(1.5), PrintEvalErrors(-1), PrintLevel(-1));
  fitFinal->Print("V");
  

  // --- draw final mass ---
  cout << "draw final mass\n";
  output_current_time();
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *fr = mass->frame(Title("")); // , Bins(nBins)
    
    // RooDataHist hProj("hProj","", RooArgList(*ctau3DErr), *redData);

    // fr->SetNbins(25);
    redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));
    finalFitModel.plotOn(fr, Name("model"), NumCPU(16));

    // finalFitModel.plotOn(fr, Name("model"), NumCPU(16), ProjWData(RooArgSet(*ctau3D, *ctau3DErr), *redData),LineColor(kBlue));
    
    // finalFitModel.plotOn(fr, Name("model"), NumCPU(16), ProjWData(RooArgSet(*ctau3D, *ctau3DErr), hProj), Normalization(hProj.sumEntries(), RooAbsReal::NumEvent));

    // finalFitModel.plotOn(fr, Name("model"), NumCPU(16), ProjWData(RooArgSet(*ctau3D, *ctau3DErr), hProj), Normalization(hProj.sumEntries(), RooAbsReal::NumEvent));
    // , Range("ctauFinal"), NormRange("ctauFinal")

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    if (ymin <= 0 || ymin >= 1e300)
      ymin = 1e-3;
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);
    // fr->SetMaximum(std::max(ymax, ymin));

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // legend
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "trueR1"))
        leg.AddEntry(o, "Decay1", "pe");
      if (auto *o = findObj(fr, "trueR2"))
        leg.AddEntry(o, "Decay2", "pe");
      // if (auto *o = findObj(fr, "bkgR"))
      //   leg.AddEntry(o, "DecayR", "pe");
      
      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
      
      // fit status
      int st = fitRes->status(); // 0 = success
      if (st != 0)
      {
        tx.DrawLatex(x, y0 + dy * k++, Form("Status=%d", st));
        std::ofstream flog(Form("logs%s/final_pT%.1f_%.1f.txt", userLabel.Data(), ptLow, ptHigh), std::ios::app);
        flog.close();
      }
    }

    // parameter latex
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;
      
      // lambda function for printing
      auto print = [&](const char *title, const char *vname)
      {
        auto *v = dynamic_cast<RooRealVar *>(ctTrueFitModel.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("f_{1}", "fTrue");
      print("#tau_{1}", "trTau1");
      print("#tau_{2}", "trTau2");
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = mass->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3Dtrue");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    // TLine line(massLow, 0.0, massHigh, 0.0);
    // line.SetLineStyle(2);
    // line.Draw("same");

    // chi2/ndf
    if (fitFinal)
    {
      int npar = fitFinal->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/final_mass_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // --- draw final ctau ---
  cout << "draw final ctau\n";
  output_current_time();
  {
    TCanvas c("c", "c", 800, 800);
    TPad pad1("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
    pad1.SetBottomMargin(0.00001);
    pad1.SetLogy();
    pad1.Draw();
    pad1.cd();

    // frame & plot
    RooPlot *fr = ctau3D->frame(Title("")); // , Bins(nBins)
    // fr->updateNormVars({*ctau3D, *mass, *ctau3DErr});

    mass->setBins(300);
    ctau3DErr->setBins(300);
    
    redData->plotOn(fr, DataError(RooAbsData::SumW2), Name("data"));

    // finalFitModel.plotOn(fr, Name("model"), NumCPU(16), Precision(1e-4));

    
    // RooDataHist hProj("hProj","", RooArgList(*ctau3DErr, *mass), *redData);
    // finalFitModel.plotOn(fr, Name("model"), NumCPU(16), ProjWData(RooArgSet(*mass, *ctau3DErr), hProj),  Normalization(hProj.sumEntries(), RooAbsReal::NumEvent));

    RooDataHist hProj("hProj","", RooArgList(*mass, *ctau3DErr), *redData);
    finalFitModel.plotOn(fr, Name("model"));
    // , Range("ctauFinal"), NormRange("ctauFinal")

    // fr->updateNormVars({*ctau3D});
    // fr->updateNormVars({*mass, *ctau3D});
    // fr->updateNormVars({*ctau3DErr, *mass, *ctau3D});
    

    // dynamic y-range for log scale
    double ymin = 1e300, ymax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(fr->getHist("data")))
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
    if (ymin <= 0 || ymin >= 1e300)
      ymin = 1e-3;
    fr->SetMinimum(ymin * 0.5);
    fr->SetMaximum(std::max(ymax, ymin) * 1e4);
    // fr->SetMaximum(std::max(ymax, ymin));

    fr->GetYaxis()->SetTitle("Events");
    fr->GetXaxis()->SetTitle("");
    fr->Draw("e");

    // legend
    TLegend leg(0.49, 0.66, 0.70, 0.94);
    {
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      leg.SetTextSize(0.03);
      if (auto *o = findObj(fr, "data"))
        leg.AddEntry(o, "Data", "lep");
      if (auto *o = findObj(fr, "model"))
        leg.AddEntry(o, "Model", "pe");
      if (auto *o = findObj(fr, "trueR1"))
        leg.AddEntry(o, "Decay1", "pe");
      if (auto *o = findObj(fr, "trueR2"))
        leg.AddEntry(o, "Decay2", "pe");
      // if (auto *o = findObj(fr, "bkgR"))
      //   leg.AddEntry(o, "DecayR", "pe");
      
      leg.Draw("same");
    }

    // CMS/info latex
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double x = 0.19, y0 = 0.90, dy = -0.06;
      int k = 0;
      tx.DrawLatex(x, y0 + dy * k++, "CMS pp Ref. #sqrt{s_{NN}} = 5.02 TeV");
      tx.DrawLatex(x, y0 + dy * k++, "Data, J/#psi #rightarrow #mu^{+}#mu^{-}");
      if (yLow == 0)
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(x, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < y < %.1f", ptLow, ptHigh, yLow, yHigh));
    }

    // parameter latex
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double x = 0.71, y0 = 0.91, dy = -0.04;
      int k = 0;
      
      // lambda function for printing
      auto print = [&](const char *title, const char *vname)
      {
        auto *v = dynamic_cast<RooRealVar *>(ctTrueFitModel.getVariables()->find(vname));
        if (!v)
          return;

        const double val = v->getVal(), err = v->getError();
        const double eps = 1e-9 * (1.0 + std::fabs(val));
        const bool fixed = v->isConstant();
        const bool atMin = v->hasMin() && std::fabs(val - v->getMin()) <= eps;
        const bool atMax = v->hasMax() && std::fabs(val - v->getMax()) <= eps;
        const bool atBound = atMin || atMax;

        TString note;
        if (fixed)
          note += "(fixed)";
        if (atBound)
          note += note.IsNull() ? (atMin ? "(at min)" : "(at max)")
                                : (atMin ? ", (at min)" : ", (at max)");

        const Int_t oldColor = tp.GetTextColor();
        if (atBound)
          tp.SetTextColor(kRed + 1);

        if (fixed)
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g %s", title, val, note.Data()));
        else
          tp.DrawLatex(x, y0 + dy * k++, Form("%s = %.4g #pm %.3g %s", title, val, err, note.Data()));

        tp.SetTextColor(oldColor);
      };

      // print
      print("f_{1}", "fTrue");
      print("#tau_{1}", "trTau1");
      print("#tau_{2}", "trTau2");
    }

    // pull pad
    c.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.4);
    pad2.Draw();
    pad2.cd();

    RooHist *hpull = fr->pullHist("data", "model");
    RooPlot *fpull = ctau3D->frame(Title(""));
    fpull->addPlotable(hpull, "P");
    fpull->GetYaxis()->SetTitle("Pull");
    fpull->GetXaxis()->SetTitle("ctau3Dtrue");
    fpull->GetXaxis()->CenterTitle();
    fpull->SetMinimum(-8);
    fpull->SetMaximum(8);
    fpull->GetYaxis()->SetNdivisions(505);
    fpull->GetYaxis()->SetTitleSize(0.12);
    fpull->GetYaxis()->SetLabelSize(0.10);
    fpull->GetXaxis()->SetTitleSize(0.15);
    fpull->GetXaxis()->SetLabelSize(0.10);
    fpull->Draw();

    // TLine line(massLow, 0.0, massHigh, 0.0);
    // line.SetLineStyle(2);
    // line.Draw("same");

    // chi2/ndf
    if (fitFinal)
    {
      int npar = fitFinal->floatParsFinal().getSize();
      double chi2ndf = fr->chiSquare("model", "data", npar);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.10);
      tc.DrawLatex(0.82, 0.86, Form("#chi^{2}/ndf = %.2f", chi2ndf));
      cout << "\nchi2/ndf = " << chi2ndf << "\n\n";
    }
    c.SaveAs(Form("%s/final_ctau_pT%.1f_%.1f.png", figDir.Data(), ptLow, ptHigh));
  }

  // --- save results ---
  TFile outFinal(Form("%s/final_pT%.1f_%.1f.root", rootDir.Data(), ptLow, ptHigh), "RECREATE");
  fitFinal->Write("fitFinal");
  outFinal.Close();

  cout << "\n=== finish fit2D_step_by_step() ===\n";
  time.Stop();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n", time.RealTime(), time.CpuTime());
}