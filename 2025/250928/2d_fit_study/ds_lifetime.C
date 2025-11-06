////////-----------////////


#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooCurve.h"

#include <iostream> 
#include <fstream> 
#include <string> 
#include <algorithm> 
#include <vector> 
#include <iomanip> 
#include <cmath> 

//#include <cstdlib>
#include "RooFit.h"// 
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h" 
#include "RooDataSet.h" 
#include "RooDataHist.h" 
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooBifurGauss.h"
#include "RooLandau.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
// #include "RooNLLVar.h"
// #include "RooMinuit.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooHist.h"
//#include "RooDstarBG.h"     // for the threshold function
#include "TLine.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TString.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TPad.h" 
//#include "RooFitModels.h"
#include "RooDstD0BG.h"  // will be used for fitting deltam of combinatoric background 

using namespace RooFit ; 
using namespace std ;

int ds_lifetime()
{
	         
  RooRealVar t("t","t_{rec} (ns)",-0.002,0.004); //decay time t
  //  RooRealVar dt("dt","dt",0,0.0009); // error in decay time t
  RooRealVar dt("dt","dt",0,0.0009); 

  RooDataSet* data1 = (RooDataSet*)RooDataSet::read("dsm_new.txt",RooArgSet(t,dt),"Q");    
  data1->Print("V");

  RooPlot *frame = t.frame(Title("#tau D_{s}^{+}"), Bins(100)); // frame to plot lifetime t
  // data1->plotOn(frame,Binning(100), Name("data"));// plotting all the data on the frame*/	
	
//----------------------BLOCK 4------------------------//
//--------here start constructing the models to be used for fitting -----------//

  //-----------pdf for signal------------//

  //-------------decay time pdf--------------------//
//----------Resolution function for signal--------------//

//--------------gauss1---------------//
  RooRealVar mean("mean","mean",0.0001,-2,0.007); // 0.00000033
  RooRealVar sc1("sc1","sc1",1.0,0,4); // scaling factor 
  RooGaussModel gauss1("gauss1","gauss",t,mean,sc1,dt); // scaling/
  RooGaussModel gauss2("gauss2", "gauss", t, mean, sc1, dt); // scaling/

  //-------lifetime variable---------//
  RooRealVar tau("tau","tau",0.000104,0,0.001);
  RooRealVar tau2("tau2", "tau", 0.00, 0, 0.01);
  RooRealVar f1("f1", "tau", 0.5, 0, 1);

  //-------------final convolved pdf for signal--------//
  RooDecay decay1("decay1","lifetime",t,tau,gauss1,RooDecay::SingleSided); // only 1G as
  RooDecay decay2("decay2", "lifetime", t, tau2, gauss2, RooDecay::SingleSided); // only 1G as resolution function

  RooAddPdf decay("decay", "", {decay1, decay2}, {f1});
  // RooDecay decay("decay", "lifetime", t, tau, gauss1, RooDecay::SingleSided); // only 1G as

  //--------------pdf for dt signal--------//
  
  //--------------Johnson Su function -----------------//                                                                                                              
  RooRealVar xi("xi","xi",0.00003609);//0.00002,0,0.00006);
  RooRealVar lambda("lambda","lambda",0.00000938);//0.00003,0,0.0003);
  RooRealVar delta("delta","delta",1.154);//1,0,5);
  RooRealVar gamma("gamma","gamma",-2.3059);//-2,-20,2);
  RooJohnson jsu("jsu","jsu",dt,xi,lambda,gamma,delta);

  //---------gaussian---------//
  RooRealVar mean_dt("mean_dt","mean",0.0000682);//0.0001,0,0.0006);
  RooRealVar sigma_dt("sigma_dt","sigma",0.00001333);//0.00004,0,0.0008);
  RooGaussian gauss1_dt("gauss1_dt","gauss1",dt,mean_dt,sigma_dt);
  RooRealVar fr_g_dt("fr_g_dt","fr_g_dt",0.093);//0.5,0,1);
  RooAddPdf gjsu("gjsu","gjsu",RooArgList(gauss1_dt,jsu),RooArgList(fr_g_dt));


  //------final signal pdf --------------//

  RooProdPdf sig_3d("sig_3d","sig_3d",gjsu,Conditional(decay,t));

  RooProdPdf sig_3dV2("sig_3dV2", "sig_3d", gjsu, Conditional(decay1, t));
  RooProdPdf sig_3dV3("sig_3dV3", "sig_3d", gjsu, Conditional(decay2, t));


  RooFitResult *fitres = sig_3d.fitTo(*data1, Minos(0), Save(), Timer(kTRUE), Hesse(0), RooFit::Verbose(true), ConditionalObservables(RooArgSet(dt)));

  // sig_3dV3.plotOn(frame, ProjWData(RooArgSet(dt), *data1), LineColor(kBlue));
  // sig_3dV2.plotOn(frame, ProjWData(RooArgSet(dt), *data1), LineColor(kRed));
  // sig_3d.plotOn(frame, LineColor(kBlue));
  sig_3d.paramOn(frame);
// ------------

  

  // --------------
  sig_3dV2.plotOn(frame, ProjWData(RooArgSet(dt), *data1), LineColor(kBlue), LineStyle(kDashed),Normalization(data1->sumEntries() * f1.getVal(), RooAbsReal::NumEvent));
  RooCurve *curve1 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  sig_3dV3.plotOn(frame, ProjWData(RooArgSet(dt), *data1), LineColor(kGreen), LineStyle(kDashed), Normalization(data1->sumEntries() * (1 - f1.getVal()), RooAbsReal::NumEvent));
  RooCurve *curve2 = (RooCurve *)frame->getObject(frame->numItems() - 1);

  std::vector<double> xvals, yvals;
  for (int i = 0; i < curve1->GetN(); i++)
  {
    double x, y1, y2;
    curve1->GetPoint(i, x, y1);
    y2 = curve2->interpolate(x);
    xvals.push_back(x);
    yvals.push_back(y1 + y2);
  }


  RooCurve *sum_curve = new RooCurve("sum_curve", "Sum of curves", *curve1, *curve2);

  sum_curve->SetLineColor(kRed);
  sum_curve->SetLineStyle(kSolid);

  // 프레임에 추가
  frame->addObject(sum_curve);
  data1->plotOn(frame,Binning(100), Name("data"));// plotting all the data on the frame*/

  frame->Draw();
  // -------------

  //-------plotting dt --------//
  RooPlot* frame1 = dt.frame();// frame to plot lifetime t
  data1->plotOn(frame1,Binning(100));// plotting all the data on the frame*/
  sig_3d.plotOn(frame1,LineColor(kBlue));

  
  
  //--------------getting # of parameters-------------//
  RooArgSet observables(t,dt);
  RooArgSet *flparams = sig_3d.getParameters(observables);
  int nparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout<<" ndf are "<<nparam<<endl<<endl;
  

  double chi2_ndf =  frame->chiSquare(nparam);//
  cout<<"chi square = "<<chi2_ndf<<endl<<endl<<endl<<endl<<endl;
  TPaveText *box= new TPaveText(0.3, 0.85, 0.6, 0.9,"BRNDC");
  box->SetFillColor(10);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextSize(0.04F);
  box->SetFillStyle(1001);
  box->SetFillColor(10);
  TText *text = 0;
  Char_t buf[30];
  sprintf( buf,  "#chi^{2}/ndf = %f", chi2_ndf );
  text = box->AddText( buf );
  frame->addObject(box) ; 
  /////////////pull distribution////////////////
		
  //****************************** 
  RooPlot* z1frame = t.frame(Title("Pull Distribution"));
  RooHist *hpull1 = frame->pullHist("data", "sum_curve");
  hpull1->SetFillColor(4);
  hpull1->SetLineColor(0);
  z1frame->addPlotable( hpull1, "B" );
  
  
  TCanvas* c1 = new TCanvas("c1","c1",800, 800) ;
  c1->Divide(1,2) ; // column, row
  double xmin = -0.002; 
  double xmax = 0.004;
  TLine *line = new TLine(xmin,0.0,xmax,0.0);
  
  TLine *line1 = new TLine(xmin,-3.0,xmax,-3.0);
  TLine *line2 = new TLine(xmin,3.0,xmax,3.0);
  
  c1->cd(2)->SetPad(0.005,0.005,0.995,0.2525); line->SetLineColor(kRed); line->SetLineWidth(1); 
  //gPad->SetLeftMargin(0.15); 
  //z1frame->GetYaxis()->SetTitleOffset(1.45); 
  z1frame->GetYaxis()->SetTitle("Pull(#sigma)");
  z1frame->GetYaxis()->SetTitleSize(0.10);
  z1frame->GetYaxis()->SetTitleOffset(0.37);
  z1frame->GetYaxis()->SetLabelSize(0.13);
  z1frame->GetXaxis()->CenterTitle();
  z1frame->Draw();
  line1->SetLineColor(kRed); line1->SetLineWidth(1); line1->SetLineStyle(2);
  line2->SetLineColor(kRed); line2->SetLineWidth(1); line2->SetLineStyle(2);
  line->Draw("SAME"); //z1frame->GetYaxis()->SetRangeUser(-3.5, 3.5);
  line1->Draw("SAME"); 
  line2->Draw("SAME");
  c1->cd(1)->SetPad(0.005,0.2525,0.995,0.995); //gPad->SetLeftMargin(0.15); 
  //frame->GetYaxis()->SetTitleOffset(1.45);
  
  
  frame->GetXaxis()->CenterTitle();
  //	sig_3d.plotOn(frame);
  //	gaussm.plotOn(frame,LineStyle(kDashed));
  //	sig_3d.plotOn(frame,Components(RooArgList(gaussmb)),LineColor(kRed),LineStyle(7), LineWidth(2));
  frame->Draw();
	
  
  c1->Draw();
  c1->SaveAs("test_t2.pdf");

  TCanvas *cdt = new TCanvas("cdt","cdt",550,400);
  frame1->Draw();
  cdt->Draw();
  cdt->SaveAs("test_dt2.pdf");
  fitres->Print();
  //	params->printLatex() ;
  
  //   cout<<"the total no. of events are "<<nentries<<endl;
  RooPlot *frame2 = t.frame(Title("#tau D_{s}^{+}")); // frame to plot lifetime t
  data1->plotOn(frame2, Binning(100), Name("data"));  // plotting all the data on the frame*/
  // sig_3d.plotOn(frame2, LineColor(kBlue), ProjWData(RooArgSet(dt), *data1), Normalization(data1->sumEntries(), RooAbsReal::NumEvent));
  sig_3d.plotOn(frame2, LineColor(kBlue), ProjWData(RooArgSet(dt), *data1));
  TCanvas c11("c11", "", 600, 600);
  frame2->Draw();
  c11.SaveAs("projected_total_model_with_dataset.pdf");

  RooPlot *frame3 = t.frame(Title("#tau D_{s}^{+}")); // frame to plot lifetime t
  sig_3d.plotOn(frame3, LineColor(kBlue), ProjWData(RooArgSet(dt), *data1), Normalization(data1->sumEntries(), RooAbsReal::NumEvent));
  data1->plotOn(frame3, Binning(100), Name("data"));  // plotting all the data on the frame*/
  frame3->Draw();
  c11.SaveAs("projected_total_model_without_dataset.pdf");

  RooPlot *frame4 = t.frame(Title("#tau D_{s}^{+}"), Bins(100)); // frame to plot lifetime t
  sig_3d.plotOn(frame4, LineColor(kBlue), ProjWData(RooArgSet(dt), *data1), Normalization(data1->sumEntries(), RooAbsReal::NumEvent));
  // data1->plotOn(frame4, Binning(100), Name("data")); // plotting all the data on the frame*/
  frame4->addObject(sum_curve);
  frame4->Draw();
  c11.SaveAs("sum_up_vs_pdf_before_data.pdf");

  return 1;
}


