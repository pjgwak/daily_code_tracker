//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "TH1.h"
#include "TH2.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TLine.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
// #include "RooNLLVar.h"
#include "TFile.h"
#include "RooFit.h"
#include "RooGaussModel.h"
#include "RooHistPdf.h"
#include "RooPolyVar.h"

#include "myRooGaussian.cpp"

using namespace RooFit ;
using namespace std;



void unbinned_Kpipi0_RS_2D_fit()  
{
	
  // get the data in a TTree format
  TFile *f = new TFile("DtoKpipi0_RS_mc14ri_a_2Dfit.root");
  TTree *tree = (TTree*)f->Get("DstD0PiKPiPi0RS_withcut");
  
  // check we have the tree
  cout << "Entries = " << tree->GetEntries() << endl;

  //create both fit observables and get the data in RooFit land
  RooRealVar m("pi0_m_prefit","#it{m}(#it{#pi}^{0})", 0.1, 0.18, "GeV/#it{c}^{2}");
  RooRealVar errm("pi0_errm","#it{m}(#it{#pi}^{0}) error", 0, 0.015, "GeV/#it{c}^{2}");
  RooDataSet* dataxy = new RooDataSet("data","data",RooArgSet(m, errm), Import(*tree));

  
  //*  S i n g l e    G a u s s i o n     R e s o l u t i o n    f u n c t i o n
  //  ------------------------------------------------------------------------

  //Create function f(y) = a0 + a1*y (consider the variation of the mean as a function of mass error)
  RooRealVar p0("p0", "p0", 0.131985, 0, 10); //0.131985  
  RooRealVar p1("p1", "p1", 0.346866, 0, 10); //0.346866
  RooPolyVar mu("mu", "mean", errm, RooArgSet(p0, p1));

  //Create function f(y) = 1 + a/x + b/x^2 
  RooRealVar a("a", "a", -5.71988e-03, -1e-01, 1e-01); //-2.71988e-03  
  RooRealVar b("b", "b", 4.42226e-05, -1e-03, 1e-03);  //1.42226e-05
  // RooFormulaVar s("s", "1.0 + (@1/@0) + (@2/(TMath::Power(@0, 2.0)))", RooArgList(errm, a, b));

  // RooFormulaVar sigma("sigma", "s*pi0_errm", RooArgList(s,errm));

  RooFormulaVar sigma("sigma", "@0+@1+@2/@0", RooArgList(errm,a,b));

  myRooGaussian gauss("gauss", "1st gaussian PDF", m, mu, sigma);
  //RooGaussian gauss("gauss", "1st gaussian PDF", m, mu, sigma);


  // C o n s t r u c t  p d f   f o r   e r r o r   m a s s
  // -----------------------------------------------------------------
  
  // Construct a histogram pdf to describe the shape of the ErrM distribution
  RooDataSet datamerr("datamerr","datamerr", RooArgSet(errm), Import(*tree));

  RooDataHist *expHistDterr = datamerr.binnedClone();
  RooHistPdf pdfErr("pdfErr", "pdfErr", errm, *expHistDterr);
 

  // C o n s t r u c t   c o n d i t i o n a l   p r o d u c t   d e c a y 
   // -----------------------------------------------------------------------------
  
  // Construct production of conditional decay_dm(dt|dterr) with empirical pdfErr(dterr)
  RooProdPdf pdf("pdf", "pdf(m|sigma_m)*pdf(sigma_m)", pdfErr, RooFit::Conditional(gauss, m));

  m.setRange("SB1",0.1,0.18) ;
  errm.setRange("SB2",0,0.0031) ;
  
  //perform the fit
  pdf.fitTo(*dataxy, RooFit::Save(true), RooFit::Strategy(2));
  
 /*-----------------------------------PLOT THE RESULTS-------------------------------------*/

  //-----------------------Canvas and pad
  TCanvas *canvas = new TCanvas("canvas","canvas", 700, 800);  //  700, 800
 
  Double_t xlow, ylow, xup, yup;
  canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);
  canvas->Divide(1,2);

  TVirtualPad *upPad = canvas->GetPad(1);
  upPad->SetPad(xlow,ylow+0.25*(yup-ylow),xup,yup);
  
  TVirtualPad *dwPad = canvas->GetPad(2);
  dwPad->SetPad(xlow,ylow,xup,ylow+0.25*(yup-ylow));


   
  //------------------------Plot pi0_M_prefit projections
  RooPlot *xframe = m.frame();
  dataxy->plotOn(xframe, CutRange("SB2"));
  const double nData = dataxy->sumEntries("", "SB2");
  cout<<" nentries in data are = "<<nData<<endl;
  // pdf.plotOn(xframe, ProjWData(*dataxy), ProjectionRange("SB2"));
  //pdf.plotOn(xframe);
  //pdf.paramOn(xframe);
  canvas->cd(1);
  pdf.plotOn(xframe,Components(gauss),LineStyle(kDashed),LineColor(kRed), ProjWData(RooArgList(errm), *dataxy), Normalization(nData, RooAbsReal::NumEvent));
  xframe->Draw();



   //----------------------- pull distribution of pi0_M_prefit
  RooHist* hpull = xframe->pullHist();
  hpull->SetFillStyle(1001);
  hpull->SetFillColor(1);
  for(int i=0;i<hpull->GetN();++i) hpull->SetPointError(i,0.0,0.0,0.0,0.0);
  RooPlot* pullplot = m.frame(Title(" "));
  pullplot->addPlotable(hpull,"B");
  pullplot->SetYTitle("Pull mass");
  pullplot->SetMinimum(-4.);
  pullplot->SetMaximum(4.);
  pullplot->GetXaxis()->SetLabelSize(0.1);
  pullplot->GetYaxis()->SetLabelSize(0.07);
  canvas->cd(2);
  pullplot->Draw();


  //------------------------- adding reference line in pull distribution
  double xmin1 = 0.1; 
  double xmax1 = 0.18;
  TLine *line = new TLine(xmin1,0.0,xmax1,0.0);
  TLine *line1 = new TLine(xmin1,3.0,xmax1,3.0);
  TLine *line2 = new TLine(xmin1,-3.0,xmax1,-3.0);
 
  line->SetLineColor(kRed); 
  line->SetLineWidth(3); 
  line1->SetLineColor(kRed);
  line2->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line2->SetLineStyle(2);
  line->Draw("SAME"); 
  line1->Draw("SAME");
  line2->Draw("SAME");
  
  

  //------------------------- Plot pi0_ErrM projections
  RooPlot *yframe = errm.frame();
  dataxy->plotOn(yframe);
  pdf.plotOn(yframe);

  canvas->Update();
  canvas->SaveAs("my_test.png");

  // TCanvas *c = new TCanvas("c", "c", 700, 800);
  // //canvas->cd(2);
  // gPad->SetLeftMargin(0.15);
  // yframe->GetYaxis()->SetTitleOffset(1.6);
  // yframe->Draw();
 

 
}

  //  RooPlot *xframe = m.frame();
  // dataxy->plotOn(xframe, CutRange("SB2"));
  // const double nData = dataxy->sumEntries("", "SB2");
  // cout<<" nentries in data are = "<<nData<<endl;
  // pdf.plotOn(xframe, ProjWData(*dataxy), CutRange("SB2"));
  // //pdf.plotOn(xframe);
  // //pdf.paramOn(xframe);
  // canvas->cd(1);
  // pdf.plotOn(xframe,Components(gauss),LineStyle(kDashed),LineColor(kRed), ProjWData(*dataxy), Normalization(nData, RooAbsReal::NumEvent), ProjectionRange("SB2"));
  // xframe->Draw();


 

  //RooRealVar constant("constant", "",1);
  
  // combine the parameters, as they enter in different combinations in the fit pdf
  // RooFormulaVar a_inv("a_inv", "@0/@1", RooArgList(a, pi0_ErrM));
  // RooFormulaVar b_inv("b_inv", "@0/(TMath::Power(@1, 2.0))", RooArgList(b, pi0_ErrM));
  // RooFormulaVar s("s", "@0+@1+@2", RooArgList(constant, a_inv, b_inv));



  //Double Gaussisn
  // RooRealVar a("a", "a", 2.71150e+00, -5, 5); //2.71150e+00
  // RooRealVar b("b", "b", -4.98568e+02, -1e+04, 1e+04);  //-4.98568e+02
  // RooRealVar c("c", "c", 3.30340e+04, -1e+06, 1e+06);  //-4.98568e+02
  // RooPolyVar s("s", "s", errm, RooArgSet(a, b, c));




 // // create the different fit parameters to be used in the fit pdf
 //  RooRealVar f_G1("f_{G1}", "frac_gauss1", 0.2, 0, 1);  //0.233
 //  RooRealVar mpi0("mpi0", "mpi0", 0.1349768); // PDG pi0 mass
 //  RooRealVar b1("b1", "b1", 0.0092, -5, 5);  
 //  RooRealVar b2("b2", "b2", 0.0081, -5, 5); 
 //  // RooRealVar s1("s1", "s1", 2.327, 0, 10); 
 //  // RooRealVar s2("s2", "s2", 0.970, 0, 10);  
 
 //  //Create function f(y) = p0 + p1*y + p2*y^2 
 //  RooRealVar p0_1("p0_1", "p0_1", 2.71150);
 //  RooRealVar p1_1("p1_1", "p1_1", -4.98568e+02);
 //  RooRealVar p2_1("p2_1", "p2_1", 3.30340e+04);
 //  RooPolyVar s1("s1", "s1", pi0_ErrM, RooArgSet(p0_1, p1_1, p2_1));

 //  RooRealVar p0_2("p0_2", "p0_2", 2.71150);
 //  RooRealVar p1_2("p1_2", "p1_2", -4.98568e+02);
 //  RooRealVar p2_2("p2_2", "p2_2", 3.30340e+04);
 //  RooPolyVar s2("s2", "s2", pi0_ErrM, RooArgSet(p0_2, p1_2, p2_2));

  
 //  //combine the parameters, as they enter in different combinations in the fit pdf
 //  RooFormulaVar mu_G1("mu_G1", "@0+@1", RooArgList(mpi0,b1)); //@mpi0+@b1
 //  RooFormulaVar sigma_G1("sigma_G1", "s1*pi0_errm", RooArgList(s1,pi0_ErrM));
  
 //  RooFormulaVar mu_G2("mu_G2", "@0+@1", RooArgList(mpi0,b2));
 //  RooFormulaVar sigma_G2("sigma_G2", "s2*pi0_errm", RooArgList(s2,pi0_ErrM));
  
 //  // Resolution function (Double Gaussion)
 //  // myRooGaussian gauss1("gauss1", "1st gaussian PDF", pi0_M_prefit, mu_G1, sigma_G1);
 //  // myRooGaussian gauss2("gauss2", "2nd gaussian PDF", pi0_M_prefit, mu_G2, sigma_G2);
 //  RooGaussModel gauss1("gauss1", "1st gaussian PDF", pi0_M_prefit, mu_G1, sigma_G1);
 //  RooGaussModel gauss2("gauss2", "2nd gaussian PDF", pi0_M_prefit, mu_G2, sigma_G2);
 //  RooAddPdf doublegaussian("doublegaussian", "doublegaussian_sig", gauss1, gauss2, f_G1);
