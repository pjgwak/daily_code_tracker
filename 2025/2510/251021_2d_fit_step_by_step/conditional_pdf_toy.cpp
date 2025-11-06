#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooAbsReal.h"
#include "RooProdPdf.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <string>

#include "RooCmdArg.h"

void conditional_pdf_toy(){

  // RooNumIntConfig* intConfig = RooAbsReal::defaultIntegratorConfig();
  // intConfig->setEpsRel(1.E-5);
  // intConfig->getConfigSection("RooAdaptiveIntegratorND").setRealValue("maxEval2D", 5000);

  cerr << "ping 1\n";
  RooWorkspace w("w");
  //peaks in x and y modelled as convoluted Hypatia2 functions
  w.factory("FCONV::hx(x[2258,2318],Hypatia2(x,lambda_x[-1,-10,0],zeta[0],beta_x[-1e-2,-1,1],sigmaSP_x[8,0,30],mu_x[2286,2285,2290],"\
                                              "a_x[1,0,50],n_x[10,0,50],a2_x[1,0,50],n2_x[10,0,50]),Gaussian(x,mg[0],sigmaMS_x[8,0,30]))");
  w.factory("FCONV::hy(y[1935,2005],Hypatia2(y,lambda_y[-1,-10,0],zeta,beta_y[-1e-2,-1,1],sigmaSP_y[8,0,30],mu_y[1970,1966,1975],"\
                                              "a_y[1,0,50],n_y[10,0,50],a2_y[1,0,50],n2_y[10,0,50]),Gaussian(y,mg,sigmaMS_y[8,0,30]))");
  //let the peak position of the 3rd dimension be a function of x and y pdfs
  w.factory("expr::effmu_z('@0*(@1-@2)+@3*(@4-@5)+@6',{rho_zx[1,-10,10],x,mu_x,rho_zy[1,-10,10],y,mu_y,mu_z[5620,5617,5630]})");
  w.factory("FCONV::hz(z[5550,5750],Hypatia2(z,lambda_z[-4,-25,0],zeta,beta_z[-5e-3,-1,1],sigmaSP_z[15,0,24],effmu_z,"\
                                              "a_z[2.5,0,50],n_z[10,0,50],a2_z[2.5,0,50],n2_z[10,0,50]),Gaussian(z,mg,sigmaMS_z[8,0,30]))");
  //make another pdf in z for fitting
  w.factory("FCONV::hz_simple(z,Hypatia2(z,lambda_z,zeta,beta_z,sigmaSP_z_simple[15,0,24],mu_z_simple[5620,5617,5630],"\
                             "a_z_simple[2.5,0,50],n_z_simple[10,0,50],a2_z_simple[2.5,0,50],n2_z_simple[10,0,50]),Gaussian(z,mg,sigmaMS_z_simple[8,0,30]))");

  cerr << "ping 2\n";
  //for doing some non-factory gymnastics, we need to get some things from the workspace that have been created by factory commands
  auto const& x = *(w.var("x")), y = *(w.var("y")), z = *(w.var("z"));
  auto hx = w.pdf("hx"), hy = w.pdf("hy"), hz = w.pdf("hz");
  // RooArgSet dsoifgn(x,y);
  // RooArgSet dwrg(*hz);
  // RooProdPdf model("model","",RooArgSet(*hx,*hy),RooCmdArg("Conditional",true,1,1,1,0,0,0,0,0,0,&dwrg,&dsoifgn));
  
  RooProdPdf model("model","",RooArgSet(*hx,*hy),RooFit::Conditional(*hz,RooArgSet(x,y),true)); // <-- what does depsAreCond mean?
  
  // looks ok: RooProdPdf::model[ hx * hy * hz|x,y ] = 1.04083e+09
  w.import(model);
  cerr << "ping 3\n";
  //make 3D pdfs out of it
  // w.factory("PROD::factory_model(hx,hy,hz|{x,y})"); //<-- looks fishy: RooProdPdf::factory_model[ hx * hy * y * hz|x ] = 2.05044e+12
  w.factory("PROD::factory_model(hx,hy,hz|~z)");
  cerr << "ping 4\n";
  w.factory("PROD::simple_model(hx,hy,hz_simple)");
  cerr << "ping 5\n";
  w.Print();
  cerr << "ping 6\n";
  

  //prepare plotting
  auto /*model = w.pdf("model"),*/ simple_model = w.pdf("simple_model");
  constexpr auto binning = 50, c_w /*width of canvas*/ = 1280, c_h /*height of canvas*/ = 720;
  TCanvas c("some_canvas", "", c_w, c_h);
  c.SetRightMargin(0.16);
  auto plot_and_save = [&c] (auto const hist) -> void {hist->Draw("cont4z"); c.SaveAs((static_cast<std::string>(hist->GetName())+".pdf").data()); };

  TStopwatch timer;
  auto const model_xy_hist_simple = simple_model->createHistogram("model_xy_hist_simple", x, RooFit::Binning(binning), RooFit::YVar(y, RooFit::Binning(binning)));
  std::cout << ">>>>>>>>  created hist for simple model" << std::endl;
  auto const model_xz_hist = model.createHistogram("model_xz_hist", x, RooFit::Binning(binning), RooFit::YVar(z, RooFit::Binning(binning)));
  plot_and_save(model_xz_hist);
  std::cout << ">>>>>>>>  created xz hist" << std::endl;
  auto const model_yz_hist = model.createHistogram("model_yz_hist", y, RooFit::Binning(binning), RooFit::YVar(z, RooFit::Binning(binning)));
  plot_and_save(model_yz_hist);
  std::cout << ">>>>>>>>  created yz hist" << std::endl;
  auto const model_xy_hist = model.createHistogram("model_xy_hist", x, RooFit::Binning(binning), RooFit::YVar(y, RooFit::Binning(binning)));
  plot_and_save(model_xy_hist);
  timer.Stop(); timer.Print();

  // 45 = 10min -> 

  timer.Start(true);//resets timer
  //generate data from the conditional model
  auto const data = model.generate(RooArgSet(x,y,z), 1e+2);
  //fit data with simple model and conditional model. measure time
  simple_model->fitTo(*data, RooFit::NumCPU(32));
  timer.Stop(); timer.Print();
  timer.Start(true);//resets timer
  model.fitTo(*data);
  timer.Stop(); timer.Print();
}
