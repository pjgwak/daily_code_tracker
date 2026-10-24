#include "projection_diagnostics.h"
#include "saved_fit_helpers.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "TMath.h"
#include <sstream>
#include <cassert>
void test_projection_diagnostics() {
    RooRealVar x("x","x",-10.,10.), mean("mean","mean",0.),sigma("sigma","sigma",.3),yield("yield","yield",1000.);
    RooGaussian g("g","g",x,mean,sigma);
    RooExtendPdf p("p","p",g,yield);
    const double exact=.5*(std::erf(.7/(.3*std::sqrt(2.)))-std::erf(-.2/(.3*std::sqrt(2.))));
    if(std::abs(oo_diagnostics::integral(p,x,-.2,.7)-exact)>1e-10)throw std::runtime_error("Gaussian bin integral mismatch");
    const double q=oo_diagnostics::quantile(p,x,.995);
    if(std::abs(q-.3*2.5758293035489)>1e-9)throw std::runtime_error("Gaussian quantile mismatch");
    // A deterministic normal sample with empty tails exercises the distinct
    // RooPlot nonzero-bin and Poisson all-bin denominator conventions.
    RooDataSet data("data","data",RooArgSet(x));
    for(int k=0;k<1000;++k){x.setVal(.3*std::sqrt(2.)*TMath::ErfInverse(2*(k+.5)/1000.-1));data.add(RooArgSet(x));}
    RooDataSet large("large","large",RooArgSet(x));
    x.setVal(0.);
    for(int k=0;k<100000;++k)large.add(RooArgSet(x));
    if(oo_diagnostics::display_bins(large,x,-1.,1.)!=650 || oo_diagnostics::display_bins(data,x,-1.,1.)!=400)
        throw std::runtime_error("Display resolution does not grow with statistics");
    x.setVal(.123);
    gSystem->Setenv("OO_DIAGNOSTICS_DIR","drawing/projection_diagnostics_tests");
    oo_diagnostics::run(p,data,x,2,100,"gaussian",true);
    if(x.getVal()!=.123)throw std::runtime_error("Diagnostic changed observable state");
    std::ifstream in("drawing/projection_diagnostics_tests/gaussian.csv");std::string line;std::getline(in,line);int rows=0;
    while(std::getline(in,line)){
        std::vector<std::string> f;std::stringstream stream(line);std::string field;
        while(std::getline(stream,field,','))f.push_back(field);
        if(f.size()!=25)throw std::runtime_error("Incomplete diagnostic CSV");
        const double retained=std::stod(f[19]);
        if(std::abs(std::stod(f[18])-1000*retained)>1e-6)throw std::runtime_error("Restricted view incorrectly renormalized");
        if(f[0]=="native" && (std::stoi(f[6])>=100 || std::stoi(f[11])!=98))throw std::runtime_error("Empty-bin denominator mismatch");
        if(f[0]=="signal" && std::abs(retained-.99)>1e-9)throw std::runtime_error("Signal coverage mismatch");
        ++rows;
    }
    if(rows!=12)throw std::runtime_error("Missing binning scenarios");
    gSystem->Unsetenv("OO_DIAGNOSTICS_DIR");
    std::cout<<"PASS: analytic integrals, 99% quantiles, restricted normalization, empty bins, statistics-dependent resolution, all twelve views, state restoration"<<std::endl;
}
