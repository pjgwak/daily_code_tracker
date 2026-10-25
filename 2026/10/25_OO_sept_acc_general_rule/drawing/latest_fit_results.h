#ifndef OO_DRAWING_LATEST_FIT_RESULTS_H
#define OO_DRAWING_LATEST_FIT_RESULTS_H
#include "../utils/nominal_bfraction.h"
#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace oo_drawing {
struct FitValues { double b, bError, yield, yieldError; };
inline FitValues readFit(const TString &path) {
    TFile file(path, "READ");
    auto *result=file.Get<RooFitResult>("modelResult");
    auto *backend=file.Get<TNamed>("configuredMinimizer");
    auto *migrad=file.Get<TParameter<int>>("minimizeStatus");
    auto *hesse=file.Get<TParameter<int>>("hesseStatus");
    if(file.IsZombie() || !result || result->status()!=0 || result->covQual()!=3 ||
       !std::isfinite(result->edm()) || !backend || std::string(backend->GetTitle())!="Minuit2" ||
       !migrad || migrad->GetVal()!=0 || !hesse || hesse->GetVal()!=0)
        throw std::runtime_error("Invalid or incompatible final fit: "+std::string(path.Data()));
    auto *b=dynamic_cast<const RooRealVar*>(result->floatParsFinal().find("bFraction"));
    auto *n=dynamic_cast<const RooRealVar*>(result->floatParsFinal().find("Nsig"));
    if(!b || !n || !std::isfinite(b->getVal()) || !std::isfinite(b->getError()) ||
       !std::isfinite(n->getVal()) || !std::isfinite(n->getError()) || b->getVal()<0 || b->getVal()>1 ||
       b->getError()<=0 || n->getVal()<=0 || n->getError()<=0)
        throw std::runtime_error("Missing or invalid B/Nsig: "+std::string(path.Data()));
    return {b->getVal(),b->getError(),n->getVal(),n->getError()};
}
template<int N> void loadSeries(const TString &roots, bool mid, const double (&low)[N],
    const double (&high)[N], double (&values)[N], double (&errors)[N], bool yield, const TString &csv) {
    std::ofstream out(csv.Data());
    if(!out)throw std::runtime_error("Cannot write plot inputs: "+std::string(csv.Data()));
    out<<std::setprecision(17)<<"pt_low,pt_high,value,error,observable,minimizer,source\n";
    for(int i=0;i<N;++i) {
        const TString tag=Form("y%.2f_%.2f",mid?0.:1.6,mid?1.6:2.4);
        const TString path=Form("%s/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root",roots.Data(),tag.Data(),tag.Data(),low[i],high[i]);
        const auto fit=readFit(path);
        values[i]=yield?fit.yield:fit.b;errors[i]=yield?fit.yieldError:fit.bError;
        out<<low[i]<<','<<high[i]<<','<<values[i]<<','<<errors[i]<<','<<(yield?"Nsig":"bFraction")<<','<<"Minuit2"<<','<<path<<'\n';
    }
}
template<int N> void loadNominal(bool mid, const double (&low)[N], const double (&high)[N],
    double (&values)[N], double (&errors)[N]) {
    for(int i=0;i<N;++i) {
        const auto &ref=nominal_bfraction::find(low[i],high[i],mid?0.:1.6,mid?1.6:2.4);
        values[i]=ref.value;errors[i]=ref.error;
    }
}
}
#endif
