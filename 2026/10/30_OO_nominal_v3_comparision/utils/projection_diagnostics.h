#ifndef OO_PROJECTION_DIAGNOSTICS_H
#define OO_PROJECTION_DIAGNOSTICS_H
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

// Independent of minimization and parameter policy. Expected counts are PDF bin
// integrals with the FULL fit normalization, including in the restricted views.
// Projection NDF is descriptive: constraints and projections do not in general
// have an exact chi-square reference distribution. No p-value is inferred here.
namespace oo_diagnostics {
struct Component { RooAbsPdf *pdf; double yield; const char *label; int color; };
inline double integral(RooAbsPdf &pdf, RooRealVar &x, double lo, double hi) {
    x.setRange("ooDiagnosticBin", lo, hi);
    RooArgSet norm(x);
    std::unique_ptr<RooAbsReal> a(pdf.createIntegral(norm, RooFit::NormSet(norm), RooFit::Range("ooDiagnosticBin")));
    const double v=a->getVal();
    if (!std::isfinite(v) || v < -1e-9) throw std::runtime_error("Invalid diagnostic PDF integral");
    return std::max(0.,v);
}
inline double quantile(RooAbsPdf &pdf, RooRealVar &x, double q) {
    double lo=x.getMin(),hi=x.getMax();
    for(int k=0;k<42;++k){ double m=(lo+hi)/2;
        if(integral(pdf,x,x.getMin(),m)<q) lo=m; else hi=m; }
    return (lo+hi)/2;
}
// Statistics-based rendering policy, independent of goodness-of-fit values.
inline int display_bins(RooAbsData &data, RooRealVar &x, double lo, double hi) {
    const double original=x.getVal();
    double sumw=0, sumw2=0;
    for(int k=0;k<data.numEntries();++k){const auto *row=data.get(k);
        const double value=row->getRealValue(x.GetName());
        if(value>=lo && value<hi){const double w=data.weight();sumw+=w;sumw2+=w*w;}}
    x.setVal(original);
    const double effective=sumw2>0?sumw*sumw/sumw2:0;
    return std::max(400,50*int(std::ceil(2*std::sqrt(effective)/50)));
}
inline void run(RooAbsPdf &pdf, RooAbsData &data, RooRealVar &x, int npar,
                int nativeBins, const TString &tag, bool lifetime,
                std::vector<Component> components={}) {
    const char *dir=gSystem->Getenv("OO_DIAGNOSTICS_DIR");
    if(!dir || !*dir) return;
    gSystem->mkdir(dir,true);
    const double original=x.getVal();
    const double total=pdf.canBeExtended()?pdf.expectedEvents(RooArgSet(x)):data.sumEntries();
    if(!(total>0)) throw std::runtime_error("Invalid projection normalization");
    const double fullLo=x.getMin(),fullHi=x.getMax();
    RooAbsPdf &core=components.empty()?pdf:*components.front().pdf;
    double focusLo=fullLo,focusHi=fullHi,coreLo=fullLo,coreHi=fullHi;
    if(lifetime){
        coreLo=quantile(core,x,.005); coreHi=quantile(core,x,.995);
        focusLo=coreLo; focusHi=coreHi;
        // Only PR/NP determine the signal window; background is displayed but
        // deliberately does not expand it back to the background's long tails.
        for(size_t c=0;c<std::min(size_t(2),components.size());++c){
            focusLo=std::min(focusLo,quantile(*components[c].pdf,x,.005));
            focusHi=std::max(focusHi,quantile(*components[c].pdf,x,.995));
        }
    }
    const int focusBins=display_bins(data,x,focusLo,focusHi), coreBins=display_bins(data,x,coreLo,coreHi);
    struct View {const char *name; double lo,hi; int bins;};
    std::vector<View> views{{"native",fullLo,fullHi,nativeBins}};
    for(int bins:{60,120,240}) views.push_back({"full",fullLo,fullHi,bins});
    if(lifetime){
        for(int bins:{60,120,240,480}) views.push_back({"signal",focusLo,focusHi,bins});
        views.push_back({"core",coreLo,coreHi,120});
    }
    views.push_back({"display_full",fullLo,fullHi,std::max(nativeBins,display_bins(data,x,fullLo,fullHi))});
    if(lifetime){
        views.push_back({"signal_display",focusLo,focusHi,focusBins});
        views.push_back({"signal_dense",focusLo,focusHi,2*focusBins});
        if(components.size()>1)views.push_back({"core_display",coreLo,coreHi,coreBins});
    }
    const TString stem=TString(dir)+"/"+tag;
    std::ofstream csv(stem+".csv");
    csv << std::setprecision(17)
        << "view,low,high,bins,width,npar,nonzero,roo_ndf,roo_chi2_ndf,pearson,deviance,all_ndf,expected_lt5,grouped_bins,grouped_ndf,grouped_pearson,grouped_deviance,observed,expected,pr_retained,np_retained,weighted,max_center_bias_pull,max_abs_pull,pull_overflow\n";
    TCanvas canvas("ooDiagnosticCanvas","Integrated projection diagnostics",1500,850);
    const int columns=lifetime?(components.size()>1?3:2):1;
    canvas.Divide(columns,2,.001,.001);
    std::vector<std::unique_ptr<TH1D>> keep;
    int panel=0;
    for(const auto &v:views){
        std::unique_ptr<TH1D> h(new TH1D("ooObserved","",v.bins,v.lo,v.hi)); h->SetDirectory(nullptr);h->Sumw2();
        for(int k=0;k<data.numEntries();++k){const auto *row=data.get(k); double value=row->getRealValue(x.GetName());
            if(value>=v.lo && value<v.hi) h->Fill(value,data.weight());}
        auto expected=std::make_unique<TH1D>("ooExpected","",v.bins,v.lo,v.hi);expected->SetDirectory(nullptr);
        auto pull=std::make_unique<TH1D>("ooResidual","",v.bins,v.lo,v.hi);pull->SetDirectory(nullptr);
        double pearson=0,deviance=0,observed=0,predicted=0,maxCenterBias=0,maxAbsPull=0;int nz=0,sparse=0,pullOverflow=0;
        std::vector<std::pair<double,double>> groups; double gn=0,gm=0;
        for(int b=1;b<=v.bins;++b){
            double mu=total*integral(pdf,x,h->GetBinLowEdge(b),h->GetBinLowEdge(b+1));
            double n=h->GetBinContent(b); expected->SetBinContent(b,mu);
            if(!(mu>0)) throw std::runtime_error("Nonpositive expected diagnostic bin");
            x.setVal(h->GetBinCenter(b));
            const RooArgSet norm(x);
            const double center=total*pdf.getVal(&norm)*h->GetBinWidth(b);
            maxCenterBias=std::max(maxCenterBias,std::abs(center-mu)/std::sqrt(mu));
            observed+=n;predicted+=mu;if(n>0)++nz;if(mu<5)++sparse;
            pearson+=(n-mu)*(n-mu)/mu;deviance+=2*(mu-n+(n>0?n*std::log(n/mu):0));
            const double residual=(n-mu)/std::sqrt(mu);
            pull->SetBinContent(b,residual);maxAbsPull=std::max(maxAbsPull,std::abs(residual));
            if(std::abs(residual)>8)++pullOverflow;
            gn+=n;gm+=mu;
            if(gm>=5){groups.push_back({gn,gm});gn=gm=0;}
        }
        if(gm>0){if(groups.empty())groups.push_back({gn,gm});else{groups.back().first+=gn;groups.back().second+=gm;}}
        double gp=0,gd=0;for(auto g:groups){gp+=(g.first-g.second)*(g.first-g.second)/g.second;
            gd+=2*(g.second-g.first+(g.first>0?g.first*std::log(g.first/g.second):0));}
        std::unique_ptr<RooPlot> frame(x.frame(RooFit::Range(v.lo,v.hi),RooFit::Bins(v.bins)));
        data.plotOn(frame.get(),RooFit::Binning(v.bins,v.lo,v.hi),RooFit::Name("data"),
                    RooFit::DataError(data.isWeighted()?RooAbsData::SumW2:RooAbsData::Poisson));
        x.setRange("ooDiagnosticFull",fullLo,fullHi);
        pdf.plotOn(frame.get(),RooFit::Range(v.lo,v.hi),RooFit::NormRange("ooDiagnosticFull"),
                   RooFit::Normalization(total,RooAbsReal::NumEvent),RooFit::Precision(1e-6),RooFit::Name("model"));
        const double nan=std::numeric_limits<double>::quiet_NaN();
        csv<<v.name<<','<<v.lo<<','<<v.hi<<','<<v.bins<<','<<(v.hi-v.lo)/v.bins<<','<<npar<<','<<nz<<','<<nz-npar<<','
           <<frame->chiSquare("model","data",npar)<<','<<(data.isWeighted()?nan:pearson)<<','<<(data.isWeighted()?nan:deviance)<<','<<v.bins-npar<<','<<sparse<<','
           <<groups.size()<<','<<int(groups.size())-npar<<','<<(data.isWeighted()?nan:gp)<<','<<(data.isWeighted()?nan:gd)<<','<<observed<<','<<predicted<<','
           <<integral(core,x,v.lo,v.hi)<<','<<(components.size()>1?integral(*components[1].pdf,x,v.lo,v.hi):nan)<<','<<data.isWeighted()<<','<<maxCenterBias<<','<<maxAbsPull<<','<<pullOverflow<<'\n';
        const std::string name=v.name;
        if(name!="display_full" && name!="signal_display" && name!="core_display") continue;
        ++panel;canvas.cd(panel);gPad->SetLogy(lifetime);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.24);
        h->SetTitle(Form("%s: %s; %s;Events / %.4g %s",tag.Data(),v.name,x.GetName(),(v.hi-v.lo)/v.bins,x.getUnit()));
        h->SetStats(false);h->SetMarkerStyle(20);h->SetMarkerSize(.35);h->SetMinimum(lifetime?.2:0);h->SetMaximum(h->GetMaximum()*(lifetime?4:1.3));h->Draw("E");
        expected->SetLineColor(kBlack);expected->SetLineWidth(2);expected->Draw("HIST SAME");
        for(auto c:components){auto ch=std::make_unique<TH1D>("ooComponent","",v.bins,v.lo,v.hi);ch->SetDirectory(nullptr);
            for(int b=1;b<=v.bins;++b)ch->SetBinContent(b,c.yield*integral(*c.pdf,x,ch->GetBinLowEdge(b),ch->GetBinLowEdge(b+1)));
            ch->SetLineColor(c.color);ch->SetLineStyle(2);ch->Draw("HIST SAME");keep.push_back(std::move(ch));}
        TLatex text;text.SetNDC();text.SetTextSize(.028);
        text.SetTextSize(.021);text.DrawLatex(.16,.96,tag);text.SetTextSize(.028);
        text.DrawLatex(.19,.92,Form("%s: %d bins",v.name,v.bins));
        text.DrawLatex(.19,.88,Form("Roo #chi^{2}/ndf = %.3f (%d ndf)",frame->chiSquare("model","data",npar),nz-npar));
        text.DrawLatex(.19,.84,Form("Grouped Pearson/ndf = %.3f (%d ndf)",gp/(int(groups.size())-npar),int(groups.size())-npar));
        if(components.size()>1)text.DrawLatex(.19,.80,Form("PR %.2f%%, NP %.2f%% retained; red PR, blue NP",100*integral(core,x,v.lo,v.hi),100*integral(*components[1].pdf,x,v.lo,v.hi)));
        canvas.cd(panel+columns);gPad->SetBottomMargin(.15);pull->SetTitle(Form("Integrated expected counts; %s; (data-model)/#sqrt{model}",x.GetName()));
        pull->SetStats(false);pull->SetMinimum(-8);pull->SetMaximum(8);pull->Draw("HIST");
        if(pullOverflow){text.DrawLatex(.19,.92,Form("%d pulls outside #pm8; max |pull| = %.1f",pullOverflow,maxAbsPull));}
        keep.push_back(std::move(h));keep.push_back(std::move(expected));keep.push_back(std::move(pull));
    }
    canvas.cd();
    canvas.Print(stem+".pdf");
    x.setVal(original);
    std::cout<<"[DIAGNOSTICS] "<<stem<<std::endl;
}
}
#endif
