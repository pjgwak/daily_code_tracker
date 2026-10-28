#ifndef OO_FIT_CONFIG_H
#define OO_FIT_CONFIG_H

#include "TFile.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TString.h"
#include <cmath>

// --- Fit policy ---

// Final parameter treatment. PDF counts and calibration values come from validated prefits.
// No nominal B value enters this policy or the likelihood.
struct OOFinalFitPolicy {
    const char *group;
    bool highPt, fullPromptCovariance, floatPromptMean, profileMass;
    int widthPolicy; // 0: all ordered widths; 1: last two; 2: last one; 3: full MC covariance
    double fractionCovarianceScale, backgroundPenaltyScale;
    double functionPrecision = -1.;
    bool calibrationSteps = false;
    static constexpr const char *minimizer = "Minuit2";
};
inline OOFinalFitPolicy oo_final_policy(bool mid, double ptLow) {
    const bool high = ptLow >= (mid ? 12. : 10.);
    // Validated common mid policy over 7 <= pT < 30; forward remains unchanged.
    if (mid) {
        OOFinalFitPolicy p{"M_7_30_tailfree_fixed", false, false, false, false, 0, 4., 1.};
        p.calibrationSteps = true;
        return p;
    }
    OOFinalFitPolicy p = ptLow<3
        ? OOFinalFitPolicy{"F_1_3_bounded_np",false,false,true,false,0,ptLow<2?16.:8.,4.}
        : ptLow<5 ? OOFinalFitPolicy{"F_3_5_joint",false,true,true,false,3,4.,1.}
        : ptLow<10 ? OOFinalFitPolicy{"F_5_10_bounded",true,false,true,false,0,4.,1.}
        : OOFinalFitPolicy{"F_10_20_np2_mean8",true,false,true,false,1,4.,1.};
    p.functionPrecision=1e-12;p.calibrationSteps=false;return p;
}
inline double oo_prompt_penalty_scale(bool mid, double ptLow, int coordinate) {
    if (mid) return 1.;
    if (ptLow == 2) return 100.;
    // One width-prior scale across the F_10_20 policy group.
    if (ptLow >= 10) return 0.1;
    return 1.;
}

inline double oo_integrated_core_strength(double pt){
 const int edges[]={5,6,7,8,9};const double factors[]={.7,.5,.6,2.4,1.8};
 for(int i=0;i<5;++i)if(pt==edges[i])return factors[i];return 0.;
}
inline double oo_integrated_width_sigma(double pt,int k,int ng){
 if(pt>=5 && pt<10){if(k==0)return .0004*oo_integrated_core_strength(pt);if(k==ng-1)return .08;if(k==1)return .0008*oo_integrated_core_strength(pt);return .004;}
 const double base[]={.010,.020,.040,.080};return base[k]*(pt<2?2.:pt<3?4.:pt<5?1.:.1);
}

// --- Initial values ---

// Numerical starts from a preliminary minimization of the same F1 likelihood.
// They do not alter MC constraint centers/covariances or fix any parameter.
struct OOInitialValue { const char *name; double value, error; };
inline const OOInitialValue *oo_f1_initial_values(unsigned &size) {
    static const OOInitialValue values[] = {
        {"Nbkg", 424081.9258028944, 746.08447836560663},
        {"Nds1", 6577.6551608981363, 3442.1357628012224},
        {"Nflip1", 33994.599241146083, 17244.023439899982},
        {"Nflip2", 13573.290525957644, 7028.0237487505638},
        {"Nplb", 173589.71694911539, 87981.466087603461},
        {"Nsig", 24421.074993430128, 395.72027801237164},
        {"NsignalSS1", 92364.287378172317, 4833.6648231800791},
        {"NsignalSS2", 64859.102807863426, 6942.6765279946521},
        {"NsignalSS3", 17868.638417832735, 3231.0176068627452},
        {"Nss1", 33599.614221887998, 17059.78136835427},
        {"Nss2", 12141.620584145354, 6270.0200094648289},
        {"Nss3", 608.77566532061553, 338.03694679205552},
        {"bFraction", 0.096057784022513248, 0.013407007280479583},
        {"bkg_flip_log_lifetime2_offset", -0.34491804947495841, 0.061151470252990836},
        {"bkg_flip_log_lifetime_offset", -2.9681290856535942, 0.14740708233614908},
        {"bkg_log_lifetime2_offset", -0.14931800817648855, 0.05963535875029552},
        {"bkg_log_lifetime3_offset", 2.9228281358889494, 0.42009305996728519},
        {"bkg_log_lifetime_offset", -5.1648077264462922, 1.4027176887248445},
        {"bkg_sym_log_lifetime_offset", 1.4697345951539826, 0.045279486222459253},
        {"ctauRecursiveFrac1", 0.36638925302200132, 0.050258120200576212},
        {"ctauRecursiveFrac2", 0.76515903610732017, 0.029837863823626076},
        {"ctauRecursiveFrac3", 0.8795913395422067, 0.023510458162060321},
        {"ctauTime1Mean", 0.0021455749876300293, 0.0025008963449622366},
        {"ctauTime1Scale", 0.081989742325220091, 0.0044997969759000347},
        {"ctauTime2Delta", 0.069563396785670667, 0.0086419062831168299},
        {"ctauTime3Delta", 0.15940150322166821, 0.025572759445899138},
        {"ctauTime4Delta", 0.52634093539988225, 0.078407479609111624},
        {"signal2_log_lifetime_offset", -0.34431323278073966, 0.027622712163227009},
        {"signal3_log_lifetime_offset", -2.3081118397811458, 0.10563541864837167},
        {"signal_log_lifetime_offset", -1.0575469913698159, 0.058049262184450834},
    };
    size = sizeof(values)/sizeof(values[0]);
    return values;
}

// --- Policy metadata ---

inline double oo_prompt_fraction_error_scale(const OOFinalFitPolicy &p, bool mid, double pt) {
    return p.fullPromptCovariance ? 10. : (!mid && pt < 3 ? std::sqrt(p.fractionCovarianceScale) : 0.);
}
inline void write_oo_forward_policy_metadata(TFile &file, const OOFinalFitPolicy &p, bool mid, double pt,
                                    int ng, int ss, int ds, int flip, int bgng) {
    TDirectory::TContext context(&file);
    auto d=[](const char *name,double v){TParameter<double>(name,v).Write("",TObject::kOverwrite);};
    auto i=[](const char *name,int v){TParameter<int>(name,v).Write("",TObject::kOverwrite);};
    d("promptFractionMCErrorScale",oo_prompt_fraction_error_scale(p,mid,pt));
    d("promptJointMCCovarianceScale",p.fullPromptCovariance?100.:0.);
    d("promptFractionMCCovarianceScale",!p.fullPromptCovariance&&!mid&&pt<5?p.fractionCovarianceScale:0.);
    d("npJointMCCovarianceScale",pt<2?175.:pt<3?700.:pt<10?1.:2.);
    d("backgroundLifetimePenaltyMultiplier",p.highPt?0.:p.backgroundPenaltyScale);
    TNamed("configuredMinimizer",p.minimizer).Write("",TObject::kOverwrite);
    d("configuredFunctionPrecision",p.functionPrecision);
    d("configuredMinimizerTolerance",1e-3);
    const double base[]={.010,.020,.040,.080};double broadest=0;
    for(int k=0;k<4;++k){
        const bool active=k<ng&&!p.fullPromptCovariance&&!(p.widthPolicy==1&&k<ng-2)&&!(p.widthPolicy==2&&k!=ng-1);
        const double sigma=active?oo_integrated_width_sigma(pt,k,ng):0.;
        d(Form("promptWidthPriorSigma%d",k+1),sigma);
        if(k==ng-1)broadest=sigma;
    }
    d("promptBroadestAbsolutePenalty",broadest);
    i("massYieldConstraintCount",0);
    i("backgroundShapeProfiled",p.highPt?1:0);
    i("activeBackgroundLifetimeConstraintCount",p.highPt?0:ss+ds+flip);
    i("nBkgSSComponents",ss);i("nBkgDSComponents",ds);i("nBkgFlipComponents",flip);i("nBkgResolutionComponents",bgng);
    i("fitPolicyMetadataVersion",1);
}

inline void write_oo_policy_metadata(TFile &file, const OOFinalFitPolicy &p, bool mid, double pt,
                                    int ng, int ss, int ds, int flip, int bgng) {
    if (!mid) { write_oo_forward_policy_metadata(file,p,mid,pt,ng,ss,ds,flip,bgng); return; }
    TDirectory::TContext context(&file);
    auto d=[](const char *name,double v){TParameter<double>(name,v).Write("",TObject::kOverwrite);};
    auto i=[](const char *name,int v){TParameter<int>(name,v).Write("",TObject::kOverwrite);};
    d("promptFractionMCErrorScale",oo_prompt_fraction_error_scale(p,mid,pt));
    d("promptJointMCCovarianceScale",p.fullPromptCovariance?100.:0.);
    d("promptFractionMCCovarianceScale",!p.fullPromptCovariance&&!mid&&pt<3?p.fractionCovarianceScale:0.);
    d("npJointMCCovarianceScale",p.highPt?0.:(pt==6?100.:1.));
    d("backgroundLifetimePenaltyMultiplier",(mid||p.highPt)?0.:p.backgroundPenaltyScale);
    TNamed("configuredMinimizer",p.minimizer).Write("",TObject::kOverwrite);
    d("configuredFunctionPrecision",p.functionPrecision);
    d("configuredMinimizerTolerance",p.functionPrecision>0?1e-1:1e-3);
    const double base[]={.010,.020,.040,.080};double broadest=0;
    for(int k=0;k<4;++k){
        const bool active=k<ng&&!p.fullPromptCovariance&&!(p.widthPolicy==1&&k<ng-2)&&!(p.widthPolicy==2&&k!=ng-1);
        const double sigma=active?base[k]*oo_prompt_penalty_scale(mid,pt,k):0.;
        d(Form("promptWidthPriorSigma%d",k+1),sigma);
        if(k==ng-1)broadest=sigma;
    }
    d("promptBroadestAbsolutePenalty",broadest);
    i("massYieldConstraintCount",0);
    i("backgroundShapeProfiled",(mid||p.highPt)?1:0);
    i("activeBackgroundLifetimeConstraintCount",(mid||p.highPt)?0:ss+ds+flip);
    i("nBkgSSComponents",ss);i("nBkgDSComponents",ds);i("nBkgFlipComponents",flip);i("nBkgResolutionComponents",bgng);
    i("fitPolicyMetadataVersion",1);
}

#endif
