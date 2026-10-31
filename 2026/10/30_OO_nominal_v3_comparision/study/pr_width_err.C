// Dump MC PR prefit (ctau resolution) parameter errors per bin: config/pr_width_err.txt
#include <cstdio>
void pr_width_err() {
    const char* names[] = {"ctauTime1Scale","ctauTime2Delta","ctauTime3Delta","ctauTime4Delta","ctauFrac1Var","ctauFrac2Var","ctauFrac3Var","ctauTime1Mean"};
    FILE* out = fopen("config/pr_width_err.txt","w");
    fprintf(out,"bin\tparameter\tvalue\terror\n");
    struct R { const char* reg; const char* y; std::vector<double> e; };
    std::vector<R> regs = {{"fwd","y1.60_2.40",{1,2,3,4,5,6,7,8,9,10,12,14,20}},{"mid","y0.00_1.60",{7,8,9,10,12,14,16,20,30}}};
    for (auto& r : regs) for (size_t i=0;i+1<r.e.size();++i) {
        TString tag = Form("%s_pt%.2f_%.2f", r.y, r.e[i], r.e[i+1]);
        TFile f(Form("roots/%s/ctau_pr/ctau_resolution_%s.root", r.y, tag.Data()));
        auto* res = (RooFitResult*)f.Get("timeResult"); if(!res){printf("no result %s\n",tag.Data());continue;}
        for (auto n : names) { auto* p=(RooRealVar*)res->floatParsFinal().find(n); if(p) fprintf(out,"%s_%g_%g\t%s\t%.8g\t%.8g\n",r.reg,r.e[i],r.e[i+1],n,p->getVal(),p->getError()); }
    }
    fclose(out);
}
