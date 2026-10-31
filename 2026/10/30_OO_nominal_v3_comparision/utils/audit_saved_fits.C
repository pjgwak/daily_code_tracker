#include "nominal_bfraction.h"
#include "TFile.h"
#include "TParameter.h"
#include "TNamed.h"
#include "RooFitResult.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <limits>

template<class T> T audit_value(TFile &file, const char *name, T fallback) {
 auto *p=file.Get<TParameter<T>>(name); return p ? p->GetVal() : fallback;
}
void audit_saved_fits(const char *rootDir="roots", const char *csv="fit_audit.csv") {
 std::ofstream out(csv);
 out << std::setprecision(17) << "file,stage,yLow,yHigh,ptLow,ptHigh,status,covQual,edm,nfloat,chi2Ndf,massChi2Ndf,timeChi2Ndf,quality,B,Berr,nominal,nominalError,differencePercent,differenceSigma,tier,policyGroup,nGauss,nNP,nBkgSS,nBkgDS,nBkgFlip,minimizer\n";
 const double nan=std::numeric_limits<double>::quiet_NaN();
 std::regex tag("y([0-9.]+)_([0-9.]+)_pt([0-9.]+)_([0-9.]+)\\.root$");
 int n=0,good=0;
 for (const auto &entry:std::filesystem::recursive_directory_iterator(rootDir)) {
  if (!entry.is_regular_file() || entry.path().extension()!=".root") continue;
  const std::string path=entry.path().string(),stage=entry.path().parent_path().filename().string();
  std::smatch match;if(!std::regex_search(path,match,tag))continue;
  const double yl=std::stod(match[1]),yh=std::stod(match[2]),lo=std::stod(match[3]),hi=std::stod(match[4]);
  TFile f(path.c_str(),"READ");
  const bool final=stage=="fit2d";
  auto *r=f.Get<RooFitResult>(final?"modelResult":(stage=="mass"||stage=="mc_mass")?"fit_result":"timeResult");
  const double chi=audit_value(f,"prefitChi2NdfValue",nan);
  const int ndm=audit_value(f,"massNdf",0),ndt=audit_value(f,"timeNdf",0);
  const double cm=ndm>0?audit_value(f,"massChi2",nan)/ndm:nan,ct=ndt>0?audit_value(f,"timeChi2",nan)/ndt:nan;
  auto *backend=f.Get<TNamed>("configuredMinimizer");
  const bool commonMinimizer=backend&&std::string(backend->GetTitle())=="Minuit2";
  bool quality=(!final||commonMinimizer)&&r&&r->status()==0&&r->covQual()==3&&std::isfinite(r->edm())&&r->edm()>=0&&std::isfinite(r->minNll());
  quality &= final?(std::isfinite(cm)&&std::isfinite(ct)&&cm>=0&&ct>=0&&cm<2&&ct<2&&audit_value(f,"hesseStatus",-1)==0&&audit_value(f,"minimizeStatus",-1)==0):(std::isfinite(chi)&&chi>=0&&chi<2&&audit_value(f,"prefitChi2Ndf",0)>0);
  out<<path<<','<<stage<<','<<yl<<','<<yh<<','<<lo<<','<<hi<<','<<(r?r->status():-99)<<','<<(r?r->covQual():-99)<<','<<(r?r->edm():nan)<<','<<(r?r->floatParsFinal().getSize():0)<<','<<chi<<','<<cm<<','<<ct<<','<<(quality?"PASS":"FAIL");
  if(final){const auto &ref=nominal_bfraction::find(lo,hi,yl,yh);double b=audit_value(f,"bFraction",nan);out<<','<<b<<','<<audit_value(f,"bFractionErr",nan)<<','<<ref.value<<','<<ref.error<<','<<nominal_bfraction::differencePercent(b,ref)<<','<<nominal_bfraction::differenceSigma(b,ref)<<','<<nominal_bfraction::tier(b,ref);}
  else out<<",,,,,,,";
  auto *group=f.Get<TNamed>("fitPolicyGroup");
  out<<','<<(group?group->GetTitle():"")<<','<<audit_value(f,"nResolutionComponents",-1)
     <<','<<audit_value(f,"nSignalSSComponents",-1)<<','<<audit_value(f,"nBkgSSComponents",-1)
     <<','<<audit_value(f,"nBkgDSComponents",-1)<<','<<audit_value(f,"nBkgFlipComponents",-1)<<','<<(backend?backend->GetTitle():"");
  out<<'\n';++n;if(quality)++good;
 }
 std::cout<<"[AUDIT] "<<csv<<" files="<<n<<" quality="<<good<<"/"<<n<<std::endl;
}
