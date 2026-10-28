#ifndef OO_NP_MOMENT_CONSTRAINT_H
#define OO_NP_MOMENT_CONSTRAINT_H
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooWrapperPdf.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace oo_np_moment {
inline std::string number(double x){std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();}
// E[t^k]/k! for a normalized exponential in the native truth interval.
inline double componentMoment(double tau,int k,double low,double high){
  if(!(tau>0&&high>low&&low>=0&&k>=0))throw std::runtime_error("Invalid NP moment domain");
  auto sum=[&](double x){double s=1,t=1;for(int j=1;j<=k;++j){t*=x/j;s+=t;}return s;};
  return std::pow(tau,k)*(std::exp(-low/tau)*sum(low/tau)-std::exp(-high/tau)*sum(high/tau))
       /(std::exp(-low/tau)-std::exp(-high/tau));
}
inline std::string componentFormula(int k,double low,double high){
  auto series=[&](double edge){std::string s="1";double fact=1;for(int j=1;j<=k;++j){fact*=j;s+="+pow("+number(edge)+"/@0,"+std::to_string(j)+")/"+number(fact);}return "("+s+")";};
  return "pow(@0,"+std::to_string(k)+")*(exp(-"+number(low)+"/@0)*"+series(low)+"-exp(-"+number(high)+"/@0)*"+series(high)+")/(exp(-"+number(low)+"/@0)-exp(-"+number(high)+"/@0))";
}
inline std::vector<double> logCoordinates(const std::vector<double>&v,const std::vector<double>&floors,int order,double low,double high){
  const int n=floors.size();double total=0;for(int i=0;i<n;++i)total+=v[2*i+1];
  if(!(total>0))throw std::runtime_error("Invalid NP total normalization");
  std::vector<double> result{std::log(total)};
  for(int k=1;k<=order;++k){double m=0;for(int i=0;i<n;++i)m+=v[2*i+1]*componentMoment(floors[i]+std::exp(v[2*i]),k,low,high);result.push_back(std::log(m/total));}
  return result;
}
struct Constraint {
  std::vector<std::unique_ptr<RooFormulaVar>> terms;
  std::vector<RooAbsReal*> coordinates;
  std::unique_ptr<RooWrapperPdf> pdf;
  TMatrixDSym covariance,correlationInverse;
  std::vector<double> centers,errors;
  std::map<std::string,double> calibrationStarts;
  double initialValue(const std::string&name)const{
    auto it=calibrationStarts.find(name);if(it==calibrationStarts.end())throw std::runtime_error("Missing NP moment initial coordinate: "+name);return it->second;
  }
  RooFormulaVar&add(const std::string&name,const std::string&formula,const RooArgList&args){
    terms.emplace_back(new RooFormulaVar(name.c_str(),formula.c_str(),args));return *terms.back();
  }
  Constraint(const RooFitResult&mc,const std::vector<RooRealVar*>&logs,const std::vector<RooRealVar*>&yields,
             const std::vector<double>&floors,int order,double scale,double low,double high,const std::string&auditPrefix)
    :covariance(order+1),correlationInverse(order+1){
    const int n=logs.size(),dim=order+1;
    if(n<2||n>3||int(yields.size())!=n||int(floors.size())!=n||order<1||order>3||order>=2*n||!(scale>0)||mc.status()!=0||mc.covQual()!=3)
      throw std::runtime_error("Invalid NP moment constraint configuration");
    std::vector<int> indices;std::vector<double> native;std::vector<RooRealVar*> variables;
    for(int i=0;i<n;++i)for(auto*v:{logs[i],yields[i]}){
      int index=mc.floatParsFinal().index(v->GetName());if(index<0)throw std::runtime_error("Missing NP native covariance coordinate");
      auto *saved=dynamic_cast<const RooRealVar*>(mc.floatParsFinal().at(index));if(!saved)throw std::runtime_error("Invalid NP native variable");
      v->setVal(saved->getVal());v->setError(saved->getError());v->setConstant(false);
      indices.push_back(index);native.push_back(saved->getVal());variables.push_back(v);
      calibrationStarts[v->GetName()]=saved->getVal();
    }
    centers=logCoordinates(native,floors,order,low,high);
    double total=std::exp(centers[0]);TMatrixD jac(dim,2*n);
    for(int i=0;i<n;++i){
      const double tau=floors[i]+std::exp(native[2*i]),yield=native[2*i+1];
      jac(0,2*i)=0;jac(0,2*i+1)=1/total;
      for(int k=1;k<=order;++k){
        double mi=componentMoment(tau,k,low,high),mean=componentMoment(tau,1,low,high),next=componentMoment(tau,k+1,low,high),mix=std::exp(centers[k]);
        jac(k,2*i)=yield/total*((k+1)*next-mi*mean)/(tau*tau)*std::exp(native[2*i])/mix;
        jac(k,2*i+1)=(mi-mix)/total/mix;
      }
    }
    std::ofstream ja(auditPrefix+"_JACOBIAN.csv"),ca(auditPrefix+"_COVARIANCE.csv"),na(auditPrefix+"_NATIVE.csv"),me(auditPrefix+"_CENTERS.csv"),nc(auditPrefix+"_NATIVE_COVARIANCE.csv");
    ja<<std::setprecision(17)<<"coordinate,native,analytic,numeric\n";
    ca<<std::setprecision(17)<<"i,j,covariance,correlation\n";
    na<<std::setprecision(17)<<"component,log_lifetime,yield,floor,truth_low,truth_high\n";
    me<<std::setprecision(17)<<"coordinate,center,error,covariance_scale\n";
    nc<<std::setprecision(17)<<"i,j,covariance\n";
    for(int i=0;i<2*n;++i)for(int j=0;j<2*n;++j)nc<<i<<','<<j<<','<<mc.covarianceMatrix()(indices[i],indices[j])<<'\n';
    for(int i=0;i<n;++i)na<<i+1<<','<<native[2*i]<<','<<native[2*i+1]<<','<<floors[i]<<','<<low<<','<<high<<'\n';
    for(int j=0;j<2*n;++j){
      double h=1e-5*std::max(1.,std::abs(native[j]));auto plus=native,minus=native;plus[j]+=h;minus[j]-=h;
      auto a=logCoordinates(plus,floors,order,low,high),b=logCoordinates(minus,floors,order,low,high);
      for(int i=0;i<dim;++i){double numeric=(a[i]-b[i])/(2*h);ja<<i<<','<<variables[j]->GetName()<<','<<jac(i,j)<<','<<numeric<<'\n';
        if(std::abs(jac(i,j)-numeric)>1e-9+1e-5*std::abs(numeric))throw std::runtime_error("NP moment Jacobian verification failed");}
    }
    for(int i=0;i<dim;++i)for(int j=0;j<dim;++j){
      covariance(i,j)=0;for(int a=0;a<2*n;++a)for(int b=0;b<2*n;++b)covariance(i,j)+=scale*jac(i,a)*mc.covarianceMatrix()(indices[a],indices[b])*jac(j,b);
    }
    for(int i=0;i<dim;++i){if(!(covariance(i,i)>0))throw std::runtime_error("Invalid NP moment variance");errors.push_back(std::sqrt(covariance(i,i)));me<<i<<','<<centers[i]<<','<<errors[i]<<','<<scale<<'\n';}
    TMatrixDSym corr(dim);
    for(int i=0;i<dim;++i)for(int j=0;j<dim;++j){corr(i,j)=covariance(i,j)/(errors[i]*errors[j]);ca<<i<<','<<j<<','<<covariance(i,j)<<','<<corr(i,j)<<'\n';}
    TDecompChol check(corr);if(!check.Decompose())throw std::runtime_error("NP moment covariance is not positive definite");correlationInverse=check.Invert();
    RooArgList yieldArgs;std::string sum;
    for(int i=0;i<n;++i){yieldArgs.add(*yields[i]);sum+=(i?"+":"")+std::string("@")+std::to_string(i);}
    auto&norm=add("npMomentTotalYield",sum,yieldArgs);coordinates.push_back(&add("npMomentLogYield","log(@0)",RooArgList(norm)));
    std::vector<RooAbsReal*>taus;
    for(int i=0;i<n;++i)taus.push_back(&add("npMomentTau"+std::to_string(i+1),number(floors[i])+"+exp(@0)",RooArgList(*logs[i])));
    for(int k=1;k<=order;++k){
      RooArgList args(norm);std::string formula="log((";
      for(int i=0;i<n;++i){
        auto&mi=add("npMomentComponent"+std::to_string(i+1)+"Order"+std::to_string(k),componentFormula(k,low,high),RooArgList(*taus[i]));
        args.add(*yields[i]);args.add(mi);formula+=(i?"+":"")+std::string("@")+std::to_string(2*i+1)+"*@"+std::to_string(2*i+2);
      }
      formula+=")/@0)";coordinates.push_back(&add("npMomentLogOrder"+std::to_string(k),formula,args));
    }
    RooArgList residuals;
    for(int i=0;i<dim;++i){
      if(std::abs(coordinates[i]->getVal()-centers[i])>1e-10)throw std::runtime_error("NP moment formula mismatch");
      residuals.add(add("npMomentResidual"+std::to_string(i),"(@0-("+number(centers[i])+"))/"+number(errors[i]),RooArgList(*coordinates[i])));
    }
    std::string penalty="0";
    for(int i=0;i<dim;++i)for(int j=0;j<dim;++j)penalty+="+("+number(.5*correlationInverse(i,j))+")*@"+std::to_string(i)+"*@"+std::to_string(j);
    auto&p=add("npMomentPenalty",penalty,residuals);auto&w=add("npMomentWeight","exp(-@0)",RooArgList(p));
    pdf.reset(new RooWrapperPdf("npMomentConstraint","NP native MC moment covariance constraint",w,true));
    RooArgSet normSet;for(auto*v:variables)normSet.add(*v);
    std::ofstream validation(auditPrefix+"_PENALTY_CHECK.csv");validation<<std::setprecision(17)<<"variable,shift_sigma,expected,actual\n";
    for(int j=0;j<2*n;++j)for(double shift:{-.25,0.,.25}){
      double value=native[j]+shift*variables[j]->getError();variables[j]->setVal(value);
      auto v=native;v[j]=value;auto c=logCoordinates(v,floors,order,low,high);double expected=0;
      for(int i=0;i<dim;++i)for(int k=0;k<dim;++k)expected+=.5*(c[i]-centers[i])/errors[i]*correlationInverse(i,k)*(c[k]-centers[k])/errors[k];
      double actual=-pdf->getLogVal(&normSet);validation<<variables[j]->GetName()<<','<<shift<<','<<expected<<','<<actual<<'\n';
      if(!std::isfinite(actual)||std::abs(actual-expected)>1e-8*(1+std::abs(expected)))throw std::runtime_error("NP moment penalty normalization mismatch");
      variables[j]->setVal(native[j]);
    }
  }
  // A single effective DATA exponential is calibrated to the mean of the full
  // native MC mixture. MC model selection still uses the original valid fit.
  Constraint(const RooFitResult&mc,RooRealVar&logTau,RooRealVar&yield,const std::vector<double>&nativeFloors,
             double dataFloor,double scale,double low,double high,const std::string&auditPrefix)
    :covariance(2),correlationInverse(2){
    std::vector<std::unique_ptr<RooRealVar>> clones;std::vector<RooRealVar*>logs,yields;
    for(unsigned i=0;i<nativeFloors.size();++i){
      for(const std::string&name:{i?"signal"+std::to_string(i+1)+"_log_lifetime_offset":std::string("signal_log_lifetime_offset"),"NsignalSS"+std::to_string(i+1)}){
        auto*v=dynamic_cast<const RooRealVar*>(mc.floatParsFinal().find(name.c_str()));if(!v)throw std::runtime_error("Missing native NP mixture coordinate");
        clones.emplace_back(static_cast<RooRealVar*>(v->Clone()));
      }
      logs.push_back(clones[2*i].get());yields.push_back(clones[2*i+1].get());
    }
    Constraint native(mc,logs,yields,nativeFloors,1,scale,low,high,auditPrefix+"_MC");
    covariance=native.covariance;correlationInverse=native.correlationInverse;centers=native.centers;errors=native.errors;
    const double mean=std::exp(centers[1]),total=std::exp(centers[0]);
    double left=dataFloor+std::exp(logTau.getMin()),right=dataFloor+std::exp(logTau.getMax());
    if(!(componentMoment(left,1,low,high)<mean&&componentMoment(right,1,low,high)>mean))throw std::runtime_error("Effective NP mean lies outside lifetime domain");
    for(int i=0;i<100;++i){double middle=.5*(left+right);if(componentMoment(middle,1,low,high)>mean)right=middle;else left=middle;}
    const double tau=.5*(left+right),initial=std::log(tau-dataFloor);
    logTau.setVal(initial);logTau.setConstant(false);yield.setVal(total);yield.setConstant(false);
    double derivative=(2*componentMoment(tau,2,low,high)-mean*mean)/(tau*tau)*(tau-dataFloor)/mean;
    logTau.setError(errors[1]/std::abs(derivative));yield.setError(total*errors[0]);
    calibrationStarts[logTau.GetName()]=initial;calibrationStarts[yield.GetName()]=total;
    auto&t=add("npMomentEffectiveTau",number(dataFloor)+"+exp(@0)",RooArgList(logTau));
    auto&m=add("npMomentEffectiveMean",componentFormula(1,low,high),RooArgList(t));
    coordinates.push_back(&add("npMomentLogYield","log(@0)",RooArgList(yield)));
    coordinates.push_back(&add("npMomentLogOrder1","log(@0)",RooArgList(m)));
    RooArgList residuals;
    for(int i=0;i<2;++i){
      if(std::abs(coordinates[i]->getVal()-centers[i])>1e-10)throw std::runtime_error("Effective NP moment mismatch");
      residuals.add(add("npMomentResidual"+std::to_string(i),"(@0-("+number(centers[i])+"))/"+number(errors[i]),RooArgList(*coordinates[i])));
    }
    std::string penalty="0";
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)penalty+="+("+number(.5*correlationInverse(i,j))+")*@"+std::to_string(i)+"*@"+std::to_string(j);
    auto&p=add("npMomentPenalty",penalty,residuals);auto&w=add("npMomentWeight","exp(-@0)",RooArgList(p));
    pdf.reset(new RooWrapperPdf("npMomentConstraint","NP effective mean from full native MC mixture",w,true));
    std::ofstream audit(auditPrefix+"_EFFECTIVE.csv"),check(auditPrefix+"_EFFECTIVE_CHECK.csv");
    audit<<std::setprecision(17)<<"native_components,truth_low,truth_high,mean,effective_tau,log_lifetime,total_yield,log_mean_error,log_total_error\n";
    audit<<nativeFloors.size()<<','<<low<<','<<high<<','<<mean<<','<<tau<<','<<initial<<','<<total<<','<<errors[1]<<','<<errors[0]<<'\n';
    check<<std::setprecision(17)<<"shift_sigma,mean,expected,actual\n";RooArgSet norm(logTau,yield);
    for(double shift:{-1.,0.,1.}){
      logTau.setVal(initial+shift*logTau.getError());double x=std::log(componentMoment(dataFloor+std::exp(logTau.getVal()),1,low,high));
      double expected=.5*std::pow((x-centers[1])/errors[1],2)*correlationInverse(1,1),actual=-pdf->getLogVal(&norm);
      check<<shift<<','<<std::exp(x)<<','<<expected<<','<<actual<<'\n';
      if(std::abs(actual-expected)>1e-8)throw std::runtime_error("Effective NP penalty mismatch");
    }
    logTau.setVal(initial);
  }
};
}
#endif
