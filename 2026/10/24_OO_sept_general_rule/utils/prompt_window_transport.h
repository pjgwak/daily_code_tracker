#ifndef OO_PROMPT_WINDOW_TRANSPORT_H
#define OO_PROMPT_WINDOW_TRANSPORT_H
#include <RooAddModel.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooRealVar.h>
#include <TMD5.h>
#include <TSystem.h>
#include <cmath>
#include <TString.h>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace oo_prompt_transport {
struct CalibrationWindow {double low,high;std::string md5;};
inline CalibrationWindow calibrationWindow(const std::string&tag,const std::string&prefit){
 const char*overridePath=gSystem->Getenv("OO_PR_CALIBRATION_WINDOWS");const std::string path=overridePath?overridePath:"config/pr_calibration_windows.csv";
 std::ifstream in(path);if(!in)throw std::runtime_error("Cannot open PR calibration domain map: "+path);std::string line;std::getline(in,line);
 while(std::getline(in,line)){std::stringstream ss(line);std::vector<std::string> col;std::string s;while(std::getline(ss,s,','))col.push_back(s);if(col.size()<5||col[0]!=tag)continue;CalibrationWindow r{std::stod(col[1]),std::stod(col[2]),col[4]};if(!(std::isfinite(r.low)&&std::isfinite(r.high)&&r.high>r.low))throw std::runtime_error("Invalid PR calibration domain");std::unique_ptr<TMD5> md5(TMD5::FileChecksum(prefit.c_str()));if(!md5||r.md5!=md5->AsString())throw std::runtime_error("PR calibration map/pre-fit hash mismatch: "+prefit);return r;}
 throw std::runtime_error("PR calibration domain missing for "+tag);
}
// Native fractions remain fractions in the MC calibration window. All mappings
// depend on current native shape coordinates, leaving their constraints intact.
struct Model {
 std::vector<std::unique_ptr<RooFormulaVar>> terms;
 std::unique_ptr<RooAddModel> resolution;
 std::unique_ptr<RooAddPdf> prompt;
 std::vector<RooAbsReal*> nativeWeights,calProbability,windowProbability,fullWeights,windowWeights;
 RooAbsReal*fullSum=nullptr;RooAbsReal*windowSum=nullptr;
 RooFormulaVar&add(const std::string&n,const std::string&formula,const RooArgList&args){terms.emplace_back(new RooFormulaVar(n.c_str(),formula.c_str(),args));return *terms.back();}
 Model(RooRealVar&time,const RooArgList&gaussians,const RooArgList&means,const RooArgList&widths,const RooArgList&nativePhysicalFirst,const CalibrationWindow&cal,double low,double high){
  const int n=gaussians.size();if(n<1||n>4||means.size()!=n||widths.size()!=n||nativePhysicalFirst.size()!=n-1)throw std::runtime_error("Invalid PR transport component layout");
  RooArgList nativeArgs,fullArgs,windowArgs,fullCoef,windowCoef;
  std::string last="1";for(int j=0;j<n-1;j++){nativeArgs.add(*nativePhysicalFirst.at(j));last+="-@"+std::to_string(j);nativeWeights.push_back(static_cast<RooAbsReal*>(nativePhysicalFirst.at(j)));}
  nativeWeights.push_back(&add("prNativeLastWeight",last,nativeArgs));
  for(int j=0;j<n;j++){
   RooArgList shape(*means.at(j),*widths.at(j));
   auto&pc=add(Form("prCalibrationProbability%d",j+1),Form("0.5*(erf((%.17g-@0)/(sqrt(2)*@1))-erf((%.17g-@0)/(sqrt(2)*@1)))",cal.high,cal.low),shape);calProbability.push_back(&pc);
   auto&pw=add(Form("prWindowProbability%d",j+1),Form("0.5*(erf((%.17g-@0)/(sqrt(2)*@1))-erf((%.17g-@0)/(sqrt(2)*@1)))",high,low),shape);windowProbability.push_back(&pw);
  }

  for(int j=0;j<n;j++){
   auto&u=add(Form("prPhysicalNumerator%d",j+1),"@0/@1",RooArgList(*nativeWeights[j],*calProbability[j]));fullWeights.push_back(&u);fullArgs.add(u);
   auto&q=add(Form("prTargetNumerator%d",j+1),"@0*@1",RooArgList(u,*windowProbability[j]));windowWeights.push_back(&q);windowArgs.add(q);
  }
  std::string sum="@0";for(int j=1;j<n;j++)sum+="+@"+std::to_string(j);
  fullSum=&add("prFullWeightSum",sum,fullArgs);windowSum=&add("prTargetWeightSum",sum,windowArgs);
  for(int j=0;j<n-1;j++){fullCoef.add(add(Form("prPhysicalFraction%d",j+1),"@0/@1",RooArgList(*fullWeights[j],*fullSum)));windowCoef.add(add(Form("prTargetFraction%d",j+1),"@0/@1",RooArgList(*windowWeights[j],*windowSum)));}
  resolution.reset(new RooAddModel("promptResolutionModel","Physical full-space PR resolution for convolution",gaussians,fullCoef));
  prompt.reset(new RooAddPdf("promptWindowPdf","PR mixture conditioned on target ctau window",gaussians,windowCoef));
 }
};
}
#endif
