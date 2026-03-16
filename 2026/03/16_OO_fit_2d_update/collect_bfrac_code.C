/*  Use this command line
// Fwd
root -l -b -q 'collect_bfrac_code.C(1.6,2.4,"roots","oo_y16to24","NooY16to24",true)'

// Mid
root -l -b -q 'collect_bfrac_code.C(0.0,1.6,"roots","oo_y16","NooY16",true)'
*/ 



#include "TFile.h"
#include "TParameter.h"
#include "TString.h"
#include "TSystem.h"
#include "RooFitResult.h"
#include "RooRealVar.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace
{
struct BinResult
{
  double yLow;
  double yHigh;
  double ptLow;
  double ptHigh;
  double value;
  double error;
  TString source;
};

template <typename T>
bool readParam(TFile &f, const char *name, T &out)
{
  auto *param = dynamic_cast<TParameter<T> *>(f.Get(name));
  if (!param)
    return false;
  out = param->GetVal();
  return true;
}

bool readResultFile(const TString &path, BinResult &out)
{
  static const std::regex kNameRe(
      R"(fit2d_result_y([0-9]+\.[0-9]+)_([0-9]+\.[0-9]+)_pt([0-9]+\.[0-9]+)_([0-9]+\.[0-9]+)\.root)");

  std::cmatch match;
  const std::string name = gSystem->BaseName(path.Data());
  if (!std::regex_match(name.c_str(), match, kNameRe))
    return false;

  std::unique_ptr<TFile> f(TFile::Open(path));
  if (!f || f->IsZombie())
    return false;

  double value = 0.0;
  double error = 0.0;
  const bool hasValue = readParam(*f, "bFraction", value);
  const bool hasError = readParam(*f, "bFractionErr", error);

  if (!(hasValue && hasError))
  {
    auto *modelResult = dynamic_cast<RooFitResult *>(f->Get("modelResult"));
    if (!modelResult)
      return false;
    auto *var = dynamic_cast<RooRealVar *>(modelResult->floatParsFinal().find("bFraction"));
    if (!var)
      return false;
    value = var->getVal();
    error = var->getError();
  }

  out.yLow = std::stod(match[1].str());
  out.yHigh = std::stod(match[2].str());
  out.ptLow = std::stod(match[3].str());
  out.ptHigh = std::stod(match[4].str());
  out.value = value;
  out.error = error;
  out.source = path;
  return true;
}

void printArray(const char *label, const std::vector<double> &vals, int decimals)
{
  std::cout << label << "{";
  for (size_t i = 0; i < vals.size(); ++i)
  {
    if (i)
      std::cout << ", ";
    std::cout << Form(("%." + std::to_string(decimals) + "f").c_str(), vals[i]);
  }
  std::cout << "};\n";
}
} // namespace

void collect_bfrac_code(double yLow = 1.6,
                        double yHigh = 2.4,
                        TString rootsDir = "roots",
                        TString prefix = "pp",
                        TString countName = "Npp",
                        bool showSources = false)
{
  void *dirp = gSystem->OpenDirectory(rootsDir);
  if (!dirp)
  {
    std::cerr << "[ERROR] roots directory not found: " << rootsDir << std::endl;
    return;
  }
  gSystem->FreeDirectory(dirp);

  std::vector<TString> files;
  const TString cmd = Form("find %s -type f -name 'fit2d_result_*.root'", rootsDir.Data());
  std::unique_ptr<TObjArray> lines(gSystem->GetFromPipe(cmd).Tokenize("\n"));
  for (int i = 0; lines && i < lines->GetEntriesFast(); ++i)
  {
    auto *obj = lines->At(i);
    if (!obj)
      continue;
    const TString path = obj->GetName();
    if (!path.IsNull())
      files.push_back(path);
  }

  std::vector<BinResult> selected;
  for (const auto &path : files)
  {
    BinResult result{};
    if (!readResultFile(path, result))
      continue;
    if (std::abs(result.yLow - yLow) > 1e-9 || std::abs(result.yHigh - yHigh) > 1e-9)
      continue;
    selected.push_back(result);
  }

  std::sort(selected.begin(), selected.end(), [](const BinResult &a, const BinResult &b) {
    if (a.ptLow != b.ptLow)
      return a.ptLow < b.ptLow;
    return a.ptHigh < b.ptHigh;
  });

  if (selected.empty())
  {
    std::cerr << Form("[ERROR] no fit2d ROOT results found for y=%.2f-%.2f under %s",
                      yLow, yHigh, rootsDir.Data())
              << std::endl;
    return;
  }

  std::vector<double> ptLowVals, ptHighVals, bfracVals, bfracErrs;
  for (const auto &item : selected)
  {
    ptLowVals.push_back(item.ptLow);
    ptHighVals.push_back(item.ptHigh);
    bfracVals.push_back(item.value);
    bfracErrs.push_back(item.error);
    if (showSources)
    {
      std::cerr << Form("[SRC] y=%.2f-%.2f pt=%.2f-%.2f <- %s",
                        item.yLow, item.yHigh, item.ptLow, item.ptHigh, item.source.Data())
                << std::endl;
    }
  }

  std::cout << Form("// y = %.2f-%.2f\n", yLow, yHigh);
  std::cout << Form("const int %s = %zu;\n", countName.Data(), selected.size());
  printArray(Form("double pt_low_%s[%s] = ", prefix.Data(), countName.Data()), ptLowVals, 1);
  printArray(Form("double pt_high_%s[%s] = ", prefix.Data(), countName.Data()), ptHighVals, 1);
  std::cout << "\n";
  printArray(Form("double val_%s[%s] = ", prefix.Data(), countName.Data()), bfracVals, 5);
  printArray(Form("double stat_%s[%s] = ", prefix.Data(), countName.Data()), bfracErrs, 5);
  std::cout << "\n";
  std::cout << Form("double sys_%s[%s];\n", prefix.Data(), countName.Data());
  std::cout << Form("for (int i = 0; i < %s; ++i)\n", countName.Data());
  std::cout << Form("  sys_%s[i] = 0.0; // sys=0\n\n", prefix.Data());
  std::cout << Form("double x_%s[%s], ex_%s[%s];\n", prefix.Data(), countName.Data(), prefix.Data(), countName.Data());
  std::cout << Form("for (int i = 0; i < %s; ++i)\n", countName.Data());
  std::cout << "{\n";
  std::cout << Form("  x_%s[i] = 0.5 * (pt_low_%s[i] + pt_high_%s[i]);\n", prefix.Data(), prefix.Data(), prefix.Data());
  std::cout << Form("  ex_%s[i] = 0.5 * (pt_high_%s[i] - pt_low_%s[i]);\n", prefix.Data(), prefix.Data(), prefix.Data());
  std::cout << "}\n";
}
