/*
Usage from this directory:

  root -l -b -q 'check_roofit_results.C("../../roots",1e-3,true)'

Arguments:
  rootsDir : directory containing the rapidity-bin output subdirectories
  edmMax   : maximum accepted RooFitResult::edm()
  showAll  : true prints every result; false prints only failed/missing results
  csvPath  : CSV output. A relative path is resolved from this macro directory.
             Pass "" to disable CSV output.

PASS requires all of:
  RooFitResult::status() == 0
  every status-history code == 0
  RooFitResult::covQual() == 3
  finite RooFitResult::edm() <= edmMax

The checker reads RooFitResult objects from the ROOT files, not text logs.
Older ctau_pr/ctau_np files made before timeResult was added are reported as
NO_ROOFITRESULT and must be regenerated before they can be checked.
*/

#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include "RooFitResult.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace
{
TString macroDir()
{
  TString path = __FILE__;
  if (!gSystem->IsAbsoluteFileName(path))
    path = TString(gSystem->WorkingDirectory()) + "/" + path;
  return gSystem->DirName(path);
}

TString existingDir(const TString &path)
{
  void *dirp = gSystem->OpenDirectory(path);
  if (!dirp)
    return "";
  gSystem->FreeDirectory(dirp);
  return path;
}

TString resolveRootsDir(const TString &rootsDir)
{
  TString dir = existingDir(rootsDir);
  if (!dir.IsNull())
    return dir;

  if (!gSystem->IsAbsoluteFileName(rootsDir))
  {
    dir = existingDir(macroDir() + "/" + rootsDir);
    if (!dir.IsNull())
      return dir;
  }

  if (rootsDir == "roots")
  {
    dir = existingDir(macroDir() + "/../../roots");
    if (!dir.IsNull())
      return dir;
  }

  return "";
}

TString resolveOutputPath(const TString &path)
{
  if (path.IsNull() || gSystem->IsAbsoluteFileName(path))
    return path;
  return macroDir() + "/" + path;
}

struct MacroSpec
{
  TString name;
  TString resultKey;
  int order = 999;
};

bool identifyMacro(const TString &path, MacroSpec &out)
{
  if (path.Contains("/mc_mass/"))
    out = {"mc_mass", "fit_result", 0};
  else if (path.Contains("/mass/"))
    out = {"mass", "fit_result", 1};
  else if (path.Contains("/ctau_pr/"))
    out = {"ctau_pr", "timeResult", 2};
  else if (path.Contains("/ctau_np/"))
    out = {"ctau_np", "timeResult", 3};
  else if (path.Contains("/ctau_bkg/"))
    out = {"ctau_bkg", "timeResult", 4};
  else if (path.Contains("/fit2d/"))
    out = {"fit2d", "modelResult", 5};
  else
    return false;
  return true;
}

struct CheckResult
{
  MacroSpec macro;
  double yLow = 0.0;
  double yHigh = 0.0;
  double ptLow = 0.0;
  double ptHigh = 0.0;
  int status = -999;
  int covQual = -999;
  double edm = std::numeric_limits<double>::quiet_NaN();
  double minNll = std::numeric_limits<double>::quiet_NaN();
  TString history;
  TString verdict;
  TString reason;
  TString source;
};

bool parseBin(const TString &path, CheckResult &out)
{
  static const std::regex kNameRe(
      R"(.*_y(-?[0-9]+\.[0-9]+)_(-?[0-9]+\.[0-9]+)_pt(-?[0-9]+\.[0-9]+)_(-?[0-9]+\.[0-9]+)\.root)");

  std::cmatch match;
  const std::string name = gSystem->BaseName(path.Data());
  if (!std::regex_match(name.c_str(), match, kNameRe))
    return false;

  out.yLow = std::stod(match[1].str());
  out.yHigh = std::stod(match[2].str());
  out.ptLow = std::stod(match[3].str());
  out.ptHigh = std::stod(match[4].str());
  return true;
}

RooFitResult *findFitResult(TFile &file, const TString &preferredKey, TString &usedKey)
{
  if (!preferredKey.IsNull())
  {
    auto *result = dynamic_cast<RooFitResult *>(file.Get(preferredKey));
    if (result)
    {
      usedKey = preferredKey;
      return result;
    }
  }

  TIter next(file.GetListOfKeys());
  while (auto *key = dynamic_cast<TKey *>(next()))
  {
    if (TString(key->GetClassName()) != "RooFitResult")
      continue;
    auto *result = dynamic_cast<RooFitResult *>(file.Get(key->GetName()));
    if (!result)
      continue;
    usedKey = key->GetName();
    return result;
  }
  return nullptr;
}

TString statusHistory(const RooFitResult &result, bool &allZero)
{
  allZero = true;
  std::ostringstream text;
  for (UInt_t i = 0; i < result.numStatusHistory(); ++i)
  {
    if (i)
      text << ";";
    const char *label = result.statusLabelHistory(i);
    const int code = result.statusCodeHistory(i);
    text << (label ? label : "UNKNOWN") << "=" << code;
    if (code != 0)
      allZero = false;
  }
  return text.str();
}

void classify(CheckResult &item, double edmMax, bool historyAllZero)
{
  std::vector<std::string> failures;
  if (item.status != 0)
    failures.push_back(Form("status=%d", item.status));
  if (!historyAllZero)
    failures.push_back("history_nonzero");
  if (item.covQual != 3)
    failures.push_back(Form("covQual=%d", item.covQual));
  if (!std::isfinite(item.edm))
    failures.push_back("edm_nonfinite");
  else if (item.edm > edmMax)
    failures.push_back(Form("edm=%.6g>%.6g", item.edm, edmMax));

  item.verdict = failures.empty() ? "PASS" : "FAIL";
  std::ostringstream reason;
  for (size_t i = 0; i < failures.size(); ++i)
  {
    if (i)
      reason << ";";
    reason << failures[i];
  }
  item.reason = failures.empty() ? "-" : reason.str();
}

CheckResult inspectFile(const TString &path, double edmMax)
{
  CheckResult out;
  out.source = path;
  if (!identifyMacro(path, out.macro) || !parseBin(path, out))
  {
    out.verdict = "SKIP";
    out.reason = "unrecognized_output_name";
    return out;
  }

  std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
  if (!file || file->IsZombie())
  {
    out.verdict = "FILE_ERROR";
    out.reason = "cannot_open";
    return out;
  }

  TString usedKey;
  RooFitResult *result = findFitResult(*file, out.macro.resultKey, usedKey);
  if (!result)
  {
    out.verdict = "NO_ROOFITRESULT";
    out.reason = "expected_key=" + out.macro.resultKey;
    return out;
  }

  out.status = result->status();
  out.covQual = result->covQual();
  out.edm = result->edm();
  out.minNll = result->minNll();
  bool historyAllZero = true;
  out.history = statusHistory(*result, historyAllZero);
  classify(out, edmMax, historyAllZero);
  return out;
}

TString csvEscape(const TString &input)
{
  TString escaped = input;
  escaped.ReplaceAll("\"", "\"\"");
  return "\"" + escaped + "\"";
}
} // namespace

void check_roofit_results(TString rootsDir = "../../roots",
                          double edmMax = 1e-3,
                          bool showAll = true,
                          TString csvPath = "fit_result_summary.csv")
{
  const TString resolvedRootsDir = resolveRootsDir(rootsDir);
  if (resolvedRootsDir.IsNull())
  {
    std::cerr << "[ERROR] roots directory not found: " << rootsDir << std::endl;
    return;
  }

  std::vector<TString> files;
  const TString command = Form("find %s -type f -name '*.root'", resolvedRootsDir.Data());
  std::unique_ptr<TObjArray> lines(gSystem->GetFromPipe(command).Tokenize("\n"));
  for (int i = 0; lines && i < lines->GetEntriesFast(); ++i)
  {
    auto *line = dynamic_cast<TObjString *>(lines->At(i));
    if (!line || line->GetString().IsNull())
      continue;
    MacroSpec spec;
    if (identifyMacro(line->GetString(), spec))
      files.push_back(line->GetString());
  }

  std::vector<CheckResult> results;
  results.reserve(files.size());
  for (const auto &path : files)
  {
    CheckResult item = inspectFile(path, edmMax);
    if (item.verdict != "SKIP")
      results.push_back(item);
  }

  std::sort(results.begin(), results.end(), [](const CheckResult &a, const CheckResult &b) {
    if (a.yLow != b.yLow)
      return a.yLow < b.yLow;
    if (a.yHigh != b.yHigh)
      return a.yHigh < b.yHigh;
    if (a.ptLow != b.ptLow)
      return a.ptLow < b.ptLow;
    if (a.ptHigh != b.ptHigh)
      return a.ptHigh < b.ptHigh;
    return a.macro.order < b.macro.order;
  });

  std::cout << "[INFO] roots=" << resolvedRootsDir << ", edmMax=" << edmMax << "\n";
  std::cout << "verdict  macro      y-bin        pt-bin       status covQual edm          history\n";
  std::cout << "--------------------------------------------------------------------------------------\n";

  std::map<std::string, std::map<std::string, int>> counts;
  for (const auto &item : results)
  {
    counts[item.macro.name.Data()][item.verdict.Data()]++;
    if (!showAll && item.verdict == "PASS")
      continue;

    std::cout << Form("%-8s %-10s %4.1f-%-4.1f  %5.1f-%-5.1f  %6d %7d %-12.5g %s",
                      item.verdict.Data(), item.macro.name.Data(),
                      item.yLow, item.yHigh, item.ptLow, item.ptHigh,
                      item.status, item.covQual, item.edm, item.history.Data())
              << "\n";
    if (item.verdict != "PASS")
      std::cout << "         reason: " << item.reason << "\n";
  }

  std::cout << "\n[SUMMARY]\n";
  const std::vector<std::string> macroOrder = {
      "mc_mass", "mass", "ctau_pr", "ctau_np", "ctau_bkg", "fit2d"};
  int totalPass = 0;
  int totalFail = 0;
  int totalMissing = 0;
  int totalFileError = 0;
  for (const auto &name : macroOrder)
  {
    const auto &c = counts[name];
    const int pass = c.count("PASS") ? c.at("PASS") : 0;
    const int fail = c.count("FAIL") ? c.at("FAIL") : 0;
    const int missing = c.count("NO_ROOFITRESULT") ? c.at("NO_ROOFITRESULT") : 0;
    const int fileError = c.count("FILE_ERROR") ? c.at("FILE_ERROR") : 0;
    totalPass += pass;
    totalFail += fail;
    totalMissing += missing;
    totalFileError += fileError;
    std::cout << Form("%-10s PASS=%2d FAIL=%2d NO_ROOFITRESULT=%2d FILE_ERROR=%2d",
                      name.c_str(), pass, fail, missing, fileError)
              << "\n";
  }
  std::cout << Form("TOTAL      PASS=%2d FAIL=%2d NO_ROOFITRESULT=%2d FILE_ERROR=%2d",
                    totalPass, totalFail, totalMissing, totalFileError)
            << "\n";

  const TString resolvedCsvPath = resolveOutputPath(csvPath);
  if (resolvedCsvPath.IsNull())
    return;

  std::ofstream csv(resolvedCsvPath.Data());
  if (!csv)
  {
    std::cerr << "[ERROR] cannot create CSV: " << resolvedCsvPath << std::endl;
    return;
  }

  csv << "verdict,macro,yLow,yHigh,ptLow,ptHigh,status,covQual,edm,minNll,history,reason,source\n";
  csv << std::setprecision(12);
  for (const auto &item : results)
  {
    csv << item.verdict << ","
        << item.macro.name << ","
        << item.yLow << ","
        << item.yHigh << ","
        << item.ptLow << ","
        << item.ptHigh << ","
        << item.status << ","
        << item.covQual << ","
        << item.edm << ","
        << item.minNll << ","
        << csvEscape(item.history) << ","
        << csvEscape(item.reason) << ","
        << csvEscape(item.source) << "\n";
  }
  csv.close();
  std::cout << "[OUTPUT] " << resolvedCsvPath << std::endl;
}
