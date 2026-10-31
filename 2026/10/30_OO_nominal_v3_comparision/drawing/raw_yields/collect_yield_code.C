/*  Use this command line
// Fwd (the defaults already produce the variable names used by raw_yield_fwd.C)
root -l -b -q 'collect_yield_code.C(1.6,2.4)'
root -l -b -q 'collect_yield_code.C(1.6,2.4,"../../roots","fit","fitBinCount",false)'

// Mid
root -l -b -q 'collect_yield_code.C(0.0,1.6)'
root -l -b -q 'collect_yield_code.C(0.0,1.6,"../../roots","fit","fitBinCount",false)'
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
TString macroDirectory()
{
  TString macroPath = __FILE__;
  if (!gSystem->IsAbsoluteFileName(macroPath))
    macroPath = TString(gSystem->WorkingDirectory()) + "/" + macroPath;
  return gSystem->DirName(macroPath);
}

TString findExistingDirectory(const TString &directoryPath)
{
  void *directoryHandle = gSystem->OpenDirectory(directoryPath);
  if (!directoryHandle)
    return "";
  gSystem->FreeDirectory(directoryHandle);
  return directoryPath;
}

TString resolveResultsDirectory(const TString &resultsDirectory)
{
  TString resolvedDirectory = findExistingDirectory(resultsDirectory);
  if (!resolvedDirectory.IsNull())
    return resolvedDirectory;

  if (!gSystem->IsAbsoluteFileName(resultsDirectory))
  {
    const TString pathRelativeToMacro = macroDirectory() + "/" + resultsDirectory;
    resolvedDirectory = findExistingDirectory(pathRelativeToMacro);
    if (!resolvedDirectory.IsNull())
      return resolvedDirectory;
  }

  if (resultsDirectory == "roots")
  {
    const TString defaultResultsDirectory = macroDirectory() + "/../../roots";
    resolvedDirectory = findExistingDirectory(defaultResultsDirectory);
    if (!resolvedDirectory.IsNull())
      return resolvedDirectory;
  }

  return "";
}

struct BinResult
{
  double rapidityLow;
  double rapidityHigh;
  double ptLow;
  double ptHigh;
  double yieldValue;
  double yieldError;
  TString sourcePath;
};

template <typename T>
bool readParameter(TFile &inputFile, const char *parameterName, T &parameterValue)
{
  auto *parameter = dynamic_cast<TParameter<T> *>(inputFile.Get(parameterName));
  if (!parameter)
    return false;
  parameterValue = parameter->GetVal();
  return true;
}

bool readFitResultFile(const TString &filePath, BinResult &result)
{
  static const std::regex kResultFilePattern(
      R"(fit2d_result_y([0-9]+\.[0-9]+)_([0-9]+\.[0-9]+)_pt([0-9]+\.[0-9]+)_([0-9]+\.[0-9]+)\.root)");

  std::cmatch fileNameMatch;
  const std::string fileName = gSystem->BaseName(filePath.Data());
  if (!std::regex_match(fileName.c_str(), fileNameMatch, kResultFilePattern))
    return false;

  std::unique_ptr<TFile> inputFile(TFile::Open(filePath));
  if (!inputFile || inputFile->IsZombie())
    return false;

  double yieldValue = 0.0;
  double yieldError = 0.0;
  const bool hasStoredYield = readParameter(*inputFile, "Nsig", yieldValue);
  const bool hasStoredError = readParameter(*inputFile, "NsigErr", yieldError);

  if (!(hasStoredYield && hasStoredError))
  {
    auto *fitResult = dynamic_cast<RooFitResult *>(inputFile->Get("modelResult"));
    if (!fitResult)
      return false;
    auto *signalYieldVariable = dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find("Nsig"));
    if (!signalYieldVariable)
      return false;
    yieldValue = signalYieldVariable->getVal();
    yieldError = signalYieldVariable->getError();
  }

  result.rapidityLow = std::stod(fileNameMatch[1].str());
  result.rapidityHigh = std::stod(fileNameMatch[2].str());
  result.ptLow = std::stod(fileNameMatch[3].str());
  result.ptHigh = std::stod(fileNameMatch[4].str());
  result.yieldValue = yieldValue;
  result.yieldError = yieldError;
  result.sourcePath = filePath;
  return true;
}

void printArray(const char *declaration, const std::vector<double> &values, int decimalPlaces)
{
  std::cout << declaration << "{";
  for (size_t valueIndex = 0; valueIndex < values.size(); ++valueIndex)
  {
    if (valueIndex)
      std::cout << ", ";
    std::cout << Form(("%." + std::to_string(decimalPlaces) + "f").c_str(), values[valueIndex]);
  }
  std::cout << "};\n";
}
} // namespace

void collect_yield_code(double rapidityLow = 1.6,
                        double rapidityHigh = 2.4,
                        TString resultsDirectory = "roots",
                        TString variablePrefix = "fit",
                        TString binCountVariable = "fitBinCount",
                        bool showSources = false)
{
  const TString resolvedResultsDirectory = resolveResultsDirectory(resultsDirectory);
  if (resolvedResultsDirectory.IsNull())
  {
    std::cerr << "[ERROR] results directory not found: " << resultsDirectory << std::endl;
    return;
  }

  std::vector<TString> resultFiles;
  const TString findCommand = Form("find %s -type f -name 'fit2d_result_*.root'", resolvedResultsDirectory.Data());
  std::unique_ptr<TObjArray> outputLines(gSystem->GetFromPipe(findCommand).Tokenize("\n"));
  for (int lineIndex = 0; outputLines && lineIndex < outputLines->GetEntriesFast(); ++lineIndex)
  {
    auto *line = outputLines->At(lineIndex);
    if (!line)
      continue;
    const TString resultFilePath = line->GetName();
    if (!resultFilePath.IsNull())
      resultFiles.push_back(resultFilePath);
  }

  std::vector<BinResult> selectedResults;
  for (const auto &resultFilePath : resultFiles)
  {
    BinResult result{};
    if (!readFitResultFile(resultFilePath, result))
      continue;
    if (std::abs(result.rapidityLow - rapidityLow) > 1e-9 || std::abs(result.rapidityHigh - rapidityHigh) > 1e-9)
      continue;
    selectedResults.push_back(result);
  }

  std::sort(selectedResults.begin(), selectedResults.end(), [](const BinResult &left, const BinResult &right) {
    if (left.ptLow != right.ptLow)
      return left.ptLow < right.ptLow;
    return left.ptHigh < right.ptHigh;
  });

  if (selectedResults.empty())
  {
    std::cerr << Form("[ERROR] no fit2d ROOT results found for y=%.2f-%.2f under %s",
                      rapidityLow, rapidityHigh, resolvedResultsDirectory.Data())
              << std::endl;
    return;
  }

  std::vector<double> ptLowValues, ptHighValues, yieldValues, yieldErrors;
  for (const auto &result : selectedResults)
  {
    ptLowValues.push_back(result.ptLow);
    ptHighValues.push_back(result.ptHigh);
    yieldValues.push_back(result.yieldValue);
    yieldErrors.push_back(result.yieldError);
    if (showSources)
    {
      std::cerr << Form("[SRC] y=%.2f-%.2f pt=%.2f-%.2f <- %s",
                        result.rapidityLow, result.rapidityHigh, result.ptLow, result.ptHigh, result.sourcePath.Data())
                << std::endl;
    }
  }

  std::cout << Form("// y = %.2f-%.2f\n", rapidityLow, rapidityHigh);
  std::cout << Form("const int %s = %zu;\n", binCountVariable.Data(), selectedResults.size());
  printArray(Form("double %sPtLow[%s] = ", variablePrefix.Data(), binCountVariable.Data()), ptLowValues, 1);
  printArray(Form("double %sPtHigh[%s] = ", variablePrefix.Data(), binCountVariable.Data()), ptHighValues, 1);
  std::cout << "\n";
  printArray(Form("double %sYields[%s] = ", variablePrefix.Data(), binCountVariable.Data()), yieldValues, 5);
  printArray(Form("double %sStatErrors[%s] = ", variablePrefix.Data(), binCountVariable.Data()), yieldErrors, 5);
  std::cout << "\n";
  std::cout << Form("double %sSysErrors[%s];\n", variablePrefix.Data(), binCountVariable.Data());
  std::cout << Form("for (int binIndex = 0; binIndex < %s; ++binIndex)\n", binCountVariable.Data());
  std::cout << Form("  %sSysErrors[binIndex] = 0.0; // sys=0\n\n", variablePrefix.Data());
  std::cout << Form("double %sPtCenters[%s], %sPtHalfWidths[%s];\n", variablePrefix.Data(), binCountVariable.Data(), variablePrefix.Data(), binCountVariable.Data());
  std::cout << Form("for (int binIndex = 0; binIndex < %s; ++binIndex)\n", binCountVariable.Data());
  std::cout << "{\n";
  std::cout << Form("  const double binWidth = %sPtHigh[binIndex] - %sPtLow[binIndex];\n", variablePrefix.Data(), variablePrefix.Data());
  std::cout << Form("  %sPtCenters[binIndex] = 0.5 * (%sPtLow[binIndex] + %sPtHigh[binIndex]);\n", variablePrefix.Data(), variablePrefix.Data(), variablePrefix.Data());
  std::cout << Form("  %sPtHalfWidths[binIndex] = 0.5 * (%sPtHigh[binIndex] - %sPtLow[binIndex]);\n", variablePrefix.Data(), variablePrefix.Data(), variablePrefix.Data());
  std::cout << Form("  %sYields[binIndex] /= binWidth;\n", variablePrefix.Data());
  std::cout << Form("  %sStatErrors[binIndex] /= binWidth;\n", variablePrefix.Data());
  std::cout << Form("  %sSysErrors[binIndex] /= binWidth;\n", variablePrefix.Data());
  std::cout << "}\n";
}
