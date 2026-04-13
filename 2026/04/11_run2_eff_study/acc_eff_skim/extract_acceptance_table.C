#include "TFile.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
TString resolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}

struct Query
{
  double ptLow;
  double ptHigh;
  double centLowPct;
  double centHighPct;
  std::string absY;
  std::string source;
};

struct SumResult
{
  double num = 0.0;
  double numErr2 = 0.0;
  double den = 0.0;
  double denErr2 = 0.0;
};

constexpr double kTol = 1e-6;

bool NearlyEqual(double a, double b)
{
  return std::abs(a - b) < kTol;
}

std::pair<int, int> FindBinRange(const TAxis *axis, double low, double high)
{
  if (!axis)
    throw std::runtime_error("null axis");

  int first = -1;
  int last = -1;
  for (int ibin = 1; ibin <= axis->GetNbins(); ++ibin)
  {
    const double binLow = axis->GetBinLowEdge(ibin);
    const double binHigh = axis->GetBinUpEdge(ibin);
    if (NearlyEqual(binLow, low))
      first = ibin;
    if (NearlyEqual(binHigh, high))
      last = ibin;
  }

  if (first < 0 || last < 0 || first > last)
    throw std::runtime_error(Form("requested range [%g, %g] does not match histogram bin edges", low, high));

  return {first, last};
}

SumResult Integrate1D(const TH1D *num, const TH1D *den, double low, double high)
{
  if (!num || !den)
    throw std::runtime_error("null histogram pointer");

  const auto bins = FindBinRange(num->GetXaxis(), low, high);
  SumResult out;
  for (int i = bins.first; i <= bins.second; ++i)
  {
    out.num += num->GetBinContent(i);
    out.numErr2 += std::pow(num->GetBinError(i), 2);
    out.den += den->GetBinContent(i);
    out.denErr2 += std::pow(den->GetBinError(i), 2);
  }
  return out;
}

void WriteResult(std::ostream &os, const Query &q, const SumResult &s)
{
  double acc = std::numeric_limits<double>::quiet_NaN();
  double err = std::numeric_limits<double>::quiet_NaN();
  if (s.den > 0.0)
  {
    acc = s.num / s.den;
    if (s.num > 0.0)
      err = acc * std::sqrt(s.numErr2 / (s.num * s.num) + s.denErr2 / (s.den * s.den));
    else
      err = 0.0;
  }

  os << std::fixed << std::setprecision(1)
     << q.ptLow << "~" << q.ptHigh << "\t"
     << q.centLowPct << "~" << q.centHighPct << "\t"
     << q.absY << "\t";
  os << std::setprecision(6)
     << acc << "\t" << err << "\t"
     << s.num << "\t" << s.den << "\n";
}
} // namespace

void extract_acceptance_table(const char *inputPath =
    "",
    const char *outputPath =
    "",
    bool verbose = false)
{
  const TString baseDir = resolveBaseDir();
  const TString tableDir = baseDir + "/outputs/tables";
  const TString inputDir = baseDir + "/skim_roots";
  gSystem->mkdir(tableDir, true);
  const TString inputPathResolved =
      (inputPath && inputPath[0] != '\0')
          ? TString(inputPath)
          : inputDir + "/acc_PbPb2018_ppInput_isMC1_PR_ncollW0_genW1_ptW1.root";
  const TString outputPathResolved =
      (outputPath && outputPath[0] != '\0')
          ? TString(outputPath)
          : tableDir + "/acceptance_table_pbpb.tsv";

  std::unique_ptr<TFile> fin(TFile::Open(inputPathResolved, "READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] failed to open " << inputPathResolved << "\n";
    return;
  }

  auto *numMid = dynamic_cast<TH1D *>(fin->Get("hist_acc_num_mid"));
  auto *denMid = dynamic_cast<TH1D *>(fin->Get("hist_acc_den_mid"));
  auto *numFwd = dynamic_cast<TH1D *>(fin->Get("hist_acc_num_fwd"));
  auto *denFwd = dynamic_cast<TH1D *>(fin->Get("hist_acc_den_fwd"));
  if (!numMid || !denMid || !numFwd || !denFwd)
  {
    std::cerr << "[ERROR] required acceptance histograms are missing in " << inputPathResolved << "\n";
    return;
  }

  const std::vector<Query> queries = {
      {3.5, 6.5, 0.0, 90.0, "1.6~2.4", "fwd"},
      {6.5, 9.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {9.0, 12.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {12.0, 40.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 0.0, 10.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 10.0, 30.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 30.0, 50.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 50.0, 90.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {6.5, 9.0, 0.0, 90.0, "0~1.6", "mid"},
      {9.0, 12.0, 0.0, 90.0, "0~1.6", "mid"},
      {12.0, 15.0, 0.0, 90.0, "0~1.6", "mid"},
      {15.0, 20.0, 0.0, 90.0, "0~1.6", "mid"},
      {20.0, 25.0, 0.0, 90.0, "0~1.6", "mid"},
      {25.0, 40.0, 0.0, 90.0, "0~1.6", "mid"},
      {6.5, 40.0, 0.0, 10.0, "0~1.6", "mid"},
      {6.5, 40.0, 10.0, 20.0, "0~1.6", "mid"},
      {6.5, 40.0, 20.0, 30.0, "0~1.6", "mid"},
      {6.5, 40.0, 30.0, 40.0, "0~1.6", "mid"},
      {6.5, 40.0, 40.0, 50.0, "0~1.6", "mid"},
      {6.5, 40.0, 50.0, 90.0, "0~1.6", "mid"},
      {6.5, 40.0, 0.0, 90.0, "0~1.6", "mid"},
  };

  std::ofstream fout(outputPathResolved.Data());
  if (!fout)
  {
    std::cerr << "[ERROR] failed to open output tsv: " << outputPathResolved << "\n";
    return;
  }

  const char *header = "pt\tcentrality_percent\tabsy\tacc\terr\tnum\tden\n";
  fout << header;
  if (verbose)
    std::cout << header;

  for (const auto &q : queries)
  {
    try
    {
      SumResult s;
      const TH1D *num = (q.source == "mid") ? numMid : numFwd;
      const TH1D *den = (q.source == "mid") ? denMid : denFwd;

      // The current acceptance inputs are derived from pp-based generator
      // samples, so the centrality-sliced acceptance histograms are not
      // meaningful for PbPb rows. Use the integrated acceptance value for
      // every centrality query instead of propagating empty cent bins.
      s = Integrate1D(num, den, q.ptLow, q.ptHigh);

      WriteResult(fout, q, s);
      if (verbose)
        WriteResult(std::cout, q, s);
    }
    catch (const std::exception &ex)
    {
      std::cerr << "[WARN] skipped "
                << q.ptLow << "~" << q.ptHigh << ", "
                << q.centLowPct << "~" << q.centHighPct << ", |y|=" << q.absY
                << " : " << ex.what() << "\n";
    }
  }

  if (verbose)
    std::cout << "[INFO] wrote TSV: " << outputPathResolved << "\n";
}
