#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
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

SumResult Integrate2D(const TH2D *num, const TH2D *den,
                      double ptLow, double ptHigh,
                      double centLowPct, double centHighPct)
{
  if (!num || !den)
    throw std::runtime_error("null histogram pointer");

  const double centLow = 2.0 * centLowPct;
  const double centHigh = 2.0 * centHighPct;

  const auto ptBins = FindBinRange(num->GetXaxis(), ptLow, ptHigh);
  const auto centBins = FindBinRange(num->GetYaxis(), centLow, centHigh);

  SumResult out;
  for (int ix = ptBins.first; ix <= ptBins.second; ++ix)
  {
    for (int iy = centBins.first; iy <= centBins.second; ++iy)
    {
      out.num += num->GetBinContent(ix, iy);
      out.numErr2 += std::pow(num->GetBinError(ix, iy), 2);
      out.den += den->GetBinContent(ix, iy);
      out.denErr2 += std::pow(den->GetBinError(ix, iy), 2);
    }
  }
  return out;
}

SumResult Integrate1D(const TH1D *num, const TH1D *den,
                     double ptLow, double ptHigh)
{
  if (!num || !den)
    throw std::runtime_error("null histogram pointer");

  const auto ptBins = FindBinRange(num->GetXaxis(), ptLow, ptHigh);

  SumResult out;
  for (int ix = ptBins.first; ix <= ptBins.second; ++ix)
  {
    out.num += num->GetBinContent(ix);
    out.numErr2 += std::pow(num->GetBinError(ix), 2);
    out.den += den->GetBinContent(ix);
    out.denErr2 += std::pow(den->GetBinError(ix), 2);
  }
  return out;
}

void WriteResult(std::ostream &os, const Query &q, const SumResult &s)
{
  double eff = 0.0;
  double err = 0.0;
  if (s.den > 0.0)
  {
    eff = s.num / s.den;
    if (s.num > 0.0)
      err = eff * std::sqrt(s.numErr2 / (s.num * s.num) + s.denErr2 / (s.den * s.den));
    else
      err = 0.0;
  }
  else
  {
    eff = std::numeric_limits<double>::quiet_NaN();
    err = std::numeric_limits<double>::quiet_NaN();
  }

  os << std::fixed << std::setprecision(1)
     << q.ptLow << "~" << q.ptHigh << "\t"
     << q.centLowPct << "~" << q.centHighPct << "\t"
     << q.absY << "\t";
  os << std::setprecision(6)
     << eff << "\t" << err << "\t"
     << s.num << "\t" << s.den << "\n";
}
} // namespace

void extract_efficiency_table(const char *inputPath =
    "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/skim_roots/eff_PbPb2018_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
    const char *outputPath =
    "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/efficiency_table.tsv",
    bool verbose = false)
{
  std::unique_ptr<TFile> fin(TFile::Open(inputPath, "READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] failed to open " << inputPath << "\n";
    return;
  }

  auto *numMid = dynamic_cast<TH2D *>(fin->Get("hist_eff2d_num_mid_fine_ptw"));
  auto *denMid = dynamic_cast<TH2D *>(fin->Get("hist_eff2d_den_mid_fine_ptw"));
  auto *numFwd = dynamic_cast<TH2D *>(fin->Get("hist_eff2d_num_fwd_fine_ptw"));
  auto *denFwd = dynamic_cast<TH2D *>(fin->Get("hist_eff2d_den_fwd_fine_ptw"));
  auto *numMid1D = dynamic_cast<TH1D *>(fin->Get("hist_eff_num_mid"));
  auto *denMid1D = dynamic_cast<TH1D *>(fin->Get("hist_eff_den_mid"));
  auto *numFwd1D = dynamic_cast<TH1D *>(fin->Get("hist_eff_num_fwd"));
  auto *denFwd1D = dynamic_cast<TH1D *>(fin->Get("hist_eff_den_fwd"));

  const bool hasFineCent = (numMid && denMid && numFwd && denFwd);
  const bool has1D = (numMid1D && denMid1D && numFwd1D && denFwd1D);
  if (!hasFineCent && !has1D)
  {
    std::cerr << "[ERROR] required histograms are missing in " << inputPath << "\n";
    return;
  }

  const std::vector<Query> queries = {
      {3.5, 6.5, 0.0, 90.0, "1.6~2.4", "fwd"},
      {6.5, 9.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {9.0, 12.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {12.0, 40.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 0.0, 10.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 10.0, 20.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 20.0, 30.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 30.0, 40.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 40.0, 50.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 10.0, 30.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 30.0, 50.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 50.0, 90.0, "1.6~2.4", "fwd"},
      {3.5, 40.0, 0.0, 90.0, "1.6~2.4", "fwd"},
      {6.5, 9.0, 0.0, 90.0, "0~1.6", "mid"},
      {9.0, 12.0, 0.0, 90.0, "0~1.6", "mid"},
      {12.0, 15.0, 0.0, 90.0, "0~1.6", "mid"},
      {15.0, 20.0, 0.0, 90.0, "0~1.6", "mid"},
      {20.0, 25.0, 0.0, 90.0, "0~1.6", "mid"},
      {25.0, 30.0, 0.0, 90.0, "0~1.6", "mid"},
      {30.0, 40.0, 0.0, 90.0, "0~1.6", "mid"},
      {25.0, 40.0, 0.0, 90.0, "0~1.6", "mid"},
      {6.5, 40.0, 0.0, 5.0, "0~1.6", "mid"},
      {6.5, 40.0, 5.0, 10.0, "0~1.6", "mid"},
      {6.5, 40.0, 10.0, 15.0, "0~1.6", "mid"},
      {6.5, 40.0, 15.0, 20.0, "0~1.6", "mid"},
      {6.5, 40.0, 0.0, 10.0, "0~1.6", "mid"},
      {6.5, 40.0, 10.0, 20.0, "0~1.6", "mid"},
      {6.5, 40.0, 20.0, 30.0, "0~1.6", "mid"},
      {6.5, 40.0, 30.0, 40.0, "0~1.6", "mid"},
      {6.5, 40.0, 40.0, 50.0, "0~1.6", "mid"},
      {6.5, 40.0, 50.0, 90.0, "0~1.6", "mid"},
      {6.5, 40.0, 0.0, 90.0, "0~1.6", "mid"},
  };

  std::ofstream fout(outputPath);
  if (!fout)
  {
    std::cerr << "[ERROR] failed to open output tsv: " << outputPath << "\n";
    return;
  }

  const char *header = "pt\tcentrality_percent\tabsy\teff\terr\tnum\tden\n";
  fout << header;
  if (verbose)
    std::cout << header;
  for (const auto &q : queries)
  {
    const TH2D *num = (q.source == "mid") ? numMid : numFwd;
    const TH2D *den = (q.source == "mid") ? denMid : denFwd;
    try
    {
      SumResult s;
      if (hasFineCent)
      {
        s = Integrate2D(num, den, q.ptLow, q.ptHigh, q.centLowPct, q.centHighPct);
      }
      else
      {
        if (!NearlyEqual(q.centLowPct, 0.0) || !NearlyEqual(q.centHighPct, 90.0))
          continue;
        const TH1D *num1D = (q.source == "mid") ? numMid1D : numFwd1D;
        const TH1D *den1D = (q.source == "mid") ? denMid1D : denFwd1D;
        s = Integrate1D(num1D, den1D, q.ptLow, q.ptHigh);
      }
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
    std::cout << "[INFO] wrote TSV: " << outputPath << "\n";
  if (!hasFineCent)
  {
    if (verbose)
      std::cout << "[INFO] pp-style fallback used: only centrality 0~90 rows are exported.\n";
  }
}

void extract_efficiency_table_pbpb_all()
{
  extract_efficiency_table(
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/skim_roots/eff_PbPb2018_isMC1_PR_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/efficiency_table_pbpb_PR.tsv",
      false);
  extract_efficiency_table(
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/skim_roots/eff_PbPb2018_isMC1_NP_ncollW1_genW1_ptW1_tnpW1_Dimuon_MiniAOD.root",
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/efficiency_table_pbpb_NP.tsv",
      false);
}
