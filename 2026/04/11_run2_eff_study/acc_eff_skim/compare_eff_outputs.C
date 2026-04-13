#include "TFile.h"
#include "TH1.h"
#include <cmath>
#include <iostream>
#include <string>

void compare_one_hist(TFile *fold, TFile *fnew, const char *name)
{
  auto *hold = dynamic_cast<TH1 *>(fold->Get(name));
  auto *hnew = dynamic_cast<TH1 *>(fnew->Get(name));
  if (!hold || !hnew)
  {
    std::cout << "[ERROR] missing histogram: " << name << "\n";
    return;
  }

  const int nbins = std::min(hold->GetNbinsX(), hnew->GetNbinsX());
  double maxAbsDiff = -1.0;
  int maxBin = -1;
  std::cout << "[HIST] " << name << "\n";
  for (int i = 1; i <= nbins; ++i)
  {
    const double oldv = hold->GetBinContent(i);
    const double newv = hnew->GetBinContent(i);
    const double diff = newv - oldv;
    if (std::fabs(diff) > maxAbsDiff)
    {
      maxAbsDiff = std::fabs(diff);
      maxBin = i;
    }
    std::cout << "  bin " << i
              << ": old=" << oldv
              << ", new=" << newv
              << ", diff=" << diff << "\n";
  }
  std::cout << "  max_abs_diff_bin=" << maxBin
            << ", max_abs_diff=" << maxAbsDiff << "\n";
}

void compare_eff_outputs(const char *oldPath, const char *newPath)
{
  TFile *fold = TFile::Open(oldPath, "READ");
  TFile *fnew = TFile::Open(newPath, "READ");
  if (!fold || fold->IsZombie() || !fnew || fnew->IsZombie())
  {
    std::cout << "[ERROR] failed to open input files\n";
    return;
  }

  std::cout << "[OLD] " << oldPath << "\n";
  std::cout << "[NEW] " << newPath << "\n";
  compare_one_hist(fold, fnew, "hist_eff_mid");
  compare_one_hist(fold, fnew, "hist_eff_fwd");
  compare_one_hist(fold, fnew, "hist_eff_num_mid");
  compare_one_hist(fold, fnew, "hist_eff_num_fwd");

  fold->Close();
  fnew->Close();
}
