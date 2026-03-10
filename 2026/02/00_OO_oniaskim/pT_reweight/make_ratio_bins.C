#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"

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
// ----------------------------------------------------------------------
// Helper container for one (y, pt) bin from the external b-fraction fit.
// ----------------------------------------------------------------------
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

// ----------------------------------------------------------------------
// Read one fit result file and recover the bin edges plus b-fraction value.
// ----------------------------------------------------------------------
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

// ----------------------------------------------------------------------
// Collect all b-fraction fit outputs for one rapidity region.
// ----------------------------------------------------------------------
std::vector<BinResult> collectBFrac(const TString &rootsDir, double yLow, double yHigh)
{
  std::vector<BinResult> selected;
  const TString cmd = Form("find %s -type f -name 'fit2d_result_*.root'", rootsDir.Data());
  std::unique_ptr<TObjArray> lines(gSystem->GetFromPipe(cmd).Tokenize("\n"));
  for (int i = 0; lines && i < lines->GetEntriesFast(); ++i)
  {
    auto *obj = lines->At(i);
    if (!obj)
      continue;

    BinResult result{};
    const TString path = obj->GetName();
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
  return selected;
}

// ----------------------------------------------------------------------
// Map the collected b-fraction values onto the data-count histogram binning.
// ----------------------------------------------------------------------
TH1D *buildBFracHist(const char *name, const char *title, const TH1D *templ, const std::vector<BinResult> &bins)
{
  auto *hist = dynamic_cast<TH1D *>(templ->Clone(name));
  hist->Reset("ICE");
  hist->SetTitle(title);
  hist->SetDirectory(nullptr);

  for (const auto &bin : bins)
  {
    const double center = 0.5 * (bin.ptLow + bin.ptHigh);
    const int idx = hist->FindBin(center);
    if (idx < 1 || idx > hist->GetNbinsX())
      continue;
    hist->SetBinContent(idx, bin.value);
    hist->SetBinError(idx, bin.error);
  }

  return hist;
}

// ----------------------------------------------------------------------
// Build Data NP/PR histograms from the inclusive data histogram and b-fraction.
//
// Error propagation assumes N and b are independent:
//   Y = N * b  ->  sigma_Y^2 = (b sigma_N)^2 + (N sigma_b)^2
//   Y = N * (1-b) -> sigma_Y^2 = ((1-b) sigma_N)^2 + (N sigma_b)^2
//
// If b is extracted from the same data yield, there can be a correlation term:
//   sigma_Y^2 = ... + 2 N b Cov(N, b)
// This covariance is not available here, so the code uses the standard
// zero-covariance approximation.
// ----------------------------------------------------------------------
TH1D *multiplyDataByFraction(const char *name, const char *title, const TH1D *dataHist, const TH1D *fracHist, bool usePromptFraction)
{
  auto *out = dynamic_cast<TH1D *>(dataHist->Clone(name));
  out->SetTitle(title);
  out->SetDirectory(nullptr);

  for (int i = 1; i <= out->GetNbinsX(); ++i)
  {
    const double n = dataHist->GetBinContent(i);
    const double dn = dataHist->GetBinError(i);
    const double b = fracHist->GetBinContent(i);
    const double db = fracHist->GetBinError(i);

    const double frac = usePromptFraction ? (1.0 - b) : b;
    const double dFdB = usePromptFraction ? -1.0 : 1.0;
    const double value = n * frac;
    const double error2 = frac * frac * dn * dn + n * n * dFdB * dFdB * db * db;

    out->SetBinContent(i, value);
    out->SetBinError(i, std::sqrt(std::max(0.0, error2)));
  }

  return out;
}

// ----------------------------------------------------------------------
// Build Data/MC ratio with standard error propagation for a quotient.
// ----------------------------------------------------------------------
TH1D *buildRatioHist(const char *name, const char *title, const TH1D *num, const TH1D *den)
{
  auto *out = dynamic_cast<TH1D *>(num->Clone(name));
  out->SetTitle(title);
  out->SetDirectory(nullptr);
  out->Reset("ICE");

  for (int i = 1; i <= out->GetNbinsX(); ++i)
  {
    const double a = num->GetBinContent(i);
    const double da = num->GetBinError(i);
    const double b = den->GetBinContent(i);
    const double db = den->GetBinError(i);

    if (b <= 0.0 || a <= 0.0)
    {
      out->SetBinContent(i, 0.0);
      out->SetBinError(i, 0.0);
      continue;
    }

    const double ratio = a / b;
    const double rel2 = (da / a) * (da / a) + (db / b) * (db / b);
    out->SetBinContent(i, ratio);
    out->SetBinError(i, ratio * std::sqrt(std::max(0.0, rel2)));
  }

  return out;
}

// ----------------------------------------------------------------------
// Normalize a histogram to unit area in bin-width convention.
// ----------------------------------------------------------------------
TH1D *makeNormalizedHist(const TH1D *input, const char *name, const char *title)
{
  auto *out = dynamic_cast<TH1D *>(input->Clone(name));
  out->SetTitle(title);
  out->SetDirectory(nullptr);
  const double integral = out->Integral(1, out->GetNbinsX());
  if (integral > 0.0)
    out->Scale(1.0 / integral, "width");
  return out;
}

// ----------------------------------------------------------------------
// Apply a common drawing style to all histograms used in the output plots.
// ----------------------------------------------------------------------
void styleHist(TH1D *hist, int color, int marker, double size)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(size);
  hist->SetLineWidth(3);
}

// ----------------------------------------------------------------------
// Draw CMS-style labels following the layout used in the reference macro.
// ----------------------------------------------------------------------
void drawCmsLabels(const TString &kinLabel)
{
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextFont(42);
    tx.SetTextAlign(31);
    tx.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.040);
    tx.SetTextFont(72);
    tx.DrawLatex(0.18, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.030);
    tx.SetTextFont(42);
    tx.DrawLatex(0.18, 0.865, kinLabel);
  }
}

// ----------------------------------------------------------------------
// Save the main comparison figure:
//   top    : Data vs MC
//   bottom : Data/MC ratio
// ----------------------------------------------------------------------
void saveComparisonPlot(const TString &outPath,
                        TH1D *dataHist,
                        TH1D *mcHist,
                        TH1D *ratioHist,
                        const TString &label,
                        const TString &species)
{
  TCanvas canv("c", "", 800, 800);
  canv.cd();

  TPad padTop("padTop", "", 0.0, 0.25, 1.0, 1.0);
  TPad padBot("padBot", "", 0.0, 0.0, 1.0, 0.25);
  padTop.SetTopMargin(0.08);
  padTop.SetBottomMargin(0.00001);
  padTop.SetLeftMargin(0.15);
  padTop.SetRightMargin(0.03);
  padBot.SetTopMargin(0.00001);
  padBot.SetBottomMargin(0.40);
  padBot.SetLeftMargin(0.15);
  padBot.SetRightMargin(0.03);
  padTop.Draw();
  padBot.Draw();

  padTop.cd();
  gPad->SetTicks(1, 1);

  styleHist(dataHist, kBlack, 20, 1.45);
  styleHist(mcHist, kAzure + 2, 24, 1.30);
  dataHist->SetTitle("");
  dataHist->GetYaxis()->SetTitle("Normalized counts");
  dataHist->GetYaxis()->CenterTitle();
  dataHist->GetYaxis()->SetTitleSize(0.055);
  dataHist->GetYaxis()->SetLabelSize(0.042);
  dataHist->GetYaxis()->SetTitleOffset(1.20);
  dataHist->GetXaxis()->SetTitle("");
  dataHist->GetXaxis()->SetLabelSize(0.0);
  dataHist->GetXaxis()->SetTitleSize(0.0);
  double maxY = std::max(dataHist->GetMaximum(), mcHist->GetMaximum());
  if (maxY <= 0.0)
    maxY = 1.0;
  dataHist->SetMinimum(0.0);
  dataHist->SetMaximum(maxY * 1.55);
  dataHist->Draw("E1");
  mcHist->Draw("E1 SAME");

  drawCmsLabels(label);

  TLegend leg(0.58, 0.74, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.028);
  leg.SetEntrySeparation(0.02);
  leg.AddEntry(dataHist, Form("Data %s", species.Data()), "lp");
  leg.AddEntry(mcHist, Form("MC %s", species.Data()), "lp");
  leg.Draw();

  padBot.cd();
  gPad->SetTicks(1, 1);
  styleHist(ratioHist, kBlack, 20, 1.20);
  ratioHist->SetTitle("");
  ratioHist->GetYaxis()->SetTitle("Data/MC");
  ratioHist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ratioHist->GetXaxis()->CenterTitle();
  ratioHist->GetYaxis()->CenterTitle();
  ratioHist->GetYaxis()->SetNdivisions(505);
  ratioHist->GetYaxis()->SetTitleSize(0.12);
  ratioHist->GetYaxis()->SetTitleOffset(0.50);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetLabelSize(0.10);
  double maxRatio = 0.0;
  for (int i = 1; i <= ratioHist->GetNbinsX(); ++i)
  {
    maxRatio = std::max(maxRatio, ratioHist->GetBinContent(i) + ratioHist->GetBinError(i));
  }
  ratioHist->SetMinimum(0.0);
  ratioHist->SetMaximum(maxRatio > 0.0 ? 1.15 * maxRatio : 2.0);
  ratioHist->Draw("E1");

  TLine unity(ratioHist->GetXaxis()->GetXmin(), 1.0, ratioHist->GetXaxis()->GetXmax(), 1.0);
  unity.SetLineColor(kRed + 1);
  unity.SetLineStyle(2);
  unity.SetLineWidth(2);
  unity.Draw("SAME");
  ratioHist->Draw("E1 SAME");

  canv.SaveAs(outPath + ".png");
  canv.SaveAs(outPath + ".pdf");
}

// ----------------------------------------------------------------------
// Check that the fit-result pt binning matches the histogram binning.
// ----------------------------------------------------------------------
bool validateMatching(const TH1D *templ, const std::vector<BinResult> &bins, const TString &tag)
{
  for (const auto &bin : bins)
  {
    const double center = 0.5 * (bin.ptLow + bin.ptHigh);
    const int idx = templ->GetXaxis()->FindBin(center);
    if (idx < 1 || idx > templ->GetNbinsX())
    {
      std::cerr << Form("[ERROR] %s pt bin center %.2f is outside histogram range",
                        tag.Data(), center)
                << std::endl;
      return false;
    }
    const double low = templ->GetBinLowEdge(idx);
    const double high = low + templ->GetBinWidth(idx);
    if (std::abs(low - bin.ptLow) > 1e-6 || std::abs(high - bin.ptHigh) > 1e-6)
    {
      std::cerr << Form("[ERROR] %s pt bin mismatch: bFrac %.2f-%.2f vs hist %.2f-%.2f",
                        tag.Data(), bin.ptLow, bin.ptHigh, low, high)
                << std::endl;
      return false;
    }
  }
  return true;
}
} // namespace

void make_ratio_bins()
{
  // --------------------------------------------------------------------
  // Global ROOT drawing setup.
  // --------------------------------------------------------------------
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  // --------------------------------------------------------------------
  // Input locations.
  // --------------------------------------------------------------------
  const TString baseDir =
      "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/pT_reweight";
  const TString fillRootsDir = baseDir + "/fill_bins_outputs/roots";
  const TString dataCountPath =
      fillRootsDir + "/pt_bins_OO2025_isMC0_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root";
  const TString mcPrPath =
      fillRootsDir + "/pt_bins_OO2025_isMC1_PR_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root";
  const TString mcNpPath =
      fillRootsDir + "/pt_bins_OO2025_isMC1_NP_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root";
  const TString bFracRootsDir =
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/09_OO_fit_2d/roots";

  // --------------------------------------------------------------------
  // Open the inclusive-data and MC template histograms.
  // --------------------------------------------------------------------
  std::unique_ptr<TFile> fData(TFile::Open(dataCountPath));
  std::unique_ptr<TFile> fMcPr(TFile::Open(mcPrPath));
  std::unique_ptr<TFile> fMcNp(TFile::Open(mcNpPath));
  if (!fData || fData->IsZombie() || !fMcPr || fMcPr->IsZombie() || !fMcNp || fMcNp->IsZombie())
  {
    std::cerr << "[ERROR] failed to open one of the input ROOT files" << std::endl;
    return;
  }

  auto *dataMidRaw = dynamic_cast<TH1D *>(fData->Get("hist_raw_mid"));
  auto *dataFwdRaw = dynamic_cast<TH1D *>(fData->Get("hist_raw_fwd"));
  auto *mcPrMidRaw = dynamic_cast<TH1D *>(fMcPr->Get("hist_raw_mid"));
  auto *mcPrFwdRaw = dynamic_cast<TH1D *>(fMcPr->Get("hist_raw_fwd"));
  auto *mcNpMidRaw = dynamic_cast<TH1D *>(fMcNp->Get("hist_raw_mid"));
  auto *mcNpFwdRaw = dynamic_cast<TH1D *>(fMcNp->Get("hist_raw_fwd"));
  if (!dataMidRaw || !dataFwdRaw || !mcPrMidRaw || !mcPrFwdRaw || !mcNpMidRaw || !mcNpFwdRaw)
  {
    std::cerr << "[ERROR] missing raw histograms in input ROOT files. Re-run make_pt_bins.C first." << std::endl;
    return;
  }

  // Detach histograms from input files so they remain valid after file close.
  dataMidRaw->SetDirectory(nullptr);
  dataFwdRaw->SetDirectory(nullptr);
  mcPrMidRaw->SetDirectory(nullptr);
  mcPrFwdRaw->SetDirectory(nullptr);
  mcNpMidRaw->SetDirectory(nullptr);
  mcNpFwdRaw->SetDirectory(nullptr);

  // --------------------------------------------------------------------
  // Collect fitted b-fraction values for mid-rapidity and forward-rapidity.
  // --------------------------------------------------------------------
  const auto bFracMidBins = collectBFrac(bFracRootsDir, 0.0, 1.6);
  const auto bFracFwdBins = collectBFrac(bFracRootsDir, 1.6, 2.4);
  if (bFracMidBins.empty() || bFracFwdBins.empty())
  {
    std::cerr << "[ERROR] failed to collect bFraction results for all rapidity bins" << std::endl;
    return;
  }

  // Make sure the fit bin definitions still line up with the counting histograms.
  if (!validateMatching(dataMidRaw, bFracMidBins, "mid"))
    return;
  if (!validateMatching(dataFwdRaw, bFracFwdBins, "fwd"))
    return;

  // --------------------------------------------------------------------
  // Build b-fraction histograms on the same binning as the counting inputs.
  // --------------------------------------------------------------------
  auto *bFracMid = buildBFracHist("bFrac_mid", ";p_{T} (GeV/c);b fraction", dataMidRaw, bFracMidBins);
  auto *bFracFwd = buildBFracHist("bFrac_fwd", ";p_{T} (GeV/c);b fraction", dataFwdRaw, bFracFwdBins);

  // --------------------------------------------------------------------
  // Split inclusive raw-count histograms into non-prompt and prompt components,
  // then build normalized copies for plotting and Data/MC ratios.
  // --------------------------------------------------------------------
  auto *dataNpMidRaw = multiplyDataByFraction("data_np_mid_raw", ";p_{T} (GeV/c);Counts", dataMidRaw, bFracMid, false);
  auto *dataNpFwdRaw = multiplyDataByFraction("data_np_fwd_raw", ";p_{T} (GeV/c);Counts", dataFwdRaw, bFracFwd, false);
  auto *dataPrMidRaw = multiplyDataByFraction("data_pr_mid_raw", ";p_{T} (GeV/c);Counts", dataMidRaw, bFracMid, true);
  auto *dataPrFwdRaw = multiplyDataByFraction("data_pr_fwd_raw", ";p_{T} (GeV/c);Counts", dataFwdRaw, bFracFwd, true);

  auto *dataNpMid = makeNormalizedHist(dataNpMidRaw, "data_np_mid", ";p_{T} (GeV/c);Normalized counts");
  auto *dataNpFwd = makeNormalizedHist(dataNpFwdRaw, "data_np_fwd", ";p_{T} (GeV/c);Normalized counts");
  auto *dataPrMid = makeNormalizedHist(dataPrMidRaw, "data_pr_mid", ";p_{T} (GeV/c);Normalized counts");
  auto *dataPrFwd = makeNormalizedHist(dataPrFwdRaw, "data_pr_fwd", ";p_{T} (GeV/c);Normalized counts");

  auto *mcNpMid = makeNormalizedHist(mcNpMidRaw, "mc_np_mid", ";p_{T} (GeV/c);Normalized counts");
  auto *mcNpFwd = makeNormalizedHist(mcNpFwdRaw, "mc_np_fwd", ";p_{T} (GeV/c);Normalized counts");
  auto *mcPrMid = makeNormalizedHist(mcPrMidRaw, "mc_pr_mid", ";p_{T} (GeV/c);Normalized counts");
  auto *mcPrFwd = makeNormalizedHist(mcPrFwdRaw, "mc_pr_fwd", ";p_{T} (GeV/c);Normalized counts");

  // --------------------------------------------------------------------
  // Build Data/MC ratios for each species and rapidity region.
  // --------------------------------------------------------------------
  auto *ratioNpMid = buildRatioHist("ratio_np_mid", ";p_{T} (GeV/c);Data/MC", dataNpMid, mcNpMid);
  auto *ratioNpFwd = buildRatioHist("ratio_np_fwd", ";p_{T} (GeV/c);Data/MC", dataNpFwd, mcNpFwd);
  auto *ratioPrMid = buildRatioHist("ratio_pr_mid", ";p_{T} (GeV/c);Data/MC", dataPrMid, mcPrMid);
  auto *ratioPrFwd = buildRatioHist("ratio_pr_fwd", ";p_{T} (GeV/c);Data/MC", dataPrFwd, mcPrFwd);

  // --------------------------------------------------------------------
  // Save histograms to an output ROOT file for later reuse.
  // --------------------------------------------------------------------
  const TString outDir = baseDir + "/ratio_bins_outputs";
  const TString figsDir = outDir + "/figs";
  const TString rootsDir = outDir + "/roots";
  gSystem->mkdir(figsDir, true);
  gSystem->mkdir(rootsDir, true);

  std::unique_ptr<TFile> fout(TFile::Open(rootsDir + "/ratio_bins.root", "RECREATE"));
  bFracMid->Write();
  bFracFwd->Write();
  dataNpMidRaw->Write();
  dataNpFwdRaw->Write();
  dataPrMidRaw->Write();
  dataPrFwdRaw->Write();
  dataNpMid->Write();
  dataNpFwd->Write();
  dataPrMid->Write();
  dataPrFwd->Write();
  mcNpMidRaw->Write("mc_np_mid_raw");
  mcNpFwdRaw->Write("mc_np_fwd_raw");
  mcPrMidRaw->Write("mc_pr_mid_raw");
  mcPrFwdRaw->Write("mc_pr_fwd_raw");
  mcNpMid->Write();
  mcNpFwd->Write();
  mcPrMid->Write();
  mcPrFwd->Write();
  ratioNpMid->Write();
  ratioNpFwd->Write();
  ratioPrMid->Write();
  ratioPrFwd->Write();
  fout->Close();

  // --------------------------------------------------------------------
  // Save final comparison plots for each rapidity/species combination.
  // --------------------------------------------------------------------
  saveComparisonPlot(figsDir + "/np_mid", dataNpMid, mcNpMid, ratioNpMid, "|y| < 1.6", "NP");
  saveComparisonPlot(figsDir + "/np_fwd", dataNpFwd, mcNpFwd, ratioNpFwd, "1.6 < |y| < 2.4", "NP");
  saveComparisonPlot(figsDir + "/pr_mid", dataPrMid, mcPrMid, ratioPrMid, "|y| < 1.6", "PR");
  saveComparisonPlot(figsDir + "/pr_fwd", dataPrFwd, mcPrFwd, ratioPrFwd, "1.6 < |y| < 2.4", "PR");

  std::cout << "[INFO] saved output ROOT file: " << rootsDir + "/ratio_bins.root" << std::endl;
  std::cout << "[INFO] saved plots under: " << figsDir << std::endl;
}
