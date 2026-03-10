#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1D.h"
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
struct FitBin
{
  double yLow = 0.0;
  double yHigh = 0.0;
  double ptLow = 0.0;
  double ptHigh = 0.0;
  double nSig = 0.0;
  double nSigErr = 0.0;
  double bFrac = 0.0;
  double bFracErr = 0.0;
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

bool readFitResultFile(const TString &path, FitBin &out)
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

  double nSig = 0.0;
  double nSigErr = 0.0;
  double bFrac = 0.0;
  double bFracErr = 0.0;

  bool hasNSig = readParam(*f, "Nsig", nSig);
  bool hasNSigErr = readParam(*f, "Nsig_err", nSigErr);
  const bool hasBFrac = readParam(*f, "bFraction", bFrac);
  const bool hasBFracErr = readParam(*f, "bFractionErr", bFracErr);

  if (!hasNSigErr)
  {
    TString massPath = path;
    massPath.ReplaceAll("/fit2d/fit2d_result_", "/mass/mass_model_");
    std::unique_ptr<TFile> fMass(TFile::Open(massPath));
    if (fMass && !fMass->IsZombie())
    {
      hasNSig = hasNSig || readParam(*fMass, "Nsig", nSig);
      hasNSigErr = readParam(*fMass, "Nsig_err", nSigErr);
    }
  }

  if (!(hasNSig && hasNSigErr && hasBFrac && hasBFracErr))
  {
    auto *fitResult = dynamic_cast<RooFitResult *>(f->Get("modelResult"));
    if (!fitResult)
      return false;

    auto *nSigVar = dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find("nSig"));
    auto *bFracVar = dynamic_cast<RooRealVar *>(fitResult->floatParsFinal().find("bFraction"));
    if (!nSigVar || !bFracVar)
      return false;

    nSig = nSigVar->getVal();
    nSigErr = nSigVar->getError();
    bFrac = bFracVar->getVal();
    bFracErr = bFracVar->getError();
  }

  out.yLow = std::stod(match[1].str());
  out.yHigh = std::stod(match[2].str());
  out.ptLow = std::stod(match[3].str());
  out.ptHigh = std::stod(match[4].str());
  out.nSig = nSig;
  out.nSigErr = nSigErr;
  out.bFrac = bFrac;
  out.bFracErr = bFracErr;
  out.source = path;
  return true;
}

std::vector<FitBin> collectFitBins(const TString &rootsDir, double yLow, double yHigh)
{
  std::vector<FitBin> selected;
  const TString cmd = Form("find %s -type f -name 'fit2d_result_*.root'", rootsDir.Data());
  std::unique_ptr<TObjArray> lines(gSystem->GetFromPipe(cmd).Tokenize("\n"));
  for (int i = 0; lines && i < lines->GetEntriesFast(); ++i)
  {
    auto *obj = lines->At(i);
    if (!obj)
      continue;

    FitBin bin;
    const TString path = obj->GetName();
    if (!readFitResultFile(path, bin))
      continue;
    if (std::abs(bin.yLow - yLow) > 1e-9 || std::abs(bin.yHigh - yHigh) > 1e-9)
      continue;
    selected.push_back(bin);
  }

  std::sort(selected.begin(), selected.end(), [](const FitBin &a, const FitBin &b) {
    if (a.ptLow != b.ptLow)
      return a.ptLow < b.ptLow;
    return a.ptHigh < b.ptHigh;
  });
  return selected;
}

TH1D *cloneHist(const TH1D *input, const char *name, const char *title)
{
  auto *hist = dynamic_cast<TH1D *>(input->Clone(name));
  hist->SetDirectory(nullptr);
  hist->SetTitle(title);
  return hist;
}

TH1D *buildFrameHist(const TH1D *templ, const char *name, const char *title, double xmin, double xmax)
{
  auto *frame = cloneHist(templ, name, title);
  frame->Reset("ICE");
  frame->SetBins(templ->GetNbinsX(), xmin, xmax);
  return frame;
}

bool validateMatching(const TH1D *templ, const std::vector<FitBin> &bins, const TString &tag)
{
  for (const auto &bin : bins)
  {
    const double center = 0.5 * (bin.ptLow + bin.ptHigh);
    const int idx = templ->GetXaxis()->FindBin(center);
    if (idx < 1 || idx > templ->GetNbinsX())
    {
      std::cerr << Form("[ERROR] %s pt center %.2f outside histogram range", tag.Data(), center) << std::endl;
      return false;
    }

    const double low = templ->GetBinLowEdge(idx);
    const double high = low + templ->GetBinWidth(idx);
    if (std::abs(low - bin.ptLow) > 1e-6 || std::abs(high - bin.ptHigh) > 1e-6)
    {
      std::cerr << Form("[ERROR] %s pt bin mismatch: fit %.2f-%.2f vs corr %.2f-%.2f",
                        tag.Data(), bin.ptLow, bin.ptHigh, low, high)
                << std::endl;
      return false;
    }
  }
  return true;
}

TH1D *buildRawSpeciesHist(const char *name,
                          const char *title,
                          const TH1D *templ,
                          const std::vector<FitBin> &bins,
                          bool isPrompt)
{
  auto *hist = cloneHist(templ, name, title);
  hist->Reset("ICE");

  for (const auto &bin : bins)
  {
    const double center = 0.5 * (bin.ptLow + bin.ptHigh);
    const int idx = hist->FindBin(center);
    if (idx < 1 || idx > hist->GetNbinsX())
      continue;

    const double width = hist->GetBinWidth(idx);
    const double frac = isPrompt ? (1.0 - bin.bFrac) : bin.bFrac;
    const double value = bin.nSig * frac;
    // Y = N * f, with f = b for NP and f = (1-b) for PR.
    // Assuming Cov(N, b) = 0:
    // sigma_Y^2 = (f sigma_N)^2 + (N sigma_f)^2, and sigma_f = sigma_b.
    const double error2 = frac * frac * bin.nSigErr * bin.nSigErr + bin.nSig * bin.nSig * bin.bFracErr * bin.bFracErr;

    hist->SetBinContent(idx, value / width);
    hist->SetBinError(idx, std::sqrt(std::max(0.0, error2)) / width);
  }

  return hist;
}

TH1D *buildCorrectionHist(const char *name,
                          const char *title,
                          const TH1D *accHist,
                          const TH1D *effHist)
{
  auto *hist = cloneHist(accHist, name, title);
  hist->Reset("ICE");

  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    const double acc = accHist->GetBinContent(i);
    const double dAcc = accHist->GetBinError(i);
    const double eff = effHist->GetBinContent(i);
    const double dEff = effHist->GetBinError(i);

    const double corr = acc * eff;
    double err = 0.0;
    if (corr > 0.0 && acc > 0.0 && eff > 0.0)
    {
      // C = A * E, assuming Cov(A, E) = 0:
      // sigma_C^2 = (E sigma_A)^2 + (A sigma_E)^2.
      const double rel2 = (dAcc / acc) * (dAcc / acc) + (dEff / eff) * (dEff / eff);
      err = corr * std::sqrt(std::max(0.0, rel2));
    }

    hist->SetBinContent(i, corr);
    hist->SetBinError(i, err);
  }

  return hist;
}

TH1D *applyCorrection(const char *name,
                      const char *title,
                      const TH1D *rawHist,
                      const TH1D *corrHist)
{
  auto *hist = cloneHist(rawHist, name, title);
  hist->Reset("ICE");

  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    const double raw = rawHist->GetBinContent(i);
    const double dRaw = rawHist->GetBinError(i);
    const double corr = corrHist->GetBinContent(i);
    const double dCorr = corrHist->GetBinError(i);

    if (corr <= 0.0 || raw <= 0.0)
    {
      hist->SetBinContent(i, 0.0);
      hist->SetBinError(i, 0.0);
      continue;
    }

    const double value = raw / corr;
    // Z = R / C, assuming Cov(R, C) = 0:
    // sigma_Z^2 = (sigma_R / C)^2 + (R sigma_C / C^2)^2.
    const double rel2 = (dRaw / raw) * (dRaw / raw) + (dCorr / corr) * (dCorr / corr);
    hist->SetBinContent(i, value);
    hist->SetBinError(i, value * std::sqrt(std::max(0.0, rel2)));
  }

  return hist;
}

void styleHist(TH1D *hist, int color, int marker, double size)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(size);
  hist->SetLineWidth(2);
}

void drawCmsLabels(const TString &rapidityLabel, const TString &species)
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
    tx.SetTextSize(0.04);
    tx.SetTextFont(72);
    tx.DrawLatex(0.19, 0.935, "CMS Internal");
  }
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.03);
    tx.SetTextFont(42);
    tx.DrawLatex(0.19, 0.865, Form("%s, %s", species.Data(), rapidityLabel.Data()));
  }
}

void saveSpeciesPlot(const TString &outPath,
                     TH1D *rawHist,
                     TH1D *corrHist,
                     const TString &rapidityLabel,
                     const TString &species)
{
  TCanvas canv("c_dn_dpt_single", "", 800, 800);
  canv.SetTitle("");

  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.SetLeftMargin(0.13);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  pad1.cd();

  gPad->SetTicks(1, 1);
  gPad->SetLogy();
  gStyle->SetOptStat(0);

  std::unique_ptr<TH1D> frame(buildFrameHist(rawHist, Form("%s_frame", rawHist->GetName()), "", 0.0, rawHist->GetXaxis()->GetXmax()));

  styleHist(rawHist, kBlack, 24, 1.5);
  styleHist(corrHist, kAzure + 2, 20, 1.4);

  frame->SetTitle("");
  frame->GetYaxis()->SetTitle("dN/dp_{T}");
  frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleSize(0.050);
  frame->GetYaxis()->SetTitleSize(0.050);
  frame->GetXaxis()->SetLabelSize(0.042);
  frame->GetYaxis()->SetLabelSize(0.042);
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.25);

  double maxY = 0.0;
  double minPositive = 1.0e30;
  for (TH1D *hist : {rawHist, corrHist})
  {
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
      const double value = hist->GetBinContent(i);
      const double err = hist->GetBinError(i);
      if (value > 0.0)
      {
        maxY = std::max(maxY, value + err);
        minPositive = std::min(minPositive, std::max(1.0e-12, value - err));
        minPositive = std::min(minPositive, value);
      }
    }
  }
  if (maxY <= 0.0)
    maxY = 1.0;
  if (minPositive >= 1.0e29)
    minPositive = 1.0e-3;

  frame->SetMinimum(minPositive * 0.5);
  frame->SetMaximum(maxY * 5.0);
  frame->Draw();
  rawHist->Draw("E1 SAME");
  corrHist->Draw("E1 SAME");

  drawCmsLabels(rapidityLabel, species);

  TLegend leg(0.60, 0.74, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.028);
  leg.SetEntrySeparation(0.015);
  leg.AddEntry(rawHist, "Raw", "lp");
  leg.AddEntry(corrHist, "Corrected", "lp");
  leg.Draw();
  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");
}
} // namespace

void dn_dpt()
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  const TString baseDir =
      "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim";
  const TString fitRootsDir =
      "/data/users/pjgwak/work/daily_code_tracker/2026/03/09_OO_fit_2d/roots";
  const TString corrDir = baseDir + "/acc_eff_skim/skim_roots";
  const TString outDir = baseDir + "/dn_dpt/outputs";
  const TString figsDir = outDir + "/figs";
  const TString rootsDir = outDir + "/roots";

  const TString accPrPath = corrDir + "/acc_OO2025_isMC1_PR_Dimuon_MiniAOD_Private_MC.root";
  const TString accNpPath = corrDir + "/acc_OO2025_isMC1_NP_Dimuon_MiniAOD_Private_MC.root";
  const TString effPrPath = corrDir + "/eff_OO2025_isMC1_PR_ptW0_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root";
  const TString effNpPath = corrDir + "/eff_OO2025_isMC1_NP_ptW0_Dimuon_MiniAOD_PromptReco_v1_Jul12_merged.root";

  gSystem->mkdir(figsDir, true);
  gSystem->mkdir(rootsDir, true);

  std::unique_ptr<TFile> fAccPr(TFile::Open(accPrPath));
  std::unique_ptr<TFile> fAccNp(TFile::Open(accNpPath));
  std::unique_ptr<TFile> fEffPr(TFile::Open(effPrPath));
  std::unique_ptr<TFile> fEffNp(TFile::Open(effNpPath));
  if (!fAccPr || fAccPr->IsZombie() || !fAccNp || fAccNp->IsZombie() ||
      !fEffPr || fEffPr->IsZombie() || !fEffNp || fEffNp->IsZombie())
  {
    std::cerr << "[ERROR] failed to open one of the acc/eff input ROOT files" << std::endl;
    return;
  }

  auto *accPrMidIn = dynamic_cast<TH1D *>(fAccPr->Get("hist_acc_mid"));
  auto *accPrFwdIn = dynamic_cast<TH1D *>(fAccPr->Get("hist_acc_fwd"));
  auto *accNpMidIn = dynamic_cast<TH1D *>(fAccNp->Get("hist_acc_mid"));
  auto *accNpFwdIn = dynamic_cast<TH1D *>(fAccNp->Get("hist_acc_fwd"));
  auto *effPrMidIn = dynamic_cast<TH1D *>(fEffPr->Get("hist_eff_mid"));
  auto *effPrFwdIn = dynamic_cast<TH1D *>(fEffPr->Get("hist_eff_fwd"));
  auto *effNpMidIn = dynamic_cast<TH1D *>(fEffNp->Get("hist_eff_mid"));
  auto *effNpFwdIn = dynamic_cast<TH1D *>(fEffNp->Get("hist_eff_fwd"));
  if (!accPrMidIn || !accPrFwdIn || !accNpMidIn || !accNpFwdIn ||
      !effPrMidIn || !effPrFwdIn || !effNpMidIn || !effNpFwdIn)
  {
    std::cerr << "[ERROR] missing histograms in acc/eff input files" << std::endl;
    return;
  }

  auto *accPrMid = cloneHist(accPrMidIn, "acc_pr_mid", ";p_{T} (GeV/c);Acceptance");
  auto *accPrFwd = cloneHist(accPrFwdIn, "acc_pr_fwd", ";p_{T} (GeV/c);Acceptance");
  auto *accNpMid = cloneHist(accNpMidIn, "acc_np_mid", ";p_{T} (GeV/c);Acceptance");
  auto *accNpFwd = cloneHist(accNpFwdIn, "acc_np_fwd", ";p_{T} (GeV/c);Acceptance");
  auto *effPrMid = cloneHist(effPrMidIn, "eff_pr_mid", ";p_{T} (GeV/c);Efficiency");
  auto *effPrFwd = cloneHist(effPrFwdIn, "eff_pr_fwd", ";p_{T} (GeV/c);Efficiency");
  auto *effNpMid = cloneHist(effNpMidIn, "eff_np_mid", ";p_{T} (GeV/c);Efficiency");
  auto *effNpFwd = cloneHist(effNpFwdIn, "eff_np_fwd", ";p_{T} (GeV/c);Efficiency");

  const auto midBins = collectFitBins(fitRootsDir, 0.0, 1.6);
  const auto fwdBins = collectFitBins(fitRootsDir, 1.6, 2.4);
  if (midBins.empty() || fwdBins.empty())
  {
    std::cerr << "[ERROR] failed to collect fit2d result files" << std::endl;
    return;
  }

  if (!validateMatching(accPrMid, midBins, "mid"))
    return;
  if (!validateMatching(accPrFwd, fwdBins, "fwd"))
    return;

  auto *rawPrMid = buildRawSpeciesHist("raw_pr_mid", ";p_{T} (GeV/c);dN/dp_{T}", accPrMid, midBins, true);
  auto *rawPrFwd = buildRawSpeciesHist("raw_pr_fwd", ";p_{T} (GeV/c);dN/dp_{T}", accPrFwd, fwdBins, true);
  auto *rawNpMid = buildRawSpeciesHist("raw_np_mid", ";p_{T} (GeV/c);dN/dp_{T}", accNpMid, midBins, false);
  auto *rawNpFwd = buildRawSpeciesHist("raw_np_fwd", ";p_{T} (GeV/c);dN/dp_{T}", accNpFwd, fwdBins, false);

  auto *corrPrMid = buildCorrectionHist("corr_pr_mid", ";p_{T} (GeV/c);Acc#timesEff", accPrMid, effPrMid);
  auto *corrPrFwd = buildCorrectionHist("corr_pr_fwd", ";p_{T} (GeV/c);Acc#timesEff", accPrFwd, effPrFwd);
  auto *corrNpMid = buildCorrectionHist("corr_np_mid", ";p_{T} (GeV/c);Acc#timesEff", accNpMid, effNpMid);
  auto *corrNpFwd = buildCorrectionHist("corr_np_fwd", ";p_{T} (GeV/c);Acc#timesEff", accNpFwd, effNpFwd);

  auto *corrYieldPrMid = applyCorrection("corr_yield_pr_mid", ";p_{T} (GeV/c);dN/dp_{T}", rawPrMid, corrPrMid);
  auto *corrYieldPrFwd = applyCorrection("corr_yield_pr_fwd", ";p_{T} (GeV/c);dN/dp_{T}", rawPrFwd, corrPrFwd);
  auto *corrYieldNpMid = applyCorrection("corr_yield_np_mid", ";p_{T} (GeV/c);dN/dp_{T}", rawNpMid, corrNpMid);
  auto *corrYieldNpFwd = applyCorrection("corr_yield_np_fwd", ";p_{T} (GeV/c);dN/dp_{T}", rawNpFwd, corrNpFwd);

  saveSpeciesPlot(figsDir + "/dn_dpt_pr_mid", rawPrMid, corrYieldPrMid, "|y| < 1.6", "Prompt");
  saveSpeciesPlot(figsDir + "/dn_dpt_np_mid", rawNpMid, corrYieldNpMid, "|y| < 1.6", "Nonprompt");
  saveSpeciesPlot(figsDir + "/dn_dpt_pr_fwd", rawPrFwd, corrYieldPrFwd, "1.6 < |y| < 2.4", "Prompt");
  saveSpeciesPlot(figsDir + "/dn_dpt_np_fwd", rawNpFwd, corrYieldNpFwd, "1.6 < |y| < 2.4", "Nonprompt");

  std::unique_ptr<TFile> fout(TFile::Open(rootsDir + "/dn_dpt.root", "RECREATE"));
  if (!fout || fout->IsZombie())
  {
    std::cerr << "[ERROR] failed to create output ROOT file" << std::endl;
    return;
  }

  rawPrMid->Write();
  rawPrFwd->Write();
  rawNpMid->Write();
  rawNpFwd->Write();
  corrPrMid->Write();
  corrPrFwd->Write();
  corrNpMid->Write();
  corrNpFwd->Write();
  corrYieldPrMid->Write();
  corrYieldPrFwd->Write();
  corrYieldNpMid->Write();
  corrYieldNpFwd->Write();

  accPrMid->Write();
  accPrFwd->Write();
  accNpMid->Write();
  accNpFwd->Write();
  effPrMid->Write();
  effPrFwd->Write();
  effNpMid->Write();
  effNpFwd->Write();

  std::cout << "[INFO] saved output ROOT file: " << rootsDir + "/dn_dpt.root" << std::endl;
  std::cout << "[INFO] saved figures under: " << figsDir << std::endl;
  std::cout << "[INFO] note: error propagation assumes zero covariance between Nsig and bFraction, and between acceptance and efficiency." << std::endl;
}
