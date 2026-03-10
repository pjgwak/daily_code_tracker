#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TPad.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace
{
struct FitConfig
{
  TString histName;
  TString tag;
  TString species;
  TString rapLabel;
};

struct FitSummary
{
  int status = 999;
  int hesse = -1;
  double chi2 = 0.0;
  int ndf = 0;
};

double singleExpWeight(double *x, double *par)
{
  const double pt = x[0];
  const double A = par[0];
  const double B = par[1];
  const double C = par[2];

  if (A < 0.0 || B < 0.0 || C < 0.0)
    return 0.0;

  const double val = A * std::exp(-B * pt) + C;
  return val > 0.0 ? val : 0.0;
}

double doubleExpWeight(double *x, double *par)
{
  const double pt = x[0];
  const double A = par[0];
  const double B = par[1];
  const double D = par[2];
  const double E = par[3];
  const double C = par[4];

  if (A < 0.0 || B < 0.0 || D < 0.0 || E < 0.0 || C < 0.0)
    return 0.0;

  const double val = A * std::exp(-B * pt) + D * std::exp(-E * pt) + C;
  return val > 0.0 ? val : 0.0;
}

void styleHist(TH1D *hist, int color, int marker)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.25);
  hist->SetLineWidth(2);
}

void styleFunc(TF1 *func, int color, int style)
{
  func->SetLineColor(color);
  func->SetLineWidth(3);
  func->SetLineStyle(style);
  func->SetNpx(600);
}

void drawCmsLabels(const TString &kinLabel)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(72);
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.19, 0.93, "CMS Internal");

  tx.SetTextFont(42);
  tx.SetTextSize(0.033);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.93, "OO 5.36 TeV (9 nb^{-1})");
  tx.SetTextAlign(11);
}

TF1 *makeFitFunction(const TString &name, const TH1D *hist, int expTerms)
{
  TF1 *func = nullptr;
  const double ymax = std::max(1.0, hist->GetMaximum());
  const double tail = std::max(0.0, hist->GetBinContent(hist->GetNbinsX()));

  if (expTerms == 1)
  {
    func = new TF1(name, singleExpWeight, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 3);
    func->SetParNames("A", "B", "C");
    func->SetParameters(std::max(0.1, ymax - tail), 0.2, std::max(0.1, tail));
    func->SetParLimits(0, 0.0, std::max(10.0, ymax * 10.0));
    func->SetParLimits(1, 1.0e-4, 10.0);
    func->SetParLimits(2, 0.0, std::max(10.0, ymax * 5.0));
  }
  else
  {
    func = new TF1(name, doubleExpWeight, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 5);
    func->SetParNames("A", "B", "D", "E", "C");
    func->SetParameters(std::max(0.1, 0.6 * (ymax - tail)),
                        0.35,
                        std::max(0.1, 0.3 * (ymax - tail)),
                        0.07,
                        std::max(0.1, tail));
    func->SetParLimits(0, 0.0, std::max(10.0, ymax * 10.0));
    func->SetParLimits(1, 1.0e-4, 10.0);
    func->SetParLimits(2, 0.0, std::max(10.0, ymax * 10.0));
    func->SetParLimits(3, 1.0e-4, 10.0);
    func->SetParLimits(4, 0.0, std::max(10.0, ymax * 5.0));
  }

  styleFunc(func, kRed + 1, 1);
  return func;
}

TF1 *fitWithRetries(TH1D *hist, const TString &name, int expTerms, FitSummary &summary)
{
  struct Seed1
  {
    double A;
    double B;
    double C;
  };

  struct Seed2
  {
    double A;
    double B;
    double D;
    double E;
    double C;
  };

  const double ymax = std::max(1.0, hist->GetMaximum());
  const double tail = std::max(0.0, hist->GetBinContent(hist->GetNbinsX()));
  const std::vector<Seed1> seeds1 = {
      {std::max(0.1, 0.4 * (ymax - tail)), 0.10, std::max(0.1, tail)},
      {std::max(0.1, 0.7 * (ymax - tail)), 0.20, std::max(0.1, 0.9 * tail)},
      {std::max(0.1, 1.0 * (ymax - tail)), 0.35, std::max(0.1, 1.1 * tail)},
      {std::max(0.1, 1.3 * (ymax - tail)), 0.60, std::max(0.1, 0.7 * tail)},
  };
  const std::vector<Seed2> seeds2 = {
      {std::max(0.1, 0.50 * (ymax - tail)), 0.35, std::max(0.1, 0.25 * (ymax - tail)), 0.07, std::max(0.1, tail)},
      {std::max(0.1, 0.80 * (ymax - tail)), 0.50, std::max(0.1, 0.30 * (ymax - tail)), 0.10, std::max(0.1, 0.8 * tail)},
      {std::max(0.1, 0.35 * (ymax - tail)), 0.20, std::max(0.1, 0.45 * (ymax - tail)), 0.03, std::max(0.1, 1.1 * tail)},
      {std::max(0.1, 1.00 * (ymax - tail)), 0.80, std::max(0.1, 0.20 * (ymax - tail)), 0.05, std::max(0.1, 0.6 * tail)},
  };

  std::unique_ptr<TF1> best;
  double bestMetric = std::numeric_limits<double>::infinity();

  const size_t nSeeds = (expTerms == 1) ? seeds1.size() : seeds2.size();
  for (size_t i = 0; i < nSeeds; ++i)
  {
    std::unique_ptr<TF1> trial(makeFitFunction(Form("%s_try%zu", name.Data(), i), hist, expTerms));
    if (expTerms == 1)
    {
      trial->SetParameter(0, seeds1[i].A);
      trial->SetParameter(1, seeds1[i].B);
      trial->SetParameter(2, seeds1[i].C);
    }
    else
    {
      trial->SetParameter(0, seeds2[i].A);
      trial->SetParameter(1, seeds2[i].B);
      trial->SetParameter(2, seeds2[i].D);
      trial->SetParameter(3, seeds2[i].E);
      trial->SetParameter(4, seeds2[i].C);
    }

    const TFitResultPtr res = hist->Fit(trial.get(), "ISRQ0");
    const int status = static_cast<int>(res);
    const int hesse = res.Get() ? res->CovMatrixStatus() : -1;
    const int ndf = trial->GetNDF();
    const double chi2 = trial->GetChisquare();
    const double metric = (status == 0 && hesse >= 2 && ndf > 0) ? chi2 / ndf : 1.0e12 + i;

    if (metric < bestMetric)
    {
      bestMetric = metric;
      best.reset(static_cast<TF1 *>(trial->Clone(name)));
      summary.status = status;
      summary.hesse = hesse;
      summary.chi2 = chi2;
      summary.ndf = ndf;
    }
  }

  if (!best)
  {
    best.reset(makeFitFunction(name, hist, expTerms));
    summary.status = 999;
    summary.hesse = -1;
    summary.chi2 = 0.0;
    summary.ndf = 0;
  }

  styleFunc(best.get(), kRed + 1, 1);
  return static_cast<TF1 *>(best.release());
}

TH1D *buildFrameHist(const TH1D *templ, const TString &name)
{
  auto *frame = static_cast<TH1D *>(templ->Clone(name));
  frame->Reset("ICES");
  frame->SetDirectory(nullptr);
  frame->SetTitle("");
  return frame;
}

void saveFitPlot(const TString &outPath,
                 TH1D *hist,
                 TF1 *func,
                 const FitConfig &cfg,
                 const FitSummary &summary,
                 int expTerms)
{
  TCanvas canv("c_fit", "", 800, 820);
  canv.cd();

  TPad padTop("padTop", "", 0.0, 0.25, 1.0, 1.0);
  TPad padBot("padBot", "", 0.0, 0.0, 1.0, 0.25);
  padTop.SetLeftMargin(0.15);
  padTop.SetRightMargin(0.04);
  padTop.SetTopMargin(0.08);
  padTop.SetBottomMargin(0.00001);
  padBot.SetLeftMargin(0.15);
  padBot.SetRightMargin(0.04);
  padBot.SetTopMargin(0.00001);
  padBot.SetBottomMargin(0.4);
  padTop.Draw();
  padBot.Draw();

  padTop.cd();
  gPad->SetTicks(1, 1);

  styleHist(hist, kBlack, 20);
  hist->GetYaxis()->SetTitle("Data/MC");
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->GetXaxis()->SetLabelSize(0.0);
  hist->GetXaxis()->SetTitleSize(0.0);
  hist->SetMinimum(0.0);
  hist->SetMaximum(std::max(1.4, hist->GetMaximum() * 1.35));
  hist->Draw("E1");
  func->Draw("SAME");

  drawCmsLabels(cfg.rapLabel);

  TLegend leg(0.50, 0.70, 0.72, 0.89);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030);
  leg.AddEntry(hist, Form("%s ratio", cfg.species.Data()), "lp");
  leg.AddEntry(func, expTerms == 1 ? "1-exp fit" : "2-exp fit", "l");
  leg.Draw();

  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.028);
    tx.SetTextFont(42);
    double xtext = 0.205, y0 = 0.865, dy = -0.05;
    int k = 0;
    tx.DrawLatex(xtext, y0 + dy * k++, Form("%s %s", cfg.species.Data(), cfg.rapLabel.Data()));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("Status : MINIMIZE=%d HESSE=%d", summary.status, summary.hesse));
    tx.DrawLatex(xtext, y0 + dy * k++, Form("#chi^{2}/ndf = %.2f / %d", summary.chi2, summary.ndf));
    tx.DrawLatex(xtext, y0 + dy * k++, expTerms == 1 ? "f = A e^{-Bp_{T}} + C" : "f = A e^{-Bp_{T}} + D e^{-Ep_{T}} + C");
  }
  {
    TLatex tp;
    tp.SetNDC();
    tp.SetTextSize(0.028);
    tp.SetTextFont(42);
    double xtext = 0.70, y0 = 0.865, dy = -0.040;
    int k = 0;
    for (int ipar = 0; ipar < func->GetNpar(); ++ipar)
    {
      tp.DrawLatex(xtext,
                   y0 + dy * k++,
                   Form("%s = %.3f #pm %.3f",
                        func->GetParName(ipar),
                        func->GetParameter(ipar),
                        func->GetParError(ipar)));
    }
  }

  padBot.cd();
  gPad->SetTicks(1, 1);

  auto *ratioHist = buildFrameHist(hist, Form("%s_fitOverData", hist->GetName()));
  ratioHist->SetTitle("");
  ratioHist->GetYaxis()->SetTitle("Fit/Data");
  ratioHist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ratioHist->GetYaxis()->CenterTitle();
  ratioHist->GetXaxis()->CenterTitle();
  ratioHist->GetYaxis()->SetNdivisions(505);
  ratioHist->GetYaxis()->SetTitleSize(0.12);
  ratioHist->GetYaxis()->SetTitleOffset(0.4);
  ratioHist->GetYaxis()->SetLabelSize(0.10);
  ratioHist->GetXaxis()->SetTitleSize(0.15);
  ratioHist->GetXaxis()->SetLabelSize(0.10);

  double maxDev = 1.1;
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    const double x = hist->GetBinCenter(i);
    const double y = hist->GetBinContent(i);
    const double ey = hist->GetBinError(i);
    if (y <= 0.0)
      continue;
    const double fitVal = func->Eval(x);
    const double ratio = fitVal / y;
    const double eratio = ey > 0.0 ? ratio * ey / y : 0.0;
    ratioHist->SetBinContent(i, ratio);
    ratioHist->SetBinError(i, eratio);
    maxDev = std::max(maxDev, ratio + eratio);
  }
  ratioHist->SetMinimum(0.6);
  ratioHist->SetMaximum(std::max(1.4, maxDev * 1.05));
  styleHist(ratioHist, kBlack, 20);
  ratioHist->Draw("E1");

  TLine line(hist->GetXaxis()->GetXmin(), 1.0, hist->GetXaxis()->GetXmax(), 1.0);
  line.SetLineColor(kBlue + 1);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw("SAME");
  ratioHist->Draw("E1 SAME");

  canv.SaveAs(outPath + ".pdf");
  canv.SaveAs(outPath + ".png");

  delete ratioHist;
}
} // namespace

void fit_ratio(int expTerms = 2)
{
  gROOT->SetBatch(kTRUE);
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if (expTerms != 1 && expTerms != 2)
  {
    std::cerr << "[ERROR] expTerms must be 1 or 2, got " << expTerms << std::endl;
    return;
  }

  const TString baseDir =
      "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/pT_reweight";
  const TString inPath = baseDir + "/ratio_bins_outputs/roots/ratio_bins.root";
  const TString outDir = baseDir + "/fit_ratio_outputs";
  const TString figsDir = outDir + "/figs";
  const TString rootsDir = outDir + "/roots";

  gSystem->mkdir(figsDir, true);
  gSystem->mkdir(rootsDir, true);

  std::unique_ptr<TFile> fin(TFile::Open(inPath));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] failed to open input ROOT file: " << inPath << std::endl;
    return;
  }

  const std::vector<FitConfig> configs = {
      {"ratio_np_mid", "np_mid", "NP", "|y| < 1.6"},
      {"ratio_np_fwd", "np_fwd", "NP", "1.6 < |y| < 2.4"},
      {"ratio_pr_mid", "pr_mid", "PR", "|y| < 1.6"},
      {"ratio_pr_fwd", "pr_fwd", "PR", "1.6 < |y| < 2.4"},
  };

  const TString expTag = Form("exp%d", expTerms);
  std::unique_ptr<TFile> fout(TFile::Open(rootsDir + Form("/fit_ratio_%s.root", expTag.Data()), "RECREATE"));
  if (!fout || fout->IsZombie())
  {
    std::cerr << "[ERROR] failed to create output ROOT file" << std::endl;
    return;
  }

  for (const auto &cfg : configs)
  {
    auto *histIn = dynamic_cast<TH1D *>(fin->Get(cfg.histName));
    if (!histIn)
    {
      std::cerr << "[ERROR] missing histogram: " << cfg.histName << std::endl;
      continue;
    }

    std::unique_ptr<TH1D> hist(static_cast<TH1D *>(histIn->Clone(cfg.histName + "_clone")));
    hist->SetDirectory(nullptr);

    FitSummary summary;
    std::unique_ptr<TF1> func(fitWithRetries(hist.get(), Form("fit_%s_exp%d", cfg.tag.Data(), expTerms), expTerms, summary));

    saveFitPlot(figsDir + "/" + cfg.tag + "_" + expTag, hist.get(), func.get(), cfg, summary, expTerms);

    fout->cd();
    hist->Write(cfg.histName + "_fit_input");
    func->Write(Form("exp%d_%s", expTerms, cfg.tag.Data()));
    TParameter<int>(cfg.tag + "_status", summary.status).Write();
    TParameter<int>(cfg.tag + "_hesse", summary.hesse).Write();
    TParameter<double>(cfg.tag + "_chi2", summary.chi2).Write();
    TParameter<int>(cfg.tag + "_ndf", summary.ndf).Write();

    if (expTerms == 1)
    {
      std::cout << Form("[INFO] %s exp1: MINIMIZE=%d HESSE=%d chi2/ndf=%.3f/%d A=%.5f B=%.5f C=%.5f",
                        cfg.tag.Data(),
                        summary.status,
                        summary.hesse,
                        summary.chi2,
                        summary.ndf,
                        func->GetParameter(0),
                        func->GetParameter(1),
                        func->GetParameter(2))
                << std::endl;
    }
    else
    {
      std::cout << Form("[INFO] %s exp2: MINIMIZE=%d HESSE=%d chi2/ndf=%.3f/%d A=%.5f B=%.5f D=%.5f E=%.5f C=%.5f",
                        cfg.tag.Data(),
                        summary.status,
                        summary.hesse,
                        summary.chi2,
                        summary.ndf,
                        func->GetParameter(0),
                        func->GetParameter(1),
                        func->GetParameter(2),
                        func->GetParameter(3),
                        func->GetParameter(4))
                << std::endl;
    }
  }

  fout->Close();

  std::cout << "[INFO] saved fit outputs to: " << rootsDir + Form("/fit_ratio_%s.root", expTag.Data()) << std::endl;
  std::cout << "[INFO] saved fit figures under: " << figsDir << std::endl;
}
