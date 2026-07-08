#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <vector>

namespace DrawCentLocal
{
const int kColorHFp = kRed + 1;
const int kColorHFm = kBlue + 1;
const int kColorHF = kGreen + 2;
const int kColorIntegrated = kOrange + 7;
const int kColorFine = kAzure + 2;

struct InputConfig
{
  TString sample;
  TString label;
  TString path;
};

TString ProjectDir()
{
  // Start from this macro path.
  TString path = __FILE__;
  // Make relative macro paths usable from the current working directory.
  if (!gSystem->IsAbsoluteFileName(path))
    path = TString(gSystem->WorkingDirectory()) + "/" + path;

  // drawing/centrality -> pilot_calcTrack_3subEvent.
  TString project = TString(gSystem->DirName(path)) + "/../..";
  // Canonicalize symlinks and relative path components.
  TString out = gSystem->GetFromPipe(Form("cd %s && pwd -P", project.Data()));
  out = out.Strip(TString::kBoth);
  return out.IsNull() ? project : out;
}

TString OutputDir(const TString &sample, const TString &category)
{
  TString out = ProjectDir() + "/drawing/" + category + "/figs/" + sample;
  gSystem->mkdir(out, true);
  return out;
}

void ApplyStyle()
{
  gROOT->SetBatch(kTRUE);
  const char *rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  if (!gSystem->AccessPathName(rootlogon))
    gROOT->ProcessLine(Form(".X %s", rootlogon));

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetEndErrorSize(4);
}

void SetupPad()
{
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);
  gPad->SetTicks(1, 1);
}

void FormatAxis(TAxis *axis, bool largePanel = false)
{
  if (!axis)
    return;
  axis->CenterTitle(true);
  axis->SetTitleSize(largePanel ? 0.050 : 0.046);
  axis->SetLabelSize(largePanel ? 0.042 : 0.040);
}

TLegend *MakeLegend(double x1, double y1, double x2, double y2, double textSize = 0.030)
{
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(textSize);
  return leg;
}

void DrawCmsLabels(const TString &sample, const char *detail = "",
                   const char *collision = "PbPb 5.36 TeV")
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.95, 0.93, collision);

  TLatex cms;
  cms.SetNDC();
  cms.SetTextSize(0.040);
  cms.SetTextFont(62);
  cms.SetTextAlign(11);
  cms.DrawLatex(0.205, 0.93, "CMS");

  TLatex extra;
  extra.SetNDC();
  extra.SetTextSize(0.040);
  extra.SetTextFont(52);
  extra.SetTextAlign(11);
  extra.DrawLatex(0.289, 0.93, "Internal");

  TLatex sub;
  sub.SetNDC();
  sub.SetTextFont(42);
  sub.SetTextSize(0.030);
  sub.SetTextAlign(11);
  const TString label = detail && detail[0] ? Form("%s, %s", sample.Data(), detail) : sample;
  sub.DrawLatex(0.18, 0.865, label);
}

void DrawPanelLabel(const char *text, double x = 0.18, double y = 0.86,
                    double size = 0.034)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(size);
  latex.DrawLatex(x, y, text);
}

void SaveCanvas(TCanvas *c, const TString &sample, const TString &category, const TString &name)
{
  if (!c)
    return;
  const TString outDir = OutputDir(sample, category);
  c->SaveAs(Form("%s/%s.pdf", outDir.Data(), name.Data()));
}

TString NormalizeSampleTag(const TString &tag)
{
  TString sample(tag);
  sample.ToLower();
  if (sample == "prompt" || sample == "pr")
    return "pr";
  if (sample == "nonprompt" || sample == "np")
    return "np";
  if (sample == "data" || sample == "data23" || sample == "pbpb2023")
    return "data";
  return sample;
}

TString ExistingFile(const TString &path)
{
  if (!path.IsNull() && !gSystem->AccessPathName(path))
    return path;
  return "";
}

TString FirstExisting(std::initializer_list<TString> paths)
{
  for (const auto &path : paths)
  {
    TString existing = ExistingFile(path);
    if (!existing.IsNull())
      return existing;
  }
  return "";
}

TString InputFileForSample(const TString &sample)
{
  const TString project = ProjectDir();

  if (sample == "pr")
    return FirstExisting({project + "/roots/pr/skim_PbPb2023_isMC1_PR_EP_3Sub_job.root"});

  if (sample == "np")
    return FirstExisting({project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test40k.root",
                          // project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test10k.root",
                          // project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test1k.root"
                        });

  if (sample == "data")
    return FirstExisting({project + "/roots/data/Oniatree_MC0_PbPb_2023_EP_3Sub_job.root",
                          // project + "/roots/data/Oniatree_MC0_PbPb_2023_EP_angles_3Sub_job.root"
                        });

  return "";
}

InputConfig GetInputConfig(const char *sampleTag)
{
  InputConfig cfg;
  // Normalize aliases like prompt/nonprompt/data23.
  cfg.sample = NormalizeSampleTag(sampleTag && sampleTag[0] ? sampleTag : "pr");
  // Pick the plot label for this sample.
  cfg.label = cfg.sample == "pr" ? "MC PR" : (cfg.sample == "np" ? "MC NP" : (cfg.sample == "data" ? "Data" : cfg.sample));
  // Pick the ROOT input file under roots/<sample>.
  cfg.path = InputFileForSample(cfg.sample);

  if (cfg.sample != "pr" && cfg.sample != "np" && cfg.sample != "data")
  {
    std::cout << "[ERROR] unknown sample tag: " << cfg.sample << " (use pr, np, or data)" << std::endl;
    return cfg;
  }

  if (cfg.path.IsNull())
    std::cout << "[ERROR] no ROOT input configured/found for sample tag '" << cfg.sample
              << "' under " << ProjectDir() << "/roots/" << cfg.sample << std::endl;
  else
    std::cout << "[INFO] sample = " << cfg.sample << ", input = " << cfg.path << std::endl;
  return cfg;
}
}


namespace
{
const int kNFineCentBins = 18;
const int kNIntegratedCentBins = 6;
const int kNAllCentBins = kNFineCentBins + kNIntegratedCentBins;

const double kCentLow[kNAllCentBins] = {
    0, 5, 10, 15, 20, 25, 30, 35, 40,
    45, 50, 55, 60, 65, 70, 75, 80, 85,
    0, 0, 10, 20, 10, 40};
const double kCentHigh[kNAllCentBins] = {
    5, 10, 15, 20, 25, 30, 35, 40, 45,
    50, 55, 60, 65, 70, 75, 80, 85, 90,
    10, 20, 20, 40, 40, 90};
enum Metric
{
  kResolution = 0,
  kTrkHFp = 1,
  kTrkHFm = 2,
  kHFpHFm = 3,
  kEntries = 4
};

struct CentPoint
{
  double low;
  double high;
  double trkHFp;
  double trkHFpErr;
  double trkHFm;
  double trkHFmErr;
  double hfpHFm;
  double hfpHFmErr;
  double resolution;
  double resolutionErr;
  double entries;
  bool valid;
};

double SafeResolution(double numerator, double denominator)
{
  if (denominator <= 0.0 || numerator / denominator < 0.0)
    return -1.0;
  return std::sqrt(numerator / denominator);
}

double ResolutionError(double a, double ea, double b, double eb, double c, double ec)
{
  const double r = SafeResolution(a * b, c);
  if (r <= 0.0 || a == 0.0 || b == 0.0 || c == 0.0)
    return 0.0;
  const double rel2 = (ea / a) * (ea / a) + (eb / b) * (eb / b) + (ec / c) * (ec / c);
  return 0.5 * r * std::sqrt(std::max(0.0, rel2));
}

CentPoint ReadPoint(TFile *f, int idx)
{
  CentPoint p{};
  p.low = kCentLow[idx];
  p.high = kCentHigh[idx];

  TH1 *hTrkHFp = f ? (TH1 *)f->Get(Form("hCos2TrkHFp_cent%d", idx)) : nullptr;
  TH1 *hTrkHFm = f ? (TH1 *)f->Get(Form("hCos2TrkHFm_cent%d", idx)) : nullptr;
  TH1 *hHFpHFm = f ? (TH1 *)f->Get(Form("hCos2HFpHFm_cent%d", idx)) : nullptr;
  if (!hTrkHFp || !hTrkHFm || !hHFpHFm)
    return p;

  p.trkHFp = hTrkHFp->GetMean(1);
  p.trkHFpErr = hTrkHFp->GetMeanError(1);
  p.trkHFm = hTrkHFm->GetMean(1);
  p.trkHFmErr = hTrkHFm->GetMeanError(1);
  p.hfpHFm = hHFpHFm->GetMean(1);
  p.hfpHFmErr = hHFpHFm->GetMeanError(1);
  p.resolution = SafeResolution(p.trkHFp * p.trkHFm, p.hfpHFm);
  p.resolutionErr = ResolutionError(p.trkHFp, p.trkHFpErr, p.trkHFm, p.trkHFmErr, p.hfpHFm, p.hfpHFmErr);
  p.entries = hTrkHFp->GetEntries();
  p.valid = p.entries > 0;
  return p;
}

double Value(const CentPoint &p, Metric metric)
{
  if (metric == kResolution)
    return p.resolution;
  if (metric == kTrkHFp)
    return p.trkHFp;
  if (metric == kTrkHFm)
    return p.trkHFm;
  if (metric == kHFpHFm)
    return p.hfpHFm;
  return p.entries;
}

double Error(const CentPoint &p, Metric metric)
{
  if (metric == kResolution)
    return p.resolutionErr;
  if (metric == kTrkHFp)
    return p.trkHFpErr;
  if (metric == kTrkHFm)
    return p.trkHFmErr;
  if (metric == kHFpHFm)
    return p.hfpHFmErr;
  return 0.0;
}

TGraphErrors *MakeGraph(const std::vector<CentPoint> &points, Metric metric, int color, int marker, const char *name)
{
  TGraphErrors *g = new TGraphErrors();
  g->SetName(name);
  int ip = 0;
  for (const auto &p : points)
  {
    if (!p.valid)
      continue;
    g->SetPoint(ip, 0.5 * (p.low + p.high), Value(p, metric));
    g->SetPointError(ip, 0.5 * (p.high - p.low), Error(p, metric));
    ++ip;
  }
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.20);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  return g;
}

void DrawIntegratedBoxes(const std::vector<CentPoint> &points, Metric metric, int color)
{
  for (const auto &p : points)
  {
    if (!p.valid)
      continue;
    const double y = Value(p, metric);
    const double ey = Error(p, metric);
    TBox *box = new TBox(p.low, y - ey, p.high, y + ey);
    box->SetFillColorAlpha(color, 0.18);
    box->SetFillStyle(1001);
    box->SetLineColor(color);
    box->SetLineWidth(2);
    box->Draw("same");
  }
}

double FindMaxY(const std::vector<CentPoint> &a, const std::vector<CentPoint> &b, Metric metric)
{
  double ymax = 0.0;
  for (const auto &p : a)
    if (p.valid)
      ymax = std::max(ymax, Value(p, metric) + Error(p, metric));
  for (const auto &p : b)
    if (p.valid)
      ymax = std::max(ymax, Value(p, metric) + Error(p, metric));
  return ymax > 0.0 ? ymax : 1.0;
}

void DrawMetricFrame(const char *yTitle, double ymax)
{
  TH1F *frame = new TH1F(Form("frame_%p", gPad), "", 90, 0, 90);
  frame->SetStats(0);
  frame->SetTitle(Form(";Centrality (%%);%s", yTitle));
  DrawCentLocal::FormatAxis(frame->GetXaxis(), true);
  DrawCentLocal::FormatAxis(frame->GetYaxis(), true);
  frame->GetYaxis()->SetTitleOffset(1.20);
  frame->SetMinimum(0.0);
  frame->SetMaximum(ymax);
  frame->Draw("axis");
}
}

void draw_cent(const char *sampleTag = "pr")
{
  using namespace DrawCentLocal; // local namespace
  TStopwatch timer;
  timer.Start();

  ApplyStyle(); // cosmetics
  InputConfig input = GetInputConfig(sampleTag);
  TString sample(input.sample); // tagging
  TString sampleLabel(input.label);

  // input file check
  if (input.path.IsNull())
    return;
  TFile *f = TFile::Open(input.path, "READ");
  if (!f || f->IsZombie())
  {
    std::cout << "[ERROR] cannot open input file: " << input.path << std::endl;
    return;
  }

  // decide number of points
  std::vector<CentPoint> finePoints;
  std::vector<CentPoint> integratedPoints;
  for (int i = 0; i < kNFineCentBins; ++i)
    finePoints.push_back(ReadPoint(f, i));
  for (int i = kNFineCentBins; i < kNAllCentBins; ++i)
    integratedPoints.push_back(ReadPoint(f, i));

  // draw EP resolution depending on centrality
  TCanvas *cRes = new TCanvas("cCentResolutionDetail", "centrality resolution", 850, 800);
  SetupPad();
  DrawMetricFrame("R_{trk}", std::max(1.05, 1.25 * FindMaxY(finePoints, integratedPoints, kResolution)));
  DrawIntegratedBoxes(integratedPoints, kResolution, kColorIntegrated);
  TGraphErrors *gFineRes = MakeGraph(finePoints, kResolution, kColorFine, 20, "gFineRes");
  TGraphErrors *gIntegratedRes = MakeGraph(integratedPoints, kResolution, kColorIntegrated, 24, "gIntegratedRes");
  gFineRes->Draw("P SAME");
  gIntegratedRes->Draw("P SAME");
  DrawCmsLabels(sampleLabel, "");
  TLegend *legRes = MakeLegend(0.7, 0.8, 0.88, 0.9);
  legRes->AddEntry(gFineRes, "5% bins", "lp");
  legRes->AddEntry(gIntegratedRes, "Integrated bins", "lp");
  legRes->Draw();
  SaveCanvas(cRes, sample, "centrality", "cent_resolution");

  // draw difference between detectors (these values were used for EP resolution calculation)
  TCanvas *cCos2 = new TCanvas("cCentCos2Detail", "centrality cos2", 1800, 600);
  cCos2->Divide(3, 1);
  const Metric metrics[3] = {kTrkHFp, kTrkHFm, kHFpHFm};
  const char *titles[3] = {"<cos 2(#Psi_{trk}-#Psi_{HF_{+}})>", "<cos 2(#Psi_{trk}-#Psi_{HF_{-}})>", "<cos 2(#Psi_{HF_{+}}-#Psi_{HF_{-}})>"};
  const char *names[3] = {"Trk - HF_{+}", "Trk - HF_{-}", "HF_{+} - HF_{-}"};
  const int colors[3] = {kColorHFp, kColorHFm, kColorHF};
  // draw plot in three pads
  for (int i = 0; i < 3; ++i)
  {
    cCos2->cd(i + 1);
    SetupPad();
    DrawMetricFrame(titles[i], std::max(1.0, 1.25 * FindMaxY(finePoints, integratedPoints, metrics[i])));
    DrawIntegratedBoxes(integratedPoints, metrics[i], colors[i]);
    TGraphErrors *gFine = MakeGraph(finePoints, metrics[i], colors[i], 20, Form("gFineCos2_%d", i));
    TGraphErrors *gIntegrated = MakeGraph(integratedPoints, metrics[i], kBlack, 24, Form("gIntegratedCos2_%d", i));
    gFine->Draw("P SAME");
    gIntegrated->Draw("P SAME");
    if (i == 0)
      DrawCmsLabels(sampleLabel, "");
    DrawPanelLabel(names[i], 0.18, i == 0 ? 0.82 : 0.86, 0.034);
  }
  SaveCanvas(cCos2, sample, "centrality", "cent_cos2");

  // check number of events in centrality bins
  TCanvas *cEntries = new TCanvas("cCentEntriesDetail", "centrality entries", 850, 800);
  SetupPad();
  DrawMetricFrame("Entries", 1.35 * FindMaxY(finePoints, integratedPoints, kEntries));
  TGraphErrors *gFineEntries = MakeGraph(finePoints, kEntries, kColorFine, 20, "gFineEntries");
  TGraphErrors *gIntegratedEntries = MakeGraph(integratedPoints, kEntries, kColorIntegrated, 24, "gIntegratedEntries");
  gFineEntries->Draw("P SAME");
  gIntegratedEntries->Draw("P SAME");
  DrawCmsLabels(sampleLabel, "");
  TLegend *legEntries = MakeLegend(0.7, 0.8, 0.88, 0.9);
  legEntries->AddEntry(gFineEntries, "5% bins", "lp");
  legEntries->AddEntry(gIntegratedEntries, "Integrated bins", "lp");
  legEntries->Draw();
  SaveCanvas(cEntries, sample, "centrality", "cent_entries");

  // save output
  f->Close();

  timer.Stop();
  timer.Print();
}
