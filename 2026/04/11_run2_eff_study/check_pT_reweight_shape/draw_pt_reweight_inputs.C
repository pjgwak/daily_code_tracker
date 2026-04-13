#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace
{
const char *kRootlogonPath = "/data/users/pjgwak/input_files/rootlogon.C";

TString resolveBaseDir()
{
  return gSystem->DirName(__FILE__);
}

struct InputSpec
{
  TString key;
  TString label;
  TString headerLabel;
  TString collisionLabel;
  TString path;
  Color_t color;
  Style_t markerStyle;
};

struct PlotContent
{
  InputSpec spec;
  std::unique_ptr<TH1D> hist;
  std::unique_ptr<TF1> func;
  std::unique_ptr<TGraphErrors> graph;
};

TH1D *loadHist(TFile *file, const char *name, const char *cloneName)
{
  TH1D *hist = dynamic_cast<TH1D *>(file->Get(name));
  if (!hist)
    return nullptr;

  TH1D *cloned = static_cast<TH1D *>(hist->Clone(cloneName));
  cloned->SetDirectory(nullptr);
  cloned->SetStats(0);
  if (cloned->GetListOfFunctions())
    cloned->GetListOfFunctions()->Clear();
  return cloned;
}

TF1 *loadFunc(TFile *file, const char *name, const char *cloneName)
{
  TF1 *func = dynamic_cast<TF1 *>(file->Get(name));
  if (!func)
    return nullptr;

  return static_cast<TF1 *>(func->Clone(cloneName));
}

TGraphErrors *makeGraphFromHist(const TH1D *hist, const char *name)
{
  auto *graph = new TGraphErrors();
  graph->SetName(name);
  int ip = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    const double y = hist->GetBinContent(i);
    const double yErr = hist->GetBinError(i);
    if (y == 0.0 && yErr == 0.0)
      continue;
    graph->SetPoint(ip, hist->GetBinCenter(i), y);
    graph->SetPointError(ip, 0.0, yErr);
    ++ip;
  }
  return graph;
}

bool loadContent(const InputSpec &spec, PlotContent &content)
{
  std::unique_ptr<TFile> fin(TFile::Open(spec.path, "READ"));
  if (!fin || fin->IsZombie())
  {
    std::cerr << "[ERROR] failed to open " << spec.path << "\n";
    return false;
  }

  std::unique_ptr<TH1D> hist(loadHist(fin.get(), "WeightFactor", spec.key + "_hist"));
  std::unique_ptr<TF1> func(loadFunc(fin.get(), "dataMC_Ratio1", spec.key + "_func"));
  if (!hist || !func)
  {
    std::cerr << "[ERROR] missing WeightFactor or dataMC_Ratio1 in " << spec.path << "\n";
    return false;
  }

  content.spec = spec;
  content.hist = std::move(hist);
  content.func = std::move(func);
  content.graph.reset(makeGraphFromHist(content.hist.get(), spec.key + "_graph"));
  return true;
}

void styleHist(TGraphErrors *graph, Color_t color, Style_t markerStyle)
{
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(1.45);
  graph->SetLineWidth(0);
}

void styleFunc(TF1 *func, Color_t color)
{
  func->SetLineColor(color);
  func->SetLineWidth(3);
}

void styleFrame(TH1D *hist, const char *xTitle, const char *yTitle)
{
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetTitleOffset(1.20);
}

void makeCanvasPad(TCanvas &canv, TPad *&pad)
{
  canv.cd();
  pad = new TPad(Form("%s_pad", canv.GetName()), "", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.08);
  pad->SetBottomMargin(0.13);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->SetTicks(1, 1);
  pad->Draw();
  pad->cd();
}

void drawCmsText(const TString &collisionLabel)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(31);
  latex.SetTextFont(42);
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.96, 0.935, Form("%s 5.02 TeV", collisionLabel.Data()));
  latex.SetTextAlign(11);
  latex.SetTextFont(72);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.19, 0.935, "CMS");
  latex.SetTextFont(42);
  latex.SetTextSize(0.032);
  latex.DrawLatex(0.285, 0.935, "Internal");
}

void drawLabel(const TString &label)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.17, 0.86, label);
}

double histMinWithError(const TH1D *hist)
{
  double minVal = std::numeric_limits<double>::max();
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
    minVal = std::min(minVal, hist->GetBinContent(i) - hist->GetBinError(i));
  return minVal;
}

double histMaxWithError(const TH1D *hist)
{
  double maxVal = -std::numeric_limits<double>::max();
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
    maxVal = std::max(maxVal, hist->GetBinContent(i) + hist->GetBinError(i));
  return maxVal;
}

double funcMinInRange(TF1 *func, double xmin, double xmax)
{
  return func->GetMinimum(xmin, xmax);
}

double funcMaxInRange(TF1 *func, double xmin, double xmax)
{
  return func->GetMaximum(xmin, xmax);
}

void setRangeFromContent(TH1D *frame, const TH1D *hist, TF1 *func, double yFloor = 0.0)
{
  const double xmin = hist->GetXaxis()->GetXmin();
  const double xmax = hist->GetXaxis()->GetXmax();
  const double rawMin = std::min(histMinWithError(hist), funcMinInRange(func, xmin, xmax));
  const double rawMax = std::max(histMaxWithError(hist), funcMaxInRange(func, xmin, xmax));
  const double span = std::max(0.15, rawMax - rawMin);
  frame->SetMinimum(std::min(yFloor, rawMin - 0.08 * span));
  frame->SetMaximum(rawMax + 0.22 * span);
}

void drawReferenceLine(double xmin, double xmax)
{
  TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
  line->SetLineStyle(7);
  line->SetLineColor(kGray + 2);
  line->Draw();
}

void drawSingle(const PlotContent &content, const TString &outDir)
{
  TCanvas canv(Form("c_%s", content.spec.key.Data()), "", 800, 800);
  TPad *pad = nullptr;
  makeCanvasPad(canv, pad);

  TH1D *frame = static_cast<TH1D *>(content.hist->Clone(content.spec.key + "_frame"));
  frame->Reset("ICES");
  if (frame->GetListOfFunctions())
    frame->GetListOfFunctions()->Clear();
  styleFrame(frame, "p_{T} (GeV/c)", "Weight");
  setRangeFromContent(frame, content.hist.get(), content.func.get());
  frame->Draw();
  drawReferenceLine(frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
  content.func->Draw("SAME");
  content.graph->Draw("PZ SAME");

  drawCmsText(content.spec.collisionLabel);
  drawLabel(content.spec.headerLabel + ", |y| < 2.4");

  TLegend leg(0.50, 0.72, 0.90, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.028);
  leg.AddEntry(content.graph.get(), content.spec.label, "p");
  leg.Draw();

  const TString outBase = outDir + "/" + content.spec.key;
  canv.SaveAs(outBase + ".pdf");
  canv.SaveAs(outBase + ".png");
}

void drawSummary(const std::vector<PlotContent> &contents, const TString &outDir)
{
  if (contents.empty())
    return;

  TCanvas canv("c_ptw_summary", "", 800, 800);
  TPad *pad = nullptr;
  makeCanvasPad(canv, pad);

  TH1D *frame = static_cast<TH1D *>(contents.front().hist->Clone("summary_frame"));
  frame->Reset("ICES");
  if (frame->GetListOfFunctions())
    frame->GetListOfFunctions()->Clear();
  styleFrame(frame, "p_{T} (GeV/c)", "Weight");

  double yMin = std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  double xMin = contents.front().hist->GetXaxis()->GetXmin();
  double xMax = contents.front().hist->GetXaxis()->GetXmax();
  for (const auto &content : contents)
  {
    xMin = std::min(xMin, content.hist->GetXaxis()->GetXmin());
    xMax = std::max(xMax, content.hist->GetXaxis()->GetXmax());
    yMin = std::min(yMin, histMinWithError(content.hist.get()));
    yMin = std::min(yMin, funcMinInRange(content.func.get(), content.hist->GetXaxis()->GetXmin(), content.hist->GetXaxis()->GetXmax()));
    yMax = std::max(yMax, histMaxWithError(content.hist.get()));
    yMax = std::max(yMax, funcMaxInRange(content.func.get(), content.hist->GetXaxis()->GetXmin(), content.hist->GetXaxis()->GetXmax()));
  }
  const double span = std::max(0.15, yMax - yMin);
  frame->GetXaxis()->SetLimits(xMin, xMax);
  frame->SetMinimum(std::min(0.0, yMin - 0.10 * span));
  frame->SetMaximum(yMax + 0.24 * span);
  frame->Draw();
  drawReferenceLine(xMin, xMax);

  TLegend leg(0.34, 0.70, 0.90, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.028);
  leg.SetNColumns(2);

  for (const auto &content : contents)
  {
    content.func->Draw("SAME");
    content.graph->Draw("PZ SAME");
    leg.AddEntry(content.graph.get(), content.spec.label, "p");
  }

  drawCmsText("PbPb, pp");
  drawLabel("PbPb, pp, |y| < 2.4");
  leg.Draw();

  canv.SaveAs(outDir + "/pt_reweight_inputs_summary.pdf");
  canv.SaveAs(outDir + "/pt_reweight_inputs_summary.png");
}
} // namespace

void draw_pt_reweight_inputs()
{
  gROOT->Macro(kRootlogonPath);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);

  const TString baseDir = resolveBaseDir();
  const TString outDir = baseDir + "/outputs";
  gSystem->mkdir(outDir, true);

  const std::vector<InputSpec> specs = {
      {"pbpb_pr", "PbPb prompt", "PbPb", "PbPb", "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_2p4_260310.root", kAzure + 2, 20},
      {"pbpb_np", "PbPb nonprompt", "PbPb", "PbPb", "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_2p4_260310.root", kRed + 1, 21},
      {"pp_pr", "pp prompt", "pp", "pp", "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_2p4_260310.root", kGreen + 2, 24},
      {"pp_np", "pp nonprompt", "pp", "pp", "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_2p4_260310_2exp.root", kOrange + 7, 25},
  };

  std::vector<PlotContent> contents;
  contents.reserve(specs.size());
  for (const auto &spec : specs)
  {
    PlotContent content;
    if (!loadContent(spec, content))
      continue;
    styleHist(content.graph.get(), spec.color, spec.markerStyle);
    styleFunc(content.func.get(), spec.color);
    drawSingle(content, outDir);
    contents.push_back(std::move(content));
  }

  drawSummary(contents, outDir);
  std::cout << "[INFO] wrote plots under " << outDir << "\n";
}
