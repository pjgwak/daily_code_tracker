/*
  Plot the number of model components used in each kinematic bin.

  Run from any directory:
    root -l -b -q plot_model_counts.C

  Optional arguments:
    root -l -b -q 'plot_model_counts.C("../../roots",".")'

  Rapidity regions are distinguished by color and model types by marker
  style.  One PDF is written for every pre-final-2D fit stage.
*/

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TList.h"
#include "TMarker.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TLatex.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace
{
struct ModelType
{
  TString key;
  TString label;
};

struct Stage
{
  TString directory;
  TString filePrefix;
  TString title;
  std::vector<ModelType> models;
};

struct RapidityBin
{
  TString directory;
  double low = 0.;
  double high = 0.;
};

struct BinValue
{
  double ptLow = 0.;
  double ptHigh = 0.;
  int count = 0;
};

bool parseRapidityDirectory(const TString &name, double &low, double &high)
{
  return std::sscanf(name.Data(), "y%lf_%lf", &low, &high) == 2;
}

bool parsePtRange(const TString &name, double &low, double &high)
{
  const Ssiz_t pos = name.Index("_pt");
  if (pos == kNPOS)
    return false;
  return std::sscanf(name.Data() + pos + 1, "pt%lf_%lf.root", &low, &high) == 2;
}

std::vector<TString> listDirectories(const TString &path)
{
  std::vector<TString> result;
  TSystemDirectory directory("directory", path);
  std::unique_ptr<TList> entries(directory.GetListOfFiles());
  if (!entries)
    return result;

  TIter next(entries.get());
  while (auto *entry = dynamic_cast<TSystemFile *>(next()))
  {
    const TString name = entry->GetName();
    if (entry->IsDirectory() && name != "." && name != "..")
      result.push_back(name);
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<TString> listRootFiles(const TString &path, const TString &prefix)
{
  std::vector<TString> result;
  TSystemDirectory directory("directory", path);
  std::unique_ptr<TList> entries(directory.GetListOfFiles());
  if (!entries)
    return result;

  TIter next(entries.get());
  while (auto *entry = dynamic_cast<TSystemFile *>(next()))
  {
    const TString name = entry->GetName();
    if (!entry->IsDirectory() && name.BeginsWith(prefix) && name.EndsWith(".root"))
      result.push_back(name);
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::vector<RapidityBin> findRapidityBins(const TString &rootsDir)
{
  std::vector<RapidityBin> bins;
  for (const TString &name : listDirectories(rootsDir))
  {
    RapidityBin bin;
    bin.directory = name;
    if (parseRapidityDirectory(name, bin.low, bin.high))
      bins.push_back(bin);
  }
  std::sort(bins.begin(), bins.end(), [](const RapidityBin &a, const RapidityBin &b) {
    return a.low == b.low ? a.high < b.high : a.low < b.low;
  });
  return bins;
}

TString rapidityLabel(const RapidityBin &bin)
{
  if (std::abs(bin.low) < 1.e-9)
    return TString::Format("|y| < %.1f", bin.high);
  return TString::Format("%.1f < |y| < %.1f", bin.low, bin.high);
}

std::vector<BinValue> readModelCounts(const TString &directory,
                                      const TString &filePrefix,
                                      const TString &parameterKey,
                                      int &missingObjects)
{
  std::vector<BinValue> values;
  for (const TString &fileName : listRootFiles(directory, filePrefix))
  {
    BinValue point;
    if (!parsePtRange(fileName, point.ptLow, point.ptHigh))
    {
      std::cerr << "[skip] Cannot parse pT range from " << fileName << std::endl;
      continue;
    }

    const TString path = directory + "/" + fileName;
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    if (!file || file->IsZombie())
    {
      std::cerr << "[skip] Cannot open " << path << std::endl;
      continue;
    }

    auto *parameter = dynamic_cast<TParameter<int> *>(file->Get(parameterKey));
    if (!parameter)
    {
      ++missingObjects;
      std::cerr << "[skip] " << parameterKey << " is missing in " << path << std::endl;
      continue;
    }
    point.count = parameter->GetVal();
    values.push_back(point);
  }

  std::sort(values.begin(), values.end(), [](const BinValue &a, const BinValue &b) {
    return a.ptLow == b.ptLow ? a.ptHigh < b.ptHigh : a.ptLow < b.ptLow;
  });
  return values;
}

void drawStage(const Stage &stage,
               const std::vector<RapidityBin> &rapidityBins,
               const TString &rootsDir,
               const TString &outputDir)
{
  const int colors[] = {kBlue + 1, kRed + 1, kAzure + 2, kGreen + 2,
                        kOrange + 7, kCyan + 2};
  const int markers[] = {24, 25, 26, 32, 27, 28, 30};
  const int lineStyles[] = {1, 2, 3, 4, 5, 6, 7, 8};

  struct Series
  {
    std::unique_ptr<TGraphErrors> graph;
    size_t rapidityIndex = 0;
    size_t modelIndex = 0;
  };

  std::vector<Series> series;
  double xMin = std::numeric_limits<double>::infinity();
  double xMax = -std::numeric_limits<double>::infinity();
  int largestCount = 0;
  int missingObjects = 0;

  for (size_t iy = 0; iy < rapidityBins.size(); ++iy)
  {
    const TString stageDir = rootsDir + "/" + rapidityBins[iy].directory + "/" + stage.directory;
    if (gSystem->AccessPathName(stageDir))
      continue;

    for (size_t im = 0; im < stage.models.size(); ++im)
    {
      const auto points = readModelCounts(stageDir, stage.filePrefix,
                                          stage.models[im].key, missingObjects);
      if (points.empty())
        continue;

      Series item;
      item.graph = std::make_unique<TGraphErrors>(points.size());
      item.rapidityIndex = iy;
      item.modelIndex = im;
      item.graph->SetName(TString::Format("g_%s_y%zu_m%zu", stage.directory.Data(), iy, im));
      item.graph->SetMarkerColor(colors[iy % (sizeof(colors) / sizeof(colors[0]))]);
      item.graph->SetLineColor(colors[iy % (sizeof(colors) / sizeof(colors[0]))]);
      item.graph->SetMarkerStyle(markers[im % (sizeof(markers) / sizeof(markers[0]))]);
      item.graph->SetLineStyle(lineStyles[im % (sizeof(lineStyles) / sizeof(lineStyles[0]))]);
      item.graph->SetMarkerSize(1.35);
      item.graph->SetLineWidth(2);

      for (size_t ip = 0; ip < points.size(); ++ip)
      {
        const double rapidityOffset =
            (static_cast<double>(iy) - 0.5 * (rapidityBins.size() - 1)) *
            0.08 * (points[ip].ptHigh - points[ip].ptLow);
        const double ptCenter =
            0.5 * (points[ip].ptLow + points[ip].ptHigh) + rapidityOffset;
        item.graph->SetPoint(ip, ptCenter, points[ip].count);
        item.graph->SetPointError(ip, 0.5 * (points[ip].ptHigh - points[ip].ptLow), 0.0);
        xMin = std::min(xMin, points[ip].ptLow);
        xMax = std::max(xMax, points[ip].ptHigh);
        largestCount = std::max(largestCount, points[ip].count);
      }
      series.push_back(std::move(item));
    }
  }

  if (series.empty())
  {
    std::cerr << "[skip] No model-count values found for stage "
              << stage.directory << std::endl;
    return;
  }

  const TString cmsStyle = "/data/users/pjgwak/input_files/rootlogon.C";
  if (!gSystem->AccessPathName(cmsStyle))
    gROOT->Macro(cmsStyle);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);
  gStyle->SetLineWidth(2);

  const TString canvasName = "c_model_count_" + stage.directory;
  auto canvas = std::make_unique<TCanvas>(canvasName, canvasName, 800, 800);
  canvas->SetTitle("");
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.08);
  canvas->SetBottomMargin(0.13);
  canvas->SetTicks(1, 1);

  const int yMaximum = std::max(3, largestCount + 2);
  const TString frameName = "frame_" + stage.directory;
  TH1F frame(frameName, "", 100, 0.0, xMax);
  frame.SetDirectory(nullptr);
  frame.SetMinimum(0.0);
  frame.SetMaximum(yMaximum);
  frame.SetLineWidth(2);
  frame.GetXaxis()->SetTitle("p_{T} (GeV/c)");
  frame.GetYaxis()->SetTitle("Number of model components");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.GetXaxis()->SetTitleFont(42);
  frame.GetYaxis()->SetTitleFont(42);
  frame.GetXaxis()->SetLabelFont(42);
  frame.GetYaxis()->SetLabelFont(42);
  frame.GetXaxis()->SetTitleSize(0.050);
  frame.GetYaxis()->SetTitleSize(0.050);
  frame.GetXaxis()->SetLabelSize(0.042);
  frame.GetYaxis()->SetLabelSize(0.042);
  frame.GetXaxis()->SetTitleOffset(1.05);
  frame.GetYaxis()->SetTitleOffset(1.30);
  frame.GetXaxis()->SetNdivisions(510);
  frame.GetYaxis()->SetNdivisions(yMaximum, false);
  frame.Draw("AXIS");

  for (auto &item : series)
    item.graph->Draw("P SAME");

  std::vector<std::unique_ptr<TGraph>> legendSamples;
  legendSamples.reserve(rapidityBins.size() + stage.models.size());
  std::set<size_t> usedRapidity;
  std::set<size_t> usedModels;
  for (const auto &item : series)
  {
    usedRapidity.insert(item.rapidityIndex);
    usedModels.insert(item.modelIndex);
  }

  const double legendTop = 0.87;
  const double legendRowHeight = 0.047;
  const double rapidityBottom =
      legendTop - legendRowHeight * (usedRapidity.size() + 1);
  TLegend rapidityLegend(0.18, rapidityBottom, 0.44, legendTop);
  rapidityLegend.SetBorderSize(0);
  rapidityLegend.SetFillStyle(0);
  rapidityLegend.SetTextFont(42);
  rapidityLegend.SetTextSize(0.026);
  rapidityLegend.SetEntrySeparation(0.01);
  rapidityLegend.SetHeader("Rapidity", "L");
  for (size_t iy : usedRapidity)
  {
    auto sample = std::make_unique<TGraph>();
    sample->SetMarkerColor(colors[iy % (sizeof(colors) / sizeof(colors[0]))]);
    sample->SetLineColor(colors[iy % (sizeof(colors) / sizeof(colors[0]))]);
    sample->SetMarkerStyle(24);
    sample->SetMarkerSize(1.35);
    sample->SetLineWidth(2);
    rapidityLegend.AddEntry(sample.get(), rapidityLabel(rapidityBins[iy]), "p");
    legendSamples.push_back(std::move(sample));
  }
  rapidityLegend.Draw();

  const double modelBottom =
      legendTop - legendRowHeight * (usedModels.size() + 1);
  TLegend modelLegend(0.48, modelBottom, 0.92, legendTop);
  modelLegend.SetBorderSize(0);
  modelLegend.SetFillStyle(0);
  modelLegend.SetTextFont(42);
  modelLegend.SetTextSize(0.026);
  modelLegend.SetEntrySeparation(0.01);
  modelLegend.SetHeader(stage.title, "L");
  for (size_t im : usedModels)
  {
    auto sample = std::make_unique<TGraph>();
    sample->SetMarkerColor(kBlack);
    sample->SetLineColor(kBlack);
    sample->SetMarkerStyle(markers[im % (sizeof(markers) / sizeof(markers[0]))]);
    sample->SetMarkerSize(1.35);
    sample->SetLineWidth(2);
    modelLegend.AddEntry(sample.get(), stage.models[im].label, "p");
    legendSamples.push_back(std::move(sample));
  }
  modelLegend.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.030);
  label.SetTextAlign(31);
  label.DrawLatex(0.96, 0.935, "OO 5.36 TeV (9 nb^{-1})");
  label.SetTextAlign(11);
  label.SetTextFont(72);
  label.SetTextSize(0.040);
  label.DrawLatex(0.15, 0.935, "CMS");
  label.SetTextFont(52);
  label.SetTextSize(0.032);
  label.DrawLatex(0.265, 0.935, "Preliminary");

  canvas->RedrawAxis();
  const TString baseName = outputDir + "/model_count_" + stage.directory;
  canvas->SaveAs(baseName + ".pdf");
  std::cout << "[done] " << baseName << ".pdf";
  if (missingObjects > 0)
    std::cout << " (missing parameters: " << missingObjects << ")";
  std::cout << std::endl;
}
} // namespace

void plot_model_counts(TString rootsDir = "", TString outputDir = "")
{
  gROOT->SetBatch(kTRUE);

  const TString macroDir = gSystem->DirName(__FILE__);
  if (rootsDir.IsNull())
    rootsDir = macroDir + "/../../roots";
  if (outputDir.IsNull())
    outputDir = macroDir;

  if (gSystem->AccessPathName(rootsDir))
  {
    std::cerr << "ERROR: ROOT result directory does not exist: " << rootsDir << std::endl;
    return;
  }
  gSystem->mkdir(outputDir, true);

  const std::vector<RapidityBin> rapidityBins = findRapidityBins(rootsDir);
  if (rapidityBins.empty())
  {
    std::cerr << "ERROR: no y<low>_<high> directories found below " << rootsDir << std::endl;
    return;
  }

  const std::vector<Stage> stages = {
      {"mc_mass", "mc_mass_model_", "Mass signal",
       {{"nSignalGaussComponents", "Signal Gaussian"},
        {"nSignalCBComponents", "Signal Crystal Ball"}}},
      {"mass", "mass_model_", "Mass background",
       {{"nBkgExpComponents", "Background exponential"},
        {"nBkgChebyOrder", "Background Chebyshev order"}}},
      {"ctau_pr", "ctau_resolution_", "c#tau resolution",
       {{"nResolutionComponents", "Resolution Gaussian"}}},
      {"ctau_np", "ctau_np_model_", "Nonprompt c#tau",
       {{"nSignalSSComponents", "Nonprompt single-sided"}}},
      {"ctau_bkg", "ctau_bkg_fitresult_", "Background c#tau",
       {{"nSignalSSComponents", "Single-sided"},
        {"nSignalFlipComponents", "Flipped single-sided"},
        {"nSignalDSComponents", "Double-sided"}}}};

  // fit2d is intentionally excluded: the request only needs the pre-final-2D
  // model choices, one plot per fit stage.
  for (const Stage &stage : stages)
    drawStage(stage, rapidityBins, rootsDir, outputDir);
}
