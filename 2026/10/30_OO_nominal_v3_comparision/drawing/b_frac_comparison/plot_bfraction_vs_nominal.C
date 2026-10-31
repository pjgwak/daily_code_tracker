#include "../latest_fit_results.h"
#include "../../utils/nominal_bfraction.h"
#include <iomanip>
// KST - Read saved mid or forward ROOT fits and compare with their regional nominal references.
// Mid: root -n -l -b -q 'drawing/b_frac_comparison/plot_bfraction_vs_nominal.C("","",2.5,"mid")'
// Nominal values come from the shared ../../utils/nominal_bfraction.h (see readBins).
//
// From this folder:
//   root -l -b -q 'plot_bfraction_vs_nominal.C()'
// From the analysis directory:
//   root -l -b -q 'drawing/b_frac_comparison/plot_bfraction_vs_nominal.C()'
// Optional ROOT input directory and PDF output directory:
//   root -l -b -q 'plot_bfraction_vs_nominal.C("/path/to/fit2d", "plots")'
// Defaults resolve relative to this macro, independently of the working directory.
// Output: bfraction_production_vs_nominal.pdf beside this macro.

#include "TCanvas.h"
#include "TROOT.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TParameter.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bfraction_comparison
{
struct Bin
{
  double low, high, nominal, nominalError, fit, fitError;
};

TString macroDirectory()
{
  TString path = __FILE__;
  if (!gSystem->IsAbsoluteFileName(path))
    path = TString(gSystem->WorkingDirectory()) + "/" + path;
  return gSystem->DirName(path);
}

// Nominal reference is read from ../../utils/nominal_bfraction.h: nominal ver2 = AN-26-048 nominal ctau template fit (2026-10-28); update the shared header, not this macro.
// KST - Share plotting code while reading the correct regional nominal table and ROOT tags.
std::vector<Bin> readBins(const TString &rootDir, bool mid = false)
{
  std::vector<Bin> bins;
  if (mid) {
    for (const auto &reference : nominal_bfraction::midBins())
      bins.push_back({reference.low, reference.high, reference.value, reference.error, 0, 0});
  } else {
    for (const auto &reference : nominal_bfraction::bins())
      bins.push_back({reference.low, reference.high, reference.value, reference.error, 0, 0});
  }
  for (auto &bin : bins)
  {
    const TString path = Form("%s/fit2d_result_%s_pt%.2f_%.2f.root",
                              rootDir.Data(), mid ? "y0.00_1.60" : "y1.60_2.40", bin.low, bin.high);
    const auto latest=oo_drawing::readFit(path);
    bin.fit=latest.b;bin.fitError=latest.bError;
    std::cout << Form("pT %g-%g: ROOT %.8f +/- %.8f; nominal %.8f +/- %.8f",
                      bin.low, bin.high, bin.fit, bin.fitError, bin.nominal, bin.nominalError)
              << '\n';
  }
  return bins;
}

void styleGraph(TGraphErrors &graph, Color_t color, Style_t marker)
{
  graph.SetLineColor(color);
  graph.SetMarkerColor(color);
  graph.SetMarkerStyle(marker);
  graph.SetMarkerSize(1.15);
  graph.SetLineWidth(2);
}

void drawGrid(const TH1D &frame, double step)
{
  TLine line;
  line.SetLineColor(TColor::GetColor("#e8e8e8"));
  line.SetLineWidth(1);
  for (double y = step * std::ceil(frame.GetMinimum() / step); y < frame.GetMaximum(); y += step)
    line.DrawLine(frame.GetXaxis()->GetXmin(), y, frame.GetXaxis()->GetXmax(), y);
}

void styleFrame(TH1 &frame, const char *yTitle)
{
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetTitle("");
  frame.GetYaxis()->SetTitle(yTitle);
  for (TAxis *axis : {frame.GetXaxis(), frame.GetYaxis()})
  {
    axis->SetTitleFont(43);
    axis->SetLabelFont(43);
    axis->SetTitleSize(23);
    axis->SetLabelSize(20);
    axis->CenterTitle();
  }
  frame.GetXaxis()->SetTickLength(0.015);
  frame.GetYaxis()->SetTickLength(0.008);
  frame.GetYaxis()->SetTitleOffset(1.7);
  frame.GetYaxis()->SetNdivisions(505);
}
} // namespace bfraction_comparison

void plot_bfraction_vs_nominal(const char *rootDir = "", const char *outputDir = "", double targetPercent = nominal_bfraction::targetPercent, const char *region = "fwd")
{
  using namespace bfraction_comparison;
  // KST - Existing calls remain forward; pass fourth argument "mid" for promoted mid results.
  const bool mid = TString(region) == "mid";
  if (!mid && TString(region) != "fwd") throw std::runtime_error("Region must be mid or fwd");
  if (!std::isfinite(targetPercent) || targetPercent <= 0)
    throw std::runtime_error("Agreement target must be positive and finite");
  // Use the common ROOT style, as in drawing/model_count.
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  const TString directory = macroDirectory();
  const TString input = rootDir[0] ? TString(rootDir)
      : directory + (mid ? "/../../roots/y0.00_1.60/fit2d" : "/../../roots/y1.60_2.40/fit2d");
  const TString output = outputDir[0] ? TString(outputDir) : directory;
  const auto bins = readBins(input, mid);
  if (gSystem->AccessPathName(output) && gSystem->mkdir(output, true) != 0)
    throw std::runtime_error("Cannot create output directory: " + std::string(output.Data()));

  void *outputDirectory = gSystem->OpenDirectory(output);
  if (!outputDirectory)
    throw std::runtime_error("Output path is not a directory: " + std::string(output.Data()));
  gSystem->FreeDirectory(outputDirectory);
  if (gSystem->AccessPathName(output, kWritePermission))
    throw std::runtime_error("Output directory is not writable: " + std::string(output.Data()));

  const int count = static_cast<int>(bins.size());
  const Color_t blue = TColor::GetColor("#226fba");
  const Color_t orange = TColor::GetColor("#d95319");
  const Color_t green = TColor::GetColor("#d9eee2");
  TGraphErrors nominal(count), band(count), within, outside, relativeWithin, relativeOutside;
  styleGraph(nominal, kBlack, 21);
  styleGraph(within, blue, 20);
  styleGraph(outside, orange, 5);
  styleGraph(relativeWithin, blue, 20);
  styleGraph(relativeOutside, orange, 5);
  band.SetFillColor(green);
  band.SetLineColor(green);
  band.SetFillStyle(1001);

  double yMin = 1, yMax = 0, relativeMin = -5, relativeMax = 5;
  int withinCount = 0;
  for (int i = 0; i < count; ++i)
  {
    const auto &bin = bins[i];
    // Equal-width categorical pT bins and small marker offsets match the reference layout.
    nominal.SetPoint(i, i - 0.06, bin.nominal);
    nominal.SetPointError(i, 0, bin.nominalError);
    band.SetPoint(i, i, bin.nominal);
    band.SetPointError(i, 0.5, bin.nominal * targetPercent / 100.0);
    const bool inTarget = std::abs(100 * (bin.fit / bin.nominal - 1)) <= targetPercent + 1e-12;
    withinCount += inTarget;
    auto &fitGraph = inTarget ? within : outside;
    auto &relativeGraph = inTarget ? relativeWithin : relativeOutside;
    const int point = fitGraph.GetN();
    fitGraph.SetPoint(point, i + 0.06, bin.fit);
    fitGraph.SetPointError(point, 0, bin.fitError);
    // Fit-only error divided by nominal; no independent-error assumption for this shared-data comparison.
    const double difference = 100 * (bin.fit / bin.nominal - 1);
    const double error = 100 * bin.fitError / bin.nominal;
    relativeGraph.SetPoint(point, i, difference);
    relativeGraph.SetPointError(point, 0, error);
    yMin = std::min(yMin, std::min(bin.nominal - bin.nominalError, bin.fit - bin.fitError));
    yMax = std::max(yMax, std::max(bin.nominal + bin.nominalError, bin.fit + bin.fitError));
    const double relativeSigma = 100 * bin.nominalError / bin.nominal;
    relativeMin = std::min(relativeMin, std::min(difference - error, -relativeSigma));
    relativeMax = std::max(relativeMax, std::max(difference + error, relativeSigma));
  }

  // Keep every uncertainty visible, including the low-pT relative errors.
  yMin = std::max(0.0, yMin - 0.02);
  yMax += 0.035;
  relativeMin = 5 * std::floor((relativeMin - 2) / 5);
  relativeMax = 5 * std::ceil((relativeMax + 2) / 5);
  TH1D topFrame("bfraction_top_frame", "", count, -0.5, count - 0.5);
  TH1D bottomFrame("bfraction_relative_frame", "", count, -0.5, count - 0.5);
  styleFrame(topFrame, "Nonprompt fraction B");
  styleFrame(bottomFrame, "Difference / nominal (%)");
  bottomFrame.SetMinimum(relativeMin);
  bottomFrame.SetMaximum(relativeMax);
  topFrame.SetMinimum(yMin);
  topFrame.SetMaximum(yMax);
  topFrame.GetXaxis()->SetLabelSize(0);
  bottomFrame.GetXaxis()->SetTitle("p_{T} (GeV/c)");
  bottomFrame.GetXaxis()->SetTitleOffset(1.25);
  for (int i = 0; i < count; ++i)
  {
    topFrame.GetXaxis()->SetBinLabel(i + 1, Form("%g-%g", bins[i].low, bins[i].high));
    bottomFrame.GetXaxis()->SetBinLabel(i + 1, Form("%g-%g", bins[i].low, bins[i].high));
  }
  topFrame.GetXaxis()->LabelsOption("h");
  bottomFrame.GetXaxis()->LabelsOption("h");
  // Lower panel only: nominal statistical +/-1 sigma replaces the dashed limits.
  // Draw this band behind the original green +/-targetPercent band.
  TGraphErrors relativeSigmaBand(count);
  for (int i = 0; i < count; ++i)
  {
    relativeSigmaBand.SetPoint(i, i, 0);
    relativeSigmaBand.SetPointError(i, 0.5, 100.0 * bins[i].nominalError / bins[i].nominal);
  }
  relativeSigmaBand.SetFillColor(TColor::GetColor("#f6b26b"));
  relativeSigmaBand.SetLineColor(TColor::GetColor("#f6b26b"));
  relativeSigmaBand.SetFillStyle(1001);

  TCanvas canvas("bfraction_comparison_canvas", "Production b-fraction vs nominal", 1300, 1000);
  canvas.SetCanvasSize(1300, 1000);
  canvas.SetFillColor(kWhite);
  // Match archived_260908/manual/drawing/b_frac_comparison: upper/lower = 0.51/0.30.
  TPad upper("bfraction_upper", "", 0, 0.44, 1, 0.95);
  TPad lower("bfraction_lower", "", 0, 0.14, 1, 0.44);
  for (TPad *pad : {&upper, &lower})
  {
    pad->SetFillColor(kWhite);
    pad->SetLeftMargin(0.09);
    pad->SetRightMargin(0.025);
    pad->SetTopMargin(0.02);
    pad->SetTicks(1, 1);
    pad->Draw();
  }
  upper.SetBottomMargin(0.025);
  lower.SetBottomMargin(0.26);
  upper.cd();
  topFrame.Draw("AXIS");
  band.Draw("2 SAME");
  drawGrid(topFrame, 0.05);
  nominal.Draw("P SAME");
  if (within.GetN()) within.Draw("P SAME");
  if (outside.GetN()) outside.Draw("P SAME");
  TLegend legend(0.11, 0.67, 0.57, 0.96);
  legend.SetBorderSize(0);
  legend.SetFillColor(kWhite);
  legend.SetTextFont(43);
  legend.SetTextSize(19);
  legend.AddEntry(&band, Form("Nominal #pm%g%% target", targetPercent), "f");
  legend.AddEntry(&relativeSigmaBand, "Lower panel: nominal #pm1#sigma statistical band", "f");
  legend.AddEntry(&nominal, mid ? "Nominal ver3 #pm statistical error" : "Nominal ver3 #pm statistical error", "pe");
  if (within.GetN()) legend.AddEntry(&within, "Within target #pm fit error", "pe");
  if (outside.GetN()) legend.AddEntry(&outside, "Outside target #pm fit error", "pe");
  legend.Draw();
  TLatex rapidity;
  rapidity.SetNDC();
  rapidity.SetTextFont(43);
  rapidity.SetTextSize(20);
  rapidity.SetTextAlign(31);
  // KST - Labels follow the selected region; neither nominal drawing reference is overwritten.
  rapidity.DrawLatex(0.96, 0.925, mid ? "OO, |y| < 1.6" : "OO, 1.6 < |y| < 2.4");
  upper.RedrawAxis();

  lower.cd();
  bottomFrame.Draw("AXIS");
  relativeSigmaBand.Draw("2 SAME");
  // Relative agreement target; nominal statistical errors appear above.
  TGraphErrors relativeBand(count);
  for (int i = 0; i < count; ++i)
  {
    relativeBand.SetPoint(i, i, 0);
    relativeBand.SetPointError(i, 0.5, targetPercent);
  }
  relativeBand.SetFillColor(green);
  relativeBand.SetLineColor(green);
  relativeBand.Draw("2 SAME");
  drawGrid(bottomFrame, 5);
  TLine zero(-0.5, 0, count - 0.5, 0);
  zero.SetLineColor(kBlack);
  zero.Draw();
  if (relativeWithin.GetN()) relativeWithin.Draw("P SAME");
  if (relativeOutside.GetN()) relativeOutside.Draw("P SAME");
  TLatex label;
  label.SetTextFont(43);
  label.SetTextSize(17);
  label.SetTextAlign(21);
  for (int i = 0; i < count; ++i)
  {
    const double difference = 100 * (bins[i].fit / bins[i].nominal - 1);
    label.DrawLatex(i, difference + 100 * bins[i].fitError / bins[i].nominal + 1,
                    Form("%+.2f%%", difference));
  }
  lower.RedrawAxis();

  canvas.cd();
  TLatex heading;
  heading.SetNDC();
  heading.SetTextFont(43);
  heading.SetTextSize(25);
  heading.SetTextAlign(23);
  heading.DrawLatex(0.53, 0.985, Form("%s fit: %g%% agreement target", mid ? "Mid" : "Forward", targetPercent));
  TLatex note;
  note.SetNDC();
  note.SetTextFont(43);
  note.SetTextSize(16);
  note.DrawLatex(0.035, 0.105, "B and fit errors are read directly from ROOT files; lower errors = 100 #times fit error / nominal.");
  note.DrawLatex(0.035, 0.080, Form("Target: 100 #times |B/B_{nominal} - 1| #leq %g%% in every bin. Green bands retain this target; orange (lower panel) shows nominal #pm1#sigma_{stat}.", targetPercent));
  note.DrawLatex(0.035, 0.055, mid ? "Model: fit2d.C. Reference: nominal ver3, bfraction_ctau_fit.xlsx (2026-10-31) (nominal_bfraction.h). Colors show agreement; validate fit quality separately." : "Model and ctau settings: fit2d.C. Reference: nominal ver3, bfraction_ctau_fit.xlsx (2026-10-31) (nominal_bfraction.h). Colors show agreement only; validate fit quality separately.");
  canvas.Modified();
  canvas.Update();
  // KST - Separate mid filenames preserve the existing forward plot and agreement table.
  const TString stem = output + (mid ? "/bfraction_mid_production_vs_nominal" : "/bfraction_production_vs_nominal");
  canvas.Print(stem + ".pdf");
  std::ifstream saved((stem + ".pdf").Data(), std::ios::binary | std::ios::ate);
  if (!saved || saved.tellg() <= 0)
    throw std::runtime_error("Failed to write PDF: " + std::string(stem.Data()) + ".pdf");
  std::ofstream table((output + (mid ? "/bfraction_mid_agreement.csv" : "/bfraction_agreement.csv")).Data());
  if (!table) throw std::runtime_error("Cannot write agreement table");
  table << std::setprecision(17)
        << "pt_low,pt_high,nominal,nominal_error,fit,fit_error,difference_percent,target_percent,agreement_pass,difference_sigma,tier,minimizer\n";
  for (const auto &bin : bins) {
    const double difference = 100 * (bin.fit/bin.nominal - 1);
    table << bin.low << ',' << bin.high << ',' << bin.nominal << ',' << bin.nominalError
          << ',' << bin.fit << ',' << bin.fitError << ',' << difference << ',' << targetPercent
          << ',' << (std::abs(difference) <= targetPercent + 1e-12 ? 1 : 0);
    const auto &ref=nominal_bfraction::find(bin.low,bin.high,mid?0.:1.6,mid?1.6:2.4);
    table << ',' << nominal_bfraction::differenceSigma(bin.fit,ref) << ',' << nominal_bfraction::tier(bin.fit,ref) << ',' << "Minuit2" << '\n';
  }
  std::cout << "Input: " << input << '\n'
            << withinCount << "/" << count << " bins within " << targetPercent << "% of nominal\n"
            << "Saved: " << stem << ".pdf\n";
}
