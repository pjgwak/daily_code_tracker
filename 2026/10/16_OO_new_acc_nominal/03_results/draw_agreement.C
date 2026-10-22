#include "plot_support.h"
#include <map>
#include <tuple>
#include <string>
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <array>
#include <memory>
#include <vector>

// Run through tools/results.py; the input table binds current ROOTs to frozen nominal values.
// Raw B from selected saved ROOTs only; no MC window transport or fitted-B error bars.
void draw_agreement(const char *inputPath, const char *outputDirectory)
{
  TString source = __FILE__;
  if (!gSystem->IsAbsoluteFileName(source)) source = TString(gSystem->WorkingDirectory()) + "/" + source;
  const TString stylePath = TString(gSystem->DirName(source)) + "/style/rootlogon.C";
  const char *style = stylePath.Data();
  if (gSystem->AccessPathName(style)) throw std::runtime_error("Missing rootlogon.C");
  Int_t error = 0;
  gROOT->Macro(style, &error);
  if (error) throw std::runtime_error("Failed to apply rootlogon.C");
  gROOT->SetBatch(kTRUE);
  gStyle->SetPaperSize(28, 17.2);
  const TString base = outputDirectory;
  std::ifstream input(inputPath);
  if (!input) throw std::runtime_error("Cannot read selected ROOT input table");
  struct Input { TString path; double nominal, error; };
  std::map<std::tuple<int,double,double>, Input> selected;
  int region; double low, high, nominal, nominalError; std::string path;
  while (input >> region >> low >> high >> nominal >> nominalError >> std::quoted(path)) {
    if ((region!=0 && region!=1) || !std::isfinite(nominal) || nominal<=0 ||
        !std::isfinite(nominalError) || nominalError<=0 || !(high>low))
      throw std::runtime_error("Invalid selected bin or nominal reference");
    if (!selected.emplace(std::make_tuple(region,low,high),Input{TString(path),nominal,nominalError}).second)
      throw std::runtime_error("Duplicate selected bin");
  }
  if (!input.eof()) throw std::runtime_error("Malformed selected ROOT input table");
  if (selected.size()!=20) throw std::runtime_error("Expected exactly20 selected ROOTs");
  std::ofstream table((base+"/agreement.csv").Data());
  table << std::setprecision(17) << "region,pt_low,pt_high,B_raw,B_raw_fit_error,B_nominal,nominal_stat_error,signed_nominal_stat_units,signed_difference_percent,tier,source_root\n";
  int count29=0,countSept=0,countSigma=0; double maximumPercent=0.,maximumSigma=0.;
  const double septemberMaximumPercent=2.91656031733;
  struct Point { double low, high, sigma, percent; Color_t color; };
  std::array<std::vector<Point>, 2> data;
  for (int mid = 0; mid < 2; ++mid) {
    for (const auto &entry : selected) {
      if (std::get<0>(entry.first)!=mid) continue;
      const nominal_bfraction::Bin ref{std::get<1>(entry.first),std::get<2>(entry.first),
          entry.second.nominal,entry.second.error,static_cast<bool>(mid)};
      const TString fitPath=entry.second.path;
      TFile policyFile(fitPath, "READ");
      auto *policyResult = policyFile.Get<RooFitResult>("modelResult");
      if (!policyResult || policyResult->floatParsFinal().find("promptResolutionScale"))
        throw std::runtime_error("Forbidden common resolution scale in selected ROOT");
      auto *commonRule = policyFile.Get<TParameter<double>>("rule_common_scale");
      if (commonRule && commonRule->GetVal()!=0.)
        throw std::runtime_error("Forbidden common-scale rule in selected ROOT");
      policyFile.Close();
      const auto fit = oo_drawing::readFit(fitPath);
      const double sigma = (fit.b - ref.value) / ref.error;
      if (!std::isfinite(sigma) || std::abs(sigma) > 1.1)
        throw std::runtime_error("Agreement point outside displayed range");
      // Preserve the original agreement-tier colors for markers and percentages.
      const std::string tier = nominal_bfraction::tier(fit.b, ref);
      const Color_t color = tier == "BEST" ? kGreen+2 : tier == "MAGINOT" ? kAzure+2
          : tier.find("PERCENT") == 0 ? kOrange+7 : kRed+1;
      data[mid].push_back({ref.low, ref.high, sigma, nominal_bfraction::differencePercent(fit.b, ref), color});
      const double percent=nominal_bfraction::differencePercent(fit.b,ref);
      countSigma += std::abs(sigma)<=0.414;
      count29 += std::abs(percent)<2.9;
      countSept += std::abs(percent)<septemberMaximumPercent;
      maximumPercent=std::max(maximumPercent,std::abs(percent));
      maximumSigma=std::max(maximumSigma,std::abs(sigma));
      table << (mid?"mid":"fwd") << ',' << ref.low << ',' << ref.high << ',' << fit.b << ',' << fit.bError << ',' << ref.value << ',' << ref.error << ',' << sigma << ',' << percent << ',' << tier << ',' << fitPath.Data() << '\n';
    }
  }
  if (data[0].size()!=12 || data[1].size()!=8) throw std::runtime_error("Expected 12 forward and 8 mid bins");
  const std::array<Int_t, 3> fills = {{TColor::GetColor("#eeeeee"),
      TColor::GetColor("#dce9f5"), TColor::GetColor("#c7e3d3")}};
  const std::array<double, 3> widths = {{1.0, 0.414, 0.138}};
  // Keep primitives alive through PDF export.
  std::vector<std::unique_ptr<TObject>> objects;
  TCanvas canvas("nominal_agreement", "Nominal agreement", 1400, 860);
  canvas.SetCanvasSize(1400, 860);
  canvas.SetFillColor(kWhite);
  for (int mid = 0; mid < 2; ++mid) {
    canvas.cd();
    auto *pad = new TPad(Form("agreement_%d", mid), "", mid ? .505 : .005, .10,
                         mid ? .995 : .495, .915);
    pad->SetLeftMargin(.17); pad->SetRightMargin(.06);
    pad->SetBottomMargin(.13); pad->SetTopMargin(.06);
    pad->SetTicks(1, 1); pad->Draw(); pad->cd();
    const int count = data[mid].size();
    auto frame = std::make_unique<TH2D>(Form("agreement_frame_%d", mid), "", 10, -1.1, 1.1, count, 0, count);
    frame->SetDirectory(nullptr); frame->SetStats(false);
    frame->GetXaxis()->SetTitle("(B_{raw} - B_{nom}) / #sigma_{stat,nom}");
    frame->GetYaxis()->SetTitle("");
    for (auto *axis : {frame->GetXaxis(), frame->GetYaxis()}) {
      axis->SetTitleFont(43); axis->SetLabelFont(43);
      axis->SetTitleSize(26); axis->SetLabelSize(22); axis->CenterTitle();
    }
    frame->GetXaxis()->SetTitleOffset(1.3); frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetXaxis()->SetNdivisions(505); frame->GetXaxis()->SetTickLength(.018);
    frame->GetYaxis()->SetTickLength(0);
    for (int i = 0; i < count; ++i)
      frame->GetYaxis()->SetBinLabel(count-i, Form("%g#font[122]{-}%g", data[mid][i].low, data[mid][i].high));
    frame->Draw("AXIS");
    for (int j = 0; j < 3; ++j) {
      auto box = std::make_unique<TBox>(-widths[j], 0, widths[j], count);
      box->SetFillColor(fills[j]); box->SetFillStyle(1001); box->SetLineWidth(0); box->Draw();
      objects.push_back(std::move(box));
    }
    TLine grid; grid.SetLineColor(kWhite);
    for (int i = 1; i < count; ++i) grid.DrawLine(-1.1, i, 1.1, i);
    TLine zero; zero.SetLineColor(TColor::GetColor("#808080")); zero.SetLineStyle(2);
    zero.DrawLine(0, 0, 0, count);
    TLatex numbers; numbers.SetTextFont(43);
    numbers.SetTextSize(21); numbers.SetTextAlign(12);
    for (int i = 0; i < count; ++i) {
      const double y = count - i - .5;
      auto point = std::make_unique<TGraph>();
      point->SetPoint(0, data[mid][i].sigma, y);
      point->SetMarkerColor(data[mid][i].color);
      point->SetMarkerStyle((mid && data[mid][i].low!=8) ? 33 : 20); point->SetMarkerSize(1.25);
      point->Draw("P SAME");
      objects.push_back(std::move(point));
      numbers.SetTextColor(data[mid][i].color);
      // Keep labels inside the frame and use the same visual offset in both panels.
      const bool labelOnLeft = data[mid][i].sigma > .60;
      numbers.SetTextAlign(labelOnLeft ? 32 : 12);
      numbers.DrawLatex(data[mid][i].sigma + (labelOnLeft ? -.055 : .055),
                        y + .012*count, Form("%+.2f%%", data[mid][i].percent));
    }
    TLatex heading; heading.SetNDC(); heading.SetTextFont(43);
    heading.SetTextSize(25); heading.SetTextAlign(12);
    heading.DrawLatex(.17, .972, mid ? "|y| < 1.6" : "1.6 < |y| < 2.4");
    TLatex ylabel; ylabel.SetNDC(); ylabel.SetTextFont(43);
    ylabel.SetTextSize(26); ylabel.SetTextAlign(22); ylabel.SetTextAngle(90);
    ylabel.DrawLatex(.040, .535, "p_{T} (GeV/c)");
    pad->RedrawAxis();
    objects.push_back(std::move(frame));
  }
  canvas.cd();
  TLatex header; header.SetNDC(); header.SetTextAlign(11);
  header.SetTextFont(63); header.SetTextSize(34); header.DrawLatex(.0883, .945, "CMS");
  header.SetTextFont(53); header.SetTextSize(26); header.DrawLatex(.161, .945, "Preliminary");
  header.SetTextAlign(31); header.SetTextFont(43); header.SetTextSize(25);
  header.DrawLatex(.9656, .945, "OO  #sqrt{s_{NN}} = 5.36 TeV  (9 nb^{-1})");
  TLatex context; context.SetNDC(); context.SetTextFont(43); context.SetTextSize(18); context.SetTextAlign(22);
  context.DrawLatex(.50,.919,Form("Raw B: %d/20 within 0.414 #sigma_{stat,nom}; no fitted common resolution scale.",countSigma));
  TLatex note; note.SetNDC(); note.SetTextFont(43); note.SetTextSize(15); note.SetTextAlign(22);
  note.DrawLatex(.50,.014,"Diamonds: 7 corrected mid fits; windows changed only for mid 7-8 and 9-10 in this revision. No postfit B transport.");
  std::array<TBox, 3> samples;
  TLegend legend(.29, .040, .73, .079);
  legend.SetNColumns(3); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.SetTextFont(43); legend.SetTextSize(21);
  legend.SetMargin(.18); legend.SetColumnSeparation(.025);
  for (int j = 0; j < 3; ++j) {
    const int k = 2-j;
    samples[j].SetFillColor(fills[k]); samples[j].SetLineWidth(0);
    legend.AddEntry(&samples[j], Form("#pm%g #sigma_{stat,nom}", widths[k]), "f");
  }
  legend.Draw();
  canvas.Modified(); canvas.Update();
  canvas.Print(base + "/nominal_agreement.pdf");
  // The CLI renders the PDF with pdftoppm; ROOT image plugins may exit early in batch mode.
  std::ofstream summary((base+"/agreement_summary.json").Data());
  summary << std::setprecision(17) << "{\"count_under_2p9_percent\":" << count29 << ",\"count_under_september_max\":" << countSept << ",\"september_max_percent\":" << septemberMaximumPercent << ",\"max_abs_raw_percent\":" << maximumPercent << ",\"max_abs_nominal_stat_units\":" << maximumSigma << ",\"MC_window_transport_applied\":false,\"significance_interpretation\":false}\n";
}
