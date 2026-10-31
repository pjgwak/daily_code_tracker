#include "../latest_fit_results.h"
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

// Run from the analysis directory:
// root -n -l -b -q 'drawing/b_frac_comparison/plot_nominal_agreement.C()'
// PDF-only export; use the shared reader to validate current fit results.
void plot_nominal_agreement()
{
  const char *style = "/data/users/pjgwak/input_files/rootlogon.C";
  if (gSystem->AccessPathName(style)) throw std::runtime_error("Missing rootlogon.C");
  Int_t error = 0;
  gROOT->Macro(style, &error);
  if (error) throw std::runtime_error("Failed to apply rootlogon.C");
  gROOT->SetBatch(kTRUE);
  gStyle->SetPaperSize(28, 17.2);
  TString macro = __FILE__;
  if (!gSystem->IsAbsoluteFileName(macro)) macro = TString(gSystem->WorkingDirectory()) + "/" + macro;
  const TString base = TString(gSystem->DirName(macro)) + "/../..";
  struct Point { double low, high, sigma, percent; Color_t color; };
  std::array<std::vector<Point>, 2> data;
  for (int mid = 0; mid < 2; ++mid) {
    for (int i = 0; i < (mid ? 8 : 12); ++i) {
      const auto &ref = mid ? nominal_bfraction::midBins()[i] : nominal_bfraction::bins()[i];
      const TString tag = mid ? "y0.00_1.60" : "y1.60_2.40";
      const auto fit = oo_drawing::readFit(Form("%s/roots/%s/fit2d/fit2d_result_%s_pt%.2f_%.2f.root",
          base.Data(), tag.Data(), tag.Data(), ref.low, ref.high));
      const double sigma = (fit.b - ref.value) / ref.error;
      if (!std::isfinite(sigma) || std::abs(sigma) > 1.1)
        throw std::runtime_error("Agreement point outside displayed range");
      // Preserve the original agreement-tier colors for markers and percentages.
      const std::string tier = nominal_bfraction::tier(fit.b, ref);
      const Color_t color = tier == "BEST" ? kGreen+2 : tier == "MAGINOT" ? kAzure+2
          : tier.find("PERCENT") == 0 ? kOrange+7 : kRed+1;
      data[mid].push_back({ref.low, ref.high, sigma, nominal_bfraction::differencePercent(fit.b, ref), color});
    }
  }
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
    frame->GetXaxis()->SetTitle("(B - B_{nom}) / #sigma_{stat,nom}");
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
      point->SetMarkerStyle(20); point->SetMarkerSize(1.25);
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
}
