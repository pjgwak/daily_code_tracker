// From this folder: root -l -b -q 'plot_bfraction_nominal_with_fwd_syst.C()'
// Optional arguments: ("/path/to/fit2d", "/path/to/output")
// Shared plot implementation: regional nominal values and errors come from nominal_bfraction.h.
// Assign the full current forward shift as a symmetric systematic uncertainty:
//   syst = abs(B_fwd - B_nominal).
// Statistical errors are shown separately; the forward fit error is not added
// to the systematic shift. Boxes span 80% of each physical pT bin.
#include "plot_bfraction_vs_nominal.C"

void plot_bfraction_nominal_with_syst(const char *rootDir, const char *outputDir, bool mid)
{
  using namespace bfraction_comparison;
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  const TString directory = macroDirectory();
  const TString input = rootDir[0] ? TString(rootDir)
      : directory + (mid ? "/../../roots/y0.00_1.60/fit2d" : "/../../roots/y1.60_2.40/fit2d");
  const TString output = outputDir[0] ? TString(outputDir) : directory;
  const auto bins = readBins(input, mid);
  const char *region = mid ? "mid" : "fwd";
  if (gSystem->AccessPathName(output) && gSystem->mkdir(output, true) != 0)
    throw std::runtime_error("Cannot create output directory: " + std::string(output.Data()));
  if (gSystem->AccessPathName(output, kWritePermission))
    throw std::runtime_error("Output directory is not writable: " + std::string(output.Data()));

  const int count = static_cast<int>(bins.size());
  TGraphErrors nominal(count), systematic(count);
  styleGraph(nominal, kBlack, 20);
  const Color_t blue = TColor::GetColor("#bfd7ef");
  systematic.SetFillColor(blue);
  systematic.SetFillStyle(1001);
  systematic.SetLineColor(TColor::GetColor("#226fba"));
  double yMin = 1, yMax = 0;
  for (int i = 0; i < count; ++i)
  {
    const auto &bin = bins[i];
    const double x = 0.5 * (bin.low + bin.high);
    const double syst = std::abs(bin.fit - bin.nominal);
    nominal.SetPoint(i, x, bin.nominal);
    nominal.SetPointError(i, 0, bin.nominalError);
    systematic.SetPoint(i, x, bin.nominal);
    systematic.SetPointError(i, 0.4 * (bin.high - bin.low), syst);
    const double extent = std::max(syst, bin.nominalError);
    yMin = std::min(yMin, bin.nominal - extent);
    yMax = std::max(yMax, bin.nominal + extent);
    std::cout << Form("pT %g-%g: nominal %.8f +/- %.8f (stat) +/- %.8f (%s syst), relative syst %.2f%%",
                      bin.low, bin.high, bin.nominal, bin.nominalError,
                      syst, region, 100 * syst / bin.nominal) << '\n';
  }

  TCanvas canvas(Form("bfraction_nominal_%s_syst_canvas", region), "Nominal with fit systematic", 1000, 750);
  canvas.SetCanvasSize(1000, 750);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.035);
  canvas.SetBottomMargin(0.13);
  canvas.SetTopMargin(0.08);
  canvas.SetTicks(1, 1);
  TH1D frame(Form("bfraction_nominal_%s_syst_frame", region), "", 100, bins.front().low, bins.back().high);
  styleFrame(frame, "Nonprompt fraction B");
  frame.GetXaxis()->SetTitle("p_{T} (GeV/c)");
  frame.GetXaxis()->SetTitleOffset(1.1);
  frame.GetYaxis()->SetTitleOffset(1.3);
  frame.SetMinimum(std::max(0.0, yMin - 0.035));
  frame.SetMaximum(yMax + 0.13);
  frame.Draw("AXIS");
  systematic.Draw("2 SAME");
  nominal.Draw("PZ SAME");
  TLegend legend(0.15, 0.72, 0.68, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextFont(43);
  legend.SetTextSize(22);
  legend.AddEntry(&nominal, "Nominal #pm statistical uncertainty", "pe");
  legend.AddEntry(&systematic, mid ? "Mid systematic: |B_{mid} - B_{nominal}|"
                                   : "Fwd systematic: |B_{fwd} - B_{nominal}|", "f");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextFont(43);
  label.SetTextSize(23);
  label.SetTextAlign(31);
  label.DrawLatex(0.96, 0.945, mid ? "OO, |y| < 1.6" : "OO, 1.6 < |y| < 2.4");
  canvas.RedrawAxis();
  canvas.Modified();
  canvas.Update();
  const TString pdf = output + Form("/bfraction_nominal_with_%s_syst.pdf", region);
  canvas.Print(pdf);
  std::ifstream saved(pdf.Data(), std::ios::binary | std::ios::ate);
  if (!saved || saved.tellg() <= 0)
    throw std::runtime_error("Failed to write PDF: " + std::string(pdf.Data()));
  std::cout << "Input: " << input << "\nSaved: " << pdf << '\n';
}

// Keep the original forward entry point and output name.
void plot_bfraction_nominal_with_fwd_syst(const char *rootDir = "", const char *outputDir = "")
{
  plot_bfraction_nominal_with_syst(rootDir, outputDir, false);
}
