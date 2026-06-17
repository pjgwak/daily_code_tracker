#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TBox.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"

namespace
{
TString macroDir()
{
  TString path = __FILE__;
  if (!gSystem->IsAbsoluteFileName(path))
    path = TString(gSystem->WorkingDirectory()) + "/" + path;
  return gSystem->DirName(path);
}

TString outputPath(const char *fileName)
{
  return macroDir() + "/" + fileName;
}
} // namespace

void raw_yield_fwd()
{
  TStopwatch __execTimer;
  __execTimer.Start();
  struct __ExecTimerGuard
  {
    TStopwatch *sw;
    explicit __ExecTimerGuard(TStopwatch *timer) : sw(timer) {}
    ~__ExecTimerGuard()
    {
      sw->Stop();
      sw->Print();
    }
  } __execTimerGuard(&__execTimer);

  // =====  CMS plot style =====
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // ===== OO raw yields (12 bins) =====
  // ctau subrange method
  const int n_old = 12;

  double pt_low_old[n_old] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                              7.0, 8.0, 9.0, 10.0, 12.0, 14.0};

  double pt_high_old[n_old] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                               8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  // nominal raw yields
  double val_pp_old[n_old] = {
      25932.8, 24533.1, 19598.2, 13779.8,
      10321.4, 6640.5, 4467.7, 2796.3,
      1771.4, 2060.3, 943.6, 952.2};

  // statistical errors
  double stat_pp_old[n_old] = {
      660.2, 474.7, 363.4, 260.2,
      191.1, 134.8, 99.4, 72.3,
      56.6, 59.0, 39.3, 38.0};

  // systematics not provided
  double sys_pp_old[n_old] = {
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0};

  double x_old[n_old], ex_old[n_old];
  for (int i = 0; i < n_old; ++i)
  {
    const double width = pt_high_old[i] - pt_low_old[i];
    x_old[i] = 0.5 * (pt_low_old[i] + pt_high_old[i]);
    ex_old[i] = 0.5 * width;
    val_pp_old[i] /= width;
    stat_pp_old[i] /= width;
    sys_pp_old[i] /= width;
  }

  // ===== OO 1.6 < |y| < 2.4 (2d fit) =====
  // y = 1.60-2.40
  const int NooY16to24 = 12;
  double pt_low_oo_y16to24[NooY16to24] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0};
  double pt_high_oo_y16to24[NooY16to24] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  double val_oo_y16to24[NooY16to24] = {24610.13871, 23941.55577, 19177.37525, 14026.88170, 10062.95819, 6557.66420, 4394.41580, 2777.29863, 1731.29077, 1991.43129, 889.11726, 921.71206};
  double stat_oo_y16to24[NooY16to24] = {333.77081, 258.94772, 202.38711, 157.12752, 115.80535, 84.44122, 65.06869, 49.94252, 38.26061, 39.18428, 26.13233, 25.60973};

  double sys_oo_y16to24[NooY16to24];
  for (int i = 0; i < NooY16to24; ++i)
    sys_oo_y16to24[i] = 0.0; // sys=0

  double x_oo_y16to24[NooY16to24], ex_oo_y16to24[NooY16to24];
  for (int i = 0; i < NooY16to24; ++i)
  {
    const double width = pt_high_oo_y16to24[i] - pt_low_oo_y16to24[i];
    x_oo_y16to24[i] = 0.5 * (pt_low_oo_y16to24[i] + pt_high_oo_y16to24[i]);
    ex_oo_y16to24[i] = 0.5 * width;
    val_oo_y16to24[i] /= width;
    stat_oo_y16to24[i] /= width;
    sys_oo_y16to24[i] /= width;
  }

  double yMin = 1.0e30;
  double yMax = 0.0;
  for (int i = 0; i < n_old; ++i)
  {
    const double yLow = val_pp_old[i] - stat_pp_old[i];
    const double yHigh = val_pp_old[i] + stat_pp_old[i];
    if (yLow > 0.0 && yLow < yMin)
      yMin = yLow;
    if (yHigh > yMax)
      yMax = yHigh;
  }
  for (int i = 0; i < NooY16to24; ++i)
  {
    const double yLow = val_oo_y16to24[i] - stat_oo_y16to24[i];
    const double yHigh = val_oo_y16to24[i] + stat_oo_y16to24[i];
    if (yLow > 0.0 && yLow < yMin)
      yMin = yLow;
    if (yHigh > yMax)
      yMax = yHigh;
  }

  // ===== Graphs =====
  // pp_old: sys band
  TGraphAsymmErrors *gSyspp_old = new TGraphAsymmErrors(n_old, x_old, val_pp_old, ex_old, ex_old, sys_pp_old, sys_pp_old);
  gSyspp_old->SetFillStyle(0);
  gSyspp_old->SetLineColor(kRed + 1);
  gSyspp_old->SetLineWidth(2);

  // pp_old: stat
  TGraphErrors *gStatpp_old = new TGraphErrors(n_old, x_old, val_pp_old, ex_old, stat_pp_old);
  gStatpp_old->SetMarkerStyle(24);
  gStatpp_old->SetMarkerSize(1.5);
  gStatpp_old->SetMarkerColor(kRed + 1);
  gStatpp_old->SetLineColor(kRed + 1);
  gStatpp_old->SetLineWidth(2);

  // OO 1.6 < |y| < 2.4: sys (currently 0)
  TGraphAsymmErrors *gSysOOY16to24 = new TGraphAsymmErrors(NooY16to24, x_oo_y16to24, val_oo_y16to24, ex_oo_y16to24, ex_oo_y16to24, sys_oo_y16to24, sys_oo_y16to24);
  gSysOOY16to24->SetFillStyle(0); // no fill
  gSysOOY16to24->SetLineColor(kAzure + 2);
  gSysOOY16to24->SetLineWidth(2);

  // OO 1.6 < |y| < 2.4: stat
  TGraphErrors *gStatOOY16to24 = new TGraphErrors(NooY16to24, x_oo_y16to24, val_oo_y16to24, ex_oo_y16to24, stat_oo_y16to24);
  gStatOOY16to24->SetMarkerStyle(20);
  gStatOOY16to24->SetMarkerSize(1.4);
  gStatOOY16to24->SetMarkerColor(kAzure + 2);
  gStatOOY16to24->SetLineColor(kAzure + 2);
  gStatOOY16to24->SetLineWidth(2);

  // ===== Draw =====
  TCanvas *c1 = new TCanvas("c1", "pp vs pp_old", 800, 800);
  c1->SetTitle("");
  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.SetLogy();
  pad1.Draw();
  pad1.cd();
  gStyle->SetOptStat(0);
  gSyspp_old->SetTitle("");
  gSyspp_old->Draw("A2");
  gSyspp_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gSyspp_old->GetYaxis()->SetTitle("dN/dp_{T}");
  gSyspp_old->GetXaxis()->CenterTitle();
  gSyspp_old->GetYaxis()->CenterTitle();
  gSyspp_old->GetXaxis()->SetTitleSize(0.05);
  gSyspp_old->GetYaxis()->SetTitleSize(0.05);
  gSyspp_old->GetXaxis()->SetLabelSize(0.042);
  gSyspp_old->GetYaxis()->SetLabelSize(0.042);
  gSyspp_old->GetYaxis()->SetTitleOffset(1.25);
  gSyspp_old->GetXaxis()->SetLimits(pt_low_old[0], pt_high_old[n_old - 1]);
  gSyspp_old->SetMinimum(0.8 * yMin);
  gSyspp_old->SetMaximum(4 * yMax);

  gSyspp_old->Draw("A2");
  gStatpp_old->Draw("P SAME");
  gStatOOY16to24->Draw("P SAME");

  // ------------------------------------------------------------------
  // labels (match ctau_pr.C style/positions)
  // ------------------------------------------------------------------
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
    // tx.SetTextColor(kRed + 1);
    tx.DrawLatex(0.19, 0.865, "1.6 < |y| < 2.4");
  }

  TLegend *leg = new TLegend(0.7, 0.8, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.024);
  leg->SetEntrySeparation(0.01);
  leg->AddEntry(gStatOOY16to24, "OO (2D fit)", "lp");
  leg->AddEntry(gStatpp_old, "OO (Nominal)", "lp");
  leg->Draw();

  c1->SaveAs(outputPath("raw_yield_OO_fwd.pdf"));
  // c1->SaveAs(outputPath("plot_OO_fwd.png"));
}