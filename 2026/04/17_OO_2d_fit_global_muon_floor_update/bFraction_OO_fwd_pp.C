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

void bFraction_OO_fwd_pp()
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

  // ===== pp_bFraction (12 bins) =====
  // not pp but OO, ctau subrange method
  const int n_old = 12;

  double pt_low_old[n_old] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                              7.0, 8.0, 9.0, 10.0, 12.0, 14.0};

  double pt_high_old[n_old] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                               8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  // f_B = 1 - f_prompt
  double val_pp_old[n_old] = {
      0.103,
      0.143,
      0.134,
      0.186,
      0.171,
      0.193,

      0.255,
      0.265,
      0.304,
      0.317,
      0.374,
      0.400 
  };

  // statistical errors: unchanged
  double stat_pp_old[n_old] = {
      0.021, 0.012, 0.009, 0.008,
      0.007, 0.008,
      0.019, 0.013, 0.017, 0.020, 0.022, 0.023};

  // systematics not provided
  double sys_pp_old[n_old] = {
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0};

  double x_old[n_old], ex_old[n_old];
  for (int i = 0; i < n_old; ++i)
  {
    x_old[i] = 0.5 * (pt_low_old[i] + pt_high_old[i]);
    ex_old[i] = 0.5 * (pt_high_old[i] - pt_low_old[i]);
  }

  // ===== OO 1.6 < |y| < 2.4 (2d fit) =====
  // y = 1.60-2.40
  const int NooY16to24 = 12;
  double pt_low_oo_y16to24[NooY16to24] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0};
  double pt_high_oo_y16to24[NooY16to24] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 20.0};

  double val_oo_y16to24[NooY16to24] = {0.12834, 0.11775, 0.14149, 0.18799, 0.18458, 0.20018, 0.24011, 0.27390, 0.28881, 0.31005, 0.37460, 0.38531};
  double stat_oo_y16to24[NooY16to24] = {0.01688, 0.01056, 0.00816, 0.00776, 0.00793, 0.00891, 0.01082, 0.01329, 0.01639, 0.01495, 0.02346, 0.02234};

  double sys_oo_y16to24[NooY16to24];
  for (int i = 0; i < NooY16to24; ++i)
    sys_oo_y16to24[i] = 0.0; // sys=0

  double x_oo_y16to24[NooY16to24], ex_oo_y16to24[NooY16to24];
  for (int i = 0; i < NooY16to24; ++i)
  {
    x_oo_y16to24[i] = 0.5 * (pt_low_oo_y16to24[i] + pt_high_oo_y16to24[i]);
    ex_oo_y16to24[i] = 0.5 * (pt_high_oo_y16to24[i] - pt_low_oo_y16to24[i]);
  }

  // ===== BPH pp (8 bins) =====
  // pp BPH results
  // from https://doi.org/10.17182/hepdata.57532.v1/t7
  const int Nbph = 8;

  double pt_low_bph[Nbph] = {0.0, 1.25, 2.0, 2.75, 3.5, 4.5, 6.5, 10.0};
  double pt_high_bph[Nbph] = {1.25, 2.0, 2.75, 3.5, 4.5, 6.5, 10.0, 30.0};

  double val_bph[Nbph] = {0.057, 0.087, 0.113, 0.139, 0.160, 0.177, 0.235, 0.374};
  double stat_bph[Nbph] = {0.021, 0.014, 0.013, 0.014, 0.014, 0.012, 0.016, 0.031};
  double sys_bph[Nbph] = {0.042, 0.022, 0.020, 0.010, 0.013, 0.012, 0.012, 0.008};

  double x_bph[Nbph], ex_bph[Nbph];
  for (int i = 0; i < Nbph; ++i)
  {
    x_bph[i] = 0.5 * (pt_low_bph[i] + pt_high_bph[i]);
    ex_bph[i] = 0.5 * (pt_high_bph[i] - pt_low_bph[i]);
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

  // BPH: sys band
  TGraphAsymmErrors *gSysBPH = new TGraphAsymmErrors(
      Nbph, x_bph, val_bph, ex_bph, ex_bph, sys_bph, sys_bph);
  gSysBPH->SetFillColor(kGray + 1);
  gSysBPH->SetFillStyle(3354);
  gSysBPH->SetLineColor(kBlack);
  gSysBPH->SetLineStyle(2);
  gSysBPH->SetLineWidth(2);

  // BPH: stat points
  TGraphErrors *gStatBPH = new TGraphErrors(Nbph, x_bph, val_bph, ex_bph, stat_bph);
  gStatBPH->SetMarkerStyle(25);
  gStatBPH->SetMarkerSize(1.5);
  gStatBPH->SetMarkerColor(kBlack);
  gStatBPH->SetLineColor(kBlack);
  gStatBPH->SetLineWidth(2);

  // ===== Draw =====
  TCanvas *c1 = new TCanvas("c1", "pp vs pp_old", 800, 800);
  c1->SetTitle("");
  TPad pad1("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
  pad1.SetTopMargin(0.08);
  pad1.SetBottomMargin(0.13);
  pad1.Draw();
  pad1.cd();
  gStyle->SetOptStat(0);
  gSyspp_old->SetTitle("");
  gSyspp_old->Draw("A2");
  gSyspp_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gSyspp_old->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
  gSyspp_old->GetXaxis()->CenterTitle();
  gSyspp_old->GetYaxis()->CenterTitle();
  gSyspp_old->GetXaxis()->SetTitleSize(0.05);
  gSyspp_old->GetYaxis()->SetTitleSize(0.05);
  gSyspp_old->GetXaxis()->SetLabelSize(0.042);
  gSyspp_old->GetYaxis()->SetLabelSize(0.042);
  gSyspp_old->GetYaxis()->SetTitleOffset(1.25);
  gSyspp_old->GetXaxis()->SetLimits(0, 30);
  gSyspp_old->SetMinimum(0);
  gSyspp_old->SetMaximum(0.9);

  gSyspp_old->Draw("A2");
  for (int i = 0; i < Nbph; ++i)
  {
    TBox *b = new TBox(
        x_bph[i] - ex_bph[i], val_bph[i] - sys_bph[i],
        x_bph[i] + ex_bph[i], val_bph[i] + sys_bph[i]);
    b->SetFillColorAlpha(kGray + 1, 0.25);
    b->SetFillStyle(1001);
    b->SetLineColor(kBlack);
    b->SetLineStyle(2);
    b->SetLineWidth(2);
    b->Draw("same");
  }

  gStatpp_old->Draw("P SAME");
  gStatBPH->Draw("P SAME");
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

  TLegend *leg = new TLegend(0.52, 0.70, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.024);
  leg->SetEntrySeparation(0.015);
  leg->AddEntry(gStatOOY16to24, "OO (2D fit)", "lp");
  leg->AddEntry(gStatpp_old, "OO (Ctau subrange)", "lp");
  leg->AddEntry(gStatBPH, "pp (BPH-10-002)", "lp");
  leg->Draw();

  c1->SaveAs("plot_OO_fwd.pdf");
  c1->SaveAs("plot_OO_fwd.png");
}
