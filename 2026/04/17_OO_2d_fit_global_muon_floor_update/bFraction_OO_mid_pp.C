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

void bFraction_OO_mid_pp()
{
  TStopwatch __execTimer;
  __execTimer.Start();
  struct __ExecTimerGuard {
    TStopwatch *sw;
    explicit __ExecTimerGuard(TStopwatch *timer) : sw(timer) {}
    ~__ExecTimerGuard() {
      sw->Stop();
      sw->Print();
    }
  } __execTimerGuard(&__execTimer);

  // =====  CMS plot style =====
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // ===== pp_bFraction (12 bins) =====
  // OO, ctau subrange method (updated bins)
  const int n_old = 8;

  double pt_low_old[n_old] = {7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  double pt_high_old[n_old] = {8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 35.0};

  // f_B = 1 - f_prompt
  double val_pp_old[n_old] = {
      0.230,
      0.259,
      0.285,
      0.314,
      0.376,
      0.402,
      0.485,
      0.546 
  };

  // statistical errors: same as prompt fraction
  double stat_pp_old[n_old] = {
      0.008,
      0.009,
      0.010,
      0.009,
      0.013,
      0.019,
      0.019,
      0.024};

  // systematics not provided
  double sys_pp_old[n_old];
  for (int i = 0; i < n_old; ++i)
    sys_pp_old[i] = 0.0; // sys=0

  double x_old[n_old], ex_old[n_old];
  for (int i = 0; i < n_old; ++i)
  {
    x_old[i] = 0.5 * (pt_low_old[i] + pt_high_old[i]);
    ex_old[i] = 0.5 * (pt_high_old[i] - pt_low_old[i]);
  }

  // ===== OO |y|<1.6 (2d fit) =====
  // y = 0.00-1.60
  const int NooY16 = 8;
  double pt_low_oo_y16[NooY16] = {7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0};
  double pt_high_oo_y16[NooY16] = {8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 35.0};

  double val_oo_y16[NooY16] = {0.23147, 0.27556, 0.28637, 0.32747, 0.39717, 0.40085, 0.48678, 0.54760};
  double stat_oo_y16[NooY16] = {0.00867, 0.00939, 0.01083, 0.00985, 0.01361, 0.01854, 0.01916, 0.02389};

  double sys_oo_y16[NooY16];
  for (int i = 0; i < NooY16; ++i)
    sys_oo_y16[i] = 0.0; // sys=0

  double x_oo_y16[NooY16], ex_oo_y16[NooY16];
  for (int i = 0; i < NooY16; ++i)
  {
    x_oo_y16[i] = 0.5 * (pt_low_oo_y16[i] + pt_high_oo_y16[i]);
    ex_oo_y16[i] = 0.5 * (pt_high_oo_y16[i] - pt_low_oo_y16[i]);
  }

  // ===== BPH pp (8 bins) =====
  // pp BPH results
  // from https://doi.org/10.17182/hepdata.57532.v1/t7 
  const int Nbph = 8;

  double pt_low_bph[Nbph]  = {0.0, 1.25, 2.0, 2.75, 3.5, 4.5, 6.5, 10.0};
  double pt_high_bph[Nbph] = {1.25, 2.0, 2.75, 3.5, 4.5, 6.5, 10.0, 30.0};

  double val_bph[Nbph]  = {0.057, 0.087, 0.113, 0.139, 0.160, 0.177, 0.235, 0.374};
  double stat_bph[Nbph] = {0.021, 0.014, 0.013, 0.014, 0.014, 0.012, 0.016, 0.031};
  double sys_bph[Nbph]  = {0.042, 0.022, 0.020, 0.010, 0.013, 0.012, 0.012, 0.008};

  double x_bph[Nbph], ex_bph[Nbph];
  for (int i = 0; i < Nbph; ++i)
  {
    x_bph[i]  = 0.5 * (pt_low_bph[i] + pt_high_bph[i]);
    ex_bph[i] = 0.5 * (pt_high_bph[i] - pt_low_bph[i]);
  }

  // ===== pp (|y| < 1.2, 2 bins; bFraction input) =====
  const int NppMid = 2;
  double pt_low_pp_mid[NppMid] = {6.5, 10.0};
  double pt_high_pp_mid[NppMid] = {10.0, 30.0};
  double val_pp_mid[NppMid] = {0.257, 0.395};
  double stat_pp_mid[NppMid] = {0.015, 0.018};
  double sys_pp_mid[NppMid] = {0.014, 0.005};


  double x_pp_mid[NppMid], ex_pp_mid[NppMid];
  for (int i = 0; i < NppMid; ++i)
  {
    x_pp_mid[i] = 0.5 * (pt_low_pp_mid[i] + pt_high_pp_mid[i]);
    ex_pp_mid[i] = 0.5 * (pt_high_pp_mid[i] - pt_low_pp_mid[i]);
  }

  // ===== pp (1.2 < |y| < 1.6, 4 bins; bFraction input) =====
  const int NppY12to16 = 4;
  double pt_low_pp_y12to16[NppY12to16] = {2.0, 4.5, 6.5, 10.0};
  double pt_high_pp_y12to16[NppY12to16] = {4.5, 6.5, 10.0, 30.0};
  double val_pp_y12to16[NppY12to16] = {0.146, 0.180, 0.203, 0.360};
  double stat_pp_y12to16[NppY12to16] = {0.021, 0.017, 0.017, 0.031};
  double sys_pp_y12to16[NppY12to16] = {0.028, 0.019, 0.014, 0.016};

  double x_pp_y12to16[NppY12to16], ex_pp_y12to16[NppY12to16];
  for (int i = 0; i < NppY12to16; ++i)
  {
    x_pp_y12to16[i] = 0.5 * (pt_low_pp_y12to16[i] + pt_high_pp_y12to16[i]);
    ex_pp_y12to16[i] = 0.5 * (pt_high_pp_y12to16[i] - pt_low_pp_y12to16[i]);
  }


  // ===== pp (|y| < 2.4, CMS 5.02 TeV) =====
  // https://www.hepdata.net/record/83746
  const int NppY24 = 12;

  double pt_low_pp_y24[NppY24] = {
      6.5, 7.5, 8.5, 9.5,
      11.0, 13.0, 15.0, 17.5,
      20.0, 25.0, 30.0, 35.0};

  double pt_high_pp_y24[NppY24] = {
      7.5, 8.5, 9.5, 11.0,
      13.0, 15.0, 17.5, 20.0,
      25.0, 30.0, 35.0, 50.0};

  double val_pp_y24[NppY24] = { 0.205, 0.228, 0.262, 0.295, 0.343, 0.391, 0.442, 0.484, 0.515, 0.574, 0.555, 0.597};

  double stat_pp_y24[NppY24] = { 0.002, 0.002, 0.002, 0.002, 0.003, 0.004, 0.005, 0.007, 0.007, 0.012, 0.018, 0.021};

  double sys_pp_y24[NppY24] = { 0.018, 0.015, 0.018, 0.015, 0.015, 0.017, 0.017, 0.020, 0.019, 0.027, 0.039, 0.034};

  double x_pp_y24[NppY24], ex_pp_y24[NppY24];
  for (int i = 0; i < NppY24; ++i)
  {
    x_pp_y24[i] = 0.5 * (pt_low_pp_y24[i] + pt_high_pp_y24[i]);
    ex_pp_y24[i] = 0.5 * (pt_high_pp_y24[i] - pt_low_pp_y24[i]);
  }

  // ===== PbPb (|y| < 2.4, CMS 5.02 TeV) =====
  // https://www.hepdata.net/record/83746
  const int NpbpbY24 = 12;

  double pt_low_pbpb_y24[NpbpbY24] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0};
  double pt_high_pbpb_y24[NpbpbY24] = {7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};
  double val_pbpb_y24[NpbpbY24] = {0.272, 0.288, 0.329, 0.356, 0.396, 0.450, 0.483, 0.537, 0.525, 0.547, 0.602, 0.574};
  double stat_pbpb_y24[NpbpbY24] = {0.008, 0.007, 0.007, 0.007, 0.008, 0.010, 0.012, 0.017, 0.018, 0.026, 0.038, 0.046};
  double sys_pbpb_y24[NpbpbY24] = {0.032, 0.021, 0.019, 0.019, 0.020, 0.020, 0.020, 0.030, 0.020, 0.028, 0.052, 0.039};

  double x_pbpb_y24[NpbpbY24], ex_pbpb_y24[NpbpbY24];
  for (int i = 0; i < NpbpbY24; ++i)
  {
    x_pbpb_y24[i] = 0.5 * (pt_low_pbpb_y24[i] + pt_high_pbpb_y24[i]);
    ex_pbpb_y24[i] = 0.5 * (pt_high_pbpb_y24[i] - pt_low_pbpb_y24[i]);
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
  gStatpp_old->SetMarkerSize(1.7);
  gStatpp_old->SetMarkerColor(kRed + 1);
  gStatpp_old->SetLineColor(kRed + 1);
  gStatpp_old->SetLineWidth(3);

  // OO |y|<1.6: sys (currently 0)
  TGraphAsymmErrors *gSysOOY16 = new TGraphAsymmErrors(NooY16, x_oo_y16, val_oo_y16, ex_oo_y16, ex_oo_y16, sys_oo_y16, sys_oo_y16);
  gSysOOY16->SetFillColorAlpha(kAzure + 2, 0.20);
  gSysOOY16->SetFillStyle(1001);
  gSysOOY16->SetLineColor(kAzure + 2);
  gSysOOY16->SetLineWidth(3);

  // OO |y|<1.6: stat
  TGraphErrors *gStatOOY16 = new TGraphErrors(NooY16, x_oo_y16, val_oo_y16, ex_oo_y16, stat_oo_y16);
  gStatOOY16->SetMarkerStyle(20);
  gStatOOY16->SetMarkerSize(1.7);
  gStatOOY16->SetMarkerColor(kAzure + 2);
  gStatOOY16->SetLineColor(kAzure + 2);
  gStatOOY16->SetLineWidth(3);

  // OO |y|<1.6 with PR MC: sys (currently 0)

  // OO |y|<1.6 with PR MC: stat

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

  // pp (|y| < 1.2): stat points
  TGraphAsymmErrors *gSysPPMid = new TGraphAsymmErrors(NppMid, x_pp_mid, val_pp_mid, ex_pp_mid, ex_pp_mid, sys_pp_mid, sys_pp_mid);
  gSysPPMid->SetFillColorAlpha(kGreen + 2, 0.20);
  gSysPPMid->SetFillStyle(1001);
  gSysPPMid->SetLineColor(kGreen + 2);
  gSysPPMid->SetLineWidth(2);

  // pp (|y| < 1.2): stat points
  TGraphErrors *gStatPPMid = new TGraphErrors(NppMid, x_pp_mid, val_pp_mid, ex_pp_mid, stat_pp_mid);
  gStatPPMid->SetMarkerStyle(33);
  gStatPPMid->SetMarkerSize(1.7);
  gStatPPMid->SetMarkerColor(kGreen + 2);
  gStatPPMid->SetLineColor(kGreen + 2);
  gStatPPMid->SetLineWidth(1);

  // pp (1.2 < |y| < 1.6): stat points
  TGraphAsymmErrors *gSysPPY12to16 = new TGraphAsymmErrors(NppY12to16, x_pp_y12to16, val_pp_y12to16, ex_pp_y12to16, ex_pp_y12to16, sys_pp_y12to16, sys_pp_y12to16);
  gSysPPY12to16->SetFillColorAlpha(kOrange + 7, 0.20);
  gSysPPY12to16->SetFillStyle(1001);
  gSysPPY12to16->SetLineColor(kOrange + 7);
  gSysPPY12to16->SetLineWidth(2);

  // pp (1.2 < |y| < 1.6): stat points
  TGraphErrors *gStatPPY12to16 = new TGraphErrors(NppY12to16, x_pp_y12to16, val_pp_y12to16, ex_pp_y12to16, stat_pp_y12to16);
  gStatPPY12to16->SetMarkerStyle(34);
  gStatPPY12to16->SetMarkerSize(1.7);
  gStatPPY12to16->SetMarkerColor(kOrange + 7);
  gStatPPY12to16->SetLineColor(kOrange + 7);
  gStatPPY12to16->SetLineWidth(1);

  // pp (|y| < 2.4): sys/stat points
  TGraphAsymmErrors *gSysPPY24 = new TGraphAsymmErrors(NppY24, x_pp_y24, val_pp_y24, ex_pp_y24, ex_pp_y24, sys_pp_y24, sys_pp_y24);
  gSysPPY24->SetFillColorAlpha(kViolet + 1, 0.16);
  gSysPPY24->SetFillStyle(1001);
  gSysPPY24->SetLineColor(kViolet + 1);
  gSysPPY24->SetLineWidth(2);

  TGraphErrors *gStatPPY24 = new TGraphErrors(NppY24, x_pp_y24, val_pp_y24, ex_pp_y24, stat_pp_y24);
  gStatPPY24->SetMarkerStyle(22);
  gStatPPY24->SetMarkerSize(1.1);
  gStatPPY24->SetMarkerColor(kViolet + 1);
  gStatPPY24->SetLineColor(kViolet + 1);
  gStatPPY24->SetLineWidth(1);

  // PbPb (|y| < 2.4): sys/stat points
  TGraphAsymmErrors *gSysPbPbY24 = new TGraphAsymmErrors(NpbpbY24, x_pbpb_y24, val_pbpb_y24, ex_pbpb_y24, ex_pbpb_y24, sys_pbpb_y24, sys_pbpb_y24);
  gSysPbPbY24->SetFillColorAlpha(kTeal + 2, 0.16);
  gSysPbPbY24->SetFillStyle(1001);
  gSysPbPbY24->SetLineColor(kTeal + 2);
  gSysPbPbY24->SetLineWidth(2);

  TGraphErrors *gStatPbPbY24 = new TGraphErrors(NpbpbY24, x_pbpb_y24, val_pbpb_y24, ex_pbpb_y24, stat_pbpb_y24);
  gStatPbPbY24->SetMarkerStyle(23);
  gStatPbPbY24->SetMarkerSize(1.1);
  gStatPbPbY24->SetMarkerColor(kTeal + 2);
  gStatPbPbY24->SetLineColor(kTeal + 2);
  gStatPbPbY24->SetLineWidth(1);

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
  gSyspp_old->GetXaxis()->SetLimits(0, 50);
  gSyspp_old->SetMinimum(0);
  gSyspp_old->SetMaximum(1.2);

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

  gStatBPH->Draw("P SAME");
  gSysPPY24->Draw("2 SAME");
  gSysPbPbY24->Draw("2 SAME");
  gSysPPMid->Draw("2 SAME");
  gSysPPY12to16->Draw("2 SAME");
  gStatPPY24->Draw("P SAME");
  gStatPbPbY24->Draw("P SAME");
  gStatPPMid->Draw("P SAME");
  gStatPPY12to16->Draw("P SAME");
  gSysOOY16->Draw("2 SAME");
  gStatpp_old->Draw("P SAME");
  gStatOOY16->Draw("P SAME");

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
    tx.DrawLatex(0.19, 0.865, "|y| < 1.6");
  }

  TLegend *leg = new TLegend(0.46, 0.60, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.017);
  leg->SetEntrySeparation(0.01);
  leg->AddEntry(gStatOOY16, "OO (|y| < 1.6, 2D fit)", "lp");
  leg->AddEntry(gStatpp_old, "OO (|y| < 1.6, Ctau subrange)", "lp");
  leg->AddEntry(gStatPPY24, "pp (CMS 5.02 TeV, |y| < 2.4)", "lp");
  leg->AddEntry(gStatPbPbY24, "PbPb (CMS 5.02 TeV, |y| < 2.4)", "lp");
  leg->AddEntry(gStatPPMid, "pp (CMS 7 TeV, |y| < 1.2)", "lp");
  leg->AddEntry(gStatPPY12to16, "pp (CMS 7 TeV, 1.2 < |y| < 1.6)", "lp");
  leg->AddEntry(gStatBPH, "pp (CMS 7 TeV, 1.6 < |y| < 2.4)", "lp");
  leg->Draw();

  c1->SaveAs("plot_OO_mid.pdf");
  c1->SaveAs("plot_OO_mid.png");
}
