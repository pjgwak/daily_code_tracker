#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"

void bFraction_comparison()
{
  const int N = 12;
  double pt_low[N] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0};
  double pt_high[N] = {7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};

  // pp_old - 2018 Raa (https://www.hepdata.net/record/ins1644903?version=1&table=Table%207)
  double val_pp_old[N] = {0.205, 0.228, 0.262, 0.295, 0.343,
      0.391, 0.442, 0.484, 0.515, 0.574,
      0.555, 0.597};
  double stat_pp_old[N] = {0.002, 0.002, 0.002, 0.002, 0.003,
      0.004, 0.005, 0.007, 0.007, 0.012,
      0.018, 0.021};
  double sys_pp_old[N] = {0.018, 0.015, 0.018, 0.015, 0.015, 0.017,
                        0.017, 0.020, 0.019, 0.027, 0.039, 0.034};

  // pp_new
  double val_pp[N] = {2.3286e-01, 0., 9, 0., 3.4129e-01,
      0., 5.0174e-01, 5.3474e-01, 0., 5.6880e-01,
      5.8908e-01, 6.1739e-01};
  double stat_pp[N] = {2.72e-02, 0.00, 0, 0.0, 8.49e-04,
      0.0, 9.28e-03, 4.79e-03, 0.00, 3.80e-03,
      7.38e-03, 6.66e-03};
  double sys_pp[N] = {0.0}; // systematic = 0

  double x[N], ex[N];
  for (int i = 0; i < N; i++)
  {
    x[i] = 0.5 * (pt_low[i] + pt_high[i]);
    ex[i] = 0.5 * (pt_high[i] - pt_low[i]);
  }

  // pp_old: sys band
  TGraphAsymmErrors *gSyspp_old = new TGraphAsymmErrors(N, x, val_pp_old, ex, ex, sys_pp_old, sys_pp_old);
  // gSyspp_old->SetFillColor(kRed - 7);
  gSyspp_old->SetFillColorAlpha(kRed - 7, 0.3);
  gSyspp_old->SetFillStyle(1001);
  // gSyspp_old->SetFillStyle(0); // no filling
  gSyspp_old->SetLineColor(kRed + 2);
  // gSyspp_old->SetLineWidth(2);

  // pp_old: stat error
  TGraphErrors *gStatpp_old = new TGraphErrors(N, x, val_pp_old, ex, stat_pp_old);
  gStatpp_old->SetMarkerStyle(24); // 21:square
  gStatpp_old->SetMarkerSize(1.2);
  gStatpp_old->SetMarkerColor(kRed + 2);
  gStatpp_old->SetLineColor(kRed + 2);

  // pp: sys=0 
  TGraphAsymmErrors *gSysPP = new TGraphAsymmErrors(N, x, val_pp, ex, ex, sys_pp, sys_pp);
  gSysPP->SetFillColor(0);
  gSysPP->SetLineColor(kBlue + 2);

  // pp: stat error
  TGraphErrors *gStatPP = new TGraphErrors(N, x, val_pp, ex, stat_pp);
  gStatPP->SetMarkerStyle(33); // 24:open circle, 33:diamond
  gStatPP->SetMarkerSize(1.6); // size for 1.2:circle, 1.6: diamond
  gStatPP->SetMarkerColor(kBlue + 2);
  gStatPP->SetLineColor(kBlue + 2);

  TCanvas *c1 = new TCanvas("c1", "pp vs pp_old", 800, 600);

  gSyspp_old->Draw("A2");
  gSyspp_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gSyspp_old->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
  gSyspp_old->GetXaxis()->SetLimits(0, 50);
  gSyspp_old->SetMinimum(0.1);
  gSyspp_old->SetMaximum(0.9);

  gStatpp_old->Draw("P SAME");
  gStatPP->Draw("P SAME");

  TLegend *leg = new TLegend(0.5, 0.2, 0.85, 0.35);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gStatPP, "pp_new", "p");
  leg->AddEntry(gStatpp_old, "pp_old", "p");
  leg->Draw();

  // c1->SaveAs("plot_pp_pp_old.pdf");
  c1->SaveAs("plot_pp_pp_old.png");
}
