void result_JPsi_XS_PPandPbPb_pt_nonprompt_2_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-10,-5.004845,52.5,3.368867);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.04);
   c1->SetTopMargin(0.08);
   c1->SetBottomMargin(0.12);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *haxes__75 = new TH1F("haxes__75","haxes",1,0,50);
   haxes__75->SetMinimum(0.0001);
   haxes__75->SetMaximum(500);
   haxes__75->SetLineStyle(0);
   haxes__75->SetMarkerStyle(20);
   haxes__75->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__75->GetXaxis()->CenterTitle(true);
   haxes__75->GetXaxis()->SetLabelFont(42);
   haxes__75->GetXaxis()->SetLabelOffset(0.007);
   haxes__75->GetXaxis()->SetLabelSize(0.05);
   haxes__75->GetXaxis()->SetTitleSize(0.06);
   haxes__75->GetXaxis()->SetTitleOffset(0.9);
   haxes__75->GetXaxis()->SetTitleFont(42);
   haxes__75->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/dp_{T} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/dp_{T} (nb / GeV/#it{c})");
   haxes__75->GetYaxis()->SetLabelFont(42);
   haxes__75->GetYaxis()->SetLabelOffset(0.007);
   haxes__75->GetYaxis()->SetLabelSize(0.05);
   haxes__75->GetYaxis()->SetTitleSize(0.043);
   haxes__75->GetYaxis()->SetTitleOffset(1.67);
   haxes__75->GetYaxis()->SetTitleFont(42);
   haxes__75->GetZaxis()->SetLabelFont(42);
   haxes__75->GetZaxis()->SetLabelOffset(0.007);
   haxes__75->GetZaxis()->SetLabelSize(0.05);
   haxes__75->GetZaxis()->SetTitleSize(0.06);
   haxes__75->GetZaxis()->SetTitleFont(42);
   haxes__75->Draw("");
   
   Double_t bin_0_syst_pp_fx3163[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_syst_pp_fy3163[12] = {
   2.596368,
   1.538445,
   1.003599,
   0.6027631,
   0.3107142,
   0.1584082,
   0.08201487,
   0.04160491,
   0.01725289,
   0.006657393,
   0.002556503,
   0.0006876037};
   Double_t bin_0_syst_pp_felx3163[12] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_pp_fely3163[12] = {
   0.167091,
   0.08360357,
   0.05058936,
   0.02403745,
   0.01057471,
   0.005512916,
   0.002893433,
   0.001575601,
   0.0006574182,
   0.0002736219,
   0.0001663665,
   3.950547e-05};
   Double_t bin_0_syst_pp_fehx3163[12] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_pp_fehy3163[12] = {
   0.167091,
   0.08360357,
   0.05058936,
   0.02403745,
   0.01057471,
   0.005512916,
   0.002893433,
   0.001575601,
   0.0006574182,
   0.0002736219,
   0.0001663665,
   3.950547e-05};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(12,bin_0_syst_pp_fx3163,bin_0_syst_pp_fy3163,bin_0_syst_pp_felx3163,bin_0_syst_pp_fehx3163,bin_0_syst_pp_fely3163,bin_0_syst_pp_fehy3163);
   grae->SetName("bin_0_syst_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_pp31613163 = new TH1F("Graph_Graph_bin_0_syst_pp31613163","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_pp31613163->SetMinimum(0.0005832884);
   Graph_Graph_bin_0_syst_pp31613163->SetMaximum(3.03974);
   Graph_Graph_bin_0_syst_pp31613163->SetDirectory(0);
   Graph_Graph_bin_0_syst_pp31613163->SetStats(0);
   Graph_Graph_bin_0_syst_pp31613163->SetLineStyle(0);
   Graph_Graph_bin_0_syst_pp31613163->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_pp31613163->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_pp31613163->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp31613163->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31613163->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31613163->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31613163->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31613163->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_pp31613163);
   
   grae->Draw("5");
   
   Double_t bin_0_pp_fx3164[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_pp_fy3164[12] = {
   2.596368,
   1.538445,
   1.003599,
   0.6027631,
   0.3107142,
   0.1584082,
   0.08201487,
   0.04160491,
   0.01725289,
   0.006657393,
   0.002556503,
   0.0006876037};
   Double_t bin_0_pp_felx3164[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_pp_fely3164[12] = {
   0.03098285,
   0.01604794,
   0.01051309,
   0.005556864,
   0.00297998,
   0.001952737,
   0.00116067,
   0.0007944869,
   0.0003446421,
   0.0002064181,
   0.0001370955,
   3.683226e-05};
   Double_t bin_0_pp_fehx3164[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_pp_fehy3164[12] = {
   0.03098285,
   0.01604794,
   0.01051309,
   0.005556864,
   0.00297998,
   0.001952737,
   0.00116067,
   0.0007944869,
   0.0003446421,
   0.0002064181,
   0.0001370955,
   3.683226e-05};
   grae = new TGraphAsymmErrors(12,bin_0_pp_fx3164,bin_0_pp_fy3164,bin_0_pp_felx3164,bin_0_pp_fehx3164,bin_0_pp_fely3164,bin_0_pp_fehy3164);
   grae->SetName("bin_0_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(4);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph_bin_0_pp31623164 = new TH1F("Graph_Graph_bin_0_pp31623164","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_pp31623164->SetMinimum(0.0005856943);
   Graph_Graph_bin_0_pp31623164->SetMaximum(2.890021);
   Graph_Graph_bin_0_pp31623164->SetDirectory(0);
   Graph_Graph_bin_0_pp31623164->SetStats(0);
   Graph_Graph_bin_0_pp31623164->SetLineStyle(0);
   Graph_Graph_bin_0_pp31623164->SetMarkerStyle(20);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_pp31623164->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_pp31623164->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp31623164->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31623164->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31623164->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31623164->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31623164->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_pp31623164);
   
   grae->Draw("p");
   
   Double_t bin_0_syst_aa_fx3165[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_syst_aa_fy3165[12] = {
   1.35157,
   0.7355143,
   0.4633089,
   0.2635224,
   0.1324728,
   0.07263043,
   0.03473078,
   0.01873368,
   0.007604387,
   0.00310657,
   0.001538558,
   0.0003258575};
   Double_t bin_0_syst_aa_felx3165[12] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_aa_fely3165[12] = {
   0.1210645,
   0.04730661,
   0.02515945,
   0.01286488,
   0.006144952,
   0.003249876,
   0.001548266,
   0.000972,
   0.0003334431,
   0.0001696841,
   0.0001062404,
   2.421701e-05};
   Double_t bin_0_syst_aa_fehx3165[12] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_aa_fehy3165[12] = {
   0.1210645,
   0.04730661,
   0.02515945,
   0.01286488,
   0.006144952,
   0.003249876,
   0.001548266,
   0.000972,
   0.0003334431,
   0.0001696841,
   0.0001062404,
   2.421701e-05};
   grae = new TGraphAsymmErrors(12,bin_0_syst_aa_fx3165,bin_0_syst_aa_fy3165,bin_0_syst_aa_felx3165,bin_0_syst_aa_fehx3165,bin_0_syst_aa_fely3165,bin_0_syst_aa_fehy3165);
   grae->SetName("bin_0_syst_aa");
   grae->SetTitle("Graph");

   ci = 2053;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_aa31593165 = new TH1F("Graph_Graph_bin_0_syst_aa31593165","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_aa31593165->SetMinimum(0.0002714764);
   Graph_Graph_bin_0_syst_aa31593165->SetMaximum(1.619868);
   Graph_Graph_bin_0_syst_aa31593165->SetDirectory(0);
   Graph_Graph_bin_0_syst_aa31593165->SetStats(0);
   Graph_Graph_bin_0_syst_aa31593165->SetLineStyle(0);
   Graph_Graph_bin_0_syst_aa31593165->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_aa31593165->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_aa31593165->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa31593165->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31593165->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31593165->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31593165->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31593165->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_aa31593165);
   
   grae->Draw("5");
   
   Double_t bin_0_aa_fx3166[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_aa_fy3166[12] = {
   1.35157,
   0.7355143,
   0.4633089,
   0.2635224,
   0.1324728,
   0.07263043,
   0.03473078,
   0.01873368,
   0.007604387,
   0.00310657,
   0.001538558,
   0.0003258575};
   Double_t bin_0_aa_felx3166[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_aa_fely3166[12] = {
   0.04334302,
   0.02008135,
   0.01195873,
   0.006027645,
   0.003126063,
   0.002060846,
   0.001176543,
   0.0008150887,
   0.0003519571,
   0.0002144256,
   0.0001493859,
   3.862374e-05};
   Double_t bin_0_aa_fehx3166[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_aa_fehy3166[12] = {
   0.04334302,
   0.02008135,
   0.01195873,
   0.006027645,
   0.003126063,
   0.002060846,
   0.001176543,
   0.0008150887,
   0.0003519571,
   0.0002144256,
   0.0001493859,
   3.862374e-05};
   grae = new TGraphAsymmErrors(12,bin_0_aa_fx3166,bin_0_aa_fy3166,bin_0_aa_felx3166,bin_0_aa_fehx3166,bin_0_aa_fely3166,bin_0_aa_fehy3166);
   grae->SetName("bin_0_aa");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(4);

   ci = TColor::GetColor("#990000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph_bin_0_aa31603166 = new TH1F("Graph_Graph_bin_0_aa31603166","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_aa31603166->SetMinimum(0.0002585104);
   Graph_Graph_bin_0_aa31603166->SetMaximum(1.534375);
   Graph_Graph_bin_0_aa31603166->SetDirectory(0);
   Graph_Graph_bin_0_aa31603166->SetStats(0);
   Graph_Graph_bin_0_aa31603166->SetLineStyle(0);
   Graph_Graph_bin_0_aa31603166->SetMarkerStyle(20);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_aa31603166->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_aa31603166->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa31603166->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31603166->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31603166->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31603166->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31603166->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_aa31603166);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.2,0.17,0.5,0.33,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.037);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("bin_0_pp","pp","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("bin_0_aa","PbPb, Cent. 0-100%","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#990000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.865,0.67,"pp global unc. 2.3 %");
tex->SetNDC();
   tex->SetTextAlign(32);
   tex->SetTextFont(42);
   tex->SetTextSize(0.034);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.92,0.6,"PbPb global unc. #splitline{+3.4 %}{#minus3.9 %}");
tex->SetNDC();
   tex->SetTextAlign(32);
   tex->SetTextFont(42);
   tex->SetTextSize(0.034);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.85,"J/#psi from b hadrons");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.055);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8,"J/#psi #rightarrow #mu^{+}#mu^{-}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.055);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.75,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.037);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.96,0.9424,"PbPb 368 #mub^{-1}, pp 28.0 pb^{-1} (5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9424,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.914,0.874,"CMS");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH1F *haxes_copy__76 = new TH1F("haxes_copy__76","haxes",1,0,50);
   haxes_copy__76->SetMinimum(0.0001);
   haxes_copy__76->SetMaximum(500);
   haxes_copy__76->SetDirectory(0);
   haxes_copy__76->SetLineStyle(0);
   haxes_copy__76->SetMarkerStyle(20);
   haxes_copy__76->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__76->GetXaxis()->CenterTitle(true);
   haxes_copy__76->GetXaxis()->SetLabelFont(42);
   haxes_copy__76->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__76->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__76->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__76->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__76->GetXaxis()->SetTitleFont(42);
   haxes_copy__76->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/dp_{T} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/dp_{T} (nb / GeV/#it{c})");
   haxes_copy__76->GetYaxis()->SetLabelFont(42);
   haxes_copy__76->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__76->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__76->GetYaxis()->SetTitleSize(0.043);
   haxes_copy__76->GetYaxis()->SetTitleOffset(1.67);
   haxes_copy__76->GetYaxis()->SetTitleFont(42);
   haxes_copy__76->GetZaxis()->SetLabelFont(42);
   haxes_copy__76->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__76->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__76->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__76->GetZaxis()->SetTitleFont(42);
   haxes_copy__76->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
