void result_JPsi_BF_PPandPbPb_pt_2_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:56 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-10,-0.135,52.5,0.99);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
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
   
   TH1F *haxes__109 = new TH1F("haxes__109","haxes",1,0,50);
   haxes__109->SetMinimum(0);
   haxes__109->SetMaximum(0.9);
   haxes__109->SetLineStyle(0);
   haxes__109->SetMarkerStyle(20);
   haxes__109->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__109->GetXaxis()->CenterTitle(true);
   haxes__109->GetXaxis()->SetLabelFont(42);
   haxes__109->GetXaxis()->SetLabelOffset(0.007);
   haxes__109->GetXaxis()->SetLabelSize(0.05);
   haxes__109->GetXaxis()->SetTitleSize(0.06);
   haxes__109->GetXaxis()->SetTitleOffset(0.9);
   haxes__109->GetXaxis()->SetTitleFont(42);
   haxes__109->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
   haxes__109->GetYaxis()->SetLabelFont(42);
   haxes__109->GetYaxis()->SetLabelOffset(0.007);
   haxes__109->GetYaxis()->SetLabelSize(0.05);
   haxes__109->GetYaxis()->SetTitleSize(0.06);
   haxes__109->GetYaxis()->SetTitleOffset(1.25);
   haxes__109->GetYaxis()->SetTitleFont(42);
   haxes__109->GetZaxis()->SetLabelFont(42);
   haxes__109->GetZaxis()->SetLabelOffset(0.007);
   haxes__109->GetZaxis()->SetLabelSize(0.05);
   haxes__109->GetZaxis()->SetTitleSize(0.06);
   haxes__109->GetZaxis()->SetTitleFont(42);
   haxes__109->Draw("");
   
   Double_t bin_0_syst_pp_fx3229[12] = {
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
   Double_t bin_0_syst_pp_fy3229[12] = {
   0.2054889,
   0.2275798,
   0.2617972,
   0.2946761,
   0.3429313,
   0.3907621,
   0.442085,
   0.4837158,
   0.5154892,
   0.5736369,
   0.5546741,
   0.5974332};
   Double_t bin_0_syst_pp_felx3229[12] = {
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
   Double_t bin_0_syst_pp_fely3229[12] = {
   0.01779859,
   0.0154695,
   0.01763396,
   0.01529129,
   0.01481817,
   0.01658948,
   0.01729382,
   0.01997755,
   0.01862013,
   0.02722514,
   0.03865391,
   0.03373631};
   Double_t bin_0_syst_pp_fehx3229[12] = {
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
   Double_t bin_0_syst_pp_fehy3229[12] = {
   0.01779859,
   0.0154695,
   0.01763396,
   0.01529129,
   0.01481817,
   0.01658948,
   0.01729382,
   0.01997755,
   0.01862013,
   0.02722514,
   0.03865391,
   0.03373631};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(12,bin_0_syst_pp_fx3229,bin_0_syst_pp_fy3229,bin_0_syst_pp_felx3229,bin_0_syst_pp_fehx3229,bin_0_syst_pp_fely3229,bin_0_syst_pp_fehy3229);
   grae->SetName("bin_0_syst_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_pp32273229 = new TH1F("Graph_Graph_bin_0_syst_pp32273229","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_pp32273229->SetMinimum(0.1433424);
   Graph_Graph_bin_0_syst_pp32273229->SetMaximum(0.6755174);
   Graph_Graph_bin_0_syst_pp32273229->SetDirectory(0);
   Graph_Graph_bin_0_syst_pp32273229->SetStats(0);
   Graph_Graph_bin_0_syst_pp32273229->SetLineStyle(0);
   Graph_Graph_bin_0_syst_pp32273229->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_pp32273229->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_pp32273229->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp32273229->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32273229->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32273229->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32273229->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32273229->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_pp32273229);
   
   grae->Draw("5");
   
   Double_t bin_0_pp_fx3230[12] = {
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
   Double_t bin_0_pp_fy3230[12] = {
   0.2054889,
   0.2275798,
   0.2617972,
   0.2946761,
   0.3429313,
   0.3907621,
   0.442085,
   0.4837158,
   0.5154892,
   0.5736369,
   0.5546741,
   0.5974332};
   Double_t bin_0_pp_felx3230[12] = {
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
   Double_t bin_0_pp_fely3230[12] = {
   0.002241082,
   0.002127785,
   0.00237449,
   0.00233039,
   0.002728743,
   0.003858328,
   0.004786639,
   0.006841942,
   0.00729952,
   0.0118172,
   0.01768031,
   0.02079184};
   Double_t bin_0_pp_fehx3230[12] = {
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
   Double_t bin_0_pp_fehy3230[12] = {
   0.002241082,
   0.002127785,
   0.00237449,
   0.00233039,
   0.002728743,
   0.003858328,
   0.004786639,
   0.006841942,
   0.00729952,
   0.0118172,
   0.01768031,
   0.02079184};
   grae = new TGraphAsymmErrors(12,bin_0_pp_fx3230,bin_0_pp_fy3230,bin_0_pp_felx3230,bin_0_pp_fehx3230,bin_0_pp_fely3230,bin_0_pp_fehy3230);
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
   
   TH1F *Graph_Graph_bin_0_pp32283230 = new TH1F("Graph_Graph_bin_0_pp32283230","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_pp32283230->SetMinimum(0.1617501);
   Graph_Graph_bin_0_pp32283230->SetMaximum(0.6597227);
   Graph_Graph_bin_0_pp32283230->SetDirectory(0);
   Graph_Graph_bin_0_pp32283230->SetStats(0);
   Graph_Graph_bin_0_pp32283230->SetLineStyle(0);
   Graph_Graph_bin_0_pp32283230->SetMarkerStyle(20);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_pp32283230->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_pp32283230->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp32283230->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32283230->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32283230->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32283230->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32283230->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_pp32283230);
   
   grae->Draw("p");
   
   Double_t bin_0_syst_aa_fx3231[12] = {
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
   Double_t bin_0_syst_aa_fy3231[12] = {
   0.2716523,
   0.2880834,
   0.3290758,
   0.3560706,
   0.3959586,
   0.4495076,
   0.4827531,
   0.5368658,
   0.5248507,
   0.5467119,
   0.6021392,
   0.574469};
   Double_t bin_0_syst_aa_felx3231[12] = {
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
   Double_t bin_0_syst_aa_fely3231[12] = {
   0.03205797,
   0.02069105,
   0.01895576,
   0.01882433,
   0.02023045,
   0.01992033,
   0.01987494,
   0.02961108,
   0.01971684,
   0.0280248,
   0.05236426,
   0.03854471};
   Double_t bin_0_syst_aa_fehx3231[12] = {
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
   Double_t bin_0_syst_aa_fehy3231[12] = {
   0.03205797,
   0.02069105,
   0.01895576,
   0.01882433,
   0.02023045,
   0.01992033,
   0.01987494,
   0.02961108,
   0.01971684,
   0.0280248,
   0.05236426,
   0.03854471};
   grae = new TGraphAsymmErrors(12,bin_0_syst_aa_fx3231,bin_0_syst_aa_fy3231,bin_0_syst_aa_felx3231,bin_0_syst_aa_fehx3231,bin_0_syst_aa_fely3231,bin_0_syst_aa_fehy3231);
   grae->SetName("bin_0_syst_aa");
   grae->SetTitle("Graph");

   ci = 2085;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_aa32253231 = new TH1F("Graph_Graph_bin_0_syst_aa32253231","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_aa32253231->SetMinimum(0.1981034);
   Graph_Graph_bin_0_syst_aa32253231->SetMaximum(0.6959944);
   Graph_Graph_bin_0_syst_aa32253231->SetDirectory(0);
   Graph_Graph_bin_0_syst_aa32253231->SetStats(0);
   Graph_Graph_bin_0_syst_aa32253231->SetLineStyle(0);
   Graph_Graph_bin_0_syst_aa32253231->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_aa32253231->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_aa32253231->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa32253231->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32253231->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32253231->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32253231->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32253231->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_aa32253231);
   
   grae->Draw("5");
   
   Double_t bin_0_aa_fx3232[12] = {
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
   Double_t bin_0_aa_fy3232[12] = {
   0.2716523,
   0.2880834,
   0.3290758,
   0.3560706,
   0.3959586,
   0.4495076,
   0.4827531,
   0.5368658,
   0.5248507,
   0.5467119,
   0.6021392,
   0.574469};
   Double_t bin_0_aa_felx3232[12] = {
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
   Double_t bin_0_aa_fely3232[12] = {
   0.007609592,
   0.006753502,
   0.007196596,
   0.00678885,
   0.007592882,
   0.009960738,
   0.01220233,
   0.01669277,
   0.01774508,
   0.02636693,
   0.03825734,
   0.04636884};
   Double_t bin_0_aa_fehx3232[12] = {
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
   Double_t bin_0_aa_fehy3232[12] = {
   0.007609592,
   0.006753502,
   0.007196596,
   0.00678885,
   0.007592882,
   0.009960738,
   0.01220233,
   0.01669277,
   0.01774508,
   0.02636693,
   0.03825734,
   0.04636884};
   grae = new TGraphAsymmErrors(12,bin_0_aa_fx3232,bin_0_aa_fy3232,bin_0_aa_felx3232,bin_0_aa_fehx3232,bin_0_aa_fely3232,bin_0_aa_fehy3232);
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
   
   TH1F *Graph_Graph_bin_0_aa32263232 = new TH1F("Graph_Graph_bin_0_aa32263232","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_aa32263232->SetMinimum(0.2264073);
   Graph_Graph_bin_0_aa32263232->SetMaximum(0.678032);
   Graph_Graph_bin_0_aa32263232->SetDirectory(0);
   Graph_Graph_bin_0_aa32263232->SetStats(0);
   Graph_Graph_bin_0_aa32263232->SetLineStyle(0);
   Graph_Graph_bin_0_aa32263232->SetMarkerStyle(20);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_aa32263232->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_aa32263232->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa32263232->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32263232->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32263232->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32263232->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32263232->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_aa32263232);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.185,0.67,0.5,0.83,NULL,"brNDC");
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
   TLatex *   tex = new TLatex(0.2,0.86,"|y| < 2.4");
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
   
   TH1F *haxes_copy__110 = new TH1F("haxes_copy__110","haxes",1,0,50);
   haxes_copy__110->SetMinimum(0);
   haxes_copy__110->SetMaximum(0.9);
   haxes_copy__110->SetDirectory(0);
   haxes_copy__110->SetLineStyle(0);
   haxes_copy__110->SetMarkerStyle(20);
   haxes_copy__110->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__110->GetXaxis()->CenterTitle(true);
   haxes_copy__110->GetXaxis()->SetLabelFont(42);
   haxes_copy__110->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__110->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__110->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__110->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__110->GetXaxis()->SetTitleFont(42);
   haxes_copy__110->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
   haxes_copy__110->GetYaxis()->SetLabelFont(42);
   haxes_copy__110->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__110->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__110->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__110->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__110->GetYaxis()->SetTitleFont(42);
   haxes_copy__110->GetZaxis()->SetLabelFont(42);
   haxes_copy__110->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__110->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__110->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__110->GetZaxis()->SetTitleFont(42);
   haxes_copy__110->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
