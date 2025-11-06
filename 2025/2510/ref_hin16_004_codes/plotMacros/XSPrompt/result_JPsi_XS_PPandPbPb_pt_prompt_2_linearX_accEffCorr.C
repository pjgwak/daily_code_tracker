void result_JPsi_XS_PPandPbPb_pt_prompt_2_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:51 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__27 = new TH1F("haxes__27","haxes",1,0,50);
   haxes__27->SetMinimum(0.0001);
   haxes__27->SetMaximum(500);
   haxes__27->SetLineStyle(0);
   haxes__27->SetMarkerStyle(20);
   haxes__27->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__27->GetXaxis()->CenterTitle(true);
   haxes__27->GetXaxis()->SetLabelFont(42);
   haxes__27->GetXaxis()->SetLabelOffset(0.007);
   haxes__27->GetXaxis()->SetLabelSize(0.05);
   haxes__27->GetXaxis()->SetTitleSize(0.06);
   haxes__27->GetXaxis()->SetTitleOffset(0.9);
   haxes__27->GetXaxis()->SetTitleFont(42);
   haxes__27->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/dp_{T} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/dp_{T} (nb / GeV/#it{c})");
   haxes__27->GetYaxis()->SetLabelFont(42);
   haxes__27->GetYaxis()->SetLabelOffset(0.007);
   haxes__27->GetYaxis()->SetLabelSize(0.05);
   haxes__27->GetYaxis()->SetTitleSize(0.043);
   haxes__27->GetYaxis()->SetTitleOffset(1.67);
   haxes__27->GetYaxis()->SetTitleFont(42);
   haxes__27->GetZaxis()->SetLabelFont(42);
   haxes__27->GetZaxis()->SetLabelOffset(0.007);
   haxes__27->GetZaxis()->SetLabelSize(0.05);
   haxes__27->GetZaxis()->SetTitleSize(0.06);
   haxes__27->GetZaxis()->SetTitleFont(42);
   haxes__27->Draw("");
   
   Double_t bin_0_syst_pp_fx3065[12] = {
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
   Double_t bin_0_syst_pp_fy3065[12] = {
   10.03871,
   5.22158,
   2.829898,
   1.442747,
   0.5953395,
   0.2469745,
   0.1035035,
   0.04440615,
   0.01621608,
   0.004948193,
   0.002052515,
   0.0004633262};
   Double_t bin_0_syst_pp_felx3065[12] = {
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
   Double_t bin_0_syst_pp_fely3065[12] = {
   0.5932276,
   0.2638729,
   0.1243144,
   0.05606081,
   0.02175842,
   0.009274044,
   0.004011535,
   0.001889757,
   0.0007026316,
   0.0003268274,
   0.0002078395,
   3.797501e-05};
   Double_t bin_0_syst_pp_fehx3065[12] = {
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
   Double_t bin_0_syst_pp_fehy3065[12] = {
   0.5932276,
   0.2638729,
   0.1243144,
   0.05606081,
   0.02175842,
   0.009274044,
   0.004011535,
   0.001889757,
   0.0007026316,
   0.0003268274,
   0.0002078395,
   3.797501e-05};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(12,bin_0_syst_pp_fx3065,bin_0_syst_pp_fy3065,bin_0_syst_pp_felx3065,bin_0_syst_pp_fehx3065,bin_0_syst_pp_fely3065,bin_0_syst_pp_fehy3065);
   grae->SetName("bin_0_syst_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_pp30633065 = new TH1F("Graph_Graph_bin_0_syst_pp30633065","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_pp30633065->SetMinimum(0.0003828161);
   Graph_Graph_bin_0_syst_pp30633065->SetMaximum(11.69509);
   Graph_Graph_bin_0_syst_pp30633065->SetDirectory(0);
   Graph_Graph_bin_0_syst_pp30633065->SetStats(0);
   Graph_Graph_bin_0_syst_pp30633065->SetLineStyle(0);
   Graph_Graph_bin_0_syst_pp30633065->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_pp30633065->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_pp30633065->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp30633065->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp30633065->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp30633065->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp30633065->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp30633065->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_pp30633065);
   
   grae->Draw("5");
   
   Double_t bin_0_pp_fx3066[12] = {
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
   Double_t bin_0_pp_fy3066[12] = {
   10.03871,
   5.22158,
   2.829898,
   1.442747,
   0.5953395,
   0.2469745,
   0.1035035,
   0.04440615,
   0.01621608,
   0.004948193,
   0.002052515,
   0.0004633262};
   Double_t bin_0_pp_felx3066[12] = {
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
   Double_t bin_0_pp_fely3066[12] = {
   0.1171035,
   0.05319822,
   0.0288862,
   0.01296806,
   0.005564225,
   0.002964048,
   0.001425661,
   0.0008283631,
   0.0003155892,
   0.0001502629,
   0.0001081051,
   2.46316e-05};
   Double_t bin_0_pp_fehx3066[12] = {
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
   Double_t bin_0_pp_fehy3066[12] = {
   0.1171035,
   0.05319822,
   0.0288862,
   0.01296806,
   0.005564225,
   0.002964048,
   0.001425661,
   0.0008283631,
   0.0003155892,
   0.0001502629,
   0.0001081051,
   2.46316e-05};
   grae = new TGraphAsymmErrors(12,bin_0_pp_fx3066,bin_0_pp_fy3066,bin_0_pp_felx3066,bin_0_pp_fehx3066,bin_0_pp_fely3066,bin_0_pp_fehy3066);
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
   
   TH1F *Graph_Graph_bin_0_pp30643066 = new TH1F("Graph_Graph_bin_0_pp30643066","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_pp30643066->SetMinimum(0.0003948251);
   Graph_Graph_bin_0_pp30643066->SetMaximum(11.17135);
   Graph_Graph_bin_0_pp30643066->SetDirectory(0);
   Graph_Graph_bin_0_pp30643066->SetStats(0);
   Graph_Graph_bin_0_pp30643066->SetLineStyle(0);
   Graph_Graph_bin_0_pp30643066->SetMarkerStyle(20);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_pp30643066->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_pp30643066->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp30643066->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp30643066->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp30643066->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp30643066->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp30643066->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_pp30643066);
   
   grae->Draw("p");
   
   Double_t bin_0_syst_aa_fx3067[12] = {
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
   Double_t bin_0_syst_aa_fy3067[12] = {
   3.623797,
   1.817616,
   0.9446006,
   0.4765622,
   0.2020894,
   0.08894734,
   0.03721237,
   0.01616085,
   0.00688428,
   0.00257571,
   0.001016595,
   0.000241375};
   Double_t bin_0_syst_aa_felx3067[12] = {
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
   Double_t bin_0_syst_aa_fely3067[12] = {
   0.2807501,
   0.1113268,
   0.04919691,
   0.02226501,
   0.009035507,
   0.004031356,
   0.001725644,
   0.000948767,
   0.0003431601,
   0.0001809504,
   0.0001269035,
   2.305294e-05};
   Double_t bin_0_syst_aa_fehx3067[12] = {
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
   Double_t bin_0_syst_aa_fehy3067[12] = {
   0.2807501,
   0.1113268,
   0.04919691,
   0.02226501,
   0.009035507,
   0.004031356,
   0.001725644,
   0.000948767,
   0.0003431601,
   0.0001809504,
   0.0001269035,
   2.305294e-05};
   grae = new TGraphAsymmErrors(12,bin_0_syst_aa_fx3067,bin_0_syst_aa_fy3067,bin_0_syst_aa_felx3067,bin_0_syst_aa_fehx3067,bin_0_syst_aa_fely3067,bin_0_syst_aa_fehy3067);
   grae->SetName("bin_0_syst_aa");
   grae->SetTitle("Graph");

   ci = 1991;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_aa30613067 = new TH1F("Graph_Graph_bin_0_syst_aa30613067","Graph",100,2.85,46.65);
   Graph_Graph_bin_0_syst_aa30613067->SetMinimum(0.0001964899);
   Graph_Graph_bin_0_syst_aa30613067->SetMaximum(4.29498);
   Graph_Graph_bin_0_syst_aa30613067->SetDirectory(0);
   Graph_Graph_bin_0_syst_aa30613067->SetStats(0);
   Graph_Graph_bin_0_syst_aa30613067->SetLineStyle(0);
   Graph_Graph_bin_0_syst_aa30613067->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_aa30613067->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_aa30613067->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa30613067->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa30613067->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa30613067->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa30613067->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa30613067->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_aa30613067);
   
   grae->Draw("5");
   
   Double_t bin_0_aa_fx3068[12] = {
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
   Double_t bin_0_aa_fy3068[12] = {
   3.623797,
   1.817616,
   0.9446006,
   0.4765622,
   0.2020894,
   0.08894734,
   0.03721237,
   0.01616085,
   0.00688428,
   0.00257571,
   0.001016595,
   0.000241375};
   Double_t bin_0_aa_felx3068[12] = {
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
   Double_t bin_0_aa_fely3068[12] = {
   0.1177039,
   0.04980035,
   0.02422505,
   0.01081128,
   0.004702762,
   0.002491556,
   0.001239323,
   0.0006939758,
   0.0003142609,
   0.0001760927,
   9.849308e-05,
   2.858609e-05};
   Double_t bin_0_aa_fehx3068[12] = {
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
   Double_t bin_0_aa_fehy3068[12] = {
   0.1177039,
   0.04980035,
   0.02422505,
   0.01081128,
   0.004702762,
   0.002491556,
   0.001239323,
   0.0006939758,
   0.0003142609,
   0.0001760927,
   9.849308e-05,
   2.858609e-05};
   grae = new TGraphAsymmErrors(12,bin_0_aa_fx3068,bin_0_aa_fy3068,bin_0_aa_felx3068,bin_0_aa_fehx3068,bin_0_aa_fely3068,bin_0_aa_fehy3068);
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
   
   TH1F *Graph_Graph_bin_0_aa30623068 = new TH1F("Graph_Graph_bin_0_aa30623068","Graph",100,2.15,54.35);
   Graph_Graph_bin_0_aa30623068->SetMinimum(0.00019151);
   Graph_Graph_bin_0_aa30623068->SetMaximum(4.11563);
   Graph_Graph_bin_0_aa30623068->SetDirectory(0);
   Graph_Graph_bin_0_aa30623068->SetStats(0);
   Graph_Graph_bin_0_aa30623068->SetLineStyle(0);
   Graph_Graph_bin_0_aa30623068->SetMarkerStyle(20);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_aa30623068->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_aa30623068->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa30623068->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa30623068->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa30623068->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa30623068->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa30623068->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_aa30623068);
   
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
      tex = new TLatex(0.2,0.85,"Prompt J/#psi");
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
   
   TH1F *haxes_copy__28 = new TH1F("haxes_copy__28","haxes",1,0,50);
   haxes_copy__28->SetMinimum(0.0001);
   haxes_copy__28->SetMaximum(500);
   haxes_copy__28->SetDirectory(0);
   haxes_copy__28->SetLineStyle(0);
   haxes_copy__28->SetMarkerStyle(20);
   haxes_copy__28->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__28->GetXaxis()->CenterTitle(true);
   haxes_copy__28->GetXaxis()->SetLabelFont(42);
   haxes_copy__28->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__28->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__28->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__28->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__28->GetXaxis()->SetTitleFont(42);
   haxes_copy__28->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/dp_{T} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/dp_{T} (nb / GeV/#it{c})");
   haxes_copy__28->GetYaxis()->SetLabelFont(42);
   haxes_copy__28->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__28->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__28->GetYaxis()->SetTitleSize(0.043);
   haxes_copy__28->GetYaxis()->SetTitleOffset(1.67);
   haxes_copy__28->GetYaxis()->SetTitleFont(42);
   haxes_copy__28->GetZaxis()->SetLabelFont(42);
   haxes_copy__28->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__28->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__28->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__28->GetZaxis()->SetTitleFont(42);
   haxes_copy__28->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
