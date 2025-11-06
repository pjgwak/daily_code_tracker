void result_JPsi_XS_PPandPbPb_rap_nonprompt_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.48,-0.17885,2.52,1.3199);
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
   
   TH1F *haxes__81 = new TH1F("haxes__81","haxes",1,0,2.4);
   haxes__81->SetMinimum(0.001);
   haxes__81->SetMaximum(1.2);
   haxes__81->SetLineStyle(0);
   haxes__81->SetMarkerStyle(20);
   haxes__81->GetXaxis()->SetTitle("|y|");
   haxes__81->GetXaxis()->CenterTitle(true);
   haxes__81->GetXaxis()->SetNdivisions(-306);
   haxes__81->GetXaxis()->SetLabelFont(42);
   haxes__81->GetXaxis()->SetLabelOffset(0.007);
   haxes__81->GetXaxis()->SetLabelSize(0.05);
   haxes__81->GetXaxis()->SetTitleSize(0.06);
   haxes__81->GetXaxis()->SetTitleOffset(0.9);
   haxes__81->GetXaxis()->SetTitleFont(42);
   haxes__81->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/d#it{y} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/d#it{y} (nb)");
   haxes__81->GetYaxis()->SetLabelFont(42);
   haxes__81->GetYaxis()->SetLabelOffset(0.007);
   haxes__81->GetYaxis()->SetLabelSize(0.05);
   haxes__81->GetYaxis()->SetTitleSize(0.05);
   haxes__81->GetYaxis()->SetTitleOffset(1.35);
   haxes__81->GetYaxis()->SetTitleFont(42);
   haxes__81->GetZaxis()->SetLabelFont(42);
   haxes__81->GetZaxis()->SetLabelOffset(0.007);
   haxes__81->GetZaxis()->SetLabelSize(0.05);
   haxes__81->GetZaxis()->SetTitleSize(0.06);
   haxes__81->GetZaxis()->SetTitleFont(42);
   haxes__81->Draw("");
   
   Double_t bin_0_syst_pp_fx3171[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_pp_fy3171[6] = {
   0.2108696,
   0.2061499,
   0.199426,
   0.1713254,
   0.1469834,
   0.1156854};
   Double_t bin_0_syst_pp_felx3171[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_pp_fely3171[6] = {
   0.008442404,
   0.01000076,
   0.009194349,
   0.008217691,
   0.00687868,
   0.005148265};
   Double_t bin_0_syst_pp_fehx3171[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_pp_fehy3171[6] = {
   0.008442404,
   0.01000076,
   0.009194349,
   0.008217691,
   0.00687868,
   0.005148265};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_pp_fx3171,bin_0_syst_pp_fy3171,bin_0_syst_pp_felx3171,bin_0_syst_pp_fehx3171,bin_0_syst_pp_fely3171,bin_0_syst_pp_fehy3171);
   grae->SetName("bin_0_syst_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_pp31693171 = new TH1F("Graph_Graph_bin_0_syst_pp31693171","Graph",100,0,2.46);
   Graph_Graph_bin_0_syst_pp31693171->SetMinimum(0.09965959);
   Graph_Graph_bin_0_syst_pp31693171->SetMaximum(0.2301895);
   Graph_Graph_bin_0_syst_pp31693171->SetDirectory(0);
   Graph_Graph_bin_0_syst_pp31693171->SetStats(0);
   Graph_Graph_bin_0_syst_pp31693171->SetLineStyle(0);
   Graph_Graph_bin_0_syst_pp31693171->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_pp31693171->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_pp31693171->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp31693171->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp31693171->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp31693171->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp31693171->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp31693171->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_pp31693171);
   
   grae->Draw("5");
   
   Double_t bin_0_pp_fx3172[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_pp_fy3172[6] = {
   0.2108696,
   0.2061499,
   0.199426,
   0.1713254,
   0.1469834,
   0.1156854};
   Double_t bin_0_pp_felx3172[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_pp_fely3172[6] = {
   0.001716515,
   0.001901429,
   0.001775926,
   0.001453643,
   0.001349139,
   0.001607842};
   Double_t bin_0_pp_fehx3172[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_pp_fehy3172[6] = {
   0.001716515,
   0.001901429,
   0.001775926,
   0.001453643,
   0.001349139,
   0.001607842};
   grae = new TGraphAsymmErrors(6,bin_0_pp_fx3172,bin_0_pp_fy3172,bin_0_pp_felx3172,bin_0_pp_fehx3172,bin_0_pp_fely3172,bin_0_pp_fehy3172);
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
   
   TH1F *Graph_Graph_bin_0_pp31703172 = new TH1F("Graph_Graph_bin_0_pp31703172","Graph",100,0,2.64);
   Graph_Graph_bin_0_pp31703172->SetMinimum(0.1042266);
   Graph_Graph_bin_0_pp31703172->SetMaximum(0.222437);
   Graph_Graph_bin_0_pp31703172->SetDirectory(0);
   Graph_Graph_bin_0_pp31703172->SetStats(0);
   Graph_Graph_bin_0_pp31703172->SetLineStyle(0);
   Graph_Graph_bin_0_pp31703172->SetMarkerStyle(20);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_pp31703172->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_pp31703172->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp31703172->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp31703172->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp31703172->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp31703172->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp31703172->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_pp31703172);
   
   grae->Draw("p");
   
   Double_t bin_0_syst_aa_fx3173[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_aa_fy3173[6] = {
   0.09437293,
   0.09828211,
   0.08716533,
   0.07684412,
   0.06895147,
   0.04966052};
   Double_t bin_0_syst_aa_felx3173[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_aa_fely3173[6] = {
   0.004523822,
   0.004955713,
   0.004996261,
   0.005846902,
   0.005453679,
   0.003918568};
   Double_t bin_0_syst_aa_fehx3173[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_aa_fehy3173[6] = {
   0.004523822,
   0.004955713,
   0.004996261,
   0.005846902,
   0.005453679,
   0.003918568};
   grae = new TGraphAsymmErrors(6,bin_0_syst_aa_fx3173,bin_0_syst_aa_fy3173,bin_0_syst_aa_felx3173,bin_0_syst_aa_fehx3173,bin_0_syst_aa_fely3173,bin_0_syst_aa_fehy3173);
   grae->SetName("bin_0_syst_aa");
   grae->SetTitle("Graph");

   ci = 2056;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_aa31673173 = new TH1F("Graph_Graph_bin_0_syst_aa31673173","Graph",100,0,2.46);
   Graph_Graph_bin_0_syst_aa31673173->SetMinimum(0.03999237);
   Graph_Graph_bin_0_syst_aa31673173->SetMaximum(0.1089874);
   Graph_Graph_bin_0_syst_aa31673173->SetDirectory(0);
   Graph_Graph_bin_0_syst_aa31673173->SetStats(0);
   Graph_Graph_bin_0_syst_aa31673173->SetLineStyle(0);
   Graph_Graph_bin_0_syst_aa31673173->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_aa31673173->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_aa31673173->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa31673173->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa31673173->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa31673173->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa31673173->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa31673173->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_aa31673173);
   
   grae->Draw("5");
   
   Double_t bin_0_aa_fx3174[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_aa_fy3174[6] = {
   0.09437293,
   0.09828211,
   0.08716533,
   0.07684412,
   0.06895147,
   0.04966052};
   Double_t bin_0_aa_felx3174[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_aa_fely3174[6] = {
   0.002076375,
   0.002051138,
   0.0018696,
   0.001671622,
   0.001603945,
   0.001907233};
   Double_t bin_0_aa_fehx3174[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_aa_fehy3174[6] = {
   0.002076375,
   0.002051138,
   0.0018696,
   0.001671622,
   0.001603945,
   0.001907233};
   grae = new TGraphAsymmErrors(6,bin_0_aa_fx3174,bin_0_aa_fy3174,bin_0_aa_felx3174,bin_0_aa_fehx3174,bin_0_aa_fely3174,bin_0_aa_fehy3174);
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
   
   TH1F *Graph_Graph_bin_0_aa31683174 = new TH1F("Graph_Graph_bin_0_aa31683174","Graph",100,0,2.64);
   Graph_Graph_bin_0_aa31683174->SetMinimum(0.04249529);
   Graph_Graph_bin_0_aa31683174->SetMaximum(0.1055912);
   Graph_Graph_bin_0_aa31683174->SetDirectory(0);
   Graph_Graph_bin_0_aa31683174->SetStats(0);
   Graph_Graph_bin_0_aa31683174->SetLineStyle(0);
   Graph_Graph_bin_0_aa31683174->SetMarkerStyle(20);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_aa31683174->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_aa31683174->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa31683174->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa31683174->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa31683174->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa31683174->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa31683174->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_aa31683174);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.16,0.55,0.46,0.71,NULL,"brNDC");
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
      tex = new TLatex(0.2,0.75,"6.5 < p_{T} < 50 GeV/c");
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
   
   TH1F *haxes_copy__82 = new TH1F("haxes_copy__82","haxes",1,0,2.4);
   haxes_copy__82->SetMinimum(0.001);
   haxes_copy__82->SetMaximum(1.2);
   haxes_copy__82->SetDirectory(0);
   haxes_copy__82->SetLineStyle(0);
   haxes_copy__82->SetMarkerStyle(20);
   haxes_copy__82->GetXaxis()->SetTitle("|y|");
   haxes_copy__82->GetXaxis()->CenterTitle(true);
   haxes_copy__82->GetXaxis()->SetNdivisions(-306);
   haxes_copy__82->GetXaxis()->SetLabelFont(42);
   haxes_copy__82->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__82->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__82->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__82->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__82->GetXaxis()->SetTitleFont(42);
   haxes_copy__82->GetYaxis()->SetTitle("#bf{#it{#Beta}} #times d#sigma/d#it{y} or #bf{#it{#Beta}} #times (1/T_{AA} N_{MB}) dN/d#it{y} (nb)");
   haxes_copy__82->GetYaxis()->SetLabelFont(42);
   haxes_copy__82->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__82->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__82->GetYaxis()->SetTitleSize(0.05);
   haxes_copy__82->GetYaxis()->SetTitleOffset(1.35);
   haxes_copy__82->GetYaxis()->SetTitleFont(42);
   haxes_copy__82->GetZaxis()->SetLabelFont(42);
   haxes_copy__82->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__82->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__82->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__82->GetZaxis()->SetTitleFont(42);
   haxes_copy__82->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
