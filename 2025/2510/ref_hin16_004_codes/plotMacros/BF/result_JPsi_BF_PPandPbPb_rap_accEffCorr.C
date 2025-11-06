void result_JPsi_BF_PPandPbPb_rap_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:56 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.48,-0.135,2.52,0.99);
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
   
   TH1F *haxes__115 = new TH1F("haxes__115","haxes",1,0,2.4);
   haxes__115->SetMinimum(0);
   haxes__115->SetMaximum(0.9);
   haxes__115->SetLineStyle(0);
   haxes__115->SetMarkerStyle(20);
   haxes__115->GetXaxis()->SetTitle("|y|");
   haxes__115->GetXaxis()->CenterTitle(true);
   haxes__115->GetXaxis()->SetNdivisions(-306);
   haxes__115->GetXaxis()->SetLabelFont(42);
   haxes__115->GetXaxis()->SetLabelOffset(0.007);
   haxes__115->GetXaxis()->SetLabelSize(0.05);
   haxes__115->GetXaxis()->SetTitleSize(0.06);
   haxes__115->GetXaxis()->SetTitleOffset(0.9);
   haxes__115->GetXaxis()->SetTitleFont(42);
   haxes__115->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
   haxes__115->GetYaxis()->SetLabelFont(42);
   haxes__115->GetYaxis()->SetLabelOffset(0.007);
   haxes__115->GetYaxis()->SetLabelSize(0.05);
   haxes__115->GetYaxis()->SetTitleSize(0.06);
   haxes__115->GetYaxis()->SetTitleOffset(1.25);
   haxes__115->GetYaxis()->SetTitleFont(42);
   haxes__115->GetZaxis()->SetLabelFont(42);
   haxes__115->GetZaxis()->SetLabelOffset(0.007);
   haxes__115->GetZaxis()->SetLabelSize(0.05);
   haxes__115->GetZaxis()->SetTitleSize(0.06);
   haxes__115->GetZaxis()->SetTitleFont(42);
   haxes__115->Draw("");
   
   Double_t bin_0_syst_pp_fx3237[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_pp_fy3237[6] = {
   0.2637829,
   0.2603897,
   0.2584116,
   0.2486694,
   0.2425548,
   0.2326316};
   Double_t bin_0_syst_pp_felx3237[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_pp_fely3237[6] = {
   0.01462651,
   0.01674748,
   0.01604399,
   0.01733656,
   0.01645004,
   0.01690512};
   Double_t bin_0_syst_pp_fehx3237[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_pp_fehy3237[6] = {
   0.01462651,
   0.01674748,
   0.01604399,
   0.01733656,
   0.01645004,
   0.01690512};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_pp_fx3237,bin_0_syst_pp_fy3237,bin_0_syst_pp_felx3237,bin_0_syst_pp_fehx3237,bin_0_syst_pp_fely3237,bin_0_syst_pp_fehy3237);
   grae->SetName("bin_0_syst_pp");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_pp32353237 = new TH1F("Graph_Graph_bin_0_syst_pp32353237","Graph",100,0,2.46);
   Graph_Graph_bin_0_syst_pp32353237->SetMinimum(0.2094582);
   Graph_Graph_bin_0_syst_pp32353237->SetMaximum(0.2846777);
   Graph_Graph_bin_0_syst_pp32353237->SetDirectory(0);
   Graph_Graph_bin_0_syst_pp32353237->SetStats(0);
   Graph_Graph_bin_0_syst_pp32353237->SetLineStyle(0);
   Graph_Graph_bin_0_syst_pp32353237->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_pp32353237->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_pp32353237->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_pp32353237->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_pp32353237->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_pp32353237->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_pp32353237->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_pp32353237->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_pp32353237);
   
   grae->Draw("5");
   
   Double_t bin_0_pp_fx3238[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_pp_fy3238[6] = {
   0.2637829,
   0.2603897,
   0.2584116,
   0.2486694,
   0.2425548,
   0.2326316};
   Double_t bin_0_pp_felx3238[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_pp_fely3238[6] = {
   0.002121966,
   0.002017173,
   0.001928439,
   0.001827103,
   0.001878562,
   0.002823068};
   Double_t bin_0_pp_fehx3238[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_pp_fehy3238[6] = {
   0.002121966,
   0.002017173,
   0.001928439,
   0.001827103,
   0.001878562,
   0.002823068};
   grae = new TGraphAsymmErrors(6,bin_0_pp_fx3238,bin_0_pp_fy3238,bin_0_pp_felx3238,bin_0_pp_fehx3238,bin_0_pp_fely3238,bin_0_pp_fehy3238);
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
   
   TH1F *Graph_Graph_bin_0_pp32363238 = new TH1F("Graph_Graph_bin_0_pp32363238","Graph",100,0,2.64);
   Graph_Graph_bin_0_pp32363238->SetMinimum(0.2261989);
   Graph_Graph_bin_0_pp32363238->SetMaximum(0.2695145);
   Graph_Graph_bin_0_pp32363238->SetDirectory(0);
   Graph_Graph_bin_0_pp32363238->SetStats(0);
   Graph_Graph_bin_0_pp32363238->SetLineStyle(0);
   Graph_Graph_bin_0_pp32363238->SetMarkerStyle(20);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_pp32363238->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_pp32363238->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_pp32363238->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_pp32363238->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_pp32363238->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_pp32363238->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_pp32363238->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_pp32363238);
   
   grae->Draw("p");
   
   Double_t bin_0_syst_aa_fx3239[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_aa_fy3239[6] = {
   0.3104231,
   0.3236807,
   0.306487,
   0.3073475,
   0.3167214,
   0.2902728};
   Double_t bin_0_syst_aa_felx3239[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_aa_fely3239[6] = {
   0.01744004,
   0.0201911,
   0.02230256,
   0.02822874,
   0.03224202,
   0.02742674};
   Double_t bin_0_syst_aa_fehx3239[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_aa_fehy3239[6] = {
   0.01744004,
   0.0201911,
   0.02230256,
   0.02822874,
   0.03224202,
   0.02742674};
   grae = new TGraphAsymmErrors(6,bin_0_syst_aa_fx3239,bin_0_syst_aa_fy3239,bin_0_syst_aa_felx3239,bin_0_syst_aa_fehx3239,bin_0_syst_aa_fely3239,bin_0_syst_aa_fehy3239);
   grae->SetName("bin_0_syst_aa");
   grae->SetTitle("Graph");

   ci = 2088;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#990000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_bin_0_syst_aa32333239 = new TH1F("Graph_Graph_bin_0_syst_aa32333239","Graph",100,0,2.46);
   Graph_Graph_bin_0_syst_aa32333239->SetMinimum(0.2542343);
   Graph_Graph_bin_0_syst_aa32333239->SetMaximum(0.3575751);
   Graph_Graph_bin_0_syst_aa32333239->SetDirectory(0);
   Graph_Graph_bin_0_syst_aa32333239->SetStats(0);
   Graph_Graph_bin_0_syst_aa32333239->SetLineStyle(0);
   Graph_Graph_bin_0_syst_aa32333239->SetMarkerStyle(20);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_syst_aa32333239->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_syst_aa32333239->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_syst_aa32333239->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_syst_aa32333239->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_syst_aa32333239->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_syst_aa32333239->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_syst_aa32333239->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_syst_aa32333239);
   
   grae->Draw("5");
   
   Double_t bin_0_aa_fx3240[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_aa_fy3240[6] = {
   0.3104231,
   0.3236807,
   0.306487,
   0.3073475,
   0.3167214,
   0.2902728};
   Double_t bin_0_aa_felx3240[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_aa_fely3240[6] = {
   0.005608597,
   0.005472332,
   0.005383786,
   0.005481929,
   0.006148722,
   0.009516474};
   Double_t bin_0_aa_fehx3240[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_aa_fehy3240[6] = {
   0.005608597,
   0.005472332,
   0.005383786,
   0.005481929,
   0.006148722,
   0.009516474};
   grae = new TGraphAsymmErrors(6,bin_0_aa_fx3240,bin_0_aa_fy3240,bin_0_aa_felx3240,bin_0_aa_fehx3240,bin_0_aa_fely3240,bin_0_aa_fehy3240);
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
   
   TH1F *Graph_Graph_bin_0_aa32343240 = new TH1F("Graph_Graph_bin_0_aa32343240","Graph",100,0,2.64);
   Graph_Graph_bin_0_aa32343240->SetMinimum(0.2759167);
   Graph_Graph_bin_0_aa32343240->SetMaximum(0.3339927);
   Graph_Graph_bin_0_aa32343240->SetDirectory(0);
   Graph_Graph_bin_0_aa32343240->SetStats(0);
   Graph_Graph_bin_0_aa32343240->SetLineStyle(0);
   Graph_Graph_bin_0_aa32343240->SetMarkerStyle(20);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_bin_0_aa32343240->GetXaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_bin_0_aa32343240->GetYaxis()->SetTitleFont(42);
   Graph_Graph_bin_0_aa32343240->GetZaxis()->SetLabelFont(42);
   Graph_Graph_bin_0_aa32343240->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_bin_0_aa32343240->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_bin_0_aa32343240->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_bin_0_aa32343240->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_bin_0_aa32343240);
   
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
   TLatex *   tex = new TLatex(0.2,0.86,"6.5 < p_{T} < 50 GeV/c");
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
   
   TH1F *haxes_copy__116 = new TH1F("haxes_copy__116","haxes",1,0,2.4);
   haxes_copy__116->SetMinimum(0);
   haxes_copy__116->SetMaximum(0.9);
   haxes_copy__116->SetDirectory(0);
   haxes_copy__116->SetLineStyle(0);
   haxes_copy__116->SetMarkerStyle(20);
   haxes_copy__116->GetXaxis()->SetTitle("|y|");
   haxes_copy__116->GetXaxis()->CenterTitle(true);
   haxes_copy__116->GetXaxis()->SetNdivisions(-306);
   haxes_copy__116->GetXaxis()->SetLabelFont(42);
   haxes_copy__116->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__116->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__116->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__116->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__116->GetXaxis()->SetTitleFont(42);
   haxes_copy__116->GetYaxis()->SetTitle("Nonprompt J/#psi fraction");
   haxes_copy__116->GetYaxis()->SetLabelFont(42);
   haxes_copy__116->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__116->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__116->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__116->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__116->GetYaxis()->SetTitleFont(42);
   haxes_copy__116->GetZaxis()->SetLabelFont(42);
   haxes_copy__116->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__116->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__116->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__116->GetZaxis()->SetTitleFont(42);
   haxes_copy__116->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
