void result_JPsi_RAA_pt_nonprompt_1_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:52 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-6,-0.225,31.5,1.65);
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
   
   TH1F *haxes__51 = new TH1F("haxes__51","haxes",1,0,30);
   haxes__51->SetMinimum(0);
   haxes__51->SetMaximum(1.5);
   haxes__51->SetLineStyle(0);
   haxes__51->SetMarkerStyle(20);
   haxes__51->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__51->GetXaxis()->CenterTitle(true);
   haxes__51->GetXaxis()->SetLabelFont(42);
   haxes__51->GetXaxis()->SetLabelOffset(0.007);
   haxes__51->GetXaxis()->SetLabelSize(0.05);
   haxes__51->GetXaxis()->SetTitleSize(0.06);
   haxes__51->GetXaxis()->SetTitleOffset(0.9);
   haxes__51->GetXaxis()->SetTitleFont(42);
   haxes__51->GetYaxis()->SetTitle("R_{AA}");
   haxes__51->GetYaxis()->SetLabelFont(42);
   haxes__51->GetYaxis()->SetLabelOffset(0.007);
   haxes__51->GetYaxis()->SetLabelSize(0.05);
   haxes__51->GetYaxis()->SetTitleSize(0.06);
   haxes__51->GetYaxis()->SetTitleOffset(1.25);
   haxes__51->GetYaxis()->SetTitleFont(42);
   haxes__51->GetZaxis()->SetLabelFont(42);
   haxes__51->GetZaxis()->SetLabelOffset(0.007);
   haxes__51->GetZaxis()->SetLabelSize(0.05);
   haxes__51->GetZaxis()->SetTitleSize(0.06);
   haxes__51->GetZaxis()->SetTitleFont(42);
   haxes__51->Draw("");
   TBox *box = new TBox(28.7,0.9572332,30,1.035917);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff9999");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#cc0000");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(28.7,0.9572332,30,1.035917);

   ci = TColor::GetColor("#ff9999");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#cc0000");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_0_syst_fx3107[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_0_syst_fy3107[9] = {
   0.3972783,
   0.3901201,
   0.343072,
   0.3645595,
   0.3570936,
   0.3474473,
   0.3522349,
   0.3885836,
   0.38451};
   Double_t bin_0_syst_felx3107[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fely3107[9] = {
   0.05615721,
   0.04135274,
   0.02652424,
   0.02309402,
   0.02614554,
   0.02058318,
   0.01808767,
   0.0190539,
   0.009381019};
   Double_t bin_0_syst_fehx3107[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fehy3107[9] = {
   0.05615721,
   0.04135274,
   0.02652424,
   0.02309402,
   0.02614554,
   0.02058318,
   0.01808767,
   0.0190539,
   0.009381019};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9,bin_0_syst_fx3107,bin_0_syst_fy3107,bin_0_syst_felx3107,bin_0_syst_fehx3107,bin_0_syst_fely3107,bin_0_syst_fehy3107);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 2015;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3107 = new TH1F("Graph_bin_0_syst3107","Graph",100,3.1,43.9);
   Graph_bin_0_syst3107->SetMinimum(0.302859);
   Graph_bin_0_syst3107->SetMaximum(0.4671243);
   Graph_bin_0_syst3107->SetDirectory(0);
   Graph_bin_0_syst3107->SetStats(0);
   Graph_bin_0_syst3107->SetLineStyle(0);
   Graph_bin_0_syst3107->SetMarkerStyle(20);
   Graph_bin_0_syst3107->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3107->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3107->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3107->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3107->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3107->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3107->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3107->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3107->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3107->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3107->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3107->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3107->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3107->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3107->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3107->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3107->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3107);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3108[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_0_fy3108[9] = {
   0.3972783,
   0.3901201,
   0.343072,
   0.3645595,
   0.3570936,
   0.3474473,
   0.3522349,
   0.3885836,
   0.38451};
   Double_t bin_0_felx3108[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_0_fely3108[9] = {
   0.02723873,
   0.02072604,
   0.01724477,
   0.01550451,
   0.01519375,
   0.01878054,
   0.0174931,
   0.02635767,
   0.05712663};
   Double_t bin_0_fehx3108[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_0_fehy3108[9] = {
   0.02723873,
   0.02072604,
   0.01724477,
   0.01550451,
   0.01519375,
   0.01878054,
   0.0174931,
   0.02635767,
   0.05712663};
   grae = new TGraphAsymmErrors(9,bin_0_fx3108,bin_0_fy3108,bin_0_felx3108,bin_0_fehx3108,bin_0_fely3108,bin_0_fehy3108);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#cc0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03108 = new TH1F("Graph_bin_03108","Graph",100,2.15,54.35);
   Graph_bin_03108->SetMinimum(0.3142463);
   Graph_bin_03108->SetMaximum(0.4532176);
   Graph_bin_03108->SetDirectory(0);
   Graph_bin_03108->SetStats(0);
   Graph_bin_03108->SetLineStyle(0);
   Graph_bin_03108->SetMarkerStyle(20);
   Graph_bin_03108->GetXaxis()->SetLabelFont(42);
   Graph_bin_03108->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03108->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03108->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03108->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03108->GetXaxis()->SetTitleFont(42);
   Graph_bin_03108->GetYaxis()->SetLabelFont(42);
   Graph_bin_03108->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03108->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03108->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03108->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03108->GetYaxis()->SetTitleFont(42);
   Graph_bin_03108->GetZaxis()->SetLabelFont(42);
   Graph_bin_03108->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03108->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03108->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03108->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03108);
   
   grae->Draw("p");
   box = new TBox(27.4,0.9543382,28.7,1.040062);

   ci = TColor::GetColor("#ffff99");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#cccc00");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(27.4,0.9543382,28.7,1.040062);

   ci = TColor::GetColor("#ffff99");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#cccc00");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_1_syst_fx3109[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_1_syst_fy3109[9] = {
   0.4816346,
   0.4988355,
   0.4747708,
   0.4404718,
   0.436888,
   0.4944584,
   0.4500249,
   0.474224,
   0.5962497};
   Double_t bin_1_syst_felx3109[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_1_syst_fely3109[9] = {
   0.05523483,
   0.04549913,
   0.03468886,
   0.02923866,
   0.02357309,
   0.02796476,
   0.02611533,
   0.02293423,
   0.01265982};
   Double_t bin_1_syst_fehx3109[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_1_syst_fehy3109[9] = {
   0.05523483,
   0.04549913,
   0.03468886,
   0.02923866,
   0.02357309,
   0.02796476,
   0.02611533,
   0.02293423,
   0.01265982};
   grae = new TGraphAsymmErrors(9,bin_1_syst_fx3109,bin_1_syst_fy3109,bin_1_syst_felx3109,bin_1_syst_fehx3109,bin_1_syst_fely3109,bin_1_syst_fehy3109);
   grae->SetName("bin_1_syst");
   grae->SetTitle("Graph");

   ci = 2017;
   color = new TColor(ci, 1, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cccc00");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_1_syst3109 = new TH1F("Graph_bin_1_syst3109","Graph",100,3.1,43.9);
   Graph_bin_1_syst3109->SetMinimum(0.3914655);
   Graph_bin_1_syst3109->SetMaximum(0.6286772);
   Graph_bin_1_syst3109->SetDirectory(0);
   Graph_bin_1_syst3109->SetStats(0);
   Graph_bin_1_syst3109->SetLineStyle(0);
   Graph_bin_1_syst3109->SetMarkerStyle(20);
   Graph_bin_1_syst3109->GetXaxis()->SetLabelFont(42);
   Graph_bin_1_syst3109->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3109->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3109->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3109->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_1_syst3109->GetXaxis()->SetTitleFont(42);
   Graph_bin_1_syst3109->GetYaxis()->SetLabelFont(42);
   Graph_bin_1_syst3109->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3109->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3109->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3109->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_1_syst3109->GetYaxis()->SetTitleFont(42);
   Graph_bin_1_syst3109->GetZaxis()->SetLabelFont(42);
   Graph_bin_1_syst3109->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3109->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3109->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3109->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_1_syst3109);
   
   grae->Draw("5");
   
   Double_t bin_1_fx3110[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_1_fy3110[9] = {
   0.4816346,
   0.4988355,
   0.4747708,
   0.4404718,
   0.436888,
   0.4944584,
   0.4500249,
   0.474224,
   0.5962497};
   Double_t bin_1_felx3110[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_1_fely3110[9] = {
   0.02485999,
   0.02090995,
   0.01929581,
   0.01584346,
   0.01632404,
   0.02159524,
   0.01878024,
   0.02891143,
   0.06985855};
   Double_t bin_1_fehx3110[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_1_fehy3110[9] = {
   0.02485999,
   0.02090995,
   0.01929581,
   0.01584346,
   0.01632404,
   0.02159524,
   0.01878024,
   0.02891143,
   0.06985855};
   grae = new TGraphAsymmErrors(9,bin_1_fx3110,bin_1_fy3110,bin_1_felx3110,bin_1_fehx3110,bin_1_fely3110,bin_1_fehy3110);
   grae->SetName("bin_1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#cccc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#cccc00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(29);
   grae->SetMarkerSize(2.3);
   
   TH1F *Graph_bin_13110 = new TH1F("Graph_bin_13110","Graph",100,2.15,54.35);
   Graph_bin_13110->SetMinimum(0.3960095);
   Graph_bin_13110->SetMaximum(0.6906627);
   Graph_bin_13110->SetDirectory(0);
   Graph_bin_13110->SetStats(0);
   Graph_bin_13110->SetLineStyle(0);
   Graph_bin_13110->SetMarkerStyle(20);
   Graph_bin_13110->GetXaxis()->SetLabelFont(42);
   Graph_bin_13110->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_13110->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_13110->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_13110->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_13110->GetXaxis()->SetTitleFont(42);
   Graph_bin_13110->GetYaxis()->SetLabelFont(42);
   Graph_bin_13110->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_13110->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_13110->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_13110->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_13110->GetYaxis()->SetTitleFont(42);
   Graph_bin_13110->GetZaxis()->SetLabelFont(42);
   Graph_bin_13110->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_13110->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_13110->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_13110->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_13110);
   
   grae->Draw("p");
   box = new TBox(26.1,0.9464743,27.4,1.073607);

   ci = TColor::GetColor("#99ccff");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#3399ff");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(26.1,0.9464743,27.4,1.073607);

   ci = TColor::GetColor("#99ccff");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#3399ff");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_2_syst_fx3111[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_2_syst_fy3111[9] = {
   0.7746424,
   0.6354407,
   0.6862662,
   0.5937709,
   0.5619589,
   0.6425758,
   0.5779454,
   0.5505761,
   0.6833567};
   Double_t bin_2_syst_felx3111[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_2_syst_fely3111[9] = {
   0.07954339,
   0.05056118,
   0.05417303,
   0.03981346,
   0.03162947,
   0.0384276,
   0.02454458,
   0.02345916,
   0.01336635};
   Double_t bin_2_syst_fehx3111[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_2_syst_fehy3111[9] = {
   0.07954339,
   0.05056118,
   0.05417303,
   0.03981346,
   0.03162947,
   0.0384276,
   0.02454458,
   0.02345916,
   0.01336635};
   grae = new TGraphAsymmErrors(9,bin_2_syst_fx3111,bin_2_syst_fy3111,bin_2_syst_felx3111,bin_2_syst_fehx3111,bin_2_syst_fely3111,bin_2_syst_fehy3111);
   grae->SetName("bin_2_syst");
   grae->SetTitle("Graph");

   ci = 2019;
   color = new TColor(ci, 0.6, 0.8, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#3399ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_2_syst3111 = new TH1F("Graph_bin_2_syst3111","Graph",100,3.1,43.9);
   Graph_bin_2_syst3111->SetMinimum(0.4944101);
   Graph_bin_2_syst3111->SetMaximum(0.8868927);
   Graph_bin_2_syst3111->SetDirectory(0);
   Graph_bin_2_syst3111->SetStats(0);
   Graph_bin_2_syst3111->SetLineStyle(0);
   Graph_bin_2_syst3111->SetMarkerStyle(20);
   Graph_bin_2_syst3111->GetXaxis()->SetLabelFont(42);
   Graph_bin_2_syst3111->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3111->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3111->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3111->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_2_syst3111->GetXaxis()->SetTitleFont(42);
   Graph_bin_2_syst3111->GetYaxis()->SetLabelFont(42);
   Graph_bin_2_syst3111->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3111->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3111->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3111->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_2_syst3111->GetYaxis()->SetTitleFont(42);
   Graph_bin_2_syst3111->GetZaxis()->SetLabelFont(42);
   Graph_bin_2_syst3111->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3111->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3111->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3111->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_2_syst3111);
   
   grae->Draw("5");
   
   Double_t bin_2_fx3112[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_2_fy3112[9] = {
   0.7746424,
   0.6354407,
   0.6862662,
   0.5937709,
   0.5619589,
   0.6425758,
   0.5779454,
   0.5505761,
   0.6833567};
   Double_t bin_2_felx3112[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_2_fely3112[9] = {
   0.03673768,
   0.02861854,
   0.0296127,
   0.02379338,
   0.02380384,
   0.03226085,
   0.02804164,
   0.03991335,
   0.09656166};
   Double_t bin_2_fehx3112[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_2_fehy3112[9] = {
   0.03673768,
   0.02861854,
   0.0296127,
   0.02379338,
   0.02380384,
   0.03226085,
   0.02804164,
   0.03991335,
   0.09656166};
   grae = new TGraphAsymmErrors(9,bin_2_fx3112,bin_2_fy3112,bin_2_felx3112,bin_2_fehx3112,bin_2_fely3112,bin_2_fehy3112);
   grae->SetName("bin_2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#3399ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#3399ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(34);
   grae->SetMarkerSize(2);
   
   TH1F *Graph_bin_23112 = new TH1F("Graph_bin_23112","Graph",100,2.15,54.35);
   Graph_bin_23112->SetMinimum(0.480591);
   Graph_bin_23112->SetMaximum(0.8414518);
   Graph_bin_23112->SetDirectory(0);
   Graph_bin_23112->SetStats(0);
   Graph_bin_23112->SetLineStyle(0);
   Graph_bin_23112->SetMarkerStyle(20);
   Graph_bin_23112->GetXaxis()->SetLabelFont(42);
   Graph_bin_23112->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_23112->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_23112->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_23112->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_23112->GetXaxis()->SetTitleFont(42);
   Graph_bin_23112->GetYaxis()->SetLabelFont(42);
   Graph_bin_23112->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_23112->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_23112->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_23112->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_23112->GetYaxis()->SetTitleFont(42);
   Graph_bin_23112->GetZaxis()->SetLabelFont(42);
   Graph_bin_23112->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_23112->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_23112->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_23112->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_23112);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.56,0.47,0.88,0.62,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("bin_0","Cent. 0-10%","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("bin_1","Cent. 10-30%","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cccc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(2.3);
   entry->SetTextFont(42);
   entry=leg->AddEntry("bin_2","Cent. 30-100%","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3399ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(2);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,30,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.044);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.85,"J/#psi from b hadrons");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.057);
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
   
   TH1F *haxes_copy__52 = new TH1F("haxes_copy__52","haxes",1,0,30);
   haxes_copy__52->SetMinimum(0);
   haxes_copy__52->SetMaximum(1.5);
   haxes_copy__52->SetDirectory(0);
   haxes_copy__52->SetLineStyle(0);
   haxes_copy__52->SetMarkerStyle(20);
   haxes_copy__52->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__52->GetXaxis()->CenterTitle(true);
   haxes_copy__52->GetXaxis()->SetLabelFont(42);
   haxes_copy__52->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__52->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__52->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__52->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__52->GetXaxis()->SetTitleFont(42);
   haxes_copy__52->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__52->GetYaxis()->SetLabelFont(42);
   haxes_copy__52->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__52->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__52->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__52->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__52->GetYaxis()->SetTitleFont(42);
   haxes_copy__52->GetZaxis()->SetLabelFont(42);
   haxes_copy__52->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__52->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__52->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__52->GetZaxis()->SetTitleFont(42);
   haxes_copy__52->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
