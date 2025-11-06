void result_JPsi_RAA_pt_prompt_1_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:50 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__3 = new TH1F("haxes__3","haxes",1,0,30);
   haxes__3->SetMinimum(0);
   haxes__3->SetMaximum(1.5);
   haxes__3->SetLineStyle(0);
   haxes__3->SetMarkerStyle(20);
   haxes__3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__3->GetXaxis()->CenterTitle(true);
   haxes__3->GetXaxis()->SetLabelFont(42);
   haxes__3->GetXaxis()->SetLabelOffset(0.007);
   haxes__3->GetXaxis()->SetLabelSize(0.05);
   haxes__3->GetXaxis()->SetTitleSize(0.06);
   haxes__3->GetXaxis()->SetTitleOffset(0.9);
   haxes__3->GetXaxis()->SetTitleFont(42);
   haxes__3->GetYaxis()->SetTitle("R_{AA}");
   haxes__3->GetYaxis()->SetLabelFont(42);
   haxes__3->GetYaxis()->SetLabelOffset(0.007);
   haxes__3->GetYaxis()->SetLabelSize(0.05);
   haxes__3->GetYaxis()->SetTitleSize(0.06);
   haxes__3->GetYaxis()->SetTitleOffset(1.25);
   haxes__3->GetYaxis()->SetTitleFont(42);
   haxes__3->GetZaxis()->SetLabelFont(42);
   haxes__3->GetZaxis()->SetLabelOffset(0.007);
   haxes__3->GetZaxis()->SetLabelSize(0.05);
   haxes__3->GetZaxis()->SetTitleSize(0.06);
   haxes__3->GetZaxis()->SetTitleFont(42);
   haxes__3->Draw("");
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
   
   Double_t bin_0_syst_fx3009[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_0_syst_fy3009[9] = {
   0.2629595,
   0.2362154,
   0.24172,
   0.2329079,
   0.2352974,
   0.2584777,
   0.2777025,
   0.308975,
   0.3676};
   Double_t bin_0_syst_felx3009[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fely3009[9] = {
   0.02857466,
   0.02092909,
   0.01669648,
   0.01411494,
   0.01560396,
   0.01470632,
   0.01330478,
   0.01738706,
   0.009285893};
   Double_t bin_0_syst_fehx3009[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fehy3009[9] = {
   0.02857466,
   0.02092909,
   0.01669648,
   0.01411494,
   0.01560396,
   0.01470632,
   0.01330478,
   0.01738706,
   0.009285893};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(9,bin_0_syst_fx3009,bin_0_syst_fy3009,bin_0_syst_felx3009,bin_0_syst_fehx3009,bin_0_syst_fely3009,bin_0_syst_fehy3009);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 1953;
   color = new TColor(ci, 1, 0.6, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cc0000");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3009 = new TH1F("Graph_bin_0_syst3009","Graph",100,3.1,43.9);
   Graph_bin_0_syst3009->SetMinimum(0.1991264);
   Graph_bin_0_syst3009->SetMaximum(0.3930458);
   Graph_bin_0_syst3009->SetDirectory(0);
   Graph_bin_0_syst3009->SetStats(0);
   Graph_bin_0_syst3009->SetLineStyle(0);
   Graph_bin_0_syst3009->SetMarkerStyle(20);
   Graph_bin_0_syst3009->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3009->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3009->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3009->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3009->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3009->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3009->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3009->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3009->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3009->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3009->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3009->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3009->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3009->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3009->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3009->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3009->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3009);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3010[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_0_fy3010[9] = {
   0.2629595,
   0.2362154,
   0.24172,
   0.2329079,
   0.2352974,
   0.2584777,
   0.2777025,
   0.308975,
   0.3676};
   Double_t bin_0_felx3010[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_0_fely3010[9] = {
   0.01858628,
   0.01279181,
   0.01212817,
   0.009938315,
   0.009937075,
   0.0137461,
   0.0135818,
   0.02091388,
   0.05498765};
   Double_t bin_0_fehx3010[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_0_fehy3010[9] = {
   0.01858628,
   0.01279181,
   0.01212817,
   0.009938315,
   0.009937075,
   0.0137461,
   0.0135818,
   0.02091388,
   0.05498765};
   grae = new TGraphAsymmErrors(9,bin_0_fx3010,bin_0_fy3010,bin_0_felx3010,bin_0_fehx3010,bin_0_fely3010,bin_0_fehy3010);
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
   
   TH1F *Graph_bin_03010 = new TH1F("Graph_bin_03010","Graph",100,2.15,54.35);
   Graph_bin_03010->SetMinimum(0.2030078);
   Graph_bin_03010->SetMaximum(0.4425494);
   Graph_bin_03010->SetDirectory(0);
   Graph_bin_03010->SetStats(0);
   Graph_bin_03010->SetLineStyle(0);
   Graph_bin_03010->SetMarkerStyle(20);
   Graph_bin_03010->GetXaxis()->SetLabelFont(42);
   Graph_bin_03010->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03010->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03010->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03010->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03010->GetXaxis()->SetTitleFont(42);
   Graph_bin_03010->GetYaxis()->SetLabelFont(42);
   Graph_bin_03010->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03010->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03010->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03010->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03010->GetYaxis()->SetTitleFont(42);
   Graph_bin_03010->GetZaxis()->SetLabelFont(42);
   Graph_bin_03010->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03010->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03010->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03010->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03010);
   
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
   
   Double_t bin_1_syst_fx3011[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_1_syst_fy3011[9] = {
   0.3655674,
   0.3430961,
   0.3304756,
   0.3302569,
   0.3395999,
   0.3551809,
   0.34159,
   0.4664485,
   0.5095057};
   Double_t bin_1_syst_felx3011[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_1_syst_fely3011[9] = {
   0.03583843,
   0.02715358,
   0.02181061,
   0.01973079,
   0.01823648,
   0.01977635,
   0.02096462,
   0.02144922,
   0.01172149};
   Double_t bin_1_syst_fehx3011[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_1_syst_fehy3011[9] = {
   0.03583843,
   0.02715358,
   0.02181061,
   0.01973079,
   0.01823648,
   0.01977635,
   0.02096462,
   0.02144922,
   0.01172149};
   grae = new TGraphAsymmErrors(9,bin_1_syst_fx3011,bin_1_syst_fy3011,bin_1_syst_felx3011,bin_1_syst_fehx3011,bin_1_syst_fely3011,bin_1_syst_fehy3011);
   grae->SetName("bin_1_syst");
   grae->SetTitle("Graph");

   ci = 1955;
   color = new TColor(ci, 1, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#cccc00");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_1_syst3011 = new TH1F("Graph_bin_1_syst3011","Graph",100,3.1,43.9);
   Graph_bin_1_syst3011->SetMinimum(0.2874087);
   Graph_bin_1_syst3011->SetMaximum(0.5424834);
   Graph_bin_1_syst3011->SetDirectory(0);
   Graph_bin_1_syst3011->SetStats(0);
   Graph_bin_1_syst3011->SetLineStyle(0);
   Graph_bin_1_syst3011->SetMarkerStyle(20);
   Graph_bin_1_syst3011->GetXaxis()->SetLabelFont(42);
   Graph_bin_1_syst3011->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3011->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3011->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3011->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_1_syst3011->GetXaxis()->SetTitleFont(42);
   Graph_bin_1_syst3011->GetYaxis()->SetLabelFont(42);
   Graph_bin_1_syst3011->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3011->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3011->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3011->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_1_syst3011->GetYaxis()->SetTitleFont(42);
   Graph_bin_1_syst3011->GetZaxis()->SetLabelFont(42);
   Graph_bin_1_syst3011->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3011->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3011->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3011->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_1_syst3011);
   
   grae->Draw("5");
   
   Double_t bin_1_fx3012[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_1_fy3012[9] = {
   0.3655674,
   0.3430961,
   0.3304756,
   0.3302569,
   0.3395999,
   0.3551809,
   0.34159,
   0.4664485,
   0.5095057};
   Double_t bin_1_felx3012[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_1_fely3012[9] = {
   0.01900666,
   0.0143484,
   0.01335744,
   0.01177337,
   0.01250223,
   0.01536474,
   0.0140314,
   0.02806655,
   0.06082041};
   Double_t bin_1_fehx3012[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_1_fehy3012[9] = {
   0.01900666,
   0.0143484,
   0.01335744,
   0.01177337,
   0.01250223,
   0.01536474,
   0.0140314,
   0.02806655,
   0.06082041};
   grae = new TGraphAsymmErrors(9,bin_1_fx3012,bin_1_fy3012,bin_1_felx3012,bin_1_fehx3012,bin_1_fely3012,bin_1_fehy3012);
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
   
   TH1F *Graph_bin_13012 = new TH1F("Graph_bin_13012","Graph",100,2.15,54.35);
   Graph_bin_13012->SetMinimum(0.2917973);
   Graph_bin_13012->SetMaximum(0.5956469);
   Graph_bin_13012->SetDirectory(0);
   Graph_bin_13012->SetStats(0);
   Graph_bin_13012->SetLineStyle(0);
   Graph_bin_13012->SetMarkerStyle(20);
   Graph_bin_13012->GetXaxis()->SetLabelFont(42);
   Graph_bin_13012->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_13012->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_13012->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_13012->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_13012->GetXaxis()->SetTitleFont(42);
   Graph_bin_13012->GetYaxis()->SetLabelFont(42);
   Graph_bin_13012->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_13012->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_13012->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_13012->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_13012->GetYaxis()->SetTitleFont(42);
   Graph_bin_13012->GetZaxis()->SetLabelFont(42);
   Graph_bin_13012->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_13012->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_13012->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_13012->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_13012);
   
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
   
   Double_t bin_2_syst_fx3013[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_2_syst_fy3013[9] = {
   0.5365502,
   0.5699709,
   0.5133281,
   0.5362739,
   0.5616396,
   0.5793436,
   0.586059,
   0.702904,
   0.8988294};
   Double_t bin_2_syst_felx3013[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_2_syst_fely3013[9] = {
   0.05019049,
   0.04347909,
   0.03427159,
   0.03163308,
   0.03032763,
   0.03189337,
   0.02408806,
   0.02939232,
   0.01990119};
   Double_t bin_2_syst_fehx3013[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_2_syst_fehy3013[9] = {
   0.05019049,
   0.04347909,
   0.03427159,
   0.03163308,
   0.03032763,
   0.03189337,
   0.02408806,
   0.02939232,
   0.01990119};
   grae = new TGraphAsymmErrors(9,bin_2_syst_fx3013,bin_2_syst_fy3013,bin_2_syst_felx3013,bin_2_syst_fehx3013,bin_2_syst_fely3013,bin_2_syst_fehy3013);
   grae->SetName("bin_2_syst");
   grae->SetTitle("Graph");

   ci = 1957;
   color = new TColor(ci, 0.6, 0.8, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#3399ff");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_2_syst3013 = new TH1F("Graph_bin_2_syst3013","Graph",100,3.1,43.9);
   Graph_bin_2_syst3013->SetMinimum(0.4350891);
   Graph_bin_2_syst3013->SetMaximum(0.9626979);
   Graph_bin_2_syst3013->SetDirectory(0);
   Graph_bin_2_syst3013->SetStats(0);
   Graph_bin_2_syst3013->SetLineStyle(0);
   Graph_bin_2_syst3013->SetMarkerStyle(20);
   Graph_bin_2_syst3013->GetXaxis()->SetLabelFont(42);
   Graph_bin_2_syst3013->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3013->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3013->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3013->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_2_syst3013->GetXaxis()->SetTitleFont(42);
   Graph_bin_2_syst3013->GetYaxis()->SetLabelFont(42);
   Graph_bin_2_syst3013->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3013->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3013->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3013->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_2_syst3013->GetYaxis()->SetTitleFont(42);
   Graph_bin_2_syst3013->GetZaxis()->SetLabelFont(42);
   Graph_bin_2_syst3013->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_2_syst3013->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_2_syst3013->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_2_syst3013->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_2_syst3013);
   
   grae->Draw("5");
   
   Double_t bin_2_fx3014[9] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   17.5,
   25,
   40};
   Double_t bin_2_fy3014[9] = {
   0.5365502,
   0.5699709,
   0.5133281,
   0.5362739,
   0.5616396,
   0.5793436,
   0.586059,
   0.702904,
   0.8988294};
   Double_t bin_2_felx3014[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_2_fely3014[9] = {
   0.02534482,
   0.02552572,
   0.02177531,
   0.02105761,
   0.02322469,
   0.02847286,
   0.02787471,
   0.04968165,
   0.1221111};
   Double_t bin_2_fehx3014[9] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   2.5,
   5,
   10};
   Double_t bin_2_fehy3014[9] = {
   0.02534482,
   0.02552572,
   0.02177531,
   0.02105761,
   0.02322469,
   0.02847286,
   0.02787471,
   0.04968165,
   0.1221111};
   grae = new TGraphAsymmErrors(9,bin_2_fx3014,bin_2_fy3014,bin_2_felx3014,bin_2_fehx3014,bin_2_fely3014,bin_2_fehy3014);
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
   
   TH1F *Graph_bin_23014 = new TH1F("Graph_bin_23014","Graph",100,2.15,54.35);
   Graph_bin_23014->SetMinimum(0.438614);
   Graph_bin_23014->SetMaximum(1.073879);
   Graph_bin_23014->SetDirectory(0);
   Graph_bin_23014->SetStats(0);
   Graph_bin_23014->SetLineStyle(0);
   Graph_bin_23014->SetMarkerStyle(20);
   Graph_bin_23014->GetXaxis()->SetLabelFont(42);
   Graph_bin_23014->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_23014->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_23014->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_23014->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_23014->GetXaxis()->SetTitleFont(42);
   Graph_bin_23014->GetYaxis()->SetLabelFont(42);
   Graph_bin_23014->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_23014->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_23014->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_23014->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_23014->GetYaxis()->SetTitleFont(42);
   Graph_bin_23014->GetZaxis()->SetLabelFont(42);
   Graph_bin_23014->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_23014->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_23014->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_23014->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_23014);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.19,0.49,0.51,0.64,NULL,"brNDC");
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
      tex = new TLatex(0.2,0.85,"Prompt J/#psi");
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
   
   TH1F *haxes_copy__4 = new TH1F("haxes_copy__4","haxes",1,0,30);
   haxes_copy__4->SetMinimum(0);
   haxes_copy__4->SetMaximum(1.5);
   haxes_copy__4->SetDirectory(0);
   haxes_copy__4->SetLineStyle(0);
   haxes_copy__4->SetMarkerStyle(20);
   haxes_copy__4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__4->GetXaxis()->CenterTitle(true);
   haxes_copy__4->GetXaxis()->SetLabelFont(42);
   haxes_copy__4->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__4->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__4->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__4->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__4->GetXaxis()->SetTitleFont(42);
   haxes_copy__4->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__4->GetYaxis()->SetLabelFont(42);
   haxes_copy__4->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__4->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__4->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__4->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__4->GetYaxis()->SetTitleFont(42);
   haxes_copy__4->GetZaxis()->SetLabelFont(42);
   haxes_copy__4->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__4->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__4->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__4->GetZaxis()->SetTitleFont(42);
   haxes_copy__4->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
