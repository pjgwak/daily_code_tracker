void result_JPsi_RAA_pt_prompt_2RapRanges_0_linearX_accEffCorr()
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
   
   TH1F *haxes__41 = new TH1F("haxes__41","haxes",1,0,30);
   haxes__41->SetMinimum(0);
   haxes__41->SetMaximum(1.5);
   haxes__41->SetLineStyle(0);
   haxes__41->SetMarkerStyle(20);
   haxes__41->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__41->GetXaxis()->CenterTitle(true);
   haxes__41->GetXaxis()->SetLabelFont(42);
   haxes__41->GetXaxis()->SetLabelOffset(0.007);
   haxes__41->GetXaxis()->SetLabelSize(0.05);
   haxes__41->GetXaxis()->SetTitleSize(0.06);
   haxes__41->GetXaxis()->SetTitleOffset(0.9);
   haxes__41->GetXaxis()->SetTitleFont(42);
   haxes__41->GetYaxis()->SetTitle("R_{AA}");
   haxes__41->GetYaxis()->SetLabelFont(42);
   haxes__41->GetYaxis()->SetLabelOffset(0.007);
   haxes__41->GetYaxis()->SetLabelSize(0.05);
   haxes__41->GetYaxis()->SetTitleSize(0.06);
   haxes__41->GetYaxis()->SetTitleOffset(1.25);
   haxes__41->GetYaxis()->SetTitleFont(42);
   haxes__41->GetZaxis()->SetLabelFont(42);
   haxes__41->GetZaxis()->SetLabelOffset(0.007);
   haxes__41->GetZaxis()->SetLabelSize(0.05);
   haxes__41->GetZaxis()->SetTitleSize(0.06);
   haxes__41->GetZaxis()->SetTitleFont(42);
   haxes__41->Draw("");
   TBox *box = new TBox(28.7,0.9543382,30,1.041388);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(28.7,0.9543382,30,1.041388);

   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_0_syst_fx3083[7] = {
   7.5,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_0_syst_fy3083[7] = {
   0.3351518,
   0.3470662,
   0.3301725,
   0.3502744,
   0.3585818,
   0.4272142,
   0.4723112};
   Double_t bin_0_syst_felx3083[7] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fely3083[7] = {
   0.03755993,
   0.02501591,
   0.01900603,
   0.01793556,
   0.01748639,
   0.02209087,
   0.01153057};
   Double_t bin_0_syst_fehx3083[7] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fehy3083[7] = {
   0.03755993,
   0.02501591,
   0.01900603,
   0.01793556,
   0.01748639,
   0.02209087,
   0.01153057};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(7,bin_0_syst_fx3083,bin_0_syst_fy3083,bin_0_syst_felx3083,bin_0_syst_fehx3083,bin_0_syst_fely3083,bin_0_syst_fehy3083);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 1995;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3083 = new TH1F("Graph_bin_0_syst3083","Graph",100,3.65,43.85);
   Graph_bin_0_syst3083->SetMinimum(0.2789669);
   Graph_bin_0_syst3083->SetMaximum(0.5024668);
   Graph_bin_0_syst3083->SetDirectory(0);
   Graph_bin_0_syst3083->SetStats(0);
   Graph_bin_0_syst3083->SetLineStyle(0);
   Graph_bin_0_syst3083->SetMarkerStyle(20);
   Graph_bin_0_syst3083->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3083->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3083->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3083->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3083->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3083->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3083->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3083->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3083->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3083->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3083->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3083->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3083->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3083->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3083->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3083->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3083->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3083);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3084[7] = {
   7.5,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_0_fy3084[7] = {
   0.3351518,
   0.3470662,
   0.3301725,
   0.3502744,
   0.3585818,
   0.4272142,
   0.4723112};
   Double_t bin_0_felx3084[7] = {
   1,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_0_fely3084[7] = {
   0.01877118,
   0.01796832,
   0.01525493,
   0.01188461,
   0.01744414,
   0.02880854,
   0.06535358};
   Double_t bin_0_fehx3084[7] = {
   1,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_0_fehy3084[7] = {
   0.01877118,
   0.01796832,
   0.01525493,
   0.01188461,
   0.01744414,
   0.02880854,
   0.06535358};
   grae = new TGraphAsymmErrors(7,bin_0_fx3084,bin_0_fy3084,bin_0_felx3084,bin_0_fehx3084,bin_0_fely3084,bin_0_fehy3084);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#000099");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03084 = new TH1F("Graph_bin_03084","Graph",100,2.15,54.35);
   Graph_bin_03084->SetMinimum(0.2926428);
   Graph_bin_03084->SetMaximum(0.5599395);
   Graph_bin_03084->SetDirectory(0);
   Graph_bin_03084->SetStats(0);
   Graph_bin_03084->SetLineStyle(0);
   Graph_bin_03084->SetMarkerStyle(20);
   Graph_bin_03084->GetXaxis()->SetLabelFont(42);
   Graph_bin_03084->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03084->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03084->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03084->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03084->GetXaxis()->SetTitleFont(42);
   Graph_bin_03084->GetYaxis()->SetLabelFont(42);
   Graph_bin_03084->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03084->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03084->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03084->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03084->GetYaxis()->SetTitleFont(42);
   Graph_bin_03084->GetZaxis()->SetLabelFont(42);
   Graph_bin_03084->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03084->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03084->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03084->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03084);
   
   grae->Draw("p");
   box = new TBox(27.4,0.9543382,28.7,1.041388);

   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(27.4,0.9543382,28.7,1.041388);

   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_3_syst_fx3085[11] = {
   3.75,
   5,
   6,
   7,
   8,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_3_syst_fy3085[11] = {
   0.5533585,
   0.3652605,
   0.3276063,
   0.3335032,
   0.3011363,
   0.3547365,
   0.2906677,
   0.3173846,
   0.358335,
   0.2983698,
   0.33122};
   Double_t bin_3_syst_felx3085[11] = {
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
   Double_t bin_3_syst_fely3085[11] = {
   0.1132607,
   0.04955064,
   0.03594794,
   0.02754959,
   0.02375441,
   0.02744593,
   0.01963667,
   0.01947041,
   0.02954741,
   0.03365004,
   0.02719949};
   Double_t bin_3_syst_fehx3085[11] = {
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
   Double_t bin_3_syst_fehy3085[11] = {
   0.1132607,
   0.04955064,
   0.03594794,
   0.02754959,
   0.02375441,
   0.02744593,
   0.01963667,
   0.01947041,
   0.02954741,
   0.03365004,
   0.02719949};
   grae = new TGraphAsymmErrors(11,bin_3_syst_fx3085,bin_3_syst_fy3085,bin_3_syst_felx3085,bin_3_syst_fehx3085,bin_3_syst_fely3085,bin_3_syst_fehy3085);
   grae->SetName("bin_3_syst");
   grae->SetTitle("Graph");

   ci = 1999;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_3_syst3085 = new TH1F("Graph_bin_3_syst3085","Graph",100,0,44.225);
   Graph_bin_3_syst3085->SetMinimum(0.2245298);
   Graph_bin_3_syst3085->SetMaximum(0.7068092);
   Graph_bin_3_syst3085->SetDirectory(0);
   Graph_bin_3_syst3085->SetStats(0);
   Graph_bin_3_syst3085->SetLineStyle(0);
   Graph_bin_3_syst3085->SetMarkerStyle(20);
   Graph_bin_3_syst3085->GetXaxis()->SetLabelFont(42);
   Graph_bin_3_syst3085->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3085->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3085->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3085->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_3_syst3085->GetXaxis()->SetTitleFont(42);
   Graph_bin_3_syst3085->GetYaxis()->SetLabelFont(42);
   Graph_bin_3_syst3085->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3085->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3085->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3085->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_3_syst3085->GetYaxis()->SetTitleFont(42);
   Graph_bin_3_syst3085->GetZaxis()->SetLabelFont(42);
   Graph_bin_3_syst3085->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3085->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3085->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3085->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_3_syst3085);
   
   grae->Draw("5");
   
   Double_t bin_3_fx3086[11] = {
   3.75,
   5,
   6,
   7,
   8,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_3_fy3086[11] = {
   0.5533585,
   0.3652605,
   0.3276063,
   0.3335032,
   0.3011363,
   0.3547365,
   0.2906677,
   0.3173846,
   0.358335,
   0.2983698,
   0.33122};
   Double_t bin_3_felx3086[11] = {
   0.75,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_3_fely3086[11] = {
   0.04663241,
   0.02444664,
   0.02203498,
   0.02022131,
   0.02105545,
   0.02627704,
   0.02038991,
   0.01825855,
   0.03277489,
   0.04209021,
   0.09205794};
   Double_t bin_3_fehx3086[11] = {
   0.75,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_3_fehy3086[11] = {
   0.04663241,
   0.02444664,
   0.02203498,
   0.02022131,
   0.02105545,
   0.02627704,
   0.02038991,
   0.01825855,
   0.03277489,
   0.04209021,
   0.09205794};
   grae = new TGraphAsymmErrors(11,bin_3_fx3086,bin_3_fy3086,bin_3_felx3086,bin_3_fehx3086,bin_3_fely3086,bin_3_fehy3086);
   grae->SetName("bin_3");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#009900");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(2.2);
   
   TH1F *Graph_bin_33086 = new TH1F("Graph_bin_33086","Graph",100,0,54.7);
   Graph_bin_33086->SetMinimum(0.2030792);
   Graph_bin_33086->SetMaximum(0.6360738);
   Graph_bin_33086->SetDirectory(0);
   Graph_bin_33086->SetStats(0);
   Graph_bin_33086->SetLineStyle(0);
   Graph_bin_33086->SetMarkerStyle(20);
   Graph_bin_33086->GetXaxis()->SetLabelFont(42);
   Graph_bin_33086->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_33086->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_33086->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_33086->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_33086->GetXaxis()->SetTitleFont(42);
   Graph_bin_33086->GetYaxis()->SetLabelFont(42);
   Graph_bin_33086->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_33086->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_33086->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_33086->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_33086->GetYaxis()->SetTitleFont(42);
   Graph_bin_33086->GetZaxis()->SetLabelFont(42);
   Graph_bin_33086->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_33086->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_33086->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_33086->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_33086);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.56,0.47,0.88,0.62,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("bin_0","|y| < 0.6","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#000099");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("bin_3","1.8 < |y| < 2.4","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(2.2);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,30,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"Cent. 0-100%");
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
   
   TH1F *haxes_copy__42 = new TH1F("haxes_copy__42","haxes",1,0,30);
   haxes_copy__42->SetMinimum(0);
   haxes_copy__42->SetMaximum(1.5);
   haxes_copy__42->SetDirectory(0);
   haxes_copy__42->SetLineStyle(0);
   haxes_copy__42->SetMarkerStyle(20);
   haxes_copy__42->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__42->GetXaxis()->CenterTitle(true);
   haxes_copy__42->GetXaxis()->SetLabelFont(42);
   haxes_copy__42->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__42->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__42->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__42->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__42->GetXaxis()->SetTitleFont(42);
   haxes_copy__42->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__42->GetYaxis()->SetLabelFont(42);
   haxes_copy__42->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__42->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__42->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__42->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__42->GetYaxis()->SetTitleFont(42);
   haxes_copy__42->GetZaxis()->SetLabelFont(42);
   haxes_copy__42->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__42->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__42->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__42->GetZaxis()->SetTitleFont(42);
   haxes_copy__42->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
