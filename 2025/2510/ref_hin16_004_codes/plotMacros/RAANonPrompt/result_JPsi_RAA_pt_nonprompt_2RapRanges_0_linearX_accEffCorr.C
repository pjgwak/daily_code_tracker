void result_JPsi_RAA_pt_nonprompt_2RapRanges_0_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__89 = new TH1F("haxes__89","haxes",1,0,30);
   haxes__89->SetMinimum(0);
   haxes__89->SetMaximum(1.5);
   haxes__89->SetLineStyle(0);
   haxes__89->SetMarkerStyle(20);
   haxes__89->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__89->GetXaxis()->CenterTitle(true);
   haxes__89->GetXaxis()->SetLabelFont(42);
   haxes__89->GetXaxis()->SetLabelOffset(0.007);
   haxes__89->GetXaxis()->SetLabelSize(0.05);
   haxes__89->GetXaxis()->SetTitleSize(0.06);
   haxes__89->GetXaxis()->SetTitleOffset(0.9);
   haxes__89->GetXaxis()->SetTitleFont(42);
   haxes__89->GetYaxis()->SetTitle("R_{AA}");
   haxes__89->GetYaxis()->SetLabelFont(42);
   haxes__89->GetYaxis()->SetLabelOffset(0.007);
   haxes__89->GetYaxis()->SetLabelSize(0.05);
   haxes__89->GetYaxis()->SetTitleSize(0.06);
   haxes__89->GetYaxis()->SetTitleOffset(1.25);
   haxes__89->GetYaxis()->SetTitleFont(42);
   haxes__89->GetZaxis()->SetLabelFont(42);
   haxes__89->GetZaxis()->SetLabelOffset(0.007);
   haxes__89->GetZaxis()->SetLabelSize(0.05);
   haxes__89->GetZaxis()->SetTitleSize(0.06);
   haxes__89->GetZaxis()->SetTitleFont(42);
   haxes__89->Draw("");
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
   
   Double_t bin_0_syst_fx3181[7] = {
   7.5,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_0_syst_fy3181[7] = {
   0.4758639,
   0.5337355,
   0.4322224,
   0.4661941,
   0.455847,
   0.4938115,
   0.4271174};
   Double_t bin_0_syst_felx3181[7] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fely3181[7] = {
   0.05689149,
   0.04293819,
   0.02777018,
   0.02439754,
   0.02199404,
   0.02462086,
   0.00849685};
   Double_t bin_0_syst_fehx3181[7] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t bin_0_syst_fehy3181[7] = {
   0.05689149,
   0.04293819,
   0.02777018,
   0.02439754,
   0.02199404,
   0.02462086,
   0.00849685};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(7,bin_0_syst_fx3181,bin_0_syst_fy3181,bin_0_syst_felx3181,bin_0_syst_fehx3181,bin_0_syst_fely3181,bin_0_syst_fehy3181);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 2057;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3181 = new TH1F("Graph_bin_0_syst3181","Graph",100,3.65,43.85);
   Graph_bin_0_syst3181->SetMinimum(0.3872301);
   Graph_bin_0_syst3181->SetMaximum(0.5938958);
   Graph_bin_0_syst3181->SetDirectory(0);
   Graph_bin_0_syst3181->SetStats(0);
   Graph_bin_0_syst3181->SetLineStyle(0);
   Graph_bin_0_syst3181->SetMarkerStyle(20);
   Graph_bin_0_syst3181->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3181->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3181->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3181->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3181->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3181->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3181->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3181->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3181->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3181->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3181->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3181->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3181->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3181->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3181->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3181->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3181->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3181);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3182[7] = {
   7.5,
   9,
   10.25,
   13,
   17.5,
   25,
   40};
   Double_t bin_0_fy3182[7] = {
   0.4758639,
   0.5337355,
   0.4322224,
   0.4661941,
   0.455847,
   0.4938115,
   0.4271174};
   Double_t bin_0_felx3182[7] = {
   1,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_0_fely3182[7] = {
   0.02661231,
   0.02779628,
   0.02014286,
   0.01599972,
   0.02244974,
   0.03362696,
   0.05959911};
   Double_t bin_0_fehx3182[7] = {
   1,
   0.5,
   0.75,
   2,
   2.5,
   5,
   10};
   Double_t bin_0_fehy3182[7] = {
   0.02661231,
   0.02779628,
   0.02014286,
   0.01599972,
   0.02244974,
   0.03362696,
   0.05959911};
   grae = new TGraphAsymmErrors(7,bin_0_fx3182,bin_0_fy3182,bin_0_felx3182,bin_0_fehx3182,bin_0_fely3182,bin_0_fehy3182);
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
   
   TH1F *Graph_bin_03182 = new TH1F("Graph_bin_03182","Graph",100,2.15,54.35);
   Graph_bin_03182->SetMinimum(0.3481169);
   Graph_bin_03182->SetMaximum(0.5809331);
   Graph_bin_03182->SetDirectory(0);
   Graph_bin_03182->SetStats(0);
   Graph_bin_03182->SetLineStyle(0);
   Graph_bin_03182->SetMarkerStyle(20);
   Graph_bin_03182->GetXaxis()->SetLabelFont(42);
   Graph_bin_03182->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03182->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03182->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03182->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03182->GetXaxis()->SetTitleFont(42);
   Graph_bin_03182->GetYaxis()->SetLabelFont(42);
   Graph_bin_03182->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03182->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03182->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03182->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03182->GetYaxis()->SetTitleFont(42);
   Graph_bin_03182->GetZaxis()->SetLabelFont(42);
   Graph_bin_03182->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03182->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03182->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03182->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03182);
   
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
   
   Double_t bin_3_syst_fx3183[11] = {
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
   Double_t bin_3_syst_fy3183[11] = {
   0.7389477,
   0.6345531,
   0.5167278,
   0.575668,
   0.4829829,
   0.3932381,
   0.4196652,
   0.425233,
   0.3568581,
   0.3908795,
   0.8555776};
   Double_t bin_3_syst_felx3183[11] = {
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
   Double_t bin_3_syst_fely3183[11] = {
   0.2613933,
   0.1109919,
   0.05985485,
   0.06045722,
   0.0547964,
   0.03650961,
   0.03325148,
   0.02828386,
   0.03287887,
   0.03599493,
   0.06545542};
   Double_t bin_3_syst_fehx3183[11] = {
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
   Double_t bin_3_syst_fehy3183[11] = {
   0.2613933,
   0.1109919,
   0.05985485,
   0.06045722,
   0.0547964,
   0.03650961,
   0.03325148,
   0.02828386,
   0.03287887,
   0.03599493,
   0.06545542};
   grae = new TGraphAsymmErrors(11,bin_3_syst_fx3183,bin_3_syst_fy3183,bin_3_syst_felx3183,bin_3_syst_fehx3183,bin_3_syst_fely3183,bin_3_syst_fehy3183);
   grae->SetName("bin_3_syst");
   grae->SetTitle("Graph");

   ci = 2061;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_3_syst3183 = new TH1F("Graph_bin_3_syst3183","Graph",100,0,44.225);
   Graph_bin_3_syst3183->SetMinimum(0.256343);
   Graph_bin_3_syst3183->SetMaximum(1.067977);
   Graph_bin_3_syst3183->SetDirectory(0);
   Graph_bin_3_syst3183->SetStats(0);
   Graph_bin_3_syst3183->SetLineStyle(0);
   Graph_bin_3_syst3183->SetMarkerStyle(20);
   Graph_bin_3_syst3183->GetXaxis()->SetLabelFont(42);
   Graph_bin_3_syst3183->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3183->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3183->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3183->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_3_syst3183->GetXaxis()->SetTitleFont(42);
   Graph_bin_3_syst3183->GetYaxis()->SetLabelFont(42);
   Graph_bin_3_syst3183->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3183->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3183->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3183->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_3_syst3183->GetYaxis()->SetTitleFont(42);
   Graph_bin_3_syst3183->GetZaxis()->SetLabelFont(42);
   Graph_bin_3_syst3183->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3183->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3183->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3183->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_3_syst3183);
   
   grae->Draw("5");
   
   Double_t bin_3_fx3184[11] = {
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
   Double_t bin_3_fy3184[11] = {
   0.7389477,
   0.6345531,
   0.5167278,
   0.575668,
   0.4829829,
   0.3932381,
   0.4196652,
   0.425233,
   0.3568581,
   0.3908795,
   0.8555776};
   Double_t bin_3_felx3184[11] = {
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
   Double_t bin_3_fely3184[11] = {
   0.0634266,
   0.04242643,
   0.0343517,
   0.03442388,
   0.0337569,
   0.02945255,
   0.02988956,
   0.0248935,
   0.0334153,
   0.0558244,
   0.2383666};
   Double_t bin_3_fehx3184[11] = {
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
   Double_t bin_3_fehy3184[11] = {
   0.0634266,
   0.04242643,
   0.0343517,
   0.03442388,
   0.0337569,
   0.02945255,
   0.02988956,
   0.0248935,
   0.0334153,
   0.0558244,
   0.2383666};
   grae = new TGraphAsymmErrors(11,bin_3_fx3184,bin_3_fy3184,bin_3_felx3184,bin_3_fehx3184,bin_3_fely3184,bin_3_fehy3184);
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
   
   TH1F *Graph_bin_33184 = new TH1F("Graph_bin_33184","Graph",100,0,54.7);
   Graph_bin_33184->SetMinimum(0.2463927);
   Graph_bin_33184->SetMaximum(1.170994);
   Graph_bin_33184->SetDirectory(0);
   Graph_bin_33184->SetStats(0);
   Graph_bin_33184->SetLineStyle(0);
   Graph_bin_33184->SetMarkerStyle(20);
   Graph_bin_33184->GetXaxis()->SetLabelFont(42);
   Graph_bin_33184->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_33184->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_33184->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_33184->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_33184->GetXaxis()->SetTitleFont(42);
   Graph_bin_33184->GetYaxis()->SetLabelFont(42);
   Graph_bin_33184->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_33184->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_33184->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_33184->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_33184->GetYaxis()->SetTitleFont(42);
   Graph_bin_33184->GetZaxis()->SetLabelFont(42);
   Graph_bin_33184->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_33184->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_33184->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_33184->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_33184);
   
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
   
   TH1F *haxes_copy__90 = new TH1F("haxes_copy__90","haxes",1,0,30);
   haxes_copy__90->SetMinimum(0);
   haxes_copy__90->SetMaximum(1.5);
   haxes_copy__90->SetDirectory(0);
   haxes_copy__90->SetLineStyle(0);
   haxes_copy__90->SetMarkerStyle(20);
   haxes_copy__90->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__90->GetXaxis()->CenterTitle(true);
   haxes_copy__90->GetXaxis()->SetLabelFont(42);
   haxes_copy__90->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__90->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__90->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__90->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__90->GetXaxis()->SetTitleFont(42);
   haxes_copy__90->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__90->GetYaxis()->SetLabelFont(42);
   haxes_copy__90->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__90->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__90->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__90->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__90->GetYaxis()->SetTitleFont(42);
   haxes_copy__90->GetZaxis()->SetLabelFont(42);
   haxes_copy__90->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__90->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__90->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__90->GetZaxis()->SetTitleFont(42);
   haxes_copy__90->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
