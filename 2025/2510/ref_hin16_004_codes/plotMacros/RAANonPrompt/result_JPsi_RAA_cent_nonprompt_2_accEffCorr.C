void result_JPsi_RAA_cent_nonprompt_2_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:53 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-84,-0.225,441,1.65);
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
   
   TH1F *haxes__59 = new TH1F("haxes__59","haxes",1,0,420);
   haxes__59->SetMinimum(0);
   haxes__59->SetMaximum(1.5);
   haxes__59->SetLineStyle(0);
   haxes__59->SetMarkerStyle(20);
   haxes__59->GetXaxis()->SetTitle("N_{part}");
   haxes__59->GetXaxis()->CenterTitle(true);
   haxes__59->GetXaxis()->SetLabelFont(42);
   haxes__59->GetXaxis()->SetLabelOffset(0.007);
   haxes__59->GetXaxis()->SetLabelSize(0.05);
   haxes__59->GetXaxis()->SetTitleSize(0.06);
   haxes__59->GetXaxis()->SetTitleOffset(0.9);
   haxes__59->GetXaxis()->SetTitleFont(42);
   haxes__59->GetYaxis()->SetTitle("R_{AA}");
   haxes__59->GetYaxis()->SetLabelFont(42);
   haxes__59->GetYaxis()->SetLabelOffset(0.007);
   haxes__59->GetYaxis()->SetLabelSize(0.05);
   haxes__59->GetYaxis()->SetTitleSize(0.06);
   haxes__59->GetYaxis()->SetTitleOffset(1.25);
   haxes__59->GetYaxis()->SetTitleFont(42);
   haxes__59->GetZaxis()->SetLabelFont(42);
   haxes__59->GetZaxis()->SetLabelOffset(0.007);
   haxes__59->GetZaxis()->SetLabelSize(0.05);
   haxes__59->GetZaxis()->SetTitleSize(0.06);
   haxes__59->GetZaxis()->SetTitleFont(42);
   haxes__59->Draw("");
   TBox *box = new TBox(400,0.911473,420,1.088527);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(400,0.911473,420,1.088527);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3125[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_syst_fy3125[6] = {
   0.5595308,
   0.6380381,
   0.6258448,
   0.8091553,
   0.8987208,
   0.8938506};
   Double_t bin_0_syst_felx3125[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fely3125[6] = {
   0.1701887,
   0.1290773,
   0.1255065,
   0.1328108,
   0.135971,
   0.1181526};
   Double_t bin_0_syst_fehx3125[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fehy3125[6] = {
   0.1696922,
   0.1282943,
   0.1249513,
   0.1323367,
   0.1367647,
   0.1404443};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3125,bin_0_syst_fy3125,bin_0_syst_felx3125,bin_0_syst_fehx3125,bin_0_syst_fely3125,bin_0_syst_fehy3125);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3125 = new TH1F("Graph_bin_0_syst3125","Graph",100,0,398.493);
   Graph_bin_0_syst3125->SetMinimum(0.3247277);
   Graph_bin_0_syst3125->SetMaximum(1.1001);
   Graph_bin_0_syst3125->SetDirectory(0);
   Graph_bin_0_syst3125->SetStats(0);
   Graph_bin_0_syst3125->SetLineStyle(0);
   Graph_bin_0_syst3125->SetMarkerStyle(20);
   Graph_bin_0_syst3125->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3125->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3125->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3125->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3125->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3125->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3125->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3125->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3125->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3125->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3125->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3125->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3125->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3125->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3125->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3125->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3125->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3125);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3126[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_fy3126[6] = {
   0.5595308,
   0.6380381,
   0.6258448,
   0.8091553,
   0.8987208,
   0.8938506};
   Double_t bin_0_felx3126[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_0_fely3126[6] = {
   0.05811291,
   0.0524912,
   0.05694825,
   0.06658004,
   0.08426589,
   0.08178288};
   Double_t bin_0_fehx3126[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_0_fehy3126[6] = {
   0.05811291,
   0.0524912,
   0.05694825,
   0.06658004,
   0.08426589,
   0.08178288};
   grae = new TGraphAsymmErrors(6,bin_0_fx3126,bin_0_fy3126,bin_0_felx3126,bin_0_fehx3126,bin_0_fely3126,bin_0_fehy3126);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(2.2);
   
   TH1F *Graph_bin_03126 = new TH1F("Graph_bin_03126","Graph",100,0,395.233);
   Graph_bin_03126->SetMinimum(0.453261);
   Graph_bin_03126->SetMaximum(1.031144);
   Graph_bin_03126->SetDirectory(0);
   Graph_bin_03126->SetStats(0);
   Graph_bin_03126->SetLineStyle(0);
   Graph_bin_03126->SetMarkerStyle(20);
   Graph_bin_03126->GetXaxis()->SetLabelFont(42);
   Graph_bin_03126->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03126->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03126->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03126->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03126->GetXaxis()->SetTitleFont(42);
   Graph_bin_03126->GetYaxis()->SetLabelFont(42);
   Graph_bin_03126->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03126->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03126->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03126->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03126->GetYaxis()->SetTitleFont(42);
   Graph_bin_03126->GetZaxis()->SetLabelFont(42);
   Graph_bin_03126->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03126->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03126->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03126->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03126);
   
   grae->Draw("p");
   box = new TBox(380,0.947524,400,1.052476);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(380,0.947524,400,1.052476);

   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_1_syst_fx3127[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_1_syst_fy3127[6] = {
   0.3578631,
   0.4240683,
   0.4896907,
   0.6003001,
   0.6365092,
   0.7365134};
   Double_t bin_1_syst_felx3127[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_1_syst_fely3127[6] = {
   0.03023398,
   0.03548467,
   0.04020204,
   0.04933663,
   0.07203595,
   0.07008365};
   Double_t bin_1_syst_fehx3127[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_1_syst_fehy3127[6] = {
   0.02907002,
   0.03420736,
   0.03912897,
   0.04863038,
   0.07278569,
   0.09394418};
   grae = new TGraphAsymmErrors(6,bin_1_syst_fx3127,bin_1_syst_fy3127,bin_1_syst_felx3127,bin_1_syst_fehx3127,bin_1_syst_fely3127,bin_1_syst_fehy3127);
   grae->SetName("bin_1_syst");
   grae->SetTitle("Graph");

   ci = 2033;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_1_syst3127 = new TH1F("Graph_bin_1_syst3127","Graph",100,0,398.493);
   Graph_bin_1_syst3127->SetMinimum(0.2773463);
   Graph_bin_1_syst3127->SetMaximum(0.8807405);
   Graph_bin_1_syst3127->SetDirectory(0);
   Graph_bin_1_syst3127->SetStats(0);
   Graph_bin_1_syst3127->SetLineStyle(0);
   Graph_bin_1_syst3127->SetMarkerStyle(20);
   Graph_bin_1_syst3127->GetXaxis()->SetLabelFont(42);
   Graph_bin_1_syst3127->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3127->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3127->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3127->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_1_syst3127->GetXaxis()->SetTitleFont(42);
   Graph_bin_1_syst3127->GetYaxis()->SetLabelFont(42);
   Graph_bin_1_syst3127->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3127->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3127->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3127->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_1_syst3127->GetYaxis()->SetTitleFont(42);
   Graph_bin_1_syst3127->GetZaxis()->SetLabelFont(42);
   Graph_bin_1_syst3127->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3127->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3127->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3127->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_1_syst3127);
   
   grae->Draw("5");
   
   Double_t bin_1_fx3128[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_1_fy3128[6] = {
   0.3578631,
   0.4240683,
   0.4896907,
   0.6003001,
   0.6365092,
   0.7365134};
   Double_t bin_1_felx3128[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_1_fely3128[6] = {
   0.01896369,
   0.02257364,
   0.0282742,
   0.03542143,
   0.04849357,
   0.05466269};
   Double_t bin_1_fehx3128[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_1_fehy3128[6] = {
   0.01896369,
   0.02257364,
   0.0282742,
   0.03542143,
   0.04849357,
   0.05466269};
   grae = new TGraphAsymmErrors(6,bin_1_fx3128,bin_1_fy3128,bin_1_felx3128,bin_1_fehx3128,bin_1_fely3128,bin_1_fehy3128);
   grae->SetName("bin_1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#009900");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_13128 = new TH1F("Graph_bin_13128","Graph",100,0,395.233);
   Graph_bin_13128->SetMinimum(0.2936718);
   Graph_bin_13128->SetMaximum(0.8364038);
   Graph_bin_13128->SetDirectory(0);
   Graph_bin_13128->SetStats(0);
   Graph_bin_13128->SetLineStyle(0);
   Graph_bin_13128->SetMarkerStyle(20);
   Graph_bin_13128->GetXaxis()->SetLabelFont(42);
   Graph_bin_13128->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_13128->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_13128->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_13128->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_13128->GetXaxis()->SetTitleFont(42);
   Graph_bin_13128->GetYaxis()->SetLabelFont(42);
   Graph_bin_13128->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_13128->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_13128->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_13128->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_13128->GetYaxis()->SetTitleFont(42);
   Graph_bin_13128->GetZaxis()->SetLabelFont(42);
   Graph_bin_13128->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_13128->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_13128->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_13128->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_13128);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.19,0.16,0.51,0.28,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("bin_0","3 < p_{T} < 6.5 GeV/c","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(33);
   entry->SetMarkerSize(2.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("bin_1","6.5 < p_{T} < 50 GeV/c","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,420,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"1.8 < |y| < 2.4");
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
      tex = new TLatex(0.96,0.9424,"PbPb 368 (<30%) / 464 (>30%) #mub^{-1}, pp 28.0 pb^{-1} (5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.03333333);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9424,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03333333);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.914,0.874,"CMS");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH1F *haxes_copy__60 = new TH1F("haxes_copy__60","haxes",1,0,420);
   haxes_copy__60->SetMinimum(0);
   haxes_copy__60->SetMaximum(1.5);
   haxes_copy__60->SetDirectory(0);
   haxes_copy__60->SetLineStyle(0);
   haxes_copy__60->SetMarkerStyle(20);
   haxes_copy__60->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__60->GetXaxis()->CenterTitle(true);
   haxes_copy__60->GetXaxis()->SetLabelFont(42);
   haxes_copy__60->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__60->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__60->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__60->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__60->GetXaxis()->SetTitleFont(42);
   haxes_copy__60->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__60->GetYaxis()->SetLabelFont(42);
   haxes_copy__60->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__60->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__60->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__60->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__60->GetYaxis()->SetTitleFont(42);
   haxes_copy__60->GetZaxis()->SetLabelFont(42);
   haxes_copy__60->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__60->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__60->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__60->GetZaxis()->SetTitleFont(42);
   haxes_copy__60->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
