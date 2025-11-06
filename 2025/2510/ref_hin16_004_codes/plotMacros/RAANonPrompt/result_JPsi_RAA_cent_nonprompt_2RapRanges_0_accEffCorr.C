void result_JPsi_RAA_cent_nonprompt_2RapRanges_0_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:55 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__93 = new TH1F("haxes__93","haxes",1,0,420);
   haxes__93->SetMinimum(0);
   haxes__93->SetMaximum(1.5);
   haxes__93->SetLineStyle(0);
   haxes__93->SetMarkerStyle(20);
   haxes__93->GetXaxis()->SetTitle("N_{part}");
   haxes__93->GetXaxis()->CenterTitle(true);
   haxes__93->GetXaxis()->SetLabelFont(42);
   haxes__93->GetXaxis()->SetLabelOffset(0.007);
   haxes__93->GetXaxis()->SetLabelSize(0.05);
   haxes__93->GetXaxis()->SetTitleSize(0.06);
   haxes__93->GetXaxis()->SetTitleOffset(0.9);
   haxes__93->GetXaxis()->SetTitleFont(42);
   haxes__93->GetYaxis()->SetTitle("R_{AA}");
   haxes__93->GetYaxis()->SetLabelFont(42);
   haxes__93->GetYaxis()->SetLabelOffset(0.007);
   haxes__93->GetYaxis()->SetLabelSize(0.05);
   haxes__93->GetYaxis()->SetTitleSize(0.06);
   haxes__93->GetYaxis()->SetTitleOffset(1.25);
   haxes__93->GetYaxis()->SetTitleFont(42);
   haxes__93->GetZaxis()->SetLabelFont(42);
   haxes__93->GetZaxis()->SetLabelOffset(0.007);
   haxes__93->GetZaxis()->SetLabelSize(0.05);
   haxes__93->GetZaxis()->SetTitleSize(0.06);
   haxes__93->GetZaxis()->SetTitleFont(42);
   haxes__93->Draw("");
   TBox *box = new TBox(400,0.9510149,420,1.048985);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(400,0.9510149,420,1.048985);

   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_0_syst_fx3189[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_syst_fy3189[6] = {
   0.3564768,
   0.4559871,
   0.5305988,
   0.5375295,
   0.6208261,
   0.7049625};
   Double_t bin_0_syst_felx3189[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fely3189[6] = {
   0.02182821,
   0.02736947,
   0.037082,
   0.03780665,
   0.05113049,
   0.06387481};
   Double_t bin_0_syst_fehx3189[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fehy3189[6] = {
   0.02019844,
   0.02541978,
   0.03570896,
   0.03706571,
   0.0521308,
   0.08755367};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3189,bin_0_syst_fy3189,bin_0_syst_felx3189,bin_0_syst_fehx3189,bin_0_syst_fely3189,bin_0_syst_fehy3189);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 2063;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3189 = new TH1F("Graph_bin_0_syst3189","Graph",100,0,398.493);
   Graph_bin_0_syst3189->SetMinimum(0.2888618);
   Graph_bin_0_syst3189->SetMaximum(0.838303);
   Graph_bin_0_syst3189->SetDirectory(0);
   Graph_bin_0_syst3189->SetStats(0);
   Graph_bin_0_syst3189->SetLineStyle(0);
   Graph_bin_0_syst3189->SetMarkerStyle(20);
   Graph_bin_0_syst3189->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3189->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3189->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3189->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3189->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3189->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3189->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3189->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3189->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3189->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3189->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3189->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3189->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3189->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3189->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3189->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3189->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3189);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3190[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_fy3190[6] = {
   0.3564768,
   0.4559871,
   0.5305988,
   0.5375295,
   0.6208261,
   0.7049625};
   Double_t bin_0_felx3190[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_0_fely3190[6] = {
   0.01138299,
   0.01608901,
   0.0212067,
   0.02495564,
   0.03612058,
   0.04212051};
   Double_t bin_0_fehx3190[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_0_fehy3190[6] = {
   0.01138299,
   0.01608901,
   0.0212067,
   0.02495564,
   0.03612058,
   0.04212051};
   grae = new TGraphAsymmErrors(6,bin_0_fx3190,bin_0_fy3190,bin_0_felx3190,bin_0_fehx3190,bin_0_fely3190,bin_0_fehy3190);
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
   
   TH1F *Graph_bin_03190 = new TH1F("Graph_bin_03190","Graph",100,0,395.233);
   Graph_bin_03190->SetMinimum(0.3048949);
   Graph_bin_03190->SetMaximum(0.7872819);
   Graph_bin_03190->SetDirectory(0);
   Graph_bin_03190->SetStats(0);
   Graph_bin_03190->SetLineStyle(0);
   Graph_bin_03190->SetMarkerStyle(20);
   Graph_bin_03190->GetXaxis()->SetLabelFont(42);
   Graph_bin_03190->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03190->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03190->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03190->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03190->GetXaxis()->SetTitleFont(42);
   Graph_bin_03190->GetYaxis()->SetLabelFont(42);
   Graph_bin_03190->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03190->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03190->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03190->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03190->GetYaxis()->SetTitleFont(42);
   Graph_bin_03190->GetZaxis()->SetLabelFont(42);
   Graph_bin_03190->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03190->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03190->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03190->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03190);
   
   grae->Draw("p");
   box = new TBox(380,0.947524,400,1.052476);

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
   
   Double_t bin_3_syst_fx3191[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_3_syst_fy3191[6] = {
   0.3578631,
   0.4240683,
   0.4896907,
   0.6003001,
   0.6365092,
   0.7365134};
   Double_t bin_3_syst_felx3191[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_3_syst_fely3191[6] = {
   0.03023398,
   0.03548467,
   0.04020204,
   0.04933663,
   0.07203595,
   0.07008365};
   Double_t bin_3_syst_fehx3191[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_3_syst_fehy3191[6] = {
   0.02907002,
   0.03420736,
   0.03912897,
   0.04863038,
   0.07278569,
   0.09394418};
   grae = new TGraphAsymmErrors(6,bin_3_syst_fx3191,bin_3_syst_fy3191,bin_3_syst_felx3191,bin_3_syst_fehx3191,bin_3_syst_fely3191,bin_3_syst_fehy3191);
   grae->SetName("bin_3_syst");
   grae->SetTitle("Graph");

   ci = 2067;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_3_syst3191 = new TH1F("Graph_bin_3_syst3191","Graph",100,0,398.493);
   Graph_bin_3_syst3191->SetMinimum(0.2773463);
   Graph_bin_3_syst3191->SetMaximum(0.8807405);
   Graph_bin_3_syst3191->SetDirectory(0);
   Graph_bin_3_syst3191->SetStats(0);
   Graph_bin_3_syst3191->SetLineStyle(0);
   Graph_bin_3_syst3191->SetMarkerStyle(20);
   Graph_bin_3_syst3191->GetXaxis()->SetLabelFont(42);
   Graph_bin_3_syst3191->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3191->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3191->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3191->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_3_syst3191->GetXaxis()->SetTitleFont(42);
   Graph_bin_3_syst3191->GetYaxis()->SetLabelFont(42);
   Graph_bin_3_syst3191->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3191->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3191->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3191->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_3_syst3191->GetYaxis()->SetTitleFont(42);
   Graph_bin_3_syst3191->GetZaxis()->SetLabelFont(42);
   Graph_bin_3_syst3191->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3191->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3191->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3191->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_3_syst3191);
   
   grae->Draw("5");
   
   Double_t bin_3_fx3192[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_3_fy3192[6] = {
   0.3578631,
   0.4240683,
   0.4896907,
   0.6003001,
   0.6365092,
   0.7365134};
   Double_t bin_3_felx3192[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_3_fely3192[6] = {
   0.01896369,
   0.02257364,
   0.0282742,
   0.03542143,
   0.04849357,
   0.05466269};
   Double_t bin_3_fehx3192[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_3_fehy3192[6] = {
   0.01896369,
   0.02257364,
   0.0282742,
   0.03542143,
   0.04849357,
   0.05466269};
   grae = new TGraphAsymmErrors(6,bin_3_fx3192,bin_3_fy3192,bin_3_felx3192,bin_3_fehx3192,bin_3_fely3192,bin_3_fehy3192);
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
   
   TH1F *Graph_bin_33192 = new TH1F("Graph_bin_33192","Graph",100,0,395.233);
   Graph_bin_33192->SetMinimum(0.2936718);
   Graph_bin_33192->SetMaximum(0.8364038);
   Graph_bin_33192->SetDirectory(0);
   Graph_bin_33192->SetStats(0);
   Graph_bin_33192->SetLineStyle(0);
   Graph_bin_33192->SetMarkerStyle(20);
   Graph_bin_33192->GetXaxis()->SetLabelFont(42);
   Graph_bin_33192->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_33192->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_33192->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_33192->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_33192->GetXaxis()->SetTitleFont(42);
   Graph_bin_33192->GetYaxis()->SetLabelFont(42);
   Graph_bin_33192->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_33192->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_33192->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_33192->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_33192->GetYaxis()->SetTitleFont(42);
   Graph_bin_33192->GetZaxis()->SetLabelFont(42);
   Graph_bin_33192->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_33192->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_33192->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_33192->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_33192);
   
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
   TLine *line = new TLine(0,1,420,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
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
   
   TH1F *haxes_copy__94 = new TH1F("haxes_copy__94","haxes",1,0,420);
   haxes_copy__94->SetMinimum(0);
   haxes_copy__94->SetMaximum(1.5);
   haxes_copy__94->SetDirectory(0);
   haxes_copy__94->SetLineStyle(0);
   haxes_copy__94->SetMarkerStyle(20);
   haxes_copy__94->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__94->GetXaxis()->CenterTitle(true);
   haxes_copy__94->GetXaxis()->SetLabelFont(42);
   haxes_copy__94->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__94->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__94->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__94->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__94->GetXaxis()->SetTitleFont(42);
   haxes_copy__94->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__94->GetYaxis()->SetLabelFont(42);
   haxes_copy__94->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__94->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__94->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__94->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__94->GetYaxis()->SetTitleFont(42);
   haxes_copy__94->GetZaxis()->SetLabelFont(42);
   haxes_copy__94->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__94->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__94->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__94->GetZaxis()->SetTitleFont(42);
   haxes_copy__94->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
