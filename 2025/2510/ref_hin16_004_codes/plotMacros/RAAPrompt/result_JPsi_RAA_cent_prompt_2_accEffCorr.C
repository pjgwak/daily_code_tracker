void result_JPsi_RAA_cent_prompt_2_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:50 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__11 = new TH1F("haxes__11","haxes",1,0,420);
   haxes__11->SetMinimum(0);
   haxes__11->SetMaximum(1.5);
   haxes__11->SetLineStyle(0);
   haxes__11->SetMarkerStyle(20);
   haxes__11->GetXaxis()->SetTitle("N_{part}");
   haxes__11->GetXaxis()->CenterTitle(true);
   haxes__11->GetXaxis()->SetLabelFont(42);
   haxes__11->GetXaxis()->SetLabelOffset(0.007);
   haxes__11->GetXaxis()->SetLabelSize(0.05);
   haxes__11->GetXaxis()->SetTitleSize(0.06);
   haxes__11->GetXaxis()->SetTitleOffset(0.9);
   haxes__11->GetXaxis()->SetTitleFont(42);
   haxes__11->GetYaxis()->SetTitle("R_{AA}");
   haxes__11->GetYaxis()->SetLabelFont(42);
   haxes__11->GetYaxis()->SetLabelOffset(0.007);
   haxes__11->GetYaxis()->SetLabelSize(0.05);
   haxes__11->GetYaxis()->SetTitleSize(0.06);
   haxes__11->GetYaxis()->SetTitleOffset(1.25);
   haxes__11->GetYaxis()->SetTitleFont(42);
   haxes__11->GetZaxis()->SetLabelFont(42);
   haxes__11->GetZaxis()->SetLabelOffset(0.007);
   haxes__11->GetZaxis()->SetLabelSize(0.05);
   haxes__11->GetZaxis()->SetTitleSize(0.06);
   haxes__11->GetZaxis()->SetTitleFont(42);
   haxes__11->Draw("");
   TBox *box = new TBox(400,0.9178617,420,1.082138);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(400,0.9178617,420,1.082138);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3027[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_syst_fy3027[6] = {
   0.3279326,
   0.4433551,
   0.5643141,
   0.5139786,
   0.565645,
   0.6967748};
   Double_t bin_0_syst_felx3027[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fely3027[6] = {
   0.05970064,
   0.07069683,
   0.08239869,
   0.06954239,
   0.07957294,
   0.09351972};
   Double_t bin_0_syst_fehx3027[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fehy3027[6] = {
   0.05921319,
   0.0700053,
   0.08170982,
   0.06917675,
   0.08010993,
   0.1106742};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3027,bin_0_syst_fy3027,bin_0_syst_felx3027,bin_0_syst_fehx3027,bin_0_syst_fely3027,bin_0_syst_fehy3027);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3027 = new TH1F("Graph_bin_0_syst3027","Graph",100,0,398.493);
   Graph_bin_0_syst3027->SetMinimum(0.2143102);
   Graph_bin_0_syst3027->SetMaximum(0.8613707);
   Graph_bin_0_syst3027->SetDirectory(0);
   Graph_bin_0_syst3027->SetStats(0);
   Graph_bin_0_syst3027->SetLineStyle(0);
   Graph_bin_0_syst3027->SetMarkerStyle(20);
   Graph_bin_0_syst3027->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3027->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3027->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3027->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3027->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3027->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3027->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3027->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3027->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3027->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3027->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3027->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3027->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3027->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3027->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3027->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3027->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3027);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3028[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_fy3028[6] = {
   0.3279326,
   0.4433551,
   0.5643141,
   0.5139786,
   0.565645,
   0.6967748};
   Double_t bin_0_felx3028[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_0_fely3028[6] = {
   0.03421559,
   0.03664711,
   0.05136604,
   0.04226348,
   0.05333248,
   0.06381975};
   Double_t bin_0_fehx3028[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_0_fehy3028[6] = {
   0.03421559,
   0.03664711,
   0.05136604,
   0.04226348,
   0.05333248,
   0.06381975};
   grae = new TGraphAsymmErrors(6,bin_0_fx3028,bin_0_fy3028,bin_0_felx3028,bin_0_fehx3028,bin_0_fely3028,bin_0_fehy3028);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(33);
   grae->SetMarkerSize(2.2);
   
   TH1F *Graph_bin_03028 = new TH1F("Graph_bin_03028","Graph",100,0,395.233);
   Graph_bin_03028->SetMinimum(0.2470292);
   Graph_bin_03028->SetMaximum(0.8072823);
   Graph_bin_03028->SetDirectory(0);
   Graph_bin_03028->SetStats(0);
   Graph_bin_03028->SetLineStyle(0);
   Graph_bin_03028->SetMarkerStyle(20);
   Graph_bin_03028->GetXaxis()->SetLabelFont(42);
   Graph_bin_03028->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03028->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03028->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03028->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03028->GetXaxis()->SetTitleFont(42);
   Graph_bin_03028->GetYaxis()->SetLabelFont(42);
   Graph_bin_03028->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03028->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03028->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03028->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03028->GetYaxis()->SetTitleFont(42);
   Graph_bin_03028->GetZaxis()->SetLabelFont(42);
   Graph_bin_03028->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03028->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03028->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03028->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03028);
   
   grae->Draw("p");
   box = new TBox(380,0.9518584,400,1.048142);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(380,0.9518584,400,1.048142);

   ci = TColor::GetColor("#99ff99");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#009900");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_1_syst_fx3029[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_1_syst_fy3029[6] = {
   0.2123862,
   0.2700907,
   0.3489137,
   0.4844321,
   0.5143154,
   0.6722098};
   Double_t bin_1_syst_felx3029[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_1_syst_fely3029[6] = {
   0.01609166,
   0.02000518,
   0.02551192,
   0.03774194,
   0.04777974,
   0.06176561};
   Double_t bin_1_syst_fehx3029[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_1_syst_fehy3029[6] = {
   0.01531758,
   0.01908134,
   0.02465036,
   0.03714023,
   0.04851593,
   0.08411424};
   grae = new TGraphAsymmErrors(6,bin_1_syst_fx3029,bin_1_syst_fy3029,bin_1_syst_felx3029,bin_1_syst_fehx3029,bin_1_syst_fely3029,bin_1_syst_fehy3029);
   grae->SetName("bin_1_syst");
   grae->SetTitle("Graph");

   ci = 1971;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_1_syst3029 = new TH1F("Graph_bin_1_syst3029","Graph",100,0,398.493);
   Graph_bin_1_syst3029->SetMinimum(0.1402916);
   Graph_bin_1_syst3029->SetMaximum(0.812327);
   Graph_bin_1_syst3029->SetDirectory(0);
   Graph_bin_1_syst3029->SetStats(0);
   Graph_bin_1_syst3029->SetLineStyle(0);
   Graph_bin_1_syst3029->SetMarkerStyle(20);
   Graph_bin_1_syst3029->GetXaxis()->SetLabelFont(42);
   Graph_bin_1_syst3029->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3029->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3029->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3029->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_1_syst3029->GetXaxis()->SetTitleFont(42);
   Graph_bin_1_syst3029->GetYaxis()->SetLabelFont(42);
   Graph_bin_1_syst3029->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3029->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3029->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3029->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_1_syst3029->GetYaxis()->SetTitleFont(42);
   Graph_bin_1_syst3029->GetZaxis()->SetLabelFont(42);
   Graph_bin_1_syst3029->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_1_syst3029->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_1_syst3029->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_1_syst3029->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_1_syst3029);
   
   grae->Draw("5");
   
   Double_t bin_1_fx3030[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_1_fy3030[6] = {
   0.2123862,
   0.2700907,
   0.3489137,
   0.4844321,
   0.5143154,
   0.6722098};
   Double_t bin_1_felx3030[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_1_fely3030[6] = {
   0.01163154,
   0.01463353,
   0.02026573,
   0.02851247,
   0.03879597,
   0.04866477};
   Double_t bin_1_fehx3030[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_1_fehy3030[6] = {
   0.01163154,
   0.01463353,
   0.02026573,
   0.02851247,
   0.03879597,
   0.04866477};
   grae = new TGraphAsymmErrors(6,bin_1_fx3030,bin_1_fy3030,bin_1_felx3030,bin_1_fehx3030,bin_1_fely3030,bin_1_fehy3030);
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
   
   TH1F *Graph_bin_13030 = new TH1F("Graph_bin_13030","Graph",100,0,395.233);
   Graph_bin_13030->SetMinimum(0.1487426);
   Graph_bin_13030->SetMaximum(0.7728866);
   Graph_bin_13030->SetDirectory(0);
   Graph_bin_13030->SetStats(0);
   Graph_bin_13030->SetLineStyle(0);
   Graph_bin_13030->SetMarkerStyle(20);
   Graph_bin_13030->GetXaxis()->SetLabelFont(42);
   Graph_bin_13030->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_13030->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_13030->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_13030->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_13030->GetXaxis()->SetTitleFont(42);
   Graph_bin_13030->GetYaxis()->SetLabelFont(42);
   Graph_bin_13030->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_13030->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_13030->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_13030->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_13030->GetYaxis()->SetTitleFont(42);
   Graph_bin_13030->GetZaxis()->SetLabelFont(42);
   Graph_bin_13030->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_13030->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_13030->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_13030->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_13030);
   
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
      tex = new TLatex(0.2,0.85,"Prompt J/#psi");
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
   
   TH1F *haxes_copy__12 = new TH1F("haxes_copy__12","haxes",1,0,420);
   haxes_copy__12->SetMinimum(0);
   haxes_copy__12->SetMaximum(1.5);
   haxes_copy__12->SetDirectory(0);
   haxes_copy__12->SetLineStyle(0);
   haxes_copy__12->SetMarkerStyle(20);
   haxes_copy__12->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__12->GetXaxis()->CenterTitle(true);
   haxes_copy__12->GetXaxis()->SetLabelFont(42);
   haxes_copy__12->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__12->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__12->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__12->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__12->GetXaxis()->SetTitleFont(42);
   haxes_copy__12->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__12->GetYaxis()->SetLabelFont(42);
   haxes_copy__12->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__12->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__12->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__12->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__12->GetYaxis()->SetTitleFont(42);
   haxes_copy__12->GetZaxis()->SetLabelFont(42);
   haxes_copy__12->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__12->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__12->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__12->GetZaxis()->SetTitleFont(42);
   haxes_copy__12->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
