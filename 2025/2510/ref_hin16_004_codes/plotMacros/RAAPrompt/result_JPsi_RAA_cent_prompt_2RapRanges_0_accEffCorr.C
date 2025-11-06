void result_JPsi_RAA_cent_prompt_2RapRanges_0_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:52 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__45 = new TH1F("haxes__45","haxes",1,0,420);
   haxes__45->SetMinimum(0);
   haxes__45->SetMaximum(1.5);
   haxes__45->SetLineStyle(0);
   haxes__45->SetMarkerStyle(20);
   haxes__45->GetXaxis()->SetTitle("N_{part}");
   haxes__45->GetXaxis()->CenterTitle(true);
   haxes__45->GetXaxis()->SetLabelFont(42);
   haxes__45->GetXaxis()->SetLabelOffset(0.007);
   haxes__45->GetXaxis()->SetLabelSize(0.05);
   haxes__45->GetXaxis()->SetTitleSize(0.06);
   haxes__45->GetXaxis()->SetTitleOffset(0.9);
   haxes__45->GetXaxis()->SetTitleFont(42);
   haxes__45->GetYaxis()->SetTitle("R_{AA}");
   haxes__45->GetYaxis()->SetLabelFont(42);
   haxes__45->GetYaxis()->SetLabelOffset(0.007);
   haxes__45->GetYaxis()->SetLabelSize(0.05);
   haxes__45->GetYaxis()->SetTitleSize(0.06);
   haxes__45->GetYaxis()->SetTitleOffset(1.25);
   haxes__45->GetYaxis()->SetTitleFont(42);
   haxes__45->GetZaxis()->SetLabelFont(42);
   haxes__45->GetZaxis()->SetLabelOffset(0.007);
   haxes__45->GetZaxis()->SetLabelSize(0.05);
   haxes__45->GetZaxis()->SetTitleSize(0.06);
   haxes__45->GetZaxis()->SetTitleFont(42);
   haxes__45->Draw("");
   TBox *box = new TBox(400,0.946786,420,1.053214);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   box = new TBox(400,0.946786,420,1.053214);

   ci = TColor::GetColor("#9999ff");
   box->SetFillColor(ci);
   box->SetFillStyle(0);

   ci = TColor::GetColor("#000099");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t bin_0_syst_fx3091[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_syst_fy3091[6] = {
   0.2413696,
   0.331585,
   0.4115941,
   0.4559406,
   0.5993138,
   0.7357528};
   Double_t bin_0_syst_felx3091[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fely3091[6] = {
   0.01517889,
   0.02073442,
   0.02870561,
   0.03287619,
   0.05076787,
   0.06708877};
   Double_t bin_0_syst_fehx3091[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fehy3091[6] = {
   0.01410663,
   0.0193776,
   0.02763823,
   0.03226345,
   0.05170722,
   0.09168759};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3091,bin_0_syst_fy3091,bin_0_syst_felx3091,bin_0_syst_fehx3091,bin_0_syst_fely3091,bin_0_syst_fehy3091);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");

   ci = 2001;
   color = new TColor(ci, 0.6, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3091 = new TH1F("Graph_bin_0_syst3091","Graph",100,0,398.493);
   Graph_bin_0_syst3091->SetMinimum(0.1660658);
   Graph_bin_0_syst3091->SetMaximum(0.8875653);
   Graph_bin_0_syst3091->SetDirectory(0);
   Graph_bin_0_syst3091->SetStats(0);
   Graph_bin_0_syst3091->SetLineStyle(0);
   Graph_bin_0_syst3091->SetMarkerStyle(20);
   Graph_bin_0_syst3091->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3091->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3091->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3091->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3091->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3091->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3091->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3091->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3091->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3091->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3091->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3091->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3091->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3091->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3091->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3091->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3091->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3091);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3092[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_0_fy3092[6] = {
   0.2413696,
   0.331585,
   0.4115941,
   0.4559406,
   0.5993138,
   0.7357528};
   Double_t bin_0_felx3092[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_0_fely3092[6] = {
   0.007752963,
   0.01171436,
   0.01635506,
   0.02091978,
   0.03446636,
   0.04335165};
   Double_t bin_0_fehx3092[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_0_fehy3092[6] = {
   0.007752963,
   0.01171436,
   0.01635506,
   0.02091978,
   0.03446636,
   0.04335165};
   grae = new TGraphAsymmErrors(6,bin_0_fx3092,bin_0_fy3092,bin_0_felx3092,bin_0_fehx3092,bin_0_fely3092,bin_0_fehy3092);
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
   
   TH1F *Graph_bin_03092 = new TH1F("Graph_bin_03092","Graph",100,0,395.233);
   Graph_bin_03092->SetMinimum(0.1790679);
   Graph_bin_03092->SetMaximum(0.8336532);
   Graph_bin_03092->SetDirectory(0);
   Graph_bin_03092->SetStats(0);
   Graph_bin_03092->SetLineStyle(0);
   Graph_bin_03092->SetMarkerStyle(20);
   Graph_bin_03092->GetXaxis()->SetLabelFont(42);
   Graph_bin_03092->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03092->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03092->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03092->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03092->GetXaxis()->SetTitleFont(42);
   Graph_bin_03092->GetYaxis()->SetLabelFont(42);
   Graph_bin_03092->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03092->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03092->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03092->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03092->GetYaxis()->SetTitleFont(42);
   Graph_bin_03092->GetZaxis()->SetLabelFont(42);
   Graph_bin_03092->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03092->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03092->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03092->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03092);
   
   grae->Draw("p");
   box = new TBox(380,0.9518584,400,1.048142);

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
   
   Double_t bin_3_syst_fx3093[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_3_syst_fy3093[6] = {
   0.2123862,
   0.2700907,
   0.3489137,
   0.4844321,
   0.5143154,
   0.6722098};
   Double_t bin_3_syst_felx3093[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_3_syst_fely3093[6] = {
   0.01609166,
   0.02000518,
   0.02551192,
   0.03774194,
   0.04777974,
   0.06176561};
   Double_t bin_3_syst_fehx3093[6] = {
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_3_syst_fehy3093[6] = {
   0.01531758,
   0.01908134,
   0.02465036,
   0.03714023,
   0.04851593,
   0.08411424};
   grae = new TGraphAsymmErrors(6,bin_3_syst_fx3093,bin_3_syst_fy3093,bin_3_syst_felx3093,bin_3_syst_fehx3093,bin_3_syst_fely3093,bin_3_syst_fehy3093);
   grae->SetName("bin_3_syst");
   grae->SetTitle("Graph");

   ci = 2005;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_3_syst3093 = new TH1F("Graph_bin_3_syst3093","Graph",100,0,398.493);
   Graph_bin_3_syst3093->SetMinimum(0.1402916);
   Graph_bin_3_syst3093->SetMaximum(0.812327);
   Graph_bin_3_syst3093->SetDirectory(0);
   Graph_bin_3_syst3093->SetStats(0);
   Graph_bin_3_syst3093->SetLineStyle(0);
   Graph_bin_3_syst3093->SetMarkerStyle(20);
   Graph_bin_3_syst3093->GetXaxis()->SetLabelFont(42);
   Graph_bin_3_syst3093->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3093->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3093->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3093->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_3_syst3093->GetXaxis()->SetTitleFont(42);
   Graph_bin_3_syst3093->GetYaxis()->SetLabelFont(42);
   Graph_bin_3_syst3093->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3093->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3093->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3093->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_3_syst3093->GetYaxis()->SetTitleFont(42);
   Graph_bin_3_syst3093->GetZaxis()->SetLabelFont(42);
   Graph_bin_3_syst3093->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_3_syst3093->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_3_syst3093->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_3_syst3093->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_3_syst3093);
   
   grae->Draw("5");
   
   Double_t bin_3_fx3094[6] = {
   358.8,
   264.2,
   189.2,
   131.4,
   86.95,
   21.87};
   Double_t bin_3_fy3094[6] = {
   0.2123862,
   0.2700907,
   0.3489137,
   0.4844321,
   0.5143154,
   0.6722098};
   Double_t bin_3_felx3094[6] = {
   2.6,
   3.8,
   4.1,
   4,
   3.7,
   1};
   Double_t bin_3_fely3094[6] = {
   0.01163154,
   0.01463353,
   0.02026573,
   0.02851247,
   0.03879597,
   0.04866477};
   Double_t bin_3_fehx3094[6] = {
   2.4,
   3.6,
   4,
   4,
   3.7,
   1.8};
   Double_t bin_3_fehy3094[6] = {
   0.01163154,
   0.01463353,
   0.02026573,
   0.02851247,
   0.03879597,
   0.04866477};
   grae = new TGraphAsymmErrors(6,bin_3_fx3094,bin_3_fy3094,bin_3_felx3094,bin_3_fehx3094,bin_3_fely3094,bin_3_fehy3094);
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
   
   TH1F *Graph_bin_33094 = new TH1F("Graph_bin_33094","Graph",100,0,395.233);
   Graph_bin_33094->SetMinimum(0.1487426);
   Graph_bin_33094->SetMaximum(0.7728866);
   Graph_bin_33094->SetDirectory(0);
   Graph_bin_33094->SetStats(0);
   Graph_bin_33094->SetLineStyle(0);
   Graph_bin_33094->SetMarkerStyle(20);
   Graph_bin_33094->GetXaxis()->SetLabelFont(42);
   Graph_bin_33094->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_33094->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_33094->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_33094->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_33094->GetXaxis()->SetTitleFont(42);
   Graph_bin_33094->GetYaxis()->SetLabelFont(42);
   Graph_bin_33094->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_33094->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_33094->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_33094->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_33094->GetYaxis()->SetTitleFont(42);
   Graph_bin_33094->GetZaxis()->SetLabelFont(42);
   Graph_bin_33094->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_33094->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_33094->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_33094->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_33094);
   
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
   
   TH1F *haxes_copy__46 = new TH1F("haxes_copy__46","haxes",1,0,420);
   haxes_copy__46->SetMinimum(0);
   haxes_copy__46->SetMaximum(1.5);
   haxes_copy__46->SetDirectory(0);
   haxes_copy__46->SetLineStyle(0);
   haxes_copy__46->SetMarkerStyle(20);
   haxes_copy__46->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__46->GetXaxis()->CenterTitle(true);
   haxes_copy__46->GetXaxis()->SetLabelFont(42);
   haxes_copy__46->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__46->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__46->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__46->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__46->GetXaxis()->SetTitleFont(42);
   haxes_copy__46->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__46->GetYaxis()->SetLabelFont(42);
   haxes_copy__46->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__46->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__46->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__46->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__46->GetYaxis()->SetTitleFont(42);
   haxes_copy__46->GetZaxis()->SetLabelFont(42);
   haxes_copy__46->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__46->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__46->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__46->GetZaxis()->SetTitleFont(42);
   haxes_copy__46->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
