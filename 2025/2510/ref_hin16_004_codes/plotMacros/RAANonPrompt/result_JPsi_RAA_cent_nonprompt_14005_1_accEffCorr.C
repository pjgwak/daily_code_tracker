void result_JPsi_RAA_cent_nonprompt_14005_1_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__85 = new TH1F("haxes__85","haxes",1,0,420);
   haxes__85->SetMinimum(0);
   haxes__85->SetMaximum(1.5);
   haxes__85->SetLineStyle(0);
   haxes__85->SetMarkerStyle(20);
   haxes__85->GetXaxis()->SetTitle("N_{part}");
   haxes__85->GetXaxis()->CenterTitle(true);
   haxes__85->GetXaxis()->SetLabelFont(42);
   haxes__85->GetXaxis()->SetLabelOffset(0.007);
   haxes__85->GetXaxis()->SetLabelSize(0.05);
   haxes__85->GetXaxis()->SetTitleSize(0.06);
   haxes__85->GetXaxis()->SetTitleOffset(0.9);
   haxes__85->GetXaxis()->SetTitleFont(42);
   haxes__85->GetYaxis()->SetTitle("R_{AA}");
   haxes__85->GetYaxis()->SetLabelFont(42);
   haxes__85->GetYaxis()->SetLabelOffset(0.007);
   haxes__85->GetYaxis()->SetLabelSize(0.05);
   haxes__85->GetYaxis()->SetTitleSize(0.06);
   haxes__85->GetYaxis()->SetTitleOffset(1.25);
   haxes__85->GetYaxis()->SetTitleFont(42);
   haxes__85->GetZaxis()->SetLabelFont(42);
   haxes__85->GetZaxis()->SetLabelOffset(0.007);
   haxes__85->GetZaxis()->SetLabelSize(0.05);
   haxes__85->GetZaxis()->SetTitleSize(0.06);
   haxes__85->GetZaxis()->SetTitleFont(42);
   haxes__85->Draw("");
   TBox *box = new TBox(400,0.9417953,420,1.058205);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(400,0.9417953,420,1.058205);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3177[16] = {
   384.3,
   358.8,
   333.3,
   285.4,
   242.9,
   205.7,
   226.7,
   172.7,
   144.1,
   118.7,
   96.51,
   77.39,
   53.86,
   30.57,
   46.81,
   8.297};
   Double_t bin_0_syst_fy3177[16] = {
   0.3648739,
   -0.349557,
   0.3753215,
   0.438496,
   0.4549329,
   0.4958641,
   -0.4559476,
   0.5161031,
   0.5714465,
   0.5965705,
   0.6006051,
   0.6831161,
   0.7207478,
   0.762282,
   -0.6091155,
   0.6150468};
   Double_t bin_0_syst_felx3177[16] = {
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fely3177[16] = {
   0.02428601,
   0.01048671,
   0.02559262,
   0.02827473,
   0.03146923,
   0.03567581,
   0.01550222,
   0.03864165,
   0.04188248,
   0.04640663,
   0.05251743,
   0.06232247,
   0.07357787,
   0.09615804,
   0.02680108,
   0.07264768};
   Double_t bin_0_syst_fehx3177[16] = {
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5,
   5};
   Double_t bin_0_syst_fehy3177[16] = {
   0.02282499,
   0.006641584,
   0.02406357,
   0.02660345,
   0.02967178,
   0.03567581,
   0.01185464,
   0.03779485,
   0.04151834,
   0.04555522,
   0.05293817,
   0.06339152,
   0.07607711,
   0.1043646,
   0.04081074,
   0.1073141};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(16,bin_0_syst_fx3177,bin_0_syst_fy3177,bin_0_syst_felx3177,bin_0_syst_fehx3177,bin_0_syst_fely3177,bin_0_syst_fehy3177);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3177 = new TH1F("Graph_bin_0_syst3177","Graph",100,0,427.9003);
   Graph_bin_0_syst3177->SetMinimum(-0.7861729);
   Graph_bin_0_syst3177->SetMaximum(1.016903);
   Graph_bin_0_syst3177->SetDirectory(0);
   Graph_bin_0_syst3177->SetStats(0);
   Graph_bin_0_syst3177->SetLineStyle(0);
   Graph_bin_0_syst3177->SetMarkerStyle(20);
   Graph_bin_0_syst3177->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3177->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3177->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3177->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3177->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3177->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3177->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3177->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3177->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3177->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3177->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3177->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3177->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3177->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3177->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3177->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3177->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3177);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3178[16] = {
   384.3,
   358.8,
   333.3,
   285.4,
   242.9,
   205.7,
   226.7,
   172.7,
   144.1,
   118.7,
   96.51,
   77.39,
   53.86,
   30.57,
   46.81,
   8.297};
   Double_t bin_0_fy3178[16] = {
   0.3648739,
   -0.349557,
   0.3753215,
   0.438496,
   0.4549329,
   0.4958641,
   -0.4559476,
   0.5161031,
   0.5714465,
   0.5965705,
   0.6006051,
   0.6831161,
   0.7207478,
   0.762282,
   -0.6091155,
   0.6150468};
   Double_t bin_0_felx3178[16] = {
   2,
   2.6,
   3.2,
   3.7,
   3.9,
   4.1,
   3.9,
   4,
   4,
   4,
   3.8,
   3.6,
   3.1,
   2.4,
   1.2,
   0.6};
   Double_t bin_0_fely3178[16] = {
   0.008959297,
   0.006499447,
   0.009540147,
   0.01161586,
   0.01257871,
   0.01465464,
   0.006585959,
   0.01661951,
   0.01800687,
   0.0211116,
   0.02454827,
   0.03023048,
   0.0287454,
   0.04349814,
   0.009679101,
   0.05221886};
   Double_t bin_0_fehx3178[16] = {
   1.8,
   2.4,
   3,
   3.5,
   3.8,
   3.9,
   3.7,
   4,
   4,
   4,
   3.8,
   3.7,
   3.2,
   2.6,
   2.4,
   1};
   Double_t bin_0_fehy3178[16] = {
   0.008959297,
   0.006499447,
   0.009540147,
   0.01161586,
   0.01257871,
   0.01465464,
   0.006585959,
   0.01661951,
   0.01800687,
   0.0211116,
   0.02454827,
   0.03023048,
   0.0287454,
   0.04349814,
   0.009679101,
   0.05221886};
   grae = new TGraphAsymmErrors(16,bin_0_fx3178,bin_0_fy3178,bin_0_felx3178,bin_0_fehx3178,bin_0_fely3178,bin_0_fehy3178);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03178 = new TH1F("Graph_bin_03178","Graph",100,0,423.9403);
   Graph_bin_03178->SetMinimum(-0.7612521);
   Graph_bin_03178->SetMaximum(0.9482376);
   Graph_bin_03178->SetDirectory(0);
   Graph_bin_03178->SetStats(0);
   Graph_bin_03178->SetLineStyle(0);
   Graph_bin_03178->SetMarkerStyle(20);
   Graph_bin_03178->GetXaxis()->SetLabelFont(42);
   Graph_bin_03178->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03178->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03178->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03178->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03178->GetXaxis()->SetTitleFont(42);
   Graph_bin_03178->GetYaxis()->SetLabelFont(42);
   Graph_bin_03178->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03178->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03178->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03178->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03178->GetYaxis()->SetTitleFont(42);
   Graph_bin_03178->GetZaxis()->SetLabelFont(42);
   Graph_bin_03178->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03178->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03178->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03178->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03178);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.51,0.47,0.83,0.62,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("bin_0","#sqrt[]{s_{NN}} = 5.02 TeV","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(28);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("NULL","#sqrt[]{s_{NN}} = 2.76 TeV","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLine *line = new TLine(0,1,420,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.044);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.73,"6.5 < p_{T} < 50 GeV/c");
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
   
   TH1F *haxes_copy__86 = new TH1F("haxes_copy__86","haxes",1,0,420);
   haxes_copy__86->SetMinimum(0);
   haxes_copy__86->SetMaximum(1.5);
   haxes_copy__86->SetDirectory(0);
   haxes_copy__86->SetLineStyle(0);
   haxes_copy__86->SetMarkerStyle(20);
   haxes_copy__86->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__86->GetXaxis()->CenterTitle(true);
   haxes_copy__86->GetXaxis()->SetLabelFont(42);
   haxes_copy__86->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__86->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__86->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__86->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__86->GetXaxis()->SetTitleFont(42);
   haxes_copy__86->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__86->GetYaxis()->SetLabelFont(42);
   haxes_copy__86->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__86->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__86->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__86->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__86->GetYaxis()->SetTitleFont(42);
   haxes_copy__86->GetZaxis()->SetLabelFont(42);
   haxes_copy__86->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__86->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__86->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__86->GetZaxis()->SetTitleFont(42);
   haxes_copy__86->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
