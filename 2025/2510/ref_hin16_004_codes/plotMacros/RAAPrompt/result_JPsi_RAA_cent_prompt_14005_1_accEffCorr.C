void result_JPsi_RAA_cent_prompt_14005_1_accEffCorr()
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
   
   TH1F *haxes__37 = new TH1F("haxes__37","haxes",1,0,420);
   haxes__37->SetMinimum(0);
   haxes__37->SetMaximum(1.5);
   haxes__37->SetLineStyle(0);
   haxes__37->SetMarkerStyle(20);
   haxes__37->GetXaxis()->SetTitle("N_{part}");
   haxes__37->GetXaxis()->CenterTitle(true);
   haxes__37->GetXaxis()->SetLabelFont(42);
   haxes__37->GetXaxis()->SetLabelOffset(0.007);
   haxes__37->GetXaxis()->SetLabelSize(0.05);
   haxes__37->GetXaxis()->SetTitleSize(0.06);
   haxes__37->GetXaxis()->SetTitleOffset(0.9);
   haxes__37->GetXaxis()->SetTitleFont(42);
   haxes__37->GetYaxis()->SetTitle("R_{AA}");
   haxes__37->GetYaxis()->SetLabelFont(42);
   haxes__37->GetYaxis()->SetLabelOffset(0.007);
   haxes__37->GetYaxis()->SetLabelSize(0.05);
   haxes__37->GetYaxis()->SetTitleSize(0.06);
   haxes__37->GetYaxis()->SetTitleOffset(1.25);
   haxes__37->GetYaxis()->SetTitleFont(42);
   haxes__37->GetZaxis()->SetLabelFont(42);
   haxes__37->GetZaxis()->SetLabelOffset(0.007);
   haxes__37->GetZaxis()->SetLabelSize(0.05);
   haxes__37->GetZaxis()->SetTitleSize(0.06);
   haxes__37->GetZaxis()->SetTitleFont(42);
   haxes__37->Draw("");
   TBox *box = new TBox(400,0.9468748,420,1.053125);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(400,0.9468748,420,1.053125);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3079[16] = {
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
   Double_t bin_0_syst_fy3079[16] = {
   0.218933,
   -0.2269219,
   0.2666771,
   0.3007051,
   0.3263139,
   0.3731811,
   -0.3339235,
   0.3942181,
   0.4531517,
   0.4776882,
   0.5641537,
   0.5994745,
   0.6490349,
   0.7336864,
   -0.5242164,
   0.7960212};
   Double_t bin_0_syst_felx3079[16] = {
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
   Double_t bin_0_syst_fely3079[16] = {
   0.01471854,
   0.006807658,
   0.0179837,
   0.01962043,
   0.02219877,
   0.02560975,
   0.0113534,
   0.02793975,
   0.03287492,
   0.03720098,
   0.04641678,
   0.05350703,
   0.0663091,
   0.09093968,
   0.02306552,
   0.0907996};
   Double_t bin_0_syst_fehx3079[16] = {
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
   Double_t bin_0_syst_fehy3079[16] = {
   0.01385116,
   0.004311517,
   0.01688436,
   0.01848863,
   0.02088646,
   0.02560975,
   0.008682011,
   0.02725556,
   0.03258318,
   0.03652002,
   0.04683658,
   0.0544656,
   0.06855796,
   0.09896706,
   0.0351225,
   0.1367286};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(16,bin_0_syst_fx3079,bin_0_syst_fy3079,bin_0_syst_felx3079,bin_0_syst_fehx3079,bin_0_syst_fely3079,bin_0_syst_fehy3079);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3079 = new TH1F("Graph_bin_0_syst3079","Graph",100,0,427.9003);
   Graph_bin_0_syst3079->SetMinimum(-0.6952851);
   Graph_bin_0_syst3079->SetMaximum(1.080753);
   Graph_bin_0_syst3079->SetDirectory(0);
   Graph_bin_0_syst3079->SetStats(0);
   Graph_bin_0_syst3079->SetLineStyle(0);
   Graph_bin_0_syst3079->SetMarkerStyle(20);
   Graph_bin_0_syst3079->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3079->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3079->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3079->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3079->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3079->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3079->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3079->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3079->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3079->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3079->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3079->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3079->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3079->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3079->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3079->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3079->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3079);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3080[16] = {
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
   Double_t bin_0_fy3080[16] = {
   0.218933,
   -0.2269219,
   0.2666771,
   0.3007051,
   0.3263139,
   0.3731811,
   -0.3339235,
   0.3942181,
   0.4531517,
   0.4776882,
   0.5641537,
   0.5994745,
   0.6490349,
   0.7336864,
   -0.5242164,
   0.7960212};
   Double_t bin_0_felx3080[16] = {
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
   Double_t bin_0_fely3080[16] = {
   0.005480322,
   0.004295478,
   0.006866621,
   0.008050311,
   0.009113084,
   0.01103663,
   0.004851973,
   0.01266949,
   0.01425902,
   0.01673658,
   0.02277298,
   0.02619523,
   0.02552463,
   0.04092289,
   0.008241533,
   0.06630805};
   Double_t bin_0_fehx3080[16] = {
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
   Double_t bin_0_fehy3080[16] = {
   0.005480322,
   0.004295478,
   0.006866621,
   0.008050311,
   0.009113084,
   0.01103663,
   0.004851973,
   0.01266949,
   0.01425902,
   0.01673658,
   0.02277298,
   0.02619523,
   0.02552463,
   0.04092289,
   0.008241533,
   0.06630805};
   grae = new TGraphAsymmErrors(16,bin_0_fx3080,bin_0_fy3080,bin_0_felx3080,bin_0_fehx3080,bin_0_fely3080,bin_0_fehy3080);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03080 = new TH1F("Graph_bin_03080","Graph",100,0,423.9403);
   Graph_bin_03080->SetMinimum(-0.6719366);
   Graph_bin_03080->SetMaximum(1.001808);
   Graph_bin_03080->SetDirectory(0);
   Graph_bin_03080->SetStats(0);
   Graph_bin_03080->SetLineStyle(0);
   Graph_bin_03080->SetMarkerStyle(20);
   Graph_bin_03080->GetXaxis()->SetLabelFont(42);
   Graph_bin_03080->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03080->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03080->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03080->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03080->GetXaxis()->SetTitleFont(42);
   Graph_bin_03080->GetYaxis()->SetLabelFont(42);
   Graph_bin_03080->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03080->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03080->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03080->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03080->GetYaxis()->SetTitleFont(42);
   Graph_bin_03080->GetZaxis()->SetLabelFont(42);
   Graph_bin_03080->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03080->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03080->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03080->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03080);
   
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
   
   TH1F *haxes_copy__38 = new TH1F("haxes_copy__38","haxes",1,0,420);
   haxes_copy__38->SetMinimum(0);
   haxes_copy__38->SetMaximum(1.5);
   haxes_copy__38->SetDirectory(0);
   haxes_copy__38->SetLineStyle(0);
   haxes_copy__38->SetMarkerStyle(20);
   haxes_copy__38->GetXaxis()->SetTitle("N_{part}");
   haxes_copy__38->GetXaxis()->CenterTitle(true);
   haxes_copy__38->GetXaxis()->SetLabelFont(42);
   haxes_copy__38->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__38->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__38->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__38->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__38->GetXaxis()->SetTitleFont(42);
   haxes_copy__38->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__38->GetYaxis()->SetLabelFont(42);
   haxes_copy__38->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__38->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__38->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__38->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__38->GetYaxis()->SetTitleFont(42);
   haxes_copy__38->GetZaxis()->SetLabelFont(42);
   haxes_copy__38->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__38->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__38->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__38->GetZaxis()->SetTitleFont(42);
   haxes_copy__38->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
