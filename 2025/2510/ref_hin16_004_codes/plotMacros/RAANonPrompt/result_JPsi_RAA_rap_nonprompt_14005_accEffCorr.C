void result_JPsi_RAA_rap_nonprompt_14005_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.48,-0.225,2.52,1.65);
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
   
   TH1F *haxes__87 = new TH1F("haxes__87","haxes",1,0,2.4);
   haxes__87->SetMinimum(0);
   haxes__87->SetMaximum(1.5);
   haxes__87->SetLineStyle(0);
   haxes__87->SetMarkerStyle(20);
   haxes__87->GetXaxis()->SetTitle("|y|");
   haxes__87->GetXaxis()->CenterTitle(true);
   haxes__87->GetXaxis()->SetNdivisions(-306);
   haxes__87->GetXaxis()->SetLabelFont(42);
   haxes__87->GetXaxis()->SetLabelOffset(0.007);
   haxes__87->GetXaxis()->SetLabelSize(0.05);
   haxes__87->GetXaxis()->SetTitleSize(0.06);
   haxes__87->GetXaxis()->SetTitleOffset(0.9);
   haxes__87->GetXaxis()->SetTitleFont(42);
   haxes__87->GetYaxis()->SetTitle("R_{AA}");
   haxes__87->GetYaxis()->SetLabelFont(42);
   haxes__87->GetYaxis()->SetLabelOffset(0.007);
   haxes__87->GetYaxis()->SetLabelSize(0.05);
   haxes__87->GetYaxis()->SetTitleSize(0.06);
   haxes__87->GetYaxis()->SetTitleOffset(1.25);
   haxes__87->GetYaxis()->SetTitleFont(42);
   haxes__87->GetZaxis()->SetLabelFont(42);
   haxes__87->GetZaxis()->SetLabelOffset(0.007);
   haxes__87->GetZaxis()->SetLabelSize(0.05);
   haxes__87->GetZaxis()->SetTitleSize(0.06);
   haxes__87->GetZaxis()->SetTitleFont(42);
   haxes__87->Draw("");
   TBox *box = new TBox(2.28,0.9543382,2.4,1.041388);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(2.28,0.9543382,2.4,1.041388);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3179[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_fy3179[6] = {
   0.4475415,
   0.4767508,
   0.437081,
   0.4485272,
   0.4691107,
   0.4292723};
   Double_t bin_0_syst_felx3179[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_fely3179[6] = {
   0.02733338,
   0.03275594,
   0.03166435,
   0.04008393,
   0.04291875,
   0.03839543};
   Double_t bin_0_syst_fehx3179[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_fehy3179[6] = {
   0.02733338,
   0.03275594,
   0.03166435,
   0.04008393,
   0.04291875,
   0.03839543};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3179,bin_0_syst_fy3179,bin_0_syst_felx3179,bin_0_syst_fehx3179,bin_0_syst_fely3179,bin_0_syst_fehy3179);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3179 = new TH1F("Graph_bin_0_syst3179","Graph",100,0,2.46);
   Graph_bin_0_syst3179->SetMinimum(0.3787616);
   Graph_bin_0_syst3179->SetMaximum(0.5241447);
   Graph_bin_0_syst3179->SetDirectory(0);
   Graph_bin_0_syst3179->SetStats(0);
   Graph_bin_0_syst3179->SetLineStyle(0);
   Graph_bin_0_syst3179->SetMarkerStyle(20);
   Graph_bin_0_syst3179->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3179->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3179->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3179->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3179->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3179->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3179->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3179->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3179->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3179->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3179->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3179->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3179->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3179->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3179->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3179->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3179->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3179);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3180[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_fy3180[6] = {
   0.4475415,
   0.4767508,
   0.437081,
   0.4485272,
   0.4691107,
   0.4292723};
   Double_t bin_0_felx3180[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_fely3180[6] = {
   0.01049904,
   0.01087813,
   0.0101508,
   0.01047291,
   0.01173123,
   0.01753272};
   Double_t bin_0_fehx3180[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_fehy3180[6] = {
   0.01049904,
   0.01087813,
   0.0101508,
   0.01047291,
   0.01173123,
   0.01753272};
   grae = new TGraphAsymmErrors(6,bin_0_fx3180,bin_0_fy3180,bin_0_felx3180,bin_0_fehx3180,bin_0_fely3180,bin_0_fehy3180);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03180 = new TH1F("Graph_bin_03180","Graph",100,0,2.64);
   Graph_bin_03180->SetMinimum(0.4041507);
   Graph_bin_03180->SetMaximum(0.4952179);
   Graph_bin_03180->SetDirectory(0);
   Graph_bin_03180->SetStats(0);
   Graph_bin_03180->SetLineStyle(0);
   Graph_bin_03180->SetMarkerStyle(20);
   Graph_bin_03180->GetXaxis()->SetLabelFont(42);
   Graph_bin_03180->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03180->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03180->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03180->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03180->GetXaxis()->SetTitleFont(42);
   Graph_bin_03180->GetYaxis()->SetLabelFont(42);
   Graph_bin_03180->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03180->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03180->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03180->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03180->GetYaxis()->SetTitleFont(42);
   Graph_bin_03180->GetZaxis()->SetLabelFont(42);
   Graph_bin_03180->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03180->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03180->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03180->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03180);
   
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
   TLine *line = new TLine(0,1,2.4,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.044);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.73,"Cent. 0-100%");
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
   
   TH1F *haxes_copy__88 = new TH1F("haxes_copy__88","haxes",1,0,2.4);
   haxes_copy__88->SetMinimum(0);
   haxes_copy__88->SetMaximum(1.5);
   haxes_copy__88->SetDirectory(0);
   haxes_copy__88->SetLineStyle(0);
   haxes_copy__88->SetMarkerStyle(20);
   haxes_copy__88->GetXaxis()->SetTitle("|y|");
   haxes_copy__88->GetXaxis()->CenterTitle(true);
   haxes_copy__88->GetXaxis()->SetNdivisions(-306);
   haxes_copy__88->GetXaxis()->SetLabelFont(42);
   haxes_copy__88->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__88->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__88->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__88->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__88->GetXaxis()->SetTitleFont(42);
   haxes_copy__88->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__88->GetYaxis()->SetLabelFont(42);
   haxes_copy__88->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__88->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__88->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__88->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__88->GetYaxis()->SetTitleFont(42);
   haxes_copy__88->GetZaxis()->SetLabelFont(42);
   haxes_copy__88->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__88->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__88->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__88->GetZaxis()->SetTitleFont(42);
   haxes_copy__88->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
