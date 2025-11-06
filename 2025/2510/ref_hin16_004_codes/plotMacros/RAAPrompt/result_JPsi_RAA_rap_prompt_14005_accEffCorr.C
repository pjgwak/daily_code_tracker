void result_JPsi_RAA_rap_prompt_14005_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:52 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__39 = new TH1F("haxes__39","haxes",1,0,2.4);
   haxes__39->SetMinimum(0);
   haxes__39->SetMaximum(1.5);
   haxes__39->SetLineStyle(0);
   haxes__39->SetMarkerStyle(20);
   haxes__39->GetXaxis()->SetTitle("|y|");
   haxes__39->GetXaxis()->CenterTitle(true);
   haxes__39->GetXaxis()->SetNdivisions(-306);
   haxes__39->GetXaxis()->SetLabelFont(42);
   haxes__39->GetXaxis()->SetLabelOffset(0.007);
   haxes__39->GetXaxis()->SetLabelSize(0.05);
   haxes__39->GetXaxis()->SetTitleSize(0.06);
   haxes__39->GetXaxis()->SetTitleOffset(0.9);
   haxes__39->GetXaxis()->SetTitleFont(42);
   haxes__39->GetYaxis()->SetTitle("R_{AA}");
   haxes__39->GetYaxis()->SetLabelFont(42);
   haxes__39->GetYaxis()->SetLabelOffset(0.007);
   haxes__39->GetYaxis()->SetLabelSize(0.05);
   haxes__39->GetYaxis()->SetTitleSize(0.06);
   haxes__39->GetYaxis()->SetTitleOffset(1.25);
   haxes__39->GetYaxis()->SetTitleFont(42);
   haxes__39->GetZaxis()->SetLabelFont(42);
   haxes__39->GetZaxis()->SetLabelOffset(0.007);
   haxes__39->GetZaxis()->SetLabelSize(0.05);
   haxes__39->GetZaxis()->SetTitleSize(0.06);
   haxes__39->GetZaxis()->SetTitleFont(42);
   haxes__39->Draw("");
   TBox *box = new TBox(2.28,0.9543382,2.4,1.041388);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(2.28,0.9543382,2.4,1.041388);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3081[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_syst_fy3081[6] = {
   0.3562073,
   0.3507093,
   0.3446304,
   0.3345537,
   0.3240817,
   0.3181873};
   Double_t bin_0_syst_felx3081[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_fely3081[6] = {
   0.02364558,
   0.02406192,
   0.02412092,
   0.02675539,
   0.02635942,
   0.02404753};
   Double_t bin_0_syst_fehx3081[6] = {
   0.05,
   0.05,
   0.05,
   0.05,
   0.05,
   0.05};
   Double_t bin_0_syst_fehy3081[6] = {
   0.02364558,
   0.02406192,
   0.02412092,
   0.02675539,
   0.02635942,
   0.02404753};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,bin_0_syst_fx3081,bin_0_syst_fy3081,bin_0_syst_felx3081,bin_0_syst_fehx3081,bin_0_syst_fely3081,bin_0_syst_fehy3081);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3081 = new TH1F("Graph_bin_0_syst3081","Graph",100,0,2.46);
   Graph_bin_0_syst3081->SetMinimum(0.2855684);
   Graph_bin_0_syst3081->SetMaximum(0.3884242);
   Graph_bin_0_syst3081->SetDirectory(0);
   Graph_bin_0_syst3081->SetStats(0);
   Graph_bin_0_syst3081->SetLineStyle(0);
   Graph_bin_0_syst3081->SetMarkerStyle(20);
   Graph_bin_0_syst3081->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3081->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3081->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3081->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3081->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3081->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3081->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3081->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3081->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3081->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3081->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3081->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3081->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3081->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3081->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3081->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3081->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3081);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3082[6] = {
   0.2,
   0.6,
   1,
   1.4,
   1.8,
   2.2};
   Double_t bin_0_fy3082[6] = {
   0.3562073,
   0.3507093,
   0.3446304,
   0.3345537,
   0.3240817,
   0.3181873};
   Double_t bin_0_felx3082[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_fely3082[6] = {
   0.008322276,
   0.007955105,
   0.007930847,
   0.007839923,
   0.008178086,
   0.01294015};
   Double_t bin_0_fehx3082[6] = {
   0.2,
   0.2,
   0.2,
   0.2,
   0.2,
   0.2};
   Double_t bin_0_fehy3082[6] = {
   0.008322276,
   0.007955105,
   0.007930847,
   0.007839923,
   0.008178086,
   0.01294015};
   grae = new TGraphAsymmErrors(6,bin_0_fx3082,bin_0_fy3082,bin_0_felx3082,bin_0_fehx3082,bin_0_fely3082,bin_0_fehy3082);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03082 = new TH1F("Graph_bin_03082","Graph",100,0,2.64);
   Graph_bin_03082->SetMinimum(0.2993189);
   Graph_bin_03082->SetMaximum(0.3704578);
   Graph_bin_03082->SetDirectory(0);
   Graph_bin_03082->SetStats(0);
   Graph_bin_03082->SetLineStyle(0);
   Graph_bin_03082->SetMarkerStyle(20);
   Graph_bin_03082->GetXaxis()->SetLabelFont(42);
   Graph_bin_03082->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03082->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03082->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03082->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03082->GetXaxis()->SetTitleFont(42);
   Graph_bin_03082->GetYaxis()->SetLabelFont(42);
   Graph_bin_03082->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03082->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03082->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03082->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03082->GetYaxis()->SetTitleFont(42);
   Graph_bin_03082->GetZaxis()->SetLabelFont(42);
   Graph_bin_03082->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03082->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03082->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03082->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03082);
   
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
   
   TH1F *haxes_copy__40 = new TH1F("haxes_copy__40","haxes",1,0,2.4);
   haxes_copy__40->SetMinimum(0);
   haxes_copy__40->SetMaximum(1.5);
   haxes_copy__40->SetDirectory(0);
   haxes_copy__40->SetLineStyle(0);
   haxes_copy__40->SetMarkerStyle(20);
   haxes_copy__40->GetXaxis()->SetTitle("|y|");
   haxes_copy__40->GetXaxis()->CenterTitle(true);
   haxes_copy__40->GetXaxis()->SetNdivisions(-306);
   haxes_copy__40->GetXaxis()->SetLabelFont(42);
   haxes_copy__40->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__40->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__40->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__40->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__40->GetXaxis()->SetTitleFont(42);
   haxes_copy__40->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__40->GetYaxis()->SetLabelFont(42);
   haxes_copy__40->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__40->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__40->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__40->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__40->GetYaxis()->SetTitleFont(42);
   haxes_copy__40->GetZaxis()->SetLabelFont(42);
   haxes_copy__40->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__40->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__40->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__40->GetZaxis()->SetTitleFont(42);
   haxes_copy__40->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
