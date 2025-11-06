void result_JPsi_RAA_pt_nonprompt_14005_2_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:54 2018) by ROOT version6.06/00
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-10,-0.225,52.5,1.65);
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
   
   TH1F *haxes__83 = new TH1F("haxes__83","haxes",1,0,50);
   haxes__83->SetMinimum(0);
   haxes__83->SetMaximum(1.5);
   haxes__83->SetLineStyle(0);
   haxes__83->SetMarkerStyle(20);
   haxes__83->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__83->GetXaxis()->CenterTitle(true);
   haxes__83->GetXaxis()->SetLabelFont(42);
   haxes__83->GetXaxis()->SetLabelOffset(0.007);
   haxes__83->GetXaxis()->SetLabelSize(0.05);
   haxes__83->GetXaxis()->SetTitleSize(0.06);
   haxes__83->GetXaxis()->SetTitleOffset(0.9);
   haxes__83->GetXaxis()->SetTitleFont(42);
   haxes__83->GetYaxis()->SetTitle("R_{AA}");
   haxes__83->GetYaxis()->SetLabelFont(42);
   haxes__83->GetYaxis()->SetLabelOffset(0.007);
   haxes__83->GetYaxis()->SetLabelSize(0.05);
   haxes__83->GetYaxis()->SetTitleSize(0.06);
   haxes__83->GetYaxis()->SetTitleOffset(1.25);
   haxes__83->GetYaxis()->SetTitleFont(42);
   haxes__83->GetZaxis()->SetLabelFont(42);
   haxes__83->GetZaxis()->SetLabelOffset(0.007);
   haxes__83->GetZaxis()->SetLabelSize(0.05);
   haxes__83->GetZaxis()->SetTitleSize(0.06);
   haxes__83->GetZaxis()->SetTitleFont(42);
   haxes__83->Draw("");
   TBox *box = new TBox(47.5,0.9543382,50,1.041388);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(47.5,0.9543382,50,1.041388);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3175[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_syst_fy3175[12] = {
   0.5205617,
   0.4780894,
   0.4616475,
   0.4371907,
   0.4263492,
   0.4585017,
   0.4234693,
   0.4502758,
   0.4407601,
   0.4666347,
   0.6018213,
   0.473903};
   Double_t bin_0_syst_felx3175[12] = {
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
   0.5,
   0.5};
   Double_t bin_0_syst_fely3175[12] = {
   0.05704963,
   0.03986511,
   0.03372379,
   0.02713588,
   0.02404351,
   0.02512273,
   0.02306843,
   0.02731351,
   0.02398772,
   0.02865439,
   0.05034567,
   0.04050557};
   Double_t bin_0_syst_fehx3175[12] = {
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
   0.5,
   0.5};
   Double_t bin_0_syst_fehy3175[12] = {
   0.05704963,
   0.03986511,
   0.03372379,
   0.02713588,
   0.02404351,
   0.02512273,
   0.02306843,
   0.02731351,
   0.02398772,
   0.02865439,
   0.05034567,
   0.04050557};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(12,bin_0_syst_fx3175,bin_0_syst_fy3175,bin_0_syst_felx3175,bin_0_syst_fehx3175,bin_0_syst_fely3175,bin_0_syst_fehy3175);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3175 = new TH1F("Graph_bin_0_syst3175","Graph",100,2.85,46.65);
   Graph_bin_0_syst3175->SetMinimum(0.3752242);
   Graph_bin_0_syst3175->SetMaximum(0.6773436);
   Graph_bin_0_syst3175->SetDirectory(0);
   Graph_bin_0_syst3175->SetStats(0);
   Graph_bin_0_syst3175->SetLineStyle(0);
   Graph_bin_0_syst3175->SetMarkerStyle(20);
   Graph_bin_0_syst3175->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3175->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3175->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3175->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3175->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3175->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3175->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3175->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3175->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3175->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3175->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3175->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3175->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3175->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3175->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3175->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3175->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3175);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3176[12] = {
   7,
   8,
   9,
   10.25,
   12,
   14,
   16.25,
   18.75,
   22.5,
   27.5,
   32.5,
   42.5};
   Double_t bin_0_fy3176[12] = {
   0.5205617,
   0.4780894,
   0.4616475,
   0.4371907,
   0.4263492,
   0.4585017,
   0.4234693,
   0.4502758,
   0.4407601,
   0.4666347,
   0.6018213,
   0.473903};
   Double_t bin_0_felx3176[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_fely3176[12] = {
   0.01781202,
   0.01397327,
   0.01285977,
   0.0107817,
   0.01086009,
   0.01418445,
   0.01554696,
   0.02139503,
   0.02221882,
   0.03530909,
   0.06675377,
   0.06164125};
   Double_t bin_0_fehx3176[12] = {
   0.5,
   0.5,
   0.5,
   0.75,
   1,
   1,
   1.25,
   1.25,
   2.5,
   2.5,
   2.5,
   7.5};
   Double_t bin_0_fehy3176[12] = {
   0.01781202,
   0.01397327,
   0.01285977,
   0.0107817,
   0.01086009,
   0.01418445,
   0.01554696,
   0.02139503,
   0.02221882,
   0.03530909,
   0.06675377,
   0.06164125};
   grae = new TGraphAsymmErrors(12,bin_0_fx3176,bin_0_fy3176,bin_0_felx3176,bin_0_fehx3176,bin_0_fely3176,bin_0_fehy3176);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03176 = new TH1F("Graph_bin_03176","Graph",100,2.15,54.35);
   Graph_bin_03176->SetMinimum(0.381857);
   Graph_bin_03176->SetMaximum(0.6946404);
   Graph_bin_03176->SetDirectory(0);
   Graph_bin_03176->SetStats(0);
   Graph_bin_03176->SetLineStyle(0);
   Graph_bin_03176->SetMarkerStyle(20);
   Graph_bin_03176->GetXaxis()->SetLabelFont(42);
   Graph_bin_03176->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03176->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03176->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03176->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03176->GetXaxis()->SetTitleFont(42);
   Graph_bin_03176->GetYaxis()->SetLabelFont(42);
   Graph_bin_03176->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03176->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03176->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03176->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03176->GetYaxis()->SetTitleFont(42);
   Graph_bin_03176->GetZaxis()->SetLabelFont(42);
   Graph_bin_03176->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03176->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03176->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03176->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03176);
   
   grae->Draw("p");
   
   TLegend *leg = new TLegend(0.19,0.49,0.51,0.64,NULL,"brNDC");
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
   TLine *line = new TLine(0,1,50,1);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.78,"|y| < 2.4");
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
   
   TH1F *haxes_copy__84 = new TH1F("haxes_copy__84","haxes",1,0,50);
   haxes_copy__84->SetMinimum(0);
   haxes_copy__84->SetMaximum(1.5);
   haxes_copy__84->SetDirectory(0);
   haxes_copy__84->SetLineStyle(0);
   haxes_copy__84->SetMarkerStyle(20);
   haxes_copy__84->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__84->GetXaxis()->CenterTitle(true);
   haxes_copy__84->GetXaxis()->SetLabelFont(42);
   haxes_copy__84->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__84->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__84->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__84->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__84->GetXaxis()->SetTitleFont(42);
   haxes_copy__84->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__84->GetYaxis()->SetLabelFont(42);
   haxes_copy__84->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__84->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__84->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__84->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__84->GetYaxis()->SetTitleFont(42);
   haxes_copy__84->GetZaxis()->SetLabelFont(42);
   haxes_copy__84->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__84->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__84->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__84->GetZaxis()->SetTitleFont(42);
   haxes_copy__84->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
