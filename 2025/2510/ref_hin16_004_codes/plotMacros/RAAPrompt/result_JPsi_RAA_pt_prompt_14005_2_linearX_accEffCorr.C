void result_JPsi_RAA_pt_prompt_14005_2_linearX_accEffCorr()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Nov  5 20:50:51 2018) by ROOT version6.06/00
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
   
   TH1F *haxes__35 = new TH1F("haxes__35","haxes",1,0,50);
   haxes__35->SetMinimum(0);
   haxes__35->SetMaximum(1.5);
   haxes__35->SetLineStyle(0);
   haxes__35->SetMarkerStyle(20);
   haxes__35->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes__35->GetXaxis()->CenterTitle(true);
   haxes__35->GetXaxis()->SetLabelFont(42);
   haxes__35->GetXaxis()->SetLabelOffset(0.007);
   haxes__35->GetXaxis()->SetLabelSize(0.05);
   haxes__35->GetXaxis()->SetTitleSize(0.06);
   haxes__35->GetXaxis()->SetTitleOffset(0.9);
   haxes__35->GetXaxis()->SetTitleFont(42);
   haxes__35->GetYaxis()->SetTitle("R_{AA}");
   haxes__35->GetYaxis()->SetLabelFont(42);
   haxes__35->GetYaxis()->SetLabelOffset(0.007);
   haxes__35->GetYaxis()->SetLabelSize(0.05);
   haxes__35->GetYaxis()->SetTitleSize(0.06);
   haxes__35->GetYaxis()->SetTitleOffset(1.25);
   haxes__35->GetYaxis()->SetTitleFont(42);
   haxes__35->GetZaxis()->SetLabelFont(42);
   haxes__35->GetZaxis()->SetLabelOffset(0.007);
   haxes__35->GetZaxis()->SetLabelSize(0.05);
   haxes__35->GetZaxis()->SetTitleSize(0.06);
   haxes__35->GetZaxis()->SetTitleFont(42);
   haxes__35->Draw("");
   TBox *box = new TBox(47.5,0.9543382,50,1.041388);
   box->SetFillColor(-1);
   box->Draw();
   box = new TBox(47.5,0.9543382,50,1.041388);
   box->SetFillColor(-1);
   box->SetFillStyle(0);
   box->Draw();
   
   Double_t bin_0_syst_fx3077[12] = {
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
   Double_t bin_0_syst_fy3077[12] = {
   0.3609823,
   0.3480969,
   0.3337932,
   0.3303158,
   0.3394523,
   0.3601479,
   0.3595278,
   0.3639327,
   0.4245342,
   0.5205356,
   0.4952924,
   0.5209613};
   Double_t bin_0_syst_felx3077[12] = {
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
   Double_t bin_0_syst_fely3077[12] = {
   0.03460031,
   0.02695614,
   0.02188608,
   0.01922661,
   0.01852416,
   0.01922774,
   0.0190938,
   0.0225332,
   0.02257819,
   0.03536442,
   0.04784256,
   0.04748044};
   Double_t bin_0_syst_fehx3077[12] = {
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
   Double_t bin_0_syst_fehy3077[12] = {
   0.03460031,
   0.02695614,
   0.02188608,
   0.01922661,
   0.01852416,
   0.01922774,
   0.0190938,
   0.0225332,
   0.02257819,
   0.03536442,
   0.04784256,
   0.04748044};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(12,bin_0_syst_fx3077,bin_0_syst_fy3077,bin_0_syst_felx3077,bin_0_syst_fehx3077,bin_0_syst_fely3077,bin_0_syst_fehy3077);
   grae->SetName("bin_0_syst");
   grae->SetTitle("Graph");
   grae->SetFillColor(-1);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_bin_0_syst3077 = new TH1F("Graph_bin_0_syst3077","Graph",100,2.85,46.65);
   Graph_bin_0_syst3077->SetMinimum(0.285354);
   Graph_bin_0_syst3077->SetMaximum(0.5941769);
   Graph_bin_0_syst3077->SetDirectory(0);
   Graph_bin_0_syst3077->SetStats(0);
   Graph_bin_0_syst3077->SetLineStyle(0);
   Graph_bin_0_syst3077->SetMarkerStyle(20);
   Graph_bin_0_syst3077->GetXaxis()->SetLabelFont(42);
   Graph_bin_0_syst3077->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3077->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3077->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3077->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_0_syst3077->GetXaxis()->SetTitleFont(42);
   Graph_bin_0_syst3077->GetYaxis()->SetLabelFont(42);
   Graph_bin_0_syst3077->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3077->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3077->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3077->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_0_syst3077->GetYaxis()->SetTitleFont(42);
   Graph_bin_0_syst3077->GetZaxis()->SetLabelFont(42);
   Graph_bin_0_syst3077->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_0_syst3077->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_0_syst3077->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_0_syst3077->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_0_syst3077);
   
   grae->Draw("5");
   
   Double_t bin_0_fx3078[12] = {
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
   Double_t bin_0_fy3078[12] = {
   0.3609823,
   0.3480969,
   0.3337932,
   0.3303158,
   0.3394523,
   0.3601479,
   0.3595278,
   0.3639327,
   0.4245342,
   0.5205356,
   0.4952924,
   0.5209613};
   Double_t bin_0_felx3078[12] = {
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
   Double_t bin_0_fely3078[12] = {
   0.01245823,
   0.01017544,
   0.009213545,
   0.00806029,
   0.008512602,
   0.01097526,
   0.0129574,
   0.0170388,
   0.02106728,
   0.03893998,
   0.05461896,
   0.06762865};
   Double_t bin_0_fehx3078[12] = {
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
   Double_t bin_0_fehy3078[12] = {
   0.01245823,
   0.01017544,
   0.009213545,
   0.00806029,
   0.008512602,
   0.01097526,
   0.0129574,
   0.0170388,
   0.02106728,
   0.03893998,
   0.05461896,
   0.06762865};
   grae = new TGraphAsymmErrors(12,bin_0_fx3078,bin_0_fy3078,bin_0_felx3078,bin_0_fehx3078,bin_0_fely3078,bin_0_fehy3078);
   grae->SetName("bin_0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(28);
   grae->SetMarkerSize(1.5);
   
   TH1F *Graph_bin_03078 = new TH1F("Graph_bin_03078","Graph",100,2.15,54.35);
   Graph_bin_03078->SetMinimum(0.2956221);
   Graph_bin_03078->SetMaximum(0.6152233);
   Graph_bin_03078->SetDirectory(0);
   Graph_bin_03078->SetStats(0);
   Graph_bin_03078->SetLineStyle(0);
   Graph_bin_03078->SetMarkerStyle(20);
   Graph_bin_03078->GetXaxis()->SetLabelFont(42);
   Graph_bin_03078->GetXaxis()->SetLabelOffset(0.007);
   Graph_bin_03078->GetXaxis()->SetLabelSize(0.05);
   Graph_bin_03078->GetXaxis()->SetTitleSize(0.06);
   Graph_bin_03078->GetXaxis()->SetTitleOffset(0.9);
   Graph_bin_03078->GetXaxis()->SetTitleFont(42);
   Graph_bin_03078->GetYaxis()->SetLabelFont(42);
   Graph_bin_03078->GetYaxis()->SetLabelOffset(0.007);
   Graph_bin_03078->GetYaxis()->SetLabelSize(0.05);
   Graph_bin_03078->GetYaxis()->SetTitleSize(0.06);
   Graph_bin_03078->GetYaxis()->SetTitleOffset(1.25);
   Graph_bin_03078->GetYaxis()->SetTitleFont(42);
   Graph_bin_03078->GetZaxis()->SetLabelFont(42);
   Graph_bin_03078->GetZaxis()->SetLabelOffset(0.007);
   Graph_bin_03078->GetZaxis()->SetLabelSize(0.05);
   Graph_bin_03078->GetZaxis()->SetTitleSize(0.06);
   Graph_bin_03078->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_bin_03078);
   
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
   
   TH1F *haxes_copy__36 = new TH1F("haxes_copy__36","haxes",1,0,50);
   haxes_copy__36->SetMinimum(0);
   haxes_copy__36->SetMaximum(1.5);
   haxes_copy__36->SetDirectory(0);
   haxes_copy__36->SetLineStyle(0);
   haxes_copy__36->SetMarkerStyle(20);
   haxes_copy__36->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   haxes_copy__36->GetXaxis()->CenterTitle(true);
   haxes_copy__36->GetXaxis()->SetLabelFont(42);
   haxes_copy__36->GetXaxis()->SetLabelOffset(0.007);
   haxes_copy__36->GetXaxis()->SetLabelSize(0.05);
   haxes_copy__36->GetXaxis()->SetTitleSize(0.06);
   haxes_copy__36->GetXaxis()->SetTitleOffset(0.9);
   haxes_copy__36->GetXaxis()->SetTitleFont(42);
   haxes_copy__36->GetYaxis()->SetTitle("R_{AA}");
   haxes_copy__36->GetYaxis()->SetLabelFont(42);
   haxes_copy__36->GetYaxis()->SetLabelOffset(0.007);
   haxes_copy__36->GetYaxis()->SetLabelSize(0.05);
   haxes_copy__36->GetYaxis()->SetTitleSize(0.06);
   haxes_copy__36->GetYaxis()->SetTitleOffset(1.25);
   haxes_copy__36->GetYaxis()->SetTitleFont(42);
   haxes_copy__36->GetZaxis()->SetLabelFont(42);
   haxes_copy__36->GetZaxis()->SetLabelOffset(0.007);
   haxes_copy__36->GetZaxis()->SetLabelSize(0.05);
   haxes_copy__36->GetZaxis()->SetTitleSize(0.06);
   haxes_copy__36->GetZaxis()->SetTitleFont(42);
   haxes_copy__36->Draw("sameaxis");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
