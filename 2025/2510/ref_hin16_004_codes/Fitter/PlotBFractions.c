#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "Riostream.h"
#include "Riostream.h"
#include <math.h>
#include "TStyle.h"
#include "TCanvas.h"
#include <TROOT.h>

void PlotBFractions() //Last updated 1 June 2016
{
	gStyle->SetOptStat(000000000);
	gStyle->SetOptFit(0);
	gStyle->SetEndErrorSize(5);
	gStyle->SetLineWidth(1);
	gStyle->SetErrorX(0);

	gStyle->SetLabelFont(62, "xyz");
	gStyle->SetTitleFont(62, "xyzt");
	gStyle->SetLabelSize(0.04, "xyz");
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetPalette(1, 0);
	gStyle->SetTitleSize(0.07, "t");
	gStyle->SetTitleSize(0.04, "xyz");
	gStyle->SetMarkerSize(1.3);


	const int nbins_ptmid = 5;
	const double bins_ptmid[nbins_ptmid + 1] = { 6.5, 9, 12, 15, 20, 30 };
	const int nbins_ptfwd = 3;
	const double bins_ptfwd[nbins_ptfwd + 1] = { 3, 6.5, 12, 30 };

	const int nbins_centmid = 6;
	const double bins_centmid[nbins_centmid + 1] = { 0, 10, 20, 30, 40, 50, 100 };
	const int nbins_centfwd = 3;
	const double bins_centfwd[nbins_centfwd + 1] = { 0, 20, 40, 100 };

	double zero_xerrminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0 }; //to hide errors along x axis (setErrorX(0) didn't work for TGraphAsymmErrors)
	double zero_xerrplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	//TGaxis::SetMaxDigits(3);

	/////pp plots

	// 1S (pp)

	TH1D* hJpsiN_pp_rap0016 = new TH1D("hJpsiN_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* hJpsiN_pp_rap1624 = new TH1D("hJpsiN_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* hFitAndre_pp_rap0016 = new TH1D("hFitAndre_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* hFitAndre_pp_rap1624 = new TH1D("hFitAndre_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* hFitMihee_pp_rap0016 = new TH1D("hFitMihee_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* hFitMihee_pp_rap1624 = new TH1D("hFitMihee_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);


	// filling pp histograms (by hand)

	hJpsiN_pp_rap0016->SetBinContent(1, 0.243);//6.5-9
	hJpsiN_pp_rap0016->SetBinContent(2, 0.318);//9-12
	hJpsiN_pp_rap0016->SetBinContent(3, 0.396);//12-15
	hJpsiN_pp_rap0016->SetBinContent(4, 0.481);//15-20
	hJpsiN_pp_rap0016->SetBinContent(5, 0.556);//20-30

	hJpsiN_pp_rap0016->SetBinError(1, 0.005);//6.5-9
	hJpsiN_pp_rap0016->SetBinError(2, 0.004);//9-12
	hJpsiN_pp_rap0016->SetBinError(3, 0.007);//12-15
	hJpsiN_pp_rap0016->SetBinError(4, 0.008);//15-20
	hJpsiN_pp_rap0016->SetBinError(5, 0.012);//20-30

	hJpsiN_pp_rap1624->SetBinContent(1, 0.177);//3-6.5
	hJpsiN_pp_rap1624->SetBinContent(2, 0.259);//6.5-12
	hJpsiN_pp_rap1624->SetBinContent(3, 0.431);//12-30

	hJpsiN_pp_rap1624->SetBinError(1, 0.005);//3-6.5
	hJpsiN_pp_rap1624->SetBinError(2, 0.005);//6.5-12
	hJpsiN_pp_rap1624->SetBinError(3, 0.009);//12-30

	hFitAndre_pp_rap0016->SetBinContent(1, 0.234);//6.5-9
	hFitAndre_pp_rap0016->SetBinContent(2, 0.305);//9-12
	hFitAndre_pp_rap0016->SetBinContent(3, 0.385);//12-15
	hFitAndre_pp_rap0016->SetBinContent(4, 0.462);//15-20
	hFitAndre_pp_rap0016->SetBinContent(5, 0.543);//20-30

	hFitAndre_pp_rap0016->SetBinError(1, 0.002);//6.5-9
	hFitAndre_pp_rap0016->SetBinError(2, 0.002);//9-12
	hFitAndre_pp_rap0016->SetBinError(3, 0.003);//12-15
	hFitAndre_pp_rap0016->SetBinError(4, 0.004);//15-20
	hFitAndre_pp_rap0016->SetBinError(5, 0.007);//20-30

	hFitAndre_pp_rap1624->SetBinContent(1, 0.161);//3-6.5
	hFitAndre_pp_rap1624->SetBinContent(2, 0.236);//6.5-12
	hFitAndre_pp_rap1624->SetBinContent(3, 0.409);//12-30

	hFitAndre_pp_rap1624->SetBinError(1, 0.002);//3-6.5
	hFitAndre_pp_rap1624->SetBinError(2, 0.002);//6.5-12
	hFitAndre_pp_rap1624->SetBinError(3, 0.005);//12-30

	hFitMihee_pp_rap0016->SetBinContent(1, 0.24);//6.5-9
	hFitMihee_pp_rap0016->SetBinContent(2, 0.31);//9-12
	hFitMihee_pp_rap0016->SetBinContent(3, 0.39);//12-15
	hFitMihee_pp_rap0016->SetBinContent(4, 0.46);//15-20
	hFitMihee_pp_rap0016->SetBinContent(5, 0.54);//20-30

	hFitMihee_pp_rap1624->SetBinContent(1, 0.14);//3-6.5
	hFitMihee_pp_rap1624->SetBinContent(2, 0.24);//6.5-12
	hFitMihee_pp_rap1624->SetBinContent(3, 0.41);//12-30


	TCanvas*  c1 = new TCanvas("c1", "pp data Jpsi", 800, 600);

	TH1I* hframe1 = new TH1I("hframe1", "hframe1", 1, 0, 30); //just holder to set the axis properly
	hframe1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe1->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe1->GetYaxis()->SetTitle("non-prompt fraction");
	hframe1->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe1->Draw();
	hframe1->SetTitle("");
	
	hJpsiN_pp_rap0016->Draw("sameP");
	hJpsiN_pp_rap0016->SetLineColor(kBlue);
	hJpsiN_pp_rap0016->SetTitle("");
	hJpsiN_pp_rap0016->SetMarkerStyle(29);
	hJpsiN_pp_rap0016->SetMarkerColor(kBlue);

	hFitAndre_pp_rap0016->SetLineColor(kRed);
	hFitAndre_pp_rap0016->SetMarkerStyle(20);
	hFitAndre_pp_rap0016->SetMarkerColor(kRed);
	hFitAndre_pp_rap0016->Draw("sameP");

	//hFitMihee_pp_rap0016->SetLineColor(kGreen);
	//hFitMihee_pp_rap0016->SetMarkerStyle(21);
	//hFitMihee_pp_rap0016->SetMarkerColor(kGreen);
	//hFitMihee_pp_rap0016->Draw("sameP");


	hJpsiN_pp_rap1624->SetLineColor(kBlue);
	hJpsiN_pp_rap1624->SetMarkerStyle(30);
	hJpsiN_pp_rap1624->SetMarkerColor(kBlue);
	hJpsiN_pp_rap1624->Draw("sameP");

	hFitAndre_pp_rap1624->SetLineColor(kRed);
	hFitAndre_pp_rap1624->SetMarkerStyle(24);
	hFitAndre_pp_rap1624->SetMarkerColor(kRed);
	hFitAndre_pp_rap1624->Draw("sameP");

	//hFitMihee_pp_rap1624->SetLineColor(kGreen);
	//hFitMihee_pp_rap1624->SetMarkerStyle(25);
	//hFitMihee_pp_rap1624->SetMarkerColor(kGreen);
	//hFitMihee_pp_rap1624->Draw("sameP");

	TLegend*leg1 = new TLegend(0.42, 0.15, 0.85, 0.38, "");
	leg1->SetFillColor(kWhite);
	leg1->SetBorderSize(0);
	leg1->SetTextSize(0.035);
	leg1->AddEntry(hJpsiN_pp_rap0016, "pp 5.02 TeV, |y|<1.6", "p");
	leg1->AddEntry(hFitAndre_pp_rap0016, "pp 5.02 TeV fit Andre, |y|<1.6", "p");
	//leg1->AddEntry(hFitMihee_pp_rap0016, "pp 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg1->AddEntry(hJpsiN_pp_rap1624, "pp 5.02 TeV, 1.6<|y|<2.4", "p");
	leg1->AddEntry(hFitAndre_pp_rap1624, "pp 5.02 TeV fit Andre, 1.6<|y|<2.4", "p");
	//leg1->AddEntry(hFitMihee_pp_rap1624, "pp 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg1->Draw("same");

	c1->SaveAs("bfrac_plots/bfraction_1S_pp_pT.png");
	c1->SaveAs("bfrac_plots/bfraction_1S_pp_pT.pdf");


	// splitting the canvas 

	TCanvas*  c1_mid = new TCanvas("c1_mid", "pp data Jpsi", 800, 600);

	TH1I* hframe1_mid = new TH1I("hframe1_mid", "hframe1_mid", 1, 0, 30); //just holder to set the axis properly
	hframe1_mid->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe1_mid->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe1_mid->GetYaxis()->SetTitle("non-prompt fraction");
	hframe1_mid->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe1_mid->Draw();
	hframe1_mid->SetTitle("");


	double p8319_d21x1y1_xval[] = { 8.5, 9.5, 10.5, 11.5, 12.75, 14.25, 16.5, 24.0, 37.5, 57.5 };
	double p8319_d21x1y1_xerrminus[] = { 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 1.5, 6.0, 7.5, 12.5 };
	double p8319_d21x1y1_xerrplus[] = { 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 1.5, 6.0, 7.5, 12.5 };
	double p8319_d21x1y1_yval[] = { 0.271, 0.285, 0.32, 0.345, 0.373, 0.417, 0.454, 0.535, 0.633, 0.646 };
	double p8319_d21x1y1_yerrminus[] = { 0.00670820393249937, 0.006403124237432849, 0.006403124237432849, 0.008602325267042627, 0.0058309518948453, 0.0078102496759066544, 0.007211102550927979, 0.007211102550927979, 0.018027756377319945, 0.044944410108488465 };
	double p8319_d21x1y1_yerrplus[] = { 0.00670820393249937, 0.006403124237432849, 0.006403124237432849, 0.008602325267042627, 0.0058309518948453, 0.0078102496759066544, 0.007211102550927979, 0.007211102550927979, 0.018027756377319945,	0.044944410108488465 };
	double p8319_d21x1y1_ystatminus[] = { 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.015, 0.038 };
	double p8319_d21x1y1_ystatplus[] = { 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.015, 0.038 };
	int p8319_d21x1y1_numpoints = 10;
	
	TGraphAsymmErrors* p8319_d21x1y1 = new  TGraphAsymmErrors(p8319_d21x1y1_numpoints, p8319_d21x1y1_xval, p8319_d21x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d21x1y1_yerrminus, p8319_d21x1y1_yerrplus);
	p8319_d21x1y1->SetName("/HepData/8319/d21x1y1");
	p8319_d21x1y1->SetTitle("/HepData/8319/d21x1y1");

	double p8319_d22x1y1_xval[] = { 8.5, 9.5, 11.0, 13.5, 22.5, 37.5 };
	double p8319_d22x1y1_xerrminus[] = { 0.5, 0.5, 1.0, 1.5, 7.5, 7.5 };
	double p8319_d22x1y1_xerrplus[] = { 0.5, 0.5, 1.0, 1.5, 7.5, 7.5 };
	double p8319_d22x1y1_yval[] = { 0.265, 0.281, 0.324, 0.38, 0.481, 0.616 };
	double p8319_d22x1y1_yerrminus[] = { 0.009899494936611667, 0.009899494936611667, 0.009219544457292887, 0.01, 0.009899494936611667, 0.04031128874149275 };
	double p8319_d22x1y1_yerrplus[] = { 0.009899494936611667, 0.009899494936611667, 0.009219544457292887, 0.01, 0.009899494936611667, 0.04031128874149275 };
	double p8319_d22x1y1_ystatminus[] = { 0.007, 0.007, 0.006, 0.006, 0.007, 0.029 };
	double p8319_d22x1y1_ystatplus[] = { 0.007, 0.007, 0.006, 0.006, 0.007, 0.029 };
	int p8319_d22x1y1_numpoints = 6;
	TGraphAsymmErrors* p8319_d22x1y1 = new TGraphAsymmErrors(p8319_d22x1y1_numpoints, p8319_d22x1y1_xval, p8319_d22x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d22x1y1_yerrminus, p8319_d22x1y1_yerrplus);
	p8319_d22x1y1->SetName("/HepData/8319/d22x1y1");
	p8319_d22x1y1->SetTitle("/HepData/8319/d22x1y1");

	double p8319_d23x1y1_xval[] = { 6.75, 7.25, 7.75, 8.25, 8.75, 9.5, 10.5, 11.5, 13.5, 22.5, 37.5 };
	double p8319_d23x1y1_xerrminus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 1.5, 7.5, 7.5 };
	double p8319_d23x1y1_xerrplus[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 1.5, 7.5, 7.5 };
	double p8319_d23x1y1_yval[] = { 0.207, 0.221, 0.219, 0.238, 0.253, 0.275, 0.299, 0.335, 0.376, 0.466, 0.583 };
	double p8319_d23x1y1_yerrminus[] = { 0.0116619037896906, 0.010816653826391968, 0.013892443989449804, 0.013038404810405297, 0.011401754250991379, 0.01, 0.010816653826391968, 0.011401754250991379, 0.012529964086141668, 0.0147648230602334, 0.03538361202590826 };
	double p8319_d23x1y1_yerrplus[] = { 0.0116619037896906, 0.010816653826391968, 0.013892443989449804, 0.013038404810405297, 0.011401754250991379, 0.01, 0.010816653826391968, 0.011401754250991379, 0.012529964086141668,	0.0147648230602334, 0.03538361202590826 };
	double p8319_d23x1y1_ystatminus[] = { 0.01, 0.009, 0.007, 0.007, 0.007, 0.006, 0.006, 0.007, 0.006,	0.007, 0.026 };
	double p8319_d23x1y1_ystatplus[] = { 0.01, 0.009, 0.007, 0.007, 0.007, 0.006, 0.006, 0.007, 0.006,	0.007, 0.026 };
	int p8319_d23x1y1_numpoints = 11;
	TGraphAsymmErrors* p8319_d23x1y1 = new TGraphAsymmErrors(p8319_d23x1y1_numpoints, p8319_d23x1y1_xval, p8319_d23x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d23x1y1_yerrminus, p8319_d23x1y1_yerrplus);
	p8319_d23x1y1->SetName("/HepData/8319/d23x1y1");
	p8319_d23x1y1->SetTitle("/HepData/8319/d23x1y1");

	hJpsiN_pp_rap0016->Draw("sameP");
	hFitAndre_pp_rap0016->Draw("sameP");
	//hFitMihee_pp_rap0016->Draw("sameP");

	p8319_d21x1y1->SetLineColor(kMagenta);
	p8319_d21x1y1->SetMarkerStyle(30);
	p8319_d21x1y1->SetMarkerColor(kMagenta);

	p8319_d22x1y1->SetLineColor(kCyan);
	p8319_d22x1y1->SetMarkerStyle(24);
	p8319_d22x1y1->SetMarkerColor(kCyan);

	p8319_d23x1y1->SetLineColor(kOrange);
	p8319_d23x1y1->SetMarkerStyle(25);
	p8319_d23x1y1->SetMarkerColor(kOrange);

	p8319_d21x1y1->Draw("ZP");
	p8319_d22x1y1->Draw("ZP");
	p8319_d23x1y1->Draw("ZP");

	TLegend*leg1_mid = new TLegend(0.40, 0.15, 0.85, 0.38, "");
	leg1_mid->SetFillColor(kWhite);
	leg1_mid->SetBorderSize(0);
	leg1_mid->SetTextSize(0.035);
	leg1_mid->AddEntry(hJpsiN_pp_rap0016, "pp 5.02 TeV, |y|<1.6", "p");
	leg1_mid->AddEntry(hFitAndre_pp_rap0016, "pp 5.02 TeV fit Andre, |y|<1.6", "p");
	//leg1_mid->AddEntry(hFitMihee_pp_rap0016, "pp 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg1_mid->AddEntry(p8319_d21x1y1, "pp 7 TeV, BPH 10-014, |y|<0.9", "p");
	leg1_mid->AddEntry(p8319_d22x1y1, "pp 7 TeV, BPH 10-014, 0.9<|y|<1.2", "p");
	leg1_mid->AddEntry(p8319_d23x1y1, "pp 7 TeV, BPH 10-014, 1.2<|y|<1.6", "p");
	leg1_mid->Draw("same");//*/

	c1_mid->SaveAs("bfrac_plots/bfraction_1S_pp_pT_mid.png");
	c1_mid->SaveAs("bfrac_plots/bfraction_1S_pp_pT_mid.pdf");



	TCanvas*  c1_fwd = new TCanvas("c1_fwd", "pp data Jpsi", 800, 600);

	TH1I* hframe1_fwd = new TH1I("hframe1_fwd", "hframe1_fwd", 1, 0, 30); //just holder to set the axis properly
	hframe1_fwd->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe1_fwd->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe1_fwd->GetYaxis()->SetTitle("non-prompt fraction");
	hframe1_fwd->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe1_fwd->Draw();
	hframe1_fwd->SetTitle("");


	double p8319_d24x1y1_xval[] = { 6.75, 7.125, 7.375, 7.75, 8.25, 8.75, 9.5, 10.5, 11.5, 13.5, 22.5 };
	double p8319_d24x1y1_xerrminus[] = { 0.25, 0.125, 0.125, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 1.5, 7.5 };
	double p8319_d24x1y1_xerrplus[] = { 0.25, 0.125, 0.125, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 1.5, 7.5 };
	double p8319_d24x1y1_yval[] = { 0.199, 0.225, 0.239, 0.233, 0.237, 0.25, 0.258, 0.299, 0.321, 0.356, 0.458 };
	double p8319_d24x1y1_yerrminus[] = { 0.015620499351813309, 0.016401219466856725, 0.016401219466856725, 0.015811388300841896, 0.017492855684535902, 0.016643316977093238, 0.01565247584249853, 0.012806248474865698, 0.018027756377319945, 0.018788294228055936, 0.02102379604162864 };
	double p8319_d24x1y1_yerrplus[] = { 0.015620499351813309, 0.016401219466856725, 0.016401219466856725, 0.015811388300841896, 0.017492855684535902, 0.016643316977093238, 0.01565247584249853, 0.012806248474865698, 0.018027756377319945, 0.018788294228055936, 0.02102379604162864 };
	double p8319_d24x1y1_ystatminus[] = { 0.01, 0.013, 0.013, 0.009, 0.009, 0.009, 0.007, 0.008, 0.01, 0.008, 0.009 };
	double p8319_d24x1y1_ystatplus[] = { 0.01, 0.013, 0.013, 0.009, 0.009, 0.009, 0.007, 0.008, 0.01, 0.008, 0.009 };
	int p8319_d24x1y1_numpoints = 11;
	TGraphAsymmErrors* p8319_d24x1y1 = new TGraphAsymmErrors(p8319_d24x1y1_numpoints, p8319_d24x1y1_xval, p8319_d24x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d24x1y1_yerrminus, p8319_d24x1y1_yerrplus);
	p8319_d24x1y1->SetName("/HepData/8319/d24x1y1");
	p8319_d24x1y1->SetTitle("/HepData/8319/d24x1y1");

	double p8319_d25x1y1_xval[] = { 6.75, 8.5, 9.5, 11.0, 13.5, 22.5 };
	double p8319_d25x1y1_xerrminus[] = { 1.25, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d25x1y1_xerrplus[] = { 1.25, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d25x1y1_yval[] = { 0.214, 0.226, 0.275, 0.294, 0.368, 0.432 };
	double p8319_d25x1y1_yerrminus[] = { 0.01769180601295413, 0.025612496949731396, 0.02690724809414742, 0.022803508501982758, 0.023345235059857503, 0.03047950130825634 };
	double p8319_d25x1y1_yerrplus[] = { 0.01769180601295413, 0.025612496949731396, 0.02690724809414742, 0.022803508501982758, 0.023345235059857503, 0.03047950130825634 };
	double p8319_d25x1y1_ystatminus[] = { 0.012, 0.016, 0.018, 0.014, 0.016, 0.02 };
	double p8319_d25x1y1_ystatplus[] = { 0.012, 0.016, 0.018, 0.014, 0.016, 0.02 };
	int p8319_d25x1y1_numpoints = 6;
	TGraphAsymmErrors* p8319_d25x1y1 = new TGraphAsymmErrors(p8319_d25x1y1_numpoints, p8319_d25x1y1_xval, p8319_d25x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d25x1y1_yerrminus, p8319_d25x1y1_yerrplus);
	p8319_d25x1y1->SetName("/HepData/8319/d25x1y1");
	p8319_d25x1y1->SetTitle("/HepData/8319/d25x1y1");




	hJpsiN_pp_rap1624->SetTitle("");
	hJpsiN_pp_rap1624->Draw("sameP");
	hFitAndre_pp_rap1624->Draw("sameP");
	//hFitMihee_pp_rap1624->Draw("sameP");

	p8319_d24x1y1->SetLineColor(kMagenta);
	p8319_d24x1y1->SetMarkerStyle(29);
	p8319_d24x1y1->SetMarkerColor(kMagenta);

	p8319_d25x1y1->SetLineColor(kCyan);
	p8319_d25x1y1->SetMarkerStyle(20);
	p8319_d25x1y1->SetMarkerColor(kCyan);

	p8319_d24x1y1->Draw("ZP");
	p8319_d25x1y1->Draw("ZP");

	TLegend*leg1_fwd = new TLegend(0.42, 0.15, 0.85, 0.38, "");
	leg1_fwd->SetFillColor(kWhite);
	leg1_fwd->SetBorderSize(0);
	leg1_fwd->SetTextSize(0.035);
	leg1_fwd->AddEntry(hJpsiN_pp_rap1624, "pp 5.02 TeV, 1.6<|y|<2.4", "p");
	leg1_fwd->AddEntry(hFitAndre_pp_rap1624, "pp 5.02 TeV fit Andre, 1.6<|y|<2.4", "p");
	//leg1_fwd->AddEntry(hFitMihee_pp_rap1624, "pp 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg1_fwd->AddEntry(p8319_d24x1y1, "pp 7 TeV, BPH 10-014, 1.6<|y|<2.1", "p");
	leg1_fwd->AddEntry(p8319_d25x1y1, "pp 7 TeV, BPH 10-014, 2.1<|y|<2.4", "p");
	leg1_fwd->Draw("same");

	c1_fwd->SaveAs("bfrac_plots/bfraction_1S_pp_pT_fwd.png");
	c1_fwd->SaveAs("bfrac_plots/bfraction_1S_pp_pT_fwd.pdf");






	// 2S (pp)

	TH1D* h2sJpsiN_pp_rap0016 = new TH1D("h2sJpsiN_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* h2sJpsiN_pp_rap1624 = new TH1D("h2sJpsiN_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* h2sFitAndre_pp_rap0016 = new TH1D("h2sFitAndre_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* h2sFitAndre_pp_rap1624 = new TH1D("h2sFitAndre_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* h2sFitMihee_pp_rap0016 = new TH1D("h2sFitMihee_pp_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* h2sFitMihee_pp_rap1624 = new TH1D("h2sFitMihee_pp_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);


	// filling pp 2S histograms (by hand)

	h2sJpsiN_pp_rap0016->SetBinContent(1, 0.327);//6.5-9
	h2sJpsiN_pp_rap0016->SetBinContent(2, 0.354);//9-12
	h2sJpsiN_pp_rap0016->SetBinContent(3, 0.437);//12-15
	h2sJpsiN_pp_rap0016->SetBinContent(4, 0.531);//15-20
	h2sJpsiN_pp_rap0016->SetBinContent(5, 0.563);//20-30

	h2sJpsiN_pp_rap0016->SetBinError(1, 0.014);//6.5-9
	h2sJpsiN_pp_rap0016->SetBinError(2, 0.013);//9-12
	h2sJpsiN_pp_rap0016->SetBinError(3, 0.020);//12-15
	h2sJpsiN_pp_rap0016->SetBinError(4, 0.027);//15-20
	h2sJpsiN_pp_rap0016->SetBinError(5, 0.042);//20-30

	h2sJpsiN_pp_rap1624->SetBinContent(1, 0.231);//3-6.5
	h2sJpsiN_pp_rap1624->SetBinContent(2, 0.333);//6.5-12
	h2sJpsiN_pp_rap1624->SetBinContent(3, 0.438);//12-30

	h2sJpsiN_pp_rap1624->SetBinError(1, 0.018);//3-6.5
	h2sJpsiN_pp_rap1624->SetBinError(2, 0.016);//6.5-12
	h2sJpsiN_pp_rap1624->SetBinError(3, 0.035);//12-30


	h2sFitAndre_pp_rap0016->SetBinContent(1, 0.339);//6.5-9
	h2sFitAndre_pp_rap0016->SetBinContent(2, 0.371);//9-12
	h2sFitAndre_pp_rap0016->SetBinContent(3, 0.450);//12-15
	h2sFitAndre_pp_rap0016->SetBinContent(4, 0.540);//15-20
	h2sFitAndre_pp_rap0016->SetBinContent(5, 0.630);//20-30

	h2sFitAndre_pp_rap0016->SetBinError(1, 0.011);//6.5-9
	h2sFitAndre_pp_rap0016->SetBinError(2, 0.011);//9-12
	h2sFitAndre_pp_rap0016->SetBinError(3, 0.017);//12-15
	h2sFitAndre_pp_rap0016->SetBinError(4, 0.022);//15-20
	h2sFitAndre_pp_rap0016->SetBinError(5, 0.035);//20-30

	h2sFitAndre_pp_rap1624->SetBinContent(1, 0.228);//3-6.5
	h2sFitAndre_pp_rap1624->SetBinContent(2, 0.330);//6.5-12
	h2sFitAndre_pp_rap1624->SetBinContent(3, 0.453);//12-30

	h2sFitAndre_pp_rap1624->SetBinError(1, 0.014);//3-6.5
	h2sFitAndre_pp_rap1624->SetBinError(2, 0.012);//6.5-12
	h2sFitAndre_pp_rap1624->SetBinError(3, 0.028);//12-30

	h2sFitMihee_pp_rap0016->SetBinContent(1, 0.59);//6.5-9
	h2sFitMihee_pp_rap0016->SetBinContent(2, 0.44);//9-12
	h2sFitMihee_pp_rap0016->SetBinContent(3, 0.50);//12-15
	h2sFitMihee_pp_rap0016->SetBinContent(4, 0.63);//15-20
	h2sFitMihee_pp_rap0016->SetBinContent(5, 0.79);//20-30
	h2sFitMihee_pp_rap0016->SetBinError(1, 0.0);//6.5-9   // NONEXISTENT ERROR
	h2sFitMihee_pp_rap0016->SetBinError(2, 0.0);//9-12   // NONEXISTENT ERROR
	h2sFitMihee_pp_rap0016->SetBinError(3, 0.03);//12-15
	h2sFitMihee_pp_rap0016->SetBinError(4, 0.03);//15-20
	h2sFitMihee_pp_rap0016->SetBinError(5, 0.0);//20-30   // NONEXISTENT ERROR

	h2sFitMihee_pp_rap1624->SetBinContent(1, 1.00);//3-6.5
	h2sFitMihee_pp_rap1624->SetBinContent(2, 0.62);//6.5-12
	h2sFitMihee_pp_rap1624->SetBinContent(3, 0.72);//12-30
	//h2sFitMihee_pp_rap1624->SetBinError(1, 0.14);//3-6.5
	//h2sFitMihee_pp_rap1624->SetBinError(2, 0.24);//6.5-12
	//h2sFitMihee_pp_rap1624->SetBinError(3, 0.41);//12-30


	TCanvas*  c2 = new TCanvas("c2", "pp data2S", 800, 600);

	TH1I* hframe2 = new TH1I("hframe2", "hframe2", 1, 0, 30); //just holder to set the axis properly
	hframe2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe2->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe2->GetYaxis()->SetTitle("non-prompt fraction");
	hframe2->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe2->Draw();
	hframe2->SetTitle("");

	h2sJpsiN_pp_rap0016->Draw("sameP");
	h2sJpsiN_pp_rap0016->SetLineColor(kBlue);
	h2sJpsiN_pp_rap0016->SetTitle("");
	h2sJpsiN_pp_rap0016->SetMarkerStyle(29);
	h2sJpsiN_pp_rap0016->SetMarkerColor(kBlue);

	h2sFitAndre_pp_rap0016->SetLineColor(kRed);
	h2sFitAndre_pp_rap0016->SetMarkerStyle(20);
	h2sFitAndre_pp_rap0016->SetMarkerColor(kRed);
	h2sFitAndre_pp_rap0016->Draw("sameP");

	//h2sFitMihee_pp_rap0016->SetLineColor(kGreen);
	//h2sFitMihee_pp_rap0016->SetMarkerStyle(21);
	//h2sFitMihee_pp_rap0016->SetMarkerColor(kGreen);
	//h2sFitMihee_pp_rap0016->Draw("sameP");


	h2sJpsiN_pp_rap1624->SetLineColor(kBlue);
	h2sJpsiN_pp_rap1624->SetMarkerStyle(30);
	h2sJpsiN_pp_rap1624->SetMarkerColor(kBlue);
	h2sJpsiN_pp_rap1624->Draw("sameP");

	h2sFitAndre_pp_rap1624->SetLineColor(kRed);
	h2sFitAndre_pp_rap1624->SetMarkerStyle(24);
	h2sFitAndre_pp_rap1624->SetMarkerColor(kRed);
	h2sFitAndre_pp_rap1624->Draw("sameP");

	//h2sFitMihee_pp_rap1624->SetLineColor(kGreen);
	//h2sFitMihee_pp_rap1624->SetMarkerStyle(25);
	//h2sFitMihee_pp_rap1624->SetMarkerColor(kGreen);
	//h2sFitMihee_pp_rap1624->Draw("sameP");

	TLegend*leg2 = new TLegend(0.38, 0.15, 0.85, 0.38, "");
	leg2->SetFillColor(kWhite);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.035);
	leg2->AddEntry(h2sJpsiN_pp_rap0016, "2S pp 5.02 TeV, |y|<1.6", "p");
	leg2->AddEntry(h2sFitAndre_pp_rap0016, "2S pp 5.02 TeV fit Andre, |y|<1.6", "p");
	//leg2->AddEntry(h2sFitMihee_pp_rap0016, "2S pp 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg2->AddEntry(h2sJpsiN_pp_rap1624, "2S pp 5.02 TeV, 1.6<|y|<2.4", "p");
	leg2->AddEntry(h2sFitAndre_pp_rap1624, "2S pp 5.02 TeV fit Andre, 1.6<|y|<2.4", "p");
	//leg2->AddEntry(h2sFitMihee_pp_rap1624, "2S pp 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg2->Draw("same");

	c2->SaveAs("bfrac_plots/bfraction_2S_pp_pT.png");
	c2->SaveAs("bfrac_plots/bfraction_2S_pp_pT.pdf");


	// splitting the canvas 
	
	TCanvas*  c2_mid = new TCanvas("c2_mid", "pp data 2S", 800, 600);

	TH1I* hframe2_mid = new TH1I("hframe2_mid", "hframe2_mid", 1, 0, 30); //just holder to set the axis properly
	hframe2_mid->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe2_mid->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe2_mid->GetYaxis()->SetTitle("non-prompt fraction");
	hframe2_mid->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe2_mid->Draw();
	hframe2_mid->SetTitle("");

	double p8319_d26x1y1_xval[] = { 7.25, 8.5, 9.5, 10.5, 11.5, 12.75, 14.25, 16.5, 24.0 };
	double p8319_d26x1y1_xerrminus[] = { 0.75, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 1.5, 6.0 };
	double p8319_d26x1y1_xerrplus[] = { 0.75, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 1.5, 6.0 };
	double p8319_d26x1y1_yval[] = { 0.326, 0.35, 0.335, 0.378, 0.394, 0.431, 0.456, 0.487, 0.572 };
	double p8319_d26x1y1_yerrminus[] = { 0.04244997055358225, 0.028653097563788806, 0.029, 0.0318276609256791, 0.03687817782917155, 0.02973213749463701, 0.03753664875824692, 0.029410882339705485, 0.03753664875824692 };
	double p8319_d26x1y1_yerrplus[] = { 0.04244997055358225, 0.028653097563788806, 0.029, 0.0318276609256791, 0.03687817782917155, 0.02973213749463701, 0.03753664875824692, 0.029410882339705485, 0.03753664875824692 };
	double p8319_d26x1y1_ystatminus[] = { 0.041, 0.025, 0.021, 0.023, 0.024, 0.022, 0.025, 0.024, 0.025 };
	double p8319_d26x1y1_ystatplus[] = { 0.041, 0.025, 0.021, 0.023, 0.024, 0.022, 0.025, 0.024, 0.025 };
	int p8319_d26x1y1_numpoints = 9;
	TGraphAsymmErrors* p8319_d26x1y1 = new TGraphAsymmErrors(p8319_d26x1y1_numpoints, p8319_d26x1y1_xval, p8319_d26x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d26x1y1_yerrminus, p8319_d26x1y1_yerrplus);
	p8319_d26x1y1->SetName("/HepData/8319/d26x1y1");
	p8319_d26x1y1->SetTitle("/HepData/8319/d26x1y1");

	double p8319_d27x1y1_xval[] = { 6.0, 7.25, 8.5, 9.5, 11.0, 13.5, 22.5 };
	double p8319_d27x1y1_xerrminus[] = { 0.5, 0.75, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d27x1y1_xerrplus[] = { 0.5, 0.75, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d27x1y1_yval[] = { 0.246, 0.308, 0.309, 0.32, 0.366, 0.375, 0.51 };
	double p8319_d27x1y1_yerrminus[] = { 0.09741663102366044, 0.03162277660168379, 0.041340053217188776, 0.05, 0.04031128874149275, 0.0425205832509386, 0.05 };
	double p8319_d27x1y1_yerrplus[] = { 0.09741663102366044, 0.03162277660168379, 0.041340053217188776, 0.05, 0.04031128874149275, 0.0425205832509386, 0.05 };
	double p8319_d27x1y1_ystatminus[] = { 0.093, 0.03, 0.035, 0.03, 0.029, 0.032, 0.03 };
	double p8319_d27x1y1_ystatplus[] = { 0.093, 0.03, 0.035, 0.03, 0.029, 0.032, 0.03 };
	int p8319_d27x1y1_numpoints = 7;
	TGraphAsymmErrors* p8319_d27x1y1 = new TGraphAsymmErrors(p8319_d27x1y1_numpoints, p8319_d27x1y1_xval, p8319_d27x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d27x1y1_yerrminus, p8319_d27x1y1_yerrplus);
	p8319_d27x1y1->SetName("/HepData/8319/d27x1y1");
	p8319_d27x1y1->SetTitle("/HepData/8319/d27x1y1");


	h2sJpsiN_pp_rap0016->Draw("sameP");
	h2sFitAndre_pp_rap0016->Draw("sameP");
	//h2sFitMihee_pp_rap0016->Draw("sameP");

	p8319_d26x1y1->SetLineColor(kMagenta);
	p8319_d26x1y1->SetMarkerStyle(30);
	p8319_d26x1y1->SetMarkerColor(kMagenta);

	p8319_d27x1y1->SetLineColor(kCyan);
	p8319_d27x1y1->SetMarkerStyle(24);
	p8319_d27x1y1->SetMarkerColor(kCyan);

	p8319_d26x1y1->Draw("sameZP");
	p8319_d27x1y1->Draw("sameZP");

	TLegend*leg2_mid = new TLegend(0.35, 0.15, 0.85, 0.38, "");
	leg2_mid->SetFillColor(kWhite);
	leg2_mid->SetBorderSize(0);
	leg2_mid->SetTextSize(0.035);
	leg2_mid->AddEntry(h2sJpsiN_pp_rap0016, "2S pp 5.02 TeV, |y|<1.6", "p");
	leg2_mid->AddEntry(h2sFitAndre_pp_rap0016, "2S pp 5.02 TeV fit Andre, |y|<1.6", "p");
	//leg2_mid->AddEntry(h2sFitMihee_pp_rap0016, "2S pp 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg2_mid->AddEntry(p8319_d26x1y1, "2S pp 7 TeV, BPH 10-014, |y|<1.2", "p");
	leg2_mid->AddEntry(p8319_d27x1y1, "2S pp 7 TeV, BPH 10-014, 1.2<|y|<1.6", "p");
	leg2_mid->Draw("same");

	c2_mid->SaveAs("bfrac_plots/bfraction_2S_pp_pT_mid.png");
	c2_mid->SaveAs("bfrac_plots/bfraction_2S_pp_pT_mid.pdf");



	TCanvas*  c2_fwd = new TCanvas("c2_fwd", "pp data 2S", 800, 600);

	TH1I* hframe2_fwd = new TH1I("hframe2_fwd", "hframe2_fwd", 1, 0, 30); //just holder to set the axis properly
	hframe2_fwd->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe2_fwd->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe2_fwd->GetYaxis()->SetTitle("non-prompt fraction");
	hframe2_fwd->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe2_fwd->Draw();
	hframe2_fwd->SetTitle("");

	double p8319_d28x1y1_xval[] = { 6.0, 7.25, 8.5, 9.5, 11.0, 13.5, 22.5 };
	double p8319_d28x1y1_xerrminus[] = { 0.5, 0.75, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d28x1y1_xerrplus[] = { 0.5, 0.75, 0.5, 0.5, 1.0, 1.5, 7.5 };
	double p8319_d28x1y1_yval[] = { 0.21, 0.19, 0.3, 0.3, 0.37, 0.372, 0.45 };
	double p8319_d28x1y1_yerrminus[] = { 0.09219544457292887, 0.042426406871192854, 0.0565685424949238, 0.06403124237432849, 0.05, 0.03807886552931954, 0.0565685424949238 };
	double p8319_d28x1y1_yerrplus[] = { 0.09219544457292887, 0.042426406871192854, 0.0565685424949238, 0.06403124237432849, 0.05, 0.03807886552931954, 0.0565685424949238 };
	double p8319_d28x1y1_ystatminus[] = { 0.06, 0.03, 0.04, 0.04, 0.03, 0.035, 0.04 };
	double p8319_d28x1y1_ystatplus[] = { 0.06, 0.03, 0.04, 0.04, 0.03, 0.035, 0.04 };
	int p8319_d28x1y1_numpoints = 7;
	TGraphAsymmErrors* p8319_d28x1y1 = new TGraphAsymmErrors(p8319_d28x1y1_numpoints, p8319_d28x1y1_xval, p8319_d28x1y1_yval, zero_xerrminus, zero_xerrplus, p8319_d28x1y1_yerrminus, p8319_d28x1y1_yerrplus);
	p8319_d28x1y1->SetName("/HepData/8319/d28x1y1");
	p8319_d28x1y1->SetTitle("/HepData/8319/d28x1y1");


	h2sJpsiN_pp_rap1624->Draw("sameP");
	h2sFitAndre_pp_rap1624->Draw("sameP");
	//h2sFitMihee_pp_rap1624->Draw("sameP");

	p8319_d28x1y1->SetLineColor(kMagenta);
	p8319_d28x1y1->SetMarkerStyle(29);
	p8319_d28x1y1->SetMarkerColor(kMagenta);

	p8319_d28x1y1->Draw("sameZP");

	TLegend*leg2_fwd = new TLegend(0.15, 0.65, 0.55, 0.88, "");
	leg2_fwd->SetFillColor(kWhite);
	leg2_fwd->SetBorderSize(0);
	leg2_fwd->SetTextSize(0.035);
	leg2_fwd->AddEntry(h2sJpsiN_pp_rap1624, "2S pp 5.02 TeV, 1.6<|y|<2.4", "p");
	leg2_fwd->AddEntry(h2sFitAndre_pp_rap1624, "2S pp 5.02 TeV fit Andre, 1.6<|y|<2.4", "p");
	//leg2_fwd->AddEntry(h2sFitMihee_pp_rap1624, "2S pp 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg2_fwd->AddEntry(p8319_d28x1y1, "2S pp 7 TeV, BPH 10-014, 1.6<|y|<2.4", "p");
	leg2_fwd->Draw("same");

	c2_fwd->SaveAs("bfrac_plots/bfraction_2S_pp_pT_fwd.png");
	c2_fwd->SaveAs("bfrac_plots/bfraction_2S_pp_pT_fwd.pdf");
	

	/////////////////////////////////////////////////

	/////////              PbPb plots          ////////

	//////////////////////////////////////////////////


	//1S

	//centrality
	TH1D* hJpsiN_cent_rap0016 = new TH1D("hJpsiN_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	TH1D* hJpsiN_cent_rap1624 = new TH1D("hJpsiN_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);

	TH1D* hFitAndre_cent_rap0016 = new TH1D("hFitAndre_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	TH1D* hFitAndre_cent_rap1624 = new TH1D("hFitAndre_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);

	//TH1D* hFitMihee_cent_rap0016 = new TH1D("hFitMihee_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	//TH1D* hFitMihee_cent_rap1624 = new TH1D("hFitMihee_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);


	//pt
	TH1D* hJpsiN_pbpb_rap0016 = new TH1D("hJpsiN_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* hJpsiN_pbpb_rap1624 = new TH1D("hJpsiN_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* hFitAndre_pbpb_rap0016 = new TH1D("hFitAndre_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* hFitAndre_pbpb_rap1624 = new TH1D("hFitAndre_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	//TH1D* hFitMihee_pbpb_rap0016 = new TH1D("hFitMihee_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	//TH1D* hFitMihee_pbpb_rap1624 = new TH1D("hFitMihee_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);


	// filling cent histograms (by hand)

	hJpsiN_cent_rap0016->SetBinContent(1, 0.405);//0-10
	hJpsiN_cent_rap0016->SetBinContent(2, 0.370);//10-20
	hJpsiN_cent_rap0016->SetBinContent(3, 0.370);//20-30
	hJpsiN_cent_rap0016->SetBinContent(4, 0.347);//30-40 // !!
	hJpsiN_cent_rap0016->SetBinContent(5, 0.320);//40-50 // !!
	hJpsiN_cent_rap0016->SetBinContent(6, 0.305);//50-100 // !!

	hJpsiN_cent_rap0016->SetBinError(1, 0.009);//0-10
	hJpsiN_cent_rap0016->SetBinError(2, 0.009);//10-20
	hJpsiN_cent_rap0016->SetBinError(3, 0.010);//20-30
	hJpsiN_cent_rap0016->SetBinError(4, 0.011);//30-40
	hJpsiN_cent_rap0016->SetBinError(5, 0.013);//40-50
	hJpsiN_cent_rap0016->SetBinError(6, 0.013);//50-100

	hJpsiN_cent_rap1624->SetBinContent(1, 0.296);//0-20
	hJpsiN_cent_rap1624->SetBinContent(2, 0.264);//20-40
	hJpsiN_cent_rap1624->SetBinContent(3, 0.255);//40-100 // !!

	hJpsiN_cent_rap1624->SetBinError(1, 0.013);//0-20
	hJpsiN_cent_rap1624->SetBinError(2, 0.014);//20-40
	hJpsiN_cent_rap1624->SetBinError(3, 0.014);//40-100


	hFitAndre_cent_rap0016->SetBinContent(1, 0.399);//0-10
	hFitAndre_cent_rap0016->SetBinContent(2, 0.370);//10-20
	hFitAndre_cent_rap0016->SetBinContent(3, 0.356);//20-30
	hFitAndre_cent_rap0016->SetBinContent(4, 0.362);//30-40
	hFitAndre_cent_rap0016->SetBinContent(5, 0.311);//40-50
	hFitAndre_cent_rap0016->SetBinContent(6, 0.296);//50-100

	hFitAndre_cent_rap0016->SetBinError(1, 0.007);//0-10
	hFitAndre_cent_rap0016->SetBinError(2, 0.007);//10-20
	hFitAndre_cent_rap0016->SetBinError(3, 0.008);//20-30
	hFitAndre_cent_rap0016->SetBinError(4, 0.010);//30-40
	hFitAndre_cent_rap0016->SetBinError(5, 0.011);//40-50
	hFitAndre_cent_rap0016->SetBinError(6, 0.011);//50-100

	hFitAndre_cent_rap1624->SetBinContent(1, 0.314);//0-20
	hFitAndre_cent_rap1624->SetBinContent(2, 0.261);//20-40
	hFitAndre_cent_rap1624->SetBinContent(3, 0.262);//40-100

	hFitAndre_cent_rap1624->SetBinError(1, 0.007);//0-20
	hFitAndre_cent_rap1624->SetBinError(2, 0.008);//20-40
	hFitAndre_cent_rap1624->SetBinError(3, 0.009);//40-100


	TCanvas*  c3 = new TCanvas("c3", "PbPb cent data Jpsi", 800, 600);

	hJpsiN_cent_rap0016->GetXaxis()->SetTitle("Centrality bin");
	hJpsiN_cent_rap0016->GetXaxis()->SetRangeUser(0, 100);
	hJpsiN_cent_rap0016->GetYaxis()->SetTitle("non-prompt fraction");
	hJpsiN_cent_rap0016->GetYaxis()->SetRangeUser(0.0, 0.7);
	hJpsiN_cent_rap0016->Draw("P");
	hJpsiN_cent_rap0016->SetLineColor(kBlue);
	hJpsiN_cent_rap0016->SetTitle("");
	hJpsiN_cent_rap0016->SetMarkerStyle(29);
	hJpsiN_cent_rap0016->SetMarkerColor(kBlue);

	hFitAndre_cent_rap0016->SetLineColor(kRed);
	hFitAndre_cent_rap0016->SetMarkerStyle(20);
	hFitAndre_cent_rap0016->SetMarkerColor(kRed);
	hFitAndre_cent_rap0016->Draw("sameP");

	/*hFitMihee_cent_rap0016->SetLineColor(kGreen);
	hFitMihee_cent_rap0016->SetMarkerStyle(21);
	hFitMihee_cent_rap0016->SetMarkerColor(kGreen);
	hFitMihee_cent_rap0016->Draw("sameP");*/


	hJpsiN_cent_rap1624->SetLineColor(kBlue);
	hJpsiN_cent_rap1624->SetMarkerStyle(30);
	hJpsiN_cent_rap1624->SetMarkerColor(kBlue);
	hJpsiN_cent_rap1624->Draw("sameP");

	hFitAndre_cent_rap1624->SetLineColor(kRed);
	hFitAndre_cent_rap1624->SetMarkerStyle(24);
	hFitAndre_cent_rap1624->SetMarkerColor(kRed);
	hFitAndre_cent_rap1624->Draw("sameP");

	/*hFitMihee_cent_rap1624->SetLineColor(kGreen);
	hFitMihee_cent_rap1624->SetMarkerStyle(25);
	hFitMihee_cent_rap1624->SetMarkerColor(kGreen);
	hFitMihee_cent_rap1624->Draw("sameP");*/

	TLegend*leg3 = new TLegend(0.40, 0.65, 0.80, 0.88, "");
	leg3->SetFillColor(kWhite);
	leg3->SetBorderSize(0);
	leg3->SetTextSize(0.035);
	leg3->AddEntry(hJpsiN_cent_rap0016, "PbPb 6.5<p_{T}<30, |y|<1.6", "p");
	leg3->AddEntry(hFitAndre_cent_rap0016, "PbPb 6.5<p_{T}<30 fit Andre, |y|<1.6", "p");
	//leg3->AddEntry(hFitMihee_cent_rap0016, "PbPb 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg3->AddEntry(hJpsiN_cent_rap1624, "PbPb 3<p_{T}<30, 1.6<|y|<2.4", "p");
	leg3->AddEntry(hFitAndre_cent_rap1624, "PbPb 3<p_{T}<30 fit Andre, 1.6<|y|<2.4", "p");
	//leg3->AddEntry(hFitMihee_cent_rap1624, "PbPb 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg3->Draw("same");

	c3->SaveAs("bfrac_plots/bfraction_1S_pbpb_cent.png");
	c3->SaveAs("bfrac_plots/bfraction_1S_pbpb_cent.pdf");



	// filling pbpb histograms pt

	hJpsiN_pbpb_rap0016->SetBinContent(1, 0.299);//6.5-9
	hJpsiN_pbpb_rap0016->SetBinContent(2, 0.363);//9-12
	hJpsiN_pbpb_rap0016->SetBinContent(3, 0.417);//12-15
	hJpsiN_pbpb_rap0016->SetBinContent(4, 0.498);//15-20
	hJpsiN_pbpb_rap0016->SetBinContent(5, 0.521);//20-30

	hJpsiN_pbpb_rap0016->SetBinError(1, 0.009);//6.5-9
	hJpsiN_pbpb_rap0016->SetBinError(2, 0.007);//9-12
	hJpsiN_pbpb_rap0016->SetBinError(3, 0.009);//12-15
	hJpsiN_pbpb_rap0016->SetBinError(4, 0.012);//15-20
	hJpsiN_pbpb_rap0016->SetBinError(5, 0.018);//20-30


	hJpsiN_pbpb_rap1624->SetBinContent(1, 0.230);//3-6.5
	hJpsiN_pbpb_rap1624->SetBinContent(2, 0.306);//6.5-12
	hJpsiN_pbpb_rap1624->SetBinContent(3, 0.441);//12-30

	hJpsiN_pbpb_rap1624->SetBinError(1, 0.019);//3-6.5
	hJpsiN_pbpb_rap1624->SetBinError(2, 0.009);//6.5-12
	hJpsiN_pbpb_rap1624->SetBinError(3, 0.015);//12-30

	hFitAndre_pbpb_rap0016->SetBinContent(1, 0.296);//6.5-9
	hFitAndre_pbpb_rap0016->SetBinContent(2, 0.366);//9-12
	hFitAndre_pbpb_rap0016->SetBinContent(3, 0.410);//12-15
	hFitAndre_pbpb_rap0016->SetBinContent(4, 0.496);//15-20
	hFitAndre_pbpb_rap0016->SetBinContent(5, 0.501);//20-30

	hFitAndre_pbpb_rap0016->SetBinError(1, 0.006);//6.5-9
	hFitAndre_pbpb_rap0016->SetBinError(2, 0.006);//9-12
	hFitAndre_pbpb_rap0016->SetBinError(3, 0.009);//12-15
	hFitAndre_pbpb_rap0016->SetBinError(4, 0.011);//15-20
	hFitAndre_pbpb_rap0016->SetBinError(5, 0.016);//20-30

	hFitAndre_pbpb_rap1624->SetBinContent(1, 0.222);//3-6.5
	hFitAndre_pbpb_rap1624->SetBinContent(2, 0.294);//6.5-12
	hFitAndre_pbpb_rap1624->SetBinContent(3, 0.433);//12-30

	hFitAndre_pbpb_rap1624->SetBinError(1, 0.007);//3-6.5
	hFitAndre_pbpb_rap1624->SetBinError(2, 0.006);//6.5-12
	hFitAndre_pbpb_rap1624->SetBinError(3, 0.014);//12-30




	TCanvas*  c4 = new TCanvas("c4", "pbpb data Jpsi", 800, 600);

	TH1I* hframe4 = new TH1I("hframe4", "hframe4", 1, 0, 30); //just holder to set the axis properly
	hframe4->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe4->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe4->GetYaxis()->SetTitle("non-prompt fraction");
	hframe4->GetYaxis()->SetRangeUser(0.0, 0.7);
	hframe4->Draw();
	hframe4->SetTitle("");

	hJpsiN_pbpb_rap0016->Draw("sameP");
	hJpsiN_pbpb_rap0016->SetLineColor(kBlue);
	hJpsiN_pbpb_rap0016->SetTitle("");
	hJpsiN_pbpb_rap0016->SetMarkerStyle(29);
	hJpsiN_pbpb_rap0016->SetMarkerColor(kBlue);

	hFitAndre_pbpb_rap0016->SetLineColor(kRed);
	hFitAndre_pbpb_rap0016->SetMarkerStyle(20);
	hFitAndre_pbpb_rap0016->SetMarkerColor(kRed);
	hFitAndre_pbpb_rap0016->Draw("sameP");

	hJpsiN_pbpb_rap1624->SetLineColor(kBlue);
	hJpsiN_pbpb_rap1624->SetMarkerStyle(30);
	hJpsiN_pbpb_rap1624->SetMarkerColor(kBlue);
	hJpsiN_pbpb_rap1624->Draw("sameP");

	hFitAndre_pbpb_rap1624->SetLineColor(kRed);
	hFitAndre_pbpb_rap1624->SetMarkerStyle(24);
	hFitAndre_pbpb_rap1624->SetMarkerColor(kRed);
	hFitAndre_pbpb_rap1624->Draw("sameP");


	TLegend*leg4 = new TLegend(0.38, 0.15, 0.85, 0.38, "");
	leg4->SetFillColor(kWhite);
	leg4->SetBorderSize(0);
	leg4->SetTextSize(0.035);
	leg4->AddEntry(hJpsiN_pbpb_rap0016, "PbPb 0-100%, |y|<1.6", "p");
	leg4->AddEntry(hFitAndre_pbpb_rap0016, "PbPb 0-100% fit Andre, |y|<1.6", "p");
	leg4->AddEntry(hJpsiN_pbpb_rap1624, "PbPb 0-100%, 1.6<|y|<2.4", "p");
	leg4->AddEntry(hFitAndre_pbpb_rap1624, "PbPb 0-100% fit Andre, 1.6<|y|<2.4", "p");
	leg4->Draw("same");

	c4->SaveAs("bfrac_plots/bfraction_1S_pbpb_pT.png");
	c4->SaveAs("bfrac_plots/bfraction_1S_pbpb_pT.pdf");




	//2S

	//centrality
	TH1D* h2sJpsiN_cent_rap0016 = new TH1D("h2sJpsiN_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	TH1D* h2sJpsiN_cent_rap1624 = new TH1D("h2sJpsiN_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);

	TH1D* h2sFitAndre_cent_rap0016 = new TH1D("h2sFitAndre_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	TH1D* h2sFitAndre_cent_rap1624 = new TH1D("h2sFitAndre_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);

	//TH1D* h2sFitMihee_cent_rap0016 = new TH1D("h2sFitMihee_cent_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_centmid, bins_centmid);
	//TH1D* h2sFitMihee_cent_rap1624 = new TH1D("h2sFitMihee_cent_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_centfwd, bins_centfwd);


	//pt
	TH1D* h2sJpsiN_pbpb_rap0016 = new TH1D("h2sJpsiN_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* h2sJpsiN_pbpb_rap1624 = new TH1D("h2sJpsiN_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	TH1D* h2sFitAndre_pbpb_rap0016 = new TH1D("h2sFitAndre_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	TH1D* h2sFitAndre_pbpb_rap1624 = new TH1D("h2sFitAndre_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);

	//TH1D* h2sFitMihee_pbpb_rap0016 = new TH1D("h2sFitMihee_pbpb_rap0016", "b contribution for rapidity 0.0 - 1.6", nbins_ptmid, bins_ptmid);
	//TH1D* h2sFitMihee_pbpb_rap1624 = new TH1D("h2sFitMihee_pbpb_rap1624", "b contribution for rapidity 1.6 - 2.4", nbins_ptfwd, bins_ptfwd);


	// filling cent histograms (by hand)

	h2sJpsiN_cent_rap0016->SetBinContent(1, 0.839);//0-10
	h2sJpsiN_cent_rap0016->SetBinContent(2, 0.635);//10-20
	h2sJpsiN_cent_rap0016->SetBinContent(3, 0.572);//20-30
	h2sJpsiN_cent_rap0016->SetBinContent(4, 0.529);//30-40 // !!
	h2sJpsiN_cent_rap0016->SetBinContent(5, 0.586);//40-50 // !!
	h2sJpsiN_cent_rap0016->SetBinContent(6, 0.653);//50-100

	h2sJpsiN_cent_rap0016->SetBinError(1, 0.161);//0-10
	h2sJpsiN_cent_rap0016->SetBinError(2, 0.126);//10-20
	h2sJpsiN_cent_rap0016->SetBinError(3, 0.102);//20-30
	h2sJpsiN_cent_rap0016->SetBinError(4, 0.111);//30-40
	h2sJpsiN_cent_rap0016->SetBinError(5, 0.146);//40-50
	h2sJpsiN_cent_rap0016->SetBinError(6, 0.104);//50-100

	h2sJpsiN_cent_rap1624->SetBinContent(1, 0.589);//0-20
	h2sJpsiN_cent_rap1624->SetBinContent(2, 1.429);//20-40
	h2sJpsiN_cent_rap1624->SetBinContent(3, 0.740);//40-100 // !!

	h2sJpsiN_cent_rap1624->SetBinError(1, 0.213);//0-20
	h2sJpsiN_cent_rap1624->SetBinError(2, 1.431);//20-40
	h2sJpsiN_cent_rap1624->SetBinError(3, 0.379);//40-100

	h2sFitAndre_cent_rap0016->SetBinContent(1, 1.000);//0-10
	h2sFitAndre_cent_rap0016->SetBinContent(2, 0.766);//10-20
	h2sFitAndre_cent_rap0016->SetBinContent(3, 0.499);//20-30
	h2sFitAndre_cent_rap0016->SetBinContent(4, 0.639);//30-40
	h2sFitAndre_cent_rap0016->SetBinContent(5, 0.478);//40-50
	h2sFitAndre_cent_rap0016->SetBinContent(6, 0.608);//50-100

	h2sFitAndre_cent_rap0016->SetBinError(1, 0.133);//0-10
	h2sFitAndre_cent_rap0016->SetBinError(2, 0.107);//10-20
	h2sFitAndre_cent_rap0016->SetBinError(3, 0.077);//20-30
	h2sFitAndre_cent_rap0016->SetBinError(4, 0.094);//30-40
	h2sFitAndre_cent_rap0016->SetBinError(5, 0.104);//40-50
	h2sFitAndre_cent_rap0016->SetBinError(6, 0.093);//50-100

	h2sFitAndre_cent_rap1624->SetBinContent(1, 0.570);//0-20
	h2sFitAndre_cent_rap1624->SetBinContent(2, 1.000);//20-40
	h2sFitAndre_cent_rap1624->SetBinContent(3, 0.689);//40-100
	
	h2sFitAndre_cent_rap1624->SetBinError(1, 0.165);//0-20
	h2sFitAndre_cent_rap1624->SetBinError(2, 0.176);//20-40
	h2sFitAndre_cent_rap1624->SetBinError(3, 0.232);//40-100



	TCanvas*  c5 = new TCanvas("c5", "PbPb cent data Jpsi", 800, 600);

	h2sJpsiN_cent_rap0016->GetXaxis()->SetTitle("Centrality bin");
	h2sJpsiN_cent_rap0016->GetXaxis()->SetRangeUser(0, 100);
	h2sJpsiN_cent_rap0016->GetYaxis()->SetTitle("non-prompt fraction");
	h2sJpsiN_cent_rap0016->GetYaxis()->SetRangeUser(0.0, 1.5);
	h2sJpsiN_cent_rap0016->Draw("P");
	h2sJpsiN_cent_rap0016->SetLineColor(kBlue);
	h2sJpsiN_cent_rap0016->SetTitle("");
	h2sJpsiN_cent_rap0016->SetMarkerStyle(29);
	h2sJpsiN_cent_rap0016->SetMarkerColor(kBlue);

	h2sFitAndre_cent_rap0016->SetLineColor(kRed);
	h2sFitAndre_cent_rap0016->SetMarkerStyle(20);
	h2sFitAndre_cent_rap0016->SetMarkerColor(kRed);
	h2sFitAndre_cent_rap0016->Draw("sameP");

	/*h2sFitMihee_cent_rap0016->SetLineColor(kGreen);
	h2sFitMihee_cent_rap0016->SetMarkerStyle(21);
	h2sFitMihee_cent_rap0016->SetMarkerColor(kGreen);
	h2sFitMihee_cent_rap0016->Draw("sameP");*/


	h2sJpsiN_cent_rap1624->SetLineColor(kBlue);
	h2sJpsiN_cent_rap1624->SetMarkerStyle(30);
	h2sJpsiN_cent_rap1624->SetMarkerColor(kBlue);
	h2sJpsiN_cent_rap1624->Draw("sameP");

	h2sFitAndre_cent_rap1624->SetLineColor(kRed);
	h2sFitAndre_cent_rap1624->SetMarkerStyle(24);
	h2sFitAndre_cent_rap1624->SetMarkerColor(kRed);
	h2sFitAndre_cent_rap1624->Draw("sameP");

	/*h2sFitMihee_cent_rap1624->SetLineColor(kGreen);
	h2sFitMihee_cent_rap1624->SetMarkerStyle(25);
	h2sFitMihee_cent_rap1624->SetMarkerColor(kGreen);
	h2sFitMihee_cent_rap1624->Draw("sameP");*/

	TLegend*leg5 = new TLegend(0.35, 0.65, 0.80, 0.88, "");
	leg5->SetFillColor(kWhite);
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.035);
	leg5->AddEntry(h2sJpsiN_cent_rap0016, "2S PbPb 6.5<p_{T}<30, |y|<1.6", "p");
	leg5->AddEntry(h2sFitAndre_cent_rap0016, "2S PbPb 6.5<p_{T}<30 fit Andre, |y|<1.6", "p");
	//leg5->AddEntry(h2sFitMihee_cent_rap0016, "PbPb 5.02 TeV fit Mihee, |y|<1.6", "p");
	leg5->AddEntry(h2sJpsiN_cent_rap1624, "2S PbPb 3<p_{T}<30, 1.6<|y|<2.4", "p");
	leg5->AddEntry(h2sFitAndre_cent_rap1624, "2S PbPb 3<p_{T}<30 fit Andre, 1.6<|y|<2.4", "p");
	//leg5->AddEntry(h2sFitMihee_cent_rap1624, "PbPb 5.02 TeV fit Mihee, 1.6<|y|<2.4", "p");
	leg5->Draw("same");

	c5->SaveAs("bfrac_plots/bfraction_2S_pbpb_cent.png");
	c5->SaveAs("bfrac_plots/bfraction_2S_pbpb_cent.pdf");



	// filling pbpb histograms pt

	h2sJpsiN_pbpb_rap0016->SetBinContent(1, 0.686);//6.5-9
	h2sJpsiN_pbpb_rap0016->SetBinContent(2, 0.616);//9-12
	h2sJpsiN_pbpb_rap0016->SetBinContent(3, 0.606);//12-15
	h2sJpsiN_pbpb_rap0016->SetBinContent(4, 0.738);//15-20
	h2sJpsiN_pbpb_rap0016->SetBinContent(5, 0.828);//20-30

	h2sJpsiN_pbpb_rap0016->SetBinError(1, 0.139);//6.5-9
	h2sJpsiN_pbpb_rap0016->SetBinError(2, 0.077);//9-12
	h2sJpsiN_pbpb_rap0016->SetBinError(3, 0.080);//12-15
	h2sJpsiN_pbpb_rap0016->SetBinError(4, 0.090);//15-20
	h2sJpsiN_pbpb_rap0016->SetBinError(5, 0.130);//20-30

	h2sJpsiN_pbpb_rap1624->SetBinContent(1, 0.886);//3-6.5
	h2sJpsiN_pbpb_rap1624->SetBinContent(2, 0.592);//6.5-12
	h2sJpsiN_pbpb_rap1624->SetBinContent(3, 0.461);//12-30

	h2sJpsiN_pbpb_rap1624->SetBinError(1, 0.565);//3-6.5
	h2sJpsiN_pbpb_rap1624->SetBinError(2, 0.163);//6.5-12
	h2sJpsiN_pbpb_rap1624->SetBinError(3, 0.149);//12-30

	h2sFitAndre_pbpb_rap0016->SetBinContent(1, 0.765);//6.5-9
	h2sFitAndre_pbpb_rap0016->SetBinContent(2, 0.715);//9-12
	h2sFitAndre_pbpb_rap0016->SetBinContent(3, 0.644);//12-15
	h2sFitAndre_pbpb_rap0016->SetBinContent(4, 0.745);//15-20
	h2sFitAndre_pbpb_rap0016->SetBinContent(5, 0.797);//20-30

	h2sFitAndre_pbpb_rap0016->SetBinError(1, 0.133);//6.5-9
	h2sFitAndre_pbpb_rap0016->SetBinError(2, 0.063);//9-12
	h2sFitAndre_pbpb_rap0016->SetBinError(3, 0.071);//12-15
	h2sFitAndre_pbpb_rap0016->SetBinError(4, 0.080);//15-20
	h2sFitAndre_pbpb_rap0016->SetBinError(5, 0.109);//20-30

	h2sFitAndre_pbpb_rap1624->SetBinContent(1, 0.846);//3-6.5
	h2sFitAndre_pbpb_rap1624->SetBinContent(2, 0.746);//6.5-12
	h2sFitAndre_pbpb_rap1624->SetBinContent(3, 0.712);//12-30

	h2sFitAndre_pbpb_rap1624->SetBinError(1, 0.290);//3-6.5
	h2sFitAndre_pbpb_rap1624->SetBinError(2, 0.149);//6.5-12
	h2sFitAndre_pbpb_rap1624->SetBinError(3, 0.138);//12-30




	TCanvas*  c6 = new TCanvas("c6", "pbpb data Jpsi", 800, 600);

	TH1I* hframe6 = new TH1I("hframe6", "hframe6", 1, 0, 30); //just holder to set the axis properly
	hframe6->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hframe6->GetXaxis()->SetRangeUser(0.0, 30.0);
	hframe6->GetYaxis()->SetTitle("non-prompt fraction");
	hframe6->GetYaxis()->SetRangeUser(0.0, 1.0);
	hframe6->Draw();
	hframe6->SetTitle("");

	//h2sJpsiN_pbpb_rap0016->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	//h2sJpsiN_pbpb_rap0016->GetXaxis()->SetRangeUser(0.0, 30.0);
	//h2sJpsiN_pbpb_rap0016->GetYaxis()->SetTitle("non-prompt fraction");
	//h2sJpsiN_pbpb_rap0016->GetYaxis()->SetRangeUser(0.0, 1.0);
	h2sJpsiN_pbpb_rap0016->Draw("sameP");
	h2sJpsiN_pbpb_rap0016->SetLineColor(kBlue);
	//h2sJpsiN_pbpb_rap0016->SetTitle("");
	h2sJpsiN_pbpb_rap0016->SetMarkerStyle(29);
	h2sJpsiN_pbpb_rap0016->SetMarkerColor(kBlue);

	h2sFitAndre_pbpb_rap0016->SetLineColor(kRed);
	h2sFitAndre_pbpb_rap0016->SetMarkerStyle(20);
	h2sFitAndre_pbpb_rap0016->SetMarkerColor(kRed);
	h2sFitAndre_pbpb_rap0016->Draw("sameP");

	h2sJpsiN_pbpb_rap1624->SetLineColor(kBlue);
	h2sJpsiN_pbpb_rap1624->SetMarkerStyle(30);
	h2sJpsiN_pbpb_rap1624->SetMarkerColor(kBlue);
	h2sJpsiN_pbpb_rap1624->Draw("sameP");

	h2sFitAndre_pbpb_rap1624->SetLineColor(kRed);
	h2sFitAndre_pbpb_rap1624->SetMarkerStyle(24);
	h2sFitAndre_pbpb_rap1624->SetMarkerColor(kRed);
	h2sFitAndre_pbpb_rap1624->Draw("sameP");


	TLegend*leg6 = new TLegend(0.35, 0.15, 0.85, 0.38, "");
	leg6->SetFillColor(kWhite);
	leg6->SetBorderSize(0);
	leg6->SetTextSize(0.035);
	leg6->AddEntry(h2sJpsiN_pbpb_rap0016, "2S PbPb 0-100%, |y|<1.6", "p");
	leg6->AddEntry(h2sFitAndre_pbpb_rap0016, "2S PbPb 0-100% fit Andre, |y|<1.6", "p");
	leg6->AddEntry(h2sJpsiN_pbpb_rap1624, "2S PbPb 0-100%, 1.6<|y|<2.4", "p");
	leg6->AddEntry(h2sFitAndre_pbpb_rap1624, "2S PbPb 0-100% fit Andre, 1.6<|y|<2.4", "p");
	leg6->Draw("same");

	c6->SaveAs("bfrac_plots/bfraction_2S_pbpb_pT.png");
	c6->SaveAs("bfrac_plots/bfraction_2S_pbpb_pT.pdf");

	



}
