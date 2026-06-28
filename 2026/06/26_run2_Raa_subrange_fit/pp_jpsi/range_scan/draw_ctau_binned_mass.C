/*
 * Draw dimuon mass spectra in ctau3D slices for a selected pp J/psi
 * kinematic bin.
 *
 * Workflow:
 *   - Load the pp J/psi RooDataSet from kDataRoot.
 *   - Apply the mass, pT, |y|, opposite-sign, muon acceptance, and finite
 *     ctau3D/ctau3DRes selections.
 *   - Split the selected sample into fixed ctau3D bins from -20 to 20 mm.
 *   - For each ctau3D slice, draw the m_mumu spectrum, mark the J/psi mass,
 *     and count entries in the signal-window and sideband regions.
 *
 * Outputs:
 *   - Per-slice mass PDFs under figs/<y-bin>/<pt-bin>/.
 *   - Per-slice histograms in roots/<y-bin>/<pt-bin>/ctau_binned_mass_*.root.
 *   - A TSV summary with entries, mass binning, signal-window counts,
 *     sideband counts, and J/side ratios.
 *
 * Example:
 *   root -l -b -q "draw_ctau_binned_mass.C(3.5,6.5,1.6,2.4,false)"
 */
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

namespace {

// ------------------------------------------------------------------
// input configuration and helper utilities
// ------------------------------------------------------------------
const char *kDataRoot = "/data/users/pjgwak/work/raa_pb18/run2_raa_pbpb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_260303_root634.root";
const char *kDatasetName = "dataset";
const char *kAccCut = "(((abs(eta1) <= 1.2) && (pt1 >= 3.5)) || ((abs(eta2) <= 1.2) && (pt2 >= 3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)))";

TString tag_value(double x)
{
	TString s = TString::Format("%.2f", x);
	s.ReplaceAll("-", "m");
	s.ReplaceAll(".", "p");
	return s;
}

TString y_tag(double yLow, double yHigh)
{
	return TString::Format("y%s_%s", tag_value(yLow).Data(), tag_value(yHigh).Data());
}

TString pt_tag(double ptLow, double ptHigh)
{
	return TString::Format("pt%s_%s", tag_value(ptLow).Data(), tag_value(ptHigh).Data());
}

TString bin_tag(double ptLow, double ptHigh, double yLow, double yHigh)
{
	return TString::Format("%s_%s", y_tag(yLow, yHigh).Data(), pt_tag(ptLow, ptHigh).Data());
}

int adaptive_mass_bins(int entries)
{
	if (entries <= 0)
		return 0;
	const int bins = static_cast<int>(std::lround(static_cast<double>(entries) / 45.0));
	return std::clamp(bins, 25, 100);
}

struct SliceCount {
	int entries = 0;
	int jpsi = 0;
	int leftSide = 0;
	int rightSide = 0;
};

SliceCount fill_mass_hist(RooDataSet &ds, TH1D &h)
{
	SliceCount count;
	for (int i = 0; i < ds.numEntries(); ++i) {
		const RooArgSet *row = ds.get(i);
		const auto *mass = dynamic_cast<const RooRealVar *>(row->find("mass"));
		if (!mass)
			continue;
		const double m = mass->getVal();
		h.Fill(m);
		++count.entries;
		if (m >= 2.9 && m <= 3.2)
			++count.jpsi;
		if (m >= 2.60 && m < 2.90)
			++count.leftSide;
		if (m > 3.30 && m <= 3.50)
			++count.rightSide;
	}
	return count;
}
}

// ------------------------------------------------------------------
// main macro
// ------------------------------------------------------------------
void draw_ctau_binned_mass(double ptLow = 3.5, double ptHigh = 6.5,
	double yLow = 1.6, double yHigh = 2.4,
	bool useLogY = false)
{
	// ------------------------------------------------------------------
	// ROOT style setup
	// ------------------------------------------------------------------
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
	gStyle->SetOptStat(0);

	// ------------------------------------------------------------------
	// load input dataset
	// ------------------------------------------------------------------
	TFile input(kDataRoot, "READ");
	if (input.IsZombie()) {
		std::cerr << "ERROR: cannot open input data file: " << kDataRoot << std::endl;
		return;
	}
	auto *data = dynamic_cast<RooDataSet *>(input.Get(kDatasetName));
	if (!data) {
		std::cerr << "ERROR: dataset not found: " << kDatasetName << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// output directory and file-name tags
	// ------------------------------------------------------------------
	const TString yTag = y_tag(yLow, yHigh);
	const TString pTag = pt_tag(ptLow, ptHigh);
	const TString tag = bin_tag(ptLow, ptHigh, yLow, yHigh);
	const TString outFigDir = TString::Format("figs/%s/%s", yTag.Data(), pTag.Data());
	const TString outRootDir = TString::Format("roots/%s/%s", yTag.Data(), pTag.Data());
	gSystem->mkdir(outFigDir, true);
	gSystem->mkdir(outRootDir, true);

	// ------------------------------------------------------------------
	// event selection
	// ------------------------------------------------------------------
	const TString baseCut = Form(
		"(mass>2.6 && mass<3.5) && (pt > %.12g && pt < %.12g) && "
		"(fabs(y) > %.12g && fabs(y) < %.12g) && (recoQQsign==0) && %s && "
		"!TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes) && "
		"(ctau3D >= -20.0 && ctau3D < 20.0)",
		ptLow, ptHigh, yLow, yHigh, kAccCut);
	std::unique_ptr<RooDataSet> selected(static_cast<RooDataSet *>(data->reduce(baseCut)));
	if (!selected || selected->numEntries() <= 0) {
		std::cerr << "ERROR: no entries after selection: " << baseCut << std::endl;
		return;
	}

	// ------------------------------------------------------------------
	// ctau3D binning and output handles
	// ------------------------------------------------------------------
	const std::vector<double> edges = {
		-20.0, -10.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, -0.10,
		0.0, 0.10, 0.25, 0.50, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0
	};
	const TString rootName = TString::Format("%s/ctau_binned_mass_%s.root", outRootDir.Data(), tag.Data());
	const TString countsName = TString::Format("%s/ctau_binned_counts_%s.tsv", outFigDir.Data(), tag.Data());

	TFile outFile(rootName, "RECREATE");
	if (outFile.IsZombie()) {
		std::cerr << "ERROR: cannot create output ROOT file: " << rootName << std::endl;
		return;
	}
	FILE *counts = fopen(countsName.Data(), "w");
	if (!counts) {
		std::cerr << "ERROR: cannot create counts file: " << countsName << std::endl;
		return;
	}
	fprintf(counts, "idx\tctLow\tctHigh\tentries\tmassBins\tjpsiWindow_2p976_3p216\tleftSide_2p60_2p90\trightSide_3p30_3p50\tjpsiOverSide\n");

	std::cout << Form("[INFO] selected entries: %d", selected->numEntries()) << std::endl;
	std::cout << Form("[INFO] ctau binning: %zu bins from %.1f to %.1f mm", edges.size() - 1, edges.front(), edges.back()) << std::endl;

	// ------------------------------------------------------------------
	// loop over ctau3D slices
	// ------------------------------------------------------------------
	TCanvas c("c_ctau_binned_mass", "c_ctau_binned_mass", 800, 800);

	for (size_t i = 0; i + 1 < edges.size(); ++i) {
		const double ctLow = edges[i];
		const double ctHigh = edges[i + 1];
		const TString cut = Form("ctau3D >= %.12g && ctau3D < %.12g", ctLow, ctHigh);
		std::unique_ptr<RooDataSet> slice(static_cast<RooDataSet *>(selected->reduce(cut)));
		const int entries = slice ? slice->numEntries() : 0;
		const int massBins = adaptive_mass_bins(entries);
		if (!slice || entries <= 0) {
			fprintf(counts, "%zu\t%.8g\t%.8g\t0\t0\t0\t0\t0\t0\n", i, ctLow, ctHigh);
			continue;
		}

		// Build the mass histogram and count signal-window/sideband entries.
		TH1D h(Form("hMass_ct_%02zu", i),
			Form(";m_{#mu#mu} [GeV/c^{2}];Events / %.1f MeV", 900.0 / massBins),
			massBins, 2.6, 3.5);
		h.Sumw2();
		const SliceCount count = fill_mass_hist(*slice, h);
		const int side = count.leftSide + count.rightSide;
		const double jpsiOverSide = side > 0 ? static_cast<double>(count.jpsi) / static_cast<double>(side) : 0.0;

		fprintf(counts, "%zu\t%.8g\t%.8g\t%d\t%d\t%d\t%d\t%d\t%.8g\n",
			i, ctLow, ctHigh, count.entries, massBins, count.jpsi, count.leftSide, count.rightSide, jpsiOverSide);
		std::cout << Form("[SLICE %02zu] ctau [%.3g, %.3g): entries=%d, J/psi=%d, side=%d, J/side=%.3f",
			i, ctLow, ctHigh, count.entries, count.jpsi, side, jpsiOverSide) << std::endl;

		outFile.cd();
		h.Write();

		// Draw the per-slice mass spectrum and mark the nominal J/psi mass.
		c.Clear();
		c.SetFillStyle(4000);
		c.SetFrameFillStyle(4000);
		c.SetLeftMargin(0.16);
		c.SetRightMargin(0.04);
		c.SetTopMargin(0.08);
		c.SetBottomMargin(0.13);
		c.SetTicks(1, 1);
		c.SetLogy(useLogY);
		h.SetTitle("");
		h.SetMarkerStyle(20);
		h.SetMarkerSize(0.80);
		h.SetLineColor(kBlack);
		h.SetLineWidth(2);
		h.SetMarkerColor(kBlack);
		h.GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
		h.GetYaxis()->SetTitle("Events");
		h.GetYaxis()->SetTitleOffset(1.45);
		h.GetXaxis()->SetTitleSize(0.050);
		h.GetYaxis()->SetTitleSize(0.050);
		h.GetXaxis()->SetLabelSize(0.042);
		h.GetYaxis()->SetLabelSize(0.036);
		h.GetXaxis()->CenterTitle();
		h.GetYaxis()->CenterTitle();
		const double ymax = std::max(1.0, h.GetMaximum());
		h.SetMinimum(useLogY ? 0.5 : 0.0);
		const double yAxisMax = useLogY ? ymax * 16.0 : ymax * 1.65;
		const double jpsiLineMax = useLogY ? ymax * 7.0 : ymax * 1.05;
		h.SetMaximum(yAxisMax);
		h.Draw("E");

		TLine jpsiLine(3.096, useLogY ? 0.5 : 0.0, 3.096, jpsiLineMax);
		jpsiLine.SetLineStyle(2);
		jpsiLine.SetLineColor(kBlue + 1);
		jpsiLine.SetLineWidth(2);
		jpsiLine.Draw("same");

		TLegend leg(0.56, 0.74, 0.80, 0.89);
		leg.SetBorderSize(0);
		leg.SetFillStyle(0);
		leg.SetTextSize(0.024);
		leg.SetEntrySeparation(0.015);
		leg.AddEntry(&h, "Data", "lep");
		leg.AddEntry(&jpsiLine, "J/#psi mean", "l");
		leg.Draw("same");

		// Add CMS-style labels and per-slice bookkeeping text.
		{
			TLatex tx;
			tx.SetNDC();
			tx.SetTextSize(0.032);
			tx.SetTextFont(42);
			tx.SetTextAlign(31);
			tx.DrawLatex(0.96, 0.935, "pp #sqrt{s} = 5.02 TeV (28.0 pb^{-1})");
		}
		{
			TLatex tx;
			tx.SetNDC();
			tx.SetTextAlign(11);
			tx.SetTextSize(0.042);
			tx.SetTextFont(62);
			tx.DrawLatex(0.23, 0.930, "CMS");
			tx.SetTextFont(52);
			tx.SetTextSize(0.032);
			tx.DrawLatex(0.325, 0.930, "Internal");
		}
		{
			TLatex tx;
			tx.SetNDC();
			tx.SetTextSize(0.026);
			tx.SetTextFont(42);
			double xtext = 0.19;
			double y0 = 0.865;
			double dy = -0.048;
			int k = 0;
			tx.DrawLatex(xtext, y0 + dy * k++, "J/#psi data");
			if (yLow == 0.0)
				tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
			else
				tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
			tx.DrawLatex(xtext, y0 + dy * k++, Form("%.3f < ctau3D < %.3f mm", ctLow, ctHigh));
		}
		{
			TLatex tp;
			tp.SetNDC();
			tp.SetTextSize(0.024);
			tp.SetTextFont(42);
			double xtext = 0.74;
			double y0 = 0.87;
			double dy = -0.045;
			int k = 0;
			tp.DrawLatex(xtext, y0 + dy * k++, Form("Entries = %d", count.entries));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("J/#psi window = %d", count.jpsi));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("Sideband = %d", side));
			tp.DrawLatex(xtext, y0 + dy * k++, Form("J/side = %.3f", jpsiOverSide));
		}

		// Save the current ctau3D slice as an individual PDF.
		TString slicePdf = TString::Format("%s/mass_%s_ct%s_%s.pdf",
			outFigDir.Data(), tag.Data(), tag_value(ctLow).Data(), tag_value(ctHigh).Data());
		c.Print(slicePdf);
	}

	// ------------------------------------------------------------------
	// finalize outputs
	// ------------------------------------------------------------------
	fclose(counts);
	outFile.Write();
	outFile.Close();
	std::cout << "[DONE] wrote slice PDFs in " << outFigDir << std::endl;
	std::cout << "[DONE] wrote " << rootName << std::endl;
	std::cout << "[DONE] wrote " << countsName << std::endl;
}
