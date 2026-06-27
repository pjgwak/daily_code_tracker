#include <iostream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#include "../headers/rootFitHeaders.h"
#include "../headers/commonUtility.h"
#include "../headers/JpsiUtility.h"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooCrystalBall.h>
#include "TStopwatch.h"
#include <RooWorkspace.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TFile.h"
#include "TParameter.h"
#include "../headers/cutsAndBin.h"
#include "../headers/CMS_lumi_v2mass.C"
#include "../headers/tdrstyle.C"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "saved_fit_helpers.h"

using namespace std;
using namespace RooFit;

struct ScopedMacroTimer
{
    TStopwatch sw;
    TString label;

    ScopedMacroTimer(const char *macroName, float ptLow, float ptHigh, float yLow, float yHigh)
        : label(TString::Format("%s(pt=%.2f-%.2f, |y|=%.2f-%.2f)", macroName, ptLow, ptHigh, yLow, yHigh))
    {
        sw.Start();
        std::cout << "[TIMER] start " << label << std::endl;
    }

    ~ScopedMacroTimer()
    {
        sw.Stop();
        std::cout << Form("[TIMER] done %s | real %.2f s | cpu %.2f s", label.Data(), sw.RealTime(), sw.CpuTime()) << std::endl;
    }
};

static std::pair<double, int> chi2_from_pull(RooHist *hpull)
{
    double chi2 = 0.0;
    int n = 0;
    if (!hpull) return {0.0, 0};
    double x = 0.0, y = 0.0;
    for (int i = 0; i < hpull->GetN(); ++i) {
        hpull->GetPoint(i, x, y);
        if (!std::isfinite(y)) continue;
        chi2 += y * y;
        ++n;
    }
    return {chi2, n};
}

void mc_mass(
		float ptLow=6.5, float ptHigh=9.0,
		float yLow=0, float yHigh=1.6,
		int PR=0, //0=PR, 1=NP, 2=Inc.
		int PRw=1, bool fEffW = false, bool fAccW = false, bool isPtW = false, bool isTnP = false,
		int nSignalCBComponents = 2,    // 0~2
		int nSignalGaussComponents = 1, // 0~2
		bool drawFromSavedFit = false,
		bool publish = false
		)
{
    ScopedMacroTimer timer("mc_mass", ptLow, ptHigh, yLow, yHigh);
    if (publish) drawFromSavedFit = true;
    gStyle->SetEndErrorSize(0);

	// ------------------------------------------------------------------
	// output and model control
	// ------------------------------------------------------------------

	TString bCont;
	if(PR==0) bCont="Prompt";
	else if(PR==1) bCont="NonPrompt";
	else if(PR==2) bCont="Inclusive";

	TString fname;
	if (PRw==1) fname="PR";
	else if (PRw==2) fname="NP";
	else fname="Inc";

	nSignalCBComponents = std::clamp(nSignalCBComponents, 0, 2);
	nSignalGaussComponents = std::clamp(nSignalGaussComponents, 0, 2);
	if (nSignalCBComponents + nSignalGaussComponents <= 0) {
		cerr << "ERROR: at least one signal component is required." << endl;
		return;
	}
	const TString modelSelTag = TString::Format("CB%d_G%d", nSignalCBComponents, nSignalGaussComponents);
	auto formatTag = [](double value) { return TString::Format("%.2f", value); };
	const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
	const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
	const TString figBaseDir = publish ? "figs_publish" : "figs";
	const TString figDir = TString::Format("%s/%s/mc_mass", figBaseDir.Data(), yTag.Data());
	const TString resultDir = TString::Format("roots/%s/mc_mass", yTag.Data());
	const TString figTag = yTag + "_" + ptTag + "_bkgOff";
	const TString modelTag = yTag + "_" + ptTag;
	auto figName = [&](const char *name) {
		return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
	};
	const TString modelFileName = TString::Format("%s/mc_mass_model_%s.root", resultDir.Data(), modelTag.Data());
	gSystem->mkdir(figDir, kTRUE);
	gSystem->mkdir(resultDir, kTRUE);
	gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

	std::unique_ptr<TFile> savedFitFile;
	if (drawFromSavedFit && !load_saved_fit_file(savedFitFile, modelFileName, "MC mass")) return;
	if (drawFromSavedFit) {
		nSignalGaussComponents = std::clamp(read_saved_int_param(savedFitFile.get(), "nSignalGaussComponents", nSignalGaussComponents), 0, 2);
		nSignalCBComponents = std::clamp(read_saved_int_param(savedFitFile.get(), "nSignalCBComponents", nSignalCBComponents), 0, 2);
	}
	
	// RooMsgService::instance().getStream(0).removeTopic(Caching);
	// RooMsgService::instance().getStream(1).removeTopic(Caching);
	// RooMsgService::instance().getStream(0).removeTopic(Plotting);
	// RooMsgService::instance().getStream(1).removeTopic(Plotting);
	// RooMsgService::instance().getStream(0).removeTopic(Integration);
	// RooMsgService::instance().getStream(1).removeTopic(Integration);
	// RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
	// RooMsgService::instance().getStream(1).removeTopic(Fitting);
	// RooMsgService::instance().getStream(1).removeTopic(Minimization);
	// RooMsgService::instance().getStream(1).removeTopic(InputArguments);
	// RooMsgService::instance().getStream(1).removeTopic(Eval);
	// RooMsgService::instance().getStream(1).removeTopic(DataHandling);
	// // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
	// RooMsgService::instance().setGlobalKillBelow(ERROR);
	// RooMsgService::instance().setSilentMode(true);

	// ------------------------------------------------------------------
	// load dataset
	// ------------------------------------------------------------------
    TFile* f1 = new TFile("/data/users/pjgwak/work/raa_pb18/run2_raa_pbpb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root", "read");

	massLow=2.6;
	massHigh=3.5;

	// ------------------------------------------------------------------
	// event selection
	// ------------------------------------------------------------------
	TString kineCut;
	kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && mass>%.2f && mass<%.2f",ptLow, ptHigh, yLow, yHigh, massLow, massHigh);

	TString accCut = "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) ) &&";//2018 acceptance cut

	TString OS="recoQQsign==0 &&";

	kineCut = OS+accCut+kineCut;


	RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
	if (!dataset) {
		cerr << "ERROR: dataset not found in input file: " << f1->GetName() << endl;
		f1->ls();
		return;
	}
	RooWorkspace *ws = new RooWorkspace("workspace");
	ws->import(*dataset);
	ws->data("dataset")->Print();
	cout << "pt: "<<ptLow<<"-"<<ptHigh<<", y: "<<yLow<<"-"<<yHigh<<endl;
	cout << "####################################" << endl;
	RooDataSet *datasetW = new RooDataSet("datasetW","A sample",*dataset->get(),Import(*dataset),WeightVar(*ws->var("weight")));
	RooDataSet *dsAB = (RooDataSet*)datasetW->reduce(RooArgSet(*(ws->var("ctau3DRes")),*(ws->var("ctau3D")), *(ws->var("ctau3DErr")), *(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );

    cout << "======= dataset ======\n";
    dataset->Print("V");
    cout << endl << endl;

    cout << "======= datasetW ======\n";
    datasetW->Print("V");
    cout << endl << endl;

    cout << "======= dsAB ======\n";
    dsAB->Print("V");
    cout << endl << endl;

	cout << "******** New Combined Dataset ***********" << endl;
	dsAB->SetName("dsAB");
	ws->import(*dsAB);
	ws->var("mass")->setRange(massLow, massHigh);
	ws->var("mass")->Print();
	// ------------------------------------------------------------------
	// mass fit model setup
	// ------------------------------------------------------------------

    RooRealVar signal_mass_mean("signal_mass_mean", "signal_mass_mean", pdgMass.JPsi, pdgMass.JPsi-0.1, pdgMass.JPsi+0.1, "GeV/c^{2}");
    RooRealVar signal_mass_sigma("signal_mass_sigma", "signal_mass_sigma", 0.03, 0.001, 0.080, "GeV/c^{2}");
    const double sigmaFloor = 1e-4;
    const double signal_mass_sigma2_init = std::max(signal_mass_sigma.getVal() + 0.03, sigmaFloor + 1e-4);
    RooRealVar signal_mass_sigma2_log_offset(
            "signal_mass_sigma2_log_offset", "log(signal_mass_sigma2 - floor)",
            std::log(std::max(signal_mass_sigma2_init - sigmaFloor, 1e-4)),
            std::log(1e-4),
            std::log(1.0));
    RooFormulaVar signal_mass_sigma2("signal_mass_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_sigma2_log_offset));
    RooFormulaVar signal_mass_sigma_delta2("signal_mass_sigma_delta2", "@0-@1", RooArgList(signal_mass_sigma2, signal_mass_sigma));

    RooRealVar signal_mass_cb_sigma_base("signal_mass_cb_sigma_base", "signal_mass_cb_sigma_base", 0.035, 0.008, 0.080, "GeV/c^{2}");
    std::unique_ptr<RooRealVar> signal_mass_cb_sigma_log_offset;
    std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_delta_from_gaus;
    std::unique_ptr<RooFormulaVar> signal_mass_cb_sigma_from_gaus;
    RooAbsReal *signal_mass_cb_sigma = &signal_mass_cb_sigma_base;
    RooAbsReal *signal_mass_cb_sigma_delta = nullptr;
    if (nSignalGaussComponents >= 1)
    {
        const double signal_mass_cb_sigma_init = std::max(signal_mass_sigma.getVal() + 0.005, sigmaFloor + 1e-4);
        signal_mass_cb_sigma_log_offset = std::make_unique<RooRealVar>(
                "signal_mass_cb_sigma_log_offset", "log(signal_mass_cb_sigma - floor)",
                std::log(std::max(signal_mass_cb_sigma_init - sigmaFloor, 1e-4)),
                std::log(1e-4),
                std::log(1.0));
        signal_mass_cb_sigma_from_gaus = std::make_unique<RooFormulaVar>(
                "signal_mass_cb_sigma", Form("%g + exp(@0)", sigmaFloor), RooArgList(*signal_mass_cb_sigma_log_offset));
        signal_mass_cb_sigma_delta_from_gaus = std::make_unique<RooFormulaVar>(
                "signal_mass_cb_sigma_delta", "@0-@1", RooArgList(*signal_mass_cb_sigma_from_gaus, signal_mass_sigma));
        signal_mass_cb_sigma = signal_mass_cb_sigma_from_gaus.get();
        signal_mass_cb_sigma_delta = signal_mass_cb_sigma_delta_from_gaus.get();
    }
    else
    {
        signal_mass_cb_sigma_delta = &signal_mass_cb_sigma_base;
    }

    const double signal_mass_cb_sigma2_init = std::max(signal_mass_cb_sigma->getVal() + 0.02, sigmaFloor + 1e-4);
    RooRealVar signal_mass_cb_sigma2_log_offset(
            "signal_mass_cb_sigma2_log_offset", "log(signal_mass_cb_sigma2 - floor)",
            std::log(std::max(signal_mass_cb_sigma2_init - sigmaFloor, 1e-4)),
            std::log(1e-4),
            std::log(1.0));
    RooFormulaVar signal_mass_cb_sigma2("signal_mass_cb_sigma2", Form("%g + exp(@0)", sigmaFloor), RooArgList(signal_mass_cb_sigma2_log_offset));
    RooFormulaVar signal_mass_cb_sigma_delta2("signal_mass_cb_sigma_delta2", "@0-@1", RooArgList(signal_mass_cb_sigma2, *signal_mass_cb_sigma));
    RooRealVar signal_mass_cb_alpha("signal_mass_cb_alpha", "signal_mass_cb_alpha", 1.5, 0.001, 30.0);
    RooRealVar signal_mass_cb_alpha2("signal_mass_cb_alpha2", "signal_mass_cb_alpha2", 2.0, 0.001, 30.0);
    RooRealVar signal_mass_cb_n("signal_mass_cb_n", "signal_mass_cb_n", 3.0, 0.001, 20.0);
    RooRealVar signal_mass_cb_n2("signal_mass_cb_n2", "signal_mass_cb_n2", 4.0, 0.001, 50.0);

    RooRealVar signal_mass_frac_ratio1("signal_mass_frac_ratio1", "signal_mass_frac_ratio1", 1.86, 1e-3, 1e3);
    RooFormulaVar signal_mass_frac1("signal_mass_frac1", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio1));
    RooRealVar signal_mass_frac_ratio2("signal_mass_frac_ratio2", "signal_mass_frac_ratio2", 0.25, 1e-3, 1e3);
    RooFormulaVar signal_mass_frac2("signal_mass_frac2", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio2));
    RooRealVar signal_mass_frac_ratio3("signal_mass_frac_ratio3", "signal_mass_frac_ratio3", 0.5, 1e-3, 1e3);
    RooFormulaVar signal_mass_frac3("signal_mass_frac3", "@0/(1.0+@0)", RooArgList(signal_mass_frac_ratio3));

    RooGaussian signal_mass_gaus("signal_mass_gaus", "signal_mass_gaus", *(ws->var("mass")), signal_mass_mean, signal_mass_sigma);
    RooGaussian signal_mass_gaus2("signal_mass_gaus2", "signal_mass_gaus2", *(ws->var("mass")), signal_mass_mean, signal_mass_sigma2);
    std::unique_ptr<RooAbsPdf> signal_mass_cb_owned;
    if (nSignalCBComponents >= 2)
    {
        signal_mass_cb_owned = std::make_unique<RooCrystalBall>(
                "signal_mass_cb", "signal_mass_cb",
                *(ws->var("mass")), signal_mass_mean,
                *signal_mass_cb_sigma, signal_mass_cb_sigma2,
                signal_mass_cb_alpha, signal_mass_cb_n,
                signal_mass_cb_alpha2, signal_mass_cb_n2);
    }
    else if (nSignalCBComponents >= 1)
    {
        signal_mass_cb_owned = std::make_unique<RooCBShape>(
                "signal_mass_cb", "signal_mass_cb",
                *(ws->var("mass")), signal_mass_mean, *signal_mass_cb_sigma, signal_mass_cb_alpha, signal_mass_cb_n);
    }
    RooAbsPdf *signal_mass_cb = signal_mass_cb_owned.get();

    std::vector<RooAbsPdf *> signalComponents;
    if (nSignalGaussComponents >= 1)
        signalComponents.push_back(&signal_mass_gaus);
    if (nSignalCBComponents >= 1)
        signalComponents.push_back(signal_mass_cb);
    if (nSignalGaussComponents >= 2)
        signalComponents.push_back(&signal_mass_gaus2);

    RooArgList signalPdfList;
    for (RooAbsPdf *pdf : signalComponents)
        signalPdfList.add(*pdf);
    RooArgList signalFracList;
    if (signalComponents.size() >= 2)
        signalFracList.add(signal_mass_frac1);
    if (signalComponents.size() >= 3)
        signalFracList.add(signal_mass_frac2);
    if (signalComponents.size() >= 4)
        signalFracList.add(signal_mass_frac3);

    std::unique_ptr<RooAbsPdf> signal_mass_owned;
    RooAbsPdf *signal_mass_pdf = nullptr;
    if (signalComponents.size() == 1)
    {
        signal_mass_pdf = signalComponents.front();
    }
    else
    {
        signal_mass_owned = std::make_unique<RooAddPdf>("signal_mass", "signal_mass", signalPdfList, signalFracList, true);
        signal_mass_pdf = signal_mass_owned.get();
    }

    RooRealVar Nsig("Nsig", "Nsig", std::max(1.0, 0.5 * dsAB->sumEntries()), 0.0, std::max(1.0, 2.0 * dsAB->sumEntries()));
    std::unique_ptr<RooAbsPdf> mass_pdf = std::make_unique<RooAddPdf>("mass_pdf", "mass_pdf", RooArgList(*signal_mass_pdf), RooArgList(Nsig));
    ws->import(*mass_pdf);

	// ------------------------------------------------------------------
	// draw mass fit
	// ------------------------------------------------------------------
    TCanvas* c_A =  new TCanvas("canvas_A","canvas_A",800,800);
    c_A->cd();
    TPad *pad_A_1 = new TPad("pad_A_1", "pad_A_1", 0.0, 0.25, 1.0, 1.0);
    pad_A_1->SetBottomMargin(0.00001);
    pad_A_1->SetTopMargin(0.08);
    pad_A_1->SetTicks(1,1);
    pad_A_1->Draw(); pad_A_1->cd();
    const int massPlotBins = 200;
    RooPlot* myPlot_A = ws->var("mass")->frame(Bins(massPlotBins), Title(""));
    myPlot_A->SetTitle("");
    ws->data("dsAB")->plotOn(myPlot_A, Binning(massPlotBins), Name("dataHist_A"));

    pad_A_1->cd();
    bool logY_flag = false;
    if (logY_flag == true) {
        gPad->SetLogy();
    }
    RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
    dsAB->plotOn(myPlot2_A, Binning(massPlotBins), Name("dataOS"), MarkerSize(.8));
    bool isWeighted = ws->data("dsAB")->isWeighted();
    std::unique_ptr<RooFitResult> mass_result;
    if (drawFromSavedFit) {
        mass_result = clone_saved_fit_result(savedFitFile.get(), "fit_result");
        if (!mass_result) {
            std::cerr << "ERROR: fit_result not found in saved MC mass file: " << modelFileName << std::endl;
            return;
        }
        apply_saved_fit_result(mass_result.get(), *mass_pdf, RooArgSet(*(ws->var("mass"))));
        std::cout << "[PlotOnly] Loaded saved MC mass fit and left ROOT file unchanged: " << modelFileName << std::endl;
    }
    else {
        cout << endl << "********* Starting Mass Dist. Fit **************" << endl << endl;
        mass_result.reset(mass_pdf->fitTo(*dsAB, Extended(), Save(), Strategy(2), PrintLevel(0), SumW2Error(isWeighted), Minimizer("Minuit2", "migrad")));
        if (mass_result && (mass_result->status() != 0 || mass_result->covQual() < 2)) {
            std::cout << "[Retry] MC mass fit status=" << mass_result->status() << " covQual=" << mass_result->covQual() << "; retry default minimizer" << std::endl;
            mass_result.reset(mass_pdf->fitTo(*dsAB, Extended(), Save(), Strategy(2), PrintLevel(0), SumW2Error(isWeighted)));
        }
        cout << endl << "********* Finished Mass Dist. Fit **************" << endl << endl;
    }
    RooFitResult* fitMass = mass_result.get();
    if (!fitMass) {
        std::cerr << "ERROR: mass fit result is null." << std::endl;
        return;
    }
    
    // Check and get fitted parameters. Use named variables so model selection does not change indexing assumptions.
    const RooArgList & fitParams = fitMass->floatParsFinal();
    
    /*
    for ( int i = 0; i < fitParams.getSize(); ++i)
    {
      auto & fitPar = (RooRealVar &) fitParams[i];
      std::cout << fitPar.GetName() << " " << fitPar.getVal() << " +- " << fitPar.getError() << std::endl;
    }
    */
    
    mass_pdf->plotOn(myPlot2_A, Name("pdfMASS_tot"), LineColor(kBlack), LineWidth(2));
    if (nSignalGaussComponents >= 1) {
        mass_pdf->plotOn(myPlot2_A, Components(signal_mass_gaus), LineColor(kBlue+1), LineStyle(kDashed), LineWidth(2), Name("signal_gaus_component"));
    }
    if (nSignalCBComponents >= 1) {
        mass_pdf->plotOn(myPlot2_A, Components(*signal_mass_cb), LineColor(kRed+1), LineStyle(kDashed), LineWidth(2), Name("signal_cb_component"));
    }
    if (nSignalGaussComponents >= 2) {
        mass_pdf->plotOn(myPlot2_A, Components(signal_mass_gaus2), LineColor(kMagenta+2), LineStyle(kDashed), LineWidth(2), Name("signal_gaus2_component"));
    }

    // Determine y-axis range from the plotted data points and their errors.
    myPlot2_A->SetFillStyle(4000);
    myPlot2_A->GetYaxis()->SetTitle("Events");
    myPlot2_A->GetYaxis()->SetTitleOffset(1.6);

    double yMax = -1.0;
    double yMinPositive = 1e99;
    RooHist *dataHistForRange = dynamic_cast<RooHist*>(myPlot2_A->getHist("dataOS"));
    if (dataHistForRange) {
        double xPoint = 0.0, yPoint = 0.0;
        for (int i = 0; i < dataHistForRange->GetN(); ++i) {
            dataHistForRange->GetPoint(i, xPoint, yPoint);
            const double yHigh = yPoint + dataHistForRange->GetErrorYhigh(i);
            const double yLow = yPoint - dataHistForRange->GetErrorYlow(i);
            if (std::isfinite(yHigh)) yMax = std::max(yMax, yHigh);
            if (std::isfinite(yLow) && yLow > 0.0) yMinPositive = std::min(yMinPositive, yLow);
            else if (std::isfinite(yPoint) && yPoint > 0.0) yMinPositive = std::min(yMinPositive, yPoint);
        }
    }
    if (yMax <= 0.0 || yMinPositive >= 1e99) {
        TH1* hRange = ws->data("dsAB")->createHistogram("hist_yaxis_range", *ws->var("mass"), Binning(massPlotBins,myPlot_A->GetXaxis()->GetXmin(),myPlot_A->GetXaxis()->GetXmax()));
        if (hRange) {
            yMax = hRange->GetBinContent(hRange->GetMaximumBin());
            for (int i=1; i<=hRange->GetNbinsX(); ++i) {
                const double yBin = hRange->GetBinContent(i);
                if (yBin > 0.0) yMinPositive = std::min(yMinPositive, yBin);
            }
        }
    }
    if (yMax <= 0.0) yMax = 1.0;
    if (yMinPositive >= 1e99) yMinPositive = 0.5;

    if (logY_flag) {
        const double yLow = std::max(0.1, 0.45 * yMinPositive);
        const double logSpan = std::log10(std::max(10.0, yMax / yLow));
        const double yHigh = yLow * std::pow(10.0, logSpan * 1.35);
        myPlot2_A->GetYaxis()->SetRangeUser(yLow, yHigh);
    }
    else {
        myPlot2_A->GetYaxis()->SetRangeUser(0.0, yMax * 1.65);
    }

    myPlot2_A->GetXaxis()->SetLabelSize(0);
    myPlot2_A->GetXaxis()->SetTitleSize(0);
    myPlot2_A->GetXaxis()->CenterTitle();
    
    myPlot2_A->GetXaxis()->SetTitle("");
    myPlot2_A->Draw("e");
    
    auto findObj = [&](RooPlot *fr, const char *n) -> TObject * { return fr ? fr->findObject(n) : nullptr; };
    TLegend leg(0.50, 0.70, 0.74, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    if (auto *o = findObj(myPlot2_A, "dataOS")) leg.AddEntry(o, "Data", "lep");
    if (auto *o = findObj(myPlot2_A, "pdfMASS_tot")) leg.AddEntry(o, "Fit", "l");
    if (nSignalGaussComponents >= 1) if (auto *o = findObj(myPlot2_A, "signal_gaus_component")) leg.AddEntry(o, "G1", "l");
    if (nSignalCBComponents >= 1) if (auto *o = findObj(myPlot2_A, "signal_cb_component")) leg.AddEntry(o, nSignalCBComponents >= 2 ? "DSCB" : "CB1", "l");
    if (nSignalGaussComponents >= 2) if (auto *o = findObj(myPlot2_A, "signal_gaus2_component")) leg.AddEntry(o, "G2", "l");
    leg.Draw("same");

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
        tx.SetTextSize(0.04);
        tx.SetTextFont(72);
        tx.DrawLatex(0.19, 0.935, "CMS Internal");
    }
    {
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.03);
        tx.SetTextFont(42);
        double xtext = 0.19, y0 = 0.865, dy = -0.05;
        int k = 0;
        tx.DrawLatex(xtext, y0 + dy * k++, Form("%s J/#psi MC", bCont.Data()));
        if (yLow == 0) tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
        else tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
    }
    if (!publish) {
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.03);
        tx.SetTextFont(42);
        const int status = fitMass ? fitMass->status() : -1;
        int hesse = -1;
        if (fitMass) {
            for (UInt_t i = 0, n = fitMass->numStatusHistory(); i < n; ++i) {
                const char *lab = fitMass->statusLabelHistory(i);
                if (lab && TString(lab) == "HESSE") {
                    hesse = fitMass->statusCodeHistory(i);
                    break;
                }
            }
        }
        if (status != 0 || hesse != 0) {
            tx.DrawLatex(0.19, 0.765, Form("Status : MINIMIZE=%d HESSE=%d", status, hesse));
        }
    }
    if (!publish) {
        TLatex tp;
        tp.SetNDC();
        tp.SetTextSize(0.024);
        tp.SetTextFont(42);
        double xtext = 0.74, y0 = 0.87, dy = -0.045;
        int k = 0;
        auto printVar = [&](const char *title, const RooAbsReal &var) {
            const RooRealVar *rrv = dynamic_cast<const RooRealVar *>(&var);
            if (rrv && rrv->isConstant()) {
                tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, rrv->getVal()));
            }
            else if (rrv) {
                tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, rrv->getVal(), rrv->getError()));
            }
            else {
                double err = fitMass ? var.getPropagatedError(*fitMass) : 0.0;
                if (std::isfinite(err) && err > 0.0)
                    tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, var.getVal(), err));
                else
                    tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, var.getVal()));
            }
        };
        printVar("N_{sig}", Nsig);
        printVar("#mu", signal_mass_mean);
        if (nSignalGaussComponents >= 1) printVar("#sigma_{G1}", signal_mass_sigma);
        if (nSignalCBComponents >= 1) {
            printVar("#sigma_{CB1}", *signal_mass_cb_sigma);
            printVar("#alpha_{CB1}", signal_mass_cb_alpha);
            printVar("n_{CB1}", signal_mass_cb_n);
        }
        if (nSignalGaussComponents >= 2) printVar("#sigma_{G2}", signal_mass_sigma2);
        if (nSignalCBComponents >= 2) {
            printVar("#sigma_{CB2}", signal_mass_cb_sigma2);
            printVar("#alpha_{CB2}", signal_mass_cb_alpha2);
            printVar("n_{CB2}", signal_mass_cb_n2);
        }
        if (signalComponents.size() >= 2) printVar("r_{f1}", signal_mass_frac_ratio1);
        if (signalComponents.size() >= 3) printVar("r_{f2}", signal_mass_frac_ratio2);
        if (signalComponents.size() >= 4) printVar("r_{f3}", signal_mass_frac_ratio3);
    }
    
	// ------------------------------------------------------------------
	// draw pull and print summary
	// ------------------------------------------------------------------
    TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0.0, 0.0, 1.0, 0.25);
    c_A->cd();
    pad_A_2->Draw();
    pad_A_2->cd();
    pad_A_2->SetTopMargin(0.00001);
    pad_A_2->SetBottomMargin(0.4);
    pad_A_2->SetFillStyle(4000);
    pad_A_2->SetFrameFillStyle(4000);
    pad_A_2->SetTicks(1,1);

    RooPlot* frameTMP = (RooPlot*)myPlot2_A->Clone("TMP");
    RooHist* hpull_A = frameTMP->pullHist("dataOS","pdfMASS_tot", true);
    hpull_A->SetMarkerSize(0.8);
    RooPlot* pullFrame_A = ws->var("mass")->frame(Title(""));
    pullFrame_A->addPlotable(hpull_A,"P");
    pullFrame_A->GetYaxis()->SetTitle("Pull");
    pullFrame_A->GetXaxis()->SetTitle("m_{#mu#mu} [GeV/c^{2}]");
    pullFrame_A->GetXaxis()->CenterTitle();
    pullFrame_A->SetMinimum(-8);
    pullFrame_A->SetMaximum(8);
    pullFrame_A->GetYaxis()->SetNdivisions(505);
    pullFrame_A->GetYaxis()->SetTitleSize(0.12);
    pullFrame_A->GetYaxis()->SetLabelSize(0.10);
    pullFrame_A->GetXaxis()->SetTitleSize(0.15);
    pullFrame_A->GetXaxis()->SetLabelSize(0.10);
    pullFrame_A->Draw();

    TLine *l1 = new TLine(massLow,0,massHigh,0);
    l1->SetLineStyle(2);
    l1->Draw("same");

    auto chiM = chi2_from_pull(hpull_A);
    int npar = fitMass ? fitMass->floatParsFinal().getSize() : 0;
    int ndf = std::max(1, chiM.second - npar);
    double pvalue = TMath::Prob(chiM.first, ndf);
    {
        TLatex tc;
        tc.SetNDC();
        tc.SetTextSize(0.085);
        tc.SetTextFont(42);
        tc.SetTextAlign(33);
        tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d", chiM.first, ndf));
    }

    TH1D* outh = new TH1D("fitResults","fit result",20,0,20);
    outh->GetXaxis()->SetBinLabel(1,"Jpsi");

    float temp1 = Nsig.getVal();
    float temp1err = Nsig.getError();

    outh->SetBinContent(1,temp1);
    outh->SetBinError(1,temp1err);

    Double_t theNLL = fitMass->minNll();
    cout << " *** NLL : " << theNLL << endl;
    //ws->Print();
    RooArgSet* fitargs = new RooArgSet();
    fitargs->add(fitMass->floatParsFinal());
    RooDataSet *datasetMass = new RooDataSet("datasetMass","dataset with Mass Fit result", *fitargs );
    datasetMass->add(*fitargs);
    c_A->Update();
    c_A->Print(figName("mass_fit"));

    TFile* outFile = nullptr;
    if (!drawFromSavedFit) {
	outFile = new TFile(modelFileName, "RECREATE");
    outFile->cd();
    TParameter<int>("nSignalCBComponents", nSignalCBComponents).Write();
    TParameter<int>("nSignalGaussComponents", nSignalGaussComponents).Write();
    TParameter<int>("fit_status", fitMass ? fitMass->status() : -1).Write();
    TParameter<int>("fit_covQual", fitMass ? fitMass->covQual() : -1).Write();
    TParameter<double>("fit_minNll", fitMass ? fitMass->minNll() : 0.0).Write();
    auto writeVar = [&](const char *name, const RooRealVar &var) {
        TParameter<double>(name, var.getVal()).Write();
        TParameter<double>(TString(name) + "_err", std::max(var.getError(), 1e-9)).Write();
    };
    auto writeAbsReal = [&](const char *name, const RooAbsReal &var) {
        double err = fitMass ? std::max(var.getPropagatedError(*fitMass), 1e-9) : 1e-9;
        TParameter<double>(name, var.getVal()).Write();
        TParameter<double>(TString(name) + "_err", err).Write();
    };
    writeVar("Nsig", Nsig);
    writeVar("signal_mass_mean", signal_mass_mean);
    if (nSignalGaussComponents >= 1)
        writeVar("signal_mass_sigma", signal_mass_sigma);
    if (nSignalCBComponents >= 1) {
        if (nSignalGaussComponents >= 1)
            writeAbsReal("signal_mass_cb_sigma_delta", *signal_mass_cb_sigma_delta);
        writeAbsReal("signal_mass_cb_sigma", *signal_mass_cb_sigma);
        writeVar("signal_mass_cb_alpha", signal_mass_cb_alpha);
        writeVar("signal_mass_cb_n", signal_mass_cb_n);
    }
    if (nSignalGaussComponents >= 2) {
        writeAbsReal("signal_mass_sigma_delta2", signal_mass_sigma_delta2);
        writeAbsReal("signal_mass_sigma2", signal_mass_sigma2);
    }
    if (nSignalCBComponents >= 2) {
        writeAbsReal("signal_mass_cb_sigma_delta2", signal_mass_cb_sigma_delta2);
        writeAbsReal("signal_mass_cb_sigma2", signal_mass_cb_sigma2);
        writeVar("signal_mass_cb_alpha2", signal_mass_cb_alpha2);
        writeVar("signal_mass_cb_n2", signal_mass_cb_n2);
    }
    if (signalComponents.size() >= 2) {
        writeVar("signal_mass_frac_ratio1", signal_mass_frac_ratio1);
        writeAbsReal("signal_mass_frac1", signal_mass_frac1);
    }
    if (signalComponents.size() >= 3) {
        writeVar("signal_mass_frac_ratio2", signal_mass_frac_ratio2);
        writeAbsReal("signal_mass_frac2", signal_mass_frac2);
    }
    if (signalComponents.size() >= 4) {
        writeVar("signal_mass_frac_ratio3", signal_mass_frac_ratio3);
        writeAbsReal("signal_mass_frac3", signal_mass_frac3);
    }
    if (fitMass) fitMass->Write("fit_result");
    mass_pdf->Write("mass_pdf");
    signal_mass_pdf->Write("signal_pdf");
    datasetMass->Write("fit_parameters_dataset");
    outh->Write("fitResults");
    ws->Write("workspace");

    outFile->Close();
    }
    else {
        std::cout << "[PlotOnly] Left ROOT model file unchanged: " << modelFileName << std::endl;
    }
    std::cout << Form("MC mass chi2/ndf: %.1f/%d = %.3f, p=%.3g", chiM.first, ndf, ndf > 0 ? chiM.first / ndf : 0.0, pvalue) << std::endl;
    std::cout << "------------------ FIT RESULT FOR MASS ONLY --------------" << std::endl;
    fitMass->Print("v");
}