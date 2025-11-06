#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Utilities/resultUtils.h"

using namespace std;
using namespace RooFit;

void DiagnosisMacro_all(const char* workDirName, const char* DSTag="DATA", int Nbins = 10, int Nsigma = 10, int CPUused = 1);

int DiagnosisMacro(int Nbins = 10, int Nsigma = 10, int CPUused = 1, TString Filename = "FIT_DATA_Psi2SJpsi_PPPrompt_Bkg_SecondOrderChebychev_pt65300_rap016_cent0200_262620_263757.root", TString Outputdir = "./") 
//Nbins: Number of points for which to calculate profile likelihood. Time required is about (1/CPU) minutes per point per parameter. 0 means do plain likelihood only
//Nsigma: The range in which the scan is performed (value-Nsigma*sigma, value+Nsigma*sigma)
//CPUused: anything larger than 1 causes weird fit results on my laptop, runs fine on lxplus with more (16)

{
	// R e a d   w o r k s p a c e   f r o m   f i l e
	// -----------------------------------------------
	// Open input file with workspace

	//Filename = "FIT_DATA_Psi2SJpsi_PP_Jpsi_DoubleCrystalBall_Psi2S_DoubleCrystalBall_Bkg_Chebychev2_pt6590_rap016_cent0200.root";
	//Filename = "FIT_DATA_Psi2SJpsi_PbPb_Jpsi_DoubleCrystalBall_Psi2S_DoubleCrystalBall_Bkg_Chebychev1_pt6590_rap016_cent0200.root";

	TFile *f = new TFile(Filename);
	// Retrieve workspace from file
	RooWorkspace* w = (RooWorkspace*)f->Get("workspace");

	// Retrieve x,model and data from workspace
	RooRealVar* x = w->var("invMass");
	RooAbsPdf* model = w->pdf("simPdf_syst");
	if (model == 0){ model = w->pdf("simPdf"); }
	if (model == 0){ model = w->pdf("pdfMASS_Tot_PP"); }
	if (model == 0){ model = w->pdf("pdfMASS_Tot_PbPb"); }
	if (model == 0){
		cout << "[ERROR] pdf failed to load from the workspace" << endl; return false;
	}

	RooAbsData* data = w->data("dOS_DATA");
	if (data == 0){ data = w->data("dOS_DATA_PP"); }
	if (data == 0){ data = w->data("dOS_DATA_PbPb"); }
	if (data == 0){
		cout << "[ERROR] data failed to load from the workspace" << endl; return false;
	}

	// Print structure of composite p.d.f.
	model->Print("t");

	/*
	// P l o t   m o d e l 
	// ---------------------------------------------------------
	// Plot data and PDF overlaid
	RooPlot* xframe = x->frame(Title("J/psi Model and Data"));
	data->plotOn(xframe);
	model->plotOn(xframe);

	// Draw the frame on the canvas
	TCanvas* c2 = new TCanvas("PlotModel", "PlotModel", 1000, 1000);
	gPad->SetLeftMargin(0.15);
	xframe->GetYaxis()->SetTitleOffset(2.0);
	xframe->Draw();//*/

	///// Check parameters

	RooArgSet* paramSet1 = model->getDependents(data);
	paramSet1->Print("v");  // Just check
	RooArgSet* paramSet2 = model->getParameters(data);
	paramSet2->Print("v");
	int Nparams = paramSet2->getSize();
	cout << "Number of parameters: " << Nparams<<endl<<endl;


	// C o n s t r u c t   p l a i n   l i k e l i h o o d
	// ---------------------------------------------------
	// Construct unbinned likelihood
	RooAbsReal* nll = model->createNLL(*data, NumCPU(CPUused));
	// Minimize likelihood w.r.t all parameters before making plots
	RooMinuit(*nll).migrad();


	//////////////////////////////////////////////////////

	///////////////////   L O O P   O V E R   P A R A M E T E R S

	/////////////////////////////////////////////////////
	
	/// Set up loop over parameters
	TString ParamName;
	double ParamValue;
	double ParamError;
	double ParamLimitLow;
	double ParamLimitHigh;
	double FitRangeLow;
	double FitRangeHigh;
	RooRealVar* vParam;
	int counter = 0;

	// Loop start
	TIterator* iter = paramSet2->createIterator();
	TObject* var = iter->Next();
	while (var != 0) {
		counter++;
		ParamName = var->GetName();
		vParam = w->var(ParamName);
		ParamValue = vParam->getVal();
		ParamError = vParam->getError();
		ParamLimitLow = vParam->getMin();
		ParamLimitHigh = vParam->getMax();
		cout << ParamName << " has value " << ParamValue << " with error: " << ParamError << " and limits: " << ParamLimitLow << " to " << ParamLimitHigh << endl << endl;

		if (ParamError == 0){   //Skipping fixed parameters
			cout << "Parameter was fixed, skipping its fitting" << endl;
			cout << endl << "DONE WITH " << counter << " PARAMETER OUT OF " << Nparams << endl << endl;
			var = iter->Next();
			continue;
		}

		// determining fit range: Nsigma sigma on each side unless it would be outside of parameter limits
		if ((ParamValue - Nsigma * ParamError) > ParamLimitLow){ FitRangeLow = (ParamValue - Nsigma * ParamError); }
		else { FitRangeLow = ParamLimitLow; }

		if ((ParamValue + Nsigma * ParamError) < ParamLimitHigh){ FitRangeHigh = (ParamValue + Nsigma * ParamError); }
		else { FitRangeHigh = ParamLimitHigh; }


		// P l o t    p l a i n   l i k e l i h o o d   a n d   C o n s t r u c t   p r o f i l e   l i k e l i h o o d
		// ---------------------------------------------------
		RooPlot* frame1;
		RooAbsReal* pll=NULL;

		if (Nbins != 0){
			frame1 = vParam->frame(Bins(Nbins), Range(FitRangeLow, FitRangeHigh), Title(TString::Format("LL and profileLL in %s", ParamName.Data())));
			nll->plotOn(frame1, ShiftToZero());

			pll = nll->createProfile(*vParam);
			// Plot the profile likelihood
			pll->plotOn(frame1, LineColor(kRed), RooFit::Precision(-1));
		}
		else{  //Skip profile likelihood
			frame1 = vParam->frame(Bins(10), Range(FitRangeLow, FitRangeHigh), Title(TString::Format("LL and profileLL in %s", ParamName.Data())));
			nll->plotOn(frame1, ShiftToZero());
		}

		// D r a w   a n d   s a v e   p l o t s
		// -----------------------------------------------------------------------

		// Adjust frame maximum for visual clarity
		frame1->SetMinimum(0);
		frame1->SetMaximum(20);

		TCanvas* c = new TCanvas("CLikelihoodResult", "CLikelihoodResult", 800, 600);
		c->cd(1);
		gPad->SetLeftMargin(0.15); 
		frame1->GetYaxis()->SetTitleOffset(1.4);
		frame1->Draw();
		TLegend* leg = new TLegend(0.70, 0.70, 0.95, 0.88, "");
		leg->SetFillColor(kWhite);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.035);
		TLegendEntry *le1 = leg->AddEntry(nll, "Plain likelihood", "l");
		le1->SetLineColor(kBlue);
		le1->SetLineWidth(3);
		TLegendEntry *le2 = leg->AddEntry(pll, "Profile likelihood", "l");
		le2->SetLineColor(kRed);
		le2->SetLineWidth(3);
		leg->Draw("same");

      //Save plot
      TString StrippedName = TString(Filename(Filename.Last('/')+1,Filename.Length()));
      StrippedName = StrippedName.ReplaceAll(".root","");
      cout << StrippedName << endl;
      gSystem->mkdir(Form("%s/root/%s", Outputdir.Data(), StrippedName.Data()), kTRUE); 
      c->SaveAs(Form("%s/root/%s/Likelihood_scan_%s.root", Outputdir.Data(), StrippedName.Data(), ParamName.Data()));
      gSystem->mkdir(Form("%s/pdf/%s", Outputdir.Data(), StrippedName.Data()), kTRUE); 
      c->SaveAs(Form("%s/pdf/%s/Likelihood_scan_%s.pdf", Outputdir.Data(), StrippedName.Data(), ParamName.Data()));
      gSystem->mkdir(Form("%s/png/%s", Outputdir.Data(), StrippedName.Data()), kTRUE); 
      c->SaveAs(Form("%s/png/%s/Likelihood_scan_%s.png", Outputdir.Data(), StrippedName.Data(), ParamName.Data()));


		delete c;
		delete frame1;
		if (pll) delete pll;

		cout << endl << "DONE WITH " << counter << " PARAMETER OUT OF " << Nparams << endl << endl;
		//if (counter == 2){ break; } //Exit - for testing
		var = iter->Next();
	}  // End of the loop

	return true;
}

void DiagnosisMacro_all(const char* workDirName, const char* DSTag, int Nbins, int Nsigma, int CPUused) {
   vector<TString> theFiles = fileList(workDirName,"",DSTag);
   TString Outputdir = Form("Output/%s/plot/RESULT/", workDirName); 

   for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
      DiagnosisMacro(Nbins, Nsigma, CPUused, *it, Outputdir);
   }
}

