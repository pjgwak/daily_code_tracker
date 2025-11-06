#ifndef check1SLimits_C
#define check1SLimits_C

#include "TString.h"
#include "TH1D.h"
#include "TLine.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"

#include <vector>
#include <map>

#include "../Fitter/Macros/Utilities/bin.h"
#include "../Fitter/Macros/Utilities/resultUtils.h"
#include "../Fitter/Systematics/syst.h"
#include "limits.h"

using namespace std;

void check1SLimits(
                   const char* workDir, // workDir: usual tag where to look for files in Output
                   const char* lFileName="cLimits_683_NominalABCD_Asym_2SPL_woSyst.csv", // file name to save limits results
                   bool dosyst = false,
                   int mode = 1, // mode=0 -> pass, mode=1 -> prompt, mode=2 -> nonprompt
                   const char* workDirFail=""
)
{
  TString slFileName(lFileName);
  if ( dosyst && !slFileName.Contains("wSys") )
  {
    cout << "Comparison requires systematics but limits file does not contain them" << endl;
    return;
  }
  
  // list of files
  set<anabin> thebins = allbins();
  const char* ppp = "../Fitter";
  
  // systematic uncertainties for fit
  map<anabin, syst> syst_All;
  if ( dosyst )
  {
     if (mode==0) syst_All = readSyst_all_pass("",ppp,workDir);
     if (mode==1) syst_All = readSyst_all_prompt("",ppp,workDir,workDirFail);
     if (mode==2) syst_All = readSyst_all_nonprompt("",ppp,workDir,workDirFail);
  }
  
  // bin edges
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  
  // histo for 1sigma limits checks
  TH1* hCL = new TH1D("hOneSigmaCLComparison","",thebins.size(),0,thebins.size());
  hCL->GetYaxis()->SetTitle("CL_{1#sigma}/#sigma");
  hCL->GetYaxis()->SetTitleOffset(1.15);
  hCL->SetStats(0);
  hCL->SetDirectory(0);
  hCL->SetMarkerColor(1);
  hCL->SetMarkerStyle(20);
  hCL->SetMarkerSize(1);
  hCL->SetLineColor(1);
  
  TLine* l1 = new TLine(0.,1.,hCL->GetXaxis()->GetXmax(),1.);
  l1->SetLineWidth(3);
  
  hCL->GetListOfFunctions()->Add(l1);

  map<anabin,limits> maplim = readLimits(Form("csv/%s",slFileName.Data()));
  
  int cnt=1;
  for (set<anabin>::const_iterator it=thebins.begin(); it!=thebins.end(); it++)
  {
     cout << "Checking 1 sigma limits for analysis bin " << cnt << endl;

     anabin thebin = *it;
     ptmin = thebin.ptbin().low();
     ptmax = thebin.ptbin().high();
     ymin = thebin.rapbin().low();
     ymax = thebin.rapbin().high();
     centmin = thebin.centbin().low();
     centmax = thebin.centbin().high();

     double sigmaDoubleR = 0;
     double doubleR = 0;
     if (mode==0) {
        doubleR = doubleratio_pass_nominal(workDir,thebin,ppp);
        sigmaDoubleR = doubleratio_pass_stat(workDir,thebin,ppp);
     }
     if (mode==1) {
        doubleR = doubleratio_prompt_nominal(workDir,workDirFail,thebin,ppp);
        sigmaDoubleR = doubleratio_prompt_stat(workDir,workDirFail,thebin,ppp);
     }
     if (mode==2) {
        doubleR = doubleratio_nonprompt_nominal(workDir,workDirFail,thebin,ppp);
        sigmaDoubleR = doubleratio_nonprompt_stat(workDir,workDirFail,thebin,ppp);
     }

     double systAll=0;
     if ( dosyst )
     {
        systAll = syst_All[thebin].value_dR;
        sigmaDoubleR = sqrt(pow(sigmaDoubleR,2)+pow(systAll,2));
     }

     limits lim = maplim[thebin];

     TString binName(Form("Pt[%.1f,%.1f]-Y[%.1f,%.1f]-C[%.1f,%.1f]",ptmin,ptmax,ymin,ymax,centmin,centmax));

     double comp = -1.;
     if ( sigmaDoubleR != 0 ) comp = (lim.val.second-lim.val.first)/(2.*sigmaDoubleR);
     hCL->SetBinContent(cnt,comp);
     hCL->GetXaxis()->SetBinLabel(cnt,binName.Data());

     cnt++;
  } // loop on the files

  TFile* fSave = new TFile("oneSigmaCLComparison.root","RECREATE");
  
  TCanvas* c = new TCanvas("cOneSigmaCLComparison","",90,116,1265,535);
  c->Range(-3.690909,-0.01066472,33.30606,0.01252061);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetRightMargin(0.1163896);
  c->SetTopMargin(0.03732809);
  c->SetBottomMargin(0.1630648);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  gPad->SetGridx();
  gPad->SetGridy();
  hCL->Draw("p");
  
  c->Write("cOneSigmaCLComparison", TObject::kOverwrite | TObject::kSingleKey);
  fSave->Close(); delete fSave;
  
}

#endif // #ifndef check1SLimits_C
