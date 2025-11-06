#ifndef compareCtau_C
#define compareCtau_C


#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"

#include <vector>
#include <map>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>

#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/bfrac.h"
#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Systematics/syst.h"

using namespace std;

void setCustomTDRStyle();

void compareCtau(TString fileNameMC, // Has to be of the form "dir1"
                 TString fileNameData, // Has to be of the form "dir1"
                 const char* prependPath="" //="",
)
{
   // Examples of usage:
   //      1) compareCtau("MCFitsCent_nonprompt_ctauRes/DataFits_006,MCFitsCent_prompt_ctauRes/DataFits_006,DataFitsCent_ctauRes_nominalErr/DataFits_006","pdfCTAUNRES_Tot","J/#psi nonprompt MC,J/#psi prompt MC, J/#psi DATA","ctauNRes","MCJPSINOPR,MCJPSIPR,DATA","","ctauRes",0)
   //      2) compareCtau("DataFitsCent_ctauErr_nominal/DataFits_006,DataFitsCent_ctauErr_CBG/DataFits_006","pdfCTAUERR_Jpsi","J/#psi CB+CB MC,J/#psi CB+G","ctauErr","DATA,DATA","","ctauErr",0)

   bool setLogScale = false;

   TString svarName("ctau");

   double maxYrange(1e5);

   TString sTitle("");
   TString sXaxis("");
   double minRange, maxRange;
   minRange = -1;
   maxRange = 4;
   sTitle = "l_{J/#psi} resolution PDFs comparison";
   //    maxYrange = 0.06;
   maxYrange = 1e5;
   sXaxis = "#frac{#font[12]{l}_{J/#psi}(-#font[12]{l}_{J/#psi}^{true})}{#font[12]{#sigma}_{J/#psi}}";
   setLogScale = true;

   map<anabin,TString> theFilesData,theFilesMC;
   if (fileNameMC.Contains(".root")) {
      anabin thebin = binFromFile(fileNameMC);
      theFilesData[thebin] = fileNameData;
      theFilesMC[thebin] = fileNameMC;
   } else {
      vector<TString> theFiles = fileList(fileNameData,"","DATA","","ctauMass");

      for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
         anabin thebin = binFromFile(*it);
         theFilesData[thebin] = *it;
      }

      theFiles = fileList(fileNameMC,"","MCJPSINOPR","","ctauReco");

      for (vector<TString>::const_iterator it=theFiles.begin(); it!=theFiles.end(); it++) {
         anabin thebin = binFromFile(*it);
         theFilesMC[thebin] = *it;
      }
   }


   // loop on the bins
   for (map<anabin,TString>::const_iterator it=theFilesData.begin(); it!=theFilesData.end(); it++) {
      anabin thebin = it->first;
      TString itData = it->second;
      TString itMC = theFilesMC[thebin];

      thebin.print();
      cout << itData << " " << itMC << endl;

      bool isPbPb      = itData.Contains("PbPb");
      TString collSystem = isPbPb ? "PbPb" : "PP";
      TString spdfName_ctauRes  = "pdfCTAURES_Jpsi_" + collSystem;
      TString spdfName_ctauTrue = "pdfCTAUCOND_JpsiNoPR_" + collSystem;
      TString sdataName         = "dOS_MCJPSINOPR_" + collSystem + "_NoBkg";

      // bin edges
      float ptmin, ptmax, ymin, ymax, centmin, centmax;


      // parse the file name to get info
      ptmin = thebin.ptbin().low();
      ptmax = thebin.ptbin().high();
      ymin = thebin.rapbin().low();
      ymax = thebin.rapbin().high();
      centmin = thebin.centbin().low();
      centmax = thebin.centbin().high();

      if (!isPbPb)
      {
         centmin = 0;
         centmax = 200;
         isPbPb = false;
      }


      RooPlot* frame(NULL);

      // plot the reconstructed ctau first (from nonprompt MC)
      TFile *f = TFile::Open(itMC); RooWorkspace *ws = NULL;
      if (!f)
      {
         cout << "[ERROR] file " << itMC << " does not exist." << endl;
      }
      else
      {
         ws = (RooWorkspace*) f->Get("workspace");
         if (!ws)
         {
            cout << "[ERROR] workspace not found in " << itMC << "." << endl;
         }
      }

      if (f && ws)
      {
         // get the dataset
         if (!ws->data(sdataName))
         {
            cout << "[ERROR] The dataset " << sdataName << " does not exist in the file " << itMC << endl;
            return;
         }

         frame = ws->var("ctau")->frame(Bins(100), Range(minRange, maxRange));

         ws->data(sdataName)->plotOn(frame,Name("ctau_reco"),LineColor(kBlack));
      }

      // now plot the convoluted pdf, ctau_true by the resolution (from the data fit)
      f = TFile::Open(itData); ws = NULL;
      if (!f)
      {
         cout << "[ERROR] file " << itData << " does not exist." << endl;
      }
      else
      {
         ws = (RooWorkspace*) f->Get("workspace");
         if (!ws)
         {
            cout << "[ERROR] workspace not found in " << itData << "." << endl;
         }
      }

      if (f && ws)
      {
         // do not forget to read the initial parameters here
         ws->loadSnapshot(Form("pdfMASS_Tot_%s_parIni",collSystem.Data()));

         RooRealVar *ctau = ws->var("ctau");
         ctau->setBins(10000,"cache") ;
         // RooAbsPdf *pCtauRes = ws->pdf(spdfName_ctauRes);
         RooAbsPdf *pCtauTrue = ws->pdf(spdfName_ctauTrue);
         if (!pCtauTrue)
         {
            cout << "[ERROR] The pdf does not exist in the file " << itData << endl;
            return;
         }
         ws->loadSnapshot(Form("pdfCTAUMASS_Tot_%s_parIni",collSystem.Data()));

         pCtauTrue->plotOn(frame,Name("ctau_reco"),LineColor(kBlack), Precision(1e-5));
      }

      // set the CMS custom style
      setTDRStyle();

      frame->SetTitle("");
      frame->GetXaxis()->SetTitle(sXaxis.Data());
      frame->GetXaxis()->CenterTitle(kTRUE);
      frame->GetXaxis()->SetLabelSize(0.04);
      frame->GetXaxis()->SetTitleSize(0.035);
      frame->GetXaxis()->SetTitleOffset(1.3);
      frame->GetXaxis()->SetTitleFont(42);
      frame->GetYaxis()->SetTitle("Prob. distribution");
      frame->GetYaxis()->SetLabelSize(0.04);
      frame->GetYaxis()->SetTitleSize(0.04);
      frame->GetYaxis()->SetTitleOffset(1.8);
      frame->GetYaxis()->SetTitleFont(42);
      frame->GetYaxis()->SetRangeUser(1E-3,maxYrange);

      TCanvas cFig("cFig","cFig",800,800);
      cFig.SetLogy(setLogScale);
      frame->Draw();

      // Drawing the text in the plot
      TLatex t; t.SetNDC(); t.SetTextSize(0.032);
      float dy = 0;

      t.SetTextSize(0.03);
      t.DrawLatex(0.35, 0.95, sTitle.Data());
      t.DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",ptmin,ptmax)); dy+=0.045;
      t.DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",ymin,ymax)); dy+=0.045;
      if (isPbPb) {t.DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(centmin/2), (int)(centmax/2))); dy+=0.045;}
      t.DrawLatex(0.21, 0.86-dy, isPbPb ? "PbPb" : "PP"); dy+=0.045;

      // Drawing the Legend
      // TLegend* leg = new TLegend(0.6303, 0.7338, 0.8308, 0.8695); leg->SetTextSize(0.03);
      // for (int i=0; i<dirCnt; i++) {
      //   if(frame->findObject(Form("PDF_%d",i))) { leg->AddEntry(frame->findObject(Form("PDF_%d",i)),pdfTitle.at(i).Data(),"l"); }
      //   else if (frame2)
      //   {
      //     if (frame2->findObject(Form("PDF_%d",i))) { leg->AddEntry(frame2->findObject(Form("PDF_%d",i)),pdfTitle.at(i).Data(),"l"); }
      //   }
      // }
      // leg->Draw("same");

      // cFig->SaveAs("ctauPlot.pdf");

      // Save the plot in different formats
      TString svarName = "ctau";
      TString saveDir = "NPResCheck";

      gSystem->mkdir(Form("pdfComparison/%s/%s/root/", svarName.Data(), saveDir.Data()), kTRUE);
      cFig.SaveAs(Form("pdfComparison/%s/%s/root/PDFComparison_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", svarName.Data(), saveDir.Data(), svarName.Data(), spdfName_ctauTrue.Data(), collSystem.Data(), (ptmin*10.0), (ptmax*10.0), (ymin*10.0), (ymax*10.0), (int)centmin, (int)centmax));

      gSystem->mkdir(Form("pdfComparison/%s/%s/pdf/", svarName.Data(), saveDir.Data()), kTRUE);
      cFig.SaveAs(Form("pdfComparison/%s/%s/pdf/PDFComparison_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", svarName.Data(), saveDir.Data(), svarName.Data(), spdfName_ctauTrue.Data(), collSystem.Data(), (ptmin*10.0), (ptmax*10.0), (ymin*10.0), (ymax*10.0), (int)centmin, (int)centmax));
   } // loop on bins

}

#endif // #define compareCtau_C
