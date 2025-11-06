#ifndef comparePDFs_C
#define comparePDFs_C


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

void comparePDFs(const char* workDirNames, // Has to be of the form "dir1,dir2,dir3,..."
                 const char* pdfName, // Name of the pdf to plot
                 const char* pdfTitles, // Titles of the PDFs for the legend
                 const char* varName, // Name of the x-axis var of the form
                 const char* DSTags, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ... // Has to be of the form "DSTagdir1,DSTagdir2,DSTagdir3,..."
                 const char* prependPath, //="",
                 const char* fitType, // "mass", "ctau"...
                 bool wantPureSMC //=false
)
{
  // Examples of usage:
  //      1) comparePDFs("MCFitsCent_nonprompt_ctauRes/DataFits_006,MCFitsCent_prompt_ctauRes/DataFits_006,DataFitsCent_ctauRes_nominalErr/DataFits_006","pdfCTAUNRES_Tot","J/#psi nonprompt MC,J/#psi prompt MC, J/#psi DATA","ctauNRes","MCJPSINOPR,MCJPSIPR,DATA","","ctauRes",0)
  //      2) comparePDFs("DataFitsCent_ctauErr_nominal/DataFits_006,DataFitsCent_ctauErr_CBG/DataFits_006","pdfCTAUERR_Jpsi","J/#psi CB+CB MC,J/#psi CB+G","ctauErr","DATA,DATA","","ctauErr",0)
  
  bool isPbPb = false;
  bool setLogScale = false;
  
  TString spdfName(pdfName);
  TString svarName(varName);
  
  double maxYrange(1.);
  
  TString sTitle("");
  TString sXaxis("");
  double minRange, maxRange;
  if (svarName.Contains("mass"))
  {
    minRange = 2.6;
    maxRange = 3.5;
    sTitle = "Invariant mass PDFs comparison";
    sXaxis = "m_{#mu#mu} (GeV/c)";
    setLogScale = true;
  }
  else if (svarName.Contains("ctauErr"))
  {
    minRange = 0.;
    maxRange = 0.15;
    sTitle = "l_{J/#psi} error PDFs comparison";
//    maxYrange = 0.1;
    maxYrange = 10.;
    sXaxis = "l_{J/#psi} Error";
    setLogScale = true;
  }
  else if (svarName.Contains("ctau") || svarName.Contains("ctauN") || svarName.Contains("ctauTrue") || svarName.Contains("ctauReco"))
  {
    minRange = -10;
    maxRange = 10;
    sTitle = "l_{J/#psi} resolution PDFs comparison";
//    maxYrange = 0.06;
    maxYrange = 1.;
    sXaxis = "#frac{#font[12]{l}_{J/#psi}(-#font[12]{l}_{J/#psi}^{true})}{#font[12]{#sigma}_{J/#psi}}";
    setLogScale = true;
  }

  
  // bin edges
  float ptmin, ptmax, ymin, ymax, centmin, centmax;
  // collision system
  Char_t collSystem[8];
  
  // workDirName: usual tag where to look for files in Output
  TString workDirNamesStr(workDirNames);
  TString workDirName; Int_t from = 0;
  vector<TString> theDirs;
  while (workDirNamesStr.Tokenize(workDirName, from , ","))
  {
    theDirs.push_back(workDirName);
  }
  
  TString DSTagNamesStr(DSTags);
  TString DSTag; Int_t fromDS = 0;
  int cntData = 0; int nData = -1; int nDirsData = 0;
  vector<TString> theDSTags;
  while (DSTagNamesStr.Tokenize(DSTag, fromDS , ","))
  {
    theDSTags.push_back(DSTag);
    if (DSTag.Contains("DATA"))
    {
      nData = cntData;
      nDirsData++;
    }
    cntData++;
  }
  
  int nDirs = theDirs.size();
  int nDSTags = theDSTags.size();
  
  if (nDirs != nDSTags)
  {
    cout << "[ERROR] number of tags is different from number of dirs" << endl;
    return;
  }
  
  bool allData(false);
  if (nDirs == nDirsData) allData = true;
  
  TString pdfTitlessStr(pdfTitles);
  TString pdfTitlesTemp; Int_t fromPDF = 0; int j = 0;
  vector<TString> pdfTitle;
  while (pdfTitlessStr.Tokenize(pdfTitlesTemp, fromPDF , ","))
  {
    pdfTitle.push_back(pdfTitlesTemp);
  }
  
  map<TString, vector<TString>> filesMap;
  int sizeFiles = 0;
  for (int i = 0 ; i < nDirs ; i++)
  {
    vector<TString> theFiles = fileList(theDirs.at(i),"",theDSTags.at(i),"",fitType);
    cout << "[INFO] Number of files in " << (TString)(theDirs.at(i)).Data() <<  " = " << theFiles.size() << endl;
    
    if ( i != 0 && sizeFiles != (int)theFiles.size())
    {
      cout << "[ERROR] number of files in directories is different" << endl;
      return;
    }
    
    sizeFiles = theFiles.size();
    filesMap[theDirs.at(i)] = theFiles;
  }
  
  vector<TString> files = filesMap[theDirs.at(0)];
  int cnt=0;
  for (vector<TString>::const_iterator it=files.begin(); it!=files.end(); it++)
  {
    vector<TString> currentFiles;
    
    TString mainFile = *it;
    currentFiles.push_back(*it);
    mainFile.Remove(0,mainFile.Last('/')+1);
    mainFile.Remove(0,mainFile.Index("_pt"));
    
    // parse the file name to get info
    anabin thebin = binFromFile(*it);
    ptmin = thebin.ptbin().low();
    ptmax = thebin.ptbin().high();
    ymin = thebin.rapbin().low();
    ymax = thebin.rapbin().high();
    centmin = thebin.centbin().low();
    centmax = thebin.centbin().high();
    strcpy(collSystem, (it->Contains("PbPb")) ? "PbPb" : "PP");
    
    if (!strcmp(collSystem,"PP"))
    {
      centmin = 0;
      centmax = 200;
      isPbPb = false;
    }
    else isPbPb = true;
    
    
    for (int i = 1 ; i < nDirs ; i++)
    {
      vector<TString> filesN = filesMap[theDirs.at(i)];
      for (vector<TString>::const_iterator it2=filesN.begin(); it2!=filesN.end(); it2++)
      {
        TString secFile = *it2;
        TString secFile2 = *it2;
        if (secFile.Contains(collSystem))
        {
          secFile.Remove(0,secFile.Last('/'));
          secFile.Remove(0,secFile.Index("_pt"));
          if(!secFile.CompareTo(mainFile))
          {
            currentFiles.push_back(*it2);
            break;
          }
        }
      }
    }
    
    if ((int)currentFiles.size() != nDirs)
    {
      cout << "[ERROR] I got a different number of files (" << currentFiles.size() << ") than the number of directories(" << nDirs << ")" << endl;
      return;
    }
    else
    {
      cout << "[INFO] I'll plot the PDFs from the following files: " << endl;
      for (int k = 0 ; k < nDirs ; k++)
      {
        cout << currentFiles.at(k).Data() << endl;
      }
    }
    
    RooPlot* frame(0x0);
    RooPlot* frame2(0x0);
    
    int dirCnt = 0;
    for (vector<TString>::const_iterator itWs=currentFiles.begin(); itWs!=currentFiles.end(); itWs++)
    {
      TFile *f = TFile::Open(*itWs); RooWorkspace *ws = NULL;
      if (!f)
      {
        cout << "[ERROR] file " << *itWs << " does not exist." << endl;
      }
      else
      {
        ws = (RooWorkspace*) f->Get("workspace");
        if (!ws)
        {
          cout << "[ERROR] workspace not found in " << *it << "." << endl;
        }
      }
      
      if (f && ws)
      {
        // get the model
        if (!ws->pdf(Form("%s_%s",spdfName.Data(),collSystem)))
        {
          cout << "[ERROR] The pdf " << Form("%s_%s",spdfName.Data(),collSystem) << " does not exist in the file " << *itWs << endl;
          return;
        }
        
        if (dirCnt != nData || allData)
        {
          if (dirCnt ==0)
          {
            frame = ws->var(svarName.Data())->frame(Bins(100), Range(minRange, maxRange));
          }
          
          ws->pdf(Form("%s_%s",spdfName.Data(),collSystem))->plotOn(frame,Name(Form("PDF_%d",dirCnt)),LineColor(kBlack+dirCnt), Precision(1e-5));
          //        ws->pdf(Form("%s%s",spdfName.Data(),collSystem))->plotOn(frame,Name("PDF"), Normalization(myws.data(dsNameCut.c_str())->sumEntries(), RooAbsReal::NumEvent),
          //                                             LineColor(kBlack), Precision(1e-5), NormRange("CtauNWindow"));
        }
        else if (!strcmp(svarName.Data(),"ctauNRes"))
        {
          frame2 = ws->var("ctauN")->frame(Bins(100), Range(minRange, maxRange));
          ws->pdf(Form("%s_%s",spdfName.Data(),collSystem))->plotOn(frame2,Name(Form("PDF_%d",dirCnt)),LineColor(kBlack+dirCnt), Precision(1e-5));
        }
      }
      dirCnt++;
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
    frame->GetYaxis()->SetRangeUser(1E-7,maxYrange);
    
    if (frame2)
    {
      frame2->SetTitle("");
      frame2->GetXaxis()->SetTitle(sXaxis.Data());
      frame2->GetXaxis()->CenterTitle(kTRUE);
      frame2->GetXaxis()->SetLabelSize(0.04);
      frame2->GetXaxis()->SetTitleSize(0.035);
      frame2->GetXaxis()->SetTitleOffset(1.3);
      frame2->GetXaxis()->SetTitleFont(42);
      frame2->GetYaxis()->SetTitle("Prob. distribution");
      frame2->GetYaxis()->SetLabelSize(0.04);
      frame2->GetYaxis()->SetTitleSize(0.04);
      frame2->GetYaxis()->SetTitleOffset(1.8);
      frame2->GetYaxis()->SetTitleFont(42);
      frame2->GetYaxis()->SetRangeUser(1E-7,maxYrange);
    }
    
    TCanvas* cFig = new TCanvas(Form("comparison_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", svarName.Data(), spdfName.Data(), collSystem, (ptmin*10.0), (ptmax*10.0), (ymin*10.0), (ymax*10.0), (int)centmin, (int)centmax),sTitle.Data(),800,800);
    cFig->SetLogy(setLogScale);
    frame->Draw();
    if (frame2)frame2->Draw("same");
    
    // Drawing the text in the plot
    TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
    float dy = 0;
    
    t->SetTextSize(0.03);
    t->DrawLatex(0.35, 0.95, sTitle.Data());
    t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",ptmin,ptmax)); dy+=0.045;
    t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",ymin,ymax)); dy+=0.045;
    if (isPbPb) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(centmin/2), (int)(centmax/2))); dy+=0.045;}
    t->DrawLatex(0.21, 0.86-dy, isPbPb ? "PbPb" : "PP"); dy+=0.045;
    
    // Drawing the Legend
    TLegend* leg = new TLegend(0.6303, 0.7338, 0.8308, 0.8695); leg->SetTextSize(0.03);
    for (int i=0; i<dirCnt; i++) {
      if(frame->findObject(Form("PDF_%d",i))) { leg->AddEntry(frame->findObject(Form("PDF_%d",i)),pdfTitle.at(i).Data(),"l"); }
      else if (frame2)
      {
        if (frame2->findObject(Form("PDF_%d",i))) { leg->AddEntry(frame2->findObject(Form("PDF_%d",i)),pdfTitle.at(i).Data(),"l"); }
      }
    }
    leg->Draw("same");
    
    
    // Save the plot in different formats
    TString saveDir(theDirs.at(0));
    saveDir.Remove(0,saveDir.Last('/')+1);
    
    gSystem->mkdir(Form("pdfComparison/%s/%s/root/", svarName.Data(), saveDir.Data()), kTRUE);
    cFig->SaveAs(Form("pdfComparison/%s/%s/root/PDFComparison_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", svarName.Data(), saveDir.Data(), svarName.Data(), spdfName.Data(), collSystem, (ptmin*10.0), (ptmax*10.0), (ymin*10.0), (ymax*10.0), (int)centmin, (int)centmax));
    
    gSystem->mkdir(Form("pdfComparison/%s/%s/pdf/", svarName.Data(), saveDir.Data()), kTRUE);
    cFig->SaveAs(Form("pdfComparison/%s/%s/pdf/PDFComparison_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", svarName.Data(), saveDir.Data(), svarName.Data(), spdfName.Data(), collSystem, (ptmin*10.0), (ptmax*10.0), (ymin*10.0), (ymax*10.0), (int)centmin, (int)centmax));
    
    cout << endl;
  }
}

#endif // #define comparePDFs_C
