#ifndef drawCtauTruePlot_C
#define drawCtauTruePlot_C

#include "Utilities/initClasses.h"

void setCtauTrueRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale);
void printCtauTrueParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted);

void drawCtauTruePlot(RooWorkspace& myws,   // Local workspace
                      string outputDir,     // Output directory
                      struct InputOpt opt,  // Variable with run information (kept for legacy purpose)
                      struct KinCuts cut,   // Variable with current kinematic cuts
                      map<string, string>  parIni,   // Variable containing all initial parameters
                      string plotLabel,     // The label used to define the output file name
                      // Select the type of datasets to fit
                      string DSTAG,         // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                      bool isPbPb,          // Define if it is PbPb (True) or PP (False)
                      // Select the type of object to fit
                      bool incJpsi,         // Includes Jpsi model
                      bool incPsi2S,        // Includes Psi(2S) model
                      bool plotPureSMC,     // Flag to indicate if we want to fit pure signal MC
                      // Select the drawing options
                      bool setLogScale,     // Draw plot with log scale
                      bool incSS,           // Include Same Sign data
                      int  nBins            // Number of bins used for plotting
                      ) 
{

  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  bool isMC = false;
  bool isNPrompt = false;
  if (DSTAG.find("MC")!=std::string::npos)
  {
    isMC = true;
    if(DSTAG.find("NOPR")!=std::string::npos) isNPrompt = true;
  }

  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string pdfName  = "pdfCTAUTRUE";
  string pdfTotName  = Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP"));
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
    
  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();

  // Create the main plot of the fit
  RooPlot*   frame = myws.var("ctauTrue")->frame(Bins(nBins), Range(cut.dMuon.ctauTrue.Min, cut.dMuon.ctauTrue.Max));

  if (incSS) {
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(pdfTotName.c_str())->plotOn(frame, Name("PDF"),  Normalization(myws.data(dsOSName.c_str())->sumEntries(), RooAbsReal::NumEvent),
                                       LineColor(kRed+2), LineStyle(1), Precision(1e-4), NormRange("CtauTrueWindow"), Range("CtauTrueWindow")
                                       );
  
  // Create the pull distribution of the fit
  RooHist *hpull = frame->pullHist(0, 0, true);
  hpull->SetName("hpull");
  RooPlot* frame2 = myws.var("ctauTrue")->frame(Title("Pull Distribution"), Range(cut.dMuon.ctauTrue.Min, cut.dMuon.ctauTrue.Max));
  frame2->addPlotable(hpull, "PX");	

  // set the CMS style
  setTDRStyle();
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cCtauTrueFig_%s", (isPbPb?"PbPb":"PP")), "cCtauTrueFig",800,800);
  TPad    *pad1  = new TPad(Form("pad1_%s", (isPbPb?"PbPb":"PP")),"",0,0.23,1,1);
  TPad    *pad2  = new TPad(Form("pad2_%s", (isPbPb?"PbPb":"PP")),"",0,0,1,.228);
  TLine   *pline = new TLine(cut.dMuon.ctauTrue.Min, 0.0, cut.dMuon.ctauTrue.Max, 0.0);

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleFont(42);
  frame->GetXaxis()->SetTitleOffset(3);
  frame->GetXaxis()->SetLabelOffset(3);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.7);
  frame->GetYaxis()->SetTitleFont(42);
  setCtauTrueRange(myws, frame, dsOSName, setLogScale);
 
  cFig->cd();
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000); 
  pad1->SetBottomMargin(0.015); 
  //plot fit
  pad1->Draw();
  pad1->cd(); 
  frame->Draw();

  printCtauTrueParameters(myws, pad1, isPbPb, pdfTotName, isWeighted);
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, "2015 HI Soft Muon ID"); dy+=0.045;
  if (isPbPb) {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } else {
    t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
  } 
  if (isPbPb) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  t->DrawLatex(0.21, 0.86-dy, Form("%.1f #leq |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=1.5*0.045;

  // Drawing the Legend
  double ymin = 0.7802;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  const char* dataName = isMC ? "MC data" : "Data";
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  if (frame->findObject("dOS")) { leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":dataName),"pe"); }
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if(frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),"Total fit","l"); }
  leg->Draw("same");

  //Drawing the title
  TString label;
  if (isPbPb) {
    if (opt.PbPb.RunNb.Start==opt.PbPb.RunNb.End){
      label = Form("PbPb Run %d", opt.PbPb.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PbPb", "HIOniaL1DoubleMu0", opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End);
    }
  } else {
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  }
  
  //CMS_lumi(pad1, isPbPb ? 105 : 104, 33, label);
  int fc = isMC ? -1 : 1;
  TString lumiLabel("");
  if (isMC)
  {
    if (isNPrompt) lumiLabel += "nonprompt";
    else lumiLabel += "prompt";
    
    if (incJpsi) lumiLabel += " J/#psi";
    else lumiLabel += " #psi(2S)";
  }
  CMS_lumi(pad1, isPbPb ? fc*108 : fc*107, 33, lumiLabel.Data());
  //CMS_lumi(pad1, isPbPb ? 108 : 107, 33, "");
  gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  //---plot pull
  pad2->Draw();
  pad2->cd();
    
  frame2->SetTitle("");
  frame2->GetYaxis()->CenterTitle(kTRUE);
  frame2->GetYaxis()->SetTitleOffset(0.4);
  frame2->GetYaxis()->SetTitleSize(0.1);
  frame2->GetYaxis()->SetLabelSize(0.1);
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetXaxis()->CenterTitle(kTRUE);
  frame2->GetXaxis()->SetTitleOffset(1);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetXaxis()->SetLabelSize(0.1);
  frame2->GetXaxis()->SetTitle("#font[12]{l}_{J/#psi} MC True (mm)");
  frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

  frame2->Draw(); 
  
  // *** Print chi2/ndof 
  printChi2(myws, pad2, frame, "ctauTrue", dsOSName.c_str(), pdfTotName.c_str(), nBins);
  
  pline->Draw("same");
  pad2->Update();
  
  // Save the plot in different formats
  gSystem->mkdir(Form("%sctauTrue/%s/plot/root/", outputDir.c_str(), DSTAG.c_str()), kTRUE); 
  cFig->SaveAs(Form("%sctauTrue/%s/plot/root/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), DSTAG.c_str(), "CTAUTRUE", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctauTrue/%s/plot/png/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctauTrue/%s/plot/png/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), DSTAG.c_str(), "CTAUTRUE", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctauTrue/%s/plot/pdf/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctauTrue/%s/plot/pdf/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), DSTAG.c_str(), "CTAUTRUE", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));

  cFig->Clear();
  cFig->Close();

}

void setCtauTrueRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("ctauTrue"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) isMC = true;
    
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/0.6)));
    Yup = YMax*TMath::Power((YMax/YMin), (0.3/0.6));
  }
  else
  {
    Ydown = max(YMin-(YMax-YMin)*0.2,0.0);
    Yup = YMax+(YMax-YMin)*0.5;
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;

}


void printCtauTrueParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.025; 
  RooArgSet* Parameters =  myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauTrue")));
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="ctauTrue"){continue;}
    else if(s1=="One"){continue;}
    else if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){
      s1=Form("#%s",s1.c_str());
    }

    if(s2=="CtauRes")  { s2="Res";   } 
    else if(s2=="JpsiNoPR")  { s2="J/#psi";   }
    else if(s2=="Jpsi" && s1=="N")  { s2="J/#psi";  } 
    else if(s2=="Psi2SNoPR") { s2="#psi(2S)"; }
    else if(s2=="Psi2S" && s1=="N")  { s2="#psi(2S)";   }
    else {continue;}
    if(s3!=""){
      label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());
    }
    // Print the parameter's results
    if(s1=="N"){ 
      t->DrawLatex(0.69, 0.70-dy, Form((isWeighted?"%s = %.6f#pm%.6f ":"%s = %.0f#pm%.0f "), label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("sigma")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.70-dy, Form("%s = %.4f#pm%.4f #mum", label.c_str(), it->getValV()*1000., it->getError()*1000.)); dy+=0.045; 
    }
    else if(s1.find("lambda")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.70-dy, Form("%s = %.4f#pm%.4f mm", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("ctau")!=std::string::npos){ 
      t->DrawLatex(0.69, 0.70-dy, Form("%s = %.4f#pm%.4f mm", (label.insert(1, string("#"))).c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else { 
      t->DrawLatex(0.69, 0.70-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
  }
};


#endif // #ifndef drawCtauTruePlot_C
