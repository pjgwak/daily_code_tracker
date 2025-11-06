#ifndef drawMassFrom2DPlot_C
#define drawMassFrom2DPlot_C

#include "Utilities/initClasses.h"
#include "TGaxis.h"

void setMassFrom2DRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale);
void printMassFrom2DParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted);

void drawMassFrom2DPlot(RooWorkspace& myws,   // Local workspace
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
                  bool incBkg,          // Includes Background model                  
                  // Select the fitting options
                  // Select the drawing options
                  bool setLogScale,     // Draw plot with log scale
                  bool incSS,           // Include Same Sign data
                  double  binWidth,     // Bin width
                  bool paperStyle=false // if true, print less info
                  ) 
{

  RooMsgService::instance().getStream(0).removeTopic(Caching);  
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  
  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));
  int nBins = min(int( round((cut.dMuon.M.Max - cut.dMuon.M.Min)/binWidth) ), 1000);
  
  string pdfTotName  = Form("pdfCTAUMASS_Tot_%s", (isPbPb?"PbPb":"PP"));
  string pdfJpsiPRName  = Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP"));
  string pdfJpsiNoPRName  = Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"));
  string pdfPsi2SPRName  = Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP"));
  string pdfPsi2SNoPRName  = Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"));
  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string dsOSNameCut = dsOSName+"_CTAUCUT";
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));

  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();
  bool isMC = (DSTAG.find("MC")!=std::string::npos);

  double normDSTot   = 1.0;  if (myws.data(dsOSNameCut.c_str()))  { normDSTot   = myws.data(dsOSName.c_str())->sumEntries()/myws.data(dsOSNameCut.c_str())->sumEntries();  }
  
  // Create the main plot of the fit
  RooPlot*   frame     = myws.var("invMass")->frame(Bins(nBins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));

 
  if (paperStyle) TGaxis::SetMaxDigits(3); // to display powers of 10
 
  myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                       FillStyle(paperStyle ? 0 : 1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed)
                                       );
  if (!paperStyle) {
    if (incJpsi) {
      if ( myws.pdf(Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP"))) ) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSIPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP"))), *myws.pdf(Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSName.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kRed+3), LineStyle(1), Precision(1e-4), NumCPU(32)
                                             );
      }
      if ( myws.pdf(Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))) ) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("JPSINOPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))), *myws.pdf(Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSName.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4), NumCPU(32)
                                             );
      }
    }
    if (incPsi2S) {
      if ( myws.pdf(Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))) ) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSName.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kRed+3), LineStyle(1), Precision(1e-4), NumCPU(32)
                                             );
      }
      if ( myws.pdf(Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))) ) {
        myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PSI2SNOPR"),Components(RooArgSet(*myws.pdf(Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))))),
                                             ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSName.c_str()), kTRUE),
                                             Normalization(normDSTot, RooAbsReal::NumEvent),
                                             LineColor(kGreen+3), LineStyle(1), Precision(1e-4), NumCPU(32)
                                             );
      }      
    } 
  }
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(pdfTotName.c_str())->plotOn(frame,Name("PDF"),
                                       ProjWData(RooArgSet(*myws.var("ctauErr")), *myws.data(dsOSName.c_str()), kTRUE),
                                       Normalization(normDSTot, RooAbsReal::NumEvent),
                                       LineColor(kBlack), NumCPU(32)
                                       );
  
  // Create the pull distribution of the fit 
  RooPlot* frameTMP = (RooPlot*)frame->Clone("TMP");
  int nBinsTMP = nBins;
  RooHist *hpull = frameTMP->pullHist(0, 0, true);
  hpull->SetName("hpull");
  RooPlot* frame2 = myws.var("invMass")->frame(Title("Pull Distribution"), Bins(nBins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  frame2->addPlotable(hpull, "PX"); 
  
  // set the CMS style
  setTDRStyle();
  
  // Create the main canvas
  TCanvas *cFig  = new TCanvas(Form("cMassFig_%s", (isPbPb?"PbPb":"PP")), "cMassFig",800,800);
  TPad    *pad1  = new TPad(Form("pad1_%s", (isPbPb?"PbPb":"PP")),"",0,paperStyle ? 0 : 0.23,1,1);
  TPad    *pad2  = new TPad(Form("pad2_%s", (isPbPb?"PbPb":"PP")),"",0,0,1,.228);
  TLine   *pline = new TLine(cut.dMuon.M.Min, 0.0, cut.dMuon.M.Max, 0.0);
  
  // TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,paperStyle ? 0.29 : 0.36,0.97,paperStyle ? 0.70 : 0.77);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  frame->SetTitle("");
  frame->GetXaxis()->CenterTitle(kTRUE);
  if (!paperStyle) {
     frame->GetXaxis()->SetTitle("");
     frame->GetXaxis()->SetTitleSize(0.045);
     frame->GetXaxis()->SetTitleFont(42);
     frame->GetXaxis()->SetTitleOffset(3);
     frame->GetXaxis()->SetLabelOffset(3);
     frame->GetYaxis()->SetLabelSize(0.04);
     frame->GetYaxis()->SetTitleSize(0.04);
     frame->GetYaxis()->SetTitleOffset(1.7);
     frame->GetYaxis()->SetTitleFont(42);
  } else {
     frame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
     frame->GetXaxis()->SetTitleOffset(1.1);
     frame->GetYaxis()->SetTitleOffset(1.45);
     frame->GetXaxis()->SetTitleSize(0.05);
     frame->GetYaxis()->SetTitleSize(0.05);
  }
  setMassFrom2DRange(myws, frame, dsOSName, setLogScale);
  if (paperStyle) {
     double Ydown = 0.;//frame->GetMinimum();
     double Yup = 0.9*frame->GetMaximum();
     frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  }
 
  cFig->cd();
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000); 
  if (!paperStyle) pad1->SetBottomMargin(0.015); 
  //plot fit
  pad1->Draw();
  pad1->cd(); 
  frame->Draw();

  printMassFrom2DParameters(myws, pad1, isPbPb, pdfTotName, isWeighted);
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  if (!paperStyle) { // do not print selection details for paper style
     t->DrawLatex(0.20, 0.86-dy, "2015 HI Soft Muon ID"); dy+=0.045;
     if (isPbPb) {
        t->DrawLatex(0.20, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=2.0*0.045;
     } else {
        t->DrawLatex(0.20, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=2.0*0.045;
     } 
  }
  if (cut.dMuon.AbsRap.Min>0.1) {t->DrawLatex(0.5175, 0.86-dy, Form("%.1f < |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=0.045;}
  else {t->DrawLatex(0.5175, 0.86-dy, Form("|y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Max)); dy+=0.045;}
  t->DrawLatex(0.5175, 0.86-dy, Form("%g < p_{T}^{#mu#mu} < %g GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
  if (isPbPb) {t->DrawLatex(0.5175, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}

  // Drawing the Legend
  double ymin = 0.7602;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  if (paperStyle) { ymin = 0.72; }
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  if (frame->findObject("dOS")) { leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe"); }
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if (frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),"Total fit","l"); }
  if (frame->findObject("JPSIPR")) { leg->AddEntry(frame->findObject("JPSIPR"),"Prompt J/#psi","l"); }
  if (frame->findObject("JPSINOPR")) { leg->AddEntry(frame->findObject("JPSINOPR"),"Non-Prompt J/#psi","l"); }
  if (incBkg && frame->findObject("BKG")) { leg->AddEntry(frame->findObject("BKG"),"Background",paperStyle ? "l" : "fl"); }
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
  
  // CMS_lumi(pad1, isPbPb ? 105 : 104, 33, label);
  CMS_lumi(pad1, isPbPb ? 108 : 107, 33, "");
  if (!paperStyle) gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  if (!paperStyle) {
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
     frame2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
     frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

     frame2->Draw(); 

     // *** Print chi2/ndof 
     printChi2(myws, pad2, frameTMP, "invMass", dsOSName.c_str(), pdfTotName.c_str(), nBinsTMP, false);

     pline->Draw("same");
     pad2->Update();
  }

  // Save the plot in different formats
  gSystem->mkdir(Form("%sctauMass/%s/plot/root/", outputDir.c_str(), DSTAG.c_str()), kTRUE); 
  cFig->SaveAs(Form("%sctauMass/%s/plot/root/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctauMass/%s/plot/png/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctauMass/%s/plot/png/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%sctauMass/%s/plot/pdf/", outputDir.c_str(), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%sctauMass/%s/plot/pdf/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  
  cFig->Clear();
  cFig->Close();
};


void setMassFrom2DRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("invMass"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
    Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
  }
  else
  {
    Ydown = max(YMin-(YMax-YMin)*(0.1/(1.0-0.1-0.4)),0.0);
    Yup = YMax+(YMax-YMin)*(0.4/(1.0-0.1-0.4));
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;
 
};


void printMassFrom2DParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.025; 
  RooArgSet* Parameters = (RooArgSet*)myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("invMass"), *myws.var("ctau"), *myws.var("ctauErr")))->selectByAttrib("Constant",kFALSE);
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="invMass" || s1=="ctauErr" || s1=="ctau"){continue;} else if(s1=="MassRatio"){continue;} 
    else if(s1=="One"){continue;} else if(s1=="mMin"){continue;} else if(s1=="mMax"){continue;}
    if(s1=="RFrac2Svs1S"){ s1="R_{#psi(2S)/J/#psi}"; } 
    else if(s1=="rSigma21"){ s1="(#sigma_{2}/#sigma_{1})"; } 
    else if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){
      s1=Form("#%s",s1.c_str());
    }
    if(s2=="PbPbvsPP")   { s2="PbPb/PP";  }
    else if(s2=="Jpsi")  { s2="J/#psi";   } 
    else if(s2=="Psi2S") { s2="#psi(2S)"; }
    else if(s2=="Bkg")   { s2="bkg";      }
    else if(s2=="CtauRes")  { continue; }
    else if(s2=="JpsiNoPR") { continue; }
    else if(s2=="JpsiPR")   { continue; }
    else if(s2=="Psi2SNoPR"){ continue; }
    else if(s2=="Psi2SPR")  { continue; }
    else if(s2=="BkgNoPR")  { continue; }
    else if(s2=="BkgPR")    { continue; }
    else if(s2=="Bkg" && (s1=="N" || s1=="b")) { continue; }
    else {continue;}
    if(s3!=""){
      label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());
    } 
    else {
      label=Form("%s^{%s}", s1.c_str(), s2.c_str());
    }
    // Print the parameter's results
    if(s1=="N"){ 
      t->DrawLatex(0.20, 0.76-dy, Form((isWeighted?"%s = %.6f#pm%.6f ":"%s = %.0f#pm%.0f "), label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("#sigma_{2}/#sigma_{1}")!=std::string::npos){ 
      t->DrawLatex(0.20, 0.76-dy, Form("%s = %.3f#pm%.3f ", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("sigma")!=std::string::npos){ 
      t->DrawLatex(0.20, 0.76-dy, Form("%s = %.2f#pm%.2f MeV/c^{2}", label.c_str(), it->getValV()*1000., it->getError()*1000.)); dy+=0.045; 
    }
    else if(s1.find("lambda")!=std::string::npos){ 
      t->DrawLatex(0.20, 0.76-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else if(s1.find("m")!=std::string::npos){ 
      t->DrawLatex(0.20, 0.76-dy, Form("%s = %.5f#pm%.5f GeV/c^{2}", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
    else { 
      t->DrawLatex(0.20, 0.76-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; 
    }
  }
};


#endif // #ifndef drawMassFrom2DPlot_C
