#ifndef drawMassPlot_C
#define drawMassPlot_C

#include "Utilities/initClasses.h"
#include "TGaxis.h"

void setMassRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, double dMuonYmin = -1., double dMuonYmax = -1., bool plotPureSMC=false, bool cutSideBand=false);
void printMassParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted);

void drawMassPlot(RooWorkspace& myws,   // Local workspace
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
                  bool cutCtau,         // Apply prompt ctau cuts
                  bool doSimulFit,      // Do simultaneous fit
                  bool plotPureSMC,     // Flag to indicate if we want to fit pure signal MC
                  // Select the drawing options
                  bool setLogScale,     // Draw plot with log scale
                  bool incSS,           // Include Same Sign data
                  bool zoomPsi,         // Zoom Psi(2S) peak on extra pad
                  int  nBins,           // Number of bins used for plotting
                  bool getMeanPT=false, // Compute the mean PT (NEED TO FIX)
                  bool paperStyle=false // if true, print less info
                  ) 
{

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));
  
  bool isMC = false;
  bool isNPrompt = false;
  if (DSTAG.find("MC")!=std::string::npos)
  {
    isMC = true;
    if(DSTAG.find("NOPR")!=std::string::npos) isNPrompt = true;
  }

  bool applyWeight_Corr = false;
  if ( (plotLabel.find("AccEff")!=std::string::npos) || (plotLabel.find("_lJpsiEff")!=std::string::npos) ) applyWeight_Corr = true;
  else applyWeight_Corr = false;
  bool SB = (incBkg&&(!incPsi2S&&!incJpsi));
  
  string dsOSName = Form("dOS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  string dsSSName = Form("dSS_%s_%s", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
  if(applyWeight_Corr)
  {
    TString corrName = "";
    if (plotLabel.find("AccEff")!=std::string::npos) corrName = "AccEff";
    else if (plotLabel.find("_lJpsiEff")!=std::string::npos) corrName = "lJpsiEff";
    dsOSName = Form("dOS_%s_%s_%s", DSTAG.c_str(),(isPbPb?"PbPb":"PP"),corrName.Data());
  }
//  if(applyWeight_Corr) dsSSName = Form("dSS_%s", DSTAG.c_str());

  string pdfName  = Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP"));
  if (plotPureSMC) dsOSName = Form("dOS_%s_%s_NoBkg", DSTAG.c_str(), (isPbPb?"PbPb":"PP"));
    
  bool isWeighted = myws.data(dsOSName.c_str())->isWeighted();
  string cutSB = parIni["BkgMassRange_FULL_Cut"];
  string cutSBLabel = parIni["BkgMassRange_FULL_Label"];
 
  // Create the main plot of the fit
  RooPlot*   frame     = myws.var("invMass")->frame(Bins(nBins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));

  Double_t norm = myws.data(dsOSName.c_str())->sumEntries();
  if (plotPureSMC)
  {
    if (incJpsi && !incPsi2S)
    {
      if ( (cut.dMuon.AbsRap.Min >= 1.6) || (cut.dMuon.AbsRap.Max > 1.6) ) norm = myws.data(dsOSName.c_str())->reduce("invMass<3.32")->sumEntries();
      else norm = myws.data(dsOSName.c_str())->reduce("invMass<3.26")->sumEntries();
    }
    else if (!incJpsi && incPsi2S)
    {
      if (cut.dMuon.AbsRap.Min >= 1.6) norm = myws.data(dsOSName.c_str())->reduce("invMass<3.95")->sumEntries();
      else norm = myws.data(dsOSName.c_str())->reduce("invMass<3.85")->sumEntries();
    }
  }

  if (paperStyle) TGaxis::SetMaxDigits(3); // to display powers of 10
    
  if (!paperStyle) {
     if (incJpsi) {
        if (incBkg) {
           if ( myws.pdf(Form("pdfMASS2_Jpsi_%s", (isPbPb?"PbPb":"PP"))) ) {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI1"), Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Jpsi_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Jpsi_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                    );
           } else {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Jpsi_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
           }
        } else {
           if ( myws.pdf(Form("pdfMASS2_Jpsi_%s", (isPbPb?"PbPb":"PP"))) ) {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI1"),Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Jpsi_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Jpsi_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                    );
           } else {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("JPSI1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Jpsi_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
           }
        }
     }
     if (incPsi2S) {
        if (incBkg) {
           if ( myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                    LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                    );
           } else {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
           }
        } else {
           if ( myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kOrange+2), LineStyle(1), Precision(1e-4)
                    );
           } else {
              myws.pdf(pdfName.c_str())->plotOn(frame,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                    Normalization(norm, RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"),
                    LineColor(kGreen+3), LineStyle(1), Precision(1e-4)
                    );
           }
        }
     }
  }
  if (incBkg && (!incJpsi && !incPsi2S)) {
    if (!paperStyle) {
       ((RooDataSet*)myws.data(dsOSName.c_str())->reduce(cutSB.c_str()))->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
       myws.pdf(pdfName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                      Normalization(myws.data(dsOSName.c_str())->reduce(cutSB.c_str())->sumEntries(), RooAbsReal::NumEvent), NormRange(cutSBLabel.c_str()),
                                      FillStyle(paperStyle ? 0 : 1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed)
                                      );
    } else {
       myws.pdf(pdfName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                      Normalization(myws.data(dsOSName.c_str())->sumEntries(), RooAbsReal::NumEvent), NormRange(cutSBLabel.c_str()),
                                      FillStyle(paperStyle ? 0 : 1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed)
                                      );
    }
  } 
  if (incBkg && (incJpsi || incPsi2S)) {
    myws.pdf(pdfName.c_str())->plotOn(frame,Name("BKG"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"),
                                      FillStyle(paperStyle ? 0 : 1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed)
                                      );
  } 
  if (incSS) { 
    myws.data(dsSSName.c_str())->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); 
  }
  myws.data(dsOSName.c_str())->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  if (incBkg && (!incJpsi && !incPsi2S)) {
    if (!paperStyle) {
       ((RooDataSet*)myws.data(dsOSName.c_str())->reduce(cutSB.c_str()))->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
       myws.pdf(pdfName.c_str())->plotOn(frame,Name("PDF"),  Normalization(myws.data(dsOSName.c_str())->reduce(cutSB.c_str())->sumEntries(), RooAbsReal::NumEvent), 
             LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange(cutSBLabel.c_str()));
    } else {
       myws.pdf(pdfName.c_str())->plotOn(frame,Name("PDF"),  Normalization(myws.data(dsOSName.c_str())->sumEntries(), RooAbsReal::NumEvent), 
                                      LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange(cutSBLabel.c_str()));
    }
  } else {
    myws.pdf(pdfName.c_str())->plotOn(frame,Name("PDF"),  Normalization(norm, RooAbsReal::NumEvent), 
                                      LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"));
  }
  
  // Create the pull distribution of the fit 
  RooPlot* frameTMP = (RooPlot*)frame->Clone("TMP");
//  int nBinsTMP = nBins;
//  if (plotPureSMC) {
//    double max = cut.dMuon.M.Max, min = cut.dMuon.M.Min;
//    if (incJpsi && !incPsi2S) {
//      if ( (cut.dMuon.AbsRap.Min >= 1.6) || (cut.dMuon.AbsRap.Max > 1.6) ) max = 3.32;
//      else  max = 3.26;
//    }
//    if (!incJpsi && incPsi2S) {
//      if (cut.dMuon.AbsRap.Min >= 1.6)  max = 3.95;
//      else  max = 3.85;
//    }
//    nBinsTMP = (int)(floor( ((max-min)/(cut.dMuon.M.Max-cut.dMuon.M.Min))*nBins ));
//    max = cut.dMuon.M.Min + ((cut.dMuon.M.Max-cut.dMuon.M.Min)/((double)nBins))*((double)nBinsTMP);
//    frameTMP = myws.var("invMass")->frame(Bins(nBinsTMP), Range(min, max));
//    myws.data(dsOSName.c_str())->plotOn(frameTMP, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
//    myws.pdf(pdfName.c_str())->plotOn(frameTMP,Name("PDF"),  Normalization(norm, RooAbsReal::NumEvent), 
//                                      LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), NormRange("MassWindow"));
//  }
  RooHist *hpull = frameTMP->pullHist(0, 0, true);
  hpull->SetName("hpull");
  RooPlot* frame2 = myws.var("invMass")->frame(Title("Pull Distribution"), Bins(nBins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  frame2->addPlotable(hpull, "PX"); 

  // Create the extra PAD for the Psi(2S) zoom
  RooPlot*   framezoom = NULL;
  if(zoomPsi) {  
    setLogScale=false;
    framezoom = myws.var("invMass")->frame(Bins(22), Range(Mass.Psi2S-0.275, Mass.Psi2S+0.275));
    myws.data(dsOSName.c_str())->plotOn(framezoom, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2)); 
    if (!paperStyle) {
       if (incPsi2S) {
          if (incBkg) {
             if ( myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) {
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kGreen+3), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kOrange+2), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
             } else {
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Psi2S_%s", (isPbPb?"PbPb":"PP"))),*myws.pdf(Form("pdfMASSTot_Bkg_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kGreen+3), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
             }
          } else {
             if ( myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ) {
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASS1_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kGreen+3), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S2"),Components(RooArgSet(*myws.pdf(Form("pdfMASS2_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kOrange+2), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
             } else {
                myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PSI2S1"),Components(RooArgSet(*myws.pdf(Form("pdfMASSTot_Psi2S_%s", (isPbPb?"PbPb":"PP"))))),
                      Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                      LineColor(kGreen+3), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)
                      );
             }
          }
       }
    }
    if (incBkg) {
      myws.pdf(pdfName.c_str())->plotOn(framezoom, Name("BKG"),Components(*myws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP")))),
                                        Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"), 
                                        FillStyle(paperStyle ? 0 : 1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed));
    }
    myws.data(dsOSName.c_str())->plotOn(framezoom, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
    myws.pdf(pdfName.c_str())->plotOn(framezoom,Name("PDF"), Normalization(norm, RooAbsReal::NumEvent), NormRange("MassWindow"),
                                      LineColor(kBlack), LineStyle(1),Precision(1e-4));
  }			
  
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
  setMassRange(myws, frame, dsOSName, setLogScale, cut.dMuon.AbsRap.Min, cut.dMuon.AbsRap.Max, plotPureSMC, SB);
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

  if (!zoomPsi) { printMassParameters(myws, pad1, isPbPb, pdfName, isWeighted); }
  pad1->SetLogy(setLogScale);

  // Drawing the text in the plot
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.032);
  float dy = 0; 
  
  t->SetTextSize(0.03);
  if (!paperStyle) { // do not print selection details for paper style
     t->DrawLatex(0.20, 0.86-dy, "2015 HI Soft Muon ID"); dy+=0.045;
     if (cutCtau) { t->DrawLatex(0.21, 0.86-dy, "#font[12]{l}_{J/#psi} cuts applied"); dy+=0.045; }
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
  if (getMeanPT){
    if (incJpsi) {
      t->DrawLatex(0.20, 0.86-dy, Form("<pt_{J/#psi}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptJpsi%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptJpsi%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
    }
    if (incPsi2S) {
      t->DrawLatex(0.20, 0.86-dy, Form("<pt_{#psi(2S)}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptPsi2S%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptPsi2S%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
    }
    if (incBkg) {
      t->DrawLatex(0.20, 0.86-dy, Form("<pt_{bkg}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptBkg%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptBkg%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
    }
  }

  // Drawing the Legend
  double ymin = 0.7802;
  if (incPsi2S && incJpsi && incSS)  { ymin = 0.7202; } 
  if (incPsi2S && incJpsi && !incSS) { ymin = 0.7452; }
  if (paperStyle) { ymin = 0.72; }
  TLegend* leg = new TLegend(0.5175, ymin, 0.7180, 0.8809); leg->SetTextSize(0.03);
  const char* dataName = isMC ? "MC data" : "Data";
  if (frame->findObject("dOS")) { leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":dataName),"pe"); }
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  if (frame->findObject("PDF")) { leg->AddEntry(frame->findObject("PDF"),"Total fit","l"); }
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
  if (!paperStyle) gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 

  // Draw the zoom framed
  if (zoomPsi) {
    framezoom->SetName("zoom_frame_PbPb");
    framezoom->SetTitle("");
    // framezoom->GetYaxis()->SetTitle(frame->GetYaxis()->GetTitle());
    framezoom->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    // if (!paperStyle) {
       framezoom->GetXaxis()->SetTitleOffset(1.0);
       framezoom->GetYaxis()->SetTitleOffset(1.4);
       framezoom->GetXaxis()->SetLabelOffset(0.012);
       framezoom->GetYaxis()->SetLabelSize(0.06);
       framezoom->GetXaxis()->SetLabelSize(0.06);
       framezoom->GetYaxis()->SetTitleSize(0.072);
       framezoom->GetXaxis()->SetTitleSize(0.072);
    // }
       setMassRange(myws, framezoom, dsOSName, setLogScale, cut.dMuon.AbsRap.Min, cut.dMuon.AbsRap.Max, plotPureSMC);

    pad4->Draw();
    pad4->cd();
    
    framezoom->Draw();
    
    cFig->cd(); 
  }

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
     printChi2(myws, pad2, frameTMP, "invMass", dsOSName.c_str(), pdfName.c_str(), nBins/*nBinsTMP*/);

     pline->Draw("same");
     pad2->Update();
  }

  // Save the plot in different formats
  gSystem->mkdir(Form("%smass%s/%s/plot/root/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE); 
  cFig->SaveAs(Form("%smass%s/%s/plot/root/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%smass%s/%s/plot/png/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%smass%s/%s/plot/png/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.png", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  gSystem->mkdir(Form("%smass%s/%s/plot/pdf/", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str()), kTRUE);
  cFig->SaveAs(Form("%smass%s/%s/plot/pdf/PLOT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.pdf", outputDir.c_str(), (SB?"SB":""), DSTAG.c_str(), "MASS", DSTAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  
  cFig->Clear();
  cFig->Close();
};


void setMassRange(RooWorkspace& myws, RooPlot* frame, string dsName, bool setLogScale, double dMuonYmin, double dMuonYmax, bool plotPureSMC, bool cutSideBand)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1* h = myws.data(dsName.c_str())->createHistogram("hist", *myws.var("invMass"), Binning(frame->GetNbinsX(),frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  // Double_t YMin = min( h->GetBinContent(h->FindFirstBinAbove(0.0)), h->GetBinContent(h->FindLastBinAbove(0.0)) );
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  
  bool isMC = false;
  if (dsName.find("MC")!=std::string::npos) isMC = true;
    
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    if (cutSideBand) {
      Yup = YMax*1000.0;
      Ydown = YMin*0.3;
    } else {
      if (isMC) Ydown = YMin*0.3;
      else Ydown = YMin/(TMath::Power((YMax/YMin), (0.1/(1.0-0.1-0.4))));
      Yup = YMax*TMath::Power((YMax/YMin), (0.4/(1.0-0.1-0.4)));
    }
  }
  else
  {
    Ydown = max(YMin-(YMax-YMin)*(0.1/(1.0-0.1-0.4)),0.0);
    Yup = YMax+(YMax-YMin)*(0.4/(1.0-0.1-0.4));
  }
  frame->GetYaxis()->SetRangeUser(Ydown,Yup);
  delete h;
  
  // Create line to indicate upper fitting range for MC
  if (plotPureSMC)
  {
    if (dsName.find("JPSI")!=std::string::npos)
    {
      TLine* lineHigh(0x0);
      if ( (dMuonYmin >= 1.6) || (dMuonYmax > 1.6) ) lineHigh = new TLine(3.32,Ydown,3.32,Yup);
      else lineHigh = new TLine(3.26,Ydown,3.26,Yup);
      lineHigh->SetLineStyle(2);
      lineHigh->SetLineColor(1);
      lineHigh->SetLineWidth(3);
      
//      TLine* lineLow(0x0);
//      lineLow = new TLine(2.2,Ydown,2.2,Yup);
//      lineLow->SetLineStyle(2);
//      lineLow->SetLineColor(1);
//      lineLow->SetLineWidth(3);
      
      frame->addObject(lineHigh);
//      frame->addObject(lineLow);
    }
    else if (dsName.find("PSI2S")!=std::string::npos)
    {
      TLine* line(0x0);
      if (dMuonYmin >= 1.6) line = new TLine(3.95,Ydown,3.95,Yup);
      else line = new TLine(3.85,Ydown,3.85,Yup);
      line->SetLineStyle(2);
      line->SetLineColor(1);
      line->SetLineWidth(3);
      
      frame->addObject(line);
    }

  }
 
};


void printMassParameters(RooWorkspace myws, TPad* Pad, bool isPbPb, string pdfName, bool isWeighted)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.025; 
  RooArgSet* Parameters =  myws.pdf(pdfName.c_str())->getParameters(*myws.var("invMass"));
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    // Parse the parameter's labels
    if(s1=="invMass"){continue;} else if(s1=="MassRatio"){continue;} 
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


#endif // #ifndef drawMassPlot_C
