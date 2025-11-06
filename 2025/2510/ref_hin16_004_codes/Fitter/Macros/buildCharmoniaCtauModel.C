#ifndef buildCharmoniaCtauModel_C
#define buildCharmoniaCtauModel_C

#include "Utilities/initClasses.h"

void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbPb);
void setCtauDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries);
bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb);
bool createCtauBkgTemplate(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth);
bool ctauBkgHistToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, string dsName, vector<double> range, bool useDataSet=true);
TH1* rebinctauBkghist(RooWorkspace& ws, TH1 *hist, double xmin=1e99, double xmax=-1e99);
void findRangeFromPlot(vector<double>& range, struct KinCuts cut, double binWidth);

bool buildCharmoniaCtauModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, string dsName,
                             struct KinCuts cut,          // Variable containing all kinematic cuts
                             bool isPbPb,                 // Determine if we are working with PbPb (True) or PP (False)
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             bool incPsi2S,               // Include Psi(2S) model
                             bool incPrompt,              // Include Prompt models
                             bool incNonPrompt,           // Include NonPrompt models
                             bool useTotctauErrPdf,       // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                             bool usectauBkgTemplate,     // If yes use a template for Bkg ctau instead of the fitted Pdf
                             double  binWidth,            // Bin width
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{
  // If the initial parameters are empty, set defaul parameter values
  setCtauDefaultParameters(parIni, isPbPb, numEntries);
  
  // Fix all psi2S parameters to jpsi
  if (incJpsi && incPsi2S) {
    fixCtauParPsi2StoJpsi(parIni, isPbPb);
  }

  // C r e a t e   m o d e l 

  bool fitMass = false;
  if ( ws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP"))) && incBkg && (incJpsi || incPsi2S) ) { fitMass = true; } 

  string pdfName     = "pdfCTAU";
  if (fitMass) { pdfName = "pdfCTAUMASS"; }
  bool isMC = (dsName.find("MC")!=std::string::npos);
  bool incCtauErrPDF = true;
  
  
  if (incJpsi) {
    if(!defineCtauResolModel(ws, "Jpsi", model.CtauRes, parIni, isPbPb)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
    if (incPrompt) {
      if(!addSignalCtauModel(ws, "JpsiPR", model.Jpsi.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Jpsi Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!addSignalCtauModel(ws, "JpsiNoPR", model.Jpsi.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Jpsi Ctau Model failed" << endl; return false; }
    }
    if (incPrompt && !incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsi(Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Jpsi_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfJpsi);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP")),
                       Form("%s_JpsiPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                       Form("pdfMASS_Jpsi_%s",(isPbPb?"PbPb":"PP"))
                       ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_JpsiPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
    if (!incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsi(Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Jpsi_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfJpsi);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                         Form("%s_JpsiNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                         Form("pdfMASS_Jpsi_%s",(isPbPb?"PbPb":"PP"))
                         ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_JpsiNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_JpsiNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
    if (incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfJpsiPR(Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Jpsi_%s"), (isPbPb?"PbPb":"PP"))),
                             Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                             );
        ws.import(pdfJpsiPR);
        if (!ws.pdf(Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")))) {
          RooProdPdf pdfJpsiNoPR(Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Jpsi_%s"), (isPbPb?"PbPb":"PP"))),
                                 Conditional( *ws.pdf(Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                                 ); 
          ws.import(pdfJpsiNoPR);
        }
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      ws.factory( parIni[Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() );
      if (incCtauErrPDF) {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      } else {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_JpsiPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Jpsi_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_JpsiNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Jpsi_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_JpsiPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Jpsi_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_Jpsi_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
  }
  if (incPsi2S) {
    if(!defineCtauResolModel(ws, "Psi2S", model.CtauRes, parIni, isPbPb)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
    if (incPrompt) {
      if(!addSignalCtauModel(ws, "Psi2SPR", model.Psi2S.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Psi2S Ctau Model failed" << endl; return false; }
    }
    if (incNonPrompt) {
      if(!addSignalCtauModel(ws, "Psi2SNoPR", model.Psi2S.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Psi2S Ctau Model failed" << endl; return false; }
    }
    if (incPrompt && !incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2S(Form("pdfCTAU_Psi2SPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Psi2S_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2S);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ){
         ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP")),
                         Form("%s_Psi2SPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                         Form("pdfMASS_Psi2S_%s",(isPbPb?"PbPb":"PP"))
                         ));
       }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_Psi2SPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
    if (!incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2S(Form("pdfCTAU_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Psi2S_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2S);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_Psi2SNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbPb?"PbPb":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2SNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_Psi2SNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
    if (incPrompt && incNonPrompt) {
      if (incCtauErrPDF) {
        RooProdPdf pdfPsi2SPR(Form("pdfCTAU_Psi2SPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Psi2S_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2SPR);
        RooProdPdf pdfPsi2SNoPR(Form("pdfCTAU_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Psi2S_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       ); 
        ws.import(pdfPsi2SNoPR);
      } else {
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      ws.factory( parIni[Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str() );
      if (incCtauErrPDF) {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      } else {
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUCOND_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if ( ws.pdf(Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP"))) ){
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_Psi2SPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_Psi2SNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Psi2S_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_Psi2SPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Psi2S_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("%s_Psi2S_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))
                      ));
    }
  }  
  if (incBkg) {
    if (!usectauBkgTemplate)
    {
      if(!defineCtauResolModel(ws, "Bkg", model.CtauRes, parIni, isPbPb)) { cout << "[ERROR] Defining the Prompt Ctau Resolution Model failed" << endl; return false; }
      if (incPrompt) {
        if(!addBackgroundCtauModel(ws, "BkgPR", model.Bkg.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Bkg Ctau Model failed" << endl; return false; }
      }
      if (incNonPrompt) {
        if(!addBackgroundCtauModel(ws, "BkgNoPR", model.Bkg.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Bkg Ctau Model failed" << endl; return false; }
      }
      if (incPrompt && !incNonPrompt) {
        if (incCtauErrPDF) {
          RooProdPdf pdf(Form("pdfCTAU_BkgPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                         Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                         );
          ws.import(pdf);
        }
        if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))) ){
          ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgPR_%s", (isPbPb?"PbPb":"PP")),
                          Form("%s_BkgPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                          Form("pdfMASS_Bkg_%s",(isPbPb?"PbPb":"PP"))
                          ));
        }
        ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_BkgPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("%s_BkgPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if (!incPrompt && incNonPrompt) {
        if (incCtauErrPDF) {
          RooProdPdf pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                         Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                         );
          ws.import(pdf);
        }
        if ( ws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))) ){
          ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")),
                          Form("%s_BkgNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                          Form("pdfMASS_Bkg_%s",(isPbPb?"PbPb":"PP"))
                          ));
        }
        ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_BkgNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("%s_BkgNoPR_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))
                        ));
      }
      if (incPrompt && incNonPrompt) {
        if (incCtauErrPDF) {
          if ( !ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbPb?"PbPb":"PP"))) ) { return false; }
          RooProdPdf pdfPR(Form("pdfCTAU_BkgPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                           Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                           );
          ws.import(pdfPR);
          RooProdPdf pdfNoPR(Form("pdfCTAU_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                             Conditional( *ws.pdf(Form("pdfCTAUCOND_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                             );
          ws.import(pdfNoPR);
        }
        ws.factory( parIni[Form("b_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str() );
        if (incCtauErrPDF) {
          ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")),
                          Form("b_Bkg_%s", (isPbPb?"PbPb":"PP")),
                          Form("pdfCTAU_BkgNoPR_%s", (isPbPb?"PbPb":"PP")),
                          Form("pdfCTAU_BkgPR_%s", (isPbPb?"PbPb":"PP"))
                          ));
        } else {
          ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_Bkg_%s", (isPbPb?"PbPb":"PP")),
                          Form("b_Bkg_%s", (isPbPb?"PbPb":"PP")),
                          Form("pdfCTAUCOND_BkgNoPR_%s", (isPbPb?"PbPb":"PP")),
                          Form("pdfCTAUCOND_BkgPR_%s", (isPbPb?"PbPb":"PP"))
                          ));
        }
      }
    }
    else if (!fitMass)
    {
      if (usectauBkgTemplate) {
        if (incPrompt && incNonPrompt ) { if( !createCtauBkgTemplate(ws, dsName, "pdfCTAUCOND", cut, incJpsi, incPsi2S, binWidth)) {cout << "[ERROR] Creating the Bkg Ctau Template failed" << endl; return false;}}
        else {cout << "[ERROR] To create the Bkg Ctau Template we need to activate prompt and nonpromt Jpsi" << endl; return false;}
      }
      if ( !ws.pdf(Form("pdfCTAUERR_Bkg_%s", (isPbPb?"PbPb":"PP"))) ) { return false; }
      RooProdPdf pdf(Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                     Conditional( *ws.pdf(Form("pdfCTAUCOND_Bkg_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       );
      ws.import(pdf);
    }
    if ( fitMass && ws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))) ){
      if (!usectauBkgTemplate) {
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_BkgPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_BkgNoPR_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP")),
                        Form("b_Bkg_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUMASS_BkgPR_%s", (isPbPb?"PbPb":"PP"))
                        ));
      } 
      else {
        RooProdPdf pdf(Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")), "", *ws.pdf(Form((useTotctauErrPdf?"pdfCTAUERRTot_Tot_%s":"pdfCTAUERR_Bkg_%s"), (isPbPb?"PbPb":"PP"))),
                       Conditional( *ws.pdf(Form("pdfCTAUCOND_Bkg_%s", (isPbPb?"PbPb":"PP"))), RooArgList(*ws.var("ctau")) )
                       );
        ws.import(pdf);
        ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP")),
                        Form("%s_Bkg_%s", (incCtauErrPDF ? "pdfCTAU" : "pdfCTAUCOND"), (isPbPb?"PbPb":"PP")),
                        Form("pdfMASS_Bkg_%s",(isPbPb?"PbPb":"PP"))
                        ));
      }
    }
    if ( !ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))) ){ ws.factory( parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str() ); } 
    ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("%sTot_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                    Form("%s_Bkg_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")),
                    Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))
                    ));
  }
  

  // Total PDF
  RooAbsPdf *themodel = NULL;
  string tag = "";
  if (incPrompt  && !incNonPrompt) { tag = "PR"; }
  if (!incPrompt && incNonPrompt)  { tag = "NoPR";   }

  if (incJpsi && incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))) 
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))), 
                                        (ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))),
                                        *ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))) 
                                         )
                             );
  }
  if (incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))), 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))),
                                        (ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))))
                                         )
                             );
  }
  if (incJpsi && !incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))), 
                                        *ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))
                                         )
                             );
  }
  if (!incJpsi && incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP"))), 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        (ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))),
                                        *ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))
                                         )
                             );
  }
  if (incJpsi && !incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Jpsi%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))
                                         )
                             );
  }
  if (!incJpsi && incPsi2S && !incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Psi2S%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        (ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))?*ws.var(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))):*ws.function(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && incBkg) {
    themodel = new RooAddPdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")), 
                             RooArgList( 
                                        *ws.pdf(Form("%s_Bkg%s_%s", pdfName.c_str(), tag.c_str(), (isPbPb?"PbPb":"PP")))
                                         ),
                             RooArgList( 
                                        *ws.var(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))
                                         )
                             );
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }
  ws.import(*themodel);
  //ws.pdf(Form("%s_Tot_%s", pdfName.c_str(), (isPbPb?"PbPb":"PP")))->setNormRange("CtauWindow");
  
  return true;
};


void fixCtauParPsi2StoJpsi(map<string, string>& parIni, bool isPbPb)
{
  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi" << endl;
  parIni[Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("sigmaMC_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), Form("sigmaMC_JpsiNoPR_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("lambdaDSS_Psi2SNoPR_PbPb")] = Form("RooFormulaVar::%s('@0',{%s})", Form("lambdaDSS_Psi2SNoPR_PbPb"), Form("lambdaDSS_JpsiNoPR_PbPb"));
  //parIni[Form("b_Psi2S_PbPb")] = Form("RooFormulaVar::%s('@0',{%s})", Form("b_Psi2S_PbPb"), Form("b_Jpsi_PbPb"));
};


bool defineCtauResolModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb)
{ 
  if (ws.pdf(Form("pdfCTAURES_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[WARNING] The %s Ctau Resolution Model has already been implemented!", object.c_str()) << endl;
    return true;
  }
  if (ws.pdf(Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[WARNING] The %s Ctau Model has already been implemented!", object.c_str()) << endl;
    return true; 
  }

  string pdfType = "pdfCTAURES";
  string varName = "ctau";
  bool usePerEventError = true;

  cout << Form("[INFO] Implementing %s Ctau Resolution Model", object.c_str()) << endl;
  
  switch(model) 
    {
    case (CtauModel::SingleGaussianResolution):  
      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) 
             )) {
 	cout << Form("[ERROR] Initial parameters where not found for Single Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }

      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );  }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")))) { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }

      // create the PDF
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
 		      Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      
      cout << Form("[INFO] %s Single Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;

    case (CtauModel::DoubleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Double Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );      }
      
      // create the two PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
 		      Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));

      cout << Form("[INFO] %s Double Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;


    case (CtauModel::TripleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP")))
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Triple Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );      }
      if (!ws.var(Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))      { ws.factory( parIni[Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );     }
      
      // create the three PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
 		      Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s23_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s23_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));

      cout << Form("[INFO] %s Triple Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;

    case (CtauModel::QuadrupleGaussianResolution):  

      if (!( 
            parIni.count(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))) &&
            parIni.count(Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP")))
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for Triple Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      
      // create the variables for this model  
      if (!ws.var("One")) { ws.factory("One[1.0]"); }
      if (!ws.var(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))))   { ws.factory( parIni[Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()  ); }
      if (!ws.var(Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))))  { ws.factory( parIni[Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() ); }
      if (!ws.var(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))))       { ws.factory( parIni[Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );      }
      if (!ws.var(Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))))      { ws.factory( parIni[Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );     }
      if (!ws.var(Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP"))))      { ws.factory( parIni[Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str() );     }
      
      // create the four PDFs
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
 		      Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, One, %s)", Form("%s4_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), varName.c_str(),
 		      Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP")), 
 		      Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                      (usePerEventError?"ctauErr":"One")
 		      ));

      // combine the two PDFs
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s34_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("%s3_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s4_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s234_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("%s2_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s34_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("%s_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("%s1_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("%s234_%s_%s", pdfType.c_str(), object.c_str(), (isPbPb?"PbPb":"PP")),  
 		      Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))
 		      ));

      cout << Form("[INFO] %s Quadruple Gaussian Ctau Resolution PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
  
    default :

      cout<< "[ERROR] Selected Ctau Resolution Model has not been implemented"<< endl; return false;

    }
 
  return true;

};


bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb) 
{
  if (ws.pdf(Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[WARNING] The %s Background Ctau Model has already been implemented!", object.c_str()) << endl;
    return true; 
  }
  if (ws.pdf(Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[INFO] The %s Background Ctau Conditional Model has already been implemented!", object.c_str()) << endl;
    return true; 
  }

  cout << Form("[INFO] Implementing %s Background Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   }

  switch(model) 
    {  
    case (CtauModel::TripleDecay): 
      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Background Triple Decay Ctau Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      } 
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      ws.factory( parIni[Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()  );
      ws.factory( parIni[Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      ws.factory( parIni[Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()     ); 
      ws.factory( parIni[Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()     );

      // create the three PDFs
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
 		      Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
 		      Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));

      // combine the three PDFs
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
 
      cout << Form("[INFO] %s Background Triple Decay Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
        
        
      case (CtauModel::QuadrupleDecay):
        if (!(
              parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("lambdaDSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("fDFSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
              )) {
          cout << Form("[ERROR] Initial parameters where not found for %s Background Quadruple Decay Ctau Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
        }
        
        // create the variables for this model
        if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))) ){
          ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
        }
        ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
        ws.factory( parIni[Form("lambdaDSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
        ws.factory( parIni[Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()  );
        ws.factory( parIni[Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
        ws.factory( parIni[Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()     );
        ws.factory( parIni[Form("fDFSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()     );
        ws.factory( parIni[Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()     );
        
        // create the four PDFs
        ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau",
                        Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau",
                        Form("lambdaDSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau",
                        Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau",
                        Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        
        // combine the three PDFs
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("fDFSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUDSS2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                        Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                        ));
        
        cout << Form("[INFO] %s Background Quadruple Decay Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
 
    case (CtauModel::Delta): 
      
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      }

      // create the PDF
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));

      cout << Form("[INFO] %s Background Delta Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 
      
    default :

      cout<< "[ERROR] Selected Background Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};


bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb) 
{
  if (ws.pdf(Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[WARNING] The %s Signal Ctau Model has already been implemented!", object.c_str()) << endl;
    return true; 
  }
  if (ws.pdf(Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))) { 
    cout << Form("[INFO] The %s Signal Ctau Conditional Model has already been implemented!", object.c_str()) << endl;
    return true; 
  }

  cout << Form("[INFO] Implementing %s Signal Ctau Model", object.c_str()) << endl;
   
  string objectInc = object; 
  if (objectInc.find("NoPR")!=std::string::npos) { objectInc.erase(objectInc.find("NoPR"), objectInc.length()); }
  if (objectInc.find("PR")!=std::string::npos)   { objectInc.erase(objectInc.find("PR"), objectInc.length());   } 

  switch(model) 
    {  
    case (CtauModel::SingleSidedDecay): 

      if (!( 
            parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
 	cout << Form("[ERROR] Initial parameters where not found for %s Signal Single Sided Decay Ctau Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      
      // create the variables for this model 
      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      }
      ws.factory( parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
            
      // create the PDF
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
 		      Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
 
      cout << Form("[INFO] %s Signal Single Sided Decay Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 
 
    case (CtauModel::Delta):

      if ( !ws.var(Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))) ){ 
        ws.factory( parIni[Form("N_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      }

      // create the three PDFs 
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAUCOND_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
 		      Form("pdfCTAURES_%s_%s", objectInc.c_str(), (isPbPb?"PbPb":"PP"))
 		      ));
 		      
      cout << Form("[INFO] %s Signal Delta Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 
 
    default :

      cout<< "[ERROR] Selected Signal Ctau Model has not been implemented"<< endl; return false;

    }
   
  return true;
};
  

void setCtauDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")), numEntries, 0.0, numEntries*2.0);
  }
  if (parIni.count(Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.2, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0.2, 0.0, 1.0);
  }
  if (parIni.count(Form("b_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("b_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("b_Bkg_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.12f,%.12f,%.12f]", Form("b_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.2, 0.0, 1.0);
  }

 // CTAU FIT PARAMETERS

 // Resolution Ctau Model
  if (parIni.count(Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    
  }
  if (parIni.count(Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")), Form("ctau1_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if (parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP")), Form("ctau2_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if (parIni[Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0',{%s})", Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP")), Form("ctau3_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if (parIni[Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("ctau4_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
  }
  if (parIni.count(Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.0, 0.01, 2.0);
  }
  if (parIni.count(Form("rS21_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("rS21_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("rS21_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS21_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("rS32_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("rS32_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("rS32_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS32_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("rS43_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("rS43_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("rS43_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rS43_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 1.0, 3.0);
  }
  if (parIni.count(Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                                                                   parIni[Form("rS21_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), Form("s1_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if ( parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                                                                   parIni[Form("rS32_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), Form("s2_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if ( parIni[Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0) {
    parIni[Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP")),
                                                                   parIni[Form("rS43_CtauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), Form("s3_CtauRes_%s", (isPbPb?"PbPb":"PP") ));
  } else if ( parIni[Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("s4_CtauRes_%s", (isPbPb?"PbPb":"PP")), 1.2, 0.001, 6.0);
  }
  if (parIni.count(Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.4, 1.0);
  }
  if (parIni.count(Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f2_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.4, 1.0);
  }
  if (parIni.count(Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("f3_CtauRes_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.4, 1.0);
  }

  // Signal Ctau Model
  if (parIni.count(Form("lambdaDSS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDSS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_JpsiNoPR_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.01, 2.0);
  }
  if (parIni.count(Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_Psi2SNoPR_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.01, 2.0);
  }

  // Background Ctau Model
  if (parIni.count(Form("f_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("f_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("fDFSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.8, 0., 1.);
  }
  if (parIni.count(Form("fDFSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("fDFSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("fDFSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.8, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("fDLIV_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("fDLIV_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.5, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.45, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("lambdaDSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDSS2_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.45, 0.0001, 5.0);
  }
  if (parIni.count(Form("lambdaDF_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDF_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("lambdaDF_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDF_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.30, 0.00001, 5.0);
  }
  if (parIni.count(Form("lambdaDDS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgNoPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.12f,%.12f,%.12f]", Form("lambdaDDS_BkgNoPR_%s", (isPbPb?"PbPb":"PP")), 0.06, 0.0001, 5.0);
  }
  if (parIni.count(Form("f_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("f_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.3, 0., 1.);
  }
  if (parIni.count(Form("fDFSS_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("fDFSS_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("fDFSS_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("fDLIV_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("fDLIV_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("fDLIV_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.9, 0., 1.);
  }
  if (parIni.count(Form("lambdaDSS_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDSS_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDSS_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.42, 0.05, 1.5);
  }
  if (parIni.count(Form("lambdaDF_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDF_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDF_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDF_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.79, 0.001, 1.5);
  }
  if (parIni.count(Form("lambdaDDS_BkgPR_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambdaDDS_BkgPR_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("lambdaDDS_BkgPR_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDDS_BkgPR_%s", (isPbPb?"PbPb":"PP")), 0.69, 0.02, 5.0);
  }

};


bool createCtauBkgTemplate(RooWorkspace& ws, string dsName, string pdfType, struct KinCuts cut, bool incJpsi, bool incPsi2S, double binWidth)
{
  string hType = pdfType;
  hType.replace(hType.find("pdf"), string("pdf").length(), "h");
  
  bool isPbPb = false;
  if (dsName.find("PbPb")!=std::string::npos) { isPbPb = true; }
  if (dsName.find("MC")!=std::string::npos)   { return false;  }  // Only accept data
  
  // create weighted data sets
  vector<double> range; range.push_back(0); range.push_back(0);
  findRangeFromPlot(range,cut,binWidth);
  
  double ctauMax = range[1], ctauMin = range[0];
  int nBins = min(int( round((ctauMax - ctauMin)/binWidth) ), 1000);
  
  TH1D* hTot = (TH1D*)ws.data(dsName.c_str())->createHistogram(Form("%s_Bkg_%s", hType.c_str(), (isPbPb?"PbPb":"PP")), *ws.var("ctau"), Binning(nBins, ctauMin, ctauMax));
  RooDataHist* dataHist = new RooDataHist("tmp", "", *ws.var("ctau"), hTot); delete dataHist; // KEEP THIS LINE, IT IS A FIX FOR ROOKEYS
  if ( !ctauBkgHistToPdf(ws, hTot, Form("%s_Bkg_%s", pdfType.c_str(), (isPbPb?"PbPb":"PP")), dsName, range)) { return false; }
  hTot->Delete();

  return true;
};

bool ctauBkgHistToPdf(RooWorkspace& ws, TH1D* hist, string pdfName, string dsName, vector<double> range, bool useDataSet)
{
  if (ws.pdf(pdfName.c_str())) {
    cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << endl;
    return true;
  }
  
  cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << endl;
  
  if (useDataSet) {
    if (ws.data(dsName.c_str())==NULL) { cout << "[ERROR] DataSet " << dsName << " for Bkg Template was not found!" << endl; return false; }
    RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataSet*)ws.data(dsName.c_str())), RooKeysPdf::MirrorAsymBoth);
    pdf->setNormRange("CtauFullWindow");
    ws.import(*pdf);
    delete pdf;
  }
  else {
    // Cleaning the input histogram
    // 1) Remove the Under and Overflow bins
    hist->ClearUnderflowAndOverflow();
    // 2) Set negative bin content to zero
    for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0) { hist->SetBinContent(i, 0.0000000001); } }
    // 3) Reduce the range of histogram and rebin it
    TH1* hClean = rebinctauBkghist(ws, hist, range[0], range[1]);
    // 4) Create the RooDataHist for the Template
    string dataName = pdfName;
    dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
    RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *ws.var("ctau"), hClean);
    if (dataHist==NULL) { cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << endl; return false; }
    if (dataHist->sumEntries()==0) { cout << "[ERROR] DataHist used to create " << pdfName << " is empty!" << endl; return false; }
    if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { cout << "[ERROR] DataHist used to create " << pdfName << "  " << " is invalid!  " << endl; return false; }
    ws.import(*dataHist);
    // 4) Make the Template
    RooHistPdf* pdf = new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataHist*)ws.data(dataName.c_str())));
    //RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctau"), *((RooDataSet*)ws.data(dataName.c_str())),RooKeysPdf::NoMirror, 0.4)
    pdf->setNormRange("CtauFullWindow");
    ws.import(*pdf);
    delete pdf;
    delete dataHist;
  }
  
  return true;
};


TH1* rebinctauBkghist(RooWorkspace& ws, TH1 *hist, double xmin, double xmax)
{
  TH1 *hcopy = (TH1*) hist->Clone("hcopy");
  
  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();
  
  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.1) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.00000001) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          hcopy->SetBinContent(i,hcopy->GetBinContent(i)/nrebin);
          break;
        }
        nrebin++;
      }
    }
  }
  
  if (xmin < newbins[1]) newbins[0] = xmin;
  if (xmax > newbins[newbins.size()-2]) newbins[newbins.size()-1] = xmax;
  
  TH1 *ans = hcopy->Rebin(newbins.size()-1,"hnew",newbins.data());
  
// //Functionality not working yet
//  vector<double> fullBins;
//  fullBins.push_back(-3.0); //hardcoded
//  for (int i=0 ; i < (int)newbins.size() ; i++)
//  {
//    fullBins.push_back(newbins[i]);
//  }
//  fullBins.push_back(5.0);//hardcoded
//  
//  RooBinning* tempBin = new RooBinning(fullBins.size()-1,fullBins.data(),"TemplateBinning");
//  ws.import(*tempBin);
  
  delete hcopy;
  return ans;
};

void findRangeFromPlot(vector<double>& range, struct KinCuts cut, double binWidth)
{
  double minRange = -4.0; //hardcoded (must be the same as in drawCtauPlot.C)
  double maxRange = 7.0;
  
  double ctauMax = cut.dMuon.ctau.Max, ctauMin = cut.dMuon.ctau.Min;
  if (ctauMax > maxRange) range[1] = maxRange;
  if (ctauMin < minRange) range[0] = minRange;
  
  int nBins = min(int( round((maxRange - minRange)/binWidth) ), 1000);
  
  double binLeftBoundLow(minRange), binBoundLeftHigh(minRange+binWidth);
  double binRighttBoundLow(maxRange-binWidth), binBoundRightHigh(maxRange);
  for (int i = 1 ; i <= nBins ; i++)
  {
    if ( (ctauMin >= binLeftBoundLow) && (ctauMin < binBoundLeftHigh)) range[0] = binLeftBoundLow;
    else
    {
      binLeftBoundLow += binWidth;
      binBoundLeftHigh += binWidth;
    }
    
    if ( (ctauMax >= binRighttBoundLow) && (ctauMax < binBoundRightHigh)) range[1] = binBoundRightHigh;
    else
    {
      binRighttBoundLow -= binWidth;
      binBoundRightHigh -= binWidth;
    }
  }
  
  return;
};

#endif // #ifndef buildCharmoniaCtauModel_C
