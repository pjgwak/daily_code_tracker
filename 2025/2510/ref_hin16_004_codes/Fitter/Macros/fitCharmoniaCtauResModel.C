#ifndef fitCharmoniaCtauResModel_C
#define fitCharmoniaCtauResModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauResModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "drawCtauResPlot.C"


void setCtauResCutParameters(struct KinCuts& cut, bool incNonPrompt);
void setCtauResFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb);
void setCtauResGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth);
bool setCtauResModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbPb);


bool fitCharmoniaCtauResModel( RooWorkspace& myws,             // Local Workspace
                               const RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                               struct KinCuts& cut,            // Variable containing all kinematic cuts
                               map<string, string>&  parIni,   // Variable containing all initial parameters
                               struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                               string outputDir,               // Path to output directory
                               // Select the type of datasets to fit
                               string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                               bool isPbPb        = false,     // isPbPb = false for pp, true for PbPb
                               bool importDS      = true,      // Select if the dataset is imported in the local workspace
                               // Select the type of object to fit
                               bool incJpsi       = true,      // Includes Jpsi model
                               bool incPsi2S      = true,      // Includes Psi(2S) model
                               // Select the fitting options
                               bool useTotctauErrPdf = false,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                               bool doFit         = true,      // Flag to indicate if we want to perform the fit
                               bool wantPureSMC   = false,     // Flag to indicate if we want to fit pure signal MC
                               bool loadFitResult = false,     // Load previous fit results
                               map<string, string> inputFitDir={},// User-defined Location of the fit results
                               int  numCores      = 2,         // Number of cores used for fitting
                               // Select the drawing options
                               bool setLogScale   = true,      // Draw plot with log scale
                               bool incSS         = false,     // Include Same Sign data
                               map<string, double> binWidth={} // User-defined Location of the fit results
                               )  
{

  if (DSTAG.find("_")!=std::string::npos) DSTAG.erase(DSTAG.find("_"));

  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    if (incJpsi && incPsi2S) { 
      cout << "[ERROR] We can only fit one type of signal using MC" << endl; return false; 
    }
    isMC = true;
  }
  if (!isMC)  { cout << "[ERROR] The ctau resolution fit can only be performed in MC!" << endl; return false; }
  bool incNonPrompt = (DSTAG.find("NOPR")!=std::string::npos);

  string COLL = (isPbPb ? "PbPb" : "PP" );
  bool usePerEventError = true;
  bool incBkg = false;

  if (importDS) {
    setCtauResCutParameters(cut, incNonPrompt);
    if (usePerEventError) {
      // check if we have already done the ctauErr fits. If yes, load their parameters
      string FileName = "";
      string pdfName = Form("pdfCTAUERR_Tot_%s", COLL.c_str());
      string plotLabel = "";
      bool incJpsi = true;
      bool fitSideBand = false;
      if (incJpsi)  { plotLabel = plotLabel + "_Jpsi";     }
      if (incPsi2S) { plotLabel = plotLabel + "_Psi2S";    }
      plotLabel = plotLabel + "_Bkg";
      setCtauErrFileName(FileName, (inputFitDir["CTAUERR"]=="" ? outputDir : inputFitDir["CTAUERR"]), "DATA", plotLabel, cut, isPbPb, fitSideBand);
      bool foundFit = false;
      if ( loadCtauErrRange(FileName, cut) ) { foundFit = true; }
      if (foundFit) { cout << "[INFO] The ctauErr fit was found and I'll load the ctau Error range used." << endl; }
      else { cout << "[ERROR] The ctauErr fit was not found!" << endl; return false; }
    }
    setCtauErrCutParameters(cut);
  }

  // Import the local datasets
  double numEntries = 1000000;
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  if (wantPureSMC) label = Form("%s_NoBkg", label.c_str());
  string dsName = Form("dOS_%s", label.c_str());
  if (importDS) {
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset(myws, inputWorkspace, cut, label);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doFit = false; }
  }
  else if (doFit && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }

  if (importDS) { 
    // Set global parameters
    setCtauResGlobalParameterRange(myws, parIni, cut, label, binWidth["CTAURES"]);
    RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName.c_str())->reduce(parIni["CtauNResRange_Cut"].c_str()))->Clone((dsName+"_CTAUNRESCUT").c_str());
    myws.import(*dataToFit);
  }
  string dsNameCut = dsName+"_CTAUNRESCUT";
  string pdfName = Form("pdfCTAURES_Tot_%s", (isPbPb?"PbPb":"PP"));

  if (!loadFitResult) {
    // Set models based on initial parameters
    struct OniaModel model;
    if (!setCtauResModel(model, parIni, isPbPb)) { return false; }
    //// LOAD CTAU ERROR PDF
    if (usePerEventError) {
      // Setting extra input information needed by each fitter
      bool loadCtauErrFitResult = true;
      bool doCtauErrFit = true;
      bool importDS = isMC;
      bool incJpsi = true;
      string DSTAG = Form("DATA_%s", (isPbPb?"PbPb":"PP"));
  
      if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                      DSTAG, isPbPb, importDS, 
                                      incJpsi, incPsi2S, incBkg, 
                                      doCtauErrFit, wantPureSMC, loadCtauErrFitResult, inputFitDir, numCores, 
                                      setLogScale, incSS, binWidth
                                      ) 
           ) { return false; }
    }

    // Build the Fit Model
    if (!buildCharmoniaCtauResModel(myws, (isPbPb ? model.PbPb : model.PP), parIni, dsName, "ctauRes", "pdfCTAURES", isPbPb, usePerEventError, useTotctauErrPdf, numEntries))  { return false; }
    if (!buildCharmoniaCtauResModel(myws, (isPbPb ? model.PbPb : model.PP), parIni, dsName, "ctauNRes", "pdfCTAUNRES", isPbPb, false, useTotctauErrPdf, numEntries))  { return false; }

    // save the initial values of the model we've just created
    RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauRes"), *myws.var("ctauNRes"), *myws.var("ctauErr"), *myws.var("ctau")));
    myws.saveSnapshot(Form("%s_parIni", pdfName.c_str()),*params,kTRUE);
    delete params;
  }

  // Define pdf and plot names
  string plotLabel = "";
  if (incJpsi || incPsi2S) { plotLabel = plotLabel + Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str()); }
  if (wantPureSMC)         { plotLabel = plotLabel + "_NoBkg"; }

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut, isPbPb);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    if (loadFitResult) return false;
    setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
  }
  bool found =  true; bool skipFit = !doFit;
  RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"), *myws.var("ctauNRes"), *myws.var("ctauRes")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (loadFitResult) {
    if ( loadPreviousFitResult(myws, FileName, DSTAG, isPbPb, false, false) ) { skipFit = true; } else  { skipFit = false; }
    if (skipFit) { cout << "[INFO] This ctau fit was already done, so I'll load the fit results." << endl; }
    myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()),*newpars,kTRUE);
  } else if (found) {
    cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
    return true;
  }

  // Fit the Datasets
  if (skipFit==false) {
    bool isWeighted = myws.data(dsName.c_str())->isWeighted();
    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsNameCut.c_str()), Extended(kTRUE), NumCPU(numCores), SumW2Error(isWeighted), Save());
    fitResult->Print("v");
    myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
    // Draw the mass plot
    drawCtauResPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbPb, wantPureSMC, setLogScale, incSS, binWidth["CTAURES"]);
    // Save the results
    string FileName = ""; setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctauRes/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);
  }
  
  return true;
};


bool setCtauResModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbPb)
{
  if (isPbPb) {
    if (parIni.count("Model_CtauRes_PbPb")>0) {
      model.PbPb.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PbPb"]];
      if (model.PbPb.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PbPb was not found in the initial parameters!" << endl; return false;
    }
  } else {
    if (parIni.count("Model_CtauRes_PP")>0) {
      model.PP.CtauRes = CtauModelDictionary[parIni["Model_CtauRes_PP"]];
      if (model.PP.CtauRes==CtauModel(0)) {
        cout << "[ERROR] The ctau resolution model: " << parIni["Model_CtauRes_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Resolution model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};


void setCtauResGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label, double binWidth)
{
  Double_t ctauNResMax; Double_t ctauNResMin;
  myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctauNRes"), ctauNResMin, ctauNResMax);
  ctauNResMin -= 0.00001;  ctauNResMax += 0.00001;
  cout << "[INFO] Range from data: ctauNResMin: " << ctauNResMin << "  ctauNResMax: " << ctauNResMax << endl;
  int nBins = min(int( round((ctauNResMax - ctauNResMin)/binWidth) ), 1000);
  TH1D* hTot = (TH1D*)myws.data(Form("dOS_%s", label.c_str()))->createHistogram("TMP", *myws.var("ctauNRes"), Binning(nBins, ctauNResMin, ctauNResMax));
  vector<double> rangeCtauNRes;
  getRange(hTot, (int)(ceil(5)), rangeCtauNRes);
  hTot->Delete();
  ctauNResMin = rangeCtauNRes[0];
  ctauNResMax = rangeCtauNRes[1];
  if (abs(ctauNResMax)>abs(ctauNResMin)) { ctauNResMax = abs(ctauNResMin); } else { ctauNResMin = -1.0*abs(ctauNResMax); }
  if (ctauNResMin<cut.dMuon.ctauNRes.Min) { ctauNResMin = cut.dMuon.ctauNRes.Min; }
  if (ctauNResMax>cut.dMuon.ctauNRes.Max) { ctauNResMax = cut.dMuon.ctauNRes.Max; }
  if (ctauNResMin<-10.0){ ctauNResMin = -10.0; }
  if (ctauNResMax>10.0) { ctauNResMax = 10.0;  }
  if (parIni.count("ctauResCut")>0 && parIni["ctauResCut"]!="") {
    parIni["ctauResCut"].erase(parIni["ctauResCut"].find("["), string("[").length());
    parIni["ctauResCut"].erase(parIni["ctauResCut"].find("]"), string("]").length());
    double ctauResCut = 0.0;
    try {
      ctauResCut = std::stod(parIni["ctauResCut"].c_str());
    } catch (const std::invalid_argument&) {
      cout << "[WARNING] ctauResCut is not a number, will ignore it!" << endl;
    }
    ctauResCut = abs(ctauResCut); if(ctauResCut>0.0) { ctauNResMax = abs(ctauResCut); ctauNResMin = -1.0*abs(ctauResCut); }
  }
  cout << "[INFO] Selected range: ctauNResMin: " << ctauNResMin << "  ctauNResMax: " << ctauNResMax << endl;
  myws.var("ctauNRes")->setRange("CtauNResWindow", ctauNResMin, ctauNResMax);
  parIni["CtauNResRange_Cut"]   = Form("(%.12f <= ctauNRes && ctauNRes <= %.12f)", ctauNResMin, ctauNResMax);
  cut.dMuon.ctauNRes.Max = ctauNResMax;
  cut.dMuon.ctauNRes.Min = ctauNResMin;
  myws.var("ctauNRes")->setRange("FullWindow", cut.dMuon.ctauNRes.Min, cut.dMuon.ctauNRes.Max);
  Double_t ctauResMax; Double_t ctauResMin;
  myws.data(Form("dOS_%s", label.c_str()))->reduce(parIni["CtauNResRange_Cut"].c_str())->getRange(*myws.var("ctauRes"), ctauResMin, ctauResMax);
  ctauResMin -= 0.00001;  ctauResMax += 0.00001;
  cout << "[INFO] Range from data: ctauResMin: " << ctauResMin << "  ctauResMax: " << ctauResMax << endl;
  myws.var("ctauRes")->setRange("CtauResWindow", ctauResMin, ctauResMax);
  parIni["CtauResRange_Cut"]   = Form("(%.12f <= ctauRes && ctauRes <= %.12f)", ctauResMin, ctauResMax);
  cut.dMuon.ctauRes.Max = ctauResMax;
  cut.dMuon.ctauRes.Min = ctauResMin;
  myws.var("ctauRes")->setRange("FullWindow", cut.dMuon.ctauRes.Min, cut.dMuon.ctauRes.Max);

  return;
};


void setCtauResFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauRes/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), "CTAURES", TAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  
  return;
};
 

void setCtauResCutParameters(struct KinCuts& cut, bool incNonPrompt)
{
  if (cut.dMuon.ctauNRes.Min==-100000.0 && cut.dMuon.ctauNRes.Max==100000.0) {
    // Default ctau values, means that the user did not specify a ctau Normalized Resolution range
    cut.dMuon.ctauNRes.Min = -30.0;
    cut.dMuon.ctauNRes.Max = 30.0;
  }
  // Define the ctau true range
  if (cut.dMuon.ctauTrue.Min==-1000.0) { cut.dMuon.ctauTrue.Min = (incNonPrompt ? 0.1 : -5.0); } // Removes cases with ctauTrue = -99
  if (cut.dMuon.ctau.Min==-1000.0)     { cut.dMuon.ctau.Min = -5.0;     } // Removes cases with ctau = -99
  
  return;
};


#endif // #ifndef fitCharmoniaCtauResModel_C
