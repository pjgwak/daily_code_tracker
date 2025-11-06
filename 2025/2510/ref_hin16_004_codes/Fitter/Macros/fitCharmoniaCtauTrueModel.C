#ifndef fitCharmoniaCtauTrueModel_C
#define fitCharmoniaCtauTrueModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauTrueModel.C"
#include "drawCtauTruePlot.C"

bool setCtauTrueModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbPb, bool incResol=false);
void setCtauTrueFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb);
void setCtauTrueGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label);
void setCtauTrueCutParameters(struct KinCuts& cut);


bool fitCharmoniaCtauTrueModel( RooWorkspace& myws,             // Local Workspace
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
                                bool incResol      = true,      // Includes Ctau True Resolution model
                                // Select the fitting options
                                bool doFit         = true,      // Flag to indicate if we want to perform the fit
                                bool wantPureSMC   = false,     // Flag to indicate if we want to fit pure signal MC
                                bool loadFitResult = false,     // Load previous fit results
                                string inputFitDir = "",        // Location of the fit results
                                int  numCores      = 2,         // Number of cores used for fitting
                                // Select the drawing options
                                bool setLogScale   = true,      // Draw plot with log scale
                                bool incSS         = false,     // Include Same Sign data
                                double binWidth    = 0.05       // Bin width used for plotting
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
  if (!isMC) { wantPureSMC=false; }

  setCtauTrueCutParameters(cut);

  string COLL = (isPbPb ? "PbPb" : "PP" );

  // Set models based on initial parameters
  struct OniaModel model;
  if (!setCtauTrueModel(model, parIni, isPbPb)) { return false; }

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

  // Set global parameters
  setCtauTrueGlobalParameterRange(myws, parIni, cut, label);

  // Build the Fit Model     
  if (!buildCharmoniaCtauTrueModel(myws, (isPbPb ? model.PbPb : model.PP), parIni, isPbPb, incJpsi, incPsi2S, numEntries))  { return false; }

  // Define pdf and plot names
  string pdfName = Form("pdfCTAUTRUE_Tot_%s", COLL.c_str());
  string plotLabel = "";
  if (incJpsi || incPsi2S) { plotLabel = plotLabel + Form("_CtauTrue_%s", parIni[Form("Model_CtauTrue_%s", COLL.c_str())].c_str());        }
  if (incResol)            { plotLabel = plotLabel + Form("_CtauTrueRes_%s", parIni[Form("Model_CtauTrueRes_%s", COLL.c_str())].c_str()) ; }
  if (wantPureSMC)         { plotLabel = plotLabel + "_NoBkg"; }

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauTrueFileName(FileName, (inputFitDir=="" ? outputDir : inputFitDir), DSTAG, plotLabel, cut, isPbPb);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    if (loadFitResult) return false;
    setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
  }
  bool found =  true; bool skipFit = !doFit;
  RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctauTrue")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (loadFitResult) {
    if ( loadPreviousFitResult(myws, FileName, DSTAG, isPbPb, false, false) ) { skipFit = true; } else  { skipFit = false; }
    if (skipFit) { cout << "[INFO] This ctauTrue fit was already done, so I'll load the fit results." << endl; }
    myws.saveSnapshot(Form("%s_parLoad", pdfName.c_str()),*newpars,kTRUE);
  } else if (found) {
    cout << "[INFO] This ctauTrue fit was already done, so I'll just go to the next one." << endl;
    return true;
  }

  // Fit the Datasets
  if (skipFit==false) {
    bool isWeighted = myws.data(dsName.c_str())->isWeighted();
    RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(isWeighted), Range("CtauTrueWindow"), NumCPU(numCores), Save());
    fitResult->Print("v");
    myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
    // Draw the mass plot
    int nBins = min(int( round((cut.dMuon.ctauTrue.Max - cut.dMuon.ctauTrue.Min)/binWidth) ), 1000);
    drawCtauTruePlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbPb, incJpsi, incPsi2S, wantPureSMC, setLogScale, incSS, nBins);
    // Save the results
    string FileName = ""; setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctauTrue/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);
  }
  
  return true;
};


bool setCtauTrueModel( struct OniaModel& model, map<string, string>&  parIni, bool isPbPb, bool incResol )
{
  if (isPbPb) {
    if (incResol) {
      if (parIni.count("Model_CtauTrueRes_PbPb")>0) {
        model.PbPb.CtauTrueRes = CtauModelDictionary[parIni["Model_CtauTrueRes_PbPb"]];
        if (model.PbPb.CtauTrueRes==CtauModel(0)) {
          cout << "[ERROR] The ctau truth resolution model: " << parIni["Model_CtauTrueRes_PbPb"] << " is invalid" << endl; return false;
        }
      } else { 
        cout << "[ERROR] Ctau Truth resolution model for PbPb was not found in the initial parameters!" << endl; return false;
      }
    } else {
      model.PbPb.CtauTrueRes = CtauModelDictionary["Delta"];
    }
    if (parIni.count("Model_CtauTrue_PbPb")>0) {
      model.PbPb.CtauTrue = CtauModelDictionary[parIni["Model_CtauTrue_PbPb"]];
      if (model.PbPb.CtauTrue==CtauModel(0)) {
        cout << "[ERROR] The ctau truth model: " << parIni["Model_CtauTrue_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau Truth model for PbPb was not found in the initial parameters!" << endl; return false;
    }
  } else {
    if (incResol) {
      if (parIni.count("Model_CtauTrueRes_PP")>0) {
        model.PP.CtauTrueRes = CtauModelDictionary[parIni["Model_CtauTrueRes_PP"]];
        if (model.PP.CtauTrueRes==CtauModel(0)) {
          cout << "[ERROR] The ctau truth resolution model: " << parIni["Model_CtauTrue_PP"] << " is invalid" << endl; return false;
        }
      } else { 
        cout << "[ERROR] Ctau truth resolution model for PP was not found in the initial parameters!" << endl; return false;
      }
    } else {
      model.PP.CtauTrueRes = CtauModelDictionary["Delta"];
    }
    if (parIni.count("Model_CtauTrue_PP")>0) {
      model.PP.CtauTrue = CtauModelDictionary[parIni["Model_CtauTrue_PP"]];
      if (model.PP.CtauTrue==CtauModel(0)) {
        cout << "[ERROR] The ctau truth model: " << parIni["Model_CtauTrue_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Ctau truth model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};


void setCtauTrueGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label)
{
  Double_t ctauTrueMax; Double_t ctauTrueMin;
  myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctauTrue"), ctauTrueMin, ctauTrueMax);
  ctauTrueMin -= 0.00001;  ctauTrueMax += 0.00001;
  if (ctauTrueMin<cut.dMuon.ctauTrue.Min) { ctauTrueMin = cut.dMuon.ctauTrue.Min; }
  if (ctauTrueMax>cut.dMuon.ctauTrue.Max) { ctauTrueMax = cut.dMuon.ctauTrue.Max; }
  cout << "Range from data: ctauTrueMin: " << ctauTrueMin << "  ctauTrueMax: " << ctauTrueMax << endl;
  myws.var("ctauTrue")->setRange("CtauTrueWindow", ctauTrueMin, ctauTrueMax);
  parIni["CtauTrueRange_Cut"]   = Form("(%.12f <= ctauTrue && ctauTrue < %.12f)", ctauTrueMin, ctauTrueMax);
  cut.dMuon.ctauTrue.Max = (double)(ceil(ctauTrueMax));
  cut.dMuon.ctauTrue.Min = (double)(floor(ctauTrueMin));

  return;
};


void setCtauTrueFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauTrue/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), "CTAUTRUE", TAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  
  return;
};
 

void setCtauTrueCutParameters(struct KinCuts& cut)
{
  // Define the ctau true range
  if (cut.dMuon.ctauTrue.Min==-1000.0 && cut.dMuon.ctauTrue.Max==1000.0) { 
    // Default ctau values, means that the user did not specify a ctau True range
    cut.dMuon.ctauTrue.Min = -10.0;
    cut.dMuon.ctauTrue.Max = 10.0;
  }
  
  return;
};


#endif // #ifndef fitCharmoniaCtauTrueModel_C
