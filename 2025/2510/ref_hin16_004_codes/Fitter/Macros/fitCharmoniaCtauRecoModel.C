#ifndef fitCharmoniaCtauRecoModel_C
#define fitCharmoniaCtauRecoModel_C

#include "Utilities/initClasses.h"
#include "buildCharmoniaCtauRecoModel.C"
#include "drawCtauRecoPlot.C"

void setCtauRecoFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb);
void setCtauRecoGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label);
void setCtauRecoCutParameters(struct KinCuts& cut);
bool isCtauRecoPdfAlreadyFound(RooWorkspace& myws, string FileName, string pdfName, bool loadCtauRecoPdf=false, string newPdfName="");


bool fitCharmoniaCtauRecoModel( RooWorkspace& myws,             // Local Workspace
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
                                bool doCtauRecoPdf = true,      // Flag to indicate if we want to perform the fit
                                bool wantPureSMC   = false,     // Flag to indicate if we want to fit pure signal MC
                                bool loadCtauRecoPdf = false,     // Load previous fit results
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

  setCtauRecoCutParameters(cut);

  string COLL = (isPbPb ? "PbPb" : "PP" );

  // Import the local datasets
  double numEntries = 1000000;
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  if (wantPureSMC) label = Form("%s_NoBkg", label.c_str());
  string dsName = Form("dOS_%s", label.c_str());
  if (importDS) {
    if ( !(myws.data(dsName.c_str())) ) {
      int importID = importDataset(myws, inputWorkspace, cut, label);
      if (importID<0) { return false; }
      else if (importID==0) { doCtauRecoPdf = false; }
    }
    numEntries = myws.data(dsName.c_str())->sumEntries(); if (numEntries<=0) { doCtauRecoPdf = false; }
  }
  else if (doCtauRecoPdf && !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the template!" << endl; return false; }

  // Set global parameters
  setCtauRecoGlobalParameterRange(myws, parIni, cut, label);
  
  // Define pdf and plot names
  string pdfName = Form("pdfCTAU_%s_%s", incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP"));
  string plotLabel = "";
  if (incJpsi || incPsi2S) { plotLabel = plotLabel + "_CtauReco";}// Form("_CtauReco_%s", parIni[Form("Model_CtauReco_%s", COLL.c_str())].c_str());        }
  if (wantPureSMC)         { plotLabel = plotLabel + "_NoBkg"; }

  // check if we have already done this fit. If yes, do nothing and return true.
  string FileName = "";
  setCtauRecoFileName(FileName, (inputFitDir=="" ? outputDir : inputFitDir), DSTAG, plotLabel, cut, isPbPb);
  if (gSystem->AccessPathName(FileName.c_str()) && inputFitDir!="") {
    cout << "[WARNING] User Input File : " << FileName << " was not found!" << endl;
    if (loadCtauRecoPdf) return false;
    setCtauRecoFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
  }
  bool found =  true; bool skipCtauRecoPdf = !doCtauRecoPdf;
  found = found && isCtauRecoPdfAlreadyFound(myws, FileName, pdfName, loadCtauRecoPdf);
  if (found) {
    if (loadCtauRecoPdf) {
      cout << "[INFO] This ctauReco Pdf was already made, so I'll load the pdf." << endl;
    } else {
      cout << "[INFO] This ctauReco Pdf was already made, so I'll just go to the next one." << endl;
    }
    return true;
  }

  // Fit the Datasets
  if (skipCtauRecoPdf==false) {
    // Create the ctau Reco Pdf
    // Build the Ctau Reco Template
    if (!buildCharmoniaCtauRecoModel(myws, parIni, cut, dsName.c_str(), incJpsi, incPsi2S, binWidth, numEntries))  { return false; }

    int nBins = min(int( round((cut.dMuon.ctauTrue.Max - cut.dMuon.ctauTrue.Min)/binWidth) ), 1000);
    // Draw the ctau reco plot
    drawCtauRecoPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbPb, incJpsi, incPsi2S, wantPureSMC, setLogScale, incSS, nBins);

    // Save the results
    string FileName = ""; setCtauRecoFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
    RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("ctau")));
    myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
    saveWorkSpace(myws, Form("%sctauReco/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);

  }
  
  return true;
};


void setCtauRecoGlobalParameterRange(RooWorkspace& myws, map<string, string>& parIni, struct KinCuts& cut, string label)
{
  Double_t ctauRecoMax; Double_t ctauRecoMin;
  myws.data(Form("dOS_%s", label.c_str()))->getRange(*myws.var("ctau"), ctauRecoMin, ctauRecoMax);
  ctauRecoMin -= 0.00001;  ctauRecoMax += 0.00001;
  if (ctauRecoMin<cut.dMuon.ctauTrue.Min) { ctauRecoMin = cut.dMuon.ctauTrue.Min; }
  if (ctauRecoMax>cut.dMuon.ctauTrue.Max) { ctauRecoMax = cut.dMuon.ctauTrue.Max; }
  cout << "[INFO] Range from data: ctauRecoMin: " << ctauRecoMin << "  ctauRecoMax: " << ctauRecoMax << endl;
  myws.var("ctau")->setRange("CtauRecoWindow", ctauRecoMin, ctauRecoMax);
  parIni["CtauRecoRange_Cut"]   = Form("(%.12f <= ctau && ctau < %.12f)", ctauRecoMin, ctauRecoMax);
  cut.dMuon.ctauTrue.Max = (double)(ceil(ctauRecoMax));
  cut.dMuon.ctauTrue.Min = (double)(floor(ctauRecoMin));

  return;
};


void setCtauRecoFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauReco/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), "CTAURECO", TAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  
  return;
};
 

void setCtauRecoCutParameters(struct KinCuts& cut)
{
  // Define the ctau true range
  if (cut.dMuon.ctauTrue.Min==-1000.0 && cut.dMuon.ctauTrue.Max==1000.0) { 
    // Default ctau values, means that the user did not specify a ctau True range
    cut.dMuon.ctauTrue.Min = -10.0;
    cut.dMuon.ctauTrue.Max = 10.0;
  }
  
  return;
};

bool isCtauRecoPdfAlreadyFound(RooWorkspace& myws, string FileName, string pdfName, bool loadCtauRecoPdf, string newPdfName)
{
  if (gSystem->AccessPathName(FileName.c_str())) {
    cout << "[INFO] Results not found for: " << FileName << endl;
    return false; // File was not found
  }
  TFile *file = new TFile(FileName.c_str());
  if (!file) return false;
  
  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    cout << "[INFO] Workspace not found in: " << FileName << endl;
    file->Close(); delete file;
    return false;
  }
  
  bool found = true;
  string dataName = pdfName;
  dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
  if ( !(ws->pdf(pdfName.c_str())) ) {
    cout << "[INFO] " << pdfName << " was not found in: " << FileName << endl; found = false;
  }
  if (loadCtauRecoPdf && found) {
    if (newPdfName!="") {
      myws.import(*(ws->pdf(pdfName.c_str())), RooFit::RenameVariable(pdfName.c_str(),newPdfName.c_str()));
      if (ws->data(dataName.c_str())) myws.import(*(ws->data(dataName.c_str())));
      pdfName = newPdfName;
    }
    if (myws.pdf(pdfName.c_str()))   { cout << "[INFO] Pdf " << pdfName << " succesfully imported!" << endl;       }
    else {  cout << "[ERROR] Pdf " << pdfName << " import failed!" << endl; found = false; }
    if (myws.data(dataName.c_str())) { cout << "[INFO] DataHist " << dataName << " succesfully imported!" << endl; }
  }

  delete ws;
  file->Close(); delete file;
  
  return found;
};

#endif // #ifndef fitCharmoniaCtauRecoModel_C
