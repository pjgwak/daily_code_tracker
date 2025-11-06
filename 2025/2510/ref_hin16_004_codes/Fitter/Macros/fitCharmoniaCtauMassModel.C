#ifndef fitCharmoniaCtauMassModel_C
#define fitCharmoniaCtauMassModel_C

#include "Utilities/initClasses.h"
#include "fitCharmoniaMassModel.C"
#include "fitCharmoniaCtauModel.C"
#include "fitCharmoniaCtauErrModel.C"
#include "fitCharmoniaCtauTrueModel.C"
#include "fitCharmoniaCtauRecoModel.C"
#include "drawMassFrom2DPlot.C"
#include "drawCtauFrom2DPlot.C"
#include "drawCtauMass2DPlot.C"

void setCtauMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb);


bool fitCharmoniaCtauMassModel( RooWorkspace& myws,             // Local Workspace
                                const RooWorkspace& inputWorkspace,   // Workspace with all the input RooDatasets
                                struct KinCuts& cut,            // Variable containing all kinematic cuts
                                map<string, string>&  parIni,   // Variable containing all initial parameters
                                struct InputOpt& opt,           // Variable with run information (kept for legacy purpose)
                                string outputDir,               // Path to output directory
                                // Select the type of datasets to fit
                                string DSTAG,                   // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
                                bool isPbPb        = false,     // isPbPb = false for pp, true for PbPb
                                // Select the type of object to fit
                                bool incJpsi       = true,      // Includes Jpsi model
                                bool incPsi2S      = false,     // Includes Psi(2S) model
                                // Select the fitting options
                                bool useTotctauErrPdf = false,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                                bool usectauBkgTemplate = false,// If yes use a template for Bkg ctau instead of the fitted Pdf
                                bool useCtauRecoPdf = false,    // If yes use the ctauReco PDF (template) instead of ctauTrue one
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
  if (DSTAG.find("MC")!=std::string::npos) {
    cout << "[ERROR] We can not make 2D fits on MC, change your input settings!" << endl; return false;
  }
  // Check if user has selected an object
  if (!incJpsi && !incPsi2S) {
    cout << "[ERROR] We need to include signal to perform 2D fits, change your input settings!" << endl; return false;
  }

  string pdfType = "pdfCTAUMASS";
  string COLL = (isPbPb ? "PbPb" : "PP" );
  string pdfName = Form("%s_Tot_%s", pdfType.c_str(), COLL.c_str());

  bool isMC          = false;
  bool incBkg        = true;      // Includes Background model
  bool incPrompt     = true;      // Includes Prompt ctau model
  bool incNonPrompt  = true;     // Includes NonPrompt ctau model
  bool fitSideBand   = false;
  bool usePerEventError = true;
  bool wantPureSMC = false;
  bool fitMass = true;
  bool fitCtau = true;

  // Set all the cuts on data
  setMassCutParameters(cut, incJpsi, incPsi2S, isMC, true);
  if (usePerEventError) {
    // check if we have already done the ctauErr fits. If yes, load their parameters
    string FileName = "";
    string pdfName = Form("pdfCTAUERR_Tot_%s", COLL.c_str());
    string plotLabel = "";
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
  setCtauCutParameters(cut, incNonPrompt);

  // Import the local datasets
  string label = ((DSTAG.find(COLL.c_str())!=std::string::npos) ? DSTAG.c_str() : Form("%s_%s", DSTAG.c_str(), COLL.c_str()));
  string dsName = Form("dOS_%s", label.c_str());
  if ( !(myws.data(dsName.c_str())) ) {
    int importID = importDataset(myws, inputWorkspace, cut, label, fitSideBand);
    if (importID<0) { return false; }
    else if (importID==0) { cout << "[ERROR] Local dataset failed to import!" << endl; return false; }
  }
  double numEntries = myws.data(dsName.c_str())->sumEntries(); 
  if ((numEntries<=0) || !(myws.data(dsName.c_str()))) { cout << "[ERROR] No local dataset was found to perform the fit!" << endl; return false; }
 
  // Set global parameters
  setMassGlobalParameterRange(myws, parIni, cut, incJpsi, incPsi2S, incBkg, false);
  setCtauErrGlobalParameterRange(myws, parIni, cut, "", binWidth["CTAUERR"], true);
  setCtauGlobalParameterRange(myws, parIni, cut, label, binWidth["CTAU"], (usectauBkgTemplate&&!isPbPb),true);

  // Cut the RooDataSet
  string dsNameCut = dsName+"_CTAUCUT";
  RooDataSet* dataToFit = (RooDataSet*)(myws.data(dsName.c_str())->reduce(parIni["CtauRange_Cut"].c_str()))->Clone(dsNameCut.c_str());
  myws.import(*dataToFit, Rename(dsNameCut.c_str()));
  double lostEvents = (myws.data(dsName.c_str())->numEntries()-myws.data(dsNameCut.c_str())->numEntries());
  double lostPerc = lostEvents*100./(myws.data(dsName.c_str())->numEntries());
  cout << "[INFO] After applying the cut: " << parIni["CtauRange_Cut"] << " , we lost " << lostEvents << " events ( " << lostPerc << " % ) " << endl;
  
  // Set models based on initial parameters
  struct OniaModel model;
  if (!setCtauModel(model, parIni, isPbPb, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt)) { return false; }

  string plotLabel = "";
  map<string, bool> plotLabels = {{"Jpsi",      (incJpsi)},
                                  {"Psi2S",     (incPsi2S)},
                                  {"Bkg",       (incBkg)},
                                  {"JpsiNoPR",  (incJpsi&&incNonPrompt)}, 
                                  {"Psi2SNoPR", (incPsi2S&&incNonPrompt)}, 
                                  {"BkgNoPR",   (incBkg&&incNonPrompt)}, 
                                  {"CtauRes",   (true)}};
  for (map<string, bool>::iterator iter = plotLabels.begin(); iter!=plotLabels.end(); iter++) {
    string obj = iter->first;
    bool cond = iter->second;
    if (cond && parIni.count(Form("Model_%s_%s", obj.c_str(), COLL.c_str()))>0) { 
      plotLabel = plotLabel + Form("_%s_%s", obj.c_str(), parIni[Form("Model_%s_%s", obj.c_str(), COLL.c_str())].c_str()); 
    }
  }

  string FileName = "";
  setCtauMassFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);

  //// LOAD MASS PDF
  if (fitMass) {
    // Setting extra input information needed by each fitter
    string iMassFitDir = inputFitDir["MASS"];
    double ibWidth = binWidth["MASS"];
    uint loadMassFitResult = (loadYields(myws, FileName, dsName, pdfName) ? 2 : 1);
    bool doMassFit = true;
    bool importDS = false;
    bool getMeanPT = false;
    bool zoomPsi = false;
    const char* applyCorr = "";
    bool doSimulFit = false;
    bool cutCtau = false;
    bool doConstrFit = false;

    if ( !fitCharmoniaMassModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                 DSTAG, isPbPb, importDS,
                                 incJpsi, incPsi2S, true, 
                                 doMassFit, cutCtau, doConstrFit, doSimulFit, wantPureSMC, applyCorr, loadMassFitResult, iMassFitDir, numCores, 
                                 setLogScale, incSS, zoomPsi, ibWidth, getMeanPT
                                 ) 
         ) { return false; }

    // Let's set all mass parameters to constant except the yields
    if (myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))) {
      cout << "[INFO] Setting mass parameters to constant!" << endl;
      myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("invMass")))->setAttribAll("Constant", kTRUE); 
    } else { cout << "[ERROR] Mass PDF was not found!" << endl; return false; }
    std::vector< std::string > objs = {"Bkg", "Jpsi", "Psi2S"};
    for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP"))))  setConstant( myws, Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")), false); }
    // Let's set the minimum value of the yields to zero
    for (auto obj : objs) { if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))) myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")))->setMin(0.0); }
    // Let's constrain the Number of Signal events to the updated mass fit results
    RooArgSet pdfList = RooArgSet("ConstrainPdfList");
    for (auto obj : objs) {
      if (myws.var(Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP"))))  {
        myws.factory(Form("Gaussian::%s_Gauss(%s,%s_Mean[%f],%s_Sigma[%f])", Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")), Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")),
                          Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")), myws.var(Form("N_%s_%s", obj.c_str(), isPbPb?"PbPb":"PP"))->getValV(),
                          Form("N_%s_%s", obj.c_str(), (isPbPb?"PbPb":"PP")) , myws.var(Form("N_%s_%s", obj.c_str(), isPbPb?"PbPb":"PP"))->getError()));
        pdfList.add(*myws.pdf(Form("N_%s_%s_Gauss", obj.c_str(), (isPbPb?"PbPb":"PP"))), kFALSE);
        std::cout << "[INFO] Constraining N_" << obj << " with Mean : " << myws.var(Form("N_%s_%s_Mean", obj.c_str(), (isPbPb?"PbPb":"PP")))->getVal()
                  << " and Sigma: " << myws.var(Form("N_%s_%s_Sigma", obj.c_str(), (isPbPb?"PbPb":"PP")))->getVal() << std::endl;
      }
    }
    myws.defineSet("ConstrainPdfList", pdfList);
    if (!myws.set("ConstrainPdfList")) { std::cout << "[ERROR] No constrain pdf list was found" << std::endl; return false; }
  }
  
  //// LOAD CTAU ERROR PDF
  if (usePerEventError) {
    // Setting extra input information needed by each fitter
    bool loadCtauErrFitResult = true;
    bool doCtauErrFit = false;
    bool importDS = false;
    bool wantPureSMC = false;
    
    if ( !fitCharmoniaCtauErrModel( myws, inputWorkspace, cut, parIni, opt, outputDir, 
                                    DSTAG, isPbPb, importDS, 
                                    incJpsi, incPsi2S, incBkg, 
                                    doCtauErrFit, wantPureSMC, loadCtauErrFitResult, inputFitDir, numCores, 
                                    setLogScale, incSS, binWidth
                                    ) 
         ) { return false; }
  }

  //// LOAD CTAU SIDEBAND PDF
  if (fitCtau) {
    // check if we have already done the sideband fits. If yes, load their results
    string FileName = "";
    bool fitSB = true;
    string plotLabel = Form("_BkgNoPR_%s_CtauRes_%s", parIni[Form("Model_BkgNoPR_%s", COLL.c_str())].c_str(), parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
    string DSTAG = Form("DATA_%s", (isPbPb?"PbPb":"PP"));
    setCtauFileName(FileName, (inputFitDir["CTAUSB"]=="" ? outputDir : inputFitDir["CTAUSB"]), DSTAG, plotLabel, cut, isPbPb, fitSB, usectauBkgTemplate);
    bool found = false;
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUSB"]!="") {
      setCtauFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb, fitSB, usectauBkgTemplate);
    } else if (inputFitDir["CTAUSB"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
      return false;
    }
    if (usectauBkgTemplate) {
      vector<string> pdfNames; pdfNames.push_back(Form("pdfCTAUCOND_Bkg_%s", COLL.c_str()));
      if ( !isPdfAlreadyFound(myws, FileName, pdfNames, true) ) {
        cout << "[ERROR] The ctau sideband fit results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau sideband fits were found, so I'll load the fit results." << endl; 
      }
    }
  }

  //// LOAD CTAU RECO PDF RESULTS
  if (fitCtau) {
    if (useCtauRecoPdf) {
      // check if we have already done the ctau Reco fits. If yes, load their results
      string FileName = "";
      string plotLabel = "_CtauReco";
      string DSTAG = Form("MCJPSINOPR_%s", (isPbPb?"PbPb":"PP"));
      string pdfName = Form("pdfCTAU_%s_%s", incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP"));
      string newPdfName = Form("pdfCTAUCOND_%s_%s", incJpsi?"JpsiNoPR":"Psi2SNoPR", (isPbPb?"PbPb":"PP"));
      setCtauRecoFileName(FileName, (inputFitDir["CTAURECO"]=="" ? outputDir : inputFitDir["CTAURECO"]), DSTAG, plotLabel, cut, isPbPb);
      bool found = false;
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURECO"]!="") {
        plotLabel = "_CtauReco_NoBkg";
        setCtauRecoFileName(FileName, (inputFitDir["CTAURECO"]=="" ? outputDir : inputFitDir["CTAURECO"]), DSTAG, plotLabel, cut, isPbPb);
      } else if (inputFitDir["CTAURECO"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURECO"]!="") {
        plotLabel = "_CtauReco";
        setCtauRecoFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
      } else if (inputFitDir["CTAURECO"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        plotLabel = "_CtauReco_NoBkg";
        setCtauRecoFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
      } else { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
        return false;
      }
      if ( !isCtauRecoPdfAlreadyFound(myws, FileName, pdfName, true, newPdfName ) ) {
        cout << "[ERROR] The ctau reco template results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau reco templates were found, so I'll load them." << endl; 
      }
    }
  }

  // Build the Fit Model
  if (!buildCharmoniaCtauModel(myws, (isPbPb ? model.PbPb : model.PP), parIni, dsName, cut, isPbPb, incBkg, incJpsi, incPsi2S, incPrompt, incNonPrompt, useTotctauErrPdf, usectauBkgTemplate, binWidth["CTAUSB"], numEntries))  { return false; }

  //// LOAD CTAU RESOLUTION PDF
  if (fitCtau && usectauBkgTemplate) {
    // check if we have already done the resolution fits. If yes, load their results
    string FileName = "";
    string plotLabel = Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
    string DSTAG = Form("MCJPSIPR_%s", (isPbPb?"PbPb":"PP"));
    if (inputFitDir["CTAURES"].find("nonPrompt")!=std::string::npos) DSTAG = Form("MCJPSINOPR_%s", (isPbPb?"PbPb":"PP"));
    if (inputFitDir["CTAURES"].find("DataFits")!=std::string::npos) DSTAG = Form("DATA_%s", (isPbPb?"PbPb":"PP"));
    setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut, isPbPb);
    if (wantPureSMC) { plotLabel = plotLabel + "_NoBkg"; }
    bool found = false;
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
      plotLabel = string(Form("_CtauRes_%s_NoBkg", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str()));
      setCtauResFileName(FileName, (inputFitDir["CTAURES"]=="" ? outputDir : inputFitDir["CTAURES"]), DSTAG, plotLabel, cut, isPbPb);
    } else if (inputFitDir["CTAURES"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAURES"]!="") {
      plotLabel = Form("_CtauRes_%s", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
      setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
    } else if (inputFitDir["CTAURES"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      plotLabel = string(Form("_CtauRes_%s_NoBkg", parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str()));
      setCtauResFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
    } else { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
      return false;
    }
    if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbPb, false, false) ) {
      cout << "[ERROR] The ctau resolution fit results were not loaded!" << endl;
      return false;
    } else { 
      cout << "[INFO] The ctau resolution fits were found, so I'll load the fit results." << endl; 
    }
    myws.pdf(Form("pdfCTAU_JpsiPR_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE);
  }

  //// LOAD CTAU SIDEBAND PDF
  if (fitCtau) {
    // check if we have already done the sideband fits. If yes, load their results
    string FileName = "";
    bool fitSB = true;
    string plotLabel = Form("_BkgNoPR_%s_CtauRes_%s", parIni[Form("Model_BkgNoPR_%s", COLL.c_str())].c_str(), parIni[Form("Model_CtauRes_%s", COLL.c_str())].c_str());
    string DSTAG = Form("DATA_%s", (isPbPb?"PbPb":"PP"));
    setCtauFileName(FileName, (inputFitDir["CTAUSB"]=="" ? outputDir : inputFitDir["CTAUSB"]), DSTAG, plotLabel, cut, isPbPb, fitSB, usectauBkgTemplate);
    bool found = false;
    if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUSB"]!="") {
      setCtauFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb, fitSB, usectauBkgTemplate);
    } else if (inputFitDir["CTAUSB"]!="") { found = true; }
    if (!found && gSystem->AccessPathName(FileName.c_str())) {
      cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
      return false;
    }
    if (!usectauBkgTemplate) {
      if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbPb, false, false) ) {
        cout << "[ERROR] The ctau sideband fit results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau sideband fits were found, so I'll load the fit results." << endl; 
      }
      if (myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbPb?"PbPb":"PP")))) {
        cout << "[INFO] Setting Prompt Background parameters to constant!" << endl;
        myws.pdf(Form("pdfCTAU_BkgPR_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
      } else { cout << "[ERROR] Prompt Background PDF was not found!" << endl; return false; }
      if (myws.pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))) {
        cout << "[INFO] Setting NonPrompt Background parameters to constant!" << endl;
        myws.pdf(Form("pdfCTAU_BkgNoPR_%s", (isPbPb?"PbPb":"PP")))->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr")))->setAttribAll("Constant", kTRUE); 
      } else { cout << "[ERROR] NonPrompt Background PDF was not found!" << endl; return false; }
      if (!setConstant(myws, Form("b_Bkg_%s", COLL.c_str()), true)) { return false; }
      if (!setConstant(myws, Form("s1_CtauRes_%s", COLL.c_str()), true)) { return false; }
    }
  }

  //// LOAD CTAU TRUE PDF RESULTS
  if (fitCtau) {
    if (!useCtauRecoPdf) {
      // check if we have already done the ctau true fits. If yes, load their results
      string FileName = "";
      string ModelName = Form("Model_JpsiNoPR_%s", COLL.c_str());
      string plotLabel = Form("_CtauTrue_%s", parIni[ModelName.c_str()].c_str());
      string DSTAG = Form("MCJPSINOPR_%s", (isPbPb?"PbPb":"PP"));
      setCtauTrueFileName(FileName, (inputFitDir["CTAUTRUE"]=="" ? outputDir : inputFitDir["CTAUTRUE"]), DSTAG, plotLabel, cut, isPbPb);
      bool found = false;
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUTRUE"]!="") {
        plotLabel = string(Form("_CtauTrue_%s_NoBkg", parIni[ModelName.c_str()].c_str()));
        setCtauTrueFileName(FileName, (inputFitDir["CTAUTRUE"]=="" ? outputDir : inputFitDir["CTAUTRUE"]), DSTAG, plotLabel, cut, isPbPb);
      } else if (inputFitDir["CTAUTRUE"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str()) && inputFitDir["CTAUTRUE"]!="") {
        plotLabel = string(Form("_CtauTrue_%s", parIni[ModelName.c_str()].c_str()));
        setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
      } else if (inputFitDir["CTAUTRUE"]!="") { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        plotLabel = string(Form("_CtauTrue_%s_NoBkg", parIni[ModelName.c_str()].c_str()));
        setCtauTrueFileName(FileName, outputDir, DSTAG, plotLabel, cut, isPbPb);
      } else { found = true; }
      if (!found && gSystem->AccessPathName(FileName.c_str())) {
        cout << "[ERROR] User Input File : " << FileName << " was not found!" << endl;
        return false;
      }
      if ( !loadPreviousFitResult(myws, FileName, DSTAG, isPbPb, false, false) ) {
        cout << "[ERROR] The ctau true fit results were not loaded!" << endl;
        return false;
      } else { 
        cout << "[INFO] The ctau true fits were found, so I'll load the fit results." << endl; 
      }
    }
  }

  // save the initial values of the model we've just created
  RooArgSet* params = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("invMass"), *myws.var("ctauErr")));
  myws.saveSnapshot((pdfName+"_parIni").c_str(),*params,kTRUE);
  delete params;

  // check if we have already done this fit. If yes, do nothing and return true.
  bool found =  true;
  RooArgSet *newpars = (RooArgSet*) myws.pdf(pdfName.c_str())->getParameters(RooArgSet(*myws.var("ctau"), *myws.var("ctauErr"), *myws.var("invMass")));
  found = found && isFitAlreadyFound(newpars, FileName, pdfName.c_str());
  if (found) {
    cout << "[INFO] This ctau fit was already done, so I'll just go to the next one." << endl;
    return true;
  }

  // Fit the Datasets
  bool isWeighted = myws.data(dsName.c_str())->isWeighted();
  RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsNameCut.c_str()), Extended(kTRUE), ExternalConstraints(*myws.set("ConstrainPdfList")), NumCPU(numCores), SumW2Error(isWeighted), Save());
  fitResult->Print("v");
  myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
  // Draw the mass plot
  drawCtauFrom2DPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbPb, incJpsi, incPsi2S, incBkg, setLogScale, incSS, binWidth["CTAU"]);
  drawMassFrom2DPlot(myws, outputDir, opt, cut, parIni, plotLabel, DSTAG, isPbPb, incJpsi, incPsi2S, incBkg, setLogScale, incSS, binWidth["MASS"]);
  drawCtauMass2DPlot(myws, outputDir, cut, plotLabel, DSTAG, isPbPb, binWidth);
  // Save the results
  myws.saveSnapshot(Form("%s_parFit", pdfName.c_str()),*newpars,kTRUE);
  saveWorkSpace(myws, Form("%sctauMass/%s/result", outputDir.c_str(), DSTAG.c_str()), FileName);

  return true;
};


void setCtauMassFileName(string& FileName, string outputDir, string TAG, string plotLabel, struct KinCuts cut, bool isPbPb)
{
  if (TAG.find("_")!=std::string::npos) TAG.erase(TAG.find("_"));
  FileName = Form("%sctauMass/%s/result/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), "CTAUMASS", TAG.c_str(), (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);

  return;
};


#endif // #ifndef fitCharmoniaCtauMassModel_C
