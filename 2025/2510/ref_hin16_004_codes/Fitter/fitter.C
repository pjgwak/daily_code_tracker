#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/fitCharmonia.C"

bool checkSettings(bool fitData, bool fitMC, bool fitPbPb, bool fitPP, bool fitMass, bool fitCtau, bool fitCtauTrue, bool fitCtauReco, bool fitRes, bool doCtauErrPDF, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt, bool cutCtau, bool doConstrFit, bool doSimulFit, bool wantPureSMC, const char* applyCorr, bool setLogScale, bool zoomPsi, bool incSS, int numCores);

bool parseFile(string FileName, vector< map<string, string> >& data);
bool parseString(string input, string delimiter, vector<double>& output);

bool iniWorkEnv( map< string, vector<string> >& DIR, const string workDirName);
void findSubDir(vector<string>& dirlist, string dirname);
bool existDir(string dir);
bool readFile(string FileName, vector< vector<string> >& content, const int nCol=-1, int nRow=-1);
bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection);
bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni, bool isPbPb, bool doConstrFit);
bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector, bool isPbPb, bool doConstrFit);


void fitter(
            const string workDirName="Test", // Working directoryi
            bool useExtFiles  = true, // Use external fit files as input
            bool useExtDS     = false, // Use external data/mc DataSets
            bool usePeriPD    = false, // If yes, use the PERIPHERAL PD provided by the user
            // Select the type of datasets to fit
            bool fitData      = true,        // Fits Data datasets
            bool fitMC        = false,         // Fits MC datasets
            bool fitPbPb      = true,         // Fits PbPb datasets
            bool fitPP        = true,        // Fits PP datasets
            bool fitMass      = true,       // Fits invariant mass distribution
            bool fitCtau      = true,       // Fits ctau distribution
            bool fitCtauTrue  = false,         // Fits ctau true MC distribution
            bool fitCtauReco  = false,      // Fit ctau reco MC distribution
            bool doCtauErrPDF = false,         // If yes, it builds the Ctau Error PDFs from data
            bool fitRes       = false,         // If yes fits the resolution from Data or MC
            // Select the type of object to fit
            bool incJpsi      = true,          // Includes Jpsi model
            bool incPsi2S     = false,         // Includes Psi(2S) model
            bool incBkg       = true,         // Includes Background model
            bool incPrompt    = true,         // Includes Prompt ctau model
            bool incNonPrompt = true,          // Includes Non Prompt ctau model 
            // Select the fitting options
            bool useTotctauErrPdf = false,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
            bool usectauBkgTemplate = true,// If yes use a template for Bkg ctau instead of the fitted Pdf
            bool useCtauRecoPdf = false,     // If yes use the ctauReco PDF (template) instead of ctauTrue one
            bool cutCtau      = false,        // Apply prompt ctau cuts
            bool doConstrFit   = false,        // Do constrained fit
            bool doSimulFit   = false,        // Do simultaneous fit
            bool wantPureSMC  = false,        // Flag to indicate if we want to fit pure signal MC
            const char* applyCorr  = "",     // Apply weight to data for correction (Acceptance & Ef , l_J/psi eff...). No correction if empty variable.
            int  numCores     = 32,           // Number of cores used for fitting
            // Select the drawing options
            bool  setLogScale  = true,         // Draw plot with log scale
            bool  incSS        = false,        // Include Same Sign data
            bool  zoomPsi      = false         // Zoom Psi(2S) peak on extra pad
            )
{
  // -------------------------------------------------------------------------------
  // STEP 0: INITIALIZE THE FITTER WORK ENVIROMENT
  // The work enviroment is divided as follows:
  /*
    main |-> Macros: Contain all the macros
         |-> Input   |-> <WorkDir> : Contain Input File, Bin and Parameter List for a given work directory (e.g. 20160201)
	 |-> Output  |-> <WorkDir> : Contain Output Plots and Results for a given work directory (e.g. 20160201)
	 |-> DataSet : Contain all the datasets (MC and Data)
  */


  gROOT->ProcessLine(".L ./Macros/Utilities/RooExtCBShape.cxx+");

  map<string, double> binWidth;
  binWidth["MASS"]     = 0.025;
  binWidth["CTAU"]     = 0.100;
  binWidth["CTAUERR"]  = 0.0025;
  binWidth["CTAUTRUE"] = 0.025;
  binWidth["CTAURECO"] = 0.100;
  binWidth["CTAURES"]  = 0.25;
  binWidth["CTAUSB"]   = 0.0150;

  if (workDirName.find("Peri")!=std::string::npos) { usePeriPD = true; }

  map<string, string> inputFitDir;
  inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_nominal/" : "DataFitsMassCent_2CB_polBkg_nominal/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_minvRange/" : "DataFitsMassCent_2CB_polBkg_minvRange/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_LLR10/" : "DataFitsMassCent_2CB_polBkg_LLR10/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_LLR25/" : "DataFitsMassCent_2CB_polBkg_LLR25/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_constrained/" : "DataFitsMassCent_2CB_polBkg_constrained/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_expBkg/" : "DataFitsMassCent_2CB_expBkg/");
  //inputFitDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/MassFits/") + (usePeriPD ? "DataFitsMassPeri_CBG_polBkg/" : "DataFitsMassCent_CBG_polBkg/");
  
  inputFitDir["CTAUERR"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauErrFits/") + (usePeriPD ? "DataFitsPeri_ctauErr_nominal/" :  "DataFitsCent_ctauErr_nominal/");
  inputFitDir["CTAUTRUE"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauTrueFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauTrue/" :  "MCFitsCent_nonPrompt_ctauTrue/");
  inputFitDir["CTAURECO"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauRecoFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauReco/" :  "MCFitsCent_nonPrompt_ctauReco/");
  
  inputFitDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauResFits/") + (usePeriPD ? "DataFitsPeri_ctauRes_nominalErr/" :  "DataFitsCent_ctauRes_nominalErr/");
  //inputFitDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauResFits/") + (usePeriPD ? "MCFitsPeri_prompt_ctauRes/" :  "MCFitsCent_prompt_ctauRes/");
  //inputFitDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauResFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauRes/" :  "MCFitsCent_nonPrompt_ctauRes/");  
  
  inputFitDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_dataCTauRes/" : "DataFitsCent_ctauBkg_dataCTauRes/");
  //inputFitDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_promptCTauRes/" : "DataFitsCent_ctauBkg_promptCTauRes/");
  //inputFitDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Output/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_nonPromptCTauRes/" : "DataFitsCent_ctauBkg_nonPromptCTauRes/");  

  map<string, string> inputInitialFilesDir;
  inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_nominal/" :  "DataFitsMassCent_2CB_polBkg_nominal/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_minvRange/" : "DataFitsMassCent_2CB_polBkg_minvRange/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_LLR10/" : "DataFitsMassCent_2CB_polBkg_LLR10/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_LLR25/" : "DataFitsMassCent_2CB_polBkg_LLR25/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_constrained/" : "DataFitsMassCent_2CB_polBkg_constrained/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_expBkg/" : "DataFitsMassCent_2CB_expBkg/");
  //inputInitialFilesDir["MASS"]     = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_CBG_polBkg/" : "DataFitsMassCent_CBG_polBkg/");
  
  inputInitialFilesDir["CTAUTRUE"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauTrueFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauTrue/" :  "MCFitsCent_nonPrompt_ctauTrue/");
  //inputInitialFilesDir["CTAUTRUE"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauTrueFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauTrue_minvRange/" :  "MCFitsCent_nonPrompt_ctauTrue_minvRange/");  

  inputInitialFilesDir["CTAURECO"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauRecoFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauReco/" :  "MCFitsCent_nonPrompt_ctauReco/");
  //inputInitialFilesDir["CTAURECO"] = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauRecoFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauReco_minvRange/" :  "MCFitsCent_nonPrompt_ctauReco_minvRange/");  

  inputInitialFilesDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauResFits/") + (usePeriPD ? "DataFitsPeri_ctauRes_nominalErr/" :  "DataFitsCent_ctauRes_nominalErr/");
  //inputInitialFilesDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauResFits/") + (usePeriPD ? "MCFitsPeri_prompt_ctauRes/" :  "MCFitsCent_prompt_ctauRes/");
  //inputInitialFilesDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauResFits/") + (usePeriPD ? "MCFitsPeri_nonPrompt_ctauRes/" :  "MCFitsCent_nonPrompt_ctauRes/");  
  //inputInitialFilesDir["CTAURES"]  = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauResFits/") + (usePeriPD ? "DataFitsPeri_ctauRes_nominalErr_minvRange/" :  "DataFitsCent_ctauRes_nominalErr_minvRange/");  

  inputInitialFilesDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_dataCTauRes/" :  "DataFitsCent_ctauBkg_dataCTauRes/");
  //inputInitialFilesDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_promptCTauRes/" : "DataFitsCent_ctauBkg_promptCTauRes/");
  //inputInitialFilesDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_nonPromptCTauRes/" : "DataFitsCent_ctauBkg_nonPromptCTauRes/");
  //inputInitialFilesDir["CTAUSB"]   = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/CtauBkgFits/") + (usePeriPD ? "DataFitsPeri_ctauBkg_dataCTauRes_minvRange/" :  "DataFitsCent_ctauBkg_dataCTauRes_minvRange/");  

  inputInitialFilesDir["CTAU"]     = string("");//string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/") + (usePeriPD ? "DataFitsPeri_2D_nominal/" :  "DataFitsCent_2D_nominal/");
  inputInitialFilesDir["FILES"]    = string("/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/Input/MassFits/") + (usePeriPD ? "DataFitsMassPeri_2CB_polBkg_nominal/" : "DataFitsMassCent_2CB_polBkg_nominal/");
  
  
  map<string, string> inputDataSet;
  inputDataSet["DOUBLEMUON"] = "/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/DataSets/";
  inputDataSet["PERIPHERAL"] = "/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/DataSets/";
  inputDataSet["MONTECARLO"] = "/afs/cern.ch/work/j/jmartinb/public/JpsiRAA/DataSets/";

  if (workDirName.find("Peri")!=std::string::npos) { usePeriPD = true; }

  bool fitTest = (workDirName=="Test");
  if (doCtauErrPDF || fitTest) { inputFitDir["CTAUERR"] = ""; inputInitialFilesDir["CTAUERR"] = "";}
  if ((fitMass && !fitCtau) || fitTest) { inputFitDir["MASS"] = ""; inputInitialFilesDir["MASS"] = "";}
  if (fitCtauTrue || fitTest) { inputFitDir["CTAUTRUE"] = ""; inputInitialFilesDir["CTAUTRUE"] = ""; }
  if (fitCtauReco || fitTest) { inputFitDir["CTAURECO"] = ""; inputInitialFilesDir["CTAURECO"] = ""; }
  if (fitRes || fitTest) { inputFitDir["CTAURES"] = ""; inputInitialFilesDir["CTAURES"] = ""; }
  if ((fitCtau && fitData && incBkg && !incJpsi) || fitTest) { inputFitDir["CTAUSB"] = ""; inputInitialFilesDir["CTAUSB"] = ""; }

  if (!checkSettings(fitData, fitMC, fitPbPb, fitPP, fitMass, fitCtau, fitCtauTrue, fitCtauReco, fitRes, doCtauErrPDF, incJpsi, incPsi2S, incBkg, incPrompt, incNonPrompt, cutCtau, doConstrFit, doSimulFit, wantPureSMC, applyCorr, setLogScale, zoomPsi, incSS, numCores)) { return; }

  map< string, vector<string> > DIR;
  if(!iniWorkEnv(DIR, workDirName)){ return; }

  vector< map<string, string> > inputFitDirs;
  inputFitDirs.push_back(inputFitDir);
  for(uint i=1; i<DIR["input"].size(); i++) {
    inputFitDirs.push_back(inputFitDir);
    for(map<string, string>::iterator iter=inputFitDirs[i].begin(); iter!=inputFitDirs[i].end(); iter++) {
      string key = iter->first;
      if (inputFitDirs[i][key]!="") {
        inputFitDirs[i][key] = DIR["input"][i];
        inputFitDirs[i][key].replace(inputFitDirs[i][key].find(DIR["input"][0]), DIR["input"][0].length(), inputFitDirs[0][key]);
      }
    }
  }

  vector< map<string, string> > inputInitialFilesDirs;
  inputInitialFilesDirs.push_back(inputInitialFilesDir);
  for(uint i=1; i<DIR["input"].size(); i++) {
    inputInitialFilesDirs.push_back(inputInitialFilesDir);
    for(map<string, string>::iterator iter=inputInitialFilesDirs[i].begin(); iter!=inputInitialFilesDirs[i].end(); iter++) {
      string key = iter->first;
      if (inputInitialFilesDirs[i][key]!="") {
        inputInitialFilesDirs[i][key] = DIR["input"][i];
        inputInitialFilesDirs[i][key].replace(inputInitialFilesDirs[i][key].find(DIR["input"][0]), DIR["input"][0].length(), inputInitialFilesDirs[0][key]);
      }
    }
  }

  // -------------------------------------------------------------------------------
  // STEP 1: CREATE/LOAD THE ROODATASETS
  /*
    Input : List of TTrees with format:  TAG <tab> FILE_NAME
    Output: Collection of RooDataSets splitted by tag name, including OS and SS dimuons.
  */
  
  const string InputTrees = (useExtDS ? inputInitialFilesDir["FILES"] : DIR["input"][0]) + "InputTrees.txt";
  map<string, vector<string> > InputFileCollection;
  if(!getInputFileNames(InputTrees, InputFileCollection)){ return; }
  
  TObjArray* aDSTAG = new TObjArray(); // Array to store the different tags in the list of trees
  aDSTAG->SetOwner(true);
  map<string, RooWorkspace> Workspace;
  
  for(map<string, vector<string> >::iterator FileCollection=InputFileCollection.begin(); FileCollection!=InputFileCollection.end(); ++FileCollection) {
    // Get the file tag which has the following format: DSTAG_COLL , i.e. DATA_PP 
    string FILETAG = FileCollection->first;  
    string DSTAG   = FILETAG;
    if (FILETAG.size()) {
      if (doSimulFit) DSTAG.erase(DSTAG.find("_")); // if doSimulFit=true the workspace will be labeled symply as DSTAG, otherwise as DSTAG_COLL (to have PP and PbPb separated ws)
    } else {
      cout << "[ERROR] FILETAG is empty!" << endl;
    }
    // Extract the filenames
    vector<string> InputFileNames = FileCollection->second; 
    string         OutputFileName;
    // If we have data, check if the user wants to fit data
    bool checkData = (fitCtau && fitMC && (existDir(inputFitDirs[inputFitDirs.size()-1]["CTAUERR"]+"ctauErr/")==false) );
    if ( (FILETAG.find("DATA")!=std::string::npos) && (fitData==true || checkData==true) ) {
      if ( (FILETAG.find("PP")!=std::string::npos)   && !fitPP   ) continue; // If we find PP, check if the user wants PP
      if ( (FILETAG.find("PbPb")!=std::string::npos) && !fitPbPb ) continue; // If we find PbPb, check if the user wants PbPb
      string dir = DIR["dataset"][0];
      
      bool isPbPb = false, isPeriPD = false;
      if (FILETAG.find("PbPb")!=std::string::npos) { isPbPb = true; }
      if (InputFileNames[0].find("HIOniaPeripheral30100")!=std::string::npos) { isPeriPD = true; }
      
      if (isPbPb && usePeriPD && !isPeriPD){
        cout << "[ERROR] Requested to use Peripheral PD but InputTress.txt does not contain the corresponding trees" << endl;
        return;
      }
      
      if (useExtDS==true && usePeriPD==true  && inputDataSet["PERIPHERAL"]!="" && (existDir(inputDataSet["PERIPHERAL"])==true)) { dir = inputDataSet["PERIPHERAL"]; }
      if (useExtDS==true && usePeriPD==false && inputDataSet["DOUBLEMUON"]!="" && (existDir(inputDataSet["DOUBLEMUON"])==true)) { dir = inputDataSet["DOUBLEMUON"]; }
      if(strcmp(applyCorr,"")){
        OutputFileName = dir + "DATASET_" + FILETAG + "_" + string(applyCorr) + (isPbPb?(isPeriPD?"_PERI":"_CENT"):"") + ".root";
        if(gSystem->AccessPathName(OutputFileName.c_str())) { OutputFileName = DIR["dataset"][0] + "DATASET_" + FILETAG + "_" + string(applyCorr) + ".root"; }
        if(!tree2DataSet(Workspace[Form("%s_%s",DSTAG.c_str(),applyCorr)], InputFileNames, FILETAG, OutputFileName)){ return; }
      }
      else {
        OutputFileName = dir + "DATASET_" + FILETAG + (isPbPb?(isPeriPD?"_PERI":"_CENT"):"") + ".root";
        if(gSystem->AccessPathName(OutputFileName.c_str())) { OutputFileName = DIR["dataset"][0] + "DATASET_" + FILETAG + (isPbPb?(isPeriPD?"_PERI":"_CENT"):"") + ".root"; }
        string NAMETAG = DSTAG;
        if (checkData) { NAMETAG = string("MC")+(incJpsi?"JPSI":"PSI2S")+(incNonPrompt?"NOPR":"PR")+"_"+(fitPP?"PP":"PbPb"); }
        if(!tree2DataSet(Workspace[NAMETAG], InputFileNames, FILETAG, OutputFileName)){ return; }
      }
      if (fitData && !aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str()));
    }
    
    // If we find MC, check if the user wants to fit MC
    if ( (FILETAG.find("MC")!=std::string::npos) && fitMC==true  ) {
      if ( (FILETAG.find("PP")!=std::string::npos)    && !fitPP     ) continue; // If we find PP, check if the user wants PP
      if ( (FILETAG.find("PbPb")!=std::string::npos)  && !fitPbPb   ) continue; // If we find PbPb, check if the user wants PbPb
      if ( (FILETAG.find("JPSI")!=std::string::npos)  && !incJpsi   ) continue; // If we find Jpsi MC, check if the user wants to include Jpsi
      if ( (FILETAG.find("PSI2S")!=std::string::npos) && !incPsi2S  ) continue; // If we find Psi2S MC, check if the user wants to include Psi2S
      if ( (FILETAG.find("NOPR")!=std::string::npos) ) { if (!incNonPrompt) continue; } // If we find Non-Prompt MC, check if the user wants to include Non-Prompt
      else if ( (FILETAG.find("PR")!=std::string::npos) && !incPrompt ) continue; // If we find Prompt MC, check if the user wants to include Prompt
      string dir = DIR["dataset"][0];
      if (useExtDS==true && inputDataSet["MONTECARLO"]!="" && (existDir(inputDataSet["MONTECARLO"])==true)) { dir = inputDataSet["MONTECARLO"]; }
      OutputFileName = dir + "DATASET_" + FILETAG + ".root";
      if(gSystem->AccessPathName(OutputFileName.c_str())) { OutputFileName = DIR["dataset"][0] + "DATASET_" + FILETAG + ".root"; }
      if(!tree2DataSet(Workspace[DSTAG], InputFileNames, FILETAG, OutputFileName)){ return; }
      if (fitMC && !aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str()));
      if (wantPureSMC)
	{
	  OutputFileName = dir + "DATASET_" + FILETAG + "_PureS" + ".root";
	  if(!tree2DataSet(Workspace[Form("%s_PureS",DSTAG.c_str())], InputFileNames, FILETAG, OutputFileName)){ return; }
	}
    }
  }
  if (Workspace.size()==0) {
    cout << "[ERROR] No onia tree files were found matching the user's input settings!" << endl; return;
  }

  // -------------------------------------------------------------------------------
  // STEP 2: LOAD THE INITIAL PARAMETERS
  /*
    Input : List of initial parameters with format PT <tab> RAP <tab> CEN <tab> iniPar ... 
    Output: two vectors with one entry per kinematic bin filled with the cuts and initial parameters
  */
  
  string InputFile;
  vector< vector< struct KinCuts > >       cutVectors;
  vector< vector< map<string, string> > >  parIniVectors;
 
  map<string, map<string, bool>> VARMAP = {
    {"MASS", 
     {
       {"BKG",   ((fitMass  && incBkg) || (fitCtau || doCtauErrPDF || (fitRes && fitData)))}, 
       {"JPSI",  (((fitMass && incJpsi) || (fitCtau || doCtauErrPDF || (fitRes && fitData))) && incJpsi) || (fitCtau && incBkg)}, 
       {"PSI2S", ((fitMass && incPsi2S) || (fitCtau || doCtauErrPDF || (fitRes && fitData))) && incPsi2S}
     }
    },
    {"CTAU", 
     {
       {"BKG",   fitCtau && incBkg && incNonPrompt}, 
       {"JPSI",  fitCtau && incJpsi && incNonPrompt},
       {"PSI2S", fitCtau && incPsi2S && incNonPrompt},
       {"RES",   (fitCtau || fitRes)},
       {"TRUE",  (fitCtauTrue || fitCtauReco)},
     }
    }
  };
  map<string, bool> COLMAP = {{"PbPb", fitPbPb}, {"PP", fitPP}};
 
  typedef map<string, map<string, bool>>::iterator var_type;
  typedef map<string, bool>::iterator it_type;
  for(uint j = 0; j < DIR["input"].size(); j++) {
    if (DIR["input"].size()>1 && j==0) continue; // First entry is always the main input directory
    vector< struct KinCuts >       cutVector;
    vector< map<string, string> >  parIniVector;
    for(var_type VAR = VARMAP.begin(); VAR != VARMAP.end(); VAR++) {
      map<string, bool> PARMAP = VAR->second;
      string name1 = "InitialParam_" + VAR->first + "_";
      for(it_type PAR = PARMAP.begin(); PAR != PARMAP.end(); PAR++) {
        if (PAR->second) {
          string name2 = name1 + PAR->first + "_";
          for(it_type COL = COLMAP.begin(); COL != COLMAP.end(); COL++) {
            if(COL->second) {
              string dir = DIR["input"][j];
              if (VAR->first=="MASS" && inputInitialFilesDirs[j]["MASS"]!="") { dir = inputInitialFilesDirs[j]["MASS"]; }
              if (VAR->first=="CTAU" && PAR->first=="TRUE" && inputInitialFilesDirs[j]["CTAUTRUE"]!="") { dir = inputInitialFilesDirs[j]["CTAUTRUE"]; }
              if (VAR->first=="CTAU" && PAR->first=="RECO" && inputInitialFilesDirs[j]["CTAURECO"]!="") { dir = inputInitialFilesDirs[j]["CTAURECO"]; }
              if (VAR->first=="CTAU" && PAR->first=="RES"  && inputInitialFilesDirs[j]["CTAURES"]!="" ) { dir = inputInitialFilesDirs[j]["CTAURES"];  }
              if (VAR->first=="CTAU" && PAR->first=="BKG"  && inputInitialFilesDirs[j]["CTAUSB"]!=""  ) { dir = inputInitialFilesDirs[j]["CTAUSB"];   }
              if (VAR->first=="CTAU" && PAR->first=="JPSI" && inputInitialFilesDirs[j]["CTAU"]!=""    ) { dir = inputInitialFilesDirs[j]["CTAU"];     }
              string name3 = name2 + COL->first + ".csv";
              InputFile = (dir + name3);
              if (!addParameters(InputFile, cutVector, parIniVector, true, doConstrFit)) { return; }
            }
          }
        }
      }
    }
    cutVectors.push_back(cutVector);
    parIniVectors.push_back(parIniVector);
  }
  
  // -------------------------------------------------------------------------------  
  // STEP 3: FIT THE DATASETS
  /*
    Input : 
              -> The cuts and initial parameters per kinematic bin
	      -> The workspace with the full datasets included.
    Output: 
              -> Plots (png, pdf and root format) of each fit.
	      -> The local workspace used for each fit.
  */
  TIter nextDSTAG(aDSTAG);
  for(uint j = 0; j < cutVectors.size(); j++) {
    int index = ( DIR["output"].size()>1 ? j+1 : j ); // First entry is always the main output directory
    string outputDir = DIR["output"][index];

    for (unsigned int i=0; i < cutVectors[j].size(); i++) {
      nextDSTAG.Reset();
      TObjString* soDSTAG(0x0);
      while ( (soDSTAG = static_cast<TObjString*>(nextDSTAG.Next())) )
        {
          TString DSTAG = (TString)(soDSTAG->GetString());
          if ( fitMass && fitCtau && !DSTAG.Contains("DATA") ) continue;

          TString wsName = "";
          if (DSTAG.Contains("MC") && wantPureSMC) wsName = Form("%s_PureS",DSTAG.Data());
          else if (DSTAG.Contains("DATA") && strcmp(applyCorr,"")) wsName = Form("%s_%s",DSTAG.Data(),applyCorr);
          else wsName = DSTAG;
          
          if (Workspace.count(wsName.Data())>0) {
            
            // DATA/MC datasets were loaded
            if (doSimulFit) {
              // If do simultaneous fits, then just fits once
              if (!fitCharmonia( Workspace[DSTAG.Data()], cutVectors[j].at(i), parIniVectors[j].at(i), outputDir,
                                 // Select the type of datasets to fit
                                 DSTAG.Data(),
                                 false,           // dummy flag when fitting simultaneously since both PP and PbPb are used
                                 // Select the type of object to fit
                                 fitMass,         // Fit mass distribution
                                 fitCtau,         // Fit ctau distribution
                                 fitCtauTrue,     // Fits ctau true MC distribution
                                 fitCtauReco,     // Fit ctau reco MC distribution
                                 incJpsi,         // Includes Jpsi model
                                 incPsi2S,        // Includes Psi(2S) model
                                 incBkg,          // Includes Background model
                                 incPrompt,       // Includes Prompt ctau model
                                 incNonPrompt,    // Includes NonPrompt ctau model
                                 doCtauErrPDF,    // If yes, it builds the Ctau Error PDFs from data
                                 fitRes,          // If yes fits the resolution from Data or MC
                                 // Select the fitting options
                                 useTotctauErrPdf,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                                 usectauBkgTemplate,// If yes use a template for Bkg ctau instead of the fitted Pdf
                                 useCtauRecoPdf,  // If yes use the ctauReco PDF (template) instead of ctauTrue one
                                 cutCtau,         // Apply prompt ctau cuts
                                 doConstrFit,     // Do constrained fit
                                 doSimulFit,      // Do simultaneous fitC
                                 wantPureSMC,     // Flag to indicate if we want to fit pure signal MC
                                 inputFitDirs[index], // Map of user-defined directory paths of previous fit results
                                 applyCorr,       // Flag to indicate if we want corrected dataset and which correction
                                 numCores,        // Number of cores used for fitting
                                 // Select the drawing options
                                 setLogScale,     // Draw plot with log scale
                                 incSS,           // Include Same Sign data
                                 zoomPsi,         // Zoom Psi(2S) peak on extra pad
                                 binWidth,        // Bin width used for plotting
                                 false            // Compute the mean PT (NEED TO FIX)
                                 )
                  ) { return; }
            } else {
              // If don't want simultaneous fits, then fit PbPb or PP separately
              if ( DSTAG.Contains("MCJPSIPR")    ) { incJpsi = true;  incPsi2S = false; incPrompt = true;  incNonPrompt = false; }
              if ( DSTAG.Contains("MCJPSINOPR")  ) { incJpsi = true;  incPsi2S = false; incPrompt = false; incNonPrompt = true;  }
              if ( DSTAG.Contains("MCPSI2SPR")   ) { incJpsi = false; incPsi2S = true;  incPrompt = true;  incNonPrompt = false; }
              if ( DSTAG.Contains("MCPSI2SNOPR") ) { incJpsi = false; incPsi2S = true;  incPrompt = false; incNonPrompt = true;  }
              
              bool isPbPb = true, proceed = false;
              if (fitPbPb && DSTAG.Contains("PbPb")) { isPbPb = true;  proceed = true; }
              if (fitPP   && DSTAG.Contains("PP"))   { isPbPb = false; proceed = true; }
              if (proceed && !fitCharmonia( Workspace[wsName.Data()], cutVectors[j].at(i), parIniVectors[j].at(i), outputDir,
                                            // Select the type of datasets to fit
                                            DSTAG.Data(),
                                            isPbPb,          // In this case we are fitting PbPb
                                            // Select the type of object to fit
                                            fitMass,         // Fit mass distribution
                                            fitCtau,         // Fit ctau distribution
                                            fitCtauTrue,     // Fits ctau true MC distribution
                                            fitCtauReco,     // Fit ctau reco MC distribution
                                            incJpsi,         // Includes Jpsi model
                                            incPsi2S,        // Includes Psi(2S) model
                                            incBkg,          // Includes Background model
                                            incPrompt,       // Includes Prompt ctau model
                                            incNonPrompt,    // Includes NonPrompt ctau model
                                            doCtauErrPDF,    // If yes, it builds the Ctau Error PDFs from data
                                            fitRes,          // If yes fits the resolution from Data or MC
                                            // Select the fitting options
                                            useTotctauErrPdf,  // If yes use the total ctauErr PDF instead of Jpsi and bkg ones
                                            usectauBkgTemplate,// If yes use a template for Bkg ctau instead of the fitted Pdf
                                            useCtauRecoPdf,  // If yes use the ctauReco PDF (template) instead of ctauTrue one
                                            cutCtau,         // Apply prompt ctau cuts
                                            doConstrFit,     // Do constrained fit
                                            doSimulFit,      // Do simultaneous fit
                                            wantPureSMC,           // Flag to indicate if we want to fit pure signal MC
                                            inputFitDirs[index],// Map of user-defined directory paths of previous fit results
                                            applyCorr,       // Flag to indicate if we want corrected dataset and which correction
                                            numCores,        // Number of cores used for fitting
                                            // Select the drawing options
                                            setLogScale,     // Draw plot with log scale
                                            incSS,           // Include Same Sign data
                                            zoomPsi,         // Zoom Psi(2S) peak on extra pad
                                            binWidth,        // Bin width used for plotting
                                            false            // Compute the mean PT (NEED TO FIX)
                                            )
                  ) { return; }
            }
          } else {
            cout << "[ERROR] The workspace for " << wsName.Data() << " was not found!" << endl; return;
          }
        }
    }
  }

  aDSTAG->Delete();
  delete aDSTAG;
};
  

bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector, bool isPbPb, bool doConstrFit)
{
  vector< map<string, string> >  data;
  if(!parseFile(InputFile, data)) { return false; }
  
  if (cutVector.size()==0) {
    for(vector< map<string, string> >::iterator row=data.begin(); row!=data.end(); ++row) {
      struct KinCuts cut; map<string, string> parIni;
      if(!setParameters(*row, cut, parIni, isPbPb, doConstrFit)) { return false; }
      cutVector.push_back(cut);  parIniVector.push_back(parIni);
    }
  }
  else {
    if (data.size()!=cutVector.size()) { cout << "[ERROR] The initial parameters in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
    for (unsigned int i=0; i<data.size(); i++) {
      struct KinCuts cut;
      if (!setParameters(data.at(i), cut, parIniVector.at(i), isPbPb, doConstrFit)) { return false; };
      if (!isEqualKinCuts(cut, cutVector.at(i), isPbPb)) { cout << "[ERROR] The bins in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
    }
  }

  return true;
};


bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni, bool isPbPb, bool doConstrFit)
{

  // set initial parameters
  cut.sMuon.Pt.Min  =  0.0;
  cut.sMuon.Pt.Max  = 100000.0;
  cut.sMuon.Eta.Min = -2.4;
  cut.sMuon.Eta.Max = 2.4;
  cut.dMuon.ctauErr.Min = -1000.0;
  cut.dMuon.ctauErr.Max = 1000.0;
  cut.dMuon.ctau.Min = -1000.0;
  cut.dMuon.ctau.Max = 1000.0;
  cut.dMuon.ctauNRes.Min = -100000.0;
  cut.dMuon.ctauNRes.Max = 100000.0;
  cut.dMuon.ctauN.Min = -100000.0;
  cut.dMuon.ctauN.Max = 100000.0;
  cut.dMuon.ctauRes.Min = -1000.0;
  cut.dMuon.ctauRes.Max = 1000.0;
  cut.dMuon.ctauTrue.Min = -1000.0; 
  cut.dMuon.ctauTrue.Max = 1000.0;
  cut.dMuon.ctauCut = "";   
  cut.dMuon.M.Min = 2.0; 
  cut.dMuon.M.Max = 5.0;  
  cut.dMuon.AbsRap.Min = 0.0;
  cut.dMuon.AbsRap.Max = 2.4;
  cut.dMuon.Pt.Min  =  0.0;
  cut.dMuon.Pt.Max  =  1000.0;
  cut.Centrality.Start = 0;
  cut.Centrality.End = 200;

  // set parameters from file
  for(map<string, string>::iterator col=row.begin(); col!=row.end(); ++col) {
    string label = col->first;
    if (label=="rap") {
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'rap' has invalid value: " << col->second << endl; return false;
      }  
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'rap' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.AbsRap.Min = v.at(0); 
      cut.dMuon.AbsRap.Max = v.at(1);
    } 
    else if (label=="pt"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'pt' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'pt' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.Pt.Min = v.at(0); 
      cut.dMuon.Pt.Max = v.at(1);
    } 
    else if (label=="cent"){
      if (isPbPb) {
        if (col->second=="" || col->second.find("-")==std::string::npos) {
          cout << "[ERROR] Input column 'cent' has invalid value: " << col->second << endl; return false;
        }
        std::vector<double> v;
        if(!parseString(col->second, "-", v)) { return false; }
        if (v.size()!=2) {
          cout << "[ERROR] Input column 'cent' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
        }
        cut.Centrality.Start = (int) (2.0*v.at(0));
        cut.Centrality.End = (int) (2.0*v.at(1));
      }
    } 
    else if (label=="mass"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'mass' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'mass' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.M.Min = v.at(0); 
      cut.dMuon.M.Max = v.at(1);
    } 
    else if (label=="ctau"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctau' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctau' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctau.Min = v.at(0); 
      cut.dMuon.ctau.Max = v.at(1);
    }  
    else if (label=="ctauErr"){
      if (col->second=="" || col->second.find("-")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctauErr' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctauErr' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctauErr.Min = v.at(0); 
      cut.dMuon.ctauErr.Max = v.at(1);
    }
    else if (label=="ctauTrue"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctauTrue' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctauTrue' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctauTrue.Min = v.at(0); 
      cut.dMuon.ctauTrue.Max = v.at(1);
    }
    else if (label=="ctauRes"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctauRes' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctauRes' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctauRes.Min = v.at(0); 
      cut.dMuon.ctauRes.Max = v.at(1);
    }
    else if (label=="ctauNRes"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctauNRes' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctauNRes' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctauNRes.Min = v.at(0); 
      cut.dMuon.ctauNRes.Max = v.at(1);
    }
    else if (label=="ctauResCut"){ 
      parIni[col->first] = col->second; 
    }
    else if (label=="ctauN"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctauN' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctauN' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctauN.Min = v.at(0); 
      cut.dMuon.ctauN.Max = v.at(1);
    }
    else if (label=="ctauNCut"){ 
      parIni[col->first] = col->second; 
    }
    else if (label.find("Model")!=std::string::npos){
      if (col->second=="") {
        cout << "[ERROR] Input column 'Model' has empty value" << endl; return false;
      }
      parIni[col->first] = col->second;
    } 
    else {
      if (col->second != "") {
	string value = col->second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple

	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << endl; return false;
	} else {
	  value.erase(value.find("["), string("[").length());
	  value.erase(value.find("]"), string("]").length());
	}
	std::vector<double> v; 
	if(!parseString(value, ",", v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          cout << "[ERROR] Initial parameter " << col->first << " has incorrect number of values, it has: " << v.size() << endl; return false;
        } 
	// everything seems alright, then proceed to save the values
	if (v.size()==1)
  {
	  // if only one value is given i.e. [ num ], consider it a constant value if doConstrFit = false. If doConstrFit = true, the parameter range will be +- 100% of the value
    if (doConstrFit) parIni[col->first] = Form("%s[ %.6f, %.6f, %.6f]", col->first.c_str(), v.at(0), v.at(0) - abs(v.at(0)), v.at(0) + abs(v.at(0)));
	  else parIni[col->first] = Form("%s[ %.6f, %.6f, %.6f]", col->first.c_str(), v.at(0), v.at(0), v.at(0));
	}
  else
  {
	  parIni[col->first] = col->first + col->second;
	}
      } else {
        parIni[col->first] = "";
      }
    }
  }

  return true;
};


bool parseString(string input, string delimiter, vector<double>& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    } else {
      cout << "[ERROR] The conversion from string to double failed!"; return false;
    }
    input = end; 
    if(input.find(delimiter.c_str())!= std::string::npos){ input.erase(input.find(delimiter.c_str()), delimiter.length()); }
  }
  return true;
};


bool parseFile(string FileName, vector< map<string, string> >& data)
{
  vector< vector<string> > content, tmp; 
  if(!readFile(FileName, tmp, -1, 1)){ return false; }
  vector<string> header = tmp.at(0);
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
  for(vector<string>::iterator rHeader=header.begin(); rHeader!=header.end(); ++rHeader) {
    if (*rHeader=="") { cout << "[ERROR] A column has no label!" << endl; return false; }
  }

  for(vector< vector<string> >::iterator row=content.begin()+1; row!=content.end(); ++row) {
    map<string, string> col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row->size()) {
	col[header.at(i)] = row->at(i);
      } else {
	col[header.at(i)] = "";
      }
    }
    data.push_back(col);
  }

  return true;
};					


bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection)
{
  vector< vector<string> > content; 
  if(!readFile(InputTrees, content, 2)){ return false; }
  for(vector< vector<string> >::iterator row=content.begin(); row!=content.end(); ++row) {
    // remove spaces
    row->at(0).erase(std::remove(row->at(0).begin(), row->at(0).end(), ' '), row->at(0).end());
    row->at(1).erase(std::remove(row->at(1).begin(), row->at(1).end(), ' '), row->at(1).end());
    // remove tabs
    row->at(0).erase(std::remove(row->at(0).begin(), row->at(0).end(), '\t'), row->at(0).end());
    row->at(1).erase(std::remove(row->at(1).begin(), row->at(1).end(), '\t'), row->at(1).end());
    if (row->at(0)!="" && row->at(1)=="") { cout << "[ERROR] There is an empty file name in your InputTrees.txt, please fix it" << endl; return false; }
    if (row->at(0)=="" && row->at(1)!="") { cout << "[ERROR] There is an empty file tag in your InputTrees.txt, please fix it" << endl; return false; }
    if (row->at(0)!="" && row->at(1)!="") {
      // store the filenames mapped by the tag
      InputFileCollection[row->at(0)].push_back(row->at(1));
    }
  }
  return true;
};


bool readFile(string FileName, vector< vector<string> >& content, const int nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Reading file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      if (nRow==0) break; else {nRow=nRow-1;}
      stringstream row(line);
      vector<string> cols; int i=0;
      while (true){
	string col; getline(row, col, ';');
	if ( (nCol>=0) ? (i>=nCol) : (col=="") ){ break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  } else {
    cout << "[ERROR] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};


bool iniWorkEnv( map< string, vector<string> >& DIR, const string workDirName)
{
  cout << "[INFO] Initializing the work enviroment" << endl;
  DIR["main"].push_back(gSystem->ExpandPathName(gSystem->pwd()));
  DIR["macros"].push_back(DIR["main"][0] + "/Macros/");
  if (existDir(DIR["macros"][0].c_str())==false){ 
    cout << "[ERROR] Input directory: " << DIR["macros"][0] << " does not exist!" << endl; 
    return false; 
  }
  DIR["input"].push_back(DIR["main"][0] + "/Input/" + workDirName + "/");
  if (existDir(DIR["input"][0])==false){ 
    cout << "[ERROR] Input directory: " << DIR["input"][0] << " does not exist!" << endl; 
    return false; 
  } else {
    findSubDir(DIR["input"], DIR["input"][0]);
  }
  DIR["output"].push_back(DIR["main"][0] + "/Output/" + workDirName + "/");
  if (existDir(DIR["output"][0].c_str())==false){ 
    cout << "[INFO] Output directory: " << DIR["output"][0] << " does not exist, will create it!" << endl;  
    gSystem->mkdir(DIR["output"][0].c_str(), kTRUE);
  }
  for(uint j = 1; j < DIR["input"].size(); j++) {
    string subdir = DIR["input"][j];
    subdir.replace(subdir.find("/Input/"), std::string("/Input/").length(), "/Output/");
    if (existDir(subdir.c_str())==false){  
      gSystem->mkdir(subdir.c_str(), kTRUE);
      cout << "[INFO] Output subdirectory: " << subdir << " created!" << endl;
    }
    DIR["output"].push_back(subdir);
  } 
  DIR["dataset"].push_back(DIR["main"][0] + "/DataSet/");
  if (existDir(DIR["dataset"][0])==false){ 
    cout << "[INFO] DataSet directory: " << DIR["dataset"][0] << " does not exist, will create it!" << endl;  
    gSystem->mkdir(DIR["dataset"][0].c_str(), kTRUE);
  }
  return true;
};


void findSubDir(vector<string>& dirlist, string dirname)
{
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  TList *subdirs = dir.GetListOfFiles();
  if (subdirs) {
    TSystemFile *subdir;
    TIter next(subdirs);
    while ((subdir=(TSystemFile*)next())) {
      if (subdir->IsDirectory() && string(subdir->GetName())!="." && string(subdir->GetName())!="..") {
        dirlist.push_back(dirname + subdir->GetName() + "/");
        cout << "[INFO] Input subdirectory: " << dirname + subdir->GetName() + "/" << " found!" << endl;
      }
    }
  }
  delete subdirs;
  return;
};


bool existDir(string dir)
{
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp){
    gSystem->FreeDirectory(dirp);
    exist = true;
  }
  return exist;
};


bool checkSettings(bool fitData, bool fitMC, bool fitPbPb, bool fitPP, bool fitMass, bool fitCtau, bool fitCtauTrue, bool fitCtauReco, bool fitRes, bool doCtauErrPDF, bool incJpsi, bool incPsi2S, bool incBkg, bool incPrompt, bool incNonPrompt, bool cutCtau, bool doConstrFit, bool doSimulFit, bool wantPureSMC, const char* applyCorr, bool setLogScale, bool zoomPsi, bool incSS, int numCores)
{ 
  cout << "[INFO] Checking user settings " << endl;

  if (!fitMass && !fitCtau && !fitCtauTrue && !doCtauErrPDF && !fitCtauReco && !fitRes) {
    cout << "[ERROR] At least one distribution has to be selected for fitting, please select either Mass, CtauTrue, CtauErr, Ctau, CtauN or CtauReco!" << endl; return false;
  }
  if (fitCtau && fitRes) {
    cout << "[ERROR] We can not fit Ctau and CtauRes at the same time, please fix your input settings!" << endl; return false;
  }
  if (doConstrFit && (fitCtau || fitCtauTrue || fitCtauReco|| fitRes)) {
    cout << "[ERROR] Constrained fits only implemented for invariant mass!" << endl; return false;
  }
  if (doConstrFit && doSimulFit) {
    cout << "[ERROR] Constrained fits not implemented with simultaneous fits!" << endl; return false;
  }
  if (fitCtauTrue && fitData) {
    cout << "[ERROR] We can not fit the truth ctau distribution on data, please select only MC!" << endl; return false;
  }
  if (fitCtauTrue && !incNonPrompt) {
    cout << "[ERROR] The ctau truth distribution is only fitted on Non-Prompt MC, please include the Non-Prompt components!" << endl; return false;
  }
  if (!fitData && !fitMC) {
    cout << "[ERROR] At least one type of dataset has to be selected for fitting, please select either Data or MC!" << endl; return false;
  }
  if (!fitMass && doSimulFit) {
    cout << "[ERROR] Can't perform simultanoeus fit without fitting the mass spectrum, please fix your fitter settings!" << endl; return false;
  }
  if (fitCtau && (!incPrompt && !incNonPrompt)) {
    cout << "[ERROR] No ctau models were included, please include either Prompt or NonPrompt when fitting ctau!" << endl; return false;
  }
  if (fitMC && (!incPrompt && !incNonPrompt)) {
    cout << "[ERROR] No MC datasets were included, please include either Prompt or Non-Prompt MC!" << endl; return false;
  }
  if (!fitPbPb && !fitPP) {
    cout << "[ERROR] No datasets were selected, please select either PbPb or PP!" << endl; return false;
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] No models were included for fitting, please include either Jpsi, Psi2S or Bkg" << endl; return false;
  }
  if ((!fitData || fitMC) && (!incJpsi && !incPsi2S)) {
    cout << "[ERROR] We can't fit only background on MC, please include also either Jpsi or Psi2S!" << endl; return false;
  }
  if (fitData && wantPureSMC) {
    cout << "[WARNING] Pure signal MC will not be fitted since you have selected DATA!" << endl;
  }
  if (!fitData && strcmp(applyCorr,"")) {
    cout << "[WARNING] Correction only to Data, please make fitData = true" << endl; return false;
  }
  if (wantPureSMC && (!incPsi2S && !incJpsi)) {
    cout << "[ERROR] Pure signal MC can only be fitted if at least one signal model is included, please fix your input settings!" << endl;
  }
  if (numCores<=0 || numCores>32) {
    cout << "[ERROR] Invalid number of cores: " << numCores << " , please select an interger number between 1 and 32." << endl; return false;
  }
  if (doSimulFit && (!fitPbPb || !fitPP)) {
    cout << "[ERROR] Can't perform simultaneous fit without both PP and PbPb, please fix your fitter settings!" << endl; return false;
  }
  if (doSimulFit && (!incPsi2S || !incJpsi)) {
    cout << "[ERROR] Simultaneous fits require both Psi2S and Jpsi, please fix your fitter settings!" << endl; return false;
  }
  if (doSimulFit && !fitData) {
    cout << "[ERROR] Can't perform simultaneous fits on MC, please fix your fitter settings!" << endl; return false;
  }
  if (setLogScale && zoomPsi) {
    cout << "[ERROR] Zooming on Psi(2S) is currently not supported with logarithmic scale." << endl;
  }

  cout << "[INFO] All user setting are correct " << endl;
  return true;
};

