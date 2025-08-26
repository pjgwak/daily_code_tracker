#pragma once
#include "TString.h"

class TFile;
class RooFitResult;
class RooWorkspace;

class McMassFit {
public:
  McMassFit(float ptLow, float ptHigh,
          float yLow, float yHigh,
          int cLow, int cHigh,
          float cosLow, float cosHigh,
          int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE);
  ~McMassFit();

  void init();
  void run();

  TString DATE = "default";
  TString inputFilePath;

  bool isWeighted = false;

  double fitLimit = 3.21;

  bool isLogon = true;

  void initVar(const std::string &varName, double init, double low, double high);

private:
  void setLabels();
  void makeOutputFolder();
  void turnOffRooFitMessage();
  void openInput();
  void processDataset();
  void setVariableRanges();
  void buildModel();
  void doFit();
  void drawPlot();
  void saveOutput();

  // --- kinematic bins ---
  float ptLow, ptHigh;
  float yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR = 1;
  int PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  // --- input files ---
  TFile *fInputMc = nullptr;

  // --- fit result
  RooFitResult *fitResult = nullptr;

  // --- internal labels ---
  TString kineLabel;
  TString bCont; // pr, np tag

  // --- workspace ---
  RooWorkspace *ws = nullptr;
};