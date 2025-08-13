#pragma once
#include "TString.h"

class TFile;
class RooFitResult;
class RooWorkspace;

class CtauResFit {
public:
  CtauResFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE);
  ~CtauResFit();
  int nCtauBins = 200; //

  void init();
  void run();

  TString DATE = "default";
  TString inputFilePath;

  bool isWeighted = false;
  int nGauss = 2;

  double ctauMin, ctauMax;
  bool isLogon = true;

  double ctauResMin = -10;
  double ctauResMax = 0;

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
  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;

  // --- fit result
  RooFitResult *fitResult = nullptr;

  // --- internal labels ---
  TString kineLabel;
  TString bCont; // pr, np tag

  // --- workspace ---
  RooWorkspace *ws = nullptr;
};

