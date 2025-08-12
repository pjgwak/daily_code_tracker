#pragma once
#include "TString.h"

class TFile;
class RooFitResult;
class RooWorkspace;

class CtauBkgFit {
public:
  CtauBkgFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE);
  ~CtauBkgFit();
  int nCtauBins = 200; //

  void init();
  void run();

  TString DATE = "default";
  TString inputFilePath;

  bool isWeighted = false;
  int nGauss, nExpL, nExpR, nExpC;

  double ctauMin, ctauMax;
  bool isLogon = true;

private:
  void setLabels();
  void makeOutputFolder();
  void turnOffRooFitMessage();
  void openInput();
  void processDataset();
  void setVariableRanges();
  void fixParameters();
  void buildCtauFitModel();
  void buildCtauCondModel();
  void doFit();
  void drawCtauPull();
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
  TFile *fCRes = nullptr;

  // --- fit result
  RooFitResult *fitResult = nullptr;

  // --- internal labels ---
  TString kineLabel;
  TString fname; // pr, np tag

  // --- workspace ---
  RooWorkspace *ws = nullptr;
};