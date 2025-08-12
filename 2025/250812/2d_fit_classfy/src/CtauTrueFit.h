#pragma once
#include "TString.h"

class TFile;
class RooFitResult;
class RooWorkspace;

class CtauTrueFit {
public:
  CtauTrueFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE);
  ~CtauTrueFit();
  int nCtauBins = 200; //

  void init();
  void run();

  TString DATE = "default";
  TString inputFilePath;

  bool isWeighted = false;
  int nExp = 2;

  double ctauMin, ctauMax;
  bool isLogon = true;

  // Todo: user custom in script
  Double_t ctau3DMin = 0.5;
  Double_t ctau3DMax = 7;
  Double_t ctau3DHigh = 7;

  int nBins = 100; // plotting bins

private:
  void setLabels();
  void makeOutputFolder();
  void turnOffRooFitMessage();
  void openInput();
  void processDataset();
  void setVariableRanges();
  void fixParameters();
  void buildCtauTrueModel();
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
  TFile *fInputTrue = nullptr;

  // --- fit result
  RooFitResult *fitResult = nullptr;

  // --- internal labels ---
  TString kineLabel;
  TString bCont; // pr, np tag

  // --- workspace ---
  RooWorkspace *ws = nullptr;
};

