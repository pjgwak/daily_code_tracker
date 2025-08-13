#pragma once
#include "TString.h"

class TFile;
class RooFitResult;
class RooWorkspace;

class Final2DFit
{
public:
  Final2DFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP, TString DATE);
  ~Final2DFit();
  int nCtauBins = 200; //

  void init();
  void run();

  TString DATE = "default";
  TString inputFilePath;

  bool isWeighted = false;
  double nGauss, nExpTrue; //Todo: 이것도 trueFit에서 만들고 가져오는 게 더 나을지도?
  bool isLogon = true;

  void initVar(const std::string &varName, double init, double low, double high);

private:
  void setLabels();
  void makeOutputFolder();
  void turnOffRooFitMessage();
  void openInput();
  void processDataset();
  void setVariableRanges();
  void fixParameters();
  void buildCtauNpModel();
  void buildMassCtauModel();
  void do2dFit();
  void drawCtauPull();
  void drawMass();
  void drawCtauRatio();
  void saveOutput();

  // --- kinematic bins ---
  float ptLow, ptHigh;
  float yLow, yHigh;
  int cLow, cHigh;
  float cosLow, cosHigh;
  int PR=1;
  int PRw;
  bool fEffW, fAccW, isPtW, isTnP;

  // --- input files ---
  TFile *fInputData = nullptr;
  TFile *fMass = nullptr;
  TFile *fCErr = nullptr;
  TFile *fCRes = nullptr;
  TFile *fCBkg = nullptr;
  TFile *fCTrue = nullptr;
  TFile *fMcParams = nullptr;

  // --- fit result
  RooFitResult *fitResult = nullptr;

  // --- internal labels ---
  TString kineLabel;
  TString fname; // pr, np tag

  // --- workspace ---
  RooWorkspace *ws = nullptr;
};