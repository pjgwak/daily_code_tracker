#pragma once
#include "TString.h"

class Final2DFit {
public:
  Final2DFit(float ptLow, float ptHigh,
             float yLow, float yHigh,
             int cLow, int cHigh,
             float cosLow, float cosHigh,
             int PR, int PRw, bool fEffW, bool fAccW, bool isPtW, bool isTnP);
  ~Final2DFit();
  int nCtauBins = 200; //

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
};