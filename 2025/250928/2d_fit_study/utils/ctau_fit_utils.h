#pragma once

#include "mass_fit_utils.h"
#include <string>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>

struct FitConfigCtauRes
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
  float ctResLow, ctResHigh;

  std::string inputFile;
  std::string figDir;
  std::string rootDir;

  std::string rootlogon;
  int numCPU;

  // parameters
  std::string ctauResPdf;
  par3 mu0, sigma1, r21, r32;
  par3 fg1, fg2;
};

struct FitConfigCtauPrMc
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
  float ctLow, ctHigh;

  std::string inputFile;
  std::string figDir;
  std::string rootDir;

  std::string rootlogon;
  int numCPU;

  // parameters
  std::string ctauResPdf;
  par3 mu0, sigma1, r21, r32;
  par3 fg1, fg23;
};

struct FitConfigCtauBkg
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
  float ctLow, ctHigh;

  std::string inputFile;
  std::string figDir;
  std::string rootDir;

  std::string rootlogon;
  int numCPU;

  // parameters
  std::string ctauResPdf;
  par3 mu0, sigma1, r21, r32;
  par3 fg1, fg23;

  std::string ctauBkgPdf;
  par3 fdecayM, fdecayLR;
  par3 tauMid, tauRatio, tauL;

  std::vector<std::string> fixNames;
};

struct FitConfigCtauTrue
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
  float ctLow, ctHigh;
  std::string inputFile;
  std::string figDir;
  std::string rootDir;

  std::string rootlogon;
  int numCPU;

  // parameters
  std::string massSigPdf;
  par3 mean, sigmaL, alphaL, nL;
  par3 sigmaRatio, alphaRatio, nRatio, sigmaGRatio;
  par3 f_cbR, f_gaus;

  std::string massBkgPdf;
  par3 tauMass; // exponential
  par3 s1, s2, s3; // cheby1 ~6
  par3 s4, s5, s6; 

  // yield
  par3 nSig, nBkg;

  // fixed parameters - sigmaL and mean should be free
  std::vector<std::string> fixNames;
};