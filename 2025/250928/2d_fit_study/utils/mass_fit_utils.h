#pragma once

#include <string>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>

struct par3 
{
  double init, min, max;
}; // {init, min, max}

struct FitConfigMcMass
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
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
};

struct FitConfigMass
{
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
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