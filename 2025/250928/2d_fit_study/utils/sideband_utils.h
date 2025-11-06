#pragma once

#include "mass_fit_utils.h"

struct FitConfigSB
{
  //  ctauErr distribution
  std::string mode;
  std::string system;
  bool isPrompt;

  float ptLow, ptHigh;
  float yLow, yHigh;
  float massLow, massHigh;
  float cLow, cHigh;
  float ctErrLow, ctErrHigh;

  std::string inputFile;
  std::string figDir;
  std::string rootDir;

  std::string rootlogon;
  int numCPU;

  // parameters
  std::string massBkgPdf;
  par3 mu0, sigma1, r21, r32;
  par3 fg1, fg2;

	float sblLow, sblHigh;
	float sigLow, sigHigh;
	float sbrLow, sbrHigh;
};