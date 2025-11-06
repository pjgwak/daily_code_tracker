#pragma once

#include "../utils/mass_fit_utils.h"
#include "../utils/ctau_fit_utils.h"
#include "../utils/sideband_utils.h"

FitConfigMcMass getConfigMcMass() {
  FitConfigMcMass cfg;

  cfg.mode = "MC";
  cfg.system = "pp";
  cfg.isPrompt = true;

  cfg.inputFile = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";

  cfg.figDir = "figs/pp_run2_pT6p5_9_y0_1p6";
  cfg.rootDir = "roots/pp_run2_pT6p5_9_y0_1p6";
  
  cfg.ptLow = 6.5; cfg.ptHigh = 9.0;
  cfg.yLow = 0; cfg.yHigh = 1.6;
  cfg.cLow = 0; cfg.cHigh = 180;

  cfg.rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  cfg.numCPU = 32;
  cfg.massLow = 2.6; cfg.massHigh = 3.5;


  // fit parameters
  cfg.massSigPdf = "DCBGauss";
  cfg.mean = {3.0969, 3.0469, 3.1469};
  cfg.sigmaL = {0.025, 0.001, 0.050};
  cfg.alphaL = {0.60, 0.20, 2.00};
  cfg.nL = {3.00, 1.00, 5.00};

  cfg.sigmaRatio = {1.00, 0.50, 5.00};
  cfg.alphaRatio = {3.60, 0.50, 10.0};
  cfg.nRatio = {0.60, 0.10, 5.00};

  cfg.sigmaGRatio = {1.00, 0.50, 5.00};

  cfg.f_cbR = {0.50, 0.00, 1.00};
  cfg.f_gaus = {0.05, 0.00, 1.00};

  return cfg;
}

FitConfigMass getConfigMass()
{
  FitConfigMass cfg;

  cfg.mode = "Data";
  cfg.system = "pp";
  cfg.isPrompt = false;

  cfg.inputFile = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";

  cfg.figDir = "figs/pp_run2_pT6p5_9_y0_1p6";
  cfg.rootDir = "roots/pp_run2_pT6p5_9_y0_1p6";

  cfg.ptLow = 6.5;
  cfg.ptHigh = 9.0;
  cfg.yLow = 0;
  cfg.yHigh = 1.6;
  cfg.cLow = 0;
  cfg.cHigh = 180;

  cfg.rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  cfg.numCPU = 32;
  cfg.massLow = 2.6;
  cfg.massHigh = 3.5;

  // fit parameters
  cfg.massSigPdf = "DCBGauss";
  cfg.mean = {3.0969, 3.0469, 3.1469};
  cfg.sigmaL = {0.025, 0.001, 0.050};
  cfg.f_cbR = {0.50, 0.00, 1.00};

  // --- fixed parameters ---
  // sigmaL and mean should be free
  cfg.fixNames = {"alphaL", "alphaRatio", "f_gaus", "nL", "nRatio", "sigmaGRatio", "sigmaRatio"};
  
  // these are need for initialization
  // but will be overwritten by MC result - don't touch
  cfg.alphaL = {0.60, 0.20, 2.00};
  cfg.nL = {3.00, 1.00, 5.00};
  cfg.sigmaRatio = {1.00, 0.50, 5.00};
  cfg.alphaRatio = {3.60, 0.50, 10.0};
  cfg.nRatio = {0.60, 0.10, 5.00};
  cfg.sigmaGRatio = {1.00, 0.50, 5.00};
  cfg.f_gaus = {0.05, 0.00, 1.00};
  // ------------------------------------------

  std::string massBkgPdf;
  cfg.massBkgPdf = "expo";
  cfg.tauMass = {0.1, -5, 1};
  cfg.s1 = {0.01, -1, 1};
  cfg.s2 = {0.01, -1, 1};

  // yields
  cfg.nSig = {10000, 1, 1000000};
  cfg.nBkg = {10000, 1, 300000};

  return cfg;
}

FitConfigSB getConfigSB()
{
  FitConfigSB cfg;

  cfg.mode = "Data";
  cfg.system = "pp";
  cfg.isPrompt = true;

  cfg.inputFile = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC0_JPsi_pp_y0.00_2.40_Effw0_Accw0_PtW1_TnP1_230221.root";

  cfg.figDir = "figs/pp_run2_pT6p5_9_y0_1p6";
  cfg.rootDir = "roots/pp_run2_pT6p5_9_y0_1p6";

  cfg.ptLow = 6.5;
  cfg.ptHigh = 9.0;
  cfg.yLow = 0;
  cfg.yHigh = 1.6;
  cfg.cLow = 0;
  cfg.cHigh = 180;

  cfg.ctErrLow = 0.0001;
  cfg.ctErrHigh = 0.3;

  cfg.sblLow = 2.6; cfg.sblHigh = 2.9;
	cfg.sigLow = 2.9; cfg.sigHigh = 3.3;
	cfg.sbrLow = 3.3; cfg.sbrHigh = 3.5;

  cfg.rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  cfg.numCPU = 32;
  cfg.massLow = 2.6;
  cfg.massHigh = 3.5;

  // fit parameters
  cfg.massBkgPdf = "expo";

  cfg.mu0 = {0.009, -0.02, 0.02};
  cfg.sigma1 = {0.5, 0.01, 5};
  cfg.r21 = {1.9, 1, 10};
  cfg.r32 = {2.9, 1, 10};

  cfg.fg1 = {0.72, 0, 1};
  cfg.fg2 = {0.26, 0, 0.4};

  return cfg;
}

FitConfigCtauRes getConfigCtauRes()
{
  FitConfigCtauRes cfg;

  cfg.mode = "MC";
  cfg.system = "pp";
  cfg.isPrompt = true;

  cfg.inputFile = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";

  cfg.figDir = "figs/pp_run2_pT6p5_9_y0_1p6";
  cfg.rootDir = "roots/pp_run2_pT6p5_9_y0_1p6";

  cfg.ptLow = 6.5;
  cfg.ptHigh = 9.0;
  cfg.yLow = 0;
  cfg.yHigh = 1.6;
  cfg.cLow = 0;
  cfg.cHigh = 180;

  cfg.ctResLow = -10;
  cfg.ctResHigh = 10;

  cfg.rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  cfg.numCPU = 32;
  cfg.massLow = 2.6;
  cfg.massHigh = 3.5;

  // fit parameters
  cfg.ctauResPdf = "Gauss3";

  cfg.mu0 = {0.009, -0.02, 0.02};
  cfg.sigma1 = {0.5, 0.01, 5};
  cfg.r21 = {1.9, 1, 10};
  cfg.r32 = {2.9, 1, 10};

  cfg.fg1 = {0.72, 0, 1};
  cfg.fg2 = {0.26, 0, 0.4};

  return cfg;
}

FitConfigCtauPrMc getConfigCtauPrMc()
{
  FitConfigCtauPrMc cfg;

  cfg.mode = "MC";
  cfg.system = "pp";
  cfg.isPrompt = true;

  cfg.inputFile = "/data/hwan/psi2S_RAA_PbPb2018/skimmedFiles/OniaRooDataSet_isMC1_Jpsi_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230717.root";

  cfg.figDir = "figs/pp_run2_pT6p5_9_y0_1p6";
  cfg.rootDir = "roots/pp_run2_pT6p5_9_y0_1p6";

  cfg.ptLow = 6.5;
  cfg.ptHigh = 9.0;
  cfg.yLow = 0;
  cfg.yHigh = 1.6;
  cfg.cLow = 0;
  cfg.cHigh = 180;

  cfg.ctResLow = -1;
  cfg.ctResHigh = 1;

  cfg.rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  cfg.numCPU = 32;
  cfg.massLow = 2.6;
  cfg.massHigh = 3.5;

  // fit parameters
  cfg.ctauResPdf = "Gauss2";

  cfg.mu0 = {0.009, -0.02, 0.02};
  cfg.sigma1 = {0.02, 0.01, 0.1};
  cfg.r21 = {1.9, 1, 50};
  cfg.r32 = {2.9, 1, 50};

  cfg.fg1 = {0.72, 0, 1};
  cfg.fg2 = {0.26, 0, 0.4};

  return cfg;
}
