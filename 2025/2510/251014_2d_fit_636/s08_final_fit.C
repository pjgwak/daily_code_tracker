// File: fit2D_with_punizi.C
// root -l -q 'fit2D_with_punizi.C("input.root","tree","mass","ctau3D", 2.6,3.5, -0.1,0.2)'

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooProdPdf.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
using namespace RooFit;

void s08_final_fit(const char *fname = "root.root",
                   const char *tname = "tree",
                   const char *massVar = "mass",
                   const char *ctauVar = "ctau3D",
                   double mmin = 2.6, double mmax = 3.5,
                   double tmin = -0.10, double tmax = 0.20,
                   int nb_m = 80, int nb_t = 80)
{
  // === 데이터 로드 ===
  TFile fin(fname, "READ");
  if (fin.IsZombie())
  {
    ::Error("fit2D", "Cannot open file");
    return;
  }
  TTree *tr = (TTree *)fin.Get(tname);
  if (!tr)
  {
    ::Error("fit2D", "Tree not found");
    return;
  }

  RooRealVar m(massVar, massVar, mmin, mmax);
  RooRealVar t(ctauVar, ctauVar, tmin, tmax);
  RooArgSet obs(m, t);
  RooDataSet data("data", "data", tr, obs); // 전부 읽기

  // === Mass: Signal (DCB + Gauss) ===
  RooRealVar mean("mean", "mean", 3.096, mmin, mmax);
  RooRealVar sigma("sigma", "CB sigma", 0.030, 0.005, 0.100);
  RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.2, 5.0);
  RooRealVar nL("nL", "nL", 3.0, 1.0, 20.0);
  RooRealVar alphaR("alphaR", "alphaR", -1.5, -5.0, -0.2);
  RooRealVar nR("nR", "nR", 3.0, 1.0, 20.0);
  RooRealVar sigmaG("sigmaG", "Gauss sigma", 0.020, 0.002, 0.100);
  RooRealVar fRight("fRight", "frac(CBR in DCB)", 0.50, 0.00, 1.00);
  RooRealVar fG("fG", "frac(Gauss in Sig)", 0.20, 0.00, 1.00);

  RooCBShape CBL("CBL", "left-tail CB", m, mean, sigma, alphaL, nL);
  RooCBShape CBR("CBR", "right-tail CB", m, mean, sigma, alphaR, nR);
  RooAddPdf DCB("DCB", "Double-CB approx", RooArgList(CBR, CBL), RooArgList(fRight), true);
  RooGaussian mG("mG", "core Gauss", m, mean, sigmaG);
  RooAddPdf MassSig("MassSig", "DCB+G", RooArgList(mG, DCB), RooArgList(fG), true);

  // === Mass: Background (Cheby 2) ===
  RooRealVar sl1("sl1", "Cheb c1", 0.0, -1.0, 1.0);
  RooRealVar sl2("sl2", "Cheb c2", 0.0, -1.0, 1.0);
  RooChebychev MassBkg("MassBkg", "cheb2", m, RooArgList(sl1, sl2));

  // === ctau: Signal (GaussRes ⊗ RooDecay SingleSided) ===
  RooRealVar muResS("muResS", "res mean S", 0.0);
  muResS.setConstant(true);
  RooRealVar sigResS("sigResS", "res sigma S", 0.020, 0.002, 0.080);
  RooGaussModel tResS("tResS", "res S", t, muResS, sigResS);
  RooRealVar tauS("tauS", "tauS", 0.030, 0.001, 0.300);
  RooDecay CtSig("CtSig", "exp⊗gauss S", t, tauS, tResS, RooDecay::SingleSided);

  // === ctau: Background (Left/Flipped, Mid/DoubleSided, Right/SingleSided) + 공통 res ===
  RooRealVar muResB("muResB", "res mean B", 0.0);
  muResB.setConstant(true);
  RooRealVar sigResB("sigResB", "res sigma B", 0.030, 0.003, 0.120);
  RooGaussModel tResB("tResB", "res B", t, muResB, sigResB);

  RooRealVar tauL("tauL", "tau left", 0.040, 0.002, 0.500);
  RooRealVar tauM("tauM", "tau mid", 0.060, 0.002, 0.800);
  RooRealVar tauR("tauR", "tau right", 0.040, 0.002, 0.500);

  RooDecay CtLeft("CtLeft", "left (flipped)", t, tauL, tResB, RooDecay::Flipped);
  RooDecay CtMid("CtMid", "double-sided mid", t, tauM, tResB, RooDecay::DoubleSided);
  RooDecay CtRight("CtRight", "right (single)", t, tauR, tResB, RooDecay::SingleSided);

  // 혼합비: fL, fM (fR = 1 - fL - fM)
  RooRealVar fL("fL", "frac left", 0.20, 0.00, 1.00);
  RooRealVar fM("fM", "frac mid", 0.50, 0.00, 1.00);
  RooFormulaVar fR("fR", "1.0 - fL - fM", RooArgList(fL, fM)); // 자동으로 1-fL-fM

  RooAddPdf CtBkg("CtBkg", "L/M/R mix",
                  RooArgList(CtRight, CtMid, CtLeft),
                  RooArgList(fR, fM)); // 3성분: fR, fM, (1-fR-fM)

  // === 2D 팩터라이즈 == //
  RooProdPdf Sig2D("Sig2D", "MassSig × CtSig", RooArgSet(MassSig, CtSig));
  RooProdPdf Bkg2D("Bkg2D", "MassBkg × CtBkg", RooArgSet(MassBkg, CtBkg));

  RooRealVar Nsig("Nsig", "Nsig", data.numEntries() * 0.6, 0.0, 1e9);
  RooRealVar Nbkg("Nbkg", "Nbkg", data.numEntries() * 0.4, 0.0, 1e9);
  RooAddPdf Tot2D("Tot2D", "Sig+Bkg", RooArgList(Sig2D, Bkg2D), RooArgList(Nsig, Nbkg))
}
