// Ref: https://root.cern/doc/v632/rf307__fullpereventerrors_8C.html
// https://root.cern.ch/doc/v636/rf303__conditional_8C_source.html

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TFile.h"
#include "TLine.h"
#include "TROOT.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "TSystem.h"
#include "RooChebychev.h"
#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooAddModel.h"
#include "RooHist.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "RooResolutionModel.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "TParameter.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <cstdlib>

using namespace RooFit;

// ROOT::EnableImplicitMT(); // Don't use it in CNU server.

// ---------- chi2 from pull points ----------
static std::pair<double, int> chi2_from_pull(RooHist *hp)
{
  double chi2 = 0.0;
  int n = 0;
  if (!hp)
    return {0.0, 0};
  for (int i = 0; i < hp->GetN(); ++i)
  {
    double x, y;
    hp->GetPoint(i, x, y);
    if (std::isfinite(y))
    {
      chi2 += y * y;
      ++n;
    }
  }
  return {chi2, n};
}

static int count_float_pars(const RooAbsPdf &pdf, const RooArgSet &obs, const RooFitResult *fr)
{
  RooArgSet params;
  pdf.getParameters(&obs, params);
  if (!fr)
    return params.getSize();
  std::unique_ptr<RooAbsCollection> common(params.selectCommon(fr->floatParsFinal()));
  return common ? common->getSize() : 0;
}

static int count_float_pars(const RooArgSet &params, const RooFitResult *fr)
{
  if (!fr)
    return params.getSize();
  std::unique_ptr<RooAbsCollection> common(params.selectCommon(fr->floatParsFinal()));
  return common ? common->getSize() : 0;
}

static void apply_err_param(RooRealVar &dst, const RooArgSet &srcPars)
{
  auto *src = dynamic_cast<RooRealVar *>(srcPars.find(dst.GetName()));
  if (!src)
  {
    printf("[WARN] errPdf param '%s' missing in file\n", dst.GetName());
    return;
  }
  dst.setVal(src->getVal());
  dst.setError(src->getError());
  dst.setConstant(true);
}

static bool load_errpdf_params(const TString &filePath,
                               RooArgSet &targetParsBkg, RooArgSet &targetParsSig,
                               RooFitResult *&fitErrBkgPtr, RooFitResult *&fitErrSigPtr)
{
  if (gSystem->AccessPathName(filePath))
  {
    printf("[ERROR] errPdf file not found: %s\n", filePath.Data());
    return false;
  }
  TFile f(filePath, "READ");
  auto *parsBkg = dynamic_cast<RooArgSet *>(f.Get("errParsBkg"));
  auto *parsSig = dynamic_cast<RooArgSet *>(f.Get("errParsSig"));
  auto *fitErrBkg = dynamic_cast<RooFitResult *>(f.Get("fitErrBkg"));
  auto *fitErrSig = dynamic_cast<RooFitResult *>(f.Get("fitErrSig"));
  if (!parsBkg || !parsSig)
  {
    printf("[ERROR] errPdf params missing in %s\n", filePath.Data());
    return false;
  }
  for (auto obj : targetParsBkg)
    if (auto *arg = dynamic_cast<RooRealVar *>(obj))
      apply_err_param(*arg, *parsBkg);
  for (auto obj : targetParsSig)
    if (auto *arg = dynamic_cast<RooRealVar *>(obj))
      apply_err_param(*arg, *parsSig);
  fitErrBkgPtr = nullptr;
  fitErrSigPtr = nullptr;
  printf("[OK] loaded errPdf params from %s\n", filePath.Data());
  return true;
}

namespace err_range
{
struct RangeStats
{
  std::size_t n = 0;
  double rawMin = 0.0;
  double rawMax = 0.0;
};

struct Ctau3DRangeParams
{
  double qLow = 0.001;
  double qHigh = 0.999;
  bool clampNonNegative = false;
};

struct ErrRangeParams
{
  double qLow = 0.001;
  double qHigh = 0.995;
  bool clampNonNegative = true;
};

constexpr Ctau3DRangeParams kCtau3D_Mass{0.001, 0.999, false};
constexpr ErrRangeParams kErr_Default{0.001, 0.995, true};

template <typename Config>
static bool compute_quantile_range(RooDataSet &ds,
                                   RooRealVar &var,
                                   const Config &cfg,
                                   double &outLow,
                                   double &outHigh,
                                   RangeStats &stats)
{
  std::vector<double> vals;
  vals.reserve(ds.numEntries());
  stats.n = 0;
  stats.rawMin = 0.0;
  stats.rawMax = 0.0;
  bool first = true;
  for (int i = 0; i < ds.numEntries(); ++i)
  {
    ds.get(i);
    const double v = var.getVal();
    if (!std::isfinite(v))
      continue;
    if (cfg.clampNonNegative && v < 0.0)
      continue;
    vals.push_back(v);
    if (first)
    {
      stats.rawMin = v;
      stats.rawMax = v;
      first = false;
    }
    else
    {
      stats.rawMin = std::min(stats.rawMin, v);
      stats.rawMax = std::max(stats.rawMax, v);
    }
  }
  stats.n = vals.size();
  if (vals.empty())
    return false;
  std::sort(vals.begin(), vals.end());
  auto qAt = [&](double q)
  {
    long long idx = llround(q * (vals.size() - 1));
    if (idx < 0)
      idx = 0;
    if (idx >= static_cast<long long>(vals.size()))
      idx = static_cast<long long>(vals.size()) - 1;
    return vals[static_cast<size_t>(idx)];
  };
  outLow = qAt(cfg.qLow);
  outHigh = qAt(cfg.qHigh);
  return std::isfinite(outLow) && std::isfinite(outHigh) && outLow < outHigh;
}

static bool compute_ctau3d_range(RooDataSet &ds,
                                 RooRealVar &var,
                                 const Ctau3DRangeParams &cfg,
                                 double /*ptHigh*/,
                                 double &outLow,
                                 double &outHigh,
                                 RangeStats &stats)
{
  return compute_quantile_range(ds, var, cfg, outLow, outHigh, stats);
}

static bool compute_err_range(RooDataSet &ds,
                              RooRealVar &var,
                              const ErrRangeParams &cfg,
                              double &outLow,
                              double &outHigh,
                              RangeStats &stats)
{
  return compute_quantile_range(ds, var, cfg, outLow, outHigh, stats);
}
} // namespace err_range

void err(float ptLow = 12, float ptHigh = 14, float yLow = 1.6, float yHigh = 2.4, bool doErrPdfFit = true)
{
  
  TStopwatch timer;
  timer.Start();

  enum ConstraintMode
  {
    kFree = 0,
    kFix = 1,
    kConstrain = 2,
    kSet = 3
  };

  // ------------------------------------------------------------------
  // user variables
  // ------------------------------------------------------------------
  // float xLow = -6, xHigh = 6; // ctau3D
  // err-PDF component controls
  // Gauss: 0~3, Landau: 0~3, Exp: 0~1
  int nErrBkgGaussUser = 1;
  int nErrBkgLandauUser = 1;
  const int nErrBkgExpUser = 0; // dont use

  int nErrSigGaussUser = 1;
  int nErrSigLandauUser = 1;
  const int nErrSigExpUser = 0; // dont use
  const bool isPt23 = (ptLow == 2.0 && ptHigh == 3.0);

  // === custom model combination ===
  // sig
  if (ptLow==1.0&&ptHigh==2.0) {
    nErrSigGaussUser = 2;
    nErrSigLandauUser = 2;
  }
  else if (ptLow == 2.0 && ptHigh == 3.0)
  {
    nErrSigGaussUser = 1;
    nErrSigLandauUser = 2;
  }

  if (ptLow == 8.0)
  {
    nErrSigGaussUser = 2;
    nErrSigLandauUser = 2;
    nErrBkgGaussUser = 2;
    nErrBkgLandauUser = 2;
  }

  // bkg
  if (ptLow == 1.0 && ptHigh == 2.0)
  {
    nErrBkgGaussUser = 3;
    nErrBkgLandauUser = 2;
  }
  else if (isPt23)
  {
    // Keep pt 2-3 bkg errPdf compact to avoid unstable over-parameterized fits.
    nErrBkgGaussUser = 1;
    nErrBkgLandauUser = 2;
  }
  else if (ptLow == 4.0 && ptHigh == 5.0)
  {
    nErrBkgGaussUser = 2;
    nErrBkgLandauUser = 3;
  }

  // Gauss/Landau parameter-range scale:
  // 1.0 = auto, <1 tighter, >1 looser.
  double errShapeRangeScale = 2.0;
  if (const char *s = gSystem->Getenv("ERR_SHAPE_RANGE_SCALE"))
  {
    const double v = std::atof(s);
    if (std::isfinite(v) && v > 0.0)
      errShapeRangeScale = std::max(0.4, std::min(2.5, v));
  }
  printf("[INFO] ERR_SHAPE_RANGE_SCALE=%.3f\n", errShapeRangeScale);



  // constraint controls
  const ConstraintMode massConstraintMode = kConstrain;
  const double massConstraintSigmaScale = 2;
  TString baseDir = gSystem->DirName(__FILE__);
  TString outFigDir = baseDir + "/figs";
  TString outRootDir = baseDir + "/roots";

  const char *DATA_ROOT = "/data/users/pjgwak/work/daily_code_tracker/2026/02/00_OO_oniaskim/onia_to_skim/roodataset_roots/RooDataSet_miniAOD_isMC0_Jpsi_cent0_200_Effw0_Accw0_PtW0_TnP0_OO25_HLT_OxyL1SingleMuOpen_v1.root";
  const char *DSET_NAME = "dataset";
  const char *massVarName = "mass";
  const char *xVarName = "ctau3D";
  const char *ptVarName = "pt";
  const char *yVarName = "y";
  const char *xErrVarName = "ctau3DErr";
  const double M_MIN = 2.60, M_MAX = 3.5;
  auto formatTag = [](double value) { return TString::Format("%.2f", value); };
  const TString yTag = TString::Format("y%s_%s", formatTag(yLow).Data(), formatTag(yHigh).Data());
  const TString ptTag = TString::Format("pt%s_%s", formatTag(ptLow).Data(), formatTag(ptHigh).Data());
  const TString figTag = yTag + "_" + ptTag;
  const TString figDir = TString::Format("%s/%s/err", outFigDir.Data(), yTag.Data());
  const TString resultDir = TString::Format("%s/%s/err", outRootDir.Data(), yTag.Data());
  const TString errPdfFile = TString::Format("%s/errpdf_%s.root", resultDir.Data(), figTag.Data());
  auto figName = [&](const char *name)
  {
    return TString::Format("%s/%s_%s.pdf", figDir.Data(), name, figTag.Data());
  };

  // ------------------------------------------------------------------
  // miscellaneous
  // ------------------------------------------------------------------
  // be quiet please
  RooMsgService &msg = RooMsgService::instance();
  msg.setGlobalKillBelow(RooFit::WARNING);
  for (int i = 0; i < msg.numStreams(); ++i)
  {
    msg.getStream(i).removeTopic(RooFit::Tracing);
    msg.getStream(i).removeTopic(RooFit::Eval);
    msg.getStream(i).removeTopic(RooFit::Integration);
    msg.getStream(i).removeTopic(RooFit::Plotting);
  }

  // CMS plot style
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");

  // ------------------------------------------------------------------
  // load and reduce dataset
  // ------------------------------------------------------------------
  auto fData = TFile::Open(DATA_ROOT);
  if (!fData || fData->IsZombie()) {
    std::cerr << "Can't open the input data\n";
    return;
  }

  auto dsTmp1 = (RooDataSet *)fData->Get("dataset");
  if (!dsTmp1)
  {
    std::cerr << "Can't read the dataset\n";
    return;
  }
  TString cutBasic = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g)", ptLow, ptHigh, yLow, yHigh);

  auto dsData = std::unique_ptr<RooDataSet>((RooDataSet *)dsTmp1->reduce(cutBasic));
  auto ctau3DTmp = (RooRealVar *)dsData->get()->find("ctau3D");
  auto ctau3DTmpErr = (RooRealVar *)dsData->get()->find("ctau3DErr");

  // ------------------------------------------------------------------
  // find ctau3DTmp range
  // ------------------------------------------------------------------
  double ctLow, ctHigh;
  double errLow, errHigh;
  {
    printf("ctau3DTmp: [%6.g, %.6g]\n", ctau3DTmp->getMin(), ctau3DTmp->getMax());
    err_range::RangeStats stats;
    if (!err_range::compute_ctau3d_range(*dsData, *ctau3DTmp, err_range::kCtau3D_Mass, ptHigh, ctLow, ctHigh, stats))
    {
      printf("ERROR: ctau3D range empty\n");
      return;
    }
    printf("collected %zu values (v>0)\n", stats.n);
    printf("raw min/max: [%.6g, %.6f]\n", stats.rawMin, stats.rawMax);
    printf("new ctau3DTmp range: [%.6g, %.6g]\n", ctLow, ctHigh);
  }

  // ------------------------------------------------------------------
  // find ctau3DTmpErr range
  // ------------------------------------------------------------------
  {
    printf("ctau3DTmpErr: [%6.g, %.6g]\n", ctau3DTmpErr->getMin(), ctau3DTmpErr->getMax());
    err_range::RangeStats stats;
    if (!err_range::compute_err_range(*dsData, *ctau3DTmpErr, err_range::kErr_Default, errLow, errHigh, stats))
    {
      printf("ERROR: ctau3DErr range empty\n");
      return;
    }
    printf("collected %zu values (v>0)\n", stats.n);
    printf("raw min/max: [%.6g, %.6g]\n", stats.rawMin, stats.rawMax);
    printf("new ctau3DTmpErr range: [%.6g, %.6g]\n", errLow, errHigh);
  }

  TString cutAll = Form("(mass>2.6 && mass<3.5) && (pt > %g && pt < %g) && (fabs(y) > %g && fabs(y) < %g) && (ctau3D >= %g && ctau3D <= %g) && (ctau3DErr >= %g && ctau3DErr <= %g)", ptLow, ptHigh, yLow, yHigh, ctLow, ctHigh, errLow, errHigh);

  auto dsFinal = std::unique_ptr<RooDataSet>((RooDataSet *)dsData->reduce(cutAll));
  dsFinal->Print();

  RooRealVar *mass = new RooRealVar("mass", "", 2.6, 3.5);
  // auto mass = (RooRealVar *)dsFinal->get()->find("mass");
  auto pt = (RooRealVar *)dsFinal->get()->find("pt");
  auto rapidity = (RooRealVar *)dsFinal->get()->find("y");
  auto ctau3D = (RooRealVar *)dsFinal->get()->find("ctau3D");
  RooRealVar *ctau3DErr = new RooRealVar("ctau3DErr", "", errLow, errHigh);
  // auto ctau3DErr = (RooRealVar *)dsFinal->get()->find("ctau3DErr");
  
  ctau3D->setRange("fit", ctLow, ctHigh);
  ctau3DErr->setRange("errFit", errLow, errHigh);
  ctau3D->setRange(ctLow, ctHigh);
  ctau3DErr->setRange(errLow, errHigh);


  // ------------------------------------------------------------------
  // Mass PDFs (signal = 2x CrystalBall + 2x Gaussian; bkg = exp or Chebychev<=6)
  // ------------------------------------------------------------------
  TString massFitFile = Form("%s/%s/mass/mass_model_%s.root", outRootDir.Data(), yTag.Data(), figTag.Data());
  int bkgModelType = 1; // 0=exp, 1=cheby
  int chebyOrder = 3;
  TFile *massCfgFile = nullptr;
  auto readMassInt = [&](const char *name, int fallback)
  {
    if (!massCfgFile)
      return fallback;
    auto *param = dynamic_cast<TParameter<int> *>(massCfgFile->Get(name));
    return param ? param->GetVal() : fallback;
  };
  auto readMassDouble = [&](const char *name, double fallback)
  {
    if (!massCfgFile)
      return fallback;
    auto *param = dynamic_cast<TParameter<double> *>(massCfgFile->Get(name));
    return param ? param->GetVal() : fallback;
  };
  if (!gSystem->AccessPathName(massFitFile))
  {
    massCfgFile = TFile::Open(massFitFile, "READ");
    if (massCfgFile && !massCfgFile->IsZombie())
    {
      const int nBkgExpComponents = readMassInt("nBkgExpComponents", 0);
      const int nBkgChebyOrder = readMassInt("nBkgChebyOrder", chebyOrder);
      if (nBkgExpComponents == 1)
        bkgModelType = 0;
      else if (nBkgChebyOrder > 0)
      {
        bkgModelType = 1;
        chebyOrder = std::max(1, std::min(6, nBkgChebyOrder));
      }
      TString bkgLabel = (bkgModelType == 0) ? "exp" : Form("cheby%d", chebyOrder);
      printf("[INFO] loaded mass bkg model from mass.C: %s\n", bkgLabel.Data());
    }
    else
    {
      printf("[WARN] failed to open mass model file: %s\n", massFitFile.Data());
      massCfgFile = nullptr;
    }
  }
  RooRealVar m0("m0", "J/psi mean", 3.0949, 3.05, 3.15);

  RooRealVar msig1("msig1", "CB1 sigma", 0.012924, 0.010, 0.180);
  RooRealVar alpha1("alpha1", "CB1 alpha", 0.48388, 0.1, 3.5);
  RooRealVar n1("n1", "CB1 n", 0.10909, 0.1, 200.0);
  RooCBShape CB1("CB1", "CB1", *mass, m0, msig1, alpha1, n1);

  RooRealVar msig2("msig2", "CB2 sigma", 0.040733, 0.010, 0.2);
  RooRealVar alpha2("alpha2", "CB2 alpha", -1.7065, -10, 5.0);
  RooRealVar n2("n2", "CB2 n", 1.2461, 0.1, 200.0);
  RooCBShape CB2("CB2", "CB2", *mass, m0, msig2, alpha2, n2);

  RooRealVar gsig("gsig", "Gauss sigma", 0.023331, 0.005, 0.4);
  RooGaussian Gm("Gm", "Gm", *mass, m0, gsig);
  RooRealVar gsig2("gsig2", "Gauss2 sigma", 0.14936, 0.010, 0.150);
  RooGaussian Gm2("Gm2", "Gm2", *mass, m0, gsig2);

  RooRealVar fCB1("fCB1", "fCB1", 0.15258, 0.0, 1.0);
  RooRealVar fCB2("fCB2", "fCB2", 0.36694, 0.0, 1.0);
  RooRealVar fG1("fG1", "fG1", 0.76131, 0.0, 1.0);
  RooAddPdf MSig("MSig", "MSig", RooArgList(CB1, CB2, Gm, Gm2), RooArgList(fCB1, fCB2, fG1), true);

  RooRealVar lambda("lambda", "exp slope", -1.0, -10.0, 0.0);
  RooRealVar c1("c1", "c1", -0.60224, -1.0, 1.0);
  RooRealVar c2("c2", "c2", -0.013589, -1.0, 1.0);
  RooRealVar c3("c3", "c3", 0.0081393, -1.0, 1.0);
  RooRealVar c4("c4", "c4", 0.0, -1.0, 1.0);
  RooRealVar c5("c5", "c5", 0.0, -1.0, 1.0);
  RooRealVar c6("c6", "c6", 0.0, -1.0, 1.0);
  std::unique_ptr<RooChebychev> MBkgCheby;
  std::unique_ptr<RooExponential> MBkgExp;
  RooAbsPdf *MBkg = nullptr;
  if (bkgModelType == 0)
  {
    MBkgExp = std::make_unique<RooExponential>("MBkg", "MBkg", *mass, lambda);
    MBkg = MBkgExp.get();
  }
  else
  {
    chebyOrder = std::max(1, std::min(6, chebyOrder));
    RooArgList chebyCoeffs;
    RooRealVar *coeffs[] = {&c1, &c2, &c3, &c4, &c5, &c6};
    for (int i = 0; i < chebyOrder; ++i)
      chebyCoeffs.add(*coeffs[i]);
    MBkgCheby = std::make_unique<RooChebychev>("MBkg", "MBkg", *mass, chebyCoeffs);
    MBkg = MBkgCheby.get();
  }

  const double nData = static_cast<double>(dsFinal->numEntries());
  RooRealVar nSigMass("nSigMass", "nSigMass", 3.4370e4, 0.0, nData * 1.5);
  RooRealVar nBkgMass("nBkgMass", "nBkgMass", 2.0480, 0.0, nData * 1.5);
  RooAddPdf massModel("massModel", "massModel", RooArgList(MSig, *MBkg), RooArgList(nSigMass, nBkgMass));
  printf("[OK] built mass model (bkg=%s%s)\n",
         (bkgModelType == 0 ? "exp" : "cheby"),
         (bkgModelType == 0 ? "" : Form("%d", chebyOrder)));

  if (massCfgFile)
  {
    if (bkgModelType == 0)
    {
      lambda.setVal(readMassDouble("bkg_mass_lambda", lambda.getVal()));
      lambda.setConstant(true);
    }
    else
    {
      RooRealVar *coeffs[] = {&c1, &c2, &c3, &c4, &c5, &c6};
      const char *names[] = {"bkg_mass_p1", "bkg_mass_p2", "bkg_mass_p3", "bkg_mass_p4", "bkg_mass_p5", "bkg_mass_p6"};
      for (int i = 0; i < chebyOrder; ++i)
      {
        coeffs[i]->setVal(readMassDouble(names[i], coeffs[i]->getVal()));
        coeffs[i]->setConstant(true);
      }
    }
  }

  // ------------------------------------------------------------------
  // constraint overrides (mass constraints only)
  // ------------------------------------------------------------------
  // OverrideMode overrides ConstraintMode per RooRealVar name
  const std::unordered_map<std::string, ConstraintMode> constraintOverrideModes = {
    // // mass
    // {"msig2", kFix},
    // {"gsig", kFix},
    // {"m0", kFix},
    // {"fCB2", kFix},
    // {"nSigMass", kFix},
    // {"nBkgMass", kFix},

  };
  const std::unordered_map<std::string, double> constraintOverrideValues = {
      // {"biasPr", 0.0},
      // {"biasBkg", 0.0},
  };

  RooArgSet mcConstraints;
  if (gSystem->AccessPathName(massFitFile))
  {
    printf("[WARN] mass fit file not found: %s\n", massFitFile.Data());
  }
  else
  {
    auto fixFromMass = [&](RooRealVar &dst, const char *massName)
    {
      const double val = readMassDouble(massName, std::numeric_limits<double>::quiet_NaN());
      if (!std::isfinite(val))
      {
        printf("[WARN] mass param '%s' not found in %s\n", massName, massFitFile.Data());
        return;
      }
      dst.setVal(val);
      dst.setConstant(true);
    };

    if (bkgModelType == 0)
    {
      fixFromMass(lambda, "bkg_mass_lambda");
    }
    else
    {
      fixFromMass(c1, "bkg_mass_p1");
      if (chebyOrder >= 2)
        fixFromMass(c2, "bkg_mass_p2");
      if (chebyOrder >= 3)
        fixFromMass(c3, "bkg_mass_p3");
      if (chebyOrder >= 4)
        fixFromMass(c4, "bkg_mass_p4");
      if (chebyOrder >= 5)
        fixFromMass(c5, "bkg_mass_p5");
      if (chebyOrder >= 6)
        fixFromMass(c6, "bkg_mass_p6");
    }
  }
  printf("[OK] loaded background model choice and fixed shape from mass.C output\n");
  
  // ------------------------------------------------------------------
  // make punzi term - bkg (mass sidebands)
  // ------------------------------------------------------------------
  std::unique_ptr<RooDataSet> dsSideband;
  std::unique_ptr<RooDataSet> dsSigRegion;
  std::unique_ptr<RooDataSet> dsSidebandErr;
  std::unique_ptr<RooDataSet> dsSigRegionErr;
  std::unique_ptr<RooAddPdf> pdfErrComboSig;
  RooFitResult *fitErrBkgPtr = nullptr;
  RooFitResult *fitErrSigPtr = nullptr;
  RooAbsPdf *pdfErrSig = nullptr;

  const double massSb1Low = 2.6, massSb1High = 2.8;
  const double massSb2Low = 3.3, massSb2High = 3.5;
  const double massSigLow = 2.8, massSigHigh = 3.3;
  const int errBins = std::max(2, ctau3DErr->getBins() * 2);
  double errLowBkg = errLow;
  double errHighBkg = errHigh;
  double errLowSig = errLow;
  double errHighSig = errHigh;
  int bkgBins = errBins;
  int sigBins = errBins;
  int bkgFinalBins = -1;
  auto clamp_comp = [](int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); };
  const int nErrBkgGauss = clamp_comp(nErrBkgGaussUser, 0, 3);
  const int nErrBkgLandau = clamp_comp(nErrBkgLandauUser, 0, 3);
  const int nErrBkgExp = clamp_comp(nErrBkgExpUser, 0, 1);
  const int nErrSigGauss = clamp_comp(nErrSigGaussUser, 0, 3);
  const int nErrSigLandau = clamp_comp(nErrSigLandauUser, 0, 3);
  const int nErrSigExp = clamp_comp(nErrSigExpUser, 0, 1);
  if (nErrBkgGauss + nErrBkgLandau + nErrBkgExp == 0 || nErrSigGauss + nErrSigLandau + nErrSigExp == 0)
  {
    printf("[ERROR] invalid err-PDF setup: at least one component is required for both bkg and sig.\n");
    return;
  }
  printf("[INFO] errPdf components: bkg(G=%d,L=%d,E=%d), sig(G=%d,L=%d,E=%d)\n",
         nErrBkgGauss, nErrBkgLandau, nErrBkgExp, nErrSigGauss, nErrSigLandau, nErrSigExp);

  ctau3DErr->setBins(errBins);
  RooArgSet errOnly(*ctau3DErr);
  TString cutSideband = Form("((mass>%g && mass<%g) || (mass>%g && mass<%g))", massSb1Low, massSb1High, massSb2Low, massSb2High);
  TString cutSignal = Form("(mass>%g && mass<%g)", massSigLow, massSigHigh);

  if (doErrPdfFit)
  {
    dsSideband.reset((RooDataSet *)dsFinal->reduce(cutSideband));
    dsSigRegion.reset((RooDataSet *)dsFinal->reduce(cutSignal));

    TString cutErr = Form("(ctau3DErr >= %g && ctau3DErr <= %g)", errLow, errHigh);
    dsSidebandErr.reset((RooDataSet *)dsSideband->reduce(errOnly, cutErr.Data()));
    dsSigRegionErr.reset((RooDataSet *)dsSigRegion->reduce(errOnly, cutErr.Data()));

    if (!dsSidebandErr || !dsSigRegionErr)
    {
      printf("ERROR: err-range reduced dataset missing\n");
      return;
    }
  }

  if (doErrPdfFit)
  {
    errLowBkg = errLow;
    errHighBkg = errHigh;
    bkgBins = errBins;
    printf("[INFO] err range (bkg): [%.6g, %.6g], bins=%d\n", errLowBkg, errHighBkg, bkgBins);

    errLowSig = errLowBkg;
    errHighSig = errHighBkg;
    sigBins = bkgBins;
    printf("[INFO] err range (sig): [%.6g, %.6g], bins=%d\n", errLowSig, errHighSig, sigBins);
  }

  std::unique_ptr<TH1> hErrBkgSB;
  if (doErrPdfFit)
  {
    hErrBkgSB.reset(dsSidebandErr->createHistogram("hErrBkgSB", *ctau3DErr, Binning(bkgBins, errLowBkg, errHighBkg)));
    if (!hErrBkgSB)
    {
      printf("ERROR: failed to create err histogram (bkg)\n");
      return;
    }
  }

  auto countNegBins = [&](const TH1 *h) -> int
  {
    int n = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      if (h->GetBinContent(i) < 0.0)
        n++;
    }
    return n;
  };

  auto clampTinyBins = [&](TH1 *h, double frac, const char *tag) -> void
  {
    if (!h)
      return;
    double maxv = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      maxv = std::max(maxv, h->GetBinContent(i));
    if (!(maxv > 0.0))
      return;
    const double thr = maxv * frac;
    int nTiny = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      const double v = h->GetBinContent(i);
      if (v > 0.0 && v < thr)
      {
        h->SetBinContent(i, 0.1);
        nTiny++;
      }
    }
    if (nTiny > 0)
      printf("[INFO] %s: clamped %d tiny bins (< %.3g) to 0.1\n", tag, nTiny, thr);
  };
  std::unique_ptr<RooAddPdf> pdfErrCombo;
  RooAbsPdf *pdfErr = nullptr;
  const double bkgSpan = std::max(1e-4, errHighBkg - errLowBkg);
  auto q_hist = [](TH1 *h, double q, double fallback) -> double
  {
    if (!h || h->Integral() <= 0.0)
      return fallback;
    double xq[1] = {q};
    double yq[1] = {fallback};
    h->GetQuantiles(1, yq, xq);
    return std::isfinite(yq[0]) ? yq[0] : fallback;
  };
  auto clamp_range = [](double v, double lo, double hi) { return std::max(lo, std::min(v, hi)); };
  auto widen_range = [&](std::pair<double, double> r, double lo, double hi) -> std::pair<double, double>
  {
    const double c = 0.5 * (r.first + r.second);
    const double w = std::max(1e-6, 0.5 * (r.second - r.first) * errShapeRangeScale);
    double a = std::max(lo, c - w);
    double b = std::min(hi, c + w);
    if (!(b > a + 1e-6))
      return {lo, hi};
    return {a, b};
  };
  auto make_range = [&](double lo, double hi) -> std::pair<double, double>
  {
    lo = std::max(lo, errLowBkg);
    hi = std::min(hi, errHighBkg);
    if (!(hi > lo + 1e-6))
      return {errLowBkg, errHighBkg};
    return {lo, hi};
  };
  const double bkgQ05 = q_hist(hErrBkgSB.get(), 0.05, errLowBkg + 0.05 * bkgSpan);
  const double bkgQ20 = q_hist(hErrBkgSB.get(), 0.20, errLowBkg + 0.20 * bkgSpan);
  const double bkgQ35 = q_hist(hErrBkgSB.get(), 0.35, errLowBkg + 0.35 * bkgSpan);
  const double bkgQ50 = q_hist(hErrBkgSB.get(), 0.50, errLowBkg + 0.50 * bkgSpan);
  const double bkgQ65 = q_hist(hErrBkgSB.get(), 0.65, errLowBkg + 0.65 * bkgSpan);
  const double bkgQ85 = q_hist(hErrBkgSB.get(), 0.85, errLowBkg + 0.85 * bkgSpan);
  const double bkgQ80 = q_hist(hErrBkgSB.get(), 0.80, errLowBkg + 0.80 * bkgSpan);
  const double bkgQ95 = q_hist(hErrBkgSB.get(), 0.95, errLowBkg + 0.95 * bkgSpan);
  const auto bkgR1 = make_range(bkgQ05, bkgQ50);
  const auto bkgR2 = make_range(bkgQ20, bkgQ80);
  const auto bkgR3 = make_range(bkgQ50, bkgQ95);
  const auto bkgRL1 = make_range(bkgQ05, bkgQ65);
  const auto bkgRL2 = make_range(bkgQ35, bkgQ85);
  const auto bkgRL3 = make_range(bkgQ65, errHighBkg);
  const auto bkgR1A = widen_range(bkgR1, errLowBkg, errHighBkg);
  const auto bkgR2A = widen_range(bkgR2, errLowBkg, errHighBkg);
  const auto bkgR3A = widen_range(bkgR3, errLowBkg, errHighBkg);
  const auto bkgRL1A = widen_range(bkgRL1, errLowBkg, errHighBkg);
  const auto bkgRL2A = widen_range(bkgRL2, errLowBkg, errHighBkg);
  const auto bkgRL3A = widen_range(bkgRL3, errLowBkg, errHighBkg);
  const double bkgBinW = bkgSpan / std::max(1, bkgBins);
  const double bkgS1Min = std::max(1e-4, 0.25 * bkgBinW);
  const double bkgS2Min = std::max(2e-4, 0.35 * bkgBinW);
  const double bkgS3Min = std::max(3e-4, 0.45 * bkgBinW);
  const double bkgS1Max = std::max(1.5 * bkgS1Min, std::min(0.03, 0.90 * (bkgR1A.second - bkgR1A.first)));
  const double bkgS2Max = std::max(1.5 * bkgS2Min, std::min(0.05, 0.90 * (bkgR2A.second - bkgR2A.first)));
  const double bkgS3Max = std::max(1.5 * bkgS3Min, std::min(0.09, 0.90 * (bkgR3A.second - bkgR3A.first)));
  const double bkgSigma1Init = clamp_range(0.50 * bkgS1Max, bkgS1Min, bkgS1Max);
  const double bkgSigma2Init = clamp_range(0.60 * bkgS2Max, bkgS2Min, bkgS2Max);
  const double bkgSigma3Init = clamp_range(0.70 * bkgS3Max, bkgS3Min, bkgS3Max);
  const double bkgMean1Init = clamp_range(bkgQ20, bkgR1A.first, bkgR1A.second);
  const double bkgMean2Init = clamp_range(bkgQ50, bkgR2A.first, bkgR2A.second);
  const double bkgMean3Init = clamp_range(bkgQ80, bkgR3A.first, bkgR3A.second);
  const double bkgMpv1Init = clamp_range(bkgQ35, bkgRL1A.first, bkgRL1A.second);
  const double bkgMpv2Init = clamp_range(bkgQ65, bkgRL2A.first, bkgRL2A.second);
  const double bkgMpv3Init = clamp_range(bkgQ80, bkgRL3A.first, bkgRL3A.second);
  RooRealVar meanGBkg("meanGBkg", "", bkgMean1Init, bkgR1A.first, bkgR1A.second);
  RooRealVar sigmaGBkg("sigmaGBkg", "", bkgSigma1Init, bkgS1Min, bkgS1Max);
  RooGaussian gausBkg("gausBkg", "", *ctau3DErr, meanGBkg, sigmaGBkg);
  RooRealVar meanGBkg2("meanGBkg2", "", bkgMean2Init, bkgR2A.first, bkgR2A.second);
  RooRealVar sigmaGBkg2("sigmaGBkg2", "", bkgSigma2Init, bkgS2Min, bkgS2Max);
  RooGaussian gausBkg2("gausBkg2", "", *ctau3DErr, meanGBkg2, sigmaGBkg2);
  RooRealVar meanGBkg3("meanGBkg3", "", bkgMean3Init, bkgR3A.first, bkgR3A.second);
  RooRealVar sigmaGBkg3("sigmaGBkg3", "", bkgSigma3Init, bkgS3Min, bkgS3Max);
  RooGaussian gausBkg3("gausBkg3", "", *ctau3DErr, meanGBkg3, sigmaGBkg3);
  RooRealVar mpvBkg("mpvBkg", "", bkgMpv1Init, bkgRL1A.first, bkgRL1A.second);
  const double bkgLandauSigmaMin = std::max((isPt23 ? 8e-4 : 2e-4), bkgS2Min);
  RooRealVar sigmaLBkg("sigmaLBkg", "", bkgSigma2Init, bkgLandauSigmaMin, bkgS2Max);
  RooLandau landauBkg("landauBkg", "", *ctau3DErr, mpvBkg, sigmaLBkg);
  RooRealVar mpvBkg2("mpvBkg2", "", bkgMpv2Init, bkgRL2A.first, bkgRL2A.second);
  RooRealVar sigmaLBkg2("sigmaLBkg2", "", bkgSigma2Init, bkgLandauSigmaMin, bkgS2Max);
  RooLandau landauBkg2("landauBkg2", "", *ctau3DErr, mpvBkg2, sigmaLBkg2);
  RooRealVar mpvBkg3("mpvBkg3", "", bkgMpv3Init, bkgRL3A.first, bkgRL3A.second);
  RooRealVar sigmaLBkg3("sigmaLBkg3", "", bkgSigma3Init, bkgS3Min, bkgS3Max);
  RooLandau landauBkg3("landauBkg3", "", *ctau3DErr, mpvBkg3, sigmaLBkg3);
  RooRealVar lambdaExpBkg("lambdaExpBkg", "", -1.0 / std::max(1e-4, bkgQ65 - bkgQ20), -5000.0, -1e-3);
  RooExponential expBkg("expBkg", "", *ctau3DErr, lambdaExpBkg);
  RooRealVar fBkgG1("fBkgG1", "", 0.55, 0.0, 1.0);
  RooRealVar fBkgG2("fBkgG2", "", 0.25, 0.0, 1.0);
  RooRealVar fBkg3("fBkg3", "", 0.2, 0.0, 1.0);
  RooRealVar fBkg4("fBkg4", "", 0.2, 0.0, 1.0);
  RooRealVar fBkg5("fBkg5", "", 0.2, 0.0, 1.0);
  RooRealVar fBkg6("fBkg6", "", 0.2, 0.0, 1.0);
  RooArgList bkgCompList;
  if (nErrBkgGauss >= 1) bkgCompList.add(gausBkg);
  if (nErrBkgGauss >= 2) bkgCompList.add(gausBkg2);
  if (nErrBkgGauss >= 3) bkgCompList.add(gausBkg3);
  if (nErrBkgLandau >= 1) bkgCompList.add(landauBkg);
  if (nErrBkgLandau >= 2) bkgCompList.add(landauBkg2);
  if (nErrBkgLandau >= 3) bkgCompList.add(landauBkg3);
  if (nErrBkgExp >= 1) bkgCompList.add(expBkg);
  RooArgList bkgFracList;
  if (bkgCompList.getSize() >= 2) bkgFracList.add(fBkgG1);
  if (bkgCompList.getSize() >= 3) bkgFracList.add(fBkgG2);
  if (bkgCompList.getSize() >= 4) bkgFracList.add(fBkg3);
  if (bkgCompList.getSize() >= 5) bkgFracList.add(fBkg4);
  if (bkgCompList.getSize() >= 6) bkgFracList.add(fBkg5);
  if (bkgCompList.getSize() >= 7) bkgFracList.add(fBkg6);
  pdfErrCombo = std::make_unique<RooAddPdf>("pdfErr", "", bkgCompList, bkgFracList);

  std::unique_ptr<RooDataHist> hErrBkg;
  std::unique_ptr<RooFitResult> fitErrBkg;
  if (doErrPdfFit && hErrBkgSB)
  {
    int nNegBkg = countNegBins(hErrBkgSB.get());
    bkgFinalBins = hErrBkgSB->GetNbinsX();

    if (nNegBkg > 0)
    {
      for (int i = 1; i <= hErrBkgSB->GetNbinsX(); ++i)
      {
        if (hErrBkgSB->GetBinContent(i) < 0.0)
          hErrBkgSB->SetBinContent(i, 0.1);
      }
      printf("[INFO] bkg err negative bins -> clamped to 0.1\n");
    }

    clampTinyBins(hErrBkgSB.get(), 1e-4, "hErrBkgSB");
    hErrBkg = std::make_unique<RooDataHist>("hErrBkg", "", *ctau3DErr, hErrBkgSB.get());
    if (!hErrBkg || hErrBkg->numEntries() == 0)
    {
      printf("ERROR: hErrBkg empty\n");
      return;
    }
  }
  if (doErrPdfFit && bkgFinalBins > 0)
  {
    errLowSig = errLowBkg;
    errHighSig = errHighBkg;
    sigBins = bkgFinalBins;
    printf("[INFO] align sig err hist/pdf to bkg: range[%.6g, %.6g], bins=%d\n", errLowSig, errHighSig, sigBins);
  }

  auto prefitFraction = [&](double n) -> float
  {
    float frac = 0.1f;
    if (n >= 200000)
      frac = 0.05f;
    if (n >= 500000)
      frac = 0.025f;
    return frac;
  };
  auto fit_failed = [&](RooFitResult *fr) -> bool
  {
    if (!fr)
      return true;
    const int status = fr->status();
    int hesse = 0;
    bool hasHesse = false;
    for (UInt_t i = 0, n = fr->numStatusHistory(); i < n; ++i)
    {
      auto lab = fr->statusLabelHistory(i);
      if (lab && !strcmp(lab, "HESSE"))
      {
        hesse = fr->statusCodeHistory(i);
        hasHesse = true;
        break;
      }
    }
    return (status != 0) || (hasHesse && hesse != 0);
  };

  auto refitBkg = [&]() -> void
  {
    if (!doErrPdfFit || !hErrBkg)
      return;
    const float frac = prefitFraction(hErrBkg->sumEntries());
    auto run_bkg_fit = [&](RooArgSet *ext) -> std::unique_ptr<RooFitResult>
    {
      if (ext && ext->getSize() > 0)
      {
        return std::unique_ptr<RooFitResult>(
            pdfErrCombo->fitTo(*hErrBkg, Save(), PrintLevel(-1), PrintEvalErrors(-1), SumW2Error(false),
                               PrefitDataFraction(frac), ExternalConstraints(*ext)));
      }
      return std::unique_ptr<RooFitResult>(
          pdfErrCombo->fitTo(*hErrBkg, Save(), PrintLevel(-1), PrintEvalErrors(-1), SumW2Error(false),
                             PrefitDataFraction(frac)));
    };
    fitErrBkg = run_bkg_fit(nullptr);
    if (fitErrBkg)
    {
      fitErrBkg->Print("V");
      fitErrBkgPtr = fitErrBkg.get();
    }
  };

  if (doErrPdfFit && hErrBkgSB)
  {
    meanGBkg.setVal(bkgMean1Init);
    sigmaGBkg.setVal(bkgSigma1Init);
    meanGBkg2.setVal(bkgMean2Init);
    sigmaGBkg2.setVal(bkgSigma2Init);
    meanGBkg3.setVal(bkgMean3Init);
    sigmaGBkg3.setVal(bkgSigma3Init);
    mpvBkg.setVal(bkgMpv1Init);
    sigmaLBkg.setVal(bkgSigma2Init);
    mpvBkg2.setVal(bkgMpv2Init);
    sigmaLBkg2.setVal(bkgSigma2Init);
    mpvBkg3.setVal(bkgMpv3Init);
    sigmaLBkg3.setVal(bkgSigma3Init);
    lambdaExpBkg.setVal(-1.0 / std::max(1e-4, bkgQ65 - bkgQ20));
    if (isPt23)
    {
      // Use stable seeds/ranges tuned for the pt 2-3 sideband err shape.
      fBkgG1.setVal(0.35);
      fBkgG2.setVal(0.60);
      sigmaLBkg.setVal(std::max(0.0035, bkgLandauSigmaMin));
      sigmaLBkg2.setVal(std::max(0.0020, bkgLandauSigmaMin));
    }
  }
  if (doErrPdfFit)
    refitBkg();
  pdfErr = pdfErrCombo.get();

  // ------------------------------------------------------------------
  // make punzi term - sig (signal region - scaled sideband)
  // 
  // ------------------------------------------------------------------
  if (doErrPdfFit)
  {
    if (!dsSideband || !dsSigRegion)
    {
      printf("ERROR: sideband or signal dataset missing\n");
      return;
    }
  }

  // errhist - SR, SB
  std::unique_ptr<TH1> hErrSigSR;
  std::unique_ptr<TH1> hErrBkgSB2;
  if (doErrPdfFit)
  {
    hErrSigSR.reset(dsSigRegion->createHistogram("hErrSigSR", *ctau3DErr, Binning(sigBins, errLowSig, errHighSig)));
    hErrBkgSB2.reset(dsSideband->createHistogram("hErrBkgSB2", *ctau3DErr, Binning(sigBins, errLowSig, errHighSig)));
    if (!hErrSigSR || !hErrBkgSB2)
    {
      printf("ERROR: failed to create err histograms\n");
      return;
    }
  }

  // mass window range
  mass->setRange("SR", massSigLow, massSigHigh);
  mass->setRange("SB1", massSb1Low, massSb1High);
  mass->setRange("SB2", massSb2Low, massSb2High);

  // MBkg ration: SR/SB
  auto iBkgSR = std::unique_ptr<RooAbsReal>(MBkg->createIntegral(*mass, Range("SR")));
  auto iBkgSB1 = std::unique_ptr<RooAbsReal>(MBkg->createIntegral(*mass, Range("SB1")));
  auto iBkgSB2 = std::unique_ptr<RooAbsReal>(MBkg->createIntegral(*mass, Range("SB2")));

  double fracBkgSR = (iBkgSR ? iBkgSR->getVal() : 0.0);
  double fracBkgSB = (iBkgSB1 ? iBkgSB1->getVal() : 0.0) + (iBkgSB2 ? iBkgSB2->getVal() : 0.0);

  // NBkg_XX = nBkgMass * fracBkgXX
  double NBkg_SR = nBkgMass.getVal() * fracBkgSR;
  double NBkg_SB = nBkgMass.getVal() * fracBkgSB;

  double scaleBkg = 0.0;
  if (NBkg_SB > 0.0)
    scaleBkg = NBkg_SR / NBkg_SB;
  else
    printf("[WARN] NBkg_SB <= 0. Cannot scale sideband -> signal region.\n");

  // subtraction: Err(SR) - scale * Err(SB)
  printf("[MassBkgScale] nBkgMass=%.1f frac(SR)=%.6f frac(SB)=%.6f => NBkg_SR=%.1f NBkg_SB=%.1f scale=%.6f\n",
         nBkgMass.getVal(), fracBkgSR, fracBkgSB, NBkg_SR, NBkg_SB, scaleBkg);

  std::unique_ptr<RooDataHist> hErrSig;
  if (doErrPdfFit && hErrSigSR && hErrBkgSB2)
  {
    hErrSigSR->Add(hErrBkgSB2.get(), -scaleBkg);

    int nNeg = countNegBins(hErrSigSR.get());
    if (nNeg > 0)
      printf("[WARN] sig err has %d negative bins (no rebin requested)\n", nNeg);

    // Keep SR subtraction bins physical for log-y plots: negative -> zero (empty bin).
    if (nNeg > 0)
    {
      for (int i = 1; i <= hErrSigSR->GetNbinsX(); ++i)
      {
        if (hErrSigSR->GetBinContent(i) < 0.0)
        {
          hErrSigSR->SetBinContent(i, 0.0);
          hErrSigSR->SetBinError(i, 0.0);
        }
      }
      printf("[INFO] subtraction produced %d negative bins -> clamped to 0 (err=0)\n", nNeg);
    }

    // Do not inflate tiny/zero SR bins to positive constants.
    hErrSig = std::make_unique<RooDataHist>("hErrSig", "", *ctau3DErr, hErrSigSR.get());
    if (!hErrSig || hErrSig->sumEntries() <= 0.0)
    {
      printf("ERROR: hErrSig empty after subtraction\n");
      return;
    }
  }

  const double sigSpan = std::max(1e-4, errHighSig - errLowSig);
  auto make_sig_range = [&](double lo, double hi) -> std::pair<double, double>
  {
    lo = std::max(lo, errLowSig);
    hi = std::min(hi, errHighSig);
    if (!(hi > lo + 1e-6))
      return {errLowSig, errHighSig};
    return {lo, hi};
  };
  const double sigQ05 = q_hist(hErrSigSR.get(), 0.05, errLowSig + 0.05 * sigSpan);
  const double sigQ20 = q_hist(hErrSigSR.get(), 0.20, errLowSig + 0.20 * sigSpan);
  const double sigQ35 = q_hist(hErrSigSR.get(), 0.35, errLowSig + 0.35 * sigSpan);
  const double sigQ50 = q_hist(hErrSigSR.get(), 0.50, errLowSig + 0.50 * sigSpan);
  const double sigQ65 = q_hist(hErrSigSR.get(), 0.65, errLowSig + 0.65 * sigSpan);
  const double sigQ80 = q_hist(hErrSigSR.get(), 0.80, errLowSig + 0.80 * sigSpan);
  const double sigQ85 = q_hist(hErrSigSR.get(), 0.85, errLowSig + 0.85 * sigSpan);
  const double sigQ95 = q_hist(hErrSigSR.get(), 0.95, errLowSig + 0.95 * sigSpan);
  const auto sigR1 = make_sig_range(sigQ05, sigQ50);
  const auto sigR2 = make_sig_range(sigQ20, sigQ80);
  const auto sigR3 = make_sig_range(sigQ50, sigQ95);
  const auto sigRL1 = make_sig_range(sigQ05, sigQ65);
  const auto sigRL2 = make_sig_range(sigQ35, sigQ85);
  const auto sigRL3 = make_sig_range(sigQ65, errHighSig);
  const auto sigR1A = widen_range(sigR1, errLowSig, errHighSig);
  const auto sigR2A = widen_range(sigR2, errLowSig, errHighSig);
  const auto sigR3A = widen_range(sigR3, errLowSig, errHighSig);
  const auto sigRL1A = widen_range(sigRL1, errLowSig, errHighSig);
  const auto sigRL2A = widen_range(sigRL2, errLowSig, errHighSig);
  const auto sigRL3A = widen_range(sigRL3, errLowSig, errHighSig);
  const double sigBinW = sigSpan / std::max(1, sigBins);
  const double sigS1Min = std::max(1e-4, 0.25 * sigBinW);
  const double sigS2Min = std::max(2e-4, 0.35 * sigBinW);
  const double sigS3Min = std::max(3e-4, 0.45 * sigBinW);
  const double sigS1Max = std::max(1.5 * sigS1Min, std::min(0.03, 0.90 * (sigR1A.second - sigR1A.first)));
  const double sigS2Max = std::max(1.5 * sigS2Min, std::min(0.05, 0.90 * (sigR2A.second - sigR2A.first)));
  const double sigS3Max = std::max(1.5 * sigS3Min, std::min(0.09, 0.90 * (sigR3A.second - sigR3A.first)));
  const double sigSigma1Init = clamp_range(0.5 * sigS1Max, sigS1Min, sigS1Max);
  const double sigSigma2Init = clamp_range(0.6 * sigS2Max, sigS2Min, sigS2Max);
  const double sigSigma3Init = clamp_range(0.7 * sigS3Max, sigS3Min, sigS3Max);
  const double sigMean1Init = clamp_range(sigQ20, sigR1A.first, sigR1A.second);
  const double sigMean2Init = clamp_range(sigQ50, sigR2A.first, sigR2A.second);
  const double sigMean3Init = clamp_range(sigQ80, sigR3A.first, sigR3A.second);
  const double sigMpv1Init = clamp_range(sigQ35, sigRL1A.first, sigRL1A.second);
  const double sigMpv2Init = clamp_range(sigQ65, sigRL2A.first, sigRL2A.second);
  const double sigMpv3Init = clamp_range(sigQ80, sigRL3A.first, sigRL3A.second);
  RooRealVar meanGSig("meanGSig", "", sigMean1Init, sigR1A.first, sigR1A.second);
  RooRealVar sigmaGSig("sigmaGSig", "", sigSigma1Init, sigS1Min, sigS1Max);
  RooGaussian gausSig("gausSig", "", *ctau3DErr, meanGSig, sigmaGSig);
  RooRealVar meanGSig2("meanGSig2", "", sigMean2Init, sigR2A.first, sigR2A.second);
  RooRealVar sigmaGSig2("sigmaGSig2", "", sigSigma2Init, sigS2Min, sigS2Max);
  RooGaussian gausSig2("gausSig2", "", *ctau3DErr, meanGSig2, sigmaGSig2);
  RooRealVar meanGSig3("meanGSig3", "", sigMean3Init, sigR3A.first, sigR3A.second);
  RooRealVar sigmaGSig3("sigmaGSig3", "", sigSigma3Init, sigS3Min, sigS3Max);
  RooGaussian gausSig3("gausSig3", "", *ctau3DErr, meanGSig3, sigmaGSig3);
  RooRealVar mpvSig("mpvSig", "", sigMpv1Init, sigRL1A.first, sigRL1A.second);
  RooRealVar sigmaLSig("sigmaLSig", "", sigSigma2Init, sigS2Min, sigS2Max);
  RooLandau landauSig("landauSig", "", *ctau3DErr, mpvSig, sigmaLSig);
  RooRealVar mpvSig2("mpvSig2", "", sigMpv2Init, sigRL2A.first, sigRL2A.second);
  RooRealVar sigmaLSig2("sigmaLSig2", "", sigSigma2Init, sigS2Min, sigS2Max);
  RooLandau landauSig2("landauSig2", "", *ctau3DErr, mpvSig2, sigmaLSig2);
  RooRealVar mpvSig3("mpvSig3", "", sigMpv3Init, sigRL3A.first, sigRL3A.second);
  RooRealVar sigmaLSig3("sigmaLSig3", "", sigSigma3Init, sigS3Min, sigS3Max);
  RooLandau landauSig3("landauSig3", "", *ctau3DErr, mpvSig3, sigmaLSig3);
  RooRealVar lambdaExpSig("lambdaExpSig", "", -1.0 / std::max(1e-4, sigQ65 - sigQ20), -5000.0, -1e-3);
  RooExponential expSig("expSig", "", *ctau3DErr, lambdaExpSig);
  RooRealVar fSigG1("fSigG1", "", 0.6, 0.0, 1.0);
  RooRealVar fSigG2("fSigG2", "", 0.2, 0.0, 1.0);
  RooRealVar fSig3("fSig3", "", 0.2, 0.0, 1.0);
  RooRealVar fSig4("fSig4", "", 0.2, 0.0, 1.0);
  RooRealVar fSig5("fSig5", "", 0.2, 0.0, 1.0);
  RooRealVar fSig6("fSig6", "", 0.2, 0.0, 1.0);
  RooArgList sigCompList;
  if (nErrSigGauss >= 1) sigCompList.add(gausSig);
  if (nErrSigGauss >= 2) sigCompList.add(gausSig2);
  if (nErrSigGauss >= 3) sigCompList.add(gausSig3);
  if (nErrSigLandau >= 1) sigCompList.add(landauSig);
  if (nErrSigLandau >= 2) sigCompList.add(landauSig2);
  if (nErrSigLandau >= 3) sigCompList.add(landauSig3);
  if (nErrSigExp >= 1) sigCompList.add(expSig);
  RooArgList sigFracList;
  if (sigCompList.getSize() >= 2) sigFracList.add(fSigG1);
  if (sigCompList.getSize() >= 3) sigFracList.add(fSigG2);
  if (sigCompList.getSize() >= 4) sigFracList.add(fSig3);
  if (sigCompList.getSize() >= 5) sigFracList.add(fSig4);
  if (sigCompList.getSize() >= 6) sigFracList.add(fSig5);
  if (sigCompList.getSize() >= 7) sigFracList.add(fSig6);
  pdfErrComboSig = std::make_unique<RooAddPdf>("pdfErrSig", "", sigCompList, sigFracList);
  if (doErrPdfFit && hErrSig)
  {
    meanGSig.setVal(sigMean1Init);
    sigmaGSig.setVal(sigSigma1Init);
    meanGSig2.setVal(sigMean2Init);
    sigmaGSig2.setVal(sigSigma2Init);
    meanGSig3.setVal(sigMean3Init);
    sigmaGSig3.setVal(sigSigma3Init);
    mpvSig.setVal(sigMpv1Init);
    sigmaLSig.setVal(sigSigma2Init);
    mpvSig2.setVal(sigMpv2Init);
    sigmaLSig2.setVal(sigSigma2Init);
    mpvSig3.setVal(sigMpv3Init);
    sigmaLSig3.setVal(sigSigma3Init);
    lambdaExpSig.setVal(-1.0 / std::max(1e-4, sigQ65 - sigQ20));
  }
  std::unique_ptr<RooFitResult> fitErrSig;
  if (doErrPdfFit)
  {
    const float frac = hErrSig ? prefitFraction(hErrSig->sumEntries()) : 0.1f;
    auto run_sig_fit = [&](RooArgSet *ext) -> std::unique_ptr<RooFitResult>
    {
      if (ext && ext->getSize() > 0)
      {
        return std::unique_ptr<RooFitResult>(
            pdfErrComboSig->fitTo(*hErrSig, Save(), PrintLevel(-1), PrintEvalErrors(-1), SumW2Error(false),
                                  PrefitDataFraction(frac), ExternalConstraints(*ext)));
      }
      return std::unique_ptr<RooFitResult>(
          pdfErrComboSig->fitTo(*hErrSig, Save(), PrintLevel(-1), PrintEvalErrors(-1), SumW2Error(false)));
    };
    fitErrSig = run_sig_fit(nullptr);
    if (fitErrSig)
    {
      fitErrSig->Print("V");
      fitErrSigPtr = fitErrSig.get();
    }
  }
  pdfErrSig = pdfErrComboSig.get();

  RooArgSet parsBkg;
  parsBkg.add(meanGBkg); parsBkg.add(sigmaGBkg); parsBkg.add(meanGBkg2); parsBkg.add(sigmaGBkg2); parsBkg.add(meanGBkg3); parsBkg.add(sigmaGBkg3);
  parsBkg.add(mpvBkg); parsBkg.add(sigmaLBkg); parsBkg.add(mpvBkg2); parsBkg.add(sigmaLBkg2); parsBkg.add(mpvBkg3); parsBkg.add(sigmaLBkg3);
  parsBkg.add(lambdaExpBkg); parsBkg.add(fBkgG1); parsBkg.add(fBkgG2); parsBkg.add(fBkg3); parsBkg.add(fBkg4); parsBkg.add(fBkg5); parsBkg.add(fBkg6);
  RooArgSet parsSig;
  parsSig.add(meanGSig); parsSig.add(sigmaGSig); parsSig.add(meanGSig2); parsSig.add(sigmaGSig2); parsSig.add(meanGSig3); parsSig.add(sigmaGSig3);
  parsSig.add(mpvSig); parsSig.add(sigmaLSig); parsSig.add(mpvSig2); parsSig.add(sigmaLSig2); parsSig.add(mpvSig3); parsSig.add(sigmaLSig3);
  parsSig.add(lambdaExpSig); parsSig.add(fSigG1); parsSig.add(fSigG2); parsSig.add(fSig3); parsSig.add(fSig4); parsSig.add(fSig5); parsSig.add(fSig6);

  if (doErrPdfFit)
  {
    gSystem->mkdir(resultDir, true);
    TFile foutErr(errPdfFile, "RECREATE");
    parsBkg.Write("errParsBkg");
    parsSig.Write("errParsSig");
    TParameter<int>("nErrBkgGauss", nErrBkgGauss).Write();
    TParameter<int>("nErrBkgLandau", nErrBkgLandau).Write();
    TParameter<int>("nErrBkgExp", nErrBkgExp).Write();
    TParameter<int>("nErrSigGauss", nErrSigGauss).Write();
    TParameter<int>("nErrSigLandau", nErrSigLandau).Write();
    TParameter<int>("nErrSigExp", nErrSigExp).Write();
    TParameter<double>("errLow", errLow).Write();
    TParameter<double>("errHigh", errHigh).Write();
    if (fitErrBkgPtr)
      fitErrBkgPtr->Write("fitErrBkg");
    if (fitErrSigPtr)
      fitErrSigPtr->Write("fitErrSig");
    foutErr.Close();
    printf("[OK] saved errPdf params: %s\n", errPdfFile.Data());
  }
  else
  {
    if (!load_errpdf_params(errPdfFile, parsBkg, parsSig, fitErrBkgPtr, fitErrSigPtr))
      return;
  }
  meanGBkg.setConstant(true);
  sigmaGBkg.setConstant(true);
  meanGBkg2.setConstant(true);
  sigmaGBkg2.setConstant(true);
  meanGBkg3.setConstant(true);
  sigmaGBkg3.setConstant(true);
  mpvBkg.setConstant(true);
  sigmaLBkg.setConstant(true);
  mpvBkg2.setConstant(true);
  sigmaLBkg2.setConstant(true);
  mpvBkg3.setConstant(true);
  sigmaLBkg3.setConstant(true);
  lambdaExpBkg.setConstant(true);
  fBkgG1.setConstant(true);
  fBkgG2.setConstant(true);
  fBkg3.setConstant(true);
  fBkg4.setConstant(true);
  fBkg5.setConstant(true);
  fBkg6.setConstant(true);
  meanGSig.setConstant(true);
  sigmaGSig.setConstant(true);
  meanGSig2.setConstant(true);
  sigmaGSig2.setConstant(true);
  meanGSig3.setConstant(true);
  sigmaGSig3.setConstant(true);
  mpvSig.setConstant(true);
  sigmaLSig.setConstant(true);
  mpvSig2.setConstant(true);
  sigmaLSig2.setConstant(true);
  mpvSig3.setConstant(true);
  sigmaLSig3.setConstant(true);
  lambdaExpSig.setConstant(true);
  fSigG1.setConstant(true);
  fSigG2.setConstant(true);
  fSig3.setConstant(true);
  fSig4.setConstant(true);
  fSig5.setConstant(true);
  fSig6.setConstant(true);
  printf("[OK] errPdf ready (bkg: G%d+L%d+E%d, sig: G%d+L%d+E%d)\n",
         nErrBkgGauss, nErrBkgLandau, nErrBkgExp, nErrSigGauss, nErrSigLandau, nErrSigExp);

  // ------------------------------------------------------------------
  // plot errPdf / errPdfSig with the datasets used to build them
  // ------------------------------------------------------------------
  if (doErrPdfFit)
  {
    const double nErrBkg = hErrBkg->sumEntries();
    const double nErrSig = hErrSig->sumEntries();

    auto findObj = [&](RooPlot *fr, const char *n) -> TObject *
    { return fr ? fr->findObject(n) : nullptr; };
    auto plotComp = [&](RooAbsPdf *pdf, RooPlot *fr, const char *compName, const char *objName, Color_t col, Style_t style)
    {
      if (!pdf || !fr)
        return;
      pdf->plotOn(fr, Components(compName), Name(objName), LineColor(col), LineStyle(style));
    };

    RooPlot *frErrBkg = ctau3DErr->frame(Title(""));
    hErrBkg->plotOn(frErrBkg, DataError(RooAbsData::SumW2), Name("data_bkg"));
    pdfErr->plotOn(frErrBkg, Name("pdf_bkg"), LineColor(kBlack), LineWidth(2));
    if (nErrBkgLandau >= 1) plotComp(pdfErr, frErrBkg, "landauBkg", "pdf_bkg_landau1", kBlue + 2, kDashed);
    if (nErrBkgLandau >= 2) plotComp(pdfErr, frErrBkg, "landauBkg2", "pdf_bkg_landau2", kOrange + 7, kDashed);
    if (nErrBkgLandau >= 3) plotComp(pdfErr, frErrBkg, "landauBkg3", "pdf_bkg_landau3", kMagenta + 1, kDashed);
    if (nErrBkgGauss >= 1) plotComp(pdfErr, frErrBkg, "gausBkg", "pdf_bkg_gaus1", kGreen + 2, kDashed);
    if (nErrBkgGauss >= 2) plotComp(pdfErr, frErrBkg, "gausBkg2", "pdf_bkg_gaus2", kRed + 1, kDashed);
    if (nErrBkgGauss >= 3) plotComp(pdfErr, frErrBkg, "gausBkg3", "pdf_bkg_gaus3", kAzure + 2, kDashed);
    if (nErrBkgExp >= 1) plotComp(pdfErr, frErrBkg, "expBkg", "pdf_bkg_exp1", kViolet + 1, kDashed);
    frErrBkg->GetYaxis()->SetTitle("Events");
    frErrBkg->GetYaxis()->SetTitleOffset(1.6);
    frErrBkg->GetXaxis()->SetTitle("");
    double errBkgYMax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frErrBkg->getHist("data_bkg")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double xval = 0.0, yval = 0.0;
        hdata->GetPoint(i, xval, yval);
        if (yval > errBkgYMax)
          errBkgYMax = yval;
      }
    }
    if (errBkgYMax > 0.0 && errBkgYMax < 1e300)
      frErrBkg->SetMaximum(errBkgYMax * 1.8);

    RooPlot *frErrSig = ctau3DErr->frame(Title(""));
    hErrSig->plotOn(frErrSig, DataError(RooAbsData::SumW2), Name("data_sig"));
    pdfErrSig->plotOn(frErrSig, Name("pdf_sig"), LineColor(kBlack), LineWidth(2));
    if (nErrSigLandau >= 1) plotComp(pdfErrSig, frErrSig, "landauSig", "pdf_sig_landau1", kBlue + 2, kDashed);
    if (nErrSigLandau >= 2) plotComp(pdfErrSig, frErrSig, "landauSig2", "pdf_sig_landau2", kOrange + 7, kDashed);
    if (nErrSigLandau >= 3) plotComp(pdfErrSig, frErrSig, "landauSig3", "pdf_sig_landau3", kMagenta + 1, kDashed);
    if (nErrSigGauss >= 1) plotComp(pdfErrSig, frErrSig, "gausSig", "pdf_sig_gaus1", kGreen + 2, kDashed);
    if (nErrSigGauss >= 2) plotComp(pdfErrSig, frErrSig, "gausSig2", "pdf_sig_gaus2", kRed + 1, kDashed);
    if (nErrSigGauss >= 3) plotComp(pdfErrSig, frErrSig, "gausSig3", "pdf_sig_gaus3", kAzure + 2, kDashed);
    if (nErrSigExp >= 1) plotComp(pdfErrSig, frErrSig, "expSig", "pdf_sig_exp1", kViolet + 1, kDashed);
    frErrSig->GetYaxis()->SetTitle("Events");
    frErrSig->GetYaxis()->SetTitleOffset(1.6);
    frErrSig->GetXaxis()->SetTitle("");
    double errSigYMax = -1e300;
    if (auto *hdata = dynamic_cast<RooHist *>(frErrSig->getHist("data_sig")))
    {
      for (int i = 0; i < hdata->GetN(); ++i)
      {
        double xval = 0.0, yval = 0.0;
        hdata->GetPoint(i, xval, yval);
        if (yval > errSigYMax)
          errSigYMax = yval;
      }
    }
    if (errSigYMax > 0.0 && errSigYMax < 1e300)
      frErrSig->SetMaximum(errSigYMax * 1.8);

    gSystem->mkdir(figDir, true);
    TCanvas cErrBkg("cErrBkg", "", 800, 800);
    TPad pad1B("pad1B", "pad1B", 0.0, 0.25, 1.0, 1.0);
    pad1B.SetBottomMargin(0.00001);
    pad1B.SetTopMargin(0.08);
    pad1B.Draw();
    pad1B.cd();
    frErrBkg->Draw("e");
    {
      auto *leg = new TLegend(0.50, 0.66, 0.74, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.03);
      if (auto *o = findObj(frErrBkg, "data_bkg"))
        leg->AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frErrBkg, "pdf_bkg"))
        leg->AddEntry(o, "Fit", "l");
      if (auto *o = findObj(frErrBkg, "pdf_bkg_gaus1"))
        leg->AddEntry(o, "Gauss 1", "l");
      if (nErrBkgGauss >= 2)
        if (auto *o = findObj(frErrBkg, "pdf_bkg_gaus2"))
        leg->AddEntry(o, "Gauss 2", "l");
      if (nErrBkgGauss >= 3)
        if (auto *o = findObj(frErrBkg, "pdf_bkg_gaus3"))
        leg->AddEntry(o, "Gauss 3", "l");
      if (auto *o = findObj(frErrBkg, "pdf_bkg_landau1"))
        leg->AddEntry(o, "Landau tail", "l");
      if (nErrBkgLandau >= 2)
        if (auto *o = findObj(frErrBkg, "pdf_bkg_landau2"))
        leg->AddEntry(o, "Landau 2", "l");
      if (nErrBkgLandau >= 3)
        if (auto *o = findObj(frErrBkg, "pdf_bkg_landau3"))
        leg->AddEntry(o, "Landau 3", "l");
      if (nErrBkgExp >= 1)
        if (auto *o = findObj(frErrBkg, "pdf_bkg_exp1"))
        leg->AddEntry(o, "Exp", "l");
      leg->Draw("same");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.032);
      tx.SetTextFont(42);
      tx.SetTextAlign(31);
      tx.DrawLatex(0.96, 0.935, "OO #sqrt{s_{NN}} = 5.36 TeV (9 nb^{-1})");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.04);
      tx.SetTextFont(72);
      tx.DrawLatex(0.19, 0.935, "CMS Internal");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double xtext = 0.19, y0 = 0.865, dy = -0.05;
      int k = 0;
      tx.DrawLatex(xtext, y0 + dy * k++, "Mass sideband data");
      if (yLow == 0)
        tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
    }
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double xtext = 0.75, y0 = 0.87, dy = -0.04;
      int k = 0;
      auto print = [&](const char *title, const RooRealVar &v)
      {
        if (v.isConstant())
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, v.getVal()));
        else
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, v.getVal(), v.getError()));
      };
      print("fG1", fBkgG1);
      print("fG2", fBkgG2);
      print("meanG1", meanGBkg);
      print("sigmaG1", sigmaGBkg);
      print("meanG2", meanGBkg2);
      print("sigmaG2", sigmaGBkg2);
      print("mpvL", mpvBkg);
      print("sigmaL", sigmaLBkg);
    }
    cErrBkg.cd();
    TPad pad2B("pad2B", "pad2B", 0.0, 0.0, 1.0, 0.25);
    pad2B.SetTopMargin(0.00001);
    pad2B.SetBottomMargin(0.4);
    pad2B.Draw();
    pad2B.cd();
    RooPlot *frPullB = ctau3DErr->frame(Title(""));
    RooHist *hpullB = frErrBkg->pullHist("data_bkg", "pdf_bkg");
    if (hpullB)
      frPullB->addPlotable(hpullB, "P");
    frPullB->GetYaxis()->SetTitle("Pull");
    frPullB->GetXaxis()->SetTitle("c#tau_{3D} err [mm]");
    frPullB->GetXaxis()->CenterTitle();
    frPullB->SetMinimum(-8);
    frPullB->SetMaximum(8);
    frPullB->GetYaxis()->SetNdivisions(505);
    frPullB->GetYaxis()->SetTitleSize(0.12);
    frPullB->GetYaxis()->SetLabelSize(0.10);
    frPullB->GetXaxis()->SetTitleSize(0.15);
    frPullB->GetXaxis()->SetLabelSize(0.10);
    frPullB->Draw();
    {
      auto chiM = chi2_from_pull(hpullB);
      RooArgSet obsErr(*ctau3DErr);
      const int npar = count_float_pars(*pdfErr, obsErr, fitErrBkgPtr);
      const int ndf = std::max(1, chiM.second - npar);
      const double pvalue = TMath::Prob(chiM.first, ndf);
      printf("[PULL][BKG] chi2=%.3f npts=%d npar=%d ndf=%d chi2/ndf=%.4f\n",
             chiM.first, chiM.second, npar, ndf, chiM.first / ndf);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.085);
      tc.SetTextFont(42);
      tc.SetTextAlign(33);
      tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", chiM.first, ndf, pvalue));
    }
    cErrBkg.SaveAs(figName("errpdf_bkg"));

    TCanvas cErrSig("cErrSig", "", 800, 800);
    TPad pad1S("pad1S", "pad1S", 0.0, 0.25, 1.0, 1.0);
    pad1S.SetBottomMargin(0.00001);
    pad1S.SetTopMargin(0.08);
    pad1S.Draw();
    pad1S.cd();
    frErrSig->Draw("e");
    {
      auto *leg = new TLegend(0.50, 0.66, 0.74, 0.89);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.03);
      if (auto *o = findObj(frErrSig, "data_sig"))
        leg->AddEntry(o, "Data", "lep");
      if (auto *o = findObj(frErrSig, "pdf_sig"))
        leg->AddEntry(o, "Fit", "l");
      if (auto *o = findObj(frErrSig, "pdf_sig_gaus1"))
        leg->AddEntry(o, "Gauss 1", "l");
      if (nErrSigGauss >= 2)
        if (auto *o = findObj(frErrSig, "pdf_sig_gaus2"))
        leg->AddEntry(o, "Gauss 2", "l");
      if (nErrSigGauss >= 3)
        if (auto *o = findObj(frErrSig, "pdf_sig_gaus3"))
        leg->AddEntry(o, "Gauss 3", "l");
      if (auto *o = findObj(frErrSig, "pdf_sig_landau1"))
        leg->AddEntry(o, "Landau tail", "l");
      if (nErrSigLandau >= 2)
        if (auto *o = findObj(frErrSig, "pdf_sig_landau2"))
        leg->AddEntry(o, "Landau 2", "l");
      if (nErrSigLandau >= 3)
        if (auto *o = findObj(frErrSig, "pdf_sig_landau3"))
        leg->AddEntry(o, "Landau 3", "l");
      if (nErrSigExp >= 1)
        if (auto *o = findObj(frErrSig, "pdf_sig_exp1"))
        leg->AddEntry(o, "Exp", "l");
      leg->Draw("same");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.032);
      tx.SetTextFont(42);
      tx.SetTextAlign(31);
      tx.DrawLatex(0.96, 0.935, "OO #sqrt{s_{NN}} = 5.36 TeV (9 nb^{-1})");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.04);
      tx.SetTextFont(72);
      tx.DrawLatex(0.19, 0.935, "CMS Internal");
    }
    {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.03);
      tx.SetTextFont(42);
      double xtext = 0.19, y0 = 0.865, dy = -0.05;
      int k = 0;
      tx.DrawLatex(xtext, y0 + dy * k++, "Signal region - scaled sideband");
      if (yLow == 0)
        tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, |y| < %.1f", ptLow, ptHigh, yHigh));
      else
        tx.DrawLatex(xtext, y0 + dy * k++, Form("%.1f < p_{T} < %.1f, %.1f < |y| < %.1f", ptLow, ptHigh, yLow, yHigh));
    }
    {
      TLatex tp;
      tp.SetNDC();
      tp.SetTextSize(0.024);
      tp.SetTextFont(42);
      double xtext = 0.75, y0 = 0.87, dy = -0.04;
      int k = 0;
      auto print = [&](const char *title, const RooRealVar &v)
      {
        if (v.isConstant())
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g (fixed)", title, v.getVal()));
        else
          tp.DrawLatex(xtext, y0 + dy * k++, Form("%s = %.4g #pm %.3g", title, v.getVal(), v.getError()));
      };
      print("fG1", fSigG1);
      print("fG2", fSigG2);
      print("meanG", meanGSig);
      print("sigmaG", sigmaGSig);
      print("mpvL", mpvSig);
      print("sigmaL", sigmaLSig);
    }
    cErrSig.cd();
    TPad pad2S("pad2S", "pad2S", 0.0, 0.0, 1.0, 0.25);
    pad2S.SetTopMargin(0.00001);
    pad2S.SetBottomMargin(0.4);
    pad2S.Draw();
    pad2S.cd();
    RooPlot *frPullS = ctau3DErr->frame(Title(""));
    RooHist *hpullS = frErrSig->pullHist("data_sig", "pdf_sig");
    if (hpullS)
      frPullS->addPlotable(hpullS, "P");
    frPullS->GetYaxis()->SetTitle("Pull");
    frPullS->GetXaxis()->SetTitle("c#tau_{3D} err [mm]");
    frPullS->GetXaxis()->CenterTitle();
    frPullS->SetMinimum(-8);
    frPullS->SetMaximum(8);
    frPullS->GetYaxis()->SetNdivisions(505);
    frPullS->GetYaxis()->SetTitleSize(0.12);
    frPullS->GetYaxis()->SetLabelSize(0.10);
    frPullS->GetXaxis()->SetTitleSize(0.15);
    frPullS->GetXaxis()->SetLabelSize(0.10);
    frPullS->Draw();
    {
      auto chiM = chi2_from_pull(hpullS);
      RooArgSet obsErr(*ctau3DErr);
      const int npar = count_float_pars(*pdfErrSig, obsErr, fitErrSigPtr);
      const int ndf = std::max(1, chiM.second - npar);
      const double pvalue = TMath::Prob(chiM.first, ndf);
      printf("[PULL][SIG] chi2=%.3f npts=%d npar=%d ndf=%d chi2/ndf=%.4f\n",
             chiM.first, chiM.second, npar, ndf, chiM.first / ndf);
      TLatex tc;
      tc.SetNDC();
      tc.SetTextSize(0.085);
      tc.SetTextFont(42);
      tc.SetTextAlign(33);
      tc.DrawLatex(0.88, 0.96, Form("#chi^{2}/ndf = %.1f/%d (%.3g)", chiM.first, ndf, pvalue));
    }
    cErrSig.SaveAs(figName("errpdf_sig"));
  }

  timer.Stop();
  char realBuf[32];
  char cpuBuf[32];
  auto formatTime = [](double sec, char *buf, size_t len)
  {
    int total = static_cast<int>(sec + 0.5);
    int h = total / 3600;
    int m = (total % 3600) / 60;
    int s = total % 60;
    if (h > 0)
      snprintf(buf, len, "%dh %dm %ds", h, m, s);
    else
      snprintf(buf, len, "%dm %ds", m, s);
  };
  formatTime(timer.RealTime(), realBuf, sizeof(realBuf));
  formatTime(timer.CpuTime(), cpuBuf, sizeof(cpuBuf));
  // printf("[OK] macro time: real %s (cpu %d)\n", realBuf, number_of_workers);
  printf("[OK] macro time: real %s\n", realBuf);
  auto extract_status_hesse = [](RooFitResult *fr, int &status, int &hesse)
  {
    status = -1;
    hesse = -1;
    if (!fr)
      return;
    status = fr->status();
    for (UInt_t i = 0, n = fr->numStatusHistory(); i < n; ++i)
    {
      const char *lab = fr->statusLabelHistory(i);
      if (lab && !strcmp(lab, "HESSE"))
      {
        hesse = fr->statusCodeHistory(i);
        break;
      }
    }
  };
  int bkgStatus, bkgHesse, sigStatus, sigHesse;
  extract_status_hesse(fitErrBkgPtr, bkgStatus, bkgHesse);
  extract_status_hesse(fitErrSigPtr, sigStatus, sigHesse);
  const bool bkgConv = (bkgStatus == 0 && bkgHesse == 0);
  const bool sigConv = (sigStatus == 0 && sigHesse == 0);
  printf("[SUMMARY] errPdf converged=%d (bkg=%d, sig=%d) bkg(MINIMIZE=%d,HESSE=%d) sig(MINIMIZE=%d,HESSE=%d)\n",
         (bkgConv && sigConv) ? 1 : 0, bkgConv ? 1 : 0, sigConv ? 1 : 0,
         bkgStatus, bkgHesse, sigStatus, sigHesse);
}

// ------------------------------------------------------------------
//
// ------------------------------------------------------------------

// ------------------------------------------------------------------
//
// ------------------------------------------------------------------
