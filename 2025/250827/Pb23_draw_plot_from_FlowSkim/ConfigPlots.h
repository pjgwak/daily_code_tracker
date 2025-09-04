// list of FlowSkim plots

#pragma once
#include "TCut.h"
#include <string>
#include <vector>
#include <algorithm>

// ----- preview options -----
// turn it on when you want to test new functions or rangeds
static const bool PREVIEW_ON = true;
static const Long64_t PREVIEW_NENTRIES = 100000; // how many entries to use?
static const Long64_t PREVIEW_FIRSTENTRY = 0;


// ----- plot structure -----
struct Plot1D {
  std::string name; // default histogram name
  std::string expr; // branch name, e.g. "mass"
  int nbins;
  double xmin, xmax;
  std::string title;
  TCut baseCut;
  std::string sliceTag; // choose slice list
};

struct Plot2D
{
  std::string name;
  std::string exprYX; // e.g. "mass:y"
  int nxbins; double xmin, xmax;
  int nybins; double ymin, ymax;
  std::string title; // ";X title;Y title"
  TCut baseCut;
  bool save_profileX = false;
  bool save_profileY = false;
  std::string sliceTag;
};

struct Slice {
  std::string suffix; // same bracnh but differnt cuts -> to handle many pT bins
  TCut cut;
};

//  ----- IO paths -----
static const std::string INPUT_FILE = "/data/users/pjgwak/work/daily_code_tracker/2025/250826/Pb23_flowSkim/FlowSkim_2023PbPbPromptRecoData_132X_miniAOD_noTrack_250826.root";
static const std::string TREE_NAME = "myTree";
static const std::string OUT_DIR = "figs";


// ===== cuts =====
// ----- some basic cuts -----
static const TCut CUT_JPSI_MASS = "mass>2.6 && mass<3.5"; // It's applied at draw_from_config.C
static const TCut CUT_FWD = "pt>3 && fabs(y)>1.6 && fabs(y)<2.4";
static const TCut CUT_MID = "pt>6.5 && fabs(y)<1.6";
static const TCut CUT_ALL = "pt>6.5 && fabs(y)<2.4";

// ----- ctau cuts -----
// ctau PR and NP
// PR: |ctau3D| < 0.05
// NP:  0.1 <= ctau3D <= 0.8
static const TCut CUT_PR_CTAU = "fabs(ctau3D) < 0.05";
static const TCut CUT_NP_CTAU = "ctau3D >= 0.1 && ctau3D <= 0.8";

// mass region (SR/LSB/RSB)
// SR : 3.00 ~ 3.20
// LSB: 2.60 ~ 2.95
// RSB: 3.21 ~ 3.50
static const TCut CUT_SR_MASS = "mass > 3.0  && mass < 3.2";
static const TCut CUT_LSB_MASS = "mass > 2.6  && mass < 2.95";
static const TCut CUT_RSB_MASS = "mass > 3.21 && mass < 3.5";

// make mass-ctau slices
inline std::vector<Slice> BuildRegion6Slices()
{
  return {
      {"_PR_SR", CUT_PR_CTAU && CUT_SR_MASS},
      {"_PR_LSB", CUT_PR_CTAU && CUT_LSB_MASS},
      {"_PR_RSB", CUT_PR_CTAU && CUT_RSB_MASS},
      {"_NP_SR", CUT_NP_CTAU && CUT_SR_MASS},
      {"_NP_LSB", CUT_NP_CTAU && CUT_LSB_MASS},
      {"_NP_RSB", CUT_NP_CTAU && CUT_RSB_MASS},
  };
}

// ----- helper functions -----
inline void trimTrailingZeros(std::string &s)
{
  const auto pos = s.find('.');
  if (pos == std::string::npos)
    return;
  while (!s.empty() && s.back() == '0')
    s.pop_back();
  if (!s.empty() && s.back() == '.')
    s.pop_back();
}

inline std::string toStr(double x)
{
  std::string s = std::to_string(x); // "6.500000"
  trimTrailingZeros(s);              // "6.5"
  return s;
}

inline std::string tagFromDouble(double x)
{
  std::string s = toStr(x);                   // "6.5"
  std::replace(s.begin(), s.end(), '.', 'p'); // "6p5"
  return s;
}

inline std::vector<Slice> BuildPtSlicesFromEdges(const std::vector<double> &edges)
{
  std::vector<Slice> out;
  if (edges.size() < 2)
    return out;
  for (size_t i = 0; i + 1 < edges.size(); ++i)
  {
    const double lo = edges[i], hi = edges[i + 1];
    const std::string loTag = tagFromDouble(lo);
    const std::string hiTag = tagFromDouble(hi);
    Slice s;
    s.suffix = std::string("_pt") + loTag + "to" + hiTag;
    s.cut = TCut((std::string("pt>") + toStr(lo) + " && pt<" + toStr(hi)).c_str());
    out.push_back(s);
  }
  return out;
}

inline std::vector<Slice> BuildRapSlices()
{
  return {
      {"_fwd", CUT_FWD},
      {"_mid", CUT_MID},
      {"_all", CUT_ALL}};
}

// mass: pT (fwd, mid+All)
inline const std::vector<double> &EdgesPtMassFwd()
{
  static const std::vector<double> e = {3, 6.5, 9, 12, 15, 20, 30, 50};
  return e;
}
inline const std::vector<double> &EdgesPtMassMidAll()
{
  static const std::vector<double> e = {6.5, 9, 12, 15, 20, 30, 50};
  return e;
}

// mass: rap × pT
inline std::vector<Slice> BuildRapxPtMassSlices()
{
  const auto rap = BuildRapSlices();
  std::vector<Slice> out;
  for (const auto &r : rap)
  {
    const bool isFwd = (r.suffix == "_fwd");
    const auto ptSlices = BuildPtSlicesFromEdges(isFwd ? EdgesPtMassFwd() : EdgesPtMassMidAll());
    for (const auto &p : ptSlices)
    {
      out.push_back({r.suffix + p.suffix, r.cut && p.cut});
    }
  }
  return out;
}

inline std::vector<Slice> BuildRegion6_RapxPtMassSlices()
{
  const auto reg6 = BuildRegion6Slices();
  const auto rap = BuildRapSlices();
  std::vector<Slice> out;
  out.reserve(reg6.size() * rap.size() * 16);

  for (const auto &r6 : reg6)
  {
    for (const auto &r : rap)
    {
      const bool isFwd = (r.suffix == "_fwd");
      const auto ptSlices = BuildPtSlicesFromEdges(isFwd ? EdgesPtMassFwd()
                                                         : EdgesPtMassMidAll());
      for (const auto &p : ptSlices)
      {
        // eg: _NP_LSB + _fwd + _pt3to6p5 -> _NP_LSB_fwd_pt3to6p5
        out.push_back({r6.suffix + r.suffix + p.suffix, r6.cut && r.cut && p.cut});
      }
    }
  }
  return out;
}

inline std::vector<Slice> BuildRapSlices_MassPtU50()
{
  return {
      {"_fwd_pt3to50", CUT_FWD && "pt<=50"},
      {"_mid_pt6p5to50", CUT_MID && "pt<=50"},
      {"_all_pt6p5to50", CUT_ALL && "pt<=50"}};
}

inline std::vector<Slice> SlicesForTag(const std::string &tag)
{
  if (tag=="NONE")  return { {"",TCut()} };
  if (tag=="RAP") return BuildRapSlices();
  if (tag=="RAPxPT_MASS") return BuildRapxPtMassSlices();
  if (tag=="RAP_MASS1D_PTUP50") return BuildRapSlices_MassPtU50();
  if (tag == "REGION6") return BuildRegion6Slices();
    if (tag=="REGION6xRAPxPT_MASS") return BuildRegion6_RapxPtMassSlices();
  return BuildRapSlices(); // default: RAP
}

// ----- plot list -----
// 1D
static const std::vector<Plot1D> CFG_1D = {
    // y (no cut)
    {"y_nocut", "y", 48, -2.4, 2.4, "rapidity;y;Events", TCut(), "NONE"},

    // y — rapidity
    {"y_by_rap", "y", 48, -2.4, 2.4, "rapidity;y;Events", TCut(), "RAP"},

    // pT (no cut)
    {"pt_nocut", "pt", 50, 0, 50, "p_{T};p_{T} [GeV];Events", TCut(), "NONE"},

    // pT — rapidity
    {"pt_by_rap", "pt", 50, 0, 50, "p_{T};p_{T} [GeV];Events", TCut(), "RAP"},

    // mass — rapidity
    {"mass_by_rap", "mass", 100, 2.6, 3.5, "mass;m_{#mu#mu} [GeV];Events", TCut(), "RAP"},

    // mass - integrated pT
    {"mass_by_rap_ptU50", "mass", 100, 2.6, 3.5, "mass;m_{#mu#mu} [GeV];Events", TCut(), "RAP_MASS1D_PTUP50"},

    // mass — rapidty + pT regions
    {"mass_by_ptbins_rap", "mass", 100, 2.6, 3.5, "mass;m_{#mu#mu} [GeV];Events", TCut(), "RAPxPT_MASS"},

    // PR, NP region: mass -> Don't use it. We need to distinguish by rapidity at least
    // {"mass_region6_only", "mass", 100, 2.6, 3.5, "mass;m_{#mu#mu} [GeV];Events", TCut(), "REGION6"},

    // PR, NP region: ctau3D -> Don't use it. We need to distinguish by rapidity at least
    // {"ctau3D_region6_only", "ctau3D", 180, -0.1, 0.8, "ctau3D;ctau3D;Events", TCut(), "REGION6"},

    // region6 × RAP × pT: mass
    {"mass_region6", "mass", 100, 2.6, 3.5, "mass;m_{#mu#mu} [GeV];Events", TCut(), "REGION6xRAPxPT_MASS"},

    // region6 × RAP × pT: ctau3D
    {"ctau3D_region6", "ctau3D", 180, -0.1, 0.8, "ctau3D;ctau3D;Events", TCut(), "REGION6xRAPxPT_MASS"},
};

// 2D
static const std::vector<Plot2D> CFG_2D = {
    // pT vs y (no cut)
    {"pt_vs_y_nocut", "pt:y", 48, -2.4, 2.4, 50, 0, 50, ";y;p_{T} [GeV]", TCut(), false, false, "NONE"},
    // pT vs y — rapidity
    {"pt_vs_y_by_rap", "pt:y", 48, -2.4, 2.4, 50, 0, 50, ";y;p_{T} [GeV]", TCut(), false, false, "RAP"},
    // pT va mass — rapidity
    {"pt_vs_mass_by_rap", "pt:mass", 100, 2.6, 3.5, 50, 0, 50, ";m_{#mu#mu} [GeV];p_{T} [GeV]", TCut(), false, false, "RAP"},

    // ctau3D vs mass - rapidity
    {"ctau3D_vs_mass_by_rap", "ctau3D:mass", 100, 2.6, 3.5, 180, -0.1, 0.8,
     ";m_{#mu#mu} [GeV];ctau3D", TCut(), false, false, "RAP"},
};

// ----- 1D overlay plot (for rapidity) -----
struct Overlay1D {
  std::string name;
  std::string expr; // "y"
  int nbins; double xmin, xmax;
  std::string title;
  std::vector<std::pair<std::string, TCut>> series; // labels, cut
  bool normalize = true;
};

static const std::vector<Overlay1D> CFG_OVERLAY_1D = {
    {"y_overlay_rap", "y", 48, -2.4, 2.4, "rapidity;y;Events", {{"fwd", CUT_FWD}, {"mid", CUT_MID}, {"all", CUT_ALL}}}
};