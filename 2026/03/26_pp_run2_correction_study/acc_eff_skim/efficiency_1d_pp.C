// ===== Weight summary =====
// GEN denominator: NColl * TreeWeight * PtWeight
// Reco numerator: NColl * TreeWeight * PtWeight * TnP
// Current status:
//   - PbPb TnP weight follows makeMuonSkimTree_Jpsi.C logic.


/*
# pp PR Reco - for Eff
/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root

# pp NP Reco - for Eff
/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root
const std::string ptWeightMidPath = isPr
          ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251125.root"
          : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251125.root";
      const std::string ptWeightFwdPath = isPr
          ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251125.root"
          : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251125.root";

// PbPb
/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251118.root
/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root

/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251118.root
/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root

*/

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/tnp_weight_lowptPbPb_num_den_new.h"
#include "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/tnp_weight_pp.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
// #include "cutsAndBins.h"
using std::cout; using std::string;

static const long MAXTREESIZE = 1000000000000;

// 2018
bool IsAcceptanceQQ(double pt, double eta)
{
  return ((fabs(eta) < 0.3 && pt > 3.4) ||
          (0.3 < fabs(eta) && fabs(eta) < 1.1 && pt > 3.3) ||
          (1.1 < fabs(eta) && fabs(eta) < 1.5 && pt > 9.08 - 5.25*fabs(eta)) ||
          (1.5 < fabs(eta) && fabs(eta) < 2.4 && pt > 0.8 && pt>2.4-0.8*fabs(eta)));
}

// Run3 - Eff numerator and standard of skim
// bool IsAcceptanceQQ(double pt, double eta)
// {
//   const double aeta = std::abs(eta);

//   if (aeta > 2.4)
//     return false;
//   if (aeta <= 1.0)
//     return pt >= 3.3;
//   if (aeta <= 1.3)
//     return pt >= ((2.1 - 3.3) / (1.3 - 1.0)) * (aeta - 1.0) + 3.3;
//   if (aeta <= 1.7)
//     return pt >= ((1.0 - 2.1) / (1.7 - 1.3)) * (aeta - 1.3) + 2.1;

//   return pt >= 1.0;
// }

// Run3 OO MuAccCut - from Jun Chen
// Acc Numerator and Eff denominator
bool IsAcceptanceQQJunChenOO(double pt, double eta)
{
  const double aeta = std::abs(eta);
  const double p = pt * std::cosh(eta);

  return (pt >= 0.5 && p >= 2.5 && aeta <= 2.5);
}

bool IsGlobalMuonInAccept2015(double pt, double eta)
{
  const double aeta = std::abs(eta);

  return (aeta < 2.4 &&
          ((aeta < 1.2 && pt >= 3.5) ||
           (1.2 <= aeta && aeta < 2.1 && pt >= 5.77 - 1.89 * aeta) ||
           (2.1 <= aeta && pt >= 1.8)));
}

double findNcoll(int hiBin)
{
  static const int nbins = 200;
  static const Double_t Ncoll[nbins] = {
      1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3,
      1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98,
      1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96,
      980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163,
      765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456,
      593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372,
      451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239,
      340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625,
      249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644,
      180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365,
      127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729,
      87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859,
      58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988,
      38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965,
      24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323,
      14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648,
      8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056,
      5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282,
      2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916,
      1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};

  if (hiBin < 0)
    return Ncoll[0];
  if (hiBin >= nbins)
    return Ncoll[nbins - 1];
  return Ncoll[hiBin];
}

double GetPtWeight(TF1 *fW, double pt, double y)
{
  (void)y;
  if (!fW)
    return 1.0;
  return fW->Eval(pt);
}

bool PassMuonSoftCutFallback(Int_t selectionType, Int_t nTrkWMea, Int_t nPixWMea,
                             Float_t dxy, Float_t dz)
{
  const bool passMuonType =
      ((selectionType & (1 << 1)) != 0) &&
      ((selectionType & (1 << 3)) != 0);

  return (nTrkWMea > 5) &&
         (nPixWMea > 0) &&
         (std::fabs(dxy) < 0.3f) &&
         (std::fabs(dz) < 20.f) &&
         passMuonType;
}

bool PassMuonSelectionTypeLikeSkim(Int_t selectionType)
{
  const bool passBit1 = ((selectionType & (1 << 1)) != 0);
  const bool passBit3 = ((selectionType & (1 << 3)) != 0);
  return passBit1 && passBit3;
}

bool ComputePbPbTnPWeight(double muPlPt, double muPlEta, ULong64_t muPlTrig,
                          double muMiPt, double muMiEta, ULong64_t muMiTrig,
                          double &tnpWeight)
{
  // Match the nominal PbPb TnP treatment used in makeMuonSkimTree_Jpsi.C.
  const int tnp_ids = 0;
  const int tnp_trks = 0;
  const int tnp_trgs = 0;
  const int kL2filter = 16;
  const int kL3filter = 17;

  tnpWeight = 1.0;
  double tnpTrigWeightMuPl = -1.0;
  double tnpTrigWeightMuMi = -1.0;

  tnpWeight *= tnp_weight_muid_pbpb(muPlPt, muPlEta, tnp_ids);
  tnpWeight *= tnp_weight_muid_pbpb(muMiPt, muMiEta, tnp_ids);
  tnpWeight *= tnp_weight_trk_pbpb(muPlEta, tnp_trks);
  tnpWeight *= tnp_weight_trk_pbpb(muMiEta, tnp_trks);

  const ULong64_t l2Mask = (ULong64_t(1) << kL2filter);
  const ULong64_t l3Mask = (ULong64_t(1) << kL3filter);
  const bool muPlL2 = (muPlTrig & l2Mask) == l2Mask;
  const bool muPlL3 = (muPlTrig & l3Mask) == l3Mask;
  const bool muMiL2 = (muMiTrig & l2Mask) == l2Mask;
  const bool muMiL3 = (muMiTrig & l3Mask) == l3Mask;

  if (!muPlL2 || !muMiL2)
    return false;

  const bool muPlIsL2 = muPlL2 && !muPlL3;
  const bool muPlIsL3 = muPlL2 && muPlL3;
  const bool muMiIsL2 = muMiL2 && !muMiL3;
  const bool muMiIsL3 = muMiL2 && muMiL3;

  if (muPlIsL2 && muMiIsL3)
  {
    tnpTrigWeightMuPl = tnp_weight_trg_pbpb(muPlPt, muPlEta, 2, tnp_trgs);
    tnpTrigWeightMuMi = tnp_weight_trg_pbpb(muMiPt, muMiEta, 3, tnp_trgs);
  }
  else if (muPlIsL3 && muMiIsL2)
  {
    tnpTrigWeightMuPl = tnp_weight_trg_pbpb(muPlPt, muPlEta, 3, tnp_trgs);
    tnpTrigWeightMuMi = tnp_weight_trg_pbpb(muMiPt, muMiEta, 2, tnp_trgs);
  }
  else if (muPlIsL3 && muMiIsL3)
  {
    int t[2] = {-1, 1};
    int l = rand() % 2;
    if (t[l] == -1)
    {
      tnpTrigWeightMuPl = tnp_weight_trg_pbpb(muPlPt, muPlEta, 2, tnp_trgs);
      tnpTrigWeightMuMi = tnp_weight_trg_pbpb(muMiPt, muMiEta, 3, tnp_trgs);
    }
    else
    {
      tnpTrigWeightMuPl = tnp_weight_trg_pbpb(muPlPt, muPlEta, 3, tnp_trgs);
      tnpTrigWeightMuMi = tnp_weight_trg_pbpb(muMiPt, muMiEta, 2, tnp_trgs);
    }
  }

  if (tnpTrigWeightMuPl < 0.0 || tnpTrigWeightMuMi < 0.0)
    return false;

  tnpWeight *= tnpTrigWeightMuPl;
  tnpWeight *= tnpTrigWeightMuMi;
  return true;
}

bool ComputePpTnPWeight(double muPlPt, double muPlEta,
                        ULong64_t muPlTrig,
                        double muMiPt, double muMiEta,
                        ULong64_t muMiTrig,
                        double &tnpWeight)
{
  const int kL2filter = 4;
  const ULong64_t l2Mask = (ULong64_t(1) << kL2filter);

  tnpWeight = 1.0;
  tnpWeight *= std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(muPlPt, muPlEta));
  tnpWeight *= std::get<0>(tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(muMiPt, muMiEta));

  if ((muPlTrig & l2Mask) != l2Mask || (muMiTrig & l2Mask) != l2Mask)
    return false;

  return true;
}

bool PassMuonQualityLikeOniaEff(bool hasIsGoodMuonBranch, Bool_t isGoodMuon, Bool_t isSoftCutBased,
                                Int_t selectionType, Int_t nTrkWMea, Int_t nPixWMea,
                                Float_t dxy, Float_t dz)
{
  const bool passMuonId = hasIsGoodMuonBranch
                              ? (isGoodMuon == 1)
                              : (isSoftCutBased || PassMuonSoftCutFallback(selectionType, nTrkWMea, nPixWMea, dxy, dz));

  return passMuonId &&
         (nTrkWMea > 5) &&
         (nPixWMea > 0) &&
         (std::fabs(dxy) < 0.3f) &&
         (std::fabs(dz) < 20.f);
}

bool BranchUsesShort(TTree *tree, const char *branchName)
{
  if (!tree)
    return false;
  TBranch *branch = tree->GetBranch(branchName);
  if (!branch)
    return false;
  TLeaf *leaf = branch->GetLeaf(branchName);
  if (!leaf && branch->GetListOfLeaves() && branch->GetListOfLeaves()->GetEntries() > 0)
    leaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->At(0));
  if (!leaf)
    return false;
  const std::string typeName = leaf->GetTypeName();
  return typeName == "Short_t" || typeName == "UShort_t";
}

TH1D *BuildWeightedEfficiencyHist(TH1D *num, TH1D *den, const char *name, const char *title)
{
  if (!num || !den)
    return nullptr;

  TH1D *eff = static_cast<TH1D *>(num->Clone(name));
  eff->SetTitle(title);
  // Use weighted ratio propagation. TnP makes the numerator non-binomial,
  // so the binomial "B" option is intentionally not used here.
  eff->Divide(num, den, 1.0, 1.0);
  return eff;
}

TH2D *BuildWeightedEfficiencyHist(TH2D *num, TH2D *den, const char *name, const char *title)
{
  if (!num || !den)
    return nullptr;

  TH2D *eff = static_cast<TH2D *>(num->Clone(name));
  eff->SetTitle(title);
  eff->Divide(num, den, 1.0, 1.0);
  return eff;
}

TH1D *ProjectCentHistFrom2D(TH2D *src, const char *name, const char *title)
{
  if (!src)
    return nullptr;

  TH1D *hist = src->ProjectionY(name, 1, src->GetNbinsX(), "e");
  hist->SetDirectory(nullptr);
  hist->SetTitle(title);
  return hist;
}

TH1D *CollapseToSingleBinHist(TH1D *src, const char *name, const char *title,
                              double xLow, double xHigh)
{
  if (!src)
    return nullptr;

  TH1D *hist = new TH1D(name, title, 1, xLow, xHigh);
  hist->SetDirectory(nullptr);

  double sum = 0.0;
  double err2 = 0.0;
  for (int ibin = 1; ibin <= src->GetNbinsX(); ++ibin)
  {
    sum += src->GetBinContent(ibin);
    err2 += std::pow(src->GetBinError(ibin), 2);
  }

  hist->SetBinContent(1, sum);
  hist->SetBinError(1, std::sqrt(err2));
  return hist;
}

TH1D *BuildCorrelatedEfficiencyHist(TH1D *num, TH1D *den, TH1D *cross,
                                    const char *name, const char *title)
{
  if (!num || !den || !cross)
    return nullptr;

  TH1D *eff = static_cast<TH1D *>(num->Clone(name));
  eff->SetTitle(title);
  eff->Reset("ICES");

  for (int i = 1; i <= num->GetNbinsX(); ++i)
  {
    const double x = num->GetBinContent(i);
    const double y = den->GetBinContent(i);
    if (y <= 0.0)
      continue;

    const double varX = std::pow(num->GetBinError(i), 2);
    const double varY = std::pow(den->GetBinError(i), 2);
    const double covXY = cross->GetBinContent(i);
    const double value = x / y;
    double variance = varX / (y * y) + (x * x * varY) / std::pow(y, 4) - (2.0 * x * covXY) / std::pow(y, 3);
    if (variance < 0.0)
      variance = 0.0;

    eff->SetBinContent(i, value);
    eff->SetBinError(i, std::sqrt(variance));
  }

  return eff;
}

void efficiency_1d_pp(long nEvt = 1000, bool isPr = true,
                   bool isNCollW = true, bool isTnPW = true, bool isGenW = true,
                   bool isPtW = true, bool useLegacy4PtWInputs = false,
                   bool isPbPb = true)
{
  cout << "Start onia_to_skim_data()\n";
  // Enable Sumw2 before any fill so weighted histograms keep statistical errors.
  TH1::SetDefaultSumw2(kTRUE);

  // ===== read input =====
  const char *treePath = "hionia/myTree";
  TChain *oniaTree = new TChain(treePath);

  const std::vector<std::string> inputFiles =
      isPbPb
          ? (isPr
                 ? std::vector<std::string>{"/data/Oniatree/miniAOD/OniaTreeMC_miniAOD_Jpsi_HydjetMB_5p02TeV_merged.root"}
                 : std::vector<std::string>{"/data/Oniatree/miniAOD/Oniatree_MC_BToJpsi_pThat-2_TuneCP5-EvtGen_HydjetDrumMB_5p02TeV_pythia8_miniAOD.root"})
          : (isPr
                 ? std::vector<std::string>{"/data/Oniatree/Jpsi/OniaTreeMC_Jpsi_Pythia8_nonUL_5p02TeV_merged.root"}
                 : std::vector<std::string>{"/data/Oniatree/Jpsi/BtoJpsiJpsi_miniAOD94X_20Jun_v1.root"});

  for (const auto &inputFile : inputFiles)
  {
    cout << "[INFO] adding input file: " << inputFile << "\n";
    oniaTree->Add(inputFile.c_str());
  }


  // ===== read pT reweight input =====
  TF1 *fPtW_mid = nullptr;
  TF1 *fPtW_fwd = nullptr;
  if (isPtW)
  {
    const bool useFitOutputMode = (!useLegacy4PtWInputs && !isPbPb);
    if (!useFitOutputMode)
    {
      if (!useLegacy4PtWInputs && isPbPb)
        cout << "[INFO] fit-output pT weights are only configured for pp here. Falling back to explicit PbPb ratioDataMC inputs.\n";

      const std::string ptWeightMidPath =
          isPbPb
              ? (isPr
                     ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y0_1p6_251118.root"
                     : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y0_1p6_251118.root")
              : (isPr
                     ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y0_1p6_251125.root"
                     : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y0_1p6_251125.root");
      const std::string ptWeightFwdPath =
          isPbPb
              ? (isPr
                     ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_Jpsi_DATA_ctauCut_y1p6_2p4_251118.root"
                     : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_AA_BtoJpsi_DATA_ctauCut_y1p6_2p4_251118.root")
              : (isPr
                     ? "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_Jpsi_DATA_ctauCut_y1p6_2p4_251125.root"
                     : "/data/users/pjgwak/work/raa_pb18/psi2S_RAA_PbPb2018/compareDataToMC/ratioDataMC_pp_BtoJpsi_DATA_ctauCut_y1p6_2p4_251125.root");
      const std::string ptWeightMidName = "dataMC_Ratio1";
      const std::string ptWeightFwdName = "dataMC_Ratio1";

      TFile *fPtWMidFile = TFile::Open(ptWeightMidPath.c_str(), "READ");
      if (!fPtWMidFile || fPtWMidFile->IsZombie())
      {
        cout << "[ERROR] failed to open mid-rapidity pT weight file: " << ptWeightMidPath << "\n";
        return;
      }

      TFile *fPtWFwdFile = TFile::Open(ptWeightFwdPath.c_str(), "READ");
      if (!fPtWFwdFile || fPtWFwdFile->IsZombie())
      {
        cout << "[ERROR] failed to open forward-rapidity pT weight file: " << ptWeightFwdPath << "\n";
        fPtWMidFile->Close();
        delete fPtWMidFile;
        return;
      }

      TF1 *fPtWMidIn = dynamic_cast<TF1 *>(fPtWMidFile->Get(ptWeightMidName.c_str()));
      TF1 *fPtWFwdIn = dynamic_cast<TF1 *>(fPtWFwdFile->Get(ptWeightFwdName.c_str()));
      if (!fPtWMidIn || !fPtWFwdIn)
      {
        cout << "[ERROR] failed to load pT weight functions: "
             << ptWeightMidPath << "::" << ptWeightMidName << ", "
             << ptWeightFwdPath << "::" << ptWeightFwdName << "\n";
        fPtWMidFile->Close();
        fPtWFwdFile->Close();
        delete fPtWMidFile;
        delete fPtWFwdFile;
        return;
      }

      fPtW_mid = static_cast<TF1 *>(fPtWMidIn->Clone(Form("%s_clone", ptWeightMidName.c_str())));
      fPtW_fwd = static_cast<TF1 *>(fPtWFwdIn->Clone(Form("%s_clone", ptWeightFwdName.c_str())));
      fPtWMidFile->Close();
      fPtWFwdFile->Close();
      delete fPtWMidFile;
      delete fPtWFwdFile;

      cout << "[INFO] pT reweight inputs (explicit ratioDataMC mode): "
           << ptWeightMidPath << "::" << ptWeightMidName << ", "
           << ptWeightFwdPath << "::" << ptWeightFwdName << "\n";
    }
    else
    {
      const std::string ptWeightFilePath =
          "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/pT_reweight/fit_ratio_outputs/roots/fit_ratio_exp2.root";
      const std::string ptWeightMidName = isPr ? "exp2_pr_mid" : "exp2_np_mid";
      const std::string ptWeightFwdName = isPr ? "exp2_pr_fwd" : "exp2_np_fwd";

      TFile *fPtWFile = TFile::Open(ptWeightFilePath.c_str(), "READ");
      if (!fPtWFile || fPtWFile->IsZombie())
      {
        cout << "[ERROR] failed to open pT weight file: " << ptWeightFilePath << "\n";
        return;
      }

      TF1 *fPtWMidIn = dynamic_cast<TF1 *>(fPtWFile->Get(ptWeightMidName.c_str()));
      TF1 *fPtWFwdIn = dynamic_cast<TF1 *>(fPtWFile->Get(ptWeightFwdName.c_str()));
      if (!fPtWMidIn || !fPtWFwdIn)
      {
        cout << "[ERROR] failed to load pT weight functions: "
             << ptWeightFilePath << "::" << ptWeightMidName << ", "
             << ptWeightFilePath << "::" << ptWeightFwdName << "\n";
        fPtWFile->Close();
        delete fPtWFile;
        return;
      }

      fPtW_mid = static_cast<TF1 *>(fPtWMidIn->Clone(Form("%s_clone", ptWeightMidName.c_str())));
      fPtW_fwd = static_cast<TF1 *>(fPtWFwdIn->Clone(Form("%s_clone", ptWeightFwdName.c_str())));
      fPtWFile->Close();
      delete fPtWFile;

      cout << "[INFO] pT reweight inputs (fit-output mode): "
           << ptWeightFilePath << "::" << ptWeightMidName << ", "
           << ptWeightFilePath << "::" << ptWeightFwdName << "\n";
    }
  }

  // user defined kinematic cuts
  int cLow = 0, cHigh = 180;
  double massLow = 2.5, massHigh = 3.5;

  // ===== labeling =====
  std::string mcLabel = isPr ? "PR" : "NP";
  const std::string collLabel = isPbPb ? "PbPb2018" : "pp5p02TeV";
  const std::string weightLabel = Form("_ncollW%d_genW%d_ptW%d_tnpW%d",
                                       isNCollW, isGenW, isPtW, isTnPW);
  if (!isPbPb && isNCollW)
    cout << "[INFO] isNCollW is ignored for pp inputs.\n";

  // ===== set Oniatree branch address (input) =====
  const long int maxBranchSize = 1000;
  const bool usesSplitMomentumBranches = (oniaTree->GetBranch("Reco_QQ_4mom_pt") != nullptr);
  const bool usesClonesMomentumBranches = (oniaTree->GetBranch("Reco_QQ_4mom") != nullptr);
  const bool useShortSizeBranches = BranchUsesShort(oniaTree, "Gen_QQ_size");
  const bool useShortRecoIndexBranches = BranchUsesShort(oniaTree, "Reco_QQ_mupl_idx");
  const bool useShortGenIndexBranches = BranchUsesShort(oniaTree, "Gen_QQ_mupl_idx");
  const bool useShortRecoMuWhichGenBranch = BranchUsesShort(oniaTree, "Reco_mu_whichGen");
  const bool useShortRecoMuTypeBranch = BranchUsesShort(oniaTree, "Reco_mu_type");
  const bool hasGenWeightBranch = (oniaTree->GetBranch("Gen_weight") != nullptr);
  const bool hasRecoMuTypeBranch = (oniaTree->GetBranch("Reco_mu_type") != nullptr);
  const bool hasRecoMuIsGlobalBranch = (oniaTree->GetBranch("Reco_mu_isGlobal") != nullptr);

  std::vector<float> *Reco_QQ_4mom_m = nullptr;
  std::vector<float> *Reco_QQ_4mom_pt = nullptr;
  std::vector<float> *Reco_QQ_4mom_phi = nullptr;
  std::vector<float> *Reco_QQ_4mom_eta = nullptr;
  
  std::vector<float> *Reco_mu_4mom_pt = nullptr;
  std::vector<float> *Reco_mu_4mom_phi = nullptr;
  std::vector<float> *Reco_mu_4mom_eta = nullptr;
  std::vector<float> *Reco_mu_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_m = nullptr;
  std::vector<float> *Gen_QQ_4mom_pt = nullptr;
  std::vector<float> *Gen_QQ_4mom_phi = nullptr;
  std::vector<float> *Gen_QQ_4mom_eta = nullptr;
  std::vector<float> *Gen_mu_4mom_pt = nullptr;
  std::vector<float> *Gen_mu_4mom_phi = nullptr;
  std::vector<float> *Gen_mu_4mom_eta = nullptr;
  std::vector<float> *Gen_mu_4mom_m = nullptr;
  Float_t Gen_weight = 1.0f;
  TClonesArray *Reco_QQ_4mom = nullptr;
  TClonesArray *Reco_mu_4mom = nullptr;
  TClonesArray *Gen_QQ_4mom = nullptr;
  TClonesArray *Gen_mu_4mom = nullptr;

  if (hasGenWeightBranch)
    oniaTree->SetBranchAddress("Gen_weight", &Gen_weight);
  else
    cout << "[WARN] Gen_weight branch not found. Falling back to unit GEN weight." << std::endl;

  if (usesSplitMomentumBranches)
  {
    oniaTree->SetBranchAddress("Reco_QQ_4mom_m", &Reco_QQ_4mom_m);
    oniaTree->SetBranchAddress("Reco_QQ_4mom_pt", &Reco_QQ_4mom_pt);
    oniaTree->SetBranchAddress("Reco_QQ_4mom_phi", &Reco_QQ_4mom_phi);
    oniaTree->SetBranchAddress("Reco_QQ_4mom_eta", &Reco_QQ_4mom_eta);

    oniaTree->SetBranchAddress("Reco_mu_4mom_pt", &Reco_mu_4mom_pt);
    oniaTree->SetBranchAddress("Reco_mu_4mom_phi", &Reco_mu_4mom_phi);
    oniaTree->SetBranchAddress("Reco_mu_4mom_eta", &Reco_mu_4mom_eta);
    oniaTree->SetBranchAddress("Reco_mu_4mom_m", &Reco_mu_4mom_m);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_m", &Gen_QQ_4mom_m);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_pt", &Gen_QQ_4mom_pt);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_phi", &Gen_QQ_4mom_phi);
    oniaTree->SetBranchAddress("Gen_QQ_4mom_eta", &Gen_QQ_4mom_eta);
    oniaTree->SetBranchAddress("Gen_mu_4mom_pt", &Gen_mu_4mom_pt);
    oniaTree->SetBranchAddress("Gen_mu_4mom_phi", &Gen_mu_4mom_phi);
    oniaTree->SetBranchAddress("Gen_mu_4mom_eta", &Gen_mu_4mom_eta);
    oniaTree->SetBranchAddress("Gen_mu_4mom_m", &Gen_mu_4mom_m);
  }
  else if (usesClonesMomentumBranches)
  {
    oniaTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
    oniaTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
    oniaTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom);
    oniaTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom);
  }
  else
  {
    cout << "[ERROR] unsupported onia tree schema: momentum branches missing\n";
    return;
  }


  // event-level scalars
  UInt_t runNb = 0;
  UInt_t eventNb = 0;
  UInt_t LS = 0;
  float zVtx = 0.f;
  Int_t Centrality = 0;
  ULong64_t HLTriggers;
  // // collection sizes
  Short_t Gen_QQ_size_short = 0;
  Short_t Reco_QQ_size_short = 0;
  Int_t Gen_QQ_size_int = 0;
  Int_t Reco_QQ_size_int = 0;

  // // per-dimuon (size = Reco_QQ_size)
  ULong64_t Reco_QQ_trig[maxBranchSize];
  Float_t Reco_QQ_VtxProb[maxBranchSize];
  Short_t Reco_QQ_mupl_idx_short[maxBranchSize];
  Short_t Reco_QQ_mumi_idx_short[maxBranchSize];
  Short_t Reco_QQ_sign_short[maxBranchSize];
  Int_t Reco_QQ_mupl_idx[maxBranchSize];
  Int_t Reco_QQ_mumi_idx[maxBranchSize];
  Int_t Reco_QQ_sign[maxBranchSize];
  Float_t Reco_QQ_ctau3D[maxBranchSize];
  Float_t Reco_QQ_ctauErr3D[maxBranchSize];
  Float_t Reco_QQ_ctau[maxBranchSize];
  Float_t Reco_QQ_ctauErr[maxBranchSize];
  Short_t Gen_QQ_mupl_idx_short[maxBranchSize];
  Short_t Gen_QQ_mumi_idx_short[maxBranchSize];
  Short_t Gen_mu_charge_short[maxBranchSize];
  Int_t Gen_QQ_mupl_idx[maxBranchSize];
  Int_t Gen_QQ_mumi_idx[maxBranchSize];
  Int_t Gen_mu_charge[maxBranchSize];

  // // per-muon (size = Reco_mu_size)
  ULong64_t Reco_mu_trig[maxBranchSize];
  Bool_t Reco_mu_isGoodMuon[maxBranchSize];
  Bool_t Reco_mu_highPurity[maxBranchSize];
  Int_t Reco_mu_nTrkHits[maxBranchSize];
  // Float_t Reco_mu_normChi2_global[maxBranchSize];
  Int_t Reco_mu_nMuValHits[maxBranchSize];
  // Int_t Reco_mu_StationsMatched[maxBranchSize];
  Float_t Reco_mu_dxy[maxBranchSize];
  // Float_t Reco_mu_dxyErr[maxBranchSize];
  Float_t Reco_mu_dz[maxBranchSize];
  Float_t Reco_mu_dzErr[maxBranchSize];
  Int_t Reco_mu_nTrkWMea[maxBranchSize];
  Int_t Reco_mu_nPixWMea[maxBranchSize];
  Int_t Reco_mu_nPixValHits[maxBranchSize];
  Int_t Reco_mu_SelectionType[maxBranchSize];
  Short_t Reco_mu_type_short[maxBranchSize];
  Int_t Reco_mu_type[maxBranchSize];
  Short_t Reco_mu_whichGen_short[maxBranchSize];
  Int_t Reco_mu_whichGen[maxBranchSize];
  Bool_t Reco_mu_isSoftCutBased[maxBranchSize];
  Bool_t Reco_mu_isGlobal[maxBranchSize];
  Float_t Reco_QQ_dca[maxBranchSize];

  // // ----- SetBranchAddress -----
  if (oniaTree->GetBranch("runNb"))
    oniaTree->SetBranchAddress("runNb", &runNb);
  if (oniaTree->GetBranch("eventNb"))
    oniaTree->SetBranchAddress("eventNb", &eventNb);
  if (oniaTree->GetBranch("LS"))
    oniaTree->SetBranchAddress("LS", &LS);
  if (oniaTree->GetBranch("zVtx"))
    oniaTree->SetBranchAddress("zVtx", &zVtx);
  if (oniaTree->GetBranch("Centrality"))
    oniaTree->SetBranchAddress("Centrality", &Centrality);
  oniaTree->SetBranchAddress("HLTriggers", &HLTriggers);

  // // sizes
  if (useShortSizeBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_short);
    oniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size_short);
  }
  else
  {
    oniaTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size_int);
    oniaTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size_int);
  }

  // // dimuon
  oniaTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig);
  oniaTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb);
  if (useShortRecoIndexBranches)
  {
    oniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx_short);
    oniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx_short);
    oniaTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign_short);
  }
  else
  {
    oniaTree->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
    oniaTree->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
    oniaTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  }
  if (useShortGenIndexBranches)
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx_short);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx_short);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge_short);
  }
  else
  {
    oniaTree->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx);
    oniaTree->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx);
    oniaTree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge);
  }
  oniaTree->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D);
  oniaTree->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau);
  oniaTree->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr);

  // single-muon
  oniaTree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig);
  if (oniaTree->GetBranch("Reco_mu_isGoodMuon"))
    oniaTree->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon);
  oniaTree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity);
  if (!usesSplitMomentumBranches && oniaTree->GetBranch("Reco_mu_type"))
  {
    if (useShortRecoMuTypeBranch)
      oniaTree->SetBranchAddress("Reco_mu_type", Reco_mu_type_short);
    else
      oniaTree->SetBranchAddress("Reco_mu_type", Reco_mu_type);
  }
  oniaTree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits);
  // oniaTree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global);
  oniaTree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits);
  // oniaTree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched);
  if (oniaTree->GetBranch("Reco_mu_dxy"))
    oniaTree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy);
  // oniaTree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr);
  if (oniaTree->GetBranch("Reco_mu_dz"))
    oniaTree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz);
  // oniaTree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr);
  oniaTree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea);
  oniaTree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits);
  if (oniaTree->GetBranch("Reco_mu_SelectionType"))
    oniaTree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType);
  if (oniaTree->GetBranch("Reco_mu_isSoftCutBased"))
    oniaTree->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased);
  if (oniaTree->GetBranch("Reco_mu_isGlobal"))
    oniaTree->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal);
  oniaTree->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca);

  if (oniaTree->GetBranch("Reco_mu_whichGen"))
  {
    if (useShortRecoMuWhichGenBranch)
      oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen_short);
    else
      oniaTree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen);
  }

  // ----- declare local variables ----
  const static long long int nMaxDimu = 1000;
  
  // event-level scalar
  Int_t event = 0;
  Int_t runN = 0;
  Int_t lumi = 0;
  Float_t vz = 0;

  // dimuon size
  Int_t nDimu = 0;

  // dimuon candidate-level variables (size = nDimu)
  Float_t mass[nMaxDimu];
  Float_t y[nMaxDimu];
  Float_t pt[nMaxDimu], pt1[nMaxDimu], pt2[nMaxDimu];
  Float_t eta[nMaxDimu], eta1[nMaxDimu], eta2[nMaxDimu];
  Float_t phi[nMaxDimu], phi1[nMaxDimu], phi2[nMaxDimu];
  Int_t recoQQsign[nMaxDimu], cBin[nMaxDimu];
  Float_t ctau3D[nMaxDimu], ctau3DErr[nMaxDimu], ctau3DRes[nMaxDimu];
  Float_t ctau[nMaxDimu], ctauErr[nMaxDimu], ctauRes[nMaxDimu];
  Float_t recoQQdca[nMaxDimu];

  // weights
  Double_t weight = 1.0; // event-level weight
  Double_t TnPweight[nMaxDimu]; // per-candidate level weight

  const char *outDir = "/data/users/pjgwak/work/daily_code_tracker/2026/03/26_pp_run2_correction_study/acc_eff_skim/skim_roots";
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] failed to create output directory: " << outDir << "\n";
    return;
  }

  std::string outFilePath = Form(
      "%s/eff_%s_isMC1%s%s_Dimuon_MiniAOD.root",
      outDir, collLabel.c_str(), mcLabel.empty() ? "" : ("_" + mcLabel).c_str(), weightLabel.c_str());
  TFile *fFlowSkim = TFile::Open(outFilePath.c_str(), "RECREATE");
  if (!fFlowSkim || fFlowSkim->IsZombie() || !fFlowSkim->IsWritable())
  {
    cout << "[ERROR] cannot open output file for writing: " << outFilePath << "\n";
    return;
  }

  const double ptBinsMid[] = {6.5, 9, 12, 15, 20, 25, 30, 40};
  const double ptBinsFwd[] = {3.5, 5, 6.5, 9, 12, 15, 20, 25, 30, 40};
  const double ptBinsRap[] = {6.5, 9, 12, 15, 20, 25, 30, 40};
  const double ptBinsRap1624[] = {3.5, 5, 6.5, 9, 12, 15, 20, 25, 30, 40};
  const double ptBinsFwdLow[] = {3.5, 5, 6.5};
  const double centBins[] = {0, 20, 40, 60, 80, 100, 180};
  const double fineCentBins[] = {0, 10, 20, 30, 40, 50, 60, 80, 100, 180};
  const double yBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  const int nPtBinsMid = sizeof(ptBinsMid) / sizeof(double) - 1;
  const int nPtBinsFwd = sizeof(ptBinsFwd) / sizeof(double) - 1;
  const int nPtBinsRap = sizeof(ptBinsRap) / sizeof(double) - 1;
  const int nPtBinsRap1624 = sizeof(ptBinsRap1624) / sizeof(double) - 1;
  const int nPtBinsFwdLow = sizeof(ptBinsFwdLow) / sizeof(double) - 1;
  const int nCentBins = sizeof(centBins) / sizeof(double) - 1;
  const int nFineCentBins = sizeof(fineCentBins) / sizeof(double) - 1;
  const int nYBins = sizeof(yBins) / sizeof(double) - 1;

  TH1D *hist_eff_den_mid = new TH1D("hist_eff_den_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_den_fwd = new TH1D("hist_eff_den_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_num_mid = new TH1D("hist_eff_num_mid", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_num_fwd = new TH1D("hist_eff_num_fwd", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_den_cent_mid = new TH1D("hist_eff_den_cent_mid", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_fwd = new TH1D("hist_eff_den_cent_fwd", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_mid = new TH1D("hist_eff_num_cent_mid", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_fwd = new TH1D("hist_eff_num_cent_fwd", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_pt_cent0010 = new TH1D("hist_eff_den_pt_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_pt_cent1030 = new TH1D("hist_eff_den_pt_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_pt_cent3090 = new TH1D("hist_eff_den_pt_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_pt_mid_cent0010 = new TH1D("hist_eff_den_pt_mid_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_den_pt_mid_cent1030 = new TH1D("hist_eff_den_pt_mid_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_den_pt_mid_cent3090 = new TH1D("hist_eff_den_pt_mid_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_den_pt_fwd_cent0010 = new TH1D("hist_eff_den_pt_fwd_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_den_pt_fwd_cent1030 = new TH1D("hist_eff_den_pt_fwd_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_den_pt_fwd_cent3090 = new TH1D("hist_eff_den_pt_fwd_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_num_pt_cent0010 = new TH1D("hist_eff_num_pt_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_pt_cent1030 = new TH1D("hist_eff_num_pt_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_pt_cent3090 = new TH1D("hist_eff_num_pt_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_pt_mid_cent0010 = new TH1D("hist_eff_num_pt_mid_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_num_pt_mid_cent1030 = new TH1D("hist_eff_num_pt_mid_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_num_pt_mid_cent3090 = new TH1D("hist_eff_num_pt_mid_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsMid, ptBinsMid);
  TH1D *hist_eff_num_pt_fwd_cent0010 = new TH1D("hist_eff_num_pt_fwd_cent0010", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_num_pt_fwd_cent1030 = new TH1D("hist_eff_num_pt_fwd_cent1030", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_num_pt_fwd_cent3090 = new TH1D("hist_eff_num_pt_fwd_cent3090", ";p_{T} (GeV/c);Counts", nPtBinsFwd, ptBinsFwd);
  TH1D *hist_eff_den_rap0024 = new TH1D("hist_eff_den_rap0024", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_rap0006 = new TH1D("hist_eff_den_rap0006", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_rap0612 = new TH1D("hist_eff_den_rap0612", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_rap1216 = new TH1D("hist_eff_den_rap1216", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_den_rap1624 = new TH1D("hist_eff_den_rap1624", ";p_{T} (GeV/c);Counts", nPtBinsRap1624, ptBinsRap1624);
  TH1D *hist_eff_den_cent_rap0024 = new TH1D("hist_eff_den_cent_rap0024", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_rap0006 = new TH1D("hist_eff_den_cent_rap0006", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_rap0612 = new TH1D("hist_eff_den_cent_rap0612", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_rap1216 = new TH1D("hist_eff_den_cent_rap1216", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_rap1624 = new TH1D("hist_eff_den_cent_rap1624", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_cent_rap1624_pt3565 = new TH1D("hist_eff_den_cent_rap1624_pt3565", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_rap0024 = new TH1D("hist_eff_num_rap0024", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_rap0006 = new TH1D("hist_eff_num_rap0006", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_rap0612 = new TH1D("hist_eff_num_rap0612", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_rap1216 = new TH1D("hist_eff_num_rap1216", ";p_{T} (GeV/c);Counts", nPtBinsRap, ptBinsRap);
  TH1D *hist_eff_num_rap1624 = new TH1D("hist_eff_num_rap1624", ";p_{T} (GeV/c);Counts", nPtBinsRap1624, ptBinsRap1624);
  TH1D *hist_eff_num_cent_rap0024 = new TH1D("hist_eff_num_cent_rap0024", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_rap0006 = new TH1D("hist_eff_num_cent_rap0006", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_rap0612 = new TH1D("hist_eff_num_cent_rap0612", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_rap1216 = new TH1D("hist_eff_num_cent_rap1216", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_rap1624 = new TH1D("hist_eff_num_cent_rap1624", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_num_cent_rap1624_pt3565 = new TH1D("hist_eff_num_cent_rap1624_pt3565", ";Centrality;Counts", nCentBins, centBins);
  TH1D *hist_eff_den_y = new TH1D("hist_eff_den_y", ";|y|;Counts", nYBins, yBins);
  TH1D *hist_eff_num_y = new TH1D("hist_eff_num_y", ";|y|;Counts", nYBins, yBins);
  TH2D *hist_eff2d_den_mid_ptw = new TH2D("hist_eff2d_den_mid_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nCentBins, centBins);
  TH2D *hist_eff2d_den_mid_noptw = new TH2D("hist_eff2d_den_mid_noptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nCentBins, centBins);
  TH2D *hist_eff2d_num_mid_ptw = new TH2D("hist_eff2d_num_mid_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nCentBins, centBins);
  TH2D *hist_eff2d_num_mid_noptw = new TH2D("hist_eff2d_num_mid_noptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nCentBins, centBins);
  TH2D *hist_eff2d_den_fwd_ptw = new TH2D("hist_eff2d_den_fwd_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nCentBins, centBins);
  TH2D *hist_eff2d_den_fwd_noptw = new TH2D("hist_eff2d_den_fwd_noptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nCentBins, centBins);
  TH2D *hist_eff2d_num_fwd_ptw = new TH2D("hist_eff2d_num_fwd_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nCentBins, centBins);
  TH2D *hist_eff2d_num_fwd_noptw = new TH2D("hist_eff2d_num_fwd_noptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nCentBins, centBins);
  TH2D *hist_eff2d_den_mid_fine_ptw = new TH2D("hist_eff2d_den_mid_fine_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nFineCentBins, fineCentBins);
  TH2D *hist_eff2d_num_mid_fine_ptw = new TH2D("hist_eff2d_num_mid_fine_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsMid, ptBinsMid, nFineCentBins, fineCentBins);
  TH2D *hist_eff2d_den_fwd_fine_ptw = new TH2D("hist_eff2d_den_fwd_fine_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nFineCentBins, fineCentBins);
  TH2D *hist_eff2d_num_fwd_fine_ptw = new TH2D("hist_eff2d_num_fwd_fine_ptw", ";p_{T} (GeV/c);Centrality", nPtBinsFwd, ptBinsFwd, nFineCentBins, fineCentBins);

  auto makeCrossHist = [](TH1D *src, const char *name) {
    TH1D *hist = static_cast<TH1D *>(src->Clone(name));
    hist->Reset("ICES");
    return hist;
  };

  TH1D *hist_eff_cross_mid = makeCrossHist(hist_eff_num_mid, "hist_eff_cross_mid");
  TH1D *hist_eff_cross_fwd = makeCrossHist(hist_eff_num_fwd, "hist_eff_cross_fwd");
  TH1D *hist_eff_cross_cent_mid = makeCrossHist(hist_eff_num_cent_mid, "hist_eff_cross_cent_mid");
  TH1D *hist_eff_cross_cent_fwd = makeCrossHist(hist_eff_num_cent_fwd, "hist_eff_cross_cent_fwd");
  TH1D *hist_eff_cross_pt_cent0010 = makeCrossHist(hist_eff_num_pt_cent0010, "hist_eff_cross_pt_cent0010");
  TH1D *hist_eff_cross_pt_cent1030 = makeCrossHist(hist_eff_num_pt_cent1030, "hist_eff_cross_pt_cent1030");
  TH1D *hist_eff_cross_pt_cent3090 = makeCrossHist(hist_eff_num_pt_cent3090, "hist_eff_cross_pt_cent3090");
  TH1D *hist_eff_cross_pt_mid_cent0010 = makeCrossHist(hist_eff_num_pt_mid_cent0010, "hist_eff_cross_pt_mid_cent0010");
  TH1D *hist_eff_cross_pt_mid_cent1030 = makeCrossHist(hist_eff_num_pt_mid_cent1030, "hist_eff_cross_pt_mid_cent1030");
  TH1D *hist_eff_cross_pt_mid_cent3090 = makeCrossHist(hist_eff_num_pt_mid_cent3090, "hist_eff_cross_pt_mid_cent3090");
  TH1D *hist_eff_cross_pt_fwd_cent0010 = makeCrossHist(hist_eff_num_pt_fwd_cent0010, "hist_eff_cross_pt_fwd_cent0010");
  TH1D *hist_eff_cross_pt_fwd_cent1030 = makeCrossHist(hist_eff_num_pt_fwd_cent1030, "hist_eff_cross_pt_fwd_cent1030");
  TH1D *hist_eff_cross_pt_fwd_cent3090 = makeCrossHist(hist_eff_num_pt_fwd_cent3090, "hist_eff_cross_pt_fwd_cent3090");
  TH1D *hist_eff_cross_rap0024 = makeCrossHist(hist_eff_num_rap0024, "hist_eff_cross_rap0024");
  TH1D *hist_eff_cross_rap0006 = makeCrossHist(hist_eff_num_rap0006, "hist_eff_cross_rap0006");
  TH1D *hist_eff_cross_rap0612 = makeCrossHist(hist_eff_num_rap0612, "hist_eff_cross_rap0612");
  TH1D *hist_eff_cross_rap1216 = makeCrossHist(hist_eff_num_rap1216, "hist_eff_cross_rap1216");
  TH1D *hist_eff_cross_rap1624 = makeCrossHist(hist_eff_num_rap1624, "hist_eff_cross_rap1624");
  TH1D *hist_eff_cross_cent_rap0024 = makeCrossHist(hist_eff_num_cent_rap0024, "hist_eff_cross_cent_rap0024");
  TH1D *hist_eff_cross_cent_rap0006 = makeCrossHist(hist_eff_num_cent_rap0006, "hist_eff_cross_cent_rap0006");
  TH1D *hist_eff_cross_cent_rap0612 = makeCrossHist(hist_eff_num_cent_rap0612, "hist_eff_cross_cent_rap0612");
  TH1D *hist_eff_cross_cent_rap1216 = makeCrossHist(hist_eff_num_cent_rap1216, "hist_eff_cross_cent_rap1216");
  TH1D *hist_eff_cross_cent_rap1624 = makeCrossHist(hist_eff_num_cent_rap1624, "hist_eff_cross_cent_rap1624");
  TH1D *hist_eff_cross_cent_rap1624_pt3565 = makeCrossHist(hist_eff_num_cent_rap1624_pt3565, "hist_eff_cross_cent_rap1624_pt3565");
  TH1D *hist_eff_cross_y = makeCrossHist(hist_eff_num_y, "hist_eff_cross_y");
  

  // fFlowSkim->SetCompressionSettings(207); // LZMA: 207, ZLIB: 1xx
  TTree *flowTree = new TTree("myTree", "");
  flowTree->SetMaxTreeSize(MAXTREESIZE);
  // flowTree->SetAutoSave(50 * 1024 * 1024); // save every 50 MB

  // ----- SetBranchAddress -----
  // event values
  flowTree->Branch("eventNb", &event, "eventNb/I");
  flowTree->Branch("runNb", &runN, "runNb/I");
  flowTree->Branch("LS", &lumi, "LS/I");
  flowTree->Branch("zVtx", &vz, "zVtx/F");

  // dimuon size
  flowTree->Branch("nDimu", &nDimu, "nDimu/I");

  // per-candidate (size = [nDimu])
  flowTree->Branch("mass", mass, "mass[nDimu]/F");
  flowTree->Branch("y", y, "y[nDimu]/F");
  flowTree->Branch("cBin", cBin, "cBin[nDimu]/I");
  flowTree->Branch("pt", pt, "pt[nDimu]/F");
  flowTree->Branch("pt1", pt1, "pt1[nDimu]/F");
  flowTree->Branch("pt2", pt2, "pt2[nDimu]/F");
  flowTree->Branch("eta", eta, "eta[nDimu]/F");
  flowTree->Branch("eta1", eta1, "eta1[nDimu]/F");
  flowTree->Branch("eta2", eta2, "eta2[nDimu]/F");
  flowTree->Branch("phi", phi, "phi[nDimu]/F");
  flowTree->Branch("phi1", phi1, "phi1[nDimu]/F");
  flowTree->Branch("phi2", phi2, "phi2[nDimu]/F");
  flowTree->Branch("recoQQsign", recoQQsign, "recoQQsign[nDimu]/I");
  flowTree->Branch("ctau3D", ctau3D, "ctau3D[nDimu]/F");
  flowTree->Branch("ctau3DErr", ctau3DErr, "ctau3DErr[nDimu]/F");
  flowTree->Branch("ctau3DRes", ctau3DRes, "ctau3DRes[nDimu]/F");
  flowTree->Branch("ctau", ctau, "ctau[nDimu]/F");
  flowTree->Branch("ctauErr", ctauErr, "ctauErr[nDimu]/F");
  flowTree->Branch("ctauRes", ctauRes, "ctauRes[nDimu]/F");
  flowTree->Branch("recoQQdca", recoQQdca, "recoQQdca[nDimu]/F");

  // weights
  flowTree->Branch("weight", &weight, "weight/D");
  flowTree->Branch("TnPweight", TnPweight, "TnPweight[nDimu]/D");

  // ===== event loop =====
  // counters
  long count_gen_dimuon = 0;
  long count_reco_dimuon = 0;
  long recoCandTotal = 0;
  long recoPassMass = 0;
  long recoPassAcceptance = 0;
  long recoPassBin = 0;
  long recoPassSign = 0;
  long recoPassMuonType = 0;
  long recoPassMuonQuality = 0;
  long recoPassVtx = 0;
  long recoPassTrig = 0;
  long recoPassMatched = 0;

  // ---- preaparation before loop -----
  const Long64_t nTot = oniaTree->GetEntries();
  if (nEvt == -1)
    nEvt = nTot;
  nEvt = std::min<Long64_t>(nEvt, nTot);

  // ----- verbose option -----
  const bool VERBOSE_EVENT = false;  // event-level
  const bool VERBOSE_CANDID = false; // dimuon-level

  // ----- start event loop -----
  for (long iev = 0; iev < nEvt; ++iev)
  {
    // print progress
    if ((iev % 100000 == 0) && (iev != 0))
    {
      cout << ">>>>> EVENT " << iev << " / " << nTot
                << " (" << int(100. * iev / std::max<Long64_t>(1, nTot)) << "%)\n";
      cout << "GEN accepted dimuon: " << count_gen_dimuon << "\n";
      cout << "Reco selected dimuon: " << count_reco_dimuon << "\n";
    }
    if (VERBOSE_EVENT) {
      cout << "Total Event: " << iev << "\n";
      cout << "Total GEN accepted dimuon: " << count_gen_dimuon << "\n";
      cout << "Total Reco selected dimuon: " << count_reco_dimuon << "\n";
    }

    // read event in Oniatree
    oniaTree->GetEntry(iev);

    // cout << "is Reco_QQ_4mom_pt->size() == Reco_QQ_size: " << (n == qq_size) << endl;
    
    // event-level values
    runN = runNb;
    event = eventNb;
    lumi = LS;
    vz = zVtx;

    double Ncoll_weight = 1.0;
    if (isNCollW && isPbPb)
    {
      Ncoll_weight = findNcoll(Centrality);
    }
    double treeWeight = 1.0;
    if (isGenW)
    {
      treeWeight = hasGenWeightBranch ? static_cast<double>(Gen_weight) : 1.0;
    }
    const double eventWeightBase = Ncoll_weight * treeWeight;

    if (isPbPb && (Centrality < cLow || Centrality >= cHigh))
      continue;

    const Int_t Gen_QQ_size = useShortSizeBranches
                                  ? static_cast<Int_t>(Gen_QQ_size_short)
                                  : Gen_QQ_size_int;
    const Int_t Reco_QQ_size = useShortSizeBranches
                                   ? static_cast<Int_t>(Reco_QQ_size_short)
                                   : Reco_QQ_size_int;
    auto GetRecoQQMuplIdx = [&](Int_t idx) -> Int_t {
      return useShortRecoIndexBranches
                 ? static_cast<Int_t>(Reco_QQ_mupl_idx_short[idx])
                 : Reco_QQ_mupl_idx[idx];
    };
    auto GetRecoQQMumiIdx = [&](Int_t idx) -> Int_t {
      return useShortRecoIndexBranches
                 ? static_cast<Int_t>(Reco_QQ_mumi_idx_short[idx])
                 : Reco_QQ_mumi_idx[idx];
    };
    auto GetRecoQQSign = [&](Int_t idx) -> Int_t {
      return useShortRecoIndexBranches
                 ? static_cast<Int_t>(Reco_QQ_sign_short[idx])
                 : Reco_QQ_sign[idx];
    };
    auto GetGenQQMuplIdx = [&](Int_t idx) -> Int_t {
      return useShortGenIndexBranches
                 ? static_cast<Int_t>(Gen_QQ_mupl_idx_short[idx])
                 : Gen_QQ_mupl_idx[idx];
    };
    auto GetGenQQMumiIdx = [&](Int_t idx) -> Int_t {
      return useShortGenIndexBranches
                 ? static_cast<Int_t>(Gen_QQ_mumi_idx_short[idx])
                 : Gen_QQ_mumi_idx[idx];
    };
    auto GetGenMuCharge = [&](Int_t idx) -> Int_t {
      return useShortGenIndexBranches
                 ? static_cast<Int_t>(Gen_mu_charge_short[idx])
                 : Gen_mu_charge[idx];
    };
    auto GetRecoMuWhichGen = [&](Int_t idx) -> Int_t {
      return useShortRecoMuWhichGenBranch
                 ? static_cast<Int_t>(Reco_mu_whichGen_short[idx])
                 : Reco_mu_whichGen[idx];
    };
    auto GetRecoMuType = [&](Int_t idx) -> Int_t {
      return useShortRecoMuTypeBranch
                 ? static_cast<Int_t>(Reco_mu_type_short[idx])
                 : Reco_mu_type[idx];
    };

    // check dimuon number
    if (Gen_QQ_size < 0) continue;
    if (Reco_QQ_size<0) continue;
    if (Gen_QQ_size != 1) continue;

    const int kTrigSel = isPbPb ? 12 : 3;
    const bool HLTPass =
        ((HLTriggers & ((ULong64_t)std::pow(2, kTrigSel))) == ((ULong64_t)std::pow(2, kTrigSel)));

    // if (TMath::Abs(vz) > 15) continue; // not included in the Jun Chen's slide

    // ===== start GEN dimuon loop =====
    // reset number of GEN-dimuon
    int nGenDimu = 0;
    bool hasGenReference = false;
    float refGenPt = 0.f;
    float refGenAbsY = 0.f;
    double refGenPtWeight = 1.0;
    TLorentzVector refGenMuPl;
    TLorentzVector refGenMuMi;

    for (Int_t irqq = 0; irqq < Gen_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID)
        cout << "  irqq: " << irqq << "\n";

      // check maxNDimuon
      if (nGenDimu >= nMaxDimu)
      {
        std::cerr << "[ERROR] nGenDimu reached nMaxDimu at event " << iev << "\n";
        break;
      }

      // check TLorentzVector pointers
      const int iMuPl = GetGenQQMuplIdx(irqq);
      const int iMuMi = GetGenQQMumiIdx(irqq);
      if (iMuPl < 0 || iMuMi < 0) continue;

      // user-defined kinematic cuts
      float gen_mass = 0.f;
      float gen_eta = 0.f;
      float gen_eta1 = 0.f;
      float gen_eta2 = 0.f;
      float gen_pt = 0.f;
      float gen_pt1 = 0.f;
      float gen_pt2 = 0.f;
      float gen_phi = 0.f;
      float gen_phi1 = 0.f;
      float gen_phi2 = 0.f;
      float gen_m1 = 0.105658f;
      float gen_m2 = 0.105658f;
      if (usesSplitMomentumBranches)
      {
        if (!Gen_QQ_4mom_m || !Gen_QQ_4mom_eta || !Gen_QQ_4mom_pt || !Gen_QQ_4mom_phi ||
            !Gen_mu_4mom_eta || !Gen_mu_4mom_pt)
          continue;
        gen_mass = Gen_QQ_4mom_m->at(irqq);
        gen_eta = Gen_QQ_4mom_eta->at(irqq);
        gen_eta1 = Gen_mu_4mom_eta->at(iMuPl);
        gen_eta2 = Gen_mu_4mom_eta->at(iMuMi);
        gen_pt = Gen_QQ_4mom_pt->at(irqq);
        gen_pt1 = Gen_mu_4mom_pt->at(iMuPl);
        gen_pt2 = Gen_mu_4mom_pt->at(iMuMi);
        gen_phi = Gen_QQ_4mom_phi->at(irqq);
        gen_phi1 = Gen_mu_4mom_phi->at(iMuPl);
        gen_phi2 = Gen_mu_4mom_phi->at(iMuMi);
        if (Gen_mu_4mom_m)
        {
          gen_m1 = Gen_mu_4mom_m->at(iMuPl);
          gen_m2 = Gen_mu_4mom_m->at(iMuMi);
        }
      }
      else
      {
        if (!Gen_QQ_4mom || !Gen_mu_4mom) continue;
        auto *genQQ = static_cast<TLorentzVector *>(Gen_QQ_4mom->At(irqq));
        auto *genMuPl = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuPl));
        auto *genMuMi = static_cast<TLorentzVector *>(Gen_mu_4mom->At(iMuMi));
        if (!genQQ || !genMuPl || !genMuMi) continue;
        gen_mass = genQQ->M();
        gen_eta = genQQ->Eta();
        gen_eta1 = genMuPl->Eta();
        gen_eta2 = genMuMi->Eta();
        gen_pt = genQQ->Pt();
        gen_pt1 = genMuPl->Pt();
        gen_pt2 = genMuMi->Pt();
        gen_phi = genQQ->Phi();
        gen_phi1 = genMuPl->Phi();
        gen_phi2 = genMuMi->Phi();
        gen_m1 = genMuPl->M();
        gen_m2 = genMuMi->M();
      }

      // reconstrcut TLorentz vector to get the rapidityu
      TLorentzVector JP;
      JP.SetPtEtaPhiM(gen_pt, gen_eta, gen_phi, gen_mass);
      float gen_y = JP.Rapidity();
      const double absGenY = std::abs(gen_y);

      if (!hasGenReference)
      {
        refGenPt = gen_pt;
        refGenAbsY = absGenY;
        refGenMuPl.SetPtEtaPhiM(gen_pt1, gen_eta1, gen_phi1, gen_m1);
        refGenMuMi.SetPtEtaPhiM(gen_pt2, gen_eta2, gen_phi2, gen_m2);
        refGenPtWeight = 1.0;
        if (isPtW)
        {
          TF1 *fPtWRef = (absGenY < 1.6) ? fPtW_mid : fPtW_fwd;
          if (!fPtWRef)
          {
            cout << "[ERROR] pT weight function is null\n";
            return;
          }
          refGenPtWeight = fPtWRef->Eval(gen_pt);
        }
        hasGenReference = true;
      }

      if (gen_mass < massLow || gen_mass > massHigh)
        continue;
      if (GetGenMuCharge(iMuPl) * GetGenMuCharge(iMuMi) > 0)
        continue;

      // ===== Efficiency denominator =====
      if (!(TMath::Abs(gen_y) < 2.4) ||
          !IsAcceptanceQQ(gen_pt1, gen_eta1) ||
          !IsAcceptanceQQ(gen_pt2, gen_eta2))
        continue;

      bool gen_inbin = false;
      if (std::abs(gen_y) < 1.6 && gen_pt >= 6.5f && gen_pt < 40.f) gen_inbin = true;
      if (std::abs(gen_y) >= 1.6 && std::abs(gen_y) < 2.4 && gen_pt >= 3.5f && gen_pt < 40.f) gen_inbin = true;
      if (!gen_inbin)
        continue;

      double Pt_weight_ = 1.0;
      if (isPtW)
      {
        TF1 *fPtW = (std::abs(gen_y) < 1.6) ? fPtW_mid : fPtW_fwd;
        if (!fPtW)
        {
          cout << "[ERROR] pT weight function is null\n";
          return;
        }
        Pt_weight_ = fPtW->Eval(gen_pt);
      }

      const double denWeightNoPt = eventWeightBase;
      const double denWeightPt = eventWeightBase * Pt_weight_;

      if (gen_pt >= 6.5f && gen_pt < 40.f)
      {
        hist_eff_den_rap0024->Fill(gen_pt, denWeightPt);
        if (absGenY < 0.6)
          hist_eff_den_rap0006->Fill(gen_pt, denWeightPt);
        else if (absGenY < 1.2)
          hist_eff_den_rap0612->Fill(gen_pt, denWeightPt);
        else if (absGenY < 1.6)
          hist_eff_den_rap1216->Fill(gen_pt, denWeightPt);
        else if (absGenY < 2.4)
          hist_eff_den_rap1624->Fill(gen_pt, denWeightPt);

        hist_eff_den_y->Fill(absGenY, denWeightPt);
        if (isPbPb)
        {
          if (Centrality < 20)
          {
            hist_eff_den_pt_cent0010->Fill(gen_pt, denWeightPt);
            if (absGenY < 1.6)
              hist_eff_den_pt_mid_cent0010->Fill(gen_pt, denWeightPt);
          }
          else if (Centrality < 60)
          {
            hist_eff_den_pt_cent1030->Fill(gen_pt, denWeightPt);
            if (absGenY < 1.6)
              hist_eff_den_pt_mid_cent1030->Fill(gen_pt, denWeightPt);
          }
          else if (Centrality < 180)
          {
            hist_eff_den_pt_cent3090->Fill(gen_pt, denWeightPt);
            if (absGenY < 1.6)
              hist_eff_den_pt_mid_cent3090->Fill(gen_pt, denWeightPt);
          }

          hist_eff_den_cent_rap0024->Fill(Centrality, denWeightPt);
          if (absGenY < 0.6)
            hist_eff_den_cent_rap0006->Fill(Centrality, denWeightPt);
          else if (absGenY < 1.2)
            hist_eff_den_cent_rap0612->Fill(Centrality, denWeightPt);
          else if (absGenY < 1.6)
            hist_eff_den_cent_rap1216->Fill(Centrality, denWeightPt);
          else if (absGenY < 2.4)
            hist_eff_den_cent_rap1624->Fill(Centrality, denWeightPt);
        }
      }
      else if (absGenY >= 1.6 && absGenY < 2.4 && gen_pt >= 3.5f)
      {
        hist_eff_den_rap1624->Fill(gen_pt, denWeightPt);
      }
      if (isPbPb && absGenY >= 1.6 && absGenY < 2.4 && gen_pt >= 3.5f && gen_pt < 40.f)
      {
        if (Centrality < 20)
          hist_eff_den_pt_fwd_cent0010->Fill(gen_pt, denWeightPt);
        else if (Centrality < 60)
          hist_eff_den_pt_fwd_cent1030->Fill(gen_pt, denWeightPt);
        else if (Centrality < 180)
          hist_eff_den_pt_fwd_cent3090->Fill(gen_pt, denWeightPt);
      }
      if (isPbPb && absGenY >= 1.6 && absGenY < 2.4 && gen_pt >= 3.5f && gen_pt < 6.5f)
        hist_eff_den_cent_rap1624_pt3565->Fill(Centrality, denWeightPt);

      if (std::abs(gen_y) < 1.6)
      {
        const int ptBin = hist_eff_den_mid->FindBin(gen_pt);
        if (ptBin > 0 && ptBin <= hist_eff_den_mid->GetNbinsX())
          hist_eff_den_mid->Fill(gen_pt, denWeightPt);
        if (isPbPb)
        {
          const int centBin = hist_eff_den_cent_mid->FindBin(Centrality);
          if (centBin > 0 && centBin <= hist_eff_den_cent_mid->GetNbinsX())
            hist_eff_den_cent_mid->Fill(Centrality, denWeightPt);
          hist_eff2d_den_mid_ptw->Fill(gen_pt, Centrality, denWeightPt);
          hist_eff2d_den_mid_noptw->Fill(gen_pt, Centrality, denWeightNoPt);
          hist_eff2d_den_mid_fine_ptw->Fill(gen_pt, Centrality, denWeightPt);
        }
      }
      else
      {
        const int ptBin = hist_eff_den_fwd->FindBin(gen_pt);
        if (ptBin > 0 && ptBin <= hist_eff_den_fwd->GetNbinsX())
          hist_eff_den_fwd->Fill(gen_pt, denWeightPt);
        if (isPbPb)
        {
          const int centBin = hist_eff_den_cent_fwd->FindBin(Centrality);
          if (centBin > 0 && centBin <= hist_eff_den_cent_fwd->GetNbinsX())
            hist_eff_den_cent_fwd->Fill(Centrality, denWeightPt);
          hist_eff2d_den_fwd_ptw->Fill(gen_pt, Centrality, denWeightPt);
          hist_eff2d_den_fwd_noptw->Fill(gen_pt, Centrality, denWeightNoPt);
          hist_eff2d_den_fwd_fine_ptw->Fill(gen_pt, Centrality, denWeightPt);
        }
      }

      // init weights
      weight = eventWeightBase;
      Double_t tnp_weight = 1.0;
      Double_t tnp_trig_w_pl = -1.0;
      Double_t tnp_trig_w_mi = -1.0;

      ++nGenDimu;
    } // end of dimuon loop

    if (nGenDimu > 0)
    {
      count_gen_dimuon += nGenDimu;
      // flowTree->Fill();
      if (VERBOSE_EVENT)
        cout << "  -> GEN accepted (" << nGenDimu << " dimuons)\n";
    }

    // ===== Reco dimuon loop =====
    nDimu = 0;
    int bestRecoIdx = -1;
    int bestRecoMuPlIdx = -1;
    int bestRecoMuMiIdx = -1;
    float bestRecoPt = 0.f;
    float bestRecoAbsY = 0.f;
    float bestRecoPt1 = 0.f;
    float bestRecoPt2 = 0.f;
    float bestRecoEta1 = 0.f;
    float bestRecoEta2 = 0.f;
    double bestRecoDr = 1e9;
    for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
    {
      if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";
      ++recoCandTotal;

      const int iMuPl = GetRecoQQMuplIdx(irqq);
      const int iMuMi = GetRecoQQMumiIdx(irqq);
      if (iMuPl < 0 || iMuMi < 0) continue;

      float mass_ = 0.f, eta_ = 0.f, eta1_ = 0.f, eta2_ = 0.f;
      float pt_ = 0.f, pt1_ = 0.f, pt2_ = 0.f;
      float phi_ = 0.f, phi1_ = 0.f, phi2_ = 0.f;
      float m1_ = 0.105658f, m2_ = 0.105658f;
      if (usesSplitMomentumBranches)
      {
        if (!Reco_QQ_4mom_m || !Reco_QQ_4mom_eta || !Reco_QQ_4mom_pt || !Reco_QQ_4mom_phi ||
            !Reco_mu_4mom_eta || !Reco_mu_4mom_pt || !Reco_mu_4mom_phi)
          continue;
        mass_ = Reco_QQ_4mom_m->at(irqq);
        eta_ = Reco_QQ_4mom_eta->at(irqq);
        eta1_ = Reco_mu_4mom_eta->at(iMuPl);
        eta2_ = Reco_mu_4mom_eta->at(iMuMi);
        pt_ = Reco_QQ_4mom_pt->at(irqq);
        pt1_ = Reco_mu_4mom_pt->at(iMuPl);
        pt2_ = Reco_mu_4mom_pt->at(iMuMi);
        phi_ = Reco_QQ_4mom_phi->at(irqq);
        phi1_ = Reco_mu_4mom_phi->at(iMuPl);
        phi2_ = Reco_mu_4mom_phi->at(iMuMi);
        if (Reco_mu_4mom_m)
        {
          m1_ = Reco_mu_4mom_m->at(iMuPl);
          m2_ = Reco_mu_4mom_m->at(iMuMi);
        }
      }
      else
      {
        if (!Reco_QQ_4mom || !Reco_mu_4mom) continue;
        auto *recoQQ = static_cast<TLorentzVector *>(Reco_QQ_4mom->At(irqq));
        auto *recoMuPl = static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuPl));
        auto *recoMuMi = static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuMi));
        if (!recoQQ || !recoMuPl || !recoMuMi) continue;
        mass_ = recoQQ->M();
        eta_ = recoQQ->Eta();
        eta1_ = recoMuPl->Eta();
        eta2_ = recoMuMi->Eta();
        pt_ = recoQQ->Pt();
        pt1_ = recoMuPl->Pt();
        pt2_ = recoMuMi->Pt();
        phi_ = recoQQ->Phi();
        phi1_ = recoMuPl->Phi();
        phi2_ = recoMuMi->Phi();
        m1_ = recoMuPl->M();
        m2_ = recoMuMi->M();
      }

      TLorentzVector JP;
      JP.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
      const float y_ = JP.Rapidity();
      const double absRecoY = std::abs(y_);

      if (mass_ < massLow || mass_ > massHigh) continue;
      ++recoPassMass;
      if (!(TMath::Abs(y_) < 2.4) ||
          !IsAcceptanceQQ(pt1_, eta1_) ||
          !IsAcceptanceQQ(pt2_, eta2_)) continue;
      ++recoPassAcceptance;
      bool rec_inbin = false;
      if (absRecoY < 1.6 && pt_ >= 6.5f && pt_ < 40.f) rec_inbin = true;
      if (absRecoY >= 1.6 && absRecoY < 2.4 && pt_ >= 3.5f && pt_ < 40.f) rec_inbin = true;
      if (!rec_inbin) continue;
      ++recoPassBin;
      if (GetRecoQQSign(irqq) != 0) continue;
      ++recoPassSign;

      const bool hasIsGoodMuonBranch = (oniaTree->GetBranch("Reco_mu_isGoodMuon") != nullptr);
      const bool hasSoftCutBranch = (oniaTree->GetBranch("Reco_mu_isSoftCutBased") != nullptr);
      const Bool_t muPlIsSoft = hasSoftCutBranch ? Reco_mu_isSoftCutBased[iMuPl] : kFALSE;
      const Bool_t muMiIsSoft = hasSoftCutBranch ? Reco_mu_isSoftCutBased[iMuMi] : kFALSE;
      const bool muPlPassType = PassMuonSelectionTypeLikeSkim(Reco_mu_SelectionType[iMuPl]);
      const bool muMiPassType = PassMuonSelectionTypeLikeSkim(Reco_mu_SelectionType[iMuMi]);
      const bool muPlQuality = PassMuonQualityLikeOniaEff(
          hasIsGoodMuonBranch, Reco_mu_isGoodMuon[iMuPl], muPlIsSoft,
          Reco_mu_SelectionType[iMuPl], Reco_mu_nTrkWMea[iMuPl],
          Reco_mu_nPixWMea[iMuPl], Reco_mu_dxy[iMuPl], Reco_mu_dz[iMuPl]);
      const bool muMiQuality = PassMuonQualityLikeOniaEff(
          hasIsGoodMuonBranch, Reco_mu_isGoodMuon[iMuMi], muMiIsSoft,
          Reco_mu_SelectionType[iMuMi], Reco_mu_nTrkWMea[iMuMi],
          Reco_mu_nPixWMea[iMuMi], Reco_mu_dxy[iMuMi], Reco_mu_dz[iMuMi]);
      if (!muPlPassType || !muMiPassType) continue;
      ++recoPassMuonType;
      if (!muPlQuality || !muMiQuality) continue;
      ++recoPassMuonQuality;
      if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;
      ++recoPassVtx;

      const bool HLTFilterPass =
          ((Reco_QQ_trig[irqq] & ((ULong64_t)std::pow(2, kTrigSel))) == ((ULong64_t)std::pow(2, kTrigSel)));
      if (!(HLTPass && HLTFilterPass)) continue;
      ++recoPassTrig;
    }

    if (hasGenReference)
    {
      for (Int_t irqq = 0; irqq < Reco_QQ_size; ++irqq)
      {
        if (VERBOSE_CANDID) cout << "  irqq: " << irqq << "\n";

        const int iMuPl = GetRecoQQMuplIdx(irqq);
        const int iMuMi = GetRecoQQMumiIdx(irqq);
        if (iMuPl < 0 || iMuMi < 0) continue;

        float mass_ = 0.f, eta_ = 0.f, eta1_ = 0.f, eta2_ = 0.f;
        float pt_ = 0.f, pt1_ = 0.f, pt2_ = 0.f;
        float phi_ = 0.f, phi1_ = 0.f, phi2_ = 0.f;
        float m1_ = 0.105658f, m2_ = 0.105658f;
        if (usesSplitMomentumBranches)
        {
          if (!Reco_QQ_4mom_m || !Reco_QQ_4mom_eta || !Reco_QQ_4mom_pt || !Reco_QQ_4mom_phi ||
              !Reco_mu_4mom_eta || !Reco_mu_4mom_pt || !Reco_mu_4mom_phi)
            continue;
          mass_ = Reco_QQ_4mom_m->at(irqq);
          eta_ = Reco_QQ_4mom_eta->at(irqq);
          eta1_ = Reco_mu_4mom_eta->at(iMuPl);
          eta2_ = Reco_mu_4mom_eta->at(iMuMi);
          pt_ = Reco_QQ_4mom_pt->at(irqq);
          pt1_ = Reco_mu_4mom_pt->at(iMuPl);
          pt2_ = Reco_mu_4mom_pt->at(iMuMi);
          phi_ = Reco_QQ_4mom_phi->at(irqq);
          phi1_ = Reco_mu_4mom_phi->at(iMuPl);
          phi2_ = Reco_mu_4mom_phi->at(iMuMi);
          if (Reco_mu_4mom_m)
          {
            m1_ = Reco_mu_4mom_m->at(iMuPl);
            m2_ = Reco_mu_4mom_m->at(iMuMi);
          }
        }
        else
        {
          if (!Reco_QQ_4mom || !Reco_mu_4mom) continue;
          auto *recoQQ = static_cast<TLorentzVector *>(Reco_QQ_4mom->At(irqq));
          auto *recoMuPl = static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuPl));
          auto *recoMuMi = static_cast<TLorentzVector *>(Reco_mu_4mom->At(iMuMi));
          if (!recoQQ || !recoMuPl || !recoMuMi) continue;
          mass_ = recoQQ->M();
          eta_ = recoQQ->Eta();
          eta1_ = recoMuPl->Eta();
          eta2_ = recoMuMi->Eta();
          pt_ = recoQQ->Pt();
          pt1_ = recoMuPl->Pt();
          pt2_ = recoMuMi->Pt();
          phi_ = recoQQ->Phi();
          phi1_ = recoMuPl->Phi();
          phi2_ = recoMuMi->Phi();
          m1_ = recoMuPl->M();
          m2_ = recoMuMi->M();
        }

        TLorentzVector JP;
        JP.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
        const float y_ = JP.Rapidity();
        const double absRecoY = std::abs(y_);

        if (mass_ < massLow || mass_ > massHigh) continue;
        if (!(TMath::Abs(y_) < 2.4) ||
            !IsAcceptanceQQ(pt1_, eta1_) ||
            !IsAcceptanceQQ(pt2_, eta2_)) continue;
        bool rec_inbin = false;
        if (absRecoY < 1.6 && pt_ >= 6.5f && pt_ < 40.f) rec_inbin = true;
        if (absRecoY >= 1.6 && absRecoY < 2.4 && pt_ >= 3.5f && pt_ < 40.f) rec_inbin = true;
        if (!rec_inbin) continue;
        if (GetRecoQQSign(irqq) != 0) continue;

        const bool hasIsGoodMuonBranch = (oniaTree->GetBranch("Reco_mu_isGoodMuon") != nullptr);
        const bool hasSoftCutBranch = (oniaTree->GetBranch("Reco_mu_isSoftCutBased") != nullptr);
        const Bool_t muPlIsSoft = hasSoftCutBranch ? Reco_mu_isSoftCutBased[iMuPl] : kFALSE;
        const Bool_t muMiIsSoft = hasSoftCutBranch ? Reco_mu_isSoftCutBased[iMuMi] : kFALSE;
        const bool muPlPassType = PassMuonSelectionTypeLikeSkim(Reco_mu_SelectionType[iMuPl]);
        const bool muMiPassType = PassMuonSelectionTypeLikeSkim(Reco_mu_SelectionType[iMuMi]);
        const bool muPlQuality = PassMuonQualityLikeOniaEff(
            hasIsGoodMuonBranch, Reco_mu_isGoodMuon[iMuPl], muPlIsSoft,
            Reco_mu_SelectionType[iMuPl], Reco_mu_nTrkWMea[iMuPl],
            Reco_mu_nPixWMea[iMuPl], Reco_mu_dxy[iMuPl], Reco_mu_dz[iMuPl]);
        const bool muMiQuality = PassMuonQualityLikeOniaEff(
            hasIsGoodMuonBranch, Reco_mu_isGoodMuon[iMuMi], muMiIsSoft,
            Reco_mu_SelectionType[iMuMi], Reco_mu_nTrkWMea[iMuMi],
            Reco_mu_nPixWMea[iMuMi], Reco_mu_dxy[iMuMi], Reco_mu_dz[iMuMi]);
        if (!muPlPassType || !muMiPassType) continue;
        if (!muPlQuality || !muMiQuality) continue;
        if (Reco_QQ_VtxProb[irqq] < 0.01f) continue;

        const bool HLTFilterPass =
            ((Reco_QQ_trig[irqq] & ((ULong64_t)std::pow(2, kTrigSel))) == ((ULong64_t)std::pow(2, kTrigSel)));
        if (!(HLTPass && HLTFilterPass)) continue;

        TLorentzVector recoMuPlVec;
        TLorentzVector recoMuMiVec;
        recoMuPlVec.SetPtEtaPhiM(pt1_, eta1_, phi1_, m1_);
        recoMuMiVec.SetPtEtaPhiM(pt2_, eta2_, phi2_, m2_);
        const double dr = std::max(recoMuPlVec.DeltaR(refGenMuPl), recoMuMiVec.DeltaR(refGenMuMi));
        if (dr >= bestRecoDr) continue;

        bestRecoIdx = irqq;
        bestRecoMuPlIdx = iMuPl;
        bestRecoMuMiIdx = iMuMi;
        bestRecoPt = pt_;
        bestRecoAbsY = absRecoY;
        bestRecoPt1 = pt1_;
        bestRecoPt2 = pt2_;
        bestRecoEta1 = eta1_;
        bestRecoEta2 = eta2_;
        bestRecoDr = dr;
      }
    }

    if (bestRecoIdx >= 0)
    {
      ++recoPassMatched;
      weight = eventWeightBase;
      Double_t tnp_weight = 1.0;
      if (isTnPW)
      {
        const bool tnpOk =
            isPbPb
                ? ComputePbPbTnPWeight(bestRecoPt1, bestRecoEta1, Reco_mu_trig[bestRecoMuPlIdx],
                                       bestRecoPt2, bestRecoEta2, Reco_mu_trig[bestRecoMuMiIdx],
                                       tnp_weight)
                : ComputePpTnPWeight(bestRecoPt1, bestRecoEta1, Reco_mu_trig[bestRecoMuPlIdx],
                                     bestRecoPt2, bestRecoEta2, Reco_mu_trig[bestRecoMuMiIdx],
                                     tnp_weight);
        if (!tnpOk)
          continue;
      }

      const double recoWeightBase = weight * refGenPtWeight;
      const double recoWeightBaseNoPt = weight;
      const double recoWeightNumNoPt = recoWeightBaseNoPt * tnp_weight;
      const double recoWeightNum = recoWeightBase * tnp_weight;
      const double recoCrossWeight = recoWeightBase * recoWeightNum;
      const bool refPassRap0024 = (refGenAbsY < 2.4f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassRap0006 = (refGenAbsY < 0.6f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassRap0612 = (refGenAbsY >= 0.6f && refGenAbsY < 1.2f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassRap1216 = (refGenAbsY >= 1.2f && refGenAbsY < 1.6f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassRap1624High = (refGenAbsY >= 1.6f && refGenAbsY < 2.4f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassRap1624Low = (refGenAbsY >= 1.6f && refGenAbsY < 2.4f && refGenPt >= 3.5f && refGenPt < 6.5f);
      const bool refPassMid = (refGenAbsY < 1.6f && refGenPt >= 6.5f && refGenPt < 40.f);
      const bool refPassFwd = (refGenAbsY >= 1.6f && refGenAbsY < 2.4f && refGenPt >= 3.5f && refGenPt < 40.f);
      if (bestRecoPt >= 6.5f && bestRecoPt < 40.f)
      {
        hist_eff_num_rap0024->Fill(refGenPt, recoWeightNum);
        hist_eff_cross_rap0024->Fill(refGenPt, recoCrossWeight);
        if (bestRecoAbsY < 0.6)
        {
          hist_eff_num_rap0006->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_rap0006->Fill(refGenPt, recoCrossWeight);
        }
        else if (bestRecoAbsY < 1.2)
        {
          hist_eff_num_rap0612->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_rap0612->Fill(refGenPt, recoCrossWeight);
        }
        else if (bestRecoAbsY < 1.6)
        {
          hist_eff_num_rap1216->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_rap1216->Fill(refGenPt, recoCrossWeight);
        }
        else if (bestRecoAbsY < 2.4)
        {
          hist_eff_num_rap1624->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_rap1624->Fill(refGenPt, recoCrossWeight);
        }

        if (refPassRap0024)
        {
          hist_eff_num_y->Fill(refGenAbsY, recoWeightNum);
          hist_eff_cross_y->Fill(refGenAbsY, recoCrossWeight);
        }
        if (isPbPb)
        {
          if (Centrality < 20)
          {
            hist_eff_num_pt_cent0010->Fill(refGenPt, recoWeightNum);
            hist_eff_cross_pt_cent0010->Fill(refGenPt, recoCrossWeight);
            if (bestRecoAbsY < 1.6)
            {
              hist_eff_num_pt_mid_cent0010->Fill(refGenPt, recoWeightNum);
              hist_eff_cross_pt_mid_cent0010->Fill(refGenPt, recoCrossWeight);
            }
          }
          else if (Centrality < 60)
          {
            hist_eff_num_pt_cent1030->Fill(refGenPt, recoWeightNum);
            hist_eff_cross_pt_cent1030->Fill(refGenPt, recoCrossWeight);
            if (bestRecoAbsY < 1.6)
            {
              hist_eff_num_pt_mid_cent1030->Fill(refGenPt, recoWeightNum);
              hist_eff_cross_pt_mid_cent1030->Fill(refGenPt, recoCrossWeight);
            }
          }
          else if (Centrality < 180)
          {
            hist_eff_num_pt_cent3090->Fill(refGenPt, recoWeightNum);
            hist_eff_cross_pt_cent3090->Fill(refGenPt, recoCrossWeight);
            if (bestRecoAbsY < 1.6)
            {
              hist_eff_num_pt_mid_cent3090->Fill(refGenPt, recoWeightNum);
              hist_eff_cross_pt_mid_cent3090->Fill(refGenPt, recoCrossWeight);
            }
          }

          if (refPassRap0024)
          {
            hist_eff_num_cent_rap0024->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_rap0024->Fill(Centrality, recoCrossWeight);
          }
          if (refPassRap0006)
          {
            hist_eff_num_cent_rap0006->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_rap0006->Fill(Centrality, recoCrossWeight);
          }
          else if (refPassRap0612)
          {
            hist_eff_num_cent_rap0612->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_rap0612->Fill(Centrality, recoCrossWeight);
          }
          else if (refPassRap1216)
          {
            hist_eff_num_cent_rap1216->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_rap1216->Fill(Centrality, recoCrossWeight);
          }
          else if (refPassRap1624High)
          {
            hist_eff_num_cent_rap1624->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_rap1624->Fill(Centrality, recoCrossWeight);
          }
        }
      }
      else if (bestRecoAbsY >= 1.6 && bestRecoAbsY < 2.4 && bestRecoPt >= 3.5f)
      {
        hist_eff_num_rap1624->Fill(refGenPt, recoWeightNum);
        hist_eff_cross_rap1624->Fill(refGenPt, recoCrossWeight);
      }
      if (isPbPb && bestRecoAbsY >= 1.6 && bestRecoAbsY < 2.4 && bestRecoPt >= 3.5f && bestRecoPt < 40.f)
      {
        if (Centrality < 20)
        {
          hist_eff_num_pt_fwd_cent0010->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_pt_fwd_cent0010->Fill(refGenPt, recoCrossWeight);
        }
        else if (Centrality < 60)
        {
          hist_eff_num_pt_fwd_cent1030->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_pt_fwd_cent1030->Fill(refGenPt, recoCrossWeight);
        }
        else if (Centrality < 180)
        {
          hist_eff_num_pt_fwd_cent3090->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_pt_fwd_cent3090->Fill(refGenPt, recoCrossWeight);
        }
      }
      if (isPbPb && bestRecoAbsY >= 1.6 && bestRecoAbsY < 2.4 && bestRecoPt >= 3.5f && bestRecoPt < 6.5f && refPassRap1624Low)
      {
        hist_eff_num_cent_rap1624_pt3565->Fill(Centrality, recoWeightNum);
        hist_eff_cross_cent_rap1624_pt3565->Fill(Centrality, recoCrossWeight);
      }

      if (bestRecoAbsY < 1.6)
      {
        const int ptBin = hist_eff_num_mid->FindBin(refGenPt);
        if (ptBin > 0 && ptBin <= hist_eff_num_mid->GetNbinsX())
        {
          hist_eff_num_mid->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_mid->Fill(refGenPt, recoCrossWeight);
        }
        if (isPbPb)
        {
          const int centBin = hist_eff_num_cent_mid->FindBin(Centrality);
          if (refPassMid && centBin > 0 && centBin <= hist_eff_num_cent_mid->GetNbinsX())
          {
            hist_eff_num_cent_mid->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_mid->Fill(Centrality, recoCrossWeight);
          }
          hist_eff2d_num_mid_ptw->Fill(refGenPt, Centrality, recoWeightNum);
          hist_eff2d_num_mid_noptw->Fill(refGenPt, Centrality, recoWeightNumNoPt);
          hist_eff2d_num_mid_fine_ptw->Fill(refGenPt, Centrality, recoWeightNum);
        }
      }
      else if (bestRecoAbsY < 2.4)
      {
        const int ptBin = hist_eff_num_fwd->FindBin(refGenPt);
        if (ptBin > 0 && ptBin <= hist_eff_num_fwd->GetNbinsX())
        {
          hist_eff_num_fwd->Fill(refGenPt, recoWeightNum);
          hist_eff_cross_fwd->Fill(refGenPt, recoCrossWeight);
        }
        if (isPbPb)
        {
          const int centBin = hist_eff_num_cent_fwd->FindBin(Centrality);
          if (refPassFwd && centBin > 0 && centBin <= hist_eff_num_cent_fwd->GetNbinsX())
          {
            hist_eff_num_cent_fwd->Fill(Centrality, recoWeightNum);
            hist_eff_cross_cent_fwd->Fill(Centrality, recoCrossWeight);
          }
          hist_eff2d_num_fwd_ptw->Fill(refGenPt, Centrality, recoWeightNum);
          hist_eff2d_num_fwd_noptw->Fill(refGenPt, Centrality, recoWeightNumNoPt);
          hist_eff2d_num_fwd_fine_ptw->Fill(refGenPt, Centrality, recoWeightNum);
        }
      }

      ++nDimu;
    }

    if (nDimu > 0)
    {
      count_reco_dimuon += nDimu;
      // flowTree->Fill();
      if (VERBOSE_EVENT) cout << "  -> Fill (" << nDimu << " dimuons)\n";
    }
  }

  // ===== compute efficiency =====
  hist_eff_den_mid->Scale(1.0, "width");
  hist_eff_den_fwd->Scale(1.0, "width");
  hist_eff_num_mid->Scale(1.0, "width");
  hist_eff_num_fwd->Scale(1.0, "width");
  hist_eff_den_cent_mid->Scale(1.0, "width");
  hist_eff_den_cent_fwd->Scale(1.0, "width");
  hist_eff_num_cent_mid->Scale(1.0, "width");
  hist_eff_num_cent_fwd->Scale(1.0, "width");
  hist_eff_den_pt_cent0010->Scale(1.0, "width");
  hist_eff_den_pt_cent1030->Scale(1.0, "width");
  hist_eff_den_pt_cent3090->Scale(1.0, "width");
  hist_eff_den_pt_mid_cent0010->Scale(1.0, "width");
  hist_eff_den_pt_mid_cent1030->Scale(1.0, "width");
  hist_eff_den_pt_mid_cent3090->Scale(1.0, "width");
  hist_eff_den_pt_fwd_cent0010->Scale(1.0, "width");
  hist_eff_den_pt_fwd_cent1030->Scale(1.0, "width");
  hist_eff_den_pt_fwd_cent3090->Scale(1.0, "width");
  hist_eff_num_pt_cent0010->Scale(1.0, "width");
  hist_eff_num_pt_cent1030->Scale(1.0, "width");
  hist_eff_num_pt_cent3090->Scale(1.0, "width");
  hist_eff_num_pt_mid_cent0010->Scale(1.0, "width");
  hist_eff_num_pt_mid_cent1030->Scale(1.0, "width");
  hist_eff_num_pt_mid_cent3090->Scale(1.0, "width");
  hist_eff_num_pt_fwd_cent0010->Scale(1.0, "width");
  hist_eff_num_pt_fwd_cent1030->Scale(1.0, "width");
  hist_eff_num_pt_fwd_cent3090->Scale(1.0, "width");
  hist_eff_den_rap0024->Scale(1.0, "width");
  hist_eff_den_rap0006->Scale(1.0, "width");
  hist_eff_den_rap0612->Scale(1.0, "width");
  hist_eff_den_rap1216->Scale(1.0, "width");
  hist_eff_den_rap1624->Scale(1.0, "width");
  hist_eff_den_cent_rap0024->Scale(1.0, "width");
  hist_eff_den_cent_rap0006->Scale(1.0, "width");
  hist_eff_den_cent_rap0612->Scale(1.0, "width");
  hist_eff_den_cent_rap1216->Scale(1.0, "width");
  hist_eff_den_cent_rap1624->Scale(1.0, "width");
  hist_eff_den_cent_rap1624_pt3565->Scale(1.0, "width");
  hist_eff_num_rap0024->Scale(1.0, "width");
  hist_eff_num_rap0006->Scale(1.0, "width");
  hist_eff_num_rap0612->Scale(1.0, "width");
  hist_eff_num_rap1216->Scale(1.0, "width");
  hist_eff_num_rap1624->Scale(1.0, "width");
  hist_eff_num_cent_rap0024->Scale(1.0, "width");
  hist_eff_num_cent_rap0006->Scale(1.0, "width");
  hist_eff_num_cent_rap0612->Scale(1.0, "width");
  hist_eff_num_cent_rap1216->Scale(1.0, "width");
  hist_eff_num_cent_rap1624->Scale(1.0, "width");
  hist_eff_num_cent_rap1624_pt3565->Scale(1.0, "width");
  hist_eff_den_y->Scale(1.0, "width");
  hist_eff_num_y->Scale(1.0, "width");

  // Apply identical width normalization to numerator and denominator before
  // building the efficiency histograms.
  TH1D *hist_eff_mid = BuildWeightedEfficiencyHist(hist_eff_num_mid, hist_eff_den_mid, "hist_eff_mid", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_fwd = BuildWeightedEfficiencyHist(hist_eff_num_fwd, hist_eff_den_fwd, "hist_eff_fwd", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_cent_mid = BuildWeightedEfficiencyHist(hist_eff_num_cent_mid, hist_eff_den_cent_mid, "hist_eff_cent_mid", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_fwd = BuildWeightedEfficiencyHist(hist_eff_num_cent_fwd, hist_eff_den_cent_fwd, "hist_eff_cent_fwd", ";Centrality;Efficiency");
  TH1D *hist_eff_pt_cent0010 = BuildWeightedEfficiencyHist(hist_eff_num_pt_cent0010, hist_eff_den_pt_cent0010, "hist_eff_pt_cent0010", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_cent1030 = BuildWeightedEfficiencyHist(hist_eff_num_pt_cent1030, hist_eff_den_pt_cent1030, "hist_eff_pt_cent1030", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_cent3090 = BuildWeightedEfficiencyHist(hist_eff_num_pt_cent3090, hist_eff_den_pt_cent3090, "hist_eff_pt_cent3090", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_mid_cent0010 = BuildWeightedEfficiencyHist(hist_eff_num_pt_mid_cent0010, hist_eff_den_pt_mid_cent0010, "hist_eff_pt_mid_cent0010", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_mid_cent1030 = BuildWeightedEfficiencyHist(hist_eff_num_pt_mid_cent1030, hist_eff_den_pt_mid_cent1030, "hist_eff_pt_mid_cent1030", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_mid_cent3090 = BuildWeightedEfficiencyHist(hist_eff_num_pt_mid_cent3090, hist_eff_den_pt_mid_cent3090, "hist_eff_pt_mid_cent3090", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_fwd_cent0010 = BuildWeightedEfficiencyHist(hist_eff_num_pt_fwd_cent0010, hist_eff_den_pt_fwd_cent0010, "hist_eff_pt_fwd_cent0010", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_fwd_cent1030 = BuildWeightedEfficiencyHist(hist_eff_num_pt_fwd_cent1030, hist_eff_den_pt_fwd_cent1030, "hist_eff_pt_fwd_cent1030", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_pt_fwd_cent3090 = BuildWeightedEfficiencyHist(hist_eff_num_pt_fwd_cent3090, hist_eff_den_pt_fwd_cent3090, "hist_eff_pt_fwd_cent3090", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_rap0024 = BuildWeightedEfficiencyHist(hist_eff_num_rap0024, hist_eff_den_rap0024, "hist_eff_rap0024", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_rap0006 = BuildWeightedEfficiencyHist(hist_eff_num_rap0006, hist_eff_den_rap0006, "hist_eff_rap0006", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_rap0612 = BuildWeightedEfficiencyHist(hist_eff_num_rap0612, hist_eff_den_rap0612, "hist_eff_rap0612", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_rap1216 = BuildWeightedEfficiencyHist(hist_eff_num_rap1216, hist_eff_den_rap1216, "hist_eff_rap1216", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_rap1624 = BuildWeightedEfficiencyHist(hist_eff_num_rap1624, hist_eff_den_rap1624, "hist_eff_rap1624", ";p_{T} (GeV/c);Efficiency");
  TH1D *hist_eff_cent_rap0024 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap0024, hist_eff_den_cent_rap0024, "hist_eff_cent_rap0024", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_rap0006 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap0006, hist_eff_den_cent_rap0006, "hist_eff_cent_rap0006", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_rap0612 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap0612, hist_eff_den_cent_rap0612, "hist_eff_cent_rap0612", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_rap1216 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap1216, hist_eff_den_cent_rap1216, "hist_eff_cent_rap1216", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_rap1624 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap1624, hist_eff_den_cent_rap1624, "hist_eff_cent_rap1624", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_rap1624_pt3565 = BuildWeightedEfficiencyHist(hist_eff_num_cent_rap1624_pt3565, hist_eff_den_cent_rap1624_pt3565, "hist_eff_cent_rap1624_pt3565", ";Centrality;Efficiency");
  TH1D *hist_eff_y = BuildCorrelatedEfficiencyHist(hist_eff_num_y, hist_eff_den_y, hist_eff_cross_y, "hist_eff_y", ";|y|;Efficiency");
  TH1D *hist_eff_den_cent_mid_stdpt = ProjectCentHistFrom2D(hist_eff2d_den_mid_ptw, "hist_eff_den_cent_mid_stdpt", ";Centrality;Counts");
  TH1D *hist_eff_num_cent_mid_stdpt = ProjectCentHistFrom2D(hist_eff2d_num_mid_ptw, "hist_eff_num_cent_mid_stdpt", ";Centrality;Counts");
  TH1D *hist_eff_den_cent_fwd_stdpt = ProjectCentHistFrom2D(hist_eff2d_den_fwd_ptw, "hist_eff_den_cent_fwd_stdpt", ";Centrality;Counts");
  TH1D *hist_eff_num_cent_fwd_stdpt = ProjectCentHistFrom2D(hist_eff2d_num_fwd_ptw, "hist_eff_num_cent_fwd_stdpt", ";Centrality;Counts");
  TH1D *hist_eff_cent_mid_stdpt = BuildWeightedEfficiencyHist(hist_eff_num_cent_mid_stdpt, hist_eff_den_cent_mid_stdpt, "hist_eff_cent_mid_stdpt", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_fwd_stdpt = BuildWeightedEfficiencyHist(hist_eff_num_cent_fwd_stdpt, hist_eff_den_cent_fwd_stdpt, "hist_eff_cent_fwd_stdpt", ";Centrality;Efficiency");
  TH1D *hist_eff_den_cent_mid_0180 = CollapseToSingleBinHist(hist_eff_den_cent_mid_stdpt, "hist_eff_den_cent_mid_0180", ";Centrality;Counts", 0.0, 180.0);
  TH1D *hist_eff_num_cent_mid_0180 = CollapseToSingleBinHist(hist_eff_num_cent_mid_stdpt, "hist_eff_num_cent_mid_0180", ";Centrality;Counts", 0.0, 180.0);
  TH1D *hist_eff_den_cent_fwd_0180 = CollapseToSingleBinHist(hist_eff_den_cent_fwd_stdpt, "hist_eff_den_cent_fwd_0180", ";Centrality;Counts", 0.0, 180.0);
  TH1D *hist_eff_num_cent_fwd_0180 = CollapseToSingleBinHist(hist_eff_num_cent_fwd_stdpt, "hist_eff_num_cent_fwd_0180", ";Centrality;Counts", 0.0, 180.0);
  TH1D *hist_eff_cent_mid_0180 = BuildWeightedEfficiencyHist(hist_eff_num_cent_mid_0180, hist_eff_den_cent_mid_0180, "hist_eff_cent_mid_0180", ";Centrality;Efficiency");
  TH1D *hist_eff_cent_fwd_0180 = BuildWeightedEfficiencyHist(hist_eff_num_cent_fwd_0180, hist_eff_den_cent_fwd_0180, "hist_eff_cent_fwd_0180", ";Centrality;Efficiency");
  TH2D *hist_eff2d_mid_ptw = BuildWeightedEfficiencyHist(hist_eff2d_num_mid_ptw, hist_eff2d_den_mid_ptw, "hist_eff2d_mid_ptw", ";p_{T} (GeV/c);Centrality;Efficiency");
  TH2D *hist_eff2d_mid_noptw = BuildWeightedEfficiencyHist(hist_eff2d_num_mid_noptw, hist_eff2d_den_mid_noptw, "hist_eff2d_mid_noptw", ";p_{T} (GeV/c);Centrality;Efficiency");
  TH2D *hist_eff2d_fwd_ptw = BuildWeightedEfficiencyHist(hist_eff2d_num_fwd_ptw, hist_eff2d_den_fwd_ptw, "hist_eff2d_fwd_ptw", ";p_{T} (GeV/c);Centrality;Efficiency");
  TH2D *hist_eff2d_fwd_noptw = BuildWeightedEfficiencyHist(hist_eff2d_num_fwd_noptw, hist_eff2d_den_fwd_noptw, "hist_eff2d_fwd_noptw", ";p_{T} (GeV/c);Centrality;Efficiency");
  TH2D *hist_eff2d_mid_fine_ptw = BuildWeightedEfficiencyHist(hist_eff2d_num_mid_fine_ptw, hist_eff2d_den_mid_fine_ptw, "hist_eff2d_mid_fine_ptw", ";p_{T} (GeV/c);Centrality;Efficiency");
  TH2D *hist_eff2d_fwd_fine_ptw = BuildWeightedEfficiencyHist(hist_eff2d_num_fwd_fine_ptw, hist_eff2d_den_fwd_fine_ptw, "hist_eff2d_fwd_fine_ptw", ";p_{T} (GeV/c);Centrality;Efficiency");

  // ===== save results =====
  fFlowSkim->cd();
  hist_eff_den_mid->Write();
  hist_eff_den_fwd->Write();
  hist_eff_num_mid->Write();
  hist_eff_num_fwd->Write();
  hist_eff_den_cent_mid->Write();
  hist_eff_den_cent_fwd->Write();
  hist_eff_num_cent_mid->Write();
  hist_eff_num_cent_fwd->Write();
  hist_eff_den_pt_cent0010->Write();
  hist_eff_den_pt_cent1030->Write();
  hist_eff_den_pt_cent3090->Write();
  hist_eff_den_pt_mid_cent0010->Write();
  hist_eff_den_pt_mid_cent1030->Write();
  hist_eff_den_pt_mid_cent3090->Write();
  hist_eff_den_pt_fwd_cent0010->Write();
  hist_eff_den_pt_fwd_cent1030->Write();
  hist_eff_den_pt_fwd_cent3090->Write();
  hist_eff_num_pt_cent0010->Write();
  hist_eff_num_pt_cent1030->Write();
  hist_eff_num_pt_cent3090->Write();
  hist_eff_num_pt_mid_cent0010->Write();
  hist_eff_num_pt_mid_cent1030->Write();
  hist_eff_num_pt_mid_cent3090->Write();
  hist_eff_num_pt_fwd_cent0010->Write();
  hist_eff_num_pt_fwd_cent1030->Write();
  hist_eff_num_pt_fwd_cent3090->Write();
  hist_eff_den_rap0024->Write();
  hist_eff_den_rap0006->Write();
  hist_eff_den_rap0612->Write();
  hist_eff_den_rap1216->Write();
  hist_eff_den_rap1624->Write();
  hist_eff_den_cent_rap0024->Write();
  hist_eff_den_cent_rap0006->Write();
  hist_eff_den_cent_rap0612->Write();
  hist_eff_den_cent_rap1216->Write();
  hist_eff_den_cent_rap1624->Write();
  hist_eff_den_cent_rap1624_pt3565->Write();
  hist_eff_num_rap0024->Write();
  hist_eff_num_rap0006->Write();
  hist_eff_num_rap0612->Write();
  hist_eff_num_rap1216->Write();
  hist_eff_num_rap1624->Write();
  hist_eff_num_cent_rap0024->Write();
  hist_eff_num_cent_rap0006->Write();
  hist_eff_num_cent_rap0612->Write();
  hist_eff_num_cent_rap1216->Write();
  hist_eff_num_cent_rap1624->Write();
  hist_eff_num_cent_rap1624_pt3565->Write();
  hist_eff_den_y->Write();
  hist_eff_num_y->Write();
  hist_eff_den_cent_mid_stdpt->Write();
  hist_eff_num_cent_mid_stdpt->Write();
  hist_eff_den_cent_fwd_stdpt->Write();
  hist_eff_num_cent_fwd_stdpt->Write();
  hist_eff_den_cent_mid_0180->Write();
  hist_eff_num_cent_mid_0180->Write();
  hist_eff_den_cent_fwd_0180->Write();
  hist_eff_num_cent_fwd_0180->Write();
  hist_eff2d_den_mid_ptw->Write();
  hist_eff2d_den_mid_noptw->Write();
  hist_eff2d_num_mid_ptw->Write();
  hist_eff2d_num_mid_noptw->Write();
  hist_eff2d_den_fwd_ptw->Write();
  hist_eff2d_den_fwd_noptw->Write();
  hist_eff2d_num_fwd_ptw->Write();
  hist_eff2d_num_fwd_noptw->Write();
  hist_eff2d_den_mid_fine_ptw->Write();
  hist_eff2d_num_mid_fine_ptw->Write();
  hist_eff2d_den_fwd_fine_ptw->Write();
  hist_eff2d_num_fwd_fine_ptw->Write();
  hist_eff_mid->Write();
  hist_eff_fwd->Write();
  hist_eff_cent_mid->Write();
  hist_eff_cent_fwd->Write();
  hist_eff_pt_cent0010->Write();
  hist_eff_pt_cent1030->Write();
  hist_eff_pt_cent3090->Write();
  hist_eff_pt_mid_cent0010->Write();
  hist_eff_pt_mid_cent1030->Write();
  hist_eff_pt_mid_cent3090->Write();
  hist_eff_pt_fwd_cent0010->Write();
  hist_eff_pt_fwd_cent1030->Write();
  hist_eff_pt_fwd_cent3090->Write();
  hist_eff_rap0024->Write();
  hist_eff_rap0006->Write();
  hist_eff_rap0612->Write();
  hist_eff_rap1216->Write();
  hist_eff_rap1624->Write();
  hist_eff_cent_rap0024->Write();
  hist_eff_cent_rap0006->Write();
  hist_eff_cent_rap0612->Write();
  hist_eff_cent_rap1216->Write();
  hist_eff_cent_rap1624->Write();
  hist_eff_cent_rap1624_pt3565->Write();
  hist_eff_y->Write();
  hist_eff_cent_mid_stdpt->Write();
  hist_eff_cent_fwd_stdpt->Write();
  hist_eff_cent_mid_0180->Write();
  hist_eff_cent_fwd_0180->Write();
  hist_eff2d_mid_ptw->Write();
  hist_eff2d_mid_noptw->Write();
  hist_eff2d_fwd_ptw->Write();
  hist_eff2d_fwd_noptw->Write();
  hist_eff2d_mid_fine_ptw->Write();
  hist_eff2d_fwd_fine_ptw->Write();
  // flowTree->Write("myTree");
  fFlowSkim->Close();

  // ----- print -----
  cout << "Total GEN accepted dimuon: " << count_gen_dimuon << "\n";
  cout << "Total Reco selected dimuon: " << count_reco_dimuon << "\n";
  cout << "Reco candidate cutflow: total=" << recoCandTotal
       << ", mass=" << recoPassMass
       << ", acceptance=" << recoPassAcceptance
       << ", inbin=" << recoPassBin
       << ", sign=" << recoPassSign
       << ", muonType=" << recoPassMuonType
       << ", muonQuality=" << recoPassMuonQuality
       << ", vtx=" << recoPassVtx
       << ", trig=" << recoPassTrig
       << ", matched=" << recoPassMatched << "\n";
  cout << "Output ROOT: " << outFilePath << "\n";
  // cout << "Total dimuon (GEN+Reco): " << (count_gen_dimuon + count_reco_dimuon) << "\n";

  cout << "Finish onia_to_skim_data()\n";
}
