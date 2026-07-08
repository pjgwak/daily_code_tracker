/* check 할 거.
--- Current cut ---
opposite-sign dimuon:
Reco_QQ_sign[iQQ] != 0 continue
dimuon |eta| < 2.4
daughter track 제거:
  |deta| < 0.03, |dphi| < 0.03 -> 0.001로 수정?
tracker Q-vector track:
|eta| <= 2.4, 0.3 < pt < 3


--- Exist in the referecne but not this code ---
trigger cut:csTree.trigHLT()[10] -> 이건 MinBias에 맞게 수정해서 apply?
event selection:csTree.evtSel()[2] -> Oniatree EventSelection에 있으니까 생략?
csTree.evtSel()[3]

centrality:100 <= centrality < 160 -> Need

J/psi candidate condition:softCand -> cutBasedSoftMuonID? -> Track base?
|y| < 2.4 -> basic
2.5 < mass < 4.3 -> Maybe...
0.2 <= pt < 20 -> Why do we need it? We want to see up to 50 GeV

track quality: highPurity
ptError/pt < 0.10
|dz/dzerr| <= 3
|dxy/dxyerr| <= 3
*/

// ===== epAng array information =====
// Index     Name   Detector Order hmin1 hmax1 hmin2 hmax2 minpt maxpt nsub mcw    rmate1    rmate2
//     0      HFm2        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp2 trackmid2
//     1      HFp2        HF     2  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2 trackmid2
//     2       HF2        HF     2 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm2   trackp2
//     3 trackmid2   Tracker     2 -0.50  0.50  0.00  0.00  0.50  3.00 3sub  no      HFm2      HFp2
//     4   trackm2   Tracker     2 -2.00 -1.00  0.00  0.00  0.50  3.00 3sub  no      HFp2      HFm2
//     5   trackp2   Tracker     2  1.00  2.00  0.00  0.00  0.50  3.00 3sub  no      HFp2      HFm2
//     6      HFm3        HF     3 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp3 trackmid3
//     7      HFp3        HF     3  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm3 trackmid3
//     8       HF3        HF     3 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm3   trackp3
//     9 trackmid3   Tracker     3 -0.50  0.50  0.00  0.00  0.50  3.00 3sub  no      HFm3      HFp3
//    10   trackm3   Tracker     3 -2.00 -1.00  0.00  0.00  0.50  3.00 3sub  no      HFp3      HFm3
//    11   trackp3   Tracker     3  1.00  2.00  0.00  0.00  0.50  3.00 3sub  no      HFp3      HFm3

// ===== PbPb2023 Data(MC) triggers =====
// 'DoubleMuonTrigger' : cms.vstring(
//   "HLT_HIL1DoubleMu0_MaxDr3p5_Open_v",#0
//   "HLT_HIL1DoubleMu0_v",#1
//   "HLT_HIL1DoubleMu0_SQ_v",#2
//   "HLT_HIL2DoubleMu0_Open_v",#3
//   "HLT_HIL2DoubleMu0_M1p5to6_Open_v",#4
//   "HLT_HIL2DoubleMu2p8_M1p5to6_Open_v",#5
//   "HLT_HIL2DoubleMu0_M7to15_Open_v",#6
//   "HLT_HIL2DoubleMu3_M7to15_Open_v",#7
//   "HLT_HIL3DoubleMu0_M0toInf_Open_v",#8
//   "HLT_HIL3DoubleMu0_Quarkonia_Open_v",#9
//   "HLT_HIL3DoubleMu2_Quarkonia_Open_v",#10
//   "HLT_HIL3DoubleMu0_M2to4p5_Open_v",#11
//   "HLT_HIL3DoubleMu2_M2to4p5_Open_v",#12
//   "HLT_HIL3DoubleMu0_M7to15_Open_v",#13
//   "HLT_HIL3DoubleMu2_M7to15_Open_v",#14
// ),
// # Single Muon Trigger List
//   'SingleMuonTrigger' : cms.vstring(
//   "HLT_HIL1SingleMu0_Open_v",#15
//   "HLT_HIL1SingleMu0_v",#16
//   "HLT_HIL2SingleMu3_Open_v",#17
//   "HLT_HIL2SingleMu5_v",#18
//   "HLT_HIL2SingleMu7_v",#19
//   "HLT_HIL3SingleMu3_Open_v",#20
//   "HLT_HIL3SingleMu5_v",#21
//   "HLT_HIL3SingleMu7_v",#22
//   "HLT_HIL3SingleMu12_v",#23
//   "HLT_HIMinimumBiasHF1AND_v", #24
//   "HLT_HIMinimumBiasHF1ANDZDC2nOR_v", #25
//   "HLT_HIMinimumBiasHF1ANDZDC1nOR_v", #26
// )

// header
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranchObject.h"
#include "math.h"
#include <vector>
#include <TVector3.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TParticlePDG.h>
#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include <typeinfo>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TSystem.h> // gSystem
#include <algorithm>
#include <cctype>
#include <cstdio> // std::remove
#include <cstring>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

using std::cout;
using std::endl;
using std::string;
using std::vector;

#define MAXCAN 10000

// Return true only for ROOT files. This lets the directory scan ignore
// production subfolders such as 0000/0001 when checking file names.
bool HasRootSuffix(const char *fileName)
{
  const size_t len = std::strlen(fileName);
  return len > 5 && std::strcmp(fileName + len - 5, ".root") == 0;
}

// Add all .root files directly inside one directory to inputFiles.
// This function is intentionally non-recursive; GetInputRootFiles controls
// exactly which directory levels are scanned.
void AddRootFilesInDir(const TString &dirName, vector<TString> &inputFiles)
{
  DIR *dir = opendir(dirName.Data());
  if (!dir)
    return;

  dirent *entry = nullptr;
  while ((entry = readdir(dir)) != nullptr)
  {
    const char *fileName = entry->d_name;
    if (HasRootSuffix(fileName))
      inputFiles.push_back(TString::Format("%s/%s", dirName.Data(), fileName));
  }
  closedir(dir);
}

// Build the real input file list from the production base folder.
// The onia output can be split as base_folder/0000, base_folder/0001, ...;
// scanning actual files avoids assuming continuous Oniatree_MC_miniAOD_N.root indices.
vector<TString> GetInputRootFiles(const TString &base_folder)
{
  vector<TString> inputFiles;
  DIR *dir = opendir(base_folder.Data());
  if (!dir)
  {
    std::cerr << "Cannot open base folder: " << base_folder << std::endl;
    return inputFiles;
  }

  // Also support the case where ROOT files are directly under base_folder.
  AddRootFilesInDir(base_folder, inputFiles);

  dirent *entry = nullptr;
  while ((entry = readdir(dir)) != nullptr)
  {
    const char *dirEntryName = entry->d_name;
    if (std::strcmp(dirEntryName, ".") == 0 || std::strcmp(dirEntryName, "..") == 0)
      continue;
    if (HasRootSuffix(dirEntryName))
      continue;

    // Scan one production subfolder level, e.g. 250828_074840/0000 and /0001.
    TString childDir = TString::Format("%s/%s", base_folder.Data(), dirEntryName);
    AddRootFilesInDir(childDir, inputFiles);
  }
  closedir(dir);

  std::sort(inputFiles.begin(), inputFiles.end());
  return inputFiles;
}

void PrintStopwatch(const char *label, TStopwatch &timer)
{
  timer.Stop();
  cout << "[Time] " << label
       << ": real " << timer.RealTime() << " s"
       << ", cpu " << timer.CpuTime() << " s" << endl;
}

std::string SecondsToHms(double sec)
{
  int h = int(sec / 3600);
  sec -= 3600 * h;
  int m = int(sec / 60);
  sec -= 60 * m;
  int s = int(sec + 0.5);
  std::ostringstream os;
  os << std::setfill('0') << std::setw(2) << h << ":"
     << std::setw(2) << m << ":"
     << std::setw(2) << s;
  return os.str();
}

Long64_t GetProgressReportEvery(Long64_t totalEvents)
{
  if (totalEvents >= 10000000)
    return 200000;
  if (totalEvents >= 1000000)
    return 100000;
  if (totalEvents >= 100000)
    return 20000;
  return 10000;
}

void PrintLoopProgress(const char *loopLabel, Long64_t index, Long64_t totalEvents, Long64_t reportEvery, TStopwatch &timer)
{
  if ((index % reportEvery) != 0 && (index + 1) != totalEvents)
    return;

  const double elapsed = timer.RealTime();
  timer.Continue();
  const double frac = (totalEvents > 0) ? double(index + 1) / double(totalEvents) : 0.0;
  const double estTotal = (frac > 0) ? (elapsed / frac) : 0.0;
  const double eta = estTotal - elapsed;

  cout << "[" << loopLabel << " "
       << std::fixed << std::setprecision(1)
       << (frac * 100.0) << "%] "
       << (index + 1) << " / " << totalEvents
       << "  | elapsed " << SecondsToHms(elapsed)
       << "  | ETA " << SecondsToHms(std::max(0.0, eta))
       << endl;
}

TString MakeSafeNameTag(const char *nameTag)
{
  TString tag = (nameTag && std::strlen(nameTag) > 0) ? nameTag : "job";
  std::string safeTag;
  safeTag.reserve(tag.Length());
  for (int i = 0; i < tag.Length(); ++i)
  {
    const unsigned char c = tag[i];
    if (std::isalnum(c) || c == '_' || c == '-' || c == '.')
      safeTag.push_back(c);
    else
      safeTag.push_back('_');
  }
  return TString(safeTag);
}

double DeltaPhiPeriodic(double phi1, double phi2)
{
  double dphi = phi1 - phi2;
  while (dphi > TMath::Pi())
    dphi -= 2.0 * TMath::Pi();
  while (dphi <= -TMath::Pi())
    dphi += 2.0 * TMath::Pi();
  return dphi;
}

double SafeResolution(double numerator, double denominator)
{
  if (denominator <= 0.0 || numerator / denominator < 0.0)
    return -1.0;
  return std::sqrt(numerator / denominator);
}

void AccumulateTrackQ(float pt, float phi, float &qxNum, float &qyNum, float &qWeight)
{
  qxNum += pt * cos(2.0 * phi);
  qyNum += pt * sin(2.0 * phi);
  qWeight += pt;
}

void FinalizeTrackQ(float qxNum, float qyNum, float qWeight,
                    Float_t &qxTree, Float_t &qyTree, Float_t &epOut,
                    Int_t &hasEp, TH2F *qHist)
{
  if (qWeight > 0)
  {
    hasEp = 1;
    qxTree = qxNum / qWeight;
    qyTree = qyNum / qWeight;
    epOut = 0.5 * atan2(qyTree, qxTree);
    if (qHist)
      qHist->Fill(qxTree, qyTree);
  }
  else
  {
    hasEp = 0;
    qxTree = 0.f;
    qyTree = 0.f;
    epOut = -99.f;
  }
}

Float_t RecenterTrackEp(Float_t qxTree, Float_t qyTree, double meanQx, double meanQy)
{
  return 0.5 * atan2(qyTree - meanQy, qxTree - meanQx);
}

void FillFlatteningMoments(Float_t epRec, TH1F *meansSin[], TH1F *meansCos[], int nMax)
{
  for (int n = 0; n < nMax; n++)
  {
    const int j = n + 1;
    meansSin[n]->Fill(sin(2.0 * j * epRec));
    meansCos[n]->Fill(cos(2.0 * j * epRec));
  }
}

void ReadFlatteningMeans(TH1F *meansSin[], TH1F *meansCos[], double meanSin[], double meanCos[], int nMax)
{
  for (int n = 0; n < nMax; n++)
  {
    meanSin[n] = meansSin[n]->GetMean();
    meanCos[n] = meansCos[n]->GetMean();
  }
}

Float_t FlattenTrackEp(Float_t epRec, const double meanSin[], const double meanCos[], int nMax)
{
  double deltaPsi2 = 0.0;
  for (int n = 0; n < nMax; n++)
  {
    const int j = n + 1;
    deltaPsi2 += (2.0 / j) * (-meanSin[n] * cos(2.0 * j * epRec) + meanCos[n] * sin(2.0 * j * epRec));
  }
  return epRec + 0.5 * deltaPsi2;
}

const int kNCentFineBins = 18;
const int kNCentIntegratedBins = 6;
const int kNCentBins = kNCentFineBins + kNCentIntegratedBins;
const int kCentLow[kNCentBins] = {
    0, 10, 20, 30, 40, 50, 60, 70, 80,
    90, 100, 110, 120, 130, 140, 150, 160, 170,
    0, 0, 20, 40, 20, 80};
const int kCentHigh[kNCentBins] = {
    10, 20, 30, 40, 50, 60, 70, 80, 90,
    100, 110, 120, 130, 140, 150, 160, 170, 180,
    20, 40, 40, 80, 80, 180};
const char *kCentLabel[kNCentBins] = {
    "0-5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%", "30-35%", "35-40%", "40-45%",
    "45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%", "75-80%", "80-85%", "85-90%",
    "0-10%", "0-20%", "10-20%", "20-40%", "10-40%", "40-90%"};

bool IsInCentBin(Int_t centrality, Int_t centBin)
{
  if (centBin < 0 || centBin >= kNCentBins)
    return false;
  return centrality >= kCentLow[centBin] && centrality < kCentHigh[centBin];
}

Int_t GetCentBin(Int_t centrality)
{
  for (int ic = 0; ic < kNCentBins; ++ic)
    if (IsInCentBin(centrality, ic))
      return ic;
  return -1;
}

TString GetCentBinLabel(Int_t centBin)
{
  if (centBin < 0 || centBin >= kNCentBins)
    return "invalid";
  return kCentLabel[centBin];
}

void SetupCentralityHist(TH1F *h)
{
  if (!h)
    return;
  for (int i = 0; i < kNCentBins; ++i)
    h->GetXaxis()->SetBinLabel(i + 1, GetCentBinLabel(i));
}

void getEpAngles_calcTrack_3subEvent(bool isMC = true, bool isPr = true, Long64_t maxEvents = 1000000, const char *outDirOverride = "", const char *nameTag = "", bool useSingleFileInput = true)
{
  // ===== input files ===== //
  // gErrorIgnoreLevel = kError;

  TStopwatch totalTimer;
  totalTimer.Start();

  // inputs
  TString prompt_base_folder = "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Prompt_with_Track/OniaTree_MC1_Run3_PbPb_track_vector_PR/PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/OniaTree_MC1_Run3_PbPb_track_vector_PR/250828_074840";
  TString nonprompt_base_folder = "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Nonprompt_with_Track/OniaTree_MC1_Run3_PbPb_track_vector_NP/NonPromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8-evtgen/OniaTree_MC1_Run3_PbPb_track_vector_NP/250829_010938";
  TString prompt_input_file = "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Prompt_with_Track/Oniatree_MC_miniAOD_PbPb23_Prompt_with_Track.root";
  TString nonprompt_input_file = "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Nonprompt_with_Track/Oniatree_MC_miniAOD_PbPb23_Nonprompt_with_Track.root";

  // mc name tag
  std::string mcLabel = "";
  if (isMC)
    mcLabel = isPr ? "PR" : "NP";
  TString safeNameTag = MakeSafeNameTag(nameTag);

  // we create and remove a buffer output file to avoid
  TString bufferFileName = TString::Format(
      "out_buffer_very_long_strange_file_name_never_be_used_MC_%s_pid%d.root",
      safeNameTag.Data(),
      gSystem->GetPid());

  // prepare to connect to input TTree
  TString nameOldTree = "hionia/myTree";
  TChain myChain(nameOldTree);
  TChain trkChain(nameOldTree);

  cout << "Input tag: " << mcLabel << endl;
  cout << "Use single input file: " << useSingleFileInput << endl;
  if (useSingleFileInput)
  {
    // use one merged root file
    TString inputFile = isPr ? prompt_input_file : nonprompt_input_file;
    cout << "Input file: " << inputFile << endl;
    if (gSystem->AccessPathName(inputFile))
    {
      std::cerr << "Input ROOT file is not accessible: " << inputFile << std::endl;
      return;
    }
    myChain.Add(inputFile);
    trkChain.Add(inputFile);
  }
  else
  {
    // check unmerged many root files
    TString base_folder = isPr ? prompt_base_folder : nonprompt_base_folder;
    vector<TString> inputFiles = GetInputRootFiles(base_folder);
    int nInputFiles = inputFiles.size();
    cout << "Base folder: " << base_folder << endl;
    cout << "Number of .root files: " << nInputFiles << endl;
    if (nInputFiles == 0)
    {
      std::cerr << "No .root input files found in " << base_folder << std::endl;
      return;
    }

    for (int i = 0; i < nInputFiles; ++i)
    {
      myChain.Add(inputFiles[i]);
      trkChain.Add(inputFiles[i]);
    }
  }

  // ===== branches ===== //
  // variables in input tree. They will be used for calculatin and will be saved
  Short_t Reco_QQ_size;
  Short_t Reco_QQ_sign[1000]; //[Reco_QQ_size] //
  TClonesArray *Reco_QQ_4mom = nullptr;
  Short_t Reco_QQ_mupl_idx[1000];
  Short_t Reco_QQ_mumi_idx[1000];
  TClonesArray *Reco_mu_4mom = nullptr;
  Short_t Reco_trk_size;
  Int_t nEP = 0;
  Int_t Centrality = -1;
  Int_t centBin = -1;
  Float_t rpAng[12];

  // variables in input tree. They will be used for calculatin but will not saved
  vector<float> *Reco_trk_4mom_pt = nullptr;
  vector<float> *Reco_trk_4mom_eta = nullptr;
  vector<float> *Reco_trk_4mom_phi = nullptr;

  // new variables
  Float_t Q2x_tree, Q2y_tree, ep, ep_rec, ep_flat;
  Float_t Q2x_trackmid_tree, Q2y_trackmid_tree, ep_trackmid, ep_trackmid_rec, ep_trackmid_flat;
  Float_t Q2x_trackm_tree, Q2y_trackm_tree, ep_trackm, ep_trackm_rec, ep_trackm_flat;
  Float_t Q2x_trackp_tree, Q2y_trackp_tree, ep_trackp, ep_trackp_rec, ep_trackp_flat;
  Int_t has_trk_ep, has_trk_mid_ep, has_trk_m_ep, has_trk_p_ep;
  Float_t ep_hfm2, ep_hfp2;
  Float_t cos2_trk_hfp, cos2_trk_hfm, cos2_hfp_hfm;
  const int kHFm2 = 0;
  const int kHFp2 = 1;

  // ===== make output file ===== //
  TFile *outTemp = TFile::Open(bufferFileName.Data(), "RECREATE");
  if (!outTemp || outTemp->IsZombie())
  {
    cout << "[ERROR] cannot open temporary buffer file: " << bufferFileName << "\n";
    return;
  }
  outTemp->cd();

  // ===== set branches from old tree ===== //
  myChain.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size);
  myChain.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign);
  myChain.SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx);
  myChain.SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx);
  myChain.SetBranchAddress("Reco_trk_size", &Reco_trk_size);

  // TClonesArray Branches
  myChain.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom);
  myChain.SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom);
  if (myChain.GetBranch("nEP"))
    myChain.SetBranchAddress("nEP", &nEP);
  if (myChain.GetBranch("Centrality"))
    myChain.SetBranchAddress("Centrality", &Centrality);
  else
    cout << "[WARN] Centrality branch is missing; centrality-binned split will be skipped" << endl;
  if (myChain.GetBranch("rpAng"))
    myChain.SetBranchAddress("rpAng", rpAng);

  // Keep large track vectors out of the cloned output tree; read them only
  // from trkChain for the Q-vector calculation.
  trkChain.SetBranchStatus("*", 0);
  trkChain.SetBranchStatus("Reco_trk_4mom_pt", 1);
  trkChain.SetBranchStatus("Reco_trk_4mom_eta", 1);
  trkChain.SetBranchStatus("Reco_trk_4mom_phi", 1);
  trkChain.SetBranchAddress("Reco_trk_4mom_pt", &Reco_trk_4mom_pt);
  trkChain.SetBranchAddress("Reco_trk_4mom_eta", &Reco_trk_4mom_eta);
  trkChain.SetBranchAddress("Reco_trk_4mom_phi", &Reco_trk_4mom_phi);

  std::set<std::string> dropBranches = {"Reco_trk_4mom_pt", "Reco_trk_4mom_eta", "Reco_trk_4mom_phi"};

  // turn off branches before cloning
  for (const std::string &branch : dropBranches)
    if (myChain.GetBranch(branch.c_str()))
      myChain.SetBranchStatus(branch.c_str(), 0);

  // ===== set new TTree ===== //
  TTree *tempTree = myChain.CloneTree(0);
  tempTree->SetName("tempTree");

  tempTree->Branch("Q2x_tree", &Q2x_tree, "Q2x_tree/F");
  tempTree->Branch("Q2y_tree", &Q2y_tree, "Q2y_tree/F");
  tempTree->Branch("ep", &ep, "ep/F");
  tempTree->Branch("has_trk_ep", &has_trk_ep, "has_trk_ep/I");
  tempTree->Branch("centBin", &centBin, "centBin/I");
  tempTree->Branch("Q2x_trackmid_tree", &Q2x_trackmid_tree, "Q2x_trackmid_tree/F");
  tempTree->Branch("Q2y_trackmid_tree", &Q2y_trackmid_tree, "Q2y_trackmid_tree/F");
  tempTree->Branch("ep_trackmid", &ep_trackmid, "ep_trackmid/F");
  tempTree->Branch("has_trk_mid_ep", &has_trk_mid_ep, "has_trk_mid_ep/I");
  tempTree->Branch("Q2x_trackm_tree", &Q2x_trackm_tree, "Q2x_trackm_tree/F");
  tempTree->Branch("Q2y_trackm_tree", &Q2y_trackm_tree, "Q2y_trackm_tree/F");
  tempTree->Branch("ep_trackm", &ep_trackm, "ep_trackm/F");
  tempTree->Branch("has_trk_m_ep", &has_trk_m_ep, "has_trk_m_ep/I");
  tempTree->Branch("Q2x_trackp_tree", &Q2x_trackp_tree, "Q2x_trackp_tree/F");
  tempTree->Branch("Q2y_trackp_tree", &Q2y_trackp_tree, "Q2y_trackp_tree/F");
  tempTree->Branch("ep_trackp", &ep_trackp, "ep_trackp/F");
  tempTree->Branch("has_trk_p_ep", &has_trk_p_ep, "has_trk_p_ep/I");

  // local TLorentz vectors
  TLorentzVector *qq4mom = nullptr;
  TLorentzVector *Reco_mupl_4mom = nullptr;
  TLorentzVector *Reco_mumi_4mom = nullptr;

  // ===== additional variables ===== //
  TH2F *hQ2xy = new TH2F("hQ2xy", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);
  TH2F *hQ2xyTrackMid = new TH2F("hQ2xyTrackMid", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);
  TH2F *hQ2xyTrackM = new TH2F("hQ2xyTrackM", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);
  TH2F *hQ2xyTrackP = new TH2F("hQ2xyTrackP", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);

  // Dimuon daughter coordinates used to remove autocorrelation tracks.
  float d1Phi[MAXCAN], d1Eta[MAXCAN];
  float d2Phi[MAXCAN], d2Eta[MAXCAN];

  // ===== main loop 1 - raw tracker Q-vector ===== //
  Long64_t nEvent = myChain.GetEntries();
  if (maxEvents > 0 && maxEvents < nEvent)
    nEvent = maxEvents;
  cout << "Total events: " << nEvent << "\n\n";
  const Long64_t reportEvery = GetProgressReportEvery(nEvent);
  TStopwatch loop1Timer;
  loop1Timer.Start();
  cout << "Start loop 1\n\n";

  for (Long64_t ievt = 0; ievt < nEvent; ++ievt)
  {
    PrintLoopProgress("loop 1", ievt, nEvent, reportEvery, loop1Timer);

    myChain.GetEntry(ievt);
    trkChain.GetEntry(ievt);
    centBin = GetCentBin(Centrality);

    Long64_t N = 0;
    Long64_t nsize = Reco_QQ_size;
    // cout << ievt << endl ;
    float Q2x_num = 0;
    float Q2y_num = 0;
    float Qw_pt_num = 0;
    float Q2x_trackmid_num = 0;
    float Q2y_trackmid_num = 0;
    float Qw_trackmid_num = 0;
    float Q2x_trackm_num = 0;
    float Q2y_trackm_num = 0;
    float Qw_trackm_num = 0;
    float Q2x_trackp_num = 0;
    float Q2y_trackp_num = 0;
    float Qw_trackp_num = 0;

    // Collect dimuon daughter tracks for autocorrelation removal.
    for (int iQQ = 0; iQQ < nsize; iQQ++)
    {
      if (Reco_QQ_sign[iQQ] != 0)
        continue; // only opposite-sign muon pairs
      //  if (Reco_QQ_VtxProb[iQQ] < 0.01) continue; // good common vertex proba

      qq4mom = (TLorentzVector *)Reco_QQ_4mom->At(iQQ);

      if (fabs(qq4mom->Eta()) > 2.4)
        continue;
      int iMuPlus = (Reco_QQ_mupl_idx)[iQQ];
      int iMuMinus = (Reco_QQ_mumi_idx)[iQQ];
      Reco_mupl_4mom = (TLorentzVector *)Reco_mu_4mom->At(iMuPlus);
      Reco_mumi_4mom = (TLorentzVector *)Reco_mu_4mom->At(iMuMinus);

      d1Phi[N] = Reco_mupl_4mom->Phi();
      d1Eta[N] = Reco_mupl_4mom->Eta();
      d2Phi[N] = Reco_mumi_4mom->Phi();
      d2Eta[N] = Reco_mumi_4mom->Eta();
      N++;
    }

    // Build the raw tracker Q-vector from selected tracks.
    nsize = Reco_trk_size;
    for (Long64_t it = 0; it < nsize; it++)
    {
      // TLorentzVector *trk4mom = (TLorentzVector *)Reco_trk_4mom->At(it);
      if (!Reco_trk_4mom_eta || !Reco_trk_4mom_pt || !Reco_trk_4mom_phi ||
          it >= (Long64_t)Reco_trk_4mom_eta->size() ||
          it >= (Long64_t)Reco_trk_4mom_pt->size() ||
          it >= (Long64_t)Reco_trk_4mom_phi->size())
        continue;

      float eta = Reco_trk_4mom_eta->at(it);
      float pt = Reco_trk_4mom_pt->at(it);
      float phi = Reco_trk_4mom_phi->at(it);

      bool reject = false;
      for (Long64_t i = 0; i < N; ++i)
      {
        if (fabs(eta - d1Eta[i]) < 1.E-3 && fabs(DeltaPhiPeriodic(phi, d1Phi[i])) < 1.E-3)
          reject = true;
        if (fabs(eta - d2Eta[i]) < 1.E-3 && fabs(DeltaPhiPeriodic(phi, d2Phi[i])) < 1.E-3)
          reject = true;
      }

      if (true == reject)
      {
        // cout << "rejected" << endl ;
        continue;
      }

      // cout << "\n\n ===========: " <<it << endl;
      if (fabs(eta) <= 2.4 && pt > 0.3 && pt < 3)
        AccumulateTrackQ(pt, phi, Q2x_num, Q2y_num, Qw_pt_num);

      // Track subevents follow the rpAng[3/4/5] definitions in epAngle_info.h.
      if (pt > 0.5 && pt < 3)
      {
        if (eta > -0.5 && eta < 0.5)
          AccumulateTrackQ(pt, phi, Q2x_trackmid_num, Q2y_trackmid_num, Qw_trackmid_num);
        if (eta > -2.0 && eta < -1.0)
          AccumulateTrackQ(pt, phi, Q2x_trackm_num, Q2y_trackm_num, Qw_trackm_num);
        if (eta > 1.0 && eta < 2.0)
          AccumulateTrackQ(pt, phi, Q2x_trackp_num, Q2y_trackp_num, Qw_trackp_num);
      }
    }

    // To keep event order; invalid tracker EP entries will be skipped in later calculations
    FinalizeTrackQ(Q2x_num, Q2y_num, Qw_pt_num, Q2x_tree, Q2y_tree, ep, has_trk_ep, hQ2xy);
    FinalizeTrackQ(Q2x_trackmid_num, Q2y_trackmid_num, Qw_trackmid_num, Q2x_trackmid_tree, Q2y_trackmid_tree, ep_trackmid, has_trk_mid_ep, hQ2xyTrackMid);
    FinalizeTrackQ(Q2x_trackm_num, Q2y_trackm_num, Qw_trackm_num, Q2x_trackm_tree, Q2y_trackm_tree, ep_trackm, has_trk_m_ep, hQ2xyTrackM);
    FinalizeTrackQ(Q2x_trackp_num, Q2y_trackp_num, Qw_trackp_num, Q2x_trackp_tree, Q2y_trackp_tree, ep_trackp, has_trk_p_ep, hQ2xyTrackP);
    // cout << Q2x_tree << "  " << Q2y_tree << endl ;

    tempTree->Fill();
  }
  cout << "loop 1 done (raw tracker Q-vector)" << endl;
  PrintStopwatch("loop 1", loop1Timer);

  // ===== main loop 2 - Recentering and flattening inputs ===== //
  // Use each valid track-Q distribution from loop 1 to recenter the matching event-plane angle.
  float mean_Q2x = hQ2xy->GetMean(1);
  float mean_Q2y = hQ2xy->GetMean(2);
  float mean_Q2x_trackmid = hQ2xyTrackMid->GetMean(1);
  float mean_Q2y_trackmid = hQ2xyTrackMid->GetMean(2);
  float mean_Q2x_trackm = hQ2xyTrackM->GetMean(1);
  float mean_Q2y_trackm = hQ2xyTrackM->GetMean(2);
  float mean_Q2x_trackp = hQ2xyTrackP->GetMean(1);
  float mean_Q2y_trackp = hQ2xyTrackP->GetMean(2);
  cout << std::defaultfloat << std::setprecision(6);
  cout << "full tracker mean Q: " << mean_Q2x << "  " << mean_Q2y << endl;
  cout << "trackmid mean Q: " << mean_Q2x_trackmid << "  " << mean_Q2y_trackmid << endl;
  cout << "trackm mean Q: " << mean_Q2x_trackm << "  " << mean_Q2y_trackm << endl;
  cout << "trackp mean Q: " << mean_Q2x_trackp << "  " << mean_Q2y_trackp << endl;

  tempTree->SetBranchAddress("Centrality", &Centrality);
  tempTree->SetBranchAddress("centBin", &centBin);
  tempTree->SetBranchAddress("Q2x_tree", &Q2x_tree);
  tempTree->SetBranchAddress("Q2y_tree", &Q2y_tree);
  tempTree->SetBranchAddress("has_trk_ep", &has_trk_ep);
  tempTree->SetBranchAddress("Q2x_trackmid_tree", &Q2x_trackmid_tree);
  tempTree->SetBranchAddress("Q2y_trackmid_tree", &Q2y_trackmid_tree);
  tempTree->SetBranchAddress("has_trk_mid_ep", &has_trk_mid_ep);
  tempTree->SetBranchAddress("Q2x_trackm_tree", &Q2x_trackm_tree);
  tempTree->SetBranchAddress("Q2y_trackm_tree", &Q2y_trackm_tree);
  tempTree->SetBranchAddress("has_trk_m_ep", &has_trk_m_ep);
  tempTree->SetBranchAddress("Q2x_trackp_tree", &Q2x_trackp_tree);
  tempTree->SetBranchAddress("Q2y_trackp_tree", &Q2y_trackp_tree);
  tempTree->SetBranchAddress("has_trk_p_ep", &has_trk_p_ep);

  TTree *recTree = tempTree->CloneTree(0);
  recTree->SetName("recTree");
  recTree->Branch("ep_rec", &ep_rec, "ep_rec/F");
  recTree->Branch("ep_trackmid_rec", &ep_trackmid_rec, "ep_trackmid_rec/F");
  recTree->Branch("ep_trackm_rec", &ep_trackm_rec, "ep_trackm_rec/F");
  recTree->Branch("ep_trackp_rec", &ep_trackp_rec, "ep_trackp_rec/F");

  const int N = 10;   // Flattening uses harmonics j = 1..10.
  TH1F *means_sin[N]; // Distributions used to extract <sin(2*j*Psi_rec)>.
  TH1F *means_cos[N]; // Distributions used to extract <cos(2*j*Psi_rec)>.
  TH1F *means_sin_trackmid[N];
  TH1F *means_cos_trackmid[N];
  TH1F *means_sin_trackm[N];
  TH1F *means_cos_trackm[N];
  TH1F *means_sin_trackp[N];
  TH1F *means_cos_trackp[N];

  for (int i = 0; i < N; i++)
  {
    means_sin[i] = new TH1F(Form("means_sin_%d", i), Form("<sin(2*%d*Psi_rec)>", i + 1), 20, -1, 1);
    means_cos[i] = new TH1F(Form("means_cos_%d", i), Form("<cos(2*%d*Psi_rec)>", i + 1), 20, -1, 1);
    means_sin_trackmid[i] = new TH1F(Form("means_sin_trackmid_%d", i), Form("<sin(2*%d*Psi_trackmid_rec)>", i + 1), 20, -1, 1);
    means_cos_trackmid[i] = new TH1F(Form("means_cos_trackmid_%d", i), Form("<cos(2*%d*Psi_trackmid_rec)>", i + 1), 20, -1, 1);
    means_sin_trackm[i] = new TH1F(Form("means_sin_trackm_%d", i), Form("<sin(2*%d*Psi_trackm_rec)>", i + 1), 20, -1, 1);
    means_cos_trackm[i] = new TH1F(Form("means_cos_trackm_%d", i), Form("<cos(2*%d*Psi_trackm_rec)>", i + 1), 20, -1, 1);
    means_sin_trackp[i] = new TH1F(Form("means_sin_trackp_%d", i), Form("<sin(2*%d*Psi_trackp_rec)>", i + 1), 20, -1, 1);
    means_cos_trackp[i] = new TH1F(Form("means_cos_trackp_%d", i), Form("<cos(2*%d*Psi_trackp_rec)>", i + 1), 20, -1, 1);
  }

  TStopwatch loop2Timer;
  loop2Timer.Start();
  cout << "\n\nStart loop 2\n\n";
  // Recenter each calculated tracker angle, then accumulate Fourier moments for flattening.
  for (Int_t i = 0; i < nEvent; i++)
  {
    PrintLoopProgress("loop 2", i, nEvent, reportEvery, loop2Timer);
    tempTree->GetEntry(i);
    recTree->GetEntry(i);

    ep_rec = has_trk_ep ? RecenterTrackEp(Q2x_tree, Q2y_tree, mean_Q2x, mean_Q2y) : -99.f;
    ep_trackmid_rec = has_trk_mid_ep ? RecenterTrackEp(Q2x_trackmid_tree, Q2y_trackmid_tree, mean_Q2x_trackmid, mean_Q2y_trackmid) : -99.f;
    ep_trackm_rec = has_trk_m_ep ? RecenterTrackEp(Q2x_trackm_tree, Q2y_trackm_tree, mean_Q2x_trackm, mean_Q2y_trackm) : -99.f;
    ep_trackp_rec = has_trk_p_ep ? RecenterTrackEp(Q2x_trackp_tree, Q2y_trackp_tree, mean_Q2x_trackp, mean_Q2y_trackp) : -99.f;

    if (has_trk_ep)
      FillFlatteningMoments(ep_rec, means_sin, means_cos, N);
    if (has_trk_mid_ep)
      FillFlatteningMoments(ep_trackmid_rec, means_sin_trackmid, means_cos_trackmid, N);
    if (has_trk_m_ep)
      FillFlatteningMoments(ep_trackm_rec, means_sin_trackm, means_cos_trackm, N);
    if (has_trk_p_ep)
      FillFlatteningMoments(ep_trackp_rec, means_sin_trackp, means_cos_trackp, N);

    recTree->Fill();
  }
  cout << "loop 2 done" << endl;
  PrintStopwatch("loop 2", loop2Timer);

  // ===== main loop 3 - Flattening ===== //
  // prepare output - folder, name tag
  TString outBaseDir = "/data/users/pjgwak/work/daily_code_tracker/2026/06/29_ep_resolution_study/pilot_calcTrack_3subEvent/roots";
  TString sampleTag = isPr ? "pr" : "np";
  TString outDir = Form("%s/%s", outBaseDir.Data(), sampleTag.Data());
  if (outDirOverride && std::strlen(outDirOverride) > 0)
    outDir = outDirOverride;
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] failed to create output directory: " << outDir << "\n";
    outTemp->Close();
    std::remove(bufferFileName.Data());
    return;
  }

  string outname = Form(
      "%s/skim_PbPb2023_isMC%d%s_EP_3Sub_%s.root",
      outDir.Data(),
      isMC,
      mcLabel.empty() ? "" : ("_" + mcLabel).c_str(),
      safeNameTag.Data());

  // output file
  TFile *outFile = TFile::Open(outname.c_str(), "RECREATE");
  if (!outFile || outFile->IsZombie() || !outFile->IsWritable())
  {
    cout << "[ERROR] cannot open output file for writing: " << outname << "\n";
    outTemp->Close();
    std::remove(bufferFileName.Data());
    return;
  }
  outFile->cd();
  TDirectory *hioniaDir = outFile->mkdir("hionia");
  outFile->cd("hionia");

  TTree *myTree = recTree->CloneTree(0);
  myTree->SetName("myTree");
  myTree->Branch("ep_flat", &ep_flat, "ep_flat/F");
  myTree->Branch("ep_trackmid_flat", &ep_trackmid_flat, "ep_trackmid_flat/F");
  myTree->Branch("ep_trackm_flat", &ep_trackm_flat, "ep_trackm_flat/F");
  myTree->Branch("ep_trackp_flat", &ep_trackp_flat, "ep_trackp_flat/F");
  myTree->Branch("ep_hfm2", &ep_hfm2, "ep_hfm2/F");
  myTree->Branch("ep_hfp2", &ep_hfp2, "ep_hfp2/F");
  myTree->Branch("cos2_trk_hfp", &cos2_trk_hfp, "cos2_trk_hfp/F");
  myTree->Branch("cos2_trk_hfm", &cos2_trk_hfm, "cos2_trk_hfm/F");
  myTree->Branch("cos2_hfp_hfm", &cos2_hfp_hfm, "cos2_hfp_hfm/F");
  recTree->SetBranchAddress("Centrality", &Centrality);
  recTree->SetBranchAddress("centBin", &centBin);
  recTree->SetBranchAddress("Q2x_tree", &Q2x_tree);
  recTree->SetBranchAddress("Q2y_tree", &Q2y_tree);
  recTree->SetBranchAddress("Q2x_trackmid_tree", &Q2x_trackmid_tree);
  recTree->SetBranchAddress("Q2y_trackmid_tree", &Q2y_trackmid_tree);
  recTree->SetBranchAddress("Q2x_trackm_tree", &Q2x_trackm_tree);
  recTree->SetBranchAddress("Q2y_trackm_tree", &Q2y_trackm_tree);
  recTree->SetBranchAddress("Q2x_trackp_tree", &Q2x_trackp_tree);
  recTree->SetBranchAddress("Q2y_trackp_tree", &Q2y_trackp_tree);
  recTree->SetBranchAddress("ep_rec", &ep_rec);
  recTree->SetBranchAddress("ep_trackmid_rec", &ep_trackmid_rec);
  recTree->SetBranchAddress("ep_trackm_rec", &ep_trackm_rec);
  recTree->SetBranchAddress("ep_trackp_rec", &ep_trackp_rec);
  recTree->SetBranchAddress("has_trk_ep", &has_trk_ep);
  recTree->SetBranchAddress("has_trk_mid_ep", &has_trk_mid_ep);
  recTree->SetBranchAddress("has_trk_m_ep", &has_trk_m_ep);
  recTree->SetBranchAddress("has_trk_p_ep", &has_trk_p_ep);
  if (recTree->GetBranch("nEP"))
    recTree->SetBranchAddress("nEP", &nEP);
  if (recTree->GetBranch("rpAng"))
    recTree->SetBranchAddress("rpAng", rpAng);

  // Flattening constants extracted from loop 2 Fourier-moment histograms.
  double mean_sin[N];
  double mean_cos[N];
  double mean_sin_trackmid[N];
  double mean_cos_trackmid[N];
  double mean_sin_trackm[N];
  double mean_cos_trackm[N];
  double mean_sin_trackp[N];
  double mean_cos_trackp[N];

  TH1F *hCos2TrkHFp = new TH1F("hCos2TrkHFp", ";cos 2(#Psi_{trk}-#Psi_{HF+});Entries", 200, -1, 1);
  TH1F *hCos2TrkHFm = new TH1F("hCos2TrkHFm", ";cos 2(#Psi_{trk}-#Psi_{HF-});Entries", 200, -1, 1);
  TH1F *hCos2HFpHFm = new TH1F("hCos2HFpHFm", ";cos 2(#Psi_{HF+}-#Psi_{HF-});Entries", 200, -1, 1);
  TH1F *hResolution3Sub = new TH1F("hResolution3Sub", "3-subevent resolution summary", 7, 0.5, 7.5);
  TH1F *hCos2TrkHFpCent[kNCentBins];
  TH1F *hCos2TrkHFmCent[kNCentBins];
  TH1F *hCos2HFpHFmCent[kNCentBins];
  for (int ic = 0; ic < kNCentBins; ++ic)
  {
    hCos2TrkHFpCent[ic] = new TH1F(Form("hCos2TrkHFp_cent%d", ic), Form("%s;cos 2(#Psi_{trk}-#Psi_{HF+});Entries", GetCentBinLabel(ic).Data()), 200, -1, 1);
    hCos2TrkHFmCent[ic] = new TH1F(Form("hCos2TrkHFm_cent%d", ic), Form("%s;cos 2(#Psi_{trk}-#Psi_{HF-});Entries", GetCentBinLabel(ic).Data()), 200, -1, 1);
    hCos2HFpHFmCent[ic] = new TH1F(Form("hCos2HFpHFm_cent%d", ic), Form("%s;cos 2(#Psi_{HF+}-#Psi_{HF-});Entries", GetCentBinLabel(ic).Data()), 200, -1, 1);
  }
  TH1F *hResolution3SubCent = new TH1F("hResolution3SubCent", ";Centrality;R_{trk}", kNCentBins, 0.5, kNCentBins + 0.5);
  TH1F *hCos2TrkHFpCentMean = new TH1F("hCos2TrkHFpCentMean", ";Centrality;<cos2 Trk-HF+>", kNCentBins, 0.5, kNCentBins + 0.5);
  TH1F *hCos2TrkHFmCentMean = new TH1F("hCos2TrkHFmCentMean", ";Centrality;<cos2 Trk-HF->", kNCentBins, 0.5, kNCentBins + 0.5);
  TH1F *hCos2HFpHFmCentMean = new TH1F("hCos2HFpHFmCentMean", ";Centrality;<cos2 HF+-HF->", kNCentBins, 0.5, kNCentBins + 0.5);
  TH1F *hCos2CentEntries = new TH1F("hCos2CentEntries", ";Centrality;Entries", kNCentBins, 0.5, kNCentBins + 0.5);
  SetupCentralityHist(hResolution3SubCent);
  SetupCentralityHist(hCos2TrkHFpCentMean);
  SetupCentralityHist(hCos2TrkHFmCentMean);
  SetupCentralityHist(hCos2HFpHFmCentMean);
  SetupCentralityHist(hCos2CentEntries);

  TStopwatch loop3Timer;
  loop3Timer.Start();
  cout << "\n\nStart loop 3\n\n";
  ReadFlatteningMeans(means_sin, means_cos, mean_sin, mean_cos, N);
  ReadFlatteningMeans(means_sin_trackmid, means_cos_trackmid, mean_sin_trackmid, mean_cos_trackmid, N);
  ReadFlatteningMeans(means_sin_trackm, means_cos_trackm, mean_sin_trackm, mean_cos_trackm, N);
  ReadFlatteningMeans(means_sin_trackp, means_cos_trackp, mean_sin_trackp, mean_cos_trackp, N);
  // Apply tracker flattening, read HF subevent angles, and fill 3-subevent correlations.
  for (Int_t i = 0; i < nEvent; i++)
  {
    PrintLoopProgress("loop 3", i, nEvent, reportEvery, loop3Timer);

    recTree->GetEntry(i);
    // cout << i << endl ;
    ep_flat = -99.f;
    ep_trackmid_flat = -99.f;
    ep_trackm_flat = -99.f;
    ep_trackp_flat = -99.f;
    ep_hfm2 = -99.f;
    ep_hfp2 = -99.f;
    cos2_trk_hfp = -99.f;
    cos2_trk_hfm = -99.f;
    cos2_hfp_hfm = -99.f;

    if (has_trk_ep)
      ep_flat = FlattenTrackEp(ep_rec, mean_sin, mean_cos, N);
    if (has_trk_mid_ep)
      ep_trackmid_flat = FlattenTrackEp(ep_trackmid_rec, mean_sin_trackmid, mean_cos_trackmid, N);
    if (has_trk_m_ep)
      ep_trackm_flat = FlattenTrackEp(ep_trackm_rec, mean_sin_trackm, mean_cos_trackm, N);
    if (has_trk_p_ep)
      ep_trackp_flat = FlattenTrackEp(ep_trackp_rec, mean_sin_trackp, mean_cos_trackp, N);

    // HF angles come from the input EP array; require indices 0 and 1 to exist.
    if (nEP > kHFp2)
    {
      ep_hfm2 = rpAng[kHFm2];
      ep_hfp2 = rpAng[kHFp2];
    }

    // Fill observed 3-subevent correlations only when all three angles are valid.
    if (has_trk_ep && ep_hfp2 > -90.f && ep_hfm2 > -90.f)
    {
      cos2_trk_hfp = cos(2.0 * (ep_flat - ep_hfp2));
      cos2_trk_hfm = cos(2.0 * (ep_flat - ep_hfm2));
      cos2_hfp_hfm = cos(2.0 * (ep_hfp2 - ep_hfm2));
      hCos2TrkHFp->Fill(cos2_trk_hfp);
      hCos2TrkHFm->Fill(cos2_trk_hfm);
      hCos2HFpHFm->Fill(cos2_hfp_hfm);

      for (int ic = 0; ic < kNCentBins; ++ic)
      {
        if (!IsInCentBin(Centrality, ic))
          continue;
        hCos2TrkHFpCent[ic]->Fill(cos2_trk_hfp);
        hCos2TrkHFmCent[ic]->Fill(cos2_trk_hfm);
        hCos2HFpHFmCent[ic]->Fill(cos2_hfp_hfm);
      }
    }
    myTree->Fill();
  }
  cout << "loop 3 done" << endl;
  PrintStopwatch("loop 3", loop3Timer);

  // Convert correlation distributions into the tracker EP resolution summary.
  Double_t Cos2TrkHFp_Mean = hCos2TrkHFp->GetMean(1);
  Double_t Cos2TrkHFm_Mean = hCos2TrkHFm->GetMean(1);
  Double_t Cos2HFpHFm_Mean = hCos2HFpHFm->GetMean(1);
  Double_t Cos2TrkR_Mean = SafeResolution(Cos2TrkHFp_Mean * Cos2TrkHFm_Mean, Cos2HFpHFm_Mean);
  cout << std::setprecision(12)
       << "Cos2TrkR_Mean = sqrt(" << Cos2TrkHFp_Mean << " * " << Cos2TrkHFm_Mean << " / " << Cos2HFpHFm_Mean << ") = " << Cos2TrkR_Mean << endl;

  hResolution3Sub->GetXaxis()->SetBinLabel(1, "<TrkHFp>");
  hResolution3Sub->GetXaxis()->SetBinLabel(2, "<TrkHFm>");
  hResolution3Sub->GetXaxis()->SetBinLabel(3, "<HFpHFm>");
  hResolution3Sub->GetXaxis()->SetBinLabel(4, "R_Trk");
  hResolution3Sub->GetXaxis()->SetBinLabel(5, "entries TrkHFp");
  hResolution3Sub->GetXaxis()->SetBinLabel(6, "entries TrkHFm");
  hResolution3Sub->GetXaxis()->SetBinLabel(7, "entries HFpHFm");
  hResolution3Sub->SetBinContent(1, Cos2TrkHFp_Mean);
  hResolution3Sub->SetBinContent(2, Cos2TrkHFm_Mean);
  hResolution3Sub->SetBinContent(3, Cos2HFpHFm_Mean);
  hResolution3Sub->SetBinContent(4, Cos2TrkR_Mean);
  hResolution3Sub->SetBinContent(5, hCos2TrkHFp->GetEntries());
  hResolution3Sub->SetBinContent(6, hCos2TrkHFm->GetEntries());
  hResolution3Sub->SetBinContent(7, hCos2HFpHFm->GetEntries());

  for (int ic = 0; ic < kNCentBins; ++ic)
  {
    const Double_t meanTrkHFpCent = hCos2TrkHFpCent[ic]->GetMean(1);
    const Double_t meanTrkHFmCent = hCos2TrkHFmCent[ic]->GetMean(1);
    const Double_t meanHFpHFmCent = hCos2HFpHFmCent[ic]->GetMean(1);
    const Double_t rTrkCent = SafeResolution(meanTrkHFpCent * meanTrkHFmCent, meanHFpHFmCent);
    hResolution3SubCent->SetBinContent(ic + 1, rTrkCent);
    hCos2TrkHFpCentMean->SetBinContent(ic + 1, meanTrkHFpCent);
    hCos2TrkHFmCentMean->SetBinContent(ic + 1, meanTrkHFmCent);
    hCos2HFpHFmCentMean->SetBinContent(ic + 1, meanHFpHFmCent);
    hCos2CentEntries->SetBinContent(ic + 1, hCos2TrkHFpCent[ic]->GetEntries());
    cout << "R_Trk " << GetCentBinLabel(ic) << " = sqrt(" << meanTrkHFpCent << " * " << meanTrkHFmCent << " / " << meanHFpHFmCent << ") = " << rTrkCent << endl;
  }

  // ===== saving and cleaning ===== //
  myTree->Write();
  outFile->cd();
  hQ2xy->Write();
  hQ2xyTrackMid->Write();
  hQ2xyTrackM->Write();
  hQ2xyTrackP->Write();
  hCos2TrkHFp->Write();
  hCos2TrkHFm->Write();
  hCos2HFpHFm->Write();
  hResolution3Sub->Write();
  for (int ic = 0; ic < kNCentBins; ++ic)
  {
    hCos2TrkHFpCent[ic]->Write();
    hCos2TrkHFmCent[ic]->Write();
    hCos2HFpHFmCent[ic]->Write();
  }
  hResolution3SubCent->Write();
  hCos2TrkHFpCentMean->Write();
  hCos2TrkHFmCentMean->Write();
  hCos2HFpHFmCentMean->Write();
  hCos2CentEntries->Write();

  outTemp->Close();
  outFile->Close();

  // remove buffer output file
  std::remove(bufferFileName.Data());

  PrintStopwatch("total", totalTimer);
}