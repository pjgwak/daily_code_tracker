/* check
--- Current cut ---
opposite-sign dimuon:
Reco_QQ_sign[iQQ] != 0 continue
dimuon |eta| < 2.4
daughter track 제거:
|deta| < 0.03, |dphi| < 0.03
tracker Q-vector track:
|eta| <= 2.4, 0.3 < pt < 3


--- Exist in the referecne but not this code ---
trigger cut:csTree.trigHLT()[10]

event selection:csTree.evtSel()[2]
csTree.evtSel()[3]

centrality:100 <= centrality < 160

J/psi candidate condition:softCand
|y| < 2.4
2.5 < mass < 4.3
0.2 <= pt < 20

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

void getEpAngles_calcTrack_3subEvent(bool isMC = true, bool isPr = true, Long64_t maxEvents = 300000, const char *outDirOverride = "", const char *nameTag = "", bool useSingleFileInput = true)
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
  Float_t rpAng[12];

  // variables in input tree. They will be used for calculatin but will not saved
  vector<float> *Reco_trk_4mom_pt = nullptr;
  vector<float> *Reco_trk_4mom_eta = nullptr;
  vector<float> *Reco_trk_4mom_phi = nullptr;

  // new variables
  Float_t Q2x_tree, Q2y_tree, ep, ep_rec, ep_flat;
  Int_t has_trk_ep;
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

  // local TLorentz vectors
  TLorentzVector *qq4mom = nullptr;
  TLorentzVector *Reco_mupl_4mom = nullptr;
  TLorentzVector *Reco_mumi_4mom = nullptr;

  // ===== additional variables ===== //
  TH2F *hQ2xy = new TH2F("hQ2xy", " ", 80, -4.0, 4.0, 80, -4.0, 4.0);

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

    Long64_t N = 0;
    Long64_t nsize = Reco_QQ_size;
    // cout << ievt << endl ;
    float Q2x_num = 0;
    float Q2y_num = 0;
    float Qw_pt_num = 0;

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
      {
        Q2x_num += pt * cos(2 * phi);
        Q2y_num += pt * sin(2 * phi);
        Qw_pt_num += pt;
      }
    }

    if (Qw_pt_num > 0)
    {
      has_trk_ep = 1;
      Q2x_tree = Q2x_num / Qw_pt_num;
      Q2y_tree = Q2y_num / Qw_pt_num;
      ep = 0.5 * atan2(Q2y_tree, Q2x_tree);
      hQ2xy->Fill(Q2x_tree, Q2y_tree);
    }
    else
    {
      // Keep event order; invalid tracker EP entries are skipped in later calculations.
      has_trk_ep = 0;
      Q2x_tree = 0;
      Q2y_tree = 0;
      ep = -99.f;
    }

    // cout << Q2x_tree << "  " << Q2y_tree << endl ;

    tempTree->Fill();
  }
  cout << "loop 1 done (raw tracker Q-vector)" << endl;
  PrintStopwatch("loop 1", loop1Timer);

  // ===== main loop 2 - Recentering and flattening inputs ===== //
  // Use the valid track-Q distribution from loop 1 to recenter each event.
  float mean_Q2x = hQ2xy->GetMean(1);
  float mean_Q2y = hQ2xy->GetMean(2);
  cout << mean_Q2x << "  " << mean_Q2y << endl;

  tempTree->SetBranchAddress("Q2x_tree", &Q2x_tree);
  tempTree->SetBranchAddress("Q2y_tree", &Q2y_tree);
  tempTree->SetBranchAddress("has_trk_ep", &has_trk_ep);

  TTree *recTree = tempTree->CloneTree(0);
  recTree->SetName("recTree");
  recTree->Branch("ep_rec", &ep_rec, "ep_rec/F");

  const int N = 10;   // Flattening uses harmonics j = 1..10.
  TH1F *means_sin[N]; // Distributions used to extract <sin(2*j*Psi_rec)>.
  TH1F *means_cos[N]; // Distributions used to extract <cos(2*j*Psi_rec)>.

  for (int i = 0; i < N; i++)
  {
    means_sin[i] = new TH1F(Form("means_sin_%d", i), Form("<sin(2*%d*Psi_rec)>", i + 1), 20, -1, 1);
    means_cos[i] = new TH1F(Form("means_cos_%d", i), Form("<cos(2*%d*Psi_rec)>", i + 1), 20, -1, 1);
  }

  TStopwatch loop2Timer;
  loop2Timer.Start();
  cout << "\n\nStart loop 2\n\n";
  // Apply tracker flattening, read HF subevent angles, and fill 3-subevent correlations.
  for (Int_t i = 0; i < nEvent; i++)
  {
    PrintLoopProgress("loop 2", i, nEvent, reportEvery, loop2Timer);
    tempTree->GetEntry(i);
    recTree->GetEntry(i);
    if (has_trk_ep)
    {
      float ep1 = 0.5 * atan2(Q2y_tree - mean_Q2y, Q2x_tree - mean_Q2x);
      ep_rec = ep1;
      float j = 1.;

      // Store j = 1..10 Fourier moments; their means become flattening constants.
      for (int i = 0; i < N; i++)
      {
        means_sin[i]->Fill(sin(2 * j * ep_rec));
        means_cos[i]->Fill(cos(2 * j * ep_rec));
        j += 1;
      }
    }
    else
    {
      ep_rec = -99.f;
    }
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
  myTree->Branch("ep_hfm2", &ep_hfm2, "ep_hfm2/F");
  myTree->Branch("ep_hfp2", &ep_hfp2, "ep_hfp2/F");
  myTree->Branch("cos2_trk_hfp", &cos2_trk_hfp, "cos2_trk_hfp/F");
  myTree->Branch("cos2_trk_hfm", &cos2_trk_hfm, "cos2_trk_hfm/F");
  myTree->Branch("cos2_hfp_hfm", &cos2_hfp_hfm, "cos2_hfp_hfm/F");
  recTree->SetBranchAddress("ep_rec", &ep_rec);
  recTree->SetBranchAddress("has_trk_ep", &has_trk_ep);
  if (recTree->GetBranch("nEP"))
    recTree->SetBranchAddress("nEP", &nEP);
  if (recTree->GetBranch("rpAng"))
    recTree->SetBranchAddress("rpAng", rpAng);

  // Flattening constants extracted from loop 2 Fourier-moment histograms.
  double mean_sin[N];
  double mean_cos[N];

  TH1F *hCos2TrkHFp = new TH1F("hCos2TrkHFp", ";cos 2(#Psi_{trk}-#Psi_{HF+});Entries", 200, -1, 1);
  TH1F *hCos2TrkHFm = new TH1F("hCos2TrkHFm", ";cos 2(#Psi_{trk}-#Psi_{HF-});Entries", 200, -1, 1);
  TH1F *hCos2HFpHFm = new TH1F("hCos2HFpHFm", ";cos 2(#Psi_{HF+}-#Psi_{HF-});Entries", 200, -1, 1);
  TH1F *hResolution3Sub = new TH1F("hResolution3Sub", "3-subevent resolution summary", 7, 0.5, 7.5);

  TStopwatch loop3Timer;
  loop3Timer.Start();
  cout << "\n\nStart loop 3\n\n";
  for (int i = 0; i < N; i++)
  {
    mean_sin[i] = means_sin[i]->GetMean();
    mean_cos[i] = means_cos[i]->GetMean();
  }

  // Apply tracker flattening, read HF subevent angles, and fill 3-subevent correlations.
  for (Int_t i = 0; i < nEvent; i++)
  {
    PrintLoopProgress("loop 3", i, nEvent, reportEvery, loop3Timer);

    recTree->GetEntry(i);
    // cout << i << endl ;
    ep_flat = -99.f;
    ep_hfm2 = -99.f;
    ep_hfp2 = -99.f;
    cos2_trk_hfp = -99.f;
    cos2_trk_hfm = -99.f;
    cos2_hfp_hfm = -99.f;

    if (has_trk_ep)
    {
      // Flatten calculated tracker EP using j = 1..10 Fourier terms.
      double deltaPsi2 = 0.0;
      for (int n = 0; n < N; n++)
      {
        const int j = n + 1;
        double term = (2.0 / j) * (-mean_sin[n] * cos(2.0 * j * ep_rec) + mean_cos[n] * sin(2.0 * j * ep_rec));
        deltaPsi2 += term;
      }
      ep_flat = ep_rec + 0.5 * deltaPsi2;
    }

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

  // ===== saving and cleaning ===== //
  myTree->Write();
  outFile->cd();
  hQ2xy->Write();
  hCos2TrkHFp->Write();
  hCos2TrkHFm->Write();
  hCos2HFpHFm->Write();
  hResolution3Sub->Write();

  outTemp->Close();
  outFile->Close();

  // remove buffer output file
  std::remove(bufferFileName.Data());

  PrintStopwatch("total", totalTimer);
}
