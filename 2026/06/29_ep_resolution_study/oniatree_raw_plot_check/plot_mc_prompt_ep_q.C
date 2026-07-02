// headers
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <vector>

// std namespace call
using std::cout;
using std::endl;
using std::vector;

// custom namespace to save the setting
namespace
{
  const char *kTreeName = "hionia/myTree";
  const int kMaxEpIndex = 5;
  const int kNPlotIndex = kMaxEpIndex + 1;

  const char *kEpName[kNPlotIndex] = {
      "HFm2", "HFp2", "HF2", "trackmid2", "trackm2", "trackp2"};

  const char *kBranchNames[] = {
      "rpAng_origin", "rpSin_origin", "rpCos_origin",
      "rpAng", "rpSin", "rpCos"};
  const int kNBranch = sizeof(kBranchNames) / sizeof(kBranchNames[0]);

  // input path
  const char *kSinglePromptInput =
      "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Prompt_with_Track/"
      "Oniatree_MC_miniAOD_PbPb23_Prompt_with_Track.root";

  const char *kPromptBaseFolder =
      "/data/Oniatree/polarization/oniatree_5p36/PbPb23_MC_Prompt_with_Track/"
      "OniaTree_MC1_Run3_PbPb_track_vector_PR/"
      "PromptJPsiToMuMu_Pthat2_TuneCP5_HydjetDrumMB_5p36TeV_pythia8/"
      "OniaTree_MC1_Run3_PbPb_track_vector_PR/250828_074840";

  const char *kRootlogon = "/data/users/pjgwak/input_files/rootlogon.C";

  bool HasRootSuffix(const char *fileName)
  {
    const size_t len = std::strlen(fileName);
    return len > 5 && std::strcmp(fileName + len - 5, ".root") == 0;
  }

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

  vector<TString> GetInputRootFiles(const TString &baseFolder)
  {
    vector<TString> inputFiles;
    DIR *dir = opendir(baseFolder.Data());
    if (!dir)
      return inputFiles;

    AddRootFilesInDir(baseFolder, inputFiles);

    dirent *entry = nullptr;
    while ((entry = readdir(dir)) != nullptr)
    {
      const char *name = entry->d_name;
      if (std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
        continue;
      if (HasRootSuffix(name))
        continue;

      AddRootFilesInDir(TString::Format("%s/%s", baseFolder.Data(), name), inputFiles);
    }
    closedir(dir);

    std::sort(inputFiles.begin(), inputFiles.end());
    return inputFiles;
  }
  
  // check existance of branch
  bool HasBranch(TTree *tree, const char *name)
  {
    return tree && tree->GetBranch(name);
  }

  TString IndexSelection(TTree *tree, int index)
  {
    if (HasBranch(tree, "nEP"))
      return TString::Format("nEP>%d", index);
    return "";
  }

  // simple histogram drawing commands with TTree
  Long64_t DrawInto(TChain &chain, const TString &expr, TH1 *hist, const TString &selection, Long64_t maxEvents)
  {
    const TString drawExpr = TString::Format("%s>>%s", expr.Data(), hist->GetName());
    const Long64_t nToRead = (maxEvents > 0) ? maxEvents : TTree::kMaxEntries;
    return chain.Draw(drawExpr.Data(), selection.Data(), "goff", nToRead);
  }

  void SaveCanvas(TCanvas &canvas, const TString &outDir, const TString &name)
  {
    // canvas.SaveAs(TString::Format("%s/%s.png", outDir.Data(), name.Data()));
    canvas.SaveAs(TString::Format("%s/%s.pdf", outDir.Data(), name.Data()));
  }

  void SetupStyle()
  {
    gROOT->SetBatch(kTRUE);
    if (!gSystem->AccessPathName(kRootlogon))
      gROOT->Macro(kRootlogon);

    gStyle->SetEndErrorSize(0);
    gStyle->SetOptStat(1110);
    gStyle->SetTitleFontSize(0.04);
    TGaxis::SetMaxDigits(3);
  }

  void StylePad()
  {
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.13);
    gPad->SetTicks(1, 1);
  }

  void StyleHist(TH1 *hist, int color)
  {
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
    hist->SetStats(true);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetXaxis()->SetLabelSize(0.040);
    hist->GetYaxis()->SetLabelSize(0.040);
    if (hist->GetMaximum() > 0)
      hist->SetMaximum(hist->GetMaximum() * 1.6);
  }

  void GetHistBinning(const TString &branchName, int epIndex, int &nBins, double &xmin, double &xmax)
  {
    nBins = 100;
    if (branchName.Contains("Ang"))
    {
      xmin = -1.8;
      xmax = 1.8;
      return;
    }

    // rpSin/rpCos are Q-vector components, not unit sin/cos values.
    if (epIndex <= 2)
    {
      xmin = -500.0;
      xmax = 500.0;
    }
    else
    {
      xmin = -150.0;
      xmax = 150.0;
    }
  }

  int BranchColor(const TString &branchName)
  {
    if (branchName.Contains("Ang"))
      return kBlack;
    if (branchName.Contains("Sin"))
      return kRed + 1;
    return kBlue + 1;
  }
} // namespace

void plot_mc_prompt_ep_q(Long64_t maxEvents = 200000,
                         bool useSingleFileInput = true,
                         const char *inputOverride = "")
{
  SetupStyle();

  const char *outDir = "figs";
  const char *rootDir = "roots";

  // bring inputs
  TChain chain(kTreeName);
  if (inputOverride && std::strlen(inputOverride) > 0)
  {
    chain.Add(inputOverride);
    cout << "Input override: " << inputOverride << endl;
  }
  else if (useSingleFileInput)
  {
    chain.Add(kSinglePromptInput);
    cout << "Input file: " << kSinglePromptInput << endl;
  }
  else
  {
    const vector<TString> inputFiles = GetInputRootFiles(kPromptBaseFolder);
    cout << "Input base folder: " << kPromptBaseFolder << endl;
    cout << "Number of .root files: " << inputFiles.size() << endl;
    for (const TString &inputFile : inputFiles)
      chain.Add(inputFile);
  }

  // check: is input empty?
  if (chain.GetEntries() <= 0)
  {
    cout << "[ERROR] No entries found in " << kTreeName << endl;
    return;
  }

  // connect trees
  chain.LoadTree(0);
  TTree *tree = chain.GetTree();
  for (int ib = 0; ib < kNBranch; ++ib)
  {
    if (!HasBranch(tree, kBranchNames[ib]))
    {
      cout << "[ERROR] Missing branch: " << kBranchNames[ib] << endl;
      return;
    }
  }

  // prepare output folders
  if (gSystem->AccessPathName(outDir) && gSystem->mkdir(outDir, true) != 0)
  {
    cout << "[ERROR] Failed to create output directory: " << outDir << endl;
    return;
  }

  if (gSystem->AccessPathName(rootDir) && gSystem->mkdir(rootDir, true) != 0)
  {
    cout << "[ERROR] Failed to create output directory: " << rootDir << endl;
    return;
  }

  const Long64_t nEntries = chain.GetEntries();
  const Long64_t nToRead = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;
  cout << "Entries in chain: " << nEntries << endl;
  cout << "Entries to draw: " << nToRead << endl;

  // prepare output file
  TFile outFile(TString::Format("%s/mc_prompt_rp_array_hists.root", rootDir), "RECREATE");
  if (outFile.IsZombie())
  {
    cout << "[ERROR] Failed to create ROOT output file in " << rootDir << endl;
    return;
  }

  // main loop
  for (int ib = 0; ib < kNBranch; ++ib)
  {
    const TString branchName = kBranchNames[ib];
    const TString branchOutDir = TString::Format("%s/%s", outDir, branchName.Data());

    if (gSystem->AccessPathName(branchOutDir) && gSystem->mkdir(branchOutDir, true) != 0)
    {
      cout << "[ERROR] Failed to create output directory: " << branchOutDir << endl;
      return;
    }

    TH1D *hBranch[kNPlotIndex] = {nullptr};
    for (int i = 0; i < kNPlotIndex; ++i)
    {
      const TString expr = TString::Format("%s[%d]", branchName.Data(), i);
      const TString selection = IndexSelection(tree, i);
      const TString tag = TString::Format("%s_idx%d_%s", branchName.Data(), i, kEpName[i]);

      int nBins = 100;
      double xmin = 0.0;
      double xmax = 0.0;
      GetHistBinning(branchName, i, nBins, xmin, xmax);

      hBranch[i] = new TH1D(TString::Format("h_%s", tag.Data()),
                            TString::Format("MC prompt %s[%d] (%s);%s[%d];Events",
                                            branchName.Data(), i, kEpName[i],
                                            branchName.Data(), i),
                            nBins, xmin, xmax);
      DrawInto(chain, expr, hBranch[i], selection, nToRead);
      StyleHist(hBranch[i], BranchColor(branchName));

      TCanvas canvas(TString::Format("c_%s", tag.Data()), "", 800, 700);
      canvas.cd();
      StylePad();
      hBranch[i]->Draw("hist");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextFont(42);
      latex.SetTextSize(0.040);
      latex.DrawLatex(0.62, 0.84, kEpName[i]);
      SaveCanvas(canvas, branchOutDir, tag);

      outFile.cd();
      hBranch[i]->Write();
    }

    TCanvas cSummary(TString::Format("c_summary_%s", branchName.Data()), "", 1800, 1200);
    cSummary.Divide(3, 2);
    for (int i = 0; i < kNPlotIndex; ++i)
    {
      cSummary.cd(i + 1);
      StylePad();
      hBranch[i]->Draw("hist");
      TLatex latex;
      latex.SetNDC();
      latex.SetTextFont(42);
      latex.SetTextSize(0.055);
      latex.DrawLatex(0.58, 0.84, kEpName[i]);
    }
    SaveCanvas(cSummary, branchOutDir, TString::Format("%s_summary_idx0to5", branchName.Data()));
  }

  outFile.Close();
}