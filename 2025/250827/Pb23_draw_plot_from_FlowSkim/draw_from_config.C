#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TDirectory.h"
#include "ConfigPlots.h"

// preview setting -> Use small portion of events
const Long64_t nent = PREVIEW_ON ? PREVIEW_NENTRIES : TTree::kMaxEntries;
const Long64_t first = PREVIEW_ON ? PREVIEW_FIRSTENTRY : 0;

// ===== helper functions =====
// prepare output root file
TFile *gOutFile = nullptr;
TDirectory *gDirH1 = nullptr, *gDirH2 = nullptr, *gDirOv = nullptr;
std::ofstream gManifest;

void OpenOutputs()
{
  gOutFile = new TFile((OUT_DIR + "/hists.root").c_str(), "RECREATE");
  gDirH1 = gOutFile->mkdir("h1");
  gDirH2 = gOutFile->mkdir("h2");
  gDirOv = gOutFile->mkdir("overlay");
  gManifest.open((OUT_DIR + "/manifest.tsv").c_str());
  gManifest << "type\tname\texpr\tcut\n";
}

void CloseOutputs()
{
  if (gOutFile)
  {
    gOutFile->Write();
    gOutFile->Close();
    delete gOutFile;
    gOutFile = nullptr;
  }
  if (gManifest.is_open())
    gManifest.close();
}

// --- save histograms ---
void Save1D(TH1 *h, const std::string &expr, const TCut &cut)
{
  if (!h)
    return;
  gDirH1->cd();
  h->Write("", TObject::kOverwrite);
  if (gManifest.is_open())
    gManifest << "1D\t" << h->GetName() << "\t" << expr << "\t" << cut.GetTitle() << "\n";
}

void Save2D(TH2 *h2, const std::string &expr, const TCut &cut)
{
  if (!h2)
    return;
  gDirH2->cd();
  h2->Write("", TObject::kOverwrite);
  if (gManifest.is_open())
    gManifest << "2D\t" << h2->GetName() << "\t" << expr << "\t" << cut.GetTitle() << "\n";
}

// ===== main macro =====
void draw_from_config()
{
  gStyle->SetOptStat(1110);
  gSystem->mkdir(OUT_DIR.c_str(), true);

  OpenOutputs();

  // ----- Inputs -----
  TFile *f = TFile::Open(INPUT_FILE.c_str());
  if (!f || f->IsZombie())
  {
    std::cerr << "Can't open the input root file\n";
    return;
  }
  TTree *tree = (TTree *)f->Get(TREE_NAME.c_str());
  if (!tree)
  {
    std::cerr << "No tree: " << TREE_NAME << "\n";
    f->Close();
    delete f;
    return;
  }

  // ----- canvas -----
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  TCanvas *c2 = new TCanvas("c2", "", 800, 800);

  // -------- 1D plots --------
  for (const auto &p : CFG_1D)
  {
    const auto slices = SlicesForTag(p.sliceTag);
    for (const auto &s : slices)
    {
      const auto name = p.name + s.suffix;
      const auto cut = p.baseCut && s.cut && CUT_JPSI_MASS;

      TString hdef = Form("%s(%d,%g,%g)", name.c_str(), p.nbins, p.xmin, p.xmax);
      tree->Draw(Form("%s >> %s", p.expr.c_str(), hdef.Data()), cut, "goff", nent, first);
      TH1 *h = (TH1 *)gDirectory->Get(name.c_str());
      if (!h)
      {
        std::cerr << "Can't make a histogram from the branch: " << name << "\n";
        continue;
      }
      h->SetTitle(p.title.c_str());
      c1->cd();
      h->Draw("E");
      c1->SaveAs(Form("%s/%s.png", OUT_DIR.c_str(), name.c_str()));
      Save1D(h, p.expr, cut);
    }
  }

  // -------- 2D plots --------
  for (const auto &p : CFG_2D)
  {
    const auto slices = SlicesForTag(p.sliceTag);
    for (const auto &s : slices)
    {
      const auto name = p.name + s.suffix;
      const auto cut = p.baseCut && s.cut && CUT_JPSI_MASS;

      TString hdef = Form("%s(%d,%g,%g,%d,%g,%g)",
                          name.c_str(), p.nxbins, p.xmin, p.xmax,
                          p.nybins, p.ymin, p.ymax);
      tree->Draw(Form("%s >> %s", p.exprYX.c_str(), hdef.Data()), cut, "goff", nent, first);
      TH2 *h2 = (TH2 *)gDirectory->Get(name.c_str());
      if (!h2)
      {
        std::cerr << "Can't make a 2D histogram from the branchs " << name << "\n";
        continue;
      }

      h2->SetTitle(p.title.c_str());
      c2->cd();
      h2->Draw("COLZ");
      c2->SaveAs(Form("%s/%s_colz.png", OUT_DIR.c_str(), name.c_str()));
      Save2D(h2, p.exprYX, cut);
    }
  }  
  CloseOutputs();

  delete c1;
  delete c2;
  f->Close();
  delete f;
}
