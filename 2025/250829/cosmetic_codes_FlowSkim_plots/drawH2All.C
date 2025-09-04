#include "PlotCore.C"
#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <iostream>

// TH2D container
struct H2Spec
{
  TString subdir;
  TString xtitle;
  TString ytitle;
};

// decide output sub directories, x and y titles
static H2Spec deduceH2Spec(const TString &n)
{
  // branch types - rapidities
  H2Spec s;
  if (n.Contains("pt_vs_y"))
  {
    s.subdir = "h2_pt_vs_y";
    s.xtitle = "y";
    s.ytitle = "p_{T} (GeV/c)";
  }
  else if (n.Contains("pt_vs_mass"))
  {
    s.subdir = "h2_pt_vs_mass";
    s.xtitle = "M_{#mu#mu} (GeV/c^{2})";
    s.ytitle = "p_{T} (GeV/c)";
  }
  else if (n.Contains("ctau3D_vs_mass"))
  {
    s.subdir = "h2_ctau3D_vs_mass";
    s.xtitle = "M_{#mu#mu} (GeV/c^{2})";
    s.ytitle = "c#tau_{3D} (mm)";
  }
  else
  {
    s.subdir = "h2_others";
    s.xtitle = "X";
    s.ytitle = "Y";
  }
  return s;
}

// ----- make sub-directories -----
static TString rapFolder(const TString &n)
{
  if (n.Contains("fwd"))
    return "fwd";
  if (n.Contains("mid"))
    return "mid";
  if (n.Contains("all"))
    return "all";
  return "anyrap";
}
static TString rapLabel2D(const TString &n)
{
  if (n.Contains("fwd"))
    return Form("%.1f < |y| < %.1f", kRapAbsMin_Fwd, kRapAbsMax_Fwd);
  if (n.Contains("mid"))
    return Form("|y| < %.1f", kRapAbsMax_Mid);
  if (n.Contains("all"))
    return Form("|y| < %.1f", kRapAbsMax_All);
  return "";
}

// main macro
void drawH2All(const char *file,
               const char *dir = "h2",
               bool savePdf = false,
               bool logz = false,
               const char *outroot = "figs")
{
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  // ----- inputs -----
  std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERROR] File open error\n";
    return;
  }

  TDirectory *d = dynamic_cast<TDirectory *>(f->Get(dir));
  if (!d)
  {
    std::cerr << "[ERROR] No dir: " << dir << "\n";
    return;
  }

  size_t nScanned = 0, nDrawn = 0;

  // main loop - iterate all histograms
  TIter nextkey(d->GetListOfKeys());
  while (auto *key = (TKey *)nextkey())
  {
    ++nScanned;
    TString name = key->GetName();
    TString full = Form("%s/%s", dir, name.Data());

    TObject *obj = f->Get(full);
    TH2 *h2 = dynamic_cast<TH2 *>(obj);
    if (!h2)
      continue;

    // decide labels - directory, x-y axes 
    H2Spec spec = deduceH2Spec(name);
    TString rapDir = rapFolder(name);

    // make output folders
    TString outdir = Form("%s/%s/%s", outroot, spec.subdir.Data(), rapDir.Data());
    gSystem->mkdir(outdir, true);

    TCanvas *c = new TCanvas(Form("c2_%s", name.Data()), "", 800, 800);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.15); // space for color bar
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    
    if (logz)
      c->SetLogz();

    // titles
    h2->GetXaxis()->SetTitle(spec.xtitle);
    h2->GetYaxis()->SetTitle(spec.ytitle);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();

    // set low limit of Z axis (if logz, lowest > 0)
    if (logz) 
    {
      double minpos = 0;
      const int nx = h2->GetNbinsX(), ny = h2->GetNbinsY();
      for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
        {
          double v = h2->GetBinContent(ix, iy);
          if (v > 0.0 && (minpos == 0 || v < minpos))
            minpos = v;
        }
      double zfloor = 1e-6;
      if (minpos > 0 && minpos < zfloor)
        zfloor = minpos * 0.5;
      h2->SetMinimum(zfloor);
    }
    else
    {
      h2->SetMinimum(0.0000001); // not to draw empty bins
      // h2->SetMinimum(0);
    }

    h2->Draw("COLZ");

    // latex - Sample name, rapidities
    TLatex tx;
    tx.SetNDC();
    tx.SetTextAlign(13);
    tx.SetTextSize(0.032);
    tx.DrawLatex(0.14, 0.92, kSampleName);
    TString rlab = rapLabel2D(name);
    if (!rlab.IsNull())
    {
      tx.SetTextSize(0.030);
      tx.DrawLatex(0.14, 0.87, rlab);
    }

    // save
    c->Modified();
    c->Update();
    c->SaveAs(Form("%s/%s.png", outdir.Data(), name.Data()));
    if (savePdf)
      c->SaveAs(Form("%s/%s.pdf", outdir.Data(), name.Data()));
    c->Close();
    delete c;

    ++nDrawn;
  }

  std::cout << "[drawH2All] scanned=" << nScanned
            << " drawn=" << nDrawn
            << " (savePdf=" << (savePdf ? "true" : "false")
            << ", logz=" << (logz ? "true" : "false")
            << ", outroot=" << outroot << ")\n";
}
