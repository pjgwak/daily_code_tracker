#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <initializer_list>
#include <iostream>

namespace DrawEpOverlayLocal
{
const int kColorTracker = kBlack;
const int kColorHFp = kRed + 1;
const int kColorHFm = kBlue + 1;
const int kColorHF = kGreen + 2;

struct InputConfig
{
  TString sample;
  TString label;
  TString path;
};

TString ProjectDir()
{
  // Start from this macro path.
  TString path = __FILE__;
  // Make relative macro paths usable from the current working directory.
  if (!gSystem->IsAbsoluteFileName(path))
    path = TString(gSystem->WorkingDirectory()) + "/" + path;

  // drawing/event_plane -> pilot_calcTrack_3subEvent.
  TString project = TString(gSystem->DirName(path)) + "/../..";
  // Canonicalize symlinks and relative path components.
  TString out = gSystem->GetFromPipe(Form("cd %s && pwd -P", project.Data()));
  out = out.Strip(TString::kBoth);
  return out.IsNull() ? project : out;
}

TString OutputDir(const TString &sample, const TString &category)
{
  TString out = ProjectDir() + "/drawing/" + category + "/figs/" + sample;
  gSystem->mkdir(out, true);
  return out;
}

TString NormalizeSampleTag(const TString &tag)
{
  TString sample(tag);
  sample.ToLower();
  if (sample == "prompt" || sample == "pr")
    return "pr";
  if (sample == "nonprompt" || sample == "np")
    return "np";
  if (sample == "data" || sample == "data23" || sample == "pbpb2023")
    return "data";
  return sample;
}

TString ExistingFile(const TString &path)
{
  if (!path.IsNull() && !gSystem->AccessPathName(path))
    return path;
  return "";
}

TString FirstExisting(std::initializer_list<TString> paths)
{
  for (const auto &path : paths)
  {
    TString existing = ExistingFile(path);
    if (!existing.IsNull())
      return existing;
  }
  return "";
}

TString InputFileForSample(const TString &sample)
{
  const TString project = ProjectDir();

  if (sample == "pr")
    return FirstExisting({project + "/roots/pr/skim_PbPb2023_isMC1_PR_EP_3Sub_job.root"});

  if (sample == "np")
    return FirstExisting({project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test40k.root",
                          // project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test10k.root",
                          // project + "/roots/np/skim_PbPb2023_isMC1_NP_EP_3Sub_test1k.root"
                        });

  if (sample == "data")
    return FirstExisting({project + "/roots/data/Oniatree_MC0_PbPb_2023_EP_3Sub_job.root",
                          // project + "/roots/data/Oniatree_MC0_PbPb_2023_EP_angles_3Sub_job.root"
                        });

  return "";
}

InputConfig GetInputConfig(const char *sampleTag)
{
  InputConfig cfg;
  // Normalize aliases like prompt/nonprompt/data23.
  cfg.sample = NormalizeSampleTag(sampleTag && sampleTag[0] ? sampleTag : "pr");
  // Pick the plot label for this sample.
  cfg.label = cfg.sample == "pr" ? "MC PR" : (cfg.sample == "np" ? "MC NP" : (cfg.sample == "data" ? "Data" : cfg.sample));
  // Pick the ROOT input file under roots/<sample>.
  cfg.path = InputFileForSample(cfg.sample);

  if (cfg.sample != "pr" && cfg.sample != "np" && cfg.sample != "data")
  {
    std::cout << "[ERROR] unknown sample tag: " << cfg.sample << " (use pr, np, or data)" << std::endl;
    return cfg;
  }

  if (cfg.path.IsNull())
    std::cout << "[ERROR] no ROOT input configured/found for sample tag '" << cfg.sample
              << "' under " << ProjectDir() << "/roots/" << cfg.sample << std::endl;
  else
    std::cout << "[INFO] sample = " << cfg.sample << ", input = " << cfg.path << std::endl;
  return cfg;
}

void ApplyStyle()
{
  gROOT->SetBatch(kTRUE);
  const char *rootlogon = "/data/users/pjgwak/input_files/rootlogon.C";
  if (!gSystem->AccessPathName(rootlogon))
    gROOT->ProcessLine(Form(".X %s", rootlogon));

  gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetEndErrorSize(4);
}

void SetupPad(bool is2D = false)
{
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(is2D ? 0.16 : 0.05);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);
  gPad->SetTicks(1, 1);
}

void FormatAxis(TAxis *axis, bool largePanel = false)
{
  if (!axis)
    return;
  axis->CenterTitle(true);
  axis->SetTitleSize(largePanel ? 0.050 : 0.046);
  axis->SetLabelSize(largePanel ? 0.042 : 0.040);
}

void FormatHist(TH1 *h, const char *xTitle, const char *yTitle = "Entries",
                int color = kBlack, int marker = 20)
{
  if (!h)
    return;
  h->SetStats(0);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(1.15);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  FormatAxis(h->GetXaxis());
  FormatAxis(h->GetYaxis());
  h->GetYaxis()->SetTitleOffset(1.35);
}

void SetYHeadroom(TH1 *h, double scale = 1.30, double floor = 1.0)
{
  if (!h)
    return;
  const double maxY = h->GetMaximum();
  h->SetMaximum(maxY > 0.0 ? maxY * scale : floor);
}

TLegend *MakeLegend(double x1, double y1, double x2, double y2, double textSize = 0.030)
{
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(textSize);
  return leg;
}

void DrawCmsLabels(const TString &sample, const char *detail = "",
                   const char *collision = "PbPb 5.36 TeV")
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.032);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.95, 0.93, collision);

  TLatex cms;
  cms.SetNDC();
  cms.SetTextSize(0.040);
  cms.SetTextFont(62);
  cms.SetTextAlign(11);
  cms.DrawLatex(0.205, 0.93, "CMS");

  TLatex extra;
  extra.SetNDC();
  extra.SetTextSize(0.040);
  extra.SetTextFont(52);
  extra.SetTextAlign(11);
  extra.DrawLatex(0.289, 0.93, "Internal");

  TLatex sub;
  sub.SetNDC();
  sub.SetTextFont(42);
  sub.SetTextSize(0.030);
  sub.SetTextAlign(11);
  const TString label = detail && detail[0] ? Form("%s, %s", sample.Data(), detail) : sample;
  sub.DrawLatex(0.18, 0.865, label);
}

void SaveCanvas(TCanvas *c, const TString &sample, const TString &category, const TString &name)
{
  if (!c)
    return;
  const TString outDir = OutputDir(sample, category);
  c->SaveAs(Form("%s/%s.pdf", outDir.Data(), name.Data()));
}

TFile *OpenInput(const TString &inputFile, TTree *&tree)
{
  tree = nullptr;
  if (inputFile.IsNull())
    return nullptr;
  TFile *f = TFile::Open(inputFile, "READ");
  if (!f || f->IsZombie())
  {
    std::cout << "[ERROR] cannot open input file: " << inputFile << std::endl;
    return nullptr;
  }
  tree = (TTree *)f->Get("hionia/myTree");
  if (!tree)
  {
    std::cout << "[ERROR] hionia/myTree not found in " << inputFile << std::endl;
    f->Close();
    return nullptr;
  }
  return f;
}

// selection cut in TTree
const char *ValidTrk()
{
  return "has_trk_ep && ep>-90 && ep_rec>-90 && ep_flat>-90";
}

const char *ValidHF()
{
  return "ep_hfp2>-90 && ep_hfm2>-90";
}

TString HistNameFromSpec(const char *histSpec)
{
  TString histName(histSpec);
  const Ssiz_t pos = histName.First('(');
  if (pos != kNPOS)
    histName.Resize(pos);
  return histName;
}

TH1 *MakeTreeHist(TTree *t, const char *expr, const char *histSpec, const char *cut,
                  const char *xTitle, int color, int marker = 20)
{
  if (!t)
    return nullptr;
  t->Draw(Form("%s>>%s", expr, histSpec), cut, "goff");
  const TString histName = HistNameFromSpec(histSpec);
  TH1 *h = (TH1 *)gDirectory->Get(histName);
  if (!h)
  {
    std::cout << "[WARN] missing histogram: " << histName << std::endl;
    return nullptr;
  }
  FormatHist(h, xTitle, "Entries", color, marker);
  SetYHeadroom(h, 1.25);
  return h;
}

double MaxHistY(std::initializer_list<TH1 *> hists)
{
  double ymax = 0.0;
  for (TH1 *h : hists)
    if (h)
      ymax = std::max(ymax, h->GetMaximum());
  return ymax;
}

}



void draw_ep_overlay(const char *sampleTag = "pr")
{
  using namespace DrawEpOverlayLocal;
  ApplyStyle();

  InputConfig input = GetInputConfig(sampleTag);
  TString sample(input.sample);
  TString sampleLabel(input.label);
  TTree *t = nullptr;
  TFile *f = OpenInput(input.path, t);
  if (!f)
    return;

  TCanvas *c = new TCanvas("cEpOverlay", "EP angle overlay", 800, 800);
  SetupPad();
  TH1 *HFpm = MakeTreeHist(t, "rpAng[2]", "h_overlay_epang(100,-1.5,1.5)", "nEP>2", "#Psi_{2} (rad)", kColorHF, 22);
  TH1 *hHFp = MakeTreeHist(t, "ep_hfp2", "h_overlay_hfp(100,-1.5,1.5)", ValidHF(), "#Psi_{2} (rad)", kColorHFp, 20);
  TH1 *hHFm = MakeTreeHist(t, "ep_hfm2", "h_overlay_hfm(100,-1.5,1.5)", ValidHF(), "#Psi_{2} (rad)", kColorHFm, 21);
  TH1 *hCalTrack = MakeTreeHist(t, "ep_flat", "h_overlay_caltrk(100,-1.5,1.5)", ValidTrk(), "#Psi_{2} (rad)", kColorTracker, 24);

  if (!HFpm || !hHFp || !hHFm || !hCalTrack)
  {
    std::cout << "[ERROR] missing histogram for EP overlay" << std::endl;
    f->Close();
    return;
  }

  const double ymax = MaxHistY({HFpm, hHFp, hHFm, hCalTrack});
  HFpm->SetTitle(";#Psi_{2} (rad);Entries");
  HFpm->SetMaximum(ymax > 0.0 ? ymax * 1.25 : 1.0);
  HFpm->Draw("hist");
  hHFp->Draw("hist same");
  hHFm->Draw("hist same");
  hCalTrack->Draw("hist same");

  TLegend *leg = MakeLegend(0.72, 0.74, 0.88, 0.90);
  leg->AddEntry(HFpm, "HF_{+,-}", "l");
  leg->AddEntry(hHFp, "HF_{+}", "l");
  leg->AddEntry(hHFm, "HF_{-}", "l");
  leg->AddEntry(hCalTrack, "CalTrack", "l");
  leg->Draw();
  DrawCmsLabels(sampleLabel, "");

  SaveCanvas(c, sample, "event_plane", "ep_overlay");
  f->Close();
}
