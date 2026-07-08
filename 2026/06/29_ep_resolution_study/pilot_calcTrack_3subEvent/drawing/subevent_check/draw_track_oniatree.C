#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <initializer_list>
#include <iostream>

namespace DrawTrackOniatreeLocal
{
const int kColorTrackMid = kMagenta + 1;
const int kColorTrackM = kBlue + 1;
const int kColorTrackP = kRed + 1;

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

  // drawing/subevent_check -> pilot_calcTrack_3subEvent.
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

void DrawPanelLabel(const char *text, double x = 0.18, double y = 0.86,
                    double size = 0.034)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(size);
  latex.DrawLatex(x, y, text);
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

TString HistNameFromSpec(const char *histSpec)
{
  TString histName(histSpec);
  const Ssiz_t pos = histName.First('(');
  if (pos != kNPOS)
    histName.Resize(pos);
  return histName;
}

TH1 *DrawTreeHist(TTree *t, const char *expr, const char *histSpec, const char *cut,
                  const char *xTitle, const char *panelLabel, int color = kBlack,
                  double labelX = 0.18, double labelY = 0.86)
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
  FormatHist(h, xTitle, "Entries", color);
  SetYHeadroom(h, 1.35);
  h->Draw("hist");
  DrawPanelLabel(panelLabel, labelX, labelY);
  return h;
}

}



void draw_track_oniatree(const char *sampleTag = "pr")
{
  using namespace DrawTrackOniatreeLocal;
  ApplyStyle();

  InputConfig input = GetInputConfig(sampleTag);
  TString sample(input.sample);
  TString sampleLabel(input.label);
  TTree *t = nullptr;
  TFile *f = OpenInput(input.path, t);
  if (!f)
    return;

  TCanvas *c = new TCanvas("cTrackOniatree", "tracker EP angles from Oniatree", 1800, 600);
  c->Divide(3, 1);
  c->cd(1);
  SetupPad();
  DrawTreeHist(t, "rpAng[3]", "h_ep_trackmid_oniatree(100,-1.5,1.5)", "nEP>3", "#Psi_{2,trackmid} (rad)", "Track_{mid} (Oniatree)", kColorTrackMid, 0.18, 0.82);
  DrawCmsLabels(sampleLabel, "");
  c->cd(2);
  SetupPad();
  DrawTreeHist(t, "rpAng[4]", "h_ep_trackm_oniatree(100,-1.5,1.5)", "nEP>4", "#Psi_{2,track-} (rad)", "Track_{-} (Oniatree)", kColorTrackM, 0.18, 0.86);
  c->cd(3);
  SetupPad();
  DrawTreeHist(t, "rpAng[5]", "h_ep_trackp_oniatree(100,-1.5,1.5)", "nEP>5", "#Psi_{2,track+} (rad)", "Track_{+} (Oniatree)", kColorTrackP, 0.18, 0.86);

  SaveCanvas(c, sample, "subevent_check", "track_oniatree");
  f->Close();
}
