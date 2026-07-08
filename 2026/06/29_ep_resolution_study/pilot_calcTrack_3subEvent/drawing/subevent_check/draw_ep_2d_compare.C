#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <initializer_list>
#include <iostream>

namespace DrawEp2DCompareLocal
{
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

const char *ValidAll()
{
  return "has_trk_ep && ep_flat>-90 && ep_hfp2>-90 && ep_hfm2>-90";
}

TString HistNameFromSpec(const char *histSpec)
{
  TString histName(histSpec);
  const Ssiz_t pos = histName.First('(');
  if (pos != kNPOS)
    histName.Resize(pos);
  return histName;
}

TH2 *DrawTreeHist2D(TTree *t, const char *expr, const char *histSpec, const char *cut,
                    const char *xTitle, const char *yTitle, const char *panelLabel)
{
  if (!t)
    return nullptr;
  gPad->SetRightMargin(0.16);
  t->Draw(Form("%s>>%s", expr, histSpec), cut, "goff");
  const TString histName = HistNameFromSpec(histSpec);
  TH2 *h = (TH2 *)gDirectory->Get(histName);
  if (!h)
  {
    std::cout << "[WARN] missing histogram: " << histName << std::endl;
    return nullptr;
  }
  h->SetStats(0);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  FormatAxis(h->GetXaxis());
  FormatAxis(h->GetYaxis());
  h->GetYaxis()->SetTitleOffset(1.35);
  h->Draw("colz");
  DrawPanelLabel(panelLabel, 0.18, 0.86);
  return h;
}

}



void draw_ep_2d_compare(const char *sampleTag = "pr")
{
  using namespace DrawEp2DCompareLocal;
  ApplyStyle();

  InputConfig input = GetInputConfig(sampleTag);
  TString sample(input.sample);
  TString sampleLabel(input.label);
  TTree *t = nullptr;
  TFile *f = OpenInput(input.path, t);
  if (!f)
    return;

  TCanvas *c = new TCanvas("cEp2DCompare", "2D EP angle comparisons", 1800, 600);
  c->Divide(3, 1);
  c->cd(1);
  SetupPad(true);
  DrawTreeHist2D(t, "ep_hfp2:ep_flat", "h2_trk_hfp(80,-1.5,1.5,80,-1.5,1.5)", ValidAll(), "#Psi_{2,trk}^{calc}", "#Psi_{2,HF_{+}}", "trk vs HF_{+}");
  DrawCmsLabels(sampleLabel, "2D EP comparison");
  c->cd(2);
  SetupPad(true);
  DrawTreeHist2D(t, "ep_hfm2:ep_flat", "h2_trk_hfm(80,-1.5,1.5,80,-1.5,1.5)", ValidAll(), "#Psi_{2,trk}^{calc}", "#Psi_{2,HF_{-}}", "trk vs HF_{-}");
  c->cd(3);
  SetupPad(true);
  DrawTreeHist2D(t, "rpAng[3]:ep_flat", "h2_trk_input(80,-1.5,1.5,80,-1.5,1.5)", "has_trk_ep && ep_flat>-90 && nEP>3", "#Psi_{2,trk}^{calc}", "#Psi_{2,trk}^{input}", "calc vs input trk");

  SaveCanvas(c, sample, "subevent_check", "ep_2d_detector_comparison");
  f->Close();
}
