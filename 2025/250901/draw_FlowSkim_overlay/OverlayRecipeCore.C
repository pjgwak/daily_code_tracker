// #include "PlotCore.C"
#include "/data/users/pjgwak/work/daily_code_tracker/2025/250829/cosmetic_codes_FlowSkim_plots/PlotCore.C"
#include <TROOT.h>
#include <TStyle.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>


struct RecipeItem
{
  TString path;          // "h1/mass_region6_PR_SR_all_pt9to12"
  TString label;         // legend
  int color = 0;         // 0: automatic change
  int marker = 0;        // 0: automatic change
  int lstyle = 1;        // 1: solid, 2: dotted, 3: solid-dot etc
  int lwidth = 2;        // line width
  int rebin = 1;         // 1: no rebin
  double scale = 1.0;    // 1: no scaling
  TString drawOpt = "E"; // "E", "hist", "E1" etc
  bool visible = true;
};

// === core functions ===
void runOverlayRecipe(const char *file,
                      const std::vector<RecipeItem> &items,
                      const char *outroot = "figs_overlay",
                      const char *outname = "overlay_custom",
                      bool savePdf = false,
                      bool normalize = false, // Integral("width")=1
                      bool logy = false,
                      double xMin = std::numeric_limits<double>::quiet_NaN(),
                      double xMax = std::numeric_limits<double>::quiet_NaN(),
                      const char *title1 = nullptr, // latex 1,2,3
                      const char *title2 = nullptr,
                      const char *title3 = nullptr)
{
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // open file
  std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERROR] File open error: " << file << "\n";
    return;
  }

  // collect items
  std::vector<TH1 *> hs;
  hs.reserve(items.size());
  std::vector<DrawOpt> opts;
  opts.reserve(items.size());
  std::vector<TString> names;
  names.reserve(items.size());

  // choose color and marker
  auto autoColor = [&](int i)
  {
    static const int cols[] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                               kOrange + 7, kAzure + 1, kViolet + 1, kCyan + 2, kPink + 9, kGray + 2};
    // if i > array, go back to first
    return cols[i % (int)(sizeof(cols) / sizeof(cols[0]))];
  };
  auto autoMarker = [&](int i)
  {
    static const int mks[] = {20, 21, 22, 23, 24, 25, 26, 33, 27, 28};
    
    // if i > array, go back to first
    return mks[i % (int)(sizeof(mks) / sizeof(mks[0]))];
  };

  TH1 *href = nullptr;
  DrawOpt oref;
  TString xtitle_ref, ytitle_ref;

  // === set drawing style ===
  for (size_t i = 0; i < items.size(); ++i)
  {
    if (!items[i].visible) continue;
    
    TObject *obj = f->Get(items[i].path);
    TH1 *src = dynamic_cast<TH1 *>(obj);
    if (!src)
    {
      std::cerr << "[WARN] not TH1: " << items[i].path << "\n";
      continue;
    }

    // clone histogram
    TH1 *h = (TH1 *)src->Clone(Form("ovR_%zu", i)); // %zu: positive integer
    h->SetDirectory(nullptr); // detach histgoram from the TFile
    h->SetStats(0);

    // rebin and scaling - before normalization
    if (items[i].rebin > 1) h->Rebin(items[i].rebin);
    if (items[i].scale != 1.0) h->Scale(items[i].scale);

    // --- use axis of first histogram ---
    if (!href)
    {
      href = h;
      PlotCtx c0 = analyzeName(src->GetName()); // functions: refer to PlotCore.C
      oref = deduceOpt(c0);
      xtitle_ref = oref.xtitle;
      ytitle_ref = makeYTitle(oref, href);
    }

    // --- normalization ---
    if (normalize)
    {
      int nb = h->GetNbinsX();
      double integ = h->Integral(1, nb, "width");
      if (integ > 0)
        h->Scale(1.0 / integ);
      ytitle_ref = "Normalized";
    }

    // --- cosmetics ---
    int col = items[i].color ? items[i].color : autoColor((int)i);
    int mks = items[i].marker ? items[i].marker : autoMarker((int)i);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(mks);
    h->SetLineStyle(items[i].lstyle);
    h->SetLineWidth(items[i].lwidth);
    h->SetMarkerSize(0.9);

    hs.push_back(h);
    opts.push_back(oref);
    names.push_back(items[i].label.IsNull() ? TString(src->GetName()) : items[i].label);
  }

  if (hs.empty())
  {
    std::cerr << "[ERROR] No histogram chosen\n";
    return;
  }

  // === canvans ===
  TCanvas *c = new TCanvas("c_recipe", "c_recipe", 800, 600);
  c->SetMargin(0.12, 0.05, 0.12, 0.05);
  if (logy) c->SetLogy();

  // --- x axis: range, title
  if (oref.setX && oref.xmax > oref.xmin)
    href->GetXaxis()->SetRangeUser(oref.xmin, oref.xmax);
  if (xMin == xMin && xMax == xMax)
    href->GetXaxis()->SetRangeUser(xMin, xMax);
  if (!xtitle_ref.IsNull())
    href->GetXaxis()->SetTitle(xtitle_ref);
  href->GetYaxis()->SetTitle(ytitle_ref);
  href->GetXaxis()->CenterTitle();
  href->GetYaxis()->CenterTitle();

  // --- y axis: logY style
  double ymax = 0, minpos = 1e300; // ymax and ymin candidate
  for (auto *h : hs) // check all hist
  {
    // ymax
    ymax = std::max(ymax, h->GetMaximum());
    
    // ymin (should > 0)
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
    {
      double v = h->GetBinContent(ib);
      if (v > 0.0 && v < minpos)
        minpos = v;
    }
  }
  if (logy)
  {
    double floor = 1e-3;
    double ymin = (minpos < 1e300) ? std::max(minpos * 0.5, floor) : floor;
    href->SetMinimum(ymin);
    href->SetMaximum((ymax > 0 ? ymax * 10. : ymin * 100.));
  }
  else
  {
    href->SetMinimum(0.0);
    href->SetMaximum(ymax > 0 ? ymax * 1.5 : 1.0);
  }

  // --- draw histograms ---
  hs.front()->Draw(items[0].drawOpt.IsNull() ? "E" : items[0].drawOpt.Data());
  for (size_t i = 1; i < hs.size(); ++i)
  {
    TString opt = (i < items.size() && !items[i].drawOpt.IsNull()) ? items[i].drawOpt : "E SAME";
    if (!opt.Contains("SAME"))
      opt += " SAME";
    hs[i]->Draw(opt);
  }

  // --- legend ---
  TLegend *leg = new TLegend(0.15, 0.70, 0.3, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.025);
  for (size_t i = 0; i < hs.size(); ++i)
    leg->AddEntry(hs[i], names[i], "pe"); // point + vertical error, lpe: point+x+y error
    // leg->AddEntry(hs[i], names[i], "pe");
  leg->Draw();

  // --- latex 1,2,3, ---
  TLatex tx;
  tx.SetNDC();
  tx.SetTextAlign(33);
  double X = 0.93, Y = 0.91, DY = 0.05;
  tx.SetTextSize(0.032);

  // draw latex
  tx.DrawLatex(X, Y, title1 ? title1 : kSampleName);
  if (title2)
  {
    tx.SetTextSize(0.030);
    tx.DrawLatex(X, Y - DY, title2);
  }
  if (title3)
  {
    tx.SetTextSize(0.030);
    tx.DrawLatex(X, Y - 2 * DY, title3);
  }

  // --- save ---
  gSystem->mkdir(outroot, true);
  TString png = Form("%s/%s.png", outroot, outname);
  c->SaveAs(png);
  if (savePdf) c->SaveAs(Form("%s/%s.pdf", outroot, outname));

  // --- clean the code ---
  c->Close();
  delete c;
  for (auto *h : hs)
    delete h; // close all histograms
}
