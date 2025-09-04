#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <regex>

// ===== basic setting - sample name, rapidity ranges =====
static const char *kSampleName = "PbPb2023 Data";
static const double kRapAbsMax_All = 2.4;
static const double kRapAbsMax_Mid = 1.6;
static const double kRapAbsMin_Fwd = 1.6;
static const double kRapAbsMax_Fwd = 2.4;

// ===== plot information container - save meta informations =====
struct PlotCtx
{
  TString name; // histogram name
  TString kind; // branch name: mass, ctau, pt, y, etc
  bool isRegion6 = false; // PR RSB, NP SR, etc
  TString comp; // B components: "PR", "NP", ""(empty)
  TString band; // "SR", "LSB", "RSB", ""
  TString rap;  // "all", "mid", "fwd", ""
  bool hasPtBin = false;
  double ptLo = 0, ptHi = 0;
};

// ===== draw option: x-title, range, isLog etc =====
struct DrawOpt
{
  TString xtitle = "", ytitle = "Events";
  TString xunit = "";
  bool logy = false, setX = false;
  double xmin = 0, xmax = 0;
};

// isUniformBin - if yes, print bin width on y-title
inline double uniformBinWidth(const TH1 *h)
{
  int nb = h->GetNbinsX();
  if (nb <= 0)
    return -1;
  double w0 = h->GetXaxis()->GetBinWidth(1);
  for (int i = 2; i <= nb; i++)
    if (fabs(h->GetXaxis()->GetBinWidth(i) - w0) > 1e-9)
      return -1;
  return w0;
}

// set y-title: uniform bin: Events /(bin width), if not Events / bin
inline TString makeYTitle(const DrawOpt &opt, const TH1 *h)
{
  double w = uniformBinWidth(h);
  if (w <= 0)
    return "Events / bin";
  return opt.xunit.Length() ? Form("Events / (%.3g %s)", w, opt.xunit.Data())
                            : Form("Events / (%.3g)", w);
}

// ===== Histogram naem parsing -> fill the PlotCtx =====
inline PlotCtx analyzeName(const TString &n)
{
  PlotCtx c;
  c.name = n;

  // kind
  if (n.Contains("pt_vs_y"))
    c.kind = "h2";
  else if (n.BeginsWith("mass_") || n.Contains("_mass_"))
    c.kind = "mass";
  else if (n.BeginsWith("ctau3D_") || n.Contains("_ctau3D_"))
    c.kind = "ctau";
  else if (n.BeginsWith("pt_") || n.Contains("_pt"))
    c.kind = "pt";
  else if (n.BeginsWith("y_") || n.Contains("_y_"))
    c.kind = "y";
  else
    c.kind = "other";

  // region6
  c.isRegion6 = n.Contains("region6");

  // comp
  if (n.Contains("_PR_"))
    c.comp = "PR";
  else if (n.Contains("_NP_"))
    c.comp = "NP";

  // band
  if (n.EndsWith("_SR") || n.Contains("_SR_"))
    c.band = "SR";
  else if (n.Contains("_LSB"))
    c.band = "LSB";
  else if (n.Contains("_RSB"))
    c.band = "RSB";

  // rap
  if (n.Contains("_fwd"))
    c.rap = "fwd";
  else if (n.Contains("_mid"))
    c.rap = "mid";
  else if (n.Contains("_all"))
    c.rap = "all";

  // pT bin: _pt6p5to9
  std::regex re("_pt([0-9]+p?[0-9]*)to([0-9]+p?[0-9]*)");
  std::cmatch m;
  if (std::regex_search(n.Data(), m, re))
  {
    auto cvt = [](std::string s)
    { for(auto&ch:s) if(ch=='p') ch='.'; return s; };
    c.hasPtBin = true;
    c.ptLo = atof(cvt(m[1].str()).c_str());
    c.ptHi = atof(cvt(m[2].str()).c_str());
  }
  return c;
}

// ===== set drawOpt according to branch type =====
inline DrawOpt deduceOpt(const PlotCtx &c)
{
  DrawOpt o;
  if (c.kind == "mass")
  {
    o.xtitle = "M_{#mu#mu} (GeV/c^{2})";
    o.xunit = "GeV/c^{2}";
    o.setX = true;
    o.xmin = 2.6;
    o.xmax = 3.5;
  }
  else if (c.kind == "ctau")
  {
    o.xtitle = "c#tau_{3D} (mm)";
    o.xunit = "mm";
    o.logy = true;
  }
  else if (c.kind == "pt")
  {
    o.xtitle = "p_{T} (GeV/c)";
    o.xunit = "GeV/c";
  }
  else if (c.kind == "y")
  {
    o.xtitle = "y";
  }
  return o;
}

// ===== make a labels =====
inline TString rapLabel(const PlotCtx &c)
{
  if (c.rap == "all")
    return Form("|y| < %.1f", kRapAbsMax_All);
  if (c.rap == "mid")
    return Form("|y| < %.1f", kRapAbsMax_Mid);
  if (c.rap == "fwd")
    return Form("%.1f < |y| < %.1f", kRapAbsMin_Fwd, kRapAbsMax_Fwd);
  return "";
}

inline TString ptLabel(const PlotCtx &c)
{
  if (c.hasPtBin)
    return Form("%.3g < p_{T} < %.3g GeV/c", c.ptLo, c.ptHi);
  // fallback by rap
  if (c.rap == "fwd")
    return "3 < p_{T} < 50 GeV/c";
  if (c.rap == "mid")
    return "6.5 < p_{T} < 50 GeV/c";
  if (c.rap == "all")
    return "6.5 < p_{T} < 50 GeV/c";
  return "";
}

// combine rap, pT label
inline TString infoLine(const PlotCtx &c)
{
  TString a = rapLabel(c), b = ptLabel(c);
  if (!a.IsNull() && !b.IsNull())
    return Form("%s, %s", a.Data(), b.Data());
  if (!a.IsNull())
    return a;
  if (!b.IsNull())
    return b;
  return "";
}

// ===== common cosmetics =====
inline void applyStyle1D(TH1 *h)
{
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.9);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

// ===== kind/조합별 커스텀(확장 포인트) =====
inline void applyCustomCosmetics(const PlotCtx &c, TH1 *h)
{
  // if (c.kind == "mass" && c.isRegion6)
  // {
  //   h->SetLineColor(kBlue + 1);
  // }
  // if (c.kind == "ctau")
  // {
  //   h->SetLineColor(kRed + 1);
  // }
}

// ===== set yMax: peak * 1.5 (1.2) =====
inline void setYMaxAuto(const PlotCtx &c, TH1 *h, bool logy)
{
  if (!h) return;
  
  const double maxv = h->GetMaximum();

  if (!logy)
  {
    h->SetMinimum(0.0);
    double scale = 1.5;
   
    if (c.band == "LSB" || c.band == "RSB")
      scale = 1.2;
   
    h->SetMaximum(maxv > 0 ? maxv * scale : 1.0);
    return;
  }

  // ----- log scale (ymin > 0) -----
  const double floor = 1e-3;
  double minpos = 1e300;
  bool foundPos = false;
  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    double v = h->GetBinContent(i);
    if (v > 0.0)
    {
      foundPos = true;
      if (v < minpos)
        minpos = v;
    }
  }

  double ymin = floor;
  if (foundPos)
    ymin = std::max(minpos * 0.5, floor);
  h->SetMinimum(ymin);

  double scale = 10.0;
  if (c.band == "LSB" || c.band == "RSB")
    scale = 5.0;

  double ymax = (maxv > 0) ? maxv * scale : ymin * 100.;
  if (ymax <= ymin)
    ymax = ymin * 100.;
  h->SetMaximum(ymax);
}

// ===== label =====
// PR, NP, RSB, SR, etc
inline TString compBandLabel(const PlotCtx &c)
{
  if (!c.isRegion6)
    return "";
  TString s;
  if (c.comp.Length())
    s += c.comp; // PR, NP
  if (c.band.Length())
  {
    if (!s.IsNull())
      s += " ";
    s += c.band; // SR, LSB, RSB
  }
  return s;
}

// draw labels
inline void drawLabelBlock(const PlotCtx &c)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextAlign(33);
  double x = 0.93, y = 0.91, dy = 0.05;
  tx.SetTextSize(0.032);

  // sample name
  tx.DrawLatex(x, y, kSampleName);
  
  // |y|, pT in same line
  TString info = infoLine(c);
  if (!info.IsNull())
  {
    tx.SetTextSize(0.030);
    tx.DrawLatex(x, y - dy, info);
  }

  // region6
  TString cb = compBandLabel(c);
  if (c.isRegion6 && !cb.IsNull())
  {
    tx.SetTextSize(0.030);
    tx.DrawLatex(x, y - 2 * dy, cb);
  }
}

// ===== set directory =====
// top directory
inline TString subdirFor(const PlotCtx &c)
{
  if (c.kind == "h2")
    return "h2_pt_vs_y";
  if (c.kind == "mass" && c.isRegion6)
    return "mass_region6";
  if (c.kind == "ctau" && c.isRegion6)
    return "ctau3D_region6";
  if (c.kind == "mass")
    return "mass";
  if (c.kind == "ctau")
    return "ctau3D";
  if (c.kind == "pt")
    return "pt";
  if (c.kind == "y")
    return "y";
  return "others";
}

// sub-directory
inline TString subdirDeep(const PlotCtx &c)
{
  TString base = subdirFor(c); // mass_region6, ctau3D_region6, ...
  
  // More directories for region6 (mass/ctau)
  if (c.isRegion6 && (c.kind == "mass" || c.kind == "ctau"))
  {
    TString comp = c.comp.Length() ? c.comp : "ANY";
    TString band = c.band.Length() ? c.band : "ANY";
    base += Form("/%s_%s", comp.Data(), band.Data()); // PR_SR, PR_LSB, etc
    TString rap = c.rap.Length() ? c.rap : "anyrap";  // fwd/mid/all/anyrap
    base += Form("/%s", rap.Data());
  }
  return base;
}

// ===== draw one plot =====
inline void plotOne(TFile *f, const char *path,
                    const char *outroot = "figs", bool savePdf = false)
{
  // get hist from the file
  TObject *obj = f->Get(path);
  if (!obj)
  {
    std::cerr << "No object " << path << "\n";
    return;
  }

  // name and labels
  TString name = obj->GetName();
  PlotCtx c = analyzeName(name);

  // make output folders
  TString outdir = Form("%s/%s", outroot, subdirDeep(c).Data());
  gSystem->mkdir(outdir, true);

  // skip TH2D
  if (c.kind == "h2")
    return;

  TH1 *h = dynamic_cast<TH1 *>(obj);
  if (!h)
    return;

  // set draw option
  DrawOpt o = deduceOpt(c);

  TCanvas *can = new TCanvas(Form("c_%s", name.Data()), "", 900, 700);
  can->SetMargin(0.12, 0.05, 0.12, 0.05);
  if (o.logy)
    can->SetLogy();

  // set titles
  if (o.setX)
    h->GetXaxis()->SetRangeUser(o.xmin, o.xmax);
  if (!o.xtitle.IsNull())
    h->GetXaxis()->SetTitle(o.xtitle);
  h->GetYaxis()->SetTitle(makeYTitle(o, h));

  // apply cosmetcis
  applyStyle1D(h);
  applyCustomCosmetics(c, h);
  setYMaxAuto(c, h, o.logy);

  // draw hist
  h->SetStats(0);
  h->Draw("E");

  // draw label
  drawLabelBlock(c);

  // draw canvas
  can->Modified();
  can->Update();

  // save canvas
  can->SaveAs(Form("%s/%s.png", outdir.Data(), name.Data()));
  if (savePdf)
    can->SaveAs(Form("%s/%s.pdf", outdir.Data(), name.Data()));
  can->Close();

  delete can;
}