#include "PlotCore.C"
#include <vector>
#include <iostream>
#include <TObjArray.h>
#include <TObjString.h>


// "a,b,c" -> vector = {"a","b","c"}
static std::vector<TString> splitCSV(const char *csv)
{
  std::vector<TString> out;
  if (!csv) return out;
  
  TString s(csv);
  s.ReplaceAll(" ", ""); // remove white space
  if (s.IsNull()) return out;

  TObjArray *arr = s.Tokenize(","); // distinguish word by comma
  if (!arr) return out;

  // casting token to TString and push to vector
  for (int i = 0; i < arr->GetEntries(); ++i)
  {
    auto *tok = dynamic_cast<TObjString *>(arr->At(i));
    if (!tok) continue;
    
    TString t = tok->GetString();

    if (!t.IsNull())
      out.push_back(t);
  }
  delete arr;
  return out;
}

// pass: have all name token in include and don't have any name token in exclude
static bool passFilterName(const TString &n,
                           const std::vector<TString> &inc,
                           const std::vector<TString> &exc)
{
  for (const auto &s : inc)
    if (!n.Contains(s))
      return false;
  for (const auto &s : exc)
    if (n.Contains(s))
      return false;
  return true;
}

// root -l -b -q 'drawH1Filter.C+("/path/hists.root","h1","mass_,PR,SR","LSB,RSB",false,false,"figs")'
void drawH1Filter(const char *file,
                  const char *dir = "h1",
                  const char *includeCSV = "",  // ex: "mass_,PR,SR"
                  const char *excludeCSV = "",  // ex: "LSB,RSB"
                  bool savePdf = false,         
                  bool dryRun = false,          // true: print target histogram list
                  const char *outroot = "figs") // output folder
{
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  // open file
  std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[ERROR] File open error\n";
    return;
  }

  // open TDirectory
  TDirectory *d = dynamic_cast<TDirectory *>(f->Get(dir));
  if (!d)
  {
    std::cerr << "[ERROR] No dir: " << dir << "\n";
    return;
  }

  // make name tokens for filtering
  auto inc = splitCSV(includeCSV);
  auto exc = splitCSV(excludeCSV);

  size_t nScanned = 0, nMatched = 0, nDrawn = 0;

  // main loop - interate all histograms
  TIter nextkey(d->GetListOfKeys());
  while (auto *key = (TKey *)nextkey())
  {
    ++nScanned;
    TString name = key->GetName();
    TString full = Form("%s/%s", dir, name.Data());

    TObject *obj = f->Get(full);
    if (!obj)
      continue;

    // skip TH2D -> use drawH2All
    if (dynamic_cast<TH2 *>(obj))
      continue;
    
    TH1 *h1 = dynamic_cast<TH1 *>(obj);
    if (!h1)
      continue;

    // check pattern
    if (!passFilterName(name, inc, exc))
      continue;
    ++nMatched;

    if (dryRun)
    {
      std::cout << "[MATCH] " << name << "\n";
      continue;
    }

    // draw
    plotOne(f.get(), full, outroot, savePdf);
    ++nDrawn;
  }

  std::cout << "[drawH1Filter] scanned=" << nScanned
            << " matched=" << nMatched
            << " drawn=" << nDrawn
            << " | savePdf=" << (savePdf ? "true" : "false")
            << " dryRun=" << (dryRun ? "true" : "false")
            << " outroot=" << outroot << "\n";
}