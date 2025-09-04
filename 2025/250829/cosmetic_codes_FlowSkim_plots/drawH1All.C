#include "PlotCore.C"

void drawH1All(const char *file, const char *dir = "h1", bool savePdf = false)
{
  gROOT->Macro("/data/users/pjgwak/input_files/rootlogon.C");
  gStyle->SetOptStat(0);

  // open file
  std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "File open error\n";
    return;
  }

  // open TDirectory
  TDirectory *d = dynamic_cast<TDirectory *>(f->Get(dir));
  if (!d)
  {
    std::cerr << "no dir: " << dir << "\n";
    return;
  }

  // iterate all histogram in the directory
  TIter nextkey(d->GetListOfKeys());
  while (auto key = (TKey *)nextkey())
  {
    TString name = key->GetName();
    TString fullpath = Form("%s/%s", dir, name.Data());
    
    // draw
    TObject *obj = f->Get(fullpath);
    if (dynamic_cast<TH2 *>(obj))
      continue; // don't draw 2D plot
    if (dynamic_cast<TH1 *>(obj))
    {
      plotOne(f.get(), fullpath, "figs", savePdf); // 3rd parameter: output folder name
    }
  }
}