#include "PlotCore.C"

void drawOne(const char *file, const char *path, bool savePdf = false)
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
  
  plotOne(f.get(), path, "figs_test", savePdf);
}