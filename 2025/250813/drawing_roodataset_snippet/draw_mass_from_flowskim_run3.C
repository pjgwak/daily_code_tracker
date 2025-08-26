void draw_mass_from_flowskim_run3()
{
  TFile f("/work/pjgwak/pol24/input_flowskim/roots/OniaFlowSkim_isMC0_Run3PbPb_Minbias.root");
  auto t = (TTree *)f.Get("myTree"); // tree name
  const char *mass = "mass";         // branch name
  const int nb = 80;
  double xmin = 2.6, xmax = 3.5;
  const double lo = 2.6, hi = 3.5;

  t->Draw(Form("%s>>h(%d,%g,%g)", mass, nb, xmin, xmax), "", "hist");
  auto h = (TH1D *)gDirectory->Get("h");

  // draw
  TCanvas c("c_mass", "mass", 800, 600);
  h->Draw("hist");
  c.SaveAs("mass_in_flowskim_Run3.png");
  printf("[Saved] %s\n", "mass_in_flowskim_Run3.png");

  // count
  Long64_t n = t->GetEntries(Form("%s>=%.9g && %s<=%.9g", mass, lo, mass, hi));
  printf("Entries in [%.3f, %.3f]: %lld\n", lo, hi, n);
}