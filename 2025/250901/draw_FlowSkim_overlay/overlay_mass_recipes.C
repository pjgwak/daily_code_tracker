#include "OverlayRecipeCore.C"

// make new function and fill the list of histograms liek below
// e.g. compare many pT bins
void mass_PR_SR_all_byPt(const char *file,
                         bool savePdf = false, bool normalize = false, bool logy = false)
{
  // histogram name, legend, color, marker style, line style, line width, rebin=1, scale=1.0, drawOpt

  std::vector<RecipeItem> v = {
      {"h1/mass_region6_PR_SR_all_pt6p5to9", "6.5 < p_{T} < 9 GeV", kBlue + 1, 20, 1, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt9to12", "9 < p_{T} < 12 GeV", kRed + 1, 21, 2, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt12to15", "12 < p_{T} < 15 GeV", kGreen + 2, 22, 3, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt15to20", "15 < p_{T} < 20 GeV", kMagenta + 1, 23, 2, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt20to30", "20 < p_{T} < 30 GeV", kOrange + 7, 24, 3, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt30to50", "30 < p_{T} < 50 GeV", kAzure + 1, 25, 1, 2, 1, 1.0, "E"},
  };

  runOverlayRecipe(file, v,
                   /*outroot=*/"figs_overlay/mass",
                   /*outname=*/"mass_region6_PR_SR_all_byPt",
                   savePdf, normalize, logy,
                   /*xMin=*/2.6, /*xMax=*/3.5,
                   /*legend1=*/"PbPb2023 Data",
                   /*legend2=*/"|y| < 2.4",
                   /*legend3=*/"PR SR");
}

// mass: compare different pT bins in same rapidity
void mass_PR_SR_rap_pt9to12(const char *file, bool savePdf = false, bool normalize = false, bool logy = false)
{
  std::vector<RecipeItem> v = {
      {"h1/mass_region6_PR_SR_fwd_pt9to12", "1.6 < |y| < 2.4, 9 < p_{T} < 12 GeV", kRed + 1, 21, 1, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_mid_pt9to12", "|y| < 1.6", kBlue + 1, 20, 2, 2, 1, 1.0, "E"},
      {"h1/mass_region6_PR_SR_all_pt9to12", "|y| < 2.4", kGreen + 2, 22, 3, 2, 1, 1.0, "E"},
  };
  runOverlayRecipe(file, v, "figs_overlay/mass",
                   "mass_region6_PR_SR_rap_pt9to12",
                   savePdf, normalize, logy,
                   2.6, 3.5, "PbPb2023 Data", "9 < p_{T} < 12 GeV", "PR SR");
}

// === main macro ===
void overlay_mass_recipes(const char *file, bool savePdf = false, bool normalize = false, bool logy = false)
{
  mass_PR_SR_all_byPt(file, savePdf, normalize, logy);
  mass_PR_SR_rap_pt9to12(file, savePdf, normalize, logy);
}
