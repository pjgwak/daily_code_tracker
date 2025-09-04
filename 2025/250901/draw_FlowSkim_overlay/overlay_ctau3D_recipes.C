// overlay_ctau3D_recipes.C
#include "OverlayRecipeCore.C"

// ctau3D_region6: PR RSB, all rap — compare many pT bins
void ctau3D_PR_RSB_all_byPt(const char *file,
                            bool savePdf = false, bool normalize = false, bool logy = true)
{
  std::vector<RecipeItem> v = {
      {"h1/ctau3D_region6_PR_RSB_all_pt6p5to9", "6.5 < p_{T} <9 GeV", kBlue + 1, 20, 1, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_RSB_all_pt9to12", "9 < p_{T} <12 GeV", kRed + 1, 21, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_RSB_all_pt12to15", "12 < p_{T} <15 GeV", kGreen + 2, 22, 3, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_RSB_all_pt15to20", "15 < p_{T} <20 GeV", kMagenta + 1, 23, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_RSB_all_pt20to30", "20 < p_{T} <30 GeV", kOrange + 7, 24, 3, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_RSB_all_pt30to50", "30 < p_{T} <50 GeV", kAzure + 1, 25, 1, 2, 1, 1.0, "e"},
  };

  runOverlayRecipe(file, v,
                   /*outroot=*/"figs_overlay/ctau3D",
                   /*outname=*/"ctau3D_region6_PR_RSB_all_byPt",
                   savePdf, normalize, /*logy=*/true,
                   /*xMin=*/std::numeric_limits<double>::quiet_NaN(),
                   /*xMax=*/std::numeric_limits<double>::quiet_NaN(),
                   "PbPb2023 Data", "|y| < 2.4", "PR RSB");
}

// ctau3D_region6: PR RSB, all rap — compare many pT bins
void ctau3D_NP_RSB_all_byPt(const char *file,
                            bool savePdf = false, bool normalize = false, bool logy = true)
{
  std::vector<RecipeItem> v = {
      {"h1/ctau3D_region6_NP_RSB_all_pt6p5to9", "6.5 < p_{T} <9 GeV", kBlue + 1, 20, 1, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt9to12", "9 < p_{T} <12 GeV", kRed + 1, 21, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt12to15", "12 < p_{T} <15 GeV", kGreen + 2, 22, 3, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt15to20", "15 < p_{T} <20 GeV", kMagenta + 1, 23, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt20to30", "20 < p_{T} <30 GeV", kOrange + 7, 24, 3, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt30to50", "30 < p_{T} <50 GeV", kAzure + 1, 25, 1, 2, 1, 1.0, "e"},
  };

  runOverlayRecipe(file, v,
                   /*outroot=*/"figs_overlay/ctau3D",
                   /*outname=*/"ctau3D_region6_NP_RSB_all_byPt",
                   savePdf, normalize, /*logy=*/true,
                   /*xMin=*/std::numeric_limits<double>::quiet_NaN(),
                   /*xMax=*/std::numeric_limits<double>::quiet_NaN(),
                   "PbPb2023 Data", "|y| < 2.4", "NP RSB");
}

// ctau3D_region6: PR RSB, all rap — compare many pT bins
void ctau3D_NP_RSB_upTo20_byPt(const char *file,
                            bool savePdf = false, bool normalize = false, bool logy = true)
{
  std::vector<RecipeItem> v = {
      {"h1/ctau3D_region6_NP_RSB_all_pt6p5to9", "6.5 < p_{T} <9 GeV", kBlue + 1, 20, 1, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt9to12", "9 < p_{T} <12 GeV", kRed + 1, 21, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt12to15", "12 < p_{T} <15 GeV", kGreen + 2, 22, 3, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_NP_RSB_all_pt15to20", "15 < p_{T} <20 GeV", kMagenta + 1, 23, 2, 2, 1, 1.0, "e"},
  };

  runOverlayRecipe(file, v,
                   /*outroot=*/"figs_overlay/ctau3D",
                   /*outname=*/"ctau3D_region6_NP_RSB_upTo20_byPt",
                   savePdf, normalize, /*logy=*/true,
                   /*xMin=*/std::numeric_limits<double>::quiet_NaN(),
                   /*xMax=*/std::numeric_limits<double>::quiet_NaN(),
                   "PbPb2023 Data", "|y| < 2.4", "NP RSB");
}

// ctau3D: compare rapiditys in the same pT region
void ctau3D_PR_SR_rap_pt20to30(const char *file, bool savePdf = false, bool normalize = false, bool logy = true)
{
  std::vector<RecipeItem> v = {
      {"h1/ctau3D_region6_PR_SR_fwd_pt20to30", "1.6 < |y| < 2.4, 20 < p_{T} < 30 GeV", kRed + 1, 21, 1, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_SR_mid_pt20to30", "|y| < 1.6", kBlue + 1, 20, 2, 2, 1, 1.0, "e"},
      {"h1/ctau3D_region6_PR_SR_all_pt20to30", "|y| < 2.4", kGreen + 2, 22, 3, 2, 1, 1.0, "e"},
  };
  runOverlayRecipe(file, v, "figs_overlay/ctau3D",
                   "ctau3D_region6_PR_SR_rap_pt20to30",
                   savePdf, normalize, logy,
                   std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                   "PbPb2023 Data", "20 < p_{T} < 30 GeV", "PR SR");
}

// === main macro ===
void overlay_ctau3D_recipes(const char *file, bool savePdf = false, bool normalize = false, bool logy = false)
{
  ctau3D_PR_RSB_all_byPt(file, savePdf, normalize, logy);
  ctau3D_NP_RSB_all_byPt(file, savePdf, normalize, logy);
  ctau3D_NP_RSB_upTo20_byPt(file, savePdf, normalize, logy);
  ctau3D_PR_SR_rap_pt20to30(file, savePdf, normalize, logy);
}
