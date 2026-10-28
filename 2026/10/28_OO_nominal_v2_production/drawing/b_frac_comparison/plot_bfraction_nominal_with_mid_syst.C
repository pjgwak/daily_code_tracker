// From this folder: root -n -l -b -q 'plot_bfraction_nominal_with_mid_syst.C()'
// Mid nominal central values/statistical errors; syst = abs(B_mid - B_nominal).
// Reuses the forward layout and rootlogon; reads saved fits only.
#include "plot_bfraction_nominal_with_fwd_syst.C"

void plot_bfraction_nominal_with_mid_syst(const char *rootDir = "", const char *outputDir = "")
{
  plot_bfraction_nominal_with_syst(rootDir, outputDir, true);
}
