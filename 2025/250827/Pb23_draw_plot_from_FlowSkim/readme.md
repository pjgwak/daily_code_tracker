### Usage
ConfigPlot.h: core function and recipes
* PREVIEW_ON: true (use a portion of events. To test macro), false (use all events)
* INPUT_FILE: input FLowSkim path
* To draw other branches
  * Add new array in CFG_1D (TH1D) or CFG_2D (TH2D)
  * Please follow the example arrays inside code to make your array
* You can combine new specil cuts
  * make new basic cuts: refer to CUT_JPSI_MASS
  * combine basic cuts: refer to BuildRapSlices_MassPtU50()
  * Add combined cut to SlicesForTag()

draw_from_config.C: draw very crude plots from TTree and save histograms for cosmetcis -> use the cosmetics codes such as drawH1All.C

It will take about 3 hours (for 300 plots)
```bash
# No parameter
root -b -q draw_from_config.C
```

draw_flowSkim_mass.C: crude macro to draw a plot form TTree (rarely used)