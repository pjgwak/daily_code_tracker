### Usage
PlotCore.C: It has the core functions to draw the plots. Don't need to touch.

drawH1All.C: draw all TH1D plots
```bash
# Parameters: root file path, TDirectory name in the root file, savePDF
root -b -q 'drawH1All.C("/data/users/pjgwak/work/daily_code_tracker/2025/250827/Pb23_draw_plot_from_FlowSkim/figs/hists.root", "h1", false)'
```

drawh1Filter.C: draw all TH1D plots passing name filter. Filter: all patterns in "include" and no pattern in "exclude".
```bash
# Parameters: root file path, Inclusive pattern, Exclusive pattern, TDirectory name in the root file, savePDF, isDryRun, output directory
root -l -b -q 'drawH1Filter.C("/data/users/pjgwak/work/daily_code_tracker/2025/250827/Pb23_draw_plot_from_FlowSkim/figs/hists.root", "h1", "mass_,PR,SR", "LSB,RSB", false, false, "figs")'
```

drawH2All.C: draw all TH2D plots
```bash
# Parameters: root file path, TDirectory name in the root file, savePDF, isLogZ, output folder
root -l -b -q 'drawH2All.C("/data/users/pjgwak/work/daily_code_tracker/2025/250827/Pb23_draw_plot_from_FlowSkim/figs/hists.root", "h2", false, false, "figs")'
```

drawOne.C: draw one target histogram -> default output folder: figs_test
```bash
# Parameters: root file path, histogram path including a directory, savePDF
root -b -q 'drawOne.C("/data/users/pjgwak/work/daily_code_tracker/2025/250827/Pb23_draw_plot_from_FlowSkim/figs/hists.root", "h1/ctau3D_region6_NP_RSB_all_pt20to30", false)'
```

drawHistBasic.C: very crude macro to draw one histogram for quick check (rarely used)