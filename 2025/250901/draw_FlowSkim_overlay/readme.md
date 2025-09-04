## Overlay Recipe
Draw overlayed plots.

Inputs: histgorams made from "draw_flowSkim" codes 

---

### Usage
Options are
```bash
const char *file, bool savePdf = false, bool normalize = false, bool logy = false
```

E.g.
```cpp
root -b -q 'overlay_mass_recipes.C("/data/users/pjgwak/work/daily_code_tracker/2025/250829/cosmetic_codes_FlowSkim_plots/hists.root", false, false, true)'

root -b -q 'overlay_ctau3D_recipes.C("/data/users/pjgwak/work/daily_code_tracker/2025/250829/cosmetic_codes_FlowSkim_plots/hists.root", false, false, true)'
```