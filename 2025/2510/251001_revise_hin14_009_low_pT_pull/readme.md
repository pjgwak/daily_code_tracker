### How to run

In short,

```bash
 bash run_script.sh
```



### Fit order

1. root -l -b -q reduce_mc_dataset.C
1. root -l -b -q mc_mass_fit.C
1. root -l -b -q mass_fit_raw.C
1. root -l -b -q reduce_data_with_err_cut.C
1. root -l -b -q mass_fit_redData.C
1. root -l -b -q make_ctauErr_pdf.C
1. root -l -b -q ctau_pr_mc_fit.C
1. root -l -b -q ctau_np_mc_fit.C
1. root -l -b -q final_2d_fit.C


### Explnations
1. reduce_mc_dataset.C: apply basic cut to MC dataset
1. mc_mass_fit.C: MC mass fitting
1. mass_fit_raw.C: Data mass fitting with dataset. It's a temporal fit.
1. reduce_data_with_err_cut.C: Decide ctauErr range. It requires mass fit result. Apply ctauErr cut. from now on, we will use this cut dataset (redData) for Data fitting.
1. mass_fit_redData.C: Data mass fit with redData
1. make_err_pdf.C: Make ctau3DErr pdfs as Punzi term
1. ctau_pr_mc_fit.C: To study reslution pdf's shape
1. ctau_np_mc_fit.C: To study signal nomprompt decay shape
1. final_2d_fit.C: Extract b-fraction