#!/bin/bash

# ---------------------
root -l -b -q reduce_mc_dataset.C
root -l -b -q mc_mass_fit.C
root -l -b -q mass_fit_raw.C
root -l -b -q reduce_data_with_err_cut.C
root -l -b -q mass_fit_redData.C
root -l -b -q make_ctauErr_pdf.C
root -l -b -q ctau_pr_mc_fit.C
root -l -b -q ctau_np_mc_fit.C
root -l -b -q ctau_bkg_fit.C
root -l -b -q final_2d_fit.C





#-----------------------------
# # config file name. No dir path
# CONFIG_NAME=pp_run2_pt12p15_y0_1p6.C


# # BE CAREFUL HERE
# # It directly change the contents in the macros - Be careful.
# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" mc_mass_fit.C

# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" mass_fit.C
# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" sideband_extraction.C
# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_mc_res_fit.C
# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_mc_fit.C
# sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_bkg_fit.C


# # run macros
# # root -l -b -q "mc_mass_fit.C(\"cfgs/$CONFIG_NAME\")"
# # root -l -b -q "mass_fit.C(\"cfgs/$CONFIG_NAME\")"
# # root -l -b -q "sideband_extraction.C(\"cfgs/$CONFIG_NAME\")" # make RooHistPdf from sideband and signal region
# # root -l -b -q "ctau_mc_fit.C(\"cfgs/$CONFIG_NAME\")" # MC prompt ctau fit
# root -l -b -q "ctau_bkg_fit.C(\"cfgs/$CONFIG_NAME\")" # ctau data bkg fit
