#!/bin/bash

# config file name. No dir path
CONFIG_NAME=pp_run2_pt12p15_y0_1p6.C


# BE CAREFUL HERE
# It directly change the contents in the macros - Be careful.
sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" mc_mass_fit.C

sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" mass_fit.C
sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" sideband_extraction.C
sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_mc_res_fit.C
sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_mc_fit.C
sed -i "s|#include \"cfgs/.*\.C\"|#include \"cfgs/$CONFIG_NAME\"|" ctau_bkg_fit.C


# run macros
# root -l -b -q "mc_mass_fit.C(\"cfgs/$CONFIG_NAME\")"
# root -l -b -q "mass_fit.C(\"cfgs/$CONFIG_NAME\")"
# root -l -b -q "sideband_extraction.C(\"cfgs/$CONFIG_NAME\")" # make RooHistPdf from sideband and signal region
# root -l -b -q "ctau_mc_fit.C(\"cfgs/$CONFIG_NAME\")" # MC prompt ctau fit
root -l -b -q "ctau_bkg_fit.C(\"cfgs/$CONFIG_NAME\")" # ctau data bkg fit



# ---------------------
# root -l -b -q ctau_mc_res_fit.C # MC prompt Res fit - not used