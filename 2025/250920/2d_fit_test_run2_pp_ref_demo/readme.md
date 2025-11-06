## Fit Order
root -b -q mc_mass_center_fit.C

root -b -q mass_center_fit.C

root -b -q 'ctau_side_fit.C("LSB")'
* option1: RSB, LSB

<!-- root -b -q ctau_mc_res_fit.C -->

root -b -q ctau_np_fit.C
