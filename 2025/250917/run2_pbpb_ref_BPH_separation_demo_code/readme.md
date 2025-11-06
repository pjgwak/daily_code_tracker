## Fit Order
root -b -q mc_mass_center_fit.C

root -b -q 'mass_side_fit.C("PR", "RSB")'
* option1: PR, NP
* option2: RSB, LSB

root -b -q 'mass_center_fit.C("PR")'
* option1: PR, NP

root -b -q 'ctau_side_fit.C("RSB")'
* option1: RSB, LSB

root -b -q ctau_mc_res_fit.C

root -b -q ctau_pr_fit.C

root -b -q get_pr_np_numbers.C