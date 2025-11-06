#!/bin/bash

source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh

# acceptance
root -l -b -q 'makeEffs_pbpb_pr.C+(1,true)'

# efficiency (syst uncertainties)
root -l -b -q 'makeEffs_pbpb_pr.C+(1)'
root -l -b -q 'makeEffs_pbpb_pr.C+(10)'
root -l -b -q 'makeEffs_pbpb_pr.C+(11)'
root -l -b -q 'makeEffs_pbpb_pr.C+(12)'
root -l -b -q 'makeEffs_pbpb_pr.C+(13)'
root -l -b -q 'makeEffs_pbpb_pr.C+(14)'
root -l -b -q 'makeEffs_pbpb_pr.C+(15)'
root -l -b -q 'makeEffs_pbpb_pr.C+(16)'
root -l -b -q 'makeEffs_pbpb_pr.C+(17)'
root -l -b -q 'makeEffs_pbpb_pr.C+(18)'

# efficiency (stat uncertainties)
root -l -b -q 'makeEffs_pbpb_pr.C+(5)'
root -l -b -q 'makeEffs_pbpb_pr.C+(6)'
root -l -b -q 'makeEffs_pbpb_pr.C+(7)'
root -l -b -q 'makeEffs_pbpb_pr.C+(8)'
