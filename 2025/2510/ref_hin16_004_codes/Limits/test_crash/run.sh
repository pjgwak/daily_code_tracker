#!/bin/bash

root -l -b -q gen1.C'("uncut_PbPb.root","PbPb")'
sleep 1 
root -l -b -q gen1.C'("uncut_PP.root","PP")'


root -l -b -q gen2.C'("uncut_PbPb.root","cut_PbPb.root","PbPb")'
root -l -b -q gen2.C'("uncut_PP.root","cut_PP.root","PP")'

root -l -b -q combine.C'("cut_PbPb.root","cut_PP.root")'
