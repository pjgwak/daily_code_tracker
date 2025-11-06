#!/bin/bash

for file in $HOME/limits*.o*; do
   bin=`grep 'it=' $file | awk -F '_wSyst_' '{print $2}' | awk -F '.root' '{print $1}'`
   lim=`grep 'The limits are' $file | awk -F ', ' '{print $2}'  | awk -F ']' '{print $1}'`
   echo $bin $lim
done
