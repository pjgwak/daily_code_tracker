#!/bin/bash

for file in csv/*csv; do
   root -l -b -q expandsyst.C+'("'$file'","'${file}2'",true,true,true)'
   mv ${file}2 $file
done
