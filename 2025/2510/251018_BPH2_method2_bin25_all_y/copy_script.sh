#!/usr/bin/env bash
TEMPLATE="ctau_fit_pT6p5_7p5.C"

PT=(6.5 7.5 8.5 9.5 11.0 13.0 15.0 17.5 20.0 25.0 30.0 35.0 50.0)

# 6.5 -> 6p5, 11.0 -> 11p0
toP(){ echo "$1" | sed 's/\./p/g'; }

mkdir -p copies
for ((i=0; i<${#PT[@]}-1; ++i)); do
  lo=${PT[i]} ; hi=${PT[i+1]}
  out="copies/ctau_fit_pT$(toP "$lo")_$(toP "$hi").C"
  cp -n "$TEMPLATE" "$out"

  sed -i -E "s/(ptLow\s*=\s*)[0-9.]+/\1$lo/; s/(ptHigh\s*=\s*)[0-9.]+/\1$hi/" "$out"

  echo "made $out"
done
