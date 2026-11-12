#!/usr/bin/env bash
# Task D with the standard CIAO recipe ("true colour" image thread):
# smooth each exposure-corrected band image with aconvolve, then combine
# red = 0.5-1.5 keV, green = 1.5-3 keV, blue = 3-8 keV with dmimg2jpg.
#   bash scripts/11_ciao_rgb.sh 12075
set -euo pipefail
OBSID=${1:?ObsID}
TOP=$(cd "$(dirname "$0")/.." && pwd)
set +u; source "$TOP/scripts/setup_ciao.sh"; set -u
I=$TOP/products/$OBSID/img
read X Y < "$TOP/products/$OBSID/pulsar_xy.txt"
# +-64" box around the pulsar (1 px = 0.492")
BOX="[sky=box($X,$Y,260,260,0)]"
for b in soft medium hard; do
  aconvolve infile="$I/band_${b}_flux.img$BOX" outfile="$I/rgb_${b}_sm.img" \
    kernelspec="lib:gaus(2,5,1,2,2)" method=fft clobber=yes verbose=0
done
# asinh scaling, limits set from the data percentiles
lims() { python -c "
import numpy as np; from astropy.io import fits
d=fits.getdata('$1'); d=d[np.isfinite(d)&(d>0)]
print('%g %g'%(np.percentile(d,1),np.percentile(d,99.8)))"; }
read rl rh < <(lims "$I/rgb_soft_sm.img")
read gl gh < <(lims "$I/rgb_medium_sm.img")
read bl bh < <(lims "$I/rgb_hard_sm.img")
dmimg2jpg infile="$I/rgb_soft_sm.img" greenfile="$I/rgb_medium_sm.img" bluefile="$I/rgb_hard_sm.img" \
  outfile="$TOP/figures/F4b_ciao_rgb_${OBSID}.jpg" \
  scalefunction=asinh \
  minred=$rl maxred=$rh mingreen=$gl maxgreen=$gh minblue=$bl maxblue=$bh \
  showaimpoint=no showgrid=no clobber=yes
echo "wrote figures/F4b_ciao_rgb_${OBSID}.jpg"
