#!/usr/bin/env bash
# Project 7 - Vela pulsar with Chandra
# Tasks A and C: broad-band and energy-resolved images.
#
#   bash scripts/02_images.sh 128
#
# Products (in products/<obsid>/img/):
#   band_{soft,medium,hard,broad}_thresh.img   counts images (1 sky px = 0.492")
#   band_*_thresh.expmap / band_*_flux.img      exposure maps / exposure-corrected
#   fine_*.img                                   counts images binned at 0.5 px
#                                                (0.246") of the inner nebula
set -euo pipefail
OBSID=${1:?ObsID}
TOP=$(cd "$(dirname "$0")/.." && pwd)
set +u; source "$TOP/scripts/setup_ciao.sh"; set -u

P=$TOP/products/$OBSID
REPRO=$TOP/data/$OBSID/repro
EVT=$P/evt2_clean.fits
mkdir -p "$P/img"
read X Y < "$P/pulsar_xy.txt"

# Image grid centred on the pulsar: +-256 sky pixels (+-2.1') at 1 px = 0.492"
x0=$(python -c "print(round($X)-256.5)"); x1=$(python -c "print(round($X)+255.5)")
y0=$(python -c "print(round($Y)-256.5)"); y1=$(python -c "print(round($Y)+255.5)")
GRID="$x0:$x1:1,$y0:$y1:1"
echo "grid: $GRID"

ASOL=$(ls "$REPRO"/pcadf*_asol1.fits | tr '\n' ',' | sed 's/,$//')
BPIX=$(ls "$REPRO"/acisf*_repro_bpix1.fits)
MSK=$(ls "$REPRO"/acisf*_msk1.fits)

# Bands: soft 0.5-1.5, medium 1.5-3, hard 3-8, broad 0.5-8 keV
# (third number = monoenergy in keV used to weight the exposure map)
punlearn fluximage
fluximage infile="$EVT" outroot="$P/img/band" \
  bands="0.5:1.5:1.0,1.5:3:2.2,3:8:4.5,0.5:8:2.3" \
  xygrid="$GRID" asolfile="$ASOL" badpixfile="$BPIX" maskfile="$MSK" \
  units=time expmapthresh=1% parallel=yes clobber=yes verbose=1

# fluximage names outputs by band limits; give them readable names
cd "$P/img"
for pair in "0.5-1.5:soft" "1.5-3:medium" "3-8:hard" "0.5-8:broad"; do
  b=${pair%%:*}; n=${pair##*:}
  for ext in thresh.img thresh.expmap flux.img; do
    [ -f "band_${b}_${ext}" ] && mv -f "band_${b}_${ext}" "band_${n}_${ext}"
  done
done

# Finer (0.5 sky px = 0.246") counts images of the inner +-64" for morphology;
# chandra_repro's default EDSER sub-pixel event repositioning makes this useful.
fx0=$(python -c "print(round($X)-130.5)"); fx1=$(python -c "print(round($X)+129.5)")
fy0=$(python -c "print(round($Y)-130.5)"); fy1=$(python -c "print(round($Y)+129.5)")
for pair in "500:1500:soft" "1500:3000:medium" "3000:8000:hard" "500:8000:broad"; do
  IFS=: read lo hi n <<< "$pair"
  dmcopy "$EVT[energy=$lo:$hi][bin x=$fx0:$fx1:0.5,y=$fy0:$fy1:0.5]" "fine_$n.img" clobber=yes
done

# Count totals per band for the log
for n in soft medium hard broad; do
  printf "%-7s counts in image: %s\n" $n "$(dmstat band_${n}_thresh.img centroid- verbose=0 >/dev/null; pget dmstat out_sum)"
done
ls -la
