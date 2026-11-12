#!/usr/bin/env bash
# Instrument responses (ARF + RMF) at the pulsar position, used to convert
# the band count ratios (hardness ratios) into power-law photon indices.
# (specextract would do this for an extended region, but it needs Sherpa,
#  which was not installed; the PWN lies within ~1' of the pulsar, where the
#  effective area changes by only a few per cent.)
#   bash scripts/05_responses.sh 128
set -euo pipefail
OBSID=${1:?ObsID}
TOP=$(cd "$(dirname "$0")/.." && pwd)
set +u; source "$TOP/scripts/setup_ciao.sh"; set -u
P=$TOP/products/$OBSID; S=$P/spec; R=$TOP/data/$OBSID/repro
mkdir -p "$S"
read X Y < "$P/pulsar_xy.txt"
EVT=$P/evt2_clean.fits
ASOL=$(ls "$R"/pcadf*asol1.fits | head -1)
MSK=$(ls "$R"/*msk1.fits)
PBK=$(ls "$R"/*pbk0.fits 2>/dev/null | head -1 || true)

punlearn asphist
asphist infile="$ASOL" outfile="$S/asp.hist" evtfile="${EVT}[ccd_id=7]" \
  dtffile="" clobber=yes verbose=0 mode=h

dmcoords "$EVT" asol="$ASOL" option=sky x=$X y=$Y verbose=0 mode=h
CX=$(pget dmcoords chipx); CY=$(pget dmcoords chipy)
echo "ObsID $OBSID: pulsar at S3 chip ($CX, $CY)"

punlearn mkarf
mkarf asphistfile="$S/asp.hist" outfile="$S/psr.arf" sourcepixelx=$X sourcepixely=$Y \
  engrid="0.3:11.0:0.01" obsfile="$EVT" detsubsys=ACIS-S3 grating=NONE \
  maskfile="$MSK" pbkfile="${PBK:-NONE}" dafile=CALDB clobber=yes verbose=0 mode=h

punlearn mkacisrmf
mkacisrmf infile=CALDB outfile="$S/psr.rmf" wmap=none energy="0.3:11.0:0.01" channel="1:1024:1" \
  chantype=PI ccd_id=7 chipx=$(printf "%.0f" $CX) chipy=$(printf "%.0f" $CY) \
  gain=CALDB obsfile="$EVT" clobber=yes verbose=0 mode=h
ls -la "$S"
