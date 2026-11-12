#!/usr/bin/env bash
# Project 7 - Vela pulsar with Chandra
# Task A: reprocess the archival data and produce a clean level-2 event file.
#
# Usage (inside the ciao-4.18 conda environment):
#     bash scripts/01_reprocess.sh 128
#
# Steps
#   1. download_obsid_caldb : fetch only the calibration files this ObsID needs
#   2. chandra_repro        : reprocess with the latest calibration -> new evt2
#   3. acisreadcorr         : flag the readout streak of the bright pulsar
#   4. light curve + deflare: remove background-flare intervals (GTI filter)
set -euo pipefail

OBSID=${1:?give an ObsID, e.g. 128}
TOP=$(cd "$(dirname "$0")/.." && pwd)
set +u; source "$TOP/scripts/setup_ciao.sh"; set -u
DATA=$TOP/data
OBS=$DATA/$OBSID
REPRO=$OBS/repro
OUT=$TOP/products/$OBSID
mkdir -p "$OUT"

# Vela pulsar (PSR B0833-45), ICRS J2000 position from VLBI (Dodson et al. 2003)
PSR_RA=128.8358       # 08h35m20.61s
PSR_DEC=-45.17635     # -45d10m34.9s

# ---------------------------------------------------------------- CALDB
# Only the files needed for this observation (instead of the 2.3 GB full CALDB).
# They go to $ASCDS_INSTALL/CALDB, which the ciao-4.18 environment picks up
# automatically when it is activated.
EVT2_ARCHIVE=$(ls "$OBS"/primary/acisf*_evt2.fits*)
download_obsid_caldb infile="$EVT2_ARCHIVE" outdir='${ASCDS_INSTALL}/CALDB' clobber=no
if [ -z "${CALDB:-}" ]; then
  export CALDB=$ASCDS_INSTALL/CALDB
  export CALDBCONFIG=$CALDB/software/tools/caldb.config
  export CALDBALIAS=$CALDB/software/tools/alias_config.fits
fi
echo "CALDB: $CALDB"

# ---------------------------------------------------------------- reprocess
punlearn chandra_repro
chandra_repro indir="$OBS" outdir="$REPRO" cleanup=yes clobber=yes verbose=1

EVT=$(ls "$REPRO"/acisf*_repro_evt2.fits)
BPIX=$(ls "$REPRO"/acisf*_repro_bpix1.fits)
ASOL=$(ls "$REPRO"/pcadf*_asol1.fits | tr '\n' ',' | sed 's/,$//')
MSK=$(ls "$REPRO"/acisf*_msk1.fits 2>/dev/null || ls "$OBS"/secondary/acisf*_msk1.fits*)

# ---------------------------------------------------------------- readout streak
# A bright source leaves a faint readout ("out-of-time") streak along the CCD
# columns. For ObsID 128 the fraction is only ~ readout time / frame time
# = (512 rows x 40 us) / 1.54 s ~ 1.3% of the ~2850 pulsar counts, i.e. ~40
# counts spread over the whole column, which is negligible next to the PWN.
# acisreadcorr without a streak spectrum removes *all* events in the strip
# (including real nebula photons) and leaves a dark stripe, so we do not apply
# it; the streak is simply noted in the analysis.
dmcoords "$EVT" asol="$ASOL" option=cel ra=$PSR_RA dec=$PSR_DEC celfmt=deg verbose=0
X=$(pget dmcoords x); Y=$(pget dmcoords y)
echo "Pulsar sky position: x=$X y=$Y"
cp "$EVT" "$OUT/evt2_repro.fits"   # (uncorrected; name kept for the steps below)

# ---------------------------------------------------------------- flares
# Light curve of a source-free region on the same chip (0.5-7 keV), 3.2*80 s bins
CHIP=$(dmkeypar "$EVT" DETNAM echo+)
echo "Detector: $CHIP"
cat > "$OUT/bkg_lc.reg" <<EOF
annulus($X,$Y,250,330)
EOF
punlearn dmextract
dmextract infile="$OUT/evt2_repro.fits[energy=500:7000][sky=region($OUT/bkg_lc.reg)][bin time=::259.28]" \
  outfile="$OUT/bkg_lc.fits" opt=ltc1 clobber=yes
punlearn deflare
deflare infile="$OUT/bkg_lc.fits" outfile="$OUT/bkg_lc.gti" method=clean \
  plot=no save="$OUT/bkg_lc" verbose=1 || true

dmcopy "$OUT/evt2_repro.fits[@$OUT/bkg_lc.gti][status=0]" "$OUT/evt2_clean.fits" clobber=yes

echo "Exposure before/after flare filtering:"
dmkeypar "$EVT" EXPOSURE echo+
dmkeypar "$OUT/evt2_clean.fits" EXPOSURE echo+
echo "Clean event file: $OUT/evt2_clean.fits"
echo "$X $Y" > "$OUT/pulsar_xy.txt"
