#!/usr/bin/env bash
# Fetch the CALDB files an ObsID needs, several at a time.
# (download_obsid_caldb fetches one file at a time and can stall on a slow
#  connection; this uses it only to list what is missing, then downloads the
#  files in parallel with curl, resuming partial files.)
#   bash scripts/fetch_caldb_parallel.sh data/12075/primary/acisf12075N004_evt2.fits.gz
set -uo pipefail
EVT=${1:?event file}
TOP=$(cd "$(dirname "$0")/.." && pwd)
set +u; source "$TOP/scripts/setup_ciao.sh"; set -u
B=https://cxc.cfa.harvard.edu/cdaftp/arcftp/ChandraCalDB
TMP=$(mktemp -d)

# index files first (cheap), via the official tool with missing=yes
missing=$(download_obsid_caldb infile="$EVT" outdir="$CALDB" missing=yes verbose=2 2>&1 |
          awk '/MISSING/{print $1}' | sed 's/\.\.\.*$//')
echo "missing: $(echo $missing | wc -w) files"

fetch() {
  f=$1
  rel=""
  for idx in "$CALDB"/data/chandra/*/caldb.indx; do
    rel=$(dmlist "${idx}[cols cal_dir,cal_file]" data,clean 2>/dev/null |
          awk -v f="$f" '$2==f{print $1"/"$2; exit}')
    [ -n "$rel" ] && break
  done
  [ -z "$rel" ] && { echo "NOT IN INDEX $f"; return; }
  mkdir -p "$CALDB/$(dirname "$rel")"
  for try in 1 2 3 4 5; do
    curl -sf -C - --speed-limit 1000 --speed-time 60 -o "$TMP/$f" "$B/$rel" && break
  done
  if [ -s "$TMP/$f" ]; then mv "$TMP/$f" "$CALDB/$rel"; echo "got $rel"; else echo "FAIL $rel"; fi
}

for f in $missing; do fetch "$f" & done
wait
download_obsid_caldb infile="$EVT" outdir="$CALDB" missing=yes 2>&1 | grep -E 'MISSING' || echo "all CALDB files present"
rm -rf "$TMP"
