# Project 7: the Vela pulsar through Chandra's eyes

This is our group's (Group 4) code for **Project 7: Vela Pulsar, X-ray Imaging with Chandra**. We take archival Chandra observations of the Vela pulsar and its wind nebula, clean them up, measure the structures you can see (the torus, the arcs, the jets and the diffuse glow), and look at how the picture changes between soft and hard X-rays.

The short version of what we found: the nebula gets more compact at higher energies, the jets are the hardest (most energetic) features, and the emission gets steadily softer as you move away from the pulsar. That's what you'd expect if particles lose energy to synchrotron radiation as they flow outward. The full story, with all the numbers and figures, is in the report.

📄 **Report:** [`Group4_FinalProjectReport.pdf`](Group4_FinalProjectReport.pdf) (18 pages)

## The data

We used two observations from the Chandra archive:

| ObsID | Date | Instrument | Clean exposure |
|---|---|---|---|
| 128 | 30 Apr 2000 | ACIS-S3, ½ subarray (1.54 s frames) | 10.6 ks |
| 12075 | 4 Sep 2010 | ACIS-S, full frame (3.24 s frames) | 41.0 ks |

ObsID 128 is our main dataset; it's the imaging observation used in the reference paper (Pavlov et al. 2001). Because it's short, we repeated everything on the much deeper ObsID 12075 as a cross-check. The results agree nicely.

The raw data aren't in this repo (they're a few hundred MB), but they're public and easy to fetch. See below.

## What's in here

```
scripts/      the whole analysis, numbered in the order you run it
figures/      every figure that appears in the report
report/       the HTML source of the report and a one-page contact sheet
products/     the measured numbers (JSON) that the report is built from
```

The scripts, in order:

| Script | What it does |
|---|---|
| `setup_ciao.sh` | Activates CIAO and points it at the calibration files |
| `01_reprocess.sh` | Downloads the calibration files, runs `chandra_repro`, checks for background flares (Task A) |
| `02_images.sh` | Makes the soft / medium / hard / broad band images and exposure maps |
| `03_measure.py` | Fits the arcs as tilted rings, measures the jets and the size of the nebula (Task B) |
| `04_energy.py` | Hardness ratios, profiles, the hardness map and colour images (Tasks C and D) |
| `05_responses.sh`, `06_photon_index.py` | Builds the instrument response and turns hardness ratios into photon indices |
| `07_doppler.py` | Uses the front/back brightness of the torus to put a limit on the flow speed |
| `08_outer_jet.py` | Traces the faint, curvy outer jet |
| `10_r50_errors.py`, `12_luminosity.py` | Error bars on the band sizes; X-ray flux and luminosity |
| `09_figures.py`, `11_ciao_rgb.sh` | The final figures, including the standard CIAO three-colour image |
| `13_make_report.py`, `14_render_report.sh` | Builds the report (HTML) and renders it to PDF |

## Running it yourself

You'll need [CIAO](https://cxc.cfa.harvard.edu/ciao/download/conda.html). We installed it with conda, without Sherpa and without the full 2.3 GB calibration database, because our internet was slow:

```bash
conda create -n ciao-4.18 --override-channels \
  -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge \
  ciao ciao-contrib ds9 scipy astropy reportlab numpy=2.3.5
```

Then grab the data and run the pipeline from the repo folder:

```bash
source scripts/setup_ciao.sh
mkdir -p data && (cd data && download_chandra_obsid 128,12075)

for o in 128 12075; do
  bash scripts/01_reprocess.sh $o
  bash scripts/02_images.sh $o
  python scripts/03_measure.py $o
  python scripts/04_energy.py $o
  bash scripts/05_responses.sh $o
  python scripts/06_photon_index.py $o
  python scripts/07_doppler.py $o
  python scripts/08_outer_jet.py $o
  python scripts/10_r50_errors.py $o
  python scripts/12_luminosity.py $o
  bash scripts/11_ciao_rgb.sh $o
done
python scripts/09_figures.py
python scripts/13_make_report.py
```

A few practical notes:

- **Calibration files.** `01_reprocess.sh` only downloads the calibration files each observation needs (about 0.5 GB in total). `download_obsid_caldb` stalled a couple of times for us, so `scripts/fetch_caldb_parallel.sh` does the same job with parallel downloads.
- **ObsID 128 has to be processed before 12075.** Script 03 borrows the PSF width measured in ObsID 128, because the pulsar is piled up in the full-frame observation.
- **Rendering the PDF** uses [WeasyPrint](https://weasyprint.org) in a separate small conda environment, plus the Space Grotesk and IBM Plex Sans fonts:

  ```bash
  conda create -n pdfdesign -c conda-forge weasyprint fonttools
  bash scripts/14_render_report.sh
  ```

- **Spectral fitting.** If you want proper spectral fits later (Sherpa / `specextract`), add Sherpa with:

  ```bash
  conda install -n ciao-4.18 -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge sherpa numpy=2.3.5
  ```

## References

The main paper is Pavlov, Zavlin, Sanwal, Burwitz & Garmire (2001), *The X-ray Spectrum of the Vela Pulsar Resolved with Chandra*, ApJ 552, L129 ([astro-ph/0103171](https://arxiv.org/abs/astro-ph/0103171)). The rest are listed at the end of the report.
