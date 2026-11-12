"""Poisson (parametric bootstrap) uncertainties on the per-band half-light radii.
  python scripts/10_r50_errors.py 128
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
OBSID = sys.argv[1]; P = products(OBSID)
wide = {b: load_image(f"{P}/img/band_{b}_thresh.img")[0] for b in BANDS}
wexp = {b: load_image(f"{P}/img/band_{b}_thresh.expmap")[0] for b in BANDS}
_, wh = load_image(f"{P}/img/band_broad_thresh.img")
WE, WN = offsets_arcsec(wh, wide["soft"].shape, OBSID); WR = np.hypot(WE, WN)
cov = load_image(f"{P}/img/band_broad_thresh.expmap")[0]
cov = cov > 0.5 * cov.max()
outer = (WR > 90) & (WR < 120) & cov
rr = np.arange(2, 101, 0.5)
ring_idx = np.digitize(WR, rr)
def r50(img, b):
    rate = np.where(cov, img / np.maximum(wexp[b], 1), 0)
    net = np.where(cov, rate - rate[outer].mean(), 0)
    sums = np.bincount(ring_idx.ravel(), weights=net.ravel(), minlength=len(rr) + 1)[1:len(rr)]
    cum = np.cumsum(sums)
    return float(np.interp(0.5 * cum[-1], cum, rr[1:]))
rng = np.random.default_rng(3)
ER = json.load(open(f"{P}/energy_results.json"))
out = {}
sims = {b: [r50(rng.poisson(wide[b]).astype(float), b) for _ in range(200)] for b in BANDS}
for b in BANDS:
    out[b] = dict(r50=r50(wide[b], b), err=float(np.std(sims[b])))
d = np.array(sims["soft"]) - np.array(sims["hard"])
out["soft_minus_hard"] = dict(value=out["soft"]["r50"] - out["hard"]["r50"], err=float(np.std(d)))
ER["band_r50_bootstrap"] = out
json.dump(ER, open(f"{P}/energy_results.json", "w"), indent=1, default=float)
print(OBSID, {k: (round(v["r50"] if "r50" in v else v["value"], 2), round(v["err"], 2)) for k, v in out.items()})
