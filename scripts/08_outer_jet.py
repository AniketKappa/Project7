"""Trace the faint, curved NW outer jet (clearly visible in the deep ObsID 12075).

  python scripts/08_outer_jet.py 12075

The jet is followed from a hand-placed guide polyline (read off the
smoothed 0.5-8 keV image); every guide point is moved to the brightness
ridge within +-5" perpendicular to the path, and each 4"-long segment is
tested: counts in a 3"-wide strip on the ridge versus two 3"-wide strips
6-9" either side.  Segments with S/N >= 3 are kept; the jet length is the
path length between the first and last significant segments.
"""
import sys, os, json
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

OBSID = sys.argv[1] if len(sys.argv) > 1 else "12075"
P = products(OBSID)
GUIDE = {  # (east, north) arcsec offsets from the pulsar, read off the image
    "12075": [(-28, 20), (-36, 27), (-44, 34), (-50, 39), (-58, 43), (-66, 45), (-74, 47),
              (-82, 50), (-88, 51), (-93, 47), (-96, 42)],
    "128": [(-44, 52), (-48, 57), (-51, 61), (-53, 64)],
}[OBSID]

img, hdr = load_image(f"{P}/img/band_broad_thresh.img")
E, N = offsets_arcsec(hdr, img.shape, OBSID)
dEdj = E[0, 1] - E[0, 0]; dNdi = N[1, 0] - N[0, 0]
sm = smooth(img, 1.5 / PIX)


def samp(im, e, n):
    return map_coordinates(im, [np.atleast_1d((np.asarray(n) - N[0, 0]) / dNdi),
                                np.atleast_1d((np.asarray(e) - E[0, 0]) / dEdj)], order=1)


# resample the guide every 2"
g = np.array(GUIDE, float)
seg = np.hypot(*np.diff(g, axis=0).T); s = np.r_[0, np.cumsum(seg)]
ss = np.arange(0, s[-1] + 1e-6, 2.0)
pts = np.c_[np.interp(ss, s, g[:, 0]), np.interp(ss, s, g[:, 1])]
# move every point to the ridge
ridge = []
for k, p in enumerate(pts):
    t = pts[min(k + 1, len(pts) - 1)] - pts[max(k - 1, 0)]; t /= np.hypot(*t)
    nrm = np.array([-t[1], t[0]])
    off = np.arange(-5, 5.01, 0.25)
    v = samp(sm, p[0] + off * nrm[0], p[1] + off * nrm[1])
    ridge.append(p + off[np.argmax(v)] * nrm)
ridge = np.array(ridge)

# significance of 4" segments (counts in strips of the unsmoothed image)
res_seg = []
for k in range(0, len(ridge) - 2, 2):
    a, b = ridge[k], ridge[k + 2]
    t = (b - a) / np.hypot(*(b - a)); nrm = np.array([-t[1], t[0]])
    rel_e, rel_n = E - a[0], N - a[1]
    along = rel_e * t[0] + rel_n * t[1]; perp = rel_e * nrm[0] + rel_n * nrm[1]
    inseg = (along >= 0) & (along < np.hypot(*(b - a)))
    on = inseg & (np.abs(perp) < 1.5)
    off = inseg & (np.abs(perp) >= 6) & (np.abs(perp) < 9)
    c_on, c_off = img[on].sum(), img[off].sum()
    scale = on.sum() / max(off.sum(), 1)
    net = c_on - scale * c_off
    err = np.sqrt(c_on + scale**2 * c_off)
    res_seg.append(dict(start=a.tolist(), end=b.tolist(), on=float(c_on), net=float(net),
                        snr=float(net / err), dist=float(np.hypot(*(0.5 * (a + b))))))
sig = [r for r in res_seg if r["snr"] >= 3]
fwhm = np.nan
if sig:
    i0 = res_seg.index(sig[0]); i1 = res_seg.index(sig[-1])
    path = np.array([res_seg[i]["start"] for i in range(i0, i1 + 1)] + [res_seg[i1]["end"]])
    L = float(np.sum(np.hypot(*np.diff(path, axis=0).T)))
    r_start, r_end = float(np.hypot(*path[0])), float(np.hypot(*path[-1]))
    pa_end = float(np.degrees(np.arctan2(path[-1][0], path[-1][1])) % 360)
    # width from the averaged perpendicular profile of the (unsmoothed) image
    sm_w = smooth(img, 0.5 / PIX)
    off = np.arange(-8, 8.01, 0.5)
    profs = []
    for k in range(len(path) - 1):
        a, b = path[k], path[k + 1]
        t = (b - a) / np.hypot(*(b - a)); nrm = np.array([-t[1], t[0]])
        for f in np.linspace(0, 1, 5):
            c = a + f * (b - a)
            profs.append(samp(sm_w, c[0] + off * nrm[0], c[1] + off * nrm[1]))
    prof = np.mean(profs, axis=0)
    fg = lambda x, A, x0, s, c: A * np.exp(-(x - x0)**2 / 2 / s**2) + c
    try:
        pp, _ = curve_fit(fg, off, prof, p0=[prof.max() - prof.min(), 0, 2, prof.min()])
        fwhm = float(np.sqrt(max((2.3548 * abs(pp[2]))**2 - (2.3548 * 0.5)**2, 0)))
    except Exception:
        pass
else:
    path = np.empty((0, 2)); L = 0.0; r_start = r_end = pa_end = np.nan
out = dict(obsid=OBSID, ridge=ridge.tolist(), segments=res_seg, path=path.tolist(),
           path_length_arcsec=L, start_distance_arcsec=r_start, end_distance_arcsec=r_end,
           end_pa_deg=pa_end, width_fwhm_arcsec=fwhm, n_sig_segments=len(sig), n_segments=len(res_seg))
json.dump(out, open(f"{P}/outer_jet.json", "w"), indent=1)
print(f"ObsID {OBSID}: {len(sig)}/{len(res_seg)} segments with S/N>=3; S/N per segment: "
      + " ".join(f"{r['snr']:.1f}" for r in res_seg))
print(f"  path length {L:.0f}\" ({arcsec_to_pc(L):.3f} pc), from r={r_start:.0f}\" to r={r_end:.0f}\" "
      f"({arcsec_to_pc(r_end):.3f} pc) at PA {pa_end:.0f} deg; width (FWHM) ~{fwhm:.1f}\"")
