"""Task B: measure the sizes of the main structures of the Vela PWN.

  python scripts/03_measure.py 128

Method
  * Image: 0.5-8 keV counts, 0.246" bins, Gaussian-smoothed.
  * Arcs: for each position angle (PA, N through E) find the radius of the
    brightness ridge inside a search window; fit an ellipse to the ridge
    points (least-squares conic fit).  Semi-axes -> size; axis ratio ->
    inclination of a circular ring, cos i = b/a.
  * Jets: brightness profile along the symmetry axis in a strip, and
    transverse profiles fitted with a Gaussian + linear background; the
    PSF width (measured on the pulsar) is removed in quadrature.
  * Bright PWN: elliptical extent of the region whose smoothed surface
    brightness exceeds the outer (r = 90-120") level by a factor of 2 and 4,
    plus a half-light radius of the nebula (pulsar excluded).
All angular sizes are converted to physical ones at d = 287 pc.
"""
import sys, json, os
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from astropy.io import fits

OBSID = sys.argv[1] if len(sys.argv) > 1 else "128"
P = products(OBSID)
OUT = {}

img, hdr = load_image(f"{P}/img/fine_broad.img")
E, N = offsets_arcsec(hdr, img.shape, OBSID)
dEdj = E[0, 1] - E[0, 0]; dNdi = N[1, 0] - N[0, 0]
BIN = abs(dNdi)                                   # arcsec per image pixel


def sample(im, e, n):
    """Bilinear sample of image `im` at offsets (e, n) arcsec."""
    j = (np.asarray(e) - E[0, 0]) / dEdj
    i = (np.asarray(n) - N[0, 0]) / dNdi
    return map_coordinates(im, [i, j], order=1, mode="nearest")


def pa_vec(pa_deg):
    t = np.radians(pa_deg)
    return np.sin(t), np.cos(t)                  # (east, north) unit vector


# ------------------------------------------------------------------ PSF
# Radial profile of the pulsar in the hard band (least nebula contamination).
# Fit a Gaussian core; FWHM is used as the effective resolution.
r_ev = np.hypot(E, N)
hard, _ = load_image(f"{P}/img/fine_hard.img")
rb = np.arange(0, 4.01, BIN)
prof = np.array([hard[(r_ev >= a) & (r_ev < a + BIN)].mean() for a in rb[:-1]])
rc = rb[:-1] + BIN / 2
g = lambda r, A, s, c: A * np.exp(-r**2 / (2 * s**2)) + c
(pA, ps, pc), _ = curve_fit(g, rc, prof, p0=[prof[0], 0.5, prof[-1]])
PSF_FWHM = 2.3548 * abs(ps)
OUT["psf_fwhm_pulsar_arcsec"] = PSF_FWHM
# In a full-frame (3.2 s) observation the pulsar core is piled up, which
# flattens and broadens it; use the value measured in the 1.54-s subarray
# observation ObsID 128 (0.73") as the resolution in that case.
if float(fits.getval(f"{P}/evt2_clean.fits", "TIMEDEL", 1)) > 2.0:
    PSF_FWHM = json.load(open(f"{products(128)}/measurements.json"))["psf_fwhm_arcsec"]
OUT["psf_fwhm_arcsec"] = PSF_FWHM


# ------------------------------------------------------------------ arcs
def ridge(sm, pa_list, rmin, rmax, dr=0.1):
    pts = []
    rr = np.arange(rmin, rmax, dr)
    for pa in pa_list:
        ue, un = pa_vec(pa)
        prof = sample(sm, rr * ue, rr * un)
        k = np.argmax(prof)
        # keep only genuine interior maxima clearly above the window edges
        if 2 < k < len(rr) - 3 and prof[k] > 1.10 * min(prof[0], prof[-1]):
            pts.append((pa, rr[k], prof[k]))
    return np.array(pts)


def _unused_fit_ellipse(e, n):
    """Direct least-squares ellipse fit (Fitzgibbon et al. 1999).
    Returns centre (e0, n0), semi-axes a >= b, PA of the major axis."""
    D = np.vstack([e * e, e * n, n * n, e, n, np.ones_like(e)]).T
    S = D.T @ D
    C = np.zeros((6, 6)); C[0, 2] = C[2, 0] = 2; C[1, 1] = -1
    w, v = np.linalg.eig(np.linalg.solve(S, C))
    k = np.argmax(np.real(w))
    A, B, Cc, Dd, Ee, F = np.real(v[:, k])
    den = B * B - 4 * A * Cc
    e0 = (2 * Cc * Dd - B * Ee) / den
    n0 = (2 * A * Ee - B * Dd) / den
    num = 2 * (A * Ee**2 + Cc * Dd**2 - B * Dd * Ee + den * F)
    s = np.sqrt((A - Cc)**2 + B**2)
    a1 = -np.sqrt(num * (A + Cc + s)) / den
    a2 = -np.sqrt(num * (A + Cc - s)) / den
    # angle of the axis a1 measured from the +e axis toward +n
    th = 0.5 * np.arctan2(-B, Cc - A)
    a, b = (a1, a2) if a1 >= a2 else (a2, a1)
    if a1 < a2:
        th += np.pi / 2
    # convert to PA (from north through east) of the major axis
    pa = (np.degrees(np.arctan2(np.cos(th), np.sin(th)))) % 180
    return dict(e0=e0, n0=n0, a=a, b=b, pa_major=pa)


def bootstrap_ellipse(pts, nboot=300, seed=1):
    rng = np.random.default_rng(seed)
    e = pts[:, 1] * np.sin(np.radians(pts[:, 0]))
    n = pts[:, 1] * np.cos(np.radians(pts[:, 0]))
    best = fit_ellipse(e, n)
    sims = []
    for _ in range(nboot):
        k = rng.integers(0, len(e), len(e))
        try:
            sims.append(fit_ellipse(e[k], n[k]))
        except Exception:
            pass
    err = {key: float(np.nanstd([s[key] for s in sims])) for key in best}
    return best, err


def ring_r(pa_deg, a, b, d, axis):
    """Radius (from the pulsar) at position angle `pa_deg` of the projected
    front half of a circular ring: an ellipse with semi-axes a (perpendicular
    to the symmetry axis) and b (along it), centred a distance d from the
    pulsar along the axis (PA `axis`)."""
    t = np.radians(pa_deg); ta = np.radians(axis)
    de, dn = np.sin(t), np.cos(t)
    ue, un = np.sin(ta), np.cos(ta)              # along axis
    q = de * ue + dn * un                        # component along axis
    p = de * un - dn * ue                        # perpendicular component
    A = p**2 / a**2 + q**2 / b**2
    B = -2 * q * d / b**2
    C = d**2 / b**2 - 1
    disc = np.maximum(B**2 - 4 * A * C, 0)
    return (-B + np.sqrt(disc)) / (2 * A)


def fit_ring(pts, p0):
    pa, r = pts[:, 0], pts[:, 1]
    lo = [1, 0.5, -30, p0[3] - 30]; hi = [60, 60, 30, p0[3] + 30]
    p, cov = curve_fit(ring_r, pa, r, p0=p0, bounds=(lo, hi), maxfev=20000)
    return p


def bootstrap_ring(pts, p0, nboot=300, seed=1):
    rng = np.random.default_rng(seed)
    best = fit_ring(pts, p0)
    sims = []
    for _ in range(nboot):
        k = rng.integers(0, len(pts), len(pts))
        try:
            sims.append(fit_ring(pts[k], best))
        except Exception:
            pass
    sims = np.array(sims)
    return best, sims.std(axis=0), sims


def ridge_guided(sm, pa_list, guess, half=3.0, dr=0.1, rmin=1.5):
    """Ridge radius at each PA, searched within +-half of a guide curve."""
    pts = []
    for pa in pa_list:
        g0 = guess(pa)
        rr = np.arange(max(rmin, g0 - half), g0 + half, dr)
        ue, un = pa_vec(pa)
        prof = sample(sm, rr * ue, rr * un)
        k = np.argmax(prof)
        if 2 < k < len(rr) - 3:
            pts.append((pa, rr[k], prof[k]))
    return np.array(pts)


def wrap(pa, axis):
    return (pa - axis + 180) % 360 - 180      # PA relative to axis in (-180,180]


AXIS0 = 310.0                                  # first guess of the symmetry axis
sm = smooth(img, 1.0 / BIN)                    # sigma = 1.0"
# Guide curves from the unwrapped (r, PA) image: the outer arc is ~14" from the
# pulsar along the axis and ~30" at 90 deg from it; the inner arc is a fainter,
# flatter arc ~7" from the pulsar along the axis and ~10" at the sides.
# (The bright blob 1.5-3.5" NW of the pulsar is the compact inner knot / jet
#  base, not the inner arc.)
g_out = lambda pa: 14 + 16 * (abs(wrap(pa, AXIS0)) / 90) ** 1.6
g_in = lambda pa: 7.0 + 3.0 * (abs(wrap(pa, AXIS0)) / 90) ** 1.5
pa_out = np.arange(AXIS0 - 75, AXIS0 + 76, 2) % 360
pa_in = np.arange(AXIS0 - 85, AXIS0 + 86, 2) % 360
# skip the jet direction itself for the inner arc (jet dominates there)
pa_in = pa_in[np.abs(wrap(pa_in, AXIS0)) > 12]
outer = ridge_guided(sm, pa_out, g_out, half=5.0)
inner = ridge_guided(smooth(img, 0.75 / BIN), pa_in, g_in, half=2.0)

for name, pts, p0 in [("outer_arc", outer, [30, 14, 0, AXIS0]),
                      ("inner_arc", inner, [10, 4, 3, AXIS0])]:
    best, err, sims = bootstrap_ring(pts, p0)
    a, b, d, ax = best
    cosi = np.clip(b / a, 0, 1); incl = np.degrees(np.arccos(cosi))
    incl_s = np.degrees(np.arccos(np.clip(sims[:, 1] / sims[:, 0], 0, 1)))
    res = pts[:, 1] - ring_r(pts[:, 0], *best)
    OUT[name] = dict(a=float(a), b=float(b), d=float(d), axis_pa=float(ax % 360),
                     err=dict(a=float(err[0]), b=float(err[1]), d=float(err[2]), axis_pa=float(err[3])),
                     inclination_deg=float(incl), inclination_err=float(np.std(incl_s)),
                     vertex_dist=float(d + b), rms_resid=float(np.sqrt(np.mean(res**2))),
                     npts=int(len(pts)), ridge=pts.tolist())

AXIS_PA = OUT["outer_arc"]["axis_pa"]
OUT["axis_pa_deg"] = float(AXIS_PA)

# Inner arc with the outer ring's inclination and axis imposed (coplanar /
# parallel rings, as in the double-torus picture); only a and d are free.
ci = OUT["outer_arc"]["b"] / OUT["outer_arc"]["a"]
f_in = lambda pa, a, d: ring_r(pa, a, a * ci, d, AXIS_PA)
rng = np.random.default_rng(2)
pc, _ = curve_fit(f_in, inner[:, 0], inner[:, 1], p0=[10, 3])
sims = []
for _ in range(300):
    k = rng.integers(0, len(inner), len(inner))
    try:
        sims.append(curve_fit(f_in, inner[k, 0], inner[k, 1], p0=pc)[0])
    except Exception:
        pass
sims = np.array(sims)
OUT["inner_arc_fixed_incl"] = dict(a=float(pc[0]), b=float(pc[0] * ci), d=float(pc[1]),
                                   err=dict(a=float(sims[:, 0].std()), b=float(sims[:, 0].std() * ci),
                                            d=float(sims[:, 1].std())),
                                   vertex_dist=float(pc[1] + pc[0] * ci),
                                   rms_resid=float(np.sqrt(np.mean((inner[:, 1] - f_in(inner[:, 0], *pc))**2))))


# ------------------------------------------------------------------ jets
def strip_profile(im, pa, rmax, half_width, step=0.25):
    ue, un = pa_vec(pa)
    pe, pn = un, -ue                          # perpendicular
    s = np.arange(0, rmax, step)
    w = np.arange(-half_width, half_width + 1e-6, BIN)
    vals = np.array([sample(im, si * ue + w * pe, si * un + w * pn).mean() for si in s])
    return s, vals


def transverse(im, pa, s0, s1, half=8.0, step=0.125):
    """Average transverse profile between distances s0..s1 along PA."""
    ue, un = pa_vec(pa)
    pe, pn = un, -ue
    t = np.arange(-half, half + 1e-6, step)
    ss = np.arange(s0, s1, BIN)
    vals = np.array([sample(im, ss * ue + ti * pe, ss * un + ti * pn).mean() for ti in t])
    return t, vals


def fit_width(t, v):
    f = lambda t, A, t0, s, c0, c1: A * np.exp(-(t - t0)**2 / (2 * s**2)) + c0 + c1 * t
    p0 = [v.max() - np.median(v), t[np.argmax(v)], 1.0, np.median(v), 0]
    p, cov = curve_fit(f, t, v, p0=p0, maxfev=20000)
    fw = 2.3548 * abs(p[2]); fw_err = 2.3548 * np.sqrt(cov[2, 2])
    intrinsic = np.sqrt(max(fw**2 - PSF_FWHM**2, 0))
    return dict(fwhm_obs=float(fw), fwhm_err=float(fw_err), fwhm_intrinsic=float(intrinsic),
                offset=float(p[1]), params=[float(x) for x in p])


sm_j = smooth(img, 0.75 / BIN)
jets = {}
for name, pa in [("NW_jet", AXIS_PA), ("SE_counterjet", (AXIS_PA + 180) % 360)]:
    # refine the jet direction: PA (within +-25 deg of the axis) maximising
    # the brightness integrated over 2-10" from the pulsar
    trial = np.arange(pa - 25, pa + 25.1, 1)
    score = [strip_profile(sm_j, p, 10, 0.75)[1][8:].mean() for p in trial]
    pa_j = float(trial[int(np.argmax(score))])
    s, prof = strip_profile(sm_j, pa_j, 30, 1.0)
    # local background: same distances, +-(6-10)" off the jet axis
    s_b1, b1 = strip_profile(sm_j, pa_j + 25, 30, 1.0)
    s_b2, b2 = strip_profile(sm_j, pa_j - 25, 30, 1.0)
    bkg = 0.5 * (b1 + b2)
    contrast = prof / bkg
    # length: last distance (beyond the PSF wings, s > 1.5") where the jet is
    # >= 20% brighter than the neighbouring nebula, for 3 consecutive samples
    good = (contrast >= 1.2) & (s > 1.5)
    if name == "NW_jet":
        # beyond ~2" inside the outer-arc vertex the profile is dominated by
        # the arc, not the jet: stop the search there
        good &= s < OUT["outer_arc"]["vertex_dist"] - 2.0
    L = 0.0
    for k in range(len(s) - 3):
        if good[k:k + 3].all():
            L = s[k + 2]
    # width: transverse profile over the inner part of the jet
    seg = (1.5, max(2.5, min(L, 8.0)))
    t, v = transverse(img if False else smooth(img, 0.25 / BIN), pa_j, *seg)
    try:
        wd = fit_width(t, v)
    except Exception as ex:
        wd = dict(error=str(ex))
    jets[name] = dict(pa_deg=pa_j, length_arcsec=float(L), width_segment=seg, **wd,
                      profile_s=s.tolist(), profile=prof.tolist(), profile_bkg=bkg.tolist(),
                      transverse_t=t.tolist(), transverse_v=v.tolist())
OUT["jets"] = jets


# ------------------------------------------------------------------ bright PWN
# Wider 1-px image to reach the outer nebula
wide, whdr = load_image(f"{P}/img/band_broad_thresh.img")
wexp, _ = load_image(f"{P}/img/band_broad_thresh.expmap")
We, Wn = offsets_arcsec(whdr, wide.shape, OBSID)
wr = np.hypot(We, Wn)
cov_ok = wexp > 0.5 * np.nanmax(wexp)
rate = np.where(cov_ok, wide / np.maximum(wexp, 1), np.nan)     # counts/s/pixel
from scipy.ndimage import gaussian_filter
w = gaussian_filter(np.nan_to_num(rate), 4) / np.maximum(gaussian_filter(cov_ok.astype(float), 4), 1e-3)
w[~cov_ok] = np.nan
bg_ann = (wr > 90) & (wr < 120) & cov_ok
bg = float(np.nanmedian(w[bg_ann]))
OUT["outer_level_cts_s_pix"] = bg
levels = {}
for fac in (2, 4):
    m = (w > fac * bg) & cov_ok & (wr < 125)
    # moments of the mask -> equivalent ellipse
    ee, nn = We[m], Wn[m]
    c = np.cov(np.vstack([ee, nn]))
    lam, vec = np.linalg.eigh(c)
    a = 2 * np.sqrt(lam[1]); b = 2 * np.sqrt(lam[0])       # uniform-ellipse semi-axes
    pa = np.degrees(np.arctan2(vec[0, 1], vec[1, 1])) % 180
    area = m.sum() * PIX**2
    levels[f"x{fac}"] = dict(a=float(a), b=float(b), pa_major=float(pa),
                             centre_e=float(ee.mean()), centre_n=float(nn.mean()),
                             area_arcsec2=float(area), eq_radius=float(np.sqrt(area / np.pi)))
OUT["bright_pwn"] = levels

# Half-light radius of the PWN (pulsar r<2" removed, outer level subtracted)
rr = np.arange(2, 120, 1.0)
net = np.nan_to_num(rate) - bg
net[~cov_ok] = 0
cum = np.array([net[(wr >= 2) & (wr < r)].sum() for r in rr])
tot = cum[np.searchsorted(rr, 100)]
OUT["half_light_radius"] = float(rr[np.searchsorted(cum, 0.5 * tot)])
OUT["r90_radius"] = float(rr[np.searchsorted(cum, 0.9 * tot)])

# ------------------------------------------------------------------ physical
def phys(theta):
    return dict(arcsec=float(theta), pc=float(arcsec_to_pc(theta)),
                cm=float(arcsec_to_cm(theta)), au=float(arcsec_to_au(theta)),
                light_days=float(light_days(theta)))

tab = []
oa, ia = OUT["outer_arc"], OUT["inner_arc"]
tab += [("Outer arc / torus: semi-major axis a", oa["a"], oa["err"]["a"]),
        ("Outer arc / torus: semi-minor axis b", oa["b"], oa["err"]["b"]),
        ("Outer arc / torus: diameter 2a", 2 * oa["a"], 2 * oa["err"]["a"]),
        ("Inner arc: semi-major axis a", ia["a"], ia["err"]["a"]),
        ("Inner arc: semi-minor axis b", ia["b"], ia["err"]["b"]),
        ("NW jet: length", jets["NW_jet"]["length_arcsec"], 1.0),
        ("NW jet: width (FWHM, PSF-corrected)", jets["NW_jet"].get("fwhm_intrinsic", np.nan), jets["NW_jet"].get("fwhm_err", np.nan)),
        ("SE counter-jet: length", jets["SE_counterjet"]["length_arcsec"], 1.0),
        ("SE counter-jet: width (FWHM, PSF-corrected)", jets["SE_counterjet"].get("fwhm_intrinsic", np.nan), jets["SE_counterjet"].get("fwhm_err", np.nan)),
        ("Bright PWN (2x outer level): semi-major axis", levels["x2"]["a"], np.nan),
        ("Bright PWN (2x outer level): semi-minor axis", levels["x2"]["b"], np.nan),
        ("Bright PWN (4x outer level): semi-major axis", levels["x4"]["a"], np.nan),
        ("Bright PWN (4x outer level): semi-minor axis", levels["x4"]["b"], np.nan),
        ("PWN half-light radius", OUT["half_light_radius"], 1.0),
        ("PWN 90%-light radius", OUT["r90_radius"], 1.0)]
OUT["table"] = [dict(name=n, err=float(e) if np.isfinite(e) else None, **phys(v)) for n, v, e in tab]

with open(f"{P}/measurements.json", "w") as f:
    json.dump(OUT, f, indent=1)

print(f"PSF FWHM (pulsar, hard band): {PSF_FWHM:.2f}\"")
print(f"Symmetry axis PA: {AXIS_PA:.1f} deg")
for k in ("outer_arc", "inner_arc"):
    d = OUT[k]
    print(f"{k}: a={d['a']:.1f}+-{d['err']['a']:.1f}\" b={d['b']:.1f}+-{d['err']['b']:.1f}\" "
          f"d={d['d']:+.1f}+-{d['err']['d']:.1f}\" axis={d['axis_pa']:.0f}+-{d['err']['axis_pa']:.0f} "
          f"vertex={d['vertex_dist']:.1f}\" i={d['inclination_deg']:.0f}+-{d['inclination_err']:.0f} deg "
          f"rms={d['rms_resid']:.2f}\" n={d['npts']}")
for k, d in jets.items():
    print(f"{k}: PA={d['pa_deg']:.0f} L={d['length_arcsec']:.1f}\" "
          f"FWHM_obs={d.get('fwhm_obs', np.nan):.2f}+-{d.get('fwhm_err', np.nan):.2f}\" "
          f"intrinsic={d.get('fwhm_intrinsic', np.nan):.2f}\"")
for k, d in levels.items():
    print(f"bright PWN {k}: a={d['a']:.0f}\" b={d['b']:.0f}\" PA={d['pa_major']:.0f} "
          f"centre=({d['centre_e']:+.0f},{d['centre_n']:+.0f}) r_eq={d['eq_radius']:.0f}\"")
print(f"half-light r = {OUT['half_light_radius']:.0f}\", r90 = {OUT['r90_radius']:.0f}\"")
print()
for row in OUT["table"]:
    print(f"{row['name']:<48s} {row['arcsec']:6.1f}\"  {row['pc']:.4f} pc  {row['cm']:.2e} cm  {row['light_days']:6.1f} lt-d")
