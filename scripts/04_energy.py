"""Tasks C and D: energy-resolved imaging, hardness ratios and colour images.

  python scripts/04_energy.py 128

Bands: soft S = 0.5-1.5 keV, medium M = 1.5-3 keV, hard H = 3-8 keV.
Colours (counts based, the standard X-ray "hardness ratio"):
    HR1 = (M - S) / (M + S)      HR2 = (H - M) / (H + M)
    HR  = (H - S) / (H + S)
Larger HR -> harder (flatter) spectrum.  Errors from Poisson statistics.
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

OBSID = sys.argv[1] if len(sys.argv) > 1 else "128"
P = products(OBSID)
FIG = os.path.join(TOP, "figures")
os.makedirs(FIG, exist_ok=True)
M = json.load(open(f"{P}/measurements.json"))
AX = M["axis_pa_deg"]
RES = {}

# ------------------------------------------------------------------ data
fine = {b: load_image(f"{P}/img/fine_{b}.img")[0] for b in BANDS + ["broad"]}
_, fh = load_image(f"{P}/img/fine_broad.img")
E, N = offsets_arcsec(fh, fine["broad"].shape, OBSID)
R = np.hypot(E, N)
PAimg = np.degrees(np.arctan2(E, N)) % 360
FBIN = abs(N[1, 0] - N[0, 0])
ext_f = extent_arcsec(fh, fine["broad"].shape, OBSID)

wide = {b: load_image(f"{P}/img/band_{b}_thresh.img")[0] for b in BANDS + ["broad"]}
wexp = {b: load_image(f"{P}/img/band_{b}_thresh.expmap")[0] for b in BANDS + ["broad"]}
wflux = {b: load_image(f"{P}/img/band_{b}_flux.img")[0] for b in BANDS + ["broad"]}
_, wh = load_image(f"{P}/img/band_broad_thresh.img")
WE, WN = offsets_arcsec(wh, wide["broad"].shape, OBSID)
WR = np.hypot(WE, WN)
WPA = np.degrees(np.arctan2(WE, WN)) % 360
ext_w = extent_arcsec(wh, wide["broad"].shape, OBSID)
cov = wexp["broad"] > 0.5 * wexp["broad"].max()


def dpa(pa, ref):
    return (pa - ref + 180) % 360 - 180


def hr(a, b):
    """(b - a)/(b + a) with Poisson error."""
    a = float(a); b = float(b)
    t = a + b
    if t <= 0:
        return np.nan, np.nan
    v = (b - a) / t
    e = 2.0 / t**2 * np.sqrt(b**2 * a + a**2 * b)
    return v, e


def counts(mask, imgs=fine):
    return {b: float(imgs[b][mask].sum()) for b in BANDS}


# ------------------------------------------------------------------ regions
ring = M["outer_arc"]


def ring_mask(a, b, d, axis, width, Ee=E, Nn=N, front_only=True):
    """Pixels within `width` arcsec of the projected ring ellipse."""
    ta = np.radians(axis)
    ue, un = np.sin(ta), np.cos(ta)
    along = Ee * ue + Nn * un - d
    perp = Ee * un - Nn * ue
    # approximate distance to the ellipse via normalised radius
    rho = np.sqrt((perp / a)**2 + (along / b)**2)
    dist = np.abs(rho - 1) * np.sqrt((a * b))      # rough (fine for width ~ 2")
    m = dist < width
    if front_only:
        m &= along > -0.3 * b
    return m


psr = R < 1.5
knot = (R >= 1.5) & (R < 3.5) & (np.abs(dpa(PAimg, AX)) < 35)
jet_nw = (R >= 3.5) & (R < min(M["jets"]["NW_jet"]["length_arcsec"], ring["vertex_dist"] - 2.5)) & \
         (np.abs((R * np.sin(np.radians(dpa(PAimg, M["jets"]["NW_jet"]["pa_deg"]))))) < 1.5) & \
         (np.abs(dpa(PAimg, M["jets"]["NW_jet"]["pa_deg"])) < 90)
cj_pa = M["jets"]["SE_counterjet"]["pa_deg"]
cjet = (R >= 2.0) & (R < M["jets"]["SE_counterjet"]["length_arcsec"]) & \
       (np.abs(R * np.sin(np.radians(dpa(PAimg, cj_pa)))) < 1.5) & (np.abs(dpa(PAimg, cj_pa)) < 90)
outer_arc = ring_mask(ring["a"], ring["b"], ring["d"], ring["axis_pa"], 1.5) & (R > 3) & \
            (np.abs(dpa(PAimg, AX)) < 75)
inner_ends = ((np.abs(R - 8.3) < 1.5) & ((np.abs(dpa(PAimg, AX)) > 60) & (np.abs(dpa(PAimg, AX)) < 90)))
used = psr | knot | jet_nw | cjet | outer_arc | inner_ends
# diffuse nebula in annuli, excluding the structures above
diffuse = {}
for r0, r1 in [(3, 15), (15, 30), (30, 45), (45, 60)]:
    diffuse[f"diffuse {r0}-{r1}\""] = (R >= r0) & (R < r1) & ~used & ~ring_mask(ring["a"], ring["b"], ring["d"], ring["axis_pa"], 3.0)

regions = {"pulsar (r<1.5\")": psr, "inner knot (1.5-3.5\" NW)": knot, "NW jet": jet_nw,
           "SE counter-jet": cjet, "inner arc (ends)": inner_ends, "outer arc": outer_arc}
regions.update(diffuse)


# local background for each compact structure: the same radii, rotated +-40 deg
def local_bkg(mask, rot=40):
    out = np.zeros_like(mask)
    for sgn in (+1, -1):
        pa_r = (PAimg - sgn * rot) % 360
        # pixel at (R, PA) is background if (R, PA - rot) lies in the region
        th = np.radians(pa_r)
        e2, n2 = R * np.sin(th), R * np.cos(th)
        j = np.round((e2 - E[0, 0]) / (E[0, 1] - E[0, 0])).astype(int)
        i = np.round((n2 - N[0, 0]) / (N[1, 0] - N[0, 0])).astype(int)
        ok = (i >= 0) & (i < mask.shape[0]) & (j >= 0) & (j < mask.shape[1])
        tmp = np.zeros_like(mask)
        tmp[ok] = mask[i[ok], j[ok]]
        out |= tmp
    return out & ~used


table = []
for name, m in regions.items():
    c = counts(m)
    row = dict(name=name, area=float(m.sum() * FBIN**2), **{k: c[k] for k in BANDS})
    row["HR1"], row["HR1_err"] = hr(c["soft"], c["medium"])
    row["HR2"], row["HR2_err"] = hr(c["medium"], c["hard"])
    row["HR"], row["HR_err"] = hr(c["soft"], c["hard"])
    # median photon energy is not available from images; use the mean band
    if not name.startswith("diffuse") and not name.startswith("pulsar"):
        bm = local_bkg(m)
        cb = counts(bm)
        scale = m.sum() / max(bm.sum(), 1)
        net = {k: c[k] - scale * cb[k] for k in BANDS}
        row.update({f"net_{k}": net[k] for k in BANDS})
        row["net_HR"], row["net_HR_err"] = hr(max(net["soft"], 1e-3), max(net["hard"], 1e-3))
        # errors including the background
        tS = c["soft"] + scale**2 * cb["soft"]; tH = c["hard"] + scale**2 * cb["hard"]
        s_, h_ = net["soft"], net["hard"]
        if s_ + h_ > 0:
            row["net_HR_err"] = float(2 / (s_ + h_)**2 * np.sqrt(h_**2 * tS + s_**2 * tH))
        row["bkg_HR"], row["bkg_HR_err"] = hr(cb["soft"], cb["hard"])
    table.append(row)
RES["regions"] = table

# ------------------------------------------------------------------ radial hardness profiles in sectors
sectors = {"NW (along jet axis)": AX, "SE (counter-jet side)": (AX + 180) % 360,
           "NE (along torus major axis)": (AX + 90) % 360, "SW (along torus major axis)": (AX - 90) % 360}
edges = np.array([2, 6, 10, 15, 20, 25, 30, 40, 50, 60, 75, 90, 110])
prof = {}
for sname, spa in sectors.items():
    rows = []
    for r0, r1 in zip(edges[:-1], edges[1:]):
        m = (WR >= r0) & (WR < r1) & (np.abs(dpa(WPA, spa)) < 45) & cov
        c = {b: float(wide[b][m].sum()) for b in BANDS}
        v, e = hr(c["soft"], c["hard"])
        v1, e1 = hr(c["soft"], c["medium"])
        rows.append(dict(r=float(0.5 * (r0 + r1)), r0=float(r0), r1=float(r1), HR=v, HR_err=e,
                         HR1=v1, HR1_err=e1, **c,
                         sb={b: float(wide[b][m].sum() / max(wexp[b][m].mean(), 1) / (m.sum() * PIX**2))
                             for b in BANDS}))
    prof[sname] = rows
# all directions
rows = []
for r0, r1 in zip(edges[:-1], edges[1:]):
    m = (WR >= r0) & (WR < r1) & cov
    c = {b: float(wide[b][m].sum()) for b in BANDS}
    v, e = hr(c["soft"], c["hard"]); v1, e1 = hr(c["soft"], c["medium"])
    rows.append(dict(r=float(0.5 * (r0 + r1)), r0=float(r0), r1=float(r1), HR=v, HR_err=e, HR1=v1, HR1_err=e1, **c,
                     sb={b: float(wide[b][m].sum() / max(wexp[b][m].mean(), 1) / (m.sum() * PIX**2)) for b in BANDS}))
prof["all directions"] = rows
RES["radial_profiles"] = prof

# ------------------------------------------------------------------ size of the nebula per band
outer_lvl = (WR > 90) & (WR < 120) & cov
sizes = {}
rr = np.arange(2, 101, 1.0)
for b in BANDS + ["broad"]:
    rate = np.where(cov, wide[b] / np.maximum(wexp[b], 1), 0)
    bg = np.mean(rate[outer_lvl])
    net = np.where(cov, rate - bg, 0)
    cum = np.array([net[(WR >= 2) & (WR < r)].sum() for r in rr])
    tot = cum[-1]
    r50 = float(np.interp(0.5 * tot, cum, rr)); r80 = float(np.interp(0.8 * tot, cum, rr))
    # extent where the smoothed band image exceeds 2x its outer level
    sm = smooth(rate, 4) / np.maximum(smooth(cov.astype(float), 4), 1e-3)
    m2 = (sm > 2 * bg) & cov & (WR < 125)
    sizes[b] = dict(r50=r50, r80=r80, eq_radius_2x=float(np.sqrt(m2.sum() * PIX**2 / np.pi)),
                    pwn_fraction=float(tot / max(np.where(cov, rate, 0)[(WR >= 2) & (WR < 100)].sum(), 1e-9)))
RES["band_sizes"] = sizes

json.dump(RES, open(f"{P}/energy_results.json", "w"), indent=1, default=float)

# ================================================================== figures
def style(ax, lim, label=True):
    ax.set_xlim(lim, -lim); ax.set_ylim(-lim, lim)
    if label:
        ax.set_xlabel("ΔRA (arcsec; east ←)"); ax.set_ylabel("ΔDec (arcsec)")


def asinh_norm(img, lo_pct=1, hi_pct=99.8, soft=0.05):
    v = img[np.isfinite(img)]
    lo, hi = np.percentile(v, lo_pct), np.percentile(v, hi_pct)
    x = np.clip((img - lo) / (hi - lo), 0, None)
    return np.clip(np.arcsinh(x / soft) / np.arcsinh(1 / soft), 0, 1)


def ring_patch(ax, d, **kw):
    t = np.linspace(0, 2 * np.pi, 400)
    ta = np.radians(d["axis_pa"]); ue, un = np.sin(ta), np.cos(ta)
    pe, pn = un, -ue
    e = d["d"] * ue + d["a"] * np.cos(t) * pe + d["b"] * np.sin(t) * ue
    n = d["d"] * un + d["a"] * np.cos(t) * pn + d["b"] * np.sin(t) * un
    ax.plot(e, n, **kw)


# --- Fig: three bands, exposure-corrected, matched smoothing, two zooms
fig, axs = plt.subplots(2, 3, figsize=(15, 10.2))
for k, b in enumerate(BANDS):
    fl = fine[b] / np.maximum(np.interp(0, [0], [1]), 1)   # counts; exposure ~flat inside 1'
    sm = smooth(fine[b], 1.0 / FBIN)
    ax = axs[0, k]
    ax.imshow(asinh_norm(sm, 1, 99.9, 0.03), origin="lower", extent=ext_f, cmap="inferno")
    style(ax, 35, label=(k == 0)); ax.set_title(f"{b.capitalize()}  {BAND_LABEL[b]}  ({int(fine[b].sum()):,} counts in view)")
    ax = axs[1, k]
    fl = np.where(cov, wflux[b], np.nan)
    sm = smooth(np.nan_to_num(fl), 3) / np.maximum(smooth(cov.astype(float), 3), 1e-3)
    sm[~cov] = np.nan
    ax.imshow(asinh_norm(sm, 2, 99.7, 0.05), origin="lower", extent=ext_w, cmap="inferno")
    style(ax, 120, label=(k == 0))
    ax.set_title(f"{b.capitalize()} – exposure-corrected, σ = 1.5\"")
for ax in axs[0]:
    ring_patch(ax, ring, color="c", lw=0.6, ls="--", alpha=0.6)
plt.tight_layout()
plt.savefig(f"{FIG}/fig_bands_{OBSID}.png", dpi=150)
plt.close()

# --- Fig: RGB three-colour images (R = soft, G = medium, B = hard)
def rgb_from(imgs, sigma, lo, hi, soft):
    chans = []
    for b in BANDS:
        sm = smooth(np.nan_to_num(imgs[b]), sigma)
        chans.append(asinh_norm(sm, lo, hi, soft))
    return np.dstack(chans)


fig, axs = plt.subplots(1, 2, figsize=(15, 7.4))
rgb_in = rgb_from(fine, 0.9 / FBIN, 1, 99.85, 0.04)
axs[0].imshow(rgb_in, origin="lower", extent=ext_f)
style(axs[0], 32); axs[0].set_title("Inner nebula (0.25\" pixels, σ = 0.9\")")
wfl = {b: np.where(cov, wflux[b], 0) for b in BANDS}
rgb_out = rgb_from(wfl, 2.5, 2, 99.8, 0.06)
rgb_out[~cov] = 0
axs[1].imshow(rgb_out, origin="lower", extent=ext_w)
style(axs[1], 120); axs[1].set_title("Whole field (exposure-corrected, σ = 1.2\")")
for ax in axs:
    ax.text(0.02, 0.97, "R: 0.5–1.5 keV   G: 1.5–3 keV   B: 3–8 keV", transform=ax.transAxes,
            color="w", va="top", fontsize=10, bbox=dict(fc="k", alpha=0.5, lw=0))
plt.tight_layout()
plt.savefig(f"{FIG}/fig_rgb_{OBSID}.png", dpi=160)
plt.close()

# --- Fig: hardness map HR = (H-S)/(H+S), 1" pixels, Gaussian sigma 2"
def rebin(a, f):
    n0, n1 = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    return a[:n0, :n1].reshape(n0 // f, f, n1 // f, f).sum(axis=(1, 3))


S = smooth(rebin(fine["soft"], 4), 2.0); H = smooth(rebin(fine["hard"], 4), 2.0)
Mm = smooth(rebin(fine["medium"], 4), 2.0)
ksum = 2 * np.pi * 2.0**2                 # counts in the Gaussian kernel ~ value x 25
hmap = (H - S) / np.maximum(H + S, 1e-9)
err = 2 / (H + S) * np.sqrt(H * S / (H + S)) / np.sqrt(ksum) * np.sqrt(ksum)  # approx
hmap[(H + S) * ksum < 150] = np.nan
hm1 = (Mm - S) / np.maximum(Mm + S, 1e-9); hm1[(Mm + S) * ksum < 150] = np.nan
fig, axs = plt.subplots(1, 2, figsize=(14, 6.4))
for ax, mp, t in zip(axs, [hmap, hm1], ["HR = (H − S)/(H + S)", "HR1 = (M − S)/(M + S)"]):
    im = ax.imshow(mp, origin="lower", extent=ext_f, cmap="coolwarm_r" if False else "RdYlBu",
                   vmin=np.nanpercentile(mp, 3), vmax=np.nanpercentile(mp, 97))
    style(ax, 45)
    ring_patch(ax, ring, color="k", lw=0.6, ls="--")
    cb = plt.colorbar(im, ax=ax, shrink=0.85); cb.set_label(t + "   (blue = harder)")
    ax.contour(np.linspace(ext_f[0], ext_f[1], fine["broad"].shape[1]),
               np.linspace(ext_f[2], ext_f[3], fine["broad"].shape[0]),
               smooth(fine["broad"], 1.2 / FBIN), levels=np.percentile(smooth(fine["broad"], 1.2 / FBIN), [80, 93, 98]),
               colors="k", linewidths=0.4, alpha=0.6)
    ax.set_title(t + "  (1\" pixels, σ = 2\")")
plt.tight_layout()
plt.savefig(f"{FIG}/fig_hardness_map_{OBSID}.png", dpi=150)
plt.close()

# --- Fig: hardness vs distance, and normalised surface-brightness profiles per band
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4))
cols = {"NW (along jet axis)": "tab:red", "SE (counter-jet side)": "tab:blue",
        "NE (along torus major axis)": "tab:green", "SW (along torus major axis)": "tab:orange",
        "all directions": "k"}
for s, rows in prof.items():
    r = [x["r"] for x in rows]; v = [x["HR"] for x in rows]; e = [x["HR_err"] for x in rows]
    axs[0].errorbar(np.array(r) * (1 + 0.012 * list(cols).index(s)), v, yerr=e, fmt="o-", ms=3, lw=1,
                    color=cols[s], label=s, capsize=2)
axs[0].set_xscale("log"); axs[0].set_xlabel("distance from pulsar (arcsec)")
axs[0].set_ylabel("HR = (H − S)/(H + S)"); axs[0].legend(fontsize=8); axs[0].grid(alpha=.3)
secax = axs[0].secondary_xaxis("top", functions=(lambda x: arcsec_to_pc(x), lambda p: p * 206265 / DIST_PC))
secax.set_xlabel("projected distance (pc, d = 287 pc)")
rows = prof["all directions"]
r = np.array([x["r"] for x in rows])
for b, c in zip(BANDS, ["tab:red", "tab:green", "tab:blue"]):
    sb = np.array([x["sb"][b] for x in rows])
    axs[1].loglog(r, sb / sb[np.argmin(np.abs(r - 12.5))], "o-", color=c, ms=3, label=f"{b} {BAND_LABEL[b]}")
axs[1].set_xlabel("distance from pulsar (arcsec)"); axs[1].set_ylabel("surface brightness (normalised at 12.5\")")
axs[1].legend(); axs[1].grid(alpha=.3, which="both")
plt.tight_layout()
plt.savefig(f"{FIG}/fig_hardness_profiles_{OBSID}.png", dpi=150)
plt.close()

# ------------------------------------------------------------------ print summary
print(f"{'region':<28s} {'S':>7s} {'M':>7s} {'H':>7s}   {'HR1':>13s} {'HR2':>13s} {'HR':>13s}  {'net HR':>13s} {'bkg HR':>13s}")
for r in table:
    line = (f"{r['name']:<28s} {r['soft']:7.0f} {r['medium']:7.0f} {r['hard']:7.0f}   "
            f"{r['HR1']:+.3f}±{r['HR1_err']:.3f} {r['HR2']:+.3f}±{r['HR2_err']:.3f} {r['HR']:+.3f}±{r['HR_err']:.3f}")
    if "net_HR" in r:
        line += f"  {r['net_HR']:+.3f}±{r['net_HR_err']:.3f} {r['bkg_HR']:+.3f}±{r['bkg_HR_err']:.3f}"
    print(line)
print()
for s, rows in prof.items():
    print(s, " ".join(f"{x['r']:.0f}:{x['HR']:+.2f}±{x['HR_err']:.2f}" for x in rows))
print()
for b, d in sizes.items():
    print(f"{b:<7s} r50={d['r50']:.1f}\" r80={d['r80']:.1f}\" r_eq(2x)={d['eq_radius_2x']:.1f}\" PWN frac={d['pwn_fraction']:.2f}")
