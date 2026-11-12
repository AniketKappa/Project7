"""Final, annotated figures for the report.

  python scripts/09_figures.py
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

FIG = os.path.join(TOP, "figures")
OBS = ["128", "12075"]
LAB = {"128": "ObsID 128 (Apr 2000, 10.6 ks)", "12075": "ObsID 12075 (Sep 2010, 41.0 ks)"}
J = {o: json.load(open(f"{products(o)}/measurements.json")) for o in OBS}
ER = {o: json.load(open(f"{products(o)}/energy_results.json")) for o in OBS}
GI = {o: json.load(open(f"{products(o)}/photon_index.json")) for o in OBS}
OJ = {o: json.load(open(f"{products(o)}/outer_jet.json")) for o in OBS}
DOP = {o: json.load(open(f"{products(o)}/doppler.json")) for o in OBS}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5})


def asinh_norm(img, lo=1, hi=99.8, soft=0.05):
    v = img[np.isfinite(img)]
    a, b = np.percentile(v, lo), np.percentile(v, hi)
    x = np.clip((img - a) / (b - a), 0, None)
    return np.clip(np.arcsinh(x / soft) / np.arcsinh(1 / soft), 0, 1)


def ring_xy(d, full=True):
    t = np.linspace(0, 2 * np.pi, 400) if full else np.linspace(0, np.pi, 200)
    ta = np.radians(d["axis_pa"] if "axis_pa" in d else d["axis"]); ue, un = np.sin(ta), np.cos(ta)
    pe, pn = un, -ue
    e = d["d"] * ue + d["a"] * np.cos(t) * pe + d["b"] * np.sin(t) * ue
    n = d["d"] * un + d["a"] * np.cos(t) * pn + d["b"] * np.sin(t) * un
    return e, n


def compass(ax, x, y, L, color="w"):
    # east is to the LEFT: x axis decreases to the right
    ax.annotate("", xy=(x, y + L), xytext=(x, y), arrowprops=dict(arrowstyle="->", color=color, lw=1))
    ax.annotate("", xy=(x + L, y), xytext=(x, y), arrowprops=dict(arrowstyle="->", color=color, lw=1))
    ax.text(x, y + L * 1.15, "N", color=color, ha="center", va="bottom", fontsize=8)
    ax.text(x + L * 1.18, y, "E", color=color, ha="right", va="center", fontsize=8)


def scalebar(ax, x, y, L, color="w"):
    ax.plot([x, x - L], [y, y], color=color, lw=2)
    ax.text(x - L / 2, y + abs(L) * 0.12, f"{L:.0f}\" = {arcsec_to_pc(L):.3f} pc", color=color,
            ha="center", va="bottom", fontsize=7.5)


def lab(ax, text, xy, xytext, color="w"):
    ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=8, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.7, shrinkA=0, shrinkB=2),
                bbox=dict(fc="k", alpha=0.45, lw=0, pad=1.2))


# ------------------------------------------------------------------ Fig 1: overview (Task A)
fig, axs = plt.subplots(2, 2, figsize=(11, 11))
for row, o in enumerate(OBS):
    P = products(o)
    wide, wh = load_image(f"{P}/img/band_broad_flux.img")
    wexp, _ = load_image(f"{P}/img/band_broad_thresh.expmap")
    cov = wexp > 0.3 * wexp.max()
    w = smooth(np.where(cov, wide, 0), 1.2 / PIX) / np.maximum(smooth(cov.astype(float), 1.2 / PIX), 1e-3)
    w[~cov] = np.nan
    ext = extent_arcsec(wh, wide.shape, o)
    ax = axs[row, 0]
    ax.imshow(asinh_norm(w, 5, 99.95, 0.05), origin="lower", extent=ext, cmap="afmhot")
    ax.set_facecolor("k")
    ax.set_xlim(110, -125); ax.set_ylim(-110, 125)
    ax.add_patch(Rectangle((-32, -32), 64, 64, fill=False, ec="c", lw=0.8, ls="--"))
    lvl = J[o]["bright_pwn"]["x2"]
    ax.add_patch(Circle((0, 0), J[o]["r90_radius"], fill=False, ec="w", lw=0.6, ls=":"))
    p = np.array(OJ[o]["path"])
    if len(p):
        mid = p[len(p) // 2]
        lab(ax, "outer jet" if o == "12075" else "outer-jet blob", tuple(mid), (mid[0] + 10, mid[1] + 30))
    lab(ax, "diffuse PWN", (-20, -45), (30, -85))
    ax.text(105, -104, "dotted circle: radius enclosing 90% of the PWN light", color="w", fontsize=7, ha="left")
    compass(ax, 75, 88, 18); scalebar(ax, -70, -92, 30)
    ax.set_title(f"{LAB[o]} — 0.5–8 keV, exposure corrected")
    ax.set_xlabel("ΔRA (arcsec; east to the left)"); ax.set_ylabel("ΔDec (arcsec)")

    fine, fh = load_image(f"{P}/img/fine_broad.img")
    ext = extent_arcsec(fh, fine.shape, o)
    ax = axs[row, 1]
    ax.imshow(asinh_norm(smooth(fine, 0.5 / 0.246), 0.5, 99.99, 0.25), origin="lower", extent=ext, cmap="afmhot")
    ax.set_xlim(32, -32); ax.set_ylim(-32, 32)
    M = J[o]; AX = np.radians(M["axis_pa_deg"])
    oa = M["outer_arc"]; ia = M["inner_arc"]
    v = oa["vertex_dist"]
    lab(ax, "outer arc", (v * np.sin(AX) * 0.95 + 4, v * np.cos(AX) + 3), (-24, 27))
    t = np.radians(M["axis_pa_deg"] - 80); r = 8.6
    lab(ax, "inner arc", (r * np.sin(t), r * np.cos(t)), (-20, -14))
    jn = M["jets"]["NW_jet"]; t = np.radians(jn["pa_deg"])
    lab(ax, "NW jet", (6 * np.sin(t), 6 * np.cos(t)), (6, 20))
    jc = M["jets"]["SE_counterjet"]; t = np.radians(jc["pa_deg"])
    lab(ax, "SE counter-jet", (7 * np.sin(t), 7 * np.cos(t)), (12, -17))
    lab(ax, "pulsar", (0.3, -0.3), (16, 6))
    lab(ax, "diffuse PWN", (18, -6), (24, -27))
    compass(ax, 22, 20, 6); scalebar(ax, -16, -29, 10)
    ax.set_title(f"{LAB[o]} — inner nebula (0.25\" pixels)")
    ax.set_xlabel("ΔRA (arcsec)"); ax.set_ylabel("ΔDec (arcsec)")
plt.tight_layout()
plt.savefig(f"{FIG}/F1_overview.png", dpi=170)
plt.close()

# ------------------------------------------------------------------ Fig 2: measurements (Task B)
fig = plt.figure(figsize=(11.5, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])
for k, o in enumerate(OBS):
    P = products(o); M = J[o]
    fine, fh = load_image(f"{P}/img/fine_broad.img")
    ext = extent_arcsec(fh, fine.shape, o)
    ax = fig.add_subplot(gs[0, k])
    ax.imshow(asinh_norm(smooth(fine, 0.75 / 0.246), 2, 99.95, 0.03), origin="lower", extent=ext, cmap="gray")
    for name, c in [("outer_arc", "tab:orange"), ("inner_arc", "tab:cyan")]:
        d = M[name]; r = np.array(d["ridge"])
        ax.plot(r[:, 1] * np.sin(np.radians(r[:, 0])), r[:, 1] * np.cos(np.radians(r[:, 0])), ".", ms=2.5, color=c)
        if d["npts"] >= 25:          # free fit only where the arc is traced well
            e, n = ring_xy(d); ax.plot(e, n, "-", color=c, lw=0.8, alpha=0.9)
    ax.plot(*ring_xy(dict(M["inner_arc_fixed_incl"], axis_pa=M["axis_pa_deg"])), ":", color="tab:cyan", lw=0.8)
    for jn, c in [("NW_jet", "tab:red"), ("SE_counterjet", "tab:green")]:
        j = M["jets"][jn]; t = np.radians(j["pa_deg"]); L = j["length_arcsec"]
        w = j.get("fwhm_intrinsic", 2) / 2
        pe, pn = np.cos(t), -np.sin(t)
        for sgn in (-1, 1):
            ax.plot([1.5 * np.sin(t) + sgn * w * pe, L * np.sin(t) + sgn * w * pe],
                    [1.5 * np.cos(t) + sgn * w * pn, L * np.cos(t) + sgn * w * pn], "-", color=c, lw=0.9)
    ta = np.radians(M["axis_pa_deg"])
    ax.plot([-40 * np.sin(ta), 40 * np.sin(ta)], [-40 * np.cos(ta), 40 * np.cos(ta)], "w--", lw=0.4)
    ax.set_xlim(30, -30); ax.set_ylim(-30, 30)
    ax.set_title(f"{LAB[o]}\nring fits (outer: orange, inner: cyan; dotted = inclination fixed to outer ring)")
    ax.set_xlabel("ΔRA (arcsec)"); ax.set_ylabel("ΔDec (arcsec)")
ax = fig.add_subplot(gs[1, 0])
for o, ls in zip(OBS, ["-", "--"]):
    for jn, c, sgn in [("NW_jet", "tab:red", 1), ("SE_counterjet", "tab:green", -1)]:
        j = J[o]["jets"][jn]; s = np.array(j["profile_s"]); p = np.array(j["profile"]); b = np.array(j["profile_bkg"])
        nrm = p[0]
        ax.plot(sgn * s, p / nrm, ls, color=c, lw=1, label=f"{jn.replace('_', ' ')} – {o}")
        ax.plot(sgn * s, b / nrm, ls, color="0.5", lw=0.7)
    v = J[o]["outer_arc"]["vertex_dist"]
    ax.axvline(v, color="tab:orange", lw=0.6, ls=ls)
ax.set_yscale("log"); ax.set_xlabel("distance along the symmetry axis (arcsec; + = NW, − = SE)")
ax.set_ylabel("surface brightness (normalised)")
ax.set_title("Brightness along the jet axis (grey: same distance, ±25° off axis)")
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3); ax.set_xlim(-25, 25)
ax.text(J["128"]["outer_arc"]["vertex_dist"] + 0.4, 0.55, "outer-arc\nvertex", color="tab:orange", fontsize=7)
ax = fig.add_subplot(gs[1, 1])
for o, ls in zip(OBS, ["-", "--"]):
    for jn, c in [("NW_jet", "tab:red"), ("SE_counterjet", "tab:green")]:
        j = J[o]["jets"][jn]; t = np.array(j["transverse_t"]); v = np.array(j["transverse_v"])
        A, t0, s, c0, c1 = j["params"]
        base = c0 + c1 * t
        ax.plot(t, (v - base) / A, ls, color=c, lw=1,
                label=f"{jn.replace('_', ' ')} – {o}: FWHM {j['fwhm_obs']:.1f}\" (intr. {j['fwhm_intrinsic']:.1f}\")")
psf = J["128"]["psf_fwhm_arcsec"]; sg = psf / 2.3548
tt = np.linspace(-8, 8, 300); ax.plot(tt, np.exp(-tt**2 / 2 / sg**2), ":", color="k", lw=1, label=f"PSF (pulsar) {psf:.2f}\" FWHM")
ax.set_xlabel("distance across the jet (arcsec)"); ax.set_ylabel("background-subtracted profile (norm.)")
ax.set_title("Transverse jet profiles with Gaussian + linear fits")
ax.legend(fontsize=6.5); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIG}/F2_measurements.png", dpi=170)
plt.close()

# ------------------------------------------------------------------ Fig 3: photon index vs distance (Task D)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
cols = {"NW (along jet axis)": "tab:red", "SE (counter-jet side)": "tab:green",
        "NE (along torus major axis)": "tab:blue", "SW (along torus major axis)": "tab:orange"}
for ax, o in zip(axs[:2], OBS):
    for s, c in cols.items():
        pr = GI[o]["radial"][s]
        r = np.array([x["r"] for x in pr]); g = np.array([x["gamma"] for x in pr])
        e1 = np.array([x["em"] for x in pr]); e2 = np.array([x["ep"] for x in pr])
        ax.errorbar(r * (1 + 0.015 * list(cols).index(s)), g, yerr=[e1, e2], fmt="o-", ms=3, lw=0.9, color=c, label=s, capsize=1.5)
    pr = GI[o]["radial"]["all directions"]
    ax.plot([x["r"] for x in pr], [x["gamma"] for x in pr], "k-", lw=2, alpha=0.6, label="all directions")
    ax.axvspan(60, 120, color="0.85", zorder=0)
    ax.text(62, 1.18, "PWN faint:\nbackground\nmatters", fontsize=7, color="0.3")
    ax.set_xscale("log"); ax.set_xlim(3, 115); ax.set_ylim(1.1, 1.9)
    ax.set_xlabel("distance from pulsar (arcsec)"); ax.set_ylabel("photon index Γ (from HR, nH = 3×10²⁰)")
    ax.set_title(LAB[o]); ax.grid(alpha=0.3, which="both")
    sec = ax.secondary_xaxis("top", functions=(lambda x: arcsec_to_pc(x), lambda p: p * 206265 / DIST_PC))
    sec.set_xlabel("projected distance (pc)")
axs[0].legend(fontsize=7, loc="upper left")
ax = axs[2]
x = np.arange(3)
for k, o in enumerate(OBS):
    sz = ER[o]["band_sizes"]
    ax.bar(x + (k - 0.5) * 0.36, [sz[b]["r50"] for b in BANDS], 0.34, color=["tab:red", "tab:green", "tab:blue"],
           alpha=0.55 if k == 0 else 0.9, edgecolor="k", lw=0.5, label=LAB[o].split(" (")[0])
    for i, b in enumerate(BANDS):
        ax.text(i + (k - 0.5) * 0.36, sz[b]["r50"] + 0.2, f"{sz[b]['r50']:.1f}", ha="center", fontsize=7)
ax.set_xticks(x, [f"{b}\n{BAND_LABEL[b]}" for b in BANDS]); ax.set_ylim(15, 25)
ax.set_ylabel("half-light radius of the PWN (arcsec)")
ax.set_title("Nebula size per band (pale: ObsID 128, dark: 12075)")
ax.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(f"{FIG}/F5_gamma_profiles.png", dpi=170)
plt.close()

# ------------------------------------------------------------------ Fig: light curves (cleaning)
fig, axs = plt.subplots(1, 2, figsize=(12, 3.3))
for ax, o in zip(axs, OBS):
    with fits.open(f"{products(o)}/bkg_lc.fits") as h:
        d = h["LIGHTCURVE"].data if "LIGHTCURVE" in h else h[1].data
        t = d["TIME"] - d["TIME"][0]; rate = d["COUNT_RATE"]; err = d["COUNT_RATE_ERR"]
    with fits.open(f"{products(o)}/bkg_lc.gti") as h:
        g = h["GTI"].data
    t0 = d["TIME"][0]
    ax.errorbar(t / 1e3, rate, yerr=err, fmt=".", ms=3, color="k", lw=0.6)
    good = np.zeros_like(rate, bool)
    for a, b in zip(g["START"], g["STOP"]):
        good |= (d["TIME"] >= a) & (d["TIME"] < b)
    ax.plot(t[~good] / 1e3, rate[~good], "rx", ms=5, label="removed")
    m = np.mean(rate[good & (rate > 0)])
    ax.axhline(m, color="tab:blue", lw=0.8); ax.axhspan(m / 1.2, m * 1.2, color="tab:blue", alpha=0.1)
    ax.set_xlabel("time (ks)"); ax.set_ylabel("rate (counts/s), 0.5–7 keV")
    ax.set_title(f"{LAB[o]}: annulus 2.0–2.7′ around the pulsar")
    ax.set_ylim(0, max(rate) * 1.2); ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig(f"{FIG}/F0_lightcurves.png", dpi=170)
plt.close()

# ------------------------------------------------------------------ Fig: ring brightness vs azimuth (Doppler)
fig, ax = plt.subplots(figsize=(6.2, 3.4))
for o, c in zip(OBS, ["tab:gray", "k"]):
    d = DOP[o]; t = np.array(d["profile_t"]); v = np.array(d["profile_excess"])
    order = np.argsort(((t + 90) % 360))
    tt = ((t + 90) % 360) - 90
    srt = np.argsort(tt)
    ax.plot(tt[srt], v[srt] / d["front_excess"], "-", color=c, lw=1, label=LAB[o])
ax.axvspan(30, 65, color="tab:red", alpha=0.08); ax.axvspan(115, 150, color="tab:red", alpha=0.08)
ax.axvspan(-150, -115, color="tab:blue", alpha=0.08); ax.axvspan(-65, -30, color="tab:blue", alpha=0.08)
ax.axhline(0, color="k", lw=0.5)
ax.set_xlim(-90, 270)
ax.set_xticks([-90, 0, 90, 180, 270], ["−90\n(back,\nSE)", "0\n(NE end)", "90\n(front,\nNW vertex)", "180\n(SW end)", "270\n(back)"])
ax.set_xlabel("azimuth along the fitted outer ring (deg)")
ax.set_ylabel("ring brightness above\nlocal nebula (÷ front mean)")
ax.set_title("Front (red bands) vs back (blue bands) of the outer ring")
ax.legend(fontsize=7); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIG}/F7_ring_doppler.png", dpi=170)
plt.close()
print("figures written to", FIG)

# ------------------------------------------------------------------ Fig 4: colour-balanced RGB (Task D)
# Each band is divided by its mean surface brightness over the bright PWN
# (3" < r < 45"), then one common asinh stretch is applied to all three
# channels.  A pixel with the average PWN spectrum is therefore grey/white;
# red = relatively soft, blue = relatively hard.
def balanced_rgb(imgs, E_, N_, sigma, stretch, cov=None):
    R_ = np.hypot(E_, N_)
    sel = (R_ > 3) & (R_ < 45)
    if cov is not None:
        sel &= cov
    ch = []
    for b in BANDS:
        sm = smooth(np.nan_to_num(imgs[b]), sigma)
        ch.append(sm / sm[sel].mean())
    ch = np.dstack(ch)
    lum = ch.mean(axis=2, keepdims=True)
    out = np.arcsinh(lum / stretch) / np.arcsinh(lum.max() / stretch * 0.6) * ch / np.maximum(lum, 1e-9)
    return np.clip(out, 0, 1)


fig, axs = plt.subplots(2, 2, figsize=(11, 11))
for row, o in enumerate(["12075", "128"]):
    P = products(o)
    fine = {b: load_image(f"{P}/img/fine_{b}.img")[0] for b in BANDS}
    _, fh = load_image(f"{P}/img/fine_broad.img")
    Ef, Nf = offsets_arcsec(fh, fine["soft"].shape, o)
    ext = extent_arcsec(fh, fine["soft"].shape, o)
    ax = axs[row, 0]
    ax.imshow(balanced_rgb(fine, Ef, Nf, 0.9 / 0.246, 0.6), origin="lower", extent=ext)
    ax.set_xlim(32, -32); ax.set_ylim(-32, 32)
    ax.set_title(f"{LAB[o]} — inner nebula, σ = 0.9\"")
    flux = {b: load_image(f"{P}/img/band_{b}_flux.img")[0] for b in BANDS}
    ex, wh = load_image(f"{P}/img/band_broad_thresh.expmap")
    cov = ex > 0.5 * ex.max()
    for b in BANDS:
        flux[b] = np.where(cov, flux[b], 0)
    Ew, Nw = offsets_arcsec(wh, ex.shape, o)
    ax = axs[row, 1]
    img = balanced_rgb(flux, Ew, Nw, 2.0 / PIX, 0.25, cov)
    img[~cov] = 0
    ax.imshow(img, origin="lower", extent=extent_arcsec(wh, ex.shape, o))
    ax.set_xlim(100, -120); ax.set_ylim(-100, 120)
    ax.set_title(f"{LAB[o]} — whole nebula, σ = 2\"")
    for a in axs[row]:
        a.set_xlabel("ΔRA (arcsec; east to the left)"); a.set_ylabel("ΔDec (arcsec)")
        a.text(0.02, 0.97, "R 0.5–1.5   G 1.5–3   B 3–8 keV\n(white = mean PWN colour)", transform=a.transAxes,
               color="w", va="top", fontsize=8, bbox=dict(fc="k", alpha=0.5, lw=0))
plt.tight_layout()
plt.savefig(f"{FIG}/F4_rgb.png", dpi=170)
plt.close()
print("F4 written")

# ------------------------------------------------------------------ split Fig 1 into wide (1a) and inner (1b) rows for the report
from PIL import Image as _PIL
_im = _PIL.open(f"{FIG}/F1_overview.png"); _W, _H = _im.size
_tl = _im.crop((0, 0, _W // 2, _H // 2)); _bl = _im.crop((0, _H // 2, _W // 2, _H))
_tr = _im.crop((_W // 2, 0, _W, _H // 2)); _br = _im.crop((_W // 2, _H // 2, _W, _H))
for name, (a, b) in {"F1a_wide.png": (_tl, _bl), "F1b_inner.png": (_tr, _br)}.items():
    out = _PIL.new("RGB", (a.width + b.width, a.height), "white"); out.paste(a, (0, 0)); out.paste(b, (a.width, 0))
    out.save(f"{FIG}/{name}")
print("F1a/F1b written")
