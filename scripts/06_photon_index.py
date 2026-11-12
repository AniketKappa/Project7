"""Convert band count ratios into power-law photon indices.

  python scripts/06_photon_index.py 128

For an absorbed power law N(E) = K E^-Gamma exp(-nH sigma(E)) folded through
the observation's ARF and RMF, compute the expected counts in the soft
(0.5-1.5), medium (1.5-3) and hard (3-8 keV) bands as a function of Gamma,
and invert the observed hardness ratios.  nH = 3e20 cm^-2 (Pavlov et al. 2001).
Cross-sections: Morrison & McCammon (1983).
"""
import sys, os, json
import numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

OBSID = sys.argv[1] if len(sys.argv) > 1 else "128"
P = products(OBSID)
NH = 3.0e20

# Morrison & McCammon (1983) table: E_lo, E_hi (keV), c0, c1, c2
MM83 = np.array([
    [0.030, 0.100, 17.3, 608.1, -2150.0], [0.100, 0.284, 34.6, 267.9, -476.1],
    [0.284, 0.400, 78.1, 18.8, 4.3], [0.400, 0.532, 71.4, 66.8, -51.4],
    [0.532, 0.707, 95.5, 145.8, -61.1], [0.707, 0.867, 308.9, -380.6, 294.0],
    [0.867, 1.303, 120.6, 169.3, -47.7], [1.303, 1.840, 141.3, 146.8, -31.5],
    [1.840, 2.471, 202.7, 104.7, -17.0], [2.471, 3.210, 342.7, 18.7, 0.0],
    [3.210, 4.038, 352.2, 18.7, 0.0], [4.038, 7.111, 433.9, -2.4, 0.75],
    [7.111, 8.331, 629.0, 30.9, 0.0], [8.331, 10.00, 701.2, 25.2, 0.0],
    [10.00, 20.00, 701.2, 25.2, 0.0]])


def sigma_mm83(E):
    s = np.zeros_like(E)
    for lo, hi, c0, c1, c2 in MM83:
        m = (E >= lo) & (E < hi)
        s[m] = (c0 + c1 * E[m] + c2 * E[m]**2) / E[m]**3 * 1e-24
    return s


# ------------------------------------------------------------------ response
with fits.open(f"{P}/spec/psr.arf") as h:
    a = h["SPECRESP"].data
    elo, ehi, area = a["ENERG_LO"], a["ENERG_HI"], a["SPECRESP"]
with fits.open(f"{P}/spec/psr.rmf") as h:
    mat = h["MATRIX"].data
    eb = h["EBOUNDS"].data
    nchan = len(eb)
    R = np.zeros((len(mat), nchan))
    for i, row in enumerate(mat):
        fch = np.atleast_1d(row["F_CHAN"]); nch = np.atleast_1d(row["N_CHAN"])
        m = np.asarray(row["MATRIX"]); k = 0
        for f, n in zip(fch, nch):
            f0 = int(f) - int(eb["CHANNEL"][0])
            R[i, f0:f0 + int(n)] = m[k:k + int(n)]; k += int(n)
    cen = 0.5 * (eb["E_MIN"] + eb["E_MAX"])
assert len(elo) == R.shape[0]
Ec = 0.5 * (elo + ehi); dE = ehi - elo
band_ch = {"soft": (cen >= 0.5) & (cen < 1.5), "medium": (cen >= 1.5) & (cen < 3.0),
           "hard": (cen >= 3.0) & (cen < 8.0)}


def band_counts(gamma, nh=NH):
    ph = Ec**(-gamma) * np.exp(-nh * sigma_mm83(Ec)) * dE      # photons/cm2/s per bin
    cts = (ph * area) @ R
    return {b: cts[m].sum() for b, m in band_ch.items()}


G = np.linspace(0.3, 4.0, 371)
tab = np.array([[band_counts(g)[b] for b in BANDS] for g in G])
S, Mm, H = tab.T
HR_of_G = (H - S) / (H + S)
HR1_of_G = (Mm - S) / (Mm + S)
HR2_of_G = (H - Mm) / (H + Mm)


def invert(hr, curve):
    if not np.isfinite(hr):
        return np.nan
    # curve decreases with Gamma
    return float(np.interp(-hr, -curve, G, left=np.nan, right=np.nan))


def gamma_with_err(v, e, curve=HR_of_G):
    g = invert(v, curve)
    lo, hi = invert(v + e, curve), invert(v - e, curve)
    return g, (hi - g if np.isfinite(hi) else np.nan), (g - lo if np.isfinite(lo) else np.nan)


ER = json.load(open(f"{P}/energy_results.json"))
out = {"nH": NH, "curve": dict(gamma=G.tolist(), HR=HR_of_G.tolist(), HR1=HR1_of_G.tolist(),
                               HR2=HR2_of_G.tolist())}
rows = []
for r in ER["regions"]:
    d = dict(name=r["name"])
    d["gamma"], d["gamma_ep"], d["gamma_em"] = gamma_with_err(r["HR"], r["HR_err"])
    d["gamma_HR1"] = invert(r["HR1"], HR1_of_G)
    d["gamma_HR2"] = invert(r["HR2"], HR2_of_G)
    if "net_HR" in r:
        d["net_gamma"], d["net_gamma_ep"], d["net_gamma_em"] = gamma_with_err(r["net_HR"], r["net_HR_err"])
        d["bkg_gamma"], _, _ = gamma_with_err(r["bkg_HR"], r["bkg_HR_err"])
    rows.append(d)
out["regions"] = rows
prof = {}
for s, pr in ER["radial_profiles"].items():
    prof[s] = [dict(r=x["r"], **dict(zip(["gamma", "ep", "em"], gamma_with_err(x["HR"], x["HR_err"]))))
               for x in pr]
out["radial"] = prof
json.dump(out, open(f"{P}/photon_index.json", "w"), indent=1)

print(f"ObsID {OBSID}: HR(Gamma=1.0, 1.5, 2.0) = "
      + ", ".join(f"{np.interp(g, G, HR_of_G):+.3f}" for g in (1.0, 1.5, 2.0)))
print(f"{'region':<28s} {'Gamma(HR)':>16s} {'G(HR1)':>7s} {'G(HR2)':>7s}   {'net Gamma':>16s} {'surround':>8s}")
for d in rows:
    s = (f"{d['name']:<28s} {d['gamma']:5.2f} +{d['gamma_ep']:.2f}/-{d['gamma_em']:.2f}  "
         f"{d['gamma_HR1']:7.2f} {d['gamma_HR2']:7.2f}")
    if "net_gamma" in d:
        s += f"   {d['net_gamma']:5.2f} +{d['net_gamma_ep']:.2f}/-{d['net_gamma_em']:.2f}  {d['bkg_gamma']:6.2f}"
    print(s)
print()
for s, pr in prof.items():
    print(f"{s:<28s}", " ".join(f"{x['r']:.0f}:{x['gamma']:.2f}" for x in pr))
