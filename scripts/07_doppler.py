"""Front/back brightness contrast of the outer ring -> Doppler-boosting flow speed.

  python scripts/07_doppler.py 12075

The fitted ring (circular ring of radius a, inclined so that b = a cos i) is
sampled along its full circumference.  For an equatorial outflow of speed
beta*c the line-of-sight velocity at ring azimuth t is beta sin(i) sin(t)
(t = +90 deg is the near/front side, t = -90 deg the far/back side), and the
surface brightness scales as D^(Gamma+1) with D = 1/[gamma (1 - beta sin i sin t)].
The ratio of front to back brightness, R = [(1 + beta sin i)/(1 - beta sin i)]^(Gamma+1),
gives beta.  Brightness on the ring is measured above the local nebula
(annuli 3-5" inside and outside the ring) and the near-axis parts
(where the jets / pulsar cross the ring) are excluded.
"""
import sys, os, json
import numpy as np
from scipy.ndimage import map_coordinates
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

OBSID = sys.argv[1] if len(sys.argv) > 1 else "12075"
P = products(OBSID)
M = json.load(open(f"{P}/measurements.json"))
ring = M["outer_arc"]
GI = json.load(open(f"{P}/photon_index.json"))
gamma_arc = [r for r in GI["regions"] if r["name"] == "outer arc"][0]["gamma"]

img, hdr = load_image(f"{P}/img/fine_broad.img")
E, N = offsets_arcsec(hdr, img.shape, OBSID)
dEdj = E[0, 1] - E[0, 0]; dNdi = N[1, 0] - N[0, 0]
sm = smooth(img, 0.75 / abs(dNdi))


def samp(e, n):
    return map_coordinates(sm, [(n - N[0, 0]) / dNdi, (e - E[0, 0]) / dEdj], order=1)


ta = np.radians(ring["axis_pa"]); ue, un = np.sin(ta), np.cos(ta); pe, pn = un, -ue
a, b, d = ring["a"], ring["b"], ring["d"]
t = np.radians(np.arange(0, 360, 2.0))


def ellipse(scale):
    e = d * ue + scale * (a * np.cos(t) * pe + b * np.sin(t) * ue)
    n = d * un + scale * (a * np.cos(t) * pn + b * np.sin(t) * un)
    return e, n


on = np.mean([samp(*ellipse(1 + s / a)) for s in np.linspace(-0.75, 0.75, 7)], axis=0)
off_in = np.mean([samp(*ellipse(1 - s / a)) for s in (3.0, 4.0, 5.0)], axis=0)
off_out = np.mean([samp(*ellipse(1 + s / a)) for s in (3.0, 4.0, 5.0)], axis=0)
excess = on - 0.5 * (off_in + off_out)
tdeg = np.degrees(t)
# exclude +-25 deg around the axis crossings (jets, pulsar, counter-jet)
front = (np.sin(t) > np.sin(np.radians(30))) & (np.abs(np.cos(t)) > np.sin(np.radians(25)))
back = (np.sin(t) < -np.sin(np.radians(30))) & (np.abs(np.cos(t)) > np.sin(np.radians(25)))
# Poisson-ish error of the excess per sample from the raw counts in the ring strip
pix_area = abs(dNdi)**2
fe, be = excess[front].mean(), excess[back].mean()
fe_err = excess[front].std() / np.sqrt(front.sum() / 4)      # samples are correlated (~4 per beam)
be_err = excess[back].std() / np.sqrt(back.sum() / 4)
incl = np.radians(ring["inclination_deg"])
res = dict(obsid=OBSID, gamma=gamma_arc, inclination_deg=ring["inclination_deg"],
           front_excess=float(fe), back_excess=float(be), front_err=float(fe_err), back_err=float(be_err))
k = gamma_arc + 1


def beta_from_ratio(R):
    x = R**(1 / k)
    return (x - 1) / (x + 1) / np.sin(incl)


if be > 2 * be_err:
    R = fe / be
    R_err = R * np.hypot(fe_err / fe, be_err / be)
    res.update(ratio=float(R), ratio_err=float(R_err), beta=float(beta_from_ratio(R)),
               beta_range=[float(beta_from_ratio(max(R - R_err, 1.0001))), float(beta_from_ratio(R + R_err))])
else:
    # back side not detected: 2-sigma upper limit on its brightness
    R_lo = fe / (max(be, 0.0) + 2 * be_err)
    res["back_upper_limit_2sigma"] = float(max(be, 0.0) + 2 * be_err)
    res["beta_lower_limit_incl64"] = float((R_lo**(1 / k) - 1) / (R_lo**(1 / k) + 1) / np.sin(np.radians(63.6)))
    res.update(ratio_lower_limit=float(R_lo), beta_lower_limit=float(beta_from_ratio(R_lo)))
res["profile_t"] = tdeg.tolist(); res["profile_excess"] = excess.tolist()
json.dump(res, open(f"{P}/doppler.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if not k.startswith("profile")}, indent=1))
