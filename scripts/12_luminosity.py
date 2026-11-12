"""PWN X-ray flux and luminosity (0.5-8 keV) from net counts + response.
  python scripts/12_luminosity.py 12075
Net counts: r = 2.5-100" minus the mean level at 90-120" (scaled by area),
converted with an absorbed power law (Gamma from the colours, nH = 3e20)
folded through the ARF/RMF.  Unabsorbed flux, L = 4 pi d^2 F at 287 pc.
"""
import sys, os, json
import numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
import importlib.util
OBSID = sys.argv[1]; P = products(OBSID)
spec = importlib.util.spec_from_file_location("pi", os.path.join(TOP, "scripts", "06_photon_index.py"))
sys.argv = [sys.argv[0], OBSID]
pi = importlib.util.module_from_spec(spec)
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(pi)
img, h = load_image(f"{P}/img/band_broad_thresh.img")
ex, _ = load_image(f"{P}/img/band_broad_thresh.expmap")
E, N = offsets_arcsec(h, img.shape, OBSID); R = np.hypot(E, N)
cov = ex > 0.5 * ex.max()
texp = float(fits.getval(f"{P}/evt2_clean.fits", "EXPOSURE", 1))
src = (R >= 2.5) & (R < 100) & cov
bkg = (R >= 90) & (R < 120) & cov
# counts in source minus outer level per unit exposure*area
rate_b = img[bkg].sum() / ex[bkg].sum()            # counts per (cm2 s) per pixel ... relative
net = img[src].sum() - rate_b * ex[src].sum()
# correct for the vignetting / mask difference between the region and the on-axis ARF:
eff = ex[src].mean() / ex[cov].max()
G = 1.5
GI = json.load(open(f"{P}/photon_index.json"))
G = [r for r in GI["regions"] if r["name"] == "diffuse 15-30\""][0]["gamma"]
Ec, dE, area = pi.Ec, pi.dE, pi.area
ph = Ec**(-G) * np.exp(-pi.NH * pi.sigma_mm83(Ec)) * dE
cts = (ph * area) @ pi.R
cen = pi.cen
band = (cen >= 0.5) & (cen < 8.0)
K = net / texp / eff / cts[band].sum()             # normalisation (photons/cm2/s/keV at 1 keV)
sel = (Ec >= 0.5) & (Ec < 8.0)
F_unabs = K * np.sum(Ec[sel]**(1 - G) * dE[sel]) * 1.602e-9      # erg/cm2/s
d_cm = DIST_PC * PC_CM
L = 4 * np.pi * d_cm**2 * F_unabs
Edot = 6.9e36
out = dict(net_counts=float(net), exposure=texp, gamma=float(G), flux_unabs_0p5_8=float(F_unabs),
           L_0p5_8=float(L), L_over_Edot=float(L / Edot))
json.dump(out, open(f"{P}/luminosity.json", "w"), indent=1)
print(OBSID, {k: (f"{v:.3g}") for k, v in out.items()})
