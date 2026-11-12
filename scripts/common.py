"""Shared helpers for the Vela PWN analysis (Project 7)."""
import os
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIX = 0.492            # ACIS sky pixel, arcsec
DIST_PC = 287.0        # VLBI parallax distance (Dodson et al. 2003)
DIST_ERR = (17.0, 19.0)
PC_CM = 3.0857e18
AU_CM = 1.496e13
C_CM_S = 2.998e10

BANDS = ["soft", "medium", "hard"]
BAND_LABEL = {"soft": "0.5–1.5 keV", "medium": "1.5–3 keV",
              "hard": "3–8 keV", "broad": "0.5–8 keV"}


def products(obsid):
    return os.path.join(TOP, "products", str(obsid))


def pulsar_xy(obsid):
    x, y = np.loadtxt(os.path.join(products(obsid), "pulsar_xy.txt"))
    return float(x), float(y)


def load_image(path):
    """Return (data, header) of a FITS image (primary HDU)."""
    with fits.open(path) as h:
        return h[0].data.astype(float), h[0].header


def physical_grid(hdr, shape):
    """Sky (physical) x, y coordinate of every image pixel centre."""
    # Chandra images carry the image->physical transform as LTM/LTV keywords
    ltm1 = hdr.get("LTM1_1", 1.0); ltm2 = hdr.get("LTM2_2", 1.0)
    ltv1 = hdr.get("LTV1", 0.0); ltv2 = hdr.get("LTV2", 0.0)
    iy, ix = np.indices(shape) + 1.0          # FITS 1-based image pixels
    x = (ix - ltv1) / ltm1
    y = (iy - ltv2) / ltm2
    return x, y


def offsets_arcsec(hdr, shape, obsid):
    """East and North offsets (arcsec) of each pixel from the pulsar.

    Chandra sky coordinates are a tangent plane with +y = North and
    +x = West (RA decreasing), so East offset = -(x - x_psr) * 0.492".
    """
    x, y = physical_grid(hdr, shape)
    x0, y0 = pulsar_xy(obsid)
    return -(x - x0) * PIX, (y - y0) * PIX


def extent_arcsec(hdr, shape, obsid):
    """matplotlib imshow extent [left, right, bottom, top] in offset arcsec
    (East to the left, i.e. left > right)."""
    e, n = offsets_arcsec(hdr, shape, obsid)
    dx = abs(e[0, 1] - e[0, 0]); dy = abs(n[1, 0] - n[0, 0])
    return [e[0, 0] + dx / 2, e[0, -1] - dx / 2, n[0, 0] - dy / 2, n[-1, 0] + dy / 2]


def smooth(img, sigma_pix):
    return gaussian_filter(img, sigma_pix)


def arcsec_to_pc(theta, d=DIST_PC):
    return theta / 206265.0 * d


def arcsec_to_cm(theta, d=DIST_PC):
    return arcsec_to_pc(theta, d) * PC_CM


def arcsec_to_au(theta, d=DIST_PC):
    return theta * d                     # 1" at 1 pc = 1 AU


def light_days(theta, d=DIST_PC):
    return arcsec_to_cm(theta, d) / C_CM_S / 86400.0
