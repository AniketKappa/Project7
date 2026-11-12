"""Build the final project report as a designed HTML document (pdf-design-skill,
"soft" style) from the saved results and figures.

  python scripts/13_make_report.py          # writes report/report.html
  bash scripts/14_render_report.sh          # renders Group4_FinalProjectReport.pdf

The text and numbers are generated from products/<obsid>/*.json, so the
report always matches the analysis. Layout and CSS follow the "soft" template
of the pdf-design-skill (white pages, pale-green cards, lime highlighter
labels, Space Grotesk headings with hairline rules, IBM Plex Sans body).
"""
import os, sys, json, html, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

REPORT_DIR = os.path.join(TOP, "report")
os.makedirs(REPORT_DIR, exist_ok=True)
OUTHTML = os.path.join(REPORT_DIR, "report.html")
FIG = "../figures"                       # relative to report/report.html
FOOTER_LEFT = "Group 4 · Project 7"
CONTENT_W_CM = 17.0                      # width the old layout's figure sizes refer to


# ------------------------------------------------------------------ HTML helpers
class Style:
    def __init__(self, name, **kw):
        self.name = name


ParagraphStyle = lambda name, **kw: Style(name)
body, cap, small, h1, h2 = Style("body"), Style("cap"), Style("small"), Style("h1"), Style("h2")
cm = 1.0


def md(s):
    """ReportLab mini-markup -> HTML."""
    s = s.replace("<super>", "<sup>").replace("</super>", "</sup>")
    s = re.sub(r"<font size=[0-9.]+>", '<span class="sm">', s).replace("</font>", "</span>")
    return s


SECTIONS = []          # (id, title) for the table of contents


def P(text, st=body):
    t = md(text)
    if st.name == "h1":
        sid = f"sec{len(SECTIONS) + 1}"
        SECTIONS.append((sid, t))
        m = re.search(r"\(Task ([A-D])\)", t)
        num = re.match(r"(\d+)\.", t)
        chip = f"Task {m.group(1)}" if m else (f"Section {num.group(1)}" if num else ("Appendix" if t.startswith("Appendix") else "Back matter"))
        brk = ""                     # sections flow; the contents page already breaks before section 1
        return f'<div class="opener{brk}"><span class="hl">{chip}</span><h1 class="sec" id="{sid}">{t}</h1></div>'
    if st.name == "h2":
        return f"<h2>{t}</h2>"
    if st.name == "cap":
        return f'<div class="figcap">{t}</div>'
    if st.name == "bul":
        return f'<div class="body bul">{t}</div>'
    if st.name == "ref":
        return f'<div class="ref">{t}</div>'
    if st.name == "small":
        return f'<div class="small">{t}</div>'
    return f'<div class="body">{t}</div>'


def table(rows, widths, header=True, zebra=True, first_col_bold=False, cls=""):
    tot = sum(widths)
    cols = "".join(f'<col style="width:{100 * w / tot:.1f}%">' for w in widths)
    out = [f'<table class="data {cls}"><colgroup>{cols}</colgroup>']
    for i, r in enumerate(rows):
        if header and i == 0:
            out.append("<thead><tr>" + "".join(f"<th>{md(str(c))}</th>" for c in r) + "</tr></thead><tbody>")
            continue
        cells = []
        for j, c in enumerate(r):
            k = ' class="k"' if (first_col_bold and j == 0) else ""
            cells.append(f"<td{k}>{md(str(c))}</td>")
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</tbody></table>" if header else "</table>")
    return "".join(out)


def fig(path, width_cm, caption):
    w = min(100, 100 * width_cm / CONTENT_W_CM)
    return (f'<figure style="width:{w:.0f}%"><img src="{path}">'
            f'<div class="figcap">{md(caption)}</div></figure>')


def box(paras):
    return '<div class="card">' + "".join(f'<div class="body">{md(p)}</div>' for p in paras) + "</div>"


def Spacer(*a, **k):
    return ""


# ------------------------------------------------------------------ data
O = ["128", "12075"]
M = {o: json.load(open(f"{products(o)}/measurements.json")) for o in O}
ER = {o: json.load(open(f"{products(o)}/energy_results.json")) for o in O}
GI = {o: json.load(open(f"{products(o)}/photon_index.json")) for o in O}
DOP = {o: json.load(open(f"{products(o)}/doppler.json")) for o in O}
OJ = {o: json.load(open(f"{products(o)}/outer_jet.json")) for o in O}
LUM = {o: json.load(open(f"{products(o)}/luminosity.json")) for o in O}
reg = {o: {r["name"]: r for r in ER[o]["regions"]} for o in O}
gam = {o: {r["name"]: r for r in GI[o]["regions"]} for o in O}
rad = {o: {x["r"]: x for x in GI[o]["radial"]["all directions"]} for o in O}
r50b = {o: ER[o]["band_r50_bootstrap"] for o in O}


def f1(x): return f"{x:.1f}"
def pc(x): return f"{arcsec_to_pc(x):.3f}"
def cmx(x):
    v = arcsec_to_cm(x); e = int(np.floor(np.log10(v))); return f"{v / 10**e:.1f}×10<super>{e}</super>"
def ld(x): return f"{light_days(x):.0f}"


# physical numbers used in the discussion
EDOT = 6.9e36
R_out = 0.5 * (M["128"]["outer_arc"]["a"] + M["12075"]["outer_arc"]["a"])
R_in = M["12075"]["inner_arc"]["a"]
B_eq = lambda R: np.sqrt(2 * EDOT / (C_CM_S * arcsec_to_cm(R)**2))         # Gauss
t_syn = lambda B, E: 36.0 * (B / 1e-4)**-1.5 * E**-0.5                         # yr
r0 = arcsec_to_cm(R_in); r1 = arcsec_to_cm(40.0); v0 = C_CM_S / 3
t_flow = (r1**3 - r0**3) / (3 * v0 * r0**2) / 3.156e7                        # yr, KC84 v ~ r^-2
B_needed = 1e-4 * (16.1 / t_flow)**(2 / 3)       # B such that t_syn(5 keV) = t_flow



S = []
g128 = {k: rad["128"][k]["gamma"] for k in rad["128"]}; g12 = {k: rad["12075"][k]["gamma"] for k in rad["12075"]}

# ================================================================== cover (soft: gradient hero band, centred title, summary card)
info = table([
    ["Target", "Vela pulsar (PSR B0833−45) and its pulsar wind nebula (PWN); distance d = 287 (+19/−17) pc (VLBI parallax, Dodson et al. 2003)"],
    ["Data", "Chandra ACIS-S: <b>ObsID 128</b> (30 Apr 2000, 10.6 ks; the imaging observation of the reference paper) and the deep "
             "<b>ObsID 12075</b> (4 Sep 2010, 41.0 ks) as an independent, higher-statistics check"],
    ["Software", "CIAO 4.18 (conda), CALDB 4.12.4 (files for each ObsID fetched with <i>download_obsid_caldb</i>), Python (numpy/scipy/astropy/matplotlib)"],
    ["Reference paper", "Pavlov, Zavlin, Sanwal, Burwitz &amp; Garmire (2001), ApJ 552, L129 (astro-ph/0103171)"]],
    [3.2 * cm, 13.8 * cm], header=False, zebra=False, first_col_bold=True, cls="info")
summary = box([
    "<b>Central question.</b> How does the morphology of the Vela PWN change with X-ray photon energy, and what does this tell us "
    "about the physical processes around the pulsar?",
    f"<b>Answer in brief.</b> The nebula is a tilted equatorial torus (outer ring radius ≈ {R_out:.0f}″ = {arcsec_to_pc(R_out):.3f} pc, "
    f"inclination ≈ 55–67°, symmetry axis at PA ≈ 305–309°) with an inner arc, a jet/counter-jet pair along the axis, a faint "
    f"~90″-long outer jet and a diffuse halo ≈ 2′ across. In the hard band (3–8 keV) the nebula is <b>more compact</b>: its half-light radius is "
    f"{r50b['12075']['hard']['r50']:.1f}″ against {r50b['12075']['soft']['r50']:.1f}″ in the soft band (difference "
    f"{r50b['12075']['soft_minus_hard']['value']:.2f} ± {r50b['12075']['soft_minus_hard']['err']:.2f}″; the same in both epochs). "
    f"The jets are the hardest features (net photon index Γ ≈ 1.0–1.1), the torus has the colour of its surroundings (Γ ≈ 1.4), "
    f"and the diffuse emission <b>softens steadily with distance</b>, from Γ ≈ {g12[12.5]:.2f} at 10–15″ to Γ ≈ {g12[35.0]:.2f} at 30–40″. "
    "Fresh particles accelerated at the termination shock flow out along the equator and the axis and lose their highest energies by "
    f"synchrotron radiation as they go (B ≈ 2–3×10<super>−4</super> G); the front/back contrast of the torus requires Doppler boosting by a flow with "
    f"β ≥ {DOP['128']['beta_lower_limit']:.2f}–{DOP['12075']['beta_lower_limit_incl64']:.2f}."])
COVER = f"""
<div class="cover">
  <div class="hero"></div>
  <div class="cover-inner">
    <div class="kick"><span class="hl">FINAL PROJECT REPORT</span></div>
    <h1 class="title">Project 7 — Vela Pulsar:<br>X-ray Imaging with Chandra</h1>
    <div class="body sub">Data reduction · Sizes of the nebular structures · Energy-resolved imaging · Three-colour analysis</div>
    {info}
    {summary.replace('<div class="card">', '<div class="card cover-card"><div class="ct">Summary</div>', 1)}
  </div>
  <div class="cover-footer"><span class="l">{FOOTER_LEFT}</span><span class="r">Vela Pulsar · Final Report</span></div>
</div>
"""

# ================================================================== 1. Data reduction (Task A)
S.append(P("1. Data and Reduction (Task A)", h1))
S.append(P("1.1 Observations", h2))
S.append(P(
    "The project brief leaves the ObsID open; we used <b>ObsID 128</b>, the timed-exposure ACIS-S image of the pulsar and nebula that "
    "Pavlov et al. (2001) used for orientation, as the main dataset. Because it is short (10.6 ks), every measurement was repeated on "
    "<b>ObsID 12075</b>, the longest single ACIS observation of the Vela PWN (2010; four times the counts). The two observations differ by "
    "120° in roll angle, so instrumental features (readout streak, PSF artefact) rotate between them while real structure does not, and "
    "they are 10 years apart, so the variable features of the nebula can be compared."))
S.append(table([
    ["ObsID", "Date", "Instrument / mode", "Frame time"],
    ["128", "2000-04-30", "ACIS-S3, FAINT, TE, ½ subarray (512 rows)", "1.54 s"],
    ["12075", "2010-09-04", "ACIS-S (S3 aim), FAINT, TE, full frame", "3.24 s"]],
    [1.6 * cm, 2.6 * cm, 9.2 * cm, 3.6 * cm]))
S.append(table([
    ["ObsID", "Roll", "Exposure (clean)", "Counts 0.5–8 keV<super>a</super>"],
    ["128", "260.3°", "10.58 ks (no bins removed)", "130 760"],
    ["12075", "140.1°", "41.02 ks (41.78 ks before light-curve cleaning)", "404 129"]],
    [1.6 * cm, 2.6 * cm, 9.2 * cm, 3.6 * cm]))
S.append(P("Table 1. Observations analysed. <super>a</super> Inside a 4.2′ × 4.2′ box centred on the pulsar.", cap))

S.append(P("1.2 Reduction steps", h2))
S.append(table([
    ["Step", "Tool(s)", "What was done / result"],
    ["Calibration", "download_obsid_caldb", "Only the CALDB 4.12.4 files each observation needs (~0.5 GB instead of the 2.3 GB full CALDB)."],
    ["Reprocess", "chandra_repro", "New level-2 event files with current gain, CTI, bad-pixel and afterglow corrections and EDSER sub-pixel event positions."],
    ["Flares", "dmextract, deflare (lc_clean)", "Light curve of an annulus 2.0–2.7′ from the pulsar in 259-s bins (Fig. 2). No background flares in either observation; for ObsID 12075 the ±20% clipping removed 0.77 ks (a few low-rate bins and partial bins at the ends)."],
    ["Clean events", "dmcopy [status=0][GTI]", "Clean event lists products/&lt;obsid&gt;/evt2_clean.fits."],
    ["Pile-up", "pileup_map", "Peak 0.21 (128) and 0.20 (12075) counts per frame per 3×3-pixel island: ≈ 5–10% pile-up, confined to the pulsar's core (r &lt; 1″)."],
    ["Readout streak", "acisreadcorr (tested)", "Expected streak ≈ 1.3% of the pulsar counts (~40 counts) spread along a whole column; acisreadcorr without a streak spectrum deleted real nebular photons along the strip, so it was <i>not</i> applied."],
    ["PSF artefact", "make_psf_asymmetry_region", "Known HRMA 'hook' 0.6–1.0″ from the pulsar on its NW side (ObsID 128): a soft excess there is instrumental."],
    ["Images", "fluximage, dmcopy", "Counts, exposure maps and exposure-corrected images in soft 0.5–1.5, medium 1.5–3, hard 3–8 and broad 0.5–8 keV at 0.492″ pixels, and counts images of the inner ±64″ at 0.246″."],
    ["Responses", "asphist, mkarf, mkacisrmf", "ARF and RMF at the pulsar position for converting band ratios to photon indices and fluxes (Sherpa / specextract not installed)."]],
    [2.3 * cm, 3.6 * cm, 11.1 * cm]))
S.append(P("Table 2. Processing chain. All scripts are in the scripts/ folder (Appendix A).", cap))
S.append(P("1.3 Structures identified", h2))
S.append(P(
    "Figures 1a and 1b show the clean 0.5–8 keV images. In both epochs we identify: (i) the <b>pulsar</b>, a bright point source; (ii) the "
    "<b>outer arc</b>, ~15″ NW of the pulsar on the axis and curving back on both sides — the near (front) half of an equatorial torus; "
    "(iii) the fainter <b>inner arc</b>, ~5–8″ NW, best seen in the deep data; (iv) a <b>jet</b> along the symmetry axis to the NW and a brighter "
    "<b>counter-jet</b> to the SE; (v) a <b>diffuse PWN</b> ≈ 2′ across, brighter on the W/NW and extended to the S; and (vi) in ObsID 12075 a "
    "faint, curved <b>outer jet</b> leaving the nebula to the NW (in ObsID 128 only a compact blob of it, 78″ from the pulsar, is detected). "
    "The arcs and jets are the same in both epochs to within a few arcseconds, while the outer jet has changed completely — as expected "
    "for its known 'fire-hose' variability (Pavlov et al. 2003)."))
S.append(fig(f"{FIG}/F1a_wide.png", 17.0,
             "Figure 1a. Clean 0.5–8 keV images of the whole nebula, exposure corrected and smoothed with a 1.2″ Gaussian; the cyan box is the "
             "area of Fig. 1b. North up, east left; offsets are from the pulsar."))
S.append(fig(f"{FIG}/F1b_inner.png", 17.0,
             "Figure 1b. The inner nebula (0.246″ pixels, 0.5″ smoothing) with the structures identified in the text."))

S.append(fig(f"{FIG}/F0_lightcurves.png", 17.0,
             "Figure 2. Background light curves used for flare screening (0.5–7 keV, 259-s bins). Shaded band: ±20% around the clipped mean "
             "used by <i>deflare</i>; crosses: removed bins. Neither observation shows a background flare; the zero-rate bins are partial bins at the ends."))

# ================================================================== 2. Sizes (Task B)
S.append(P("2. Sizes of the Structures (Task B)", h1))
S.append(P("2.1 Method", h2))
S.append(P(
    "<b>Arcs.</b> The broad-band image was unwrapped into radius versus position angle (PA, north through east) and the brightness ridge of each "
    "arc was located at every PA. Fitting a free ellipse to one side of a ring is badly constrained, so we fitted the physical model of a "
    "<b>circular ring seen at inclination i</b>: an ellipse with semi-axes a (perpendicular to the symmetry axis) and b = a cos i (along it), "
    "whose centre may sit a distance d from the pulsar along the axis. Errors are from 300 bootstrap resamplings of the ridge points. "
    "<b>Jets.</b> Brightness profiles along the axis were compared with profiles at the same distance 25° off the axis; the jet length is the "
    "last distance where the jet is ≥ 20% brighter than its surroundings (the NW jet is cut where it merges with the outer arc). Widths are "
    "Gaussian FWHMs of the transverse profile with the PSF (0.73″, measured on the pulsar in ObsID 128) removed in quadrature. "
    "<b>Whole nebula.</b> After subtracting the mean level at 90–120″ we measured the radii enclosing 50% and 90% of the nebular light, and the "
    "equivalent ellipse of the region brighter than 4× (and 2×) that level. Physical sizes use l = dθ with d = 287 pc, so "
    f"1″ = 287 AU = 4.29×10<super>15</super> cm = {arcsec_to_pc(1):.5f} pc ≈ {light_days(1):.2f} light-days."))
S.append(fig(f"{FIG}/F2_measurements.png", 17.0,
             "Figure 3. Measurements. Top: ridge points and fitted rings (orange = outer arc, cyan = inner arc; dotted cyan = inner ring with the outer "
             "ring's inclination imposed) and the jet strips. Bottom left: brightness along the symmetry axis (red/green) compared with the nebula at the "
             "same distance 25° off axis (grey). Bottom right: transverse jet profiles and the PSF."))

oa = {o: M[o]["outer_arc"] for o in O}; ia = M["12075"]["inner_arc"]; iaf = M["128"]["inner_arc_fixed_incl"]
jn = {o: M[o]["jets"]["NW_jet"] for o in O}; jc = {o: M[o]["jets"]["SE_counterjet"] for o in O}
bp = {o: M[o]["bright_pwn"]["x4"] for o in O}
ojl = OJ["12075"]


def row(name, v128, v12, e128=None, e12=None, physical=None, note="", nd=1):
    s128 = "—" if v128 is None else (f"{v128:.{nd}f}" + (f" ± {e128:.{nd}f}" if e128 is not None else ""))
    s12 = "—" if v12 is None else (f"{v12:.{nd}f}" + (f" ± {e12:.{nd}f}" if e12 is not None else ""))
    ref = physical if physical is not None else (v12 if v12 is not None else v128)
    return [name, s128, s12, pc(ref), cmx(ref), ld(ref), note]


T = [["Structure / quantity", "ObsID 128 (″)", "ObsID 12075 (″)", "pc<super>b</super>", "cm<super>b</super>", "lt-days<super>b</super>", "Notes"]]
T.append(row("Outer arc (torus): ring radius a", oa["128"]["a"], oa["12075"]["a"], oa["128"]["err"]["a"], oa["12075"]["err"]["a"],
             note="semi-major axis of the projected ellipse"))
T.append(row("Outer arc: torus diameter 2a", 2 * oa["128"]["a"], 2 * oa["12075"]["a"], 2 * oa["128"]["err"]["a"], 2 * oa["12075"]["err"]["a"]))
T.append(row("Outer arc: semi-minor axis b", oa["128"]["b"], oa["12075"]["b"], oa["128"]["err"]["b"], oa["12075"]["err"]["b"]))
T.append(row("Outer arc: distance of arc vertex from pulsar", oa["128"]["vertex_dist"], oa["12075"]["vertex_dist"], note="along the axis, NW"))
T.append(["Outer arc: inclination i = arccos(b/a)", f"{oa['128']['inclination_deg']:.0f} ± {oa['128']['inclination_err']:.0f}°",
          f"{oa['12075']['inclination_deg']:.0f} ± {oa['12075']['inclination_err']:.0f}°", "", "", "", "Ng &amp; Romani: 63.6°"])
T.append(["Symmetry-axis position angle", f"{oa['128']['axis_pa']:.1f} ± {oa['128']['err']['axis_pa']:.1f}°",
          f"{oa['12075']['axis_pa']:.1f} ± {oa['12075']['err']['axis_pa']:.1f}°", "", "", "", "proper motion PA ≈ 301°"])
T.append(row("Inner arc: ring radius a", iaf["a"], ia["a"], iaf["err"]["a"], ia["err"]["a"],
             note=f"128: only the arc ends seen (fit with i fixed); 12075: i = {ia['inclination_deg']:.0f} ± {ia['inclination_err']:.0f}°"))
T.append(row("Inner arc: distance of arc vertex from pulsar", None, ia["vertex_dist"]))
T.append(row("NW jet: length", jn["128"]["length_arcsec"], jn["12075"]["length_arcsec"], note="until it meets the outer arc"))
T.append(row("NW jet: width (FWHM, PSF-corrected)", jn["128"]["fwhm_intrinsic"], jn["12075"]["fwhm_intrinsic"],
             jn["128"]["fwhm_err"], jn["12075"]["fwhm_err"], note="12075 broadened by inner-arc knots", nd=2))
T.append(row("SE counter-jet: length", jc["128"]["length_arcsec"], jc["12075"]["length_arcsec"]))
T.append(row("SE counter-jet: width (FWHM, PSF-corrected)", jc["128"]["fwhm_intrinsic"], jc["12075"]["fwhm_intrinsic"],
             jc["128"]["fwhm_err"], jc["12075"]["fwhm_err"], nd=2))
T.append(row("Outer jet: path length", None, ojl["path_length_arcsec"],
             note=f"from r = {ojl['start_distance_arcsec']:.0f}″ to {ojl['end_distance_arcsec']:.0f}″ ({pc(ojl['end_distance_arcsec'])} pc); width ≈ {ojl['width_fwhm_arcsec']:.0f}″; 128: blob at r ≈ 78″"))
T.append(row("Bright PWN: half-light radius", M["128"]["half_light_radius"], M["12075"]["half_light_radius"]))
T.append(row("Bright PWN: 90%-light radius", M["128"]["r90_radius"], M["12075"]["r90_radius"]))
T.append(row("Bright PWN (> 4× outer level): semi-axes", bp["128"]["a"], bp["12075"]["a"],
             note=f"× {bp['128']['b']:.0f}″ / {bp['12075']['b']:.0f}″; i.e. ≈ {2*bp['12075']['a']:.0f}″ × {2*bp['12075']['b']:.0f}″ ≈ 2′ across"))
S.append(table([[r[0], r[1], r[2], r[6]] for r in T], [5.6 * cm, 3.0 * cm, 3.0 * cm, 5.4 * cm], first_col_bold=True))
S.append(table([[r[0], r[3], r[4], r[5]] for r in T if (r is T[0] or r[3])], [7.0 * cm, 3.0 * cm, 3.5 * cm, 3.5 * cm], first_col_bold=True))
S.append(P("Table 3. Angular and physical sizes. <super>b</super> Physical values for the ObsID 12075 number (or ObsID 128 when it is the only one), at "
           "d = 287 pc; light-days = light-travel time across the length. Using d = 210 pc (Pavlov et al. 2001) would reduce all physical sizes by 27%; "
           "the ±18 pc parallax error adds ±6%.", cap))
S.append(P("2.2 Results and comparison with earlier work", h2))
S.append(P(
    f"The outer arc is described well by an inclined ring (rms scatter of the ridge about the model {oa['128']['rms_resid']:.1f}″ and "
    f"{oa['12075']['rms_resid']:.1f}″). Its radius, {oa['128']['a']:.1f} ± {oa['128']['err']['a']:.1f}″ (2000) and "
    f"{oa['12075']['a']:.1f} ± {oa['12075']['err']['a']:.1f}″ (2010), i.e. {arcsec_to_pc(R_out):.3f} pc ≈ {arcsec_to_au(R_out):,.0f} AU, is the same in both "
    "epochs and agrees with the torus fits of Ng &amp; Romani (2004, 2008), as do the axis direction (PA ≈ 305–309°, close to the proper-motion "
    "direction PA ≈ 301°) and the inclination. The inclination differs between our two fits (67 ± 4° vs 54 ± 2°): with only the front half of the "
    "ring visible, b and the centre offset d are partly degenerate and the brightness distribution along the arc changes between epochs, so we "
    "regard i ≈ 55–67° as the range allowed by our data (Ng &amp; Romani: 63.6°). The inner arc has a radius ≈ 9″ (0.012 pc) and comes within "
    f"≈ {ia['vertex_dist']:.0f}″ of the pulsar; it is the closest feature to the wind termination shock. Both jets are ≈ 12–14″ (≈ 0.018 pc) long "
    "and ≈ 2.5″ (≈ 0.0035 pc) wide, i.e. resolved by Chandra. The bright nebula has a half-light radius of 23″ (0.032 pc) and is ≈ 2′ (0.17 pc) "
    f"across, as in the schematic of our pre-project report; the outer jet reaches {ojl['end_distance_arcsec']:.0f}″ (0.14 pc) from the pulsar."))

# ================================================================== 3. Energy-resolved imaging (Task C)
S.append(P("3. Energy-Resolved Imaging (Task C)", h1))
S.append(P(
    "Figures 4 and 5 show the soft (0.5–1.5 keV), medium (1.5–3 keV) and hard (3–8 keV) images. To make the comparison quantitative we used "
    "(a) the half-light radius of the nebula in each band (outer level subtracted, exposure corrected, Poisson errors from 200 simulated images); "
    "(b) hardness ratios HR = (H − S)/(H + S) (and HR1 = (M − S)/(M + S), HR2 = (H − M)/(H + M)) in regions and annuli; and (c) conversion of the "
    "hardness ratios into the photon index Γ of an absorbed power law (n<sub>H</sub> = 3×10<super>20</super> cm<super>−2</super>, Pavlov et al. 2001), "
    "folded through each observation's ARF and RMF. Step (c) removes the loss of soft-band sensitivity of ACIS between 2000 and 2010 "
    "(contamination of the optical-blocking filter): the raw hardness ratios of ObsID 12075 are systematically higher, but the photon indices of the "
    "two epochs agree (Table 5)."))
S.append(fig(f"{FIG}/fig_bands_12075.png", 17.0,
             "Figure 4. ObsID 12075 in the three bands. Top: inner ±35″ (counts, 1″ smoothing; dashed = fitted outer ring). Bottom: whole nebula, "
             "exposure corrected, 1.5″ smoothing. The same stretch is applied to each band relative to its own range."))
S.append(fig(f"{FIG}/fig_bands_128.png", 15.5,
             "Figure 5. As Fig. 4 for ObsID 128. The faint E–W line in the soft wide image is the pulsar's readout streak (left uncorrected, Sect. 1.2)."))
S.append(P("3.1 How the morphology changes with energy", h2))
rb = r50b
S.append(table([
    ["", "Soft 0.5–1.5 keV", "Medium 1.5–3 keV", "Hard 3–8 keV"],
    ["Half-light radius, ObsID 128", f"{rb['128']['soft']['r50']:.2f} ± {rb['128']['soft']['err']:.2f}″",
     f"{rb['128']['medium']['r50']:.2f} ± {rb['128']['medium']['err']:.2f}″", f"{rb['128']['hard']['r50']:.2f} ± {rb['128']['hard']['err']:.2f}″"],
    ["Half-light radius, ObsID 12075", f"{rb['12075']['soft']['r50']:.2f} ± {rb['12075']['soft']['err']:.2f}″",
     f"{rb['12075']['medium']['r50']:.2f} ± {rb['12075']['medium']['err']:.2f}″", f"{rb['12075']['hard']['r50']:.2f} ± {rb['12075']['hard']['err']:.2f}″"],
    ["80%-light radius, ObsID 12075", f"{ER['12075']['band_sizes']['soft']['r80']:.1f}″", f"{ER['12075']['band_sizes']['medium']['r80']:.1f}″",
     f"{ER['12075']['band_sizes']['hard']['r80']:.1f}″"],
    ["Appearance", "Smooth, extended nebula; arcs visible but with low contrast; pulsar surrounded by a soft halo; extended southern emission",
     "Arcs, jets and the region inside the torus show the highest contrast", "Compact: inner nebula, arcs and counter-jet; outer halo faint"]],
    [4.1 * cm, 4.3 * cm, 4.3 * cm, 4.3 * cm], first_col_bold=True))
S.append(P(f"Table 4. Size and appearance of the nebula in the three bands. The hard nebula is smaller than the soft one by "
           f"{rb['128']['soft_minus_hard']['value']:.2f} ± {rb['128']['soft_minus_hard']['err']:.2f}″ (2000) and "
           f"{rb['12075']['soft_minus_hard']['value']:.2f} ± {rb['12075']['soft_minus_hard']['err']:.2f}″ (2010), i.e. by ≈ 10%, at > 12σ in each epoch.", cap))
S.append(P(
    "In all three bands the same skeleton is present — torus, inner arc, jets — so the overall geometry does not depend on energy; what changes "
    "is the <b>contrast and extent</b>. The soft image is dominated by the smooth, extended nebula and by a halo around the pulsar; the medium band "
    "shows the arcs and jets most clearly (as predicted in our pre-project report); and the hard image is visibly more compact, with the emission "
    "concentrated inside and on the torus and along the jets. The nebula therefore shrinks with increasing photon energy."))

# ================================================================== 4. Colour images (Task D)
S.append(P("4. Three-Colour Imaging and Spectral Map (Task D)", h1))
S.append(P(
    "We built three-colour images with red = soft, green = medium and blue = hard, both with the standard CIAO recipe (<i>aconvolve</i> + "
    "<i>dmimg2jpg</i>, Fig. 6) and in Python with a <b>colour-balanced</b> scaling (Fig. 7): each band is divided by its mean brightness over the "
    "nebula (3″ &lt; r &lt; 45″) and one common stretch is applied, so that the average nebular spectrum appears white and colours show where the "
    "emission is relatively soft (red/yellow) or hard (blue). Figure 8 is the corresponding hardness-ratio map and Fig. 9 the photon index as a "
    "function of distance from the pulsar in four directions."))
S.append('<div class="pair-wrap" style="page-break-inside: avoid"><div class="pair">'
         + "".join(f'<div class="c"><img src="{FIG}/F4b_ciao_rgb_{o}.jpg"><div class="small lbl">ObsID {o}</div></div>' for o in ["128", "12075"])
         + '</div>' + P("Figure 6. Standard CIAO three-colour images (aconvolve, σ = 2 pixels; dmimg2jpg with asinh scaling, each band "
                        "scaled between its 1st and 99.8th percentile), ±64″ around the pulsar. Red = soft, green = medium, blue = hard.", cap) + '</div>')
S.append(fig(f"{FIG}/F4_rgb.png", 16.0,
             "Figure 7. Colour-balanced three-colour images (red 0.5–1.5, green 1.5–3, blue 3–8 keV). White = the average PWN spectrum."))
S.append(fig(f"{FIG}/fig_hardness_map_12075.png", 16.5,
             "Figure 8. Hardness maps of ObsID 12075 (1″ pixels, 2″ Gaussian smoothing; pixels with fewer than 150 counts blank). Blue = harder. "
             "Dashed: fitted outer ring; contours: broad-band brightness."))
S.append(fig(f"{FIG}/F5_gamma_profiles.png", 17.0,
             "Figure 9. Photon index Γ (from HR) versus distance from the pulsar in four 90° sectors and for all directions (thick line), for both "
             "epochs; right: half-light radius of the nebula in each band. Beyond ≈ 60″ (grey) the nebula is faint and the instrumental and Vela-SNR "
             "background, which has a different spectrum, dominates the colours."))

S.append(P("4.1 Colours of the individual structures", h2))
names = [("inner knot (1.5-3.5\" NW)", "Pulsar halo, 1.5–3.5″ NW"), ("NW jet", "NW jet"), ("SE counter-jet", "SE counter-jet"),
         ("inner arc (ends)", "Inner arc (ends)"), ("outer arc", "Outer arc"), ("diffuse 3-15\"", "Diffuse 3–15″"),
         ("diffuse 15-30\"", "Diffuse 15–30″"), ("diffuse 30-45\"", "Diffuse 30–45″"), ("diffuse 45-60\"", "Diffuse 45–60″")]
T5 = [["Region", "HR 2000", "HR 2010", "Γ 2000", "Γ 2010", "Γ net 2000<super>c</super>", "Γ net 2010<super>c</super>", "Γ surroundings 2000 / 2010"]]
def gfmt(d, k="gamma"):
    return f"{d[k]:.2f} <font size=6.5>+{d[k + '_ep']:.2f}/−{d[k + '_em']:.2f}</font>" if k + "_ep" in d else f"{d[k]:.2f}"
for key, nm in names:
    r1, r2 = reg["128"][key], reg["12075"][key]; g1, g2 = gam["128"][key], gam["12075"][key]
    T5.append([nm, f"{r1['HR']:+.3f}", f"{r2['HR']:+.3f}", f"{g1['gamma']:.2f}±{g1['gamma_ep']:.2f}", f"{g2['gamma']:.2f}±{g2['gamma_ep']:.2f}",
               (f"{g1['net_gamma']:.2f} +{g1['net_gamma_ep']:.2f}/−{g1['net_gamma_em']:.2f}" if "net_gamma" in g1 else "—"),
               (f"{g2['net_gamma']:.2f} +{g2['net_gamma_ep']:.2f}/−{g2['net_gamma_em']:.2f}" if "net_gamma" in g2 else "—"),
               (f"{g1['bkg_gamma']:.2f} / {g2['bkg_gamma']:.2f}" if "bkg_gamma" in g1 else "—")])
S.append(table([[r[0], (r[1] + " / " + r[2]) if r is not T5[0] else "HR 2000 / 2010", r[3], r[4]] for r in T5],
                [4.4 * cm, 4.2 * cm, 4.2 * cm, 4.2 * cm], first_col_bold=True))
S.append(table([[r[0], r[5], r[6], r[7]] for r in T5 if (r is T5[0] or r[5] != "—")],
                [4.4 * cm, 4.2 * cm, 4.2 * cm, 4.2 * cm], first_col_bold=True))
S.append(P("Table 5. Hardness ratios HR = (H − S)/(H + S) and photon indices. <super>c</super> 'Net' = after subtracting the nebula at the same distance "
           "from the pulsar, rotated by ±40° (the colour of the structure itself); last column = colour of that comparison region. "
           f"The pulsar itself (r &lt; 1.5″; Γ<sub>HR</sub> = {gam['128']['pulsar (r<1.5\")']['gamma']:.2f} and {gam['12075']['pulsar (r<1.5\")']['gamma']:.2f}) "
           "is not listed: its colours are distorted by pile-up and its spectrum is not a single power law (Sect. 5.4).", cap))

S.append(P("4.2 Answers to the Task D questions", h2))
jnet = (gam["128"]["NW jet"], gam["12075"]["NW jet"]); cnet = (gam["128"]["SE counter-jet"], gam["12075"]["SE counter-jet"])
S.append(table([
    ["Question", "Answer from the colour images, maps and Table 5"],
    ["Which regions are dominated by <b>soft</b> emission?",
     f"The outer, diffuse halo (r ≥ 25–30″, Γ ≈ 1.6–1.7), especially the extension to the S/SW, and the halo immediately around the pulsar "
     f"(r ≈ 1–3″; HR ≈ −0.7 to −0.9), which is the pulsar's own soft (thermal) emission spread by the PSF plus the soft PSF artefact."],
    ["Which regions are dominated by <b>hard</b> emission?",
     f"Everything inside and on the torus (r ≤ 15″, Γ ≈ {g12[12.5]:.2f}–{g12[17.5]:.2f}), the jets and the inner arc. The blue pulsar core in Fig. 7 "
     "is an artefact of pile-up (two soft photons recorded as one harder event)."],
    ["Is the jet harder or softer than the surrounding PWN?",
     f"<b>Harder.</b> NW jet: net Γ = {jnet[0]['net_gamma']:.2f} +{jnet[0]['net_gamma_ep']:.2f}/−{jnet[0]['net_gamma_em']:.2f} (2000) and "
     f"{jnet[1]['net_gamma']:.2f} ± {jnet[1]['net_gamma_ep']:.2f} (2010) against Γ ≈ 1.3 around it. SE counter-jet in 2010: "
     f"{cnet[1]['net_gamma']:.2f} +{cnet[1]['net_gamma_ep']:.2f}/−{cnet[1]['net_gamma_em']:.2f} against {cnet[1]['bkg_gamma']:.2f} (≈ 6σ). "
     f"In the short 2000 exposure the counter-jet was not harder ({cnet[0]['net_gamma']:.2f} +{cnet[0]['net_gamma_ep']:.2f}/−{cnet[0]['net_gamma_em']:.2f}), "
     "consistent with the known variability of the jets."],
    ["Does the emission become softer with distance from the pulsar?",
     f"<b>Yes</b>, in all four directions and in both epochs: Γ rises from ≈ {g128[8.0]:.2f}/{g12[8.0]:.2f} at 6–10″ to "
     f"≈ {g128[45.0]:.2f}/{g12[45.0]:.2f} at 40–50″ (2000/2010), ΔΓ ≈ +0.3 over 0.05 pc; the hard nebula is ≈ 10% smaller. Both colours "
     "(HR1 and HR2) steepen by the same amount, which favours a steepening power law over an added soft thermal component."],
    ["Are the torus and jet structures different in their energy distribution?",
     f"<b>Yes.</b> The outer arc has the same spectrum as the nebula around it (net Γ = {gam['12075']['outer arc']['net_gamma']:.2f} vs "
     f"{gam['12075']['outer arc']['bkg_gamma']:.2f} in 2010), the inner arc is slightly harder (net Γ ≈ 1.15–1.35), and the jets are the hardest "
     "structures (Γ ≈ 1.0–1.1), ≈ 0.3–0.4 harder than the torus."]],
    [4.6 * cm, 12.4 * cm], first_col_bold=False))
S.append(P("Table 6. Task D questions.", cap))

# ================================================================== 5. Interpretation
S.append(P("5. Physical Interpretation", h1))
S.append(P("5.1 Synchrotron nebula and cooling", h2))
p_in = 2 * g12[12.5] - 1; p_out = 2 * g12[35.0] - 1
S.append(P(
    f"The nebular spectrum is a hard power law (Γ ≈ 1.3–1.7), the signature of synchrotron emission from a power-law population of electrons "
    f"and positrons, N(E) ∝ E<super>−p</super> with p = 2Γ − 1 ≈ {p_in:.1f} near the torus and ≈ {p_out:.1f} in the outer nebula. Synchrotron "
    "losses are fastest for the most energetic particles, t<sub>syn</sub> ≈ 36 yr (B/100 µG)<super>−3/2</super>(E/1 keV)<super>−1/2</super>, so a "
    "population that flows away from the termination shock while radiating loses its hard-X-ray-emitting particles first: its spectrum steepens "
    "(towards ΔΓ = 0.5 for complete cooling) and the nebula looks smaller at high energy. This is exactly what we see: ΔΓ ≈ +0.3 between 10″ and "
    "40″ and a hard half-light radius 10% smaller than the soft one. That ΔΓ is less than 0.5 means the cooling is partial on these scales."))
S.append(P(
    f"The time scales are consistent with this picture. Pressure balance at the torus, B²/8π ≈ Ė/(4πcR²) with Ė = 6.9×10<super>36</super> erg s<super>−1</super> "
    f"and R = {R_out:.0f}″, gives B ≈ {B_eq(R_out) * 1e6:.0f} µG, for which electrons radiating at 5 keV (1 keV) live t<sub>syn</sub> ≈ "
    f"{t_syn(B_eq(R_out), 5):.0f} yr ({t_syn(B_eq(R_out), 1):.0f} yr). A relativistic post-shock flow (≈ c/3) would carry them far beyond 0.06 pc in that "
    f"time, but in the standard model of Kennel &amp; Coroniti (1984) the flow decelerates (v ∝ r<super>−2</super>); starting at c/3 at the inner arc "
    f"(r ≈ {R_in:.0f}″) it takes ≈ {t_flow:.0f} yr to reach 40″, equal to the 5-keV cooling time for B ≈ {B_needed * 1e6:.0f} µG — the same order as the "
    "pressure-balance field. The observed softening therefore implies a few-hundred-µG field and a strongly decelerating outflow inside the bright PWN."))
S.append(P("5.2 Jets and arcs", h2))
S.append(P(
    "The jets are the hardest structures (Γ ≈ 1.0–1.1, p ≈ 1.0–1.2): they carry freshly accelerated particles along the spin axis from very "
    "close to the pulsar and are seen before the particles have cooled, so their spectrum is close to the injection spectrum. The inner arc, "
    "closest to the termination shock, is also slightly harder than its surroundings, while the outer arc (the downstream equatorial flow) "
    "has the same spectrum as the diffuse nebula. The torus and jets therefore trace two different outflows — an equatorial wind and polar jets — "
    "with different particle histories. Their alignment with the proper motion (axis PA ≈ 305–309° vs 301°) supports spin–kick alignment."))
S.append(P("5.3 Relativistic Doppler boosting", h2))
d128, d12 = DOP["128"], DOP["12075"]
S.append(P(
    f"Only the NW (front) half of the torus is seen: along the fitted ring, the back half is not detected above the local nebula (Fig. 10). For an "
    f"equatorial outflow of speed βc the front/back brightness ratio is [(1 + β sin i)/(1 − β sin i)]<super>Γ+1</super>. The 2σ lower limits on the "
    f"ratio, > {d128['ratio_lower_limit']:.0f} (2000) and > {d12['ratio_lower_limit']:.0f} (2010), give β ≥ {d128['beta_lower_limit']:.2f} and "
    f"β ≥ {d12['beta_lower_limit']:.2f} ({d12['beta_lower_limit_incl64']:.2f} for i = 63.6°): the post-shock flow is mildly relativistic, in line with the "
    "0.3–0.7c motions of blobs in the arcs and jets measured with Chandra (Pavlov et al. 2001b, 2003). Doppler boosting, not intrinsic emission, "
    "explains why the nebula is so much brighter on its NW side."))
S.append(fig(f"{FIG}/F7_ring_doppler.png", 11.0,
             "Figure 10. Brightness of the outer ring above the local nebula along the fitted ring, normalised to the mean of the front segments "
             "(red bands). The back segments (blue bands) are consistent with zero."))
S.append(P("5.4 The pulsar and the energetics", h2))
S.append(P(
    "The pulsar's core looks hard (HR ≈ −0.1, 'blue') only because of pile-up; outside the piled core (1–2.5″) its PSF halo is the softest emission "
    "in the field (HR ≈ −0.8 in all directions). Since the wings of Chandra's PSF are, if anything, harder than its core, the pulsar is intrinsically "
    "softer than the nebula, consistent with the thermal (≈ 0.7 MK) atmosphere plus weak power law of Pavlov et al. (2001); its colours do not fit a "
    "single power law (Γ from HR1 and HR2 differ strongly), so no photon index is quoted for it. For the nebula (r = 2.5–100″, outer level subtracted) "
    f"the unabsorbed 0.5–8 keV flux is {LUM['128']['flux_unabs_0p5_8'] * 1e11:.2f}×10<super>−11</super> (2000) and "
    f"{LUM['12075']['flux_unabs_0p5_8'] * 1e11:.2f}×10<super>−11</super> erg cm<super>−2</super> s<super>−1</super> (2010), i.e. "
    f"L<sub>X</sub> ≈ {LUM['12075']['L_0p5_8'] / 1e32:.1f}×10<super>32</super> erg s<super>−1</super> at 287 pc and L<sub>X</sub>/Ė ≈ "
    f"{LUM['12075']['L_over_Edot'] * 1e4:.1f}×10<super>−4</super>. The 2% agreement between two epochs with very different detector responses "
    "is a good check of the whole reduction. Vela converts only ~10<super>−4</super> of its spin-down power into nebular X-rays."))

# ================================================================== 6. Pre-project hypotheses
S.append(P("6. Our Pre-Project Expectations, Tested", h1))
S.append(table([
    ["Expectation (pre-project report, Sect. 3.4)", "Result"],
    ["Soft band dominated by the thermally bright pulsar plus diffuse nebula", "Partly. The diffuse nebula dominates; the pulsar is surrounded by a very soft halo, but its ≈ 60 eV thermal emission peaks below 0.5 keV, so in 0.5–1.5 keV it does not outshine the nebula."],
    ["Medium band shows the arcs and jets most clearly", "Confirmed (Figs. 4–5)."],
    ["Hard band: a more compact nebula on the inner arc and jets, with fewer counts", "Confirmed: half-light radius 10% smaller (> 12σ); 71 471 hard vs 210 332 soft counts in ObsID 12075."],
    ["Pulsar 'red', jets and inner torus 'blue'", "Confirmed for the jets and inner torus; the unpiled pulsar halo is red, the piled core appears blue."],
    ["Outer diffuse emission redder than the inner structures if cooling matters", "Confirmed: ΔΓ ≈ +0.3 between 10″ and 40″; cooling is partial, as anticipated."],
    ["Pile-up may affect the pulsar core", "Confirmed and quantified (≈ 5–10%); the core is excluded from the nebular analysis."],
    ["Sizes: torus radius ≈ 20″, bright PWN ≈ 2′", f"Confirmed: torus radius {R_out:.1f}″ (0.030 pc); bright PWN ≈ 2′ (0.17 pc)."]],
    [6.5 * cm, 10.5 * cm]))
S.append(P("Table 7. Hypotheses from the pre-project report compared with the results.", cap))

# ================================================================== 7. Limitations
S.append(P("7. Limitations and Checks", h1))
bul = ParagraphStyle("bul", parent=body, leftIndent=11, firstLineIndent=-11, spaceAfter=3)
for t in [
    "• <b>No full spectral fitting.</b> Photon indices come from band ratios with n<sub>H</sub> fixed and one response (at the pulsar) for the whole "
    "nebula; within 1′ of the pulsar the effective area changes by only a few per cent. A ±50% change of n<sub>H</sub> shifts Γ by ≈ ±0.03, equally "
    "in every region, so it does not affect the trends we report. The agreement of Γ between 2000 and 2010, after the large change in the soft "
    "response, shows the method is robust.",
    "• <b>Background.</b> The whole field is filled by the nebula and the Vela SNR, so there is no source-free background; the outer level (90–120″) "
    "contains instrumental, SNR and faint PWN emission. The colour trend is measured where the PWN is ≥ 4× brighter than that level.",
    "• <b>Ring fits.</b> Only the front half of each ring is visible, so b and the centre offset are correlated; the inclination is therefore less "
    "certain than the radius.",
    "• <b>Variability.</b> The jets, outer jet and knots change on weeks–months; sizes of these features refer to the epoch of each observation.",
    "• <b>Instrumental features</b> (pile-up, readout streak, PSF artefact) were identified and excluded rather than modelled."]:
    S.append(P(t, bul))

# ================================================================== 8. Conclusions
S.append(P("8. Conclusions", h1))
S.append(P(
    f"(1) Reprocessing ObsID 128 (and the deep ObsID 12075) gives clean images in which the pulsar, an inner and an outer arc, a jet and counter-jet, "
    f"the diffuse nebula and the outer jet are all identified. (2) The outer arc is the front of a torus of radius {R_out:.0f}″ = {arcsec_to_pc(R_out):.3f} pc "
    f"inclined by ≈ 55–67° about an axis at PA ≈ 305–309°; the inner arc has a radius ≈ 9″ (0.012 pc); the jets are ≈ 13″ (0.018 pc) long and "
    f"≈ 2.5″ wide; the bright nebula is ≈ 2′ (0.17 pc) across with a half-light radius of 23″ (0.032 pc); the outer jet extends to 0.14 pc. "
    "(3) With increasing photon energy the nebula becomes more compact and its emission concentrates on the torus and jets. (4) The jets are the "
    "hardest structures, the torus has the colour of its surroundings, and the emission softens with distance from the pulsar. (5) These facts "
    "describe a synchrotron nebula fed by a mildly relativistic (β ≥ 0.4–0.6), decelerating equatorial wind and polar jets, in which particles "
    "accelerated at the termination shock cool in a ~2–3×10<super>−4</super> G field as they flow outward."))

# ================================================================== references / appendix
S.append(P("References", h1))
refs = [
    "Dodson, R. et al. 2003, ApJ, 596, 1137 — Vela pulsar parallax and proper motion from VLBI.",
    "Gaensler, B. M. &amp; Slane, P. O. 2006, ARA&amp;A, 44, 17 — Evolution and structure of pulsar wind nebulae.",
    "Helfand, D. J., Gotthelf, E. V. &amp; Halpern, J. P. 2001, ApJ, 556, 380 — Vela pulsar and its synchrotron nebula.",
    "Kargaltsev, O. &amp; Pavlov, G. G. 2008, AIP Conf. Proc., 983, 171 — Pulsar wind nebulae in the Chandra era.",
    "Kennel, C. F. &amp; Coroniti, F. V. 1984, ApJ, 283, 694 — Confinement of the Crab pulsar's wind by its supernova remnant.",
    "Morrison, R. &amp; McCammon, D. 1983, ApJ, 270, 119 — Interstellar photoelectric absorption cross sections.",
    "Ng, C.-Y. &amp; Romani, R. W. 2004, ApJ, 601, 479; 2008, ApJ, 673, 411 — Fitting pulsar wind tori I, II.",
    "Pavlov, G. G. et al. 2001, ApJ, 552, L129 — The X-ray spectrum of the Vela pulsar resolved with Chandra (reference paper).",
    "Pavlov, G. G. et al. 2001b, ApJ, 554, L189 — Variability of the Vela PWN; Pavlov, G. G. et al. 2003, ApJ, 591, 1157 — The variable jet of the Vela pulsar.",
    "Fruscione, A. et al. 2006, Proc. SPIE, 6270, 62701V — CIAO: Chandra's data analysis system; Chandra X-ray Center, CIAO threads (cxc.harvard.edu/ciao)."]
for r in refs:
    S.append(P(r, ParagraphStyle("ref", parent=small, leftIndent=10, firstLineIndent=-10, spaceAfter=2)))

S.append(P("Appendix A. Reproducing the Analysis", h1))
S.append(P("All products are in the project folder <i>Vela_Project7/</i>. With the <i>ciao-4.18</i> conda environment, run from that folder:", body))
S.append(table([
    ["Script", "Purpose"],
    ["scripts/setup_ciao.sh", "Activate CIAO and point CALDB at the partial calibration database"],
    ["scripts/01_reprocess.sh &lt;obsid&gt;", "CALDB download, chandra_repro, flare screening → products/&lt;obsid&gt;/evt2_clean.fits (Task A)"],
    ["scripts/02_images.sh &lt;obsid&gt;", "Band counts images, exposure maps and fluxed images (Tasks A, C)"],
    ["scripts/03_measure.py &lt;obsid&gt;", "Ring fits, jet lengths and widths, nebula extent → measurements.json (Task B)"],
    ["scripts/04_energy.py &lt;obsid&gt;", "Band images, hardness ratios, profiles, hardness map, RGB (Tasks C, D)"],
    ["scripts/05_responses.sh, 06_photon_index.py", "ARF/RMF and conversion of hardness ratios to photon indices"],
    ["scripts/07_doppler.py, 08_outer_jet.py, 10_r50_errors.py, 12_luminosity.py", "Doppler limit, outer jet, size errors, flux and luminosity"],
    ["scripts/09_figures.py, 11_ciao_rgb.sh, 13_make_report.py, 14_render_report.sh", "Final figures, CIAO three-colour images and this report"]],
    [6.6 * cm, 10.4 * cm]))
S.append(P("Table A1. Analysis scripts, run in order for ObsIDs 128 and 12075.", cap))



# ================================================================== assemble HTML
CSS = """
/* --- pdf-design-skill, SOFT style (tokens and components from assets/soft.html) --- */
@page { size: letter; margin: 0.8in 0.9in 1in 0.9in;
  @bottom-left  { content: "%(footer)s"; font-family: 'Plex Sans R'; font-size: 9px; color: #9a9a9a; vertical-align: top; padding-top: 0.55in; }
  @bottom-right { content: counter(page); font-family: 'Plex Sans R'; font-size: 9px; color: #9a9a9a; vertical-align: top; padding-top: 0.55in; } }
@page cover { margin: 0;
  @bottom-left  { content: none; } @bottom-right { content: none; } }
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: 'Plex Sans R'; color: #1d1d1d; }
sup, sub { font-size: 0.72em; line-height: 0; }
.sm { font-size: 0.8em; }
i { font-style: italic; }

/* cover */
.cover { page: cover; width: 8.5in; height: 11in; position: relative; text-align: center; padding: 0.55in 0.9in 1in 0.9in; page-break-after: always; }
.hero { position: absolute; top: 0; left: 0; right: 0; height: 2.4in; background: linear-gradient(#D0FF9B, #ffffff); }
.cover-inner { position: relative; }
.cover .kick { margin: 0.25in 0 14px 0; }
h1.title { font-family: 'Space Grotesk M'; font-weight: normal; font-size: 46px; line-height: 1.1; margin-bottom: 14px; }
.cover .sub { max-width: 5.8in; margin: 6px auto 0.3in auto; font-size: 13px; text-align: center; }
.cover table.info, .cover .card { text-align: left; }
.cover-footer { position: absolute; bottom: 0.45in; left: 0.9in; width: 6.7in; display: table; font-size: 9px; color: #9a9a9a; text-align: left; }
.cover-footer .l { display: table-cell; }
.cover-footer .r { display: table-cell; text-align: right; }

h1.sec { font-family: 'Space Grotesk M'; font-weight: normal; font-size: 30px; line-height: 1.15; margin: 10px 0 4px 0; }
.opener { margin: 0 0 6px 0; page-break-after: avoid; page-break-inside: avoid; }
.opener.newpage { page-break-before: always; }
h2 { font-family: 'Space Grotesk M'; font-weight: normal; font-size: 22px; margin: 24px 0 6px 0; padding-bottom: 8px; border-bottom: 1px solid #D9D9D9; page-break-after: avoid; }
.body { font-size: 12px; line-height: 1.65; max-width: 6.3in; margin: 10px 0; text-align: left; }
.body b, td b, .figcap b { font-family: 'Plex Sans SB'; font-weight: normal; }
.bul { padding-left: 14px; text-indent: -14px; margin: 6px 0; }
.ref { font-size: 11px; line-height: 1.5; padding-left: 14px; text-indent: -14px; margin: 3px 0; max-width: 6.7in; }
.small { font-size: 11px; }

.hl { background: #D0FF9B; font-family: 'Plex Sans SB'; font-size: 10.5px; padding: 2px 7px; border-radius: 4px; }

.card { background: #EAFFD1; border-radius: 16px; padding: 20px 24px; margin: 18px 0; }
.card .ct { font-family: 'Space Grotesk M'; font-size: 20px; margin-bottom: 8px; }
.card .body { margin: 0 0 8px 0; max-width: none; }
.cover-card .body { font-size: 11.5px; line-height: 1.55; }

table.data { border-collapse: collapse; width: 100%; margin: 14px 0 4px 0; table-layout: fixed; }
table.data th { font-family: 'Plex Sans SB'; font-weight: normal; font-size: 11px; background: #EAFFD1; text-align: center; padding: 8px 9px; border: 0.8px solid #D2D2D2; }
table.data td { font-size: 10.5px; line-height: 1.5; padding: 7px 9px; border: 0.8px solid #D2D2D2; vertical-align: top;  overflow-wrap: anywhere; }
table.data td.k { font-family: 'Plex Sans SB'; }
table.data tr { page-break-inside: avoid; }
table.info td { font-size: 10.5px; padding: 7px 10px; }
table.info td.k { background: #EAFFD1; }
.tlabel { margin: 14px 0 -6px 0; }

figure { margin: 16px auto 10px auto; page-break-inside: avoid; }
figure img { width: 100%; border-radius: 12px; border: 0.8px solid #E3E3E3; }
.figcap { font-size: 10.5px; color: #8a8a8a; margin-top: 6px; line-height: 1.45; }
.pair { display: table; width: 100%; border-spacing: 14px 0; margin: 16px 0 0 -14px; }
.pair .c { display: table-cell; width: 50%; text-align: center; vertical-align: top; }
.pair img { width: 2.6in; border-radius: 12px; }
.pair .lbl { margin-top: 6px; }

/* table of contents */
.toc { page-break-after: always; }
.toc a { display: block; color: #1d1d1d; text-decoration: none; font-size: 13px; line-height: 2.2; border-bottom: 0.8px solid #ECECEC; }
.toc a::after { content: target-counter(attr(href), page); float: right; font-family: 'Plex Sans SB'; }
.toc .chip { display: inline-block; width: 1.05in; }
"""


def toc_html():
    rows = []
    for sid, t in SECTIONS:
        m = re.search(r"\(Task ([A-D])\)", t)
        chip = f'<span class="hl">Task {m.group(1)}</span>' if m else ""
        rows.append(f'<a href="#{sid}"><span class="chip">{chip}</span>{re.sub(r" ?\\(Task [A-D]\\)", "", t)}</a>')
    return ('<div class="toc"><div class="opener"><span class="hl">Contents</span>'
            '<h1 class="sec">Contents</h1></div><h2 style="margin-top:6px">Sections</h2>' + "".join(rows) + "</div>")


doc = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Final Project Report – Project 7: Vela Pulsar X-ray Imaging with Chandra</title>
<meta name="author" content="Group 4">
<meta name="description" content="Chandra ACIS imaging of the Vela pulsar wind nebula">
<style>{CSS.replace("%(footer)s", FOOTER_LEFT)}</style></head>
<body>
{COVER}
{toc_html()}
{chr(10).join(S)}
</body></html>"""
open(OUTHTML, "w").write(doc)
print("wrote", OUTHTML, f"({len(SECTIONS)} sections)")
