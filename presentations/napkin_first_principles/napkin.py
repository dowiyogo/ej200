#!/usr/bin/env python3
"""
napkin.py  —  First-principles photon budget for the SHiP timing bar napkin deck.

Model:  All photons are eventually guided to one of the two ends by Vikuiti reflectors
        (R_nonTIR=0.95) + TIR (R=1). Losses: (a) bulk attenuation, (b) Vikuiti 5% /bounce.
        We report f_eff = Npe_sim / N_produced to capture the combined collection×PDE factor.
        Material-to-material differences are dominated by yield × bulk_survival(L/2).

Outputs (in values/):
  napkin_values.csv   — one row per quantity: name, value, unit, description
  napkin_macros.tex   — LaTeX \\newcommand macros, never hardcoded in .tex

Physics:
  Knoll §8.3 Eq. 8.15  —  f_TIR = 2·cos(θ_c) − 1  for slab geometry
  θ_c = arcsin(1/n)    with n = 1.58 (PVT)
  v_group = c / n
  σ_timing napkin      ≈  τ_decay / √(N_pe)
  σ_x                  =  σ_END · v_group / √2
"""

import csv
import math
import numpy as np
from pathlib import Path

DECK   = Path(__file__).parent
VALUES = DECK / "values"
VALUES.mkdir(exist_ok=True)

# ── Geometry & universal constants ─────────────────────────────────────────
n          = 1.58          # PVT refractive index
L_mm       = 1400.0        # bar length  (mm)
W_mm       = 60.0          # bar width   (mm)
H_mm       = 10.0          # bar height  (mm)
c_mm_ns    = 299.792458    # speed of light (mm/ns)
R_vik      = 0.95          # Vikuiti non-TIR reflectivity  (Materials.cc, Geant4)
R_vik_nap  = 0.98          # Vikuiti reflectivity in napkin slides (manufacturer spec)

# ── Muon beam ──────────────────────────────────────────────────────────────
track_mm     = H_mm           # 1 GeV μ⁻ vertical through H
dEdx_MeV_cm  = 2.0            # MIP in PVT  (MeV/cm)
dE_MeV       = dEdx_MeV_cm * (track_mm / 10.0)   # 0.2 MeV per mm × 10 mm = 2.0 MeV

# ── SiPM geometry at END face ──────────────────────────────────────────────
n_sipm_per_end  = 8           # SiPMs per END face
sipm_side_mm    = 6.0         # SiPM active square side (mm)
sipm_area_mm2   = sipm_side_mm**2          # 36 mm²
end_face_mm2    = W_mm * H_mm              # 600 mm²
sipm_active_mm2 = n_sipm_per_end * sipm_area_mm2   # 288 mm²
f_coverage      = sipm_active_mm2 / end_face_mm2   # 0.48
PDE_napkin      = 0.40        # assumed PDE for napkin budget

# ── Angle at critical ray (just outside TIR, β measured from bar axis) ──────
# β_c = 90° − θ_c  (angle from bar axis = complement of critical angle from surface)
# tan(β_c) = tan(90° − θ_c) = 1/tan(θ_c) ... actually β_c is the angle FROM bar axis
# The ray just outside TIR makes angle β_c = 90° - θ_c with the bar surface normal,
# equivalently it makes angle α_c = θ_c from the surface.
# For bounce counting: angle from bar axis = 90° - θ_c ≈ 50.7°
alpha_c_deg = 90.0 - math.degrees(math.asin(1.0 / n))   # ≈ 50.7° (from surface)
tan_alpha_c = math.tan(math.radians(alpha_c_deg))        # ≈ 1.22

# ── Material table  (from src/Materials.cc) ────────────────────────────────
materials = {
    "EJ-200": dict(yield_per_MeV=10000, tau_rise_ns=0.9,  tau_decay_ns=2.1, latt_cm=380),
    "EJ-204": dict(yield_per_MeV=10400, tau_rise_ns=0.7,  tau_decay_ns=1.8, latt_cm=160),
    "EJ-230": dict(yield_per_MeV=9700,  tau_rise_ns=0.5,  tau_decay_ns=1.5, latt_cm=120),
}

# ── Simulation results at x=0  (from convention.json + phase_ab_optimal.csv) ──
sim_npe   = {"EJ-200": 1310.8, "EJ-204": 941.0, "EJ-230": 701.3}   # PE/end
sim_sigma = {"EJ-200": 47.6,   "EJ-204": 52.1,  "EJ-230": 49.6}    # ps
# Effective propagation velocity and timing (EJ-230 END-only, from Geant4 analysis)
veff_end_mm_ns  = 184.6       # mm/ns  (END SiPMs, EJ-230)
sigma_Dt_end_ps = 111.5       # ps     (σ of t_R − t_L for END, EJ-230)
sigma_x_end_mm  = (veff_end_mm_ns / 2.0) * (sigma_Dt_end_ps * 1e-3)   # mm = 10.3
# BLUE best result (EJ-230, END + 20 TOP SiPMs, x=0)
sigma_END_best_ps  = 53.68
sigma_TOP_best_ps  = 15.20
sigma_BLUE_best_ps = 15.21
w_TOP_best         = 0.987
w_END_best         = 0.013

# ══════════════════════════════════════════════════════════════════════════
# 1.  Critical angle and TIR fraction  (Knoll 8.15)
# ══════════════════════════════════════════════════════════════════════════
theta_c_rad   = math.asin(1.0 / n)
theta_c_deg   = math.degrees(theta_c_rad)
f_TIR         = 2.0 * math.cos(theta_c_rad) - 1.0   # ≈ 0.548  (trapped in thin slab)
f_TIR_per_end = f_TIR / 2.0                          # ≈ 0.274  (one end by symmetry)

# ══════════════════════════════════════════════════════════════════════════
# 2.  Group velocity
# ══════════════════════════════════════════════════════════════════════════
v_group_mm_ns = c_mm_ns / n    # ≈ 189.7 mm/ns

# ══════════════════════════════════════════════════════════════════════════
# 3.  Vikuiti bounce geometry  (photons zig-zagging to one end)
#
#  A trapped photon bouncing between the H-faces makes one bounce every 2·H.
#  To travel L/2 in the x-direction it needs N_bounce = (L/2) / (2·H) bounces.
#  Each bounce has survival R_vik.  Bulk: exp(–L/2 / λ_att).
#
#  Note: non-trapped photons also reach the end via Vikuiti reflection on the
#  W-faces.  The full 3D problem is complex; the napkin quotes N_bounce only
#  for the TIR-guided fraction.  Absolute Npe is taken from simulation.
# ══════════════════════════════════════════════════════════════════════════
N_bounce          = (L_mm / 2.0) / (2.0 * H_mm)    # 35 for our bar
R_surv_bounce     = R_vik ** N_bounce                # Vikuiti-only survival

# ══════════════════════════════════════════════════════════════════════════
# 4.  Per-material derived quantities
# ══════════════════════════════════════════════════════════════════════════
path_mm = L_mm / 2.0   # center → one end

bulk_surv  = {}   # exp(–L/2 / λ_att)
N_produced = {}   # yield × dE
f_eff      = {}   # Npe_sim / N_produced  (collection × PDE combined)
f_eff_norm = {}   # f_eff normalised to EJ-200  (material comparison)
sigma_x    = {}   # position resolution (mm)

for mat, p in materials.items():
    latt_mm  = p["latt_cm"] * 10.0
    bulk_surv[mat]  = math.exp(-path_mm / latt_mm)
    N_produced[mat] = p["yield_per_MeV"] * dE_MeV
    f_eff[mat]      = sim_npe[mat] / N_produced[mat]
    sigma_x[mat]    = sim_sigma[mat] * 1e-3 * v_group_mm_ns / math.sqrt(2.0)

ej200_f = f_eff["EJ-200"]
for mat in materials:
    f_eff_norm[mat] = f_eff[mat] / ej200_f

# Napkin cross-check: ratio EJ-204/EJ-200 and EJ-230/EJ-200 from first principles
# Predicted ratio: (yield_B / yield_A) × (bulk_surv_B / bulk_surv_A)
npe_ratio_pred = {}
npe_ratio_obs  = {}
for mat in ["EJ-204", "EJ-230"]:
    npe_ratio_obs[mat]  = sim_npe[mat] / sim_npe["EJ-200"]
    npe_ratio_pred[mat] = ((materials[mat]["yield_per_MeV"] / materials["EJ-200"]["yield_per_MeV"])
                           * (bulk_surv[mat] / bulk_surv["EJ-200"]))

# ══════════════════════════════════════════════════════════════════════════
# 5.  Build value table
# ══════════════════════════════════════════════════════════════════════════
rows = []

def add(name, value, unit, desc):
    rows.append({"name": name, "value": float(value), "unit": unit, "description": desc})

# Universal
add("n_pvt",          n,              "",        "PVT refractive index")
add("theta_c_deg",    theta_c_deg,    "deg",     "Critical angle: arcsin(1/n)")
add("f_TIR",          f_TIR,          "",        "TIR fraction (slab): 2cos(θ_c)−1  [Knoll 8.15]")
add("f_TIR_per_end",  f_TIR_per_end,  "",        "TIR fraction toward one end: f_TIR/2")
add("v_group_mm_ns",  v_group_mm_ns,  "mm/ns",   "Group velocity: c/n")
add("dE_MeV",         dE_MeV,         "MeV",     "Energy deposited: 1-GeV μ⁻ through H=10 mm")
add("R_vikuiti",      R_vik,          "",        "Vikuiti non-TIR reflectivity (Materials.cc)")
add("N_bounce",       N_bounce,       "",        "Bounces to traverse L/2 along H-faces")
add("R_surv_bounce",  R_surv_bounce,  "",        "Vikuiti bounce survival: R^N_bounce")
add("dEdx_MeV_cm",    dEdx_MeV_cm,    "MeV/cm",  "MIP dE/dx in PVT (minimum ionising)")
add("L_mm",           L_mm,           "mm",      "Bar length")
add("W_mm",           W_mm,           "mm",      "Bar width")
add("H_mm",           H_mm,           "mm",      "Bar height (thin direction)")
add("n_sipm_per_end", n_sipm_per_end, "",        "SiPM per END face")

# Per-material
for mat, p in materials.items():
    tag = mat.lower().replace("-", "")
    add(f"yield_{tag}",      p["yield_per_MeV"],  "ph/MeV",  f"{mat} scintillation yield")
    add(f"tau_rise_{tag}",   p["tau_rise_ns"],     "ns",      f"{mat} rise time")
    add(f"tau_decay_{tag}",  p["tau_decay_ns"],    "ns",      f"{mat} decay time")
    add(f"latt_mm_{tag}",    p["latt_cm"] * 10.0, "mm",      f"{mat} attenuation length (mm)")
    add(f"latt_cm_{tag}",    p["latt_cm"],         "cm",      f"{mat} attenuation length (cm)")
    add(f"bulk_surv_{tag}",  bulk_surv[mat],       "",        f"{mat} bulk survival: exp(−700/{p['latt_cm']*10:.0f})")
    add(f"N_prod_{tag}",     N_produced[mat],      "photons", f"{mat} photons produced: yield×dE")
    add(f"Npe_sim_{tag}",    sim_npe[mat],         "pe",      f"{mat} simulated Npe/end at x=0")
    add(f"f_eff_{tag}",      f_eff[mat],           "",        f"{mat} collection×PDE: Npe_sim/N_produced")
    add(f"sigma_sim_{tag}",  sim_sigma[mat],       "ps",      f"{mat} simulated σ(T0) at x=0")
    add(f"sigma_x_{tag}",    sigma_x[mat],         "mm",      f"{mat} position resolution: σ·v_g/√2")

# Napkin cross-check ratios
for mat in ["EJ-204", "EJ-230"]:
    tag = mat.lower().replace("-", "")
    add(f"ratio_pred_{tag}",  npe_ratio_pred[mat],  "",  f"{mat}/EJ-200 Npe ratio (napkin prediction)")
    add(f"ratio_obs_{tag}",   npe_ratio_obs[mat],   "",  f"{mat}/EJ-200 Npe ratio (simulation)")

# ── Slide 7 & 14: velocity and timing ──────────────────────────────────────
add("veff_end_mm_ns",   veff_end_mm_ns,       "mm/ns",  "Effective propagation velocity (END, EJ-230)")
add("sigma_Dt_end_ps",  sigma_Dt_end_ps,      "ps",     "σ of (t_R−t_L) for END estimator")
add("sigma_x_end_mm",   sigma_x_end_mm,       "mm",     "Position resolution: (veff/2)·σΔt")
beta_eff_end_deg = math.degrees(math.acos(veff_end_mm_ns / v_group_mm_ns))
add("beta_eff_end_deg", beta_eff_end_deg,     "deg",    "Effective axial angle END: arccos(veff/vg)")

# ── Slide 8: bounce counting at α_c ────────────────────────────────────────
d_half   = L_mm / 2.0      # = 700 mm
N_bH = d_half * tan_alpha_c / H_mm    # thin face (10 mm)
N_bW = d_half * tan_alpha_c / W_mm    # wide face (60 mm)
add("alpha_c_deg",    alpha_c_deg,    "deg",   "Angle from bar surface for ray just outside TIR")
add("tan_alpha_c",    tan_alpha_c,    "",      "tan(α_c) for bounce counting formula")
add("N_bounce_H",     N_bH,          "",      "Bounces on H=10 mm face to reach L/2")
add("N_bounce_W",     N_bW,          "",      "Bounces on W=60 mm face to reach L/2")

# ── Slide 9: reflector attenuation length ───────────────────────────────────
# Λ_refl = −D / (tan α_c · ln R)    [mm]
Lambda_refl_H = -H_mm / (tan_alpha_c * math.log(R_vik_nap))   # thin face, mm
Lambda_refl_W = -W_mm / (tan_alpha_c * math.log(R_vik_nap))   # wide face, mm
add("R_vik_nap",        R_vik_nap,        "",   "Vikuiti reflectivity (napkin, manufacturer spec)")
add("Lambda_refl_H_mm", Lambda_refl_H,    "mm", "Reflector attenuation length, H=10 mm face")
add("Lambda_refl_W_mm", Lambda_refl_W,    "mm", "Reflector attenuation length, W=60 mm face")

# ── Slide 12: photon budget chain (EJ-230 at x=0) ──────────────────────────
N_prod_230   = materials["EJ-230"]["yield_per_MeV"] * dE_MeV   # 19400
N_tir_end    = N_prod_230 * f_TIR_per_end                       # × 0.274
N_bulk       = N_tir_end  * bulk_surv["EJ-230"]                 # × 0.558
N_at_sipm    = N_bulk     * f_coverage                          # × 0.48
N_pe_napkin  = N_at_sipm  * PDE_napkin                          # × 0.4
N_pe_LR      = 2.0 * N_pe_napkin
add("sipm_active_mm2",  sipm_active_mm2,  "mm2",     "Total SiPM active area per END face (8×36)")
add("end_face_mm2",     end_face_mm2,     "mm2",     "END face area: W×H = 60×10")
add("f_coverage",       f_coverage,       "",        "SiPM coverage fraction: 288/600")
add("PDE_napkin",       PDE_napkin,       "",        "Assumed PDE for napkin budget")
add("N_tir_end_230",    N_tir_end,        "photons", "EJ-230: TIR-guided toward one end")
add("N_bulk_surv_230",  N_bulk,           "photons", "EJ-230: after bulk attenuation (70 cm)")
add("N_at_sipm_230",    N_at_sipm,        "photons", "EJ-230: reaching SiPM surface")
add("N_pe_napkin_230",  N_pe_napkin,      "pe",      "EJ-230: napkin Npe/end (TIR-only path)")
add("N_pe_LR_230",      N_pe_LR,          "pe",      "EJ-230: napkin Npe L+R total")

# ── Slide 13: audit table ────────────────────────────────────────────────────
# Index-matched coupling: all TIR-guided escape → 0.274 acceptance
f_idx_match  = f_TIR_per_end   # = 0.274
npe_idx_match = N_prod_230 * f_idx_match * bulk_surv["EJ-230"] * f_coverage * PDE_napkin
# Air gap: only escape cone → (1 − cos θ_c)/2
f_escape_cone = (1.0 - math.cos(theta_c_rad)) / 2.0   # ≈ 0.113
npe_air_gap   = N_prod_230 * f_escape_cone * bulk_surv["EJ-230"] * f_coverage * PDE_napkin
ratio_lo = npe_idx_match / npe_idx_match    # = 1.0  (normalised)
ratio_hi = npe_idx_match / npe_air_gap      # ratio between the two scenarios
add("f_escape_cone",   f_escape_cone,   "",   "Escape cone fraction at END: (1−cosθ_c)/2")
add("npe_idx_match",   npe_idx_match,   "pe", "EJ-230 napkin Npe/end (index-matched coupling)")
add("npe_air_gap",     npe_air_gap,     "pe", "EJ-230 napkin Npe/end (air-gap escape cone only)")
add("ratio_vik_TIR_hi", ratio_hi,       "",   "Napkin max Vikuiti/TIR ratio (air-gap case ≈ 2.4)")

# ── Slide 2: best simulation result (BLUE, EJ-230 + END + 20 TOP SiPMs) ────
add("sigma_END_best_ps",  sigma_END_best_ps,  "ps",  "Best σ_END (m*=8): EJ-230 + Vikuiti, x=0")
add("sigma_TOP_best_ps",  sigma_TOP_best_ps,  "ps",  "Best σ_TOP (k*=3): EJ-230 + 20 TOP SiPMs")
add("sigma_BLUE_best_ps", sigma_BLUE_best_ps, "ps",  "Best σ_BLUE: BLUE combination at x=0")
add("w_TOP_best",         w_TOP_best,         "",    "BLUE weight for TOP estimator")
add("w_END_best",         w_END_best,         "",    "BLUE weight for END estimator")

# ══════════════════════════════════════════════════════════════════════════
# 6.  Write CSV
# ══════════════════════════════════════════════════════════════════════════
csv_path = VALUES / "napkin_values.csv"
with open(csv_path, "w", newline="") as fout:
    writer = csv.DictWriter(fout, fieldnames=["name", "value", "unit", "description"])
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {csv_path}  ({len(rows)} rows)")

# ══════════════════════════════════════════════════════════════════════════
# 7.  Write LaTeX macros
# ══════════════════════════════════════════════════════════════════════════
def fmt(v):
    if not math.isfinite(v):
        return r"\infty"
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1000:
        return f"{v:.0f}"
    elif av >= 100:
        return f"{v:.1f}"
    elif av >= 10:
        return f"{v:.2f}"
    elif av >= 1:
        return f"{v:.3f}"
    elif av >= 0.01:
        return f"{v:.4f}"
    else:
        exp = int(math.floor(math.log10(av)))
        mant = v / 10**exp
        return rf"{mant:.2f}\times10^{{{exp}}}"

import re as _re

def _latexify_name(name):
    """Convert CSV name to a letter-only LaTeX macro suffix.
    Digits are not letters in TeX, so a macro like \\NapNpeLR230 tokenises
    as \\NapNpeLR + 230.  Map material suffixes to letter codes:
      ej200 → ejA,  ej204 → ejB,  ej230 → ejC
      bare 230/204/200 tail → C/B/A
    """
    s = name.replace("_", "").replace("-", "")
    s = _re.sub(r'ej200$', 'ejA', s)
    s = _re.sub(r'ej204$', 'ejB', s)
    s = _re.sub(r'ej230$', 'ejC', s)
    s = _re.sub(r'200$',   'A',   s)
    s = _re.sub(r'204$',   'B',   s)
    s = _re.sub(r'230$',   'C',   s)
    s = _re.sub(r'mm2$',   'mmsq', s)   # mm² suffix (not letter-safe)
    return s

tex_path = VALUES / "napkin_macros.tex"
with open(tex_path, "w") as fout:
    fout.write("% napkin_macros.tex — AUTO-GENERATED by napkin.py, do not edit\n\n")
    for row in rows:
        name_tex  = _latexify_name(row["name"])
        val_str   = fmt(row["value"])
        unit      = row["unit"]
        if unit:
            macro_val = rf"{val_str}\,\text{{{unit}}}"
        else:
            macro_val = val_str
        fout.write(f"\\newcommand{{\\Nap{name_tex}}}{{{macro_val}}}  % {row['description']}\n")
print(f"Wrote {tex_path}  ({len(rows)} macros)")

# ══════════════════════════════════════════════════════════════════════════
# 8.  Summary
# ══════════════════════════════════════════════════════════════════════════
W  = 52   # column width for table
print(f"\n{'─'*W}")
print(f"θ_c = {theta_c_deg:.2f}°   f_TIR = {f_TIR:.3f}   f_TIR/end = {f_TIR_per_end:.3f}")
print(f"v_group = {v_group_mm_ns:.2f} mm/ns   N_bounce = {N_bounce:.0f}   R^N = {R_surv_bounce:.4f}")
print(f"dE = {dE_MeV:.1f} MeV  (1-GeV μ⁻ through {H_mm:.0f} mm)\n")
hdr = f"{'Mat':<8} {'yield':>7} {'τ_d':>5} {'λ(mm)':>7} {'bulk':>6} {'N_prod':>7} {'Npe_sim':>8} {'f_eff':>7} {'σ(ps)':>7} {'σ_x(mm)':>9}"
print(hdr)
print("─"*len(hdr))
for mat in materials:
    p   = materials[mat]
    tag = mat.lower().replace("-", "")
    print(f"{mat:<8} {p['yield_per_MeV']:>7d} {p['tau_decay_ns']:>5.1f} "
          f"{p['latt_cm']*10:>7.0f} {bulk_surv[mat]:>6.3f} "
          f"{N_produced[mat]:>7.0f} {sim_npe[mat]:>8.1f} "
          f"{f_eff[mat]:>7.4f} {sim_sigma[mat]:>7.1f} {sigma_x[mat]:>9.2f}")
print()
for mat in ["EJ-204", "EJ-230"]:
    print(f"  {mat}/EJ-200 Npe: napkin={npe_ratio_pred[mat]:.3f}   sim={npe_ratio_obs[mat]:.3f}")

# ══════════════════════════════════════════════════════════════════════════
# Appended: extra constants for slides 12, 13, 7, 14
# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    pass  # already executed above
