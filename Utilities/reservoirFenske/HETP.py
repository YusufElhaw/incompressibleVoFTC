#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HETP.py — Complete postprocessing for incompressibleVoFTC rectification cases.

Data sources:
  postProcessing/<func-name>/0/values.dat          (primary source, default: reservoirFenske)
  postProcessing/apparatusBalance/0/<patch>.dat    (one file per patch, summed over group)
  postProcessing/rawPatchMolarFlux/0/<patch>.dat   (optional, one file per patch)

Output: plotsHETP/  with numbered PNG files + console summary

Usage:
    python3 HETP.py
    python3 HETP.py --case /path/to/case --out results
    python3 HETP.py --top-patches patch1 patch2 --bottom-patches patch3 patch4
    python3 HETP.py --func-name reservoirFenske --median-window 51
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Settings  — edit here, then run:  python3 HETP.py
# ---------------------------------------------------------------------------
SPECIES1_NAME = "Cyclohexan"
SPECIES2_NAME = "n-Heptan"
COL_HEIGHT_M  = 0.05       # column height [m] — blockMeshDict: y=50mm; used only for axis annotation

FUNC_NAME = "reservoirFenske"   # postProcessing subfolder written by the function object

CONDENSER_PATCHES = [           # top-apparatus patches (inletOutletCondenser BC)
    "bottom_Top_on_top",
]

BOILER_PATCHES = [              # bottom-apparatus patches (inletOutletBoiler BC)
    "bottom_Top_on_bottom",
]


# ---------------------------------------------------------------------------
# Helper functions: reading
# ---------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    """Read OpenFOAM postProcessing .dat file with '# Time ...' header."""
    header = None
    rows = []
    with path.open("r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                parts = s[1:].strip().split()
                if parts and parts[0] in ("Time", "time"):
                    header = parts
                continue
            if header is not None:
                vals = s.split()
                if len(vals) == len(header):
                    rows.append(vals)
    if header is None:
        raise ValueError(f"No '# Time ...' header found in {path}")
    df = pd.DataFrame(rows, columns=header)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["Time"]).sort_values("Time").reset_index(drop=True)


def find_pp_file(case: Path, func: str, filename: str = "values.dat") -> Path | None:
    """Find postProcessing/<func>/*/<filename>. Prefers subfolder 0."""
    base = case / "postProcessing" / func
    if not base.exists():
        return None
    p0 = base / "0" / filename
    if p0.exists():
        return p0
    cands = []
    for p in base.glob(f"*/{filename}"):
        try:
            t = float(p.parent.name)
        except ValueError:
            t = -1e300
        cands.append((t, p))
    return sorted(cands)[-1][1] if cands else None


def load_apparatus_group(case: Path, folder: str, patches: list[str]) -> pd.DataFrame | None:
    """Load and aggregate multiple patch .dat files from postProcessing/<folder>/0/.

    Flux columns (N_*, L_*, G_*, dN_*) are summed across patches.
    Composition columns (x_*) are averaged — all coupled patches share the same
    reservoir composition, so the mean equals any individual patch value.
    """
    base = case / "postProcessing" / folder / "0"
    frames = [read_table(base / f"{p}.dat")
              for p in patches
              if (base / f"{p}.dat").exists()]
    if not frames:
        return None
    if len(frames) == 1:
        return frames[0]
    aligned = [f.set_index("Time") for f in frames]
    total = aligned[0].copy()
    for f in aligned[1:]:
        total = total.add(f, fill_value=0)
    # Composition columns were summed; divide back to get the mean.
    x_cols = [c for c in total.columns if c.startswith("x")]
    if x_cols:
        total[x_cols] = total[x_cols] / len(frames)
    return total.reset_index()


def load_apparatus_pair(
    case: Path, tops: list[str], bottoms: list[str]
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Load and aggregate condenser and boiler apparatusBalance .dat files."""
    return (
        load_apparatus_group(case, "apparatusBalance", tops),
        load_apparatus_group(case, "apparatusBalance", bottoms),
    )


def load_raw_pair(
    case: Path, tops: list[str], bottoms: list[str]
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Load and aggregate condenser and boiler rawPatchMolarFlux .dat files."""
    return (
        load_apparatus_group(case, "rawPatchMolarFlux", tops),
        load_apparatus_group(case, "rawPatchMolarFlux", bottoms),
    )


# ---------------------------------------------------------------------------
# Helper functions: plotting
# ---------------------------------------------------------------------------

def med(df: pd.DataFrame, window: int) -> pd.DataFrame:
    """Apply rolling median to all numeric columns."""
    out = df.copy()
    if window <= 1:
        return out
    for c in out.columns:
        if c != "Time":
            out[c] = out[c].rolling(window, center=True, min_periods=1).median()
    return out


def savefig(path: Path, dpi: int = 180) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi)
    plt.close()


def smart_ylim(series: pd.Series, frac_last: float = 0.30,
               margin_lo: float = 0.0, margin_hi: float = 1.5) -> tuple[float, float] | None:
    """
    Compute Y-axis limits from the last frac_last fraction of the data.
    Returns None if too few valid values are available.
    """
    s = series.dropna()
    if len(s) < 5:
        return None
    n_tail = max(5, int(len(s) * frac_last))
    tail = s.iloc[-n_tail:]
    lo = float(tail.quantile(0.05))
    hi = float(tail.quantile(0.95))
    span = hi - lo if hi > lo else abs(hi) * 0.1 + 1e-30
    return (lo - margin_lo * span, hi + margin_hi * span)


def apply_smart_ylim(ax, series: pd.Series, frac_last: float = 0.30,
                     hard_min: float | None = None) -> None:
    lims = smart_ylim(series, frac_last)
    if lims is None:
        return
    lo, hi = lims
    if hard_min is not None:
        lo = max(lo, hard_min)
    if hi > lo:
        ax.set_ylim(lo, hi)


# ---------------------------------------------------------------------------
# Global conservation calculation
# ---------------------------------------------------------------------------

def compute_global_conservation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute drift of NSystem, N1System, N2System relative to t=0.
    Assumes these columns already contain column + apparatus totals.
    """
    out = df.copy()
    for col, label in [("NSystem", "N"), ("N1System", "N1"), ("N2System", "N2")]:
        if col not in df.columns:
            continue
        n0 = float(df[col].iloc[0])
        out[f"drift_{label}_mol"] = df[col] - n0
        out[f"drift_{label}_pct"] = (df[col] - n0) / max(abs(n0), 1e-300) * 100.0
    return out


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_01_konzentrationen(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    _, ax = plt.subplots(figsize=(12, 5))
    for col, label, ls in [
        ("xD", f"xD Kondensator ({SPECIES1_NAME})", "-"),
        ("xB", f"xB Boiler ({SPECIES1_NAME})", "-"),
        ("yB", f"yB Boiler-Dampf ({SPECIES1_NAME})", "--"),
    ]:
        if col in df:
            ax.plot(df["Time"], df[col], alpha=0.18, color=f"C{['xD','xB','yB'].index(col)}")
            ax.plot(m["Time"], m[col], label=label, ls=ls,
                    color=f"C{['xD','xB','yB'].index(col)}")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel(f"Molbruch {SPECIES1_NAME} [-]")
    ax.set_ylim(0, 1)
    ax.set_title("Kopf- und Sumpfkonzentration über Zeit")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "01_konzentrationen_xD_xB.png")


def plot_02_HETP(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    if "HETPColumnGeo" not in df.columns:
        return
    hetp_mm     = df["HETPColumnGeo"] * 1000.0
    hetp_mm_med = m["HETPColumnGeo"] * 1000.0

    _, ax = plt.subplots(figsize=(12, 5))
    ax.plot(df["Time"], hetp_mm, alpha=0.18, color="C0")
    ax.plot(m["Time"], hetp_mm_med, color="C0", label="HETP [mm]")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("HETP [mm]")
    ax.set_title("HETP über Zeit")
    ax.grid(True, alpha=0.4)

    # Y-axis: lower bound = 0; upper bound = max(5 * final value, column height).
    # Clips the initial transient (HETP → ∞) while showing the full convergence curve.
    s = hetp_mm_med.dropna()
    if len(s) >= 5:
        n_tail    = max(5, len(s) // 5)
        final_val = float(s.iloc[-n_tail:].median())
        y_max     = max(final_val * 5.0, COL_HEIGHT_M * 1000 * 1.1)
        ax.set_ylim(0, y_max)

    # Column height as reference line
    ax.axhline(COL_HEIGHT_M * 1000.0, color="gray", ls=":", lw=1,
               label=f"Kolonnenhöhe {COL_HEIGHT_M*1000:.0f} mm")
    ax.legend()
    savefig(out / "02_HETP.png")


def plot_03_Nmin(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    if "NminColumnGeo" not in df.columns:
        return
    _, ax = plt.subplots(figsize=(12, 5))
    ax.plot(df["Time"], df["NminColumnGeo"], alpha=0.18, color="C1")
    ax.plot(m["Time"], m["NminColumnGeo"], color="C1", label="Nmin (Fenske, inkl. Boiler)")
    ax.plot(m["Time"], m["NminColumnGeo"] - 1, color="C2", ls="--",
            label="Npacking = Nmin − 1 (Kolonnenstufen)")
    ax.axhline(1, color="gray", ls=":", lw=1, label="Nmin = 1 (nur Boiler-VLE)")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Stufenzahl [-]")
    ax.set_title("Fenske-Stufenzahl über Zeit\n"
                 "Nmin inkl. Boiler als Gleichgewichtsstufe; Npacking = Nmin − 1")
    ax.grid(True, alpha=0.4)
    # Y-axis from 0 to 1.3 * final value — shows the full rise from 0.
    s = m["NminColumnGeo"].dropna()
    if len(s) >= 5:
        n_tail    = max(5, len(s) // 5)
        final_val = float(s.iloc[-n_tail:].median())
        ax.set_ylim(0, max(final_val * 1.3, 0.5))
    ax.legend()
    savefig(out / "03_Nmin.png")


def plot_04_massenbilanz_fehler_rel(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    cols_pct  = [("drift_N_pct",  "N gesamt",             "C0", "-"),
                 ("drift_N1_pct", f"N1 ({SPECIES1_NAME})", "C1", "--"),
                 ("drift_N2_pct", f"N2 ({SPECIES2_NAME})", "C2", ":")]
    available = [(c, l, clr, ls) for c, l, clr, ls in cols_pct if c in df.columns]
    if not available:
        return
    _, ax = plt.subplots(figsize=(12, 5))
    for col, label, clr, ls in available:
        ax.plot(df["Time"], df[col], alpha=0.18, color=clr, ls=ls)
        ax.plot(m["Time"], m[col], label=label, color=clr, ls=ls)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Globale Drift [%]  (N(t)–N(0)) / N(0) · 100")
    ax.set_title("Globale Massenbilanz: relativer Fehler [%]\n"
                 "(N_Kolonne + N_Kondensator + N_Boiler)")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "04_massenbilanz_fehler_rel.png")


def plot_05_massenbilanz_fehler_abs(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    cols_mol  = [("drift_N_mol",  "ΔN gesamt",             "C0", "-"),
                 ("drift_N1_mol", f"ΔN1 ({SPECIES1_NAME})", "C1", "--"),
                 ("drift_N2_mol", f"ΔN2 ({SPECIES2_NAME})", "C2", ":")]
    available = [(c, l, clr, ls) for c, l, clr, ls in cols_mol if c in df.columns]
    if not available:
        return
    _, ax = plt.subplots(figsize=(12, 5))
    for col, label, clr, ls in available:
        ax.plot(df["Time"], df[col] * 1e6, alpha=0.18, color=clr, ls=ls)
        ax.plot(m["Time"], m[col] * 1e6, label=label, color=clr, ls=ls)
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Absolute Drift [µmol]")
    ax.set_title("Globale Massenbilanz: absoluter Fehler")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "05_massenbilanz_fehler_abs.png")


def plot_06_reservoir_fuellstand(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    cols      = [("ND",      "N Kondensator", "C0"),
                 ("NB",      "N Boiler",      "C1"),
                 ("NColumn", "N Kolonne",     "C2")]
    available = [(c, l, clr) for c, l, clr in cols if c in df.columns]
    if not available:
        return
    _, ax = plt.subplots(figsize=(12, 5))
    for col, label, clr in available:
        ax.plot(df["Time"], df[col] * 1000, alpha=0.18, color=clr)
        ax.plot(m["Time"],  m[col]  * 1000, label=label, color=clr)
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Stoffmenge [mmol]")
    ax.set_title("Füllstände: Kondensator, Boiler, Kolonne")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "06_reservoir_fuellstand.png")


def plot_07_A12(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    if "A12ColumnGeo" not in df.columns:
        return
    _, ax = plt.subplots(figsize=(12, 4))
    ax.plot(df["Time"], df["A12ColumnGeo"], alpha=0.18, color="C3")
    ax.plot(m["Time"],  m["A12ColumnGeo"],  color="C3", label="A12 (geometrisch)")
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Relative Flüchtigkeit A12 [-]")
    ax.set_title("VLE relative Flüchtigkeit A12")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "07_relative_fluechtigkeit.png")


def plot_08_apparate_molarstroeme(cond: pd.DataFrame | None,
                                   boil: pd.DataFrame | None,
                                   out: Path, window: int) -> None:
    if cond is None and boil is None:
        return
    _, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    for ax, df_app, title in [
        (axes[0], cond, "Kondensator"),
        (axes[1], boil, "Boiler"),
    ]:
        if df_app is None:
            ax.set_title(f"{title}: keine Daten")
            continue
        dm = med(df_app, window)
        for col, label, clr in [
            ("N_to",   "N → Apparat [mol/s]",   "C0"),
            ("N_from", "N ← Apparat [mol/s]",   "C1"),
            ("N1_to",  "N1 → Apparat [mol/s]",  "C0"),
            ("N1_from","N1 ← Apparat [mol/s]",  "C1"),
        ]:
            if col not in df_app.columns:
                continue
            ls = "-" if "N_" in col else "--"
            ax.plot(df_app["Time"], df_app[col], alpha=0.15, color=clr)
            ax.plot(dm["Time"],     dm[col],     label=label, color=clr, ls=ls)
        ax.set_ylabel("Molarstrom [mol/s]")
        ax.set_title(f"{title}: Molarströme N / N1 rein/raus")
        ax.grid(True, alpha=0.4)
        ax.legend(ncol=2, fontsize=8)
    axes[-1].set_xlabel("Zeit [s]")
    savefig(out / "08_apparate_molarstroeme.png")


def plot_09_phasenstroeme(cond: pd.DataFrame | None,
                           boil: pd.DataFrame | None,
                           out: Path, window: int) -> None:
    if cond is None and boil is None:
        return
    _, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    for ax, df_app, title in [
        (axes[0], cond, "Kondensator"),
        (axes[1], boil, "Boiler"),
    ]:
        if df_app is None:
            ax.set_title(f"{title}: keine Daten")
            continue
        dm = med(df_app, window)
        for col, label, clr in [
            ("L_to",   "Flüssig → Apparat", "C0"),
            ("L_from", "Flüssig ← Apparat", "C0"),
            ("G_to",   "Gas → Apparat",      "C1"),
            ("G_from", "Gas ← Apparat",      "C1"),
        ]:
            if col not in df_app.columns:
                continue
            ls = "-" if "_to" in col else "--"
            ax.plot(df_app["Time"], df_app[col], alpha=0.15, color=clr)
            ax.plot(dm["Time"],     dm[col],     label=label, color=clr, ls=ls)
        ax.set_ylabel("Molarstrom [mol/s]")
        ax.set_title(f"{title}: Phasenströme L / G")
        ax.grid(True, alpha=0.4)
        ax.legend(ncol=2, fontsize=8)
    axes[-1].set_xlabel("Zeit [s]")
    savefig(out / "09_phasenstroeme.png")


def plot_10_apparate_bilanzfehler(cond: pd.DataFrame | None,
                                   boil: pd.DataFrame | None,
                                   out: Path) -> None:
    if cond is None and boil is None:
        return
    _, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    for ax, df_app, title in [
        (axes[0], cond, "Kondensator"),
        (axes[1], boil, "Boiler"),
    ]:
        if df_app is None:
            ax.set_title(f"{title}: keine Daten")
            continue
        for col, label, clr in [
            ("dN_error",  "dN Fehler",  "C0"),
            ("dN1_error", "dN1 Fehler", "C1"),
            ("dN2_error", "dN2 Fehler", "C2"),
        ]:
            if col not in df_app.columns:
                continue
            ax.plot(df_app["Time"], df_app[col], label=label, color=clr, lw=0.8)
        ax.axhline(0, color="k", lw=0.5, ls="--")
        ax.set_ylabel("Fehler pro Zeitschritt [mol]")
        ax.set_title(f"{title}: Reservoir Soll/Ist-Bilanzfehler (sollte ≈ 0)")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("Zeit [s]")
    savefig(out / "10_apparate_bilanzfehler.png")


def plot_11_x_to_from(cond: pd.DataFrame | None,
                       boil: pd.DataFrame | None,
                       out: Path, window: int) -> None:
    """Composition of streams entering and leaving each apparatus."""
    if cond is None and boil is None:
        return
    _, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    for ax, df_app, title in [
        (axes[0], cond, "Kondensator"),
        (axes[1], boil, "Boiler"),
    ]:
        if df_app is None:
            ax.set_title(f"{title}: keine Daten")
            continue
        dm = med(df_app, window)
        for col, label, clr in [
            ("x_to",   f"x → Apparat ({SPECIES1_NAME})", "C0"),
            ("x_from", f"x ← Apparat ({SPECIES1_NAME})", "C1"),
        ]:
            if col not in df_app.columns:
                continue
            ax.plot(df_app["Time"], df_app[col], alpha=0.15, color=clr)
            ax.plot(dm["Time"],     dm[col],     label=label, color=clr)
        if "x_old" in df_app.columns:
            ax.plot(dm["Time"], dm["x_old"], "k--", lw=1, label="x Reservoir", alpha=0.7)
        ax.set_ylim(0, 1)
        ax.set_ylabel(f"Molbruch {SPECIES1_NAME} [-]")
        ax.set_title(f"{title}: Zusammensetzung rein/raus")
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel("Zeit [s]")
    savefig(out / "11_zusammensetzung_apparate.png")


def plot_12_raw_vs_ist(cond_raw: pd.DataFrame | None, boil_raw: pd.DataFrame | None,
                        cond_app: pd.DataFrame | None, boil_app: pd.DataFrame | None,
                        out: Path, window: int) -> None:
    """Compare raw rho/Y/W/U/Sf fluxes against apparatusBC tracked fluxes."""
    if (cond_raw is None and boil_raw is None) or (cond_app is None and boil_app is None):
        return

    _, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    for ax, raw, app, title in [
        (axes[0], cond_raw, cond_app, "Kondensator"),
        (axes[1], boil_raw, boil_app, "Boiler"),
    ]:
        if raw is None or app is None:
            ax.set_title(f"{title}: keine Daten")
            continue
        dm_raw = med(raw, window)
        dm_app = med(app, window)
        for col, label, clr in [
            ("N_to",   "N → raw",  "C0"),
            ("N_from", "N ← raw",  "C1"),
            ("N1_to",  "N1 → raw", "C2"),
            ("N1_from","N1 ← raw", "C3"),
        ]:
            if col in raw.columns:
                ax.plot(dm_raw["Time"], dm_raw[col], label=label, color=clr)
        for col, label, clr, ls in [
            ("N_to",   "N → apparatusBC",  "C0", "--"),
            ("N_from", "N ← apparatusBC",  "C1", "--"),
            ("N1_to",  "N1 → apparatusBC", "C2", "--"),
            ("N1_from","N1 ← apparatusBC", "C3", "--"),
        ]:
            if col in app.columns:
                ax.plot(dm_app["Time"], dm_app[col], label=label, color=clr, ls=ls)
        ax.set_ylabel("Molarstrom [mol/s]")
        ax.set_title(f"{title}: Raw rho/Y/W vs. apparatusBC Flüsse")
        ax.grid(True, alpha=0.4)
        ax.legend(ncol=2, fontsize=7)
    axes[-1].set_xlabel("Zeit [s]")
    savefig(out / "12_raw_vs_apparatusBC.png")


def plot_13_global_system(df: pd.DataFrame, m: pd.DataFrame, out: Path) -> None:
    """NSystem, N1System, N2System over time — the true global conservation check."""
    cols      = [("NSystem",  "N gesamt",             "C0", "-"),
                 ("N1System", f"N1 ({SPECIES1_NAME})", "C1", "--"),
                 ("N2System", f"N2 ({SPECIES2_NAME})", "C2", ":")]
    available = [(c, l, clr, ls) for c, l, clr, ls in cols if c in df.columns]
    if not available:
        return
    _, ax = plt.subplots(figsize=(12, 5))
    for col, label, clr, ls in available:
        ax.plot(df["Time"], df[col] * 1000, alpha=0.18, color=clr, ls=ls)
        ax.plot(m["Time"],  m[col]  * 1000, label=label, color=clr, ls=ls)
    ax.set_xlabel("Zeit [s]")
    ax.set_ylabel("Stoffmenge [mmol]")
    ax.set_title("Gesamtsystem N(t) = N_Kolonne + N_Kondensator + N_Boiler\n"
                 "(echter Erhaltungscheck — konstant = perfekte Erhaltung)")
    ax.grid(True, alpha=0.4)
    ax.legend()
    savefig(out / "13_gesamtsystem_N.png")


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

def print_summary(df: pd.DataFrame,
                  cond: pd.DataFrame | None, boil: pd.DataFrame | None) -> None:
    if df.empty:
        return

    last    = df.iloc[-1]
    t_end   = float(last["Time"])
    t_start = float(df.iloc[0]["Time"])

    print(f"\n{'='*60}")
    print(f"  HETP Postprocessing Summary")
    print(f"  Time range: {t_start:.4g} ... {t_end:.4g} s  ({len(df)} rows)")
    print(f"{'='*60}")

    print("\n--- Concentrations (final) ---")
    for c, label in [("xD", "xD condenser"), ("xB", "xB boiler"), ("yB", "yB boiler vapour")]:
        if c in df.columns:
            print(f"  {label:25s} = {float(last[c]):.6f}")

    print("\n--- Separation result (final) ---")
    if "xD" in df.columns and "xB" in df.columns:
        xD = float(last["xD"])
        xB = float(last["xB"])
        print(f"  dx = xD - xB              = {xD - xB:.6f}")
    if "A12ColumnGeo" in df.columns:
        print(f"  A12 (geometric)            = {float(last['A12ColumnGeo']):.4f}")
    if "NminColumnGeo" in df.columns:
        Nmin = float(last["NminColumnGeo"])
        print(f"  Nmin (Fenske)              = {Nmin:.4f}  (theoretical stages)")
    if "HETPColumnGeo" in df.columns:
        hetp = float(last["HETPColumnGeo"])
        print(f"  HETP                       = {hetp:.4f} m  =  {hetp*1000:.2f} mm")
        print(f"  Column height              = {COL_HEIGHT_M*1000:.0f} mm")
        if hetp > 0:
            print(f"  Stages (H/HETP)            = {COL_HEIGHT_M/hetp:.2f}")

    print("\n--- Global mass balance ---")
    print("  (N_column + N_condenser + N_boiler vs. t=0 — independent of flux measurement)")
    for label, col_mol, col_pct in [
        ("N  total",            "drift_N_mol",  "drift_N_pct"),
        (f"N1 ({SPECIES1_NAME[:6]})", "drift_N1_mol", "drift_N1_pct"),
        (f"N2 ({SPECIES2_NAME[:6]})", "drift_N2_mol", "drift_N2_pct"),
    ]:
        if col_mol not in df.columns:
            continue
        drift_mol = float(last[col_mol])
        drift_pct = float(last[col_pct])
        max_drift = float(df[col_mol].abs().max())
        max_pct   = float(df[col_pct].abs().max())
        print(f"  {label:12s}: drift(T)={drift_mol:+.3e} mol ({drift_pct:+.4f}%)  "
              f"max|drift|={max_drift:.3e} mol ({max_pct:.4f}%)")

    # Apparatus summary
    for name, app in [("Condenser", cond), ("Boiler", boil)]:
        if app is None or app.empty:
            continue
        last_a = app.iloc[-1]
        print(f"\n--- {name} (final) ---")
        for c, label in [("N_to","N→"), ("N_from","N←"), ("N1_to","N1→"), ("N1_from","N1←"),
                          ("L_to","L→"), ("G_to","G→"),  ("L_from","L←"), ("G_from","G←")]:
            if c in app.columns:
                print(f"  {label:6s} = {float(last_a[c]):.4e} mol/s", end="   ")
        print()
        for c, label in [("dN_error","max|dN_error|"), ("dN1_error","max|dN1_error|"),
                          ("dN2_error","max|dN2_error|")]:
            if c in app.columns:
                print(f"  {label:16s} = {float(app[c].abs().max()):.3e} mol")

    print(f"\n{'='*60}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Complete HETP postprocessing for incompressibleVoFTC rectification cases."
    )
    ap.add_argument("--case",          default=".",         help="OpenFOAM case directory.")
    ap.add_argument("--out",           default="plotsHETP", help="Output directory.")
    ap.add_argument("--func-name",     default=FUNC_NAME,
                    help="Function-object name (postProcessing subfolder). "
                         f"Default from Settings: {FUNC_NAME}")
    ap.add_argument("--top-patches",    nargs="+", default=CONDENSER_PATCHES,
                    help="One or more top-apparatus patch names (condenser side). "
                         "Multiple patches are aggregated (fluxes summed, "
                         "compositions averaged). "
                         f"Default from Settings: {CONDENSER_PATCHES}")
    ap.add_argument("--bottom-patches", nargs="+", default=BOILER_PATCHES,
                    help="One or more bottom-apparatus patch names (boiler side). "
                         f"Default from Settings: {BOILER_PATCHES}")
    ap.add_argument("--median-window", type=int, default=101,
                    help="Rolling median window width (in time steps).")
    args = ap.parse_args()

    case    = Path(args.case).resolve()
    out_dir = case / args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    window = args.median_window

    # 1) Primary data: reservoirFenske (or user-specified func-name)
    res_path = find_pp_file(case, args.func_name)
    if res_path is None:
        print(f"ERROR: postProcessing/{args.func_name} not found.")
        return 1

    df = read_table(res_path)
    df = compute_global_conservation(df)
    df.to_csv(out_dir / f"{args.func_name}.csv", index=False)
    m = med(df, window)

    # 2) Apparatus balance (condenser + boiler) — aggregate over all listed patches
    cond_app, boil_app = load_apparatus_pair(case, args.top_patches, args.bottom_patches)

    # 3) Raw patch flux (optional, for cross-check) — same patch lists
    cond_raw, boil_raw = load_raw_pair(case, args.top_patches, args.bottom_patches)

    # 4) Generate all plots
    print(f"Generating plots in {out_dir} ...")

    plot_01_konzentrationen(df, m, out_dir)
    plot_02_HETP(df, m, out_dir)
    plot_03_Nmin(df, m, out_dir)
    plot_04_massenbilanz_fehler_rel(df, m, out_dir)
    plot_05_massenbilanz_fehler_abs(df, m, out_dir)
    plot_06_reservoir_fuellstand(df, m, out_dir)
    plot_07_A12(df, m, out_dir)
    plot_08_apparate_molarstroeme(cond_app, boil_app, out_dir, window)
    plot_09_phasenstroeme(cond_app, boil_app, out_dir, window)
    plot_10_apparate_bilanzfehler(cond_app, boil_app, out_dir)
    plot_11_x_to_from(cond_app, boil_app, out_dir, window)
    plot_12_raw_vs_ist(cond_raw, boil_raw, cond_app, boil_app, out_dir, window)
    plot_13_global_system(df, m, out_dir)

    # 5) Console summary
    print_summary(df, cond_app, boil_app)

    print(f"Output: {out_dir}")
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
