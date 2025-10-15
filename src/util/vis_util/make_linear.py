from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter, ScalarFormatter
import sys
import os

# Add the parent directory of 'vis_util' to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import feds_util
import general_util as gen_util

LABELS = {
    # fire_spread
    "fline":"Active fire line",
    "fperim":"Fire perimeter",
    "nfp":"New fire pixels",

    # fuel_topo
    "adj":"Surface spread rate adjustment",
    "asp":"Topographic aspect (°)",
    "cbd":"Canopy bulk density (×100 kg m⁻³)",
    "cbh":"Canopy base height (×10 m)",
    "cc":"Canopy cover (%)",
    "ch":"Canopy height (×10 m)",
    "dem":"Digital elevation model (m)",
    "fbfm40":"Fire behavior fuel model",
    "ignition_mask":"Ignition mask",
    "phi":"Initial φ (level set variable)",
    "slp":"Topographic slope (°)",

    # high_res_climate
    "lh":"Live herbaceous fuel moisture (%)",
    "lw":"Live woody fuel moisture (%)",
    "m1":"1-h dead fuel moisture (%)",
    "m10":"10-h dead fuel moisture (%)",
    "m100":"100-h dead fuel moisture (%)",
    "wd":"20-ft wind direction (°)",
    "ws":"20-ft wind speed (mph)",

    # FRP (values are in CSV as given)
    "frp_density":"Fire radiative power density (MW m⁻²)",
    "frp_total_mw":"Total fire radiative power (MW)",
    "frp_on_fline_mw":"FRP on active fire line (MW)",

    # landfire
    "200evt":"Existing Vegetation Type",
    "200f13_20":"13 Anderson Fire Behavior Fuel Models 2020",
    "200f40_20":"40 Scott and Burgan Fire Behavior Fuel Models 2020",
    "asp2020":"Topographic Aspect (°)",
    "elev2020":"Topographic Elevation (m)",
    "slpd2020":"Topographic Slope Degrees (°)",

    # low_res_climate
    "d2m":"2-m dew-point temperature (°C)",
    "sp":"Surface pressure (Pa)",
    "t2m":"2-m temperature (°C)",
    "tp":"Total precipitation (mm)",

    # derived
    "vpd":"Vapor-pressure deficit (kPa)",
    "mean_rh":"Mean relative humidity (%)",
    "fline_mean_elev_m":"Mean elevation over fline (m)",
}

FEDS_VARS = ("fline","fperim","nfp")
SKIP_ALWAYS = {
    "UTC","Local","fline","fperim","nfp",
    "fline_area_km2","fperim_area_km2","nfp_area_km2"
}

def _slugify(name: str) -> str:
    s = (name or "").strip()
    s = re.sub(r"[^\w]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:80]

def _make_axis_converters(local_series: pd.Series, utc_series: pd.Series):
    xL = mdates.date2num(pd.to_datetime(local_series))
    xU = mdates.date2num(pd.to_datetime(utc_series))

    # Sort by Local and drop duplicate Local values
    order = np.argsort(xL)
    xL, xU = xL[order], xU[order]
    maskL = np.concatenate([[True], np.diff(xL) > 0])
    xL, xU = xL[maskL], xU[maskL]

    def local2utc(x):
        return np.interp(x, xL, xU, left=xU[0], right=xU[-1])

    orderU = np.argsort(xU)
    xU2, xL2 = xU[orderU], xL[orderU]
    maskU = np.concatenate([[True], np.diff(xU2) > 0])
    xU2, xL2 = xU2[maskU], xL2[maskU]

    def utc2local(x):
        return np.interp(x, xU2, xL2, left=xL2[0], right=xL2[-1])

    return local2utc, utc2local

def _span_days(lo, hi) -> float:
    dt = pd.to_datetime(hi) - pd.to_datetime(lo)
    return dt.total_seconds() / 86400.0

def _best_loc_and_fmt(lo, hi):
    """
    Choose a locator/formatter based on span (in days).
    Smallest interval is 1 day (no sub-daily ticks).
    """
    d = _span_days(lo, hi)
    if d <= 10:
        loc = mdates.DayLocator(interval=1)
        fmt = mdates.DateFormatter("%b %d")
    elif d <= 35:
        loc = mdates.WeekdayLocator(byweekday=mdates.MO, interval=1)
        fmt = mdates.DateFormatter("%b %d")
    elif d <= 120:
        loc = mdates.MonthLocator(interval=1)
        fmt = mdates.DateFormatter("%b %Y")
    else:
        loc = mdates.MonthLocator(interval=2)
        fmt = mdates.DateFormatter("%b %Y")
    return loc, fmt

def _grid_freq(lo, hi) -> str:
    d = _span_days(lo, hi)
    if d <= 10:   return "D"
    if d <= 35:   return "W-MON"
    if d <= 120:  return "MS"
    return "2MS"

def _iter_ticks(start, end, freq):
    s = pd.to_datetime(start)
    e = pd.to_datetime(end)
    rng = pd.date_range(start=s.floor("D"), end=e.ceil("D"), freq=freq)
    return rng[(rng >= s) & (rng <= e)]

def _maybe_wider_fig(lo, hi, base=(12, 4)):
    d = _span_days(lo, hi)
    width = max(base[0], min(22, 12 + 0.12 * d))
    return (width, base[1])

def _frp_density_display_unit(max_mw_per_m2: float):
    m = max_mw_per_m2
    if not np.isfinite(m) or m == 0.0:
        return 1e6, "W m⁻²"
    if m >= 1.0:
        return 1.0, "MW m⁻²"
    if m >= 1e-3:
        return 1e3, "kW m⁻²"
    if m >= 1e-6:
        return 1e6, "W m⁻²"
    if m >= 1e-9:
        return 1e9, "mW m⁻²"
    if m >= 1e-12:
        return 1e12, "µW m⁻²"
    return 1e15, "nW m⁻²"

def make_plots(fid: str) -> None:
    # Fire metadata
    firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)
    mask = firelist["Event_ID"].astype(str) == str(fid)
    if not mask.any():
        raise KeyError(f"Event_ID {fid} not found in firelist: {feds_util.feds_firelist}")
    row = firelist.loc[mask].iloc[0]

    fire_name = str(row.get("Incid_Name", "")).strip()
    fire_slug = _slugify(fire_name)

    # Output dirs
    dir_suffix = f"{fid}_{fire_slug}" if fire_slug else fid
    base_dir   = Path(gen_util.dir_output) / "sum_vis" / dir_suffix
    base_dir.mkdir(parents=True, exist_ok=True)
    linear_dir = base_dir / "linear_plots"
    linear_dir.mkdir(parents=True, exist_ok=True)

    csv_path = base_dir / f"data_table_{dir_suffix}.csv"
    if not csv_path.exists():
        matches = list(base_dir.glob(f"data_table_{fid}_*.csv"))
        if not matches:
            raise FileNotFoundError(f"No CSV found in {base_dir} for {fid}")
        csv_path = matches[0]

    # Load & sort
    df = (
        pd.read_csv(csv_path, sep="\t", parse_dates=["UTC","Local"])
          .sort_values("Local")
          .reset_index(drop=True)
    )

    # full window from CSV
    LOCAL_START = df["Local"].iloc[0]
    LOCAL_END   = df["Local"].iloc[-1]
    UTC_START   = df["UTC"].iloc[0]
    UTC_END     = df["UTC"].iloc[-1]

    # mapping functions for a true UTC top axis + gridline positions
    local2utc, utc2local = _make_axis_converters(df["Local"], df["UTC"])

    def _title_prefix() -> str:
        return f"{fid} – {fire_name}" if fire_name else f"{fid}"

    def _add_synced_utc_axis(ax, loU: pd.Timestamp, hiU: pd.Timestamp):
        secax = ax.secondary_xaxis("top", functions=(local2utc, utc2local))
        locU, fmtU = _best_loc_and_fmt(loU, hiU)
        secax.xaxis.set_major_locator(locU)
        secax.xaxis.set_major_formatter(fmtU)
        secax.set_xlabel(f"UTC time\n{loU:%b %d %H:%M}  →  {hiU:%b %d %H:%M}")
        return secax

    def _draw_utc_gridlines_on_local(ax):
        gfreq = _grid_freq(UTC_START, UTC_END)
        for ut in _iter_ticks(UTC_START, UTC_END, gfreq):
            num = utc2local(mdates.date2num(pd.to_datetime(ut)))
            ax.axvline(mdates.num2date(num),
                       linestyle=":", linewidth=0.7, alpha=0.45, color="grey")

    def _draw_local_guides(ax):
        gfreq = _grid_freq(LOCAL_START, LOCAL_END)
        for dt in _iter_ticks(LOCAL_START, LOCAL_END, gfreq):
            ax.axvline(dt, linestyle=":", linewidth=0.8, alpha=0.6, color="grey")

    def plot_feds(var: str) -> None:
        area_col = f"{var}_area_km2"
        cols = ["UTC","Local",var] + ([area_col] if area_col in df.columns else [])
        d = df[cols].dropna(subset=[var])
        if d.empty:
            print(f"[skip] {var}: only NaNs")
            return

        # captions from the data actually plotted
        loL, hiL = d["Local"].iloc[[0, -1]]
        loU, hiU = d["UTC"].iloc[[0, -1]]

        fig, ax = plt.subplots(figsize=_maybe_wider_fig(LOCAL_START, LOCAL_END))
        ax.set_axisbelow(True)

        ax.step(d["Local"], d[var], where="post", label="pixels")
        changed = d[var].ne(d[var].shift())
        ax.scatter(d.loc[changed,"Local"], d.loc[changed,var], zorder=10, s=30, label="Value changed")

        if area_col in d and d[area_col].notna().any():
            axR = ax.twinx()
            axR.step(d["Local"], d[area_col], where="post", linestyle="--")
            axR.set_ylabel("Area (km²)")

        ax.set_ylabel("Number of pixels")
        ax.set_title(f"{_title_prefix()} – {var}")
        ax.set_xlim(LOCAL_START, LOCAL_END)

        locL, fmtL = _best_loc_and_fmt(LOCAL_START, LOCAL_END)
        ax.xaxis.set_major_locator(locL)
        ax.xaxis.set_major_formatter(fmtL)
        ax.set_xlabel(f"Local time\n{loL:%b %d %H:%M}  →  {hiL:%b %d %H:%M}")

        _add_synced_utc_axis(ax, loU, hiU)
        _draw_local_guides(ax)
        _draw_utc_gridlines_on_local(ax)

        fig.tight_layout()
        fig.subplots_adjust(bottom=0.18)
        out_png = linear_dir / f"{var}_{dir_suffix}_linear_plot.png"
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        print("Saved:", out_png)

    for v in FEDS_VARS:
        if v in df.columns:
            plot_feds(v)

    def plot_column(col: str) -> None:
        area_col = f"{col}_area_km2"
        keep = ["UTC","Local",col] + ([area_col] if area_col in df.columns else [])
        d = df[keep].dropna(subset=[col])
        if d.empty:
            print(f"[skip] {col}: no data in window")
            return

        fig, ax = plt.subplots(figsize=_maybe_wider_fig(LOCAL_START, LOCAL_END))
        ax.set_axisbelow(True)

        y = d[col]

        if col == "frp_density":
            max_mw = float(np.nanmax(pd.to_numeric(y, errors="coerce"))) if y.notna().any() else 0.0
            disp_factor, disp_unit = _frp_density_display_unit(max_mw)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v * disp_factor:.3g}"))
            ylabel = f"Fire radiative power density ({disp_unit})"
        elif col in ("frp_total_mw", "frp_on_fline_mw"):
            ylabel = LABELS.get(col, col)
        else:
            sf = ScalarFormatter(useMathText=True)
            sf.set_powerlimits((-3, 3))
            ax.yaxis.set_major_formatter(sf)
            ylabel = LABELS.get(col, col)

        is_constant = y.nunique(dropna=True) == 1

        if is_constant:
            v = float(y.dropna().iloc[0])

            if col == "frp_density":
                max_mw = float(np.nanmax(pd.to_numeric(y, errors="coerce"))) if y.notna().any() else 0.0
                disp_factor, _ = _frp_density_display_unit(max_mw)
                min_span_display = 0.1                 # 0.1 in chosen unit (e.g., 0.1 W m^-2)
                span_raw = min_span_display / disp_factor  # convert back to MW m^-2
                top = max(abs(v), span_raw) * 1.2
                ax.hlines(v, LOCAL_START, LOCAL_END, linewidth=1.8)
                ax.set_ylim(-0.02 * top, top)
                loL, hiL = LOCAL_START, LOCAL_END
                loU, hiU = UTC_START,   UTC_END
            else:
                # default behavior for other constant series
                ax.hlines(v, LOCAL_START, LOCAL_END, linewidth=1.8)
                loL, hiL = LOCAL_START, LOCAL_END
                loU, hiU = UTC_START,   UTC_END
                pad = max(0.08 * (abs(v) if v != 0 else 1.0), 0.05)
                ymin, ymax = v - pad, v + pad
                if ymax - ymin < 1e-9:
                    ymin, ymax = v - 1e-6, v + 1e-6
                ax.set_ylim(ymin, ymax)

        else:
            ax.step(d["Local"], y, where="post", linewidth=1.8)
            changed = y.ne(y.shift())
            ax.scatter(d.loc[changed, "Local"], y.loc[changed],
                    zorder=10, s=28, linewidths=0.0)

            loL, hiL = d["Local"].iloc[[0, -1]]
            loU, hiU = d["UTC"].iloc[[0, -1]]

            if col == "frp_density":
                y_min_raw = float(np.nanmin(y))
                y_max_raw = float(np.nanmax(y))
                if not np.isfinite(y_max_raw) or y_max_raw <= 0:
                    y_max_raw = 1e-9
                pad_raw = 0.15 * y_max_raw
                ax.set_ylim(-0.02 * (y_max_raw + pad_raw), y_max_raw + pad_raw)
            else:
                # default behavior for other time series
                y_min, y_max = float(np.nanmin(y)), float(np.nanmax(y))
                span = y_max - y_min
                if not np.isfinite(span) or span == 0:
                    span = max(abs(y_max), abs(y_min), 1.0) * 0.1
                pad = max(0.08 * span, 0.05)
                ax.set_ylim(y_min - pad, y_max + pad)

        ax.set_ylabel(ylabel)

        if area_col in d and d[area_col].notna().any():
            axR = ax.twinx()
            axR.step(d["Local"], d[area_col], where="post", linestyle="--")
            axR.set_ylabel("Area (km²)")

        ax.set_title(f"{_title_prefix()} – {col}")
        ax.set_xlim(LOCAL_START, LOCAL_END)

        locL, fmtL = _best_loc_and_fmt(LOCAL_START, LOCAL_END)
        ax.xaxis.set_major_locator(locL)
        ax.xaxis.set_major_formatter(fmtL)

        ax.set_xlabel(f"Local time\n{loL:%b %d %H:%M}  →  {hiL:%b %d %H:%M}")
        _add_synced_utc_axis(ax, loU, hiU)

        _draw_local_guides(ax)
        _draw_utc_gridlines_on_local(ax)

        fig.tight_layout()
        fig.subplots_adjust(bottom=0.18, right=0.96)
        out_path = linear_dir / f"{col}_{dir_suffix}_linear_plot.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print("Saved:", out_path)

    for col in (c for c in df.columns if c not in SKIP_ALWAYS):
        plot_column(col)

if __name__ == "__main__":
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'
    fid_to_use = zogg_id

    if len(sys.argv) == 1:
        make_plots(fid_to_use)
    elif len(sys.argv) != 2:
        print("Usage: python src/make_linear.py <FIRE_EVENT_ID>")
        raise SystemExit(2)
    else:
        make_plots(sys.argv[1])
