from pathlib import Path
import warnings
import re
import numpy as np
import pandas as pd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
import sys
import os
from typing import Optional

# Optional timezone libs (script still works without them via fixed-offset fallback)
try:
    from timezonefinder import TimezoneFinder
    import pytz
    _TZF = TimezoneFinder()
except Exception:
    _TZF = None
    pytz = None
    warnings.filterwarnings("ignore")

# Allow `import feds_util` etc. from the parent of vis_util
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import feds_util
import general_util as gen_util
import processing_util as proc_util

# Heuristics for mapping variable stems to subfolders
SUBFOLDER_HINT = {
    "fline": "fire_spread", "fperim": "fire_spread", "nfp": "fire_spread",
    "frp": "frp",
    "lh": "high_res_climate", "lw": "high_res_climate",
    "m1": "high_res_climate", "m10": "high_res_climate", "m100": "high_res_climate",
    "wd": "high_res_climate", "ws": "high_res_climate",
    "d2m": "low_res_climate", "sp": "low_res_climate",
    "t2m": "low_res_climate", "tp": "low_res_climate",
    "adj": "fuel_topo", "asp": "fuel_topo", "cbd": "fuel_topo", "cbh": "fuel_topo",
    "cc": "fuel_topo", "ch": "fuel_topo", "dem": "fuel_topo", "fbfm40": "fuel_topo",
    "ignition_mask": "fuel_topo", "phi": "fuel_topo", "slp": "fuel_topo",
    "200evt": "landfire", "200f13_20": "landfire", "200f40_20": "landfire",
    "asp2020": "landfire", "elev2020": "landfire", "slpd2020": "landfire",
}

# Folders scanned to auto-discover available rasters
SCAN_SUBFOLDERS = [
    "fire_spread",
    "frp",
    "high_res_climate",
    "low_res_climate",
    "fuel_topo",
    "landfire",
]


def _slugify(name: str) -> str:
    """
    Create a filesystem-safe slug from a fire name.
    """
    if not name or not isinstance(name, str):
        return ""
    s = name.strip()
    s = re.sub(r"[^\w]+", "_", s, flags=re.UNICODE)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:80]


def _read_tif_bands(tif_path: Path):
    """
    Read bands and pixel area; keep other metadata for callers that need it.
    """
    with rasterio.open(tif_path) as ds:
        pix_km2 = abs(ds.transform.a * ds.transform.e) / 1_000_000  # area per pixel in km²
        bands = [ds.read(i) for i in range(1, ds.count + 1)]
        nodata = ds.nodata
        crs = ds.crs
        transform = ds.transform
        width, height = ds.width, ds.height
    return bands, pix_km2, nodata, crs, transform, width, height


def _sat_vp(temp_c):
    """
    Saturation vapor pressure (kPa) from °C.
    """
    return np.exp(temp_c * 17.27 / (temp_c + 237.3)) * 0.6108


def _mean_rh(temp_c, dew_c):
    """
    Relative humidity (%) from air temp and dew point (°C).
    """
    return _sat_vp(dew_c) / _sat_vp(temp_c) * 100.0


def _mean_elevation_series(fline_tif: Path, dem_tif: Path,
                           start_time_utc: pd.Timestamp, freq="h") -> pd.Series:
    """
    Mean elevation (m) over active fireline pixels for each FEDS band.
    """
    from numpy import nanmean
    with rasterio.open(fline_tif) as fsrc, rasterio.open(dem_tif) as dsrc:
        # Reproject DEM on the fly to fireline grid for exact masking/averaging
        with WarpedVRT(
            dsrc, crs=fsrc.crs, transform=fsrc.transform,
            width=fsrc.width, height=fsrc.height, resampling=Resampling.bilinear
        ) as dem_vrt:
            dem = dem_vrt.read(1).astype("float32")
            if dsrc.nodata is not None:
                dem = np.where(dem == dsrc.nodata, np.nan, dem)
            vals = []
            for i in range(1, fsrc.count + 1):
                mask = fsrc.read(i) > 0  # active fireline
                vals.append(nanmean(dem[mask]) if mask.any() else np.nan)
    idx = pd.date_range(start_time_utc, periods=len(vals), freq=freq, tz="UTC")
    return pd.Series(vals, index=idx)


def _timezone_of(lat: float, lon: float):
    """
    Best-effort local timezone (pytz if available; else fixed offset by longitude).
    """
    if _TZF is not None and pytz is not None:
        try:
            tzname = _TZF.timezone_at(lat=lat, lng=lon)
            if tzname:
                return pytz.timezone(tzname)
        except Exception:
            pass
    hours = int(round(lon / 15.0))

    class _Fixed:
        def __init__(self, h): self._h = h
        def utcoffset(self, _): return pd.Timedelta(hours=self._h)
        def tzname(self, _): return f"UTC{self._h:+d}"

    return _Fixed(hours)


def _local_series_to_utc(s: pd.Series, tz_local):
    """
    Convert naive/local times to timezone-aware UTC times.
    """
    s = pd.to_datetime(s)
    if pytz is not None and hasattr(tz_local, "localize"):
        if getattr(s.dt, "tz", None) is not None:
            return s.dt.tz_convert("UTC")
        return s.dt.tz_localize(tz_local, ambiguous="infer", nonexistent="shift_forward").dt.tz_convert("UTC")
    # Fallback: fixed offset shift
    offset = tz_local.utcoffset(None)
    return (s - offset).dt.tz_localize("UTC")


def _has_time_token(s: str) -> bool:
    """
    Heuristic: detect hh:mm or AM/PM in a human-entered string.
    """
    s = str(s).strip().lower()
    return (":" in s) or ("am" in s) or ("pm" in s)


def _extract_firelist_window(row: pd.Series):
    """
    Find start/end columns across multiple possible firelist schemas.
    """
    candidates = [
        ("tst", "ted"),
        ("start_local", "end_local"), ("Start_local", "End_local"),
        ("start", "end"), ("Start", "End"),
        ("t_start", "t_end"), ("t0", "t1"),
        ("start_time", "end_time"), ("Start Time", "End Time"),
    ]
    for a, b in candidates:
        if a in row and b in row and pd.notna(row[a]) and pd.notna(row[b]):
            s_raw, e_raw = str(row[a]), str(row[b])
            s = pd.to_datetime(s_raw, errors="coerce")
            e = pd.to_datetime(e_raw, errors="coerce")
            # If date only, pin start to 01:30 local to line up with FEDS t0+30m
            if not _has_time_token(s_raw):
                s = s + pd.Timedelta(hours=1, minutes=30)
            # If end has time token "12:00" without AM/PM, interpret as midnight
            if _has_time_token(e_raw):
                low = e_raw.strip().lower()
                if ("12:00" in low) and ("am" not in low) and ("pm" not in low):
                    e = e.replace(hour=0, minute=0, second=0)
            else:
                e = e.replace(hour=0, minute=0, second=0)
            return s, e, (a, b), s_raw, e_raw
    return None, None, (None, None), None, None


def _frp_series(cubes_dir: Path,
                fire_start_utc: pd.Timestamp) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """
    FRP time series (density, total, restricted to fireline).
    """
    frp_tif = cubes_dir / "frp" / "frp.tif"
    if not frp_tif.exists():
        return None, None, None, None

    fline_tif = cubes_dir / "fire_spread" / "fline.tif"
    with rasterio.open(frp_tif) as frp_ds:
        nb = frp_ds.count
        pix_m2 = abs(frp_ds.transform.a * frp_ds.transform.e)  # for total FRP in MW

        # Optional VRT to the FEDS grid for intersecting FRP and fireline per band
        fline_vrt = None
        if fline_tif.exists():
            fsrc = rasterio.open(fline_tif)
            fline_vrt = WarpedVRT(
                fsrc,
                crs=frp_ds.crs,
                transform=frp_ds.transform,
                width=frp_ds.width,
                height=frp_ds.height,
                resampling=Resampling.nearest,
            )

        dens_vals, total_vals, on_fline_vals, dens_on_fline_vals = [], [], [], []
        for i in range(1, nb + 1):
            frp_band = frp_ds.read(i).astype("float32")
            if frp_ds.nodata is not None:
                frp_band = np.where(frp_band == frp_ds.nodata, np.nan, frp_band)

            # Mean FRP density over positive pixels
            dens_vals.append(np.nanmean(np.where(frp_band > 0, frp_band, np.nan)))

            # Total FRP over positive pixels (density * area)
            vals = np.where(frp_band > 0, frp_band, np.nan)
            total = np.nansum(vals) * pix_m2 if np.any(~np.isnan(vals)) else np.nan
            total_vals.append(total)

            # Restrict to fireline pixels if FEDS exists
            on_fline = np.nan
            dens_on_fline = np.nan
            if fline_vrt is not None and i <= fline_vrt.count:
                fl_band = fline_vrt.read(i)
                mask = (fl_band > 0) & (frp_band > 0)
                if np.any(mask):
                    on_fline = np.nansum(np.where(mask, frp_band, np.nan)) * pix_m2
                    dens_on_fline = np.nanmean(np.where(mask, frp_band, np.nan))
            on_fline_vals.append(on_fline)
            dens_on_fline_vals.append(dens_on_fline)

        # FRP bands are hourly aligned; anchor to FEDS start + 30 min
        start = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
        idx = pd.date_range(start, periods=nb, freq="h", tz="UTC")
        if fline_vrt is not None:
            fline_vrt.close()
            fsrc.close()

        s_dens    = pd.Series(dens_vals, index=idx, name="frp_density")
        s_tot     = pd.Series(total_vals, index=idx, name="frp_total")
        s_fl      = pd.Series(on_fline_vals, index=idx, name="frp_on_fline")
        s_dens_fl = pd.Series(dens_on_fline_vals, index=idx, name="frp_density_on_fline")
        return s_dens, s_tot, s_fl, s_dens_fl


def _circular_mean_deg(deg_array, weights=None):
    """
    Circular mean on unit circle; robust to wrap-around (359° and 1° are close).
    """
    arr = np.asarray(deg_array, dtype="float64")
    a = np.deg2rad(arr)
    if weights is None:
        x = np.nanmean(np.cos(a))
        y = np.nanmean(np.sin(a))
    else:
        w = np.asarray(weights, dtype="float64")
        mask = ~np.isnan(arr) & ~np.isnan(w)
        if not np.any(mask):
            return np.nan
        x = np.nansum(w[mask] * np.cos(a[mask])) / np.nansum(w[mask])
        y = np.nansum(w[mask] * np.sin(a[mask])) / np.nansum(w[mask])
    ang = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0
    return ang


def _wind_series(cubes_dir: Path, grid_start_utc: pd.Timestamp) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    """
    Vector-averaged wind: compute u,v per band then aggregate to ws_mean, wd_mean.
    u_mean/v_mean stay available for vector plots or diagnostics.

    Notes:
      - wd follows meteorological FROM-direction (0° = from North, 90° = from East).
      - u positive eastward (m/s); v positive northward (m/s).
    """
    ws_tif = cubes_dir / "high_res_climate" / "ws.tif"
    wd_tif = cubes_dir / "high_res_climate" / "wd.tif"
    if not (ws_tif.exists() and wd_tif.exists()):
        return None, None, None, None

    with rasterio.open(ws_tif) as ws_ds, rasterio.open(wd_tif) as wd_ds:
        n = min(ws_ds.count, wd_ds.count)
        u_vals, v_vals, ws_vals, wd_vals = [], [], [], []
        for i in range(1, n + 1):
            ws_band = ws_ds.read(i).astype("float32")
            wd_band = wd_ds.read(i).astype("float32")
            if ws_ds.nodata is not None:
                ws_band = np.where(ws_band == ws_ds.nodata, np.nan, ws_band)
            if wd_ds.nodata is not None:
                wd_band = np.where(wd_band == wd_ds.nodata, np.nan, wd_band)

            # Convert FROM-direction and scalar speed to components
            u = -ws_band * np.sin(np.deg2rad(wd_band))  # eastward
            v = -ws_band * np.cos(np.deg2rad(wd_band))  # northward

            # Spatial mean components (vector average)
            u_mean = np.nanmean(u)
            v_mean = np.nanmean(v)

            # Back to scalar speed and FROM-direction
            ws_mean = np.hypot(u_mean, v_mean)
            wd_mean = (270.0 - np.degrees(np.arctan2(v_mean, u_mean))) % 360.0

            u_vals.append(u_mean)
            v_vals.append(v_mean)
            ws_vals.append(ws_mean)
            wd_vals.append(wd_mean)

    idx = pd.date_range(grid_start_utc, periods=len(ws_vals), freq="h", tz="UTC")
    return (
        pd.Series(u_vals,  index=idx, name="u_mean"),
        pd.Series(v_vals,  index=idx, name="v_mean"),
        pd.Series(ws_vals, index=idx, name="ws_mean"),
        pd.Series(wd_vals, index=idx, name="wd_mean"),
    )


def _aspect_circ_mean_from(tif_path: Path) -> Optional[float]:
    """
    Circular-mean aspect (deg) from an aspect GeoTIFF (single- or multi-band).
    """
    if not tif_path.exists():
        return None
    with rasterio.open(tif_path) as ds:
        band_means = []
        for i in range(1, ds.count + 1):
            a = ds.read(i).astype("float32")
            if ds.nodata is not None:
                a = np.where(a == ds.nodata, np.nan, a)
            band_means.append(_circular_mean_deg(a))
        if len(band_means) == 0:
            return None
        return float(np.nanmean(band_means))


def build_table(fid: str) -> Path:
    """
    Build a time-aligned table for a fire:
      - FEDS footprints & areas,
      - FRP series (overall and restricted to fireline),
      - Vector-averaged wind (u_mean, v_mean, ws_mean, wd_mean),
      - Circular-mean aspects per source (asp, asp2020),
      - Common derived vars (VPD, mean RH), and
      - Elevation over active fireline.

    Returns path to the saved TSV.
    """
    firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)
    row = firelist[firelist["Event_ID"] == fid].iloc[0]
    fire_name_raw = str(row.get("Incid_Name", "")).strip()
    fire_name_slug = _slugify(fire_name_raw)
    center_lat = float((row.get("lat0") + row.get("lat1")) / 2.0)
    center_lon = float((row.get("lon0") + row.get("lon1")) / 2.0)
    tz_local = _timezone_of(center_lat, center_lon)

    s_local, e_local, used_cols, s_raw, e_raw = _extract_firelist_window(row)
    if s_local is None:
        gdf_fperim_rd, _, _ = feds_util.read_1fire(fid)
        if gdf_fperim_rd is None:
            raise FileNotFoundError(f"FEDS file for fire {fid} does not exist and firelist has no times.")
        times_local = pd.to_datetime(gdf_fperim_rd.t)
        local_buf = proc_util.add_time_buffers(times_local)
        buf_min_local, buf_max_local = local_buf.min(), local_buf.max()
        fire_start_utc, fire_end_utc = _local_series_to_utc(
            pd.Series([buf_min_local, buf_max_local]), tz_local
        )
    else:
        local_buf = proc_util.add_time_buffers(pd.Series([s_local, e_local]))
        buf_min_local, buf_max_local = local_buf.min(), local_buf.max()
        fire_start_utc, fire_end_utc = _local_series_to_utc(
            pd.Series([buf_min_local, buf_max_local]), tz_local
        )

    cubes_dir = Path(gen_util.dir_output) / gen_util.dir_cubes / fid

    # Ensure timeline covers full FEDS fireline stack if present
    fline_tif = cubes_dir / "fire_spread" / "fline.tif"
    if fline_tif.exists():
        with rasterio.open(fline_tif) as ds:
            n_fline = ds.count
        fline_start = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
        last_needed = fline_start + pd.Timedelta(hours=n_fline - 1)
        if last_needed > fire_end_utc:
            fire_end_utc = last_needed

    grid_start_utc = fire_start_utc.floor("D")
    grid_end_utc   = fire_end_utc.ceil("30min")

    dir_suffix = f"{fid}_{fire_name_slug}" if fire_name_slug else fid
    out_dir = Path(gen_util.dir_output) / "sum_vis" / dir_suffix
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"data_table_{dir_suffix}.csv"

    full_utc = pd.date_range(grid_start_utc, grid_end_utc, freq="30min", tz="UTC")
    if pytz is not None and hasattr(tz_local, "localize"):
        full_local = [pd.Timestamp(t).tz_convert(tz_local).tz_localize(None) for t in full_utc]
    else:
        full_local = [(pd.Timestamp(t) + tz_local.utcoffset(None)).tz_localize(None) for t in full_utc]
    df = pd.DataFrame({"UTC": full_utc, "Local": full_local}, index=full_utc)

    # Discover available rasters
    discovered = {}
    for sub in SCAN_SUBFOLDERS:
        subdir = cubes_dir / sub
        if not subdir.exists():
            continue
        for tif in sorted(subdir.glob("*.tif")):
            var = tif.stem
            folder = SUBFOLDER_HINT.get(var, sub)
            discovered[var] = folder

    # Ensure canonical keys are considered even if not discovered by glob order
    for k in ["fline", "fperim", "nfp", "frp"]:
        folder = SUBFOLDER_HINT.get(k)
        if folder and (cubes_dir / folder / f"{k}.tif").exists():
            discovered[k] = folder

    def _series_from_stack(var: str, folder: str):
        """
        Convert a stack to a time series (means), or counts/areas for FEDS masks.
        """
        tif = cubes_dir / folder / f"{var}.tif"
        if not tif.exists():
            return None, None
        bands, pix_km2, nodata, _, _, _, _ = _read_tif_bands(tif)
        n = len(bands)

        # FEDS stacks are hourly, aligned to start+30min
        if var in {"fline", "fperim", "nfp"}:
            start = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
            idx = pd.date_range(start, periods=n, freq="h", tz="UTC")
            pix = []
            for b in bands:
                if nodata is not None:
                    b = np.where(b == nodata, 0, b)
                pix.append((b > 0).sum())
            area = [p * pix_km2 for p in pix]
            return pd.Series(pix, index=idx), pd.Series(area, index=idx)

        # Generic numeric mean per band
        vals = []
        for b in bands:
            bb = b.astype("float32")
            if nodata is not None:
                bb = np.where(bb == nodata, np.nan, bb)
            vals.append(np.nanmean(bb))

        # Unit conversions for specific ERA vars
        if var in {"t2m", "d2m"}:
            vals = [v - 273.15 if v is not None and not np.isnan(v) else v for v in vals]  # K→°C
        elif var == "tp":
            vals = [v * 100.0 * 10.0 if v is not None and not np.isnan(v) else v for v in vals]  # m→mm

        # One-band static rasters → broadcast a single value across the timeline
        is_static = (folder in {"fuel_topo", "landfire"}) or (n == 1)
        if is_static:
            value = float(vals[0]) if len(vals) else np.nan
            s_single = pd.Series([value], index=[df.index[0]])
            return s_single.reindex(df.index), None

        start = grid_start_utc
        idx = pd.date_range(start, periods=n, freq="h", tz="UTC")
        return pd.Series(vals, index=idx), None

    # FEDS counts/areas
    for var in ("fline", "fperim", "nfp"):
        if discovered.get(var) == "fire_spread":
            s, s_area = _series_from_stack(var, "fire_spread")
            if s is not None:
                df[var] = s.reindex(df.index)
            if s_area is not None:
                df[f"{var}_area_km2"] = s_area.reindex(df.index)

    # FRP series
    s_frp_dens, s_frp_tot, s_frp_fline, s_frp_dens_fline = _frp_series(cubes_dir, fire_start_utc)
    if s_frp_dens is not None:
        df["frp_density"] = s_frp_dens.reindex(df.index)                # MW/m²
    if s_frp_tot is not None:
        df["frp_total"] = s_frp_tot.reindex(df.index)                   # MW
    if s_frp_fline is not None:
        df["frp_on_fline"] = s_frp_fline.reindex(df.index)              # MW
    if s_frp_dens_fline is not None:
        df["frp_density_on_fline"] = s_frp_dens_fline.reindex(df.index) # MW/m²

    # Vector-averaged wind (also keeps components for plots/debug)
    u_s, v_s, ws_s, wd_s = _wind_series(cubes_dir, grid_start_utc)
    if u_s is not None:
        df["u_mean"]  = u_s.reindex(df.index)       # m/s, +east
        df["v_mean"]  = v_s.reindex(df.index)       # m/s, +north
        df["ws_mean"] = ws_s.reindex(df.index)      # m/s
        df["wd_mean"] = wd_s.reindex(df.index)      # deg (FROM)

    # Circular-mean aspects per source
    asp_fp     = cubes_dir / "fuel_topo" / "asp.tif"
    asp2020_fp = cubes_dir / "landfire"  / "asp2020.tif"

    asp_deg = _aspect_circ_mean_from(asp_fp)
    if asp_deg is not None:
        df["asp"] = pd.Series([asp_deg], index=[df.index[0]]).reindex(df.index)

    asp2020_deg = _aspect_circ_mean_from(asp2020_fp)
    if asp2020_deg is not None:
        df["asp2020"] = pd.Series([asp2020_deg], index=[df.index[0]]).reindex(df.index)

    # Generic ingestion for everything else
    # Skip those handled specially above to avoid duplicates/overwrite.
    skip_vars = {"fline", "fperim", "nfp", "frp", "ws", "wd", "asp", "asp2020"}
    for var, folder in discovered.items():
        if var in skip_vars:
            continue
        s, s_area = _series_from_stack(var, folder)
        if s is not None:
            df[var] = s.reindex(df.index)
        if s_area is not None:
            df[f"{var}_area_km2"] = s_area.reindex(df.index)

    # Derived variables from climate
    if {"t2m", "d2m"}.issubset(df.columns):
        df["vpd"] = _sat_vp(df["t2m"]) - _sat_vp(df["d2m"])
        df["mean_rh"] = _mean_rh(df["t2m"], df["d2m"])

    # Elevation along fireline over time
    dem_tif = cubes_dir / "fuel_topo" / "dem.tif"
    if fline_tif.exists() and dem_tif.exists():
        fline_start_for_elev = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
        df["fline_mean_elev_m"] = _mean_elevation_series(
            fline_tif, dem_tif, fline_start_for_elev
        ).reindex(df.index)

    # Column ordering (cosmetic)
    static_vars, dynamic_vars = [], []
    for var, folder in discovered.items():
        if var in {"fline", "fperim", "nfp"} or var == "frp":
            continue
        tif = cubes_dir / folder / f"{var}.tif"
        if tif.exists():
            try:
                with rasterio.open(tif) as ds:
                    is_static = (folder in {"fuel_topo", "landfire"}) or (ds.count == 1)
            except Exception:
                is_static = (folder in {"fuel_topo", "landfire"})
            (static_vars if is_static else dynamic_vars).append(var)

    feds_area_cols = [c for c in (f"{v}_area_km2" for v in ("fline","fperim","nfp")) if c in df.columns]
    feds_count_cols = [c for c in ("fline","fperim","nfp") if c in df.columns]
    frp_cols = [c for c in ("frp_density","frp_total","frp_on_fline","frp_density_on_fline") if c in df.columns]
    wind_cols = [c for c in ("u_mean","v_mean","ws_mean","wd_mean") if c in df.columns]

    ordered = (
        ["UTC","Local"]
        + feds_area_cols
        + feds_count_cols
        + frp_cols
        + wind_cols
        + sorted([v for v in dynamic_vars if v in df.columns])
        + sorted([v for v in static_vars if v in df.columns])
    )

    remaining = [c for c in df.columns if c not in set(ordered)]
    ordered += sorted(remaining)

    df = df[[c for c in ordered if c in df.columns]]

    # Save with UTC naive timestamps in the file for easy CSV tooling
    df_out = df.copy()
    df_out["UTC"] = df_out["UTC"].dt.tz_convert("UTC").dt.tz_localize(None)
    out_dir.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_csv, sep="\t", index=False, float_format="%.6f", date_format="%Y-%m-%d %H:%M:%S")
    print(f"Saved ({fire_name_raw or 'Unknown fire'}):", out_csv)
    return out_csv


if __name__ == "__main__":
    """
    Run as:
        python -m src.util.vis_util.make_table <FIRE_EVENT_ID>
        # or (defaults to zogg_id if no arg)
        python src/util/vis_util/make_table.py
    """
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'
    fid_to_use = zogg_id

    if len(sys.argv) == 1:
        build_table(fid_to_use)
    elif len(sys.argv) != 2:
        print("Usage: python src/make_table.py <FIRE_EVENT_ID>")
        raise SystemExit(2)
    else:
        build_table(sys.argv[1])
