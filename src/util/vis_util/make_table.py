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

try:
    from timezonefinder import TimezoneFinder
    import pytz
    _TZF = TimezoneFinder()
except Exception:
    _TZF = None
    pytz = None
    warnings.filterwarnings("ignore")

# Add the parent directory of 'vis_util' to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import feds_util
import general_util as gen_util
import processing_util as proc_util

SUBFOLDER_HINT = {
    # fire_spread
    "fline": "fire_spread", "fperim": "fire_spread", "nfp": "fire_spread",
    # FRP
    "frp": "frp",
    # high_res_climate
    "lh": "high_res_climate", "lw": "high_res_climate",
    "m1": "high_res_climate", "m10": "high_res_climate", "m100": "high_res_climate",
    "wd": "high_res_climate", "ws": "high_res_climate",
    # low_res_climate
    "d2m": "low_res_climate", "sp": "low_res_climate",
    "t2m": "low_res_climate", "tp": "low_res_climate",
    # fuel_topo
    "adj": "fuel_topo", "asp": "fuel_topo", "cbd": "fuel_topo", "cbh": "fuel_topo",
    "cc": "fuel_topo", "ch": "fuel_topo", "dem": "fuel_topo", "fbfm40": "fuel_topo",
    "ignition_mask": "fuel_topo", "phi": "fuel_topo", "slp": "fuel_topo",
    # landfire (example names; will also accept ANY name we discover there)
    "200evt": "landfire", "200f13_20": "landfire", "200f40_20": "landfire",
    "asp2020": "landfire", "elev2020": "landfire", "slpd2020": "landfire",
}

# Which subfolders to scan for stacks
SCAN_SUBFOLDERS = [
    "fire_spread",
    "frp",
    "high_res_climate",
    "low_res_climate",
    "fuel_topo",
    "landfire",
]

def _slugify(name: str) -> str:
    """Make a filesystem-safe slug from a fire name."""
    if not name or not isinstance(name, str):
        return ""
    s = name.strip()
    s = re.sub(r"[^\w]+", "_", s, flags=re.UNICODE)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:80]

def _read_tif_bands(tif_path: Path):
    with rasterio.open(tif_path) as ds:
        pix_km2 = abs(ds.transform.a * ds.transform.e) / 1_000_000
        bands = [ds.read(i) for i in range(1, ds.count + 1)]
        nodata = ds.nodata
        crs = ds.crs
        transform = ds.transform
        width, height = ds.width, ds.height
    return bands, pix_km2, nodata, crs, transform, width, height

def _sat_vp(temp_c):
    return np.exp(temp_c * 17.27 / (temp_c + 237.3)) * 0.6108

def _mean_rh(temp_c, dew_c):
    return _sat_vp(dew_c) / _sat_vp(temp_c) * 100.0

def _mean_elevation_series(fline_tif: Path, dem_tif: Path,
                           start_time_utc: pd.Timestamp, freq="h") -> pd.Series:
    from numpy import nanmean
    with rasterio.open(fline_tif) as fsrc, rasterio.open(dem_tif) as dsrc:
        with WarpedVRT(
            dsrc, crs=fsrc.crs, transform=fsrc.transform,
            width=fsrc.width, height=fsrc.height, resampling=Resampling.bilinear
        ) as dem_vrt:
            dem = dem_vrt.read(1).astype("float32")
            if dsrc.nodata is not None:
                dem = np.where(dem == dsrc.nodata, np.nan, dem)
            vals = []
            for i in range(1, fsrc.count + 1):
                mask = fsrc.read(i) > 0
                vals.append(nanmean(dem[mask]) if mask.any() else np.nan)
    idx = pd.date_range(start_time_utc, periods=len(vals), freq=freq, tz="UTC")
    return pd.Series(vals, index=idx)

def _timezone_of(lat: float, lon: float):
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
    s = pd.to_datetime(s)
    if pytz is not None and hasattr(tz_local, "localize"):
        if getattr(s.dt, "tz", None) is not None:
            return s.dt.tz_convert("UTC")
        return s.dt.tz_localize(tz_local, ambiguous="infer", nonexistent="shift_forward").dt.tz_convert("UTC")
    offset = tz_local.utcoffset(None)
    return (s - offset).dt.tz_localize("UTC")

def _has_time_token(s: str) -> bool:
    s = str(s).strip().lower()
    return (":" in s) or ("am" in s) or ("pm" in s)

def _extract_firelist_window(row: pd.Series):
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
            if not _has_time_token(s_raw):
                s = s + pd.Timedelta(hours=1, minutes=30)  # 01:30 local
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
    Returns (frp_density, frp_total, frp_on_fline, frp_density_on_fline)

    - frp_density: mean FRP density (MW/m^2) over VALID and NON-ZERO pixels
    - frp_total: sum over pixels where FRP > 0 (density * pixel_area_m2) → MW
    - frp_on_fline: sum over pixels where (fline>0 and FRP>0) → MW
    - frp_density_on_fline: mean FRP density where (fline>0 and FRP>0) → MW/m^2
    """
    frp_tif = cubes_dir / "frp" / "frp.tif"
    if not frp_tif.exists():
        return None, None, None, None

    fline_tif = cubes_dir / "fire_spread" / "fline.tif"
    with rasterio.open(frp_tif) as frp_ds:
        nb = frp_ds.count
        pix_m2 = abs(frp_ds.transform.a * frp_ds.transform.e)

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

            dens_vals.append(np.nanmean(np.where(frp_band > 0, frp_band, np.nan)))

            vals = np.where(frp_band > 0, frp_band, np.nan)
            total = np.nansum(vals) * pix_m2 if np.any(~np.isnan(vals)) else np.nan
            total_vals.append(total)

            on_fline = np.nan
            dens_on_fline = np.nan
            if fline_vrt is not None and i <= fline_vrt.count:
                fl_band = fline_vrt.read(i)
                mask = (fl_band > 0) & (frp_band > 0)
                if np.any(mask):
                    # total FRP restricted to fireline, excluding 0s
                    on_fline = np.nansum(np.where(mask, frp_band, np.nan)) * pix_m2
                    # density restricted to fireline, excluding 0s
                    dens_on_fline = np.nanmean(np.where(mask, frp_band, np.nan))
            on_fline_vals.append(on_fline)
            dens_on_fline_vals.append(dens_on_fline)

        start = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
        idx = pd.date_range(start, periods=nb, freq="h", tz="UTC")
        if fline_vrt is not None:
            fline_vrt.close()
            fsrc.close()

        s_dens       = pd.Series(dens_vals, index=idx, name="frp_density")
        s_tot        = pd.Series(total_vals, index=idx, name="frp_total")
        s_fl         = pd.Series(on_fline_vals, index=idx, name="frp_on_fline")
        s_dens_fl    = pd.Series(dens_on_fline_vals, index=idx, name="frp_density_on_fline")
        return s_dens, s_tot, s_fl, s_dens_fl

def build_table(fid: str) -> Path:
    # firelist & timezone
    firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)
    row = firelist[firelist["Event_ID"] == fid].iloc[0]
    fire_name_raw = str(row.get("Incid_Name", "")).strip()
    fire_name_slug = _slugify(fire_name_raw)
    center_lat = float((row.get("lat0") + row.get("lat1")) / 2.0)
    center_lon = float((row.get("lon0") + row.get("lon1")) / 2.0)
    tz_local = _timezone_of(center_lat, center_lon)

    # window from firelist (preferred), else FEDS perimeters
    s_local, e_local, used_cols, s_raw, e_raw = _extract_firelist_window(row)
    if s_local is None:
        gdf_fperim_rd, _, _ = feds_util.read_1fire(fid)
        if gdf_fperim_rd is None:
            raise FileNotFoundError(f"FEDS file for fire {fid} does not exist and firelist has no times.")
        times_local = pd.to_datetime(gdf_fperim_rd.t)                     # LOCAL
        local_buf = proc_util.add_time_buffers(times_local)               # buffer on LOCAL timeline
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

    # paths
    cubes_dir = Path(gen_util.dir_output) / gen_util.dir_cubes / fid

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

    full_utc = pd.date_range(
        grid_start_utc,
        grid_end_utc,
        freq="30min",
        tz="UTC",
    )
    if pytz is not None and hasattr(tz_local, "localize"):
        full_local = [pd.Timestamp(t).tz_convert(tz_local).tz_localize(None) for t in full_utc]
    else:
        full_local = [(pd.Timestamp(t) + tz_local.utcoffset(None)).tz_localize(None) for t in full_utc]
    df = pd.DataFrame({"UTC": full_utc, "Local": full_local}, index=full_utc)

    discovered = {}
    for sub in SCAN_SUBFOLDERS:
        subdir = cubes_dir / sub
        if not subdir.exists():
            continue
        for tif in sorted(subdir.glob("*.tif")):
            var = tif.stem 
            folder = SUBFOLDER_HINT.get(var, sub)
            discovered[var] = folder

    # always ensure canonical FEDS & frp keys are included if present
    for k in ["fline", "fperim", "nfp", "frp"]:
        folder = SUBFOLDER_HINT.get(k)
        if folder and (cubes_dir / folder / f"{k}.tif").exists():
            discovered[k] = folder

    def _series_from_stack(var: str, folder: str):
        tif = cubes_dir / folder / f"{var}.tif"
        if not tif.exists():
            return None, None
        bands, pix_km2, nodata, _, _, _, _ = _read_tif_bands(tif)
        n = len(bands)

        # FEDS stacks: index relative to fire start (t0 + 30min)
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

        # values for band means
        vals = []
        for b in bands:
            bb = b.astype("float32")
            if nodata is not None:
                bb = np.where(bb == nodata, np.nan, bb)
            vals.append(np.nanmean(bb))

        # unit conversions
        if var in {"t2m", "d2m"}:
            vals = [v - 273.15 if v is not None and not np.isnan(v) else v for v in vals]  # K -> °C
        elif var == "tp":
            vals = [v * 1000.0 if v is not None and not np.isnan(v) else v for v in vals]  # m -> mm

        # Static layers (fuel_topo/landfire) or single-band rasters
        is_static = (folder in {"fuel_topo", "landfire"}) or (n == 1)
        if is_static:
            value = float(vals[0]) if len(vals) else np.nan
            s_single = pd.Series([value], index=[df.index[0]])
            return s_single.reindex(df.index), None

        start = grid_start_utc
        idx = pd.date_range(start, periods=n, freq="h", tz="UTC")
        return pd.Series(vals, index=idx), None

    # core FEDS counts + areas
    for var in ("fline", "fperim", "nfp"):
        if discovered.get(var) == "fire_spread":
            s, s_area = _series_from_stack(var, "fire_spread")
            if s is not None:
                df[var] = s.reindex(df.index)
            if s_area is not None:
                df[f"{var}_area_km2"] = s_area.reindex(df.index)

    s_frp_dens, s_frp_tot, s_frp_fline, s_frp_dens_fline = _frp_series(cubes_dir, fire_start_utc)
    if s_frp_dens is not None:
        df["frp_density"] = s_frp_dens.reindex(df.index)                # MW/m^2
    if s_frp_tot is not None:
        df["frp_total"] = s_frp_tot.reindex(df.index)                   # MW
    if s_frp_fline is not None:
        df["frp_on_fline"] = s_frp_fline.reindex(df.index)              # MW
    if s_frp_dens_fline is not None:
        df["frp_density_on_fline"] = s_frp_dens_fline.reindex(df.index) # MW/m^2

    skip_vars = {"fline", "fperim", "nfp", "frp"}
    for var, folder in discovered.items():
        if var in skip_vars:
            continue
        s, s_area = _series_from_stack(var, folder)
        if s is not None:
            df[var] = s.reindex(df.index)
        if s_area is not None:
            df[f"{var}_area_km2"] = s_area.reindex(df.index)

    # derived variables
    if {"t2m", "d2m"}.issubset(df.columns):
        df["vpd"] = _sat_vp(df["t2m"]) - _sat_vp(df["d2m"])
        df["mean_rh"] = _mean_rh(df["t2m"], df["d2m"])

    dem_tif   = cubes_dir / "fuel_topo"   / "dem.tif"
    if fline_tif.exists() and dem_tif.exists():
        fline_start_for_elev = fire_start_utc.floor("h") + pd.Timedelta(minutes=30)
        df["fline_mean_elev_m"] = _mean_elevation_series(
            fline_tif, dem_tif, fline_start_for_elev
        ).reindex(df.index)

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

    # keep areas for FEDS up front if present
    feds_area_cols = [c for c in (f"{v}_area_km2" for v in ("fline","fperim","nfp")) if c in df.columns]
    feds_count_cols = [c for c in ("fline","fperim","nfp") if c in df.columns]
    frp_cols = [c for c in ("frp_density","frp_total","frp_on_fline","frp_density_on_fline") if c in df.columns]

    ordered = (
        ["UTC","Local"]
        + feds_area_cols
        + feds_count_cols
        + frp_cols
        + sorted([v for v in dynamic_vars if v in df.columns])
        + sorted([v for v in static_vars if v in df.columns])
    )

    remaining = [c for c in df.columns if c not in set(ordered)]
    ordered += sorted(remaining)

    df = df[[c for c in ordered if c in df.columns]]

    df_out = df.copy()
    df_out["UTC"] = df_out["UTC"].dt.tz_convert("UTC").dt.tz_localize(None)
    out_dir.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_csv, sep="\t", index=False, float_format="%.6f", date_format="%Y-%m-%d %H:%M:%S")
    print(f"Saved ({fire_name_raw or 'Unknown fire'}):", out_csv)
    return out_csv

# CLI
if __name__ == "__main__":
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
