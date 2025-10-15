from pathlib import Path
import re
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.plot import plotting_extent
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import ListedColormap
import matplotlib.dates as mdates
import os
import util.feds_util as feds_util
import util.general_util as gen_util

def _span_days(lo, hi) -> float:
    dt = pd.to_datetime(hi) - pd.to_datetime(lo)
    return dt.total_seconds() / 86400.0

def _best_loc_and_fmt(lo, hi):
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
    if d <= 10:  return "D"
    if d <= 35:  return "W-MON"
    if d <= 120: return "MS"
    return "2MS"

def _iter_ticks(start, end, freq):
    s = pd.to_datetime(start)
    e = pd.to_datetime(end)
    rng = pd.date_range(start=s.floor("D"), end=e.ceil("D"), freq=freq)
    return rng[(rng >= s) & (rng <= e)]

CURSOR_WIDTH = 2.0

def _slugify(name: str) -> str:
    s = (name or "").strip()
    s = re.sub(r"[^\w]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:80]

def _find_evt_crosswalk(fname: str = "LF2024_EVT.csv") -> Path:
    env = os.getenv("FIRECUBE_EVT_CROSSWALK")
    if env:
        p = Path(env).expanduser().resolve()
        if p.exists():
            return p
    try:
        proj_root = Path(__file__).resolve().parents[1]
    except Exception:
        proj_root = Path.cwd()
    candidates = [
        proj_root / "inputData" / fname,
        proj_root / fname,
        Path.cwd() / fname,
    ]
    for p in candidates:
        if p.exists():
            return p.resolve()
    raise FileNotFoundError(
        f"Could not find {fname}. Tried: "
        + ", ".join(str(c) for c in candidates)
        + " or set FIRECUBE_EVT_CROSSWALK=/full/path/to/LF2024_EVT.csv"
    )

def _make_axis_converters(local_series: pd.Series, utc_series: pd.Series):
    xL = mdates.date2num(pd.to_datetime(local_series))
    xU = mdates.date2num(pd.to_datetime(utc_series))
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

def _read_table_for_fid(fid: str):
    firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)
    row = firelist[firelist["Event_ID"] == fid].iloc[0]
    fire_name = str(row.get("Incid_Name", "")).strip()
    fire_slug = _slugify(fire_name)
    dir_suffix = f"{fid}_{fire_slug}" if fire_slug else fid
    sum_vis_dir = Path(gen_util.dir_output) / "sum_vis" / dir_suffix
    sum_vis_dir.mkdir(parents=True, exist_ok=True)
    csv_path = sum_vis_dir / f"data_table_{dir_suffix}.csv"
    if not csv_path.exists():
        matches = list(sum_vis_dir.glob(f"data_table_{fid}_*.csv"))
        if not matches:
            raise FileNotFoundError(f"No CSV found in {sum_vis_dir} for {fid}")
        csv_path = matches[0]
    df = (pd.read_csv(csv_path, sep="\t", parse_dates=["UTC", "Local"])
            .sort_values("Local")
            .reset_index(drop=True))
    local2utc, utc2local = _make_axis_converters(df["Local"], df["UTC"])
    return df, sum_vis_dir, dir_suffix, fire_name, local2utc, utc2local

def _cube_dirs_for_fid(fid: str):
    root = Path(gen_util.dir_output) / "cubes" / fid
    return {
        "root": root,
        "fuel_topo": root / "fuel_topo",
        "fire_spread": root / "fire_spread",
        "landfire": root / "landfire",
    }

def _ensure_movies_dir(sum_vis_dir: Path):
    out = sum_vis_dir / "movies"
    out.mkdir(parents=True, exist_ok=True)
    return out

def read_raster(path: Path, *, single_band: bool = True):
    with rasterio.open(path) as ds:
        arr = ds.read(1 if single_band else None).astype("float32")
        left, bottom, right, top = ds.bounds
        dx, dy = ds.res
        nodata = ds.nodata
    if nodata is not None:
        arr[arr == nodata] = np.nan
    return arr, (left, bottom, right, top), (dx, dy)

def get_mask_resampled(ds, band_idx, out_h, out_w):
    band = ds.read(
        band_idx,
        out_shape=(1, out_h, out_w),
        resampling=Resampling.nearest
    )
    return (band > 0).astype("uint8")

def infer_km2_per_pixel(dfp: pd.DataFrame, var: str, area_col: str, dx=None, dy=None):
    if area_col in dfp and dfp[area_col].notna().any():
        d_pix  = dfp[var].diff()
        d_area = dfp[area_col].diff()
        m = d_pix.notna() & d_area.notna() & (d_pix != 0)
        r = (d_area[m] / d_pix[m]).to_numpy()
        r = r[np.isfinite(r) & (r > 0)]
        if r.size:
            return float(np.median(r))
        lvl = (dfp[area_col] / dfp[var].replace(0, np.nan)).to_numpy()
        lvl = lvl[np.isfinite(lvl) & (lvl > 0)]
        if lvl.size:
            return float(np.median(lvl))
    if dx is not None and dy is not None:
        return abs(dx * dy) / 1e6
    return None

def red_mask_cmap(alpha: float = 0.9) -> ListedColormap:
    return ListedColormap(["none", (1, 0, 0, alpha)])

def black_mask_cmap(alpha: float = 0.6) -> ListedColormap:
    return ListedColormap(["none", (0, 0, 0, alpha)])

def _find_evt_raster(landfire_dir: Path) -> tuple[Path, str]:
    if not landfire_dir.exists():
        raise FileNotFoundError(f"landfire dir not found: {landfire_dir}")

    # 1) exact
    p = landfire_dir / "200evt.tif"
    if p.exists():
        return p, p.stem

    # 2) any *evt*.tif
    candidates = sorted(landfire_dir.glob("*evt*.tif"))
    if candidates:
        return candidates[0], candidates[0].stem

    # 3) fallback to any tif
    any_tifs = sorted(landfire_dir.glob("*.tif"))
    if any_tifs:
        return any_tifs[0], any_tifs[0].stem

    raise FileNotFoundError(f"No EVT-like tif found in {landfire_dir}")

def build_evt_colormap(evt_arr, xwalk_csv: Path):
    df = pd.read_csv(xwalk_csv, dtype={"VALUE": "int64"}, engine="python")
    df = df[df["VALUE"] > 0].set_index("VALUE")
    lifeforms = df["EVT_LF"].unique()
    def to_hex(r, g, b):
        return f"#{int(r):02x}{int(g):02x}{int(b):02x}"
    palette = {lf: to_hex(*df.loc[df["EVT_LF"] == lf, ["R", "G", "B"]].median()) for lf in lifeforms}
    palette["Unknown"] = "#000000"
    code2lf  = df["EVT_LF"].to_dict()
    lf_array = np.vectorize(lambda c: code2lf.get(int(c), "Unknown"))(evt_arr)
    lf_labels, lf_index = np.unique(lf_array, return_inverse=True)
    lf_index = lf_index.reshape(evt_arr.shape)
    colors   = [palette.get(lf, "#000000") for lf in lf_labels]
    cmap_evt = ListedColormap(colors, name="evt_lf")
    return lf_labels, colors, lf_index, cmap_evt

def _derive_time_vectors(dfp: pd.DataFrame, n_frames: int, local2utc):
    start_local, end_local = dfp["Local"].iloc[[0, -1]]
    start_num = mdates.date2num(start_local)
    end_num   = mdates.date2num(end_local)
    times_local_num = np.linspace(start_num, end_num, n_frames)
    times_local = mdates.num2date(times_local_num)
    times_utc_num = local2utc(times_local_num)
    times_utc = mdates.num2date(times_utc_num)
    times_utc_str = [pd.to_datetime(t).strftime("%b %d %H:%M UTC") for t in times_utc]
    return (start_local, end_local, times_local, times_utc, times_utc_str)

def _pad_xlim(lo_local, hi_local):
    lo = pd.to_datetime(lo_local); hi = pd.to_datetime(hi_local)
    if lo == hi: hi = lo + pd.Timedelta(minutes=30)
    return lo, hi

def _secondary_time_axis(ax, local2utc, utc2local, lo_local, hi_local):
    lo, hi = _pad_xlim(lo_local, hi_local)
    ax.set_xlim(lo, hi)
    locL, fmtL = _best_loc_and_fmt(lo, hi)
    ax.xaxis.set_major_locator(locL)
    ax.xaxis.set_major_formatter(fmtL)
    secax = ax.secondary_xaxis("top", functions=(local2utc, utc2local))
    locU, fmtU = _best_loc_and_fmt(lo, hi)
    secax.xaxis.set_major_locator(locU)
    secax.xaxis.set_major_formatter(fmtU)
    return secax

def make_trio_over_dem(fid: str, df: pd.DataFrame, local2utc, utc2local,
                       movies_dir: Path, fire_name: str, cube_dirs: dict):
    dem_path = cube_dirs["fuel_topo"] / "dem.tif"
    spread_dir = cube_dirs["fire_spread"]
    fpaths = {
        "fline":  spread_dir / "fline.tif",
        "fperim": spread_dir / "fperim.tif",
        "nfp":    spread_dir / "nfp.tif",
    }
    if not dem_path.exists():
        print(f"[skip] DEM not found: {dem_path}")
        return
    dem, (l, b, r, t), (dx, dy) = read_raster(dem_path, single_band=True)
    dx_half, dy_half = dx / 2.0, abs(dy) / 2.0
    dem_extent  = [l - dx_half, r + dx_half, b - dy_half, t + dy_half]

    for var, raster_path in fpaths.items():
        if not raster_path.exists():
            print(f"[skip] {var}: raster not found → {raster_path}")
            continue
        stack, _, _ = read_raster(raster_path, single_band=False)  # (T,H,W)
        n_frames = stack.shape[0]
        area_col = f"{var}_area_km2"
        cols = ["UTC","Local",var] + ([area_col] if area_col in df.columns else [])
        dfp = df[cols].dropna(subset=[var]).copy()
        if dfp.empty:
            print(f"[skip] {var}: only NaNs in time series")
            continue

        start_local, end_local, times_local, times_utc, times_utc_str = \
            _derive_time_vectors(dfp, n_frames, local2utc)

        fig = plt.figure(figsize=(13.5, 10), constrained_layout=True)
        gs  = fig.add_gridspec(nrows=2, ncols=21, height_ratios=[3.1, 2.4])
        ax_map   = fig.add_subplot(gs[0, 0:20])
        cax      = fig.add_subplot(gs[0, 20])
        ax_ts    = fig.add_subplot(gs[1, :])

        dem_im = ax_map.imshow(dem, cmap="terrain", extent=dem_extent, origin="upper")
        cbar = plt.colorbar(dem_im, cax=cax); cbar.set_label("Elevation (m)")
        ax_map.set_aspect("equal")
        ax_map.set_xlabel("Easting (map units)")
        ax_map.set_ylabel("Northing (map units)")

        overlay = ax_map.imshow(
            stack[0], cmap=red_mask_cmap(0.6), vmin=0, vmax=1,
            extent=dem_extent, origin="upper",
        )
        title = ax_map.set_title(
            f"{fid} – {fire_name} — {var} on DEM — {times_utc_str[0]}",
            fontsize=16, pad=10
        )

        ax_ts.set_axisbelow(True); ax_ts.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_ts.step(dfp["Local"], dfp[var], where="post", linewidth=2.0, zorder=3)
        value_changes = dfp[var].ne(dfp[var].shift())
        ax_ts.scatter(
            dfp.loc[value_changes, "Local"], dfp.loc[value_changes, var],
            color="black", zorder=4, s=30
        )

        km2_per_pixel = infer_km2_per_pixel(dfp, var, area_col, dx, dy)
        if km2_per_pixel and np.isfinite(km2_per_pixel) and km2_per_pixel > 0:
            sec_y = ax_ts.secondary_yaxis('right',
                functions=(lambda y: y * km2_per_pixel, lambda y: y / km2_per_pixel))
            sec_y.set_ylabel("Area (km²)", labelpad=8)

        ax_ts.set_ylabel("Number of pixels", labelpad=6)
        ax_ts.set_title(f"{fire_name} – {var}", pad=8)

        secax = _secondary_time_axis(ax_ts, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        ax_ts.set_xlabel(
            f"Local — {start_local:%b %d %H:%M} → {end_local:%b %d %H:%M}",
            labelpad=8
        )

        gfreq = _grid_freq(df["Local"].iloc[0], df["Local"].iloc[-1])
        for dt in _iter_ticks(df["Local"].iloc[0], df["Local"].iloc[-1], gfreq):
            ax_ts.axvline(dt, linestyle=":", color="grey", linewidth=0.8, zorder=1)

        fig.canvas.draw()
        ymin, ymax = ax_ts.get_ylim()
        cursor, = ax_ts.plot([times_local[0], times_local[0]], [ymin, ymax],
                             linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)

        def init():
            overlay.set_data(stack[0])
            y0, y1 = ax_ts.get_ylim()
            cursor.set_data([times_local[0], times_local[0]], [y0, y1])
            title.set_text(f"{fid} – {fire_name} — {var} on DEM — {times_utc_str[0]}")
            return overlay, cursor, title

        def update(i: int):
            overlay.set_data(stack[i])
            title.set_text(f"{fid} – {fire_name} — {var} on DEM — {times_utc_str[i]}")
            y0, y1 = ax_ts.get_ylim()
            cursor.set_data([times_local[i], times_local[i]], [y0, y1])
            return overlay, cursor, title

        ani = animation.FuncAnimation(fig, update, init_func=init, frames=n_frames, interval=300, blit=True)

        out_mp4 = movies_dir / f"{var}_on_dem_with_linear.mp4"
        try:
            writer = animation.FFMpegWriter(fps=10, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
            ani.save(out_mp4, writer=writer)
        except (RuntimeError, FileNotFoundError):
            ani.save(out_mp4, writer=animation.FFMpegWriter(fps=10, codec="mpeg4"))
        plt.close(fig)
        print("MP4 saved:", out_mp4)

def make_trio_over_evt(fid: str, df: pd.DataFrame, local2utc, utc2local,
                       movies_dir: Path, fire_name: str, cube_dirs: dict):
    landfire_dir = cube_dirs["landfire"]
    try:
        evt_path, evt_stem = _find_evt_raster(landfire_dir)
    except FileNotFoundError as e:
        print(f"[skip] EVT not found: {e}")
        return

    xwalk_csv = _find_evt_crosswalk()
    spread_dir = cube_dirs["fire_spread"]
    fpaths = {
        "fline":  spread_dir / "fline.tif",
        "fperim": spread_dir / "fperim.tif",
        "nfp":    spread_dir / "nfp.tif",
    }

    with rasterio.open(evt_path) as ds_evt:
        evt = ds_evt.read(1)
        evt_nodata = ds_evt.nodata
        if evt_nodata is not None:
            evt = np.where(evt == evt_nodata, 0, evt)
        extent = plotting_extent(ds_evt)
        evt_dx, evt_dy = ds_evt.res
        H, W = evt.shape

    try:
        lf_labels, lf_colors, lf_index, cmap_evt = build_evt_colormap(evt, xwalk_csv)
    except Exception as e:
        print(f"[warn] Crosswalk not applied for {evt_stem} ({xwalk_csv}): {e}")
        lf_index = np.clip(evt, 0, np.nanmax(evt)).astype(int)
        cmap_evt, lf_labels, lf_colors = None, None, None

    var_to_cmap = {
        "fperim": black_mask_cmap(0.6),
        "fline":  red_mask_cmap(0.9),
        "nfp":    red_mask_cmap(0.6),
    }

    for var, raster_path in fpaths.items():
        if not raster_path.exists():
            print(f"[skip] {var}: raster not found → {raster_path}")
            continue

        ds_var = rasterio.open(raster_path)
        n_frames = ds_var.count

        area_col = f"{var}_area_km2"
        cols = ["UTC", "Local", var] + ([area_col] if area_col in df.columns else [])
        dfp = df[cols].dropna(subset=[var]).copy()
        if dfp.empty:
            print(f"[skip] {var}: only NaNs in time series")
            ds_var.close()
            continue

        start_local, end_local, times_local, times_utc, times_utc_str = \
            _derive_time_vectors(dfp, n_frames, local2utc)

        fig = plt.figure(figsize=(14.5, 10.5), constrained_layout=True)
        gs  = fig.add_gridspec(nrows=2, ncols=25, height_ratios=[3.0, 2.4])

        ax_map  = fig.add_subplot(gs[0, 0:20])
        ax_leg  = fig.add_subplot(gs[0, 20:25])
        ax_ts   = fig.add_subplot(gs[1, :])

        im_evt = ax_map.imshow(
            lf_index, cmap=cmap_evt, extent=extent,
            vmin=0, vmax=(np.nanmax(lf_index) if np.isfinite(lf_index).any() else 1),
            alpha=1.0, origin="upper"
        )

        ax_leg.axis("off")
        if lf_labels is not None and lf_colors is not None:
            from matplotlib.patches import Patch
            patches = [Patch(color=lf_colors[i], label=lf_labels[i]) for i in range(len(lf_labels))]
            leg = ax_leg.legend(handles=patches, loc="upper left", frameon=True, fontsize=9)
            ax_leg.add_artist(leg)

        var_cmap = var_to_cmap.get(var, red_mask_cmap(0.6))
        im_overlay = ax_map.imshow(
            np.zeros_like(evt), cmap=var_cmap, vmin=0, vmax=1,
            alpha=0.9, extent=extent, origin="upper"
        )

        ax_map.set_title(f"{fid} – {fire_name} — {var} on {evt_stem} — {times_utc_str[0]}", fontsize=15, pad=8)
        ax_map.set_xlabel("Easting")
        ax_map.set_ylabel("Northing")

        ax_ts.set_axisbelow(True)
        ax_ts.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_ts.step(dfp["Local"], dfp[var], where="post", linewidth=2.0, zorder=3)
        value_changes = dfp[var].ne(dfp[var].shift())
        ax_ts.scatter(dfp.loc[value_changes, "Local"], dfp.loc[value_changes, var],
                      color="black", zorder=4, s=30)

        km2_per_pixel = infer_km2_per_pixel(dfp, var, area_col, evt_dx, evt_dy)
        if km2_per_pixel and np.isfinite(km2_per_pixel) and km2_per_pixel > 0:
            sec_y = ax_ts.secondary_yaxis(
                'right',
                functions=(lambda y: y * km2_per_pixel, lambda y: y / km2_per_pixel)
            )
            sec_y.set_ylabel("Area (km²)", labelpad=8)

        ax_ts.set_ylabel("Number of pixels", labelpad=6)
        ax_ts.set_title(f"{fire_name} – {var}", pad=8)

        secax = _secondary_time_axis(ax_ts, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        ax_ts.set_xlabel(
            f"Local — {start_local:%b %d %H:%M} → {end_local:%b %d %H:%M}",
            labelpad=8
        )

        gfreq = _grid_freq(df["Local"].iloc[0], df["Local"].iloc[-1])
        for dt in _iter_ticks(df["Local"].iloc[0], df["Local"].iloc[-1], gfreq):
            ax_ts.axvline(dt, linestyle=":", color="grey", linewidth=0.8, zorder=1)

        fig.canvas.draw()
        ymin, ymax = ax_ts.get_ylim()
        cursor, = ax_ts.plot([times_local[0], times_local[0]], [ymin, ymax],
                             linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)

        def init():
            m0 = get_mask_resampled(rasterio.open(raster_path), 1, H, W)
            im_overlay.set_data(m0)
            y0, y1 = ax_ts.get_ylim()
            cursor.set_data([times_local[0], times_local[0]], [y0, y1])
            ax_map.set_title(f"{fid} – {fire_name} — {var} on {evt_stem} — {times_utc_str[0]}")
            return im_overlay, cursor

        def update(i: int):
            m = get_mask_resampled(rasterio.open(raster_path), i + 1, H, W)
            im_overlay.set_data(m)
            ax_map.set_title(f"{fid} – {fire_name} — {var} on {evt_stem} — {times_utc_str[i]}")
            y0, y1 = ax_ts.get_ylim()
            cursor.set_data([times_local[i], times_local[i]], [y0, y1])
            return im_overlay, cursor

        ani = animation.FuncAnimation(fig, update, init_func=init, frames=n_frames, interval=300, blit=True)

        out_mp4 = movies_dir / f"{var}_on_{evt_stem}_with_linear.mp4"
        try:
            writer = animation.FFMpegWriter(fps=10, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
            ani.save(out_mp4, writer=writer)
        except (RuntimeError, FileNotFoundError):
            ani.save(out_mp4, writer=animation.FFMpegWriter(fps=10, codec="mpeg4"))
        plt.close(fig)
        ds_var.close()
        print("MP4 saved:", out_mp4)

def make_combo_over_evt(fid: str, df: pd.DataFrame, local2utc, utc2local,
                        movies_dir: Path, fire_name: str, cube_dirs: dict):
    landfire_dir = cube_dirs["landfire"]
    try:
        evt_path, evt_stem = _find_evt_raster(landfire_dir)
    except FileNotFoundError as e:
        print("[skip] combo EVT:", e)
        return

    xwalk_csv = _find_evt_crosswalk()
    spread_dir = cube_dirs["fire_spread"]
    fline_path = spread_dir / "fline.tif"
    fperim_path = spread_dir / "fperim.tif"

    if not (evt_path.exists() and fline_path.exists() and fperim_path.exists()):
        print("[skip] combo EVT: missing one of EVT/fline/fperim rasters")
        return

    with rasterio.open(evt_path) as ds_evt:
        evt = ds_evt.read(1)
        evt_nodata = ds_evt.nodata
        if evt_nodata is not None:
            evt = np.where(evt == evt_nodata, 0, evt)
        extent = plotting_extent(ds_evt)
        evt_dx, evt_dy = ds_evt.res
        H, W = evt.shape

    try:
        lf_labels, lf_colors, lf_index, cmap_evt = build_evt_colormap(evt, xwalk_csv)
    except Exception as e:
        print(f"[warn] Crosswalk not applied for {evt_stem} ({xwalk_csv}): {e}")
        lf_index = np.clip(evt, 0, np.nanmax(evt)).astype(int)
        cmap_evt = None
        lf_labels, lf_colors = None, None

    fline_ds  = rasterio.open(fline_path)
    fperim_ds = rasterio.open(fperim_path)
    n_frames  = min(fline_ds.count, fperim_ds.count)

    def prep_df(var: str):
        area_col = f"{var}_area_km2"
        cols = ["UTC", "Local", var] + ([area_col] if area_col in df.columns else [])
        d = df[cols].dropna(subset=[var]).copy()
        return d, area_col

    dfp_fp, area_fp = prep_df("fperim")
    dfp_fl, area_fl = prep_df("fline")
    if dfp_fp.empty and dfp_fl.empty:
        print("[skip] combo EVT: both fperim and fline time series empty")
        fline_ds.close(); fperim_ds.close()
        return

    spans = []
    if not dfp_fp.empty: spans.append(dfp_fp["Local"].iloc[[0, -1]].to_list())
    if not dfp_fl.empty: spans.append(dfp_fl["Local"].iloc[[0, -1]].to_list())
    start_local = min(s for s, _ in spans)
    end_local   = max(e for _, e in spans)

    start_num = mdates.date2num(start_local)
    end_num   = mdates.date2num(end_local)
    times_local_num = np.linspace(start_num, end_num, n_frames)
    times_local = mdates.num2date(times_local_num)
    times_utc_num = local2utc(times_local_num)
    times_utc = mdates.num2date(times_utc_num)
    times_utc_str = [pd.to_datetime(t).strftime("%b %d %H:%M UTC") for t in times_utc]

    fig = plt.figure(figsize=(14.8, 12.0), constrained_layout=True)
    gs  = fig.add_gridspec(nrows=3, ncols=25, height_ratios=[3.1, 2.2, 2.2])

    ax_map = fig.add_subplot(gs[0, 0:21])
    ax_leg = fig.add_subplot(gs[0, 21:25])
    ax_fp  = fig.add_subplot(gs[1, :])
    ax_fl  = fig.add_subplot(gs[2, :])

    ax_map.imshow(lf_index, cmap=cmap_evt, extent=extent,
                  vmin=0, vmax=(np.nanmax(lf_index) if np.isfinite(lf_index).any() else 1),
                  alpha=1.0, origin="upper")

    ax_leg.axis("off")
    if lf_labels is not None and lf_colors is not None:
        from matplotlib.patches import Patch
        patches = [Patch(color=lf_colors[i], label=lf_labels[i]) for i in range(len(lf_labels))]
        leg = ax_leg.legend(handles=patches, loc="upper left", frameon=True, fontsize=9)
        ax_leg.add_artist(leg)

    im_fperim = ax_map.imshow(np.zeros_like(evt), cmap=black_mask_cmap(0.6), vmin=0, vmax=1,
                              alpha=0.8, extent=extent, origin="upper")
    im_fline  = ax_map.imshow(np.zeros_like(evt), cmap=red_mask_cmap(0.9), vmin=0, vmax=1,
                              alpha=1.0, extent=extent, origin="upper")

    title = ax_map.set_title(
        f"{fid} – {fire_name} — fperim & fline on {evt_stem} — {times_utc_str[0]}",
        fontsize=16, pad=10
    )

    # fperim TS
    if not dfp_fp.empty:
        ax_fp.set_axisbelow(True); ax_fp.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_fp.step(dfp_fp["Local"], dfp_fp["fperim"], where="post", linewidth=2.0, zorder=3)
        vc = dfp_fp["fperim"].ne(dfp_fp["fperim"].shift())
        ax_fp.scatter(dfp_fp.loc[vc, "Local"], dfp_fp.loc[vc, "fperim"], color="black", zorder=4, s=30)
        k_fp = infer_km2_per_pixel(dfp_fp, "fperim", area_fp, evt_dx, evt_dy)
        if k_fp and np.isfinite(k_fp) and k_fp > 0:
            y2 = ax_fp.secondary_yaxis('right', functions=(lambda y: y * k_fp, lambda y: y / k_fp))
            y2.set_ylabel("Area (km²)", labelpad=8)
        secax_fp = _secondary_time_axis(ax_fp, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax_fp.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        sL, eL = dfp_fp["Local"].iloc[[0, -1]]
        ax_fp.set_xlabel(f"Local — {sL:%b %d %H:%M} → {eL:%b %d %H:%M}", labelpad=8)
        ax_fp.set_ylabel("Number of pixels", labelpad=6)
        ax_fp.set_title(f"{fire_name} – fperim", pad=8)

    # fline TS
    if not dfp_fl.empty:
        ax_fl.set_axisbelow(True); ax_fl.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_fl.step(dfp_fl["Local"], dfp_fl["fline"], where="post", linewidth=2.0, zorder=3)
        vc = dfp_fl["fline"].ne(dfp_fl["fline"].shift())
        ax_fl.scatter(dfp_fl.loc[vc, "Local"], dfp_fl.loc[vc, "fline"], color="black", zorder=4, s=30)
        k_fl = infer_km2_per_pixel(dfp_fl, "fline", area_fl, evt_dx, evt_dy)
        if k_fl and np.isfinite(k_fl) and k_fl > 0:
            y2 = ax_fl.secondary_yaxis('right', functions=(lambda y: y * k_fl, lambda y: y / k_fl))
            y2.set_ylabel("Area (km²)", labelpad=8)
        secax_fl = _secondary_time_axis(ax_fl, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax_fl.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        sL, eL = dfp_fl["Local"].iloc[[0, -1]]
        ax_fl.set_xlabel(f"Local — {sL:%b %d %H:%M} → {eL:%b %d %H:%M}", labelpad=8)
        ax_fl.set_ylabel("Number of pixels", labelpad=6)
        ax_fl.set_title(f"{fire_name} – fline", pad=8)

    fig.canvas.draw()
    ymin1, ymax1 = ax_fp.get_ylim()
    ymin2, ymax2 = ax_fl.get_ylim()
    cursor_fp, = ax_fp.plot([times_local[0], times_local[0]], [ymin1, ymax1],
                            linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)
    cursor_fl, = ax_fl.plot([times_local[0], times_local[0]], [ymin2, ymax2],
                            linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)

    def init():
        im_fperim.set_data(get_mask_resampled(fperim_ds, 1, H, W))
        im_fline.set_data(get_mask_resampled(fline_ds,  1, H, W))
        y0a, y1a = ax_fp.get_ylim(); cursor_fp.set_data([times_local[0], times_local[0]], [y0a, y1a])
        y0b, y1b = ax_fl.get_ylim(); cursor_fl.set_data([times_local[0], times_local[0]], [y0b, y1b])
        title.set_text(f"{fid} – {fire_name} — fperim & fline on {evt_stem} — {times_utc_str[0]}")
        return im_fperim, im_fline, cursor_fp, cursor_fl, title

    def update(i: int):
        im_fperim.set_data(get_mask_resampled(fperim_ds, i + 1, H, W))
        im_fline.set_data(get_mask_resampled(fline_ds,  i + 1, H, W))
        y0a, y1a = ax_fp.get_ylim(); cursor_fp.set_data([times_local[i], times_local[i]], [y0a, y1a])
        y0b, y1b = ax_fl.get_ylim(); cursor_fl.set_data([times_local[i], times_local[i]], [y0b, y1b])
        title.set_text(f"{fid} – {fire_name} — fperim & fline on {evt_stem} — {times_utc_str[i]}")
        return im_fperim, im_fline, cursor_fp, cursor_fl, title

    ani = animation.FuncAnimation(fig, update, init_func=init, frames=n_frames, interval=300, blit=True)

    out_mp4 = movies_dir / f"fperim_fline_on_{evt_stem}_with_linear.mp4"
    try:
        writer = animation.FFMpegWriter(fps=10, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
        ani.save(out_mp4, writer=writer)
    except (RuntimeError, FileNotFoundError):
        ani.save(out_mp4, writer=animation.FFMpegWriter(fps=10, codec="mpeg4"))
    plt.close(fig)
    fline_ds.close()
    fperim_ds.close()
    print("MP4 saved:", out_mp4)

def make_combo_over_dem(fid: str, df: pd.DataFrame, local2utc, utc2local,
                        movies_dir: Path, fire_name: str, cube_dirs: dict):
    dem_path = cube_dirs["fuel_topo"] / "dem.tif"
    spread_dir = cube_dirs["fire_spread"]
    fline_path = spread_dir / "fline.tif"
    fperim_path = spread_dir / "fperim.tif"

    if not (dem_path.exists() and fline_path.exists() and fperim_path.exists()):
        print("[skip] combo DEM: missing one of DEM/fline/fperim rasters")
        return

    dem, (l, b, r, t), (dx, dy) = read_raster(dem_path, single_band=True)
    dx_half, dy_half = dx / 2.0, abs(dy) / 2.0
    dem_extent  = [l - dx_half, r + dx_half, b - dy_half, t + dy_half]

    fline_stack,  _, _ = read_raster(fline_path,  single_band=False)
    fperim_stack, _, _ = read_raster(fperim_path, single_band=False)
    n_frames = int(min(fline_stack.shape[0], fperim_stack.shape[0]))

    def prep_df(var: str):
        area_col = f"{var}_area_km2"
        cols = ["UTC", "Local", var] + ([area_col] if area_col in df.columns else [])
        d = df[cols].dropna(subset=[var]).copy()
        return d, area_col

    dfp_fp, area_fp = prep_df("fperim")
    dfp_fl, area_fl = prep_df("fline")
    if dfp_fp.empty and dfp_fl.empty:
        print("[skip] combo DEM: both fperim and fline time series empty")
        return

    spans = []
    if not dfp_fp.empty: spans.append(dfp_fp["Local"].iloc[[0, -1]].to_list())
    if not dfp_fl.empty: spans.append(dfp_fl["Local"].iloc[[0, -1]].to_list())
    start_local = min(s for s, _ in spans)
    end_local   = max(e for _, e in spans)

    start_num = mdates.date2num(start_local)
    end_num   = mdates.date2num(end_local)
    times_local_num = np.linspace(start_num, end_num, n_frames)
    times_local = mdates.num2date(times_local_num)
    times_utc_num = local2utc(times_local_num)
    times_utc = mdates.num2date(times_utc_num)
    times_utc_str = [pd.to_datetime(t).strftime("%b %d %H:%M UTC") for t in times_utc]

    fig = plt.figure(figsize=(14.2, 12.0), constrained_layout=True)
    gs  = fig.add_gridspec(nrows=3, ncols=21, height_ratios=[3.2, 2.2, 2.2])

    ax_map   = fig.add_subplot(gs[0, 0:20])
    cax      = fig.add_subplot(gs[0, 20])
    ax_fp_ts = fig.add_subplot(gs[1, :])
    ax_fl_ts = fig.add_subplot(gs[2, :])

    dem_im = ax_map.imshow(dem, cmap="terrain", extent=dem_extent, origin="upper")
    cbar = plt.colorbar(dem_im, cax=cax); cbar.set_label("Elevation (m)")
    ax_map.set_aspect("equal")
    ax_map.set_xlabel("Easting (map units)")
    ax_map.set_ylabel("Northing (map units)")

    im_fperim = ax_map.imshow(fperim_stack[0], cmap=black_mask_cmap(0.6), vmin=0, vmax=1,
                              extent=dem_extent, origin="upper")
    im_fline = ax_map.imshow(fline_stack[0], cmap=red_mask_cmap(0.9), vmin=0, vmax=1,
                             extent=dem_extent, origin="upper")

    title = ax_map.set_title(
        f"{fid} – {fire_name} — fperim & fline on DEM — {times_utc_str[0]}",
        fontsize=16, pad=10
    )

    # fperim TS
    if not dfp_fp.empty:
        ax_fp_ts.set_axisbelow(True); ax_fp_ts.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_fp_ts.step(dfp_fp["Local"], dfp_fp["fperim"], where="post", linewidth=2.0, zorder=3)
        vc = dfp_fp["fperim"].ne(dfp_fp["fperim"].shift())
        ax_fp_ts.scatter(dfp_fp.loc[vc, "Local"], dfp_fp.loc[vc, "fperim"], color="black", zorder=4, s=30)
        k_fp = infer_km2_per_pixel(dfp_fp, "fperim", area_fp, dx, dy)
        if k_fp and np.isfinite(k_fp) and k_fp > 0:
            y2 = ax_fp_ts.secondary_yaxis('right', functions=(lambda y: y * k_fp, lambda y: y / k_fp))
            y2.set_ylabel("Area (km²)", labelpad=8)
        secax_fp = _secondary_time_axis(ax_fp_ts, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax_fp.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        sL, eL = dfp_fp["Local"].iloc[[0, -1]]
        ax_fp_ts.set_xlabel(f"Local — {sL:%b %d %H:%M} → {eL:%b %d %H:%M}", labelpad=8)
        ax_fp_ts.set_ylabel("Number of pixels", labelpad=6)
        ax_fp_ts.set_title(f"{fire_name} – fperim", pad=8)

    # fline TS
    if not dfp_fl.empty:
        ax_fl_ts.set_axisbelow(True); ax_fl_ts.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax_fl_ts.step(dfp_fl["Local"], dfp_fl["fline"], where="post", linewidth=2.0, zorder=3)
        vc = dfp_fl["fline"].ne(dfp_fl["fline"].shift())
        ax_fl_ts.scatter(dfp_fl.loc[vc, "Local"], dfp_fl.loc[vc, "fline"], color="black", zorder=4, s=30)
        k_fl = infer_km2_per_pixel(dfp_fl, "fline", area_fl, dx, dy)
        if k_fl and np.isfinite(k_fl) and k_fl > 0:
            y2 = ax_fl_ts.secondary_yaxis('right', functions=(lambda y: y * k_fl, lambda y: y / k_fl))
            y2.set_ylabel("Area (km²)", labelpad=8)
        secax_fl = _secondary_time_axis(ax_fl_ts, local2utc, utc2local, df["Local"].iloc[0], df["Local"].iloc[-1])
        secax_fl.set_xlabel(
            f"UTC — {times_utc[0].strftime('%b %d %H:%M')} → {times_utc[-1].strftime('%b %d %H:%M')}",
            labelpad=8
        )
        sL, eL = dfp_fl["Local"].iloc[[0, -1]]
        ax_fl_ts.set_xlabel(f"Local — {sL:%b %d %H:%M} → {eL:%b %d %H:%M}", labelpad=8)
        ax_fl_ts.set_ylabel("Number of pixels", labelpad=6)
        ax_fl_ts.set_title(f"{fire_name} – fline", pad=8)

    fig.canvas.draw()
    ymin1, ymax1 = ax_fp_ts.get_ylim()
    ymin2, ymax2 = ax_fl_ts.get_ylim()
    cursor_fp, = ax_fp_ts.plot([times_local[0], times_local[0]], [ymin1, ymax1],
                               linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)
    cursor_fl, = ax_fl_ts.plot([times_local[0], times_local[0]], [ymin2, ymax2],
                               linewidth=CURSOR_WIDTH, color="tab:orange", alpha=0.95, zorder=5)

    def init():
        im_fperim.set_data(fperim_stack[0])
        im_fline.set_data(fline_stack[0])
        y0a, y1a = ax_fp_ts.get_ylim(); cursor_fp.set_data([times_local[0], times_local[0]], [y0a, y1a])
        y0b, y1b = ax_fl_ts.get_ylim(); cursor_fl.set_data([times_local[0], times_local[0]], [y0b, y1b])
        title.set_text(f"{fid} – {fire_name} — fperim & fline on DEM — {times_utc_str[0]}")
        return im_fperim, im_fline, cursor_fp, cursor_fl, title

    def update(i: int):
        im_fperim.set_data(fperim_stack[i])
        im_fline.set_data(fline_stack[i])
        y0a, y1a = ax_fp_ts.get_ylim(); cursor_fp.set_data([times_local[i], times_local[i]], [y0a, y1a])
        y0b, y1b = ax_fl_ts.get_ylim(); cursor_fl.set_data([times_local[i], times_local[i]], [y0b, y1b])
        title.set_text(f"{fid} – {fire_name} — fperim & fline on DEM — {times_utc_str[i]}")
        return im_fperim, im_fline, cursor_fp, cursor_fl, title

    ani = animation.FuncAnimation(fig, update, init_func=init, frames=n_frames, interval=300, blit=True)

    out_mp4 = movies_dir / "fperim_fline_on_dem_with_linear.mp4"
    try:
        writer = animation.FFMpegWriter(fps=10, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
        ani.save(out_mp4, writer=writer)
    except (RuntimeError, FileNotFoundError):
        ani.save(out_mp4, writer=animation.FFMpegWriter(fps=10, codec="mpeg4"))
    plt.close(fig)
    print("MP4 saved:", out_mp4)

def make_movies(fid: str):
    df, sum_vis_dir, dir_suffix, fire_name, local2utc, utc2local = _read_table_for_fid(fid)
    movies_dir = _ensure_movies_dir(sum_vis_dir)
    cube_dirs  = _cube_dirs_for_fid(fid)

    print(f"Fire: {fid} – {fire_name or '(unknown)'}")
    print(f"CSV:  {sum_vis_dir}/data_table_{dir_suffix}.csv")
    print(f"Cube: {cube_dirs['root']}")
    print(f"Out:  {movies_dir}")

    # A) Three movies over DEM
    make_trio_over_dem(fid, df, local2utc, utc2local, movies_dir, fire_name, cube_dirs)

    # B) Combined DEM movie
    make_combo_over_dem(fid, df, local2utc, utc2local, movies_dir, fire_name, cube_dirs)

    # C) Three movies over EVT (auto-discovered)
    make_trio_over_evt(fid, df, local2utc, utc2local, movies_dir, fire_name, cube_dirs)

    # D) Combined EVT movie (auto-discovered)
    make_combo_over_evt(fid, df, local2utc, utc2local, movies_dir, fire_name, cube_dirs)

# CLI
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python src/make_movies.py <FIRE_EVENT_ID>")
        raise SystemExit(2)
    make_movies(sys.argv[1])
