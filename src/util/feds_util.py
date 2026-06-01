from datetime import timedelta
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import rasterio
import shutil
from rasterio.features import rasterize
from rasterio.transform import from_origin

import util.general_util as gen_util

# Supress SettingWithCopyWarning
pd.set_option('mode.chained_assignment', None)

dir_feds25 = os.path.join('inputData', 'Full_FEDS')
dir_firepix = os.path.join('inputData', 'firepix')
feds_firelist = os.path.join('inputData', 'fireslist2012-2024_withtype_no_fline_hawaii.csv')

def set_gdffile(event_id):
    """
    Find the GeoPackage file corresponding to the given fire event ID in the FEDS input folder.

    Args:
        event_id (str): Fire event ID to search for.

    Returns:
        str: Path to the GeoPackage file for the desired fire if it exists; None otherwise.
    """
    target_file = f'{event_id}.gpkg'

    # For all years under FEDS folder, try to find file corresponding to given event ID
    for root, dirs, files in os.walk(dir_feds25):
        if target_file in files:
            return os.path.join(root, target_file)
    return None

def check_gdffile(event_id):
    """
    Checks whether a GeoPackage file for the given fire event ID exists.

    Args:
        event_id (str): Fire event ID whose file to search for.

    Returns:
        bool: True if a GeoPackage file exists; False otherwise.
    """
    fnm = set_gdffile(event_id)
    return os.path.exists(os.path.expanduser(fnm))

def set_dd(layer="perimeter"):
    """
    Set the diagnostic data for the given layer. This is in addition to the geometries supplied in a
    GeoPackage file or GeoDataFrame.

    Args:
        layer (str, optional): Layer for which diagnostic data needs to be set. Defaults to "perimeter".

    Returns:
        dict: Dictionary mapping diagnostic data names to their types.
    """

    # diagnostic data name and types saved in GeoJSON files (in addition to geometries)
    if layer == "perimeter":
        dd = {
            "t": "datetime64[ns]",
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "sumFRP": "float",  # mean FRP of the new fire pixels
            "flinelen": "float",  # active fire front line length
        }
    elif layer == "fireline":
        dd = {
            "t": "datetime64[ns]",
        }
    elif layer == "newfirepix":
        dd = {
            "t": "datetime64[ns]",
        }

    return dd

def read_gdf_fire(event_id, layer='perimeter'):
    """
    Read one layer for of a output file for the given single fire event.

    Args:
        event_id (str): Fire event ID.
        layer (str, optional): Layer to which data needs to be read for the given fire. Defaults to 'perimeter'.

    Returns:
        GeoDataFrame: GeoDataFrame with data for the specified layer for the specified fire.
    """
    
    dd = set_dd(layer=layer)
    
    fnm = set_gdffile(event_id)
    if check_gdffile(event_id):
        gdf = gpd.read_file(fnm, layer=layer)
        # In addition to the GeoDataFrame that is read, add the diagnostic data entries
        for k, tp in dd.items():
            gdf[k] = gdf[k].astype(tp)        
    else:
        gdf = None
        print('File does not exist')
    return gdf

def read_1fire(Event_ID):
    """
    Read all three layers (fire perimeter, fire line, new fire pixels) of a output file for the given single fire event.

    Args:
        Event_ID (str): Fire event ID.

    Returns:
        tuple: Tuple of three GeoDataFrames (one each for fire perimeter, fire line, new fire pixels).
    """
    # Layer name is perimeter, but we are actually interested in the stored area data
    gdf_farea_rd = read_gdf_fire(Event_ID, layer='perimeter')
    gdf_fline_rd = read_gdf_fire(Event_ID, layer='fireline')
    gdf_nfp_rd = read_gdf_fire(Event_ID, layer='newfirepix')
    return gdf_farea_rd, gdf_fline_rd, gdf_nfp_rd

def get_gdf_times_and_time_range(gdf, conv_delta=None, start_time=None, end_time=None, num_hours=None, save_csv=None):
    """
    Get the set of times corresponding to the geometries in the given GeoDataFrame, and the time range for rasterization.
    The set of times is the set of unique times corresponding to the geometries in the given GeoDataFrame, with conv_delta applied if given.
    The caveat is that each of the times should be spaced 12 hours apart and if the geometry at that time is null, we still include
    it as a time step for rasterization (with null geometry). This is because FEDS data is recorded every 12 hours, but we want to rasterize 
    at a 1-hour time step, so we need to include all the 1-hour time steps even if there is no data for them.
    The time range is a list of times starting from start_time (or minimum time in the GeoDataFrame if start_time is not given) and ending 
    at end_time (or maximum time in the GeoDataFrame if end_time is not given), with a step of 1 hour. If num_hours is given, it is used to 
    calculate end_time based on start_time instead of using the maximum time in the GeoDataFrame.

    Args:
        gdf (GeoDataFrame): GeoDataFrame with geometries of fire data.
        conv_delta (pandas.Timedelta, optional): Time delta to apply to FEDS data times for conversion from LST to UTC. Defaults to None.
        start_time (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        end_time (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
        save_csv(str, optional): Path to save the CSV file with the times corresponding to fire duration and whether there are
         FEDS observations for that time. Do not save CSV if the given path is None. Defaults to None.

    Returns:
        tuple: Tuple of (set of times corresponding to the geometries in the given GeoDataFrame, time range for rasterization).
    """
    start_time = gdf['t'].min() if start_time is None else start_time
    # If num_hours is given, use that to generate times; otherwise, calculate num_hours based on end_time
    if num_hours is None or num_hours == 0:
        end_time = gdf['t'].max() if end_time is None else end_time
        num_hours = int( (end_time - start_time).total_seconds() / 3600)
    gdf_times = { t + conv_delta if conv_delta is not None else t for t in gdf.loc[gdf.geometry.notnull(), 't'].unique() }
    # Fill in times with 12-hr time intervals
    gdf_times_12hr_grid = pd.date_range(start=min(gdf_times), end=max(gdf_times), freq="12h")
    time_range = [start_time + timedelta(hours=i) for i in range(num_hours+1)]

    if save_csv is not None:
        df = pd.DataFrame(time_range, columns=['time'])
        df['feds'] = df['time'].isin(gdf_times_12hr_grid)
        df.to_csv(save_csv, index=False)

    return gdf_times_12hr_grid, time_range


def get_canonical_grid_from_bounds(bounds, resolution):
    """
    Build a canonical north-up raster grid from bounds and resolution.

    Args:
        bounds (tuple/list): Bounds in (minx, miny, maxx, maxy) order.
        resolution (float): Pixel resolution in meters.

    Returns:
        tuple: Tuple (transform, width, height).
    """
    minx, miny, maxx, maxy = [float(v) for v in bounds]
    width = int(np.round((maxx - minx) / float(resolution)))
    height = int(np.round((maxy - miny) / float(resolution)))

    if width <= 0 or height <= 0:
        raise ValueError(f'Invalid target grid for bounds {bounds} and resolution {resolution}')

    transform = from_origin(minx, maxy, float(resolution), float(resolution))
    return transform, width, height


def rasterize_gdf_and_save_as_tif(
    gdf, out_tif, gdf_times=None, time_range=None, crs='EPSG:5070', start_time=None, end_time=None,
    num_hours=None, use_prev=False, conv_delta=None, *, out_transform, out_width, out_height
):
    """
    Rasterize the given GeoDataFrame with FEDS data and save it as a TIF file at the specified output location. This process 
    may include projection to a new CRS (coordinate reference system) and/or resolution.

    Args:
        gdf (GeoDataFrame): GeoDataFrame with geometries of fire data.
        out_tif (str): Name of the output TIF file to where rasterized data should be stored.
        crs (str, optional): Coordinate reference system code (according to EPSG) in which the rasterized data
         should be stored. Defaults to 'EPSG:5070'.
        start_time (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        end_time (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
        use_prev (bool, optional): Whether to use previous non-null data for rasterization if no data exists. Defaults to False.
        conv_delta (pandas.Timedelta, optional): Time delta to apply to FEDS data times for conversion from LST to UTC. Defaults to None.
        out_transform (Affine): Target transform for writing to the canonical final grid.
        out_width (int): Target width for writing to the canonical final grid.
        out_height (int): Target height for writing to the canonical final grid.
    """
    # Convert gdf (GeoDataFrame) to correct CRS and get width, height based on desired resolution
    gdf = gdf.to_crs(crs)
    gdf['t_shift'] = gdf['t'] + conv_delta if conv_delta is not None else gdf['t']
    if gdf_times is None or time_range is None:
        gdf_times, time_range = get_gdf_times_and_time_range(
            gdf, conv_delta=conv_delta, start_time=start_time, end_time=end_time, num_hours=num_hours
        )

    transform = out_transform
    width = int(out_width)
    height = int(out_height)

    prev_non_null_df_time = None
    rasterized_bands = []

    # For each time, if there is FEDS data, use it for rasterization; otherwise, use previously used data (if use_prev is True)
    # This is needed because the time step in the time range is 1-hour, but FEDS data is recorded every 12 hours
    for t in time_range:
        if t in gdf_times:
            df = gdf[(gdf['t_shift'] == t) & (gdf.geometry.notna())]
            prev_non_null_df_time = t
        elif use_prev and prev_non_null_df_time is not None:
            df = gdf[gdf['t_shift'] == prev_non_null_df_time].dropna(subset=['geometry'])
        else:
            if use_prev:
                df = None
            else:
                continue

        # Rasterize the selected geometries; if no data, create an empty raster (with NaN values)
        if df is not None and not df.empty:
            shapes = ((geom, 1) for geom in df.geometry if geom is not None)
            raster = rasterize(
                shapes,
                out_shape=(height, width),
                transform=transform,
                fill=0,
                all_touched=True,
                dtype='uint8'
            )
        else:
            raster = np.full((height, width), np.nan, dtype='float32')
        
        rasterized_bands.append(raster)
    
    if len(rasterized_bands) == 0:
        raise ValueError('No raster bands were produced for the given time range')

    # Merge all rasters into a single array, and then a single TIF file
    stacked_array = np.stack(rasterized_bands)
    with rasterio.open(
        out_tif,
        'w',
        driver='GTiff',
        height=stacked_array.shape[1],
        width=stacked_array.shape[2],
        count=stacked_array.shape[0],  # number of bands
        dtype=stacked_array.dtype,
        crs=crs,
        transform=transform
    ) as dst:
        for i in range(stacked_array.shape[0]):
            dst.write(stacked_array[i], i + 1)  # rasterio bands are 1-based
            
def driver_feds(
    fid, final_bounds, res=300.0, fire_start=None, fire_end=None, num_hours=None, plot_orig=False, use_prev=False,
    conv_delta=None
):
    """
    Driver function for obtaining, cropping, resampling, and plotting FEDS data.

    Args:
        fid (str): Fire event ID.
        final_bounds (tuple): Tuple of the bounding box to which the output FEDS data should be padded. The tuple must
         contain four values, corresponding to the minimum x-coordinate, minimum y-coordinate, maximum x-coordinate,
         and maximum y-coordinate, in that order.
        res (float, optional): Resolution in which FEDS data should be rasterized. Defaults to 300.0.
        fire_start (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        fire_end (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
        plot_orig (bool, optional): True if FEDS data should be plotted; False otherwise. Defaults to False.
        use_prev (bool, optional): Whether to use previous non-null data for rasterization if no data exists. Defaults to False.
        conv_delta (pandas.Timedelta, optional): Time delta to apply to FEDS data times for conversion from LST to UTC. Defaults to None.
    Raises:
        ValueError: Occurs when fire area or new fire pixel data is empty.
    """
    gdf_farea_rd, gdf_fline_rd, gdf_nfp_rd = read_1fire(fid)
    gdfs = {
        "farea" : gdf_farea_rd,   # fire area
        "fline" : gdf_fline_rd,     # active fire line
        "nfp" : gdf_nfp_rd          # new fire pixels
    }

    times, trange = get_gdf_times_and_time_range(
        gdf_farea_rd, conv_delta=conv_delta, 
        start_time=fire_start, end_time=fire_end, num_hours=num_hours,
        save_csv=os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, f'fire_times.csv')
    )
    final_transform, final_width, final_height = get_canonical_grid_from_bounds(final_bounds, res)

    for var in gdfs:
        # For each "variable" (farea, fline, nfp), get the temporary data/video file names
        var_tif = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )

        # Some fires do not have fire line data (because they are too small or satellite observations are sporadic), so
        # just skip them. However, if fires do not have fire area or new fire pixel data, this is a problem, so
        # raise an error.
        if gdfs[var].empty:
            if var == 'fline':
                continue
            else:
                raise ValueError(f'GDF for variable {var} is empty')

        final_out_tif = gen_util.get_output_data_filename(fid, var, gen_util.subdir_firespread)

        rasterize_gdf_and_save_as_tif(
            gdfs[var],
            out_tif=var_tif,
            gdf_times=times,
            time_range=trange,
            use_prev=use_prev,
            conv_delta=conv_delta,
            out_transform=final_transform,
            out_width=final_width,
            out_height=final_height
        )
        shutil.copyfile(var_tif, final_out_tif)

        # Plot FEDS data if needed
        if plot_orig:
            gen_util.create_animation_plot_from_tif(
                var_tif, var_vid, start_time=fire_start
            )

def read_firepix_1fire(yr, fid):
    """
    Read the fire pixel data for a specific year and fire from the FEDS25 dataset.

    Args:
        yr (int): Year in which the desired fire occurred.
        fid (str): Fire event ID.

    Returns:
        DataFrame: DataFrame with fire pixel data corresponding to the given fire.
    """

    fnm = os.path.join(dir_firepix, 'Firepix_' + str(yr) + '.csv')
    df = pd.read_csv(fnm, index_col=0, parse_dates=['t'])
    
    # Get sub-DataFrame corresponding only to the given fire (should only be one row)
    df_fp = df[(df.Event_ID == fid)]

    return df_fp

def get_gdf_firepix_t(df_fp, t, out_crs="epsg:4326"):
    """
    Get the GeoDataFrame of fire pixel data corresponding to the given time step using the given DataFrame.

    Args:
        df_fp (DataFrame): DataFrame with fire pixel data corresponding to a single fire.
        t (pandas.Timestamp): Time at which fire pixel data should be extracted.
        out_crs (str, optional): Output coordinate reference system (CRS) EPSG code in which the GeoDataFrame should be. 
         Defaults to "epsg:4326".

    Returns:
        GeoDataFrame: GeoDataFrame of fire pixel data at the given time step using data from the given DataFrame.
    """
    # Extract FRP data for the current time step
    df_fp_t = df_fp[(df_fp.t == t)]
    # FRP density in MW/m^2
    df_fp_t['FRPden'] = df_fp_t['FRP'] / (df_fp_t['DT'] * df_fp_t['DS']) / 1e4

    # Convert to GeoDataFrame in CRS EPSG:4326
    gdf_fp_t = gpd.GeoDataFrame(df_fp_t, geometry=gpd.points_from_xy(df_fp_t.Lon, df_fp_t.Lat, crs="epsg:4326"))

    # Create shapes for rasterization in the output CRS
    gdf_fp_t = gdf_fp_t.to_crs(out_crs)

    return gdf_fp_t

def rasterize_frp_and_save_as_tif(
    gdf_farea_rd, df_fp, out_tif, crs='EPSG:5070', start_time=None, end_time=None, num_hours=None,
    use_prev=False, conv_delta=None, *, out_transform, out_width, out_height
):
    """
    Rasterize the given GeoDataFrame with FRP data and save it as a TIF file at the specified output location. This process 
    may include projection to a new CRS (coordinate reference system) and/or resolution.

    Args:
        gdf_farea_rd (GeoDataFrame): GeoDataFrame with geometries of fire area data.
        df_fp (DataFrame): DataFrame with fire pixel data corresponding to a single fire.
        out_tif (str): Name of the output TIF file to where rasterized data should be stored.
        crs (str, optional): Coordinate reference system code (according to EPSG) in which the rasterized data
         should be stored. Defaults to 'EPSG:5070'.
        start_time (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        end_time (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
        use_prev (bool, optional): Whether to use previous non-null data for rasterization if no data exists. Defaults to False.
        conv_delta (pandas.Timedelta, optional): Time delta to apply to FEDS data times for conversion from LST to UTC. Defaults to None.
        out_transform (Affine): Target transform for writing to the canonical final grid.
        out_width (int): Target width for writing to the canonical final grid.
        out_height (int): Target height for writing to the canonical final grid.
    """
    # Convert gdf_farea_rd (GeoDataFrame) to correct CRS and get width, height based on desired resolution
    gdf_farea_rd = gdf_farea_rd.to_crs(crs)
    transform = out_transform
    width = int(out_width)
    height = int(out_height)

    start_time = gdf_farea_rd['t'].min() if start_time is None else start_time
    # If num_hours is given, use that to generate times; otherwise, calculate num_hours based on end_time
    if num_hours is None or num_hours == 0:
        end_time = gdf_farea_rd['t'].max() if end_time is None else end_time
        num_hours = int( (end_time - start_time).total_seconds() / 3600)
    gdf_times = { t + conv_delta if conv_delta is not None else t for t in gdf_farea_rd.loc[gdf_farea_rd.geometry.notnull(), 't'].unique() }
    time_range = [start_time + timedelta(hours=i) for i in range(num_hours+1)]

    # Add conv_delta to the times in df_fp for consistency with gdf_farea_rd times if conv_delta is given
    df_fp_local = df_fp.copy()
    df_fp_local['t'] = df_fp_local['t'] + conv_delta if conv_delta is not None else df_fp_local['t']

    prev_non_null_df_time = None
    rasterized_bands = []

    # For each time, if there is FEDS/FRP data, use it for rasterization; otherwise, use previously used data (if use_prev is True)
    # This is needed because the time step in the time range is 1-hour, but FEDS/FRP data is recorded every 12 hours
    for t in time_range:
        if t in gdf_times:
            gdf_fp_t = get_gdf_firepix_t(df_fp_local, t, out_crs=crs)
            prev_non_null_df_time = t
        elif use_prev and prev_non_null_df_time is not None:
            gdf_fp_t = get_gdf_firepix_t(df_fp_local, prev_non_null_df_time, out_crs=crs)
        else:
            if use_prev:
                gdf_fp_t = None
            else:
                continue

        # Rasterize the selected geometries; if no data, create an empty raster
        if gdf_fp_t is not None and not gdf_fp_t.empty:
            shapes = ((geom, value) for geom, value in zip(gdf_fp_t.geometry, gdf_fp_t['FRPden']))
            raster = rasterize(
                shapes,
                out_shape=(height, width),
                transform=transform,
                fill=0,
                all_touched=True,
                dtype='float'
            )
        else:
            raster = np.zeros((height, width), dtype='float')

        rasterized_bands.append(raster)
    
    if len(rasterized_bands) == 0:
        raise ValueError('No FRP raster bands were produced for the given time range')

    # Merge all rasters into a single array, and then a single TIF file
    stacked_array = np.stack(rasterized_bands)
    with rasterio.open(
        out_tif,
        'w',
        driver='GTiff',
        height=stacked_array.shape[1],
        width=stacked_array.shape[2],
        count=stacked_array.shape[0],  # number of bands
        dtype=stacked_array.dtype,
        crs=crs,
        transform=transform
    ) as dst:
        for i in range(stacked_array.shape[0]):
            dst.write(stacked_array[i], i + 1)  # rasterio bands are 1-based


def driver_frp(
    fid, final_bounds, res=300.0, fire_start=None, fire_end=None, num_hours=None, plot_orig=False, use_prev=False,
    conv_delta=None
):
    """
    Driver function for obtaining, cropping, resampling, and plotting FRP data.

    Args:
        fid (str): Fire event ID.
        final_bounds (tuple): Tuple of the bounding box to which the output FRP data should be padded. The tuple must
         contain four values, corresponding to the minimum x-coordinate, minimum y-coordinate, maximum x-coordinate,
         and maximum y-coordinate, in that order.
        res (float, optional): Resolution in which FRP data should be rasterized. Defaults to 300.0.
        fire_start (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        fire_end (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
        plot_orig (bool, optional): True if FRP data should be plotted; False otherwise. Defaults to False.
        use_prev (bool, optional): Whether to use previous non-null data for rasterization if no data exists. Defaults to False.
        conv_delta (pandas.Timedelta, optional): Time delta to apply to FEDS data times for conversion from LST to UTC. Defaults to None.
    """
    # Read FEDS2.5 MTBS fire area data for the given fire
    gdf_farea_rd, _, _ = read_1fire(fid) 
    if gdf_farea_rd is None or gdf_farea_rd.empty:
        print(f"No fire area data for {fid}, cannot process FRP.")
        return
    # Read FEDS2.5 MTBS fire pixel data for the given specific fire (data at all time steps)
    df_fp = read_firepix_1fire(gdf_farea_rd.t.min().year, fid) 
    var = 'frp'

    # Get the name for the temporary FRP data file
    var_tif = gen_util.get_temp_data_video_filename(
        fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
    )

    final_out_tif = gen_util.get_output_data_filename(fid, var, gen_util.subdir_firespread)

    final_transform, final_width, final_height = get_canonical_grid_from_bounds(final_bounds, res)
    rasterize_frp_and_save_as_tif(
        gdf_farea_rd,
        df_fp,
        out_tif=var_tif,
        start_time=fire_start,
        end_time=fire_end,
        num_hours=num_hours,
        use_prev=use_prev,
        conv_delta=conv_delta,
        out_transform=final_transform,
        out_width=final_width,
        out_height=final_height
    )
    shutil.copyfile(var_tif, final_out_tif)

    # Plot FRP data if needed
    if plot_orig:
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )
        gen_util.create_animation_plot_from_tif(
            var_tif, var_vid, start_time=fire_start
        )
