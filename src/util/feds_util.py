from datetime import timedelta
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_bounds

import util.general_util as gen_util
import util.processing_util as proc_util

# Supress SettingWithCopyWarning
pd.set_option('mode.chained_assignment', None)

dir_feds25 = os.path.join('inputData', 'Full_FEDS')
dir_firepix = os.path.join('inputData', 'firepix')
feds_firelist = os.path.join('inputData', 'fireslist2012-2024.csv')

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
    gdf_fperim_rd = read_gdf_fire(Event_ID, layer='perimeter')
    gdf_fline_rd = read_gdf_fire(Event_ID, layer='fireline')
    gdf_nfp_rd = read_gdf_fire(Event_ID, layer='newfirepix')
    return gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd

def rasterize_gdf_and_save_as_tif(gdf, out_tif, resolution, crs='EPSG:5070', start_time=None, end_time=None, num_hours=None):
    """
    Rasterize the given GeoDataFrame with FEDS data and save it as a TIF file at the specified output location. This process 
    may include projection to a new CRS (coordinate reference system) and/or resolution.

    Args:
        gdf (GeoDataFrame): GeoDataFrame with geometries of fire data.
        out_tif (str): Name of the output TIF file to where rasterized data should be stored.
        resolution (float): Output resolution in which the rasterized data should be stored.
        crs (str, optional): Coordinate reference system code (according to EPSG) in which the rasterized data
         should be stored. Defaults to 'EPSG:5070'.
        start_time (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        end_time (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
    """
    # Convert gdf (GeoDataFrame) to correct CRS and get width, height based on desired resolution
    gdf = gdf.to_crs(crs)
    minx, miny, maxx, maxy = gdf.total_bounds
    width = int((maxx - minx) / resolution)
    height = int((maxy - miny) / resolution)
    transform = from_bounds(minx, miny, maxx, maxy, width, height)

    start_time = gdf['t'].min() if start_time is None else start_time
    # If num_hours is given, use that to generate times; otherwise, calculate num_hours based on end_time
    if num_hours is None or num_hours == 0:
        end_time = gdf['t'].max() if end_time is None else end_time
        num_hours = int( (end_time - start_time).total_seconds() / 3600)
    gdf_times = set(gdf.loc[gdf.geometry.notnull(), 't'].unique())
    time_range = [start_time + timedelta(hours=i) for i in range(num_hours+1)]

    prev_non_null_df_time = None
    rasterized_bands = []

    # For each time, if there is FEDS data, use it for rasterization; otherwise, use previously used data (if it exists)
    # This is needed because the time step in the time range is 1-hour, but FEDS data is recorded every 12 hours
    for t in time_range:
        if t in gdf_times:
            df = gdf[gdf['t'] == t].dropna(subset=['geometry'])
            prev_non_null_df_time = t
        elif prev_non_null_df_time is not None:
            df = gdf[gdf['t'] == prev_non_null_df_time].dropna(subset=['geometry'])
        else:
            df = None

        # Rasterize the selected geometries; if no data, create an empty raster
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
            raster = np.zeros((height, width), dtype='uint8')
        
        rasterized_bands.append(raster)
    
    # Merge all rasters into a single array, and then a single TIF file
    stacked_array = np.stack(rasterized_bands)
    temp_tif_file = 'temp_rasterized_gdf.tif'
    with rasterio.open(
        temp_tif_file,
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

    # Resample the tif to correct resolution in case the above procedure slightly shifts resolution
    proc_util.resample_tif(temp_tif_file, out_tif, target_res=resolution)
    os.remove(temp_tif_file)
            
def driver_feds(fid, final_bounds, res=300.0, fire_start=None, fire_end=None, num_hours=None, plot_orig=False):
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
    """
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = read_1fire(fid)
    gdfs = {
        "fperim" : gdf_fperim_rd,   # fire perimeter
        "fline" : gdf_fline_rd,     # active fire line
        "nfp" : gdf_nfp_rd          # new fire pixels
    }

    for var in gdfs:
        # For each "variable" (fperim, fline, nfp), get the temporary data/video file names
        var_tif = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )

        # Save rasterized data to the temporary TIF file
        rasterize_gdf_and_save_as_tif(
            gdfs[var], out_tif=var_tif, resolution=res, start_time=fire_start, end_time=fire_end, num_hours=num_hours
        )

        out_batch = gen_util.get_out_batch_for_tif(var_tif)
        final_out_tif = gen_util.get_output_data_filename(fid, var, out_batch)

        # Add padding to the temporary TIF file based on the final_bounds argument
        proc_util.pad_tif_to_bounds(var_tif, final_out_tif, final_bounds)

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

def rasterize_frp_and_save_as_tif(gdf_fperim_rd, df_fp, out_tif, resolution, crs='EPSG:5070', start_time=None, end_time=None, num_hours=None):
    """
    Rasterize the given GeoDataFrame with FRP data and save it as a TIF file at the specified output location. This process 
    may include projection to a new CRS (coordinate reference system) and/or resolution.

    Args:
        gdf_fperim_rd (GeoDataFrame): GeoDataFrame with geometries of fire perimeter data.
        df_fp (DataFrame): DataFrame with fire pixel data corresponding to a single fire.
        out_tif (str): Name of the output TIF file to where rasterized data should be stored.
        resolution (float): Output resolution in which the rasterized data should be stored.
        crs (str, optional): Coordinate reference system code (according to EPSG) in which the rasterized data
         should be stored. Defaults to 'EPSG:5070'.
        start_time (pandas.Timestamp, optional): Start time from which data should be rasterized. Defaults to None.
        end_time (pandas.Timestamp, optional): End time until which data should be rasterized. Defaults to None.
        num_hours (int, optional): Number of hours for which data should be rasterized. Defaults to None.
    """
    # Convert gdf_fperim_rd (GeoDataFrame) to correct CRS and get width, height based on desired resolution
    gdf_fperim_rd = gdf_fperim_rd.to_crs(crs)
    minx, miny, maxx, maxy = gdf_fperim_rd.total_bounds
    width = int((maxx - minx) / resolution)
    height = int((maxy - miny) / resolution)
    transform = from_bounds(minx, miny, maxx, maxy, width, height)

    start_time = gdf_fperim_rd['t'].min() if start_time is None else start_time
    # If num_hours is given, use that to generate times; otherwise, calculate num_hours based on end_time
    if num_hours is None or num_hours == 0:
        end_time = gdf_fperim_rd['t'].max() if end_time is None else end_time
        num_hours = int( (end_time - start_time).total_seconds() / 3600)
    gdf_times = set(gdf_fperim_rd.loc[gdf_fperim_rd.geometry.notnull(), 't'].unique())
    time_range = [start_time + timedelta(hours=i) for i in range(num_hours+1)]

    prev_non_null_df_time = None
    rasterized_bands = []

    # For each time, if there is FEDS/FRP data, use it for rasterization; otherwise, use previously used data (if it exists)
    # This is needed because the time step in the time range is 1-hour, but FEDS/FRP data is recorded every 12 hours
    for t in time_range:
        if t in gdf_times:
            gdf_fp_t = get_gdf_firepix_t(df_fp, t, out_crs=crs)
            prev_non_null_df_time = t
        elif prev_non_null_df_time is not None:
            gdf_fp_t = get_gdf_firepix_t(df_fp, prev_non_null_df_time, out_crs=crs)
        else:
            gdf_fp_t = None

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
    
    # Merge all rasters into a single array, and then a single TIF file
    stacked_array = np.stack(rasterized_bands)
    temp_tif_file = 'temp_rasterized_gdf.tif'
    with rasterio.open(
        temp_tif_file,
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

    # Resample the tif to correct resolution in case the above procedure slightly shifts resolution
    proc_util.resample_tif(temp_tif_file, out_tif, target_res=resolution)
    os.remove(temp_tif_file)

def driver_frp(fid, final_bounds, res=300.0, fire_start=None, fire_end=None, num_hours=None, plot_orig=False):
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
    """
    # Read FEDS2.5 MTBS fire perimeter data for the given fire
    gdf_fperim_rd, _, _ = read_1fire(fid) 
    if gdf_fperim_rd is None or gdf_fperim_rd.empty:
        print(f"No fire perimeter data for {fid}, cannot process FRP.")
        return
    # Read FEDS2.5 MTBS fire pixel data for the given specific fire (data at all time steps)
    df_fp = read_firepix_1fire(gdf_fperim_rd.t.min().year, fid) 
    var = 'frp'

    # Get the name for the temporary FRP data file
    var_tif = gen_util.get_temp_data_video_filename(
        fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_frp, var_type=gen_util.subdir_type_resample
    )

    # Save rasterized data to the temporary TIF file
    rasterize_frp_and_save_as_tif(
        gdf_fperim_rd, df_fp, out_tif=var_tif, resolution=res, start_time=fire_start, end_time=fire_end, num_hours=num_hours
    )

    out_batch = gen_util.get_out_batch_for_tif(var_tif)
    final_out_tif = gen_util.get_output_data_filename(fid, var, out_batch)

    # Add padding to the temporary TIF file based on the final_bounds argument
    proc_util.pad_tif_to_bounds(var_tif, final_out_tif, final_bounds)

    # Plot FRP data if needed
    if plot_orig:
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_frp, var_type=gen_util.subdir_type_resample
        )
        gen_util.create_animation_plot_from_tif(
            var_tif, var_vid, start_time=fire_start
        )