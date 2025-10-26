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
feds_firelist = os.path.join('inputData', 'fireslist2012-2023_MTBSfull.csv')

def set_gdffile( Event_ID):
    """
    Constructs the file path for a GeoPackage file based on the given year and Event ID.

    Args:
        yr (int): The year associated with the file.
        Event_ID (str): The unique identifier for the event.

    Returns:
        str: The constructed file path for the GeoPackage file.
    """
    # for all years under FEDS folder, try to find file corresponding to given event ID
    target_file = f'{Event_ID}.gpkg'

    for root, dirs, files in os.walk(dir_feds25):
        if target_file in files:
            return os.path.join(root, target_file)
    return None

    # fnm = os.path.join(dir_feds25, Event_ID + '.gpkg')
    # return fnm

def check_gdffile( Event_ID):
    """
    Checks if a GDF file exists for a given year and event ID.

    Args:
        yr (int): The year of the event.
        Event_ID (str): The unique identifier for the event.

    Returns:
        bool: True if the GDF file exists, False otherwise.
    """
    fnm = set_gdffile( Event_ID)
    return os.path.exists(os.path.expanduser(fnm))

def set_dd(layer="perimeter"):

    # diagnostic data name and types saved in geojson files (in addition to geometries)
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

def read_gdf_fire(Event_ID,layer='perimeter'):
    """ read 1 layer of a output file for a single fire event
    """
    dd = set_dd(layer=layer)
    
    fnm = set_gdffile( Event_ID)
    if check_gdffile( Event_ID):
        gdf = gpd.read_file(fnm, layer=layer)
        for k, tp in dd.items():
            gdf[k] = gdf[k].astype(tp)        
    else:
        gdf = None
        print('File does not exist')
    return gdf

def read_1fire(Event_ID):
    """ read all three layers of a output file for a single fire event
    """
    gdf_fperim_rd = read_gdf_fire(Event_ID, layer='perimeter')
    gdf_fline_rd = read_gdf_fire(Event_ID, layer='fireline')
    gdf_nfp_rd = read_gdf_fire(Event_ID, layer='newfirepix')
    return gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd

def rasterize_gdf_and_save_as_tif(gdf, out_tif, resolution, crs='EPSG:5070', start_time=None, end_time=None, num_hours=None):
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
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = read_1fire(fid)
    gdfs = {
        "fperim" : gdf_fperim_rd,   # fire perimeter
        "fline" : gdf_fline_rd,     # active fire line
        "nfp" : gdf_nfp_rd          # new fire pixels
    }

    for var in gdfs:
        var_tif = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_feds, var_type=gen_util.subdir_type_resample
        )

        rasterize_gdf_and_save_as_tif(
            gdfs[var], out_tif=var_tif, resolution=res, start_time=fire_start, end_time=fire_end, num_hours=num_hours
        )

        out_batch = gen_util.get_out_batch_for_tif(var_tif)
        final_out_tif = gen_util.get_output_data_filename(fid, var, out_batch)

        proc_util.pad_tif_to_bounds(var_tif, final_out_tif, final_bounds)

        if plot_orig:
            gen_util.create_animation_plot_from_tif(
                var_tif, var_vid, start_time=fire_start
            )

def read_firepix_1fire(yr, fid):
    """ Reads the fire pixel data for a specific year from the FEDS25 dataset. """

    fnm = os.path.join(dir_firepix, 'Firepix_' + str(yr) + '.csv')
    df = pd.read_csv(fnm, index_col=0, parse_dates=['t'])
    
    df_fp = df[(df.Event_ID == fid)]

    return df_fp

def get_gdf_firepix_t(df_fp, t, out_crs="epsg:4326"):
        # extract FRP data for the current time step
        df_fp_t = df_fp[(df_fp.t == t)]
        df_fp_t['FRPden'] = df_fp_t['FRP'] / (df_fp_t['DT'] * df_fp_t['DS']) / 1e4  # FRP density in MW/m^2

        # convert to GeoDataFrame (epsg:4326)
        gdf_fp_t = gpd.GeoDataFrame(df_fp_t, geometry=gpd.points_from_xy(df_fp_t.Lon, df_fp_t.Lat, crs="epsg:4326"))

        # create shapes for rasterization
        gdf_fp_t = gdf_fp_t.to_crs(out_crs)

        return gdf_fp_t

def rasterize_frp_and_save_as_tif(gdf_fperim_rd, df_fp, out_tif, resolution, crs='EPSG:5070', start_time=None, end_time=None, num_hours=None):
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

    # For each time, if there is FEDS data, use it for rasterization; otherwise, use previously used data (if it exists)
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
            # shapes = ((geom, 1) for geom in df.geometry if geom is not None)
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
    # read FEDS25MTBS
    gdf_fperim_rd, _, _ = read_1fire(fid) # read FEDS25MTBS fire perimeter for a specific fire
    if gdf_fperim_rd is None or gdf_fperim_rd.empty:
        print(f"No fire perimeter data for {fid}, cannot process FRP.")
        return
    df_fp = read_firepix_1fire(gdf_fperim_rd.t.min().year, fid) # read FEDS25MTBS firepix data for a specific fire (all time steps)
    var = 'frp'

    var_tif = gen_util.get_temp_data_video_filename(
        fid, var, dir_type=gen_util.dir_data, data_source=gen_util.subdir_frp, var_type=gen_util.subdir_type_resample
    )

    rasterize_frp_and_save_as_tif(
        gdf_fperim_rd, df_fp, out_tif=var_tif, resolution=res, start_time=fire_start, end_time=fire_end, num_hours=num_hours
    )

    out_batch = gen_util.get_out_batch_for_tif(var_tif)
    final_out_tif = gen_util.get_output_data_filename(fid, var, out_batch)

    proc_util.pad_tif_to_bounds(var_tif, final_out_tif, final_bounds)

    if plot_orig:
        var_vid = gen_util.get_temp_data_video_filename(
            fid, var, dir_type=gen_util.dir_videos, data_source=gen_util.subdir_frp, var_type=gen_util.subdir_type_resample
        )
        gen_util.create_animation_plot_from_tif(
            var_tif, var_vid, start_time=fire_start
        )