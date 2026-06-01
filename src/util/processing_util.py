from datetime import timedelta
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.windows import from_bounds, Window
import xarray as xr
import multiprocessing as mpc
import os
import psutil

from util.validation_util import GLOBAL_NULL_VALUE

TOT_THREADS = mpc.cpu_count()
USE_THREADS = max(1, TOT_THREADS // 4)

os.environ["GDAL_NUM_THREADS"] = str(USE_THREADS)
# Get available system memory (bytes)
available_ram = psutil.virtual_memory().available
# Use 20% of available RAM
gdal_cache_mb = int(available_ram * 0.2 / 1024 / 1024)
# Clamp to reasonable bounds
gdal_cache_mb = max(1024, min(gdal_cache_mb, 16384))
os.environ["GDAL_CACHEMAX"] = str(gdal_cache_mb)


def bufferbnds(bnds, res=0.005, bufgd=1):
    """
    Create new bounds with the given resolution (in degrees) and grid buffer that enclose the original bounds.

    Args:
        bnds (list): Original bounds, listed in the order of minimum longitude, minimum latitude, maximum
         longitude, and maximum latitude.
        res (float, optional): Grid resolution in degrees. Defaults to 0.005.
        bufgd (int, optional): Buffer pixels (in terms of resolution). Defaults to 1.

    Returns:
        list: New bounds, listed in the order of minimum longitude, minimum latitude, maximum longitude, 
        and maximum latitude.
    """
    # Derive grid bounds
    lon0 = np.floor(bnds[0] / res - bufgd) * res 
    lat0 = np.floor(bnds[1] / res - bufgd) * res
    lon1 = np.ceil(bnds[2] / res + bufgd) * res 
    lat1 = np.ceil(bnds[3] / res + bufgd) * res 

    bufferedbnds = [lon0,lat0,lon1,lat1]
    return bufferedbnds

def add_time_buffers(df_t):
    """
    Add a one-day buffer to the beginning and end of the given time Series. Four new timestamps will be added: 24 hours
    before the beginning, 12 hours before the beginning, 12 hours after the end, 24 hours after the end.

    Args:
        df_t (pandas.Series): Series of timestamps

    Returns:
        pandas.Series: Series of timestamps with a one day buffer
    """
    # Add 1-day buffer at the beginning and the end
    pds1 = pd.Series([df_t.min()+timedelta(days=-1),df_t.min()+timedelta(hours=-12)], index=[-1,-2])
    pds2 = pd.Series([df_t.max()+timedelta(hours=12),df_t.max()+timedelta(days=1)], index=[-3,-4])
    return pds1._append(df_t)._append(pds2)

def clean_xr_dataset_by_times(ds, start_time, end_time, replace_nan=False):
    """
    Clean the given xarray Dataset by removing times outside the given start to end timeframe and dropping entries
    with duplicate timestamps.

    Args:
        ds (xarray.Dataset): Multi-layer dataset with data recorded at several time-steps.
        start_time (pandas.Timestamp): Timestamp of the start time of the valid timeframe
        end_time (pandas.Timestamp): Timestamp of the end time of the valid timeframe
        replace_nan (bool, optional): True if NaN values in the dataset should be replaced by a global null
         value (-9999); False otherwise. Defaults to False.

    Returns:
        xarray.Dataset: Cleaned dataset with unique data-time pairings within the given timeframe
    """
    time_vals = ds['time'].values
    step_vals = ds['step'].values
    valid_times = ds['valid_time'].values

    records = []
    for i in range(len(time_vals)):
        for j in range(len(step_vals)):
            vt = valid_times[i, j]
            # Filter to only valid_time within [start_time, end_time]
            if start_time <= vt <= end_time:
                is_valid = True
                for var in ds.data_vars:
                    val = ds[var].isel(time=i, step=j).values
                    # Even if there are NaN values, if we wish to replace them (replace_nan=True), then 
                    # treat the record as valid
                    if np.isnan(val).any() and not replace_nan:
                        is_valid = False
                        break
                if is_valid:
                    records.append({
                        'time_idx': i,
                        'step_idx': j,
                        'time': time_vals[i],
                        'step': step_vals[j],
                        'valid_time': vt
                    })

    df = pd.DataFrame(records)

    # Drop duplicates, keeping the first with data (i.e., earliest forecast without NaN data)
    df = df.sort_values(['valid_time', 'time'])  # ensures earliest forecast is kept
    df = df.drop_duplicates(subset='valid_time', keep='first')

    # Create an empty mask of shape (time, step)
    mask = np.full(ds['valid_time'].shape, False)

    # Set mask to True for selected (time_idx, step_idx) pairs
    for row in df.itertuples():
        mask[row.time_idx, row.step_idx] = True

    # Build a new dimension to hold filtered valid_times
    valid_time_vals = df['valid_time'].values

    # Prepare filtered dataset
    ds_filtered = xr.Dataset(coords={'valid_time': valid_time_vals,
                                    'latitude': ds.latitude,
                                    'longitude': ds.longitude})

    # Add each variable
    for var in ds.data_vars:
        data = ds[var].values  # shape: (time, step, lat, lon)
        # Get the data where mask is True, and reshape to (valid_time, lat, lon)
        filtered_data = data[mask]  # shape: (valid_time, lat, lon)
        if replace_nan:
            filtered_data = np.nan_to_num(filtered_data, nan=GLOBAL_NULL_VALUE)
        ds_filtered[var] = xr.DataArray(filtered_data, dims=('valid_time', 'latitude', 'longitude'),
                                    coords={'valid_time': valid_time_vals,
                                            'latitude': ds.latitude,
                                            'longitude': ds.longitude})

    return ds_filtered

def change_tif_crs(in_tif_file, out_tif_file, new_crs):
    """
    Change the coordinate reference system (CRS) of a given TIF file to a (new) given CRS.

    Args:
        in_tif_file (str): Path to the original TIF file.
        out_tif_file (str): Path to the new TIF file (with the new CRS).
        new_crs (str): EPSG code of the new CRS (e.g. "EPSG:5070")
    """
    with rasterio.open(in_tif_file) as src:
        transform, width, height = calculate_default_transform(
            src.crs, new_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': new_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(out_tif_file, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):  # assuming each time-band is a separate band
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=new_crs,
                    resampling=Resampling.nearest,  # or bilinear/cubic etc.
                    num_threads=USE_THREADS
                )

def resample_tif(in_tif, out_tif, target_res=None, catype=False):
    """
    Resample the given TIF to the given target resolution using either nearest-neighbor or bilinear resampling.
    If no target resolution is given, then the default resolution used is the closest multiple of 30 to the
    current TIF resolution.

    Args:
        in_tif (str): Path to the input TIF that needs to be resampled.
        out_tif (str): Path to the output TIF where resampled data will be written.
        target_res (float, optional): Target resolution for resampling. If None, then the target resolution
         will be the closest multiple of 30 to the original resolution. Defaults to None.
        catype (bool, optional): True if the resampling method should be nearest-neighbor; False if the resampling
         method should be bilinear. Defaults to False.

    Raises:
        ValueError: Occurs if the input TIF file has a non-projected CRS
    """
    # print(f"Resampling {in_tif} to {out_tif} with target resolution {target_res} and catype={catype}")
    resampling = Resampling.nearest if catype else Resampling.bilinear

    with rasterio.open(in_tif) as src:
        if not src.crs or not src.crs.is_projected:
            raise ValueError("CRS must be projected")

        # Original resolution
        xres = src.transform.a

        # Default target resolution
        if target_res is None:
            target_res = round(xres / 30) * 30

        # Compute new transform and dimensions
        transform, width, height = calculate_default_transform(
            src.crs,
            src.crs,
            src.width,
            src.height,
            *src.bounds,
            resolution=target_res
        )

        profile = src.profile.copy()
        profile.update(
            transform=transform,
            width=width,
            height=height,
            tiled=True,
            blockxsize=512,
            blockysize=512,
            interleave="band",
            compress="lzw",
        )

        with rasterio.open(out_tif, "w", **profile) as dst:
            chunk_size = 32  # tune this

            for b in range(1, src.count + 1, chunk_size):
                band_range = list(range(b, min(b + chunk_size, src.count + 1)))

                reproject(
                    source=rasterio.band(src, band_range),
                    destination=rasterio.band(dst, band_range),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=src.crs,
                    resampling=resampling
                )

def crop_tif_based_on_area(in_tif, out_tif, bounds):
    """
    Crop the given TIF to the specified bounds.

    Args:
        in_tif (str): Path of the input TIF file to be cropped.
        out_tif (str): Path where the output (cropped) TIF file will be saved.
        bounds (list): List of bounding box coordinates in the same CRS as the input TIF. There should be
         four coordinates: x_min, y_min, x_max, y_max.
    """
    with rasterio.open(in_tif) as src:
        # Create a window from the bounding box
        x1,y1,x2,y2 = bounds
        window = from_bounds(x1, y1, x2, y2, transform=src.transform)
        # Round window bounds up to make sure we get the entire fire area
        window = Window(window.col_off, window.row_off, int(round(window.width)), int(round(window.height)))

        # Read the data within the window
        data = src.read(window=window)

        # Compute the new transform for the cropped window
        new_transform = src.window_transform(window)

        # Define new metadata
        new_meta = src.meta.copy()
        new_meta.update({
            "height": data.shape[1],
            "width": data.shape[2],
            "transform": new_transform
        })

        with rasterio.open(out_tif, 'w', **new_meta) as dst:
            dst.write(data)

def pad_bounds_to_resolution_multiple(bounds, coarsest_res=9000):
    """
    Pad the given bounds to be a multiple of the given resolution (typically the coarsest resolution amongst a set
    of TIFs). If the required padding is less than 10% of the given resolution, add an additional full resolution buffer.

    Args:
        bounds (list): List of bounding box coordinates in the same CRS as the input TIF. There should be
         four coordinates: x_min, y_min, x_max, y_max.
        coarsest_res (float, optional): Coarsest resolution in meters. Defaults to 9000.

    Returns:
        list: Padded bounds.
    """
    x1,y1,x2,y2 = bounds
    center_x, center_y = (x1 + x2) / 2, (y1 + y2) / 2

    # Use coarsest tif to determine number of pixels required to cover bounds
    width_m = x2 - x1   # width in meters
    height_m = y2 - y1  # height in meters

    # Padding needed in each direction to get dimensions that are a multiple of coarsest_res
    pad_x = coarsest_res - (width_m % coarsest_res) if width_m % coarsest_res != 0 else 0
    pad_y = coarsest_res - (height_m % coarsest_res) if height_m % coarsest_res != 0 else 0

    # If the previously derived padding is smaller than 10% of the coarsest resolution, add a 100% coarsest 
    # resolution buffer in that direction.
    if pad_x < coarsest_res / 10:   
        pad_x += coarsest_res
    if pad_y < coarsest_res / 10:
        pad_y += coarsest_res

    # Dimensions for output in meters (multiple of coarsest_res)
    final_width_m = width_m + pad_x
    final_height_m = height_m + pad_y

    # Bounding box for output
    half_width = final_width_m / 2
    half_height = final_height_m / 2
    final_x1, final_x2 = center_x - half_width, center_x + half_width
    final_y1, final_y2 = center_y - half_height, center_y + half_height

    return [final_x1, final_y1, final_x2, final_y2]

def center_and_crop_tifs_to_same_area(in_tifs, out_tifs, bounds):
    """
    Centers and crops the given TIFs to the same spatial area (given by bounds).

    Args:
        in_tifs (list): List of paths to input TIF files.
        out_tifs (list): List of paths to output (cropped and centered) TIF files.
        bounds (list): List of bounding box coordinates in the same CRS as the input TIF. There should be
         four coordinates: x_min, y_min, x_max, y_max.
    """
    # Basic assertions
    assert len(in_tifs) == len(out_tifs) and len(bounds) == 4

    final_bounds = pad_bounds_to_resolution_multiple(bounds, 9000)

    for (in_tif, out_tif) in zip(in_tifs, out_tifs):
        crop_tif_based_on_area(in_tif, out_tif, final_bounds)

def get_tif_bounds(tif):
    """
    Get the bounding box for the given TIF.

    Args:
        tif (str): Path to the input TIF.

    Returns:
        np.array: NumPy array containing raster bounds [x_min, y_min, x_max, y_max]
    """
    with rasterio.open(tif) as src:
        bounds = src.bounds
        return np.array([bounds.left, bounds.bottom, bounds.right, bounds.top])
    
def convert_null_values_to_nan(in_tif, out_tif=None, null_values=[GLOBAL_NULL_VALUE]):
    """
    Convert all specified null values in the given input TIF file to np.nan and write the
    modified data to the given output TIF file. If no output TIF is given, write the
    modified data back to the input TIF.

    Args:
        in_tif (str): Path to input TIF.
        out_tif (str, optional): Path to output TIF. Defaults to None.
        null_values (list, optional): List of values to be treated as null. Defaults to [GLOBAL_NULL_VALUE].

    Returns:
        str: Path to output TIF.
    """
    # Open input
    with rasterio.open(in_tif) as src:
        data = src.read()  # Read as (bands, rows, cols)
        meta = src.meta.copy()

        # print(in_tif, data.shape)

        # Include raster's built-in nodata value if present
        if src.nodata is not None:
            null_values = list(null_values) + [src.nodata]

    # Convert to float so nan can exist
    data = data.astype(float)

    # Convert all specified null values to nan
    for nv in null_values:
        data[data == nv] = np.nan

    # Output filename
    if out_tif is None:
        out_tif = in_tif

    # Write output
    meta.update(dtype='float32', nodata=np.nan)

    with rasterio.open(out_tif, 'w', **meta) as dst:
        dst.write(data.astype("float32"))

    return out_tif
