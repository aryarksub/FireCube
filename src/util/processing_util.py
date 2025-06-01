from datetime import timedelta
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine
from rasterio.warp import calculate_default_transform, reproject, Resampling
import xarray as xr

def bufferbnds(bnds,res=0.005,bufgd=1):
    ''' create new bounds with res that enclose the original bounds

    Parameters
    ----------
    bnds : list
        [lon0,lat0,lon1,lat1] of the original bounds
    res : float
        grid resolution
    bufgd : int 
        the buffer (in number of grids with res)

    Returns
    -------
    bufferedbnds : list
        [lon0,lat0,lon1,lat1] of the new bounds
    '''

    # derive grid bounds
    lat0 = np.floor(bnds[1]/res-bufgd)*res 
    lon0 = np.floor(bnds[0]/res-bufgd)*res 
    lat1 = np.ceil(bnds[3]/res+bufgd)*res 
    lon1 = np.ceil(bnds[2]/res+bufgd)*res

    bufferedbnds = [lon0,lat0,lon1,lat1]
    return bufferedbnds

def add_time_buffers(df_t):
    # add 1-day buffer at the beginning and the end
    pds1 = pd.Series([df_t.min()+timedelta(days=-1),df_t.min()+timedelta(hours=-12)], index=[-1,-2])
    pds2 = pd.Series([df_t.max()+timedelta(hours=12),df_t.max()+timedelta(days=1)], index=[-3,-4])
    return pds1._append(df_t)._append(pds2)

def clean_xr_dataset_by_times(ds, start_time, end_time):
    time_vals = ds['time'].values
    step_vals = ds['step'].values
    valid_times = ds['valid_time'].values

    records = [
        {'time_idx': i, 'step_idx': j,
        'time': time_vals[i], 'step': step_vals[j], 'valid_time': valid_times[i, j]}
        for i in range(len(time_vals)) for j in range(len(step_vals))
    ]

    df = pd.DataFrame(records)

    # Filter to only valid_time within [start_time, end_time]
    df = df[(df['valid_time'] >= start_time) & (df['valid_time'] <= end_time)]

    # Drop duplicates, keeping the first (i.e., earliest forecast)
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
        ds_filtered[var] = xr.DataArray(filtered_data, dims=('valid_time', 'latitude', 'longitude'),
                                    coords={'valid_time': valid_time_vals,
                                            'latitude': ds.latitude,
                                            'longitude': ds.longitude})

    return ds_filtered

def change_tif_crs(in_tif_file, out_tif_file, new_crs):
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
                    resampling=Resampling.nearest  # or bilinear/cubic etc.
                )

def resample_tif(in_tif, out_tif, target_res=None):
    with rasterio.open(in_tif) as src:
        transform = src.transform
        crs = src.crs
        width = src.width
        height = src.height

        # Current pixel size
        pixel_width = transform.a
        pixel_height = -transform.e

        # Ensure CRS is in meters
        if not crs.is_projected:
            raise ValueError("CRS must be projected (e.g., UTM) to use meters.")
        
        if target_res is None:
            target_res = round(pixel_width / 30) * 30

        # Compute new width/height to avoid padding (crop if needed)
        new_width = int((width * pixel_width) // target_res)
        new_height = int((height * pixel_height) // target_res)

        # New transform — same origin, new resolution
        new_transform = Affine(
            target_res, 0, transform.c,
            0, -target_res, transform.f
        )

        # Write to new file
        profile = src.profile
        profile.update({
            "height": new_height,
            "width": new_width,
            "transform": new_transform
        })

        with rasterio.open(out_tif, "w", **profile) as dst:
            for i in range(1, src.count+1):
                band_data = src.read(
                    i,
                    out_shape=(new_height, new_width),
                    resampling=Resampling.bilinear
                )
                dst.write(band_data, i)