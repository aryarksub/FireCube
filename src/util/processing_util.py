from datetime import timedelta
import numpy as np
import pandas as pd
import rasterio
from rasterio.coords import BoundingBox
from rasterio.transform import Affine, from_origin
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.windows import from_bounds, Window
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

    records = []
    for i in range(len(time_vals)):
        for j in range(len(step_vals)):
            vt = valid_times[i, j]
            # Filter to only valid_time within [start_time, end_time]
            if start_time <= vt <= end_time:
                is_valid = True
                for var in ds.data_vars:
                    val = ds[var].isel(time=i, step=j).values
                    if np.isnan(val).any():
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

def crop_tif_based_on_area(in_tif, out_tif, bounds):
    with rasterio.open(in_tif) as src:
        # Create a window from the bounding box
        x1,y1,x2,y2 = bounds
        window = from_bounds(x1, y1, x2, y2, transform=src.transform)
        # Round window bounds up to make sure we get the entire fire area
        window = Window(window.col_off, window.row_off, int(np.ceil(window.width)), int(np.ceil(window.height)))

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

def get_coarsest_resolution(tifs):
    max_res = 0
    for tif in tifs:
        with rasterio.open(tif) as src:
            res_x, res_y = src.res
            assert res_x == res_y
            max_res = max(max_res, res_x)
    return max_res

def center_and_crop_tifs_to_same_area(in_tifs, out_tifs, bounds):
    assert len(in_tifs) == len(out_tifs) and len(bounds) == 4

    x1,y1,x2,y2 = bounds
    center_x, center_y = (x1 + x2) / 2, (y1 + y2) / 2
    coarsest_res = get_coarsest_resolution(in_tifs)

    # Use coarsest tif to determine number of pixels required to cover bounds
    width_m = x2 - x1   # width in meters
    height_m = y2 - y1  # height in meters

    # Padding needed in each direction to get dimensions that are a multiple of coarsest_res
    pad_x = coarsest_res - (width_m % coarsest_res) if width_m % coarsest_res != 0 else 0
    pad_y = coarsest_res - (height_m % coarsest_res) if height_m % coarsest_res != 0 else 0

    # Dimensions for output in meters (multiple of coarsest_res)
    final_width_m = width_m + pad_x
    final_height_m = height_m + pad_y

    # Bounding box for output
    half_width = final_width_m / 2
    half_height = final_height_m / 2
    final_x1, final_x2 = center_x - half_width, center_x + half_width
    final_y1, final_y2 = center_y - half_height, center_y + half_height

    for (in_tif, out_tif) in zip(in_tifs, out_tifs):
        crop_tif_based_on_area(in_tif, out_tif, (final_x1, final_y1, final_x2, final_y2))

def pad_tif_to_bounds(in_tif, out_tif, bounds):
    with rasterio.open(in_tif) as src:
        original_data = src.read()
        original_bounds = src.bounds
        dtype = src.dtypes[0]
        count = src.count

        xres, yres = src.res
        assert xres == yres

        xmin, ymin, xmax, ymax = bounds
        bounds = BoundingBox(left=xmin, bottom=ymin, right=xmax, top=ymax)

        # Compute new width and height in pixels
        new_width = int(np.ceil((bounds.right - bounds.left) / xres))
        new_height = int(np.ceil((bounds.top - bounds.bottom) / yres))

        # Compute new transform for the new bounding box
        new_transform = from_origin(bounds.left, bounds.top, xres, yres)

        # Create empty array filled with zeros
        padded_data = np.zeros((count, new_height, new_width), dtype=dtype)

        # Compute pixel offsets to place original data
        col_offset = int((original_bounds.left - bounds.left) / xres)
        row_offset = int((bounds.top - original_bounds.top) / yres)

        # Copy original data into the padded array
        padded_data[:, row_offset:row_offset + src.height, col_offset:col_offset + src.width] = original_data

        # Update metadata
        profile = src.profile.copy()
        profile.update({
            'height': new_height,
            'width': new_width,
            'transform': new_transform
        })

        with rasterio.open(out_tif, 'w', **profile) as dst:
            dst.write(padded_data)

def get_tif_bounds(tif):
    with rasterio.open(tif) as src:
        bounds = src.bounds
        return np.array([bounds.left, bounds.bottom, bounds.right, bounds.top])