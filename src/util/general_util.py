from datetime import datetime, timedelta
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
import rasterio
from tqdm import tqdm

dir_temp = 'temp'
dir_data = 'data'
dir_videos = 'videos'
subdir_era5 = 'era5'
subdir_pyr = 'pyr'
subdir_lf = 'lf'
subdir_type_original = 'original'
subdir_type_converted = 'converted'
subdir_type_resample = 'resample'

dir_output = 'output'
dir_sum_vis = 'sum_vis' # summary visualizations
dir_cubes = 'cubes' # top-level dir for cube outputs
subdir_vis = 'vis' # visualizations
subdir_lrc = 'low_res_climate' # era5 9000m
subdir_hrc = 'high_res_climate' # pyregence 600m
subdir_fuel_topo = 'fuel_topo' # pyregence 30m
subdir_landfire = 'landfire' # landfire 30m
subdir_firespread = 'fire_spread' # FEDS rasters 30m

data_sources = [subdir_era5, subdir_pyr, subdir_lf]
var_types = [subdir_type_original, subdir_type_converted, subdir_type_resample]

data_batches = [subdir_vis, subdir_lrc, subdir_hrc, subdir_fuel_topo, subdir_landfire, subdir_firespread]

def create_dirs_for_fire(fid):
    # create folders in temporary (local) directory
    for top_dir in [dir_data, dir_videos]:
        top_path = os.path.join(dir_temp, top_dir)
        if not os.path.exists(top_path):
            os.makedirs(top_path)
        
        for data_src_dir in data_sources:
            full_path = os.path.join(top_path, fid, data_src_dir)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            
            for data_type_dir in var_types:
                full_path2 = os.path.join(full_path, data_type_dir)
                if not os.path.exists(full_path2):
                    os.makedirs(full_path2)

    # create folders in output directory
    for top_dir in [dir_sum_vis, dir_cubes]:
        path = os.path.join(dir_output, top_dir)
        if not os.path.exists(path):
            os.makedirs(path)
    for batch in data_batches:
        batch_path = os.path.join(dir_output, dir_cubes, fid, batch)
        if not os.path.exists(batch_path):
            os.makedirs(batch_path)

def get_era5_nc_filename(fid):
    return os.path.join(dir_temp, dir_data, fid, subdir_era5, f'{fid}_era5_original.nc')

def get_pyr_tar_filename(fid):
    return os.path.join(dir_temp, dir_data, fid, subdir_pyr, f'{fid}_pyr_original.tar')

def get_lf_zip_filename(fid):
    return os.path.join(dir_temp, dir_data, fid, subdir_lf, f'{fid}_lf_original.zip')

def get_temp_data_video_filename(fid, var, dir_type='data', data_source='era5', var_type='original'):
    assert dir_type in [dir_data, dir_videos]
    assert data_source in data_sources
    assert var_type in var_types
    dir = os.path.join(dir_temp, dir_type, fid, data_source, var_type)
    filename = f"{var}.{'tif' if dir_type == dir_data else 'mp4'}"
    return os.path.join(dir, filename)

def get_out_batch_for_tif(tif):
    with rasterio.open(tif) as src:
        res_x, res_y = src.transform.a, -src.transform.e
        assert res_x == res_y
    
    if res_x == 9000 and subdir_era5 in tif:
        return subdir_lrc
    elif res_x == 600 and subdir_pyr in tif:
        return subdir_hrc
    elif res_x == 30 and subdir_pyr in tif:
        return subdir_fuel_topo
    elif subdir_lf in tif:
        return subdir_landfire
    else:
        print(f'Resolution + Var combination for file {tif} is not supported')

def get_output_data_filename(fid, var, batch_dir):
    assert batch_dir in data_batches
    dir = os.path.join(dir_output, dir_cubes, fid, batch_dir)
    filename = f"{var}.{'mp4' if batch_dir == subdir_vis else 'tif'}"
    return os.path.join(dir, filename)

def to_datetime(date):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return datetime.utcfromtimestamp(timestamp)

def get_fire_start_from_ds(ds):
    return to_datetime(ds['valid_time'].values[0])

## Code to get the extent of a raster
def get_extent(profile):
    transform = profile['transform']
    width, height = profile['width'], profile['height']
    left, top = transform * (0, 0)
    right, bottom = transform * (width, height)
    return (left, right, bottom, top)

def create_animation_plot_from_tif(in_tif, out_file, start_time, mask=False, ignore_small_neg=False):
    with rasterio.open(in_tif) as src:
        data = src.read()
        transform = src.transform
        crs = src.crs
        extent = get_extent(src.profile)

    # If mask, then do not include zero and NaN values in min/max computation
    # If ignore_small_negs, then do not include small values (<-500) in min/max computation
    condition = np.ones_like(data, dtype=bool)
    if mask:
        condition &= ~np.isnan(data)
        condition &= data != 0
    if ignore_small_neg:
        condition &= data >= -500
    masked_data = data[condition]

    vmin = np.min(masked_data) if len(masked_data) != 0 else 0
    vmax = np.max(masked_data) if len(masked_data) != 0 else 1

    base  = start_time
    times = [base + timedelta(hours=i) for i in range(data.shape[0])]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))
    img = ax.imshow(data[0], cmap='viridis', vmin=vmin, vmax=vmax, extent=extent, origin='upper')
    title = ax.set_title(f"Band 0 - Time {times[0]}")

    pbar = tqdm(total=data.shape[0], desc='Animating', unit='frames', leave=False)
    def update(frame):
        pbar.update(1)
        pbar.set_postfix_str(f"Frame {frame+1}/{data.shape[0]}")
        img.set_array(data[frame])
        title.set_text(f"Band {frame} - Times {times[frame]}")
        return img, title
    
    ani = animation.FuncAnimation(
        fig, update, frames=data.shape[0], interval=300, blit=False
    )
    ani.save(out_file, writer='ffmpeg', fps=10)
    plt.close()

def get_tif_vars_in_dir(dir):
    if type(dir) is str:
        dir = Path(dir)
    return [f.stem for f in dir.glob('*.tif')]