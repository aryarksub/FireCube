import os

dir_temp = 'temp'
dir_data = 'data'
dir_videos = 'videos'
subdir_era5 = 'era5'
subdir_pyr = 'pyr'
subdir_type_original = 'original'
subdir_type_converted = 'converted'
subdir_type_resample = 'resample'

def create_dirs_for_fire(fid):
    for top_dir in [dir_data, dir_videos]:
        top_path = os.path.join(dir_temp, top_dir)
        if not os.path.exists(top_path):
            os.makedirs(top_path)
        
        for data_src_dir in [subdir_era5, subdir_pyr]:
            full_path = os.path.join(top_path, fid, data_src_dir)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            
            for data_type_dir in [subdir_type_original, subdir_type_converted, subdir_type_resample]:
                full_path2 = os.path.join(full_path, data_type_dir)
                if not os.path.exists(full_path2):
                    os.makedirs(full_path2)

def get_era5_nc_filename(fid):
    return os.path.join(dir_temp, dir_data, fid, subdir_era5, f'{fid}_era5_original.nc')

def get_temp_data_video_filename(fid, var, dir_type='data', data_source='era5', var_type='original'):
    assert dir_type in [dir_data, dir_videos]
    assert data_source in [subdir_era5, subdir_pyr]
    assert var_type in [subdir_type_original, subdir_type_converted, subdir_type_resample]
    dir = os.path.join(dir_temp, dir_type, fid, data_source, var_type)
    filename = f"{var}.{'tif' if dir_type == dir_data else '.mp4'}"
    return os.path.join(dir, filename)