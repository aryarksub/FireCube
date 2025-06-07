from io import BytesIO
import os
import rasterio
import requests
import time
import zipfile

import util.general_util as gen_util
import util.processing_util as proc_util

LF_BASE_URL = "https://lfps.usgs.gov"
LF_JOB_URL = f"{LF_BASE_URL}/api/job/submit"
LF_STATUS_URL = f"{LF_BASE_URL}/api/job/status"

def download_landfire_data(bounds, out_file='temp.zip', email="a@a.com", layers=['ASP2020'], out_proj='102003', resample_res=300, redo=False):
    if not redo and os.path.exists(out_file):
        return

    # set the LF bounds for the region (use 0.1 degree resolution and 1-pixel buffer at the boundary)
    LF_bounds = proc_util.bufferbnds(bounds, res=0.1, bufgd=1)

    job_payload = {
        "Email": email,
        "Layer_List": ';'.join(layers),
        "Area_of_Interest": " ".join([str(bnd) for bnd in LF_bounds]), # west, south, east, north
        "Output_Projection": out_proj,
        "Resample_Resolution": resample_res
    }

    submit_response = requests.post(LF_JOB_URL, json=job_payload)
    submit_data = submit_response.json()
    job_id = submit_data.get("jobId")
    status_payload = {
        "JobId": job_id
    }

    status = "Executing"
    while status == "Executing":
        time.sleep(10)  # wait 10 seconds between polls
        status_response = requests.post(LF_STATUS_URL, json=status_payload)
        status_data = status_response.json()
        status = status_data.get("status")

    if status == "Succeeded":
        download_url = status_data['outputFile']
        if download_url:
            result_response = requests.get(download_url)
            with open(out_file, "wb") as f:
                f.write(result_response.content)
    else:
        print('Data download failed')

def split_tifs_in_zip(zip_path, fid):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        file_list = zip_ref.namelist()
        tif_files = [file for file in file_list if file.endswith('.tif')]
        assert len(tif_files) == 1
        
        with zip_ref.open(tif_files[0]) as tif_file:
            with rasterio.open(BytesIO(tif_file.read())) as src:
                profile = src.profile.copy()
                for i in range(1, src.count+1):
                    band_data = src.read(i)
                    band_name = src.descriptions[i-1] if src.descriptions[i-1] else f'band_{i}'
                    var_name = band_name.replace('US_', '').lower()

                    band_profile = profile.copy()
                    band_profile.update({
                        "count": 1,
                        "dtype": band_data.dtype,
                        "driver": "GTiff"
                    })

                    output_path = gen_util.get_temp_data_video_filename(
                        fid, var_name, dir_type=gen_util.dir_data,
                        data_source=gen_util.subdir_lf, var_type=gen_util.subdir_type_original
                    )
                    with rasterio.open(output_path, "w", **band_profile) as dst:
                        dst.write(band_data, 1)
                        dst.set_band_description(1, band_name)

def get_lf_layers_given_vars_and_year(var_names, year):
    year_cutpoints = [2014,2020,2022,2023,2024]
    version_lists = {
        'EVT': [str(num) + 'EVT' for num in [105,200,220,230,240]],
        'FBFM13': ['105FBFM13', '200F13_20', '220F13_22', '230FBFM13', '240FBFM13'],
        'FBFM40': ['105FBFM40', '200F40_20', '220F40_22', '230FBFM40', '240FBFM40'],
    }
    latest_versions = {
        'EVT': '250EVT',
        'FBFM13': '250FBFM13',
        'FBFM40': '250FBFM40'
    }
    layers = []

    for var in var_names:
        if var in ['ASP', 'ELEV', 'SLPD']:
            layer_name = f'{var}2020'
        elif var in ['EVT', 'FBFM13', 'FBFM40']:
            versions = version_lists[var]
            for ind in range(len(year_cutpoints)):
                if year <= year_cutpoints[ind]:
                    layer_name = versions[ind]
                    break
            else:
                layer_name = latest_versions[var]
        else:
            print(f'Layer for {var} not in existing cases; will not be included in data download')
        layers.append(layer_name)
    return layers


def driver_landfire(fid, var_names, bounds, fire_start, plot_types=[]):
    # bounds should be in W,S,E,N format
    layers = get_lf_layers_given_vars_and_year(var_names, int(fire_start.year))
    lf_zip_path = gen_util.get_lf_zip_filename(fid)

    download_landfire_data(
        bounds=bounds, out_file=lf_zip_path,
        layers=layers, out_proj="102003", resample_res=300, redo=False
    )
    split_tifs_in_zip(lf_zip_path, fid)

    data_vars = gen_util.get_tif_vars_in_dir(
        os.path.join(gen_util.dir_temp, gen_util.dir_data, fid, gen_util.subdir_lf, gen_util.subdir_type_original)
    )

    for var in data_vars:
        lf_var_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_data,
                data_source=gen_util.subdir_lf, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]
        lf_plot_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_videos,
                data_source=gen_util.subdir_lf, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]

        # Change CRS of original data to EPSG:5070
        proc_util.change_tif_crs(lf_var_fnames[0], lf_var_fnames[1], 'EPSG:5070')
        # Resample CRS-converted data to resolution defined in resample_tif (closest multiple of 30)
        proc_util.resample_tif(lf_var_fnames[1], lf_var_fnames[2], target_res=None)

        for plot_type in set(plot_types):
            if plot_type in gen_util.var_types:
                index = gen_util.var_types.index(plot_type)
                gen_util.create_animation_plot_from_tif(
                    in_tif=lf_var_fnames[index],
                    out_file=lf_plot_fnames[index],
                    start_time=fire_start,
                    mask=True,
                    ignore_small_neg=True
                )



