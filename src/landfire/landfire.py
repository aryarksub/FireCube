from io import BytesIO
import os
import rasterio
import requests
import time
import zipfile
from datetime import datetime

import util.general_util as gen_util
import util.processing_util as proc_util

LF_BASE_URL = "https://lfps.usgs.gov"
LF_JOB_URL = f"{LF_BASE_URL}/api/job/submit"
LF_STATUS_URL = f"{LF_BASE_URL}/api/job/status"
LF_HEALTH_URL = f"{LF_BASE_URL}/api/healthCheck"

def download_landfire_data(bounds, out_file='temp.zip', email="a@a.com", layers=['ASP2020'], out_proj='102003', resample_res=30, redo=False):
    """
    Download the LANDFIRE data given the specified parameters. The downloaded data will be stored as a ZIP file.

    Args:
        bounds (list): List of lat/lon coordinates specified in terms of West/South/East/North coordinates.
        out_file (str, optional): Name of output file to which data will be downloaded. Defaults to 'temp.zip'.
        email (str, optional): Email to which data will be sent. Defaults to "a@a.com".
        layers (list, optional): List of LANDFIRE layers to download. Full list of available layers can be found at
         https://www.landfire.gov/sites/default/files/documents/LF_Data_Dictionary.pdf in Section 4. Defaults to ['ASP2020'].
        out_proj (str, optional): Output projection (ESRI code). Defaults to '102003'.
        resample_res (int, optional): Resolution to which data should be resampled (in meters). Defaults to 30.
        redo (bool, optional): True if data should be re-downloaded even if the temporary ZIP file exists; False otherwise. 
         Defaults to False.

    Returns:
        bool: True if data download succeeded; False otherwise.
    """
    if not redo and os.path.exists(out_file):
        return True
    
    # If LANDFIRE server is not healthy, download will fail (indicate as such and end process)
    health_response = requests.get(LF_HEALTH_URL, timeout=120)
    health_data = health_response.json()
    healthy = health_data['success']
    if not healthy:
        print('LANDFIRE API not healthy. Retry later')
        return False

    # Set the LF bounds for the region (use 0.1 degree resolution and 2-pixel buffer at the boundary)
    # 2-pixel buffer ensures we have at least 9km buffer even at high latitude parts of CONUS
    LF_bounds = proc_util.bufferbnds(bounds, res=0.1, bufgd=2)

    job_payload = {
        "Email": email,
        "Layer_List": ';'.join(layers),
        "Area_of_Interest": " ".join([str(bnd) for bnd in LF_bounds]), # west, south, east, north
        "Output_Projection": out_proj
    }

    # If resample resolution is in the range [31,9999], include it in the API request
    if isinstance(resample_res, int):
        if (resample_res >= 31) & (resample_res <= 9999):
            job_payload["Resample_Resolution"] = resample_res

    # Submit request to get data and retrieve corresponding job ID
    submit_response = requests.post(LF_JOB_URL, json=job_payload, timeout=120)
    submit_data = submit_response.json()
    job_id = submit_data.get("jobId")
    status_payload = {
        "JobId": job_id
    }

    # While the job is executing/pending, check status every 10 seconds
    status = "Executing"
    while status == "Executing" or status == "Pending":
        time.sleep(10) 
        status_response = requests.post(LF_STATUS_URL, json=status_payload, timeout=120)
        status_data = status_response.json()
        status = status_data.get("status")

    # If job succeeds, download the data from the URL specified in job response
    if status == "Succeeded":
        download_url = status_data['outputFile']
        if download_url:
            result_response = requests.get(download_url, timeout=600)
            with open(out_file, "wb") as f:
                f.write(result_response.content)
            return True
        return False
    else:
        print('Data download failed')
        return False

def split_tifs_in_zip(zip_path, fid):
    """
    Split the ZIP file with all LANDFIRE data into separate TIF files for each layer.

    Args:
        zip_path (str): Path corresponding to ZIP file with all LANDFIRE data.
        fid (str): Fire event ID.
    """
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # There should be one TIF file in the ZIP file with all data
        file_list = zip_ref.namelist()
        tif_files = [file for file in file_list if file.endswith('.tif')]
        assert len(tif_files) == 1
        
        with zip_ref.open(tif_files[0]) as tif_file:
            with rasterio.open(BytesIO(tif_file.read())) as src:
                profile = src.profile.copy()
                # Each band corresponds to a single LANDFIRE layer
                for i in range(1, src.count+1):
                    band_data = src.read(i)
                    # Get band name based on description and normalize it
                    band_name = src.descriptions[i-1] if src.descriptions[i-1] else f'band_{i}'
                    var_name = band_name.replace('US_', '').lower()

                    # Create profile specific to the band (along with generic profile)
                    band_profile = profile.copy()
                    band_profile.update({
                        "count": 1,
                        "dtype": band_data.dtype,
                        "driver": "GTiff"
                    })

                    # Get name of TIF file data should be written to, and write data there
                    output_path = gen_util.get_temp_data_video_filename(
                        fid, var_name, dir_type=gen_util.dir_data,
                        data_source=gen_util.subdir_lf, var_type=gen_util.subdir_type_original
                    )
                    with rasterio.open(output_path, "w", **band_profile) as dst:
                        dst.write(band_data, 1)
                        dst.set_band_description(1, band_name)

def get_lf_layers_given_vars_and_start_time(var_names, start_time, state):
    """
    Get LANDFIRE layer name for the given variables based on the fire year and state.

    Args:
        var_names (list): List of LANDFIRE variables.
        year (int): Year to use when determining version.
        state (str): Abbreviation of the state for which LANDFIRE layer names need to be obtained.

    Returns:
        list: List of LANDFIRE layer names (one for each given variable).
    """
    # Yearly cutpoints (one for each time LANDFIRE released a new version, based on most recent LANDFIRE reconfiguration)
    cutoffs = [
        (datetime(2016, 12, 31), 2016),
        (datetime(2022, 12, 31), 2022),
        (datetime(2023, 10, 1), 2023),
        (datetime(2024, 10, 1), 2024),
    ]
    # Map of variables EVT/FBFM13/FBFM40 to layer names based on each yearly cutpoint
    version_lists = {
        'EVT': [str(num) + 'EVT' for num in [200,230,240]],
        'FBFM13': ['200F13_20', '230FBFM13', '240FBFM13'],
        'FBFM40': ['200F40_20', '230FBFM40', '240FBFM40'],
    }
    latest_versions = {
        'EVT': '250EVT',
        'FBFM13': '250FBFM13',
        'FBFM40': '250FBFM40'
    }
    layers = []

    for var in var_names:
        # ASP/ELEV/SLPD only have one layer name
        if var in ['ASP', 'ELEV']:
            layer_name = f'LF2020_{var.upper()[0]+var[1:].lower()}'
        elif var in ['SLPD']:
            layer_name = f'LF2020_SlpD'
        # EVT/FBFM13/FBFM40 have layer names as defined in maps above
        elif var in ['EVT', 'FBFM13', 'FBFM40']:
            versions = version_lists[var]
            layer_name = None
            for cutoff, year in cutoffs:
                if start_time >= cutoff:
                    layer_name = f'LF{year}_{var}' #versions[ind]
            if layer_name is None:
                #NOTE: pass error here if we want to strictly enforce the landfire version to precede the fire
                layer_name = f'LF{2016}_{var}'
        # ROADS only has layer name 240ROADS_23 for CONUS
        # Currently, the supposed layer for Alaska/Hawaii is 220ROADS_20, but this doesn't work, so skip for now
        elif var in ['ROADS']:
            if state in {'AK', 'HI'}:
                continue
            else:
                layer_name = 'LF2023_Roads'
        # Other layer names are not supported by our code
        else:
            print(f'Layer for {var} not in existing cases; will not be included in data download')
        layers.append(layer_name)
    return layers


def driver_landfire(fid, var_names, bounds, fire_start, plot_types=[]):
    """
    Driver function for downloading, cropping, resampling, and plotting LANDFIRE data.

    Args:
        fid (str): Fire event ID.
        var_names (list): List of LANDFIRE layers to download. Full list of available layers can be found at
         https://www.landfire.gov/sites/default/files/documents/LF_Data_Dictionary.pdf in Section 4.
        bounds (list): List of lat/lon coordinates specified in terms of West/South/East/North coordinates.
        fire_start (pandas.TimeStamp): Time of fire start (including buffer).
        plot_types (list, optional): List of types of data plots should be created for. Available types of data are listed
         in gen_util.var_types. Defaults to [].

    Raises:
        RuntimeError: Occurs when LANDFIRE download failed.
    """
    # bounds should be in W,S,E,N format
    layers = get_lf_layers_given_vars_and_start_time(var_names, fire_start, fid[:2])
    lf_zip_path = gen_util.get_lf_zip_filename(fid)

    download_success = download_landfire_data(
        bounds=bounds, out_file=lf_zip_path,
        layers=layers, out_proj="102003", resample_res=30, redo=False
    )
    if not download_success:
        raise RuntimeError('LANDFIRE data download failed')
    
    split_tifs_in_zip(lf_zip_path, fid)

    # Get all LANDFIRE data variables
    data_vars = gen_util.get_tif_vars_in_dir(
        os.path.join(gen_util.dir_temp, gen_util.dir_data, fid, gen_util.subdir_lf, gen_util.subdir_type_original)
    )

    for var in data_vars:
        # Get data/plot file names
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
        proc_util.resample_tif(lf_var_fnames[1], lf_var_fnames[2], target_res=None, catype=True)

        # Plot the specified types of data (original/converted/resampled)
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
