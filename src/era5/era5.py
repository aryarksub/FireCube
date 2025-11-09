import cdsapi as cds
import numpy as np
import os
import pandas as pd
import rasterio
from rasterio.transform import from_origin
import xarray as xr

import util.general_util as gen_util
import util.processing_util as proc_util

def download_ERA5_reg_1m(y,m,days,bnds, varERA5 = ['2m_temperature'], tmpfile='ERA5.grib'):
    """
    Download the ERA5 data given the specified parameters. The downloaded data will be in GRIB format.

    Args:
        y (int): Year from which data should be downloaded.
        m (int): Month from which data should be downloaded (January = 1, February = 2, ..., December = 12).
        days (list): List of days from which data should be downloaded (1, 2, ..., 31).
        bnds (list): List of lat/lon coordinates specified in terms of West/South/East/North coordinates.
        varERA5 (list, optional): List of ERA5 layers to download. Full list of available layers can be found 
         here: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Parameterlistings.
         Defaults to ['2m_temperature'].
        tmpfile (str, optional): Name of output file to which data will be downloaded. Defaults to 'ERA5.grib'.
    """

    # Data will be retrieved from ERA5-Land Reanalysis dataset
    # Details: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview
    dataset = 'reanalysis-era5-land'
    
    # ERA5 area is defined as [N,W,S,E]; bnds is defined as [W,S,E,N]
    # Convert from bnds [lon0,lat0,lon1,lat1] to area [lat1,lon0,lat0,lon1]
    area = [bnds[i] for i in [3,0,1,2]]   
    
    # API request with variables (layers), days/month/year, times (all hours of the day), area, data/download format
    request = {
        'variable': varERA5,
        'year': str(y),
        'month': str(m).zfill(2),
        'day': [str(d) for d in days],
        'time': [
            '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'
        ],
        'area': area,
        'data_format': 'grib',
        'download_format': 'unarchived',
    }
    
    # Use cdsapi to fetch the ERA5 data for multiple days in a month and save to a temporary file.
    # Note that to use cdsapi, it must be locally set up. See this guide: https://cds.climate.copernicus.eu/how-to-api
    c = cds.Client()
    c.retrieve(dataset, request).download(target=tmpfile)
    
def download_ERA5_reg(fid, df_t, bnds, varERA5=['2m_temperature'], fnmERA5='./ERA5.nc', redo=False):
    """
    Download and extract the rasterized data for ERA5 layers. The final data will be stored in NetCDF format at the
    file specified by fnmERA5.

    Args:
        fid (str): Fire event ID
        df_t (DataFrame): Series of fire start and end times, with a 1-day buffer.
        bnds (list): List of coordinates in the order - min longitude (west), min latitude (south), max longitude (east),
         max latitude (north). All are represented in degrees.
        varERA5 (list, optional): List of ERA5 layers to download. Full list of available layers can be found 
         here: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Parameterlistings.
         Defaults to ['2m_temperature'].
        fnmERA5 (str, optional): Name of file to which final data will be stored. Defaults to './ERA5.nc'.
        redo (bool, optional): True if data should be re-downloaded even if the temporary NetCDF file exists; False otherwise. 
         Defaults to False.
    """
    # Stop process if we do not need to redo the download
    if not redo and os.path.exists(fnmERA5):
        return

    # Set the ERA5 bounds for the region (use 0.1 degree resolution and 1-pixel buffer at the boundary)
    ERA5bnds = proc_util.bufferbnds(bnds, res=0.1, bufgd=1)

    # Name all temporary files with "tmp_F" prefix so they can be easily found and deleted later
    temp_file_prefix = os.path.join(gen_util.dir_temp, f'tmp_F{fid}_')

    # For each unique month-year pair in the time series, download data if it does not already exist
    year_months = pd.to_datetime(df_t).dt.to_period('M')
    unique_months_sorted = sorted(year_months.unique())
    for ym in unique_months_sorted:
        temp_file = f'{temp_file_prefix}{ym.year}_{ym.month}.grib'
        if not os.path.exists(temp_file):
            days = list(np.unique(df_t[df_t.dt.month == ym.month].dt.day.values))
            download_ERA5_reg_1m(ym.year, ym.month, days, ERA5bnds, varERA5, temp_file)

    # Merge all GRIB files into one dataset
    ds = xr.open_mfdataset(temp_file_prefix + '*.grib', engine='cfgrib')
    # Clean the dataset: Only keep times within the desired period and remove data for duplicated times
    ds_clean = proc_util.clean_xr_dataset_by_times(ds, df_t.min(), df_t.max(), replace_nan=True)
    # Convert dataset to NetCDF
    ds_clean.to_netcdf(fnmERA5)

    # Remove all temporary files
    file_list = [f for f in os.listdir(gen_util.dir_temp) if f.startswith('tmp_F'+fid) and '.grib' in f]
    for file_path in file_list:
        os.remove(os.path.join(gen_util.dir_temp, file_path))

def get_data_vars_from_era5_dataset(ds):
    """
    Get all data variables/layers from the given ERA5 dataset. These are all the listed variables that are not also
    "coordinates" (examples of coordinates are "longitude", "latitude", "valid_time").

    Args:
        ds (xarray.Dataset): Dataset of ERA5 data.

    Returns:
        set: Set of all data variables.
    """
    return set(ds.variables.keys()) - set(ds.coords)

def convert_era5_nc_to_tif(ds, fid, data_variables):
    """
    Convert the given ERA5 dataset (taken from NetCDF file and stored in xarray Dataset) to TIF files.

    Args:
        ds (xarray.Dataset): Dataset of ERA5 data.
        fid (str): Fire event ID.
        data_variables (set): Set of ERA5 data variables/layers.
    """
    for var in data_variables:
        da = ds[var]
        lat, lon = ds['latitude'].values, ds['longitude'].values

        # Flip data if lat is ascending
        if lat[0] < lat[-1]:
            da = da.sel(latitude=lat[::-1])

        data = da.transpose("valid_time", "latitude", "longitude").values

        # Compute transform (assumes regular grid)
        res_lat = abs(lat[1] - lat[0])
        res_lon = abs(lon[1] - lon[0])
        transform = from_origin(west=lon[0], north=lat[0], xsize=res_lon, ysize=res_lat)

        # Write to multi-band GeoTIFF
        with rasterio.open(
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_data, 
                data_source=gen_util.subdir_era5, var_type=gen_util.subdir_type_original
            ),
            "w",
            driver="GTiff",
            height=data.shape[1],
            width=data.shape[2],
            count=data.shape[0],  # number of bands = time * step
            dtype=data.dtype,
            crs="EPSG:4326",
            transform=transform,
        ) as dst:
            for i in range(data.shape[0]):
                dst.write(data[i, :, :], i + 1)

def driver_era5(fid, vars, df_t, bounds, out_nc_file, plot_types=[]):
    """
    Driver function for downloading, cropping, resampling, and plotting ERA5 data.

    Args:
        fid (str): Fire event ID.
        vars (list): List of ERA5 variables to download. Full list of available layers can be found 
         here: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Parameterlistings.
        df_t (DataFrame): Series of fire start and end times, with a 1-day buffer.
        bounds (list): List of coordinates in the order - min longitude (west), min latitude (south), max longitude (east),
         max latitude (north). All are represented in degrees.
        out_nc_file (str): Name of file to which final data will be stored (should have .nc suffix for NetCDF format).
        plot_types (list, optional): List of types of data plots should be created for. Available types of data are listed
         in gen_util.var_types. Defaults to [].
    """

    # Download ERA5 data
    download_ERA5_reg(fid, df_t, bounds, vars, out_nc_file)
    ds = xr.open_dataset(gen_util.get_era5_nc_filename(fid), engine='netcdf4')
    data_vars = get_data_vars_from_era5_dataset(ds)
    # Convert data for TIF files
    convert_era5_nc_to_tif(ds, fid, data_vars)

    for var in data_vars:
        # Get data/plot file names
        era5_var_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_data,
                data_source=gen_util.subdir_era5, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]
        era5_plot_fnames = [
            gen_util.get_temp_data_video_filename(
                fid, var, dir_type=gen_util.dir_videos,
                data_source=gen_util.subdir_era5, var_type=vtype
            )
            for vtype in gen_util.var_types
        ]

        # Change CRS of original data to EPSG:5070
        proc_util.change_tif_crs(era5_var_fnames[0], era5_var_fnames[1], 'EPSG:5070')
        # Resample CRS-converted data to resolution of 9000m (9km)
        proc_util.resample_tif(era5_var_fnames[1], era5_var_fnames[2], 9000)

        # Plot the specified types of data (original/converted/resampled)
        for plot_type in set(plot_types):
            if plot_type in gen_util.var_types:
                index = gen_util.var_types.index(plot_type)
                gen_util.create_animation_plot_from_tif(
                    in_tif=era5_var_fnames[index],
                    out_file=era5_plot_fnames[index],
                    start_time=gen_util.get_fire_start_from_ds(ds),
                    mask=True,
                    ignore_small_neg=False
                )