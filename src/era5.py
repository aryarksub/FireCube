import cdsapi as cds
import numpy as np
import os
import pandas as pd
import xarray as xr

import util.general_util as gen_util
import util.processing_util as proc_util

def download_ERA5_reg_1m(y,m,days,bnds, varERA5 = ['2m_temperature'], tmpfile='ERA5.grib'):
    '''  use cdsapi to fetch the ERA5 data for multiple days in a month and save to a temporary file
    for each day, the 24 hourly data are retrieved.
    '''
    dataset = 'reanalysis-era5-land'
    # ERA5 area is defined as [N,W,S,E]; bnds is defined as [W,S,E,N]
    area = [bnds[i] for i in [3,0,1,2]]   # convert from bnds [lon0,lat0,lon1,lat1] to area [lat1,lon0,lat0,lon1]
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
    c = cds.Client()
    c.retrieve(dataset, request).download(target=tmpfile)
    
def download_ERA5_reg(fid, df_t, bnds, varERA5=['2m_temperature'], fnmERA5='./ERA5.nc', redo=False):
    ''' extract the ERA5 raster for all days of a fire. The output file `fnmERA5` is in netcdf format
    
    Parameters
    ----------
    df_t : pandas.DataFrame
        time series of fire start and end times
    bnds : list
        [latmin,latmax,lonmin,lonmax]
    varERA5 : str
        variable name in ERA5 data
    fnmERA5 : str
        name of netcdf file to save ERA5 data to
    
    Returns
    -------
    None    
    '''

    if not redo and os.path.exists(fnmERA5):
        return

    # set the ERA5 bounds for the region (use 0.1 degree resolution and 1-pixel buffer at the boundary)
    ERA5bnds = proc_util.bufferbnds(bnds, res=0.1, bufgd=1)

    temp_file_prefix = os.path.join(gen_util.dir_temp, f'tmp_F{fid}_')

    year_months = pd.to_datetime(df_t).dt.to_period('M')
    unique_months_sorted = sorted(year_months.unique())
    for ym in unique_months_sorted:
        temp_file = f'{temp_file_prefix}{ym.year}_{ym.month}.grib'
        if not os.path.exists(temp_file):
            days = list(np.unique(df_t[df_t.dt.month == ym.month].dt.day.values))
            download_ERA5_reg_1m(ym.year, ym.month, days, ERA5bnds, varERA5, temp_file)

    # merge all files into one netcdf file in the pre-set directory
    ds = xr.open_mfdataset(temp_file_prefix + '*.grib', engine='cfgrib')
    ds_clean = proc_util.clean_xr_dataset_by_times(ds, np.datetime64('2020-09-27T00:00:00'), np.datetime64('2020-10-03T23:00:00'))
    ds_clean.to_netcdf(fnmERA5)

    # remove all temporary files
    file_list = [f for f in os.listdir(gen_util.dir_temp) if f.startswith('tmp_F'+fid) and '.grib' in f]
    for file_path in file_list:
        os.remove(os.path.join(gen_util.dir_temp, file_path))

def driver_era5(fid, vars, df_t, bounds, out_file):
    download_ERA5_reg(fid, df_t, bounds, vars, out_file)