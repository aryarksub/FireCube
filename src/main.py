import pandas as pd

from era5 import driver_era5
from pyregence.pyregence import driver_pyregence
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util

firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)

def process_single_fire(fid, era5_vars=[], do_pyr=True, verbose=False, plot=False):
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    bnds = proc_util.bufferbnds(gdf_fperim_rd.total_bounds, res=0.005, bufgd=1)
    df_t = pd.to_datetime(gdf_fperim_rd.t)
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    fire_row = firelist[firelist['Event_ID'] == fid]
    center_lat = round((fire_row['lat0'].values[0] + fire_row['lat1'].values[0]) / 2, 2)
    center_lon = round((fire_row['lon0'].values[0] + fire_row['lon1'].values[0]) / 2, 2)
    fire_start = pd.Timestamp(df_t_with_buffer.min())
    fire_hours = int( (df_t_with_buffer.max() - df_t_with_buffer.min()).total_seconds() / 3600) + 12

    if len(era5_vars) != 0:
        if verbose: print(f'Getting ERA5 data for fire {fid}')
        
        driver_era5(
            fid, era5_vars, df_t_with_buffer, bnds, 
            gen_util.get_era5_nc_filename(fid), 
            plot_types=gen_util.var_types if plot else []
        )
    else:
        if verbose: print(f'No ERA5 variables specified; not getting ERA5 data for fire {fid}')

    if do_pyr:
        if verbose: print(f'Getting Pyregence data for fire {fid}')
        driver_pyregence(fid, (center_lat, center_lon), fire_start, fire_hours)
    else:
        if verbose: print(f'Skipping Pyregence data for fire {fid}')

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'

    era5_vars = [] #['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']

    fid_to_use = zogg_id
    gen_util.create_dirs_for_fire(fid_to_use)
    process_single_fire(fid_to_use, era5_vars, do_pyr=True, verbose=True, plot=True)