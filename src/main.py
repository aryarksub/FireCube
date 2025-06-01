import os
import pandas as pd

from era5 import driver_era5
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util

def process_single_fire(fid, era5_vars=[], verbose=False):
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    bnds = proc_util.bufferbnds(gdf_fperim_rd.total_bounds, res=0.005, bufgd=1)
    df_t = pd.to_datetime(gdf_fperim_rd.t)
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    if len(era5_vars) != 0:
        if verbose: print(f'Getting ERA5 data for fire {fid}')
        
        driver_era5(fid, era5_vars, df_t_with_buffer, bnds, gen_util.get_era5_nc_filename(fid))
    else:
        if verbose: print(f'No ERA5 variables specified; not getting ERA5 data for fire {fid}')
    

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'

    era5_vars = ['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']

    fid_to_use = zogg_id
    gen_util.create_dirs_for_fire(fid_to_use)
    process_single_fire(fid_to_use, era5_vars, verbose=True)