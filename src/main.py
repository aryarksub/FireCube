import os
import pandas as pd

from era5 import driver_era5
from landfire import driver_landfire
from pyregence.pyregence import driver_pyregence
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util

firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)

def process_single_fire(fid, era5_vars=[], do_pyr=True, lf_vars=[], verbose=False, plot=[], multi_plot=False):
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    bnds = proc_util.bufferbnds(gdf_fperim_rd.total_bounds, res=0.005, bufgd=1) # W,S,E,N
    df_t = pd.to_datetime(gdf_fperim_rd.t)
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    fire_row = firelist[firelist['Event_ID'] == fid]
    center_lat = round((fire_row['lat0'].values[0] + fire_row['lat1'].values[0]) / 2, 2)
    center_lon = round((fire_row['lon0'].values[0] + fire_row['lon1'].values[0]) / 2, 2)
    fire_start = pd.Timestamp(df_t_with_buffer.min())
    # Add 11 to go from mid-day to end of last day (hour 23:00 is the last time)
    fire_hours = int( (df_t_with_buffer.max() - df_t_with_buffer.min()).total_seconds() / 3600) + 11

    if len(era5_vars) != 0:
        if verbose: print(f'Getting ERA5 data for fire {fid}')
        
        driver_era5(
            fid, era5_vars, df_t_with_buffer, bnds, 
            gen_util.get_era5_nc_filename(fid), 
            plot_types=gen_util.var_types if gen_util.subdir_era5 in plot else []
        )
    else:
        if verbose: print(f'No ERA5 variables specified; not getting ERA5 data for fire {fid}')

    if do_pyr:
        if verbose: print(f'Getting Pyregence data for fire {fid}')
        driver_pyregence(
            fid, (center_lat, center_lon), fire_start, fire_hours, 
            plot_types=gen_util.var_types if gen_util.subdir_pyr in plot else []
        )
    else:
        if verbose: print(f'Skipping Pyregence data for fire {fid}')

    if len(lf_vars) != 0:
        if verbose: print(f'Getting LANDFIRE data for fire {fid}')

        driver_landfire(
            fid, lf_vars, bnds, fire_start, 
            plot_types=gen_util.var_types if gen_util.subdir_lf in plot else []
        )
    else:
        if verbose: print(f'No LANDFIRE variables specified; not getting LANDFIRE data for fire {fid}')

    # Convert gdf_perim to EPSG:5070
    gdf_fperim_5070 = gdf_fperim_rd.to_crs('EPSG:5070') 
    bounds_5070 = gdf_fperim_5070.total_bounds       

    # Crop era5/pyr/lf tifs to just surround the fire perim
    all_variable_input_tifs, all_variable_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(fid)
    proc_util.center_and_crop_tifs_to_same_area(all_variable_input_tifs, all_variable_output_tifs, bounds_5070)

    if multi_plot:
        for batch in gen_util.data_batches:
            if batch != gen_util.subdir_vis and batch != gen_util.subdir_firespread: # fire spread rasters (FEDS) not supported yet
                if verbose:
                    print(f'Generating plot for batch {batch} - fire {fid}')
                gen_util.create_multi_animation_for_dir(
                    os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, batch),
                    gen_util.get_output_data_filename(fid, batch, gen_util.subdir_vis),
                    fire_start
                )
        
        if verbose:
            print(f'Generating plot for all variables - fire {fid}')
        gen_util.create_multi_animation_from_tifs(
            all_variable_output_tifs,
            os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, gen_util.subdir_vis, 'all.mp4'),
            start_time=fire_start
        )

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'

    era5_vars = []#['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']
    get_pyr_data = False
    lf_vars = []#['ASP', 'ELEV', 'SLPD', 'EVT', 'FBFM13', 'FBFM40']
    plot_sources = []

    fid_to_use = zogg_id
    gen_util.create_dirs_for_fire(fid_to_use)
    process_single_fire(
        fid_to_use, era5_vars, do_pyr=get_pyr_data, lf_vars=lf_vars, verbose=True, plot=plot_sources, multi_plot=False
    )