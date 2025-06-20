import os
import pandas as pd

from era5 import driver_era5
from landfire import driver_landfire
from pyregence.pyregence import driver_pyregence
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util

firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)

def process_single_fire(fid, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, plot=[], batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False):
    if verbose:
        print(f'Processing fire {fid}')
    
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    if gdf_fperim_rd is None:
        if verbose:
            print(f'FEDS file for fire {fid} does not exist - no processing will be done')
        return
    
    gen_util.create_dirs_for_fire(fid)
    bnds = proc_util.bufferbnds(gdf_fperim_rd.total_bounds, res=0.005, bufgd=1) # W,S,E,N
    df_t = pd.to_datetime(gdf_fperim_rd.t)
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    # Convert gdf_perim to EPSG:5070
    gdf_fperim_5070 = gdf_fperim_rd.to_crs('EPSG:5070') 
    bounds_5070 = gdf_fperim_5070.total_bounds

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

    # Crop era5/pyr/lf tifs to just surround the fire perim (ignore FEDS data for cropping since it needs to be padded)
    non_feds_input_tifs, non_feds_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(
        fid, exclude=[gen_util.subdir_feds]
    )
    # If there are no output tifs, then stop processing
    if len(non_feds_output_tifs) == 0:
        if verbose:
            print(f'No output TIFs were created - stopping processing for fire {fid}')
        return

    proc_util.center_and_crop_tifs_to_same_area(non_feds_input_tifs, non_feds_output_tifs, bounds_5070)

    # Bounding box for all variable/layer tifs is the same, so we can just take the box for the first tif
    largest_var_tif_bounds = proc_util.get_tif_bounds(non_feds_output_tifs[0])

    if do_feds:
        if verbose: print(f'Getting/Rasterizing FEDS data for fire {fid}')
        feds_util.driver_feds(
            fid, largest_var_tif_bounds, res=300.0, fire_start=fire_start, num_hours=fire_hours, plot_orig=False
        )
    else:
        if verbose: print(f'Skipping FEDS data for fire {fid}')

    all_variable_input_tifs, all_variable_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(fid)

    if batch_plot:
        for batch in gen_util.data_batches:
            if batch != gen_util.subdir_vis:
                if verbose:
                    print(f'Generating plot for batch {batch} - fire {fid}')
                gen_util.create_multi_animation_for_dir(
                    os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, batch),
                    gen_util.get_output_data_filename(fid, batch, gen_util.subdir_vis),
                    fire_start
                )

    if all_plot: 
        if verbose:
            print(f'Generating plot for all variables - fire {fid}')
        gen_util.create_multi_animation_from_tifs(
            all_variable_output_tifs,
            os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, gen_util.subdir_vis, 'all.mp4'),
            start_time=fire_start
        )

    gen_util.remove_temp_dir_files(
        fid, 
        del_dir_types=gen_util.dir_types,
        del_data_sources=del_sources,
        del_var_types=gen_util.var_types,
        remove_intermediate=del_intermediate,
        verbose=verbose
    )

def process_multiple_fires(fid_list=[], fid_file=None, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, plot=[], batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False):
    if len(fid_list) == 0 and fid_file is None:
        if verbose:
            print('No FIDs given - no processing will be done')
        return
    
    if len(fid_list) > 0:
        if verbose:
            print('Processing fires given in list argument')
        for fid in fid_list:
            process_single_fire(fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate)
    else:
        if verbose:
            print('Processing fires given in file-path argument')
        try:
            with open(fid_file, 'r', encoding='utf-8') as file:
                for line in file:
                    fid = line.strip()
                    process_single_fire(fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate)
        except:
            if verbose:
                print(f'Error when reading file {fid_file} - no processing will be done')

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'

    era5_vars = ['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']
    get_pyr_data = True
    lf_vars = ['ASP', 'ELEV', 'SLPD', 'EVT', 'FBFM13', 'FBFM40']
    rasterize_feds = True
    plot_sources = []

    fids_to_use = [zogg_id]
    process_multiple_fires(
        fid_list=fids_to_use, era5_vars=era5_vars, do_pyr=get_pyr_data, lf_vars=lf_vars, do_feds=rasterize_feds,
        verbose=True, plot=plot_sources, batch_plot=False, all_plot=False, del_sources=gen_util.data_sources,
        del_intermediate=False
    )