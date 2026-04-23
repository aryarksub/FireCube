import os
import pandas as pd
import geopandas as gpd
import time

from era5.era5 import driver_era5
from landfire.landfire import driver_landfire
from pyregence.pyregence import driver_pyregence
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util
import util.validation_util as valid_util

# DataFrame of all fires stored by FEDS/MTBS
firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)

def process_single_fire(fid, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, plot=[], 
                        batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False,
                        feds_direct_final_grid=False, remove_old_date=time.mktime((2025, 1, 1, 0, 0, 0, 0, 0, -1))):
    """
    Driver function to process a single fire. Processing includes downloading data from different sources (ERA5, Pyregence,
    LANDFIRE, FEDS), cropping/rasterizing/resampling as needed, and creating plots of the layers.

    Args:
        fid (str): Event ID of the fire to process.
        era5_vars (list, optional): List of layers to download from ERA5. Full list of available layers
         can be found at https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Parameterlistings. 
         Defaults to [].
        do_pyr (bool, optional): True if Pyregence layers should be downloaded; False otherwise. Defaults to True.
        lf_vars (list, optional): List of layers to download from LANDFIRE. Full list of available layers can be found at
         https://www.landfire.gov/sites/default/files/documents/LF_Data_Dictionary.pdf in Section 4. Defaults to [].
        do_feds (bool, optional): True if FEDS fireline data should be obtained; False otherwise. Defaults to True.
        verbose (bool, optional): True if descriptive messages should be printed; False otherwise. Defaults to False.
        plot (list, optional): List of sources for which plots should be generated (see gen_util.data_sources). Defaults to [].
        batch_plot (bool, optional): True if plots should be generated for layer categories; False otherwise. Defaults to False.
        all_plot (bool, optional): True if a plot should be generated with all layers; False otherwise. Defaults to False.
        del_sources (list, optional): List of sources for which temporarily created data should be deleted. 
         Defaults to gen_util.data_sources.
        del_intermediate (bool, optional): True if intermediate downloaded data (e.g. zip/tar files) should be deleted;
         False otherwise. Defaults to False.
        feds_direct_final_grid (bool, optional): True to run FEDS/FRP rasterization directly on the final
         canonical grid for A/B testing; False for legacy rasterize+resample+pad. Defaults to False.
        remove_old_date (float, optional): Unix timestamp (seconds since epoch) before which files should be removed. Defaults to 1/1/2025.

    Raises:
        ValueError: Occurs when the downloaded data is invalid (e.g. layer value is outside of min/max range).
    """
    if verbose:
        print(f'Processing fire {fid}')
    
    # Only process a fire if it has FEDS data, specifically for fire area
    gdf_farea_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    if gdf_farea_rd is None:
        if verbose:
            print(f'FEDS file for fire {fid} does not exist - no processing will be done')
        return
    
    gen_util.create_dirs_for_fire(fid)
    # West/South/East/North bounds for the fire
    bnds = proc_util.bufferbnds(gdf_farea_rd.total_bounds, res=0.005, bufgd=1)
    df_t = pd.to_datetime(gdf_farea_rd.t)
    # Add time-buffer to DataFrame
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    # print('Original fire start/end', df_t.min(), df_t.max())
    # print('Buffered fire start/end', df_t_with_buffer.min(), df_t_with_buffer.max())

    # Convert gdf_farea to EPSG:5070
    gdf_farea_5070 = gdf_farea_rd.to_crs('EPSG:5070') 
    bounds_5070 = gdf_farea_5070.total_bounds

    # Get fire center
    fire_row = firelist[firelist['Event_ID'] == fid]
    center_lat = round((fire_row['lat0'].values[0] + fire_row['lat1'].values[0]) / 2, 2)
    center_lon = round((fire_row['lon0'].values[0] + fire_row['lon1'].values[0]) / 2, 2)
   
    # Convert from LST to UTC (inverting FEDS method for computing local time)
    # Add one hour for alignnment from LST to UTC, but this is a rough estimate since the actual time shift depends on the 
    # longitude of the fire and the time of year (daylight savings time).
    conversion_delta = pd.to_timedelta(1, unit="hours") - pd.to_timedelta(round(center_lon / 15), unit="hours")
    df_t_with_buffer = df_t_with_buffer + conversion_delta

    # print('Buffered fire start/end after LST to UTC conversion', df_t_with_buffer.min(), df_t_with_buffer.max())

    fire_start = pd.Timestamp(df_t_with_buffer.min().normalize())
    fire_end = pd.Timestamp(df_t_with_buffer.max().normalize()) + pd.Timedelta(hours=23)
    fire_hours = int( (fire_end - fire_start).total_seconds() / 3600)

    # print('Fire start/end after normalization + num hours', fire_start, fire_end, fire_hours)

    # ERA5 download
    if len(era5_vars) != 0:
        if verbose: print(f'Getting ERA5 data for fire {fid}')
        
        driver_era5(
            fid, era5_vars, df_t_with_buffer, bnds, 
            gen_util.get_era5_nc_filename(fid), 
            plot_types=gen_util.var_types if gen_util.subdir_era5 in plot else []
        )
    else:
        if verbose: print(f'No ERA5 variables specified; not getting ERA5 data for fire {fid}')

    # Pyregence download
    if do_pyr:
        if verbose: print(f'Getting Pyregence data for fire {fid}')
        driver_pyregence(
            fid, (center_lat, center_lon), fire_start, fire_hours, 
            plot_types=gen_util.var_types if gen_util.subdir_pyr in plot else []
        )
    else:
        if verbose: print(f'Skipping Pyregence data for fire {fid}')

    # LANDFIRE download
    if len(lf_vars) != 0:
        if verbose: print(f'Getting LANDFIRE data for fire {fid}')

        driver_landfire(
            fid, lf_vars, bnds, fire_start, 
            plot_types=gen_util.var_types if gen_util.subdir_lf in plot else []
        )
    else:
        if verbose: print(f'No LANDFIRE variables specified; not getting LANDFIRE data for fire {fid}')       

    # Crop era5/pyr/lf tifs to just surround the fire perim (ignore FEDS data for cropping since it needs to be padded)
    # The output TIFs here are pre-processing and may not exist yet
    non_feds_input_tifs, pre_proc_non_feds_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(
        fid, exclude=[gen_util.subdir_feds]
    )
    # The output TIFs here are post-processing and may not exist if the fire is being processed for the first time
    post_proc_non_feds_output_tifs = gen_util.get_all_tifs_in_output_dir_for_fire(fid)

    # print(len(non_feds_input_tifs), len(pre_proc_non_feds_output_tifs), len(post_proc_non_feds_output_tifs))
    # print('INPUT', non_feds_input_tifs)
    # print('PRE-PROC', pre_proc_non_feds_output_tifs)
    # print('POST-PROC', post_proc_non_feds_output_tifs)

    # If no output TIFs do not exist, stop processing (there is an issue)
    if len(post_proc_non_feds_output_tifs) == 0 and len(pre_proc_non_feds_output_tifs) == 0:
        if verbose:
            print(f'No output TIFs were created - stopping processing for fire {fid}')
        return
    # If we have pre-processing output TIF file names, then do the necessary cropping
    elif len(pre_proc_non_feds_output_tifs) != 0:
        non_feds_output_tifs = pre_proc_non_feds_output_tifs
        # Make all non-FEDS TIFs centered and cropped to the same region
        proc_util.center_and_crop_tifs_to_same_area(non_feds_input_tifs, non_feds_output_tifs, bounds_5070)
    # Otherwise, we have post-processing output TIFs, so continue with those
    else:
        non_feds_output_tifs = post_proc_non_feds_output_tifs
    
    non_feds_output_tifs = list(set(non_feds_output_tifs).union(set(post_proc_non_feds_output_tifs)))

    # print(len(non_feds_output_tifs), non_feds_output_tifs)

    # Bounding box for all variable/layer tifs is the same, so we can just take the box for the first tif
    largest_var_tif_bounds = proc_util.get_tif_bounds(non_feds_output_tifs[0])

    # print('BOUNDS 5070:', bounds_5070)
    # for var_tif in non_feds_output_tifs:
    #     print(proc_util.get_tif_bounds(var_tif))

    # FEDS + FRP download/rasterization
    if do_feds:
        if verbose: print(f'Getting/Rasterizing FEDS data for fire {fid}')
        feds_util.driver_feds(
            fid, largest_var_tif_bounds, res=300.0, fire_start=fire_start, num_hours=fire_hours, 
            plot_orig=True if gen_util.subdir_feds in plot else False,
            use_prev=False, conv_delta=conversion_delta,
            direct_to_final_grid=feds_direct_final_grid
        )
        feds_util.driver_frp(
            fid, largest_var_tif_bounds, res=300.0, fire_start=fire_start, num_hours=fire_hours, 
            plot_orig=True if gen_util.subdir_feds in plot else False,
            use_prev=False, conv_delta=conversion_delta,
            direct_to_final_grid=feds_direct_final_grid
        )
    else:
        if verbose: print(f'Skipping FEDS data for fire {fid}')

    # Get output TIF names (as stored in output/cubes directory)
    all_variable_input_tifs, all_variable_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(fid)

    # Convert all null values (-1/-9999) to np.nan
    for var_tif in all_variable_output_tifs:
        # print(f'Processing variable TIF: {var_tif}')
        # If the output file does not exist, just skip it
        try:
            # phi values are allowed to be -1, but this is considered null data for other variables
            if 'phi' in var_tif:
                proc_util.convert_null_values_to_nan(in_tif=var_tif, out_tif=var_tif)
            else:
                proc_util.convert_null_values_to_nan(in_tif=var_tif, out_tif=var_tif, null_values=[valid_util.GLOBAL_NULL_VALUE, -1])
        except:
            pass

    # # Validate all downloaded data
    # valid_data, invalid_layers, invalid_spatial_layers, missing_layers = valid_util.validate_one_fire_data(fid)
    # if not valid_data:
    #     if verbose:
    #     #     print('Layers with invalid data:')
    #     #     for layer in invalid_layers:
    #     #         print(layer)
    #     #     print('Layers with invalid spatial info:')
    #     #     for layer in invalid_spatial_layers:
    #     #         print(layer)
    #     #     print('Missing layers:')
    #     #     for layer in missing_layers:
    #     #         print(layer)
    #     # raise ValueError(f'Fire {fid} has invalid data')
    #         print(f'Fire {fid} has invalid data')

    # Create batch plots
    if batch_plot:
        for batch in gen_util.data_batches:
            if batch != gen_util.subdir_vis:
                if verbose:
                    print(f'Generating plot for batch {batch} - fire {fid}')
                try:
                    # Add 1 to fire_hours since we want to go up to midnight of next day
                    gen_util.create_multi_animation_for_dir(
                        os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, batch),
                        gen_util.get_output_data_filename(fid, batch, gen_util.subdir_vis),
                        fire_start, num_steps=fire_hours+1
                    )
                except:
                    pass

    # Create one large plot with all variables
    if all_plot: 
        if verbose:
            print(f'Generating plot for all variables - fire {fid}')
        gen_util.create_multi_animation_from_tifs(
            all_variable_output_tifs,
            os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, gen_util.subdir_vis, 'all.mp4'),
            start_time=fire_start
        )

    # Remove temporary data files
    gen_util.remove_temp_dir_files(
        fid, 
        del_dir_types=gen_util.dir_types,
        del_data_sources=del_sources,
        del_var_types=gen_util.var_types,
        remove_intermediate=del_intermediate,
        verbose=verbose
    )

    for subdir in [gen_util.subdir_hrc, gen_util.subdir_fuel_topo, gen_util.subdir_landfire]:
        data_dir = os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, subdir)
        gen_util.remove_old_files(data_dir, cutoff_timestamp=remove_old_date)

    # gen_util.add_to_new_fires_list(fid, '2018UTVA')

def process_multiple_fires(fid_list=[], fid_file=None, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, 
                           plot=[], batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False, 
                           feds_direct_final_grid=False, remove_old_date=time.mktime((2025, 1, 1, 0, 0, 0, 0, 0, -1))):
    """
    Driver function to process multiple fire. Processing includes downloading data from different sources (ERA5, Pyregence,
    LANDFIRE, FEDS), cropping/rasterizing/resampling as needed, and creating plots of the layers. The multiple fires to process
    can be given as either a list of strings (fire event IDs) or a file where each fire event ID is on a separate line.

    Args:
        fid_list (list, optional): List of fire event IDs. Defaults to [].
        fid_file (str, optional): Name of file with fire event IDs. Defaults to None.
        era5_vars (list, optional): List of layers to download from ERA5. Full list of available layers
         can be found at https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#heading-Parameterlistings. 
         Defaults to [].
        do_pyr (bool, optional): True if Pyregence layers should be downloaded; False otherwise. Defaults to True.
        lf_vars (list, optional): List of layers to download from LANDFIRE. Full list of available layers can be found at
         https://www.landfire.gov/sites/default/files/documents/LF_Data_Dictionary.pdf in Section 4. Defaults to [].
        do_feds (bool, optional): True if FEDS fireline data should be obtained; False otherwise. Defaults to True.
        verbose (bool, optional): True if descriptive messages should be printed; False otherwise. Defaults to False.
        plot (list, optional): List of sources for which plots should be generated (see gen_util.data_sources). Defaults to [].
        batch_plot (bool, optional): True if plots should be generated for layer categories; False otherwise. Defaults to False.
        all_plot (bool, optional): True if a plot should be generated with all layers; False otherwise. Defaults to False.
        del_sources (list, optional): List of sources for which temporarily created data should be deleted. 
         Defaults to gen_util.data_sources.
        del_intermediate (bool, optional): True if intermediate downloaded data (e.g. zip/tar files) should be deleted;
         False otherwise. Defaults to False.
        feds_direct_final_grid (bool, optional): True to run FEDS/FRP rasterization directly on the final canonical grid for A/B testing
        remove_old_date (float, optional): Unix timestamp (seconds since epoch) before which files should be removed. Defaults to 1/1/2025.
    """
    if len(fid_list) == 0 and fid_file is None:
        if verbose:
            print('No FIDs given - no processing will be done')
        return
    
    count = 0
    if len(fid_list) > 0:
        if verbose:
            print('Processing fires given in list argument')
        for fid in fid_list:
            process_single_fire(
                fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate,
                feds_direct_final_grid=feds_direct_final_grid, remove_old_date=remove_old_date
            )
            count += 1
            if count % 10 == 0 or count == len(fid_list):
                print(f'******* STATUS CHECK: Processed {count} fires *******')
    else:
        if verbose:
            print('Processing fires given in file-path argument')
        try:
            with open(fid_file, 'r', encoding='utf-8') as file:
                for line in file:
                    fid = line.strip()
                    process_single_fire(
                        fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate,
                        feds_direct_final_grid=feds_direct_final_grid, remove_old_date=remove_old_date
                    )
                    count += 1
                    if count % 20 == 0:
                        print(f'Processed {count} fires')
        except:
            if verbose:
                print(f'Error when reading file {fid_file} - no processing will be done')

def random_select_fids(n=5, size_threshold=None, min_size=1000, duration_threshold=None, year=None, method='random', skip_exist=True):
    """
    Select n fire event IDs from the entire list of fires based on a given size and duration threshold. Selection
    can be done either randomly (method='random') or by selecting the largest fires (method='size'). If the number
    of available fires to download is less than n, then the output size will be the number of available fires
    instead of n.

    Args:
        n (int, optional): Number of fire event IDs to return. Defaults to 10.
        size_threshold (int, optional): Upper bound of fire size in terms of burned area. Defaults to None.
        min_size (int, optional): Lower bound of fire size in terms of burned area. Defaults to 1000 acres.
        duration_threshold (int, optional): Upper bound of fire duration in days. Defaults to None.
        year (int, optional): Year for which to select fires. Defaults to None.
        method (str, optional): Method by which FIDs should be selected (random, size). Defaults to 'random'.
        skip_exist (bool, optional): True if FIDs that already have downloaded data should be exempt from 
         selection; False otherwise. Defaults to True.

    Raises:
        ValueError: Occurs when the size threshold is less than 1000 or the duration threshold is less than 1.
        NotImplementedError: Occurs when the given selection method is not supported.

    Returns:
        list: List of FIDs that meet the given size/duration thresholds.
    """
    existing_fids = [fid for fid in os.listdir(os.path.join(gen_util.dir_output, gen_util.dir_cubes))]
    remote = []
    try:
        with open(os.path.join('src', f'temp_new_fids{"_" + str(year) if year is not None else ""}.txt'), 'r', encoding='utf-8') as f:
            for line in f:
                remote.append(line.strip())
    except:
        pass
    existing_fids = list(set(existing_fids).union(set(remote)))

    fires_df = firelist[(~firelist['Event_ID'].str.contains('HI', na=False))]

    # If we should avoid selecting FIDs for fires whose data already exists, then mask these FIDs out
    if skip_exist:
        fires_df = fires_df[~fires_df['Event_ID'].isin(existing_fids)]

    # Filter based on size threshold if it is provided
    if size_threshold is not None:
        if size_threshold >= 1000:
            fires_df = fires_df[fires_df['BurnBndAc'] < size_threshold]
        else:
            raise ValueError('Size threshold for filtering must be at least 1000 acres')
    
    fires_df = fires_df[fires_df['BurnBndAc'] >= min_size]
        
    # Filter based on duration threshold if it is provided
    if duration_threshold is not None:
        if duration_threshold >= 1:
            fires_df = fires_df[
                (pd.to_datetime(fires_df['ted']).dt.normalize() - pd.to_datetime(fires_df['tst'])).dt.days.between(0, duration_threshold)
            ]
        else:
            raise ValueError('Duration threshold for filtering must be at least 1 day')
        
    # Filter based on year if it is provided
    if year is not None:
        fires_df = fires_df[fires_df['Year'] == year]

    print(fires_df.shape[0], 'fires meet the given size/duration thresholds and do not already have downloaded data')
    # Select n FIDs if there are that many; otherwise, select all available FIDs
    num_to_select = min(n, fires_df.shape[0])
    
    if method == 'random':
        sample_fires = fires_df.sample(n=num_to_select)
        sample_fids = list(sample_fires['Event_ID'])
    elif method == 'size':
        sorted_df = fires_df.sort_values(by='BurnBndAc', ascending=False)
        top_n = sorted_df.head(num_to_select)
        sample_fids = list(top_n['Event_ID'])
    elif method == 'size_reverse':
        sorted_df = fires_df.sort_values(by='BurnBndAc', ascending=True)
        top_n = sorted_df.head(num_to_select)
        sample_fids = list(top_n['Event_ID'])
    else:
        raise NotImplementedError(f'The given selection method {method} is not supported')

    return sample_fids

def get_all_existing_fids():
    """
    Get list of all fire event IDs for which data has already been downloaded and processed 
    (i.e. there is a folder for the fire in output/cubes).

    Returns:
        list: List of fire event IDs for which data has already been downloaded and processed.
    """
    existing_fids = [fid for fid in os.listdir(os.path.join(gen_util.dir_output, gen_util.dir_cubes))]
    return existing_fids

def get_fids_from_file(file_path):
    with open(file_path, 'r') as f:
        fids = [line.strip() for line in f if line.strip()]
    return fids

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'
    caldor_id = 'CA3858612053820210815'
    temp_id = 'CA3338011743020180706'

    # To skip ERA5 download, set era5_vars = []
    era5_vars = ['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']
    # To skip Pyregence download, set do_pyr = False
    get_pyr_data = True
    # To skip LANDFIRE download, set lf_vars = []
    lf_vars = ['ASP', 'ELEV', 'SLPD', 'EVT', 'FBFM13', 'FBFM40', 'ROADS']
    # To skip FEDS rasterization, set do_feds = False
    rasterize_feds = True
    # To skip plotting, set plot_sources = []
    plot_sources = []
    existing_fids = get_all_existing_fids()
    # for x in existing_fids:
    #     print(x)
    
    # True : FIDs should be randomly selected
    # False: Use hard-coded FID(s)
    do_sample_fids = True
    feds_direct_final_grid = True
    cutoff_date = time.mktime((2026, 4, 1, 0, 0, 0, 0, 0, -1))

    if do_sample_fids:
        fids_to_use = random_select_fids(n=1, min_size=1000, size_threshold=2000000, duration_threshold=150, method='random')
    else:
        fids_to_use = [temp_id] #existing_fids
    
    # Standard procedure is to use del_sources=gen_util.data_sources to delete temporary created data files
    process_multiple_fires(
        fid_list=fids_to_use, era5_vars=era5_vars, do_pyr=get_pyr_data, lf_vars=lf_vars, do_feds=rasterize_feds,
        verbose=True, plot=plot_sources, batch_plot=False, all_plot=False, del_sources=gen_util.data_sources,
        del_intermediate=False, feds_direct_final_grid=feds_direct_final_grid,
        remove_old_date=cutoff_date
    )
