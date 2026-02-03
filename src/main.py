import os
import pandas as pd
import geopandas as gpd

from era5.era5 import driver_era5
from landfire.landfire import driver_landfire
from pyregence.pyregence import driver_pyregence
import util.feds_util as feds_util
import util.general_util as gen_util
import util.processing_util as proc_util
import util.validation_util as valid_util

# DataFrame of all fires stored by FEDS/MTBS
firelist = pd.read_csv(feds_util.feds_firelist, index_col=0)

def process_single_fire(fid, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, plot=[], batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False):
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

    Raises:
        ValueError: Occurs when the downloaded data is invalid (e.g. layer value is outside of min/max range).
    """
    if verbose:
        print(f'Processing fire {fid}')
    
    # Only process a fire if it has FEDS data, specifically for fire perimeter
    gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd = feds_util.read_1fire(fid)
    if gdf_fperim_rd is None:
        if verbose:
            print(f'FEDS file for fire {fid} does not exist - no processing will be done')
        return
    
    gen_util.create_dirs_for_fire(fid)
    # West/South/East/North bounds for the fire
    bnds = proc_util.bufferbnds(gdf_fperim_rd.total_bounds, res=0.005, bufgd=1)
    df_t = pd.to_datetime(gdf_fperim_rd.t)
    # Add time-buffer to DataFrame
    df_t_with_buffer = proc_util.add_time_buffers(df_t)

    # Convert gdf_perim to EPSG:5070
    gdf_fperim_5070 = gdf_fperim_rd.to_crs('EPSG:5070') 
    bounds_5070 = gdf_fperim_5070.total_bounds

    # Get fire center
    fire_row = firelist[firelist['Event_ID'] == fid]
    center_lat = round((fire_row['lat0'].values[0] + fire_row['lat1'].values[0]) / 2, 2)
    center_lon = round((fire_row['lon0'].values[0] + fire_row['lon1'].values[0]) / 2, 2)
   
    #Convert from LST to UTC (inverting FEDS method for computing local time)
    df_t_with_buffer = df_t_with_buffer - pd.to_timedelta(round(center_lon / 15), unit="hours")

    fire_start = pd.Timestamp(df_t_with_buffer.min().normalize())
    fire_end = pd.Timestamp(df_t_with_buffer.max().normalize()) + pd.Timedelta(hours=23)
    fire_hours = int( (fire_end - fire_start).total_seconds() / 3600)

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
    non_feds_input_tifs, non_feds_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(
        fid, exclude=[gen_util.subdir_feds]
    )
    # If there are no output tifs, then stop processing
    if len(non_feds_output_tifs) == 0:
        if verbose:
            print(f'No output TIFs were created - stopping processing for fire {fid}')
        return

    # Make all non-FEDS TIFs centered and cropped to the same region
    proc_util.center_and_crop_tifs_to_same_area(non_feds_input_tifs, non_feds_output_tifs, bounds_5070)

    # Bounding box for all variable/layer tifs is the same, so we can just take the box for the first tif
    largest_var_tif_bounds = proc_util.get_tif_bounds(non_feds_output_tifs[0])

    # FEDS + FRP download/rasterization
    if do_feds:
        if verbose: print(f'Getting/Rasterizing FEDS data for fire {fid}')
        feds_util.driver_feds(
            fid, largest_var_tif_bounds, res=300.0, fire_start=fire_start, num_hours=fire_hours, 
            plot_orig=True if gen_util.subdir_feds in plot else False
        )
        feds_util.driver_frp(
            fid, largest_var_tif_bounds, res=300.0, fire_start=fire_start, num_hours=fire_hours, 
            plot_orig=True if gen_util.subdir_feds in plot else False
        )
    else:
        if verbose: print(f'Skipping FEDS data for fire {fid}')

    # Get output TIF names (as stored in output/cubes directory)
    all_variable_input_tifs, all_variable_output_tifs = gen_util.get_all_var_and_output_tifs_for_fire(fid)

    # Convert all null values (-1/-9999) to np.nan
    for var_tif in all_variable_output_tifs:
        # phi values are allowed to be -1, but this is considered null data for other variables
        if 'phi' in var_tif:
            proc_util.convert_null_values_to_nan(in_tif=var_tif, out_tif=var_tif)
        else:
            proc_util.convert_null_values_to_nan(in_tif=var_tif, out_tif=var_tif, null_values=[valid_util.GLOBAL_NULL_VALUE, -1])

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
                    gen_util.create_multi_animation_for_dir(
                        os.path.join(gen_util.dir_output, gen_util.dir_cubes, fid, batch),
                        gen_util.get_output_data_filename(fid, batch, gen_util.subdir_vis),
                        fire_start
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

def process_multiple_fires(fid_list=[], fid_file=None, era5_vars=[], do_pyr=True, lf_vars=[], do_feds=True, verbose=False, plot=[], batch_plot=False, all_plot=False, del_sources=gen_util.data_sources, del_intermediate=False):
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
            process_single_fire(fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate)
            count += 1
            if count % 20 == 0:
                print(f'Processed {count} fires')
    else:
        if verbose:
            print('Processing fires given in file-path argument')
        try:
            with open(fid_file, 'r', encoding='utf-8') as file:
                for line in file:
                    fid = line.strip()
                    process_single_fire(fid, era5_vars, do_pyr, lf_vars, do_feds, verbose, plot, batch_plot, all_plot, del_sources, del_intermediate)
                    count += 1
                    if count % 20 == 0:
                        print(f'Processed {count} fires')
        except:
            if verbose:
                print(f'Error when reading file {fid_file} - no processing will be done')

def random_select_fids(n=5, size_threshold=None, duration_threshold=None, method='random', skip_exist=True):
    """
    Select n fire event IDs from the entire list of fires based on a given size and duration threshold. Selection
    can be done either randomly (method='random') or by selecting the largest fires (method='size'). If the number
    of available fires to download is less than n, then the output size will be the number of available fires
    instead of n.

    Args:
        n (int, optional): Number of fire event IDs to return. Defaults to 10.
        size_threshold (int, optional): Upper bound of fire size in terms of burned area. Defaults to None.
        duration_threshold (int, optional): Upper bound of fire duration in days. Defaults to None.
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

    # TODO: Revert this after LANDFIRE versioning has been fixed
    fires_df = firelist[(~firelist['Event_ID'].str.contains('HI|AK', na=False)) & (firelist['Year'] >= 2020) & (firelist['Year'] <= 2022)]

    # If we should avoid selecting FIDs for fires whose data already exists, then mask these FIDs out
    if skip_exist:
        fires_df = fires_df[~fires_df['Event_ID'].isin(existing_fids)]

    # Filter based on size threshold if it is provided
    if size_threshold is not None:
        if size_threshold >= 1000:
            fires_df = fires_df[fires_df['BurnBndAc'] < size_threshold]
        else:
            raise ValueError('Size threshold for filtering must be at least 1000 acres')
        
    # Filter based on duration threshold if it is provided
    if duration_threshold is not None:
        if duration_threshold >= 1:
            fires_df = fires_df[
                (pd.to_datetime(fires_df['ted']).dt.normalize() - pd.to_datetime(fires_df['tst'])).dt.days.between(1, duration_threshold)
            ]
        else:
            raise ValueError('Duration threshold for filtering must be at least 1 day')
        
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

if __name__=='__main__':
    creek_id = 'CA3720111927220200905'
    zogg_id = 'CA4054112256820200927'
    temp_id = 'AZ3338011031820240727'

    era5_vars = ['surface_pressure', 'total_precipitation', '2m_temperature', '2m_dewpoint_temperature']
    get_pyr_data = True
    lf_vars = ['ASP', 'ELEV', 'SLPD', 'EVT', 'FBFM13', 'FBFM40', 'ROADS']
    rasterize_feds = True
    plot_sources = []
    
    # True : FIDs should be randomly selected
    # False: Use hard-coded FID(s)
    do_sample_fids = False

    if do_sample_fids:
        fids_to_use = random_select_fids(n=1, size_threshold=100000, duration_threshold=28, method='random')
    else:
        fids_to_use = [temp_id]
    
    # Standard procedure is to use del_sources=gen_util.data_sources to delete temporary created data files
    process_multiple_fires(
        fid_list=fids_to_use, era5_vars=era5_vars, do_pyr=get_pyr_data, lf_vars=lf_vars, do_feds=rasterize_feds,
        verbose=True, plot=plot_sources, batch_plot=False, all_plot=False, del_sources=gen_util.data_sources,
        del_intermediate=False
    )
