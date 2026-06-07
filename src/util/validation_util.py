import pandas as pd
import numpy as np
import os
import rasterio

output_cubes_dir = os.path.join('output', 'cubes')
layer_range_file = os.path.join('inputData', 'layer_ranges.csv')
GLOBAL_NULL_VALUE = -9999
BATCH_RESOLUTION_MAP = {
    "fire_spread" : 300,
    "fuel_structure" : 30,
    "high_res_climate" : 600,
    "veg_fm_topo" : 30,
    "low_res_climate" : 9000
}

def get_layer_range_data():
    """
    Get dictionary mapping layers to important information. Keys are layer names. Values are dictionaries that have
    keys "category", "min", "max", "nan". These represent the category of the layer, min/max values, and a value
    that is treated as a null value for the layer.

    Returns:
        dict: Dictionary of layer information.
    """
    layer_range_df = pd.read_csv(layer_range_file)
    key = 'layer'
    return layer_range_df.set_index(key).to_dict(orient='index')

def find_layer_file_containing(abbrev, layer_dir):
    """
    Find the layer file corresponding to the given abbreviation in a directory.

    Args:
        abbrev (str): Abbreviation of the layer to search for.
        layer_dir (str): Path of the directory in which we should search.

    Returns:
        str: Name of the corresponding layer file.
    """
    for file in sorted(os.listdir(layer_dir)):
        if abbrev in file and os.path.isfile(os.path.join(layer_dir, file)):
            return file
    return None

def get_tif_spatial_properties(src):
    """
    Given an opened TIF file (dataset reader), i.e. the output of rasterio.open(...), get the spatial area properties
    of the TIF. This includes the coordinates of the top-left pixel, the spatial width, and spatial height.
    
    Args:
        src (rasterio.DatasetReader): Dataset reader of an opened TIF file

    Returns:
        tuple: Tuple (x, y, w, h) where (x,y) is the coordinates of the top-left corner of the top-left pixel,
        w is the spatial width of the TIF, and h is the spatial height of the TIF
    """

    # transform = Affine object whose matrix is represented by 6 parameters:
    # a = x-direction pixel width
    # b = row rotation (typically 0)
    # c = x-coordinate of upper-left corner of upper-left pixel
    # d = column rotation (typically 0)
    # e = y-direction pixel height (negative value for north-up images)
    # f = y-coordinate of upper-left corner of upper-left pixel
    transform = src.transform
    width = src.width
    height = src.height

    top_left = (transform.c, transform.f)
    x_max = transform.c + width * transform.a
    y_min = transform.f + height * transform.e

    spatial_width = abs(x_max - transform.c)
    spatial_height = abs(transform.f - y_min)

    return *top_left, spatial_width, spatial_height

def check_band_count(layer_tif1, layer_tif2):
    """
    Determine whether the band count of layer_tif2 is either 1 or equal to the band count of layer_tif1. Due to changes
    so that fire_spread layers are recorded every 12 hours, we also all band count of layer_tif1 to be either 4 or 5
    less than the band count of layer_tif2 (assuming layer_tif1 is a fire_spread layer).

    Args:
        layer_tif1 (str): Path to the first TIF file.
        layer_tif2 (str): Path to the second TIF file.

    Returns:
        tuple: Tuple of (is_valid, tif1_band_count, tif2_band_count) where is_valid is a boolean indicating whether
         the band counts are valid (i.e. either 1 or equal or correct relative to hourly vs 12-hourly temporal difference).
    """
    with rasterio.open(layer_tif1) as src_tif1:
        tif1_band_count = src_tif1.count

        with rasterio.open(layer_tif2) as src_tif2:
            tif2_band_count = src_tif2.count

            return (
                tif2_band_count in {1, tif1_band_count} or tif2_band_count // 12 - tif1_band_count in {4,5}, 
                tif1_band_count, 
                tif2_band_count
            )

def validate_one_fire_data(fid, layer_data_dict=None):
    """
    Determine whether the downloaded data for the given fire (fid) is valid. Validation is done by checking that the data
    for each layer is either within the min/max range defined in layer_data_dict or the corresponding null value (also
    defined in layer_data_dict).

    Args:
        fid (str): Fire event ID for which to validate data.
        layer_data_dict (dict, optional): Dictionary of layer information (see get_layer_range_data() for dictionary
         format). Defaults to None.

    Raises:
        NotADirectoryError: Occurs when there is no directory for the given fire (fid) or there is no directory for
         a certain layer category.

    Returns:
        tuple: Tuple of a boolean and four lists. The boolean takes a True value if the fire data is valid; False otherwise.
         The first list contains tuples of the form (layer, band, indices). These store the invalid layers, the time band at
         which the data is invalid, and the indices in the 2D TIF data array that are invalid. The second list contains
         the names of layers that have unexpected spatial properties as compared to other layers. It also contains the 
         names of layers whose correspond TIF has unexpected resolution. The third list contains the names of layers that
         are missing (i.e. no TIF files). The fourth list contains layers with invalid band counts.
    """

    # If no layer data dictionary is specified, then obtain it with get_layer_range_data()
    layer_data_dict = layer_data_dict if layer_data_dict is not None else get_layer_range_data()

    # Ensure there is an output folder for the given FID
    fire_folder = os.path.join('output', 'cubes', fid)
    if not os.path.exists(fire_folder):
        raise NotADirectoryError(f'No folder for fire {fid} exists')
    
    invalid_layers = []
    evt_found = False

    expected_spatial_properties = None
    invalid_layers_spatial = []
    missing_layers = []
    invalid_layers_band_count = []

    farea_tif = os.path.join(fire_folder, 'fire_spread', 'farea.tif')
    
    for layer, layer_data in layer_data_dict.items():
        # Keep support for layer range files with multiple EVT aliases.
        if 'evt' in layer and evt_found:
            continue

        category = layer_data['category']
        min_val, max_val, null_val = layer_data['min'], layer_data['max'], layer_data['nan']
        
        # Ensure there is a folder for the layer's category
        layer_folder = os.path.join(fire_folder, category)
        if not os.path.exists(layer_folder):
            raise NotADirectoryError(f'No folder for layer category {category} and fire {fid} exists')
        
        # VFMT layers may include prefixes/suffixes in their TIF names.
        if category == 'veg_fm_topo':
            layer_file = find_layer_file_containing(layer, layer_folder)
            if layer_file is None and layer in {'fbfm13', 'fbfm40'}:
                layer_file = find_layer_file_containing(layer.replace('fbfm', 'f'), layer_folder)
            if layer_file is None:
                layer_tif = None
            else:
                layer_tif = os.path.join(layer_folder, layer_file)
        else:
            layer_tif = os.path.join(layer_folder, f'{layer}.tif')
        
        if layer_tif is None or not os.path.exists(layer_tif):
            missing_layers.append(layer)
            continue
        
        if 'evt' in layer:
            evt_found = True

        is_valid, farea_band_count, layer_band_count = check_band_count(farea_tif, layer_tif)
        if not is_valid:
            invalid_layers_band_count.append((layer, farea_band_count, layer_band_count))
        
        # Data validation for each band in TIF
        with rasterio.open(layer_tif) as src:
            # Store spatial properties of first layer TIF
            if expected_spatial_properties is None:
                expected_spatial_properties = get_tif_spatial_properties(src=src)
            # Validate all other TIF spatial properties with the stored properties
            else:
                spatial_properties = get_tif_spatial_properties(src=src)
                if not np.allclose(
                    np.array(expected_spatial_properties, dtype=float).ravel(), 
                    np.array(spatial_properties, dtype=float).ravel()
                ):
                    invalid_layers_spatial.append((layer, expected_spatial_properties, spatial_properties))

            # Make sure layer TIF has the expected resolution based on the layer category (batch)
            expected_res = BATCH_RESOLUTION_MAP[category]
            actual_res = src.transform.a
            if expected_res != actual_res:
                invalid_layers_spatial.append((layer, expected_res, actual_res))

            for band in range(src.count):
                # Band counting is 1-indexed
                data = src.read(band + 1)
                # Allow null value if specified
                if null_val:
                    # TODO: remove null/global null checks, if null: then allow data == np.nan
                    valid_data_mask = ((data >= min_val) & (data <= max_val)) | (np.isnan(data))
                else:
                    valid_data_mask = (data >= min_val) & (data <= max_val)
                
                # Store the first occurrence of invalid data
                if not np.all(valid_data_mask):
                    invalid_indices = np.where(~valid_data_mask)
                    invalid_locations = list(zip(invalid_indices[0], invalid_indices[1]))
                    invalid_layers.append((layer, band+1, invalid_locations))
                    break
        
    return (
        len(invalid_layers) == 0 and len(invalid_layers_spatial) == 0 and len(invalid_layers_band_count) == 0, 
        invalid_layers,
        invalid_layers_spatial, 
        missing_layers, 
        invalid_layers_band_count
    )


def validate_existing_fires(layer_data_dict=None, verbose=False):
    """
    Perform data validation for all fires that are stored in the output directory. This is done by simply iterating
    through the directory and running validate_one_fire_data().

    Args:
        layer_data_dict (dict, optional): Dictionary of layer information (see get_layer_range_data() for dictionary
         format). Defaults to None.
        verbose (bool, optional): True if descriptive messages should be printed; False otherwise. Defaults to False.

    Returns:
        tuple: Tuple with three lists: fires with invalid data, fires with missing data, and fires with invalid
        band counts.
    """
    # If no layer data dictionary is specified, then obtain it with get_layer_range_data()
    layer_data_dict = layer_data_dict if layer_data_dict is not None else get_layer_range_data()

    invalid_fires = []
    fires_with_missing_data = []
    invalid_fires_band_count = []
    processed_count = 0

    # Note that folders are named with the fire event ID (so folder = FID)
    for folder in os.listdir(output_cubes_dir):
        if os.path.isdir(os.path.join(output_cubes_dir, folder)):
            valid_data, invalid_layers, invalid_layers_spatial, missing_layers, invalid_layers_band_count = validate_one_fire_data(folder, layer_data_dict)
            if not valid_data:
                invalid_fires.append((folder, invalid_layers, invalid_layers_spatial, missing_layers, invalid_layers_band_count))

                if verbose:
                    print(f'Fire {folder} has invalid data')
            
            if len(missing_layers) > 0:
                fires_with_missing_data.append((folder, missing_layers))
            if len(invalid_layers_band_count) > 0:
                invalid_fires_band_count.append((folder, invalid_layers_band_count))
            
            processed_count += 1

            if verbose and processed_count % 10 == 0:
                print(f'Validated {processed_count} fires')
    
    if verbose:
        print(f'Done with validation ... Validated {processed_count} fires')
    
    return invalid_fires, fires_with_missing_data, invalid_fires_band_count

if __name__ == '__main__':
    fire_id = 'ID4308111519720240903'

    layer_data = get_layer_range_data()

    # True if only one fire needs to be validated (fire_id); False if all stored fires need to be validated
    single_fire = False

    # To validate one fire specified here
    if single_fire:
        valid_data, invalid_layers, invalid_layers_spatial, missing_layers, invalid_layers_band_count = validate_one_fire_data(fire_id, layer_data)
        for layer in invalid_layers:
            print('Invalid layer (data):', layer)
        for layer in invalid_layers_spatial:
            print('Invalid layer (spatial):', layer)
        for layer in invalid_layers_band_count:
            print('Invalid layer (band count):', layer)
        print(missing_layers)
        print(f'Fire {fire_id} has {"" if valid_data else "in"}valid data')
    else:
        invalid_fires, fires_with_missing_data, invalid_fires_band_count = validate_existing_fires(layer_data_dict=layer_data, verbose=True)
        print('Fires with invalid data:')
        for fire in invalid_fires:
            print(fire)
        print('Fires with missing data:')
        for fire, missing_layers in fires_with_missing_data:
            print(fire, missing_layers)
        print('Fires with invalid band count:')
        for fire, invalid_layers_band_count in invalid_fires_band_count:
            print(fire, invalid_layers_band_count)
        
