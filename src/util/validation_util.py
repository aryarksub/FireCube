import pandas as pd
import numpy as np
import os
import rasterio

output_cubes_dir = os.path.join('output', 'cubes')

def get_layer_range_data():
    """
    Get dictionary mapping layers to important information. Keys are layer names. Values are dictionaries that have
    keys "category", "min", "max", "nan". These represent the category of the layer, min/max values, and a value
    that is treated as a null value for the layer.

    Returns:
        dict: Dictionary of layer information.
    """
    layer_range_df = pd.read_csv('inputData\layer_ranges.csv')
    key = 'layer'
    return layer_range_df.set_index(key).to_dict(orient='index')

def find_landfire_layer_file(abbrev, lf_dir):
    """
    Find the layer file corresponding to the given abbreviation in the given directory of LANDFIRE data.

    Args:
        abbrev (str): Abbreviation of the LANDFIRE layer to search for.
        lf_dir (str): Path of the LANDFIRE directory in which we should search.

    Returns:
        str: Name of the corresponding LANDFIRE layer file.
    """
    for file in os.listdir(lf_dir):
        if abbrev in file and os.path.isfile(os.path.join(lf_dir, file)):
            return file
    return None

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
        FileNotFoundError: Occurs when there is no TIF file for a certain layer.

    Returns:
        tuple: Tuple of a boolean and a list. The boolean takes a True value if the fire data is valid; False otherwise.
         The list contains tuples of the form (layer, band, indices). These store the invalid layers, the time band at
         which the data is invalid, and the indices in the 2D TIF data array that are invalid.
    """

    # If no layer data dictionary is specified, then obtain it with get_layer_range_data()
    layer_data_dict = layer_data_dict if layer_data_dict is not None else get_layer_range_data()

    # Ensure there is an output folder for the given FID
    fire_folder = os.path.join('output', 'cubes', fid)
    if not os.path.exists(fire_folder):
        raise NotADirectoryError(f'No folder for fire {fid} exists')
    
    invalid_layers = []
    evt_found = False
    
    for layer, layer_data in layer_data_dict.items():
        # Since there are multiple potential EVT layers (LANDFIRE), if we have already found data for one, skip the others
        if 'evt' in layer and evt_found:
            continue

        category = layer_data['category']
        min_val, max_val, null_val = layer_data['min'], layer_data['max'], layer_data['nan']
        
        # Ensure there is a folder for the layer's category
        layer_folder = os.path.join(fire_folder, category)
        if not os.path.exists(layer_folder):
            raise NotADirectoryError(f'No folder for layer category {category} and fire {fid} exists')
        
        # LANDFIRE layers are named differently than just "{layer}.tif"
        if category == 'landfire':
            # Get the TIF file for the LANDFIRE layer
            lf_layer_file = find_landfire_layer_file(layer, layer_folder)

            if lf_layer_file is None:
                # Fuel behavior model layers (f13/f40) can also be named (fbfm13/fbfm40)
                if layer == 'f13' or layer == 'f40':
                    lf_layer_file = find_landfire_layer_file(layer.replace('f', 'fbfm'), layer_folder)
                    layer_tif = os.path.join(layer_folder, lf_layer_file)
                # If layer files for 105/200evt cannot be found, skip since layer may be 220evt
                elif layer == '105evt' or layer == '200evt':
                    continue
                # Since 220evt is the last possible EVT layer, if it cannot be found, then the layer does not exist
                elif layer == '220evt':
                    layer_tif = None
                else:
                    layer_tif = None
            else:
                layer_tif = os.path.join(layer_folder, lf_layer_file)
        else:
            layer_tif = os.path.join(layer_folder, f'{layer}.tif')
        
        if layer_tif is None or not os.path.exists(layer_tif):
            raise FileNotFoundError(f'No file for layer {layer} and fire {fid} exists')
        
        if 'evt' in layer:
            evt_found = True
        
        # Data validation for each band in TIF
        with rasterio.open(layer_tif) as src:
            for band in range(src.count):
                # Band counting is 1-indexed
                data = src.read(band + 1)
                # Allow null value if specified
                if null_val:
                    valid_data_mask = ((data >= min_val) & (data <= max_val)) | (data == null_val)
                else:
                    valid_data_mask = (data >= min_val) & (data <= max_val)
                
                # Store the first occurrence of invalid data
                if not np.all(valid_data_mask):
                    invalid_indices = np.where(~valid_data_mask)
                    invalid_locations = list(zip(invalid_indices[0], invalid_indices[1]))
                    invalid_layers.append((layer, band+1, invalid_locations))
                    break
        
    return len(invalid_layers) == 0, invalid_layers


def validate_existing_fires(layer_data_dict=None, verbose=False):
    """
    Perform data validation for all fires that are stored in the output directory. This is done by simply iterating
    through the directory and running validate_one_fire_data().

    Args:
        layer_data_dict (dict, optional): Dictionary of layer information (see get_layer_range_data() for dictionary
         format). Defaults to None.
        verbose (bool, optional): True if descriptive messages should be printed; False otherwise. Defaults to False.

    Returns:
        list: List of tuples of the form (fid, invalid_layers), where fid corresponds to a fire event ID with invalid
         data and invalid_layers is the list returned by validate_one_fire_data().
    """
    # If no layer data dictionary is specified, then obtain it with get_layer_range_data()
    layer_data_dict = layer_data_dict if layer_data_dict is not None else get_layer_range_data()

    invalid_fires = []
    processed_count = 0

    # Note that folders are named with the fire event ID (so folder = FID)
    for folder in os.listdir(output_cubes_dir):
        if os.path.isdir(os.path.join(output_cubes_dir, folder)):
            valid_data, invalid_layers = validate_one_fire_data(folder, layer_data_dict)
            if not valid_data:
                invalid_fires.append((folder, invalid_layers))

                if verbose:
                    print(f'Fire {folder} has invalid data')
            
            processed_count += 1

            if verbose and processed_count % 10 == 0:
                print(f'Validated {processed_count} fires')
    
    return invalid_fires

if __name__ == '__main__':
    fire_id = 'AL3105608678420130203'

    layer_data = get_layer_range_data()

    # True if only one fire needs to be validated (fire_id); False if all stored fires need to be validated
    single_fire = False

    # To validate one fire specified here
    if single_fire:
        valid_data, invalid_layers = validate_one_fire_data(fire_id, layer_data)
        for layer in invalid_layers:
            print(layer)
        print(f'Fire {fire_id} has {"" if valid_data else "in"}valid data')
    else:
        invalid_fires = validate_existing_fires(layer_data_dict=layer_data, verbose=True)
        for fire in invalid_fires:
            print(fire)

