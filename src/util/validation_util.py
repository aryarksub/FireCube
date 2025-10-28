import pandas as pd
import numpy as np
import os
import rasterio

output_cubes_dir = os.path.join('output', 'cubes')

def get_layer_range_data():
    layer_range_df = pd.read_csv('inputData\layer_ranges.csv')
    key = 'layer'
    return layer_range_df.set_index(key).to_dict(orient='index')

def find_landfire_layer_file(abbrev, lf_dir):
    for file in os.listdir(lf_dir):
        if abbrev in file and os.path.isfile(os.path.join(lf_dir, file)):
            return file
    return None

def validate_one_fire_data(fid, layer_data_dict):
    fire_folder = os.path.join('output', 'cubes', fid)
    if not os.path.exists(fire_folder):
        raise NotADirectoryError(f'No folder for fire {fid} exists')
    
    invalid_layers = []
    evt_found = False
    
    for layer, layer_data in layer_data_dict.items():
        if 'evt' in layer and evt_found:
            continue

        category = layer_data['category']

        min_val, max_val, null_val = layer_data['min'], layer_data['max'], layer_data['nan']
        
        layer_folder = os.path.join(fire_folder, category)
        if not os.path.exists(layer_folder):
            raise NotADirectoryError(f'No folder for layer category {category} and fire {fid} exists')
        
        # LANDFIRE layers are named differently than just "{layer}.tif"
        if category == 'landfire':
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
        
        with rasterio.open(layer_tif) as src:
            for band in range(src.count):
                data = src.read(band + 1)
                if null_val:
                    valid_data_mask = ((data >= min_val) & (data <= max_val)) | (data == null_val)
                else:
                    valid_data_mask = (data >= min_val) & (data <= max_val)
                
                if not np.all(valid_data_mask):
                    invalid_indices = np.where(~valid_data_mask)
                    invalid_locations = list(zip(invalid_indices[0], invalid_indices[1]))
                    invalid_layers.append((layer, band+1, invalid_locations))
                    break
        
    return len(invalid_layers) == 0, invalid_layers


def validate_existing_fires(layer_data_dict=None, verbose=False):
    layer_data_dict = layer_data_dict if layer_data_dict is not None else get_layer_range_data()

    invalid_fires = []
    processed_count = 0

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

