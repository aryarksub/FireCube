import os
import argparse
import zipfile
import pandas as pd
import rasterio
from rasterio.windows import from_bounds, Window
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.transform import from_origin
import numpy as np
import time
import tempfile
import shutil
import math

def get_zip_path_from_LF_version_and_variable(landfire_dir, LF_year, variable, region):
    # Map version year to code (a unique identifier for each LF version filename)
    LF_version_year_to_code = {
        2010 : '_120_',
        2012 : '_130_',
        2014 : '_140_',
        2016 : '2016_',
        2020 : '_220_',
        2022 : '2022_',
        2023 : '2023_',
        2024 : '2024_'
    }

    try:
        LF_code = LF_version_year_to_code[LF_year]
    except:
        raise ValueError(f"Provided LF version {LF_year} is not valid")

    # Find zipfile with code and variable name
    zip_paths = os.listdir(landfire_dir)
    candidates = [z for z in zip_paths if '.zip' in z and LF_code in z and variable in z and (region == 'AK') == ('AK' in z)]

    if len(candidates) == 1:
        zip_path = os.path.join(landfire_dir, candidates[0])
        return zip_path
    elif len(candidates) == 0:
        raise ValueError(f"Provided LF Version {LF_year} does not have variable {variable} in {landfire_dir}")
    else:
        print(candidates)
        raise ValueError(f"Multiple valid zipfiles for {LF_year} and variable {variable}")

def process_landfire_variable(landfire_dir, landfire_version, variable, csv_path, output_dir, bilinear=False, region='CONUS'):
    os.makedirs(output_dir, exist_ok=True)
    
    # Retrieve fire list corresponding to LF version
    df = pd.read_csv(csv_path)
    df = df[df['LF_version'] == landfire_version]
    print(f"Processing {len(df)} fires...")
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Retrieve TIF corresponding to LF version and variable
        load_tif_start = time.perf_counter()
        zip_path = get_zip_path_from_LF_version_and_variable(landfire_dir, landfire_version, variable, region)
        print(f"Extracting tif into temporary directory...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            tif_files = [f for f in zip_ref.namelist() if f.endswith('.tif')]
            if not tif_files:
                raise ValueError(f"No {tif_name} found in the provided zip archive.")
            tif_name = tif_files[0]

            zip_ref.extract(tif_name, tmp_dir)
            full_tif_path = os.path.join(tmp_dir, tif_name)
    
        # Load the entire TIF into memory
        print(f"Loading {tif_name} into memory...")
        with rasterio.open(full_tif_path) as src:
            data = src.read() 
            profile = src.profile
            transform = src.transform
            max_rows = data.shape[1]
            max_cols = data.shape[2]
        load_tif_stop = time.perf_counter()
        runtime_min = (load_tif_stop - load_tif_start)/60
        print(f"TIF loaded successfully. Runtime = {runtime_min:.2f} minutes")
    
        # Set up tif storage directory
        tif_save_path = os.path.join(tmp_dir, variable)
        os.mkdir(tif_save_path)

        # Iterate through the fires and slice the array
        if bilinear:
            print("Generating tifs with bilinear resampling...")
        else:
            print("Generating tifs with nearest-neighbor resampling...") 
        fires_start = time.perf_counter()
        for index, row in df.iterrows():
            fire_id = row['Event_ID']
            minx, miny, maxx, maxy = row['minx'], row['miny'], row['maxx'], row['maxy']

            if region == 'AK':
                minx_5070, miny_5070, maxx_5070, maxy_5070 = minx, miny, maxx, maxy
                minx_3338, miny_3338, maxx_3338, maxy_3338 = transform_bounds(
                    "EPSG:5070", "EPSG:3338", 
                    minx_5070, miny_5070, maxx_5070, maxy_5070
                )
                # Add a 5,000 meter buffer to the source bounds to account for rotation
                buffer_m = 5000
                minx_3338 -= buffer_m
                miny_3338 -= buffer_m
                maxx_3338 += buffer_m
                maxy_3338 += buffer_m
                window = from_bounds(minx_3338, miny_3338, maxx_3338, maxy_3338, transform)
            else:
                window = from_bounds(minx, miny, maxx, maxy, transform)
                window = Window(window.col_off, window.row_off, int(round(window.width)), int(round(window.height)))
            
            h = int(window.height)
            w = int(window.width)
            
            if bilinear:
                # Compute bilinear weights based on window offset
                r_exact, c_exact = window.row_off, window.col_off
                r0, c0 = int(math.floor(r_exact)), int(math.floor(c_exact))
                
                dy = r_exact - r0
                dx = c_exact - c0
                
                w_tl = (1 - dx) * (1 - dy)  # Top-Left weight
                w_tr = dx * (1 - dy)        # Top-Right weight
                w_bl = (1 - dx) * dy        # Bottom-Left weight
                w_br = dx * dy              # Bottom-Right weight

                r1 = r0 + 1
                c1 = c0 + 1
                
                # Handle potential out-of-bounds errors
                if min(r0, c0) < 0 or (r1 + h) > max_rows or (c1 + w) > max_cols:
                    print(f"{fire_id} exceeded tif bounds. Skipping")
                    continue
                
                # Extract the four slices
                slice_tl = data[:, r0 : r0 + h, c0 : c0 + w]
                slice_tr = data[:, r0 : r0 + h, c1 : c1 + w]
                slice_bl = data[:, r1 : r1 + h, c0 : c0 + w]
                slice_br = data[:, r1 : r1 + h, c1 : c1 + w]
                
                # Apply weights
                sliced_data = (
                    w_tl * slice_tl + 
                    w_tr * slice_tr + 
                    w_bl * slice_bl + 
                    w_br * slice_br
                ).astype(data.dtype)
            else:
                # Mimic nearest neighbor resampling by snapping to the nearest gridcell
                r0 = round(window.row_off)
                c0 = round(window.col_off)
                
                # Handle potential out-of-bounds errors
                if min(r0, c0) < 0 or (r0 + h) > max_rows or (c0 + w) > max_cols:
                    print(f"{fire_id} exceeded tif bounds. Skipping")
                    continue
                
                # Slice the array directly from memory
                sliced_data = data[:, r0 : r0 + h, c0 : c0 + w]

            if region == 'AK':
                slice_transform = rasterio.windows.transform(window, transform)

                dst_width = int(round((maxx_5070 - minx_5070) / 30))
                dst_height = int(round((maxy_5070 - miny_5070) / 30))
                dst_transform = from_origin(minx_5070, maxy_5070, 30, 30)
                
                slice_5070 = np.zeros((data.shape[0], dst_height, dst_width), dtype=data.dtype)

                resampling = Resampling.bilinear if bilinear else Resampling.nearest
                
                reproject(
                    source=sliced_data,
                    destination=slice_5070,
                    src_transform=slice_transform,
                    src_crs="EPSG:3338",
                    dst_transform=dst_transform,
                    dst_crs="EPSG:5070",
                    resampling=resampling
                )

                sliced_data = slice_5070

                # Update metadata
                out_profile = profile.copy()
                out_profile.update({
                    "crs": "EPSG:5070",
                    "transform": dst_transform,
                    "width": dst_width,
                    "height": dst_height,
                    "compress": "lzw"
                })
                
            else:
                # Update metadata
                new_transform = rasterio.windows.transform(window, transform)
                out_profile = profile.copy()
                out_profile.update({
                    "height": sliced_data.shape[1],
                    "width": sliced_data.shape[2],
                    "transform": new_transform,
                    "compress": "lzw"
                })
            
            # Save to temporary directory
            out_path = os.path.join(tif_save_path, f"{fire_id}.tif")
            with rasterio.open(out_path, "w", **out_profile) as dest:
                dest.write(sliced_data)
            print(fire_id)
        
        fires_stop = time.perf_counter()
        print(f"Successfully processed {len(df)} fires. Runtime = {fires_stop - fires_start:.2f} seconds")

        # Zip tifs and save
        destination_zip = os.path.join(output_dir, f'LF{landfire_version}_{region}_{variable}')
        shutil.make_archive(destination_zip, 'zip', tif_save_path)
        print(f"Saved tifs to {destination_zip}.zip")

def main():
    parser = argparse.ArgumentParser(
        description="Retrieve tifs for a provided LF version and variable.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--landfire-path',
        dest='lf_path',
        required=True,
        help='Path to directory containing LANDFIRE archives.'
    )
    parser.add_argument(
        '--landfire-version',
        dest='version',
        required=True,
        help='LF Version to retrieve'
    )
    parser.add_argument(
        '--landfire-var',
        dest='var',
        required=True,
        help='LF variable to retrieve'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Path to the csv containing FIDs, bounds, and LF versions.'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the destination directory to save the zipped tifs.'
    )
    parser.add_argument(
        '--bilinear',
        required=False,
        default="False",
        help='Whether to use bilinear resampling.'
    )
    parser.add_argument(
        '--region',
        required=False,
        default="CONUS",
        help='CONUS or AK.'
    )

    args = parser.parse_args()
    
    process_landfire_variable(args.lf_path,
                              int(args.version),
                              args.var,
                              args.input,
                              args.output,
                              args.bilinear.lower() == 'true',
                              args.region)

if __name__ == "__main__":
    main()