import os
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import geopandas as gpd

def get_datetime_cutoffs_for_LF_version(LFv):
    #define the earliest possible date for a LF version to be the optimal version
    cutoffs = [
        (datetime(2010, 12, 31), 2010),
        (datetime(2012, 12, 31), 2012),
        (datetime(2014, 12, 31), 2014),
        (datetime(2016, 12, 31), 2016),
        (datetime(2020, 12, 31), 2020),
        (datetime(2022, 12, 31), 2022),
        (datetime(2023, 10, 1), 2023),
        (datetime(2024, 10, 1), 2024),
    ]
    
    lower = None
    upper = datetime(2025, 1, 1)

    # seek the lower and upper bound for the date range of a given LF version
    for cutoff, v in cutoffs:
        if LFv < v:
            upper = cutoff
            break
        else:
            lower = cutoff

    if lower is None:
        raise ValueError("LANDFIRE Version is too early. Must be 2010 or later")
    
    return lower, upper

def get_fires_for_LF_versions(firelist_csv, LF_versions, region='CONUS', subset=None):
    # open full PyroStack firelist
    df = pd.read_csv(firelist_csv)
    df['tst'] = pd.to_datetime(df['tst'])

    # filter to only include fires such that the provided LF version is the most recent LF version
    dfs_by_version = []
    for LFv in LF_versions:
        lower_cutoff, upper_cutoff = get_datetime_cutoffs_for_LF_version(LFv)
        df_sub = df[(df['tst'] >= lower_cutoff) & (df['tst'] < upper_cutoff)].reset_index(drop=True)
        dfs_by_version.append(df_sub)

    df_all_versions = pd.concat(
        [df.assign(LF_version=label) for df, label in zip(dfs_by_version, LF_versions)], 
        ignore_index=True
    )

    # filter by region of US
    if region == 'CONUS':
        df_all_versions = df_all_versions[~df_all_versions['Event_ID'].str.contains('AK')]
    elif region == 'AK':
        df_all_versions = df_all_versions[df_all_versions['Event_ID'].str.contains('AK')]
    else:
        raise ValueError(f"Region {region} is not supported. Must be CONUS or AK")

    # filter using unique FID substring
    if subset is None:
        return df_all_versions[['Year', 'Event_ID', 'LF_version']]
    else:
        df_subset = df_all_versions[np.isin(df_all_versions['Event_ID'], subset)]
        if len(df_subset.index) == 0:
            raise ValueError("Provided FID subset does not match any fires from the specified LF version")
        return df_subset[['Year', 'Event_ID', 'LF_version']]

def pad_bounds_to_resolution_multiple(bounds, coarsest_res=9000):
    """
    Pad the given bounds to be a multiple of the given resolution (typically the coarsest resolution amongst a set
    of TIFs). If the required padding is less than 10% of the given resolution, add an additional full resolution buffer.

    Args:
        bounds (list): List of bounding box coordinates in the same CRS as the input TIF. There should be
         four coordinates: x_min, y_min, x_max, y_max.
        coarsest_res (float, optional): Coarsest resolution in meters. Defaults to 9000.

    Returns:
        list: Padded bounds.
    """
    x1,y1,x2,y2 = bounds
    center_x, center_y = (x1 + x2) / 2, (y1 + y2) / 2

    # Use coarsest tif to determine number of pixels required to cover bounds
    width_m = x2 - x1   # width in meters
    height_m = y2 - y1  # height in meters

    # Padding needed in each direction to get dimensions that are a multiple of coarsest_res
    pad_x = coarsest_res - (width_m % coarsest_res) if width_m % coarsest_res != 0 else 0
    pad_y = coarsest_res - (height_m % coarsest_res) if height_m % coarsest_res != 0 else 0

    # If the previously derived padding is smaller than 10% of the coarsest resolution, add a 100% coarsest 
    # resolution buffer in that direction.
    if pad_x < coarsest_res / 10:   
        pad_x += coarsest_res
    if pad_y < coarsest_res / 10:
        pad_y += coarsest_res

    # Dimensions for output in meters (multiple of coarsest_res)
    final_width_m = width_m + pad_x
    final_height_m = height_m + pad_y

    # Bounding box for output
    half_width = final_width_m / 2
    half_height = final_height_m / 2
    final_x1, final_x2 = center_x - half_width, center_x + half_width
    final_y1, final_y2 = center_y - half_height, center_y + half_height

    return [final_x1, final_y1, final_x2, final_y2]

def get_bounds_for_fires(args, df_fires):
    # generate a list of geopackage paths
    years = df_fires['Year'].unique()

    gpkgs = []
    for y in years:
        year_path = os.path.join(args.pyrostack_path, 'inputData', 'Full_FEDS', str(y))
        gpkgs.extend(os.listdir(year_path))

    # retrieve the reprojected, padded bounds for each fire
    bounds = []
    for i, row in df_fires.iterrows():
        year = row['Year']
        fid = row['Event_ID']

        #Retrieve gpkg path
        matches = [m for m in gpkgs if fid in m]
    
        if len(matches) == 0:
            raise RuntimeError("No FID match for " + fid + " in " + args.pyrostack_path)
        elif len(matches) > 1:
            raise RuntimeError("Multiple FID matches for " + fid + " in " + args.pyrostack_path)
        
        gpkg_path = os.path.join(args.pyrostack_path, 'inputData', 'Full_FEDS', str(year), matches[0])

        #Reproject fire to EPSG:5070 and get total bounds
        gdf_fperim_rd = gpd.read_file(gpkg_path, layer='perimeter')
        gdf_fperim_5070 = gdf_fperim_rd.to_crs('EPSG:5070') 
        bounds_5070 = gdf_fperim_5070.total_bounds

        padded_bounds_5070 = pad_bounds_to_resolution_multiple(bounds_5070)

        bounds.append(padded_bounds_5070)

    df_fires[['minx', 'miny', 'maxx', 'maxy']] = pd.DataFrame(bounds, index=df_fires.index)

    return df_fires.drop(columns='Year')

def main():
    parser = argparse.ArgumentParser(
        description="Retrieve tifs for a provided LF version and variable.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--pyrostack-path',
        required=True,
        help='Path to PyroStack base directory'
    )
    parser.add_argument(
        '--landfire-version',
        dest='version',
        required=True,
        help='LF Version to retrieve'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to save the csv containing FIDs, bounds, and LF versions.'
    )
    parser.add_argument(
        '--region',
        required=False,
        default="CONUS",
        help='CONUS or AK.'
    )
    parser.add_argument(
        '--subset',
        required=False,
        default=None,
        help='Substring unique to the subset of Fire IDs that you seek to extract (e.g., `CA` to only retrieve California fires)'
    )

    args = parser.parse_args()

    firelist_csv = os.path.join(args.pyrostack_path, 'inputData', 'fireslist2012-2024_withtype_no_fline_hawaii.csv')
    
    df = get_fires_for_LF_versions(firelist_csv, [int(args.version)], args.region, args.subset)
    df_bounds = get_bounds_for_fires(args, df)
    df_bounds.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()