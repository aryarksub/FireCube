"""
run_pyretechnics.py

Purpose:
    For an individual fire interval:
        Retrieves the environmental fields for PyroStack and pre-processes them to become Pyretechnics inputs
        Runs a single Pyretechnics simulation
        Plots an example output figure

Usage:
    python src/run_pyretechnics.py 
        --fire-id FIRE_ID
        --start-day START_DAY
        --duration DURATION
        --pyrostack-path /path/to/pyrostack_dir
        --output-path /path/to/output_dir
    
Arguments:
    --fire-id: Fire ID from FEDS database.
    --start-day: Day since first FEDS retrieval ("duration" column of FEDS gpkg).
    --duration: Length of interval in days.
    --pyrostack-path: Path to the PyroStack directory.
    --output-path: Path to the output directory to save figure.
"""

# pre-process / simulation requirements
import os
import time
import argparse
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
import pyretechnics.eulerian_level_set as els
import pyretechnics.load_landfire as lf
from pyretechnics.space_time_cube import SpaceTimeCube

# plotting requirements
import cartopy.crs as ccrs
from rasterio.plot import show
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import matplotlib.axes as maxes
from mpl_toolkits.axes_grid1 import make_axes_locatable

def parse_args(args=None):
    def str2float(v):
        return float(v)
    
    parser = argparse.ArgumentParser(
        description="Example Pyretechnics Simulation Script",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--fire-id',
        dest='fid',
        required=True,
        help='Fire ID'
    )
    parser.add_argument(
        '--pyrostack-path',
        required=True,
        help='Path to the PyroStack base directory.'
    )
    parser.add_argument(
        '--start-day',
        type=str2float,
        required=False,
        default=None,
        help='Start day of the desired interval to plot.'
    )
    parser.add_argument(
        '--duration',
        type=str2float,
        required=False,
        default=None,
        help='Duration in days of the desired interval to plot.'
    )
    parser.add_argument(
        '--output-path',
        required=True,
        default=True,
        help='Where to save output figure'
    )

    args = parser.parse_args(args)
    return args

def get_gpkg_path(args):
    year = args.fid[-8:-4]
    pyrostack_dir = os.path.join(args.pyrostack_path, 'inputData', 'Full_FEDS', year)
    
    #find gpkg file that contians FID in PyroStack inputs
    if os.path.isdir(pyrostack_dir):
        matches = []
        for fname in os.listdir(pyrostack_dir):
            if args.fid in fname:
                matches.append(os.path.join(pyrostack_dir, fname))
    else:
        raise RuntimeError("Directory does not exist: " + pyrostack_dir)

    if len(matches) == 0:
        raise RuntimeError("No FID match for " + args.fid + " in " + pyrostack_dir)
    elif len(matches) > 1:
        raise RuntimeError("Multiple FID matches for " + args.fid + " in " + pyrostack_dir)
    
    return matches[0]

def get_pyrostack_indices(args, gpkg_path):
    gdf = gpd.read_file(gpkg_path, layer='perimeter')
    gdf["duration"] = gdf["duration"].apply(pd.to_numeric, errors='coerce')
    
    #1) map start day to a timestamp
    current_day = args.start_day
    
    diff = np.abs(gdf.duration.values - current_day)
    idx = np.argmin(diff)
        
    if diff[idx] > 1e-5:
        raise ValueError(f"Day {current_day} not found in geopackage!")
    
    t_start = pd.to_datetime(gdf.loc[idx, 't']) + pd.Timedelta(hours=1, minutes=30)
    
    #2) determine the first timestamp of the geopackage entries
    t0_gpkg = t_start - pd.Timedelta(days=current_day)
    
    #3) determine the timestamp of the first band of tifs
    initial_geometry = gpd.GeoSeries(gdf.iloc[idx].geometry, crs=gdf.crs)
    projected_geometry = initial_geometry.to_crs("EPSG:4326").iloc[0]
    longitude = projected_geometry.centroid.x
    
    t_start_utc = t_start - pd.to_timedelta(round(longitude / 15), unit="hours")
    t0_tif = (t0_gpkg - pd.to_timedelta(round(longitude / 15), unit="hours") - pd.Timedelta(days=1)).normalize()
    
    #4) use difference between timestamps to find the start and end index
    timesteps = int(args.duration * 24)
    start_index = int((t_start_utc - t0_tif).total_seconds() // 3600)
    end_index = start_index + timesteps + 1

    return start_index, end_index, timesteps

def get_landfire_tif_path(fire_path, raster):
    landfire_tifs = os.listdir(os.path.join(fire_path, "veg_fm_topo"))
    if raster == "fbfm40.tif":
        tif = [f for f in landfire_tifs if ("f40" in f or "m40" in f) and ".tif" in f]
    elif raster == "slp.tif":
        tif = [f for f in landfire_tifs if "slp" in f and ".tif" in f]
    elif raster == "asp.tif":
        tif = [f for f in landfire_tifs if "asp" in f and ".tif" in f]
    
    if len(tif) == 1:
        return os.path.join(fire_path, "veg_fm_topo", tif[0])
    else:
        raise RuntimeError(f"Could not find landfire tif: {raster}")

def rename_keys(the_dict, keys_dict):
    for (f, t) in keys_dict.items():
        the_dict[t] = the_dict.pop(f)

def load_partial_raster(file_path, band_slice, dtype=None, cube_shape_divisors=(1,1,1),
                        resampling_policy="nearest_match",
                        resampling_method=Resampling.nearest):
    """
    Adapted from: https://github.com/pyregence/pyretechnics
    resampling_policy: "always_upsample" or "nearest_match"
    resampling_method: any rasterio.enums.Resampling method
    """
    with rasterio.open(file_path, "r") as input_raster:
        metadata  = lf.raster_metadata(input_raster)

        if band_slice:
            start, stop = band_slice
            target_indexes = list(range(start + 1, stop + 1))
            bands = len(target_indexes)
        else:
            target_indexes = None # None implies "read all bands"
            bands = metadata["bands"]
        
        rows      = metadata["rows"]
        cols      = metadata["cols"]
        (b, r, c) = cube_shape_divisors
        new_bands = lf.maybe_resample_resolution(bands, b, resampling_policy)
        new_rows  = lf.maybe_resample_resolution(rows, r, resampling_policy)
        new_cols  = lf.maybe_resample_resolution(cols, c, resampling_policy)
        if new_bands == bands and new_rows == rows and new_cols == cols:
            return {
                "array"   : input_raster.read(indexes=target_indexes, out_dtype=dtype),
                "metadata": metadata,
            }
        else:
            metadata["bands"]     = new_bands
            metadata["rows"]      = new_rows
            metadata["cols"]      = new_cols
            metadata["transform"] = (input_raster.transform
                                     * input_raster.transform.scale(
                                         cols / new_cols,
                                         rows / new_rows,
                                     ))
            array = input_raster.read(
                indexes=target_indexes,
                out_dtype=dtype,
                out_shape=(new_bands, new_rows, new_cols),
                resampling=resampling_method,
            )
            return {
                "array"   : array,
                "metadata": metadata,
            }

def feds_ignition_field(args):
    gpkg_path = get_gpkg_path(args)

    perimeters = gpd.read_file(gpkg_path, layer="perimeter")
    perimeters["duration"] = perimeters["duration"].apply(pd.to_numeric, errors='coerce')
    flines = gpd.read_file(gpkg_path, layer="fireline")
    
    perim_init = perimeters[perimeters["duration"] == args.start_day]
    perim_final = perimeters[perimeters["duration"] == args.start_day + args.duration]
    fline_init = flines[flines["t"] == perim_init["t"].iloc[0]]
    
    perimeter_reprojected = gpd.GeoDataFrame(geometry=perim_init["geometry"],
                                             crs="EPSG:4326").to_crs("EPSG:5070")
    perimeter_final_reprojected = gpd.GeoDataFrame(geometry=perim_final["geometry"],
                                                   crs="EPSG:4326").to_crs("EPSG:5070")
    fline_reprojected = gpd.GeoDataFrame(geometry=fline_init["geometry"],
                                         crs="EPSG:4326").to_crs("EPSG:5070")

    ignition_hour = fline_init["t"].iloc[0].hour
    fuel_tif = os.path.join(args.pyrostack_path, 'output', 'cubes', args.fid, 'fuel_structure', 'cc.tif')

    with rasterio.open(fuel_tif) as src:
        perimeter_shapes = ((geom, 1) for geom in perimeter_reprojected.geometry if geom is not None)
        perimeter_final_shapes = ((geom, 1) for geom in perimeter_final_reprojected.geometry if geom is not None)
        fline_shapes = ((geom, 1) for geom in fline_reprojected.geometry if geom is not None)

        cold_burn_raster = rasterize(
            perimeter_shapes,
            out_shape=src.shape,
            transform=src.transform,
            fill=0,
            all_touched=True,
            dtype='uint8'
        )
        active_burn_raster = rasterize(
            fline_shapes,
            out_shape=src.shape,
            transform=src.transform,
            fill=0,
            all_touched=True,
            dtype='uint8'
        )
        
        pre_burned = np.where((cold_burn_raster == 1) & (active_burn_raster == 0), 1, 0)
        ignition = active_burn_raster.astype('uint8')
        
        return pre_burned, ignition, ignition_hour

def load_and_transform_rasters(args, band_slice):
    fire_path = os.path.join(args.pyrostack_path, "output", "cubes", args.fid)

    #load rasters
    pyregence_rasters = {
        raster: lf.load_raster(os.path.join(fire_path, "fuel_structure", raster), dtype="float32")
        for raster in [
            "cc.tif",
            "ch.tif",
            "cbh.tif",
            "cbd.tif",
        ]
    }
    landfire_rasters = {
        raster: lf.load_raster(get_landfire_tif_path(fire_path, raster), dtype="float32")
        for raster in [
            "fbfm40.tif",
            "slp.tif",
            "asp.tif",
        ]
    }
    non_landfire_rasters = {
        raster: load_partial_raster(os.path.join(fire_path, "high_res_climate", raster),
                                                    band_slice,
                                                    dtype="float32")
        for raster in [
            "ws.tif",
            "wd.tif",
            "m1.tif",
            "m10.tif",
            "m100.tif",
            "lh.tif",
            "lw.tif",
        ]
    }
    temp_raster = {"t2m.tif" : load_partial_raster(os.path.join(fire_path, "low_res_climate", "t2m.tif"), band_slice, dtype="float32")}
    rasters_dict = pyregence_rasters | landfire_rasters | non_landfire_rasters | temp_raster
    rename_keys(
        rasters_dict,
        {
            "slp.tif": "slope",  # magnitude
            "asp.tif": "aspect",  # direction
            "fbfm40.tif": "fuel_model",  # surface vegetation
            "cc.tif": "canopy_cover",  # how close trees are 0% (no trees) to 100% (all trees)
            "ch.tif": "canopy_height",  # how tall the trees are
            "cbh.tif": "canopy_base_height",  # how high to the lowest branches
            "cbd.tif": "canopy_bulk_density",  # amount of material in tree branches
            "ws.tif": "wind_speed_10m",  # magnitude
            "wd.tif": "upwind_direction",  # direction
            "m1.tif": "fuel_moisture_dead_1hr",  # thinner
            "m10.tif": "fuel_moisture_dead_10hr",  # 0.5 cm to 2cm
            "m100.tif": "fuel_moisture_dead_100hr",  # 2cm to 6.5cm
            "lh.tif": "fuel_moisture_live_herbaceous",  # leaves, grasses etc.
            "lw.tif": "fuel_moisture_live_woody",  # branches, trunks etc.
            "t2m.tif": "temperature", #2 meter temperature
        },
    )

    # SLOPE
    slope_array = rasters_dict["slope"]["array"]
    np.tan(np.radians(slope_array, out=slope_array), out=slope_array) #convert degrees to ratio

    # ASPECT
    aspect_array = rasters_dict["aspect"]["array"]
    np.subtract(180, aspect_array, out=aspect_array)
    np.mod(aspect_array, 360, out=aspect_array)
    aspect_array[np.isnan(aspect_array)] = -1

    # FUEL_MODEL
    fuel_model_array = rasters_dict["fuel_model"]["array"]
    fuel_model_array[np.isnan(fuel_model_array)] = -9999

    # CANOPY_COVER
    canopy_cover_array = rasters_dict["canopy_cover"]["array"]
    np.multiply(canopy_cover_array, 0.01, out=canopy_cover_array)
    canopy_cover_array[np.isnan(canopy_cover_array)] = -9999

    # CANOPY_HEIGHT
    canopy_height_array = rasters_dict["canopy_height"]["array"]
    np.multiply(canopy_height_array, 0.1, out=canopy_height_array)
    canopy_height_array[np.isnan(canopy_height_array)] = -9999

    # CANOPY_BASE_HEIGHT
    canopy_base_height_array = rasters_dict["canopy_base_height"]["array"]
    np.multiply(canopy_base_height_array, 0.1, out=canopy_base_height_array)
    canopy_base_height_array[np.isnan(canopy_base_height_array)] = -9999

    # CANOPY_BULK_DENSITY
    canopy_bulk_density_array = rasters_dict["canopy_bulk_density"]["array"]
    np.multiply(canopy_bulk_density_array, 0.01, out=canopy_bulk_density_array)
    canopy_bulk_density_array[np.isnan(canopy_bulk_density_array)] = -9999

    # WIND_SPEED
    wind_speed_10m_array = rasters_dict["wind_speed_10m"]["array"] # miles/hour at 20 ft
    np.multiply(wind_speed_10m_array, 1.609344 * 1.15, out=wind_speed_10m_array) # convert to km/h at 10 meters

    # UPWIND_DIRECTION
    upwind_direction_array = rasters_dict["upwind_direction"]["array"]
    np.subtract(180, upwind_direction_array, out=upwind_direction_array)
    np.mod(upwind_direction_array, 360, out=upwind_direction_array)

    # FUEL_MOISTURE_DEAD_1HR
    fuel_moisture_dead_1hr_array = rasters_dict["fuel_moisture_dead_1hr"]["array"]
    np.multiply(fuel_moisture_dead_1hr_array, 0.01, out=fuel_moisture_dead_1hr_array)
    fuel_moisture_dead_1hr_array[np.isnan(fuel_moisture_dead_1hr_array)] = -9999

    # FUEL_MOISTURE_DEAD_10HR
    fuel_moisture_dead_10hr_array = rasters_dict["fuel_moisture_dead_10hr"]["array"]
    np.multiply(fuel_moisture_dead_10hr_array, 0.01, out=fuel_moisture_dead_10hr_array)
    fuel_moisture_dead_10hr_array[np.isnan(fuel_moisture_dead_10hr_array)] = -9999

    # FUEL_MOISTURE_DEAD_100HR
    fuel_moisture_dead_100hr_array = rasters_dict["fuel_moisture_dead_100hr"]["array"]
    np.multiply(fuel_moisture_dead_100hr_array, 0.01, out=fuel_moisture_dead_100hr_array)
    fuel_moisture_dead_100hr_array[np.isnan(fuel_moisture_dead_100hr_array)] = -9999

    # FUEL_MOISTURE_LIVE_HERBACEOUS
    fuel_moisture_live_herbaceous_array = rasters_dict["fuel_moisture_live_herbaceous"]["array"]
    np.multiply(fuel_moisture_live_herbaceous_array, 0.01, out=fuel_moisture_live_herbaceous_array)

    # FUEL_MOISTURE_LIVE_WOODY
    fuel_moisture_live_woody_array = rasters_dict["fuel_moisture_live_woody"]["array"]
    np.multiply(fuel_moisture_live_woody_array, 0.01, out=fuel_moisture_live_woody_array)

    # TEMPERATURE
    temperature_array = rasters_dict["temperature"]["array"]
    np.subtract(temperature_array, 273.15, out=temperature_array)

    # GENERATE IGNITION FIELD
    pre_burned, ignition, ignition_hour = feds_ignition_field(args)
    # set pre-burned fuel to unburnable
    fuel_model_array[:, pre_burned == 1] = 99.0
    
    return rasters_dict, pre_burned, ignition

def make_simple_plot(args, pre_burn, ignition, toa_matrix, gpkg_path):
    # retrieve fire name
    firelist = pd.read_csv(os.path.join(args.pyrostack_path, "inputData", "fireslist2012-2024.csv"))
    firename = firelist.loc[firelist[firelist["Event_ID"] == args.fid].index[0]]["Incid_Name"]
    sup_title = f"{firename} (ID: {args.fid})"

    # retrieve tif data
    pyrostack_output_dir = os.path.join(args.pyrostack_path, "output", "cubes", args.fid)
    prefix = 'ak_' if 'AK' in args.fid else ''
    with rasterio.open(os.path.join(pyrostack_output_dir, 'veg_fm_topo', prefix + 'elev2020.tif')) as dem_src, \
         rasterio.open(os.path.join(pyrostack_output_dir, 'veg_fm_topo', prefix + 'slpd2020.tif')) as slp_src, \
         rasterio.open(os.path.join(pyrostack_output_dir, 'veg_fm_topo', prefix + 'asp2020.tif')) as asp_src:

        dem = dem_src.read(1)
        slp = slp_src.read(1)
        asp = asp_src.read(1)
        asp[np.isnan(asp)] = 0

        # compute hillshade
        slp_rad = np.radians(slp)
        asp_rad = np.radians(asp)
        
        azimuth_sun = 315
        altitude_sun = 45
        
        zenith_rad = np.radians(90 - altitude_sun)
        azimuth_rad = np.radians(azimuth_sun)
            
        hillshade = (np.cos(zenith_rad) * np.cos(slp_rad)) + \
                    (np.sin(zenith_rad) * np.sin(slp_rad) * np.cos(azimuth_rad - asp_rad))
        hillshade = np.maximum(hillshade, 0)

        #retrieve metadata
        map_proj = ccrs.epsg(dem_src.crs.to_epsg())
        transform = dem_src.transform
        b = dem_src.bounds
        extent = [b.left, b.right, b.bottom, b.top]

    # retrieve perimeter data
    perimeters = gpd.read_file(gpkg_path, layer="perimeter")
    perimeters["duration"] = perimeters["duration"].apply(pd.to_numeric, errors='coerce')
    flines = gpd.read_file(gpkg_path, layer="fireline")
    
    perims = perimeters[(perimeters["duration"] == args.start_day) | (perimeters["duration"] == args.start_day + args.duration)]
    fline_init = flines[flines["t"] == perims["t"].iloc[0]]
    perimeters_reprojected = perims.to_crs("EPSG:5070")
    fline_reprojected = gpd.GeoDataFrame(geometry=fline_init["geometry"],
                                         crs="EPSG:4326").to_crs("EPSG:5070")

    # crop to the region where the burn occurs
    rows, cols = np.where((toa_matrix > 0) | (pre_burn))
    buffer = 10
    perim_bounds = perimeters_reprojected.total_bounds
    perim_xmin, perim_ymin, perim_xmax, perim_ymax = perim_bounds
    perim_buffer = 300
    
    r0, r1 = max(0, rows.min() - buffer), min(toa_matrix.shape[0], rows.max() + buffer) 
    c0, c1 = max(0, cols.min() - buffer), min(toa_matrix.shape[1], cols.max() + buffer)
    
    x_at_c0, y_at_r0 = transform * (c0, r0)
    x_at_c1, y_at_r1 = transform * (c1, r1)
    
    left = min(x_at_c0, x_at_c1, perim_xmin - perim_buffer)
    right = max(x_at_c0, x_at_c1, perim_xmax + perim_buffer)
    bottom = min(y_at_r0, y_at_r1, perim_ymin - perim_buffer)
    top = max(y_at_r0, y_at_r1, perim_ymax + perim_buffer)

    # set up figure
    fig_width = 7
    fig_height = (fig_width * (top - bottom) / (right - left)) + 1

    current_width = right - left
    zoom_scale = min((toa_matrix.shape[1] * transform[0]) / current_width, 4)
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), subplot_kw={'projection': map_proj})

    # plot hillshade layer
    data_crs = ccrs.epsg(5070)
    ax.imshow(hillshade, cmap='gray', extent=extent, 
              origin='upper', transform=data_crs, interpolation='none', alpha=1.0)

    # plot time of arrival field
    toa_img = ax.imshow(toa_matrix, cmap='Oranges', extent=extent, 
                       origin='upper', transform=data_crs, interpolation='none', alpha=0.5)

    # plot final observed perimeter
    perimeters_reprojected.iloc[[1]].plot(
        ax=ax,
        facecolor="none",
        edgecolor="red",
        transform=data_crs,
        linewidth=1 * zoom_scale,
        zorder=3
    )

    # plot initial cold perimeter
    perimeters_reprojected.iloc[[0]].plot(
        ax=ax,
        facecolor="none",
        edgecolor="blue",
        transform=data_crs,
        linewidth=0.5 * zoom_scale,
        zorder=4
    )

    # plot initial active fireline
    fline_reprojected.plot(
        ax=ax,
        facecolor="none",
        edgecolor="lime",
        transform=data_crs,
        linewidth=0.7 * zoom_scale,
        zorder=5
    )

    # crop
    ax.set_extent([left, right, bottom, top], crs=data_crs)

    # make colorbar and bottom legend
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1, axes_class=maxes.Axes)
    cbar_mappable = cm.ScalarMappable(cmap=plt.colormaps['Oranges'])
    cbar_mappable.set_array(toa_matrix)
    fig.colorbar(cbar_mappable, cax=cax, label="Time of Arrival (hours)")

    lax = divider.append_axes("bottom", size="5%", pad=0.6, axes_class=maxes.Axes)
    lax.axis("off")
    proxy_handles = [
        Line2D([0], [0], color='blue', lw=1, label='Initial Cold Perimeter'),
        Line2D([0], [0], color='lime', lw=1, label='Initial Fireline'),
        Line2D([0], [0], color='red', lw=1, label='Final Perimeter')
    ]
    labels = [h.get_label() for h in proxy_handles]
    legend = lax.legend(
        proxy_handles, 
        labels,
        loc='center',
        ncol=4,
        fontsize=8,
        frameon=True,
        bbox_to_anchor=(0.5, 0.5)
    )

    # add titles
    fig.suptitle(
        "Time of Arrival Field of " + sup_title,
        fontsize=10,
        fontweight='bold'
    )

    # add lat/lon grid
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), 
                  linestyle='--', color='black', alpha=0.3)
    gl.top_labels = False
    gl.right_labels = False

    fname = args.fid + "_interval_" + str(args.start_day).replace(".", "_") + "_test.png"
    output_file = os.path.join(args.output_path, fname)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def spread_one_fire(args):
    # retrieve unique fire observation intervals from FEDS geopackages
    gpkg_path = get_gpkg_path(args)
    start_index, end_index, timesteps = get_pyrostack_indices(args, gpkg_path)
    
    # preprocess input data
    rasters_dict, pre_burned, ignition = load_and_transform_rasters(args, (start_index, end_index))
    
    # set foliar moisture (no observation available)
    rasters_dict["foliar_moisture"] = {"array" : 1.0}

    # set spread adjustment factor (e.g., prescribed burn window during daytime)
    rasters_dict["weather_spread_adjustment"] = {"array" : np.ones(end_index - start_index)} #NOTE: here we prescribe no adjustment
    
    # set up dimensions
    sample_shape = rasters_dict["slope"]["array"].shape
    cube_shape = (timesteps + 1, sample_shape[1], sample_shape[2])
    
    cube_resolution = ( 
        60, # band_duration: minutes
        30, # cell_height:   meters
        30, # cell_width:    meters
    )
    
    # specify cube refresh rates
    cube_refresh_rates = {
        "wind_speed_10m"           : 1.0 / 15.0,
        "upwind_direction"         : 1.0 / 15.0,
        "fuel_moisture_dead_1hr"   : 1.0 / 30.0,
        "temperature"              : 1.0 / 30.0,
        "fuel_spread_adjustment"   : 0.0,
        "weather_spread_adjustment": 1.0 / 30.0,
    }
    
    spot_config = {
        "random_seed"                 : 1234567890,
        "firebrands_per_unit_heat"    : 1e-9,       # firebrands/kJ
        "downwind_distance_mean"      : 10.0,       # meters
        "fireline_intensity_exponent" : 0.3,        # downwind_distance_mean multiplier [I^fireline_intensity_exponent]
        "wind_speed_exponent"         : 0.55,       # downwind_distance_mean multiplier [U^wind_speed_exponent]
        "downwind_variance_mean_ratio": 425.0,      # meters^2 / meter [downwind_variance_mean_ratio = Var(X) / E(X)]
        "crosswind_distance_stdev"    : 100.0,      # meters
        "decay_distance"              : 200.0,      # meters
    }
    
    # create SpaceTimeCube objects from input fields
    space_time_cubes = {name: SpaceTimeCube(cube_shape, array["array"])
                        for (name, array) in rasters_dict.items()}
    
    # create a SpreadState object and specify ignition from the FEDS fireline
    phi = np.where(ignition == 1, -1, 1).astype(np.float32)
    spread_state = els.SpreadState(cube_shape).ignite_cells(lower_left_corner=(0, 0),
                                                            ignition_matrix=phi)

    # run simulation
    runtime_start       = time.perf_counter()
    fire_spread_results = els.spread_fire_with_phi_field(space_time_cubes,
                                                         spread_state,
                                                         cube_resolution,
                                                         30,              #account for 30 min offset between VIIRS overpass and PyroStack layers
                                                         timesteps * 60,
                                                         spot_config=spot_config,
                                                         surface_lw_ratio_model="rothermel",
                                                         cube_refresh_rates=cube_refresh_rates)
    runtime_stop        = time.perf_counter()
    print("Simulation finished. Runtime = " + str(runtime_stop - runtime_start) + " seconds")

    # extract output matrices
    
    spread_state           = fire_spread_results["spread_state"]   # updated SpreadState object (mutated from inputs)
    output_matrices        = spread_state.get_full_matrices(layers=["fire_type", "spread_rate", "fireline_intensity", "flame_length", "time_of_arrival", "phi"])
    fire_type_matrix       = output_matrices["fire_type"] # 0 no fire 1 surface fire 2 passive crown fire 3 active crown fire
    spread_rate_matrix     = output_matrices["spread_rate"] # in m/s
    fireline_intensity_matrix = output_matrices["fireline_intensity"] # kW/m
    flame_length_matrix    = output_matrices["flame_length"] # how tall in m
    time_of_arrival_matrix = output_matrices["time_of_arrival"] - 30 # minutes (subtract out offset)
    np.divide(time_of_arrival_matrix, 60.0, out=time_of_arrival_matrix)   # hours
    
    # generate simple example plot
    make_simple_plot(args, pre_burned, ignition, time_of_arrival_matrix, gpkg_path)

def main():
    args = parse_args()

    spread_one_fire(args)

if __name__ == "__main__":
    main()
