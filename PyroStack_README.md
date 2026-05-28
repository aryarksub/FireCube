# PyroStack

PyroStack is a multi-band spatiotemporal dataset recording fire spread and covariate information for fires across the contiguous United States (CONUS) and Alaska between 2012-2024. In total, it contains data for 6994 fires.

## Dataset Structure

The top-level directory is named `pyrostack`. Within this directory, fire information is stored in ZIP files corresponding to the fire year. All fires in the year `YYYY` are stored in the `fires_{YYYY}.zip` file. We also supply a `manifest.csv` file which contains identifying information for each fire, including its name, unique MTBS identifier, size, duration, and ignition latitude and longitude. 

Unzipping a particular ZIP file will result in a set of directories, each of which correspond to a single fire. Within this fire directory, there is a `fire_times.csv` file, as well as five subdirectories corresponding to particular data categories:
- `fire_spread`: Fire observations and characteristics layers, originally sourced using the FEDS algorithm (Chen et al., 2022) on MTBS data (Chen et al., 2026)
- `low_res_climate`: Low resolution climate data, originally sourced from ERA5-Land (Mu&ntilde;oz Sabater et al., 2021)
- `high_res_climate`: High resolution climate data, originally sourced from RTMA via the CloudFire Inc. Worldgen Server
- `fuel_structure`: Fuel structure (canopy) data, originally sourced from LANDFIRE via the CloudFire Inc. Worldgen Server
- `veg_fm_topo`: Vegetation, fuel model, and topography data, originally sourced from LANDFIRE (Ryan and Opperman, 2013)

All of the data in these subdirectories are stored as TIF files. The `fire_times.csv` file contains two columns: `time` and `feds`. The `time` column indicate the date and hour corresponding to data observations, and the `feds` column is a boolean flag indicating whether the corresponding time coincides with a FEDS fire observation. 

TIFs in the `fire_spread` subdirectory have multiple bands, each of which correspond to a 12-hour fire observation at times with a `True` value in the `fire_times.csv` file. TIFs in the `low_res_climate` and `high_res_climate` subdirectories also have multiple bands, each of which correspond to the times in each row in the `fire_times.csv` file. TIFs in the `fuel_structure` and `veg_fm_topo` subdirectories are single band and correspond to static data layers.

## Data Overview

The below table shows the spatial and temporal resolutions of data stored in each subdirectory

| Subdirectory | Spatial Resolution | Temporal Resolution |
| :---: | :---: | :---: |
| `fire_spread` | 300m | 12-hr |
| `low_res_climate` | 9000m | 1-hr |
| `high_res_climate` | 600m | 1-hr |
| `fuel_structure` | 30m | Static |
| `veg_fm_topo` | 30m | Static |

The below table summarizes the variables found in each subdirectory.

| Layer Name | Description | Subdirectory | Units |
| :---: | :---: | :---: | :---: |
| `farea` | Fire area | `fire_spread` | Binary |
| `fline` | Active fire line | `fire_spread` | Binary |
| `nfp` | New fire pixels | `fire_spread` | Binary |
| `frp` | Fire radiative power | `fire_spread` | W/m<sup>2</sup> |
| `d2m` | 2-meter dewpoint temperature | `low_res_climate` | Kelvin (K) |
| `sp` | Surface pressure | `low_res_climate` | Pascal (Pa) |
| `t2m` | 2-meter temperature | `low_res_climate` | Kelvin (K) |
| `tp` | Total precipitation | `low_res_climate` | Meters/hr (m/hr) |
| `lh` | Live herbaceous fuel moisture | `high_res_climate` | Percent (%) |
| `lw` | Live woody fuel moisture | `high_res_climate` | Percent (%) |
| `m1` | 1-hour dead fuel moisture | `high_res_climate` | Percent (%) |
| `m10` | 10-hour dead fuel moisture | `high_res_climate` | Percent (%) |
| `m100` | 100-hour dead fuel moisture | `high_res_climate` | Percent (%) |
| `wd` | 20-foot wind direction | `high_res_climate` | Degrees (°) |
| `ws` | 20-foot wind speed | `high_res_climate` | Miles per hour (mph) |
| `cbd` | Canopy bulk density | `fuel_structure` | 100 kg/m<sup>3</sup> |
| `cbh` | Canopy base height | `fuel_structure` | 10-m |
| `cc` | Canopy cover | `fuel_structure` | Percent (%) |
| `ch` | Canopy height | `fuel_structure` | 10-m |
| `*evt*` | Existing vegetation type | `veg_fm_topo` | Categorical |
| `*f[bfm]13*` | 13 Anderson Fire Behavior Fuel Models (2020) | `veg_fm_topo` | Categorical |
| `*f[bfm]40*` | 40 Scott and Burgan Fire Behavior Fuel Models (2020) | `veg_fm_topo` | Categorical |
| `*roads*` | Operational roads | `veg_fm_topo` | Categorical |
| `*asp*` | Topographic aspect | `veg_fm_topo` | Degrees (°) |
| `*elev*` | Topographic elevation | `veg_fm_topo` | Meters (m) |
| `*slpd*` | Topographic slope | `veg_fm_topo` | Degrees (°) |

Asterisks in the layer name indicate that there may prefixes and/or suffixes in the name of the corresponding TIF files for some of the fires. The `fbfm13` and `fbfm40` layers may also be recorded as `f13` and `f40`, respectively. 

All data is recorded using the EPSG:5070 Coordinate Reference System. For a single fire, all data is also recorded within the same spatial bounding box, as defined by the fire's final perimeter, along with some additional padding. 