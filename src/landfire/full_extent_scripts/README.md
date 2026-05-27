# PyroStack Layer Production from LANDFIRE Full-Extent Archives

This directory contains scripts designed to extract PyroStack LANDFIRE layers efficiently from full-extent LANDFIRE archives.

## Overview

The workflow automates the process of categorizing fire events by their optimal LANDFIRE version and slicing the full-extent datasets into fire-specific domains.

## Script Reference

* `build_fire_list.py`  
    Retrieves all fires from the PyroStack fire list that correspond with a specified LANDFIRE version, compiling their metadata used for extraction.
* `extract_tifs.py`  
    Loads full-extent `.tif` files from the local LANDFIRE archives (which must be pre-downloaded from [landfire.gov](https://www.landfire.gov/)) into memory and rapidly slices them into each fire's respective spatial domain.
* `process_lf_version.sh`  
    A shell script that streamlines the pipeline, allowing you to iteratively process all variables for a single LANDFIRE version in a single execution.