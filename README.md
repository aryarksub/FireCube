# PyroStack

Code used to generate the PyroStack wildfire trajectory dataset. For a detailed description of the PyroStack dataset's contents, view `Dataset_README.md`.

## Setup

Install the core datacube-generation environment with uv:

```bash
uv sync
```

Install the optional dependencies for the machine learning training code:

```bash
uv sync --extra ml
```

Run commands through uv so they use the project environment:

```bash
uv run python -c "import main"
```

## Repository Structure

We separate the input files used as an initial point for PyroStack generation from the actual code files that house scripts to run the entire data acquisition and cleaning process. 

### Input Files

Input files with broader metadata information about fires are stored in the `inputData` directory. There are four such files:
- `LF2024_EVT.csv`: Contains integer codes for existing vegetation type (EVT) as used by the 2024 LANDFIRE version. This file is primarily used as an aid for external visualization code and is not part of the data acquisition pipeline.
- `fireslist2012-2024.csv`: Complete list of 7,739 fires as noted in the FEDS-MTBS dataset. Information for each fire includes its year (`Year`), MTBS incident name (`Incid_Name`), MTBS unique identifier (`Event_ID`), final burned area in acres (`BurnBndAc`), start date (`tst`), end date (`ted`), minimum longitude (`lon0`), maximum longitude (`lon1`), minimum latitude (`lat0`), and maximum latitude (`lat1`).
  - The latitude-longitude pairs (`lon0`, `lat0`) and (`lon1`, `lat1`) correspond to the bottom-left and top-right spatial coordinates for a bounding box of the fire's final burned area.
- `fireslist2012-2024_withtype_no_fline_hawaii.csv`: Complete list of 6,994 fires recorded in PyroStack. This file is a subset of the above with fires without any active fireline data removed, along with fires occurring in Hawaii. The schema is the same, with the addition of the `IncidType` column, which indicates the type of fire ("Prescribed Fire" versus "Wildfire").
- `layer_ranges.csv`: Contains mimimum and maximum possible values for each data layer stored in PyroStack, along with how missing data (NaN) values are represented. This file is primarily used for validating downloaded data after the data acquisition process.

### Source Code

All source code for data acquisition, dataset loading, and model runs can be found in the `src` directory. The pipeline to generate PyroStack by downloading data from external sources is stored in `main.py`. We also store several subdirectories that serve separate purposes:
- `era5`: Stores code for downloading low-resolution climate data from the ERA5-Land dataset via the Copernicus Climate Data Store API (see `era5.py`).
- `landfire`: Stores code for downloading vegetation, fuel model, and topography data from the LANDFIRE API (see `landfire.py`). Also contains a subdirectory `archive_processing` that stores code for obtaining data from locally downloaded full-extent LANDFIRE datasets when the API does not provide data from updated versions. This subdirectory contains its own README file with more details.
- `ml`: Stores code related to loading the PyroStack dataset (see `data.py`) as well as other relevant code for training and evaluating a basic Vision Transformer model for one-step-ahead prediction.
- `pyregence`: Stores code for downloading fuel structure and high-resolution climate data from the CloudFire Inc. Worldgen server (see `pyregence.py`) as well as other helper files.
- `pyretechnics`: Stores code for running the Pyretechnics physical fire behavior/spread model (see `run_pyretechnics.py`).
- `util`: Stores code for utility functions that are used across several other code files. Notably, `feds_util.py` is used to rasterize original satellite observations to obtain fire spread data layers in PyroStack. Also, visualization is handled via code stored under the `vis_util` subdirectory. 

