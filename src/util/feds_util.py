import geopandas as gpd
import os

dir_feds25 = os.path.join('inputData', 'FEDS2.5')
feds_firelist = os.path.join('inputData', 'feds2.5_firelist.csv')

def set_gdffile( Event_ID):
    """
    Constructs the file path for a GeoPackage file based on the given year and Event ID.

    Args:
        yr (int): The year associated with the file.
        Event_ID (str): The unique identifier for the event.

    Returns:
        str: The constructed file path for the GeoPackage file.
    """
    fnm = os.path.join(dir_feds25, Event_ID + '.gpkg')
    return fnm

def check_gdffile( Event_ID):
    """
    Checks if a GDF file exists for a given year and event ID.

    Args:
        yr (int): The year of the event.
        Event_ID (str): The unique identifier for the event.

    Returns:
        bool: True if the GDF file exists, False otherwise.
    """
    fnm = set_gdffile( Event_ID)
    return os.path.exists(os.path.expanduser(fnm))

def set_dd(layer="perimeter"):

    # diagnostic data name and types saved in geojson files (in addition to geometries)
    if layer == "perimeter":
        dd = {
            "t": "datetime64[ns]",
            "n_pixels": "int",  # number of total pixels
            "n_newpixels": "int",  # number of new pixels
            "farea": "float",  # fire size
            "fperim": "float",  # fire perimeter length
            "duration": "float",  # fire duration
            "pixden": "float",  # fire pixel density
            "sumFRP": "float",  # mean FRP of the new fire pixels
            "flinelen": "float",  # active fire front line length
        }
    elif layer == "fireline":
        dd = {
            "t": "datetime64[ns]",
        }
    elif layer == "newfirepix":
        dd = {
            "t": "datetime64[ns]",
        }

    return dd

def read_gdf_fire(Event_ID,layer='perimeter'):
    """ read 1 layer of a output file for a single fire event
    """
    dd = set_dd(layer=layer)
    
    fnm = set_gdffile( Event_ID)
    if check_gdffile( Event_ID):
        gdf = gpd.read_file(fnm, layer=layer)
        for k, tp in dd.items():
            gdf[k] = gdf[k].astype(tp)        
    else:
        gdf = None
        print('File does not exist')
    return gdf

def read_1fire(Event_ID):
    """ read all three layers of a output file for a single fire event
    """
    gdf_fperim_rd = read_gdf_fire(Event_ID, layer='perimeter')
    gdf_fline_rd = read_gdf_fire(Event_ID, layer='fireline')
    gdf_nfp_rd = read_gdf_fire(Event_ID, layer='newfirepix')
    return gdf_fperim_rd, gdf_fline_rd, gdf_nfp_rd