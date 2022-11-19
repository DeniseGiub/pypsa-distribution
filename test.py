#%%
#Load yaml files
import os
import yaml

with open(r'config.yaml') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
# %%
#This script calculates the population in the mini grid 
import geopandas as gpd
from tqdm import tqdm
microgrid_shape = gpd.read_file(f"C://Users//denis//OneDrive//Desktop//Mini grids//pypsa-distribution//mydata.geojson") 
microgrid_shape=microgrid_shape.rename(columns={'subType': 'Geometry'}).rename(columns={'geometry' : 'Centre'})
#TODO add a column with the minigrid name
microgrid_name=config["microgrids_list"]["micA"]["name"]
country_codes=["BW"]
# %%
def _init_process_pop(df_gadm_, year_, worldpop_method_):
    global df_gadm, year, worldpop_method
    df_gadm, year, worldpop_method = df_gadm_, year_, worldpop_method_
#%%
import requests
import shutil
def download_WorldPop_standard(
    country_codes,
    wordlpop_method,
    year=2020,
    update=True,
    out_logging=False,
    size_min=300,
):
    WorldPop_filename = f"srb_ppp_{year}_UNadj_constrained.tif"
    WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/SRB/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/SRB/{WorldPop_filename}",
        ]

    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    WorldPop_filename : str
        Name of the file
    """

    if country_codes == "XK":
        WorldPop_filename = f"srb_ppp_{year}_UNadj_constrained.tif"
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/SRB/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/SRB/{WorldPop_filename}",
        ]

    WorldPop_inputfile = os.path.join(
        os.getcwd(), "data", "WorldPop"
    )  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
    
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break

    return WorldPop_inputfile

download_WorldPop_standard(
    country_codes,
    config["build_shape_options"]["worldpop_method"],
    config["build_shape_options"]["worldpop_method"],
    update=True,
    out_logging=False,
    size_min=300,
)

def _init_process_pop(df_gadm_, year_, worldpop_method_):
    global df_gadm, year, worldpop_method
    df_gadm, year, worldpop_method = df_gadm_, year_, worldpop_method_

#%%
def generalized_mask(src, geom, **kwargs):
    "Generalize mask function to account for Polygon and MultiPolygon"
    if geom.geom_type == "Polygon":
        return mask(src, [geom], **kwargs)
    elif geom.geom_type == "MultiPolygon":
        return mask(src, geom.geoms, **kwargs)
    else:
        return mask(src, geom, **kwargs)

#%%
import numpy as np
def _sum_raster_over_mask(shape, img):
    """
    Function to sum the raster value within a shape
    """
    # select the desired area of the raster corresponding to each polygon
    # Approximation: the population is measured including the pixels
    #   where the border of the shape lays. This leads to slightly overestimate
    #   the output, but the error is limited and it enables halving the
    #   computational time
    out_image, out_transform = generalized_mask(
        img, shape, all_touched=True, invert=False, nodata=0.0
    )
    # calculate total output in the selected geometry
    out_image[np.isnan(out_image)] = 0
    out_sum = out_image.sum()
    # out_sum = out_image.sum()/2 + out_image_int.sum()/2

    return out_sum

def _process_func_pop(c_code):

    # get subset by country code
    country_rows = df_gadm.loc[df_gadm["country"] == c_code].copy()

    # get worldpop image
    WorldPop_inputfile, WorldPop_filename = download_WorldPop_standard(
        c_code, worldpop_method, year, False, False
    )

    with rasterio.open(WorldPop_inputfile) as src:

        for i in country_rows.index:
            country_rows.loc[i, "pop"] = _sum_raster_over_mask(
                country_rows.geometry.loc[i], src
            )

    return country_rows

import multiprocessing as mp
import rasterio

def add_population_data(                        
    microgrid_shape,
    country_codes,
    worldpop_method,
    year=2020,
    update=False,
    out_logging=False,
    nprocesses=2,
    disable_progressbar=False,
   
):

    """
    Function to add population data to arbitrary number of shapes in a country

    Inputs:
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Outputs:
    --------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["pop"]
    """

    # initialize new population column
    microgrid_shape["pop"] = 0.0

    tqdm_kwargs = dict(
        ascii=False,
        unit=" countries",
        # total=len(country_codes),
        desc="Compute population ",
    )

    kwargs = {
            "initializer": _init_process_pop,
            "initargs": (microgrid_shape, year, worldpop_method),
            "processes": nprocesses,
        }
    with mp.get_context("spawn").Pool(**kwargs) as pool:
            if disable_progressbar:
                _ = list(pool.map(_process_func_pop, country_codes))
                for elem in _:
                    microgrid_shape.loc[elem.index, "pop"] = elem["pop"]
            else:
                _ = list(
                    tqdm(
                        pool.imap(_process_func_pop, country_codes),
                        total=len(country_codes),
                        **tqdm_kwargs,
                    )
                )
                for elem in _:
                    microgrid_shape.loc[elem.index, "pop"] = elem["pop"]

#%%
add_population_data(
            microgrid_shape,
            country_codes,
            config["build_shape_options"]["worldpop_method"],
            config["build_shape_options"]["year"],
            config["build_shape_options"]["update_file"],
            config["build_shape_options"]["out_logging"],
            config["build_shape_options"]["nprocesses"],
        )