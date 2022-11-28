
import geopandas as gpd
import fiona
import rasterio
import rasterio.mask
import georasters as gr
import geojson
import pandas as pd
import yaml
from geojson import FeatureCollection
from create_network import n

with open(f"C://Users//denis//OneDrive//Desktop//Mini grids//pypsa-distribution//Scripts//config.yaml") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

# #This script downloads WorldPop data for SL country. 2019 data are taken, since 2020 data are not available for Sierra Leone

# import requests 
# import shutil
# import os

# def download_WorldPop_standard(
#     country_code,
#     year=2019,
#     update=False,
#     size_min=300
# ):
#     """
#     Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

#     Parameters
#     ----------
#     country_code : str
#         Two letter country codes of the downloaded files.
        # Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
#     year : int
#         Year of the data to download
#     update : bool
#         Update = true, forces re-download of files
#     size_min : int
#         Minimum size of each file to download
#     Returns
#     -------
#     WorldPop_inputfile : str
#         Path of the file
#     WorldPop_filename : str
#         Name of the file
#     """
#     WorldPop_filename = f"sle_ppp_{year}_constrained.tif"
#     WorldPop_urls = [
#             f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2019/BSGM/SLE/{WorldPop_filename}"]
#     WorldPop_inputfile = os.path.join(
#         os.getcwd(), "data", "WorldPop", WorldPop_filename
#     )  # Input filepath tif

#     if not os.path.exists(WorldPop_inputfile) or update is True:
    
#         #  create data/osm directory
#         os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

#         loaded = False
#         for WorldPop_url in WorldPop_urls:
#             with requests.get(WorldPop_url, stream=True) as r:
#                 with open(WorldPop_inputfile, "wb") as f:
#                     if float(r.headers["Content-length"]) > size_min:
#                         shutil.copyfileobj(r.raw, f)
#                         loaded = True
#                         break

#     return WorldPop_inputfile, WorldPop_filename


# download_WorldPop_standard(
#                     config["countries"], config["year"], False, 300)

#%%
my_feature=[]

#I create a rectangle
my_feature={
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates":  [
        [
          [config["microgrids_list"]["micA"]["Coordinates"]["Point1"]["lon1"], config["microgrids_list"]["micA"]["Coordinates"]["Point1"]["lat1"]],
          [config["microgrids_list"]["micA"]["Coordinates"]["Point2"]["lon2"], config["microgrids_list"]["micA"]["Coordinates"]["Point2"]["lat2"]],
          [config["microgrids_list"]["micA"]["Coordinates"]["Point3"]["lon3"], config["microgrids_list"]["micA"]["Coordinates"]["Point3"]["lat3"]],         
          [config["microgrids_list"]["micA"]["Coordinates"]["Point4"]["lon4"], config["microgrids_list"]["micA"]["Coordinates"]["Point4"]["lat4"]],
        ]
    ]
  },

  "properties": {
    "Microgrid": config["microgrids_list"]["micA"]["name"]
  }
}

FeatureCollection = FeatureCollection([my_feature])

with open('microgrid_shape.geojson', 'w') as f:
        f.write(geojson.dumps(FeatureCollection))

gdf = gpd.read_file('microgrid_shape.geojson')
gdf.to_file('microgrid_shape.shp')

with fiona.open("microgrid_shape.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

with rasterio.open("sle_ppp_2019_constrained.tif") as src:
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
    out_meta = src.meta

out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

with rasterio.open("SL.masked.tif", "w", **out_meta) as dest:
    dest.write(out_image)

myRaster = 'sle_ppp_2019_constrained.tif'
total_pop= gr.from_file(myRaster)
    
total_pop=total_pop.to_geopandas() 

total_pop=(total_pop['value'].sum()) #Total SL population

myRaster = 'SL.masked.tif'
pop_microgrid = gr.from_file(myRaster)
    
pop_microgrid=pop_microgrid.to_geopandas() 

pop_microgrid=(pop_microgrid['value'].sum()) #Microgrid population

#I import the dataframe of electricity demand for Africa
import pandas as pd
df_demand=pd.read_excel(r'C:\Users\denis\OneDrive\Desktop\Mini grids\pypsa-distribution\Africa.xlsx', index_col = None)

#I select the rows related to Benin (since there are no data for SL)
df_demand_SL=df_demand.loc[26280:35039, :]

# I select the column "electricity demand"
df_demand_SL=df_demand_SL["Electricity demand"]

df_demand_SL=pd.DataFrame(df_demand_SL)

p=(pop_microgrid/total_pop)*100 #Coefficient

electric_load=df_demand_SL/p #Electric load of the minigrid

electric_load_xlsx=electric_load.to_excel('electric_load_1.xlsx', index=False) 
