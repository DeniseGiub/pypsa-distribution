#%%
import geopandas as gpd
import fiona
import rasterio
import rasterio.mask
import georasters as gr
import geojson
from geojson import FeatureCollection
import yaml

with open(r'config.yaml') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

my_feature=[]

#I create a rectangle
my_feature={
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates":  [
        [
          [-12.023283,7.951067],
          [-11.450478, 7.951067],
          [-11.450478, 7.714312],         
          [-12.023283, 7.714312]
        ]
    ]
  },
  "properties": {
    "Microgrid": config["microgrids_list"]["micA"]["name"]
  }
}

FeatureCollection = FeatureCollection([my_feature])

with open('mydata.geojson', 'w') as f:
        f.write(geojson.dumps(FeatureCollection))

gdf = gpd.read_file('mydata.geojson')
gdf.to_file('file.shp')

with fiona.open("file.shp", "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

with rasterio.open("SL.tif") as src:
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
    out_meta = src.meta

out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

with rasterio.open("sle.masked.tif", "w", **out_meta) as dest:
    dest.write(out_image)

myRaster = 'SL.tif'
total_pop= gr.from_file(myRaster)
    
total_pop=total_pop.to_geopandas()

total_pop=(total_pop['value'].sum())

myRaster = 'sle.masked.tif'
pop_microgrid = gr.from_file(myRaster)
    
pop_microgrid=pop_microgrid.to_geopandas()

pop_microgrid=(pop_microgrid['value'].sum())

#Percentage of population inside the microgrid 

p=(pop_microgrid/total_pop)*100


# %%
