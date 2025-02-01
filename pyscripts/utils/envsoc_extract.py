"""

Select areas defined by the World Database on Proctected Areas (WDPA) and Landscan

Description:
This script performs the following tasks:
1. Read shapefiles from the WDPA
2. Read shapefiles from Landscan
2. Combine and select regions within specified geographical boundaries to support the creation of clusters

Usage:
1. Configure parameters in the header of extract_shp_envsoc.py
2. Run the current script (check the instances below)
3. The new shapefile will be support in the folder specified in the header

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
import os
import time
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import matplotlib.pyplot as plt

# get existing path
script_dir = os.path.dirname(__file__)  # absolute directory

# HEADER
lat_min, lat_max = -36, 7
lon_min, lon_max = -74, -30
out_shp_pol=os.path.join(script_dir,'../../support/shapefile/envsoc/envsoc_poly.shp')
out_shp_pnt=os.path.join(script_dir,'../../support/shapefile/envsoc/envsoc_point.shp')
wdpa_pol1=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_0/WDPA_Feb2024_Public_shp-polygons.shp')
wdpa_pol2=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_1/WDPA_Feb2024_Public_shp-polygons.shp')
wdpa_pol3=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_2/WDPA_Feb2024_Public_shp-polygons.shp')
wdpa_pnt1=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_0/WDPA_Feb2024_Public_shp-points.shp')
wdpa_pnt2=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_1/WDPA_Feb2024_Public_shp-points.shp')
wdpa_pnt3=os.path.join(script_dir,'../../support/shapefile/WDPA/WDPA_Feb2024_Public_shp (1)/WDPA_Feb2024_Public_shp_2/WDPA_Feb2024_Public_shp-points.shp')

# read shapefiles
shp_pol1 = gpd.read_file(wdpa_pol1)
shp_pol2 = gpd.read_file(wdpa_pol2)
shp_pol3 = gpd.read_file(wdpa_pol3)
shp_pnt1 = gpd.read_file(wdpa_pnt1)
shp_pnt2 = gpd.read_file(wdpa_pnt2)
shp_pnt3 = gpd.read_file(wdpa_pnt3)

# select boundaries for all
sel_boundaries = box(lon_min, lat_min, lon_max, lat_max)

# clip the polygons with the box
clipped_pol = gpd.clip(pd.concat([shp_pol1, shp_pol2, shp_pol3]), sel_boundaries)
clipped_pnt = gpd.clip(pd.concat([shp_pnt1, shp_pnt2, shp_pnt3]), sel_boundaries)

# save clipped shapefiles
clipped_pol.to_file(out_shp_pol)
clipped_pnt.to_file(out_shp_pnt)
