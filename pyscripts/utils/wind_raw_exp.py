"""

Extract data from shapefile provided by the regulator into a Pandas DataFrame

Description:
This script performs the following tasks:
1. Read shapefile and extract 

Usage:
1. Configure the parameters of simulation using the pv_config.py
2. Run the current script (check the instances below)
3. Outputs will be used in the workbook_wind.xlsx

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com


"""
from wind_stat_config import shp_path_in, shp_path_out
import geopandas as gpd
import pandas as pd

# begin
print("[+] Begin")
try:
   # create a geodataframe
   df = gpd.read_file(shp_path_in, encoding='UTF-8')
except FileNotFoundError as e:
   print(f"[!] Error loading files: {e}")
   exit(1)

# extract the points from the current geometry
df['lon'] = df['geometry'].x
df['lat'] = df['geometry'].y

# convert the geodataframe to a pandas dataframe
dfp = pd.DataFrame(df)

# select only existing projects
dfp = dfp.loc[dfp['OPERACAO']=='Sim']

# select the columns to be stored
dfsel = dfp[['POT_MW', 'ALT_TOTAL', 'ALT_TORRE', 'DIAM_ROTOR', 'DATA_ATUAL', 'EOL_VERSAO', 'NOME_EOL', 'DEN_AEG', 'OPERACAO', 'UF', 'CEG', 'lon', 'lat']]
dfsel = dfsel.sort_values('NOME_EOL')
dfsel.to_csv(shp_path_out, index=False, encoding="utf-8-sig")
print("[+] Done")
