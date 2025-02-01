"""

Extract data from shapefile provided by the regulator into a Pandas DataFrame

Description:
This script performs the following tasks:
1. Read shapefile and extract 

Usage:
1. Configure the parameters of simulation using the pv_config.py
2. Run the current script (check the instances below)
3. Outputs will be used in the workbook_solar.xlsx

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
from pv_stat_config import shp_path_in, shp_path_out
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
dfp = dfp.loc[dfp['FASE_USINA']=='Operação']

# select the columns to be stored
dfsel = dfp[['NOME', 'CEG', 'UF', 'INIC_OPER', 'POT_KW', 'POT_FISC_K', 'ID_EMPREEN', 'FASE_USINA', 'lon', 'lat']]
dfsel = dfsel.sort_values('NOME')
dfsel.to_csv(shp_path_out, index=False, encoding="utf-8-sig")
print("[+] Done")
