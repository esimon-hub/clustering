"""

Benchmark simulations against data (capacity factors) provided by the system operator (hourly, daily, and monthly)

Description:
This script performs the following tasks:
1. Read raw data downloaded and preprocessed from the ISO (ONS)
2. Read simulated data from the files (match the files provided by the simulation)

Usage:
1. Configure the parameters of simulation using the pv_stat_config.py
2. Create the folder "Simulated" within "/data/ACT/Solar/Simulated"
3. Run the current script (check the instances below)
4. Outputs will be stored in within the following folder:
   "$CURRENT_PATH$/Solar/benchmark/"

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
import pandas as pd
import os
import numpy as np
from pv_stat_config import yr_beg, yr_end
from pv_stat_config import in_cf_act, in_asset, in_folder_cf_model, in_cap_model
from pv_stat_config import out_cl_iso_h, out_cl_sim_h, out_sin_iso_h, out_sin_sim_h
from pv_stat_config import out_cl_iso_d, out_cl_sim_d, out_sin_iso_d, out_sin_sim_d
from pv_stat_config import out_cl_iso_m, out_cl_sim_m, out_sin_iso_m, out_sin_sim_m
from datetime import datetime

# begin
print("[+] Begin")
# read dataframes needed to perform benchmark analysis
df_asset = pd.read_csv(in_asset, index_col=False)
df_iso_h = pd.read_csv(in_cf_act, index_col=False)
df_model_cap = pd.read_csv(in_cap_model, index_col=False)

# read multiple files from prior simulation
dataframes = []

# iterate through the years, read each file and append to the list
for year in range(yr_beg, yr_end + 1):
    file_name = f'solar_cf_{year}.csv'
    file_path = os.path.join(in_folder_cf_model, file_name)

    # check if the file exists before attempting to read
    if os.path.exists(file_path):
        try:
            df = pd.read_csv(file_path)
            dataframes.append(df)
        except Exception as e:
            # if there is error reading the file, print an error message
            print(f'[+] Error reading {file_name}: {e}')
    else:
        # if the files does not exist, print a warning message
        print(f'[+] File does not exist: {file_path}')

# concatenates all the dataframes in the list into a single dataframe, if any
if dataframes:
    df_model_cf_h = pd.concat(dataframes, ignore_index=True)
else:
    print('No files were found or successfully read. No combined DataFrame created.')

# adjust column for datetime
df_iso_h['datetime'] = pd.to_datetime(df_iso_h['din_instante'])

# create a time series to make sure all hours are available in the dataframe (ISO data has a few missing hours since 2010)
dtime1=datetime(yr_beg,1,1)
dtime2=datetime(yr_end,12,31,23,0,0)
date_range_h = pd.date_range(start=dtime1, end=dtime2, freq='H')

# time series to be used along with merger
ts_h=pd.Series(date_range_h, name='datetime')

# reframe the format of the datasets to become rows (hours) and columns (clusters)
col_names=df_iso_h['nom_usina'].drop_duplicates().to_list()

# create the dataframes for benchmark (ISO)
bench_iso_cf_h=pd.DataFrame(columns=col_names, index=pd.to_datetime(date_range_h))
bench_iso_gen_h=pd.DataFrame(columns=col_names, index=pd.to_datetime(date_range_h))
bench_iso_cap_h=pd.DataFrame(columns=col_names, index=pd.to_datetime(date_range_h))

# create the dataframes for benchmarks (Simulated)
bench_model_cf_h=pd.DataFrame(0.0, columns=col_names, index=pd.to_datetime(date_range_h))
bench_model_gen_h=pd.DataFrame(0.0, columns=col_names, index=pd.to_datetime(date_range_h))
bench_model_cap_h=pd.DataFrame(0.0, columns=col_names, index=pd.to_datetime(date_range_h))  #This is an auxiliary set of dataframes to compute capacity factors for clusters (col_names)

# adjust timezone from simulations to GMT -3 (shift all entries to account for the change in timezone)
df_model_cf_h['datetime'] =  pd.to_datetime(df_model_cf_h['datetime'])
df_model_cf_h['datetime'] = df_model_cf_h['datetime'] - pd.Timedelta('3H')
df_model_cf_h = df_model_cf_h.iloc[3:]
new_rows = df_model_cf_h[-24:-21].copy()
new_rows['datetime'] = new_rows['datetime'] + pd.Timedelta(hours=24)
df_model_cf_h = pd.concat([df_model_cf_h,new_rows])
df_model_cf_h.reset_index(inplace=True, drop=True)

# create dataframe for installed capacity to deal with commissioning dates (simulated data)
full_cap = df_model_cap.loc[df_model_cap.index.repeat(len(df_model_cf_h.index))].reset_index(drop=True)
df_model_cap_h = pd.DataFrame(data=full_cap.values, index=df_model_cf_h.index, columns=df_model_cap.columns)
df_asset['COMMISSIONING'] = pd.to_datetime(df_asset['COMMISSIONING'], format='%m/%d/%Y %H:%M')

# cross-check to make sure that the number of data points is enough to perform tasks
min_sim = df_model_cf_h['datetime'].min().year
min_iso = dtime1.year
#if min_sim >= min_iso:
#    print("[!] Limited availability of simulated data!")
#    exit(1)

#
# HOURLY SIMULATIONS VS ACTUALS
#
print("[+] Hourly at cluster level")
# ISO data
for irow in col_names:
    # select cluster from the dataframe
    iso_iter=df_iso_h.loc[df_iso_h['nom_usina']==irow][['datetime','Fator de Carga Comiss','val_capacinstaladacomissconmw','val_geraenergiaconmwmed']]
    if(not iso_iter.empty):
        # create the new structure
        tmp = pd.merge(ts_h, iso_iter, left_on='datetime', right_on='datetime', how='left')
        bench_iso_cf_h[irow]=tmp['Fator de Carga Comiss'].values                     # dataframe used for hourly analysis
bench_iso_cf_h.index.name = 'datetime'
# Simulation data
# perform element-wise multiplication on the other columns
df_model_gen_h = df_model_cf_h.drop(columns=['datetime']) * df_model_cap_h  # datetime column is available only in df_model_cf_h
df_model_gen_h['datetime'] = df_model_cf_h['datetime']
df_model_cap_h['datetime'] = df_model_cf_h['datetime']
for _,irow in df_asset.iterrows():
    # zero data before commissioning period (important to calculate cluster equivalent generation)
    idproj=irow['ID']
    date_cut=irow['COMMISSIONING'].date()
    capacity=irow['POT_KW']
    # store hourly generation
    df_model_gen_h.loc[df_model_gen_h['datetime'].dt.date < date_cut, idproj] = 0.0
    # store hourly capacity
    df_model_cap_h.loc[df_model_cap_h['datetime'].dt.date < date_cut, idproj] = 0.0
# populate the generation and capacity dataframes for the same number of rows
df_model_gen_h = df_model_gen_h.loc[(df_model_gen_h['datetime'].dt.date >= datetime.date(dtime1)) & (df_model_gen_h['datetime'].dt.date <= datetime.date(dtime2))]
df_model_cap_h = df_model_cap_h.loc[(df_model_cap_h['datetime'].dt.date >= datetime.date(dtime1)) & (df_model_cap_h['datetime'].dt.date <= datetime.date(dtime2))]
# reindex to perform operations and interpolate values to full hour
df_model_gen_h.drop(columns=['datetime'], inplace=True)
df_model_gen_h['datetime'] = date_range_h
df_model_gen_h.set_index('datetime', inplace=True)
df_model_gen_h = df_model_gen_h.resample('30T').asfreq()
df_model_gen_h = df_model_gen_h.interpolate(method='linear')
df_model_gen_h = df_model_gen_h.resample('H').first()
df_model_cap_h.drop(columns=['datetime'], inplace=True)
df_model_cap_h['datetime'] = date_range_h
df_model_cap_h.set_index('datetime', inplace=True)
df_model_cap_h = df_model_cap_h.resample('30T').asfreq()
df_model_cap_h = df_model_cap_h.interpolate(method='linear')
df_model_cap_h = df_model_cap_h.resample('H').first()

# calculate adjusted capacity factors
for irow in col_names:
    prjclstr=df_asset.loc[df_asset['CLUSTER']==irow]
    # look into project s belonging to the current cluster
    for _,line in prjclstr.iterrows():
        n_reduced = line['ID']
        bench_model_gen_h[irow] = bench_model_gen_h[irow] + df_model_gen_h[n_reduced]
        bench_model_cap_h[irow] = bench_model_cap_h[irow] + df_model_cap_h[n_reduced]
    # use the same bench_model_h to estimate the cf for the current cluster
    bench_model_cf_h[irow]=bench_model_gen_h[irow]/bench_model_cap_h[irow]
# Store capacity factors for clusters
bench_iso_cf_h.to_csv(out_cl_iso_h, encoding='utf-8-sig')
bench_model_cf_h.to_csv(out_cl_sim_h, encoding='utf-8-sig')

#
# HOURLY SIN LEVELS
#
print("[+] Hourly at SIN level")
# create a structure to store the selected clusters from the ISO
df_iso_h.dropna(subset=['val_capacinstaladacomissconmw'], inplace=True)
sin_iso_gen = df_iso_h[['datetime','val_geraenergiaconmwmed']]
sin_iso_cap = df_iso_h[['datetime','val_capacinstaladacomissconmw']]
# calculate the sum for each hour
sin_iso_gen = sin_iso_gen.groupby(['datetime']).sum()
sin_iso_cap = sin_iso_cap.groupby(['datetime']).sum()
sin_iso_gen.reset_index(inplace=True)                   # drop the datetime for index and keep it as a column
sin_iso_cap.reset_index(inplace=True)
# create the dataframe for hourly cf
sin_iso_gen['cf']=np.where(sin_iso_cap['val_capacinstaladacomissconmw'] != 0,sin_iso_gen['val_geraenergiaconmwmed'] / sin_iso_cap['val_capacinstaladacomissconmw'],np.nan)
# match the entire time series (ts_merge)
sin_iso_h = pd.merge(ts_h, sin_iso_gen[['datetime','cf']], left_on='datetime', right_on='datetime', how='left')
sin_iso_h.set_index('datetime', inplace=True)
# SIMULATION (not just clusters)
# create a structure to store the selected data from the SIMULATION for the entire system (SIN)
sin_sim_gen_h = pd.DataFrame()
#sin_sim_gen['datetime'] = pd.to_datetime(date_range_h)
sin_sim_gen_h['gen'] = df_model_gen_h.sum(axis=1)
sin_sim_cap_h = pd.DataFrame()
#sin_sim_cap['datetime'] = pd.to_datetime(date_range_h)
sin_sim_cap_h['cap'] = df_model_cap_h.sum(axis=1)
sin_sim_h = pd.DataFrame()
#sin_sim_h['datetime'] = pd.to_datetime(date_range_h)
sin_sim_h['cf'] = sin_sim_gen_h['gen']/sin_sim_cap_h['cap']
# Store capacity factors for clusters
sin_iso_h.to_csv(out_sin_iso_h, encoding='utf-8-sig')
sin_sim_h.to_csv(out_sin_sim_h, encoding='utf-8-sig')

#
# DAILY SIMULATIONS VS ACTUALS
#
print("[+] Daily at cluster level")
# ISO
# create a structure to store the selected clusters from the ISO
bench_iso_cf_d = bench_iso_cf_h.resample('D').mean()
bench_iso_cf_d = bench_iso_cf_d.where(bench_iso_cf_d >= 0.001)
# SIMULATION
bench_sim_cf_d = bench_model_cf_h.resample('D').mean()
bench_sim_cf_d = bench_sim_cf_d.where(bench_sim_cf_d >= 0.001)
# Store capacity factors for clusters
bench_iso_cf_d.to_csv(out_cl_iso_d, encoding='utf-8-sig')
bench_sim_cf_d.to_csv(out_cl_sim_d, encoding='utf-8-sig')

#
# DAILY SIN LEVELS
#
print("[+] Daily at SIN level")
# ISO (combined clusters)
sin_iso_d = sin_iso_h.resample('D').mean()
# Model (combined clusters)
sin_sim_d = sin_sim_h.resample('D').mean()
# Store capacity factors for clusters
sin_iso_d.to_csv(out_sin_iso_d, encoding='utf-8-sig')
sin_sim_d.to_csv(out_sin_sim_d, encoding='utf-8-sig')

#
# MONTHLY SIMULATIONS VS ACTUALS
#
print("[+] Monthly at cluster level")
# ISO
# create a structure to store the selected clusters from the ISO
bench_iso_cf_m = bench_iso_cf_h.resample('M').mean()
bench_iso_cf_m = bench_iso_cf_m.where(bench_iso_cf_m >= 0.001)
# SIMULATION
bench_sim_cf_m = bench_model_cf_h.resample('M').mean()
bench_sim_cf_m = bench_sim_cf_m.where(bench_sim_cf_m >= 0.001)
# Store capacity factors for clusters
bench_iso_cf_m.to_csv(out_cl_iso_m, encoding='utf-8-sig')
bench_sim_cf_m.to_csv(out_cl_sim_m, encoding='utf-8-sig')

#
# MONTHLY SIN LEVELS
#
print("[+] Monthly at SIN level")
# ISO (combined clusters)
sin_iso_m = sin_iso_h.resample('M').mean()
# Model (combined clusters)
sin_sim_m = sin_sim_h.resample('M').mean()
# Store capacity factors for clusters
sin_iso_m.to_csv(out_sin_iso_m, encoding='utf-8-sig')
sin_sim_m.to_csv(out_sin_sim_m, encoding='utf-8-sig')


print("[+] Done")

