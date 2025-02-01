"""
Solar PV Parameters

Description:
This script sets the overall assumptions to perform the benchmark analysis:

Usage:
1. Set all folders for actuals and simulated

For support or inquiries, contact emanueel@gmail.com

"""

# Set definitions for solar PV
import os

# get existing path
script_dir = os.path.dirname(__file__)  # absolute directory

#
# downscaling_map.py
#
yr_down1=1999
yr_down2=2018
era5_src=os.path.join(script_dir,'../../data/ERA5/')
gsa_src_ghi=os.path.join(script_dir,'../../data/GSA/GHI.tif')
map_path_folder=os.path.join(script_dir,'figures/maps/')

#
# pv_map.py
#
yr_res1=1994
yr_res2=2023
era5_src=os.path.join(script_dir,'../../data/ERA5/')
map_path_folder=os.path.join(script_dir,'figures/maps/')
shapeon=os.path.join(script_dir,'../../support/shapefile/qgis/onshore_only.shp')

#
# pv_iso_exp.py
#
# select mode to execute
# 1 - process raw data (create month and year columns; reduce the number of columns to relevant data; adjust timezone)
# 2 - calculate averages from step 1 and remove outliers
# 3 - plot and inspection of datasets
# 4 - plot results for all projects
mode_exec=2
# inputs
in_iso_path1=os.path.join(script_dir,'data/Solar/ISO/solar_data_ons_cf_2015_2022.csv')
in_iso_path2_dir=os.path.join(script_dir,'data/Solar/ISO/2023/')
in_window=os.path.join(script_dir,'data/Solar/ISO/iso_window_out.csv')
in_assetexc=os.path.join(script_dir,'data/Solar/ISO/iso_asset_out.csv')
# outputs (location output for raw data processing mode 1)
out_path=os.path.join(script_dir,'data/Solar/ISO/output1_solar.csv')
# output  (location output for raw data processing mode 2)
out_cf_path0=os.path.join(script_dir,'data/Solar/ISO/output2_calc_solar_cluster_H.csv')
# plot selected asset (mode 3)
sel_asset="Conj. Boa Hora"
# output path for figures (mode 4). inspect first, then create in_window and in_assetexc
inspect_figure=os.path.join(script_dir,'figures/inspect/')

#
# pv_raw_exp.py
# 
shp_path_in=os.path.join(script_dir,'data/Solar/_ags_output1/zipfolder/Centrais_Geradoras__Solares_Fotovoltaicas_-_UFV.shp')
shp_path_out=os.path.join(script_dir,'data/Solar/shp_export.csv')

#
# pv_benchmark.py
#
yr_beg=2010
yr_end=2023
in_asset=os.path.join(script_dir,'../../data/ACT/Solar/Asset/pv_assets.csv')
in_folder_cf_model=os.path.join(script_dir,'../../data/ACT/Solar/Simulated/')
in_cap_model=os.path.join(script_dir,'../../data/ACT/Solar/Simulated/solar_capacity.csv')
in_cf_act=os.path.join(script_dir,'data/Solar/ISO/output2_calc_solar_cluster_H.csv')
# outputs (hourly)
out_cl_iso_h=os.path.join(script_dir,'data/Solar/benchmark/h_cl_iso.csv')
out_cl_sim_h=os.path.join(script_dir,'data/Solar/benchmark/h_cl_sim.csv')
out_sin_iso_h=os.path.join(script_dir,'data/Solar/benchmark/h_sin_iso.csv')
out_sin_sim_h=os.path.join(script_dir,'data/Solar/benchmark/h_sin_sim.csv')
# outputs (daily)
out_cl_iso_d=os.path.join(script_dir,'data/Solar/benchmark/d_cl_iso.csv')
out_cl_sim_d=os.path.join(script_dir,'data/Solar/benchmark/d_cl_sim.csv')
out_sin_iso_d=os.path.join(script_dir,'data/Solar/benchmark/d_sin_iso.csv')
out_sin_sim_d=os.path.join(script_dir,'data/Solar/benchmark/d_sin_sim.csv')
# outputs (monthly)
out_cl_iso_m=os.path.join(script_dir,'data/Solar/benchmark/m_cl_iso.csv')
out_cl_sim_m=os.path.join(script_dir,'data/Solar/benchmark/m_cl_sim.csv')
out_sin_iso_m=os.path.join(script_dir,'data/Solar/benchmark/m_sin_iso.csv')
out_sin_sim_m=os.path.join(script_dir,'data/Solar/benchmark/m_sin_sim.csv')

#
# pv_plot.py
#
# selection mode
# 1: Scatter chart
# 2: Time series for clusters and combined clusters
# 3: Selected time series (sin or any other asset)
# 4: Boxplot for the current technology
# 5: Boxplot for different subgrid configurations
# 6: Produce statistics from simulated and observed
modeplt=5
# dataset selection (1: BC; 2: NBC)
modebias=1
# exceptions (projects that do not have enough data)
exceptions=[]
# time series plot (mode 3)
asset='sin'         # sin or any other asset name
beg='2022-01-01'
end='2022-12-31'
# boxplot analysis
minpts=30
# define parameters to plot statistics
stat_min=36         # minimum number of data points  for the benchmark (36 months for at least 3 years)
# folder to save figures from the code
fig_path_folder=os.path.join(script_dir,'figures/plot/')
stats_path_folder=os.path.join(script_dir,'data/Solar/benchmark/stats/')
# hourly ISO
cl_iso_h=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_cl_iso.csv')
sin_iso_h=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_sin_iso.csv')
# hourly BC
bc_cl_sim_h=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_cl_sim.csv')
bc_sin_sim_h=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_sin_sim.csv')
# hourly BNC
nbc_cl_sim_h=os.path.join(script_dir,'data/Solar/benchmark/NBC/h_cl_sim.csv')
nbc_sin_sim_h=os.path.join(script_dir,'data/Solar/benchmark/NBC/h_sin_sim.csv')
# daily ISO
cl_iso_d=os.path.join(script_dir,'data/Solar/benchmark/BC_1/d_cl_iso.csv')
sin_iso_d=os.path.join(script_dir,'data/Solar/benchmark/BC_1/d_sin_iso.csv')
# daily BC
bc_cl_sim_d=os.path.join(script_dir,'data/Solar/benchmark/BC_1/d_cl_sim.csv')
bc_sin_sim_d=os.path.join(script_dir,'data/Solar/benchmark/BC_1/d_sin_sim.csv')
# daily NBC
nbc_cl_sim_d=os.path.join(script_dir,'data/Solar/benchmark/NBC/d_cl_sim.csv')
nbc_sin_sim_d=os.path.join(script_dir,'data/Solar/benchmark/NBC/d_sin_sim.csv')
# monthly ISO
cl_iso_m=os.path.join(script_dir,'data/Solar/benchmark/BC_1/m_cl_iso.csv')
sin_iso_m=os.path.join(script_dir,'data/Solar/benchmark/BC_1/m_sin_iso.csv')
# monthly BC
bc_cl_sim_m=os.path.join(script_dir,'data/Solar/benchmark/BC_1/m_cl_sim.csv')
bc_sin_sim_m=os.path.join(script_dir,'data/Solar/benchmark/BC_1/m_sin_sim.csv')
# monthly NBC
nbc_cl_sim_m=os.path.join(script_dir,'data/Solar/benchmark/NBC/m_cl_sim.csv')
nbc_sin_sim_m=os.path.join(script_dir,'data/Solar/benchmark/NBC/m_sin_sim.csv')