"""
Wind Parameters

Description:
This script sets the overall assumptions to perform simulations:

Usage:
1. Set all parameters below before running the simulation

For support or inquiries, contact emanueel@gmail.com

"""

# Set definitions for wind
import os

# get existing path
script_dir = os.path.dirname(__file__)  # absolute directory

# COMMON & DATES
yr_simcl1=1994           	# wind_simulation.py (cluster)
yr_simcl2=2023
yr_simex1=1994           	# wind_simulation.py (asset)
yr_simex2=2023
bc_act=True                 # wind_simulation.py (apply bias correction to actuals)
yr_gwa1=2008                # years used in GWA
yr_gwa2=2017
method_sel=2                # selection method of the subgrid - 1: average; 2: centroid
subdim=30                   # dimension subgrid NxN
nth_perc=0.75               # select the top 1-nth_perc  
# PATHS
era_src=os.path.join(script_dir,'../../data/ERA5/')
gwa_src_A_10=os.path.join(script_dir,'../../data/GWA/BRA_combined-Weibull-A_10m.tif')
gwa_src_A_100=os.path.join(script_dir,'../../data/GWA/BRA_combined-Weibull-A_100m.tif')
gwa_src_k_10=os.path.join(script_dir,'../../data/GWA/BRA_combined-Weibull-k_10m.tif')
gwa_src_k_100=os.path.join(script_dir,'../../data/GWA/BRA_combined-Weibull-k_100m.tif')
gwa_src_v_10=os.path.join(script_dir,'../../data/GWA/BRA_wind-speed_10m.tif')
gwa_src_v_100=os.path.join(script_dir,'../../data/GWA/BRA_wind-speed_100m.tif')
file_red_act=os.path.join(script_dir,'../../data/ACT/Wind/Asset/wind_assets.csv')           # unit-level assets
file_red_cls=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wind_subgrid_sim.csv')    # post-processed SVD clusters
folder_out_act=os.path.join(script_dir,'../../data/ACT/Wind/Simulated/')                    # store simulation at the unit-level for existing assets
folder_out_cls=os.path.join(script_dir,'../../data/CLS/Wind/Simulated/')                    # store simulation at the cluster level
svd_wind=os.path.join(script_dir,'../../data/CLS/Wind/SVD/wind_db_out.csv')
svd_wind_label=os.path.join(script_dir,'../../data/CLS/Wind/SVD/wind_db_out_label.csv')
wind_onwpc_file=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wpc_onshore.csv')      # save power curves
wind_offwpc_file=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wpc_offshore.csv')
