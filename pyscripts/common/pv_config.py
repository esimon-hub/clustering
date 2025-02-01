"""
Solar PV Parameters

Description:
This script sets the overall assumptions to perform simulations:

Usage:
1. Set all parameters below before running the simulation
2. Make sure nc files follow the pattern CCC_YYYY_MM.nc (CCC=country; YYYY=year; MM=month)

For support or inquiries, contact emanueel@gmail.com

"""

# Set definitions for solar PV
import os

# get existing path
script_dir = os.path.dirname(__file__)  # absolute directory

# COMMON & DATES
yr_simcl1=1994           	# pv_simulation.py (cluster)
yr_simcl2=2023
yr_simex1=1994           	# pv_simulation.py (asset)
yr_simex2=2023
bc_act=True                 # pv_simulation.py (apply bias correction to actuals)
yr_gsa1=1999                # years used in GSA
yr_gsa2=2018
method_sel=2                # selection method of the subgrid - 1: average; 2: centroid
subdim=10                   # dimension subgrid NxN
cl_tech=1                   # cluster technology - 0: fixed tilt; 1: single-axis tracker
nth_perc=0.75               # select the top 1-nth_perc         
# PATHS
era_src=os.path.join(script_dir,'../../data/ERA5/')
gsa_src=os.path.join(script_dir,'../../data/GSA/GHI.tif')
file_red_act=os.path.join(script_dir,'../../data/ACT/Solar/Asset/pv_assets.csv')            # project-level assets
file_red_cls=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/pv_subgrid_sim.csv')     # post-processed SVD clusters
folder_out_act=os.path.join(script_dir,'../../data/ACT/Solar/Simulated/')                   # store simulation at the unit-level for existing assets
folder_out_cls=os.path.join(script_dir,'../../data/CLS/Solar/Simulated/')                   # store simulation at the cluster level
svd_solar=os.path.join(script_dir,'../../data/CLS/Solar/SVD/solar_db_out.csv')
svd_solar_label=os.path.join(script_dir,'../../data/CLS/Solar/SVD/solar_db_out_label.csv')