"""
Wind simulation

Description:
This script performs the following tasks:
1. Read data ERA5 and GWA to perform simulations for either existing assets or subgrids (clusters)
2. Perform simulations with and without bias corrected data using GWA
3. Store outputs at either the asset level or subgrid

Usage:
1. Configure the parameters of simulation using the wind_config.py
2. Run the current script (check the instances below)
3. Outputs will be stored in following folders (make sure folders are created):
   "/data/ACT/Wind/Simulated/" (actuals)
   "/data/CLS/Wind/Simulated/" (subgrid)

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
# current version
version = "1.0.0"

# importing libraries
import argparse
import time
import datetime
import os

# Importing custom classes and functions
from classes.wind_simulation import WindAsset
from classes.wind_simulation import WindMeshAsset
from classes.wind_simulation import WindCluster

# begin
def main():
    parser = argparse.ArgumentParser(description="Run simulations")
    parser.add_argument("--WindAsset", action="store_true", help="Run simulations for existing assets (with/without SD-BC) without subgrids")
    parser.add_argument("--WindMeshAsset", action="store_true", help="Run simulations for existing assets using subgrids")
    parser.add_argument("--WindCluster", action="store_true", help="Run simulations for clusters using subgrids")

    # parse the arguments
    args = parser.parse_args()

    # set up logging and time tracking
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.time()
    print("# Wind simulation")
    print(f"Current time: {current_time}")
    print(f"Version: {version}")
    print('#')

    # conditional execution based on arguments
    if args.WindAsset:
        windsim_asset = WindAsset()             # create instance of WindSimAsset (with/without BC using unit-level data)
        windsim_asset.run_simulation() 
    if args.WindMeshAsset:
        windsim_asset = WindMeshAsset()         # create instance of WindMesh for existing assets using subgrids
        windsim_asset.run_simulation()
    if args.WindCluster:
        windsim_cluster = WindCluster()         # create instance of WindMesh for clusters using subgridss
        windsim_cluster.run_simulation()

    # end of script runtime information
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()
