"""
Solar PV simulation 

Description:
This script performs the following tasks:
1. Read data ERA5 and GSA to perform simulations for either existing assets or subgrids
2. Perform simulations with and without bias corrected data from GSA
3. Store outputs at either the asset level or subgrid

Usage:
1. Configure the parameters of simulation using the pv_config.py
2. Run the current script (check the instances below)
3. Outputs will be stored in following folders (make sure folders are created):
   "/data/ACT/Solar/Simulated/" (actuals)
   "/data/CLS/Solar/Simulated/" (subgrid)

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
from classes.pv_simulation import PVAsset
from classes.pv_simulation import PVMeshAsset
from classes.pv_simulation import PVCluster

# begin
def main():
    parser = argparse.ArgumentParser(description="Run simulations")
    parser.add_argument("--PVAsset", action="store_true", help="Run simulations for existing assets (with/without SD-BC) without subgrids")
    parser.add_argument("--PVMeshAsset", action="store_true", help="Run simulations for existing assets using subgrids")
    parser.add_argument("--PVCluster", action="store_true", help="Run simulations for clusters using subgrids")

    # parse the arguments
    args = parser.parse_args()

    # set up logging and time tracking
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.time()
    print("# Solar PV simulation")
    print(f"Current time: {current_time}")
    print(f"Version: {version}")
    print('#')

    # conditional execution based on arguments
    if args.PVAsset:
        pvsim_asset = PVAsset()           # create instance of PVAsset (with/without BC using project-level data)
        pvsim_asset.run_simulation() 
    if args.PVMeshAsset:
        pvsim_asset = PVMeshAsset()       # create instance of PVMeshAsset for existing assets using subgrids
        pvsim_asset.run_simulation()
    if args.PVCluster:
        pvsim_cluster = PVCluster()       # # create instance of PVCluster for clusters using subgrids
        pvsim_cluster.run_simulation()

    # end of script runtime information
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()
