"""
Create clusters from singular value decomposition and post-process information

Description:
This script performs the following tasks:
1. Read results from the singular value decomposition analysis
2. Create clusters using K-Means
3. Post-process clusters to eliminate unwanted clusters (environmental constraints or urban areas)

Usage:
1. Configure the parameters of simulation using the main_config.py
2. Run the current script (check the instances below)
3. Outputs will be stored within the following folders:
   "/data/CLS/Solar/Cluster/pv_cluster.csv"
   "/data/CLS/Wind/Cluster/wind_cluster.csv"

For support or inquiries, contact emanueel@gmail.com

"""
# current version
version = "1.0.0"

# importing libraries
import argparse
import time
import datetime
import os

# importing custom classes and functions
from classes.clusters import Clusters

# begin
def main():
    parser = argparse.ArgumentParser(description="Run clustering operations.")
    parser.add_argument("--calc_wcss", action="store_true", help="Calculate WCSS")
    parser.add_argument("--ksubplots", action="store_true", help="Plot clusters as a function of k")
    parser.add_argument("--singlek_pre", action="store_true", help="Create dataframe and plot map for a given k")
    parser.add_argument("--opassets", action="store_true", help="Plot a map with existing assets")
    parser.add_argument("--postproc", action="store_true", help="Create subgrids and eliminate exclusion zones")
    parser.add_argument("--seltop", action="store_true", help="Select the Nth best resources for each subgrid")
    parser.add_argument("--simsubgrid_view", action="store_true", help="Plot a map with simulated subgrids")
    parser.add_argument("--simgrid_view", action="store_true", help="Plot a map with clusters and labels")

    # parse the arguments
    args = parser.parse_args()

    # set up logging and time tracking
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.time()
    print("# Clustering module")
    print(f"Current time: {current_time}")
    print(f"Version: {version}")
    print('#')

    # create instance of Clusters
    cluster_vre = Clusters()

    # conditional execution based on arguments
    if args.calc_wcss:
        cluster_vre.calc_wcss()
    if args.ksubplots:
        cluster_vre.ksubplots()
    if args.singlek_pre:
        cluster_vre.singlek_pre()
    if args.opassets:
        cluster_vre.opassets()
    if args.postproc:
        cluster_vre.postproc()
    if args.seltop:
        cluster_vre.seltop()
    if args.simsubgrid_view:
        cluster_vre.simsubgrid_view()
    if args.simgrid_view:
        cluster_vre.simgrid_view()    

    # end of script runtime information
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()
