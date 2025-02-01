"""
Create charts and tables for solar pv and wind resource assessment

Description:
This script performs the following tasks:
1. Read results from the cluster simulation
2. Adjust solar pv and wind generation profiles to the same format and timezone
3. Files are saved in the specified folders

Usage:
1. Configure the parameters of simulation using the main_config.py
2. Run the current script (check the instances below)
3. Outputs will be stored within the following folders:
   "/data/CLS/Solar/Cluster/"
   "/data/CLS/Wind/Cluster/"

For support or inquiries, contact emanueel@gmail.com

"""
# current version
version = "1.0.0"

# importing libraries
import argparse
import time
import datetime
import os

#   Importing custom classes and functions
from classes.resource_assessment import Resource

# begin
def main():
    parser = argparse.ArgumentParser(description="Run SVD.")
    parser.add_argument("--tables", action="store_true", help="Create table with geographical summary")
    parser.add_argument("--assessment", action="store_true", help="Perform geographical assessment for each cluster")
    parser.add_argument("--correl", action="store_true", help="Create a correlation matrix with all clusters")

    # parse the arguments
    args = parser.parse_args()

    # set up logging and time tracking
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.time()
    print("# Renewable resource assessment module")
    print(f"Current time: {current_time}")
    print(f"Version: {version}")
    print('#')

    # create instance of Resource
    res_assess = Resource()

    # conditional execution based on arguments
    if args.tables:
        res_assess.table_summary()
    if args.assessment:
        res_assess.resource_assessment()
    if args.correl:
        res_assess.corr()

    # end of script runtime information
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()














