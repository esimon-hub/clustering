"""
Perform singular value decomposition for solar and wind resources 

Description:
This script performs the following tasks:
1. Read parameters defined on main_config.py
2. Read ERA5 climate variables
3. Perform singular value decomposition

Usage:
1. Configure the parameters of simulation using the main_config.py
2. Run the current script (check the instances below)
3. Outputs will be stored in within the following folders:
   "/data/SVD/Solar/"
   "/data/SVD/Wind/"

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
# current version
version = "1.0.0"

# importing libraries
import argparse
import time
import datetime

# importing custom classes and functions
from classes.svd import PerformSVD

# begin
def main():
    parser = argparse.ArgumentParser(description="Run SVD.")
    parser.add_argument("--preset", action="store_true", help="Preprocess data before running decomposition")
    parser.add_argument("--run_svd", action="store_true", help="Run Singular Value decomposition")
    parser.add_argument("--plot", action="store_true", help="Create plots with results")

    # parse the arguments
    args = parser.parse_args()

    # set up logging and time tracking
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    start_time = time.time()
    print("# Singular Value Decomposition")
    print(f"Current time: {current_time}")
    print(f"Version: {version}")
    print('#')

    # create instance of Clusters
    svd_perform = PerformSVD()

    # conditional execution based on arguments
    if args.preset:
        svd_perform.preset()
    if args.run_svd:
        svd_perform.run_svd()
    if args.plot:
        svd_perform.plot()

    # end of script runtime information
    print(f"Execution time: {time.time() - start_time} seconds")

if __name__ == "__main__":
    main()

