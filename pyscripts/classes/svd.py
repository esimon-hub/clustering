"""
SVD Simulation Class 

Description:
This class includes the methods needed to perform the SVD for solar and wind resources.

Usage:
Import this into the svd_exe.py to create instances and perform decomposition.

For support or inquiries, contact emanueel@gmail.com

"""
import xarray as xr
import glob
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import gc
import os
import time
from shapely.geometry import Point
from sklearn.preprocessing import MinMaxScaler
# This algorithm (randomized_svd) provides an approximation of the SVD, but it can be 
# significantly faster than the deterministic algorithm used by TruncatedSVD.     
from sklearn.utils.extmath import randomized_svd

#   Importing libraries and settings
from common.main_config import climatevar, minvals, shape, threshold                                            # general parameters
from common.main_config import yr_svd1, yr_svd2                                                                 # years to perform analysis
from common.main_config import shapeon, shapeonoff                                                              # shapefiles
from common.main_config import era_src                                                                          # ERA5 source
from common.main_config import solar_geocoord, solar_geoindex, solar_space, solar_time, solar_avrg, solar_sigma # solar files
from common.main_config import wind_geocoord, wind_geoindex, wind_space, wind_time, wind_avrg, wind_sigma       # wind files
from common.main_config import solar_svd_fig1, solar_svd_fig2                                                   # solar png
from common.main_config import solar_svd_pdf1, solar_svd_pdf2                                                   # solar pdf
from common.main_config import wind_svd_fig1, wind_svd_fig2                                                     # wind png
from common.main_config import wind_svd_pdf1, wind_svd_pdf2                                                     # wind pdf
from common.main_config import solar_geoindex_onshore, solar_geoindex_offshore                                  # dataframe with indexes
from common.main_config import wind_geoindex_onshore, wind_geoindex_offshore                                    # dataframe with indexes

# class to perform simulations for existing solar PV assets
class PerformSVD():
    def __init__(self):
        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        seltf = [f'{yr}{mn:02}' for yr in range(yr_svd1,yr_svd2+1) for mn in range(1,13)]
        file_path = []
        for pattern in seltf:
            file_path.append(glob.glob(era_src + '*_' + str(pattern) + '*.nc')[0])
        # check whether the content was loaded or not
        if not file_path:
            print("[!] Error. The list of file paths is empty. Stopping execution.")
            sys.exit(1)
        else:
            try:
                with xr.open_mfdataset(file_path, combine='nested', concat_dim='time', chunks={'lon': 200, 'lat':200}) as self.era5:
                    num_el=str(self.era5.ssrd.size)                             # grid points (lat.lon) x 8760
                    memory_usage_bytes = self.era5.ssrd.nbytes
                    memory_usage_megabytes = memory_usage_bytes / (1024**2)     # Convert bytes to megabytes
                    num_lat = self.era5.dims['latitude']
                    num_lon = self.era5.dims['longitude']
                    print(f"[+] Total number of elements (grid points x hours year): {num_el}")
                    print(f"[+] Total number of grid cells (original): {num_lat*num_lon}")
                    print(f"[+] Memory usage raw: {memory_usage_bytes}")
                    print(f"[+] Memory usage raw: {memory_usage_megabytes} MB")
                    print("[+] ERA5 WGS84")
                    # loaded ok
                    print("[+] ERA5 for BC successfully loaded.")
            except FileNotFoundError:
                # print error if exception occurs
                print(f"[!] Error. File not found. Stopping execution.")
                sys.exit(1)               
            except Exception as e:
                # print error if exception occurs
                print(f"[!] Error. Unable to open datasets. Stopping execution.")
                sys.exit(1)

    # method to preprocess the data before running decomposition
    def preset(self):
        # select the relevant data for each resource
        if climatevar==1:
            print("[+] SDV for solar resources")
            # solar surface radiation downwards
            xres=self.era5.ssrd.sel()
            var_size = xres.nbytes / (1024**3)      # give the size in GB
            print(f"[+] Memory size of the reanalysis variables (total): {var_size:.2f}")            
            poly_shape=shapeon                      # prefilter data using shapefiles (select onshore grid cells within the shapefiles)
            # output
            self.out_geocoords=solar_geocoord
            self.out_geoindex=solar_geoindex
            self.out_space=solar_space
            self.out_time=solar_time
            self.out_avrg=solar_avrg
            self.out_sigma=solar_sigma
            self.geoindex_onshore=solar_geoindex_onshore
            self.geoindex_offshore=solar_geoindex_offshore
        elif climatevar==2:
            print("[+] SDV for wind resources")
            # wind components u and v at 100 meters
            u100_data=self.era5.u100.sel()
            v100_data=self.era5.v100.sel()
            # calculate the absolute wind speed for each hour individually
            xres = ((u100_data**2+v100_data**2)**0.5).compute()
            var_size = xres.nbytes / (1024**3)      # give the size in GB
            print(f"[+] Memory size of the reanalysis variables (total): {var_size:.2f}")
            del u100_data
            del v100_data
            gc.collect()
            poly_shape=shapeonoff                   # prefilter data using shapefiles (select onshore grid cells within the shapefiles)
            # output
            self.out_geocoords=wind_geocoord
            self.out_geoindex=wind_geoindex
            self.out_space=wind_space
            self.out_time=wind_time
            self.out_avrg=wind_avrg
            self.out_sigma=wind_sigma
            self.geoindex_onshore=wind_geoindex_onshore
            self.geoindex_offshore=wind_geoindex_offshore
        else:
            print("[!] Climate variable not available.")
            quit()
        # close nc
        self.era5.close()
        # provide information about the size of dimensions and memory usage
        var_dim = xres.dims
        var_shape = xres.shape
        var_size = xres.nbytes / (1024**3)  # give the size in GB
        print("[+] Dimensions of the reanalysis variable: ", var_dim)
        print("[+] Size of dimensions of the reanalysis variable: ", var_shape)
        print(f"[+] Memory size of the reanalysis variable: {var_size:.2f}")
        # get dimensions of the xarray
        n_time, n_lat, n_lon = xres.shape
        # get the components of the xarray individually
        xsts = xres.coords['time'].values
        xlat = xres.coords['latitude'].values
        xlon = xres.coords['longitude'].values
        xs = xres.data.reshape(n_time, n_lat*n_lon)
        # eliminate unneeded variables before perfoming operations
        del xres
        gc.collect()

        # reduce the values to show just resources above a certain threshold
        if minvals:
            xs_avrg = np.mean(xs, axis=0)
            if climatevar==1:
                xs_avrg /= 3600000
                keep_index_1 = np.where(xs_avrg>0.1826)[0]  # 4.38 kWh/m2/d (1600kWh/8760h) 0.1826
            elif climatevar==2:
                keep_index_1 = np.where(xs_avrg>=5)[0]       # 5 m/s min to make sure values with high potential within the cells are captured
        # stack locations to perform the geo join (keep points within the borders)
        self.df_points = pd.DataFrame({
            'lat': np.repeat(xlat, len(xlon)),
            'lon': np.tile(xlon, len(xlat))
        })
        geometry = [Point(xy) for xy in zip(self.df_points['lon'], self.df_points['lat'])]
        # GeoDataFrame is endogenously modifying self.df_points by addint a new column "geometry"
        gdf_points = gpd.GeoDataFrame(self.df_points, geometry=geometry, crs='EPSG:4326')

        # perform the spatial merge to get the relevant climate variables
        poly_data=gpd.read_file(poly_shape, encoding='utf-8')
        geo_join = gpd.sjoin(gdf_points, poly_data, how='left', op='within')   # spatial join (points within a polygon)
        keep_index_2 = geo_join[geo_join['index_right'].notnull()].index
        # onshore
        poly_data_on=gpd.read_file(shapeon, encoding='utf-8')
        geo_join = gpd.sjoin(gdf_points, poly_data_on, how='left', op='within')   # spatial join (points within a polygon)
        onshore_index = geo_join[geo_join['index_right'].notnull()].index
        onshore_index = pd.DataFrame({'onshore_index': onshore_index})
        onshore_index.to_csv(self.geoindex_onshore, index=False)
        # offshore
        poly_data_onoff=gpd.read_file(shapeonoff, encoding='utf-8')
        geo_join = gpd.sjoin(gdf_points, poly_data_onoff, how='left', op='within')
        offshore_index = geo_join[geo_join['index_right'].notnull()].index
        offshore_index = pd.DataFrame({'offshore_index': offshore_index})
        offshore_index = offshore_index[~offshore_index['offshore_index'].isin(onshore_index['onshore_index'])]
        offshore_index.to_csv(self.geoindex_offshore, index=False)

        # combine indexes that are needed to perform the analysis (goal is to reduce memory size for SVD and also improve clusters)
        if minvals:
            keep_index = np.intersect1d(keep_index_1, keep_index_2.values)
        else:
            keep_index = keep_index_2

        # keep selected indexes
        self.keep_index = pd.DataFrame({'keep_index': keep_index})

        # keep just the selected elements and print the new size in memory
        self.xs = xs[:,keep_index]

        var_shape = self.xs.shape
        var_size = self.xs.nbytes / (1024**3)  # give the size in GB
        print("[+] Reduced memory dimensions: ", var_shape)
        print(f"[+] Reduced memory size of the reanalysis variable: {var_size:.2f}")

        # storage average resource values (select k percentile for clustering)
        res_avrg = np.mean(xs,axis=0)
        res_avrg = pd.DataFrame({'average': res_avrg})
        res_avrg.to_csv(self.out_avrg, index=False)

        del res_avrg

    # SVD
    def run_svd(self):
        # perform the SVD with just the shape of the resources (variance = 1)
        if shape:
            scaler = MinMaxScaler()
            xs = scaler.fit_transform(self.xs)

        # calculate singular value decomposition
        n_comp=min(xs.shape)
        start_time = time.time()
        U, sigma, VT = randomized_svd(xs, n_components=n_comp, random_state=None)
        print("[+] Total execution time (SVD): ", time.time()-start_time)
        S = np.diag(sigma)
        # calculate the cumuulative variance
        exp_var = (sigma ** 2) / np.sum(sigma ** 2)
        cum_exp_var = np.cumsum(exp_var)
        # truncate the matrix to meet the required threshold
        num_sing_vals = np.argmax(cum_exp_var>=threshold) + 1
        # print relevant information to the user and truncated the matrices
        print("[+] Number of singular values (truncated): ", num_sing_vals)
        U = U[:,:num_sing_vals]
        Msigma = np.diag(sigma[:num_sing_vals])
        VT = VT[:num_sing_vals,:]
        # create dataframe to store data
        print("[+] Storing sigma values from singular value decomposition")
        df_sigma = pd.DataFrame({
            'sigma': sigma,
            'cum_exp_var': cum_exp_var,
            'truncated': [num_sing_vals] + [np.nan] * (len(sigma) - 1), 
            'n_comp': [n_comp] + [np.nan] * (len(sigma) - 1)
        })
        df_sigma.to_csv(self.out_sigma, index=False)
        # calculate spatial and temporal weights
        Usigma = np.matmul(U,Msigma)
        sigmaVT = np.matmul(Msigma,VT)
        # create a dataframe for weights
        timeprod = pd.DataFrame(Usigma)
        spaceprod = pd.DataFrame(sigmaVT)
        # store the weighting factors and geographical coordinates to perform clustering
        timeprod.to_csv(self.out_time, index=False)
        spaceprod.to_csv(self.out_space, index=False)
        self.df_points.to_csv(self.out_geocoords, index=False)
        self.keep_index.to_csv(self.out_geoindex, index=False)

        # print execution time
        print("[+] Done.")
    
    def plot_xy(self, x, y_singular, y_cumulative, n_comp, title):
        # Initialize a figure with subplots
        fig, axs = plt.subplots(1, 2, figsize=(7, 3))  # Adjusted the overall figsize for better display
        plt.rcParams['font.size'] = 9
        plt.subplots_adjust(left=0.1, right=0.97, bottom=0.15, top=None, wspace=0.3, hspace=None)
        fig.suptitle(title, fontsize=10, fontweight='bold')

        # Singular values plot
        axs[0].semilogy(range(1, len(x) + 1), y_singular, 'ko', markersize=3, markerfacecolor='black', markeredgewidth=0)
        axs[0].set_xlabel('d', fontsize=9)
        axs[0].set_ylabel('Singular value (\u03C3)', fontsize=9)
        axs[0].grid(True, linestyle='--', linewidth=0.5, color='grey')

        # Cumulative energy plot
        rel_error = 1 - np.array(y_cumulative)
        axs[1].plot(range(1, len(x) + 1), rel_error, 'ko', markersize=3, markerfacecolor='black', markeredgewidth=0)
        axs[1].set_xlabel('d', fontsize=9)
        axs[1].set_ylabel('Relative error (\u03B5)', fontsize=9)
        axs[1].grid(True, linestyle='--', linewidth=0.5, color='grey')

        plt.savefig(self.out_fig_1, dpi=1000, format='jpg')  # Save the figure as a single combined image
        plt.savefig(self.out_pdf_1, format='pdf')

        # Annotate specific components if n_comp > 500
        if n_comp > 500:
            high_comps = np.array([1, 41])
            singular_high_vals = y_singular[high_comps - 1]
            error_high_vals = rel_error[high_comps - 1]

            # singular values annotations
            axs[0].scatter(high_comps, singular_high_vals, color='red', s=3, zorder=3)
            for comp, value in zip(high_comps, singular_high_vals):
                xan = 7 if comp < 50 else 7
                axs[0].annotate(f'd={comp}', (comp, value), textcoords="offset points", xytext=(xan, -2),
                                ha='left', fontsize=10, color='red')
            
            # relative error annotations
            axs[1].scatter(high_comps, error_high_vals, color='red', s=3, zorder=3)
            for comp, value in zip(high_comps, error_high_vals):
                xan = 7 if comp < 50 else 7
                axs[1].annotate(f'd={comp}', (comp, value), textcoords="offset points", xytext=(xan, -2),
                                ha='left', fontsize=10, color='red')

        plt.savefig(self.out_fig_2, dpi=1000, format='jpg')  # Save the figure as a single combined image
        plt.savefig(self.out_pdf_2, format='pdf')


    # function to plot figures (singular values and cumulative values)
    def plot(self):
        # create dataframe to store data
        print("[+] Ploting sigma values from singular value decomposition")
        # select the relevant data for each resource
        if climatevar==1:
            print("[+] Plot solar singular values")
            self.out_sigma=solar_sigma
            self.out_fig_1=solar_svd_fig1
            self.out_pdf_1=solar_svd_pdf1
            self.out_fig_2=solar_svd_fig2
            self.out_pdf_2=solar_svd_pdf2
            title="Solar radiation"
        elif climatevar==2:
            print("[+] Plot wind singular values")
            self.out_sigma=wind_sigma
            self.out_fig_1=wind_svd_fig1
            self.out_pdf_1=wind_svd_pdf1
            self.out_fig_2=wind_svd_fig2
            self.out_pdf_2=wind_svd_pdf2
            title="Wind speeds"
        else:
            print("[!] Climate variable not available.")
            quit()
        # read stored csv
        df_sigma = pd.read_csv(self.out_sigma)
        # Load metadata from the first row
        sigma = df_sigma['sigma']
        cum_exp_var = df_sigma['cum_exp_var']
        n_comp = int(df_sigma['n_comp'].dropna().iloc[0])  # Assuming scalar metadata is stored in the first row
        # plot singular values and cumulative energy
        self.plot_xy(sigma, sigma, cum_exp_var, n_comp, title)