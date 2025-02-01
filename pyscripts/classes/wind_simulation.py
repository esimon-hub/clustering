"""
Wind Simulation Class 

Description:
This class includes the methods needed to perform the simulation of clusters and existing wind assets.

Usage:
Import this into the wind_simulation_exec.py to create instances and perform simulations.

For support or inquiries, contact emanueel@gmail.com

"""
# Importing libraries and settings
from common.wind_config import era_src                                                                                      # common settings
from common.wind_config import yr_simex1, yr_simex2, bc_act, file_red_act, folder_out_act                                   # actual settings
from common.wind_config import yr_simcl1, yr_simcl2, method_sel, file_red_cls, folder_out_cls                               # cluster settings
from common.wind_config import gwa_src_A_10, gwa_src_k_10, gwa_src_v_10                                                     # GWA files (10)
from common.wind_config import gwa_src_A_100, gwa_src_k_100, gwa_src_v_100                                                  # GWA files (100)
from common.wind_config import yr_gwa1, yr_gwa2, subdim                                                                     # GWA years and dimension
from common.wind_config import wind_onwpc_file, wind_offwpc_file                                                            # power curves for clusters 
from scipy.special import gamma
from scipy.spatial import cKDTree
from datetime import datetime
from multiprocessing import Pool
from rasterio.windows import Window
from .wind_shared import calc_on_cf, calc_of_cf                                                                             # parametric simulation
from .wind_shared import calc_cf, wpc_onshore, wpc_offshore                                                                 # power curves for clusters 
import xarray as xr
import glob
import sys
import os
import pandas as pd
import numpy as np
import rasterio
import gc           # garbage collector to release unreferenced memory back to OS
import psutil       # measure memory usage

# class to perform simulations for existing wind assets
class WindAsset():
    def __init__(self):
        print("[+] Begin")
        # read turbine information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and wind turbine data')
            self.turbine_level = pd.read_csv(file_red_act, index_col=False)
            # --- GeoTIFF Factor A, k, v at 100-meter high ---
            self.gwa_A10 = rasterio.open(gwa_src_A_10,'r')
            self.gwa_A100 = rasterio.open(gwa_src_A_100,'r')
            self.gwa_k10 = rasterio.open(gwa_src_k_10,'r')
            self.gwa_k100 = rasterio.open(gwa_src_k_100,'r')
            #self.gwa_v10 = rasterio.open(gwa_src_v_10,'r')
            #self.gwa_v100 = rasterio.open(gwa_src_v_100,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        # plot tiff info
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gwa_A10: {self.gwa_A10.width} x {self.gwa_A10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A10: {self.gwa_A10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_A100: {self.gwa_A100.width} x {self.gwa_A100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A100: {self.gwa_A100.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k10: {self.gwa_k10.width} x {self.gwa_k10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k10: {self.gwa_k10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k100: {self.gwa_k100.width} x {self.gwa_k100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k100: {self.gwa_k100.bounds}")
        print("[+] GeoTIFF files successfully opened")
        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simex1,yr_gwa1)
        yr_max=max(yr_simex2,yr_gwa2)
        seltf = [f'{yr}{mn:02}' for yr in range(yr_min,yr_max+1) for mn in range(1,13)]
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
                    print(f"[+] Memory usage raw: {memory_usage_megabytes} MB")
                    print("[+] ERA5 for BC successfully loaded.")
            except FileNotFoundError:
                # print error if exception occurs
                print(f"[!] Error. File not found. Stopping execution.")
                sys.exit(1)               
            except Exception as e:
                # print error if exception occurs
                print(f"[!] Error. Unable to open datasets. Stopping execution.")
                sys.exit(1)

    # check whether coordinates are out of bound or not
    def check_bounds(self, py, px, width, height):
        if not (0 <= px <= width and 0 <= py <= height):
            print(['[!] Error. Coordinates are out of bounds!'])
            sys.exit(1)

    # method to extract the value from a GeoTIFF file at the specified geographical coordinate
    def get_val(self, lon, lat, param_A, param_k):
        # Convert the geographical coordinates to dataset's pixel coordinates
        Apy, Apx = param_A.index(lon, lat)
        self.check_bounds(Apy, Apx, param_A.width, param_A.height)
        kpy, kpx = param_k.index(lon, lat)
        self.check_bounds(kpy, kpx, param_k.width, param_k.height)
        #vpy, vpx = param_v.index(lon, lat)
        #self.check_bounds(vpy, vpx, param_v.width, param_v.height)
        # window defines the data window we want to read (in this case, a single pixel)
        windowA = rasterio.windows.Window(Apx, Apy, 1, 1)
        windowk = rasterio.windows.Window(kpx, kpy, 1, 1)
        #windowv = rasterio.windows.Window(vpx, vpy, 1, 1)
        # Read the data in the window, assuming the raster has only one band
        bandA = param_A.read(1, window=windowA)
        bandk = param_k.read(1, window=windowk)
        #bandv = param_v.read(1, window=windowv)
        # Extract the value from the array
        A = bandA[0, 0]
        k = bandk[0, 0]
        #v = bandv[0, 0]
        # return values
        return A, k

    # execute simulation for existing assets
    def run_simulation(self):
        # select data for existing assets
        col_names = self.turbine_level['ID'].to_list()
        latlon=self.turbine_level[['lat','lon']].values

        # *************** Selection for GWA timeframe *************** 
        # retrieve weibull parameters from GWA to perform simulation with BC
        if bc_act:
            # calculate weibull parameters using selected ERA5 timeframe
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Calculating Weibull parameters for grid-cells at {current_time}")               
            era5 = self.era5.sel(time=slice(f'{yr_gwa1}-01-01',f'{yr_gwa2}-12-31'))

            # calculate absolute wind speeds (10m)
            sel_wnd_u10 = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v10 = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
            sel_wnd_10.name = 'ws10'
            df_wnd10_m = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_names)
            # calculate absolute wind speeds (100m)
            sel_wnd_u100 = era5.u100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                    longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v100 = era5.v100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
            sel_wnd_100.name = 'ws100'
            df_wnd100_m = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_names)
            # calculate mean, std deviation, A, and k for each location (ERA5)
            wnd10_m_mean = df_wnd10_m.mean()
            wnd10_m_std = df_wnd10_m.std()
            k10m = (wnd10_m_std / wnd10_m_mean) ** (-1.086)
            A10m = wnd10_m_mean / gamma(1+1/k10m)
            wnd100_m_mean = df_wnd100_m.mean()
            wnd100_m_std = df_wnd100_m.std()
            k100m = (wnd100_m_std / wnd100_m_mean) ** (-1.086)
            A100m = wnd100_m_mean / gamma(1+1/k100m)
            
            # select Weibull parameters using GWA
            A10_arr = []
            A100_arr = []
            k10_arr = []
            k100_arr = []
            #v10set = []
            #v100set = []
            # select the coordinates for the current project 
            for _, row in self.turbine_level.iterrows():
                # select data for the current project
                asset_id = row['ID']
                lat, lon = row['lat'], row['lon']
                print(f"[+] Retrieving parameters for asset: {asset_id}")
                # select Weibull parameters (A and k) from GWA
                Au10_val, ku10_val = self.get_val(lon,lat,self.gwa_A10,self.gwa_k10)
                Au100_val, ku100_val = self.get_val(lon,lat,self.gwa_A100,self.gwa_k100)
                # record parameters
                A10_arr.append(Au10_val)
                A100_arr.append(Au100_val)
                k10_arr.append(ku10_val)
                k100_arr.append(ku100_val)
            # create a dataframe that includes the name of the assets and the new values
            A10u = pd.Series(A10_arr, index=col_names)
            A100u = pd.Series(A100_arr, index=col_names)
            k10u = pd.Series(k10_arr, index=col_names)
            k100u = pd.Series(k100_arr, index=col_names)
            # print status of the dataframe creation
            memuse = psutil.Process(os.getpid())
            print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
            # store dataframe for future simulations
            #df_bc.to_csv(file_bcf_act, index=False, encoding="utf-8-sig")
            # delete temporary variables
            del sel_wnd_u10, sel_wnd_v10, sel_wnd_10, df_wnd10_m
            del sel_wnd_u100, sel_wnd_v100, sel_wnd_100, df_wnd100_m
            gc.collect()

        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for existing assets at {current_time}")        
        # execute simulations with or without downscaling
        for year in range(yr_simex1,yr_simex2+1):
            # select current year before the nearest method
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # calculate absolute wind speeds (10m)
            sel_wnd_u10 = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v10 = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
            sel_wnd_10.name = 'ws10'
            wnd10m = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_names)
            # calculate absolute wind speeds (100m)
            sel_wnd_u100 = era5.u100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v100 = era5.v100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
            sel_wnd_100.name = 'ws100'
            wnd100m = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_names)
            # downscale datasets
            if bc_act:
                wnd_10 = A10u*((wnd10m/A10m)**(k10m/k10u))
                wnd_100 = A100u*((wnd100m/A100m)**(k100m/k100u))
                alpha = np.log(wnd_100 / wnd_10) / np.log(100 / 10)
            else:
                wnd_10 = wnd10m
                wnd_100 = wnd100m
                alpha = np.log(wnd_100 / wnd_10) / np.log(100 / 10)
            # prepare for multithreading
            date_index = wnd_10.index
            subtasks = []
            # iterate over the dataframe to obtain the time series
            for _, row in self.turbine_level.iterrows():
                id = row['ID']
                sp = row['SPOW']
                hheight = row['ALT_TORRE']
                sel_wnd100 = wnd_100[id]
                sel_alpha = alpha[id]
                wnd_adj = (sel_wnd100 * (hheight/100)**sel_alpha)
                wnd_adj = np.where(wnd_adj<=25,wnd_adj,0)
                # create tuple to perform multithreading
                asset_tuple = (wnd_adj,sp)
                subtasks.append(asset_tuple)
            # plot status before multithreading
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Running simulation for year {year} with multithreading at {current_time}")
            # run multithreading with tasks
            with Pool(processes=4) as pool:
                # run calc_on_cf to each set of tasks
                results_yr = pool.starmap(calc_on_cf,subtasks)
            # plot status (completed)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Storing simulation for year {year} at {current_time}")
            # create a new dataframe for each year and store into a csv
            df_results_yr = pd.DataFrame(results_yr).T
            df_results_yr.columns = col_names
            df_results_yr.index = date_index
            df_results_yr.index.name = 'datetime'
            df_results_yr = df_results_yr.round(4)      # reduce size in disk
            csv_filename_cf = f'{folder_out_act}wind_cf_{year}.csv'
            df_results_yr.to_csv(csv_filename_cf, encoding="utf-8-sig")
        # print status
        df_capacity = pd.DataFrame(self.turbine_level.set_index('ID')['POT_MW']).transpose()
        csv_filename_cap=f'{folder_out_act}wind_capacity.csv'
        df_capacity.to_csv(csv_filename_cap, index=False, encoding="utf-8-sig")
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')

#
# class to perform simulations for existing wind assets using mesh
#
class WindMeshAsset():
    def __init__(self):
        print("[+] Begin")
        # read turbine information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and wind turbine data')
            self.turbine_level = pd.read_csv(file_red_act, index_col=False)
            # --- GeoTIFF Factor A, k, v at 100-meter high ---
            self.gwa_A10 = rasterio.open(gwa_src_A_10,'r')
            self.gwa_A100 = rasterio.open(gwa_src_A_100,'r')
            self.gwa_k10 = rasterio.open(gwa_src_k_10,'r')
            self.gwa_k100 = rasterio.open(gwa_src_k_100,'r')
            #self.gwa_v10 = rasterio.open(gwa_src_v_10,'r')
            #self.gwa_v100 = rasterio.open(gwa_src_v_100,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        # plot tiff info
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gwa_A10: {self.gwa_A10.width} x {self.gwa_A10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A10: {self.gwa_A10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_A100: {self.gwa_A100.width} x {self.gwa_A100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A100: {self.gwa_A100.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k10: {self.gwa_k10.width} x {self.gwa_k10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k10: {self.gwa_k10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k100: {self.gwa_k100.width} x {self.gwa_k100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k100: {self.gwa_k100.bounds}")
        print("[+] GeoTIFF files successfully opened")
        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simex1,yr_gwa1)
        yr_max=max(yr_simex2,yr_gwa2)
        seltf = [f'{yr}{mn:02}' for yr in range(yr_min,yr_max+1) for mn in range(1,13)]
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
                    print(f"[+] Memory usage raw: {memory_usage_megabytes} MB")
                    print("[+] ERA5 for BC successfully loaded.")
            except FileNotFoundError:
                # print error if exception occurs
                print(f"[!] Error. File not found. Stopping execution.")
                sys.exit(1)               
            except Exception as e:
                # print error if exception occurs
                print(f"[!] Error. Unable to open datasets. Stopping execution.")
                sys.exit(1)

    # create a mesh with the subgrids for all locations
    def mesh_generation(self, subdim):
        # ERA5 grid size is fixed at .25 degree for both latitude and longitude
        era5_size = 0.25
        # calculate the size of each subgrid based on the subdim
        subsize = era5_size / subdim
        # structure to store the mesh information
        mesh_info = []
        # create an extended array to allow for the geographical boundaries
        ext_lat = [self.era5.latitude.values[0] + era5_size] + \
            list(self.era5.latitude.values) + \
                [self.era5.latitude.values[-1] - era5_size]
        ext_lon = [self.era5.longitude.values[0] - era5_size] + \
            list(self.era5.longitude.values) + \
                [self.era5.longitude.values[-1] + era5_size]   
        for lat in ext_lat:
            for lon in ext_lon:
                # generate the subgrid for the current cell
                for i in range(subdim):
                    for j in range(subdim):
                        # calculate the center (centroid) of each subgrid cell
                        sub_center_lat = lat - (i + 0.5) * subsize
                        sub_center_lon = lon + (j + 0.5) * subsize
                        # store the information in the mesh_info list
                        mesh_info.append({
                            'Grid lat': lat,
                            'Grid lon': lon,
                            'Sub lat': sub_center_lat,  # centroid
                            'Sub lon': sub_center_lon,  # centroid
                        })
        # convert the mesh list into a dataframe
        mesh_df = pd.DataFrame(mesh_info)
        return mesh_df

    # function to retrieve the geographical limits (coordinates) for each subgrid
    def get_bounds(self, slat, slon, sdim):
        # calculate bounds of the subgrid
        subdelta = 0.25 / (sdim*2)
        min_lat = slat - subdelta
        max_lat = slat + subdelta
        min_lon = slon - subdelta
        max_lon = slon + subdelta
        # returns a dictionary instead of individual values
        bounds = {
            'min_lat': min_lat,
            'max_lat': max_lat,
            'min_lon': min_lon,
            'max_lon': max_lon
        }
        return bounds

    # check whether coordinates are out of bound or not
    def check_bounds(self, py, px, width, height):
        if not (0 <= px <= width and 0 <= py <= height):
            print(['[!] Error. Coordinates are out of bounds!'])
            sys.exit(1)

    # convert bounds to geotiff pixel coordinates for each of A, k, and v
    def get_geo_zone(self, bounds, geotiff):
        min_lat = bounds['min_lat']
        max_lat = bounds['max_lat']
        min_lon = bounds['min_lon']
        max_lon = bounds['max_lon']
        inv_transform = ~geotiff.transform
        left, top = inv_transform * (min_lon, max_lat)
        right, bottom = inv_transform * (max_lon, min_lat)
        left, bottom, right, top = map(int, [left, bottom, right, top])
        window = Window.from_slices((top, bottom), (left, right))
        # read the data within the window
        data = geotiff.read(1, window=window)
        # filter out any nodata values
        nodata = geotiff.nodata
        if nodata is not None:
            data = data[data != nodata]
        return data

    # method to extract the value from a GeoTIFF file at the specified geographical coordinate
    def get_geo_centroid(self, slon, slat, geotiff):
        # Convert the geographical coordinates to dataset's pixel coordinates
        py, px = geotiff.index(slon, slat)
        self.check_bounds(py, px, geotiff.width, geotiff.height)
        # window defines the data window we want to read (in this case, a single pixel)
        window_val = rasterio.windows.Window(px, py, 1, 1)
        # Read the data in the window, assuming the raster has only one band
        data_val = geotiff.read(1, window=window_val)
        # Extract the value from the array
        data = data_val[0,0]
        # return single element list to keep consistency with get_geo_zone
        return data
    
    # calculate the average values using downscaling
    def down_grid_zone(self, up_wnd, Am, Ak, km, kk):
        # perform downscaling using the given parameters (m - macro; u - micro)
        nk = len(Ak)
        wnd_k=0
        for ik in range(nk):
            Au = Ak[ik]
            ku = kk[ik]
            wnd_k += Au*((up_wnd/Am)**(km/ku))
        # return the average time series
        wnd_k = wnd_k / nk
        return wnd_k

    # calculate the average values using downscaling
    def down_grid_centroid(self, up_wnd, Am, Ak, km, kk):
        # perform downscaling using the given parameters (m - macro; u - micro)
        wnd_k = Ak*((up_wnd/Am)**(km/kk))
        return wnd_k

    # execute simulation for existing assets
    def run_simulation(self):
        # Grid-level data
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Select nearest ERA5 grid-cells for existing assets at {current_time}")
        era5_lat = []
        era5_lon = []
        for _, row in self.turbine_level.iterrows():
            era5_loc = self.era5.sel(latitude=row['lat'], longitude=row['lon'], method='nearest')
            n_lat = era5_loc.latitude.values.item()
            n_lon = era5_loc.longitude.values.item()
            era5_lat.append(n_lat)
            era5_lon.append(n_lon)
        self.turbine_level['Grid lat'] = era5_lat
        self.turbine_level['Grid lon'] = era5_lon

        #  Subgrid-level data
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f'[+] Creating meshes for subgrids at {current_time}')
        submesh = self.mesh_generation(subdim)
        sub_points = submesh[['Sub lat','Sub lon']].values
        turbine_points = self.turbine_level[['lat','lon']].values
        tree = cKDTree(sub_points)                              # create a KDTree for efficient nearest neighbor search
        distances, indices = tree.query(turbine_points)         # query the nearest subgrid for each wind turbine
        self.turbine_level['Sub lat'] = submesh.iloc[indices]['Sub lat'].values
        self.turbine_level['Sub lon'] = submesh.iloc[indices]['Sub lon'].values

        # Drop duplicates and create specific dataframes (up ERA5, down GWA) to improve the efficiency of the code
        sel_up = self.turbine_level[['Grid lat','Grid lon']].drop_duplicates().reset_index(drop=True)
        sel_up['id'] = range(len(sel_up))
        up_map = {(lat,lon): idx for lat, lon, idx in sel_up[['Grid lat','Grid lon','id']].to_numpy()}
        self.turbine_level['id up'] = self.turbine_level.apply(lambda row: up_map[(row['Grid lat'],row['Grid lon'])], axis=1)
        sel_down = self.turbine_level[['Grid lat','Grid lon','Sub lat','Sub lon']].drop_duplicates().reset_index(drop=True)
        sel_down['id'] = range(len(sel_down))
        down_map = {(lat,lon,slat,slon): idx for lat, lon, slat, slon, idx in sel_down[['Grid lat','Grid lon','Sub lat','Sub lon','id']].to_numpy()}
        self.turbine_level['id down'] = self.turbine_level.apply(lambda row: down_map[(row['Grid lat'],row['Grid lon'],row['Sub lat'],row['Sub lon'])], axis=1)

        # *************** Selection for GWA timeframe *************** 
        col_up = sel_up['id'].to_list()                             # ERA5 columns
        col_down = sel_down['id'].to_list()                         # Subgrid columns
        latlon_up=sel_up[['Grid lat','Grid lon']].values            # ERA5 geographical coordinates
        # retrieve Weibull parameters from GWA from 'sel_down' to perform simulations
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating Weibull parameters for grid-cells at {current_time}")
        era5 = self.era5.sel(time=slice(f'{yr_gwa1}-01-01',f'{yr_gwa2}-12-31'))
        # calculate absolute wind speeds (10m)
        sel_wnd_u10 = era5.u10.sel(
            latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
            longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
            method='nearest').compute()
        sel_wnd_v10 = era5.v10.sel(
            latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
            longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
            method='nearest').compute()
        sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
        sel_wnd_10.name = 'ws10'
        df_wnd10_m = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_up)
        # calculate absolute wind speeds (100m)
        sel_wnd_u100 = era5.u100.sel(
            latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
            longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
            method='nearest').compute()
        sel_wnd_v100 = era5.v100.sel(
            latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
            longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
            method='nearest').compute()
        sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
        sel_wnd_100.name = 'ws100'
        df_wnd100_m = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_up)
        # calculate mean, std deviation, A, and k for each location (ERA5)
        wnd10_m_mean = df_wnd10_m.mean()
        wnd10_m_std = df_wnd10_m.std()
        wnd100_m_mean = df_wnd100_m.mean()
        wnd100_m_std = df_wnd100_m.std()
        # calculate the shape (A) and scale (k) parameters for each location
        k10 = (wnd10_m_std / wnd10_m_mean) ** (-1.086)
        A10 = wnd10_m_mean / gamma(1+1/k10)
        k100 = (wnd100_m_std / wnd100_m_mean) ** (-1.086)
        A100 = wnd100_m_mean / gamma(1+1/k100)
        # print status of the dataframe creation
        memuse = psutil.Process(os.getpid())
        print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
        # store dataframe for future simulations
        #df_bc.to_csv(file_bcf_act, index=False, encoding="utf-8-sig")
        # delete temporary variables
        del sel_wnd_u10, sel_wnd_v10, sel_wnd_10, df_wnd10_m
        del sel_wnd_u100, sel_wnd_v100, sel_wnd_100, df_wnd100_m
        gc.collect()
        # create a dictionary with GWA parameters
        k10_dict = {}
        A10_dict = {}
        k100_dict = {}
        A100_dict = {}
        # iterate through each row in the DataFrame
        for _, row in sel_down.iterrows():
            # select data for the current subgrid
            sid = row['id']
            slat, slon = row['Sub lat'], row['Sub lon']
            print(f"[+] Retrieving parameters for lat{slat} and lon{slon}")
            # methods to downscale are average and centroid
            if method_sel==1:
                # calculate bounds to retrieve local data
                bounds = self.get_bounds(slat, slon, subdim)

                # select Weibull parameters (A and k) from GWA (returns a list)
                kk10 = self.get_geo_zone(bounds, self.gwa_k10)
                Ak10 = self.get_geo_zone(bounds, self.gwa_A10)
                kk100 = self.get_geo_zone(bounds, self.gwa_k100)
                Ak100 = self.get_geo_zone(bounds, self.gwa_A100)

                # store outputs in dictionaries
                k10_dict[sid] = kk10
                A10_dict[sid] = Ak10
                k100_dict[sid] = kk100
                A100_dict[sid] = Ak100
            elif method_sel==2:
                # select Weibull parameters (A and k) from GWA (returns a list)
                kk10 = self.get_geo_centroid(slon, slat, self.gwa_k10)
                Ak10 = self.get_geo_centroid(slon, slat, self.gwa_A10)
                kk100 = self.get_geo_centroid(slon, slat, self.gwa_k100)
                Ak100 = self.get_geo_centroid(slon, slat, self.gwa_A100)
                # store outputs in dictionaries
                k10_dict[sid] = kk10
                A10_dict[sid] = Ak10
                k100_dict[sid] = kk100
                A100_dict[sid] = Ak100

        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for existing assets at {current_time}")        
        # execute simulations with or without downscaling
        for year in range(yr_simex1,yr_simex2+1):
            # select current year before the nearest method
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # calculate absolute wind speeds (10m)
            sel_wnd_u10 = era5.u10.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v10 = era5.v10.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
            sel_wnd_10.name = 'ws10'
            wnd10m = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_up)
            # calculate absolute wind speeds (100m)
            sel_wnd_u100 = era5.u100.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v100 = era5.v100.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
            sel_wnd_100.name = 'ws100'
            wnd100m = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_up)
            # create dataframe to store downscale
            wnd10_sub = pd.DataFrame(index=sel_wnd_10.time.values, columns=col_down)
            wnd100_sub = pd.DataFrame(index=sel_wnd_100.time.values, columns=col_down)
            shearu_sub = pd.DataFrame(index=sel_wnd_100.time.values, columns=col_down)
            # downscale data for the projects
            for _, row in sel_down.iterrows():
                # get column for ERA5 (up) and GWA (down)
                id = row['id']
                lat = row['Grid lat']
                lon = row['Grid lon']
                slat = row['Sub lat']
                slon = row['Sub lon']
                id_up = up_map[(lat,lon)]
                id_down = down_map[(lat,lon,slat,slon)]
                # fetch wind speed and Weibull parameters for the corresponding ERA5 grid cell
                up_wnd10 = wnd10m.loc[:,id_up]
                Am10 = A10[id_up]
                km10 = k10[id_up]
                up_wnd100 = wnd100m.loc[:,id_up]
                Am100 = A100[id_up]
                km100 = k100[id_up]
                # downscale using either zone or centroid
                if method_sel==1:
                    # fetch Weibull parameters for the corresponding subgrid
                    Ak = A10_dict[id]
                    kk = k10_dict[id]
                    down_wnd10 = self.down_grid_zone(up_wnd10, Am10, Ak, km10, kk)
                    Ak = A100_dict[id]
                    kk = k100_dict[id]
                    down_wnd100 = self.down_grid_zone(up_wnd100, Am100, Ak, km100, kk)
                    down_shear = np.log(down_wnd100 / down_wnd10) / np.log(100 / 10)
                elif method_sel==2:
                    # fetch Weibull parameters for the centroid
                    Ak = A10_dict[id]
                    kk = k10_dict[id]
                    down_wnd10 = self.down_grid_centroid(up_wnd10, Am10, Ak, km10, kk)
                    Ak = A100_dict[id]
                    kk = k100_dict[id]
                    down_wnd100 = self.down_grid_centroid(up_wnd100, Am100, Ak, km100, kk)
                    down_shear = np.log(down_wnd100 / down_wnd10) / np.log(100 / 10)
                # store outputs using dataframes
                wnd10_sub[id_down] = down_wnd10
                wnd100_sub[id_down] = down_wnd100
                shearu_sub[id_down] = down_shear
            # Running simulations
            subtasks = []
            for _, row in self.turbine_level.iterrows():
                id = row['ID']
                sp = row['SPOW']
                hheight = row['ALT_TORRE']
                id_down = row['id down']
                # create the wind series for the current asset
                wnd_adj = wnd100_sub[id_down]
                alpha = shearu_sub[id_down]
                wnd_adj = (wnd_adj * (hheight/100)**alpha)
                wnd_adj = np.where(wnd_adj<=25,wnd_adj,0)
                # create tuple to perform multithreading
                asset_tuple = (wnd_adj,sp)
                subtasks.append(asset_tuple)
            # plot status before multithreading
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Running simulation for year {year} with multithreading at {current_time}")
            # run multithreading with tasks
            with Pool(processes=4) as pool:
                # run calc_on_cf to each set of tasks
                res_pool_yr = pool.starmap(calc_on_cf,subtasks)
            # plot status (completed)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Storing simulation for year {year} at {current_time}")
            # create a new dataframe for each year and store into a csv
            results = pd.DataFrame(res_pool_yr).T
            results.columns = self.turbine_level['ID']
            results.index = sel_wnd_100.time.values
            results.index.name = 'datetime'
            results = results.round(4)      # reduce size in disk
            csv_filename_cf = f'{folder_out_act}wind_cf_{year}.csv'
            results.to_csv(csv_filename_cf, encoding="utf-8-sig")
        # print status
        df_capacity = pd.DataFrame(self.turbine_level.set_index('ID')['POT_MW']).transpose()
        csv_filename_cap=f'{folder_out_act}wind_capacity.csv'
        df_capacity.to_csv(csv_filename_cap, index=False, encoding="utf-8-sig")
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')


# class to perform simulations for subgrids defined by the clustering module
class WindCluster():
    def __init__(self):
        print("[+] Begin")
        # read turbine information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and wind turbine data')
            self.wt_data = pd.read_csv(file_red_cls, index_col=False)
            # --- GeoTIFF Factor A, k, v at 100-meter high ---
            self.gwa_A10 = rasterio.open(gwa_src_A_10,'r')
            self.gwa_A100 = rasterio.open(gwa_src_A_100,'r')
            self.gwa_k10 = rasterio.open(gwa_src_k_10,'r')
            self.gwa_k100 = rasterio.open(gwa_src_k_100,'r')
            #self.gwa_v10 = rasterio.open(gwa_src_v_10,'r')
            #self.gwa_v100 = rasterio.open(gwa_src_v_100,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        # plot tiff info
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gwa_A10: {self.gwa_A10.width} x {self.gwa_A10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A10: {self.gwa_A10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_A100: {self.gwa_A100.width} x {self.gwa_A100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_A100: {self.gwa_A100.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k10: {self.gwa_k10.width} x {self.gwa_k10.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k10: {self.gwa_k10.bounds}")
        print(f"[+] Shape of the file (width x height) for gwa_k100: {self.gwa_k100.width} x {self.gwa_k100.height}")
        print(f"[+] Geographic extent (bounds) for gwa_k100: {self.gwa_k100.bounds}")
        print("[+] GeoTIFF files successfully opened")
        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simcl1,yr_gwa1)
        yr_max=max(yr_simcl2,yr_gwa2)
        seltf = [f'{yr}{mn:02}' for yr in range(yr_min,yr_max+1) for mn in range(1,13)]
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
                    print(f"[+] Memory usage raw: {memory_usage_megabytes} MB")
                    print("[+] ERA5 for BC successfully loaded.")
            except FileNotFoundError:
                # print error if exception occurs
                print(f"[!] Error. File not found. Stopping execution.")
                sys.exit(1)               
            except Exception as e:
                # print error if exception occurs
                print(f"[!] Error. Unable to open datasets. Stopping execution.")
                sys.exit(1)

    # check whether coordinates are out of bound or not
    def check_bounds(self, py, px, width, height):
        if not (0 <= px <= width and 0 <= py <= height):
            print(['[!] Error. Coordinates are out of bounds! Check boundaries of the geographical set!'])
            sys.exit(1)

    # method to extract the value from a GeoTIFF file at the specified geographical coordinate
    def get_val(self, lon, lat, param_A, param_k):
        # Convert the geographical coordinates to dataset's pixel coordinates
        Apy, Apx = param_A.index(lon, lat)
        self.check_bounds(Apy, Apx, param_A.width, param_A.height)
        kpy, kpx = param_k.index(lon, lat)
        self.check_bounds(kpy, kpx, param_k.width, param_k.height)
        #vpy, vpx = param_v.index(lon, lat)
        #self.check_bounds(vpy, vpx, param_v.width, param_v.height)
        # window defines the data window we want to read (in this case, a single pixel)
        windowA = rasterio.windows.Window(Apx, Apy, 1, 1)
        windowk = rasterio.windows.Window(kpx, kpy, 1, 1)
        #windowv = rasterio.windows.Window(vpx, vpy, 1, 1)
        # Read the data in the window, assuming the raster has only one band
        bandA = param_A.read(1, window=windowA)
        bandk = param_k.read(1, window=windowk)
        #bandv = param_v.read(1, window=windowv)
        # Extract the value from the array
        A = bandA[0, 0]
        k = bandk[0, 0]
        #v = bandv[0, 0]
        # return values
        return A, k

    # execute simulation for existing clusters
    def run_simulation(self):
        # *************** Wind power curves ***************
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Creating wind power curves with Gaussian filtering: {current_time}")               
        # onshore
        df_wpc_on, interp_wpc_on = wpc_onshore()
        df_wpc_on.to_csv(wind_onwpc_file, encoding="utf-8-sig")
        # offshore
        df_wpc_off, interp_wpc_off = wpc_offshore()
        df_wpc_off.to_csv(wind_offwpc_file, encoding="utf-8-sig")

        # *************** Prepare inputs ***************
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Reading simulation inputs {current_time}")
        clu_top = self.wt_data[self.wt_data['Top'] == True]
        if clu_top.empty:
            print("[!] No subgrids available for the simulation!")
            sys.exit(1)
        col_names = clu_top['id'].to_list()
        latlon=clu_top[['Sub lat','Sub lon']].values

                
        # *************** Selection for GWA timeframe *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating Weibull parameters for subgrids: {current_time}")               
        # the operation will be performed on an yearly basis to reduce memory usage
        cum_sums_10m = pd.Series(0, index=col_names)
        cum_sums_100m = pd.Series(0, index=col_names)
        cum_sums_sq_10m = pd.Series(0, index=col_names)
        cum_sums_sq_100m = pd.Series(0, index=col_names)
        total_hours = 0
        for year in range(yr_gwa1,yr_gwa2+1):
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # section to calculate parameters at 10m
            sel_wnd_u10 = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v10 = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
            d_wnd10_y = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_names)            
            # section to calculate parameters at 100m
            sel_wnd_u100 = era5.u100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v100 = era5.v100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
            d_wnd100_y = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_names)
            # a numerical stability analysis was performed before implementing the approach below to calculate the standard deviation for wind speeds
            total_hours+=d_wnd10_y.shape[0]
            hourly_wnd10=d_wnd10_y.sum()
            cum_sums_10m+=hourly_wnd10
            hourly_sq_wnd10=d_wnd10_y **2
            hourly_sq_wnd10=hourly_sq_wnd10.sum()
            cum_sums_sq_10m+=hourly_sq_wnd10
            # same for 100m
            hourly_wnd100=d_wnd100_y.sum()      # after sum it becomes a Series (instead of DataFrame)
            cum_sums_100m+=hourly_wnd100
            hourly_sq_wnd100=d_wnd100_y **2
            hourly_sq_wnd100=hourly_sq_wnd100.sum()
            cum_sums_sq_100m+=hourly_sq_wnd100
        # calculations for mean and stdev were confirmed by exogenous calculations (all set)
        wnd10_m_mean = cum_sums_10m / total_hours
        wnd10_m_std = np.sqrt((cum_sums_sq_10m-(cum_sums_10m**2 / total_hours)) / total_hours)
        # same for 100m
        wnd100_m_mean = cum_sums_100m / total_hours
        wnd100_m_std = np.sqrt((cum_sums_sq_100m-(cum_sums_100m**2 / total_hours)) / total_hours)
        # calculations for mean and stdev all confirmed with exogenous analysis (all set)
        k10m = (wnd10_m_std / wnd10_m_mean) ** (-1.086)
        A10m = wnd10_m_mean / gamma(1+1/k10m)
        k100m = (wnd100_m_std / wnd100_m_mean) ** (-1.086)
        A100m = wnd100_m_mean / gamma(1+1/k100m)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculation completed at {current_time}")
        # select Weibull parameters using GWA
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating GWA parameters for subgrids at {current_time}")
        A10_arr = []
        A100_arr = []
        k10_arr = []
        k100_arr = []
        # select the coordinates for the current cluster 
        for _, row in clu_top.iterrows():
            # select data for the current project
            lat, lon = row['Sub lat'], row['Sub lon']
            # select Weibull parameters (A and k) from GWA
            Au10_val, ku10_val = self.get_val(lon,lat,self.gwa_A10,self.gwa_k10)
            Au100_val, ku100_val = self.get_val(lon,lat,self.gwa_A100,self.gwa_k100)
            A10_arr.append(Au10_val)
            A100_arr.append(Au100_val)
            k10_arr.append(ku10_val)
            k100_arr.append(ku100_val)
        # create a dataframe that includes the name of the assets and the new values ('u' stands for micro)
        A10u = pd.Series(A10_arr, index=col_names)
        A100u = pd.Series(A100_arr, index=col_names)
        k10u = pd.Series(k10_arr, index=col_names)
        k100u = pd.Series(k100_arr, index=col_names)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculation completed at {current_time}")
        # print status of the dataframe creation
        memuse = psutil.Process(os.getpid())
        print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
        # delete temporary variables
        del sel_wnd_u10, sel_wnd_v10, sel_wnd_10
        del sel_wnd_u100, sel_wnd_v100, sel_wnd_100
        gc.collect()

        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for clusters at {current_time}")
        yrs_count = yr_simcl2 - yr_simcl1 + 1
        # onshore wind
        clu_on_res_yr = []
        clu_on_result = clu_top[clu_top['location']==0]
        clu_on_names = clu_on_result['Label'].to_list()
        clu_on_res_assess = pd.Series(0,index=clu_on_names)    
        # offshore wind
        clu_of_res_yr = []
        clu_of_result = clu_top[clu_top['location']==1]
        clu_of_names = clu_of_result['Label'].to_list()
        clu_of_res_assess = pd.Series(0,index=clu_of_names)
        for year in range(yr_simcl1,yr_simcl2+1):
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # wind speeds 10m
            sel_wnd_u10 = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v10 = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_10 = np.sqrt(sel_wnd_u10**2 + sel_wnd_v10**2)
            wnd10m = pd.DataFrame(sel_wnd_10.values, index=sel_wnd_10.time.values, columns=col_names)
            # wind speeds 100m
            sel_wnd_u100 = era5.u100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v100 = era5.v100.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_100 = np.sqrt(sel_wnd_u100**2 + sel_wnd_v100**2)
            wnd100m = pd.DataFrame(sel_wnd_100.values, index=sel_wnd_100.time.values, columns=col_names)
            # downscale datasets (array division)
            wnd_10 = A10u*((wnd10m/A10m)**(k10m/k10u))
            wnd_100 = A100u*((wnd100m/A100m)**(k100m/k100u))
            alpha = np.log(wnd_100 / wnd_10) / np.log(100 / 10)
            # prepare for multithreading (subtasks for onshore and offshore projects to optimize use of processors)
            date_index = wnd_10.index

            # onshore wind
            subtasks_onw = []
            for _, row in clu_on_result.iterrows():
                id = row['id']
                sel_wnd100 = wnd_100[id]
                sel_alpha = alpha[id]
                hheight = 125
                wnd_adj = (sel_wnd100 * (hheight/100)**sel_alpha)
                wnd_adj = np.where(wnd_adj<=25,wnd_adj,0)
                asset_tuple = (wnd_adj, interp_wpc_on)
                subtasks_onw.append(asset_tuple)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] {current_time}: Running simulation for onshore wind in {year}")
            # run multithreading with tasks
            if len(subtasks_onw)>0:
                with Pool(processes=8) as pool:
                    # run cluster_on_cf to each set of tasks (Gaussian filtering)
                    results_on_yr = pool.starmap(calc_cf,subtasks_onw)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] {current_time}: Simulation completed for onshore wind")
            if len(results_on_yr)>0:
                current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[+] {current_time}: Storing simulation for onshore wind")
                df_results_on_yr = pd.DataFrame(results_on_yr).T
                df_results_on_yr.columns = clu_on_names
                df_results_on_yr.index = date_index
                df_results_on_yr.index.name = 'datetime'
                df_clu_on_avrg_yr = df_results_on_yr.groupby(by=df_results_on_yr.columns, axis=1).mean()  # force by using columns, otherwise pandas will not execute
                clu_on_res_yr.append(df_clu_on_avrg_yr)
                # populate dataframe for resource assessment
                res_assess = df_results_on_yr.sum()
                clu_on_res_assess += res_assess

            # offshore wind
            subtasks_ofw = []
            for _, row in clu_of_result.iterrows():
                id = row['id']
                sel_wnd100 = wnd_100[id]
                sel_alpha = alpha[id]
                hheight = 150
                wnd_adj = (sel_wnd100 * (hheight/100)**sel_alpha)
                wnd_adj = np.where(wnd_adj<=25,wnd_adj,0)
                asset_tuple = (wnd_adj, interp_wpc_off)
                subtasks_ofw.append(asset_tuple)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] {current_time}: Running simulation for offshore wind in {year}")            
            if len(subtasks_ofw)>0:
                with Pool(processes=8) as pool:
                    results_of_yr = pool.starmap(calc_cf,subtasks_ofw)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] {current_time}: Simulation completed for offshore wind")            
            if len(results_of_yr)>0:
                current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[+] {current_time}: Storing simulation for offshore wind")
                df_results_of_yr = pd.DataFrame(results_of_yr).T
                df_results_of_yr.columns = clu_of_names
                df_results_of_yr.index = date_index
                df_results_of_yr.index.name = 'datetime'
                df_clu_of_avrg_yr = df_results_of_yr.groupby(by=df_results_of_yr.columns, axis=1).mean()  # force by using columns, otherwise pandas will not execute
                clu_of_res_yr.append(df_clu_of_avrg_yr)
                # populate dataframe for resource assessment
                res_assess = df_results_of_yr.sum()
                clu_of_res_assess += res_assess                   
        # concatenate results and produce resource assessment
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] {current_time}: Storing simulation results")            
        if len(clu_on_res_yr)>0:
            df_cluster_onshore=pd.concat(clu_on_res_yr)            # concat lists
            clu_on_res_assess.index.name = 'Cluster'
            clu_on_res_assess.name = 'CF'
            clu_on_res_assess = clu_on_res_assess / (yrs_count*8760)        
            csv_filename_cf = f'{folder_out_cls}wind_onshore.csv'
            csv_filename_rass = f'{folder_out_cls}wind_onshore_rass.csv'
            df_cluster_onshore.to_csv(csv_filename_cf, encoding="utf-8-sig")
            clu_on_res_assess.to_csv(csv_filename_rass, encoding="utf-8-sig")
        if len(clu_of_res_yr)>0:
            df_cluster_offshore=pd.concat(clu_of_res_yr)            # concat lists
            clu_of_res_assess.index.name = 'Cluster'
            clu_of_res_assess.name = 'CF'
            clu_of_res_assess = clu_of_res_assess / (yrs_count*8760)        
            csv_filename_cf = f'{folder_out_cls}wind_offshore.csv'
            csv_filename_rass = f'{folder_out_cls}wind_offshore_rass.csv'
            df_cluster_offshore.to_csv(csv_filename_cf, encoding="utf-8-sig")
            clu_of_res_assess.to_csv(csv_filename_rass, encoding="utf-8-sig")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')
