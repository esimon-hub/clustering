"""
Solar PV Simulation Class 

Description:
This class includes the methods needed to perform the simulation of clusters and existing solar PV assets.

Usage:
Import this into the pv_simulation_exec.py to create instances and perform simulations.

For support or inquiries, contact emanueel@gmail.com

"""
#   Importing libraries and settings
from common.pv_config import era_src                                                                                # common settings
from common.pv_config import yr_simex1, yr_simex2, bc_act, file_red_act, folder_out_act                             # actual settings
from common.pv_config import yr_simcl1, yr_simcl2, cl_tech, nth_perc, method_sel, file_red_cls, folder_out_cls      # cluster settings
from common.pv_config import gsa_src                                                                                # GSA file
from common.pv_config import yr_gsa1, yr_gsa2, subdim                                                               # GWA years and dimension
from multiprocessing import Pool
from datetime import datetime
from scipy.spatial import cKDTree
from multiprocessing import Pool
from rasterio.windows import Window
from .pv_shared import complement_halfhour
from .pv_shared import interp_halfhour
from .pv_shared import pv_sim_multi
import xarray as xr
import glob
import sys
import os
import pandas as pd
import numpy as np
import rasterio
import gc           # garbage collector to release unreferenced memory back to OS
import psutil       # measure memory usage

# class to perform simulations for existing solar PV assets
class PVAsset():
    def __init__(self):
        print("[+] Begin")
        # read solar PV plant information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and solar PV data')
            self.pv_data = pd.read_csv(file_red_act, index_col=False)
            self.gsa_fl = rasterio.open(gsa_src,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        # plot tiff info
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gsa: {self.gsa_fl.width} x {self.gsa_fl.height}")
        print(f"[+] Geographic extent (bounds) for gsa: {self.gsa_fl.bounds}")
        print("[+] GeoTIFF files successfully opened")

        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simex1,yr_gsa1)
        yr_max=max(yr_simex2,yr_gsa2)
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
                    print("[+] ERA5 successfully loaded.")
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
    def get_val(self, lon, lat, param):
        # Convert the geographical coordinates to dataset's pixel coordinates
        py, px = param.index(lon, lat)
        self.check_bounds(py, px, param.width, param.height)
        # window defines the data window we want to read (in this case, a single pixel)
        window = rasterio.windows.Window(px, py, 1, 1)
        # Read the data in the window, assuming the raster has only one band
        band = param.read(1, window=window)
        value = band[0, 0]
        return value

    # execute simulation for existing assets
    def run_simulation(self):
        # select data for existing assets
        col_names = self.pv_data['ID'].to_list()
        latlon=self.pv_data[['lat','lon']].values

        # *************** Selection for GSA timeframe *************** 
        if bc_act:
            # calculate BC factors using selected ERA5 timeframe
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Calculating radiation parameters for grid-cells at {current_time}")            
            rtime1=datetime(yr_gsa1,1,1)
            rtime2=datetime(yr_gsa2,12,31,23,0,0)
            r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
            r_date_range_30m_up = r_date_range + pd.Timedelta('30min')                  # generation shifts +30m            
            era5 = self.era5.sel(time=slice(f'{yr_gsa1}-01-01',f'{yr_gsa2}-12-31'))
            # select relevant variables from ERA5
            sel_ghi = era5.ssrd.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            # calculate daily mean for each location (ERA5)
            d_ghim = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_names)
            d_ghim = d_ghim.clip(lower=0) / 3600
            d_ghim = complement_halfhour(d_ghim)            # delete 00:00 & generate 00:30 from 01:00 | add 00:00 & generate 23:30 from 00:00
            d_ghim.index = r_date_range_30m_up
            d_ghim = d_ghim.resample('D').sum()/1000        #kWh/m2
            d_ghim = d_ghim.mean()
            # select GHI estimates using GSA
            GHI_arr = []
            # select the coordinates for the current project 
            for _, row in self.pv_data.iterrows():
                # select data for the current project
                asset_id = row['ID']
                lat, lon = row['lat'], row['lon']
                print(f"[+] Retrieving parameters for asset: {asset_id}")
                # select GHI from GSA
                ghi_val = self.get_val(lon,lat,self.gsa_fl)
                # record values
                GHI_arr.append(ghi_val)
            # create a dataframe that includes the name of the assets and the new values
            d_ghiu = pd.Series(GHI_arr, index=col_names)
            # print status of the dataframe creation
            memuse = psutil.Process(os.getpid())
            print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
            # delete temporary variables
            del sel_ghi
            gc.collect()

        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for existing assets at {current_time}")        
        # execute simulations with or without downscaling
        for year in range(yr_simex1,yr_simex2+1):
            # create date range customized for each year (consider the +- 30m)
            rtime1=datetime(year,1,1)
            rtime2=datetime(year,12,31,23,0,0)
            r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
            r_date_range_30m_up = r_date_range + pd.Timedelta('30min')      # generation shifts +30m
            # select current year before the nearest method
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # calculate GHI
            sel_ghi = era5.ssrd.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            ghim = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_names)
            ghim = ghim.clip(lower=0) / 3600
            ghim = complement_halfhour(ghim)            # delete 00:00 & generate 00:30 from 01:00 | add 00:00 & generate 23:30 from 00:00
            ghim.index = r_date_range_30m_up
            # calculate absolute wind speeds (10m)
            sel_wnd_u = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd = np.sqrt(sel_wnd_u**2 + sel_wnd_v**2)
            sel_wnd = sel_wnd*(np.log(1/0.25)/np.log(10/0.25))  # adjust wind speed to 1 meter high using 0.25 of roughness
            wndm = pd.DataFrame(sel_wnd.values, index=sel_wnd.time.values, columns=col_names)
            wndm = interp_halfhour(wndm)
            # select temperature
            sel_tmp = (era5.t2m.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest') - 273.15).compute()
            tmpm = pd.DataFrame(sel_tmp.values, index=sel_tmp.time.values, columns=col_names)
            tmpm = interp_halfhour(tmpm)
            # downscale datasets
            if bc_act:
                # multiplicative factor (array division)
                ghi_adj = d_ghiu/d_ghim*ghim
            else:
                ghi_adj = ghim
            # prepare for multithreading
            date_index = ghim.index
            subtasks = []
            # iterate over the dataframe to obtain the time series
            for _, row in self.pv_data.iterrows():
                id = row['ID']
                lat, lon = row['lat'], row['lon']
                trx=row['TRACKER']
                ghi = ghi_adj[id]
                wnd = wndm[id]
                tmp = tmpm[id]
                asset_tuple = (ghi, tmp, wnd, lat, lon, trx, r_date_range_30m_up)
                subtasks.append(asset_tuple)
            # plot status before multithreading
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Running simulation for year {year} with multithreading at {current_time}")
            # run multithreading with tasks
            with Pool(processes=4) as pool:
                # run calc_on_cf to each set of tasks
                results_yr = pool.starmap(pv_sim_multi,subtasks)
            # plot status (completed)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Storing simulation for year {year} at {current_time}")
            # create a new dataframe for each year and store into a csv
            df_results_yr = pd.DataFrame(results_yr).T
            df_results_yr.columns = col_names
            df_results_yr.index = date_index
            df_results_yr.index.name = 'datetime'
            df_results_yr = df_results_yr.round(4)      # reduce size in disk
            csv_filename_cf = f'{folder_out_act}solar_cf_{year}.csv'
            df_results_yr.to_csv(csv_filename_cf, encoding="utf-8-sig")
        # print status
        df_capacity = pd.DataFrame(self.pv_data.set_index('ID')['POT_KW']).transpose()
        csv_filename_cap=f'{folder_out_act}solar_capacity.csv'
        df_capacity.to_csv(csv_filename_cap, index=False, encoding="utf-8-sig")
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')

#
# class to perform simulations for existing solar assets using mesh
#
class PVMeshAsset():
    def __init__(self):
        print("[+] Begin")
        # read solar PV plant information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and solar PV data')
            self.pv_data = pd.read_csv(file_red_act, index_col=False)
            self.gsa_fl = rasterio.open(gsa_src,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gsa: {self.gsa_fl.width} x {self.gsa_fl.height}")
        print(f"[+] Geographic extent (bounds) for gsa: {self.gsa_fl.bounds}")
        print("[+] GeoTIFF files successfully opened")

        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simex1,yr_gsa1)
        yr_max=max(yr_simex2,yr_gsa2)
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
                    print("[+] ERA5 successfully loaded.")
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

    # convert bounds to geotiff pixel coordinates for each GHI
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
    def down_grid_zone(self, up_ghi, gm, ghiu):
        # perform downscaling using the given parameters (m - macro; u - micro)
        nk = len(ghiu)
        ghi_k=0
        for ik in range(nk):
            gk = ghiu[ik]
            ghi_k += gk/gm*up_ghi
        # return the average time series
        ghi_k = ghi_k / nk
        return ghi_k

    # calculate the average values using downscaling
    def down_grid_centroid(self, up_ghi, gm, ghiu):
        # perform downscaling using the given parameters (m - macro; u - micro)
        ghi_k = ghiu/gm*up_ghi
        return ghi_k

    # execute simulation for existing assets
    def run_simulation(self):
        # Grid-level data
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Select nearest ERA5 grid-cells for existing assets at {current_time}")
        era5_lat = []
        era5_lon = []
        for _, row in self.pv_data.iterrows():
            # select the nearest point in ERA5 data for the current turbine location
            era5_loc = self.era5.sel(latitude=row['lat'], longitude=row['lon'], method='nearest')
            n_lat = era5_loc.latitude.values.item()
            n_lon = era5_loc.longitude.values.item()
            era5_lat.append(n_lat)
            era5_lon.append(n_lon)
        self.pv_data['Grid lat'] = era5_lat
        self.pv_data['Grid lon'] = era5_lon
        #  Subgrid-level data
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f'[+] Creating meshes for subgrids at {current_time}')
        submesh = self.mesh_generation(subdim)
        sub_points = submesh[['Sub lat','Sub lon']].values
        asset_points = self.pv_data[['lat','lon']].values
        tree = cKDTree(sub_points)                              # create a KDTree for efficient nearest neighbor search
        distances, indices = tree.query(asset_points)         # query the nearest subgrid for each wind turbine
        self.pv_data['Sub lat'] = submesh.iloc[indices]['Sub lat'].values
        self.pv_data['Sub lon'] = submesh.iloc[indices]['Sub lon'].values
        # Drop duplicates and create specific dataframes (up ERA5, down GWA) to improve the efficiency of the code
        sel_up = self.pv_data[['Grid lat','Grid lon']].drop_duplicates().reset_index(drop=True)
        sel_up['id'] = range(len(sel_up))
        up_map = {(lat,lon): idx for lat, lon, idx in sel_up[['Grid lat','Grid lon','id']].to_numpy()}
        self.pv_data['id up'] = self.pv_data.apply(lambda row: up_map[(row['Grid lat'],row['Grid lon'])], axis=1)
        sel_down = self.pv_data[['Grid lat','Grid lon','Sub lat','Sub lon']].drop_duplicates().reset_index(drop=True)
        sel_down['id'] = range(len(sel_down))
        down_map = {(lat,lon,slat,slon): idx for lat, lon, slat, slon, idx in sel_down[['Grid lat','Grid lon','Sub lat','Sub lon','id']].to_numpy()}
        self.pv_data['id down'] = self.pv_data.apply(lambda row: down_map[(row['Grid lat'],row['Grid lon'],row['Sub lat'],row['Sub lon'])], axis=1)

        # *************** Selection for GSA timeframe *************** 
        col_up = sel_up['id'].to_list()                                 # ERA5 columns
        latlon_up = sel_up[['Grid lat','Grid lon']].values              # ERA5 geographical coordinates
        col_down = sel_down['id'].to_list()                             # Subgrid columns
        # retrieve solar irradiation from GSA to perform simulation with BC
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating radiation parameters for grid-cells at {current_time}")
        rtime1=datetime(yr_gsa1,1,1)
        rtime2=datetime(yr_gsa2,12,31,23,0,0)
        r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
        r_date_range_30m_up = r_date_range + pd.Timedelta('30min')                  # generation shifts +30m            
        era5 = self.era5.sel(time=slice(f'{yr_gsa1}-01-01',f'{yr_gsa2}-12-31'))        
        # select relevant variables from ERA5
        sel_ghi = era5.ssrd.sel(
            latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
            longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
            method='nearest').compute()
        # calculate daily mean for each location (ERA5)
        d_ghim = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_up)
        d_ghim = d_ghim.clip(lower=0) / 3600
        d_ghim = complement_halfhour(d_ghim)            # eliminate the first row (23:30) and add the last row (23:30)
        d_ghim.index = r_date_range_30m_up
        d_ghim = d_ghim.resample('D').sum()/1000        #kWh/m2
        d_ghim = d_ghim.mean()
        # create a dictionary with GWA parameters
        ghi_dict = {}
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
                # select GHI from GSA (returns a list)
                ghi_gsa = self.get_geo_zone(bounds, self.gsa_fl)
                # store outputs in dictionaries
                ghi_dict[sid] = ghi_gsa
            elif method_sel==2:
                # select GHI from GSA (returns a list)
                ghi_gsa = self.get_geo_centroid(slon, slat, self.gsa_fl)
                # store outputs in dictionaries
                ghi_dict[sid] = ghi_gsa
        # print status of the dataframe creation
        memuse = psutil.Process(os.getpid())
        print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
        # delete temporary variables
        del sel_ghi
        gc.collect()


        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for existing assets at {current_time}")        
        # execute simulations with or without downscaling
        for year in range(yr_simex1,yr_simex2+1):
            # create date range customized for each year (consider the +- 30m)
            rtime1=datetime(year,1,1)
            rtime2=datetime(year,12,31,23,0,0)
            r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
            r_date_range_30m_up = r_date_range + pd.Timedelta('30min')      # generation shifts +30m
            # select current year before the nearest method
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            # calculate GHI
            sel_ghi = era5.ssrd.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            ghim = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_up)
            ghim = ghim.clip(lower=0) / 3600
            ghim = complement_halfhour(ghim)    # delete 00:00 & generate 00:30 from 01:00 | add 00:00 & generate 23:30 from 00:00
            ghim.index = r_date_range_30m_up
            # calculate absolute wind speeds (10m)
            sel_wnd_u = era5.u10.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v = era5.v10.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd = np.sqrt(sel_wnd_u**2 + sel_wnd_v**2)
            sel_wnd = sel_wnd*(np.log(1/0.25)/np.log(10/0.25))  # adjust wind speed to 1 meter high using 0.25 of roughness
            wndm = pd.DataFrame(sel_wnd.values, index=sel_wnd.time.values, columns=col_up)
            wndm = interp_halfhour(wndm)
            # select temperature
            sel_tmp = (era5.t2m.sel(
                latitude=xr.DataArray(latlon_up[:, 0], dims="points"),
                longitude=xr.DataArray(latlon_up[:, 1], dims="points"),
                method='nearest') - 273.15).compute()
            tmpm = pd.DataFrame(sel_tmp.values, index=sel_tmp.time.values, columns=col_up)
            tmpm = interp_halfhour(tmpm)
            # create dataframe to store downscale
            ghi_sub = pd.DataFrame(index=ghim.index, columns=col_down)
            wnd_sub = pd.DataFrame(index=wndm.index, columns=col_down)
            tmp_sub = pd.DataFrame(index=tmpm.index, columns=col_down)
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
                # fetch data (ghi, wind, temperature) parameters for the corresponding ERA5 grid cell
                up_ghi = ghim.loc[:,id_up]
                up_wnd = wndm.loc[:,id_up]
                up_tmp = tmpm.loc[:,id_up]
                dghiup = d_ghim[id_up]
                dghik = ghi_dict[id]
                # downscale using either zone or centroid
                if method_sel==1:
                    down_ghi = self.down_grid_zone(up_ghi, dghiup, dghik)
                elif method_sel==2:
                    down_ghi = self.down_grid_centroid(up_ghi, dghiup, dghik)
                # store outputs using dataframes
                ghi_sub[id_down] = down_ghi     # down for GHI
                wnd_sub[id_down] = up_wnd       # up for wind
                tmp_sub[id_down] = up_tmp       # up for temperature
            # Running simulations
            subtasks = []
            for _, row in self.pv_data.iterrows():
                id = row['ID']
                id_down = row['id down']
                lat, lon = row['lat'], row['lon']
                trx=row['TRACKER']
                ghi = ghi_sub[id_down]
                wnd = wnd_sub[id_down]
                tmp = tmp_sub[id_down]
                asset_tuple = (ghi, tmp, wnd, lat, lon, trx, r_date_range_30m_up)
                subtasks.append(asset_tuple)
            # plot status before multithreading
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Running simulation for year {year} with multithreading at {current_time}")
            # run multithreading with tasks
            with Pool(processes=4) as pool:
                # run calc_on_cf to each set of tasks
                res_pool_yr = pool.starmap(pv_sim_multi,subtasks)
            # plot status (completed)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Storing simulation for year {year} at {current_time}")
            # create a new dataframe for each year and store into a csv
            results = pd.DataFrame(res_pool_yr).T
            results.columns = self.pv_data['ID']
            results.index = ghim.index      # CHANGE LOG MODIFICATION
            results.index.name = 'datetime'
            results = results.round(4)      # reduce size in disk
            csv_filename_cf = f'{folder_out_act}solar_cf_{year}.csv'
            results.to_csv(csv_filename_cf, encoding="utf-8-sig")
        # print status
        df_capacity = pd.DataFrame(self.pv_data.set_index('ID')['POT_KW']).transpose()
        csv_filename_cap=f'{folder_out_act}solar_capacity.csv'
        df_capacity.to_csv(csv_filename_cap, index=False, encoding="utf-8-sig")
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')

#
# class to perform simulations for subgrids defined by the clustering module
#
class PVCluster():
    def __init__(self):
        print("[+] Begin cluster generation")
        # read solar PV plant information and geotiff data
        try:
            print('[+] Reading GeoTIFF files and solar PV data')
            self.pv_data = pd.read_csv(file_red_cls, index_col=False)
            self.gsa_fl = rasterio.open(gsa_src,'r')
        except FileNotFoundError as e:
            print(f"[!] Error loading files: {e}")
            sys.exit(1)
        except rasterio.errors.RasterioIOError as e:
            print("Unable to open GeoTIFF datasets.")
            sys.exit(1)
        print(f"[+] Initialization completed")
        print(f"[+] Shape of the file (width x height) for gsa: {self.gsa_fl.width} x {self.gsa_fl.height}")
        print(f"[+] Geographic extent (bounds) for gsa: {self.gsa_fl.bounds}")
        print("[+] GeoTIFF files successfully opened")
        # --- ERA5 ---
        print("[+] Reading ERA5 files")
        # selection of files between begin and end year
        yr_min=min(yr_simcl1,yr_gsa1)
        yr_max=max(yr_simcl2,yr_gsa2)
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
                    print("[+] ERA5 successfully loaded.")
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
            print('[!] Error. Coordinates are out of bounds! Check boundaries of the geographical set!')
            sys.exit(1)

    # method to extract the value from a GeoTIFF file at the specified geographical coordinate
    def get_val(self, lon, lat, param):
        # Convert the geographical coordinates to dataset's pixel coordinates
        py, px = param.index(lon, lat)
        self.check_bounds(py, px, param.width, param.height)
        window = rasterio.windows.Window(px, py, 1, 1)
        band = param.read(1, window=window)
        value = band[0, 0]
        return value
    
    # execute simulation for existing clusters
    def run_simulation(self):
        # *************** Prepare inputs ***************
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Reading simulation inputs {current_time}")
        clu_top = self.pv_data[self.pv_data['Top'] == True]
        if clu_top.empty:
            print("[!] No subgrids available for the simulation!")
            sys.exit(1)        
        col_names = clu_top['id'].to_list()
        clu_names = clu_top['Label'].to_list()
        latlon=clu_top[['Sub lat','Sub lon']].values

        # *************** Selection for GSA timeframe *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating ERA5 parameters for subgrids at {current_time}")
        # the operation will be performed on an yearly basis to reduce memory usage
        cum_sums = pd.Series(0, index=col_names)
        total_days=0
        for year in range(yr_gsa1,yr_gsa2+1):
            rtime1=datetime(year,1,1)
            rtime2=datetime(year,12,31,23,0,0)
            r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
            r_date_range_30m_up = r_date_range + pd.Timedelta('30min')                  # generation shifts +30m            
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            sel_ghi = era5.ssrd.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            # calculate daily mean for each location (ERA5)
            d_ghim_y = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_names)
            d_ghim_y = d_ghim_y.clip(lower=0) / 3600
            d_ghim_y = complement_halfhour(d_ghim_y)            # delete 00:00 & generate 00:30 from 01:00 | add 00:00 & generate 23:30 from 00:00
            d_ghim_y.index = r_date_range_30m_up
            d_ghim_y = d_ghim_y.resample('D').sum()/1000        #kWh/m2
            total_days+=d_ghim_y.shape[0]                       # d_ghim_y has now daily entries
            d_ghim_y = d_ghim_y.sum()
            cum_sums+=d_ghim_y
        d_ghim = cum_sums/total_days
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculation completed at {current_time}")
        # select GHI estimates using GSA
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculating GSA parameters for subgrids at {current_time}")
        GHI_arr = []
        for _, row in clu_top.iterrows():
            # select data for the current project
            lat, lon = row['Sub lat'], row['Sub lon']
            ghi_val = self.get_val(lon,lat,self.gsa_fl)
            GHI_arr.append(ghi_val)
        # create a dataframe that includes the name of the clusters and the new values('d' stands for downscale)
        d_ghiu = pd.Series(GHI_arr, index=col_names)
        # adjust values that are located close to borders (water bodies or ocean)
        d_ghiu = d_ghiu.astype(float)
        mean_nonzero = d_ghiu[d_ghiu > 0].mean()
        d_ghiu.fillna(mean_nonzero, inplace=True)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Calculation completed at {current_time}")
        # print status of the dataframe creation
        memuse = psutil.Process(os.getpid())
        print(f"[+] Memory usage: {memuse.memory_info().rss / 1024 ** 2:.2f} MB")
        # delete temporary variables
        del sel_ghi
        gc.collect()

        # *************** Downscaling and Simulation *************** 
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Starting simulations for clusters at {current_time}")       
        clu_res_yr = []
        clu_res_assess = pd.Series(0, index=clu_names)
        yrs_count = yr_simcl2 - yr_simcl1 + 1
        # execute simulations with or without downscaling
        for year in range(yr_simcl1,yr_simcl2+1):
            # create date range customized for each year (consider the +- 30m)
            rtime1=datetime(year,1,1)
            rtime2=datetime(year,12,31,23,0,0)
            r_date_range = pd.date_range(start=rtime1, end=rtime2, freq='H', tz=None)
            r_date_range_30m_up = r_date_range + pd.Timedelta('30min')      # generation shifts +30m
            era5 = self.era5.sel(time=slice(f'{year}-01-01',f'{year}-12-31'))
            sel_ghi = era5.ssrd.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            ghim = pd.DataFrame(sel_ghi.values, index=sel_ghi.time.values, columns=col_names)
            ghim = ghim.clip(lower=0) / 3600
            ghim = complement_halfhour(ghim)            # delete 00:00 & generate 00:30 from 01:00 | add 00:00 & generate 23:30 from 00:00
            ghim.index = r_date_range_30m_up
            # calculate absolute wind speeds (10m)
            sel_wnd_u = era5.u10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd_v = era5.v10.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest').compute()
            sel_wnd = np.sqrt(sel_wnd_u**2 + sel_wnd_v**2)
            sel_wnd = sel_wnd*(np.log(1/0.25)/np.log(10/0.25))  # adjust wind speed to 1 meter high using 0.25 of roughness
            wndm = pd.DataFrame(sel_wnd.values, index=sel_wnd.time.values, columns=col_names)
            wndm = interp_halfhour(wndm)
            # select temperature
            sel_tmp = (era5.t2m.sel(
                latitude=xr.DataArray(latlon[:, 0], dims="points"),
                longitude=xr.DataArray(latlon[:, 1], dims="points"),
                method='nearest') - 273.15).compute()
            tmpm = pd.DataFrame(sel_tmp.values, index=sel_tmp.time.values, columns=col_names)
            tmpm = interp_halfhour(tmpm)
            # downscale datasets (array division)
            ghi_adj = d_ghiu/d_ghim*ghim
            # prepare for multithreading
            date_index = ghim.index
            subtasks = []
            # iterate over the dataframe to obtain the time series
            for _, row in clu_top.iterrows():
                id = row['id']
                lat, lon = row['Sub lat'], row['Sub lon']
                ghi = ghi_adj[id]
                wnd = wndm[id]
                tmp = tmpm[id]
                asset_tuple = (ghi, tmp, wnd, lat, lon, cl_tech, r_date_range_30m_up)
                subtasks.append(asset_tuple)
            # plot status before multithreading
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Running simulation for year {year} with multithreading at {current_time}")
            # run multithreading with tasks
            with Pool(processes=10) as pool:
                # run calc_on_cf to each set of tasks
                results_yr = pool.starmap(pv_sim_multi,subtasks)
            # plot status (completed)
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[+] Storing simulation for year {year} at {current_time}")
            # create a new dataframe for each year and store into a csv
            df_results_yr = pd.DataFrame(results_yr).T
            df_results_yr.columns = clu_names
            df_results_yr.index = date_index
            df_results_yr.index.name = 'datetime'
            df_cluster_average_yr = df_results_yr.groupby(by=df_results_yr.columns, axis=1).mean()  # force by using columns, otherwise pandas will not execute
            clu_res_yr.append(df_cluster_average_yr)
            # populate dataframe for resource assessment
            res_assess = df_results_yr.sum()
            clu_res_assess += res_assess
        # concatenate results and produce resource assessment
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        df_cluster_average=pd.concat(clu_res_yr)
        clu_res_assess.index.name = 'Cluster'
        clu_res_assess.name = 'CF'
        clu_res_assess = clu_res_assess / (yrs_count*8760)
        # store files
        csv_filename_cf = f'{folder_out_cls}solar_tech_{cl_tech}_cf.csv'
        csv_filename_rass = f'{folder_out_cls}solar_tech_{cl_tech}_rass.csv'
        df_cluster_average.to_csv(csv_filename_cf, encoding="utf-8-sig")
        clu_res_assess.to_csv(csv_filename_rass, encoding="utf-8-sig")
        print(f"[+] Process completed at {current_time}")
        print('[+] Done.')
