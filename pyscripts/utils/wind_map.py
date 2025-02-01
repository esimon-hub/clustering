"""

Plot wind resource map

Description:
This script performs the following tasks:
1. Read ERA5 raw data
2. Plot resource maps for both solar PV and wind

Usage:
1. Configure parameters for wind_stat_config.py
2. Run the current script (check the instances below)
3. Outputs will are displayed during the script execution and saved in the determined folder

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import geopandas as gpd
import sys
import glob
import os
import time
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from shapely.geometry import mapping
from wind_stat_config import era5_src, map_path_folder
from wind_stat_config import yr_res1, yr_res2
from wind_stat_config import shapeonoff
from matplotlib.ticker import FormatStrFormatter

# begin
print("[+] Begin")
start_time = time.time()

# --- SHP ---
onshore_data = gpd.read_file(shapeonoff)

# --- ERA5 ---
print("[+] Reading ERA5 files")
# selection of files between begin and end year
seltf = [f'{yr}{mn:02}' for yr in range(yr_res1,yr_res2+1) for mn in range(1,13)]
file_path = []
for pattern in seltf:
    file_path.append(glob.glob(era5_src + '*_' + str(pattern) + '*.nc')[0])

# check whether the content was loaded or not
if not file_path:
    print("[!] Error. The list of file paths is empty. Stopping execution.")
    sys.exit(1)
else:
    try:
        with xr.open_mfdataset(file_path, combine='nested', concat_dim='time', chunks={'lon': 200, 'lat':200}) as era5:
            num_el=str(era5.ssrd.size)                                  # grid points (lat.lon) x 8760
            memory_usage_bytes = era5.ssrd.nbytes
            memory_usage_megabytes = memory_usage_bytes / (1024**2)     # Convert bytes to megabytes
            num_lat = era5.dims['latitude']
            num_lon = era5.dims['longitude']                    
            print(f"[+] Total number of elements (grid points x hours year): {num_el}")
            print(f"[+] Total number of grid cells (original): {num_lat*num_lon}")
            print(f"[+] Memory usage raw: {memory_usage_bytes}")
            print(f"[+] Memory usage raw: {memory_usage_megabytes} MB")
            print("[+] ERA5 WGS84")
            # process ssrd
            ssrd_data=era5.ssrd.sel()
            ssrd_hourly=ssrd_data*(1/3600000)
            ssrd_daily=ssrd_hourly.resample(time='1D').sum()
            ds_weighted = ssrd_daily.groupby("time.season").mean()
            ds_weighted = ds_weighted.rio.write_crs("epsg:4326")
            clipped = ds_weighted.rio.clip(onshore_data.geometry.apply(mapping), onshore_data.crs, drop=False)
            print("[+] ERA5 for BC successfully loaded.")
            # process wind speeds
            u100_data=era5.u100.sel()
            v100_data=era5.v100.sel()
            w_abs = ((u100_data**2+v100_data**2)**0.5).compute()
            ds_weighted = w_abs.groupby("time.season").mean()
            ds_weighted = ds_weighted.rio.write_crs("epsg:4326")
            clipped = ds_weighted.rio.clip(onshore_data.geometry.apply(mapping), onshore_data.crs, drop=False)
            print("[+] ERA5 for BC successfully loaded.")
    except FileNotFoundError:
        # print error if exception occurs
        print(f"[!] Error. File not found. Stopping execution.")
        sys.exit(1)               
    except Exception as e:
        # print error if exception occurs
        print(f"[!] Error. Unable to open datasets. Stopping execution.")
        sys.exit(1)

# calculate averages by quarter DJF / JJA / MAM / SON
lats = ds_weighted.coords['latitude']
lons = ds_weighted.coords['longitude']
# define seasons to be used in charts
seasons=['DJF','MAM','JJA','SON']

# execution time
end_time = time.time()
print(f'[+] Execution time: {end_time-start_time}')

# create figures based on preset parameters (2x2)
n_rows, n_cols = 2, 2
vmin, vmax = 2, 10
levels = np.linspace(vmin, vmax, num=100)
ticks = np.arange(vmin, vmax + 0.5, 0.5)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(7,8.3))
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.11, top=0.99, wspace=0.2, hspace=0.03)
cmap_sel=plt.get_cmap('Spectral_r')
cmap_sel.set_bad('lightgrey')
norm = Normalize(vmin=vmin, vmax=vmax)
mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap_sel)
fig.suptitle('Wind speeds', fontsize=12, fontweight='bold')
# populate subplots
for i in range(n_rows):
    for j in range(n_cols):
        ax = axes[i,j]
        iseason = seasons[i*n_cols+j]
        ax.set_title(f'Season = {iseason}')
        map = Basemap(projection='cyl',llcrnrlat=-33.6,urcrnrlat=5.6, llcrnrlon=-74.2,urcrnrlon=-32.5,resolution='h', ax=ax)
        map.drawcountries(color='black', linewidth=0.2)
        map.drawcoastlines(color='black', linewidth=0.2)
        map.fillcontinents(color='#ffe2ab', lake_color='#b0e0e6')  # Filling lakes with the same color as the ocean for consistency
        map.drawmapboundary(fill_color='#b0e0e6')  # Filling the ocean        
        map.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0])        # draw parallels and meridians.
        map.drawmeridians(np.arange(-180.,180.,15.),labels=[0,0,0,1])
        llons, llats = np.meshgrid(lons, lats)
        x,y = map(llons,llats)
        ws=clipped.sel(season=iseason).values
        avrg = np.nanmean(ws)
        std_dev = np.nanstd(ws)
        # add text below each figure
        text_str = f'\u03BC: {avrg: .2f} m/s\n\u03C3: ±{std_dev: .2f} m/s'
        ax.text(0.5,-0.07, text_str, transform=ax.transAxes, ha='center', va='top', fontsize=10)
        temp_plot = map.contourf(x, y, ws, levels=levels, vmin=vmin, vmax=vmax, cmap=cmap_sel, shading='interp', extend='both')
# create a single colorbar
cbar_ax = fig.add_axes([0.1, 0.06, 0.8, 0.02])  # Adjust these values to fit the colorbar in your layout
cb = ColorbarBase(cbar_ax, cmap=cmap_sel, norm=norm, orientation='horizontal', extend='both', ticks=ticks)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
cb.set_label('Long-term average of wind speeds (m/s)', fontsize=10)
cb.ax.tick_params(labelsize=10)
figure_name = f"wind_map_matrix.png"
save_path = os.path.join(map_path_folder, figure_name)
plt.savefig(save_path, dpi=1000, format='jpg')

# create figures based on preset parameters (1x4)
n_rows, n_cols = 1, 4
vmin, vmax = 2, 10
levels = np.linspace(vmin, vmax, num=100)
ticks = np.arange(vmin, vmax + 0.5, 1)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(7,2.9))
plt.rcParams['font.size'] = 9
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.3, top=0.95, wspace=0.25)
cmap_sel=plt.get_cmap('Spectral_r')
norm = Normalize(vmin=vmin, vmax=vmax)
mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap_sel)
fig.suptitle('Wind speeds', fontsize=10, fontweight='bold')
# populate subplots
for i in range(n_cols):
    ax = axes[i]
    iseason = seasons[i]
    ax.set_title(f'Season = {iseason}', fontsize=9)
    map = Basemap(projection='cyl',llcrnrlat=-33.6,urcrnrlat=5.6, llcrnrlon=-74.2,urcrnrlon=-32.5,resolution='h', ax=ax)
    map.drawcountries(color='black', linewidth=0.2)
    map.drawcoastlines(color='black', linewidth=0.2)
    map.fillcontinents(color='#ffe2ab', lake_color='#b0e0e6')   # Filling lakes with the same color as the ocean for consistency
    map.drawmapboundary(fill_color='#b0e0e6')                   # Filling the ocean
    map.drawparallels(np.arange(-90.,120.,15.),labels=[1,0,0,0], fontsize=9)        # draw parallels and meridians.
    map.drawmeridians(np.arange(-180.,180.,15.),labels=[0,0,0,1], fontsize=9)
    llons, llats = np.meshgrid(lons, lats)
    x,y = map(llons,llats)
    ws=clipped.sel(season=iseason).values
    avrg = np.nanmean(ws)
    std_dev = np.nanstd(ws)
    # add text below each figure
    text_str = f'\u03BC: {avrg: .2f} m/s\n\u03C3: ±{std_dev: .2f} m/s'
    ax.text(0.5,-0.12, text_str, transform=ax.transAxes, ha='center', va='top', fontsize=9)
    temp_plot = map.contourf(x, y, ws, levels=levels, vmin=vmin, vmax=vmax, cmap=cmap_sel, shading='interp', extend='both')
# create a single colorbar
cbar_ax = fig.add_axes([0.15, 0.17, 0.7, 0.03])  # Adjust these values to fit the colorbar in your layout
cb = ColorbarBase(cbar_ax, cmap=cmap_sel, norm=norm, orientation='horizontal', extend='both', ticks=ticks)
cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
cb.set_label('Long-term average of wind speeds (m/s)', fontsize=9)
cb.ax.tick_params(labelsize=9)  # Adjust tick label size for clarity
# save plot into a file
figure_name_png = f"wind_map_uni.png"  # Construct the filename
figure_name_pdf = f"wind_map_uni.pdf"
save_path_png = os.path.join(map_path_folder, figure_name_png)  # Combine the directory with the filename
save_path_pdf = os.path.join(map_path_folder, figure_name_pdf)
plt.savefig(save_path_png, dpi=1000, format='jpg')   # Save the figure
plt.savefig(save_path_pdf, format='pdf')

# close nc
era5.close()

