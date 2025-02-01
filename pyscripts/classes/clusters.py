"""
Clustering and Cell Selection Class 

Description:
This class includes the methods needed to perform the clustering of solar and wind resources.

Usage:
Import this into the clustering_exe.py to create instances and create clusters.

For support or inquiries, contact emanueel@gmail.com

"""
import pandas as pd
import numpy as np
import os
import rasterio.windows
import shapefile
import rasterio
from rasterio.windows import Window
import sys
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
import math
from shapely.geometry import Point, box
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans

#   Importing libraries and settings
from common.main_config import climatevar, max_clusters, subdim, exop                                                       # general parameters
from common.main_config import transon, genon, labels, k_single                                                             # flags and labels
from common.main_config import solar_geocoord, solar_geoindex, solar_space, solar_avrg, solar_db_out, solar_db_out_label    # solar files
from common.main_config import wind_geocoord, wind_geoindex, wind_space, wind_avrg, wind_db_out, wind_db_out_label          # wind files
from common.main_config import solar_geoindex_onshore, solar_geoindex_offshore                                              # dataframe with indexes
from common.main_config import wind_geoindex_onshore, wind_geoindex_offshore                                                # dataframe with indexes
from common.main_config import solar_k_subimg                                                                               # solar selected k (subplots)
from common.main_config import wind_k_subimg                                                                                # wind selected k (subplots)
from common.main_config import shp_regional_arg1, shp_regional_arg2                                                         # regional shapefiles
from common.main_config import shp_ex_arg1, shp_ex_arg2, shp_plan_arg1, shp_plan_arg2                                       # transmission shapefiles
from common.main_config import shp_solar_arg1, shp_solar_arg2                                                               # solar shapefiles
from common.main_config import shp_wind_arg1, shp_wind_arg2                                                                 # wind shapefiles
from common.main_config import solar_cluster_fig1, solar_cluster_fig2, solar_cluster_fig3a, solar_cluster_fig3b             # solar figures
from common.main_config import solar_cluster_pdf1, solar_cluster_pdf2, solar_cluster_pdf3b                                  # pdf
from common.main_config import wind_cluster_fig1, wind_cluster_fig2, wind_cluster_fig3a, wind_cluster_fig3b                 # wind figures
from common.main_config import wind_cluster_pdf1, wind_cluster_pdf2, wind_cluster_pdf3b                                     # pdf
from common.main_config import solar_cluster_roc                                                                            # solar roc (csv)
from common.main_config import wind_cluster_roc                                                                             # wind roc (csv)
from common.main_config import solar_cluster_fig4a, solar_cluster_fig4b, solar_cluster_fig4_label                           # final cluster with labels
from common.main_config import solar_cluster_pdf4b                                                                          # pdf
from common.main_config import wind_cluster_fig4a, wind_cluster_fig4b, wind_cluster_fig4_label                              # final cluster with labels
from common.main_config import wind_cluster_pdf4b                                                                           # pdf
from common.main_config import solar_singlek, solar_subgrid, solar_subgrid_sim                                              # dimension of the subgrid 
from common.main_config import wind_singlek, wind_subgrid, wind_subgrid_sim                                                 # dimension of the subgrid
from common.main_config import shp_wdpa_poly, shp_wdpa_pnt, gpw_tiff, globcover, gebco, srtm                                # subgrid constraints
from common.main_config import solar_assets_exop, wind_assets_exop                                                          # subgrid constraints for existing resources
from common.main_config import gsa_src, gwa_src

# class to perform simulations for existing solar PV assets
class Clusters():
    def __init__(self):
        # read coords and space structure
        if climatevar==1:
            print("[+] Reading solar coordinates and indices")
            self.geocoords = pd.read_csv(solar_geocoord)
            keep_index = pd.read_csv(solar_geoindex)
            self.keep_index = keep_index['keep_index'].values
            print("[+] Reading solar space vector")
            self.space = pd.read_csv(solar_space)
            self.average_res = pd.read_csv(solar_avrg)
            # geographical coordinates for onshore and offshore locations
            onshore_index = pd.read_csv(solar_geoindex_onshore)
            self.onshore_index = onshore_index['onshore_index'].values
            offshore_index = pd.read_csv(solar_geoindex_offshore)
            self.offshore_index = offshore_index['offshore_index'].values
            # shapefile for generation
            self.shapegen_arg1 = shp_solar_arg1
            self.shapegen_arg2 = shp_solar_arg2
            # figures
            self.fig1 = solar_cluster_fig1
            self.pdf1 = solar_cluster_pdf1
            self.fig2 = solar_cluster_fig2
            self.pdf2 = solar_cluster_pdf2
            self.fig3a = solar_cluster_fig3a
            self.fig3b = solar_cluster_fig3b
            self.pdf3b = solar_cluster_pdf3b
            self.fig4a = solar_cluster_fig4a
            self.fig4b = solar_cluster_fig4b
            self.pdf4b = solar_cluster_pdf4b
            self.fig4_label = solar_cluster_fig4_label
            self.roc = solar_cluster_roc
            self.sel_clusters = solar_k_subimg
            # singlek (cluster csv) and subgrid (csv to perform cluster simulations)
            self.clusterk = solar_singlek
            self.cluster_subgrid = solar_subgrid
            self.cluster_subgrid_sim = solar_subgrid_sim
            self.cluster_act = solar_assets_exop
            # db out (locations to perform analysis)
            self.db_out_path = solar_db_out
            self.db_out_path_label = solar_db_out_label
            # atlas
            self.atlas = gsa_src
        elif climatevar==2:
            print("[+] Reading wind coordinates and indices")
            self.geocoords = pd.read_csv(wind_geocoord)
            keep_index = pd.read_csv(wind_geoindex)
            self.keep_index = keep_index['keep_index'].values
            print("[+] Reading wind space vector")
            self.space = pd.read_csv(wind_space)
            self.average_res = pd.read_csv(wind_avrg)
            # geographical coordinates for onshore and offshore locations
            onshore_index = pd.read_csv(wind_geoindex_onshore)
            self.onshore_index = onshore_index['onshore_index'].values
            offshore_index = pd.read_csv(wind_geoindex_offshore)
            self.offshore_index = offshore_index['offshore_index'].values
            # shapefile for generation
            self.shapegen_arg1 = shp_wind_arg1
            self.shapegen_arg2 = shp_wind_arg2
            # figures
            self.fig1 = wind_cluster_fig1
            self.pdf1 = wind_cluster_pdf1
            self.fig2 = wind_cluster_fig2
            self.pdf2 = wind_cluster_pdf2
            self.fig3a = wind_cluster_fig3a
            self.fig3b = wind_cluster_fig3b
            self.pdf3b = wind_cluster_pdf3b
            self.fig4a = wind_cluster_fig4a
            self.fig4b = wind_cluster_fig4b
            self.pdf4b = wind_cluster_pdf4b
            self.fig4_label = wind_cluster_fig4_label
            self.roc = wind_cluster_roc
            self.sel_clusters = wind_k_subimg
            # singlek (cluster csv) and subgrid (csv to perform cluster simulations)
            self.clusterk = wind_singlek
            self.cluster_subgrid = wind_subgrid
            self.cluster_subgrid_sim = wind_subgrid_sim
            self.cluster_act = wind_assets_exop
            # db out (locations to perform analysis)
            self.db_out_path = wind_db_out
            self.db_out_path_label = wind_db_out_label
            # atlas
            self.atlas = gwa_src
        else:
            print("[!] Climate data not found")
        # legend for k_plots and wcss
        self.leg_index = ['a','b','c','d','e','f','g','h']
        self.era5_size = 0.25


    # calculate the within-cluster sum of squares (Elbow)
    def calc_wcss(self):

        print('hello')
        quit()
        # Climate variable selection
        if climatevar==1:
            title = 'Solar radiation'
            cutoff=15
            xoff1=-11
            yoff1=-12
            xoff2=-11
            yoff2=15
        elif climatevar==2:
            title = 'Wind speeds'
            cutoff=16
            xoff1=-11
            yoff1=-13
            xoff2=-11
            yoff2=15
        else:
            print("[!] Climate data not found")
        
        # calculate the sum of the squared distance between each point and the centroid in a cluster
        print("[+] Calculate the WCSS (within cluster sum of squares)")
        wcss = []
        var = []
        # WCSS
        for n_clusters in range(1,max_clusters+1):
            kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(self.space.T)
            wcss.append(kmeans.inertia_)
            if n_clusters==1:
                initial_wcss = kmeans.inertia_
            var.append((kmeans.inertia_ / initial_wcss) * 100)
        # Calculate rate of change of explained variance
        roc = [var[n - 1] - var[n] for n in range(1, len(var))]
        roc.insert(0, None)
        
        # map details
        plt.rcParams['font.size'] = 9
        fig, axs = plt.subplots(1, 2, figsize=(7, 3))  # 7/3 Adjusted the overall figsize for better display
        fig.suptitle(title, fontsize=10, fontweight='bold')
        plt.subplots_adjust(left=0.1, right=0.97, bottom=0.15, top=None, wspace=0.3, hspace=None)
        
        # WCSS
        axs[0].plot(range(1, len(wcss)+1), wcss, 'ko', markersize=3, markerfacecolor='black', markeredgewidth=0)
        axs[0].set_xlabel('Number of clusters (k)', fontsize=9)
        axs[0].set_ylabel('WCSS', fontsize=9)
        axs[0].grid(True, linestyle='--', linewidth=0.5, color='grey')
        
        # Rate of change
        axs[1].plot(range(1, len(var)+1), roc, 'ko', markersize=3, markerfacecolor='black', markeredgewidth=0)
        axs[1].axhline(y=1, color='blue', linestyle='--', linewidth=1)   
        axs[1].annotate('1%', xy=(1, 1), xytext=(1, 1.5),
                           textcoords="data", ha='left', fontsize=8, color='blue')
        axs[1].set_xlabel('Number of clusters (k)', fontsize=9)
        axs[1].set_ylabel('Rate of change (%)', fontsize=9)
        axs[1].grid(True, linestyle='--', linewidth=0.5, color='grey')

        # Set the y-axis to scientific notation
        max_wcss = max(wcss)
        max_roc = max(roc[1:])
        axs[0].set_ylim(0, max_wcss * 1.1)      # Add 10% buffer to the top  
        axs[1].set_ylim(0, max_roc * 1.1)  

        # highlight selected clusters WCSS
        for i, cluster in enumerate(self.sel_clusters):
            if cluster <= len(wcss):
                axs[0].scatter(cluster, wcss[cluster-1], color='yellow', s=3, label=self.leg_index[i], zorder=3)
                # Annotating selected clusters
                axs[0].annotate(self.leg_index[i], (cluster, wcss[cluster-1]), textcoords="offset points", xytext=(3, 5), ha='right', fontsize=8, color='black')

        # highlight selected clusters RoC
        self.sel_clusters = self.sel_clusters[1:]
        self.leg_index = self.leg_index[1:]
        for i, cluster in enumerate(self.sel_clusters):
            if cluster <= len(wcss):
                axs[1].scatter(cluster, roc[cluster-1], color='yellow', s=3, label=self.leg_index[i], zorder=3)
                # Annotating selected clusters
                axs[1].annotate(self.leg_index[i], (cluster, roc[cluster-1]), textcoords="offset points", xytext=(3, 5.5), ha='right', fontsize=8, color='black')                

        # cluster selection annotations (WCSS)
        wcss_val = wcss[cutoff-1]
        axs[0].scatter(cutoff, wcss_val, color='red', s=3, zorder=3)
        axs[0].annotate(f'k={cutoff}', (cutoff, wcss_val), textcoords="offset points", xytext=(xoff1, yoff1),
                        ha='left', fontsize=8, color='red')
        
        # cluster selection annotations (variance)
        roc_val = roc[cutoff-1]
        axs[1].scatter(cutoff, roc_val, color='red', s=3, zorder=3)
        axs[1].annotate(f'k={cutoff}', (cutoff, roc_val), textcoords="offset points", xytext=(xoff2, yoff2),
                        ha='left', fontsize=8, color='red')
        
        # store figure and save dataframe
        plt.grid(True)
        plt.savefig(self.fig1, dpi=1000, format='jpg')   # Save the figure
        plt.savefig(self.pdf1, format='pdf')
        df_roc = pd.DataFrame({
            'Cluster': range(1, max_clusters+1),
            'WCSS': wcss,
            'Rate of change': roc
        })
        df_roc.to_csv(self.roc,index=False)
        print("[+] WCSS done.")


    # plot multiple cluster configurations (varying number of clusters)
    def ksubplots(self):
        # Climate variable selection
        if climatevar==1:
            title = 'Solar radiation'
        elif climatevar==2:
            title = 'Wind speeds'
        else:
            print("[!] Climate data not found")

        # create subplots for clustering results
        n_rows, n_cols = 2, 4
        leg_index = ['a','b','c','d','e','f','g','h']
        plt.rcParams['font.size'] = 9
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(7, 3.6))  # Adjusted the overall figsize for better display
        fig.suptitle(title, fontsize=10, fontweight='bold')
        plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.93, wspace=0.04, hspace=0.04)
        # prepare to rebuild dimensions
        geo_dim=self.geocoords.shape[0]
        keep_index = self.keep_index
        for i, ax in enumerate(axes.flat):
            # perform K-Means clustering
            kmeans = KMeans(n_clusters=self.sel_clusters[i], random_state=42).fit(self.space.T)
            cluster_vals=kmeans.labels_
            # rebuild the dataframe with all elements (lat / lon)
            rebuilt_geo = np.full(geo_dim, np.nan)
            rebuilt_geo[keep_index] = cluster_vals
            # create a basemap for each subplot
            bmap = Basemap(ax=ax, projection='merc',llcrnrlat=-35,urcrnrlat=7, llcrnrlon=-77,urcrnrlon=-31,resolution='h')
            # add details to basemap
            bmap.readshapefile(shp_regional_arg1, shp_regional_arg2, drawbounds=True)
            bmap.shadedrelief() # can also use map.etopo()
            bmap.drawmapscale(-37, -31.5, -37, -31.5, 1000, barstyle='simple', units='km', fontsize=8, labelstyle='simple', fillcolor1='b', fillcolor2='white')
            # add scatter
            val_lat = self.geocoords['lat'].values
            val_lon = self.geocoords['lon'].values
            val_cls = rebuilt_geo
            bmap.scatter(val_lon, val_lat, c=val_cls, cmap='tab20c', latlon=True, marker='s', edgecolors='none', s=1)  # Square markers
            # add details about the cluster being displayed
            ax.text(0.95, 0.95, f'k={self.sel_clusters[i]}', color='yellow', transform=ax.transAxes, fontsize=8, ha='right', va='top', fontweight='bold')
            circle = Circle((0.08, 0.92), 0.05, color='yellow', transform=ax.transAxes, clip_on=False)
            ax.add_patch(circle)
            ax.text(0.08, 0.92, self.leg_index[i], color='black', fontsize=8, ha='center', va='center', transform=ax.transAxes)
        # plot the final version of the chart
        plt.savefig(self.fig2, dpi=1000, format='jpg')   # Save the figure
        plt.savefig(self.pdf2, format='pdf')
        print("[+] Clusters k done.")

    # set the number of clusters and plot the chart
    def singlek_pre(self):
        # prepare to rebuild dimensions
        geo_dim=self.geocoords.shape[0]
        keep_index = self.keep_index
        onshore_index = self.onshore_index
        offshore_index = self.offshore_index
        
        # create plot for a single k
        plt.figure(figsize=(3.5,3.4))   # single column
        kmeans = KMeans(n_clusters=k_single, random_state=42).fit(self.space.T)
        cluster_vals=kmeans.labels_
        
        # rebuild the spatial structure with onshore and offshore locations
        rebuilt_geo = np.full(geo_dim, np.nan)
        rebuilt_geo[onshore_index] = 0          # 0 for onshore cells and 1 for offshore cells
        rebuilt_geo[offshore_index] = 1
        # rebuild the cluster relationship with lat and lon
        rebuilt_clu = np.full(geo_dim, np.nan)
        rebuilt_clu[keep_index] = cluster_vals
        
        # create a basemap for each subplot
        bmap = Basemap(projection='merc',llcrnrlat=-35,urcrnrlat=7, llcrnrlon=-77,urcrnrlon=-31,resolution='h')
        val_lat = self.geocoords['lat'].values
        val_lon = self.geocoords['lon'].values
        val_cls = rebuilt_clu
        bmap.scatter(val_lon, val_lat, c=val_cls, label=val_cls, cmap='tab20c', latlon=True, marker='s', edgecolors='none', s=1.3)    # single column
        bmap.readshapefile(shp_regional_arg1, shp_regional_arg2, drawbounds=True)
        bmap.shadedrelief() # can also use map.etopo()
        bmap.drawmapscale(-36, -32, -36, -32, 700, barstyle='simple', units='km', fontsize=8, labelstyle='simple', fillcolor1='b', fillcolor2='white')

        # add labels for clusters
        if labels:
            # draw meridians
            bmap.drawmeridians(np.arange(round(val_lon.min()), round(val_lon.max()), 2), labels=[0,0,0,1], fontsize=2)  # Adjust interval as necessary
            bmap.drawparallels(np.arange(round(val_lat.min()), round(val_lat.max()), 2), labels=[1,0,0,0], fontsize=2)  # Adjust interval as necessary
            # add labels for grids
            data_labels = list(zip(val_lon,val_lat,val_cls))    # create a list of tuples
            flip=0
            for lon,lat,cluster in data_labels:
                flip+=1
                if not math.isnan(cluster):
                    if flip>5:
                        flip=0
                        x,y = bmap(lon,lat)
                        plt.text(x,y,int(cluster),fontsize=4, color='red')

        print("[+] Saving file with clusters and locations")
        clusterk = pd.DataFrame({
            'lat': val_lat,
            'lon': val_lon,
            'location': rebuilt_geo,
            'cluster': val_cls
        })
        clusterk.to_csv(self.clusterk, index=False)
        # plot the final version of the chart
        print("[+] Saving figure with clusters")
        #plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.tight_layout
        plt.savefig(self.fig3a, dpi=1000, format='jpg')   # Save the figure
        print("[+] Selected cluster done.")

    # plot map with existing assets (appendix)
    def opassets(self):
        # create plot for a single k
        plt.figure(figsize=(3.5,3.4))   # single column

        # create a basemap for each subplot
        bmap = Basemap(projection='merc',llcrnrlat=-35,urcrnrlat=7, llcrnrlon=-77,urcrnrlon=-31,resolution='h')
        bmap.readshapefile(shp_regional_arg1, shp_regional_arg2, drawbounds=True)
        bmap.shadedrelief() # can also use map.etopo()
        bmap.drawmapscale(-36, -32, -36, -32, 700, barstyle='simple', units='km', fontsize=8, labelstyle='simple', fillcolor1='b', fillcolor2='white')

        # add chart legend
        if climatevar==1:
            leg_label = 'Solar PV'
        elif climatevar==2:
            leg_label = 'Onshore wind'
        # add generation layer
        if genon and not transon:
            bmap.readshapefile(self.shapegen_arg1, self.shapegen_arg2)
            for shape in getattr(bmap,self.shapegen_arg2):
                bmap.plot(shape[0],shape[1],marker='.',alpha=1,c='black',linestyle='None',markersize=2,linewidth=0.2,markeredgewidth=0.2)
            leg1 = mlines.Line2D([], [], color='black', marker='.', markersize=4, label=leg_label, linestyle='None', markeredgewidth=0.3)
            plt.legend(handles=[leg1], loc='lower left', fontsize=8, frameon=True, title_fontsize=8)
        # add transmission layer
        if not genon and transon:
            colormap=['blue','red']
            bmap.readshapefile(shp_ex_arg1, shp_ex_arg2, linewidth=0.3, color=colormap[0])
            bmap.readshapefile(shp_plan_arg1, shp_plan_arg2, linewidth=0.3, color=colormap[1], zorder=1)
            legend=['Existing line', 'Planned line']
            leg1 = mlines.Line2D([0], [0], color=colormap[0], label=legend[0])
            leg2 = mlines.Line2D([0], [0], color=colormap[1], label=legend[1])
            plt.legend(handles=[leg1, leg2], loc='lower left', fontsize=8, frameon=True, title_fontsize=8)
        if genon and transon:
            bmap.readshapefile(self.shapegen_arg1, self.shapegen_arg2)
            for shape in getattr(bmap,self.shapegen_arg2):
                bmap.plot(shape[0],shape[1],marker='.',alpha=1,c='black',linestyle='None',markersize=2,linewidth=0.2,markeredgewidth=0.2, zorder=2)
            leg1 = mlines.Line2D([], [], color='black', marker='.', markersize=4, label=leg_label, linestyle='None', markeredgewidth=0.3)
            colormap=['blue','red']
            bmap.readshapefile(shp_ex_arg1, shp_ex_arg2, linewidth=0.3, color=colormap[0], zorder=1)
            bmap.readshapefile(shp_plan_arg1, shp_plan_arg2, linewidth=0.3, color=colormap[1], zorder=1)
            legend=['Existing line', 'Planned line']
            leg2 = mlines.Line2D([0], [0], color=colormap[0], label=legend[0])
            leg3 = mlines.Line2D([0], [0], color=colormap[1], label=legend[1])
            plt.legend(handles=[leg1, leg2, leg3], loc='lower left', fontsize=8, frameon=True, title_fontsize=8)
            
        # plot the final version of the chart
        print("[+] Saving figure with operating assets")
        #plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.tight_layout()
        plt.savefig(self.fig3b, dpi=1000, format='jpg')   # Save the figure
        plt.savefig(self.pdf3b, format='pdf')
        print("[+] Operating assets done.")


    # create a mesh with the subgrids for all locations (used within filtered)
    def mesh_generation(self, subdim):
        # calculate the size of each subgrid based on the subdim
        subsize = self.era5_size / subdim
        # read the geographical coordinates of the SVD analysis (get unique and sorted lat/lon values)
        arr_lat = np.unique(self.geocoords['lat'].values)
        arr_lon = np.unique(self.geocoords['lon'].values)
        # structure to store the mesh information
        mesh_info = []
        # create an extended array to allow for the geographical boundaries
        for lat in arr_lat:
            for lon in arr_lon:
                # generate the subgrid for the current cell
                for i in range(subdim):
                    for j in range(subdim):
                        # calculate the center (centroid) of each subgrid cell
                        sub_center_lat = lat - (i + 0.5) * subsize + self.era5_size/2
                        sub_center_lon = lon + (j + 0.5) * subsize - self.era5_size/2

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

    # create function to calculate overlaps and return a new dataframe
    def exclusion_zone(self, mesh_df):
        # SVD file (lat,lon,location,cluster)
        if os.path.exists(self.clusterk):
            svd_ck = pd.read_csv(self.clusterk)
        else:
            print(f"File not found: {self.clusterk}")
            sys.exit(1)
        # WDPA poly
        if os.path.exists(shp_wdpa_poly):
            shp_poly = gpd.read_file(shp_wdpa_poly)
        else:
            print(f"File not found: {shp_wdpa_poly}")
            sys.exit(1)
        # WDPA points
        if os.path.exists(shp_wdpa_pnt):
            shp_point = gpd.read_file(shp_wdpa_pnt)
        else:
            print(f"File not found: {shp_wdpa_pnt}")
            sys.exit(1)
        # GlobCover data
        if os.path.exists(globcover):
            globcover_pixel = rasterio.open(globcover)
        else:
            print(f"File not found: {globcover}")
            sys.exit(1)
        # GPWV4
        if os.path.exists(gpw_tiff):
            gpw_pixel = rasterio.open(gpw_tiff)
        else:
            print(f"File not found: {gpw_tiff}")
            sys.exit(1)
        # GEBCO
        if os.path.exists(gebco):
            gebco_pixel = rasterio.open(gebco)
        else:
            print(f"File not found: {gebco}")
            sys.exit(1)
        # SRTM-3
        if os.path.exists(srtm):
            srtm_pixel = rasterio.open(srtm)
        else:
            print(f"File not found: {srtm}")
            sys.exit(1)

        # match the corresponding clusters from SVD in mesh (create a new column with the cluster)
        mesh_df = mesh_df.rename(columns={'Grid lat': 'lat', 'Grid lon': 'lon'})
        mesh = pd.merge(mesh_df, svd_ck, on=['lat', 'lon'], how='left')
        mesh = mesh.dropna(subset=['cluster','location'])                       # keep just subgrids that have an associated cluster and location
        mesh = mesh.rename(columns={'lat': 'Grid lat', 'lon': 'Grid lon'})

        # eliminate polygons within WDPA
        gdf_mesh = gpd.GeoDataFrame(mesh, geometry=gpd.points_from_xy(mesh['Sub lon'], mesh['Sub lat']))
        gdf_mesh.set_crs('EPSG:4326', inplace=True)
        gdf_mesh = gpd.sjoin(gdf_mesh, shp_poly, how='left', op='intersects')
        gdf_mesh = gdf_mesh[gdf_mesh['index_right'].isnull()]
        mesh = gdf_mesh[['Grid lat', 'Grid lon', 'Sub lat', 'Sub lon', 'cluster', 'location']]    # now it becomes pandas

        # eliminate points within WDPA
        subsize = 0.25 / subdim
        geo_poly = mesh.apply(lambda row: box(
            row['Sub lon'] - subsize / 2, row['Sub lat'] - subsize / 2,
            row['Sub lon'] + subsize / 2, row['Sub lat'] + subsize / 2), axis=1)
        gdf_mesh = gpd.GeoDataFrame(mesh, geometry=geo_poly)
        gdf_mesh.set_crs('EPSG:4326', inplace=True)
        gdf_mesh = gpd.sjoin(gdf_mesh, shp_point, how='left', op='contains')
        gdf_mesh = gdf_mesh[gdf_mesh['index_right'].isnull()]
        mesh = gdf_mesh[['Grid lat', 'Grid lon', 'Sub lat', 'Sub lon', 'cluster', 'location']]    # now it becomes pandas

        # eliminate unwanted values/classes of landcover (globcover)
        drop_classes = [40, 50, 60, 70, 90, 100, 160, 170, 190, 210, 220]
        drop_indices = []
        for index, row in mesh.iterrows():
            # Only process rows where location is onshore
            if row['location'] == 0:
                # Access geographical coordinates
                sub_lat = row['Sub lat']
                sub_lon = row['Sub lon']
                # Convert the geographical coordinates to dataset's pixel coordinates
                py, px = globcover_pixel.index(sub_lon, sub_lat)
                if not (0 <= px < globcover_pixel.width and 0 <= py < globcover_pixel.height):
                    print(f'[!] Error (GlobCover). Coordinates are out of bounds! Lat: {sub_lat} Lon: {sub_lon}')
                    continue
                # window defines the data window we want to read (in this case, a single pixel)
                window = rasterio.windows.Window(px, py, 1, 1)
                # Read the data in the window, assuming the raster has only one band
                band = globcover_pixel.read(1, window=window)
                value = band[0, 0]
                # Check if the land cover class is one of the unwanted classes
                if value in drop_classes:
                    drop_indices.append(index)
        # drop rows from Dataframe and reset index
        mesh = mesh.drop(drop_indices).reset_index(drop=True)

        # eliminate subgrids above defined threshold from GPWv4
        threshold_pop = 150
        drop_indices = []
        for index, row in mesh.iterrows():
            # Only process rows where location is onshore
            if row['location'] == 0:
                # Access geographical coordinates
                sub_lat = row['Sub lat']
                sub_lon = row['Sub lon']
                # Convert the geographical coordinates to dataset's pixel coordinates
                py, px = gpw_pixel.index(sub_lon, sub_lat)
                if not (0 <= px < gpw_pixel.width and 0 <= py < gpw_pixel.height):
                    print(f'[!] Error (GPWv4). Coordinates are out of bounds! Lat: {sub_lat} Lon: {sub_lon}')
                    continue
                # window defines the data window we want to read (in this case, a single pixel)
                window = rasterio.windows.Window(px, py, 1, 1)
                # Read the data in the window, assuming the raster has only one band
                band = gpw_pixel.read(1, window=window)
                value = band[0, 0]
                # eliminate subgrids with population density higher than threshold
                if value >= threshold_pop:
                    drop_indices.append(index)
        # drop rows from Dataframe and reset index
        mesh = mesh.drop(drop_indices).reset_index(drop=True)

        # eliminate subgrids below water depth of 50m (GEBCO)
        threshold_depth = -50
        drop_indices = []
        for index, row in mesh.iterrows():
            # Only process rows where location is offshore
            if row['location'] == 1:
                # Access geographical coordinates
                sub_lat = row['Sub lat']
                sub_lon = row['Sub lon']
                # Convert the geographical coordinates to dataset's pixel coordinates
                py, px = gebco_pixel.index(sub_lon, sub_lat)
                if not (0 <= px < gebco_pixel.width and 0 <= py < gebco_pixel.height):
                    print(f'[!] Error (GEBCO). Coordinates are out of bounds! Lat: {sub_lat} Lon: {sub_lon}')
                    continue
                # window defines the data window we want to read (in this case, a single pixel)
                window = rasterio.windows.Window(px, py, 1, 1)
                # Read the data in the window, assuming the raster has only one band
                band = gebco_pixel.read(1, window=window)
                value = band[0, 0]
                # eliminate subgrids with water depth below than threshold
                if value < threshold_depth:
                    drop_indices.append(index)
        # drop rows from Dataframe and reset index
        mesh = mesh.drop(drop_indices).reset_index(drop=True)
       
        # eliminate subgrids with steep slope and altitude (SRTM-3)
        drop_indices = []
        max_slope = np.arctan(0.20)     # output in radians
        # iterate through rows to evaluate subgrids
        for index, row in mesh.iterrows():
            # Only process rows where location is onshore
            if row['location'] == 0:
                # Access geographical coordinates
                sub_lat = row['Sub lat']
                sub_lon = row['Sub lon']
                # Calculate corners for each subgrid
                halfsub = self.era5_size / (subdim*2)
                corners = [
                    (sub_lat - halfsub, sub_lon - halfsub),
                    (sub_lat - halfsub, sub_lon + halfsub),
                    (sub_lat + halfsub, sub_lon - halfsub),
                    (sub_lat + halfsub, sub_lon + halfsub)
                ]
                # calculate distance between centroid and corners using Haversine formula
                lat1, lon1, lat2, lon2 = map(np.radians, [sub_lat, sub_lon, sub_lat - halfsub, sub_lon - halfsub])
                dlat = lat2 - lat1
                dlon = lon2 - lon1
                a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
                c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 
                distance = 6371 * c * 1000  # distance in meters
                # seek elevation at the centroid
                py, px = srtm_pixel.index(sub_lon, sub_lat)
                if not (0 <= px < srtm_pixel.width and 0 <= py < srtm_pixel.height):
                    print(f'[!] Error (SRMT-3). Coordinates are out of bounds! Lat: {sub_lat} Lon: {sub_lon}')
                    continue
                centroid_el = srtm_pixel.read(1, window=rasterio.windows.Window(px, py, 1, 1))[0, 0]
                if centroid_el > 2000:
                    drop_indices.append(index)
                    continue
                # check elevation at corners and compute slope
                drop_row = False
                for corner in corners:
                    py, px = srtm_pixel.index(corner[1], corner[0])
                    if not (0 <= px < srtm_pixel.width and 0 <= py < srtm_pixel.height):
                        print(f'[!] Error (SRMT-3). Coordinates are out of bounds! Lat: {corner[0]} Lon: {corner[1]}')
                        continue
                    corner_el = srtm_pixel.read(1, window=rasterio.windows.Window(px, py, 1, 1))[0, 0]
                    diff_el = abs(corner_el - centroid_el)
                    slope = diff_el / distance
                    if slope > max_slope:
                        drop_row = True
                        break
                if drop_row:
                    drop_indices.append(index)
                                    
        # drop rows from Dataframe and reset index
        mesh = mesh.drop(drop_indices).reset_index(drop=True)        

        # create an index for the dataframe
        mesh['id'] = range(1, len(mesh) + 1)
        return mesh


    # function to execute resource assessment and select the best available sites
    def res_assessment(self, net_mesh):
        # calculate areas for each subgrid
        areas = []
        for index, row in net_mesh.iterrows():
            # access geographical coordinates
            sub_lat = row['Sub lat']
            sub_lon = row['Sub lon']
            halfsub = self.era5_size / (subdim*2)
            # calculate distances using Haversine formula - lat
            lat1, lat2 = map(np.radians, [sub_lat + halfsub, sub_lat - halfsub])
            dlat = lat2 - lat1
            a_lat = np.sin(dlat/2)**2
            distance_lat = 2* 6371 * np.arctan2(np.sqrt(a_lat), np.sqrt(1-a_lat))       # distance in km
            # lon
            lon1, lon2 = map(np.radians, [sub_lon - halfsub, sub_lon + halfsub])
            lat_rad = np.radians(sub_lat)
            dlon = lon2 - lon1
            a_lon = np.cos(lat_rad)**2*np.sin(dlon/2)**2
            distance_lon = 2 * 6371 * np.arctan2(np.sqrt(a_lon), np.sqrt(1 - a_lon))    # distance in km
            # calculate the area of the subgrid
            area = distance_lat * distance_lon
            areas.append(area)
        # add the area column to the dataframe
        net_mesh['Area'] = areas

        return net_mesh

    # function to post process the gross net including exclusion zones
    def postproc(self):
        # create mesh using existing coordinates
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Creating mesh for subgrids at {current_time}")
        submesh = self.mesh_generation(subdim)

        # Exclude zones with geographical constraints WDPA, GlobCover, GPWv4, GEBCO, SRTM-3)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Eliminating areas not available for power generation {current_time}")
        net_mesh = self.exclusion_zone(submesh)

        # Select the subgrids eligible for resource assessment
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Select subgrids eligible for resource assessment {current_time}")
        top_mesh = self.res_assessment(net_mesh)
        
        # store output for reference
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Saving subgrids to perform simulation {current_time}")
        top_mesh.to_csv(self.cluster_subgrid, index=False)

    # function to select the top highest resources
    def seltop(self):
        # load net mesh with subgrids
        if os.path.exists(self.cluster_subgrid):
            net_mesh = pd.read_csv(self.cluster_subgrid)
        else:
            print(f"File not found: {self.cluster_subgrid}")
            sys.exit(1)

        # load global atlases to rank resources
        if os.path.exists(self.atlas):
            atlas_pixel = rasterio.open(self.atlas)
        else:
            print(f"File not found: {self.atlas}")
            sys.exit(1)

        # eliminate subgrids containing operating assets 
        if exop:
            # load operating assets
            if os.path.exists(self.cluster_act):
                act_asset = pd.read_csv(self.cluster_act, index_col=False)
            else:
                print(f"File not found: {self.cluster_act}")
                sys.exit(1)

            # eliminate subgrids with operating assets
            subsize = 0.25 / subdim
            geo_poly = net_mesh.apply(lambda row: box(
                row['Sub lon'] - subsize / 2, row['Sub lat'] - subsize / 2,
                row['Sub lon'] + subsize / 2, row['Sub lat'] + subsize / 2), axis=1)
            gdf_mesh = gpd.GeoDataFrame(net_mesh, geometry=geo_poly, crs='EPSG:4326')
            gdf_asset = gpd.GeoDataFrame(act_asset, geometry=gpd.points_from_xy(act_asset['lon'], act_asset['lat']), crs='EPSG:4326')
            gdf_mesh_exop = gpd.sjoin(gdf_mesh, gdf_asset, how='left', op='contains')
            gdf_mesh_exop = gdf_mesh_exop[gdf_mesh_exop['index_right'].isnull()]
            net_mesh = gdf_mesh_exop[['id','Grid lat', 'Grid lon', 'Sub lat', 'Sub lon', 'cluster', 'location','Area']]    # now it becomes pandas
            net_mesh.reset_index(drop=True, inplace=True)

        # iterate through the subgrids to define the average resource
        res_value = [np.nan] * len(net_mesh)
        for index, row in net_mesh.iterrows():
            # access geographical coordinates
            sub_lat = row['Sub lat']
            sub_lon = row['Sub lon']
            # Convert the geographical coordinates to dataset's pixel coordinates
            py, px = atlas_pixel.index(sub_lon, sub_lat)
            if not (0 <= px < atlas_pixel.width and 0 <= py < atlas_pixel.height):
                print(f'[!] Error (Atlas). Coordinates are out of bounds! Lat: {sub_lat} Lon: {sub_lon}')
                continue
            # window defines the data window we want to read (in this case, a single pixel)
            window = rasterio.windows.Window(px, py, 1, 1)
            band = atlas_pixel.read(1, window=window)
            value = band[0, 0]
            res_value[index] = value
        net_mesh['Resource'] = res_value

        # calculate top rankings for each cluster
        if climatevar==1:
            # define the first quartile of the resources
            for cluster_id, group in net_mesh.groupby('cluster'):
                cutoff = int(len(group) * 0.25)
                top_indices = group['Resource'].nlargest(cutoff).index
                net_mesh.loc[top_indices, 'Top'] = True
            
            # calculate area for top subgrids within each cluster
            top_clusters = net_mesh[net_mesh['Top'] == True]
            cluster_area = top_clusters.groupby('cluster')['Area'].sum().sort_values(ascending=False)
            cluster_labels = {cluster: f'SOL-{str(index+1).zfill(2)}' for index, cluster in enumerate(cluster_area.index)}
            net_mesh['Label'] = net_mesh['cluster'].map(cluster_labels)
        elif climatevar==2:
            # define the first quartile of the resources
            for (cluster_id, location), group in net_mesh.groupby(['cluster', 'location']):
                cutoff = int(len(group) * 0.25)
                top_indices = group['Resource'].nlargest(cutoff).index
                net_mesh.loc[top_indices, 'Top'] = True

            # calculate area for top subgrids within each cluster (segregating onshore and offshore)
            top_clusters = net_mesh[net_mesh['Top'] == True]
            cluster_area = top_clusters.groupby('cluster')['Area'].sum().sort_values(ascending=False)
            cluster_ranks = {cluster: idx + 1 for idx, cluster in enumerate(cluster_area.index)}
            net_mesh['Label'] = net_mesh.apply(lambda row: f"{'ONW' if row['location'] == 0 else 'OFW'}-{str(cluster_ranks[row['cluster']]).zfill(2)}", axis=1)

        # store output with new columns
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[+] Saving subgrids with selected subgrids {current_time}")
        net_mesh.to_csv(self.cluster_subgrid_sim, index=False)

    # plot a map with simulated subgrids
    def simsubgrid_view(self):
        # read data structure from the original clustering process
        if os.path.exists(self.clusterk):
            grid = pd.read_csv(self.clusterk)
        else:
            print(f"File not found: {self.clusterk}")
            sys.exit(1)

        # read data structure after pos processing
        if os.path.exists(self.cluster_subgrid_sim):
            subgrid = pd.read_csv(self.cluster_subgrid_sim)
        else:
            print(f"File not found: {self.cluster_subgrid_sim}")
            sys.exit(1)

        # create plot for a single k
        plt.figure(figsize=(3.5,3.4))   # single column

        # create a basemap for each subplot
        bmap = Basemap(projection='merc',llcrnrlat=-35,urcrnrlat=7, llcrnrlon=-77,urcrnrlon=-31,resolution='l')
        bmap.shadedrelief() # can also use map.etopo()

        # plotting clusters
        val_lat = grid['lat'].values
        val_lon = grid['lon'].values
        val_cls = grid['cluster'].values
        bmap.scatter(val_lon, val_lat, c=val_cls, label=val_cls, cmap='tab20c', latlon=True, marker='s', edgecolors='none', s=2.3, zorder=1)    # single column

        # plotting top subgrids
        top_areas = subgrid[subgrid['Top'] == True]
        top_lat = top_areas['Sub lat'].values
        top_lon = top_areas['Sub lon'].values
        bmap.scatter(top_lon, top_lat, color='black', latlon=True, marker='.', edgecolors='none', s=0.1, alpha=0.6, zorder=2, label='Selected subgrid')

        # adding map scale and legend
        bmap.readshapefile(shp_regional_arg1, shp_regional_arg2, drawbounds=True)
        for collection in plt.gca().collections:
            collection.set_linewidth(0.3)
        bmap.drawmapscale(-36, -32, -36, -32, 700, barstyle='simple', units='km', fontsize=8, labelstyle='simple', fillcolor1='b', fillcolor2='white')

        # add details about the cluster being displayed
        #plt.show()
        print("[+] Saving figure with clusters and subgrids")
        plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.savefig(self.fig4a, dpi=1000, format='jpg')   # Save the figure
        print("[+] Selected cluster done.")

    # plot a map with simulated subgrids
    def simgrid_view(self):
        # read data structure from the original clustering process
        if os.path.exists(self.clusterk):
            grid = pd.read_csv(self.clusterk)
        else:
            print(f"File not found: {self.clusterk}")
            sys.exit(1)

        # read cluster labels
        if os.path.exists(self.fig4_label):
            fig_labels = pd.read_csv(self.fig4_label)
        else:
            print(f"File not found: {self.fig4_label}")
            sys.exit(1)

        # create plot for a single k
        plt.figure(figsize=(3.5,3.4))   # single column

        # create a basemap for each subplot
        bmap = Basemap(projection='merc',llcrnrlat=-35,urcrnrlat=7, llcrnrlon=-77,urcrnrlon=-31,resolution='l')
        bmap.shadedrelief() # can also use map.etopo()

        # plotting clusters
        val_lat = grid['lat'].values
        val_lon = grid['lon'].values
        val_cls = grid['cluster'].values
        bmap.scatter(val_lon, val_lat, c=val_cls, label=val_cls, cmap='tab20c', latlon=True, marker='s', edgecolors='none', s=2.3, zorder=1)    # single column

        # add labels
        for _, row in fig_labels.iterrows():
            xpt, ypt = bmap(row['lon1'], row['lat1'])
            bmap.plot([xpt, bmap(row['lon2'], row['lat2'])[0]], [ypt, bmap(row['lon2'], row['lat2'])[1]], color='black', linewidth=0.5)
            plt.text(xpt, ypt, row['label'], fontsize=5, ha='center', color='black', bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.5'))

        # adding map scale and legend
        bmap.readshapefile(shp_regional_arg1, shp_regional_arg2, drawbounds=True)
        for collection in plt.gca().collections:
            collection.set_linewidth(0.3)
        bmap.drawmapscale(-36, -32, -36, -32, 700, barstyle='simple', units='km', fontsize=8, labelstyle='simple', fillcolor1='b', fillcolor2='white')

        # add details about the cluster being displayed
        #plt.show()
        print("[+] Saving figure with clusters with labels")
        plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        plt.savefig(self.fig4b, dpi=1000, format='jpg')   # Save the figure
        plt.savefig(self.pdf4b, format='pdf')
        print("[+] Selected cluster done.")

