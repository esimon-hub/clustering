"""
Renewable Resource Assessment 

Description:
This class includes the methods needed to perform the solar pv and wind resource assessment.

Usage:
Import this into the resource_assessment_exe.py to create instances.

For support or inquiries, contact emanueel@gmail.com

"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as patches
import numpy as np
import seaborn as sns
import networkx as nx
from matplotlib.ticker import FuncFormatter
from matplotlib import cm
from matplotlib import colors

#   Importing libraries and settings
from common.main_config import solar_subgrid_cf, wind_subgrid_on_cf, wind_subgrid_off_cf                        # input data for cf
from common.main_config import solar_subgrid_resass, wind_subgrid_on_resass, wind_subgrid_off_resass            # input data for cf
from common.main_config import solar_kna, onshore_kna, offshore_kna                                             # exclude these clusters from the assessement
from common.main_config import solar_cluster_fig5, wind_cluster_fig5a, wind_cluster_fig5b                       # subplot (cf) for capacity factors
from common.main_config import solar_cluster_pdf5, wind_cluster_pdf5a, wind_cluster_pdf5b                       # pdf
from common.main_config import solar_cluster_fig6, wind_cluster_fig6a, wind_cluster_fig6b                       # subplot for resource assessment
from common.main_config import solar_cluster_pdf6, wind_cluster_pdf6a, wind_cluster_pdf6b                       # pdf
from common.main_config import solar_table_res, wind_table_onshore_res, wind_table_offshore_res                 # summary table with resources
from common.main_config import solar_dist_res, wind_onshore_dist_res, wind_offshore_dist_res                    # summary table with resources
from common.main_config import solar_density, onwind_density, offwind_density                                   # define installed capacity density
from common.main_config import solar_subgrid_sim, wind_subgrid_sim                                              # calculate subgrid area by cluster
from common.main_config import solar_singlek, wind_singlek                                                      # calculate total (gross) area by cluster 
from common.main_config import assess_summary_table_res, assess_corr_table, assess_corr_chart, assess_corr_pdf  # summary table and correlation table
from common.main_config import ncols_boxplot                                                                    # number of columns for boxplot

# class to perform solar pv and wind resource assessment
class Resource():
    def __init__(self):
        print("[+] Reading input files for resource assessment")
        # load files with capacity factors for solar pv
        solar_cf_ts = pd.read_csv(solar_subgrid_cf, parse_dates=['datetime'])
        solar_cf_ts.set_index('datetime', inplace=True)
        solar_cf_ts.index = solar_cf_ts.index - pd.Timedelta(minutes=30)
        self.solar_cf = solar_cf_ts.rolling(window=2).mean()
        self.solar_cf = self.solar_cf[self.solar_cf.index.minute == 0]
        self.solar_cf.dropna(inplace=True)
        nrow = self.solar_cf.iloc[[0]].copy()
        nrow_index = self.solar_cf.index[0] - pd.Timedelta(hours=1)
        nrow.index = [nrow_index]
        self.solar_cf = pd.concat([nrow, self.solar_cf])
        self.solar_cf.index = self.solar_cf.index - pd.Timedelta(hours=3)
        self.solar_cf.index.name = 'datetime'
        if solar_kna:
            self.solar_cf.drop(columns=solar_kna, inplace=True)

        # load files with capacity factors for onshore wind
        self.onshore_wind_cf = pd.read_csv(wind_subgrid_on_cf, parse_dates=['datetime'])
        self.onshore_wind_cf.set_index('datetime', inplace=True)
        self.onshore_wind_cf.index = self.onshore_wind_cf.index - pd.Timedelta(hours=3)
        if onshore_kna:
            self.onshore_wind_cf.drop(columns=onshore_kna, inplace=True)

        # load files with capacity factors for offshore wind
        self.offshore_wind_cf = pd.read_csv(wind_subgrid_off_cf, parse_dates=['datetime'])
        self.offshore_wind_cf.set_index('datetime', inplace=True)
        self.offshore_wind_cf.index = self.offshore_wind_cf.index - pd.Timedelta(hours=3)
        if offshore_kna:
            self.offshore_wind_cf.drop(columns=offshore_kna, inplace=True)

        # load files with resource mapping for solar pv
        self.solar_res = pd.read_csv(solar_subgrid_resass, index_col=False)
        if solar_kna:
            self.solar_res = self.solar_res[~self.solar_res['Cluster'].isin(solar_kna)]

        # load files with resource mapping for onshore wind
        self.onshore_wind_res = pd.read_csv(wind_subgrid_on_resass, index_col=False)
        if onshore_kna:
            self.onshore_wind_res = self.onshore_wind_res[~self.onshore_wind_res['Cluster'].isin(onshore_kna)]

        # load files with resource mapping for offshore wind
        self.offshore_wind_res = pd.read_csv(wind_subgrid_off_resass, index_col=False)
        if offshore_kna:
            self.offshore_wind_res = self.offshore_wind_res[~self.offshore_wind_res['Cluster'].isin(offshore_kna)]        

        # load geographical stats (area) for solar pv (clusters and total)
        # calculate area for all clusters
        solar_tot_cl_areas = pd.read_csv(solar_singlek)
        solar_tot_cl_areas.dropna(subset=['cluster'], inplace=True)
        area_tot = 0
        for _, row in solar_tot_cl_areas.iterrows():
            lat = row['lat']
            lon = row['lon']
            halfsub = 0.25 / 2
            # calculate distances using Haversine formula - lat
            lat1, lat2 = map(np.radians, [lat + halfsub, lat - halfsub])
            dlat = lat2 - lat1
            a_lat = np.sin(dlat/2)**2
            distance_lat = 2* 6371 * np.arctan2(np.sqrt(a_lat), np.sqrt(1-a_lat))       # distance in km
            lon1, lon2 = map(np.radians, [lon - halfsub, lon + halfsub])
            lat_rad = np.radians(lat)
            dlon = lon2 - lon1
            a_lon = np.cos(lat_rad)**2*np.sin(dlon/2)**2
            distance_lon = 2 * 6371 * np.arctan2(np.sqrt(a_lon), np.sqrt(1 - a_lon))    # distance in km
            # calculate area of the grid
            area = distance_lat * distance_lon
            area_tot+=area
        # calculate area for subgrids
        solar_areas = pd.read_csv(solar_subgrid_sim, index_col=False)
        # stats for all areas within each cluster
        if solar_areas.empty:
            print("[!] No subgrids available for the simulation (all areas)!")
            sys.exit(1)
        all_areas_avrg = solar_areas.groupby('Label')['Area'].mean()
        all_counts = solar_areas.groupby('Label').size()
        all_area = all_areas_avrg * all_counts
        # stats for top areas within each cluster
        areas_top = solar_areas[solar_areas['Top'] == True]
        if areas_top.empty:
            print("[!] No subgrids available for the simulation (selected subgrids)!")
            sys.exit(1)
        top_areas_avrg = areas_top.groupby('Label')['Area'].mean()
        top_counts = areas_top.groupby('Label').size()
        top_area = top_areas_avrg * top_counts
        # create area dataframe for charts and tables
        self.solar_cluster_stats = pd.DataFrame({
            'Total_Clusters': area_tot,
            'All_Areas_Avrg': all_areas_avrg,
            'All_Counts': all_counts,
            'All_Area': all_area,
            'Top_Areas_Avrg': top_areas_avrg,
            'Top_Counts': top_counts,
            'Top_Area': top_area
        }).reset_index()

        # load geographical stats (area) for onshore wind (clusters and total)
        # calculate area for all clusters
        onshore_tot_cl_areas = pd.read_csv(wind_singlek)
        onshore_tot_cl_areas.dropna(subset=['cluster'], inplace=True)
        onshore_tot_cl_areas = onshore_tot_cl_areas[onshore_tot_cl_areas['location']==0]
        area_tot = 0
        for _, row in onshore_tot_cl_areas.iterrows():
            lat = row['lat']
            lon = row['lon']
            halfsub = 0.25 / 2
            # calculate distances using Haversine formula - lat
            lat1, lat2 = map(np.radians, [lat + halfsub, lat - halfsub])
            dlat = lat2 - lat1
            a_lat = np.sin(dlat/2)**2
            distance_lat = 2* 6371 * np.arctan2(np.sqrt(a_lat), np.sqrt(1-a_lat))       # distance in km
            lon1, lon2 = map(np.radians, [lon - halfsub, lon + halfsub])
            lat_rad = np.radians(lat)
            dlon = lon2 - lon1
            a_lon = np.cos(lat_rad)**2*np.sin(dlon/2)**2
            distance_lon = 2 * 6371 * np.arctan2(np.sqrt(a_lon), np.sqrt(1 - a_lon))    # distance in km
            # calculate area of the grid
            area = distance_lat * distance_lon
            area_tot+=area
        # calculate area for subgrids
        onshore_wind_areas = pd.read_csv(wind_subgrid_sim, index_col=False)
        onshore_wind_areas = onshore_wind_areas[onshore_wind_areas['location']==0.0]
        # stats for all areas within each cluster
        if onshore_wind_areas.empty:
            print("[!] No subgrids available for the simulation (all areas)!")
            sys.exit(1)
        all_areas_avrg = onshore_wind_areas.groupby('Label')['Area'].mean()
        all_counts = onshore_wind_areas.groupby('Label').size()
        all_area = all_areas_avrg * all_counts
        # stats for top areas within each cluster
        areas_top = onshore_wind_areas[onshore_wind_areas['Top'] == True]
        if areas_top.empty:
            print("[!] No subgrids available for the simulation (selected subgrids)!")
            sys.exit(1)
        top_areas_avrg = areas_top.groupby('Label')['Area'].mean()
        top_counts = areas_top.groupby('Label').size()
        top_area = top_areas_avrg * top_counts
        # create area dataframe for charts and tables
        self.onshore_wind_cluster_stats = pd.DataFrame({
            'Total_Clusters': area_tot,
            'All_Areas_Avrg': all_areas_avrg,
            'All_Counts': all_counts,
            'All_Area': all_area,
            'Top_Areas_Avrg': top_areas_avrg,
            'Top_Counts': top_counts,
            'Top_Area': top_area
        }).reset_index()

        # load geographical stats (area) for offshore wind (clusters and total)
        # calculate area for all clusters
        offshore_tot_cl_areas = pd.read_csv(wind_singlek)
        offshore_tot_cl_areas.dropna(subset=['cluster'], inplace=True)
        offshore_tot_cl_areas = offshore_tot_cl_areas[offshore_tot_cl_areas['location']==1]
        area_tot = 0
        for _, row in offshore_tot_cl_areas.iterrows():
            lat = row['lat']
            lon = row['lon']
            halfsub = 0.25 / 2
            # calculate distances using Haversine formula - lat
            lat1, lat2 = map(np.radians, [lat + halfsub, lat - halfsub])
            dlat = lat2 - lat1
            a_lat = np.sin(dlat/2)**2
            distance_lat = 2* 6371 * np.arctan2(np.sqrt(a_lat), np.sqrt(1-a_lat))       # distance in km
            lon1, lon2 = map(np.radians, [lon - halfsub, lon + halfsub])
            lat_rad = np.radians(lat)
            dlon = lon2 - lon1
            a_lon = np.cos(lat_rad)**2*np.sin(dlon/2)**2
            distance_lon = 2 * 6371 * np.arctan2(np.sqrt(a_lon), np.sqrt(1 - a_lon))    # distance in km
            # calculate area of the grid
            area = distance_lat * distance_lon
            area_tot+=area
        # calculate area for subgrids
        offshore_wind_areas = pd.read_csv(wind_subgrid_sim, index_col=False)
        offshore_wind_areas = offshore_wind_areas[offshore_wind_areas['location']==1.0]
        # stats for all areas within each cluster
        if offshore_wind_areas.empty:
            print("[!] No subgrids available for the simulation (all areas)!")
            sys.exit(1)
        all_areas_avrg = offshore_wind_areas.groupby('Label')['Area'].mean()
        all_counts = offshore_wind_areas.groupby('Label').size()
        all_area = all_areas_avrg * all_counts
        # stats for top areas within each cluster
        areas_top = offshore_wind_areas[offshore_wind_areas['Top'] == True]
        if areas_top.empty:
            print("[!] No subgrids available for the simulation (selected subgrids)!")
            sys.exit(1)
        top_areas_avrg = areas_top.groupby('Label')['Area'].mean()
        top_counts = areas_top.groupby('Label').size()
        top_area = top_areas_avrg * top_counts
        # create area dataframe for charts and tables
        self.offshore_wind_cluster_stats = pd.DataFrame({
            'Total_Clusters': area_tot,
            'All_Areas_Avrg': all_areas_avrg,
            'All_Counts': all_counts,
            'All_Area': all_area,
            'Top_Areas_Avrg': top_areas_avrg,
            'Top_Counts': top_counts,
            'Top_Area': top_area
        }).reset_index()

    # function to generate boxplots for solar pv and wind
    def plot_boxplots(self, res_cf, output, output_pdf):
        # group data by each hour of the day
        res_cf['hour'] = res_cf.index.hour
        grouped = res_cf.groupby('hour')

        # determine key parameters for subplot
        total_plots = len(res_cf.columns) - 1
        n_rows = (total_plots + ncols_boxplot - 1) // ncols_boxplot
        
        # plot setup
        fig_height = 1.5 * n_rows
        fig, axs = plt.subplots(n_rows, ncols_boxplot, figsize=(7, fig_height), sharex=False, sharey=True)
        axs = axs.ravel()

        # define reference hours to display on the x-axis
        reference_hours = [0, 6, 12, 18]
        reference_labels = [''] * 24  # Initialize all labels as empty
        for hour in reference_hours:
            reference_labels[hour] = str(hour)

        # adjust chart details
        for i, col in enumerate(res_cf.columns[:-1]): 
            data_to_plot = [group[col].dropna().values for name, group in grouped]
            box = axs[i].boxplot(data_to_plot, widths=0.5, showfliers=False, patch_artist=True, capprops={'visible': False}, medianprops={'linewidth': 0})
            for patch in box['boxes']:
                patch.set(facecolor='#1f77b4')  # Change to desired color
            axs[i].set_title(col, fontsize=8)

            # Plot median points (instead of using lines from boxplot)
            medians = [np.median(group[col].dropna().values) if not group[col].dropna().empty else 0 for name, group in grouped]
            axs[i].scatter(range(1, len(medians)+1), medians, color='Khaki', zorder=3, s=3)            

            # Y labels only on the first column
            if i % ncols_boxplot == 0:
                axs[i].set_ylabel('Capacity Factor', fontsize=8)

            # set X-axis Y-axis ticks and limits
            axs[i].set_xticks(range(1, 25))
            axs[i].set_xticklabels(reference_labels, fontsize=8, rotation=90)
            axs[i].set_ylim(0, 1.0)
            axs[i].set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            axs[i].tick_params(axis='y', labelsize=8)
        
        # hide unused axes if the number of projects is less than n_cols * n_rows
        for j in range(total_plots, len(axs)):
            axs[j].set_visible(False)        

        # store figure
        plt.tight_layout()
        plt.savefig(output, dpi=1000, format='jpg')
        plt.savefig(output_pdf, format='pdf')
        plt.show()

    # create bar chart with the distribution of resources by cluster and cf
    def resource_distribution(self, geo_res, geo_stats, output1, output2, output3, caparea, min, max, delta):
        # create bins to associate with the ranges
        num_steps = int((max - min) / delta) + 1
        bins = np.linspace(min, max, num_steps).tolist()
        bins[-1] = np.inf
        labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins) - 2)]
        labels.append(f">={bins[-2]:.2f}")
        geo_res['range'] = pd.cut(geo_res['CF'], bins=bins, labels=labels, right=False) # assign data to bins

        # group by clusters and cf range
        grouped_cf = geo_res.groupby(['Cluster', 'range']).size().unstack(fill_value=0)

        # create column with the average area by subgrid
        subgrid_area = geo_stats[['Label', 'Top_Areas_Avrg']]
        subgrid_area = subgrid_area.rename(columns={'Label': 'Cluster'})

        # merge area with grouped data and multiply each column by its area
        grouped_cf = grouped_cf.merge(subgrid_area, on='Cluster', how='left')
        for col in labels:
            grouped_cf[col] = grouped_cf[col] * grouped_cf['Top_Areas_Avrg'] * caparea * 10**-3
        grouped_cf = grouped_cf.set_index('Cluster').drop(columns=['Top_Areas_Avrg'])

        # calculate total potential for each cluster
        grouped_cf['total'] = grouped_cf.sum(axis=1)

        # sort clusters by descending order
        grouped_cf = grouped_cf.sort_values('total', ascending=False).drop(columns=['total'])
        
        # plot as stacked bar chart
        plt.figure(figsize=(3.5,5))
        ax = grouped_cf.plot(kind='bar', stacked=True, figsize=(7, 5), colormap='cividis')

        # save distribution to csv
        grouped_cf.set_index(labels)
        grouped_cf.to_csv(output1)

        # create plot with the data
        plt.xlabel('')
        plt.ylabel('Potential capacity (GW)', fontsize=12)
        plt.xticks(rotation=45)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x:,.0f}'))

        # reserve the order of the legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title='Capacity factor', fontsize=12, loc='upper right')

        # save plot
        plt.tight_layout()
        plt.savefig(output2, dpi=1000, format='jpg')
        plt.savefig(output3, format='pdf')

    # create table with the distribution of resources by cluster
    def table_distribution(self, res_cf, res_geo_subgrid, output, caparea):
        # define a list with the clusters included for the analysis
        clusters_all = res_cf.columns.tolist()
        clusters = [col for col in clusters_all if col not in solar_kna]

        # calculate potential areas and capacity for each cluster (suitable and top)
        pot_area = []
        top_area = []
        cap_area = []
        cum_top = 0
        for cluster in clusters:
            sel_data = res_geo_subgrid[res_geo_subgrid['Label'] == cluster]
            if not sel_data.empty:
                pot_area_val = sel_data['All_Area'].values[0]
                top_area_val = sel_data['Top_Area'].values[0]
                cap_area_val = top_area_val * caparea * 10**-3
                cum_top+=top_area_val
            else:
                pot_area_val = np.nan
                top_area_val = np.nan
                cap_area_val = np.nan
            pot_area.append(pot_area_val)
            top_area.append(top_area_val)
            cap_area.append(cap_area_val)

        # calculate share of the cluster
        cnt_area = []
        for top in top_area:
            if cum_top > 0 and not np.isnan(top):
                cnt_area_val = top / cum_top * 100
            else:
                cnt_area.append = np.nan
            cnt_area.append(cnt_area_val)

        # calculate average capacity factor and standard deviation
        avrg = []
        stdev = []
        coeff_var = []
        for cluster in clusters:
            if cluster in res_cf.columns:
                avrg_val = res_cf[cluster].mean()
                stdev_val = res_cf[cluster].std()
            else:
                avrg_val = np.nan
                stdev_val = np.nan
            avrg.append(avrg_val)
            stdev.append(stdev_val)
            if avrg_val > 0:
                coeff_var_val = stdev_val/avrg_val
            else:
                coeff_var_val = np.nan
            coeff_var.append(coeff_var_val)

        # calculate generation potential
        avrg = np.array(avrg)
        top_gen = cap_area * avrg * 8760 * 10**-3       # output in TWh

        # create and save dataframe with the calculated data
        cluster_distribution = pd.DataFrame({
            'Cluster': clusters,
            'Pot_Area (km2)': pot_area,
            'Top_Area (km2)': top_area,
            'Share of the cluster (%)': cnt_area,
            'Potential capacity (GW)': cap_area,
            'Potential generation (TWh/yr)': top_gen,
            'Average (CF)': avrg,
            'Standard Deviation': stdev,
            'Coefficient of Variation': coeff_var
        })
        cluster_distribution.to_csv(output, index=False)


    # create table with summary for all technologies
    def table_summary(self):
        # create table with summary for all clusters
        print("[+] Create summary table for all clusters")

        # select geographical stats for all clusters
        total_solar = self.solar_cluster_stats['Total_Clusters'][0]
        total_onshore = self.onshore_wind_cluster_stats['Total_Clusters'][0]
        total_offshore = self.offshore_wind_cluster_stats['Total_Clusters'][0]

        # create geographical stats for solar pv
        total_all_solar = self.solar_cluster_stats['All_Area'].sum()
        total_top_solar = self.solar_cluster_stats['Top_Area'].sum()

        # create geographical stats for onshore wind
        total_all_onshore_wind = self.onshore_wind_cluster_stats['All_Area'].sum()
        total_top_onshore_wind = self.onshore_wind_cluster_stats['Top_Area'].sum()

        # create geographical stats for offshore wind
        total_all_offshore_wind = self.offshore_wind_cluster_stats['All_Area'].sum()
        total_top_offshore_wind = self.offshore_wind_cluster_stats['Top_Area'].sum()

        # calculate suitable areas and prime areas in relation to the overall country's area
        label = ['Solar PV', 'Onshore wind', 'Offshore wind']
        total_cluster = [total_solar, total_onshore, total_offshore]
        total_area = [total_all_solar, total_all_onshore_wind, total_all_offshore_wind]
        total_area_pct = [total_all_solar/total_solar, total_all_onshore_wind/total_onshore, total_all_offshore_wind/total_offshore]
        top_area = [total_top_solar, total_top_onshore_wind, total_top_offshore_wind]
        top_area_pct = [total_top_solar/total_solar, total_top_onshore_wind/total_onshore, total_top_offshore_wind/total_offshore]

        # create dataframe and store it as a csv file
        table_summary = pd.DataFrame({
            'Technology': label,
            'Total cluster (km2)': total_cluster,
            'Suitable Areas (km2)': total_area,
            'Suitable Area Coverage (%)': total_area_pct,
            'Selected Areas (km2)': top_area,
            'Selected Area Coverage (%)': top_area_pct
        })
        table_summary.to_csv(assess_summary_table_res, index=False)
        

    # create chart and table for the correlation matrix
    def corr(self):
        #print("[+] Create correlation matrix for selected clusters")
        comb_clusters = pd.concat([self.solar_cf, self.onshore_wind_cf, self.offshore_wind_cf], axis=1)
        
        # correlation matrix
        corr_matrix = comb_clusters.corr()

        # determining the size of each category
        num_solar = self.solar_cf.shape[1]
        num_onshore_wind = self.onshore_wind_cf.shape[1]
        num_offshore_wind = self.offshore_wind_cf.shape[1]

        # Indices for boxes
        solar_index = (0, num_solar)
        onshore_wind_index = (num_solar, num_solar + num_onshore_wind)
        offshore_wind_index = (num_solar + num_onshore_wind, num_solar + num_onshore_wind + num_offshore_wind)        

        # set parameters for seaborn
        sns.set_theme()
        plt.rcParams['axes.labelsize'] = 8  # Adjusts the font size of the labels
        plt.rcParams['xtick.labelsize'] = 8  # Adjusts the font size of the x-tick labels
        plt.rcParams['ytick.labelsize'] = 8  # Adjusts the font size of the y-tick labels

        # define figure parameters
        plt.figure(figsize=(7, 8))
        ax = sns.heatmap(corr_matrix, annot=False, center=0, cmap='seismic', 
                    linewidths=.75, cbar=True, cbar_kws={'shrink': 0.5, 'location': 'bottom', 'pad': 0.09},
                    xticklabels=True, yticklabels=True)
        ax.set_xlim([-0.1, corr_matrix.shape[1]+0.1])
        ax.set_ylim([-0.1, corr_matrix.shape[0]+0.1])
        ax.invert_yaxis()        
        # solar pv
        ax.add_patch(patches.Rectangle((solar_index[0], solar_index[0]), 
                                    solar_index[1] - solar_index[0], solar_index[1] - solar_index[0], 
                                    fill=False, edgecolor='black', lw=1.5, linestyle='dashed'))
        # Onshore Wind
        ax.add_patch(patches.Rectangle((onshore_wind_index[0], onshore_wind_index[0]), 
                                    onshore_wind_index[1] - onshore_wind_index[0], onshore_wind_index[1] - onshore_wind_index[0], 
                                    fill=False, edgecolor='black', lw=1.5, linestyle='dashed'))
        # Offshore Wind
        ax.add_patch(patches.Rectangle((offshore_wind_index[0], offshore_wind_index[0]), 
                                    offshore_wind_index[1] - offshore_wind_index[0], offshore_wind_index[1] - offshore_wind_index[0], 
                                    fill=False, edgecolor='black', lw=1.5, linestyle='dashed'))

        plt.tight_layout()
        plt.savefig(assess_corr_chart, dpi=1000, format='jpg')   # Save the figure
        plt.savefig(assess_corr_pdf, format='pdf')
        # save correlation matrix to a csv file
        corr_matrix.to_csv(assess_corr_table, index=False)
        print("[+] Correlation matrix and chart created.")


    # produce resource assessment (boxplot and tables) for all technologies
    def resource_assessment(self):
        # boxplots for hourly generation
        print("[+] Plotting boxplots for clusters (capacity factors)")
        self.plot_boxplots(self.solar_cf, solar_cluster_fig5, solar_cluster_pdf5)
        self.plot_boxplots(self.onshore_wind_cf, wind_cluster_fig5a, wind_cluster_pdf5a)
        self.plot_boxplots(self.offshore_wind_cf, wind_cluster_fig5b, wind_cluster_pdf5b)

        # plot resource distribution by clusters
        print("[+] Plotting resource distribution")
        # solar pv
        min_cf = 0.22
        max_cf = 0.32
        del_cf = 0.02
        self.resource_distribution(self.solar_res, self.solar_cluster_stats, solar_dist_res, solar_cluster_fig6, solar_cluster_pdf6, solar_density, min_cf, max_cf, del_cf)
        # onshore wind
        min_cf = 0.25
        max_cf = 0.75
        del_cf = 0.1
        self.resource_distribution(self.onshore_wind_res, self.onshore_wind_cluster_stats, wind_onshore_dist_res, wind_cluster_fig6a, wind_cluster_pdf6a, onwind_density, min_cf, max_cf, del_cf)
        # offshore wind
        min_cf = 0.25
        max_cf = 0.75
        del_cf = 0.1
        self.resource_distribution(self.offshore_wind_res, self.offshore_wind_cluster_stats, wind_offshore_dist_res, wind_cluster_fig6b, wind_cluster_pdf6b, offwind_density, min_cf, max_cf, del_cf)
        
        # create table with summary for solar pv clusters
        print("[+] Create summary table for solar clusters")
        self.table_distribution(self.solar_cf, self.solar_cluster_stats, solar_table_res, solar_density)
        # create table with summary for onshore wind clusters
        print("[+] Create summary table for onshore wind clusters")
        self.table_distribution(self.onshore_wind_cf, self.onshore_wind_cluster_stats, wind_table_onshore_res, onwind_density)
        # create table with summary for onshore wind clusters
        print("[+] Create summary table for offshore wind clusters")
        self.table_distribution(self.offshore_wind_cf, self.offshore_wind_cluster_stats, wind_table_offshore_res, offwind_density)
