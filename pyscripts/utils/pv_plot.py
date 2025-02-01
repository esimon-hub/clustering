"""

Plot comparisons between actual (data from the ISO) and simulation for each asset (hourly, daily, and monthly)

Description:
This script performs the following tasks:
1. Read data preprocessed with adjusted time series for actuals and simulated
2. Plot the data using a preset of charts
3. Plot statistics of main statistical parameters for hourly, daily and monthly timeframes

Usage:
1. Configure the parameters of simulation using the pv_stat_config.py
2. Run the current script (check the instances below)
3. Outputs will are displayed during the script execution

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
import numpy as np
import pandas as pd
import os
from pv_stat_config import exceptions, modebias, modeplt
from pv_stat_config import cl_iso_h, sin_iso_h                                          # ISO (Hourly)
from pv_stat_config import bc_cl_sim_h, bc_sin_sim_h                                    # BC  (Hourly)
from pv_stat_config import nbc_cl_sim_h, nbc_sin_sim_h                                  # NBC (Hourly)
from pv_stat_config import cl_iso_d, sin_iso_d                                          # ISO (Daily)
from pv_stat_config import bc_cl_sim_d, bc_sin_sim_d                                    # BC (Daily)
from pv_stat_config import nbc_cl_sim_d, nbc_sin_sim_d                                  # NBC (Daily)
from pv_stat_config import cl_iso_m, sin_iso_m                                          # ISO (Monthly)
from pv_stat_config import bc_cl_sim_m, bc_sin_sim_m                                    # BC (Monthly)
from pv_stat_config import nbc_cl_sim_m, nbc_sin_sim_m                                  # NBC (Monthly)
from pv_stat_config import fig_path_folder
from pv_stat_config import asset, beg, end
from pv_stat_config import minpts, stat_min
from pv_stat_config import stats_path_folder                                            # statistical analysis for subgrids
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
from calendar import monthrange
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
import math
from datetime import datetime

# begin
print("[+] Begin")

# Define the base directory where the figures will be saved
if not os.path.exists(fig_path_folder):
    os.makedirs(fig_path_folder)  # Create the directory if it doesn't exist

# read hourly dataframes needed to perform benchmark analysis
if modebias==1:
    # read hourly
    cl_iso_h=pd.read_csv(cl_iso_h, index_col=0, parse_dates=True)
    cl_sim_h=pd.read_csv(bc_cl_sim_h, index_col=0, parse_dates=True)
    sin_iso_h=pd.read_csv(sin_iso_h, index_col=0, parse_dates=True)
    sin_sim_h=pd.read_csv(bc_sin_sim_h, index_col=0, parse_dates=True)
    # read daily
    cl_iso_d=pd.read_csv(cl_iso_d, index_col=0, parse_dates=True)
    cl_sim_d=pd.read_csv(bc_cl_sim_d, index_col=0, parse_dates=True)
    sin_iso_d=pd.read_csv(sin_iso_d, index_col=0, parse_dates=True)
    sin_sim_d=pd.read_csv(bc_sin_sim_d, index_col=0, parse_dates=True)
    # read monthly
    cl_iso_m=pd.read_csv(cl_iso_m, index_col=0, parse_dates=True)
    cl_sim_m=pd.read_csv(bc_cl_sim_m, index_col=0, parse_dates=True)
    sin_iso_m=pd.read_csv(sin_iso_m, index_col=0, parse_dates=True)
    sin_sim_m=pd.read_csv(bc_sin_sim_m, index_col=0, parse_dates=True)
else:
    # read hourly
    cl_iso_h=pd.read_csv(cl_iso_h, index_col=0, parse_dates=True)
    cl_sim_h=pd.read_csv(nbc_cl_sim_h, index_col=0, parse_dates=True)
    sin_iso_h=pd.read_csv(sin_iso_h, index_col=0, parse_dates=True)
    sin_sim_h=pd.read_csv(nbc_sin_sim_h, index_col=0, parse_dates=True)
    # read daily
    cl_iso_d=pd.read_csv(cl_iso_d, index_col=0, parse_dates=True)
    cl_sim_d=pd.read_csv(nbc_cl_sim_d, index_col=0, parse_dates=True)
    sin_iso_d=pd.read_csv(sin_iso_d, index_col=0, parse_dates=True)
    sin_sim_d=pd.read_csv(nbc_sin_sim_d, index_col=0, parse_dates=True)
    # read monthly
    cl_iso_m=pd.read_csv(cl_iso_m, index_col=0, parse_dates=True)
    cl_sim_m=pd.read_csv(nbc_cl_sim_m, index_col=0, parse_dates=True)
    sin_iso_m=pd.read_csv(sin_iso_m, index_col=0, parse_dates=True)
    sin_sim_m=pd.read_csv(nbc_sin_sim_m, index_col=0, parse_dates=True)

# drop clusters with limited data availability and/or spurious data
exceptions_all = [column for column in cl_iso_h if cl_iso_h[column].count() < 1000]
exceptions_all.extend(exceptions)
for icol in exceptions_all:
    cl_iso_h.drop(columns=icol, inplace=True)
    cl_sim_h.drop(columns=icol, inplace=True)
    cl_iso_d.drop(columns=icol, inplace=True)
    cl_sim_d.drop(columns=icol, inplace=True)
    cl_iso_m.drop(columns=icol, inplace=True)
    cl_sim_m.drop(columns=icol, inplace=True)

# function to plot observed and simulated capacity factors
def plot_xy(cl_iso, cl_sim, tf):
    # remove outliers from the current cluster and plot a scatter
    clist=list(cl_iso.columns)
    nlist=len(clist)
    npage = True
    page = 0
    ielem = 0
    cols = 6
    rows = 3
    subplot_width = 7 / cols    # Width of each subplot
    subplot_height = 1.3       # Height of each subplot
    font_size_axis_labels = 7
    font_size_ticks = 7
    font_size_annotation = 7
    # create a dictionary to use an alias for each cluster and fit into the subplot
    n_map = {
        'Assú V': 'Assú V',
        'Conj. Alex': 'Alex',
        'Conj. BJL': 'BJL',
        'Conj. Boa Hora': 'Boa Hora',
        'Conj. Bom Jesus': 'Bom Jesus',
        'Conj. Brígida': 'Brígida',
        'Conj. Calcário': 'Calcário',
        'Conj. Dracena': 'Dracena',
        'Conj. FV SJP': 'FV SJP',
        'Conj. Floresta': 'Floresta',
        'Conj. Guaimbê': 'Guaimbê',
        'Conj. Horizonte': 'Horizonte',
        'Conj. Ituverava': 'Ituverava',
        'Conj. Jaíba': 'Jaíba',
        'Conj. Juazeiro Solar': 'Juazeiro S.',
        'Conj. Juazeiro Solar 2': 'Juazeiro S. 2',
        'Conj. Lapa': 'Lapa',
        'Conj. Lavras': 'Lavras',
        'Conj. Nova Olinda': 'N. Olinda',
        'Conj. Paracatu': 'Paracatu',
        'Conj. Pereira Barreto': 'P. Barreto',
        'Conj. Pirapora': 'Pirapora',
        'Conj. Rio Alto': 'Rio Alto',
        'Conj. Sertão Solar Barreiras': 'S. Barreiras',
        'Conj. Sol do Futuro': 'S. Futuro',
        'Conj. Sol do Sertão': 'S. Sertão',
        'Conj. São Gonçalo': 'S. Gonçalo',
        'Conj. São Pedro': 'S. Pedro'
    }
    while ielem < nlist:
        if npage:  # Starting a new page
            # Define the total number of plots for this page, considering remaining items
            fig_width = cols * subplot_width
            fig_height = rows * subplot_height  # Always use fixed number of rows for consistency
            fig, axs = plt.subplots(rows, cols, figsize=(fig_width, fig_height), squeeze=False, sharey=True)
            i, j = 0, 0  # Reset subplot indices
            npage = False  # Reset page flag

        while ielem < nlist and i < rows:  # Ensure we don't go beyond the designated rows and columns
            if j >= cols:  # If we've filled the row, move to the next
                i += 1
                j = 0
            if i >= rows:  # If we've filled the page, break to create a new page
                break
            
            iclst = clist[ielem]  # Cluster index
            alias = n_map.get(iclst,iclst)
            df_sct = pd.DataFrame({'obs': cl_iso[iclst], 'sim': cl_sim[iclst]}).dropna()
            filtered_df = df_sct[(df_sct['obs'] > 0.001) & (df_sct['sim'] > 0.001)]
            
            if not filtered_df.empty:  # Check if there is data
                # Plot only if there is enough data
                if tf==1:
                    ext_max=1
                elif tf==2:
                    ext_max=0.5
                elif tf==3:
                    ext_max=0.5


                # prepare configuration to plot
                colors = [(0, "white"), 
                          (0.3, "gray"), 
                          (1, "black")]
                cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

                hb = axs[i, j].hexbin(filtered_df['obs'], 
                                      filtered_df['sim'], 
                                      gridsize=15, 
                                      cmap=cmap, 
                                      mincnt=1, 
                                      extent=[0, ext_max, 0, ext_max]
                                      )
                axs[i, j].set_aspect('equal', adjustable='box')
                axs[i, j].plot([0, ext_max], [0, ext_max], 'r--')  # Plot y=x line
                axs[i, j].text(0.05, 0.95, f"{alias}", transform=axs[i, j].transAxes, verticalalignment='top', fontsize=font_size_annotation)
                axs[i, j].set_xticks(np.arange(0, ext_max+0.25, 0.5))
                axs[i, j].set_yticks(np.arange(0, ext_max+0.25, 0.5))
                axs[i, j].set_xlabel('Observed C.F.', fontsize=font_size_axis_labels)
                axs[i, j].set_ylabel('Simulated C.F.', fontsize=font_size_axis_labels)
                axs[i, j].tick_params(axis='both', which='major', labelsize=font_size_ticks)
                axs[i, j].grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
                j += 1  # Move to the next subplot in the row
            ielem += 1  # Move to the next data element

        # Hide any unused subplots on the current page
        for x in range(i, rows):  # Start from the last used row
            for y in range(j if x == i else 0, cols):  # Start from the next column if in the used row, else start from first column
                axs[x, y].axis('off')  # Hide the unused subplot
        npage = True  # Indicate new page setup for the next iteration
        page+=1
        plt.subplots_adjust(left=0.08, right=0.99, bottom=0.05, top=0.98, wspace=0.6, hspace=0.13)

        # save the files with higher resolution
        if tf==1:
            figure_name = f"pv_1_hourly_bench_{page}.png"  # Construct the filename
        elif tf==2:
            figure_name = f"pv_1_daily_bench_{page}.png"
        elif tf==3:
            figure_name = f"pv_1_monthly_bench_{page}.png"
        save_path = os.path.join(fig_path_folder, figure_name)  # Combine the directory with the filename
        plt.savefig(save_path, dpi=300)  # Save the figure

# function to plot observed and simulated capacity factors
def plot_xy_comb(sin_iso, sin_sim, tf):
    # plot definition
    width = 3.5                 # Width
    height = 3.0                # Height
    font_size_axis_labels = 7
    font_size_ticks = 7

    # clear the current figure
    plt.clf()  
    
    # plot only if there is enough data
    plt.figure(figsize=(width,height))

    # plot following tf
    if tf==1:
        ext_max=1
    elif tf==2:
        ext_max=0.5
    elif tf==3:
        ext_max=0.5

    df_sct = pd.DataFrame({'obs': sin_iso['cf'], 'sim': sin_sim['cf']}).dropna()
    filtered_df = df_sct[(df_sct['obs'] > 0.01) & (df_sct['sim'] > 0.01)]


    # prepare configuration to plot
    colors = [(0, "white"), 
                (0.3, "gray"), 
                (1, "black")]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

    if not filtered_df.empty:  # Check if there is data
        plt.hexbin(filtered_df['obs'], 
                   filtered_df['sim'], 
                   gridsize=15, 
                   cmap=cmap,
                   mincnt=1, 
                   extent=[0, ext_max, 0, ext_max])  # Adjust 'gridsize' and 'cmap'
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot([0,ext_max], [0,ext_max], 'r--')
        plt.xticks(np.arange(0,ext_max+0.25,0.25), fontsize=font_size_ticks)
        plt.yticks(np.arange(0,ext_max+0.25,0.25), fontsize=font_size_ticks)
        
        # Adding labels to both axes
        plt.xlabel('Observed C.F.', fontsize=font_size_axis_labels)  # X-axis label
        plt.ylabel('Simulated C.F.', fontsize=font_size_axis_labels)  # Y-axis label
        plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')     

        # adjust layout to avoid cutting off labels
        plt.tight_layout()

        # save the files with higher resolution
        if tf==1:
            figure_name = f"pv_2_hourly_combined.png"  # Construct the filename
        elif tf==2:
            figure_name = f"pv_2_daily_combined.png"
        elif tf==3:
            figure_name = f"pv_2_monthly_combined.png"
        save_path = os.path.join(fig_path_folder, figure_name)  # Combine the directory with the filename
        plt.savefig(save_path, dpi=300)  # Save the figure        
        # plot to check consistency
        #plt.show()  # Show the current page of plots

# function to plot observed and simulated capacity factors for clusters
def plot_ts(iso, sim, index, col):
    plt.clf()  # Clear the current figure
    width = 15              # width
    height = 7.0            # height    
    fig = plt.figure(figsize=(width, height))
    plt.text(0.05, 0.95, f"{col}", transform=plt.gca().transAxes, verticalalignment='top', fontsize=12)
    plt.plot(iso.index, iso, label='Observed')
    plt.plot(sim.index, sim, label='Simulated')
    plt.ylabel('Capacity Factor')    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')   
    plt.legend(loc='lower right')

    # adjust layout to avoid cutting off labels
    plt.tight_layout()

    figure_name = f"pv_3_cl_hourly_{index+1}"  # Construct the filename
    # combined
    if col=='sin':
        figure_name = f"pv_3_combined_hourly_{index+1}"
    # cluster
    save_path = os.path.join(fig_path_folder, figure_name)  # Combine the directory with the filename
    plt.savefig(save_path, dpi=300)  # Save the figure        
    #plt.show()  # Show the current page of plots
    plt.close(fig)

# function to plot observed and simulated capacity factors for clusters (side by side)
def plot_ts_2(sin_iso, sin_sim, index, col):
    plt.clf()  # Clear the current figure
    width = 15              # width
    height = 7.0            # height    
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))  # Create a figure and two side-by-side subplots

    # Plot the observed data in the first subplot
    axs[0].plot(sin_iso.index, sin_iso , label='Observed')
    axs[0].set_title(f"{col} - Observed")  # You can adjust the title as needed
    axs[0].set_ylabel('Capacity Factor')    
    axs[0].grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')

    # Plot the simulated data in the second subplot
    axs[1].plot(sin_sim.index, sin_sim, label='Simulated', color='orange')
    axs[1].set_title(f"{col} - Simulated")  # You can adjust the title as needed
    axs[1].set_ylabel('Capacity Factor')    
    axs[1].grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')

    # Adjust layout for better structure and saving
    plt.tight_layout()
    # save plot into a file
    figure_name = f"pv_3_cl_hourly_{index + 1}_2.png"  # Construct the filename
    save_path = os.path.join(fig_path_folder, figure_name)  # Combine the directory with the filename
    plt.savefig(save_path, dpi=300)  # Save the figure
    # plt.show()  # Uncomment if you want to display the figure as well
    plt.close(fig)

# function to plot observed and simulated capacity factors for clusters
def plot_inspect_ts(iso, sim, begin_date, end_date):
    plt.clf()  # Clear the current figure
    fig = plt.figure(figsize=(15, 7))
    date_range=pd.date_range(start=begin_date, end=end_date, freq='H')
    plt.plot(date_range, iso.reindex(date_range), label='Observed')
    plt.plot(date_range, sim.reindex(date_range), label='Simulated')
    plt.ylabel('Capacity Factor')    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')   
    plt.legend()
    plt.show()
    plt.close(fig)

# function to create boxplot for time series
def boxplot(iso, sim_bc, sim_nbc):
    # plot definition
    width = 3.5                 # width
    height = 3.0                # height
    font_size_axis_labels = 7
    font_size_ticks = 7

    # initialize dataframes to store evaluation metrics for each cluster
    metrics = ['Correlation','RMSE','MBE']
    results = {metric: pd.DataFrame() for metric in metrics}

    # calculate statistics for each cluster
    for cluster in iso.columns:
        print(cluster)
        # prepare data by removing NA values and ensuring sufficient data points
        valid_data = pd.DataFrame({
            'ISO': iso[cluster],
            'BC': sim_bc[cluster],
            'NBC': sim_nbc[cluster]
        }).dropna()

        if len(valid_data)>minpts:
            # compute evaluation metrics
            results['Correlation'][cluster] = [
                valid_data['ISO'].corr(valid_data['BC']),
                valid_data['ISO'].corr(valid_data['NBC'])
            ]
            results['RMSE'][cluster] = [
                np.sqrt(((valid_data['ISO'] - valid_data['BC']) ** 2).mean()),
                np.sqrt(((valid_data['ISO'] - valid_data['NBC']) ** 2).mean())
            ]
            results['MBE'][cluster] = [
                (valid_data['ISO'] - valid_data['BC']).mean(),
                (valid_data['ISO'] - valid_data['NBC']).mean()
            ]
        
    # plotting setup
    fig, axes = plt.subplots(nrows=3, figsize=(width,height*3), constrained_layout=True)
    y_labels = ['Correlation Coefficient', 'RMSE', 'MBE']
    # boxplot settings and labels
    for ax, metric, y_label in zip(axes, metrics, y_labels):
        ax.boxplot(results[metric].transpose(), widths=0.7, patch_artist=True, medianprops=dict(color='black'))
        ax.set(ylabel=y_label)
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=font_size_ticks)
        ax.set_aspect('auto', adjustable='box')

    # adjust bottom plot labels
    axes[-1].set_xticklabels(['BC', 'NBC'], fontsize=font_size_ticks)
    #axes[0].set_xticklabels(['', ''])
    #axes[1].set_xticklabels(['', ''])
    #axes[2].set_xticklabels(['BC', 'NBC'], fontsize=10)

    # save plot into a file
    figure_name = f"pv_4_boxplot.png"  # Construct the filename
    save_path = os.path.join(fig_path_folder, figure_name)  # Combine the directory with the filename
    plt.savefig(save_path, dpi=300)  # Save the figure
    # plt.show()  # Uncomment if you want to display the figure as well
    plt.close(fig)
    # display plots
    #plt.show()

# function to create boxplot for time series
def boxplot_all(labels, iso, sim1, sim2, sim3, sim4, sim5):
    # plot definition
    width = 3.5                     # width
    height = 1.5                    # height
    font_size_ticks = 8
    font_size_labels = 8

    # initialize dataframes to store evaluation metrics for each cluster
    metrics = ['Correlation','RMSE','MBE']
    results = {metric: pd.DataFrame() for metric in metrics}

    # calculate statistics for each cluster
    for cluster in iso.columns:
        print(cluster)
        # prepare data by removing NA values and ensuring sufficient data points
        valid_data = pd.DataFrame({
            'ISO': iso[cluster],
            'BC1': sim1[cluster],
            'BC2': sim2[cluster],
            'BC3': sim3[cluster],
            'BC4': sim4[cluster],
            'BC5': sim5[cluster]
        }).dropna()

        if len(valid_data)>minpts:
            # compute evaluation metrics
            results['Correlation'][cluster] = [
                valid_data['ISO'].corr(valid_data['BC1']),
                valid_data['ISO'].corr(valid_data['BC2']),
                valid_data['ISO'].corr(valid_data['BC3']),
                valid_data['ISO'].corr(valid_data['BC4']),
                valid_data['ISO'].corr(valid_data['BC5'])
            ]
            results['RMSE'][cluster] = [
                np.sqrt(((valid_data['BC1'] - valid_data['ISO']) ** 2).mean()),
                np.sqrt(((valid_data['BC2'] - valid_data['ISO']) ** 2).mean()),
                np.sqrt(((valid_data['BC3'] - valid_data['ISO']) ** 2).mean()),
                np.sqrt(((valid_data['BC4'] - valid_data['ISO']) ** 2).mean()),
                np.sqrt(((valid_data['BC5'] - valid_data['ISO']) ** 2).mean())
            ]
            results['MBE'][cluster] = [
                (valid_data['BC1'] - valid_data['ISO']).mean(),
                (valid_data['BC2'] - valid_data['ISO']).mean(),
                (valid_data['BC3'] - valid_data['ISO']).mean(),
                (valid_data['BC4'] - valid_data['ISO']).mean(),
                (valid_data['BC5'] - valid_data['ISO']).mean()
            ]
        
    # plotting setup
    fig, axes = plt.subplots(nrows=3, figsize=(width,height*3), constrained_layout=True, sharex=True)
    y_labels = ['Correlation Coefficient', 'RMSE', 'MBE']
    y_limits = [(0.77,1.0), (0.05,0.27),(-0.05,0.08)]

    # boxplot settings and labels
    for ax, metric, y_label, y_lim in zip(axes, metrics, y_labels, y_limits):
        box = ax.boxplot(results[metric].transpose(), widths=0.7, patch_artist=True, 
                         medianprops=dict(color='Khaki', linewidth=2),
                         flierprops={'marker': 'o', 'color': 'red', 'markersize': 4})
        for patch in box['boxes']:
            patch.set(facecolor='#1f77b4')  # Blue color for the box
        ax.set_ylabel(y_label, fontsize=font_size_labels, labelpad=8)
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=font_size_ticks)
        ax.set_ylim(y_lim)
        ax.set_aspect('auto', adjustable='box')

    # adjust bottom plot labels
    tick_positions = np.arange(1, len(labels) + 1)
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(labels, fontsize=font_size_ticks)

    # save plot into a file
    figure_name_png = f"solar_5_boxplot_all.png"  # Construct the filename
    figure_name_pdf = f"solar_5_boxplot_all.pdf"
    save_path_png = os.path.join(fig_path_folder, figure_name_png)  # Combine the directory with the filename
    save_path_pdf = os.path.join(fig_path_folder, figure_name_pdf)
    plt.savefig(save_path_png, dpi=1000)
    plt.savefig(save_path_pdf, format='pdf')
    plt.close(fig)

    # save csv files
    save_path1 = os.path.join(stats_path_folder, 'boxplot_corr.csv')
    save_path2 = os.path.join(stats_path_folder, 'boxplot_rmse.csv')
    save_path3 = os.path.join(stats_path_folder, 'boxplot_mbe.csv')
    results['Correlation'].to_csv(save_path1, index=False, encoding="utf-8-sig")
    results['RMSE'].to_csv(save_path2, index=False, encoding="utf-8-sig")
    results['MBE'].to_csv(save_path3, index=False, encoding="utf-8-sig")

# Function to calculate statistics from arguments
def stats_calc(tf, iso,sim_h1,sim_h2,sim_h3,sim_h4,sim_h5):
    # create dataframe from function arguments
    df_eval = pd.DataFrame({
        'iso': iso,
        'sim_h1': sim_h1,
        'sim_h2': sim_h2,
        'sim_h3': sim_h3,
        'sim_h4': sim_h4,
        'sim_h5': sim_h5
    })

    # drop NA elements to create statistics and avoid error messages; merge clean dataframe
    df_eval.dropna(inplace=True)
    # initialize an empty dataframe for results
    sin_stats = pd.DataFrame()
    # ensure there is a minimum amount of datapoints
    countval = len(df_eval[df_eval['iso'] > 0.01])
    if countval > stat_min:  # it has at least a minimum amount of datapoints
        # Calculate correlations
        correl_nbc_1 = df_eval['iso'].corr(df_eval['sim_h1'])
        correl_bc_2 = df_eval['iso'].corr(df_eval['sim_h2'])
        correl_bc_3 = df_eval['iso'].corr(df_eval['sim_h3'])
        correl_bc_4 = df_eval['iso'].corr(df_eval['sim_h4'])
        correl_bc_5 = df_eval['iso'].corr(df_eval['sim_h5'])

        # Calculate RMSE
        rmse_nbc_1 = ((df_eval['sim_h1'] - df_eval['iso'])**2).mean()**0.5
        rmse_bc_2 = ((df_eval['sim_h2'] - df_eval['iso'])**2).mean()**0.5
        rmse_bc_3 = ((df_eval['sim_h3'] - df_eval['iso'])**2).mean()**0.5
        rmse_bc_4 = ((df_eval['sim_h4'] - df_eval['iso'])**2).mean()**0.5
        rmse_bc_5 = ((df_eval['sim_h5'] - df_eval['iso'])**2).mean()**0.5
        
        # Calculate MBE
        mbe_nbc_1 = (df_eval['sim_h1'] - df_eval['iso']).mean()
        mbe_bc_2 = (df_eval['sim_h2'] - df_eval['iso']).mean()
        mbe_bc_3 = (df_eval['sim_h3'] - df_eval['iso']).mean()
        mbe_bc_4 = (df_eval['sim_h4'] - df_eval['iso']).mean()
        mbe_bc_5 = (df_eval['sim_h5'] - df_eval['iso']).mean()
        
        # Store the values of correlation, RMSE, MBE in a DataFrame
        sin_stats = pd.DataFrame({
            'parameter': ['corr_1', 'corr_2', 'corr_3', 'corr_4', 'corr_5',
                          'R_1', 'R_2', 'R_3', 'R_4', 'R_5',
                          'rmse_1', 'rmse_2', 'rmse_3', 'rmse_4', 'rmse_5',
                          'mbe_1', 'mbe_2', 'mbe_3', 'mbe_4', 'mbe_5'],
            'values': [correl_nbc_1, correl_bc_2, correl_bc_3, correl_bc_4, correl_bc_5,
                       correl_nbc_1**2, correl_bc_2**2, correl_bc_3**2, correl_bc_4**2, correl_bc_5**2, 
                       rmse_nbc_1, rmse_bc_2, rmse_bc_3, rmse_bc_4, rmse_bc_5,  
                       mbe_nbc_1, mbe_bc_2, mbe_bc_3, mbe_bc_4, mbe_bc_5]
        })

        # save csv files
        if tf==1:
            csv_name = f"pv_hourly_stats.csv"  # Construct the filename
        elif tf==2:
            csv_name = f"pv_daily_stats.csv"
        elif tf==3:
            csv_name = f"pv_monthly_stats.csv"
        save_path = os.path.join(stats_path_folder, csv_name)  # Combine the directory with the filename
        sin_stats.to_csv(save_path, index=False)

# Function to print statistics previously calculated
def print_stats(sin_iso,bc_sin_sim,stats):
    size1=len(sin_iso.index)
    size2=len(bc_sin_sim.index)
    print("\n\n+Number of elements in the series 1:"+str(size1)+" & 2:"+str(size2))
    # https://medium.com/@ottaviocalzone/mae-mse-rmse-and-f1-score-in-time-series-forecasting-d04021ffa7ce
    correl = str(stats[stats['parameter']=='correl_bc']['values'].iloc[0])
    rmse = str(stats[stats['parameter']=='correl_bc']['values'].iloc[0])
    mbe = str(stats[stats['parameter']=='correl_bc']['values'].iloc[0])
    rsquared = correl**2
    print("[+] correlation bc (Pearson): "+str(correl))
    print("[+] RMSE bc: "+str(rmse))
    print("[+] MBE bc: "+str(mbe))
    print("[+] r-squared bc: "+str(rsquared**2))

#
#  EXECUTION
#
# Hexbin
if modeplt==1:
    # plot xy charts
    print("[+] Hourly XY charts")
    plot_xy(cl_iso_h,cl_sim_h,1)
    print("[+] Daily XY charts")
    plot_xy(cl_iso_d,cl_sim_d,2)
    print("[+] Monthly XY charts")
    plot_xy(cl_iso_m,cl_sim_m,3)
    # plot xy combined charts
    print("[+] Hourly combined XY charts")
    plot_xy_comb(sin_iso_h,sin_sim_h,1)
    print("[+] Daily combined XY charts")
    plot_xy_comb(sin_iso_d,sin_sim_d,2)
    print("[+] Monthly combined XY charts")
    plot_xy_comb(sin_iso_m,sin_sim_m,3)

# Hourly time series (clusters and combined output)
if modeplt==2:
    for index, col in enumerate(cl_iso_h.columns):
        print(f"[+] Hourly time series: {col}")
        plot_ts(cl_iso_h[col],cl_sim_h[col],index,col)
        plot_ts_2(cl_iso_h[col],cl_sim_h[col],index,col)
    # plot sin
    plot_ts(sin_iso_h['cf'],sin_sim_h['cf'],0,'sin')

# Selected hourly time series (clusters and combined output)
if modeplt==3:
    if asset=='sin':
        iso=sin_iso_h['cf']
        sim=sin_sim_h['cf']
    else:
        iso=cl_iso_h[asset]
        sim=cl_sim_h[asset]
    # plot the selected time frame
    plot_inspect_ts(iso, sim, beg, end)

# Create a boxplot for the projects
if modeplt==4:
    # this module needs to read both biased and unbiased data (iso is already loaded)
    iso_h = cl_iso_h
    bc_sim_h=pd.read_csv(bc_cl_sim_h, index_col=0, parse_dates=True)
    nbc_sim_h=pd.read_csv(nbc_cl_sim_h, index_col=0, parse_dates=True)
    # function calculates statistics and create boxplot
    boxplot(cl_iso_h, bc_sim_h, nbc_sim_h)

# Create a boxplot for the projects
if modeplt==5:
    # this module needs to read both biased and unbiased data (iso is already loaded)
    script_dir = os.path.dirname(__file__)

    cl_iso_h_file=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_cl_iso.csv')
    cl_sim_h1_file=os.path.join(script_dir,'data/Solar/benchmark/NBC/h_cl_sim.csv')
    cl_sim_h2_file=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_cl_sim.csv')
    cl_sim_h3_file=os.path.join(script_dir,'data/Solar/benchmark/BC_30_C/h_cl_sim.csv')
    cl_sim_h4_file=os.path.join(script_dir,'data/Solar/benchmark/BC_10_C/h_cl_sim.csv')
    cl_sim_h5_file=os.path.join(script_dir,'data/Solar/benchmark/BC_5_C/h_cl_sim.csv')

    cl_iso_h=pd.read_csv(cl_iso_h_file, index_col=0, parse_dates=True)
    nbc_sim_h1=pd.read_csv(cl_sim_h1_file, index_col=0, parse_dates=True)
    bc_sim_h2=pd.read_csv(cl_sim_h2_file, index_col=0, parse_dates=True)
    bc_sim_h3=pd.read_csv(cl_sim_h3_file, index_col=0, parse_dates=True)
    bc_sim_h4=pd.read_csv(cl_sim_h4_file, index_col=0, parse_dates=True)
    bc_sim_h5=pd.read_csv(cl_sim_h5_file, index_col=0, parse_dates=True)

    labels=['NBC','Asset','1:900','1:100','1:25']

    # function calculates statistics and create boxplot
    boxplot_all(labels, cl_iso_h, nbc_sim_h1, bc_sim_h2, bc_sim_h3, bc_sim_h4, bc_sim_h5)

# Create and plot statistics for BC and NBC datasets (combined assets)
if modeplt==6:
    # define path for hourly, daily, and monthly
    script_dir = os.path.dirname(__file__)

    # read and pre-process country-level files
    sys_iso_h_file=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_sin_iso.csv')
    sys_sim_h1_file=os.path.join(script_dir,'data/Solar/benchmark/NBC/h_sin_sim.csv')
    sys_sim_h2_file=os.path.join(script_dir,'data/Solar/benchmark/BC_1/h_sin_sim.csv')
    sys_sim_h3_file=os.path.join(script_dir,'data/Solar/benchmark/BC_30_C/h_sin_sim.csv')
    sys_sim_h4_file=os.path.join(script_dir,'data/Solar/benchmark/BC_10_C/h_sin_sim.csv')
    sys_sim_h5_file=os.path.join(script_dir,'data/Solar/benchmark/BC_5_C/h_sin_sim.csv')

    sys_iso_h=pd.read_csv(sys_iso_h_file, index_col=0, parse_dates=True)
    nbc_sim_h1=pd.read_csv(sys_sim_h1_file, index_col=0, parse_dates=True)
    bc_sim_h2=pd.read_csv(sys_sim_h2_file, index_col=0, parse_dates=True)
    bc_sim_h3=pd.read_csv(sys_sim_h3_file, index_col=0, parse_dates=True)
    bc_sim_h4=pd.read_csv(sys_sim_h4_file, index_col=0, parse_dates=True)
    bc_sim_h5=pd.read_csv(sys_sim_h5_file, index_col=0, parse_dates=True)

    # call function to calculate statistics
    print('[+] Plotting statistics for hourly benchmarks')
    sin_hourly_stats = stats_calc(1, sys_iso_h['cf'], nbc_sim_h1['cf'], bc_sim_h2['cf'], bc_sim_h3['cf'], bc_sim_h4['cf'], bc_sim_h5['cf'])

print("[+] Done")

