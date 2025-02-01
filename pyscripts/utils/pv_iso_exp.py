"""

Extract power generation data provided by the ISO into a Pandas DataFrame

Description:
This script performs the following tasks:
1. Read raw data downloaded from the ISO (ONS)
2. Replace data gaps and issues with NaN (maintaining a single index for all assets)

Usage:
1. Configure the parameters of simulation using the pv_stat_config.py
2. Run the current script (check the instances below)

Note: Ensure that all required dependencies are installed before running

For support or inquiries, contact emanueel@gmail.com

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob

from pv_stat_config import in_iso_path1, in_iso_path2_dir, in_window, in_assetexc, mode_exec 
from pv_stat_config import out_path, out_cf_path0
from pv_stat_config import sel_asset
from pv_stat_config import inspect_figure

# begin
print("[+] Begin")
# process raw data from ONS (adjust for timezone)
if mode_exec==1:
    print("[+] Executing mode 1")
    try:
        # read data and create dataframes
        parse_dts=['din_instante']
        iso_data1 = pd.read_csv(in_iso_path1, dtype={'val_capacinstaladacomissconmw': float, 'val_geraenergiaconmwmed': float}, parse_dates=parse_dts, index_col=False)
        iso_data1 = iso_data1[['din_instante','nom_usina','dsc_estado','Fator de Carga Comiss','val_capacinstaladacomissconmw','val_geraenergiaconmwmed']]
        files_path=glob.glob(os.path.join(in_iso_path2_dir,'*.csv'))
        iso_data2_list = [pd.read_csv(file, dtype={'val_capacinstaladacomissconmw': float, 'val_geraenergiaconmwmed': float}, delimiter=';', parse_dates=parse_dts, index_col=False) for file in files_path]
        iso_data2_yr = pd.concat(iso_data2_list, ignore_index=True)
        iso_data2_yr.rename(columns={'nom_estado': 'dsc_estado', 'nom_usina_conjunto': 'nom_usina', 'val_fatorcapacidade': 'Fator de Carga Comiss', 'val_capacidadeinstalada': 'val_capacinstaladacomissconmw', 'val_geracaoverificada': 'val_geraenergiaconmwmed'}, inplace=True)
        iso_data2_sel = iso_data2_yr[iso_data2_yr['nom_tipousina']=='Solar']
        iso_data2 = iso_data2_sel[['din_instante','nom_usina','dsc_estado','Fator de Carga Comiss','val_capacinstaladacomissconmw','val_geraenergiaconmwmed']]
        iso_data3 = pd.concat([iso_data1, iso_data2], ignore_index=True)
    except FileNotFoundError as e:
        print(f"[!] Error loading files: {e}")
        exit(1)

    # create new columns with day, month, and year to group entries
    iso_data3['day'] = pd.DatetimeIndex(iso_data3['din_instante']).day
    iso_data3['month'] = pd.DatetimeIndex(iso_data3['din_instante']).month
    iso_data3['year'] = pd.DatetimeIndex(iso_data3['din_instante']).year

    # select subset to store
    df1 = iso_data3[['din_instante', 'nom_usina', 'dsc_estado', 'day', 'month', 'year', 'Fator de Carga Comiss', 'val_capacinstaladacomissconmw', 'val_geraenergiaconmwmed']].copy()

    # ONS has cells without data (generation/cf) let's rid the dataset of empty cells
    df1['Fator de Carga Comiss'].replace('', np.nan, inplace=True)
    df1.dropna(subset=['Fator de Carga Comiss'], inplace=True)

    # DAYLIGHT SAVING TIME
    # ONS provides the historical generation figures using DST in a format that is not recognized by
    # pyzt. The data from the ISO overwrites the transition between DST and regular hours
    # The code below checks whether the time is DST because it is never DST in UTC
    # select the years where data will be benchmarked (years where ISO has data)
    # https://pt.wikipedia.org/wiki/Lista_de_per%C3%ADodos_em_que_vigorou_o_hor%C3%A1rio_de_ver%C3%A3o_no_Brasil
    # 1 datapoint will be missing 23 hour
    df1['din_instante'] = pd.to_datetime(df1['din_instante'])
    df1['din_instante_dst'] = pd.to_datetime(df1['din_instante']).dt.tz_localize('America/Sao_Paulo', ambiguous=False)
    df1['flag'] = df1['din_instante_dst'].map(lambda x : x.dst().total_seconds()!=0)
    for idx in df1.index:
         if df1.loc[idx,'flag'] == True:
            df1.loc[idx,'din_instante'] = df1.loc[idx,'din_instante'] - pd.Timedelta(hours=1)

    # store in the same folder
    df1.to_csv(out_path,index=False,encoding="utf-8-sig")

# preprocess entries (eliminate 1: commissioning periods; 2: CF greater than 1; 3: windows of spurious or unreliable data) and calculate means
if mode_exec==2:
    print("[+] Executing mode 2")
    try:
        # open the file from step 1 read data and create dataframes
        pos_data1 = pd.read_csv(out_path, index_col=False)
        window_out = pd.read_csv(in_window, index_col=False)
        asset_exc = pd.read_csv(in_assetexc, index_col=False)
    except FileNotFoundError as e:
        print(f"[!] Error loading files: {e}")
        exit(1)

    # ensure datetime column has the correct data type
    pos_data1['din_instante'] = pd.to_datetime(pos_data1['din_instante'])
    window_out['Begin'] = pd.to_datetime(window_out['Begin'])
    window_out['End'] = pd.to_datetime(window_out['End'])

    # eliminate commissioning period after shifts in installed capacity
    def extra_commissioning(group):
        # ensure the group is sorted by datetime
        group = group.sort_values('din_instante')
        # detect where installed capacity changes
        capchg = group['val_capacinstaladacomissconmw'] != group['val_capacinstaladacomissconmw'].shift()
        # for each change, identify the date 2 months after the change
        for datechg in group.loc[capchg,'din_instante']:
            # set the timeframe between step change and commissioning period
            cutoff_date = datechg + pd.DateOffset(months=2)
            # set cf to NaN for the calculated timeframe
            group.loc[(group['din_instante'] >= datechg) & (group['din_instante'] <= cutoff_date), 'Fator de Carga Comiss'] = np.nan
            group.loc[(group['din_instante'] >= datechg) & (group['din_instante'] <= cutoff_date), 'val_capacinstaladacomissconmw'] = np.nan
        return group
    # apply the function to each project group
    pos_data1 = pos_data1.groupby('nom_usina').apply(extra_commissioning).reset_index(drop=True)

    # eliminate entries with capacity factors greater than 1 (become NaN)
    pos_data1['Fator de Carga Comiss'] = np.where(pos_data1['Fator de Carga Comiss']>1,np.nan,pos_data1['Fator de Carga Comiss'])
    pos_data1['val_capacinstaladacomissconmw'] = np.where(pos_data1['Fator de Carga Comiss']>1,np.nan,pos_data1['val_capacinstaladacomissconmw'])

    # eliminate entries where cf remain <=0 for more than 48 consecutive hours
    def filter_cf(group):
        group = group.sort_values('din_instante')
        # identify rows where capacity factors is <= 0
        group['cf_spur'] = group['Fator de Carga Comiss'] <= 0
        # use cumsum to identify consecutive periods
        group['block'] = (group['cf_spur'] != group['cf_spur'].shift()).cumsum()
        # count consecutive hours in each block
        block_counts = group.groupby('block').size()
        # identify consecutive blocks of cf <= 0
        block_del = block_counts[block_counts > 48].index
        # filter out blocks with more than 48 hours of cf <= 0
        group.loc[group['block'].isin(block_del), 'Fator de Carga Comiss'] = np.nan
        group.loc[group['block'].isin(block_del), 'val_capacinstaladacomissconmw'] = np.nan
        group = group.drop(columns=['cf_spur', 'block'])
        return group
    pos_data1 = pos_data1.groupby('nom_usina').apply(filter_cf).reset_index(drop=True)

    # eliminate projects with no reported "Fator de Carga Comiss"
    def remove_unreported_cf(group):
        # If all 'Fator de Carga Comiss' values in the group are NaN, return None, else return the group
        if group['Fator de Carga Comiss'].isna().all():
            return None
        else:
            return group    
    pos_data1 = pos_data1.groupby('nom_usina').apply(remove_unreported_cf).reset_index(drop=True).dropna(subset=['nom_usina'])

    # eliminate windows of entries with inconsistent data
    for _, row in window_out.iterrows():
        asset_name = row['nom_usina']  # look for individual entries to eliminate outliers
        start_date, end_date = row['Begin'], row['End']  # begin and end dates
        # logic to set entries to NaN (cfs)
        pos_data1.loc[(pos_data1['nom_usina']==asset_name) &
                      (pos_data1['din_instante']>=start_date) &
                      (pos_data1['din_instante']<=end_date),
                      'Fator de Carga Comiss'] = np.nan
        # logic to set entries to NaN (capacity)
        pos_data1.loc[(pos_data1['nom_usina']==asset_name) &
                      (pos_data1['din_instante']>=start_date) &
                      (pos_data1['din_instante']<=end_date),
                      'val_capacinstaladacomissconmw'] = np.nan

    # exclude projects that include that have unreliable entries
    asset_sel = set(asset_exc['nom_usina'].unique())
    pos_data1 = pos_data1[~pos_data1['nom_usina'].isin(asset_sel)]

    # store output in the same folder
    pos_data1.to_csv(out_cf_path0,index=False,encoding="utf-8-sig")

# plot results to identify outliers
if mode_exec==3:
    print("[+] Executing mode 3") 
    try:
        # open the file from step 1 read data and create dataframes
        pos_data1 = pd.read_csv(out_path, index_col=False)
        pos_data2 = pd.read_csv(out_cf_path0, index_col=False)
    except FileNotFoundError as e:
        print(f"[!] Error loading files: {e}")
        exit(1)

    # ensure datetime
    pos_data1['din_instante'] = pd.to_datetime(pos_data1['din_instante'])
    pos_data2['din_instante'] = pd.to_datetime(pos_data2['din_instante'])

    # function to calculate the moving average
    def plot_comparison(project_name):
        # filter each dataframe by project name
        original_project = pos_data1[pos_data1['nom_usina'] == project_name]
        modified_project = pos_data2[pos_data2['nom_usina'] == project_name]

        # adjust to datetime format
        original_project['din_instante'] = pd.to_datetime(original_project['din_instante'])
        modified_project['din_instante'] = pd.to_datetime(modified_project['din_instante'])

        # ensure they are sorted by time and set din_instante as index for resampling
        original_project.sort_values('din_instante', inplace=True)
        modified_project.sort_values('din_instante', inplace=True)
        original_project.set_index('din_instante', inplace=True)
        modified_project.set_index('din_instante', inplace=True)

        # plot hourly data
        # Create a subplot with 2 rows and 1 column
        fig, axs = plt.subplots(2, 2, figsize=(15, 10), sharex=True)

        # Plot original data on the first subplot
        axs[0, 0].plot(original_project.index, original_project['Fator de Carga Comiss'], label='Original', alpha=0.5)
        axs[0, 0].set_title(f'Original CF time series for {project_name}')
        axs[0, 0].set_ylabel('Capacity Factor')
        axs[0, 0].legend()
        # Original installed capacity
        axs[1, 0].plot(original_project.index, original_project['val_capacinstaladacomissconmw'], label='Original', alpha=0.5)
        axs[1, 0].set_title(f'Original installed capacity for {project_name}')
        axs[1, 0].set_ylabel('MW')
        axs[1, 0].set_xlabel('Time')
        axs[1, 0].legend()

        # Plot modified data
        axs[0, 1].plot(modified_project.index, modified_project['Fator de Carga Comiss'], label='Modified', alpha=0.5)
        axs[0, 1].set_title(f'Modified CF time series for {project_name}')
        axs[0, 1].set_ylabel('Capacity Factor')
        axs[0, 1].legend()
        # Plot modified installed capacity
        axs[1, 1].plot(modified_project.index, modified_project['val_capacinstaladacomissconmw'], label='Modified installed capacity', alpha=0.5)
        axs[1, 1].set_title(f'Modified installed capacity for {project_name}')
        axs[1, 1].set_ylabel('MW')
        axs[1, 1].set_xlabel('Time')
        axs[1, 1].legend()

        # Automatically adjust subplot params so that the subplot(s) fits in to the figure area
        plt.tight_layout()
        plt.show()

    # plot selected projects
    plot_comparison(sel_asset)

# plot results for all projects 
if mode_exec==4:
    print("[+] Executing mode 4")
    try:
        # open the file from step 1 read data and create dataframes
        pos_data1 = pd.read_csv(out_path, index_col=False)
        pos_data2 = pd.read_csv(out_cf_path0, index_col=False)
    except FileNotFoundError as e:
        print(f"[!] Error loading files: {e}")
        exit(1)

    # ensure datetime
    pos_data1['din_instante'] = pd.to_datetime(pos_data1['din_instante'])
    pos_data2['din_instante'] = pd.to_datetime(pos_data2['din_instante'])

    # function to calculate the moving average
    def plot_comparison(project_name, index):
        # filter each dataframe by project name
        original_project = pos_data1[pos_data1['nom_usina'] == project_name]
        modified_project = pos_data2[pos_data2['nom_usina'] == project_name]

        # adjust to datetime format
        original_project['din_instante'] = pd.to_datetime(original_project['din_instante'])
        modified_project['din_instante'] = pd.to_datetime(modified_project['din_instante'])

        # ensure they are sorted by time and set din_instante as index for resampling
        original_project.sort_values('din_instante', inplace=True)
        modified_project.sort_values('din_instante', inplace=True)
        original_project.set_index('din_instante', inplace=True)
        modified_project.set_index('din_instante', inplace=True)

        # plot hourly data
        # Create a subplot with 2 rows and 1 column
        fig, axs = plt.subplots(2, 2, figsize=(15, 10), sharex=True)

        # Plot original data on the first subplot
        axs[0, 0].plot(original_project.index, original_project['Fator de Carga Comiss'], label='Original', alpha=0.5)
        axs[0, 0].set_title(f'Original CF time series for {project_name}')
        axs[0, 0].set_ylabel('Capacity Factor')
        axs[0, 0].legend()
        # Original installed capacity
        axs[1, 0].plot(original_project.index, original_project['val_capacinstaladacomissconmw'], label='Original', alpha=0.5)
        axs[1, 0].set_title(f'Original installed capacity for {project_name}')
        axs[1, 0].set_ylabel('MW')
        axs[1, 0].set_xlabel('Time')
        axs[1, 0].legend()

        # Plot modified data
        axs[0, 1].plot(modified_project.index, modified_project['Fator de Carga Comiss'], label='Modified', alpha=0.5)
        axs[0, 1].set_title(f'Modified CF time series for {project_name}')
        axs[0, 1].set_ylabel('Capacity Factor')
        axs[0, 1].legend()
        # Plot modified installed capacity
        axs[1, 1].plot(modified_project.index, modified_project['val_capacinstaladacomissconmw'], label='Modified installed capacity', alpha=0.5)
        axs[1, 1].set_title(f'Modified installed capacity for {project_name}')
        axs[1, 1].set_ylabel('MW')
        axs[1, 1].set_xlabel('Time')
        axs[1, 1].legend()

        # Automatically adjust subplot params so that the subplot(s) fits in to the figure area
        plt.tight_layout()
        figure_name = f"pv_inspect_{index}.png"   # Construct the filename
        save_path = os.path.join(inspect_figure, figure_name)  # Combine the directory with the filename
        plt.savefig(save_path, dpi=100)             # Save the figure        

    # plot for all projects
    proj_list = pos_data2['nom_usina'].unique()

    # iterate over these unique values and call prot_comparison
    for index, project in enumerate(proj_list):
        print(f'[+] Plotting project: {project}')
        plot_comparison(project, index)

print("[+] Done")
