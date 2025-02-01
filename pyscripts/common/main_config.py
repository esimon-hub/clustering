"""
SVD Parameters

Description:
This script sets the overall assumptions to perform SVD and clustering techniques:

Usage:
1. Set all parameters below before running the simulation
2. Make sure nc files follow the pattern CCC_YYYY_MM.nc (CCC=country; YYYY=year; MM=month)

For support or inquiries, contact emanueel@gmail.com

"""

import os
script_dir = os.path.dirname(__file__)

#
# svd_exe.py & clusters_exe.py
#
climatevar=2                                       # 1:solar; 2:wind
yr_svd1=1994                                        # read ERA5 files between these dates
yr_svd2=2023
max_clusters=40                                     # elbow chart
threshold=0.99                                      # required threshold for the cumulative energy
solar_k_subimg=[1,5,10,15,20,25,30,35]              # number of clusters in the subplot (ready for solar radiation)
wind_k_subimg=[1,5,10,15,20,25,30,35]               # number of clusters in the subplot (ready for wind speeds)
k_single=16                                         # single plot for a given k
subdim=10                                           # dimension subgrid NxN (needs to be consistent with simulation)
nupper=75                                           # calculate the Nth percentile
exop=True                                          # exclude operating assets from clustering process
minvals=True                                        # set minimum values for the climate variable to run the SVD
shape=True                                          # shape (adimensional) or Absolute values (m/2)
genon=True                                         # add generation shapefile
transon=True                                       # add transmission shapefile for the clusters
labels=False                                        # add labels for the clusters

#
# resource_assessment_exe.py
#
ncols_boxplot = 4                                   # columns required for the subplot (boxplot)
solar_density = 32                                  # installed capacity per squared kilometer
onwind_density = 5.33
offwind_density = 7.44
solar_subplot = [3,4]                               # define number of rows and cols for subplots
onshore_subplot = [3,4]
offshore_subplot = [3,4]
solar_kna = []                                      # remove this clusters from the analysis
onshore_kna = ['ONW-05','ONW-10','ONW-12','ONW-15']
offshore_kna = []

# paths
era_src=os.path.join(script_dir,'../../data/ERA5/')
shapeon=os.path.join(script_dir,'../../support/shapefile/qgis/onshore_only.shp')                # shapefile onshore borders (for solar) & onshore wind borders for Nth percentile
shapeonoff=os.path.join(script_dir,'../../support/shapefile/qgis/onshore_offshore_final.shp')   # shapefile onshore and offshore borders (for wind)
solar_geocoord=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_geocoord.csv')
solar_geoindex=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_geo_keep_index.csv')
solar_space=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_space.csv')
solar_time=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_time.csv')
solar_avrg=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_avrg.csv')
solar_sigma=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_sigma.csv')
solar_svd_fig1=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_plot.jpg')
solar_svd_pdf1=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_plot.pdf')
solar_svd_fig2=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_plot_labels.jpg')
solar_svd_pdf2=os.path.join(script_dir,'../../data/CLS/Solar/SVD/svd_plot_labels.pdf')
solar_geoindex_onshore=os.path.join(script_dir,'../../data/CLS/Solar/SVD/geoindex_onshore.csv')
solar_geoindex_offshore=os.path.join(script_dir,'../../data/CLS/Solar/SVD/geoindex_offshore.csv')
solar_db_out=os.path.join(script_dir,'../../data/OUT/Solar/solar_db_out.csv')
solar_db_out_label=os.path.join(script_dir,'../../data/OUT/Solar/solar_db_out_label.csv')
solar_cluster_fig1=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster1.jpg')
solar_cluster_pdf1=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster1.pdf')
solar_cluster_fig2=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster2.jpg')
solar_cluster_pdf2=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster2.pdf')
solar_cluster_fig3a=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster3a.jpg')
solar_cluster_fig3b=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster3b.jpg')
solar_cluster_pdf3b=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster3b.pdf')
solar_cluster_fig4a=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster4a.jpg')
solar_cluster_fig4b=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster4b.jpg')
solar_cluster_pdf4b=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster4b.pdf')
solar_cluster_fig4_label=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig4_label.csv')
solar_cluster_fig5=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster5.jpg')
solar_cluster_pdf5=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster5.pdf')
solar_cluster_fig6=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster6.jpg')
solar_cluster_pdf6=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/fig_cluster6.pdf')
solar_cluster_roc=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/roc.csv')
solar_singlek=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/pv_singlek.csv')
solar_subgrid=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/pv_subgrid.csv')
solar_assets_exop=os.path.join(script_dir,'../../data/ACT/Solar/Asset/pv_assets.csv')                                       # active assets
solar_subgrid_cf=os.path.join(script_dir,'../../data/CLS/Solar/Simulated/solar_tech_1_cf.csv')                              # simulated outputs cf (to provide path)
solar_subgrid_resass=os.path.join(script_dir,'../../data/CLS/Solar/Simulated/solar_tech_1_rass.csv')                        # resource assessment
solar_subgrid_sim=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/pv_subgrid_sim.csv')
solar_table_res=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/table_assessment.csv')
solar_dist_res=os.path.join(script_dir,'../../data/CLS/Solar/Cluster/table_distribution.csv')                               # cf distribution by cluster
gsa_src=os.path.join(script_dir,'../../data/GSA/GHI.tif')
wind_geocoord=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_geocoord.csv')
wind_geoindex=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_geo_keep_index.csv')
wind_space=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_space.csv')
wind_time=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_time.csv')
wind_avrg=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_avrg.csv')
wind_sigma=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_sigma.csv')
wind_svd_fig1=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_plot.jpg')
wind_svd_pdf1=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_plot.pdf')
wind_svd_fig2=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_plot_labels.jpg')
wind_svd_pdf2=os.path.join(script_dir,'../../data/CLS/Wind/SVD/svd_plot_labels.pdf')
wind_geoindex_onshore=os.path.join(script_dir,'../../data/CLS/Wind/SVD/geoindex_onshore.csv')
wind_geoindex_offshore=os.path.join(script_dir,'../../data/CLS/Wind/SVD/geoindex_offshore.csv')
wind_db_out=os.path.join(script_dir,'../../data/OUT/Wind/wind_db_out.csv')
wind_db_out_label=os.path.join(script_dir,'../../data/OUT/Wind/wind_db_out_label.csv')
wind_cluster_fig1=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster1.jpg')
wind_cluster_pdf1=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster1.pdf')
wind_cluster_fig2=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster2.jpg')
wind_cluster_pdf2=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster2.pdf')
wind_cluster_fig3a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster3a.jpg')
wind_cluster_fig3b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster3b.jpg')
wind_cluster_pdf3b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster3b.pdf')
wind_cluster_fig4a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster4a.jpg')
wind_cluster_fig4b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster4b.jpg')
wind_cluster_pdf4b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster4b.pdf')
wind_cluster_fig4_label=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig4_label.csv')
wind_cluster_fig5a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster5a.jpg')
wind_cluster_pdf5a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster5a.pdf')
wind_cluster_fig5b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster5b.jpg')
wind_cluster_pdf5b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster5b.pdf')
wind_cluster_fig6a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster6a.jpg')
wind_cluster_pdf6a=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster6a.pdf')
wind_cluster_fig6b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster6b.jpg')
wind_cluster_pdf6b=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/fig_cluster6b.pdf')
wind_cluster_roc=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/roc.csv')
wind_singlek=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wind_singlek.csv')
wind_subgrid=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wind_subgrid.csv')
wind_assets_exop=os.path.join(script_dir,'../../data/ACT/Wind/Asset/wind_assets.csv')                                       # active wind assets
wind_subgrid_on_cf=os.path.join(script_dir,'../../data/CLS/Wind/Simulated/wind_onshore.csv')                                # simulated outputs onshore cf (to provide path)
wind_subgrid_on_resass=os.path.join(script_dir,'../../data/CLS/Wind/Simulated/wind_onshore_rass.csv')                       # resource assessment onshore
wind_subgrid_off_cf=os.path.join(script_dir,'../../data/CLS/Wind/Simulated/wind_offshore.csv')                              # simulated outputs onshore cf (to provide path)
wind_subgrid_off_resass=os.path.join(script_dir,'../../data/CLS/Wind/Simulated/wind_offshore_rass.csv')                     # resource assessment onshore
wind_subgrid_sim=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/wind_subgrid_sim.csv')
wind_table_onshore_res=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/table_onshore_assessment.csv')                  # summary table for each cluster
wind_onshore_dist_res=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/table_onshore_distribution.csv')                 # cf distribution by cluster
wind_table_offshore_res=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/table_offshore_assessment.csv')
wind_offshore_dist_res=os.path.join(script_dir,'../../data/CLS/Wind/Cluster/table_offshore_distribution.csv')               # cf distribution by cluster
assess_summary_table_res=os.path.join(script_dir,'../../data/CLS/Assessment/table_summary.csv')                             # summary for all technologies
assess_corr_table=os.path.join(script_dir,'../../data/CLS/Assessment/table_correlation.csv')                                # correlation table      
assess_corr_chart=os.path.join(script_dir,'../../data/CLS/Assessment/correlation.jpg')                                      # correlation matrix
assess_corr_pdf=os.path.join(script_dir,'../../data/CLS/Assessment/correlation.pdf')                                        # correlation matrix pdf
gwa_src=os.path.join(script_dir,'../../data/GWA/BRA_wind-speed_100m.tif')
shp_regional_arg1=os.path.join(script_dir,'../../support/shapefile/regioes_2010/regioes_2010')                              # regions 
shp_regional_arg2='regioes_2010'                                                              
shp_ex_arg1=os.path.join(script_dir,'../../support/transmission/Linhas_de_Transmissão_-_Base_Existente')                    # existing transmission
shp_ex_arg2=os.path.join(script_dir,'Linhas_de_Transmissão_-_Base_Existente')
shp_plan_arg1=os.path.join(script_dir,'../../support/transmission/Linhas_de_Transmissão_-_Expansão_Planejada')              # planned transmission
shp_plan_arg2=os.path.join(script_dir,'Linhas_de_Transmissão_-_Expansão_Planejada')
shp_wind_arg1=os.path.join(script_dir,'../../support/generation/wind/wind')                                                 # wind assets
shp_wind_arg2='wind'
shp_solar_arg1=os.path.join(script_dir,'../../support/generation/solar/solar')                                              # solar assets
shp_solar_arg2='solar'
shp_wdpa_poly=os.path.join(script_dir,'../../support/shapefile/envsoc/envsoc_poly.shp')                                     # WDPA areas
shp_wdpa_pnt=os.path.join(script_dir,'../../support/shapefile/envsoc/envsoc_point.shp')                                     # WDPA points
gpw_tiff=os.path.join(script_dir,'../../support/GPW/gpw-v4/gpw_v4_population_density_rev11_2020_30_sec.tif')                # population density
globcover=os.path.join(script_dir,'../../support/globcover/Globcover2009_V2.3_Global_/GLOBCOVER_L4_200901_200912_V2.3.tif') # global land use
gebco=os.path.join(script_dir,'../../support/GEBCO/gebco_2023_n10.0_s-38.0_w-53.0_e-30.0.tif')                              # bathymetric data
srtm=os.path.join(script_dir,'../../support/srtm/full_srtm3_geotiff.tif')                                                   # shuttle radar topography mission

