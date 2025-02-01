# A solar and wind clustering framework with downscaling and bias correction of reanalysis data using singular value decomposition

A Python-based framework for systematically identifying and evaluating spatial clusters of solar and wind resources using multi-decade reanalysis data.

## Table of Contents
- [Description](#description)
- [Methodology](#methodology)
  - [1. Spatial Clustering](#1-spatial-clustering)
  - [2. Statistical Downscaling and Bias Correction](#2-statistical-downscaling-and-bias-correction)
  - [3. Simulation and Resource Potential Assessment](#3-simulation-and-resource-potential-assessment)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Description
This repository provides a framework designed to identify and evaluate spatial clusters for variable renewable energy (VRE) resources—specifically solar PV, onshore wind, and offshore wind. By leveraging decades of hourly data, the methodology uncovers the most relevant hourly patterns for each resource type and efficiently simulates their potential for energy generation. These clusters are statistically downscaled and bias-corrected with the Global Solar Atlas and the Global Wind Atlas, incorporating a new mesh of subgrids defined by the user. Subsequently, the framework provides a detailed assessment of the resources available within each cluster, encompassing similar hourly generation profiles from subgrids suitable for VRE.

## Methodology

### 1. Spatial Clustering
The methodology begins with the **spatial clustering** of VRE resources, focusing on solar radiation and wind speed. To manage and analyze these variables, singular value decomposition (SVD) is employed to derive a low-rank approximation of the solar radiation and wind speed matrices. Following the dimensionality reduction, the k-means clustering algorithm is applied to these lower-dimensional data structures.

### 2. Statistical Downscaling and Bias Correction
To enhance resolution and reduce biases for the the resource assessment, the framework implements **statistical downscaling and bias correction methods** with the support of a mesh of subgrids created by the algorithm. It refines coarse-resolution reanalysis data to a finer resolution using the Global Solar Atlas and the Global Wind Atlas. This  process establishes a relationship between the coarse-scale inputs and local-scale conditions.

### 3. Simulation and Resource Potential Assessment
Following the downscaling and bias-correction through subgrids, the following steps are performed:
1. **Simulations** are run for **solar PV**, **onshore wind**, and **offshore wind** at the subgrid level after excluding areas not suitable for VRE development. 
2. A **multiprocessing strategy** is employed to distribute the computation of these time series across multiple cores, significantly reducing execution time.
3. Once the simulations are complete, a **resource potential assessment** quantifies the VRE resources available within each spatial cluster.

## Datasets
Datasets utilized in the framework are available for download directly from the sources listed below.

### 1. ERA5
ERA5 is the fifth-generation of atmospheric reanalysis from the Copernicus Climate Change Service (C3S) at the European Centre for Medium-Range Weather Forecasts (ECMWF). It provides hourly estimates of atmospheric, ocean-wave, and land-surface variables from 1940 onwards. The data is structured in a gridded format, conveying representative values for each cell, as opposed to representing a particular point-measurement in time. Grid-cells are provided in a regular 0.25-degree latitude-longitude grid (approximately 28 km by 28 km), derived from the native reduced-Gaussian grid N320. Reanalysis products are accessed via the Copernicus Data Store (CDS) using its Application Programming Interface.

This framework relies on six reanalysis datasets to simulate solar PV and wind power generation. Orthogonal components at 10 and 100 meters are used to calculate the absolute wind speeds at their respective elevations. Additionally, hourly totals of downward solar radiation (SSRD) are converted to power per square meter (W.m-2) by dividing by 3600. Temperatures from the original reanalysis product are converted to Celsius.
Source: https://cds.climate.copernicus.eu/

### 2. Global Solar Atlas
Global Solar Atlas (GSA) version 2.0 provides a dataset of solar resources, spanning land areas from 60º North to 55º South, except for Latin America where the coverage extends to 45º South [73]. It provides key solar parameters on a grid with a spatial resolution of 9 arcseconds. GSA is produced by numerical models using satellite imagery, atmospheric, and meteorological data to estimate factors such as cloud transmittance, atmospheric conditions, and terrain characteristics. Data spans from 1999 to 2018 across the Americas and is validated by over 220 weather stations and measurement campaigns globally. Global Horizontal Irradiance (GHI) are extracted  from a raster file, representing the long-term annual averages of daily totals in kWh/m2.
Source: https://globalsolaratlas.info/map

### 3. Global Wind Atlas
Global Wind Atlas (GWA) version 3.3 provides a dataset for all terrestrial locations and marine areas within 200 kilometers of coastlines, excluding Antarctica. It offers wind resource parameters a spatial resolution of 9 arcseconds. GWA results from advances in global reanalysis products, improved geographical granularity of topographic data, and the development of numerical wind atlas. The dataset relies on a downscaling methodology employing mesoscale and microscale models to estimate wind speeds at various heights. Atmospheric reanalysis data spanning from 2008 to 2017 are used as inputs to the mesoscale model. Validation efforts of GWA include 35 sites across six countries. Weibull parameters (A and k) are selected from two raster files for 10 and 100 meters of altitude.
Source: https://globalwindatlas.info/em

### 4. Geographical datasets
World Database on Protected Areas (WDPA) is a database encompassing marine and terrestrial regions subject to restrictions which could impede the development of solar PV, onshore and offshore wind projects. This dataset provides detailed polygons and points that outline the boundaries and locations of protected areas.
Source: https://www.protectedplanet.net

GlobCover is an automated global-scale land mapping largely used for climate modeling, including land cover and land use dynamics. It includes 22 land cover types such as croplands, forests, wetlands, water bodies, artificial surfaces, and permanent snow and ice. GlobCover is available with a resolution of approximately 10 arcseconds.
Source: https://due.esrin.esa.int/page_globcover.php

Gridded Population of the World version 4 (GPWv4) provides detailed global population density estimates based on national censuses and population registers for 2020. GPWv4 uses proportional allocation gridding algorithms across nearly 13.5 million administrative units with a geographical resolution of 30 arcseconds. Population density estimates are available from raster format by dividing the population counts by the land area for the grid cells.
Source: https://doi.org/10.7927/H49C6VHW

Suttle Radar Topography Mission (SRTM-3) dataset is derived from a global digital elevation model of the Earth. Elevation is provided for terrestrial areas covering 60º North to 56º South, representing nearly 80% of the world’s land surface. SRTM-3 includes elevation data in a gridded format with a resolution of 3 arcseconds available.
Source: https://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL3.003

GEBCO_2023 is a dataset that contains a global elevation model for ocean and land built upon SRTM15+ [81]. It is derived from on a global bathymetry and topography grid with resolution of 15 arcseconds with coverage between latitudes of 60º North and 50º South. SRTM15+ was originated from shipboard soundings and estimates of depth using satellite altimetry.
Source: https://www.gebco.net/data_and_products/gridded_bathymetry_data/

---
## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/<YourUserName>/<YourRepoName>.git
    cd <YourRepoName>
    ```

2. **Create and activate a virtual environment** (recommended):
    ```bash
    # For Windows
    python -m venv venv
    venv\Scripts\activate
    
    # For macOS/Linux
    python3 -m venv venv
    source venv/bin/activate
    ```

3. **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
---
## Data Paths and Directory Structure

### Overview
Due to the substantial size of the datasets used in this project, they are not hosted directly within the GitHub repository. It is recommended that users independently download the necessary data using the links provided previously. Following a consistent directory structure is crucial for ensuring that the scripts function correctly as they rely on specific paths to access input data and to save outputs.

### Recommended Directory Structure
To facilitate the correct functioning of the framework, a brief description of the folders is provided below.

- `data/`: Root directory for all data files.
  - `ERA5/`: Contains all reanalysis files (*.nc).
  - `ACT/`: Data specific to asset classes.
    - `Solar/`: Solar-related data.
      - `Asset/`: Inputs for simulation using existing solar PV assets.
      - `Simulated/`: Outputs from the hourly simulation for solar PV assets.
    - `Wind/`: Wind-related data.
      - `Asset/`: Inputs for simulation using existing wind assets.
      - `Simulated/`: Outputs from the hourly simulation for wind assets.
  - `CLS/`: Cluster analysis outputs and assessments.
    - `Assessment/`: Correlation outputs among all clusters (solar PV, onshore wind, and offshore wind).
    - `Solar/`: Solar PV clustering data.
      - `Cluster/`: Output data for solar PV clusters.
      - `Simulated/`: Hourly data outputs for solar PV clusters.
      - `SVD/`: Outputs from the singular value decomposition for solar.
    - `Wind/`: Wind clustering data.
      - `Cluster/`: Output data for wind clusters.
      - `Simulated/`: Hourly data outputs for wind clusters.
      - `SVD/`: Outputs from the singular value decomposition for wind.
  - `GSA/`: Global Solar Atlas data.
    - `GHI.tif`: GeoTIFF file for solar radiation.
  - `GWA/`: Global Wind Atlas data.
    - `BRA_combined-Weibull-A_10m.tif`: Weibull distribution parameter 'A' at 10m.
    - `BRA_combined-Weibull-A_100m.tif`: Weibull distribution parameter 'A' at 100m.
    - `BRA_combined-Weibull-k_10m.tif`: Weibull distribution parameter 'k' at 10m.
    - `BRA_combined-Weibull-k_100m.tif`: Weibull distribution parameter 'k' at 100m.
- `pyscripts/`: Root directory for the Python scripts.
  - `classes/`: Classes for the framework.
    - `clusters.py`: Methods required to perform the clustering techniques.
    - `pv_shared.pv`: Simulation engine for solar PV and time series treatment.
    - `pv_simulation.pv`: Methods defined to perform simulations either at the plant or subgrid level.
    - `resource_assessment.pv`: Methods implemented to evaluate solar and wind resources.
    - `svd.pv`: Perform SVD and determine the truncated version of the datasets.
    - `wind_shared.pv`: Shared methods to perform onshore and offshore wind simulations.
    - `wind_simulation.pv`: Methods defined to perform simulations either at the plant or subgrid level.    
  - `common/`: Definitions of the framework. 
    - `main_config.py`: Configure main parameters to perform SVD and clustering, including file location.
    - `pv_config.pv`: Configure parameters to perform solar PV simulation, including file location.
    - `wind_config.pv`: Configure parameters to perform wind simulation, including file location.	
  - `utils/`: Root directory with several scripts to support further analysis and map generation.
  - `clusters_exe`: Execution of K-Means clustering based on lower-dimensional data.
  - `resource_assessment_exe`: Perform resource assessment after simulations.
  - `simulation_pv_exe`: Perform solar PV simulation.
  - `simulation_wind_exe`: Perform onshore and offshore wind simulation.
  - `svd_exe`: Perform singular value decomposition.
- `support/`: Supporting geographical datasets.
  - `GEBCO/`: Bathymetric data.
    - `gebco_2023_n10.0_s-38.0_w-53.0_e-30.0.tif`: GEBCO dataset for water depth.
  - `srtm/`: Radar topography data.
    - `full_srtm3_geotiff.tif`: SRTM3 - GeoTIFF.
  - `globcover/`: Land cover data.
    - `GLOBCOVER_L4_200901_200912_V2.3.tif`: GeoTIFF with GlobCover data.
  - `GPW/`: Population data.
    - `gpw_v4_population_density_rev11_2020_30_sec.tif`: GeoTIFF with population density.
  - `shapefile/envsoc/`: Shapefiles for environmental and social data.
    - `envsoc_poly.shp`: WDPA shapefile for polygons.
    - `envsoc_point.shp`: WDPA shapefile for points.
	

### Preparation and Verification
**Before running any scripts, ensure that all required files are correctly placed within the designated directories.** Please check the configuration files described below where more information about the paths are provided.

---

## Requirements

Below are key Python packages required for running the scripts. See `requirements.txt` for the complete list with pinned versions.

- `pandas==1.5.3`
- `geopandas==0.12.2`
- `matplotlib==3.5.3`
- `numpy==1.23.5`
- `shapely==1.8.2`
- `xarray==2023.2.0`
- `rasterio==1.3.9`
- `scipy==1.10.1`
- `sklearn==0.0.post1`
- `pvlib==0.9.5`
- `psutil==5.9.8`
- `seaborn==0.12.2`
- `networkx==3.3`
- `cdsapi==0.6.1`

---

## Usage
Please edit the following configuration files to perform the tasks.

1. **Perform Singular Value Decomposition**  
   - Configure parameters for the SVD by editing parameters in the following file: **common/main_config.py**.
   - Perform singular value deccomposition for solar and wind resources
   - Example command:  
     ```bash
     python svd_exe.py [options]
     ```
   - [Options]
   - Below are the options available, each triggering specific functionalities within the script:
      - **`--preset`**: Preprocess the solar and wind datasets to reduce size (e.g.: excluding areas with relatively lower sources and geographical boundaries). This is a required step.
      - **`--run_svd`**: Run the Singular Value Decomposition module.
      - **`--plot`**: Plot the results of the decomposition (e.g.: singular values and relative errors)

2. **Clustering renewables with truncated SVD datasets**  
   - Configure parameters for clustering in the following file: **common/main_config.py**.
   - This script performs various clustering operations. It allows for specific functions to be executed via command-line arguments. 
   - Example command:  
     ```bash
     python clusters_exe.py [options]
     ```
   - [Options]
   - Below are the options available, each triggering specific functionalities within the script:
      - **`--calc_wcss`**: Calculate the Within Cluster Sum of Squares (WCSS) to evaluate the optimal number of clusters.
      - **`--ksubplots`**: Generate plots of clusters varying with the number `k` to visualize the effect of different cluster sizes.
      - **`--singlek_pre`**: Create a dataframe and plot a map for a predetermined `k` value.
      - **`--opassets`**: Plot a map that includes existing assets and transmission lines.
      - **`--postproc`**: Perform post-processing operations such as creating subgrids and removing exclusion zones.
      - **`--seltop`**: Select the top `N` resources for each subgrid, aiding in resource allocation and planning.
      - **`--simsubgrid_view`**: Visualize simulated subgrids on a map to assess their practical layout and feasibility.
      - **`--simgrid_view`**: Display a map with clusters and labels for easy identification and analysis.

3. **Run Solar PV simulations**  
   - Configure parameters for solar PV simulations in the following file: **common/pv_config.py**.
   - Perform simulations for solar PV at either the plant or subgrid level.  
   - Example command:  
     ```bash
     python simulation_pv_exe.py [options]
     ```
   - [Options]
   - Below are the options available, each triggering specific functionalities within the script:
      - **`--PVAsset`**: Run solar PV simulations for a list of existing assets. This does not use the cluster definition and is merely used for benchmark purposes. The option allows for simulation at the plant level with or without SD/BC.
      - **`--PVMeshAsset`**: Run solar PV simulations for a list of existing assets. This does not use the cluster definition and is merely used for benchmark purposes. The option allows for simulation at the subgrid level.
      - **`--PVCluster`**: Run solar PV simulation for each subgrid individually. This option provides the results for the resource assessment considering the clustering methodology.

3. **Run onshore (and offshore, if applicable) simulations**  
   - Configure parameters for onshore simulations in the following file: **common/wind_config.py**.
   - Perform simulations for wind at either the plant or subgrid level.  
   - Example command:  
     ```bash
     python simulation_wind_exe.py [options]
     ```
   - [Options]
   - Below are the options available, each triggering specific functionalities within the script:
      -  **`--WindAsset`**: Run wind simulations for a list of existing assets. This does not use the cluster definition and is merely used for benchmark purposes. The option allows for simulation at the plant level with or without SD/BC.
      -  **`--WindMeshAsset`**: Run wind simulations for a list of existing assets. This does not use the cluster definition and is merely used for benchmark purposes. The option allows for simulation at the subgrid level.
      -  **`--WindCluster`**: Run simulations for each subgrid individually. This option provides the results for the resource assessment considering the clustering methodology.

4. **Resource Potential Assessment**  
   - Configure parameters for onshore simulations in the following file: **common/main_config.py**.
   - After simulations, run the resource assessment script to compute energy potential metrics for each cluster.  
   - Example command:  
     ```bash
     python resource_assessment_exe.py
     ```
   - [Options]
   - Below are the options available, each triggering specific functionalities within the script:
      - **`--tables`**: Create table with geographical summary.
      - **`--assessment`**: Perform resource potential assessment for each cluster.
      - **`--correl`**: Create a correlation matrix with all cluster.

---

## Contributing

1. **Fork the Project**  
2. **Create a Feature Branch** (`git checkout -b feature/YourFeature`)
3. **Commit Your Changes** (`git commit -m 'Add YourFeature'`)
4. **Push to the Branch** (`git push origin feature/YourFeature`)
5. **Open a Pull Request** describing your feature in detail

We welcome contributions that enhance the methodology, add new functionalities, or improve documentation.

---

## License

MIT License

Copyright (c) 2025 esimon-hub

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
