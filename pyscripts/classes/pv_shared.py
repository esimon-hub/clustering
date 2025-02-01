"""
Solar PV Shared Functions

Author: Emanuel Simon
Date: 10/03/2023

Description:
This class includes several functions that are used to perform simulations and analysis

Usage:
Import this class into the pv_simulation.py

"""
#   Importing libraries and settings
import xarray as xr
import pandas as pd
import numpy as np
import pvlib
from datetime import datetime
from pvlib.location import Location
from pvlib.tools import cosd
from pvlib.pvsystem import PVSystem
from pvlib.modelchain import ModelChain
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS
from multiprocessing import Pool

# create the last entry for the radiation dataframe (23:30)
def complement_halfhour(timeseries):
    # replicate the last row of the dataframe (-24 means 00:00 of the previous day)
    last_row = timeseries.iloc[-24]
    last_timestamp = timeseries.index[-1] + pd.Timedelta(minutes=60)
    timeseries.loc[last_timestamp] = last_row
    # eliminate the first row of the dataframe (as it refers to the previous year of the simulation)
    timeseries = timeseries[1:]
    return timeseries

# adjust temperature and wind speed data to mid-hour (13:00 *13:30* 14:00)
def interp_halfhour(timeseries):
    # resample to consider 30T
    ts_resampled = timeseries.resample('30T').asfreq()
    # interpolate to fill in the NaN values
    ts_half_hour = ts_resampled.interpolate(method='linear')
    # replicate the last row of the dataframe
    last_row = ts_half_hour.iloc[-1]
    last_timestamp = ts_half_hour.index[-1] + pd.Timedelta(minutes=30)
    ts_half_hour.loc[last_timestamp] = last_row
    # get every other row and eliminate the instantaneous values at the full hour
    ts_half_hour = ts_half_hour[1::2]
    return ts_half_hour

# method to calculate PV output using PVLib (inclination equal to latitude or single-axis)
def pv_sim_multi(radiation, temperature, wind_speed, lat, lon, tracking, ts_up):
    # set day from index
    dayfromindex = ts_up.dayofyear
    # calculates location with a reference to UTC
    location = Location(latitude=lat,longitude=lon,tz='UTC')
    # temperature model
    temperature_model_parameters = TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']
    # get solar position
    ephem_pv_installation = pvlib.solarposition.pyephem(ts_up,location.latitude,location.longitude,temperature=np.mean(temperature))
    zenith_ang = ephem_pv_installation['zenith']
    if tracking==0:
        slope=lat
        if lat>=0:
            aspect = 0
        elif lat < 0:
            aspect = 180
    elif tracking == 1:
        tracker_data = pvlib.tracking.singleaxis(ephem_pv_installation['apparent_zenith'], \
            ephem_pv_installation['azimuth'],
            axis_tilt=0,
            axis_azimuth=0,
            max_angle=90,
            backtrack=True,
            gcr=2.0/7.0)
        slope  = tracker_data['surface_tilt']
        aspect = tracker_data['surface_azimuth']        
    # calculate DNI and DHI
    dni_pre = pvlib.irradiance.disc(radiation,zenith_ang,dayfromindex)['dni']
    dhi_pre = radiation - dni_pre*cosd(zenith_ang)  # radiation means GHI
    # create a dataframe to store the weather variables
    weather = pd.DataFrame({'ghi': radiation, 
                        'dni': dni_pre, 
                        'dhi': dhi_pre, 
                        'temp_air': temperature, 
                        'wind_speed': wind_speed},
                    index=ts_up)
    # select modules and inverters from PVLIB
    sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
    cec_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')
    # updated repo https://github.com/NREL/SAM/tree/develop/deploy/libraries
    pv_module = sandia_modules['Canadian_Solar_CS6X_300M__2013_']
    pv_inverter = cec_inverters['Sungrow_Power_Supply_Co___Ltd___SG125HV__600V_']
    output_inversor=125178
    #area_1kWp = number_of_panels_1kWp * sandia_module['Area']
    # prepare the system based on the inputs above
    system = PVSystem(surface_tilt=slope, surface_azimuth=aspect, \
        module_parameters=pv_module,
        inverter_parameters=pv_inverter,
        temperature_model_parameters=temperature_model_parameters,
        modules_per_string=27,
        strings_per_inverter=18)
    mc = ModelChain(system, location)
    mc.run_model(weather)
    pv_output = (mc.results.ac * 1/output_inversor).fillna(0)
    return pv_output

# eliminate the first row because it refers to the previous day (ssrd is accumulated)
# duplicate the last entry to represent the last element of the array
def cf_ts_adj(ts_cf):
    #  eliminate the first row of the dataframe
    ts_cf = ts_cf.iloc[1:]
    #  duplicate the last row of the dataframe
    last_row = ts_cf.iloc[[-1]]
    last_row = last_row.shift(freq='1H')
    ts_cf=pd.concat([ts_cf,last_row])
    return ts_cf

# Transform output data into capacity factors
def capacity_factors(ts_power):
    cf_1 =  ts_power/ts_power.max().copy()
    cf = cf_1.where(ts_power > 0,0).copy()
    return cf

