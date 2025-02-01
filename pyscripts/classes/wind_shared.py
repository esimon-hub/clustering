"""
Wind Shared Functions

Author: Emanuel Simon
Date: 10/03/2023

Description:
This class includes several functions that are used to perform simulations and analysis

Usage:
Import this class into the wind_simulation.py

"""
#   Importing libraries and settings
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

# parameters for onshore wind using parametric approach by Ryberg
capFac = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 
          18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 
          36.0, 37.0, 38.0, 39.0,  40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 
          54.0, 55.0, 56.0, 57.0, 58.0,  59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 
          72.0, 73.0, 74.0, 75.0, 76.0, 77.0,  78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 
          90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0,  97.0, 98.0, 99.0, 100.0,])

constA = np.array([-0.8691217362419719, -0.5846306874157908, -0.4199402034594731, -0.39102684258488274, 
          -0.3894573187489953, -0.4058020013186425, -0.40432027949816374, -0.39236913666698187, -0.3819249132228733, 
          -0.36330046354548506, -0.3435178362137303, -0.3178321425403834, -0.29117508148616195, -0.269508217230584, 
          -0.2532127879218888, -0.23716397286478894, -0.22262634723157404, -0.20608305884473196, -0.18294597412929683, 
          -0.1595198473797326, -0.1370173698072664, -0.11673980125939647, -0.10109357095443004, -0.08657611801800887, 
          -0.07475367908452243, -0.0632268060157558, -0.05125428172949209, -0.039703953493727594, -0.028335883786178517, 
          -0.01558658032251571, -0.0023729354603777156, 0.012792316480288236, 0.034139993622620955, 0.054711366860818865, 
          0.07301254668013582, 0.0898846001854323, 0.10482533277884826, 0.1150798118926572, 0.11964163861256258, 
          0.1242663158085397, 0.1294060065593988, 0.1352503099574173, 0.14152170196334551, 0.14683148290674713, 
          0.15136974087405913, 0.15579167325105514, 0.16027451616989363, 0.16455629622664028, 0.16836020011407535, 
          0.17250900734008895, 0.1765350074146051, 0.18057485786938174, 0.18426444139957254, 0.1885531927827845, 
          0.19215673631235866, 0.19482283957199606, 0.19705093291250078, 0.19823204136166156, 0.1986947366716821, 
          0.19857859581037052, 0.19826130240180936, 0.19787106118220874, 0.19697825851443274, 0.19541633763382596, 
          0.19326578494880048, 0.19090088568221777, 0.18833523184930082, 0.1856944901502562, 0.1823018438171243, 
          0.17902157722026388, 0.17583645944964005, 0.17433875596971168, 0.17282643301966263, 0.17183139919421217, 
          0.1712745290384288, 0.17100886657189338, 0.17155504072621522, 0.17181169091183046, 0.1717475206313266, 
          0.17205953689104467, 0.17203759637770671, 0.17135373589368072, 0.17038596526121402, 0.1696531890637046, 
          0.16777914643443562, 0.16631979884256623, 0.1650077001337763, 0.1638990474023653, 0.16392813072956272, 
          0.16592132918666147, 0.16818170112634986, 0.16826496645049946, 0.1685499709652857, 0.17168741672630294, 
          0.1785778934680895, 0.18353861618665934, 0.17884824961686804, 0.17285863610537955, 0.15069868246914123, 
          0.13128579344491695, 0.11012455979004396,])

constB = np.array([0.30864428729945903, 0.298395231937226, 0.28855482910643393, 0.29670835655698713, 
         0.30694659174340955, 0.3186606646987472, 0.3262109237633102, 0.3309590104414416, 0.3353076855262877, 
         0.3377929657612922, 0.33963490414680153, 0.33994744034478136, 0.3397628567355643, 0.3400709179760609, 
         0.34104589967775417, 0.34182215842728025, 0.342677119769971, 0.3430452815466755, 0.34212097609298914, 
         0.3409467860008212, 0.3397912532745987, 0.33891784784556933, 0.33873460341449246, 0.3386213565615072, 
         0.33886827438947414, 0.3390797757839704, 0.3391432173371577, 0.33922232252559076, 0.33926121476773186, 
         0.33900163826510576, 0.33859625081614964, 0.33780110674821306, 0.33589696154420867, 0.33409700749274673, 
         0.33264489345560033, 0.33139670650839426, 0.3304333161479266, 0.33022289054218695, 0.33093232690690216, 
         0.3315942385260518, 0.33213539920570395, 0.3325279563119545, 0.3328180841186797, 0.33324903030112274, 
         0.3337880560119182, 0.3343176761761118, 0.33481440494506853, 0.33532147113139826, 0.3358949127330813, 
         0.3363941497376186, 0.3369001413399903, 0.3373854094451347, 0.33791086689131156, 0.33830763653262885, 
         0.3388061917885345, 0.33945446631863424, 0.3401692545183817, 0.3410558652808998, 0.3420604313985917, 
         0.3431598218970714, 0.34428530273814234, 0.34541702411566705, 0.34662664250084896, 0.34794309830105724, 
         0.3493553281647887, 0.3507948657524085, 0.3522646588773047, 0.3537436504328669, 0.35535160917143216, 
         0.3569385284006994, 0.35850750514375196, 0.35978319892006216, 0.3610620326969219, 0.3622481306704468, 
         0.3633573244159112, 0.3644204304211945, 0.36534953287098765, 0.3663298056163852, 0.36736857090133734, 
         0.36834900739477233, 0.36940125738628526, 0.370576568712213, 0.3718120905892753, 0.3730195804206644, 
         0.3744312805803984, 0.37580157824845994, 0.37720392901934174, 0.3786249969503181, 0.37989766263735397, 
         0.38086578501576873, 0.3818189748119938, 0.38319957281022093, 0.38462943240495706, 0.38567391597445144, 
         0.38618181662154005, 0.38720768445078335, 0.3902082649368858, 0.39378413492769676, 0.4008966514755601, 
         0.4087302008778277, 0.4212662545480378,])

# NREL IEA-15-240-RWT derived @ https://github.com/IEAWindTask37/IEA-15-240-RWT/blob/master/Documentation/IEA-15-240-RWT_tabular.xlsx
offshore_ws = np.array([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 
                        11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 
                        19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0])

offshore_output = np.array([0.00283308303523257, 0.0177530081392936, 0.0374416715958697, 0.0623630027020519, 0.092998591126361, 0.12985207435673, 
                            0.173436565753144, 0.224264919995628, 0.282813047483084, 0.34972365750845, 0.425739429296416, 0.511373003897932, 
                            0.605711074255507, 0.711292914042876, 0.822190477797176, 0.976635181139334, 0.999617322939717, 1.0, 0.999956530286412, 
                            0.999998491198691, 0.999984043990525, 0.999984740291669, 0.999992283899712, 0.999992931574012, 0.999990528158469, 
                            0.999990162360238, 0.99998969913173, 0.999989266994872, 0.999989023940737, 0.9999888046739, 0.999988606040099, 
                            0.999988457625655, 0.999988341909147, 0.999988116255081, 0.99998806504247, 0.999988425546058, 0.999987324389552, 
                            0.999986599978043, 0.999993648207722, 1.0, 1.0, 0.999996515354873, 0.999997166608538, 1.0, 1.0])

# Vestas V150-4.2 MW
onshore_ws = np.array([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 
                       14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5])

onshore_output = np.array([0.0185714285714286, 0.040952380952381, 0.0683333333333333,	0.101428571428571, 0.143095238095238, 0.193809523809524, 
                           0.25452380952381, 0.32547619047619, 0.408809523809524, 0.502380952380952, 0.606190476190476, 0.714761904761905, 
                           0.816190476190476, 0.898333333333333, 0.955238095238095, 0.983571428571429, 0.996666666666667, 0.999523809523809, 
                           1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])


# 
def syntpc(specp, cutout=25.00):
    # create the curve using the coeficients proposed by Ryberg
    specp=int(specp)
    windspeed = [0,]
    windspeed.extend(np.exp(constA + constB*np.log(specp)))
    ratedWS = np.arange(windspeed[-1] + 0.25, cutout + 0.25, 0.25)
    windspeed.extend(ratedWS)
    windspeed = np.array(windspeed)
    cfactors = [0,]
    cfactors.extend(capFac/100)
    cfactors.extend([1]*len(ratedWS))
    cfactors = np.array(cfactors)
    # create a dictionary with ws and cf
    pc = {
        'ws': windspeed,
        'cf': cfactors
    }
    return pc

# calculate onshore wind capacity factors (Ryberg)
def calc_on_cf(wind_series, spower):
    # get a dictionary with the wind speed and cf for a given specific power
    ws_cf=syntpc(specp=spower)
    iws=ws_cf['ws']
    icf=ws_cf['cf']
    calc_inter = interp1d(iws,icf)(wind_series) # interp1d is a class that returns an object and we give the second paranthesis as argument
    return pd.Series(calc_inter)
    
# calculate offshore wind capacity factors (NREL WPC without filtering)
def calc_of_cf(wind_series):
    interp_func=interp1d(offshore_ws, offshore_output, fill_value="extrapolate")
    cf_out=interp_func(wind_series)
    return pd.Series(cf_out)

# apply gaussian filtering in the wind power curve
def wpc_gaussian(vmin, vmax, ws, wpc):
    # set parameters for the filter
    res = 0.01
    a = 0.15
    b = 0.6

    # create a high-resolution dataset
    ws_full = np.arange(0,40 + res, res)        # full means adding data below cut-in and above cut-off
    output_full = np.zeros_like(ws_full)
    output_gaussian = np.zeros_like(ws_full)
    valid_range = (ws_full >= vmin) & (ws_full <= vmax)
    spline = CubicSpline(ws, wpc)
    output_full[valid_range] = spline(ws_full[valid_range])

    # zero-padding to perform convolution
    sigma_max = b + a*max(ws_full)
    padding_size = int(4*sigma_max/res)
    output_full_padded = np.pad(output_full, (padding_size, padding_size), mode='constant', constant_values=0)

    # manual implementation of the gaussian filter through convolution
    for i in range(len(ws_full)):
        # calculate sigma for a given wind speed
        sigma = b + a*ws_full[i]
        lb = int(-4*sigma/res)
        up = int(4*sigma/res)
        conv_pc = 0
        for j in range(lb,up+1):
            wj = j * res
            window_pc = output_full_padded[padding_size + i - j]
            gamma = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-(wj ** 2) / (2 * sigma ** 2))
            conv_pc+=window_pc*gamma*res
        output_gaussian[i] = conv_pc
    
    # return both power curves in a dictionary
    wpc = pd.DataFrame({
        'ws': ws_full,
        'wpc_full': output_full,
        'wpc_gaussian': output_gaussian 
    })
    return wpc

# create onshore wind power curve for clusters with gaussian filtering
def wpc_onshore():
    vmin = min(onshore_ws)
    vmax = max(onshore_ws)
    df_wpc = wpc_gaussian(vmin, vmax, onshore_ws, onshore_output)
    interp_wpc=interp1d(df_wpc['ws'], df_wpc['wpc_gaussian'], fill_value="extrapolate")
    return df_wpc, interp_wpc

# create offshore wind power curve for clusters with gaussian filtering
def wpc_offshore():
    vmin = min(offshore_ws)
    vmax = max(offshore_ws)
    df_wpc = wpc_gaussian(vmin, vmax, offshore_ws, offshore_output)
    interp_wpc=interp1d(df_wpc['ws'], df_wpc['wpc_gaussian'], fill_value="extrapolate")
    return df_wpc, interp_wpc

# perform simulation using the adjusted power curve for both onshore and offshore projects
def calc_cf(ts, interp):
    cf_out=interp(ts)
    return pd.Series(cf_out)