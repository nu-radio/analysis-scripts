import numpy as np


def efficiency_sigmoid(x, eff_low, eff_high, loge_turn, loge_halfmax):
    # sigmoid function in logE for efficiency between 0 and 1
    logx = np.log10(x)
    # choose factors conveniently
    # loge_halfmax should correspond to units in logE from turnover, where 0.25/0.75 of max are reached
    # = number of orders of magnitude in x between 0.25..0.75*(max-min) range
    b = np.log(3) / loge_halfmax

    eff = ((eff_low - eff_high) / (1 + (np.exp(b * (logx - loge_turn))))) + eff_high
    # do not allow below 0
    eff = np.maximum(0, eff)
    return  eff


def bound_efficiency_sigmoid(x, eff_low, eff_high, loge_turn, loge_halfmax):
    # sigmoid function in logE for efficiency between 0 and 1
    # hard limits between 0 and 1
    eff = efficiency_sigmoid(x, eff_low, eff_high, loge_turn, loge_halfmax)
    # limit to range between 0 and 1
    eff = np.maximum(0, eff)
    eff = np.minimum(1, eff)
    return  eff


def rnog_analysis_efficiency(E, minval, maxval, log_turnon_gev, log_turnon_width):
    # any3_gtr3
    # [-0.19848465  0.92543898  7.42347294  1.2133977 ]
    # 3phased_2support_gtr3
    # [-0.06101333  0.89062991  8.50399113  0.93591249]
    # 3power_gtr3
    # [-0.16480194  0.76853897  8.46903659  1.03517252]
    # 3power_2support_gtr3
    # [-0.0923469   0.73836631  8.72327879  0.85703575]
    # return bound_efficiency_sigmoid(E, -0.16480194,  0.76853897,  8.46903659,  1.03517252)
    return bound_efficiency_sigmoid(E, minval, maxval, log_turnon_gev, log_turnon_width)


def rnog_angular_reco_efficiency(logE):
    from scipy.interpolate import interp1d
    E = [16, 17, 18, 19, 20]
    e = [0, 0.03, 0.15, 0.38, 0.7]
    return interp1d(E, e, bounds_error=False, fill_value=(0, 1))(logE)


def PA_efficiency(logE):
    # Kaeli's analysis
    from scipy.interpolate import interp1d
    E = [16, 17, 18, 19]
    #e = [0.15, 0.25, 0.55, 0.65]
    e = [0.576 * 0.85, 0.688 * 0.85, 0.86 * 0.85, 0.94 * 0.85] 

    return interp1d(E, e, bounds_error=False, fill_value="extrapolate")(logE)

def PA_efficiency_pess(logE):
    # Kaeli's analysis
    from scipy.interpolate import interp1d
    E = [16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5]
    #e = np.array([0.55, 0.66, 0.672, 0.765, 0.851, 0.908, 0.939, 0.946]) * 0.85

    e = np.array([0.527, 0.646, 0.654, 0.753, 0.84, 0.904, 0.933, 0.951]) * 0.85
    #e = e / e

    return interp1d(E, e, bounds_error=False, fill_value="extrapolate")(logE)

def PA_efficiency_opp(logE):
    # Kaeli's analysis
    from scipy.interpolate import interp1d
    E = [16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5]
    e = np.array([0.624, 0.723, 0.722, 0.793, 0.868, 0.916, 0.941, 0.936])

    return interp1d(E, e, bounds_error=False, fill_value="extrapolate")(logE)


def Hpol3sigma(logE):
    from scipy.interpolate import interp1d
    E = [17, 18.6]
    e = [0.2, 0.27]
    return interp1d(E, e, bounds_error=False, fill_value="extrapolate")(logE)
