import numpy as np
import h5py
from six import iteritems

def read_input_hdf5_file(filename):
    fin = h5py.File(filename, 'r')
    _fin = {}
    _fin_stations = {}
    _fin_attrs = {}
    for key, value in iteritems(fin):
        if isinstance(value, h5py._hl.group.Group):
            _fin_stations[key] = {}
            for key2, value2 in iteritems(value):
                _fin_stations[key][key2] = np.array(value2)
        _fin[key] = np.array(value)
    for key, value in iteritems(fin.attrs):
        _fin_attrs[key] = value
    fin.close()
    return _fin, _fin_stations, _fin_attrs

def get_coszenbins():
    return np.linspace(-1,1,21)

def get_logEs():
    # return np.arange(18.5, 20.1, 0.5)
    return np.arange(19.0, 19.1, 0.5)

def get_number_of_parts_and_events(flavor, logE, czmin):
    num_parts=1
    num_events=1

    if flavor=='e':
        # the electron flavor is always fast, and can always take 1000 events
        if czmin < -0.5:
            num_events = int(10000)
            num_parts = int(10)
        else:
            num_events = int(5000)
            num_parts = int(20)
    elif flavor=='mu' or flavor=='tau':
        # otherwise, be clever

        if czmin < -0.3:
            num_events = int(10000)
            num_parts = int(10)
        
        elif (czmin < 0.4 and czmin >=-0.3 ):
            if logE<19:
                num_events = int(250)
                num_parts = int(400)
            elif (logE<19.5 and logE>=19):
                num_events = int(100)
                num_parts = int(1000)
            elif (logE<20 and logE>=19.5):
                num_events = int(50)
                num_parts = int(2000)
            elif (logE>=20):
                num_events = int(25)
                num_parts = int(4000)
        
        elif (czmin >= 0.4 ):
            if logE<19:
                num_events = int(1000)
                num_parts = int(100)
            elif (logE<19.5 and logE>=19):
                num_events = int(500)
                num_parts = int(200)
            elif (logE<20 and logE>=19.5):
                num_events = int(250)
                num_parts = int(400)
            elif (logE>=20):
                num_events = int(100)
                num_parts = int(1000)

    return num_parts, num_events

def get_number_of_parts_and_events_grid(flavor, logE, czmin):
    # this is to get number of parts (files) and events per file
    # for the grid, where Brian ran the very highest energies
    # these are tuned to finish in roughly 2 to 2.5 hours
    # assuming a relatively "moderate" array, like the baseline or hex hybrid
    num_parts = 1
    num_events = 1

    if logE>=19:
        # right now, Brian is only doing E==19 (other energies to come later)
        if flavor=='e':
            if czmin < -0.3:
                num_parts = int(10)
                num_events = int(5000)
            else:
                num_parts = int(250)
                num_events = int(200)
        if flavor=='mu':
            if czmin < -0.3:
                num_parts = int(10)
                num_events = int(5000)
            elif (czmin < 0.5 and czmin >=-0.3):
                num_parts = int(1000)
                num_events = int(50)
            elif (czmin >= 0.5):
                num_parts = int(500)
                num_events = int(100)
        if flavor=='tau':
            if czmin < -0.3:
                num_parts = int(10)
                num_events = int(5000)
            elif (czmin < 0.5 and czmin >=-0.3):
                num_parts = int(1000)
                num_events = int(50)
            elif (czmin >= 0.5):
                num_parts = int(250)
                num_events = int(200)
    
    return num_parts, num_events  

def get_file_pattern(flavor, logE, cosz_bin):
    coszenbins = get_coszenbins()
    czen1 = coszenbins[cosz_bin]
    czen2 = coszenbins[cosz_bin + 1]
    pattern = f"{flavor}_{logE:.2f}eV_{czen1:.1f}_{czen2:.1f}"
    return pattern

def get_gen_filename(flavor, logE, cosz_bin, part):
    pattern = get_file_pattern(flavor, logE, cosz_bin)
    out_filename = "in_" + f"{pattern}" + f".part{part:06}" + ".hdf5"
    return out_filename

def full_zmin_below_ice(logE, cz0, ice_depth=-2.7, default_dzmin=-0.6):
    # maximum additional z extension needed in km
    # calculated using 95% quantile of tau length for a given energy times the cz0 (or minimum cosz for which PREM transmisssion prob is > 1e-5)
    # values for ice have been divided by 2.5 (for rough / conservative estimate on rock)

    dzmin = {
        15.0: {-1.0: -0.1, -0.9: -0.1, -0.8: -0.0, -0.7: -0.0, -0.6: -0.0, -0.5: -0.0, -0.4: -0.0, -0.3: -0.0, -0.2: -0.0, -0.1: -0.0},
        15.5: {-1.0: -0.2, -0.9: -0.2, -0.8: -0.1, -0.7: -0.1, -0.6: -0.1, -0.5: -0.1, -0.4: -0.1, -0.3: -0.1, -0.2: -0.0, -0.1: -0.0},
        16.0: {-1.0: -0.5, -0.9: -0.5, -0.8: -0.5, -0.7: -0.4, -0.6: -0.4, -0.5: -0.3, -0.4: -0.2, -0.3: -0.2, -0.2: -0.1, -0.1: -0.1},
        16.5: {-1.0: -1.4, -0.9: -1.4, -0.8: -1.4, -0.7: -1.2, -0.6: -1.1, -0.5: -0.9, -0.4: -0.7, -0.3: -0.5, -0.2: -0.4, -0.1: -0.2},
        17.0: {-1.0: -2.6, -0.9: -2.6, -0.8: -2.6, -0.7: -2.6, -0.6: -2.6, -0.5: -2.3, -0.4: -1.8, -0.3: -1.4, -0.2: -0.9, -0.1: -0.5},
        17.5: {-1.0: -4.2, -0.9: -4.2, -0.8: -4.2, -0.7: -4.2, -0.6: -4.2, -0.5: -4.2, -0.4: -4.0, -0.3: -3.0, -0.2: -2.0, -0.1: -1.0},
        18.0: {-1.0: -5.4, -0.9: -5.4, -0.8: -5.4, -0.7: -5.4, -0.6: -5.4, -0.5: -5.4, -0.4: -5.4, -0.3: -5.2, -0.2: -3.5, -0.1: -1.7},
        18.5: {-1.0: -5.5, -0.9: -5.5, -0.8: -5.5, -0.7: -5.5, -0.6: -5.5, -0.5: -5.5, -0.4: -5.5, -0.3: -5.5, -0.2: -5.0, -0.1: -2.5},
        19.0: {-1.0: -4.8, -0.9: -4.8, -0.8: -4.8, -0.7: -4.8, -0.6: -4.8, -0.5: -4.8, -0.4: -4.8, -0.3: -4.8, -0.2: -4.8, -0.1: -3.0},
        19.5: {-1.0: -3.9, -0.9: -3.9, -0.8: -3.9, -0.7: -3.9, -0.6: -3.9, -0.5: -3.9, -0.4: -3.9, -0.3: -3.9, -0.2: -3.9, -0.1: -3.4},
        20.0: {-1.0: -3.4, -0.9: -3.4, -0.8: -3.4, -0.7: -3.4, -0.6: -3.4, -0.5: -3.4, -0.4: -3.4, -0.3: -3.4, -0.2: -3.4, -0.1: -3.4}}
    return ice_depth + min(dzmin[np.round(logE, 1)][np.round(cz0, 1)], default_dzmin)
