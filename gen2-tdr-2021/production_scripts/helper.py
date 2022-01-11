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
