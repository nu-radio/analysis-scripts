from tabnanny import check
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import h5py
import glob
import os
import sys
import pickle
from functools import partial
from multiprocessing import Pool
from NuRadioReco.utilities import units
from NuRadioMC.utilities.Veff import get_Veff_water_equivalent as gwe
from NuRadioMC.utilities import cross_sections

def get_list_of_advanced_options():
    options = ['stations_and_triggers', 'triggers_coincidences']
    return options

def check_requested_mode(mode):
    if mode not in get_list_of_advanced_options():
        raise NotImplementedError('Advanced mode {} is not supported.'.format(mode))

def get_partial_func(mode, hybrid_list, shallow_list, deep_trigger, shallow_trigger):

    check_requested_mode(mode)

    if mode=='stations_and_triggers':
        local_func = partial(tmp_by_stations_and_triggers, 
            hybrid_list=hybrid_list, 
            shallow_list=shallow_list,
            deep_trigger=deep_trigger,
            shallow_trigger=shallow_trigger,
        )
        return local_func
    elif mode=='triggers_coincidences':
        local_func = partial(tmp_by_trigger_coincidences, 
            hybrid_list=hybrid_list, 
            shallow_list=shallow_list,
            deep_trigger=deep_trigger,
            shallow_trigger=shallow_trigger,
        )
        return local_func

def get_list_and_labels(mode):

    check_requested_mode(mode)

    the_list = None
    labels_dict = None


    if mode=='stations_and_triggers':
        
        """
        for breaking out by veff by stations AND triggers
        hd = hybrid, deep component
        hs = hybrid, shallow component
        ss = surface, shallow component
        hd_hs = hybrid-deep and hybrid-surface coincidence
        hs_ss = hybrid-deep and surface-shallow coincidnece
        hs_ss = hybrid-shallow and surface-shallow coincidence
        hd_hs_ss = all component coincience
        """
        the_list = ['total', 'hd', 'hs', 'ss', 'hd_hs', 'hd_ss', 'hs_ss', 'hd_hs_ss']
        labels_dict = {
            'total': 'Total',
            'hd': 'HD: Hybrid-Deep',
            'hs': 'HS: Hybrid-Shallow',
            'ss': 'SS: Surface-Shallow',
            'hd_hs': 'HD-HS: Hybrid-Deep + Hybrid-Shallow',
            'hd_ss': 'HD-SS: Hybrid-Deep + Surface-Shallow',
            'hs_ss': 'HS-SS: Hybrid-Shallow + Surface-Shallow',
            'hd_hs_ss': 'HD-HS-SS: All Components'
        }

    elif mode=='triggers_coincidences':

        """
        for breaking out by veff by trigger (deep vs shallow) coincidences
        se = shallow-exclusive, seen in only one shallow station
        ss = shallow-exclusive, seen in multiple shallow stations
        de = deep-exclusive, seen in only one deep station
        dd = deep-exclusive, seen in multiple deep stations
        sd = shallow-deep coincidences
        """
        the_list = ['total', 'se', 'ss', 'de', 'dd', 'ds']
        labels_dict = {
            'total':    'Total',
            'se':       '1 Shallow Trigger (shallow-only)',
            'ss':       '>1 Shallow Trigger (shallow-only)',
            'de':       '1 Deep Trigger (deep-only)',
            'dd':       '>1 Deep Trigger (deep-only)',
            'ds':       'Deep-Shallow Trigger Coincidence',
        }
    
    return the_list, labels_dict

def get_markers():
    markers = ['o', 's', '^', 'v', 'd', '+', 'x', '>']
    return markers

detectors = [
        "baseline_array",
        "hex_hybrid_only_array",
        "hex_shallow_array",
        "hex_shallowheavy_array"
    ]
def get_detectors():
    return detectors

flavors = [
    "e",
    "mu",
    "tau"
]
def get_flavors():
    return flavors


def tmp(filename, hybrid_list, shallow_list, deep_trigger, shallow_trigger,
    by_deepshallow_or_hybridhshallow='deepshallow'):

    """Summary or Description of the Function

    Parameters:
    filename: the name of the merged hdf5 file to be analyzed
    hybrid_list: an array/list of all the hybrid stations in the array
    shallow_list: an array/list of all the shallow-nly stations in the array
    deep_trigger: the string of the deep trigger used
    shallow_trigger: the string of the shallow trigger used
    by_deepshallow_or_hybridshallow: whether to break out the effective volumes
        by "deep vs shallow" or by "hybrid vs shallow-only"

    Returns:
    dictionary: dictionary of the energy, czmin, cmax, etc, along with the weights

   """

    assert by_deepshallow_or_hybridhshallow in ['deepshallow', 'hybridshallow'], 'Mode not supported'

    how_to_count = None
    if by_deepshallow_or_hybridhshallow == 'deepshallow':
        # this means the shallow component of the hybrid station should
        # count towards the shallow count
        how_to_count = 2
    elif by_deepshallow_or_hybridhshallow == 'hybridshallow':
        # this means the shallow component of the hybrid station should
        # count towards the hybrid count
        how_to_count = 1


    print("Filename {}".format(filename))

    fin = h5py.File(filename, "r")
    summary = {}
    summary['n_events'] = fin.attrs['n_events']
    summary['volume'] = fin.attrs['volume']
    summary['czmin'] = np.cos(fin.attrs['thetamax'])
    summary['czmax'] = np.cos(fin.attrs['thetamin'])
    summary['thetamin'] = fin.attrs['thetamin']
    summary['thetamax'] = fin.attrs['thetamax']
    summary['tot_weight'] = 0 

    # our goal is to figure out
    # first, how many events triggered in either the deep or the surface array
    # second, how many events triggered in just the deep array, just the surface array, or both
    # so, we will make a dictionary  of information
    # for every event (key), we will store an array (values)
    # the array will contain a weight
    # if it triggers at all
    # if it triggers deep, triggers shallow, or triggers both
    event_information = {}

    # get group ids of events that triggered
    ugids, uindex = np.unique(np.array(fin['event_group_ids']), 
        return_index=True)
    weights = np.array(fin['weights'])[uindex]
    weights_dict = {A:B for A,B in zip(ugids, weights)}
    try:
        tnames = fin.attrs['trigger_names']
    except:
        return None

    # tnames should be:
    # 0 = 'LPDA_2of4_100Hz'
    # 1 = 'LPDA_2of4_10mHz'
    # 2 = 'PA_8channel_100Hz'
    # 3 = 'PA_4channel_100Hz'
    # 4 = 'PA_8channel_1mHz'
    # 5 = 'PA_4channel_1mHz'
    # we want 0 and 3

    tname_to_index = {}
    for i, key in enumerate(tnames):
        tname_to_index[key] = i
    
    if deep_trigger not in tname_to_index:
        raise KeyError(f'Deep trigger ({deep_trigger}) is not in the hdf5 file list ({tnames})')
    if shallow_trigger not in tname_to_index:
        raise KeyError(f'Shallow trigger ({shallow_trigger}) is not in the hdf5 file list ({tnames})')


    # first, the hybrid-only part
    for key in fin.keys():
        if(key.startswith("station")):
            hybrid_id = int(key.split('_')[1])
            if hybrid_id not in hybrid_list:
                continue
            s_hybrid = fin[key]
            if 'multiple_triggers_per_event' in s_hybrid:

                # this is the deep component of the hybrid array
                hybrid_triggers_deep = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[deep_trigger]]
                evs_triggers_deep_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_deep])
                
                for ev in evs_triggers_deep_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][1] = 1 # seen in SOME deep component
                    elif ev in event_information:
                        event_information[ev][1] +=1 # increment the number of deep components where this is seen

                # now the shallow component of the hybrid array
                hybrid_triggers_shallow = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_shallow])
                
                for ev in evs_triggers_shallow_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][how_to_count] = 1 # seen in SOME shallow component
                    elif ev in event_information:
                        event_information[ev][how_to_count] +=1 # increment the number of shallow or hybrid components where this is seen

                
    # second, the shallow-only part
    for key in fin.keys():
        if(key.startswith("station")):
            shallow_id = int(key.split('_')[1])
            if shallow_id not in shallow_list:
                continue
            s_shallow = fin[key]
            if 'multiple_triggers_per_event' in s_shallow:
                triggers_shallow = s_shallow['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow = np.unique(np.array(s_shallow['event_group_ids'])[triggers_shallow])
                
                for ev in evs_triggers_shallow:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][2] = 1 # seen in SOME shallow station
                    elif ev in event_information:
                        event_information[ev][2] +=1 # increment the number of hybrid stations where this is seen

    tot_weight = 0
    deeporhybrid_only_weight = 0
    shallow_only_weight = 0
    dual_weight = 0

    for ev in event_information:
        weight = event_information[ev][0]
        trig_deeporhybrid = event_information[ev][1]
        trig_shal = event_information[ev][2]

        num_trig_total = trig_deeporhybrid + trig_shal
        if num_trig_total > 0:
            tot_weight+=weight
        
        if trig_deeporhybrid and trig_shal:
            dual_weight+=weight
        elif trig_deeporhybrid and not trig_shal:
            deeporhybrid_only_weight += weight
        elif trig_shal and not trig_deeporhybrid:
            shallow_only_weight += weight
        
    summary['tot_weight'] = tot_weight
    summary['deeporhybrid_only_weight']  = deeporhybrid_only_weight
    summary['shallow_only_weight'] = shallow_only_weight
    summary['dual_weight'] = dual_weight

    return summary

def tmp_by_stations_and_triggers(filename, hybrid_list, shallow_list, deep_trigger, shallow_trigger):

    """Summary or Description of the Function

    Parameters:
    filename: the name of the merged hdf5 file to be analyzed
    hybrid_list: an array/list of all the hybrid stations in the array
    shallow_list: an array/list of all the shallow-nly stations in the array
    deep_trigger: the string of the deep trigger used
    shallow_trigger: the string of the shallow trigger used

    Returns:
    dictionary: dictionary of the energy, czmin, cmax, etc, along with the weights

   """

    print("Filename {}".format(filename))

    fin = h5py.File(filename, "r")
    summary = {}
    summary['n_events'] = fin.attrs['n_events']
    summary['volume'] = fin.attrs['volume']
    summary['czmin'] = np.cos(fin.attrs['thetamax'])
    summary['czmax'] = np.cos(fin.attrs['thetamin'])
    summary['thetamin'] = fin.attrs['thetamin']
    summary['thetamax'] = fin.attrs['thetamax']
    summary['tot_weight'] = 0 

    # this version tries to break out the different part of the array volume more aggressively
    # there are three detector components in this model
    # the hybrid-deep, hybrid-shallow, and shallow-standadlone
    # so, we now need three indicators
    # make a dictionary of information
    # for every event (key), store an array (values)
    # the array contains a weight
    # and if it triggers at all, if it triggers hybrid-deep, hybrid-shallow, or shallow-standalone
    event_information = {}

    # get group ids of events that triggered
    ugids, uindex = np.unique(np.array(fin['event_group_ids']), 
        return_index=True)
    weights = np.array(fin['weights'])[uindex]
    weights_dict = {A:B for A,B in zip(ugids, weights)}
    try:
        tnames = fin.attrs['trigger_names']
    except:
        return None

    # tnames should be:
    # 0 = 'LPDA_2of4_100Hz'
    # 1 = 'LPDA_2of4_10mHz'
    # 2 = 'PA_8channel_100Hz'
    # 3 = 'PA_4channel_100Hz'
    # 4 = 'PA_8channel_1mHz'
    # 5 = 'PA_4channel_1mHz'
    # we want 0 and 3

    tname_to_index = {}
    for i, key in enumerate(tnames):
        tname_to_index[key] = i
    
    if deep_trigger not in tname_to_index:
        raise KeyError(f'Deep trigger ({deep_trigger}) is not in the hdf5 file list ({tnames})')
    if shallow_trigger not in tname_to_index:
        raise KeyError(f'Shallow trigger ({shallow_trigger}) is not in the hdf5 file list ({tnames})')


    # first, the hybrid-only part
    for key in fin.keys():
        if(key.startswith("station")):
            hybrid_id = int(key.split('_')[1])
            if hybrid_id not in hybrid_list:
                continue
            s_hybrid = fin[key]
            if 'multiple_triggers_per_event' in s_hybrid:

                # this is the deep component of the hybrid array
                hybrid_triggers_deep = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[deep_trigger]]
                evs_triggers_deep_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_deep])
                
                for ev in evs_triggers_deep_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0, 0]
                        event_information[ev][1] = 1 # seen in SOME hybrid-deep
                    elif ev in event_information:
                        event_information[ev][1] +=1 # increment the number of hybrid-deep components where this is seen

                # now the shallow component of the hybrid array
                hybrid_triggers_shallow = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_shallow])
                
                for ev in evs_triggers_shallow_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0, 0]
                        event_information[ev][2] = 1 # seen in SOME hybrid-shallow component
                    elif ev in event_information:
                        event_information[ev][2] +=1 # increment the number of hybrid-shallow components where this is seen

                
    # second, the shallow-only part
    for key in fin.keys():
        if(key.startswith("station")):
            shallow_id = int(key.split('_')[1])
            if shallow_id not in shallow_list:
                continue
            s_shallow = fin[key]
            if 'multiple_triggers_per_event' in s_shallow:
                triggers_shallow = s_shallow['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow = np.unique(np.array(s_shallow['event_group_ids'])[triggers_shallow])
                
                for ev in evs_triggers_shallow:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0, 0]
                        event_information[ev][3] = 1 # seen in SOME shallow-only station
                    elif ev in event_information:
                        event_information[ev][3] +=1 # increment the number of shallow-only stations where this is seen

    tot_weight = 0
    hd = 0 # hybrid-deep only #
    hs = 0 # hybrid-shallow only #
    ss = 0 # shallow-only #
    hd_hs = 0 # hybrid-deep + hybrid-shallow
    hd_ss = 0 # hybrid-deep + shallow-only
    hs_ss = 0 # hybrid-shallow + shallow-only
    all_evs = 0 # hd + hs + ss (all at the same time) #

    for ev in event_information:
        weight = event_information[ev][0]
        trig_hd = event_information[ev][1]
        trig_hs = event_information[ev][2]
        trig_ss = event_information[ev][3]


        num_trig_total = trig_hd + trig_hs + trig_ss
        if num_trig_total > 0:
            tot_weight+=weight

        if trig_hd and not trig_hs and not trig_ss:
            # print("Condition is hd {}, hs {}, so {}".format(trig_hd, trig_hs, trig_so))
            # only hybrid deep
            hd+=weight

        if trig_hs and not trig_hd and not trig_ss:
            # only hybrid shallow
            hs+=weight

        if trig_ss and not trig_hd and not trig_hs:
            # only shallow only
            ss+=weight
        
        if (trig_hd and trig_hs) and not trig_ss:
            # hybrid deep + hybrid shallow
            hd_hs+=weight

        if (trig_hd and trig_ss) and not trig_hs:
            # hybrid deep + shallow only
            hd_ss+=weight

        if (trig_hs and trig_ss) and not trig_hd:
            # hybrid shallow + shallow only
            hs_ss+=weight

        if trig_hd and trig_hs and trig_ss:
            # all three components
            all_evs+=weight

    summary['total_weight'] = tot_weight
    summary['hd_weight']  = hd
    summary['hs_weight'] = hs
    summary['ss_weight'] = ss
    summary['hd_hs_weight']  = hd_hs
    summary['hd_ss_weight'] = hd_ss
    summary['hs_ss_weight'] =  hs_ss
    summary['hd_hs_ss_weight'] = all_evs

    return summary


def tmp_by_trigger_coincidences(filename, hybrid_list, shallow_list, deep_trigger, shallow_trigger):

    """Summary or Description of the Function

    Parameters:
    filename: the name of the merged hdf5 file to be analyzed
    hybrid_list: an array/list of all the hybrid stations in the array
    shallow_list: an array/list of all the shallow-nly stations in the array
    deep_trigger: the string of the deep trigger used
    shallow_trigger: the string of the shallow trigger used

    Returns:
    dictionary: dictionary of the energy, czmin, cmax, etc, along with the weights

   """

    print("Filename {}".format(filename))

    fin = h5py.File(filename, "r")
    summary = {}
    summary['n_events'] = fin.attrs['n_events']
    summary['volume'] = fin.attrs['volume']
    summary['czmin'] = np.cos(fin.attrs['thetamax'])
    summary['czmax'] = np.cos(fin.attrs['thetamin'])
    summary['thetamin'] = fin.attrs['thetamin']
    summary['thetamax'] = fin.attrs['thetamax']
    summary['tot_weight'] = 0 

    # this version tries to break out the veff by triggers
    # care about
    # events seen only in 1 shallow trigger
    # events seen in >1 shallow trigger
    # events seen only in 1 deep trigger
    # events seen only >1 deep trigger
    # events een in both a shallow and deep trigger
    # there are three detector components in this model
    # only need two indicators: # deep, # shallow
    # make a dictionary of information
    # for every event (key), store an array (values)
    # the array contains a weight
    # if it triggers deep, and if it triggers shallow
    event_information = {}

    # get group ids of events that triggered
    ugids, uindex = np.unique(np.array(fin['event_group_ids']), 
        return_index=True)
    weights = np.array(fin['weights'])[uindex]
    weights_dict = {A:B for A,B in zip(ugids, weights)}
    try:
        tnames = fin.attrs['trigger_names']
    except:
        return None

    tname_to_index = {}
    for i, key in enumerate(tnames):
        tname_to_index[key] = i
    
    if deep_trigger not in tname_to_index:
        raise KeyError(f'Deep trigger ({deep_trigger}) is not in the hdf5 file list ({tnames})')
    if shallow_trigger not in tname_to_index:
        raise KeyError(f'Shallow trigger ({shallow_trigger}) is not in the hdf5 file list ({tnames})')

    # first, the hybrid-only part
    for key in fin.keys():
        if(key.startswith("station")):
            hybrid_id = int(key.split('_')[1])
            if hybrid_id not in hybrid_list:
                continue
            s_hybrid = fin[key]
            if 'multiple_triggers_per_event' in s_hybrid:

                # this is the deep component of the hybrid array
                hybrid_triggers_deep = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[deep_trigger]]
                evs_triggers_deep_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_deep])
                
                for ev in evs_triggers_deep_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][1] = 1 # seen in SOME deep component
                    elif ev in event_information:
                        event_information[ev][1] +=1 # increment the number of deep components where this is seen

                # now the shallow component of the hybrid array
                hybrid_triggers_shallow = s_hybrid['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow_hybrid = np.unique(np.array(s_hybrid['event_group_ids'])[hybrid_triggers_shallow])
                
                for ev in evs_triggers_shallow_hybrid:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][2] = 1 # seen in SOME shallow component
                    elif ev in event_information:
                        event_information[ev][2] +=1 # increment the number of shallow or hybrid components where this is seen

                
    # second, the shallow-only part
    for key in fin.keys():
        if(key.startswith("station")):
            shallow_id = int(key.split('_')[1])
            if shallow_id not in shallow_list:
                continue
            s_shallow = fin[key]
            if 'multiple_triggers_per_event' in s_shallow:
                triggers_shallow = s_shallow['multiple_triggers_per_event'][:, tname_to_index[shallow_trigger]]
                evs_triggers_shallow = np.unique(np.array(s_shallow['event_group_ids'])[triggers_shallow])
                
                for ev in evs_triggers_shallow:
                    if ev not in event_information:
                        event_information[ev] = [weights_dict[ev], 0, 0]
                        event_information[ev][2] = 1 # seen in SOME shallow station
                    elif ev in event_information:
                        event_information[ev][2] +=1 # increment the number of hybrid stations where this is seen

    tot_weight = 0
    se = 0 # only one shallow station
    ss = 0 # >1 shallow station
    de = 0 # only 1 deep station
    dd = 0 # >1 deep station
    ds = 0 # seen in both shallow and deep

    for ev in event_information:
        weight = event_information[ev][0]
        trig_deep = event_information[ev][1]
        trig_shallow = event_information[ev][2]

        num_trig_total = trig_deep + trig_shallow
        if num_trig_total > 0:
            tot_weight+=weight

        if trig_deep and trig_shallow:
            # trig both a deep and shallow trigger
            ds += weight
        
        if trig_deep == 0:
            # no deep trigger at all
            if trig_shallow == 1:
                # seen in only one shallow station
                se += weight
            elif trig_shallow > 1:
                # seen in multiple shallow stations
                ss += weight
        
        if trig_shallow == 0:
            # no shallow trigger at all
            if trig_deep == 1:
                # seen in only 1 deep station
                de += weight
            elif trig_deep > 1:
                # seen in multiple deep stations
                dd += weight

    summary['total_weight'] = tot_weight
    summary['se_weight']  = se
    summary['ss_weight'] = ss
    summary['de_weight'] = de
    summary['dd_weight']  = dd
    summary['ds_weight'] = ds

    return summary

