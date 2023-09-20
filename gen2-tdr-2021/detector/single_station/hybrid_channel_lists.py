
# List of channel types and detector configurations for the Gen2_hybrid.json detector
import numpy as np
ch_ids = {}
ch_ids['lpda_ids_all'] = np.array([0, 1, 2, 3])
ch_ids['hpol_ids_all'] = np.array([12, 13, 14, 15, 18, 19, 20, 21, 24, 37, 40, 43, 46])
ch_ids['vpol_ids_all'] = np.array([4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 41, 42, 44, 45])
ch_ids['lpda_ids_baseline'] = np.array([0, 1, 2, 3])
ch_ids['vpol_ids_baseline'] = np.array([8, 9,10,11, 16, 17, 22, 23, 26, 27, 29, 31])
ch_ids['hpol_ids_baseline'] = np.array([12, 13, 18, 24])
ch_ids['channel_ids_power_string'] = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 25, 26, 27, 28, 29, 30, 31, 32])
ch_ids['PA_4ch'] = [8, 9, 10, 11]
ch_ids['vpol_bottom_baseline'] = [8, 16, 22]

ch_ids['vpol_ids_baseline_power'] = np.intersect1d(ch_ids['channel_ids_power_string'], ch_ids['vpol_ids_baseline'])
ch_ids['vpol_ids_baseline_support'] = np.setdiff1d(ch_ids['vpol_ids_baseline'], ch_ids['channel_ids_power_string'])
ch_ids['vpol_ids_baseline_power_noPA'] = np.setdiff1d(ch_ids['vpol_ids_baseline_power'], ch_ids['PA_4ch'])