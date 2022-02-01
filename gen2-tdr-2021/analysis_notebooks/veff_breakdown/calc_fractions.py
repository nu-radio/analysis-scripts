import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import h5py
import glob
import os
import sys
import pickle
from multiprocessing import Pool, pool
from functools import partial
from NuRadioReco.utilities import units
from NuRadioMC.utilities.Veff import get_Veff_water_equivalent as gwe
import helper
import argparse

"""
This script caluculates the fractions of the effective volume in
various detector components. It has two modes.
A "deepshallow" mode, where the coincidences are calculated between
the deep and shallow components of the array.
And, a "hybridshallow" mode, where the coincidences are calculated
between the hybrid stations and surface stations.
So they are slightly different ways of slicing things up.

The output is a pkl file with the triggering stats and the veffs.

"""


parser = argparse.ArgumentParser()
parser.add_argument('--det', required=True, help='which detector index')
parser.add_argument('--mode', required=True, help='which mode: deepshallow or hybridshallow')
args = parser.parse_args()

the_det = int(args.det)
the_mode = args.mode

top = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3"
detectors = [
    "baseline_array",
    "hex_hybrid_only_array",
    "hex_shallow_array",
    "hex_shallowheavy_array"
]

detector = detectors[the_det]
config = "config_ARZ2020_noise"
detsim = "D01detector_sim"
path=os.path.join(top, detector, config, detsim)

deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'
mode = the_mode
hybrid_list = np.genfromtxt(f"../../detector/station_lists/stations_{detector}_hybrid.csv")
shallow_list = np.genfromtxt(f"../../detector/station_lists/stations_{detector}_shallow.csv")

flavors = [
    "e",
    "mu",
    "tau"
]

zen_bins = np.linspace(-1,1,21)
coszen_to_bin = {}
for iZ, czen1 in enumerate(zen_bins[:-1]):
    coszen_to_bin[f"{czen1:.1f}"] = iZ
logEs = np.arange(16, 20.1, 0.5)
n_cores = 10

local_func = partial(helper.tmp, 
    hybrid_list=hybrid_list, 
    shallow_list=shallow_list,
    deep_trigger=deep_trigger,
    shallow_trigger=shallow_trigger,
    by_deepshallow_or_hybridhshallow=mode
    )

for flavor in flavors:
    total_veff = {}
    for logE in logEs:
        total_veff[f"{logE:.1f}"] = {}
        
        combined_total_veff = 0
        combined_deeporhybrid_only_veff = 0
        combined_shallow_only_veff = 0
        combined_dual_veff = 0
        num_zen_bins = 20
        combined_total_veff_bins = np.zeros(20)
        combined_deeporhybrid_only_veff_bins = np.zeros(20)
        combined_shallow_only_veff_bins = np.zeros(20)
        combined_dual_veff_bins = np.zeros(20)
        n_thrown = 0
        trig_weight = 0
        
        print("Flavor {}, logE {}".format(flavor, logE))

        local_path = os.path.join(path, f"{flavor}/{flavor}_{logE:.2f}*.hdf5")
        filenames = glob.glob(local_path)
        print("Filenames {}".format(filenames))

        with Pool(n_cores) as p:
            pool_result = p.map(local_func, filenames)
        
        for result in pool_result:
            if result is None:
                continue
                
            czmin = str('{:.1f}'.format(result['czmin']))
            zen_bin = coszen_to_bin[czmin]
            combined_total_veff_bins[zen_bin] += result['tot_weight']/result['n_events'] * result['volume']
            combined_deeporhybrid_only_veff_bins[zen_bin] += result['deeporhybrid_only_weight']/result['n_events'] * result['volume']
            combined_shallow_only_veff_bins[zen_bin] += result['shallow_only_weight']/result['n_events'] * result['volume']
            combined_dual_veff_bins[zen_bin] += result['dual_weight']/result['n_events'] * result['volume']

            combined_total_veff += result['tot_weight']/result['n_events'] * result['volume']
            combined_deeporhybrid_only_veff += result['deeporhybrid_only_weight']/result['n_events'] * result['volume']
            combined_shallow_only_veff += result['shallow_only_weight']/result['n_events'] * result['volume']
            combined_dual_veff += result['dual_weight']/result['n_events'] * result['volume']
            n_thrown += result['n_events']
            trig_weight += result['tot_weight']

        total_veff[f"{logE:.1f}"]['total_veff'] = gwe(combined_total_veff/num_zen_bins)
        total_veff[f"{logE:.1f}"]['deeporhybrid_only_veff'] = gwe(combined_deeporhybrid_only_veff/num_zen_bins)
        total_veff[f"{logE:.1f}"]['shallow_only_veff'] = gwe(combined_shallow_only_veff/num_zen_bins)
        total_veff[f"{logE:.1f}"]['dual_veff'] = gwe(combined_dual_veff/num_zen_bins)
        total_veff[f"{logE:.1f}"]['n_thrown'] = n_thrown
        total_veff[f"{logE:.1f}"]['trig_weight'] = trig_weight

        # binned by veff
        total_veff[f"{logE:.1f}"]['czmins'] = zen_bins[:-1]
        total_veff[f"{logE:.1f}"]['total_veff_zen_bins'] = gwe(combined_total_veff_bins)
        total_veff[f"{logE:.1f}"]['deeporhybrid_only_veff_zen_bins'] = gwe(combined_deeporhybrid_only_veff_bins)
        total_veff[f"{logE:.1f}"]['shallow_only_veff_zen_bins'] = gwe(combined_shallow_only_veff_bins)
        total_veff[f"{logE:.1f}"]['dual_veff_zen_bins'] = gwe(combined_dual_veff_bins)

        print("total veff at 1 EeV for {} is {}".format(total_veff[f"{logE:.1f}"]['total_veff']/units.km**3 * 4 * np.pi, flavor))

    pkl_file_name = os.path.join('data/overlap_' + detector + "_deeptrig_" + deep_trigger + "_shallowtrig_" + shallow_trigger + "_mode_"+ mode+'_' + flavor + ".pkl")
    with open(pkl_file_name, "wb") as fout:
        pickle.dump(total_veff, fout, protocol=4)