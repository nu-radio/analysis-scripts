import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import h5py
import glob
import json
import os
import sys
import pickle
from NuRadioReco.utilities import units
from NuRadioMC.utilities import cross_sections
import helper
import argparse

"""
This script caluculates takes the output of the calc_fractions_advanced.py,
and makes (1) a csv file of the resultant veffs and fractions,
as well as a plot of the Veff and fractions.

"""


parser = argparse.ArgumentParser()
parser.add_argument('--det', required=True, help='which detector index')
parser.add_argument('--mode', required=True, help='which mode: deepshallow or hybridshallow')
args = parser.parse_args()

the_det = int(args.det)
the_mode = args.mode

detectors = [
    "baseline_array",
    "hex_hybrid_only_array",
    "hex_shallow_array",
    "hex_shallowheavy_array"
]

top = "/data/sim/Gen2/radio/2020/gen2-tdr-2021/simulation_output/secondaries_1700km2/step3"
detector = detectors[the_det]
config = "config_ARZ2020_noise"
detsim = "D01detector_sim"
deep_trigger = 'PA_4channel_100Hz'
shallow_trigger = 'LPDA_2of4_100Hz'
mode = the_mode

flavors = [
    "e",
    "mu",
    "tau"
]

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

result = {}
for flavor in flavors:
    result[flavor] = {}
    result[flavor]['E'] = []
    for l in the_list:
        result[flavor][l+'_veff'] = []
        result[flavor][l+'_aeff'] = []
    
    pkl_file_name = os.path.join('data/overlap_' + detector + "_deeptrig_" + deep_trigger + "_shallowtrig_" + shallow_trigger + "_mode_"+ mode+'_' + flavor + ".pkl")
    with open(pkl_file_name, "rb") as fin:
        coincidences = pickle.load(fin)
    
    for lgE, value in coincidences.items():
        result[flavor]['E'].append(lgE)
        lint = cross_sections.get_interaction_length(np.power(10., float(lgE)))

        for l in the_list:
            v_temp = value[l] # the veff water equivalent
            a_temp = v_temp/lint # get the aeff water equivalent
            v_temp = v_temp * 4 * np.pi / units.km**3 # convert to km^3 steradians
            a_temp = a_temp * 4 * np.pi / units.km**2 # convert to km^2 steradians

            result[flavor][l+'_veff'].append(v_temp)
            result[flavor][l+'_aeff'].append(a_temp)

energies = result['e']['E']
result_average = {}
for l in the_list:
    result_average[l+'_veff'] = np.zeros(9)
    result_average[l+'_aeff'] = np.zeros(9)
for l in the_list:
    for iflavor, flavor in enumerate(flavors):
        for iE, energy in enumerate(energies):
            result_average[l+'_veff'][iE] += result[flavor][l+'_veff'][iE]/3
            result_average[l+'_aeff'][iE] += result[flavor][l+'_aeff'][iE]/3

result_fractions = {}
for l in the_list:
    result_fractions[l] = result_average[l+'_veff']/result_average['total_veff']

fig, axs = plt.subplots(1, 2, figsize=(12,6))
# colors = ['C0', 'C1', 'C2']
# markers = ['o', 's', '^']
# styles = ['C0o-', 'C1s--', 'C2^-.', 'C3v:']
xx = result[flavor]['E']
fig.suptitle(f"{detector}")
# for l in [the_list[0]]:
for l in the_list:
    axs[0].plot(xx, result_average[l+'_aeff'], label=l)
    axs[1].plot(xx, result_fractions[l])

axs[0].set_yscale('log')
axs[0].set_xlabel("log10(energy [eV])")
axs[0].set_ylabel(r"[km$^2$ * str]")
axs[0].set_title("Water Equivalent Effective Area ")
axs[0].set_ylim([1E-5,1E2])
axs[0].legend(loc='lower right')


axs[1].set_ylim([0, 1.1])
axs[1].set_xlabel("log10(energy [eV])")
axs[1].set_ylabel("Fraction")
axs[1].set_title("Fraction of the All-Sky Effective Area")
fig.tight_layout(pad=2)
outfilename = f'plots/veff_fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.png'
fig.savefig(outfilename)

del fig, axs

import matplotlib as mpl

plt.rcParams.update({
    'xtick.labelsize': 14, 
    'ytick.labelsize': 14,
    'xtick.major.size': 5, 
    'ytick.major.size': 5,
    'axes.titlesize': 14,
    'axes.labelsize': 14
})
mpl.rc('font', size=12)
mpl.rc('axes', titlesize=14)

markers = ['o', 's', '^', 'v', 'd', '+', 'x', '>']

# plot the just the fractions in an easier to read format
fig, axs = plt.subplots(1, 1, figsize=(12,6))
for i, l in enumerate(the_list):
    axs.plot(xx, result_fractions[l], label=labels_dict[l], marker=markers[i])

axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs.set_ylim([0, 1.1])
axs.set_xlabel("log10(energy [eV])")
axs.set_ylabel("Fraction")
axs.set_title("Fraction of the All-Sky Effective Area")
fig.suptitle(f"{detector}")
fig.tight_layout()
outfilename = f'plots/fractions_{detector}_deeptrig_{deep_trigger}_shallowtrig_{shallow_trigger}_mode_{mode}.png'
fig.savefig(outfilename)

# now, we need a sensible way to accumulate some things
# for example, we might want to plot 
# 1) the hybrid-deep only
# 2) the hybrid-shallow only
# 3) the surface-only only
# 4) all cross-terms

the_list = ['total', 'hd', 'hs', 'ss', 'hd_hs', 'hd_ss', 'hs_ss', 'hd_hs_ss']
